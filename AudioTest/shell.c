/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Switches      : <COMPILER SWITCH HERE>
 *                 <add description here>
 *
 * Description   : Delaying fullband micsignal outside filterbank to be used
 *                 as clean mic signal when there is near-end talk only.
 *                 It improves the quality of the micsignal sent to far-end.
 *                 When echo, double-talk or noise situation, the output
 *                 from the filterbank must be used for best possible cancellation.
 *
 *                 The control parameter is pShell->gain,
 *                    if 0-ish: use filterbank signal
 *                    if 1-ish: use shell signal
 *
  * ------------------------------------------------------------------------
 * Major changes made
 * ------------------------------------------------------------------------
 * 2009-08-30    : JPS
 *                 Ported from calypso branch and reduced complexity.
 *                 Removed some unnecessary hifi-stuff because we now have 48kHz filterbank.
 **************************************************************************/
#include "shell_priv.h"
#include "shell.h"
#include "mathfun.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "noisereduction_priv.h"
#include "analyse_priv.h"
#include "lsprocess_priv.h"

#ifdef __ARM_NEON__
#include <arm_neon.h>
//#pragma thumb           // for compatibility with gcc
#endif

#define ONE_DIV_FRAMESIZE 0.0020833333f
#define SHELLGVAR_LIMIT  1.33f
#define SHELLGVAR_INCSPEED 0.2f  // up-speed
#define SHELLGVAR_DECSPEED 0.8f  // down-speed
#define SHELL_LSLEVSPEED 0.125f
#define SHELL_SUBUSED_START 4   //avoid using lowest subbands in summation of levels

SHELL_PTR shell_create(void)
{
    SHELL_PTR pShell =(SHELL_PTR) malloc(sizeof(SHELL));
    if ( pShell == NULL )
    {
        fprintf(stderr, "shell_create: Could not allocate pShell\n");
        return 0;
    }

    return pShell;
}

void shell_destroy(SHELL_PTR pShell)
{
    free(pShell);
}

void shell_init(SHELL_PTR pShell, int channels)
{
    int i, m;

    for (i = 0; i < SHELL_TOTALDELAY; i++)
    {
        pShell->delayline[i] = 0.0f;
    }

    pShell->variable_gain = 0.0f;
    pShell->gain          = 0.0f;
    pShell->prevgain      = 0.0f;
    pShell->writeIx       = 0;
    pShell->rampIx        = FFTSIZE;
    pShell->outIx         = FRAMESIZE;
    pShell->debug         = 0;
    pShell->lslevAdjust   = 3.0f;
    pShell->noisegnAdjust = 1.6f;
    pShell->nlpgAdjust    = 1.0f;
    
    pShell->varIncSpeed   = SHELLGVAR_INCSPEED;
    pShell->varDecSpeed   = SHELLGVAR_DECSPEED;

    for (m = 0; m < MAX_NUM_CHANNELS; m++)
    {
        for (i = 0; i < SUBUSED16; i++)
        {
            pShell->lslevel[m][i] = 0.0f;
        }
    }

    if (channels >1) 
    {
        pShell->shellOn   = false;
    }
    else
    {
        pShell->shellOn   = true;
    }

}



/*  Fills in new micdata into delayline;
*   starting at writeIx, emulating circular buffer.
*
*   pShell->delayline is in SDRAM, so we have to make sure the loop is interruptable.
*   Thus, adding pragma 
*/
void shell_input(SHELL_PTR pShell, float *  inputbuf)
{
    int wIx = pShell->writeIx;
    float * buf = pShell->delayline;
    int loopend = FRAMESIZE;

    if( wIx + FRAMESIZE > SHELL_TOTALDELAY)
        loopend = SHELL_TOTALDELAY - wIx;
    
    memcpy(&buf[wIx], inputbuf, loopend*sizeof(float));
    
    if(loopend == FRAMESIZE)
    {
        wIx += FRAMESIZE;
        if(wIx==SHELL_TOTALDELAY)
          wIx = 0;
    }
    else
    {
        // After process one part of memory, process another part of memory which rollback.
        wIx = FRAMESIZE-loopend;
        memcpy(buf, &inputbuf[loopend], wIx*sizeof(float)); 
    }
    
    pShell->writeIx = wIx;
}


/*  Adds output from shell to synth_out.
*   starting at outIx, emulating circular buffer */
void shell_output(SHELL_PTR pShell, float *  outBuf)
{
    int i;
    int outIx = pShell->outIx;
    float * buf = pShell->delayline;

#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
    int loopend = FRAMESIZE;
    float32x4x2_t tmpd0,tmpd1,tmpd2,tmpd3, tmpd4;

    /*
    Firstly process 16x in loop. 
    If all of outIx are NOT beyond SHELL_TOTALDELAY, the routine will return.
    Otherwise, it continues to do the remain. 
      a) do 16 element to assure the OutIx rolls back to header of buffer.
      b) process the remain with 16x.
    */
    if( outIx + FRAMESIZE > SHELL_TOTALDELAY)
    {
        loopend = (SHELL_TOTALDELAY - outIx) & ~15;
    }
            
    for (i = 0; i < loopend; i+=16, outIx+=16)
    {
        tmpd0 = vld2q_f32(&buf[outIx]);
        tmpd1 = vld2q_f32(&outBuf[i]);            
        tmpd2 = vld2q_f32(&buf[outIx+8]);
        tmpd3 = vld2q_f32(&outBuf[i+8]);
        tmpd4.val[0] = vaddq_f32(tmpd0.val[0], tmpd1.val[0]);
        tmpd4.val[1] = vaddq_f32(tmpd0.val[1], tmpd1.val[1]);
        tmpd2.val[0] = vaddq_f32(tmpd2.val[0], tmpd3.val[0]);
        tmpd2.val[1] = vaddq_f32(tmpd2.val[1], tmpd3.val[1]);

        vst2q_f32(&outBuf[i],   tmpd4);
        vst2q_f32(&outBuf[i+8], tmpd2);
    }
    
    if(outIx==SHELL_TOTALDELAY)
        outIx = 0; 

    if( i < FRAMESIZE)
    {
        loopend = i + 16;
        for ( ; i < loopend; i++)
        {
            outBuf[i] += buf[outIx];
            outIx++;
            if(outIx == SHELL_TOTALDELAY)
                outIx = 0;
        } 
        for ( ; i < FRAMESIZE; i+=16, outIx+=16)
        {
            tmpd0 = vld2q_f32(&buf[outIx]);
            tmpd1 = vld2q_f32(&outBuf[i]);            
            tmpd2 = vld2q_f32(&buf[outIx+8]);
            tmpd3 = vld2q_f32(&outBuf[i+8]);
            tmpd4.val[0] = vaddq_f32(tmpd0.val[0], tmpd1.val[0]);
            tmpd4.val[1] = vaddq_f32(tmpd0.val[1], tmpd1.val[1]);
            tmpd2.val[0] = vaddq_f32(tmpd2.val[0], tmpd3.val[0]);
            tmpd2.val[1] = vaddq_f32(tmpd2.val[1], tmpd3.val[1]);

            vst2q_f32(&outBuf[i],   tmpd4);
            vst2q_f32(&outBuf[i+8], tmpd2);
        }       
    }

#else
    for (i = 0; i < FRAMESIZE; i++)
    {
        outBuf[i] += buf[outIx];
        outIx = (outIx == ((SHELL_TOTALDELAY) - 1)) ? 0 : outIx+1;
    }
#endif

    pShell->outIx = outIx;
}

        
/* multiplies the part of the delayline corresponding to the present gain, by gain */
static void shell_rampgain(SHELL_PTR pShell)
{
    int i;
    float gainslope;
    float invBufLength = ONE_DIV_FRAMESIZE;
    float *  buf = pShell->delayline;
    float prevgain = pShell->prevgain;
    float a = (pShell->gain - prevgain) * invBufLength; //stigningstall for gainslope
    int rampIx = pShell->rampIx;

#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
    float         a_mul_i[4] = {0.0f, a, 2*a, 3*a};
    float32x4_t   a_mul_4_n = vdupq_n_f32(4*a);
    float32x4_t   prevgain_n = vdupq_n_f32(prevgain);
    float32x4_t   gainslope_n = vld1q_f32(a_mul_i);
    float32x4_t   gainslope_plus_n;

    float32x4_t   tmpd0, tmpd1;
    float32x4_t   tmp0, tmp1;
    int loopend = FRAMESIZE;
    float k; 

    /* --- the pseudocode ---
    if( rampIx + FRAMESIZE > SHELL_TOTALDELAY)
    {
        loopend = (SHELL_TOTALDELAY - rampIx) & ~7;
    }
    
    for (i = 0; i < loopend; i+=8, rampIx+=8)
    {
        loop 8 {
            gainslope = a * i + prevgain;
            buf[rampIx] *= gainslope;
        }
    }
    
    if(loopend < FRAMESIZE)
    {
        loopend = i+8;
        for ( ; i < loopend; i++, rampIx++)
        {
            gainslope = a * i + prevgain;
            buf[rampIx] *= gainslope;
            rampIx++;
            if(rampIx == SHELL_TOTALDELAY)
                rampIx = 0;
        }
        
        for ( ; i < FRAMESIZE; i+=8, rampIx+=8)
        {
            loop 8 {
                gainslope = a * i + prevgain;
                buf[rampIx] *= gainslope;
            }
        } 
    }
    */

    gainslope_n = gainslope_plus_n = vaddq_f32(prevgain_n, gainslope_n);

    if( rampIx + FRAMESIZE > SHELL_TOTALDELAY)
    {
        loopend = (SHELL_TOTALDELAY - rampIx) & ~7;
    }
    
    tmpd0 = vld1q_f32(&buf[rampIx]);
    tmp0  = vmulq_f32(gainslope_n, tmpd0);

    k = 8.0f;
    for (i = 0; i < loopend; i+=8, rampIx+=8, k+=8.0f)
    {
        float a_mul_i = k * a;
        gainslope_n = vaddq_f32(gainslope_n, a_mul_4_n);
        tmpd1 = vld1q_f32(&buf[rampIx+4]);
        tmp1  = vmulq_f32(gainslope_n, tmpd1);

        vst1q_f32(&buf[rampIx], tmp0);
        gainslope_n = vdupq_n_f32(a_mul_i);
        tmpd0 = vld1q_f32(&buf[rampIx+8]);
        /* 
        gainslope_n is assigned ((k*a)<<96)|((k*a)<<64)|((k*a)<<32)|(k*a), and it need to 
        add ((3*a+prevgain)<<96)| ((2*a+prevgain)<<64) | ((a+prevgain)<<32) | (prevgain) 
        which equals original prevgain_n plus original gainslpoe_n.
        */
        gainslope_n = vaddq_f32(gainslope_n, gainslope_plus_n);
        vst1q_f32(&buf[rampIx+4], tmp1); 
        tmp0  = vmulq_f32(gainslope_n, tmpd0);             
    }

    if(rampIx==SHELL_TOTALDELAY)
        rampIx = 0; 

    if(loopend < FRAMESIZE)
    {
        if(rampIx&7)
        {
            gainslope = i * a + prevgain;
            loopend = 8+i;
            
            for ( ; i < loopend; i++)
            {
                buf[rampIx] *= gainslope;
                gainslope += a;
                rampIx++;
                if(rampIx == SHELL_TOTALDELAY)
                    rampIx = 0;    
            }
            gainslope_n = vdupq_n_f32(gainslope);
        }

        tmpd0 = vld1q_f32(&buf[rampIx]);
        tmp0  = vmulq_f32(gainslope_n, tmpd0);

        k = i+8.0f;
        for ( ; i < FRAMESIZE; i+=8, rampIx+=8, k+=8.0f)
        {
            float a_mul_i = k * a;
            gainslope_n = vaddq_f32(gainslope_n, a_mul_4_n);
            tmpd1 = vld1q_f32(&buf[rampIx+4]);
            tmp1  = vmulq_f32(gainslope_n, tmpd1);

            vst1q_f32(&buf[rampIx], tmp0);
            gainslope_n = vdupq_n_f32(a_mul_i);
            tmpd0 = vld1q_f32(&buf[rampIx+8]);
            gainslope_n = vaddq_f32(gainslope_n, gainslope_plus_n);
            vst1q_f32(&buf[rampIx+4], tmp1); 
            tmp0  = vmulq_f32(gainslope_n, tmpd0);             
        }
    }
#else
    for (i=0; i<FRAMESIZE; i++)
    {
        gainslope = a * i + prevgain;
        buf[rampIx] *= gainslope;
        rampIx = (rampIx == ((SHELL_TOTALDELAY) - 1)) ? 0 : rampIx+1;
    }
#endif

    pShell->rampIx = rampIx;
    pShell->prevgain = pShell->gain;
    return;
}

static float shell_lslevel(float * lev, COMPLEX32*  lsBuf)
{
    int m;
    float sp = SHELL_LSLEVSPEED;
    float sum = 0.0f;

#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now 
    float32x4x2_t tmpd0;
    float32x4_t   tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
    float32x4_t   sp_half_n = vdupq_n_f32(0.5f*sp);
    float32x4_t   sp_rsb1_n = vdupq_n_f32(1.0f-sp);
    float32x4_t   sum_n = vdupq_n_f32(0.0f);
    float32x2_t   tmpX;

    tmpd0 = vld2q_f32((float*)&lsBuf[SHELL_SUBUSED_START]);
    tmp0 = vabsq_f32(tmpd0.val[0]);
    tmp1 = vabsq_f32(tmpd0.val[1]);
    tmp2 = vld1q_f32(&lev[SHELL_SUBUSED_START]);
    tmp3 = vaddq_f32(tmp0, tmp1);
    tmp4 = vmulq_f32(tmp2, sp_rsb1_n);
    tmp5 = vmulq_f32(tmp3, sp_half_n);
    for (m = SHELL_SUBUSED_START; m < SUBUSED16; m+=4)
    {
        tmpd0 = vld2q_f32((float*)&lsBuf[m+4]);
        tmp0 = vabsq_f32(tmpd0.val[0]);
        tmp1 = vabsq_f32(tmpd0.val[1]);
        tmp6 = vaddq_f32(tmp4, tmp5);
        tmp2 = vld1q_f32(&lev[m+4]);   
        tmp3 = vaddq_f32(tmp0, tmp1);
        sum_n = vaddq_f32(sum_n, tmp6);
        tmp4 = vmulq_f32(tmp2, sp_rsb1_n);
        tmp5 = vmulq_f32(tmp3, sp_half_n);  
        vst1q_f32(&lev[m], tmp6);
    }
    
    tmpX = vpadd_f32(vget_high_f32(sum_n), vget_low_f32(sum_n));
    sum = vget_lane_f32(tmpX, 1) + vget_lane_f32(tmpX, 0);    
#else
    for (m = SHELL_SUBUSED_START; m < SUBUSED16; m++)
    {
        float re = lsBuf[m].re;
        float im = lsBuf[m].im;

        float abs_re = _fabsf(re);
        float abs_im = _fabsf(im);
        lev[m] = lev[m] * (1.0f - sp) + 0.5f * sp * (abs_re + abs_im);
        sum += lev[m];
    }
#endif
    return (sum);
}

static void shell_calcgain(SHELL_PTR pShell, NOIRED_PTR pNoired, LSPROCESS_PTR pLs)
{
  float lsParam = 1.0f, lsAdjust = 0.0f, gaintmp = 0.0f, noisgn = 0.0f;
  float aerlInv = 1.0f; //dummyvariable - should aerl be implemented?
  float nlpg = pNoired->shellnlp;
  float sum_lslevel = 0.0f, tmp = 0.0f;
  float sum_totlev = 0.0f, sum_beflev = 0.0f, sum_noilev = 0.0f;
  int i, m;

  for (m = SHELL_SUBUSED_START; m < SUBUSED16; m++)
  {
    sum_totlev += pNoired->totlev[m];
    sum_beflev += pNoired->beflev[m];
    sum_noilev += pNoired->noilev[m];
  }

  for (i = 0; i < pLs->nChannels; i++)
  {
      LSPROCESS_CHANNEL *channel = &pLs->channel[i];

      tmp = shell_lslevel(pShell->lslevel[i], channel->pAnalyseSdram->fftout);
      sum_lslevel = max(sum_lslevel, tmp);
  }

  noisgn = min(1.0f,max(0.0f, pShell->noisegnAdjust *(1.0f - 1.25f * fastdiv(sum_noilev, (sum_beflev+ EPSILON)) ))) ;

  lsParam = pShell->lslevAdjust * aerlInv * fastdiv(sum_lslevel, (sum_totlev + EPSILON));
  lsAdjust = min(1.0f,max(0.0f, 1.1f*(1.0f - lsParam)));

  nlpg = min(1.0f, nlpg*pShell->nlpgAdjust);

  gaintmp = nlpg * lsAdjust * noisgn;
  gaintmp = max(0.0f, gaintmp);     //*1.5; %Strekking: HifiSize

  if (pShell->variable_gain <= gaintmp)
  {
      pShell->variable_gain = ( (1.0f-pShell->varIncSpeed) * pShell->variable_gain ) + ( pShell->varIncSpeed * gaintmp ) ;
  }
  else
  {
      pShell->variable_gain = ( (1.0f-pShell->varDecSpeed) * pShell->variable_gain ) + ( pShell->varDecSpeed * gaintmp ) ;
  }
  pShell->variable_gain = (min(SHELLGVAR_LIMIT,(pShell->variable_gain - 0.05f)*1.1f ));

  if (pShell->shellOn)
  {
      pShell->gain = min(1.0f,max(0.0f, pShell->variable_gain));
  }
  else
  {
      pShell->gain = 0.0f;
  }

  if (pShell->debug == 67)
  {
      static int dbgcounter = 0;
      static float sum1 = 0.0f, sum2 = 0.0f, sum3 = 0.0f, sum4 = 0.0f, sum5 = 0.0f, sum6 = 0.0f, sum7 = 0.0f, sum8 = 0.0f, sum9 = 0.0f, sum10 = 0.0f;

      /*
       1. sum_beflev
       2. sum_noilev
       3. sum_totlev
       4. sum_lslevel
       5. nlpg
       6. lsAdjust
       7. noisgn
       8. gaintmp
       9. variable_gain
       10. gain
      */

      sum1 += sum_beflev;
      sum2 += sum_noilev;
      sum3 += sum_lslevel;
      sum4 += sum_totlev;
      sum5 += nlpg;
      sum6 += lsAdjust;
      sum7 += noisgn;
      sum8 += gaintmp;
      sum9 += pShell->variable_gain;
      sum10 += pShell->gain;

      if (++dbgcounter == 20)
      {
          printf("  %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n",
                 sum1*0.05, sum2*0.05, sum3*0.05, sum4*0.05, sum5*0.05, sum6*0.05, sum7*0.05, sum8*0.05, sum9*0.05, sum10*0.05);
          sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = 0.0f, sum9 = 0.0f, sum10 = 0.0f;
          dbgcounter = 0;
      }
  }
return;
}



SHELL_PTR shell_load(SHELL_PTR pShellSdram)
{
    /* MEMXFERMUSTDIE */
    return pShellSdram;
}

void shell_flush(SHELL_PTR pShellSdram, SHELL_PTR pShellIram)
{
    /* MEMXFERMUSTDIE */
    (void)pShellSdram;
    (void)pShellIram;
}



/* Ramp current micbuffer
 * Mixes outbuf signal with delayed input.
 * Starting at writeIx, emulating circular buffer and updated writeIx
 */
float shell_process(SHELL_PTR pShell, NOIRED_PTR pNoired, LSPROCESS_PTR pLs)
{
    shell_calcgain(pShell, pNoired, pLs);
    shell_rampgain(pShell);
    return(pShell->gain);
}

void shell_setShellOn(SHELL_PTR pShell, bool onoff)
{
    const char *status[] = {"off", "on"};
    pShell->shellOn = onoff;
    printf("shell is set %s\n\n", (char*)status[pShell->shellOn]);
    return;
}

void shell_setlslevAdjust(SHELL_PTR pShell, float value)
{
    pShell->lslevAdjust = value;
    printf("pShell->lslevAdjust is set to %4.2f \n", pShell->lslevAdjust);
    return;
}

void shell_setnoisegnAdjust(SHELL_PTR pShell, float value)
{
    pShell->noisegnAdjust = value;
    printf("pShell->noisegnAdjust is set to %4.2f \n", pShell->noisegnAdjust);
    return;
}

void shell_setnlpgAdjust(SHELL_PTR pShell, float value)
{
    pShell->nlpgAdjust = value;
    printf("pShell->nlpgAdjust is set to %4.2f \n", pShell->nlpgAdjust);
    return;
}

void shell_setvarIncSpeed(SHELL_PTR pShell, float value)
{
    pShell->varIncSpeed = value;
    printf("pShell->varIncSpeed is set to %4.2f \n", pShell->varIncSpeed);
    return;
}

void shell_setvarDecSpeed(SHELL_PTR pShell, float value)
{
    pShell->varDecSpeed = value;
    printf("pShell->varDecSpeed is set to %4.2f \n", pShell->varDecSpeed);
    return;
}

void shell_setDebug(SHELL_PTR pShell, int value)
{
    pShell->debug = value;
    printf("pShell->debug is set to %d \n", pShell->debug);
}

void shell_status(SHELL_PTR pShell)
{
    const char *status[] = {"off", "on"};
    printf("\rStatus: shell");
    printf("\r\n   shellOn           - %s",    (char*)status[pShell->shellOn]);
    printf("\r\n   lslevAdjust       - %4.2f",  pShell->lslevAdjust);
    printf("\r\n   noisegnAdjust     - %4.2f",  pShell->noisegnAdjust);
    printf("\r\n   nlpgAdjust        - %4.2f",  pShell->nlpgAdjust);
    printf("\r\n   varIncSpeed       - %4.2f",  pShell->varIncSpeed);
    printf("\r\n   varDecSpeed       - %4.2f",  pShell->varDecSpeed);
    printf("\r\n   debug             - %d",     pShell->debug);
    printf("\r\n");
}

#ifdef UNITTEST
/************************************************************************
 * UNITTEST
 * Desc.: run-trhough of shell
 ************************************************************************/
#include "unittest.h"
#include "testdata/shelltestdata.h"      /* input and output data defined here */
#include <math.h>

#define nWRITE_UNITTEST_OUTPUT

/* Private testdata */
static float shellTestInput[]      = {SHELL_TESTINPUT_DEFINE};
static float shellTestOutput[]     = {SHELL_TESTOUTPUT_DEFINE};
static float shellTestOutputGain[] = {SHELL_TESTOUTPUTGAIN_DEFINE};

/*************************************************************************
 * COMPARERESULT
 *      Desc.: Compares two float buffers
 *************************************************************************/
static int compareResults(float *inp1, float *inp2, int len)
{
    int differ = 0;
    int i;
    float epsilon = 3.0e-7f;

    for (i = 0; i < len && differ != 1; i++)
    {
        if ((fabsf(inp1[i] - inp2[i])) > epsilon)
        {
            differ = 1;
        }
    }
    return differ;
}
/*************************************************************************
 * UNITTEST_SHELL
 ************************************************************************/
void unittest_shell()
{
    SHELL_PTR pShell;
    NOIRED_PTR pNoired;
    LSPROCESS_PTR pLs;
    int i;
    int numFrames = 10;
    int differ = 0;
    int differGain = 0;
    const char  * testName = "Shell";
    float * buf;
    FILTERBANK_USE_TYPE filterbankUse = FILTERBANK_EC;
    float gain;

#ifdef WRITE_UNITTEST_OUTPUT
    int cnt;
    FILE *fid;
    const char *Filename = "/home/jps/src/testdata/shelltestoutput.dat";
    FILE *fid2;
    const char *FilenameGain = "/home/jps/src/testdata/shelltestoutputgain.dat";

#endif


    /* Initialization*/
    unittest_context(testName);

    pShell = shell_create();
    if ( (pShell == NULL) )
    {
        fprintf(stderr, "shell unittest: Could not allocate shell struct \n");
        unittest_assert(0);
        return;
    }
    shell_init(pShell, 1);

    pNoired = noisereduction_create();
    noisereduction_init(pNoired, filterbankUse);

    pLs = lsprocess_create(filterbankUse, 1);
    lsprocess_init(pLs, filterbankUse);
    memset(pLs->pAnalyseSdram[0]->fftout, 0, sizeof(float) * FFTSIZE);

    buf = (float *) calloc(FRAMESIZE, sizeof(float));
    pNoired->shellnlp = 1.0f;

#ifdef WRITE_UNITTEST_OUTPUT
        printf("write shell outputs to files \n");

        fid = (void*)fopen(Filename,"wb");
        if (fid == NULL)
        {
            fprintf(stderr, "unittest_shell: Could not open file %s \n", Filename);
        }

        fid2 = (void*)fopen(FilenameGain,"wb");
        if (fid2 == NULL)
        {
            fprintf(stderr, "unittest_shell: Could not open file %s \n", FilenameGain);
        }

#endif

    /* Run process and verify result*/
    for (i = 0; i < numFrames; i++)
    {
        SHELL_PTR pShellIram;
        scratchmem_lock();

        memcpy(buf, shellTestInput + i * FRAMESIZE, sizeof(float) * FRAMESIZE);

        shell_input(pShell, buf);

        pShellIram = shell_load(pShell);
        memxfer_waitAddress(pShellIram);

        gain = shell_process(pShellIram, pNoired, pLs);
        memset(buf, 0, sizeof(float)*FRAMESIZE);
        shell_output(pShellIram, buf);

        shell_flush(pShell, pShellIram);

        differ += compareResults(buf, shellTestOutput + i * FRAMESIZE, FRAMESIZE);
        differGain += compareResults(&gain, shellTestOutputGain + i, 1);

        // dummy update of noisereduction levels independent of shell-functions
        pNoired->sum_beflev += 0.02f;
        pNoired->sum_noilev += 0.001f;
        pNoired->sum_totlev += 0.01f;
        pNoired->shellnlp   -= 0.001f;

        scratchmem_unlock();

#ifdef WRITE_UNITTEST_OUTPUT
        cnt = fwrite(buf, sizeof(float), FRAMESIZE, fid);
        printf("Wrote %d samples to file %s \n", cnt, Filename);
        cnt = fwrite(&gain, sizeof(float), 1, fid2);
        printf("Wrote %d samples to file %s \n", cnt, FilenameGain);
    }
    fclose(fid);
#else
    }
#endif

   unittest_assert(differ == 0);
   unittest_assert(differGain == 0);

   noisereduction_destroy(pNoired);
   lsprocess_destroy(pLs, filterbankUse);
   shell_destroy(pShell);
   free(buf);
}

#endif /* UNITTEST */

