/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2009
 *==========================================================================
 *
 * Author        : Ole J Gauteplass, Tandberg Telecom AS.
 * Co-author     : Ingvar Flaten Aarnes, Tandberg Telecom AS.
 *
 * Switches      : -
 *
 * Description   : Part of echocanceller. By estimating the aerl-level and
 *                 signal level from far-end+vcr etc. destined for loudspeaker,
 *                 the loudspeaker signal is adjusted to avoid saturation of
 *                 the (local) mic signal. (If the micsignal is distorted,
 *                 the echocanceller filters will not converge.
 *
 * Known issues:   If the 1/aerl is small and mic is heavily distorted
 *                 by a signal with more or less constant high level, the
 *                 echo-filters will not converge, and lslimiter will not
 *                 attenuate loudspeaker. (This is a rare situation).
 *
 * Note          : <notes here>
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *                 Modification made.
 **************************************************************************/

#include "lslimiter_priv.h"
#include "lslimiter.h"
#include "audec_defs.h"
#include "mathfun.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


LSLIMITER_PTR lslimiter_create(void)
{
    LSLIMITER_PTR pLslimiter;
    pLslimiter = (LSLIMITER_PTR) malloc(sizeof(LSLIMITER));
    if(pLslimiter == NULL)
    {
        fprintf(stderr, "Could not allocate Lslimiter struct\n");
    }
    return pLslimiter;
}


void lslimiter_destroy(LSLIMITER_PTR pLslimiter)
{
    free(pLslimiter);
}


void lslimiter_init(LSLIMITER_PTR pLslimiter, int numAcousticInputs)
{
    (void) numAcousticInputs;

    pLslimiter->gain_current_vector[0] = pLslimiter->gain_current_vector[1] = 1.0f;
    pLslimiter->stereoWeight      = 1.0f;
    pLslimiter->stereoWeightOld   = 1.0f;
    pLslimiter->forceMonoCntr     = 500;
    pLslimiter->limiter_threshold = 1.0f;

    pLslimiter->oldsp=0;
    pLslimiter->limcount=0;
    pLslimiter->limoldg = 1.0f;
    pLslimiter->lslimDelta = 0.6f;

    pLslimiter->lslimC = 1.0f;
    pLslimiter->lslimAdjust = 1.0f;
    pLslimiter->lslimDebug = false;

}

/* loadLsLimiter **********************************************************
 * Auth.: GVA
 * Desc.: loads pEq struct into iram
 * Returns: ptr to eq struct in iram
 ************************************************************************/
LSLIMITER_PTR loadLsLimiter(LSLIMITER_PTR pLsLimiterIn)
{
    /* MEMXFERMUSTDIE */
    return pLsLimiterIn;
}

/* flushEq **********************************************************
 * Auth.: GVA
 * Desc.: flush pEq struct to sdram
 ************************************************************************/
void flushLsLimiter(LSLIMITER_PTR pLsLimiterSdram, LSLIMITER_PTR pLsLimiterIram)
{
    /* MEMXFERMUSTDIE */
    (void)pLsLimiterSdram;
    (void)pLsLimiterIram;
}

/* STEREO COLLAPSING *********************************************************/
/* Auth.: GEO                                                                */
/*****************************************************************************/
void lslimiter_collapse_stereo_detect(float bef_lvl,
                                      float aft_lvl,
                                      float est_lvl,
                                      float noi_lvl,
                                      float stQty,
                                      LSLIMITER_PTR pLslimiter
                                      )

{
#define STATT            0.707946f
#define TOTAMP_TO_SATANP sqrtf(1.0f+1.0f/(2.0f*STATT*STATT))
    /* estimated amplitude of signal from all speakers divided by
       estimated amplitude of signal from satelite speakers */

    float b2alim   = 1.41f; /* float b2alim = 3.16227766f;  +10dB */
    float stonfact = 1.41f; /* 2.0f; // +6dB */
    int stfreeze   = 100;
    float stpow    = 1.25f;
    float pow_b2alim;
    float pow_stonfact;
    float beflev,org_beflev;
    float aftlev,org_aftlev;
    float estlev;
    float beta;

    (void)aft_lvl;
    (void)noi_lvl;
    (void)stQty;

    pLslimiter->stereoWeightOld = pLslimiter->stereoWeight;
    pow_b2alim   = powf(b2alim,stpow);
    pow_stonfact = powf(stonfact,stpow);
    beflev = powf(bef_lvl,stpow);
    estlev = powf(est_lvl,stpow);
    aftlev = beflev - estlev;
    if(aftlev < 0.0f)
    {
        aftlev = 0.0f;
    }
    if (stonfact>0.0f)
    {
        float a, b;

        org_beflev = beflev;
        org_aftlev = aftlev;
        beta = 0.0f;
        if (stfreeze<0)
        {
            beflev *= beflev;
            aftlev *= estlev;
        }

        a = beflev * TOTAMP_TO_SATANP;
        b = (aftlev * (beta + TOTAMP_TO_SATANP) - beta * beflev) * pow_b2alim;
        
        if (stfreeze==0)
        {
            if((int)(pLslimiter->stereoWeight * 100.0f) != 100)
            {
                if (a > b * pow_stonfact)
                    pLslimiter->stereoWeight = 1.0f;
            }
            else
            {
                if (a < b)
                    pLslimiter->stereoWeight = 0.0f;
            }
        }
        else
        {
            if((int)(pLslimiter->stereoWeight * 100.0f) != 100)
            {
                if (a > b * pow_stonfact)
                    pLslimiter->stereoWeight = (float) !(--(pLslimiter->forceMonoCntr));
            }
            else
            {
                if (a < b)
                {
                    pLslimiter->forceMonoCntr = abs(stfreeze);
                    pLslimiter->stereoWeight = 0.0f;
                }
            }
        }


        aftlev = org_aftlev;
        beflev = org_beflev;
    }
    pLslimiter->stereoWeight = 1.0f;  /*  no collapsing will be done here  */
}


/* LSLIMITER PROCESS ROUTINE ************************************************* */
/* Auth.: OGA                                                                  */
/* Limiting mono or stereo loudspeaker-signal and create stereo or             */
/* mono-signal by the same amount.                                             */
/* Performs loudspeaker limiter (to avoid overdrive in microphone),            */
/* and limiter (to avoid overdrive in loudspeaker).                            */
/* This is a soft-knee lookahead limiter with 5ms delay.                       */
/* *************************************************************************** */
void lslimiter_process(LSLIMITER_PTR lslimPtr,
                       float *inbuf[],
                       float *outbuf[],
                       float *delaylines[],
                       int nchannels,
                       float aerlInv_current,
                       float * limgain
                       )
{
    int n,i,j;
    float g,tmp;
    float tbuf[FRAMESIZE * 2];
    const float aerl_threshold = lslimPtr->lslimC * lslimPtr->lslimDelta * (aerlInv_current > EPSILON ? _rcpsp(aerlInv_current) : (1.0f/EPSILON));
    const float limthresh =  min(aerl_threshold, lslimPtr->limiter_threshold * lslimPtr->lslimAdjust * LSLIM_CLIP_HEADROOM);
    const float cliplevpos = lslimPtr->limiter_threshold;
    const float cliplevneg = -lslimPtr->limiter_threshold;
    float limth1=limthresh;
    float limth2=limthresh;
    float limratio=0;

    //if (aerl_threshold < limthresh) 
    {
        float knee=0.7f; //Soft compression starts 3dB below threshold
        limth1=limthresh*knee*0.9999f;
        limth2=fastdiv(limthresh,knee);
        limratio=fastdiv((knee-1.0f),(limth2-limth1));
    }

    if ((nchannels<=0) || (nchannels>LSLIMITER_MAXCHANNELS)) {
        return;
    }

    if (limgain==NULL)
        limgain = tbuf + FRAMESIZE;

    /* Calculate max(max(abs(x))) */
    {
    float *  inbufn = inbuf[0];
    for (i=0; i<FRAMESIZE; i++) {
        limgain[i] = fabsf(inbufn[i]);
    }
    }
    for (n=1; n<nchannels; n++) {
        float *  inbufn = inbuf[n];
        for (i=0; i<FRAMESIZE; i++) {
            limgain[i] = max(limgain[i],fabsf(inbufn[i]));
        }
    }

    /* Limiter */
    {
    float s,tmp2,tmp3;
    int itmp;
    tmp2=lslimPtr->oldsp;
    tmp3=lslimPtr->limoldg;
    itmp=lslimPtr->limcount;
    for (i=0; i<FRAMESIZE; i++) //17 cycles
    {
        s=limgain[i];
        if (s>tmp2) {
            tmp2=s;
            itmp=LSLIM_DELAY*5;
        }
        else if (itmp<=0) {
            tmp2=tmp2*(1.0f-LSLIM_REL)+s*LSLIM_REL;
            itmp=0;
        }
        else {
            itmp--;
        }

        if (tmp2>limth2) {
            g=fastdiv(limthresh,tmp2);  /* Hard knee */
        }
        else if (tmp2>limth1) {
            g=1.0f+limratio*(tmp2-limth1);  /* Soft knee */
        }
        else {
            g=1.0f;
        }

        tmp=g-tmp3;
        if (g>tmp3) {
            tmp3 += tmp*LSLIM_REL;
        }
        else {
            tmp3 += tmp*LSLIM_ATT;
        }
        limgain[i]=tmp3;
    }
    lslimPtr->limoldg=tmp3;
    lslimPtr->oldsp=tmp2;
    lslimPtr->limcount=itmp;
    lslimPtr->gain_current_vector[0] = lslimPtr->gain_current_vector[1] = lslimPtr->limoldg; /* Read by dbgdraw in aud_ectask.c */
    }

    for (n=0; n<nchannels; n++) {
        float *  delayline = delaylines[n];
        float * outbufn = outbuf[n];
        float * inbufn = inbuf[n];

        for (i=0; i<LSLIM_DELAY; i++)
        {
            tbuf[i] = inbufn[i]; /* Make a copy of inbuf sample in case inbuf == outbuf */
        }
        for (i=0,j=FRAMESIZE-LSLIM_DELAY; i<LSLIM_DELAY; i++,j++)
        {
            outbufn[i] = max(min(limgain[i]*delayline[i], cliplevpos), cliplevneg);
            delayline[i] = inbufn[j];  /* Copy most recent samples to delay-buffer */
        }
        for (i=LSLIM_DELAY,j=0; i<FRAMESIZE; i++,j++)
        {
            outbufn[i] = max(min(limgain[i]*tbuf[j], cliplevpos), cliplevneg);
        }
    }

    /* debug */
    if( lslimPtr->lslimDebug )
    {
        /*
        1. min limgain
        2. max inlevel channel 0
        3. max outlevel channel 0
        4. min aerl_threshold
        5. min limthresh
        */
        static int counter = 0;
        static float limgain_min=1e38;
        static float aerl_threshold_min=1e38, limthresh_min=1e38;
        static float inlevel_max=0, outlevel_max=0;

        for (n=0; n<nchannels; n++) {
            float *  inbufn = inbuf[n];
            float *  outbufn = outbuf[n];
            for (i=0; i<FRAMESIZE; i++) {
                inlevel_max = max(inlevel_max, fabsf(inbufn[i]));
                outlevel_max = max(outlevel_max, fabsf(outbufn[i]));
            }
        }
        for (i=0; i<FRAMESIZE; i++) {
            limgain_min = min(limgain_min,limgain[i]);
        }
        aerl_threshold_min = min(aerl_threshold_min, aerl_threshold);
        limthresh_min = min(limthresh_min, limthresh);

        if( ++counter == 20 )
        {
            printf("%8.2e,%8.2e,%8.2e,%8.2e,%8.2e  \n",
                   limgain_min,
                   inlevel_max, outlevel_max,
                   aerl_threshold_min, limthresh_min );
            limgain_min=1e38;
            inlevel_max=0, outlevel_max=0;
            aerl_threshold_min=1e38, limthresh_min=1e38;
            counter = 0;
        }
    }
}


/***************************************************************************
 *  Author:     : JPS
 *  Description : Set gatedependent loudspeaker treshold
 ***************************************************************************/
void lslimiter_setLslimAdjust(LSLIMITER_PTR lslimPtr, float val)
{
    lslimPtr->lslimAdjust = val;
    return;
}

/***************************************************************************
 *  Author:     : JPS
 *  Description : Set gatedependent aerl multiplication factor C
 ***************************************************************************/
void lslimiter_setLslimC(LSLIMITER_PTR lslimPtr, float val)
{
    lslimPtr->lslimC = val;
    return;
}

void lslimiter_status(LSLIMITER_PTR lslimPtr)
{
    printf("\rStatus: Lslimiter");
    printf("\r\n   lslimAdjust    - %4.2ef", lslimPtr->lslimAdjust);
    printf("\r\n   lslimC         - %4.2ef", lslimPtr->lslimC);
    printf("\r\n   lslimDebug     - %s", lslimPtr->lslimDebug ? "true" : "false");
    printf("\r\n");
    return;
}

/***************************************************************************
 *  Author:     : OJG
 *  Description : Turn dbgdraw on|off
 ***************************************************************************/
void lslimiter_setdebug(LSLIMITER_PTR lslimPtr, bool val)
{
    lslimPtr->lslimDebug = val;
    return;
}

#ifdef UNITTEST

#include "unittest.h"
#include "testdata/lslimiter_testdata.h"
#include <stdio.h>


/******************************************************************************
 * UNITTEST
 *
 * Description: Unittest of lslimiter. 1 s of audio data is run through the
 *              limiter using aerl_current = 1 and threshold = 1. The result
 *              is compared with the result of a MatLab simulation. If for each
 *              sample there is less than 0.17 dB difference the test passes.
 *
 *****************************************************************************/
static float lslimiterTestInput[]    = {LSLIMITER_TESTINPUT_DEFINE};
static float lslimTestOutput[]   = {LSLIMITER_TESTOUTPUT_DEFINE};
static float lslimTestOutputStereo[]   = {LSLIMITER_TESTOUTPUTSTEREO_DEFINE};

#define LSLIMITER_TESTLIM_LOW 0.91f
#define LSLIMITER_TESTLIM_HIGH 1.1f

void unittest_lslimiter()
{

    int i;
    int n;
    //FILE* filepointer = fopen("testresult.txt", "wb");
    float temp;

    LSLIMITER_PTR lslimiter;
    float *testResult = (float*)malloc(sizeof(float)*4800);
    float *testInputRight = (float*)malloc(sizeof(float)*4800);
    float *testResultMulti = (float*)malloc(sizeof(float)*4800*2);
    float *delay = (float *)malloc(sizeof(float)*LSLIM_DELAY*2);
    int differ = 0;

    lslimiter = lslimiter_create();

    /* Mono test */
    lslimiter_init(lslimiter,2);
    for(i=0; i<LSLIM_DELAY*2; i++)
    {
        delay[i] = 0.0f;
    }

    lslimiter->limiter_threshold = 1.0f;

    unittest_context("Lslimiter");

    //fprintf(filepointer,"Input     MatlabFasit    CModul \n");
    for (i=0; i<10; i++)
    {
        float *input[1];
        float *output[1];
        float *delaylines[1];
        input[0] = lslimiterTestInput+(i*FRAMESIZE);
        output[0] = testResult+(i*FRAMESIZE);
        delaylines[0] = delay;
        lslimiter_process(lslimiter,
                            input,
                            output,
                            delaylines,
                            1,
                            1,
                            NULL);

       if(i== 0)
       {
           for( n=0; n<LSLIM_DELAY; n++ )
           {
               if(testResult[i*FRAMESIZE+n] != 0)
               {
                   differ++;
               }
               //fprintf(filepointer,"%f,  %f,   %f,  0,   %d, \n", lslimiterTestInput[i*FRAMESIZE48 + n],lslimTestOutput[i*FRAMESIZE48 + n],testResult[i*FRAMESIZE48 + n],differ);
           }

           for(n=LSLIM_DELAY; n<FRAMESIZE; n++)
           {
               temp = testResult[i*FRAMESIZE+n]/lslimTestOutput[i*FRAMESIZE+n];
               if(temp < LSLIMITER_TESTLIM_LOW || temp > LSLIMITER_TESTLIM_HIGH )
               {
                   differ++;
               }
               //fprintf(filepointer,"%f,  %f,  %f,  %f,  %d \n", lslimiterTestInput[i*FRAMESIZE48 + n],lslimTestOutput[i*FRAMESIZE48 + n],testResult[i*FRAMESIZE48 + n],temp,differ);
           }
       }
       else
       {
           for(n = 0; n < FRAMESIZE; n++)
           {
               temp = testResult[i*FRAMESIZE+n]/lslimTestOutput[i*FRAMESIZE+n];
               if(temp < LSLIMITER_TESTLIM_LOW ||
                  temp > LSLIMITER_TESTLIM_HIGH )
               {
                   differ++;
               }

               //fprintf(filepointer,"%f,  %f,  %f,  %f,  %d\n", lslimiterTestInput[i*FRAMESIZE48 + n],lslimTestOutput[i*FRAMESIZE48 + n],testResult[i*FRAMESIZE48 + n],temp, differ);
           }
       }

    }

    /* Stereo test */
    lslimiter_init(lslimiter,2);
    for(i=0; i<LSLIM_DELAY*2; i++)
    {
        delay[i] = 0.0f;
    }

    lslimiter->limiter_threshold = 1.0f;

    for (i=0; i<4800; i++)
    {
      testInputRight[i] = lslimiterTestInput[i]*0.7;
    }

    //fprintf(filepointer,"Input     MatlabFasit    CModul \n");
    for (i=0; i<10; i++)
    {
        float *input[2];
        float *output[2];
        float *delaylines[2];
        
        input[0] = lslimiterTestInput+(i*FRAMESIZE);
        input[1] = testInputRight+(i*FRAMESIZE);
        output[0] = testResultMulti+(i*FRAMESIZE);
        output[1] = testResultMulti+4800+(i*FRAMESIZE);
        delaylines[0] = delay;
        delaylines[1] = delay + LSLIM_DELAY;
        lslimiter_process(lslimiter,
                          input,
                          output,
                          delaylines,
                          2,
                          1,
                          NULL);

       if(i== 0)
       {
           for( n=0; n<LSLIM_DELAY; n++ )
           {
               if((testResultMulti[i*FRAMESIZE+n] != 0) || (testResultMulti[i*FRAMESIZE+n+4800] != 0))
               {
                   differ++;
               }
               //fprintf(filepointer,"%f,  %f,  %f,  %f,  %f,  %f,  0,  %d\n", lslimiterTestInput[i*FRAMESIZE48 + n],testInputRight[i*FRAMESIZE48 + n],lslimTestOutputStereo[i*FRAMESIZE48 + n],lslimTestOutputStereo[i*FRAMESIZE48 + n+4800],testResultMulti[i*FRAMESIZE48 + n],testResultMulti[i*FRAMESIZE48 + n+4800], differ);
           }

           for(n=LSLIM_DELAY; n<FRAMESIZE; n++)
           {
               temp = testResultMulti[i*FRAMESIZE+n]/lslimTestOutputStereo[i*FRAMESIZE+n];
               if(temp < LSLIMITER_TESTLIM_LOW || temp > LSLIMITER_TESTLIM_HIGH )
               {
                   differ++;
               }
               temp = testResultMulti[i*FRAMESIZE+n+4800]/lslimTestOutputStereo[i*FRAMESIZE+n+4800];
               if(temp < LSLIMITER_TESTLIM_LOW || temp > LSLIMITER_TESTLIM_HIGH )
               {
                   differ++;
               }
               //fprintf(filepointer,"%f,  %f,  %f,  %f,  %f,  %f,  %f,  %d\n", lslimiterTestInput[i*FRAMESIZE48 + n],testInputRight[i*FRAMESIZE48 + n],lslimTestOutputStereo[i*FRAMESIZE48 + n],lslimTestOutputStereo[i*FRAMESIZE48 + n+4800],testResultMulti[i*FRAMESIZE48 + n],testResultMulti[i*FRAMESIZE48 + n+4800], temp, differ);
           }
       }
       else
       {
           for(n = 0; n < FRAMESIZE; n++)
           {
               temp = testResultMulti[i*FRAMESIZE+n]/lslimTestOutputStereo[i*FRAMESIZE+n];
               if(temp < LSLIMITER_TESTLIM_LOW ||
                  temp > LSLIMITER_TESTLIM_HIGH )
               {
                   differ++;
               }
               temp = testResultMulti[i*FRAMESIZE+n+4800]/lslimTestOutputStereo[i*FRAMESIZE+n+4800];
               if(temp < LSLIMITER_TESTLIM_LOW ||
                  temp > LSLIMITER_TESTLIM_HIGH )
               {
                   differ++;
               }

               //fprintf(filepointer,"%f,  %f,  %f,  %f,  %f,  %f,  %f,  %d\n", lslimiterTestInput[i*FRAMESIZE48 + n],testInputRight[i*FRAMESIZE48 + n],lslimTestOutputStereo[i*FRAMESIZE48 + n],lslimTestOutputStereo[i*FRAMESIZE48 + n+4800],testResultMulti[i*FRAMESIZE48 + n],testResultMulti[i*FRAMESIZE48 + n+4800], temp, differ);
           }
       }

    }
    //fclose(filepointer);

    unittest_assert(differ == 0);

    free(testResultMulti);
    free(testInputRight);
    lslimiter_destroy(lslimiter);
    free(testResult);
}

#endif
