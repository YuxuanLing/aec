/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 * Co-author     : Geir Ole Øverby (GEO), Tandberg Telecom AS
 *
 * Switches      : <COMPILER SWITCH HERE>
 *           <add description here>
 *
 * Description   : <add description here>
 *
 * Note      : <notes here>
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *         Modification made.
 **************************************************************************/

#include "lsexchndl_priv.h"
#include "micexchndl_priv.h"
#include "micexchndl.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static float micexchndl_power(float power, float * subbuf);

static void micexchndl_xcorr(LSEXCHNDL_PTR const  pLsexchndl[],
                             MICEXCHNDL_PTR pMicexchndl,
                             float * subbuf);

static void micexchndl_sumxc(MICEXCHNDL_PTR pMicexchndl, float * sumxc);

static void micexchndl_detect(LSEXCHNDL_PTR pLsexchndl[],
                              MICEXCHNDL_PTR pMicexchndl,
                              float * sumxc);

static void micexchndl_filtgain(MICEXCHNDL_PTR pMicexchndl);


/***************************************************************************
 MICEXCHNDL_CREATE
 *   Description: Allocates memory for struct MICEXCHNDL.
 *
 *   Parameters:   -
 *
 *   Returns:      pointer to struct MICEXCHNDL.
 *
 *   Note:         -
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
MICEXCHNDL_PTR micexchndl_create(void)
{

    MICEXCHNDL_PTR pMicexchndl;

    pMicexchndl = (MICEXCHNDL_PTR ) malloc(sizeof(MICEXCHNDL));

    if(pMicexchndl == NULL)
    {
        fprintf(stderr,"micexchndl_create: Could not allocate micexchndl buffer\n");
    }

    return pMicexchndl;
}


/***************************************************************************
 * MICEXCHNDL_DESTROY
 *   Description:  Frees allocated memory pointed by input parameter.
 *
 *   Parameters:   pMicexchndl
 *
 *   Returns:      -
 *
 *   Note:         -
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
void micexchndl_destroy(MICEXCHNDL_PTR pMicexchndl)
{
    free(pMicexchndl);
}


/***************************************************************************
 * MICEXCHNDL_INIT
 *  Description:  Initialises members in struct micexchndl.
 *
 *  Parameters:   pMicexchndl
 *
 *  Returns:      -
 *
 *  Note:         -
 *
 *  Globals used: -
 *
 *  Restrictions: Memory pointed by pMicexchndl shall be allocated by
 *                calling function micexchndl_create
 ***************************************************************************/
void micexchndl_init(MICEXCHNDL_PTR  pMixexchndl, int channels)
{
    int i, j, ch;

    pMixexchndl->channels = channels;
    pMixexchndl->excgain=(float)1.0;
    pMixexchndl->excbefbef = 0.0f;
    pMixexchndl->excaftaft = 0.0f;
    pMixexchndl->excgainu  = 1.0f;
    pMixexchndl->excgainf  = 1.0f;
    pMixexchndl->dlreadix = 0; /* where to write present data in excaftls (which row) */

    for (ch = 0; ch < channels; ch++)
    {
        for(i = 0; i < EXCLAGCNT; i++)
        {
            for(j = 0; j < EXCSUBUSED * 2; j++)
            {
                pMixexchndl->excaftls[ch][i][j]=(float)0.0;
            }
        }
    }

    pMixexchndl->ExptHndlOn = true;
}


/***************************************************************************
 * MICEXCHNDL_PROCESS
 *   Description:   Main micexchndl routine
 *
 *   Parameters:    pLsexchndl     - pointer to LSEXCHNDL struct
 *                  pMicexchndl    - pointer to MICEXCHNDL struct.
 *                  micfft         - output from analyse filter
 *                  estecho        - estimated echo from echocomp
 *
 *   Returns:       -
 *
 *   Note:          -
 *
 *   Globals used:  -
 *
 *   Restrictions:
 ***************************************************************************/
void micexchndl_process(LSEXCHNDL_PTR pLsexchndl[],
                        MICEXCHNDL_PTR pMicexchndl,
                        float *  micfft,
                        float *  estecho)
{
  int i, j;
  float diff_micest[EXCSUBUSED * 2];
  float sumxc[MAX_NUM_CHANNELS][EXCLAGCNT];

 //float subused[EXCSUBUSED * 2];
 //float sumxcLSBEF[EXCLAGCNT];
 //float sumxcAttack[EXC_FIR_LAG];

  /* calculates the difference between micfft and estecho */
  j = (SUBUSED_START * 2) + (EXCSUBUSED_START * 2);
  for (i = 0; i < EXCSUBUSED * 2; i++)
  {
      diff_micest[i] = micfft[j] - estecho[j];
      //subused[i] = micfft[j];
      j++;
  }

  pMicexchndl->excbefbef = micexchndl_power(pMicexchndl->excbefbef, &micfft[(SUBUSED_START + EXCSUBUSED_START) * 2]);
  pMicexchndl->excaftaft = micexchndl_power(pMicexchndl->excaftaft, &diff_micest[0]);

  micexchndl_xcorr(pLsexchndl, pMicexchndl, &diff_micest[0]);
  micexchndl_sumxc(pMicexchndl, sumxc[0]);
  micexchndl_detect(pLsexchndl, pMicexchndl, sumxc[0]);
  micexchndl_filtgain(pMicexchndl);

/*************************************
 * //Needs FIR implementation
 * if(ExcDelayCalc)
 * {
 *     float dummyVar;
 *     dummyVar = micexchndl_powerBef(pMicexchndl->excbefbef, subused, pMicexchndl);
 *     (void)dummyVar;
 *
 *     //-- cross correlate the loudspeaker signal with the diff_micest
 *     micexchndl_xcorrFir(pLsexchndl, pMicexchndl, diff_micest, subused);
 *
 *     //shift delayline of aft and bef signals
 *     micexchndl_shift(pMicexchndl,diff_micest,subused);
 *     micexchndl_SNR(estecho, pMicexchndl);
 *     micexchndl_sumxcFirDelay(pMicexchndl, sumxcLSBEF);
 *     micexchndl_CalcDelay(pLsexchndl,   // (I) LS Exception handler
 *                          pMicexchndl, // (I/O) Mic Exception handler
 *                          sumxcLSBEF);
 *
 * }
 *************************************/



  if (pMicexchndl->ExptHndlOn) /* testing for debug */
  {
      pMicexchndl->excgain = pMicexchndl->excgainf; // excgain doesn't really need to be in the struct...but it's simpler...!
  }
  else
  {
      pMicexchndl->excgain = 1;
  }

}

/*************************************************************************************
 * MICEXCHNDL_POWER
 *   Auth.: IFA
 *   Desc.: mic power calculation routine. Calculates an approximation
 *          of the power of micfft, based on the most significant subbands for speech.
 *          subbuf must point to the first "SUBUSED" band in micfft
 ************************************************************************************/
float micexchndl_power(float power, float * subbuf)
{
    int i;
    float tmp = 0.0f;
#ifdef __INTEL_COMPILER
#pragma novector
#endif
    for (i = 0; i < (EXCSUBUSED) * 2; i++)
    { /* squared l2-norm of exch-subbands of micfft */
        tmp += subbuf[i] * subbuf[i]; // sum(buf(params.excsubix).*conj(buf(params.excsubix)))
    }

    return( ((1.0f - EXCPSPD) * power) + (EXCPSPD * tmp) ) ; /* returns *new* power */
}

/*************************************************************************************
 * MICEXCHNDL_XCORR
 *  Auth.:  IFA
 *  Desc.:  Cross correlation calculation.
 *          Calculates the correlation of each row of the ls-signal-delayline
 *          and present echo-cancelled mic-signal.
 *          Calculates column-wise, to avoid reading subbuf several times...
 *          subbuf = subused bands of micfft-estecho.
 ************************************************************************************/
void micexchndl_xcorr(LSEXCHNDL_PTR const  pLsexchndl[], MICEXCHNDL_PTR pMicexchndl, float * subbuf)
{
    int i, j, k, ch;
    float tmpRe = 0.0f;
    float tmpIm = 0.0f;
    float a, b, c, d; /* to avoid extra reading */
    for (j = 0; j < EXCSUBUSED * 2; j += 2)
    {
        a = subbuf[j]; /* Re1 in */
        b = subbuf[j+1]; /* Im1 in */
        k = pMicexchndl->dlreadix; /* excdline is not rotated, dlwriteix-1 denotes newest row. */
        for (i = 0; i < EXCLAGCNT; i++)
        {
            for (ch = 0; ch < pMicexchndl->channels; ch++)
            {
                c = pLsexchndl[ch]->excdline[k][j];  /* Re2 in */
                d = pLsexchndl[ch]->excdline[k][j+1]; /* Im2 in */

                tmpRe = a * c + b * d;
                tmpIm = b * c - a * d;

                pMicexchndl->excaftls[ch][i][j]   = (1.0f - EXCPSPD) * (pMicexchndl->excaftls[ch][i][j])  +  (EXCPSPD * tmpRe);
                pMicexchndl->excaftls[ch][i][j+1] = (1.0f - EXCPSPD) * (pMicexchndl->excaftls[ch][i][j+1]) + (EXCPSPD * tmpIm);
            }
            ((k == 0) ? (k = (EXCLAGCNT - 1)):(k--)); /* next row, modulo EXCLAGCNT, (counts backwards) */
            /* note that excdline has opposite orientation than in the MatLab model.
               model: [new, 2., 3., ...] here: (ex:) [4., 3., 2, new, 10., 9., ...]
               (Still, excaftls has the same orientation as ML-model.)*/
        }
    }

    pMicexchndl->dlreadix++;
    if (pMicexchndl->dlreadix == EXCLAGCNT)
    {
        pMicexchndl->dlreadix = 0;
    }
    /* dlreadix rotates modulo EXCLAGCNT, it points to where 'newest' data should be read from excdline. */

}
/******************************************************************************
 * MICEXCHNDL_SUMXC
 *    Auth.: IFA
 *    Desc.: Cross correlation sum
 *           Calculates the niveau of each (delayed) rom of excaftls.
******************************************************************************/
void micexchndl_sumxc(MICEXCHNDL_PTR pMicexchndl, float *  sumxc)
{
    int i, j, ch;
    float tmp;

    for (ch = 0; ch < pMicexchndl->channels; ch++)
    {
        float tmp1 = 0.0f;
        float tmp2 = 0.0f;

        for (i = 0; i < EXCLAGCNT; i++)
        {
#ifdef __INTEL_COMPILER
#pragma novector
#endif
            for (j = 0; j < EXCSUBUSED * 2; j += 2)
            {
                tmp1 += fabsf(pMicexchndl->excaftls[ch][i][j]);
                tmp2 += fabsf(pMicexchndl->excaftls[ch][i][j+1]);
            }
            tmp                   = tmp1 + tmp2;
            sumxc[i+ch*EXCLAGCNT] = tmp * tmp;
            tmp1                  = 0.0f;
            tmp2                  = 0.0f;
        }
    }
}

/******************************************************************************
 * MICEXCHNDL_DETECT
 *   Auth.: IFA
 *   Desc.: calculates whether an exception has occured: excgainu
 ****************************************************************************/
void micexchndl_detect(LSEXCHNDL_PTR pLsexchndl[], MICEXCHNDL_PTR pMicexchndl, float * sumxc)
{
    int i, j, ch;
    float caft, cbef;
    float tmp;
    float exccond1[MAX_NUM_CHANNELS];
    float exccond2[MAX_NUM_CHANNELS];

    for (ch = 0; ch < pMicexchndl->channels; ch++)
    {
        exccond1[ch] = 0.0f;
        exccond2[ch] = 0.0f;
    }

    if ((int)(pMicexchndl->excgainu * 100.0f) == 100)
    { // exception last buffer
        caft = EXCCAFTEXCLAST * pMicexchndl->excaftaft;
        cbef = EXCCBEFEXCLAST * pMicexchndl->excbefbef;
    }
    else
    {                       // no exception last buffer
        caft = EXCCAFTNOEXCLAST * pMicexchndl->excaftaft;
        cbef = EXCCBEFNOEXCLAST * pMicexchndl->excbefbef;
    }

    for (ch = 0; ch < pMicexchndl->channels; ch++)
    {
        j = pLsexchndl[ch]->lsLastIx;
        exccond1[ch] = pLsexchndl[ch]->exclsls[j] - (EXCMINLSPOW);
        for (i = 0; i < EXCLAGCNT; i++)
        {
            float exclsls = pLsexchndl[ch]->exclsls[j];

            tmp = sumxc[i+ch*EXCLAGCNT] - caft * exclsls - cbef * exclsls;
            if (tmp > 0)
            {
                exccond2[ch] += tmp;
            } // sum(max(0,tmp))
            if (j == 0)
            {
                j = EXCLAGCNT; // rotates modulo EXCLAGCNT; backwards.
            }
            j--;
        }
    }

    pMicexchndl->excgainu = 1;
    for (ch = 0; ch < pMicexchndl->channels; ch++)
    {
        if ((exccond1[ch] > 0.0f) && (exccond2[ch] > 0.0f))
        {
            pMicexchndl->excgainu = 0;
            break;
        }
    }
}


/******************************************************************************
 * MICEXCHNDL_FILTGAIN
 *   Auth.: IFA
 *   Desc.: Gain filtering,
 *          filters (updating slowly form excgainu) the exceptionhandler gain : excgainf
 *          matlab: state.excgainf = state.excgainf + EXCGAINFILTSP * (state.excgainu - state.excgainf);
*****************************************************************************************/
void micexchndl_filtgain(MICEXCHNDL_PTR pMicexchndl)
{
    pMicexchndl->excgainf = (1 - EXCGAINFILTSP) * pMicexchndl->excgainf + EXCGAINFILTSP * pMicexchndl->excgainu;
}



/* micexchndl_load **********************************************************
 * Auth.: JPS
 * Desc.: loads pMicexchndl struct into iram
 * Returns: ptr to pMicexchndl struct in iram
 ************************************************************************/
MICEXCHNDL_PTR micexchndl_load(MICEXCHNDL_PTR pMicexchndlIn)
{
    /* MEMXFERMUSTDIE */
    return pMicexchndlIn;
}

/* micexchndl_flush **********************************************************
 * Auth.: JPS
 * Desc.: flush pMicexchndl struct to sdram
 ************************************************************************/
void micexchndl_flush(MICEXCHNDL_PTR pMicexchndlSdram, MICEXCHNDL_PTR pMicexchndlIram)
{
    /* MEMXFERMUSTDIE */
    (void)pMicexchndlSdram;
    (void)pMicexchndlIram;
}

/* micexchndl_getExcgain *************************************************
 * Auth.: JPS
 * Desc.: returns pMicexchndl->excgain
 ************************************************************************/
float micexchndl_getExcgain(MICEXCHNDL_PTR pMicexchndl)
{
    return pMicexchndl->excgain;
}



/****************************************************************************
 * Functions for setting testvariables inside noisereduction
 ***************************************************************************/
void micexchndl_set(MICEXCHNDL_PTR pMicexchndl, bool onoff)
{
    const char *status[] = {"off", "on"};
    pMicexchndl->ExptHndlOn = onoff;
    printf("MicExchndl is set %s \n\n", (char*)status[pMicexchndl->ExptHndlOn]);
    return;
}

void micexchndl_statusMicexchndl(MICEXCHNDL_PTR pMicexchndl)
{
    const char *status[] = {"off", "on"};

    printf("\rStatus: Micexchndl");
    printf("\r\n   micexchndl     - %s", (char*)status[pMicexchndl->ExptHndlOn]);
    printf("\r\n");
}








#ifdef UNITTEST
/**************************************************************************
 * UNITTEST_MICEXCHNDL
 *   Description.: simple run-trhough of micexchndl with random input
 *                 and check against testvector verified in matlab.
 *                 No difference for 16kHz or 48kHz
 **************************************************************************/
#include "unittest.h"
#include "testdata/micexchndltestdata.h"  /* testvectors defined here */

/* Private testdata */
static float micfftTestInput[]     = {MICFFT_TESTINPUT_DEFINE};
static float estechoTestInput[]    = {ESTECHO_TESTINPUT_DEFINE};
static float lsexchndlTestInput1[] = {LSEXCHNDL_TESTINPUT_CH1_DEFINE};
static float lsexchndlTestInput2[] = {LSEXCHNDL_TESTINPUT_CH2_DEFINE};
static float excbefbefTestOutput[] = {EXCBEFBEF_TESTOUTPUT_DEFINE};
static float excaftaftTestOutput[] = {EXCAFTAFT_TESTOUTPUT_DEFINE};
static float excaftlsTestOutput1[] = {EXCAFTLS_TESTOUTPUT_CH1_DEFINE};
static float excaftlsTestOutput2[] = {EXCAFTLS_TESTOUTPUT_CH2_DEFINE};
/**************************************************************************
 * COMPARERESULTS
 *    Desc.: Compares two float buffers. Due to compilator optimization, a
 *           a small difference 'epsilon' is allowed when comparing with the
 *           testvector
 **************************************************************************/
static int compareResults(float * inp1, float * inp2, int len)
{
    int differ = 0;
    int i;
    float epsilon = 1e-10f;

    for(i = 0; i < len; i++)
    {
        if((fabsf(inp1[i] - inp2[i])) > epsilon)
        {
            differ = 1;
            //printf("%d \t %23.18f \t %23.18f \t %23.18f \n ",i,inp1[i],inp2[i], fabsf(inp1[i] - inp2[i]));
        }

    }

    return differ;
}

/**************************************************************************
 * UNITTEST_MICEXCHNDL
 **************************************************************************/
#include <string.h>
void unittest_micexchndl()
{
    MICEXCHNDL_PTR pMicexchndl;
    LSEXCHNDL_PTR pLsexchndl[2];
    float * estecho;
    float * pexcaftls1;       //pointer to output buffer
    float * pexcaftls2;       //pointer to output buffer
    float * pbef;
    float * paft;
    int i, numFrames;
    int excbefbefDiffer, excaftaftDiffer, excaftlsDiffer;
	int size = 7680;     //number of samples in input
    int subused = 116;
    int fftsize = 256;

    /*FILE *fid;
    char *Filename ="d:\\tmp\\micexchndltestoutput.dat"; */

    //criteria flag to be returned by function. 0: OK, 1: ERROR, 2:Nothing done, hence ERROR
    excbefbefDiffer = excaftaftDiffer = excaftlsDiffer = 2; 

    pMicexchndl = micexchndl_create();
    if(pMicexchndl == NULL)
    {
        unittest_assert(0);
        return;
    }
    micexchndl_init(pMicexchndl, 2);

    pLsexchndl[0] = lsexchndl_create();
    lsexchndl_init(pLsexchndl[0]);
    pLsexchndl[1] = lsexchndl_create();
    lsexchndl_init(pLsexchndl[1]);

    numFrames = size/fftsize;

    //allocating memory for output buffer and initializing to zero values
    pexcaftls1 = (float *) calloc(EXCLAGCNT * EXCSUBUSED * 2 * numFrames, sizeof(float));
    pexcaftls2 = (float *) calloc(EXCLAGCNT * EXCSUBUSED * 2 * numFrames, sizeof(float));
    pbef = (float *) calloc(numFrames, sizeof(float));
    paft = (float *) calloc(numFrames, sizeof(float));

    estecho =(float *) malloc(sizeof(float) * fftsize);
    for(i = 0; i < fftsize; i++)
    {
        estecho[i] = 0.0f;
    }

    for (i = 0; i < numFrames; i++)
    {
        memcpy(estecho + 2 * SUBUSED_START, estechoTestInput + i * subused * 2, sizeof(float) * subused * 2);
        scratchmem_lock();
        lsexchndl_process(pLsexchndl[0], lsexchndlTestInput1 + i * fftsize, 1.0f);
        lsexchndl_process(pLsexchndl[1], lsexchndlTestInput2 + i * fftsize, 1.0f);
        micexchndl_process(pLsexchndl,pMicexchndl, micfftTestInput + i * fftsize, estecho);
        scratchmem_unlock();
        // printf("%23.18f %23.18f %23.18f \n", pMicexchndl->excgain, pMicexchndl->excbefbef, pMicexchndl->excaftaft);

        pbef[i] = pMicexchndl->excbefbef;
        paft[i] = pMicexchndl->excaftaft;
        memcpy(pexcaftls1 + i * EXCLAGCNT * EXCSUBUSED * 2, pMicexchndl->excaftls[0], sizeof(float) * EXCLAGCNT * EXCSUBUSED * 2);
        memcpy(pexcaftls2 + i * EXCLAGCNT * EXCSUBUSED * 2, pMicexchndl->excaftls[1], sizeof(float) * EXCLAGCNT * EXCSUBUSED * 2);
    }

    /*
      printf("write micexchndl output to file \n");
      fid = (void*)fopen(Filename,"wb");
      //fwrite(pexcaftls,sizeof(float),numFrames*EXCLAGCNT*EXCSUBUSED*2,fid);
      fwrite(pbef,sizeof(float),numFrames,fid);
      //fwrite(paft,sizeof(float),numFrames,fid);
      fclose(fid);
      free(fid);
      free(Filename);*/


    excbefbefDiffer = compareResults(pbef, excbefbefTestOutput, numFrames);
    excaftaftDiffer = compareResults(paft, excaftaftTestOutput, numFrames);
    excaftlsDiffer  = compareResults(pexcaftls1, excaftlsTestOutput1, numFrames * EXCLAGCNT * EXCSUBUSED * 2);
    excaftlsDiffer  = compareResults(pexcaftls2, excaftlsTestOutput2, numFrames * EXCLAGCNT * EXCSUBUSED * 2);



    free(pexcaftls1);
    free(pexcaftls2);
    free(pbef);
    free(paft);
    free(estecho);
    lsexchndl_destroy(pLsexchndl[0]);
    lsexchndl_destroy(pLsexchndl[1]);
    micexchndl_destroy(pMicexchndl);

    unittest_context("excbef");
    unittest_assert(excbefbefDiffer == 0);
    unittest_context("excaft");
    unittest_assert(excaftaftDiffer == 0);
    unittest_context("excaftls");
    unittest_assert(excaftlsDiffer == 0);
}

#endif /* UNITTEST */
