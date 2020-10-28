/***************************************************************************
*                              A U D I O
*--------------------------------------------------------------------------
*                  (C) Copyright Tandberg Telecom AS 2004
*==========================================================================
* Revision      : $Revision: 1.18 $
*
* Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
* Co-author     : -
*
* Switches      : <COMPILER SWITCH HERE>
*         <add description here>
*
* Description   : Part of echocanceller. Main objectives:
*                 1) Reduce noise by decreasing the gain in subbands not
*          containing talk.
*                 2) Supress echoes by attenuating signal when far-end-talk.
*                 3) Mute signal when echocanceller is not working (exception
*          occurs)
*                 4) Fill in synthetic noise when 2) or 3) mutes signal.
*                 5) Produce some input to mixer and hifigain.
*
* Note    : <notes here>
*
* Documentation : <Xxxxxx-Document.doc>
*
**************************************************************************/
#include "dereverb.h"
#include "noisereduction_priv.h"
#include "noisereduction.h"
#include "keyclickremoval.h"
#include "mathfun.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define SUBUSED16 116              /* used for estimates without high-frequencies */

#ifdef WRITE_NOIREDDEBUG_TO_FILE
	FILE * fidNoired;
#endif

float noisereduction_measureLevel1(float speed,
                                   float * const  levelPtr,
                                   const COMPLEX32 *bufPtr);

float noisereduction_measureLevel2(float speed,
                                   float * const  levelPtr,
                                   const COMPLEX32 *bufPtrA,
                                   const COMPLEX32 *bufPtrB);

short noisereduction_checknoise(NOIRED_PTR pNoired);

float noisereduction_estimatenoise(float decspeed,
                                   float incspeed,
                                   const float * const  totlev,
                                   float *noilev);

void noisereduction_calcnoisegainlimit(NOIRED_PTR pNoired, bool noisereductionOn);

void noisereduction_calcnoisegain(float *  noisgn,
                                  const float *noilev,
                                  const float *totlev,
                                  float nrglimit);

void noisereduction_calcnlpgain(NOIRED_PTR pNoired,
                                float excgain,
                                float *nlpgain,
                                float nlpconst);

void noisereduction_selectgain(const float *noisegain,
                               const float *dereverbgain,
                               const float *nlpgain,
                               float * const  gain,
                               bool dereverbOn);

void noisereduction_calccomfortnoise(NOIRED_PTR pNoired,
                                            const float *gain,
                                            COMPLEX32 *noise,
                                            float *nfgn);

/***************************************************************************
NOISEREDUCTION_CREATE
*   Description: Allocates memory for struct noired.
*
*   Parameters:   -
*
*   Returns:      struct noired.
*
*   Note:         -
*
*   Globals used: -
*
*   Restrictions: -
***************************************************************************/
NOIRED_PTR noisereduction_create(void)
{
    NOIRED_PTR pNoired;

    pNoired = (NOIRED_PTR) malloc(sizeof(NOIRED));

    if( pNoired == NULL )
    {
        fprintf(stderr, "noisereduction_create: Could not allocate noisereductions buffer\n");
        return 0;
    }

#ifdef WRITE_NOIREDDEBUG_TO_FILE
	fidNoired = fopen("noiredDebug.txt", "w");
#endif
    return pNoired;
}

/***************************************************************************
* NOISEREDUCTION_DESTROY
*   Description: Frees allocated memory pointed by ptr.
*
*   Parameters:   pNoired
*
*   Returns:      -
*
*   Note:         -
*
*   Globals used: -
*
*   Restrictions: -
***************************************************************************/
void noisereduction_destroy(NOIRED_PTR pNoired)
{
    free(pNoired);
	#ifdef WRITE_NOIREDDEBUG_TO_FILE
	fclose(fidNoired);
	#endif

}

/***************************************************************************
* NOISEREDUCTION_INIT
*  Description:  Initialises members in struct noired.
*
*  Parameters:   pNoired
*
*  Restrictions: Memory pointed by pNoired shall be allocated by
*                calling function noisereduction_create.
***************************************************************************/
void noisereduction_init(NOIRED_PTR pNoired, FILTERBANK_USE_TYPE filterbankUse)
{
    int i = 0;

    for( i = 0; i < SUBUSED; i++ )
    {
        pNoired->beflev[i] = 0.0f;
        pNoired->aftlev[i] = 0.0f;
        pNoired->estlev[i] = 0.0f;
        pNoired->totlev[i] = 0.0f;
        pNoired->noilev[i] = 0.0f;
        pNoired->noisgn[i] = 1.0f;
    }

    pNoired->sum_beflev = 0.0f;
    pNoired->sum_aftlev = 0.0f;
    pNoired->sum_estlev = 0.0f;
    pNoired->sum_totlev = 0.0f;
    pNoired->sum_noilev = 0.0f;
    pNoired->fullgPrev  = 0.0f;
    pNoired->fullgPrevHi  = 0.0f;

    for( i = 0; i < 2; i++ )
    {
        pNoired->beflevPrev[i] = 0.0f;
        pNoired->beflevPrevHi[i] = 0.0f;
    }

    for( i = 0; i < SUBUSED; i++ )
    {
        pNoired->subgPrev[i] = 0.0f;
    }
    pNoired->nrgoptimal = NRGLMAX;
    pNoired->noiseamp = (float)0.0;

    //if( variant == AL_CONNECTORVAR_HANDSET || variant == AL_CONNECTORVAR_HEADSET )
    if(filterbankUse == FILTERBANK_EC)
        pNoired->comNoiseAmp = 1.0f;
    else
        pNoired->comNoiseAmp = 1.12f;   /* tuned to zerodelaynlp-routine */

    pNoired->ncnt = 0;
    pNoired->ncnt_clim = NOISECLIM;
    pNoired->ncnt_boot = 2 * NOISECLIM;
    pNoired->seed = 0x10001;


    pNoired->NlpGainOn    = true;
    pNoired->dereverbOn   = false;
    pNoired->ComNoiseOn   = true;
    pNoired->NoiseRedOn   = true;
    pNoired->NLPSubOn     = true;
    pNoired->NlpHiOn      = true;
    pNoired->fullgNlpOn   = true;
    pNoired->fullgNlpHiOn = true;
    pNoired->extragainOn  = true;
    pNoired->fullgNlpTransition= 10;
    pNoired->debug        = 0;
    pNoired->nlpSubusedEndHi = SUBUSED16;

    pNoired->shellnlp = 0.0f;
    pNoired->nlpspeed = NLPSPEED;

//==========================================================
// Code to detect CPU capabilities on IA processors
//==========================================================
#ifdef __INTEL_COMPILER
    pNoired->MMX_PRESENT   = 0;
    pNoired->SSE_PRESENT   = 0;
    pNoired->SSE2_PRESENT  = 0;
    pNoired->SSE3_PRESENT  = 0;
    pNoired->SSSE3_PRESENT = 0;
    pNoired->SSE4_PRESENT  = 0;
    {
#define MMX_FLAG 0x0800000
#define SSE_FLAG 0x2000000   // Safe to use -QxK/-xK
#define SSE2_FLAG 0x4000000  // Safe to use -QxN/-xN
#define SSE3_FLAG 0x0000001  // Safe to use -QxP/-xP
#define SSSE3_FLAG 0x0000200 // Safe to use -QxT/-xT
        unsigned int cpeinfo;
        unsigned int cpsse3;
        unsigned long vender[3];
        __asm {
            mov eax, 0
            cpuid
            mov long ptr [vender+0], ebx
            mov long ptr [vender+4], edx
            mov long ptr [vender+8], ecx
        }
        __asm {
            mov eax, 01h
            cpuid
            mov cpeinfo, edx
            mov cpsse3, ecx
        }
        if( strncmp((char *)vender,"GenuineIntel",12) ) cpsse3 = 0;
        pNoired->MMX_PRESENT   = (cpeinfo & MMX_FLAG   )>0;
        pNoired->SSE_PRESENT   = (cpeinfo & SSE_FLAG   )>0;
        pNoired->SSE2_PRESENT  = (cpeinfo & SSE2_FLAG  )>0;
        pNoired->SSE3_PRESENT  = (cpsse3  & SSE3_FLAG  )>0;
        pNoired->SSSE3_PRESENT = (cpsse3  & SSSE3_FLAG )>0;
    }
#endif
}


/***************************************************************************
NOISEREDUCTION_LOADNOIRED
*   Description: copies the noisereduction struct to iram  using dmax.
*
*   Parameters:   - pNoisereduction, pointer to NOIRED struct
*
*   Returns:      - pointer to noired struct in iram.
*
*   Note:         - should use memxfer_waitAddress on the returned pointer
*                   before use in order to confirm completed dmax transfer
***************************************************************************/
NOIRED_PTR noisereduction_loadNoired(NOIRED_PTR pNoisereduction)
{
    /* MEMXFERMUSTDIE */
    return pNoisereduction;
}

/***************************************************************************
NOISEREDUCTION_FLUSHNOIRED
*   Description: flushes the noisereduction struct to sdram  using dmax.
*
*   Parameters:   - pNoiredSdram
*                 - pNoiredIram
***************************************************************************/
void noisereduction_flushNoired(NOIRED_PTR pNoiredSdram, NOIRED_PTR pNoiredIram)
{
    /* MEMXFERMUSTDIE */
    (void)pNoiredSdram;
    (void)pNoiredIram;
}



/***************************************************************************
NOISEREDUCTION_LOADGAIN
*   Description: copies the gain vector to iram  using dmax.
*
*   Parameters:   - gainSdram, pointer to gain vector
*
*   Returns:      - pointer to gain vector in iram.
*
*   Note:         - should use memxfer_waitAddress on the returned pointer
*                   before use in order to confirm completed dmax transfer
***************************************************************************/
float* noisereduction_loadGain(float *gainSdram)
{
    /* MEMXFERMUSTDIE */
    return gainSdram;
}


/***************************************************************************
* NOISEREDUCTION_PROCESS
*  Description: Main noisereduction routine
*               Determines if the input signals are noisy or not.
*
*   Parameters: fftout            - pointer to fftout from ANALYSE struct of mic
*               excgain           - exceptionhandler gain from MICEXCHNDL struct.
*               pNoired           - pointer to NOIRED struct containing variables
*                                   defining noise reduction filters.
*               estecho           - Echo estimation from echo compensation.
*               gain              - gain output from noise reduction.
*               noise             - noise output from noise reduction.
*               nlpconst          - non linear processing const.
*               noisereductionOn  - Noisereduction turned on/off
*               GSM_noiseDetected -
*
*   Returns:    1 - There is noise
*               0 - No noise
*
*   Note:       -
*
*   Globals used: -
*
*   Restrictions: Memory pointed by pointers shall be allocated by
*                 calling function.
***************************************************************************/
short noisereduction_process(COMPLEX32 *  fftout,
                             float excgain,
                             NOIRED_PTR pNoired,
                             DEREVERB_PTR pDereverb,
                             COMPLEX32 *  estecho,
                             float *  gain,
                             COMPLEX32 *  noise,
                             float nlpconst,
                             bool  noisereductionOn,
                             float *decay)
{

    DEREVERB_PROCESS dereverbProcess;
    short isnoise = 0;
    float nlpgain[SUBUSED];
    float nfgn[SUBUSED];

    /* Decay buffer and DEREVERB_PTR to IRAM */
    if (pNoired->dereverbOn){    
        dereverb_loadDereverb(pDereverb, &dereverbProcess);
    }

    /* ---- normal process routine: ------------------------------------*/
    /* update level vectors: */
    pNoired->sum_beflev = noisereduction_measureLevel1(pNoired->nlpspeed,
                                                       pNoired->beflev,
                                                       fftout); // beflev

    pNoired->sum_estlev = noisereduction_measureLevel1(pNoired->nlpspeed,
                                                       pNoired->estlev,
                                                       estecho); // estlev

    pNoired->sum_aftlev = noisereduction_measureLevel2(pNoired->nlpspeed,
                                                       pNoired->aftlev,
                                                       fftout,
                                                       estecho); // aftlev

    pNoired->sum_totlev = noisereduction_measureLevel2(NRSPEED,
                                                       pNoired->totlev,
                                                       fftout,
                                                       estecho); // totlev

    /* estimate noise level (based on noisedetection) */
    isnoise = noisereduction_checknoise(pNoired);

    if( isnoise == 1 )
    {
        pNoired->sum_noilev = noisereduction_estimatenoise(ISNOISE_DECSPEED,
                                                           ISNOISE_INCSPEED,
                                                           pNoired->totlev,
                                                           pNoired->noilev);
    }
    else
    {
        pNoired->sum_noilev = noisereduction_estimatenoise(NONOISE_DECSPEED,
                                                           NONOISE_INCSPEED,
                                                           pNoired->totlev,
                                                           pNoired->noilev);
    }

    /* calculation of noise gain limit */
    noisereduction_calcnoisegainlimit(pNoired, noisereductionOn);
    /* calculation of noise gain */
    noisereduction_calcnoisegain(pNoired->noisgn,
                                 pNoired->noilev,
                                 pNoired->totlev,
                                 pNoired->nrglimit);
    /* calculation of nlp gain */
    noisereduction_calcnlpgain(pNoired,
                               excgain,
                               nlpgain,
                               nlpconst);

    /* Calculation of dereverb gain */
    if (pNoired->dereverbOn) {
        
        dereverb_process(&dereverbProcess,
                         (float *)fftout,
                         decay,
                         pNoired->nrglimit);
    }

    noisereduction_selectgain(pNoired->noisgn, 
                              decay, 
                              nlpgain, 
                              gain,
                              pNoired->dereverbOn);

    if (pNoired->dereverbOn) {
        dereverb_flushDereverb(pDereverb, &dereverbProcess);
    }

    /* calculation of comfort noise */
    noisereduction_calccomfortnoise(pNoired, gain, noise, nfgn);
#ifdef WRITE_NOIREDDEBUG_TO_FILE
	fprintf(fidNoired, "%13.11e %13.11e %13.11e %13.11e\n", pNoired->sum_beflev, pNoired->sum_estlev, pNoired->sum_aftlev, pNoired->sum_noilev);
#endif

    return isnoise;
}

void noisereduction_setNcntBoot(NOIRED_PTR pNoired, int iCntBoot)
{
    pNoired->ncnt_boot = (short)iCntBoot;
}

/****************************************************************************
* Functions for setting testvariables inside noisereduction
***************************************************************************/
void noisereduction_setNlp(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->NlpGainOn = onoff;
    printf("Nlp is set %s\n\n", status[pNoired->NlpGainOn]);
    return;
}

void noisereduction_setDereverb(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->dereverbOn = onoff;
    printf("Dereverb is set %s\n\n", status[pNoired->dereverbOn]);
    return;
}

void noisereduction_setComNoise(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->ComNoiseOn = onoff;
    printf("ComNoise is set %s\n\n", status[pNoired->ComNoiseOn]);
    return;
}

void noisereduction_setComNoiseAmp(NOIRED_PTR pNoired, float value)
{
    pNoired->comNoiseAmp = value;
    printf("pNoired->comNoiseAmp is set %5.3f \n", pNoired->comNoiseAmp);
}


void noisereduction_setNoiseRed(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->NoiseRedOn = onoff;
    printf("NoiseRed is set %s\n\n", status[pNoired->NoiseRedOn]);
    return;
}

void noisereduction_setNlpSub(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->NLPSubOn = onoff;
    printf("NlpSub is set %s \n\n", status[pNoired->NLPSubOn]);
    return;
}

void noisereduction_setNlpHi(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->NlpHiOn = onoff;
    printf("NlpHi is set %s\n\n", status[pNoired->NlpHiOn]);
    return;
}

void noisereduction_setNlpfullg(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->fullgNlpOn = onoff;
    printf("fullgNlp is set %s\n\n", status[pNoired->fullgNlpOn]);
    return;
}

void noisereduction_setNlpfullgHi(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->fullgNlpHiOn = onoff;
    printf("fullgNlpHi is set %s\n\n", status[pNoired->fullgNlpHiOn]);
    return;
}


void noisereduction_setNlpExtragain(NOIRED_PTR pNoired, bool onoff)
{
    const char *status[] = {"off", "on"};
    pNoired->extragainOn = onoff;
    printf("pNoired->extragain is set %s \n\n", status[pNoired->extragainOn]);
    return;
}

void noisereduction_setNlpSubusedEndHi(NOIRED_PTR pNoired, int value)
{
    pNoired->nlpSubusedEndHi = value;
    printf("pNoired->nlpSubusedEndHi is set to %d \n\n", pNoired->nlpSubusedEndHi);
    return;
}

void noisereduction_setNlpTransition(NOIRED_PTR pNoired, int value)
{
    pNoired->fullgNlpTransition = value;
    printf("pNoired->fullgNlpTransition is set to %d \n\n", pNoired->fullgNlpTransition);
    return;
}

void noisereduction_setDebug(NOIRED_PTR pNoired, int value)
{
    pNoired->debug = value;
    printf("Noired->debug is set to %d \n", pNoired->debug);
    return;
}

void noisereduction_status(NOIRED_PTR pNoired)
{
    const char *status[] = {"off", "on"};

    printf("\rStatus: Noired");
    printf("\r\n   nlp               :  %s", status[pNoired->NlpGainOn]);
    printf("\r\n   comNoise          :  %s", status[pNoired->ComNoiseOn]);
    printf("\r\n   comNoiseAmp       :  %f", pNoired->comNoiseAmp);
    printf("\r\n   noiseRed          :  %s", status[pNoired->NoiseRedOn]);
    printf("\r\n   nlpSubband        :  %s", status[pNoired->NLPSubOn]);
    printf("\r\n   nlpHi             :  %s", status[pNoired->NlpHiOn]);
    printf("\r\n   nlpSubusedEndHi   :  %d", pNoired->nlpSubusedEndHi);
    printf("\r\n   fullgNlp          :  %s", status[pNoired->fullgNlpOn]);
    printf("\r\n   fullgNlpHi        :  %s", status[pNoired->fullgNlpHiOn]);
    printf("\r\n   fullgNlpTransition- %d", pNoired->fullgNlpTransition);
    printf("\r\n   nlpExtragain      :  %s", status[pNoired->extragainOn]);
    printf("\r\n   debug             :  %d", pNoired->debug);
    printf("\r\n");
}

void noisereduction_ZerodelayStatus(NOIRED_PTR pNoired)
{
    const char *status[] = {"off", "on"};
    printf("\r\n   comNoise          :  %s", status[pNoired->ComNoiseOn]);
    printf("\r\n   comNoiseAmp       :  %f", pNoired->comNoiseAmp);
    printf("\r\n");
}





#ifdef UNITTEST
/***************************************************************************
* UNITTEST_NOISEREDUCTION
*    Desc.: run-trhough of noisereduction with 10 blocks equal random input on
*           louds and micinput.
*           The unittest checks mic.gain and noisevector output.
*           The noise output is not equal to matlab, but the nfgn has been
*           verified manually.
***************************************************************************/
#include "unittest.h"
#include "testdata/noisereductiontestdata.h"  /* testvectors defined here */

/* Private testdata */
static float noiredTestMicInput48[]     = {NOIRED48_MICTESTINPUT_DEFINE};
static float noiredTestEstmicInput48[]  = {NOIRED48_ESTMICTESTINPUT_DEFINE};
static float noiredTestGainOutput48[]   = {NOIRED48_GAINTESTOUTPUT_DEFINE};
static float noiredTestNoiseOutput48[]  = {NOIRED48_NOISETESTOUTPUT_DEFINE};

#define nWRITE2FILE

/***************************************************************************
* COMPARERESULTS
*    Desc.: Compares two float buffers. Due to compilator optimization, a
*           a small difference 'epsilon' is allowed when comparing with the
*           testvector
***************************************************************************/
static int compareResults(float *inp1, float *inp2, int len)
{
    int differ = 0;
    int i;
    float epsilon = 9.0e-7f;

    for( i = 0; i < len; i++ )
    {

        if( (fabsf(inp1[i] - inp2[i])) > epsilon )
        {
            differ = 1;
        }
    }
    return differ;
}

/***************************************************************************
* UNITTEST_NOISEREDUCTION
**************************************************************************/
void unittest_noisereduction()
{
    NOIRED_PTR pNoired; //pointer to noired struct
    DEREVERB_PTR pDereverb = NULL; //pointer to noired struct
    float *fftout;
    float excgain = 1.0f;
    float *estecho;
    float *gain;
    float *noise;
    float nlpconst = NLPC_MONO;
    bool  noisereductionOn = true;
    char const * testName = {"Noisereduction-48kHz"};
    int i, numFrames;
    int differGain; //criteria flag to be returned by function. 0: OK, 1: ERROR, 2:Nothing done, hence ERROR
    int differNoise; //criteria flag to be returned by function. 0: OK, 1: ERROR, 2:Nothing done, hence ERROR
    int size = 7680;
    FILTERBANK_USE_TYPE filterbankUse = FILTERBANK_EC;
    float *decay = NULL;

#ifdef WRITE2FILE
    float *outbuf;
    int cnt;
    FILE * fid;
    const char * Filename = "noiredtestoutput_host.dat";
#endif

    pNoired = noisereduction_create();
    if( pNoired == NULL )
    {
        fprintf(stderr, "unittest_noisereduction: memory allocation failed\n");
        unittest_assert(0);
        return;
    }

    unittest_context(testName);
    differGain = 2;
    differNoise = 2;

    noisereduction_init(pNoired, filterbankUse);
    pNoired->ncnt_boot = 0;   //to avoid boot behaviour, for comparison with matlab

    numFrames = size / FFTSIZE;

    //Allocating memory buffers
    estecho =(float *) calloc(FFTSIZE, sizeof(float));
    noise   =(float *) calloc(SUBUSED * 2, sizeof(float));
    gain    =(float *) malloc(SUBUSED * sizeof(float));
    fftout  =(float *) malloc(FFTSIZE * sizeof(float));
    for( i = 0; i < SUBUSED; i++ )
    {
        gain[i] = 1.0f;
    }
#ifdef WRITE2FILE
    outbuf =(float *) calloc(SUBUSED*numFrames,sizeof(float));
#endif
    for( i = 0; i < numFrames && differNoise != 1 && differGain != 1; i++ )
    {
        scratchmem_lock();
        memcpy(fftout, noiredTestMicInput48 + i * FFTSIZE, sizeof(float) * FFTSIZE);
        memcpy(estecho + 2 * SUBUSED_START, noiredTestEstmicInput48 + i * SUBUSED * 2, sizeof(float) * SUBUSED * 2);

        noisereduction_process(fftout,
                               excgain,
                               pNoired,
                               pDereverb,
                               estecho,
                               gain,
                               noise,
                               nlpconst,
                               noisereductionOn,
                               decay);

        differGain  = compareResults(gain,
                                     noiredTestGainOutput48 + i * SUBUSED,
                                     SUBUSED);
        differNoise = compareResults(noise,
                                     noiredTestNoiseOutput48 + i * SUBUSED * 2,
                                     SUBUSED*2);

#ifdef WRITE2FILE
        memcpy(outbuf+i*SUBUSED, gain, SUBUSED*sizeof(float));
#endif
        scratchmem_unlock();
    }


#ifdef WRITE2FILE
    printf("write noisereduction output to file \n");
    fid = (void*)fopen(Filename,"wb");
    cnt = fwrite(outbuf,sizeof(float),SUBUSED*numFrames,fid);
    fclose(fid);
    free(outbuf);
    (void)cnt;
#endif

    free(gain);
    free(noise);
    free(estecho);
    free(fftout);
    unittest_assert(differNoise == 0);
    unittest_assert(differGain == 0);

    noisereduction_destroy(pNoired);
}

#endif /* UNITTEST */
