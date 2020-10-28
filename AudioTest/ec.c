#include "ec_priv.h"
#include "ec.h"
#include "adjustments.h"
#include "audec_defs.h"
#include "fft.h"
#include "mathfun.h"
#include "deesser.h"
#include "shell.h"
#include "delayEstimation.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "mm_malloc.h"

#include "analyse_priv.h"
#include "lsprocess_priv.h"
#include "keyclickremoval.h"

/* EC_CREATE */
EC_PTR ec_create(FILTERBANK_USE_TYPE filterbankUse, int channels)
{
    EC_PTR pEc;

    if (channels < 1 || channels > MAX_NUM_CHANNELS)
    {
        /* Illegal number of channels */
        abort();
    }

    pEc =(EC_PTR) _mm_malloc(sizeof(EC), 64);
    if( pEc == NULL )
    {
        /* Could not allocate ec struct */
        abort();
    }
    pEc->pAnalyseSdram  = analyse_create(filterbankUse);
    pEc->pEchocomp   = echocomp_create(channels);
    pEc->pNoired     = noisereduction_create();
    pEc->pDereverb     = dereverb_create();
    pEc->pMicexchndl = micexchndl_create();
    pEc->pSynth      = synth_create();
    pEc->pShell      = shell_create();
    pEc->pGain       = (float *) malloc(SUBUSED * sizeof(float));
	pEc->delayEstimation = delayEstimation_create();

    if( pEc->pGain == NULL )
    {
        /* Could not allocate pEc->pGain */
        abort();
    }

    if(pEc->pAnalyseSdram   == NULL ||
       pEc->pEchocomp       == NULL ||
       pEc->pNoired         == NULL ||
       pEc->pMicexchndl     == NULL ||
       pEc->pSynth          == NULL ||
       pEc->pShell          == NULL ||
       pEc->pGain           == NULL ||
	   pEc->delayEstimation == NULL)
    {
        /* failed allocation */
        abort();
    }

    pEc->pKeyclickRemoval = keyclickremoval_create();
    pEc->channels = channels;

    return pEc;
}


/* EC_DESTROY */
void ec_destroy(EC_PTR pEc, FILTERBANK_USE_TYPE filterbankUse)
{
    analyse_destroy(pEc->pAnalyseSdram, filterbankUse);
    echocomp_destroy(pEc->pEchocomp);
    noisereduction_destroy(pEc->pNoired);
    dereverb_destroy(pEc->pDereverb);
    micexchndl_destroy(pEc->pMicexchndl);
    synth_destroy(pEc->pSynth);
    shell_destroy(pEc->pShell);
	delayEstimation_destroy(pEc->delayEstimation);
    free(pEc->pGain);
    keyclickremoval_destroy(pEc->pKeyclickRemoval);
    _mm_free(pEc);
}


/* EC_INIT */
void ec_init(EC_PTR pEc, FILTERBANK_USE_TYPE filterbankUse, LSPROCESS_PTR pLs)
{
    int i;
    (void)pLs;  // no longer used

    analyse_init(pEc->pAnalyseSdram, filterbankUse);
    echocomp_init(pEc->pEchocomp);
    noisereduction_init(pEc->pNoired, filterbankUse);
    dereverb_init(pEc->pDereverb);
    micexchndl_init(pEc->pMicexchndl, pEc->channels);
    synth_init(pEc->pSynth);
    shell_init(pEc->pShell, pEc->channels);
	delayEstimation_init(pEc->delayEstimation);
    fft_init();

    pEc->noiHPdelayline[0] = 0.0f;
    pEc->noiHPdelayline[1] = 0.0f;
    //pEc->pAfir = NULL;

    for( i = 0; i<SUBUSED; i++ )
    {
        pEc->pGain[i] = 1.0f;
    }
    pEc->debug = 0;
    pEc->noisereductionOn = true;
    pEc->pAfir = NULL;

#ifdef PLATFORM_SNOOPY
    pEc->keyprocessOn = true;
#else
    pEc->keyprocessOn = false;
#endif

    pEc->noiredHPfiltOn = true;
    keyclickremoval_init(pEc->pKeyclickRemoval);
}

/* EC_LOAD
 * Starts dmax transfers on analyse and lsexchndl data.
 * Requires a later call to memxfer_wait with buffer-pointers
 *  to ensure successful transfers.
 */
void ec_load(EC_PTR pEc)
{
    analyse_loadInbuf(pEc->pAnalyseSdram);
    pEc->pAfir = analyse_loadAfir(pEc->pAnalyseSdram);
}


/* EC_PROCESS
 * Filterbank and Echo cancellator
 */
void ec_process(EC_PTR pEc,
                LSPROCESS_PTR pLs,
                float *micBuf,
                float *outBuf,
                bool runLsprocess)
{
    /* scratch variables */
    float mictemp[FRAMESIZE];
    LSEXCHNDL_PTR lsexchndl[MAX_NUM_CHANNELS];
    float *micinput;
    DEESSER_PTR pDeEsser;
    float* decay = NULL; //FIXME! Where is the original?

    float nlpconst = NLPC_MONO;
    float *gain;
    __declspec(align(64)) COMPLEX32 noise[SUBUSED], echoest[FFTSIZE/2], synthIn[FFTSIZE/2];
    MICEXCHNDL_PTR  pMicexchndl;
    NOIRED_PTR      pNoired;
    SYNTH_PTR       pSynth;
    SHELL_PTR       pShell;
    int m, n, ch;
    float loudsGain_tmp = pLs->loudsGain; /* only used for lsexchndl when runLsprocess is true */
    float shellgain;
	int  movedDelay = 0;

    /* MEMXFERMUSTDIE */

    if( pEc->noisereductionOn && pEc->noiredHPfiltOn )
    {
        micinput   = mictemp; //output from noiredHPfilt
        noiredHPfilt(&(pEc->noiHPdelayline[0]), micBuf, micinput);
    }
    else
    {
        micinput = micBuf;
    }

    shell_input(pEc->pShell, micinput);  //make copy of micinput

    pEc->pAnalyse = analyse_loadAnalyse(pEc->pAnalyseSdram);    /* uses memcpy */
    analyse_process(pEc->pAnalyse, micinput, pEc->pAfir);
    analyse_flushAnalyse(pEc->pAnalyseSdram, pEc->pAnalyse);      /* uses memcpy */

    pDeEsser = deesser_getDeEsser();
    if (pDeEsser->deEsserOn)
    {
        deesser_process(pEc->pAnalyse->fftout);
    }

    if(runLsprocess)
    {
        lsprocess_process(pLs, pEc->pAfir);
    }

    movedDelay = delayEstimation_estimateDelay(pEc->delayEstimation, pEc->pAnalyseSdram->fftout, pLs->channel[0].pAnalyseSdram->fftout);
    echocomp_moveFilters(pEc->pEchocomp, pEc->delayEstimation, movedDelay); //move taps in adapt and fixed filters according to delay

    if (pEc->channels == 1)
    {
        echocomp_process(pEc->pEchocomp,
                         pLs->channel[0].pAnalyseSdram->fftout,
                         echoest + SUBUSED_START,
                         pEc->pAnalyseSdram->fftout + SUBUSED_START,
                         micexchndl_getExcgain(pEc->pMicexchndl),
                         runLsprocess);
    }
    else /* channels == 2 */
    {
        echocomp_processStereo(pEc->pEchocomp,
                               pLs->channel[0].pAnalyseSdram->fftout,
                               pLs->channel[1].pAnalyseSdram->fftout,
                               echoest + SUBUSED_START,
                               pEc->pAnalyseSdram->fftout + SUBUSED_START,
                               micexchndl_getExcgain(pEc->pMicexchndl),
                               runLsprocess);
    }

    for (ch = 0; ch < pEc->channels; ch++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[ch];
        channel->pLsexchndl  = lsexchndl_load(channel->pLsexchndlSdram);
    }
    pMicexchndl     = micexchndl_load(pEc->pMicexchndl);
    pNoired         = noisereduction_loadNoired(pEc->pNoired);
    gain            = noisereduction_loadGain(pEc->pGain);
    pShell          = shell_load(pEc->pShell);
    pSynth          = synth_loadSynth(pEc->pSynth); /* uses memcpy */
    synth_loadfftst(pSynth);
    synth_loadSfir(pSynth);

    for (ch = 0; ch < pEc->channels; ch++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[ch];
        if (runLsprocess)
        {
            lsexchndl_process(channel->pLsexchndl, (float *) channel->pAnalyseSdram->fftout, loudsGain_tmp);
        }
        lsexchndl[ch] = channel->pLsexchndl;
    }
    micexchndl_process(lsexchndl, pMicexchndl, (float *) pEc->pAnalyseSdram->fftout, (float *) echoest);
    micexchndl_flush(pEc->pMicexchndl, pMicexchndl);

    if (runLsprocess)
    {
        for (ch = 0; ch < pEc->channels; ch++)
        {
            LSPROCESS_CHANNEL *channel = &pLs->channel[ch];
            lsexchndl_flush(channel->pLsexchndlSdram, channel->pLsexchndl);
        }
    }

    noisereduction_process(pEc->pAnalyseSdram->fftout,
                           micexchndl_getExcgain(pEc->pMicexchndl),
                           pNoired,
                           pEc->pDereverb,
                           echoest,
                           gain,
                           noise,
                           nlpconst,
                           pEc->noisereductionOn,
                           decay); // FIXME: decay, where does it come from?
    
    shellgain = shell_process(pShell, pNoired, pLs);
    noisereduction_flushNoired(pEc->pNoired, pNoired);

    // TODO: Move these loops into a separate function to reduce the size of this one
    for( m = 0; m < SUBUSED_START; m++ )
    {
        synthIn[m].re = 0.0f;
        synthIn[m].im = 0.0f;
    }

    for( m = SUBUSED_START, n=0; m < SUBUSED_START + SUBUSED; m++, n++ )
    {
        synthIn[m].re = (pEc->pAnalyseSdram->fftout[m].re - echoest[m].re) * gain[n]
                          + noise[m - SUBUSED_START].re
                          - shellgain * pEc->pAnalyseSdram->fftout[m].re;
        synthIn[m].im = (pEc->pAnalyseSdram->fftout[m].im - echoest[m].im) * gain[n]
                          + noise[m - SUBUSED_START].im
                          - shellgain * pEc->pAnalyseSdram->fftout[m].im;
    }

    for( m = SUBUSED_START + SUBUSED; m < FFTSIZE/2; m++ )
    {
        synthIn[m].re = 0.0f;
        synthIn[m].im = 0.0f;
    }

    synth_process(pSynth, synthIn, outBuf);
    synth_flushSynth(pEc->pSynth, pSynth); /* uses memcpy */

    shell_output(pShell, outBuf);
    shell_flush(pEc->pShell, pShell);

    if( pEc->keyprocessOn )   /* ON for snoopy only */
    {
        keyclickremoval_limiter_process(pEc->pKeyclickRemoval, outBuf);
    }

}


/* ec_getAerl */
float ec_getAerl(EC_PTR pEc, LSPROCESS_PTR pLs,  float miclevel)
{
    return echocomp_getAerl(pEc->pEchocomp, miclevel, pLs->loudsGain);
}

void ec_updateAfirIndex(EC_PTR pEc, int index)
{
    analyse_updateInbufStart(pEc->pAnalyseSdram, index);
}

int ec_getAfirIndex(EC_PTR pEc)
{
    return(analyse_getInbufStart(pEc->pAnalyseSdram));
}

void ec_keyclickEvent(EC_PTR pEc, int *inbuf, int sampleindex)
{
    keyclickremoval_verifyEvent(pEc->pKeyclickRemoval, inbuf, sampleindex);
}

/* enable or disable noiredHPfilter and subband-noisereduction  */
void ec_setNoisereduction(EC_PTR pEc, bool onoff)
{
    pEc->noisereductionOn = onoff;
    return;
}


/*****************************************************
 * Test functions
 ****************************************************/

/* ec_setcalcNewDelta */
void ec_setCalcNewDelta(EC_PTR pEc, bool onoff)
{
    echocomp_setCalcNewDelta(pEc->pEchocomp, onoff);
}

/* ec_setMinDelta */
void ec_setMinDelta(EC_PTR pEc, float value)
{
    echocomp_setMinDelta(pEc->pEchocomp, value);
}

/* ec_setMaxDelta */
void ec_setMaxDelta(EC_PTR pEc, float value)
{
    echocomp_setMaxDelta(pEc->pEchocomp, value);
}

/* ec_setFilter */
void ec_setFilter(EC_PTR pEc, int value)
{
    echocomp_setFilter(pEc->pEchocomp, value);
}
/* ec_setFilterZero */
void ec_setFilterZero(EC_PTR pEc, int value)
{
    echocomp_setZero(pEc->pEchocomp, value);
}

/* ec_setNlp */
void ec_setNlp(EC_PTR pEc, bool onoff)
{
    noisereduction_setNlp(pEc->pNoired, onoff);
}

/* ec_setNlp */
void ec_setDereverb(EC_PTR pEc, bool onoff)
{
    noisereduction_setDereverb(pEc->pNoired, onoff);
}

/* ec_setComNoise */
void ec_setComNoise(EC_PTR pEc, bool onoff)
{
    noisereduction_setComNoise(pEc->pNoired, onoff);
}

/* ec_setNoiseRed
 * enable or disable subband-noisereduction  */
void ec_setSubbandNoiseRed(EC_PTR pEc, bool onoff)
{
    noisereduction_setNoiseRed(pEc->pNoired, onoff);
}

/* ec_setNlpSub */
void ec_setNlpSub(EC_PTR pEc, bool onoff)
{
    noisereduction_setNlpSub(pEc->pNoired, onoff);
}

/* ec_setNlpHi */
void ec_setNlpHi(EC_PTR pEc, bool onoff)
{
    noisereduction_setNlpHi(pEc->pNoired, onoff);
}

/* ec_setNlpfullg */
void ec_setNlpfullg(EC_PTR pEc, bool onoff)
{
    noisereduction_setNlpfullg(pEc->pNoired, onoff);
}

/* ec_setNlpfullgHi */
void ec_setNlpfullgHi(EC_PTR pEc, bool onoff)
{
    noisereduction_setNlpfullgHi(pEc->pNoired, onoff);
}

/* ec_setNlpExtragain */
void ec_setNlpExtragain(EC_PTR pEc, bool onoff)
{
    noisereduction_setNlpExtragain(pEc->pNoired, onoff);
}

/* ec_setNlpSubusedEndHi */
void ec_setNlpSubusedEndHi(EC_PTR pEc, float value)
{
    noisereduction_setNlpSubusedEndHi(pEc->pNoired, value);
}

/* ec_setNlpTransition */
void ec_setNlpTransition(EC_PTR pEc, float value)
{
    noisereduction_setNlpTransition(pEc->pNoired, value);
}

void ec_setKeyclickprocess(EC_PTR pEc, bool onoff)
{
    const char *status[] = {"off", "on"};
    pEc->keyprocessOn = onoff;
    printf("pEc->keyprocessOn is set %s \n", status[pEc->keyprocessOn]);

    keyclickremoval_setKeyclickprocess(pEc->pKeyclickRemoval, onoff);
}


void ec_setKeyclickmaxlimitcount(EC_PTR pEc, int value)
{
    keyclickremoval_setmaxlimitcount(pEc->pKeyclickRemoval, value);
}

void ec_setKeyclickgi(EC_PTR pEc, float value)
{
    keyclickremoval_setgi(pEc->pKeyclickRemoval, value);
}

void ec_setKeyclickgd(EC_PTR pEc, float value)
{
    keyclickremoval_setgd(pEc->pKeyclickRemoval, value);
}


void ec_setShellOn(EC_PTR pEc, bool onoff)
{
    shell_setShellOn(pEc->pShell, onoff);
}

void ec_setShellLsLevAdjust(EC_PTR pEc, float value)
{
    shell_setlslevAdjust(pEc->pShell, value);
}

void ec_setShellnoisgnAdjust(EC_PTR pEc, float value)
{
    shell_setnoisegnAdjust(pEc->pShell, value);
}

void ec_setShellnlpgAdjust(EC_PTR pEc, float value)
{
    shell_setnlpgAdjust(pEc->pShell, value);
}

void ec_setShellvarIncspeed(EC_PTR pEc, float value)
{
    shell_setvarIncSpeed(pEc->pShell, value);
}

void ec_setShellvarDecspeed(EC_PTR pEc, float value)
{
    shell_setvarDecSpeed(pEc->pShell, value);
}


/* ec_setFbPrototypefilt */
void ec_setFbPrototypefilt(EC_PTR pEc, LSPROCESS_PTR pLs, int type, FILTERBANK_USE_TYPE filterbankUse)
{
    int i;

    /* Call init functions to clear state */
    analyse_init(pEc->pAnalyseSdram, filterbankUse);
    for (i = 0; i < pEc->channels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];
        analyse_init(channel->pAnalyseSdram, filterbankUse);
    }
    synth_init(pEc->pSynth);

    analyse_setFbPrototypefilt(pEc->pAnalyseSdram, type);
    for (i = 0; i < pEc->channels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];
        analyse_setFbPrototypefilt(channel->pAnalyseSdram, type);
    }
    synth_setFbPrototypefilt(pEc->pSynth, type);
}

/* ec_setMicexchndl */
void ec_setMicexchndl(EC_PTR pEc, bool onoff)
{
    micexchndl_set(pEc->pMicexchndl, onoff);
}

void ec_setComNoiseAmp(EC_PTR pEc, float value)
{
    noisereduction_setComNoiseAmp(pEc->pNoired, value);
}

void ec_setNoiseredHpFilt(EC_PTR pEc, bool onoff)
{
    pEc->noiredHPfiltOn = onoff;
    return;
}

void ec_setStereoDelta(EC_PTR pEc, float value)
{
    echocomp_setStereoDelta(pEc->pEchocomp, value);
}

/* ec_setDebug */
void ec_setDebug(EC_PTR pEc, int value)
{
    pEc->debug = value;
    printf("pEc->debug is set to %d \n", pEc->debug);
    noisereduction_setDebug(pEc->pNoired, value);
    dereverb_setDebug(pEc->pDereverb, value);
    echocomp_setDebug(pEc->pEchocomp, value);
    shell_setDebug(pEc->pShell, value);
    return;
}

void ec_printNumChannels(EC_PTR pEc)
{
    printf("ec_channels = %d\n", pEc->channels);
}

/* ec_status */
void ec_status(EC_PTR pEc)
{
    const char *status[] = {"off", "on"};

    printf("\rStatus: Ec");
    printf("\r\n   noiseredHpFilt   - %s",     status[pEc->noiredHPfiltOn]);
    printf("\r\n   channels         - %d",     pEc->channels);
    printf("\r\n   debug            - %d",     pEc->debug);
    printf("\r\n");
    micexchndl_statusMicexchndl(pEc->pMicexchndl);
    noisereduction_status(pEc->pNoired);
    echocomp_status(pEc->pEchocomp);
    synth_status(pEc->pSynth);
    analyse_status(pEc->pAnalyseSdram);
    keyclickremoval_status(pEc->pKeyclickRemoval);
    shell_status(pEc->pShell);
}

/* Find and print max sample value of loudspeaker signal. Before and after volume adjustment */
void ec_debug12print(float * pLsBuf, float loudsGain)
{
    static float maxSample = 0.0f;
    static int cnt = 0;
    int i;
    float loudsGainInv = 1.0f;

    if (loudsGain > 1.0e-7f)
    {
        loudsGainInv = fastdiv(1.0f, loudsGain);     // loudsGainInv = 1/loudsGain
    }
    else
    {
        loudsGainInv = 1.0f;
        //printf("analyse_gainAdjust: loudsGain < EPSILON \n");
    }
#ifdef __INTEL_COMPILER
#pragma novector
#endif
    for( i=0; i<AUD_INT_BUFSIZE;i++ )
    {
        float tmp = fabsf(pLsBuf[i]);
        if( tmp > maxSample )
        {
            maxSample = tmp;
        }
    }

    if( ++cnt == 20 )
    {
        printf("%11.5e, %11.5e \n", maxSample * loudsGainInv, maxSample);
        maxSample = 0.0f;
        cnt = 0;
    }
}

void ec_setWindowsDelay(EC_PTR pEc, float value)
{
	delayEstimation_setWindowsDelay(pEc->delayEstimation, value);
}

void ec_setDelayEstimationOnOff(EC_PTR pEc, int onoff)
{
	delayEstimation_setDelayingOn(pEc->delayEstimation, onoff);
}

#ifdef UNITTEST
/******************************************************************************
 * UNITTEST
 *   Description:
 *       Run-through of echocanceller.
 *       The test input and output signals to this unittest are read/written
 *       to/from file.
 *****************************************************************************/
#include "unittest.h"


/* compareResults ******************************************************* */
/* Auth.:                                                                 */
/* Desc.: Compares two float buffers. Due to compilator optimization, a   */
/*        a small difference 'epsilon' is allowed when comparing with the */
/*        testvector                                                     */
/* ********************************************************************** */
static int compareResults(float *inp1, float *inp2, int len)
{
    int differ = 0;
    int i;
    //const float epsilon = 1.5e-8f;
    const float epsilon = 1.0e-3f;

    for(i = 0; i < len; i++)
    {
        if(fabsf(inp1[i] - inp2[i]) > epsilon)
        {
            differ = 1;
            //printf("%d \t %23.18f \t %23.18f \t %23.18f \n ",i,inp1[i],inp2[i], fabsf(inp1[i] - inp2[i]));
            break;
        }
    }

    return differ;
}

// TODO: Use portable clock wrappers
#ifdef CHIP_674X
#include <std.h>
#include <clk.h>
#include <sts.h>
#include "timing.h"
extern STS_Obj sts_ec_unittest;
#endif

void run_unittest_ec(const char fileMic[], const char fileNetworkIn[],
                     const char fileNetworkOutReference[], const char fileNetworkOut[],
                     int numFrames)
{
    enum { channels = 1 };
    
    extern int IRAM;

    EC_PTR pEc;
    LSPROCESS_PTR pLs;
    FILTERBANK_USE_TYPE filterbankUse = FILTERBANK_EC;

    FILE *fidMic;
    FILE *fidNetworkIn;
    FILE *fidNetworkOutReference;
    FILE *fidNetworkOut;

    float *micBuf;
    float *lsBuf;
    float *outBufReference;
    float *outBuf;

    int i;
    int cnt = 0;
    int differ = 0;
    int minFrames = 10;

    #ifdef CHIP_674X
    float t;
    #endif

    unittest_context("echo canceller");

    pLs = lsprocess_create(filterbankUse, channels);
    lsprocess_init(pLs, filterbankUse);

    pEc = ec_create(filterbankUse, channels);
    ec_init(pEc, filterbankUse, pLs);

    micBuf = (float*) dsp_memAllocSeg(IRAM, FRAMESIZE*sizeof(float), 128); 
    lsBuf = (float*) dsp_memAllocSeg(IRAM, FRAMESIZE*sizeof(float), 128);
    outBufReference = (float*) dsp_memAllocSeg(IRAM, FRAMESIZE*sizeof(float), 128);
    outBuf = (float*) dsp_memAllocSeg(IRAM, FRAMESIZE*sizeof(float), 128);

    memset(micBuf, 0, FRAMESIZE*sizeof(float));
    memset(lsBuf, 0, FRAMESIZE*sizeof(float));
    memset(outBufReference, 0, FRAMESIZE*sizeof(float));
    memset(outBuf, 0, FRAMESIZE*sizeof(float));

    if(micBuf == NULL || lsBuf == NULL || outBuf == NULL || outBufReference == NULL)
    {
        printf("ERROR: Could not allocate test buffers\n");
        exit(2);
    }

    fidMic = fopen(fileMic, "rb");
    if(fidMic == NULL)
    {
        printf("ERROR: Could not open file %s\n", fileMic);
        exit(2);
    }

    fidNetworkIn = fopen(fileNetworkIn, "rb");
    if(fidNetworkIn == NULL)
    {
        printf("ERROR: Could not open file %s\n", fileNetworkIn);
        exit(2);
    }

    fidNetworkOutReference = fopen(fileNetworkOutReference, "rb");
    if(fidNetworkOutReference == NULL)
    {
        printf("ERROR: Could not open file %s\n", fileNetworkOutReference);
        exit(2);
    }

    fidNetworkOut = fopen(fileNetworkOut, "wb");
    if(fidNetworkOut == NULL)
    {
        printf("ERROR: Could not open file %s\n", fileNetworkOut);
        exit(2);
    }

    #ifdef CHIP_674X
    STS_reset(&sts_ec_unittest);
    #endif

    printf("START ec_loop\n");
    for(i = 0; i < numFrames; i++)
    {
        if(i % 10 == 0)
        {
            printf("Reading ec unittest input signals, frame %d\n", i);
        }
        cnt = fread(micBuf, sizeof(float), FRAMESIZE, fidMic);
        cnt = fread(lsBuf, sizeof(float), FRAMESIZE, fidNetworkIn);
        cnt = fread(outBufReference, sizeof(float), FRAMESIZE, fidNetworkOutReference);

        scratchmem_lock();

        lsprocess_load(pLs, lsBuf);
        ec_load(pEc);

        #ifdef CHIP_674X
        STS_set(&sts_ec_unittest, CLK_gethtime());
        tic();
        #endif

        ec_process(pEc, pLs, micBuf, outBuf, true);

        #ifdef CHIP_674X
        t = toc();
        STS_delta(&sts_ec_unittest, CLK_gethtime());
        printf("one frame: %.2e\n", t);
        #endif

        scratchmem_unlock();

        cnt = fwrite(outBuf, sizeof(float), FRAMESIZE, fidNetworkOut);

        differ = compareResults(outBuf, outBufReference, FRAMESIZE);

        if(differ && i > minFrames)
            break;
    }
    printf("END ec_loop after %d frames\n", i);

    #ifdef CHIP_674X
    printf("Timing statistics:\n");
    printf("    num: %d\n", sts_ec_unittest.num);
    printf("    acc: %d\n", sts_ec_unittest.acc);
    printf("    max: %d\n", sts_ec_unittest.max);
    printf("    avg: %.2e\n", ((float)sts_ec_unittest.acc) / sts_ec_unittest.num);
    #endif

    // TODO: Use combined assertion/reporting function from test framework 
    if(differ != 0)
    {
        printf("Computed network output differs from reference data.\n");
        unittest_assert(differ == 0);
    }

    fclose(fidMic);
    fclose(fidNetworkIn);
    fclose(fidNetworkOutReference);
    fclose(fidNetworkOut);

    dsp_memFreeSeg(IRAM, micBuf, FRAMESIZE*sizeof(float));
    dsp_memFreeSeg(IRAM, lsBuf, FRAMESIZE*sizeof(float));
    dsp_memFreeSeg(IRAM, outBufReference, FRAMESIZE*sizeof(float));
    dsp_memFreeSeg(IRAM, outBuf, FRAMESIZE*sizeof(float));

    ec_destroy(pEc, filterbankUse);

    printf("END ec_unittest\n");

    (void)cnt;
}

void unittest_ec(void)
{
    // TODO: Place test data somewhere central
    const char fileMic[]        = "/home/jps/src/testdata/ecMic48.dat";
    const char fileNetworkIn[]  = "/home/jps/src/testdata/ecNetworkIn48.dat";
    const char fileNetworkOutReference[] = "/home/jps/src/testdata/ecNetworkOutReference48.dat";
    const char fileNetworkOut[] = "/home/jps/src/testdata/ecNetworkOutComputed48.dat";
    int numFrames = 1145;
    run_unittest_ec(fileMic, fileNetworkIn,
                    fileNetworkOutReference, fileNetworkOut,
                    numFrames);
}

#endif /* UNITTEST */
