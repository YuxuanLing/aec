#ifndef EC_PRIV_H
#define EC_PRIV_H

#include "analyse.h"
#include "lsexchndl.h"  //not really necessary because it is included in micexchndl.h
#include "dereverb.h"
#include "echocomp.h"
#include "micexchndl.h"
#include "noisereduction.h"
#include "synth.h"
#include "shell.h"

typedef struct EC {
    ANALYSE_PTR    pAnalyse;
    ANALYSE_PTR    pAnalyseSdram;
    ECHOCOMP_PTR   pEchocomp;
    MICEXCHNDL_PTR pMicexchndl;
    NOIRED_PTR     pNoired;
    DEREVERB_PTR   pDereverb;
    SYNTH_PTR      pSynth;
    SHELL_PTR      pShell;
	struct DELAY_ESTIMATION * delayEstimation;
    float *pAfir;
    float *pSynthIn;
    float *pGain;
    float *pNoise;
    float noiHPdelayline[2];      /* delayline for noisereducing highpassfilters (adsjustments.c: noiredHPfilt()) */
    bool noisereductionOn;
    int debug;
    KEYCLICKREMOVAL_PTR pKeyclickRemoval;
    bool keyprocessOn;
    bool noiredHPfiltOn;
    int channels;
} EC;

#endif /* EC_PRIV_H */
