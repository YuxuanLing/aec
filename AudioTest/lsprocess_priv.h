#ifndef LSPROCESS_PRIV_H
#define LSPROCESS_PRIV_H

#include "analyse.h"
#include "lsexchndl.h"
#include "audec_defs.h"
#include <stdbool.h>

typedef struct LSPROCESS_CHANNEL {
    float          *lsinput;
    ANALYSE_PTR    pAnalyse;
    ANALYSE_PTR    pAnalyseSdram;
    LSEXCHNDL_PTR  pLsexchndl;
    LSEXCHNDL_PTR  pLsexchndlSdram;
} LSPROCESS_CHANNEL;

typedef struct LSPROCESS {
    LSPROCESS_CHANNEL channel[MAX_NUM_CHANNELS];
    int nChannels;

    float loudsGain;
    bool  lsGainAdjustOn;
#ifdef PLATFORM_SNOOPY
    short filtlen[SUBUSED];
#endif
    int debug;
} LSPROCESS;


#endif
