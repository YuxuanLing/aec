#include "lsprocess_priv.h"
#include "lsprocess.h"
#include "analyse.h"
#include "lsexchndl.h"
#include "ec.h"
#include "audec_defs.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "analyse_priv.h"


LSPROCESS_PTR lsprocess_create(FILTERBANK_USE_TYPE filterbankUse, int channels)
{
    LSPROCESS_PTR pLs = (LSPROCESS_PTR) malloc(sizeof(LSPROCESS));
    int i;

    if( pLs == NULL )
    {
        fprintf(stderr, "lsprocess_create: Could not allocate pLs\n");
        return 0;
    }

    pLs->nChannels = channels;

    for (i = 0; i < channels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];

        channel->pLsexchndlSdram = lsexchndl_create();
        channel->pAnalyseSdram = analyse_create(filterbankUse);
    }
    return pLs;
}


void lsprocess_destroy(LSPROCESS_PTR pLs, FILTERBANK_USE_TYPE filterbankUse)
{
    int i;

    for (i = 0; i < pLs->nChannels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];

        lsexchndl_destroy(channel->pLsexchndlSdram);
        analyse_destroy(channel->pAnalyseSdram, filterbankUse);
    }
    free(pLs);
}

void lsprocess_init(LSPROCESS_PTR pLs, FILTERBANK_USE_TYPE filterbankUse)
{
    int i;

    for (i = 0; i < pLs->nChannels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];

        lsexchndl_init(channel->pLsexchndlSdram);
        analyse_init(channel->pAnalyseSdram, filterbankUse);
        channel->lsinput = NULL;
    }

    pLs->loudsGain = 1.0f;
    pLs->lsGainAdjustOn = true;
    pLs->debug = 0;
}

void lsprocess_load(LSPROCESS_PTR pLs, float * lsBuf)
{
    int i;

    for (i = 0; i < pLs->nChannels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];

        channel->lsinput = &lsBuf[i * FRAMESIZE];
        analyse_loadInbuf(channel->pAnalyseSdram);
    }
}

void lsprocess_process(LSPROCESS_PTR pLs, float *afir)
{
    int i;

    for (i = 0; i < pLs->nChannels; i++)
    {
        LSPROCESS_CHANNEL *channel = &pLs->channel[i];

        channel->pAnalyse = analyse_loadAnalyse(channel->pAnalyseSdram);  /* uses memcpy */

        if(pLs->debug == 12)
        {
            ec_debug12print(channel->lsinput, pLs->loudsGain);  /* debug print need to be here because analyse prefilter modifies lsinput */
        }

        analyse_process(channel->pAnalyse, channel->lsinput, afir);  //runs analysefilter on louds buffer

        if( pLs->lsGainAdjustOn)
        {
            analyse_gainAdjust(channel->pAnalyse, pLs->loudsGain);
        }
        else
        {
            pLs->loudsGain = 1.0f;
        }

        analyse_flushAnalyse(channel->pAnalyseSdram, channel->pAnalyse);  /* uses memcpy */
    }
}


void lsprocess_setloudsGain(LSPROCESS_PTR pLs, float value)
{
    pLs->loudsGain = value;
}

void lsprocess_setlsGainAdjust(LSPROCESS_PTR pLs, bool onoff)
{
    const char *status[] = {"off", "on"};
    pLs->lsGainAdjustOn = onoff;
    printf("lsGainAdjust is set %s\n\n", (char*)status[pLs->lsGainAdjustOn]);
}

void lsprocess_setDebug(LSPROCESS_PTR pLs, int value)
{
    pLs->debug = value;
    printf("pLs->debug is set to %d \n", pLs->debug);
}

void lsprocess_status(LSPROCESS_PTR pLs)
{
    const char *status[] = {"off", "on"};
    printf("\rStatus: Ec");
    printf("\r\n   lsGainAdjust     - %s",     (char*)status[pLs->lsGainAdjustOn]);
    printf("\r\n   debug            - %d",     pLs->debug);
    printf("\r\n");
}
