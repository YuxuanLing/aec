#ifndef DEREVERB_H
#define DEREVERB_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Gisle Enstad (GEN), Tandberg Telecom AS.
 *
 **************************************************************************/

#include "audtypes.h"
#include <stdbool.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DEREVERB * DEREVERB_PTR;

typedef struct DEREVERB_PROCESS
{
    float * delaylineBufN;
    float * delaylineBufNmin1;
    float * delaylineBufEnd;
    int delaylineIndex;
}DEREVERB_PROCESS;

DEREVERB_PTR dereverb_create(void);

void dereverb_destroy(DEREVERB_PTR pDereverb);

void dereverb_init(DEREVERB_PTR pDereverb);

void dereverb_setDebug(DEREVERB_PTR pDereverb, int value);

void dereverb_process(DEREVERB_PROCESS * pDereverb,
                      float *  fftout,
                      float *  decayInGainOut,
                      float mingain);

void dereverb_loadDereverb(DEREVERB_PTR pDereverb, DEREVERB_PROCESS * pDereverbProcess);

void dereverb_flushDereverb(DEREVERB_PTR pDereverbSdram, DEREVERB_PROCESS * pDereverbProcess);

#ifdef __cplusplus
}
#endif


#endif /* DEREVERB_H */


