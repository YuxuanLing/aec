/***************************************************************************
*                              A U D I O
*--------------------------------------------------------------------------
*                  (C) Copyright Tandberg Telecom AS 2009
*==========================================================================
*
* Author        : Gisle Enstad (GEN), Tandberg Telecom AS.
* Co-author     : -
*
* Switches      : <COMPILER SWITCH HERE>
*         <add description here>
*
* Description   : Dereverberation by spectral subtraction
*
* Note    : <notes here>
*
* Documentation : <Xxxxxx-Document.doc>
*
**************************************************************************/
#include "audec_defs.h"
#include "dereverb_priv.h"
#include "dereverb.h"
#include "mathfun.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


DEREVERB_PTR dereverb_create(void)
{
    DEREVERB_PTR pDereverb;

    pDereverb = (DEREVERB_PTR) malloc(sizeof(DEREVERB));

    if( pDereverb == NULL )
    {
        fprintf(stderr, "dereverb_create: Could not allocate dereverbs buffer\n");
        return 0;
    }

    return pDereverb;
}

void dereverb_destroy(DEREVERB_PTR pDereverb)
{
    free(pDereverb);
}

void dereverb_init(DEREVERB_PTR pDereverb)
{
    int i, j;

    for(i = 0; i < DEREVERB_DELAYLINE; i++) {
        for(j = 0; j < SUBUSED; j++) {
            pDereverb->siglevDelayline[i][j] = 0.0f;
        }
    }

    pDereverb->delaylineIndex = 0;
    pDereverb->debug = 0;
}


void dereverb_setDebug(DEREVERB_PTR pDereverb, int value)
{
    pDereverb->debug = value;
    printf("pDereverb->debug is set to %d \n", pDereverb->debug);
    return;
}

void dereverb_loadDereverb(DEREVERB_PTR pDereverb, DEREVERB_PROCESS * pDereverbProcess)
{
    int Nmin1;
    int Nend;

    pDereverbProcess->delaylineIndex = pDereverb->delaylineIndex;

    if(pDereverbProcess->delaylineIndex == (DEREVERB_DELAYLINE-1)){
        Nmin1 = pDereverbProcess->delaylineIndex - 1;
        Nend = 0;
    }
    else if(pDereverbProcess->delaylineIndex == 0){
        Nmin1 = (DEREVERB_DELAYLINE-1);
        Nend = pDereverbProcess->delaylineIndex + 1;
    }
    else{
        Nmin1 = pDereverbProcess->delaylineIndex - 1;
        Nend = pDereverbProcess->delaylineIndex + 1;
    }
    
    /* MEMXFERMUSTDIE */
    pDereverbProcess->delaylineBufN     = &pDereverb->siglevDelayline[pDereverbProcess->delaylineIndex][0];
    pDereverbProcess->delaylineBufNmin1 = &pDereverb->siglevDelayline[Nmin1][0];
    pDereverbProcess->delaylineBufEnd   = &pDereverb->siglevDelayline[Nend ][0];
}

void dereverb_flushDereverb(DEREVERB_PTR pDereverbSdram, DEREVERB_PROCESS * pDereverbProcess)
{
    /* MEMXFERMUSTDIE */
    (void) pDereverbProcess;

    pDereverbSdram->delaylineIndex++;
    if(pDereverbSdram->delaylineIndex == DEREVERB_DELAYLINE) {
        pDereverbSdram->delaylineIndex = 0;
    }

}

void dereverb_process(DEREVERB_PROCESS * pDereverb,
                      float *  fftout,
                      float *  decayInGainOut,
                      float mingain)
{
    int n, m;
    float alpha = 0.9f;
    float maxgain = 1.0f;
    float tmp;
    float RSR;
    float *  delaylineBufN = pDereverb->delaylineBufN;
    float *  delaylineBufNmin1 = pDereverb->delaylineBufNmin1;
    float *  delaylineBufEnd = pDereverb->delaylineBufEnd;

    /* Level delayline index */
    n = SUBUSED_START * 2;

    /* Level estimate delayline */
    for( m = 0; m < SUBUSED; m++, n += 2 )
    {
        tmp = fftout[n] * fftout[n] + fftout[n + 1] * fftout[n + 1];
        
        delaylineBufN[m] = tmp + alpha*(delaylineBufNmin1[m] - tmp);
    }

    /* Signal to reverberation level estimate, SRR */

    for( m = 0; m < SUBUSED; m++){

        tmp = 1.0f;

        for(n = 0; n < DEREVERB_DELAYLINE; n++) {
            tmp*=decayInGainOut[m];
        }

        tmp*=delaylineBufEnd[m];

        RSR = fasterdiv(tmp, delaylineBufN[m]); // Reverb-to-Signal Ratio

        decayInGainOut[m] = min(maxgain, max(mingain, (1.0f-RSR))); // Wiener gain
    }
}


#ifdef UNITTEST
/***************************************************************************
* UNITTEST_DEREVERB
*    Desc.: run-trhough of dereverb with 10 blocks equal random input on
*           louds and micinput.
*           The unittest checks mic.gain and noisevector output.
*           The noise output is not equal to matlab, but the nfgn has been
*           verified manually.
***************************************************************************/
#include "unittest.h"

/***************************************************************************
* UNITTEST_DEREVERB
**************************************************************************/
void unittest_dereverb()
{
}

#endif /* UNITTEST */
