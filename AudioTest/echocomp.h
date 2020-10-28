#ifndef ECHOCOMP_H
#define ECHOCOMP_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 *         Geir Ole Øverby (GEO), Tandberg Telecom AS.
 * Co-author     : Torgeir Grothe Lien (TGL), Tandberg Telecom AS.
 *
 **************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#include "delayEstimation.h"
#include <stdbool.h>

typedef struct ECHOCOMP * ECHOCOMP_PTR;
typedef struct DECAY_ATTACK DECAY_ATTACK;

/***************************************************************************
 * ECHOCOMP_CREATE
 *  Description: Claims memory for the echo compensation filters, returns pointer.
 *
 *   Parameters: -
 *
 *   Returns:   A pointer to an ECHOCOMP structure
 *
 *   Note:      Memory for the ECHOCOMP structure is allocated here. To free
 *              this memory, echocomp_destroy has to be called.
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
ECHOCOMP_PTR echocomp_create(int channels);


/***************************************************************************
 * ECHOCOMP_DESTROY
 *  Description: Frees memory. Should not be necessary.
 *
 *   Parameters: pEchocomp - Pointer to the ECHOCOMP structure that shall be
 *              freed
 *
 *      Returns: -
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: echocomp_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void echocomp_destroy(ECHOCOMP_PTR pEchocomp);

/***************************************************************************
 * ECHOCOMP_INIT
 *  Description: Initialization of the struct
 *
 *   Parameters: pEchocomp - Pointer to the ANALYSE structure to be initialized
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: echocomp_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void echocomp_init(ECHOCOMP_PTR pEchocomp);


/***************************************************************************
 * ECHOCOMP PROCESS
 *  Description: Main routine for echocpmpansation processing.
 *
 *   Parameters: pEchocomp  - Pointer to the ECHOCOMP structure
 *               lsBuf      - Pointer to the loudspeaker fftout, must be double-word aligned
 *               micin      - Pointer to fftsize new samples
 *               echoest    - Pointer to estimated echo, subused*2 samples
 *               systemGain -
 *
 *   Returns:    0
 *
 *   Note:
 *
 *   Globals used: -
 *
 *   Restrictions: - lsBuf must be double-word aligned
 *                 - This function is dependent on apa.c
 ***************************************************************************/
int echocomp_process(ECHOCOMP_PTR pEchocomp,
                     const COMPLEX32 *lsBuffer,
                     COMPLEX32 *echoest,
                     const COMPLEX32 *micin,
                     const float excgain,
                     bool runShift);

int echocomp_processStereo(ECHOCOMP_PTR pEcho,
                           const COMPLEX32 *lsBuf,
                           const COMPLEX32 *lsBuf2,
                           COMPLEX32 *echoest,
                           const COMPLEX32 *micin,
                           const float excgain,
                           bool runShift);

/***************************************************************************
 * ECHOCOMP CALC NEW DELTA
 *   Auth.: EHO(GEO)
 *   Desc.: This routine calculates the delta parameter which is used to
 *          control the speed vs stability tradeoff in the filter adaption.
 *          Large delta improves stability but will limit convergence speed if
 *          chosen too large. This is often refered to as regularization.
 ***************************************************************************/
void echocomp_calcNewDelta(ECHOCOMP_PTR pEchocomp,
                           const float excgain,
                           const int subband_start,
                           const int subband_end);


/* ECHOCOMP_loadAdapttaps************************************************
 * Auth.: JPS
 * Desc.: loads adapt filtertaps into iram
 ************************************************************************/
void echocomp_loadAdaptTaps(ECHOCOMP_PTR pEchocomp);

/* ECHOCOMP_loadFixedTaps************************************************
 * Auth.: JPS
 * Desc.: loads fixed filtertaps into iram
 ************************************************************************/
void echocomp_loadFixedTaps(ECHOCOMP_PTR pEchocomp);

/* ECHOCOMP_getAerl
 * Calculating and returning the estimated aerlInverse.
 * It is basically the maxTap in the fixed filter, with a time smoothing.
 * The aerlInvers estimate is dependent on the miclevel setting and the loudspeaker volume.
 * Returns: the estimated aerlInvers
 */
float echocomp_getAerl(ECHOCOMP_PTR pEchocomp, float miclevel, float loudsGain);


/* ECHOCOMP_getDecay
 * Returns: pointer to decay
 */
//float* echocomp_loadDecay(ECHOCOMP_PTR pEchocomp);

void echocomp_setCalcNewDelta(ECHOCOMP_PTR pEchocomp, bool onoff);
void echocomp_setMinDelta(ECHOCOMP_PTR pEchocomp, float value);
void echocomp_setMaxDelta(ECHOCOMP_PTR pEchocomp, float value);
void echocomp_setFilter(ECHOCOMP_PTR pEchocomp, int value);
void echocomp_setZero(ECHOCOMP_PTR pEchocomp, int value);
void echocomp_setDebug(ECHOCOMP_PTR pEchocomp, int value);
void echocomp_setStereoDelta(ECHOCOMP_PTR pEchocomp, float value);
void echocomp_status(ECHOCOMP_PTR pEchocomp);
void echocomp_moveFilters(ECHOCOMP_PTR pEchocomp, const struct DELAY_ESTIMATION * delayEstimation, int delay_diff);

#ifdef UNITTEST
/* unittest_echocomp **************************************************** */
/* Auth.: Jens Petter Stang (JPS)                                         */
/* Desc.: simple run-trhough of echocomp                                  */
/* ********************************************************************** */
void unittest_echocomp(void);
#endif /* UNITTEST */

#ifdef __cplusplus
}
#endif

#endif
