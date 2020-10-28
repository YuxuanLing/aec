#ifndef NOISEREDUCTION_H
#define NOISEREDUCTION_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 * Co-author     : Olof Lindström (OLI), S&T ES
 *
 * Description   : Part of echocanceller. Main objectives:
 *                 1) Reduce noise by decreasing the gain in subbands not
 *                    containing talk.
 *                 2) Supress echoes by attenuating signal when far-end-talk.
 *                 3) Mute signal when echocanceller is not working (exception
 *                    occurs)
 *                 4) Fill in synthetic noise when 2) or 3) mutes signal.
 *                 5) Produce some input to mixer and hifigain.
 *
 **************************************************************************/

#include "audec_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*-------------------------------------------------------------------------
 * -- INCLUDE FILES AREA
 *-------------------------------------------------------------------------*/

#include "keyclickremoval.h"
#include "audtypes.h"
#include "dereverb.h"
#include <stdbool.h>

/*-------------------------------------------------------------------------
 * -- PREPROCESSOR CONSTANTS AREA
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * -- EXPORTED TYPES AREA
 *-------------------------------------------------------------------------*/
typedef struct NOIRED * NOIRED_PTR;

/*-------------------------------------------------------------------------
 * -- EXPORTED DATA AREA
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * -- EXPORTED MACROS AREA
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * -- EXPORTED FUNCTION PROTOTYPES AREA
 *-------------------------------------------------------------------------*/

/***************************************************************************
 NOISEREDUCTION_CREATE
 *   Description: Allocates memory for struct noired.
 *
 *   Parameters:   -
 *
 *   Returns:      pointer to NOIRED struct.
 *
 *   Note:         -
 *
 *   Globals used: -
 *
 *   Restrictions: -
 *
 ***************************************************************************/
NOIRED_PTR noisereduction_create(void);

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
 *
 ***************************************************************************/
void noisereduction_destroy(NOIRED_PTR pNoired);

/***************************************************************************
 * NOISEREDUCTION_INIT
 *  Description: Initialises members in struct noired.
 *
 *  Parameters:   pNoired
 *
 *  Restrictions: Memory pointed by pNoired shall be allocated by
 *                calling function noisereduction_create.
 *
 ***************************************************************************/
void noisereduction_init(NOIRED_PTR pNoired, FILTERBANK_USE_TYPE filterbankUse);


/***************************************************************************
 * NOISEREDUCTION_PROCESS
 *  Description: Main noisereduction routine
 *               Determines if the input signals are noisy or not.
 *
 *   Parameters: pAnalyse       - pointer to ANALYSE struct of mic
 *               excgain        - exceptionhandler gain from MICEXCHNDL struct.
 *               pNoired        - pointer to NOIRED struct containing variables
 *                                defining noise reduction filters.
 *               estecho        - Echo estimation from echo compensation.
 *               gain           - gain output from noise reduction.
 *               noise          - noise output from noise reduction.
 *               fftSize        - length of estecho array
 *               nlpconst       - non linear processing const.
 *               noiseredOn     - Noisereduction turned on/off
 *             GSM_noiseDetected-
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
 *
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
                             float *decay);

/***************************************************************************
 NOISEREDUCTION_LOADGAIN
 *   Description: copies the gain vector to iram  using dmax.
 *
 *   Parameters:   - gainSdram, pointer to gain vector
 *                 - subused, number of subbands
 *
 *   Returns:      - pointer to gain vector in iram.
 *
 *   Note:         - should use memxfer_waitAddress on the returned pointer
 *                   before use in order to confirm completed dmax transfer
 ***************************************************************************/
float* noisereduction_loadGain(float *gainSdram);


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
NOIRED_PTR noisereduction_loadNoired(NOIRED_PTR pNoisereduction);

/***************************************************************************
 NOISEREDUCTION_FLUSHNOIRED
 *   Description: flushes the noisereduction struct to sdram  using dmax.
 *
 *   Parameters:   - pNoiredSdram
 *                 - pNoiredIram
  ***************************************************************************/
void noisereduction_flushNoired(NOIRED_PTR pNoiredSdram, NOIRED_PTR pNoiredIram);

void noisereduction_getComfortNoise(NOIRED_PTR pNoired,
                                   COMPLEX32 *micfft,      /* subband micsignal (i) */
                                   COMPLEX32 *noise,       /* subband noise (o) */
                                   float nlpgain);     /* zero delay nlpgain (i) */


/****************************************************************************
 * Functions for setting testvariables inside noisereduction
 ***************************************************************************/
void noisereduction_setDereverb(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlp(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setComNoise(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setComNoiseAmp(NOIRED_PTR pNoired, float value);
void noisereduction_setNoiseRed(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlpSub(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlpHi(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlpExtragain(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlpfullg(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlpfullgHi(NOIRED_PTR pNoired, bool onoff);
void noisereduction_setNlpSubusedEndHi(NOIRED_PTR pNoired, int value);
void noisereduction_setNlpTransition(NOIRED_PTR pNoired, int value);
void noisereduction_status(NOIRED_PTR pNoired);
void noisereduction_setDebug(NOIRED_PTR pNoired, int value);
void noisereduction_setNcntBoot(NOIRED_PTR pNoired, int iCntBoot);
void noisereduction_ZerodelayStatus(NOIRED_PTR pNoired);


#ifdef UNITTEST
/***************************************************************************
 * UNITTEST_NOISEREDUCTION
 *    Desc.: run-trhough of noisereduction with 10 blocks equal random input on
 *           louds and micinput.
 *           The unittest checks mic.gain and noisevector output.
 *           The noise output is not equal to matlab, but the nfgn has been
 *           verified manually.
 ***************************************************************************/
void unittest_noisereduction(void);
#endif /* UNITTEST */

#ifdef __cplusplus
}
#endif

#endif /* NOISEREDUCTION_H */
/*----------------------------- End of File -------------------------------*/
