#ifndef ANALYSE_H
#define ANALYSE_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes, Tandberg Telecom AS
 * Co-author     : Johan Malmstigen, S&T ES
 *
 * Description   : Interface for the analyse filter in the Echo Canceller
 *         module.
 *
 *         To create an instance of the analyse filter, call
 *          analyse_create(), followed by
 *          analyseInit()
 *
 *         To transform a n frames of n*160 or n*480 sound samples from
 *         time-domain to frequency domain, call analyse_process()
 *
 ***************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/*-------------------------------------------------------------------------
 * -- EXPORTED TYPES AREA
 *-------------------------------------------------------------------------*/
typedef struct ANALYSE * ANALYSE_PTR;

/*-------------------------------------------------------------------------
 * -- EXPORTED FUNCTIONS AREA
 *-------------------------------------------------------------------------*/

#include "audtypes.h"

/***************************************************************************
 * ANALYSE CREATE
 *  Description: Claims memory for the analyze filter, returns pointer.
 *
 *   Parameters: -
 *
 *   Returns:   A pointer to a ANALYSE structure used by the analyse filter
 *
 *   Note:      Memory for the ANALYSE structure is allocated here. To free
 *              this memory, analyse_destroy has to be called.
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
ANALYSE_PTR analyse_create(FILTERBANK_USE_TYPE filterbankUse);

/***************************************************************************
 * ANALYSE DESTROY
 *  Description: Frees memory. Should not be necessary.
 *
 *   Parameters: pAnalyse - Pointer to the ANALYSE structure that shall be
 *              freed
 *
 *      Returns: -
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: analyse_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void analyse_destroy(ANALYSE_PTR pAnalyse, FILTERBANK_USE_TYPE filterbankUse);

/***************************************************************************
 * ANALYSE INIT
 *  Description: Initialization of the struct that contains the necessary
 *               data for each signal to be analysed
 *
 *   Parameters: pAnalyse - Pointer to the ANALYSE structure to be initialized
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: analyse_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void analyse_init(ANALYSE_PTR pAnalyse, FILTERBANK_USE_TYPE filterbankUse);

/***************************************************************************
 * ANALYZE SET FILTERBANK PROTOTYPE FILTER
 *  Description: select analysis filter
 *
 *   Parameters: pAnalyse - Pointer to the ANALYSE structure
  *               filtertype - integer
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: analyse_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void analyse_setFbPrototypefilt(ANALYSE_PTR pAnalyse, int filtertype);

void analyse_status(ANALYSE_PTR pAnalyse);

/***************************************************************************
 * ANALYSE PROCESS
 *  Description: Filterbank: Takes a sampled signal from the time domain
 *               and puts out the corresponding signal in the frequency-
 *               domain, after some filtering (PRE- and DC-filters).
 *
 *   Parameters: pAnalyse  - Pointer to the ANALYSE structure
 *               inbuf     - Pointer to n*480 new samples
 *               afir      - Pointer to analyse filter
 *
 *   Returns:    0
 *
 *   Note:       The buffersize is set in the analyse_create operation. The
 *               filter output is stored in member fftout in ANALYSE struct.
 *               The output consists of n*384 complex values,
 *               dependent of the sampling frequency. In the Texas
 *               implementation, n=2 and the second part is the conjugate
 *               of the first part.
 *
 *   Remark:     FFTSIZE is 256 or 768. 1 complex value is made up by two parts,
 *               one real part and one imaginary part. The filter output
 *               does not include the conjugated part.
 *
 *   Globals used: -
 *
 *   Restrictions: This function is dependent on the fft-module, which has
 *                 to be initialized before this function is called.
 ***************************************************************************/
void analyse_process(ANALYSE_PTR pAnalyse, float * inbuf, float * afir);


/* ANALYSE_GAINADJUST
 * Auth.: JPS
 * Desc.: Multiply output from analyse filter with gain from loudspeaker volume
 *        to keep the loudspeaker samples in range for the delta-calculation in
 *        echocomp.
 *        Only use on lsfft.
 */
void analyse_gainAdjust(ANALYSE_PTR pAnalyse, float loudsGain);
void analyse_updateInbufStart(ANALYSE_PTR pAnalyse, int index);
int analyse_getInbufStart(ANALYSE_PTR pAnalyse);

/* ANALYSE_LOADINBUF *****************************************************
 * Auth.: JPS
 * Desc.: loads pAnalyse->inbuf into iram
 ************************************************************************/
void analyse_loadInbuf(ANALYSE_PTR pAnalyse);

/* ANALYSE_LOADAFIR ******************************************************
 * Auth.: JPS
 * Desc.: loads afir into iram
 ************************************************************************/
float * analyse_loadAfir(ANALYSE_PTR pAnalyse);

/* analyse_loadAnalyse ***************************************************
 * Auth.: JPS
 * Desc.: loads pAnalyse struct into iram
 * Returns: ptr to analyse struct in iram
 ************************************************************************/
ANALYSE_PTR analyse_loadAnalyse(ANALYSE_PTR pAnalyseIn);

/* analyse_flushAnalyse ***************************************************
 * Auth.: JPS
 * Desc.: flush pAnalyse struct to sdram
 ************************************************************************/
void analyse_flushAnalyse(ANALYSE_PTR pAnalyseSdram, ANALYSE_PTR pAnalyseIram);

#ifdef UNITTEST
/************************************************************************
 * UNITTEST
 * Desc.: run-trhough of analyse filter for both 16kHz and 48kHz
 *        with 0.3sec random input.
 *        Check against testvector verified in matlab
 *        Pre-filter is active
 *        All sub-bands to FFTSIZE/2
 ************************************************************************/
void unittest_analyse(void);
#endif /* UNITTEST */


#ifdef __cplusplus
}
#endif

#endif /* ANALYSE_H */
