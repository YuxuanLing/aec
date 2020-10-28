#ifndef SYNTH_H
#define SYNTH_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Geir Ole �verby, Tandberg Telecom AS
 * Co-author     : Johan Malmstigen, S&T ES
 *
 * Switches      : -
 *
 * Description   : Interface for the synthese filter in the Echo Canceller
 *         module.
 *
 *         To create an instance of the synthese filter, call
 *          synth_create(), followed by
 *          synth_init()
 *
 *         To transform n frames of FFTSIZE subbands from frequency-
 *         domain to time-domain, call
 *          synth_process()
 *
 **************************************************************************/

#include "audec_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*-------------------------------------------------------------------------
 * -- EXPORTED TYPES AREA
 *-------------------------------------------------------------------------*/
typedef struct SYNTH * SYNTH_PTR;

/*-------------------------------------------------------------------------
 * -- EXPORTED FUNCTIONS AREA
 *-------------------------------------------------------------------------*/


/******************************************************************
 * SYNTH_CREATE
 *  Description: Allocate memory for the synthese filter.
 *
 *  Parameters:  -
 *
 *  Returns:     Pointer to SYNTH structure used by the synthese filter
 *
 *  Note:        Memory for the SYNTH structure is allocated here. To free
 *               this memory, synth_destroy has to be called.
 *****************************************************************/
SYNTH_PTR synth_create(void);


/******************************************************************
 * SYNTH_DESTROY
 *  Description:  Frees memory for fftst matrix and synth struct
 *
 *  Parameters:   pSynth - Pointer to the SYNTH structure that shall be freed
 *
 *  Restrictions: synth_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void synth_destroy(SYNTH_PTR pSynth);


/******************************************************************
 * SYNTH_INIT
 *  Description:    Initialization of the struct that contains the necessary
 *                  data for each signal to be analysed
 *
 *  Parameters:     pSynth - Pointer to the SYNTH structure to be initialized
 *                  sampleFreqKhz - integer, 16 or 48 kHz
 *
 *  Returns:        None
 *
 *  Restrictions:   synth_create has to be called once before this function
 *                  is called.
 *****************************************************************/
void synth_init(SYNTH_PTR pSynth);


/***************************************************************************
 * SYNTH SET FILTERBANK PROTOTYPE FILTER
 *  Description: select synthesis filter
 *
 *   Parameters: pSynth - Pointer to the SYNTH structure
 *               sampleFreqKhz - integer, 16 or 48
 *               filtertype - integer
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: synth_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void synth_setFbPrototypefilt(SYNTH_PTR pSynth, int filtertype);


void synth_status(SYNTH_PTR pSynth);

/******************************************************************************
 * SYNTH_PROCESS
 *
 *  Description: Filterbank: Takes a signal from the frequency-domain and
 *       puts the corresponding signal in the time-domain, after
 *       some filtering.
 *
 *       The follwoing C code is used to synthezie the in signal
 *       containing a set of complex frequecy bands to the out
 *       signal. Only half of the subbands are represnted in the
 *       insignal, the rest are calculated based on the fact that
 *       the signal the in signal was analyzed from should be real.
 *       The filtering strongly depends on the fact that the length
 *       of the synthezie filter h, (not including the last part
 *       filled with only zeros) is a multiplum of the number of
 *       subbands, FFTSTLY. This makes it possible to reduce the
 *       number of filtertaps used in computing each sample out to
 *       (SFIRLEN/FRAMESIZE) (rounded up) by computing an ifft of
 *       the subbands first. Note also that this routine will
 *       overwrite the last(FFTSTLY) real subbands with 0, because
 *       we do not want these to contribute. h is padded with zeroes
 *       until it is a multiplum of FRAMESIZE to acheive a simpler
 *       synth_build.
 *
 *   Parameters: pSynth        - Pointer to the SYNTH structure
 *               inbuf         - Pointer to buffer in frequency-domain (subbands)
 *               outbuf        - Pointer to buffer in time-domain (sound samples)
 *               bDoPostfilter - Applies postfiltering if true.
 *
 *   Returns: None
 *
 *   Note: The buffersize is set in the synth_create operation. The
 *         output consists of FRAMESIZE real values.
 *
 *   Note (TMS): It is required (by the fft) that in-vector lies on a
 *               double-boundary
 *
 *   Globals used: -
 *
 *   Restrictions: This function is dependent on the fft-module, which has
 *                 to be initialized before this function is called.
 ***************************************************************************/
void synth_process(SYNTH_PTR pSynth, COMPLEX32 * inbuf, float * outbuf);


/* SYNTH_LOADFFTST ******************************************************
 * Auth.: JPS
 * Desc.: loads pSynth->fftst into iram
 ************************************************************************/
void synth_loadfftst(SYNTH_PTR pSynth);


/* SYNTH_LOADSFIR ********************************************************
 * Auth.: JPS
 * Desc.: loads sfir into iram
 ************************************************************************/
void synth_loadSfir(SYNTH_PTR pSynth);


/******************************************************************************
 * SYNTH_PROCESSLOOPDETECT
 *****************************************************************************/
void synth_processLoopDetect(SYNTH_PTR pSynth,
                             float *fftst,
                             float *in,
                             float * out);

/* SYNTH_LOADSYNTH *******************************************************
 * Auth.: JPS
 * Desc.: loads synthese struct into iram
 * Returns: pointer to synthese struct in iram
 ************************************************************************/
SYNTH_PTR synth_loadSynth(SYNTH_PTR pSynthIn);

/* SYNTH_FLUSHSYNTH ******************************************************
 * Auth.: JPS
 * Desc.: flush synthese struct to sdram
 ************************************************************************/
void synth_flushSynth(SYNTH_PTR pSynthSdram, SYNTH_PTR pSynthIram);


#ifdef UNITTEST
/******************************************************************************
 * UNITTEST
 *   Description:
 *       Run-trhough of synthese filter with 10 blocks of input signal.
 *       The test input signal to this unittest is the fftout (mic.micfft)
 *       of analysefilter in the matlab model when the input to analyse is random noise
 *       The result is verified with matlab.
 *       Both 16 kHz and 48 kHz implementation is verified.
 *****************************************************************************/
void unittest_synthese(void);
#endif /* UNITTEST */

#ifdef __cplusplus
}
#endif

#endif /* SYNTH_H */
