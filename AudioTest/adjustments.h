#ifndef ADJUSTMENTS_H
#define ADJUSTMENTS_H

/* pre-, post-filters, powercalculation 16kHz, highpass-iir-filter. */

/* PREFILTER ********************************************************** */
/* Auth.: IFA                                                           */
/* Desc.: Lowering low frequencies while rasing the level of high freqs.*/
/* Note:  Works only when the filterlength==2!!                         */
/*        framesize mod FRAMESIZE16 should be 0                         */
/*        For Texas, framesize has to be FRAMESIZE16 or 2*FRAMESIZE16   */
/* ******************************************************************** */
void ecprefilt(float * buf, float * prefiltState, float prefiltercoeff, float prefiltergain);

/* POSTFILTER ********************************************************* */
/* Auth.: IFA                                                           */
/* Desc.: Raises low frequencies while lowering the level of high       */
/* frequencies.(Inverse of PREFILTER)                                   */
/* Note:  Works only when the filterlength==2!!                         */
/*        framesize mod FRAMESIZE16 should be 0                         */
/*        For Texas, framesize has to be FRAMESIZE16 or 2*FRAMESIZE16   */
/* ******************************************************************** */
void ecpostfilt(float * buf, float * postfiltState, float postfiltercoeff, float postfiltergain);

/* DC-FILTER ********************************************************** */
/* Auth.: GEO                                                           */
/* Desc.: Removes (lowers) constant and almost constant frequencies.    */
/* Note:  This is the ecdcfilt-routine from the "old" MatLab-model.     */
/* ******************************************************************** */
float dcremove(float * in, float dc_offset, float tap);

/* POW16 ************************************************************* */
/* Auth.: IFA                                                          */
/* calculation of the power of the lower frequencies, for use in       */
/* the calculation of the hifigain. It is separated from the rest of   */
/* hifigain_noise for optimisation (to avoid superfluous cache misses).*/
/* ******************************************************************* */
float pow16_calc(float * inbuf16);

/* LEV16 ************************************************************* */
/* Auth.: IFA                                                          */
/* calculation of the level of the lower frequencies, for use in       */
/* the calculation of the hifigain. It is separated from the rest of   */
/* hifigain_noise for optimisation (to avoid superfluous cache misses).*/
/* ******************************************************************* */
#ifdef _TMS320C6X
float lev16_calc(float * inbuf16);
#endif

/* noiredHPfilt *********************************************************** */
/* Auth.: IFA                                                           */
/* Highpass IIR-filter attenuating frequencies lower than 185Hz - thus not  */
/* deteriorating speech-intelligibility. The motivation for this is to      */
/* avoid too loud sounds in the loweer frequncies. Should be by-passed      */
/* when noisereduction is turned off.                                       */
/* ************************************************************************ */
void noiredHPfilt(float * delayLine, float * const  inbuf, float * outbuf);


/* noiredHPfilt *********************************************************** */
/* Auth.: IFA, modified EHO                                                 */
/* Highpass IIR-filter attenuating frequencies higher than 6800Hz-thus not  */
/* deteriorating speech-intelligibility. The motivation for this is to      */
/* avoid echo in full band EC, in bands where echo is not cancelled.        */
/*                                                                          */
/* ************************************************************************ */
void noiredLPfilt(float * delayLine, float * const  inbuf, float * outbuf);

#endif
