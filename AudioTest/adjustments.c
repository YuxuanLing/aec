#include "adjustments.h"
#include "audec_defs.h"

/* CALYPSO-filter framesize=FRAMESIZE48
 *   filter constructed in MatLab with two consecutive bilinear IIR-filters
 *   (to obtain high "q-value" in the resulting bi-quadratic filter):
 *   f6db=10^2.07453;
 *   [b,a] = butter(1,f6db/24000,'high');
 *   a=conv(a,a); b=conv(b,b);
 *
 *   This filter has -6dB at 119 Hz (original calypso-filter)
 */
#define NOIREDHP_b0 0.98463833886716f
#define NOIREDHP_b1 -1.96927667773431f
#define NOIREDHP_b2 0.98463833886716f
#define NOIREDHP_a0 1.00000000000000f /* is thus removed from routine ..*/
#define NOIREDHP_a1 -1.96915777235858f
#define NOIREDHP_a2 0.96939558311005f

/*****************************************************************************
 * PREFILTER
 *   Auth.: IFA/GEO
 *   Desc.: Lowering low frequencies while rasing the level of high freqs.
 *   Note:  Works only when the filterlength==2!!
 ****************************************************************************/
void ecprefilt(float * buf, float * prefiltStatePtr, float prefiltercoeff, float prefiltergain)
{
    int i;
    float tmp;
    float prefiltState = *prefiltStatePtr;

    for( i=0; i<FRAMESIZE; i++ )
    {
        tmp=buf[i];
        buf[i] = (tmp + prefiltercoeff * prefiltState) * prefiltergain;
        /* (a[0]*y[n] = b[0]*x[n]+b[1]*x[n-1]) , inplace */
        prefiltState=tmp;
    }

    *prefiltStatePtr = prefiltState;
}

/************************************************************************
 * POSTFILTER
 *   Auth.: IFA/GEO
 *   Desc.: Raises low frequencies while lowering the level of high
 *          frequencies.(Inverse of PREFILTER)
 *   Note:  Works only when the filterlength==2!!
 ************************************************************************/
void ecpostfilt(float * buf, float * postfiltStatePtr, float postfiltercoeff, float postfiltergain)
{
    int i;
    float postfiltState = *postfiltStatePtr;

    for( i=0; i<FRAMESIZE; i++ )
    {
        postfiltState = buf[i]*postfiltergain - postfiltercoeff * postfiltState;
        /* (a[0]*y[n] = b[0]*x[n]-a[1]*y[n-1]) , inplace */
        buf[i]=postfiltState;
    }
    *postfiltStatePtr = postfiltState;
    /* this is different compared to Matlab, but gives the same filtering. */
}

/* DC-FILTER ********************************************************** */
/* Auth.: GEO                                                           */
/* Desc.: Removes (lowers) constant and almost constant frequencies.    */
/* Note:  This is the ecdcfilt-routine from the "old" MatLab-model.     */
/* ******************************************************************** */

float dcremove(float * in, float dc_offset, float tap)
{
  int i;
  float s=0;
  float a;

  for(i=0;i<FRAMESIZE;i++)
  {
    a=in[i];
    in[i] = a - dc_offset;
    s += a;
  }
  s -= dc_offset*FRAMESIZE;
  dc_offset+=s*tap;
  return dc_offset;
}

 /***************************************************************************
 * NOIREDHPFILT
 *   Auth.: IFA
 *   Desc.: Highpass IIR-filter attenuating frequencies lower than 185Hz - thus not
 *          deteriorating speech-intelligibility. The motivation for this is to
 *          avoid too loud sounds in the loweer frequncies. Should be by-passed
 *          when noisereduction is turned off.
 ****************************************************************************/
void noiredHPfilt(float * delayLine, float * const  inbuf, float * outbuf)
{
  /* Temporarily disabled this function for 674X evaluation project,
     need properly built library like mentioned in comment below */
#if defined(_TMS320C6X) && !defined(CHIP_674X)

  /* ------ */
  /* using DSPLIB routine - see spru657 and spra947 */
  /*void DSPF_sp_biquad(float *x,      float pointer to in-buffer  */
  /*                    float *b,      fir-filter part coefficients  */
  /*                    float *a,      iir-filter part coefficients  */
  /*                    float *delay,  delayline  */
  /*                    float *r,      float pointer to out-buffer  */
  /*                    int nx )       inbuffer framesize  */
  /* note: in order to run correctly, the interrupt threshold must be larger than 4*480+76 = 1996  */
  /* ------ */
  DSPF_sp_biquad(inbuf, noiredHP_B, noiredHP_A, delayLine, outbuf, FRAMESIZE);

#else

  int i;

  /*  float a0Inv = 1/NOIREDHP_a0;  (removed since a0=1)*/

  /* first two samples are exceptions from for-loop */
  outbuf[0] =   ( NOIREDHP_b0*inbuf[0] ) + delayLine[0];

  outbuf[1] =   ( NOIREDHP_b0*inbuf[1]  + NOIREDHP_b1*inbuf[0] )
    - ( NOIREDHP_a1*outbuf[0] ) + delayLine[1] ;

  /* next two samples are removed from for-loop (optimization) */
  outbuf[2] =   ( NOIREDHP_b0*inbuf[2]  + NOIREDHP_b1*inbuf[2-1] + NOIREDHP_b2*inbuf[2-2] )
    - ( NOIREDHP_a1*outbuf[2-1] + NOIREDHP_a2*outbuf[2-2] ) ;
  outbuf[3] =   ( NOIREDHP_b0*inbuf[3]  + NOIREDHP_b1*inbuf[3-1] + NOIREDHP_b2*inbuf[3-2] )
    - ( NOIREDHP_a1*outbuf[3-1] + NOIREDHP_a2*outbuf[3-2] ) ;

  for(i=4; i<FRAMESIZE; i++)
    {
      outbuf[i] =   ( NOIREDHP_b0*inbuf[i]  + NOIREDHP_b1*inbuf[i-1] + NOIREDHP_b2*inbuf[i-2] )
  - ( NOIREDHP_a1*outbuf[i-1] + NOIREDHP_a2*outbuf[i-2] ) ;
    }

  /*update delayline*/
  delayLine[0] =  ( NOIREDHP_b1*inbuf[FRAMESIZE-1] + NOIREDHP_b2*inbuf[FRAMESIZE-2] )
    - ( NOIREDHP_a1*outbuf[FRAMESIZE-1] + NOIREDHP_a2*outbuf[FRAMESIZE-2] ) ;
  delayLine[1] =  ( NOIREDHP_b2*inbuf[FRAMESIZE-1] )
    - ( NOIREDHP_a2*outbuf[FRAMESIZE-1] ) ;
#endif

}
