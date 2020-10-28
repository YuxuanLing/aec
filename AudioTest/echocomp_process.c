/***************************************************************************
*                              A U D I O
*--------------------------------------------------------------------------
*                  (C) Copyright Tandberg Telecom AS 2004
*==========================================================================
*
* Author        : Espen Holmbakken (EHO)
*                 Geir Ole �verby  (GEO)
*
* Co-author     : -
*
* Switches      : <COMPILER SWITCH HERE>
*           <add description here>
*
* Description   : <add description here>
*
* Note      : <notes here>
*
* Documentation : <Xxxxxx-Document.doc>
*
* ------------------------------------------------------------------------
* Major changes made (complement to cvs log)
* ------------------------------------------------------------------------
* yyyy-mm-dd    : <signature>
*         Modification made.
*
* yyyy-mm-dd    : <signature>
*         File created.
*
**************************************************************************/
#include "apa.h"
#include "echocomp_priv.h"
#include <math.h>
#ifdef __ARM_NEON__
#include <arm_neon.h>
//#pragma thumb
#define fabsf fabs      // compiler workaround
#else
#include <string.h>
#endif
#define fastdiv(a,b) ((a)/(b))
#define fasterdiv(a,b) ((a)/(b))
#define _rcpsp(x) (1.0f/(x))
#define _rsqrsp(x) (1.0f/sqrtf(x))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

/* ECHOCOMP SHIFT *******************************************************
* Auth.: IFA, later modified by OBI
* Desc.: Delays loudspeaker subbands.
* Each subband has a different delay-length, given by filtlen vector.
************************************************************************/
void echocomp_shift(ECHOCOMP_PTR pEchocomp,
                    int ch, // channel
                    const COMPLEX32 *  lsBuf,
                    const int subband_start,
                    const int subband_end)
{
  ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];
  int i;
#if defined(ENV_IOS)
  lsBuf += SUBUSED_START+subband_start;
#endif

  for ( i = subband_start; i < subband_end; i++ )
  {  
    ECHOCOMP_SUBBAND *subband = &channel->subband[i];
    COMPLEX32 * dl = subband->delayline;
    int len = subband->filtlen;

#if defined(ENV_IOS)
    if (dl > subband->delayline_base)
    {
      dl[len].re = 0.0;
      dl[len].im = 0.0;
      subband->delayline = --dl;
    }
    else
    {     
#if defined(__ARM_NEON__)
      /*Since subband->delayline is allocated with ((len+3)&3)+4, there is enough buffer to do the following safely*/
      len = (len+3)&~3; 
      dl = subband->delayline_base+len;
      if (len&4) 
      {
         dl  -=4;
         len -=4;
         vst2q_f32((float*)(dl+9), vld2q_f32((float*)dl));
      }
      
      for ( ; len>0; len-=8)
      {
         dl -= 8;
         vst4q_f32((float*)(dl+9), vld4q_f32((float*)dl));
      }
#else
      dl = subband->delayline_base;
      memmove(dl+9, dl, len*sizeof(COMPLEX32)); 
#endif
      subband->delayline = dl = dl+8;
    }
    dl->re = lsBuf->re;
    dl->im = lsBuf->im;
    lsBuf++;
#else
    memmove(&dl[1], dl, len*sizeof(COMPLEX32));
    dl->re = lsBuf[SUBUSED_START+i].re;
    dl->im = lsBuf[SUBUSED_START+i].im;
#endif
  }
}

/* ECHOCOMP SUBTRACT ************************************************** */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine simply computes (mic[] - est[]) and places the   */
/* the answer in canc[].                                                */
/* ******************************************************************** */
void echocomp_subtract(COMPLEX32 * canc,
                       const COMPLEX32 *  mic,
                       const COMPLEX32 *  est,
                       const int subband_start,
                       const int subband_end)
{
  int i = subband_start;

#if defined(__ARM_NEON__) && defined(ENV_IOS)
  {
    int end = (subband_end-subband_start) &~3;
    end += subband_start;
    for( ; i < end; i+=4 )
    {
      float32x4x2_t val_mic = vld2q_f32((float32_t*)&mic[i]);
      float32x4x2_t val_est = vld2q_f32((float32_t*)&est[i]);
      float32x4x2_t val_can;
      float32x4_t re = vsubq_f32(val_mic.val[0], val_est.val[0]);
      float32x4_t im = vsubq_f32(val_mic.val[1], val_est.val[1]);
      val_can.val[0] = re;
      val_can.val[1] = im;
      vst2q_f32((float32_t*)&canc[i], val_can);
    }
  }  
#endif

  for ( ; i < subband_end; i++)
  {
    canc[i].re = mic[i].re - est[i].re;
    canc[i].im = mic[i].im - est[i].im;
  }
}

/* ECHOCOMP EVALUATION ************************************************ */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine calculates the amount of signal in buf. This is  */
/* done by summing the absolute values of all the real and imaginary    */
/* parts, and then filter the value obtained.                           */
/* ******************************************************************** */
float echocomp_eval(float *level,
                    const COMPLEX32 *  buf,
                    const int subband_start,
                    const int subband_end)
{
  int i = subband_start;
  float sum_abs_real = 0.0f, sum_abs_imag = 0.0f;
  float sum=0.0;

#ifdef PRESERVE_ORIGINAL_BUGS
  buf = (COMPLEX32 *)((float *)buf + subband_start); // this looks wrong, but aligns with original code
#else
  buf += subband_start;
#endif

#if defined(__ARM_NEON__) && defined(ENV_IOS)
  {
    float32x4_t sum_re;
    float32x4_t sum_im;
    float32x2_t sum_2f;
    float32x4_t abs0, abs1;
    float32x4x2_t tmpd;
    int end;

    sum_re = vdupq_n_f32(0.0f);
    sum_im = vdupq_n_f32(0.0f);
    end = subband_start + ((subband_end-subband_start)&~3);
    for( ; i < end; i+=4 )
    {
      tmpd = vld2q_f32((float32_t*)buf);
      buf += 4;
      abs0 = vabsq_f32(tmpd.val[0]);
      abs1 = vabsq_f32(tmpd.val[1]);
      sum_re = vaddq_f32(sum_re, abs0);
      sum_im = vaddq_f32(sum_im, abs1);
    }
     
    sum_re = vaddq_f32(sum_re, sum_im);
    sum_2f = vpadd_f32(vget_high_f32(sum_re), vget_low_f32(sum_re));
    sum = vget_lane_f32(sum_2f, 0) + vget_lane_f32(sum_2f, 1);
  }  
#endif

  for ( ; i < subband_end; i++)
  {
    sum_abs_real += fabsf(buf->re);
    sum_abs_imag += fabsf(buf->im);
    buf++;
  }
  sum += (sum_abs_real + sum_abs_imag);

  *level += (sum - *level) * SSBUSP;

  return sum;
}

/* ECHOCOMP WEIGHT **************************************************** */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine simply mixes the two signals in fxbuf and adbuf  */
/* into the signal in mixedbuf. The weights used in the mixing is       */
/* fxweight and (1-fxweight).                                           */
/* ******************************************************************** */
void echocomp_weight(COMPLEX32 *  mixedbuf,
                     const COMPLEX32 *  fxbuf,
                     const COMPLEX32 *  adbuf,
                     const float fxweight,
                     const int subband_start,
                     const int subband_end)
{
  int i = subband_start;

#if defined(__ARM_NEON__) && defined(ENV_IOS)
  {
    float32x4_t val_fxw = vdupq_n_f32(fxweight);
    int end = (subband_end-subband_start) &~3;
    end += subband_start;
    for( ; i < end; i+=4 )
    {
      float32x4x2_t val_fx = vld2q_f32((float32_t*)&fxbuf[i]);
      float32x4x2_t val_ad = vld2q_f32((float32_t*)&adbuf[i]);
      float32x4x2_t val_mx;
      float32x4_t re = vsubq_f32(val_fx.val[0], val_ad.val[0]);
      float32x4_t im = vsubq_f32(val_fx.val[1], val_ad.val[1]);
      re = vaddq_f32(val_ad.val[0], vmulq_f32(val_fxw, re));
      im = vaddq_f32(val_ad.val[1], vmulq_f32(val_fxw, im));
      val_mx.val[0] = re;
      val_mx.val[1] = im;
      vst2q_f32((float32_t*)&mixedbuf[i], val_mx);
    }
  }  
#endif

  for ( ; i < subband_end; i++)
  {
    // The usual "(1-a)x + ay" form requires two multiplications;
    // rewriting it as "a(y-x) + x" requires only one.
    mixedbuf[i].re = fxweight * (fxbuf[i].re - adbuf[i].re) + adbuf[i].re;
    mixedbuf[i].im = fxweight * (fxbuf[i].im - adbuf[i].im) + adbuf[i].im;
  }
}

/* f raised to the power of n ***************************************** */
/* Auth.: EHO(GEO)                                                      */
/* This routine will compute f^n, when n is a non-negative integer,     */
/* fast implementation, for instance: x^320, requiers 9 multiplications */
/* Not paricularly efficient for very large exponents!                  */
/* ******************************************************************** */
float powf_to_n(float x, int n)
{
  float y;
  if ( n&0x1 )
  {
    y=x;
  }
  else
  {
    y=1.0f;
  }
  n /= 2;
  while( n!=0 )
  {
    x *= x;
    if ( n&0x1 )
    {
      y*=x;
    }
    n /= 2;
  }
  return y;
}

/* ECHOCOMP FIND DECAY ************************************************ */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine computes the dynamic adation weights a[]. This   */
/* is split in two parts: one early part common for all subands called  */
/* attack_taps, and one late part different for each subband_group      */
/* called decay. We only need one value for the later part in a subband */
/* as it is assumed to be exponential. When both these have been found  */
/* a scaling factor is computed to ensure the norm of a[] is constant.  */
/* ******************************************************************** */
void echocomp_find_decay(ECHOCOMP_PTR pEchocomp)
{
  int i, ch;
  float decay_16;
  float decay_1est, decay_2est, decay_4est, decay_8est, decay_16est;
  float new_attack_gain, new_decay;

  int ix = pEchocomp->decAttData.decay_ix;

  if (pEchocomp->decAttData.decay_compute != 1 && ix != -1)
  {
    float early_sum = pEchocomp->decAttData.decay_early_power;
    float late_sum  = pEchocomp->decAttData.decay_late_power;

    for (ch = 0; ch < pEchocomp->nChannels; ch++)
    {
      ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];
      COMPLEX32 *fixed = channel->subband[ix].fixed;
      int flen = channel->subband[ix].filtlen;

      for (i = 0; i < flen; i++)
      {
        float real_tap = fixed->re;
        float imag_tap = fixed->im;
        float power_of_tap = real_tap * real_tap + imag_tap * imag_tap;

        if ( power_of_tap > EPSILON )
        {
          channel->filter_power[i] += sqrtf(power_of_tap);
        }
        if ( i >= DECAY_START_TAP && i < DECAY_START_TAP + DECAY_SUBBAND_GROUPLENGTH )
        {
          early_sum += power_of_tap;
        }
        else if ( i >= DECAY_START_TAP + DECAY_SUBBAND_GROUPLENGTH && i < DECAY_START_TAP + 2*DECAY_SUBBAND_GROUPLENGTH )
        {
          late_sum += power_of_tap;
        }
        fixed++;
      }
    }
    pEchocomp->decAttData.decay_early_power = early_sum;
    pEchocomp->decAttData.decay_late_power  = late_sum;
  }
  else if ( ix != -1 && pEchocomp->decAttData.decay_compute == 1 )
  {
    float denom;

    if (pEchocomp->decAttData.decay_early_power > EPSILON)
    {
      decay_16 = pEchocomp->decAttData.decay_late_power / pEchocomp->decAttData.decay_early_power;

      if ( decay_16 < DECAY_16_MIN )
        decay_16 = DECAY_16_MIN; /* DECAY_MIN; */
      else if ( decay_16 > DECAY_16_MAX )
        decay_16 = DECAY_16_MAX; /* DECAY_MAX; */

      decay_8est = 1.0f / sqrtf(decay_16);   /* approximately decay_16 ^(-1/2 ) */
      decay_4est = 1.0f / sqrtf(decay_8est); /* approximately decay_16 ^( 1/4 ) */
      decay_2est = 1.0f / sqrtf(decay_4est); /* approximately decay_16 ^(-1/8 ) */
      decay_1est = 1.0f / sqrtf(decay_2est); /* approximately decay_16 ^( 1/16), our first guess for the decay */

      /* will now compute better value for decay using Newtons method: */
      /* x[n] = x[n+1]- f(x[n])/f'(x[n])   where   f(t) = x^16 - decay_16 */
      /* we do only one iteration which amounts to */
      /* new_decay  =  decay_1est  -  (decay_1est^16 - decay_16) / (16 * decay_1est^15) */
      /* we start off with computing decay_1est^16 and 16 * decay_1est^15 */
      decay_2est   = decay_1est * decay_1est;
      denom        = decay_1est * 16.0f;
      decay_4est   = decay_2est * decay_2est;
      denom       *= decay_2est;
      decay_8est   = decay_4est * decay_4est;
      denom       *= decay_4est;
      decay_16est  = decay_8est * decay_8est;     /*      decay_1est^16 */
      denom       *= decay_8est;                  /* 16 * decay_1est^15 */
      new_decay    = decay_1est - ((decay_16est - decay_16) / denom);
    }
    else
    {
      new_decay = pEchocomp->channel[0].subband[ix].decay_default;
    }

    /* now its time to ensure sum(mu[time]) = filtlen to guaranty maximum convergence speed while keeping stability.
    The first few mu ([ 0 - DECAY_START_SUBBAND ]) are the same shape for all subbands (except an individual gain)
    while the rest mu([ DECAY_START_SUBBAND - filtlen]) is shaped with different decay for each DECAY_SUBBAND_GROUP */

    for (ch = 0; ch < pEchocomp->nChannels; ch++)
    {
      ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];
      int filtlen = channel->subband[ix].filtlen;
      int number_of_attack_taps = channel->number_of_used_attack_taps_new;
      float this_weight = channel->attack_taps_new[number_of_attack_taps-1];
      float num = (float)filtlen;
      denom = channel->decay_attack_taps_sum_new;

      if ( 1 ) /*efficient version */
      {
        denom += this_weight / (new_decay-1) * (powf_to_n(new_decay, (filtlen-number_of_attack_taps+1))-new_decay);
      }
      //else   /* dumb version */
      //{
      //  int j;
      //  for ( j=number_of_attack_taps; j<filtlen; j++ )
      //  {
      //    this_weight *= new_decay;
      //    denom += this_weight;
      //  }
      //}
      new_attack_gain = num / denom;

      for ( i = ix; i > ix-DECAY_SUBBAND_GROUPSIZE; i-- )
      {
        ECHOCOMP_SUBBAND *subband = &channel->subband[i];

        subband->decay_attack_gain = new_attack_gain;
        subband->decay = new_decay;
        subband->weights_need_recalculate = 1;
      }
    }

    pEchocomp->decAttData.decay_early_power = 0.0f;
    pEchocomp->decAttData.decay_late_power  = 0.0f;
    /* Here decay_compute=1 and ix!= -1 */
  }
  else  /* pEchocomp->decay_compute == 1 && pEchocomp->decay_ix == -1 */
  {
    int last_subband_computed = ((pEchocomp->decAttData.subused_finddecay / DECAY_SUBBAND_GROUPSIZE) * DECAY_SUBBAND_GROUPSIZE) - 1;

    for (ch = 0; ch < pEchocomp->nChannels; ch++)
    {
      ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];
      float largest_tap_value = 0.0f;
      int largest_tap_index = 0;
      float attack_taps_sum;
      float new_attack_tap, this_tap, next_tap;
      float minimum_value_accepted;
      int number_of_used_attack_taps_new; //for loop optimization

      new_decay = channel->subband[last_subband_computed].decay;
      new_attack_gain = channel->subband[last_subband_computed].decay_attack_gain;

      for (i=last_subband_computed+1; i<SUBUSED; i++)
      {
        channel->subband[i].decay = new_decay;
        channel->subband[i].decay_attack_gain = new_attack_gain;
      }

      /* store old values, and find maximum tap in impulse response */
      channel->number_of_used_attack_taps_old = channel->number_of_used_attack_taps_new;
      channel->decay_attack_taps_sum_old = channel->decay_attack_taps_sum_new;

      for ( i=0; i<FILTLEN_MAX; i++ )
      {
        new_attack_tap = channel->filter_power[i];
        if ( new_attack_tap < EPSILON )
        {
          new_attack_tap = EPSILON;
        }
        channel->filter_power[i] = 0.0f;

        if ( i < NUMBER_OF_ATTACK_TAPS )
        {
          channel->attack_taps_old[i] = channel->attack_taps_new[i];
          channel->attack_taps_new[i] = new_attack_tap;
        }
        if ( largest_tap_value < new_attack_tap )
        {
          largest_tap_value = new_attack_tap;
          largest_tap_index = i;
        }
      }
      /* if delay is unreasonably high or low set attack_taps_new to default value */
      if ( largest_tap_index < (DECAY_LARGEST_TAP_AT_MIN-1) || (DECAY_LARGEST_TAP_AT_MAX-1) < largest_tap_index )
      {
        channel->attack_taps_new[0] = 0.2000f;
        channel->attack_taps_new[1] = 0.5000f;
        channel->attack_taps_new[2] = 3.0000f;
        channel->attack_taps_new[3] = 3.0000f;
        channel->attack_taps_new[4] = 3.0000f;
        channel->attack_taps_new[5] = 2.8125f;
        channel->attack_taps_new[6] = 2.6367f;
        channel->attack_taps_new[7] = 2.4719f;
        channel->attack_taps_new[8] = 2.3174f;
        channel->attack_taps_new[9] = 2.1725f;
        largest_tap_value = 0.0f;
        largest_tap_index = 0;

        for ( i=0; i<NUMBER_OF_ATTACK_TAPS; i++ )
        {
          if ( largest_tap_value < channel->attack_taps_new[i] )
          {
            largest_tap_value = channel->attack_taps_new[i];
            largest_tap_index = i;
          }
        }
      }
      /* ensure not to large differende between cosecutive attack_taps_new, and compte sum of them */
      attack_taps_sum = channel->attack_taps_new[largest_tap_index];
      channel->number_of_used_attack_taps_new = (short)(1 + largest_tap_index + DECAY_PRE_ROOM_DOMINACE_TIME);

      for ( i=largest_tap_index; i>0; i-- )
      {
        this_tap = channel->attack_taps_new[i-1];
        next_tap = channel->attack_taps_new[i];
        minimum_value_accepted = MAX_ADAPT_TAPWEIGHT_CHANGE * next_tap;
        if ( minimum_value_accepted < MAX_ADAPT_TAPWEIGHT_CHANGE_TOTAL * largest_tap_value )
        {
          minimum_value_accepted = MAX_ADAPT_TAPWEIGHT_CHANGE_TOTAL * largest_tap_value;
        }
        if ( this_tap < minimum_value_accepted )
        {
          this_tap = minimum_value_accepted;
        }
        else if ( this_tap > next_tap )
        {
          this_tap = next_tap;
        }
        attack_taps_sum += this_tap;
        channel->attack_taps_new[i-1] = this_tap;
      }

      number_of_used_attack_taps_new = channel->number_of_used_attack_taps_new;

      for ( i=largest_tap_index; i < number_of_used_attack_taps_new-1; i++ )
      {
        this_tap = channel->attack_taps_new[i];
        next_tap =channel->attack_taps_new[i+1];
        minimum_value_accepted = MAX_ADAPT_TAPWEIGHT_CHANGE * this_tap;
        if ( minimum_value_accepted<(MAX_ADAPT_TAPWEIGHT_CHANGE_TOTAL * largest_tap_value) )
        {
          minimum_value_accepted = (MAX_ADAPT_TAPWEIGHT_CHANGE_TOTAL * largest_tap_value);
        }
        if ( next_tap < minimum_value_accepted )
        {
          next_tap = minimum_value_accepted;
        }
        attack_taps_sum += next_tap;
        channel->attack_taps_new[i+1] = next_tap;
      }
      channel->decay_attack_taps_sum_new = attack_taps_sum;

      for (i=0; i<SUBUSED; i++)
      {
        channel->subband[i].weights_need_recalculate = 1;
      }
    }
    /* Here we know dacay_compute=1 and ix=-1, at exit from find_decay they will both be zero  */
  }

  if ( ( (ix+1)%DECAY_SUBBAND_GROUPSIZE != 0 ) && ( pEchocomp->decAttData.decay_compute == 0 ) )
  {
    if (pEchocomp->decAttData.decay_ix >= 0)
      pEchocomp->channel[0].subband[pEchocomp->decAttData.decay_ix].weights_need_recalculate = 1;

    pEchocomp->decAttData.decay_ix++;
  }
  else if ( pEchocomp->decAttData.decay_compute == 1 &&
            (ix+1) != (pEchocomp->decAttData.subused_finddecay / DECAY_SUBBAND_GROUPSIZE) * DECAY_SUBBAND_GROUPSIZE )
  {
    pEchocomp->decAttData.decay_compute = 0;

    if (pEchocomp->decAttData.decay_ix >= 0)
      pEchocomp->channel[0].subband[pEchocomp->decAttData.decay_ix].weights_need_recalculate = 1;

    pEchocomp->decAttData.decay_ix++;
  }
  else if ( pEchocomp->decAttData.decay_compute == 1 )
  {
    pEchocomp->decAttData.decay_ix = -1;
    for (i=0; i<SUBUSED; i++)
    {
      pEchocomp->channel[0].subband[i].weights_need_recalculate = 1;
    }
  }
  else
  {
    pEchocomp->decAttData.decay_compute = 1; /* must have come from here , ix == -1 (mod DECAY_SUBBAND_GROUPSIZE) */
  }
}

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
                           const int subband_end)
{
  const float DeltaTC = 0.97f;
  const float DeltaGain = 0.2f;
  float MinDelta = pEchocomp->minDelta; //const float MinDelta = 1e-5f;
  float MaxDelta = pEchocomp->maxDelta; //const float MaxDelta = 5e-4f;
  int i;

  float DeltaTC_excgain = DeltaTC * excgain;
  float DeltaGain_excgain = DeltaGain * (1.0f - DeltaTC) * excgain;

  for (i = subband_start ; i < subband_end; i++)
  {
    ECHOCOMP_SUBBAND *subband = &pEchocomp->channel[0].subband[i];
    float  sum = subband->est_err.re * subband->est_err.re + subband->est_err.im * subband->est_err.im;
    subband->delta = subband->delta*DeltaTC_excgain + DeltaGain_excgain * sum; 
    subband->delta = max(min(subband->delta, MaxDelta), MinDelta);
  }
  
}

