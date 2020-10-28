/****************************************************************************
 *                        AUDIO - ECHO MODULE
 * ---------------------------------------------------------------------------
 *            (C) Copyright Tandberg Telecom ASA
 * ===========================================================================
 * Revision: $Revision: 1.16 $
 * Author:   EHO(GEO)
 * Date:
 * Desc:     APA echo estimation model.
 *
 * Note:
 *---------------------------------------------------------------------------
 * Revisions:
 *---------------------------------------------------------------------------
 * dd-mm-yyyy <Auth>
 *            <desc>
 *****************************************************************************/
#include "echocomp.h"
#include "echocomp_priv.h"
#include <math.h>
#define fastdiv(a,b) ((a)/(b))
#define fasterdiv(a,b) ((a)/(b))

#define nDBG_SUBBAND_MAXTAPS
#ifdef  DBG_SUBBAND_MAXTAPS
#endif

#ifdef __ARM_NEON__
#include <arm_neon.h>
//#pragma thumb           // for compatibility with gcc
#endif

static void apa_recalculate_weights(const ECHOCOMP_PTR pEchocomp,
                                    const int i /* subband index */)
{
  ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[0]; // this function only used in mono mode
  ECHOCOMP_SUBBAND *subband = &channel->subband[i];
  const int filtlen = subband->filtlen;
  float *adaptWeights = subband->adaptWeights;
  const float subbandDecay = subband->decay;
  const float attackGain = subband->decay_attack_gain;
  float expDecay = 1.0f;
  int attack_length;
  const float *attackTaps;
  int j;

  if ( i < (pEchocomp->decAttData.decay_ix / DECAY_SUBBAND_GROUPSIZE) * DECAY_SUBBAND_GROUPSIZE )
  {
    attackTaps    = channel->attack_taps_new;
    attack_length = channel->number_of_used_attack_taps_new;
  }
  else
  {
    attackTaps    = channel->attack_taps_old;
    attack_length = channel->number_of_used_attack_taps_old;
  }
  for (j=0; j<attack_length; j++)
  {
    adaptWeights[j] = attackGain * attackTaps[j];
  }
  if (attack_length == 0)
  {
    adaptWeights[0] = expDecay;
    attack_length = 1;
  }
  else
  {
    expDecay = adaptWeights[attack_length-1];
  }
  for (j=attack_length; j<filtlen; j++)
  {
    expDecay *= subbandDecay;
    adaptWeights[j] = expDecay;
  }
  for (j=filtlen; j<((filtlen+3)& ~3); j++)
  {
    adaptWeights[j] = 0.0f;
  }
  subband->weights_need_recalculate = 0;
}

/* ECHOCOMP FILTER APA ************************************************ */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine calculates the DOT-product of filtertaps[] (z(k) */
/* and pLsprocess->delaylines[] and adds a correction given by mu*beta1*r12  */
/* and places the answer in est for each subband. Thus est[] contains   */
/* the estimated echo signal to reach the microphone.                   */
/*                                                                      */
/* Returns: The highest maxTapAdapt^2 */
/* ******************************************************************** */

float apa_filter_APA2(const ECHOCOMP_PTR pEchocomp,
                      COMPLEX32 * est,
                      COMPLEX32 * est2,
                      const int subband_start,
                      const int subband_end)
{
  float maxTapAdapt = 0;
  float absTapAdapt;
  int i, j;

  for (i = subband_start; i < subband_end; i++)
  {
    ECHOCOMP_SUBBAND *subband = &pEchocomp->channel[0].subband[i];

    float est_real  = 0.0f;
    float est_imag  = 0.0f;
    float est_real2 = 0.0f;
    float est_imag2 = 0.0f;

#if 0//defined(__ARM_NEON__) && defined(ENV_DARWIN)	// should work on Android, but leaving for now
    const int filtlen = (subband->filtlen + 3) & ~3;

    float32x4_t Max  = vdupq_n_f32(0.0f);
    float32x4_t E1re = vdupq_n_f32(0.0f);
    float32x4_t E1im = vdupq_n_f32(0.0f);
    float32x4_t E2re = vdupq_n_f32(0.0f);
    float32x4_t E2im = vdupq_n_f32(0.0f);
    float32x2_t v2; // temporary
    float32_t maxhi, maxlo;
  
    for (j=0; j<filtlen; j += 4)
    {
      // load 4 entries from each set of taps
      float32x4x2_t DL = vld2q_f32((float32_t *)&subband->delayline[j]);
      float32x4x2_t AD = vld2q_f32((float32_t *)&subband->adapt[j]);
      float32x4x2_t FX = vld2q_f32((float32_t *)&subband->fixed[j]);

      float32x4_t DLre = DL.val[0], DLim = DL.val[1];
      float32x4_t ADre = AD.val[0], ADim = AD.val[1];
      float32x4_t FXre = FX.val[0], FXim = FX.val[1];
      
      Max = vmaxq_f32(Max, vaddq_f32(vmulq_f32(ADre, ADre), vmulq_f32(ADim, ADim)));

      E1re = vmlsq_f32(vmlaq_f32(E1re, ADre, DLre), ADim, DLim);
      E1im = vmlaq_f32(vmlaq_f32(E1im, ADre, DLim), ADim, DLre);
      E2re = vmlsq_f32(vmlaq_f32(E2re, FXre, DLre), FXim, DLim);
      E2im = vmlaq_f32(vmlaq_f32(E2im, FXre, DLim), FXim, DLre);
    }
  
    // collapse accumulators
    v2 = vpadd_f32(vget_high_f32(E1re), vget_low_f32(E1re));
    est_real += vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);

    v2 = vpadd_f32(vget_high_f32(E1im), vget_low_f32(E1im));
    est_imag += vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);

    v2 = vpadd_f32(vget_high_f32(E2re), vget_low_f32(E2re));
    est_real2 += vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);

    v2 = vpadd_f32(vget_high_f32(E2im), vget_low_f32(E2im));
    est_imag2 += vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);

    v2 = vpmax_f32(vget_high_f32(Max), vget_low_f32(Max));
  
    maxhi = vget_lane_f32(v2, 1);
    maxlo = vget_lane_f32(v2, 0);
  
    absTapAdapt = (maxhi > maxlo) ? maxhi : maxlo;
  
    if(absTapAdapt > maxTapAdapt)
    {
      maxTapAdapt = absTapAdapt;
    }
#else
    const int filtlen = subband->filtlen;

    const COMPLEX32 *DL = subband->delayline;
    const COMPLEX32 *AD = subband->adapt;
    const COMPLEX32 *FX = subband->fixed;

    absTapAdapt = 0.0;

    for (j=0; j<filtlen; j++)
    {
      absTapAdapt += AD->re * AD->re + AD->im * AD->im;

      est_real  += AD->re * DL->re - AD->im * DL->im;
      est_imag  += AD->re * DL->im + AD->im * DL->re;
      est_real2 += FX->re * DL->re - FX->im * DL->im;
      est_imag2 += FX->re * DL->im + FX->im * DL->re;

      AD++; FX++; DL++;
    }

    if(absTapAdapt > maxTapAdapt * filtlen)
    {
      maxTapAdapt = absTapAdapt / filtlen;
    }

#endif
    est[i].re  = est_real + subband->FiltCorr.re;
    est[i].im  = est_imag + subband->FiltCorr.im;
    est2[i].re = est_real2;  
    est2[i].im = est_imag2;
  }
  return maxTapAdapt;
}

/* APA CALC R11 AND R12 *********************************************** */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine calculates the componenets of a matrix R, where  */
/* R=(X'*A*X) in Matlab notation, so ' means hermitian transpose. A is  */
/* the diagonal matrix which has the dynamic adaption weights on its    */
/* diagonal. Note that R will be hermitic so it is enough to calculate  */
/* the upper half of the matrix, for instance r21=conj(12). If it is    */
/* desirable to save some CPU it is possible to do this recursivly, but */
/* care must be taken to recalculate from scratch when A is updated.    */
/* This should thrn be done for one subaband at a time in find_decay    */
/* after a new decay has been computed for one subbandgroup.            */
/* ******************************************************************** */
void apa_calcR11_and_R12(ECHOCOMP_PTR pEchocomp,
                         const int subband_start,
                         const int subband_end)
{
  ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[0]; // this function only used in mono mode
  int i, j;

  for (i = subband_start; i < subband_end; i++)
  {
    ECHOCOMP_SUBBAND *subband = &channel->subband[i];
    const int filtlen = subband->filtlen;
  
    float r11    = 0.0f;
    float r22    = 0.0f;
    float r12_re = 0.0f;
    float r12_im = 0.0f;

    if (subband->weights_need_recalculate)
    {
      apa_recalculate_weights(pEchocomp, i);
    }

#if defined(__ARM_NEON__) && defined(ENV_DARWIN)	// should work on Android, but leaving for now
    {
      const int filtlen = (subband->filtlen + 3) & ~3;  // hides previous declaration
      
      float32x4_t R11   = vdupq_n_f32(0.0f);
      float32x4_t R22   = vdupq_n_f32(0.0f);
      float32x4_t R12re = vdupq_n_f32(0.0f);
      float32x4_t R12im = vdupq_n_f32(0.0f);

      float32x4x2_t L1  = vld2q_f32((float32_t*)&subband->delayline[0]);   // preload first 4 entries
      float32x4_t L1re  = L1.val[0], L1im  = L1.val[1];
      float32x2_t v2; // temporary

      for (j=0; j<filtlen; j += 4)
      {
        // load 4 entries from each set of taps
        float32x4x2_t Next = vld2q_f32((float32_t*)&subband->delayline[j+4]);
        float32x4_t Nextre = Next.val[0], Nextim  = Next.val[1];
        float32x4_t L2re   = vextq_f32(L1re, Nextre, 1);
        float32x4_t L2im   = vextq_f32(L1im, Nextim, 1);
        float32x4_t A      = vld1q_f32(&subband->adaptWeights[j]);

        R11   = vmlaq_f32(R11,   A, vaddq_f32(vmulq_f32(L1re, L1re), vmulq_f32(L1im, L1im)));
        R22   = vmlaq_f32(R22,   A, vaddq_f32(vmulq_f32(L2re, L2re), vmulq_f32(L2im, L2im)));
        R12re = vmlaq_f32(R12re, A, vaddq_f32(vmulq_f32(L1re, L2re), vmulq_f32(L1im, L2im)));
        R12im = vmlaq_f32(R12im, A, vsubq_f32(vmulq_f32(L1im, L2re), vmulq_f32(L1re, L2im)));

        L1re = Nextre; L1im = Nextim;
      }
      
      // collapse accumulators
      v2 = vpadd_f32(vget_high_f32(R11), vget_low_f32(R11));
      r11 = vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);
      
      v2 = vpadd_f32(vget_high_f32(R22), vget_low_f32(R22));
      r22 = vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);
      
      v2 = vpadd_f32(vget_high_f32(R12re), vget_low_f32(R12re));
      r12_re = vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);
      
      v2 = vpadd_f32(vget_high_f32(R12im), vget_low_f32(R12im));
      r12_im = vget_lane_f32(v2, 1) + vget_lane_f32(v2, 0);
    }
#else
    for (j=0; j<filtlen; j++)
    {
      const float a = subband->adaptWeights[j];

      const float r1 = subband->delayline[j].re;
      const float i1 = subband->delayline[j].im;
      const float r2 = subband->delayline[j+1].re;
      const float i2 = subband->delayline[j+1].im;

      r11    += a * (r1*r1 + i1*i1);
      r22    += a * (r2*r2 + i2*i2);
      r12_re += a * (r1*r2 + i1*i2);
      r12_im += a * (i1*r2 - r1*i2);
    }   
#endif
    subband->r11    = r11 + subband->delta;
    subband->r22    = r22 + subband->delta;
    subband->r12.re = r12_re;
    subband->r12.im = r12_im;
  }
}

/* ECHOCOMP CALCULATE FILTER CORRECTION FOR APA *********************** */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: This routine calculates the nessecary correction factor to be */
/* added in the filtering routine when the filters in the APA-algorithm */
/* are transformed to achieve faster performance in the filter adaption */
/* routine for APA. This correcton factor must then be added to the     */
/* result of a normal filtering (as in a straight forward APA or LMS)   */
/* ******************************************************************** */
void apa_calcFiltCorr(const ECHOCOMP_PTR pEchocomp,
                      const int subband_start,
                      const int subband_end)
{
  float mu = pEchocomp->mu;
  int i;

  for (i = subband_start; i < subband_end; i++)
  {
    ECHOCOMP_SUBBAND *subband = &pEchocomp->channel[0].subband[i];
    COMPLEX32 r12;

    r12.re = subband->r12.re;
    r12.im = subband->r12.im;

    subband->FiltCorr.re = mu * (subband->beta_1.re * r12.re -
                                 subband->beta_1.im * r12.im);
    subband->FiltCorr.im = mu * (subband->beta_1.im * r12.re +
                                 subband->beta_1.re * r12.im);
  }
}

/* ECHOCOMP ADAPTION APA ********************************************** */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: Updates filter according to APA.                              */
/* ******************************************************************** */
void apa_adapt_APA(ECHOCOMP_PTR pEchocomp,
                   const COMPLEX32 *  adbuf,
                   const COMPLEX32 *  fxbuf,
                   const int subband_start,
                   const int subband_end)
{
  ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[0]; // this function only used in mono mode
  COMPLEX32 r12, est_err2, beta_2, beta_12;
  float mu = pEchocomp->mu;
  float r11, r22, den;
  int i, j;

  if ((int) (mu * 100.0f) == 100)
  {
    for (i = subband_start; i < subband_end; i++)
    {
      ECHOCOMP_SUBBAND *subband = &channel->subband[i];
      const int filtlen = subband->filtlen;

      r11 = subband->r11;
      r22 = subband->r22;
      r12.re = subband->r12.re;
      r12.im = subband->r12.im;

      den = 1.0f / (r11*r22 - r12.re*r12.re - r12.im*r12.im);

      subband->est_err.re = adbuf[i].re;
      subband->est_err.im = adbuf[i].im;

      beta_2.re = -(r12.re * subband->est_err.re + r12.im * subband->est_err.im) * den;
      beta_2.im = -(r12.re * subband->est_err.im - r12.im * subband->est_err.re) * den;

      beta_12.re = subband->beta_1.re + beta_2.re;
      beta_12.im = subband->beta_1.im + beta_2.im;

      if (subband->weights_need_recalculate)
      {
        apa_recalculate_weights(pEchocomp, i);
      }
#if defined(__ARM_NEON__) && defined(ENV_DARWIN)	// should work on Android, but leaving for now
      {
        const int filtlen = (subband->filtlen + 3) & ~3;  // hides previous declaration
        
        for (j=0; j<filtlen; j += 4)
        {
          // load 4 entries from each set of taps
          float32x4x2_t DL = vld2q_f32((float32_t *)&subband->delayline[j+1]);
          float32x4x2_t AD = vld2q_f32((float32_t *)&subband->adapt[j]);

          float32x4_t DLre = DL.val[0], DLim = DL.val[1];
          float32x4_t ADre = AD.val[0], ADim = AD.val[1];
          float32x4_t A    = vld1q_f32(&subband->adaptWeights[j]);
          
          AD.val[0] = vaddq_f32(ADre, vmulq_f32(A, vaddq_f32(vmulq_n_f32(DLre, beta_12.re), vmulq_n_f32(DLim, beta_12.im))));
          AD.val[1] = vaddq_f32(ADim, vmulq_f32(A, vsubq_f32(vmulq_n_f32(DLre, beta_12.im), vmulq_n_f32(DLim, beta_12.re))));
          
          vst2q_f32((float32_t *)&subband->adapt[j], AD);
        }
      }
#else
      for (j=0; j<filtlen; j++)
      {
        const float a = subband->adaptWeights[j];
        const float dl_real = subband->delayline[j+1].re;
        const float dl_imag = subband->delayline[j+1].im;

        subband->adapt[j].re += a * (beta_12.re * dl_real + beta_12.im * dl_imag);
        subband->adapt[j].im += a * (beta_12.im * dl_real - beta_12.re * dl_imag);
      }
#endif
      subband->beta_1.re = subband->r22 * subband->est_err.re * den;
      subband->beta_1.im = subband->r22 * subband->est_err.im * den;

      subband->beta_fixed.re = subband->r22 * fxbuf[i].re * den;
      subband->beta_fixed.im = subband->r22 * fxbuf[i].im * den;
    }
  }
  else
  {
    for (i = subband_start; i < subband_end; i++)
    {
      ECHOCOMP_SUBBAND *subband = &channel->subband[i];
      const int filtlen = subband->filtlen;

      r11 = subband->r11;
      r22 = subband->r22;
      r12.re = subband->r12.re;
      r12.im = subband->r12.im;

      den = 1.0f /(r11*r22 - r12.re*r12.re - r12.im*r12.im);

      est_err2.re = (1.0f - mu) * subband->est_err.re;
      est_err2.im = (1.0f - mu) * subband->est_err.im;

      subband->est_err.re =  adbuf[i].re;
      subband->est_err.im =  adbuf[i].im;

      beta_2.re = r11 * est_err2.re - (r12.re*subband->est_err.re + r12.im*subband->est_err.im);
      beta_2.im = r11 * est_err2.im - (r12.re*subband->est_err.im - r12.im*subband->est_err.re);

      beta_12.re = (subband->beta_1.re + beta_2.re) * (mu * den);
      beta_12.im = (subband->beta_1.im + beta_2.im) * (mu * den);

      if (subband->weights_need_recalculate)
      {
        apa_recalculate_weights(pEchocomp, i);
      }

      for (j=0; j<filtlen; j++)
      {
        const float a = subband->adaptWeights[j];
        const float dl_real = subband->delayline[j+1].re;
        const float dl_imag = subband->delayline[j+1].im;

        subband->adapt[j].re += a * (beta_12.re * dl_real + beta_12.im * dl_imag);
        subband->adapt[j].im += a * (beta_12.im * dl_real - beta_12.re * dl_imag);
      }

      subband->beta_1.re = (subband->r22 * subband->est_err.re -
                            (subband->r12.re * est_err2.re - subband->r12.im * est_err2.im)) * den;
      subband->beta_1.im = (subband->r22 * subband->est_err.im -
                            (subband->r12.re * est_err2.im + subband->r12.im * est_err2.re)) * den;

      /* The following is not correct. To make it correct we need to store
         previous fxbuf and do as in the computation of beta_1 above */
      subband->beta_fixed.re = subband->r22 * fxbuf[i].re * den;
      subband->beta_fixed.im = subband->r22 * fxbuf[i].im * den;
    }
  }
}

/* ECHOCOMP FINDBEST APA*********************************************** */
/* Auth.: EHO(GEO)                                                      */
/* Desc.: Decides is fixed or adaptive filter is to be used, or restets */
/* adaptive filter to zero. If adapt is sufficiently good, fixed filter */
/* is replaced by it, or if fixed is better, adapt is replaced by it.   */
/* Note that a translation of the filter is needed when copying between */
/* the adaptive filter and fixed filter since the fixed filter must     */
/* contain the true impulse response while the adative is transformed.  */
/* Returns:                                                             */
/*  -3 if adapt is reset to zero                                        */
/*   1 if fixed is set to adapt                                         */
/*  -1 if adapt is reset to fixed                                       */
/*   0 if processing was normal, (no copying of filters or resetting)   */
/* ******************************************************************** */
short apa_findbest_APA(ECHOCOMP_PTR  pEchocomp,
                       int k,
                       const int subband_start,
                       const int subband_end)
{
    ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[0]; // this function only used in mono mode
    short didcopy = 0;
    int i,j;
    float Threshold;
    float expDecay;
    float attackGain;
    float subbandDecay;

    if (pEchocomp->adaptlevel[k] > 1.3333333333f * pEchocomp->level[k])
    {
        /* resets adaptive filter to zero    */
        /* printf("\nAdapt set to zero \n"); */
        for (i = subband_start; i < subband_end; i++)
        {
            ECHOCOMP_SUBBAND *subband = &channel->subband[i];
            const int filtlen = subband->filtlen;
          
            for (j=0; j<filtlen; j++)
            {
                subband->adapt[j].re = 0.0f;
                subband->adapt[j].im = 0.0f;
            }
            subband->beta_1.re = 0.0f;
            subband->beta_1.im = 0.0f;
        }
        pEchocomp->adaptlevel[k] = pEchocomp->level[k];
        didcopy = -3;
    }
    else if (pEchocomp->fixedlevel[k] >  pEchocomp->adaptlevel[k])
    {
        /* adaptive filter is best */
        /* reduce fixed weight     */
        if (pEchocomp->fixedweight[k] > 0.125f)
        {
            pEchocomp->fixedweight[k] -= 0.125f;
        }
        else
        {
            pEchocomp->fixedweight[k]  = 0.0f;
        }

        Threshold = 1.0f - 0.5f * pEchocomp->adaptlevel[k] / pEchocomp->level[k];
        if (Threshold > 0.875f)
            Threshold = 0.875f;

        if ( pEchocomp->adaptlevel[k] < pEchocomp->fixedlevel[k] * Threshold)
           {
            /* copy adaptive filter to fixed filter */
            /* printf("\nFixed set to Adapt\n");    */

                for (i = subband_start; i < subband_end; i++)
                {
                    ECHOCOMP_SUBBAND *subband = &channel->subband[i];
                    const int filtlen = subband->filtlen;
                    float adaptWeights[FILTLEN_MAX];
                    int attack_length;
                    const float *attackTaps;

                    subbandDecay = subband->decay;
                    expDecay = pEchocomp->mu;
                    attackGain = expDecay * subband->decay_attack_gain;

                    if (i < (pEchocomp->decAttData.decay_ix / DECAY_SUBBAND_GROUPSIZE) * DECAY_SUBBAND_GROUPSIZE)
                    {
                        attackTaps    = channel->attack_taps_new;
                        attack_length = channel->number_of_used_attack_taps_new;
                    }
                    else
                    {
                        attackTaps    = channel->attack_taps_old;
                        attack_length = channel->number_of_used_attack_taps_old;
                    }
                    for (j=0; j<filtlen; j++)
                    {
#ifdef PRESERVE_ORIGINAL_BUGS
                        if (j*2<attack_length)  // looks wrong, but matches original
                        {
                            expDecay = attackGain * (attackTaps[j*2]);  // looks wrong, but matches original
                        }
#else
                        if (j<attack_length)
                        {
                            expDecay = attackGain * (attackTaps[j]);
                        }
#endif
                        else
                        {
                            expDecay *= subbandDecay;
                        }
                        adaptWeights[j] = expDecay;
                    }
                    for (j=0; j<filtlen; j++)
                    {
                        expDecay = adaptWeights[j];
                        subband->fixed[j].re = subband->adapt[j].re + expDecay * (subband->beta_1.re * subband->delayline[j+1].re +
                                                                                  subband->beta_1.im * subband->delayline[j+1].im);
                        subband->fixed[j].im = subband->adapt[j].im + expDecay * (subband->beta_1.im * subband->delayline[j+1].re -
                                                                                  subband->beta_1.re * subband->delayline[j+1].im);
                    }
                }
                pEchocomp->fixedlevel[k] = pEchocomp->adaptlevel[k];
                didcopy = 1;
           }
    }
    else
    {
        /* fixed filter is best  */
        /* increase fixed weight */
        if (pEchocomp->fixedweight[k] < (1.0f - 0.125f))
        {
            pEchocomp->fixedweight[k] += 0.125f;
        }
        else
        {
            pEchocomp->fixedweight[k]  = 1.0f;
        }

        if ( 1.3333333333f*pEchocomp->fixedlevel[k] < pEchocomp->adaptlevel[k] )
        {
            /* if big, difference, copy fixed filter to adaptive filter */
            /* printf("\nAdapt set to Fixed\n");                        */
            for (i = subband_start; i < subband_end; i++)
            {
                ECHOCOMP_SUBBAND *subband = &channel->subband[i];
                const int filtlen = subband->filtlen;
                float adaptWeights[FILTLEN_MAX];
                int attack_length;
                const float *attackTaps;

                subbandDecay = subband->decay;
                expDecay = pEchocomp->mu;
                attackGain = expDecay * subband->decay_attack_gain;
                if (i < (pEchocomp->decAttData.decay_ix / DECAY_SUBBAND_GROUPSIZE) * DECAY_SUBBAND_GROUPSIZE)
                {
                    attackTaps    = channel->attack_taps_new;
                    attack_length = channel->number_of_used_attack_taps_new;
                }
                else
                {
                    attackTaps    = channel->attack_taps_old;
                    attack_length = channel->number_of_used_attack_taps_old;
                }
                for (j=0; j<filtlen; j++)
                {
#ifdef PRESERVE_ORIGINAL_BUGS
                    if (j*2<attack_length)  // looks wrong, but matches original
                    {
                        expDecay = attackGain * (attackTaps[j*2]);  // looks wrong, but matches original
                    }
#else
                    if (j<attack_length)
                    {
                        expDecay = attackGain * (attackTaps[j]);
                    }
#endif
                    else
                    {
                        expDecay *= subbandDecay;
                    }
                    adaptWeights[j] = expDecay;
                }
                for (j=0; j<filtlen; j++)
                {
                    expDecay = adaptWeights[j];
                    subband->adapt[j].re = subband->fixed[j].re - expDecay * (subband->beta_fixed.re * subband->delayline[j+1].re +
                                                                              subband->beta_fixed.im * subband->delayline[j+1].im);
                    subband->adapt[j].im = subband->fixed[j].im - expDecay * (subband->beta_fixed.im * subband->delayline[j+1].re -
                                                                              subband->beta_fixed.re * subband->delayline[j+1].im);
                }
                subband->beta_1.re = 0.0f;
                subband->beta_1.im = 0.0f;
            }

            pEchocomp->adaptlevel[k] = pEchocomp->fixedlevel[k];
            didcopy = -1;
        }
    }
    return didcopy;
}
