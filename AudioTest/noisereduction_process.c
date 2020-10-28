// armcc -c --fpmode=fast --asm --c99 -O3 -Otime -- --vectorize --cpu=Cortex-A9 --fpu=softvfp+vfpv3 --library_interface=aeabi_clib -Ono_memcpy --thumb --diag_warning=optimizations --diag_style=ide --depend_format=unix_escaped --no_depend_system_headers

#include "noisereduction.h"
#include "noisereduction_priv.h"
#include <math.h>
#include <float.h>

#define SUBUSED16 116
#define M_4_PI 4.0f/3.14159265358979323846f
#define fastlog(x) (logf(x))
#define fasterdiv(a,b) ((a)/(b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

#ifdef __ARM_NEON__
#include <arm_neon.h>
//#pragma thumb           // for compatibility with gcc
#endif

/***************************************************************************
* MEASURELEVEL 1 & 2
*  Auth.: IFA
*  Desc.:  These routines takes a complex vector (the output of an fft)
*          and calculates the level of each subband (used in noisereduction).
*          The levels are updated slower than the real level change from last
*          measurement.
*          '1' takes one input vector while '2' takes two. '2' does the
*          same as '1', but if the input vector is a difference of two vectors.
*          The bufPtr(s) is assumed to be a complex vector of length >= FFTSIZE,
*          and the output levelPtr will have length SUBUSED.
**************************************************************************/
float noisereduction_measureLevel1(float speed, float * const  levelPtr, const COMPLEX32 *bufPtr) 
{ 
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
  /*
  Xn+1 = Xn*(3-d*Xn*Xn)/2 converge to 1/sqrt(d)
  sqrt(d) ~= d * X0*(3-d*X0*X0)/2
  *//*
     tmp3 = tmp4 - tmp3 * M_4_PI; tmp4 = tmp4 - tmp3 * speed;
  => tmp4 = tmp4 - (tmp4 - tmp3 * M_4_PI)*speed = tmp4 - tmp4 *speed + tmp3*M_4_PI*speed
  */
  
  int m;
  float32x4_t   M4_4_PI_mul_speed;
  float32x4_t   flt_min_n;
  float32x4_t   speed_n;
  float32x4x2_t tmpd, tmpd2;
  float32x4_t   tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  float32x4_t   tmp11,tmp12,tmp13,tmp14,tmp15,tmp16;

  // Put the first iteration out of loop in order to reduce the stall of "vrsqrteq_f32"

  // interleaved real/imag, vld2q will de-interleave
  tmpd = vld2q_f32((float32_t *)&bufPtr[1]);
  tmp1 = tmpd.val[0];   // real
  tmp2 = tmpd.val[1];   // imag

  tmpd2 = vld2q_f32((float32_t *)&bufPtr[5]);
  tmp11 = tmpd2.val[0];   // real
  tmp12 = tmpd2.val[1];   // imag
 
  tmp1 = vmulq_f32(tmp1, tmp1);
  tmp2 = vmulq_f32(tmp2, tmp2);
  tmp11= vmulq_f32(tmp11,tmp11);
  tmp12= vmulq_f32(tmp12,tmp12);

  M4_4_PI_mul_speed = vdupq_n_f32(M_4_PI*speed);
  flt_min_n = vdupq_n_f32(FLT_MIN);
  speed_n = vdupq_n_f32(speed);
  
  tmp3 = vaddq_f32(tmp1, tmp2);
  tmp13= vaddq_f32(tmp11,tmp12);
  
  tmp3 = vmaxq_f32(flt_min_n, tmp3);   // d
  tmp13= vmaxq_f32(flt_min_n, tmp13);
  
  tmp5 = vrsqrteq_f32(tmp3);           // X0
  tmp15= vrsqrteq_f32(tmp13);   
  
  for (m = 0; m < SUBUSED; m+=8)
  {
    tmp6 = vmulq_f32(tmp5, tmp5);             // X0*X0
    tmp16= vmulq_f32(tmp15,tmp15);

    tmp4 = vld1q_f32(&levelPtr[m]);
    tmp14= vld1q_f32(&levelPtr[m+4]);

    tmpd = vld2q_f32((float32_t *)&bufPtr[m+9]);
    tmp1 = tmpd.val[0];   // real
    tmp2 = tmpd.val[1];   // imag

    tmp6 = vrsqrtsq_f32(tmp3, tmp6);          // (3-d*X0*X0)/2
    tmp16= vrsqrtsq_f32(tmp13,tmp16);

    tmp3 = vmulq_f32(tmp3, tmp5);           // X0*d
    tmp13= vmulq_f32(tmp13,tmp15);
    
    tmpd2 = vld2q_f32((float32_t *)&bufPtr[m+13]);
    tmp11 = tmpd2.val[0];   // real
    tmp12 = tmpd2.val[1];   // imag

    tmp1 = vmulq_f32(tmp1, tmp1);
    tmp2 = vmulq_f32(tmp2, tmp2);

    tmp4 = vmlsq_f32(tmp4, tmp4, speed_n);    // tmp4 - tmp4 * speed
    tmp14= vmlsq_f32(tmp14,tmp14,speed_n);

    tmp11= vmulq_f32(tmp11,tmp11);
    tmp12= vmulq_f32(tmp12,tmp12);   
    
    tmp6 = vmulq_f32(tmp3, tmp6);          // d * X0*(3-d*X0*X0)/2 = sqrt(d)
    tmp16= vmulq_f32(tmp13,tmp16);
    
    tmp3 = vaddq_f32(tmp1, tmp2);
    tmp13= vaddq_f32(tmp11,tmp12);
    
    tmp6 = vmlaq_f32(tmp4, tmp6, M4_4_PI_mul_speed);   //tmp4 + tmp3*M_4_PI*speed
    tmp16= vmlaq_f32(tmp14,tmp16,M4_4_PI_mul_speed);

    tmp3 = vmaxq_f32(flt_min_n, tmp3);
    tmp13= vmaxq_f32(flt_min_n, tmp13);
    
    tmp5 = vrsqrteq_f32(tmp3);
    tmp15= vrsqrteq_f32(tmp13);    
    
    vst1q_f32(&levelPtr[m], tmp6);
    vst1q_f32(&levelPtr[m+4], tmp16);
  }

  // Calculate to get the sum. Firstly, process 16x in the loop, then process 4x, 
  // and then process the remian. In fact, SUBUSED16 is 4x.
  float sum;
  tmp3 = vdupq_n_f32(0.0f);
  for (m = 0; m < (SUBUSED16&~15); m+=16)
  {
    tmpd  = vld2q_f32((float32_t *)&levelPtr[m]);
    tmpd2 = vld2q_f32((float32_t *)&levelPtr[m+8]);
    tmp3 = vaddq_f32(tmp3, tmpd.val[0]);
    tmp3 = vaddq_f32(tmp3, tmpd.val[1]);
    tmp3 = vaddq_f32(tmp3, tmpd2.val[0]);
    tmp3 = vaddq_f32(tmp3, tmpd2.val[1]);
  }
  for (; m < SUBUSED16; m+=4)
  {
    tmp1  = vld1q_f32((float32_t *)&levelPtr[m]);  
    tmp3 = vaddq_f32(tmp3, tmp1);
  }
  
  float32x2_t tmpX = vpadd_f32(vget_high_f32(tmp3), vget_low_f32(tmp3));
  sum = vget_lane_f32(tmpX, 1) + vget_lane_f32(tmpX, 0);
  
#else
	for (int m = 0; m < SUBUSED; m++)
	{
		float tmp1 = bufPtr[m+1].re;
		float tmp2 = bufPtr[m+1].im;
		float tmp3 = tmp1 * tmp1 + tmp2 * tmp2;
		float tmp4 = levelPtr[m];
		tmp3 = max(tmp3, FLT_MIN);
		tmp3 = sqrtf(tmp3);
		tmp3 = tmp4 - tmp3 * M_4_PI;
		tmp4 = tmp4 - tmp3 * speed;
		levelPtr[m] = tmp4;
	}

  float sum = 0;
  for (int m = 0; m < SUBUSED16; m++)
    sum += levelPtr[m];
#endif

  return sum;
}

/* LEVEL MEASUREMENT : see comment above ****************************** */
float noisereduction_measureLevel2(float speed, float * const  levelPtr,  const COMPLEX32 *bufPtrA, const COMPLEX32 *bufPtrB)
{
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
  /*
  Xn+1 = Xn*(3-d*Xn*Xn)/2 converge to 1/sqrt(d) 
  sqrt(d) ~= d * X0*(3-d*X0*X0)/2
  *//*
     tmp3 = tmp4 - tmp3 * M_4_PI; tmp4 = tmp4 - tmp3 * speed;
  => tmp4 = tmp4 - (tmp4 - tmp3 * M_4_PI)*speed = tmp4 - tmp4 *speed + tmp3*M_4_PI*speed
  */
  int m;
  float32x4_t   M4_4_PI_mul_speed;
  float32x4_t   flt_min_n;
  float32x4_t   speed_n; 
  float32x4x2_t tmpd, tmpd2;
  float32x4_t   tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  float32x4_t   tmp11,tmp12,tmp13,tmp14,tmp15,tmp16;

  // Put the first iteration out of loop in order to reduce the stall of "vrsqrteq_f32"

  // interleaved real/imag, vld2q will de-interleave
  tmpd  = vld2q_f32((float32_t *)&bufPtrA[1]);
  tmpd2 = vld2q_f32((float32_t *)&bufPtrB[1]);

  tmp1 = vsubq_f32(tmpd.val[0], tmpd2.val[0]);   // real
  tmp2 = vsubq_f32(tmpd.val[1], tmpd2.val[1]);   // imag
  
  tmpd  = vld2q_f32((float32_t *)&bufPtrA[5]);
  tmpd2 = vld2q_f32((float32_t *)&bufPtrB[5]);

  tmp11 = vsubq_f32(tmpd.val[0], tmpd2.val[0]);   // real
  tmp12 = vsubq_f32(tmpd.val[1], tmpd2.val[1]);   // imag    
  
  tmp1 = vmulq_f32(tmp1, tmp1);
  tmp2 = vmulq_f32(tmp2, tmp2);
  tmp11= vmulq_f32(tmp11,tmp11);
  tmp12= vmulq_f32(tmp12,tmp12);

  M4_4_PI_mul_speed = vdupq_n_f32(M_4_PI*speed);
  flt_min_n = vdupq_n_f32(FLT_MIN);
  speed_n = vdupq_n_f32(speed);
  
  tmp3 = vaddq_f32(tmp1, tmp2);
  tmp13= vaddq_f32(tmp11,tmp12);
  
  tmp3 = vmaxq_f32(flt_min_n, tmp3);   // d
  tmp13= vmaxq_f32(flt_min_n, tmp13);
  
  tmp5 = vrsqrteq_f32(tmp3);           // X0
  tmp15= vrsqrteq_f32(tmp13);   
  
  for (m = 0; m < SUBUSED; m+=8)
  {
    tmp6 = vmulq_f32(tmp5, tmp5);             // X0*X0
    tmp16= vmulq_f32(tmp15,tmp15);

    tmp4 = vld1q_f32(&levelPtr[m]);
    tmp14= vld1q_f32(&levelPtr[m+4]);
    
    tmp4 = vmlsq_f32(tmp4, tmp4, speed_n);    // tmp4 - tmp4 * speed
    tmp14= vmlsq_f32(tmp14,tmp14,speed_n);
    
    tmp6 = vrsqrtsq_f32(tmp3, tmp6);          // (3-d*X0*X0)/2
    tmp16= vrsqrtsq_f32(tmp13,tmp16);

    tmp3 = vmulq_f32(tmp3, tmp5);             // X0 * d
    tmp13= vmulq_f32(tmp13,tmp15);

    tmpd = vld2q_f32((float32_t *)&bufPtrA[m+9]);
    tmpd2= vld2q_f32((float32_t *)&bufPtrB[m+9]);
    tmp1 = vsubq_f32(tmpd.val[0], tmpd2.val[0]);   // real
    tmp2 = vsubq_f32(tmpd.val[1], tmpd2.val[1]);   // imag

    tmp6 = vmulq_f32(tmp3, tmp6);             // d * X0*(3-d*X0*X0)/2 = sqrt(d)
    tmp16= vmulq_f32(tmp13,tmp16);

    tmp1 = vmulq_f32(tmp1, tmp1);
    tmp2 = vmulq_f32(tmp2, tmp2);

    tmpd = vld2q_f32((float32_t *)&bufPtrA[m+13]);
    tmpd2= vld2q_f32((float32_t *)&bufPtrB[m+13]);
    tmp11= vsubq_f32(tmpd.val[0],tmpd2.val[0]);   // real
    tmp12= vsubq_f32(tmpd.val[1],tmpd2.val[1]);   // imag 
    
    tmp11= vmulq_f32(tmp11,tmp11);
    tmp12= vmulq_f32(tmp12,tmp12);
    
    tmp6 = vmlaq_f32(tmp4, tmp6, M4_4_PI_mul_speed);   //tmp4 + tmp3*M_4_PI*speed
    tmp16= vmlaq_f32(tmp14,tmp16,M4_4_PI_mul_speed);   
    
    tmp3 = vaddq_f32(tmp1, tmp2);
    tmp13= vaddq_f32(tmp11,tmp12);
   
    tmp3 = vmaxq_f32(flt_min_n, tmp3);
    tmp13= vmaxq_f32(flt_min_n, tmp13);

    vst1q_f32(&levelPtr[m],  tmp6);
    vst1q_f32(&levelPtr[m+4],tmp16);
    
    tmp5 = vrsqrteq_f32(tmp3);
    tmp15= vrsqrteq_f32(tmp13);
  }

  // Calculate to get the sum. Firstly, process 16x in the loop, then process 4x, 
  // and then process the remian. In fact, SUBUSED16 is 4x.
  float sum;
  tmp3 = vdupq_n_f32(0.0f);
  for (m = 0; m < (SUBUSED16&~15); m+=16)
  {
    tmpd  = vld2q_f32((float32_t *)&levelPtr[m]);
    tmpd2 = vld2q_f32((float32_t *)&levelPtr[m+8]);
    tmp3 = vaddq_f32(tmp3, tmpd.val[0]);
    tmp3 = vaddq_f32(tmp3, tmpd.val[1]);
    tmp3 = vaddq_f32(tmp3, tmpd2.val[0]);
    tmp3 = vaddq_f32(tmp3, tmpd2.val[1]);
  }
  for (; m < SUBUSED16; m+=4)
  {
    tmp1  = vld1q_f32((float32_t *)&levelPtr[m]);  
    tmp3 = vaddq_f32(tmp3, tmp1);
  }
  
  float32x2_t tmpX = vpadd_f32(vget_high_f32(tmp3), vget_low_f32(tmp3));
  sum = vget_lane_f32(tmpX, 1) + vget_lane_f32(tmpX, 0);
  
#else
  for(int m = 0; m < SUBUSED; m++)
  {
    float tmp1 = bufPtrA[m+1].re - bufPtrB[m+1].re;
    float tmp2 = bufPtrA[m+1].im - bufPtrB[m+1].im;
    float tmp3 = tmp1 * tmp1 + tmp2 * tmp2;
    float tmp4 = levelPtr[m];
    tmp3 = max(tmp3, FLT_MIN);
    tmp3 = sqrtf(tmp3);
    tmp3 = tmp4 - tmp3 * M_4_PI;
    tmp4 = tmp4 - tmp3 * speed;
    levelPtr[m] = tmp4;
  }
  float sum = 0;
  for(int m = 0; m < SUBUSED16; m++)     
    sum += levelPtr[m];
#endif

  return sum;
}

/***************************************************************************
* ESTIMATENOISE
*   Auth.: IFA
*   Desc.: In the subbands where the total level is lower than the
*          estimated noise level, the est. noise level i decreased by an amount
*          given by decspeed. Likewise, in the subbands where the total level is
*          higher than the est. noise level; est. noise level is increased.
*          decspeed and incspeed will (normally) vary depending on whether noise
*          is detected or not. (the output of the 'checknoise' routine.
************************************************************************ */
float noisereduction_estimatenoise(float decSpeed, float incSpeed, const float *const  totlev, float * noilev)
{
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
  int m;
  float32x4_t   decSpeed_n = vdupq_n_f32(decSpeed); 
  float32x4_t   incSpeed_n = vdupq_n_f32(incSpeed);
  float32x4x2_t tmpd0, tmpd1, tmpd2, tmpd3;
  float32x4_t   tmp0, tmp1, tmp2, tmp3;
  float32x4_t   tmp10, tmp11,tmp12,tmp13;
  uint32x4_t    cmpresult0, cmpresult1, cmpresult10, cmpresult11;

  tmpd0 = vld2q_f32(&totlev[0]);
  tmpd2 = vld2q_f32(&totlev[8]);

  for (m = 0; m < SUBUSED; m+=16)
  {
    tmpd1 = vld2q_f32(&noilev[m]);
    tmpd3 = vld2q_f32(&noilev[m+8]);

    // totlev[m] >= noilev[m] ?
    cmpresult0 = vcgeq_f32(tmpd0.val[0], tmpd1.val[0]);
    cmpresult1 = vcgeq_f32(tmpd0.val[1], tmpd1.val[1]);
    cmpresult10= vcgeq_f32(tmpd2.val[0], tmpd3.val[0]);
    cmpresult11= vcgeq_f32(tmpd2.val[1], tmpd3.val[1]);

    // totlev[m] - noilev[m]
    tmp0 = vsubq_f32(tmpd0.val[0], tmpd1.val[0]);
    tmp1 = vsubq_f32(tmpd0.val[1], tmpd1.val[1]);
    tmp10= vsubq_f32(tmpd2.val[0], tmpd3.val[0]);
    tmp11= vsubq_f32(tmpd2.val[1], tmpd3.val[1]);

    // mulor = totlev[m] >= noilev[m] ? incSpeed_n : decSpeed_n
    tmp2 = vbslq_f32(cmpresult0, incSpeed_n, decSpeed_n);
    tmp3 = vbslq_f32(cmpresult1, incSpeed_n, decSpeed_n);
    tmp12= vbslq_f32(cmpresult10,incSpeed_n, decSpeed_n);
    tmp13= vbslq_f32(cmpresult11,incSpeed_n, decSpeed_n);

    // product = (totlev[m] - noilev[m]) * mulor
    tmp2 = vmulq_f32(tmp0, tmp2);
    tmp3 = vmulq_f32(tmp1, tmp3);
    tmp12= vmulq_f32(tmp10, tmp12);
    tmp13= vmulq_f32(tmp11, tmp13);

    // preload to reduce the stall of vmul
    tmpd0 = vld2q_f32(&totlev[m+16]);
    tmpd2 = vld2q_f32(&totlev[m+24]);

    // noilev[m] += XXXSpeed * (totlev[m] - noilev[m]);
    tmp2 = vaddq_f32(tmpd1.val[0], tmp2);
    tmp3 = vaddq_f32(tmpd1.val[1], tmp3);
    tmp12= vaddq_f32(tmpd3.val[0], tmp12);
    tmp13= vaddq_f32(tmpd3.val[1], tmp13);
    
    tmpd1.val[0] = tmp2;
    tmpd1.val[1] = tmp3;
    tmpd3.val[0] = tmp12;
    tmpd3.val[1] = tmp13;

    vst2q_f32(&noilev[m],  tmpd1);
    vst2q_f32(&noilev[m+8],tmpd3);
  } 
  
  // Calculate to get the sum. Firstly, process 16x in the loop, then process 4x, 
  // and then process the remian. In fact, SUBUSED16 is 4x.
  float sum_noi;
  tmp3 = vdupq_n_f32(0.0f);
  for (m = 0; m < (SUBUSED16&~15); m+=16)
  {
    tmpd0  = vld2q_f32((float32_t *)&noilev[m]);
    tmpd2 = vld2q_f32((float32_t *)&noilev[m+8]);
    tmp3 = vaddq_f32(tmp3, tmpd0.val[0]);
    tmp3 = vaddq_f32(tmp3, tmpd0.val[1]);
    tmp3 = vaddq_f32(tmp3, tmpd2.val[0]);
    tmp3 = vaddq_f32(tmp3, tmpd2.val[1]);
  }
  for (; m < SUBUSED16; m+=4)
  {
    tmp1  = vld1q_f32((float32_t *)&noilev[m]);  
    tmp3 = vaddq_f32(tmp3, tmp1);
  }
  
  float32x2_t tmpX = vpadd_f32(vget_high_f32(tmp3), vget_low_f32(tmp3));
  sum_noi = vget_lane_f32(tmpX, 1) + vget_lane_f32(tmpX, 0);
#else
  for(int m = 0; m < SUBUSED; m++ )
  {
    if( totlev[m] < noilev[m] )
    {
      noilev[m] += decSpeed * (totlev[m] - noilev[m]);
    }
    else
    {
      noilev[m] += incSpeed * (totlev[m] - noilev[m]);
    }

  }  
  float sum_noi  = 0;
  for(int m = 0; m < SUBUSED16; m++)     
    sum_noi += noilev[m];
#endif

  return sum_noi;
}

/***************************************************************************
* SELECTGAIN
*    Auth.: IFA
*    Desc.: Choose which of the gains calculated to be used - for each subband.
***************************************************************************/
void noisereduction_selectgain(const float *noisegain,
                               const float *dereverbgain,
                               const float *nlpgain,
                               float * const  gain,
                               bool dereverbOn)
{
  if(dereverbOn) {
    for(int m = 0; m < SUBUSED; m++ )
    {
      gain[m] = min(dereverbgain[m], min(nlpgain[m], noisegain[m]));
    }
  }
  else{
    for(int m = 0; m < SUBUSED; m++ )
    {
      gain[m] = min(nlpgain[m], noisegain[m]);
    }
  }
}

/***************************************************************************
* CALCULATE NOISEGAIN
*   Auth.: IFA
*   Desc.: Optimised version uses fastdiv - less but enough precision. See
*          CCS-help and mathfun_inline.c. (_rcpsp and Newton's method)
***************************************************************************/
void noisereduction_calcnoisegain(float *  noisgn,
                                 const float *  noilev,
                                 const float *  totlev,
                                 float nrglimit)
{
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
  /*
  Xn+1 = Xn*(2-d*Xn) converge to 1/d 
  1/d ~= X0*(2-d*X0)
  */
  int m;
  float32x4_t   cst1SubNRGSP_n = vdupq_n_f32(1.0f - (NRGSP)); 
  float32x4_t   cstNRGSP_n = vdupq_n_f32(NRGSP); 
  float32x4_t   cst1dot0_n = vdupq_n_f32(1.0f);
  float32x4_t   nrglimit_n = vdupq_n_f32(nrglimit);
  float32x4_t   FLT_MIN_n = vdupq_n_f32(FLT_MIN);
  float32x4x2_t tmpd;
  float32x4_t   tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10;

  tmpd = vld2q_f32(&totlev[0]);
  tmp1 = vaddq_f32(tmpd.val[0], FLT_MIN_n);         // (totlev[m] + FLT_MIN)
  tmp2 = vaddq_f32(tmpd.val[1], FLT_MIN_n);
  tmp3 = vrecpeq_f32(tmp1);                         // X0
  tmp4 = vrecpeq_f32(tmp2);
    
  for (m = 0; m < SUBUSED; m+=8)
  {
    tmp5 = vrecpsq_f32(tmp1,tmp3);                  // (2-d*X0)
    tmp6 = vrecpsq_f32(tmp2,tmp4);

    tmpd = vld2q_f32(&noisgn[m]);
    tmp1 = vmulq_f32(tmpd.val[0], cst1SubNRGSP_n);  // (1.0f - (NRGSP)) * noisgn[m]
    tmp2 = vmulq_f32(tmpd.val[1], cst1SubNRGSP_n); 

    tmpd = vld2q_f32(&noilev[m]); 
    tmp3 = vmulq_f32(tmp3,tmpd.val[0]);             // noilev[m]*X0;
    tmp4 = vmulq_f32(tmp4,tmpd.val[1]);

    tmp9 = vaddq_f32(tmp1, cstNRGSP_n);             // (1.0f - (NRGSP)) * noisgn[m] + NRGSP
    tmp10= vaddq_f32(tmp2, cstNRGSP_n);    
    
    tmp7 = vmulq_f32(tmp3, tmp5);                   // noilev[m]*X0 *(2-d*X0) = noilev[m]/(totlev[m] + FLT_MIN);
    tmp8 = vmulq_f32(tmp4, tmp6);
 
    tmpd = vld2q_f32(&totlev[m+8]);
    tmp1 = vaddq_f32(tmpd.val[0], FLT_MIN_n);       // NEXT (totlev[m] + FLT_MIN)
    tmp2 = vaddq_f32(tmpd.val[1], FLT_MIN_n);

    tmp7 = vsubq_f32(tmp9, tmp7);                   // (1.0f - (NRGSP)) * noisgn[m] + NRGSP + noilev[m]/(totlev[m] + FLT_MIN);
    tmp8 = vsubq_f32(tmp10,tmp8);
    
    tmp7 = vminq_f32(tmp7,cst1dot0_n);
    tmp8 = vminq_f32(tmp8,cst1dot0_n);
    
    tmp7 = vmaxq_f32(tmp7,nrglimit_n);
    tmp8 = vmaxq_f32(tmp8,nrglimit_n); 

    tmp3 = vrecpeq_f32(tmp1);                       // NEXT X0
    tmp4 = vrecpeq_f32(tmp2);

    tmpd.val[0] = tmp7;
    tmpd.val[1] = tmp8;
    vst2q_f32(&noisgn[m],tmpd);
  }
    
#else
  for(int m = 0; m < SUBUSED; m++ )
  {
    float ng = (1.0f - (NRGSP)) * noisgn[m] + (NRGSP) - noilev[m]/(totlev[m] + FLT_MIN);
    if( ng > 1.0f )
    {
      ng = 1.0f;  /* "do not amplify noise" */
    }
    noisgn[m] = max(ng, nrglimit);   /* "do not try to attenuate noise more than nrglimit allows for" */
  }
#endif
}

/***************************************************************************
* CALCULATE COMFORT NOISE
*    Auth.: IFA
*    Desc.: Here the calculations are based on the effect (i.e.the square)
*           of the different signals.
*           Comfort noise is added mainly when there is little near end sound.
*           Thus the whole near end signal is damped, and the amount of comfort
*           noise added (random noise in each subband) is consequently larger. [...]
***************************************************************************/
int noisereduction_noise(
  float rmax_inv,		
  int tmp_seed,
  float nrgl2,
  float comf_noise_ampl,
  const float * gain,
  float * no,
  COMPLEX32 * noise,
  float * nfgn)
{
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
  float32x4_t   tmpd0, tmpd1;//, tmpd2, tmpd3;
  float32x4_t   tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  float32x4x2_t seed0, seed1;
  float32x4_t   nrgl2_n   = vdupq_n_f32(nrgl2);
  float32x4_t   FLT_MIN_n = vdupq_n_f32(FLT_MIN);
  float32x4_t   cna_mul_ri = vdupq_n_f32(comf_noise_ampl*rmax_inv);
  float32x4_t   cna_mul_05 = vdupq_n_f32(comf_noise_ampl*0.5f);
  
  for(int m = 0; m < SUBUSED; m+=8 )
  {
    COMPLEX32 seed[8];
    for (int n = 0; n < 8; n++)
    {
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      seed[n].re = (float) tmp_seed;
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      seed[n].im = (float) tmp_seed;
    }
    
    /*
    for (int n = 0; n < 8; n++)
    {
      float tmpdiff = nrgl2 - gain[m+n] * gain[m+n];
      tmpdiff = max(tmpdiff, FLT_MIN);
      nfgn_tmp = sqrtf(tmpdiff);
      float  comf_noise_ampl_5 = comf_noise_ampl*0.5f;
      float  comf_noise_ampl_rmax_inv = comf_noise_ampl*rmax_inv;
      float  nfgn_tmp_no = nfgn_tmp * no[m+n];
      noise[m+n].re = (comf_noise_ampl_rmax_inv * seed[n].re - comf_noise_ampl_5) * nfgn_tmp_no;
      noise[m+n].im = (comf_noise_ampl_rmax_inv * seed[n].im - comf_noise_ampl_5) * nfgn_tmp_no;
      nfgn[m+n] = nfgn_tmp;
    }
    */

    /*
    Xn+1 = Xn*(3-d*Xn*Xn)/2 converge to 1/sqrt(d) 
    sqrt(d) ~= d * X0*(3-d*X0*X0)/2
    */
    
    tmpd0  = vld1q_f32(&gain[m]);
    tmpd1  = vld1q_f32(&gain[m+4]);
    tmp0  = vmulq_f32(tmpd0, tmpd0);
    tmp1  = vmulq_f32(tmpd1, tmpd1);

    seed0 = vld2q_f32((float*)&seed[0]);
    seed1 = vld2q_f32((float*)&seed[4]); 

    tmp0  = vsubq_f32(nrgl2_n, tmp0);
    tmp1  = vsubq_f32(nrgl2_n, tmp1);
    
    tmp0  = vmaxq_f32(FLT_MIN_n, tmp0);     // d
    tmp1  = vmaxq_f32(FLT_MIN_n, tmp1);

    tmp2  = vrsqrteq_f32(tmp0);             // X0
    tmp3  = vrsqrteq_f32(tmp1);

    tmp6  = vmulq_f32(cna_mul_ri, seed0.val[0]);  // seed[n].re * comf_noise_ampl * rmax_inv
    tmp7  = vmulq_f32(cna_mul_ri, seed1.val[0]);
    tmp8  = vmulq_f32(cna_mul_ri, seed0.val[1]);  // seed[n].im * comf_noise_ampl * rmax_inv
    tmp9  = vmulq_f32(cna_mul_ri, seed1.val[1]);    

    tmp4  = vmulq_f32(tmp2, tmp2);        // X0 *X0
    tmp5  = vmulq_f32(tmp3, tmp3);

    tmp4  = vrsqrtsq_f32(tmp0, tmp4);     // (3-d*X0*X0)/2 
    tmp5  = vrsqrtsq_f32(tmp1, tmp5);   

    tmp2  = vmulq_f32(tmp0, tmp2);        // d * X0
    tmp3  = vmulq_f32(tmp1, tmp3);
 
    tmp6  = vsubq_f32(tmp6, cna_mul_05);
    tmp7  = vsubq_f32(tmp7, cna_mul_05);
    
    tmp4  = vmulq_f32(tmp2, tmp4);        // d * X0 * (3-d*X0*X0)/2 => 1/sqrt(d)
    tmp5  = vmulq_f32(tmp3, tmp5);

    tmpd0 = vld1q_f32(&no[m]);
    tmpd1 = vld1q_f32(&no[m+4]);
    
    tmp4  = vmulq_f32(tmpd0, tmp4);
    tmp5  = vmulq_f32(tmpd1, tmp5);       // nfgn_tmp * no[m+n];

    tmp8  = vsubq_f32(tmp8, cna_mul_05);
    tmp9  = vsubq_f32(tmp9, cna_mul_05);

    tmp6  = vmulq_f32(tmp6, tmp4);        // (comf_noise_ampl_rmax_inv * seed[n].re - comf_noise_ampl_5) * nfgn_tmp_no;
    tmp7  = vmulq_f32(tmp7, tmp5);
    tmp8  = vmulq_f32(tmp8, tmp4);        // (comf_noise_ampl_rmax_inv * seed[n].im - comf_noise_ampl_5) * nfgn_tmp_no;
    tmp9  = vmulq_f32(tmp9, tmp5); 

    vst1q_f32(&nfgn[m],   tmp4);          // nfgn[m+n] = nfgn_tmp;
    vst1q_f32(&nfgn[m+4], tmp5);

    seed0.val[0] = tmp6;
    seed1.val[0] = tmp7;
    seed0.val[1] = tmp8;
    seed1.val[1] = tmp9;

    vst2q_f32((float*)&noise[m],   seed0);
    vst2q_f32((float*)&noise[m+4], seed1);
  }
#else
  float nfgn_tmp, rand;
  for(int m = 0; m < SUBUSED; m+=8 )
  {
    COMPLEX32 seed[8];
    for (int n = 0; n < 8; n++)
    {
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      seed[n].re = (float) tmp_seed;
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      tmp_seed = tmp_seed >> 1 ^ (tmp_seed & 0x1 ? (MLS_POLYNOME) : 0);
      seed[n].im = (float) tmp_seed;
    }
    for (int n = 0; n < 8; n++)
    {
      float tmpdiff = nrgl2 - gain[m+n] * gain[m+n];
      tmpdiff = max(tmpdiff, FLT_MIN);
      nfgn_tmp = sqrtf(tmpdiff);
      rand = comf_noise_ampl * (seed[n].re * rmax_inv - 0.5f);
      noise[m+n].re = nfgn_tmp * no[m+n] * rand;
      rand = comf_noise_ampl * (seed[n].im * rmax_inv - 0.5f);
      noise[m+n].im = nfgn_tmp * no[m+n] * rand;
      nfgn[m+n] = nfgn_tmp;
    }
  }	
#endif
  return(tmp_seed);
}

/***************************************************************************
* CALCULATE COMFORT NOISE
*    Auth.: IFA
*    Desc.: Here the calculations are based on the effect (i.e.the square)
*           of the different signals.
*           Comfort noise is added mainly when there is little near end sound.
*           Thus the whole near end signal is damped, and the amount of comfort
*           noise added (random noise in each subband) is consequently larger. [...]
***************************************************************************/
void noisereduction_calccomfortnoise(NOIRED_PTR pNoired,
                                            const float * gain,
                                            COMPLEX32 * noise,
                                            float * nfgn)
{
  //  float nfgn[SUBUSED48]; correct
  float        nrgl2    = pNoired->nrglimit * pNoired->nrglimit; /* to avoid repeating this in the loop */
  unsigned int tmp_seed = pNoired->seed;
  float        comf_noise_ampl;
  float        y0, y1, a, b;

  /*
  * Calculate comf_noise_amp_factor.
  * (pNoired->comNoiseAmp less than 1.12 means attenuating cnoise at high noiselevels)
  * If noise is high, nrgl=-6dB, comfnoise is not further amplified
  * i.e. nrgl = nrglmin => comf_noise_ampl = 1.12
  * and  nrgl = nrglmax => comf_noise_ampl = 1.0
  */
  y1 = 1.00f * pNoired->comNoiseAmp;
  y0 = 1.12f;

  a = fasterdiv((y1 - y0), (NRGLMAX - NRGLMIN));
  b = y0 - (NRGLMIN) * a;

  comf_noise_ampl  = a * pNoired->nrglimit + b;
  comf_noise_ampl *= 2.0f;      /* compensating for the random-numbers being in [0,1] */
  if( pNoired->ComNoiseOn == false ) comf_noise_ampl = 0.0f; /* turning comfort noise off for test */

  tmp_seed = noisereduction_noise(	
    MLS_RAND_MAX_INV,		  
    tmp_seed,		
    nrgl2,
    comf_noise_ampl,
    gain,
    pNoired->noilev,
    noise,
    nfgn);

  pNoired->seed = tmp_seed;
}

/***************************************************************************
* CALCULATE NLP GAIN
*   Auth.: IFA
*   Desc.:
***************************************************************************/
void noisereduction_calcnlpgain(NOIRED_PTR pNoired,
                                float excgain,
                                float *  nlpgain,
                                float nlpconst)
{
  float fullg = 0.0f, fullgCorr = 0.0f, fullgNlp = 0.0f;
  float fullgHi = 0.0f, fullgCorrHi = 0.0f, fullgNlpHi = 0.0f;
  float extrag = 0.0f;
  float aft = 0.0f, bef = 0.0f, est = 0.0f;
  float tmp;
  float fullgn = 0.0f, fullgd = 0.0f;
  float fullgnCorr = 0.0f, fullgdCorr = 0.0f;
  short slowflag = 0, slowflagHi = 0;
  int m = 0;
  float levpow = 1.25f;
  float powerFrac = 1.60f;
//  float NLP_SLOW_SPD_ = 0xDE.01p0;

  /* Use nlpgain pointer directly in calculation of subg. It saves use of scratchmem */
  //float *  subg;
  //subg = scratchmem_alloc(sizeof(float) * SUBUSED);

  for( m = NLP_SUBUSED_START; m < NLP_SUBUSED_END; m++ )
  {
    aft += pNoired->aftlev[m];
    bef += pNoired->beflev[m];
    est += pNoired->estlev[m];
  }

  /* fullg = 2*((aft / bef)^(2)  */
  fullgn = aft;
  fullgd = bef;
  fullgn = ( (fullgn > 0.0f) ? (fullgn) : (0.0f) ); /* max(0,fullgn) */
  fullgd = ( (fullgd > 0.0f) ? (fullgd) : (0.0f) ); /* max(0,fullgd) */
  fullg = fasterdiv(fullgn, (fullgd + FLT_MIN));
  fullg = nlpconst * fullg * fullg;
  fullg = ( (fullg <= 1.0f) ? (fullg) : (1.0f) ); /* min(1,fullg) */

  /* updating fullgain slowly if beflev is falling */
  tmp = 0.5f * (pNoired->beflevPrev[0] + pNoired->beflevPrev[1]) ;
  if( tmp >= bef ) /* if beflev has a falling slope */
  {
    fullg = fullg * NLP_SLOW_SPD + (1.0f - NLP_SLOW_SPD) * pNoired->fullgPrev;
    slowflag = 1;
  }
  pNoired->fullgPrev = fullg;
  pNoired->beflevPrev[0] = pNoired->beflevPrev[1];
  pNoired->beflevPrev[1] = bef;

  /* Correction factor fullgCorr = 2*((bef^c-est^c) / bef^c)^(b/c)  */
  aft = powf(aft, levpow);
  bef = powf(bef, levpow);
  est = powf(est, levpow);

  fullgnCorr = bef - est;
  fullgdCorr = bef;
  fullgnCorr = ( (fullgnCorr > 0.0f) ? (fullgnCorr) : (0.0f) ); /* max(0,fullgn) */
  fullgdCorr = ( (fullgdCorr > 0.0f) ? (fullgdCorr) : (0.0f) ); /* max(0,fullgd) */

  fullgCorr = fasterdiv(fullgnCorr, (fullgdCorr + FLT_MIN));
  fullgCorr = nlpconst * powf(fullgCorr, powerFrac);
  fullgCorr = ( (fullgCorr <= 1.0f) ? (fullgCorr) : (1.0f) ) ; /* min(1,fullgCorr) */

  if( fullgCorr<fullg && pNoired->fullgNlpOn )
  {
    fullgNlp = fasterdiv(fullgCorr, fullg + FLT_MIN);
  }
  else
  {
    fullgNlp = 1.0f;
  }

  /* calculate single talk additional attenuation */
  if( pNoired->extragainOn )
  {
    extrag = ( (fullgCorr < SGTLIM) ? (fullgCorr/SGTLIM) : (1.0f) ); // min(1,fullg/SGTLIM) - extra atteunation if guaranteed single talk
  }
  else
  {
    extrag = 1.0f;
  }

  if( pNoired->NlpHiOn )
  {
    /****** start calculations for higher subbands ************/
    aft = 0.0f;
    bef = 0.0f;
    est = 0.0f;
    for( m = NLP_SUBUSED_END; m < pNoired->nlpSubusedEndHi; m++ )
    {
      aft += pNoired->aftlev[m];
      bef += pNoired->beflev[m];
      est += pNoired->estlev[m];
    }

    /* fullg = 2*((aft / bef)^(2)  */
    fullgn = aft;
    fullgd = bef;

    fullgn = ( (fullgn > 0.0f) ? (fullgn) : (0.0f) ); /* max(0,fullgn) */
    fullgd = ( (fullgd > 0.0f) ? (fullgd) : (0.0f) ); /* max(0,fullgd) */
    fullgHi = fasterdiv(fullgn, (fullgd + FLT_MIN));
    fullgHi = nlpconst * fullgHi* fullgHi;
    fullgHi = ( (fullgHi <= 1.0f) ? (fullgHi) : (1.0f) ); /* min(1,fullg) */

    /* updating fullgain slowly if beflev is falling */
    tmp = 0.5f * (pNoired->beflevPrevHi[0] + pNoired->beflevPrevHi[1]) ;
    if( tmp >= bef ) /* if beflev has a falling slope */
    {
      fullgHi = fullgHi * NLP_SLOW_SPD + (1.0f - NLP_SLOW_SPD) * pNoired->fullgPrevHi;
      slowflagHi = 1;
    }
    pNoired->fullgPrevHi = fullgHi;
    pNoired->beflevPrevHi[0] = pNoired->beflevPrevHi[1];
    pNoired->beflevPrevHi[1] = bef;

    /* Correction factor fullgCorr = 2*((bef^c-est^c) / bef^c)^(b/c)  */
    aft = powf(aft, levpow);
    bef = powf(bef, levpow);
    est = powf(est, levpow);

    fullgnCorr = bef - est;
    fullgdCorr = bef;
    fullgnCorr = ( (fullgnCorr > 0.0f) ? (fullgnCorr) : (0.0f) ); /* max(0,fullgn) */
    fullgdCorr = ( (fullgdCorr > 0.0f) ? (fullgdCorr) : (0.0f) ); /* max(0,fullgd) */

    fullgCorrHi = fasterdiv(fullgnCorr, (fullgdCorr + FLT_MIN));
    fullgCorrHi = nlpconst * powf(fullgCorrHi, powerFrac);
    fullgCorrHi = ( (fullgCorrHi <= 1.0f) ? (fullgCorrHi) : (1.0f) ) ; /* min(1,fullgCorr) */

    if( fullgCorrHi<fullgHi && pNoired->fullgNlpHiOn )
    {
      fullgNlpHi = fasterdiv(fullgCorrHi, fullgHi + FLT_MIN);
    }
    else
    {
      fullgNlpHi = 1.0f;
    }
  } /********** end calculations for higher subbands **************/
  else
  {
    fullgHi     = fullg;
    fullgCorrHi = fullgCorr;
    fullgNlpHi  = fullgNlp;
    slowflagHi  = slowflag;
  }

  /* calculate subband gain (vector) */
  /* 2*((aft / bef)^(2) */
  {
    float * const  subgn_tmp = pNoired->aftlev;
    float * const  subgd_tmp = pNoired->beflev;
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
    float32x4x2_t tmpd0, tmpd1, tmpd2,tmpd10, tmpd11, tmpd12;
    float32x4_t   tmp0,tmp1,tmp2,tmp3,tmp4,tmp5;
    float32x4_t   tmp10,tmp11,tmp12,tmp13,tmp14,tmp15;
    float32x4_t nlpconst_n = vdupq_n_f32(nlpconst);
    float32x4_t cst1dot0_n = vdupq_n_f32(1.0f);

    /*
    Xn+1 = Xn*(3-d*Xn*Xn)/2 converge to 1/sqrt(d) 
    sqrt(d) ~= d * X0*(3-d*X0*X0)/2
    */
    
    for( m = 0; m < SUBUSED; m+=16 )
    {
      tmpd0 = vld2q_f32(&subgd_tmp[m]);
      tmp0 = tmpd0.val[0];          // d
      tmp1 = tmpd0.val[1]; 
      tmp2 = vrecpeq_f32(tmp0);    // X0
      tmp3 = vrecpeq_f32(tmp1);
      
      tmpd10 = vld2q_f32(&subgd_tmp[m+8]);
      tmp10 = tmpd10.val[0];          // d
      tmp11 = tmpd10.val[1]; 
      tmp12 = vrecpeq_f32(tmp10);    // X0
      tmp13 = vrecpeq_f32(tmp11);

      tmp0 = vrecpsq_f32(tmp0,tmp2);   // (2-d*X0)  
      tmp1 = vrecpsq_f32(tmp1,tmp3);
      tmp10 = vrecpsq_f32(tmp10,tmp12);   // (2-d*X0)  
      tmp11 = vrecpsq_f32(tmp11,tmp13);
      
      tmpd1 = vld2q_f32(&subgn_tmp[m]);
      tmp4 = tmpd1.val[0];          // d
      tmp5 = tmpd1.val[1]; 
      tmpd11 = vld2q_f32(&subgn_tmp[m+8]);
      tmp14 = tmpd11.val[0];          // d
      tmp15 = tmpd11.val[1]; 

      tmp4 = vmulq_f32(tmp4,tmp2);      // divided_num*X0
      tmp5 = vmulq_f32(tmp5,tmp3);   
      tmp14 = vmulq_f32(tmp14,tmp12);      // divided_num*X0
      tmp15 = vmulq_f32(tmp15,tmp13);   

      tmp4 = vmulq_f32(tmp4,tmp0);      // divided_num*X0*(2-d*X0) => gain
      tmp5 = vmulq_f32(tmp5,tmp1); 
      tmp14 = vmulq_f32(tmp14,tmp10);      // divided_num*X0*(2-d*X0) => gain
      tmp15 = vmulq_f32(tmp15,tmp11); 


      tmp4 = vmulq_f32(tmp4,tmp4);      // gain * gain
      tmp5 = vmulq_f32(tmp5,tmp5); 
      tmp14 = vmulq_f32(tmp14,tmp14);      // gain * gain
      tmp15 = vmulq_f32(tmp15,tmp15); 

      tmp4 = vmulq_f32(nlpconst_n,tmp4);      // nlpconst * gain * gain
      tmp5 = vmulq_f32(nlpconst_n,tmp5); 
      
      tmp14 = vmulq_f32(nlpconst_n,tmp14);      // nlpconst * gain * gain
      tmp15 = vmulq_f32(nlpconst_n,tmp15);

      tmp4 = vminq_f32(cst1dot0_n,tmp4);      // min(1.0* nlpconst * gain * gain)
      tmp5 = vminq_f32(cst1dot0_n,tmp5); 

      tmp14 = vminq_f32(cst1dot0_n,tmp14);      // min(1.0* nlpconst * gain * gain)
      tmp15 = vminq_f32(cst1dot0_n,tmp15); 

      tmpd2.val[0] = tmp4;
      tmpd2.val[1] = tmp5;

      vst2q_f32(&nlpgain[m], tmpd2);
      
      tmpd12.val[0] = tmp14;
      tmpd12.val[1] = tmp15;
      vst2q_f32(&nlpgain[m+8], tmpd12);
      
    }
#else
    for( m = 0; m < SUBUSED; m++ )
    {
      nlpgain[m] = fasterdiv(subgn_tmp[m], subgd_tmp[m] + FLT_MIN);
      nlpgain[m] = nlpconst * nlpgain[m] * nlpgain[m];
      nlpgain[m] = min(1.0f, nlpgain[m]);
      //pNoired->subgPrev[m] = nlpgain[m];
    }
#endif
  }

  /* updating subgain slowly if beflev is falling */
  if (slowflag)
  {
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
    float32x4x2_t tmpd0, tmpd1, tmpd2;
    float32x4_t   tmp0,tmp1,tmp2,tmp3,tmp4,tmp5;
    float32x4_t   NLP_SLOW_SPD_n, NLP_SLOW_SPD_rsb1_n;

    NLP_SLOW_SPD_n      = vdupq_n_f32(NLP_SLOW_SPD);
    NLP_SLOW_SPD_rsb1_n = vdupq_n_f32(1.0f - NLP_SLOW_SPD);

    // To reduce stall of vmul, put the part of first iteration out of loop.
    tmpd0 = vld2q_f32(&nlpgain[0]);
    tmpd1 = vld2q_f32(&pNoired->subgPrev[0]);
    tmp0 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[0]);
    tmp1 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[0]);
    tmp2 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[1]);
    tmp3 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[1]); 
    // firstly do 8x.
    for( m=0; m < (NLP_SUBUSED_END&~7); m+=8 )
    {
      tmp4 = vaddq_f32(tmp0, tmp1);
      tmp5 = vaddq_f32(tmp2, tmp3);
      tmpd2.val[0] = tmp4;
      tmpd2.val[1] = tmp5;
      vst2q_f32(&nlpgain[m],tmpd2);
      tmpd0 = vld2q_f32(&nlpgain[m+8]);
      tmpd1 = vld2q_f32(&pNoired->subgPrev[m+8]);
      tmp0 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[0]);
      tmp1 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[0]);
      tmp2 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[1]);
      tmp3 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[1]);
    } 
    // then process the remain which is less than 8.
    for( ; m < NLP_SUBUSED_END; m++ )
    {
      nlpgain[m] = nlpgain[m] * NLP_SLOW_SPD + (1.0f - NLP_SLOW_SPD) * pNoired->subgPrev[m];
    }
#else
    for( m = 0; m < NLP_SUBUSED_END; m++ )
    {
      nlpgain[m] = nlpgain[m] * NLP_SLOW_SPD + (1.0f - NLP_SLOW_SPD) * pNoired->subgPrev[m];
    }
#endif
  }
  if (slowflagHi)
  {
#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
    float32x4x2_t tmpd0, tmpd1, tmpd2;
    float32x4_t   tmp0,tmp1,tmp2,tmp3,tmp4,tmp5;
    float32x4_t   NLP_SLOW_SPD_n, NLP_SLOW_SPD_rsb1_n;

    // firstly do it in order to make nlpgain aligned 32-byte.
    for( m = NLP_SUBUSED_END; m < ((NLP_SUBUSED_END+7) & ~7); m++ )
    {
      nlpgain[m] = nlpgain[m] * NLP_SLOW_SPD + (1.0f - NLP_SLOW_SPD) * pNoired->subgPrev[m];
    }

    // To reduce stall of vmul, put the part of first iteration out of main loop.
    NLP_SLOW_SPD_n      = vdupq_n_f32(NLP_SLOW_SPD);
    NLP_SLOW_SPD_rsb1_n = vdupq_n_f32(1.0f - NLP_SLOW_SPD);
    tmpd0 = vld2q_f32(&nlpgain[m]);
    tmpd1 = vld2q_f32(&pNoired->subgPrev[m]);
    tmp0 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[0]);
    tmp1 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[0]);
    tmp2 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[1]);
    tmp3 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[1]);    
    for( ; m < SUBUSED; m+=8 )
    {
      tmp4 = vaddq_f32(tmp0, tmp1);
      tmp5 = vaddq_f32(tmp2, tmp3);
      tmpd2.val[0] = tmp4;
      tmpd2.val[1] = tmp5;
      vst2q_f32(&nlpgain[m],tmpd2);
      tmpd0 = vld2q_f32(&nlpgain[m+8]);
      tmpd1 = vld2q_f32(&pNoired->subgPrev[m+8]);
      tmp0 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[0]);
      tmp1 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[0]);
      tmp2 = vmulq_f32(NLP_SLOW_SPD_n, tmpd0.val[1]);
      tmp3 = vmulq_f32(NLP_SLOW_SPD_rsb1_n, tmpd1.val[1]);
    } 
#else
    for( m = NLP_SUBUSED_END; m < SUBUSED; m++ )
    {
      nlpgain[m] = nlpgain[m] * NLP_SLOW_SPD + (1.0f - NLP_SLOW_SPD) * pNoired->subgPrev[m];
    }
#endif
  }
  
  for( m = 0; m < SUBUSED; m++ )
  {
  	pNoired->subgPrev[m] = nlpgain[m];
  } 
  
 
  /* TWISTING of subgain */
  /* Limit nlpgain per subband to avoid 'Mickey-mouse' effect
  * subbands below NLP_SUBUSED_END are limited to +/- 12 dB according to
  *   the calculated NLP_fullgain(='average' value for the subbands NLP_SUBUSED_START to NLP_SUBUSED_END)
  * subbands from and above NLP_SUBUSED_END are limited to +/- 12dB according to
  *   the calculated NLP_fullgainHi(='average' value for the subbands NLP_SUBUSED_END to pNoired->nlpSubusedEndHi)
  */
  {
    float maxg = 0.0f, ming = 0.0f;

    maxg =( (1.0f < fullg*MAXTWIST ) ? (1.0f) : (fullg*MAXTWIST) ); // min(1,params.maxtwist*fullg);
    ming = fasterdiv(fullg, MAXTWIST);

    for( m = 0; m < NLP_SUBUSED_END; m++ )      /* assumes that maxg > ming */
    {
      if( nlpgain[m] > maxg )                    /* nlpgain = min(maxg,subg) */
      {
        nlpgain[m] = maxg;
      }
      else if( nlpgain[m] < ming )
      {
        nlpgain[m] = ming;
      }
    }

    /* updating the twist-limits for higher subbands */
    maxg =( (1.0f < fullgHi*MAXTWIST ) ? (1.0f) : (fullgHi*MAXTWIST) ); // min(1,params.maxtwist*fullg);
    ming = fasterdiv(fullgHi, MAXTWIST);

    for( m = NLP_SUBUSED_END; m < SUBUSED; m++ )      /* twisting the rest of the subbands */
    {
      if( nlpgain[m] > maxg )
      {
        nlpgain[m] = maxg;
      }
      else if( nlpgain[m] < ming )
      {
        nlpgain[m] = ming;
      }
    }
  }

  /* calculating the final nlpgain */
  /* adding extra gain (single talk gain) and exception handler gain */
  if( pNoired->NlpGainOn )
  {
    if( pNoired->NLPSubOn )     /* listening to subgain */
    {
      for( m = 0; m < NLP_SUBUSED_END; m++ )
      {
        nlpgain[m] = nlpgain[m] * fullgNlp * extrag * excgain;
      }


      { /* 10 subbands used as transition area of fullgNlp between low and high frequencies */
        float diff = (fullgNlpHi - fullgNlp) * 0.1f;
        float fullgNlpTmp = fullgNlp;
        int loopend = NLP_SUBUSED_END + pNoired->fullgNlpTransition;
        for( m = NLP_SUBUSED_END; m < loopend; m++ )
        {
          fullgNlpTmp += diff;
          nlpgain[m] = nlpgain[m] * fullgNlpTmp * extrag * excgain;
        }
      }

#if defined(__ARM_NEON__) && defined(ENV_IOS)    // should work on Android, but leaving for now
      {
        float32x4x2_t tmpd0, tmpd1;
        float32x4_t   tmp0,tmp1,tmp2,tmp3;
        float32x4_t   fullgee_n;
        float fullgee = fullgNlpHi * extrag * excgain;
        int n;

        // firstly do it in order to prepare 4x process and make nlpgain aligned 16-byte.
        m = NLP_SUBUSED_END + pNoired->fullgNlpTransition;
        n = (m+3)&~3;
        for( ; m < n; m++ )                     /* twisting the rest of the subbands */
        {
          nlpgain[m] = nlpgain[m] * fullgee;
        }

        fullgee_n = vdupq_n_f32(fullgee);

        // do it 4x in order to prepare 16x process and make nlpgain aligned 64-byte.
        tmp1 = vld1q_f32(&nlpgain[m]);
        tmp0 = vmulq_f32(fullgee_n, tmp1);
        n = (m+15)&~15;
        for( ; m < n; m+=4 )                    /* twisting the rest of the subbands */
        {
           tmp1 = vld1q_f32(&nlpgain[m+4]);
           vst1q_f32(&nlpgain[m],tmp0);
           tmp0 = vmulq_f32(fullgee_n, tmp1);        
        }

        // do it 16x.
        tmpd0 = vld2q_f32(&nlpgain[m]);
        tmp0 = vmulq_f32(fullgee_n, tmpd0.val[0]);
        tmp1 = vmulq_f32(fullgee_n, tmpd0.val[1]);
        for( ; m < SUBUSED; m+=16 )             /* twisting the rest of the subbands */
        {
          tmpd0 = vld2q_f32(&nlpgain[m+8]);
          tmp2 = vmulq_f32(fullgee_n, tmpd0.val[0]);
          tmp3 = vmulq_f32(fullgee_n, tmpd0.val[1]);
          
          tmpd1.val[0] = tmp0;
          tmpd1.val[1] = tmp1;
          vst2q_f32(&nlpgain[m],tmpd1);
          
          tmpd0 = vld2q_f32(&nlpgain[m+16]);
          tmp0 = vmulq_f32(fullgee_n, tmpd0.val[0]);
          tmp1 = vmulq_f32(fullgee_n, tmpd0.val[1]);
          
          tmpd1.val[0] = tmp2;
          tmpd1.val[1] = tmp3;
          vst2q_f32(&nlpgain[m+8],tmpd1);
        }        
      }
#else  
      for( m = NLP_SUBUSED_END + pNoired->fullgNlpTransition; m < SUBUSED; m++ )
      {
        nlpgain[m] = nlpgain[m] * fullgNlpHi * extrag * excgain;
      }
#endif
    }
    else                      /* listening to fullgain */
    {
      for( m = 0; m < SUBUSED; m++ )
      {
        nlpgain[m] = fullg * extrag * excgain;
      }
    }
  }
  else                        /* turns off nlp-gain, for test purposes */
  {
    for( m = 0; m < SUBUSED; m++ )
    {
      nlpgain[m] = 1.0f;
    }
  }

  pNoired->shellnlp = fullgNlp * extrag * excgain;  
}

/***************************************************************************
* CHECKNOISE
*   Auth.: IFA
*   Desc.: Returned value indicates whether there is noise or not.
*   Modifications:
*      2003-02-07 TFM: Noise defined also if sum noilev > sum beflev
*      2003-11-05 BWI/IFA: If totlev idicates more than 10dB noise-increase
*                          the noise must prove to be stable for a longer
*                          time -see comment in code.
*      2004-?-?   IFA: if stereo, noisClim is increased - to help avoid
*                      detecting music as noise.
***************************************************************************/
short noisereduction_checknoise(NOIRED_PTR pNoired)
{
    short isnoise;
    float beflev;
    float noilev;
    float totlev;
    float tmp;
    int newclim;

    beflev = pNoired->sum_beflev;
    noilev = pNoired->sum_noilev;
    totlev = pNoired->sum_totlev;

    if( beflev > FLT_MIN )
    {
        if( pNoired->ncnt_boot > 0 ) /* just booted or plugged in mic */
        {
            pNoired->ncnt_boot--;
        }
    }
    else /* no signal */
    {
        pNoired->ncnt_boot = 2 * NOISECLIM;
    }

    if( fabsf(beflev - pNoired->noiseamp) < beflev * NOISEDLIM ) /* stable level */
    {
        /* update counter */
        if( ++pNoired->ncnt > pNoired->ncnt_clim )
        {
            pNoired->ncnt = pNoired->ncnt_clim;
            /* This (if) means: ncnt=min(ncnt+1,NOISCLIM) */
        }
        /* update counter limit: */
        tmp = fasterdiv(totlev, noilev + FLT_MIN);
        newclim = (int)((20 * fastlog(tmp)) * 3.0f); /* increase 3 iterations per dB over 10 */

        /* limiting downwards */
        if( newclim < NOISECLIM )
        {
            newclim = NOISECLIM;
        }
        /* limiting upwards MONO */
        if( newclim > (NOISECLIM + 12) )
        {
            newclim = NOISECLIM + 12;
        }
        pNoired->ncnt_clim = (short)newclim;
    }
    else
    {
        pNoired->noiseamp = beflev;
        pNoired->ncnt = 0;
        pNoired->ncnt_clim = NOISECLIM;
    }
    /* standard ISNOISE kriterion: */
    if( pNoired->ncnt >= pNoired->ncnt_clim )
    {
        isnoise = 1;
    }
    /* other ISNOISE kriteria: */
    else if( noilev > beflev )
    {
        isnoise = 1;
    }
    else if( pNoired->ncnt_boot > 0 )
    {
        isnoise = 1;
    }
    /* default - slow noilev speed: */
    else
    {
        isnoise = 0;
    }


    return isnoise;
}

/***************************************************************************
* CALCULATE NOISEGAIN LIMIT
*   Auth.: IFA
*   Desc.: Calculates the nrglimit as sum of all noiselevels in the subbands used
*          in noisereduction; then assures that it is neither too large nor
*          too small.  nrglimit is thus the same for every subband.
***************************************************************************/
void noisereduction_calcnoisegainlimit(NOIRED_PTR pNoired, bool noisereductionOn)
{
    float tmpgain = 0.0f;

    tmpgain = 16.0f *  pNoired->sum_noilev + 0.125f;      /* (1/8);  more noise => less relative noise reduction */
    /* limits the gain to a band between the MIN- and MAX-limit, max(NRGLMIN,min(NRGLMAX,tmpgain)) */
    tmpgain = min(tmpgain, NRGLMAX);
    tmpgain = max(tmpgain, NRGLMIN);

    pNoired->nrgoptimal = pNoired->nrgoptimal + (tmpgain - pNoired->nrgoptimal) * NRGL_UPDSPEED; /* updates slowly */
    if( (noisereductionOn) && (pNoired->NoiseRedOn) ) /* two ways to turn off noisered.; com audioctrl (and menu) and directly via test-variable pNoired->NoiseRedOn */
    {
        pNoired->nrglimit = pNoired->nrgoptimal;
    }
    else
    {
        pNoired->nrglimit = 1.0f; /* setting noiseregulationlimit to 1 turns off noisereduction without affecting nlpgain */
    }
}
