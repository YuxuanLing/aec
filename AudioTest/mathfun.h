/*
 * Header for mathfun.
 *
 * This header should contain the interface for all functions 
 * exported by mathfun, including inlined functions (commented
 * out, as below), as well as includes to platform-specific 
 * implementations.
 *
 */
#ifndef MATHFUN_H
#define MATHFUN_H

#include <math.h>
#include <stdint.h>
#include <limits.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define logof2 0.69314718055994530941723212145818F
#define logof2inv 1.4426950408889634073599246810019F
#define dB2LinBase 1.2589254117941672104239541063958F

#define PI16 0.19634954084936207740391521145497F
#define PI16inv 5.0929581789406507446042804279205F

#ifndef PI
#define PI  3.14159265358979F
#endif
#ifndef dPI
#define dPI  3.14159265358979
#endif
#define LN_TO_LOG10             0.2302585093
#define LOG10CONST              0.43429448190325182765112891891661F

#ifdef __INTEL_COMPILER
// #define inline __inline  // ipazarba: inline is a keyword. For your custom macros, use CAPITAL letters or a prefix or a suffix.
# define CC_INLINE __inline
#else
# define CC_INLINE inline  // in fact, gcc, too, supports __inline
#endif

//#define _abs(a)          abs(a)
//#define _lmbd(a, b)      (((b) == (0)) ? (32) : (31 - _bit_scan_reverse(b)))
//__forceinline int _norm(unsigned int x) 
//{
//  return(_lmbd(0,abs(x + (x >> 31)))-1); 
//}

static CC_INLINE int _sadd(int src1, int src2)
{
  int result = src1 + src2;
  if( ((src1 ^ src2) & INT_MIN) == 0)
    if((result ^ src1) & INT_MIN)
      result = (src1<0) ? INT_MIN : INT_MAX;
  return(result);
}

static CC_INLINE int _smpy(int src1, int src2)
{
  int result = ((short)src1 * (short)src2) << 1;
  if (result != INT_MIN)	
    return(result);
  else 
    return(INT_MAX);
}

#define _sshl(a, c) ((a)<<(c))

#define _add2(a, b)      c_add2(a, b)
#define _extu(a, b, c)   c_extu(a, b, c)
#define _mpy(a, b)       c_mpy(a, b)
#define _mpyh(a, b)      c_mpyh(a, b)
#define _mpyhl(a, b)     c_mpyhl(a, b)
#define _mpylh(a, b)     c_mpylh(a, b)
#define _smpyh(a, b)     c_smpyh(a, b)
#define _ssub(a, b)      c_ssub(a, b)


static CC_INLINE int32_t c_add2(int32_t a, int32_t b)
{
  return a + b - (((a & ~0x10000) + (b & ~0x10000)) & 0x10000);
}

static CC_INLINE uint32_t c_extu(uint32_t a, uint32_t b, uint32_t c)
{
  return (a << b ) >> c;
}

static CC_INLINE int32_t c_mpy(int32_t a, int32_t b)
{
  return (int16_t)a * (int16_t)b;
}

static CC_INLINE int32_t c_mpyh(int32_t a, int32_t b)
{
  return (a >> 16) * (b >> 16);
}

static CC_INLINE int32_t c_mpyhl(int32_t a, int32_t b)
{
  return (a >> 16) * (int16_t)b;
}

static CC_INLINE int32_t c_mpylh(int32_t a, int32_t b)
{
  return (int16_t)a * (b >> 16);
}

static CC_INLINE int32_t c_smpyh(int32_t a,int32_t b)
{
  int32_t result;

  result = ( (int16_t)((a & ~0xffff) >> 16) *
     (int16_t)((b & ~0xffff) >> 16)) << 1;

  if ( result != (int32_t)0x80000000)
    return(result);
  else {
    return(0x7fffffff);
  }
}

static CC_INLINE int32_t c_ssub(int32_t a, int32_t b)
{
  return (a ^ b) < 0 && ((a - b) ^ a) < 0 ? (int32_t)0x80000000 - (a >= 0) : a - b;
}



typedef struct
{
    int K;
    float expod[8];
    double q;
    float qinv;
    unsigned int Kmask;
} fastexpcnftype;

typedef struct
{
    int K;
    float expod[8];
    double q;
    float qinv;
    unsigned int Kmask;
} db2lincnftype;

/* Inlined functions: */
/* static inline float fastlog(float x); */
/* static inline float fastexp(float x); */
/* static inline float fastdiv(float a,float b); */
/* static inline float fasterdiv(float a,float b); */
/* static inline float fastsqrtf(float x); */
/* static inline float specExp(float x); */
/* static inline float dB2Lin(float x); */

void fastloginit(void);
void fastexpinit(void);
void dB2Lininit(void);
float reciproc_sqrtf(float x);
void triginit(void);
void pol2cart(float * polar,float * out,int num);
void pol2cartClib(float * polar,float * out,int num);
float expIntQ(int x);
void mathfunInit(void);
int floatvec_validnumbers(float *x, int num);
int floatvec_validrange(float *x, int num);
void calcAutocorr(float * X,int n0,int N,int m0,int M,float * r);


#define fastdiv(a,b) ((a)/(b))
#define fasterdiv(a,b) ((a)/(b))
#define fastlog(x) (logf(x))
#define _fabsf(x) (fabsf(x))
#define reciproc_sqrtf(x) (1.0f/sqrtf(x))
#define fastexp(x) (expf(x))

#ifdef __INTEL_COMPILER
#define _lmbd(a, b)      (((b) == (0)) ? (32) : (31 - _bit_scan_reverse(b)))
#else
#define _lmbd(a, b)      (((b) == (0)) ? (32) : (31 - __builtin_clz(b)))
#endif /* __INTEL_COMPILER */

#define _rcpsp(x) (1.0f/(x))
#define _rsqrsp(x) (1.0f/sqrtf(x))

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))

#define _lo(a) (taa_lo(a))
#define _hi(a) (taa_hi(a))
#define _itof(x) (taa_itof(x))
#define _amemd8(p) (*(uint64_t *) (p))

static CC_INLINE int32_t taa_lo(uint64_t a)
{
    union
    {
        uint64_t dval;
        int32_t uval[2];
    }u;
    
    u.dval = a;
#ifdef PROC_LITTLE_ENDIAN
    return u.uval[0];
#else
    return u.uval[1];
#endif
}

static CC_INLINE int32_t taa_hi(uint64_t a)
{
    union
    {
        uint64_t dval;
        int32_t uval[2];
    }u;
    
    u.dval = a;
    
#ifdef PROC_LITTLE_ENDIAN
    return u.uval[1];
#else
    return u.uval[0];
#endif
}

static  float taa_itof(uint32_t x)
{
    union {float f; uint32_t u;} u_src;
    
    u_src.u = x;
    
    return (u_src.f);
}

#endif
