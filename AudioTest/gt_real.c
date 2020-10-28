/******************************************************************************* */
/*                                gt_real.c                                      */
/******************************************************************************* */


/* ----------------------------------------------------------------------------  */
/*                                gt_real()                                      */
/* ----------------------------------------------------------------------------  */
/*                                                                               */
/*  Author        : MAH                                                          */
/*                                                                               */
/*  Description   : Forward Fast Fourier Transform (FFT) of pure real input      */
/*                  for some selected lengths, N.                                */
/*                                                                               */
/*  This routine computes the FFT of a pure real-valued signal. The input is     */
/*  of the form [re0, re1, ..., reN-1], thus there are no interleaved zeros      */
/*  for the non-existing imaginary data. It computes only the first half of      */
/*  the output + element number "N", as it is assumed that this is the only      */
/*  part that the caller is going to be interested in keeping in mind that       */
/*  the second half would be redundant given the real-valued input.              */
/*                                                                               */
/*  To accomplish the N-point real-valued FFT, we use a complex N/2-point        */
/*  FFT plus one additional twiddle stage. The user has to provide a pointer     */
/*  to the complex FFT as an argument into the function. This pointer argument   */
/*  is being passed to the function as void *, but is originally of either the   */
/*  type void (*) (double *, double *), or of type void (*) (double *).          */
/*                                                                               */
/*  In the first case, input is taken from the location of the first pointer     */
/*  and output might be written to back to the same location, or to the location */
/*  of the second argument (see comment about output_in_t further down).         */
/*                                                                               */
/*  In the second case the complex FFT is a true in-place routine that requires  */
/*  no other storage than the input/output array.                                */
/*                                                                               */
/*  All data is float even though the pointers are of type double *. We use      */
/*  this as a means to provide alignment info to the compiler.                   */
/*                                                                               */
/*  This routine was primarily intended to being used with Good-Thomas FFT       */
/*  implementations, however any complex fft can be used as long as it           */
/*  interfaces the same way, and as long as the requirements on the size         */
/*  N are adhered to (see notes further down).                                   */
/*                                                                               */
/*  The function is called as:                                                   */
/*                                                                               */
/*  float gt_real(      int N,                                                   */
/*                      double * x,                                      */
/*                      double * t,                                      */
/*                      double * twx,                                    */
/*                      void   *cplx_fft,                                        */
/*                      int    output_in_t                                       */
/*               )                                                               */
/*                                                                               */
/*                                                                               */
/*          N           --  Size of "full" FFT. I.e. for an M point real fft,    */
/*                          N should be set to M and an M/2 pt. complex FFT      */
/*                          routine should be provided via the function          */
/*                          pointer argument further down. It is required        */
/*                          that N be divisible by 4 and also that               */
/*                          _GT_REALFFT_MIN_N <= N <= _GT_REALFFT_MAX_N          */
/*                                                                               */
/*          x           --  Double-word aligned pointer to input or input and    */
/*                          output data, N floats                                */
/*                                                                               */
/*          t           --  DW-aligned output or intermediate array, N floats.   */
/*                          Should be NULL if cplx_fft is a true in-place        */
/*                          function.                                            */
/*                                                                               */
/*          twx         --  DW-aligned array of N/2 floats for the post-fft      */
/*                          stage. Generate with tw_gen_realfft(twx, N);         */
/*                                                                               */
/*          cplx_fft    --  pointer to an N/2 point complex FFT                  */
/*                                                                               */
/*          output_in_t --  flag to tell whether cplx_fft's output appears       */
/*                          in t[] (output_in_t=1) or in x[] (output_in_t=0)     */
/*                                                                               */
/*                                                                               */
/*  Note that t[] serves as an indicator of whether cplx_fft is an in-place      */
/*  routine ((int) t == NULL) or not. On C6x targets, because of calling         */
/*  convention details that it won't be elaborated on any further here, the      */
/*  value of t happens to be of no consequence if a true inplace routine is      */
/*  passed as cplx_fft. However for other targets such (erronous) usage          */
/*  might imply a fatal run-time error. Thus, it is strongly recommended         */
/*  that t indeed be NULL whenever cplx_fft is true in-place (i.e. takes         */
/*  one argument only).                                                          */
/*                                                                               */
/*  The requirement that N be divisible by 4 comes from i) the fact that N/2     */
/*  must be an integer (remember we want to use an N/2 pt complex FFT) and       */
/*  ii) from a manual unroll-by-2 in the post-FFT stage. This unroll could be    */
/*  removed if required to support N divisible by 2 (but not 4).                 */
/*                                                                               */
/*  Please respect the alignment restrictions listed above, as well as the       */
/*  requirements on N.                                                           */
/*                                                                               */
/*                                                                               */
/*  Returns       : "out[N]", which doesn't fit in out[].                        */
/*                                                                               */
/******************************************************************************* */

#include "gt_fft.h"

#include <stdlib.h>

#define _GT_REALFFT_MIN_N  60
#define _GT_REALFFT_MAX_N 960


float gt_real(  int    N,             /* size of real-valued FFT                 */
                double * x,   /* DW-aligned input/output array, N floats */
                double * t,   /* DW-aligned intermediate array, N floats */
                double * twx, /* DW-aligned twiddle factors, N/2 floats  */
                void   *cplx_fft,     /* pointer to complex FFT of size N/2      */
                int    output_in_t    /* flag to tell whether cplx_fft's output  */
                                      /* appears in t[] (=1) or in x[] (=0)      */
                )
{


    /* ============ Call the complex N/2 pt. FFT ====================== */
    {
        union { void *_; void (*__) (double *          ); } cplx_fft_ip;
        union { void *_; void (*__) (double *, double *); } cplx_fft_op;

        cplx_fft_ip._ = cplx_fft_op._ = cplx_fft;

        if   (t == NULL) (*cplx_fft_ip.__) (x   );
        else             (*cplx_fft_op.__) (x, t);
    }


    /* ============ Perform the post-FFT twiddle stage ================ */
    {
        int k;

        float xf_N;
        float * xf   = output_in_t ? (float *) t: (float *) x;
        float * xf_  = output_in_t ? (float *) t: (float *) x;
        float * twxf = (float *) twx;

        xf_N = xf[0] - xf[1]; xf[0] = xf[0] + xf[1]; xf[1] = 0;

        for(k = 1; k < (N/4); k++)
        {

            float tw_c = twxf [      k    ];
            float tw_s = twxf [N/4 - k    ];
            float a    = xf   [    2*k    ];
            float b    = xf   [    2*k + 1];
            float c    = xf   [N - 2*k    ];
            float d    = xf   [N - 2*k + 1];

            xf_[    2*k    ] = (a + c) * .5f + tw_c * (b + d) - tw_s * (a - c);
            xf_[    2*k + 1] = (b - d) * .5f - tw_c * (a - c) - tw_s * (b + d);
            xf_[N - 2*k    ] = (a + c) * .5f - tw_c * (b + d) + tw_s * (a - c);
            xf_[N - 2*k + 1] = (b - d) *-.5f - tw_c * (a - c) - tw_s * (b + d);
        }

        //xf[N/2] =  xf[N/2];
        xf[N/2 + 1] = -xf[N/2 + 1];

        return (xf_N);
    }
}

/* ----------------------------------------------------------------------------  */
/*                                gt_real_inv()                                  */
/* ----------------------------------------------------------------------------  */
/*                                                                               */
/*  Author        : MAH                                                          */
/*                                                                               */
/*  Description   : Inverse Forward Fast Fourier Transform (IFFT) of complex     */
/*                  input that are such that the IFFT of it is purely real-      */
/*                  valued. For some selected lengths, N.                        */
/*                                                                               */
/*  This routine computes a real IFFT. The input is of the form [re0, 0, re1,    */
/*  im2, ..., reN/2-1, imN/2-1]. In addition, reN/2 is passed outside of the     */
/*  array. It is not as complicated as it sounds - just have a look at gt_real() */
/*  and things will get clearer.                                                 */
/*                                                                               */
/*  The output is a real-valued array of size N floats.                          */
/*                                                                               */
/*  To accomplish the N-point real-valued IFFT, we use a complex N/2-point       */
/*  FFT plus a twiddle stage before and a conjugation step afterwards. The user  */
/*  has to provide a pointer to the complex FFT as an argument into the function.*/
/*  This pointer argument is being passed to the function as void *, but is      */
/*  originally of either the  type void (*) (double *, double *), or of type     */
/*  void (*) (double *).                                                         */
/*                                                                               */
/*  In the first case, input is taken from the location of the first pointer     */
/*  and output might be written to back to the same location, or to the location */
/*  of the second argument (see comment about output_in_t further down).         */
/*                                                                               */
/*  In the second case the complex FFT is a true in-place routine that requires  */
/*  no other storage than the input/output array.                                */
/*                                                                               */
/*  All data is float even though the pointers are of type double *. We use      */
/*  this as a means to provide alignment info to the compiler.                   */
/*                                                                               */
/*  This routine was primarily intended to being used with Good-Thomas FFT       */
/*  implementations, however any complex fft can be used as long as it           */
/*  interfaces the same way, and as long as the requirements on the size         */
/*  N are adhered to (see notes further down).                                   */
/*                                                                               */
/*  The function is called as:                                                   */
/*                                                                               */
/*    void gt_real_inv(                                                          */
/*                       int    N,                                               */
/*                       double * x,                                     */
/*                       double * t,                                     */
/*                       double * twx,                                   */
/*                       void   *cplx_fft,                                       */
/*                       int    output_in_t,                                     */
/*                       float xf_N                                              */
/*                    )                                                          */
/*                                                                               */
/*          N           --  Size of "full" FFT. I.e. for an M point real fft,    */
/*                          N should be set to M and an M/2 pt. complex FFT      */
/*                          routine should be provided via the function          */
/*                          pointer argument further down. It is required        */
/*                          that N be divisible by 4 and also that               */
/*                          _GT_REALFFT_MIN_N <= N <= _GT_REALFFT_MAX_N          */
/*                                                                               */
/*          x           --  Double-word aligned pointer to input or input and    */
/*                          output data, N floats                                */
/*                                                                               */
/*          t           --  DW-aligned output or intermediate array, N floats.   */
/*                          Should be NULL if cplx_fft is a true in-place        */
/*                          function.                                            */
/*                                                                               */
/*          twx         --  DW-aligned array of N/2 floats for the post-fft      */
/*                          stage. Generate with tw_gen_realfft(twx, N);         */
/*                                                                               */
/*          cplx_fft    --  pointer to an N/2 point complex FFT                  */
/*                                                                               */
/*          output_in_t --  flag to tell whether cplx_fft's output appears       */
/*                          in t[] (output_in_t=1) or in x[] (output_in_t=0)     */
/*                                                                               */
/*          float xf_N  --  "in[N]", which did not fit in x[]                    */
/*                                                                               */
/*  Note that t[] serves as an indicator of whether cplx_fft is an in-place      */
/*  routine ((int) t == NULL) or not. On C6x targets, because of calling         */
/*  convention details that it won't be elaborated on any further here, the      */
/*  value of t happens to be of no consequence if a true inplace routine is      */
/*  passed as cplx_fft. However for other targets such (erronous) usage          */
/*  might imply a fatal run-time error. Thus, it is strongly recommended         */
/*  that t indeed be NULL whenever cplx_fft is true in-place (i.e. takes         */
/*  one argument only).                                                          */
/*                                                                               */
/*  The requirement that N be divisible by 4 comes from i) the fact that N/2     */
/*  must be an integer (remember we want to use an N/2 pt complex FFT) and       */
/*  ii) from a manual unroll-by-2 in the post-FFT stage. This unroll could be    */
/*  removed if required to support N divisible by 2 (but not 4).                 */
/*                                                                               */
/*  Please respect the alignment restrictions listed above, as well as the       */
/*  requirements on N.                                                           */
/*                                                                               */
/******************************************************************************* */

void gt_real_inv(
                int    N,             /* size of real-valued FFT                 */
                double *x,   /* DW-aligned input/output array, N floats */
                double * t,   /* DW-aligned intermediate array, N floats */
                double * twx, /* DW-aligned twiddle factors, N/2 floats  */
                void   *cplx_fft,     /* pointer to complex FFT of size N/2      */
                int    output_in_t,   /* flag to tell whether cplx_fft's output  */
                                      /* appears in t[] (=1) or in x[] (=0)      */
                float xf_N            /* "in[N]"                                 */
                )
{
    int   k;
    float i0, i1;
    float * xf   = (float *) x;
    float * xf_  = (float *) x;
    float * twxf = (float *) twx;
    float twoOverN;

    /* Pre-calced 2/N for some commonly used values of N */
    switch (N)
    {
        case 640:   twoOverN = 0.00312500000f;  break;
        case 320:   twoOverN = 0.00625000000f;  break;
        case 160:   twoOverN = 0.01250000000f;  break;
        default:    twoOverN = 2.0f/(float) N;  break;
    }


    /*
     * First do the twiddle stage. We also scale by 2/N to compensate for that
     * we use a forward fft instead of an ifft (the factor 1/N) and the fact
     * that we use an fft of half the size (the factor 2). We also conjugate
     * the input, again because we use a forward fft.
     */
    i0 = (xf[0] - xf[1] + xf_N) * .5f;
    i1 = (xf[0] + xf[1] - xf_N) * .5f;

    xf[0] =  i0 * twoOverN;
    xf[1] = -i1 * twoOverN;

    for(k=1; k<N/4; k++)
    {
        float tw_c = twxf [      k    ];
        float tw_s = twxf [N/4 - k    ];
        float a    = xf   [    2*k    ];
        float b    = xf   [    2*k + 1];
        float c    = xf   [N - 2*k    ];
        float d    = xf   [N - 2*k + 1];

        xf_[    2*k    ] =  ((a + c) * .5f - tw_c * (b + d) - tw_s * (a - c)) * twoOverN;
        xf_[    2*k + 1] =  ((b - d) *-.5f - tw_c * (a - c) + tw_s * (b + d)) * twoOverN;
        xf_[N - 2*k    ] =  ((a + c) * .5f + tw_c * (b + d) + tw_s * (a - c)) * twoOverN;
        xf_[N - 2*k + 1] =  ((b - d) * .5f - tw_c * (a - c) + tw_s * (b + d)) * twoOverN;

    }

    xf[2*N/4  ] = xf[2*N/4  ] * twoOverN;
    xf[2*N/4+1] = xf[2*N/4+1] * twoOverN;


     /* ============ Call the complex N/2 pt. FFT ====================== */
     {
         union { void *_; void (*__) (double *          ); } cplx_fft_ip;
         union { void *_; void (*__) (double *, double *); } cplx_fft_op;

         cplx_fft_ip._ = cplx_fft_op._ = cplx_fft;

         if   (t == NULL) (*cplx_fft_ip.__) (x   );
         else             (*cplx_fft_op.__) (x, t);
    }


    /* ============ Perform the post-FFT conjugation =================== */
    {

        float * y   = output_in_t ? (float *) t: (float *) x;

        for(k=0; k<N/2; k++)
        {
            y[2*k + 1] = -y[2*k + 1];
        }

    }
}


/* ----------------------------------------------------------------------------  */
/*                                gt_real_inv_noDivide()                         */
/* ----------------------------------------------------------------------------  */
/*                                                                               */
/*  Author        : MAH (JPS)                                                    */
/*                                                                               */
/*  Description   : Same as gt_real_inv, but without division by fftsize N.      */
/*                  The function is used by synthesize filter                    */
/*                                                                               */
/*  The function is called as:                                                   */
/*                                                                               */
/*    void gt_real_inv_noDivide(                                                 */
/*                       int    N,                                               */
/*                       double * x,                                     */
/*                       double * t,                                     */
/*                       double * twx,                                   */
/*                       void   *cplx_fft,                                       */
/*                       int    output_in_t,                                     */
/*                       float xf_N                                              */
/*                    )                                                          */
/*                                                                               */
/*          N           --  Size of "full" FFT. I.e. for an M point real fft,    */
/*                          N should be set to M and an M/2 pt. complex FFT      */
/*                          routine should be provided via the function          */
/*                          pointer argument further down. It is required        */
/*                          that N be divisible by 4 and also that               */
/*                          _GT_REALFFT_MIN_N <= N <= _GT_REALFFT_MAX_N          */
/*                                                                               */
/*          x           --  Double-word aligned pointer to input or input and    */
/*                          output data, N floats                                */
/*                                                                               */
/*          t           --  DW-aligned output or intermediate array, N floats.   */
/*                          Should be NULL if cplx_fft is a true in-place        */
/*                          function.                                            */
/*                                                                               */
/*          twx         --  DW-aligned array of N/2 floats for the post-fft      */
/*                          stage. Generate with tw_gen_realfft(twx, N);         */
/*                                                                               */
/*          cplx_fft    --  pointer to an N/2 point complex FFT                  */
/*                                                                               */
/*          output_in_t --  flag to tell whether cplx_fft's output appears       */
/*                          in t[] (output_in_t=1) or in x[] (output_in_t=0)     */
/*                                                                               */
/*          float xf_N  --  "in[N]", which did not fit in x[]                    */
/*                                                                               */
/*  Note that t[] serves as an indicator of whether cplx_fft is an in-place      */
/*  routine ((int) t == NULL) or not. On C6x targets, because of calling         */
/*  convention details that it won't be elaborated on any further here, the      */
/*  value of t happens to be of no consequence if a true inplace routine is      */
/*  passed as cplx_fft. However for other targets such (erronous) usage          */
/*  might imply a fatal run-time error. Thus, it is strongly recommended         */
/*  that t indeed be NULL whenever cplx_fft is true in-place (i.e. takes         */
/*  one argument only).                                                          */
/*                                                                               */
/*  The requirement that N be divisible by 4 comes from i) the fact that N/2     */
/*  must be an integer (remember we want to use an N/2 pt complex FFT) and       */
/*  ii) from a manual unroll-by-2 in the post-FFT stage. This unroll could be    */
/*  removed if required to support N divisible by 2 (but not 4).                 */
/*                                                                               */
/*  Please respect the alignment restrictions listed above, as well as the       */
/*  requirements on N.                                                           */
/*                                                                               */
/******************************************************************************* */
void gt_real_inv_noDivide(
                int    N,             /* size of real-valued FFT                 */
                double *x,   /* DW-aligned input/output array, N floats */
                double * t,   /* DW-aligned intermediate array, N floats */
                double * twx, /* DW-aligned twiddle factors, N/2 floats  */
                void   *cplx_fft,     /* pointer to complex FFT of size N/2      */
                int    output_in_t,   /* flag to tell whether cplx_fft's output  */
                                      /* appears in t[] (=1) or in x[] (=0)      */
                float xf_N            /* "in[N]"                                 */
                )
{
    int   k;
    float i0, i1;
    float * xf   = (float *) x;
    float * xf_  = (float *) x;
    float * twxf = (float *) twx;
    float two = 2.0f;


     /*
     * First do the twiddle stage.
     * We use a forward fft instead of an ifft (the factor 1/N), but
     * we do not divide by N because the synthesize filter doesn't want this compensation.
     * We multiply with 2 because we use an fft of half the size. We also conjugate
     * the input, again because we use a forward fft.
     */
    i0 = (xf[0] - xf[1] + xf_N) * .5f;
    i1 = (xf[0] + xf[1] - xf_N) * .5f;

    xf[0] =  i0 * two;
    xf[1] = -i1 * two;

    for(k=1; k<N/4; k++)
    {
        float tw_c = twxf [      k    ];
        float tw_s = twxf [N/4 - k    ];
        float a    = xf   [    2*k    ];
        float b    = xf   [    2*k + 1];
        float c    = xf   [N - 2*k    ];
        float d    = xf   [N - 2*k + 1];

        xf_[    2*k    ] =  ((a + c) * .5f - tw_c * (b + d) - tw_s * (a - c)) * two;
        xf_[    2*k + 1] =  ((b - d) *-.5f - tw_c * (a - c) + tw_s * (b + d)) * two;
        xf_[N - 2*k    ] =  ((a + c) * .5f + tw_c * (b + d) + tw_s * (a - c)) * two;
        xf_[N - 2*k + 1] =  ((b - d) * .5f - tw_c * (a - c) + tw_s * (b + d)) * two;

    }

    xf[2*N/4  ] = xf[2*N/4  ] * two;
    xf[2*N/4+1] = xf[2*N/4+1] * two;


     /* ============ Call the complex N/2 pt. FFT ====================== */
     {
         union { void *_; void (*__) (double *          ); } cplx_fft_ip;
         union { void *_; void (*__) (double *, double *); } cplx_fft_op;

         cplx_fft_ip._ = cplx_fft_op._ = cplx_fft;

         if   (t == NULL) (*cplx_fft_ip.__) (x   );
         else             (*cplx_fft_op.__) (x, t);
    }


    /* ============ Perform the post-FFT conjugation =================== */
    {

        float * y   = output_in_t ? (float *) t: (float *) x;

        for(k=0; k<N/2; k++)
        {
            y[2*k + 1] = -y[2*k + 1];
        }

    }
}




/* ============================================================================= */
/*                               end of gt_real.c                                */
/* ============================================================================= */
