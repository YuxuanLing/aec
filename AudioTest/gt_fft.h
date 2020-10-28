/*****************************************************************************
 * ---------------------------------------------------------------------------
 *                                gt_fft.h
 * ---------------------------------------------------------------------------
 *
 *  Author        : MAH
 *
 *  Description   : Header file for real and complex Good-Thomas FFT
 *                  implementations for some selected lengths, N.
 *
 *****************************************************************************/
#ifndef GT_FFT_H
#define GT_FFT_H

//#ifdef _TMS320C6X
#define GT_ACCESSTYPE double
//#else
//#define GT_ACCESSTYPE float
 //#endif

/* -------- gtNNN() ---------------------------------------------+
 |                                                               |
 |  Complex fft's of size NNN.                                   |
 |                                                               |
 |  Implementations that use a temporary array (t).              |
 |                                                               |
 |                                                               |
 |  x[]    --  input and output array, containing 2*NNN  floats  |
 |  t[]    --  temporary array, also containing 2*NNN  floats    |
 |                                                               |
 +-------------------------------------------------------------- */

void gt60    (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);



void gt80    (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);
void gt160   (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);
void gt240   (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);
void gt320   (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);
void gt384   (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);
void gt480   (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);
void gt768   (GT_ACCESSTYPE *x, GT_ACCESSTYPE *t);


/* -------- gtNNNip() -------------------------------------------+
 |                                                               |
 |  Complex fft's of size NNN.                                   |
 |                                                               |
 |  Implementations that are true in-place.                      |
 |                                                               |
 |                                                               |
 |  x[]    --  input and output array, containing 2*NNN  floats  |
 |                                                               |
 +-------------------------------------------------------------- */

void gt60ip  (GT_ACCESSTYPE *x);
void gt240ip (GT_ACCESSTYPE *x);
void gt320ip (GT_ACCESSTYPE *x);
void gt480ip (GT_ACCESSTYPE *x);



/* -------- gt_real() -------------------------------------------+
 |                                                               |
 |  Wrapper for real fft of size N using compex fft of size N/2  |
 |                                                               |
 +-------------------------------------------------------------- */

float gt_real(  int    N,             /* size of real-valued FFT                 */
                GT_ACCESSTYPE * x,   /* DW-aligned input/output array, N floats */
                GT_ACCESSTYPE * t,   /* DW-aligned intermediate array, N floats */
                GT_ACCESSTYPE * twx, /* DW-aligned twiddle factors, N/2 floats  */
                void   *cplx_fft,     /* pointer to complex FFT of size N/2      */
                int    output_in_t    /* flag to tell whether cplx_fft's output  */
                                      /* appear in t[] (=1) or in x[] (=0)       */
                );

void gt_real_inv(
                int    N,             /* size of real-valued FFT                 */
                GT_ACCESSTYPE * x,   /* DW-aligned input/output array, N floats */
                GT_ACCESSTYPE * t,   /* DW-aligned intermediate array, N floats */
                GT_ACCESSTYPE * twx, /* DW-aligned twiddle factors, N/2 floats  */
                void   *cplx_fft,     /* pointer to complex FFT of size N/2      */
                int    output_in_t,   /* flag to tell whether cplx_fft's output  */
                                      /* appears in t[] (=1) or in x[] (=0)      */
                float xf_N            /* "in[N]"                                 */
                );


void gt_real_inv_noDivide(
                int    N,                    /* size of real-valued FFT                 */
                GT_ACCESSTYPE * x,   /* DW-aligned input/output array, N floats */
                GT_ACCESSTYPE * t,   /* DW-aligned intermediate array, N floats */
                GT_ACCESSTYPE * twx, /* DW-aligned twiddle factors, N/2 floats  */
                void   *cplx_fft,            /* pointer to complex FFT of size N/2      */
                int    output_in_t,          /* flag to tell whether cplx_fft's output  */
                                             /* appears in t[] (=1) or in x[] (=0)      */
                float xf_N                   /* "in[N]"                                 */
                );



/* -------- tw_gen_realfft() ------------------------------------+
 |                                                               |
 |  Initialization of extra twiddle factors for gt_real()        |
 |                                                               |
 +-------------------------------------------------------------- */

void tw_gen_realfft(
                        float *twc, /* (O) Array of N/2 floats to be initialized */
                        int N       /* (I) # of points in real FFT               */
                    );

void gt240_loadTables(void *pTableMem);
void gt240ip_loadTables(void *pTableMem);

#endif
/*****************************************************************************
 *                          end of gt_fft.h
 ****************************************************************************/
