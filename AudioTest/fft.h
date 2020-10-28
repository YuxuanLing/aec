#ifndef FFTHEADER
#define FFTHEADER

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 * Description   : This routine must be called before using fft_process
 *
 *
 * Returns       : void
 ****************************************************************************/
void fft_init(void);

/*****************************************************************************
 * Description   : fft_process performs a comlex n-point FFT on the data x.
 *                 Data are organised like this.
 *                 x[0] - real part data1
 *                 x[1] - imaginary part data1
 *                 x[2] - real part data2
 *                 x[3] - imaginary part data2
 *                     ...
 *                     ...
 *
 * Parameters    : x - Data input and output.
 *                 n - Number of complex indata (n =< 256)
 *
 * Returns       : void
 ****************************************************************************/
void fft_process128(float *x);


/*****************************************************************************
 * Description   : real fft process for fftsizes 256 and 768, when input is
 *                 real,
 *                     ...
 *                     ...
 *
 * Parameters    : x - Data input and output.
 *                 n - Number of complex indata (n =< 256)
 *
 * Returns       : float representing N/2 value of fft when n==768, else 0.0f
 ****************************************************************************/
float realfft_process(float *x, double *scratchpad, short n);


/*****************************************************************************
 * Description   : real ifft process for fftsizes 256 and 768, when input is
 *                 complex and result is real
 *                     ...
 *                     ...
 *
 * Parameters    : x       - Data input and output.
 *                 lastval - last N/2 value from realfft_process (use 0.0f if
 *                           not known)
 *                 n - Number of complex indata (n == 768 || n== 256)
 *
 * Returns       : void
 ****************************************************************************/

void realIfft_process(float *x, double *scratchpad, float lastval, short n);

#ifdef __cplusplus
}
#endif

#endif
