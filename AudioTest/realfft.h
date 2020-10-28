void tw_gen_realfft(float *twc, int N);



/* ------------------------------------------------------------------------------------------------------
			faster versions that only work for N=768
   ---------------------------------------------------------------------------------------------------- */

float fft_real768(
				const float *ptr_x,			/* pointer to a input array of size N */
				const float *ptr_w,         /* pointer to a input array of size N */
				float *ptr_y,               /* pointer to a input array of size N */
				const unsigned char *brev,	/* bit reversal look-up table, see TI DSPLIB doc. */
				const float *twc,			/* pointer to twiddles for the add'l stage. DW-aligned. */
				float *out                  /* pointer to output array */
				);


void ifft_real_noDivide768(
                float *in,                  /* DW-aligned input array, N floats */
				const float *twiddleTI,     /* DW-aligned tw. factors for TI FFT, N floats */
				float *interm,              /* DW-aligned intermediate array, N floats */
				const unsigned char *brev,  /* bit reversal table, see TI DSPLIB doc. */
				const float *twiddleXtra,   /* DW-aligned tw. factors for the add'l stage, N/2 floats */
				float *out,                 /* DW-aligned output array, size N floats */
				float lastval               /* "in[N]" */
				);
