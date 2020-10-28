/*****************************************************************************
 * ---------------------------------------------------------------------------
 *                                fft_real
 * ---------------------------------------------------------------------------
 *
 *  Author        : MAH
 *
 *  Description   : Forward Fast Fourier Transform (FFT) of pure real input.
 *
 *  This routine computes the FFT of a pure real-valued signal. The input is of
 *  the form [re0, re1, ..., reN-1], i.e there are no interleaved zeros for the
 *  non-existent imaginary data. It computes only the first half of the output
 *  + element number "N", since it is assumed that this is the only part that
 *  the caller is going to be interested in. The second half would be redundant
 *  given the real input.
 *
 *  To accomplish the N-pt. real-valued FFT, we use a complex N/2-pt. FFT plus
 *  one additional twiddle stage. Hence, two sets of twiddle factors are being
 *  used: one for the TI FFT and one for the add'l twiddle stage.
 *
 *  The TI DSPLIB routine is not an in-place routine and consequently the user
 *  has to provide three pointers for data: i) input to TI FFT ii) output from
 *  TI FFT and iii) output from the additional twiddle stage. i) and ii) must
 *  not be the same physical address, but iii) can be the same as either of i)
 *  and ii), or it can be a different buffer altogether. For best cache
 *  performance, the latter option is really only recommended if the output is
 *  not going to be read again soon after FFT routine has completed, or if the
 *  output is to be directly DMAed out wihout further modification.
 *
 *  NOTE:
 *  -----
 *    o  The input buffer in[] is overwritten with garbage data, unless
 *       out[] = in[], in which case it will be overwritten with useful data
 *
 *    o  twiddleTI[] should be initialized using sp_tw_gen(twiddleTI, N/2);
 *
 *    o  twiddleXtra[] should be initialized using
 *       tw_gen_realfft(twiddleXtra, N);
 *
 *    o  Some buffers are assumed to be double-word aligned, see parameter
 *       list.
 *
 *    o  out[] may be either *exactly* the same as in[] or interm[], or it
 *       can be a different array.
 *
 *    o  in[] and interm[] must be distinct, non-overlapping arrays.
 *
 *    o  This routine was optimized for 64 =< N =< 256, but works for
 *       _REALFFT_MIN_N =< N =< _REALFFT_MAX_N. Expect a certain amount
 *       of L1D thrashing for N >= 512.
 *
 *    o  Radix for the TI FFT call is computed at run time.
 *
 *  Returns       : "out[N]", which doesn't fit in out[].
 *
*****************************************************************************/

#include <math.h> /* For twiddle factor generation */
#include "gt_fft.h" /*  */
#include "fftutil.h"
#include "mathfun.h"
#include "realfft.h"

/* --- Defining the range of supported N-values --- */
#define _REALFFT_MIN_N  64
#define _REALFFT_MAX_N  1024

#define NOASSUME


/*****************************************************************************
 * ---------------------------------------------------------------------------
 *                       ifft_real / ifft_real_noDivide
 * ---------------------------------------------------------------------------
 *
 *  Author        : MAH
 *
 *  Description   : Inverse Fast Fourier Transform (IFFT) whose output is
 *                  known to be real beforehand.
 *
 *  This routine computes the IFFT of an FFT of a pure real-valued signal.
 *  Due to symmetry reasons, it sufficies to know re0, im0, re1, im1,..
 *  .., reN/2-1, imN/2-1, reN/2 to compute the IFFT given the result be
 *  real, and hence these are the only values that the routine take as
 *  input. The last value (reN/2) is passed as a separate argument to keep the
 *  input array even-sized. The output is on the form [re0, re1, ..., reN-1],
 *  i.e there are no interleaved zeros for the non-existent imaginary data.
 *
 *  To accomplish the N-pt. real-valued FFT, we use a complex N/2-pt. FFT plus
 *  some additional twiddle stages. Hence, two sets of twiddle factors are being
 *  used: one for the TI FFT and one for the additional twiddle stages.
 *
 *  The TI DSPLIB routine is not an in-place routine and consequently the user
 *  has to provide three pointers for data: i) input to TI FFT ii) output from
 *  TI FFT and iii) output from the additional twiddle stage. i) and ii) must
 *  not be the same physical address, but iii) can be the same as either of i)
 *  and ii), or it can be a different buffer altogether. For best cache
 *  performance, the latter option is really only recommended if the output is
 *  not going to be read again soon after FFT routine has completed , or if the
 *  output is to be directly DMAed out wihout further modification (which
 *  should be a fairly unlikely scenario). Recommended synopsis for
 *  ifft_real() is that 'out == interm', and for ifft_real_noDivide() it
 *  is recommended that either 'out == interm', or that 'out == in'.
 *
 *  NOTE:
 *  -----
 *    o  The input buffer in[] is overwritten with garbage data, unless
 *       out[] = in[], in which case it will be overwritten with useful data
 *
 *    o  twiddleTI[] should be initialized using sp_tw_gen(twiddleTI, N/2);
 *
 *    o  twiddleXtra[] should be initialized using
 *       tw_gen_realfft(twiddleXtra, N);
 *
 *    o  Some buffers are assumed to be double-word aligned, see parameter
 *       list.
 *
 *    o  out[] may be either *exactly* the same as in[] or interm[], or it
 *       can be a different array.
 *
 *    o  in[] and interm[] must be distinct, non-overlapping arrays.
 *
 *    o  This routine was optimized for 64 =< N =< 256, but works for
 *       _REALFFT_MIN_N =< N =< _REALFFT_MAX_N. Expect a certain amount
 *       of L1D thrashing for N >= 512.
 *
 *    o  Radix for the TI FFT call is computed at run time.
 *
 *****************************************************************************/
float fft_real768(
	            const float * in,          /* DW-aligned input array, N floats */
				const float *twiddleTI,            /* DW-aligned tw. factors for TI FFT., N floats */
				float * interm,            /* DW-aligned intermediate array, N floats */
				const unsigned char *brev,         /* bit reversal table, see TI DSPLIB doc. */
				const float * twiddleXtra, /* DW-aligned tw. factors for the add'l stage. */
				float * out                /* DW-aligned output array, size N floats */
				)
{
	/* int rad; */
        int k;
	float lastval;
	int N = 768;    /* hardcoded version... */

#ifndef NOASSUME
	_nassert(((int) in          & 0x7) == 0);
	_nassert(((int) twiddleTI   & 0x7) == 0);
	_nassert(((int) interm      & 0x7) == 0);
	_nassert(((int) twiddleXtra & 0x7) == 0);
	_nassert(((int) out         & 0x7) == 0);
#endif

    (void)twiddleTI;
    (void)brev;

	/* determine the radix of N/2 */
	/* rad = _lmbd(1, N) & 0x1 ? 2 : 4; */

	/* First, call the TI DSPLIB routine for forward FFT of size N/2 */
        //	DSPF_sp_fftSPxSP(N/2, in, twiddleTI, interm, brev, rad, 0, N/2);

        gt384((GT_ACCESSTYPE *) in,(GT_ACCESSTYPE *) interm);
        for(k = 0;k<384;k++)
        {
            interm[k] = in[k];
        }



	/* Then do the last twiddles */
	lastval = interm[0] -interm[1];
	out[0]  = interm[0] + interm[1];
	out[1]  = 0;

	for(k=1; k<N/4; k++)
	{
		float tw_c = twiddleXtra[    k];
        float tw_s = twiddleXtra[N/4-k];

		float a = interm[2*k  ];
		float b = interm[2*k+1];
		float c = interm[2*(N/2-k)  ];
		float d = interm[2*(N/2-k)+1];

		out[2*k  ] = (a+c)*.5f + tw_c*(b+d) - tw_s*(a-c);
		out[2*k+1] = (b-d)*.5f - tw_c*(a-c) - tw_s*(b+d);

		out[2*(N/2-k)  ] = (a+c)*.5f - tw_c*(b+d) + tw_s*(a-c);
		out[2*(N/2-k)+1] = (d-b)*.5f - tw_c*(a-c) - tw_s*(b+d);
	}

	out[2*N/4]   =  interm[2*N/4];
	out[2*N/4+1] = -interm[2*N/4+1];

	return (lastval); /* ...and by the way, we give you "out[2*N/2]" as well! */
}

void ifft_real_noDivide768(
                float * in,                /* DW-aligned input array, N floats */
				const float * twiddleTI,   /* DW-aligned tw. factors for TI FFT, N floats */
				float * interm,            /* DW-aligned intermediate array, N floats */
				const unsigned char *brev,         /* bit reversal table, see TI DSPLIB doc. */
				const float * twiddleXtra, /* DW-aligned tw. factors for the add'l stage, N/2 floats */
				float * out,               /* DW-aligned output array, size N floats */
				float lastval                      /* "in[N]" */
				)
{
	/*int rad */
        int k;
	float i0, i1;
	int N = 768;    /* hardcoded version... */
	float * _inf = &in[  2];
	float * _inr = &in[N-2];


#ifndef NOASSUME
	_nassert(((int) in          & 0x7) == 0);
	_nassert(((int) twiddleTI   & 0x7) == 0);
	_nassert(((int) interm      & 0x7) == 0);
	_nassert(((int) twiddleXtra & 0x7) == 0);
	_nassert(((int) out         & 0x7) == 0);
#endif

    (void)twiddleTI;
    (void)brev;

        //	rad = 2;
	/* rad = _lmbd(1, N) & 0x1 ? 2 : 4; */

	/*
	 * First do the twiddle stage. Similar to the corresponding
	 * stage in ifft_real() but it also conjugates the input
	 */
    i0 = (in[0] - in[1] + lastval) * .5f;
    i1 = (in[0] + in[1] - lastval) * .5f;

    in[0] =  i0;
    in[1] = -i1;

	for(k=1; k<N/4; k++)
	{
		float tw_c = twiddleXtra[    k];
		float tw_s = twiddleXtra[N/4-k];

		float a = _inf[0], b = _inf[1];
		float c = _inr[0], d = _inr[1];

		_inf[0] =   (a+c)*.5f - tw_c*(b+d) - tw_s*(a-c);
		_inf[1] = -((b-d)*.5f + tw_c*(a-c) - tw_s*(b+d));

		_inr[0] =   (a+c)*.5f + tw_c*(b+d) - tw_s*(c-a);
		_inr[1] = -((d-b)*.5f - tw_c*(c-a) - tw_s*(b+d));

		_inf +=2; _inr -=2;
	}

	//in[2*N/4  ] = in[2*N/4  ];
	//in[2*N/4+1] = in[2*N/4+1];


	/* Then call TI DSPLIB routine for forward FFT of size N/2 */
    //	DSPF_sp_fftSPxSP(N/2, in, twiddleTI, interm, brev, rad, 0, N/2);

        gt384((GT_ACCESSTYPE *) in,(GT_ACCESSTYPE *) interm);
        for(k = 0;k<384;k++)
        {
            interm[k] = in[k];
        }

	/*
	 * And lastly, conjugate the output of the FFT. We also need to
	 * multiply by a factor of 2 since we actually used an FFT of size
	 * N/2 rather than N. It would be somewhat more efficient to do
	 * the mpy-by-2 on the input side instead, but this is a bit
	 * more straight forward and the cycle difference is small
	 */
	for(k=0; k<N/2; k++)
	{
		float im = _itof(_hi(_amemd8(&interm[2*k])));
		float re = _itof(_lo(_amemd8(&interm[2*k])));

		out[2*k  ] =  2.0f * re;
		out[2*k+1] = -2.0f * im;
	}
}


/*****************************************************************************
 * Author        : MAH
 * Description   : Routine for run-time initialization of twiddle factor array
 *				   to support fft_real() described elsewhere.
 ****************************************************************************/

void tw_gen_realfft(
						float *twc,	/* (O) Array of n/2 floats to be initialized */
						int N		/* (I) # of points in real FFT */
					)
{
	int k;

	for(k=0; k<(N/2); k++)
	{
		twc[k] = (float)(cos(2.0 * 3.1415926535898 * k/(double) N) * .5);
	}
}
