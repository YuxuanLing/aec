/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Switches      : <COMPILER SWITCH HERE>
 *                    <add description here>
 *
 * Description   : <add description here>
 *
 * Note          : <notes here>
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *                 Modification made.
 *
 * yyyy-mm-dd    : <signature>
 *                 File created.
 *
 **************************************************************************/
#define _MATT_

#include "fft.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "realfft.h"
#include "gt_fft.h"
#include "fftutil.h"

#define FFTSIZE 256

static __declspec(align(64)) float tw768_real[768 / 2];

static unsigned char digit_rev_table[FFTSIZE]; /* Assumes FFTSIZE <= 256 */

/* Function prototypes */
static void fft_gen_digit_rev_table(unsigned char *digit_rev_table, short n);

void fft_init(void)
{
    /* generate digit reversal table. */
    fft_gen_digit_rev_table(digit_rev_table, FFTSIZE);
    //#ifdef _MATT_ //need those for reall fft
    tw_gen_realfft(tw768_real, 768);
    //#endif
}

/* The follwoing C code is used to generate the digit-reverse table */
/* unsigned char digit_rev_table assumes n <= 256 */

static void fft_gen_digit_rev_table(unsigned char *digit_rev_table, short n)
{
    int i, j, k, l, m;
    int cnt = 0;
    k = n;
    m = 0;
    while(k != 4)
    {
        k >>= 2;
        m += 2;
    }


    for(i = 0; i < n; i++)
    {
        j = 0;
        k = i;

        for(l = m; k != 0; l -= 2)
        {
            j += (k & 3) << l;
            k >>= 2;
            cnt++;
        }

        digit_rev_table[i] = (unsigned char)j;
    }

}

/* The follwoing C code is used to compute the fft*/
/* unsigned char digit_rev_table assumes n <= 256 */

/*****************************************************************************
 * Description   : real fft process for fftsizes 256 and 768, when input is
 *                 real,
 *                     ...
 *                     ...
 *
 * Parameters    : x          - Data input and output.
 *                 scratchpad - intermediate array used for efficient
 *                              calculation of fft
 *                 n          - Number of complex indata (n =< 256)
 *
 * Returns       : float representing N/2 value of fft
 ****************************************************************************/

float realfft_process(float *x, double *scratchpad, short n)
{
  assert(n==768);
  float last = 0.0f;

  if (n == 768)
  {
    last = gt_real(768,
      (GT_ACCESSTYPE *) x,
      (GT_ACCESSTYPE *) scratchpad,
      (GT_ACCESSTYPE *) tw768_real,
      gt384,
      0);
  }

  return last;
}


/*****************************************************************************
 * Description   : real ifft process for fftsizes 256 and 768, when input is
 *                 complex and result is real
 *                     ...
 *                     ...
 *
 * Parameters    : x          - Data input and output.
 *                 scratchpad - intermediate array used for efficient
 *                              calculation of ifft

 *                 lastval    - last N/2 value from realfft_process (use 0.0f if
 *                              not known)
 *                 n          - Number of complex indata (n == 768 || n== 256)
 *
 * Returns       : void
 ****************************************************************************/
void realIfft_process(float *x, double *scratchpad, float lastval, short n)
{
  assert(n==768);
        gt_real_inv_noDivide(n,
                             (GT_ACCESSTYPE *) x,
                             (GT_ACCESSTYPE *) scratchpad,
                             (GT_ACCESSTYPE *) tw768_real,
                             gt384,
                             0,
                             lastval);
}
