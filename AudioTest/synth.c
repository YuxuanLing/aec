/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Geir Ole �verby, Tandberg Telecom AS
 * Co-author     : -
 *
 * Switches      : -
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
 **************************************************************************/

#include "synth_priv.h"
#include "synth.h"
#include "sFirData16.h"
#include "sFirData48.h"
#include "audec_defs.h" /* contains parameters and struct definitions */
#include "fft.h"
#include "adjustments.h" /* for postfilter */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef _TMS320C6X
#define FHI(a) _itof(_hi(a))
#define FLO(a) _itof(_lo(a))
#else
#define FHI(a) *((float*)(&a) + 1)
#define FLO(a) *(((float*)(&a)))
#endif

#define NUM_SYNTH_FILTERS_48 2
#define DEFAULT_SYNTH_FILTER_48_IDX 1

static const float
sfir48[NUM_SYNTH_FILTERS_48][3360] = {{SYNTH48_FILTER_SYMM_TAPS},
		                              {SYNTH48_FILTER_ASYMD1536_TAPS}};



static void synth_copy(SYNTH_PTR pSynth, float * fftst, float * fft);
static void synth_build (SYNTH_PTR pSynth, float * fftst, const float * sfir, float * out);

/******************************************************************
 * SYNTH_CREATE
 *  Description: Allocate memory for the synthese filter.
 *
 *  Parameters:  -
 *
 *  Returns:     Pointer to SYNTH structure used by the synthese filter
 *
 *  Note:        Memory for the SYNTH structure is allocated here. To free
 *               this memory, synth_destroy has to be called.
 *****************************************************************/

SYNTH_PTR synth_create()
{
    SYNTH_PTR pSynth;

    pSynth = (SYNTH_PTR) malloc(sizeof(SYNTH));
    if(pSynth == NULL)
    {
        fprintf(stderr, "synth_create: Could not allocate synthese buffer \n");
        return 0;
    }

    pSynth->fftst = (float *) malloc(sizeof(float) * FFTSIZE * ((SFIRLEN48 + FRAMESIZE - 1) / FRAMESIZE));  //pSynth->synthDlLength
    if(pSynth->fftst == NULL)
    {
        fprintf(stderr, "synth_create: Could not allocate synthese delay line \n");
        free(pSynth);
        return 0;
    }

    return pSynth;
}

/******************************************************************
 * SYNTH_DESTROY
 *  Description:  Frees memory for fftst matrix and synth struct
 *
 *  Parameters:   pSynth - Pointer to the SYNTH structure that shall be freed
 *
 *  Restrictions: synth_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void synth_destroy(SYNTH_PTR pSynth)
{
    free(pSynth->fftst);
    free(pSynth);
}


/******************************************************************
 * SYNTH_INIT
 *  Description:    Initialization of the struct that contains the necessary
 *                  data for each signal to be analysed
 *
 *  Parameters:     pSynth - Pointer to the SYNTH structure to be initialized
 *
 *  Returns:        None
 *
 *  Restrictions:   synth_create has to be called once before this function
 *                  is called.
 *****************************************************************/
void synth_init(SYNTH_PTR pSynth)
{
    int i;

    /* fftst_ix and fftst_iy are x and y coordinates
       of where the current subband signal recieved has its start,
       ie: fftst[fftst_ix+fftst_iy*FFTSTLX] contains subband 0 in
       the last input to the synthese filter*/

    pSynth->fftst_ix = 0;
    pSynth->fftst_iy = 0;

    pSynth->postfiltState = 0.0f;
    pSynth->fftStlY = FFTSIZE;
    pSynth->fftStlX = ((SFIRLEN48 + FRAMESIZE - 1) / FRAMESIZE);
    pSynth->synthDlLength = pSynth->fftStlX * pSynth->fftStlY;
    pSynth->synthBuildLen = FRAMESIZE * pSynth->fftStlX;
    pSynth->sfirSDRAM = sfir48[DEFAULT_SYNTH_FILTER_48_IDX];
    pSynth->filtertype = DEFAULT_SYNTH_FILTER_48_IDX;
    pSynth->postfiltercoeff = PFT48_COEFF2;
    pSynth->postfiltergain = PST48_GAIN;

    for(i = 0; i < pSynth->synthDlLength; i++)
    {
        pSynth->fftst[i] = 0.0f;
    }

}

/***************************************************************************
 * SYNTH SET FILTERBANK PROTOTYPE FILTER
 *  Description: select synthesis filter
 *
 *   Parameters: pSynth - Pointer to the SYNTH structure
 *               filtertype - integer
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: synth_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void synth_setFbPrototypefilt(SYNTH_PTR pSynth, int filtertype)
{
    if((filtertype >=0) && (filtertype < NUM_SYNTH_FILTERS_48))
    {
        pSynth->sfirSDRAM = sfir48[filtertype];
        pSynth->filtertype = filtertype;
        printf("synth_setFbPrototypefilt: filtertype is set to %d \n",filtertype);
    }
    else
    {
        fprintf(stderr, "synth_setFbPrototypefilt: prototype filter index "
                        "%d is invalid\n",filtertype);
    }
}

void synth_status(SYNTH_PTR pSynth)
{
    printf("\rStatus: synthesis");
    printf("\r\n   fftst          - %p", pSynth->fftst);
    printf("\r\n   filtertype     - %d", pSynth->filtertype);
    printf("\r\n   synthBuildLen  - %d", pSynth->synthBuildLen);
    printf("\r\n");
}

/******************************************************************************
 * SYNTH_COPY
 *  Description:
 *       The follwoing C code is used to copy the FFTSTLY real frequency bands
 *       into fftst. fftst holds the last FFTSTLX sets of FFTSTLY real frequency
 *       bands. These are stored in such a way that fftst_ix decrease by one for
 *       each set of subbands and fftst_iy icrease by FRAMESIZE. No unnessecary
 *       storage(almost) is used, and no movement of data in fftst is done. This
 *       is achived by wrap-around occuring in x and y direction at position
 *       FFTSTLX and FFTSTLY. The seemingly complicated storage sequence in
 *       fftst is done to synthezie as fast as possible in synth_build.
 *       (Earlier: Note that since the data was extended in the natural way
 *       before performing the ifft) We know that the complex frequecys should
 *       all be real when they come as input to this routine. Hence we only save
 *       every other (the real part) in fftst. The other half will be almost 0
 *       (restricted to precission in fft).
 *
 *  Parameters:  pSynth - Pointer to the SYNTH structure to be initialized
 *               fft    - fft buffer
 *
 *  Returns:     -
 ******************************************************************************/
static void synth_copy(SYNTH_PTR  pSynth,
                       float *fftst,
                       float *  fft)
{
    float *p;
    short iy, py;

    pSynth->fftst_ix--;
    if (pSynth->fftst_ix < 0)
    {
        pSynth->fftst_ix = (short)(pSynth->fftStlX - 1);
    }

    pSynth->fftst_iy +=  FRAMESIZE;
    if (pSynth->fftst_iy >=  pSynth->fftStlY)
    {
        pSynth->fftst_iy = (short)(pSynth->fftst_iy - pSynth->fftStlY);
    }

    p = fftst + pSynth->fftst_ix;
    py = pSynth->fftst_iy;


    for (iy = py; iy < pSynth->fftStlY; iy++)
    {
        p[iy * pSynth->fftStlX] = *fft++;
    }

    for(iy = 0; iy < py; iy++)
    {
        p[iy * pSynth->fftStlX] = *fft++;
    }
}



/******************************************************************************
 * SYNTH_BUILD Description:
 *  The follwoing C code is used to build the
 *  synthezied signal from the last FFTSTLY real frequency bands
 *  stored in fftst. Each sample out is computed from FFTSTLX
 *  different elements in fftst. Because of the nifty way the samples
 *  are stored in fftst these are right next to each other, and to
 *  compute the next sample out one only has to move one line
 *  (y-direction) down, ie FFTSTLX higher.  The code below may seem a
 *  bit hairy, but this is only because of the wrap-around in x and y
 *  direction that has to be accounted for. The samples out are
 *  computed starting with the last one (FRAMESIZE-1) and working its
 *  way up to the first sample to go out.
 *
 *  Parameters:  pSynth - Pointer to the SYNTH structure to be initialized
 *               *sfir  - Pointer to the synthesis filter
 *               fft    - fft buffer
 *
 *  Returns:     -
******************************************************************************/
static void synth_build (SYNTH_PTR pSynth,
                         float *fftst,
                         const float *sfir,
                         float *out)
{
    short ix, iy;
    int max, j, i, k, lastSmplOut, wrapXlim, jstart, done = 0;
    float *  p, s0, s1, s2, s3;
    const float *  ph;

    ph = sfir;

    ix = pSynth->fftst_ix;
    iy = pSynth->fftst_iy;

    lastSmplOut = (iy + ( FRAMESIZE - 1)) % ( pSynth->fftStlY);

    max = lastSmplOut - ( FRAMESIZE - 1) > 0 ? lastSmplOut - ( FRAMESIZE - 1): 0;

    wrapXlim =  pSynth->fftStlX - ix;

    p = fftst + lastSmplOut * pSynth->fftStlX + ix;
    ph += pSynth->synthBuildLen -  pSynth->fftStlX;

    k = ((FRAMESIZE));
    jstart = lastSmplOut;

    do {

        for(j = jstart; j >= max; j -= 4)
        {
            float *_p0 = &p[-0 * pSynth->fftStlX];
            float *_p1 = &p[-1 * pSynth->fftStlX];
            float *_p2 = &p[-2 * pSynth->fftStlX];
            float *_p3 = &p[-3 * pSynth->fftStlX];

            s0 = s1 = s2 = s3 = 0;

            for(i = 0; i < pSynth->fftStlX; i++)
            {
                //unsigned int wrapCorr = (i == wrapXlim) ? FFTSTLX : 0;
                int wrapCorr = (i >= wrapXlim) ? pSynth->fftStlX : 0;

                s0 += _p0[i - wrapCorr] * ph[i - 0 * pSynth->fftStlX];
                s1 += _p1[i - wrapCorr] * ph[i - 1 * pSynth->fftStlX];
                s2 += _p2[i - wrapCorr] * ph[i - 2 * pSynth->fftStlX];
                s3 += _p3[i - wrapCorr] * ph[i - 3 * pSynth->fftStlX];
            }

            out[--k] = s0; out[--k] = s1;
            out[--k] = s2; out[--k] = s3;

            p -= 4 * pSynth->fftStlX; ph -= 4 * pSynth->fftStlX;
        }

        /* speculatively set up parameter's for a possible second iteration */
        p += pSynth->synthDlLength;
        jstart = (pSynth->fftStlY - 1);
        max = iy;

        /*
         * if k != 0 then wraparound in y direction occurs in fftst
         * when computing the wanted FRAMESIZE samples. This means
         * that we need to do the second iteration of the do-while loop.
         * If k==0 then we are done after the first iteration
         */
        done += !k;

    } while (!done++);
}


/******************************************************************************
 * SYNTH_PROCESS
 *
 *  Description: Filterbank: Takes a signal from the frequency-domain and
 *       puts the corresponding signal in the time-domain, after
 *       some filtering.
 *
 *       The follwoing C code is used to synthezie the in signal
 *       containing a set of complex frequecy bands to the out
 *       signal. Only half of the subbands are represnted in the
 *       insignal, the rest are calculated based on the fact that
 *       the signal the in signal was analyzed from should be real.
 *       The filtering strongly depends on the fact that the length
 *       of the synthezie filter h, (not including the last part
 *       filled with only zeros) is a multiplum of the number of
 *       subbands, FFTSTLY. This makes it possible to reduce the
 *       number of filtertaps used in computing each sample out to
 *       (SFIRLEN/FRAMESIZE) (rounded up) by computing an ifft of
 *       the subbands first. Note also that this routine will
 *       overwrite the last(FFTSTLY) real subbands with 0, because
 *       we do not want these to contribute. h is padded with zeroes
 *       until it is a multiplum of FRAMESIZE to acheive a simpler
 *       synth_build.
 *
 *   Parameters: pSynth        - Pointer to the SYNTH structure
 *               inbuf         - Pointer to buffer in frequency-domain (subbands)
 *               outbuf        - Pointer to buffer in time-domain (sound samples)
 *               bDoPostfilter - Applies postfiltering if true.
 *
 *   Returns: None
 *
 *   Note: The buffersize is set in the synth_create operation. The
 *         output consists of FRAMESIZE real values.
 *
 *   Note (TMS): It is required (by the fft) that in-vector lies on a
 *               double-boundary
 *
 *   Globals used: -
 *
 *   Restrictions: This function is dependent on the fft-module, which has
 *                 to be initialized before this function is called.
 ***************************************************************************/
void synth_process(SYNTH_PTR pSynth, COMPLEX32 *in, float *out)
{
    __declspec(align(64)) double ifftScratchPad[FFTSIZE];

    realIfft_process((float *)in, ifftScratchPad, 0.0f, FFTSIZE);

    synth_copy(pSynth, pSynth->fftst, (float *)in);

    synth_build(pSynth, pSynth->fftst, pSynth->sfir, out);

    /* postfiltering */
    ecpostfilt(out,
               &(pSynth->postfiltState),
               pSynth->postfiltercoeff,
               pSynth->postfiltergain);
}


/* SYNTH_LOADFFTST ******************************************************
 * Auth.: JPS
 * Desc.: loads pSynth->fftst into iram
 ************************************************************************/
void synth_loadfftst(SYNTH_PTR pSynth)
{
    /* MEMXFERMUSTDIE */
    (void)pSynth;
}

/* SYNTH_LOADSFIR ********************************************************
 * Auth.: JPS
 * Desc.: loads sfir into iram
 ************************************************************************/
void synth_loadSfir(SYNTH_PTR pSynth)
{
    /* MEMXFERMUSTDIE */
    pSynth->sfir = pSynth->sfirSDRAM;
}

/* SYNTH_LOADSYNTH *******************************************************
 * Auth.: JPS
 * Desc.: loads synthese struct into iram
 * Returns: pointer to synthese struct in iram
 ************************************************************************/
SYNTH_PTR synth_loadSynth(SYNTH_PTR pSynthIn)
{
    /* MEMXFERMUSTDIE */
    return pSynthIn;
}

/* SYNTH_FLUSHSYNTH *******************************************************
 * Auth.: JPS
 * Desc.: flush synthese struct to sdram
 ************************************************************************/
void synth_flushSynth(SYNTH_PTR pSynthSdram, SYNTH_PTR pSynthIram)
{
    /* MEMXFERMUSTDIE */
    (void)pSynthSdram;
    (void)pSynthIram;
}


/******************************************************************************
 * SYNTH_PROCESSLOOPDETECT
 *****************************************************************************/
void synth_processLoopDetect(SYNTH_PTR pSynth,
                             float *fftst,
                             float *in,
                             float * out)
{
    synth_copy(pSynth, fftst, in);
    synth_build(pSynth, fftst, pSynth->sfir, out);
}


#ifdef UNITTEST
/******************************************************************************
 * UNITTEST
 *   Description:
 *       Run-trhough of synthese filter with 10 blocks of input signal.
 *       The test input signal to this unittest is the fftout (mic.micfft)
 *       of analysefilter in the matlab model when the input to analyse is
 *       random noise
 *       The result is verified with matlab.
 *****************************************************************************/
#include "testdata/synthtestdata.h"  //input and output data defined here
#include "unittest.h"
#include <math.h>

#undef WRITE_UNITTEST_OUTPUT /* define to write to file for further analysis */

/* Private testdata */
static float syntheseTestInput48[]  = {SYNTHESE48_TESTINPUT_DEFINE};
static float syntheseSymmTestOutput48[] = {SYNTHESE48_SYMM_TESTOUTPUT_DEFINE};
static float syntheseAsymD1536TestOutput48[] = {SYNTHESE48_ASYMD1536_TESTOUTPUT_DEFINE};

/******************************************************************************
* COMPARERESULT
* Compares two float buffers
 *****************************************************************************/
static int compareResults(float *inp1, float *inp2, int len, float epsilon)
{
    int differ = 0;
    int i;

    for(i = 0;i < len && differ != 1; i++)
    {
        if(fabs(inp1[i] - inp2[i]) > epsilon)
        {
            differ = 1;
        }
    }

    return differ;
}

/******************************************************************************
 * UNITTEST_SYNTHESE
 * Commented lines are for writing to file
 *****************************************************************************/
void unittest_synthese()
{
    SYNTH_PTR pSynth;
    float *inbuf;
    float *outbuf;   //pointer to output buffer
    int size[] = {7680, 7680};
    char const * testName[] = {"Synth-48kHz", "Synth-48kHz(32ms delay)"};
    float *syntheseTestInputs[] = {syntheseTestInput48, syntheseTestInput48};
    float *syntheseTestOutputs[] = {syntheseSymmTestOutput48,
                                    syntheseAsymD1536TestOutput48};
    int filterIndices[] = {0,1};
    int i, m, numFrames;
    int differ;  //criteria flag to be returned by function. 0: OK, 1: ERROR, 2:Nothing done, hence ERROR
    float epsilon = 1e-5;

#ifdef WRITE_UNITTEST_OUTPUT
    int count;
    FILE *fid;
    const char *Filename[] = {"/home/jtk/tmp/synthesetestoutput_16_host.dat",
    		                  "/home/jtk/tmp/synthesetestoutput_48_host.dat",
    		                  "/home/jtk/tmp/synthesetestoutput_asymm_48_host.dat"};
#endif
    float * output;

    pSynth = synth_create();
    if(pSynth == NULL)
    {
        unittest_assert(0);
        return;
    }
    fft_init();

    for (m = 0; m < 2; m++)
    {
        unittest_context(testName[m]);
        differ = 0;

        synth_init(pSynth);
        if(filterIndices[m]!=-1)
        	synth_setFbPrototypefilt(pSynth, filterIndices[m]);

        /* allocating and initializing to zero */
        inbuf  =(float *) calloc(FFTSIZE, sizeof(float));
        outbuf =(float *) calloc(FRAMESIZE, sizeof(float));

        numFrames = size[m] / FFTSIZE;

        output =(float *) calloc(numFrames*FRAMESIZE, sizeof(float));

        for(i = 0; i < numFrames; i++)
        {
            memcpy(inbuf,
                   syntheseTestInputs[m] + i * FFTSIZE,
                   sizeof(float) * FFTSIZE);
            scratchmem_lock();
            synth_loadfftst(pSynth);
            synth_loadSfir(pSynth);
            synth_process(pSynth, inbuf, outbuf);
            memxfer_waitAddress(outbuf);
            differ |= compareResults(outbuf,
                                    syntheseTestOutputs[m] + i * FRAMESIZE,
                                    FRAMESIZE, epsilon);
            scratchmem_unlock();
            memcpy(output+i*FRAMESIZE, outbuf, sizeof(float)*FRAMESIZE);
        }

#ifdef WRITE_UNITTEST_OUTPUT
        printf("write synthese output to file \n");
		fid = (void*)fopen(Filename[m],"wb");
		if (fid == NULL)
		{
			printf("Could not open file output\n");
		}
		else
		{
			count = fwrite(output,sizeof(float),numFrames*FRAMESIZE,fid);
			printf("Wrote %d samples to file %s\n", count, Filename[m]);
			fclose(fid);
		}
#endif
		free(output);
        free(inbuf);
        free(outbuf);

        unittest_assert(differ == 0);
    }

    synth_destroy(pSynth);

}
#endif /* UNITTEST */
