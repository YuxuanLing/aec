/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 * Revision      : $Revision: 1.12 $
 *
 * Author        : Ingvar Flaten Aarnes, Tandberg Telecom AS
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

#include "analyse.h"
#include "analyse_priv.h"
#include "aFirData48.h"
#include "fft.h"
#include "adjustments.h"
#include "mathfun.h"

#include "mm_malloc.h"

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

#define NUM_ANALYSIS_FILTERS_48 2
#define DEFAULT_ANALYSIS_FILTER_48_IDX 1

static const float
afir48[NUM_ANALYSIS_FILTERS_48][4608*2] = {{ANALYSIS_FILTER_SYMM_TAPS_48,
		                                  ANALYSIS_FILTER_SYMM_TAPS_48},
		                                 {ANALYSIS_FILTER_ASYMD1536_TAPS_48,
		                                  ANALYSIS_FILTER_ASYMD1536_TAPS_48}};

/* Overlap routine: */
static int analyse_overl(float *fftin,
                         float *inbuf,
                         float *afirtaps,
                         int startIndex);
/* Importing the new samples: */
static int analyse_copy(float *newbuf,
                        float *inbuf,
                        int inbufStart);

/***************************************************************************
 * ANALYZE CREATE
 *  Description: Claims memory for the analyze filter, returns pointer.
 *
 *   Parameters: -
 *
 *   Returns:   A pointer to a ANALYSE structure used by the analyse filter
 *
 *   Note:      Memory for the ANALYSE structure is allocated here. To free
 *              this memory, analyse_destroy has to be called.
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
ANALYSE_PTR analyse_create(FILTERBANK_USE_TYPE filterbankUse)
{
    ANALYSE_PTR pAnalyse;

    pAnalyse = (ANALYSE_PTR) malloc(sizeof(ANALYSE));
    if(pAnalyse == NULL)
    {
        fprintf(stderr, "analyse_create: Could not allocate analysis struct\n");
        return 0;
    }

    pAnalyse->inbufSDRAM  = (float *) malloc(sizeof(float) * AFIRLEN48);
    if(pAnalyse->inbufSDRAM == NULL)
    {
        fprintf(stderr, "analyse_create: Could not allocate analysis inbuf\n");
        free(pAnalyse);
        return 0;
    }

    if (filterbankUse == FILTERBANK_EC)
    {
        pAnalyse->fftout =(COMPLEX32 *) _mm_malloc(sizeof(COMPLEX32) * FFTSIZE/2, 64);
        if(pAnalyse->fftout == NULL)
        {
            fprintf(stderr, "analyse_create: Could not allocate analysis fftout in iram\n");
            free(pAnalyse->inbufSDRAM);
            free(pAnalyse);
            return 0;
        }
    }
    else
    {
        pAnalyse->fftout  = NULL; /* zerodelaynlp use scratchmem */
    }

    return pAnalyse;
}


/***************************************************************************
 * ANALYZE DESTROY
 *  Description: Frees memory. Should not be necessary.
 *
 *   Parameters: pAnalyse - Pointer to the ANALYSE structure that shall be
 *              freed
 *
 *      Returns: -
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: analyse_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void analyse_destroy(ANALYSE_PTR pAnalyse, FILTERBANK_USE_TYPE filterbankUse)
{
    if (filterbankUse == FILTERBANK_EC)  /* zerodelaynlp use scratchmem */
    {
        _mm_free(pAnalyse->fftout);
    }

    free(pAnalyse->inbufSDRAM);
    free(pAnalyse);
}

/***************************************************************************
 * ANALYZE INIT
 *  Description: Initialization of the struct that contains the necessary
 *               data for each signal to be analysed
 *
 *   Parameters: pAnalyse - Pointer to the ANALYSE structure to be initialized
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: analyse_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void analyse_init(ANALYSE_PTR pAnalyse, FILTERBANK_USE_TYPE filterbankUse)
{
    int i = 0;

    pAnalyse->afir = afir48[DEFAULT_ANALYSIS_FILTER_48_IDX];
    pAnalyse->filtertype = DEFAULT_ANALYSIS_FILTER_48_IDX;
    pAnalyse->prefiltercoeff = PFT48_COEFF2;
    pAnalyse->prefiltergain = PRE48_GAIN;

    pAnalyse->inbufStart = 0;
    pAnalyse->prefiltState = (float)0.0;
    for (i = 0; i < AFIRLEN48; i++)
    {
        pAnalyse->inbufSDRAM[i] = (float)0.0;
    }
    if (filterbankUse == FILTERBANK_EC) /* zdnlp uses scratchmem */
    {
        for (i = 0; i < FFTSIZE/2; i++)
        {
            pAnalyse->fftout[i].re = 0.0f;
            pAnalyse->fftout[i].im = 0.0f;
        }
    }
}

/***************************************************************************
 * ANALYZE SET FILTERBANK PROTOTYPE FILTER
 *  Description: select analysis filter
 *
 *   Parameters: pAnalyse - Pointer to the ANALYSE structure
 *               filtertype - integer
 *
 *      Returns: void
 *
 * Restrictions: analyse_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void analyse_setFbPrototypefilt(ANALYSE_PTR pAnalyse, int filtertype)
{
    if((filtertype >=0) && (filtertype < NUM_ANALYSIS_FILTERS_48))
    {
        pAnalyse->afir = afir48[filtertype];
        pAnalyse->filtertype = filtertype;
        printf("analyse_setFbPrototypefilt: filtertype is set to %d\n",filtertype);
    }
    else
    {
        fprintf(stderr, "analyse_setFbPrototypefilt: prototype filter index %d is invalid\n",filtertype);
    }
}

void analyse_status(ANALYSE_PTR pAnalyse)
{
    printf("\rStatus: analyse");
    printf("\r\n   inbufSDRAM     - %p", pAnalyse->inbufSDRAM);
    printf("\r\n   fftout         - %p", pAnalyse->fftout);
    printf("\r\n   filtertype     - %d", pAnalyse->filtertype);
    printf("\r\n   prefiltState   - %f", pAnalyse->prefiltState);
    printf("\r\n");
}

#ifndef __ARM_NEON__
/***************************************************************************
 * ANALYSE PROCESS
 *  Description: Filterbank: Takes a sampled signal from the time domain
 *               and puts out the corresponding signal in the frequency-
 *               domain, after some filtering (PRE- and DC-filters).
 *
 *   Parameters: pAnalyse  - Pointer to the ANALYSE structure
 *               inbuf     - Pointer to n*480 new samples
 *               afir      - Pointer to analyse filter
 *
 *   Returns:    0
 *
 *   Note:       The buffersize is set in the analyse_create operation. The
 *               filter output is stored in member fftout in ANALYSE struct.
 *               The output consists of n*384 complex values,
 *               dependent of the sampling frequency. In the Texas
 *               implementation, n=2 and the second part is the conjugate
 *               of the first part.
 *
 *   Remark:     FFTSIZE is 768. 1 complex value is made up by two parts,
 *               one real part and one imaginary part. The filter output
 *               does not include the conjugated part.
 *
 *   Globals used: -
 *
 *   Restrictions: This function is dependent on the fft-module, which has
 *                 to be initialized before this function is called.
 ***************************************************************************/
void analyse_process(ANALYSE_PTR pAnalyse, float *inbuf, float * afir)
{
    float tst;/*realfft_process returns value of the N/2 as well, but will not use*/
    __declspec(align(64)) double fftScratchPad[FFTSIZE];

    /* prefiltering */
    ecprefilt(inbuf,
              &(pAnalyse->prefiltState),
              pAnalyse->prefiltercoeff,
              pAnalyse->prefiltergain);

    /* copy new data */
    analyse_copy(inbuf, pAnalyse->inbuf, pAnalyse->inbufStart);

    /* update index: */
    pAnalyse->inbufStart = (pAnalyse->inbufStart + FRAMESIZE) % AFIRLEN48;

    /* -OVERLAP- notice that the afir input is not the first element of afir[], but afir[AFIRLEN-inbufStart] */
    analyse_overl((float *)pAnalyse->fftout,
                  pAnalyse->inbuf,
                  afir,
                  pAnalyse->inbufStart);

    memcpy(pAnalyse->inbufSDRAM, pAnalyse->inbuf, sizeof(float) * AFIRLEN48);

    /* -FFT- outvektor format: [re(0), im(0), re(1), im(1), ... re(FFTSIZE-1), im(FFTSIZE-1)] */
    /* re(FFTSIZE) can be returned from this function if we want. */
    tst = realfft_process((float *)pAnalyse->fftout, fftScratchPad, FFTSIZE);
    (void) tst;
}
#endif

/***************************************************************************
 * OVERLAP
 *    Auth.: IFA (after the matlab-model)
 *
 *    Description: Takes 768th element of inbuf and multiplies it
 *                 with the corresponding element of afir; adds these numbers
 *                 (p.t. six) and puts them into fftin at places
 *                 If 48kHz: fftsize - (startIndex mod(fftsize))+ 0, 2, 4 a.s.o.
 *                           GCD(AFIRLEN48, FRAMESIZE48) = 96, so p.t.
 *                           startIndex is a multiple of 96.
 *                 The input is not in "time-order", but the output has to be.
 *                 ("oldest" samples first into fft.)
 *
 *   Parameters:   fftin      - Pointer to first element of fftout
 *                 inbuf      - Pointer to first element of inbuf
 *                 afirtaps   - pointer to the element for the afir-filter that
 *                              corresponds to fftin[0]
 *                 startIndex - where the "oldest" soundsamples are stored in
 *                              the inbuf-array
 *
 *   Return:       0
 ***********************************************************************/
#ifdef _TMS320C6X
static int analyse_overl(float * const  fftin,
                         float * const  inbuf,
                         float *afirtaps,
                         int startIndex)
{
    int m, n, cntr;
    float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    double *a1, *b1 ;

    cntr = FFTSIZE - startIndex;
    cntr += AFIRLEN48;     /* make sure cntr is always positive before calculating modulo */
    cntr %= FFTSIZE;     /* modulo FFTSIZE */


    a1 = (double *) afirtaps;
    b1 = (double *) inbuf;


    for(n = 0; n < FFTSIZE; n += 8)
    {
        tmp0 = tmp1 = tmp2 = tmp3 = 0.0f;
        tmp4 = tmp5 = tmp6 = tmp7 = 0.0f;

        for(m = n / 2; m < AFIRLEN48 / 2; m += FFTSIZE / 2)
        {
            tmp0 += FLO(a1[m])     * FLO(b1[m]);
            tmp1 += FHI(a1[m])     * FHI(b1[m]);
            tmp2 += FLO(a1[m + 1]) * FLO(b1[m + 1]);
            tmp3 += FHI(a1[m + 1]) * FHI(b1[m + 1]);
            tmp4 += FLO(a1[m + 2]) * FLO(b1[m + 2]);
            tmp5 += FHI(a1[m + 2]) * FHI(b1[m + 2]);
            tmp6 += FLO(a1[m + 3]) * FLO(b1[m + 3]);
            tmp7 += FHI(a1[m + 3]) * FHI(b1[m + 3]);
        }

        fftin[cntr++] = tmp0; /* real part */
        fftin[cntr++] = tmp1; /* real part */
        fftin[cntr++] = tmp2; /* real part */
        fftin[cntr++] = tmp3; /* real part */
        fftin[cntr++] = tmp4; /* real part */
        fftin[cntr++] = tmp5; /* real part */
        fftin[cntr++] = tmp6; /* real part */
        fftin[cntr++] = tmp7; /* real part */

        cntr %= FFTSIZE;      /* if cntr has reached FFTSIZE, it is set to 0 */
    }

    return 0;
}  /* end (analyse_overl()) */

#else

static int analyse_overl(float *fftin,
                  float *inbuf,
                  float *afirtaps,
                  int startIndex)
{
    int rc = 0, m = 0, n = 0, cntr;

    cntr = FFTSIZE - startIndex;
    cntr += AFIRLEN48; /* make sure cntr is always positive before calculating modulo */
    cntr %= FFTSIZE; /* modulo FFTSIZE */

    for (n = 0; n < FFTSIZE; n++)
    {
        float tmp = 0.0;
        m = n;
        while (m < AFIRLEN48)
        {
            tmp += afirtaps[m] * inbuf[m];
            m   += FFTSIZE;
        }
        fftin[cntr++] = tmp;   /* input is real, no imaginary part nessecarry */
        cntr %= FFTSIZE;       /* if cntr has reached FFTSIZE, it is set to 0 */
    }
    return rc;
}  /* end (analyse_overl()) */
#endif

/***************************************************************************
 * ANALYSE_COPY
 * Auth.: IFA
 *
 * Desc.: Fills in new (newbuf-)data into inbuf; (overwrites old data)
 *        starting at inbufStart
 *
 *   Parameters:   newbuf      - Pointer to first element of fftout
 *                 inbuf       - Pointer to first element of inbuf
 *                 inbufstart  - pointer to the element for the afir-filter that
 *                               corresponds to fftin[0]
 *
 *   Returns:      inbufStart
 ***********************************************************************/
static int analyse_copy(float *newbuf, float *inbuf, int inbufStart)
{
    int i = 0;
    for (i = 0; i < FRAMESIZE; i++)
    {
        inbuf[inbufStart] = newbuf[i];
        inbufStart++;
        if (inbufStart == AFIRLEN48) {inbufStart=0;} /* modulo AFIRLEN */
    }
    return inbufStart;
}

/* ANALYSE_GAINADJUST
 * Auth.: JPS
 * Desc.: Multiply output from analyse filter with loudspeaker volume.
 *        Only in use on lsfft.
 */
void analyse_gainAdjust(ANALYSE_PTR pAnalyse, float loudsGain)
{
    int m;
    float loudsGainInv;
    COMPLEX32 *  fftout = pAnalyse->fftout;

    if (loudsGain > 1.0e-7f)
    {
        loudsGainInv = fastdiv(1.0f, loudsGain);     // loudsGainInv = 1/loudsGain
    }
    else
    {
        loudsGainInv = 1.0f;
        //printf("analyse_gainAdjust: loudsGain < EPSILON \n");
    }

    for(m = 0; m < FFTSIZE/2; m++)
    {
        fftout[m].re *= loudsGainInv;
        fftout[m].im *= loudsGainInv;
    }
}

void analyse_updateInbufStart(ANALYSE_PTR pAnalyse, int index)
{
    pAnalyse->inbufStart = index; 
}

int analyse_getInbufStart(ANALYSE_PTR pAnalyse)
{
    return(pAnalyse->inbufStart);
}



/* analyse_loadAnalyse ***************************************************
 * Auth.: JPS
 * Desc.: loads pAnalyse struct into iram
 * Returns: ptr to analyse struct in iram
 ************************************************************************/
ANALYSE_PTR analyse_loadAnalyse(ANALYSE_PTR pAnalyseIn)
{
    /* MEMXFERMUSTDIE */
    return pAnalyseIn;
}

/* analyse_flushAnalyse **************************************************
 * Auth.: JPS
 * Desc.: flush pAnalyse struct to sdram
 ************************************************************************/
void analyse_flushAnalyse(ANALYSE_PTR pAnalyseSdram, ANALYSE_PTR pAnalyseIram)
{
    /* MEMXFERMUSTDIE */
    (void)pAnalyseSdram;
    (void)pAnalyseIram;
}


/* ANALYSE_LOADINBUF *****************************************************
 * Auth.: JPS
 * Desc.: loads pAnalyse->inbuf into iram
 ************************************************************************/
void analyse_loadInbuf(ANALYSE_PTR pAnalyse)
{
    /* MEMXFERMUSTDIE */
    pAnalyse->inbuf = pAnalyse->inbufSDRAM;
}

/* ANALYSE_LOADAFIR ******************************************************
 * Auth.: JPS
 * Desc.: loads afir into iram
 ************************************************************************/
float * analyse_loadAfir(ANALYSE_PTR pAnalyse)
{
    /* update startindex for afir for memxfer */
    int bufStart = (pAnalyse->inbufStart + FRAMESIZE) % AFIRLEN48;

    return (float *)(pAnalyse->afir + AFIRLEN48 - bufStart);
}



#ifdef UNITTEST
/************************************************************************
 * UNITTEST
 * Desc.: run-trhough of analyse filter for both 16kHz and 48kHz
 *        with 0.3sec random input.
 *        Check against testvector verified in matlab
 *        Pre-filter is active
 *        All sub-bands to FFTSIZE/2
 ************************************************************************/
#include "unittest.h"
#include "testdata/analysetestdata.h"      /* input and output data defined here */
#include <math.h>

#undef WRITE_UNITTEST_OUTPUT /* define to write to file for further analysis */

/* Private testdata */
static float analyseTestInput[]    = {ANALYSE_TESTINPUT_DEFINE};

static float analyseSymmTestOutput48[] = {ANALYSE48_SYMM_TESTOUTPUT_DEFINE};
static float analyseAsymD1536TestOutput48[] = {ANALYSE48_ASYMD1536_TESTOUTPUT_DEFINE};

/*************************************************************************
 * COMPARERESULT
 *      Desc.: Compares two float buffers
 *************************************************************************/
static int compareResults(float *inp1, float *inp2, int len)
{
    int differ = 0;
    int i;
    float epsilon = 1.0e-7f;

    for(i = 0; i < len && differ != 1; i++)
    {
        if((fabsf(inp1[i] - inp2[i])) > epsilon)
        {
            differ = 1;
        }
    }
    return differ;
}
/*************************************************************************
 * UNITTEST_ANALYSE
 *
 *    Note: Commented lines are for writing to file
 ************************************************************************/
void unittest_analyse()
{
    ANALYSE_PTR pAnalyse; //pointer to analyse struct
    int   i, m;
    int   numFrames;
    int   differ; /* criteria flag to be returned by function. 0: OK,
                     1: ERROR, 2:Nothing done, hence ERROR */
    int   size[] = {4800, 4800};
    const char  * testName[] = {"Analysefilter-48kHz",
    		                    "Analysefilter-48kHz (32ms delay)"};
    float * analyseTestOutputs[] = {analyseSymmTestOutput48,
    		                        analyseAsymD1536TestOutput48};
	int analyseFilterIndices[] = {0, 1};/* -1 means do not change from default.. */
    float * inbuf, *pAfir;
    int   createOk = 1;
    FILTERBANK_USE_TYPE filterbankUse = FILTERBANK_ZDNLP;

#ifdef WRITE_UNITTEST_OUTPUT
    int count;
    FILE *fid;
    const char *Filename[] = {"/home/jtk/tmp/analysetestoutput_48_host.dat",
    		                  "/home/jtk/tmp/analysetestoutput_asymm_48_host.dat"};
    //char *Filename ="analysetestoutput_host.dat";
#endif
    float *outbuf;


    /* Allocate memory
     * create analyse buffers for unittest instead of calling analyse_create which allocates in IRAM */

    pAnalyse = analyse_create(filterbankUse);
    if ( (pAnalyse != NULL) && (pAnalyse->fftout == NULL) )
    {
        pAnalyse->fftout =(float *) malloc(sizeof(float) * FFTSIZE);
        for (i = 0; i < FFTSIZE; i++)
        {
            pAnalyse->fftout[i] = 0.0f;
        }
    }

    if(pAnalyse == NULL)  //try to allocate only in sdram
    {
        createOk = 0;

        pAnalyse = (ANALYSE_PTR) malloc(sizeof(ANALYSE));
        if(pAnalyse == NULL)
        {
            fprintf(stderr, "analyse unittest: Could not allocate analysis struct");
            unittest_assert(0);
            return;
        }

        pAnalyse->inbufSDRAM  = (float *) malloc(sizeof(float) * AFIRLEN48);
        if(pAnalyse->inbufSDRAM == NULL)
        {
            fprintf(stderr, "analyse unittest: Could not allocate analysis inbuf");
            unittest_assert(0);
            free(pAnalyse);
            return;
        }

        pAnalyse->fftout =(float *) malloc(sizeof(float) * FFTSIZE);
        if(pAnalyse->fftout == NULL)
        {
            fprintf(stderr, "analyse unittest: Could not allocate analysis fftout");
            unittest_assert(0);
            free(pAnalyse->inbufSDRAM);
            free(pAnalyse);
            return;
        }
    }

    fft_init();

    for (m = 0; m < 2; m++)
    {
        /* Initialization*/
        unittest_context(testName[m]);
        differ = 0;

        //inbuf = dsp_memAllocSeg(IRAM, sizeof(float) * size[m], 8);
        inbuf = malloc(sizeof(float) * size[m]);
        //copy input because prefilter modifies input
        memcpy(inbuf, analyseTestInput, size[m] * sizeof(*inbuf));

        analyse_init(pAnalyse, filterbankUse);
        if(analyseFilterIndices[m]>=0)  /* -1 means do not change from default.. */
        	analyse_setFbPrototypefilt(pAnalyse, analyseFilterIndices[m]);

        numFrames = size[m] / FRAMESIZE;

        outbuf =(float *) calloc(numFrames*FFTSIZE,sizeof(float));

        /* Run process and verify result*/
        for(i = 0; i < numFrames; i++)
        {
            scratchmem_lock(); /* level 0 */

            pAfir = analyse_loadAfir(pAnalyse);
            analyse_loadInbuf(pAnalyse);

            analyse_process(pAnalyse, inbuf + i * FRAMESIZE, pAfir);  //runs analysefilter on mic buffer
            memcpy(outbuf+i*FFTSIZE, pAnalyse->fftout, sizeof(float)*FFTSIZE);
            differ |= compareResults(pAnalyse->fftout, analyseTestOutputs[m] + i * FFTSIZE, FFTSIZE);

            scratchmem_unlock(); /* level 0 */
         }


#ifdef WRITE_UNITTEST_OUTPUT
        printf("write analyse fftout to file \n");
        fid = (void*)fopen(Filename[m],"wb");
        if (fid == NULL)
        {
            printf("Could not open file output");
        }
        else
        {
	        count = fwrite(outbuf,sizeof(float),numFrames*FFTSIZE,fid);
			printf("Wrote %d samples to file %s\n", count, Filename[m]);
	        fclose(fid);
        }
#endif
        free(outbuf);

        unittest_assert(differ == 0);
        //dsp_memFreeSeg(IRAM, inbuf, sizeof(float) * size[m]);
        free(inbuf);
   }

   if(createOk)
   {
       analyse_destroy(pAnalyse, filterbankUse);
   }
   else
   {
       free(pAnalyse->fftout);
       free(pAnalyse->inbufSDRAM);
       free(pAnalyse);
   }
}

#endif /* UNITTEST */
