/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Espen Holmbakken (EHO)
 *                 Geir Ole �verby  (GEO)
 *
 * Co-author     : -
 *
 * Switches      : <COMPILER SWITCH HERE>
 *           <add description here>
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
 *
 * yyyy-mm-dd    : <signature>
 *         File created.
 *
 **************************************************************************/

#include "echocomp.h"
#include "echocomp_priv.h"
#include "lsprocess.h"
#include "mathfun.h"
#include "apa.h"
//#include "nlms.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lsprocess_priv.h"

#include "mm_malloc.h"

/* forward declarations of private functions */
void echocomp_shift(ECHOCOMP_PTR pEchocomp,
                    int ch, // channel
                    const COMPLEX32 *  lsBuf,
                    const int subband_start,
                    const int subband_end);

float echocomp_eval(float * level,
                    const COMPLEX32 *  buf,
                    const int subband_start,
                    const int subband_end);

void echocomp_weight(COMPLEX32 *  mixedbuf,
                     const COMPLEX32 *  fxbuf,
                     const COMPLEX32 *  adbuf,
                     const float fxweight,
                     const int subband_start,
                     const int subband_end);

void echocomp_find_decay(ECHOCOMP_PTR pEchocomp);

void echocomp_subtract(COMPLEX32 *  canc,
                       const COMPLEX32 *  mic,
                       const COMPLEX32 *  est,
                       const int subband_start,
                       const int subband_end);

float powf_to_n(float x, int n);

/***************************************************************************
 * ECHOCOMP_CREATE
 *  Description: Claims memory for the echo compensation filters, returns pointer.
 *
 *   Parameters: -
 *
 *   Returns:   A pointer to an ECHOCOMP structure
 *
 *   Note:      Memory for the ECHOCOMP structure is allocated here. To free
 *              this memory, echocomp_destroy has to be called.
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
ECHOCOMP_PTR echocomp_create(int channels)
{
    static const short filtlen48[SUBUSED] = { FILTLEN_TAPS48 };
    ECHOCOMP_PTR pEchocomp = (ECHOCOMP_PTR) malloc(sizeof(ECHOCOMP));
    int i;

    if ( pEchocomp == NULL )
    {
        fprintf(stderr, "echocomp_create: Could not allocate echocomp struct\n");
        return NULL;
    }
    memset(pEchocomp, 0, sizeof(ECHOCOMP)); /* set to zero to facilitate cleanup */

    pEchocomp->nChannels = channels;

    for (i = 0; i < channels; i++)
    {
        ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[i];
        int j;

        for (j = 0; j < SUBUSED; j++)
        {
            ECHOCOMP_SUBBAND *subband = &channel->subband[j];
            int len = filtlen48[j];

            subband->filtlen = len;
            len = ((len + 3) & ~3) + 4;  // filtlen, rounded up to multiple of 4, plus 4

#if defined(HAVE_IOS)
            subband->delayline_base = (COMPLEX32 *) _mm_malloc((len+8) * sizeof(COMPLEX32), 64);
#else
            subband->delayline = (COMPLEX32 *) _mm_malloc(len * sizeof(COMPLEX32), 64);
#endif
            subband->adapt     = (COMPLEX32 *) _mm_malloc(len * sizeof(COMPLEX32), 64);
            subband->fixed     = (COMPLEX32 *) _mm_malloc(len * sizeof(COMPLEX32), 64);

#if defined(HAVE_IOS)
            subband->delayline = subband->delayline_base + 8;
#endif
        
            if ( subband->delayline == NULL )
            {
                fprintf(stderr, "echocomp_create: Could not allocate echocomp delayline buffer\n");
                goto cleanup;
            }
            if ( subband->adapt == NULL )
            {
                fprintf(stderr, "echocomp_create: Could not allocate echocomp adapt buffer\n");
                goto cleanup;
            }
            if ( subband->fixed == NULL )
            {
                fprintf(stderr, "echocomp_create: Could not allocate echocomp fixed buffer\n");
                goto cleanup;
            }
        }
    }
    return pEchocomp;

 cleanup:
    echocomp_destroy(pEchocomp);
    return NULL;
}


/***************************************************************************
 * ECHOCOMP_DESTROY
 *  Description: Frees memory. Should not be necessary.
 *
 *   Parameters: pEchocomp - Pointer to the ECHOCOMP structure that shall be
 *              freed
 *
 *      Returns: -
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: echocomp_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void echocomp_destroy(ECHOCOMP_PTR pEchocomp )
{
    int i;

    for (i = 0; i < pEchocomp->nChannels; i++)
    {
        ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[i];
        int j;

        for (j = 0; j < SUBUSED; j++)
        {
            ECHOCOMP_SUBBAND *subband = &channel->subband[j];

#if defined(HAVE_IOS)
            if ( subband->delayline_base != NULL ) _mm_free(subband->delayline_base);
#else
            if ( subband->delayline != NULL ) _mm_free(subband->delayline);
#endif
            if ( subband->adapt != NULL )     _mm_free(subband->adapt);
            if ( subband->fixed != NULL )     _mm_free(subband->fixed);
        }
    }
    free(pEchocomp);
}

/***************************************************************************
 * ECHOCOMP_INIT
 *  Description: Initialization of the struct
 *
 *   Parameters: pEchocomp - Pointer to the ECHOCOMP structure to be initialized
 *               sampFreqKhz - integer, 16 or 48 kHz
 *
 *      Returns: void
 *
 *         Note: -
 *
 * Globals used: -
 *
 * Restrictions: echocomp_create has to be called once before this function
 *               is called.
 ***************************************************************************/
void echocomp_init(ECHOCOMP_PTR pEchocomp)
{
    static const float attackTapsInit[NUMBER_OF_ATTACK_TAPS] = { DECAY_ATTACK_TAPS_DEFAULT };
    static const float decay_default[SUBUSED] = { DECAY_DEFAULT };
    const float delta = 1.0e-5f;
    int i;

    pEchocomp->decAttData.decay_early_power = 0.0f;
    pEchocomp->decAttData.decay_late_power = 0.0f;
    pEchocomp->decAttData.decay_ix = 0;
    pEchocomp->decAttData.decay_compute = 0;
    pEchocomp->decAttData.subused_finddecay = SUBUSED_FINDDECAY48;

    for (i = 0; i < pEchocomp->nChannels; i++)
    {
        ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[i];
        float sum = 0.0f;
        int j;

        for (j = 0; j < NUMBER_OF_ATTACK_TAPS; j++)
        {
            channel->attack_taps_old[j] = attackTapsInit[j];
            channel->attack_taps_new[j] = attackTapsInit[j];
        }

        channel->number_of_used_attack_taps_new = DECAY_NUMBER_OF_USED_ATTACK_TAPS_DEFAULT;
        channel->number_of_used_attack_taps_old = DECAY_NUMBER_OF_USED_ATTACK_TAPS_DEFAULT;

        for (j = 0; j < channel->number_of_used_attack_taps_new; j++)
        {
            sum += channel->attack_taps_new[j];
        }
        channel->decay_attack_taps_sum_old = sum;
        channel->decay_attack_taps_sum_new = sum;

        for (j = 0; j < FILTLEN_MAX; j++)
        {
            channel->filter_power[j] = 0.0f;
        }

        for (j = 0; j < SUBUSED; j++)
        {
            ECHOCOMP_SUBBAND *subband = &channel->subband[j];
            float num = (float) subband->filtlen;
            float denom = channel->decay_attack_taps_sum_new;
            int number_of_attack_taps = channel->number_of_used_attack_taps_new;
            float this_weight = channel->attack_taps_new[number_of_attack_taps-1];
            int k;

            subband->delta = delta;
            subband->r11 = delta;
            subband->r22 = 0.0f;

            subband->FiltCorr.re = subband->FiltCorr.im = 0.0f;
            subband->beta_1.re = subband->beta_1.im = 0.0f;
            subband->beta_fixed.re = subband->beta_fixed.im = 0.0f;
            subband->r12.re = subband->r12.im = 0.0f;
            subband->est_err.re = subband->est_err.im = 0.0f;

            subband->decay = subband->decay_default = decay_default[j];

            for (k = number_of_attack_taps; k < subband->filtlen; k++)
            {
                this_weight *= subband->decay;
                denom += this_weight;
            }

            subband->decay_attack_gain = num / denom;
            subband->weights_need_recalculate = 1;

            for (k = 0; k < ((subband->filtlen + 3) & ~3) + 4; k++)
            {
                subband->delayline[k].re = subband->delayline[k].im = 0.0f;
                subband->adapt[k].re = subband->adapt[k].im = 0.0f;
                subband->fixed[k].re = subband->fixed[k].im = 0.0f;
            }
        }
    }

    for (i = 0; i < NUM_SUBBAND_GROUPS_STEREO; i++)
    {
        pEchocomp->level[i] = 0.0f;
        pEchocomp->fixedlevel[i] = 0.0f;
        pEchocomp->adaptlevel[i] = 0.0f;
        pEchocomp->copyd[i] = 0;
        pEchocomp->fixedweight[i] = 0.0f;
    }

    pEchocomp->adaptMaxTap = 0.0f;
    pEchocomp->fixedMaxTap = 0.0f;
    pEchocomp->aerlInverse = 0.0f;

    /* APA variables */
    pEchocomp->mu = 1.0f;

    /* debug parameters */
    pEchocomp->calcNewDeltaOn = true;
    pEchocomp->minDelta = 1e-5f;
    pEchocomp->maxDelta = 5e-4f;
    pEchocomp->debug = 0;
    pEchocomp->filter = 0;

    pEchocomp->nlmsDelta = 0.0001f;

//==========================================================
// Code to detect CPU capabilities on IA processors
//==========================================================
#ifdef __INTEL_COMPILER
    pEchocomp->MMX_PRESENT   = 0;
	pEchocomp->SSE_PRESENT   = 0;
	pEchocomp->SSE2_PRESENT  = 0;
	pEchocomp->SSE3_PRESENT  = 0;
	pEchocomp->SSSE3_PRESENT = 0;
	pEchocomp->SSE4_PRESENT  = 0;


	{
#define MMX_FLAG 0x0800000
#define SSE_FLAG 0x2000000   // Safe to use -QxK/-xK
#define SSE2_FLAG 0x4000000  // Safe to use -QxN/-xN
#define SSE3_FLAG 0x0000001  // Safe to use -QxP/-xP
#define SSSE3_FLAG 0x0000200 // Safe to use -QxT/-xT
		unsigned int cpeinfo;
		unsigned int cpsse3;
		unsigned int vender[4] = {0,0,0,0};
		__asm {
			mov eax, 0
				cpuid
				mov long ptr [vender+0], ebx
				mov long ptr [vender+4], edx
				mov long ptr [vender+8], ecx
		}
		__asm {
			mov eax, 01h
				cpuid
				mov cpeinfo, edx
				mov cpsse3, ecx
		}
    if (strncmp((char *)vender,"GenuineIntel",12) != 0 &&
        strncmp((char *)vender,"AuthenticAMD",12) != 0)
    {
        cpsse3 = 0;
    }
		pEchocomp->MMX_PRESENT   = (cpeinfo & MMX_FLAG   )>0;
		pEchocomp->SSE_PRESENT   = (cpeinfo & SSE_FLAG   )>0;
		pEchocomp->SSE2_PRESENT  = (cpeinfo & SSE2_FLAG  )>0;
		pEchocomp->SSE3_PRESENT  = (cpsse3  & SSE3_FLAG  )>0;
		pEchocomp->SSSE3_PRESENT = (cpsse3  & SSSE3_FLAG )>0;
	}
#endif
//==========================================================
}


/***************************************************************************
 * ECHOCOMP PROCESS
 *  Description: Main routine for echocpmpansation processing.
 *
 *   Parameters: pEchocomp  - Pointer to the ECHOCOMP structure
 *               lsBuf      - Pointer to the loudspeaker fftout, must be double-word aligned
 *               micin      - Pointer to fftsize new samples
 *               echoest    - Pointer to estimated echo, SUBUSED48*2 samples
 *
 *   Returns:    0
 *
 *   Note:       Processes one channel (mono only)
 *
 *   Globals used: -
 *
 *   Restrictions: - lsBuf must be double-word aligned
 *                 - This function is dependent on apa.c
 ***************************************************************************/
int echocomp_process(ECHOCOMP_PTR pEchocomp,
                     const COMPLEX32 *lsBuf,
                     COMPLEX32 *echoest,
                     const COMPLEX32 *micin,
                     const float excgain,
                     bool runShift)
{
  COMPLEX32 adest[SUBUSED];
  COMPLEX32 fxest[SUBUSED];
  COMPLEX32 adcanc[SUBUSED];
  COMPLEX32 fxcanc[SUBUSED];
  int i, k;
  int subband_start = 0;
  int subband_end = 0;
  int subband_start_next = 0;
  int subband_end_next = 0;
  short decay_computed = 0;

  bool copied = false;

  int subband_groups[NUM_SUBBAND_GROUPS] = {27, 28, 31, 34 ,39 ,45, 57, 59};  /* dynamic groupsize for memxfer optimization */

  pEchocomp->adaptMaxTap = 0.0f;
  subband_end = subband_groups[0];

  for ( k = 0; k < NUM_SUBBAND_GROUPS; k++ )
  {
    float adaptMaxTap_subband = 0.0f;

    if ( k != NUM_SUBBAND_GROUPS - 1 )
    {
      subband_start_next = subband_end;
      if ( k == NUM_SUBBAND_GROUPS - 2 )
        subband_end_next   = SUBUSED; //needed for 16kHz
      else
        subband_end_next   = subband_start_next + subband_groups[k+1];
    }

    if ( decay_computed == 0 )
    {
      if ( ((subband_start <= pEchocomp->decAttData.decay_ix) && (pEchocomp->decAttData.decay_ix < subband_end)) || 
          (pEchocomp->decAttData.decay_ix == -1) )
      {
        echocomp_find_decay(pEchocomp);
        decay_computed = 1;
      }
    }

    if ( pEchocomp->calcNewDeltaOn )
    {
      echocomp_calcNewDelta(pEchocomp,
        excgain,
        subband_start,
        subband_end);
    }

    if ( runShift )
    {
      echocomp_shift(pEchocomp, 0,
        lsBuf,
        subband_start,
        subband_end);
    }

    apa_calcR11_and_R12(pEchocomp,
      subband_start,
      subband_end);

    apa_calcFiltCorr(pEchocomp,
      subband_start,
      subband_end);

    adaptMaxTap_subband = apa_filter_APA2(pEchocomp,
      adest,
      fxest,
      subband_start,
      subband_end);

    echocomp_subtract(fxcanc,
      micin,
      fxest,
      subband_start,
      subband_end);

    echocomp_subtract(adcanc,
      micin,
      adest,
      subband_start,
      subband_end);

    echocomp_eval(&pEchocomp->fixedlevel[k],
      fxcanc,
      subband_start,
      subband_end);

    echocomp_eval(&pEchocomp->adaptlevel[k],
      adcanc,
      subband_start,
      subband_end);

    echocomp_eval(&pEchocomp->level[k],
      micin,
      subband_start,
      subband_end);

    pEchocomp->copyd[k] = apa_findbest_APA(pEchocomp,
      k,
      subband_start,
      subband_end);

    if ( pEchocomp->copyd[k] == -3 )    /* Adapt set to Zero */
    {
      for (i = subband_start; i < subband_end; i++)
      {
        adest[i].re = 0.0f;
        adest[i].im = 0.0f;
        adcanc[i] = micin[i];
      }
    }
    else if ( pEchocomp->copyd[k] == -1 )   /* Fixed copied to adapt */
    {
      for (i = subband_start; i < subband_end; i++)
      {
        adest[i] = fxest[i];
        adcanc[i] = fxcanc[i];
      }
    }
    else if ( pEchocomp->copyd[k] == 1 )    /* Adapt copied to fixed */
    {
      for (i = subband_start; i < subband_end; i++)
      {
        fxest[i] = adest[i];
        fxcanc[i] = adcanc[i];
      }
    }

    apa_adapt_APA(pEchocomp,
      adcanc,
      fxcanc,
      subband_start,
      subband_end);

    if ( pEchocomp->filter == 1 )
    {
      pEchocomp->fixedweight[k] = 1.0f;  /*listen to fixed*/
    }
    else if ( pEchocomp->filter == 2 )
    {
      pEchocomp->fixedweight[k] = 0.0f;  /*listen to adapt*/
    }

    echocomp_weight(echoest,
      fxest,
      adest,
      pEchocomp->fixedweight[k],
      subband_start,
      subband_end);

    subband_start = subband_start_next;
    subband_end = subband_end_next;

    if ( (pEchocomp->copyd[k] == 1) && ( adaptMaxTap_subband > pEchocomp->adaptMaxTap) )
    {
      pEchocomp->adaptMaxTap = adaptMaxTap_subband;
      copied = true;
    }
  }
  if ( copied )
  {
    pEchocomp->fixedMaxTap = sqrtf(pEchocomp->adaptMaxTap);
  }

  return 0;
}

/*
 * Calculating and returning the estimated aerlInverse.
 * It is basically the maxTap in the fixed filter, with a time smoothing.
 * The aerlInvers estimate is dependent on the miclevel setting and the loudspeaker volume.
  */
float echocomp_getAerl(ECHOCOMP_PTR pEchocomp, float miclevel, float loudsGain)
{
    float tmp;

    if (loudsGain < 1.0e-7f)
    {
            loudsGain = 1.0f;
    }

    tmp = pEchocomp->fixedMaxTap / (miclevel * loudsGain);
    if (tmp > 100.0f) {
        /* Safety valve in case fixedMaxTap is abnormally large, lslim should not attenuate more than 40dB */
        tmp = 100.0f - pEchocomp->aerlInverse;
    }
    else {
        tmp -= pEchocomp->aerlInverse;
    }

    if ( tmp>0 )
    {
        pEchocomp->aerlInverse += tmp * 0.08f;
    }
    else
    {
        pEchocomp->aerlInverse += tmp * 0.02f;
    }


    if ( pEchocomp->debug == 10 )
    /* Print AERL and adding this to the print from AUDLOCALIN_debug10print()*/
    {
        static int cnt = 0;
        static float Mean_AerlInverse = 0.0f, Max_AerlInverse = 0.0f;
        static float Min_AerlInverse = 1.0f;

        Mean_AerlInverse += pEchocomp->aerlInverse;
        Max_AerlInverse = max(Max_AerlInverse, pEchocomp->aerlInverse);
        Min_AerlInverse = min(Min_AerlInverse, pEchocomp->aerlInverse);

        if ( ++cnt == 20 )
        {
            printf("%11.5e, %11.5e, %11.5e %11.5e\n", Mean_AerlInverse*0.05, Min_AerlInverse, Max_AerlInverse, loudsGain);
            Mean_AerlInverse = 0.0f;
            Max_AerlInverse  = 0.0f;
            Min_AerlInverse  = 1.0f;
            cnt  = 0;
        }
    }
    return pEchocomp->aerlInverse;
}

/* Make decay estimate available */
//float* echocomp_loadDecay(ECHOCOMP_PTR pEchocomp)
//{
//    return pEchocomp->decAttData->decay;
//}

int echocomp_processStereo(ECHOCOMP_PTR pEcho,
                           const COMPLEX32 *lsBuf,
                           const COMPLEX32 *lsBuf2,
                           COMPLEX32 *echoest,
                           const COMPLEX32 *micin,
                           const float excgain,
                           bool runShift)

{
#if 0
    float tempdest[SUBUSED * 2 * 4];
    float *  adest;
    float *  fxest;
    float *  adcanc;
    float *  fxcanc;
    double *  adest_d;
    double *  fxest_d;
    double *  adcanc_d;
    double *  fxcanc_d;
    const double *  micin_d = (double *) micin;
    float * delaylines[MAX_NUM_CHANNELS];
    float * adapt[MAX_NUM_CHANNELS];
    float * fixed[MAX_NUM_CHANNELS];
    float * delaylines_next[MAX_NUM_CHANNELS];
    float * adapt_next[MAX_NUM_CHANNELS];
    float * fixed_next[MAX_NUM_CHANNELS];
    ECHOCOMP_PTR pEchocomp;
    int i, k, ch;
    int subband_start = 0;
    int subband_end = 0;
    int subband_start_next =0;
    int subband_end_next = 0;
    int NumFloatsInSubbandGroup = 0;
    int NumFloatsInSubbandGroup_next = 0;
    int SizeOfSubbandGroup = 0;
    int SizeOfSubbandGroup_next = 0;
    int FilterStartIndex_next = 0;
    short decay_computed = 0;
    int m = 0;

    void * IRAM_buf[MAX_NUM_CHANNELS];
    float * IRAM_buf_ad[MAX_NUM_CHANNELS][3];
    float * IRAM_buf_fx[MAX_NUM_CHANNELS][3];
    float * IRAM_buf_dl[MAX_NUM_CHANNELS][3];
    bool copied = false;

    /* dynamic groupsize for memxfer optimization */
    int NUM_SUBBAND_GROUPS = NUM_SUBBAND_GROUPS_STEREO; /*48kHz*/
    int subband_groups[NUM_SUBBAND_GROUPS_STEREO] = {14, 14, 15, 15, 16, 17, 17, 18, 20, 21, 23, 26, 30, 36, 38};

    for (ch = 0; ch < pEcho->channels; ch++)
    {
        delaylines_next[ch] = NULL;
        adapt_next[ch]      = NULL;
        fixed_next[ch]      = NULL;
    }

    pEchocomp = pEcho;

    adest  = (float*) tempdest;
    fxest  = (float*) adest + SUBUSED*2;
    adcanc = (float*) fxest + SUBUSED*2;
    fxcanc = (float*) adcanc + SUBUSED*2;

    adest_d = (double *) adest;
    fxest_d = (double *) fxest;
    adcanc_d = (double *) adcanc;
    fxcanc_d = (double *) fxcanc;

    pEchocomp->adaptMaxTap = 0.0f;
    pEchocomp->FilterStartIndex = 0;
    subband_end = subband_groups[0];

    for (i = subband_start; i < subband_end; i++)
    {
        NumFloatsInSubbandGroup += 2 * (pEchocomp->filtlen[i] + 1);
    }

    SizeOfSubbandGroup = sizeof(float) * NumFloatsInSubbandGroup;

    for (ch = 0; ch < pEcho->channels; ch++)
    {
        IRAM_buf[ch] = _mm_malloc(9 * SizeOfSubbandGroup, 64);

        IRAM_buf_dl[ch][0] = (float*) IRAM_buf[ch];
        IRAM_buf_dl[ch][1] = (float*) IRAM_buf_dl[ch][0] + NumFloatsInSubbandGroup;
        IRAM_buf_dl[ch][2] = (float*) IRAM_buf_dl[ch][1] + NumFloatsInSubbandGroup;
        
        IRAM_buf_ad[ch][0] = (float*) IRAM_buf_dl[ch][2] + NumFloatsInSubbandGroup;
        IRAM_buf_ad[ch][1] = (float*) IRAM_buf_ad[ch][0] + NumFloatsInSubbandGroup;
        IRAM_buf_ad[ch][2] = (float*) IRAM_buf_ad[ch][1] + NumFloatsInSubbandGroup;
        
        IRAM_buf_fx[ch][0] = (float*) IRAM_buf_ad[ch][2] + NumFloatsInSubbandGroup;
        IRAM_buf_fx[ch][1] = (float*) IRAM_buf_fx[ch][0] + NumFloatsInSubbandGroup;
        IRAM_buf_fx[ch][2] = (float*) IRAM_buf_fx[ch][1] + NumFloatsInSubbandGroup;

        delaylines[ch] = IRAM_buf_dl[ch][0];
        adapt     [ch] = IRAM_buf_ad[ch][0];
        fixed     [ch] = IRAM_buf_fx[ch][0];

        memcpy(delaylines[ch], pEchocomp->delaylines[ch], SizeOfSubbandGroup);
        memcpy(adapt     [ch], pEchocomp->adapt     [ch], SizeOfSubbandGroup);
        memcpy(fixed     [ch], pEchocomp->fixed     [ch], SizeOfSubbandGroup);
    }

    for ( k = 0; k < NUM_SUBBAND_GROUPS; k++ )
    {
        float adaptMaxTap_subband = 0.0f;
        m++;
        m %= 3;

        /* memxfer_load of buffers for k+1 */
        if ( k != NUM_SUBBAND_GROUPS - 1 )
        {
            subband_start_next = subband_end;
            if ( k == NUM_SUBBAND_GROUPS - 2 )
            {
                subband_end_next   = SUBUSED;
            }
            else
            {
                subband_end_next   = subband_start_next + subband_groups[k+1];
            }

            FilterStartIndex_next = pEchocomp->FilterStartIndex + NumFloatsInSubbandGroup;

            NumFloatsInSubbandGroup_next = 0;

            for (i = subband_start_next; i < subband_end_next; i++)
            {
                NumFloatsInSubbandGroup_next += 2 * (pEchocomp->filtlen[i] + 1);
            }

            SizeOfSubbandGroup_next = sizeof(float) * NumFloatsInSubbandGroup_next;

            for (ch = 0; ch < pEcho->channels; ch++)
            {
                delaylines_next[ch] = IRAM_buf_dl[ch][m];
                adapt_next     [ch] = IRAM_buf_ad[ch][m];
                fixed_next     [ch] = IRAM_buf_fx[ch][m];

                memcpy(delaylines_next[ch], pEchocomp->delaylines[ch] + FilterStartIndex_next, SizeOfSubbandGroup_next);
                memcpy(adapt_next     [ch], pEchocomp->adapt     [ch] + FilterStartIndex_next, SizeOfSubbandGroup_next);
                memcpy(fixed_next     [ch], pEchocomp->fixed     [ch] + FilterStartIndex_next, SizeOfSubbandGroup_next);
            }
        } /* end memxfer_load of buffers for k+1 */

        if ( decay_computed == 0 )
        {
            int decay_ix = pEchocomp->decAttData->decay_ix;

            if ( (subband_start <= decay_ix && decay_ix < subband_end) || decay_ix == -1 )
            {
                echocomp_find_decay(pEchocomp, pEchocomp->decAttData, fixed);
                decay_computed = 1;
            }
        }

        if ( pEchocomp->calcNewDeltaOn )
        {
            echocomp_calcNewDelta(pEchocomp,
                                  excgain,
                                  subband_start,
                                  subband_end);
        }

        if ( runShift )
        {
            for (ch = 0; ch < pEcho->channels; ch++)
            {
                echocomp_shift(delaylines[ch],
                               (ch == 0) ? lsBuf : lsBuf2,
                               pEchocomp->filtlen,
                               subband_start,
                               subband_end);

                memcpy(pEchocomp->delaylines[ch] + pEchocomp->FilterStartIndex,
                       delaylines[ch],
                       SizeOfSubbandGroup);
            }
        }

        adaptMaxTap_subband = nlms_filterStereo(pEchocomp,
                                                delaylines[0],
                                                delaylines[1],
                                                adapt[0],
                                                adapt[1],
                                                fixed[0],
                                                fixed[1],
                                                adest,
                                                fxest,
                                                subband_start,
                                                subband_end);

        echocomp_subtract(fxcanc,
                          micin,
                          fxest,
                          subband_start,
                          subband_end);

        echocomp_subtract(adcanc,
                          micin,
                          adest,
                          subband_start,
                          subband_end);

        echocomp_eval(&pEchocomp->fixedlevel[k],
                      fxcanc,
                      subband_start,
                      subband_end);

        echocomp_eval(&pEchocomp->adaptlevel[k],
                      adcanc,
                      subband_start,
                      subband_end);

        echocomp_eval(&pEchocomp->level[k],
                      micin,
                      subband_start,
                      subband_end);

        pEchocomp->copyd[k] = nlms_findBest(pEchocomp,
                                            adapt[0],
                                            adapt[1],
                                            fixed[0],
                                            fixed[1],
                                            k,
                                            subband_start,
                                            subband_end);

        if ( pEchocomp->copyd[k] == -3 )    /* Adapt set to Zero */
        {
            for (i = subband_start; i < subband_end; i++)

            {
                adest_d[i] = 0.0f;
                adcanc_d[i] = micin_d[i];
            }
        }
        else if ( pEchocomp->copyd[k] == -1 )   /* Fixed copied to adapt */
        {
            for (i = subband_start; i < subband_end; i++)
            {
                adest_d[i] = fxest_d[i];
                adcanc_d[i] = fxcanc_d[i];
            }
        }
        else if ( pEchocomp->copyd[k] == 1 )    /* Adapt copied to fixed */
        {

            for (i = subband_start; i < subband_end; i++)
            {
                fxest_d[i] = adest_d[i];
                fxcanc_d[i] = adcanc_d[i];
            }

        }

        for (ch = 0; ch < pEcho->channels; ch++)
        {
            memcpy(pEchocomp->fixed[ch] + pEchocomp->FilterStartIndex,
                   fixed[ch],
                   SizeOfSubbandGroup);
        }

        nlms_adaptStereo(pEchocomp,
                         pEchocomp->decAttData,
                         delaylines[0],
                         delaylines[1],
                         adapt[0],
                         adapt[1],
                         adcanc,
                         micin,
                         subband_start,
                         subband_end);

        for (ch = 0; ch < pEcho->channels; ch++)
        {
            memcpy(pEchocomp->adapt[ch] + pEchocomp->FilterStartIndex,
                   adapt[ch],
                   SizeOfSubbandGroup);
        }

        echocomp_weight(echoest,
                        fxest,
                        adest,
                        pEchocomp->fixedweight[k],
                        subband_start,
                        subband_end);

        for (ch = 0; ch < pEcho->channels; ch++)
        {
            fixed[ch] = fixed_next[ch];
            delaylines[ch] = delaylines_next[ch];
            adapt[ch] = adapt_next[ch];
        }
        subband_start = subband_start_next;
        subband_end = subband_end_next;
        pEchocomp->FilterStartIndex = FilterStartIndex_next;
        NumFloatsInSubbandGroup = NumFloatsInSubbandGroup_next;
        SizeOfSubbandGroup = sizeof(float) * NumFloatsInSubbandGroup;

        if ( (pEchocomp->copyd[k] == 1) && ( adaptMaxTap_subband > pEchocomp->adaptMaxTap) )
        {
            pEchocomp->adaptMaxTap = adaptMaxTap_subband;
            copied = true;
        }

    }

    if ( copied )
    {
        pEchocomp->fixedMaxTap = sqrtf(pEchocomp->adaptMaxTap);
    }

    if ( pEchocomp->debug == 11 )  /* print copy status for adapt */
        /* Printing a print-variable depending on what kind of copying is happening
         * Printing one value for each NUM_SUBBAND_GROUPS */
    {
        bool flag1 = true, flag2=true, flag3=true;
        static int cnt = 0;
        static int copyStatus[NUM_SUBBAND_GROUPS] = {0,0,0,0,0,0,0,0};

        for ( i=0; i<NUM_SUBBAND_GROUPS; i++ )
        {
            if ( pEchocomp->copyd[i] == -3 && flag1 )  /*Adapt to zero*/
            {
                copyStatus[i] = 4;
                flag1 = false;
            }
            if ( pEchocomp->copyd[i] == -1 && flag2 ) /*Fixed to adapt*/
            {
                copyStatus[i] = 2;
                flag2 = false;
            }
            if ( pEchocomp->copyd[i] == 1 && flag3 ) /*Adapt to fixed*/
            {
                copyStatus[i] = 1;
                flag3 = false;
            }
        }

        if ( ++cnt == 20 )
        {
            for ( i=0; i<NUM_SUBBAND_GROUPS; i++ )
            {
                printf("%3d, ",copyStatus[i]+i*10);
                copyStatus[i]=0;
            }
            printf("  \r\n");
            cnt = 0;
            flag1 = flag2 = flag3 = true;
        }
    }

    /*DEBUG prints*/
    if ( pEchocomp->debug == 41 )  /* look at fixedlevel in subbandgroups */
    {
        static float Sum0=0, Sum1=0, Sum2=0, Sum3=0, Sum4=0, Sum5=0, Sum6=0, Sum7=0;
        static int cnt = 0;
        Sum0 += pEchocomp->fixedlevel[0];
        Sum1 += pEchocomp->fixedlevel[1];
        Sum2 += pEchocomp->fixedlevel[2];
        Sum3 += pEchocomp->fixedlevel[3];
        Sum4 += pEchocomp->fixedlevel[4];
        Sum5 += pEchocomp->fixedlevel[5];
        Sum6 += pEchocomp->fixedlevel[6];
        Sum7 += pEchocomp->fixedlevel[7];

        if ( ++cnt == 20 )
        {
            printf("  %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n",
                   Sum0*0.05f, Sum1*0.05f,Sum2*0.05f,Sum3*0.05f,Sum4*0.05f, Sum5*0.05f,Sum6*0.05f,Sum7*0.05f);

            cnt = 0;
            Sum0 = 0.0f;
            Sum1 = 0.0f;
            Sum2 = 0.0f;
            Sum3 = 0.0f;
            Sum4 = 0.0f;
            Sum5 = 0.0f;
            Sum6 = 0.0f;
            Sum7 = 0.0f;
        }
    }
    if ( pEchocomp->debug == 42 )  /* look at adaptlevel in subbandgroups */
    {
        static float Sum0=0, Sum1=0, Sum2=0, Sum3=0, Sum4=0, Sum5=0, Sum6=0, Sum7=0;
        static int cnt = 0;
        Sum0 += pEchocomp->adaptlevel[0];
        Sum1 += pEchocomp->adaptlevel[1];
        Sum2 += pEchocomp->adaptlevel[2];
        Sum3 += pEchocomp->adaptlevel[3];
        Sum4 += pEchocomp->adaptlevel[4];
        Sum5 += pEchocomp->adaptlevel[5];
        Sum6 += pEchocomp->adaptlevel[6];
        Sum7 += pEchocomp->adaptlevel[7];

        if ( ++cnt == 20 )
        {
            printf("  %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e, %8.2e\n",
                   Sum0*0.05f, Sum1*0.05f,Sum2*0.05f,Sum3*0.05f,Sum4*0.05f, Sum5*0.05f,Sum6*0.05f,Sum7*0.05f);

            cnt = 0;
            Sum0 = 0.0f;
            Sum1 = 0.0f;
            Sum2 = 0.0f;
            Sum3 = 0.0f;
            Sum4 = 0.0f;
            Sum5 = 0.0f;
            Sum6 = 0.0f;
            Sum7 = 0.0f;
        }
    }

    for (ch = 0; ch < pEcho->channels; ch++)
    {
        _mm_free(IRAM_buf[ch]);
    }
#endif
    return 0;
}

/* Test functions */

void echocomp_setCalcNewDelta(ECHOCOMP_PTR pEchocomp, bool onoff)
{
    const char *status[] = {"off", "on"};
    pEchocomp->calcNewDeltaOn = onoff;
    printf("calcNewDelta is set %s \n\n", (char*)status[pEchocomp->calcNewDeltaOn]);
    return;
}

void echocomp_setMinDelta(ECHOCOMP_PTR pEchocomp, float value)
{
    pEchocomp->minDelta = value;
    printf("Echocomp->minDelta is set to %4.2ef \n\n", pEchocomp->minDelta);
    return;
}

void echocomp_setMaxDelta(ECHOCOMP_PTR pEchocomp, float value)
{
    pEchocomp->maxDelta = value;
    printf("Echocomp->maxDelta is set to %4.2ef \n\n", pEchocomp->maxDelta);
    return;
}

void echocomp_setStereoDelta(ECHOCOMP_PTR pEchocomp, float value)
{
    pEchocomp->nlmsDelta = value;
    printf("Echocomp->nlmsDelta is set to %4.2ef \n\n", pEchocomp->nlmsDelta);
    return;
}

void echocomp_setFilter(ECHOCOMP_PTR pEchocomp, int value)
{
    const char *status[] = {"default", "fixed", "adapt"};
    if ( value != 0 && value != 1 && value != 2 )
    {
        printf("value is not a valid option \n");
        value = 0;
    }
    pEchocomp->filter = value;
    printf("Echocomp->filter is set to %d: %s \n\n", pEchocomp->filter, (char*)status[pEchocomp->filter]);
    return;
}

static void echocomp_setZeroMono(ECHOCOMP_PTR pEchocomp, int value)
{
    ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[0];
    int i;

    if ( value == 1 || value == 0 )
    {
        for ( i=0; i<SUBUSED; i++ )
        {
            ECHOCOMP_SUBBAND *subband = &channel->subband[i];
            int j;

            for (j = 0; j < subband->filtlen; j++)
            {
               subband->fixed[j].re = 0.0f;
               subband->fixed[j].im = 0.0f;
            }
        }
        printf("Fixed filter set to zero \n\n");
    }

    if ( value == 2 || value == 0 )
    {
        for ( i=0; i<SUBUSED; i++ )
        {
            ECHOCOMP_SUBBAND *subband = &channel->subband[i];
            int j;

            for (j = 0; j < subband->filtlen; j++)
            {
               subband->adapt[j].re = 0.0f;
               subband->adapt[j].im = 0.0f;
            }
            subband->beta_1.re = 0.0f;
            subband->beta_1.im = 0.0f;
        }
        for ( i=0; i<NUM_SUBBAND_GROUPS; i++ )
        {
            pEchocomp->adaptlevel[i] = pEchocomp->level[i];
            pEchocomp->copyd[i] = -3;
        }
        printf("Adapt filter set to zero \n\n");
    }

    return;
}

static void echocomp_setZeroStereo(ECHOCOMP_PTR pEchocomp, int value)
{
    int i, ch;

    if (value == 0)
    {
        printf("Error: Cannot clear both fixed and adapt at the same time in stereo\n\n");
        return;
    }

    if (value == 1)
    {
        for (ch = 0; ch < pEchocomp->nChannels; ch++)
        {
            ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];

            for ( i=0; i<SUBUSED; i++ )
            {
                ECHOCOMP_SUBBAND *subband = &channel->subband[i];
                int j;

                for (j = 0; j < subband->filtlen; j++)
                {
                   subband->fixed[j].re = 0.0f;
                   subband->fixed[j].im = 0.0f;
                }
            }
        }
        printf("Fixed filter set to zero \n\n");
    }
    else if (value == 2)
    {
        for (ch = 0; ch < pEchocomp->nChannels; ch++)
        {
            ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];

            for ( i=0; i<SUBUSED; i++ )
            {
                ECHOCOMP_SUBBAND *subband = &channel->subband[i];
                int j;

                for (j = 0; j < subband->filtlen; j++)
                {
                   subband->adapt[j].re = 0.0f;
                   subband->adapt[j].im = 0.0f;
                }
            }
        }

        for (i = 0; i < NUM_SUBBAND_GROUPS; i++)
        {
            pEchocomp->adaptlevel[i] = pEchocomp->level[i];
            pEchocomp->copyd[i] = -3;
        }
        printf("Adapt filter set to zero \n\n");
    }
}

void echocomp_setZero(ECHOCOMP_PTR pEchocomp, int value)
{
    if (value != 0 && value != 1 && value != 2)
    {
        printf("No filter represented by this value \n\n");
        return;
    }

    if (pEchocomp->nChannels == 1)
    {
        echocomp_setZeroMono(pEchocomp, value);
    }
    else if (pEchocomp->nChannels == 2)
    {
        echocomp_setZeroStereo(pEchocomp, value);
    }
    else
    {
        printf("Invalid number of channels: %d\n", pEchocomp->nChannels);
    }
}

void echocomp_setDebug(ECHOCOMP_PTR pEchocomp, int value)
{
    pEchocomp->debug = value;
    printf("Echocomp->debug is set to %d \n\n", pEchocomp->debug);
    return;
}

void echocomp_status(ECHOCOMP_PTR pEchocomp)
{
    const char *status[] = {"off", "on"};
    const char *filter[] = {"default", "fixed", "adapt"};

    printf("\rStatus: Echocomp");
    printf("\r\n   calcnewdelta   - %s",     (char*)status[pEchocomp->calcNewDeltaOn]);
    printf("\r\n   minDelta       - %4.2ef", pEchocomp->minDelta);
    printf("\r\n   maxDelta       - %4.2ef", pEchocomp->maxDelta);
    printf("\r\n   nlmsDelta      - %4.2ef", pEchocomp->nlmsDelta);
    printf("\r\n   filter         - %s",     (char*)filter[pEchocomp->filter]);
    printf("\r\n   debug          - %d",     pEchocomp->debug);
    printf("\r\n");
}

/* Move tap values according to delay */
void echocomp_moveFilters(ECHOCOMP_PTR pEchocomp, const struct DELAY_ESTIMATION * delayEstimation, int delay_diff)
{
  int ch;
  int abs_delay = abs(delay_diff);

  if (delay_diff == 0)
    return;

  for (ch = 0; ch < pEchocomp->nChannels; ch++)
  {
    ECHOCOMP_CHANNEL *channel = &pEchocomp->channel[ch];
    int j;

    if ( delay_diff > 0 )   // Fill in appropriate lsbuffers in echocomp_delayline when delay decreases - "fast forward"  
    {
      int i;
      for (i = 0; i < abs_delay; i++)
      {
        __declspec(align(64)) COMPLEX32 lsFft[FFTSIZE/2];
        delayEstimation_returnIntermediateLsBuf(delayEstimation, lsFft, abs_delay, i);
        echocomp_shift(pEchocomp, ch, lsFft, 0, SUBUSED);
      }
    }

    for (j = 0; j < SUBUSED; j++)
    {
      ECHOCOMP_SUBBAND *subband = &channel->subband[j];

      int frames_to_move = subband->filtlen - abs_delay;
      int bytes_to_move = sizeof (COMPLEX32) * frames_to_move;
      int bytes_to_zero = sizeof (COMPLEX32) * abs_delay;

      if (bytes_to_move > 0)
      {
        if ( delay_diff > 0 )
        {
          memmove(subband->fixed + abs_delay, subband->fixed, bytes_to_move);
          memmove(subband->adapt + abs_delay, subband->adapt, bytes_to_move);

          memset (subband->fixed, 0, bytes_to_zero);
          memset (subband->adapt, 0, bytes_to_zero);
        }
        else if ( delay_diff < 0 )
        {
          int i;
          memmove(subband->fixed, subband->fixed + abs_delay, bytes_to_move);
          memmove(subband->adapt, subband->adapt + abs_delay, bytes_to_move);
          memmove(subband->delayline, subband->delayline + abs_delay, bytes_to_move);

          memset (subband->fixed + frames_to_move, 0, bytes_to_zero);
          memset (subband->adapt + frames_to_move, 0, bytes_to_zero);

          for ( i = 0; i < abs_delay ; i++)
          {
            COMPLEX32 lsSample;
            delayEstimation_returnHistoryLsSample(delayEstimation, &lsSample, delay_diff, subband->filtlen, j, i );
            subband->delayline[frames_to_move + i].re = lsSample.re;
          }
        }
      }
      else
      {
        bytes_to_zero = sizeof (float) * subband->filtlen;

        memset (subband->fixed, 0, bytes_to_zero);
        memset (subband->adapt, 0, bytes_to_zero);
        memset (subband->delayline, 0, bytes_to_zero);
      }

      //fixed_ptr += pEchocomp->filtlen[j]*2 + 2;
      //adapt_ptr += pEchocomp->filtlen[j]*2 + 2;
      //delayline_ptr += pEchocomp->filtlen[j]*2 + 2;
    }
  }
}


#ifdef UNITTEST
/* UNITTEST ************************************************************* */
/* Auth.: Jens Petter Stang (JPS)                                         */
/* Desc.: simple run-trhough of echocomp with 0.3sec random input on      */
/*        louds and micinput= filter(h, 1, louds); where h is a impulse   */
/*        response. The unittest verifies micest                          */
/* ********************************************************************** */
#include "unittest.h"
#include "testdata/echocomptestdata.h"  /* testvectors defined here */
#include "analyse_priv.h"

#pragma DATA_ALIGN(echocompTestMicInput48,8);
#pragma DATA_ALIGN(echocompTestLsInput48,8);

/* Private testdata */
static float echocompTestMicInput48[]  = {ECHOCOMP48_MICTESTINPUT_DEFINE};
static float echocompTestLsInput48[]  = {ECHOCOMP48_LSTESTINPUT_DEFINE};
static float echocompTestOutput48[] = {ECHOCOMP48_TESTOUTPUT_DEFINE};

/* compareResults ******************************************************* */
/* Auth.:                                                                 */
/* Desc.: Compares two float buffers. Due to compilator optimization, a   */
/*        a small difference 'epsilon' is allowed when comparing with the */
/*        testvector                                                     */
/* ********************************************************************** */
static int compareResults(float *inp1, float *inp2, int len)
{
    int differ = 0;
    int i;
    float epsilon = 1.5e-8f;

    for ( i=0;i<len && differ!=1; i++ )
    {
        if ( (fabsf(inp1[i] - inp2[i]))>epsilon )
        {
            differ = 1;
            //printf("%d \t %23.18f \t %23.18f \t %23.18f \n ",i,inp1[i],inp2[i], fabsf(inp1[i] - inp2[i]));
        }
    }

    return differ;
}

/* unittest_echocomp ****************************************************** */
/* Auth.:                                                                   */
/* Desc.:                                                                   */
/*                                                                          */
/* ********************************************************************     */
void unittest_echocomp()
{
    LSPROCESS_PTR pLs;
    ECHOCOMP_PTR pEchocomp;     //pointer to echocomp struct
    int   i;
    int   numFrames = 30;
    char  const * testName ={"Echocomp-48kHz"};
    float *pMicInput;
    float *pLsInput;
    float *pOutput;
    float *outputBuf;
    int differ = 0;         //criteria flag to be returned by function. 0: OK, 1: ERROR, 2:Nothing done, hence ERROR
    int createOk = 1;
    FILTERBANK_USE_TYPE filterbankUse = FILTERBANK_EC;

    /*
      FILE *fid;
      //const char *Filename ="d:\\tmp\\echocomptestoutput.dat";
      const char *Filename ="echocomptestoutput_host.dat";
      float * outbuf;
      int cnt;
    */

    /* Allocate memory */

    pLs = lsprocess_create(filterbankUse, 1);
    if ( pLs == NULL || pLs->echocompDelaylines[0] == NULL )
    {
            fprintf(stderr,"echocomp_unittest: Could not allocate lsprocess buffers\n");
            unittest_assert(0);
            return;
    }

    pEchocomp = echocomp_create(1);
    if ( pEchocomp==NULL )
    {
        createOk = 0;
        pEchocomp = (ECHOCOMP_PTR) malloc(sizeof(ECHOCOMP));
        if ( pEchocomp==NULL )
        {
            fprintf(stderr,"echocomp_unittest: Could not allocate echocomp struct buffer\n");
            unittest_assert(0);
            return;
        }

        pEchocomp->adapt[0] = (float*) malloc(sizeof(float)*SUBUSED * FILTLEN_MEAN * 2);
        if ( pEchocomp->adapt[0] == NULL )
        {
            fprintf(stderr,"echocomp_unittest: Could not allocate echocomp adapt buffer\n");
            unittest_assert(0);
            free(pEchocomp);
            return;
        }

        pEchocomp->fixed[0] = (float*) malloc(sizeof(float)*SUBUSED * FILTLEN_MEAN * 2);
        if ( pEchocomp->fixed[0] == NULL )
        {
            fprintf(stderr,"echocomp_unittest: Could not allocate echocomp fixed buffer\n");
            unittest_assert(0);
            //free(pEchocomp->adapt);
            free(pEchocomp);
            return;
        }

        //pEchocomp->delaylines = (float*) malloc(sizeof(float)*SUBUSED * FILTLEN_MEAN * 2);
    }

    /* Initialization */
    unittest_context(testName);
    lsprocess_init(pLs, filterbankUse);
    echocomp_init(pEchocomp, pLs->echocompDelaylines);

    // Use fixed delta = 1e-4 for now. When the delta calculation routine
    // and default values are decided upon, new testdata should be generated
    for ( i=0; i<SUBUSED; i++ )
    {
        pEchocomp->delta[i] = 1e-4;
        pEchocomp->r11[i] = 1e-4;
    }

    pEchocomp->calcNewDeltaOn = false;

    outputBuf = (float *) calloc(2*SUBUSED,sizeof(float)); //Allocating and initializing buffer
    //outbuf = (float *) calloc(numFrames*2*SUBUSED48,sizeof(float));

    differ = 0;

    for ( i=0; i<numFrames && differ==0; i++ )
    {
        pMicInput = echocompTestMicInput48 + 2*SUBUSED_START + i*FFTSIZE;
        pLsInput = echocompTestLsInput48 + i*FFTSIZE;
        pOutput = echocompTestOutput48 + i*2*SUBUSED;
    #ifdef PLATFORM_SNOOPY
        memcpy(pLs->pAnalyseSdram[0]->fftout, pLsInput, sizeof(float)*FFTSIZE);
        lsprocess_shift(pLs);   /* needing pLs->pAnalyseSdram->fftout */
        scratchmem_lock();
        echocomp_process(pEchocomp,pLs->pAnalyseSdram[0]->fftout, outputBuf, pMicInput, 1.0f, true);
    #else
        scratchmem_lock();
        echocomp_process(pEchocomp, pLsInput, outputBuf, pMicInput, 1.0f, true);
    #endif
        scratchmem_unlock();
        differ = compareResults(outputBuf, pOutput, 2*SUBUSED);


        //memcpy(outbuf+i*2*(SUBUSED48), outputBuf, sizeof(float)*2*SUBUSED48);
    }

    /*
      printf("write echoest to file \n");
      fid = (void*)fopen(Filename,"wb");
      if (fid == NULL)
      {
      tt_fatal(TTLOC,"Could not open file output");
      }
      cnt = fwrite(outbuf,sizeof(float),numFrames*2*SUBUSED48,fid);
      printf("cnt: %d \n",cnt);
      fclose(fid);
      free(outbuf);
      (void)cnt;
    */

    unittest_assert(differ == 0);
    free(outputBuf);


    if (createOk)
    {
        echocomp_destroy(pEchocomp);
        lsprocess_destroy(pLs, filterbankUse);
    }
    else
    {
        free(pEchocomp->adapt[0]);
        free(pEchocomp->fixed[0]);
        free(pEchocomp->delaylines[0]);
        free(pEchocomp);
    }


}

#endif /* UNITTEST */
