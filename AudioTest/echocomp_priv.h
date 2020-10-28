#ifndef ECHOCOMP_PRIV_H
#define ECHOCOMP_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 *         Geir Ole Øverby (GEO), Tandberg Telecom AS.
 *
 * Description   : Private data for echocomp module.
 **************************************************************************/

#include "audec_defs.h"
#include <stdbool.h>

struct DECAY_ATTACK
{
    float decay_early_power;             /* temporary sum of early impulse response power */
    float decay_late_power;              /* temporary sum of late  impulse response power */
    short decay_ix;                      /* current subband to compute (early and late) filter power of */
    char  decay_compute;                 /* 0 when summing power of filter, 1 when finished, then computes decay estimate */
    int   subused_finddecay;
};

typedef struct ECHOCOMP_SUBBAND          /* per-subband data */
{
    short filtlen;                       /* length of adaptive and fixed filter */
    float delta;                         /* delta */
    float r11;                           /* r11 */
    /* APA variables */
    float r22;                           /* r22 = previous r11 */
    COMPLEX32 r12;                       /* r12 (real and imag) */
    COMPLEX32 FiltCorr;                  /* Correction for h->z transform (real and imag) */
    COMPLEX32 beta_1;                    /* Beta1 (real and imag) */
    COMPLEX32 beta_fixed;                /* Beta_fixed (real and imag) */
    COMPLEX32 est_err;                   /* previous estimation error */

    /* pointers to per-subband filter taps */
    COMPLEX32* delayline;
    COMPLEX32* adapt;
    COMPLEX32* fixed;

    // per-subband decay attack
    float decay_attack_gain;             /* gain to adjust decay_attack_taps with */
    float decay;                         /* decay-factor between two consecutive filtertaps in fixed filter */
    float decay_default;                 /* default decay-factor between two consecutive filtertaps in fixed filter */
    float adaptWeights[FILTLEN_MAX];
    int   weights_need_recalculate;      /* optimization, don't recompute weights on every pass */

#if defined(HAVE_IOS)
    /* pointers to per-subband filter taps */
    COMPLEX32* delayline_base;
#endif

} ECHOCOMP_SUBBAND;

typedef struct ECHOCOMP_CHANNEL
{
    ECHOCOMP_SUBBAND subband[SUBUSED];   /* each channel is divided into 320 sub-bands */

    // decay attack
    float decay_attack_taps_sum_new;     /* will be set to the sum of the used part of attack_taps_new */
    float decay_attack_taps_sum_old;     /* will be set to the sum of the used part of attack_taps_old */
    float attack_taps_new[NUMBER_OF_ATTACK_TAPS]; /* array of tap weights before room dominates filter */
    float attack_taps_old[NUMBER_OF_ATTACK_TAPS]; /* array of tap weights before room dominates filter */
    int   number_of_used_attack_taps_new;/* the contents of attack_taps[] after this number are ignored */
    int   number_of_used_attack_taps_old;/* the contents of attack_taps[] after this number are ignored */
    float filter_power[FILTLEN_MAX];     /* power of taps in filter in each 10ms packet */
} ECHOCOMP_CHANNEL;

/* ----------------------------------------------------*/
/* defines 'struct echocomp' which contains            */
/* variables etc. that define the echoestimation       */
/* filters. These should be filled by an init routine. */
/* ----------------------------------------------------*/
typedef struct ECHOCOMP
{
    ECHOCOMP_CHANNEL channel[MAX_NUM_CHANNELS];
    int   nChannels;

    float fixedweight[NUM_SUBBAND_GROUPS_STEREO];   /* amount of fixed filter used (adaptweight = 1-fixedweight) */
    float level[NUM_SUBBAND_GROUPS_STEREO];         /* level of uncaceled mic signal*/
    float fixedlevel[NUM_SUBBAND_GROUPS_STEREO];    /* level of canceled mic signal if fixed filter is used */
    float adaptlevel[NUM_SUBBAND_GROUPS_STEREO];    /* level of canceled mic signal if adapt filter is used */
    float fixedMaxTap;                              /* max Tap of fixed filter */
    float adaptMaxTap;                              /* max Tap^2 of adapt filter */
    float aerlInverse;
    short copyd[NUM_SUBBAND_GROUPS_STEREO];         /* tells if copying between adapt and fixed filters or reset of adapt has occured */
    float mu;                                       /* stepsize */

    DECAY_ATTACK decAttData;

    /* debug parameters */
    bool  calcNewDeltaOn;
    float minDelta;
    float maxDelta;
    int   filter;
    int   debug;
    float nlmsDelta;

#ifdef __INTEL_COMPILER
	int MMX_PRESENT;
	int SSE_PRESENT;
	int SSE2_PRESENT;
	int SSE3_PRESENT;
	int SSSE3_PRESENT;
	int SSE4_PRESENT;
#endif
} ECHOCOMP;

#endif
