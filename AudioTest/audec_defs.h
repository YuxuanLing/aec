#ifndef AUDEC_DEFS_H
#define AUDEC_DEFS_H

/*-------------------------------------------------------------------------
 * -- PREPROCESSOR CONSTANTS AREA
 *-------------------------------------------------------------------------*/
#define nWRITE_DEBUG_TO_FILE
#define nWRITE_NOIREDDEBUG_TO_FILE

#ifdef WRITE_DEBUG_TO_FILE
#include <stdio.h>
FILE * fidEc;
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define EPSILON     1e-35f         /* to avoid division by zero, add epsilon */

/* general */
#define FFTSIZE 768              /* length of fft, 1536 real numbers */
#define FRAMESIZE 480            /* 48kHz input soundbuffer has this length */
#define SUBUSED 320              /* subbands used for 48kHz*/
#define SUBUSED_START 1            /* start-index subbands used in relation to a full subbandbuffer */

#define MAX_NUM_CHANNELS 2

#define PRESERVE_ORIGINAL_BUGS

typedef struct COMPLEX32
{
  float re; // real part
  float im; // imaginary part
} COMPLEX32;


/* size of echoestimation filter: */
/* SUBUSED elements for different filterlengths in each subband */
/* numbers in each row must be the same because members of one DECAY_SUBBAND_GROUP */
/* are assumed to have the same number of taps in computation of decay_attack_gain */
#define FILTLEN_TAPS48  34, 34, 34, 34, 34, 34, 34, 34, \
                        33, 33, 33, 33, 33, 33, 33, 33, \
                        33, 33, 33, 33, 33, 33, 33, 33, \
                        32, 32, 32, 32, 32, 32, 32, 32, \
                        31, 31, 31, 31, 31, 31, 31, 31, \
                        31, 31, 31, 31, 31, 31, 31, 31, \
                        30, 30, 30, 30, 30, 30, 30, 30, \
                        29, 29, 29, 29, 29, 29, 29, 29, \
                        29, 29, 29, 29, 29, 29, 29, 29, \
                        28, 28, 28, 28, 28, 28, 28, 28, \
                        27, 27, 27, 27, 27, 27, 27, 27, \
                        27, 27, 27, 27, 27, 27, 27, 27, \
                        26, 26, 26, 26, 26, 26, 26, 26, \
                        25, 25, 25, 25, 25, 25, 25, 25, \
                        25, 25, 25, 25, 25, 25, 25, 25, \
                        24, 24, 24, 24, 24, 24, 24, 24, \
                        23, 23, 23, 23, 23, 23, 23, 23, \
                        23, 23, 23, 23, 23, 23, 23, 23, \
                        22, 22, 22, 22, 22, 22, 22, 22, \
                        21, 21, 21, 21, 21, 21, 21, 21, \
                        21, 21, 21, 21, 21, 21, 21, 21, \
                        20, 20, 20, 20, 20, 20, 20, 20, \
                        19, 19, 19, 19, 19, 19, 19, 19, \
                        19, 19, 19, 19, 19, 19, 19, 19, \
                        18, 18, 18, 18, 18, 18, 18, 18, \
                        17, 17, 17, 17, 17, 17, 17, 17, \
                        17, 17, 17, 17, 17, 17, 17, 17, \
                        16, 16, 16, 16, 16, 16, 16, 16, \
                        15, 15, 15, 15, 15, 15, 15, 15, \
                        15, 15, 15, 15, 15, 15, 15, 15, \
                        14, 14, 14, 14, 14, 14, 14, 14, \
                        13, 13, 13, 13, 13, 13, 13, 13, \
                        13, 13, 13, 13, 13, 13, 13, 13, \
                        12, 12, 12, 12, 12, 12, 12, 12, \
                        11, 11, 11, 11, 11, 11, 11, 11, \
                        11, 11, 11, 11, 11, 11, 11, 11, \
                        11, 11, 11, 11, 11, 11, 11, 11, \
                        11, 11, 11, 11, 11, 11, 11, 11, \
                        11, 11, 11, 11, 11, 11, 11, 11, \
                        11, 11, 11, 11, 11, 11, 11, 11


#define FILTLEN_MAX 36 /* max value of FILTLEN-vector + 1, rounded up to a multiple of 4 (APA) */
#define FILTLEN_MIN 11
#define FILTLEN_MEAN 23 /* ceil(mean(FILTLEN_TAPS48) + 1)*/

/* pre- postfilter */
#define PREFIRLEN 2                      /* length of pre- and postfilter */
#define PFT48_COEFF2 -0.8125f           /* this means that the pre-filter B= /post-filter A= [1, -0.8125] */
#define PRE48_GAIN    2.0f              /* gainadjustment of pre-filter, B = [b0, b1]*gain| */
#define PST48_GAIN    0.5f              /* = 1/PRE48_GAIN, gainadjustment of post-filter, A = [a0, a1]*gain */

/* analyse */
#define AFIRLEN48 4608      /* length of the 48kHz analysefilter */

/* impulse response decay estiamtion */
#define SUBUSED_FINDDECAY48 159; /* The finddecay routine requires at least 20 filtertaps */

#define DECAY_DEFAULT 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, \
                      0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, \
                      0.94f, 0.94f, 0.94f, 0.94f, 0.94f, 0.94f, 0.94f, 0.94f, \
                      0.93f, 0.93f, 0.93f, 0.93f, 0.93f, 0.93f, 0.93f, 0.93f, \
                      0.93f, 0.93f, 0.93f, 0.93f, 0.93f, 0.93f, 0.93f, 0.93f, \
                      0.92f, 0.92f, 0.92f, 0.92f, 0.92f, 0.92f, 0.92f, 0.92f, \
                      0.91f, 0.91f, 0.91f, 0.91f, 0.91f, 0.91f, 0.91f, 0.91f, \
                      0.91f, 0.91f, 0.91f, 0.91f, 0.91f, 0.91f, 0.91f, 0.91f, \
                      0.90f, 0.90f, 0.90f, 0.90f, 0.90f, 0.90f, 0.90f, 0.90f, \
                      0.89f, 0.89f, 0.89f, 0.89f, 0.89f, 0.89f, 0.89f, 0.89f, \
                      0.89f, 0.89f, 0.89f, 0.89f, 0.89f, 0.89f, 0.89f, 0.89f, \
                      0.88f, 0.88f, 0.88f, 0.88f, 0.88f, 0.88f, 0.88f, 0.88f, \
                      0.87f, 0.87f, 0.87f, 0.87f, 0.87f, 0.87f, 0.87f, 0.87f, \
                      0.87f, 0.87f, 0.87f, 0.87f, 0.87f, 0.87f, 0.87f, 0.87f, \
                      0.86f, 0.86f, 0.86f, 0.86f, 0.86f, 0.86f, 0.86f, 0.86f, \
                      0.85f, 0.85f, 0.85f, 0.85f, 0.85f, 0.85f, 0.85f, 0.85f, \
                      0.85f, 0.85f, 0.85f, 0.85f, 0.85f, 0.85f, 0.85f, 0.85f, \
                      0.84f, 0.84f, 0.84f, 0.84f, 0.84f, 0.84f, 0.84f, 0.84f, \
                      0.83f, 0.83f, 0.83f, 0.83f, 0.83f, 0.83f, 0.83f, 0.83f, \
                      0.83f, 0.83f, 0.83f, 0.83f, 0.83f, 0.83f, 0.83f, 0.83f, \
                      0.82f, 0.82f, 0.82f, 0.82f, 0.82f, 0.82f, 0.82f, 0.82f, \
                      0.81f, 0.81f, 0.81f, 0.81f, 0.81f, 0.81f, 0.81f, 0.81f, \
                      0.81f, 0.81f, 0.81f, 0.81f, 0.81f, 0.81f, 0.81f, 0.81f, \
                      0.80f, 0.80f, 0.80f, 0.80f, 0.80f, 0.80f, 0.80f, 0.80f, \
                      0.79f, 0.79f, 0.79f, 0.79f, 0.79f, 0.79f, 0.79f, 0.79f, \
                      0.79f, 0.79f, 0.79f, 0.79f, 0.79f, 0.79f, 0.79f, 0.79f, \
                      0.78f, 0.78f, 0.78f, 0.78f, 0.78f, 0.78f, 0.78f, 0.78f, \
                      0.77f, 0.77f, 0.77f, 0.77f, 0.77f, 0.77f, 0.77f, 0.77f, \
                      0.77f, 0.77f, 0.77f, 0.77f, 0.77f, 0.77f, 0.77f, 0.77f, \
                      0.76f, 0.76f, 0.76f, 0.76f, 0.76f, 0.76f, 0.76f, 0.76f, \
                      0.75f, 0.75f, 0.75f, 0.75f, 0.75f, 0.75f, 0.75f, 0.75f, \
                      0.75f, 0.75f, 0.75f, 0.75f, 0.75f, 0.75f, 0.75f, 0.75f, \
                      0.74f, 0.74f, 0.74f, 0.74f, 0.74f, 0.74f, 0.74f, 0.74f, \
                      0.73f, 0.73f, 0.73f, 0.73f, 0.73f, 0.73f, 0.73f, 0.73f, \
                      0.73f, 0.73f, 0.73f, 0.73f, 0.73f, 0.73f, 0.73f, 0.73f, \
                      0.72f, 0.72f, 0.72f, 0.72f, 0.72f, 0.72f, 0.72f, 0.72f, \
                      0.71f, 0.71f, 0.71f, 0.71f, 0.71f, 0.71f, 0.71f, 0.71f, \
                      0.71f, 0.71f, 0.71f, 0.71f, 0.71f, 0.71f, 0.71f, 0.71f, \
                      0.70f, 0.70f, 0.70f, 0.70f, 0.70f, 0.70f, 0.70f, 0.70f, \
                      0.69f, 0.69f, 0.69f, 0.69f, 0.69f, 0.69f, 0.69f, 0.69f

#define DECAY_MIN 0.75f
#define DECAY_MAX 0.984375f
#define DECAY_16_MIN 1.00225957576185464859e-2f /* DECAY_MIN^16 */
#define DECAY_16_MAX 0.7772651709434128567135f /* DECAY_MAX^16 */
#define DECAY_START_TAP 4 /* time distance between loudspeaker and microphone must be less than DECAY_START_TAP, meaning actual distance less than 340m/s * 40ms = 13.6m */
#define DECAY_SUBBAND_GROUPSIZE 8 /* number of subbands considered as one group with one common decayestimate */
#define DECAY_SUBBAND_GROUPLENGTH 8 /* number of 10ms packets to use in computation of early and late power of impulse response*/
#define DECAY_ATTACK_TAPS_DEFAULT 0.2f, 0.5f, 3.0f, 3.0f, 3.0f, 2.8125f, 2.6367f, 2.4719f, 2.3174f, 2.1725f /* default value must be  NUMBER_OF_ATTACK_TAPS long */
#define DECAY_NUMBER_OF_USED_ATTACK_TAPS_DEFAULT 7
#define DECAY_LARGEST_TAP_AT_MAX 8 /* largest value of echo filter (total power of subbands) occur at maximum this 10ms packet */
#define DECAY_LARGEST_TAP_AT_MIN 3 /* largest value of echo filter (total power of subbands) occur at minimum this 10ms packet */
#define DECAY_PRE_ROOM_DOMINACE_TIME 2 /* number of 10ms packets from largest value of echo filter occur until the room dominates the echo filter */
#define NUMBER_OF_ATTACK_TAPS (DECAY_LARGEST_TAP_AT_MAX + DECAY_PRE_ROOM_DOMINACE_TIME)/* maximum number of individual tapweights */
#define MAX_ADAPT_TAPWEIGHT_CHANGE 0.1f /* -20dB maximum tolerated realtive difference between two consecutive attack_taps */
#define MAX_ADAPT_TAPWEIGHT_CHANGE_TOTAL 0.001f /* -60dB maximum tolerated realtive difference between smallest and largest attack_taps */

#define NUM_SUBBAND_GROUPS 8 //used in echocomp_init, the subband groups are of different sizes due to effective memxfer transfers
#define NUM_SUBBAND_GROUPS_STEREO 15

/* exceptionhandler constants */
#define EXCSUBUSED_START (6 - SUBUSED_START)             /* start-index exc-subbands used in relation to a "subused" subbandbuffer */
#define EXCSUBUSED       12               /* uses only some subbands to reduce cpu load */
#define EXCSUBUSED_END   (EXCSUBUSED_START + EXCSUBUSED) /* stop-index exc-subbands used in relation to a "subused" subbandbuffer */

#if 1   /* IIR implementation */
#define EXCLAGCNT        10        /* exception handler history */
#define EXCGAINFILTSP    0.25f     /* 1/4     gain update speed */
#define EXCPSPD          0.0625f   /* 2^(-4)  power update speed*/
#define EXCMINLSPOW      3.814697266e-6f  /* 1/(64*4096) minimum ls power to evaluate exception */
#define EXCCAFTNOEXCLAST 0.0625f   /* 1/16    after  echo level correlation correction, if exception did not occur last frame */
#define EXCCBEFNOEXCLAST 0.0f      /* 0       before echo level correlation correction, if exception did not occur last frame */
#define EXCCAFTEXCLAST   0.0f      /* 0       after  echo level correlation correction, if exception did occur last frame */
#define EXCCBEFEXCLAST   0.375f
#else   /* FIR implementation*/
#define EXC_FIR_LAG         4     /* number of different lags used in FIR correlation calcluation  */
#define EXC_LAGSHIFT        0     /* number of samples in subband ls signal is shifted before exc calc begins:4-5 optimal for FIR */
#define EXCLAGCNT           32    /* exception handler history */
#define EXCLAGCNT_ATTACK    16    /* exception handler history for short XCORR calculation: attack : must be =< EXCLAGCNT  */
#define EXC_DELAYCALC       12    /* EXC_DELAYCALC must be greater than EXC_FIR_LAG*/
#define EXC_DELAY_SUBUSED   4     /* EXC_DELAY_SUBUSED must be less or equal to EXCSUBUSED*/
#define EXC_DELAYPOWSPD     0.015625f
#define EXC_DELAY_SWITCHCNT 5     /*number of concecutive maxvalues before setting new delay*/
#define EXCGAINFILTSPATT    0.6f  /* 1/4     gain update speed */
#define EXCGAINFILTSPDEC    0.15f /* 1/4     gain update speed */
#define EXCCAFTNOEXCLAST    0.05f /* 1/16    after  echo level correlation correction, if exception did not occur last frame */
#define EXCCBEFNOEXCLAST    0.0f  /* 0       before echo level correlation correction, if exception did not occur last frame */
#define EXCCAFTEXCLAST      0.0f  /* 0       after  echo level correlation correction, if exception did occur last frame */
#define EXCCBEFEXCLAST      0.3f
#endif

/* Non-linear processing setup */
#define NLP_SUBUSED 48                /* voiceband used in echo-reduction */
#define NLP_SUBUSED_START 5           /* start-index nlp subbands used in relation to a "subused" subbandbuffer (aftlev etc.)*/
#define NLP_SUBUSED_END (NLP_SUBUSED_START + NLP_SUBUSED)  /* stop-index nlp subbands used in relation to a "subused" subbandbuffer (aftlev etc.)*/
#define NLPSPEED  0.25f              /* 2^(-2) NB! Movi specific; Original value was 2^(-3) */
#define NLPC_NOILEV 0.9375f           /* 15/16; underestimation of noilev */
#define NLPC_MONO   2.0f              /* fullg to be increased by factor: 2 (in noisereduction.c)*/
#define NLPC_STEREO 1.5f              /* fullg to be increased by factor: 1.5 (in noisereduction.c) if stereo */
#define SGTLIM    0.125f              /* 1/8;       */
#define NLP_SLOW_SPD 0.01f         /* slow update when beflev is falling */
#define MAXTWIST  4

/* noisereduction */
#define NRSPEED 0.25f              /* 2^(-2); noise reduction level calculation update speed */
#define ISNOISE_INCSPEED 0.03125f  /* 2^(-5); noise estimate increase speed, if noise detector indicates noise */
#define ISNOISE_DECSPEED 0.03125f  /* 2^(-5); noise estimate decrease speed, if noise detector indicates noise */
#define NONOISE_INCSPEED 2.384185791e-7f /* 2^(-22); noise estimate increase speed, if noise detector doesn't indicate noise */
#define NONOISE_DECSPEED 3.814697266e-6f /* 2^(-18); noise estimate decrease speed, if noise detector doesn't indicate noise */
#define NRGLMAX  0.5011872336f     /* 10^( -6/20); minimum noise reduction gain if much noise */
#define NRGLMIN  0.1778279410f     /* 10^(-15/20); minimum noise reduction gain if little noise */
#define NRGL_UPDSPEED 0.00390625f  /* (1/256); */
#define NOISE_OVEREST 1.333333333f /* 4/3; noise overestimation factor */
#define NRGSP 0.75f                /* 3/4; noise reduction gain speed */
#define NOISEDLIM 0.0625f          /* 2^(-4) noise variation limit */
#define NOISECLIM 30               /* noise duration limit */

/* mls "random" number generator */
#define MLS_POLYNOME 0x20040       /* only certain (binary-representations of) polynomes can be used, this one repeats itself after 262143 iterations */
#define MLS_RAND_MAX 262143.0f     /* defined by the polynome (i.e. the maximum value the mls_rand kan take.) */
#define MLS_RAND_MAX_INV 3.814711816e-6f                        /* 1/MLS_RAND_MAX */

/* synth */
#define SFIRLEN48 3360 /* the 48kHz-filter is 3072, but with zeropadding; 3360 */

/* echocomp */
/*#define SSBUSP 0.125f*/
#define SSBUSP 0.0625f

#define LOOPDETECT_CNT 6
#define LOOPDETECT_START 3
#define LOOPDETECT_SUBUSED 3, 7, 15, 80, 95, 114   /*ecMixPtr->micData[0]->mic_fft[m] is size FFTSIZE*/
#define LOOPDETECT_SUB_NOILEV 2, 6, 14, 79, 94, 113 /*nrPtr->noilev[m] is size SUBUSED*/
#define LOOPDETECT_MLSLEN 31
#define LOOPDETECT_MAX_IMPRES 4
#define LOOPDETECT_MAX_IMPRES_INV 1.0f/LOOPDETECT_MAX_IMPRES
#define LOOPDETECT_CONVLEN LOOPDETECT_MLSLEN*LOOPDETECT_MAX_IMPRES
#define LOOPDETECT_DELAYLINE 20

#define LOOPDETECT_COMBINATIONS 6//(LOOPDETECT_MAX_IMPRES-1)*(LOOPDETECT_MAX_IMPRES-2)
#define LOOPSUB1 24
#define LOOPSUB2 30
#define LOOPSUB3 36
#define LOOPSUB4 65
#define LOOPSUB5 76
#define LOOPSUB6 90


#define LOOPMICSPEED  0.25f
//#define CONVLENFAC = 1/LOOPDETECT_CONVLEN

#define LSLIMITER_MAXCHANNELS 8




#if 0  /* old defines */
#define FRAMESIZE16 160            /* used in adjustments.c */
#define SUBUSED16 116
#define FFTSIZE16 256              /* length of fft, 512 real numbers */

#define PFT16_COEFF2 -0.625f            /* this means that the pre-filter B= /post-filter A= [1, -0.625] */
#define PRE16_GAIN    1.0f              /* gainadjustment of pre-filter, B = [b0, b1]*gain */
#define PST16_GAIN    1.0f              /* = 1/PRE16_GAIN, gainadjustment of post-filter, A = [a0, a1]*gain */
#define AFIRLEN16 1536      /* length of the 16kHz analysefilter */
#define SUBUSED_FINDDECAY16 116; /* The finddecay routine requires at least 20 filtertaps */
#define SFIRLEN16 1120 /* the 16kHz-filter is 1024, but with zeropadding; 1120 */

/* loudspeaker delay compensation */
#define MAXEXTRALSDELAY 100 /* in milliseconds should be a muliple of 10ms */
#define LS16DELAYBUFFERLENGTH (FRAMESIZE16*(MAXEXTRALSDELAY/10)) /* 10ms * 10 = 100ms maximum delay to compensate */

/* dc removal */
/* #define AudEC_DCTAP 1.0f/(256*FRAMESIZE48) */

/* loudspeakerprocessing */
#define LSLEVSPEED 0.125f            /* 1/8 update speed used by adaption ctrl */
#define ADLIMIT    7.62939453125e-6f /* 2^(-17) subband loudspeaker adaption level limit */

/* 48kHz constants */
#define DECIM_FIRLEN 90   /* consistency with file "decimfir90.h" */
#define DELAYBEFORECANC (((AFIRLEN-1)+(PREFIRLEN-1))/2 + FRAMESIZE16)  /* consistency with file "aFirData.h" */
#define DELAYBEFORECANC48 ((DELAYBEFORECANC) * DECIMFACTOR + DECIM_FIRLEN/2 - 1) /* (2828 smpl)  decimfir even length, 1/2 sample delay
                                                                                   accounted for in "DELAYBEFORECANC". */
//#define DELAYAFTERCANC ((SFIRLEN-1)+(PREFIRLEN-1))/2       /* consistency with file "decimfir90.h" */
#define DELAYAFTERCANC (((1024-1)+(PREFIRLEN-1))/2  - FRAMESIZE16)  /* consistency with file "sFirData.h" */
#define DELAYAFTERCANC48 ((DELAYAFTERCANC) * DECIMFACTOR + DECIM_FIRLEN/2 ) /*(2061 smpl)  decimfir even length, 1/2 sample delay
                                                                              accounted for in "DELAYBEFORECANC". */
#define TOTALDELAY (DELAYBEFORECANC + DELAYAFTERCANC)
#define TOTALDELAY48 (DELAYBEFORECANC48 + DELAYAFTERCANC48) /* (4889 smpl / 101,85 ms) (in addition comes ping-pong buffer of 480 smpl/10ms.)*/
#define AUDEC_MIC_DELAY_BUF FRAMESIZE
#define AUDEC_MIC_DELAY_TOT (TOTALDELAY48+AUDEC_MIC_DELAY_BUF)
#define AUDEC_MIC_DELAY_TOT_MS (AUDEC_MIC_DELAY_TOT*1000/48000)

#endif


#endif /*AUDEC_DEFS_H*/

