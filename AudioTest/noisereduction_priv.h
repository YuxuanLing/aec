#ifndef NOISEREDUCTIONS_PRIV_H
#define NOISEREDUCTIONS_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 * Co-author     : -
 *
 * Switches      : <COMPILER SWITCH HERE>
 *         <add description here>
 *
 * Description   : <add description here>
 *
 * Note    : <notes here>
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *       Modification made.
 *
 * yyyy-mm-dd    : <signature>
 *       File created.
 *
 **************************************************************************/

#include "audec_defs.h"
#include "dereverb.h"
#include <stdbool.h>
#include "keyclickremoval.h"

/*-------------------------------------------------------------------------
 * -- EXPORTED MACROS AREA
 *-------------------------------------------------------------------------*/

/*
#define NOISEREDUCTIONS_PRIV_SET_NCNT_BOOT(noired_ptr, val)      \
    ((noired_ptr)->ncnt_boot = (val))
*/

/* ---------------------------------------------------*/
/* defines 'struct NOIRED' which contains             */
/* variables etc. that define the noisereduction      */
/* filters.                                           */
/* ---------------------------------------------------*/
typedef struct NOIRED {
    float beflev[SUBUSED];  /* level before echo compensation */
    float aftlev[SUBUSED];  /* level after echo compensation */
    float estlev[SUBUSED];  /* level of estimated echo */
    float noilev[SUBUSED];  /* est. noise level */
    float totlev[SUBUSED];  /* total level */
    float noisgn[SUBUSED];  /* noise gain vector  */

    float sum_beflev;
    float sum_aftlev;
    float sum_estlev;
    float sum_noilev;
    float sum_totlev;

    float subgPrev[SUBUSED];        /* subgain delay-line */
    float fullgPrev;                /* fullgain delay-line */
    float beflevPrev[2];            /* beflev delay-line (fullband) */
    float fullgPrevHi;              /* fullgain delay-line */
    float beflevPrevHi[2];          /* beflev delay-line (fullband) */
    float nrgoptimal;               /* optimal noise reduction limit */
    float nrglimit;                 /* noise regulation limit  */
    float noiseamp;                 /* noiseamplification */
    float comNoiseAmp;              /* comNoise amplification */
    short ncnt;                     /* noisedetectioncounter  */
    short ncnt_clim;                /* noisecount limit */
    short ncnt_boot;                /* noisecount from boot */
    unsigned int seed;              /* the random number "seed" can be set to any value except 0 */
    float nlpspeed;
    float shellnlp;

     /* test variables */
    bool NlpGainOn;
    bool dereverbOn;
    bool ComNoiseOn;
    bool NoiseRedOn;
    bool NLPSubOn;
    bool NlpHiOn;
    bool fullgNlpOn;
    bool fullgNlpHiOn;
    bool extragainOn;
    int  nlpSubusedEndHi;
    int  fullgNlpTransition;
    int  debug;

#ifdef __INTEL_COMPILER
	int MMX_PRESENT;
	int SSE_PRESENT;
	int SSE2_PRESENT;
	int SSE3_PRESENT;
	int SSSE3_PRESENT;
	int SSE4_PRESENT;
#endif



} NOIRED;


#endif /* NOISEREDUCTIONS_PRIV_H */
