#ifndef LSLIMITER_PRIV_H
#define LSLIMITER_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2009
 *==========================================================================
 *
 * Author        : Ole J Gauteplass
 * Co-author     : 
 *
 * Switches      : -
 *
 * Description   : Private data for the lslimiter module
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
 **************************************************************************/
#include "audec_defs.h"
#include <stdbool.h>

#define LSLIM_CLIP_HEADROOM 0.89f
#define LSLIM_ATT (1.0f/(LSLIM_DELAY/3.1415f)) // 0.013f /* Pi is such a nice number! */
#define LSLIM_REL 0.00041667f /* 1/(fs*.05) */
#define AUD_MAX_NUM_ACOUSTIC_INPUTS 4

typedef struct LSLIMITER {

  /* private members */
  float gain_current_vector[2]; 
  float stereoWeight;
  float stereoWeightOld;
  int forceMonoCntr;
  float limiter_threshold; /* max output to loudspeaker */

  int limcount;
  float oldsp;
  float limoldg;

  float lslimDelta; /* Productdependent aerl-constant */

  /* debug variables */
  float lslimC;
  float lslimAdjust;
  bool lslimDebug;

} LSLIMITER;

#endif /* LSLIMITER_PRIV_H */
