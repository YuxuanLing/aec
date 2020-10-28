#ifndef DEREVERB_PRIV_H
#define DEREVERB_PRIV_H
/***************************************************************************
*                              A U D I O
*--------------------------------------------------------------------------
*                  (C) Copyright Tandberg Telecom AS 2009
*==========================================================================
*
* Author        : Gisle Enstad (GEN), Tandberg Telecom AS.
* Co-author     : -
*
**************************************************************************/

#include "audec_defs.h"
#include <stdbool.h>

#define DEREVERB_DELAYLINE 6       

typedef struct DEREVERB {

    /* Dereverb buffers/variables */
    float siglevDelayline[DEREVERB_DELAYLINE][SUBUSED];
    int delaylineIndex;
    int debug;

} DEREVERB;
#endif /* DEREVERB_PRIV_H */
