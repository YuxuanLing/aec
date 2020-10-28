#ifndef SHELL_PRIV_H
#define SHELL_PRIV_H

#include "audec_defs.h"
#include <stdbool.h>

#define SHELL_TOTALDELAY (1536) // asymmetric filterbank has delay: 32ms * 48000; (symmetric filterbank has delay = afirlen/2 + sfirlen/2 = 3840)
#define SUBUSED16 116   /* used for level estimate without high-frequencies */

typedef struct SHELL
{
    float delayline[SHELL_TOTALDELAY];
    float gain;
    float prevgain;
    float variable_gain;
    int   writeIx;
    int   outIx;
    int   rampIx;
    bool  shellOn;
    float lslevAdjust;
    float noisegnAdjust;
    float nlpgAdjust;
    float varIncSpeed;
    float varDecSpeed;
    float lslevel[MAX_NUM_CHANNELS][SUBUSED16];
    int   debug;
} SHELL;

#endif

