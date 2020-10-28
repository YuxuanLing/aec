#ifndef DEESSER_H
#define DEESSER_H

#include "audec_defs.h"

//#define DEESSER_DEBUG

#define INIT_BW_THRESHOLD 0.3f
#define INIT_S_THRESHOLD  1.05f
#define INIT_GAIN         0.35f
#define INIT_FADE_STEPS   12
#define INIT_LOGGING      0
#define INIT_DEESSER_ON   0
#define INIT_SHOW_SHARPNESS 0

#define CRITICAL_BANDS    48
#define FFT_WIDTH         62.5f
#define LOUDNESS_CONST    0.08f
#define CRITICAL_BAND_PEAK_MINIMUM 28

typedef struct DEESSER {
    float bwThreshold;
    float sharpnessThreshold;
    float gain;
    int   fadeSteps;
    int   logging;
    int   deEsserOn;
    int   showSharpness;
} DEESSER;

typedef struct DEESSER * DEESSER_PTR;

void  deesser_calcCriticalBandIntensities(COMPLEX32 * micfft);
float deesser_calcSharpness(void);
void  deesser_reduceSibilants(COMPLEX32 * micfft);
void  deesser_process(COMPLEX32 * micfft);
void  deesser_status(void);
DEESSER_PTR deesser_getDeEsser(void);


#endif //DEESSER_H
