#ifndef DELAY_ESTIMATION_H
#define DELAY_ESTIMATION_H

#include "audec_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>

struct DELAY_ESTIMATION * delayEstimation_create(void);
void delayEstimation_destroy(struct DELAY_ESTIMATION * pDelayEstimation);
void delayEstimation_init(struct DELAY_ESTIMATION * pDelayEstimation);
int delayEstimation_estimateDelay(struct DELAY_ESTIMATION * delayEstimation, const COMPLEX32 *  micFft, COMPLEX32 *  lsFft);
void delayEstimation_printStatus(struct DELAY_ESTIMATION * delayEstimation);
void delayEstimation_setDebug(struct DELAY_ESTIMATION * delayEstimation, int value);
void delayEstimation_setWindowsDelay(struct DELAY_ESTIMATION * delayEstimation, float value);
void delayEstimation_setDelayingOn(struct DELAY_ESTIMATION * delayEstimation, int onoff);
void delayEstimation_returnIntermediateLsBuf(const struct DELAY_ESTIMATION * delayEstimation, COMPLEX32 *  lsFft, int changedDelay, int delay_position);
void delayEstimation_returnHistoryLsSample(const struct DELAY_ESTIMATION * delayEstimation, 
                                           COMPLEX32 *  lsComplexSample, 
                                           int movedDelay, 
                                           int filtlen, 
                                           int subband,
                                           int delay_position);

struct INDEXED_MAX_SUM
{
    int   index;
    float max_sum;
};

 void delayEstimation_addDataToLsLine(struct DELAY_ESTIMATION * delayEstimation, const COMPLEX32 *  lsFft);
 void delayEstimation_crossCorrelate(struct DELAY_ESTIMATION * delayEstimation, const COMPLEX32 *  micFft);
 struct INDEXED_MAX_SUM delayEstimation_findMaxSum(struct DELAY_ESTIMATION * delayEstimation);
 bool delayEstimation_evaluateMaxCrossCorrelation(const struct DELAY_ESTIMATION * delayEstimation, const struct INDEXED_MAX_SUM * maxCorrelationValue);
 int delayEstimation_decideUponDelayEstimate(struct DELAY_ESTIMATION * delayEstimation, int estimatedDelay);
 void delayEstimation_delayLsBuf(const struct DELAY_ESTIMATION * delayEstimation, COMPLEX32 * lsFft);
 void delayEstimation_measurePower(float * power, const COMPLEX32 * lsbuf);



#ifdef __cplusplus
}
#endif

#endif /* DELAY_ESTIMATION_H */
