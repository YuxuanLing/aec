#ifndef LSPROCESS_H
#define LSPROCESS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "audtypes.h"

typedef struct LSPROCESS * LSPROCESS_PTR;

LSPROCESS_PTR lsprocess_create(FILTERBANK_USE_TYPE filterbankUse, int channels);
void lsprocess_destroy(LSPROCESS_PTR pLs, FILTERBANK_USE_TYPE filterbankUse);
void lsprocess_init(LSPROCESS_PTR pLs, FILTERBANK_USE_TYPE filterbankUse);
void lsprocess_load(LSPROCESS_PTR, float * lsBuf);
void lsprocess_process(LSPROCESS_PTR pLs, float *afir);
void lsprocess_setlsGainAdjust(LSPROCESS_PTR pLs, bool onoff);
void lsprocess_status(LSPROCESS_PTR pLs);
void lsprocess_setloudsGain(LSPROCESS_PTR pLs, float value);
void lsprocess_setDebug(LSPROCESS_PTR pLs, int value);
void lsprocess_shift(LSPROCESS_PTR pLs);

#ifdef __cplusplus
}
#endif

#endif /* LSPROCESS_H */
