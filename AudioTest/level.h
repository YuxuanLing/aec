#ifndef _LEVEL_H
#define _LEVEL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//#if !defined(__INTEL_COMPILER) && !defined(__GNUC__)
//#define 
//#endif

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct LEVELCONTROL * LEVELCONTROL_PTR;

/**
 * Create levectrl instance
 *
 * @retval : pointer to levelctrl struct
 */
LEVELCONTROL_PTR levelctrl_create(void);

/**
 * Destroy levectrl instance
 */
void levelctrl_destroy(LEVELCONTROL_PTR levelcontrol_ptr);

/**
 * Reset levectrl instance
 */
void levelctrl_reset(LEVELCONTROL_PTR levelCtrl, int buffersize);

/**
 * Remove DC, ramp and find level. If DC is already removed, call ramp instead
 */
void levelctrl_rampAndDcRemove(LEVELCONTROL_PTR  levelCtrl,
                               int nChan,
                               float * AbsLevel,
                               uint32_t * DcRemoveDone,
                               float * Gain,
                               float *pData);

/**
 * Calculate abslevel and remove DC. If DC is already removed, update abslevel
 */
void levelctrl_findLevelAndDcRemove(LEVELCONTROL_PTR  levelCtrl,
                                    int nChan,
                                    float * AbsLevel,
                                    uint32_t * DcRemoveDone,
                                    float *pData,
                                    bool consecutiveData,
                                    float *pData_ch2);

/**
 * Get signal to noise ratio (based on absolute value calculations)
 */
float levelctrl_getSNR(LEVELCONTROL_PTR  levelCtrl);

#if 0
/**
 * Add two streams together and ramp the gain to be smooth. Call
 * levelctrl_findLevelAndDcRemove first if DC not has been removed already
 *
 * NOTE: If this function is to be used again, please avoid using
 *       AUDINT_STREAM_DATA_STRUCT, since this struct causes a dependency
 *       to an FSM-header.
 */
void levelctrl_rampAndAdd(LEVELCONTROL_PTR  levelCtrl1,
                          LEVELCONTROL_PTR  levelCtrl2,
                          AUDINT_STREAM_DATA_STRUCT * pAudBuf1,
                          AUDINT_STREAM_DATA_STRUCT * pAudBuf2);
#endif

float findAbsLevel(float *pData,
                   int buflen);

/**
 * Ramp the gain to be nice and smooth. Will call levelctrl_findLevelAndDcRemove
 * on first if DC not has been removed allready.
 */
void levelctrl_ramp(LEVELCONTROL_PTR  levelCtrl,
                    int nChan,
                    float * AbsLevel,
                    uint32_t * DcRemoveDone,
                    float * Gain,
                    float *pData);
/**
 * Get absolute level
 */
float levelctrl_getLevel(LEVELCONTROL_PTR levelCtrl);

/**
 * Ramp buffer
 */
void ramp(float *pSrc,  /* Must be 64-bits aligned */
          float oldGain,
          float newGain,
          int   buflen, /* Must be a multiple of 8 */
          int   nChan);

float levelctrl_getSignalLevel(LEVELCONTROL_PTR  levelCtrl);

float levelctrl_getNoiseLevel(LEVELCONTROL_PTR  levelCtrl);

#ifdef __cplusplus
}
#endif

#endif
