#ifndef EC_H
#define EC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "keyclickremoval.h"
#include "audtypes.h"
#include "lsprocess.h"

typedef struct EC * EC_PTR;

EC_PTR ec_create(FILTERBANK_USE_TYPE filterbankUse, int channels);
void ec_destroy(EC_PTR pEc, FILTERBANK_USE_TYPE filterbankUse);
void ec_init(EC_PTR pEc, FILTERBANK_USE_TYPE filterbankUse, LSPROCESS_PTR pLs);

/* EC_LOAD
 * Starts dmax transfers on analyse and lsexchndl data.
 * Requires a later call to memxfer_wait with buffer-pointers
 *  to ensure successful transfers.
 */
void ec_load(EC_PTR pEc);
/* EC_PROCESS
 * Filterbank and Echo cancellator
 */
void ec_process(EC_PTR pEc,
                LSPROCESS_PTR pLs,
                float * micBuf,
                float * outBuf,
                bool runLsprocess);

/* ec_getAerl */
float ec_getAerl(EC_PTR pEc, LSPROCESS_PTR pLs,  float miclevel);
void ec_updateAfirIndex(EC_PTR pEc, int index);
int ec_getAfirIndex(EC_PTR pEc);
void ec_keyclickEvent(EC_PTR pEc, int *inbuf, int sampleindex);
bool ec_getkeyclickState(EC_PTR pEc);
void ec_keyclickremovalClip(EC_PTR pEc, bool onoff);
void ec_setNoisereduction(EC_PTR pEc, bool onoff);
void ec_setWindowsDelay(EC_PTR pEc, float value);
void ec_setDelayEstimationOnOff(EC_PTR pEc, int onoff);

/* Test functions */
void ec_setlsGainAdjust(EC_PTR pEc, bool onoff);
void ec_setCalcNewDelta(EC_PTR pEc, bool onoff);
void ec_setMinDelta(EC_PTR pEc, float value);
void ec_setMaxDelta(EC_PTR pEc, float value);
void ec_setFilter(EC_PTR pEc, int value);
void ec_setFilterZero(EC_PTR pEc, int value);
void ec_setNlp(EC_PTR pEc, bool onoff);
void ec_setDereverb(EC_PTR pEc, bool onoff);
void ec_setComNoise(EC_PTR pEc, bool onoff);
void ec_setSubbandNoiseRed(EC_PTR pEc, bool onoff);
void ec_setNlpSub(EC_PTR pEc, bool onoff);
void ec_setNlpHi(EC_PTR pEc, bool onoff);
void ec_setNlpExtragain(EC_PTR pEc, bool onoff);
void ec_setNlpfullg(EC_PTR pEc, bool onoff);
void ec_setNlpfullgHi(EC_PTR pEc, bool onoff);
void ec_setNlpSubusedEndHi(EC_PTR pEc, float value);
void ec_setNlpTransition(EC_PTR pEc, float value);
void ec_setMicexchndl(EC_PTR pEc, bool onoff);
void ec_setDebug(EC_PTR pEc, int value);
void ec_status(EC_PTR pEc);
void ec_setFbPrototypefilt(EC_PTR pEc, LSPROCESS_PTR pLs, int type, FILTERBANK_USE_TYPE filterbankUse);
void ec_setKeyclickprocess(EC_PTR pEc, bool onoff);
void ec_setKeyclickmute(EC_PTR pEc, bool onoff);
void ec_setKeyclicklength(EC_PTR pEc, int value);
void ec_setKeyclickgi(EC_PTR pEc, float value);
void ec_setKeyclickgd(EC_PTR pEc, float value);
void ec_setKeyclickmaxmutecount(EC_PTR pEc, int value);
void ec_setKeyclickmaxlimitcount(EC_PTR pEc, int value);
void ec_setComNoiseAmp(EC_PTR pEc, float value);
void ec_debug12print(float * pLsBuf, float loudsGain);
void ec_setNoiseredHpFilt(EC_PTR pEc, bool onoff);
void ec_setShellOn(EC_PTR pEc, bool onoff);
void ec_setShellLsLevAdjust(EC_PTR pEc, float value);
void ec_setShellnoisgnAdjust(EC_PTR pEc, float value);
void ec_setShellnlpgAdjust(EC_PTR pEc, float value);
void ec_setShellvarIncspeed(EC_PTR pEc, float value);
void ec_setShellvarDecspeed(EC_PTR pEc, float value);
void ec_setStereoDelta(EC_PTR pEc, float value);
void ec_printNumChannels(EC_PTR pEc);

#ifdef UNITTEST
void run_unittest_ec(const char fileMic[], const char fileNetworkIn[],
                     const char fileNetworkOutReference[], const char fileNetworkOut[],
                     int numFrames);
void unittest_ec(void);
#endif /* UNITTEST */

#ifdef __cplusplus
}
#endif

#endif /* EC_H */
