#ifndef SHELL_H
#define SHELL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "noisereduction.h"
#include "lsprocess.h"
#include <stdbool.h>

typedef struct SHELL * SHELL_PTR;

SHELL_PTR shell_create(void);
void shell_destroy(SHELL_PTR pShell);
void shell_init(SHELL_PTR pShell, int channels);
float shell_process(SHELL_PTR pShell, NOIRED_PTR pNoired, LSPROCESS_PTR pLs);
SHELL_PTR shell_load(SHELL_PTR pShellSdram);
void shell_flush(SHELL_PTR pShellSdram, SHELL_PTR pShellIram);
void shell_loadDelayline(SHELL_PTR pShell);
void shell_flushDelayline(SHELL_PTR pShell);
void shell_input(SHELL_PTR pShell, float *inputbuf);
void shell_output(SHELL_PTR pShell, float *outBuf);
void shell_setlslevAdjust(SHELL_PTR pShell, float value);
void shell_setnoisegnAdjust(SHELL_PTR pShell, float value);
void shell_setnlpgAdjust(SHELL_PTR pShell, float value);
void shell_setvarIncSpeed(SHELL_PTR pShell, float value);
void shell_setvarDecSpeed(SHELL_PTR pShell, float value);
void shell_setShellOn(SHELL_PTR pShell, bool onoff);
void shell_status(SHELL_PTR pShell);
void shell_setDebug(SHELL_PTR pShell, int value);

#ifdef UNITTEST
void unittest_shell(void);
#endif

#ifdef __cplusplus
}
#endif

#endif /* SHELL_H */

