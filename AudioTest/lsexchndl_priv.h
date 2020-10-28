#ifndef LSEXCHNDL_PRIV_H
#define LSEXCHNDL_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 **************************************************************************/

#include "audec_defs.h"

#define EXC_MAXLEN_LX (EXCLAGCNT + EXC_DELAYCALC + EXC_LAGSHIFT)
/* -------------------------------------------------*/
/* defines 'struct lsexchndl' which contains       */
/* variables etc. that define the exceptionhandling */
/* -------------------------------------------------*/
typedef struct LSEXCHNDL {
    short dlwriteix; /* where to write present data in excdline (which row) */
    short lsLastIx;  /* where to write present data in exclsls */
    float exclsls[EXCLAGCNT];
    float excdline[EXCLAGCNT][EXCSUBUSED * 2];
#if 0 /*FIR implementation */
    int   dlStartIdx;
    float exclslsDelay[EXC_DELAYCALC];
    float exclslsDelayTmp[EXCLAGCNT + EXC_DELAYCALC + EXC_LAGSHIFT];    
    float excdline[(EXC_MAXLEN_LX) * 2][EXCSUBUSED * 2];
    float exclslsTmp[EXCLAGCNT + EXC_FIR_LAG + EXC_LAGSHIFT];
    float exclslsAtt[EXC_FIR_LAG]; /*lspower for short xcorr calculation: attack*/
    float exclslsTmpAtt[EXCLAGCNT_ATTACK + EXC_FIR_LAG + EXC_LAGSHIFT];
#endif
} LSEXCHNDL;


#endif /* LSEXCHNDL_PRIV_H */
