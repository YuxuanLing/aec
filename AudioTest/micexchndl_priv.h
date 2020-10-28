#ifndef MICEXCHNDL_PRIV_H
#define MICEXCHNDL_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 * Co-author     : Geir Ole Øverby (GEO), Tandberg Telecom AS
 *
 * Switches      : <COMPILER SWITCH HERE>
 *		     <add description here>
 *
 * Description   : <add description here>
 *
 * Note		 : <notes here>
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *		   Modification made.
 *
 * yyyy-mm-dd    : <signature>
 *		   File created.
 *
 **************************************************************************/

#include "audec_defs.h"
#include "audtypes.h"

#include <stdbool.h>

/* -------------------------------------------------*/
/* defines 'struct micexchndl' which contains      */
/* variables etc. that define the exceptionhandling */
/* -------------------------------------------------*/
typedef struct MICEXCHNDL {

    int   channels;
    float excgain;   /* level exceptionhandler gain (flat) */
    float excbefbef; /* */
    float excaftaft; /* */
    float excgainu;  /* */
    float excgainf;  /* */
    short dlreadix ; /* where to read last data from ls-exchandler's excdline (which row) */
    float excaftls[MAX_NUM_CHANNELS][EXCLAGCNT][EXCSUBUSED*2];//   = zeros(params.exclagcnt,params.excsubcnt);  /* */
    bool  ExptHndlOn;

#if 0 //FIR implentation
    float excbefbefAtt; /* bef power short*/
    float excaftaftAtt; /* aft power short*/
    float excbefbefDelay;
    float excbefls[EXCLAGCNT][EXC_DELAY_SUBUSED*2]; // = zeros(params.exclagcnt,params.excsubcnt);  /* */
    float excaftlsAtt[EXC_FIR_LAG][EXCSUBUSED*2];   // = zeros(params.exclagcnt,params.excsubcnt);  /*  REDUCE TO */
    float excdlinebef[EXCLAGCNT*2][EXCSUBUSED*2];   // = zeros(params.exclagcnt,params.excsubcnt);  /* */
    float excdlineaft[EXCLAGCNT*2][EXCSUBUSED*2];   // = zeros(params.exclagcnt,params.excsubcnt);  /* */
    int dlStartIdx;
    int systemDelay;
    float maxDelayCorr;
    int maxDelayCnt;
    int prevDelayIdx;
    int currentDelayIdx;
    float excSNR;
    float excNoiLev[EXCSUBUSED];
    float excPowerLong[EXC_DELAYCALC];
    float excNoiseLevel;
    float excEchoLevel;

#endif

} MICEXCHNDL;

#endif /* MICEXCHNDL_PRIV_H */
