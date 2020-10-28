#ifndef ANALYSE_PRIV_H
#define ANALYSE_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes, Tandberg Telecom AS
 * Co-author     : Johan Malmstigen, S&T ES
 *
 * Switches      : -
 *
 * Description   : Private data for the analyse module
 *
 * Note		 : -
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *		   Modification made.
 *
 * 2004-02-03    : JMA
 *	           Created from eccanc_common.h and analyse_priv.h from SOLO.
 **************************************************************************/

/*-------------------------------------------------------------------------
 * -- INCLUDE FILES AREA
 *-------------------------------------------------------------------------*/
#include "audec_defs.h"

/*-------------------------------------------------------------------------
 * -- PREPROCESSOR CONSTANTS AREA
 *-------------------------------------------------------------------------*/
/* see audec_defs.h */

/*-------------------------------------------------------------------------
 * -- PRIVATE MACROS AREA
 *-------------------------------------------------------------------------*/
//#define ANALYSE_PRIV_GET_OUTPUT(analyse_ptr)	((analyse_ptr)->fftout)

/*-------------------------------------------------------------------------
 * -- PRIVATE TYPES AREA
 *-------------------------------------------------------------------------*/
/* -------------------------------------------------*/
/* defines 'struct ANALYSE' which contains          */
/* variables etc. that define the analysis filter.  */
/* These should be filled by an init routine.       */
/* -------------------------------------------------*/
typedef struct ANALYSE {
    int   inbufStart;        /* where new data is stored in inbuf */
    float prefiltState;
    float *inbuf;          /* sound samples in (aFirLen real numbers)   */
    float *inbufSDRAM;     /* SDRAM - sound samples in (aFirLen real numbers)   */
    COMPLEX32 *fftout;     /* frequency samples out 384 (768) complex numbers) */
    const float *afir;
    int   filtertype;
    float prefiltercoeff;  /* b[1] in prefilter. b[0] assumed to be 1 */
    float prefiltergain;
} ANALYSE;



#endif /* ANALYSE_PRIV_H */
