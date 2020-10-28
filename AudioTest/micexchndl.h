#ifndef MICEXCHNDL_H
#define MICEXCHNDL_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes (IFA), Tandberg Telecom AS.
 * Co-author     : Geir Ole Øverby (GEO), Tandberg Telecom AS
 *		   Johan Malmstigen (JMA), S&T ES
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
 **************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#include "lsexchndl.h"
#include <stdbool.h>

typedef struct MICEXCHNDL * MICEXCHNDL_PTR;

MICEXCHNDL_PTR micexchndl_create(void);
void micexchndl_destroy(MICEXCHNDL_PTR pMicexchndl);
void micexchndl_init(MICEXCHNDL_PTR pMicexchndl,
                     int channels);
void micexchndl_process(LSEXCHNDL_PTR pLsexchndl[],
			MICEXCHNDL_PTR pMicexchndl,
			float *  micfft,
			float *  estecho
			);

/* micexchndl_load **********************************************************
 * Auth.: JPS
 * Desc.: loads pMicexchndl struct into iram
 * Returns: ptr to pMicexchndl struct in iram
 ************************************************************************/
MICEXCHNDL_PTR micexchndl_load(MICEXCHNDL_PTR pMicexchndlIn);

/* micexchndl_flush **********************************************************
 * Auth.: JPS
 * Desc.: flush pMicexchndl struct to sdram
 ************************************************************************/
void micexchndl_flush(MICEXCHNDL_PTR pMicexchndlSdram, MICEXCHNDL_PTR pMicexchndlIram);

/* micexchndl_getExcgain ************************************************
 * Auth.: JPS
 * Desc.: returns pMicexchndl->excgain
 ************************************************************************/
float micexchndl_getExcgain(MICEXCHNDL_PTR pMicexchndl);

void micexchndl_set(MICEXCHNDL_PTR pMicexchndl, bool onoff);
void micexchndl_statusMicexchndl(MICEXCHNDL_PTR pMicexchndl);



#ifdef UNITTEST
/* unittest_micexchndl **************************************************** */
/* Auth.: Jens Petter Stang (JPS)                                         */
/* Desc.: simple run-trhough of micexchndl with 0.3sec random input        */
/*        and check against testvector verified in matlab                 */
/* ********************************************************************** */
void unittest_micexchndl(void);
#endif /*UNITTEST*/

#ifdef __cplusplus
}
#endif

#endif /* MICEXCHNDL_H */
