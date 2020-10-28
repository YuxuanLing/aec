#ifndef LSEXCHNDL_H
#define LSEXCHNDL_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flåten Aarnes (IFA), Tandberg Telecom AS.
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

typedef struct LSEXCHNDL * LSEXCHNDL_PTR;


/***************************************************************************
 LSEXCHNDL_CREATE
 *   Description: Allocates memory for struct LSEXCHNDL.
 *
 *   Parameters:   -
 *
 *   Returns:      pointer to struct LSEXCHNDL.
 *
 *   Note:         -
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
LSEXCHNDL_PTR lsexchndl_create(void);

/***************************************************************************
 * LSEXCHNDL_DESTROY
 *   Description:  Frees allocated memory pointed by input parameter.
 *
 *   Parameters:   pLsexchndl
 *
 *   Returns:      -
 *
 *   Note:         -
 *
 *   Globals used: -
 *
 *   Restrictions: -
 ***************************************************************************/
void lsexchndl_destroy(LSEXCHNDL_PTR pLsexchndl);

/***************************************************************************
 * LSEXCHNDL_INIT
 *  Description:  Initialises members in struct lsexchndl.
 *
 *  Parameters:   pLsexchndl
 *
 *  Returns:      -
 *
 *  Note:         -
 *
 *  Globals used: -
 *
 *  Restrictions: Memory pointed by pLsexchndl shall be allocated by
 *                calling function micexchndl_create
 ***************************************************************************/
void lsexchndl_init(LSEXCHNDL_PTR pLsexchndl);

/***************************************************************************
 * LSEXCHNDL_PROCESS
 *   Description:   Main lsexchndl routine
 *
 *   Parameters:    pLsexchndl     - pointer to LSEXCHNDL struct
 *                  lsbuf          - loudspeaker output from analyse filter
 *                  loudsgain      -
 *
 *   Returns:       -
 *
 *   Note:          -
 *
 *   Globals used:  -
 *
 *   Restrictions:
 ***************************************************************************/
void lsexchndl_process(LSEXCHNDL_PTR pLsexchndl,
		       float * lsbuf,
		       float loudsGain);


/************************************************************************************
 * LSEXCHNDL_DLPOW
 *   Auth.: IFA
 *
 *   Desc.: delay line power calculation routine. Calculates an approximation
 *          of the power of lsfft, based on the most significant subbands for speech.
 *          subbuf contain the "EXCSUBUSED" bands in lsfft
 *
 *   Parameters:    pLsexchndl     - pointer to LSEXCHNDL struct
 *                  subuf          -
 ***********************************************************************************/
void lsexchndl_dlpow(LSEXCHNDL_PTR pLsexchndl, float * subbuf);

/************************************************************************************
 * LSEXCHNDL_SHIFT
 *   Auth.: IFA
 *
 *   Desc.: Delayline shift
 *          Writes new data onto delayline. The write-index indicates where
 *          present starting-point is.
 *
 *   Parameters:    pLsexchndl     - pointer to LSEXCHNDL struct
 *                  subuf          -
  ************************************************************************************/
void lsexchndl_shift(LSEXCHNDL_PTR pLsexchndl, float * subbuf);


/* lsexchndl_load **********************************************************
 * Auth.: JPS
 * Desc.: loads pLsexchndl struct into iram
 * Returns: ptr to pLsexchndl struct in iram
 ************************************************************************/
LSEXCHNDL_PTR lsexchndl_load(LSEXCHNDL_PTR pExchndlIn);


/* lsexchndl_flush **********************************************************
 * Auth.: JPS
 * Desc.: flush pLsexchndl struct to sdram
 ************************************************************************/
void lsexchndl_flush(LSEXCHNDL_PTR pLsexchndlSdram, LSEXCHNDL_PTR pLsexchndlIram);




#ifdef UNITTEST
/*************************************************************************
 * UNITTEST_LSEXCHNDL
 *   Auth.: Oystein Birkenes (OBI)
 *   Desc.: Test of lsexchndl.
 ************************************************************************/
void unittest_lsexchndl(void);
#endif /* UNITTEST */

#ifdef __cplusplus
}
#endif

#endif /* LSEXCHNDL_H */
