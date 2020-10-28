/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flåten Aarnes (IFA), Tandberg Telecom AS.
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
 **************************************************************************/

#include "lsexchndl.h"
#include "lsexchndl_priv.h"
#include <stdlib.h>
#include <stdio.h>

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
LSEXCHNDL_PTR lsexchndl_create(void)
{
    LSEXCHNDL_PTR pLsexchndl;
    pLsexchndl = (LSEXCHNDL_PTR) malloc(sizeof(LSEXCHNDL));

    if (pLsexchndl == NULL)
    {
        fprintf(stderr, "lsexchndl_create: Could not allocate lsexchndl buffer");
    }

    return pLsexchndl;
}

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
void lsexchndl_destroy(LSEXCHNDL_PTR pLsexchndl)
{
    free(pLsexchndl);
}


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
 *                calling function lsexchndl_create
 ***************************************************************************/
void lsexchndl_init(LSEXCHNDL_PTR pLsexchndl)
{
    int i, j;
    pLsexchndl->dlwriteix = 0; /* where to write present data in excdline (which row) */
    pLsexchndl->lsLastIx = (EXCLAGCNT - 1); /* where to read data in exclsls */

    for (i = 0; i < EXCLAGCNT; i++)
    {
        pLsexchndl->exclsls[i] = (float)0.0;
    }

    for (i = 0; i < EXCLAGCNT; i++)
    {
        for (j = 0; j < EXCSUBUSED * 2; j++)
        {
            pLsexchndl->excdline[i][j] = (float)0.0;
        }
    }
}

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
void lsexchndl_process(LSEXCHNDL_PTR pLsexchndl, float * lsbuf, float loudsGain)
{
    int i;
    float lsexsubbuf[EXCSUBUSED * 2];

    for (i = 0; i < EXCSUBUSED * 2; i += 2)
    {
        lsexsubbuf[i]     = lsbuf[(SUBUSED_START + EXCSUBUSED_START) * 2 + i]     * loudsGain; /* Re */
        lsexsubbuf[i + 1] = lsbuf[(SUBUSED_START + EXCSUBUSED_START) * 2 + i + 1] * loudsGain; /* Im */
    }

    lsexchndl_dlpow(pLsexchndl, lsexsubbuf); /* lsbuf's first "EXCSUBUSED" band*/
    lsexchndl_shift(pLsexchndl, lsexsubbuf); /* lsbuf's first "EXCSUBUSED" band*/

    //  lsexchndl_dlpow(pLsexchndl , &lsbuf[SUBUSED_START*2]); /* lsbuf's first "SUBUSED" band*/
    //  lsexchndl_shift(pLsexchndl , &lsbuf[SUBUSED_START*2]); /* lsbuf's first "SUBUSED" band*/
}


/************************************************************************************
 * LSEXCHNDL_DLPOW
 *   Auth.: IFA
 *
 *   Desc.: delay line power calculation routine. Calculates an approximation
 *          of the power of lsfft, based on the most significant subbands for speech.
 *          subbuf contains the "EXCSUBUSED" bands in lsfft
 *
 *   Parameters:    pLsexchndl     - pointer to LSEXCHNDL struct
 *                  subbuf         -
 ***********************************************************************************/
void lsexchndl_dlpow(LSEXCHNDL_PTR pLsexchndl, float * subbuf)
{
    int i;
    int j;
    float tmp = 0.0f;
#ifdef __INTEL_COMPILER
#pragma novector
#endif
    for (i = 0; i < (EXCSUBUSED) * 2; i++)
    { /* squared l2-norm of exch-subbands of lsfft */
        tmp += subbuf[i] * subbuf[i]; // sum(buf(params.excsubix).*conj(buf(params.excsubix)))
    }

    i = pLsexchndl->lsLastIx;
    j = (pLsexchndl->lsLastIx + 1) % EXCLAGCNT; // rotates modulo EXCLAGCNT
    pLsexchndl->exclsls[j] = ((1.0f - EXCPSPD) * pLsexchndl->exclsls[i]) + (EXCPSPD * tmp);

    pLsexchndl->lsLastIx = (short)j;
    /* lsLastIx now points to where the next element of exclsls should
       be written. (rotating modulo EXCLAGCNT)) */
}

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
void lsexchndl_shift(LSEXCHNDL_PTR pLsexchndl, float *  subbuf)
{
    int i;
    int j = pLsexchndl->dlwriteix;

    /* *subbuf = lsfft[2*SUBUSED_START+2*EXCSUBUSED_START] */
#ifdef __INTEL_COMPILER
#pragma novector
#endif
    for (i = 0; i < EXCSUBUSED * 2; i++)
    {
        pLsexchndl->excdline[j][i] = * subbuf;
        subbuf++;
    }
    ((pLsexchndl->dlwriteix == (EXCLAGCNT - 1)) ? (pLsexchndl->dlwriteix = 0) : (pLsexchndl->dlwriteix++));
    /* dlwriteix now points to where the next row of excdline should
       be written. (rotating modulo EXCLAGCNT)) */
}



/* lsexchndl_load **********************************************************
 * Auth.: JPS
 * Desc.: loads pLsexchndl struct into iram
 * Returns: ptr to pLsexchndl struct in iram
 ************************************************************************/
LSEXCHNDL_PTR lsexchndl_load(LSEXCHNDL_PTR pLsexchndlIn)
{
    /* MEMXFERMUSTDIE */
    return pLsexchndlIn;
}

/* lsexchndl_flush **********************************************************
 * Auth.: JPS
 * Desc.: flush pLsexchndl struct to sdram
 ************************************************************************/
void lsexchndl_flush(LSEXCHNDL_PTR pLsexchndlSdram, LSEXCHNDL_PTR pLsexchndlIram)
{
    /* MEMXFERMUSTDIE */
    (void)pLsexchndlSdram;
    (void)pLsexchndlIram;
}




#ifdef UNITTEST
/*************************************************************************
 * UNITTEST_LSEXCHNDL
 *   Auth.: Oystein Birkenes (OBI)
 *   Desc.: Test of lsexchndl.
 ************************************************************************/

#include "unittest.h"
#include <math.h>

static float TestInput[EXCLAGCNT][(SUBUSED_START + EXCSUBUSED_START + EXCSUBUSED) * 2];

/*************************************************************************
 * UNITTEST_LSEXCHNDL
 *************************************************************************/
void unittest_lsexchndl()
{
    float loudsGain = 1.0;
    int i, j;
    int subused_start = (SUBUSED_START + EXCSUBUSED_START) * 2; /* The index of the first subband used */
    float exclsls[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float ip[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float in_re, in_im, out_re, out_im;
    float epsilon = 1e-7;
    int differ = 0;
    int numFrames = EXCLAGCNT;
    LSEXCHNDL_PTR pLsexchndl;

    /* Allocate memory */
    pLsexchndl = lsexchndl_create();

    /* Generate test input (0.001, 0.002, 0.003,...) */
    for (i = 0; i < numFrames; i++)
    {
        for (j = 0; j < subused_start + EXCSUBUSED * 2; j++)
        {
            TestInput[i][j] = 0.001 * (i * (subused_start + EXCSUBUSED * 2) + j);
        }
    }

    /* Run lsexchndl */
    lsexchndl_init(pLsexchndl);
    for (i = 0; i < numFrames; i++)
    {
        lsexchndl_process(pLsexchndl, TestInput[i], loudsGain);
    }

    /* Generate inner products of each frame with itself */
    for (i = 0; i < numFrames; i++)
    {
        for (j = 0; j < EXCSUBUSED * 2; j++)
        {
            ip[i] += TestInput[i][subused_start + j] * TestInput[i][subused_start + j];
        }
    }

    /* Compute power */
    exclsls[0] = EXCPSPD * ip[0];
    for (i = 1; i < numFrames; i++)
    {
        exclsls[i] = (1.0f - EXCPSPD) * exclsls[i - 1] + EXCPSPD * ip[i];
    }

    /* Check if pLsexchnld->exclsls is (almost) the same as exclsls */
    for (i = 0; i < EXCLAGCNT; i++)
    {
        if (fabs(exclsls[i] - pLsexchndl->exclsls[i]) > epsilon)
        {
            differ = 1;
            break;
        }
    }


    /* Check if delaylines is correct (out == in) */
    for (i = 0; i < EXCLAGCNT && differ == 0; i++)
    {
        for (j = 0; j < EXCSUBUSED * 2 && differ == 0; j += 2)
        {
            out_re = pLsexchndl->excdline[i][j];
            out_im = pLsexchndl->excdline[i][j + 1];
            in_re = TestInput[i][subused_start + j];
            in_im = TestInput[i][subused_start + j + 1];

            differ = (out_re != in_re) || (out_im != in_im);
        }
    }

    lsexchndl_destroy(pLsexchndl);

    unittest_context("lsexchndl");
    unittest_assert(differ == 0);
}

#endif /* UNITTEST */
