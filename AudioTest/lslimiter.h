#ifndef LSLIMITER_H
#define LSLIMITER_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Ingvar Flaten Aarnes, Tandberg Telecom AS.
 * Co-author     : Johan Malmstigen, S&T ES
 *
 * Switches      : -
 *
 * Description   : Part of echocanceller. By estimating the aerl-level and
 *                 signal level from farend+vcr etc. destined for loudspeaker,
 *                 the loudspeaker signal is adjusted to avoid saturation of
 *                 the (local) mic signal. (If the micsignal is distorted,
 *                 the echocanceller filters will not converge.
 *
 * Known issues:   If the 1/aerl is small and mic is heavily distorted
 *                 by a signal with more or less constant high level, the
 *                 echo-filters will not converge, and lslimiter will not
 *                 attenuate loudspeaker. (This is a rare situation).
 *
 * Note    : <notes here>
 *
 * Documentation : <Xxxxxx-Document.doc>
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *       Modification made.
 *****************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define LSLIM_DELAY 240 /* delay */

typedef struct LSLIMITER * LSLIMITER_PTR;

/*****************************************************************************
 * Author  : IFA, JMA
 * Description   : Create routine, called once at startup
 *
 * Returns       : pointer to an lslimiter instance, NULL if no free memory
 *       left
 ****************************************************************************/
LSLIMITER_PTR lslimiter_create(void);

/*****************************************************************************
 * Author        : JMA
 * Description   : Destory routine, called once at startup
 *
 * Returns       : pointer to an lslimiter instance, NULL if no free memory
 *       left
 ****************************************************************************/
void lslimiter_destroy(LSLIMITER_PTR pLslimiter);

/*****************************************************************************
 * Author        : IFA, JMA
 * Description   : Init routine, called once at startup. Initializes the
 *       lslimiter structure pointed by pLslimiter.
 *
 * Returns       : -
 ****************************************************************************/
void lslimiter_init(LSLIMITER_PTR pLslimiter, int numAcousticInputs);

/*****************************************************************************
 * Author        : OGA, IFA, JMA
 * Description   : If loudspeaker signal is too loud, this routine tries to
 *       attenuate the signal.
 *
 *             Works only for 10ms buffers
 * Returns       : -
 ****************************************************************************/
void lslimiter_process(LSLIMITER_PTR pLslimiter,
                       float *inbuf[],
                       float *outbuf[],
                       float *delaylines[],
                       int nchannels,
                       float aerlInv_current,
                       float * limgain
                      );

/* STEREO COLLAPSING *********************************************************/
/* Auth.: GEO                                                                */
/*****************************************************************************/
void lslimiter_collapse_stereo_detect(float bef_lvl,
                                      float aft_lvl,
                                      float est_lvl,
                                      float noi_lvl,
                                      float stQty,
                                      LSLIMITER_PTR pLslimiter
                                      );


/***************************************************************************
 *  Author:     : JPS
 *  Description : Set gatedependent loudspeaker treshold
 ***************************************************************************/
void lslimiter_setLslimAdjust(LSLIMITER_PTR lslimPtr, float val);

/***************************************************************************
 *  Author:     : JPS
 *  Description : Set gatedependent aerl multiplication factor C
 ***************************************************************************/
void lslimiter_setLslimC(LSLIMITER_PTR lslimPtr, float val);

/***************************************************************************
 *  Author:     : JPS
 *  Description : print status of debug variables
 ***************************************************************************/
void lslimiter_status(LSLIMITER_PTR lslimPtr);

/***************************************************************************
 *  Author:     : JPS
 *  Description : Turn dbgdraw on|off
 ***************************************************************************/
void lslimiter_setdebug(LSLIMITER_PTR lslimPtr, bool val);

/***************************************************************************
 *  Author:     : GVA
 ***************************************************************************/
void unittest_lslimiter(void);

/***************************************************************************
 *  Author:     : GVA
 ***************************************************************************/
LSLIMITER_PTR loadLsLimiter(LSLIMITER_PTR pLsLimiterIn);

/***************************************************************************
 *  Author:     : GVA
 ***************************************************************************/
void flushLsLimiter(LSLIMITER_PTR pLsLimiterSdram, LSLIMITER_PTR pLsLimiterIram);


#ifdef __cplusplus
}
#endif

#endif /* LSLIMITER_H */
