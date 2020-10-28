#ifndef SYNTH_PRIV_H
#define SYNTH_PRIV_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : Geir Ole �verby, Tandberg Telecom AS
 * Co-author     : Johan Malmstigen, S&T ES
 *
 * Switches      : -
 *
 * Description   : Private data for the synth module
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
 **************************************************************************/

/*-------------------------------------------------------------------------
 * -- INCLUDE FILES AREA
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * -- PREPROCESSOR CONSTANTS AREA
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * -- PRIVATE MACROS AREA
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * -- PRIVATE TYPES AREA
 *-------------------------------------------------------------------------*/
/* -------------------------------------------------*/
/* defines 'struct SYNTH' which contains            */
/* variables etc. that define the synthesis filter. */
/* These should be filled by an init routine.       */
/* Auth: geo                                        */
/* -------------------------------------------------*/

typedef struct SYNTH {
    float *fftst;
    short fftst_ix;
    short fftst_iy;
    float postfiltState;
    int fftStlX;
    int fftStlY;
    int synthDlLength;
    int synthBuildLen;
    const float *sfirSDRAM;
    int filtertype;
    const float *sfir;
    float postfiltercoeff;
    float postfiltergain;
} SYNTH;


#endif /*SYNTH_PRIV_H*/
