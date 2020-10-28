#ifndef APA_H
#define APA_H
/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Author        : BWI,IFA, GEO
 *
 * Description   : Header file for APA. APA is a submodule of AEC-submodule
 *                 echocomp. (AEC=Acoustic Echo Canceller)
 *
 *
 ***************************************************************************/

#include "echocomp.h"

/* AEC-naming conventions now enforced (ifa) - names should be reviewed?  */

void apa_calcFiltCorr(ECHOCOMP_PTR pEchocomp,
                      const int subband_start,
                      const int subband_end);


void apa_calcR11_and_R12(ECHOCOMP_PTR pEchocomp,
                         const int subband_start,
                         const int subband_end);

float apa_filter_APA2(ECHOCOMP_PTR pEchocomp,
                     COMPLEX32 * est,
                     COMPLEX32 * est2,
                     const int subband_start,
                     const int subband_end);


void apa_adapt_APA(ECHOCOMP_PTR pEchocomp,
                   const COMPLEX32 *  adbuf,
                   const COMPLEX32 *  fxbuf,
                   const int subband_start,
                   const int subband_end);

short apa_findbest_APA(ECHOCOMP_PTR  pEchocomp,
                       int k,
                       const int subband_start,
                       const int subband_end);
#endif
