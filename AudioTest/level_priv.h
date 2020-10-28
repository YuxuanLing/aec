#ifndef _LEVEL_PRIV_H_
#define _LEVEL_PRIV_H_
/*****************************************************************************
 *                              A U D I O
 *
 *----------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2003
 *============================================================================
 *
 * Revision      : $Revision: 1.6 $
 *
 * Description   : Findlevel, ramp and DC-remove.
 *                 A LEVELCONTROL struct is assosiated with a spesific data stream
 * Documentation :
 ****************************************************************************/

#include <stdbool.h>

#define LEVELCTRL_NO_ADJUST_THR (1.0e-6f)
#define DCREMOVERATE         0.004f
#define LEVEL_SIGNATTACK     1.0f        /* <1.0f>  */
#define LEVEL_SIGNDECAY      0.133975f   /* <0.25f> */
#define LEVEL_NOISATTACK     2.44170e-4f /* <0.000488281f> */
#define LEVEL_NOISDECAY      0.292893f   /* <0.5f> */
#define LEVEL_NOISEMIN       5.0e-5f
#define LEVEL_SIGMIN         5.0e-5f


typedef struct LEVELCONTROL
{
  float sigAbsLev;
  float noiseAbsLev;

  float PrevGain;
  float DcOffset[2];
  float PrevDcOffset[2];
  bool  DcRemoveOn;

  int BufLength;
  float InvBufLength;
} LEVELCONTROL;

#endif /* _LEVEL_PRIV_H_ */
