/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 *
 * Switches      : <COMPILER SWITCH HERE>
 *           <add description here>
 *
 * Description   : Findlevel, ramp and DC-remove.
 *
 * Note      :
 *
 * Documentation :
 *
 * ------------------------------------------------------------------------
 * Major changes made (complement to cvs log)
 * ------------------------------------------------------------------------
 * yyyy-mm-dd    : <signature>
 *         Modification made.
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "level.h"
#include "level_priv.h"
#include "mathfun.h"

/** 
 * Finds average abslevel for a frame. Not callable from outside to ensure DC is
 * removed first. If so findLevelAndDcRemove will call this function.
 */
static void levelctrl_findLevel(LEVELCONTROL_PTR  levelCtrl,
                                float* AbsLevel,
                                int nChan,
                                float *pData,
                                bool consecutiveData,
                                float * pData_ch2);

/** 
 * Calculates abs levels of signal and noise. Must only be called once for each
 * frame otherwise time constants will be wrong. Therefore called right after
 *  DC has been removed.
 */
static void levelctrl_findSigAndNoise(LEVELCONTROL_PTR levelCtrl, float currinlevel);

/** 
 * Create levelctrl instance
 *
 * @retval : pointer to struct
 */
LEVELCONTROL_PTR levelctrl_create()
{
    LEVELCONTROL_PTR pThis;
    pThis = (LEVELCONTROL_PTR) malloc(sizeof(LEVELCONTROL));
    assert(NULL != pThis);
    return pThis;
}

/** 
 * Destroy levelctrl instance
 */
void levelctrl_destroy(LEVELCONTROL_PTR levelcontrol_ptr)
{
    free(levelcontrol_ptr);
}

/** 
 * Initialises levelcontrol
 */
void levelctrl_reset(LEVELCONTROL_PTR levelCtrl, int buffersize)
{
  int i;
  levelCtrl->PrevGain   = 0.0f;
  levelCtrl->DcRemoveOn = true;
  for (i=0; i<2; i++)
  {
    levelCtrl->DcOffset[i]     = 0.0f;
    levelCtrl->PrevDcOffset[i] = 0.0f;
  }
  levelCtrl->BufLength    = buffersize;
  levelCtrl->InvBufLength = 1.0f/((float)buffersize);
  levelCtrl->sigAbsLev    = 0.0f;
  levelCtrl->noiseAbsLev  = 0.0f;
}

/** 
 * Remove dc and adjust gain for AUDINT_STREAM_DATA_STRUCT 
 * - mono or stereo
 *   then updates power, gain, etc. if the gain is within a certain distance
 *   from 1.0 the function will return immediately
 */
void levelctrl_rampAndDcRemove(LEVELCONTROL_PTR  levelCtrl,
                               int nChan,
                               float * AbsLevel,
                               uint32_t * DcRemoveDone,
                               float * Gain,
                               float *pData)

{
    int k,ch;
    float stepsizex2,DcRemoveStepsizex2;
    float ramp1,ramp2,tmp1,tmp2,out1,out2,fm1,fm2,fp1,fp2;
    float * bufstart,*  bufptr;
    float DcOffset1,DcOffset2;
    float stepsize,OldGn,meanLevel,absLevel;
    float initialDcOffset,DcRemoveStepsize,InvBufLength;
    int bufLength;

    int AllDcRemovesFinished =
        ((!levelCtrl->DcRemoveOn) ||
         (DcRemoveDone[0] &&  DcRemoveDone[1]) ||
         (DcRemoveDone[0] && (nChan != 2)) );

    if (AllDcRemovesFinished)
    {
        levelctrl_ramp(levelCtrl,
                       nChan,
                       AbsLevel,
                       DcRemoveDone,
                       Gain,
                       pData);
        return;
    }

    if((nChan>2)||(nChan<1))
    {
      printf("rampAndDcRemove: err, channels = %d\n",
      nChan);
      nChan = 1;
    }

    /* init variables */
    bufLength    = levelCtrl->BufLength;
    InvBufLength = levelCtrl->InvBufLength;
    stepsize = (*Gain-levelCtrl->PrevGain) * InvBufLength;
    stepsizex2         = 2.0f * stepsize;

    /* ramp from OldGn */
    OldGn=levelCtrl->PrevGain;
    /* store for next frame */
    levelCtrl->PrevGain = *Gain;

    for(ch=0;ch<nChan;ch++)
    {
        /* init DC-remove variables */
        if(DcRemoveDone[ch]==1)
        {
            initialDcOffset  = 0.0f;
            DcRemoveStepsize = 0.0f;
        }
        else
        {
            initialDcOffset  = levelCtrl->PrevDcOffset[ch];
            DcRemoveStepsize = (levelCtrl->DcOffset[ch] - initialDcOffset) * InvBufLength;
        }
        DcRemoveStepsizex2 = 2.0f * DcRemoveStepsize;
        levelCtrl->PrevDcOffset[ch] = levelCtrl->DcOffset[ch];

        /* prepare loop-variables */
        bufstart = pData+ch*bufLength; /* points to start of channel-data. (L(/R)) */

        fm1 = 0.0f;
        fm2 = 0.0f;

        fp1 = 0.0f;
        fp2 = 0.0f;

        bufptr = bufstart;

        DcOffset1 = initialDcOffset + DcRemoveStepsize;
        DcOffset2 = initialDcOffset + DcRemoveStepsizex2;

        ramp1 = OldGn + stepsize;
        ramp2 = OldGn + stepsizex2;

        for (k = 0; k < bufLength; k+=2)
        {

            tmp1       = *bufptr;
            fm1       += tmp1;
            out1       = ramp1 * (tmp1 - DcOffset1);
            fp1       += fabsf(out1);
            *bufptr++  = out1;
            ramp1     += stepsizex2;
            DcOffset1 += DcRemoveStepsizex2;

            tmp2       = *bufptr;
            fm2       += tmp2;
            out2       = ramp2 * (tmp2 - DcOffset2);
            fp2       += fabsf(out2);
            *bufptr++  = out2;
            ramp2     += stepsizex2;
            DcOffset2 += DcRemoveStepsizex2;

        }
        /* After adjust */
        meanLevel = (fm1+fm2) * InvBufLength;
        DcRemoveDone[ch] = 1;
        AbsLevel[ch] = (fp1+fp2) * InvBufLength;
        /* Compute amount the signal is to be adjusted next time DcReomve runs because of DC components*/
        levelCtrl->DcOffset[ch] += (meanLevel - levelCtrl->DcOffset[ch]) * DCREMOVERATE;
    }
    /* After adjust */
    *Gain = 1.0f;

    absLevel = AbsLevel[0];  /*  mono  */
    if(nChan == 2)       /* stereo */
    {
        absLevel = AbsLevel[0] + AbsLevel[1];
    }
    levelctrl_findSigAndNoise(levelCtrl, absLevel);
}

/** 
 * Removes dc and computes findlevel. Keeps track of wether DCremove has already
 * been done, and in that case will only update abslevel.
 */
void levelctrl_findLevelAndDcRemove(LEVELCONTROL_PTR  levelCtrl,
                                    int nChan,
                                    float * AbsLevel,
                                    uint32_t * DcRemoveDone,
                                    float *pData,
                                    bool consecutiveData,
                                    float *pData_ch2)

{
    float tmp1,tmp2,out1,out2;
    float fp1,fp2,fm1,fm2;
    float * bufstart,*  bufptr;
    float DcOffset1,DcOffset2;
    float initialDcOffset,DcRemoveStepsize,DcRemoveStepsizex2,meanLevel,InvBufLength;
    float absLevel;
    int bufLength,k,ch;

    int AllDcRemovesFinished =
        ((!levelCtrl->DcRemoveOn) ||
         (DcRemoveDone[0] && DcRemoveDone[1]) ||
         (DcRemoveDone[0] && nChan != 2));

    if (AllDcRemovesFinished)
    {
        levelctrl_findLevel(levelCtrl, AbsLevel, nChan, pData, consecutiveData, pData_ch2);
        return;
    }

    if((nChan>2)||(nChan<1))
    {
      printf("findLevelAndDcRemove: err, channels = %d\n",
      nChan);
      nChan = 1;
    }

    /* init variables */
    bufLength = levelCtrl->BufLength;
    InvBufLength = levelCtrl->InvBufLength;

    for(ch=0;ch<nChan;ch++)
    {
        /* init DC-remove variables */
        if(DcRemoveDone[ch]==1)
        {
            initialDcOffset  = 0.0f;
            DcRemoveStepsize = 0.0f;
        }
        else
        {
            initialDcOffset  = levelCtrl->PrevDcOffset[ch];
            DcRemoveStepsize = (levelCtrl->DcOffset[ch] - initialDcOffset) * InvBufLength;
        }

        DcRemoveStepsizex2 = 2.0f * DcRemoveStepsize;
        levelCtrl->PrevDcOffset[ch] = levelCtrl->DcOffset[ch];

        /* prepare loop-variables */
        if(!consecutiveData) 
        {
            if(ch == 0) {bufstart = pData;    }
            else        {bufstart = pData_ch2;}

        }
        else
        {
            bufstart = pData+ch*bufLength; /* points to start of channel-data. (L(/R)) */
        }

        fm1 = 0.0f;
        fm2 = 0.0f;

        fp1 = 0.0f;
        fp2 = 0.0f;

        bufptr = bufstart;

        DcOffset1 = initialDcOffset + DcRemoveStepsize;
        DcOffset2 = initialDcOffset + DcRemoveStepsizex2;

        for (k = 0; k < bufLength; k+=2)
        {
            tmp1       = *bufptr;
            fm1       += tmp1;
            out1       = tmp1 - DcOffset1;
            *bufptr++  = out1;
            fp1       += fabsf(out1);
            DcOffset1 += DcRemoveStepsizex2;

            tmp2       = *bufptr;
            fm2       += tmp2;
            out2       = tmp2 - DcOffset2;
            *bufptr++  = out2;
            fp2       += fabsf(out2);
            DcOffset2 += DcRemoveStepsizex2;
        }
        /* After adjust */
        meanLevel = (fm1+fm2) * InvBufLength;
        DcRemoveDone[ch] = 1;
        AbsLevel[ch] = (fp1+fp2) * InvBufLength;
        /* Compute amount the signal is to be adjusted nest time DcReomve runs because of DC components*/
        levelCtrl->DcOffset[ch] += (meanLevel - levelCtrl->DcOffset[ch]) * DCREMOVERATE;

    }
    /* Compute average power */
    absLevel = AbsLevel[0];  /*  mono  */
    if (nChan == 2)               /* stereo */
    {
        absLevel = AbsLevel[0] + AbsLevel[1];
    }

    levelctrl_findSigAndNoise(levelCtrl, absLevel);

}

/** 
 * Calculates AbsLevel for mono or stereo AUDINT_STREAM_DATA_STRUCT. Can, will
 * and shall ONLY be called from levelctrl_findLevelAndDcRemove as that function
 * will keep track of wheter it is nessecary to remove DC first, if this is not
 * done AbsLevel will be rubbish.
 *
 * Modified by TVJ:
 *  Added levelctrl_findSigAndNoise to make sure it is done in
 *  levelctrl_findLevelAndDcRemove. Assuming levelctrl_findLevelAndDcRemove
 *  is called with a different LEVELCONTROL_PTR if called more than once per frame.
 */
static void levelctrl_findLevel(LEVELCONTROL_PTR  levelCtrl,
                                float* AbsLevel,
                                int nChan,
                                float *pData,
                                bool consecutiveData,
                                float * pData_ch2)
{
    int k, ch;
    float fp1,fp2,fp3,fp4;
    float *  bufptr;
    float * bufstart;
    float InvBufLength;
    float absLevel;
    int bufLength;

    if((nChan>2)||(nChan<1))
    {
        printf("findlevel: err, channels = %d\n",
                nChan);
        nChan = 1;
    }

    /* Init variables */
    bufLength=levelCtrl->BufLength;
    InvBufLength=levelCtrl->InvBufLength;

    for(ch=0;ch<nChan;ch++)
    {
        if(!consecutiveData) 
        {
            if(ch == 0) {bufstart = pData;    }
            else        {bufstart = pData_ch2;}

        }
        else
        {
            bufstart = pData+ch*bufLength; /* points to start of channel-data. (L(/R)) */
        }

        bufstart = pData+ch*bufLength; /* points to start of channel-data. (L(/R)) */

        bufptr = bufstart;
        
        fp1 = 0.0f;
        fp2 = 0.0f;
        fp3 = 0.0f;
        fp4 = 0.0f;

        for (k = 0; k < bufLength; k+=4)
        {
            fp1 += fabsf(*bufptr++);
            fp2 += fabsf(*bufptr++);
            fp3 += fabsf(*bufptr++);
            fp4 += fabsf(*bufptr++);
        }
        AbsLevel[ch] = ((fp1+fp2)+(fp3+fp4))*InvBufLength;
    }
    
    /* Compute average power.
       Assuming levelctrl_findLevelAndDcRemove is called with a different
       LEVELCONTROL_PTR if called more than once per frame.*/
    absLevel = AbsLevel[0];  /*  mono  */
    if (nChan == 2)               /* stereo */
    {
        absLevel = AbsLevel[0] + AbsLevel[1];
    }
    levelctrl_findSigAndNoise(levelCtrl, absLevel);
}

/** 
 * Adjust gain for AUDINT_STREAM_DATA_STRUCT
 * - mono or stereo then updates gain, etc.
 *   if the gain is within a certain distance from,
 *   1.0 the function will return immediately
 */
void levelctrl_ramp(LEVELCONTROL_PTR  levelCtrl,
                    int nChan,
                    float * AbsLevel,
                    uint32_t * DcRemoveDone,
                    float * Gain,
                    float *pData)
{
    int bufLength,k,ch;
    float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
    double tmp12,tmp34,tmp56,tmp78;
    float ramp1,ramp2,ramp3,ramp4,ramp5,ramp6,ramp7,ramp8;
    float * bufstart,*  bufptr;
    float stepsize,stepsizex8,OldGn,InvBufLength;
    int AllDcRemovesFinished = 
        ((!levelCtrl->DcRemoveOn) ||
         (DcRemoveDone[0] &&  DcRemoveDone[1])||
         (DcRemoveDone[0] && (nChan != 2)));

    if (!AllDcRemovesFinished)
    {
        printf("ramp: err, DC not removed!\n");
        levelctrl_rampAndDcRemove(levelCtrl,
                                  nChan,
                                  AbsLevel,
                                  DcRemoveDone,
                                  Gain,
                                  pData);
        return;
    }

    if(fabsf(*Gain-1.0)<LEVELCTRL_NO_ADJUST_THR &&
       fabsf(levelCtrl->PrevGain-1.0)<LEVELCTRL_NO_ADJUST_THR)

    {
        return;
    }

    if((nChan>2)||(nChan<1))
    {
      printf("ramp: err, channels = %d\n",
      nChan);
      nChan = 1;
    }

    /* init variables */
    bufLength    = levelCtrl->BufLength;
    InvBufLength = levelCtrl->InvBufLength;
    stepsize     = (*Gain-levelCtrl->PrevGain) * InvBufLength;
    stepsizex8   = 8.0f * stepsize;
    bufstart     = pData;

    /* ramp from OldGn */
    OldGn=levelCtrl->PrevGain;
    /* store for next frame */
    levelCtrl->PrevGain=*Gain;

    for(ch=0;ch<nChan;ch++)
    {
        bufstart = pData+ch*bufLength; /* points to start of channel-data. (L(/R)) */

        bufptr = bufstart;

        ramp1 = OldGn +        stepsize;
        ramp2 = OldGn + 2.0f * stepsize;
        ramp3 = OldGn + 3.0f * stepsize;
        ramp4 = OldGn + 4.0f * stepsize;
        ramp5 = OldGn + 5.0f * stepsize;
        ramp6 = OldGn + 6.0f * stepsize;
        ramp7 = OldGn + 7.0f * stepsize;
        ramp8 = OldGn +        stepsizex8;
        for (k = 0; k < bufLength; k+=8)
        {
            tmp12 = (double)_amemd8((double*)&(bufptr[0]));
            tmp34 = (double)_amemd8((double*)&(bufptr[2]));
            tmp56 = (double)_amemd8((double*)&(bufptr[4]));
            tmp78 = (double)_amemd8((double*)&(bufptr[6]));

            tmp1  = _itof(_lo(tmp12));
            tmp2  = _itof(_hi(tmp12));
            tmp3  = _itof(_lo(tmp34));
            tmp4  = _itof(_hi(tmp34));
            tmp5  = _itof(_lo(tmp56));
            tmp6  = _itof(_hi(tmp56));
            tmp7  = _itof(_lo(tmp78));
            tmp8  = _itof(_hi(tmp78));

            *bufptr++   = ramp1 * tmp1;
            ramp1     += stepsizex8;

            *bufptr++   = ramp2 * tmp2;
            ramp2     += stepsizex8;

            *bufptr++   = ramp3 * tmp3;
            ramp3     += stepsizex8;

            *bufptr++   = ramp4 * tmp4;
            ramp4     += stepsizex8;

            *bufptr++   = ramp5 * tmp5;
            ramp5     += stepsizex8;

            *bufptr++   = ramp6 * tmp6;
            ramp6     += stepsizex8;

            *bufptr++   = ramp7 * tmp7;
            ramp7     += stepsizex8;

            *bufptr++   = ramp8 * tmp8;
            ramp8     += stepsizex8;
        }
        /* After adjust */
        AbsLevel[ch] *= *Gain;
    }
    /* After adjust */
    *Gain = 1.0f;
}

float findAbsLevel(float *pData, int buflen)
{
    int k;
    float fp1,fp2,fp3,fp4;

    fp1 = 0.0f;
    fp2 = 0.0f;
    fp3 = 0.0f;
    fp4 = 0.0f;

    for (k = 0; k < buflen; k+=4)
    {
        fp1 += fabsf(*pData++);
        fp2 += fabsf(*pData++);
        fp3 += fabsf(*pData++);
        fp4 += fabsf(*pData++);
    }
    return ((fp1+fp2+fp3+fp4)/buflen);
}

static void levelctrl_findSigAndNoise(LEVELCONTROL_PTR levelCtrl,
                                      float currinlevel)
{
  if (levelCtrl->sigAbsLev<currinlevel)
    levelCtrl->sigAbsLev +=
      LEVEL_SIGNATTACK * (currinlevel - levelCtrl->sigAbsLev);
  else
    levelCtrl->sigAbsLev +=
      LEVEL_SIGNDECAY  * (currinlevel - levelCtrl->sigAbsLev);

  levelCtrl->sigAbsLev = max(LEVEL_NOISEMIN, levelCtrl->sigAbsLev);

  if (levelCtrl->noiseAbsLev<currinlevel)
    levelCtrl->noiseAbsLev +=
      LEVEL_NOISATTACK * (currinlevel - levelCtrl->noiseAbsLev);
  else
    levelCtrl->noiseAbsLev +=
      LEVEL_NOISDECAY  * (currinlevel - levelCtrl->noiseAbsLev);

  levelCtrl->noiseAbsLev = max(LEVEL_SIGMIN, levelCtrl->noiseAbsLev);

}

float levelctrl_getSNR(LEVELCONTROL_PTR  levelCtrl)
{
    if((int)(levelCtrl->noiseAbsLev * 1000.f) == 0)
       return 0.0f;
    else
       return fasterdiv(levelCtrl->sigAbsLev,levelCtrl->noiseAbsLev);
}

#if 0
/** 
 * This function adds the signals in pAudBuf1 and pAudBuf2 and writes the result
 * back to pAudBuf1. The signals are also faded in/out independently from
 *  previous gain held in levelCtrl1/2 to new gain set in pAudBuf1/2
 */
void levelctrl_rampAndAdd(LEVELCONTROL_PTR  levelCtrl1,
                          LEVELCONTROL_PTR  levelCtrl2,
                          AUDINT_STREAM_DATA_STRUCT *  pAudBuf1,
                          AUDINT_STREAM_DATA_STRUCT *  pAudBuf2)
{
    int k;
    float tmp11,tmp12,tmp13,tmp14,tmp21,tmp22,tmp23,tmp24;
    float ramp11,ramp12,ramp13,ramp14,ramp21,ramp22,ramp23,ramp24;
    float *  bufptr11,*  bufptr12,*  bufptr13,
        *  bufptr14;
    float *  bufptr21,*  bufptr23;
    float *  bufptr1L,*  bufptr1R,*  bufptr2L,
        *  bufptr2R;
    float stepsize1,OldGn1,stepsize2,OldGn2,stepsize1x4,stepsize2x4;
    int bufLength;
    int nChan1,nChan2;
    int AllDcRemovesFinished1,AllDcRemovesFinished2;

    AllDcRemovesFinished1 =
        ((!levelCtrl1->DcRemoveOn) ||
         (pAudBuf1->header.DcRemoveDone[0] && pAudBuf1->header.DcRemoveDone[1]) ||
         (pAudBuf1->header.DcRemoveDone[0] && (pAudBuf1->header.iChannels != 2)));

    if (!AllDcRemovesFinished1)
    {
        printf("rampAndAdd: err1, DC not removed!, %d, %x\n",
               pAudBuf1->header.iSequenceCnt, (unsigned int)pAudBuf1);
        levelctrl_findLevelAndDcRemove(
            levelCtrl1, &pAudBuf1->header, pAudBuf1->pData, true, 0);
    }

    AllDcRemovesFinished2 =
        ((!levelCtrl2->DcRemoveOn) ||
         (pAudBuf2->header.DcRemoveDone[0] &&  pAudBuf2->header.DcRemoveDone[1]) ||
         (pAudBuf2->header.DcRemoveDone[0] && (pAudBuf2->header.iChannels != 2)));

    if (!AllDcRemovesFinished2)
    {
        printf("rampAndAdd: err2, DC not removed!, %d, %x\n",
               pAudBuf2->header.iSequenceCnt,(unsigned int)pAudBuf2);
        levelctrl_findLevelAndDcRemove(
            levelCtrl2, &pAudBuf2->header, pAudBuf2->pData, true, 0);
    }

    nChan1    = pAudBuf1->header.iChannels;
    nChan2    = pAudBuf2->header.iChannels;

    if((nChan1>2)||(nChan1<1))
    {
        printf("rampAndAdd: err1, channels = %d\n",nChan1);
        nChan1 = 1;
    }

    if((nChan2>2)||(nChan2<1))
    {
        printf("rampAndAdd: err2, channels = %d\n",nChan2);
        nChan2 = 1;
    }

    /* init variables*/
    OldGn1    = levelCtrl1->PrevGain;
    OldGn2    = levelCtrl2->PrevGain;
    stepsize1 = (pAudBuf1->header.Gain - OldGn1) * levelCtrl1->InvBufLength;
    stepsize2 = (pAudBuf2->header.Gain - OldGn2) * levelCtrl2->InvBufLength;
    bufLength = levelCtrl1->BufLength;

    bufptr1L = pAudBuf1->pData;
    bufptr1R = &pAudBuf1->pData[bufLength];

    bufptr2L = pAudBuf2->pData;
    bufptr2R = &pAudBuf2->pData[bufLength];

    /* store for next frame*/
    levelCtrl1->PrevGain = pAudBuf1->header.Gain;
    levelCtrl2->PrevGain = pAudBuf2->header.Gain;

    stepsize1x4 = 4.0f * stepsize1;
    stepsize2x4 = 4.0f * stepsize2;

    ramp11 = OldGn1 +        stepsize1;
    ramp12 = OldGn1 + 2.0f * stepsize1;
    ramp13 = OldGn1 + 3.0f * stepsize1;
    ramp14 = OldGn1 +        stepsize1x4;

    ramp21 = OldGn2 +        stepsize2;
    ramp22 = OldGn2 + 2.0f * stepsize2;
    ramp23 = OldGn2 + 3.0f * stepsize2;
    ramp24 = OldGn2 +        stepsize2x4;

    bufptr11 = bufptr1L;
    bufptr12 = bufptr1L + 1;
    bufptr13 = bufptr1L + 2;
    bufptr14 = bufptr1L + 3;

    bufptr21 = bufptr2L;
    bufptr23 = bufptr2L + 2;

    if((nChan1==1)&&(nChan2==1))
    {
double bufSrc112,bufSrc134,bufSrc212,bufSrc234;

#pragma MUST_ITERATE(20, 240, 20);
        for (k = 0; k < bufLength; k+=4)
        {
            bufSrc112 = _amemd8((double*)&(bufptr11[k]));
            bufSrc134 = _amemd8((double*)&(bufptr13[k]));
            bufSrc212 = _amemd8((double*)&(bufptr21[k]));
            bufSrc234 = _amemd8((double*)&(bufptr23[k]));

            tmp11  = _itof(_lo(bufSrc112));
            tmp12  = _itof(_hi(bufSrc112));
            tmp13  = _itof(_lo(bufSrc134));
            tmp14  = _itof(_hi(bufSrc134));
            tmp21  = _itof(_lo(bufSrc212));
            tmp22  = _itof(_hi(bufSrc212));
            tmp23  = _itof(_lo(bufSrc234));
            tmp24  = _itof(_hi(bufSrc234));

            bufptr11[k] = ramp11*tmp11 + ramp21*tmp21;
            ramp11     += stepsize1x4;
            ramp21     += stepsize2x4;

            bufptr12[k] = ramp12*tmp12 + ramp22*tmp22;
            ramp12     += stepsize1x4;
            ramp22     += stepsize2x4;

            bufptr13[k] = ramp13*tmp13 + ramp23*tmp23;
            ramp13     += stepsize1x4;
            ramp23     += stepsize2x4;

            bufptr14[k] = ramp14*tmp14 + ramp24*tmp24;
            ramp14     += stepsize1x4;
            ramp24     += stepsize2x4;

        }
    }
    else if((nChan1==2)&&(nChan2==1))
    {
        float buf2L1,buf2L2,buf1L1,buf1L2,buf1R1,buf1R2;
        double buf2L12,buf1L12,buf1R12;
        float stepsize1x2 = 2.0f * stepsize1;
        float stepsize2x2 = 2.0f * stepsize2;

        ramp11 = OldGn1 + stepsize1;
        ramp12 = OldGn1 + stepsize1x2;

        ramp21 = OldGn2 + stepsize2;
        ramp22 = OldGn2 + stepsize2x2;

#pragma MUST_ITERATE(40, 240, 40);
        for (k = 0; k < bufLength; k+=2)
        {
            buf2L12 = _amemd8((double*)&(bufptr2L[k]));
            buf1L12 = _amemd8((double*)&(bufptr1L[0]));
            buf1R12 = _amemd8((double*)&(bufptr1R[0]));

            buf2L1  = _itof(_lo(buf2L12));
            buf2L2  = _itof(_hi(buf2L12));
            buf1L1  = _itof(_lo(buf1L12));
            buf1L2  = _itof(_hi(buf1L12));
            buf1R1  = _itof(_lo(buf1R12));
            buf1R2  = _itof(_hi(buf1R12));

            tmp21       = ramp21*buf2L1;
            *bufptr1L++ = ramp11*buf1L1 + tmp21;
            *bufptr1R++ = ramp11*buf1R1 + tmp21;
            tmp22       = ramp22*buf2L2;
            *bufptr1L++ = ramp12*buf1L2 + tmp22;
            *bufptr1R++ = ramp12*buf1R2 + tmp22;

            ramp11 += stepsize1x2;
            ramp21 += stepsize2x2;
            ramp12 += stepsize1x2;
            ramp22 += stepsize2x2;
        }
    }
    else if((nChan1==1)&&(nChan2==2))
    {
        float stepsize1x2 = 2.0f * stepsize1;
        float stepsize2x2 = 2.0f * stepsize2;

        ramp11 = OldGn1 + stepsize1;
        ramp12 = OldGn1 + stepsize1x2;

        ramp21 = OldGn2 + stepsize2;
        ramp22 = OldGn2 + stepsize2x2;

#pragma MUST_ITERATE(40, 240, 40);
        for (k = 0; k < bufLength; k+=2)
        {
            bufptr1L[k]   = ramp11*bufptr1L[k]   + 0.5f*ramp21*(bufptr2L[k]  +bufptr2R[k]);
            bufptr1L[k+1] = ramp12*bufptr1L[k+1] + 0.5f*ramp22*(bufptr2L[k+1]+bufptr2R[k+1]);

            ramp11 += stepsize1x2;
            ramp21 += stepsize2x2;
            ramp12 += stepsize1x2;
            ramp22 += stepsize2x2;
        }
    }
    else if((nChan1==2)&&(nChan2==2))
    {
        float stepsize1x2 = 2.0f * stepsize1;
        float stepsize2x2 = 2.0f * stepsize2;
        float bufSrc1,bufSrc2,bufSrc3,bufSrc4,bufSrc5,bufSrc6,bufSrc7,bufSrc8;
        double bufSrc12,bufSrc34,bufSrc56,bufSrc78;
        ramp11 = OldGn1 + stepsize1;
        ramp12 = OldGn1 + stepsize1x2;

        ramp21 = OldGn2 + stepsize2;
        ramp22 = OldGn2 + stepsize2x2;

#pragma MUST_ITERATE(40, 240, 40);
#pragma UNROLL(1)
        for (k = 0; k < bufLength; k+=2)
        {
            bufSrc12 = _amemd8((double*)&(bufptr1L[0]));
            bufSrc34 = _amemd8((double*)&(bufptr2L[k]));
            bufSrc56 = _amemd8((double*)&(bufptr1R[0]));
            bufSrc78 = _amemd8((double*)&(bufptr2R[k]));

            bufSrc1  = _itof(_lo(bufSrc12));
            bufSrc2  = _itof(_hi(bufSrc12));
            bufSrc3  = _itof(_lo(bufSrc34));
            bufSrc4  = _itof(_hi(bufSrc34));
            bufSrc5  = _itof(_lo(bufSrc56));
            bufSrc6  = _itof(_hi(bufSrc56));
            bufSrc7  = _itof(_lo(bufSrc78));
            bufSrc8  = _itof(_hi(bufSrc78));

            *bufptr1L++ = ramp11*bufSrc1 + ramp21*bufSrc3;
            *bufptr1R++ = ramp11*bufSrc5 + ramp21*bufSrc7;

            *bufptr1L++ = ramp12*bufSrc2 + ramp22*bufSrc4;
            *bufptr1R++ = ramp12*bufSrc6 + ramp22*bufSrc8;

            ramp11 += stepsize1x2;
            ramp21 += stepsize2x2;

            ramp12 += stepsize1x2;
            ramp22 += stepsize2x2;
        }
    }

    pAudBuf1->header.Gain = 1.0;
    pAudBuf2->header.Gain = 1.0;
}
#endif 

float levelctrl_getLevel(LEVELCONTROL_PTR levelCtrl)
{
    return levelCtrl->sigAbsLev;
}

void ramp(float *pSrc,  /* Must be 64-bits aligned */
          float oldGain,
          float newGain,
          int   buflen, /* Must be a multiple of 8 */
          int   nChan)
{
    int k,ch;
    float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
    double tmp12,tmp34,tmp56,tmp78;
    float ramp1,ramp2,ramp3,ramp4,ramp5,ramp6,ramp7,ramp8;
    float *  bufptr;
    float stepsize,stepsizex8;

    if((nChan>2)||(nChan<1))
    {
        return;
    }

    /* init variables */
    stepsize     = (newGain-oldGain)/buflen;
    stepsizex8   = 8.0f * stepsize;
    bufptr       = pSrc;

    for(ch=0;ch<nChan;ch++)
    {
        ramp1 = oldGain +        stepsize;
        ramp2 = oldGain + 2.0f * stepsize;
        ramp3 = oldGain + 3.0f * stepsize;
        ramp4 = oldGain + 4.0f * stepsize;
        ramp5 = oldGain + 5.0f * stepsize;
        ramp6 = oldGain + 6.0f * stepsize;
        ramp7 = oldGain + 7.0f * stepsize;
        ramp8 = oldGain +        stepsizex8;
        for (k = 0; k < buflen; k+=8)
        {
            tmp12 = (double)_amemd8((double*)&(bufptr[0]));
            tmp34 = (double)_amemd8((double*)&(bufptr[2]));
            tmp56 = (double)_amemd8((double*)&(bufptr[4]));
            tmp78 = (double)_amemd8((double*)&(bufptr[6]));

            tmp1  = _itof(_lo(tmp12));
            tmp2  = _itof(_hi(tmp12));
            tmp3  = _itof(_lo(tmp34));
            tmp4  = _itof(_hi(tmp34));
            tmp5  = _itof(_lo(tmp56));
            tmp6  = _itof(_hi(tmp56));
            tmp7  = _itof(_lo(tmp78));
            tmp8  = _itof(_hi(tmp78));

            *bufptr++   = ramp1 * tmp1;
            ramp1     += stepsizex8;

            *bufptr++   = ramp2 * tmp2;
            ramp2     += stepsizex8;

            *bufptr++   = ramp3 * tmp3;
            ramp3     += stepsizex8;

            *bufptr++   = ramp4 * tmp4;
            ramp4     += stepsizex8;

            *bufptr++   = ramp5 * tmp5;
            ramp5     += stepsizex8;

            *bufptr++   = ramp6 * tmp6;
            ramp6     += stepsizex8;

            *bufptr++   = ramp7 * tmp7;
            ramp7     += stepsizex8;

            *bufptr++   = ramp8 * tmp8;
            ramp8     += stepsizex8;
        }
    }
}


float levelctrl_getSignalLevel(LEVELCONTROL_PTR  levelCtrl)
{
    return levelCtrl->sigAbsLev;
}


float levelctrl_getNoiseLevel(LEVELCONTROL_PTR  levelCtrl)
{
    return levelCtrl->noiseAbsLev;
}

