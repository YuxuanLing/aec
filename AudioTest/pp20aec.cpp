/***************************************************************************
*                              A U D I O
*--------------------------------------------------------------------------
*                  (C) Copyright Tandberg Telecom AS 2007
*==========================================================================
*
* Author        : Håvard Graff
*
*
*
* Description   :
*
* Note    : <notes here>
*
* Documentation : <Xxxxxx-Document.doc>
*
**************************************************************************/

/*-------------------------------------------------------------------------
* -- INCLUDE FILES AREA
*-------------------------------------------------------------------------*/
#include "pp20aec.h"

#include "level.h"
#include "ec.h"
#include "lslimiter.h"

#include <stdlib.h>
#include <string.h>

//#ifdef WIN32
#include <malloc.h>
//#endif

#ifdef HAVE_EMMINTRIN_H
#include <emmintrin.h>
#endif

//#if defined(PME_ISIM)
#include "mm_malloc.h"
//#endif

#define nREAD_NOISE_INTO_FAR_END_SIGNAL
#ifdef READ_NOISE_INTO_FAR_END_SIGNAL
#include <stdio.h>
FILE * fidNoise;
static int counter;
#endif

struct PP20AECData
{
  LEVELCONTROL_PTR mLevelCtrlMic;
  LEVELCONTROL_PTR mLevelCtrlNet;

  LSLIMITER_PTR mLsLimiter;

  FILTERBANK_USE_TYPE mFilterbankUse;
  EC_PTR mEc;
  LSPROCESS_PTR mLsProcess;

  float * mLsLimiterDelayline;
  float * mLsProcessBuffer;

  float mAerlInverse;
  int mChannels;
};

PP20AEC::PP20AEC()
{
  pData = new PP20AECData;
  InitAECModules();
}

PP20AEC::~PP20AEC(void)
{
  DestroyAECModules();
  delete pData;
}

void PP20AEC::InitAECModules()
{
  pData->mFilterbankUse = FILTERBANK_EC;
  pData->mChannels = 1; //TODO?

  // Init DC remove modules
  pData->mLevelCtrlMic = levelctrl_create();
  pData->mLevelCtrlNet = levelctrl_create();
  levelctrl_reset (pData->mLevelCtrlMic, AEC_BUFSIZE_SAMPLES);
  levelctrl_reset (pData->mLevelCtrlNet, AEC_BUFSIZE_SAMPLES);

  pData->mLsLimiter = lslimiter_create ();
  lslimiter_init (pData->mLsLimiter, 1); // Second arg not in use...?

  pData->mLsLimiterDelayline = (float * ) _mm_malloc (
    LSLIM_DELAY * sizeof(float), 64);
  memset (pData->mLsLimiterDelayline, 0, LSLIM_DELAY * sizeof(float));

  pData->mLsProcessBuffer = (float * ) _mm_malloc (
    AEC_BUFSIZE_SAMPLES * sizeof(float), 64);
  memset (pData->mLsProcessBuffer, 0, AEC_BUFSIZE_SAMPLES * sizeof(float));

  pData->mLsProcess = lsprocess_create(pData->mFilterbankUse, pData->mChannels);
  lsprocess_init (pData->mLsProcess, pData->mFilterbankUse);

  pData->mEc = ec_create(pData->mFilterbankUse, pData->mChannels);
  ec_init(pData->mEc, pData->mFilterbankUse, pData->mLsProcess);

  pData->mAerlInverse = 0.0f;

#ifdef WRITE_AEC_AUDIO_DATA
  network_in.open           ("aec_network_in.raw"         , std::ios::binary);
  network_out.open          ("aec_network_out.raw"        , std::ios::binary);
  mic_in.open               ("aec_mic_in.raw"             , std::ios::binary);
  spkr_out.open             ("aec_spkr_out.raw"           , std::ios::binary);
#endif

#ifdef READ_NOISE_INTO_FAR_END_SIGNAL
  fidNoise = fopen("noise.bin", "rb");
  counter = 0;
#endif
}

void PP20AEC::DestroyAECModules()
{
  levelctrl_destroy(pData->mLevelCtrlMic);
  levelctrl_destroy(pData->mLevelCtrlNet);

  lslimiter_destroy (pData->mLsLimiter);

  ec_destroy (pData->mEc, pData->mFilterbankUse);
  lsprocess_destroy (pData->mLsProcess, pData->mFilterbankUse);

  _mm_free (pData->mLsLimiterDelayline);
  _mm_free (pData->mLsProcessBuffer);

#ifdef READ_NOISE_INTO_FAR_END_SIGNAL
  fclose(fidNoise);
  counter = 0;
#endif

#ifdef WRITE_AEC_AUDIO_DATA
  network_in.close();
  network_out.close();
  mic_in.close();
  spkr_out.close();
#endif
}

void PP20AEC::SetDelay(float ms)
{
  ec_setWindowsDelay(pData->mEc, ms);
}

void PP20AEC::Reset()
{
    DestroyAECModules();
    delete pData;
    pData = new PP20AECData;
    InitAECModules();
}

bool PP20AEC::Process(float *inpFromNw,  //from network
                      float *inpFromMic, //from mic
                      float *outToLs,    //to loudspeaker
                      float *outToNw,    //to network
                      int samples)       //how many samples from network and mic
{

  if (samples % AEC_BUFSIZE_SAMPLES)
  {
    //    LOG_WARNING_F("Error! Number of samples into process() of aec48 is invalid.");
    return false;
  }

  int frames = samples / AEC_BUFSIZE_SAMPLES;

  float *curNetworkInput;
  float *curMicInput;
  float *curNetworkOutput;
  float *curSpeakerOutput;

#ifdef READ_NOISE_INTO_FAR_END_SIGNAL
  if (counter > 300 && counter < (300+1500))
    fread(inpFromNw, sizeof(float), 480*2, fidNoise);
  counter++;
#endif

  for(int k = 0; k < frames; k++)
  {
    /* advance buffer pointers one frame */
    curNetworkInput   = inpFromNw  + (k * AEC_BUFSIZE_SAMPLES);
    curMicInput       = inpFromMic + (k * AEC_BUFSIZE_SAMPLES);
    curNetworkOutput  = outToNw    + (k * AEC_BUFSIZE_SAMPLES);
    curSpeakerOutput  = outToLs    + (k * AEC_BUFSIZE_SAMPLES);

    /* adjust DC and measure level on signals */
    float abs_level_net;
    float abs_level_mic;
    uint32_t DcRemoveDoneNet [] = {0, 0};
    uint32_t DcRemoveDoneMic [] = {0, 0};
    levelctrl_findLevelAndDcRemove(pData->mLevelCtrlNet, 1,
      &abs_level_net, DcRemoveDoneNet, curNetworkInput, true, 0);
    levelctrl_findLevelAndDcRemove(pData->mLevelCtrlMic, 1,
      &abs_level_mic, DcRemoveDoneMic, curMicInput, true, 0);

    float mic_gain = 1.0; //FIXME ?
    float speaker_gain = 1.0; //FIXME - this should reflect the pre-adjustments done to the volume before it arrives here...
    // If we tell the PP20AEC about our leveling in GUI, it will adapt faster!

    lslimiter_process(pData->mLsLimiter, &curNetworkInput, &curSpeakerOutput,
      &pData->mLsLimiterDelayline, pData->mChannels, pData->mAerlInverse, NULL);

    lsprocess_setloudsGain(pData->mLsProcess, speaker_gain); //Set speaker gain for optimal adaptionspeed

    // Need a separate buffer to give to lsprocess_load
    memcpy (pData->mLsProcessBuffer, curSpeakerOutput, AEC_BUFSIZE_SAMPLES * sizeof (float));
    lsprocess_load(pData->mLsProcess, pData->mLsProcessBuffer);

    ec_load(pData->mEc);
    ec_process(pData->mEc, pData->mLsProcess, curMicInput, curNetworkOutput, true);

    pData->mAerlInverse = ec_getAerl(pData->mEc, pData->mLsProcess, mic_gain);
  }


#ifdef WRITE_AEC_AUDIO_DATA
  network_in.write  (reinterpret_cast<char*>(inpFromNw ), samples * sizeof(float));
  network_out.write (reinterpret_cast<char*>(outToNw   ), samples * sizeof(float));
  mic_in.write      (reinterpret_cast<char*>(inpFromMic), samples * sizeof(float));
  spkr_out.write    (reinterpret_cast<char*>(outToLs   ), samples * sizeof(float));
#endif

  return true;
}
