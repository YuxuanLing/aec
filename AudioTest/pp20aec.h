#ifndef __PP20_AEC_H__
#define __PP20_AEC_H__

#define nWRITE_AEC_AUDIO_DATA

#ifdef WRITE_AEC_AUDIO_DATA
#include <fstream>
#endif

const int AEC_SAMPLERATE_KHZ = 48;
const int AEC_BUFSIZE_MS = 10;
const int AEC_BUFSIZE_SAMPLES = AEC_BUFSIZE_MS * AEC_SAMPLERATE_KHZ;

struct PP20AECData;

class PP20AEC
{
public:
  PP20AEC();
  ~PP20AEC();

  bool Process(float *inpFromNw,  //from network
               float *inpFromMic, //from mic
               float *outToLs,    //to loudspeaker
               float *outToNw,    //to network
               int samples);      //how many samples from network and mic

  void SetDelay(float ms);
  void Reset();

private:
  PP20AECData *pData;

  void InitAECModules();
  void DestroyAECModules();

#ifdef WRITE_AEC_AUDIO_DATA
  std::ofstream network_in;
  std::ofstream network_out;
  std::ofstream mic_in;
  std::ofstream spkr_out;
#endif

};

#endif /* __PP20_AEC_H__ */
