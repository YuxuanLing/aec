#include "aecdelayline.h"

//#include "log.h"

#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_EMMINTRIN_H
#include <emmintrin.h>
#endif

#if defined(PME_ISIM)
#include <mm_malloc.h>
#endif

using namespace std;

AECDelayLine::AECDelayLine(int delay, int samplerate_khz, int framesize)
{
  this->frameSize = framesize;

  // assuming samplerate is 48khz
  this->frameSizeMs = framesize / samplerate_khz;

  if (delay < 0) {
    //LOG_WARNING_F("AECDelayLine:: Cannot assign negative delay, defaulting to zero");
    delay = 0;
  }

  int numFrames = delay / frameSizeMs;

  if (frameSizeMs % 10 != 0) {
    //LOG_WARNING_F("AECDelayLine:: These frames each last for %d ms, not divisible by 10", frameSizeMs);
  }

  for (int i = 0; i < numFrames; i++)
  {
    float * mem = (float *)_mm_malloc (frameSize * sizeof(float), 64);
    memset (mem, 0, frameSize * sizeof(float));
    line.push_back(mem);
  }
}

AECDelayLine::~AECDelayLine() {
  std::list<float *>::iterator it;
  for (it = line.begin(); it != line.end(); it++)
    _mm_free (*it);
}

unsigned int AECDelayLine::GetDelay()
{
  return this->frameSizeMs * line.size();
}

void AECDelayLine::SetDelay(int msDelay)
{
  if (msDelay < 0) {
    //LOG_WARNING_F("AECDelayLine:: Cannot assign negative delay, defaulting to zero");
    msDelay = 0;
  }
  //int configuredDelay = line.size()*this->frameSizeMs;
  // We need to make sure that the configured delay is always
  // lower than or equal to the delay observed in the in-out
  // loop. So, always round down to the nearest multiple of
  // the frame size in the delay line.

  unsigned int newCapacity = floor((float)msDelay / (float)this->frameSizeMs);

  unsigned int currentCapacity = line.size();
  if (newCapacity < currentCapacity)
  {
    for (unsigned int i = 0; i < (currentCapacity - newCapacity); i++)
    {
      Reduce();
    }
    //LOG_INFO_F("AECDelayLine:: reducing delay to %d ms", line.size()*this->frameSizeMs);
  }
  else if (newCapacity > currentCapacity)
  {
    for (unsigned int i = 0; i < (newCapacity - currentCapacity); i++)
    {
      Expand();
    }
    //LOG_INFO_F("AECDelayLine:: increasing delay to %d ms", line.size()*this->frameSizeMs);
  }
}

bool AECDelayLine::Process(float *inBuffer, float **outBuffer)
{
  float * newbuffer;
  newbuffer = (float *) _mm_malloc(this->frameSize * sizeof(float), 64);
  memcpy(newbuffer, inBuffer, this->frameSize * sizeof(float));
  line.push_front(newbuffer);

  *outBuffer = line.back();

  line.pop_back();

  return true;
}

void AECDelayLine::Expand()
{
  if (line.size() <= 1) {
    // if we only have one or zero frames to delay, we just insert silence
    float *silence = (float *)_mm_malloc(frameSize * sizeof(float), 64);
    memset (silence, 0, frameSize * sizeof(float));
    line.push_front(silence);
  }
  else {
    list<float *>::iterator iter = line.end();
    list<float *>::iterator last = --iter;
    list<float *>::iterator secondLast = --iter;

    float *newFrame;
    float *frameA = *secondLast;
    float *frameC = *last;
    expandRegion(frameA, &newFrame, frameC);
    line.insert(secondLast, newFrame);
  }
}

void AECDelayLine::Reduce()
{
  if (line.size() <= 0) {
    //LOG_WARNING_F("AECDelayLine:: Cannot reduce the delay of a zero-delay delayLine. No change made.");
    return;
  }

  if (line.size() <= 2) {
    _mm_free(line.front());
    line.pop_front();
  }
  else {
    list<float *>::iterator iter = line.end();
    list<float *>::iterator last = --iter;
    list<float *>::iterator secondLast = --iter;
    list<float *>::iterator thirdLast = --iter;

    float **frameA = &(*thirdLast);
    float **frameC = &(*last);
    reduceRegion(frameA, *secondLast, frameC);
    _mm_free(*secondLast);
    line.erase(secondLast);
  }
}

/**
  Modify the frames around an expand operation. frameA and frameC are input frames,
  and we are left to create newFrame to be placed in-between on the delay line.
*/
void AECDelayLine::expandRegion(float *frameA, float **newFrame, float *frameC)
{
  //LOG_DEBUG_F ("frameA = %f.frameC=%f \n", frameA, frameC);
  // just make some silence for now
  *newFrame = (float *)_mm_malloc (frameSize* sizeof(float), 64);
  memset (*newFrame, 0, frameSize * sizeof(float));
}

/**
  Modify the frames around a reduce operation. We are given all three frames with the
  assumption that the middle will be removed. From here we can modify the surrounding
  frames frameA and frameC
*/
void AECDelayLine::reduceRegion(float **frameA, float *frameToBeRemoved, float **frameC)
{
  //LOG_DEBUG_F ("frameA = %f.frameToBeRemoved = %f, frameC=%f \n",
    //frameA, frameToBeRemoved, frameC);
  // here we could modify surrounding frames, but we'll leave them for now
}
