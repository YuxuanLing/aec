#ifndef _AEC_DELAYLINE_H_
#define _AEC_DELAYLINE_H_

#include <list>

class AECDelayLine
{
private:
  int frameSize;
  int frameSizeMs;
  std::list<float *> line;

  void Expand();
  void Reduce();
  void expandRegion(float *frameA, float **newFrame, float *frameC);
  void reduceRegion(float **frameA, float *frameToBeRemoved, float **frameC);

public:
  AECDelayLine(int delay, int samplerate_khz, int framesize);
  ~AECDelayLine(void);

  unsigned int GetDelay();
  void SetDelay(int msDelay);

  bool Process(float *inBuffer, float **outBuffer);

  static void SelfTest(const char *inputFile, const char *outputFile, int delay);
};

#endif /* _AEC_DELAYLINE_H_ */
