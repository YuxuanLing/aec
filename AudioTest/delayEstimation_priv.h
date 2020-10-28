#ifndef DELAY_ESTIMATION_PRIV_H
#define DELAY_ESTIMATION_PRIV_H

#define DELAY_ESTIMATION_MAX_LAG (50)
#define DELAY_ESTIMATION_SUBUSED_START (5 - SUBUSED_START) 
#define DELAY_ESTIMATION_SUBUSED (20)
#define DELAY_ESTIMATION_MAX_DECISION_COUNT (2)
#define DELAY_ESTIMATION_UPDATE_SPEED (0.0039f)   //power update speed (2^(-8))
#define DELAY_ESTIMATION_MIN_LSPOWER (3.8147e-6f) //minimum loudspeaker power to evaluate correlation: 1/(4096*64)
#define DELAY_ESTIMATION_AEC_OFFSET (5)		 // number of frames that the aec loudspeaker buffer is ahead of max correlation tap

typedef struct DELAY_ESTIMATION 
{
    COMPLEX32 lsLine[DELAY_ESTIMATION_MAX_LAG][FFTSIZE/2];
    COMPLEX32 crossCorrelation[DELAY_ESTIMATION_MAX_LAG][DELAY_ESTIMATION_SUBUSED];
    float sumCrossCorrelation[DELAY_ESTIMATION_MAX_LAG];
    float lsPower;
    int   currentFrameDelay;
    int   newFrameDelay;
    int   decisionCounter;
    int   lsLineWriteIndex;
    int   lsLineReadIndex;
    int   debug;
	float windowsdelay;
	int   offset;      //aec wants at least 20 ms before max correlation peak for best possible adaption 
	int   maxDecisionCount;
	int   doDelaying;
	int   useDelay;
} DELAY_ESTIMATION;


#endif /* DELAY_ESTIMATION_PRIV_H */
