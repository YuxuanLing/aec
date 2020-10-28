#include "delayEstimation.h"
#include "delayEstimation_priv.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef WRITE_DEBUG_TO_FILE
#include <locale.h>
#endif

#ifdef __ARM_NEON__
#include <arm_neon.h>
//#pragma thumb           // for compatibility with gcc
#endif

struct DELAY_ESTIMATION * delayEstimation_create()
{
    struct DELAY_ESTIMATION * delayEstimation = (struct DELAY_ESTIMATION *) malloc(sizeof(DELAY_ESTIMATION));
    if( delayEstimation == NULL )
    {
        fprintf(stderr, "delayEstimation_create: Could not allocate delayEstimation struct\n");
        return 0;
    }
#ifdef WRITE_DEBUG_TO_FILE
    fidEc = fopen("delay.txt", "w");
    setlocale( LC_NUMERIC, "English" );
#endif
    return delayEstimation;
}

void delayEstimation_destroy(struct DELAY_ESTIMATION * delayEstimation)
{
    free(delayEstimation);
#ifdef WRITE_DEBUG_TO_FILE
    fclose(fidEc);
#endif

}

void delayEstimation_init(struct DELAY_ESTIMATION * delayEstimation)
{
    int i,j,k;
    delayEstimation->lsPower           = 0.0f;
    delayEstimation->currentFrameDelay = 0;;
    delayEstimation->newFrameDelay     = 0;
    delayEstimation->decisionCounter   = 0;
    delayEstimation->lsLineWriteIndex  = 0;
    delayEstimation->lsLineReadIndex   = 0;
    delayEstimation->debug             = 0;
    delayEstimation->windowsdelay      = 0.0f;
    delayEstimation->maxDecisionCount  = DELAY_ESTIMATION_MAX_DECISION_COUNT;
    delayEstimation->offset            = DELAY_ESTIMATION_AEC_OFFSET;
    delayEstimation->doDelaying        = 1;
    delayEstimation->useDelay          = 0;

    for( i = 0; i < DELAY_ESTIMATION_MAX_LAG; i++ )
    {
        delayEstimation->sumCrossCorrelation[i] = 0.0f;
        for( j = 0; j < DELAY_ESTIMATION_SUBUSED; j++ )
        {
            delayEstimation->crossCorrelation[i][j].re = 0.0f;
            delayEstimation->crossCorrelation[i][j].im = 0.0f;
        }
        for( k = 0; k < FFTSIZE/2; k++ )
        {
            delayEstimation->lsLine[i][k].re = 0.0f;
            delayEstimation->lsLine[i][k].im = 0.0f;
        }
    }
}

void delayEstimation_measurePower(float * power, const COMPLEX32 * lsbuf)
{
    int m;
    int n = SUBUSED_START + DELAY_ESTIMATION_SUBUSED_START; 

    float sum = 0.0f;
    for( m = 0; m < DELAY_ESTIMATION_SUBUSED; m++, n++ )
    {
        sum += lsbuf[n].re * lsbuf[n].re + lsbuf[n].im * lsbuf[n].im;
    }
    *power = *power * (1-DELAY_ESTIMATION_UPDATE_SPEED) + sum * DELAY_ESTIMATION_UPDATE_SPEED;
}

void delayEstimation_addDataToLsLine(struct DELAY_ESTIMATION * delayEstimation, const COMPLEX32 *  lsbuf)
{
    int j = delayEstimation->lsLineWriteIndex;
    memcpy(delayEstimation->lsLine[j], lsbuf, sizeof(COMPLEX32)*FFTSIZE/2);

    delayEstimation->lsLineWriteIndex = (delayEstimation->lsLineWriteIndex + 1) % DELAY_ESTIMATION_MAX_LAG;
    /* lsLineWriteIndex now points to where the next row of delayline should be written. (rotating modulo DELAY_ESTIMATION_MAX_LAG)) */
}

void delayEstimation_crossCorrelate(struct DELAY_ESTIMATION * delayEstimation, const COMPLEX32 *  micFft)
{
    int i;
    int k = delayEstimation->lsLineReadIndex; /* lsLine is not rotated, lsLineWriteIndex-1 denotes newest row. */

#if defined(__ARM_NEON__) && defined(ENV_IOS)
    micFft += SUBUSED_START + DELAY_ESTIMATION_SUBUSED_START;
#endif 

    for ( i=0; i < DELAY_ESTIMATION_MAX_LAG; i++ )
    {
        int m;

#if defined(__ARM_NEON__) && defined(ENV_IOS)	// should work on Android, but leaving for now
        COMPLEX32 *lsLine = delayEstimation->lsLine[k] + SUBUSED_START + DELAY_ESTIMATION_SUBUSED_START;
        COMPLEX32 *crossC = delayEstimation->crossCorrelation[i];
    
        for( m = 0; m < DELAY_ESTIMATION_SUBUSED; m+=4 )
        {
            float32x4_t a, b, c, d, tmp_re, tmp_im, cc_re, cc_im;
          
            // interleaved real/imag, vld2q will de-interleave
            float32x4x2_t tmpd = vld2q_f32((float32_t *)&micFft[m]);
          
            a = tmpd.val[0];  // real
            b = tmpd.val[1];  // imag
          
            tmpd = vld2q_f32((float32_t *)&lsLine[m]);
          
            c = tmpd.val[0];  // real
            d = tmpd.val[1];  // imag

            tmp_re = vaddq_f32(vmulq_f32(a, c), vmulq_f32(b, d));
            tmp_im = vsubq_f32(vmulq_f32(b, c), vmulq_f32(a, d));

            tmpd = vld2q_f32((float32_t *)&crossC[m]);
          
            cc_re = tmpd.val[0];  // real
            cc_im = tmpd.val[1];  // imag

            tmpd.val[0] = vaddq_f32(cc_re, vmulq_n_f32(vsubq_f32(tmp_re, cc_re), DELAY_ESTIMATION_UPDATE_SPEED));
            tmpd.val[1] = vaddq_f32(cc_im, vmulq_n_f32(vsubq_f32(tmp_im, cc_im), DELAY_ESTIMATION_UPDATE_SPEED));

            vst2q_f32((float32_t *)&crossC[m], tmpd);
        }
#else
        int j = SUBUSED_START + DELAY_ESTIMATION_SUBUSED_START;
        for( m = 0; m < DELAY_ESTIMATION_SUBUSED; m++, j++ )
        {
            float a = micFft[j].re; 
            float b = micFft[j].im;
            float c = delayEstimation->lsLine[k][j].re;
            float d = delayEstimation->lsLine[k][j].im;

            float tmpRe = a * c + b * d;
            float tmpIm = b * c - a * d;

            // The usual "x = (1-a)*x + a*y" form requires two multiplications;
            // rewriting it as "x += a*(y-x)" requires only one.
            delayEstimation->crossCorrelation[i][m].re += DELAY_ESTIMATION_UPDATE_SPEED * (tmpRe - delayEstimation->crossCorrelation[i][m].re);
            delayEstimation->crossCorrelation[i][m].im += DELAY_ESTIMATION_UPDATE_SPEED * (tmpIm - delayEstimation->crossCorrelation[i][m].im);
        }
#endif
        /* next row (k - 1), modulo DELAY_ESTIMATION_MAX_LAG */
        if (--k < 0)
            k += DELAY_ESTIMATION_MAX_LAG;

        /* note that lsLine has opposite orientation than in the MatLab model.
           model: [new, 2., 3., ...] here: (lsLine:) [4., 3., 2, new, 10., 9., ...]
           (Still, crossCorrelation has the same orientation as ML-model.)*/
    }

    /* lsLineReadIndex rotates modulo DELAY_ESTIMATION_MAX_LAG, it points to where 'newest' data should be read from ecrossCorrelationdline. */
    i = delayEstimation->lsLineReadIndex+1;
    if(i >= DELAY_ESTIMATION_MAX_LAG)
        i -= DELAY_ESTIMATION_MAX_LAG;
    delayEstimation->lsLineReadIndex = i;
}


struct INDEXED_MAX_SUM delayEstimation_findMaxSum(struct DELAY_ESTIMATION * delayEstimation)
{
    struct INDEXED_MAX_SUM current = {0,0.0f};

#if defined(__ARM_NEON__) && defined(ENV_IOS)
    COMPLEX32 *pCC;
    float32x4_t sum_re;
    float32x4_t sum_im;
    float32x2_t sum_2f;
    float sum;
    int i, j;

    for( i = 0; i < DELAY_ESTIMATION_MAX_LAG; i++ )
    {
        pCC = delayEstimation->crossCorrelation[i];
        sum_re = vdupq_n_f32(0.0f);
        sum_im = vdupq_n_f32(0.0f);

        for( j = 0; j < DELAY_ESTIMATION_SUBUSED; j+=4 )
        {
            float32x4x2_t cc = vld2q_f32((float32_t*)&pCC[j]);
            sum_re = vaddq_f32(sum_re, vabsq_f32(cc.val[0]));
            sum_im = vaddq_f32(sum_im, vabsq_f32(cc.val[1]));
        }
         
        sum_re = vaddq_f32(sum_re, sum_im);
        sum_2f = vpadd_f32(vget_high_f32(sum_re), vget_low_f32(sum_re));
        sum = vget_lane_f32(sum_2f, 0) + vget_lane_f32(sum_2f, 1);
        delayEstimation->sumCrossCorrelation[i] = sum;
        if(i==0 || sum > current.max_sum)
        {
            current.index  = i;
            current.max_sum = sum;
        }
    }
#else
    int i;
    for( i = 0; i < DELAY_ESTIMATION_MAX_LAG; i++ )
    {
        float sum = 0.0f;
        int j;
        for( j = 0; j < DELAY_ESTIMATION_SUBUSED; j++ )
        {
            sum += fabsf(delayEstimation->crossCorrelation[i][j].re);
            sum += fabsf(delayEstimation->crossCorrelation[i][j].im);
        }

        if( i==0 )
        {
            current.index  = i;
            current.max_sum = sum;
        }
        else if( sum > current.max_sum )
        {
            current.index  = i;
            current.max_sum = sum;
        }
        delayEstimation->sumCrossCorrelation[i] = sum;
    }
#endif

    return current;
}


bool delayEstimation_evaluateMaxCrossCorrelation(const struct DELAY_ESTIMATION * delayEstimation, const struct INDEXED_MAX_SUM * maxCorrelationValue)
{
    int k; 

    //evaluating indexes before correlation peak
    int last_neighbour = maxCorrelationValue->index;
    for( k = maxCorrelationValue->index-1; k >= 0; k-- )
    {
        if( delayEstimation->sumCrossCorrelation[k] > (maxCorrelationValue->max_sum * 0.5) )  //only significant correlation values are evaluated
        {
            if( (last_neighbour-k) != 1 )
            {
                //fprintf(stdout, "not increasing neighbours, at k = %d \n", k);
                return false;
            }
            else
                last_neighbour = k;
            if( delayEstimation->sumCrossCorrelation[k+1] < delayEstimation->sumCrossCorrelation[k] )
            {
                //fprintf(stdout, "not increasing correlation value, at k = %d, k+1: %9.7f, k: %9.7f \n", k, delayEstimation->sumCrossCorrelation[k+1], delayEstimation->sumCrossCorrelation[k]);
                return false;
            }
            if( delayEstimation->sumCrossCorrelation[k] > maxCorrelationValue->max_sum *0.75 )  //high correlation on lower neighbours, dangerous to use new estimate
            {
                //fprintf(stdout, "high correlation on lower neighbours, dangerous to use new estimate, at k = %d \n", k);
                return false;
            }
        }
    }

    //evaluating indexes after correlation peak
    last_neighbour = maxCorrelationValue->index;
    for( k = maxCorrelationValue->index+1; k < DELAY_ESTIMATION_MAX_LAG; k++ )
    {
        if( delayEstimation->sumCrossCorrelation[k] > (maxCorrelationValue->max_sum * 0.5) )  //only significant correlation values are evaluated
        {
            if( (k-last_neighbour) != 1 )
            {
                //fprintf(stdout, "not decreasing neighbours, at k = %d \n", k);
                return false;
            }
            else
                last_neighbour = k;
            if( delayEstimation->sumCrossCorrelation[k] > delayEstimation->sumCrossCorrelation[k-1] )
            {
                //fprintf(stdout, "not decreasing correlation value, at k = %d, k+1: %9.7f, k: %9.7f \n", k, delayEstimation->sumCrossCorrelation[k], delayEstimation->sumCrossCorrelation[k-1]);
                return false;
            }
        }
    }
    return true;
}

int delayEstimation_decideUponDelayEstimate(struct DELAY_ESTIMATION * delayEstimation, int estimatedDelay)
{
	int MovedDelay = 0;
    if( estimatedDelay == delayEstimation->currentFrameDelay )
        return MovedDelay;
    else if( estimatedDelay < 0 )
        return MovedDelay;
    else if( estimatedDelay == delayEstimation->newFrameDelay )
    {
        delayEstimation->decisionCounter++;
        if( delayEstimation->decisionCounter > delayEstimation->maxDecisionCount - 1)
            delayEstimation->decisionCounter = delayEstimation->maxDecisionCount - 1;
        //fprintf(stdout,"delaycount++ = %d \n", delayEstimation->decisionCounter);
    }
    else
    {
        delayEstimation->newFrameDelay = estimatedDelay;    
        delayEstimation->decisionCounter = 0;
        //fprintf(stdout,"delayEstimation: newFrameDelay changed to: %d \n", delayEstimation->newFrameDelay);
    }

	if( delayEstimation->decisionCounter == delayEstimation->maxDecisionCount - 1 )
    {
		delayEstimation->currentFrameDelay = delayEstimation->newFrameDelay;
		
		if( delayEstimation->currentFrameDelay < delayEstimation->offset)
		{
			delayEstimation->useDelay = 0;
		}
		else
		{//only change delay we use if maxTap moves outside tap 3-8 in filter - to avoid unecessary re-adaption 
			int difference = delayEstimation->currentFrameDelay - (delayEstimation->useDelay + delayEstimation->offset);
			if( (difference < 3-delayEstimation->offset) || (difference > 8-delayEstimation->offset) ) //offset == 5 
			{
				MovedDelay = delayEstimation->useDelay - (delayEstimation->currentFrameDelay - delayEstimation->offset);
				delayEstimation->useDelay = delayEstimation->currentFrameDelay - delayEstimation->offset;
			}
		}
    }
    return MovedDelay;
}


void delayEstimation_delayLsBuf(const struct DELAY_ESTIMATION * delayEstimation, COMPLEX32 *  lsFft)
{
    int index = (delayEstimation->lsLineReadIndex - 1) - (delayEstimation->useDelay);

    while (index < 0 )
      index += DELAY_ESTIMATION_MAX_LAG;  //make sure index is positive before doing modulo
    index %= DELAY_ESTIMATION_MAX_LAG;

    memcpy(lsFft, delayEstimation->lsLine[index], sizeof(COMPLEX32)*FFTSIZE/2 );
}

int delayEstimation_estimateDelay(struct DELAY_ESTIMATION * delayEstimation, const COMPLEX32 *  micFft, COMPLEX32 *  lsFft)
{
    int movedDelay = 0;
	struct INDEXED_MAX_SUM maxCorrelationValue;
    delayEstimation_measurePower(&(delayEstimation->lsPower), lsFft);
    delayEstimation_addDataToLsLine(delayEstimation, lsFft);
    delayEstimation_crossCorrelate(delayEstimation, micFft);
    maxCorrelationValue = delayEstimation_findMaxSum(delayEstimation);

	if (maxCorrelationValue.index == 0 || maxCorrelationValue.index == DELAY_ESTIMATION_MAX_LAG)
	{
		;
		// TODO: Fix for movi
		//fprintf(stderr, "delayEstimation: estimated framedelay reached boundary: %d \n", maxCorrelationValue.index);
	}

    if( delayEstimation->lsPower > DELAY_ESTIMATION_MIN_LSPOWER && (maxCorrelationValue.index != 0 && maxCorrelationValue.index != DELAY_ESTIMATION_MAX_LAG) )
    {
        bool delayIsValid = delayEstimation_evaluateMaxCrossCorrelation(delayEstimation, &maxCorrelationValue);
        if( delayIsValid )
        {
            movedDelay = delayEstimation_decideUponDelayEstimate(delayEstimation, maxCorrelationValue.index);
        }
    }

	if(delayEstimation->currentFrameDelay == 0)
		delayEstimation->useDelay = (int) delayEstimation->windowsdelay;

	if( delayEstimation->doDelaying ) 
		delayEstimation_delayLsBuf(delayEstimation, lsFft);

    if( delayEstimation->debug == 70 )
    {
        static int counter = 0;
        static float Sum0 = 0.0f, Sum1 = 0.0f;
        Sum0 += delayEstimation->currentFrameDelay;
        Sum1 += delayEstimation->newFrameDelay;
        if( ++counter == 20 )
        {
            printf("  %8.2e, %8.2e\n", Sum0*0.05f, Sum1*0.05f);
            counter = 0;
            Sum0 = 0.0f;
            Sum1 = 0.0f;
        }
    }

#ifdef WRITE_DEBUG_TO_FILE
	fprintf(fidEc, " %4.2f %d %7.5e %d ",delayEstimation->windowsdelay, delayEstimation->currentFrameDelay, delayEstimation->lsPower, delayEstimation->useDelay);
#endif
    return movedDelay;
}


void delayEstimation_printStatus(struct DELAY_ESTIMATION * delayEstimation)
{
    printf("\rStatus: DelayEstimation");
    printf("\r\n   delay            - %d",     delayEstimation->currentFrameDelay);
    printf("\r\n   channels         - %d",     delayEstimation->newFrameDelay);
    printf("\r\n   debug            - %d",     delayEstimation->decisionCounter);
    printf("\r\n");
}

void delayEstimation_setDebug(struct DELAY_ESTIMATION * delayEstimation, int value)
{
    delayEstimation->debug = value;
}

void delayEstimation_setWindowsDelay(struct DELAY_ESTIMATION * delayEstimation, float value)
{
	delayEstimation->windowsdelay = value * 0.1f;
}

void delayEstimation_setDelayingOn(struct DELAY_ESTIMATION * delayEstimation, int onoff)
{
	delayEstimation->doDelaying = onoff;
}

void delayEstimation_returnIntermediateLsBuf(const struct DELAY_ESTIMATION * delayEstimation, COMPLEX32 *  lsFft, int movedDelay, int delay_position)
{
    int index = (delayEstimation->lsLineReadIndex - 1) - (delayEstimation->useDelay + movedDelay - delay_position);

    while (index < 0 )
      index += DELAY_ESTIMATION_MAX_LAG;  //make sure index is positive before doing modulo
    index %= DELAY_ESTIMATION_MAX_LAG;

    memcpy(lsFft, delayEstimation->lsLine[index], sizeof(COMPLEX32)*FFTSIZE/2 );
}

void delayEstimation_returnHistoryLsSample(const struct DELAY_ESTIMATION * delayEstimation, 
                                           COMPLEX32 *  lsComplexSample, 
                                           int movedDelay, 
                                           int filtlen, 
                                           int subband,
                                           int delay_position)
{
    int pre_shift_offset = 1;
    int index = (delayEstimation->lsLineReadIndex - 1) - (delayEstimation->useDelay + filtlen + movedDelay + delay_position + pre_shift_offset);

    while (index < 0 )
      index += DELAY_ESTIMATION_MAX_LAG;  //make sure index is positive before doing modulo
    index %= DELAY_ESTIMATION_MAX_LAG;

    lsComplexSample->re = delayEstimation->lsLine[index][SUBUSED_START+subband-1].re;
    lsComplexSample->im = delayEstimation->lsLine[index][SUBUSED_START+subband-1].im;
}

