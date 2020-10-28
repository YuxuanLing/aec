/***************************************************************************
 *                        D E E S S E R   M O D U L E
 *==========================================================================
 *
 * Description   : Module for detecting and reducing sibilants.
 *
 *                 Based on the article "Adaptive Algorithm for Detecting
 *                 and Reducing Sibilands in Recorded Speech" by Martin
 *                 Wolters, Markus Sapp and Joerg Becker-Schweitzer.
 *
 **************************************************************************/
#include "deesser.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

static float criticalBandIntensity[CRITICAL_BANDS];

static int   criticalBandBoundries[CRITICAL_BANDS] = {
    50, 100, 150, 200, 250, 300, 350, 400, 450, 510, 570, 630, 700, 770, 840,
    920, 1000, 1080, 1170, 1270, 1370, 1480, 1600, 1720, 1850,
    2000, 2150, 2320, 2500, 2700, 2900, 3150, 3400, 3700, 4000, 4400, 4800,
    5300, 5800, 6400, 7000, 7700, 8500, 9500, 10500, 12000, 13500, 15500
};

static float alpha[CRITICAL_BANDS] = {
    0.0128,    0.0150,    0.0176,    0.0207,    0.0245,    0.0291,    0.0347,    0.0416,
    0.0501,    0.0607,    0.0739,    0.0904,    0.1111,    0.1374,    0.1707,    0.2132,
    0.2675,    0.3373,    0.4271,    0.5428,    0.6920,    0.8840,    1.1303,    1.4442,
    1.8399,    2.3307,    2.9252,    3.6203,    4.3914,    5.1785,    5.8735,    6.3161,
    6.3147,    5.7101,    4.4797,    2.8317,    1.1745,   -0.1130,   -0.9668,   -1.6657,
   -2.5853,   -4.0800,   -6.6584,  -11.3317,  -20.3345,  -39.0856,  -82.3873, -197.7311
};

static float gain[CRITICAL_BANDS] = {
    1.0893,    1.1865,    1.2924,    1.4078,    1.5334,    1.6703,    1.8194,    1.9818,
    2.1587,    2.3514,    2.5613,    2.7899,    3.0389,    3.3102,    3.6056,    3.9275,
    4.2781,    4.6599,    5.0759,    5.5290,    6.0225,    6.5601,    7.1456,    7.7835,
    8.4782,    9.2350,   10.0593,   10.9572,   11.9353,   13.0007,   14.1611,   15.4252,
   16.8020,   18.3018,   19.9355,   21.7149,   23.6532,   25.7646,   28.0643,   30.5694,
   33.2981,   36.2703,   39.5079,   43.0344,   46.8757,   51.0599,   55.6176,   60.5821
};

static DEESSER deEsser = {
    INIT_BW_THRESHOLD,
    INIT_S_THRESHOLD,
    INIT_GAIN,
    INIT_FADE_STEPS,
    INIT_LOGGING,
    INIT_DEESSER_ON,
    INIT_SHOW_SHARPNESS
};

DEESSER_PTR deesser_getDeEsser()
{
    return &deEsser;
}


void deesser_calcCriticalBandIntensities(COMPLEX32 * micfft)
{
    int i, fftIndex;
    float relCriticalBandwidth;
    float freq, startRest, endRest;

    fftIndex = 0;
    startRest = endRest = 0.0f;
    freq = FFT_WIDTH/2;
    for (i=0;i<CRITICAL_BANDS;i++)
    {
        criticalBandIntensity[i] = 0.0f;
    }

    for (i=0; i<CRITICAL_BANDS; i++)
    {
        while ( (freq+=FFT_WIDTH) < criticalBandBoundries[i])
        {
            fftIndex++;
            criticalBandIntensity[i] += ( micfft[fftIndex].re*micfft[fftIndex].re +
                                          micfft[fftIndex].im*micfft[fftIndex].im );
        }
        fftIndex++;
        endRest = (FFT_WIDTH - (freq - criticalBandBoundries[i]))/FFT_WIDTH;
        criticalBandIntensity[i] += (micfft[fftIndex].re*micfft[fftIndex].re +
                                     micfft[fftIndex].im*micfft[fftIndex].im) * endRest;

        // Start next critical bandwidth?
        if (i<CRITICAL_BANDS-1)
        {
            startRest = 1-endRest;
            relCriticalBandwidth = (criticalBandBoundries[i+1] - criticalBandBoundries[i])/FFT_WIDTH;
            // Check if one whole critical bandwidth is within an fft bin.
            // This will only happen when i and the critical bandwidth is small. Last test is for coverity.
            if ( (startRest > relCriticalBandwidth) && (i<CRITICAL_BANDS-2) )
            {
                criticalBandIntensity[i+1] = (micfft[fftIndex].re*micfft[fftIndex].re +
                                              micfft[fftIndex].im*micfft[fftIndex].im) * relCriticalBandwidth;
                startRest -= relCriticalBandwidth;
                i++;
            }
            criticalBandIntensity[i+1] += (micfft[fftIndex].re*micfft[fftIndex].re +
                                           micfft[fftIndex].im*micfft[fftIndex].im) * startRest;
        }
    }
}


float deesser_calcSharpness()
{
    float loudness, sharpness;
    float numerator, denominator;
    int i;

    numerator = denominator = 0;

    for (i=0; i<CRITICAL_BANDS; i++)
    {
        loudness = sqrtf(sqrtf(alpha[i]*criticalBandIntensity[i]));
        numerator   += loudness*gain[i];
        denominator += loudness;
    }
    sharpness = LOUDNESS_CONST*numerator/denominator;

    if (deEsser.showSharpness)
    {
        printf("\nSharpness = %f\n",sharpness);
    }

    return sharpness;
}


void deesser_reduceSibilants(COMPLEX32 * micfft)
{
    float max = 0.0f;
    int maxIndex = 0;
    int lowCutOffIndex, highCutOffIndex, i;

    for (i=0; i<CRITICAL_BANDS; i++)
    {
        if (criticalBandIntensity[i] > max)
        {
            max = criticalBandIntensity[i];
            maxIndex = i;
        }
    }

    // Do nothing if peak is below 15 bark
    if (maxIndex < CRITICAL_BAND_PEAK_MINIMUM)
    {
        return;
    }

    // find lower limit
    i = maxIndex-1;
    while ( (criticalBandIntensity[i]>deEsser.bwThreshold*criticalBandIntensity[maxIndex]) && (i>1) )
    {
        i--;
    }
    lowCutOffIndex = (int) ( (criticalBandBoundries[i]-FFT_WIDTH/2)/FFT_WIDTH );

    // find upper limit
    i = maxIndex+1;
    while ( (criticalBandIntensity[i]>deEsser.bwThreshold*criticalBandIntensity[maxIndex]) && (i<CRITICAL_BANDS-1) )
    {
        i++;
    }
    highCutOffIndex = (int) ( (criticalBandBoundries[i]-FFT_WIDTH/2)/FFT_WIDTH );

    if (deEsser.logging)
    {
        printf("Max freq: %d\n",criticalBandBoundries[maxIndex]);
        printf("Filter bank cutoff freqs: %f, %f\n", 
               31.25f+62.5*lowCutOffIndex, 31.25f+62.5*highCutOffIndex);
    }

#ifdef DEESSER_DEBUG
    if (deEsser.logging)
    {
        printf("\npre = [");
        for (i=0;i<320;i++)
        {
            printf(" %f", micfft[i].re*micfft[i].re + micfft[i].im*micfft[i].im);
        }
        printf("];\n");
    }
#endif

    // gain fade:
    for (i=1; i<deEsser.fadeSteps; i++)
    {
        // real and imag start fade
        micfft[lowCutOffIndex-1+i].re *= (1-((1-deEsser.gain) * i/(float)deEsser.fadeSteps));
        micfft[lowCutOffIndex-1+i].im *= (1-((1-deEsser.gain) * i/(float)deEsser.fadeSteps));

        // real and imag end fade
        micfft[highCutOffIndex+1-i].re *= (1-((1-deEsser.gain) * i/(float)deEsser.fadeSteps));
        micfft[highCutOffIndex+1-i].im *= (1-((1-deEsser.gain) * i/(float)deEsser.fadeSteps));
    }

    for (i=lowCutOffIndex+(deEsser.fadeSteps-1); i<=highCutOffIndex-(deEsser.fadeSteps-1); i++)
    {
        micfft[i].re *= deEsser.gain;
        micfft[i].im *= deEsser.gain;
    }

#ifdef DEESSER_DEBUG
    if (deEsser.logging)
    {
        printf("\npost = [");
        for (i=0;i<320;i++)
        {
            printf(" %f", micfft[i].re*micfft[i].re + micfft[i].im*micfft[i].im);
        }
        printf("];\n");
    }
#endif
}


void deesser_status()
{
    printf("\n");
    printf("DeEssing is turned %s\n", deEsser.deEsserOn ? "on" : "off");
    printf("DeEsser parameters:\n");
    printf("Bandwidth threshold  = %f\n",deEsser.bwThreshold);
    printf("Sharpness threshold  = %f\n",deEsser.sharpnessThreshold);
    printf("Gain                 = %f\n",deEsser.gain);
    printf("Number of fade steps = %d\n",deEsser.fadeSteps);
    printf("Logging is turned %s\n", deEsser.logging ? "on" : "off");
}


void deesser_process(COMPLEX32 * micfft)
{
    deesser_calcCriticalBandIntensities(micfft);
    if (deesser_calcSharpness() > deEsser.sharpnessThreshold)
    {
        deesser_reduceSibilants(micfft);
    }
}

