/***************************************************************************
 *                              A U D I O
 *--------------------------------------------------------------------------
 *                  (C) Copyright Tandberg Telecom AS 2004
 *==========================================================================
 * Author     : Trygve F. Marton, see fastnoisered_mod.m
 * Implemented: Jens Petter Stang (JPS)
 *
 *
 * Desc:    Part of noisereduction. Main objectives:
 *          Remove noise from keypad use by decreasing the gain in subbands
 *              according to an a priori keyEstimate.
 *
 * Notes:   We need information about when a key is pressed. On snoopy this
 *              information is served from the fpga via audlocalin.
 *
 *          (If the mic input is overloaded at the AD, the routine can shift
 *              to mute-mode, instead of using a subband gain, not active)
 *
 *           Jan 12th 2009: only the mute-functionality is active, meaning
 *                  muting all subbands when there is a keyevent. It's a quick
 *                  and dirty implementation of a keyclickremoval, which
 *                  degrades the sound quality if the user is talking and at
 *                  the same time using the keypad.
 *
 **************************************************************************/

#include "keyclickremoval_priv.h"
#include "keyclickremoval.h"
#include "testdata/keyclickremovaltestdata.h"
#include "audec_defs.h"
#include "mathfun.h"
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define nKEYCLICKREMOVAL_DEBUGPRINT

/**
 * @file keyremoval.c
 * @ingroup audlocalin
 * @brief Remove keyclick noise from snoopy keypad in subbands based on an a priori estimate of the key noise.
 * Trigged by fpga when an keyevent occurs. A 0xACDC bitmask is inserted for identifying a keyevent, and the
 * following 48 bits identifies the keyId, one unique bit is set high for each key on the keypad.
 */

static void keyclickremoval_trigger(KEYCLICKREMOVAL_PTR pKeyclickRemoval);


/**
 *
 * @return KEYCLICKREMOVAL_PTR
 */
KEYCLICKREMOVAL_PTR keyclickremoval_create(void)
{
    KEYCLICKREMOVAL_PTR pKeyclickRemoval;
    pKeyclickRemoval = malloc(sizeof(KEYCLICKREMOVAL));
    if( pKeyclickRemoval == NULL )
    {
        fprintf(stderr, "keyclickremoval_create: Could not allocate pKeyclickRemoval \n");
    }
    return pKeyclickRemoval;
}

/**
 *
 * @param pKeyclickRemoval
 */
void keyclickremoval_destroy(KEYCLICKREMOVAL_PTR pKeyclickRemoval)
{
    free(pKeyclickRemoval);
}

/**
 *
 * @param pKeyclickRemoval
 */
void keyclickremoval_init(KEYCLICKREMOVAL_PTR pKeyclickRemoval)
{
    pKeyclickRemoval->keyId = KEY_NONE;
    pKeyclickRemoval->keyState = KEYUP;
    pKeyclickRemoval->limitcount = 0;
    pKeyclickRemoval->limitcountmax = 5;

    pKeyclickRemoval->keyclickremovalOn = true;
    pKeyclickRemoval->removalStarted    = false;

    pKeyclickRemoval->lim_gain = 1.0f;
    pKeyclickRemoval->gi = 1.05f;
    pKeyclickRemoval->gd = 0.95f;
    pKeyclickRemoval->lev = 0.01f;
    pKeyclickRemoval->tc = 0.9995f;
}


static KEYCLICKREMOVAL_PTR keyclickremoval_load(KEYCLICKREMOVAL_PTR pKeyclickRemovalSdram)
{
    /* MEMXFERMUSTDIE */
    return pKeyclickRemovalSdram;
}

/**
 *
 */
void keyclickremoval_limiter_process(KEYCLICKREMOVAL_PTR pKeyclickRemovalSdram, float *  micbuf)
{
    KEYCLICKREMOVAL_PTR pKeyclickRemoval;
    int i;
    float gain_curr, lim_gain, gi, gd, lev;
    
    /* MEMXFERMUSTDIE */

    pKeyclickRemoval = keyclickremoval_load(pKeyclickRemovalSdram);

    if( pKeyclickRemoval->removalStarted && pKeyclickRemoval->keyclickremovalOn )   /* if key is trigged */
    {
        lim_gain = pKeyclickRemoval->lim_gain;
        gi = pKeyclickRemoval->gi;
        gd = pKeyclickRemoval->gd;
        lev = pKeyclickRemoval->lev;

        for (i = 0; i < FRAMESIZE; i++)
        {
            micbuf[i] *= lim_gain;

            gain_curr = gi;

            if(fabsf(micbuf[i]) > lev)
            {
                gain_curr = gd;
            }
            lim_gain *= gain_curr;
            if (lim_gain > 1.0f)
            {
                lim_gain = 1.0f;
            }
        }
        pKeyclickRemoval->lim_gain = lim_gain;

        if( ++pKeyclickRemoval->limitcount >= pKeyclickRemoval->limitcountmax ) /* removal is finished */
        {
            pKeyclickRemoval->removalStarted = false;
            pKeyclickRemoval->keyId = KEY_NONE;
            pKeyclickRemoval->limitcount = 0;
        }
    }
    else
    {
        for (i = 0; i < FRAMESIZE; i++)
        {
            micbuf[i] *= pKeyclickRemoval->lim_gain;
            pKeyclickRemoval->lim_gain *= pKeyclickRemoval->gi;
            if (pKeyclickRemoval->lim_gain > 1.0f)
            {
                pKeyclickRemoval->lim_gain = 1.0f;
            }
            pKeyclickRemoval->lev = pKeyclickRemoval->lev*pKeyclickRemoval->tc + (1.0f-pKeyclickRemoval->tc)*fabsf(micbuf[i]);
        }
    }
}

/**
 * Sets necessary parameters for the keyremovalprocess and trigger the routine to start
 *
 * @param pKeyclickRemoval
 * @param keytrigged
 * @param timing
 * @param keyId
 */
static void keyclickremoval_trigger(KEYCLICKREMOVAL_PTR pKeyclickRemoval)
{
        if(pKeyclickRemoval->keyState == KEYUP && pKeyclickRemoval->keyId != KEY_NONE)  /* new key-down event */
        {
            #ifdef KEYCLICKREMOVAL_DEBUGPRINT
            printf("keyclickremoval_trigger: key-down \n");
            #endif
            pKeyclickRemoval->keyState   = KEYDOWN;
            pKeyclickRemoval->removalStarted = true;
            pKeyclickRemoval->limitcount = 0;
            return;
        }
        else if(pKeyclickRemoval->keyState == KEYDOWN && pKeyclickRemoval->keyId == KEY_NONE)  /* new key-up event */
        {
            #ifdef KEYCLICKREMOVAL_DEBUGPRINT
            printf("keyclickremoval_trigger: Key-up \n ");
            #endif
            pKeyclickRemoval->keyState   = KEYUP;
            pKeyclickRemoval->removalStarted = true;
            pKeyclickRemoval->limitcount = 0;
            return;
        }
        else if(pKeyclickRemoval->keyState == KEYDOWN && pKeyclickRemoval->keyId != KEY_NONE) /* error, new key-down event while key is down */
        {
            #ifdef KEYCLICKREMOVAL_DEBUGPRINT
            printf("keyclickremoval_trigger: Key-down event before key-up event \n ");
            #endif
            pKeyclickRemoval->keyState   = KEYDOWN;
            pKeyclickRemoval->removalStarted = true;
            pKeyclickRemoval->limitcount = 0;
            return;
        }
        else if(pKeyclickRemoval->keyState == KEYUP && pKeyclickRemoval->keyId == KEY_NONE)
        {
            #ifdef KEYCLICKREMOVAL_DEBUGPRINT
            printf("keyclickremoval_trigger: Two key-up events in a row \n");
            #endif
            pKeyclickRemoval->keyState = KEYUP;
            pKeyclickRemoval->limitcount = 0;
            return;
        }
        else
        {
            #ifdef KEYCLICKREMOVAL_DEBUGPRINT
            printf("keyclickremoval_trigger: Does this happen?\n");
            #endif
            return;
        }
}

/**
 * Service function that returns the status of the keyremovalroutine.
 * @param pKeyclickRemoval
 *
 * @return bool
 */
bool keyclickremoval_getState(KEYCLICKREMOVAL_PTR pKeyclickRemoval)
{
    return pKeyclickRemoval->removalStarted;
}

/**
 * Find keyId from the 48 bits following the 0xACDC mask.
 * We need to look through two uint32 samples.
 * The keyId is the one high bit found. If more high bits are found, there is an error.
 *
 * @param pKeyclickremoval
 * @param inbuf
 * @param sampleindex
 * @param bitstartindex
 */
static void keyclickremoval_findKeyid(KEYCLICKREMOVAL_PTR pKeyclickremoval, int * inbuf, int sampleindex)
{
    unsigned int buf;
    int i = 47;
    int sampleLoadindex = sampleindex + 1; /* need 2 samples to find keyId */
    int counter = 0;
    KEYCLICKREMOVAL_KeyId_t keyId = KEY_NONE;

    if( sampleLoadindex >= 480 )
    {
        #ifdef KEYCLICKREMOVAL_DEBUGPRINT
        printf("ACDC detection too late, in sample %d. Not able to find keyId if present \n", sampleindex);
        #endif
        pKeyclickremoval->timing = sampleindex;
        pKeyclickremoval->keyId  = KEY_DEFAULT;
        keyclickremoval_trigger(pKeyclickremoval);
        return;
    }
    buf = (unsigned int) inbuf[sampleLoadindex];

    while( i >= 0 )
    {
        if( (buf & 0x1) != 0x0 )
        {
            counter++;
            keyId = (KEYCLICKREMOVAL_KeyId_t) i;
        }

        buf >>= 1;

        i--;
        if( i == 15 )
        {
            buf = (unsigned int)inbuf[--sampleLoadindex];
        }
    }

//     printf("sample %d: 0x%x, sample %d: 0x%x, sample %d: 0x%x \n", sampleindex,
//            inbuf[sampleindex], sampleindex+1, inbuf[sampleindex+1], sampleindex+2, inbuf[sampleindex+2] );

    if( counter == 0 )
    {
        #ifdef KEYCLICKREMOVAL_DEBUGPRINT
        printf("keyclickremoval: KeyEvent verified in sample %d, no keyId \n", sampleindex);
        #endif
        keyId = KEY_NONE;
    }
    else if( counter > 1 )
    {
        #ifdef KEYCLICKREMOVAL_DEBUGPRINT
        printf("keyclickremoval: KeyEvent verified in sample %d, %d keyIds found, use default estimate \n ", sampleindex, counter);
        #endif
        keyId = KEY_DEFAULT;
    }
    else
    {
        #ifdef KEYCLICKREMOVAL_DEBUGPRINT
        printf("keyclickremoval: KeyEvent verified in sample %d, keyId: %d \n", sampleindex, keyId);
        #endif
    }

    pKeyclickremoval->timing = sampleindex;
    pKeyclickremoval->keyId  = keyId;
    keyclickremoval_trigger(pKeyclickremoval);

    return;
}

/**
 * A key event is identified when the fpga sends a 0xACDC bitmask.
 * Assumes alignment of the 0xACDC mask.
 *
 * @param pKeyclickremoval
 * @param inbuf
 * @param sampleindex
 */
void keyclickremoval_verifyEvent(KEYCLICKREMOVAL_PTR pKeyclickRemoval, int *inbuf, int sampleindex)
{
    #define ACDC_MASK 0xacdc
    /* check for acdc mask with start at bit position 0, from the most significant bit */

    if( pKeyclickRemoval->keyclickremovalOn ) /* No use in detecting keyEvents if not running the routine */
    {
        unsigned int buf = (unsigned int)inbuf[sampleindex];

        buf >>= 16;
        if( (buf ^ ACDC_MASK) == 0x0 )   /* bitwise exclusive OR */
        {
            //printf("keyclickremoval_verifyEvent: ACDC_mask verified in sample %d \n", sampleindex);
            keyclickremoval_findKeyid(pKeyclickRemoval, inbuf, sampleindex);
        }
        else
        {
            #ifdef KEYCLICKREMOVAL_DEBUGPRINT
            printf("keyclickremoval_verifyEvent: acdc_mask not found in sample %d \n", sampleindex);
            #endif
        }
    }
    return;
}

/********* test functions ************/

void keyclickremoval_setKeyclickprocess(KEYCLICKREMOVAL_PTR pKeyclickRemoval, bool onoff)
{
    const char *status[] = {"off", "on"};
    pKeyclickRemoval->keyclickremovalOn = onoff;
    printf("pKeyclickRemoval->keyclickremovalOn is set %s \n\n", status[pKeyclickRemoval->keyclickremovalOn]);
    return;
}

void keyclickremoval_setmaxlimitcount(KEYCLICKREMOVAL_PTR pKeyclickRemoval, int value)
{
    pKeyclickRemoval->limitcountmax = value;
    printf("pKeyclickRemoval->limitcountmax is set to %d \n\n", pKeyclickRemoval->limitcountmax);
    return;
}

void keyclickremoval_setgi(KEYCLICKREMOVAL_PTR pKeyclickRemoval, float value)
{
    pKeyclickRemoval->gi = value;
    printf("pKeyclickRemoval->gi is set to %f \n\n", pKeyclickRemoval->gi);
    return;
}

void keyclickremoval_setgd(KEYCLICKREMOVAL_PTR pKeyclickRemoval, float value)
{
    pKeyclickRemoval->gd = value;
    printf("pKeyclickRemoval->gd is set to %f \n\n", pKeyclickRemoval->gd);
    return;
}

void keyclickremoval_status(KEYCLICKREMOVAL_PTR pKeyclickRemoval)
{
    const char *status[] = {"off", "on"};

    printf("\rStatus: Keyclick");
    printf("\r\n   keyremovalOn      -  %s", status[pKeyclickRemoval->keyclickremovalOn]);
    printf("\r\n   limitcountmax     -  %d", pKeyclickRemoval->limitcountmax);
    printf("\r\n   limit gi          -  %f", pKeyclickRemoval->gi);
    printf("\r\n   limit gd          -  %f", pKeyclickRemoval->gd);
    printf("\r\n   limit lev         -  %f", pKeyclickRemoval->lev);
    printf("\r\n");
}

#ifdef UNITTEST
/***************************************************************************
 * UNITTEST_KEYCLICKREMOVAL
 *    Desc.:
 ***************************************************************************/
#include "unittest.h"
#include "testdata/keyclickremovaltestdata.h"  /* testvectors defined here */
#include <string.h>

/* Private testdata */
static float keyclickremovalTestInput[]  = {KEYCLICKREMOVAL_TESTINPUT_DEFINE};
static float keyclickremovalTestOutput[] = {KEYCLICKREMOVAL_TESTOUTPUT_DEFINE};

#define nWRITE2FILE

/***************************************************************************
 * COMPARERESULTS
 *    Desc.: Compares two float buffers. Due to compilator optimization, a
 *           a small difference 'epsilon' is allowed when comparing with the
 *           testvector
 ***************************************************************************/
static int compareResults(float *inp1, float *inp2, int len)
{
    int differ = 0;
    int i;
    float epsilon = 3.0e-7f;

    for( i = 0; i < len; i++ )
    {

        if( (fabsf(inp1[i] - inp2[i])) > epsilon )
        {
            differ = 1;
            printf("keyclickremoval_unittest: too large difference: %e \n", (fabsf(inp1[i] - inp2[i])));
        }
    }
    return differ;
}


/***************************************************************************
 * UNITTEST_KEYCLICKREMOVAL
 **************************************************************************/
void unittest_keyclickremoval(void)
{
    KEYCLICKREMOVAL_PTR pKeyclickRemoval;
    float *micbuf;
    char const testName[] = {"Keyremoval"};
    int i;
    int differ = 0; //criteria flag to be returned by function. 0: OK, 1: ERROR, 2:Nothing done, hence ERROR
    int numFrames = 7;
    int startpress = 1;


#ifdef WRITE2FILE
    float *outbuf;
    int cnt;
    FILE * fid;
    const char * Filename = "/home/jps/src/testdata/keyclickremovaltestoutput.dat";
#endif

    pKeyclickRemoval = keyclickremoval_create();
    if( pKeyclickRemoval == NULL )
    {
        fprintf(stderr, "unittest_keyclickremoval: memory allocation failed\n");
        unittest_assert(0);
        return;
    }

    unittest_context(testName);
    keyclickremoval_init(pKeyclickRemoval);

    pKeyclickRemoval->keyclickremovalOn = true;

    //Allocating memory buffers
    micbuf = (float *) calloc(FRAMESIZE, sizeof(float));

#ifdef WRITE2FILE
    outbuf =(float *) calloc(FRAMESIZE*numFrames,sizeof(float));
#endif
    for( i = 0; i < numFrames; i++ )
    {
        scratchmem_lock();
        //printf("frame: %d\n", i);
        memcpy(micbuf, keyclickremovalTestInput + i * FRAMESIZE, sizeof(float) * FRAMESIZE);

        if( i == startpress )  //key pressed at this frame
        {
            pKeyclickRemoval->timing = 0;
            pKeyclickRemoval->keyId = KEY_UNITTEST;
            keyclickremoval_trigger(pKeyclickRemoval);
        }
        keyclickremoval_limiter_process(pKeyclickRemoval, micbuf);

        differ += compareResults(micbuf, keyclickremovalTestOutput + i * FRAMESIZE, FRAMESIZE);

#ifdef WRITE2FILE
        memcpy(outbuf+i*FRAMESIZE, micbuf, FRAMESIZE*sizeof(float));
#endif
        scratchmem_unlock();
    }

#ifdef WRITE2FILE
    fid = (void*)fopen(Filename,"wb");
    cnt = fwrite(outbuf,sizeof(float),FRAMESIZE*numFrames,fid);
    printf("unittest_keyclickremoval: write output to file, %d floats \n", cnt);
    fclose(fid);
    free(outbuf);
    (void)cnt;
#endif

    unittest_assert(differ == 0);

    free(micbuf);
    keyclickremoval_destroy(pKeyclickRemoval);
}

#endif /* UNITTEST */
