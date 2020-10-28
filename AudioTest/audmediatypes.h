#ifndef __AUDMEDIATYPES_H__
#define __AUDMEDIATYPES_H__

#define TYPE_END_MARKER         0x7FFFFFFF

typedef enum
{
    AUDCSTD_OFF,               /* 00 */
    AUDCSTD_UNCOMP8,           /* 01 */     /*Uncompressed */
    AUDCSTD_G711MU_48,         /* 02 */
    AUDCSTD_G711MU_56,         /* 03 */
    AUDCSTD_G711MU_64,         /* 04 */
    AUDCSTD_G711A_48,          /* 05 */
    AUDCSTD_G711A_56,          /* 06 */
    AUDCSTD_G711A_64,          /* 07 */
    AUDCSTD_G722_48,           /* 08 */
    AUDCSTD_G722_56,           /* 09 */
    AUDCSTD_G722_64,           /* 10 */
    AUDCSTD_G722_1_16,         /* 11 */
    AUDCSTD_G722_1_24,         /* 12 */
    AUDCSTD_G722_1_32,         /* 13 */
    AUDCSTD_G722_1_24_32kHz,   /* 14 */
    AUDCSTD_G722_1_32_32kHz,   /* 15 */
    AUDCSTD_G722_1_48_32kHz,   /* 16 */
    AUDCSTD_G728,              /* 17 */
    AUDCSTD_G723_1,            /* 18 */
    AUDCSTD_G729,              /* 19 */
    AUDCSTD_G729A,             /* 20 */
    AUDCSTD_G729AB,            /* 21 */
    AUDCSTD_AACLD_48,          /* 22 */
    AUDCSTD_AACLD_56,          /* 23 */
    AUDCSTD_AACLD_64,          /* 24 */
    AUDCSTD_AACLD_128,         /* 25 */
    AUDCSTD_CISCO_PCM16_256,   /* 26 */ /* Cisco L16bit linear pcm at 16Khz sample rate */
    AUDCSTD_COMFORT_NOISE,     /* 27 */
    AUDCSTD_AACLC_48,          /* 28 */
    AUDCSTD_AACLC_56,          /* 29 */
    AUDCSTD_AACLC_64,          /* 30 */
    AUDCSTD_AACLC_96,          /* 31 */
    AUDCSTD_AACLC_128,         /* 32 */
    AUDCSTD_RFC2833EVENT,      /* 33 */
    AUDCSTD_L16_768,           /* 34 */ /* L16bit linear pcm at 48Khz sample rate */
    AUDCSTD_UNDEF,  /* undefined std */
    AUD_STD_TYPE_LAST,
    AUD_STD_TYPE_UNDEF,
    AUD_STD_TYPE_END = TYPE_END_MARKER
} AUD_STD_TYPE;


typedef enum
{
    AF_off,
    AF_G711mu,
    AF_G711A,
    AF_G722,
    AF_G722_1,
    AF_G728,
    AF_G729,
    AF_G729A,
    AF_G729AB,
    AF_PCM_Float,
    AF_PCM_Fixed,
    AF_AAC,
    AF_Last,      /* AF_Last must be the last element, OK???
                   * (It will be used to find the count of enumerations)*/

    AUD_FORMAT_TYPE_LAST,
    AUD_FORMAT_TYPE_UNDEF,
    AUD_FORMAT_TYPE_END = TYPE_END_MARKER
} AUD_FORMAT_TYPE;

#endif
