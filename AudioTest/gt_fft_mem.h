#ifndef GT_FFT_MEM_H
#define GT_FFT_MEM_H

/* WARNING: DON'T CHANGE THE LAYOUT OF THESE STRUCT TYPES 
 * WITHOUT CHANGING THE INITIALIZATION WHERE THEY ARE USED */

typedef struct
{
    unsigned char brev[64];
    float         T_16[32];
    unsigned char ktab240[240];
    unsigned char itab240[240];
    unsigned char ntab240[240];
}GT240TABLE_STRUCT;

typedef struct
{
    float         T_16[32];
    unsigned char ktabnums240[8];
    unsigned char ktabspec240[31 + 31 + 31 + 31 + 31 + 31 + 31 + 16];
    unsigned char ntabnums240[6];
    unsigned char ntabspec240[33 + 17 + 33 + 33 + 33 + 17];
    unsigned char itabnums240[3];
    unsigned char itabspec240[123 + 115 + 4];
    unsigned char nswaps240[2*38];
}GT240IPTABLE_STRUCT;

/* FIXME: twiddle factors are generated in fft.c */
extern float w_half[256];

#endif

