#ifndef FFTUTIL_H
#define FFTUTIL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __ANDROID__
# define complex     _Complex
# define _Complex_I  (__extension__ 1.0iF)
# define I           _Complex_I
#else
# include <complex.h>
#endif

//#include "complex.h"
void partial_fft(int N, int M, float z[]);
void complete_fft(const int N, float zi[], float zo[]);

#endif
