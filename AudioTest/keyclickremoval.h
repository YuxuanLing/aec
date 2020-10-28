#ifndef KEYCLICKREMOVAL_H
#define KEYCLICKREMOVAL_H

#include <stdbool.h>


typedef struct KEYCLICKREMOVAL * KEYCLICKREMOVAL_PTR;


KEYCLICKREMOVAL_PTR keyclickremoval_create(void);
void keyclickremoval_destroy(KEYCLICKREMOVAL_PTR pKeyclickRemoval);
void keyclickremoval_init(KEYCLICKREMOVAL_PTR pKeyclickRemoval);
void keyclickremoval_process(KEYCLICKREMOVAL_PTR pKeyclickRemoval, float *micfft, float * gain);
void keyclickremoval_limiter_process(KEYCLICKREMOVAL_PTR pKeyclickRemoval, float *micbuf);
void keyclickremoval_verifyEvent(KEYCLICKREMOVAL_PTR pKeyclickRemoval, int *inbuf, int sampleindex);

void keyclickremoval_setclip(KEYCLICKREMOVAL_PTR pKeyclickRemoval, bool onoff);
bool keyclickremoval_getState(KEYCLICKREMOVAL_PTR pKeyclickRemoval);

/* test functions */
void keyclickremoval_setKeyclickprocess(KEYCLICKREMOVAL_PTR pKeyclickRemoval, bool onoff);
void keyclickremoval_setmute(KEYCLICKREMOVAL_PTR pKeyclickRemoval, bool onoff);
void keyclickremoval_setlength(KEYCLICKREMOVAL_PTR pKeyclickRemoval, int value);
void keyclickremoval_setmaxmutecount(KEYCLICKREMOVAL_PTR pKeyclickRemoval, int value);
void keyclickremoval_setmaxlimitcount(KEYCLICKREMOVAL_PTR pKeyclickRemoval, int value);
void keyclickremoval_status(KEYCLICKREMOVAL_PTR pKeyclickRemoval);
void keyclickremoval_setgi(KEYCLICKREMOVAL_PTR pKeyclickRemoval, float value);
void keyclickremoval_setgd(KEYCLICKREMOVAL_PTR pKeyclickRemoval, float value);


#ifdef UNITTEST
void unittest_keyclickremoval(void);
#endif //UNITTEST

#endif //KEYCLICKREMOVAL_H
