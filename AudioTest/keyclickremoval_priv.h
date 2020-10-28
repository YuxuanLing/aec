#ifndef KEYCLICKREMOVAL_PRIV_H
#define KEYCLICKREMOVAL_PRIV_H

#include "audec_defs.h"
#include <stdbool.h>

/* number is given after bit position in the 48 bits identifying a key */
typedef enum
{
    KEY_NONE = 0,

    KEY_HANDSFREE = 2,
    KEY_HEADSET = 4,
    KEY_VIDEO = 5,
    KEY_SELFVIEW = 9,
    KEY_MIC = 12,
    KEY_VOLPLUSS = 22,
    KEY_VOLMINUS = 23,
    KEY_PC = 13,
    KEY_LETTER = 14,
    KEY_PHONEBOOK = 15,

    KEY_HOME = 25,
    KEY_ARROWDOWN = 26,
    KEY_ARROWUP = 33,
    KEY_ARROWLEFT = 34,
    KEY_ARROWRIGHT = 41,
    KEY_OK = 42,

    KEY_F1 = 7,          /* softbutton 1 */
    KEY_F2 = 3,          /* softbutton 2 */
    KEY_F3 = 1,          /* softbutton 3 */
    KEY_F4 = 17,         /* softbutton 4 */
    KEY_F5 = 6,          /* softbutton 5 */
    KEY_C = 39,
    KEY_1 = 45,
    KEY_2 = 37,
    KEY_3 = 29,
    KEY_4 = 44,
    KEY_5 = 36,
    KEY_6 = 28,
    KEY_7 = 43,
    KEY_8 = 35,
    KEY_9 = 27,
    KEY_0 = 38,
    KEY_HASH = 30,
    KEY_STAR = 46,
    KEY_RED = 31,
    KEY_GREEN = 47,

    KEY_DEFAULT,
    KEY_UNITTEST   /* for unittest */

} KEYCLICKREMOVAL_KeyId_t;


typedef enum
{
    KEYUP,
    KEYDOWN
}
KEYCLICKREMOVAL_State_t;

typedef struct KEYCLICKREMOVAL {
    KEYCLICKREMOVAL_KeyId_t keyId;
    int timing;
    KEYCLICKREMOVAL_State_t keyState;

    bool keyclickremovalOn;
    bool removalStarted;
    bool limit;
    int limitcount;
    int limitcountmax;
    float lim_gain;
    float gi;
    float gd;
    float lev;
    float tc;
} KEYCLICKREMOVAL;

#endif //KEYCLICKREMOVAL_PRIV_H
