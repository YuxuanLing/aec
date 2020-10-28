#ifndef AUDTYPES_H
#define AUDTYPES_H

#include "auddefs.h"

#include <inttypes.h>

/** 
 * Use of filterbank type
 */
typedef enum {
    FILTERBANK_NONE = 0,
    FILTERBANK_EC,
    FILTERBANK_ZDNLP
} FILTERBANK_USE_TYPE;

/** 
 * Audio physical connector type
 */
typedef enum
{
    AL_CONNECTORVAR_XLR,
    AL_CONNECTORVAR_RCA,
    AL_CONNECTORVAR_HDMI,
    AL_CONNECTORVAR_TLINK,
    AL_CONNECTORVAR_ALERTSPKR,
    AL_CONNECTORVAR_INTMICSPKR,
    AL_CONNECTORVAR_HEADSET,
    AL_CONNECTORVAR_HANDSET,
    AL_CONNECTORVAR_BLUETOOTH,
    AL_CONNECTORVAR_SNOOPYFPGA,
    AL_CONNECTORVAR_DIT,
    AL_CONNECTORVAR_INTSPKR,
    AL_CONNECTORVAR_MIC,
    AL_CONNECTORVAR_LINE,
    AL_CONNECTORVAR_UNDEF,
    AL_CONNECTORVAR_CNT=AL_CONNECTORVAR_UNDEF
} AL_ConnectorVariant_t;

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

