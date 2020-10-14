// Funtions for DSO9254A oscilliscope
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

/* Defines */
#define MAX_LENGTH 10000000
#define IO_TIMEOUT 20000

/* Type definitions */
typedef unsigned long int UINT64; /* This defines a 64-bit unsigned
                                   * integer for Microsoft platforms.
                                   */

/* Structure and Union definitions */
union DATATYPE
{
   char buffer[MAX_LENGTH]; /* Buffer for reading word format data */
   char byte[MAX_LENGTH];
   unsigned short word[MAX_LENGTH / 2];
   UINT64 longlong[MAX_LENGTH / 4];
};

typedef struct
{
   char Cookie[2];
   char Version[2];
   int FileSize;
   int NumberOfWaveforms;
} FileHeader;
FileHeader fileHeader;
   
const char COOKIE[2] = {'A', 'G'};
const char VERSION[2] = {'1', '0'};

#define DATE_TIME_STRING_LENGTH 16
#define FRAME_STRING_LENGTH 24
#define SIGNAL_STRING_LENGTH 16

typedef struct
{
   int HeaderSize;
   int WaveformType;
   int NWaveformBuffers;
   int Points;
   int Count;
   float XDisplayRange;
   double XDisplayOrigin;
   double XIncrement;
   double XOrigin;
   int XUnits;
   int YUnits;
   char Date[DATE_TIME_STRING_LENGTH];
   char Time[DATE_TIME_STRING_LENGTH];
   char Frame[FRAME_STRING_LENGTH];
   char WaveformLabel[SIGNAL_STRING_LENGTH];
   double TimeTag;
   unsigned int SegmentIndex;
} WaveformHeader;
WaveformHeader waveformHeader;

typedef struct
{
   int HeaderSize;
   short BufferType;
   short BytesPerPoint;
   int BufferSize;
} WaveformDataHeader;

WaveformDataHeader waveformDataHeader;

typedef enum
{
   PB_UNKNOWN,
   PB_NORMAL,
   PB_PEAK_DETECT,
   PB_AVERAGE,
   PB_HORZ_HISTOGRAM,
   PB_VERT_HISTOGRAM,
   PB_LOGIC
} WaveformType;

typedef enum
{
   PB_DATA_UNKNOWN,
   PB_DATA_NORMAL,
   PB_DATA_MAX,
   PB_DATA_MIN,
   PB_DATA_TIME,
   PB_DATA_COUNTS,
   PB_DATA_LOGIC
} DataType;

/* Prototypes */
void OutputNormalWaveform(WaveformHeader waveformHeader);
/* Globals */
double xOrg = 0L, xInc = 0L; /* Values necessary to create time data */
union DATATYPE WaveFormData; /* Used to input and output data */
FILE *InputFile = NULL;
char *buffer;
float Volts[MAX_LENGTH];
