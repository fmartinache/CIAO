// =========================================================================
//                 Shared memory data structure definition
// =========================================================================

#ifndef _SHM_DS_H
#define _SHM_DS_H


#define CHAR 1
#define INT 2
#define FLOAT 3
#define DOUBLE 4
#define COMPLEX_FLOAT 5
#define COMPLEX_DOUBLE 6
#define USHORT 7




typedef struct {
  float re;
  float im;
} complex_float;


typedef struct {
  double re;
  double im;
} complex_double;

typedef struct {  // ======= KEYWORD structure ==========
  char name[16];
  char type; // N: unused, L: long, D: double, S: 16-char string 
  union {
    long numl;
    double numf;
    char valstr[16];
  } value;
  char comment[80];
} SHM_KWD;


typedef struct { // ========= METADATA structure ============
  char name[80];                // image name
  long naxis;                   // number of axis
  long size[3];                 // image size 
  long nel;		        // number of elements in image
  int atype;			// data type code   

  double creation_time;	        // creation time (since program start)
  double last_access;		// last access time (since program start)
  struct timespec wtime;

  int shared;                   // 1 if in shared memory
  int write;                    // 1 if image is being written  
  int status;
  long cnt0;                    // counter (if image is updated)
  long cnt1;
  
  long nbKW;                    // number of keywords

} SHM_METADATA;

typedef struct { // ========= DATA ARRAY structure =================
  SHM_METADATA *md;            // 1. a METADATA header
  union {
    char           *C;
    int            *I;
    float          *F;
    double         *D;
    complex_float  *CF;
    complex_double *CD;
    unsigned short *U;
  } array;                      // 2. a pointer to data array
  SHM_KWD *kw;                  // 3. a pointer to a series of keywords
} SHM_DATA;

#endif
