/************************************************************************
  
     Integral Systems, Inc.
  
***********************************************************************
  
     Project   : Operations Ground Equipment for GOES
     System    : Sensors Processing System
     Source    : cilpevsc.h
     Programmer: Igor Levine
  
     Ver.    Data    By   Comment
     ----  --------  ---  ---------------------------------------------
     1     11/28/88  IL   Initial creation
     2     03/02/93  IL   Conversion from FORTRAN to C
     3     05/17/95  IL   Deleted obsolete variables from CILPEVSC.
***********************************************************************/
#include <stdio.h>
#include <math.h>

#ifndef CILPEVSC_H
#define CILPEVSC_H

#ifndef  TRUE
#define  TRUE    1
#endif

#ifndef  FALSE
#define  FALSE   0
#endif

#ifndef  SQ
#define  SQ(a)  ((a)*(a))
#endif

#ifndef  max
#define  max(a,b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef  min
#define  min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef  sign
#define  sign(a)   (((a) < 0) ? -1 : 1)
#endif

#ifndef  ABS
#define  ABS(a)  (((a) >= 0) ? (a) : -(a))    
#endif

#ifndef  NINT
#define  NINT(a)  ( (a) + (((a) < 0) ? -0.5 : 0.5))    
#endif

#ifndef  QUATER
typedef  union {
         float         r4;
         long          i4;
         unsigned long u4;
         char          c4;
         }     QUATER;
#endif


typedef struct {
        double  nomorb;    /* Nominal distance = 42164.365 km */
        double  pi;
        double  pi2;
        double  rad;       /* Degrees to radians conversion  */
        double  ae;        /* Earth equatorial radius (km) */
        double  be;        /* Earth polar radius (km) */
        double  aebe;      /*  ae/be  */
        double  aebe2;     /* (ae/be)**2 */  
        double  aebe3;     /* (ae/be)**2 - 1 */
        double  aebe4;     /* (be/ae)**4 - 1 */
        double  xs[3];     /* Normalized s/c position in ECEF coordinates */
        double  b[3][3];   /* ECEF to s/c coordinates transform. matrix */
        double  bt[3][3];  /* ECEF to instrument transformation matrix */
        double  q3;        /* used in function lpoint() */
        double  pitch;     /* pitch angle of instrument (rad) */
        double  roll;      /* roll angle of instrument (rad) */
        double  yaw;       /* yaw angle of instrument (rad) */
        double  pma;       /* pitch misalignments (rad) */
        double  rma;       /* roll misalignment (rad) */
        }    SPSAREA;

typedef struct {
        int    iymax[2];     /* Increments per cycle in elevation angle  */
        int    ixmax[2];     /* Increments per cycle in scan angle */
        float  alfmax[2];    /* Bounds in elevation (rad) */
        float  zetmax[2];    /* Bounds in scan angle (rad) */
        float  alfincr[2];   /* Change in elev. angle per increment (rad)  */
        float  zetincr[2];   /* Change in scan angle per increment (rad) */
        float  incrln[2];    /* # of increments per detector scan line */
        float  incrpx[2];    /* # of increments per pixel */
        float  alfln[2];     /* Elev. angle per detector scan line (rad)  */
        float  zetpx[2];     /* Scan angle per pixel (rad) */
        }      CILPEVSC;

#endif
