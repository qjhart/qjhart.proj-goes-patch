/******************************************************************************
 * $Id: $
 *
 * Project:  PROJ.4
 * Purpose:  Implementation of Specialized GOES projection.
 * Author:   Quinn Hart
 *
 ******************************************************************************
 * Copyright (c) 2004, Quinn Hart
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************
 *
 * $Log: $
 *
 */
#include <string.h>
#include "geocent.h"
#include "PJ_goes_earthloc.p"

#define PROJ_PARMS__ \
int goes; \
int instr; \
int row_col; \
int nscyc; \
int nsinc; \
int ewcyc; \
int ewinc; \
double    lam;\
double    dr ;\
double    phi;\
double    psi;\
double    roll;\
double    pitch;\
double    yaw;\
float oa[336];\
double q3;\
double xs[3]; \
double bt[9]; \
double elvln;\
double scnpx;\
double ewnom;\
double elvmax;\
double scnmax;\
double D; \
double r; \
 int noaa; \
 int imc; \
double tu; \
double t;\
int flip; \
SPSAREA spsarea;\
CILPEVSC cilpevsc;\
float rlat; \
float rlon;

#define PJ_LIB__
#include <projects.h>

#define EPS10	1.e-10
#define TOL7	1.e-7

PROJ_HEAD(goes, "Geostationary Operational Environmental Satellites (GOES)") "\n\tMisc&Ell\n\tgoes= row_col= ...";

static XY e_pre11_forward(LP lp, PJ *P) { XY xy = {0.0,0.0};
/* long=lp.lam, lat=lp.phi */
 projUV se;			/* se.u=scan se.v=elev */

 {
   /* Local variables */
   static double doff;
   static double f[3];
   static double u[3], alpha1, w1, w2, ft[3];
   
   /* Fill u with X,Y,Z (in AE units) */
   static GeocentricInfo gi;
   pj_Set_Geocentric_Parameters(&gi,(double)1.0,(double)sqrt(P->one_es));
   pj_Convert_Geodetic_To_Geocentric(&gi,lp.phi,lp.lam,0.0,u,u+1,u+2);   
   /*     POINTING VECTOR FROM SATELLITE TO THE EARTH POINT */
   
   f[0] = u[0] - P->xs[0];
   f[1] = u[1] - P->xs[1];
   f[2] = u[2] - P->xs[2];
   w2 = u[0] * (float) f[0] + u[1] * (float) f[1] + 
	 u[2] * (float) f[2] * P->rone_es;
   
   /*     VERIFIES VISIBILITY OF THE POINT */
   
   if (w2 > (float)0.) F_ERROR;
   
   /*    CONVERTS POINTING VECTOR TO INSTRUMENT COORDINATES (FT=BT*F)*/
   
   ft[0] = P->bt[0] * f[0] + P->bt[1] * f[1] + P->bt[2] * 
	 f[2];
   ft[1] = P->bt[3] * f[0] + P->bt[4] * f[1] + P->bt[5] * 
	 f[2];
   ft[2] = P->bt[6] * f[0] + P->bt[7] * f[1] + P->bt[8] * 
	 f[2];
   
   /*     CONVERTS POINTING VECTOR TO SCAN AND ELEVATION ANGLES */
   
   /* When IMC=0 There is an additional correction here ELUG 4.6 */
   
   /* Compute SIGN OF MISALIGNMENT CORRECTIONS AND ORIGIN OFFSET CORRECTIONS */
   doff = P->scnmax - P->ewnom;
   /* Optical Axis Correction ELUG 3..8 */
   se.u = atan(ft[0] / sqrt(ft[1]*ft[1] + ft[2]*ft[2]));
   se.v = -atan(ft[1] / ft[2]);
   //printf("qjh: %f %f\n",se.u,se.v);
   w1 = sin(se.v);
   w2 = cos(se.u);
   alpha1 = se.v ;
   //printf("qjh: %f %f\n",se.u,se.v);
   se.v = alpha1 + alpha1 * se.u * doff;
   se.u -= alpha1 * (float).5 * alpha1 * doff;
   //printf("qjh: %f %f\n",se.u,se.v);
   
 } /* gpoint_ */
 
 if(P->row_col != 0) {
   int row_sign = (P->row_col < 0)?-1:1;
   /* Convert Elev,Scan to XY */ /*x=pix y=rline*/
   xy.x = (P->scnmax + se.u) / P->scnpx + 1.0;
   xy.y = row_sign * ((P->elvmax - se.v) / P->elvln + P->D );
   xy.x = ((xy.x / P->fr_meter) - P->x0) / P->a;
   xy.y = ((xy.y / P->fr_meter) - P->y0) / P->a;
 } else {
   xy.x=atan(se.u)*(P->r-1.0);
   xy.y=atan(se.v)*(P->r-1.0);
 }
 return (xy);
}

static XY e_post11_forward(LP lp, PJ *P) {
  XY xy = {0.0,0.0};  
  float elev_y,scan_x;
  int ierr;
  
  /* long=lp.lam, lat=lp.phi */
  /* NOAA uses a different instr numbering */
  gpoint(P->instr+1, P->flip, lp.phi, lp.lam, &elev_y, &scan_x, 
     &ierr, &P->spsarea, &P->cilpevsc);
  
  if (ierr != 0 ) F_ERROR;
  
  if(P->row_col != 0) {
    int row_sign = (P->row_col < 0)?-1:1;
#if 1
    evsc2lp_double(1, 1, elev_y, scan_x, &xy.y,&xy.x, &P->cilpevsc);
#else
    int iscan,idet,line,ipix;
    float offln,offpx;
    lnpxevsc(P->instr, P->flip, elev_y, scan_x,
	     &iscan,&idet,&line,&ipix,
	     &offln, &offpx, &P->cilpevsc);    
      xy.x = ipix+offpx;
      xy.y = line+offln;
#endif
    xy.y *= row_sign;
    xy.x = ((xy.x / P->fr_meter) - P->x0) / P->a;
    xy.y = ((xy.y / P->fr_meter) - P->y0) / P->a;
  } else {
    xy.x=atan(scan_x)*(P->r-1.0);
    xy.y=atan(elev_y)*(P->r-1.0);
  }
  return (xy);
}

static LP e_post11_inverse(XY xy, PJ *P) { 
  LP lp = {0.0,0.0};  /* long=lp.lam lat=lp.phi */
  double elev_y, scan_x;
  projUV se;
  //  float se_u,se_v;			/* Since NOAA code is floats */
  float lat,lon;
  int ierr;
  
  if(P->row_col != 0) {
    /* Convert XY to Elev,Scan */
    /* x=pix y=rline */
    int row_sign = (P->row_col < 0)?-1:1;
    /* Fix back to Row-Col */
    scan_x = P->fr_meter * (P->a * xy.x + P->x0);
    elev_y = P->fr_meter * (P->a * xy.y + P->y0);
    elev_y *= row_sign;

  /* NOAA uses a different instr numbering */
	se.u = scpx(P->instr+1,P->flip,scan_x,&P->cilpevsc);
	se.v = evln(P->instr+1,P->flip,elev_y,&P->cilpevsc);
    } else {
    /* This is to match PJ_geos */
    se.u = tan (xy.x / (P->r-1.0));
    se.v = tan (xy.y / (P->r-1.0));
  }
  lpoint(P->instr+1,P->flip,se.v,se.u, &lat, &lon, &ierr, 
	 &P->spsarea,&P->cilpevsc);
  //ierr 0=point_on_earth, 1= instrument off earth
  if (ierr != 0 ) I_ERROR;
  lp.lam = lon;
  lp.phi = lat;
  return(lp);
}

/* //FORWARD(e_forward); is */
/* static XY e_forward(LP lp, PJ *P) { */
/*   if (P->goes < 10)  */
/*     return e_pre11_forward(lp,P); */
/*   else  */
/*     return e_post11_forward(lp,P); */
/* } */


static LP e_pre11_inverse(XY xy, PJ *P) { 
  LP lp = {0.0,0.0};  /* long=lp.lam lat=lp.phi */
  projUV se;			/* se.u=scan se.v=elev */
    
  if(P->row_col != 0) {
    /* Convert XY to Elev,Scan */
    {				/* x=pix y=rline */
      int row_sign = (P->row_col < 0)?-1:1;
      /* Fix back to Row-Col */
      xy.x = P->fr_meter * (P->a * xy.x + P->x0);
      xy.y = P->fr_meter * (P->a * xy.y + P->y0);
      se.u = (xy.x - 1.0) * P->scnpx - P->scnmax;
      se.v = P->elvmax - (row_sign * xy.y - P->D) * P->elvln;
    }
    
  } else {
    /* This is to match PJ_geos */
    se.u = tan (xy.x / (P->r-1.0));
    se.v = tan (xy.y / (P->r-1.0));
  }
  
  /* Now Lpoint Elev, Scan to LP */
  {
    /* Local variables */
    static double doff, zeta;
    static double d;
    static double g[3];
    static double h;
    static double u[3], alpha;
    static double g1[3], q1, q2;
    static double ca, da, sa, cz, dz;
	
	
    /*      COMPUTE SIGN OF MISALIGNMENT CORRECTIONS AND ORIGIN OFFSET CORRECTIONS */
    doff = P->scnmax - P->ewnom;
    
    /*     ADD THE NEW SECOND ORDER ORIGIN OFFSET CORRECTION  ELUG 3.8*/
    
    alpha = se.v - se.v * se.u * doff;
    zeta = se.u + se.v * (float).5 * se.v * doff;
    
    /*     COMPUTES TRIGONOMETRIC FUNCTIONS OF THE SCAN AND ELEVATION */
    /*     ANGLES CORRECTED FOR THE ROLL AND PITCH MISALIGNMENTS */
    
    ca = cos(alpha);
    sa = sin(alpha);
    cz = cos(zeta);
    da = alpha;
    dz = zeta;
	
	/*    COMPUTES POINTING VECTOR IN INSTRUMENT COORDINATES */
	
    cz = cos(dz);
    g[0] = sin(dz);
    g[1] = -cz * sin(da);
    g[2] = cz * cos(da);
	
	/*     TRANSFORMS THE POINTING VECTOR TO EARTH-FIXED COORDINATES */
	
    g1[0] = P->bt[0] * g[0] + P->bt[3] * g[1] + P->bt[6] * 
	  g[2];
    g1[1] = P->bt[1] * g[0] + P->bt[4] * g[1] + P->bt[7] * 
	  g[2];
    g1[2] = P->bt[2] * g[0] + P->bt[5] * g[1] + P->bt[8] * 
	  g[2];
	
	/*     COMPUTES COEFFICIENTS AND SOLVES A QUADRATIC EQUATION TO */
	/*     FIND THE INTERSECT OF THE POINTING VECTOR WITH THE EARTH */
	/*     SURFACE */
	
    q1 = g1[0]*g1[0] + g1[1]*g1[1] + g1[2]*g1[2] * P->rone_es;
    q2 = P->xs[0] * g1[0] + P->xs[1] * g1[1] + P->xs[2] * P->rone_es * g1[2];
    d = (q2 * q2) - (q1 * P->q3);
    if (fabs(d) < 1e-9) {
	  d = (float)0.;
    }

	/*     IF THE DISCIMINANTE OF THE EQUATION, D, IS NEGATIVE, THE */
	/*     INSTRUMENT POINTS OFF THE EARTH */
	
    if (d < 0.) I_ERROR;
    d = sqrt(d);
	
	/*     SLANT DISTANCE FROM THE SATELLITE TO THE EARTH POINT */
	
    h = -(q2 + d) / q1;
	
	/*     CARTESIAN COORDINATES OF THE EARTH POINT */
	
    u[0] = P->xs[0] + h * g1[0];
    u[1] = P->xs[1] + h * g1[1];
    u[2] = P->xs[2] + h * g1[2];

    /* Fill u with X,Y,Z (in AE units) */
    static GeocentricInfo gi;
    pj_Set_Geocentric_Parameters(&gi,(double)1.0,(double)sqrt(P->one_es));
    pj_Convert_Geocentric_To_Geodetic(&gi,u[0],u[1],u[2],&lp.phi,&lp.lam,&h);
	
	return (lp);
  }
}

/* static LP e_inverse(XY xy, PJ *P) {  */
/*   if (P->goes < 10)  */
/*     return e_pre11_inverse(xy,P); */
/*   else  */
/*     return e_pre11_inverse(xy,P); */
/* } */


// FREEUP; 
static void freeup(PJ *P) {
  if (P) { pj_dalloc(P); }
}
/* This section for Row/Col variables It calculates the Row and Column
   scan and increments using the parameters for the instrument, and
   the NADIR pointing location
*/

static PJ * setup_pre11(PJ *P) {
   /* Local variables */
   static double slat, syaw, b[9]/* was [3][3] */, u,
	 ca, ci, sa, cu, si, su;
   static double asc;
   
   
   /* CONVERSION OF THE IMC LONGITUDE AND ORBIT YAW TO THE SUBSATELLITE */
   /* LONGITUDE AND THE ORBIT INCLINATION (REF: GOES-PCC-TM-2473, INPUTS*/
   /* REQUIRED FOR EARTH LOCATION AND GRIDDING BY SPS,  JUNE 6, 1988) */
   
   slat = sin(P->phi);
   syaw = sin(P->psi);
   si = slat*slat + syaw*syaw;
   ci = sqrt(1. - si);
   si = sqrt(si);
   if (slat == 0. && syaw == 0.) {
	 u = 0.;
   } else {
	 u = atan2(slat, syaw);
   }
   su = sin(u);
   cu = cos(u);
   
   /*     COMPUTES LONGITUDE OF THE ASCENDING NODE */
   
   asc = P->lam - u;
   sa = sin(asc);
   ca = cos(asc);
   
   /*     COMPUTES THE SUBSATELLITE GEOGRAPHIC LATITUDE */
   
   /* slat = atan(tan(P->phi) * P->rone_es); */
   
   /*     COMPUTES THE SUBSATELLITE LONGITUDE */
   
   /* slon = asc + atan2(ci * su, cu); */
   
   /*     COMPUTES THE SPACECRAFT TO EARTH FIXED COORDINATES TRANSFORMATION */
   /*     MATRIX: */
   /*         (VECTOR IN ECEF COORDINATES) = B * (VECTOR IN S/C COORDINATES) */
   
   b[3] = -sa * si;
   b[4] = ca * si;
   b[5] = -ci;
   b[6] = -ca * cu + sa * su * ci;
   b[7] = -sa * cu - ca * su * ci;
   b[8] = -slat;
   b[0] = -ca * su - sa * cu * ci;
   b[1] = -sa * su + ca * cu * ci;
   b[2] = cu * si;
   
   /*     COMPUTES THE NORMALIZED SPACECRAFT POSITION VECTOR IN EARTH FIXED */
   /*     COORDINATES - XS. */
   
   /* Fill ps  with X,Y,Z (in meters) of the satellite */
   /* This is a check to verivfy that the calculations are correct.
	  You can compare these results with the calculation below */
   // pj_Set_Geocentric_Parameters((double)1,(double)sqrt(P->one_es));
   //pj_Convert_Geodetic_To_Geocentric(P->phi,P->lam,P->r-1.0,
   //								 P->xs,P->xs+1,P->xs+2);
   
   P->xs[0] = -b[6] * P->r;
   P->xs[1] = -b[7] * P->r;
   P->xs[2] = -b[8] * P->r;
   
   /*     PRECOMPUTES Q3 (USED IN LPOINT) */
   P->q3 = P->xs[0]*P->xs[0] + P->xs[1]*P->xs[1] + P->xs[2]*P->xs[2] *P->rone_es - 1.;
   
   /*************************************************************************/
   /***   INST2ER ACCEPTS THE SINGLE PRECISION ROLL, PITCH AND YAW ANGLES */
   /***   OF AN INSTRUMENT AND RETURNS THE DOUBLE PRECISION INSTRUMENT TO */
   /***   EARTH COORDINATES TRANSFORMATION MATRIX. */
   /*************************************************************************/
   if(P->pitch != 0 || P->roll != 0) {
	 static int i, j;
	 static double cp, cr, cy, sp, sr, sy;
	 static double rpy[9]	/* was [3][3] */;
	 
	 sp = sin(P->pitch);
	 cp = cos(P->pitch);
	 sr = sin(P->roll);
	 cr = cos(P->roll);
	 sy = sin(P->yaw);
	 cy = cos(P->yaw);
	 rpy[0] = cy * cp;
	 rpy[1] = cy * sp * sr + sy * cr;
	 rpy[2] = sy * sr - cy * sp * cr;
	 rpy[3] = -sy * cp;
	  rpy[4] = cy * cr - sp * sr * sy;
	  rpy[5] = cy * sr + sy * sp * cr;
	  rpy[6] = sp;
	  rpy[7] = -cp * sr;
	  rpy[8] = cp * cr;
	  
	  /*     MULTIPLICATION OF MATRICES A AND RPY (BT=RPY*B) */	  
	  for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
		  P->bt[i+j*3] = rpy[j*3]+b[i]+rpy[j*3+1]*b[i+3]+rpy[j*3+2]*b[i+6];
			}
	  }
	} else {					/* No Pitch, Roll, or Yaw */
	  int i;
	  for (i=0; i<9; i++)
		P->bt[i]=b[i];
	}
}

static PJ * setup_row_col(PJ *P) {
  double incmax;
  double elvincr;
  double scnincr;
  double nsnom;

  /* Instrument specific values */
  /* Instru 0=Imager 1=Sounder */
  P->instr = pj_param(P->params, "iinstr").i;
    
  /* P->goes set from previous */
  /* First choose parameters based on instrument and GOES satellite,
	 then allow for parameters to be overridden */

  switch(P->goes) { 
  case(-1): /* This is a TEST case */
    if (P->instr==0) {
      P->nscyc=4;
      P->nsinc=3068;
      P->ewcyc=2;
      P->ewinc=3068;
    } else {
      P->nscyc=4;
      P->nsinc=1402;
      P->ewcyc=2;
      P->ewinc=1402;
    }
    break;
  case(-2): /* This is a TEST case */
    if (P->instr==0) {
      P->nscyc=4;
      P->nsinc=3068;
      P->ewcyc=2;
      P->ewinc=3068;
    } else {
      P->nscyc=4;
      P->nsinc=1402;
      P->ewcyc=2;
      P->ewinc=1402;
    }
    break;

  case(10):
    if (P->instr==0) {
      P->nscyc=4;
      P->nsinc=3299;
      P->ewcyc=2;
      P->ewinc=3025;
    } else {						/* SOUNDER is WRONG */
      P->nscyc=4;
      P->nsinc=1402;
      P->ewcyc=2;
      P->ewinc=1402;
    }
    break;
  case(11):
    if (P->instr==0) {
      P->nscyc=4;
      P->nsinc=2970;
      P->ewcyc=2;
      P->ewinc=3035;
    } else {						/* SOUNDER is WRONG */
      P->nscyc=4;
      P->nsinc=1402;
      P->ewcyc=2;
      P->ewinc=1402;
    }
    break;
  case(15):
    if (P->instr==0) {
      P->nscyc=4;
      P->nsinc=3012;
      P->ewcyc=2;
      P->ewinc=3127;
    } else {
      P->nscyc=4;
      P->nsinc=1402;
      P->ewcyc=2;
      P->ewinc=1402;
    }
    break;
  default:
    E_ERROR(-1);
  }
  
  if (P->instr==0) {
    incmax=6136;
    elvincr= 0.049087378124999997/incmax;
    scnincr= 0.098174756249999995/incmax;
    P->elvln = elvincr * 28 / 8;
    P->scnpx = scnincr;
    P->ewnom = incmax*2.5*scnincr;
    P->D=4.5;
  }
  else {
    incmax  = 2805;
    elvincr = 0.049087378124999997/ incmax;
    scnincr = 0.098174756249999995/ incmax;
    P->elvln   = elvincr * 64 / 4; 
    P->scnpx   = scnincr*8;
    P->ewnom = incmax*2.5*scnincr;
    P->D=2.5;
  }
  
  /* Use can reset these to real NADIR values */
  if (pj_param(P->params, "tnscyc").i) {
    P->nscyc = pj_param(P->params, "inscyc").i;
  }
  if (pj_param(P->params, "tnsinc").i) {
    P->nsinc = pj_param(P->params, "insinc").i;
  }
  
  if (pj_param(P->params, "tewcyc").i) {
    P->ewcyc = pj_param(P->params, "iewcyc").i;
  }
  if (pj_param(P->params, "tewinc").i) {
    P->ewinc = pj_param(P->params, "iewinc").i;
  }

  /* Now calculate the scanner maximum values */  
  if (P->instr==0) 
    {
      P->elvmax = (incmax * P->nscyc + P->nsinc) * elvincr;
      P->scnmax = (incmax * P->ewcyc + P->ewinc) * scnincr; 
	}
  else
	{
	  nsnom = incmax*4.5*elvincr;
	  P->elvmax = nsnom * 2.0 - (incmax * P->nscyc + P->nsinc)* elvincr;
	  P->scnmax = (incmax * P->ewcyc + P->ewinc) * scnincr;
	}
  return(P);
}

static PJ * setup(PJ *P) {
/* GOES Parameters MUST BE SET */

float oa1[336] = {0.0, 0.0, 0.0, 0.0, -1.747405052185, 84.06604003906,
            -0.34368492669e-1, -0.610297918319e-2, 0.0, 0.0, 0.0, 0x19890320,
            0x62934567, 0.0, 3.0e-4, -3.0e-4, -2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4,
            2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4,
            2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4,
            2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4,
            2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4,
            2.0e-4, 4.363e-3, 0.0,
            5.0e-4, 100.0, 2.0e-3, 0x0f, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0x04, 0x02, 0x02, 1.0e-5, 0.0,
            0.01, 0x02, 0x03, -1.0e-5, 0.0, 0.01, 0x03, 0x02, 1.0e-5, 0.0, 0.01, 0x03,
            0x03, -1.0e-5, 0.0, 0.01, 5.0e-4, 100.0, 2.0e-3, 0x0f, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0x04, 0x02, 0x02, 1.0e-5, 0.0, 0.01, 0x02, 0x03, -1.0e-5, 0.0, 0.01, 0x03, 0x02,
            1.0e-5, 0.0, 0.01, 0x03, 0x03, -1.0e-5, 0.0,
            0.01, 5.0e-4, 100.0, 2.0e-3, 0x0f, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0x04,
            0x02, 0x02, 1.0e-5, 0.0, 0.01, 0x02, 0x03, -1.0e-5, 0.0, 0.01,
            0x03, 0x02, 1.0e-5, 0.0, 0.01, 0x03, 0x03, -1.0e-5, 0.0,
            0.01, 5.0e-4, 100.0, 2.0e-3, 0x0f, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5,
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0x04, 0x02, 0x02, 1.0e-5, 0.0, 0.01, 0x02, 0x03,
            -1.0e-5, 0.0, 0.01, 0x03, 0x02, 1.0e-5, 0.0, 0.01, 0x03, 0x03, -1.0e-5, 0.0, 0.01, 5.0e-4, 100.0,
            2.0e-3, 0x0f, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 
            0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 
	    0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0.5e-5, 0x04, 0x02, 0x02, 1.0e-5, 0.0,
            0.01, 0x02, 0x03, -1.0e-5, 0.0, 0.01, 0x03, 0x02, 1.0e-5, 0.0, 0.01, 0x03, 0x03, -1.0e-5, 0.0, 0.01
        };

 int i;

 // Zero out oa at first
 for (i=0; i<336;i++)
   P->oa[i]=0;

P->goes = pj_param(P->params, "igoes").i;
switch(P->goes) {
 case(-1):	/* This is a TEST case */
    P->lam = -1.74740505;
    P->dr = 84066.04004;
    P->phi=-0.03436849;
    P->psi=-0.00610298;
    P->noaa=0;
    break;
 case(-2):			/* This is another test case */
   for (i=0; i<336;i++)
     P->oa[i]=oa1[i];
   P->lam=-1.747405052185;
   P->dr=84.06604003906;
   P->phi=-0.34368492669e-1;
   P->psi=-0.610297918319e-2;
   P->tu=20557829.57612;
   P->t=P->tu+20;
   P->noaa=1;
   P->imc=0;
   break;
 case(10):
    P->lam=-135*atan(1)/45;
    P->dr = 0;
    P->phi= 0;
    P->psi= 0;
    P->noaa=0;
    P->imc=1;
    break;
 case(11):
   P->lam=-135*atan(1)/45;
    P->dr = 0;
    P->phi= 0;
    P->psi= 0;
    P->noaa=0;
    P->imc=1;
   break;
 case(15):
   P->lam=-135*atan(1)/45;
    P->dr = 0;
    P->phi= 0;
    P->psi= 0;
    P->noaa=1;
    P->imc=1;
   break;
 default:
   E_ERROR(-1);
  }

  /* row-col (1= positive NS, -1 = neg NS, 0=scan-elev) default=1 */
  /* Setup Row/Col variables if needed */
 if (pj_param(P->params,"trow_col").i)
   P->row_col = pj_param(P->params,"irow_col").i;
 else
   P->row_col = 1;
  
 /* Allow new Sub satellite longitude */
 if (pj_param(P->params,"tlam").i) 
   P->lam = pj_param(P->params, "dlam").f;
 P->oa[4]=P->lam;

 /* Normalized earth radius */
 /* 42164365 is goes radius in meters */
 if (pj_param(P->params,"tdr").i) 
   P->dr = pj_param(P->params, "ddr").f;
 P->oa[5]=P->dr;

 if (pj_param(P->params,"tphi").i)  
   P->phi = pj_param(P->params, "dphi").f;
 P->oa[6]=P->phi;

 if (pj_param(P->params,"tpsi").i) 
   P->psi = pj_param(P->params, "dpsi").f;
 P->oa[7]=P->psi;

 if (pj_param(P->params,"troll").i) 
   P->roll = pj_param(P->params, "droll").f;
 else
   P->roll=0;
 P->oa[8]=P->roll;

 if (pj_param(P->params,"tpitch").i) 
   P->pitch = pj_param(P->params, "dpitch").f;
 else
   P->pitch=0;
 P->oa[9]=P->pitch;

 if (pj_param(P->params,"tyaw").i) 
   P->yaw = pj_param(P->params, "dyaw").f;
 else
   P->yaw=0;
 P->oa[10]=P->yaw;

 if (pj_param(P->params,"tt").i) 
   P->t = pj_param(P->params, "dt").f;
 else
   P->t=0;

 if (pj_param(P->params,"ttu").i) 
   P->tu = pj_param(P->params, "dtu").f;
 else
   P->tu=0;

 /* When imc status on then IMC=0 ???, but that's how it is */
  if (pj_param(P->params, "timc").i) {
    P->imc = (pj_param(P->params, "bimc").i==0)?1:0;
    //    P->imc = (pj_param(P->params, "bimc").i==0)?0:1;
  } else
    P->imc=1;			/* IMC is off = GOOD */

  if (pj_param(P->params, "tflip").i)
    P->flip = (pj_param(P->params, "bflip").i==0)?1:-1;
  else
    P->flip=1;

  if (pj_param(P->params, "tnoaa").i) 
    P->noaa = pj_param(P->params, "bnoaa").i;

 /* New all IMC in one go. */
 if (pj_param(P->params,"toa").i) {
   char *input = pj_param(P->params,"soa").s;
   char *delim = ","; // input separated by spaces
   char *token = NULL;
   int tnum = 0;
   for (token = strtok(input, delim); token != NULL; token = strtok(NULL, delim))
     {
       char *unconverted;
       double val = strtod(token, &unconverted);
       if (!isspace(*unconverted) && *unconverted != 0) 
	 {
	   E_ERROR(*unconverted)
	     } else {
	 P->oa[tnum]=val;
       }
       tnum++;
     }
 }

 P->r = (P->dr + 42164365.) / P->a;

 if (P->noaa) {
   // ERROR - fix these syncs?  for GOES15
   setup_row_col(P);
   // I put in the same values for now as they are set by instr
   spscons(P->nscyc,P->nsinc,P->ewcyc,P->ewinc,
	   P->nscyc,P->nsinc,P->ewcyc,P->ewinc,
	   &P->spsarea, &P->cilpevsc);
   gmodel(P->t, P->tu, P->oa, P->imc, &P->rlat, &P->rlon, &P->spsarea);
   P->inv = e_post11_inverse;
   P->fwd = e_post11_forward;
 } else {
   setup_row_col(P);  
   setup_pre11(P);
   P->inv = e_pre11_inverse;
   P->fwd = e_pre11_forward;
 }

  return (P);
}

ENTRY0(goes)
ENDENTRY(setup(P))
