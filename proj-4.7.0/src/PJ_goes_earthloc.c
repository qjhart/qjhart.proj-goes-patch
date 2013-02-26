/************************************************************************

     Integral Systems, Inc.

 ***********************************************************************

     Project   : Operations Ground Equipment for GOES
     System    : Sensors Processing System
     Functions : gmodel, gatt, inst2er, spscons, time50, 
                 gpoint, lnpxevsc, scpx, ev2ci, sc2ci, pxci, 
                 lpoint, limb, clipper, sndeloc, evci
     Source    : earthloc.c
     Programmer: Igor Levine

     ver.    data    by   comment
     ----  --------  ---  ---------------------------------------------
     1     03/05/93  IL   Converted from FORTRAN to C.
     2     03/13/93  IL   Created function clipper()
     3     04/28/93  IL   Added trap for the case slat = syaw = 0
                          in function gmodel()
     4     11/29/93  IL   Added factor of 2 for two last terms in the
                          expression for lam in gmode().
     5     11/30/93  IL   Implemented new formulae for misalignment
                          corrections of scan angles.  Affected
                          functions: gpoint, lpoint, limb, clipper.
     6     04/29/94  IL   Deleted roll & pitch misalignments from the
                          expressions for roll & pitch in gmodel().
     7     06/01/94  IL   Added arguments to spscons()
     V1.6  07/20/94  ncf  mask out flywheel bit from bcd day (in time50)
     8     04/14/95  IL   Added origin offset correction for Sounder
                          (affected functions: evci, ev2ci, sndeloc).
     9     04/18/95  IL   Added functions evsc2lp, px2ci, utime50,
                          tupack50.
    10     05/01/95  IL   Added function slant_range().
    11     05/09/95  IL   In gatt() the false library function pow(x,y)
                          was replaced with macro function IPOW().
    12     05/12/95  IL   Added ev2ln(), touched lnpxevsc().
    13     05/17/95  IL   Deleted obsolete variables from spscons().
    14     12/02/97  IL   Implemented new formulae for invert operations;
                          implemented new misalignment corrections for
                          the Sounder. Affected functions: gpoint, lpoint,
                          limb, clipper, sndeloc, evci, pxci, saci, scpx,
                          ev2ci, sc2ci, px2ci, evsc2lp, evln, ev2ln,
                          lnpxevsc.
    15     22/01/98  IL   Changed direction of image rotation in
                          sndeloc(). Implemented elevation and scan
                          angle corrections due to nadir offsets;
                          affected functions: gpoint, lpoint, clipper.
 ***********************************************************************/

#include "PJ_goes_earthloc.p"
#include <math.h>
/* #include "earthloc.p"*/
#ifndef  IPOW
#define  IPOW(x, n) (((n%2) ? sign((x)) : 1))*exp(n*log(ABS((x))))
#endif

//#define TEST 

/***********************************************************************
   The function gmodel() computes the position of the satellite and 
   the attitude of the imager or sounder.  the calculations are based
   on the oats orbit and attitude model represented by the O&A
   parameter set in GVAR Block 0.
     
     Inputs:
          time, epoch time, O&A parameter set, IMC status.

     Outputs:
          the spacecraft position vector in earth fixed coordinates; 
          the earth fixed to spacecraft coordinates transformation
          matrix; the geometric roll, pitch, yaw angles and the roll,
          pitch misalignments for either the imager or the sounder;  
          the earth fixed to instrument frame transformation matrix; 
          geographic latitude and longitude at subsatellite point.

   Description:
   gmodel() accepts an input double precision time in minutes from
   1950, Jan.1.0 and an input set of O&A parameters and computes
   position of the satellite, the attitude angles and attitude
   misalignments and the instrument to earth fixed coordinates
   transformation matrix.

   Arguments:

      t    = input time in minutes from 1950, January 1.0 UTC
      tu   = input epoch time in minutes from 1950 January 1.0 UTC
      rec[] = input O&A parameter set
      imc   = input IMC status: 0 - on, 1 - off 
      rlat  = subsatellite geodetic latitude in radians
      rlon  = subsatellite longitude in radians    
      spsarea = pointer to DB structure

 ************************************************************************* */

void gmodel(double t, double tu, float rec[], int imc,
        float *rlat, float *rlon, SPSAREA *spsarea) {
    double a1, a2, gatt();
    double asc, c2w, ca, ci, cu, cw, cw1, cw3, phi, psi,
            s2w, sa, si, slat, su, sw, sw1, sw3, syaw, te, u, w, wa,
            dlat, dr, dyaw, lam, r, ts;
    void inst2er();


    /*  Assign reference values to the subsatellite longitude and
        latitude, the radial distance and the orbit yaw.   */
#ifdef TEST
    tu = t - 525.0;
#endif
    lam = (double) rec[4];
    dr = (double) rec[5];
    phi = (double) rec[6];
    psi = (double) rec[7];

    /*     Assign reference values to the attitudes and misalignments  */
    spsarea->roll = (double) rec[8];
    spsarea->pitch = (double) rec[9];
    spsarea->yaw = (double) rec[10];
    spsarea->rma = 0.;
    spsarea->pma = 0.;

    /*     If imc is off, compute changes in the satellite orbit  */
    if (imc != 0) {
        /* Set to zero reference radial distance, latitude and orbit yaw */
        dr = 0.;
        phi = 0.;
        psi = 0.;

        /*  Determine time since epoch (in minutes) */
        ts = t - tu;

        /*  Computes orbit angle and the related trigonometric functions.
            earth rotational rate=.7292115e-4 (rad/s)  */
        w = 0.7292115e-4 * ts * 60.;
        sw = sin(w);
        cw = cos(w);
        sw1 = sin(.927 * w);
        cw1 = cos(.927 * w);
        s2w = sin(2. * w);
        c2w = cos(2. * w);
        sw3 = sin(1.9268 * w);
        cw3 = cos(1.9268 * w);

	//        printf(" x= %f, sw= %f, cw= %f, sw1= %f, cw1= %f, s2w= %f, c2w= %f, sw3= %f, cw3= %f \n", w, sw, cw, sw1, cw1, s2w, c2w, sw3, cw3);
        /*  Computes change in the imc longitude from the reference */
        lam = lam + (double) rec[17] + ((double) rec[18] +
                (double) rec[19] * w) * w + ((double) rec[20] * sw +
                (double) rec[21] * cw + (double) rec[22] * s2w +
                (double) rec[23] * c2w + (double) rec[24] * sw3 +
                (double) rec[25] * cw3 + w * ((double) rec[28] * sw +
                (double) rec[29] * cw) + (double) rec[26] * sw1 +
                (double) rec[27] * cw1)*2.;

        /*  Computes change in radial distance from the reference (km) */
        dr = dr + (double) rec[30] + (double) rec[31] * cw +
                (double) rec[32] * sw + (double) rec[33] * c2w +
                (double) rec[34] * s2w + (double) rec[35] * cw3 +
                (double) rec[36] * sw3 + (double) rec[37] * cw1 +
                (double) rec[38] * sw1 + w * ((double) rec[39] * cw +
                (double) rec[40] * sw);

        /*  Computes the sinus of the change in the geocentric latitude */
        dlat = (double) rec[41] + (double) rec[42] * cw + (double) rec[43] * sw +
                (double) rec[44] * c2w + (double) rec[45] * s2w +
                w * ((double) rec[46] * cw + (double) rec[47] * sw) +
                (double) rec[48] * cw1 + (double) rec[49] * sw1;

        /*  Computes geocentric latitude by using an expansion for arcsin */
        phi = phi + dlat * (1. + dlat * dlat / 6.);

        /*  Computes sinus of the change in the orbit yaw  */
        dyaw = (double) rec[50] + (double) rec[51] * sw + (double) rec[52] * cw +
                (double) rec[53] * s2w + (double) rec[54] * c2w +
                w * ((double) rec[55] * sw + (double) rec[56] * cw) +
                (double) rec[57] * sw1 + (double) rec[58] * cw1;

        /*  Computes the orbit yaw by using an expansion for arcsinus. */
        psi = psi + dyaw * (1. + dyaw * dyaw / 6.);

	//        printf(" lam= %f, dr= %f, dlat=%f, phi= %f, dyaw= %f, psi= %f \n", lam, dr, dlat, phi, dyaw, psi);
        /*  Calculation of changes in the satellite orbit ends here */
    }

    /* Conversion of the imc longitude and orbit yaw to the subsatellite
       longitude and the orbit inclination (ref: GOES-PCC-TM-2473, inputs
       required for earth location and gridding by SPS,  June 6, 1988) */

    slat = sin(phi);
    syaw = sin(psi);
    si = SQ(slat) + SQ(syaw);
    ci = sqrt(1. - si);
    si = sqrt(si);

    /* Check if slat = syaw = 0  */
    if (slat == 0. && syaw == 0.)
        u = 0.;
    else
        u = atan2(slat, syaw);
    su = sin(u);
    cu = cos(u);

    /* Computes longitude of the ascending node */
    asc = lam - u;
    sa = sin(asc);
    ca = cos(asc);

    /* Computes the subsatellite geographic latitude */
    *rlat = (float) (atan(spsarea->aebe2 * tan(phi)));

    /*  Computes the subsatellite longitude */
    *rlon = (float) (asc + atan2(ci*su, cu));

    /*  Computes the spacecraft to earth fixed coordinates transformation
        matrix:  (vector in ECEF coord.) = b * (vector in s/c coord.) */
    //    printf("ca= %f, sa= %f, ci=%f, si= %f, cu= %f, su= %f, slat= %f\n", ca, sa, ci, si, cu, su, slat);
    spsarea->b[2][0] = sa*si;
    spsarea->b[2][1] = -ca*si;
    spsarea->b[2][2] = ci;
    spsarea->b[0][0] = -ca * cu + sa * su*ci;
    spsarea->b[0][1] = -sa * cu - ca * su*ci;
    spsarea->b[0][2] = -slat;
    spsarea->b[1][0] = ca * su + sa * cu*ci;
    spsarea->b[1][1] = sa * su - ca * cu*ci;
    spsarea->b[1][2] = -cu*si;

    /*  Computes the normalized spacecraft position vector in earth fixed
        coordinates - xs. */
    r = (spsarea->nomorb + dr) / spsarea->ae;
    spsarea->xs[0] = -spsarea->b[0][0] * r;
    spsarea->xs[1] = -spsarea->b[0][1] * r;
    spsarea->xs[2] = -spsarea->b[0][2] * r;

    /*  Precomputes q3 (used in lpoint) */
    spsarea->q3 = SQ(spsarea->xs[0]) + SQ(spsarea->xs[1]) +
            spsarea->aebe2 * SQ(spsarea->xs[2]) - 1.0;

    /*  Computes the attitudes and misalignments if IMC is off */

    if (imc != 0) {
        /* Computes the solar orbit angle */
        wa = (double) rec[59] * ts;

        /*  Computes the difference between current time, ts, and the
            exponential time, rec[60]. Note that both times are since
            epoch. */
        te = ts - (double) rec[60];

        /*  Computes roll */
        a1 = gatt(61, rec, wa, te);

        /*  Computes pitch */
        a2 = gatt(116, rec, wa, te);

        /*  Computes yaw */
        spsarea->yaw += gatt(171, rec, wa, te);

        /*  Computes roll misalignment */
        spsarea->rma = gatt(226, rec, wa, te);

        /*  Computes pitch misalignment */
        spsarea->pma = gatt(281, rec, wa, te);

	//        printf("a1= %f, a2=%f, yaw=%f, rma=%f, pma=%f", a1, a2, spsarea->yaw, spsarea->rma, spsarea->pma);
        /*  Computes the true geometrical roll and pitch */
        spsarea->roll = spsarea->roll + a1; /* corrected 4/29/94 */
        spsarea->pitch = spsarea->pitch + a2;

        /*  Apply the spacecraft compensation */
        spsarea->roll = spsarea->roll + (double) rec[14];
        spsarea->pitch = spsarea->pitch + (double) rec[15];
        spsarea->yaw = spsarea->yaw + (double) rec[16];
    }

    /* Computes the instrument to earth fixed coordinates
       transformation matrix - bt */
    inst2er(spsarea->roll, spsarea->pitch, spsarea->yaw,
            spsarea->b, spsarea->bt);
    /**
    int i, j;
    for(i=0; i<3; i++) {
        for(j=0; j<3; j++) {
            printf("i=%d, j= %d, spsarea.b= %f, spsarea.bt=%f\n", i, j, spsarea->b[i][j], spsarea->bt[i][j]);
        }
    } */
    return;
}

/************************************************************************

      The function gatt() computes an attitude/misalignment angle 
      from a given subset of the O&A parameters in GVAR Block 0.
      Argument k0 indicates the first word of the subset.

      Arguments:

         k0  = starting position of a parameter subset in the O&A set 
         rec = input O&A parameter set 
         wa  = input solar orbit angle in radians 
         te  = input exponential time delay from epoch (minutes)  

 *********************************************************************** */

double gatt(int k0, float rec[], double wa, double te) {
    int k, l, ll;
    double datt, wrk1, wrk2;
    QUATER I, J, M;

    k = k0;
    datt = (double) rec[k + 2];
    //printf("1. datt= %f\n", datt);
    /*  Computes the exponential term  */

    if (te >= 0. && (double) rec[k] > 0.0)
        datt += (double) rec[k] * exp(-te / (double) rec[k + 1]);
    //printf("2. datt= %f\n", datt);
    /*  Extracts the number of sinusoids */
    I.r4 = rec[k + 3];

    /*  Calculation of sinusoids */
    for (l = 1; l <= (int) I.r4; l++) {
        wrk1 = wa * (double) l + (double) rec[k + l * 2 + 3];
        datt += (double) rec[k + l * 2 + 2] * cos(wrk1);
    }
    //printf("3. datt= %f I.r4= %f, rec[64]=%f\n", datt, I.r4, rec[64]);
    /*  Pointer to the number of monomial sinusoids */
    k = k + 34;

    /*  Extracts number of monomial sinusoids */
    I.r4 = rec[k];

    /* Computes monomial sinusoids */
    for (l = 1; l <= (int) I.r4; l++) {
        ll = k + 5 * l;

        /* Order of sinusoid (c2) */
        J.r4 = rec[ll - 4];

        /* Order of monomial sinusoid (c5) */
        M.r4 = rec[ll - 3];
        wrk1 = wa * (double) J.r4 + (double) rec[ll - 1];
        wrk2 = pow(wa - (double) rec[ll], (double) M.r4);
#if 0   
        wrk2 = IPOW(wa - (double) rec[ll], M.i4); /* for REAL/IX */
#endif
        //printf("wrk1= %f, wrk2=%f \n", wrk1, wrk2);
        datt += (double) rec[ll - 2] * wrk2 * cos(wrk1);
    }
    //printf(" datt= %f\n", datt);
    return ( datt);
}

/***********************************************************************

     The function gpoint() computes elevation and scan angles 
     which aim the instrument at a point on the earth with given 
     geographic coordinates.

     Arguments:
        instr = instrument code (1-imager, 2-sounder)
        flip = instrument orientation flag (1-nominal, (-1)-inverted)
        rlat    =  geographic latitude in radians (input) 
        rlon    =  geographic longitude in radians (input)
        alf     =  elevation angle in radians (output) 
        gam     =  scan angle in radians (output) 
        ierr    =  output status; 0 - successful completion, 
                   1 - point with given lat/lon is invisible 
        spsarea = pointer to DB structure
        cilpevsc = pointer to DB structure
 ************************************************************************/

void gpoint(int instr, int flip, float rlat, float rlon, float *alf,
        float *gam, int *ierr, SPSAREA *spsarea,
        CILPEVSC *cilpevsc) {
    double ft[3], sing, slat, u[3], w1, w2, f[3], sa, cz, z, FF,
            Doff, alpha, alpha1;

    FF = (double) flip;
    Doff = FF * (cilpevsc->zetmax[instr - 1] -
            2.5 * (float) cilpevsc->ixmax[instr - 1]
            * cilpevsc->zetincr[instr - 1]);
    if (instr == 2) FF = -FF;

    /*  Computes sinus of geographic (geodetic) latitude */
    sing = sin((double) rlat);
    w1 = spsarea->aebe4 * sing*sing;

    /*  Sinus of the geocentric latitude */
    slat = ((0.375 * w1 - 0.5) * w1 + 1.) * sing / spsarea->aebe2;

    /*  Computes local earth radius at specified point */
    w2 = slat*slat;
    w1 = spsarea->aebe3*w2;
    w1 = (0.375 * w1 - 0.5) * w1 + 1.;

    /*  Computes Cartesian coordinates of the point */
    u[2] = slat*w1;
    w2 = w1 * sqrt(1. - w2);
    u[0] = w2 * cos((double) rlon);
    u[1] = w2 * sin((double) rlon);

    /*  Pointing vector from satellite to the earth point */
    f[0] = u[0] - spsarea->xs[0];
    f[1] = u[1] - spsarea->xs[1];
    f[2] = u[2] - spsarea->xs[2];
    w2 = u[0] * f[0] + u[1] * f[1] + u[2] * f[2] * spsarea->aebe2;

    /*  Verifies visibility of the point */
    if (w2 > 0.) {
        /* invisible point on the earth */
        *ierr = 1;
        *alf = 99999.;
        *gam = 99999.;
        return;
    }

    /* Converts pointing vector to instrument coordinates */
    ft[0] = spsarea->bt[0][0] * f[0] + spsarea->bt[0][1] * f[1] +
            spsarea->bt[0][2] * f[2];
    ft[1] = spsarea->bt[1][0] * f[0] + spsarea->bt[1][1] * f[1] +
            spsarea->bt[1][2] * f[2];
    ft[2] = spsarea->bt[2][0] * f[0] + spsarea->bt[2][1] * f[1] +
            spsarea->bt[2][2] * f[2];

    /*  Converts pointing vector to scan and elevation angles and
        correct the angles for the origin offset */

    alpha = atan(ft[2] / ft[0]); /* elevation */
    w2 = -ft[1] / sqrt(SQ(ft[0]) + SQ(ft[2])); /* tan(scan angle) */
    z = atan(w2); /* scan angle */
    //printf("noaa: %f %f\n",z,alpha);
    alpha1 = alpha + alpha * z*Doff; /* corrected elevation */
    z -= 0.5 * alpha * alpha*Doff; /* corrected scan angle */
    //printf("noaa:%f %f doff=%f\n",z,alpha1,Doff);
    sa = sin(alpha1);
    cz = cos(z);

    /*     Corrects for the roll and pitch misalignments */

    *gam = (float) (z - FF * spsarea->rma * sa);
    *alf = (float) (alpha1 + spsarea->rma * (1. - cos(alpha1) / cz)
            + spsarea->pma * sa * (tan(z) + FF / cz));
    *ierr = 0;
    //printf("noaa: %f %f\n",*gam,*alf);
    return;
}

/************************************************************************

     The function lpoint() computes geographic coordinates of the earth
     point defined by the elevation and scan angles of the instrument.

     Arguments:
        instr = instrument code (1-imager, 2-sounder)
        flip = instrument orientation flag (1-nominal, (-1)-inverted)
        alpha0    = elevation angle (rad) 
        zeta0     = scan angle (rad)
        rlat      = latitude in radians (output)
        rlon      = longitude in radians (output) 
        ierr      = output status; 0 - point on the earth,
                    1 - instrument points off earth
        spsarea = pointer to DB structure
        cilpevsc = pointer to DB structure

 ************************************************************************/

void lpoint(int instr, int flip, float alpha0, float zeta0,
        float *rlat, float *rlon, int *ierr,
        SPSAREA *spsarea, CILPEVSC *cilpevsc) {
    double ca, sa, cz, sz, da, dz, g[3], d, g1[3], h, Doff,
            q1, q2, u[3], d1, FF, alpha, zeta;

    *ierr = 1;
    FF = (double) flip;

    /*  Computes the second order corrections of the elevation
        and scan angles due to the origin offset */
    Doff = FF * (cilpevsc->zetmax[instr - 1] -
            2.5 * (float) cilpevsc->ixmax[instr - 1]
            * cilpevsc->zetincr[instr - 1]);
    alpha = alpha0 - alpha0 * zeta0*Doff;
    zeta = zeta0 + 0.5 * alpha0 * alpha0*Doff;

    //printf("alpha=%f, zeta= %f\n", alpha, zeta);
    /*  Computes the scan and elevation angles corrections
        due to the roll and pitch misalignments */
    if (instr == 2) FF = -FF;
    ca = cos(alpha);
    sa = sin(alpha);
    cz = cos(zeta);
    sz = sin(zeta);
    da = -(spsarea->pma * sa * (FF + sz)
            + spsarea->rma * (cz - ca)) / cz; /* elev. angle correction */
    dz = FF * spsarea->rma*sa; /* scan angle correction */

    /*  Computes pointing vector in instrument coordinates */
    cz = cos(zeta + dz);
    g[0] = cos(alpha + da) * cz;
    g[1] = -sin(zeta + dz);
    g[2] = sin(alpha + da) * cz;

    /*  Transforms the pointing vector to earth fixed coordinates */
    g1[0] = spsarea->bt[0][0] * g[0] + spsarea->bt[1][0] * g[1]
            + spsarea->bt[2][0] * g[2];
    g1[1] = spsarea->bt[0][1] * g[0] + spsarea->bt[1][1] * g[1]
            + spsarea->bt[2][1] * g[2];
    g1[2] = spsarea->bt[0][2] * g[0] + spsarea->bt[1][2] * g[1]
            + spsarea->bt[2][2] * g[2];

    /*  Computes coefficients and solves a quadratic equation to
        find the intersect of the pointing vector with the earth
        surface  */
    q1 = SQ(g1[0]) + SQ(g1[1]) + spsarea->aebe2 * SQ(g1[2]);
    q2 = spsarea->xs[0] * g1[0] + spsarea->xs[1] * g1[1]
            + spsarea->aebe2 * spsarea->xs[2] * g1[2];
    d = q2 * q2 - q1 * spsarea->q3;
    if (fabs(d) < 1.e-9)
        d = 0.;

    /*  If the disciminante of the equation, d, is negative, the
        instrument points off the earth */
    if (d < 0) {
        *rlat = 999999.;
        *rlon = 999999.;
        return;
    }
    d = sqrt(d);

    /*  Slant distance from the satellite to the earth point */
    h = -(q2 + d) / q1;

    /*  Cartesian coordinates of the earth point */
    u[0] = spsarea->xs[0] + h * g1[0];
    u[1] = spsarea->xs[1] + h * g1[1];
    u[2] = spsarea->xs[2] + h * g1[2];

    /*  Sinus of geocentric latitude */
    d1 = u[2] / sqrt(SQ(u[0]) + SQ(u[1]) + SQ(u[2]));

    /* int i; */
    /* for(i=0; i<3; i++) { */
    /*     printf(" i= %d, g=%f, g1=%f, u= %f, spsarea.xs=%f\n", i, g[i], g1[i], u[i], spsarea->xs[i]); */
    /* } */
    /*  Geographic (geodetic) coordinates of the point */
    *rlat = (float) (atan(spsarea->aebe2 * d1 / sqrt(1. - d1 * d1)));
    *rlon = (float) (atan2(u[1], u[0]));
    *ierr = 0;
    return;
}

/************************************************************************

     inst2er() accepts the single precision roll, pitch and yaw angles
     of an instrument and returns the double precision instrument to
     earth coordinates transformation matrix.

       Arguments:                             
          r  = roll angle in radians
          p  = pitch angle in radians
          y  = yaw angle in radians 
          a  = spacecraft to ECEF coordinates transformation matrix
          at = instrument to ECEF coordinates transformation matrix

 *********************************************************************** */

void inst2er(double r, double p, double y, double a[][3], double at[][3]) {
    int i, j;
    double rpy[3][3];

    /* Instrument to body coordinates transformation matrix is
       computed by using a small angle approximation of trigonometric
       functions of the roll, pitch and yaw. */
    rpy[0][0] = 1. - 0.5 * (p * p + r * r);
    rpy[1][0] = p - r*y;
    rpy[2][0] = -r - p*y;
    rpy[0][1] = -p;
    rpy[1][1] = 1. - 0.5 * (y * y + p * p);
    rpy[2][1] = -y;
    rpy[0][2] = r;
    rpy[1][2] = y + r*p;
    rpy[2][2] = 1. - 0.5 * (y * y + r * r);

    /*  Multiplication of matrices A and RPY */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            at[j][i] = a[0][i] * rpy[j][0] + a[1][i] * rpy[j][1]
                + a[2][i] * rpy[j][2];
    return;
}

/*********************************************************************  
   evci() converts instrument cycles/increments and
   the detector number to the elevation angle of the
   detector's optical axis.
  
   Arguments:
      instr = instrument code (1-imager, 2-sounder)
      flip = instrument orientation flag (1-nominal, (-1)-inverted)
      cy    = instrument cycles   
      incy  = instrument increments   
      idet  = index of a visible detector (1-8  for imager,    
              1-4 for sounder, 0 corresponds to the optical 
              axis of the instrument    

 *********************************************************************** */

float evci(int instr, int flip, int cy, int incy, int idet,
        CILPEVSC *cilpevsc) {
    float d;
    long iw;

    if (idet == 0)
        d = 0.;
    else {
        if (instr == 1)
            d = ((float) idet - 4.5) * cilpevsc->alfln[instr - 1];
        else
            d = (2.5 - (float) idet) * cilpevsc->alfln[instr - 1];
    }
    iw = (long) cy * (long) cilpevsc->iymax[instr - 1] + (long) incy;
    if (instr == 1)
        d = (float) flip * (cilpevsc->alfmax[0] - (float) iw *
            cilpevsc->alfincr[0]) - d;
    else
        d = (float) flip * (cilpevsc->alfmax[1] + (float) (iw -
            9 * (long) cilpevsc->iymax[1]) * cilpevsc->alfincr[1]) + d;
    return (d);
}

/***********************************************************************

   lnpxevsc() converts the elevation and scan angles
   to the nearest absolute scan, visible detector line and
   pixel numbers and the fractional detector line and pixel
   offsets.

   Arguments:
      instr = instrument code (1-imager, 2-sounder)
      flip = instrument orientation flag (1-nominal, (-1)-inverted)
      alpha = elevation angle (rad) 
      zeta  = scan angle (rad) 
      iscan = absolute scan number 
      idet  = detector line number in the scan
      line  = absolute line number 
      ipix  = pixel number 
      offln = vertical offset relative to idet 
      offpx = horizontal offset relative to ipix 

 ************************************************************************/
void lnpxevsc(int instr, int flip, float alpha, float zeta,
        int *iscan, int *idet, int *line, int *ipix,
        float *offln, float *offpx, CILPEVSC *cilpevsc) {
    int i;
    float rl, rp;
    if (instr == 1)
        rl = 4.5;
    else
        rl = 2.5;

    if (flip == 1) {
        rl = (cilpevsc->alfmax[instr - 1] - alpha) /
                cilpevsc->alfln[instr - 1] + rl;
        rp = (cilpevsc->zetmax[instr - 1] + zeta) /
                cilpevsc->zetpx[instr - 1] + 1.;
    } else {
        rl = 9. * cilpevsc->iymax[instr - 1] / cilpevsc->incrln[instr - 1]
                - (cilpevsc->alfmax[instr - 1] + alpha) /
                cilpevsc->alfln[instr - 1] + rl;
        rp = 5. * cilpevsc->ixmax[instr - 1] / cilpevsc->incrpx[instr - 1]
                - (cilpevsc->zetmax[instr - 1] - zeta) /
                cilpevsc->zetpx[instr - 1] + 1.;
    }
    *line = rl + 0.5 * sign(rl); /* integer output line # */
    *offln = rl - (float) (*line); /* vertical offset (det. line) */
    if (instr == 1) {
        i = (*line - 1) / 8;
        *idet = *line - 8 * i; /* detector line # within the scan */
    } else {
        i = (*line - 1) / 4;
        *idet = *line - 4 * i;
    }
    *iscan = i + 1;
    *ipix = rp + 0.5 * sign(rp); /* pixel number */
    *offpx = rp - (float) (*ipix); /* horizontal offset (pixel) */
    return;
}

/***********************************************************************

   pxci() converts instrument cycles/increments to
   the fractional pixel number.

   Arguments:
      instr = instrument code (1-imager, 2-sounder)
      flip = instrument orientation flag (1-nominal, (-1)-inverted)
      cx    = instrument cycles 
      incx  = instrument increments

 ************************************************************************/
float pxci(int instr, int flip, int cx, int incx, CILPEVSC *cilpevsc) {
    if (flip == 1)
        return ( 1. + (float) (cx * cilpevsc->ixmax[instr - 1]
            + incx) / cilpevsc->incrpx[instr - 1]);
    else
        return ( 1. + (float) ((5 - cx) * cilpevsc->ixmax[instr - 1]
            - incx) / cilpevsc->incrpx[instr - 1]);
}

/***********************************************************************

   saci() converts instrument cycles/increments
   to the scan angle.

   Arguments:
      instr = instrument code (1-imager, 2-sounder)
      flip = instrument orientation flag (1-nominal, (-1)-inverted)
      cx    = instrument cycles
      incx  = instrument increments 

 ************************************************************************/
float saci(int instr, int flip, int cx, float incx, CILPEVSC *cilpevsc) {
    return ( (float) flip * ((cx * cilpevsc->ixmax[instr - 1] + incx)
            * cilpevsc->zetincr[instr - 1] - cilpevsc->zetmax[instr - 1]));
}

/***********************************************************************

   scpx() converts fractional pixel number to scan angle.

   Arguments:
      instr = instrument code (1-imager, 2-sounder)
      flip = instrument orientation flag (1-nominal, (-1)-inverted)
      pix   = fractional pixel number

 ************************************************************************/
float scpx(int instr, int flip, float pix, CILPEVSC *cilpevsc) {
    if (flip == 1)
        return ( (pix - 1.) * cilpevsc->zetpx[instr - 1] -
            cilpevsc->zetmax[instr - 1]);
    else
        return ( (pix - 1.) * cilpevsc->zetpx[instr - 1] -
            5. * (float) cilpevsc->ixmax[instr - 1] * cilpevsc->zetincr[instr - 1]
            + cilpevsc->zetmax[instr - 1]);
}

/***********************************************************************

     ev2ci() accepts a single precision elevation angle of a 
     specific detector and computes the related pointing of the
     instrument expressed in cycles/increments coordinates.
     
     Arguments:
          instr = instrument code (1-imager, 2-sounder)
          flip = instrument orientation flag (1-nominal, (-1)-inverted)
          ev    = elevation angle in radians
          idet  = index of a visible detector (1-8  for imager, 
                  1-4 for sounder, 0 corresponds to the optical 
                  axis of the instrument) 
          cy    = instrument cycles
          incy  = instrument increments 

 ************************************************************************/

void ev2ci(int instr, int flip, float ev, int idet, int *cy, int *incy,
        CILPEVSC *cilpevsc) {
    int i;
    float r;

    r = cilpevsc->alfmax[instr - 1] - (float) flip*ev;
    if (idet != 0) {
        if (instr == 1)
            r -= ((float) idet - 4.5) * cilpevsc->alfln[0] * flip;
        else
            r -= ((float) idet - 2.5) * cilpevsc->alfln[1] * flip;
    }
    if (instr == 2)
        r = 9. * (float) cilpevsc->iymax[1] * cilpevsc->alfincr[1] - r;
    i = (int) (r / cilpevsc->alfincr[instr - 1] + 0.5);
    *cy = i / cilpevsc->iymax[instr - 1];
    *incy = i - *cy * cilpevsc->iymax[instr - 1];
    return;
}

/***********************************************************************

     sc2ci accepts single precision scan angle of the instrument
     and expresses it in integer cycles and increments
     
     Arguments:
          instr = instrument code (1-imager, 2-sounder)
          flip = instrument orientation flag (1-nominal, (-1)-inverted)
          sc    = scan angle in radians
          cx    = instrument cycles
          incx  = instrument increments 

 ************************************************************************/
void sc2ci(int instr, int flip, float sc, int *cx, int *incx,
        CILPEVSC *cilpevsc) {
    int i;
    i = (int) ((cilpevsc->zetmax[instr - 1] + (float) flip * sc) /
            cilpevsc->zetincr[instr - 1] + 0.5);
    *cx = i / cilpevsc->ixmax[instr - 1];
    *incx = i - *cx * cilpevsc->ixmax[instr - 1];
    return;
}

/**********************************************************************

   spscons() set constants in structures CILPEVSC and SPSAREA

 ***********************************************************************/

void spscons(int NScyc1, int NSinc1, int EWcyc1, int EWinc1,
        int NScyc2, int NSinc2, int EWcyc2, int EWinc2,
        SPSAREA *spsarea, CILPEVSC *cilpevsc)
/*
    NScyc1, NSinc1 - the Imager frame origin in N-S cycles/increments
    EWcyc1, EWinc1 - the Imager frame origin in E-W cycles/increments
    NScyc2, NSinc2 - the Sounder frame origin in N-S cycles/increments
    EWcyc2, EWinc2 - the Sounder frame origin in E-W cycles/increments
 */ {
    double fe, w, dy, dx;

    /*  Generates numerical constants */

    spsarea->pi = 3.141592653589793;
    spsarea->pi2 = spsarea->pi + spsarea->pi;
    spsarea->rad = spsarea->pi / 180.0;

    /*  Nominal orbit radius in km  */
    spsarea->nomorb = 42164.365;

    /*  Earth equatorial radius in km  */
    spsarea->ae = 6378.137;

    /*  Earth polar radius in km  */
    spsarea->be = 6356.7533;

    /*  Earth flattening  */

    /*  Note: It is assumed that the earth shape is defined via ae and be.
              alternatively it may be defined via earth equatorial radius,
              ae, and the earth flattening  fe=1.-(be/ae)
              in the last case, earth polar radius, be, has to be found as

               be=ae*(1.-fe)
     */

    /* Constants related to ae and be  */
    spsarea->aebe = spsarea->ae / spsarea->be; /* aebe = 1./(1.-fe) */
    spsarea->aebe2 = SQ(spsarea->aebe);
    spsarea->aebe3 = spsarea->aebe2 - 1.;
    w = SQ(spsarea->be / spsarea->ae);
    spsarea->aebe4 = SQ(w) - 1.; /*  aebe4 = (1.-fe)**4 - 1 */

    /*  Instrument-related constants */

    /*  Constants for conversions between different instrument coordinates */

    cilpevsc->iymax[0] = 6136;
    cilpevsc->iymax[1] = 2805;
    cilpevsc->ixmax[0] = 6136;
    cilpevsc->ixmax[1] = 2805;
    cilpevsc->alfincr[0] = 8.e-6;
    cilpevsc->alfincr[1] = 17.5e-6;
    cilpevsc->zetincr[0] = 16.e-6;
    cilpevsc->zetincr[1] = 35.e-6;
    cilpevsc->incrln[0] = 3.5;
    cilpevsc->incrln[1] = 16.;
    cilpevsc->incrpx[0] = 1.;
    cilpevsc->incrpx[1] = 8.;
    cilpevsc->alfln[0] = 28.e-6;
    cilpevsc->alfln[1] = 280.e-6;
    cilpevsc->zetpx[0] = 16.e-6;
    cilpevsc->zetpx[1] = 280.e-6;
    cilpevsc->alfmax[0] = (float) cilpevsc->iymax[0]*4.5 * cilpevsc->alfincr[0];
    cilpevsc->alfmax[1] = (float) cilpevsc->iymax[1]*4.5 * cilpevsc->alfincr[1];
    cilpevsc->zetmax[0] = (float) cilpevsc->ixmax[0]*2.5 * cilpevsc->zetincr[0];
    cilpevsc->zetmax[1] = (float) cilpevsc->ixmax[1]*2.5 * cilpevsc->zetincr[1];

    /* Set new offsets for the instruments */
    dx = 5.6 * spsarea->rad;
    dy = 3.6 * spsarea->rad;
    w = (cilpevsc->iymax[0] * NScyc1 + NSinc1) * cilpevsc->alfincr[0];
    if (fabs(w - cilpevsc->alfmax[0]) < dy)
        cilpevsc->alfmax[0] = w;
    w = (cilpevsc->ixmax[0] * EWcyc1 + EWinc1) * cilpevsc->zetincr[0];
    if (fabs(w - cilpevsc->zetmax[0]) < dx)
        cilpevsc->zetmax[0] = w;

    w = (cilpevsc->iymax[1] * NScyc2 + NSinc2) * cilpevsc->alfincr[1];
    if (fabs(w - cilpevsc->alfmax[1]) < dy)
        cilpevsc->alfmax[1] = 2. * cilpevsc->alfmax[1] - w;
    w = (cilpevsc->ixmax[1] * EWcyc2 + EWinc2) * cilpevsc->zetincr[1];
    if (fabs(w - cilpevsc->zetmax[1]) < dx)
        cilpevsc->zetmax[1] = w;
   
    return;
}

/***********************************************************************

     The function limb() defines the intersect of the scan plane 
     with the earth surface and the position of the limbic points.

     Arguments:
         instr = instrument code (1 - imager, 2 - sounder)
         flip = instrument orientation flag (1-nominal, (-1)-inverted)
         alf  = elevation angle in radians
         c    = coefficients of the equation defining the intersect
                of the scan plane with the earth surface
         st   = defines position of limbic points
         ierr = output status; 0 - successful completion,
                1 - scan plane does not intersect earth
         spsarea = pointer to DB

 **********************************************************************/

void limb(int instr, int flip, float alf, double c[][3], double *st,
        int *ierr, SPSAREA *spsarea) {
    double a0, a1, a2, a3, ab1, ab2, b1, ct, gn, l, m1, m2, w1, w2,
            g[3], g1[3], da, d, FF;

    *ierr = 1;
    FF = (double) flip;
    if (instr == 2) FF = -FF;
    /* da - mean elev. angle correction */
    da = FF * spsarea->pma * sin((double) alf);
    g[0] = -sin((double) alf + da);
    g[1] = 0;
    g[2] = cos((double) alf + da);
    g1[0] = spsarea->bt[0][0] * g[0] + spsarea->bt[1][0] * g[1] +
            spsarea->bt[2][0] * g[2];
    g1[1] = spsarea->bt[0][1] * g[0] + spsarea->bt[1][1] * g[1] +
            spsarea->bt[2][1] * g[2];
    g1[2] = spsarea->bt[0][2] * g[0] + spsarea->bt[1][2] * g[1] +
            spsarea->bt[2][2] * g[2];
    a0 = g1[0] * spsarea->xs[0] + g1[1] * spsarea->xs[1] + g1[2] * spsarea->xs[2];
    g1[2] = g1[2] / spsarea->aebe;
    gn = sqrt(SQ(g1[0]) + SQ(g1[1]) + SQ(g1[2]));
    a1 = g1[0] / gn;
    a2 = g1[1] / gn;
    a3 = g1[2] / gn;
    a0 = a0 / gn;
    d = 1. - a0*a0;
    if (d < 1.e-12)
        return; /* no intersections with earth */
    b1 = a1 * a1 + a2*a2;
    if (b1 < 1.e-7) {
        b1 = 0.;
        ab1 = 1.;
        ab2 = 0.;
    } else {
        b1 = sqrt(b1);
        ab1 = a1 / b1;
        ab2 = a2 / b1;
    }
    m1 = (ab1 * spsarea->xs[0] + ab2 * spsarea->xs[1]) * a3 -
            b1 * spsarea->xs[2] / spsarea->aebe;
    m2 = ab1 * spsarea->xs[1] - ab2 * spsarea->xs[0];
    l = sqrt(m1 * m1 + m2 * m2);
    d = sqrt(d);
    ct = d / l;
    w1 = ct*m1;
    w2 = ct*m2;
    c[0][0] = a3 * ab1 * w1 - ab2*w2;
    c[1][0] = -a3 * ab1 * w2 - ab2*w1;
    c[2][0] = a1*a0;
    c[0][1] = a3 * ab2 * w1 + ab1*w2;
    c[1][1] = -a3 * ab2 * w2 + ab1*w1;
    c[2][1] = a2*a0;
    c[0][2] = -b1 * w1 / spsarea->aebe;
    c[1][2] = b1 * w2 / spsarea->aebe;
    c[2][2] = a3 * a0 / spsarea->aebe;
    *st = sqrt(1. - ct * ct);
    *ierr = 0;
    return;
}

/************************************************************************

     The function clipper() computes the image clipping pixel boundaries.
     It returns 1, if the specified scan is off the earth, and 
     0 otherwise. 

     Arguments:
          instr     = instrument code (1 - imager, 2 - sounder)
          flip = instrument orientation flag (1-nominal, (-1)-inverted)
          cy, incy  = North-South scan address in cycles/increments 
          cw        = clipping wedge in radians 
          wpix      = western pixel number 
          epix      = eastern pixel number 
          spsarea   = pointer to DB
          cilpevsc  = pointer to DB

 *************************************************************************/

int clipper(int instr, int flip, int cy, int incy, float cw,
        float *wpix, float *epix, SPSAREA *spsarea, CILPEVSC *cilpevsc) {
    double ct, st, c[3][3], gam1, gam2, s1, u[3], f[3], fs[3], FF, Doff;
    int ierr;
    float alf;

    FF = (double) flip;
    Doff = FF * (cilpevsc->zetmax[instr - 1] -
            2.5 * (float) cilpevsc->ixmax[instr - 1]
            * cilpevsc->zetincr[instr - 1]);
    if (instr == 2) FF = -FF;
    alf = evci(instr, flip, cy, incy, 0, cilpevsc);
    limb(instr, flip, alf, c, &st, &ierr, spsarea);
    if (ierr != 0)
        return (1);

    /*  ECEF coordinates of the left limbic point */
    ct = sqrt(1. - st * st);
    u[0] = c[0][0] * ct - c[1][0] * st + c[2][0];
    u[1] = c[0][1] * ct - c[1][1] * st + c[2][1];
    u[2] = c[0][2] * ct - c[1][2] * st + c[2][2];

    /*  Scan angle at left limbic point */
    fs[0] = u[0] - spsarea->xs[0];
    fs[1] = u[1] - spsarea->xs[1];
    fs[2] = u[2] - spsarea->xs[2];
    f[0] = spsarea->bt[0][0] * fs[0] + spsarea->bt[0][1] * fs[1] +
            spsarea->bt[0][2] * fs[2];
    f[1] = spsarea->bt[1][0] * fs[0] + spsarea->bt[1][1] * fs[1] +
            spsarea->bt[1][2] * fs[2];
    f[2] = spsarea->bt[2][0] * fs[0] + spsarea->bt[2][1] * fs[1] +
            spsarea->bt[2][2] * fs[2];
    s1 = atan(f[2] / f[0]); /* elevation */
    gam1 = -asin(f[1] / sqrt(SQ(f[0]) + SQ(f[1]) + SQ(f[2])))
            - 0.5 * s1 * s1 * Doff - FF * spsarea->rma * sin(s1);

    /*   ECEF coordinates of the right limbic point */
    u[0] = c[0][0] * ct + c[1][0] * st + c[2][0];
    u[1] = c[0][1] * ct + c[1][1] * st + c[2][1];
    u[2] = c[0][2] * ct + c[1][2] * st + c[2][2];

    /*  Scan angle at right limbic point */
    fs[0] = u[0] - spsarea->xs[0];
    fs[1] = u[1] - spsarea->xs[1];
    fs[2] = u[2] - spsarea->xs[2];
    f[0] = spsarea->bt[0][0] * fs[0] + spsarea->bt[0][1] * fs[1] +
            spsarea->bt[0][2] * fs[2];
    f[1] = spsarea->bt[1][0] * fs[0] + spsarea->bt[1][1] * fs[1] +
            spsarea->bt[1][2] * fs[2];
    f[2] = spsarea->bt[2][0] * fs[0] + spsarea->bt[2][1] * fs[1] +
            spsarea->bt[2][2] * fs[2];
    s1 = atan(f[2] / f[0]); /* elevation */
    gam2 = -asin(f[1] / sqrt(SQ(f[0]) + SQ(f[1]) + SQ(f[2])))
            - 0.5 * s1 * s1 * Doff - FF * spsarea->rma * sin(s1);
    gam1 -= (double) cw;
    gam2 += (double) cw;
    *wpix = 1. + ((float) gam1 + cilpevsc->zetmax[instr - 1]) /
            cilpevsc->zetpx[instr - 1];
    *epix = 1. + ((float) gam2 + cilpevsc->zetmax[instr - 1]) /
            cilpevsc->zetpx[instr - 1];
    return (0);
}

/************************************************************************

    The function evln() converts fractional line number to 
    elevation angle in radians.
 
    Arguments:
             instr = instrument code (1-imager, 2-sounder) 
             flip = instrument orientation flag (1-nominal, (-1)-inverted)
             line  = fractional line number 

 ************************************************************************/
float evln(int instr, int flip, float rline, CILPEVSC *cilpevsc) {
    float w;
    //printf("alphaln=%f, rline= %f", cilpevsc->alfln[instr-1], rline);
    if (instr == 1)
        w = (rline - 4.5) * cilpevsc->alfln[instr - 1];
    else
        w = (rline - 2.5) * cilpevsc->alfln[instr - 1];

    if (flip == 1)
        return (cilpevsc->alfmax[instr - 1] - w);
    else
        return (9. * (float) cilpevsc->iymax[instr - 1] *
            cilpevsc->alfincr[instr - 1] - cilpevsc->alfmax[instr - 1] - w);
}

/************************************************************************

      sndeloc() accepts the mirror position in cycles and increments, 
      servo error values, and the positional offsets for four detectors 
      of a selected sounder channel and computes the detector earth 
      locations in latitude/longitude coordinates.
      
      Arguments:
         flip = instrument orientation flag (1-nominal, (-1)-inverted)
         cyew  = E-W cycles 
         incew = E-W increments   
         cyns  = N-S cycles   
         incns = N-S increments   
         svew  = E-W servo error in radians   
         svns  = N-S servo error in radians   
         doff  = offsets (rad) of 4 detectors with respect to their 
                 nominal positions (doff[0][0..3] - E-W offsets, 
                 doff[1][0..3] - N-S offsets)   
         geo   = geographic coordinates (rad) related to 4 detectors 
                 (geo[0][0..3] - latitude, geo[1][0..3] - longitude)
         cilpevsc, spsarea = pointers to DB

 *********************************************************************/
void sndeloc(int flip, int cyew, int incew, int cyns, int incns,
        float svew, float svns, float doff[][4], float geo[][4],
        CILPEVSC *cilpevsc, SPSAREA *spsarea) {
    int i, ier;
    float cose, de, ds, e, ev, h, s, sc, sine;

    /*  Convert the mirror position, given in cycles and increments,
        to elevation and scan angles. Add servo errors to obtain
        the true mirror gimbal angles */

    e = (((float) (cyns - 9) * cilpevsc->iymax[1] + incns)
            * cilpevsc->alfincr[1] + cilpevsc->alfmax[1] + svns)
            *(float) flip;
    s = ((float) (cyew * cilpevsc->ixmax[1] + incew)
            * cilpevsc->zetincr[1] - cilpevsc->zetmax[1] + svew)
            *(float) flip;
    sine = (float) sin((double) e)*(float) flip;
    cose = (float) cos((double) e);

    /* Compute positional offsets for each detector, convert them to
       angular offsets, correct elevation and scan angles, and then
       compute latitude and longitude of the corresponding points on
       the earth surface.
       Note: If a detector looks off the earth, the related latitude
             and longitude are set to 999999.  */

    h = -2. * cilpevsc->zetpx[1];
    for (i = 1; i <= 4; i++) {
        de = (2.5 - (float) i) * cilpevsc->alfln[1] + doff[1][i - 1];
        ds = h + doff[0][i - 1];
        ev = e + de * cose + ds*sine;
        sc = s - de * sine + ds*cose;
        lpoint(2, flip, ev, sc, &geo[0][i - 1], &geo[1][i - 1], &ier,
                spsarea, cilpevsc);
        h = -h;
    }
    return;
}

/**********************************************************************

   Integral Systems, Inc.

 ***********************************************************************

   Project   : Operations Ground Equipment for GOES
   System    : Sensors Processing System
   Function  : time50()
   Source    : time50.c
   Programmer: Igor Levine

   Ver.    Data    By   Comment
   ----  --------  ---  ---------------------------------------------
   1     03/02/93  IL   Initial creation

 ***********************************************************************

    time50() converts data and time given in BCD format to minutes
    from 1950 Jan. 1.0 UTC.

 **********************************************************************/

double time50(float x[]) {
    static int k, y, d, h, m, i;
    int extract_left4();
    static double s;
    static long l;
    static QUATER a[2];

    /* Convert year, day of year, hours, and minutes
       from the BCD format to integer values  */

    a[0].r4 = x[0];
    a[1].r4 = x[1];
    y = extract_left4(&a[0].u4);
    for (i = 1; i <= 3; i++)
        y = y * 10 + extract_left4(&a[0].u4);
    d = extract_left4(&a[0].u4);
    d &= 0x07;
    for (i = 1; i <= 2; i++)
        d = d * 10 + extract_left4(&a[0].u4);
    h = extract_left4(&a[0].u4);
    h = h * 10 + extract_left4(&a[1].u4);
    m = extract_left4(&a[1].u4);
    m = m * 10 + extract_left4(&a[1].u4);

    /* Convert seconds from the BCD format  */

    s = 0.;
    for (i = 1; i <= 5; i++) {
        k = extract_left4(&a[1].u4);
        s = s * 10. + (double) k;
    }

    //printf("year = %d, day = %d, hour = %d\n", y, d, h);
    /* Convert year and day of year to number of days
       from 0 hrs UTC, 1950 Jan. 1  */

    l = (long) d + 1461 * ((long) y + 4799) / 4 - 3 * (((long) y
            + 4899) / 100) / 4 - 2465022;
    //printf("jday = %d\n", l);
    /* Return time in minutes from 1950 1.0  */

    return (1440. * (double) l + 60. * (double) h + (double) m + s / 60000.);
}

int extract_left4(unsigned long *a) {
    static unsigned long Mask = 0Xf0000000;
    static int b;
    b = (int) ((*a & Mask) >> 28);
    *a = *a << 4;
    return (b);
}

/************************************************************************

     tupack50 converts date and time given in BCD format to integer
     day of year, hours, minutes, seconds, milliseconds, and minutes
     from 1950 Jan. 1.0 UTC. 

     Arguments:
          x[2] = date/time in BCD format 
                 ( x[0] = YYYY DDDH, x[1] = HMMSSLLL ) 
          yr = year
          doy = day of year
          hrs = hours
          mns = minutes
          sec = seconds
          msec = milliseconds
          t50 = time in minutes since 1950, Jan. 1.0 UTC.

 ***********************************************************************/
void tupack50(float x[], int *yr, int *doy, int *hrs, int *mns,
        int *sec, int *msec, double *t50) {
    static int k, y, d, h, m, i;
    int extract_left4();
    static double s;
    static long l;
    static QUATER a[2];

    /* Convert year, day of year, hours, and minutes
       from the BCD format to integer values  */

    a[0].r4 = x[0];
    a[1].r4 = x[1];
    y = extract_left4(&a[0].u4);
    for (i = 1; i <= 3; i++)
        y = y * 10 + extract_left4(&a[0].u4);
    d = extract_left4(&a[0].u4);
    d &= 0x07;
    for (i = 1; i <= 2; i++)
        d = d * 10 + extract_left4(&a[0].u4);
    h = extract_left4(&a[0].u4);
    h = h * 10 + extract_left4(&a[1].u4);
    m = extract_left4(&a[1].u4);
    m = m * 10 + extract_left4(&a[1].u4);

    /* Convert seconds from the BCD format  */

    s = 0.;
    for (i = 1; i <= 5; i++) {
        k = extract_left4(&a[1].u4);
        s = s * 10. + (double) k;
    }

    /* Convert year and day of year to number of days
       from 0 hrs UTC, 1950 Jan. 1  */

    l = (long) d + 1461 * ((long) y + 4799) / 4 - 3 * (((long) y
            + 4899) / 100) / 4 - 2465022;

    *t50 = 1440. * (double) l + 60. * (double) h + (double) m + s / 60000.;
    *yr = y;
    *doy = d;
    *hrs = h;
    *mns = m;
    *sec = s / 1000;
    *msec = s - *sec * 1000;
    return;
}

/************************************************************************

     utime50() converts double precision minutes from 1950 Jan. 1.0
     to integer year, day of year, hour, minutes, seconds and
     milliseconds.

     Arguments:
           t = time in minutes since Jan. 1.0, 1950.
           yr = year
           doy = day of year
           hrs = hours
           mns = minutes
           sec = seconds
           msec = milliseconds

     Note:
       This conversion is correct only for data between Jan. 1, 1900  
       and  Dec. 31, 2100.

     ver.    date    by                comment
     ----  --------  ----------------  --------------------------------
     1     12/29/88  I. Levine         Initial creation.
     2     04/17/95  I. Levine         Converted from FORTRAN to C.

 ***********************************************************************/
void utime50(double t, int *yr, int *doy, int *hrs, int *mns,
        int *sec, int *msec) {
    long jd, l, ldoy, lyr;
    float s;
    double r;

    /* Extract number of days from 0 hour UT, Jan.1, 1950 */
    jd = (long) (t / 1440.);
    r = t - (double) jd * 1440.;

    /* Extract hour of day */
    *hrs = (int) (r / 60.);
    r = r - (double) (*hrs)*60.;

    /* minutes */
    *mns = (int) r;

    /* seconds */
    s = (r - (double) (*mns))*60.;
    *sec = (int) s;

    /* millseconds */
    *msec = (int) (1000. * (s - *sec) + 0.5);

    /* Convert Julian day from 1950, Jan 1.0 to year and day of year */
    ldoy = jd + 365;
    l = ldoy / 1461;
    lyr = 1949 + l * 4;
    ldoy = ldoy - l * 1461;
    if (ldoy > 1095) l = 3;
    else l = ldoy / 365;
    *yr = lyr + l;
    *doy = ldoy - l * 365 + 1;
    return;
}

/************************************************************************

     px2ci() accepts integer pixel number and computes the related 
     E-W pointing of the instrument expressed in cycles/increments.

     Arguments:
        instr = instrument code (1-Imager, 2-Sounder)
        flip = instrument orientation flag (1-nominal, (-1)-inverted)
        px = pixel number
        cx = cycles
        incr = increments

 ************************************************************************/

void px2ci(int instr, int flip, int px, int *cx, int *incx,
        CILPEVSC *cilpevsc) {
    int i;
    if (flip == 1)
        i = (px - 1) * cilpevsc->incrpx[instr - 1] + 0.5;
    else
        i = 5 * cilpevsc->ixmax[instr - 1] -
            px * cilpevsc->incrpx[instr - 1] + 1.5;
    *cx = i / cilpevsc->ixmax[instr - 1];
    *incx = i - *cx * cilpevsc->ixmax[instr - 1];
    return;
}

/***********************************************************************

     evsc2lp() converts elevation and scan angles
     to the nearest integer line and pixel numbers.

     Arguments:
        instr = instrument code
        flip = instrument orientation flag (1-nominal, (-1)-inverted)
        alpha = elevation (N/S) angle in radians
        zeta = scan (E/W) angle in radians
        line = line number
        ipix = pixel number

 ***********************************************************************/
void evsc2lp_double(int instr, int flip, float alpha, float zeta, double *line,
        double *ipix, CILPEVSC *cilpevsc) {
    float rl, rp;
    if (instr == 1) rl = 4.5;
    else rl = 2.5;
    if (flip == 1) {
        rl += (cilpevsc->alfmax[instr - 1] - alpha)
                / cilpevsc->alfln[instr - 1];
        rp = (cilpevsc->zetmax[instr - 1] + zeta) / cilpevsc->zetpx[instr - 1] + 1.;
    } else {
        rl += 9. * (float) cilpevsc->iymax[instr - 1] / cilpevsc->incrln[instr - 1]
                - (cilpevsc->alfmax[instr - 1] + alpha) / cilpevsc->alfln[instr - 1];
        rp = 5. * (float) cilpevsc->ixmax[instr - 1] / cilpevsc->incrpx[instr - 1]
                - (cilpevsc->zetmax[instr - 1] - zeta) / cilpevsc->zetpx[instr - 1] + 1.;
    }
    *line=rl;
    *ipix=rp;
    //    *line = rl + 0.5 * sign(rl);
    //*ipix = rp + 0.5 * sign(rp);
    return;
}

/************************************************************************

     Integral Systems, Inc.

 ***********************************************************************

     Project   : Operations Ground Equipment for GOES
     System    : Sensors Processing System
     Functions : slant_range.
     Source    : slant.c
     Programmer: Igor Levine

     ver.    data    by   comment
     ----  --------  ---  ---------------------------------------------
     1     05/01/95  IL   Initial creation.
 ***********************************************************************

     The function slant_range() computes the slant distance (in km) from 
     a given point on the earth to the satellite. The distance is negative
     if the point is invisible from the satellite.

     Arguments:
        rlat    =  geographic latitude in radians (input), 
        rlon    =  geographic longitude in radians (input),
        h       =  altitude in meters (input), 
        spsarea = pointer to DB structure.

 ************************************************************************/

double slant_range(float rlat, float rlon, float h, SPSAREA *spsarea) {
    double sing, slat, u[3], w2, f[3], r, rloc, H, latc, diff, lath;

    /*  Sinus of geographic (geodetic) latitude */
    sing = sin((double) rlat);

    /*  Sinus of the geocentric latitude */
    slat = sing / (spsarea->aebe2 * sqrt(1. + spsarea->aebe4 * sing * sing));
    latc = asin(slat); /* geocentric latitude */
    diff = (double) rlat - latc; /* geographic - geocentric lat */

    /*  Computes local earth radius at specified point */
    rloc = 1. / sqrt(1. + spsarea->aebe3 * slat * slat);

    /* Compute the geocentric latitude and radius-vector of the point
       with given altitude */

    H = 0.001 * h / spsarea->ae;
    r = sqrt(H * H + rloc * rloc + 2. * H * rloc * cos(diff));
    lath = latc + asin(H * sin(diff) / r);

    /*  Computes Cartesian coordinates of the point */
    u[2] = r * sin(lath);
    w2 = r * cos(lath);
    u[0] = w2 * cos((double) rlon);
    u[1] = w2 * sin((double) rlon);

    /*  Pointing vector from satellite to the earth point */
    f[0] = u[0] - spsarea->xs[0];
    f[1] = u[1] - spsarea->xs[1];
    f[2] = u[2] - spsarea->xs[2];
    w2 = u[0] * f[0] + u[1] * f[1] + u[2] * f[2] * spsarea->aebe2;

    /*  Verify visibility of the point */
    if (w2 > 0.)
        return (-1.); /* invisible point on the earth */

    /* Compute slant distance */

    w2 = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    return (spsarea->ae * w2);
}

/***********************************************************************

     ev2ln() converts elevation angle in radians to the absolute
     line number.

     Arguments:
        instr = instrument code
        flip = instrument orientation flag (1-nominal, (-1)-inverted)
        alpha = elevation (N/S) angle in radians

 ***********************************************************************/
float ev2ln(int instr, int flip, float alpha, CILPEVSC *cilpevsc) {
    float rl;
    if (instr == 1) rl = 4.5;
    else rl = 2.5;
    if (flip == 1)
        return (rl + (cilpevsc->alfmax[instr - 1] - alpha) / cilpevsc->alfln[instr - 1]);
    else {
        rl += 9. * (float) cilpevsc->iymax[instr - 1] * cilpevsc->alfincr[instr - 1];
        return (rl - (cilpevsc->alfmax[instr - 1] + alpha) / cilpevsc->alfln[instr - 1]);
    }
}

/***********************************************************************

     retateZ() perform the coordinate rotation around z-axis.

     Arguments:
        v = the input 3-dim vector
        angle = the rotation angle.

 ***********************************************************************/

void rotateZ(double *v, double angle) {
    double p1, p2, p3;
    p1 = *v;
    p2 = *(v + 1);
    p3 = *(v + 2);
    *v = cos(angle) * p1 + sin(angle) * p2;
    *(v+1) = -sin(angle) * p1 + cos(angle) * p2;
    *(v+2) = p3;
}

/***********************************************************************

     retateY() perform the coordinate rotation around y-axis.

     Arguments:
        v = the input 3-dim vector
        angle = the rotation angle.

 ***********************************************************************/

void rotateY(double *v, double angle) {
    double p1, p2, p3;
    p1 = *v;
    p2 = *(v+1);
    p3 = *(v+2);
    *v = cos(angle) * p1 + sin(angle) * p3;
    *(v+1) = p2;
    *(v+2) = -p1 * sin(angle) + p3 * cos(angle);
}

/***********************************************************************

     retateX() perform the coordinate rotation around x-axis.

     Arguments:
        v = the input 3-dim vector
        angle = the rotation angle.

 ***********************************************************************/

void rotateX(double *v, double angle) {
    double p1, p2, p3;
    p1=*v;
    p2 = *(v+1);
    p3 = *(v+2);
    *v  = p1;
    *(v+1) = p2 * cos(angle) - p3 * sin(angle);
    *(v+2) = p2 * sin(angle) + p3 * cos(angle);
}

