//========================= GRFFDE =========================//

double  p [INA][JNA][NPR] ;  /* space for prims */
double  dq[INA][JNA][NPR] ;  /* slopes */
double  F1[INA][JNA][NPR] ;  /* fluxes */
double  F2[INA][JNA][NPR] ;  /* fluxes */
double  ph[INA][JNA][NPR] ;  /* half-step prims */

/* for debug */

/* grid functions */
double conn[INA][JNA][NDIM][NDIM][NDIM] ;
double gcon[INA][JNA][NPG][NDIM][NDIM] ;
double gcov[INA][JNA][NPG][NDIM][NDIM] ;
double gdet[INA][JNA][NPG] ;

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
double a ;

/* numerical parameters */
int i_horizon ;

double Rin,Rout,hslope,R0 ;
double cour ;
double dx[NDIM],startx[NDIM] ;
double dt ;
double t,tf ;
int istart,istop,jstart,jstop ;
double dminarg1,dminarg2 ;
int nstep ;

/* output parameters */
double DTd ;
double DTl ;
double DTb ;
int    DTr ;
int    dump_cnt ;
int    bin_cnt ;
int    rdump_cnt ;

/* global flags */
int    lim ;

/* current local position */
int icurr,jcurr ;
