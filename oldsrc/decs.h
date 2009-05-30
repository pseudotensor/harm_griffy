//====================== GRFFDE =====================//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* number of grid zones */
#define N1	64
#define N2	64

/** MNEMONICS SECTION **/

/* mnemonics for primitive vars; conserved vars */
#define M01	0	/* Maxwell components M^ti*/
#define M02	1
#define M03	2
#define F10	3       /* Faraday components F^it*/
#define F20	4
#define F30	5

/* mnemonics for dimensional indices */
#define TT	0	
#define RR	1
#define TH	2
#define PH	3

/* mnemonics for centering of grid functions */
#define FACE1	0	
#define FACE2	1
#define CORN	2
#define CENT	3

/* mnemonics for slope limiter */
#define MC	0
#define VANL	1
#define MINM	2

/* mnemonics for diagnostic calls */
#define INIT_OUT	0
#define DUMP_OUT	1
#define IMAGE_OUT	2
#define LOG_OUT		3
#define FINAL_OUT	4
#define BINARY_OUT	5

/* mnemonics for termination conditions */

#define END_BINARY_OUTPUT_OPENING_FILE 150

/** GLOBAL ARRAY SECTION **/

/* size of global arrays */
#define NPR	6	/* number of primitive variables */
#define NDIM	4	/* number of total dimensions.  Never changes */
#define NPG	4	/* number of positions on grid for grid functions */

/* indices and sizes of field arrays */
#define NGHOSTS (2)              /* number of ghost zones on each side */
#define IS      (NGHOSTS)        /* first active zone */
#define IE      (N1+IS-1)        /* last active zone */
#define ISA     (0)              /* beginning of the array in this direction */
#define IEA     (N1+2*NGHOSTS-1) /* end of the array in this direction */
#define ISG     (IS-NGHOSTS)     /* first ghost zone */
#define IEG     (IE+NGHOSTS)     /* last ghost zone */
#define IN      (IE-IS +1)       /* number of active zones in this direction */
#define INA     (IEA-ISA +1)     /* array extent in this direction */
#define ING     (IEG-ISG +1)     /* total number of zones in this direction */

#define JS      (NGHOSTS)        /* first active zone */
#define JE      (N2+JS-1)        /* last active zone */
#define JSA     (0)              /* beginning of the array in this direction */
#define JEA     (N2+2*NGHOSTS-1) /* end of the array in this direction */
#define JSG     (JS-NGHOSTS)     /* first ghost zone */
#define JEG     (JE+NGHOSTS)     /* last ghost zone */
#define JN      (JE-JS +1)       /* number of active zones in this direction */
#define JNA     (JEA-JSA +1)     /* array extent in this direction */
#define JNG     (JEG-JSG +1)     /* total number of zones in this direction */

extern double  p [INA][JNA][NPR] ;	/* space for primitive vars */
extern double  dq[INA][JNA][NPR] ;      /* slopes */
extern double  F1[INA][JNA][NPR] ;	/* fluxes */
extern double  F2[INA][JNA][NPR] ;	/* fluxes */
extern double  ph[INA][JNA][NPR] ;	/* half-step primitives */

/* for debug -- not all used right now */

/* grid functions */
extern double conn[INA][JNA][NDIM][NDIM][NDIM] ; /* connection */
extern double gcon[INA][JNA][NPG][NDIM][NDIM] ;  /* g contra */
extern double gcov[INA][JNA][NPG][NDIM][NDIM] ;  /* g cov */
extern double gdet[INA][JNA][NPG] ;              /* det g */

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters*/
extern double a ;

/* numerical parameters */

extern int i_horizon ;            /* i-value at horizon */
extern double hslope, R0 ;
extern double Rin, Rout ;         /* for KS prime coord */
extern double cour ;
extern double dx[NDIM],startx[NDIM] ;
extern double dt ;
extern double t,tf ;
extern double x1curr,x2curr ;
extern int nstep ;

/* output parameters-'restart' routine not used now */
extern double DTd ;
extern double DTl ;
extern double DTb ;
extern int    DTr ;
extern int    dump_cnt ;
extern int    bin_cnt ;
extern int    rdump_cnt ;

/* global flags */
extern int lim ;

/** MACROS **/

/* loop over all active zones */
#define ILOOP for(i=IS;i<=IE;i++)
#define JLOOP for(j=JS;j<=JE;j++)
#define ZLOOP ILOOP JLOOP
#define IJLOOP ILOOP JLOOP

#define ZLOOPG for(i=ISG;i<=IEG;i++)for(j=JSG;j<=JEG;j++)

/* specialty loops */
extern int istart,istop,jstart,jstop ;
#define ZSLOOP(istart,istop,jstart,jstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)

/* loop over Primitive variables */
#define PLOOP  for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)

/* set global variables that indicate current local metric, etc. */
extern int icurr,jcurr,pcurr ;

struct of_geom {
	double gcon[NDIM][NDIM] ;
	double gcov[NDIM][NDIM] ;
	double g ;
} ;

struct of_state {
        double Mcon[NDIM][NDIM] ;
        double T[NDIM][NDIM] ;
} ;

//=========== function declarations ==========*/

/** single use/file functions **/
void init(void) ;
void ener(FILE *fp) ;
void get_startx(void) ;
double timestep(void) ;
double max2(double xx, double yy) ;
double min2(double xx, double yy) ;
double vchar(const int i, const int j, const int dir) ;
void dump(FILE *fp) ;
void gridoutput(char *binnam) ;

void bound_prim(double pr[INA][JNA][NPR]) ;
void set_arrays(void) ;
void set_grid(void) ;
double slope_lim(const double yy1, const double yy2, const double yy3) ;

/* bdump.c */
void bdump_primitive(const char *binnam, const int field) ;
void bdump_gdet(const char *binnam, const int loc) ;
void bdump_gcon(const char *binnam, const int loc, const int j, const int k) ;
void bdump_gcov(const char *binnam, const int loc, const int j, const int k) ;
void bdump_conn(const char *binnam, const int i, const int j, const int k) ;
void bdump_array(const char *binnam, const int n, double a[n],
		 const int begin, const int end);
void bdump_field(const char *binname, double field[INA][JNA]) ;
void bdump_fluxes(const int label) ;
void bdump_source(const double dU_debug[INA][JNA][NPR] ) ;

/* coord.c */
void bl_coord(const double X[NDIM], double *r, double *th) ;
void gcov_func(const double X[NDIM], double gcov_local[NDIM][NDIM]) ;
void set_points(void) ;
void fix_flux(double F1_local[INA][JNA][NPR], double F2_local[INA][JNA][NPR]);

/* diag.c  */
void diag(const int call_code) ;

/* metric.c */
void coord(const int i, const int j, const int loc, double X[NDIM]) ;
#ifdef NUMERICAL
/* Linear Algebra */
void gaussj_invert(const int n, double array[n][n]) ; /* to get gcon */
/* to get determinant of gcov */
void ludcmp(const int n,double array[n][n],int indx[n],double *d);
#endif /*  NUMERICAL */

/* metric.c or coord.c */
double gdet_func(const double X[NDIM], double lgcov[NDIM][NDIM]) ;
void gcon_func(const double X[NDIM],
	       double lgcov[NDIM][NDIM], double lgcon[NDIM][NDIM]) ;
void conn_func(const double X[NDIM], const struct of_geom *geom,
	       double lconn[NDIM][NDIM][NDIM]) ;

/*phys.c */
void primtoU(const double pr[NPR], double U[NPR], const struct of_geom *geom) ;
void primtoflux(const double pr[NPR], double f[NPR],
		const struct of_state *state, 
		const struct of_geom *geom, const int dir) ;
void stresstensor(double Mcon[NDIM][NDIM],
		  double T[NDIM][NDIM], const struct of_geom *geom) ;
void source(const double ph_l[NPR], double dU[NPR],
	    const int ii, const int jj, 
	    const struct of_state *state, const struct of_geom *geom) ;
void UtoPrim(double pr[NPR], double U[NPR], const struct of_geom *geom) ;
void Max_con(const double pr[NPR],
	     double Mcon[NDIM][NDIM], const struct of_geom *geom) ;
void lower_v(const double vcon[NDIM], double vcov[NDIM],
	     const struct of_geom *geom) ;
void lower_T(double Tcon[NDIM][NDIM],
	     double Tcov[NDIM][NDIM], const struct of_geom *geom) ;
void get_geometry(const int ii, const int jj, const int kk,
		  struct of_geom *geom) ;
void get_state(const double pr[NPR], 
	       struct of_state *state, const struct of_geom *geom) ;

/* step_ch.c */
void step_ch(void) ;
void advance( double pi[INA][JNA][NPR], double pb[INA][JNA][NPR], 
	      const double Dt, double pf[INA][JNA][NPR]) ;
void flux_ct(double F1_l[INA][JNA][NPR],double F2_l[INA][JNA][NPR]) ;
void fluxcalc( double pr[INA][JNA][NPR], double F[INA][JNA][NPR], 
	       const int dir ) ;
