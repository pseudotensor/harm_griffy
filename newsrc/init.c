//============================== GRFFDE ==============================//

/*
 * sets initial conditions in KS' coordinates
 * for a magnetic monopole
 *
 */

#include "decs.h"

#define SMALL 1e-9

void init()
{
	int i,j ;
	double r,r_horizon,th,Bi ;
	double X[NDIM] ;

	struct of_geom geom ;

        /* some numerical parameters */
        lim = MC ;      /* choose limiter */
        cour = 1.1  ;   /* courant # */

	/* for KS prime coordinates */
	R0 = 0. ;
	hslope = 1. ;    

	a = 0.9375 ;                          /* angular momentum */
  	r_horizon = (1.+sqrt(1.-a*a)) ;    /* horizon radius */
	i_horizon = IS+4 ;                 /* i-value at horizon */

	Bi = 1. ;                          /* initial B (arbitrary) */
	Rout = 200. ;                      /* outer radius of grid */
	//Rin = 0.9 * r_horizon ;          

	/* output parameters */
	DTd = 20.0  ;	/* dumping frequency */
	DTb = 1.0   ;   /* binary dump frequency */
	DTl = 0.1   ;	/* logfile frequency */

	/* solve for Rin given choices above */
	get_startx() ;

	/* Set up all the geometry */
        set_arrays() ;
        set_grid() ;

	/* tick-tock */
        t = 0. ;                            /* initial time */
	tf = 300.0 ;                        /* final time */

	/* print info in terminal */
	fprintf(stderr,"\n\nUsing Modified Kerr-Schild Coordinates\n\n") ;

	coord(ISG,JS,CORN,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"rmin (corner): %10.5g\n",r) ;
	coord(ISG,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"rmin (center): %10.5g\n",r) ;

	coord(IS,JS,CORN,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"First active r (corner): %10.5g\n",r) ;
	coord(IS,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"First active r (center): %10.5g\n",r) ;

	coord(IE,JS,CORN,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"Last active r (corner): %10.5g\n",r) ;
	coord(IE,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"Last active r (center): %10.5g\n",r) ;

	coord(IEG,JS,CORN,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"rmax (corner): %10.5g\n",r) ;
	coord(IEG,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"rmax (center): %10.5g\n",r) ;

	coord(IS+N1/4,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"r (quater-grid fluxes): %10.5g\n",r) ;

	coord(IS+N1/2,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"r (midgrid fluxes): %10.5g\n",r) ;

	coord(IE-2*NGHOSTS,JS,CENT,X) ;
	bl_coord(X,&r,&th) ;
	fprintf(stderr,"r (outer fluxes): %10.5g\n",r) ;

	fprintf(stderr,"Rin: %10.5g  Rout:  %10.5g\n",Rin,Rout) ;
	fprintf(stderr,"r_horizon: %g\n",r_horizon) ;
	fprintf(stderr,"i_horizon-IS: %d\n",i_horizon-IS) ;
	fprintf(stderr,"rmin/r_horizon: %g\n",r/r_horizon) ;

 	dt = timestep();
	fprintf(stderr,"dt: %10.5g\n",dt) ;
	fprintf(stderr,"Number of timesteps needed to finish: %d\n",
		(int) ceil((tf-t)/dt) ) ;

	DTl = 0.1*dt   ;	/* logfile frequency (every timestep) */

	/*
	 * geting down to business: build our monopole
	 */

	ZLOOP {

	    /* load "cartesian" grid coordinates */
	    coord(i,j,CENT,X) ;
	    /* load bl coordinates */
	    bl_coord(X,&r,&th) ;

	    /* geometric business */
	    get_geometry(i,j,CENT,&geom) ;

	    //====================
	    // check to make sure grid set up properly
	    if(i == i_horizon) 
	      if(fabs(r - r_horizon)> SMALL ){
		fprintf(stderr,"grid not set up properly!\n") ;
		fprintf(stderr,"r:  %10.5g  r_horizon:  %10.5g\n",r,r_horizon);
	      }
	    //
	    //====================

	    /* the easy prims */

	    p[i][j][F10] = 0. ; /* E-field */
	    p[i][j][F20] = 0. ;
	    p[i][j][F30] = 0. ;

	    p[i][j][M02] = 0. ; /* only radial B-field */
	    p[i][j][M03] = 0. ;

	    /* radial magnetic field component */
	    p[i][j][M01] = Bi*gdet[i_horizon][j][CENT]/(gdet[i][j][CENT]) ;

	    /* done ! */

	}

	/* enforce boundary conditions */
	bound_prim(p) ;

}

#undef SMALL

/*
 * this funtion assurs that the horizon is on
 * an integer 'i' value on the centered grid -- it sets based on 
 * the desired i_horizon value set above
 */

void get_startx()
{

      double i_horlocal,r_hor,k,b ;

      i_horlocal = (double)(i_horizon-IS) ;
      r_hor = (1. + sqrt(1. - a*a)) ;
      
      k = (i_horlocal + 0.5)/N1 ;
      b = 1. - k ;


      Rin = R0 + exp( (log(r_hor-R0) - k*log(Rout-R0))/b) ;

      return ;

}

/* set up all grid functions */
void set_grid()
{
	int i,j ;
	double X[NDIM] ;
	struct of_geom geom ;

	/* set up boundaries, steps in coordinate grid */
	set_points() ;

	DLOOPA X[j] = 0. ;

	ZSLOOP(ISG,IEG,JSG,JEG) {
		
		/* zone-centered */
		coord(i,j,CENT,X) ;
		gcov_func(X, gcov[i][j][CENT]) ;
		gdet[i][j][CENT] = gdet_func(X,gcov[i][j][CENT]) ;
		gcon_func(X,gcov[i][j][CENT], gcon[i][j][CENT]) ;

		get_geometry(i,j,CENT,&geom) ;
		conn_func(X, &geom, conn[i][j]) ;

		/* corner-centered */
		coord(i,j,CORN,X) ;
		gcov_func(X,gcov[i][j][CORN]) ;
		gdet[i][j][CORN] = gdet_func(X, gcov[i][j][CORN]) ;

		gcon_func(X, gcov[i][j][CORN],gcon[i][j][CORN]) ;

		/* r-face-centered */
		coord(i,j,FACE1,X) ;
		gcov_func(X,gcov[i][j][FACE1]) ;
		gdet[i][j][FACE1] = gdet_func(X, gcov[i][j][FACE1]) ;
		gcon_func(X, gcov[i][j][FACE1],gcon[i][j][FACE1]) ;

		/* theta-face-centered */
		coord(i,j,FACE2,X) ;
		gcov_func(X,gcov[i][j][FACE2]) ;
		gdet[i][j][FACE2] = gdet_func(X, gcov[i][j][FACE2]) ;
		gcon_func(X, gcov[i][j][FACE2],gcon[i][j][FACE2]) ;
	}

	/* done! */
}

void set_arrays()
{
	int i,j,k ;

	/* everything must be initialized to zero */
	ZSLOOP(ISG,IEG,JSG,JEG) {
		PLOOP {
			p [i][j][k] = 0. ;
			ph[i][j][k] = 0. ;
			dq[i][j][k] = 0. ;
			F1[i][j][k] = 0. ;
			F2[i][j][k] = 0. ;
		}
	}
}
