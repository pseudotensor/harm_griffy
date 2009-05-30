//============================= GRFFDE ===========================//

/*
 * This steps primitive variables forward in time
 */

#include "decs.h"

void step_ch()
{

	/* halfstep */
	fprintf(stderr,"half ") ;
	advance(p, p, 0.5*dt, ph) ;

	/* full step */
	fprintf(stderr," full ") ;
	advance(p, ph, dt,    p) ;

        /* increment time */
        t += dt ;
	if(t + dt > tf) dt = tf - t ;  /* don't step beyond end of run */

        /* done! */
}

void advance(
	double pi[INA][JNA][NPR],   /* initial prims in */
	double pb[INA][JNA][NPR],   /* prims used to calc fluxes */
	const double Dt,            /* timestep for "advance" */
	double pf[INA][JNA][NPR]    /* final prims out */
	)
{
	int i,j,k ;

	double U[NPR],dU[NPR] ;

#ifdef DEBUG
	double dU_debug[INA][JNA][NPR] ;
#endif /* DEBUG */

	struct of_state state ;
	struct of_geom geom ;

	/* update fluxes using pb */
	fluxcalc(pb, F1, 1) ;
	fluxcalc(pb, F2, 2) ;

#ifdef DEBUG
	bdump_fluxes(1);
#endif /* DEBUG */

	/* impose axial BCs on the flux */
	fix_flux(F1, F2);

#ifdef DEBUG
	bdump_fluxes(2);
#endif /* DEBUG */
	/* CT constraint */
	flux_ct(F1,F2) ;

#ifdef DEBUG
	bdump_fluxes(3);
	PLOOP ZLOOP dU_debug[i][j][k]=-999.;
#endif /* DEBUG */

	/* now update pi to pf */
	ZLOOP {

	        /* get local geometry */
		get_geometry(i,j,CENT,&geom) ;

		get_state(pb[i][j],&state,&geom) ;

                /* calculate source term (stress-energy T * connection) */
		source(pb[i][j],dU,i,j,&state,&geom) ;
#ifdef DEBUG
		PLOOP dU_debug[i][j][k]=dU[k] ;
#endif /* DEBUG */
		/* pi to cons */
		primtoU(pi[i][j],U,&geom) ;

		/* update cons using new fluxes */

		PLOOP {
		       U[k] += Dt*( 
		  		   - (F1[i+1][j][k] - F1[i][j][k])/dx[1] 
		  		   - (F2[i][j+1][k] - F2[i][j][k])/dx[2] 
		  		   + dU[k]
				   ) ; 

		}  

		/* update prims using new cons */
                UtoPrim(pf[i][j],U,&geom) ; // This alters U

	}

#ifdef DEBUG
	bdump_source(dU_debug);
#endif /* DEBUG */

	/* apply boundary conditions */
        bound_prim(pf) ;

	return ;
}


void fluxcalc(
	      double pr[INA][JNA][NPR],    /* prims used to calc flux */
	      double F[INA][JNA][NPR],     /* Fluxes */
	const int dir                      /* direction? */
	)
{
	int i,j,k,idel,jdel,face ;
	double p_l[NPR],p_r[NPR],F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ;
	double ctop ;

	struct of_state state_l,state_r ;
	struct of_geom geom ;

	/* which direction? */
        if     (dir == 1) {idel = 1; jdel = 0; face = FACE1;}
	else if(dir == 2) {idel = 0; jdel = 1; face = FACE2;}
	else { exit(10); }

	/** evaluate slopes of primitive variables **/
	ZSLOOP(IS-1,IE+1,JS-1,JE+1) PLOOP {
			dq[i][j][k] = slope_lim(
					pr[i-idel][j-jdel][k],
					pr[i][j][k],
					pr[i+idel][j+jdel][k]
					) ;
	}

	/* calculate left/right prims */
        ZSLOOP(IS-jdel,IE+1,JS-idel,JE+1) {
	  /* this avoids problems on the pole */
	  if(dir == 2 && (j == JS || j == JE+1)) {
	    PLOOP F[i][j][k] = 0. ;
	  }
	  else {


                PLOOP {
                        p_l[k] = pr[i-idel][j-jdel][k] 
			       + 0.5*dq[i-idel][j-jdel][k] ;
                        p_r[k] = pr[i][j][k]   
			       - 0.5*dq[i][j][k] ;
                }

		/* get geometry depending on face */
		get_geometry(i,j,face,&geom) ;

		get_state(p_l,&state_l,&geom) ;
		get_state(p_r,&state_r,&geom) ;

		/* convert left/right prims to left/right fluxes */
		primtoflux(p_l,F_l,&state_l,&geom,dir) ;
		primtoflux(p_r,F_r,&state_r,&geom,dir) ;

		ctop=vchar(i,j,dir);

		/* convert left/right prims to left/right cons*/
		primtoU(p_l,U_l,&geom) ;
		primtoU(p_r,U_r,&geom) ;

		/* use LAX technique */
		PLOOP F[i][j][k] =0.5*(F_l[k] + F_r[k]
				     + ctop*(U_l[k] - U_r[k])) ;
	  }
	}

	return ;

}


void flux_ct(double F1_local[INA][JNA][NPR], double F2_local[INA][JNA][NPR])
{
	int i,j ;
	static double emf[INA][JNA] ;

	/* calculate EMFs */
	/* Toth approach: just average */
	ZSLOOP(IS,IE+1,JS,JE+1)
	    emf[i][j] = 0.25*(F1_local[i][j][M02] + F1_local[i][j-1][M02]
			      - F2_local[i][j][M01] - F2_local[i-1][j][M01]) ;

	/* rewrite EMFs as fluxes, after Toth */
        ZSLOOP(IS,IE+1,JS,JE) {
                F1_local[i][j][M01] = 0. ;
                F1_local[i][j][M02] =  0.5*(emf[i][j] + emf[i][j+1]) ;
        }
        ZSLOOP(IS,IE,JS,JE+1) {
                F2_local[i][j][M01] = -0.5*(emf[i][j] + emf[i+1][j]) ;
                F2_local[i][j][M02] = 0. ;
	}

}

