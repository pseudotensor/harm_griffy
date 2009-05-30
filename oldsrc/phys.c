//============================ GRFFDE =============================//

#include "decs.h"

/* 
 * convert primitives to conserved quantities 
 */

void primtoU(const double pr[NPR], double U[NPR], const struct of_geom *geom)
{

      int k ;
      double F1t,F2t,F3t,Mt1,Mt2,Mt3 ;

      /* load prims */
      F1t = pr[F10] ;   /* F^it -- E components */
      F2t = pr[F20] ;
      F3t = pr[F30] ;
      Mt1 = pr[M01] ;   /* M^ti -- B components */
      Mt2 = pr[M02] ;
      Mt3 = pr[M03] ;

      /* M^ti -- just copy!*/
      U[0] = Mt1 ;
      U[1] = Mt2 ;
      U[2] = Mt3 ;

      /* T^t_i terms( t up i down):
       *
       * T^t_i = - sqrt(-g) * [ijk] * M^tj * F^kt
       *
       */

      U[3] = - geom->g * (Mt2*F3t - Mt3*F2t) ;
      U[4] = - geom->g * (Mt3*F1t - Mt1*F3t) ;
      U[5] = - geom->g * (Mt1*F2t - Mt2*F1t) ;

      /* mult by geometric factor */
      PLOOP U[k] *= geom->g ;
}

/* 
 * convert prims to fluxes
 */
/* pr unused here.  Check if this is correct */
void primtoflux(const double pr[NPR], double f[NPR],
		const struct of_state *state, 
		const struct of_geom *geom, const int dir)
{
      int k ;

      if(dir == 1) {

	f[0] = 0 ;
	f[1] = state->Mcon[1][2] ;
	f[2] = state->Mcon[1][3] ;
	f[3] = state->T[1][1] ;
	f[4] = state->T[1][2] ;
	f[5] = state->T[1][3] ;

      }

      else if(dir == 2) {

	f[0] = state->Mcon[2][1] ;
	f[1] = 0 ;
	f[2] = state->Mcon[2][3] ;
	f[3] = state->T[2][1] ;
	f[4] = state->T[2][2] ;
	f[5] = state->T[2][3] ;

      }
      else fprintf(stderr,"Unimplemented dir value.  Check phys.c") ;

      /* multiply by geometric factor */
      PLOOP f[k] *= geom->g ;
}

/*
 * Convert conserved quantities back to primitives
 */

void UtoPrim(double pr[NPR], double U[NPR], const struct of_geom *geom)
{
      int k ;

      double F1t, F2t, F3t, Mt1, Mt2, Mt3 ;  /* cons */
      double Bsq ;

      double T0[NDIM] ;   /* stress tensor terms T^t_i */
      double Mx[NDIM] ;   /* mixed Maxwell components M^t_i */

      PLOOP  U[k] *= 1/(geom->g) ;

      Mt1 = U[0] ;  /* M^ti */
      Mt2 = U[1] ;
      Mt3 = U[2] ;
      T0[1] = U[3] ; /* T^t_i terms (t up, i down) */
      T0[2] = U[4] ;
      T0[3] = U[5] ;

      /* lower second index of Maxwell to get M^t_i */

      Mx[1] = geom->gcov[1][1] * Mt1
	    + geom->gcov[1][2] * Mt2
	    + geom->gcov[1][3] * Mt3 ;
      Mx[2] = geom->gcov[2][1] * Mt1
	    + geom->gcov[2][2] * Mt2
	    + geom->gcov[2][3] * Mt3 ;
      Mx[3] = geom->gcov[3][1] * Mt1
	    + geom->gcov[3][2] * Mt2
	    + geom->gcov[3][3] * Mt3 ;

       /* denominator */
      
      Bsq = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
      
      /* calculate prims: F^it
       *
       * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
       * cyclic 1,2,3.
       *
       */
      
      F1t = (1/geom->g)*(Mx[2]*T0[3] - Mx[3]*T0[2])/Bsq ;
      F2t = (1/geom->g)*(Mx[3]*T0[1] - Mx[1]*T0[3])/Bsq ;
      F3t = (1/geom->g)*(Mx[1]*T0[2] - Mx[2]*T0[1])/Bsq ;
      
      /* reassign prims */

      pr[F10] = F1t ; /* F^it */
      pr[F20] = F2t ;
      pr[F30] = F3t ;
      pr[M01] = Mt1 ; /* M^ti */
      pr[M02] = Mt2 ;
      pr[M03] = Mt3 ;
}      

/*
 * add in source terms (e.g. connection times stress tensor)
 */
/* ph unused here.  Check if this is correct */
void source(const double ph_l[NPR], double dU[NPR],
	    const int ii, const int jj, 
	    const struct of_state *state, const struct of_geom *geom)
{
      int j,k ;

      /* contract FFDE stress tensor with connection */
      PLOOP dU[k] = 0. ;
      DLOOP {
	     dU[3] += state->T[j][k]*conn[ii][jj][k][1][j] ;
	     dU[4] += state->T[j][k]*conn[ii][jj][k][2][j] ;
	     dU[5] += state->T[j][k]*conn[ii][jj][k][3][j] ;
      }

      /* multiply by geometric factor */
      PLOOP dU[k] *= geom->g ;
}

//======================================================================
// Business functions
//======================================================================

/*
 * Calculate ALL the components of the stress-energy tensor
 * (note specific to FFDE), first index up, second index down.
 */

void stresstensor(double Mcon[NDIM][NDIM],
		  double T[NDIM][NDIM], const struct of_geom *geom)
{
      int i,j,k ;

      double Mcov[NDIM][NDIM] ; /* covariant Maxwell */
      double Msq ;

      /* get covariant Maxwell from contravariant */
      lower_T(Mcon, Mcov, geom) ;

      /* out with the old */
      DLOOP  T[j][k] = 0. ;

      /* in with the new:
       *
       * a, b, c, d run 0 to 4:
       *
       * T^a_b = M^ac * M_bc - 1/4 KroneckerDelta[a,b] (M_cd M^cd)
       *
       */

      Msq = 0. ;
      DLOOP Msq += Mcov[j][k]*Mcon[j][k] ;

      DLOOP {
	for(i = 0; i < NDIM; i++) {
	  T[j][k] += Mcon[j][i]*Mcov[k][i] ;
	}
      }
      DLOOPA T[j][j] -= 0.25 * Msq ;
}

/*
 * load local geometry into structure geom
 */

void get_geometry(const int ii, const int jj, const int kk,
		  struct of_geom *geom)
{
	int j,k ;

	DLOOP geom->gcov[j][k] = gcov[ii][jj][kk][j][k] ;
	DLOOP geom->gcon[j][k] = gcon[ii][jj][kk][j][k] ;
	geom->g = gdet[ii][jj][kk] ;
	icurr = ii ;
	jcurr = jj ;
}

/*
 * load stress tensor and contravariant Maxwell in structure state
 */

void get_state(const double pr[NPR], 
	       struct of_state *state, const struct of_geom *geom)
{

      /* get contravariant Maxwell and stress tensor */
      Max_con(pr, state->Mcon, geom) ;
      stresstensor(state->Mcon, state->T, geom) ;
}

/*
 * calculate the contravariant Maxwell
 */
void Max_con(const double pr[NPR],
	     double Mcon[NDIM][NDIM], const struct of_geom *geom)
{
      int j ;

      double F1t, F2t, F3t, Mt1, Mt2, Mt3 ;    /* prims */
      double Econ[NDIM], Ecov[NDIM] ;          /* E four-vector */
      double beta[NDIM] ;                      /* shift  3-vector  */

      /* lapse */
      const double alpha = 1./sqrt(- geom->gcon[0][0]);  

      /* beta */
      SLOOPA beta[j] = geom->gcon[0][j]*alpha*alpha ;

      F1t = pr[F10] ; /* F^it */
      F2t = pr[F20] ;
      F3t = pr[F30] ;
      Mt1 = pr[M01] ; /* M^ti */
      Mt2 = pr[M02] ;
      Mt3 = pr[M03] ;

      /* create E vector */
      Econ[0] = 0. ;
      Econ[1] = - alpha * F1t ;
      Econ[2] = - alpha * F2t ;
      Econ[3] = - alpha * F3t ;

      lower_v(Econ, Ecov, geom) ;

      /* diagonal */
      Mcon[0][0] = Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0. ;

      /* out with the old */
      //DLOOP  Mcon[j][k] = 0. ;

      /* M^ij terms:
       *
       * M^ij = -beta^i M^tj + beta^j M^ti
       * + alpha * (1/sqrt(-g)) * Ecov[k]
       *
       * use primitive variables for M's
       *
       */

      Mcon[1][2] = (-beta[1] * Mt2 + beta[2] * Mt1
	         + alpha/(geom->g) * Ecov[3]) ;
      Mcon[2][3] = (-beta[2] * Mt3 + beta[3] * Mt2
                 + alpha/(geom->g) * Ecov[1]) ;
      Mcon[3][1] = (-beta[3] * Mt1 + beta[1] * Mt3
                 + alpha/(geom->g) * Ecov[2]) ;
	 
      /* copy remaining spacial terms */
   
      Mcon[2][1] = -Mcon[1][2] ;
      Mcon[3][2] = -Mcon[2][3] ;
      Mcon[1][3] = -Mcon[3][1] ;

      /* time terms - easy */
      Mcon[0][1] = Mt1 ;
      Mcon[0][2] = Mt2 ;
      Mcon[0][3] = Mt3 ;
      Mcon[1][0] = - Mt1 ;
      Mcon[2][0] = - Mt2 ;
      Mcon[3][0] = - Mt3 ;
}

/*
 * lower index of a four-vector
 */

void lower_v(const double vcon[NDIM], double vcov[NDIM],
	     const struct of_geom *geom)
{
	vcov[0] = geom->gcov[0][0]*vcon[0] 
		+ geom->gcov[0][1]*vcon[1] 
		+ geom->gcov[0][2]*vcon[2] 
		+ geom->gcov[0][3]*vcon[3] ;
	vcov[1] = geom->gcov[1][0]*vcon[0] 
		+ geom->gcov[1][1]*vcon[1] 
		+ geom->gcov[1][2]*vcon[2] 
		+ geom->gcov[1][3]*vcon[3] ;
	vcov[2] = geom->gcov[2][0]*vcon[0] 
		+ geom->gcov[2][1]*vcon[1] 
		+ geom->gcov[2][2]*vcon[2] 
		+ geom->gcov[2][3]*vcon[3] ;
	vcov[3] = geom->gcov[3][0]*vcon[0] 
		+ geom->gcov[3][1]*vcon[1] 
		+ geom->gcov[3][2]*vcon[2] 
		+ geom->gcov[3][3]*vcon[3] ;
}

/*
 * lower both indices of an anti-symmetric  tensor
 *  -- go from contravariant to covariant
 */

void lower_T(double Tcon[NDIM][NDIM],
	     double Tcov[NDIM][NDIM], const struct of_geom *geom)
{
        int j,k ;

	/* out with the old */
	DLOOP  Tcov[j][k] = 0. ;

	/* in with the new */
	DLOOP {
       	       Tcov[0][1] += (geom->gcov[0][j])*(geom->gcov[1][k])*Tcon[j][k] ;
       	       Tcov[0][2] += (geom->gcov[0][j])*(geom->gcov[2][k])*Tcon[j][k] ;
       	       Tcov[0][3] += (geom->gcov[0][j])*(geom->gcov[3][k])*Tcon[j][k] ;
	       Tcov[1][2] += (geom->gcov[1][j])*(geom->gcov[2][k])*Tcon[j][k] ;
	       Tcov[2][3] += (geom->gcov[2][j])*(geom->gcov[3][k])*Tcon[j][k] ;
	       Tcov[3][1] += (geom->gcov[3][j])*(geom->gcov[1][k])*Tcon[j][k] ;
	      }

	/*
	 * diagonal is already set to zero,
	 * just copy the rest
	 */

	Tcov[1][0] = - Tcov[0][1] ;
	Tcov[2][0] = - Tcov[0][2] ;
	Tcov[3][0] = - Tcov[0][3] ;
	Tcov[2][1] = - Tcov[1][2] ;
	Tcov[3][2] = - Tcov[2][3] ;
	Tcov[1][3] = - Tcov[3][1] ;
}
