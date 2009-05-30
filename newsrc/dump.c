//============================== GRFFDE ==============================//

/*
 * for gr trials
 */


#include "decs.h"

void dump(FILE *fp)
{
	int i,j,k ;
	double divb ;
	double X[NDIM] ;
	double r,th ;
	double BsqmEsqvar;

	struct of_geom geom ;
	struct of_state state ;

	/* print on first line */	
	fprintf(fp,"%10.5g %d %d %d %10.5g %10.5g %10.5g %10.5g\n",
			t,N1,N2,i_horizon-IS,startx[1],startx[2],dx[1],dx[2]) ;

	/*
	 * colums: i  j  x[1]  x[2]  r  th  b1  b2  b3  e1
	 *         e2  e3  divb  gdet  E_flux
	 */

	ZLOOP {
		coord(i,j,CENT,X) ;
		bl_coord(X,&r,&th) ;

		get_geometry(i,j,CENT,&geom) ;

		fprintf(fp,"%d    %d   ",i-IS,j-JS) ;  
		fprintf(fp,"%15.7g %15.7g",X[1],X[2]) ;
		fprintf(fp,"%15.7g %15.7g",r,th) ;
		PLOOP fprintf(fp,"%15.7g ",(p[i][j][k])) ;


                /* divb flux-ct defn; corner-centered.  Use
		   only interior corners */
		if(i > IS && j > JS && i <= IE && j <= JE ) {
			divb = fabs( 0.5*(
				p[i][j][M01]*gdet[i][j][CENT]
				+ p[i][j-1][M01]*gdet[i][j-1][CENT]
				- p[i-1][j][M01]*gdet[i-1][j][CENT]
				- p[i-1][j-1][M01]*gdet[i-1][j-1][CENT]
				)/dx[1] +
				0.5*(
				p[i][j][M02]*gdet[i][j][CENT]
				+ p[i-1][j][M02]*gdet[i-1][j][CENT]
				- p[i][j-1][M02]*gdet[i][j-1][CENT]
				- p[i-1][j-1][M02]*gdet[i-1][j-1][CENT]
				)/dx[2]) ;
		}

		else divb = 0. ;

		fprintf(fp,"%15.7g ",divb) ;
		fprintf(fp,"%15.7g",geom.g) ;

		/*
                for(int jdim=0;jdim<NDIM;jdim++)
		for(int kdim=0;kdim<NDIM;kdim++)
		  fprintf(fp," M %d %d %15.7g %15.7g %15.7g %15.7g",jdim,kdim,
			  geom.gcov[jdim][kdim], gcov[i][j][CENT][jdim][kdim],
			  geom.gcon[jdim][kdim], gcon[i][j][CENT][jdim][kdim]);
		*/

		/* for energy flux */
		get_state(p[i][j],&state,&geom) ;

		fprintf(fp,"%15.7g",-state.T[1][0]) ;

		BsqmEsq(p[i][j],&geom,&BsqmEsqvar);

		fprintf(fp,"%15.7g",BsqmEsqvar) ;


		/* since get_state has already been called, you 
		 * could also print out all of T or Mcon
		 */

		fprintf(fp,"\n") ;

	}

}
