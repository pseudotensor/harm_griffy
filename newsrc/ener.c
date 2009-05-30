//========================= GRFFDE =========================//

#include "decs.h"

/*
 * output radial energy and angular momentum flux data to "ener.out"
 *
 * could also do all kinds of other output with this . . .
 */

void ener(FILE *fp)
{
      int i,j ; 

      double ener_flux, angular_flux ;

      struct of_geom geom ;
      struct of_state state ;

      /* print out for each entry: t, E_flux, L_flux, at different radii */

      fprintf(fp,"%10.5g",t);

      /* sum to get E_tot */

      ener_flux = 0. ;
      angular_flux = 0. ;
      i = i_horizon ;
      JLOOP {
	get_geometry(i,j,CENT,&geom) ;
	get_state(p[i][j],&state,&geom) ;
	ener_flux    += (geom.g)*(state.T[1][0])*dx[2] ;
	angular_flux += (geom.g)*(state.T[1][3])*dx[2] ;
      }
      ener_flux = -ener_flux*2.*M_PI ;
      angular_flux = -angular_flux*2.*M_PI ;
      fprintf(fp," %10.5g %10.5g",ener_flux,angular_flux);

      ener_flux = 0. ;
      angular_flux = 0. ;
      i = IS+N1/4 ;
      JLOOP {
	get_geometry(i,j,CENT,&geom) ;
	get_state(p[i][j],&state,&geom) ;
	ener_flux    += (geom.g)*(state.T[1][0])*dx[2] ;
	angular_flux += (geom.g)*(state.T[1][3])*dx[2] ;
      }
      ener_flux = -ener_flux*2.*M_PI ;
      angular_flux = -angular_flux*2.*M_PI ;
      fprintf(fp," %10.5g %10.5g",ener_flux,angular_flux);

      ener_flux = 0. ;
      angular_flux = 0. ;
      i = IS+N1/2;
      JLOOP {
	get_geometry(i,j,CENT,&geom) ;
	get_state(p[i][j],&state,&geom) ;
	ener_flux    += (geom.g)*(state.T[1][0])*dx[2] ;
	angular_flux += (geom.g)*(state.T[1][3])*dx[2] ;
      }
      ener_flux = -ener_flux*2.*M_PI ;
      angular_flux = -angular_flux*2.*M_PI ;
      fprintf(fp," %10.5g %10.5g",ener_flux,angular_flux);

      ener_flux = 0. ;
      angular_flux = 0. ;
      i = IE-2*NGHOSTS;
      JLOOP {
	get_geometry(i,j,CENT,&geom) ;
	get_state(p[i][j],&state,&geom) ;
	ener_flux    += (geom.g)*(state.T[1][0])*dx[2] ;
	angular_flux += (geom.g)*(state.T[1][3])*dx[2] ;
      }
      ener_flux = -ener_flux*2.*M_PI ;
      angular_flux = -angular_flux*2.*M_PI ;
      fprintf(fp," %10.5g %10.5g",ener_flux,angular_flux);

      fprintf(fp,"\n");

      return ;
}
