//============================== GRFFDE ==============================//

/*
 * sets the timestep
 *
 */

#include "decs.h"

double timestep()
{
  double vmax, vlarger ;
  double dt1, dt2 ;
  int i, j ;
  int imax, jmax ;


  /* Calculate worst dx[1] case */
  imax = -10 ;
  jmax = -10 ;
  vmax = -100.0 ; 
  ZLOOP {
    vlarger=vchar(i,j,1);
    if (vmax < vlarger){
      vmax=vlarger ;
      imax = i ;
      jmax = j ;
    }
  }
  dt1 = cour * dx[1] / vmax ;
  fprintf(stderr,
	  "1: vmax=%10.5g at i-IS=%d j-JS=%d\n",
	  vmax,imax-IS,jmax-JS) ;
  fprintf(stderr,"dt1= %10.5g\n",dt1);

  /* Calculate worst dx[2] case */
  imax = -10 ;
  jmax = -10 ;
  vmax = -100.0 ; 
  ZLOOP {
    vlarger=vchar(i,j,2);
    if (vmax < vlarger){
      vmax=vlarger ;
      imax = i ;
      jmax = j ;
    }
  }
  dt2 = cour * dx[2] / vmax ;
  fprintf(stderr,"2: vmax=%10.5g at i-IS=%d j-JS=%d\n",
	  vmax,imax-IS,jmax-JS) ;
  fprintf(stderr,"dt2= %10.5g\n",dt2);

  return 1. / (1. /dt1 + 1. / dt2) ;
}

double max2(double xx, double yy)
{
  return (xx > yy) ? xx : yy ;
}

double min2(double xx, double yy)
{
  return (xx < yy) ? xx : yy ;
}

/* Calculates the largest phase velocity in a given
     direction, assuming MKS coordinates in the
     calculation for dir==2 */
double vchar(const int i, const int j, const int dir)
{
  double vplus, vminus, vlarger;
  double coeff_a, coeff_b, coeff_c;
 
  if (dir==1){
    coeff_c = gcon[i][j][CENT][1][1] ;
    coeff_b = -2. * gcon[i][j][CENT][1][0] ;
    coeff_a = gcon[i][j][CENT][0][0] ;
  
    vplus = -coeff_b + sqrt(coeff_b*coeff_b - 4. * coeff_a*coeff_c) ;
    vplus = fabs (vplus / (2. * coeff_a) ) ;

    vminus = -coeff_b - sqrt(coeff_b*coeff_b - 4. * coeff_a*coeff_c) ;
    vminus = fabs (vminus / (2. * coeff_a) ) ;

    vlarger = max2(vplus, vminus) ;
  }

  if (dir==2){ 
    coeff_c = gcon[i][j][CENT][2][2] ;
    coeff_a = gcon[i][j][CENT][0][0] ;

    vlarger=sqrt(fabs(coeff_c/coeff_a)) ;
  }

  return vlarger ;
}
