//============================ GRFFDE =============================//

#include "decs.h"

#define KS

#define SMALL   1.e-20 

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions.
 *
 **/

/* should return Boyer-Lindquist coordinte of point */
void bl_coord(const double X[NDIM], double *r, double *th)
{
        *r = exp(X[1]) + R0 ;
        *th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]) ;
}

/* insert metric here: */
void gcov_func(const double X[NDIM], double gcov_local[NDIM][NDIM])
{
      int j,k ;

      /*
       * chose metric/coordinates: Minkowski or Kerr-Schild
       */

#ifdef MINK
      DLOOP {
           if(j==k) {
                 if(j == 0)
                       gcov_local[j][k] = -1. ;
                 else
                       gcov_local[j][k] = 1. ;
           }
           else
                 gcov_local[j][k] = 0. ;
      }
#endif /* MINK */

#ifdef KS
        double sth,cth,s2,rho2 ;
        double r,th ;
        double tfac,rfac,hfac,pfac ;

        DLOOP gcov_local[j][k] = 0. ;

        bl_coord(X,&r,&th) ;

        cth = cos(th) ;
        sth = fabs(sin(th)) ;
        if (sth<SMALL) sth=SMALL ;
        s2 = sth*sth ;
        rho2 = r*r + a*a*cth*cth ;

        tfac = 1. ;
        rfac = r - R0 ;
        hfac = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) ;
        pfac = 1. ;

        gcov_local[TT][TT] = (-1. + 2.*r/rho2)       * tfac*tfac ;
        gcov_local[TT][RR] = (2.*r/rho2)             * tfac*rfac ;
        gcov_local[TT][TH] = 0.                      * tfac*hfac ;
        gcov_local[TT][PH] = (-2.*a*r*s2/rho2)       * tfac*pfac ;

        gcov_local[RR][TT] = gcov_local[TT][RR] ;
        gcov_local[RR][RR] = (1. + 2.*r/rho2)         * rfac*rfac ;
        gcov_local[RR][TH] = 0.                       * rfac*hfac ;
        gcov_local[RR][PH] = (-a*s2*(1. + 2.*r/rho2)) * rfac*pfac ;

        gcov_local[TH][TT] = gcov_local[TT][TH] ;
        gcov_local[TH][RR] = gcov_local[RR][TH] ;
        gcov_local[TH][TH] = rho2                     * hfac*hfac ;
        gcov_local[TH][PH] = 0.                       * hfac*pfac ;

        gcov_local[PH][TT] = gcov_local[TT][PH] ;
        gcov_local[PH][RR] = gcov_local[RR][PH] ;
        gcov_local[PH][TH] = gcov_local[TH][PH] ;
        gcov_local[PH][PH] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)) * pfac*pfac ;

#endif /* KS */
}


void fix_flux(double F1_local[INA][JNA][NPR], double F2_local[INA][JNA][NPR])
{
  /* Intended to set F2 to zero at the axis */

  int i,k ;
  for(i=IS-1;i<=IE+1;i++) {

    /* Designed to get F2[i][JS][M01]=F2[i][JE+1][M01]=0. after
       constrained transport is done with Toth's averaging approach */
    F1_local[i][JS-1][M02] = -F1_local[i][JS][M02] ;
    F1_local[i][JE+1][M02] = -F1_local[i][JE][M02] ;

    PLOOP F2_local[i][JE+1][k] = 0. ;
    PLOOP F2_local[i][JS][k]   = 0. ;
  }
}

#ifndef NUMERICAL
/* 
 * returns sqrt(-g) for the Kerr metric in the
 * modified Kerr Schild coordinates.
 * Does not assume if gcov has been set first or not.
 */
double gdet_func(const double X[NDIM], double gcov_local[NDIM][NDIM]) 
{
#ifdef MINK
  return 1.0 ;
#endif /* MINK */

#ifdef KS
        double sth,cth,s2,rho2 ;
        double r,th ;
        double tfac,rfac,hfac,pfac ;

        bl_coord(X,&r,&th) ;

        cth = cos(th) ;
        sth = fabs(sin(th)) ;
        if (sth<SMALL) sth=SMALL ;
        s2 = sth*sth ;
        rho2 = r*r + a*a*cth*cth ;

        tfac = 1. ;
        rfac = r - R0 ;
        hfac = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) ;
        pfac = 1. ;

        return rho2*sth*tfac*rfac*hfac*pfac;
#endif /* KS */
}

void gcon_func(const double X[NDIM], double gcov_local[NDIM][NDIM], 
               double gcon_local[NDIM][NDIM])
{
      int j,k ;
#ifdef MINK
      DLOOP {
           if(j==k) {
                 if(j == 0)
                       gcon_local[j][k] = -1. ;
                 else
                       gcon_local[j][k] = 1. ;
           }
           else
                 gcon_local[j][k] = 0. ;
      }
#endif /* MINK */

#ifdef KS
        double sth,cth,s2,rho2, rho2i;
        double r,th ;
        double tfaci,rfaci,hfaci,pfaci ;

        DLOOP gcon_local[j][k] = 0. ;

        bl_coord(X,&r,&th) ;

        cth = cos(th) ;
        sth = fabs(sin(th)) ;
        if (sth<SMALL) sth=SMALL ;
        s2 = sth*sth ;
        rho2 = r*r + a*a*cth*cth ;
        rho2i = 1./rho2 ;

        tfaci = 1. ;
        rfaci = 1./(r - R0) ;
        hfaci = 1./(M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2])) ;
        pfaci = 1. ;

        gcon_local[TT][TT] = -(r*(r+2.)+a*a*cth*cth) * rho2i * tfaci*tfaci ;
        gcon_local[TT][RR] = 2.*r                    * rho2i * tfaci*rfaci ;
        gcon_local[TT][TH] = 0.                      * rho2i * tfaci*hfaci ;
        gcon_local[TT][PH] = 0.                      * rho2i * tfaci*pfaci ;

        gcon_local[RR][TT] = gcon_local[TT][RR] ;
        gcon_local[RR][RR] = (r*(r-2.)+a*a)          * rho2i * rfaci*rfaci ;
        gcon_local[RR][TH] = 0.                      * rho2i * rfaci*hfaci ;
        gcon_local[RR][PH] = a                       * rho2i * rfaci*pfaci ;

        gcon_local[TH][TT] = gcon_local[TT][TH] ;
        gcon_local[TH][RR] = gcon_local[RR][TH] ;
        gcon_local[TH][TH] = 1.                      * rho2i * hfaci*hfaci ;
        gcon_local[TH][PH] = 0.                      * rho2i * hfaci*pfaci ;

        gcon_local[PH][TT] = gcon_local[TT][PH] ;
        gcon_local[PH][RR] = gcon_local[RR][PH] ;
        gcon_local[PH][TH] = gcon_local[TH][PH] ;
        gcon_local[PH][PH] = (1./s2)                 * rho2i * pfaci*pfaci ;
#endif /* KS */
}

void conn_func(const double X[NDIM], const struct of_geom *geom,
               double gamma[NDIM][NDIM][NDIM])
{
#ifdef MINK
  int i,j,k ;
  for(i=0;i<NDIM;i++)
    for(j=0;j<NDIM;j++) 
      for(k=0;k<NDIM;k++)
        gamma[i][j][k] = 0. ;
#endif /* MINK */

#ifdef KS

#define SQ(x)     ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define FOURTH(x) ((x)*(x)*(x)*(x))
#define FIFTH(x)  ((x)*(x)*(x)*(x)*(x))
#define SIXTH(x)  ((x)*(x)*(x)*(x)*(x)*(x))

  double sth,cth,s2,c2;
  double r,th ;
  double tfac,rfac,hfac,pfac ;
  double tfaci,rfaci,hfaci,pfaci ;
  double Pi, Delta ;
  double Sigma, Sigmai, Sigmai2, Sigmai3 ;
  double longcoefficient, longcoefficient2, longcoefficient3 ;
  double longcoefficient4, longcoefficient5, longcoefficient6 ;
  double a2,a3,a4,a5;
  double r2,r3,r4,r5;
  double rfacprime, hfacprime;

  bl_coord(X,&r,&th) ;

  a2=SQ(a);
  a3=CUBE(a);
  a4=FOURTH(a);
  a5=FIFTH(a);

  r2=SQ(r);
  r3=CUBE(r);
  r4=FOURTH(r);
  r5=FIFTH(r);

  cth = cos(th) ;
  sth = fabs(sin(th)) ;
  if (sth<SMALL) sth=SMALL ;

  s2 = SQ(sth) ;
  c2 = SQ(cth) ;
  Sigma = r2 + a2 * c2 ;
  Sigmai = 1./ Sigma ;
  Sigmai2 = SQ(Sigmai) ;
  Sigmai3 = CUBE(Sigmai) ;
  Pi      =  r2 - a2 * c2 ;
  Delta   = r2 +a2 -2.*r ;

  longcoefficient  = (r5+r*a4*FOURTH(cth)-a2*r2*s2+c2*(2.*a2*r3+a4*s2)) ;
  longcoefficient2 = (a4*c2*(2.*c2-1.)    +a2*r2*s2-2.*r*Pi-r4) ;
  longcoefficient4 = (a2*c2*(2*r2+a2*c2) +2.*a2*r*s2+2.*r*Sigma+r4) ;
  longcoefficient3 = (r3+r4-a2*r*c2-FOURTH(a*cth)) ;
  longcoefficient5 = SIXTH(r)+SIXTH(a*cth) + a2*r3*(4.+r)*s2 +
    2.*a4*r*FOURTH(sth)+
    FOURTH(cth)*(3.*a4*r2+SIXTH(a)*s2)+
    a2*r*c2*(3.*r3+2.*a2*(r+2.)*s2) ;
  longcoefficient6 = SQ(Sigma)/sth + a2*r*2.*sth ;

  tfac = 1. ;
  rfac = r - R0 ;
  hfac = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) ;
  pfac = 1. ;

  tfaci = 1. ;
  rfaci = 1./ rfac ;
  hfaci = 1./ hfac ;
  pfaci = 1. ;

  rfacprime = 1. ;
  hfacprime = (2. * SQ(M_PI) * (hslope - 1.) * sin(2.*M_PI*X[2]) )/ hfac ;

  gamma[TT][TT][TT] = 2.*r*Pi*Sigmai3                                  *tfac ;
  gamma[TT][TT][RR] = Pi*(Sigma+2.*r)*Sigmai3                          *rfac ;
  gamma[TT][TT][TH] = -2.*a2*r*cth*sth*Sigmai2                         *hfac ;
  gamma[TT][TT][PH] = -2.*a*r*Pi*s2*Sigmai3                            *pfac ;

  gamma[TT][RR][TT] = gamma[TT][TT][RR] ;
  gamma[TT][RR][RR] = 2.* longcoefficient3*Sigmai3          *tfaci*rfac*rfac ;
  gamma[TT][RR][TH] = -2.*a2*r*cth*sth*Sigmai2              *tfaci*rfac*hfac ;
  gamma[TT][RR][PH] = -a*Pi*(Sigma+2.*r)*s2*Sigmai3         *tfaci*rfac*pfac ;

  gamma[TT][TH][TT] = gamma[TT][TT][TH] ;
  gamma[TT][TH][RR] = gamma[TT][RR][TH] ;
  gamma[TT][TH][TH] = -2.*r2*Sigmai                         *tfaci*hfac*hfac ;
  gamma[TT][TH][PH] = 2.*CUBE(a*sth)*r*cth*Sigmai2          *tfaci*hfac*pfac ;

  gamma[TT][PH][TT] = gamma[TT][TT][PH] ;
  gamma[TT][PH][RR] = gamma[TT][RR][PH] ;
  gamma[TT][PH][TH] = gamma[TT][TH][PH] ;
  gamma[TT][PH][PH] = -2.*r*s2*longcoefficient*Sigmai3      *tfaci*pfac*pfac ; 


  gamma[RR][TT][TT] = Delta*Pi*Sigmai3                      *rfaci*tfac*tfac ;
  gamma[RR][TT][RR] = -Pi*(2.*r-a2*s2)*Sigmai3                    *tfac      ;
  gamma[RR][TT][TH] = 0. ;
  gamma[RR][TT][PH] = -a* Delta*Pi*s2 *Sigmai3              *rfaci*tfac*pfac ;

  gamma[RR][RR][TT] = gamma[RR][TT][RR] ;
  gamma[RR][RR][RR] = longcoefficient2 * Sigmai3*rfac + rfacprime ;
  gamma[RR][RR][TH] = -a2*sth*cth*Sigmai                               *hfac ;
  gamma[RR][RR][PH] = a*s2*(longcoefficient+2.*r*Pi)*Sigmai3           *pfac ;

  gamma[RR][TH][TT] = gamma[RR][TT][TH] ;
  gamma[RR][TH][RR] = gamma[RR][RR][TH] ;
  gamma[RR][TH][TH] = -r*Delta/Sigma                        *rfaci*hfac*hfac ;
  gamma[RR][TH][PH] = 0. ;

  gamma[RR][PH][TT] = gamma[RR][TT][PH] ;
  gamma[RR][PH][RR] = gamma[RR][RR][PH] ;
  gamma[RR][PH][TH] = gamma[RR][TH][PH] ;
  gamma[RR][PH][PH] = -longcoefficient*Delta*s2*Sigmai3     *rfaci*pfac*pfac ;
  

  gamma[TH][TT][TT] = -2.*a2*r*cth*sth*Sigmai3              *hfaci*tfac*tfac ;
  gamma[TH][TT][RR] = -2.*a2*r*cth*sth*Sigmai3              *hfaci*tfac*rfac ; 
  gamma[TH][TT][TH] = 0. ;
  gamma[TH][TT][PH] =  2.*a*r*(a2+r2)*cth*sth*Sigmai3       *hfaci*tfac*pfac ; 

  gamma[TH][RR][TT] = gamma[TH][TT][RR] ;
  gamma[TH][RR][RR] = -2.*a2*r*cth*sth*Sigmai3              *hfaci*rfac*rfac ;
  gamma[TH][RR][TH] = r*Sigmai                                    *rfac      ;
  gamma[TH][RR][PH] = a*cth*sth*longcoefficient4*Sigmai3    *hfaci*rfac*pfac ;

  gamma[TH][TH][TT] = gamma[TH][TT][TH] ;
  gamma[TH][TH][RR] = gamma[TH][RR][TH] ;
  gamma[TH][TH][TH] = -a2*cth*sth*Sigmai*hfac + hfacprime      ;
  gamma[TH][TH][PH] = 0. ;

  gamma[TH][PH][TT] = gamma[TH][TT][PH] ;
  gamma[TH][PH][RR] = gamma[TH][RR][PH] ;
  gamma[TH][PH][TH] = gamma[TH][TH][PH] ;
  gamma[TH][PH][PH] = -cth*sth*longcoefficient5*Sigmai3     *hfaci*pfac*pfac ;


  gamma[PH][TT][TT] = a*Pi*Sigmai3                          *pfaci*tfac*tfac ;
  gamma[PH][TT][RR] = a*Pi*Sigmai3                          *pfaci*tfac*rfac ; 
  gamma[PH][TT][TH] = -2.*a*r*cth*Sigmai2/sth               *pfaci*tfac*hfac ;
  gamma[PH][TT][PH] = -a2*Pi*s2*Sigmai3                           *tfac      ;
                                         
  gamma[PH][RR][TT] = gamma[PH][TT][RR] ;
  gamma[PH][RR][RR] = a*Pi*Sigmai3                          *pfaci*rfac*rfac ;
  gamma[PH][RR][TH] = -a*cth*(Sigma+2.*r)*Sigmai2/sth       *pfaci*rfac*hfac ;
  gamma[PH][RR][PH] = longcoefficient*Sigmai3                     *rfac      ;
                                         
  gamma[PH][TH][TT] = gamma[PH][TT][TH] ;
  gamma[PH][TH][RR] = gamma[PH][RR][TH] ;
  gamma[PH][TH][TH] = -a*r*Sigmai                           *pfaci*hfac*hfac ;
  gamma[PH][TH][PH] = longcoefficient6 *cth*Sigmai2               *hfac      ;
                                         
  gamma[PH][PH][TT] = gamma[PH][TT][PH] ;
  gamma[PH][PH][RR] = gamma[PH][RR][PH] ;
  gamma[PH][PH][TH] = gamma[PH][TH][PH] ;
  gamma[PH][PH][PH] = -longcoefficient*a*s2*Sigmai3               *pfac      ;

#undef SQ
#undef CUBE
#undef FOURTH
#undef FIFTH
#undef SIXTH

#endif /* KS */
}
#endif /* NUMERICAL */

/* some grid location, dxs */
void set_points()   /* flat space! and rectangular Cartesian grid! */
{
#ifdef MINK
         startx[1] = -0.5 ;
         startx[2] = -0.5 ;
         dx[1] = 1./N1 ;
         dx[2] = 1./N2 ;
#endif /* MINK */
#ifdef KS 
        startx[1] = log(Rin - R0) ;
        startx[2] = 0. ;
        dx[1] = log((Rout - R0)/(Rin - R0))/N1 ;
        dx[2] = 1./N2 ;
#endif /* KS */
}
#undef SMALL
#undef KS
