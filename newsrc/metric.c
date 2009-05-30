//=========================== GRFFDE ==============================//

#include "decs.h"

void coord(const int i, const int j, const int loc, double X[NDIM])
{
  int iloc=i-IS, jloc=j-JS;
  if(loc == FACE1) {
    X[1] = startx[1] + iloc*dx[1] ;
    X[2] = startx[2] + (jloc + 0.5)*dx[2] ;
  }
  else if(loc == FACE2) {
    X[1] = startx[1] + (iloc + 0.5)*dx[1] ;
    X[2] = startx[2] + jloc*dx[2] ;
  }
  else if(loc == CENT) {
    X[1] = startx[1] + (iloc + 0.5)*dx[1] ;
    X[2] = startx[2] + (jloc + 0.5)*dx[2] ;
  }
  else if(loc == CORN) {
    X[1] = startx[1] + iloc*dx[1] ;
    X[2] = startx[2] + jloc*dx[2] ;
  }
  else {
    fprintf(stderr,"Unknown grid location\n");
    exit(11);
  }
}


#ifdef NUMERICAL
/* 
 * returns sqrt(-g),
 * assumes gcov has been set first
 */

double gdet_func(const double X[NDIM], double gcov_local[NDIM][NDIM]) 
{
	double tmp[NDIM][NDIM] ;
	double d ;
	int j,k,indx[NDIM] ;

	DLOOP tmp[j][k] = gcov_local[j][k] ;

	ludcmp(NDIM,tmp,indx,&d) ;

	for(j=0;j<NDIM;j++) d *= tmp[j][j] ;
	return(sqrt(fabs(d))) ; 
}

/*
 * invert gcov to get gcon
 */

void gcon_func(const double X[NDIM], 
	       double gcov_local[NDIM][NDIM], 
	       double gcon_local[NDIM][NDIM])
{
	int j,k ;
	double tmp[NDIM][NDIM] ;

        DLOOP tmp[j][k] = gcov_local[j][k] ;

        gaussj_invert(NDIM,tmp) ;

        DLOOP gcon_local[j][k] = tmp[k][j] ;
}

/* 
 * this gives the connection coefficient
 *	\Gamma^{i}_{j,k} = conn[..][i][j][k]
 * where i = {1,2,3,4} corresponds to {t,r,theta,phi}
 */

#define DELTA 1.e-5

/* NOTE: parameter hides global variable */
void conn_func(const double X[NDIM],
	       const struct of_geom *geom,
	       double conn_local[NDIM][NDIM][NDIM])
{
	int i,j,k,l ;
	double tmp[NDIM][NDIM][NDIM] ;
	double Xh[NDIM],Xl[NDIM] ;
	double gh[NDIM][NDIM] ;
	double gl[NDIM][NDIM] ;

	for(k=0;k<NDIM;k++) {
		for(l=0;l<NDIM;l++) Xh[l] = X[l] ;
		for(l=0;l<NDIM;l++) Xl[l] = X[l] ;
		Xh[k] += DELTA ;
		Xl[k] -= DELTA ;
		gcov_func(Xh,gh) ;
		gcov_func(Xl,gl) ;

		for(i=0;i<NDIM;i++)
		for(j=0;j<NDIM;j++) 
		  conn_local[i][j][k] = (gh[i][j] - gl[i][j])/(Xh[k] - Xl[k]) ;
	}

	/* now rearrange to find \Gamma_{ijk} */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++) 
	  tmp[i][j][k] = 0.5*(conn_local[j][i][k] + 
			      conn_local[k][i][j] - conn_local[k][j][i]) ;

	/* finally, raise index */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++)  {
	  conn_local[i][j][k] = 0. ;
	  for(l=0;l<NDIM;l++) 
	    conn_local[i][j][k] += geom->gcon[i][l]*tmp[l][j][k] ;
	}

	/* done! */
}

#undef DELTA

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void gaussj_invert(const int n,double array[n][n])
{
  int indxc[n], indxr[n], ipiv[n];
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;
  int lb = 1 ;

  for (j=1;j<=n;j++) ipiv[j-lb]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j-lb] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k-lb] == 0) {
	    if (fabs(array[j-lb][k-lb]) >= big) {
	      big=fabs(array[j-lb][k-lb]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k-lb] > 1) {
	    fprintf(stderr,
		    "choke in gaussj zone= %d %d\n",
		    icurr, jcurr) ;
	    exit(2) ;
	  }
	}
    ++(ipiv[icol-lb]);
    if (irow != icol) {
      for (l=1;l<=n;l++) SWAP(array[irow-lb][l-lb],array[icol-lb][l-lb]) ;
    }
    indxr[i-lb]=irow;
    indxc[i-lb]=icol;
    if (array[icol-lb][icol-lb] == 0.0) {
      fprintf(stderr,"gaussj: Singular Matrix-2 zone = %d %d\n",
	      icurr,jcurr);
      exit(3) ;
    }
    pivinv=1.0/array[icol-lb][icol-lb];
    array[icol-lb][icol-lb]=1.0;
    for (l=1;l<=n;l++) array[icol-lb][l-lb] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=array[ll-lb][icol-lb];
	array[ll-lb][icol-lb]=0.0;
	for (l=1;l<=n;l++) array[ll-lb][l-lb] -= array[icol-lb][l-lb]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l-lb] != indxc[l-lb])
      for (k=1;k<=n;k++)
	SWAP(array[k-lb][indxr[l-lb]-lb],array[k-lb][indxc[l-lb]-lb]);
  }
}
#undef SWAP

#define TINY 1.0e-20;
void ludcmp(const int n,double array[n][n],int indx[n],double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double vv[n];
	int lb=1;

	imax = 0 ;

	*d=1.0;
	for (i=1;i<=n;i++) {
	  big=0.0;
	  for (j=1;j<=n;j++) {
	    if ((temp=fabs(array[i-lb][j-lb])) > big) big=temp;
	  }
	  if (big == 0.0) {
	    fprintf(stderr,"sing. matr. in ludcmp: %d %d\n",icurr,jcurr) ;
	  }
	  vv[i-lb]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
		  sum=array[i-lb][j-lb];
		  for (k=1;k<i;k++) sum -= array[i-lb][k-lb]*array[k-lb][j-lb];
		  array[i-lb][j-lb]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
		  sum=array[i-lb][j-lb];
		  for (k=1;k<j;k++)
		    sum -= array[i-lb][k-lb]*array[k-lb][j-lb];
		  array[i-lb][j-lb]=sum;
		  if ( (dum=vv[i-lb]*fabs(sum)) >= big) {
		    big=dum;
		    imax=i;
		  }
		}
		if (j != imax) {
		  for (k=1;k<=n;k++) {
		    dum=array[imax-lb][k-lb];
		    array[imax-lb][k-lb]=array[j-lb][k-lb];
		    array[j-lb][k-lb]=dum;
		  }
		  *d = -(*d);
		  vv[imax-lb]=vv[j-lb];
		}
		indx[j-lb]=imax;
		if (array[j-lb][j-lb] == 0.0) array[j-lb][j-lb]=TINY;
		if (j != n) {
		  dum=1.0/(array[j-lb][j-lb]);
		  for (i=j+1;i<=n;i++) array[i-lb][j-lb] *= dum;
		}
	}
}
#undef TINY

#endif /* NUMERICAL */
