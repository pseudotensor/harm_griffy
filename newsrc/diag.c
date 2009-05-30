//============================ GRFFDE ===============================//

#include "decs.h"

/* all diagnostics subroutine */

void diag(const int call_code)
{
#define FILENAMESIZE 100
  char dfnam[FILENAMESIZE],
    binnam[FILENAMESIZE];
#undef FILENAMESIZE

  int i,j,k;
  int imax, jmax ;
  int i2max, j2max ;
  int ii, jj;
  double divb,divbmax ;            /* for divB */
  double BsqmEsqvar,BsqmEsqmin;

  static int terminationflag=0;
  static int alloutputflag=0;
  static int count=0;

  static FILE *ener_file ;
  FILE *dump_file ;

  struct of_geom geom ;
  struct of_state state ;
  double temp[INA][JNA] ;
  double Mcov[NDIM][NDIM] ;        /* covariant Maxwell */

  double U[NPR] ;                  /* cons */
  double cons[NPR] ;
    
  static double cons_init[NPR] ;   /* initial cons */
  static double cons_fin[NPR] ;    /* final cons */

  static double dV ;
  dV = dx[1]*dx[2] ;

  if(call_code==INIT_OUT) {
    /*  Initialize counters */
    dump_cnt = 0 ;
    bin_cnt = 0 ;
    /* create ener.out file */
    ener_file = fopen("ener.out","a") ;
    if(ener_file==NULL) {
      fprintf(stderr,"error opening energy output file\n") ;
      exit(1) ;
    }
    fprintf(stderr," %d\n",NGHOSTS);
    fprintf(stderr," %d",IS      );
    fprintf(stderr," %d",IE      );
    fprintf(stderr," %d",ISA     );
    fprintf(stderr," %d",IEA     );
    fprintf(stderr," %d",ISG     );
    fprintf(stderr," %d",IEG     );
    fprintf(stderr," %d",IN      );
    fprintf(stderr," %d",INA     );
    fprintf(stderr," %d\n",ING     );
    fprintf(stderr," %d",JS      );
    fprintf(stderr," %d",JE      );
    fprintf(stderr," %d",JSA     );
    fprintf(stderr," %d",JEA     );
    fprintf(stderr," %d",JSG     );
    fprintf(stderr," %d",JEG     );
    fprintf(stderr," %d",JN      );
    fprintf(stderr," %d",JNA     );
    fprintf(stderr," %d\n",JNG     );
#ifdef GRIDOUTPUT
    gridoutput(binnam);
#endif /* GRIDOUTPUT */
    }


  /* calculate conserved quantities */
  if(call_code==INIT_OUT || 
     call_code==LOG_OUT ||
     call_code==FINAL_OUT ) {
    /* zero everything */
    PLOOP  cons[k] = 0. ;
    divbmax = 0. ;
    BsqmEsqmin = 1E30;

    ZLOOP  {
      get_geometry(i,j,CENT,&geom) ;
      primtoU(p[i][j],U,&geom) ;
      /* update cons */
      PLOOP cons[k] += U[k]*dV ;

      /* flux-ct defn */
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

      BsqmEsq(p[i][j],&geom,&BsqmEsqvar);


      /* calculate divbmax */
      if(divb > divbmax && i > IS && j > JS) {
	imax = i ;
	jmax = j ;
	divbmax = divb ;
      }

      /* calculate BsqmEsqmin */
      if(BsqmEsqvar < BsqmEsqmin && i > IS && j > JS) {
	i2max = i ;
	j2max = j ;
	BsqmEsqmin = BsqmEsqvar ;
      }

    }
    if (alloutputflag==0 && (fabs(divbmax) > 1.e-9 )){
      alloutputflag=1;
    }
    if (fabs(divbmax)>0.001){
      terminationflag=1;
    }
  }

  if (alloutputflag==1) count++;
  if (count>50) terminationflag=1;

  /* record initial values of cons */
  if(call_code == INIT_OUT) PLOOP cons_init[k] = cons[k] ;

  /* final aout */
  if(call_code == FINAL_OUT || alloutputflag==1) {

    /* record final values of cons */
    PLOOP cons_fin[k] = cons[k] ;
 
    /* print conserved info in terminal */
    fprintf(stderr,"\n\nCon_quantities: ini,fin,del:\n\n") ;
    PLOOP fprintf(stderr,"%d     %10.5g     %10.5g     %10.5g\n",
		  k,cons_init[k],cons_fin[k],(cons_fin[k]-cons_init[k])) ;
  }

  /* log entry */
  if(call_code == INIT_OUT || 
     call_code == LOG_OUT ||
     call_code == FINAL_OUT || alloutputflag==1) {

    ener(ener_file) ;
    fprintf(stderr,"\ndivbmax: %g  %d  %d\n\n",divbmax,imax-IS,jmax-JS) ;
    fprintf(stderr,"\nBsqmEsqmin: %g  %d  %d\n\n",BsqmEsqmin,i2max-IS,j2max-JS) ;
    //fprintf(ener_file,"%10.5g  %10.5g",t,) ;
    //fprintf(ener_file,"\n") ;
    fflush(ener_file) ;
  }


  /* dump at regular intervals */
  if(call_code == INIT_OUT || 
     call_code == DUMP_OUT ||
     call_code == FINAL_OUT  || alloutputflag==1 ) {
    /* make regular dump file */
    sprintf(dfnam,"dump%03d",dump_cnt) ;
    fprintf(stderr,"\nDUMP     file=%s\n",dfnam) ;

    dump_file = fopen(dfnam,"w") ;

    /* oops */
    if(dump_file==NULL) {
      fprintf(stderr,"error opening dump file\n") ;
      exit(2) ;
    }

    dump(dump_file) ;
    fclose(dump_file) ;

    dump_cnt++ ;
  }

  /* make binary output files */
  if(call_code == BINARY_OUT ||
     call_code == INIT_OUT ||
     call_code == FINAL_OUT  || alloutputflag==1 ) {

    sprintf(binnam,"bin_XX_%04d",bin_cnt);
    fprintf(stderr,"BINARY   files=%s time %10.5g nstep=%d\n",
	    binnam,t,nstep) ;

    sprintf(binnam,"bin_b1_%04d",bin_cnt) ; bdump_primitive(binnam,M01) ;
    sprintf(binnam,"bin_b2_%04d",bin_cnt) ; bdump_primitive(binnam,M02) ;
    sprintf(binnam,"bin_b3_%04d",bin_cnt) ; bdump_primitive(binnam,M03) ;
    sprintf(binnam,"bin_e1_%04d",bin_cnt) ; bdump_primitive(binnam,F10) ;
    sprintf(binnam,"bin_e2_%04d",bin_cnt) ; bdump_primitive(binnam,F20) ;
    sprintf(binnam,"bin_e3_%04d",bin_cnt) ; bdump_primitive(binnam,F30) ;

    /* Dump the antisymmetric contravariant Maxwell tensor */
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<ii;jj++) {
      sprintf(binnam,"bin_Mcon_%02d_%02d_%04d",ii,jj,bin_cnt) ;
      ZLOOP  {
	get_geometry(i,j,CENT,&geom);
	get_state(p[i][j],&state,&geom) ;
	temp[i][j]=state.Mcon[ii][jj] ;
      }
      bdump_field(binnam,temp) ;
    }

    /* Dump the antisymmetric covariant Maxwell tensor */
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<ii;jj++) {
      sprintf(binnam,"bin_Mcov_%02d_%02d_%04d",ii,jj,bin_cnt) ;
      ZLOOP  {
	get_geometry(i,j,CENT,&geom);
	get_state(p[i][j],&state,&geom) ;

	/* get covariant Maxwell from contravariant */
	lower_T(state.Mcon,Mcov,&geom) ;
	temp[i][j]=Mcov[ii][jj] ;
      }
      bdump_field(binnam,temp) ;
    }


    /* Dump the mixed contravariant-covariant stress tensor */
    for(ii=0;ii<NDIM;ii++) for(jj=0;jj<NDIM;jj++) {
      sprintf(binnam,"bin_T_%02d_%02d_%04d",ii,jj,bin_cnt) ;
      ZLOOP  {
	get_geometry(i,j,CENT,&geom);
	get_state(p[i][j],&state,&geom) ;
	temp[i][j]=state.T[ii][jj] ;
      }
      bdump_field(binnam,temp) ;
    }

    bin_cnt++ ;
  }
  if(terminationflag==1){exit(122);}
}

void gridoutput(char *binnam)
{
  /* Binary output of geometry data */

  double X[NDIM] ;
  double x1[INA], x2[JNA], r[INA], th[JNA], placeholder ;
  int i, j, k ;
  int ibegin, iend, jbegin, jend ;

#ifdef DEBUG
  ibegin=ISG ;
  iend  =IEG ;
  jbegin=JSG ;
  jend  =JEG ;
#else
  ibegin=IS ;
  iend  =IE ;
  jbegin=JS ;
  jend  =JE ;
#endif

  j=JS;
  for (i=ibegin;i<=iend;i++) {
    coord(i,j,CENT,X) ;
    bl_coord(X,&(r[i]),&placeholder) ;
    x1[i]=X[1] ;
  }

  i=IS;
  for (j=jbegin;j<=jend;j++) {
    coord(i,j,CENT,X) ;
    bl_coord(X,&placeholder,&(th[j])) ;
    x2[i]=X[2] ;
  }

  sprintf(binnam,"b_r")  ; bdump_array(binnam,INA,r, ibegin,iend) ;
  sprintf(binnam,"b_x1") ; bdump_array(binnam,INA,x1,ibegin,iend) ;
  sprintf(binnam,"b_th") ; bdump_array(binnam,JNA,th,jbegin,jend) ;
  sprintf(binnam,"b_x2") ; bdump_array(binnam,INA,x2,jbegin,jend) ;

  for(int loc=0;loc<NPG;loc++){
    sprintf(binnam,"b_gdet_%04d",loc) ; bdump_gdet(binnam,loc) ;
    DLOOP {
      sprintf(binnam,"b_gcov_%04d_%04d_%04d",loc,j,k) ;
      bdump_gcov(binnam,loc,j,k) ;
    }
    DLOOP {
      sprintf(binnam,"b_gcon_%04d_%04d_%04d",loc,j,k) ;
      bdump_gcon(binnam,loc,j,k) ;
    }
  }
  for(i=0;i<NDIM;i++){
    DLOOP {
      sprintf(binnam,"b_conn_%04d_%04d_%04d",i,j,k) ;
      bdump_conn(binnam,i,j,k) ;
    }
  }
}


#if 0

/** some diagnostic routines -- not currently implemented**/

/* map out region around failure point */
void area_map(int i, int j, double prim[][JNA][NPR])
{
	int k ;

	fprintf(stderr,"area map\n") ;

	PLOOP {
		fprintf(stderr,"variable %d \n",k) ;
		fprintf(stderr,"i = \t %12d %12d %12d\n",i-1,i,i+1) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j+1,
				prim[i-1][j+1][k],
				prim[i][j+1][k],
				prim[i+1][j+1][k]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j,
				prim[i-1][j][k],
				prim[i][j][k],
				prim[i+1][j][k]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j-1,
				prim[i-1][j-1][k],
				prim[i][j-1][k],
				prim[i+1][j-1][k]) ;
	}

	/* print out other diagnostics here */

}

/* evaluate fluxed based diagnostics; put results in
 * global variables */

void diag_flux(double F1[][JNA][NPR], double F2[][JNA][NPR])
{
	int j ;

        mdot = edot = ldot = 0. ;
        JLOOP {
	        mdot += F1[IS][j][RHO]*2.*M_PI*dx[2] ;
                edot -= (F1[IS][j][UU] - F1[IS][j][RHO])*2.*M_PI*dx[2] ;
                ldot += F1[IS][j][U3] *2.*M_PI*dx[2] ;
        }
}

#endif /* 0 */
