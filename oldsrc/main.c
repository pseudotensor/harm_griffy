//========================= GRFFDE =========================//

#include "decs.h"
#include "defs.h"

int main(int argc,char *argv[])
{
	double tdump,tlog,tbinary ;
	 
	nstep = 0 ;
	dump_cnt = 0. ;

	/* perform initializations */
	init() ;

	/* do initial diagnostics */
	diag(INIT_OUT) ;

	tdump = t+DTd ;
	tlog = t+DTl ;
        tbinary = t+DTb ;

	/* tell me what you're printing out! */
	fprintf(stderr,"time,   dt,  nstep\n\n") ;

	while(t < tf) {

		/* step variables forward in time */
		step_ch() ;

		nstep++ ;

		fprintf(stderr,"%10.5g %10.5g %8d \n",t,dt,nstep) ;

		/* perform diagnostics */
		if(t >= tdump) {
			diag(DUMP_OUT) ;
			tdump += DTd ;
		}
		if(t >= tlog) {
			diag(LOG_OUT) ;
			tlog += DTl ;
			}
		if(t >= tbinary) {
			diag(BINARY_OUT) ;
			tbinary += DTb ;
		}
		/* don't go on forever */
		//if(nstep >= 20000) {
		//   fprintf(stderr,"too many steps!\n\n") ;
		//   exit(0) ;
		//}
	}
	fprintf(stderr,"done!\n\n") ;

	/* do final diagnostics */
	diag(FINAL_OUT) ;

	return(0) ;
}
