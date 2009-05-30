
#include "decs.h"

double slope_lim(const double yy1, const double yy2, const double yy3) 
{
	double Dqm,Dqp,Dqc,s ;

	/* woodward, or monotonized central, slope limiter */
	if(lim == MC) {
		Dqm = 2. *(yy2 - yy1) ;
		Dqp = 2. *(yy3 - yy2) ;
		Dqc = 0.5*(yy3 - yy1) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else {
			if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
				return(Dqm) ;
			else if(fabs(Dqp) < fabs(Dqc))
				return(Dqp) ;
			else
				return(Dqc) ;
		}
	}
	/* van leer slope limiter */
	else if(lim == VANL) {
		Dqm = (yy2 - yy1) ;
		Dqp = (yy3 - yy2) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else
			return(2.*s/(Dqm+Dqp)) ;
	}
	/* minmod slope limiter (crude but robust) */
	else if(lim == MINM) {
		Dqm = (yy2 - yy1) ;
		Dqp = (yy3 - yy2) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else if(fabs(Dqm) < fabs(Dqp)) return Dqm ;
		else return Dqp ;
	}

	fprintf(stderr,"unknown slope limiter\n") ;
	exit(10) ;

	return(0.) ;
}

