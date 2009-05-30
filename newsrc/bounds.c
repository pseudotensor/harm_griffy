//========================= GRFFDE =========================//

#include "decs.h"

/* bound array containing entire set of primitive variables */

void bound_prim( double prim[INA][JNA][NPR] )
{
	int i,j,k ;

	/* inner r boundary condition */
	JLOOP {

		prim[IS-1][j][M01] = prim[IS][j][M01] *
			gdet[IS][j][CENT]/gdet[IS-1][j][CENT] ;
		prim[IS-2][j][M01] = prim[IS][j][M01] *
			gdet[IS][j][CENT]/gdet[IS-2][j][CENT] ;

		prim[IS-1][j][M02] = prim[IS][j][M02] ;
		  /* * gdet[IS][j][CENT]/gdet[IS-1][j][CENT]*/
		prim[IS-2][j][M02] = prim[IS][j][M02] ;
		  /* * gdet[IS][j][CENT]/gdet[IS-2][j][CENT]*/

		prim[IS-1][j][M03] = prim[IS][j][M03] ;
		  /* * gdet[IS][j][CENT]/gdet[IS-1][j][CENT]*/
		prim[IS-2][j][M03] = prim[IS][j][M03] ;
		  /* * gdet[IS][j][CENT]/gdet[IS-2][j][CENT]*/

		prim[IS-1][j][F10] = prim[IS][j][F10]
		  * gdet[IS][j][CENT]/gdet[IS-1][j][CENT] ;
		prim[IS-2][j][F10] = prim[IS][j][F10]
		  * gdet[IS][j][CENT]/gdet[IS-2][j][CENT] ;

		prim[IS-1][j][F20] = prim[IS][j][F20]
		  /* * gdet[IS][j][CENT]/gdet[IS-1][j][CENT]*/ ;
		prim[IS-2][j][F20] = prim[IS][j][F20]
		  /* * gdet[IS][j][CENT]/gdet[IS-2][j][CENT]*/ ;

		prim[IS-1][j][F30] = prim[IS][j][F30]
		  /* * gdet[IS][j][CENT]/gdet[IS-1][j][CENT]*/ ;
		prim[IS-2][j][F30] = prim[IS][j][F30]
		  /* * gdet[IS][j][CENT]/gdet[IS-2][j][CENT]*/ ;

	}

	/* outer BC */
	JLOOP {
		prim[IE+1][j][M01] = prim[IE][j][M01] *
			gdet[IE][j][CENT]/gdet[IE+1][j][CENT] ;
		prim[IE+2][j][M01] = prim[IE][j][M01] *
			gdet[IE][j][CENT]/gdet[IE+2][j][CENT] ;

		prim[IE+1][j][M02] = prim[IE][j][M02]
		  /* * gdet[IE][j][CENT]/gdet[IE+1][j][CENT]*/ ;
		prim[IE+2][j][M02] = prim[IE][j][M02]
		  /* * gdet[IE][j][CENT]/gdet[IE+2][j][CENT]*/ ;

		prim[IE+1][j][M03] = prim[IE][j][M03]
		  /* * gdet[IE][j][CENT]/gdet[IE+1][j][CENT]*/ ;
		prim[IE+2][j][M03] = prim[IE][j][M03]
		  /* * gdet[IE][j][CENT]/gdet[IE+2][j][CENT]*/ ;

		prim[IE+1][j][F10] = prim[IE][j][F10] *
			gdet[IE][j][CENT]/gdet[IE+1][j][CENT] ;
		prim[IE+2][j][F10] = prim[IE][j][F10] *
			gdet[IE][j][CENT]/gdet[IE+2][j][CENT] ;

		prim[IE+1][j][F20] = prim[IE][j][F20]
		  /* * gdet[IE][j][CENT]/gdet[IE+1][j][CENT]*/ ;
		prim[IE+2][j][F20] = prim[IE][j][F20]
		  /* * gdet[IE][j][CENT]/gdet[IE+2][j][CENT]*/ ;

		prim[IE+1][j][F30] = prim[IE][j][F30]
		  /* * gdet[IE][j][CENT]/gdet[IE+1][j][CENT]*/ ;
		prim[IE+2][j][F30] = prim[IE][j][F30]
		  /* * gdet[IE][j][CENT]/gdet[IE+2][j][CENT]*/ ;

	}

	/* 
	 * polar BCs  -- this is the Z-axis 
	 * just copy
	 */

	for(i=ISG;i<=IEG;i++) PLOOP {
		prim[i][JS-1][k] = prim[i][JS][k] ;
		prim[i][JE+1][k] = prim[i][JE][k] ;
		prim[i][JS-2][k] = prim[i][JS+1][k] ;
		prim[i][JE+2][k] = prim[i][JE-1][k] ;
	}

	/* 
	 *make sure b and e are antisymmetric at the poles
	 * (theta components) 
	 */

	for(i=ISG;i<=IEG;i++) {
		for(j=JSG;j<JS;j++) {
		  prim[i][j][M02] *= -1. ;
		  prim[i][j][F20] *= -1. ;
		}
		for(j=JE+1;j<=JEG;j++) {
		  prim[i][j][M02] *= -1. ;
		  prim[i][j][F20] *= -1. ;
		}
	}

}
