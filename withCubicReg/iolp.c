/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>

#ifndef PC
#include <unistd.h>
#include <fcntl.h>
#endif

#ifdef QuadPrec
#include "Quad/Quad.h"
#define double Quad
#else
#define high(x) (x)
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "conmin.h"
#include "myalloc.h"

/* Define functions */

CONMIN	*openlp(void)
{
	CONMIN *lp;

	MALLOC(lp, 1, CONMIN);

	if (lp != NULL) {
		lp->n = 0;
		lp->f = 0.0;
		lp->x = NULL;
		lp->c = NULL;
		lp->max	    = 1;     /*	max = -1,   min	= 1 */
		lp->itnlim  = 500;   /* iteration limit */
		lp->timlim  = HUGE_VAL; /* time limit */
		lp->verbose = 2;     /*	verbosity level	*/
		lp->inftol  = 1.0e-6;/*	infeasibility requested	*/
		lp->quadratic = 0; /* assert problem is QP */
		lp->init_vars	  = deflt_init;
		lp->pertlim = 10;  /* limit on number of lambda updates for cubic regularization, U in Algorithm 2 */
/*
		lp->h_init        = nl_init_mps;
		lp->h_update 	  = nl_update_dummy;
		lp->h_close  	  = nl_close_mps;
		lp->objval = objval_dummy;
		lp->objgrad = objgrad_dummy;
*/
	}

	return lp;
}

void	closelp(
	CONMIN	*lp
)
{

	FREE(lp->c);
	FREE(lp->x);
	FREE(lp);

}

#ifdef __cplusplus
}
#endif
