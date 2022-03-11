/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>


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
#include "conmin1.h"
#include "myalloc.h"

#define	MAXIT  200
#define EPS 1.0e-15

/* Table of constant values */
static int sigint = FALSE;

void handler(
    int sig
)
{
    if (sig == SIGINT) {sigint = TRUE;}
}

/* Prototype static functions */

static void error(
	int	num,
	char	*text
);

void deflt_init(
	void *vlp
);

void opt_message(
	int status
);

void deflt_hook(
	void *vlp
);

/* Define functions */

int	solvelp(
	CONMIN	*lp
)
{
	DECLAREALL

	signal(SIGINT, handler);

	/*---------------------------------------------------------+
	| transcribe information from lp			  */

	n	= lp->n;
	c	= lp->c;
	f	= lp->f;
	x	= lp->x;

	max    	   = lp->max;
	inftol	   = lp->inftol;
	verbose	   = lp->verbose;

	stopping_rule =	lp->stopping_rule;
	init_vars     =	lp->init_vars;

	h_init       = lp->h_init;
	h_update     = lp->h_update;
	h_close      = lp->h_close;

	/*---------------------------------------------------------+
	| Allocate space for unallocated arrays.		  */

	if (c == NULL){CALLOC( c, n, double );}
	lp->c = c;

	/*------------------------------------------------------+
	| Here is an initialization hook.                      */

	h_init(lp);
	COPYBACK(lp);

	if (verbose>1) {
	    printf("variables: %8d\n", n);
	}

	/*---------------------------------------------------------+
	| allocate storage.					  */

	MALLOC( dx,	n,	double );
	MALLOC( x0,	n,	double );
	MALLOC( temp,	n,	double );
	MALLOC( temp1,	n,	double );
	MALLOC( temp2,	n,	double );
	MALLOC( p, 	n, 	double );
	MALLOC( pt, 	n, 	double );
	MALLOC( y, 	n, 	double );
	MALLOC( yt, 	n, 	double );
	MALLOC( c0,	n,	double );
	MALLOC( c00,	n,	double );
	MALLOC( dx0,	n,	double );

	/*------------------------------------------------------+
	| initialize parameters.				*/

	iter=0; 
	elaptime = 0.0;
	starttime = cputimer();
	ifun = 0;		
	pertcnt = 0;
	acc = 1e-20;
	alpha = 0.0;
	dg = 0.0;
	if (verbose>1) {
		printf(
"-----------------------------------------------------\n"
"Iter |		Obj Value     	Residual	| P   \n"
"- - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
		);
	}

	/*------------------------------------------------------+
	| initialize x						*/

	init_vars(lp);
	COPYBACK(lp);
	dprat = 0.0;
	ddprat = 0.0;

restart:
	nrst = n;
	restart = 0;
	lambda = 0.0;

	while ( iter < itnlim && elaptime < timlim ) {

		xTx = dotprod(x, x, n); 
		cTc = dotprod(c, c, n);

		/* Check for convergence */
		if ( cTc <= inftol*inftol*MAX(1.0, xTx) ) { 
			status = 0;
			goto end;
		}


		/* Compute step direction */
		if ( restart == 0 ) {
			if ( verbose > 1 ) {
				printf("%4d :\t%14.6e\t %14.6e\t | %d \n", iter, f, sqrt(cTc/MAX(1.0, xTx)), pertcnt);
			}
			for (j = 0; j<n; j++) dx[j] = -c[j];
		} else {
			/* Test to see if the Powell restart criterion holds */
printf("%14.6e \t %14.6e\n", lambda, ABS(dotprod(c, c0, n)/cTc));
			if ( ABS(dotprod(c, c0, n)/cTc) > 0.2 && restart0 > 1 && nrst != n ) { 
/*
			if ( ABS(dotprod(c, c0, n)/(sqrt(dotprod(c,c,n))*sqrt(dotprod(c0,c0,n)))) > 0.5 && restart0 > 1 && nrst != n ) {
*/
printf("Powell restart\n"); 
/*
				nrst = n;
*/

if ( lambda == 0.0 ) {
				if ( dprat0 != 0.0 ) ddprat = ( ABS(dotprod(c, c0, n)/cTc) - 2*prat + prat0 ) / ( (lambda - olambda)*(oolambda-olambda) );
				/* if ( dprat0 != 0.0 ) ddprat = ( dprat - dprat0 ) / (lambda/2.0 - oolambda/2.0); */
				dprat0 = dprat;
				if ( lambda > 0 ) dprat = ( ABS(dotprod(c, c0, n)/cTc) - prat ) / ( lambda - olambda );
				prat0 = prat;
				prat = ABS(dotprod(c, c0, n)/cTc);
				alpha = alpha0;
				restart = restart0;
				for ( j=0; j<n; j++) x[j] = x0[j];
				for ( j=0; j<n; j++) dx[j] = dx0[j];
				for ( j=0; j<n; j++) c0[j] = c00[j];
				nerror = h_update(lp);
				COPYBACK(lp);
				xTx = dotprod(x, x, n);
				cTc = dotprod(c, c, n);
				oolambda = olambda;
				olambda = lambda;
				if ( lambda == 0.0 ) {
					lambda = prat/0.2; 
				} else {
					if ( pertcnt <= 2 ) {
						lambda = olambda - prat/dprat;
					} else {
						if ( dprat*dprat - 2*ddprat*prat > 0 && ABS(ddprat) > 1e-8 ) {
							lambda = olambda + MAX( (-dprat + sqrt(dprat*dprat - 2*ddprat*prat))/ddprat,
										(-dprat - sqrt(dprat*dprat - 2*ddprat*prat))/ddprat ); 
						} else {
							lambda = 2*olambda;
						}
					}
					if ( lambda < 1e-12 ) lambda = 2*olambda;
				}
} else {
	lambda = 2*lambda;
}
				pertcnt++;
				iter--;
			} else {
				for ( j=0; j<n; j++ ) dx0[j] = dx[j];
				lambda = 0.0;
				dprat = 0.0;
				dprat0 = 0.0;
				prat0 = 0.0;
				if ( verbose > 1 ) {
					printf("%4d :\t%14.6e\t %14.6e\t | %d \n", iter, f, sqrt(cTc/MAX(1.0, xTx)), pertcnt);
				}
				pertcnt = 0;
			}

			/* If performing a restart, update the Beale restart vectors */
			if ( nrst == n ) {
				for ( j=0; j<n; j++ ) {
					pt[j] = alpha*dx[j];
					yt[j] = c[j] - c0[j];
					ytTyt = dotprod(yt, yt, n);
					cTyt = dotprod(pt, yt, n);
					cTct = dotprod(pt, pt, n);
				}
			}

			for ( j=0; j<n; j++ ) {
				p[j] = alpha*dx[j];
				y[j] = c[j] - c0[j];
			}

			pTc = dotprod(pt, c, n);
			yTc = dotprod(yt, c, n);

			u1 = -dotprods(pt, c, n, ytTyt);
			u2 = 0.0;
			for (j=0; j<n; j++) u2 = u2 + 2*pt[j]*c[j]/cTyt - yt[j]*c[j]/ytTyt;
			u3 = cTyt/ytTyt;

			bracket = lambda*lambda + 2*lambda*ytTyt/cTyt + ytTyt/cTct;
			a = -cTyt/(lambda*cTyt + ytTyt);
			b = (-lambda*cTyt - 2*ytTyt)*ytTyt/(cTyt*cTct*bracket*(lambda*cTyt + ytTyt));
			d = lambda/(bracket*(lambda*cTyt + ytTyt));
			e = ytTyt/(cTct*bracket*(lambda*cTyt + ytTyt));

			bracket = lambda*lambda+lambda/u3+lambda*ytTyt/cTyt+cTyt/(u3*cTct);
			denom = u3*u3*cTct*lambda*bracket + u3*cTct*bracket;
			d = u3*u3*lambda*(cTct/cTyt) / denom;
			/*
			a = -u3/(lambda*u3+1.0);
			b = -u3*lambda/denom - 1.0/denom - u3*ytTyt/(cTyt*denom);
			e = u3 / denom;
			*/

			for ( j=0; j<n; j++ ) dx[j] = a*c[j] + b*pTc*pt[j] + d*yTc*yt[j] + e*yTc*pt[j] + e*pTc*yt[j];

			if ( nrst != n ) {
				/* inv(H+lI)*y */
				ptTy = dotprod(pt, y, n);
				ytTy = dotprod(yt, y, n);
				for ( j=0; j<n; j++ ) temp1[j] = -(a*y[j] + b*ptTy*pt[j] + d*ytTy*yt[j] + e*ytTy*pt[j] + e*ptTy*yt[j]);

				/* inv(H+lI)*invH*p */
				a2 = 1.0/(lambda*u3+1.0);
				b2 = lambda*b;
				d2 = lambda*d;
				e2 = lambda*e;
				ptTp = dotprod(pt, p, n);
				ytTp = dotprod(yt, p, n);
				for ( j=0; j<n; j++ ) temp2[j] = a2*p[j] + b2*ptTp*pt[j] + d2*ytTp*yt[j] + e2*ytTp*pt[j] + e2*ptTp*yt[j];

				u10 = dotprod(temp1, c, n); /* yT inv(H+lI) g */
				u11 = dotprod(temp2, c, n); /* pT invH inv(H+lI) g */
				u12 = dotprod(temp1, y, n); /* yT inv(H+lI) y */
				u13 = dotprod(temp2, y, n); /* pT invH inv(H+lI) y */
				u14 = dotprod(temp2, p, n); /* pT invH inv(H+lI) p */
				u15 = dotprod(p, y, n);	/* pTy */

				denom = lambda*u14*(u15 + u12) + u13*u13;
				
				/* dx = -inv(H+lI)*g --- Remember the minus sign */
				for ( j=0; j<n; j++ ) dx[j] = dx[j] + lambda*u14*u10*temp1[j]/denom - (u15 + u12)*u11*temp2[j]/denom
								+ u13*u10*temp2[j]/denom + u13*u11*temp1[j]/denom;

			}

		}

		/* Check that the search direction is a descent direction */
		dxTc = dotprod(dx, c, n);
		if ( dxTc > 0 ) {
			status = 3;
			goto end;
		}

		/* Save the current point */
		f0 = f;
		for ( j=0; j<n; j++ ) {
			x0[j] = x[j];
			c00[j] = c0[j];
			c0[j] = c[j];
		}
		alpha0 = alpha;
		restart0 = restart;

		if ( restart == 0 ) {
			restart = 1;
		} else {
			if ( pertcnt == 0 ) {
				if ( nrst == n ) nrst = 0;
				nrst++;
				restart = 2;
			}
		}

		/* Compute the steplength */
		fmin = f;
		alpha = alpha*dg/dxTc;
		if ( nrst == 1 ) alpha = 1.0;
		if ( restart <= 1 ) alpha = 1.0/sqrt(cTc);
		ap = 0;
		fp = fmin;
		dp = dxTc;
		dg = dxTc;
	
		linesearch = TRUE;
		tried = 0;

		while ( linesearch ) {

			if ( tried >= 300 ) {
				linesearch = FALSE;
				continue;
			}
		
			/* If the linear search failed, restart if possible.  Otherwise, declare error and exit. */	
			if ( alpha*dotprod(dx, dx, n) < acc ) {
				for ( j=0; j<n; j++ ) x[j] = x0[j];
				nerror = h_update(lp);
				COPYBACK(lp);
				if ( restart <= 1 ) {
					status = 2;
					goto end;
				} else {
					goto restart;
				}
			}

			for ( j=0; j<n; j++ ) x[j] = x0[j] + alpha*dx[j];
			nerror = h_update(lp);
			COPYBACK(lp);
			tried++;
			dalpha = dotprod(c, dx, n);

/* HANDE: Changed the threshold for dalpha to -1e-12 from -1e-16 to force the issue.  The signs of dalpha and dp prevent us from 
 * going to the correct part of the bounding process later.  The signs come out to be the same, even though we have bracketed the minimum. */
			if ( f > (1.0+1e-12)*fmin && dalpha < -1e-12 ) {
			/* If the function value is higher at alpha, but we have negative slope,
			 * it means we've already passed the minimum and gone too far.  Search for the min
			 * between 0 and alpha/3 */
				alpha = alpha/3.0;
				ap = 0;
				fp = fmin;
				dp = dg;
				continue;
			} 

			interp = FALSE;				

/* HANDE: What happens when dalpha and dg are so small that we keep on having problems? Clearly, we stop with errors in the next iteration
 * when we can't find argmin alpha.  But what if it takes forever to find that optimal value? */
			/* If no sufficient descent, continue with cubic interpolation */
			if ( f <= fmin + 1e-4*alpha*dg && ABS(dalpha/dg) < 0.9) {
				/* If sufficient descent and two points have been tested and/or minimum has been found, linesearch is done. */
				if ( tried >= 2 || ABS(dalpha/dg) < inftol ) {
					linesearch = FALSE;
				} else {
					interp = TRUE;
				}
			} else {
				interp = TRUE;
			}

			/* If we are still interpolating.... */
			if ( interp ) {

				/* apply Davidon's formula for the minimum of a cubic */
				u1 = dp + dalpha - 3.0*(fp - f)/(ap - alpha);
				u2 = u1*u1 - dp*dalpha;
				if ( u2 < 0 ) u2 = 0;
				u2 = sqrt(u2);
				at = alpha - (alpha - ap)*(dalpha + u2 - u1)/(dalpha - dp + 2.0*u2);

/* HANDE: ap, alpha seem switched, going past the optimum, etc. */
				if (dalpha/dp > 0) { /* minimum has not been bracketed */
					if (dalpha <= 0 || at <= 0 || at >= 0.99*MIN(alpha,ap)) {
						if (dalpha > 0 || at <= 1.01*MAX(ap, alpha)) {
							if (dalpha <= 0) {
								at = 2*MAX(ap,alpha);
							} else {
								at = MIN(ap,alpha)/2;
							}

						}	
					}
				} else { /* minimum has been bracketed */
					if ( at < 1.01*MIN(alpha,ap) || at > 0.99*MAX(alpha,ap) ) {
						at = (alpha + ap)/2;
					}

				}
				
				ap = alpha;
				fp = f;
				dp = dalpha;
				alpha = at;
			}
		}


		/* Take the step and update function value and gradient */
/*
		for ( j=0; j<n; j++ ) x[j] = x0[j] + alpha*dx[j];
		nerror = h_update(lp);
		COPYBACK(lp);
*/


		/* Housekeeping */
		iter = iter + 1;
		elaptime = cputimer() - starttime;

	}

	if ( iter == itnlim ) status = 1;
		

end:
	lp->iter = iter;
if ( verbose > 1 ) {
	printf("%4d :\t%14.6e\t %14.6e\t | %d \n", iter, f, sqrt(cTc/MAX(1.0, xTx)), pertcnt);
}
if (verbose >= 4) {
	printf("\n");
	printf("dx: \n");
	for (j=0; j<lp->n; j++) {
	    printf("%14.6e \n", dx[j]);
	}
	printf("\n");
}
	    FREE(dx);
	    FREE(x0);

	    h_close(lp);

	if (verbose>1) {
		printf("----------------------\n");
		opt_message( status );
	}

	return status;
}

/* Define static functions */

static void error(
	int	num,
	char	*text
)
{
	char str[250];

	switch (num) {
	case   2: sprintf(str, "cannot open file %s\n",text);
		  break;
	case   3: sprintf(str, "cannot read file %s\n",text);
		  break;
	case   4: sprintf(str, "cannot create file %s\n",text);
		  break;
	case   5: sprintf(str, "cannot write file %s\n",text);
		  break;
	case   6: sprintf(str, "cannot allocate space\n");
		  break;
	case   9: sprintf(str, "cannot solve dual as primal when ranges are present\n");
		  break;
	case  10: sprintf(str, "dimension conflict in %s\n",text);
		  break;
	case  11: sprintf(str, "NAME not found\n");
		  break;
	case  26: sprintf(str, "unrecognized section label: \n   %s \n",text);
		  break;
	case  31: sprintf(str, "negative ranges not allowed: %s \n",text);
		  break;
	case  40: sprintf(str, "multiple entry in matrix at %s \n"
		               "(these labels may be incorrect unless NOPREP"
			       " is set) \n", text);
		  break;
	case  50: sprintf(str, "cannot evaluate obj and/or constraint at "
			       "initial solution \n");
		  break;
	case  60: sprintf(str, "student version limited to "
			       "300 variables/constraints \n");
		  break;
	}
	printf("CONMIN ERROR(%d): %s\n", num, str); exit(1);
}

void deflt_init(
        void *vlp
)
{
        CONMIN    *lp= (CONMIN *)vlp;
        int     n;
        int     nerror;
        double  *x;

        int 	j;

        n = lp->n; 
	x = lp->x;

        if (x == NULL) {
                CALLOC(x, n, double);
		for (j=0; j<n; j++) x[j] = 0.0;
		lp->x = x;
        }

	nerror = lp->h_update(lp);
        if (nerror) { error(50,""); }

}

void opt_message( int status )
{
	switch ( status )
	{
	    case 0:
		printf("OPTIMAL SOLUTION FOUND\n");
		break;
	    case 1:
		printf("ITERATION LIMIT\n");
		break;
	    case 2:
		printf("LINEAR SEARCH FAILED TO IMPROVE THE FUNCTION VALUE.\n");
		break;
	    case 3:
		printf("SEARCH DIRECTION WAS NOT A DESCENT DIRECTION.  PLEASE INCREASE INFTOL.\n");
		break;
	    case 4:
		printf("TIME LIMIT\n");
		break;
	}
	/* fflush(stdout); */
}



#ifdef __cplusplus
}
#endif
