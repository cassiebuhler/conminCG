/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

#include <math.h>
#include <stdlib.h>

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

double dotprod(		/* inner product between n-vectors x and y */
	double	*x,
	double	*y,
	int	n
)
{
	int i; 
	double dotprod=0.0e0;
	int m = n % 5;

	/* loop-unrolling */
	/* add remainders from beginning then add the last n/5 groups of 5 */
	for (i=0; i<m; i++ ) dotprod = dotprod + x[i]*y[i];
	for (i=m; i<n; i+=5) dotprod = dotprod + x[i]*y[i]+x[i+1]*y[i+1]+x[i+2]*y[i+2]+x[i+3]*y[i+3]+x[i+4]*y[i+4];

	return (dotprod);
}

double dotprods(		/* inner product between n-vectors x and y divided by scalar s */
	double	*x,
	double	*y,
	int	n,
	double	s
)
{
	int i; 
	double dotprod=0.0e0;
	int m = n % 5;

	/* loop-unrolling */
	/* add remainders from beginning then add the last n/5 groups of 5 */
	for (i=0; i<m; i++ ) dotprod = dotprod + x[i]*y[i]/s;
	for (i=m; i<n; i+=5) dotprod = dotprod + x[i]*y[i]/s+x[i+1]*y[i+1]/s+x[i+2]*y[i+2]/s+x[i+3]*y[i+3]/s+x[i+4]*y[i+4]/s;

	return (dotprod);
}

double dotprodm(	/* inner product between n-vectors x and y */
	double	*x,
	double	*y,
	int	n,
	int     *mask,
	int	val
)
{
	int i; 
	double dotprod=0.0e0;

	for (i=0; i<n; i++) 
	    if ( mask[i] & val ) {
		dotprod += x[i]*y[i];
	    }

	return (dotprod);
}

double mincompl(	/* min-complementarity between n-vectors x and y */
	double  min,
	double	*x,
	double	*y,
	int	n,
	int     *mask,
	int	val
)
{
	int i; 
	double newmin=min;

	for (i=0; i<n; i++) 
	    if ( mask[i] & val ) {
		newmin = MIN(newmin, x[i]*y[i]);
	    }

	return (newmin);
}

void smx(		/* y = sparse matrix (A,kA,iA) times x */
	int	m,
	int	n,
	double	*A,
	int	*kA,
	int	*iA,
	double	*x,
	double	*y
)    
{
	int i,j,k;

	for (i=0; i<m; i++) y[i] = 0.0e0;
	for (j=0; j<n; j++) 
		for (k=kA[j]; k<kA[j+1]; k++)
			y[iA[k]] += A[k]*x[j];
}

void atnum(		/* (kAt,iAt,At)	= transpose of (kA,iA,A) */
	int	m,
	int	n,
	int	*kA,
	int	*iA,
	double	*A,
	int	*kAt,
	int	*iAt,
	double	*At
)  
{
	int i,j,k,row,addr;
	int *iwork;

	CALLOC(	iwork, m, int );

	for (k=0; k<kA[n]; k++)	{
		row = iA[k];
		iwork[row]++;
	}
	kAt[0] = 0;
	for (i=0; i<m; i++) {
		kAt[i+1] = kAt[i] + iwork[i];
		iwork[i] = 0;
	}
	for (j=0; j<n; j++) {
		for (k=kA[j]; k<kA[j+1]; k++) {
			row = iA[k];
			addr = kAt[row]	+iwork[row];
			iwork[row]++;
			iAt[addr] = j;
			At[addr]  = A[k];
		}
	}
	FREE( iwork );
}

double maxv(		/* compute componentwise maximum of n-vector x */
	double *x,
	int	n
)
{
	int i;
	double maxv=0.0e0;

	for (i=0; i<n; i++) maxv = MAX(maxv, fabs(x[i]));

	return (maxv);
}

#ifdef __cplusplus
}
#endif
