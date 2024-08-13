/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

#include <stdio.h>
/* #include <sigfpe.h> */

#define	TRUE 1
#define	FALSE 0

#define LP 0
#define QP 1

#define	MAX(x,y)  ((x) > (y) ? (x) : (y))
#define	MIN(x,y)  ((x) > (y) ? (y) : (x))
#define	ABS(x)	  ((x) > 0   ? (x) : -(x))
#define	SGN(x)	  ((x) > 0   ? (1.0) : (-1.0))

typedef	struct conmin {
	int n;		/* number of columns */
	double *c;	/* pointer to array containing objective function */
	double f;	/* fixed adjustment to objective function */
	double *x;	/* pointer to array containing primal solution */
	double *dx;	/* step direction for primal solution */

	int max;	/* max = -1, min = 1 */
	double inftol;	/* infeasibility tolerance */
	int itnlim;	/* iteration limit */
	double timlim;	/* time limit */
	int verbose;	/* level of verbosity */
	int quadratic;	/* assert problem is quadratic */

	int  (*stopping_rule)(void *);/* pointer to stopping	rule fcn */
	void (*init_vars)(void *);    /* pointer to initialization fcn */

	void (*h_init)(void *);       /* pointer to initialization hook fcn */
	int  (*h_update)(void *);     /* pointer to f,g,h update hook fcn */
	void (*h_close)(void *);      /* pointer to update hook fcn */
	void (*h_step)(void *);       /* pointer to step hook fcn */

	double  (*objval)       ( double * );
	void    (*objgrad)      ( double *, double * );
	void    (*hessian)      ( double *, double *, double * );
	void    (*conval)       ( double *, double * );
	void    (*congrad)      ( double *, double *, double * );

	int    iter;	    /* current iteration number	*/
	double elaptime;    /* elapsed time */
	double primal_obj;  /* primal objective	value */

	int	inpower;

} CONMIN;

CONMIN	*openlp(void);

int	solvelp(
	CONMIN	*lp
);

void	closelp(
	CONMIN	*lp
);

void	my_exit(
	int num,
	char *str
);

extern void message();

double sdotprod(
	double	*a,
	int	*ja,
	double  *x,
	int	na
);

double dotprod(		/* inner product between n-vectors x and y */
	double	*x,
	double	*y,
	int	n
);

double dotprods(		/* inner product between n-vectors x and y */
	double	*x,
	double	*y,
	int	n,
	double	s
);

double dotprodm(	/* inner product between n-vectors x and y */
	double	*x,
	double	*y,
	int	n,
	int	*mask,
	int	val
);

double compl(	/* mu-complementarity between n-vectors x and y */
	double  mu,
	double	*x,
	double	*y,
	int	n,
	int	*mask,
	int	val
);

double mincompl(	/* min-complementarity between n-vectors x and y */
	double  min,
	double	*x,
	double	*y,
	int	n,
	int     *mask,
	int	val
);

void smx(		/* y = sparse matrix (A,kA,iA) times x */
	int	m,
	int	n,
	double	*A,
	int	*kA,
	int	*iA,
	double	*x,
	double	*y
);

void smtx(		/* y = x times sparse matrix (A,kA,iA) */
	int	n,
	double	*A,
	int	*kA,
	int	*iA,
	double	*x,
	double	*y
);

void atnum(		/* (kAt,iAt,At)	= transpose of (kA,iA,A) */
	int	m,
	int	n,
	int	*kA,
	int	*iA,
	double	*A,
	int	*kAt,
	int	*iAt,
	double	*At
);

double maxv(		/* compute componentwise maximum of n-vector x */
	double *x,
	int	n
);

void deflt_init(
    void *
);

void inv_clo(void);

char *my_strdup(
	char	*s1
);

double cputimer(
);

void nlsetup(
    void *vlp
);

void nlobjterm(
        void (*func)( double *z, double *param,
		      double *pval, double *grad, double **hessian),
			   /* function that computes val, grad, hessian at z */
        int k,             /* number of arguments for func() */
        char **collabs,    /* col labels for arugments to func() */
	int np,	           /* number of parameters in param array */
	double *param	   /* parameters to pass to func() */
);

int rd_specs(char *Spec_name);
void *binsearch(char **sp);
void set_opns(CONMIN *lp);

double objval_dummy( double *x );
void objgrad_dummy( double *c, double *x );
void hessian_dummy( double *Q, double *x, double *y );
void conval_dummy( double *h, double *x );
void congrad_dummy( double *A, double *At, double *x );
