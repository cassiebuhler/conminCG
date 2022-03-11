#define DECLAREALL \
\
int n; \
double *c; \
double f; \
double *x; \
double *dx; \
double *x0; \
int max; \
double inftol; \
int itnlim; \
double timlim; \
int quadratic; \
int verbose; \
int (*stopping_rule)(void *); \
void (*init_vars)(void *); \
void (*h_init)(void *); \
int  (*h_update)(void *); \
void (*h_close)(void *); \
int iter; \
double elaptime; \
int	ifun; \
double	starttime; \
double *temp; \
int j; \
double acc; \
int nrst; \
int restart; \
int status; \
double u1, u2, u3, u4; \
double xTx, cTc, ytTyt, dxTc, cTyt; \
double *p, *pt, *y, *yt, *c0; \
double alpha, f0, fmin, dg, ap, fp, dp; \
double dalpha, at; \
int linesearch, tried, nerror, interp, pertcnt; \
double *temp1, *temp2, *c00; \
double lambda, a, b, d, e, bracket, denom; \
double a2, b2, d2, e2, u10, u11, u12, u13, u14, u15; \
double cTct, pTc, yTc, ptTy, ytTy, ptTp, ytTp; \
double alpha0, lambda0; \
double *dx0; \
int restart0; \
double prat, dprat, olambda, ddprat, oolambda, dprat0, prat0;

#define	COPYBACK(lp) \
\
n = lp->n; \
c = lp->c; \
f = lp->f; \
x = lp->x; \
max = lp->max; \
inftol = lp->inftol; \
itnlim = lp->itnlim; \
timlim = lp->timlim; \
quadratic = lp->quadratic; \
stopping_rule = lp->stopping_rule; \
init_vars = lp->init_vars; \
h_init = lp->h_init; \
h_update = lp->h_update; \
h_close = lp->h_close; 
