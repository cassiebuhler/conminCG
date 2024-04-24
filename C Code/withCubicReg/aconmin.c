/****************************************************************
Copyright (C) 1997-2000 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
#include "conmin.h"
#include "myalloc.h"
#ifdef __cplusplus
	}
#endif

#define STDIO_H_included
#include "asl_pfgh.h"
#include "getstub.h"

 static void nlclose ANSI((void*));
 static void nlinit ANSI((void*));
 static int  nlupdate ANSI((void*));
 static void qx (int,double*,int*,int*,double*,double*);
 static void getpenaltyparam ANSI((void*));

#ifdef __cplusplus
extern "C" {
#endif

 extern char conmin_bsname[], conmin_version[];

 struct options
{
	int maxflag;	/* max or min */

	int itnlim;	/* iteration limit */
	int pertlim;	/* perturbation limit */
	double timlim;	/* runtime limit */
	double inftol;	/* convergence tolerance */

	int v;		/* verbosity level */
};

 static struct options opns = { 0, 500, 10, 0. /* HUGE_VAL */, 1e-6, 0 };

 static int neval, ngeval, time_flag, wantver;
 static real ftimes[5];	/* func, grad, constr, Jacobians, Hessians */

 static I_Known maximize = { -1, &opns.maxflag };

 static keyword
keywds[] = {	/* must be alphabetical */
 KW("inftol",	D_val,	&opns.inftol,	"Convergence tolerance"),
 KW("iterlim",	I_val,	&opns.itnlim,	"Iteration Limit (ITERLIM)"),
 KW("max",	IK_val,	&maximize,	"maximize the objective"),
 KW("maximize",	IK_val,	&maximize,	"maximize the objective"),
 KW("maxit",	I_val,	&opns.itnlim,	"Iteration Limit (ITERLIM)"),
 KW("min",	IK1_val,&opns.maxflag,	"minimize the objective"),
 KW("minimize",	IK1_val,&opns.maxflag,	"minimize the objective"),
 KW("outlev",	I_val,	&opns.v,	"output level"),
 KW("pertlim",	I_val,	&opns.pertlim,	"Perturbation Limit"),
 KW("timing",	I_val,	&time_flag,	"timing destination: 1 = stdout, 2 = stderr, 3 = both"),
 KW("timlim",	D_val,	&opns.timlim,	"Time limit"),
 KW("verbose",	I_val,	&opns.v,	"synonym for outlev"),
 KW("version",	I_val,	&wantver,	"report version"),
 KW("wantsol",	WS_val, 0,		WS_desc_ASL+5),
 };

 static char *usage_msg[] = {
	"  where  stub  is from  \"ampl -obstub\"  or  \"ampl -ogstub\".",
	"  -m with no stub ==> read an MPS file on stdin.",
	0 };

 static Option_Info Oinfo = {
	"conmin", conmin_bsname, "conmin_options", keywds, nkeywds, 1,
	conmin_version, usage_msg, 0, 0, 0, 1 };

#undef INTEGER
#undef DOUBLE
#define INTEGER(x) (int *)Malloc((x)*Sizeof(int))
#define DOUBLE(x) (double *)Malloc((x)*Sizeof(double))

 static int
qtest(char *s) {
	ASL *asl;
	int flag, i;
	fint *colqp, *rowqp;
	real *delsqp;
	FILE *nl;

	asl = ASL_alloc(ASL_read_fg);
	nl = jac0dim(s, (long)strlen(s));
	flag = nlc == 0;
	if (nlo == 0 || flag == 0) {
		fclose(nl);
		goto done;
		}

	qp_read(nl,0);
	if ( nlo > 0 && nqpcheck(i, &rowqp, &colqp, &delsqp) < 0)
		flag = 0;
done:
	ASL_free(&asl);
	return flag;
	}

 int
#ifdef KR_headers
amplin(asl, argv, stub, kp) ASL *asl; char **argv, *stub; CONMIN *kp;
#else
amplin(ASL *asl, char **argv, char *stub, CONMIN *kp)
#endif
{
	fint M, N, NZ, NO, MXROW, MXCOL;
	double *c;
	cgrad *cg, **cgx;
	ograd *og;
	int j, n;
	FILE *nl;
	double *x;

	if (!stub)
		usage_ASL(&Oinfo, 1);

	/* cur_ASL is used implicitly by nlinit and the nlupdate */
	/* routines.  It was changed by qtest.  */
	cur_ASL = asl;
	want_xpi0 = 3;
	nl = jacdim(stub, &M, &N, &NO, &NZ, &MXROW, &MXCOL, (ftnlen)strlen(stub));
	Uvx = DOUBLE(N);
	pfgh_read(nl, ASL_findgroups | ASL_find_co_class);

	/*
	if (qtest(stub)) {
		kp->quadratic = 1;
		if (opns.v) {
			printf("It's a QP.\n");
			need_nl = 0;
		}
	}
	*/

	kp->c = c = DOUBLE(N);
	memset(c, 0, N*sizeof(double));
	kp->f = 0;
	if ( n_obj > 1 ) {
		printf("Cannot handle multiple objective functions.  Exiting. \n");
		return 1;
	}

	kp->n = N;
	if (j = nlogv + niv + nlvbi + nlvci + nlvoi) {
		printf("ignoring integrality of %d variables\n", j);
		need_nl = 0; 
	}

	if (X0) kp->x = X0; else CALLOC( kp->x, N, double ) /* macro gives { ... } */

	return 0;
}
	
 void
#ifdef KR_headers
amplout(asl, kp, status) ASL *asl; CONMIN *kp; int status;
#else
amplout(ASL *asl, CONMIN *kp, int status)
#endif
{
	char buf[32], hbuf[256];
	int i;
	double *y, *ye;
	typedef struct { char *msg; int code; } Sol_info;
	static Sol_info solinfo[] = {
		{ /* 0 */ "optimal solution", 0 },
		{ /* 1 */ "iteration limit", 100 },
		{ /* 2 */ "suboptimal solution", 400 },
		{ /* 9 */ "resource limit", 500 },
		{ "??? CONMIN bug", 510 }
		};

	if (status < 0 || status > 9) status = 9;
	solve_result_num = solinfo[status].code;
	if (wantver)
		i = Sprintf(hbuf, "CONMIN %s, ASL(%ld):\n", conmin_version,
			ASLdate_ASL);
	else
		i = Sprintf(hbuf, "%s: ", conmin_bsname);
		i += Sprintf(hbuf+i, "%s (%d %siterations, %d evaluations)",
			solinfo[status].msg, kp->iter,
			kp->quadratic ? "QP " : "", neval);
	if (status < 3) {
		g_fmtop(buf, kp->primal_obj);
		i += Sprintf(hbuf+i, "\n objective %s", buf);
		}
	write_sol(hbuf, kp->x, NULL, &Oinfo);
	}

 static double Times[4];

 static void
show_times(VOID)
{
	int i, j, nn[2];
	FILE *f;
	real t;
	static char *what[2] = { "function", "gradient" };

	time_flag = 1;
	Times[3] = xectim_();
	for(i = 1; i <= 2; i++)
	    if (time_flag & i) {
		f = i == 1 ? stdout : stderr;
		fprintf(f,
		"\nTimes (seconds):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);

/*
		fprintf(f, "\nnn nonlinear evaluations = t seconds:\n");
		nn[0] = neval;
		nn[1] = ngeval;
		t = 0;
		for(j = 0; j < 2; j++) {
			fprintf(f, "\t%4d %10ss = %g\n",
				nn[j], what[j], ftimes[j]);
			t += ftimes[j];
			}
		fprintf(f, "Total nonlinear evaluation seconds = %g\n", t);
*/
		}
	}

 static void
#ifdef KR_headers
namecpy(s, t) char *s; char *t;
#else
namecpy(char *s, char *t)
#endif
{
	if (t) {
		strncpy(s,t,10);
		s[10] = 0;
		}
	else
		*s = 0;
	}

 void
#ifdef KR_headers
update_kp(kp) CONMIN *kp;
#else
update_kp(CONMIN *kp)
#endif
{
	kp->inftol = opns.inftol;
	kp->verbose = opns.v;
	kp->max = opns.maxflag ? opns.maxflag : 1;
	kp->itnlim = opns.itnlim;
	kp->timlim = opns.timlim;
	kp->pertlim = opns.pertlim;
	}

#ifdef __cplusplus
}
#endif

 int
#ifdef KR_headers
main(argc, argv) char **argv;
#else
main(int argc, char **argv)
#endif
{
	ASL	*asl;
	int	namelen, status;
	char	*av[2], **av0, *stub;
	char	fname[128];	/* solution file name */
	CONMIN	*kp;
	FILE 	*f;

	asl = ASL_alloc(ASL_read_pfgh);
	kp = openlp();
	if (kp == NULL) {
		fprintf(stderr,"Bug: openlp failure!\n");
		return 1;
	}

	Times[0] = xectim_();
	opns.timlim = HUGE_VAL;
	av0 = argv;
	stub = getstub(&argv, &Oinfo);
	if (wantver) {
		printf("CONMIN %s\n", conmin_version);
		need_nl = 0;
	}
	if (getopts(argv, &Oinfo)) return 1;
	if (amplin(asl,av0,stub,kp)) return 1;
	update_kp(kp);
	kp->h_init = nlinit;
	kp->h_update = nlupdate;
	kp->h_close = nlclose;
	Times[1] = xectim_();
	status = solvelp(kp);
	Times[2] = xectim_();
	amplout(asl, kp, status);
	show_times();
	closelp(kp);
	return 0;
}

#define asl cur_ASL

static int *
#ifdef KR_headers
ficopy(n, f) int n; fint *f;
#else
ficopy(int n, fint *f)
#endif
{
	int *i, *ie, *rv;

	if (sizeof(int) == sizeof(fint))
		rv = (int*)f;
	else {	/* a good compiler will eliminate this code in most cases... */
		i = rv = (int*)Malloc(n*sizeof(int));
		ie = i + n;
		while(i < ie)
			*i++ = *f++;
		}
	return rv;
	}

 static void
#ifdef KR_headers
nlinit(lp) Char *lp;
#else
nlinit(Char *lp)
#endif
{ Not_Used(lp); }

 static void
#ifdef KR_headers
nlclose(lp) Char *lp;
#else
nlclose(Char *lp)
#endif
{ Not_Used(lp); }

 static int
#ifdef KR_headers
nlupdate(vlp) Char *vlp; 
#else
nlupdate(Char *vlp)
#endif
{
	CONMIN *lp;
	double f;
	double *c, *x;
	fint nerror = 0;
	real ft[3];

	lp = (CONMIN*)vlp;
	c = lp->c;
	x = lp->x;

	/*----------------------------------------------------+
	| Nonlinear objective stuff                          */

	ft[0] = xectim_();
	neval++;
	f = objval(0, x, &nerror);
	lp->primal_obj = f;
	lp->f = f;
			
	if (nerror) return 1;
					
	ft[1] = xectim_();
	ngeval++;
	objgrd(0, x, c, 0);
	ft[2] = xectim_();

	ftimes[0] += ft[1] - ft[0];
	ftimes[1] += ft[2] - ft[1];
	return 0;
}

