/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

#include "defaults.h"
#include "model.h"
#include "wave.h"
#include "cmat2.h"
#include "rf.h"

#define VPVS(sigma)  sqrt((1.-(sigma))/(.5-(sigma)))

char *progname;

void synrf(Wave &wave,
	int nlay,
        double *z,  double *vp, double *vs,
        double *rh, double *qp, double *qs,
        RfnParam &par, double *syn_z, double *syn_r, double *syn_rf);

extern "C"
{

int synrf_cwrap(
	int nsamp, double fsamp, double tshift, double p, double a,
	double nsv, double sigma, int waveno, int nlay,
	double  *z, double *vp, double *vs,
	double *rh, double *qp, double *qs,
	double *fz, double *fr, double *rf);

int coeff_cwrap(double p,
	double vp1, double vs1, double rh1,
	double vp2, double vs2, double rh2,
	int displacement,
	double *rd11r, double *rd12r, double *rd21r, double *rd22r,
	double *rd11i, double *rd12i, double *rd21i, double *rd22i,
	double *td11r, double *td12r, double *td21r, double *td22r, 
	double *td11i, double *td12i, double *td21i, double *td22i, 
	double *ru11r, double *ru12r, double *ru21r, double *ru22r, 
	double *ru11i, double *ru12i, double *ru21i, double *ru22i, 
	double *tu11r, double *tu12r, double *tu21r, double *tu22r, 
	double *tu11i, double *tu12i, double *tu21i, double *tu22i, 
	double *rhdr,  double *thdr,  double *rhur,  double *thur,
	double *rhdi,  double *thdi,  double *rhui,  double *thui);

int coeffs_cwrap(double p,
	double vp, double vs, double rh,
	double *ru11r, double *ru12r, double *ru21r, double *ru22r, 
	double *ru11i, double *ru12i, double *ru21i, double *ru22i, 
	double *rhur,  double *rhui);
}

#define DEGREES_PER_KM 0.00899

int
synrf_cwrap(
	int nsamp, double fsamp, double tshift, double p, double a,
	double nsv, double sigma, int waveno, int nlay,
	double  *z, double *vp, double *vs,
	double *rh, double *qp, double *qs,
	double *fz, double *fr, double *rf)
{
	RfnParam par;

	par.p = p;
	par.a = a;
	par.nsamp = nsamp;
	par.fsamp = fsamp;
	par.tshift = tshift;

	par.vptop = nsv*VPVS(sigma);
	par.vstop = nsv;

	Wave wave = { waveno, par.p*DEGREES_PER_KM }; // MESSY!!!
	synrf(wave, nlay, z, vp, vs, rh, qp, qs, par, fz, fr, rf);

	return 1;
}

void
coeffm (double  u,
        double  vp1,  double  vs1,  double  rho1,
        double  vp2,  double  vs2,  double  rho2,
        Cmat2   &rd,  Cmat2   &td,  Cmat2   &ru,  Cmat2   &tu,
        Complex &rhd, Complex &thd, Complex &rhu, Complex &thu);


int
coeff_cwrap(double p,
	    double vp1, double vs1, double rh1,
	    double vp2, double vs2, double rh2,
	    int displacement,
	    double *rd11r, double *rd12r, double *rd21r, double *rd22r,
	    double *rd11i, double *rd12i, double *rd21i, double *rd22i,
	    double *td11r, double *td12r, double *td21r, double *td22r, 
	    double *td11i, double *td12i, double *td21i, double *td22i, 
	    double *ru11r, double *ru12r, double *ru21r, double *ru22r, 
	    double *ru11i, double *ru12i, double *ru21i, double *ru22i, 
	    double *tu11r, double *tu12r, double *tu21r, double *tu22r, 
	    double *tu11i, double *tu12i, double *tu21i, double *tu22i, 
	    double *rhdr,  double *thdr,  double *rhur,  double *thur,
	    double *rhdi,  double *thdi,  double *rhui,  double *thui)
{
	Cmat2   rd,  td,  ru,  tu;
	Complex rhd, thd, rhu, thu;

	coeffm(p, vp1,vs1,rh1, vp2,vs2,rh2, rd,td,ru,tu, rhd,thd,rhu,thu);

	if (displacement) {
		ru.c12 *= vs2/vp2;
		ru.c21 *= vp2/vs2;
		tu.c11 *= vp2/vp1;
		tu.c12 *= vs2/vp1;
		tu.c21 *= vp2/vs1;
		tu.c22 *= vs2/vs1;

		rd.c12 *= vs1/vp1;
		rd.c21 *= vp1/vs1;
		td.c11 *= vp1/vp2;
		td.c12 *= vs1/vp2;
		td.c21 *= vp1/vs2;
		td.c22 *= vs1/vs2;
	}

	*rd11r = real(rd.c11);  *rd11i = imag(rd.c11);
	*rd12r = real(rd.c12);  *rd12i = imag(rd.c12);
	*rd21r = real(rd.c21);  *rd21i = imag(rd.c21);
	*rd22r = real(rd.c22);  *rd22i = imag(rd.c22);

	*td11r = real(td.c11);  *td11i = imag(td.c11);
	*td12r = real(td.c12);  *td12i = imag(td.c12);
	*td21r = real(td.c21);  *td21i = imag(td.c21);
	*td22r = real(td.c22);  *td22i = imag(td.c22);

	*ru11r = real(ru.c11);  *ru11i = imag(ru.c11);
	*ru12r = real(ru.c12);  *ru12i = imag(ru.c12);
	*ru21r = real(ru.c21);  *ru21i = imag(ru.c21);
	*ru22r = real(ru.c22);  *ru22i = imag(ru.c22);

	*tu11r = real(tu.c11);  *tu11i = imag(tu.c11);
	*tu12r = real(tu.c12);  *tu12i = imag(tu.c12);
	*tu21r = real(tu.c21);  *tu21i = imag(tu.c21);
	*tu22r = real(tu.c22);  *tu22i = imag(tu.c22);

	*rhdr  = real(rhd);     *rhdi  = imag(rhd);
	*thdr  = real(thd);     *thdi  = imag(thd);
	*rhur  = real(rhu);     *rhui  = imag(rhu);
	*thur  = real(thu);     *thui  = imag(thu);

	return 0;
}

void coeffs(double u, double vp, double vs, Cmat2 &ru, Complex &rhu);

int coeffs_cwrap(double p,
	         double vp, double vs, double rh,
	         double *ru11r, double *ru12r, double *ru21r, double *ru22r,
	         double *ru11i, double *ru12i, double *ru21i, double *ru22i,
	         double *rhur,  double *rhui)
{
	Cmat2   ru;
	Complex rhu;

	coeffs(p, vp, vs, ru, rhu);

	*ru11r = real(ru.c11);  *ru11i = imag(ru.c11);
	*ru12r = real(ru.c12);  *ru12i = imag(ru.c12);
	*ru21r = real(ru.c21);  *ru21i = imag(ru.c21);
	*ru22r = real(ru.c22);  *ru22i = imag(ru.c22);

	*rhur  = real(rhu);     *rhui  = imag(rhu);
	return 0;
}
