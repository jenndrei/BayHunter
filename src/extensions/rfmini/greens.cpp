/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

# include <math.h>
# include <vector>
# include "synrf.h"
# include "cmat2.h"
# include "model.h"
# include "wave.h"

static inline Cmat2 exe (Cmat2 &e, Cmat2 &x);

#define WTYPE_PSV 0
#define WTYPE_SH  1

void
coeffm (double  u,
        double  vp1,  double  vs1,  double  rho1,
        double  vp2,  double  vs2,  double  rho2,
        Cmat2   &rd,  Cmat2   &td,  Cmat2   &ru,  Cmat2   &tu,
        Complex &rhd, Complex &thd, Complex &rhu, Complex &thu)
{
	double	mue1 = rho1*vs1*vs1, mue2 = rho2*vs2*vs2,
		c=2.*(mue1-mue2), u2=u*u, cu2=c*u2, t1, t2, t3;
	Complex rpp, rps, rsp, rss, tpp, tps, tsp, tss,
		d1, d2, t4, t5, t7,
		a1 = conj(sqrt(Complex(1./(vp1*vp1)-u2))),
		a2 = conj(sqrt(Complex(1./(vp2*vp2)-u2))),
		b1 = conj(sqrt(Complex(1./(vs1*vs1)-u2))),
		b2 = conj(sqrt(Complex(1./(vs2*vs2)-u2)));

	// (table 1) coefficients for an incident wave travelling
	//           in the medium with index 1 (downward)
	t1  = cu2 - rho1 + rho2;
	t2  = cu2 - rho1;
	t3  = cu2 + rho2;
	t4  = t3*a1 - t2*a2;

	d1  = t1*t1*u2 + t2*t2*a2*b2 + rho1*rho2*a2*b1;
	d2  = c*c*u2*a1*a2*b1*b2 + t3*t3*a1*b1 + rho1*rho2*a1*b2;
	t5  = 1./(d1+d2);
	t7  = 2.*rho1*t5;

	rpp = (d2-d1)*t5;
	rps =-2.*u*a1*t5*(t1*t3+c*t2*a2*b2);
	tpp = a1*t7*(t3*b1-t2*b2);
	tps =-a1*t7*u*(t1+c*a2*b1);
	rss = (d2-d1-2.*rho1*rho2*(a1*b2-a2*b1))*t5;
	rsp = 2.*u*b1*t5*(t1*t3+c*t2*a2*b2);
	tss = b1*t7*t4;
	tsp = b1*t7*u*(t1+c*a1*b2);

	rd = Cmat2(rpp, rsp, rps, rss);
	td = Cmat2(tpp, tsp, tps, tss);

	// (table 2) coefficients for an incident wave travelling
	//           in the medium with index 2 (upward)
	d1  = t1*t1*u2 + t3*t3*a1*b1 + rho1*rho2*a1*b2;
	d2  = c*c*u2*a1*a2*b1*b2 + t2*t2*a2*b2 + rho1*rho2*a2*b1;
	t5  = 1./(d1+d2);
	t7  = 2.*rho2*t5;

	rpp = (d2-d1)*t5;
	rps = 2.*u*a2*t5*(t1*t2+c*t3*a1*b1);
	tpp = a2*t7*(t3*b1-t2*b2);
	tps =-a2*t7*u*(t1+c*a1*b2);
	rss = (d2-d1-2.*rho1*rho2*(a2*b1-a1*b2))*t5;
	rsp = -2.*u*b2*t5*(t1*t2+c*t3*a1*b1);
	tss = b2*t7*t4;
	tsp = b2*t7*u*(t1+c*a2*b1);

	ru = Cmat2(rpp, rsp, rps, rss);
	tu = Cmat2(tpp, tsp, tps, tss);

	// coefficients for the sh case
	Complex mb1 = mue1*b1, mb2 = mue2*b2, mmm = 1./(mb1+mb2);

	rhd = (mb1-mb2)*mmm;
	rhu = -rhd;
	thd = 2.*mb1*mmm;
	thu = 2.*mb2*mmm;
}

void
coeffs(double u, double vp, double vs, Cmat2 &ru, Complex &rhu)
{
	// computes the free-surface reflection coefficients
	double u2 = u*u;
	Complex rpp, rps, rsp, rss, a, b, t1, t2, t3, d, d1, d2;

	a = sqrt(Complex(1./(vp*vp)-u2));
	b = sqrt(Complex(1./(vs*vs)-u2));

	t1  = 2.*vs*vs;
	t2  = t1*u2-1.;
	d1  = t2*t2;
	d2  = t1*t1*u2*a*b;
	d   = d1+d2;
	t3  = 2.*t1*u*t2/d;
	rpp = (d2-d1)/d;
	rsp = -b*t3;
	rps = a*t3;
	rss = rpp;

	ru = Cmat2(rpp, rsp, rps, rss);

	// total reflexion of SH wave at the surface
	rhu = Complex (1., 0.);
}

static void
coeff (	double     slowness,	    // slowness
	FlatLayer  *upp,	    // upper layer
	FlatLayer  *low,	    // lower layer
	Cmat2   &rd,  Cmat2   &td,  Cmat2   &ru,  Cmat2   &tu,
	Complex &rhd, Complex &thd, Complex &rhu, Complex &thu)
{
	if (upp==NULL) {
		// compute free-surface coefficients
		coeffs(slowness, low->vp, low->vs, ru, rhu);
		 rd =  td =  tu = 0.;
		rhd = thd = thu = 0.;
	}
	else	coeffm(slowness,
			upp->vp, upp->vs, upp->rh,
			low->vp, low->vs, low->rh,
			rd,  td,  ru,  tu,
			rhd, thd, rhu, thu);
}

void ccfork (int n, Complex *x, int signi);

static void
iftr (int nsamp, Complex *cf, double *f)
{
	int i;
	Complex cx[nsamp];

	// transform the complex spectrum to time domain
	// by invoking the FFT routine 'ccfork'
	// the resultung time series is real

	// scale factor for fft
	double q = 1./sqrt((double)nsamp);

	for (i=0; i<nsamp/2+1; i++)
		cx[i] = cf[i];
	for (i=nsamp/2+1; i<nsamp; i++)
		cx[i] = conj (cx[nsamp-i]);

	ccfork(nsamp, cx, 1);

	for (i=0; i<nsamp; i++)
		f[i] = q*real(cx[i]);
}


static void
iftr2 (int nsamp, Complex *cf1, Complex *cf2, double *f1, double *f2)
{
	int i;
	Complex cx[nsamp], cx1[nsamp], cx2[nsamp], imi(0.,1.);

	// simultaneous transformation of two complex spectra
	// to time domain
	// the resultung time series is real

	// scale factor for fft
	double q = 1./sqrt((double)nsamp);

	// transform the r and z component *simultaneously*
	for (i=0; i<nsamp/2+1; i++) {
		cx1[i] = cf1[i];
		cx2[i] = cf2[i];
	}
	for (i=nsamp/2+1; i<nsamp; i++) {
		cx1[i] = conj (cx1[nsamp-i]);
		cx2[i] = conj (cx2[nsamp-i]);
	}
	for (i=0; i<nsamp; i++)
		cx[i] = cx1[i] + imi*cx2[i];

	ccfork(nsamp, cx, 1);

	for (i=0; i<nsamp; i++)	{
		f1[i] = q*real(cx[i]);
		f2[i] = q*imag(cx[i]);
	}

	return;
}

static void
top_down (int nlay,
	  Cmat2  e[], Cmat2 nt[], Cmat2 nb[],
	  Cmat2 ru[], Cmat2 rd[], Cmat2 tu[], Cmat2 td[],
	  Cmat2  g[], int options=0)
{
	// top down approach like in Mueller (1985)

	Cmat2 q(0.);
	//const Cmat2 I (1.,0.,0.,1.); // 2x2 unity matrix
	const Cmat2 I(1.); // 2x2 unity matrix

	for (int i=1; i<nlay; i++) {
		if (i==1) nt[i] = ru[1];  // Mueller eq. (44)
		else      nt[i] = ru[i] + td[i]*nb[i-1]*q;
  
		if (options & SUPPRESS_MULTIPLES) {
			nt[i] -= ru[i];
			q = tu[i+1];
		}
		else { // normal case
			nb[i] = exe(e[i],nt[i]);  // nb[i] = e[i]*nt[i]*e[i];
			q = inv(I-rd[i+1]*nb[i])*tu[i+1];
		}

		if (i==1) g[i] =        e[1]*q;
		else      g[i] = g[i-1]*e[i]*q;
	}
}

// This is the corresponding routine for an SH wave
static void
top_down (int nlay,
	  Complex  e[], Complex nt[], Complex nb[],
	  Complex ru[], Complex rd[], Complex tu[], Complex td[],
	  Complex  g[], int options=0)
{
	Complex q;

	for (int i=1; i<nlay; i++) {
		if (i==1) nt[i] = ru[1];
		else      nt[i] = ru[i] + td[i]*nb[i-1]*q;
  
		if (options & SUPPRESS_MULTIPLES) {
			nt[i] -= ru[i];
			q      = tu[i+1];
		}
		else { // normal case
			nb[i] = e[i]*nt[i]*e[i];
			q     = 1./(1.-rd[i+1]*nb[i])*tu[i+1];
		}

		if (i==1) g[i] =        e[1]*q;
		else      g[i] = g[i-1]*e[i]*q;
	}
}

static void
bottom_up (int nlay,
	   Cmat2  e[], Cmat2 mt[], Cmat2 mb[],
	   Cmat2 ru[], Cmat2 rd[], Cmat2 tu[], Cmat2 td[],
	   Cmat2  f[], int options=0)
{
	// bottom up approach, not described in Mueller (1985)

	Cmat2 q;
	const Cmat2 I(1.); // 2x2 unity matrix

	for (int i=nlay-1; i>=1; i--) {
		if (i==(nlay-1))
		// initial condition of the recursion: MT[n]   = 0
		//                                <=>  MB[n-1] = Rd[n]
		// i.e. at the lowermost layer boundary the reflectivity
		// is the reflection coefficient rd
			mb[i] = rd[nlay];
		else	mb[i] = rd[i+1] + tu[i+1]*q*mt[i+1]*td[i+1];
      
		// mt[i] = e[i]*mb[i]*e[i];
		mt[i] = exe (e[i],mb[i]);
		q     = inv(I-mt[i]*ru[i]);
      
		if (i==(nlay-1)) f[i] = q*e[i]*tu[i+1];
		else             f[i] = q*e[i]*tu[i+1]*f[i+1];
	}
}

static void
bottom_up (int nlay,
	   Complex  e[], Complex mt[], Complex mb[],
	   Complex ru[], Complex rd[], Complex tu[], Complex td[],
	   Complex  f[], int options=0)
{
	Complex q;

	for (int i=nlay-1; i>=1; i--) {
		if (i==(nlay-1))
		// initial condition of the recursion: MT[n]   = 0
		//                                <=>  MB[n-1] = Rd[n]
		// i.e. at the lowermost layer boundary the reflectivity
		// is the reflection coefficient rd
			mb[i] = rd[nlay];
		else	mb[i] = rd[i+1] + tu[i+1]*q*mt[i+1]*td[i+1];
      
		mt[i] = e[i]*mb[i]*e[i];
		q     = 1./(1.-mt[i]*ru[i]);
      
		if (i==(nlay-1)) f[i] = q*e[i]*tu[i+1];
		else             f[i] = q*e[i]*tu[i+1]*f[i+1];
	}
}

static void
displacement_matrix (double p, double vp, double vs, Cmat2 &m)
// computation of the matrix h for the free surface displacements
// after: Mueller (1985)
{
        double  vp2 = vp*vp, vs2 = vs*vs, p2 = p*p, x = 1.-2.*vs2*p2;
        Complex // ****** eq. (89) part 1 *************
                a1 = conj(sqrt(Complex(1./vp2-p2))),
                b1 = conj(sqrt(Complex(1./vs2-p2))),
                q  = 1./(x*x + 4.*vs2*vs2 *p2*a1*b1);

        m.c11 =  q*a1* b1*2.*vs2*p;
        m.c12 =  q*b1*(1.-2.*vs2*p2);
        m.c21 =  q*a1*(1.-2.*vs2*p2);
        m.c22 = -q*a1* b1*2.*vs2*p;
}

static void 
decomp (int n, Complex *cz, Complex *cr, double p, double vp, double vs)
{
	/* XXX WAS IST MIT p>1/vp ??? */
	double	a   = sqrt(1./(vp*vp)-p*p),
		b   = sqrt(1./(vs*vs)-p*p),
		m11 = -(2*vs*vs*p*p-1.)/(vp*a),
		m12 =   2.*p*vs*vs/vp,
		m21 =  -2.*p*vs,
		m22 =  (1.-2.*vs*vs*p*p)/(vs*b);

	for (int i=0; i<n; i++) {
		Complex cx(cz[i]*m11 + cr[i]*m12),
			cy(cz[i]*m21 + cr[i]*m22);
		cz[i] = cx;
		cr[i] = cy;
	}
}

static void
compute_rf (int wave_type, Complex *cr, Complex *cz,
       int nsamp, double fsamp, double water, double tshift, double a,
       double p, double vp0, double vs0, Complex *crf)
//     computes a receiver function in the frequency domain
//     from the radial and vertical components 'cr' and 'cz'
//
//     'nsamp' and 'fsamp' are the number of samples and
//         sampling frequency, and refer to the corresponding
//         **time domain** functions.
//     'water' is the waterlevel for the deconvolution.
//     the origin of the time axis is shifted by 'tshift'.
//     'a' is the constant for the gauss filter

//     changed for index range 0...nfreq-1 31.08.97
//     added parameter 'angle' 07.02.99
{
	double	zmax=0.,w,wa,dw = 2.0*M_PI*fsamp/nsamp,
		denom,q = sqrt(M_PI)*fsamp/a;
	int j,nfreq = nsamp/2+1;
	Complex cq;

	if (vs0 > 0.01 && fabs(p)>0.0001)
		// decompose Z/R -> P/SV (somewhat redundant here...)
		decomp (nfreq, cz, cr, p, vp0, vs0);

	if (wave_type==SV_Wave) {
		// Deconvolve P with SV, unlike P-RF's where we 
		// deconvolve SV with P. Thus, swap data pointers.
		Complex *tmp = cz; cz = cr; cr = tmp;
	}

	// XXX The waterlevel stabilization is actually obsolete for (noise free) synthetics.

	for (j=0; j<nfreq; j++)
		if (abs(cz[j])>zmax) 
			zmax = abs(cz[j]);

	for (j=0; j<nfreq; j++) {
		w = dw*j;
		denom = real(cz[j]*conj(cz[j]));
//		if (denom<wlev) denom = wlev;
		crf[j] = cr[j]*conj(cz[j])/denom;

		// application of the gauss filter and shifting
		// of the origin of the time axis by 'tshift'
		wa = w/a;
		wa = (wa>50.0)?50.0:wa; // avoid overflow at high frequencies
		cq     = q*exp(Complex(-0.25*(wa*wa),-w*tshift));
		crf[j] = crf[j]*cq;
		cr[j]  =  cr[j]*cq;
		cz[j]  =  cz[j]*cq;
	}

	return;
}

static void
calcresp_core (
	FlatModel model,                                        // model
	Wave      wave,
	double    fref,        // reference freq. for q.f.      // model
	int       nsamp,       // number of samples             // trace
	double    fsamp,       // sampling frequency            // trace
	Complex  *cz,	 Complex  *cr,	  Complex  *ct,
	Complex **czp=0, Complex **crp=0, Complex **ctp=0,
	// drdp==NULL indicates that no pd matrix
	// needs to be computed (previously: pdflag)
	int       laymin=0,    // first perturbed layer         // model
	int       laymax=0,    // last  perturbed layer         // model
	double    pert=0.001,  // perturbation
	int       options=0)   // option flags
{
	int	nlay = model.numLayers(),
		j, k=0, nfreq=nsamp/2 + 1,
		// if crp and czp are not null pointers, we
		// want to compute the partial derivatives
		compute_partial_derivatives = (crp && czp) ? -1 : 0;
	double	p = wave.slowness, p2 = p*p, dw, wref;
	Complex ii(0.,1.), dummy;
	Cmat2	// reflection and transmission coefficients for each interface
		ru[nlay+1], rd[nlay+1], tu[nlay+1], td[nlay+1],
		// reflectivities
		mb[nlay+1], mt[nlay+1], nb[nlay+1], nt[nlay+1],
		// transmissivities
		f[nlay+1], // transmissivity from the bottom up approach
		g[nlay+1], // transmissivity from the top down  approach
		e[nlay+1], // Phase matrix e
			   // only e11 and e22 are used, the others are zero
		// Perturbed reflection and transmission
		// coefficients connected to the UPPER boundary of
		// the perturbed layer:
		rup1[nlay+1], rdp1[nlay+1], tup1[nlay+1], tdp1[nlay+1],
		// ... and those connected to its LOWER boundary:
		rup2[nlay+1], rdp2[nlay+1], tup2[nlay+1], tdp2[nlay+1];

	Complex rhu[nlay+1], rhd[nlay+1], thu[nlay+1], thd[nlay+1];
	Complex rhup1[nlay+1], rhdp1[nlay+1], thup1[nlay+1], thdp1[nlay+1];
	Complex rhup2[nlay+1], rhdp2[nlay+1], thup2[nlay+1], thdp2[nlay+1];
	Complex eh[nlay+1]; // Phase matrix e
  
	Cmat2 t, h, hp, q;
	const Cmat2 I(1); // 2x2 unity matrix

	wref = 2.*M_PI*fref;

	if (compute_partial_derivatives) {
		// if we want the responses / receiver function only,
		// we can compute them for arbitrary vs and rho
		// but if we want to compute the partial derivative matrix,
		// the vs and rho of the input model must be set into a
		// relationship to vp. This relationship is kept during
		// the whole inversion.
		for (int i=1; i<=nlay; i++)
			model.layer(i).perturb(0);
	}

	// complex reflexion and transmission matrices and coefficients
	// for displacement for all boundaries
	for (int i=1; i<=nlay; i++) {
		// we compute the SH coefficients even if we don't need them
		// (doesn't cost us much CPU...)
		coeff (p, i==1 ? 0 : &model.layer(i-1), &model.layer(i),
		       rd[i],  td[i],  ru[i],  tu[i],
		       rhd[i], thd[i], rhu[i], thu[i]);
	}

	if (compute_partial_derivatives) {
		for (int i=1; i<=nlay; i++) {
			FlatLayer currentLayer = model.layer(i);

			currentLayer.perturb (pert);

			// reflection and transmission coefficients for
			// the UPPER boundary of the perturbed layer
			coeff(p, i==1 ? 0 : &model.layer(i-1), &currentLayer,
				rdp1[i],  tdp1[i],  rup1[i],  tup1[i],
				rhdp1[i], thdp1[i], rhup1[i], thup1[i]);

			// The lower halfspace has no lower
			// boundary, so we are done.
			if (i==nlay) break;

			// reflection and transmission coefficients for
			// the LOWER boundary of the perturbed layer
			coeff(p, &currentLayer, &model.layer(i+1),
				rdp2[i],  tdp2[i],  rup2[i],  tup2[i],
				rhdp2[i], thdp2[i], rhup2[i], thup2[i]);
		} // for
	} // if (compute_partial_derivatives)

	// compute the matrix h ... (for P/SV)
	displacement_matrix (p, model.layer(1).vp, model.layer(1).vs, h);

	// ... and the corresponding perturbed matrix hp (for P/SV)
	if (compute_partial_derivatives) {
		FlatLayer lay = model.layer(1);
		lay.perturb (pert);
		displacement_matrix (p, lay.vp, lay.vs, hp);
	}

	// *********** START OF FREQUENCY LOOP ***********

	// omega increment for fft
	dw = 2.0*M_PI*fsamp/nsamp;

	// determine the time the direct wave travels through the medium
	double t0 = 0.;
        switch (wave.type) {
        case P_Wave:   // for an incident P wave
		for (int i=1; i<=nlay; i++) {
			double vp = model.layer(i).vp,
				d = model.layer(i).h;
			t0 += d*sqrt(1./(vp*vp)-p2);
		}
                break;
        case SV_Wave:  // for an incident SV wave
		for (int i=1; i<=nlay; i++) {
			double vs = model.layer(i).vs,
				d = model.layer(i).h;
			t0 += d*sqrt(1./(vs*vs)-p2);
		}
                break;
        }

	for (j=0; j<nfreq; j++) {
		double	w   = dw*j,                 // angular frequency w
			lgw = j ? log(w/wref) : 0;  // =0 for w=0

		// phase matrix e[i] (eq. (23))
		for (int i=1; i<=nlay; i++) {
			double z, d, vp, vs, rho, qp, qs;
			model.layer(i).getAll(&z,&d,&vp,&vs,&rho,&qp,&qs);
			Complex	miwd (0., -w*d), // miiwd = -ii*w*d,
				// frequency dependent, complex velocities
				// after Mueller (1985) eq. (132)
				vpc = vp*(1.+lgw/(M_PI*qp)+ii/(2.*qp)),
				vsc = vs*(1.+lgw/(M_PI*qs)+ii/(2.*qs)),
				// vertical slownesses
				plc = sqrt (1./(vpc*vpc)-p2),
				slc = sqrt (1./(vsc*vsc)-p2);

			// We only need to compute the diagonal elements
			// e11 and e22 because e12 and e21 are always zero.
			e[i] = Cmat2 (exp(miwd*plc), 0, 0, exp(miwd*slc));
			eh[i] = exp(miwd*slc);
		}

		// beginning of the recursion

		{ // SH
		    std::complex<double>	ru_[nlay+1], rd_[nlay+1],
		    tu_[nlay+1], td_[nlay+1],
		    /*	mb[nlay+1], mt[nlay+1], */  nb_[nlay+1], nt_[nlay+1],
		    /*	f[nlay+1],  */ g_[nlay+1],  e_[nlay+1];
		top_down (nlay, e_, nt_, nb_, ru_, rd_, tu_, td_, g_, options);
		ct[j] = 2.* /* 2* */ g_[nlay-1];
		} // XXX XXX XXX 

		// computation of the unperturbed transmissivities
		// EITHER by the bottom-up OR top-down method

		// top down approach like in Mueller (1985)
		top_down (nlay, e, nt, nb, ru, rd, tu, td, g, options);

		if (compute_partial_derivatives)
			bottom_up (nlay, e, mt, mb, ru, rd, tu, td, f);

		// t = 2*h*f[1];    // bottom-up
		t = 2*h*g[nlay-1];  // top-down

		// copy responses for the desired wave type
		switch (wave.type) {
		case P_Wave:   // for an incident P wave
			t.get (cr+j, &dummy, cz+j, &dummy);
			break;
		case SV_Wave:  // for an incident SV wave
			t.get (&dummy, cr+j, &dummy, cz+j);
			break;
		}
		Complex qq = exp(Complex(0., w*t0));
		cr[j] *= qq;
		cz[j] *= qq;

		// if only unperturbed responses / receiver functions
		// are desired, the computation of the matrix drdp can
		// be skipped -> next frequency
		if (compute_partial_derivatives==0) continue;

		// **********************************************************
		// ** Beginning of the computation of perturbed quantities **
		// **********************************************************

		// perturb layers from laymin to laymax
		for (k=laymin; k<=laymax; k++) {
			Cmat2	mbp, ntp, // perturbed reflectivities
				ep,       // perturbed phase matrix
				x,        // matrix x for the indirect approach
				fp, gp;

			if      (k==1) t = 2*hp; // use perturbed version of h
			else if (k==2) t = 2*h;
			else           t = 2*h*g[k-2];

			if ( ! model.layer(k).isLowerHalfspace()) {
				// compute perturbed phase matrix ep[k]
				double z,d,vp,vs,rho,qp,qs;
				FlatLayer currentLayer = model.layer(k);
				currentLayer.perturb (pert);
				currentLayer.getAll(&z,&d,&vp,&vs,&rho,&qp,&qs);
				Complex	miwd (0., -w*d), // miiwd = -ii*w*d,
					// frequency dependent,
					// complex velocities after
					// Mueller (1985) eq. (132)
					vpc = vp*(1.+lgw/(M_PI*qp)+ii/(2.*qp)),
					vsc = vs*(1.+lgw/(M_PI*qs)+ii/(2.*qs)),
					// vertical slownesses
					plc = sqrt (1./(vpc*vpc)-p2),
					slc = sqrt (1./(vsc*vsc)-p2);

				// We only need to compute the
				// diagonal elements e11 and e22
				// because e12 and e21 are always
				// zero.

				ep = Cmat2 (exp(miwd*plc), 0, 0,
					    exp(miwd*slc));
			}

			// compute perturbed gp(k-1) and ntp(k)
			if (k==1) // uppermost layer: nt=ru[1]
				ntp = rup1[1];
			else {
				q  = inv(I-rdp1[k]*nb[k-1])*tup1[k];
				gp = e[k-1]*q;

				if (k<nlay) // perturbed ntp(k)
					ntp = rup1[k] + tdp1[k]*nb[k-1]*q;

				// update transmissivity t = t*gp(k-1)
				t *= gp;
			} // perturbed gp(k-1) and ntp(k) ready

			// For k=nlay all of the following computations
			// of the loop can me omitted.
			if (k<nlay) {
				// compute perturbed fp(k+1) and mbp(k)
				if (k==(nlay-1))
					mbp = rdp2[k];
				// f(k+1) is not needed for k=nlay-1
				else {
					q   = inv(I-mt[k+1]*rup2[k]);
					fp  = q*e[k+1]*tu[k+2];
					mbp = rdp2[k] +
					      tup2[k]*q*mt[k+1]*tdp2[k];
				} // perturbed fp(k+1) and mbp(k) ready

				// x = inv(I-ep*mbp*ep*ntp)*ep*tup2[k];
				x = inv(I-exe(ep,mbp)*ntp)*ep*tup2[k];

				if (k==(nlay-1))
					t *= x;      // tt=tt*x(k)
				else    t *= x*fp;   // tt=tt*x(k)*fp(k+1)
				if (k<(nlay-2))
					t *= f[k+2]; // tt=tt*f[k+2]
			}

			switch (wave.type) {
			case P_Wave:   // for an incident P wave
				t.get (crp[k]+j, &dummy, czp[k]+j, &dummy);
				break;
			case SV_Wave:  // for an incident SV wave
				t.get (&dummy, crp[k]+j, &dummy, czp[k]+j);
				break;
			}
			crp[k][j] *= qq;
			czp[k][j] *= qq;
		}
	}
	// **** end of frequency loop *****
}

void
calcresp (
	FlatModel model,                                        // model
	Wave      wave,
	double    fref,        // reference freq. for q.f.      // model
	int       nsamp,       // number of samples             // trace
	double    fsamp,       // sampling frequency            // trace
	double    water,       // waterlevel                    // RF
	double    tshift,      // time shift                    // RF
	double    a,           // Gauss parameter               // RF
	double	  vp_top,
	double	  vs_top,
	double    *zz,         // << vertical responses         // trace
	double    *rr,         // << radial responses           // trace
	double    *rf,         // << receiver function          // RF trace
	int       first,       // XXX first sample used in inversion
	int       nused,       // number of used samples 
	double    **drdp,      // << partial derivative matrix
	// drdp==NULL indicates that no pd matrix
	// needs to be computed (previously: pdflag)
	int       laymin,      // first perturbed layer         // model
	int       laymax,      // last  perturbed layer         // model
	InvParam &ipar, // XXX
	double    pert,  // perturbation
	int       options)   // option flags
{
	int	nlay  = model.numLayers(),    
		nfreq = nsamp/2 + 1;
	Complex cz[nfreq], cr[nfreq], ct[nfreq];
	Complex	**czp=0, **crp=0, **ctp=0;
//	Complex czp[nlay+1][nfreq], crp[nlay+1][nfreq], ctp[nlay+1][nfreq];
	double	*tt=rf;
tt=NULL; /* XXX */

	if (drdp) {
		czp  = new Complex* [nlay+1];
		crp  = new Complex* [nlay+1];
		ctp  = new Complex* [nlay+1];

		for (int k=laymin; k<=laymax; k++) {
			czp[k] = new Complex [nfreq];
			crp[k] = new Complex [nfreq];
			ctp[k] = new Complex [nfreq];
		}
	}

	if (drdp) calcresp_core (model, wave, fref, nsamp, fsamp,
			     cz, cr, ct, czp, crp, ctp,
			     laymin, laymax, pert);
	else      calcresp_core (model, wave, fref, nsamp, fsamp,
			     cz, cr, ct);


	{	Complex crf [nfreq];

		// compute and transform the unperturbed receiver
		// function as well as r and z components
		compute_rf (wave.type, cr, cz, nsamp, fsamp, water, tshift, a,
			wave.slowness, vp_top, vs_top, crf);

		// transform the unperturbed receiver function to
		// time domain
		iftr (nsamp, crf, rf);

		if (rr != NULL  &&  zz != NULL)
			// transform the unperturbed r and z
			// components simultaneously
			iftr2 (nsamp, cr, cz, rr, zz);

		if (tt != NULL)
			iftr (nsamp, ct, tt);
	}

	// compute and transform the PERTURBED receiver functions
	// from the transmission responses saved in crp and czp

	if (drdp) {
		double  ipert = 1./pert, rfp1[nsamp], rfp2[nsamp];
		Complex	cf1[nfreq], cf2[nfreq];

		for (int k=laymin; k<=laymax; k++) {
			int i;

			if (k&1) { // for odd layer indices
				compute_rf (wave.type,crp[k], czp[k], nsamp, fsamp, water,
				       tshift, a, wave.slowness,
					vp_top, vs_top, cf1);

				if (k==laymax) {
					// last layer in case that nlay is odd
					iftr (nsamp,cf1,rfp1);
					for (int j=1; j<=nused; j++) {
						i = j+first-1;
						drdp[j][k] =
						    ipert*(rfp1[i]-rf[i]);
					}
				}
			}
			else { // for even layer indices
				compute_rf (wave.type,crp[k], czp[k], nsamp, fsamp, water,
				       tshift, a, wave.slowness,
					vp_top, vs_top, cf2);

				if (k==laymin) {
				// first layer in case that laymin is even
					iftr (nsamp, cf2, rfp2);

					for (int j=1; j<=nused; j++) {
						i = j+first-1;
						drdp[j][k] =
						    ipert*(rfp2[i]-rf[i]);
					}
				}
				else {
					// transform the receiver functions
					// for perurbed layers k and k-1
					// simultaneously
					iftr2 (nsamp, cf1, cf2, rfp1, rfp2);

					// compute the partial derivation
					// for k and k-1
					for (int j=1; j<=nused; j++) {
						i = j+first-2;
						drdp[j][k-1] =
						    ipert*(rfp1[i]-rf[i]);
						drdp[j][k] =
						    ipert*(rfp2[i]-rf[i]);
					}
				}
			}
		} // for k

		for (int k=laymin; k<=laymax; k++) {
			delete [] crp[k]; 
			delete [] czp[k];
			delete [] ctp[k];
		}
		delete [] crp; 
		delete [] czp;
		delete [] ctp;
	}
}


static inline Cmat2
exe (Cmat2 &e, Cmat2 &x)
{
	// Computes e*x*e if e12 and e21 are zero (not checked!).
	// This way many superfluous complex multiplications are
	// avoided and the computation becomes much faster!
	Complex	x11, x12, x21, x22, e11, e12, e21, e22;

	e.get (&e11, &e12, &e21, &e22);
	x.get (&x11, &x12, &x21, &x22);

	e12  = e11*e22;
	e11 *= e11;
	e22 *= e22;

	return Cmat2 (x11*e11, x12*e12, x21*e12, x22*e22);
}

