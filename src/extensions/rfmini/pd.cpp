# include <stdio.h>
# include <stdlib.h>
# include <math.h>

static int
svbksb (double **u, double *w, double **v, int m, int n, double *b, double *x)
{
	int i, j, k;
	double s, tmp[n+1];

	for (j=1; j<=n; j++) {
		s = 0.;

		if (w[j] != 0.) {
			for (i=1; i<=m; i++)
				s += u[i][j]*b[i];
			s /= w[j];
		}

		tmp[j] = s;
	}

	for (j=1; j<=n; j++) {
		s = 0.;
		for (k=1; k<=n; k++)
			s += v[j][k]*tmp[k];
		x[j] = s;
	}

	return 0;
}


/* ARGHHHH... */
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))

static int
svdcmp (double **a, int m, int n, double *w, double **v)
{
	int flag,i,its,j,jj,k,l=0,nm=0;
	double anorm=0.,c,f,g=0.,h,s,scale=0.,x,y,z,rv1[n+1];

	for (i=1; i<=n; i++) {
		l = i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.;

		if (i <= m) {
			for (k=i; k<=m; k++)
				scale += fabs(a[k][i]);

			if (scale != 0.) {
				for (k=i; k<=m; k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}

				f = a[i][i];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][i] = f-g;

				for (j=l; j<=n; j++) {
					for (s=0.,k=i; k<=m; k++)
						s += a[k][i]*a[k][j];
					f = s/h;
					for (k=i; k<=m; k++)
						a[k][j] += f*a[k][i];
				}

				for (k=i; k<=m; k++)
					a[k][i] *= scale;
			}
		}

		w[i] = scale*g;
		g = s = scale = 0.;

		if (i<=m && i!=n) {
			for (k=l; k<=n; k++)
				scale += fabs(a[i][k]);

			if (scale!=0.) {
				double ih, isc=1./scale;

				for (k=l; k<=n; k++) {
					a[i][k] *= isc;
					s += a[i][k]*a[i][k];
				}

				f = a[i][l];
				g = -SIGN(sqrt(s),f);
				h = f*g - s; ih=1.0/h;
				a[i][l] = f-g;

				for (k=l; k<=n; k++)
					rv1[k] = a[i][k]*ih;

				for (j=l; j<=m; j++) {
					for (s=0.,k=l; k<=n; k++)
						s += a[j][k]*a[i][k];
					for (k=l; k<=n; k++)
						a[j][k] += s*rv1[k];
				}

				for (k=l; k<=n; k++)
					a[i][k] *= scale;
			}
		}

		anorm = FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}

	for (i=n; i>=1; i--) {
		if (i<n) {
			if (g!=0) {
				double ig=1./g;

				for (j=l; j<=n; j++)
					v[j][i] = (a[i][j]/a[i][l])*ig;

				for (j=l; j<=n; j++) {
					for (s=0.,k=l; k<=n; k++)
						s += a[i][k]*v[k][j];
					for (k=l; k<=n; k++)
						v[k][j] += s*v[k][i];
				}
			}

			for (j=l; j<=n; j++)
				v[i][j] = v[j][i] = 0.;
		}

		v[i][i] = 1.;
		g = rv1[i];
		l = i;
	}

	for (i=IMIN(m,n); i>=1; i--) {
		l = i+1;
		g = w[i];

		for (j=l; j<=n; j++)
			a[i][j] = 0.;

		if (g!=0) {
			g = 1./g;

			for (j=l; j<=n; j++) {
				for (s=0.,k=l; k<=m; k++)
					s += a[k][i]*a[k][j];

				f = (s/a[i][i])*g;

				for (k=i; k<=m; k++)
					a[k][j] += f*a[k][i];
			}
			for (j=i; j<=m; j++)
				a[j][i] *= g;
		}
		else	for (j=i; j<=m; j++)
				a[j][i] = 0.;

		a[i][i] += 1.;
	}

	for (k=n; k>=1; k--) {
		for (its=1; its<=30; its++){
			flag = 1;

			for (l=k; l>=1; l--) {
				nm = l-1;

				if (fabs(rv1[l])+anorm == anorm) {
					flag = 0;
					break;
				}

				if (fabs(w[nm])+anorm == anorm)
					break;
			}

			if (flag) {
				c = 0.;
				s = 1.;

				for (i=l; i<=k; i++) {
					f      = s*rv1[i];
					rv1[i] = c*rv1[i];

					if (fabs(f)+anorm == anorm)
						break;

					g = w[i];
					h = sqrt (f*f+g*g);
					w[i] = h;
					h = 1./h;
					c = g*h;
					s = -f*h;

					for (j=1; j<=m; j++) {
					    y = a[j][nm];
					    z = a[j][i];
					    a[j][nm] = y*c+z*s;
					    a[j][i]  = z*c-y*s;
					}
				}
			}

			z = w[k];

			if (l == k) {
				if (z < 0.) {
					w[k] = -z;
					for (j=1; j<=n; j++)
						v[j][k] = -v[j][k];
				}
				break;
			}

			if (its == 30) {
				fprintf (stderr,
				"no convergence in 30 svdcmp iterations\n");
				return -1;
			}

			x  = w[l];
			nm = k-1;
			y  = w[nm];
			g  = rv1[nm];
			h  = rv1[k];
			f  = ((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g  = sqrt(f*f+1.);
			f  = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.;

			for (j=l; j<=nm; j++) {
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = sqrt(f*f+h*h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;

				for (jj=1; jj<=n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}
				z = sqrt (f*f+h*h);
				w[j] = z;

				if (z!=0.) {
					z = 1./z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;

				for (jj=1; jj<=m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.;
			rv1[k] = f;
			w[k]   = x;
		}
	}

	return 0;
}

int
solvpd (double **drdp,  /* [1...m][1...n] */
	int m,          /* m=nused        */
	int n,          /* n=nlay         */
	double *dr,     /* [1...m]        */
	double *dp,     /* [1...n]        */
	double trunc,
	int   ntrunc)	
{
	double *vv[m+1], ww[n+1], wmin, wmax=0.0;
	int   i, k, err;
	/* solve the linear system using the
	 * singular value decomposition (SVD)
	 */

	for (i=0; i<m; i++) vv[i+1] = new double[n+1];
	
	err = svdcmp(drdp, m, n, ww, vv);

	if (err==0) {
		/* eigenvalue truncation */
		wmax = 0.;
		for (k=1; k<=n; k++)
			if (ww[k]>wmax) wmax = ww[k];

		wmin = trunc*wmax;
		for (k=1; k<=n; k++)
			if (ww[k]<wmin)
				ww[k] = 0.;

		svbksb(drdp, ww, vv, m, n, dr, dp);
	}

	for (i=0; i<m; i++) delete [] vv[i+1];
	
	return err;
}
