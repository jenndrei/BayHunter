/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

#include "Complex.h"
#include <math.h>

void
ccfork (int n, Complex *x, int signi)
{
	//  FFT-Routine fuer C++
	//  fuer N = 2er-Potenz
	//  entwickelt aus der Fortran-Routine FORK
	//  Indizes fuer x von 0...N-1
	//  Bedeutung von signi:
	//  signi =-1 => Fourier-Transformation 
	//  signi = 1 => inverse Fourier-Transformation 
	//
	//  22.05.96 J.Saul    saul@dkrz.de

	Complex w,tmp;
	double sc;
	int i,istep,j=0,l,m;
  
	// normalization factor
	sc = sqrt(1./(double)n);

	for (i=0; i<n; i++) {
		if (i<=j) {
			tmp  = x[j]*sc;
			x[j] = x[i]*sc;
			x[i] = tmp;
		}
		m = n>>1;

		do {
			if (j<m) break;
			j -= m;
			m >>= 1;
		} while (m>=1);
		j += m;
	}
	l = 1;

	do {
		istep = 2*l;
		for (m=0; m<l; m++) {
			w = exp(Complex(0.0,
				M_PI*(double)(signi*m)/(double)l));
			for (i=m; i<n; i+=istep) {
				tmp    = w*x[i+l];
				x[i+l] = x[i]-tmp;
				x[i]  +=      tmp;
			}
		}
		l = istep;
	} while (l<n);
}
