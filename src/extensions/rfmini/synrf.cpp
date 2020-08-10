/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "synrf.h"
# include "defaults.h"
# include "model.h"
# include "wave.h"

void
synrf (Wave &wave,
	int nlay,
	double *z,  double *vp, double *vs,
	double *rh, double *qp, double *qs,
	RfnParam &par,
	double *syn_z, double *syn_r, double *syn_rf)
{
	int i, options = 0; // currently practically a dummy
	double fref = 1.; // currently this cannot be changed
	FlatModel model;
        FlatLayer lay[nlay+1];
        for (i=0; i<nlay-1; i++)
                lay[i+1].setAll(z[i], z[i+1]-z[i], vp[i], vs[i],
			       rh[i], qp[i], qs[i]);
        lay[nlay].setAll(z[nlay-1], -1, vp[nlay-1], vs[nlay-1],
		        rh[nlay-1],     qp[nlay-1], qs[nlay-1]);
        model.setLayers(nlay, lay);
	model.flatten();
InvParam ipar;

	// compute the responses and the receiver function
	calcresp (model,
              wave,
              fref,       // reference freq. for q.f.
              par.nsamp,  // number of samples
              par.fsamp,  // sampling frequency
              par.water,  // waterlevel
              par.tshift, // time shift
              par.a,      // Gauss parameter
              par.vptop, par.vstop,
              syn_z,      // << vertical responses
              syn_r,      // << radial responses
              syn_rf,     // << receiver function
		0,	  // XXX first sample used in inversion
              0,          // number of used samples 
              NULL,       // << pd matrix
              0, 0, ipar, 0.001,
              options);
}

