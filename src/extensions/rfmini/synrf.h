/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

#ifndef _SYNRF_H_
#define _SYNRF_H_

#include <stdio.h>
#include "cmat2.h"
#include "model.h"
#include "wave.h"

#include "rf.h" // XXX

double rhovp (double vp);
double vsvp  (double vp);

void  cmprf (Complex *cr, Complex *cz, int nsamp, double fsamp,
             double water, double tshift, double a, double angle, Complex *crf);
Cmat2 cmatrh(double p, double vp, double vs);

void ztod (int nlay, double *z, double *d);

void
calcresp (FlatModel model,                                        // model
          Wave      wave,
          double     fref,        // reference freq. for q.f.      // model
          int       nsamp,       // number of samples             // record
          double     fsamp,       // sampling frequency            // record
          double     water,       // waterlevel                    // RF
          double     tshift,      // time shift                    // RF
          double     a,           // Gauss parameter               // RF
          double     vptop, double vstop,// NEW                     // RF
          double     *zz,         // << vertical responses         // record
          double     *rr,         // << radial responses           // record
          double     *rf,         // << receiver function          // RF
          int first,
	  int       nused,       // number of used samples 
          double     **drdp,      // << pd matrix
          // drdp==NULL indicates that no pd matrix
          // needs to be computed (previously: pdflag)
          int       laymin,      // first perturbed layer         // model
          int       laymax,      // last  perturbed layer         // model
InvParam &ipar,
	  double     pert=0.001,  // perturbation
          int       options=0);  // option flags



# define SUPPRESS_MULTIPLES   1
# define WITHOUT_ANELASTICITY 2

# endif
