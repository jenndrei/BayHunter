/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"

FlatLayer::FlatLayer ()
{
	z  = h  =      0.;
	vp = vs = rh = 0.;
	qp = qs =      0.;
	vpvs = 1.73; // dummy at this point
    
	flattened = 0;
}

FlatLayer::FlatLayer (const FlatLayer &other)
{
	z	= other.z;
	h	= other.h;
	vp	= other.vp;
	vs	= other.vs;
	vpvs	= other.vpvs;
	rh	= other.rh;
	qp	= other.qp;
	qs	= other.qs;

	flattened = other.flattened;
}


void
FlatLayer::setThick (double xh)
{
	h = xh;
}

void
FlatLayer::setPVel (double xvp)
{
	vp = xvp;
}

void
FlatLayer::setSVel (double xvs)
{
	// Sets the S velocity of the layer.
	//
	//  // If the argument is omitted, the S velocity is
	//  // calculated from the P velocity.

	if (xvs >= 0.)
		vs = xvs;
	else {
		//    double vpvs = vp/vs;
		//    vs = xvs;
	}
}

double rho_vp (double vp);

void
FlatLayer::setDens (double xrh)
{
	// Sets the density of the layer.
	//
	// If the argument is omitted, the density is
	// calculated from the P velocity.

	if (xrh >= 0.)
		rh = xrh;
	else	rh = rho_vp (vp);

	// The vp-rho relationship of Bertheussen (1977) is
	// only valid for high-velocity crystalline rocks. The
	// low-velocity modification approximates the vp-rho
	// relationship for sediments of Gardner et al. (1974)
	// for low velocities (vp < 5km/s).
}

void
FlatLayer::setAll (double xz, double xh, 
		   double xvp, double xvs, double xrh,
		   double xqp, double xqs)
{
	z	= xz;
	h	= xh;
	vp	= xvp;
	vs	= xvs;
	vpvs	= vp/vs;
	rh	= xrh;

	qp	= xqp;
	qs	= xqs;

//    fprintf(stderr,"%f %f %f %f %f %f %f %f\n",z,h,vp,vs,vpvs,rh,qp,qs);
}

double
FlatLayer::getDepth ()
{
	return z;
}

double
FlatLayer::getThick ()
{
	return h;
}

double
FlatLayer::getPVel ()
{
	return vp;
}

double
FlatLayer::getSVel ()
{
	return vs;
}

double
FlatLayer::getDens ()
{
	return rh;
}

void
FlatLayer::getAll (double *xz, double *xh, 
		   double *xvp, double *xvs, double *xrh,
		   double *xqp, double *xqs)
{
	if (xz)  *xz  = z;
	if (xh)  *xh  = h;
	if (xvp) *xvp = vp;
	if (xvs) *xvs = vs;
	if (xrh) *xrh = rh;

	if (xqp) *xqp = qp;
	if (xqs) *xqs = qs;
}


double
rho_vp (double vp)
{
    // Based on Bertheussen (1977): rho = 0.77 + 0.32*vp
    // modified for low velocity sediments.
    //
    // The vp-rho relationship of Bertheussen (1977) is
    // only valid for high-velocity crystalline rocks. The
    // low-velocity modification approximates the vp-rho
    // relationship for sediments of Gardner et al. (1974)
    // for low velocities (vp < 5km/s).

	return 0.77+0.32*vp                         // original Berteussen term
	  + 0.68*exp(-0.12*pow(vp-1.8,2))           // sediments
	  - 0.09*(vp-5.5)*exp(-0.7*pow(vp-5.5,2));  // transition
}

# define EarthRadius 6371.0

void
FlatLayer::perturb (double dv)
{
	// We perturb vs and adjust vp and rho accordingly.
	// Vp is adjusted by keeping the vp/vs ratio constant.
	// For adjusting rho, we first have to check whether
	// the layer is currently flattened, since we have to
	// use the unflattened vp value for settung rho.
	double q = vp/vs;

	vs += dv;    // perturb vs
	vp  = q*vs;  // adjust vp

	if (flattened) {
		double	r   = EarthRadius * exp(-z/EarthRadius),
			q   = r/EarthRadius,
			vp2 = q * vp;
			rh  = q * rho_vp (vp2);

		return;
	}

	rh = rho_vp (vp);
}

# undef EarthRadius

int
FlatLayer::isUpperHalfspace ()
{
	if (h > 0.)
		return 0;

	if (vp > 1. && rh > 0.1)
		return 0;

	return 1;
}

int
FlatLayer::isLowerHalfspace ()
{
	if (h > 0.)
		return 0;

	if (vp < 1.  &&  rh < 0.1)
		return 0;

	return 1;
}


# define EarthRadius 6371.0

void
FlatLayer::flatten ()
  // earth-flattening transformation
{
	if (flattened) return;

	double	zb = z + h, // depth of bottom
		r  = EarthRadius - z,
		q  = EarthRadius/r;
  
	z = EarthRadius * log(q);
	vp *= q;
	vs *= q;
	rh /= q;

	if( ! isLowerHalfspace()) {
		// compute transformed layer thickness
		r  = EarthRadius - zb;
		q  = EarthRadius/r;
		zb = EarthRadius * log(q);

		h  = zb - z;
	}

	// Not transforming the Q factors is certainly not
	// correct, but how are they transformed correctly???

	flattened = 1;
}

void
FlatLayer::unflatten ()
{
	if( ! flattened) return;

	double
		zb = z + h, // depth of bottom
		r  = EarthRadius*exp(-z/EarthRadius),
		q  = r/EarthRadius;
  
	z = EarthRadius - r;
	vp *= q;
	vs *= q;
	rh /= q;

	if ( ! isLowerHalfspace()) {
		// compute layer thickness
		r  = EarthRadius*exp(-zb/EarthRadius),
		zb = EarthRadius - r;

		h  = zb - z;
	}

	flattened = 0;
}

FlatModel::FlatModel()
{
	nLay      = 0;
	lay       = 0;
	flattened = 0;
}

FlatModel::FlatModel (FlatModel const &other)
{
	nLay = other.nLay;

	lay = new FlatLayer [nLay+1];

	for (int i=0; i<=nLay; i++)
		lay[i] = other.lay[i];

	flattened = other.flattened;
}

FlatModel::~FlatModel()
{
	if (lay) delete [] lay;
}

int
FlatModel::setLayers (int n, const FlatLayer *xlay)
{
	// re-allocate space for the layers
	if (lay) delete [] lay;
	lay = new FlatLayer [n+1];

	for (int i=1; i<=n; i++)
		lay[i] = xlay[i];

	nLay = n;

	return -1;
}

void
FlatModel::flatten ()
{
    if (flattened) return;

    for (int i=1; i<=nLay; i++)
        lay[i].flatten();

    flattened = -1;
}

void
FlatModel::unflatten ()
{
    if ( ! flattened) return;

    for (int i=1; i<=nLay; i++)
        lay[i].unflatten();

    flattened = 0;
}

FlatLayer&
FlatModel::layer (int index)
{
	if (index < 0  ||  index > nLay) {
		fprintf (stderr, "layer index %d out of range\n", index);
		exit(0);
	}
	return lay[index];
}


