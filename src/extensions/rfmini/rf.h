/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

# ifndef _INVRF_H_
# define _INVRF_H_

# include <string.h>
# include "defaults.h"


// parameters for receiver function computation
struct RfnParam
{
	double	p,       // horizontal slowness
		fsamp,   // sampling frequency
		water,   // waterlevel
		a,       // gauss parameter
		tshift,  // time shift  
		vptop, vstop; // rotation velocities
	int	nsamp;   // number of samples

	RfnParam () {
		// initialize with default values
		p      = PDEF_P;
		fsamp  = PDEF_FSAMP;
		water  = PDEF_WATER;
		a      = PDEF_GAUSS;
		tshift = PDEF_TSHIFT;
		vstop  = vptop = 0.0;
		nsamp  = PDEF_NSAMP;
	}
};

// parameters for receiver function inversion
struct InvParam
{
	double	smooth,  // model smoothness
		dsdz,    // increase of smoothness with depth
		speed,   // inversion speed
		length,  // length of time window for inversion
		svtf;    // fraction of max. singular value after which to truncate
    
	int	imin,    // minimum number of iterations
		imax,    // maximum number of iterations
		laymin,  // topmost variable layer
		laymax,  // lowermost variable layer
		direct,  // direct inversion
		verbose, // verbose mode
		save_pd, // save partial derivative matrix
		save_intermediate, // save intermediate models
		svtn;    // number of singular values after which to truncate

	InvParam () {
		// initialize with default values
		smooth  = PDEF_SMOOTH;
		dsdz    = 0.;
		speed   = PDEF_ISPEED;
		length  = 20.0;
		imin    = 6;
		imax    = 30;
		laymin  = 0;
		laymax  = -1;
		direct  = 0;
		verbose = 0; 
		save_pd = 0;
		save_intermediate = 0;
		svtn    = 0;
		svtf    = 0.0;
	}  
};

# endif
