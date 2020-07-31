/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

# ifndef _MODEL_H_
# define _MODEL_H_

# include "defaults.h"

# define LOWER_HALFSPACE -1

class FlatLayer
{
public:
  FlatLayer ();
  FlatLayer (const FlatLayer &);
 ~FlatLayer () {}

  void      setDepth (double);
  void      setThick (double);
  void      setPVel  (double);
  void      setSVel  (double  vs = -1.);
  void      setDens  (double rho = -1.);
  void      setAll   (double z,  double h,
		      double vp, double vs, double rho,
		      double qp=DEF_QS, double qs=DEF_QS);

  // members for data retrieval
  double     getDepth ();
  double     getThick ();
  double     getPVel  ();
  double     getSVel  ();
  double     getDens  ();
  void      getAll   (double *z,  double *h,
		      double *vp, double *vs, double *rho,
		      double *qp=0, double *qs=0);

  int       realVerticalSlownesses (double, double*, double*);
  //int       complexVerticalSlownesses (double, Complex*, Complex*);

  // earth-flattening transform
  enum      Type {Spherical=0, Flat=1};
  void      flatten ();
  void      unflatten ();

  void      perturb (double dv=0.001);

  int       isUpperHalfspace ();
  int       isLowerHalfspace ();

//private:
  Type      type;

  double      z;     // layer thickness
  double      h;     // layer thickness
  double     vp;     // P wave velocity
  double     vs;     // S wave velocity
  double     vpvs;   // vp/vs
  double     rh;     // density
  double     qp;     // P wave quality factor
  double     qs;     // S wave quality factor

  int       flattened; //
};

// The class `FlatModel' describes an earthmodel consisting of only
// flat, isotropic, homogeneous layers.
//
// The layer indexing is as follows:
// The topmost layer has the index 1
// The lower halfspace has the index n, so that if the lower
// halfspace is considered the lowermost `layer', we get a number
// of layers of n. However, the memory is physically allocated
// for n+1 layers because the layer index 0 is also allocated
// but not used. In fact, one could consider the medium with
// index 0 be the upper halfspace.

class FlatModel
{
public:
    FlatModel ();
    FlatModel (FlatModel const &);
   ~FlatModel ();

    int        isNull() { return nLay>0 ? 0 : -1; }
    int        numLayers () {return nLay;}
    FlatLayer& layer (int index);
    int        setLayers (int, const FlatLayer*);

    void       flatten ();
    void       unflatten ();

//private:
    int        nLay;      // number of layers in model
    FlatLayer *lay;       // pointer to array of layers
    int        flattened; // flag for earth-flattening transform
};

# endif
