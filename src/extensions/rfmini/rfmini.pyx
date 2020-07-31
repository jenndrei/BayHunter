from numpy cimport ndarray as array
import numpy

cdef extern from "Python.h":
    object PyComplex_FromDoubles(double, double)

cdef extern from "numpy/arrayobject.h":
    cdef enum PyArray_TYPES:
        PyArray_CHAR, PyArray_UBYTE, PyArray_SBYTE, PyArray_SHORT,
        PyArray_USHORT, PyArray_INT, PyArray_UINT, PyArray_LONG, PyArray_FLOAT,
        PyArray_DOUBLE, PyArray_CFLOAT, PyArray_CDOUBLE, PyArray_OBJECT,
        PyArray_NTYPES, PyArray_NOTYPE

# force inclusion of some header files
#cdef extern from "seispy_unique.h": pass
cdef extern from "numpy/arrayobject.h": pass

#cdef extern from "seispy/seispy.h":
#
#    ctypedef struct Model_1D:
#        int     nlay
#        double  *z
#        double  *vp
#        double  *vs
#        double  *rh
#
#    ctypedef struct ReceiverFunction:
#        int     nsamp
#        double  *data
#        double  fsamp, tshift, slow, azi
#
#    Model_1D *Model_1D_new()
#    int     Model_1D_init(Model_1D *mod, int nlay, double *z, double *vp, double *vs, double *rh)
#    void    Model_1D_delete(Model_1D *mod)
#
#   ReceiverFunction* ReceiverFunction_new()
#   int     ReceiverFunction_init(ReceiverFunction *rf, int nsamp, double *data, double fsamp, double tshift, double slow, double azi)
#   int     ReceiverFunction_moco(ReceiverFunction *rf, Model_1D *mod, int phase, double p0)
#   int     ReceiverFunction_zmigr(ReceiverFunction *rf, Model_1D *mod, int phase, ReceiverFunction *zf)
#   void    ReceiverFunction_delete(ReceiverFunction *rf)

#   int     rfzmigr(double *rf_data, int rf_nsamp, double rf_fsamp, double rf_slow, double rf_tshift, int nz1, double *z, double *vp, double *vs, double *st_data, int st_nsamp, double st_fsamp, int phase)
    

cdef extern int synrf_cwrap(
        int nsamp, double fsamp, double tshift, double p, double a,
        double nsv, double sigma, int waveno, int nlay,
        double  *z, double *vp, double *vs,
        double *rh, double *qp, double *qs,
        double *fz, double *fr, double *rf)

cdef extern int coeff_cwrap(double p,
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

cdef extern int coeffs_cwrap(double p,
            double vp, double vs, double rh,
            double *ru11r, double *ru12r, double *ru21r, double *ru22r,
            double *ru11i, double *ru12i, double *ru21i, double *ru22i,
            double *rhur,  double *rhui)


def synrf(array  z_arr, array vp_arr, array vs_arr,
          array rh_arr, array qp_arr, array qs_arr,
          double p, double a, int nsamp, double fsamp, double tshift,
          double nsv, double sigma, str wave):

    cdef array fz_arr, fr_arr, rf_arr
    cdef double *z
    cdef double *vp
    cdef double *vs
    cdef double *rh
    cdef double *qp
    cdef double *qs
    cdef double *fz
    cdef double *fr
    cdef double *rf
    cdef int nlay, waveno

    try:
        waveno = ["P", "SV", "SH"].index(wave)
    except ValueError:
        raise ValueError, "wave must be 'P', 'SV' or 'SH', not '%s'" % wave
    fz_arr = numpy.zeros(nsamp, dtype=numpy.float64)
    fr_arr = numpy.zeros(nsamp, dtype=numpy.float64)
    rf_arr = numpy.zeros(nsamp, dtype=numpy.float64)

    nlay = len(z_arr) 
    z  = <double*>  z_arr.data
    vp = <double*> vp_arr.data
    vs = <double*> vs_arr.data
    rh = <double*> rh_arr.data
    qp = <double*> qp_arr.data
    qs = <double*> qs_arr.data

    fz = <double*> fz_arr.data
    fr = <double*> fr_arr.data
    rf = <double*> rf_arr.data

    synrf_cwrap(nsamp, fsamp, tshift, p, a, nsv, sigma, waveno, nlay,
                z, vp, vs, rh, qp, qs, fz, fr, rf)

    return fz_arr,fr_arr,rf_arr


cdef d_copy_cast(array arr):
    # Casts AND copies an array to double. Array is copied
    # even if the original array already is of type double.
    return arr.astype(numpy.float)  # XXX Maybe better to use PyArray_Cast() here?

#def mocorr(array rf_arr, array  z_arr,
#           array vp_arr, array vs_arr,
#           double rf_slow, double p0, double rf_fsamp, int phase):
#    # leaves the input arrays untouched
#
#    cdef double *rf_data
#    cdef int i, err, rf_nsamp
#    cdef Model_1D *mod
#    cdef ReceiverFunction *rf
#
#    z_arr  = d_copy_cast(z_arr)
#    vp_arr = d_copy_cast(vp_arr)
#    vs_arr = d_copy_cast(vs_arr)
#    rf_arr = d_copy_cast(rf_arr)
#
#    rf_nsamp = rf_arr.dimensions[0]
#    rf_data  = <double*> rf_arr.data
#    rf_slow  = rf_slow/111.195
#    p0       = p0/111.195
#
#    mod = Model_1D_new()
#    rf  = ReceiverFunction_new()
#    if not mod or not rf:
#        raise MemoryError
#    err = Model_1D_init(mod, z_arr.dimensions[0],
#                    <double*>  z_arr.data,
#                    <double*> vp_arr.data,
#                    <double*> vs_arr.data, <double*>0)
#    if err: raise MemoryError
#    err = ReceiverFunction_init(rf, rf_nsamp, rf_data, rf_fsamp, 0, rf_slow, 0)
#    if err: raise MemoryError
#
#    ReceiverFunction_moco(rf, mod, phase, p0)
#    for i from 0 <= i < rf_nsamp:
#        rf_data[i] = rf.data[i]
#    ReceiverFunction_delete(rf)
#    Model_1D_delete(mod)
#    return rf_arr
#
#def tmigr(array rf_arr, double rf_slow, double rf_fsamp, double rf_azi,
#          double rf_sta_x, double rf_sta_y,
#          array z_arr, array vp_arr, array vs_arr,
#          double p0, int phase,
#          array st_arr, array wt_arr):
#
#    cdef double *rf_data
#    cdef int i, err, rf_nsamp
#    cdef Model_1D *mod
#    cdef ReceiverFunction *rf
#    cdef ReceiverFunction *zf
#
#    z_arr    = d_copy_cast(z_arr)
#    vp_arr   = d_copy_cast(vp_arr)
#    vs_arr   = d_copy_cast(vs_arr)
#
#    rf_arr   = d_copy_cast(rf_arr)
#    rf_nsamp = rf_arr.dimensions[0]
#    rf_data  = <double*> rf_arr.data
#    rf_slow  = rf_slow/111.195
#
#    mod = Model_1D_new()
#    rf  = ReceiverFunction_new()
#    zf  = ReceiverFunction_new()
#    if not mod or not rf:
#        raise MemoryError
#    err = Model_1D_init(mod, z_arr.dimensions[0],
#                    <double*>  z_arr.data,
#                    <double*> vp_arr.data,
#                    <double*> vs_arr.data, <double*>0)
#    if err: raise MemoryError
#    err = ReceiverFunction_init(rf, rf_nsamp, rf_data, rf_fsamp, 0, rf_slow, 0)
#    if err: raise MemoryError
#
#### ReceiverFunction_tmigr(rf, mod, phase, st_grd, wt_grd)
#### XXX
#                      
#def zmigr(array rf_arr, array  z_arr,
#          array vp_arr, array vs_arr,
#          double rf_slow, double rf_fsamp, double rf_tshift,
#          double zmax, double dz, int phase):
#
#    cdef double *rf_data
#    cdef double *zf_data
#    cdef double  zf_fsamp
#    cdef int i, rf_nsamp, zf_nsamp, nz1
#    cdef array zf_arr
#    cdef Model_1D *mod
#    cdef ReceiverFunction *zf
#    cdef ReceiverFunction *rf
#
#    z_arr    = d_copy_cast(z_arr)
#    vp_arr   = d_copy_cast(vp_arr)
#    vs_arr   = d_copy_cast(vs_arr)
#
#    rf_arr   = d_copy_cast(rf_arr)
#    rf_nsamp = rf_arr.dimensions[0]
#    rf_data  = <double*> rf_arr.data
#    rf_slow  = rf_slow/111.195
#
#    zf_nsamp = <int>(zmax/dz+0.9)
#    zf_arr   = vector_new_d(zf_nsamp)
#    zf_data  = <double*> zf_arr.data
#    zf_fsamp = 1./dz
#
#    mod = Model_1D_new()
#    rf  = ReceiverFunction_new()
#    zf  = ReceiverFunction_new()
#    if not mod or not rf:
#        raise MemoryError
#    err = Model_1D_init(mod, z_arr.dimensions[0],
#                    <double*>  z_arr.data,
#                    <double*> vp_arr.data,
#                    <double*> vs_arr.data, <double*>0)
#    if err: raise MemoryError
#    err = ReceiverFunction_init(rf, rf_nsamp, rf_data, rf_fsamp, 0, rf_slow, 0)
#    if err: raise MemoryError
#    err = ReceiverFunction_init(zf, zf_nsamp, zf_data, zf_fsamp, 0, 0, 0)
#    if err: raise MemoryError
#
#    ReceiverFunction_zmigr(rf, mod, phase, zf)
#
#    for i from 0 <= i < zf_nsamp:
#        zf_data[i] = zf.data[i]
#    ReceiverFunction_delete(rf)
#    ReceiverFunction_delete(zf)
#    Model_1D_delete(mod)
#
#    return zf_arr


def coeff(double p, double vp1, double vs1, double rh1,
                    double vp2, double vs2, double rh2, int dis):

    """
coeff(p, vp1, vs1, rh1, vp2, vs2, rh2, dis)

Computes the reflection/transmission coefficients for a plane wave with
slowness p, incident at a plane, horizontal interface between two
halfspaces with properties vp1,vs1,rh1 and vp2,vs2,rh2. If 'dis' is not
0, displacement coefficients are computed; otherwise they are potential.

The coefficients are returned as a tuples of five matrices, each
consisting of 4 complex values.
(rd11,rd12,rd21,rd22),(td11,td12,td21,td22),
(ru11,ru12,ru21,ru22),(tu11,tu12,tu21,tu22),
(rdh,tdh,ruh,tuh)
The first 4 matrices contain the P/SV coefficients, the 5th matrix
contains the corresponding SH coefficients.
"""

    cdef double rd11r, rd12r, rd21r, rd22r, \
                rd11i, rd12i, rd21i, rd22i, \
                td11r, td12r, td21r, td22r, \
                td11i, td12i, td21i, td22i, \
                ru11r, ru12r, ru21r, ru22r, \
                ru11i, ru12i, ru21i, ru22i, \
                tu11r, tu12r, tu21r, tu22r, \
                tu11i, tu12i, tu21i, tu22i, \
                rhdr,  thdr,  rhur,  thur,  \
                rhdi,  thdi,  rhui,  thui

    coeff_cwrap(p, vp1, vs1, rh1, vp2, vs2, rh2, dis,
                &rd11r, &rd12r, &rd21r, &rd22r,
                &rd11i, &rd12i, &rd21i, &rd22i,
                &td11r, &td12r, &td21r, &td22r,
                &td11i, &td12i, &td21i, &td22i,
                &ru11r, &ru12r, &ru21r, &ru22r,
                &ru11i, &ru12i, &ru21i, &ru22i,
                &tu11r, &tu12r, &tu21r, &tu22r,
                &tu11i, &tu12i, &tu21i, &tu22i,
                &rhdr,  &thdr,  &rhur,  &thur,
                &rhdi,  &thdi,  &rhui,  &thui)

    return (  ( PyComplex_FromDoubles(rd11r, rd11i),
                PyComplex_FromDoubles(rd12r, rd12i),
                PyComplex_FromDoubles(rd21r, rd21i),
                PyComplex_FromDoubles(rd22r, rd22i) ),
              ( PyComplex_FromDoubles(td11r, td11i),
                PyComplex_FromDoubles(td12r, td12i),
                PyComplex_FromDoubles(td21r, td21i),
                PyComplex_FromDoubles(td22r, td22i) ),
              ( PyComplex_FromDoubles(ru11r, ru11i),
                PyComplex_FromDoubles(ru12r, ru12i),
                PyComplex_FromDoubles(ru21r, ru21i),
                PyComplex_FromDoubles(ru22r, ru22i) ),
              ( PyComplex_FromDoubles(tu11r, tu11i),
                PyComplex_FromDoubles(tu12r, tu12i),
                PyComplex_FromDoubles(tu21r, tu21i),
                PyComplex_FromDoubles(tu22r, tu22i) ),
              ( PyComplex_FromDoubles(rhdr,  rhdi),
                PyComplex_FromDoubles(thdr,  thdi),
                PyComplex_FromDoubles(rhur,  rhui),
                PyComplex_FromDoubles(thur,  thui) )  )

def coeffs(double p, double vp, double vs, double rh):

    cdef double ru11r, ru12r, ru21r, ru22r, \
                ru11i, ru12i, ru21i, ru22i, \
                rhur, rhui

    coeffs_cwrap(p, vp, vs, rh,
                &ru11r, &ru12r, &ru21r, &ru22r,
                &ru11i, &ru12i, &ru21i, &ru22i,
                &rhur, &rhui)

    return (  ( PyComplex_FromDoubles(ru11r, ru11i),
                PyComplex_FromDoubles(ru12r, ru12i),
                PyComplex_FromDoubles(ru21r, ru21i),
                PyComplex_FromDoubles(ru22r, ru22i) ),
                PyComplex_FromDoubles(rhur,  rhui) )
