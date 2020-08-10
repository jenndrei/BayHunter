# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import numpy as np
from BayHunter.surfdisp96_ext import surfdisp96


class SurfDisp(object):
    """Forward modeling of dispersion curves based on surf96 (Rob Herrmann).

    The quick fortran routine is from Hongjian Fang:
        https://github.com/caiweicaiwei/SurfTomo

    BayHunter.SurfDisp leaning on the python wrapper of Marius Isken:
        https://github.com/miili/pysurf96
    """

    def __init__(self, obsx, ref):
        self.obsx = obsx
        self.kmax = obsx.size
        self.ref = ref

        self.modelparams = {
            'mode': 1,  # mode, 1 fundamental, 2 first higher
            'flsph': 0  # flat earth model
            }

        self.wavetype, self.veltype = self.get_surftags(ref)

        if self.kmax > 60:
            message = "Your observed data vector exceeds the maximum of 60 \
periods that is allowed in SurfDisp. For forward modeling SurfDisp will \
reduce the samples to 60 by linear interpolation within the given period \
span.\nFrom this data, the dispersion velocities to your observed periods \
will be determined. The precision of the data will depend on the distribution \
of your samples and the complexity of the input velocity-depth model."
            self.obsx_int = np.linspace(obsx.min(), obsx.max(), 60)
            print(message)

    def set_modelparams(self, **mparams):
        self.modelparams.update(mparams)

    def get_surftags(self, ref):
        if ref == 'rdispgr':
            return (2, 1)

        elif ref == 'ldispgr':
            return (1, 1)

        elif ref == 'rdispph':
            return (2, 0)

        elif ref == 'ldispph':
            return (1, 0)
        else:
            tagerror = "Reference is not available in SurfDisp. If you defined \
a user Target, assign the correct reference (target.ref) or update the \
forward modeling plugin with target.update_plugin(MyForwardClass()).\n \
* Your ref was: %s\nAvailable refs are: rdispgr, ldispgr, rdispph, ldispph\n \
(r=rayleigh, l=love, gr=group, ph=phase)" % ref
            raise ReferenceError(tagerror)

    def get_modelvectors(self, h, vp, vs, rho):
        nlayer = len(h)
        thkm = np.zeros(100)
        thkm[:nlayer] = h

        vpm = np.zeros(100)
        vpm[:nlayer] = vp

        vsm = np.zeros(100)
        vsm[:nlayer] = vs

        rhom = np.zeros(100)
        rhom[:nlayer] = rho

        return thkm, vpm, vsm, rhom

    def run_model(self, h, vp, vs, rho, **params):
        """ The forward model will be run with the parameters below.

        thkm, vpm, vsm, rhom: model for dispersion calculation
        nlayer - I4: number of layers in the model
        iflsph - I4: 0 flat earth model, 1 spherical earth model
        iwave - I4: 1 Love wave, 2 Rayleigh wave
        mode - I4: ith mode of surface wave, 1 fundamental, 2 first higher, ...
        igr - I4: 0 phase velocity, > 0 group velocity
        kmax - I4: number of periods (t) for dispersion calculation
        t - period vector (t(NP))
        cg - output phase or group velocities (vector,cg(NP))

        """
        nlayer = len(h)
        h, vp, vs, rho = self.get_modelvectors(h, vp, vs, rho)

        iflsph = self.modelparams['flsph']
        mode = self.modelparams['mode']
        iwave = self.wavetype
        igr = self.veltype

        if self.kmax > 60:
            kmax = 60
            pers = self.obsx_int

        else:
            pers = np.zeros(60)
            kmax = self.kmax
            pers[:kmax] = self.obsx

        dispvel = np.zeros(60)  # result
        error = surfdisp96(h, vp, vs, rho, nlayer, iflsph, iwave,
                           mode, igr, kmax, pers, dispvel)

        if error == 0:
            if self.kmax > 60:
                disp_int = np.interp(self.obsx, pers, dispvel)
                return self.obsx, disp_int

            return pers[:kmax], dispvel[:kmax]

        return np.nan, np.nan
