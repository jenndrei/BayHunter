# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import numpy as np
from BayHunter import rfmini


class RFminiModRF(object):
    """Forward modeling of receiver functions based on SeisPy (Joachim Saul).
    """
    def __init__(self, obsx, ref):
        self.ref = ref
        self.obsx = obsx
        self._init_obsparams()

        if self.ref in ['prf', 'seis']:
            self.modelparams = {'wtype': 'P'}
        elif self.ref in ['srf']:
            self.modelparams = {'wtype': 'SV'}

        self.modelparams.update(
            {'gauss': 1.0,
             'p': 6.4,
             'water': 0.001,
             'nsv': None
             })

        self.keys = {'z': '%.2f',
                     'vp': '%.4f',
                     'vs': '%.4f',
                     'rho': '%.4f',
                     'qp': '%.1f',
                     'qs': '%.1f',
                     'n': '%d'}

    def _init_obsparams(self):
        """Extract parameters from observed x-data (time vector).

        fsamp = sampling frequency in Hz
        tshft = time shift by which the RF is shifted to the left
        nsamp = number of samples, must be 2**x
        """

        # get fsamp
        deltas = np.round((self.obsx[1:] - self.obsx[:-1]), 4)
        if np.unique(deltas).size == 1:
            dt = float(deltas[0])
            self.fsamp = 1. / dt
        else:
            raise ValueError("Target: %s. Sampling rate must be constant."
                             % self.ref)
        # get tshft
        self.tshft = -self.obsx[0]

        # get nsamp
        ndata = self.obsx.size
        self.nsamp = 2.**int(np.ceil(np.log2(ndata * 2)))

    def write_startmodel(self, h, vp, vs, rho, modfile, **params):
        qp = params.get('qp', np.ones(h.size) * 500)
        qs = params.get('qs', np.ones(h.size) * 225)

        z = np.cumsum(h)
        z = np.concatenate(([0], z[:-1]))

        mparams = {'z': z, 'vp': vp, 'vs': vs, 'rho': rho,
                   'qp': qp, 'qs': qs}
        mparams = dict((a, b) for (a, b) in mparams.items()
                       if b is not None)
        pars = mparams.keys()

        nkey = 0
        header = []
        mline = []
        data = np.empty((len(pars), mparams[pars[0]].size))
        for key in ['z', 'vp', 'vs', 'rho', 'qp', 'qs']:
            if key in pars:
                header.append(key)
                mline.append(self.keys[key])
                data[nkey, :] = mparams[key]
                nkey += 1

        header = '\t'.join(header) + '\n'
        mline = '\t'.join(mline) + '\n'

        with open(modfile, 'w') as f:
            f.write(header)
            for i in np.arange(len(data[0])):
                f.write(mline % tuple(data.T[i]))

    def set_modelparams(self, **mparams):
        self.modelparams.update(mparams)

    def compute_rf(self, h, vp, vs, rho, **params):
        """
        Compute RF using self.modelsparams (dict) for parameters.
        e.g. usage: self.set_modelparams(gauss=1.0)

        Parameters are:
        # z  depths of the top of each layer
        gauss: Gauss parameter
        water: water level
        p: angular slowness in sec/deg
        wtype: type of incident wave; must be 'P' or 'SV'
        nsv: tuple with near-surface S velocity and Poisson's ratio
            (will be computed by input model, if None)
        """
        gauss = self.modelparams['gauss']
        water = self.modelparams['water']
        p = self.modelparams['p']
        wtype = self.modelparams['wtype']
        nsv = self.modelparams['nsv']

        qp = params.get('qp', np.ones(h.size) * 500.)
        qs = params.get('qs', np.ones(h.size) * 225.)

        z = np.cumsum(h)
        z = np.concatenate(([0], z[:-1]))

        nsvp, nsvs = float(vp[0]), float(vs[0])
        vpvs = nsvp / nsvs
        poisson = (2 - vpvs**2)/(2 - 2 * vpvs**2)

        if nsv is None:
            nsv = nsvs

        time = np.arange(self.nsamp) / self.fsamp - self.tshft

        fp, fsv, qrf = rfmini.synrf(
            z, vp, vs, rho, qp, qs,
            p, gauss, self.nsamp, self.fsamp,
            self.tshft, nsv, poisson, wtype)

        # must be converted to float64
        qrfdata = qrf.astype(float)

        return time[:self.obsx.size], qrfdata[:self.obsx.size]

    def run_model(self, h, vp, vs, rho, **params):

        assert h.size == vp.size == vs.size == rho.size

        h = h.astype(float)
        vp = vp.astype(float)
        vs = vs.astype(float)
        rho = rho.astype(float)

        time, qrf = self.compute_rf(h, vp, vs, rho, **params)
        return time, qrf
