# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import numpy as np
import quick_routine   # fictive routine


class MyForwardModel(object):
    """
    """
    def __init__(self, obsx, ref):
        self.ref = ref
        self.obsx = obsx

        # default parameters necessary for forward modeling
        # the dictionary can be updated by the user
        self.modelparams.update(
            {'test': 5,
             })

    def set_modelparams(self, **mparams):
        self.modelparams.update(mparams)

    def compute_data(self, h, vp, vs, rho, **params):
        """
        Method to compute the synthetic data. Here you probably need to
        include your quick e.g. fortran written code.
        """
        test = self.modelparams['test']

        z = np.cumsum(h)
        z = np.concatenate(([0], z[:-1]))

        xmod, ymod = quick_routine(test, z, vp, vs, rho)

    def validate(self, xmod, ymod):
        """Some condition that modeled data is valid. """
        if ymod.size == self.obsx.size:
            # xmod == xobs !!!
            return xmod, ymod
        else:
            return np.nan, np.nan

    def run_model(self, h, vp, vs, rho, **params):
        # incoming model is float32
        xmod, ymod = self.compute_rf(h, vp, vs, rho, **params)

        return self.validate(xmod, ymod)
