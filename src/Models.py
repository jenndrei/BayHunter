# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import numpy as np
import copy


class Model(object):
    """Handle interpolating methods for a single model vector."""

    @staticmethod
    def split_modelparams(model):
        model = model[~np.isnan(model)]
        n = int(model.size / 2)  # layers

        vs = model[:n]
        z_vnoi = model[-n:]

        return n, vs, z_vnoi

    @staticmethod
    def get_vp(vs, vpvs=1.73, mantle=[4.3, 1.8]):
        """Return vp from vs, based on crustal and mantle vpvs."""
        ind_m = np.where((vs >= mantle[0]))[0]  # mantle

        vp = vs * vpvs  # correct for crust
        if len(ind_m) == 0:
            return vp
        else:
            ind_m[0] == np.int
            vp[ind_m[0]:] = vs[ind_m[0]:] * mantle[1]
        return vp

    @staticmethod
    def get_vp_vs_h(model, vpvs=1.73, mantle=None):
        """Return vp, vs and h from a input model [vs, z_vnoi]"""
        n, vs, z_vnoi = Model.split_modelparams(model)
        # discontinuities:
        z_disc = (z_vnoi[:n-1] + z_vnoi[1:n]) / 2.
        h_lay = (z_disc - np.concatenate(([0], z_disc[:-1])))
        h = np.concatenate((h_lay, [0]))

        if mantle is not None:
            vp = Model.get_vp(vs, vpvs, mantle)
        else:
            vp = vs * vpvs
        return vp, vs, h

    @staticmethod
    def get_stepmodel(model, vpvs=1.73, mantle=None):
        """Return a steplike model from input model, for plotting."""
        vp, vs, h = Model.get_vp_vs_h(model, vpvs, mantle)

        dep = np.cumsum(h)

        # insert steps into velocity model
        dep = np.concatenate([(d, d) for d in dep])
        dep_step = np.concatenate([[0], dep[:-1]])
        vp_step = np.concatenate([(v, v) for v in vp])
        vs_step = np.concatenate([(v, v) for v in vs])

        dep_step[-1] = np.max([150, dep_step[-1] * 2.5])  # half space

        return vp_step, vs_step, dep_step

    @staticmethod
    def get_stepmodel_from_h(h, vs, vpvs=1.73, dep=None, vp=None, mantle=None):
        """Return a steplike model from input model."""
        # insert steps into velocity model
        if dep is None:
            dep = np.cumsum(h)

        if vp is None:
            if mantle is not None:
                vp = Model.get_vp(vs, vpvs, mantle)
            else:
                vp = vs * vpvs

        dep = np.concatenate([(d, d) for d in dep])
        dep_step = np.concatenate([[0], dep[:-1]])
        vp_step = np.concatenate([(v, v) for v in vp])
        vs_step = np.concatenate([(v, v) for v in vs])

        dep_step[-1] = dep_step[-1] * 2.5  # half space

        return vp_step, vs_step, dep_step

    @staticmethod
    def get_interpmodel(model, dep_int, vpvs=1.73, mantle=None):
        """
        Return an interpolated stepmodel, for (histogram) plotting.

        Model is a vector of the parameters.
        """
        vp_step, vs_step, dep_step = Model.get_stepmodel(model, vpvs, mantle)
        vs_int = np.interp(dep_int, dep_step, vs_step)
        vp_int = np.interp(dep_int, dep_step, vp_step)

        return vp_int, vs_int


class ModelMatrix(object):
    """
    Handle interpolating methods for a collection of single models.

    Same as the Model class, but for a matrix. Only for Plotting
    or after inversion.
    """

    @staticmethod
    def _delete_nanmodels(models):
        """Remove nan models from model-matrix."""
        cmodels = copy.copy(models)
        mean = np.nanmean(cmodels, axis=1)
        nanidx = np.where((np.isnan(mean)))[0]

        if nanidx.size == 0:
            return cmodels
        else:
            return np.delete(cmodels, nanidx, axis=0)

    @staticmethod
    def _replace_zvnoi_h(models):
        """
        Return model matrix with (vs, h) - models.

        Each model in the matrix is parametrized with (vs, z_vnoi).
        For plotting, h will be computed from z_vnoi."""
        models = ModelMatrix._delete_nanmodels(models)

        for i, model in enumerate(models):
            _, vs, h = Model.get_vp_vs_h(model)
            newmodel = np.concatenate((vs, h))
            models[i][:newmodel.size] = newmodel
        return models

    @staticmethod
    def get_interpmodels(models, dep_int):
        """Return model matrix with interpolated stepmodels.

        Each model in the matrix is parametrized with (vs, z_vnoi)."""
        models = ModelMatrix._delete_nanmodels(models)

        deps_int = np.repeat([dep_int], len(models), axis=0)
        vss_int = np.empty((len(models), dep_int.size))

        for i, model in enumerate(models):
            # for vs, dep 2D histogram
            _, vs_int = Model.get_interpmodel(model, dep_int)
            vss_int[i] = vs_int

        return vss_int, deps_int

    @staticmethod
    def get_singlemodels(models, dep_int=None, misfits=None):
        """Return specific single models from model matrix (vs, depth).
        The model is a step model for plotting.

        -- interpolated
        (1) mean
        (2) median
        (3) minmax
        (4) stdminmax

        -- binned, vs step: 0.025 km/s
                   dep step: 0.5 km or as in dep_int
        (5) mode (histogram)

        -- not interpolated
        (6) bestmisfit   - min misfit
        """
        singlemodels = dict()

        if dep_int is None:
            # interpolate depth to 0.5 km bins.
            dep_int = np.linspace(0, 100, 201)

        vss_int, deps_int = ModelMatrix.get_interpmodels(models, dep_int)

        # (1) mean, (2) median
        mean = np.mean(vss_int, axis=0)
        median = np.median(vss_int, axis=0)

        # (3) minmax
        minmax = np.array((np.min(vss_int, axis=0), np.max(vss_int, axis=0))).T

        # (4) stdminmax
        stdmodel = np.std(vss_int, axis=0)
        stdminmodel = mean - stdmodel
        stdmaxmodel = mean + stdmodel

        stdminmax = np.array((stdminmodel, stdmaxmodel)).T

        # (5) mode from histogram
        vss_flatten = vss_int.flatten()
        vsbins = int((vss_flatten.max() - vss_flatten.min()) / 0.025)
        # in PlotFromStorage posterior_models2d
        data = np.histogram2d(vss_int.flatten(), deps_int.flatten(),
                              bins=(vsbins, dep_int))
        bins, vs_bin, dep_bin = np.array(data).T
        vs_center = (vs_bin[:-1] + vs_bin[1:]) / 2.
        dep_center = (dep_bin[:-1] + dep_bin[1:]) / 2.
        vs_mode = vs_center[np.argmax(bins.T, axis=1)]
        mode = (vs_mode, dep_center)

        # (6) bestmisfit - min misfit
        if misfits is not None:
            ind = np.argmin(misfits)
            _, vs_best, dep_best = Model.get_stepmodel(models[ind])

            singlemodels['minmisfit'] = (vs_best, dep_best)

        # add models to dictionary
        singlemodels['mean'] = (mean, dep_int)
        singlemodels['median'] = (median, dep_int)
        singlemodels['minmax'] = (minmax.T, dep_int)
        singlemodels['stdminmax'] = (stdminmax.T, dep_int)
        singlemodels['mode'] = mode

        return singlemodels

    @staticmethod
    def get_weightedvalues(weights, models=None, likes=None, misfits=None,
                           noiseparams=None, vpvs=None):
        """
        Return weighted matrix of models, misfits and noiseparams, and weighted
        vectors of likelihoods.

        Basically just repeats values, as given by weights.
        """
        weights = np.array(weights, dtype=int)
        wlikes, wmisfits, wmodels, wnoise, wvpvs = (None, None, None, None, None)

        if likes is not None:
            wlikes = np.repeat(likes, weights)

        if misfits is not None:
            if type(misfits[0]) in [int, float, np.float64]:
                wmisfits = np.repeat(misfits, weights)
            else:
                wmisfits = np.ones((np.sum(weights), misfits[0].size)) * np.nan
                n = 0
                for i, misfit in enumerate(misfits):
                    for rep in range(weights[i]):
                        wmisfits[n] = misfit
                        n += 1

        if models is not None:
            wmodels = np.ones((np.sum(weights), models[0].size)) * np.nan

            n = 0
            for i, model in enumerate(models):
                for rep in range(weights[i]):
                    wmodels[n] = model
                    n += 1

        if noiseparams is not None:
            wnoise = np.ones((np.sum(weights), noiseparams[0].size)) * np.nan

            n = 0
            for i, noisepars in enumerate(noiseparams):
                for rep in range(weights[i]):
                    wnoise[n] = noisepars
                    n += 1

        if vpvs is not None:
            wvpvs = np.repeat(vpvs, weights)

        return wmodels, wlikes, wmisfits, wnoise, wvpvs
