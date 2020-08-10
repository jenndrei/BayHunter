# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import copy
import time
import numpy as np
import os.path as op

from BayHunter import Model, ModelMatrix
from BayHunter import utils

import logging
logger = logging.getLogger()


PAR_MAP = {'vsmod': 0, 'zvmod': 1, 'birth': 2, 'death': 2,
           'noise': 3, 'vpvs': 4}


class SingleChain(object):

    def __init__(self, targets, chainidx=0, initparams={}, modelpriors={},
                 sharedmodels=None, sharedmisfits=None, sharedlikes=None,
                 sharednoise=None, sharedvpvs=None, random_seed=None):
        self.chainidx = chainidx
        self.rstate = np.random.RandomState(random_seed)

        defaults = utils.get_path('defaults.ini')
        self.priors, self.initparams = utils.load_params(defaults)
        self.initparams.update(initparams)
        self.priors.update(modelpriors)
        self.dv = (self.priors['vs'][1] - self.priors['vs'][0])

        self.nchains = self.initparams['nchains']
        self.station = self.initparams['station']

        # set targets and inversion specific parameters
        self.targets = targets

        # set parameters
        self.iter_phase1 = int(self.initparams['iter_burnin'])
        self.iter_phase2 = int(self.initparams['iter_main'])
        self.iterations = self.iter_phase1 + self.iter_phase2
        self.iiter = -self.iter_phase1
        self.lastmoditer = self.iiter

        self.propdist = np.array(self.initparams['propdist'])
        self.acceptance = self.initparams['acceptance']
        self.thickmin = self.initparams['thickmin']
        self.maxlayers = int(self.priors['layers'][1]) + 1

        self.lowvelperc = self.initparams['lvz']
        self.highvelperc = self.initparams['hvz']
        self.mantle = self.priors['mantle']

        # chain models
        self._init_chainarrays(sharedmodels, sharedmisfits, sharedlikes,
                               sharednoise, sharedvpvs)

        # init model and values
        self._init_model_and_currentvalues()


# init model and misfit / likelihood

    def _init_model_and_currentvalues(self):
        ivpvs = self.draw_initvpvs()
        self.currentvpvs = ivpvs
        imodel = self.draw_initmodel()
        # self.currentmodel = imodel
        inoise, corrfix = self.draw_initnoiseparams()
        # self.currentnoise = inoise

        rcond = self.initparams['rcond']
        self.set_target_covariance(corrfix[::2], inoise[::2], rcond)

        vp, vs, h = Model.get_vp_vs_h(imodel, ivpvs, self.mantle)
        self.targets.evaluate(h=h, vp=vp, vs=vs, noise=inoise)

        # self.currentmisfits = self.targets.proposalmisfits
        # self.currentlikelihood = self.targets.proposallikelihood

        logger.debug((vs, h))

        self.n = 0  # accepted models counter
        self.accept_as_currentmodel(imodel, inoise, ivpvs)
        self.append_currentmodel()

    def draw_initmodel(self):
        keys = self.priors.keys()
        zmin, zmax = self.priors['z']
        vsmin, vsmax = self.priors['vs']
        layers = self.priors['layers'][0] + 1  # half space

        vs = self.rstate.uniform(low=vsmin, high=vsmax, size=layers)
        vs.sort()

        if (self.priors['mohoest'] is not None and layers > 1):
            mean, std = self.priors['mohoest']
            moho = self.rstate.normal(loc=mean, scale=std)
            tmp_z = self.rstate.uniform(1, np.min([5, moho]))  # 1-5
            tmp_z_vnoi = [moho-tmp_z, moho+tmp_z]

            if (layers - 2) == 0:
                z_vnoi = tmp_z_vnoi
            else:
                z_vnoi = np.concatenate((
                    tmp_z_vnoi,
                    self.rstate.uniform(low=zmin, high=zmax, size=(layers - 2))))

        else:  # no moho estimate
            z_vnoi = self.rstate.uniform(low=zmin, high=zmax, size=layers)

        z_vnoi.sort()
        model = np.concatenate((vs, z_vnoi))

        return(model if self._validmodel(model)
               else self.draw_initmodel())

    def draw_initnoiseparams(self):
        # for each target the noiseparams are (corr and sigma)
        noiserefs = ['noise_corr', 'noise_sigma']
        init_noise = np.ones(len(self.targets.targets)*2) * np.nan
        corrfix = np.zeros(len(self.targets.targets)*2, dtype=bool)

        self.noisepriors = []
        for i, target in enumerate(self.targets.targets):
            for j, noiseref in enumerate(noiserefs):
                idx = (2*i)+j
                noiseprior = self.priors[target.noiseref + noiseref]

                if type(noiseprior) in [int, float, np.float64]:
                    corrfix[idx] = True
                    init_noise[idx] = noiseprior
                else:
                    init_noise[idx] = self.rstate.uniform(
                        low=noiseprior[0], high=noiseprior[1])

                self.noisepriors.append(noiseprior)

        self.noiseinds = np.where(corrfix == 0)[0]
        if len(self.noiseinds) == 0:
            logger.warning('All your noise parameters are fixed. On Purpose?')

        return init_noise, corrfix

    def draw_initvpvs(self):
        if type(self.priors['vpvs']) == np.float:
            return self.priors['vpvs']

        vpvsmin, vpvsmax = self.priors['vpvs']
        return self.rstate.uniform(low=vpvsmin, high=vpvsmax)

    def set_target_covariance(self, corrfix, noise_corr, rcond=None):
        # SWD noise hyper-parameters: if corr is not 0, the correlation of data
        # points assumed will be exponential.
        # RF noise hyper-parameters: if corr is not 0, but fixed, the
        # correlation between data points will be assumed gaussian (realistic).
        # if the prior for RFcorr is a range, the computation switches
        # to exponential correlated noise for RF, as gaussian noise computation
        # is too time expensive because of computation of inverse and
        # determinant each time _corr is perturbed

        for i, target in enumerate(self.targets.targets):
            target_corrfix = corrfix[i]
            target_noise_corr = noise_corr[i]

            if not target_corrfix:
                # exponential for each target
                target.get_covariance = target.valuation.get_covariance_exp
                continue

            if (target_noise_corr == 0 and np.any(np.isnan(target.obsdata.yerr))):
                # diagonal for each target, corr inrelevant for likelihood, rel error
                target.get_covariance = target.valuation.get_covariance_nocorr
                continue

            elif target_noise_corr == 0:
                # diagonal for each target, corr inrelevant for likelihood
                target.get_covariance = target.valuation.get_covariance_nocorr_scalederr
                continue

            # gauss for RF
            if target.noiseref == 'rf':
                size = target.obsdata.x.size
                target.valuation.init_covariance_gauss(
                    target_noise_corr, size, rcond=rcond)
                target.get_covariance = target.valuation.get_covariance_gauss

            # exp for noise_corr
            elif target.noiseref == 'swd':
                target.get_covariance = target.valuation.get_covariance_exp

            else:
                message = 'The noise correlation automatically defaults to the \
exponential law. Explicitly state a noise reference for your user target \
(target.noiseref) if wished differently.'
                logger.info(message)
                target.noiseref == 'swd'
                target.get_covariance = target.valuation.get_covariance_exp

    def _init_chainarrays(self, sharedmodels, sharedmisfits, sharedlikes,
                          sharednoise, sharedvpvs):
        """from shared arrays"""
        ntargets = self.targets.ntargets
        chainidx = self.chainidx
        nchains = self.nchains

        accepted_models = int(self.iterations * np.max(self.acceptance) / 100.)
        self.nmodels = accepted_models  # 'iterations'

        msize = self.nmodels * self.maxlayers * 2
        nsize = self.nmodels * ntargets * 2
        missize = self.nmodels * (ntargets + 1)
        dtype = np.float32

        models = np.frombuffer(sharedmodels, dtype=dtype).\
            reshape((nchains, msize))
        misfits = np.frombuffer(sharedmisfits, dtype=dtype).\
            reshape((nchains, missize))
        likes = np.frombuffer(sharedlikes, dtype=dtype).\
            reshape((nchains, self.nmodels))
        noise = np.frombuffer(sharednoise, dtype=dtype).\
            reshape((nchains, nsize))
        vpvs = np.frombuffer(sharedvpvs, dtype=dtype).\
            reshape((nchains, self.nmodels))

        self.chainmodels = models[chainidx].reshape(
            self.nmodels, self.maxlayers*2)
        self.chainmisfits = misfits[chainidx].reshape(
            self.nmodels, ntargets+1)
        self.chainlikes = likes[chainidx]
        self.chainnoise = noise[chainidx].reshape(
            self.nmodels, ntargets*2)
        self.chainvpvs = vpvs[chainidx]
        self.chainiter = np.ones(self.chainlikes.size) * np.nan


# update current model (change layer number and values)

    def _model_layerbirth(self, model):
        """
        Draw a random voronoi nucleus depth from z and assign a new Vs.

        The new Vs is based on the before Vs value at the drawn z_vnoi
        position (self.propdist[2]).
        """
        n, vs_vnoi, z_vnoi = Model.split_modelparams(model)

        # new voronoi depth
        zmin, zmax = self.priors['z']
        z_birth = self.rstate.uniform(low=zmin, high=zmax)

        ind = np.argmin((abs(z_vnoi - z_birth)))  # closest z_vnoi
        vs_before = vs_vnoi[ind]
        vs_birth = vs_before + self.rstate.normal(0, self.propdist[2])

        z_new = np.concatenate((z_vnoi, [z_birth]))
        vs_new = np.concatenate((vs_vnoi, [vs_birth]))

        self.dvs2 = np.square(vs_birth - vs_before)
        return np.concatenate((vs_new, z_new))

    def _model_layerdeath(self, model):
        """
        Remove a random voronoi nucleus depth from model. Delete corresponding
        Vs from model.
        """
        n, vs_vnoi, z_vnoi = Model.split_modelparams(model)
        ind_death = self.rstate.randint(low=0, high=(z_vnoi.size))
        z_before = z_vnoi[ind_death]
        vs_before = vs_vnoi[ind_death]

        z_new = np.delete(z_vnoi, ind_death)
        vs_new = np.delete(vs_vnoi, ind_death)

        ind = np.argmin((abs(z_new - z_before)))
        vs_after = vs_new[ind]
        self.dvs2 = np.square(vs_after - vs_before)
        return np.concatenate((vs_new, z_new))

    def _model_vschange(self, model):
        """Randomly chose a layer to change Vs with Gauss distribution."""
        ind = self.rstate.randint(0, model.size / 2)
        vs_mod = self.rstate.normal(0, self.propdist[0])
        model[ind] = model[ind] + vs_mod
        return model

    def _model_zvnoi_move(self, model):
        """Randomly chose a layer to change z_vnoi with Gauss distribution."""
        ind = self.rstate.randint(model.size / 2, model.size)
        z_mod = self.rstate.normal(0, self.propdist[1])
        model[ind] = model[ind] + z_mod
        return model

    def _get_modelproposal(self, modify):
        model = copy.copy(self.currentmodel)

        if modify == 'vsmod':
            propmodel = self._model_vschange(model)
        elif modify == 'zvmod':
            propmodel = self._model_zvnoi_move(model)
        elif modify == 'birth':
            propmodel = self._model_layerbirth(model)
        elif modify == 'death':
            propmodel = self._model_layerdeath(model)

        return self._sort_modelproposal(propmodel)

    def _sort_modelproposal(self, model):
        """
        Return the sorted proposal model.

        This method is necessary, if the z_vnoi from the new proposal model
        are not ordered, i.e. if one z_vnoi value is added or strongly modified.
        """
        n, vs, z_vnoi = Model.split_modelparams(model)
        if np.all(np.diff(z_vnoi) > 0):   # monotone increasing
            return model
        else:
            ind = np.argsort(z_vnoi)
            model_sort = np.concatenate((vs[ind], z_vnoi[ind]))
        return model_sort

    def _validmodel(self, model):
        """
        Check model before the forward modeling.

        - The model must contain all values > 0.
        - The layer thicknesses must be at least thickmin km.
        - if lvz: low velocity zones are allowed with the deeper layer velocity
           no smaller than (1-perc) * velocity of layer above.
        - ... and some other constraints. E.g. vs boundaries (prior) given.
        """
        vp, vs, h = Model.get_vp_vs_h(model, self.currentvpvs, self.mantle)

        # check whether nlayers lies within the prior
        layermin = self.priors['layers'][0]
        layermax = self.priors['layers'][1]
        layermodel = (h.size - 1)
        if not (layermodel >= layermin and layermodel <= layermax):
            logger.debug("chain%d: model- nlayers not in prior"
                         % self.chainidx)
            return False

        # check model for layers with thicknesses of smaller thickmin
        if np.any(h[:-1] < self.thickmin):
            logger.debug("chain%d: thicknesses are not larger than thickmin"
                         % self.chainidx)
            return False

        # check whether vs lies within the prior
        vsmin = self.priors['vs'][0]
        vsmax = self.priors['vs'][1]
        if np.any(vs < vsmin) or np.any(vs > vsmax):
            logger.debug("chain%d: model- vs not in prior"
                         % self.chainidx)
            return False

        # check whether interfaces lie within prior
        zmin = self.priors['z'][0]
        zmax = self.priors['z'][1]
        z = np.cumsum(h)
        if np.any(z < zmin) or np.any(z > zmax):
            logger.debug("chain%d: model- z not in prior"
                         % self.chainidx)
            return False

        if self.lowvelperc is not None:
            # check model for low velocity zones. If larger than perc, then
            # compvels must be positive
            compvels = vs[1:] - (vs[:-1] * (1 - self.lowvelperc))
            if not compvels.size == compvels[compvels > 0].size:
                logger.debug("chain%d: low velocity zone issues"
                             % self.chainidx)
                return False

        if self.highvelperc is not None:
            # check model for high velocity zones. If larger than perc, then
            # compvels must be positive.
            compvels = (vs[:-1] * (1 + self.highvelperc)) - vs[1:]
            if not compvels.size == compvels[compvels > 0].size:
                logger.debug("chain%d: high velocity zone issues"
                             % self.chainidx)
                return False

        return True

    def _get_hyperparameter_proposal(self):
        noise = copy.copy(self.currentnoise)
        ind = self.rstate.choice(self.noiseinds)

        noise_mod = self.rstate.normal(0, self.propdist[3])
        noise[ind] = noise[ind] + noise_mod
        return noise

    def _validnoise(self, noise):
        for idx in self.noiseinds:
            if noise[idx] < self.noisepriors[idx][0] or \
                    noise[idx] > self.noisepriors[idx][1]:
                return False
        return True

    def _get_vpvs_proposal(self):
        vpvs = copy.copy(self.currentvpvs)
        vpvs_mod = self.rstate.normal(0, self.propdist[4])
        vpvs = vpvs + vpvs_mod
        return vpvs

    def _validvpvs(self, vpvs):
        # only works if vpvs-priors is a range
        if vpvs < self.priors['vpvs'][0] or \
                vpvs > self.priors['vpvs'][1]:
            return False
        return True


# accept / save current models

    def adjust_propdist(self):
        """
        Modify self.propdist to adjust acceptance rate of models to given
        percentace span: increase or decrease by five percent.
        """
        with np.errstate(invalid='ignore'):
            acceptrate = self.accepted / self.proposed * 100

        # minimum distribution width forced to be not less than 1 m/s, 1 m
        # actually only touched by vs distribution
        propdistmin = np.full(acceptrate.size, 0.001)

        for i, rate in enumerate(acceptrate):
            if np.isnan(rate):
                # only if not inverted for
                continue
            if rate < self.acceptance[0]:
                new = self.propdist[i] * 0.95
                if new < propdistmin[i]:
                    new = propdistmin[i]
                self.propdist[i] = new

            elif rate > self.acceptance[1]:
                self.propdist[i] = self.propdist[i] * 1.05
            else:
                pass

    def get_acceptance_probability(self, modify):
        """
        Acceptance probability will be computed dependent on the modification.

        Parametrization alteration (Vs or voronoi nuclei position)
            the acceptance probability is equal to likelihood ratio.

        Model dimension alteration (layer birth or death)
            the probability was computed after the formulation of Bodin et al.,
            2012: 'Transdimensional inversion of receiver functions and
            surface wave dispersion'.
        """
        if modify in ['vsmod', 'zvmod', 'noise', 'vpvs']:
            # only velocity or thickness changes are made
            # also used for noise changes
            alpha = self.targets.proposallikelihood - self.currentlikelihood

        elif modify in ['birth', ]:
            theta = self.propdist[2]  # Gaussian distribution
            # self.dvs2 = delta vs square = np.square(v'_(k+1) - v_(i))
            A = (theta * np.sqrt(2 * np.pi)) / self.dv
            B = self.dvs2 / (2. * np.square(theta))
            C = self.targets.proposallikelihood - self.currentlikelihood

            alpha = np.log(A) + B + C

        elif modify in ['death', ]:
            theta = self.propdist[2]  # Gaussian distribution
            # self.dvs2 = delta vs square = np.square(v'_(j) - v_(i))
            A = self.dv / (theta * np.sqrt(2 * np.pi))
            B = self.dvs2 / (2. * np.square(theta))
            C = self.targets.proposallikelihood - self.currentlikelihood

            alpha = np.log(A) - B + C

        return alpha

    def accept_as_currentmodel(self, model, noise, vpvs):
        """Assign currentmodel and currentvalues to self."""
        self.currentmisfits = self.targets.proposalmisfits
        self.currentlikelihood = self.targets.proposallikelihood
        self.currentmodel = model
        self.currentnoise = noise
        self.currentvpvs = vpvs
        self.lastmoditer = self.iiter

    def append_currentmodel(self):
        """Append currentmodel to chainmodels and values."""
        self.chainmodels[self.n, :self.currentmodel.size] = self.currentmodel
        self.chainmisfits[self.n, :] = self.currentmisfits
        self.chainlikes[self.n] = self.currentlikelihood
        self.chainnoise[self.n, :] = self.currentnoise
        self.chainvpvs[self.n] = self.currentvpvs

        self.chainiter[self.n] = self.iiter
        self.n += 1

# run optimization

    def iterate(self):
        if self.iiter < (-self.iter_phase1 + (self.iterations * 0.01)):
            # only allow vs and z modifications the first 1 % of iterations
            modify = self.rstate.choice(['vsmod', 'zvmod'] + self.noisemods +
                                        self.vpvsmods)
        else:
            modify = self.rstate.choice(self.modifications)

        if modify in self.modelmods:
            proposalmodel = self._get_modelproposal(modify)
            proposalnoise = self.currentnoise
            proposalvpvs = self.currentvpvs
            if not self._validmodel(proposalmodel):
                proposalmodel = None

        elif modify in self.noisemods:
            proposalmodel = self.currentmodel
            proposalnoise = self._get_hyperparameter_proposal()
            proposalvpvs = self.currentvpvs
            if not self._validnoise(proposalnoise):
                proposalmodel = None

        elif modify == 'vpvs':
            proposalmodel = self.currentmodel
            proposalnoise = self.currentnoise
            proposalvpvs = self._get_vpvs_proposal()
            if not self._validvpvs(proposalvpvs):
                proposalmodel = None

        if proposalmodel is None:
            # If not a valid proposal model and noise params are found,
            # leave self.iterate and try with another modification
            # should not occur often.
            logger.debug('Not able to find a proposal for %s' % modify)
            self.iiter += 1
            return

        # compute synthetic data and likelihood, misfit
        vp, vs, h = Model.get_vp_vs_h(proposalmodel, proposalvpvs, self.mantle)
        self.targets.evaluate(h=h, vp=vp, vs=vs, noise=proposalnoise)

        paridx = PAR_MAP[modify]
        self.proposed[paridx] += 1

        # Replace self.currentmodel with proposalmodel with acceptance
        # probability alpha. Accept candidate sample (proposalmodel)
        # with probability alpha, or reject it with probability (1 - alpha).
        # these are log values ! alpha is log.
        u = np.log(self.rstate.uniform(0, 1))
        alpha = self.get_acceptance_probability(modify)

        # #### _____________________________________________________________
        if u < alpha:
            # always the case if self.jointlike > self.bestlike (alpha>1)
            self.accept_as_currentmodel(proposalmodel, proposalnoise, proposalvpvs)
            self.append_currentmodel()
            self.accepted[paridx] += 1

        # print inversion status information
        if self.iiter % 5000 == 0:
            runtime = time.time() - self.tnull
            current_iterations = self.iiter + self.iter_phase1

            if current_iterations > 0:
                acceptrate = float(self.n) / current_iterations * 100.

                logger.info('%6d %5d + hs %8.3f\t%9d |%6.1f s  | %.1f ' % (
                    self.lastmoditer, self.currentmodel.size/2 - 1,
                    self.currentmisfits[-1], self.currentlikelihood,
                    runtime, acceptrate) + r'%')

            self.tnull = time.time()

        # stabilize model acceptance rate
        if self.iiter % 1000 == 0:
            if np.all(self.proposed) != 0:
                self.adjust_propdist()

        self.iiter += 1

    def run_chain(self):
        t0 = time.time()
        self.tnull = time.time()
        self.iiter = -self.iter_phase1

        self.modelmods = ['vsmod', 'zvmod', 'birth', 'death']
        self.noisemods = [] if len(self.noiseinds) == 0 else ['noise']
        self.vpvsmods = [] if type(self.priors['vpvs']) == np.float else ['vpvs']
        self.modifications = self.modelmods + self.noisemods + self.vpvsmods

        self.accepted = np.zeros(len(self.propdist))
        self.proposed = np.zeros(len(self.propdist))

        while self.iiter < self.iter_phase2:
            self.iterate()

        runtime = (time.time() - t0)

        # update chain values (eliminate nan rows)
        self.chainmodels = self.chainmodels[:self.n, :]
        self.chainmisfits = self.chainmisfits[:self.n, :]
        self.chainlikes = self.chainlikes[:self.n]
        self.chainnoise = self.chainnoise[:self.n, :]
        self.chainvpvs = self.chainvpvs[:self.n]
        self.chainiter = self.chainiter[:self.n]

        # only consider models after burnin phase
        p1ind = np.where(self.chainiter < 0)[0]
        p2ind = np.where(self.chainiter >= 0)[0]

        if p1ind.size != 0:
            wmodels, wlikes, wmisfits, wnoise, wvpvs = self.get_weightedvalues(
                pind=p1ind, finaliter=0)
            self.p1models = wmodels  # p1 = phase one
            self.p1misfits = wmisfits
            self.p1likes = wlikes
            self.p1noise = wnoise
            self.p1vpvs = wvpvs

        if p2ind.size != 0:
            wmodels, wlikes, wmisfits, wnoise, wvpvs = self.get_weightedvalues(
                pind=p2ind, finaliter=self.iiter)
            self.p2models = wmodels  # p2 = phase two
            self.p2misfits = wmisfits
            self.p2likes = wlikes
            self.p2noise = wnoise
            self.p2vpvs = wvpvs

        accmodels = float(self.p2likes.size)  # accepted models in p2 phase
        maxmodels = float(self.initparams['maxmodels'])  # for saving
        self.thinning = int(np.ceil(accmodels / maxmodels))
        self.save_finalmodels()

        logger.debug('time for inversion: %.2f s' % runtime)

    def get_weightedvalues(self, pind, finaliter):
        """
        Models will get repeated (weighted).

        Each iteration, if there was no model proposal accepted, the current
        model gets repeated once more. This weight is based on self.chainiter,
        which documents the iteration of the last accepted model."""
        pmodels = self.chainmodels[pind]  # p = phase (1 or 2)
        pmisfits = self.chainmisfits[pind]
        plikes = self.chainlikes[pind]
        pnoise = self.chainnoise[pind]
        pvpvs = self.chainvpvs[pind]
        pweights = np.diff(np.concatenate((self.chainiter[pind], [finaliter])))

        wmodels, wlikes, wmisfits, wnoise, wvpvs = ModelMatrix.get_weightedvalues(
            pweights, models=pmodels, likes=plikes, misfits=pmisfits,
            noiseparams=pnoise, vpvs=pvpvs)
        return wmodels, wlikes, wmisfits, wnoise, wvpvs

    def save_finalmodels(self):
        """Save chainmodels as pkl file"""
        savepath = op.join(self.initparams['savepath'], 'data')
        names = ['models', 'likes', 'misfits', 'noise', 'vpvs']

        # phase 1 -- burnin
        try:
            for i, data in enumerate([self.p1models, self.p1likes,
                                     self.p1misfits, self.p1noise,
                                     self.p1vpvs]):
                outfile = op.join(savepath, 'c%.3d_p1%s' % (self.chainidx, names[i]))
                np.save(outfile, data[::self.thinning])
        except:
            logger.info('No burnin models accepted.')

        # phase 2 -- main / posterior phase
        try:
            for i, data in enumerate([self.p2models, self.p2likes,
                                     self.p2misfits, self.p2noise,
                                     self.p2vpvs]):
                outfile = op.join(savepath, 'c%.3d_p2%s' % (self.chainidx, names[i]))
                np.save(outfile, data[::self.thinning])

            logger.info('> Saving %d models (main phase).' % len(data[::self.thinning]))
        except:
            logger.info('No main phase models accepted.')
