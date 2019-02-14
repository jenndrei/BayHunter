# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import os
# set os.environment variables to ensure that numerical computations
# do not do multiprocessing !! Essential !! Do not change !
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import os.path as op
import matplotlib
matplotlib.use('PDF')

from BayHunter import PlotFromStorage
from BayHunter import Targets
from BayHunter import utils
from BayHunter import MCMC_Optimizer
from BayHunter import ModelMatrix
from BayHunter import SynthObs
import logging


#
# console printout formatting
#
formatter = ' %(processName)-12s: %(levelname)-8s |  %(message)s'
logging.basicConfig(format=formatter, level=logging.INFO)
logger = logging.getLogger()


#
# ------------------------------------------------------------  obs SYNTH DATA
#
# Load priors and initparams from config.ini or simply create dictionaries.
initfile = 'config.ini'
priors, initparams = utils.load_params(initfile)

# Load observed data (synthetic test data)
xsw, _ysw = np.loadtxt('observed/st3_rdispph.dat', usecols=[0, 1]).T
xrf, _yrf = np.loadtxt('observed/st3_prf.dat', usecols=[0, 1]).T

# add noise to create observed data
# order of noise values (correlation, amplitude):
# noise = [corr1, sigma1, corr2, sigma2] for 2 targets
noise = [0.0, 0.012, 0.90, 0.005]
ysw_err = SynthObs.compute_expnoise(_ysw, corr=noise[0], sigma=noise[1])
ysw = _ysw + ysw_err
yrf_err = SynthObs.compute_gaussnoise(_yrf, corr=noise[2], sigma=noise[3])
yrf = _yrf + yrf_err


#
# -------------------------------------------  get refernece model for BayWatch
#
# Create truemodel only if you wish to have reference values in plots
# and BayWatch. You ONLY need to assign the values in truemodel that you
# wish to have visible.
dep, vs = np.loadtxt('observed/st3_mod.dat', usecols=[0, 2], skiprows=1).T
pdep = np.concatenate((np.repeat(dep, 2)[1:], [150]))
pvs = np.repeat(vs, 2)

truenoise = np.concatenate(([noise[0]], [np.std(ysw_err)],   # target 1
                            [noise[2]], [np.std(yrf_err)]))  # target 2

explike = SynthObs.compute_explike(yobss=[ysw, yrf], ymods=[_ysw, _yrf],
                                   noise=truenoise, gauss=[False, True],
                                   rcond=initparams['rcond'])
truemodel = {'model': (pdep, pvs),
             'nlays': 3,
             'noise': truenoise,
             'explike': explike,
             }

print truenoise, explike


#
#  -----------------------------------------------------------  DEFINE TARGETS
#
# Only pass x and y observed data to the Targets object which is matching
# the data type. You can chose for SWD any combination of Rayleigh, Love, group
# and phase velocity. Default is the fundamendal mode, but this can be updated.
# For RF chose P or S. You can also use user defined targets or replace the
# forward modeling plugin wih your own module.
target1 = Targets.RayleighDispersionPhase(xsw, ysw)
target2 = Targets.PReceiverFunction(xrf, yrf)
target2.moddata.plugin.set_modelparams(gauss=1., water=0.01, p=6.4)

# Join the targets. targets must be a list instance with all targets
# you want to use for MCMC Bayesian inversion.
targets = Targets.JointTarget(targets=[target1, target2])


#
#  ---------------------------------------------------  Quick parameter update
#
# "priors" and "initparams" from config.ini are python dictionaries. You could
# also simply define the dictionaries directly in the script, if you don't want
# to use a config.ini file. Or update the dictionaries as follows, e.g. if you
# have station specific values, etc.
# See docs/bayhunter.pdf for explanation of parameters

priors.update({'mohoest': (38, 4),  # optional, moho estimate (mean, std)
               'rfnoise_corr': 0.9,
               'swdnoise_corr': 0.
               # 'rfnoise_sigma': np.std(yrf_err),  # fixed to true value
               # 'swdnoise_sigma': np.std(ysw_err),  # fixed to true value
               })

initparams.update({'nchains': 5,
                   'iter_burnin': (2048 * 32),
                   'iter_main': (2048 * 16),
                   'propdist': (0.025, 0.025, 0.015, 0.005),
                   })


#
#  -------------------------------------------------------  MCMC BAY INVERSION
#
# Save configfile for baywatch. refmodel must not be defined.
utils.save_baywatch_config(targets, path='.', priors=priors,
                           initparams=initparams, refmodel=truemodel)
optimizer = MCMC_Optimizer(targets, initparams=initparams, priors=priors,
                           random_seed=None)
# default for the number of threads is the amount of cpus == one chain per cpu.
# if baywatch is True, inversion data is continuously send out (dtsend)
# to be received by BayWatch (see below).
optimizer.mp_inversion(nthreads=6, baywatch=True, dtsend=1)


#
# #  ---------------------------------------------- Model resaving and plotting
path = initparams['savepath']
cfile = '%s_config.pkl' % initparams['station']
configfile = op.join(path, 'data', cfile)
obj = PlotFromStorage(configfile)
# The final distributions will be saved with save_final_distribution.
# Beforehand, outlier chains will be detected and excluded.
# Outlier chains are defined as chains with a likelihood deviation
# of dev * 100 % from the median posterior likelihood of the best chain.
obj.save_final_distribution(maxmodels=100000, dev=0.05)
# Save a selection of important plots
obj.save_plots(refmodel=truemodel)
obj.merge_pdfs()

#
# If you are only interested on the mean posterior velocity model, type:
# file = op.join(initparams['savepath'], 'data/c_models.npy')
# models = np.load(file)
# singlemodels = ModelMatrix.get_singlemodels(models)
# vs, dep = singlemodels['mean']

#
# #  ---------------------------------------------- WATCH YOUR INVERSION
# if you want to use BayWatch, simply type "baywatch ." in the terminal in the
# folder you saved your baywatch configfile or type the full path instead
# of ".". Type "baywatch --help" for further options.

# if you give your public address as option (default is local address of PC),
# you can also use BayWatch via VPN from 'outside'.
# address = '139.?.?.?'  # here your complete address !!!
