.. role:: raw-latex(raw)
   :format: latex

Appendix
========

Inversion example
-----------------

::

   import numpy as np
   import os.path as op
   from BayHunter import utils
   from BayHunter import SynthObs
   from BayHunter import Targets
   from BayHunter import MCMC_Optimizer
   from BayHunter import PlotFromStorage

   # ---------------------------------------------------- obs SYNTH DATA
   # Load observed data (synthetic test data)
   xsw, _ysw = np.loadtxt('rdispph.dat').T
   xrf, _yrf = np.loadtxt('prf.dat').T

   # Create noise and add to synthetic data
   ysw = _ysw + SynthObs.compute_expnoise(_ysw, corr=0, sigma=0.012)
   yrf = _yrf + SynthObs.compute_gaussnoise(_yrf, corr=0.92, sigma=0.005)

   #  ---------------------------------------------------------- TARGETS
   # Assign data to target classes
   target1 = Targets.RayleighDispersionPhase(xsw, ysw)
   target2 = Targets.PReceiverFunction(xrf, yrf)
   target2.moddata.plugin.set_modelparams(gauss=1.0, water=0.01, p=6.4)

   # Join the targets
   targets = Targets.JointTarget(targets=[target1, target2])

   #  ------------------------------------------------------- PARAMETERS
   # Define parameters as dictionaries ...
   priors = {'vs': (2, 5),
             'layers': (1, 20),
             'vpvs': 1.73,
              'rfnoise_corr': 0.92,
              ...
             }

   initparams = {'nchains': 21,
          'iter_burnin': 100000,
          'iter_main': 50000,
          ...
                 }

   # ... or load from file
   initfile = 'config.ini'
   priors, initparams = utils.load_params(initfile)

   #  --------------------------------------------------- MCMC INVERSION
   # Save configfile for baywatch
   utils.save_baywatch_config(targets, priors=priors, initparams=initparams)
   optimizer = MCMC_Optimizer(targets, initparams=initparams, priors=priors,
                              random_seed=None)

   # start inversion, activate BayWatch
   optimizer.mp_inversion(nthreads=8, baywatch=True, dtsend=1)

   #  ------------------------------------------------ SAVING / PLOTTING
   # Initiate plotting object
   path = initparams['savepath']
   cfile = '%s_config.pkl' % initparams['station']
   configfile = op.join(path, 'data', cfile)
   obj = PlotFromStorage(configfile)

   # Save posterior distribution to combined files, incl. outlier detection
   obj.save_final_distribution(maxmodels=100000, dev=0.05)

   # Save a selection of important plots
   obj.save_plots()
   obj.merge_pdfs()
