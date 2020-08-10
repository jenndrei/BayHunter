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

   # --------------------------------------------------- obs SYNTH DATA
   # Load observed data (synthetic test data)
   xsw, _ysw = np.loadtxt('rdispph.dat').T
   xrf, _yrf = np.loadtxt('prf.dat').T

   # Create noise and add to synthetic data
   ysw = _ysw + SynthObs.compute_expnoise(_ysw, corr=0, sigma=0.012)
   yrf = _yrf + SynthObs.compute_gaussnoise(_yrf, corr=0.98, sigma=0.005)

   #  --------------------------------------------------------- TARGETS
   # Assign data to target classes
   target1 = Targets.RayleighDispersionPhase(xsw, ysw)
   target2 = Targets.PReceiverFunction(xrf, yrf)
   target2.moddata.plugin.set_modelparams(gauss=1.0, water=0.01, p=6.4)

   # Join the targets
   targets = Targets.JointTarget(targets=[target1, target2])

   #  ------------------------------------------------------ PARAMETERS
   # Define parameters as dictionaries ...
   priors = {'vs': (2, 5),
             'layers': (1, 20),
             'vpvs': 1.73,
             'rfnoise_corr': 0.98,
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

   #  -------------------------------------------------- MCMC INVERSION
   # Save configfile for baywatch
   utils.save_baywatch_config(targets, priors=priors, initparams=initparams)
   optimizer = MCMC_Optimizer(targets, initparams=initparams, priors=priors,
                              random_seed=None)

   # start inversion, activate BayWatch
   optimizer.mp_inversion(nthreads=8, baywatch=True, dtsend=1)

   #  ----------------------------------------------- SAVING / PLOTTING
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
   
  
Estimation of :math:`r_{RF}`
----------------------------

::

   import numpy as np
   import matplotlib.pyplot as plt
   from BayHunter import utils

   # ------------------------------------------------ RF and parameters
   # load observed data (synthetic RF, Gauss factor a=1)
   rfx, rfy = np.loadtxt('observed/st3_prf.dat').T
   rfa = 1  # a

   # define parameters
   dt = 0.2
   draws = 40000
   rrfs = [0.75, 0.85, 0.95, 0.97, 0.98, 0.99]

   pars = {'rfx': rfx, 'rfy': rfy, 'rfa': rfa,
           'a': rfa, 'dt': dt, 'rrfs': rrfs,
           'draws': draws}

   # ------------------------------- visualize 'raw' data and estimates
   fig = utils.plot_rrf_estimate(pars=pars)
   fig.savefig('st3_rrf_estimate.pdf', bbox_inches='tight')


   # --------------------------- return values for costum visualization
   # update parameters...
   pars['rrfs'] = np.linspace(0.9, 0.999, 25)
   pars['draws'] = 2000

   # get rrf-values with corresponding a-values
   rrf, a = utils.rrf_estimate(pars=pars)

   # custom plot example
   fig, ax = plt.subplots()
   ax.plot(rrf, a, color='k', marker='x', ls='')
   ax.axhline(rfa, color='gray', label='reference')
   ax.set_xlabel('$r_{RF}$')
   ax.set_ylabel('Gauss factor a')
   ax.grid(color='lightgray')
   ax.legend(loc=1)
   fig.savefig('rrf-a_rel.pdf', bbox_inches='tight')

