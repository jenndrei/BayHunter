.. role:: raw-latex(raw)
   :format: latex

.. _sec:tutorial:

Tutorial
========

This chapter contains the installation instructions of BayHunter,
followed by an example of how to set up and run an inversion. A
minimalistic working example is shown in the :doc:`Appendix <appendix>`. Furthermore, results from an inversion using synthetic data are shown and discussed. Be here also referred to the `full tutorial <https://github.com/jenndrei/BayHunter/tree/master/tutorial>`_ including data files.


Requirements and installation
-----------------------------

BayHunter is currently for a Python 2 environment (as of October, 2019). After installation of the required Python-modules, simply type the following to install BayHunter:

::

    sudo python setup.py install

============== ==================================
``numpy``      numerical computations
``matplotlib`` plotting library
``pyPdf``      merging PDFs
``configobj``  configuration file
``zmq``        BayWatch, inversion live-streaming
``Cython``     C-extensions for Python
============== ==================================

The forward modeling codes for surface wave dispersion and receiver functions are already included in the BayHunter package and will be compiled when installing BayHunter. BayHunter uses a Python wrapper interfacing the
``SURF96`` routine from :doc:`Herrmann and Ammon, 2002 <references>`, and  ``rfmini`` developed for BayHunter by
`Joachim Saul, GFZ <https://www.gfz-potsdam.de/en/staff/joachim-saul/>`_.

.. _sec:baysetup:

Setting up and running an inversion
-----------------------------------

Setting up the targets
~~~~~~~~~~~~~~~~~~~~~~~

As mentioned in :ref:`sec:intarg`, BayHunter provides six target classes (four SWD and two RF),
which use two types of forward modeling plugins (``SURF96``,
``rfmini``). For both targets, the user may update the default forward
modeling parameters with *set_modelparams* (see :doc:`Appendix <appendix>`). Parameters and default values are given in :numref:`Table {number} <tab:targetpars>`.

.. container::
   :name: tab:targetpars

   .. table:: Default forward modeling parameters for SWD and RF.

      +-----+----------------+----------------------------------------------+
      | SWD | mode = 1       | 1=fundamental mode, 2=1st higher mode, etc.  |
      +-----+----------------+----------------------------------------------+
      | RF  | gauss = 1.0    | Gauss factor, low pass filter                |
      +-----+----------------+----------------------------------------------+
      |     | water = 0.001  | water level stabilization                    |
      +-----+----------------+----------------------------------------------+
      |     | p = 6.4        | slowness in deg/s                            |
      +-----+----------------+----------------------------------------------+
      |     | nsv = ``None`` | near surface velocity in km/s for computation|
      |     |                | of incident angle (trace rotation). If       |
      |     |                | ``None``, nsv is taken from velocity-model.  |
      +-----+----------------+----------------------------------------------+

If the user wants to implement own forward modeling code, a new forward
modeling class for it is needed. After normally initializing a target
with BayHunter, an instance of the new forward modeling class must be
initialized and passed to the *update_plugin* method of the target. If
an additional data set is wished to be included in the inversion, i.e.,
from a non pre-defined target class, a new target class needs to be
implemented, additionally to the forward modeling class that handles the
synthetic data computation. For both, the forward modeling class and the
new target class, a template is stored on the GitHub repository. It is
important that the classes implement specifically named methods and
parameters to ensure the correct interface with BayHunter.

.. _sec:parsetup:

Setting up parameters
~~~~~~~~~~~~~~~~~~~~~~

Each chain will be initialized with the targets and with parameter
dictionaries. The model priors and inversion parameters that need to be
defined are listed with default values in :numref:`Table {number} <tab:invpars>`,
and are explained below in detail.

.. container::
   :name: tab:invpars

   .. table:: Default model priors and inversion parameters (SI-units, i.e., km, km/s, %). Model prior tuples define the limits (min, max) of a uniform distribution. ``None`` implies that the constraint is not used. Abbreviations and constraints are explained in the text.
   
      +---------------+-----------------------+------------+---------------+
      | vs            | = (1, 5)              | nchains    | = 3           |
      +---------------+-----------------------+------------+---------------+
      | z             | = (0, 60)             | :math:`i   | = 4096        |
      |               |                       | ter_       |               |
      |               |                       | {burnin}`  |               |
      +---------------+-----------------------+------------+---------------+
      | layers        | = (1, 20)             | :math:`i   | = 2048        |
      |               |                       | te         |               |
      |               |                       | r_{main}`  |               |
      +---------------+-----------------------+------------+---------------+
      | vpvs          | = (1.5, 2.1)          | acceptance | = (40,        |
      |               |                       |            | 45)           |
      +---------------+-----------------------+------------+---------------+
      | mantle        | = ``None``            | propdist   | = (0.015,     |
      | :math:`^1`    |                       | :math:`^3` | 0.015, 0.005, |
      |               |                       |            | 0.015, 0.005) |
      | mohoest       | = ``None``            |            |               |
      | :math:`^2`    |                       |            |               |
      +---------------+-----------------------+------------+---------------+
      | :math:`r_     | = (0.35, 0.75)        | thickmin   | = 0.          |
      | {RF}`         |                       |            |               |
      +---------------+-----------------------+------------+---------------+
      | :math:`\sigma | = (1e-5, 0.05)        | lvz        | = ``None``    |
      | _{RF}`        |                       |            |               |
      +---------------+-----------------------+------------+---------------+
      | :math:`r_     | = 0.                  | hvz        | = ``None``    |
      | {SWD}`        |                       |            |               |
      +---------------+-----------------------+------------+---------------+
      | :math:`\sigma | = (1e-5, 0.1)         | rcond      | = ``None``    |
      | _{SWD}`       |                       |            |               |
      +---------------+-----------------------+------------+---------------+
      | :math:`^1` i.e., (vs\ :math:`_m`,     | station    | = ’test’      |
      | vpvs\ :math:`_m`), e.g., (4.2, 1.8)   |            |               | 
      +---------------+-----------------------+------------+---------------+
      | :math:`^2` i.e., (z\ :math:`_{mean}`, | savepath   | = ’results/’  |
      | z\ :math:`_{std}`), e.g., (40, 4)     |            |               |
      +---------------+-----------------------+------------+---------------+
      | :math:`^3` i.e.,                      | maxmodels  | = 50 000      |
      | (vs, z\ :math:`_{move}`,              |            |               |
      | vs\ :math:`_{birth/death}`,           |            |               |
      | noise, vpvs)                          |            |               |
      +---------------+-----------------------+------------+---------------+

The priors for the velocity-depth structure include :math:`\mathrm{V_S}`
and depth, the number of layers, and average crustal
:math:`\mathrm{V_P}`/:math:`\mathrm{V_S}`. The ranges as given in
:numref:`Table {number} <tab:invpars>` indicate the bounds of uniform distributions.
:math:`\mathrm{V_P}`/:math:`\mathrm{V_S}` can also be given as a float
digit (e.g., 1.73), indicating a constant value during the inversion.
The parameter :math:`layers` does not include the underlying half space,
which is always added to the model. A mantle condition (vs\ :math:`_m`,
vpvs\ :math:`_m`), i.e., a vs\ :math:`_m` threshold beyond which
:math:`\mathrm{V_P}` is computed from vpvs\ :math:`_m`, can be chosen if
appropriate. There is also the option to give a single interface depth
estimate through :math:`mohoest`. It can be any interface, but the
initial idea behind was to give a Moho estimate. As explained in :ref:`sec:inmod`, this should only be
considered for testing purposes. Each noise scaling parameter
(:math:`r`, :math:`\sigma`) can be given by a range or a digit,
corresponding to the bounds of a uniform distribution (the parameter is
inverted for) or a constant value (unaltered during the inversion),
repectively.

For surface waves, the exponential correlation law
(Eq. :eq:`exp`) is a realistic estimate of the correlation
between data points and is automatically applied. For receiver
functions, the assumed correlation law should be Gaussian
(Eq. :eq:`gauss`), if the RFs are computed using a
Gaussian filter, and exponential, if the RFs are computed applying an
exponential filter. The inversion for :math:`r_{RF}` is viable for the
latter, however, not for the Gaussian correlation law as of
computational reasons (see :ref:`sec:complike`). Only if :math:`r_{RF}` is estimated by giving a
single digit, the Gaussian correlation law is considered. Otherwise, if
given a range for :math:`r_{RF}`, the exponential correlation law is
used. Note that the estimation of :math:`r_{RF}` using the exponential
law during an inversion, may not lead to correct results if the input RF
was Gaussian filtered.

Nevertheless, :math:`r_{RF}` can be estimated, as it is dependent on the
sampling rate and the applied Gaussian filter width. :numref:`Figure {number} <fig:rrf_est>` shows an application of the BayHunter implemented tool to estimate :math:`r_{RF}`. You will find a minimalistic code example in the :doc:`Appendix <appendix>` and an executable file with plenty of comments in the tutorial folder of the repository.

.. _fig:rrf_est:

.. figure :: _static/st3_rrf_est.png

	Visual estimation of :math:`r_{RF}`. Top: Synthetic receiver function from a 3-layer crustal model, applying a Gaussian low pass filter with the Gaussian factor :math:`a` =1. Bottom: Frequency spectrum of synthetic receiver function (solid black) and Gaussian filter with :math:`a` =1 (dashed black). Transparently colored areas correspond to the spectra of large sample draws of synthetic Gaussian correlated noise using different values of :math:`r_{RF}`. The solid colored lines represent the Gaussian curves matching the data 'envelope'. The legend displays corresponding :math:`r_{RF}` and :math:`a` values. For the given receiver function, a proper estimate of :math:`r_{RF}` is 0.98.

If
:math:`r_{RF}` is too large (i.e., very close to 1), :math:`R^{-1}`
becomes instable and small eigenvalues need to be suppressed. The user
can define the cutoff for small singular values by defining
:math:`rcond`. Singular values smaller than :math:`rcond` x the largest
singular value (both in modulus) are set to zero. :math:`rcond` is not
ascribed to the prior dictionary, but to the inversion parameter
dictionary (see configuration file).

The inversion parameters can be subdivided into three categories: (1)
actual inversion parameters, (2) model constraints and (3) saving
options. Parameters to constrain the inversion are the number of chains,
the number of iterations for the burn-in and the main phase, the initial
proposal distribution widths, and the acceptance rate. A large number of
chains is preferable and assures good coverage of the solution space
sampling, as each chain starts with a randomly drawn model only bound by
the priors. The number of iterations should also be set high, as it can
benefit, but not guarantee, the convergence of the chain towards the
global likelihood maximum. The total amount of iterations is
:math:`iter_{total} = iter_{burnin} + iter_{main}`. We recommend to
increase the ratio towards the iterations in the burn-in phase
(i.e., :math:`iter_{burnin}>iter_{main}`), so a chain is more likely to
have converged when entering the exploration phase for the posterior
distribution.

The initial proposal distributions, i.e., Gaussian distributions
centered at zero, for model modifications, must be given as standard
deviations according to each of the model modification methods (:ref:`sec:propmod`). The values must be given as a vector of size five, the order representing following modifications: (1) :math:`\mathrm{V_S}`, (2) depth, (3) birth/death, (4) noise, and (5) :math:`\mathrm{V_P}`/:math:`\mathrm{V_S}`. The first
three distributions represent :math:`\mathrm{V_S}`-depth model
modifications referring to alterations of :math:`\mathrm{V_S}`\ (1,3)
and z (2) of a Voronoi nucleus. There is one proposal distribution for
both noise parameters :math:`r` and :math:`\sigma` (4) and one for
:math:`\mathrm{V_P}`/:math:`\mathrm{V_S}`\ (5).

If the proposal distributions were constant, the percentage of accepted
proposal models would decrease with ongoing inversion progress, i.e.,
the acceptance rate decreases at the expense of an efficient sampling.
To efficiently sample the parameter space, an acceptance rate of
:math:`\sim`\ 40–45 % is forced for each proposal method by dynamically
adapting the width of each proposal distribution. We implemented a
minimum standard deviation of 0.001 for each proposal distribution.

The most accepted model modifications are (1) and (2); their acceptance
rates get easily forced to the desired percentage without even coming
close to the defined minimum width of a proposal distribution. Birth and
death steps, however, barely get accepted after an initial phase of high
acceptance; if not limiting the proposal distribution width to a
minimum, the standard deviations for (3) will get as small as
:math:`10^{-10}` km/s and smaller to try to keep the acceptance rate up.
However, as discussed in :ref:`sec:propmod`,
the distribution width does not in the first place influence the
model-modification, but the added or removed Voronoi nucleus. Models
modified by birth and death steps will naturally not be accepted very
often and even less the further the inversion progresses. Therefore, the
overall acceptance rate is stuck with a specific level below the forced
rate. An estimate of the actual overall acceptance rate can be made,
assuming a realistic acceptance for the birth and death steps, e.g.,
1 %. A user given target rate of 40 % for each method would give an
actual overall acceptance rate of :math:`\sim`\ 3 %.
(:math:`\rightarrow` 6 methods, 4 reach 40 %, 2 only 1 % = 30 % over
all.) The target acceptance rate must be given as an interval.

There are three additional conditions, which might be worthwhile to use
to constrain the velocity-depth model. However, using any of them could
bias the posterior distribution. The user is allowed to set a minimum
thickness of layers. Furthermore low and high velocity zones can be
suppressed. If not ``None``, the value for :math:`lvz` (or :math:`hvz`)
indicates the percentage of allowed :math:`\mathrm{V_S}` decrease (or
increase) from each layer of a model relative to the underlying layer.
For instance, if :math:`lvz`\ =0.1, then a drop of :math:`\mathrm{V_S}`
by 10 %, but not more, to the underlying layer is allowed. As
:math:`\mathrm{V_S}` naturally increases with depth, and the algorithm
only compares each layer with the layer underneath, the :math:`hvz`
criteria should only be used if observing extreme high velocity zones in
the output. Otherwise sharp (but real) discontinuities could be smoothed
out, if chosen too small. The :math:`lvz` and :math:`hvz` criteria will
be checked every time a velocity-depth model is proposed and the model
will be rejected if the constraints are not fulfilled.

The saving parameters include the :math:`station`, :math:`savepath` and
:math:`maxmodels`. The :math:`station` name is optional and is only used
as a reference for the user, for the automatically saved configuration
file after initiation of an inversion. :math:`savepath` represents the
path where all the result files will be stored. A subfolder *data* will
contain the configuration file and all the *SingleChain* output files,
the combined posterior distribution files and an outlier information
file. :math:`savepath` also serves as figure directory.
:math:`maxmodels` is the number of p2-models that will be stored from
each chain.

Running an inversion
~~~~~~~~~~~~~~~~~~~~~

The inversion will start through the *optimizer.mp_inversion* command
with the option to chose the number of threads, :math:`nthreads`, for
parallel computing. By default, :math:`nthreads` is equal to the number
of CPUs of the user’s PC. One thread is occupied if using BayWatch.
Ideally, one chain is working on one thread. If fully exhausting the
capacity of a 8 CPUs PC, give :math:`nthreads`\ =8 and
:math:`nchains`\ =multiple of :math:`nthreads` or (:math:`nthreads`-1)
if using BayWatch. This would cause :math:`nthreads`\ (-1) chains to run
parallel at all times, until :math:`nchains` are worked off.

The speed of the inversion will not increase by choosing a larger
:math:`nthreads`. In fact, the speed is determined by the number of
CPUs. If, for instance, the user doubles :math:`nthreads`, the number of
chains running parallel at once is also double, but the chains queue for
some non-threadable computations blocking one CPU at a time, so each
chain runs half the speed. To decrease :math:`nthreads` offers a
possibility to minimize the workload for a PC and that it is still
accessible for other tasks during an inversion.

Although having access to a cluster, inversions were also performed on a
single work station to determine the duration of an inversion with
standard PC equipment (e.g., Memory: 16 GB, Processor model: 3.60 GHz x
8 cores). The runtime is not only dependent on the PC model, but also on
the number of chains and iterations, and the number of layers of the
actual velocity-depth structures, which directly influences the
computational time of the forward modeling. The inversion for the
example given in :ref:`sec:testdata` with 21 chains,
150,000 iterations and models with 3–10 layers, took 20.4 minutes; so
each batch of 7 chains took 7 minutes.

Another argument to set when starting an inversion is :math:`baywatch`.
If set to True, model data will be send out with an interval of
:math:`dtsend`\ =0.5 s and can be received by BayWatch until the
inversion has finished.

.. _sec:testdata:

Testing with synthetic data
---------------------------

A set of test data was computed with the *BayHunter.SynthObs* module,
which provides methods for computing receiver functions (P, S), surface
wave dispersion curves (Love, Rayleigh, phase, group), and synthetic
noise following the exponential or the Gaussian correlation law. We
computed the P-RF and the fundamental mode SWD of the Rayleigh wave
phase velocity from a six-layer model including a low velocity zone. We
computed non-correlated noise for SWD and Gaussian correlated noise for
the RF with values for :math:`r` and :math:`\sigma` as given in :numref:`Table {number} <tab:testpars>` (*true*). Noise and synthetic data were then
added to create observed data. An example script, including these steps,
can be found in the :doc:`Appendix <appendix>` and the online repository.

.. container::
   :name: tab:testpars

   .. table:: Model priors and inversion parameters for synthetic test inversion and *true* values used for modeling of the observed data. Model prior tuples define the limits (min, max) of a uniform distribution.
   
      +-----------------+------------------+---------------+---------------+
      | vs              | = (2, 5)         | nchains       | = 21          |
      +-----------------+------------------+---------------+---------------+
      | z               | = (0, 60)        | :math:`i      | = 100,000     |
      |                 |                  | ter_{burnin}` |               |
      +-----------------+------------------+---------------+---------------+
      | layers          | = (1, 20)        | :math:`i      | = 50,000      |
      |                 |                  | ter_{main}`   |               |
      +-----------------+------------------+---------------+---------------+
      | vpvs            | = (1.5, 2.1)     | acceptance    | = (50, 55)    |
      +-----------------+------------------+---------------+---------------+
      | :math:`r_{RF}`  | = 0.92           | propdist      | = (0.005,     |
      |                 |                  |               | 0.005, 0.005, |
      |                 |                  |               | 0.005, 0.005) |
      +-----------------+------------------+---------------+---------------+
      | :math:`\sigma_  | = (1e-5,         | rcond         | = 1e-6        |
      | {RF}`           | 0.05)            |               |               |
      +-----------------+------------------+---------------+---------------+
      | :math:`r_{SWD}` | = 0.             | station       | = ’st6’       |
      |                 |                  |               |               |
      +-----------------+------------------+---------------+---------------+
      | :math:`\sigma_  | = (1e-5, 0.1)    |               |               |
      | {SWD}`          |                  |               |               |
      +-----------------+------------------+---------------+---------------+

Two targets (*PReceiverFunction*, *RayleighDispersionPhase*) were
initialized with the "observed" data and combined to a
*BayHunter.JointTarget* object. The latter and the two parameter
dictionaries of model priors and inversion parameters
(:numref:`Table {number} <tab:testpars>`) were passed to the Optimizer. Parameters
that were not defined fall back to the default values. We purposely show
a run with only 150,000 iterations to visualize the convergence of
different chains and the outlier detection method. The inversion
finished after 20 minutes, saving and plotting methods were applied
afterwards.

.. _fig:bh_iiter:

.. figure :: _static/c_iiter_likes_example.png

    Development of likelihood over iteration for all 21 chains (top) and a small selection of chains (bottom).
  
:numref:`Figure {number} <fig:bh_iiter>` shows the likelihood
development over the iterations for all and for a selection of chains. A
strong increase of likelihood can be observed at the first iterations in
the burn-in phase, converging towards a stable value with increasing
number of iteration. Some chains reached the final likelihood plateau in
the burn-in phase (e.g., :math:`c0`), some within the posterior sampling
phase (e.g., :math:`c4`), and some did not converge at all (:math:`c2`).
The chain :math:`c2` (also :math:`c1` and :math:`c3`) had a good chance
of reaching the maximum likelihood, if the small number of iterations
would not have stopped the exploration at this early stage. However, the
number of iterations cannot be eternal; in any case it is necessary to
compare the convergence level of the chains.

Here, we defined a 0.02 deviation condition for outliers. With a maximum
median posterior likelihood of 1674 (:math:`c13`), the likelihood
threshold is 1640, which declared 13 chains with deviations of
0.032–0.159 as outliers (see :numref:`Table {number} <tab:outliers>`). In a real
case inversion, the number of iterations should be much higher, and the
number of outlier chains is small. The detected outlier chains will be
excluded from the posterior distribution.

.. container::
   :name: tab:outliers

   .. table:: Deviations of each chain’s median likelihood from the maximum median likelihood of the chain ensemble. Only outlier chains with deviations :math:`>`\ 0.02 (2 %) are listed.

      == ===== == ===== === ===== === ===== === =====
      c1 0.039 c6 0.061 c9  0.059 c15 0.150 c19 0.059
      c2 0.111 c7 0.033 c10 0.150 c16 0.033     
      c5 0.032 c8 0.109 c14 0.033 c17 0.033     
      == ===== == ===== === ===== === ===== === =====

:numref:`Figure {number} <fig:bh_datafits>` shows the current :math:`\mathrm{V_S}`-depth models from different chains and corresponding data fits (same chains as in :numref:`Figure {number} <fig:bh_iiter>`, bottom). Chains :math:`c1` and :math:`c2` show the worst data fits; they were declared as outliers. The
other chains (:math:`c0`, :math:`c3`, :math:`c4`) show a reasonably good
data fit with very similar velocity models. Chains :math:`c0` and
:math:`c4` already found a six-layer model, :math:`c3` found a
five-layer model averaging the low velocity zone.

.. _fig:bh_datafits:

.. figure :: _static/c_currentmodels.png

    Current velocity-depth models and data fits of corresponding SWD and RF data from different chains with likelihoods as illustrated in :numref:`Figure {number} <fig:bh_iiter>` (bottom). The black line (left) is the synthetic :math:`\mathrm{V_S}`-depth structure.


The posterior distribution of the eight converged chains, containing 100,000 models, are illustrated in :numref:`Figure {number} <fig:bh_models>`. The mean (and mode) posterior :math:`\mathrm{V_S}`-depth structure images the true model
very well, including the low velocity zone. The number of layers is determined to be most likely six. The :math:`\sigma` distributions of both, RF and SWD show a Gaussian shape, inhering a tail of higher values from models of chains that only converged within the exploration phase of the inversion (e.g., :math:`c3` and :math:`c4`). The distribution of :math:`\sigma_{SWD}` already represents a good estimate, slightly overestimated, which falls back to the number of iterations. Tests with more iterations show that the median of :math:`\sigma_{SWD}` is in
perfect agreement with the true value.

.. _fig:bh_models:

.. figure :: _static/c_posteriormodels.png

    Recovered posterior distributions of :math:`\mathrm{V_S}`, interface depth, number of layers, and noise level for synthetic data. Red lines indicate the true model, as given in :numref:`Table {number} <tab:testpars>`. The posterior distribution is assembled by 100,000 models collected by 8 chains.

:math:`\sigma_{RF}` is underestimated, which theoretically means that
noise was interpreted as signal and receiver function data is
overfitted. The difference to SWD is the type of noise correlation (=
Gaussian) and the assumption of the correlation :math:`r` of data noise
(:math:`r\neq0`). We computed synthetic RF data applying a Gaussian
lowpass filter with a Gaussian factor of 1. Separately, noise was
generated randomly with a correlation :math:`r` estimated to represent
the applied Gauss filter, and added to the synthetic data. The random
process of generating noise does not output a noise vector which exactly
matches the given values of :math:`r` and :math:`\sigma`. If only
drawing one single realization with a determined amount of samples from
the multivariate normal distribution may produce deviations from the
targeted :math:`r` and :math:`\sigma`. From the generated noise the true
:math:`\sigma` can be computed by the standard deviation. However, the
true :math:`r` is not easy to reconstruct. Assuming a wrong :math:`r`
for the covariance matrix of noise cannot lead to the correct
:math:`\sigma`.

.. _fig:bh_corr:

.. figure :: _static/c_covfix.png

    Comparison of residuals of the best fitting RF model and one realization of noise through :math:`C_e^{RF}` for receiver functions. Both noise vectors are of coherent appearance in frequency and amplitude, hence, the estimate of :math:`r_{RF}` is appropriate.

It is possible to clarify whether the assumed correlation parameter
:math:`r` is in tendency correct. :numref:`Figure {number} <fig:bh_corr>` shows a comparison of (1) the RF data residuals of the best fitting model and (2) one realization of noise with the given correlation :math:`r` and the estimated
:math:`\sigma_{RF}`; both noise vectors should be of coherent appearance
in frequency and amplitude, if the estimate of :math:`r_{RF}` is
appropriate.
