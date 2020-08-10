.. role:: raw-latex(raw)
   :format: latex

.. _sec:faq:

FAQs
====


- **Which Python version requieres BayHunter?**

BayHunter is compatible with Python2 and Python3.


- **What noise sources are represented in the noise amplitude?**

Multiple sources of noise are included in :math:`\sigma` and consider observational errors, processing errors and theory errors (:doc:`Bodin et al., 2012 <references>`). Below is an overview of these sources.

- Observational errors from background seismic noise (micro-seisms) and instrumental noise.

- Processing errors introduced by parameter choices and generally unstable operations (e.g., the deconvolution for receiver functions).

- Theory errors from forward modeling and inversion algorithms through physical approximations of the Earth, e.g., horizontal homogeneous isotropic layers. This noise type is coherent and reproducible, and part of the signal chosen not to explain (e.g., Gaussian filtering for receiver functions with frequency band corresponding to major scale structures).

These contributions are not simply additive and are all partly reflected in the noise amplitudes :math:`\sigma_{RF}` and :math:`\sigma_{SWD}`


- **How to choose priors for the noise scaling parameters?**

*Noise amplitude* :math:`\sigma` :

The noise parameters :math:`\sigma_{RF}` and :math:`\sigma_{SWD}` represent the level or the amplitude of noise. The ranges should be given sufficiently wide, so for your first tests you can try something like (1e-5, 0.05) for :math:`\sigma_{RF}` and (1e-5, 0.1) for :math:`\sigma_{SWD}`. If a prior range was chosen too small, you will see it after the inversion, if your posterior distribution of :math:`\sigma` has settled on a limit value. Be sure to choose the prior minimum to be larger than zero (i.e., not exactly 0).

As an example: The inversions in :doc:`Dreiling et al. (2020) <references>` yielded :math:`\sigma_{SWD}` between 0.0043--0.0241 km/s (median: 0.0078 km/s) and :math:`\sigma_{RF}` between 0.0048–0.0157 (median: 0.0082).

*Noise correlation factor* :math:`r` :

Because RF are time series, adjacent data points show a correlation. If you computed RF using an exponential filter, you can give a range for :math:`r_{RF}`. The standard processing routine to construct RF, however, is to apply a Gauss filter. Therefore, you cannot invert for :math:`r_{RF}`, as this would take a tremendous amount of computation time for the inversion, but instead need to give an estimate (digit) for :math:`r_{RF}`. As commented by :doc:`Bodin et al. (2012) <references>`: Using the exponential correlation law to model RF uncertainties (of Gaussian filtered RF) would erroneously assume high frequency components in the noise which have obviously been cut by the Gaussian filter.

BayHunter provides a tool to estimate :math:`r_{RF}`, as explained in the documentation (:ref:`sec:parsetup`). For a code example, see the :doc:`Appendix <appendix>` or the tutorial folder in the repository.

Since dispersion data are not time series but travel‐time measurements at different periods, you can consider the noise independent (i.e., not correlated). So for SWD you can keep the correlation :math:`r_{SWD}` =0.

