.. role:: raw-latex(raw)
   :format: latex


Overview
========

Abstract
--------

BayHunter is an open source Python tool to perform an McMC
transdimensional Bayesian inversion of surface wave dispersion and/or
receiver functions. The algorithm follows a data-driven strategy and
solves for the velocity-depth structure, the number of layers,
:math:`\mathrm{V_P}`/:math:`\mathrm{V_S}` ratio and noise parameters,
i.e., data noise correlation and amplitude. The forward modeling codes
are provided within the package, but are easily replaceable with own
codes. It is also possible to add (completely different) data sets.

The BayWatch module can be used to live-stream the inversion while it is
running (see animation below). This makes it easy to see how each chain is exploring the
parameter space, how the data fits and models change and in which
direction the inversion progresses.

.. video:: _static/baywatch.mp4
   :width: 85%


Available on GitHub: https://github.com/jenndrei/BayHunter


Citation
--------

Dreiling, J., Tilmann, F. BayHunter – McMC transdimensional Bayesian inversion of receiver functions and surface wave dispersion. GFZ Data Services, 2019. https://doi.org/10.5880/GFZ.2.4.2019.001.


Application examples
--------------------

-  Dreiling et al., 2020: Crustal structure of Sri
   Lanka derived from joint inversion of surface wave dispersion and
   receiver functions using a Bayesian approach. Journal of Geophysical
   Research: Solid Earth. https://doi.org/10.1029/2019JB018688.

-  Green et al., 2020: Magmatic and sedimentary
   structure beneath the Klyuchevskoy Volcanic Group, Kamchatka, from
   ambient noise tomography. Journal of Geophysical Research: Solid
   Earth. https://doi.org/10.1029/2019JB018900.

-  Mauerberger et al., n.a.: The multifaceted
   Scandinavian lithosphere imaged by surface waves and ambient noise.
   In preparation.


Documentation
-------------

The documentation to BayHunter contains three chapters, accessible through the navigation menu, and downloadable as `PDF <https://github.com/jenndrei/BayHunter/blob/master/documentation/BayHunter_v2.1_documentation.pdf>`_. 

- **Introduction**: Bayes theorem, McMC inversion approach

- **Algorithm**: Python framework behind BayHunter, including optimizer and chain modules, saving and plotting options, and Baywatch

- **Tutorial**: Setting up and running an inversion, targets and parameters, example inversion using synthetic data (minimalistic working example of the code in the appendix)


Comments and Feedback
---------------------

BayHunter is ready to use. It is quick and efficient and I am happy with
the performance. Still, there are always things that can be improved to
make it even faster and more efficient, and user friendlier.

Although we tested the software with a variety of synthetic and real
data, each data set is still unique and shows own characteristics. If
you observe any unforeseen behavior, please share it with me to wipe out
possible problems we haven’t considered.

I am happy to share my experience with you and also if you share your
thoughts with me. I am looking forward to your feedback.


-------------------------------------------------

Quickstart
==========

Requirements
------------

- ``matplotlib``

-  ``numpy``

-  ``pyPdf``

-  ``configobj``

-  ``zmq``

-  ``Cython``

Installation
------------

(compatible with Python 2 and 3)

::

   git clone https://github.com/jenndrei/BayHunter.git
   cd BayHunter
   sudo python setup.py install


Resources
----------

-  Algorithm: based on the work of :doc:`Bodin et al., 2012 <references>`.

-  SWD forward modeling: ``SURF96`` from Computer Programs in Seismology (:doc:`Herrmann and Ammon, 2002 <references>`). Python wrapper using `pysurf96 <https://github.com/miili/pysurf96>`_ and `SurfTomo <https://github.com/caiweicaiwei/SurfTomo>`_.

-  RF forward modeling: ``rfmini`` (`Joachim Saul, GFZ <https://www.gfz-potsdam.de/en/staff/joachim-saul/>`_).

-------------------------------------------------

.. toctree::
    :maxdepth: 2
    :hidden:
    
    introduction
    algorithm
    tutorial
    references
    appendix
    FAQs
