# BayHunter v2.1

BayHunter is an open source Python tool to perform an McMC transdimensional Bayesian inversion of surface wave dispersion and/or receiver functions. The algorithm follows a data-driven strategy and solves for the velocity-depth structure, the number of layers, Vp/Vs ratio and noise parameters, i.e., data noise correlation and amplitude. The forward modeling codes are provided within the package, but are easily replaceable with own codes. It is also possible to add (completely different) data sets.

The BayWatch module can be used to live-stream the inversion while it is running: this makes it easy to see how each chain is exploring the parameter space, how the data fits and models change and in which direction the inversion progresses.


### Citation

> Dreiling, Jennifer; Tilmann, Frederik (2019): BayHunter - McMC transdimensional Bayesian inversion of receiver functions and surface wave dispersion. GFZ Data Services. [http://doi.org/10.5880/GFZ.2.4.2019.001](http://doi.org/10.5880/GFZ.2.4.2019.001)


### Application examples

> Dreiling et al. (2020): Crustal structure of Sri Lanka derived from joint inversion of surface wave dispersion and receiver functions using a Bayesian approach. Journal of Geophysical Research: Solid Earth. [https://doi.org/10.1029/2019JB018688](https://doi.org/10.1029/2019JB018688).

> Green et al. (2020): Magmatic and sedimentary structure beneath the Klyuchevskoy Volcanic Group, Kamchatka, from ambient noise tomography. Journal of Geophysical Research: Solid Earth. [https://doi.org/10.1029/2019JB018900](https://doi.org/10.1029/2019JB018900).

> Mauerberger et al. (n.a.): The multifaceted Scandinavian lithosphere imaged by surface waves and ambient noise. In preparation.


### Comments and Feedback

BayHunter is ready to use. It is quick and efficient and I am happy with the performance. Still,
there are always things that can be improved to make it even faster and more efficient, and user
friendlier. 

Although we tested the software with a variety of synthetic and real data, each data set is
still unique and shows own characteristics. If you observe any unforeseen behavior, please share
it with me to wipe out possible problems we haven’t considered.

I am happy to share my experience with you and also if you share your thoughts with me. I am looking forward to your feedback. 


### Who am I?

I am Jennifer Dreiling. I finished my PhD studies at the German Research Center for Geosciences (GFZ) in Potsdam, Germany. I created BayHunter in the frame of my PhD program. [Contact me](https://www.gfz-potsdam.de/en/staff/jennifer-dreiling/).



## Quick start

### Requirements
* matplotlib
* numpy
* PyPDF2
* configobj
* zmq
* Cython

### Installation (compatible with Python 2 and 3)*

*Although BayHunter is currently compatible with Python 2 and 3, we recommend you to upgrade to Python 3, as the official support for Python 2 has stopped in January 2020.

```sh
git clone https://github.com/jenndrei/BayHunter.git
cd BayHunter
sudo python setup.py install
```

### Documentation and Tutorial

The documentation to BayHunter offers background information on the inversion algorithm, the parameters and usage of BayHunter and BayWatch (tutorial). See the [documentation here](https://jenndrei.github.io/BayHunter/) or download the [PDF](https://github.com/jenndrei/BayHunter/blob/master/documentation/BayHunter_v2.1_documentation.pdf). Also check out the [FAQs](https://jenndrei.github.io/BayHunter/FAQs).

An example inversion can be found in the **tutorial folder**.
The file to be run, `tutorialhunt.py`, is spiked with comments.
You can also create your own synthetic data set with `create_testdata.py`.


### Resources

* Algorithm: based on the work of [Bodin et al., 2012](https://doi.org/10.1029/2011JB008560).
* SWD forward modeling: SURF96 from Computer Programs in Seismology ([CPS](http://www.eas.slu.edu/eqc/eqccps.html)). Python wrapper using [pysurf96](https://github.com/miili/pysurf96) and [SurfTomo](https://github.com/caiweicaiwei/SurfTomo).
* RF forward modeling: **rfmini** from [Joachim Saul, GFZ](https://www.gfz-potsdam.de/en/staff/joachim-saul/).
