# BayHunter

BayHunter is an open source Python tool to perform an McMC transdimensional Bayesian inversion of receiver functions and/ or surface wave dispersion. It is inverting for the velocity-depth structure, the number of layers, Vp/Vs ratio and noise parameters (noise correlation and amplitude). The forward modeling codes are provided within the package, but are easily replaceable with own codes. It is also possible to add (completely different) data sets.

The BayWatch module can be used to live-stream the inversion while it is running: this makes it easy to see how each chain is exploring the parameter space, how the data fits and models change and in which direction the inversion progresses.

**Citation:**

> Dreiling, Jennifer; Tilmann, Frederik (2019): BayHunter - McMC transdimensional Bayesian inversion of receiver functions and surface wave dispersion. V. 1.0. GFZ Data Services.

[http://doi.org/10.5880/GFZ.2.4.2019.001](http://doi.org/10.5880/GFZ.2.4.2019.001)

## Quick start

### Requirements
* matplotlib
* numpy
* pyPdf
* configobj
* zmq
* rfmini, only if inverting for RF (`rfmini.tar.gz`)

### Installation (python2 environment)

```sh
git clone https://github.com/jenndrei/BayHunter.git
cd BayHunter
sudo python setup.py install
```

### Tutorial

An example of how to run an inversion can be found in the **tutorial folder**.
The file to be run `tutorialhunt.py` is spiked with comments.
You can also create your own synthetic data set with `create_testdata.py`.

Use the input file `config.ini` for adjusting the inversion parameters.

More background information about how to chose the best parameters, and about BayHunter and BayWatch in general can be found in the file `docs/bayhunter.pdf`.

### References

* SWD forward modeling is based on surf96 from [CPS](http://www.eas.slu.edu/eqc/eqccps.html) from Rob Herrmann, St. Louis University: BayHunter uses the python wrapper [pysurf96](https://github.com/miili/pysurf96) from Marius Isken wrapping the
quick surf96 routine [SurfTomo](https://github.com/caiweicaiwei/SurfTomo) from Hongjian Fang.
* RF forward modeling using [rfmini](https://git.gfz-potsdam.de/saul/rfmini) from Joachim Saul, GFZ.
* Most influence offered the work from Bodin et al., 2012: *Transdimensional inversion of receiver functions and surface wave dispersion*.

## Outlook and Feedback

**BayHunter is ready to use**. It is quick and efficient and I am happy with the performance. Still, there are always things that can be improved to make it even faster and more efficient, and user friendlier.  

BayHunter was mostly tested with a joint data set of RF and SWD and depths down to 80 km. Colleagues tested BayHunter using:  
1. only one SWD with depths down to 200 km (real data)  
2. joint SWD down to 30 km including very low surface velocities (real data)  
3. RF and an additional user data set (synthetic data).

Thus, we could eliminate some problems. However, each data set and each inversion has its own characteristics. If you observe any unforeseen behavior, please share it with me to wipe out possible problems we haven't considered.

I am happy to share my experience with you and also if you share your thoughts with me. I am looking forward to your feedback. 

## Who am I?

I am Jennifer Dreiling, final sprint PhD candidate at GFZ (German Research Center for Geosciences) in Potsdam, Germany. BayHunter was created by me in the frame of my PhD program. [Contact me](https://www.gfz-potsdam.de/en/staff/jennifer-dreiling/).
