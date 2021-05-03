# *If you like it, you should `put-a-ring-on-it`*:

*Modelling photon orbits and gravitational lensing in a general relativistic framework*

## Getting started
### Pre-requisites

Keeping track of the names and versions of the python modules used.
```
argparse, scipy, numpy, matplotlib, PIL, gaussxw
```

### Installing

### Use

`initialize_photons.py` takes in an image / object and converts it into a numpy array for processing

`rk4.py` contains the rk4 solver module we can call

`integral.py` contains the integral solver module we can call

`helmholtz.py` can be called from the command line to produce images of a Gaussian object transiting behind a lens of chosen mass, using Helmholtz regression
`python helmholtz.py -s rk4 -l 1e3 --savefile` uses rk4 to output and save (`--savefile` is an optional argument) images, created with lens mass `-l` = 1e3 solar masses.

## Licensing

Georgia Institute of Technology, 2021
Computational Physics 6260

## Authors

* Peter Addison
* Mi Do
* Si Ferrel
* Tamir Gonen Cohen
* Julie Malewicz
