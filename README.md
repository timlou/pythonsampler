# pythonsampler
## Basics
sampler.sh -- bash script to run the code.

Usage: ./sampler.sh inputfile.txt outputfile.txt n_pts_generated

PointsGenerator.py -- main python code that contains the PointsGenerator() class for parsing files and generating points.
sobol_seq.py -- sobol sequence generator taken from https://github.com/naught101/sobol_seq.git

## Prerequisites
python3.5 is used, other versions may or may not work
sympy package is utilized for symbol simplications for inequalities
vegas package is used for adaptive monte carlo integration

## Algorithm
The code first parses the input file for inequalities, and attempt to simplify it using sympy. So far only the most basic univariate inequalities can be simplified and converted into basic boundaries for individual variables. Then, there are two algorithms that attemp to generate sample points:

1. Basic Algorithm: quasi-random numbers are generated on the variables boundaries via the sobol_seq.py. This offers good uniform coverage. Invalid points are simply discarded
2. Adaptive Algorithm: an effective potential is used as an integrand (which vanishes when the points are invalid). Using the vegas package, an integration is performed, and the sampled points extracted.


