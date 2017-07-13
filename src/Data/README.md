
# Synopsis
This folder contains data files for the analysis of various water models.
Below is a brief description of the data format.

## Equation of state
Here the data is of the form:
1. density = f(temperature, pressure, model)
2. pressure = f(temperature, density, model)

The files for the equation state are of the form with regex:
1. compiledEOSAll*
2. pruned*compiledEOSAll


### Pruning
The data is also pruned (see ../src/prune.py) on the basis of structural
relaxation F_s(k,t), i.e. "pruned"+xxx+"compiledEOSAll", where xxx is the upper 
bound of the relaxation time in nanoseconds.