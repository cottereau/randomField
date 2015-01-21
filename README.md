# Random Fields Generation Library


## PRESENTATION

This library aims to generate random fields with prescribed first-order marginal distribution and correlation structure over structured grids or non-structured point arrays. Currently the library offers:

1. Random field generation using the spectral method by Shinozuka and Deodatis and a variation suited for isotropic media.

This software is mainly developed at laboratoire MSSMat (Ecole Centrale Paris - CNRS).

* contact : [Luciano de Carvalho Paludo](mailto:luciano.de-carvalho@ecp.fr)
* contributors (by order of first commit): R. Cottereau, L. Paludo, V. Bouvier

It is developed in FORTRAN 90.

## REFERENCES

Original references for the theory are:

1. Shinozuka, M. and Deodatis, G., Simulation of multi-dimensional Gaussian stochastic fields by spectral representation, App. Mech. Rev. 49 (1), 1996, pp. 29-53.
1. Shinozuka, M. and Deodatis, G., Simulation of stochastic processes by spectral representation, App. Mech. Rev. 44 (4), 1991, pp. 191-205.


## INSTALLATION

### Requirements

There are two librairies needed to build the Random Fields Generation Library: 

HDF5 : www.hdfgroup.org is a library used to manage big size files
MPI  : a MPI library

The path to these libraries should be informed in the Makefile (`LIBHDF5`, `INCLUDEHDF5`, `LIBMPI` and `INCLUDEMPI`)

### Compiling

Inside the folder randomField the command `make` will build the file "librandomfield.a" that will be file of the library that your makefile should point to.


## USE

To create "N" realizations the inputs are:

* Points coordinates (xPoints)
* Correlation model
* First-order marginal law
* Correlation length
* Field average
* Field variance
* Number of realizations (N)
* Generation method

## SYNTAX

The syntax to call the function is:

`createStandardGaussianFieldUnstruct(xPoints, corrMod, margiFirst, corrL, fieldAvg, fieldVar, Nmc, method, randField)`


Where:

xPoints - is a matrix where each column contains X, Y and Z coordinates (real numbers) for each point.

          |X1  X2  X3  ...... Xn|
          |Y1  Y2  Y3  ...... Yn|
          |Z1  Z2  Z3  ...... Zn|


corrMod - string containing the name of the correlation model. Only the "gaussian" option is implemented.

margiFirst - string containing the name of the first-order marginal law. The options are "gaussian" or "lognormal"

corrL - vector of positive real numbers containing the correlation length in each direction

          |Correlation Length in X|
          |Correlation Length in Y|
          |Correlation Length in Z|

fieldAvg - real number containing the field average

fieldVar - positive real number containing the field variance

Nmc - integer representing the number of realisations

method - generation method (integer). The options are "2" for Shinozuka and Deodatis method and "1" for the isotropic spectrum optimized version

randField - matrix tha will store the results. Each column has a realization. If each realization is represented by a letter and each point by a number it would be a matrix like:

          |A1 B1 C1 ..... |
          |A2 B2 C2 ..... |
          |A3 B3 C3 ..... |

Obs: all the variables should be allocated before the call

