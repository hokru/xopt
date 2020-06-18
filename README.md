
# XOPT - an eXternal OPTimizer

## Purpose
The goal is to proving a robust optimizer for quantum chemical and semi-empirical method
that is suitable for large and complex molecules.

## Notes/Versions
The current version (2.0 beta) is under continuous development and no warranty for correctness can be given at the current development stage.

The legacy version (v1.0.1) is published in H. Kruse, J. Sponer PCCP, 2015,17, 1399-1410. 

## Manual
Execute `xopt -h` for command line options.
Online documentation (unfinished): [![Documentation Status](https://readthedocs.org/projects/xopt/badge/?version=latest)](http://xopt.readthedocs.io/en/latest/?badge=latest)
Check `Manual.pdf` for complimentary options.


## customatization
getgrad.f90 contains most of the system calls which might need adaption to your work environment.
Most system calls can also be set in `$HOME/.xoptrc`. 

See also the online documentation.

## build
You can build xopt with a standard Makefile (see configs/Makefile.xxx for examples).
Also cmake is supported:
cmake -H. -Bobjdir <flags>
cmake --build objdir

Available compiler flags are:
* -DGNU=ON (gfortran)
* -DINTEL=ON (ifort)
* -DPGI=ON (pgfortran)
* -DBLAS=MKL/OpenBLAS/Generic

One can help the BLAS/LAPACK autodetection setting the MATHROOT variable in the shell.

