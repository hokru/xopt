
# XOPT - an eXternal OPTimizer
![CI](https://github.com/hokru/xopt/workflows/CI/badge.svg)

## Purpose
The goal is to proving a robust optimizer for quantum chemical and semi-empirical method
that is suitable for large and complex molecules.

## Notes/Versions
The current version (2.0 beta) is under continuous development and no warranty for correctness can be given. It is a significant extension and re-write of the legacy version published in H. Kruse, J. Sponer PCCP, 2015,17, 1399-1410 that introduced the approach of restrained optimizations for biomolecules.

## build
Standard way of building is using cmake:
```
cmake -H. -Bobjdir <flags>
cmake --build objdir
```

Available compiler flags are:
* `-DBLAS=MKL/OpenBLAS/Generic`
 Compiler can be set via `$FC` variable `-DCmake_Fortran_COMPILER=` flag or using one of the following flags
* -DGNU=ON (gfortran)
* -DINTEL=ON (ifort)
* -DPGI=ON (pgfortran)

One can help the BLAS/LAPACK autodetection setting the MATHROOT variable in the shell.

Alternatively, building via Makefile is still possible (see configs/Makefile.xxx for examples).

## Manual
Execute `xopt -h` for command line options.
Online documentation (unfinished): [![Documentation Status](https://readthedocs.org/projects/xopt/badge/?version=latest)](http://xopt.readthedocs.io/en/latest/?badge=latest)
Check `Manual.pdf` for complimentary options.


## customatization
getgrad.f90 contains most of the system calls which might need adaption to your work environment.
Most system calls can also be set in `$HOME/.xoptrc`. 

See also the online documentation.
