# XOPT - an eXternal OPTimizer

## Purpose
The goal is to proving a robust optimizer for quantum chemical and semi-empirical method
that is suitable for large and complex molecules.

## Capabilities

Coordinate systems:
- approx. normal coordinates
- cartesian "

Hessian Updates:
- modified SG1-BFGS
- ...

Step determination:
- RFO
- SI-RFO
- conj. gradient (cg)
- RFO-cg mixture

Interfaces:
- Turbomole
- ORCA
- Amber
- mopac

Constraints/Restraints:
- cartesian space constraints
- restrained primititves (form: U(x)=k(x)^2, x=bond/angle/torsion deviation)


## USAGE

see manual.pdf

Adapt getgrad.f90 according to your own environment!


