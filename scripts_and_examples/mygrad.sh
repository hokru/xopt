#!/bin/bash
#
# mygrad.sh example Turbomole
#
# run define to prepare usual Turbomole input first

#1. make use of xopt.xyz
babel -ixyz xopt.xyz -otmol coord 

#2. execute
dscf > dscf.out

#3. write energy into file
grep "|  total energy" dscf.out | awk '{print $5}' > xopt.energy.tmp 
