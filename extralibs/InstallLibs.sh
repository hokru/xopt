#!/bin/bash

echo "Making library files:"

echo "--> d3lib"
cd d3lib
 make clean
 make &>> ../make.out
 tail -1 ../make.out
cd ..

echo "--> gcplib"
cd gcplib
 make clean
 make &>> ../make.out
 tail -1 ../make.out
cd ..

echo "--> lbfgs"
cd lbfgs
 make clean
 make &>> ../make.out
 tail -1 ../make.out
cd ..

echo " "
echo " "
echo "* Library files *"
ls *.a
echo " "

