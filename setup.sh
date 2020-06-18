rm -r objdir
#cmake -H. -Bobjdir -DCMAKE_Fortran_COMPILER=gfortran-7
#cmake -H. -Bobjdir -DINTEL=ON
cmake -H. -Bobjdir  $@
#cmake -H. -Bobjdir
