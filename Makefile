
  PROG = xopt


  OBJS1=modules.o main.o io.o  string.o eval_opt.o dist.o print.o getHess.o copt.o getgrad.o hessup.o progs.o control.o \
       anc.o ancopt.o lambda.o defaults.o constrains.o conjgrad.o TRproj.o intbond.o
  LIB = lib/timestamp.o
  OBJS = $(OBJS1) $(LIB)

#--------------------------------------------------------------------------
# Testing for 1. gfortran 2. ifort
FC := $(shell which gfortran 2>> error.dat)
FC := $(shell which ifort 2> error.dat) 
ifndef FC
#ifeq ($(FC),)
$(warning ifort not found.)
$(warning gfortran not found.)
$(warning Please adjust FC= yourself!)
$(error aborting...)
endif

#--------------------------------------------------------------------------
# BLAS library needed
#  FC = ifort 
   FC = ifort   -mkl=parallel 
#  FC = ifort  -mkl=sequential 
  LINKER = $(FC) -static 
  FFLAGS =    -O 

# MKL  (static)
#  FC = ifort  -mkl=parallel 
#  FC = ifort  -mkl=sequential 
#  LINKER = $(FC) -static 
#FFLAGS =    -O
 
        
#dynamic openblas
#OPENBLAS = /usr/qc/openblas_generic/
#FC = ifort
#LINKER = $(FC) -L$(OPENBLAS)/lib/ -lopenblas
#FFLAGS =    -O


.PHONY: all
.PHONY: clean


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@


BUILID:=$(shell date)
#
all:$(PROG)
	$(shell echo "print*,' build date: XXX'" > version.dat)
	$(shell sed -i s/"XXX"/"$(BUILID)"/ version.dat)
	@echo '*** done ***'

$(PROG):$(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

tarball:
	tar -czf xopt.tar.gz *.f90 Makefile lib/timestamp.f90 version.dat MANUAL.pdf
	zip xopt.zip *.f90 Makefile lib/timestamp.f90 version.dat MANUAL.pdf

clean:
	rm -f *.o *.mod $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))
