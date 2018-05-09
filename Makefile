# THIS IS THE MAIN MAKEFILE FOR XOPT
# NOTES:
# - BLAS/LAPACK:
#   xopt relies on various lapack and blas functions. However, speed is not crucial unless one
#   runs very large molecules. Thus any library will do the trick. 
#   OpenBLAS is an excellent and free threaded blas/lapack implementation.
#   for Intel MKL use new mkl_rt, mkl link advisor or ifort -mkl=<..>
#   MKL is available in Anaconda (miniconda). mkl and mkl-devel are needed.
#
# GLOBAL SETTINGS
USE_DEV=no
PROJDIR :=./
SOURCEDIR := src
BUILDDIR := build
MOD:= mod
#######################
# USER CONFIGURATION  #
#######################
# CHOOSE TARGET LOCATION
PROG = ~/bin/xopt

# choose compiler
USE_GNU=yes
#USE_PGI=no
#USE_INTEL=yes

# static/dynamic
DO_STATIC=yes

# set non-essential flags (at least -O2 recommened) 
#FFLAGS= -O0 -g -traceback -check bounds  -openmp -fpe0 
FFLAGS= -O2 

## set BLAS/LAPCK paths for optimized libraries
 OPENBLAS = /usr/qc/OpenBLAS/
 BLASLIB = -L$(OPENBLAS) -lopenblas -lpthread

# MKLROOT = ${HOME}/miniconda3/
# BLASLIB =  -L${MKLROOT}/lib/ -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl  # for 'new' single-dynamic MKL
# BLASLIB  = -mkl=parallel                                   # for ifort

###        examples      ###
#  BLASLIB  = -L/local/intel_mkl -<see intel link advisor>    # MKL explicit
#  BLASLIB  = /local/openblas_lib/lib/ -lopenblas             # OpenBLAS (may need -lpthread)
#  BLASLIB  = -L/usr/lib64 -llapack -lblas                    # native blas
#  BLASLIB  = -mkl=parallel                                   # for ifort

### EXTRALIBS (dftd3,gcp,lbfgs)
# set to "yes" after building them.
USE_EXTRA=yes

USE_DEV=no

#############################################################
# DO NOT MODIFY BELOW UNLESS YOU KNOW WHAT YOU ARE DOING !  #
#############################################################
#
ifeq ($(USE_DEV),yes)
USE_GNU=yes
USE_INTEL=no
FFLAGS =  -fbounds-check -Og -g  -Wunused-dummy-argument -std=f2008 -fall-intrinsics -fbacktrace
#OPTIM =  -fopt-info-vec -fopt-info-loop  
endif

EXTRAS=version.dat README.md Makefile header.tmp place.sh

FPPSOURCES:= modules.F90 main.F90 print.F90

FSOURCES:=\
 string.f90\
 Hmass.f90\
 TRproj.f90\
 anc.f90\
 ancopt.f90\
 ciopt.f90\
 conjgrad.f90\
 constrains.f90\
 control.f90\
 copt.f90\
 defaults.f90\
 eval_opt.f90\
 gdiis.f90\
 getHess.f90\
 getgrad.f90\
 hessup.f90\
 intcoords.f90\
 intopt.f90\
 io.f90\
 irc_opt.f90\
 lambda.f90\
 math.f90\
 md.f90\
 molecule.f90\
 oniom.f90\
 printmat.f90\
 progs.f90\
 random.f90\
 tools.f90\
 hbonds.f90

OFPP= $(foreach dir, $(FPPSOURCES), $(addprefix $(SOURCEDIR)/, $(dir)))
OF= $(foreach dir, $(FSOURCES), $(addprefix $(SOURCEDIR)/, $(dir)))
OBJ := $(subst $(SOURCEDIR),$(BUILDDIR),$(OFPP:%.F90=%.o)) $(subst $(SOURCEDIR),$(BUILDDIR),$(OF:%.f90=%.o))


# build and git info
BUILID:=$(shell date)
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
$(info Building: $(GIT_VERSION))


#####  INTEL COMPILER  ##########
ifeq ($(USE_INTEL),yes)
FC = ifort # -check all -fpe0 -warn -traceback -debug extended
DFLAGS+= -DINTEL
FC+= -module $(MOD)
LIBS+= $(BLASLIB)
  ifeq ($(DO_STATIC),yes)
   FC+= -static
  endif
endif
       

##### GFORTRAN #########
ifeq ($(USE_GNU),yes)
FC = gfortran 
DFLAGS+= -DGNU
FFLAGS+= -J$(MOD) -ffree-line-length-none
LIBS+= $(BLASLIB)
  ifeq ($(DO_STATIC),yes)
   FC+= -static-libgfortran -static
  endif
endif


##### PGI FORTRAN #############
# PLEASE ADAPT THIS YOURSELF
ifeq ($(USE_PGI),yes)
FC=pgfortran -Bstatic_pgi -pgf90libs -Mbackslash  -Bstatic
LIBS =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_pgi_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl 
FFLAGS= -tp=x64 -O
DFLAGS+= -DPGI
endif

ifeq ($(USE_EXTRA),yes)
LIBS+=-L$(PROJDIR)/extralibs/ -ldftd3 -lgcp -llbfgs
DFLAGS+= -DDFTD3 -DGCP -DLBFGS
endif

# make targets:
.PHONY: all
.PHONY: clean
.PHONY: version
.PHONY: tarball
.PHONY: pgrad
.PHONY: version.dat
.PHONY: extralibs


$(BUILDDIR)/%.o $(BUILDDIR)/%.mod: $(SOURCEDIR)/%.F90
	$(FC) $(DFLAGS) $(FFLAGS)  -c $< -o $@  

$(BUILDDIR)/%.o: $(SOURCEDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(PROG): version.dat $(OBJ) 
	$(FC) $(OBJ) $(LIBS) -o $(PROG)
	@echo '*** all done ***'


#%.o: %.F90
#	@echo "making $@ from $<"
#$(FC) $(DFLAGS) $(FFLAGS)  -c $< -o $@  


version.dat:
	@echo 'writing new version.dat'
	@echo     " print*,  ' build date    : $(BUILID)     '"      > version.dat
	@echo     " print*,  ' git version   : $(GIT_VERSION)'"      >> version.dat
	@echo     " print*,  ' git repo      : https://github.com/hokru/xopt.git'"      >> version.dat


# make the MPI gradient helper binary
pgrad:
	$(MAKE) -C pgrad/

# packing all good things by hand
tarball:
	tar -czf xopt.tar.gz *.f90 Makefile lib/timestamp.f90  Manual/man.pdf lib/wregrad.f90 pgrad/mpigrad.f90 pgrad/Makefile
	zip xopt.zip *.f90 Makefile lib/timestamp.f90 Manual/man.pdf lib/wregrad.f90 pgrad/mpigrad.f90 pgrad/Makefile

# zip the git
archive:
	git archive master --format=zip --output=../xopt-$(GIT_VERSION).zip

clean:
	rm -f $(BUILDDIR)/* $(MOD)/* $(PROG) xopt.tar.gz xopt.zip version.dat

############################### ENDE ##############################################
