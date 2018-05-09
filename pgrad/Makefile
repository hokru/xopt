  PROG =  xopt.pgrad
  OBJS =  mpigrad.o 

#--------------------------------------------------------------------------
# set to "yes" here or on use: make MPI=yes
MPI=yes

#--------------------------------------------------------------------------
        
# gfortran
ifeq ($(MPI),yes)
$(info ************  MPI VERSION ************)
 FC = mpif90 -ffree-line-length-none -m64 -static-libgfortran
 FC+= -Dpara
else
$(info ************  SERIAL VERSION ************)
 FC = gfortran -fopenmp -ffree-line-length-none -m64 -static 
endif


#LIBS = -L$(OPENBLAS) -lopenblas 
#LIBS= $(OPENBLAS)/libopenblas.a  
FFLAGS =  -cpp  -O
#FFLAGS =  -cpp -fbacktrace -O0 -fimplicit-none -fbounds-check


# make targets:
.PHONY: all
.PHONY: clean
.PHONY: version
.PHONY: tarball


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

BUILID:=$(shell date)
#
all:$(PROG)

$(PROG):$(OBJS) 
	$(FC) $(OBJS) $(LIBS) -o $(PROG)

clean:
	rm -f *.o *.mod $(PROG) 

