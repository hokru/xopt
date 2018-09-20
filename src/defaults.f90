subroutine loaddef
use parm
use popt
use logic
use progs
use MDdat
use internals, only : int_LJ_cut
use fiso, only: r8
implicit none

! defaults
thrR=5.0_r8

iiter=1

cell=0_r8
molvol=0_r8
debug=.false.
large=.false. 
d3hess=.false.

mode=1  ! 1= min ; >2 TS
! TS
tsopt = .false.
tsmode = 0

!conv criteria
econv = 1e-7_r8 ! SCF energy change
gconv = 1.e-3_r8 ! GNORM threshold
maxgrad = 1e-3_r8     ! MAX grad component
dconv = 0.005_r8 ! MAX displacment threshold

maxiter = 500 ! max iterations
maxd = 0.35_r8 ! max displacement


idiis=5
do_gdiis=.false.

!penalty function CI opt
ciopt=.false.
alpha=0.025_r8
sigma=3.5_r8

! opt options
iupdate=1 ! SR1-BFGS
iopt = 21 ! ANC+RFO
iconv = 4 ! all 4 conv crit must be reached for convergence

qmmm=.false.
readonly=.false.
readHess=.false.
tightopt=.false.
preopt=.false.
numgrad=.false.

freeze=.false.
restrain=.false.
restart=.false.

apbs=.false.
mopac=.false.
tm=.false.
orca=.false.
gaus=.false.
hmodel='simple'
openbabel=.false.
tmhuge=.false.
driver=.false.
psi4=.false.
gamess=.false.
riper=.false.
nwchem=.false.
xtb=.false.
GEI=.false.

Doshift=.true.

int_LJ_cut=7.0_r8

do_dvv=.false.
dvv_err=0.003_r8
dvv_step=0.04_r8 
dvv_init=1.0_r8 ! forward direction
dvv_de=1e-8_r8
dvv_up=5e-6_r8

orient=.false.

! parallel
nproc=1


usrscr="local"
scratchjob=.false.

! PROGRAM SYSTEM CALLS
scall_orca='orca40'
scall_xtb='xtb'
scall_nwchem='nwchem.x $PARNODES nw.in nw.out'
scall_gms='rungms'//' '//trim(gmsin)//' $PARNODES'
scall_mopac='mopac16'
scall_psi4='psi4'
scall_gaus='run-g09'
scall_amber='sander -O -i '//trim(ambin)// &
             ' -c amber.rst'// &
             ' -p '//trim(ambtop)// &
!             ' -o '//trim(xjob)
             ' -o xopt.job 2> xopt.err '
scall_gei='mygrad.sh'
scall_mpi='mpiexec -bynode'
prog_flags=''

geigrad='xopt.grad'

exlib=.false.

! MD section
do_md=.false.
nvar=3*nat
dt=0.5_r8
temp0=298.15_r8
maxstep=10000 ! = 5ps
iprint=1
ifile=1
thermo_nh=.true.
thermo_berend=.false.
thermo_scale=.false.


do_hmass=.false.

! catch & fix default for special case
if(nat>3000) then
call message('using unit Hessian for nat>3000, can be overwritten by request')
 hmodel='hdiag'
endif


end subroutine loaddef
