module fiso
! wrapper around iso_fortran_env
#ifdef GNU
 use, intrinsic :: iso_fortran_env, only : &
     stdin=> input_unit, stdout=>output_unit,stderr=>error_unit, & ! I/O units, replaces non-def 0,5,6
      i8=>INT64,i4=>INT32               , & 
      compiler_version                            ! nerdy info
! IEEE 754
 integer, parameter :: r4 = SELECTED_REAL_KIND(6,37)
 integer, parameter :: r8 = SELECTED_REAL_KIND(15,307)
 integer, parameter :: qp = SELECTED_REAL_KIND(33,4931)

 integer, parameter :: io_debug=1111
#elif GNU_LEGACY # 4.5
 integer, parameter :: r4 = SELECTED_REAL_KIND(6,37)
 integer, parameter :: r8 = SELECTED_REAL_KIND(15,307)
 integer, parameter :: qp = SELECTED_REAL_KIND(33,4931)
 integer, parameter :: io_debug=1111
 integer, parameter :: i4 = selected_int_kind(4)
 integer, parameter :: i8 = selected_int_kind(8)
 integer, parameter :: stdin =0
 integer, parameter :: stdout =6
 integer, parameter :: stderr= 5
 contains
 character(7) function compiler_version()
 implicit none
  compiler_version='legacy gnu'
 end function
#else 
! for cases that do not have yet compiler_version like ifort 15 
  use, intrinsic :: iso_fortran_env, only : &
       stdin=> input_unit, stdout=>output_unit,stderr=>error_unit & ! I/O units, replaces non-def 0,5,6
       ,i8=>INT64,i4=>INT32   
! IEEE 754
 integer, parameter :: r4 = SELECTED_REAL_KIND(6,37)
 integer, parameter :: r8 = SELECTED_REAL_KIND(15,307)
 integer, parameter :: qp = SELECTED_REAL_KIND(33,4931)

 integer, parameter :: io_debug=1111
 contains
 character(7) function compiler_version()
 implicit none
  compiler_version='unknown'
 end function
#endif

!real(r8), parameter :: eins=1.0_r8,nul=0.0_r8
!note:
! The iso fortran literals REAL64 etc only declared _storage size_, not precision
! Stevel Lionel (Dr Fortran) recommends selected_xx_kind
! https://software.intel.com/en-us/blogs/2017/03/27/doctor-fortran-in-it-takes-all-kinds
end module fiso




module constant
use fiso, only: r8
! in part taken from psi4
real(r8), parameter:: pi = 3.141592653589793_r8  ! cake
real(r8), parameter:: au2ang = 0.52917720859_r8
real(r8), parameter:: amu2au=1.66053886E-27_r8/9.10938215E-31_r8
real(r8), parameter:: au2cm =219474.63067_r8
real(r8), parameter:: au2kcal = 627.5095_r8
real(r8), parameter:: au2cm1 = 219474.6_r8
real(r8), parameter:: au2ev= 27.21138_r8
real(r8), parameter:: au2mhz=6.579684E9_r8

real(r8), parameter:: kb_J=1.3806504E-23_r8 ! J/K
real(r8), parameter:: kB_au=3.1668114E-6_r8 ! Eh/K
real(r8), parameter:: au2fs=0.02418884326505_r8

real(r8), parameter :: Planck= 6.62606896D-34       ! The Planck constant (Js)
real(r8), parameter :: AvoN= 6.02214179E23_r8       ! Avagadro's number
real(r8), parameter :: pc_c= 2.99792458E8_r8        ! Speed of light (ms$^{-1}$)
real(r8), parameter :: bohr2m= 0.52917720859E-10_r8 !   Bohr to meters conversion factor
real(r8), parameter :: amu2kg= 1.660538782E-27_r8   !  Atomic mass units to kg conversion factor
real(r8), parameter :: hartree2J= 4.359744E-18_r8   ! hartree to Joule
real(r8), parameter :: au2m = au2ang/(10.0_r8**10)  ! bohr to metre 

real(r8), parameter :: w8Gb = 8.0_r8/(1024.0_r8**3)  ! 64bit word to Gb
end module

module timing
use fiso, only: r8,i8
integer(i8) cr ! count rate
!integer cr ! count rate
real(r8) tstart,elapsed
real(r8) tend,timer
end module

module parm
use fiso, only: r8
! contains common arrays and variables
  integer, parameter :: maxat = 50000 ! max. atoms
  integer, parameter :: maxrestr = 500 ! max. number of restrains
  integer i,j,k,l
  real(r8) thrR        ! distance threshold, 50 bohr
  integer nat  ! number of atoms
  integer nvar  ! number of variables to opt
  integer nat3  ! number of atoms times 3
  integer npair ! number of atom paris
  real(r8) energy
  real(r8) gnorm
  integer, allocatable :: iat(:)

  real(r8), allocatable :: xyz(:,:) ! coordintes
  real(r8), allocatable :: xyz0(:,:) ! start coordintes
  real(r8), allocatable :: dist(:)  ! distances for all pairs
  integer, allocatable :: kj(:),ki(:) ! i,j index for dist vector

  ! integer :: ifrag(:) ! fragment optimization

  real(r8), allocatable ::  chess(:,:) ! cartesian hessian
  real(r8), allocatable ::  hint(:)    ! internal hessian nvar*(nvar+1)/2)
  real(r8), allocatable ::  B(:,:)     ! Wilson B-matrix
  real(r8), allocatable ::  grad(:,:) ! cartesian gradient
  real(r8), allocatable ::  gint(:) ! internal gradient
! 
  real(r8) :: cell(3),molvol ! molecular cell size and effective(=restrained) volume

! cartesian constraints
  integer, allocatable :: ifrez(:)
! primitive restraints
integer ires
real(r8), allocatable:: irest_vol(:,:)
integer, allocatable:: irest_atom(:)
integer, allocatable:: irest_bond(:,:)
integer, allocatable:: irest_ang(:,:)
integer, allocatable:: irest_dihed(:,:)
real(r8), allocatable:: irest_konst(:,:)
real(r8), allocatable :: val0(:) ! restart values for restrains
real(r8), allocatable :: eneR(:) ! restraint energies (E_h)
real(r8) eneTOT ! total energy restraints (E_h)


end module

!************************************************************

module logic
logical echo ! printout
logical debug
logical relax !no calcs, just relax structure
logical exlib ! extralibs
logical large
logical d3hess
! programs
logical readonly
logical gradonly
logical mopac
logical tm,tmri,tmcc,tmhuge,riper
logical orca,xtb
logical gaus
logical amber,apbs
logical openbabel
logical driver
logical qmmm
logical numgrad
logical ciopt
logical psi4
logical gamess
logical nwchem
logical GEI ! general external interfacea
! add-ons
logical d3,gcp,ppot
character(80) d3func,gcplevel
! model hessian
character(20) hmodel
! shift hessian
logical Doshift
! parallel
integer nproc ! number of processes for mpicall to external programs
integer nomp ! number of threads for BLAS/LAPACK
logical scratchjob
! PA frame
logical orient
! MD
logical do_md
! reduced precision for speed
logical do_float

! internals
logical int_deloc

end module logic

!************************************************************

module popt
use fiso, only: r8
! parameter for optimization and other options

integer iiter

logical readHess
logical restart
character(80) hessname

integer mode
! TS opt:
logical tsopt   ! do TS search
integer tsmode  ! eigenvector to be maximized
! re/con-straints
logical freeze
logical restrain
logical constrain
! composite keywords
logical tightopt,preopt
! conv criteria
real(r8) econv  ! SCF energy change
real(r8) gconv  ! GNORM threshold
real(r8) maxgrad     ! MAX grad component
real(r8) dconv     ! MAX displ theshold

real(r8) Hmin   ! minimal Hessian value

integer maxiter ! max iterations
real (r8) maxd  ! max displacement

integer iupdate ! 1=SR1BFGS,2=SR1,3=BFGS
integer iopt  ! 11=cart-RFO 12=cart-SI-RFO  21=ANC+RFO 22=ANC+SI-RFO
integer iconv ! number of conv criteria fullfilled =4 (all),=2 (dE+G),=1 (dE)

! RFO stuff
real(r8) sdiag  ! scaling parm for SI-RFO
real(r8) soff ! scaling parm for SI-RFO
! GDIIS
 integer idiis
logical do_gdiis
! ciopt
real(r8) alpha, sigma
! DDV (IRC via MD)
logical do_dvv
real(r8) dvv_err,dvv_step ! 0.003 and 0.04/0.08 in paper
! initial velocity direction to chose educt/product: 1 or -1
  ! /= 1 to damp/increase the initial step.
real(r8) dvv_init
real(r8) dvv_de,dvv_up

end module


module progs
! handle input names for programs and therelikes
! ToDo: allow more user changes
character(255) :: scall_orca,command_orca
character(255) :: scall_xtb,command_xtb
character(255) :: scall_nwchem,command_nwchem
character(255) :: scall_gms,command_gms
character(255) :: scall_mopac,command_mopac
character(255) :: scall_psi4,command_psi4
character(255) :: scall_gaus,command_gaus
character(255) :: scall_amber,command_amber
character(255) :: scall_gei,command_gei
character(255) :: scall_mpi,command_mpi
character(255) :: prog_flags

character(255) :: geigrad

character(3), parameter :: redir =' > '
character(6), parameter :: errout =' 2>&1 '
character(80), parameter :: xjob = 'xopt.job  2>&1 '
character(80), parameter :: orcain = 'orca.in'
character(80), parameter :: xyzfile = 'xopt.xyz'

character(80), parameter :: tmin = 'coord'
character(80), parameter :: mopin = 'mopac.in'
character(80), parameter :: gauin = 'gau.in'
character(80), parameter :: gmsin = 'gms' ! file actually gms.inp

character(80), parameter :: ambin = 'amber.in'
character(80), parameter :: ambcrd = 'amber.crd'
character(80), parameter :: ambtop = 'amber.top'

character(200) binpath

character(200) scrdir,workdir,usrscr
end module

module intcoords
use fiso
integer nints ! number of internal coordinates
integer nbond,nang,ndih ! number of bonds,angles,dihedrals
real(r8), allocatable :: bond(:)
real(r8), allocatable :: ang(:)
real(r8), allocatable :: dih(:)
end module

!module parallel
!use mpi
!integer ierr,nproc,rank,strlen
!character(MPI_MAX_PROCESSOR_NAME) procname
!end module


module mod_qmmm
integer  imod,ireal
integer, allocatable :: isys(:)
character(120) nreal,nlow,nhigh,namber
logical doEE,debug_qmmm
end module


module atomdata
use fiso, only: r8
real(r8) ams(107)
real(r8) ams2(118)
!real(8) ams(118)

data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,      &
    10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0, &
    20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0, &
    30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0, &
    40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0, &
    54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0, &
    65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0, &
    79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0, &
    91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0, &
    102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0, &
    118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0, &
    132.9054d0, 137.3300d0, 15*0.000d0, 178.4900d0, 180.9479d0, &
    183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0, &
    196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0, &
    18*0.000d0,   0.0000d0,  5*0.000d0/


! taken from PSI4
data ams2 /1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406, &
12,14.00307400478,15.99491461956,18.998403224,19.99244017542, &
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629, &
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983, &
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,  &
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,  &
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,  &
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,  &
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292 , &
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676, &
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932, &
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297, &
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757, &
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089, &
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109, &
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011, &
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271, &
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325, &
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108, &
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724, &
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500, &
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451, &
285.183698,287.191186,292.199786,291.206564,293.214670/


character(2) elem(95)
DATA elem/'h ','he',                                             &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu','xx'/ ! elem(95) being a dummy


real(r8) rvdw(94),rcov(94)
! atomic radii from Mantina, Valero, Cramer, Truhlar "Atomic radii of elements"
! Copyed by hand.
!            H       He
data rvdw /1.10d0,1.40d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.82d0,1.53d0,1.92d0,1.70d0,1.54d0,1.52d0,1.47d0,1.54d0, &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    2.27d0,1.73d0,1.84d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0, &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.75d0,2.31d0,2.15d0,2.11d0,2.07d0,2.06d0,2.05d0,2.04d0,2.00d0,1.97d0,1.96d0,2.01d0,1.87d0,2.11d0,1.85d0,1.90d0,1.85d0,2.02d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    3.03d0,2.49d0,2.26d0,2.23d0,2.18d0,2.17d0,2.16d0,2.13d0,2.10d0,2.10d0,2.11d0,2.18d0,1.93d0,2.17d0,2.06d0,2.06d0,1.98d0,2.16d0, &
    ! Cs Ba
    3.32d0,2.68d0, &
    ! La-Lu
    2.43d0,2.42d0,2.40d0,2.46d0,2.38d0,2.36d0,2.35d0,2.34d0,2.33d0,2.31d0,2.30d0,2.29d0,2.27d0,2.26d0,2.24d0, &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    2.23d0,2.22d0,2.18d0,2.16d0,2.16d0,2.13d0,2.13d0,2.23d0,2.23d0,2.11d0,2.02d0,2.07d0,1.97d0,2.02d0,2.20d0, &
    ! Fr-Pu
    3.48d0,2.83d0,2.47d0,2.45d0,2.43d0,2.41d0,2.39d0,2.43d0/

! in Angstrom
data rcov /0.32d0,0.37d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.30d0,0.99d0,0.84d0,0.75d0,0.71d0,0.64d0,0.60d0,0.62d0,  &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    1.60d0,1.40d0,1.24d0,1.14d0,1.09d0,1.04d0,1.00d0,1.01d0,  &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.00d0,1.74d0,1.59d0,1.48d0,1.44d0,1.30d0,1.29d0,1.24d0,1.18d0,1.17d0,1.22d0,1.20d0,1.23d0,1.20d0,1.20d0,1.18d0,1.17d0,1.24d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    2.15d0,1.90d0,1.78d0,1.64d0,1.56d0,1.46d0,1.38d0,1.36d0,1.34d0,1.30d0,1.36d0,1.40d0,1.42d0,1.40d0,1.40d0,1.37d0,1.32d0,1.36d0, &
    ! Cs Ba
    2.38d0,2.06d0,  &
    ! La-Lu
    1.94d0,1.84d0,1.90d0,1.73d0,1.86d0,1.85d0,1.83d0,1.82d0,1.81d0,1.80d0,1.79d0,1.77d0,1.77d0,1.78d0,1.74d0,  &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    1.64d0,1.58d0,1.50d0,1.41d0,1.36d0,1.32d0,1.30d0,1.64d0,1.88d0,1.48d0,1.45d0,1.50d0,1.42d0,1.47d0,1.46d0,  &
    ! Fr-Pu
    2.42d0,2.11d0,2.01d0,1.90d0,1.84d0,1.83d0,1.80d0,1.80d0/

! Covalent radii revisited", Dalton Trans., 2008, 2832. Cordero et al.
real(r8) rcov2(96)
data rcov2 / &
0.31, & !H
0.28, & !He
1.28, & !Li
0.96, & !Be
0.84, & !B
0.76, & !C
0.71, & !N
0.66, & !O
0.57, & !F
0.58, & !Ne
1.66, & !Na
1.41, & !Mg
1.21, & !Al
1.11, & !Si
1.07, & !P
1.05, & !S
1.02, & !Cl
1.06, & !Ar
2.03, & !K
1.76, & !Ca
1.70, & !Sc
1.60, & !Ti
1.53, & !V
1.39, & !Cr
1.61, & !Mn
1.52, & !Fe
1.50, & !Co
1.24, & !Ni
1.32, & !Cu
1.22, & !Zn
1.22, & !Ga
1.20, & !Ge
1.19, & !As
1.20, & !Se
1.20, & !Br
1.16, & !Kr
2.20, & !Rb
1.95, & !Sr
1.90, & !Y
1.75, & !Zr
1.64, & !Nb
1.54, & !Mo
1.47, & !Tc
1.46, & !Ru
1.42, & !Rh
1.39, & !Pd
1.45, & !Ag
1.44, & !Cd
1.42, & !In
1.39, & !Sn
1.39, & !Sb
1.38, & !Te
1.39, & !I
1.40, & !Xe
2.44, & !Cs
2.15, & !Ba
2.07, & !La
2.04, & !Ce
2.03, & !Pr
2.01, & !Nd
1.99, & !Pm
1.98, & !Sm
1.98, & !Eu
1.96, & !Gd
1.94, & !Tb
1.92, & !Dy
1.92, & !Ho
1.89, & !Er
1.90, & !Tm
1.87, & !Yb
1.87, & !Lu
1.75, & !Hf
1.70, & !Ta
1.62, & !W
1.51, & !Re
1.44, & !Os
1.41, & !Ir
1.36, & !Pt
1.36, & !Au
1.32, & !Hg
1.45, & !Tl
1.46, & !Pb
1.48, & !Bi
1.40, & !Po
1.50, & !At
1.50, & !Rn
2.60, & !Fr
2.21, & !Ra
2.15, & !Ac
2.06, & !Th
2.00, & !Pa
1.96, & !U
1.90, & !Np
1.87, & !Pu
1.80, & !Am
1.69/   !Cm


! CSD atomic radii + .20 (from molden)
real(r8) rcov3(100)
    data rcov3/0.430d0,0.741d0,0.880d0,0.550d0,1.030d0,0.900d0,0.880d0,&
               0.880d0,0.840d0,0.815d0,1.170d0,1.300d0,1.550d0,1.400d0,&
               1.250d0,1.220d0,1.190d0,0.995d0,1.530d0,1.190d0,1.640d0,&
               1.670d0,1.530d0,1.550d0,1.555d0,1.540d0,1.530d0,1.700d0,&
               1.720d0,1.650d0,1.420d0,1.370d0,1.410d0,1.420d0,1.410d0,&
               1.069d0,1.670d0,1.320d0,1.980d0,1.760d0,1.680d0,1.670d0,&
               1.550d0,1.600d0,1.650d0,1.700d0,1.790d0,1.890d0,1.830d0,&
               1.660d0,1.660d0,1.670d0,1.600d0,1.750d0,1.870d0,1.540d0,&
               2.070d0,2.030d0,2.020d0,2.010d0,2.000d0,2.000d0,2.190d0,&
               1.990d0,1.960d0,1.950d0,1.940d0,1.930d0,1.920d0,2.140d0,&
               1.920d0,1.770d0,1.630d0,1.570d0,1.550d0,1.570d0,1.520d0,&
               1.700d0,1.700d0,1.900d0,1.750d0,1.740d0,1.740d0,1.880d0,&
               0.200d0,0.200d0,0.200d0,2.100d0,2.080d0,1.990d0,1.810d0,&
               1.780d0,1.750d0,0.200d0,1.710d0,0.200d0,0.200d0,1.730d0,&
               0.100d0,0.200d0/


end module


module MDdat
use fiso, only: r8
real(r8) dt ! time step (will be converted to au from fs)
real(r8), allocatable:: mass(:)   ! atomic masses in au after setMass
real(r8), allocatable:: velo(:,:) ! atomic velocity
real(r8) temp0  ! target temperature in K
integer(r8) maxstep
integer ifile,iprint


!thermostats
logical  do_thermo
real(r8) tau  ! berendsen thermostat
logical thermo_nh, thermo_berend,thermo_scale
real(r8) m_nh(2) ! masses
real(r8) p_nh(2) ! positions
real(r8) v_nh(2) ! velocities
real(r8) c_nh(2) ! coupling to bath
real(r8) g_nh(2)
real(r8) r_nh(2)
end module


module internals
use fiso, only: r8
! internal arrays
integer int_nb,int_na,int_nt,nints
real(r8), allocatable :: int_bval(:)
integer, allocatable :: int_bcast(:,:)
real(r8), allocatable :: int_aval(:)
integer, allocatable :: int_acast(:,:)
real(r8), allocatable :: int_tval(:)
integer, allocatable :: int_tcast(:,:)

! LJ-type bond internals
integer :: int_nlj
real(r8) :: int_LJ_cut
real(r8), allocatable :: int_ljval(:)
integer, allocatable :: int_ljcast(:,:)

! connect fragments
integer frag_nh, frag_nb
real(r8), allocatable :: frag_hval(:) ! H-bonds
integer, allocatable :: frag_hcast(:,:)
real(r8), allocatable :: frag_bval(:) ! non-bonded distances
real(r8), allocatable :: frag_bcast(:,:)
! real(r8), allocatable :: frag_aval(:,:) ! non-bonded angles
! real(r8), allocatable :: frag_acast(:,:)

real(r8), allocatable :: hdiag(:)

! transformator and projectors
real(r8), allocatable :: Gmat(:,:)
real(r8), allocatable :: Umat(:,:)
real(r8), allocatable :: Bdeloc(:,:)
real(r8), allocatable :: Ginv(:,:) 
! real(r8), allocatable :: P_int(:,:) 
! real(r8), allocatable :: A_int(:,:)
real(r8), allocatable :: inthess(:,:)

end module internals

