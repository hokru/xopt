!**************************************************************************************
!    Copyright 2013, 2014 Holger Kruse                                                *
!                                                                                     *
!    This file is part of Xopt.                                                       *
!                                                                                     *
!    Xopt is free software: you can redistribute it and/or modify                     *
!    it under the terms of the GNU Lesser General Public License as published by      *
!    the Free Software Foundation, either version 3 of the License, or                *
!    (at your option) any later version.                                              *
!                                                                                     *
!    xopt is distributed in the hope that it will be useful,                        *
!    but WITHOUT ANY WARRANTY; without even the implied warranty of                   *
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    *
!    GNU Lesser General Public License for more details.                              *
!                                                                                     *
!    You should have received a copy of the GNU Lesser General Public License         *
!    along with Xopt.  If not, see <http://www.gnu.org/licenses/>.                    *
!                                                                                     *
!**************************************************************************************
!************************************************************
module constant
real(8), parameter:: pi = 3.141592653589793d0
real(8), parameter:: au2ang = 0.5291770d0
end module


module parm
! contains common arrays and variables
  integer, parameter :: maxat = 50000 ! max. atoms
  integer, parameter :: maxrestr = 500 ! max. number of restrains
  integer i,j,k,l
!  real(8), parameter :: thrR = 50        ! distance threshold, 50 bohr
  real(8) thrR        ! distance threshold, 50 bohr
  integer nat  ! number of atoms
  integer nvar  ! number of variables to opt
  integer nat3  ! number of atoms times 3
  integer npair ! number of atom paris
  real(8) energy 
  real(8) gnorm
  integer, allocatable :: iat(:)

  real(8), allocatable :: xyz(:,:) ! coordintes
  real(8), allocatable :: xyz0(:,:) ! start coordintes
  real(8), allocatable :: dist(:)  ! distances for all pairs
  integer, allocatable :: kj(:),ki(:) ! i,j index for dist vector

  real(8), allocatable ::  chess(:,:)
  real(8), allocatable ::  hint(:)
  real(8), allocatable ::  B(:,:)
  real(8), allocatable ::  grad(:,:) ! cartesian gradient
  real(8), allocatable ::  gint(:) ! internal gradient



! cartesian constraints
  integer, allocatable :: ifrez(:)
! primitive restraints
integer ires
integer, allocatable:: irest_bond(:,:)
integer, allocatable:: irest_ang(:,:)
integer, allocatable:: irest_dihed(:,:)
real(8), allocatable:: irest_konst(:,:)
real(8), allocatable :: val0(:) ! restart values for restrains
real(8), allocatable :: eneR(:) ! restraint energies (E_h)
real(8) eneTOT ! total energy restraints (E_h)

end module

!************************************************************

module logic
logical echo ! printout
logical relax !no calcs, just relax structure

! programs
logical readonly
logical mopac
logical tm,tmri,tmcc,tmhuge
logical orca
logical gaus
logical amber,apbs
logical openbabel
logical driver
logical qmmm
logical numgrad

logical d3,gcp,ppot
character(80) d3func,gcplevel

! model hessian
character(20) hmodel

! shift hessian
logical Doshift
end module logic

!************************************************************

module popt
! params for optimization and other options

logical readHess
logical restart

integer mode
! TS opt:
logical doTS
integer iTS  ! 


logical freeze  ! 
logical restrain
logical constrain

logical tightopt,preopt
!conv criteria
real(8) econv  ! SCF energy change
real(8) gconv  ! GNORM threshold
real(8) maxgrad     ! MAX grad component
real(8) dconv     ! MAX displ theshold

real(8) Hmin   ! minimal Hessian value

integer maxiter ! max iterations
real (8) maxd  ! max displacement

integer iupdate ! 1=SR1BFGS,2=SR1,3=BFGS
integer iopt  ! 11=cart-RFO 12=cart-SI-RFO  21=ANC+RFO 22=ANC+SI-RFO
integer iconv ! number of conv criteria fullfilled =4 (all),=2 (dE+G),=1 (dE)

! RFO stuff
real(8) sdiag  ! scaling parm for SI-RFO
real(8) soff ! scaling parm for SI-RFO
! GDIIS
 integer, parameter :: idiis = 10

end module



module progs
! handle input names for programs and therelikes

character(80), parameter :: xjob = 'xopt.job  2>&1 '
character(80), parameter :: orcain = 'orca.in'
character(80), parameter :: xyzfile = 'xopt.xyz'

character(80), parameter :: tmin = 'coord'
character(80), parameter :: mopin = 'mopac.in'
character(80), parameter :: gauin = 'gau.in'
character(80), parameter :: gammin = 'gam.in'


character(80), parameter :: ambin = 'amber.in'
character(80), parameter :: ambcrd = 'amber.crd'
character(80), parameter :: ambtop = 'amber.top'


end module

module intcoords

integer nints ! number of internal coordinates
integer nbond,nang,ndih ! number of bonds,angles,dihedrals
real(8), allocatable :: bond(:)
real(8), allocatable :: ang(:)
real(8), allocatable :: dih(:)

end module



module mod_qmmm
integer  imod,ireal
integer, allocatable :: isys(:)
character(120) nreal,nlow,nhigh,namber
logical doEE,debug
end module
