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
subroutine loaddef
use parm
use popt
use logic
implicit none

! defaults
thrR=5d0


mode=1  ! 1= min ; >2 TS
! TS
doTS = .false.
iTS = 1

!conv criteria
econv = 1e-7 ! SCF energy change
gconv = 1.e-3 ! GNORM threshold
maxgrad = 1e-3     ! MAX grad component
dconv = 0.005 ! MAX displacment threshold

maxiter = 500 ! max iterations
maxd = 0.35d0 ! max displacement


! opt options
iupdate=1 ! SR1-BFGS
iopt = 21 ! RFO
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
hmodel='hdiag'
openbabel=.false.
tmhuge=.false.
driver=.false.

Doshift=.true.


return

end subroutine loaddef
