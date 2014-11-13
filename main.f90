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
!*************************
!* An eXternal OPTimizer *
!*************************

program  xopt
use parm  ! essential parameter
use logic ! essential logic
use popt 
use progs
implicit none
integer nbond,istat,iiter
real(8) Ebond,Ebend,Edisp,Erepul,e
character(120) infile,outfile
logical debug

! load defaults
call loaddef
iiter=1
debug=.false.


! some print out
call head(echo)

call rcontrol
call eval_options(infile,outfile)

call tmolrd(infile,.false.,.true.) ! just determined nat
if(nat.gt.maxat) print*,'too many atoms, increase: maxat =',maxat
call p_info(echo)

npair=(nat*(nat-1))/2
nat3=nat*3

allocate(xyz(3,nat),stat=istat)
allocate(xyz0(3,nat),stat=istat)
allocate(grad(3,nat),stat=istat)
allocate(iat(nat),stat=istat)


if(freeze) allocate(ifrez(nat),stat=istat)
if(restrain) then
 allocate(irest_bond(maxrestr,2))
 allocate(irest_ang(maxrestr,3))
 allocate(irest_dihed(maxrestr,4)) 
 allocate(irest_konst(maxrestr,3)) 
 allocate(val0(maxrestr)) 
 allocate(eneR(maxrestr)) 
endif

allocate(chess(3*nat,3*nat))

if(istat.ne.0) stop 'allocation error @ startup'

call tmolrd(infile,.true.,.false.) ! read coordinates
if(debug) print*,'tmolrd done'
xyz0=xyz


! calculate restraining energy from restart file
if(ReadOnly) then
restart=.true.
call rcontrolrestrain
call rdhess(nat3,chess,'xopt.restart')
call addRestrainGrad
call printFinalRestrain
print*,'Reading done..stopping'
stop
endif

call rcontrolfreeze
call rcontrolrestrain



call wrxyz(trim(xyzfile))
call prepPROG
if(debug) print*,'prepPROG done'
!call get_dist(echo) ! distances of all unique pairs

! make approx. hessian
call getHess
if(debug) print*,'getHess done'
! form coordinate system

!

if(readHess) call rdhess(nat3,chess,'xopt.chess')
if(restart)  call rdhess(nat3,chess,'xopt.restart')
if(debug) print*,'rdhess done'
!  call readORCAhess(nat3,chess,'xopt.chess')

! enter TS search


! enter minimizer (default)
select case(iopt)
 case(11:19)
   if(restrain) then
    call copt
   else
    call copt
   endif
 case(21:29)
  call getanc
  if(restrain) then
    call ancopt(iiter)
  else
    call ancopt(iiter)
  endif
 case(31:39)
  nvar=npair
  call int_bonds_B
  call ancopt(iiter)
end select

call wrxyz('xopt.xyz')

call wrxyz('xopt.opt')

if(restrain) then
 open(unit=33,file='xopt.restrain.tmp',status='old')
 close(33,status='delete')
endif

end program xopt
