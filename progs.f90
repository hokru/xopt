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
subroutine prepPROG
use logic
use progs
implicit none
character(120) aa
integer i

!****************
!* ORCA         *
!****************
if(orca) then
print*, " * THIS IS AN ORCA OPTIMIZATION *"
call system('rm orca.gbw')

!i=0
!open(unit=13,file=orcain,status='old')
!do
! i=i+1
! read(13,'(a)',end=111) a(i)
!enddo
!111 continue
!close(130

!open(unit=13,file=orcain,status='old',position='append')
!write(13,'(a)')'!ENGRAD '
!write(13,'(a,a)')'*xyzfile 0 1 ', trim(xyzfile)
!close(13)
return
endif

if(openbabel) then
print*, " * THIS IS AN OPENBABEL(GAFF) OPTIMIZATION *"
! openbabel currently(?) cannot print "stand-alone" gradients
return
endif


!************************
!* QM/MM                *
!************************
if(qmmm) then
print*, " * THIS IS A QM/MM OPTIMIZATION (good luck) *"
! anything to prepare?
return
endif

if(numgrad) then
print*, " * THIS IS A NUMERICAL GRAD OPTIMIZATION (good luck) *"
return
endif

!************************
!* DRIVER (Grimme)      *
!************************
if(driver) then
print*, " * THIS IS AN DRIVER OPTIMIZATION *"
print*,' *** '
print*,' driver input: '
open(55,file='.DRIVERINPUT')
do
 read(55,'(a)',end=551) aa
 print*,trim(aa)
enddo
551 close(55)
print*,' *** '
return
endif


!****************
!* MOPAC12      *
!****************
if(mopac) then
print*, " * THIS IS AN MOPAC OPTIMIZATION *"
return
endif

!****************
!* TM           *
!****************
if(tm.or.tmhuge) then
print*, " * THIS IS A TURBOMOLE OPTIMIZATION *"
 call system('actual -r >/dev/null')
 tmri=.false.
 print*,' Checking control file for $rij'
 open(111,file='control')
  do
  read(111,'(a)',end=666) aa
    if(index(aa,'$rij').ne.0) then
      tmri=.true.
      write(6,'(5x,a)') ' -> found $rij option'
      exit
    endif
  enddo
 666 continue
 close(111)
return
endif



!****************
!* gaussian09   *
!****************
if(gaus) then
print*, " * THIS IS A GAUSSIAN OPTIMIZATION *"
return
endif


if(amber) then
!call system('cp amber.crd amber.rst')
print*, " * THIS IS AN AMBER OPTIMIZATION *"
return
endif

call error(6,'select prog for energy and gradient!')

end subroutine prepPROG
