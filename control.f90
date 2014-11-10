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
!****************************************************
!* this manages various kinds of internal log files *
!****************************************************
!
! xopt.grad.log   all gradient of the optimization numbered with iteration and energy
! xopt.control    flags for the program
!
!
!
!



!********************
!* log the gradient *
!********************
subroutine loggrad(nat,iter,grad,energy)
implicit none
integer iter,nat,i
real(8) grad(3,nat),energy

open(unit=11,file='xopt.grad.log',status='old',position='append')
write(11,'(a,2x,I5,2x,F)') '$iter',iter, energy
do i=1,nat
 write(11,'(3F18.14)') grad(1,i),grad(2,i),grad(3,i)
enddo
close(11)
end subroutine

subroutine touchlog
implicit none
open(unit=11,file='xopt.grad.log',status='replace')
write(11,'(a)') 'logfile for the xopt gradients'
close(11)

open(unit=11,file='xopt.log',status='replace')
!write(11,'(a)') 'logfile for the xopt gradients'
close(11)

end subroutine touchlog




!*************************************
!* restart file containing:          *
!* last hessian                      *
!* restraining informations          *
!*************************************
subroutine wrestart
use parm
use logic
use progs
use popt
implicit none
character(120) aa

if(restrain) open(unit=33,file='xopt.restrain.tmp',status='old')
open(44,file='xopt.hess',status='old')

open(unit=11,file='xopt.restart',status='replace')
write(11,'(a)') 'This is the restart file for xopt'
write(11,'(a)') 'you may change values by hand   '
write(11,'(a)') 

if(restrain) then
! copy xopt.restrain.tmp into the xopt.restart
write(11,'(a)') '$restr'
do
 read(33,'(a)',end=333) aa
 write(11,'(a)') aa
enddo
333 close(33)
endif

! write freeze section


! write last hessian
do
 read(44,'(a)',end=444) aa
 write(11,'(a)') aa
enddo
444 close(44)

close(11)

end subroutine wrestart



! read xopt.control
! handle logical and basic things
! detailed reading (eg for restraints) is done separately 
subroutine rcontrol
use parm
use logic
use progs
use popt
character(80) aa,bb,cc
logical da
integer ii

inquire(file='xopt.control',exist=da)
if(da) then
  print*,''
  print*, ' * found control file: xopt.control ! *'
  print*,''
else
 return
endif


open(unit=11,file='xopt.control',status='old')
do
 read(11,'(a)',end=123) aa

! convergence
 if(index(aa,'$conv').ne.0) then
  do
    read(11,'(a)',end=123) aa
  if(index(aa,'econv').ne.0) read(aa,*) bb,econv
  if(index(aa,'gnorm').ne.0)read(aa,*) bb,gnorm
  if(index(aa,'maxgrad').ne.0)read(aa,*) bb,maxgrad
  if(index(aa,'max displ').ne.0)read(aa,*) bb,dconv
  if(index(aa,'displ').ne.0)read(aa,*) bb,maxd
  if(index(aa,'maxiter').ne.0)read(aa,*) bb,maxiter
  if(index(aa,'$').ne.0) then ! next item, one line back and exit loop
   backspace(11)
   exit
  endif
  enddo
 endif

 if(index(aa,'$freeze').ne.0) freeze=.true.
 if(index(aa,'$restr').ne.0) then
   restrain=.true.
   iconv=1
 endif


 if(index(aa,'$oniom').ne.0) then
   oniom=.true.

 endif

! if(index(aa,'').ne.0) then
!  do
!    read(11,'(a)',end=123) aa
!  if(index(aa,'$').ne.0) then
!   backspace(11)
!   exit
!  endif
!  enddo
! endif


enddo
123 continue
close(11)
end subroutine rcontrol

! read restraints and write the restart file
subroutine rcontrolrestrain
use parm
use logic
use progs
use popt
implicit none
character(80) aa,bb,cc
logical da
integer ii,jj,kk,ll
integer s2i
real(8) s2r,dih,dgrad,ang,dang,bl,di360,anggrad

if(restart) then
  open(unit=11,file='xopt.restart',status='old')
else
  inquire(file='xopt.control',exist=da)
  if(.not.da)  return
  open(unit=11,file='xopt.control',status='old')
endif

! save initial data
open(unit=33,file='xopt.restrain.tmp')

do
 read(11,'(a)',end=124) aa

 if(index(aa,'$restr').ne.0) then 
  print*,''
  if(restart) then
    print*, ' * reading restraints information: xopt.restart ! *'
  else
    print*, ' * reading restraints information: xopt.control ! *'
  endif
  print*,''
 restrain=.true.
 irest_bond=0
 irest_ang=0
 irest_dihed=0
 irest_konst=0.0d0
  do
    read(11,'(a)',end=123) aa

    if(index(aa,'bond').ne.0) then 
     ires=ires+1
     call charXsplit(aa,bb,2)
      irest_bond(ires,1)=s2i(bb)
     call charXsplit(aa,bb,3)
      irest_bond(ires,2)=s2i(bb)
     call charXsplit(aa,bb,4)
      irest_konst(ires,1)=s2r(bb)
    stop 'not tested!'
    write(*,'(2x,a,2I4,2x,a,F9.5)') 'bond',irest_bond(ires,1),irest_bond(ires,2),'konst.',irest_konst(ires,1)
    write(33,'(2x,a,2I4,2x,F9.5)') 'bond',irest_bond(ires,1),irest_bond(ires,2),irest_konst(ires,1)
    endif
    if(index(aa,'ang').ne.0) then 
     ires=ires+1
     call charXsplit(aa,bb,2)
      irest_ang(ires,1)=s2i(bb)
      ii=irest_ang(ires,1)
     call charXsplit(aa,bb,3)
      irest_ang(ires,2)=s2i(bb)
      jj=irest_ang(ires,2)
     call charXsplit(aa,bb,4)
      irest_ang(ires,3)=s2i(bb)
      kk=irest_ang(ires,3)
     call charXsplit(aa,bb,5)
      irest_konst(ires,2)=s2r(bb)   
     call angle(xyz0,ii,jj,kk,ang,dang)
    val0(ires)=ang
    write(*,'(2x,a,3I4,2x,a,F9.5,a,F7.2)') 'angle', ii,jj,kk,' | restrain',irest_konst(ires,2),' | value: ',dang
    write(33,'(2x,a,3I4,2x,F9.5,2x,F)') 'angle', ii,jj,kk,irest_konst(ires,2),ang
    endif
    if(index(aa,'dihed').ne.0) then 
     ires=ires+1
     call charXsplit(aa,bb,2)
      irest_dihed(ires,1)=s2i(bb)
      ii=irest_dihed(ires,1)
     call charXsplit(aa,bb,3)
      irest_dihed(ires,2)=s2i(bb)
      jj=irest_dihed(ires,2)
     call charXsplit(aa,bb,4)
      irest_dihed(ires,3)=s2i(bb)
      kk=irest_dihed(ires,3)
     call charXsplit(aa,bb,5)
      irest_dihed(ires,4)=s2i(bb)
      ll=irest_dihed(ires,4)
     call charXsplit(aa,bb,6)
      irest_konst(ires,3)=s2r(bb) 
     if(restart) then
         call charXsplit(aa,bb,7)
         dih=s2r(bb) 
         dgrad=di360(dih)
     else
         call dihed(xyz0,ii,jj,kk,ll,dih,dgrad)
     endif
     val0(ires)=dih
    write(*,'(2x,a,4I4,2x,a,F9.5,a,F7.2)') 'torsion ',ii,jj,kk,ll,'| restrain:',irest_konst(ires,3),' | value: ',dgrad
    write(33,'(2x,a,4I4,2x,F9.5,2x,F)') 'dihed ',ii,jj,kk,ll,irest_konst(ires,3),dih
    endif
    if(index(aa,'$').ne.0) then
     backspace(11)
     exit
    endif
  enddo
 endif

enddo
123 continue
  if (irest_ang(1,1).gt.0) print*,' angle restrain gradient is done numerically'
124 continue
close(11)
close(33)
  print*,''

end subroutine



subroutine rcontrolfreeze
use parm
use logic
use progs
use popt
character(80) aa,bb,cc
logical da
integer ii

inquire(file='xopt.control',exist=da)
if(da) then
else
 return
endif

open(unit=11,file='xopt.control',status='old')
do
 read(11,'(a)',end=123) aa

 if(index(aa,'$freeze').ne.0) then ! set up ifreeze vector
  print*,''
  print*, ' * reading frozen atoms information: xopt.control ! *'
  print*,''
 freeze=.true.
 ifrez=0
  do
    read(11,'(a)',end=123) aa
  if(index(aa,'Hopt').ne.0) then ! freeze all non-H 
   do i=1,nat
    if(iat(i).ne.1) ifrez(i)=1
   enddo
  endif
  if(index(aa,'atom').ne.0) then ! freeze atom number X
   read(aa,*) bb,ii
   ifrez(ii)=1
  endif
  if(index(aa,'elem').ne.0) then ! freeze elem X
   call charXsplit(aa,bb,2)
   call elem(bb,ii)
   do i=1,nat
    if(iat(i)==ii) ifrez(i)=1
   enddo
  endif
  if(index(aa,'$').ne.0) then
   backspace(11)
   exit
  endif
  enddo
 endif

enddo
123 continue
close(11)
end subroutine


subroutine ReadControlOniom
use parm
use logic
use progs
use popt
use mod_qmmm
implicit none
character(80) aa,bb,cc
logical da
integer ii,s2i

!* ireal: cart. size of real system
!* imod : cart. size of model system
!* isys : 0=real system 1=model
!*
!* 
!* OUTPUT:
!* energy+gradient (parm module)
!*
!*
!* We do this in 3 subdirectories: R=real system, H=high level model, L=low level model
!* There is a control file which handles all stuff: xopt.oniom:
!*
!* $oniom
!*   amber <input for amber>
!*   real <name2>.crd/top   ! crd + top file
!*   low  <name3>.crd/top   ! crd + top file
!*   high <name4>.xyz
!*   EE (if found, do electronic/electrostatic embedding via point charges)
!*   TM  <name>.xyz (default)
!*   ORCA  (use orca)


inquire(file='xopt.control',exist=da)
if(da) then
else
 return
endif

print*,'Reading xopt.control $oniom'
open(unit=11,file='xopt.control',status='old')
do
 read(11,'(a)',end=123) aa

if(index(aa,'$oniom').ne.0) then ! found oniom yay
 do
  read(11,'(a)',end=123) aa
  if(index(aa,'amber').ne.0) then
  call charXsplit(aa,bb,2)
  namber=trim(bb)
  endif
  if(index(aa,'real').ne.0) then 
  call charXsplit(aa,bb,2)
  nreal=trim(bb)
  endif
  if(index(aa,'low').ne.0) then 
  call charXsplit(aa,bb,2)
  nlow=trim(bb)
  endif
  if(index(aa,'high').ne.0) then 
  call charXsplit(aa,bb,2)
  nhigh=trim(bb)
  endif
!  if(index(aa,'size').ne.0) then 
!  call charXsplit(aa,bb,2)
!  call charXsplit(aa,cc,3)
!  ireal=s2i(bb)
!  imod=s2i(cc)
!  endif
  if(index(aa,'orca').ne.0.or.index(aa,'ORCA').ne.0) orca=.true.
  if(index(aa,'ee ').ne.0.or.index(aa,'EE ').ne.0) doEE=.true.

 enddo
 

endif
call error(6,"no $oniomin xopt.control")
enddo
123 continue
close(11)



end subroutine
