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
subroutine getgrad
use parm
use logic
use progs
implicit none
character(80) aa,bb
character(120) aaa

if(orca) then
  aa='orca30 orca.in > '//xjob
  call system(trim(aa))
  call rdorcagrad(nat,grad,energy)
  return
endif

! for testing purposes
!if(qmmm)then
!call getoniom
!stop
!endif

!if(numgrad) then
!print*," './mygrad.sh' reads from energy from 'xopt.energy.tmp' "
!call get_numgrad
!return
!endif

if(gaus) then
call system('cat g.in > gtmp.in')
call system("sed '1,2'd xopt.xyz >> gtmp.in")
call system('echo >> gtmp.in')
aa='run-g09 gtmp.in xopt.job'
call system(trim(aa))
call g09grad(nat,grad,energy)
call system("cp xopt.job xopt.last")
return
endif

if(driver) then
call wtm('coord')
call system('driver > xopt.job')
call drivergrad(nat, grad, energy)
call system("cp xopt.job xopt.last")
return
endif


if(mopac) then
call system("echo 'mxyz'  > INPUT")
aa='gcx < '//trim(xyzfile)//' > tmpout'
call system(aa)
!CALL system("cp ~/.setup_mopac SETUP")
call system("mv GCOORD mopac.dat")
call system("mopac12 mopac.dat 2> /dev/null ")
call system("sed -i s/'-'/' -'/g mopac.aux ")
call system("sed -i s/'='/'= '/g mopac.aux ")
aa='cat mopac.out > '//trim(xjob)
call system(aa)

call mopacgrad(nat,grad,energy)
 if(d3) then
   aa='rm -f gradient'
   call system(trim(aa))
   aa='dftd3 '//trim(xyzfile)//' -grad -zero -func pm6 >> '//xjob
!   aa='dftd3 coord -grad -func '//trim(d3func)//' >> '//xjob
   call system(trim(aa))
   call readD3Grad(nat,grad) ! adds gradient
  endif
 call system("cp xopt.job xopt.last")
  return
endif

!if(openbabel) then
!call babelnumgrad
! call system("cp xopt.job xopt.last")
!return
!endif

if(tm.and..not.tmhuge) then
if(relax) goto 909
call wtm('coord')
 if(tmri) then
   aa='ridft > '//xjob
   bb='rdgrad >> '//xjob
   call system(aa)
   call system(bb)
  if(ppot) then 
   aa='ppot coord -tmgrad >> '//xjob
   call system(aa)
  elseif(d3) then
   aa='dftd3 coord -grad -func '//trim(d3func)//' >> '//xjob
   call system(aa)
  elseif(gcp) then
   aa='gcp coord -grad -l "'//trim(gcplevel)// '" >> '//xjob
   call system(aa)
  endif 
 elseif(tmcc) then
   aa='dscf > '//xjob
   bb='ricc2 >> '//xjob
   call system(aa)
   call system(bb)
 else
   aa='dscf > '//xjob
   bb='grad >> '//xjob
   call system(aa)
   call system(bb)
 if(ppot) then                                              
  aa='ppot coord -tmgrad >> '//xjob
  call system(aa)
 elseif(d3) then
  aa='dftd3 coord -grad -func '//trim(d3func)//' >> '//xjob
  call system(aa)
 elseif(gcp) then
  aa='gcp coord -grad -l "'//trim(gcplevel)// '" >> '//xjob
  call system(aa)
 endif 

 endif
909 continue
 call tmgrad(nat,grad,energy)
 call system("cp xopt.job xopt.last")
 return
endif

if(tmhuge) then
call wtm('coord')
 if(tmri) then
   aa='ridft_huge > '//xjob
   bb='rdgrad_huge >> '//xjob
   call system(aa)
   call system(bb)
 elseif(tmcc) then
   aa='dscf_huge > '//xjob
   bb='ricc2_huge >> '//xjob
   call system(aa)
   call system(bb)
 else
   aa='dscf_huge > '//xjob
   bb='grad_huge >> '//xjob
   call system(aa)
   call system(bb)
 endif

 if(ppot) then
  aa='ppot coord -tmgrad >> '//xjob
  call system(aa)
 elseif(d3) then
  aa='dftd3 coord -grad -func '//trim(d3func)//' >> '//xjob
  call system(aa)
 elseif(gcp) then
  aa='gcp coord -grad -l "'//trim(gcplevel)// '" >> '//xjob
  call system(aa)
 endif

 call tmgrad(nat,grad,energy)
 call system("cp xopt.job xopt.last")
 return
endif


if(amber) then
! prepare
aa='rm -f forcedump.dat force.dat > /dev/null'
call system (trim(aa))
call wamber(nat3,xyz,'amber.rst')


! run
aaa='sander -O -i '//trim(ambin)// &
             ' -o '//trim(xjob)// &
             ' -c amber.rst'// &
!             ' -c '//trim(ambcrd)// &
             ' -p '//trim(ambtop)
!print*,aaa

if(APBS) then
aaa='sander.APBS -O -i '//trim(ambin)// &
             ' -o '//trim(xjob)// &
             ' -c amber.rst'// &
             ' -p '//trim(ambtop)
endif
call system(trim(aaa))
!call wamber(nat3,xyz)
call ambergrad(nat,grad,energy,apbs)
 call system("cp xopt.job xopt.last")
return
endif


stop ' no gradient obtained! '
end subroutine getgrad



!*****************
!* AMBER routines *
!*****************
subroutine ambergrad(xnat,g,e,apbs)
implicit none
integer xnat,i,j
real*8 g(3,xnat),e,s2r,g2(3,xnat)
character(80) a,wx
logical da, da2,apbs

da=.false.
da2=.false.


if(apbs) then
open(44,file='xopt.apbs.grad')
   do i=1,xnat
     read (44,*) g(1,i), g(2,i), g(3,i)
   enddo
g=-g*0.529177260d0/627.509d0 
close(44)

open(45,file='xopt.job')
 do
  read(45,'(a)',end=669) a
  if(index(a,'NSTEP       ENERGY').ne.0) then
  read(45,'(a)',end=669) a
  call charXsplit(a,wx,2)
  e=s2r(wx)/627.509d0
  endif
 enddo
669 continue
close(45)

return
endif

inquire(file='forcedump.dat',exist=da)
if(.not.da) call error(6,'forcedump.dat is missing. Check amber.in settings')
g=0d0
open(44,file='forcedump.dat')
 do
  read(44,'(a)',end=666) a
  if(index(a,'Total Force').ne.0) then
   do i=1,xnat
     read (44,*) g(1,i), g(2,i), g(3,i)
   enddo
  endif
 enddo
666 continue
close(44)

! forces (not gradient!) in kcal/mol / angstrom
g=-g*0.529177260d0/627.509d0 


open(44,file='xopt.job')
 do
  read(44,'(a)',end=667) a
  if(index(a,'Etot').ne.0) then
  call charXsplit(a,wx,3)
  e=s2r(wx)/627.509d0
  endif
 enddo
667 continue
close(44)


! for PBSA
! electrostatic forces = qE forces + DB forces
!inquire(file='force.dat',exist=da2)
if(da2) then
open(44,file='force.dat')
 do
  read(44,'(a)',end=668) a
  if(index(a,'Atomic qE forces').ne.0) then
  g2=0d0
   do i=1,xnat
     read (44,*) g2(1,i), g2(2,i), g2(3,i)
   enddo
   g2=-g2*0.529177260d0/627.509d0
   g=g+g2
  endif
  if(index(a,'Atomic DB forces').ne.0) then
   g2=0d0
   do i=1,xnat
     read (44,*) g2(1,i), g2(2,i), g2(3,i)
   enddo
   g2=-g2*0.529177260d0/627.509d0
   g=g+g2
  endif
 enddo
668 continue
close(44)
endif

end subroutine ambergrad

subroutine wamber(nat3,xyzv,fin)
use constant
implicit none
character(10) num
integer nat3,maxcol,mincol,i,j
real(8) xyzv(nat3)
character(*) fin

write(num,'(I10)') nat3/3
open(55,file=fin)
write(55,'(a)') 'made by xopt'
write(55,'(3x,a)')  adjustl(trim(num))
!do i=1,nat3
  maxcol = 0
  do while (maxcol.lt.nat3)
    mincol = maxcol + 1
    maxcol = min(maxcol+6,nat3)
    write(55,'(6f12.7)')(xyzv(j)*au2ang,j=mincol,maxcol)
  enddo 
!enddo
close(55)

end subroutine


!*****************
!* ORCA routines *
!*****************
subroutine rdorcagrad(nat,g,e)
implicit none
real*8 g(3,nat),e
integer nat,i,j
character(80) a

g=0d0
open(unit=33,file='orca.engrad')
read(33,'(a)')a
read(33,'(a)')a
read(33,'(a)')a
read(33,*) i
read(33,'(a)')a
read(33,'(a)')a
read(33,'(a)')a
read(33,*) e
read(33,'(a)')a
read(33,'(a)')a
read(33,'(a)')a
do j=1,nat
   read(33,*)g(1,j)
   read(33,*)g(2,j)
   read(33,*)g(3,j)
enddo
close(33)
return
end




!***************
! TM routines  *
!***************
subroutine tmgrad(xnat,g,e)
implicit none
integer io,i,j,l,k,xnat
real(8), intent(out) :: g(3,xnat),e
real(8) s2r
character*80 aa,wx
io=111


g=0d0
call system('sed -i s/"=-"/"= -"/ gradient')
open (io,file='gradient')
do
 read (io,'(A)',end=666) aa
 if (index(aa,"cycle")/=0) then
  call charXsplit(aa,wx,7)
!  read(wx,'(F)') e
   e=s2r(wx)
   do i=1,xnat
     read (io,'(A)') aa
   enddo
   do i=1,xnat
     read (io,*) g(1,i), g(2,i), g(3,i)
   enddo
! print*,e
 endif
enddo
666 continue
close(io)
end subroutine tmgrad



!************
! mopac12   *
!************
subroutine mopacgrad(nat,grad,e)
implicit none
!real*8 g(3*nat)
real*8 e,grad(nat*3),s2r
integer nat, nn,j
character*80 a,wx
character*10 b,c,d
INTEGER :: i, ios, unit_num,l,k
logical da
CHARACTER(LEN=200) :: record, filename
unit_num = 45
grad=0d0


inquire(file='OLD_MOPAC',exist=da)
if(da) then
print*,'reading mopac.out | OLD_MOPAC found'
! normal output, sometimes printout in *.aux is bad
 filename ="mopac.out"
 open(unit=unit_num, file=trim(filename))
 do
  read(unit_num,'(a)') a
 if(index(a,'FINAL  POINT  AND  DERIVATIVES').ne.0)then
  read(unit_num,'(a)') a
  read(unit_num,'(a)') a
   do i=1,nat*3
    read(unit_num,'(a)') a 
    call charXsplit(a,wx,7)
    grad(i)=s2r(wx)
   enddo
   exit
 endif
 enddo
! now energy
 filename ="mopac.aux"
 open(unit=unit_num, file=trim(filename))
 do
  read(unit_num,'(a)') a
  if(index(a,'HEAT_OF_FORMATION:KCAL/MOL=').ne.0)then
    call charXsplit(a,wx,2)
    e=s2r(wx)/627.509d0
  exit
  endif
 enddo
 close(unit_num)

else ! more precise output
filename ="mopac.aux"
open(unit=unit_num, file=trim(filename))
do
 read(unit_num,'(a)') a
 if(index(a,'HEAT_OF_FORMATION:KCAL/MOL=').ne.0)then
  call charXsplit(a,wx,2)
  e=s2r(wx)/627.509d0
 endif
 if(index(a,'GRADIENTS:KCAL/MOL/ANGSTROM').ne.0)then
   do k=1,int(nat*3/10)
     read(unit_num,*) (grad(nn+(k-1)*10),nn=1,10)
   enddo
   if(int(nat*3/10).eq.0) then ! weniger als 10
    read(unit_num,*) (grad(i),i=1,mod(nat*3,10))
    exit
   endif
  if(int(nat*3/10).eq.1) k=1  ! über 10, unter 20
  if(mod(nat*3,10).gt.1) then   ! über 10, rest lesen
    read(unit_num,*) (grad(i+(k-1)*10),i=1,mod(nat*3,10))
   exit
  endif
  if(mod(nat*3,10).eq.1) read(unit_num,*) grad(nat*3)  ! 11,21,31,...
    exit
 endif
enddo
  close(unit_num)
endif
   grad=grad/(627.509541d0/0.529177260d0)

end subroutine






subroutine readD3Grad(nat,gold)
implicit none
real*8 g(3,nat),gold(3,nat)
integer nat, i
 open(unit=11, file='dftd3_gradient')
 do i=1,nat
   read(11,*)g(1,i),g(2,i),g(3,i)
 enddo
 close(11)
 do i=1,nat
   gold(1:3,i)=gold(1:3,i)+g(1:3,i)
 enddo
end




!**********************
! obenbabel routines  *
!**********************
subroutine ebabel(e)
!real(8) function ebabel
implicit none
integer i,j,n
real(8) s2r,e
character(80)aa
character(80) atmp,a,wx


aa='xopt.py > babel.tmp.out'
call system(aa)
open(321,file='babel.tmp.out')
read(321,*) e
close(321)

return

aa='obenergy -ff GAFF xopt.xyz > babel.tmp.out'
call system(aa)
open(321,file='babel.tmp.out')
do
read(321,'(a)',end=666) atmp
 if(index(atmp,'TOTAL ENERGY').ne.0) then
   call charXsplit(atmp,wx,4)
  e=s2r(wx)/627.509d0
 exit
 endif
enddo
666 continue
close(321)
end subroutine
!end function




subroutine babelnumgrad
use parm
use logic
use progs
implicit none
real(8) el,er,step
print*,'numerical GAFF gradient',nat
step=0.001d0
!energy=ebabel
call ebabel(energy)
do i=1,nat
 do j=1,3
   xyz(j,i)=xyz(j,i)+step
   call newxyz(nat,iat,xyz) ! write new xopt.xyz
   call ebabel(el)
!   el=ebabel
   xyz(j,i)=xyz(j,i)-2d0*step
   call newxyz(nat,iat,xyz)
   call ebabel(er)
!   er=ebabel
   xyz(j,i)=xyz(j,i)+step
   grad(j,i)=(el+er)/2d0
 enddo
enddo

end subroutine



subroutine GetAmberinput(aa)
implicit none
character(2) aa

open(88,file='amber.in',status='replace')
if(aa=='pb') then
write(88,'(a)')  ' '
write(88,'(a)')  '*****************'
write(88,'(a)')  '&cntrl'
write(88,'(a)')  'nstlim  =0, ntx=1,'
write(88,'(a)')  'ipb=2, ntf=1, ntc=1,'
write(88,'(a)')  'ntb=0, irest=0'
write(88,'(a)')  '/'
write(88,'(a)')  '&debugf'
write(88,'(a)')  'do_debugf=1,dumpfrc=1'
write(88,'(a)')  '/'
write(88,'(a)')  '&pb'
write(88,'(a)')  'radiopt=1,space=0.2,dbfopt=1,arcres=0.0625,sprob=1.6'
write(88,'(a)')  '/'
write(88,'(a)')  '&end'
elseif(aa=='gb') then
write(88,'(a)')  ' '
write(88,'(a)')  '*****************'
write(88,'(a)')  '&cntrl'
write(88,'(a)')  'nstlim  =0, ntx=1, '
write(88,'(a)')  'igb=1, ntf=1, ntc=1,'
write(88,'(a)')  'ntb=0, irest=0'
write(88,'(a)')  '/'
write(88,'(a)')  '&debugf'
write(88,'(a)')  'do_debugf=1,dumpfrc=1'
write(88,'(a)')  '/'
write(88,'(a)')  '&end'

endif

close(88)
end subroutine



!**************
! GAUSSIAN
!*************
subroutine g09grad(xnat,grad,e)
implicit none
integer i,j,k,l,xnat,nn
character(80) a,wx
real(8) grad(3,xnat),e,s2r

open(unit=444,file='xopt.job')
do
read(444,'(a)') a
 if(index(a,'SCF Done: ').ne.0)then
   call charXsplit(a,wx,5)
  e=s2r(wx)
 endif
 if (index(a,"Forces (Hartrees/Bohr)")/=0) then
   read(444,'(a)') a
   read(444,'(a)') a
   do i=1,xnat
    read(444,*) k,l,(grad(j,i),j=1,3)
   enddo
   exit
 endif
enddo
grad=-grad !*0.529177260d0
!print*,grad
!stop
close(444)
end subroutine



subroutine drivergrad(xnat,grad,e)
implicit none
integer i,xnat
real(8) grad(3,xnat),e

open(unit=33,file='.DRIVER')
read(33,*) e
do i=1,xnat
   read(33,*) grad(1:3,i)
enddo
close(33)

end subroutine






subroutine get_numgrad
use parm
use logic
use progs
implicit none
real(8) el,er,step
print*,'numerical gradient',nat
step=0.001d0
call system("./mygrad.sh")
call ene(energy)
do i=1,nat
 do j=1,3
   xyz(j,i)=xyz(j,i)+step
   call newxyz(nat,iat,xyz) ! write new xopt.xyz
   call system("mygrad.sh")
   call ene(er)
!   el=ebabel
   xyz(j,i)=xyz(j,i)-2d0*step
   call newxyz(nat,iat,xyz)
   call system("mygrad.sh")
   call ene(el)
   xyz(j,i)=xyz(j,i)+step
   grad(j,i)=(el+er)/2d0
 enddo
enddo

end subroutine

subroutine ene(e)
implicit none
real*8 e
open(44,file='xopt.energy.tmp')
read(44,*) e
close(44)
end subroutine
