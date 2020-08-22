!*******************************************************
!* this manages various kinds of control and log files *
!*******************************************************
!
! OVERVIEW (incomplete..)
! xopt.grad.log   all gradient of the optimization numbered with iteration and energy
! xopt.control    flags for the program
! xopt.restart
!
!
!


!**********************
!* global config file *
!**********************
subroutine read_gconfig
use fiso, only: r8, stdout
use popt
use logic
implicit none
character(200) home,config,aa,s
logical da,fstr
real(r8) s2r
integer s2i

da=.false.
call get_environment_variable('HOME', home)
config=trim(home)//'.xoptrc'

inquire(file=trim(config),exist=da)
if(.not.da) return

! process
write(stdout,'(2x,a)') ' Found global config: .xoptrc'

open(111,file=trim(config))
do
 read(111,'(a)') aa
 call  charXsplit(aa,s,2)
 if(fstr(aa,'maxiter')) maxiter=s2i(s)
 if(fstr(aa,'maxd'))    maxd=s2r(s)
 if(fstr(aa,'gconv'))   gconv=s2r(s)
 if(fstr(aa,'econv'))   econv=s2r(s)
 if(fstr(aa,'maxgrad')) econv=s2r(s)
enddo
close(111)


end subroutine




!********************
!* log the gradient *
!********************
subroutine loggrad(nat,iter,grad,energy)
use fiso, only: r8
implicit none
integer iter,nat,i
real(r8) grad(3,nat),energy

open(unit=11,file='xopt.grad.log',status='old',position='append')
write(11,'(a,2x,I5,2x,F16.8)') '$iter',iter, energy
do i=1,nat
 write(11,'(3F18.14)') grad(1,i),grad(2,i),grad(3,i)
enddo
close(11)
end subroutine

!*************************************************
!* make initial log files so we can append later *
!*************************************************
subroutine touchlog
implicit none
! gradients
open(unit=11,file='xopt.grad.log',status='replace')
  write(11,'(a)') 'logfile for the xopt gradients'
close(11)

!geometries
open(unit=11,file='xopt.log',status='replace')
!write(11,'(a)') 'logfile for the xopt gradients'
close(11)

end subroutine touchlog




!*************************************
!* restart informations
!* writes binary hessian: xopt.hess.restart
!* restraining informations
!* velocities if MD
!*************************************
subroutine wrestart(iter)
use parm
use logic
use progs
use popt
use MDdat, only: velo
implicit none
integer, intent(in) :: iter
character(120) aa

if(restrain) open(unit=33,file='xopt.restrain.tmp',status='old')

open(unit=11,file='xopt.restart',status='replace')
write(11,'(a)') 'This is the restart file for xopt'
write(11,'(a)') 'you may change values by hand   '
write(11,'(a)')

if(restrain) then
write(11,'(a)') '$restr'
do
 read(33,'(a)',end=333) aa
 write(11,'(a)') aa
enddo
333 close(33)
endif

! write freeze section
write(11,'(a)') '$popt'
write(11,'(a,2x,I5)') 'iter ',iter

if(do_md) then
write(11,'(a)') '$velo'
do i=1,nat
 write(11,'(3(F20.12,1x))') velo(1,i),velo(2,i),velo(3,i)
enddo
endif

close(11)
! write last hessian
if(.not.do_md) call writebin(nat3,chess,'xopt.hess.restart')

end subroutine wrestart



! read xopt.control
! handle logical and basic things
! detailed reading (eg for restraints) is done below
subroutine rcontrol
use fiso, only: stdout
use parm
use logic
use progs
use popt
implicit none
character(80) aa,bb
logical da

inquire(file='xopt.control',exist=da)
if(da) then
  write(stdout,'(a)')''
  write(stdout,'(a)') ' * found control file: xopt.control ! *'
  write(stdout,'(a)')''
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
    if(index(aa,'econv').ne.0)     read(aa,*) bb,econv
    if(index(aa,'gnorm').ne.0)     read(aa,*) bb,gnorm
    if(index(aa,'maxgrad').ne.0)   read(aa,*) bb,maxgrad
    if(index(aa,'max displ').ne.0) read(aa,*) bb,dconv
    if(index(aa,'displ').ne.0)     read(aa,*) bb,maxd
    if(index(aa,'maxiter').ne.0)   read(aa,*) bb,maxiter
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


! if(index(aa,'$oniom').ne.0) then
!   qmmm=.true.
! endif

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
use strings
use fiso, only: r8,stdout
use constant, only:au2ang
implicit none
character(255) aa,bb
logical da,tar
integer ii,jj,kk,ll
integer s2i
real(r8) s2r,dih,dgrad,ang,dang,di360,grad2rad
real(r8) bdist,dbond

tar=.false.
da=.false.

if(restart) then
  inquire(file='xopt.restart',exist=da)
  if(.not.da) then
    write(stdout,'(a)')'INFO: no xopt.restart file found. Just using old hessian. '
    return
  endif
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
!  print*,'If using restraints consider citing : H. Kruse, J. Sponer PCCP, 2015,17, 1399-1410'
!  print*,''
  if(restart) then
    write(stdout,'(a)') ' * reading restraints information: xopt.restart ! *'
  else
    write(stdout,'(a)') ' * reading restraints information: xopt.control ! *'
  endif
  write(stdout,'(a)')''
 restrain=.true.
 irest_atom=0
 irest_vol=0_r8
 irest_bond=0
 irest_ang=0
 irest_dihed=0
 irest_konst=0.0_r8
  do
  tar=.false.
    read(11,'(a)',end=123) aa
    if(index(aa,'#').ne.0) cycle

! volume restraints. only sphere for now
!   vol <type> <params>
    if(index(aa,'vol').ne.0) then
       ires=ires+1
       if(index(aa,'sphere').ne.0) then
         ! vol sphere <radius>/vdw
         if(index(aa,'vdw').ne.0) then ! molecular extension plus vdW radius. Scaled by +10%
             irest_vol(ires,1)=(maxval(cell)/2d0)*1.1
         else
             call str_parse(aa,3,irest_vol(ires,1))
         endif
         call str_parse(aa,4,irest_konst(ires,5))
       write(stdout,'(2x,a,2x,F7.2,a,F10.5)') 'sphere radius [A]', irest_vol(ires,1),'| restrain ',irest_konst(ires,5)
       write(33,'(2x,a,F7.2,2x,F7.2)') 'sphere',irest_vol(ires,1),irest_konst(ires,5)
       endif
    endif

! cart. space restraints
    if(index(aa,'atom').ne.0) then
!    atom <at nr> <k>
     ires=ires+1
     call charXsplit(aa,bb,2)
     irest_atom(ires)=s2i(bb)
     ii=irest_atom(ires)
     call charXsplit(aa,bb,3)
     irest_konst(ires,4)=s2r(bb)
     if(restart) then
         call charXsplit(aa,bb,4)
         bdist=s2r(bb)
     else
         bdist=0.0_r8
     endif
     val0(ires)=bdist
     write(stdout,'(2x,a,I4,2x,a,F10.5,a,F7.2)') 'atom [A]', ii,' | restrain',irest_konst(ires,4),' | value: ',bdist*au2ang
     write(33,'(2x,a,I4,2x,F9.5)') 'atom',ii,irest_konst(ires,4),bdist
    endif

!   bond restraints
    if(index(aa,'bond').ne.0) then
     ires=ires+1
     call charXsplit(aa,bb,2)
      irest_bond(ires,1)=s2i(bb)
      ii=irest_bond(ires,1)
     call charXsplit(aa,bb,3)
      irest_bond(ires,2)=s2i(bb)
      jj=irest_bond(ires,2)
     call charXsplit(aa,bb,4)
      irest_konst(ires,1)=s2r(bb)

     if(restart) then
         call charXsplit(aa,bb,5)
         bdist=s2r(bb)
     else
         bdist=dbond(xyz0(1,ii),xyz0(1,jj))
     endif

     if(index(aa,'target').ne.0) then
       tar=.true.
       call charXsplit(aa,bb,6)
       bdist=s2r(bb)/au2ang
     endif

     val0(ires)=bdist

!     print*, '** not well tested!! **'
     if(tar) then
       write(stdout,'(2x,a,2I4,2x,a,F10.5,a,F7.2)') 'bond [A]', ii,jj,' | restrain',irest_konst(ires,1),' | target: ',bdist*au2ang
     else
       write(stdout,'(2x,a,2I4,2x,a,F10.5,a,F7.2)') 'bond [A]', ii,jj,' | restrain',irest_konst(ires,1),' | value: ',bdist*au2ang
     endif
     write(33,'(2x,a,2I4,2(2x,F9.5))') 'bond',ii,jj,irest_konst(ires,1),bdist
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

     if(restart) then
         call charXsplit(aa,bb,6)
         ang=s2r(bb)
         dang=di360(ang)
     else
      call angle(xyz0,ii,jj,kk,ang,dang)
     endif

      if(index(aa,'target').ne.0) then
      tar=.true.
       call charXsplit(aa,bb,7)
       dang=s2r(bb)
       ang=grad2rad(dang)
      endif
    val0(ires)=ang
    if(tar) then
      write(stdout,'(2x,a,3I4,2x,a,F9.5,a,F7.2)') 'angle [deg]', ii,jj,kk,' | restrain',irest_konst(ires,2),' | target: ',dang
    else
      write(stdout,'(2x,a,3I4,2x,a,F9.5,a,F7.2)') 'angle [deg]', ii,jj,kk,' | restrain',irest_konst(ires,2),' | value: ',dang
    endif
    write(33,'(2x,a,3I4,2x,F9.5,2x,F12.8)') 'angle', ii,jj,kk,irest_konst(ires,2),ang
    endif
    if(index(aa,'dihed').ne.0.or.index(aa,'torsion').ne.0) then
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
      if(index(aa,'target').ne.0) then
      tar=.true.
       call charXsplit(aa,bb,8)
       dgrad=s2r(bb)
       dih=grad2rad(dgrad)
      endif
     val0(ires)=dih
   if(tar) then
    write(stdout,'(2x,a,4I4,2x,a,F9.5,a,F7.2)') 'torsion [deg]',ii,jj,kk,ll,'| restrain:',irest_konst(ires,3),' | target: ',dgrad
   else
    write(stdout,'(2x,a,4I4,2x,a,F9.5,a,F7.2)') 'torsion [deg]',ii,jj,kk,ll,'| restrain:',irest_konst(ires,3),' | value: ',dgrad
   endif
    write(33,'(2x,a,4I4,2x,F9.5,2x,F12.8)') 'dihed ',ii,jj,kk,ll,irest_konst(ires,3),dih
    endif
    if(index(aa,'$').ne.0) then
     backspace(11)
     exit
    endif
  enddo
 endif

enddo
123 continue
  ! if (irest_ang(1,1)s.gt.0) write(stdout,'(a)') ' angle restrain gradient is done numerically'
124 continue
close(11)
close(33)
write(stdout,'(a)')''

end subroutine



subroutine rcontrolfreeze
use parm
use logic
use progs
use popt
use constant, only: au2ang
use fiso, only: r8,stdout
implicit none
character(80) aa,bb
logical da,fstr
integer ii,itemp(nat),ial,s2i
real(r8) bufsize,s2r,center(3)
!integer, allocatable :: flist(:)
!allocate(flist(nat))

inquire(file='xopt.control',exist=da)
if(da) then
else
 return
endif

itemp=0
bufsize=0.0_r8

open(unit=11,file='xopt.control',status='old')
do
 read(11,'(a)',end=123) aa

 if(index(aa,'$freeze').ne.0) then ! set up ifreeze vector
  write(stdout,'(a)')''
  write(stdout,'(a)') ' * reading frozen atoms information: xopt.control ! *'
  write(stdout,'(a)')''
  freeze=.true.
  ifrez=0
  do
   read(11,'(a)',end=123) aa
   if(index(aa,'#').ne.0) cycle
   !** freeze all non-H atoms **!
   if(index(aa,'Hopt').ne.0) then
    write(stdout,'(a)') ' --> H-only optimization. Will set all non-H gradients to zero'
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
   if(index(aa,'list').ne.0) then ! freeze a integer list of atoms, range "-" allowed!
    call charXsplit(aa,bb,2)
    call atlist(trim(bb),itemp,ial)
    do i=1,nat
     do ii=1,ial
      if(i==itemp(ii)) ifrez(i)=1
     enddo
    enddo
   endif
   if(index(aa,'region').ne.0) then
    do
      read(11,'(a)',end=123) aa
      if(index(aa,'#').ne.0) cycle
      if(index(aa,'sphere').ne.0) then
        call charXsplit(aa,bb,2)
        bufsize=s2r(bb)
      endif
      if(index(aa,'box').ne.0) then
        call charXsplit(aa,bb,2)
        bufsize=s2r(bb)
        call freeze_box(bufsize,ifrez)
      endif
      if(fstr(aa,'atom').or.fstr(aa,'atom')) then
      call charXsplit(aa,bb,2)
       center(1:3)=xyz(1:3,s2i(bb))
      endif
      if(fstr(aa,'space').or.fstr(aa,'space')) then
       call charXsplit(aa,bb,2)
       center(1)=s2r(bb)
       call charXsplit(aa,bb,3)
       center(2)=s2r(bb)
       call charXsplit(aa,bb,4)
       center(3)=s2r(bb)
      endif
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

if(bufsize>0.0_r8) then
  write(stdout,'(a)')'Freezing:'
  write(stdout,'(a,F6.2)') '  buffer [A] ',bufsize
  write(stdout,'(a,3(F8.3,1x))') '  center [A] ',center*au2ang
  call freeze_sphere(bufsize,center)
endif

if(freeze) then
! print freeze info
j=0
write(stdout,'('' freeze info (atom nr. list): '')')
do i=1,nat
 if(ifrez(i)==1) j=j+1
 if(ifrez(i)==1.and.nat<10) write(stdout,'(I2)',advance="no") i
 if(ifrez(i)==1.and.nat<100.and.nat>10) write(stdout,'(I3)',advance="no") i
 if(ifrez(i)==1.and.nat<1000.and.nat>100) write(stdout,'(I4)',advance="no") i
 if(ifrez(i)==1.and.nat<10000.and.nat>1000) write(stdout,'(I5)',advance="no") i

 if(j==10) then
  write(stdout,'(a)') ''
  j=0
 endif
enddo
 write(stdout,'(a)') ''
endif

end subroutine

subroutine freeze_box(buffer,ifrez)
implicit none
real(8) buffer
integer ifrez(*)

call error('not implemented')
end subroutine


subroutine freeze_sphere(buffer,center)
! freeze all atoms outside a sphere with radius buffer
! located at center
! buffer in angstrom
use parm, only: xyz,nat,ifrez
use constant, only: au2ang
implicit none
real(8) :: buffer,center(3)
real(8) :: r
integer :: i

do i=1,nat
  call veclen2(xyz(:,i),center,r)
  r=r*au2ang
  if(r>buffer) ifrez(i)=1
! print*, i,ifrez(i),r,buffer
enddo

end subroutine


subroutine ReadControlOniom
use parm
use logic
use progs
use popt
use mod_qmmm
implicit none
character(80) aa,bb
logical da

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
call error("no $oniomin xopt.control")
enddo
123 continue
close(11)

end subroutine


subroutine MDrestart
use fiso, only: r8,stdout
use MDdat
use parm, only: nat
implicit none
integer        :: io,i
character(255) :: ch
logical        :: fstr

open(newunit=io,file='xopt.restart',status='old')
do 
 read(io,'(a)',end=666) ch
 if(fstr(ch,'$velo')) then
   write(stdout,'(a)')' Reading velocities from xopt.restart'
   do i=1,nat
      read(io,*) velo(1,i),velo(2,i),velo(3,i)
   enddo
   exit
 endif
enddo
close(io)
return 
! failing to find $velo:
666 continue
call error('failed reading velocities from restart file!')
end subroutine


subroutine CheckStop()
! if STOP is found, stop program and remove STOP file
implicit none
logical da

inquire(file='STOP',exist=da)
if(da) then
 open(unit=1234, file='STOP', status='old')
 close(1234, status='delete')
 call error('found STOP file! Removing and stopping.')
endif

end subroutine


