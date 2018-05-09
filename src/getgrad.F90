subroutine getgrad
! central gradient routine
use parm
use logic
use progs
use constant, only: au2ang
use fiso, only: r8
implicit none
character(255) aa,bb
character(255) aaa
real(r8) dummy,gdummy(3,nat)

!************
!* wrappers *
!************
if(ciopt) then
call wtm('coord')
if(.not.numgrad) then
  call wtm('stateI.xopt/coord')
  call wtm('stateJ.xopt/coord')
  call wrxyz('stateI.xopt/xopt.xyz')
  call wrxyz('stateJ.xopt/xopt.xyz')
endif
call cioptgrad
return

endif


! for testing purposes
if(qmmm)then
call error('nope')
!call getoniom
!call orcaoniom
stop
endif

!
! if(numgrad_internal) then
!
!
! endif

!*********************
!* purely standalone *
!*********************

! G E I
if(gei) then
  call system(trim(command_gei))
  call egradfile(nat,grad,energy,trim(geigrad))
  call getd3gcp(nat,energy,grad)
  return
endif

!* O R C A
if(orca) then
  call system(trim(command_orca))
  call rdorcagrad(nat,grad,energy,'orca.engrad')
  call getd3gcp(nat,energy,grad)
  call cpfile('xopt.job','xopt.last','.')
  return
endif

if(xtb) then
  call system(trim(command_xtb))
  call tmgrad(nat,grad,energy)  ! xtb writes TM gradient/energy file
  call cpfile('xopt.job','xopt.last','.')
  return
endif

!*********************
!* numerical grad    *
!* GEI interface     *
!*********************
if(numgrad) then
  if(nproc>1) then
   call get_pgrad(energy,dummy,grad,gdummy)
   return
  else
   call get_numgrad(xyz,energy,grad)
  endif
  call getd3gcp(nat,energy,grad)
  return
endif

! N W C H E M
if(nwchem) then
  call IOnwchem('nw.in')
  call system(trim(command_nwchem))
  call nwchemGrad(nat,grad,energy,'xopt.job')
  call getd3gcp(nat,energy,grad)
call cpfile('xopt.job','xopt.last','.')
endif

! G A U S S I A N   0 9
if(gaus) then
  call IOgaussian('g.in')
  call system(trim(command_gaus))
  call g09grad(nat,grad,energy)
  call getd3gcp(nat,energy,grad)
  call cpfile('xopt.job','xopt.last','.')
return
endif

! D R I V E R (Grimme)
if(driver) then
call wtm('coord')
call system('driver > xopt.job')
call drivergrad(nat, grad, energy)
call getd3gcp(nat,energy,grad)
call cpfile('xopt.job','xopt.last','.')
return
endif


! M O P A C (modern)
if(mopac) then
 call IOmopac12('mopac.dat')
 call system(trim(command_mopac))
 ! make sure large numbers are separated from "=" and "-" symbols.
 call system("sed -i s/'-'/' -'/g mopac.aux ")
 call system("sed -i s/'='/'= '/g mopac.aux ")
 call cpfile('mopac.out','xopt.job','.')
 call mopacgrad(nat,grad,energy)
 call getd3gcp(nat,energy,grad)
!  if(d3) then
!    aa='rm -f gradient'
!    call system(trim(aa))
!    aa='dftd3 '//trim(xyzfile)//' -grad -func '//trim(d3func)//' >> '//xjob
!    call system(trim(aa))
!    call readD3Grad(nat,grad,energy,'.EDISP','dftd3_gradient') ! adds gradient
!   endif
call cpfile('xopt.job','xopt.last','.')
return
endif

! O P E N B A B E L
if(openbabel) then
call babelnumgrad
 call cpfile('xopt.job','xopt.last','.')
return
endif

! T U R B O M O L E
if(tm) then
  call wtm('coord')
  aa='dscf > '//xjob
  bb='grad >> '//xjob
  if(tmhuge) aa='dscf_huge > '//xjob
  if(tmhuge) aa='grad_huge >> '//xjob
  if(tmri) then
    aa='ridft > '//xjob
    bb='rdgrad >> '//xjob
    if(tmhuge) bb='ridft_huge > '//xjob
    if(tmhuge) bb='rdgrad_huge >> '//xjob
  endif
  if(tmcc) then
    bb='ricc2 >> '//xjob
    if(tmhuge) bb='ricc2_huge >> '//xjob
  endif
  if(riper) then
   aa='riper > '//xjob
   if(tmhuge) aa='riper_huge > '//xjob
   bb='' 
  endif
  if(debug) then
    print*,'TM calls:'
    print*,trim(aa)
    print*,trim(bb)
  endif
  call system(trim(aa))
  call system(trim(bb))
  call tmgrad(nat,grad,energy)
  call getd3gcp(nat,energy,grad)
  call cpfile('xopt.job','xopt.last','.')
return
endif

! A M B E R (sander)
if(amber) then
  ! prepare
  aa='rm -f forcedump.dat force.dat > /dev/null'
  call system (trim(aa))
  call wamber(nat3,xyz,'amber.rst')
  if(APBS) then
  aaa='sander.APBS -O -i '//trim(ambin)// &
               ' -c amber.rst'// &
               ' -p '//trim(ambtop)// &
               ' -o '//trim(xjob)
  endif
  call system(trim(command_amber))
  !call wamber(nat3,xyz)
  call ambergrad(nat,grad,energy,apbs)
  call cpfile('xopt.job','xopt.last','.')
  return
endif

!  P S I 4
if(psi4) then
  call system(trim(command_psi4))
  call psi4grad(nat,grad,energy,'xopt.job')
  call cpfile('xopt.job','xopt.last','.')
  return
endif

! G A M E S S
if(gamess) then
  call IOgamess(trim(gmsin)//'.inp')
  call system(trim(command_gms))
  call gmsgrad(nat,grad,energy,'xopt.job')
  call getd3gcp(nat,energy,grad)
  call cpfile('xopt.job','xopt.last','.')
  return
endif

call error(' no gradient obtained! (ups?!) ')
end subroutine getgrad

!*******************
!* GAMESS routines *
!*******************
subroutine gmsgrad(n,g,e,filen)
use fiso, only: r8
implicit none
integer n,io,i,ii
real(r8) g(3,n),e,s2r
character(200) aa
character(*) filen
logical fstr
character(80) wx
e=9
open(newunit=io,file=filen)
 do
  read(io,'(a)',end=666) aa
  if(fstr(aa,'FINAL RDFTB ENERGY')) then
   call charXsplit(aa,wx,5)
   e=s2r(wx)
  endif
  if(fstr(aa,'GRADIENT OF THE ENERGY')) then
   read(io,'(a)') aa
   read(io,'(a)') aa
   read(io,'(a)') aa
   do i=1,n
     read(io,*) ii,aa,g(1:3,i)
   enddo
  endif
 enddo
666 continue
if(e>0) call error(' no energy in GAMESS output found')
close(io)

end subroutine


!*****************
!* PSI4 routines *
!*****************
subroutine psi4grad(n,g,e,filen)
use fiso, only: r8
implicit none
integer n,io,i
real(r8) g(3,n),e,x,s2r
character(200) aa
character(*) filen
logical fstr
character(80) wx

open(newunit=io,file=filen)
  do
   read(io,'(a)',end=666) aa
   call lower_case(aa)
   if(fstr(aa,'total energy =')) then
     call charXsplit(aa,wx,4)
     e=s2r(wx)
   endif
   if(fstr(aa,'total gradient:')) then
     read(io,'(a)') aa
     read(io,'(a)') aa
     do i=1,n
       read(io,*) x, g(1:3,i)
     enddo
   exit
   endif
  enddo
  666 continue
close(io)
  return
end subroutine

!*****************
!* AMBER routines *
!*****************
subroutine ambergrad(xnat,g,e,apbs)
use fiso, only: r8, stdout
use constant, only: au2ang, au2kcal
implicit none
integer xnat,i
real(r8) g(3,xnat),e,s2r,g2(3,xnat)
character(80) a,wx
logical da, da2,apbs

da=.false.
da2=.false.

inquire(file='amber.force',exist=da)
if(da) then
!write(stdout,'(a)') 'using amber.force'
! for the modifed amber
open(44,file='amber.force')
   read(44,*) e
   do i=1,xnat
     read (44,*) g(1,i), g(2,i), g(3,i)
   enddo
! forces (not gradient!) in kcal/mol / angstrom
!g=-g*0.529177260d0/627.509d0
!e=e/627.509d0
g=-g
close(44)
return
endif
! ignore rest for now.

!
! if(apbs) then
! open(44,file='xopt.apbs.grad')
!    do i=1,xnat
!      read (44,*) g(1,i), g(2,i), g(3,i)
!    enddo
! g=-g*0.529177260d0/627.509d0
! close(44)
!
!
! open(45,file='xopt.job')
!  do
!   read(45,'(a)',end=669) a
!   if(index(a,'NSTEP       ENERGY').ne.0) then
!   read(45,'(a)',end=669) a
!   call charXsplit(a,wx,2)
!   e=s2r(wx)/627.509d0
!   endif
!  enddo
! 669 continue
! close(45)
! return
! endif

! kept for legacy amber
inquire(file='forcedump.dat',exist=da)
if(.not.da) call error('forcedump.dat is missing. Check amber.in settings')
g=0_r8
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
! g=-g*0.529177260d0/627.509d0
g=-g*au2ang/au2kcal

open(44,file='xopt.job')
 do
  read(44,'(a)',end=667) a
  if(index(a,'Etot').ne.0) then
  call charXsplit(a,wx,3)
  e=s2r(wx)/au2kcal
  endif
 enddo
667 continue
close(44)


! for PBSA
! electrostatic forces = qE forces + DB forces
inquire(file='force.dat',exist=da2)
if(da2) then
open(44,file='force.dat')
 do
  read(44,'(a)',end=668) a
  if(index(a,'Atomic qE forces').ne.0) then
  g2=0_r8
   do i=1,xnat
     read (44,*) g2(1,i), g2(2,i), g2(3,i)
   enddo
   g2=-g2*au2ang/au2kcal !0.529177260d0/627.509d0
   g=g+g2
  endif
  if(index(a,'Atomic DB forces').ne.0) then
   g2=0_r8
   do i=1,xnat
     read (44,*) g2(1,i), g2(2,i), g2(3,i)
   enddo
   g2=-g2*au2ang/au2kcal! 0.529177260d0/627.509d0
   g=g+g2
  endif
 enddo
668 continue
close(44)
endif

end subroutine ambergrad

subroutine wamber(nat3,xyzv,fin)
use constant, only: au2ang
use fiso, only: r8
implicit none
character(10) num
integer nat3,maxcol,mincol,j
real(r8) xyzv(nat3)
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
subroutine rdorcagrad(nat,g,e,fname)
use fiso, only: r8
implicit none
integer nat,i,j
real(r8) g(3,nat),e
character(80) a
character(*) fname

g=0_r8
open(unit=33,file=fname)
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
use fiso, only: r8
implicit none
integer io,i,xnat
real(r8), intent(out) :: g(3,xnat),e
real(r8) s2r
character(80) aa,wx
io=111


g=0_r8
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
use constant, only: au2ang, au2kcal
use fiso, only: r8
implicit none
integer nat, nn
real(r8) e,grad(nat*3),s2r
character(80) a,wx
character(10) c,d
INTEGER :: i, unit_num,k
logical da
CHARACTER(LEN=200) :: filename

grad=0_r8
inquire(file='OLD_MOPAC',exist=da)
if(da) then
  print*,'reading mopac.out | OLD_MOPAC found'
  ! normal output, sometimes printout in *.aux is bad
   filename ="mopac.out"
   open(newunit=unit_num, file=trim(filename))
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
      e=s2r(wx)/au2kcal
    exit
    endif
   enddo
   close(unit_num)
else ! more precise output (use this!)
  filename ="mopac.aux"
  open(unit=unit_num, file=trim(filename))
  do
   read(unit_num,'(a)') a
   if(index(a,'HEAT_OF_FORMATION:KCAL/MOL=').ne.0)then
    call charXsplit(a,wx,2)
    e=s2r(wx)/au2kcal
   endif
   if(index(a,'GRADIENTS:KCAL/MOL/ANGSTROM').ne.0)then
     do k=1,int(nat*3/10)
       read(unit_num,*) (grad(nn+(k-1)*10),nn=1,10)
     enddo
     if(int(nat*3/10).eq.0) then ! weniger als 10
      read(unit_num,*) (grad(i),i=1,mod(nat*3,10))
      exit
     endif
    if(int(nat*3/10).eq.1) k=1  ! ueber 10, unter 20
    if(mod(nat*3,10).gt.1) then   ! ueber 10, rest lesen
      read(unit_num,*) (grad(i+(k-1)*10),i=1,mod(nat*3,10))
     exit
    endif
    if(mod(nat*3,10).eq.1) read(unit_num,*) grad(nat*3)  ! 11,21,31,...
      exit
   endif
  enddo
    close(unit_num)
endif
   grad=grad/(au2kcal/au2ang)
  !  grad=grad/(627.509541d0/0.529177260d0)
end subroutine




! read a simple energy/gradient file as written by gcp and dftd3
! updates energy and gradient
subroutine readD3Grad(nat,gold,eold,fileE,fileG)
use fiso, only: r8
implicit none
integer nat, i
real(r8) g(3,nat),gold(3,nat),e,eold
character(*) fileE,fileG

open(unit=12,file=fileE)
 read(12,*) e
close(12)
eold=eold+e

 open(unit=11, file=fileG)
 do i=1,nat
   read(11,*)g(1,i),g(2,i),g(3,i)
 enddo
 close(11)
 do i=1,nat
   gold(1:3,i)=gold(1:3,i)+g(1:3,i)
 enddo
end subroutine


!**********************
! obenbabel routines  *
!**********************
subroutine ebabel(e)
! not useful
use fiso, only: r8
use constant, only: au2kcal
implicit none
real(r8)       :: s2r,e
character(80) :: aa
character(80) :: atmp,wx

aa='obenergy -ff GAFF xopt.xyz > babel.tmp.out'
call system(aa)
open(321,file='babel.tmp.out')
do
read(321,'(a)',end=666) atmp
 if(index(atmp,'TOTAL ENERGY').ne.0) then
   call charXsplit(atmp,wx,4)
  e=s2r(wx)/au2kcal
 exit
 endif
enddo
666 continue
close(321)
end subroutine



subroutine babelnumgrad
use parm
use logic
use progs
use fiso, only: r8,stdout
implicit none
real(r8) el,er,step
write(stdout,*) 'numerical GAFF gradient',nat
step=0.001_r8

call ebabel(energy)
do i=1,nat
 do j=1,3
   xyz(j,i)=xyz(j,i)+step
   call newxyz(nat,iat,xyz)
   call ebabel(el)
   xyz(j,i)=xyz(j,i)-2_r8*step
   call newxyz(nat,iat,xyz)
   call ebabel(er)
   xyz(j,i)=xyz(j,i)+step
   grad(j,i)=(el-er)/(2_r8*step)
 enddo
enddo

end subroutine

!***************************
!* pre-defined Amber input *
!***************************
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
!write(88,'(a)')  '&debugf'
!write(88,'(a)')  'do_debugf=1,dumpfrc=1'
!write(88,'(a)')  '/'
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
!write(88,'(a)')  '&debugf'
!write(88,'(a)')  'do_debugf=1,dumpfrc=1'
!write(88,'(a)')  '/'
write(88,'(a)')  '&end'
endif

close(88)
end subroutine



!************************
! GAUSSIAN 09           *
!************************
subroutine g09grad(xnat,grad,e)
use fiso, only: r8
implicit none
integer i,j,k,l,xnat
character(80) a,wx
real(r8) grad(3,xnat),e,s2r

open(unit=444,file='xopt.job')
do
read(444,'(a)') a
 if(index(a,'SCF Done: ').ne.0)then
  call charXsplit(a,wx,5)
  e=s2r(wx)
 endif
 if(index(a,' E=  ').ne.0) then ! for CASSCF
!~   call charXsplit(a,wx,5)
  wx=a(22:40)
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
close(444)
end subroutine


!***********
!* DRIVER  *
!***********
subroutine drivergrad(xnat,grad,e)
use fiso, only: r8
implicit none
integer i,xnat
real(r8) grad(3,xnat),e

open(unit=33,file='.DRIVER')
read(33,*) e
do i=1,xnat
   read(33,*) grad(1:3,i)
enddo
close(33)
end subroutine

! general gradient file with energy
! energy
! gx1 gy1 gz1
! gx2 gy2  gz2
subroutine egradfile(xnat,grad,e,fname)
use fiso, only: r8
implicit none
integer i,xnat,io
character(*) fname
real(r8) grad(3,xnat),e

open(newunit=io,file=fname)
read(io,*) e
do i=1,xnat
   read(io,*) grad(1:3,i)
enddo
close(io)
end subroutine


!**********************
!* numerical gradient *
!**********************

! INTERNAL NUMGRAD
! not used
! maybe switch on/off flags here?
subroutine get_internal_numgrad(e,g)
use fiso, only: r8
use parm, only: i,j,nat,iat,xyz,energy
use logic, only: numgrad
implicit none
real(r8) xyz0(3,nat)
real(r8) el,er,step,e,g(3,nat)

step=0.005_r8
xyz0=xyz
numgrad=.false.
do i=1,nat
write(*,'(a,I2,a,I2,a)') 'gradient of atom [',i, ']/[', nat,']'
 do j=1,3
   xyz(j,i)=xyz(j,i)+step
   call newxyz(nat,iat,xyz)
   call getgrad
   er=energy
   xyz(j,i)=xyz(j,i)-2_r8*step
   call newxyz(nat,iat,xyz)
   call getgrad
   el=energy
   xyz(j,i)=xyz(j,i)+step
   g(j,i)=(er-el)/(step*2_r8)
 enddo
enddo

numgrad=.true.
xyz=xyz0

end subroutine





! EXTERNAL NUMGRAD
subroutine get_numgrad(xyz,e,g)
use parm, only: i,j,k,l,nat,iat
use fiso, only: r8, stdout
!use logic
use progs, only: command_gei
implicit none
real(r8), intent(in) :: xyz(3,nat)
real(r8) nxyz(3,nat)
real(r8) el,er,step,e,g(3,nat)
nxyz=xyz
step=0.005_r8
!step=0.010

! inquire(file='mygrad.sh',exist=da)
! if(.not.da) stop 'mygrad.sh script missing!'
call system(trim(command_gei))
! call system("./mygrad.sh")
call ene(e)
! if(nat>99) stop 'system to large! Comment me out if you REALLY want to do this'
do i=1,nat
write(stdout,'(a,I2,a,I2,a)') 'gradient of atom [',i, ']/[', nat,']'
 do j=1,3
   nxyz(j,i)=nxyz(j,i)+step
   call newxyz(nat,iat,nxyz)
   call system(trim(command_gei))
   call ene(er)
   nxyz(j,i)=nxyz(j,i)-2_r8*step
   call newxyz(nat,iat,nxyz)
   call system(trim(command_gei))
   call ene(el)
   nxyz(j,i)=nxyz(j,i)+step
   g(j,i)=(er-el)/(step*2_r8)
!   print*, el,er,g(j,i)
 enddo
enddo
end subroutine

subroutine ene(e)
use fiso, only: r8
implicit none
real(r8) e
  e=9e9
  open(44,file='xopt.energy.tmp')
  read(44,*,end=666) e
  666 continue
  if(e>=8e9) call error('No energy')
  close(44)
end subroutine


! form two grads at the same time, from the same output
! usefull for CIopt calculations where you can get 2 states at once
subroutine get_numgrad2(e1,e2,g1,g2)
use parm, only: xyz,i,j,k,l,nat,iat
use fiso, only: r8
!use logic
use progs, only: command_gei
implicit none
real(r8) el,er,el2,er2,step,e1,e2,g1(3,nat),g2(3,nat)
logical da
step=0.005_r8
inquire(file=trim(command_gei),exist=da)
if(.not.da) stop 'GEI script missing!'
call system(trim(command_gei))
call ene2(e1,e2)

do i=1,nat
write(*,'(a,I2,a,I2,a)') 'gradient of atom [',i, ']/[', nat,']'
 do j=1,3
   xyz(j,i)=xyz(j,i)+step
   call newxyz(nat,iat,xyz) ! write new xopt.xyz
   call system(trim(command_gei))
   call ene2(er,er2)
   xyz(j,i)=xyz(j,i)-2_r8*step
   call newxyz(nat,iat,xyz)
   call system(trim(command_gei))
   call ene2(el,el2)
   xyz(j,i)=xyz(j,i)+step
   g1(j,i)=0.5_r8*(er-el)/step
   g2(j,i)=0.5_r8*(er2-el2)/step
 enddo
enddo
end subroutine

subroutine ene2(e,e2)
use fiso, only: r8
implicit none
real(r8) e,e2
e=9e9
open(44,file='xopt.energy.tmp')
read(44,*,end=666) e
read(44,*,end=666) e2
666 continue
if(e>=8e9) call error('error. No energy!! <ene2>')
if(e2>=8e9) e2=0
close(44)
end subroutine


! parallel gradient from xopt.pgrad helper
! always returns 2 energies and 2 gradients
! in case we do a CIopt from a single calculation and
! eg read two states at once.
! only e1,g1 are non-zero otherwise.
subroutine get_pgrad(e1,e2,g1,g2)
use fiso, only: r8
use parm,  only: i,nat,iat,xyz
use logic, only: nproc,scratchjob
use progs, only: binpath,command_gei,command_mpi
use constant, only: au2ang
use progs, only: scrdir, workdir,usrscr
implicit none
real(r8) g1(3,nat),g2(3,nat),e1,e2

if(scratchjob) then
 call copy2scr()
 call chdir(trim(scrdir))
 call system(trim(command_gei))
 call ene2(e1,e2)
 call chdir(trim(workdir))
else
 call system(trim(command_gei))
 call ene2(e1,e2)
endif

 !write info
 open(222,file='xopt.para.tmp')
 write(222,'(I5)') nat
 do i=1,nat
  write(222,'(1x,3E22.13,1x,I2)')xyz(1:3,i)*au2ang,iat(i)
 enddo
 write(222,'(a)') trim(usrscr)
 write(222,'(a)') trim(command_gei)
 ! if(.not.usrscr=="0") write(222,'(a)') trim(usrscr)
 close(222)
 !call xopt.pgrad
 ! write(bb,*) nproc
 ! write(command_mpi,'(a)') trim(scall_mpi)//' -output-filename xopt.slave -n '//trim(adjustl(bb))//' '//trim(binpath)
 call system(trim(command_mpi))
 call cpfile('xopt.slave.1.0','xopt.pgrad.out','.')
 ! read gradient 1
 open(unit=33,file='xopt.pgrad.tmp')
 do i=1,nat
    read(33,*) g1(1:3,i)
 enddo
 ! read gradient 2 (might be zero)
 do i=1,nat
    read(33,*) g2(1:3,i)
 enddo
 close(33)

return
end subroutine



subroutine getd3gcp(nat,energy,grad)
! D3 and gCP
! using the library version is preferred
! as the Grimme codes might change behaviour every now and then..
use fiso, only: r8
use progs, only: xjob
use logic, only: d3func,gcplevel,debug,d3,gcp,exlib
implicit none
integer :: nat
real(r8) :: grad(3,nat),energy
character(255) :: aa

if(exlib) then
  call getd3gcpLIB(energy,grad)
  return
else
  if(d3) then
  if(debug) print*,' adding D3 from system call'
!    aa='dftd3 xopt.xyz -grad -func '//trim(d3func)//' >> '//xjob
    aa='dftd3 xopt.xyz -grad -func '//trim(d3func)//' > xopt.dftd3.out'
    call system(aa)
    call readD3Grad(nat,grad,energy,'.EDISP','dftd3_gradient') ! adds gradient
  endif
  if(gcp) then
  if(debug) print*,' adding gCP from system call'
    aa='gcp xopt.xyz -grad -l "'//trim(gcplevel)// '" >> '//xjob
    call system(aa)
    call readD3Grad(nat,grad,energy,'.CPC','gcp_gradient') ! adds gradient
  endif
endif
end subroutine

! gCP-D3 library calls
subroutine getd3gcpLIB(energy,grad)
use fiso, only: r8
use progs, only: xjob
use parm, only: nat,iat,xyz
use logic, only: d3func,gcplevel,debug,d3,gcp
implicit none
integer :: narg,i,io
real(r8) :: grad(3,nat),energy
real(r8) :: e_gcp, e_disp, g_disp(3,nat),g_gcp(3,nat)
character(255) :: arg(10)

e_gcp=0_r8
e_disp=0_r8
arg=''
open(newunit=io,file='xopt.gcpd3.out')
if(d3) then
#ifdef DFTD3
      call cstring(trim(d3func),narg)
      if(narg>9) call error('too many d3 arguments(max 10)')
      do i=1,narg
          call charXsplit(d3func,arg(i),i)
          ! print*,arg(i)
      enddo
      arg(narg+1)='-grad'
      call call_dftd3(nat,iat,xyz,g_disp,e_disp,.true.,arg,io)
      energy=energy+e_disp
      grad=grad+g_disp
#else
  call error('dftd3 library not compiled!')
#endif
endif
if(gcp) then
#ifdef GCP
    call gcp_call &
    (nat,xyz,iat,e_gcp,g_gcp,.true.,.false.,trim(gcplevel),.true.,.false.,io)
      energy=energy+e_gcp
      grad=grad+g_gcp
#else
    call error('gcp library not compiled!')
#endif
endif
close(io)
if(debug) then
print*, 'D3 energy',e_disp
print*, 'gCP energy',e_gcp
endif

!open(newunit=io,file='xopt.job',position='append')
!close(io)
end subroutine



subroutine nwchemGrad(xnat,grad,e,filen)
! should work for SCF/DFT
! Sadly, there does not seem to be a template output format for all nwchem
! modules.
use strings
use fiso, only: r8
integer xnat,io,i,ii
character(*) filen
character(255) aa
character(2) el
real(r8) e,grad(3,xnat),dummy(3) 
logical fstr
grad=0.0
e=0.0
open(newunit=io,file=filen,status='old')
do
 read(io,'(a)',end=666) aa
  if(fstr(aa,'ENERGY GRADIENTS')) then
  read(io,'(a)') aa
  read(io,'(a)') aa
  read(io,'(a)') aa
  do i=1,xnat
    read(io,'(a)') ii,el,dummy(1:3),grad(1:3,i)
  enddo
  endif
  if(fstr(aa,'Total DFT energy')) then
    call str_parse(aa,5,e)
  endif
  if(fstr(aa,'Total SCF energy')) then
    call str_parse(aa,5,e)
  endif
enddo
666 continue
close(io)

end subroutine
