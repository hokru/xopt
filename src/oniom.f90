!* This subroutine handles level-2 ONIOM-ME calculations with sander(amber suite)
!* We do this in cartesian coordinates
!* Option1: E_oniom(QM:MM)=E_MM_real + (E_QM_model-E_MM_real)
!* planned: Option2: E_oniom(QM1:QM2)=E_QM2_real + (E_QM1_model-E_QM2_model)
!* The user has to prepare all necessary input before calling xopt.
!* stage1. check if input files are there
!* stage2. real system , low level -> ere,gre
!* stage3. model system, low level -> el,gl
!* stage4. model system, high level -> eh,gh
!* stage5. finalize -> eon,gon
!*

!* INPUT: (module qmmm)
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
!*   real <name2>.crd <name2>.top  (name2 needs to be identical)
!*   low  <name3>.crd <name3>.top  (name3 needs to be identical)
!*   high <name4>.xyz
!*   EE (if found, do electronic/electrostatic embedding via point charges)
!*   TM  <name>.xyz (default)
!*   ORCA  (use orca)
!*
!*
!*
subroutine getoniom
use parm
use popt
use logic
use mod_qmmm
implicit none
integer is
real(8) ere,el,eh,eon,pc(nat),diff
real(8), allocatable:: gre(:,:),gl(:,:),gh(:,:),gon(:,:)
real(8) xgre(3,1000)
real(8), allocatable::  creal(:,:),cmod(:,:),ctmp(:)
!integer, intent(in) :: isys(ireal)
logical ex1,ex2,ex3,ex
character(80), parameter :: xjob = 'xopt.job  2>&1 '
character(120) aa,bb
character(280) pwd

ex=.false.
ex1=.false.
ex2=.false.
ex3=.false.



!***********************************************************
! stage1, check input files, determined imod,ireal,isys
if(debug) print*,'Entering... stage1(peparation)'
call getcwd(pwd) ! working directory
call ReadControlOniom ! read control file data or oniom

aa='mkdir R H L'
call system (trim(aa))

! read model data
aa=trim(nlow)//'.crd'
open(111,file=trim(aa))
 read(111,'(a)') aa ! name
 read(111,*) imod
 allocate(ctmp(3*imod),cmod(3,imod))
 read(111,*) ctmp
close(111)

k=0
do i=1,imod
 do j=1,3
 k=k+1
 cmod(j,i)=ctmp(k)
 enddo
enddo
deallocate(ctmp)

! read real system data
aa=trim(nreal)//'.crd'
open(111,file=trim(aa))
 read(111,'(a)') aa ! name
 read(111,*) ireal
 allocate(ctmp(3*ireal),creal(3,ireal),isys(ireal))
 read(111,*) ctmp
close(111)

k=0
do i=1,ireal
 do j=1,3
 k=k+1
 creal(j,i)=ctmp(k)
 enddo
enddo
deallocate(ctmp)

allocate(gre(3,ireal),gh(3,imod),gl(3,imod),gon(3,ireal))

isys=0
do i=1,ireal
 do j=1,imod
  diff=creal(1,i)-cmod(1,j)
  if(abs(diff)<=2.5d-3) isys(i)=1
 enddo
! print*,i,isys(i)
enddo


!***********************************************************

if(debug) print*,'Entering... stage2(real system)'

! stage2, calculate MM for real system using Amber
!call system('cd R')
call chdir(trim(pwd)//'/R')

! remove old amber files
aa='rm -f forcedump.dat force.dat > /dev/null'
call system (trim(aa))

! copy inputs
aa='cp ../'//trim(nreal)//'.crd .'
 call system(aa)
aa='cp ../'//trim(nreal)//'.top .'
 call system(aa)
aa='cp ../'//trim(namber)//' .'
 call system(aa)


! calculate
!call wamber(nat3,xyz,trim(nreal)//'.rst')  !  ???
aa='sander -O -i '//trim(namber)// &
             ' -o xopt.real.out -c '//trim(nreal)//'.crd'// &
             ' -p '//trim(nreal)//'.top'
call system(aa)
call ambergrad(ireal,gre,ere,.false.)
!do i=1,ireal
!write(*,*) gre(1:3,i)
!enddo

if(debug) print*,'Entering... stage3(low model)'

! stage3, calculate MM for low level model system using Amber
call chdir(trim(pwd)//'/L')
aa='rm -f forcedump.dat force.dat > /dev/null'
call system (trim(aa))

! copy inputs
aa='cp ../'//trim(nlow)//'.crd .'
 call system(aa)
aa='cp ../'//trim(nlow)//'.top .'
 call system(aa)
aa='cp ../'//trim(namber)//' .'
 call system(aa)

! calculate
!call wamber(nat3,xyz,trim(nreal)//'.rst')  !  ???
aa='sander -O -i '//trim(namber)// &
             ' -o xopt.low.out -c '//trim(nlow)//'.crd'// &
             ' -p '//trim(nlow)//'.top'
call system(aa)
call ambergrad(imod,gl,el,.false.)
!do i=1,imod
!write(*,*) gl(1:3,i)
!enddo


if(debug) print*,'Entering... stage4(high model)'

! stage4, QM calulation of high level model system
call chdir(trim(pwd)//'/H')
aa='cp ../'//trim(nhigh)//' .'
 call system(aa)

if(ORCA) then
stop 'Sorry, not done yet'

endif

if(TM) then
  print*,'test implementation, PBE/SVP'
  call wtm('coord')
!  call system('cefine -bas SVP -func pbe -noopt -chrg -1')
  if(doEE) then
  call system('kdg end')
!  call readamber(ireal,iat,xyz,pc,trim(nlow)//'.crd',trim(nlow)//'.top') ??
  open(unit=22,file='control',status='old',position='append')
   write(22,'(a)') '$point_charges'
   do i=1,nat
     write(22,'(3(F18.14,1x),F16.8)') xyz(1:3,i), pc(i)
   enddo
  close(22)
  call system('echo "$end" >> control')
  endif
  call wtm('coord')
  aa='ridft > '//xjob
  bb='rdgrad >> '//xjob
!  call system(aa)
!  call system(bb)
  call tmgrad(imod,gh,eh)
  call system("cp xopt.job xopt.last")
endif

call chdir(trim(pwd))

!**********************
!* CALCULATIONS DONE. *
!**********************


! glue energie together
eon=ere+(eh-el)

! glue gradient together
gon=0d0
do i=1,ireal
 if(isys(i)==0) then
   gon(1:3,i)=gre(1:3,i)
 else
   gon(1:3,i)=gre(1:3,i)+(gh(1:3,i)-gl(1:3,i))
 endif
enddo

if(debug) then
print*,'Energy ',eon, eon*627.51
write(*,*) 'Grad '
do i=1,ireal
write(*,*) gon(1:3,i)
enddo
endif

end subroutine



! not used
!subroutine onamber()
!implicit none
!integer info
!character(80) aa
!character(*) inp,outp,coord,top
!
!aa='sander -O -i '//trim(inp)// &
!             ' -o '//trim(outp)// &
!             ' -c '//trim(coord) &
!             ' -p '//trim(top)
!call system(aa,info)
!if(info.ne.0) then
! print*,info
! call EXECUTE_COMMAND_LINE(aa,exitstat=info)
! call error(' amber shell command failed: <onamber>')
!subroutine
