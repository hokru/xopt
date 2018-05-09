subroutine composition
use parm, only: nat,iat
use fiso, only: stdout
implicit none
integer comp(107),c,k,i,ic(107)
character(2) esym
character(5) num
character(255) string

comp=0
ic=0
do k=1,107
c=0
 do i=1,nat
  if(iat(i)==k) then
    c=c+1
    comp(k)=c
    ic(k)=i
  endif
 enddo
enddo

write(stdout,'(3x,a)',advance="no") ' element composition : '
do k=1,107
 if(ic(k)/=0) then
  write(num,'(I5)') comp(k)
  string=trim(adjustl(esym(k)))//trim(adjustl(num))//' '
  call upper_case(string)
  write(stdout,'(a)',advance="no") trim(string)
 endif
enddo
write(stdout,'(a)')''
end subroutine

subroutine getCOM()
use fiso, only: r8,stdout
use atomdata, only: ams
use parm, only: nat,xyz,iat
use constant, only: au2ang
use logic, only: orient
implicit none
real(r8) :: com(3)
real(r8) :: mmass,coma(3)
integer i,j

!xyz=xyz*au2ang
mmass=0_r8
com=0_r8
do i=1,nat
 mmass=mmass+ams(iat(i))
 do j=1,3
  com(j)=com(j)+xyz(j,i)*ams(iat(i))
 enddo
enddo

com=com/mmass
coma=com*au2ang

write(stdout,'(3x,a,F12.5)')  ' molecular mass      : ',mmass
write(stdout,'(3x,a,3(F12.5,1x))') ' center of mass (A)  : ',coma(1:3)
!write(stdout,'(3x,a,3F12.5)') 'center of mass (au)',com(1:3)

if(orient) then
! move molecule to COM
do i=1,nat
 do j=1,3
    xyz(j,i)=xyz(j,i)-com(j)
 enddo
enddo
endif
!call wrxyz('cema_xopt.xyz')
end


subroutine getIntertia
! calculate moment of inertia and principle axis
! z**2+y**2    -xy        -xz
! -xy         z**2+x**2   -yz
! -xz           -yz        x*2+y**2
use parm, only: nat,iat,xyz,i,j
use fiso, only: r8,stdout
use atomdata, only: ams
use constant, only: Planck, pc_c, pi, bohr2m, amu2kg
use logic, only: orient
implicit none
real(r8) mom(3,3),e(3),rot(3)
real(r8) x,y,z,m
real(r8) conv,s,ddot
real(r8) coord(3,nat)

! Conversion factor from moments to rotational constants.
  conv= Planck / (8.0_r8 * pi *pi * pc_c)
! Add factor to put moments into SI units - give result in wavenumbers.
  conv=conv / (bohr2m * bohr2m * amu2kg * 100.0_r8)
!  conv=conv * au2ang*au2ang

mom=0.0_r8
!xyz=xyz*au2ang
do i=1,nat
x=xyz(1,i)
y=xyz(2,i)
z=xyz(3,i)
m=ams(iat(i))
mom(1,1)=mom(1,1)+(z**2+y**2)*m
mom(1,2)=mom(1,2)-x*y*m
mom(1,3)=mom(1,3)-x*z*m
mom(2,2)=mom(2,2)+(z**2+x**2)*m
mom(2,3)=mom(2,3)-y*z*m
mom(3,3)=mom(3,3)+(x**2+y**2)*m
enddo

!print*,xyz(1:3,1)
mom(2,1)=mom(1,2)
mom(3,1)=mom(1,3)
mom(3,2)=mom(2,3)

!do i=1,3
!write(stdout,'(3F14.6)') mom(i,1:3)
!enddo

! diag, mom contains now eigenvectors
call DiagSM(3,mom,e)
!print*, mom(1:3,1)
!print*, mom(1:3,2)
!print*, mom(1:3,3)

!mom(1:3,2)=-mom(1:3,2)
!mom(1:3,3)=-mom(1:3,3)

do i=1,3
if(e(i)<1e-6_r8) then
  rot(i)=0.0_r8
else
  rot(i)=conv/e(i)
endif
enddo

if(abs(e(1))<1e-6_r8) then
 write(stdout,'(3x,a)') '  linear molecule!'
 write(stdout,'(3x,a,2(F12.5,a))') ' Rotational constants: A= ---  B= ',rot(2) ,' C= ',rot(3), ' [cm^-1]'
 rot=rot*pc_c/10000.0_r8
 write(stdout,'(3x,a,2(F12.5,a))') ' Rotational constants: A= ---  B= ',rot(2) ,' C= ',rot(3), ' [Mhz]'
else
 write(stdout,'(3x,a,3(F12.5,a))') ' Rotational constants:  A= ', rot(1),' B= ',rot(2) ,' C= ',rot(3), ' [cm^-1]'
 rot=rot*pc_c/10000.0_r8
 write(stdout,'(3x,a,3(F12.5,a))') ' Rotational constants:  A= ', rot(1),' B= ',rot(2) ,' C= ',rot(3), ' [Mhz]'
endif


if(orient) then
write(stdout,'(a)')  '   ** rotating molecule to principle axis frame  **'
! handedness
s=mom(1,1)*(mom(2,2)*mom(3,3)-mom(3,2)*mom(2,3)) +  &
  mom(1,2)*(mom(2,3)*mom(3,1)-mom(2,1)*mom(3,3)) +  &
  mom(1,3)*(mom(2,1)*mom(3,2)-mom(2,2)*mom(3,1))
! invert if left-handed
if(s<0) then
  do i=1,3
  mom(i,1)=-mom(i,1)
  enddo
endif



coord=xyz
! rotate to principle axis frame
do i=1,nat
 do j=1,3
 xyz(j,i)=ddot(3,coord(1,i),1,mom(1,j),1)
 enddo
enddo

!xyz=xyz/au2ang
write(stdout,'(a)') '   ** writing orient_xopt.xyz  **'
call wrxyz('orient_xopt.xyz')
endif


end subroutine

integer pure function period(atnum)
! return period from atomic number
implicit none
integer, intent(in) :: atnum

select case(atnum)
 case(1:2)
   period=1
 case(3:10)
   period=2
 case(11:18)
   period=3
 case(19:36)
   period=4
 case(37:)
   period=5
end select 
end function


subroutine bondmatrix(nat,iat,xyz,bond,cn)
! bond matrix and coordinate number
use fiso, only: r8
use logic, only: debug
use atomdata
use constant, only: au2ang
implicit none
integer, intent(in) :: nat,iat(nat)
real(r8), intent(in) :: xyz(3,nat)
integer :: i,j,io
real(r8) :: r1,r2,dbond,fac_bond
character(2) :: esym
integer, intent(out) :: bond(nat,nat),cn(nat)

!fac_bond=1.1d0
fac_bond=1.0d0
! full symmetric bond matrix
bond=0
do i=1,nat
 do j=1,nat
   if(i==j) cycle
!   r1=(rcov(iat(i))+rcov(iat(j))) ! in A
   r1=(rcov3(iat(i))+rcov3(iat(j))) ! in A
   r2=dbond(xyz(1,i),xyz(1,j))*au2ang
   if(r2 <= r1*fac_bond) bond(i,j)=1
 enddo
enddo



if(debug) print*,'Writing xopt.bondmat'
 open(newunit=io,file='xopt.bondmat')
  call printimat(io,nat,nat,bond,'bond matrix')
!   write(io,*) bond
 close(io)

! coordination number (integer)
cn=0
do i=1,nat
 do j=1,nat
 if(bond(j,i)==1) then
  cn(i)=cn(i)+1
 endif
 enddo
enddo

if(debug) then
  write(*,'(a)')'CN:'
  do i=1,nat
  write(*,'(2x,I5,''['',a2,'']'',2x,I3)') i,esym(iat(i)),cn(i)
  enddo
endif

end subroutine

subroutine molbox(nat,xyz,iat,cell,vdw,echo)
use atomdata, only: rvdw
use fiso, only: r8,stdout
use logic, only: debug
use constant, only: au2ang
implicit none
integer :: nat, iat(nat),i
real(r8) :: xyz(3,nat),cell(3),vol
integer :: xmax(2),ymax(2),zmax(2),xmin(2),ymin(2),zmin(2)
logical echo, vdw


xmax=-999d0
ymax=-999d0
zmax=-999d0
xmin=999d0
ymin=999d0
zmin=999d0

do i=1,nat
 if(xyz(1,i)>xmax(1)) then
  xmax(1)=xyz(1,i)
  xmax(2)=i
 endif

 if(xyz(2,i)>ymax(1)) then
  ymax(1)=xyz(2,i)
  ymax(2)=i
 endif

 if(xyz(3,i)>zmax(1)) then
   zmax(1)=xyz(3,i)
   zmax(2)=i
 endif

 if(xyz(1,i)<xmin(1)) then
   xmin(1)=xyz(1,i)
   xmin(2)=i
 endif

 if(xyz(2,i)<ymin(1))then
   ymin(1)=xyz(2,i)
   ymin(2)=i
 endif

 if(xyz(3,i)<zmin(1)) then
  zmin(1)=xyz(3,i)
  zmin(2)=i
 endif
enddo


cell(1)=xmax(1)-xmin(1)
cell(2)=ymax(1)-ymin(1)
cell(3)=zmax(1)-zmin(1)

cell=cell*au2ang
if(echo) write(stdout,'(3x,a,3F8.2)') '     cell / A   : ',cell(1:3)
if(vdw) then
  ! add vdw radii of corner atoms
  cell(1)=cell(1)+rvdw(iat(xmax(2)))+rvdw(iat(xmin(2)))
  cell(2)=cell(2)+rvdw(iat(ymax(2)))+rvdw(iat(ymin(2)))
  cell(3)=cell(3)+rvdw(iat(zmax(2)))+rvdw(iat(zmin(2)))
endif
vol=cell(1)*cell(2)*cell(3)
if(echo) write(stdout,'(3x,a,3F8.2)') ' vdW cell / A   : ',cell(1:3)
if(echo) write(stdout,'(3x,a,F12.2)') ' vdW cell volume  / A^3 : ',vol

end subroutine


