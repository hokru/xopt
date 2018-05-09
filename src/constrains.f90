
! module constr
! ! not used
! integer nfbond ! number of fixed bonds
! integer nfangle
! integer nfdihed
! real(8), allocatable :: rbond(:) ! reference values (distances)
! real(8), allocatable :: rang(:) ! reference values
! real(8), allocatable :: rdih(:) ! reference values
! real(8), allocatable :: kbond(:) ! reference values (distances)
! real(8), allocatable :: kang(:) ! reference values
! real(8), allocatable :: kdih(:) ! reference values
! integer, allocatable :: fbond(:,:) ! atom index, fixed bonds
! integer, allocatable :: fangle(:,:) ! atom index, fixed bonds
! integer, allocatable :: fdihed(:,:) ! atom index, fixed bonds
!
! end module constr

! restrain distance between atoms a and b
! returns energy and gradient in a.u.
! V(x) = 0.5x**2 ; x=r_ab-r_ref
! a index = fbond(i,1)
! b index = fbond(i,2)
! subroutine fixedbond(e,gbond)
! ! not used
! use parm, only: xyz,i
! use popt
! use constr
! implicit none
! real(8), intent(out) :: e,gbond(nfbond)
! real(8) rab,r(3)
! real(8) xx,x,y,r2,dx,dy,dz
!
! e=0
! do i=1,nfbond
!  dx=xyz(1,fbond(i,1))-xyz(1,fbond(i,2))
!  dy=xyz(2,fbond(i,1))-xyz(2,fbond(i,2))
!  dz=xyz(3,fbond(i,1))-xyz(3,fbond(i,2))
!  rab=sqrt(dx**2+dy**2+dz**2)
!  x=rab-rbond(i)
!  y =  0.5d0 * kbond(i)*x*x
!  ! print*,'e',y
!  e=e + xx
!  gbond(i)=gbond(i) + kbond(i)*x
! enddo
!
! end subroutine



! freeze all non-H atoms
!subroutine freezeHB(nat,nvar,b,it)
!implicit none
!real(8), intent(inout):: b(nat*3,nat*3)
!integer nat,nvar,it(nat),i,j,k,icount
!
!do i=1,nvar
!   icount=0
!   do j=1,nat
!      do k=1,3
!         icount=icount+1
!          if(it(j).ne.1)  b(icount,i)=0.0d0
!      enddo
!   enddo
! enddo
!end subroutine



subroutine freezeBmat
! freeze arbitrary atom numbers
use fiso, only: r8
use parm, only: nvar,ifrez,i,j,k,nat,B
implicit none
integer icount

do i=1,nvar
   icount=0
   do j=1,nat
      do k=1,3
         icount=icount+1
          if(ifrez(j)==1) b(icount,i)=0.0_r8
      enddo
   enddo
 enddo
end subroutine

subroutine freezeGrad
use fiso, only: r8
use parm, only: nat,ifrez,i,grad
implicit none
do i=1,nat
 if(ifrez(i)==1) grad(1:3,i)=0.0_r8
enddo
end subroutine


subroutine printFinalRestrain
use parm
use fiso, only: r8,stdout
use constant,only: au2ang,au2kcal
implicit none
integer aa,bb,cc,dd
real(r8) dih,dihgrad,di0,di0grad,r0(3),rab(3)
real(r8) v0,vab,ang,ang0,agrad,agrad0,di360,anggrad

write(stdout,*) ' '
write(stdout,*) ' **  Restraint Summary **'
do i=1,ires


! space restraints
if(irest_atom(i)>0) then
  aa=irest_atom(i)
  r0(1:3)=xyz0(1:3,aa)-xyz0(1:3,aa)
  rab(1:3)=xyz(1:3,aa)-xyz0(1:3,aa)
  call veclen(r0,v0)
  call veclen(rab,vab)
  write(stdout,'(a,1x,I4,8x,a,F7.2,a,F7.2,a,ES10.2)') &
  'atom ', aa, ' | start: ',v0*au2ang,' -> final: ',vab*au2ang,' angstrom | E_r ',eneR(i)*au2kcal
endif

! bond restraints
 if(irest_bond(i,1)>0) then
  aa=irest_bond(i,1)
  bb=irest_bond(i,2)
  r0(1:3)=xyz0(1:3,aa)-xyz0(1:3,bb)
  rab(1:3)=xyz(1:3,aa)-xyz(1:3,bb)
  call veclen(r0,v0)
  call veclen(rab,vab)
  write(stdout,'(a,1x,I4,I4,8x,a,F7.2,a,F7.2,a,ES10.2)') &
  'bond ', aa,bb, ' | start: ',v0*au2ang,' -> final: ',vab*au2ang,' angstrom | E_r ',eneR(i)*au2kcal
 endif

! angle restraints
 if(irest_ang(i,1)>0) then
  aa=irest_ang(i,1)
  bb=irest_ang(i,2)
  cc=irest_ang(i,3)
  call angle(xyz,aa,bb,cc,ang,agrad)
  ang0=val0(i)
  agrad0=anggrad(ang0)
  write(stdout,'(a,3(I4),3x,a,F7.2,a,F7.2,a,ES10.2)') &
  'angle ', aa,bb,cc, ' | start: ',agrad0,' -> final: ',agrad,' degree | E_r ',eneR(i)*au2kcal
 endif

! dihedral restraints
 if(irest_dihed(i,1)>0) then
  aa=irest_dihed(i,1)
  bb=irest_dihed(i,2)
  cc=irest_dihed(i,3)
  dd=irest_dihed(i,4)
  call dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
  di0=val0(i)
  di0grad=di360(di0)
  write(stdout,'(a,4(I4),a,F7.2,a,F7.2,a,ES10.2)') &
  'dihed ', aa,bb,cc,dd, ' | start: ',di0grad,' -> final: ',dihgrad,' degree | E_r ',eneR(i)*au2kcal
 endif

if(irest_vol(i,1)>0) then ! sphere volume
  write(stdout,'(a,2x,F7.2,a,F10.5)') 'sphere ',irest_vol(ires,1),' [A] | E_r ',eneR(i)*au2kcal
endif

enddo
write(stdout,*)''
write(stdout,'(a,F14.6)') 'Total E_restraint (kcal/mol): ',eneTOT*au2kcal
write(stdout,*)''
write(stdout,*)''

end



subroutine addRestrainGrad
use parm
use constant
use fiso, only: r8
implicit none
integer aa,bb,cc,dd
real(r8) r0,rab,rx,dbond,Konst,s,e(3),g(3),ee
real(r8) dih,dihgrad,di0,ang,agrad,ang0,rgrad(3,ires)
real(r8) gd(3,4),ga(3,3),va(3),vb(3),vc(3),rbc
real(r8) co,si,co2,si2,x2,x3,x4,crosa(3),crosb(3),rcd,x1,lab
real(r8) v12(3),v23(3),v34(3),v43(3),v21(3),v32(3),vd(3),di0grad
real(r8) rba,eba(3),ebc(3)
logical qI,qIV,shift
real(r8) di360
rgrad=0.0_r8

! scaling factor s is currently not used
eneTOT=0.0_r8
eneR=0.0_r8
rx=0_r8
! loop over all restrains
do i=1,ires


!********************
! volume restraints *
!********************
if(irest_vol(i,1)>0) then ! sphere volume
 r0=irest_vol(i,1)
 konst=irest_konst(i,5)
 do k=1,nat
    g=0d0; ee=0d0
    e(1:3)=xyz(1:3,k)*au2ang
    call veclen(e,lab)
    if(lab>=r0) then
      rx=(lab-r0)
      g(1)=konst*rx*2.0_r8*e(1)/lab
      g(2)=konst*rx*2.0_r8*e(2)/lab
      g(3)=konst*rx*2.0_r8*e(3)/lab
      grad(1:3,i)=grad(1:3,i)+g(1:3)
      ee=ee+konst+rx+rx
    endif
 enddo
 eneR(i)=ee
endif

!*******************
! space restraints *
!********************
if(irest_atom(i).gt.0) then
 ! e = konst*rx^2
 aa=irest_atom(i)
 S=1.0_r8
 g=0.0_r8
 r0=val0(i)
 rab=dbond(xyz(1:3,aa),xyz0(1:3,aa))
 e(1:3)=xyz(1:3,aa)-xyz0(1:3,aa)
 rx=(rab-r0)
 if(abs(rx).lt.1e-4) cycle
 konst=irest_konst(i,4)*s
 call veclen(e,lab)
 g(1)=konst*rx*2.0_r8*e(1)/lab
 g(2)=konst*rx*2.0_r8*e(2)/lab
 g(3)=konst*rx*2.0_r8*e(3)/lab
 grad(1:3,aa)=grad(1:3,aa)+g(1:3)
 eneR(i)=konst*rx*rx
! print*,aa,rx,eneR(i)
endif

!*****************
! bond restraint *
!*****************
if(irest_bond(i,1).gt.0) then
  g=0.0_r8
 aa=irest_bond(i,1)
 bb=irest_bond(i,2)
! S=irest_scale(ires,1)
 S=1.0_r8


 r0=val0(i)
 rab=dbond(xyz(1:3,aa),xyz(1:3,bb))
 rx=(rab-r0)

 if(abs(rx).lt.1e-4) cycle
 e(1:3)=xyz(1:3,aa)-xyz(1:3,bb)
 call veclen(e,lab)
 e=e/lab
 konst=irest_konst(i,1)*s
 g(1)=Konst*rx*e(1)*2.0_r8
 g(2)=Konst*rx*e(2)*2.0_r8
 g(3)=Konst*rx*e(3)*2.0_r8

 grad(1:3,aa)=grad(1:3,aa)+g(1:3)
 grad(1:3,bb)=grad(1:3,bb)-g(1:3)

 eneR(i)=konst*rx*rx
endif


!**********
!* angles *
!**********
 if(irest_ang(i,1).gt.0) then
  ga=0.0_r8
  s=1.0_r8
  aa=irest_ang(i,1)
  bb=irest_ang(i,2)
  cc=irest_ang(i,3)
  va(1:3)=xyz(1:3,aa)
  vb(1:3)=xyz(1:3,bb)
  vc(1:3)=xyz(1:3,cc)
  ang0=val0(i)
  call angle(xyz,aa,bb,cc,ang,agrad)
  rx=ang-ang0
  if(abs(rx).lt.1e-4) cycle

  konst=irest_konst(i,2)*s



  ! call numangle(nat,xyz,xyz0,ga,aa,bb,cc,konst)
  call evec(va,vb,eba)
  call evec(vc,vb,ebc)
  call veclen2(va,vb,rba)
  call veclen2(vc,vb,rbc)
  co=dot_product(eba,ebc)
  si=sqrt(1.0_r8-co**2)
  ga(1:3,1)=(co*eba(1:3)-ebc(1:3))/(rba*si)
  ga(1:3,2)=((rba-rbc*co)*eba(1:3)+(rbc-rba*co)*ebc(1:3))/(rba*rbc*si)
  ga(1:3,3)=(co*ebc(1:3)-eba(1:3))/(rbc*si)
  ga=ga*konst*2.0_r8*rx
  grad(1:3,aa)=grad(1:3,aa)+ga(1:3,1)
  grad(1:3,bb)=grad(1:3,bb)+ga(1:3,2)
  grad(1:3,cc)=grad(1:3,cc)+ga(1:3,3)

  eneR(i)=konst*rx*rx
 endif

!*******************
!* dihedral angles *
!*******************
 if(irest_dihed(i,1).gt.0) then
  gd=0.0_r8
 aa=irest_dihed(i,1)
 bb=irest_dihed(i,2)
 cc=irest_dihed(i,3)
 dd=irest_dihed(i,4)
 S=1.0_r8

 call dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
 di0=val0(i)
 di0grad=di360(di0)

! write(*,'(2x,a,4I4,2x,F9.5,2x,2(F6.2,2x))') 'dihed ',aa,bb,cc,dd,irest_konst(ires,3),dihgrad,di0grad

! Check for the case that the torsion goes from the I to the IV quadrant and adjust accordingly
! we check for "greater/less equal" since we might want to reach 0 as target value
!  IV to I
qI=.false.
qIV=.false.
shift=.false.
 if(dihgrad.ge.0.and.dihgrad.le.90) qI=.true.
 if(di0grad.ge.270.and.di0grad.le.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
  dih=dih+2.0_r8*pi
  rx=dih-di0
!  print*,'case1',dihgrad,di0grad,rx
  goto 999
 endif

 qI=.false.
 qIV=.false.

! I to VI
qI=.false.
qIV=.false.
shift=.false.
 if(di0grad.ge.0.and.di0grad.le.90) qI=.true.
 if(dihgrad.ge.270.and.dihgrad.le.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
  dih=dih-2.0_r8*pi
  rx=dih-di0
!   print*,'case2',dihgrad,di0grad,rx
  goto 999
 endif

999 continue

 if(.not.shift)  rx=dih-di0
 if(abs(rx).lt.1e-4) cycle

 konst=irest_konst(i,3)*s

! ** numerical and analytical gradient agree **
! call numdihedgrad(nat,xyz,xyz0,gd,aa,bb,cc,dd,konst)
 va(1:3)=xyz(1:3,aa)
 vb(1:3)=xyz(1:3,bb)
 vc(1:3)=xyz(1:3,cc)
 vd(1:3)=xyz(1:3,dd)

 call evec(va,vb,v12)
 call evec(vb,vc,v23)
 call evec(vc,vd,v34)
 call evec(vd,vc,v43)
 call evec(vc,vb,v32)
 call evec(vb,va,v21)
 co=DOT_PRODUCT(v21,v23)
 co2=DOT_PRODUCT(v32,v34)
 si=sqrt(1.0_r8-co**2)
 si2=sqrt(1.0_r8-co2**2)

 call cross_prod(crosa,v12,v23)
 call cross_prod(crosb,v43,v32)
 call veclen2(va,vb,lab)
 call veclen2(vb,vc,rbc)
 call veclen2(vc,vd,rcd)
 x1=(rbc-lab*co)/(rbc*lab*si**2)
 x2=co2/(rbc*si2**2)
 x3=(rbc-rcd*co2)/(rbc*rcd*si2**2)
 x4=co/(rbc*si**2)

 gd(1:3,1)=-crosa(1:3)/(lab*si**2)
 gd(1:3,2)=x1*crosa(1:3)+x2*crosb(1:3)
 gd(1:3,3)=x3*crosb(1:3)+x4*crosa(1:3)
 gd(1:3,4)=-crosb(1:3)/(rcd*si2**2)

 gd=gd*konst*2.0_r8*rx
 grad(1:3,aa)=grad(1:3,aa)+gd(1:3,1)
 grad(1:3,bb)=grad(1:3,bb)+gd(1:3,2)
 grad(1:3,cc)=grad(1:3,cc)+gd(1:3,3)
 grad(1:3,dd)=grad(1:3,dd)+gd(1:3,4)
 endif

! energy of restraints
 eneR(i)=konst*rx*rx
enddo
eneTOT=sum(eneR)
end subroutine



subroutine angle(xyz,aa,bb,cc,ang,anggrad)
use constant, only: pi
use fiso, only: r8
implicit none
integer aa,bb,cc
real(r8) xyz(3,*),v1(3),v2(3)
real(r8) ang,anggrad
real(r8) cv(3),lcv

v1=xyz(1:3,bb)-xyz(1:3,aa)
v2=xyz(1:3,bb)-xyz(1:3,cc)

call cross_prod(cv,v1,v2)
call veclen(cv,lcv)
ang=atan2(lcv, dot_product(v1,v2))
anggrad=ang*180.0_r8/pi
end subroutine

subroutine dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
use constant, only: pi
use fiso, only: r8
implicit none
integer aa,bb,cc,dd
real(r8) b1(3),b2(3),b3(3),n1(3),n2(3)
real(r8) un1(3),un2(3),ub2(3),m(3),um(3),dix,diy
real(r8) dih,dihgrad,xyz(3,*)

 b1=xyz(1:3,aa)-xyz(1:3,bb)
 b2=xyz(1:3,bb)-xyz(1:3,cc)
 b3=xyz(1:3,cc)-xyz(1:3,dd)

 ! normal of the planes
 call cross_prod(n1,b1,b2)
 call cross_prod(n2,b2,b3)

 call unitvec(n1,un1)
 call unitvec(n2,un2)
 call unitvec(b2,ub2)

 call cross_prod(m,un1,ub2)
 call unitvec(m,um)

 dix=DOT_PRODUCT(un1,un2)
 diy=DOT_PRODUCT(um,un2)

 dih=atan2(diy,dix)

!  Quadrant    Angle              sin    cos    tan
!----------------------------------------------------
!  I           0    < α < π/2     > 0    > 0    > 0
!  II          π/2  < α < π       > 0    < 0    < 0
!  III         π    < α < 3π/2    < 0    < 0    > 0
!  IV          3π/2 < α < 2π      < 0    > 0    < 0
! atan2(0,1) =   0
! atan2(1,0) =   pi/2
! atan2(-1,0) = -pi/2
! atan2(0,-1) =  pi

!  print*,atan2(1.0,1.0)*180.0d0/pi
!  print*,atan2(-1.0,0.0) *180.0d0/pi

! give results in 0 to 360 degree
 if(dih<0.0_r8) dih=dih+pi*2
 dihgrad=dih*180.0_r8/pi
end subroutine


subroutine numdihedgrad(nat,xyz,xyz0,grad,aa,bb,cc,dd,konst)
use fiso, only: r8
implicit none
integer i,j,aa,bb,cc,dd,at(4),nat
real(r8) xyz(3,nat),xyz0(3,nat),grad(3,4)
real(r8) step,dih,dummy,rx,di0,konst,er,el

step=0.001_r8
grad=0.0_r8

call dihed(xyz0,aa,bb,cc,dd,di0,dummy)

at(1)=aa
at(2)=bb
at(3)=cc
at(4)=dd

 do i=1,4
   do j=1,3
    xyz(j,at(i))=xyz(j,at(i))+step
      call dihed(xyz,aa,bb,cc,dd,dih,dummy)
      rx=dih-di0
      er=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))-2.0_r8*step
      call dihed(xyz,aa,bb,cc,dd,dih,dummy)
      rx=dih-di0
      el=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))+step
    grad(j,i)=(er-el)/(2.0_r8*step)
   enddo ! j-loop
 enddo  ! i-loop

end subroutine


subroutine numangle(nat,xyz,xyz0,grad,aa,bb,cc,konst)
use fiso, only: r8
implicit none
integer i,j,aa,bb,cc,at(3),nat
real(r8) xyz(3,nat),xyz0(3,nat),grad(3,3)
real(r8) step,dummy,ang,rx,ang0,konst,er,el

at(1)=aa
at(2)=bb
at(3)=cc
step=0.001_r8
grad=0.0_r8

call  angle(xyz0,aa,bb,cc,ang0,dummy)
do i=1,3
 do j=1,3
    xyz(j,at(i))=xyz(j,at(i))+step
    call  angle(xyz,aa,bb,cc,ang,dummy)
    rx=ang-ang0
    er=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))-2.0_r8*step
    call  angle(xyz,aa,bb,cc,ang,dummy)
    rx=ang-ang0
    el=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))+step
    grad(j,i)=(er-el)/(2.0_r8*step)
 enddo
enddo
end subroutine

real(r8) function di360(x)
use fiso, only: r8
use constant, only: pi
implicit none
real(r8) x
if(x<0.0_r8) x=x+pi*2
di360=x*180.0_r8/pi
end function

real(r8) function grad2rad(x)
use constant, only: pi
use fiso, only: r8
implicit none
real(r8) x
grad2rad=x*pi/180.0_r8
end function

real(r8) function anggrad(x)
use fiso, only: r8
use constant, only: pi
implicit none
real(r8) x
anggrad=x*180.0_r8/pi
end function
