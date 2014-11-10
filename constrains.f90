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
! handle constrains
! control file:
! $constr
! bond a b ref k
! angle a b c ref k
! dihed a b c d ref k
! 
!
!
module constr
integer nfbond ! number of fixed bonds
integer nfangle
integer nfdihed
real(8), allocatable :: rbond(:) ! reference values (distances)
real(8), allocatable :: rang(:) ! reference values 
real(8), allocatable :: rdih(:) ! reference values 
real(8), allocatable :: kbond(:) ! reference values (distances)
real(8), allocatable :: kang(:) ! reference values 
real(8), allocatable :: kdih(:) ! reference values 
integer, allocatable :: fbond(:,:) ! atom index, fixed bonds
integer, allocatable :: fangle(:,:) ! atom index, fixed bonds
integer, allocatable :: fdihed(:,:) ! atom index, fixed bonds

end module constr

! THIS IS NOT USED
subroutine XXgetconstrXX
use constr
implicit none
character(80) aa,cc
integer cb,ca,cd,n,s2i,nwords ! counters
real(8) s2r

nfbond=0
nfangle=0
nfdihed=0
! read xopt.control
open(44,file='xopt.control')

! read and count
 do
  read(44,'(a)',end=666) aa 
  if(index(aa,'$constr').ne.0) then
   do   
     read(44,'(a)',end=666) aa 
     call lower_case(aa)
     if(index(aa,'bond').ne.0)  nfbond=nfbond+1
     if(index(aa,'angle').ne.0) nfangle=nfangle+1
     if(index(aa,'dihed').ne.0) nfdihed=nfdihed+1
     exit ! break point
   enddo ! "$constr loop"
  endif
 enddo ! "xopt.control loop"
666 continue

! then allocate and rewind
rewind(44)
allocate(rbond(nfbond))
allocate(rang(nfbond))
allocate(rdih(nfbond))
allocate(kbond(nfbond))
allocate(kang(nfbond))
allocate(kdih(nfbond))
allocate(fbond(nfbond,2))
allocate(fangle(nfbond,3))
allocate(fdihed(nfbond,4))

cb=0
ca=0
cd=0
do
  read(44,'(a)',end=667) aa
  if(index(aa,'#').ne.0) cycle
  if(index(aa,'$constr').ne.0) then
   do
     read(44,'(a)',end=667) aa
     call lower_case(aa)
     if(index(aa,'bond').ne.0) then ! bonds
      cb=cb+1
      call charXsplit(aa,cc,2)
      fbond(cb,1)=s2i(cc)
      call charXsplit(aa,cc,3)
      fbond(cb,2)=s2i(cc)
      call charXsplit(aa,cc,4)
      rbond(cb)=s2i(cc)
      if(nwords(cc).gt.4) then
       call charXsplit(aa,cc,5)
       kbond(cb)=s2i(cc)
      endif
     endif
     if(index(aa,'angle').ne.0) then ! angles
      ca=ca+1 
      call charXsplit(aa,cc,2)
      fangle(ca,1)=s2i(cc)
      call charXsplit(aa,cc,3)
      fangle(ca,2)=s2i(cc)
      call charXsplit(aa,cc,4)
      fangle(ca,3)=s2i(cc)
      call charXsplit(aa,cc,5)
      rang(ca)=s2r(cc)
     endif
     if(index(aa,'dihed').ne.0) then ! dihedrals
      cd=cd+1
      call charXsplit(aa,cc,2)
      fdihed(cd,1)=s2i(cc)
      call charXsplit(aa,cc,3)
      fdihed(cd,2)=s2i(cc)
      call charXsplit(aa,cc,4)
      fdihed(cd,3)=s2i(cc)
      call charXsplit(aa,cc,5)
      fdihed(cd,4)=s2i(cc)
      call charXsplit(aa,cc,6)
      rdih(cd)=s2r(cc) 
     endif
    exit ! break point
   enddo ! "$constr loop"
  endif
 enddo ! "xopt.control loop"
667 continue



close(44)

!debug
print*,'constrains:'
print*,'bonds',fbond,rbond,kbond
print*,'angles',fangle,rang
print*,'diheds',fdihed,rdih

! setup variables


end subroutine 


! restrain distance between atoms a and b
! returns energy and gradient in a.u.
! V(x) = 0.5x**2 ; x=r_ab-r_ref
! a index = fbond(i,1)
! b index = fbond(i,2)
subroutine fixedbond(e,gbond)
use parm
use popt
use constr
implicit none
real(8), intent(out) :: e,gbond(nfbond)
real(8) rab,r(3)
real(8) xx,x,y,r2,dx,dy,dz

e=0
do i=1,nfbond
 dx=xyz(1,fbond(i,1))-xyz(1,fbond(i,2))
 dy=xyz(2,fbond(i,1))-xyz(2,fbond(i,2))
 dz=xyz(3,fbond(i,1))-xyz(3,fbond(i,2))
 rab=sqrt(dx**2+dy**2+dz**2)
 x=rab-rbond(i)
 y =  0.5d0 * kbond(i)*x*x
 print*,'e',y
 e=e + xx
 gbond(i)=gbond(i) + kbond(i)*x
enddo

end subroutine 



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


! freeze arbitrary atom numbers
subroutine freezeBmat
use parm
implicit none
integer icount

do i=1,nvar
   icount=0
   do j=1,nat
      do k=1,3 
         icount=icount+1
          if(ifrez(j)==1) b(icount,i)=0.0d0
      enddo 
   enddo 
 enddo 
end subroutine

subroutine freezeGrad
use parm
implicit none
do i=1,nat
 if(ifrez(i)==1) grad(1:3,i)=0.0d0
enddo
end subroutine


subroutine printFinalRestrain
use parm
use constant
implicit none
integer aa,bb,cc,dd
real(8) dih,dihgrad,di0,di0grad,r0(3),rab(3)
real(8) v0,vab,ang,ang0,agrad,agrad0,di360,anggrad

write(*,*) 
write(*,*) ' **  Restrain summary **'
do i=1,ires
! bond restraint
 if(irest_bond(i,1).gt.0.0d0) then
  aa=irest_bond(i,1)
  bb=irest_bond(i,2)
  r0(1:3)=xyz0(1:3,aa)-xyz0(1:3,bb)
  rab(1:3)=xyz(1:3,aa)-xyz(1:3,bb)
  call veclen(r0,v0)
  call veclen(rab,vab)
  write(*,'(a,1x,I4,I4,8x,a,F7.2,a,F7.2,a)') 'bond ', aa,bb, ' | start: ',v0*au2ang,' -> final: ',vab*au2ang,' angstrom'
 endif
 
 if(irest_ang(i,1).gt.0.0d0) then
  aa=irest_ang(i,1)
  bb=irest_ang(i,2)
  cc=irest_ang(i,3)
  call angle(xyz,aa,bb,cc,ang,agrad)
!  call angle(xyz0,aa,bb,cc,ang0,agrad0)
  ang0=val0(i)
  agrad0=anggrad(ang0)
  write(*,'(a,3(I4),3x,a,F7.2,a,F7.2,a)') 'angle ', aa,bb,cc, ' | start: ',agrad0,' -> final: ',agrad,' degree'
 endif


 if(irest_dihed(i,1).gt.0.0d0) then
 aa=irest_dihed(i,1)
 bb=irest_dihed(i,2)
 cc=irest_dihed(i,3)
 dd=irest_dihed(i,4)
 call dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
 di0=val0(i)
 di0grad=di360(di0)
! call dihed(xyz0,aa,bb,cc,dd,di0,di0grad)
!  write(*,'(a,4(I4),a,F7.2,a,F7.2,a)') 'dihed ', aa,bb,cc,dd, ' | start: ',di0grad,' -> final: ',dihgrad,' degree'
  write(*,'(a,4(I4),a,F7.2,a,F7.2,a,ES10.2)') 'dihed ', aa,bb,cc,dd, ' | start: ',di0grad,' -> final: ',dihgrad,' degree | E_r ',eneR(i)*627.509
  endif
enddo
write(*,*) 
write(*,'(a,F12.4)') 'Total E_restraint (kcal/mol): ',eneTOT*627.509
!write(*,'(a,F18.8)') 'Total E_restraint (Eh): ',eneTOT
write(*,*) 
write(*,*) 

end



subroutine addRestrainGrad
use parm
use constant
implicit none
integer nb,na,nd,aa,bb,cc,dd
real(8) r0(3),rab(3),rx,dbond,Konst,s,e(3),g(3),rgnorm,dnrm2
real(8) dih,dihgrad,di0,ang,agrad,ang0,rgrad(3,ires)
real(8) gd(3,4),ga(3,3),inv,inv2,va(3),vb(3),vc(3),eba(3),ebc(3),rba,rbc
real(8) co,si,xx,co2,si2,x2,x3,x4,crosa(3),crosb(3),rcd,x1,lab
real(8) v12(3),v23(3),v34(3),v43(3),v21(3),v32(3),vd(3),di0grad
logical qI,qIV,shift
real(8) di360
rgrad=0d0

! scaling factor s is currently not used
eneTOT=0d0
eneR=0d0
! loop over all restrains
do i=1,ires


!*****************
! bond restraint *
!*****************
 if(irest_bond(i,1).gt.0.0d0) then
  g=0d0
 aa=irest_bond(i,1)
 bb=irest_bond(i,2)
! S=irest_scale(ires,1)
 S=1.0d0


 r0(1:3)=xyz0(1:3,aa)-xyz0(1:3,bb)
 rab(1:3)=xyz(1:3,aa)-xyz(1:3,bb)
 rx=dbond(rab,r0)*0.5
 if(abs(rx).lt.1e-4) cycle
 call evec(rab,r0,e)
 konst=irest_konst(i,1)*s
 g(1)=Konst*rx*e(1)*2d0
 g(2)=Konst*rx*e(2)*2d0
 g(3)=Konst*rx*e(3)*2d0
 
 grad(1:3,aa)=grad(1:3,aa)+g(1:3)
 grad(1:3,bb)=grad(1:3,bb)-g(1:3)
 endif


!**********
!* angles *
!**********
 if(irest_ang(i,1).gt.0.0d0) then
  ga=0d0
  S=1.0d0
  aa=irest_ang(i,1)
  bb=irest_ang(i,2)
  cc=irest_ang(i,3)
  va(1:3)=xyz(1:3,aa)
  vb(1:3)=xyz(1:3,bb)
  vc(1:3)=xyz(1:3,cc)
!  call angle(xyz0,aa,bb,cc,ang0,agrad)
  ang0=val0(i)
  call angle(xyz,aa,bb,cc,ang,agrad)
  rx=ang-ang0
  if(abs(rx).lt.1e-4) cycle

  konst=irest_konst(i,2)*s



! ** numerical and analytical grad DIFFER BUT SIGN! dont know why... **
  call numangle(nat,xyz,xyz0,ga,aa,bb,cc,konst)
! goto 333
! print*,ga(1:3,1)
! ga=0
!  call evec(vb,va,eba)
!  call evec(vb,vc,ebc)
!  call veclen2(vb,va,rba)
!  call veclen2(vb,vc,rbc)
!  co=dot_product(eba,ebc)
!  si=sqrt(1.d0-co**2)
! 
!  ga(1:3,1)=(co*eba(1:3)-ebc(1:3))/(rba*si)
!  ga(1:3,2)=((rba-rbc*co)*eba(1:3)+(rbc-rba*co)*ebc(1:3))/(rba*rbc*si)
!  ga(1:3,3)=(co*ebc(1:3)-eba(1:3))/(rbc*si)
! 
!  ga=ga*konst*2d0*rx
! print*,ga(1:3,1)
! print*,'data',ang,ang0,rx,ang0
! 333 continue
  grad(1:3,aa)=grad(1:3,aa)+ga(1:3,1)
  grad(1:3,bb)=grad(1:3,bb)+ga(1:3,2)
  grad(1:3,cc)=grad(1:3,cc)+ga(1:3,3)
 endif

!*******************
!* dihedral angles *
!*******************
 if(irest_dihed(i,1).gt.0.0d0) then
  gd=0d0
 aa=irest_dihed(i,1)
 bb=irest_dihed(i,2)
 cc=irest_dihed(i,3)
 dd=irest_dihed(i,4)
 S=1.0d0

 call dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
! call dihed(xyz0,aa,bb,cc,dd,di0,di0grad)
 di0=val0(i)
 di0grad=di360(di0)

! write(33,'(2x,a,4I4,2x,F9.5,2x,F)') 'dihed ',aa,bb,cc,dd,irest_konst(ires,3),dih

! Check for the case that the torsion goes from the I to the IV quadrant and adjust accordingly

!  IV to I
qI=.false.
qIV=.false.
shift=.false.
 if(dihgrad.gt.0.and.dihgrad.lt.90) qI=.true.
 if(di0grad.gt.270.and.di0grad.lt.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
!  print*,'case1',dihgrad,di0grad
  dih=dih+2.0d0*pi
  rx=dih-di0
  goto 999
 endif

 qI=.false.
 qIV=.false.

! I to VI
qI=.false.
qIV=.false.
shift=.false.
 if(di0grad.gt.0.and.di0grad.lt.90) qI=.true.
 if(dihgrad.gt.270.and.dihgrad.lt.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
!   print*,'case2',dihgrad,di0grad
  dih=dih-2.0d0*pi
  rx=dih-di0
  goto 999
 endif

999 continue

 if(.not.shift)  rx=dih-di0
 if(abs(rx).lt.1e-4) cycle
  
 konst=irest_konst(i,3)*s


! ** numerical and analytical gradient agree **
! call numdihedgrad(nat,xyz,xyz0,gd,aa,bb,cc,dd,konst)

! print*,'num',gd
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
  si=sqrt(1.d0-co**2)
  si2=sqrt(1.d0-co2**2)

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

     gd=gd*konst*2d0*rx
!     print*,'ana',gd
!     print*,'///////'
 grad(1:3,aa)=grad(1:3,aa)+gd(1:3,1)
 grad(1:3,bb)=grad(1:3,bb)+gd(1:3,2)
 grad(1:3,cc)=grad(1:3,cc)+gd(1:3,3)
 grad(1:3,dd)=grad(1:3,dd)+gd(1:3,4)
 endif

! energy of restraints
 eneR(i)=konst*rx*rx
! print*,eneR(i),rx,konst

enddo

 eneTOT=sum(eneR)
! print*,eneTOT
return

end subroutine

subroutine unitvec(x,e)
implicit none
real(8) x(3),e(3),t(3)
t=DOT_PRODUCT(x,x)
!t=x(1)**2+x(2)**2+x(3)**2
e=x/sqrt(t)
end

subroutine veclen(x,v)
implicit none
real(8) x(3),v
v=dsqrt(dot_product(x,x))
end subroutine


subroutine veclen2(a,b,v)
implicit none
real(8) a(3),b(3),v,x(3)
x=a-b
v=dsqrt(dot_product(x,x))
end subroutine


subroutine cross_prod(y,x2,x3)
implicit none
real(8) y(3),x2(3),x3(3)
  y(1) =  x2(2)*x3(3) - x3(2)*x2(3)
  y(2) = -x2(1)*x3(3) + x3(1)*x2(3)
  y(3) =  x2(1)*x3(2) - x3(1)*x2(2)
end subroutine


subroutine angle(xyz,aa,bb,cc,ang,anggrad)
use constant
implicit none
integer aa,bb,cc
real(8) xyz(3,*), v1(3),v2(3),lv1,lv2,uv1(3),uv2(3)
real(8) t, tt,ang,anggrad
real(8) cv,lcv

v1=xyz(1:3,bb)-xyz(1:3,aa)
v2=xyz(1:3,bb)-xyz(1:3,cc)

call cross_prod(cv,v1,v2)
call veclen(cv,lcv)
ang=atan2(lcv, dot_product(v1,v2))
anggrad=ang*180.0d0/pi
!print*,anggrad
!if (ang.lt.0.0d0) ang=ang+2*pi

!t=DOT_PRODUCT(v1,v2)
!call veclen(v1,lv1)
!call veclen(v2,lv2)
!tt=lv1*lv2
!ang=acos(t/tt)
!anggrad=ang*180.0d0/pi
!print*,anggrad
end subroutine

subroutine dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
use constant
implicit none
integer aa,bb,cc,dd
real(8) b1(3),b2(3),b3(3),n1(3),n2(3)
real(8) un1(3),un2(3),ub2(2),m(3),um(3),dix,diy
real(8) dih,dihgrad,xyz(3,*)

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


!stop
!  print*,'DIX,DIY',dix,diy


! give results in 0 to 360 degree
 if(dih<0.0d0) dih=dih+pi*2
 dihgrad=dih*180.0d0/pi
end subroutine


subroutine numdihedgrad(nat,xyz,xyz0,grad,aa,bb,cc,dd,konst)
implicit none
real(8) xyz(3,nat),xyz0(3,nat),grad(3,4)
integer i,j,aa,bb,cc,dd,at(4),nat
real(8) step,dih,dummy,dih0,rx,di0,konst,er,el

step=0.001d0
grad=0d0

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
    xyz(j,at(i))=xyz(j,at(i))-2.0d0*step
      call dihed(xyz,aa,bb,cc,dd,dih,dummy)
      rx=dih-di0
      el=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))+step
    grad(j,i)=(er-el)/(2.0d0*step)
   enddo ! j-loop
 enddo  ! i-loop

end subroutine


subroutine numangle(nat,xyz,xyz0,grad,aa,bb,cc,konst)
implicit none
real(8) xyz(3,nat),xyz0(3,nat),grad(3,3)
integer i,j,aa,bb,cc,dd,at(3),nat
real(8) step,dih,dummy,ang,rx,ang0,konst,er,el

at(1)=aa
at(2)=bb
at(3)=cc
step=0.001d0
grad=0d0

call  angle(xyz0,aa,bb,cc,ang0,dummy)
do i=1,3
 do j=1,3
    xyz(j,at(i))=xyz(j,at(i))+step
      call  angle(xyz,aa,bb,cc,ang,dummy)
      rx=ang-ang0
      er=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))-2.0d0*step
      call  angle(xyz,aa,bb,cc,ang,dummy)
      rx=ang-ang0
      el=konst*rx*rx
    xyz(j,at(i))=xyz(j,at(i))+step
    grad(j,i)=(er-el)/(2.0d0*step)
 enddo
enddo
end subroutine

real*8 function di360(x)
implicit none
real(8), parameter:: pi = 3.141592653589793d0
real(8) x
if(x<0.0d0) x=x+pi*2
di360=x*180.0d0/pi
end function

real(8) function anggrad(x)
implicit none
real(8) x
real(8), parameter:: pi = 3.141592653589793d0
anggrad=x*180.0d0/pi
end function
