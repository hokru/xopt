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
!*****************************
!* stuff for hessian updates *
!*****************************
! xvar is the local nvar cause here we can do cartesian and internal coords

!**********************
!* main wrapper       *
!**********************
subroutine Hupdate(iter,chess,hs,xvar,grad,oldG,displ,switch,doit)
use popt
use logic
implicit none
integer iter,i
integer xvar,switch
real(8) chess(xvar,xvar),hs(xvar*(xvar+1)/2)
real(8) grad(xvar),oldG(xvar),displ(xvar)
logical doit

if(.not.doit) return
!if(iter==1) return
!if(iter<=3) return
!print*, "updating hessian"
!if(iter==5) then
!call getHess
!print*,'New guess Hessian obtained'
!endif

if(switch>=1)call packM(xvar,chess,hs,'pack')
!iupdate=3

select case(iupdate)
 case default
 case(1)
!  if(amber) then
!   call fastSR1BFGS(xvar,grad,oldG,displ,hs)
!  else
   call SR1BFGS(xvar,grad,oldG,displ,hs)
!  endif
 case(2)
  call SR1(xvar,grad,oldG,displ,hs)
 case(3)
  call BFGS(xvar,grad,oldG,displ,hs)
 case(4)
  call SS_BFGS(xvar,grad,oldG,displ,hs)
! case(5)
!   call SS_SR1BFGS(xvar,grad,oldG,displ,hs)
end select


if(switch>=1) call packM(xvar,chess,hs,'unpack')



end subroutine Hupdate

!*****************************
!* SR1 update                *
!*****************************
subroutine SR1(xvar,g,go,disp,h)
implicit none
integer i,j,k,xvar,ij
real(8) g(xvar),go(xvar),ddot,xx,u(xvar)
real(8) h(xvar*(xvar+1)/2),hdg(xvar)
real(8) dg(xvar),hdisp(xvar),disp(xvar)
real(8) ud

! grad=grad_old
!k=0
!do i=1,xvar/3
!   do j=1,3
!      k=k+1
!      dg(k)=g(j,i)-go(j,i)
!   enddo
!enddo

dg=g-go

! H*displ
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

u=dg-hdisp
ud=ddot(xvar,u,1,disp,1)

do i=1,xvar
  do j=1,i
   ij = i*(i-1)/2 + j  
    xx=(u(i)*u(j))/(ud)
    if(abs(xx).gt.epsilon(1d0)) h(ij)=h(ij)+xx
  enddo
enddo
end subroutine SR1


!*****************************
!* BFGS update               *
!*****************************
subroutine BFGS(xvar,g,go,disp,h)
implicit none
integer i,j,k,xvar,ij
real(8) g(xvar),go(xvar),ddot,xx,u(xvar)
real(8) h(xvar*(xvar+1)/2),hdg(xvar)
real(8) hdisp(xvar),dg(xvar),disp(xvar),dhd,dgdx

dg=g-go
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

dgdx=ddot(xvar,dg,1,disp,1)
dhd=ddot(xvar,disp,1,hdisp,1)

do i=1,xvar
  do j=1,i
   ij = i*(i-1)/2 + j
    h(ij)=h(ij) + (dg(i)*dg(j)/dgdx) - hdisp(i)*hdisp(j)/dhd
  enddo
enddo
end subroutine BFGS

!***********************************************
!* SR1+BFGS update                             *
!* following DOI: 10.1063/1.1515483 (Helgakar) *
!* Hup=(1-phi)*SR1+phi*BFGS                    *
!***********************************************
subroutine SR1BFGS(xvar,g,go,disp,h)
implicit none
integer i,j,k,xvar,ij
real(8) phi
real(8) g(xvar),go(xvar),ddot,xx
real(8) h(xvar*(xvar+1)/2),hdg(xvar)
real(8) hdisp(xvar),dg(xvar),disp(xvar)
real(8) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(8) xud,xdhd

dg=g-go
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

u=dg-hdisp
ud=ddot(xvar,u,1,disp,1)

dxdu=ddot(xvar,disp,1,u,1)
dxdx=ddot(xvar,disp,1,disp,1)
dudu=ddot(xvar,u,1,u,1)
dgdx=ddot(xvar,dg,1,disp,1)
dhd=ddot(xvar,disp,1,hdisp,1)


! get weighing factor phi
!if(abs(dxdx) > eps) phi=1d0-sqrt(dxdu*dxdu/(dxdx*dudu))
!phi=1d0-sqrt(dxdu*dxdu/(dxdx*dudu))
phi=dsqrt(dxdu*dxdu/(dxdx*dudu))

if(abs(dxdx) < epsilon(1d0)) stop 'no displacement! stoppping'

xud=1.0d0/ud
xdhd=1.0d0/dhd
xx=1.0d0/dgdx

do i=1,xvar
  do j=1,i
    ij = i*(i-1)/2 + j
!    h(ij)= h(ij) + (1d0-phi)*(u(i)*u(j))*xud  ! SR1
!    h(ij)= h(ij) + phi*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
    h(ij)= h(ij) + phi*(u(i)*u(j))*xud  ! SR1
    h(ij)= h(ij) + (1d0-phi)*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
  enddo
enddo
end subroutine SR1BFGS


subroutine packM(n,F,T,AA)
! pack/unpack square matrix to/from upper triangle
! a(i,j) is a(i+j*(j-1)/2) for (j.ge.i).
implicit none
integer i,j,k,n
real(8) F(n,n),T(n*(n+1)/2)
character(*) AA

select case(AA)
 case('pack')
k=0
do i=1,n
   do j=1,i
      k=k+1
      T(k)=F(j,i)
   enddo
enddo
case ('unpack')
k=0
do i=1,n
   do j=1,i
      k=k+1
      F(j,i)=T(k)
      F(i,j)=T(k)
   enddo
enddo
end select
end subroutine packM






!***********************************************
!* SR1+BFGS update                             *
!* following DOI: 10.1063/1.1515483 (Helgakar) *
!* Hup=(1-phi)*SR1+phi*BFGS                    *
!***********************************************
subroutine fastSR1BFGS(xvar,xg,xgo,xdisp,h)
implicit none
integer i,j,k,xvar,ij
real(8) phi
real(8) xg(xvar),xgo(xvar),xdisp(xvar)
real(4) g(xvar),go(xvar),disp(xvar)
real(4) sdot,xx
real(8) h(xvar*(xvar+1)/2)
real(4) hdg(xvar),hf(xvar*(xvar+1)/2)
real(4) hdisp(xvar),dg(xvar)
real(4) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(4) xud,xdhd


call floatVec(xvar*(xvar+1)/2,h,hf)
call floatVec(xvar,xg,g)
call floatVec(xvar,xgo,go)
call floatVec(xvar,xdisp,disp)
dg=g-go
call sspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

u=dg-hdisp
ud=sdot(xvar,u,1,disp,1)

dxdu=sdot(xvar,disp,1,u,1)
dxdx=sdot(xvar,disp,1,disp,1)
dudu=sdot(xvar,u,1,u,1)
dgdx=sdot(xvar,dg,1,disp,1)
dhd=sdot(xvar,disp,1,hdisp,1)


! get weighing factor phi
!if(abs(dxdx) > eps) phi=1d0-sqrt(dxdu*dxdu/(dxdx*dudu))
!phi=1d0-sqrt(dxdu*dxdu/(dxdx*dudu))
phi=sqrt(dxdu*dxdu/(dxdx*dudu))

if(abs(dxdx) < epsilon(1d0)) stop 'no displacement! stoppping'

xud=1.0d0/ud
xdhd=1.0d0/dhd
xx=1.0d0/dgdx

do i=1,xvar
  do j=1,i
    ij = i*(i-1)/2 + j
!    h(ij)= h(ij) + (1d0-phi)*(u(i)*u(j))*xud  ! SR1
!    h(ij)= h(ij) + phi*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
    h(ij)= h(ij) + phi*(u(i)*u(j))*xud  ! SR1
    h(ij)= h(ij) + (1d0-phi)*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
  enddo
enddo


call doubleVec(xvar*(xvar+1)/2,hf,h)
end subroutine fastSR1BFGS





!* Hessian damping
subroutine hessdamp(nvar,hess,hold,gnorm)
integer i,j
real(8) Hd(nvar*(nvar+1)/2),hess(nvar*(nvar+1)/2),hold(nvar*(nvar+1)/2)
real(8) damp,gnorm
integer ig

ig=3
if(gnorm.lt.0.1) ig=2
if(gnorm.lt.0.001) ig=1
if(gnorm.lt.0.0008) ig=3

select case(ig)
 case(1)  
 damp=0.1
 case(2)  
 damp=0.2d0
 case(3)
 damp=0d0
end select
write(*,'(''damping hessian: '',F6.3)') damp
Hd=(hess+damp*hold)/(1.0d0+damp)
hess=Hd
end subroutine 



!************************************************
!*SS-BFGS update (spectral scaling)             *
!************************************************
subroutine SS_BFGS(xvar,g,go,disp,h)
implicit none
integer i,j,k,xvar,ij
real(8) g(xvar),go(xvar),ddot,xx,u(xvar)
real(8) h(xvar*(xvar+1)/2),hdg(xvar)
real(8) hdisp(xvar),dg(xvar),disp(xvar),dhd,dgdx,sfac,dnrm2

dg=g-go
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

dgdx=ddot(xvar,dg,1,disp,1)
dhd=ddot(xvar,disp,1,hdisp,1)

sfac=dgdx/(dnrm2(xvar,dg,1)**2)
!print*,'sfac',sfac

do i=1,xvar
  do j=1,i
   ij = i*(i-1)/2 + j
    h(ij)=h(ij) + (sfac*dg(i)*dg(j)/dgdx) - hdisp(i)*hdisp(j)/dhd
  enddo
enddo
end subroutine

