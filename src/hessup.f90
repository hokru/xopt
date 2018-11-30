!*****************************
!* stuff for hessian updates *
!*****************************
! xvar is the local nvar cause here we can do cartesian and internal coords

!**********************
!* main wrapper       *
!**********************
subroutine Hupdate(iter,chess,hs,xvar,grad,oldG,displ,switch,doit)
! wrapper around varous Hessian update formulas
use fiso, only: r8
use popt
use logic
implicit none
integer iter
integer xvar,switch
real(r8) chess(xvar,xvar),hs(xvar*(xvar+1)/2)
real(r8) grad(xvar),oldG(xvar),displ(xvar)
logical doit

if(.not.doit) return

if(switch>=1)call packM(xvar,chess,hs,'pack')
! call screenMat(xvar,xvar,chess,1e-5_r8,.true.,'H(int)')
select case(iupdate)
 case default
 case(1)
   call SR1BFGS(xvar,grad,oldG,displ,hs)
 case(2)
  call SR1(xvar,grad,oldG,displ,hs)
 case(3)
  call BFGS(xvar,grad,oldG,displ,hs)
 case(4)
  call SS_BFGS(xvar,grad,oldG,displ,hs)
 case(5)
  call SS_SR1BFGS(xvar,grad,oldG,displ,hs)
 case(6)
  call SR1PSB(xvar,grad,oldG,displ,hs)
 case(7)
  call PSB(xvar,grad,oldG,displ,hs)
end select
if(switch>=1) call packM(xvar,chess,hs,'unpack')
end subroutine Hupdate

!*****************************
!* SR1 update                *
!*****************************
subroutine SR1(xvar,g,go,disp,h)
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) g(xvar),go(xvar),ddot,xx,u(xvar)
real(r8) h(xvar*(xvar+1)/2)
real(r8) dg(xvar),hdisp(xvar),disp(xvar)
real(r8) ud

dg=g-go

! H*displ
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

u=dg-hdisp
ud=ddot(xvar,u,1,disp,1)

do i=1,xvar
  do j=1,i
   ij = i*(i-1)/2 + j
    xx=(u(i)*u(j))/(ud)
    if(abs(xx).gt.epsilon(1.0d0)) h(ij)=h(ij)+xx
  enddo
enddo
end subroutine SR1


!*****************************
!* BFGS update               *
!*****************************
subroutine BFGS(xvar,g,go,disp,h)
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) g(xvar),go(xvar),ddot
real(r8) h(xvar*(xvar+1)/2)
real(r8) hdisp(xvar),dg(xvar),disp(xvar),dhd,dgdx

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
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) phi
real(r8) g(xvar),go(xvar),ddot,xx
real(r8) h(xvar*(xvar+1)/2)
real(r8) hdisp(xvar),dg(xvar),disp(xvar)
real(r8) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(r8) xud,xdhd

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
phi=sqrt(dxdu*dxdu/(dxdx*dudu))

if(abs(dxdx) < epsilon(1.0d0)) call error('no displacement! stoppping <sr1bfgs>')

xud=1.0_r8/ud
xdhd=1.0_r8/dhd
xx=1.0_r8/dgdx

do i=1,xvar
  do j=1,i
    ij = i*(i-1)/2 + j
    h(ij)= h(ij) + phi*(u(i)*u(j))*xud  ! SR1
    h(ij)= h(ij) + (1_r8-phi)*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
  enddo
enddo
end subroutine SR1BFGS

!***********************************************
!* SR1+PSB update for TS state searches        *
!* following DOI: 10.1063/1.1515483 (Helgakar) *
!* Hup=(1-phi)*SR1+phi*PSB                     *
!***********************************************
subroutine SR1PSB(xvar,g,go,disp,h)
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) phi
real(r8) g(xvar),go(xvar),ddot,xx
real(r8) h(xvar*(xvar+1)/2)
real(r8) hdisp(xvar),dg(xvar),disp(xvar)
real(r8) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(r8) xud,xdhd,term1,term2

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
phi=sqrt(dxdu*dxdu/(dxdx*dudu))

if(abs(dxdx) < epsilon(1d0)) call error('no displacement! stopping <sr1psb>')

xud=1.0_r8/ud
xdhd=1.0_r8/dhd
xx=1.0_r8/dgdx

do i=1,xvar
  do j=1,i
    ij = i*(i-1)/2 + j
    h(ij)= h(ij) + phi*(u(i)*u(j))*xud  ! SR1
    term1=(u(i)*disp(j)+disp(i)*u(j))/dxdx
    term2=(dxdu*disp(i)*disp(j))/dxdx*dxdx
    h(ij)=h(ij) + (1_r8-phi)*(term1-term2)! PSB
    ! term1=u(i)*disp(j)+disp(i)*u(j)
    ! term2=dxdu*disp(i)*disp(j)/dxdx
    ! h(ij)=h(ij) + (1_r8-phi)*(term1-term2)/dxdx! PSB
  enddo
enddo
end subroutine SR1PSB


subroutine packM(n,F,T,AA)
! pack/unpack square matrix to/from upper triangle
! a(i,j) is a(i+j*(j-1)/2) for (j.ge.i).
use fiso, only: r8
implicit none
integer i,j,k,n
real(r8) F(n,n),T(n*(n+1)/2)
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




subroutine packMF(n,F,T,AA)
! pack/unpack square matrix to/from upper triangle
! a(i,j) is a(i+j*(j-1)/2) for (j.ge.i).
use fiso, only: r4
implicit none
integer i,j,k,n
real(r4) F(n,n),T(n*(n+1)/2)
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
end subroutine packMF




!***********************************************
!* SR1+BFGS update                             *
!* following DOI: 10.1063/1.1515483 (Helgakar) *
!* Hup=(1-phi)*SR1+phi*BFGS                    *
!* single-precision version                    *
!***********************************************
subroutine fastSR1BFGS(xvar,xg,xgo,xdisp,h)
use fiso, only: r4, r8
implicit none
integer i,j,xvar,ij
real(r8) phi
real(r8) xg(xvar),xgo(xvar),xdisp(xvar)
real(r8) g(xvar),go(xvar),disp(xvar)
real(r4) sdot,xx
real(r8) h(xvar*(xvar+1)/2)
real(r4) hf(xvar*(xvar+1)/2)
real(r4) hdisp(xvar),dg(xvar)
real(r4) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(r4) xud,xdhd


call floatVec(xvar*(xvar+1)/2,h,hf)
call floatVec(xvar,xg,g)
call floatVec(xvar,xgo,go)
call floatVec(xvar,xdisp,disp)
dg=g-go
call sspmv('u',xvar,1.0,h,disp,1,0.0,hdisp,1)

u=dg-hdisp
ud=sdot(xvar,u,1,disp,1)

dxdu=sdot(xvar,disp,1,u,1)
dxdx=sdot(xvar,disp,1,disp,1)
dudu=sdot(xvar,u,1,u,1)
dgdx=sdot(xvar,dg,1,disp,1)
dhd=sdot(xvar,disp,1,hdisp,1)


! get weighing factor phi
phi=sqrt(dxdu*dxdu/(dxdx*dudu))

if(abs(dxdx) < epsilon(1.0)) call error('no displacement! stoppping')

xud=1.0_r4/ud
xdhd=1.0_r4/dhd
xx=1.0_r4/dgdx

do i=1,xvar
  do j=1,i
    ij = i*(i-1)/2 + j
    h(ij)= h(ij) + phi*(u(i)*u(j))*xud  ! SR1
    h(ij)= h(ij) + (1_r4-phi)*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
  enddo
enddo


call doubleVec(xvar*(xvar+1)/2,hf,h)
end subroutine fastSR1BFGS


subroutine hessdamp(nvar,hess,hold,gnorm)
! Hessian damping (not used)
use fiso, only: r8,stdout
implicit none
integer nvar
real(r8) Hd(nvar*(nvar+1)/2),hess(nvar*(nvar+1)/2),hold(nvar*(nvar+1)/2)
real(r8) damp,gnorm
integer ig

ig=3
if(gnorm.lt.0.1) ig=2
if(gnorm.lt.0.001) ig=1
if(gnorm.lt.0.0008) ig=3

select case(ig)
 case(1)
 damp=0.1
 case(2)
 damp=0.2_r8
 case(3)
 damp=0_r8
end select
write(stdout,'(''damping hessian: '',F6.3)') damp
Hd=(hess+damp*hold)/(1.0_r8+damp)
hess=Hd
end subroutine



!************************************************
!*SS-BFGS update (spectral scaling)             *
!************************************************
subroutine SS_BFGS(xvar,g,go,disp,h)
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) g(xvar),go(xvar),ddot
real(r8) h(xvar*(xvar+1)/2)
real(r8) hdisp(xvar),dg(xvar),disp(xvar),dhd,dgdx,sfac,dnrm2

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

!***********************************************
!* spectral scaling SR1+BFGS update            *
!* SS for phi seems to work well...            *
!***********************************************
subroutine SS_SR1BFGS(xvar,g,go,disp,h)
! experimental
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) phi
real(r8) g(xvar),go(xvar),ddot,xx
real(r8) h(xvar*(xvar+1)/2)
real(r8) hdisp(xvar),dg(xvar),disp(xvar)
real(r8) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(r8) xud,xdhd,sfac,dnrm2

dg=g-go
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

u=dg-hdisp
ud=ddot(xvar,u,1,disp,1)

dxdu=ddot(xvar,disp,1,u,1)
dxdx=ddot(xvar,disp,1,disp,1)
dudu=ddot(xvar,u,1,u,1)
dgdx=ddot(xvar,dg,1,disp,1)
dhd=ddot(xvar,disp,1,hdisp,1)

sfac=dgdx/(dnrm2(xvar,dg,1)**2)
!phi=sfac*dxdu*dxdu/(dxdx*dudu)
phi=sfac*dsqrt(dxdu*dxdu/(dxdx*dudu))

if(abs(dxdx) < epsilon(1d0)) call error('no displacement! stoppping <ss_sr1bfgs>')

xud=1.0_r8/ud
xdhd=1.0_r8/dhd
xx=1.0_r8/dgdx

do i=1,xvar
  do j=1,i
    ij = i*(i-1)/2 + j
    h(ij)= h(ij) + phi*(u(i)*u(j))*xud  ! SR1
    h(ij)= h(ij) + (1_r8-phi)*((dg(i)*dg(j)*xx) - hdisp(i)*hdisp(j)*xdhd) ! BFGS
  enddo
enddo
end subroutine


!***********************************************
!* Powell update for TS state searches         *
!***********************************************
subroutine psb(xvar,g,go,disp,h)
use fiso, only: r8
implicit none
integer i,j,xvar,ij
real(r8) phi
real(r8) g(xvar),go(xvar),ddot,xx
real(r8) h(xvar*(xvar+1)/2)
real(r8) hdisp(xvar),dg(xvar),disp(xvar)
real(r8) dxdu,dxdx,dudu,dgdx,dhd,u(xvar),ud
real(r8) xud,xdhd,term1,term2

dg=g-go
call dspmv('u',xvar,1.0d0,h,disp,1,0.0d0,hdisp,1)

u=dg-hdisp

dxdu=ddot(xvar,disp,1,u,1)
dxdx=ddot(xvar,disp,1,disp,1)

if(abs(dxdx) < epsilon(1d0)) call error('no displacement! stopping <sr1psb>')

! xud=1.0_r8/ud
! xdhd=1.0_r8/dhd
! xx=1.0_r8/dgdx

do i=1,xvar
  do j=1,i
      ij = i*(i-1)/2 + j
    !  h(ij)= h(ij) + ( u(i)*disp(j) + disp(i)*u(j) - disp(i)*dxdu*disp(j))/(dxdu/dxdx)
     h(ij)=h(ij) + (u(i)*disp(j)+disp(i)*u(j)/dxdx) - (disp(j)*u(i)*dxdx)/dxdx**2
  enddo
enddo
end subroutine 
