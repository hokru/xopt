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

!***********************
!* cartesian optimizer *
!***********************
subroutine copt
use parm
use logic
use popt
use progs
implicit none
integer iter, info
integer nat1,lwork
real(8), allocatable :: displ(:)
real(8), allocatable :: Uaug(:,:), eaug(:),aux(:)
real(8) f,dE,gmax,maxdipl,dnorm,oldE,oldG(3,nat)
character(2) esym
real(8) ddot,dnrm2
real(8) hs(nat*3*(nat*3+1)/2)
logical EE,GG,GM,DM,doupdate

EE=.false.
GG=.false.
GM=.false.

nat3=nat*3
nat1=nat3+1
lwork = 1 + 6*nat1 + 2*nat1**2

allocate(displ(nat3))
allocate(Uaug(nat1,nat1),eaug(nat1),aux(lwork))

energy=0d0
oldG=0d0
displ=0d0
Doupdate=.false.

! main loop
call wrxyz('xopt.log')
call touchlog

write(6,'(2x,a)') ''
write(6,'(2x,a)') 'Starting cartesian geometry optimization'
write(6,'(2x,a)') 'Iter          Energy         dE           gnorm         gmax           dnorm       maxdispl       lambda'
do iter=1,maxiter

! energy+gradient
oldE=energy
oldG=grad

call wrxyz(trim(xyzfile))
call getgrad
if(restrain) call addRestrainGrad
if (freeze) call freezeGrad
call loggrad(nat,iter,grad,energy)

dE=energy-oldE
gnorm=dnrm2(nat,grad(1,1),1)
gnorm=gnorm+dnrm2(nat,grad(2,1),1)
gnorm=gnorm+dnrm2(nat,grad(3,1),1)
gmax=maxval(grad)
call appxyz('xopt.log')

! update hessian
call Hupdate(iter,chess,hs,nat3,grad,oldG,displ,1,Doupdate)
!call wrhess(nat3,chess,'xopt.hess')
Doupdate=.true.
!if(iter>1) then
!  call packM(nat3,chess,hs,'pack')
! call SR1(nat3,grad,oldG,displ,hs)
! call BFGS(nat3,grad,oldG,displ,hs)
! call SR1BFGS(nat3,grad,oldG,displ,hs)
!  call bfgsms(nat3,gnorm,grad,oldG,displ,info,hs)
!  call packM(nat3,chess,hs,'unpack')
!endif


! form aug Hessian
Uaug = 0d0
do i=1,nat3
   do j=1,i
      Uaug(i,j)=chess(i,j)
      Uaug(j,i)=chess(i,j)
   enddo
enddo

k=0
do i=1,nat
 do j=1,3
  k=k+1
  Uaug(k,nat1)=grad(j,i)
  Uaug(nat1,k)=grad(j,i)
 enddo
 enddo
!print*, Uaug

call dsyev ('V','U',nat1,Uaug,nat1,eaug,aux,lwork,info)
if(info/=0) stop 'diag error'
if(abs(Uaug(nvar+1,1)).lt.1.d-10)stop'internal RF error'
displ(1:nat3)=Uaug(1:nat3,1)/Uaug(nat1,1)


! check step
do i=1,nat3
   if(abs(displ(i)).gt.maxd) then
      if(displ(i) < 0) displ(i)=-maxd
      if(displ(i) > 0) displ(i)= maxd
!      if(displ(i) > 0) then
!        print*,maxd,displ(i)
!        displ(i)= maxd
!        print*,maxd,displ(i),'X'
!      endif
   endif
enddo
dnorm=dnrm2(nat3,displ,1)
print*,dnrm2(nat3,displ,1),dnrm2(nat,displ,1)
maxdipl=maxval(displ)

k=0
do i=1,nat
 do j=1,3
 k=k+1
  xyz(j,i)=xyz(j,i)+displ(k)
 enddo
enddo



!print*,'COORDS'
!print*, xyz


! f=0.5291770d0
! open(unit=11,file='xopt.xyz')
! do i=1,nat
! write(11,'(a2,5x,3F18.14)') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
! enddo
!close(11)


if(abs(dE)<=econv) EE=.true.
if(gnorm<=gconv) GG=.true.
if(gmax<=maxgrad) GM=.true.
if(maxdipl<=dconv) DM=.true.

if(iter==1) call pciter(iter,energy,0d0,gnorm,gmax,dnorm,maxdipl,eaug(1),EE,GG,GM,DM)
if(iter>1) call pciter(iter,energy,dE,gnorm,gmax,dnorm,maxdipl,eaug(1),EE,GG,GM,DM)


if(EE.and.restrain.or.(EE.and.GG.and.GM.and.DM)) then
print*,'CONVERGED !'
write(6,'(3x,a)') '               criteria   actual value'
write(6,'(3x,a,2(ES10.3,2x))') ' energy    = ',econv,dE
write(6,'(3x,a,2(ES10.3,2x))') ' gnorm     = ',gconv,gnorm
write(6,'(3x,a,2(ES10.3,2x))') ' max grad  = ',maxgrad,gmax
write(6,'(3x,a,2(ES10.3,2x))') ' max displ = ',dconv,dnorm
if(restrain) call printFinalRestrain
  return
endif


if(iter==2000) then
  xyz0=xyz
  call getanc
  call ancopt(iter+1)
  return
endif


enddo ! end main loop
call wrhess(nat3,chess,'xopt.hess')
call printFinalRestrain

print*,'FAILED!'
return

end subroutine copt
