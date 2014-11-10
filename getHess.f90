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
subroutine getHess
use parm
use logic
implicit none
real(8) be

!print*,'BE: ', be

select case(hmodel)
 case('numerical')
 call numhess
 case('bond')
!  call bondhess
 case('lindh')
  stop 'sorry. not implemented'
 case('hdiag')
  chess=0.d0
  do i=1,3*nat
   do j=1,3*nat
    chess(i,j)=0.0
    if(i==j) chess(i,j)=1.0
   enddo
  enddo
 case('scaled')
 case default
end select 

end subroutine getHess



subroutine modlindh
use parm
implicit none
real(8) dx,dy,dz,dis

do i=1,nat-1
 do j=i+1,nat
   dx=xyz(1,i)-xyz(1,j)
   dy=xyz(2,i)-xyz(2,j)
   dz=xyz(3,i)-xyz(3,j)
   dis=sqrt(dx*dx+dy*dy+dz*dz)

   

 enddo
enddo

end subroutine

real(8) function be
use parm
implicit none
be=1
!be=0
!do i=1,npair
!" r=dist(i)
! be=be+0.1*exp(-0.5**r**2)
!enddo

end function



subroutine numhess
! numerical cartesian hessian
use parm
use logic
use progs
implicit none
integer ii,jj,hi,hj,yy
real*8  step
real*8  t1,t0,xx,thrE
real*8, allocatable ::  hess(:,:), gr(:,:),gl(:,:)


thrE=epsilon(1d0)
step=0.005d0
  print*, 'Doing Hessian numerically ..'
  allocate(hess(nat3,nat3),gl(3,nat),gr(3,nat))
  hess=0.0d0
  do i=1,nat
  if(int(mod(i,nint(nat/10.0))).eq.0) write(6,'(I3,A)') nint(100d0*i/dble(nat)),' % done'
   do j=1,3
    hi=(i-1)*3+j
    xyz(j,i)=xyz(j,i)+step
    call newxyz(nat,iat,xyz) ! write new xopt.xyz
    call getgrad
    gr=grad
    xyz(j,i)=xyz(j,i)-2d0*step
    call newxyz(nat,iat,xyz)
    call getgrad
    gl=grad
    xyz(j,i)=xyz(j,i)+step
  do ii=1,nat
   do jj=1,3
    hj=(ii-1)*3+jj
    xx=(gr(jj,ii)-gl(jj,ii))/(2d0*step)
    if(abs(xx).gt.thrE)  hess(hi,hj)=xx
   enddo ! jj-loop
  enddo  ! ii-loop

 enddo ! j-loop
 enddo  ! i-loop

xyz=xyz0
  deallocate(gl,gr)

! symmetrize
  print*, 'Symmetrizing Hessian ..'
  do i=1,nat3
     do j=1,nat3
        chess(j,i)=(hess(i,j)+hess(j,i))*0.5d0
!        chess(i,j)=chess(j,i)
     enddo
  enddo

  print*, 'xopt cartesian Hessian  written into .. ''xopt.chess'' '
 call wrhess(nat3,chess,'xopt.chess')

deallocate(hess)

!  call cpu_time(t1)
!  call prtim(6,t1-t0,'t','hessian')

end subroutine 


subroutine newxyz(nat,iat,xyz)
use progs
implicit none
integer k,nat,iat(nat)
real(8) xyz(3,nat),f
character(2) esym
    f=0.5291770d0

 open(unit=55,file=xyzfile)
 write(55,'(I5)') nat
 write(55,'(2F16.8)') 0,0
  do k=1,nat
  write(55,'(a2,5x,3(F18.14,3x))') esym(iat(k)), xyz(1,k)*f,xyz(2,k)*f,xyz(3,k)*f
  enddo
  close(55)
end subroutine
