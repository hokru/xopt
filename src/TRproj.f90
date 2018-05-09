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
!    Foobar is distributed in the hope that it will be useful,                        *
!    but WITHOUT ANY WARRANTY; without even the implied warranty of                   *
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    *
!    GNU Lesser General Public License for more details.                              *
!                                                                                     *
!    You should have received a copy of the GNU Lesser General Public License         *
!    along with Xopt.  If not, see <http://www.gnu.org/licenses/>.                    *
!                                                                                     *
!**************************************************************************************
subroutine getcema(nat,xyzcm)
use fiso
implicit none
integer i,nat
real(r8) xx,yy,zz
real(r8), intent(inout) :: xyzcm(3,nat)


xx=0d0
yy=0d0
zz=0d0


do i=1,nat
 xx=xx+xyzcm(1,i)
 yy=yy+xyzcm(2,i)
 zz=zz+xyzcm(3,i)
enddo

xx=xx/nat
yy=yy/nat
zz=zz/nat

do i=1,nat
 xyzcm(1,i)=xyzcm(1,i)-xx
 xyzcm(2,i)=xyzcm(2,i)-yy
 xyzcm(3,i)=xyzcm(3,i)-zz
enddo
end subroutine


subroutine TRpr(nat,xyz,hess)
! calculating the translational-rotational projection matrix
! orthogonalize
! project out
use fiso, only:r8
implicit none
integer i,nat,nat3
real(r8) cema(3,nat),xyz(3,nat)
real(r8) TR(nat*3,6)
real(r8) hess(nat*3*(nat*3+1)/2)

integer lwork,info
real(r8) qwork,S(nat*3),U(nat*3,6),VT(6,6)
real(r8), allocatable :: work(:)

cema=xyz
call getcema(nat,cema)
TR=0_r8

nat3=nat*3
do i=1,nat
! translation
  TR(3*(i-1)+1,1) = 1.0_r8
  TR(3*(i-1)+2,2) = 1.0_r8
  TR(3*(i-1)+3,3) = 1.0_r8

! x,y,z rotation
  TR(3*(i-1)+1,4) = 0.0_r8
  TR(3*(i-1)+2,4) = -cema(3,i)
  TR(3*(i-1)+3,4) =  cema(2,i)

  TR(3*(i-1)+1,5) = cema(3,i)
  TR(3*(i-1)+2,5) = 0.0_r8
  TR(3*(i-1)+3,5) = -cema(1,i)

  TR(3*(i-1)+1,6) = -cema(2,i)
  TR(3*(i-1)+2,6) = cema(1,i)
  TR(3*(i-1)+3,6) = 0.0_r8
enddo

! orthogonalize
! SVD -> U
! query lwork, allocate, do it
call DGESVD('S','S', nat3, 6, TR , nat3, S, U, nat3, VT, 6, QWORK, -1, INFO )
lwork=int(qwork)
allocate(work(lwork))
call DGESVD('S','S', nat3, 6, TR , nat3, S, U, nat3, VT, 6, WORK, lwork, INFO )
deallocate(work)

if(info.ne.0) stop 'SVD orthogonalization failed'
! nearest ortho mat U*Vt
call matmult('N','N',nat3,6,6,u,vt,tr)
! call DGEMM('N','N',nat3,6,6,1.0d0,U,nat3,VT,6,0.0d0,TR,nat3)

! projection
call projMat(nat3,6,TR,hess)
end subroutine



subroutine projMat(n,m,pmat,fmat)
! apply orthogonal projection pmat to fmat
use fiso, only: r8
implicit none
integer i,j,k
!real(8) scrnm(n,m),scrnn(n,n)
integer, intent(in) :: n,m
real(r8), allocatable :: scrnm(:,:),scrnn(:,:)
real(r8), intent(in):: pmat(n,m)
real(r8), intent(inout) :: fmat(n*(n+1)/2)

allocate(scrnm(n,m),scrnn(n,n))

call packM(n,scrnn,fmat,'unpack') !fmat -> scrnn
call dsymm('l','u',n,m,1.0d0,scrnn,n,pmat,n,0.0d0,scrnm,n) !scrnnm=scrnn*pmat
call dgemm('n','t',n,n,m,1.0d0,scrnm,n,pmat,n,0.0d0,scrnn,n) !scrnn=scrnm*pmat'

! fmat=fmat-scrnn-scrnn'

do i=1,n
 do j=1,i
   k = i*(i-1)/2 + j
   fmat(k) = fmat(k) - scrnn(i,j) - scrnn(j,i)
 end do
end do

  call dgemm('t','n',n,m,n,1.0d0,scrnn,n,pmat,n,0.0d0,scrnm,n)
  call dgemm('n','t',n,n,m,1.0d0,pmat,n,scrnm,n,0.0d0,scrnn,n)

! fmat=fmat+scrnn
do i=1,n
 do j=1,i
   k = i*(i-1)/2 + j
   fmat(k) = fmat(k) + scrnn(i,j)
 end do
end do


deallocate(scrnm,scrnn)
end subroutine
