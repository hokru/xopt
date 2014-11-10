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
subroutine int_bonds_B
use parm
use popt
real(8) e(3)
integer aa,bb
allocate(dist(npair),stat=istat)
allocate(kj(npair*(npair+1)/2),stat=istat)
allocate(ki(npair*(npair+1)/2),stat=istat)
if(stat.ne.0) call error(6,'internal_bonds allocation')

call get_dist(.true.)
print*,npair
nvar=npair
allocate(b(nat3,npair))
allocate(hint(nvar*(nvar+1)/2))

!make B matrix

!do i=1,nat-1
! do j=i+1,nat

!do i=1,nat-1
! do j=i+1,nat
!k=k+1
!  call evec(xyz(1,i),xyz(1,j),e)
!b(i,k)=-e(1)
!b(i+1,k)=-e(1)
!b(i+2,k)=-e(1)
! enddo
!enddo
!1 4
!7 10

do k=1,npair
! print*,k,ki(k),kj(k)
 aa=3*(ki(k)-1)
 bb=3*(kj(k)-1)
 call evec(xyz(1,ki(k)),xyz(1,kj(k)),e)
! print*,k,aa,bb
 b(aa+1,k)= -e(1)
 b(aa+2,k)= -e(2)
 b(aa+3,k)= -e(3)
 b(bb+1,k)= e(1)
 b(bb+2,k)= e(2)
 b(bb+3,k)= e(3)
enddo
!print*,B(1,1:6)
!print*,' '
print*,B(1:6,1)
!stop
call hxyz2int(nvar,nat3,hint,chess,b)

end subroutine

! get unit vector
subroutine evec(a,b,e)
implicit none
real(8), intent(in) :: a(3),b(3)
real(8), intent(out) :: e(3)
real(8) rab,dbond

rab=dbond(a,b);
e(1)=(a(1)-b(1))/rab
e(2)=(a(2)-b(2))/rab
e(3)=(a(3)-b(3))/rab

end subroutine

! distance between cartesian vectors a and b
real(8) pure function dbond(a,b)
implicit none
real(8), intent(in) :: a(3),b(3)
integer i
real(8) ab(3)
ab=a-b
dbond=sqrt( ab(1)**2 + ab(2)**2 + ab(3)**2 )
end function

! Hint=Bt*chess*B
!nvar,nat3 * nat3,nat3  * nat3,nvar
! M,K     * K*N
subroutine hxyz2int(nvar,nat3,hint,chess,b)
implicit none
real*8 hint(nvar*(nvar+1)/2),b(nat3,nvar)
real*8 chess(nat3,nat3)
integer nvar,nat3
integer i,j,k,padr
real*8 :: bh(nvar,nat3),xx
!real*8 :: bh(nat3,nvar),xx
real(8) h(nvar,nvar)
!                   M    N    K
call dgemm('T','N',nvar,nat3,nat3,1d0,b,nat3,chess,nat3,0d0,bh,nvar)  
call dgemm('N','N',nvar,nvar,nat3,1d0,bh,nvar,b,nat3,0d0,h,nvar)
call packM(nvar,h,hint,'pack')
end subroutine


