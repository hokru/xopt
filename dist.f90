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
subroutine get_dist(echo)
use parm

real(8) dx,dy,dz
logical echo

ki=0
kj=0
! loop over all unique pairs
k=1
do i=1,nat-1
 do j=i+1,nat
 dx=xyz(1,i)-xyz(1,j)
 dy=xyz(2,i)-xyz(2,j)
 dz=xyz(3,i)-xyz(3,j)
 dist(k)=sqrt(dx*dx+dy*dy+dz*dz)
 if(dist(k).gt.thrR) cycle
 kj(k)=j
 ki(k)=i
 k=k+1
 enddo
enddo
k=k-1

if(echo) then
write(*,'(2x,a,F6.2,a)') 'Keeping: ',(dble(k)/dble(npair)*100d0),' % of atom pairs after distance screening'
write(*,'(2x,a,I8)') 'npair before screening : ' ,npair
npair=k
write(*,'(2x,a,I8)') 'npair after  screening : ' ,npair
endif
end subroutine
