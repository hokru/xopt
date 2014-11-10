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
! conjugate gradient
! Polak-Ribiere beta value
! Taken from F Jensen "introduction to comp chem"
subroutine conjgrad(xvar,g,og,od,d)
implicit none
integer i,j,xvar
real(8), intent(in) :: g(xvar),og(xvar),od(xvar)
real(8), intent(out) :: d(xvar)
real(8) dg(xvar),ogog,dgg,ddot,beta

dg=g-og

dgg=ddot(xvar,g,1,dg,1)
ogog=ddot(xvar,g,1,og,1)
beta=dgg/ogog
do i=1,xvar
  d(i)=-g(i)+od(i)*beta
enddo

return
end subroutine conjgrad
