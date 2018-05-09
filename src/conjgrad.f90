subroutine conjgrad(xvar,g,og,od,d)
! conjugate gradient
! Polak-Ribiere beta value
use fiso, only: r8
implicit none
integer i,xvar
real(r8), intent(in) :: g(xvar),og(xvar),od(xvar)
real(r8), intent(out) :: d(xvar)
real(r8) dg(xvar),ogog,dgg,ddot,beta

dg=g-og
dgg=ddot(xvar,g,1,dg,1)
ogog=ddot(xvar,og,1,og,1)
beta=dgg/ogog
if(sum(og)<1e-10) beta=0.0_r8

do i=1,xvar
  d(i)=-g(i)+od(i)*beta
enddo
end subroutine conjgrad
