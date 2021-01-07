subroutine Hmass(h,string)
! take full Hessian h on input, mass-weight it, and output full Hessian again
! project trans/rot if string=project
use fiso, only: r8,stdout
use parm, only: iat,xyz,nat,nat3,i,j,k
use atomdata, only: ams
implicit none
integer ii
real(r8), allocatable :: hv(:),isqm(:)
real(r8), intent(inout):: h(nat3,nat3)
character(*), intent(in) :: string

allocate(hv(nat3*(nat3+1)/2),isqm(nat3))

do i = 1, nat
  do j = 1, 3
   ii = (i-1)*3+j
   isqm(ii)=1./sqrt(ams(iat(i)))
 enddo
enddo

call packM(nat3,h,hv,'pack')

if(string=='project') then
  call TRpr(nat,xyz,hv)
  write(stdout,'(2x,a)') 'T/R projecting hessian in mass-weight routine'
endif

write(stdout,'(2x,a)') 'mass-weighting Hessian'
!unpack + mass-weight
k=0
do i=1,nat3
    do j=1,i
       k=k+1
       h(j,i)=hv(k)*isqm(i)*isqm(j)
       h(i,j)=h(j,i)
    enddo
enddo

deallocate(hv,isqm)
end subroutine Hmass
