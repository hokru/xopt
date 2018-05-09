!program test
!implicit none
!real(8) r,mean,s
!real(8), allocatable:: data(:)
!integer(8) i,n
!call init_random_seed()
!
!n=1e7
!allocate(data(n))
!do i=1,n
!call randGaus(r)
!data(i)=r
!enddo
!
!mean=sum(data)/n
!print*, 'mean', mean
!
!s=0
!do i=1,n
!s=s+(data(i)-mean)**2
!enddo
!s=sqrt(s/n)
!print*, 's',s
!end


! normal distributed Random numbers
! polar form of Box-Mueller transformation
subroutine randGaus(y1)
use fiso, only: r8
implicit none
real(r8), intent(out) :: y1
real(r8) x1,x2,y2,w,r(2)

w=99
do while(w>=1.0_r8)
 call random_number(r)
 x1=2.0_r8*r(1)-1.0_r8
 x2=2.0_r8*r(2)-1.0_r8
 w=x1**2+x2**2
enddo
 w=sqrt( (-2.0_r8*log(w)) / w )
 y1=x1*w
 y2=x2*w ! not used
end



! GCC RANDOM_SEED EXAMPLE FROM DOCUMENATION
subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid,getpid
            integer(int64) :: t

            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
end subroutine init_random_seed
