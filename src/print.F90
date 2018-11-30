subroutine head
use fiso, only: stdout, compiler_version
! prints main header data (name, version)
use progs, only: workdir
implicit none
character(80) host
call get_environment_variable('HOSTNAME',host)
call get_environment_variable('PWD',workdir)
write(stdout,'(a)') "              "
write(stdout,'(a)') "  -------------------------- "
write(stdout,'(a)') " |                         | "
write(stdout,'(a)') " |      \_/ _  _ |_        | "
write(stdout,'(a)') " |      / \(_)|_)|_        | "
write(stdout,'(a)') " |            |            | "
write(stdout,'(a)') "  ------------------------- "
write(stdout,'(a)') "              "
! write(stdout,'(a)') ' |-----------------------------|'
! write(stdout,'(a)') ' |  __   __               _    |'
! write(stdout,'(a)') ' |  \ \ / /              | |   |'
! write(stdout,'(a)') ' |   \ V /   ___   _ __  | |_  |'
! write(stdout,'(a)') " |   /   \  / _ \ | '_ \ | __| |"
! write(stdout,'(a)') ' |  / /^\ \| (_) || |_) || |_  |'
! write(stdout,'(a)') ' |  \/   \/ \___/ | .__/  \__| |'
! write(stdout,'(a)') ' |              | |            |'
! write(stdout,'(a)') ' |              |_|            |'
! write(stdout,'(a)') ' |-----------------------------|'
write(stdout,'(a)') '  An eXternal OPTimizer by'
write(stdout,'(a)') '   Holger Kruse (mail2holger@gmail.com)'
write(stdout,'(a)') '              '
write(stdout,'(a)')' Version:  2.0 rc1      '
write(stdout,'(a)')'             '
write(stdout,'(a)',advance="no") ' Date: '
call timestamp()
write(stdout,'(2a)') ' Host:   ',trim(host)
write(stdout,'(2a)') ' working directory: ',trim(workdir)
write(stdout,'(a)') ''
include '../version.dat'
write(stdout,*)'             '
write(stdout,'(2a)') ' compiler version: ',compiler_version()
write(stdout,'(a)') ''
write(stdout,'(a)') ''
! write(stdout,'(2a)') ' compiler options: ',compiler_options() ! ugly
call p_gcpd3
end subroutine head

subroutine p_gcpd3
use fiso,only: stdout
!use logic, only: d3,gcp

write(stdout,'(2a)',advance="no") ' internal dftd3 library : '
#ifdef DFTD3
write(stdout,'(a)') 'yes'
#else
write(stdout,'(a)') 'no'
#endif

write(stdout,'(2a)',advance="no") ' internal gcp library   : '
#ifdef GCP
write(stdout,'(a)') 'yes'
#else
write(stdout,'(a)') 'no'
#endif

write(stdout,'(2a)',advance="no") ' internal syva library  : '
#ifdef SYVA
write(stdout,'(a)') 'yes'
#else
write(stdout,'(a)') 'no'
#endif


write(stdout,'(a)')'             '
write(stdout,'(a)')'             '


end subroutine p_gcpd3


subroutine p_info
! prints essential data about molecule and algorithms
use fiso, only: stdout
use popt
use logic
use parm
use MDdat
implicit none
integer io
logical da
character(120) aa

if(do_output) then
 write(stdout,'(2a)') ' output: ',trim(output_name)
endif
io=stdout

write(io,'(a)')
write(io,'(2x,a)') ' Molecular data:'
write(io,'(3x,a,I6)') ' Atoms   = ',nat
call composition
call getCOM
call getIntertia
write(io,'(3x,a)') '   '
call molbox(nat,xyz,iat,cell,.true.,.true.)
write(io,'(3x,a)') '   '

!******************************             
!*         MOLECULAR DYNAMICS *
!******************************             
if(do_md) then
  write(io,'(3x,a)') '!----------------------------! '
  write(io,'(3x,a)') '!  molecular dynamics module !'
  write(io,'(3x,a)') '!----------------------------!'
  write(io,'(3x,a)') '   '
  write(io,'(2x,a)')   ' Integrator    :      --> velocity-verlet '
  write(io,'(2x,a)',advance="no") ' Thermostat    :'
  if(thermo_NH) then
    write(io,'(3x,a)') '   --> Nose-Hoover chains'
  elseif(thermo_scale) then
    write(io,'(3x,a)') '   --> velocity rescaling '
  elseif(thermo_berend) then
    write(io,'(3x,a)') '   --> Berendsen scaling'
  else
   write(io,'(3x,a)') '   --> none'
  endif
  write(io,'(3x,a)') '   '
  write(io,'(2x,a)') ' MD settings   :'
  write(io,'(3x,a,I8)') ' # steps            = ',maxstep
  write(io,'(3x,a,F7.4)') ' time step (fs)   = ',dt
  write(io,'(3x,a,F7.2)') ' temperature (K)  = ',Temp0
  write(io,'(3x,a,I7)')   ' particles        = ',nat
  write(io,'(3x,a)') '   '
  write(io,'(2x,a)') ' print settings:'
  write(io,'(3x,a,I8)') ' file logging    = ',ifile
  write(io,'(3x,a,I8)') ' screen printing = ',iprint

!******************************             
!*          QUASI-NEWTON      *
!******************************             
else
  write(io,'(3x,a)') '!---------------------------------! '
  write(io,'(3x,a)') '!  quasi-Newton Raphson module    !'
  write(io,'(3x,a)') '!---------------------------------!'
  write(io,'(3x,a)') '   '
  write(io,'(2x,a)',advance="no")      ' run type         : '
  if(mode+tsmode>1)  write(io,'(3x,a)') '      --> TS search'
  if(mode+tsmode==1) write(io,'(3x,a)') '      --> minimization'
  write(io,'(a)')

! general
  write(io,'(2x,a)') ' Convergence criteria:'
  write(io,'(3x,a,ES10.3)') ' energy    = ',econv
  write(io,'(3x,a,ES10.3)') ' gnorm     = ', gconv
  write(io,'(3x,a,ES10.3)') ' max grad  = ',maxgrad
  write(io,'(3x,a,ES10.3)') ' max displ = ',dconv
  write(io,'(a)')
  write(io,'(2x,a)') ' Setup: '
  write(io,'(3x,a,F7.4)') ' max displacement = ', maxd
  write(io,'(3x,a,I4)')  ' max iter         = ',maxiter
  if(tsopt) write(io,'(3x,a,I4)')   ' EF mode          = ', tsmode
  write(io,'(a)')

  ! selected methods
  write(io,'(2x,a)',advance="no")      ' selected program : '
  if(tm) write(io,'(5x,a,L1)')   '    --> Turbomole'
  if(orca) write(io,'(5x,a,L1)') '    --> ORCA'
  if(mopac) write(io,'(5x,a,L1)')'    --> mopac12'
  if(gaus) write(io,'(5x,a,L1)') '    --> Gaussian09'
  if(psi4) write(io,'(5x,a,L1)') '    --> PSI4'
  if(openbabel) write(io,'(5x,a,L1  )') '    --> openbabel'
  if(GEI) write(io,'(5x,a)')   '     --> external interface'
  write(io,'(a)') ''


  write(io,'(2x,a,I2)',advance="no")   ' Hessian update   : ',iupdate
  if(iupdate==1)   write(io,'(3x,a)') '    --> SR1-BFGS'
  if(iupdate==2)   write(io,'(3x,a)') '    --> SR1'
  if(iupdate==3)   write(io,'(3x,a)') '    --> BFGS'
  if(iupdate==4)   write(io,'(3x,a)') '    --> SS-BFGS'
  if(iupdate==5)   write(io,'(3x,a)') '    --> SS-SR1-BFGS'
  if(iupdate==6)   write(io,'(3x,a)') '    --> SR1-PSB'
  if(iupdate==7)   write(io,'(3x,a)') '    --> PSB'


  write(io,'(2x,a)',advance="no")   ' initial Hessian  : '
  if(readHess) then
                              write(io,'(3x,a)') '      --> user provided'
  else
    if(hmodel=='numerical')   write(io,'(3x,a)') '      --> numerical'
    if(hmodel=='lindh')       write(io,'(3x,a)') '      --> H(int) Lindh '
    if(hmodel=='simple')      write(io,'(3x,a)') '      --> H(int) simple '
    if(hmodel=='alm')         write(io,'(3x,a)') '      --> H(int) Fischer/Almloef'
    if(hmodel=='hdiag')       write(io,'(3x,a)') '      --> H(cart) unit'
    if(hmodel=='mass')        write(io,'(3x,a)') '      --> H(cart) mass-weighted'
  endif

  write(io,'(2x,a,I2)',advance="no")   ' convergence mode : ',iconv
  if(iconv==4) write(io,'(3x,a)') '    --> all criteria'
  if(iconv==2) write(io,'(3x,a)') '    --> energy + gradient norm'
  if(iconv==1) write(io,'(3x,a)') '    --> energy only'

  write(io,'(2x,a,I3)',advance="no")   ' optimizer mode   : ',iopt
  if(iopt==11) write(io,'(3x,a)') '   --> cartesian + RFO'
  if(iopt==12) write(io,'(3x,a)') '   --> cartesian + SI-RFO'
  if(iopt==21) write(io,'(3x,a)') '   --> approx. normal coords + RFO'
  if(iopt==22) write(io,'(3x,a)') '   --> approx. normal coords + SI-RFO'
  if(iopt==23) write(io,'(3x,a)') '   --> approx. normal coords + CG'
  if(iopt==24) write(io,'(3x,a)') '   --> approx. normal coords + CG/RFO mix'
  if(iopt==25) write(io,'(3x,a)') '   --> approx. normal coords + P-RFO'
  write(io,'(a)')
endif

if(ciopt) then
  write(io,'(a)')
  write(io,'(1x,a)')      ' Requested:  smoothed penalty function CI gradient'
  write(io,'(3x,a,F8.4)') '   smoothing factor alpha: ', alpha
  write(io,'(3x,a,F8.4)') '   penalty factor sigma  : ', sigma
  write(io,'(a)')
endif


if(ppot.or.d3.or.gcp) then
write(io,'(2x,a)')   ' helper & options  : '
if(gcp.and..not.exlib) write(io,'(5x,2a)')   '   --> (external) gcp ', trim(gcplevel)
if(gcp.and.exlib)      write(io,'(5x,2a)')   '   --> (internal) gcp ', trim(gcplevel)
if(d3.and..not.exlib)  write(io,'(5x,2a)')   '   --> (external) dftd3 ', trim(d3func)
if(d3.and.exlib)       write(io,'(5x,2a)')   '   --> (internal) dftd3 ', trim(d3func)
if(ppot) write(io,'(5x,a)')   '   --> ppot '
if(mopac) then
 inquire(file='SETUP',exist=da)
 if(da) then
  open(88,file='SETUP')
  read(88,'(a)') aa
  write(io,*) 'SETUP: ',trim(aa)
 else
  call error(6,' no mopac SETUP file found!!')
endif

if(large) then
  write(io,'(2x,a)')   ' * LARGE mode activated * '
endif

endif
write(io,'(a)')
endif

end subroutine p_info

subroutine pchead()
use fiso, only: stdout
implicit none
write(stdout,'(2x,a)') 'Iter          Energy         dE           gnorm         gmax      dnorm     maxdispl   lambda     dE(pred)'
end subroutine pchead

subroutine pciter(iter,e,de,gnorm,gmax,dnorm,disp,lam,EE,GG,GM,DM,prede)
! prints values for the optimization step
use fiso, only: stdout,r8
implicit none
real(r8), intent(in) :: e,de,disp,gmax,lam,dnorm,gnorm,prede
integer, intent(in) :: iter
logical, intent(in) :: EE,GG,GM,DM

write(stdout,'(2x,i4,F18.8,2x,ES10.3,1x,L1,2x,F8.5,1x,L1,2x,ES10.3,1x,L1,2x,F7.3,2x,F9.4,1x,L1,2x,ES9.2,2x,ES9.2)') &
             iter,e,      de       ,EE,   gnorm,  GG,   gmax,     GM,    dnorm,disp,DM,lam,prede
end subroutine pciter


subroutine endopt(dE,gnorm,gmax,dnorm)
! prints final convergence criteria
use fiso, only: stdout,r8
use popt
implicit none
real(r8), intent(in) :: dE,gnorm,gmax,dnorm
write(stdout,*) 'CONVERGED !'
write(stdout,'(3x,a)') '               criteria   actual value'
write(stdout,'(3x,a,2(ES10.3,2x))') ' energy    = ',econv,dE
write(stdout,'(3x,a,2(ES10.3,2x))') ' gnorm     = ',gconv,gnorm
write(stdout,'(3x,a,2(ES10.3,2x))') ' max grad  = ',maxgrad,gmax
write(stdout,'(3x,a,2(ES10.3,2x))') ' max displ = ',dconv,dnorm
end subroutine


subroutine prtime(io)
! prints elapsed time
use timing, only: elapsed
use fiso,only: r8,stdout
implicit none
real(r8) s
integer d,h,m
integer io

call timedecomp(elapsed,s,m,h,d)

write(stdout,'("wall-time: ",1x,i3," d",1x,i3," h",1x,i3," m",f5.1," s")') d,h,m,s
end subroutine prtime

subroutine timedecomp(tin,s,m,h,d)
! input the seconds and get out the decomposition in s,m,h,d
use fiso, only: r8
real(r8), intent(in):: tin
real(r8),intent(out):: s
integer, intent(out):: d,h,m
real(r8) :: t
real(r8), parameter :: s2d = 86400_r8
real(r8), parameter :: s2h = 3600_r8
real(r8), parameter :: s2m = 60_r8

t=tin
d=int(t/s2d)!days
t=t-d*s2d
h=int(t/s2h) !hours
t=t-h*s2h
m=int(t/s2m)  !mins
t=t-m*s2m
s=t          !secs

end subroutine


subroutine progtime_start()
! obtains initial time
use timing
implicit none
integer t

call system_clock(t) !
tstart=dble(t)
end subroutine


subroutine progtime_end()
! obtains end time and provides elapsed time
use timing
implicit none
integer t
call system_clock(t) !
elapsed=(dble(t)-tstart)/dble(cr)
end subroutine


! real(8) function stime()
! ! function to obtain start clock time (unused)
! use timing
! implicit none
! integer t
! call system_clock(t) !
! stime=dble(t)
! end function
!
! real(8) function mytime(t1)
! ! unused timing function
! use timing
! implicit none
! integer t
! real(8) t1,t2
!
! call system_clock(t) !
! t2=dble(t)
! mytime=(t2-t1)/dble(cr)
!
! end function
!


subroutine debug1(string,time)
! timing function
use fiso, only: stdout,r8,i8
implicit none
real(r8),intent(out):: time ! current clock time
character(*),intent(in):: string ! custom name for the timing
integer(i8) t

call system_clock(t) !
time=dble(t)
write(stdout,'(a)',advance="no") string
end subroutine

subroutine status1(time)
! timing function
use fiso, only: stdout,r8,i8
implicit none
real(r8),intent(out):: time 
integer(i8) t

call system_clock(t) !
time=dble(t)
end subroutine



subroutine debug2(time)
! timing function
use fiso, only: stdout,r8,i8
use timing
implicit none
real(r8),intent(in):: time 
real(r8)t2,time2
integer(i8) t

call system_clock(t) !
t2=dble(t)
time2=(t2-time)/dble(cr)
if(time2<500_r8) write(stdout,'(a,F12.4,a)') '  (',time2,' s)'
if(time2>500_r8) write(stdout,'(a,F12.4,a)') '  (',time2/60_r8,' m)'
!if(time<500) write(stdout,'(a,F12.4,a)') ' done in: ',time,' s'
!if(time>500) write(stdout,'(a,F12.4,a)') ' done in: ',time/60_r8,' m'
end subroutine

subroutine status2(time)
! timing function
use fiso, only: stdout,r8
use timing
implicit none
real(r8),intent(in):: time 
real(r8)t2,time2
integer(8) t

call system_clock(t) !
t2=dble(t)
time2=(t2-time)/dble(cr)
if(time2<500_r8) write(stdout,'(2x,a,F12.4,a)') '--> (',time2,' s)'
if(time2>500_r8) write(stdout,'(2x,a,F12.4,a)') '--> (',time2/60_r8,' m)'
write(stdout,'(a)') ' '
end subroutine

subroutine message(str)
! print text line to stdout
use fiso, only: stdout
character(*), intent(in) :: str
 write(stdout,'(2x,a)') trim(str)
end subroutine

subroutine message_head(str)
! print text line to stdout
use fiso, only: stdout
character(*), intent(in) :: str
 write(stdout,'(2x,a)') trim(str)
end subroutine





subroutine timestamp ( )
!*****************************************************************************80
! original by John Burkardt, modified by H.Kruse
!  Licensing:
!  This code is distributed under the GNU LGPL license.
use fiso, only: stdout
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)
  write (stdout, '(2x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,1x)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, ' '
  return

  ! no AM/PM :-)
  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write (stdout, '(2x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, ' ', trim ( ampm )
end



subroutine p_memory(isize,title)
! prints memory of a word input with title.
! max 15 characters for pretty printing
! isize = word
! title = title/name
use fiso, only: stdout,r8
implicit none
real(r8), intent(in)  :: isize
character(*), intent(in) :: title
real(r8), parameter :: f2=8.0_r8/(1024.0_r8**2)
real(r8), parameter :: f3=8.0_r8/(1024.0_r8**3)
real(r8), parameter :: f4=8.0_r8/(1024.0_r8**4)
real(r8) :: x
integer :: l
integer, parameter:: maxl=15
character(40) :: string
character(maxl) :: title2

l=len(title)
 string='(3x,a,'': '',F9.3,a)'
if(l<=maxl) then
  title2(1:l)=title(1:l)
  title2(l+1:maxl)=' '
else
  title2(1:maxl)=title(1:maxl)
endif

x=isize*f2
 if(x<=1024) then
   write(stdout,string) title2,isize*f2,' Mib'
 elseif(x>1024.and.x<=1e6) then                          
   write(stdout,string) title2,isize*f3,' Gib'
 elseif(x>=1e6) then                                     
   write(stdout,string) title2,isize*f4,' Tib'
 endif
end 

subroutine printvib(nat3,e)
! subroutine printvib(nat3,e,max)
! print frequency informations
! maximal nmax (optional, default=18)
use fiso, only: r8,stdout
implicit none
integer, intent(in)  :: nat3
real(r8), intent(in) :: e(nat3)
integer              :: i,ii,iii
integer              :: n,loop,nmax
real(r8)             :: freqcm

! why is this not working?
! integer, optional :: max
! if(present(max)) then
!  nmax=max
! else
!  nmax=18
! endif
! print*,max,present(max)

nmax=18
loop=0
n=nat3
if(n<=6) then
   do i=1,n
      write(stdout,'(4x,I2,'':'',1x,F8.2)') i,freqcm(e(i))
   enddo
elseif(n<=12) then
   if(mod(n,2)==0) loop=n/2
   if(mod(n,2)==1) loop=(n-1)/2
   do i=1,loop
      ii=i+loop
      write(stdout,'(2x,2(2x,I2,'':'',1x,F8.2))') i,freqcm(e(i)),ii,freqcm(e(ii))
   enddo
   if(mod(n,2)==1)  write(stdout,'(14x,4x,I2,'':'',1x,F8.2)') n,freqcm(e(n))
else
   if(nat3>nmax) n=nmax
   if(mod(n,3)==0) loop=n/3
   if(mod(n,3)==1) loop=(n-1)/3
   do i=1,loop
      ii=i+loop
      iii=ii+loop
      write(stdout,'(2x,3(2x,I2,'':'',1x,F8.2))') i,freqcm(e(i)),ii,freqcm(e(ii)),iii,freqcm(e(iii))
   enddo
   if(mod(nat3,3)==1)  write(stdout,'(14x,4x,I2,'':'',1x,F8.2)') n,freqcm(e(n))
endif

end subroutine
