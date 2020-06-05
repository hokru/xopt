
subroutine cioptgrad
! simple conical section optimization using Levine/Martinez algorithm
! states  I<J, i.e. j=i+1, e.g. I=gs, J=s1 (opposite as in the paper)
! if roots flip, we simply interchange energies and gradients.
use parm
use logic
!use timing
use progs
use popt
use constant, only: au2ev
use fiso, only: r8, stdout
implicit none
character(255) aa,bb
character(255) pwd,dir1,dir2,exepath
real(r8) ei,ej,deij,meij,gij,fold,etmp,gtmp(3,nat)
real(r8) fij,term,ddot,dnrm2
! change these to allocation...
real(r8) gradij(3,nat),gipgj(3,nat),gimgj(3,nat),gi(3,nat),gj(3,nat),uvec(nat3),gvec(nat3),xvec(nat3)
logical flipped
real(r8) ida,ida2,sig_max ! to avoid division by zero

sig_max=50
fold=energy

if(tm) then
  call system('firm alles.cc')
  call system('firm restart.cc')
  call system('firm gradient')
  call system('firm energy')
endif


! for numerical gradient we can do everything in current directory!
if(numgrad) then
    if(nproc>1) then
      call get_pgrad(ei,ej,gi,gj)
    else
      call get_numgrad2(ei,ej,gi,gj)
    endif
else

  ! directory definitions. Fixed for now
  dir1='/stateI.xopt'
  dir2='/stateJ.xopt'
  !call getcwd(pwd)
  call get_environment_variable('PWD',pwd)

  ! change to dir1
  exepath=adjustl(trim(pwd))//trim(adjustl(dir1))
  call chdir(adjustl(trim(exepath)))

  !----------------------------------------
  ! It is an unfortunate design to have to repeat getgrad.
  ! Long term this should be done better
  ! for state I
  if(gei) then
    call system(trim(command_gei))
    call egradfile(nat,gi,ei,trim(geigrad))
  elseif (tm) then
      call system('actual -r &> /dev/null')
      aa='dscf > '//xjob
      bb='ricc2 >> '//xjob
      call system(aa)
      call system(bb)
      call tmgrad(nat,gi,ei)
  elseif(orca) then
      call system(trim(command_orca))
      call rdorcagrad(nat,gi,ei   ,'orca.engrad')
      call system("cp xopt.job xopt.last")
  elseif(gaus) then
      call IOgaussian('g.in')
      call system(trim(command_gaus))
      call g09grad(nat,gi,ei)
      call system("cp xopt.job xopt.last")
  elseif(amber) then
      aa='rm -f forcedump.dat force.dat > /dev/null'
      call system (trim(aa))
      call wamber(nat3,xyz,'amber.rst')
      call system(trim(command_amber))
      call ambergrad(nat,gi,ei,apbs)
      call cpfile('xopt.job','xopt.last','.')
  !    aa='sander -O -i '//trim(ambin)// &
  !             ' -o '//trim(xjob)// &
  !             ' -c amber.rst'// &
  !             ' -p '//trim(ambtop)
  !    call system(trim(aa))
  !    call ambergrad(nat,grad,energy,apbs)
  endif
  ! call getd3gcp(nat,ei,gi)
  !----------------------------------------

  ! change to dir2
  exepath=adjustl(trim(pwd))//trim(adjustl(dir2))
  call chdir(trim(exepath))

  !----------------------------------------
  ! for state J
  if(gei) then
    call system(trim(command_gei))
    call egradfile(nat,gj,ej,trim(geigrad))
  elseif (tm) then
      call system('actual -r &> /dev/null')
      aa='dscf > '//xjob
      bb='ricc2 >> '//xjob
      call system(aa)
      call system(bb)
      call tmgrad(nat,gj,ej)
  elseif(orca) then
      aa='orca30 orca.in > '//xjob
      call system(trim(aa))
      call rdorcagrad(nat,gj,ej,'orca.engrad')
      call system("cp xopt.job xopt.last")
  elseif(gaus) then
      call IOgaussian('g.in')
      aa='run-g09 g.in xopt.job'
      call system(trim(aa))
      call g09grad(nat,gj,ej)
      call system("cp xopt.job xopt.last")
  elseif(amber) then
      aa='rm -f forcedump.dat force.dat > /dev/null'
      call system (trim(aa))
      call wamber(nat3,xyz,'amber.rst')
      call system(trim(command_amber))
      call ambergrad(nat,gj,ej,apbs)
      call cpfile('xopt.job','xopt.last','.')
  !    aa='sander -O -i '//trim(ambin)// &
  !             ' -o '//trim(xjob)// &
  !             ' -c amber.rst'// &
  !             ' -p '//trim(ambtop)
  !    call system(trim(aa))
  !    call ambergrad(nat,grad,energy,apbs)
  endif
  ! call getd3gcp(nat,ej,gj)
  !----------------------------------------

  ! return to working dir
  call chdir(adjustl(trim(pwd)))
endif

flipped=.false.
if(ej<ei) then
flipped=.true.
!write(6,'(12x,a)') 'swapping I<->J states'
 etmp=ej
 ej=ei
 ei=etmp
 gtmp=gj
 gj=gi
 gi=gtmp
endif

! form intermediates
meij=(ej+ei)/2.0_r8
deij=ej-ei
ida=1.0_r8/(deij+alpha)
if(ida<1e-12_r8) ida=0.0_r8
!gij=(deij*deij)/(deij+alpha)
gij=(deij*deij)*ida
fij=meij+sigma*gij

do i=1,nat
 gipgj(1:3,i)=gj(1:3,i)+gi(1:3,i)
 gimgj(1:3,i)=gj(1:3,i)-gi(1:3,i)
enddo
gipgj=gipgj*0.5_r8

ida2=1.0_r8/((deij+alpha)**2)
if(ida2<1e-12) ida2=0
term=(deij**2+2.0_r8*alpha*deij)*ida2

! update sigma somewhat dynamically
if(abs(fold-fij)<5e-5.and.abs(deij)>1e-3) then
  write(stdout,'(12x,a,F7.3)') ' increasing sigma by',sigma*(2*deij/alpha)
  sigma=sigma+sigma*(2*deij/alpha)
  maxd=0.2
 if(sigma>sig_max) sigma=sig_max
endif


write(stdout,'(12x,a,F7.3,2x,a,F5.1,2(a,1x,F13.7),a,L1)') &
 'gap[eV]: ',deij*au2ev,' penalty: ', sigma,' E(low): ', ei, ' E(high): ', ej, ' root flip: ',flipped


gradij=gipgj+term*gimgj*sigma

! modified gradient along penalty direction (eq. 11,12)
gnorm=0.0_r8

k=0
do i=1,nat
  do j=1,3
    k=k+1
    gvec(k)=gradij(j,i)
    xvec(k)=term*gimgj(j,i) ! penalty part
  enddo
enddo

! eq 13
uvec=gvec/dnrm2(nat3,xvec,1)

! eq 11
gnorm=(1.0_r8/sigma)*ddot(nat3,gvec,1,uvec,1)
!~ print*,gnorm

! eq 12
xvec=gvec-ddot(nat3,gvec,1,uvec,1)*uvec

!write(6,'(3x,a,D10.3,2x,E10.3,x,L)') '| para gconv',gnorm,0.005,gnorm<0.005
!write(6,'(3x,a,E10.3,2x,E10.3,x,L)') '| perp gconv',dnrm2(nat3,xvec,1),0.005,dnrm2(nat3,xvec,1)<0.005

if(deij<1.e-3.and.sigma>=50) then
write(stdout,*) 'Gap small and max. penalty reached. Probably wont get better, likely worse.'
write(stdout,*) 'Gap:',deij,' penalty:',sigma
write(stdout,*) 'E(I):',ei,'E(J):',ej
write(stdout,*) 'stopping and writing final geometry'
call wrxyz(trim(xyzfile))
call error('ciopt request urgent stop')
endif

! transfer key variables back to optimizer
energy=fij
grad=gradij

end subroutine

!
! subroutine inc_sigma(sigma,sig_max)
! ! fixed, stepwise sigma update (not used)
! use fiso, only: r8, stdout
! implicit none
! real(r8) :: sigma,sig_max
!  write(stdout,*) '   --> increasing sigma'
! if(sigma<4)then
!  sigma=5
! elseif(sigma<6) then
!  sigma=8
! elseif(sigma<10) then
!  sigma=10
! else
!  sigma=sigma*2
! endif
! if(sigma>sig_max) sigma=sig_max
! end subroutine
