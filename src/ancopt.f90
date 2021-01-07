!**********************************
!* approx. normal coord optimizer *
!**********************************
subroutine ancopt()
use fiso, only: r8, stdout
use parm
use logic
use popt
use progs
implicit none
integer :: iter,c
integer :: nvar1,iflag
real(r8), allocatable :: displ(:),anc(:),oldanc(:)
real(r8) :: dE,gmax,maxdipl,dnorm,oldE,oldG(nvar),oldD(nvar)
real(r8) :: Lam,oldL
logical :: EE,GG,GM,DM,DOupdate,remake,downscale
real(r8) :: DNRM2,echange
real(r8) :: hold(nvar*(nvar+1)/2)
real(r8) :: time
! DIIS
real(r8) :: ogint(nvar,idiis),oanc(nvar,idiis)
integer :: ind
! TS
real(r8), allocatable :: hmode(:)

logical upstep

! do_gdiis=.false.
! do_gdiis=.true.

EE=.false.
GG=.false.
GM=.false.
remake=.false.

nvar1=nvar+1
allocate(gint(nvar),anc(nvar))
allocate(displ(nvar),oldanc(nvar))
if(tsopt) allocate(hmode(nvar))

energy=0.0_r8
oldG=0.0_r8
displ=0.0_r8
anc=0.0_r8
c=0
DOupdate=.false.
upstep=.false.

ogint=0.0_r8
oanc=0.0_r8
ind=0
hold=0.0_r8 ! hess damping
oldd=0.0_r8
lam=0.0_r8
oldL=0.0_r8

!if(ciopt) remake=.true.

!if(amber) upstep=.true.

! 'touch' log files
call touchlog

! main loop
write(stdout,'(2x,a)') ''
write(stdout,'(2x,a)') 'Starting approx. normal coordinates geometry optimization'
call pchead()

do iter=iiter,maxiter
call check_stop()

!if(ciopt.and.gnorm<1e-3.and..not.remake) remake=.true.
if(remake) then
call wrhess(nat3,chess,'xopt.hess')
  call int2xyz(nvar,anc,xyz)
  deallocate(b,hint)
  xyz0=xyz
  call getanc
 if(freeze)  call freezeBmat
  anc=0
  oldG=0.0_r8
  displ=0.0_r8
  c=0
 DOupdate=.false.
 remake=.false.
endif

!  upstep=.false.

oldE=energy
oldG=gint
oldD=displ
oldL=lam

if (debug) write(stdout,'(a)',advance="no") 'int2xyz ...'
! new coordinates
call int2xyz(nvar,anc,xyz)
if (debug) write(stdout,'(a)') '  done'


! new xyz file
call wrxyz(trim(xyzfile))
! energy+gradient
if(debug) call debug1('gradient',time)
call getgrad
if(debug) call debug2(time)
if (debug) write(stdout,'(a)',advance="no") 'restrain/constrain ...'
if(restrain)  call addRestrainGrad
if (freeze) call freezeGrad
if (debug) write(stdout,'(a)') '  done'

dE=energy-oldE


if (debug) write(stdout,'(a)',advance="no") 'loggrad ...'
call loggrad(nat,iter,grad,energy)
if (debug) write(stdout,'(a)') '  done'

if(debug) call debug1('gxyz2int',time)
call gxyz2gint
if(debug) call debug2(time)

if(.not.ciopt) gnorm=dnrm2(nvar,gint,1)
gmax=maxval(gint)


! log files
if (debug) write(stdout,'(a)',advance="no") 'appxyz ...'
call appxyz('xopt.log')
if (debug) write(stdout,'(a)') '  done'



if(debug) call debug1('Hessian update',time)
! update hint hessian
call Hupdate(iter,chess,hint,nvar,gint,oldG,displ,0,DOupdate)
! make new cart. chess
call hint2xyz(nvar,nat3,hint,chess,b)
DOupdate=.true.
if(debug) call debug2(time)


! call hessdamp(nvar,hint,hold,gnorm)

! GDIIS for interpolated gradient
! something is wrong
! if(do_gdiis) then
!     ind=ind+1
!     if(ind>5) ind=1
!     ogint(1:nvar,ind)=gint(1:nvar)
!     if(iter>=idiis.and.ind==5) then
!       call gdiis(idiis,iter,ogint,oanc,displ,gint)
!       gnorm=dnrm2(nvar,gint,1)
!       gmax=maxval(gint)
!     endif
! endif

if(tsopt) then
  if (debug) write(stdout,'(a)',advance="no") 'mode following ...'
  iflag=1
  if(iter==1) then
    iflag=0
    write(stdout,'(2x,''| EF optimization: Following mode'',i2)') tsmode
    call eigovl_int(hint,nvar,hmode,tsmode,iflag)
  else
    call eigovl_int(hint,nvar,hmode,tsmode,iflag)
  endif
  if (debug) write(stdout,'(a)') '  done'
endif

if(debug) call debug1('step estimate (RFO/CG)',time)
select case(iopt)
  case(21,31)
    if(tsopt) mode=tsmode+1 ! +1 because of augmented Hessian
    if(nat>2000.or.do_float) then
      call FastgetLambdaRFO(displ,mode,Lam) ! "no" gain
    else
      call getLambdaRFO(displ,mode,Lam)
    endif
 case(22,32)
    if(tsopt) mode=tsmode+1 ! +1 because of augmented Hessian
    call getLambdaSIRFO(displ,mode,Lam)
 case(23,33) ! experimental
    ! reset every 10 iters
    if(iter==1.or.mod(iter,10)==0) then
    displ=-gint
    else
      call conjgrad(nvar,gint,oldg,oldD,displ)
    endif
      DOupdate=.false.
      lam=0.0_r8
 case(24) ! experimental
    call conjgrad(nvar,gint,oldg,oldD,displ)
    call getLambdaRFO(oldD,mode,Lam)
  displ=0.9_r8*displ+oldD*0.1_r8
 case(25)
    if(tsopt) mode=tsmode
    call getLambdaPRFO(nvar,hint,gint,displ,lam,mode)
 case(26)
! idiis=iter
! call gdiis2(idiis,iter,nat3)
 case(99)
!  call splitRFO(displ,mode,Lam)
 case default
  call error('this should not happen :-)')
end select
if(debug) call debug2(time)

! check step
do i=1,nvar
   if(displ(i) > maxd) displ(i)= maxd
   if(displ(i) < -maxd) displ(i)= -maxd
enddo

downscale=.false.
dnorm=dnrm2(nvar,displ,1)
if(iter>1.and..not.tsopt) then
  if((oldL/lam)<0.1) then
  downscale=.true.
  displ=displ*0.5_r8
  endif
endif
dnorm=dnrm2(nvar,displ,1)
maxdipl=maxval(abs(displ))

oldanc=anc
anc=anc+displ

hold=hint

! open(unit=921,file="anc.debug",position='append')
! do i=1,nvar
!   write(921,'(2(F15.7,2x))') displ(i),gint(i)
! enddo
! write(921,*) ''
! close(921)

EE=.false.
GG=.false.
GM=.false.
DM=.false.
if(abs(dE)<=econv) EE=.true.
if(gnorm<=gconv) GG=.true.
if(gmax<=maxgrad) GM=.true.
if(maxdipl<=dconv) DM=.true.

if(iter==1) then
  dE=0.0_r8
  echange=0.0_r8
endif

call pciter(iter,energy,dE,gnorm,gmax,dnorm,maxdipl,Lam,EE,GG,GM,DM,echange)
if((abs(dE)<1e-12_r8).and.(.not.amber).and.(.not.iter==1)) call error('ERROR: no energy change!')

! predicted energy change
call Echange2_int(nvar,gint,hint,displ,echange)

if(downscale) then
 write(stdout,*) ' unreliable step detected! scaling down displacements'
endif

select case(iconv)
  case default
  case(4) ! all 4 thres reached (default)
    if(EE.and.GG.and.GM.and.DM) then
      goto 421 ! success
    endif
  case(2) ! only energy and gradient
    if(EE.and.GG) then
      goto 421 ! success
    endif
  case(1) ! only energy
    if(EE) then
      goto 421 ! success
    endif
end select


call wrestart(iter)

enddo ! end main loop


! final hessian in ascii (TM format)
! Takes a lot of disk space and not really needed?
!call wrhess(nat3,chess,'xopt.hess')

if(restrain) call printFinalRestrain
write(stdout,'(a)') 'FAILED!'
return ! <<<<< failed opt + exit >>>>>>>

!*************
!* success ! *
!*************
421 continue
call endopt(dE,gnorm,gmax,dnorm)
if(restrain) call printFinalRestrain
call wrestart(iter)
return

end subroutine ancopt
