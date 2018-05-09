
!**********************************
!*       internal coord optimizer *
!**********************************
subroutine intopt()
use parm, only: nat,nvar,gint,chess,grad,gnorm,hint,nat3,energy,iat,xyz
use logic
use popt
use progs, only: xyzfile
use fiso, only: r8
implicit none
integer iter, info,c,i
real(r8)  displ(nvar),dint(nvar),oldint(nvar)
real(r8) f,dE,gmax,maxdipl,dnorm,oldE,oldG(nvar),oldD(nvar)
character(2) esym
real(r8) ddot,Lam
logical EE,GG,GM,DM,DOupdate
real(r8) DNRM2
real(r8) xyznew(3,nat)

logical upstep

EE=.false.
GG=.false.
GM=.false.

allocate(gint(nvar))

energy=0d0
oldG=0d0
displ=0d0
dint=0d0
c=0
DOupdate=.false.
upstep=.false.

oldd=0d0


! 'touch' log files
call touchlog

! main loop
write(6,'(2x,a)') ''
write(6,'(2x,a)') 'Starting internal coordinates geometry optimization'
write(6,'(2x,a)') 'Iter          Energy         dE           gnorm         gmax           dnorm       maxdispl       lambda'

do iter=iiter,maxiter

oldE=energy
oldG=gint
oldD=displ

! energy+gradient

call wrxyz(trim(xyzfile))
call getgrad
call loggrad(nat,iter,grad,energy)
call grad2int(grad,gint)

dE=energy-oldE


gnorm=dnrm2(nvar,gint,1)
gmax=maxval(gint)


! log files
call appxyz('xopt.log')


upstep=.true.

! update hessian, chess not used
call Hupdate(iter,chess,hint,nvar,gint,oldG,displ,0,DOupdate)
! call hint2xyz(nvar,nat3,hint,chess,b)
call wrhess(nat3,chess,'xopt.hess')
DOupdate=.true.

! before we make the QN step we need G+H in nonredundant space
call proj_nonred(gint,hint)


select case(iopt)
 case(21,31)
  call getLambdaRFO(displ,mode,Lam)
 case(22,32)
  call getLambdaSIRFO(displ,mode,Lam)
 case(23,33)
  call conjgrad(nvar,gint,oldg,oldD,displ)
  DOupdate=.false.
  lam=0d0
 case(24)
  call conjgrad(nvar,gint,oldg,oldD,displ)
  call getLambdaRFO(oldD,mode,Lam)
  displ=0.9*displ+oldD*0.1
end select


! check step components
do i=1,nvar
   if(abs(displ(i)).gt.maxd) then
      if(displ(i) < 0) displ(i)=-maxd
      if(displ(i) > 0) displ(i)= maxd
   endif
enddo

dnorm=dnrm2(nvar,displ,1)
maxdipl=maxval(abs(displ))

! back transformation
! linear step
! call get_Amatrix(Amat)
call dint2cart(dint,xyznew)
! new internals
call get_primitives(nat,iat,xyz)
call get_Bmatrix(nat,xyz)
call get_G_matrix()
! while ()
! 1.displace internals q(new)=q(old)+dint
! 2. make new xyz: xyz(new)=xyz0+dcart
!    2a call get_Amatrix(Amat)
!    2b call dint2cart(nat3,nints,dint,dcart,Amt)
!    2c make new internals q(trial) from xyz(new) (B+G matrix)
! 3. delta_q= qnew-q(trial)
!    3a new xyz: dint2cart with delta_q
! 4. check convergence 

! debug code
!do i=1,nvar
!  write(*,'(''ANC,Grad,Displ'',3F12.6)')ANC(i),gint(i),displ(i)
!enddo



EE=.false.
GG=.false.
GM=.false.
DM=.false.
if(abs(dE)<=econv) EE=.true.
if(gnorm<=gconv) GG=.true.
if(gmax<=maxgrad) GM=.true.
if(maxdipl<=dconv) DM=.true.

if(iter==1) then
 call pciter(iter,energy,0d0,gnorm,gmax,dnorm,maxdipl,Lam,EE,GG,GM,DM)
else
 call pciter(iter,energy,dE,gnorm,gmax,dnorm,maxdipl,Lam,EE,GG,GM,DM)
 if(dE==0.and..not.amber) stop 'ERROR: no energy change!'
endif


if(dE>=0d0) print*, 'Positive step detected!'

select case(iconv)
 case default
 case(4) ! all 4 thres reached (default)
  if(EE.and.GG.and.GM.and.DM) then
   call endopt(dE,gnorm,gmax,dnorm)
   return
  endif
 case(2) ! only energy and gradient
  if(EE.and.GG) then
   call endopt(dE,gnorm,gmax,dnorm)
   return
  endif
 case(1) ! only energy (eg for pre-optimizations)
  if(EE) then
   call endopt(dE,gnorm,gmax,dnorm)
   return
  endif
end select



!if(abs(dE)<=1e-5) then
!  xyz0=xyz
!  call gethess
!  call copt
!  return
!endif


enddo ! end main loop

print*,'FAILED!'
return

end subroutine 


