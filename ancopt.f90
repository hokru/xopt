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

!**********************************
!* approx. normal coord optimizer *
!**********************************
subroutine ancopt(iiter)
use parm
use logic
use popt
use progs
implicit none
integer iter, info,c,iiter
integer nvar1,lwork
real(8), allocatable :: displ(:),anc(:),oldanc(:)
!real(8), allocatable :: Uaug(:,:), eaug(:),aux(:)
real(8) f,dE,gmax,maxdipl,dnorm,oldE,oldG(nvar),oldD(nvar)
character(2) esym
real(8) ddot,Lam
logical EE,GG,GM,DM,DOupdate
real(8) DNRM2
real(8) hold(nvar*(nvar+1)/2)

! DIIS
real(8) ogint(nvar,idiis),oanc(nvar,idiis),gdispl(1:nvar)
integer ind
logical ok

logical upstep

EE=.false.
GG=.false.
GM=.false.


nvar1=nvar+1
allocate(gint(nvar),anc(nvar))
allocate(displ(nvar),oldanc(nvar))

!lwork = 1 + 6*nvar1 + 2*nvar1**2
!allocate(Uaug(nvar1,nvar1),eaug(nvar1),aux(lwork))

energy=0d0
oldG=0d0
displ=0d0
anc=0
c=0
DOupdate=.false.
upstep=.false.

ogint=0
oanc=0
ind=0
hold=0d0 ! hess damping
oldd=0d0

!if(amber) upstep=.true.

! 'touch' log files
call touchlog

! main loop
write(6,'(2x,a)') ''
write(6,'(2x,a)') 'Starting approx. normal coordinates geometry optimization'
write(6,'(2x,a)') 'Iter          Energy         dE           gnorm         gmax           dnorm       maxdispl       lambda'

do iter=iiter,maxiter


! every X iterations make new normal coordinates
goto 321 !skip for now
c=c+1
if(iter==5.and..not.restrain.and..not.restart) then
call wrhess(nat3,chess,'xopt.hess')
  call int2xyz(nvar,anc,xyz)
  deallocate(b,hint)
  xyz0=xyz
  call getanc
 if(freeze)  call freezeBmat
  anc=0
  oldG=0d0
  displ=0d0
  c=0
 DOupdate=.false.
endif
321 continue

!  call int2xyz(nvar,oldanc,xyz)
!  upstep=.false.

oldE=energy
oldG=gint
oldD=displ

! new coordinates
call int2xyz(nvar,anc,xyz)


! new xyz file
call wrxyz(trim(xyzfile))
! energy+gradient
call getgrad
if(restrain)  call addRestrainGrad
if (freeze) call freezeGrad

dE=energy-oldE

!if(dE>0d0.and..not.upstep) then
!  print*, 'Positive step detected! Doing NR-step'
!  call int2xyz(nvar,oldanc,xyz)
!  call wrxyz(trim(xyzfile))
!  gint=oldG
!  energy=oldE
!  call getLambdaPRFO(nvar,hint,gint,displ,lam,-1)
!  print*,displ
!  upstep=.true.
!  goto 42
!endif
!upstep=.false.


call loggrad(nat,iter,grad,energy)
call gxyz2gint

gnorm=dnrm2(nvar,gint,1)
gmax=maxval(gint)


! log files
call appxyz('xopt.log')



 ! update hessian, chess not used
 call Hupdate(iter,chess,hint,nvar,gint,oldG,displ,0,DOupdate)
 !print*,'hint2xyz'
 call hint2xyz(nvar,nat3,hint,chess,b)
 !call oldhint2xyz(nvar,nat3,hint,chess,b)
 call wrhess(nat3,chess,'xopt.hess')
 DOupdate=.true.

! call hessdamp(nvar,hint,hold,gnorm)

! GDIIS, "ind" goes from 1-idiis, not working proberly
!ind=ind+1
!if(ind>idiis) ind=1
!if(ind>=1) ogint(1:nvar,ind)=gint(1:nvar)
!if(ind>=1) oanc(1:nvar,ind)=displ(1:nvar)
!if(iter>50.and.ind>=1) then
!  call gdiis(idiis,iter,ogint,oanc,displ,ok)
!   ind=-10
!  gnorm=dnrm2(nvar,gint,1)
!  gmax=maxval(gint)
!  if(ok) goto 42
!endif


!if(abs(dE)<(1e-5).and.iter>=1.and.iopt==21) then
! iopt=22
! print*,'Switching to SI-RFO'
!endif

select case(iopt)
 case(21,31)
  call getLambdaRFO(displ,mode,Lam)
! broken  call splitRFO(displ,mode,Lam)
!  if(amber) call FastgetLambdaRFO(displ,mode,Lam) ! "no" gain
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

!42 continue
! check step
do i=1,nvar
   if(abs(displ(i)).gt.maxd) then
      if(displ(i) < 0) displ(i)=-maxd
      if(displ(i) > 0) displ(i)= maxd
   endif
enddo

dnorm=dnrm2(nvar,displ,1)
maxdipl=maxval(abs(displ))

oldanc=anc
anc=anc+displ

hold=hint ! hess damping

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
 if(dE==0d0.and..not.amber) stop 'ERROR: no energy change!'
endif

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
   if(restrain) call printFinalRestrain
   return
  endif
end select


!if(abs(dE)<=1e-5) then
!  xyz0=xyz
!  call gethess
!  call copt
!  return
!endif


call wrestart
enddo ! end main loop

if(restrain) call printFinalRestrain
print*,'FAILED!'
return

end subroutine ancopt


