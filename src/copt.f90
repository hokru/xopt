
!***********************
!* cartesian optimizer *
!***********************
subroutine copt
! simple, cartesian optimizer for reference
use parm
use logic
use popt
use progs
use fiso, only: r8, stdout
implicit none
integer iter, info
integer nat1,lwork,liwork
real(r8), allocatable :: displ(:)
real(r8), allocatable :: Uaug(:,:), eaug(:),aux(:)
real(r8) f,dE,gmax,maxdipl,dnorm,oldE,oldG(3,nat),oldD(nat3)
character(2) esym
real(r8) ddot,dnrm2,rwork,dummy,time
real(r8) hs(nat*3*(nat*3+1)/2),freqcm
logical EE,GG,GM,DM,doupdate
integer, allocatable :: iwork(:),isuppz(:)
real(r8), parameter :: null=0.0_r8
!diis
integer ind
real(r8) newxyz(3,nat),q(nat3),qk(nat3),newgrad(3,nat)
real(r8) evec(idiis,nat3)

!gdiis.f
real(r8) ogint(nat3,idiis),oanc(nat3,idiis),gdispl(1:nat3)
logical ok

if(tsopt) mode=tsmode+1 ! +1 because of augmented Hessian

do_gdiis=.false. ! matrix inversion failing? hessian singular?
!do_gdiis=.true.

!debug=.false.
EE=.false.
GG=.false.
GM=.false.

nat3=nat*3
nat1=nat3+1
lwork = 1 + 6*nat1 + 2*nat1**2
liwork=3+5*nat1
info=123

allocate(displ(nat3))
allocate(Uaug(nat1,nat1),eaug(nat1))

energy=null
oldG=null
oldD=null
displ=null
Doupdate=.false.

! main loop
call wrxyz('xopt.log')
call touchlog

! project out trans/rot from input hessian
!if(restart) then
print*,' removing trans/rot from input hessian'
k=0
do i=1,nat3
   do j=1,i
      k=k+1
      hs(k)=0.5*(chess(i,j)+chess(j,i))
      if(i.ne.j .and. abs(hs(k)).lt.1d-12) hs(k)=null
   enddo
enddo
call TRpr(nat,xyz0,hs)
call packM(nat3,chess,hs,'unpack')
!endif

if(tsopt) then
  call DiagSM(nat*3,chess,gdispl)
  write(*,'(2x,a)') 'lowest projected vib. freq (cm-1) of non-weighted input Hessian after shift:'
  call printvib(nat*3,gdispl)
  print*, 'Following ',tsmode,freqcm(gdispl(tsmode)),' cm-1 '
endif

gdispl=null
evec=null
ind=0
write(stdout,'(2x,a)') ''
write(stdout,'(2x,a)') 'Starting cartesian geometry optimization'
! write(stdout,'(2x,a)') 'Iter          Energy         dE           gnorm         gmax           dnorm       maxdispl       lambda'
!write(stdout,'(2x,a)') 'Iter          Energy         dE           gnorm         gmax      dnorm     maxdispl     lambda     dE(pred)'
call pchead()
do iter=iiter,maxiter
call check_stop()

! energy+gradient
oldE=energy
oldG=grad
oldD=displ

call wrxyz(trim(xyzfile))
call getgrad
if(restrain) call addRestrainGrad
if (freeze) call freezeGrad
call loggrad(nat,iter,grad,energy)

dE=energy-oldE
call cartgnorm(nat,grad,gnorm)
gmax=maxval(grad)
call appxyz('xopt.log')

if(debug) call debug1('Hupdate',time)
! update hessian
call Hupdate(iter,chess,hs,nat3,grad,oldG,displ,1,Doupdate)
!call wrhess(nat3,chess,'xopt.hess')
Doupdate=.true.
if(debug) call debug2(time)



if(debug) call debug1('forming augmented Hessian',time)
! form aug Hessian
Uaug = null
do i=1,nat3
   do j=1,nat3
      Uaug(j,i)=chess(j,i)
   enddo
enddo

k=0
do i=1,nat
 do j=1,3
  k=k+1
  Uaug(k,nat1)=grad(j,i)
  Uaug(nat1,k)=grad(j,i)
 enddo
 enddo
if(debug) call debug2(time)


if(debug) call debug1('RFO diag',time)
!call DiagSM3_lowest(nat1,Uaug,eaug,1,1d-13)
call DiagSM(nat1,Uaug,eaug)
if(debug) call debug2(time)

if(abs(Uaug(nat1,mode)).lt.1.e-10) call error('internal RF error')
displ(1:nat3)=Uaug(1:nat3,mode)/Uaug(nat1,mode)



k=0
do i=1,nat
 do j=1,3
  k=k+1
  qk(k)=grad(j,i)
  q(k)=oldG(j,i)
 enddo
enddo

!call conjgrad(nat3,qk,q,oldD,displ)
!oldD=displ
!if(abs(eaug(1))>1e-4) then
!  f=0.2
!else
!f=0
!endif
!f=0
!displ(1:nat3)=(1-f)*(Uaug(1:nat3,1)/Uaug(nat1,1))+f*oldD


! check step
if(.not.do_gdiis) then
do i=1,nat3
   if(abs(displ(i)).gt.maxd) then
      if(displ(i) < 0.0_r8) displ(i)=-maxd
      if(displ(i) > 0.0_r8) displ(i)= maxd
   endif
enddo
endif
dnorm=dnrm2(nat3,displ,1)
maxdipl=maxval(displ)

k=0
do i=1,nat
 do j=1,3
  k=k+1
!  q(k)=xyz(j,i) ! old
  xyz(j,i)=(xyz(j,i)+displ(k))
!  qk(k)=grad(j,i) ! new
!  qk(k)=xyz(j,i) ! new
 enddo
enddo

!if(iter>idiis)   xyz=newxyz

EE=.false.
GG=.false.
GM=.false.
DM=.false.
if(abs(dE)<=econv) EE=.true.
if(gnorm<=gconv) GG=.true.
if(gmax<=maxgrad) GM=.true.
if(maxdipl<=dconv) DM=.true.

if(iter==1) then
 dE=0d0 
endif
call pciter(iter,energy,dE,gnorm,gmax,dnorm,maxdipl,eaug(mode),EE,GG,GM,DM,0d0)


if(EE.and.restrain.or.(EE.and.GG.and.GM.and.DM)) then
print*,'CONVERGED !'
write(stdout,'(3x,a)') '               criteria   actual value'
write(stdout,'(3x,a,2(ES10.3,2x))') ' energy    = ',econv,dE
write(stdout,'(3x,a,2(ES10.3,2x))') ' gnorm     = ',gconv,gnorm
write(stdout,'(3x,a,2(ES10.3,2x))') ' max grad  = ',maxgrad,gmax
write(stdout,'(3x,a,2(ES10.3,2x))') ' max displ = ',dconv,dnorm
if(restrain) call printFinalRestrain
  return
endif


!if(iter==2000) then
!  xyz0=xyz
!  call getanc
!  call ancopt(iter+1)
!  return
!endif

call wrestart(iter)
enddo ! end main loop

!if(nat<2000) call wrhess(nat3,chess,'xopt.hess')
if(restrain) call printFinalRestrain

write(stdout,*) 'FAILED!'
return

end subroutine copt


subroutine cartgnorm(nat,g,gnorm)
use fiso, only: r8
implicit none
integer nat
real(r8) g(3,nat),gnorm,dnrm2

gnorm=dnrm2(nat,g(1,1),1)
gnorm=gnorm+dnrm2(nat,g(2,1),1)
gnorm=gnorm+dnrm2(nat,g(3,1),1)
!print*, gnorm
end
