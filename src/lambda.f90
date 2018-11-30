
! single precision version
subroutine FastgetLambdaRFO(displ,mode,Lam)
use parm, only: i,j,k,hint,gint,nvar
use fiso, only: r8,r4
implicit none
integer lwork,nvar1,mode
real(r4), allocatable :: auxF (:)
real(r4), allocatable :: eaugF(:)
real(r4), allocatable :: UaugF(:,:)

real(r8) displ(nvar)
real(r8) lam,tmp2
logical debug
real(r8) time
debug=.false.

nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
allocate(UaugF(nvar+1,nvar+1),eaugF(nvar+1), auxF(lwork))

! augment hessian by grad:
if(debug) call debug1('float Uaug',time)
  UaugF = 0_r4
  k = 0
  do i=1,nvar
     do j=1,i
        k=k+1
        UaugF(i,j)=real(hint(k),4)
        UaugF(j,i)=real(hint(k),4)
     enddo
  enddo
  UaugF(1:nvar,nvar+1)=real(gint(1:nvar),4)
  UaugF(nvar+1,1:nvar)=real(gint(1:nvar),4)

if(debug) call debug2(time)

if(debug) call debug1('DiagSM3_lowest_SP',time)
call DiagSM3_lowest_SP(nvar+1,UaugF,eaugF,mode,1e-13)
if(debug) call debug2(time)


tmp2=real(UaugF(nvar+1,mode),8)
if(abs(tmp2).lt.1e-10) call error('RF failure <FastgetLambdaRFO>')

displ(1:nvar)=real(UaugF(1:nvar,mode),8)/tmp2
Lam=real(eaugF(mode),8)
deallocate(UaugF,eaugF,auxF)
end subroutine FastgetLambdaRFO





subroutine getLambdaRFO(displ,mode,Lam)
use parm, only: i,j,k,hint,gint,nvar
use fiso, only: r8
use logic, only: debug
implicit none
integer lwork,nvar1,info,mode
real(r8), allocatable :: aux (:)
real(r8), allocatable :: eaug(:)
real(r8) displ(nvar)
real(r8), allocatable :: Uaug(:,:)
real(r8) F(nvar),Q,edum,damp,hlow,ddot,temp,lam,lambda,tmp2,tmp3

! condition number
real(r8) cnr

!logical debug
real(r8) time
!debug=.false.

nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
allocate(Uaug(nvar+1,nvar+1),eaug(nvar+1),aux(lwork))


if(debug) call debug1('augHess',time)
! augment hessian by grad:
      Uaug = 0.0_r8
      k = 0
      do i=1,nvar
         do j=1,i
            k=k+1
            Uaug(i,j)=hint(k)
            Uaug(j,i)=hint(k)
         enddo
      enddo
      Uaug(1:nvar,nvar+1)=gint(1:nvar)
      Uaug(nvar+1,1:nvar)=gint(1:nvar)
if(debug) call debug2(time)

!call get_cnr(nvar1,nvar1,Uaug,cnr)

if(debug) call debug1('lambda diag',time)
! get eigenvalues
!else
if(nvar>1000) then
if(debug) print*,'fast diag(DiagSM3_lowest)'
  call DiagSM3_lowest(nvar+1,Uaug,eaug,mode,5d-10)
else
if(debug) print*,'normal diag(DiagSM)'
 ! call DiagSM2(nvar+1,Uaug,eaug)
 call DiagSM(nvar+1,Uaug,eaug)
endif
if(debug) call debug2(time)
! devide by last element
      if(abs(Uaug(nvar1,1)).lt.1.e-10_r8) call error('RF error (getLambdaRFO)')
      displ(1:nvar)=Uaug(1:nvar,mode)/Uaug(nvar+1,mode)
      Lam=eaug(mode)
return
end subroutine getLambdaRFO



subroutine getLambdaSIRFO(displ,mode,Lam)
use parm, only: i,j,k,hint,gint,nvar,nat3
use fiso, only: r8
implicit none
integer lwork,nvar1,info,mode
real(r8), allocatable :: eaug(:)
real(r8), allocatable :: aux (:)
real(r8), allocatable :: Uaug(:,:)
real(r8), allocatable :: S(:,:)
real(r8) displ(nvar),lam
real(r8) soff,sdiag,dnrm2
integer liwork
integer, allocatable :: iwork(:)

nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
liwork= 3 + 5*nvar1

allocate(Uaug(nvar+1,nvar+1),eaug(nvar+1),aux(lwork))
allocate (iwork(liwork))
allocate(S(nvar1,nvar1))


! augment hessian by grad:
      Uaug = 0.0_r8
      k = 0
      do i=1,nvar
         do j=1,i
            k=k+1
            Uaug(i,j)=hint(k)
            Uaug(j,i)=hint(k)
         enddo
      enddo
      Uaug(1:nvar,nvar+1)=gint(1:nvar)
      Uaug(nvar+1,1:nvar)=gint(1:nvar)

! form scaling matrix S
!dnrm2(nvar,gint,1)
      ! sdiag=1.0_r8/sqrt(dble(nvar))
      sdiag=sqrt(dble(nvar))
      soff=0.0_r8
! meh, a really good scaling would improve things. sqrt(natoms) is ok, but the RFO correction
! is be too small sometimes.

      S = 0.0_r8
      k = 0
      do i=1,nvar
         do j=1,i
            k=k+1
            S(i,j)=soff
            S(j,i)=soff
            if(i==j) S(i,j)=sdiag
         enddo
      enddo
      S(1:nvar,nvar1)=0.0_r8
      S(nvar1,1:nvar)=0_r8
      S(nvar1,nvar1)=1_r8

    call dsygvd(1, 'V', 'U', nvar1, Uaug, nvar1, S, nvar1, eaug, aux, lwork, iwork, liwork,info )

! get eigenvalues
!      call dsyev ('V','U',nvar+1,Uaug,nvar+1,eaug,aux,lwork,info)
! devide by last element
      if(abs(Uaug(nvar1,1)).lt.1.e-10_r8) call error('RF error')
      displ(1:nvar)=Uaug(1:nvar,mode)/Uaug(nvar+1,mode)
      Lam=eaug(mode)
!print*,lam,Uaug(nvar+1,mode)
return
end subroutine getLambdaSIRFO



!  p-rfo
subroutine getLambdaPRFO(nvar,hess,gint,displ,lambda,mode)
use fiso, only: r8,stdout
implicit none
integer k,i,j,lwork,nvar1,nvar,info,mode
real(r8), allocatable :: aux (:)
real(r8), allocatable :: eaug(:)
real(r8), allocatable :: he(:)
! real(r8), allocatable :: Uaug(:,:)
real(r8), allocatable :: Hint(:,:),U(:,:)
real(r8) hess(nvar*(nvar+1)/2)
real(r8) displ(nvar), gint(nvar)
real(r8) F(nvar),temp,lam,lambda,freqcm
nvar1=nvar+1
! lwork = 1 + 6*nvar1 + 2*nvar1**2
! allocate(eaug(nvar+1),aux(lwork))


!~ he=0d0
!aux=0d0
allocate(Hint(nvar,nvar),he(nvar),U(nvar,nvar))
! eigenvalues of internal hessian
     Hint = 0.0_r8
     k = 0
     do i=1,nvar
        do j=1,i
           k=k+1
           Hint(i,j)=hess(k)
           Hint(j,i)=hess(k)
          enddo
     enddo
! call dsyev ('V','U',nvar,Hint,nvar,he,aux,lwork,info)
call DiagSM(nvar,Hint,he)

! gradient along local eigenmode of internal hessian
do i=1,nvar
  F(i)=0.0_r8
  do j=1,nvar
    F(i)=F(i)+Hint(j,i)*gint(i)
  enddo
enddo

! P-RFO with newton-raphson lambda from mopac
call NRLamda(F,he,nvar,lambda,mode)

! lambda for maximization
! lam=0.5_r8* he(mode) + ( 0.5_r8* sqrt(he(mode)*he(mode) + 4.0_r8* F(mode)*F(mode)))
lam= he(mode) + sqrt (he(mode)*he(mode) + 4.0_r8* F(mode)*F(mode) )
lam=0.5_r8*lam
write(stdout,'(2x,''| vib. freq. of selected mode: '',F9.2)') freqcm(he(mode))

! ! check & correct lambda
do i=1,nvar
! check if lambda is less than he(i)
 if(lambda.gt.he(i)) lambda=lambda-2*(lambda-he(i))
 if(lambda.lt.0.0_r8) lambda=abs(lambda)*0.1_r8  ! positive and small
enddo


! take the step
displ=0.0_r8
do i=1,nvar
  temp=F(i)/(lambda-he(i))
  if(i==mode) temp=F(i)/(lam-he(i))
  do j=1,nvar
   displ(i)=displ(i)+temp*Hint(j,i)
  enddo
enddo
deallocate(Hint,he,U)

end subroutine getLambdaPRFO

! overlap between eigenvectors
! for TS search
subroutine eigovl(U,mode,nvar,tm)
use fiso, only: r8,stdout
implicit none
integer tm,otm,i,nvar
real(r8) mode(nvar),S,oS
real(r8) U(nvar,nvar)
real(r8) ddot,v1(nvar),v2(nvar)

otm=tm

!v1=0_r8
!v2=0_r8
v1(1:nvar)=U(1:nvar,tm)
v2(1:nvar)=mode(1:nvar)
! normalize
v1(1:nvar)=v1(1:nvar)/sqrt(ddot(nvar,v1,1,v1,1))
v2(1:nvar)=v2(1:nvar)/sqrt(ddot(nvar,v2,1,v2,1))
! overlap with current mode
oS=ddot(nvar,v1,1,v2,1)
oS=abs(oS)
write(stdout,'(2x,''| overlap (%) with current mode is: '',F6.2)') oS*100_r8
do i=1,nvar
  v1(1:nvar)=U(1:nvar,i)
  v1=v1/sqrt(ddot(nvar,v1,1,v1,1))
  S=abs(ddot(nvar,v1,1,v2,1))
  if(S>oS) then
   oS=S
   tm=i
endif

enddo

if(tm.ne.otm) then
 write(stdout,'(2x,''|new transition mode: '',i2,'' with overlap: '',F5.2,''%'')') tm,(oS*100)
 if(oS*100.le.40) then
  call warning('safety stop: please (re-)select mode by hand!')
 endif
endif

return
end subroutine eigovl


! overlap between eigenvectors
! for TS search
subroutine eigovl_int(hint,nvar,hmode,tm,iflag)
use fiso, only: r8,stdout
implicit none
integer tm,otm,i,nvar,iflag
real(r8) hmode(nvar),S,oS
real(r8) U(nvar,nvar),hint(nvar*(nvar+1)/2)
real(r8) ddot,v1(nvar),v2(nvar)

otm=tm

if(iflag==0) then
 hmode=U(1:nvar,tm)
!  print*, hmode(1:6)
 return
endif

call packM(nvar,u,hint,'unpack')

v1(1:nvar)=U(1:nvar,tm)
v2(1:nvar)=hmode(1:nvar)
! print*, tm
! print*, hmode(1:6)
! print*, U(1:6,tm)

! normalize
v1(1:nvar)=v1(1:nvar)/sqrt(ddot(nvar,v1,1,v1,1))
v2(1:nvar)=v2(1:nvar)/sqrt(ddot(nvar,v2,1,v2,1))
! overlap with current mode
oS=ddot(nvar,v1,1,v2,1)
oS=abs(oS)
write(stdout,'(2x,''| overlap (%) with current mode is: '',F6.2)') oS*100_r8
do i=1,nvar
  v1(1:nvar)=U(1:nvar,i)
  v1=v1/sqrt(ddot(nvar,v1,1,v1,1))
  S=abs(ddot(nvar,v1,1,v2,1))
  if(S>oS) then
   oS=S
   tm=i
endif

enddo

if(tm.ne.otm) then
 write(stdout,'(2x,''|new transition mode: '',i2,'' with overlap: '',F5.2,''%'')') tm,(oS*100)
 if(oS*100.le.40) then
  call warning('safety stop: please (re-)select mode by hand!')
 endif
endif

return
end subroutine eigovl_int


! basically from mopac, like everyone does it?!
! find lamda that minimizes all other modes than tsmode
subroutine NRLamda(F,he,nvar,lamda,tsmode)
use fiso, only: r8, stdout
! set tsmode < 1 to minimize all nvar
implicit none
integer nvar,i,j,tsmode,k,ts1,maxiter
real(r8) F(nvar),he(nvar),minstep,big
real(r8) lamda,fff,ff1,ff3,temp,toll,ff2,l1,l2
real(r8) xstep
logical conv,bisect,root

big=1d+3
maxiter=9999
conv=.false.
bisect=.false.
root=.false.
toll=1.0e-8_r8
minstep=0.05_r8

ts1=tsmode+1
if(ts1>2) ts1=1 

! guess: equal zero or less than he(ts1)
lamda=0.0_r8
if(he(ts1)<0.0_r8) lamda=he(ts1)-minstep

do j=1,1000
      fff = 0.0_r8
      ff1 = 0.0_r8
      do  i=1,nvar
            if(i==tsmode) cycle
            fff = fff + (f(i)*f(i))/(lamda-he(i))
            ff1 = ff1 - (f(i)*f(i))/(lamda-he(i))**2
      enddo
      fff = fff - lamda
      ff1 = ff1 - 1.0_r8
      temp = fff/ff1
      lamda = lamda - temp
      write(stdout,'(''lambda: '',d14.6,'' iter: '',i5,'' conv: '',d14.6)') lamda,j,temp
      if (abs(temp)< toll) then
        conv=.true.
        exit
      endif
enddo

return

! check if L is OK
123 continue
if(conv) then
  if(lamda>he(ts1).or.(he(ts1)>0.and.lamda>0)) then
    write(stdout,*) 'lamda, eigenvalue',lamda,he(ts1)
    goto 222 ! try bisect 1 time 
    call error('unacceptable lambda!') ! oh no.. :-(
  else
      write(stdout,*) 'final lambda', lamda
      return ! all good! :-)
  endif
endif

! try bisecting
222 continue
 
if (bisect) then
  write(stdout,*) 'bisecting lambda failed. Potentially unreliable lambda!'
return
! call error('P-RFO lamda bisec failed. Unable to continue') ! we tried before
endif 
write(stdout,*) 'trying to bisect lamda'
bisect=.true.

l1=0.0
if(he(ts1)<0) l1=he(ts1)-minstep
ff1=0.0
do i=1,nvar
      if(i==tsmode) cycle
      ff1=ff1+ (f(i)*f(i))/(l1 - he(i)) 
enddo
ff1=ff1-l1
xstep = abs(l1)/big
if(xstep < 10.0_r8) xstep=10.0_r8                                      
do k=1,maxiter
      ff2=0.0
      l2=l1-k*xstep
      do i=1,nvar
            if(i==tsmode) cycle
            ff2=ff2+ (f(i)*f(i))/(l2 - he(i))  
      enddo
      ff2=ff2-l2
      ! print*, ff1,ff2,l1,l2
      ! if(k==10) stop
      if(ff1*ff2<0.0_r8) then
        root=.true.
        exit
      endif

enddo
write(stdout,*) 'bisect boundaries:',l1,l2

if(root) then
! print*,'bisect boundaries:',l1,l2
do k=1,maxiter

temp=(l1+l2)*0.5_r8
ff3=0.0
do i=1,nvar
 if(i==tsmode) cycle
 ff3 = ff3 + (f(i)*f(i))/(temp - he(i)) 
enddo
ff3=ff3-temp
if(abs(temp-l2)< toll ) then
 lamda=temp
 conv=.true.
 goto 123 ! check conv again
endif
if(ff1*ff3<0.0) then
  l2=temp
else
  l1=temp
endif
enddo
endif

! print*,'Potentially unreliable lambda!'
! call error('no further roots found. NR-lambda failed')s



end subroutine NRLamda
