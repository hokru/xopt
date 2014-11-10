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
!* does not work
subroutine splitRFO(displ,mode,Lam)
use parm
implicit none
integer lwork,nvar1,info,mode,var1,var2,lw2,lw3
real*8, allocatable :: aux (:),aux2(:),aux3(:)
real*8, allocatable :: eaug(:),eaug2(:)
real*8, allocatable :: Uaug(:,:),Uaug2(:,:)
real*8 displ(nvar)
real*8 F(nvar) ,Q,edum,damp,hlow,ddot,temp,lam,lambda,tmp2,tmp3,he(nvar),Hint2(nvar,nvar)
nvar1=nvar+1


var1=int(nvar/2 )
var2=var1+mod(nvar,2)
var2=nvar
print*,var1,var2,nvar

lwork = 1 + 6*var1 + 2*var1**2
lw2 = 1 + 6*var2 + 2*var2**2
allocate(Uaug(var1+1,var1+1),eaug(var1+1),aux(lwork))
allocate(Uaug2(var2+1,var2+1),eaug2(var2+1),aux2(lw2))

!sort hessian (breaks stuff...)
!call bsort(nvar,hint)

! split into pieces
      Uaug = 0d0
      k = 0
      do i=1,var1
         do j=1,i
            k=k+1
            Uaug(i,j)=hint(k)
            Uaug(j,i)=hint(k)
         enddo
      enddo
      Uaug(1:var1,var1+1)=gint(1:var1)
      Uaug(var1+1,1:var1)=gint(1:var1)
! get eigenvalues
      call dsyev ('V','U',var1+1,Uaug,var1+1,eaug,aux,lwork,info)
! devide by last element
      if(abs(Uaug(var1+1,1)).lt.1.d-10)stop'internal RF error'
      displ(1:var1)=Uaug(1:var1,mode)/Uaug(var1+1,mode)
      Lam=eaug(mode)

print*,displ(1:var1)
!print*,lam,'2nd'
! augment hessian by grad:
      Uaug2 = 0d0
      do i=var1+1,nvar
         do j=1,i
            k=k+1
            Uaug2(i,j)=hint(k)
            Uaug2(j,i)=hint(k)
         enddo
      enddo
      Uaug2(var1+1:var2,var2+1)=gint(var1+1:nvar)
      Uaug2(var2+1,var1+1:var2)=gint(var1+1:nvar)
! get eigenvalues
      call dsyev ('V','U',var2+1,Uaug2,var2+1,eaug2,aux2,lw2,info)
! devide by last element
      if(abs(Uaug2(var2+1,1)).lt.1.d-10)stop'internal RF error'
      displ(var1+1:nvar)=Uaug2(var1+1:var2,mode)/Uaug2(var2+1,mode)
      Lam=eaug2(mode)
deallocate(Uaug,Uaug2,aux,aux2,eaug,eaug2)
return
end subroutine splitRFO






subroutine FastgetLambdaRFO(displ,mode,Lam)
use parm
implicit none
integer lwork,nvar1,info,mode
real*8, allocatable :: aux (:)
real*8, allocatable :: eaug(:)
real*8, allocatable :: Uaug(:,:)

real*4, allocatable :: auxF (:)
real*4, allocatable :: eaugF(:)
real*4, allocatable :: UaugF(:,:)

real*8 displ(nvar)
real*8 F(nvar) ,Q,edum,damp,hlow,ddot,temp,lam,lambda,tmp2,tmp3
nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
allocate(Uaug(nvar+1,nvar+1),eaug(nvar+1),aux(lwork))

!call system('touch FastRFO')
! augment hessian by grad:
      Uaug = 0d0
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

allocate(UaugF(nvar+1,nvar+1),eaugF(nvar+1), auxF(lwork))
call floatMat(nvar+1,nvar+1,Uaug,UaugF)
call ssyev ('V','U',nvar+1,UaugF,nvar+1,eaugF,auxF,lwork,info)
call doubleMat(nvar+1,nvar+1,UaugF,Uaug)
call doubleVec(nvar+1,eaugF,eaug)
deallocate(UaugF,eaugF,auxF)

!      call dsyev ('V','U',nvar+1,Uaug,nvar+1,eaug,aux,lwork,info)
! devide by last element
      if(abs(Uaug(nvar+1,1)).lt.1.d-10)stop'internal RF error'
      displ(1:nvar)=Uaug(1:nvar,mode)/Uaug(nvar+1,mode)
      Lam=eaug(mode)
return
end subroutine FastgetLambdaRFO





subroutine getLambdaRFO(displ,mode,Lam)
use parm
implicit none
integer lwork,nvar1,info,mode
real*8, allocatable :: aux (:)
real*8, allocatable :: eaug(:)
real*8, allocatable :: Uaug(:,:)
real*8 displ(nvar)
real*8 F(nvar) ,Q,edum,damp,hlow,ddot,temp,lam,lambda,tmp2,tmp3
nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
allocate(Uaug(nvar+1,nvar+1),eaug(nvar+1),aux(lwork))


! augment hessian by grad:
      Uaug = 0d0
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
! get eigenvalues
      call dsyev ('V','U',nvar+1,Uaug,nvar+1,eaug,aux,lwork,info)
! devide by last element
      if(abs(Uaug(nvar+1,1)).lt.1.d-10)stop'internal RF error'
      displ(1:nvar)=Uaug(1:nvar,mode)/Uaug(nvar+1,mode) 
      Lam=eaug(mode)
return
end subroutine getLambdaRFO



subroutine getLambdaSIRFO(displ,mode,Lam)
use parm
implicit none
integer lwork,nvar1,info,mode
real*8, allocatable :: aux (:)
real*8, allocatable :: eaug(:)
real*8, allocatable :: Uaug(:,:)
real*8, allocatable :: S(:,:)
real*8 displ(nvar)
real*8 F(nvar) ,Q,edum,damp,hlow,ddot,temp,lam,lambda,tmp2,tmp3

real(8) soff,sdiag,dnrm2
integer liwork
integer, allocatable :: iwork(:)

nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
liwork= 3 + 5*nvar1

allocate(Uaug(nvar+1,nvar+1),eaug(nvar+1),aux(lwork))
allocate (iwork(liwork))
allocate(S(nvar1,nvar1))


! augment hessian by grad:
      Uaug = 0d0
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
      sdiag=1d0/sqrt(dble(nvar))
      soff=sdiag*0.01
!      sdiag=(1d0-lam)/sqrt(dble(nvar))
!      sdiag=1d0/sqrt(dble(nat))
!      sdiag=1d0/sqrt(dnrm2(nvar,gint,1))

      S = 0d0
      k = 0
      do i=1,nvar
         do j=1,i
            k=k+1
            S(i,j)=soff
            S(j,i)=soff
            if(i==j) S(i,j)=sdiag
         enddo
      enddo 
      S(1:nvar,nvar1)=0d0
      S(nvar1,1:nvar)=0d0
      S(nvar1,nvar1)=1d0

    call dsygvd(1, 'V', 'U', nvar1, Uaug, nvar1, S, nvar1, eaug, aux, lwork, iwork, liwork,info )

! get eigenvalues
!      call dsyev ('V','U',nvar+1,Uaug,nvar+1,eaug,aux,lwork,info)
! devide by last element
      if(abs(Uaug(nvar+1,1)).lt.1.d-10)stop'internal RF error'
      displ(1:nvar)=Uaug(1:nvar,mode)/Uaug(nvar+1,mode)
      Lam=eaug(mode)
!print*,lam,Uaug(nvar+1,mode)
return
end subroutine getLambdaSIRFO




subroutine getLambdaPRFO(nvar,hess,gint,displ,lambda,mode)
implicit none
integer k,i,j,lwork,nvar1,nvar,info,mode
real*8, allocatable :: aux (:)
real*8, allocatable :: eaug(:)
real*8, allocatable :: he(:)
real*8, allocatable :: Uaug(:,:)
real*8, allocatable :: Hint(:,:),U(:,:)
real*8 hess(nvar*(nvar+1)/2),hm,hdum(nvar*(nvar+1)/2)
real*8 displ(nvar), gint(nvar)
real*8 F(nvar) ,Q,edum,damp,hlow,ddot,temp,lam,lambda
nvar1=nvar+1
lwork = 1 + 6*nvar1 + 2*nvar1**2
allocate(eaug(nvar+1),aux(lwork))


he=0d0
!aux=0d0
allocate(Hint(nvar,nvar),he(nvar),U(nvar,nvar))
! eigenvalues of internal hessian
     Hint = 0d0
     k = 0
     do i=1,nvar
        do j=1,i
           k=k+1
           Hint(i,j)=hess(k)
           Hint(j,i)=hess(k)
        enddo
     enddo
call dsyev ('V','U',nvar,Hint,nvar,he,aux,lwork,info)


! gradient along local eigenmode of internal hessian
do i=1,nvar
F(i)=0.0d0
do j=1,nvar
F(i)=F(i)+Hint(j,i)*gint(i)
enddo
enddo


! P-RFO
!***********************************
   call NRLamda(F,he,nvar,lambda)

! lambda for maximization
!lam=0.5d0* he(mode) + ( 0.5d0* sqrt(he(mode)*he(mode) + 4.0d0* F(mode)*F(mode)))


! check lambda
!do i=1,nvar
! check if lambda is less than and he(i)
!if(lambda.gt.he(i)) lambda=lambda-2*(lambda-he(i))
!if(lambda.lt.0.0d0) lambda=abs(lambda)*0.1  ! positive and small
!enddo


! take the step
displ=0.0d0
!lam=eaug(mode)
do i=1,nvar
 temp=F(i)/(lambda-he(i))
if(i==mode) temp=F(i)/(lam-he(i))
do j=1,nvar
displ(i)=displ(i)+temp*Hint(j,i)
enddo
enddo
!print*,displ

!displ=0
!lam=eaug(mode)
!do i=1,nvar
! temp=gint(i)/(lambda-he(i))
!if(i==mode) temp=gint(i)/(lam-he(i))
!do j=1,nvar
!displ(i)=displ(i)+temp
!enddo
!enddo

end subroutine getLambdaPRFO


subroutine eigovl(U,mode,nvar,tm)
implicit none
real*8 mode(nvar),S,oS
real*8 U(nvar,nvar)
real*8 ddot
integer tm,otm,i,nvar

if(tm==1) return
otm=tm
! write(*,'(''Following '',i2, ''. mode '')') tm

! overlap with current mode
oS=ddot(nvar,U(1,tm),1,mode,1)
oS=abs(oS)
write(*,'(''Overlap with current mode is: '',F7.2)') oS
! test 
do i=1,nvar
! write(*,*) 'hit'
S=abs(ddot(nvar,U(1,i),1,mode,1))
! write(*,*) 'S:',S
if(S.gt.oS) then
oS=S
tm=i
endif
enddo


if(tm.ne.otm) then
write(*,'(''new transition mode: '',i2,'' with overlap: '',F5.2,''%'')') tm,(oS*100)
if(oS*100.le.40) then
 write(*,*) 'safe stop: please (re-)select mode by hand'
! write(*,*) 'eht -hess > ancopt.eht.out'
! call system('eht -hess > ancopt.eht.out')
 stop 'sorry'
endif
endif
return
end subroutine eigovl

! See ef.x
subroutine NRLamda(F,he,nvar,lambda)
!              gradient, hessian eigenvalue, dimension, shift
real*8 F(nvar),he(nvar), lambda
integer nvar,i,j,k
real*8 lamda,fff,ff1, temp,toll,ff2

      TOLL=1.0d-8
      LAMDA=0d0   
!        lamda=eaug(1)                                                             
!     LETS FIRST TRY A NEWTON-RAPHSON WITH ABOVE LAMDA AS GUESS.                
!     IT WILL ALMOST ALWAYS WORK AND IT IS FAST. FFF IS FUNCTION,               
!     FF1 THE GRADIENT OF THE LAMDA FUNCTION.                                   
     do j=1,1000000
      FFF = 0d0                                                                
      FF1 = 0d0                                                                
      FF2=0d0
      DO  I=1,NVAR                                                            
!         FFF = FFF + (gint(I)*gint(I))/(LAMDA-he(I))                            
         FFF = FFF + (F(I)*F(I))/(LAMDA-he(I))                            
!          FF1 = F	F1 - (gint(I)*gint(I))/((LAMDA-he(I))**2)                       
         !FF1 = FF1 - ((gint(I))/(LAMDA-he(I)))**2                       
         FF1 = FF1 - ((F(I))/(LAMDA-he(I)))**2                       
      !ff2=ff2 + ( 2d0* (gint(i)*gint(i)) ) / ( (lamda-he(i))**3 )   
      enddo                                                                  
      FFF = FFF - LAMDA                                                         
      FF1 = FF1 - 1.d0                                                           
      TEMP = FFF/FF1                                                            
      LAMDA = LAMDA -TEMP  
! WRITE(*,'(''Lambda iter: '',D12.6,'' iter: '',I5,'' Diff: '',D12.6)') lamda,j,eaug(1)-lamda
      IF (ABS(TEMP) .LT. TOLL) GOTO 130 
      enddo
130 continue

      lambda=LAMDA
 WRITE(*,'(''NR-Lambda: '',D14.6,'' iter: '',I5)') lamda,j
return
end subroutine NRLamda


      subroutine floatMat(n,m,iM,oM)
      integer n,m,i,j
      real(8),intent(in) :: iM(n,m)
      real(4),intent(out) :: oM(n,m)
      do j=1,m
       do i=1,n
        oM(i,j)=real(iM(i,j))
       enddo
      enddo
      end subroutine

      subroutine floatVec(n,iV,oV)
      integer n,m,i
      real(8),intent(in) :: iV(n)
      real(4),intent(out) :: oV(n)
       do i=1,n
        oV(i)=real(iV(i))
       enddo
      end subroutine

      subroutine doubleMat(n,m,iM,oM)
      integer n,m,i,j
      real(4),intent(in) :: iM(n,m)
      real(8),intent(out) :: oM(n,m)
      do j=1,m
       do i=1,n
        oM(i,j)=dble(iM(i,j))
       enddo
      enddo
      end subroutine

      subroutine doubleVec(n,iV,oV)
      integer n,m,i
      real(4),intent(in) :: iV(n)
      real(8),intent(out) :: oV(n)
       do i=1,n
        oV(i)=dble(iV(i))
       enddo
      end subroutine

