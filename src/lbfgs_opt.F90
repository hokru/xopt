subroutine lbfgs_opt()
use fiso, only: stdout, r8
implicit none

#ifdef LBFGS


#else
call error('L-BFGS lib not compiled')
#endif

end subroutine



!   ierror=0
!   !lwork=2*nhist+(2*nhist+1)*npar
!   lwork=2*nhist+(2*nhist+1)*nfit
!   allocate(work(lwork))
!   allocate(diag(nfit))

!   ! init
!   call tn_sfun(fit,fitgrad)
!   CALL LBFGS(nfit,nhist,fitpar,fit,fitgrad,.false.,diag,1,xtol,1e-7,work,ierror)
!   i=1
!   do
!     call tn_sfun(fit,fitgrad)
!     CALL LBFGS(nfit,nhist,fitpar,fit,fitgrad,.false.,diag,1,xtol,1e-7,work,ierror)
!   !  call update_params()
!    call update_params(npar,params,nfit,fitpar,iopt)
!     if(ierror==0) exit
!   !  if(ierror<0.or.ierror>1) then
!   !   print*,'error code ',ierror
!   !   stop 'L-BFGS error'
!   !  endif
!     i=i+1
!     if(i>1000) then
!       print*,'max iter (1000) reached'
!       exit
!     endif
!   enddo
! endif

