! GDIIS for normal coords
! buggy??
subroutine gdiis(idiis,iter,ogint,oanc,displ,xgint,ok)
use parm
use fiso, only: r8,stdout
implicit none
integer ind,num,iswitch,idiis,info,iter
integer ipiv(idiis+1),ic
real(r8) e(nvar,idiis),a(idiis+1,idiis+1)
real(r8) x(idiis+1),displ(nvar),s
real(r8) ogint(nvar,idiis),ddot,xgint(nvar)
character(80) wx
real(r8) dnrm2,qanc(nvar),oanc(nvar,idiis)
logical ok



! calculate error vectors e
do i=1,idiis
 e(1:nvar,i)=ogint(1:nvar,i)
enddo

! form gdiis matrix a
do i=1,idiis
 do j=1,idiis
  a(i,j)=ddot(nvar,e(1,i),1,e(1,j),1)
  if(i==j) a(i,j)=a(i,j)*1.2_r8 ! stabilize
 enddo
enddo
a(1:idiis,idiis+1)=-1.0_r8
a(idiis+1,1:idiis)=-1.0_r8
a(idiis+1,idiis+1)=0.0_r8

x=0_r8
x(idiis+1)=-1.0_r8

! solve linear equation
call dgesv(idiis+1,1,a,idiis+1,ipiv,x,idiis+1,info)

if(info/=0) then
 write(stdout,*) 'info:', info
 call error('GDIIS error (dgesv failed)')
endif

!print*,'diis constraint',sum(x)
!if((sum(x)-1d0).gt.0.005) return

xgint=0d0
! form new gradient
do j=1,nvar
 do i=1,idiis
  xgint(j)=xgint(j)+x(i)*ogint(j,i)
 enddo
enddo

if( dnrm2(nvar,gint,1) .lt. dnrm2(nvar,xgint,1) ) then
 write(stdout,*)'     no gdiis step'
 return ! skip gdiis
endif

! Farkas rule c)
s=0
do i=1,idiis
 if(x(i)>0d0) s=s+x(i)
enddo
if(s>15) then
 write(stdout,*) 'no gdiis (large extrapolation)'
 return ! skip gdiis
endif

write(stdout,*)'     --> gdiis step'
return

end subroutine
