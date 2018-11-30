subroutine getHess
! provide start hessian
use fiso, only: r8,stdout
use logic, only: hmodel,large
use internals, only: hdiag
use parm, only: B,chess
implicit none
logical transform

write(stdout,'(2x,a)') 'forming model Hessian'
transform=.false.
select case(hmodel)
 case('numerical')
     call numhess
    allocate(b(1,1),hdiag(1))
 case('lindh')
    call make_primitives()
    call modlindh
    transform=.true.
 case('hdiag')
    call getunitH()
    allocate(b(1,1),hdiag(1))
 case('mass') ! no gain from first tests
    call getunitH()
    call Hmass(chess,'nope')
    allocate(b(1,1),hdiag(1))
 case('simple')
    call make_primitives()
    call modsimple
    transform=.true.
 case('alm')
    call make_primitives()
    call modalmloef
    transform=.true.
 case default
end select

if(transform) then
  write(stdout,'(2x,a)') 'transforming model H(int) to H(cart)'
  call Hdiag2Cart()
endif



end subroutine getHess


subroutine getunitH
! cart. unit hessian
use fiso, only: r8
use parm, only: nat, chess,i
implicit none

chess=0.0_r8
do i=1,3*nat
  chess(i,i)=1.0_r8
enddo
end subroutine getunitH


subroutine getscaledH
use parm, only: nat, chess,i,j,iat
use fiso, only: r8
use atomdata, only: ams2
use constant, only: amu2au
implicit none
real(r8) factor
chess=0.0_r8
i=0
do j=1,nat
    factor=1.0+1.0/sqrt(ams2(iat(j)))
    factor=1.0/(ams2(iat(j)))
    i=1+3*(j-1)
    chess(i,i)=1.0_r8*factor
    chess(i+1,i+1)=1.0_r8*factor
    chess(i+2,i+2)=1.0_r8*factor
    print*,j,i,i+1,i+2,chess(i,i)
enddo
end



subroutine modlindh
! diagonal form!
use parm, only: xyz,nat
use fiso, only: r8,stdout
use internals
use logic, only: debug
implicit none
real(r8) dbond
real(r8), parameter :: kb=0.45,ka=0.15,kt=0.005
integer i,a,b,idx,c,d
real(r8) rab,rbc,rcd,lindh_rho

hdiag=0
idx=0
do i=1,int_nb
  idx=idx+1
  a=int_bcast(i,1)
  b=int_bcast(i,2)
  rab=dbond(xyz(1,a),xyz(1,b))
  hdiag(idx)=kb*lindh_rho(a,b,rab)
enddo

do i=1,frag_nh+frag_nb
  idx=idx+1
  a=frag_bcast(i,1)
  b=frag_bcast(i,2)
  rab=dbond(xyz(1,a),xyz(1,b))
  hdiag(idx)=kb*lindh_rho(a,b,rab)
enddo

do i=1,int_na
  idx=idx+1
  a=int_acast(i,1)
  b=int_acast(i,2)
  c=int_acast(i,3)
  rab=dbond(xyz(1,a),xyz(1,b))
  rbc=dbond(xyz(1,b),xyz(1,c))
  hdiag(idx)=ka*lindh_rho(a,b,rab)*lindh_rho(b,c,rbc)
enddo

do i=1,int_nt
  idx=idx+1
  a=int_tcast(i,1)
  b=int_tcast(i,2)
  c=int_tcast(i,3)
  d=int_tcast(i,4)
  rab=dbond(xyz(1,a),xyz(1,b))
  rbc=dbond(xyz(1,b),xyz(1,c))
  rcd=dbond(xyz(1,c),xyz(1,d))
  hdiag(idx)=kt*lindh_rho(a,b,rab)*lindh_rho(b,c,rbc)*lindh_rho(c,d,rcd)
enddo

! We cannot afford that an internal is not counting (close to zero) since
! we have diagnoal non-redudant internals!
do i=1,int_nb+int_na+int_nt+frag_nh+frag_nb
if (hdiag(i)<0.05) then
  if(debug) write(stdout,*) 'small Lindh hess value! Adjusting!',hdiag(i),i
  hdiag(i)=0.05
endif
enddo

end subroutine

subroutine modsimple
! diagonal form!
use parm, only: xyz,i,j,nat
use fiso, only: r8
use internals
implicit none
integer idx
! real(r8) hdiag(nints)

idx=0
do i=1,int_nb
  idx=idx+1
  hdiag(idx)=0.5
enddo
do i=1,frag_nh+frag_nb
  idx=idx+1
  hdiag(idx)=0.2
enddo
do i=1,int_na
  idx=idx+1
  hdiag(idx)=0.2
enddo
do i=1,int_nt
  idx=idx+1
  hdiag(idx)=0.1
enddo
! this could be substituted with a D2 like expression
do i=1,int_nLJ
  idx=idx+1
  ! hdiag(idx)=0.2/(int_LJval(i))**2
   hdiag(idx)=0.05 !0.01+2.0/(int_LJval(i))
enddo
end subroutine


subroutine modalmloef
! diagonal form!
use parm, only: xyz,i,j,nat,iat
use fiso, only: r8
use internals
use atomdata, only: rcov
implicit none
integer idx,a,b,c,d
real(r8) :: rab,rbc,cab,cbc,dbond
real(r8) :: f1,f2,f3,f4,f5
integer :: nb


idx=0
do i=1,int_nb
  idx=idx+1
  a=int_bcast(i,1)
  b=int_bcast(i,2)
  cab=rcov(iat(a))+rcov(iat(b))
  rab=dbond(xyz(1,a),xyz(1,b))
  f1=0.3601
  f2=1.944
  hdiag(idx)=f1*exp(-f2*(rab-cab))
enddo

! fragment + H-bonds
do i=1,frag_nh+frag_nb
  idx=idx+1
  hdiag(idx)=0.2
enddo

do i=1,int_na
  idx=idx+1
  a=int_acast(i,1)
  b=int_acast(i,2)
  c=int_acast(i,3)
  f1=0.089
  f2= 0.11
  f3= 0.44
  f4=-0.42
  cab=rcov(iat(a))+rcov(iat(b))
  cbc=rcov(iat(b))+rcov(iat(c))
  rab=dbond(xyz(1,a),xyz(1,b))
  rbc=dbond(xyz(1,b),xyz(1,c))
  hdiag(idx)= f1+f2/((cab*cbc)**f4)*exp(-f3* (rab-cab+rbc-cbc) )
enddo

do i=1,int_nt
  idx=idx+1
  ! a=int_tcast(i,1)
  b=int_tcast(i,2)
  c=int_tcast(i,3)
  ! d=int_tcast(i,4)
  f1= 0.0015
  f2= 14.0
  f3= 2.85
  f4= 0.57
  f5= 4.00
  nb=cn(b)+cn(c)-2 ! -2 because central bond is twice in the CN
  rbc=dbond(xyz(1,b),xyz(1,c))
  cbc=rcov(iat(b))+rcov(iat(c))
  hdiag(idx)=f1+f2*nb**f4/(rbc-cbc)**f5*exp(-f3*(rbc-cbc))
enddo


end subroutine




subroutine numhess
! numerical cartesian hessian
use parm
use logic
use progs
use fiso, only: r8,stdout
implicit none
integer ii,jj,hi,hj
real(r8)  step
real(r8)  xx,thrE
real(r8), allocatable ::  hess(:,:), gr(:,:),gl(:,:)

! machine precision epsilon
thrE=epsilon(1.0_r8)

step=0.005_r8
  write(stdout,*) 'Doing Hessian numerically ..'
  write(stdout,*) '#displacement: ',6*nat
  allocate(hess(nat3,nat3),gl(3,nat),gr(3,nat))
  hess=0.0_r8
  do i=1,nat
  if(int(mod(i,nint(nat/10.0))).eq.0) write(stdout,'(I3,A)') nint(100_r8*i/dble(nat)),' % done'
   do j=1,3
    hi=(i-1)*3+j
    xyz(j,i)=xyz(j,i)+step
    call newxyz(nat,iat,xyz) ! write new xopt.xyz
    call getgrad
    gr=grad
    xyz(j,i)=xyz(j,i)-2_r8*step
    call newxyz(nat,iat,xyz)
    call getgrad
    gl=grad
    xyz(j,i)=xyz(j,i)+step
  do ii=1,nat
   do jj=1,3
    hj=(ii-1)*3+jj
    xx=(gr(jj,ii)-gl(jj,ii))/(2_r8*step)
    if(abs(xx).gt.thrE)  hess(hi,hj)=xx
   enddo ! jj-loop
  enddo  ! ii-loop

 enddo ! j-loop
 enddo  ! i-loop

xyz=xyz0
  deallocate(gl,gr)

! symmetrize
  write(stdout,*) 'Symmetrizing Hessian ..'
  do i=1,nat3
     do j=1,nat3
        chess(j,i)=(hess(i,j)+hess(j,i))*0.5_r8
     enddo
  enddo

 write(stdout,'(a)')  'xopt cartesian Hessian  written into .. ''xopt.chess'' '
 call wrhess(nat3,chess,'xopt.chess')


! do mass weighting+projection
call Hmass(chess,'project')
!~ call DiagSM(nat3,chess,eig)
! and printing of Freq+ZPVE, and g98fake output

deallocate(hess)
end subroutine





! HELPER FUNCTIONS FOR LINDH Hessian
real(r8) function lindh_REF(a,b)
use fiso, only: r8
implicit none
integer a,b,c
real(r8) val

! sawp if not b>a
if(a>b) then
 c=a
 a=b
 b=c
endif

val=0d0
if(a==1.and.b==1) then
  val=1.35_r8
elseif(a==1.and.b==2) then
 val=2.1_r8
elseif(a==1.and.b>2) then
 val=2.53_r8
elseif(a==2.and.b==2) then
  val=2.87_r8
elseif(a==2.and.b>2) then
  val=3.40_r8
elseif(a>2) then
 val=3.40_r8
endif

if(val<0.001) call error('this must not happend')
lindh_REF=val
end function

real(r8) function lindh_a(a,b)
use fiso, only: r8
implicit none
integer a,b,c
real(r8) val

if(a==1.and.b==1) then
 val=1.00_r8
elseif(a==1.and.b>1) then
 val=0.3939_r8
elseif(a>1.and.b>1) then
 val=0.280_r8
endif
if(val<0.001) call error('this must not happend')
lindh_a=val
end function

real(r8) function lindh_rho(a,b,rab)
use fiso, only: r8
implicit none
integer a,b,pa,pb,period
real(r8) rab,alpha,ref,lindh_REF,lindh_a

pa=period(a)
pb=period(b)

ref=lindh_REF(pa,pb)
alpha=lindh_a(pa,pb)

lindh_rho=exp(-alpha*(rab**2-ref**2))
end function


subroutine HcartD3(nat,iat,xyz,chess)
! adds numerical cartesian hessian for d3
!USE OMP_LIB
! use parm, only: nat,iat,xyz,nat3,nat,chess
use fiso, only: r8,stdout
use logic, only: d3func,d3
implicit none
integer nat,nat3,iat(nat),nat2
integer ii,jj,hi,hj,i,j,io,tid
real(r8), parameter ::  step=0.005_r8
real(r8), parameter ::  thrE=1e-12_r8
real(r8)  xx,xyz(3,nat),e1,e2,timer
! real(r8), allocatable ::  hess(:,:), gr(:,:),gl(:,:)
real(r8)  hess(nat*3,nat*3), gr(3,nat),gl(3,nat),chess(nat*3,nat*3)
character(255) :: arg(10)
logical switchoff

!call error('disabled!') ! re-enable d3 lib

switchoff=.false.

! #ifdef DFTD3
arg=''
arg(1)='-func'
arg(2)='hf'
arg(3)='-zero'
arg(4)='-grad'
arg(5)='-cnthr'
arg(6)='20'
arg(7)='-cutoff'
arg(8)='20'
nat3=nat*3
d3func=' -func hf -zero -cnthr 20 -cutoff 20'
if(.not.d3) then
  d3=.true.
  switchoff=.true.
endif
print *,'D3: ',trim(d3func)
! step=0.005_r8
print*, '#displacement: ',6*nat
 ! allocate(hess(nat3,nat3),gl(3,nat),gr(3,nat))
call status1(timer)
hess=0.0_r8
e1=0
e2=0
do i=1,nat
   if(int(mod(i,nint(nat/10.0))).eq.0) then
    write(stdout,'(2x,I3,A)') nint(100_r8*i/dble(nat)),' % done'
    call status2(timer) 
    call status1(timer)
   endif
   do j=1,3
      hi=(i-1)*3+j
      xyz(j,i)=xyz(j,i)+step
      call getd3gcp(nat,e1,gr)
      ! call call_dftd3(nat,iat,xyz,gr,e1,.false.,arg,6)
      xyz(j,i)=xyz(j,i)-2_r8*step
      call getd3gcp(nat,e2,gl)
      ! call call_dftd3(nat,iat,xyz,gl,e2,.false.,arg,6)
      xyz(j,i)=xyz(j,i)+step
      do ii=1,nat
        do jj=1,3
          hj=(ii-1)*3+jj
          xx=(gr(jj,ii)-gl(jj,ii))/(2_r8*step)
          if(abs(xx).gt.thrE)  hess(hi,hj)=xx
        enddo ! jj-loop
      enddo  ! ii-loop
    enddo ! j-loop
enddo  ! i-loop
!  deallocate(gl,gr)

! symmetrize
  print*, 'adding Hessians ..'
  do i=1,nat3
     do j=1,nat3
        chess(j,i)=chess(j,i)+(hess(i,j)+hess(j,i))*0.5_r8
          ! chess(j,i)=(hess(i,j)+hess(j,i))*0.5_r8
     enddo
  enddo
! call printMat(stdout,nat3,nat3,chess,"H(D3)")
 write(stdout,'(a)')  'cartesian D3-Hessian done'
! do mass weighting+projection
! call Hmass(chess,'project')
! call DiagSM(nat3,chess,eig)
! and printing of Freq+ZPVE, and g98fake output

! deallocate(hess)
! #else
!   call error('Needs internal dftd3 library!')
! #endif
if ( switchoff) d3=.false.

end subroutine


