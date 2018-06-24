!  blas/lapack wrapper and general math functions


! symmetric eigensolver, wrapper around simple lapack driver.
! slow for large matrices
subroutine DiagSM(xvar,mat,eig)
use fiso, only: r8
implicit none
real(r8), allocatable :: aux(:)
integer info,lwork,xvar
real(r8) ,intent(in) :: mat(xvar,xvar)
real(r8) xx
real(r8), intent(out) :: eig(xvar)

eig=0.0_r8
call dsyev ('V','U',xvar,mat,xvar,eig,xx,-1,info)
lwork=int(xx)
allocate(aux(lwork))
call dsyev ('V','U',xvar,mat,xvar,eig,aux,lwork,info)
if(info/=0) print*,'Diagonalization failed !!',info
end subroutine


! single-precision variant
subroutine DiagSMF(xvar,matF,eigF)
use fiso, only: r4
implicit none
real(r4), allocatable :: auxF(:)
integer info,lwork,xvar
!real(r8) ,intent(in) :: mat(xvar,xvar)
!real(r8), intent(out) :: eig(xvar)
!real(r8), allocatable :: auxF(:)
real(r4) ,intent(in) :: matF(xvar,xvar)
real(r4) xx
real(r4), intent(out) :: eigF(xvar)

eigF=0.0_r4
call ssyev ('V','U',xvar,matF,xvar,eigF,xx,-1,info)
lwork=int(xx)
allocate(auxF(lwork))
call ssyev ('V','U',xvar,matF,xvar,eigF,auxF,lwork,info)

if(info/=0) print*,'DiagSMF: Diagonalization failed !!',info

end subroutine

!******************************************************************************

! wrapper around divide-and-conquer lapack eiensolver for symmetric problems
! fast, but significant memory requirements
subroutine DiagSM2(xvar,mat,eig)
use fiso, only: r8
implicit none
real(r8), allocatable :: aux(:)
integer, allocatable :: iaux(:)
integer info,lwork,xvar,qi,liwork
real(r8) ,intent(in) :: mat(xvar,xvar)
real(r8) qr
real(r8), intent(out) :: eig(xvar)

eig=0.0_r8
call dsyevd('V','U',xvar,mat,xvar,eig,qr,-1,qi,-1,info)
lwork=int(qr)
allocate(aux(lwork))
allocate(iaux(qi))
liwork=qi
call dsyevd('V','U',xvar,mat,xvar,eig,aux,lwork,iaux,liwork,info)

if(info/=0) print*,'Diagonalization failed [subroutine DiagSM2] !!',info
end subroutine

!******************************************************************************

! wrapper around relatively-robust eigensolver
! fastest and smallest workspace
! but needs to double the input matrix!
! should this be the default? accuracy?
! FULL SOLVER = ALL EIGENVECTORS!
! provide tolerance for solver tol=0 means machine epsilon
subroutine DiagSM3(xvar,mat,eig,tol)
use fiso, only: r8
implicit none
real(r8), allocatable :: aux(:)
integer, allocatable :: iaux(:)
integer info,lwork,xvar,qi,liwork
real(r8) ,intent(inout) :: mat(xvar,xvar)
real(r8) qr
real(r8), intent(out) :: eig(xvar)
integer il,iu,neig,isuppz(2*xvar)
real(r8) vl,vu,tol
real(r8) emat(xvar,xvar)

eig=0
!tol=DLAMCH( 'S' )
!tol=0
call dsyevr('V','A','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,qr,-1,qi,-1,info)
lwork=int(qr)
allocate(aux(lwork))
allocate(iaux(qi))
liwork=qi
call dsyevr('V','A','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,aux,lwork,iaux,liwork,info)
mat=emat

if(info/=0) print*,'Diagonalization failed [subroutine DiagSM3] !!',info
end subroutine


! find lowest neig eigenvectors and eigenvalues and given tolerance
subroutine DiagSM3_lowest(xvar,mat,eig,nlow,tol)
use fiso, only: r8
implicit none
real(r8), allocatable :: aux(:)
integer, allocatable :: iaux(:)
integer info,lwork,xvar,qi,liwork
real(r8) ,intent(inout) :: mat(xvar,xvar)
real(r8) qr
real(r8), intent(out) :: eig(xvar)
integer il,iu,neig,isuppz(2*xvar),nlow
real(r8) vl,vu,tol
real(r8) emat(xvar,xvar) ! this can be smaller!

il=1
iu=nlow
eig=0
neig=0
!tol=DLAMCH('S')
call dsyevr('V','I','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,qr,-1,qi,-1,info)
lwork=int(qr)
allocate(aux(lwork))
allocate(iaux(qi))
liwork=qi
call dsyevr('V','I','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,aux,lwork,iaux,liwork,info)
mat=emat
!print*,'#eigenvalues: ',neig,' of ',nlow,' requested'
if(info/=0) print*,'Diagonalization failed [subroutine DiagSM3_lowest] !!',info
end subroutine

! single-precision
subroutine DiagSM3_lowest_SP(xvar,mat,eig,nlow,tol)
use fiso, only: r4
implicit none
real(r4), allocatable :: aux(:)
integer, allocatable :: iaux(:)
integer info,lwork,xvar,qi,liwork
real(r4) ,intent(inout) :: mat(xvar,xvar)
real(r4) qr
real(r4), intent(out) :: eig(xvar)
integer il,iu,neig,isuppz(2*xvar),nlow
real(r4) vl,vu,tol
real(r4) emat(xvar,xvar)

tol=real(tol,4) ! just to be sure
il=1
iu=nlow
eig=0
neig=0
!tol=DLAMCH('S')
call ssyevr('V','I','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,qr,-1,qi,-1,info)
lwork=int(qr)
allocate(aux(lwork))
allocate(iaux(qi))
liwork=qi
call ssyevr('V','I','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,aux,lwork,iaux,liwork,info)
mat=emat
!print*,'#eigenvalues: ',neig,' of ',nlow,' requested'
if(info/=0) print*,'Diagonalization failed [subroutine DiagSM3_lowest] !!',info
end subroutine



! RR-solver for just the eigenvalues
! with given tolerance
subroutine DiagSM3_N(xvar,mat,eig,tol)
use fiso, only: r8
implicit none
real(r8), allocatable :: aux(:)
integer, allocatable :: iaux(:)
integer info,lwork,xvar,qi,liwork
real(r8) ,intent(inout) :: mat(xvar,xvar)
real(r8) qr
real(r8), intent(out) :: eig(xvar)
integer il,iu,neig,isuppz(2*xvar)
real(r8) vl,vu,tol,dlamch
real(r8) emat(xvar,xvar)

eig=0
tol=DLAMCH( 'S' )
call dsyevr('N','A','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,qr,-1,qi,-1,info)
lwork=int(qr)
allocate(aux(lwork))
allocate(iaux(qi))
liwork=qi
call dsyevr('N','A','U',xvar,mat,xvar,vl,vu,il,iu,tol,neig,eig,emat,xvar,isuppz,aux,lwork,iaux,liwork,info)
mat=emat
if(info/=0) print*,'Diagonalization failed [subroutine DiagSM3_N] !!'
end subroutine

!**********************************************

subroutine mat_svp(A,m,n,k,U,S,VT)
! SVP of matrix A
! this uses the full dimensions!
! U(n,n) * Diag(S(n,m) * VT(m,m)
use fiso, only: r8
implicit none
integer, intent(in) :: n,m,k
real(r8),intent(in):: A(m,n)
real(r8) :: S(k),U(m,m),VT(n,n)
integer :: lwork,info,ldu,ldv
real(r8), allocatable :: work(:)
real(r8) qwork
call DGESVD('A','A', m, n, A , m, S, U, m , VT, n, qwork, -1, info )
lwork=int(qwork)
allocate(work(lwork))
call DGESVD('A','A', m, n, A , m, S, U, m, VT, n, work, lwork, info )
deallocate(work)
if(info.ne.0) call error('SVD orthogonalization failed')

end subroutine

! LAPACK-based symmetric matrix inverse (LU decomposition)
subroutine matinv(A,n)
use fiso, only: r8
 implicit none
 integer :: n,ipiv(n),info
 real(r8) :: A(n,n),work(n)
 ! LU factorization
 call DGETRF(n, n, A, n, ipiv, info)
 if (info /= 0) call error('matinv: matrix is numerically singular!')
 ! inverse through LU factorization
 call DGETRI(n, A, n, ipiv, work, n, info)
 if (info /= 0) call error('matinv: Matrix inversion failed!')
end subroutine matinv


subroutine mat_pinv(A,m,n,Ainv,cutoff)
! svd-based computation of the pseudo-inverse
use fiso, only: r8
implicit none
integer, intent(in)     :: n,m
integer k,i,j
real(r8),intent(in)  :: A(n,m)
real(r8),intent(out) :: Ainv(m,n)
real(r8), allocatable :: U(:,:),VT(:,:),sdiag(:)
real(r8), optional :: cutoff
real(r8), parameter :: alpha=1.0_r8
real(r8), parameter :: beta=0.0_r8

k=min(m,n)
allocate(U(m,k),VT(k,n),sdiag(k))

if(m<n) call error('m<n not allowed!')
if(.not.present(cutoff)) cutoff=1e-10_r8

!1) SVD  ->  A=U*S*V^t
call mat_svp(A,m,n,min(m,n),U,sdiag,VT)

!2) set small s-values to zero; invert s-values
do i=1,k
 if(sdiag(i)>cutoff) then
   sdiag(i)=1.0_r8/sdiag(i)
 else
  sdiag(i)=0.0_r8
 endif
enddo

! ainv=V*S*U^t
do j = 1, k
   call dSCAL( M, sdiag(j), U(1,j),1)
end do
call DGEMM( 't', 't', N, M, n, alpha, VT, n, U, M,beta, ainv, N)
! return transpose?

end subroutine


subroutine unitvec(x,e)
use fiso, only: r8
! unit vector of vector
implicit none
real(r8) x(3),e(3),t(3)
t=DOT_PRODUCT(x,x)
e=x/sqrt(t)
end


subroutine veclen(x,v)
use fiso, only: r8
! vector(3) norm
implicit none
real(r8) x(3),v
v=sqrt(dot_product(x,x))
end subroutine


subroutine veclen2(a,b,v)
use fiso, only: r8
! difference vector length
implicit none
real(r8) a(3),b(3),v,x(3)
x=a-b
v=sqrt(dot_product(x,x))
end subroutine


subroutine cross_prod(y,x2,x3)
use fiso, only: r8
! cross product
implicit none
real(r8) y(3),x2(3),x3(3)
  y(1) =  x2(2)*x3(3) - x3(2)*x2(3)
  y(2) = -x2(1)*x3(3) + x3(1)*x2(3)
  y(3) =  x2(1)*x3(2) - x3(1)*x2(2)
end subroutine

! floating point precision trafo
! Is this actually any good?
subroutine floatMat(n,m,iM,oM)
use fiso, only: r8,r4
integer n,m,i,j
real(r8),intent(in) :: iM(n,m)
real(r4),intent(out) :: oM(n,m)
do j=1,m
 do i=1,n
  oM(i,j)=real(iM(i,j))
 enddo
enddo
end subroutine

subroutine floatVec(n,iV,oV)
use fiso, only: r8,r4
integer n,i
real(r8),intent(in) :: iV(n)
real(r4),intent(out) :: oV(n)
 do i=1,n
  oV(i)=real(iV(i))
 enddo
end subroutine

subroutine doubleMat(n,m,iM,oM)
use fiso, only: r8,r4
integer n,m,i,j
real(r4),intent(in) :: iM(n,m)
real(r8),intent(out) :: oM(n,m)
do j=1,m
 do i=1,n
  oM(i,j)=dble(iM(i,j))
 enddo
enddo
end subroutine

subroutine doubleVec(n,iV,oV)
use fiso, only: r8,r4
integer n,i
real(r4),intent(in) :: iV(n)
real(r8),intent(out) :: oV(n)
 do i=1,n
  oV(i)=dble(iV(i))
 enddo
end subroutine

! compute condition number of a given matrix(MxN) using SVD, M>N
subroutine get_cnr(m,n,matrix,cnr)
use fiso, only: r8
implicit none
integer m,n,lwork2,info
real(r8) cnr,qwork
real(r8) matrix(m,n),backup(m,n)
real(r8), allocatable :: aux2(:)
real(r8) S(m),U(m,m),VT(n,n)

backup=matrix

! get condition number of augmented Hessian
call DGESVD('S','S', m, n, matrix , m, S, U, m, VT, n, QWORK, -1, INFO )
lwork2=int(qwork)
allocate(aux2(lwork2))
call DGESVD('S','S', m, n, matrix ,m, S, U, m, VT, n, aux2, lwork2, INFO )
!write(*,*) (S(i),i=1,nvar1)
cnr=s(1)/s(m)
write(*,'(3x,a,E018.3,2x,F8.5,2x,E10.3)')'  --> CND: ',cnr, log(cnr),1d0/cnr

matrix=backup

end subroutine


! get unit vector
subroutine evec(a,b,e)
use fiso, only: r8
implicit none
real(r8), intent(in) :: a(3),b(3)
real(r8), intent(out) :: e(3)
real(r8) rab,dbond

rab=dbond(a,b);
e(1)=(a(1)-b(1))/rab
e(2)=(a(2)-b(2))/rab
e(3)=(a(3)-b(3))/rab

end subroutine

! distance between cartesian vectors a and b
real(r8) pure function dbond(a,b)
use fiso, only: r8
implicit none
real(r8), intent(in) :: a(3),b(3)
real(r8) ab(3)
ab=a-b
dbond=sqrt( ab(1)**2 + ab(2)**2 + ab(3)**2 )
end function


subroutine matmult(ct1,ct2,m,n,k,a,b,c)
! C(m,n) = A(m,k) x B(k,n)
use fiso, only: r8
character(1), intent(in)  :: ct1,ct2 ! "char temp"
character(1)  ch1,ch2                ! actual char for transpose indication
integer :: m,n,k,lda,ldb,ldc
real(r8), intent(in) :: a(*),b(*)
real(r8), intent(inout) :: c(*)
! real(r8), parameter  :: one=1d0, zero=0d0
real(r8), parameter  :: one=1.0_r8, zero=0.0_r8

ch1=ct1
ch2=ct2
call lower_case(ch1)
call lower_case(ch2)

lda=m
ldb=k
ldc=m
! neither transposed
if(ch1=='n'.and.ch2=='n') then ! tested
 call dgemm(ch1,ch2,m,n,k,one,a,lda,b,ldb,zero,c,ldc)
elseif(ch1=='n'.and.ch2=='t') then ! tested
 ldb=n
!  print*,m,n,k,lda,ldb,ldc
 call dgemm(ch1,ch2,m,n,k,one,a,lda,b,ldb,zero,c,ldc)
elseif(ch1=='t'.and.ch2=='n') then
 lda=k
 call dgemm(ch1,ch2,m,n,k,one,a,lda,b,ldb,zero,c,ldc)
elseif(ch1=='t'.and.ch2=='t') then
 if(m/=n) call error('Wrong DGEMM dimensions!')
 lda=k
 ldb=n ! =m
 call dgemm(ch1,ch2,m,n,k,one,a,lda,b,ldb,zero,c,ldc)
else
 call error('wrong transpose informations')
endif

end subroutine


subroutine screenMat(m,n,mat,thr,do_print,title)
! remove numerical noise from Matrix and tell sparsity
use fiso, only: r8,stdout
implicit none
integer i,j,m,n
real(r8), intent(inout) :: mat(m,n)
real(r8), intent(in) :: thr
logical, intent(in) :: do_print
character(*) title
real(r8) nzero,nelem,count

count=0
nelem=n*m
nzero=0_r8
! write(stdout,'(E10.3)') epsilon(1d0)

do j=1,n
 do i=1,m
  ! print*,mat(i,j)
  if(abs(mat(i,j))<thr) then
    ! print*,mat(i,j)
    mat(i,j)=0.0_r8
    nzero=nzero+1_r8
  endif
 enddo
enddo

if(do_print) then
 write(stdout,'(a,a,2x,F7.3)') 'sparsity of ',title, (nzero/nelem)*100_r8
endif

end subroutine
