subroutine getanc
! take Hessian and make normal coordinates
use fiso, only: r8,stdout
use parm
use popt
use logic
use timing, only: timer
implicit none
integer pmax
real(r8), allocatable :: e(:),h(:,:),tr(:)
real(r8) hmax,damp,edum,freqcm
! real(8) D(nat3*(nat3+1)/2)

real(r8) time

call status1(timer)
call message_head('Generating normal coordinates')
allocate(e(nat3))
allocate(tr(nat3*(nat3+1)/2))
allocate(h(nat3,nat3))
Hmin=0.001_r8
Hmax=1e+5_r8
e=0_r8
tr=0_r8
h=0_r8


if(nat>4) then
   pmax=6
else
  pmax=int(nat3/2)
endif


if(debug) call debug1('Tr/Rot',time)
! Projecting out Tr/Rot
k=0
do i=1,nat3
   do j=1,i
      k=k+1
      tr(k)=0.5_r8*(chess(i,j)+chess(j,i))
      if(i.ne.j .and. abs(tr(k)).lt.1e-10_r8) tr(k)=0.0_r8
   enddo
enddo

! skip T/R projection for 2 atom molecules (bug work around)
if(nat>2) call TRpr(nat,xyz0,tr)
if(debug) call debug2(time)

! diagonalize
k=0
do i=1,nat3
   do j=1,i
      k=k+1
      h(i,j)=tr(k)
      h(j,i)=tr(k)
   enddo
enddo

! mass-weight 
! if(tsopt) call Hmass(h,'noproject')

if(debug) call debug1('DIAG',time)
! if(nat>2000) then
  ! call DiagSM3(nat3,h,e,1e-10_r8)
! else
  call DiagSM(nat3,h,e)
! endif
if(debug) call debug2(time)

call findlroot(nat3,e,edum)
! leave out trans/rot modes
!print*, e(1:10)
nvar=0
do i=1,nat3
   if(abs(e(i)).gt.1.e-12_r8) then
     nvar=nvar+1
   else
    e(i)=0.0_r8
   endif
enddo


call status2(timer)
write(stdout,'(2x,a,I5)') 'degrees of freedom(anc) : ',nvar
write(stdout,'(a)') ''
if(nat>2) then
 if(nvar/=nat3-6) then
  call exclaim ('possible incomplete Hessian!')
  write(*,'(2x,a,I5)') 'expected: degrees of freedom(3N-6) : ',nat3-6
  write(*,'(a)') ''
 endif
endif 

! remove negative eigenvalues by shifting
if(Doshift) then
  damp=0.0_r8
  if(.not.tsopt) then
     damp = hmin - edum
     if(damp.lt.0) damp = 0.0_r8
     do i=1,nat3
     if(abs(e(i)).gt.1.e-8_r8.and.abs(e(i)).lt.1.0_r8)then
        e(i) = e(i) + damp
     endif
     enddo
  endif
   write(*,'(2x,a,F8.2)') 'Hessian shift (cm-1) :', freqcm(damp)
   write(*,'(2x,a)') 'lowest projected vib. freq (cm-1) of non mass-weighted input Hessian after shift:'
   call printvib(nat3,e)
Doshift=.false.
endif

! Now make Bmatrix

allocate(b(nat3,nvar))
allocate(hint(nvar*(nvar+1)/2))

b = 0_r8
k = 0
hint = 0_r8
do i=nat3,1,-1
if(abs(e(i)).gt.1.e-8_r8 .and.k.le.nvar)then
  k=k+1
  b(1:nat3,k)=h(1:nat3,i)
      if(.not.tsopt)then
        !  print*,k,k+k*(k-1)/2,i
         hint(k+k*(k-1)/2)=max(e(i),hmin)
         if(hint(k+k*(k-1)/2) > hmax) hint(k+k*(k-1)/2)=hmax
      else
         hint(k+k*(k-1)/2)=e(i)
      endif
   endif
enddo

if(k/=nvar) then
 print*,nvar,k
 call error('k/=nvar! error')
endif


if(freeze) call freezeBmat

call sort(nat3,nvar,hint,B)

deallocate(e)
deallocate(h,tr)

end subroutine



subroutine gxyz2gint
! g*B*gT
use fiso, only: r8
use parm,only: nat,grad,gint,nvar,b
implicit none
integer i,j,k,l

gint = 0_r8
do l=1,nvar
 k=0
  do i=1,nat
    do j=1,3
      k=k+1
      gint(l)=gint(l)+b(k,l)*grad(j,i)
    enddo
  enddo
enddo
end subroutine gxyz2gint


subroutine int2xyz(xvar,cint,cart)
! colinear transformation of internals to cartesian
use parm
use fiso, only: r8
implicit none
integer xvar,m
real(r8) xx,cart(3,nat),cint(xvar)

cart = 0.0_r8
do i=1,xvar
   m=0
   do j=1,nat
      do k=1,3
         m=m+1
         xx=b(m,i)*cint(i)
         cart(k,j)=cart(k,j)+xx
      enddo
   enddo
enddo
cart = cart + xyz0
end subroutine int2xyz

subroutine hint2xyz(nvar,nat3,hint,chess,b)
! Hxyz=Bt*H*B
! HERE B is (nat3,nvar)
use fiso, only: r8
implicit none
integer :: nvar,nat3
real(r8) :: hint(nvar*(nvar+1)/2),b(nat3,nvar)
real(r8) :: chess(nat3,nat3)
real(r8) :: bh(nvar,nat3)
real(r8) :: h(nvar,nvar)
!  M    N    K
!  M,K * K,N
call packM(nvar,h,hint,'unpack')
! print*,h(1:5,1),b(1:5,1)
! stop
! call dgemm('N','T',nvar,nat3,nvar,1d0,h,nvar,b,nat3,0d0,bh,nvar)
! call dgemm('N','N',nat3,nat3,nvar,1d0,b,nat3,bh,nvar,0d0,chess,nat3)
! print*,chess(1:3,1)
call matmult('n','t',nvar,nat3,nvar,h,b,bh)
call matmult('n','n',nat3,nat3,nvar,b,bh,chess)
! print*,chess(1:3,1)
end subroutine


subroutine findlroot(n,e,el)
! find lowest, positive value in vector
use fiso, only: r8
implicit none
integer n,i
real(r8) e(n), el
el=HUGE(1_r8)
do i=1,n
   if(abs(e(i)) .gt. 1.e-10_r8 .and. e(i) .lt. el ) el = e(i)
enddo
end

subroutine bsort(n,e)
! bubble sort the vector e
use fiso, only: r8
implicit none
integer i,nn,n
real(r8) e(n),tt
logical order

nn=n
order=.false.
do
if(order) exit
order=.true.
 do i=1,nn-1
    if (e(i).gt.e(i+1) ) then ! swap
      tt=e(i)
      e(i)=e(i+1)
      e(i+1)=tt
      order = .false.
     endif
 enddo
nn=nn-1
enddo

end subroutine


integer function padr(i1,i2)
! address of the diagonal elements of a
! triagonal matrix (eg 1,3,6,10, etc)
integer i1,i2,idum1,idum2
 idum1=max(i1,i2)
 idum2=min(i1,i2)
 padr=idum2+idum1*(idum1-1)/2
 return
end

subroutine sort(nat3,nvar,hess,b)
! sort bmatrix
use fiso, only: r8
implicit none
integer, intent(in) :: nat3,nvar
integer :: ii,k,j,m,i
real(r8) hess(nvar*(nvar+1)/2),b(nat3,nat3),pp,sc1
real(r8) edum(nvar)

do k=1,nvar
    edum(k)=hess(k+k*(k-1)/2)
enddo
!c sort
 do ii = 2, nvar
       i = ii - 1
       k = i
       pp= edum(i)
       do j = ii, nvar
          if (edum(j) .gt. pp) cycle
          k = j
          pp= edum(j)
       enddo
       if (k .eq. i) cycle
       edum(k) = edum(i)
       edum(i) = pp
       do m=1,nat3
          sc1=b(m,i)
          b(m,i)=b(m,k)
          b(m,k)=sc1
       enddo
 enddo

 do k=1,nvar
    hess(k+k*(k-1)/2)=edum(k)
 enddo

end subroutine




real(r8) function freqcm(e)
! convert frequency from au to cm
use fiso, only: r8
use constant, only : amu2au,au2cm
implicit none
real(r8) e
freqcm=au2cm*sign(sqrt(abs(e)),e)/sqrt(amu2au)
end function



