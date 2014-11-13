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
!**
!* take updated, approx. Hessian and form from this approx. normal coordinates
subroutine getanc
use parm
use popt
use logic
implicit none
integer ij
real(8), allocatable :: e(:),h(:,:),totsym(:),tr(:)
real(8) hmax,damp,edum
real(8) D(nat3*(nat3+1)/2)

print*,' Making approx. normal coords'
allocate(e(nat3))
allocate(totsym(nat3))
allocate(tr(nat3*(nat3+1)/2))
allocate(h(nat3,nat3))
Hmin=0.001d0
Hmax=1e+5

! Projecting out Tr/Rot
k=0
do i=1,nat3
   do j=1,i
      k=k+1
      tr(k)=0.5*(chess(i,j)+chess(j,i))
      if(i.ne.j .and. abs(tr(k)).lt.1d-10) tr(k)=0.0d0
   enddo
enddo
call TRpr(nat,xyz0,tr)

! diagonalize
k=0
do i=1,nat3
   do j=1,i
      k=k+1
      h(i,j)=tr(k)
      h(j,i)=tr(k)
   enddo
enddo
call DiagSM(nat3,h,e)

call findlroot(nat3,e,edum)

! 5 to 6 modes should be zero, leave them out
totsym=0
do i=1,nat3
   if(abs(e(i)).gt.1.d-8) totsym(i)=1
!   if(abs(e(i)).gt.1.d-8) then
!     totsym(i)=1
!   else
!    e(i)=0.0d0
!   endif
enddo
nvar=sum(totsym)

write(*,'(3x,a,I4)') 'Nvar(c1) : ',nvar

if(Doshift) then
  ! find shift
!  call findlroot(nat3,e,edum)

  ! shift
  if(.not.doTS) then
     damp = hmin - edum
     if(damp.lt.0) damp = 0
     do i=1,nat3
!     if(abs(e(i)).gt.1.d-8)then
     if(abs(e(i)).gt.1.d-8.and.abs(e(i)).lt.1.0)then
        e(i) = e(i) + damp
     endif
     enddo
  endif
   write(*,*) 'Shifting diagonal of input Hessian by ', damp
   write(*,*) 'Lowest  eigenvalues of approx. Hessian'
   write(*,'(6F13.8)')(e(i),i=1,min(36,nat3))
Doshift=.false.
endif

allocate(b(nat3,nvar))
allocate(hint(nvar*(nvar+1)/2))

b = 0
k = 0
hint = 0
do i=nat3,1,-1
!if(abs(e(i)).gt.1.d-8 .and. k.le.nvar  .and.  totsym(i).eq.1)then
if(abs(e(i)).gt.1.d-8 .and.k.le.nvar)then
  k=k+1
  b(1:nat3,k)=h(1:nat3,i)
      if(.not.doTS)then
         hint(k+k*(k-1)/2)=max(e(i),hmin)
         if(hint(k+k*(k-1)/2) > hmax) hint(k+k*(k-1)/2)=hmax
      else
         hint(k+k*(k-1)/2)=e(i)
      endif
   endif
enddo

if(k/=nvar) then
 print*,nvar,k
 stop 'k/=nvar! error'
endif


if(freeze) call freezeBmat
call sort(nat3,nvar,hint,B)



deallocate(e,totsym)
deallocate(h,tr)

end subroutine



subroutine gxyz2gint
use parm
implicit none
integer m

gint = 0
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
use parm
implicit none
integer xvar,m
real(8) xx,cart(3,nat),cint(xvar)

cart = 0.0d0
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

subroutine DiagSM(xvar,mat,eig)
implicit none
integer i,j,k
real(8), allocatable :: aux(:)
integer info,lwork,xvar
real(8) ,intent(in) :: mat(xvar,xvar)
real(8) xx
real(8), intent(out) :: eig(xvar)

eig=0
call dsyev ('V','U',xvar,mat,xvar,eig,xx,-1,info)
lwork=int(xx)
allocate(aux(lwork))
call dsyev ('V','U',xvar,mat,xvar,eig,aux,lwork,info)
if(info/=0) print*,'Diagonalization failed !!'

end subroutine

! Hxyz=Bt*H*B
subroutine hint2xyz(nvar,nat3,hint,chess,b)
implicit none
real*8 hint(nvar*(nvar+1)/2),b(nat3,nvar)
real*8 chess(nat3,nat3)
integer nvar,nat3
integer i,j,k,padr
real*8 :: bh(nvar,nat3),xx
real(8) h(nvar,nvar)
!                   M    N    K
!             M,K * K,N 
call packM(nvar,h,hint,'unpack')
call dgemm('N','T',nvar,nat3,nvar,1d0,h,nvar,b,nat3,0d0,bh,nvar)
call dgemm('N','N',nat3,nat3,nvar,1d0,b,nat3,bh,nvar,0d0,chess,nat3)
end subroutine




subroutine findlroot(n,e,el)
implicit none
integer n,i
real*8 e(n), el
el=1.d+42
do i=1,n
   if(abs(e(i)) .gt. 1.d-10 .and. e(i) .lt. el ) el = e(i)
enddo
end



subroutine bsort(n,e)
! bubble sort the vector e
implicit none
integer i,j,k,l,nn,n,ii
real(8) e(n),tt
character(80) cc
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

return
end subroutine




! address in packed array
integer function padr(i1,i2)
integer i1,i2,idum1,idum2
idum1=max(i1,i2)
idum2=min(i1,i2)
padr=idum2+idum1*(idum1-1)/2
return
end

subroutine sort(nat3,nvar,hess,b)
implicit none
integer nat3,nvar,ii,k,j,m,i
real*8 hess(nvar*(nvar+1)/2),b(nat3,nat3),pp,sc1
real*8 :: edum(nvar)

   do k=1,nvar
      edum(k)=hess(k+k*(k-1)/2)
   enddo
!c sort
 do 140   ii = 2, nvar
    i = ii - 1
    k = i
    pp= edum(i)
    do 120   j = ii, nvar
       if (edum(j) .gt. pp) go to 120
       k = j
       pp= edum(j)
  120    continue
     if (k .eq. i) go to 140
     edum(k) = edum(i)
     edum(i) = pp
     do m=1,nat3
        sc1=b(m,i)
        b(m,i)=b(m,k)
        b(m,k)=sc1
     enddo
  140 continue

do k=1,nvar
   hess(k+k*(k-1)/2)=edum(k)
enddo

end




