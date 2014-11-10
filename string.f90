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

!********************************
!* convert a word to lower case *
!********************************
      subroutine lower_case(word)
      character (len=*) , intent(in out) :: word
      integer :: i,ic,nlen
      nlen = len(word)
      do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
      end do
      end subroutine lower_case


!********************************************
!* split a string s into n separate words w *
!********************************************
subroutine charsplit(s,n,w)
implicit none
integer i,n,k
character*80, intent(in) :: s
character*80, intent(out) :: w(n)
character*80   a,aa

aa=adjustl(s)
do i=1,n
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
enddo
return
end subroutine



!***************************************************
!* split a string s iand the return the x'ths word *
!***************************************************
subroutine charXsplit(s,wx,x)
implicit none
integer i,n,k,x
character*80, intent(in) :: s
character*80, intent(out) ::wx
character*80   w(20)
character*80   a,aa

aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
wx=w(x)
!print*, trim(wx)
return
end subroutine


!***********************************
!*  count words n in string s      *
!***********************************
subroutine cstring(s,n)
implicit none
integer i,n,k,x
character(*), intent(in) :: s
character*80   a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
!  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
n=i
return
end subroutine

integer function nwords(s)
implicit none
integer i,k,x
character(*), intent(in) :: s
character*80   a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
!  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
nwords=i
return
end function



! String to integer
integer pure function s2i(a)
implicit none
character(*), intent(in):: a
!read(a,'(i)') s2i
read(a,*) s2i
return
end function


! String to real(8)
real(8) pure function s2r(a)
implicit none
character(*), intent(in):: a
!read(a,'(f)') s2r
read(a,*) s2r
return
end function

