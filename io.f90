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

!******************
!* write xyz      *
!******************
subroutine wrxyz(infile)
use parm
implicit none
integer ierr
character(*) infile
character(2) esym
real(8) f

f=0.5291770d0
open(unit=55,file=infile,iostat=ierr,status='replace')
if(ierr.ne.0) stop 'cannot write xopt.xyz'
write(55,'(I5)') nat
write(55,'(2F16.8)')energy,gnorm
do i=1,nat
 write(55,'(a2,5x,3(F18.14,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(55)
end subroutine wrxyz 


!**********************
!* app xyz geoms      *
!**********************

subroutine appxyz(infile)
use parm
implicit none
character(*) infile
character(2) esym
real(8) f

f=0.5291770d0
open(unit=11,file=infile,status='old',position='append')
write(11,'(I5)') nat
write(11,'(2F16.8)') energy,gnorm
do i=1,nat
 write(11,'(a2,5x,3(F16.12,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(11)
end subroutine appxyz


!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd(infile,echo,c_nat)
use parm
use popt
implicit none
character*2 cc,ff
character*80  atmp
character*(*) infile
real(kind=8) txyz(3,maxat),xx(5)
real(kind=8) bohr
integer tiat(maxat),nn,tifrez(maxat),istat,iff
logical da,c_nat,echo
bohr=0.52917726d0
i=0
tifrez=0
iff=0

inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(echo) write(*,'('' reading...'',$)')

open(unit=33,file=infile)
! test for tmol or txyz file

 read(33,'(a)') atmp ! $coord
rewind(33)
if(index(atmp,'$coord').ne.0) then

 ! count number of atoms
 do while (da)
  read(33,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo
 nat=i
 100 continue
 if(c_nat) then
  close(33)
  return  ! just return number of atoms
 endif
 rewind(unit=33)

 ! read TMOL file
 read(33,*) atmp ! $coord
 do j=1,nat
    read(33,'(a)') atmp ! $coord
    backspace(33)
    if(index(atmp,' f ').ne.0) then
    read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc,ff
     tifrez(j)=1
    iff=iff+1
    else ! default
     read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc
   endif
   call elem(cc,tiat(j))
   txyz(1:3,j)=txyz(1:3,j)
  enddo
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(33)

else ! txyz file
       read(33,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do
            nat=nat+1
            read(33,'(a)',end=123) atmp
           enddo
            if(c_nat) then
              close(33)
              return  ! just return number of atoms
            endif
          else
            nat=idint(xx(1))
            if(c_nat) then
             close(33)
             return  ! just return number of atoms
            endif
           read(33,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(33,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,tiat(i))
            txyz(1:3,i)=xx(1:3)*1d0/bohr
       enddo
 101  close(33)
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif


if(maxval(tifrez,nat).eq.1) then
freeze=.true.
 if(echo) then
  write(*,'(a,x,I4,x,a)') '  found ',iff,' frozen cart. coordinates'
  if(iff.lt.50) then ! dont spam the output to much ...
   write(*,'($,a)') '  atom nr: '
   do i=1,nat
     if(tifrez(i)==1) write(*,'($,x,I2)') i
   enddo
   print*,''
  endif
 endif
endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select

do i=1,nat
xyz(1:3,i)=txyz(1:3,i)
iat(i)=tiat(i)
if(freeze) ifrez(i)=tifrez(i)
enddo
return

end subroutine


      SUBROUTINE READL(A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) A1
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=READAA(A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END


      FUNCTION READAA(A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 READAA
      CHARACTER*(*) A
      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')
      IEND=0
      IEND2=0
      IDIG=0
      C1=0
      C2=0
      ONE=1.D0
      X = 1.D0
      NL=LEN(A)
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J:J))
         M=ICHAR(A(J+1:J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20

   10 CONTINUE
      READAA=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1.D0
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10+N-IZERO
            X = X /10
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE
!C
!C PUT THE PIECES TOGETHER
!C
   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO 31 I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
   31 CONTINUE
   61 CONTINUE
   70 READAA=READAA*10**(ONE*C1)
      RETURN
      END
 SUBROUTINE ELEM(KEY1, NAT)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 CHARACTER*(*) KEY1
 CHARACTER*2 ELEMNT(94),E

 DATA ELEMNT/'h ','he',                                      &
 'li','be','b ','c ','n ','o ','f ','ne',                    &
 'na','mg','al','si','p ','s ','cl','ar',                    &
 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',     &
 'zn','ga','ge','as','se','br','kr',                         &
 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',     &
 'cd','in','sn','sb','te','i ','xe',                         &
 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
 'au','hg','tl','pb','bi','po','at','rn',                    &
 'fr','ra','ac','th','pa','u ','np','pu'/

 nat=0
 e='  '
 k=1
 DO J=1,len(key1)
    if (k.gt.2)exit
    N=ICHAR(key1(J:J))
    if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
       e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
       k=k+1
    endif
    if(n.ge.ichar('a') .and. n.le.ichar('z') )then
       e(k:k)=key1(j:j)
       k=k+1
    endif
 enddo
 DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

character*2 FUNCTION ESYM(I)
CHARACTER*2 ELEMNT(94)
DATA ELEMNT/'h ','he',                                           &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu'/
  ESYM=ELEMNT(I)
  RETURN
END


      subroutine ordhess(nat3,h,fname)
      implicit none
      integer nat3
      real*8 h(nat3,nat3)
      character*(*) fname
      integer iunit,i,j,mincol,maxcol
      character*5 adum
      character*80 a80

      write(*,*)
      write(*,*) 'Reading <',trim(fname),'>'
      iunit=1
      open(unit=iunit,file=fname)
   50 read(iunit,'(a)')a80
      if(index(a80,'$hessian').ne.0)then
      do 100 i=1,nat3
        maxcol = 0
  200   mincol = maxcol + 1
        maxcol = min(maxcol+5,nat3)
        read(iunit,'(a5   ,5f15.10)')adum ,(h(i,j),j=mincol,maxcol)
        if (maxcol.lt.nat3) goto 200
  100 continue
      close(iunit,status='delete')
      goto 300
      endif
      goto 50

  300 return
      end



subroutine rdhess(nat3,h,fname)
implicit none
integer nat3
real*8 h(nat3,nat3)
character*(*) fname
integer iunit,i,j,mincol,maxcol
character*5 adum
character*80 a80

write(*,*)
write(*,*) 'Reading <',trim(fname),'>'
iunit=11
open(unit=iunit,file=fname)

do
read(iunit,'(a)',end=201) a80
  if(index(a80,'$hessian').ne.0)then
   do  i=1,nat3
      maxcol = 0
 200   mincol = maxcol + 1
       maxcol = min(maxcol+5,nat3)
       read(iunit,'(a5   ,5f15.10)')adum ,(h(i,j),j=mincol,maxcol)
       if (maxcol.lt.nat3) goto 200
   enddo
   exit
 endif
!print*,a80
enddo
 close(iunit,status='keep')
  300 return
201 write(*,*) ' no TM hessian found'
 close(iunit,status='keep')

end

subroutine wrhess(nat3,h,fname)
implicit none
integer nat3
real*8 h(nat3,nat3)
character*(*) fname
integer iunit,i,j,mincol,maxcol
character*5 adum
character*80 a80

adum='   '
iunit=11
open(unit=iunit,file=fname)
a80='$hessian'
write(iunit,'(a)')a80
do i=1,nat3
  maxcol = 0
  do while (maxcol.lt.nat3)
    mincol = maxcol + 1
    maxcol = min(maxcol+5,nat3)
    write(iunit,'(a5,5f15.10)')adum,(h(i,j),j=mincol,maxcol)
  enddo !if (maxcol.lt.nat3) goto 200
enddo
write(iunit,'(''$end'')')
close(iunit)
end


subroutine readORCAhess(n3,h,filen)
implicit none
real(8) h(n3,n3),x,xx(10)
integer i,j,n3,k,d,c,l
character*(*) filen
character*80 a
logical da

inquire(file=filen,exist=da)
if(.not.da) print*,'cannot read ',filen,' (ERROR)'
open(unit=33,file=filen)
c=mod(n3,6)
  do
     read(33,'(a)',end=201)a
     if(index(a,'$hessian').ne.0) then
     read(33,*,end=201) d ! dimension
     if(d.ne.n3) stop 'wrong hessian dimension!'
     j=0
     do k=1,int(n3/6)
         read(33,'(a)',end=201)a ! legend
         do i=1,n3
           read(33,*) x,(h(i,l),l=1+j,6+j)
!        write(*,'(7f7.4)') x,(h(i,l),l=1+j,6+j)
         enddo
      j=j+6
      enddo
! rest
          read(33,'(a)',end=201)a ! legend
         do i=1,n3
          read(33,*) x,(h(i,l),l=j+1,j+c)
!          write(*,'(4f7.4)') x,(h(i,l),l=1+j,j+c)
         enddo
     endif
   enddo
201  continue
!write(*,*) h(1,1),h(7,2),h(n3,n3)
close(33)

end subroutine readORCAhess



subroutine wtm(outfile)
use parm
implicit none
real*8 bohr
character(*) outfile
character(2) esym
real(8) f

bohr=1
! bohr=0.52917726
open(unit=44,file=outfile)
write(44,'(a)')'$coord'
do i=1,nat
  write(44,'(3F18.12,2x,a2)') xyz(1:3,i)/bohr , esym(iat(i))
enddo
write(44,'(a)')'$end'
close(44)
end subroutine wtm


subroutine error(io,aa)
implicit none
integer io
character(*) aa


write(io,'(a)')
write(io,'(3x,a)')'************  ERROR  *********************'
write(io,'(4x,a)') trim(aa)
write(io,'(3x,a)')'******************************************'
write(io,'(a)')
stop

end subroutine error





subroutine readamber(nat,iat,xyz,pc,ncrd,ntop)
implicit none
integer i,j,k,l,nat,iat(nat)
real(8) xyz(3,nat),vxyz(3*nat),mass(nat),pc(nat)
integer ierr
character(*) ncrd,ntop  
character(80) F1,aa
character(2) esym
real(8) f


! read prmtop file (ntop)
print*,' reading ',trim(ntop)
iat=0
open(123,file=trim(ntop))
do
  read(123,'(a)',end=777) aa
    if(index(aa,'%FLAG ATOMIC_NUMBER').ne.0) then
     read(123,'(a)') aa !dummy
        read(123,*) iat(1:nat)
    endif
    if(index(aa,'%FLAG CHARGE').ne.0) then
     read(123,'(a)') aa !dummy
        read(123,*) pc(1:nat)
    endif
    
! sometimes ???? the FLAG ATOMIC NUMBER is done there. Maybe use the masses flag instead?
enddo
777 continue
close(123)

! from Amber documentation
!  CHARGE : the atom charges.  (Divide by 18.2223 to convert to charge
!           in units of the electron charge)

pc=pc/18.2223
print*,pc(1:5)

!if(iat(1)==0) stop 'something went wrong in readamber'

! use masses .... 
if(iat(1)==0) then
print*,'flag %FLAG ATOMIC_NUMBER not found. Guessing from masses'
open(123,file=trim(F1))
do
  read(123,'(a)',end=771) aa
    if(index(aa,'%FLAG MASS').ne.0) then
     read(123,'(a)') aa
        read(123,*) mass(1:nat)
      exit
    endif
enddo
771 continue
close(123)

do i=1,nat
if(mass(i)==1.008d0) iat(i)=1
if(mass(i)==16.d0) iat(i)=8
if(mass(i)==12.01d0) iat(i)=6
if(mass(i)==12.001d0) iat(i)=6
if(mass(i)==14.01d0) iat(i)=7
if(mass(i)==30.97d0) iat(i)=15
if(iat(i)==0) then
 print*,'Mass', mass(i)
 stop'undefined mass encoutered. stopping'
endif
enddo

endif



! read crd file (ncrd)


f=0.5291770d0
open(222,file=ncrd)
read(222,'(a)') aa ! comment
read(222,'(a)') aa ! number of atoms
read(222,*) xyz(1:3,1:nat)
close(222)
xyz=xyz/f
end subroutine
