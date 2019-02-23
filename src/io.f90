
!******************
!* write xyz      *
!******************
subroutine wrxyz(infile)
use fiso, only: r8
use constant, only: au2ang
use parm, only: xyz,i,energy,gnorm,nat,iat
implicit none
integer ierr,io
character(*) infile
character(2) esym
real(r8) f

f=au2ang
open(newunit=io,file=infile,iostat=ierr,status='replace')
if(ierr.ne.0) call error('cannot write xopt.xyz')
write(io,'(I5)') nat
write(io,'(2F16.8)')energy,gnorm
do i=1,nat
 write(io,'(a2,5x,3(F18.14,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(io)
end subroutine wrxyz


!*************************
!* write 'local xyz      *
!*************************
subroutine localxyz(infile,n,ia,nxyz)
use constant, only: au2ang
use fiso, only: r8
implicit none
integer ierr,n,ia(n),i,io
real(r8) nxyz(3,n)
character(*) infile
character(2) esym
real(r8) f

f=au2ang
open(newunit=io,file=infile,iostat=ierr,status='replace')
if(ierr.ne.0) call error('cannot write xopt.xyz')
write(io,'(I5)') n
write(io,'(a)') 'xopt xyz file'
do i=1,n
 write(io,'(a2,5x,3(F18.14,3x))') esym(ia(i)), nxyz(1,i)*f,nxyz(2,i)*f,nxyz(3,i)*f
enddo
close(io)
end subroutine localxyz


!*************************
!* read local xyz        *
!*************************
subroutine read_localcoord(infile,n,ia,nxyz,ftype)
use constant, only: au2ang
use fiso, only: r8
implicit none
integer ierr,n,ia(n),i,io
real(r8) nxyz(3,n)
character(*) infile,ftype
character(120) aa
character(2) el
real(r8) f
logical fstr

if(ftype=='xyz') then
  f=au2ang
  open(newunit=io,file=infile,iostat=ierr,status='old')
  if(ierr.ne.0) call error('error reading '//trim(infile))
  read(io,'(a)') aa !dummy
  read(io,'(a)') aa !dummy
  do i=1,n
   read(io,*) el, nxyz(1,i),nxyz(2,i),nxyz(3,i)
   call elem(el,ia(i))
  enddo
  close(io)
  nxyz=nxyz/f ! bohrs
elseif(ftype=='tm') then
  open(newunit=io,file=infile,iostat=ierr,status='old')
  if(ierr.ne.0) call error('error reading '//trim(infile))
  do
    read(io,'(a)') aa
    if(fstr(aa,'$coord')) then
      do i=1,n
       read(io,*) nxyz(1,i),nxyz(2,i),nxyz(3,i),el
       call elem(el,ia(i))
      enddo
      exit
    endif
  enddo
  close(io)
else
 call error('no file type given!')
endif
end subroutine read_localcoord

!******************************
!* find nat in tmole/xyz file *
!******************************
integer function tm_nat(infile,ftype)
use fiso, only: stdout
implicit none
character(*) infile,ftype
character(120) aa
integer io,n,iio
logical fstr

n=0

if(ftype=='tm') then
  open(newunit=io,file=infile,status='old')
  do
   read(io,'(a)',iostat=iio) aa
   if(iio/=0) exit
   if(fstr(aa,'$coord')) then
     do
       read(io,'(a)',iostat=iio) aa
       if(fstr(aa,'$')) exit
       if(iio/=0) exit
       n=n+1
     enddo
   endif
   if(fstr(aa,'$')) exit
  enddo
  close(io)
elseif(ftype=='xyz') then
  open(newunit=io,file=infile,iostat=iio,status='old')
  if(iio/=0) call error('error reading '//trim(infile))
  read(io,*) n
  close(io)
endif

write(stdout,*) '  Found ',n,' atoms in',trim(infile)
if(n==0) call error('no coordinates found! <function tm_nat>')
tm_nat=n
end function

!**********************
!* app xyz geoms      *
!**********************

subroutine appxyz(infile)
use parm,only: energy, gnorm,xyz,i,nat,iat
use constant, only: au2ang
use fiso, only: r8
implicit none
integer :: io
character(*), intent(in) :: infile
character(2) :: esym
real(r8) :: f

f=au2ang
open(newunit=io,file=infile,status='old',position='append')
write(io,'(I5)') nat
write(io,'(2F16.8)') energy,gnorm
do i=1,nat
 write(io,'(a2,5x,3(F16.12,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(io)
end subroutine appxyz


!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd(infile,echo,c_nat)
use parm,     only: nat,iat,xyz,i,j,maxat,ifrez
use popt,     only: freeze
use constant, only: au2ang
use fiso, only: r8, stdout
implicit none
character(2) cc,ff
character(255)  atmp,wx
character(*) infile
real(r8) txyz(3,maxat),xx(5)
real(r8) bohr
integer tiat(maxat),nn,tifrez(maxat),iff,s2i,io
logical da,c_nat,echo
bohr=au2ang
i=0
tifrez=0
iff=0

inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(echo) write(stdout,'('' reading...'')',advance="no")

open(newunit=io,file=infile)
! test for tmol or txyz file

read(io,'(a)') atmp ! $coord
rewind(io)
if(index(atmp,'$coord').ne.0) then
 ! count number of atoms
 do while (da)
  read(io,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo
 nat=i
 100 continue
 if(c_nat) then
  close(io)
  return  ! just return number of atoms
 endif
 rewind(unit=io)

 ! read TMOL file
 read(io,*) atmp ! $coord
 do j=1,nat
    read(io,'(a)') atmp ! $coord
    backspace(io)
    if(index(trim(atmp),' f ').ne.0) then
    read(io,*) txyz(1,j),txyz(2,j),txyz(3,j),cc,ff
     tifrez(j)=1
     iff=iff+1
    else ! default
     read(io,*) txyz(1,j),txyz(2,j),txyz(3,j),cc
   endif
   call elem(cc,tiat(j))
   txyz(1:3,j)=txyz(1:3,j)
  enddo
 if(echo) write(stdout,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(io)

else ! txyz file
       read(io,'(a)',end=101) atmp
      ! print*,atmp
! check for first two lines
      !  call readl(atmp,xx,nn)
        call cstring(atmp,nn)

        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do
            nat=nat+1
            read(io,'(a)',end=123) atmp
           enddo
            if(c_nat) then
              close(io)
              return  ! just return number of atoms
            endif
        else
           call charXsplit(atmp,wx,1)
            nat=s2i(wx)
            if(c_nat) then
             close(io)
             return  ! just return number of atoms
            endif
           read(io,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
       read(io,*) cc, xx(1),xx(2),xx(3)
!            read(33,'(a)') atmp
!            call readl(atmp,xx,nn)
!            call elem(atmp,tiat(i))
            call elem(cc,tiat(i))
            txyz(1:3,i)=xx(1:3)/au2ang
       enddo
 101  close(io)
      if(echo) write(stdout,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif


if(maxval(tifrez).eq.1) then
freeze=.true.
 if(echo) then
  write(stdout,'(a,1x,I4,1x,a)') '  found ',iff,' frozen cart. coordinates'
  if(iff.lt.50) then ! dont spam the output to much ...
   write(stdout,'(a)',advance="no") '  atom nr: '
   do i=1,nat
     if(tifrez(i)==1) write(stdout,'(1x,I2)',advance="no") i
   enddo
   write(stdout,*)''
  endif
 endif
endif

case (.false.)
  call error(' no input file <'//trim(infile) //'> found !! ')
end select

do i=1,nat
  xyz(1:3,i)=txyz(1:3,i)
  iat(i)=tiat(i)
  if(freeze) ifrez(i)=tifrez(i)
enddo

end subroutine

subroutine elem(key1, nat)
implicit none
!implicit double precision (a-h,o-z)
integer :: nat,k,n,j,i
character(*) key1
character(2) elemnt(95),e

data elemnt/'h ','he',                                      &
'li','be','b ','c ','n ','o ','f ','ne',                    &
'na','mg','al','si','p ','s ','cl','ar',                    &
'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',     &
'zn','ga','ge','as','se','br','kr',                         &
'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',     &
'cd','in','sn','sb','te','i ','xe',                         &
'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
'au','hg','tl','pb','bi','po','at','rn',                    &
'fr','ra','ac','th','pa','u ','np','pu','xx'/

nat=0
e='  '
k=1
do j=1,len(key1)
   if (k.gt.2)exit
   n=ichar(key1(j:j))
   if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
      e(k:k)=char(n+ichar('a')-ichar('A'))
      k=k+1
   endif
   if(n.ge.ichar('a') .and. n.le.ichar('z') )then
      e(k:k)=key1(j:j)
      k=k+1
   endif
enddo
do i=1,107
    if(e.eq.elemnt(i))then
       nat=i
       return
    endif
enddo

end subroutine

character(2) FUNCTION ESYM(I)
implicit none
integer :: i
! returns element symbol
CHARACTER(2) ELEMNT(95)
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
  'fr','ra','ac','th','pa','u ','np','pu','xx'/
  ESYM=ELEMNT(I)
  RETURN
END


subroutine rdhess(nat3,h,fname)
! read TM hessian from fname
use fiso, only: r8, stdout
implicit none
integer nat3
real(r8) h(nat3,nat3)
character(*) fname
integer iunit,i,j,mincol,maxcol,io
character(5) adum
character(255) a255
io=stdout

write(io,'(a)') ''
write(io,'(a)') 'Reading <'//trim(fname)//'>'
open(newunit=iunit,file=fname)
do
read(iunit,'(a)',end=201) a255
  if(index(a255,'$hessian').ne.0)then
   do  i=1,nat3
       maxcol = 0
 200   mincol = maxcol + 1
       maxcol = min(maxcol+5,nat3)
       read(iunit,'(a5   ,5f15.10)')adum ,(h(i,j),j=mincol,maxcol)
       if (maxcol.lt.nat3) goto 200
   enddo
   exit
  endif
enddo
close(iunit,status='keep')
return

201 write(io,'(a)') ' no TM hessian found'
 close(iunit,status='keep')

end

! TM style hessian
subroutine wrhess(nat3,h,fname)
use fiso, only: r8
implicit none
integer nat3
real(r8) h(nat3,nat3)
character(*) fname
integer iunit,i,j,mincol,maxcol
character(5) adum
character(80) a80

adum='   '
open(newunit=iunit,file=fname)
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


subroutine writebin(n,X,fname)
! write unformated sym. matrix X(n,n)
use fiso, only: r8
implicit none
integer, intent(in)  :: n
real(r8), intent(in) :: X(n,n)
character(*)         :: fname
integer              :: iunit
open(newunit=iunit,file=fname,form='unformatted',access='stream')
 write(iunit) X
close(iunit)
end

subroutine readbin(n,X,fname)
! read unformated sym. matrix X(n,n)
use fiso, only: r8,stdout
implicit none
integer, intent(in) ::  n
real(r8), intent(inout) :: X(n,n)
character(*) fname
integer iunit
logical da
write(stdout,*)' Reading ... ',trim(fname)
da=.false.
inquire(file=trim(fname),exist=da)
X=0.0d0
if(da) then
  open(newunit=iunit,file=fname,form='unformatted',access='stream')
   read(iunit) X
  close(iunit)
else
  write(stdout,*) trim(fname) ,'not found !'
endif
end

subroutine readPSI4hess(n3,h,filen)
!
! Also works for CFOUR hessians
!
use fiso, only: r8, stdout
implicit none
integer k,i,j,n3,io,nat,x
real(r8) h(n3,n3),vec(n3*n3)
character(*) filen

write(stdout,*) 'Reading <',trim(filen),'>'
open(newunit=io,file=filen)
  read(io,*) nat,x
  if(nat*3/=n3) call error('wrong hessian dimension!')
  k=1
  do i=1,n3*n3/3
    read(io,*) vec(k),vec(k+1),vec(k+2)
    k=k+3
  enddo
  k=0
  do i=1,n3
   do j=1,n3
    k=k+1
    h(i,j)=vec(k)
   enddo
  enddo
close(io)
call printmat(6,n3,n3,h,'H(PSI4)')
end subroutine


subroutine readORCAhess(n3,h,filen)
use fiso, only: r8,stdout
implicit none
integer i,j,n3,k,d,c,l,nval
real(r8) h(n3,n3),x
character(*) filen
character(80) a
logical da

inquire(file=filen,exist=da)
if(.not.da) write(stdout,*) 'cannot read ',filen,' (ERROR)'

write(stdout,*) 'Reading <',trim(filen),'>'
nval=6 ! 
open(unit=33,file=filen)
c=mod(n3,6)
  do
     read(33,'(a)',end=201)a
     if(index(a,'$hessian').ne.0) then
     read(33,*,end=201) d ! dimension
     ! print*, d,n3
     if(d.ne.n3) stop 'wrong hessian dimension!'
     ! check orca30 vs orca40 hessian
     read(33,'(a)') a
     read(33,'(a)') a
     print*, trim(a),len(trim(a))
     if(len(trim(a))>77) nval=5
     backspace(33)
     backspace(33)
     j=0
     do k=1,int(n3/nval)
         read(33,'(a)',end=201) a ! legend
         do i=1,n3
           read(33,*) x,(h(i,l),l=1+j,nval+j)
         enddo
      j=j+nval
      enddo
! rest
      if(.not.int(n3/nval)*nval==n3) then
          read(33,'(a)',end=201)a ! legend
         do i=1,n3
          read(33,*) x,(h(i,l),l=j+1,j+c)
         enddo
      endif
     endif
   enddo
201  continue

close(33)
end subroutine readORCAhess


subroutine readG09hess(n3,h,filen)
! read hessian from gaussian09 output.
use fiso, only: r8, stdout
implicit none
integer i,j,n3,c,l,ii,istart
real(r8) h(n3,n3),x
character(*) filen
character(120) a
logical da

inquire(file=filen,exist=da)
if(.not.da) write(stdout,*) 'cannot read ',filen,' (ERROR)'

write(stdout,*) 'Reading <',trim(filen),'>'
h=99
open(unit=33,file=filen)
c=mod(n3,5)
  do
     read(33,'(a)',end=201)a
     if(index(a,'Force constants in Cartesian coordinates:').ne.0) then
         read(33,'(a)',end=201)a ! dummy read
       istart=1
       do
       ii=0
          do i=istart,n3
              ii=ii+1
              l=min(5,ii)
              read(33,*) x,(h(i,j),j=istart,l+istart-1)
          enddo
         read(33,'(a)',end=201)a ! dummy read
          istart=istart+5
          if(istart>n3) exit
       enddo
      exit
     endif
   enddo
201  continue
close(33)

! fill full matrix
  do i=1,n3
     do j=1,i
     if(i==j) cycle
        h(j,i)=h(i,j)
     enddo
  enddo
end subroutine readG09hess


subroutine wtm(outfile)
use parm, only : xyz,nat,iat
implicit none
integer i
character(*) outfile
character(2) esym

open(unit=44,file=outfile)
write(44,'(a)')'$coord'
do i=1,nat
  write(44,'(3F18.12,2x,a2)') xyz(1:3,i) , esym(iat(i))
enddo
write(44,'(a)')'$end'
close(44)
end subroutine wtm


subroutine error(aa)
! error message +step
! io unit is optional (default=stdout)
use fiso, only: stdout
implicit none
integer io
! integer, optional :: setio
character(*), intent(in) :: aa

io=stdout
! if(present(setio)) io=setio

write(io,'(a)')
write(io,'(3x,a)')'************  ERROR  *********************'
write(io,'(4x,a)') trim(aa)
write(io,'(3x,a)')'******************************************'
write(io,'(a)')
stop
end subroutine error

subroutine warning(aa)
use fiso, only: stdout
implicit none
integer io
! integer, intent(in), optional :: setio
character(*), intent(in) :: aa

io=stdout
! if(present(setio)) io=setio

write(io,'(a)')
write(io,'(3x,a)')'!!--------------W A R N I N G------------------ !!'
write(io,'(4x,a)') trim(aa)
write(io,'(3x,a)')'!!--------------------------------------------- !!'
! stop
end subroutine warning


subroutine exclaim(aa)
use fiso, only: stdout
implicit none
integer io
character(*) aa
io=stdout
write(io,'(a)')
write(io,'(3x,a)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(io,'(4x,a)') trim(aa)
write(io,'(3x,a)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end subroutine





subroutine readamber(nat,iat,xyz,pc,ncrd,ntop)
use fiso, only: r8,stdout
use constant, only: au2ang
implicit none
integer i,nat,iat(nat)
real(r8) xyz(3,nat),mass(nat),pc(nat)
character(*) ncrd,ntop
character(80) F1,aa
real(r8) f,thres


! read prmtop file (ntop)
write(stdout,'(a)') ' reading ',trim(ntop)
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

pc=pc/18.2223_r8
! write(stdout,*) pc(1:5)

!if(iat(1)==0) stop 'something went wrong in readamber'

! use masses ....
if(iat(1)==0) then
write(stdout,'(a)') 'flag %FLAG ATOMIC_NUMBER not found. Guessing from masses'
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

thres=1e-4_r8
do i=1,nat
if( abs(mass(i)-1.008_r8)<thres ) iat(i)=1
if( abs(mass(i)-16.0_r8)<thres ) iat(i)=8
if( abs(mass(i)-12.01_r8)<thres) iat(i)=6
if( abs(mass(i)-12.001_r8)<thres) iat(i)=6
if( abs(mass(i)-14.01_r8)<thres) iat(i)=7
if( abs(mass(i)-30.97_r8)<thres) iat(i)=15
if(iat(i)==0) then
 write(stdout,*) 'Mass', mass(i)
 call error('undefined mass encoutered. stopping')
endif
enddo

endif



! read crd file (ncrd)

f=au2ang
open(222,file=ncrd)
read(222,'(a)') aa ! comment
read(222,'(a)') aa ! number of atoms
read(222,*) xyz(1:3,1:nat)
close(222)
xyz=xyz/f
end subroutine


subroutine printEIG(io,e,norbs)
use fiso, only: r8
implicit none
real(r8)  e(*)
integer n,ntimes,nrest,norbs,j,k,i,n2,io

n=6
ntimes=norbs/n
nrest=mod(norbs,n)
if(ntimes.eq.0) nrest=norbs
j=1
n2=n
do k=1,ntimes
   write(io,100)(i,i=j,n2)
   write(io,300)(e(i),i=j,n2)
   j =j +n
   n2=n2+n
enddo
if(nrest.gt.0.or.ntimes.eq.0) then
  write(io,100)(i,i=j,j+nrest-1)
  write(io,300)(e(i),i=j,j+nrest-1)
endif

 100  format(' #    : ',2x,6(4x,i4,2x))
 300  format('eigval: ',2x,6f10.3)
end subroutine


subroutine cpfile(infile,outfile,pr)
! simple file copy
use fiso, only: stdout
implicit none
character(*) infile, outfile,pr
integer i,o
character(255) line
if(pr=='print') write(stdout,'(a,a,a,a)') 'copying ',infile,' to ',outfile
open(newunit=i,file=infile)
open(newunit=o,file=outfile)
do
    read(i,'(a)',end=666) line
    write(o,'(a)') trim(line)
enddo
666 continue
close(i)
close(o)
end subroutine



!********************************************************************
!* reading gaussian input and re-writing with updated coordinates   *
!* we start processing from the first # occurance                   *
!********************************************************************
subroutine IOgaussian(gfile)
use parm, only : nat,xyz,iat,i
use constant, only: au2ang
implicit none
integer io,tmp,empty,k
character(*) gfile
character(120) aa
character(2) esym
!~ io=333
empty=0
k=0
open(newunit=io,file=gfile) ! fortran 2008 feature
open(newunit=tmp,file='xopt.tmp.IOgaussian')
do
k=k+1
 read(io,'(a)',end=666) aa
 write(tmp,'(a)') trim(aa)
   if(aa=='') empty=empty+1
   if(empty==2) then
    empty=3
    read(io,'(a)') aa ! charge,mult line
    write(tmp,'(a)') trim(aa)
    do i=1,nat
      write(tmp,'(a2,5x,3(F14.10,3x))') esym(iat(i)), xyz(1,i)*au2ang,xyz(2,i)*au2ang,xyz(3,i)*au2ang
      read(io,'(a)') aa ! to advance in the file
    enddo
   endif
enddo
666 continue
close(io)
close(tmp)
call cpfile('xopt.tmp.IOgaussian',gfile,'.')
open(newunit=tmp,file='xopt.tmp.IOgaussian',status='old')
close(tmp,status='delete')
end subroutine


subroutine newxyz(nat,iat,xyz)
use fiso, only: r8
use progs, only: xyzfile
use constant, only:au2ang
implicit none
integer k,nat,iat(nat)
real(r8) xyz(3,nat),f
character(2) esym

f=au2ang
open(unit=55,file=xyzfile)
write(55,'(I5)') nat
write(55,'(I2,2x,I2)') 0,0
  do k=1,nat
    write(55,'(a2,5x,3(F18.14,3x))') esym(iat(k)), xyz(1,k)*f,xyz(2,k)*f,xyz(3,k)*f
  enddo
close(55)
end subroutine newxyz

! read last n geometries from xopt.log
! e.g. for gdiis n=idiis
subroutine readlog_geo(nat,n,logxyz)
use fiso, only: r8
use constant, only: au2ang
implicit none
integer i,j,ic,n,num,ind,nat
character(120) aa
real(r8) qxyz(3,nat),logxyz(n,3,nat)

ic=0
! count geometries in xopt.log
open(unit=45,file='xopt.log',status='old')
do
  read(45,'(a)',end=665) aa
  call cstring(aa,num)
  if(num==1) ic=ic+1
enddo
665 close(45)
!print*,'# structures', ic

ind=ic-n
ic=0
!print*, 'reading', ind, ind+n
! read last n geometries
open(unit=45,file='xopt.log',status='old')
do
  read(45,'(a)',end=667) aa
  call cstring(aa,num)
   if(num==1) ic=ic+1
   if(ic==ind) then
    backspace(45)
     do i=1,n
       read(45,'(a)') aa ! #atoms
       read(45,'(a)') aa ! energy line
        do j=1,nat
         read(45,*) aa,qxyz(1,j),qxyz(2,j),qxyz(3,j)
!         print*, qxyz(1,j),qxyz(2,j),qxyz(3,j)
        enddo
!        print*,'-------------'
          logxyz(i,1,1:nat)=qxyz(1,1:nat)
          logxyz(i,2,1:nat)=qxyz(2,1:nat)
          logxyz(i,3,1:nat)=qxyz(3,1:nat)
      enddo
    goto 667 ! exit
   endif
enddo
667 continue ! exit point
close(45)

logxyz=logxyz/au2ang
end subroutine

subroutine readlog_grad(nat,n,loggrad)
use fiso, only: r8
implicit none
integer i,j,ic,n,num,ind,nat
character(120) aa,wx
real(r8) g(3,nat),loggrad(n,3,nat)

! get total grads in xopt.grad.log

open(unit=44,file='xopt.grad.log',status='old')
do
 read(44,'(a)',end=112) aa
 if(index(aa,'$iter').ne.0) then
   call charXsplit(aa,wx,2)
   read(wx,*) num
 endif
enddo
112 close(44)

ind=num-n
ic=0
open(unit=44,file='xopt.grad.log',status='old')
  read(44,'(a)') aa ! title line
do
  read(44,'(a)',end=666) aa
  if(index(aa,'$iter').ne.0) then
   call charXsplit(aa,wx,2)
   read(wx,*,end=111) num
!   print*,num,ind
   if(num==ind) then
     backspace(44)
     do i=1,n
       read(44,'(a)',end=111) aa ! info line
       do j=1,nat
         read(44,*) g(1,j),g(2,j),g(3,j)
!         print*,j, g(1,j),g(2,j),g(3,j)
       enddo
       loggrad(i,1,1:nat)=g(1,1:nat)
       loggrad(i,2,1:nat)=g(2,1:nat)
       loggrad(i,3,1:nat)=g(3,1:nat)
     enddo
     goto 666 ! exit
   endif
  endif
enddo
111 call error('reading error in xopt.grad.log')
666 continue ! exit point
close(44)
end subroutine

! make scrdir in specified path
subroutine para_scrdir2(basepath)
use progs, only: scrdir
use fiso, only: stdout
implicit none
integer pid,getpid
character(*) basepath
character(200) spid,pwd,aa

! process ID to string
pid=getpid()
write(spid,*) pid

! current path
call get_environment_variable('PWD',pwd)

write(scrdir,'(a)') '/'//trim(basepath)//'/xopt-scr-'//trim(adjustl(spid))
write(stdout,*) 'scratch name: ', trim(scrdir)
aa='mkdir '//trim(scrdir)
call system(aa)
end subroutine

!subroutine change_dir(dirname)
!use progs, only: workdir
!implicit none
!character(*) dirname
!character(200) pwd,path

!path=adjustl(trim(workdir))//trim(adjustl(dirname))
!call chdir(trim(path))
!end subroutine

!subroutine copy2scr(n,list)
! simple for now
subroutine copy2scr()
use progs, only: scrdir,workdir
implicit none
character(300) aa
aa='cp '//trim(workdir)//'/* '//trim(scrdir)//'/. 2> /dev/null '
call system(aa)
end subroutine

! scrdir back to workdir (make sure you actually want that)
subroutine scr2work()
use progs, only: scrdir, workdir
implicit none
character(300) aa
aa='cp '//trim(scrdir)//'/* '//trim(workdir)//'/. 2> /dev/null '
call system(aa)
end subroutine


! remove file and do selected actions
! switch=1 : be verbose
! switch=2 : stop program & print error
subroutine rmfile(string,switch)
  use fiso, only : stdout
implicit none
logical da
character(*) string
integer switch

inquire(file=trim(string),exist=da)
if(da) then
 open(unit=1234, file=trim(string), status='old')
 close(1234, status='delete')
if(switch==2) call error('user requested stop')
if(switch==1) write(stdout,*) 'removing ', trim(string)
endif


end subroutine


! set psi4 python commands
subroutine IOpsi4(filen)
implicit none
character(*) filen
character(200) aa
integer io,iot
logical fstr

open(newunit=io,file=filen,status='old',position='rewind')
 read(io,'(a)',end=666) aa
 if(fstr(aa,'xopt')) then ! already set, close and return
   close(io)
   return
 endif
 rewind(io)
open(newunit=iot,file='psi4.tmp',status='new')

write(iot,'(a)') '# set by xopt'
write(iot,'(a)') "qmol = qcdb.Molecule.init_with_xyz('xopt.xyz')"
write(iot,'(a)') 'mol = geometry(qmol.create_psi4_string_from_molecule())'
write(iot,'(a)') 'mol.update_geometry()'
write(iot,'(a)') "mol.reset_point_group('c1')"
! copy rest of the content
do
 read(io,'(a)',end=666) aa
 write(iot,'(a)') trim(aa)
enddo
666 continue
close(io)
close(iot)

call cpfile('psi4.tmp',filen,'n')
call rmfile('psi4.tmp',0)

end subroutine


subroutine IOmopac12(filen)
! write mopac.dat file with local SETUP file
use parm, only: xyz,iat,nat
use constant, only: au2ang
implicit none
character(*) filen
character(2) esym
integer i,io

open(newunit=io,file=filen)
write(io,'(a)') ' SETUP'
write(io,'(a)') ' written by xopt'
write(io,'(a)') ' '
do i=1,nat
 write(io,'(a2,5x,3(F18.14,3x))') esym(iat(i)), xyz(1,i)*au2ang,xyz(2,i)*au2ang,xyz(3,i)*au2ang
enddo
close(io)

end subroutine

subroutine IOgamess(filen)
! update DATA section of gamess input with current coordinates
use parm, only: xyz,iat,nat
use constant, only: au2ang
implicit none
character(*) filen
character(255) s
character(2) esym
integer i,wio,rio
logical fstr

open(newunit=rio,file=filen)
open(newunit=wio,file='xopt.gms.tmp')
do
 read(rio,'(a)',end=666) s
 write(wio,'(a)') trim(s)
 if(fstr(s,'$DATA')) then
  write(wio,'(a)') ' xopt xyz'
  read(rio,'(a)',end=666) s ! dummy
  write(wio,'(a)') ' C1 '
  read(rio,'(a)',end=666) s ! dummy
  do i=1,nat
    read(rio,'(a)',end=666) s ! dummy
    write(wio,'(2x,a2,1x,I3,3(F18.14,3x))') esym(iat(i)),iat(i),xyz(1,i)*au2ang,xyz(2,i)*au2ang,xyz(3,i)*au2ang
  enddo
 endif
enddo
666 continue
close(wio)
close(rio)
call cpfile('xopt.gms.tmp',trim(filen),'.')
end subroutine


subroutine IOnwchem(filen)
! exchange geometry block in nw.in
character(*), intent(in)  :: filen
character(255) :: s
integer        :: wio,rio
logical        :: fstr,wait


wait=.false.
open(newunit=rio,file=filen)
open(newunit=wio,file='xopt.nwchem.tmp')
do 
 read(rio,'(a)',end=666) s
 if(fstr(s,'GEOMETRY').or.fstr(s,'geometry')) then
   wait=.true.
   write(wio,'(a)') 'GEOMETRY units angstrom noautoz'
   write(wio,'(a)') 'LOAD xopt.xyz'
   write(wio,'(a)') 'END'
 endif
 if(.not.wait) write(wio,'(a)') trim(s)
 if(fstr(s,'END').or.fstr(s,'end')) wait=.false.
enddo
666 continue
close(wio)
close(rio)
call cpfile('xopt.nwchem.tmp',trim(filen),'.')
end subroutine


subroutine IOdftbplus(filen)
use parm, only: xyz,iat,nat
use constant, only: au2ang
implicit none
integer :: io
character(*) filen

open(newunit=io,file=filen)
write(io,'(''Geometry = GenFormat {'')')
write(io,'(i4,'' C'')') nat

end subroutine


subroutine write_genformat(io)
use parm, only: xyz,iat,nat
use constant, only: au2ang
implicit none
integer, intent(in) :: io
integer :: ic(107),i,c,k,comp(107),id(nat)

! unqiue elements in molecule
comp=0
ic=0
id=0
c=0
do i=1,nat
 do k=1,107
  if(iat(i)==k.and..not.any(id==k)) then
    c=c+1
    comp(i)=c
    id(i)=k
  endif
 enddo
enddo
write(io,*) id
stop

! write
xyz=xyz*au2ang
do i=1,nat
write(io,'(2i4,3(1x,F22.14))') i,id(i),xyz(1,i),xyz(2,i),xyz(3,i)
enddo

end subroutine

real(r8) function hubdata(iat)
use fiso, only: r8
implicit none
integer :: iat
hubdata=99
if(iat==1) hubdata=-0.1857
if(iat==2) hubdata=-0.1450
if(iat==6) hubdata=-0.1492
if(iat==7) hubdata=-0.1535
if(iat==8) hubdata=-0.1575
if(iat==9) hubdata=-0.1623
if(iat==11) hubdata=-0.0350
if(iat==12) hubdata=-0.0200
if(iat==15) hubdata=-0.1400
if(iat==16) hubdata=-0.1100
if(iat==17) hubdata=-0.0697
if(hubdata>90) call error('missing Hubbard parameters')
end function
