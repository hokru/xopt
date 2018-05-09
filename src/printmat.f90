

!program test
!implicit none
!integer m,n,i,j
!real(8) r
!real(8), allocatable :: mat(:,:)
!real(8), allocatable :: vec(:)
!integer, allocatable :: imat(:,:)
!
!n=10
!m=12
!allocate(mat(n,m),imat(n,m),vec(n))
!do i=1,n
! do j=1,m
!  call random_number(r)
!  mat(i,j)=r
!  mat(i,j)=0
!  vec(i)=r*10000
!  if(i==j) mat(i,j)=1d0
!  if(i==j) imat(i,j)=1
! enddo
!enddo
!
!call printMat(6,n,m,mat,'test')
!call printMat(6,n,1,vec,'test vec')
! call printiMat(6,n,m,imat,'itest')
! end



subroutine printMat(io,dim1,dim2,mat,string)
! print matrix mat(n,m) to output unit io
! with descriptor string
! for vectors use dim2=1
use fiso, only: r8
implicit none
character(*) string
integer dim1,dim2,io,i,j,k
integer maxcol,nparts,left
integer col,row,icol
real(r8) mat(dim1,dim2)
character(80) widefmt,prfmt

! max width of printout
maxcol=10
! dynamic print format
write(widefmt,'(a,I2,a)') "(1x,5x,",maxcol,"(3x,I5,3x))"
write(prfmt,'(a,I2,a)')   "(1x,I5,",maxcol,"(F10.4,1x))"
!write(prfmt,'(a,I2,a)')   "(1x,I5,",maxcol,"(E10.3,1x))"
! header
write(io,'(1x,5x,a)') trim(string)

! initial indices
left=0
col=min(dim2,maxcol)
row=dim1
if(dim2<=maxcol) then
  nparts=1
  col=dim2
else
 nparts=dim2/maxcol
 left=mod(dim2,maxcol)
 col=maxcol
endif

if(dim2==1) then
 left=0
 maxcol=0
 col=1
endif

!main part
icol=1
do k=1,nparts
 write(io,widefmt) (i,i=icol,col)
 write(io,*)''
 do i=1,row
  write(io,prfmt) i, (mat(i,j),j=icol,col)
 enddo
 write(io,*)''
 col=col+maxcol
 icol=icol+maxcol
enddo

if(left>0) then
!leftover
col=col-maxcol+left
write(io,widefmt) (i,i=icol,col)
 write(io,*)''
do i=1,row
 write(io,prfmt) i, (mat(i,j),j=icol,col)
enddo
 write(io,*)''
endif
end subroutine


subroutine printiMat(io,dim1,dim2,mat,string)
! print matrix mat(n,m) to output unit io
! with descriptor string
! prints integer matrix
implicit none
character(*) string
integer dim1,dim2,io,i,j,k
integer maxcol,nparts,left
integer col,row,icol
integer mat(dim1,dim2)
character(80) widefmt,prfmt

! max width of printout
maxcol=10
left=0
! dynamic print format
write(widefmt,'(a,I2,a)') "(1x,4x,",maxcol,"(1x,I4,x))"
write(prfmt,'(a,I2,a)')   "(1x,I4,",maxcol,"(1x,I4,x))"
! print*,prfmt,widefmt

! header
write(io,'(1x,5x,a)') trim(string)

! initial indices
col=min(dim2,maxcol)
row=dim1

if(dim2<=maxcol) then
  nparts=1
  col=dim2
else
 nparts=dim2/maxcol
 left=mod(dim2,maxcol)
 col=maxcol
endif


!main part
icol=1
do k=1,nparts
 write(io,widefmt) (i,i=icol,col)
 write(io,*)''
 do i=1,row
  write(io,prfmt) i, (mat(i,j),j=icol,col)
 enddo
 write(io,*)''
 col=col+maxcol
 icol=icol+maxcol
enddo

if(left>0) then
!leftover
col=col-maxcol+left
write(io,widefmt) (i,i=icol,col)
 write(io,*)''
do i=1,row
 write(io,prfmt) i, (mat(i,j),j=icol,col)
enddo
 write(io,*)''
endif
end subroutine
