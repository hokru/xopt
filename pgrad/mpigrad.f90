module parallel
use mpi
#define master rank==0
integer ierr,nproc,rank,strlen
character(MPI_MAX_PROCESSOR_NAME) procname
character(300) scrdir,usrscr,syscall
end module


program mpigrad
use mpi
use parallel
implicit none
integer nat,i,j,k,l,ntasks,ptask,taskid,rank1,leftover
real(8), allocatable :: xyz(:,:)
real(8), allocatable :: grad(:,:)
real(8), allocatable :: lgrad(:)
real(8), allocatable :: lgrad2(:)
real(8), allocatable :: grad2(:,:)
integer, allocatable :: iat(:)
integer, allocatable :: istart(:),iend(:)
logical              :: da, debug,locals
character(255)       :: aa,pwd,pwdold
integer              :: id,rtasks,namel
integer, allocatable :: rcounts(:),offset(:)

real(8) e,g,f,g2

usrscr="/scratch/"
debug=.true.
locals=.true.
f=0.5291770d0
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

rank1=rank+1

allocate(istart(nproc),iend(nproc))
allocate(rcounts(nproc),offset(nproc))
rcounts=0
offset=0

if(master) then
  write(*,'(a)')''
  write(*,'(a)')' xopt: parallel numerical gradient helper '
  write(*,'(a)')''
  write(*,'(a)')''
  write(*,'(3x,a,i3,a)') '# Running  with  :', nproc ,' threads'

  ! read nat,xyz,iat from file
  open(99,file='xopt.para.tmp')
  read(99,*) nat
endif

call MPI_BCAST(nat,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
call MPI_GET_PROCESSOR_NAME(procname,namel, ierr)


allocate(xyz(3,nat),iat(nat),grad(3,nat),grad2(3,nat))
ntasks=3*nat
ptask=int(ntasks/nproc)
leftover=mod(ntasks,nproc)

if(master) then
  do i=1,nat
     read(99,*) xyz(1:3,i),iat(i)
  enddo
  read(99,'(a)',end=123) usrscr
  read(99,'(a)',end=123) syscall
123 continue
close(99)
xyz=xyz/f
if(trim(usrscr)=='local') locals=.true.
!if(.not.trim(usrscr)=="0") print*,' custom scratch location:', trim(usrscr)
if(.not.locals) print*,' custom scratch location:', trim(usrscr)
print*,' external script: ', trim(syscall)
endif

! broadcast data
call MPI_BCAST (usrscr,len(usrscr),MPI_CHARACTER,0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (syscall,len(syscall),MPI_CHARACTER,0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (xyz,3*nat,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
call MPI_BCAST (iat,nat,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)

if(master) then
  write(*,'(3x,a,i5)') '# tasks    : ', ntasks
  write(*,'(3x,a,i5,i5)') '# task/proc + remainder    : ', ptask, mod(ntasks,nproc)
  if(ntasks<nproc) stop 'tasks<nproc! What are you doing?'
  if(nproc==1) stop ' Just 1 proc! What are you doing? '
endif

! use rank1 to offset rank index starting from 0
istart(rank1)=1+ptask*rank
iend(rank1)=ptask*(rank+1)
! last proc also does the leftovers (if there are any)
if(rank1==nproc) iend(rank1)=ntasks
rtasks=1+(iend(rank1)-istart(rank1)) 
rcounts(rank1)=rtasks


if(master) then
 offset=0
 rcounts(1)=ptask
 do i=2,nproc
   offset(i)=offset(i-1)+ptask
   if(i==nproc) ptask=ptask+leftover
   rcounts(i)=ptask
 enddo
 if(debug) then
  do i=1,nproc
     print*,'proc/offsets/tasks',i,offset(i),rcounts(i)
  enddo
 endif
endif

allocate(lgrad(rtasks),lgrad2(rtasks))



print*,'I am ', trim(procname)
print*,'TASKS INFO (proc/start/end/#tasks):',rank1,istart(rank1),iend(rank1),rtasks
!print*,'Making scratch and copying mygrad.sh to scratch dir'
print*,'Making scratch and copying all data to scratch dir'

! MAKE SCRATCH, local or custom location
call getcwd(pwdold)
!if (.not.trim(usrscr)=="0") then
if (.not.locals) then
 call para_scrdir2(usrscr)
else
 call para_scrdir('xopt')
endif
call chdir(trim(scrdir))
call getcwd(pwd)
print*,'Running in ',trim(pwd) 


! COPY
aa='cp '//trim(pwdold)//'/* '//trim(scrdir)//'/. '
!aa='cp '//trim(pwdold)//'/mygrad.sh '//trim(scrdir)//'/.'
print*, trim(aa)
call system(aa)


print*,'Starting calculations...'
grad=0d0
lgrad=0d0
k=0
do taskid=istart(rank1),iend(rank1)
k=k+1
call get_displ_grad(nat,xyz,iat,taskid,g,g2)
lgrad(k)=g
lgrad2(k)=g2
enddo
! now each threads holds a part of the gradient



!print*,'gradient',lgrad
!call MPI_GATHER(lgrad,rtasks,MPI_DOUBLE_PRECISION,fgrad,rtasks,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  ! truncated message if there are leftovers! otherwise working
call MPI_Gatherv(lgrad,rcounts(rank1),MPI_DOUBLE_PRECISION,grad,rcounts,offset,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr) 
call MPI_Gatherv(lgrad2,rcounts(rank1),MPI_DOUBLE_PRECISION,grad2,rcounts,offset,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr) 

call chdir(pwdold)
if(master) then
open(99,file='xopt.pgrad.tmp')
! write(99,'(F18.10)') e
 print*,'gradient 1'
 do i=1,nat
  write(*,'(3F18.10,2x)') grad(1:3,i)
  write(99,'(3F18.10,2x)') grad(1:3,i)
 enddo
 print*,'gradient 2'
  do i=1,nat
  write(*,'(3F18.10,2x)') grad2(1:3,i)
  write(99,'(3F18.10,2x)') grad2(1:3,i)
 enddo

close(99)
endif

if(master) then
 write(*,'(a)') ''
 write(*,'(a)') 'all done'
 write(*,'(a)') ''
endif

! delete tmp dirs
aa='rm -r '//trim(scrdir)
call system(aa)

CALL MPI_FINALIZE(ierr)

end program

!*****************************************************************



! make scrdir in current path for process <id>  and return full path
! string = custom directory name
subroutine para_scrdir(basepath)
use parallel, only: rank,scrdir
integer pid,id
character(*) basepath
character(200) spid,pwd,aa
character(80) sid

! process ID to string
pid=getpid()
write(spid,*) pid

! current path
call getcwd(pwd)

write(sid,*) rank
write(scrdir,'(a)') trim(pwd)//'/'//trim(basepath)//'-tmp-'//trim(adjustl(spid))
print*, 'scratch name: ', trim(scrdir)
aa='mkdir '//trim(scrdir)
call system(aa)
end subroutine


! make scrdir in specified path 
subroutine para_scrdir2(basepath)
use parallel, only: rank,scrdir
integer pid,id
character(*) basepath
character(200) spid,pwd,aa
character(80) sid

! process ID to string
pid=getpid()
write(spid,*) pid

! current path
call getcwd(pwd)

write(sid,*) rank
write(scrdir,'(a)') '/'//trim(basepath)//'/xopt-scr-'//trim(adjustl(spid))
!call chdir(trim(basepath))
print*, 'scratch name: ', trim(scrdir)
aa='mkdir '//trim(scrdir)
call system(aa)
end subroutine




! gradient from single displacement
subroutine get_displ_grad(nat,xyz,iat,taskid,lgrad,lgrad2)
use parallel, only: rank,syscall
implicit none
integer i,j,k,l,nat,iat(nat),jmod,taskid
real(8) xyz(3,nat)
real(8) el,er,step,e,g,e2,el2,er2,lgrad2
real(8) lgrad,f
!real(8) grad(3,nat),lgrad,f
! COPY 
character(1) ind(3)
data ind/'x','y','z'/

! CALC!
!step=0.005d0
! coords are in Angstrom
step=0.005d0
!step=0.010

! get atom (i) and x,y,z (j) index from tasklist which runs from 1 to 3*nat as follows:
! atom1: 1 2 3
! atom2: 4 5 6
! atom3: 7 8 9
! ...
! probably not very smart:
if(taskid<=3) then
 i=1
elseif(mod(taskid,3)==0) then
 i=int(taskid/3d0)
else
 i=int(taskid/3d0)+1
endif

jmod=mod(taskid,3)
if(jmod==0) then
 j=3
else
 j=jmod
endif

print*,'proc',rank+1,'atom:',i,'component ',ind(j),' task ',taskid
xyz(j,i)=xyz(j,i)+step
call newxyz(nat,iat,xyz)
!call system("./mygrad.sh")
call system(trim(syscall))
call ene(er,er2)
xyz(j,i)=xyz(j,i)-2d0*step
call newxyz(nat,iat,xyz)
!call system("./mygrad.sh")
call system(trim(syscall))
call ene(el,el2)
xyz(j,i)=xyz(j,i)+step
!grad(j,i)=(er-el)/(step*2d0)
lgrad=(er-el)/(step*2d0)
lgrad2=(er2-el2)/(step*2d0)

end subroutine

subroutine ene(e,e2)
implicit none
real*8 e,e2
e=9d9
open(44,file='xopt.energy.tmp')
read(44,*,end=666) e
read(44,*,end=666) e2
666 continue
if(e==9d9) stop 'error. No energy!!'
if(e2==9d9) e2=0
close(44)
end subroutine


subroutine newxyz(nat,iat,xyz)
implicit none
integer k,nat,iat(nat)
real(8) xyz(3,nat),f
character(2) esym
character(80) xyzfile

xyzfile='xopt.xyz'

f=0.5291770d0
open(unit=55,file=xyzfile)
write(55,'(I5)') nat
write(55,'(I2,2x,I2)') 0,0
  do k=1,nat
    write(55,'(a2,5x,3(F18.14,3x))') esym(iat(k)), xyz(1,k)*f,xyz(2,k)*f,xyz(3,k)*f
  enddo
close(55)
end subroutine newxyz

character(2) FUNCTION ESYM(I)
implicit none
integer i
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

