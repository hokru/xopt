!*************************
!* An eXternal OPTimizer *
!*************************

program  xopt
use fiso, only: stdout,i8
use parm
use logic
use popt
use progs
use timing
use internals
implicit none
integer istat
character(255) infile
logical da, newhess
real(r8) min_mem

da=.false.
newhess=.false.

! PGI has fully buffered I/O to files, but we want to see the progress step by step!
#ifdef PGI
call setvbuf(6,1,0)
#endif

! clock
call system_clock(count_rate=cr) ! count rate, once is enough
call progtime_start() ! end timing


! load defaults
call loaddef()
! load user config defaults
call loadrc()
! control options
call eval_options(infile)

! output file handling
if(do_output) then
  stdout=1234
  open(stdout,file=trim(output_name))
endif

! header
call head()

call rcontrol()

if(.not.do_md) newhess=.true.
if(do_dvv) then
 newhess=.false.
!  readHess=.true.
endif

! dimensions and memory
call tmolrd(infile,.false.,.true.) ! just determined nat
if(nat.gt.maxat) write(stdout,*) 'too many atoms, increase: maxat =',maxat
if(nat==0) call error(' no atoms found!')

npair=(nat*(nat-1))/2
nat3=nat*3

!write(stdout,'(2x,a)') ' min. memory requirements: '
min_mem=dble(nat)*9.0+(dble(nat)*3.0)**2
call p_memory(min_mem,'min. memory requirements:')
write(stdout,'(a)') ' '

allocate(xyz(3,nat),stat=istat)
allocate(xyz0(3,nat),stat=istat)
allocate(grad(3,nat),stat=istat)
allocate(iat(nat),stat=istat)

if(freeze) allocate(ifrez(nat),stat=istat)

if(restrain) then
 allocate(irest_vol(10,3),stat=istat)
 allocate(irest_atom(maxrestr),stat=istat)
 allocate(irest_bond(maxrestr,2),stat=istat)
 allocate(irest_ang(maxrestr,3),stat=istat)
 allocate(irest_dihed(maxrestr,4),stat=istat)
 allocate(irest_konst(maxrestr,5),stat=istat)
 allocate(val0(4*maxrestr),stat=istat)
 allocate(eneR(4*maxrestr),stat=istat)
endif


allocate(chess(3*nat,3*nat),stat=istat)
if(istat.ne.0) call error('allocation error: Hessian')

! read coordinates
call tmolrd(infile,.true.,.false.)
if(debug) write(stdout,*) 'tmolrd done'
xyz0=xyz

! print info
call p_info

! calculate restraining energy from restart file
if(ReadOnly) then
  restart=.true.
  call rcontrolrestrain
  call readbin(nat3,chess,'xopt.hess.restart')
  call addRestrainGrad
  call printFinalRestrain
  write(stdout,*) 'Reading done..stopping'
  stop
endif

! read control file for restraints/constraints
! ifrag=1
call rcontrolfreeze
call rcontrolrestrain

if(debug) write(stdout,*) 're/constraints done'


! write xopt.xyz
call wrxyz(trim(xyzfile))

! prepare programs
call prepPROG

write(stdout,'(a)') ''

! initial Hessians
if(restart) then
    newhess=.false.
! This is broken for some reason?!
    inquire(file='xopt.hess.restart',exist=da)
    if(freeze) then
     call warning('careful with restarts and $freeze. potentially incomplete Hessian.')
    endif
    if(da) then
        call readbin(nat3,chess,'xopt.hess.restart')
    else
        call exclaim('old hessian not found, making new.')
        newhess=.true.
    endif
endif
if(readHess) then
    newhess=.false.
    call checkhess(hessname,istat)
    if(istat==1) call rdhess(nat3,chess,hessname)       ! TM
    if(istat==2) call readORCAhess(nat3,chess,hessname) ! ORCA
    if(istat==3) call readG09hess(nat3,chess,hessname)  ! G09
    if(istat==4) call readPSI4hess(nat3,chess,hessname) ! PSI4
endif
if(newhess) then
    call status1(timer)
    call message_head('* initial Hessian *')
    call getHess
    call status2(timer)
endif


if(d3hess) then
    call status1(timer)
    call message_head('* numerical D3 Hessian *')
    call HcartD3(nat,iat,xyz,chess)
    call status2(timer)
endif

if(debug.and.nat<=50) call printmat(stdout,nat3,nat3,chess,'initial H(cart)')

! select minimzer or MD
select case(iopt) ! cartesians
 case(11:19)
    if(newhess) deallocate(hdiag,B)
    call copt
 case(21:29,99) ! normal coords
    if(newhess) deallocate(hdiag,B)
    call getanc
    call ancopt()
case(31:39)  ! internals
  nvar=0
  call intopt()
 case(666) ! gradient only
  call getgrad
  write(stdout,'(a)')' writing .XOPT '
  open(unit=331,file='.XOPT')
  write(331,*) energy
  do i=1,nat
    write(331,'(3E22.13)')grad(1,i),grad(2,i),grad(3,i)
  enddo

 ! MD section
 case(900)
   write(stdout,*) 'EXPERIMENTAL MD SECTION'
  call runMD
 ! IRC
 case(800)
 call irc_driver()
end select

! final geometries
call wrxyz('xopt.xyz')
call message('writing final structure: xopt.opt')
call wrxyz('xopt.opt')

if(restrain) then
 open(unit=34,file='xopt.restrain.tmp',status='old')
 close(34,status='delete')
endif

call progtime_end() ! end timing
call prtime(stdout) ! print timing std output


 ! delete scratch dirs for parallel gradients
if(nproc>1.and.scratchjob) call system('rm -r '//trim(scrdir))


! the end
call message('normal program termination')
call message(' ')
call message(' ')

end program xopt
