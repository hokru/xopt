subroutine loadrc()
use parm
use popt
use logic
use progs, only: &
  usrscr,scall_orca,scall_psi4,scall_gaus,scall_mopac, &
  scall_gei,scall_nwchem,scall_gms,scall_amber,scall_mpi,prog_flags
use MDdat
use strings
use fiso, only: stdout
implicit none
character(255) :: homedir,sfile,string
logical :: da,fstr
integer io
character(6), parameter :: fs='(2x,a)'

call get_environment_variable('HOME', homedir)
sfile=trim(homedir)//'/.xoptrc'
inquire(file=sfile,exist=da)
if(da) then
  write(*,fs) 'Found global config: '//trim(sfile)
endif

open(newunit=io,file=trim(sfile))
do
 read(io,'(a)',end=666) string
 if(fstr(string,'call_orca')) then
    i=index(string,'=')+1
    j=len(string)
    scall_orca=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom ORCA call:'
    write(stdout,fs)'  -> '//trim(scall_orca)
 endif
 if(fstr(string,'call_gaus')) then
    i=index(string,'=')+1
    j=len(string)
    scall_gaus=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom Gaussian call:'
    write(stdout,fs)'  -> '//trim(scall_gaus)
 endif
 if(fstr(string,'call_psi4')) then
    i=index(string,'=')+1
    j=len(string)
    scall_psi4=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom PSI4 call:'
    write(stdout,fs)'  -> '//trim(scall_psi4)
 endif
 if(fstr(string,'call_gei')) then
    i=index(string,'=')+1
    j=len(string)
    scall_gei=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom GEI call:'
    write(stdout,fs)'  -> '//trim(scall_gei)
 endif
 if(fstr(string,'call_amber')) then
    i=index(string,'=')+1
    j=len(string)
    scall_amber=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom AMBER call:'
    write(stdout,fs)'  -> '//trim(scall_amber)
 endif
 if(fstr(string,'call_gms')) then
    i=index(string,'=')+1
    j=len(string)
    scall_gms=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom GAMESS call:'
    write(stdout,fs)'  -> '//trim(scall_gms)
 endif
 if(fstr(string,'call_mpi')) then
    i=index(string,'=')+1
    j=len(string)
    scall_mpi=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom MPI(pgrag) call:'
    write(stdout,fs)'  -> '//trim(scall_mpi)
 endif
  if(fstr(string,'call_nwchem')) then
    i=index(string,'=')+1
    j=len(string)
    scall_nwchem=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom NWCHEM call:'
    write(stdout,fs)'  -> '//trim(scall_nwchem)
 endif
  if(fstr(string,'call_mopac')) then
    i=index(string,'=')+1
    j=len(string)
    scall_mopac=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'custom mopac call:'
    write(stdout,fs)'  -> '//trim(scall_mopac)
 endif
  if(fstr(string,'scratch')) then
    i=index(string,'=')+1
    j=len(string)
    usrscr=trim(adjustl(string(i:j )))
    write(stdout,fs,advance='no')'scratch '
    write(stdout,fs)'  : '//trim(usrscr)
 endif

!
  if(fstr(string,'maxiter')) call str_parse(string,2,maxiter)
  if(fstr(string,'maxd')) call str_parse(string,2,maxd)
  if(fstr(string,'gconv')) call str_parse(string,2,gconv)
  if(fstr(string,'econv'))  call str_parse(string,2,econv)
  if(fstr(string,'maxgrad')) call str_parse(string,2,econv)
  if(fstr(string,'gnorm')) call str_parse(string,2,gnorm)
enddo
666 continue
close(io)

end subroutine


subroutine eval_options (infile)
use logic
use popt
use progs, only: usrscr,prog_flags
use MDdat
use fiso, only: r8, stdout
implicit none
integer i,maxarg
character(80), allocatable :: arg(:)
character(255), intent(out):: infile
character(80) :: ftmp,narg
logical :: fstr, processed
integer :: s2i
real(r8) :: s2r

 call get_command_argument(1,infile)
 infile=trim(infile)
 if(infile=='') call error('no structure file!!')
 if(infile=='-h'.or.infile=='-help') call help
! if(infile=='restart') call set_restart()

maxarg=command_argument_count()
if(maxarg.gt.0) then

 allocate(arg(maxarg))
   do i=1,maxarg
     call get_command_argument(i,arg(i))
   enddo

 do i=2,maxarg
  processed=.false.
  ftmp=arg(i)
  if(i/=maxarg) narg=trim(arg(i+1))
  if(index(ftmp,'-h ').ne.0) then
   call help
  endif
!  if(index(ftmp,'-noprint ').ne.0) echo=.false.
  ! <<<<<         PROGS              >>>>>>
  if(index(ftmp,'-tm ').ne.0) tm=.true.
  if(index(ftmp,'-tmcc ').ne.0) then
   tm=.true.
   tmcc=.true.
  endif
  if(index(ftmp,'-riper ').ne.0) then
   tm=.true.
   riper=.true.
  endif
  if(index(ftmp,'-grad ').ne.0) then
    iopt=666 ! stop after printing the gradient
    write(stdout,'(a)') ' ** GRADIENT ONLY CALCULATION **'
    write(stdout,'(a)') ' energy/gradient in (hidden) .XOPT file '
  endif
  if(index(ftmp,'-readonly ').ne.0) readonly=.true.
  if(index(ftmp,'-tmhuge').ne.0) tmhuge=.true.
  if(index(ftmp,'-gms').ne.0) gamess=.true.
  if(index(ftmp,'-xtb ').ne.0) xtb=.true.
  if(index(ftmp,'-g09 ').ne.0) gaus=.true.
  if(index(ftmp,'-psi4 ').ne.0) psi4=.true.
  if(index(ftmp,'-mop ').ne.0) mopac=.true.
  if(index(ftmp,'-mopac ').ne.0) mopac=.true.
  if(index(ftmp,'-orca ').ne.0) orca=.true.
  if(index(ftmp,'-obabel ').ne.0) openbabel=.true.
  if(index(ftmp,'-relax ').ne.0) relax=.true.
  if(index(ftmp,'-driver ').ne.0) driver=.true.
  if(index(ftmp,'-gdiis ').ne.0) do_gdiis=.true.
  if(index(ftmp,'-gei ').ne.0) then
    gei=.true.
!    call_gei=trim(narg) ??
  endif
  if(index(ftmp,'-ppot ').ne.0) then
   ppot=.true.
   tm=.true.
  endif
  if(index(ftmp,'-apbs ').ne.0) APBS=.true.
!  if(index(ftmp,'-oniom ').ne.0) qmmm=.true.
  if(index(ftmp,'-amber ').ne.0) then
    amber=.true.
    maxd=0.1_r8
    maxiter=1000
  endif
  if(index(ftmp,'-rhess ').ne.0) then
   readHess=.true.
   hessname=trim(narg)
  endif
  if(index(ftmp,'-restart ').ne.0) then
   restart=.true.
  endif
  if(index(ftmp,'-extd3 ').ne.0) then
    d3=.true.
    d3func=trim(arg(i+1))
  endif
  if(index(ftmp,'-extgcp ').ne.0) then
    gcp=.true.
    gcplevel=trim(arg(i+1))
  endif
  ! library versions of gcp and dftd3
   if(index(ftmp,'-d3 ').ne.0) then 
    exlib=.true.
    d3=.true.
    d3func=trim(arg(i+1))
  endif
  if(index(ftmp,'-gcp ').ne.0) then
    exlib=.true.
    gcp=.true.
    gcplevel=trim(arg(i+1))
  endif
  if(index(ftmp,'-flags ').ne.0) then
    prog_flags=trim(arg(i+1))
  endif
  if(index(ftmp,'-numhess ').ne.0) hmodel='numerical'
  if(index(ftmp,'-numgrad ').ne.0) numgrad=.true.
  if(index(ftmp,'-hmodel ').ne.0) then
   
   if(s2i(narg)==1) then
     hmodel='hdiag' ! cart unit Hessian
   elseif(s2i(narg)==2) then
     hmodel='simple'! simple Hdiag (default)
   elseif(s2i(narg)==3) then
      hmodel='alm'   ! Fischer Almloef
   elseif(s2i(narg)==4) then
      hmodel='lindh' ! simplified Lindh
   elseif(s2i(narg)==5) then
      hmodel='mass' ! simplified Lindh
   else
      hmodel=trim(narg)
   endif
  endif

  ! <<<<<         Minimizer              >>>>>>
  if(index(ftmp,'-cart ').ne.0) iopt=11
  if(index(ftmp,'-cart-si ').ne.0) iopt=12
  if(index(ftmp,'-anc ').ne.0) iopt=21
  if(index(ftmp,'-anc-si ').ne.0) iopt=22
  if(index(ftmp,'-anc-cg ').ne.0) iopt=23
  if(index(ftmp,'-anc-cgrfo ').ne.0) iopt=24
  if(index(ftmp,'-anc-prfo ').ne.0) iopt=25
!  if(index(ftmp,'-anc-prfo ').ne.0) iopt=25
!  if(index(ftmp,'-anc-gdiis ').ne.0) iopt=26
  if(index(ftmp,'-devrfo ').ne.0) iopt=99
  if(index(ftmp,'-int ').ne.0) iopt=31
  if(index(ftmp,'-bfgs ').ne.0) iupdate=3
  if(index(ftmp,'-sr1 ').ne.0) iupdate=2
  if(index(ftmp,'-psb ').ne.0) iupdate=7
  if(index(ftmp,'-sr1-bfgs ').ne.0) iupdate=1
  if(index(ftmp,'-ss-bfgs ').ne.0) iupdate=4
  if(index(ftmp,'-ss-sr1-bfgs ').ne.0) iupdate=5
  if(index(ftmp,'-sr1-psb ').ne.0) iupdate=6
  



  ! <<<<<         parallel             >>>>>>
  if(index(ftmp,'-n ').ne.0) then
    nproc=s2i(narg)
  endif
  if(index(ftmp,'-omp  ').ne.0) then
    nomp=s2i(narg)
  endif

  ! <<<<<         scratch             >>>>>>
  if(index(ftmp,'-scratch ').ne.0) then
    usrscr=narg
    call para_scrdir2(trim(usrscr)) ! sets scrdir
    scratchjob=.true.
  endif


  ! <<<<<         Setup             >>>>>>
  if(fstr(ftmp,'-c ')) maxiter=s2i(narg) ! max iterations
  if(fstr(ftmp,'-iconv ')) iconv=s2i(narg) ! which thresholds are to be reached
  if(fstr(ftmp,'-dpl ')) maxd=s2r(narg) ! max displacement
  if(fstr(ftmp,'-gconv ')) gconv=s2r(narg) ! gnorm threshold
  if(fstr(ftmp,'-econv ')) econv=s2r(narg) ! energy threshold

 ! <<<<<         TS search             >>>>>>
  if(fstr(ftmp,'-ts ')) then
     tsopt=.true.
    !  iopt=25 <--- needs some re-work
    iopt=21 
    !  maxd=0.2_r8
     iupdate=6 
  endif
  if(fstr(ftmp,'-tsmode ')) then
    tsopt=.true.
    tsmode=s2i(narg)
  endif

  ! <<<<<         Compound keywords           >>>>>>
  if(fstr(ftmp,'-tight ')) tightopt=.true.
  if(fstr(ftmp,'-loose ')) preopt=.true.


  ! <<<<<         penalty function CI optimization           >>>>>>
  if(fstr(ftmp,'-ciopt ')) then
     ciopt=.true.
     econv=5e-7_r8
     gconv=1e-3_r8
     maxd=0.2_r8
     maxgrad=1e-3_r8 ! loose, unclear what is good
     dconv= 2e-2_r8 ! loose
  endif
  if(fstr(ftmp,'-sigma '))  sigma=s2r(narg) ! penalty CI opt
  if(fstr(ftmp,'-alpha '))  alpha=s2r(narg)! smoothing CI opt

  ! <<<<<         MD simulation           >>>>>>
  if(index(ftmp,'-md ').ne.0) then
    do_md=.true.
    iopt=900
  endif
  if(fstr(ftmp,'-dt '))  dt=s2r(narg)
  if(fstr(ftmp,'-steps ').or.fstr(ftmp,'-maxsteps ').or.fstr(ftmp,'-c ')) maxstep=s2i(narg)
  if(fstr(ftmp,'-temp ')) temp0=s2r(narg)
  if(fstr(ftmp,'-ifile ')) ifile=s2i(narg)
  if(fstr(ftmp,'-iprint ')) iprint=s2i(narg)
  if(fstr(ftmp,'-thermo ')) then
   thermo_nh=.false.
   if(s2i(narg)==1) thermo_nh=.true.
   if(s2i(narg)==2) thermo_berend=.true.
   if(s2i(narg)==3) thermo_scale=.true.
  endif
  
  ! <<<<<         IRC optimization           >>>>>>
  if(index(ftmp,'-dvv ').ne.0) then
    tsmode=1
    iopt=800
    do_dvv=.true.
  endif
  if(fstr(ftmp,'-dvv_err ')) dvv_err=s2r(narg)  ! error
  if(fstr(ftmp,'-dvv_step ')) dvv_step=s2r(narg)  ! step length/damping
  if(fstr(ftmp,'-dvv_init ')) dvv_init=s2r(narg)  ! direction (def: -1)
  if(fstr(ftmp,'-dvv_de ')) dvv_de=s2r(narg)  ! delta_E threshold to stop (def: 1e-8 au)
  if(fstr(ftmp,'-dvv_up ')) dvv_up=s2r(narg)  ! uphild threshold to stop (def: 5e-6 au)

  ! <<<<<         special keywords           >>>>>>
  ! move to COM and rotate to principle axis frame
  if(fstr(ftmp,'-orient')) orient=.true.

  if(fstr(ftmp,'-debug')) debug=.true.

  if(fstr(ftmp,'-exlib')) exlib=.true.
  if(fstr(ftmp,'-large')) large=.true.
  if(fstr(ftmp,'-d3hess')) d3hess=.true.

  if(index(ftmp,'-amber-pb ').ne.0) then
   write(stdout,'(a)') ' writing  amber.in for MMPB optimizations'
   call GetAmberinput('pb')
    amber=.true.
    maxd=0.01
    maxiter=1000
    iconv=2
    econv=1e-6
  endif
  if(index(ftmp,'-amber-gb ').ne.0) then
   write(stdout,'(a)') ' writing  amber.in for MMGB optimizations'
   call GetAmberinput('gb')
    amber=.true.
    maxd=0.01
    maxiter=1000
    iconv=2
    econv=1e-6
  endif
 enddo


endif
! need to safe all options into

if(tightopt) then
econv = 5e-8 ! SCF energy change
gconv = 1e-5 ! GNORM threshold
maxgrad = 1e-4     ! MAX grad component
dconv = 0.002 ! MAX displacment threshold
endif

if(preopt) then
econv = 5e-6 ! SCF energy change
gconv = 5e-3 ! GNORM threshold
maxgrad = 5e-4     ! MAX grad component
dconv = 0.02 ! MAX displacment threshold
endif


if(tsopt) then
if(tsmode<1) tsmode=1
!maxd=0.2
! iupdate=6 ! enforce sr1-psb update
endif

end subroutine



!* print dat HELP
subroutine help
use fiso, only: stdout
implicit none
integer io

! call get_environment_variable('XOPT_DIR',xoptdir)
! aa=trim(xoptdir)//'/'//"help.dat"
! ! read and print the help.dat file from the installation directory
! open(99,file=trim(aa))
! do
!   read(99,'(a)',end=123) a
!   write(6,'(a)') a
! enddo
! 123 close(99)
!

io=stdout
write(io,'(a)')'********************************'
write(io,'(a)')'**        H   E   L   P       **'
write(io,'(a)')'********************************'
write(io,'(a)')''
write(io,'(a)')'usage: xopt <coordinates> [-options]'
write(io,'(a)')''
write(io,'(a)')'coordinates can be in Turbomole format(bohrs)'
write(io,'(a)')'or in xyz coordinates (angstrom)'
write(io,'(a)')''
write(io,'(a)')''
 write(io,'(a)')'-h                 this help'
!write(io,'(a)')''-noprint
write(io,'(a)')''
write(io,'(a)')'PROGRAM INTERFACES'
 write(io,'(a)')'-tm                Turbomole(Xscf/Xgrad)'
 write(io,'(a)')'-tmcc              Turbomole(Xscf/ricc2)'
 write(io,'(a)')'-grad              stop after gradient(includes restraints!)'
 write(io,'(a)')'-readonly          (debug) stop after reading the setup'
 write(io,'(a)')'-tmhuge            switch to _huge Turbomole binaries '
 write(io,'(a)')'-g09               Gaussian 09'
 write(io,'(a)')'-mop/-mopac        mopac2012'
 write(io,'(a)')'-orca              ORCA'
 write(io,'(a)')'-gms               GAMESS'
 write(io,'(a)')'-psi4              PSI4'
 write(io,'(a)')'-nwchem            NWCHEM (under dev)'
 write(io,'(a)')'-xtb               xTB '
!  write(io,'(a)')'-obabel     ' did not work well
!  write(io,'(a)')'-apbs              ' no longer working/special version
!write(io,'(a)')''-oniom '        unfinished..
 write(io,'(a)')'-amber             sander module'
 write(io,'(a)')'-amber-pb          preset PB/MM optimization'
 write(io,'(a)')'-amber-gb          preset GB/MM optimization'
 write(io,'(a)')'-d3 <string>       dftd3 lib interface  '
 write(io,'(a)')'-gcp <string>      gcp   lib inferface  '
 write(io,'(a)')'-extd3 <string>    dftd3 system call  '
 write(io,'(a)')'-extgcp <string>   gcp   system call  '
 write(io,'(a)')''
 write(io,'(a)')'-numgrad           arbitrary numerical gradient (provide mygrad.sh)'
 write(io,'(a)')''
 write(io,'(a)')'HESSIANS'
 write(io,'(a)')'-rhess <filename>  read TM, ORCA or G09 hessian           '
 write(io,'(a)')'-restart           read xopt.restart.hess + xopt.restart if present'
 write(io,'(a)')'-numhess           calculate input hessian numerically'
 write(io,'(a)')'-hmodel <int>      select model Hessian'
 write(io,'(a)')'           1         cart. unit Hessian'
 write(io,'(a)')'           2         diagH(int) simple (default)'
 write(io,'(a)')'           3         diagH(int) Fischer/Almloef'
 write(io,'(a)')'           4         diagH(int) Lindh'
 write(io,'(a)')''
 write(io,'(a)')'PARALLEL SETTINGS'
 write(io,'(a)')'-n <int>           number of MPI threads, enables parallel mode'
 write(io,'(a)')'-scratch <string>  custom scratch directory '
 write(io,'(a)')''
 write(io,'(a)')'OPTIMIZER'
!write(io,'(a)')''-shess '
 write(io,'(a)')'-cart              cartesian coord. + RFO step'
 write(io,'(a)')'-cart-si           cartesian coord. + SI-RFO step'
 write(io,'(a)')''
 write(io,'(a)')'-anc               approx.normal coord. + RFO'
 write(io,'(a)')'-anc-si            approx.normal coord. + SI-RFO'
 write(io,'(a)')'-anc-cg            approx.normal coord. + CG'
 write(io,'(a)')'-anc-cgrfo         approx.normal coord. + RFO/CG mixture'
 write(io,'(a)')'-anc-prfo          approx.normal coord. + P-RFO (for TS only)'
 write(io,'(a)')''
 write(io,'(a)')'-ts                TS search (default: anc + EF + RFO)'
!  write(io,'(a)')'-ts                TS search (default: anc + EF + P-RFO)'
 write(io,'(a)')''
 write(io,'(a)')' CONICAL INTERSECTION SEARCH'
 write(io,'(a)')'-ciopt             penalty function conical interaction optimization'
 write(io,'(a)')'-sigma <float>     penalty for CIopt        '
 write(io,'(a)')'-alpha <float>     smoothing for CIopt        '
 write(io,'(a)')''
!write(io,'(a)')''-anc-prfo '
!write(io,'(a)')''-anc-gdiis '
!  write(io,'(a)')'-devrfo '
 write(io,'(a)')'-int               To be done...'
 write(io,'(a)')''
 write(io,'(a)')' HESSIAN UPDATE SCHEMES:'
 write(io,'(a)')'-sg1-bfgs          (default for minimization)'
 write(io,'(a)')'-sg1 '
 write(io,'(a)')'-bfgs '
 write(io,'(a)')'-ss-bfgs '
 write(io,'(a)')'-ss-sg1-bfgs       (experimental)'
 write(io,'(a)')'-sg1-psb           (default for TS search)'
 write(io,'(a)')'-psb              '
 write(io,'(a)')''

 write(io,'(a)')'MD (velocity verlet)'
 write(io,'(a)')'-md                run MD'
 write(io,'(a)')'-dt                time step [fs, default=0.5]'
 write(io,'(a)')'-temp              temperaure [K,default =300]'
 write(io,'(a)')'-steps             # steps [default =2000]'
 write(io,'(a)')'-thermo <int>      # thermostat: '
 write(io,'(a)')'                      Nose-Hoover chains   =  1'
 write(io,'(a)')'                      Berendsen rescaling  =  2'
 write(io,'(a)')'                      velocity scaling     =  3'


 write(io,'(a)')'IRC DRIVER'
 write(io,'(a)')'  (provide hessian, select eigenmode (default lowest) with -tsmode)'
 write(io,'(a)')'-dvv               IRC/DRC via DVV (damped velocity verlet)'
 write(io,'(a)')'-dvv_err  <flaot>  error (default 0.003 au)     '
 write(io,'(a)')'-dvv_step <flaot>  step length/damping (default 0.04 au/fs)       '
 write(io,'(a)')'-dvv_init <float>  directionality+damping (default= -1)'
 write(io,'(a)')'-dvv_de   <float>  dE threshold to stop (default= 1e-8 au)'
 write(io,'(a)')'-dvv_up   <float>  uphild threshold to stop (default= 5e-6 au)'
 

 write(io,'(a)')''
 write(io,'(a)')' OPTIONS'
 write(io,'(a)')'-c <int>           max. cycles'
 write(io,'(a)')'-dpl   <float>     max. displacement (default=0.35)'
 write(io,'(a)')'-gconv <float>     gradient norm threshold'
 write(io,'(a)')'-econv <float>     energy threshold'
 write(io,'(a)')'-tsmode <int>      vibrational mode to follow'
 write(io,'(a)')'-tight             econv=5e-8,gconv=1e-5'
 write(io,'(a)')'-loose             econv=5e-6,gconv=53e-3'
 write(io,'(a)')''
 write(io,'(a)')'  * penalty function CIopt*'
 write(io,'(a)')'-sigma <float>     penalty for CIopt'
 write(io,'(a)')'-alpha <float>     smoothing for CIopt'
 write(io,'(a)')''
 write(io,'(a)')'-orient            principle axis alignment'
 write(io,'(a)')''
 ! write(io,'(a)')'-ppot              ppot'
 write(io,'(a)')'-iconv <int>       <int>=1-4 to select how many thresholds to obey'
 write(io,'(a)')'-debug              switches on various debug printouts'
 write(io,'(a)')''
 write(io,'(a)')'-flags <str>        add additonal flags for: xtb'
 write(io,'(a)')''
 write(io,'(a)')'SPECIAL'
 write(io,'(a)')'-driver            driver'
 write(io,'(a)')''

write(io,'(a)')''
write(io,'(a)')" Good luck & Don't Panic"
write(io,'(a)')''

stop
end subroutine help


