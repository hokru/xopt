subroutine prepPROG
use logic
use progs
use fiso, only: stdout
implicit none
character(120) aa

write(stdout,'(a)') ''
if(ciopt) then
write(stdout,*) " *  penalty function-based conical intersection (CI) gradient       *"
write(stdout,*) " *  Levine et al.10.1021/jp0761618                                  *"
write(stdout,*) " *                  +                                               *"
    if(tm) then
         write(stdout,*) " *  Turbomole                                                  *"
    elseif(orca) then
         write(stdout,*) " *  ORCA                                                       *"
    elseif(gaus) then
         write(stdout,*) " *  Gaussian                                                   *"
    elseif(amber) then
         write(stdout,*) " *  sander                                                     *"
    elseif(numgrad) then
           write(stdout,*) " * using parallel numerical gradients                             *"
         if(nproc>1) then
           write(stdout,*) "   --> # proc: ", nproc
           call set_binpath(binpath)
         else
           write(stdout,*) " * using numerical gradients                                      *"
         endif
    else
        call error("selected internal interface not implemented")
    endif
!return
endif

!************************
!* QM/MM                *
!************************
if(qmmm) then
write(stdout,*) " * THIS IS A QM/MM OPTIMIZATION (good luck) *"
write(stdout,'(a)') ''
! anything to prepare?
return
endif


!****************
!* gamess         *
!****************
if(gamess) then
write(stdout,*) " * INTERNAL GAMESS INTERFACE *"
command_gms=trim(scall_gms)//redir//xjob
write(stdout,'(4x,a)') 'GAMESS system call:'
write(stdout,'(4x,a)') trim(command_gms)
write(stdout,'(4x,a)') ''

return
endif

!************************
!* GEI               *
!************************
if(GEI) then
write(stdout,*) " * GENERAL EXTERNAL INTERFACE *"
command_gei='./'//trim(scall_gei)
write(stdout,'(4x,a)') 'external system call:'
write(stdout,'(4x,a)') trim(command_gei)
write(stdout,'(4x,a)') ''
return
endif


!************************
!* xTB                  *
!************************
if(xtb) then
write(stdout,*) " * INTERNAL XTB INTERFACE *"
call rmfile('gradient',1)
prog_flags=' -grad '//trim(prog_flags)
command_xtb=trim(scall_xtb)//' '//trim(xyzfile)//trim(prog_flags)//redir//xjob
write(stdout,'(4x,a)') 'XTB system call:'
write(stdout,'(4x,a)') trim(command_xtb)
write(stdout,'(4x,a)') ''
return
endif


!****************
!* ORCA         *
!****************
if(orca) then
write(stdout,*) " * INTERNAL ORCA INTERFACE *"
call rmfile('orca.gbw',1)
command_orca=trim(scall_orca)//' '//trim(orcain)//redir//xjob
write(stdout,'(4x,a)') 'ORCA system call:'
write(stdout,'(4x,a)') trim(command_orca)
write(stdout,'(4x,a)') ''

return
endif

if(openbabel) then
write(stdout,*) " * INTERNAL OPENBABEL(GAFF) INTERFACE *"
write(stdout,'(a)') ''
! openbabel currently(?) cannot print "stand-alone" gradients
return
endif


! for external interface
if(numgrad) then
    write(stdout,*) " * GENERAL EXTERNAL INTERFACE *"
    command_gei='./'//trim(scall_gei)
    write(stdout,'(4x,a)') 'external system call:'
    write(stdout,'(4x,a)') trim(command_gei)
     write(stdout,'(4x,a)') ''
 if(nproc>1) then
   write(stdout,*) " * PARALLEL NUMERICAL GRADIENT REQUESTED*"
   write(stdout,*) "   --> # proc: ", nproc
   call set_binpath(binpath)
   write(aa,*) nproc
   command_mpi=trim(scall_mpi)//' -output-filename xopt.slave -n '//trim(adjustl(aa))//' '//trim(binpath)
   write(stdout,'(4x,a)') 'MPI system call:'
   write(stdout,'(4x,a)') trim(command_mpi)
   write(stdout,'(4x,a)') ''
 else
    write(stdout,*) " *  NUMERICAL GRADIENT REQEUSTED  *"
 endif
write(stdout,'(a)') ''
return
endif


!
! if(numgrad_internal) then
! write(stdout,*) " *  NUMERICAL GRADIENT REQEUSTED  *"
! write(stdout,'(a)') ''
! return
! endif


!************************
!* DRIVER (Grimme)      *
!************************
if(driver) then
write(stdout,*) " * INTERNAL DRIVER INTERFACE *"
write(stdout,'(a)') ''
write(stdout,*)' *** '
write(stdout,*)' driver input: '
open(55,file='.DRIVERINPUT')
do
 read(55,'(a)',end=551) aa
 write(stdout,*)trim(aa)
enddo
551 close(55)
write(stdout,*)' *** '
return
endif


!****************
!* MOPAC16      *
!****************
if(mopac) then
write(stdout,*) " * INTERNFAL MOPAC INTERFACE *"
write(stdout,'(4x,a)') 'MOPAC system call:'
command_mopac=trim(scall_mopac)//' mopac.dat  2> /dev/null '
write(stdout,'(4x,a)') trim(command_mopac)
write(stdout,'(4x,a)') ''
return
endif

!****************
!* TM           *
!****************
if(tm.or.tmhuge) then
write(stdout,*) " * INTERNAL TURBOMOLE INTERFACE *"
call rmfile('gradient',1)
call rmfile('energy',1)
call get_environment_variable('TURBODIR',aa)
write(stdout,'(4x,a)') 'TM install:'
write(stdout,'(4x,a)') trim(aa)
write(stdout,'(4x,a)') ''

 call system('actual -r >/dev/null')
 tmri=.false.
 write(stdout,*)' Checking control file for $rij'
 open(111,file='control')
  do
  read(111,'(a)',end=666) aa
    if(index(aa,'$rij').ne.0) then
      tmri=.true.
      write(6,'(5x,a)') ' -> found $rij option'
      exit
    endif
  enddo
 666 continue
 close(111)
 if(tmcc) write(6,'(5x,a)') ' -> ricc2 gradient'
 if(riper) write(6,'(5x,a)') ' -> riper module'
return
endif



!****************
!* gaussian09   *
!****************
if(gaus) then
  write(stdout,*) " * INTERNAL GAUSSIAN09 INTERFACE  *"
  command_gaus=trim(scall_gaus)//' g.in '//trim(xjob)
  write(stdout,'(4x,a)') 'Gaussian system call:'
  write(stdout,'(4x,a)') trim(command_gaus)
  write(stdout,'(4x,a)') ''
  return
endif

!****************
!* AMBER/SANDER *
!****************
if(amber) then
  !call system('cp amber.crd amber.rst')
  write(stdout,*) " * INTERNAL AMBER INTERFACE  *"
  command_amber=trim(scall_amber)
  write(stdout,'(4x,a)') 'Amber system call:'
  write(stdout,'(4x,a)') trim(command_amber)
  write(stdout,'(4x,a)') ''
  return
endif

!****************
!*      PSI4    *
!****************
if(psi4) then
  call IOpsi4('psi4.in')
  write(stdout,*) " * INTERNAL PSI4 INTERFACE  *"
  prog_flags=' '//trim(prog_flags)
  command_psi4=trim(scall_psi4)//trim(prog_flags)//' psi4.in '//trim(xjob)
  write(stdout,'(4x,a)') 'PSI4 system call:'
  write(stdout,'(4x,a)') trim(command_psi4)
  write(stdout,'(4x,a)') ''
  return
endif


!****************
!*  NWCHEM      *
!****************
if(nwchem) then
  call IOnwchem('nw.in')
  write(stdout,*) " * INTERNAL NWCHEM INTERFACE  *"
  command_nwchem=trim(scall_nwchem)
  write(stdout,'(4x,a)') 'NWCHEM system call:'
  write(stdout,'(4x,a)') trim(command_nwchem)
  write(stdout,'(4x,a)') ''
  return
endif

!****************
!*   GAMESS     *
!****************
if(gamess) then
  call IOgamess(trim(gmsin)//'.inp')
  write(stdout,*) " * INTERNAL GAMESS INTERFACE  *"
  command_gms=trim(scall_gms)
  write(stdout,'(4x,a)') 'GAMESS system call:'
  write(stdout,'(4x,a)') trim(command_gms)
  write(stdout,'(4x,a)') ''
  return
endif


call error('select prog for energy and gradient!')

end subroutine prepPROG


subroutine set_binpath(binpath)
use fiso, only:stdout
implicit none
character(200) binpath
 call get_environment_variable('HOME', binpath)
 binpath=trim(binpath)//'/bin/xopt.pgrad'
 write(stdout,*) "   using: ", trim(binpath)
end


subroutine checkhess(fname,hlogic)
use fiso, only:stdout
implicit none
character(*) fname
integer hlogic
character(120) aa
integer nstrings

! hlogic:
! 1 = TM (default)
! 2 = ORCA
! 3 = G09
! 4 = PSI4

hlogic=-1
open(11,file=fname)
do
  read(11,'(a)',end=999) aa
  call cstring(aa,nstrings)
  if(index(aa,'$hessian').ne.0) hlogic=1
  if(index(aa,'orca_hessian_file').ne.0) hlogic=2
  if(index(aa,'Force constants in Cartesian coordinates:').ne.0) hlogic=3
  if(nstrings==2) hlogic=4 ! not 100% safe
enddo
999 close(11)
if(hlogic==1) write(stdout,*) ' Found TM hessian'
if(hlogic==2) write(stdout,*) ' Found ORCA hessian'
if(hlogic==3) write(stdout,*) ' Found Gaussian09 hessian'
if(hlogic==4) write(stdout,*) ' Found PSI4 hessian'

if(hlogic<0) call error('Could not identify hessian type')
end
