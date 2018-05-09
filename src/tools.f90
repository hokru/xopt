! wrapper for L-BFGS





subroutine Echange2_int(nvar,gint,hint,displ,echange)
! predict energy change according to 2nd order Taylor expansion
use fiso, only: r8
implicit none
integer :: nvar
real(r8), intent(in) :: displ(nvar),gint(nvar),hint(nvar*(nvar+1)/2)
real(r8), intent(out) :: echange
real(r8) :: ddot,gdot,hdot
real(r8) :: hd(nvar)

call dspmv('u',nvar,0.5_r8,hint,displ,1,0.0_r8,hd,1)
gdot=ddot(nvar,displ,1,gint,1)
hdot=ddot(nvar,displ,1,hd,1)
echange = hdot + gdot

end subroutine


! WRAPPERS FOR COMPILER/GNU EXTENSIONS

subroutine system(input)
! wrapper for non-F2008 standard intrinsic system
! if used with intrinsics it must be declated exernal
implicit none
character(*), intent(in):: input
integer :: stat=1
character(255) :: cmdm=''

call execute_command_line(trim(input),exitstat=stat,wait=.true.,cmdmsg=cmdm)
end subroutine


!character(255) function getcwd()
! wrapper for non-F2008 extension.
! gets $PWD from console for the current working directory.
! if used with intrinsics it must be declated exernal
!implicit none
!
!call get_environment_variable('PWD',getcwd)
!end function


! how to do this with fortran??
! subroutine chdir(input)
! implicit none
! character(255), intent(in) :: input
!
! end

subroutine check_stop()
implicit none
logical da
integer io
inquire(file='STOP',exist=da)
if(da) then
 open(newunit=io, file='STOP', status='old')
 close(io, status='delete')
 call warning('STOP FILE FOUND!')
 stop
endif
end subroutine
