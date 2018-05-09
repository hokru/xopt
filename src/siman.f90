! simple simmulated annealing
subroutine siman
implicit none
real(8) Eold,Enew,energy,randr,ran
integer i,nsteps,Temp,tstep
real(8) fitold(nfit)
real(8), external :: frms
character(200) aa
logical da

call init_random_seed()

Temp=500
nsteps=250
tstep=5
!Temp=sa_temp
!nsteps=sa_mc
!tstep=sa_step

write(stdout,'(a)') ' '
write(stdout,'(a)') ' '
write(stdout,'(a)') ' ** SIMULATED ANNEALING OF NORMAL MODES **'
write(stdout,'(a)') '       H.Kruse unpublished '
write(stdout,'(a)') ' '
write(stdout,'(2x,a,I4)') 'Temperature: ',temp
write(stdout,'(2x,a,I4)') 'Temp step  : ',tstep
write(stdout,'(2x,a,I4)') 'MC steps   : ',nsteps


call status1(timer)
call message_head('calculating initial fitness ')
ancold=anc
call getenergy()
Eold=energy
call status2(timer)


write(stdout,'(2x,a)')'|   ---    MC steps done    ---        |'
!main loop (Temp)
do while (Temp>1)
  if(da('STOP')) stop '** STOP FOUND **'
  do i=1,nsteps
  if(int(mod(i,nint(nsteps/10.0))).eq.0) write(stdout,'(2x,I3,A,$)') nint(100d0*i/dble(nsteps)),' '

    ! random perturbation of normal coordinates
    call random_vib(nvar,anc)

    ! fitness
    call wrxyz(trim(xyzfile))
    call getenergy()
    if(restrain)  call addRestrainGrad
    if (freeze) call freezeGrad
    Enew=energy
    dE=energy-Eold

    call metropolis(Eold,Enew,Temp,fitold,fitpar)
  enddo

write(stdout,'(2x,a,I4,a,F18.8,a,ES10.3)') '% --> Temp:', Temp,' fitness: ',Enew ,' dE: ',dE

!logging---------------------------------------------------


!----------------------------------------------------------

  Temp=Temp-tstep
enddo

end subroutine

subroutine metropolis(Eold,Enew,temp,nvar,ancold,anc)
implicit none
integer temp
integer, intent(in) :: nvar
real(8) Eold,Enew,boltz,ran
real(8),intent(inout) anc(nvar),ancold(nvar)

call random_number(ran)

boltz=exp(-1d0*(Enew-Eold)/temp)
if(Enew>Eold) then
  ! revert step
  Enew=Eold
  anc=ancold
  if(ran<=boltz) then
  !accept
    Eold=Enew
    ancold=anc
  endif

else
! lower fitness solution, accept by default
  Eold=Enew
  ancold=anc
endif
end subroutine


! randomly perturb vibration
! change by up to +/- 50%
subroutine random_vib(nvar,anc)
implicit none
integer i,nvar
real(8) ran,randr,r1,r2
real(8), intent(inout) :: anc(nvar)

r1=-0.5d0
r2=0.5d0
do i=1,nfit
 call random_number(ran)
 ran=(r2-r1)*ran+r1
 anc(i)=anc(i)+ran*anc(i)
enddo
end subroutine

