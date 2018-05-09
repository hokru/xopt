! intrinsic-reaction coordinate driver
subroutine irc_driver()
! use parm, only: iat,nat,xyz
use popt, only: do_dvv
use fiso, only: r8,stdout
! real(r8) :: mxyz(3,nat),mgrad(3,nat)

write(stdout,'(a)') ''
write(stdout,'(a)') ' Reaction Coordinate Driver'
write(stdout,'(a)') ''



if(do_dvv) call irc_dvv()
    

! call get_mxyz(mxyz)
! call get_mgrad(mgrad)



end subroutine


subroutine irc_dvv()
! damped velocity verlect
! dynamic reaction path integrator
use parm, only: iat,nat,xyz,energy,grad,gnorm,chess
use progs, only: xyzfile
use popt, only: dvv_err, dvv_init,dvv_step,tsmode,dvv_de,dvv_up
use mddat, only: maxstep,dt
use fiso, only: r8,stdout
use constant, only: au2fs,au2kcal,au2m,amu2kg,hartree2J
! use logic, only: debug
implicit none
integer i,j,k,iter,c
real(r8) :: v(3,nat),a(3,nat)
real(r8) :: mxyz(3,nat),mgrad(3,nat),vnorm
real(r8) :: rdist,r(3,nat),rmax,de,eold
real(r8) :: am2(3,nat,2),xm2(3,nat,2),vm2(3,nat,2),xp(3,nat),dtm1,err
real(r8) :: e(nat*3),freqcm,dtnew,max_dt,min_dt

! gradient to acceleration (bohr amu**0.5/fs**2)
real(r8), parameter :: grad2acc=1.0_r8/(1d15*1d15*au2m**2*amu2kg/hartree2J)
    
write(stdout,'(a)') '  '    
write(stdout,'(a)') '  damped velocity verlet (DVV) integrator'
write(stdout,'(a)') '  Hratchian/Schlegel: 10.1021/jp012125b'
write(stdout,'(a)') '  '
write(stdout,'(2x,a)') 'DVV settings  '
write(stdout,'(2x,a,F7.3)') ' error threshold  :', dvv_err
write(stdout,'(2x,a,F7.3)') ' velocity damping :', dvv_step
write(stdout,'(2x,a,F7.3)') ' initial step     :', dvv_init
write(stdout,'(2x,a,E10.3)') ' de threshold     :', dvv_de
write(stdout,'(2x,a,E10.3)') ' up threshold     :', dvv_up
write(stdout,'(a)') '  '


! init
call wrxyz('xopt.log')
call touchlog

v=0.0_r8;vm2=0.0_r8;xm2=0.0_r8
dt=dt
max_dt=2.0_r8/au2fs
min_dt=0.025_r8/au2fs
dvv_step=dvv_step

dt=dt/au2fs

! first energy+gradient
call getgrad

! mass-weight
call get_mxyz(mxyz,xyz,nat,iat,1)
call get_mgrad(mgrad)

! initial acceleration
a=-mgrad*grad2acc
call get_acc(grad,a)

! set initial velocity from Hessian eigenvector
! project? mass-weight?
write(stdout,'(2x,a)') 'generating initial path velocity from Hessian eigenvector'
call Hmass(chess,'project')
call DiagSM(nat*3,chess,e)
write(stdout,'(2x,a)') 'vib. freq (cm-1)'
call printvib(nat*3,e)
write(stdout,'(2x,a,I2,a,F10.2,a)') 'Following ',tsmode,' with ',freqcm(e(tsmode)),' cm-1 '


k=0
do i=1,nat
 do j=1,3
    k=k+1
    v(j,i)=chess(k,tsmode)
 enddo
enddo
! call get_mxyz(v,v,nat,iat,1)  ! <== ??s

! direction 
v=v*dvv_init

! initial damping
call cartgnorm(nat,v,vnorm)
v=v*(dvv_step/vnorm)
! call cartgnorm(nat,v,vnorm)
! print*,'initial |velocity|', vnorm


 iter=0
 c=0
print*, '' 
! <<<< IRC-DVV LOOP >>>>>
write(stdout,'(a)') '    iter    dt       ddt             E            error         dE'
do
    call check_stop()
    iter=iter+1
    eold=energy
    do i =1,nat
        do j=1,3
            v(j,i)=v(j,i) + a(j,i)*0.5_r8*dt
            mxyz(j,i)=mxyz(j,i) + v(j,i)*dt + a(j,i)*0.5_r8*dt*dt
        enddo
    enddo

    !  print*, mxyz(1,1),xm2(1,1,2),xm2(1,1,1)
    ! call get_coord_dist(nat,iat,mxyz,xm2(:,:,1),rdist,0)
    ! print*, rdist
    call get_mxyz(mxyz,xyz,nat,iat,0)
    call wrxyz(trim(xyzfile))
    call getgrad()
    call loggrad(nat,iter,grad,energy)
    call cartgnorm(nat,grad,gnorm)
    call appxyz('xopt.log')
    call get_mgrad(mgrad)
    de=energy-eold

    a=-mgrad*grad2acc
    ! call get_acc(grad,a)
    do i =1,nat
        do j=1,3
            v(j,i)=v(j,i)+a(j,i)*0.5_r8*dt
        enddo
    enddo
    
    ! damp
    call cartgnorm(nat,v,vnorm)
    v=v*(dvv_step/vnorm)


    if(iter<=2) write(stdout,'(I6,2x,F8.3,2x,8x,F18.7,2x,9x,2x,ES10.3)') iter,dt*au2fs,energy,de
 
    
    if(iter>2) then
        if(de>dvv_up) then
        ! call warning('Walking uphill')
        write(stdout,'(2x, a)') 'walking uphill! stopping'
        !  call write_irc_restart(xyz,v)
        exit
        elseif(abs(de)<dvv_de) then
        write(stdout,'(2x, a)') 'very small step! stopping'
        exit
        endif


        ! error estimate
        do i=1,nat
          do j=1,3
          xp(j,i)=xm2(j,i,1)+vm2(j,i,1)*(dtm1+dt)+0.5_r8*am2(j,i,1)*(dtm1+dt)**2
          enddo
          
        enddo
        ! print*, xm2(1:3,1,1),xp(1:3,1)
        !  call get_mxyz(mxyz,xyz,nat,iat,0)
        ! call localxyz('xp.xyz',nat,iat,xyz)
        ! stop

        r=0_r8
        rmax=0_r8
        call get_coord_dist(nat,iat,xp,mxyz,rdist,0)
!         do k=1,nat
!           do j=1,3
!             r(j,k)=xp(j,k)-mxyz(j,k)
! !            if(r(j,k)>=rmax) rmax=r(j,k)
!           enddo
!         enddo
!         call get_mxyz(r,r,nat,iat,0)
!         call cartgnorm(nat,r,rdist)
        err=rdist
        
        ! dtnew=min_dt
        dtnew=dt*(dvv_err/err)**(1.0_r8/3.0_r8)
        ! limit
        if(dtnew>=max_dt) dtnew=max_dt
        if(dtnew<=min_dt) dtnew=min_dt
        write(stdout,'(I6,2x,F8.3,2x,F8.3,F18.7,2x,ES9.2,2x,ES10.3)') iter,dt*au2fs,(dt-dtnew)*au2fs,energy,err,de
        dt=dtnew
    endif
    

    
    dtm1=dt
    ! safe previous steps
    if(iter<=2) then
        vm2(1:3,1:nat,iter)=v(1:3,1:nat)
        am2(1:3,1:nat,iter)=a(1:3,1:nat)
        xm2(1:3,1:nat,iter)=mxyz(1:3,1:nat)
    endif
    if(iter>2) then
      vm2(1:3,1:nat,1)=vm2(1:3,1:nat,2)
      am2(1:3,1:nat,1)=am2(1:3,1:nat,2)
      xm2(1:3,1:nat,1)=xm2(1:3,1:nat,2)
      vm2(1:3,1:nat,2)=v(1:3,1:nat)
      am2(1:3,1:nat,2)=a(1:3,1:nat)
      xm2(1:3,1:nat,2)=mxyz(1:3,1:nat)
    endif
      
    

if(iter==maxstep) exit
enddo



end subroutine


subroutine get_coord_dist(nat,iat,xyz1,xyz2,rdist,iflag)
! distance between mass-weighted xyz1 and xyz2
use fiso, only: r8
implicit none
integer :: nat,j,k,iat(nat),iflag
real(r8) :: xyz1(3,nat),xyz2(3,nat),r(3,nat),rdist
     if(iflag>0) then
     call get_mxyz(xyz1,xyz1,nat,iat,0)
     call get_mxyz(xyz2,xyz2,nat,iat,0)
     endif
     do k=1,nat
          do j=1,3
            r(j,k)=xyz1(j,k)-xyz2(j,k)
          enddo
        !   print*,xyz1(1:3,k)
     enddo
        call get_mxyz(r,r,nat,iat,0)
        call cartgnorm(nat,r,rdist)
return
end subroutine

subroutine get_mxyz(mxyz,xyz,nat,iat,iflag)
! use parm, only: iat,nat,xyz
use fiso, only: r8,stdout
use atomdata, only: ams2
use constant, only: amu2au
integer :: i,iflag,nat,iat(nat)
real(r8),intent(inout) :: mxyz(3,nat)
real(r8),intent(inout) :: xyz(3,nat)
if(iflag==1) then
    do i=1,nat
        mxyz(1:3,i)=xyz(1:3,i)*sqrt(ams2(iat(i))*amu2au)
    enddo
elseif(iflag==0) then
    do i=1,nat
        xyz(1:3,i)=mxyz(1:3,i)/sqrt(ams2(iat(i))*amu2au)
    enddo
endif
end subroutine


subroutine get_mgrad(mgrad)
use parm, only: iat,nat,grad
use fiso, only: r8,stdout
use atomdata, only: ams2
use constant, only: amu2au
integer :: i
real(r8),intent(inout) :: mgrad(3,nat)
do i=1,nat
 mgrad(1:3,i)=grad(1:3,i)/sqrt(ams2(iat(i))*amu2au)
enddo
end subroutine


subroutine get_acc(g,a)
use parm, only: iat,nat
use fiso, only: r8,stdout
use atomdata, only: ams2
use constant, only: amu2au,au2m,amu2kg,hartree2J
real(r8), parameter :: grad2acc=1.0_r8/(1d15*1d15*au2m**2*amu2kg/hartree2J)
integer :: i
real(r8),intent(in) :: g(3,nat)
real(r8),intent(out) :: a(3,nat)
a=0d0
do i=1,nat
 a(1:3,i)=-g(1:3,i)*grad2acc/sqrt(ams2(iat(i)))
enddo
end subroutine

