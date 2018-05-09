! all things molecular dynamics


subroutine runMD
use fiso, only: r8,stdout
use MDdat
use popt, only: restart,freeze,restrain
use progs, only: xyzfile
use parm, only: xyz, energy, grad,nat,nvar,gnorm,enetot
use constant, only: au2fs,kb_au,au2kcal
implicit none
integer step,sstep
integer i,j
real(r8) etk,temp,time,tav,vol_r,p
logical da


allocate(velo(3,nat),mass(nat))
nvar=3*nat
velo=0.0_r8
p=0.0_r8

! some defaults
sstep=1
sstep=1


! thermo stuff
do_thermo=.true.
tau=1.0_r8
m_nh=nvar*temp0*kb_au*(dt*30)**2
m_nh=m_nh/au2fs !for dt
m_nh=0.1_r8*1000_r8/au2fs
!print*,'Q NHC ',m_nh
c_nh=0.0_r8
p_nh=0.0_r8
v_nh=0.0_r8



! units
dt=dt/au2fs
tau=tau/au2fs
!*************
! initialize *
!*************
call setMass
if(restart) then
 call MDrestart
else
 call randomVelo
endif
call getETK(etk)
call getTemp(etk,temp)

write(stdout,'(a)') 'initial conditions:'
write(stdout,'(1x,a,3x,F8.2,1x,F8.2)') 'Temp(instant)/Temp(target)', temp,temp0

grad=0.0_r8
! initial energy/forces
call getgrad

!****************
! START MD loop *
!****************
open(555,file='xopt.md')
call touchlog
write(stdout,'(a)') ''
write(stdout,'(a)') ''
write(stdout,'(a)') '  **  Velocity-Verlet MD  **'
!write(stdout,'(a)') '    step    time (fs)           E           Temp (avg)     '
write(stdout,'(a)') '    step    time (fs)           E        Temp    T(avg)    E_r [kcal]      P '
tav=0.0_r8
do step=sstep,maxstep
tav=tav+temp
![+periodically remove trans/rot]
![+translate to COM]

call rmfile('STOP',2)
inquire(file='STOP',exist=da)
if(da) then
 open(unit=1234, file='STOP', status='old')
 close(1234, status='delete')
 stop 'user requested stop'
endif


! apply thermostat
if(do_thermo)  call thermo_select(dt,etk)

! velo and positions
!call MD_VVstep(1,dt)
! mid-step
do i=1,nat
 do j=1,3
  velo(j,i)=velo(j,i)+(-grad(j,i)/mass(i))*dt*0.5_r8
  xyz(j,i)=xyz(j,i)+velo(j,i)*dt
 enddo
enddo

! get energy/gradient + write log
call wrxyz(trim(xyzfile))
call getgrad
if(restrain) call addRestrainGrad
if (freeze) call freezeGrad
call cartgnorm(nat,grad,gnorm)
call appxyz('xopt.log')

! final velo
!call MD_VVstep(2,dt)
do i=1,nat
 do j=1,3
  velo(j,i)=velo(j,i)+(-grad(j,i)/mass(i))*dt*0.5d0
 enddo
enddo


call getETK(etk)
call getTemp(etk,temp)

! thermostat
if(do_thermo)  call thermo_select(dt,etk)


! print runtime info
time=dble(step)*dt*au2fs
!write(stdout,'(I6,2x,F8.2,2x,F18.6,2x,F7.2,2x,F7.2)') step,time,energy+etk,temp,tav/step
write(stdout,'(I6,2x,F8.2,2x,F18.6,2x,F7.2,2x,F7.2,4x,F7.2,6x,F7.2)') step,time,energy+etk,temp,tav/step,enetot*au2kcal,p
write(555,'(I6,2x,F8.2,2x,F18.6,2x,F7.2,2x,F7.2)') step,time,energy+etk,temp,tav/step


! write restart file
call wrestart(step)

enddo
close(555)
!***************
!* END MD LOOP *
!***************
if(restrain) call printFinalRestrain

end subroutine


subroutine thermo_select(dt,etk)
use fiso, only: r8
use MDdat, only: thermo_nh, thermo_berend,thermo_scale
implicit none
real(r8), intent(in) :: dt
real(r8), intent(inout) :: etk

if(thermo_nh) then
 call thermostat_NH(dt,etk)
elseif(thermo_berend) then
  call thermostat_berendsen(dt,etk)
elseif(thermo_scale) then
 call thermostat_scale(etk)
else
 call error('no thermostat selected!')
endif

end subroutine

! set masses in au
subroutine setMass
use parm, only: i,nat,iat
use MDdat, only: mass
use atomdata, only: ams
use constant, only: amu2au
implicit none

do i=1,nat
   mass(i)=ams(iat(i))*amu2au
enddo
end subroutine

! normal distributed velocities for given temperature
subroutine randomVelo
use fiso, only: r8,stdout
use constant,only: kb_au
use MDdat,   only: Temp0,velo,mass
use parm,    only: nat,nvar
implicit none
integer i,j
real(r8) x,aveEkin
logical gauss
real(r8) t,etk

call init_random_seed()
velo=0.0_r8

gauss=.true.
!gauss=.false.

aveEkin=0.5_r8*kb_au*Temp0

if(gauss) then
  do i=1,nat
   do j=1,3
    call randGaus(x)
    velo(j,i) = sqrt(2.0_r8*aveEkin/mass(i))*x
   end do
  end do
  call getETK(etk)
  call thermostat_scale(etk)
  call getETK(etk)
  call getTemp(etk,t)
  write(stdout,'(a)') 'initial velocity: '
  write(stdout,'(2x,a,F8.2)') 'randomized & scaled gaussian distribution with T=',t
else
 do i=1,nat
  do j=1,3
   call random_number(x)
   x=x*(0.5_r8+0.5_r8)-0.5_r8
   velo(j,i) = sqrt(2*aveEkin/mass(i))*x
  end do
 end do
 call getETK(etk)
 t=2.0_r8*etk/(nvar*kb_au)
 write(stdout,'(a,F8.2)')'Generated uniform velocity distribution [-0.5/0.5] with T=',t
endif

end

! calculate instantanous temperature temp at given kinetic energy
subroutine getTemp(etk,temp)
use fiso, only: r8
use constant, only: kb_au
use parm, only: nvar
implicit none
real(r8) temp,etk
 Temp=2.0_r8*etk/(kB_au*nvar)
end subroutine


! instantaneous internal kinetic energy etk
subroutine getETK(etk)
use MDdat, only: velo, mass
use parm  ,only: nat
use fiso, only: r8
implicit none
integer i
real(r8) etk,vv1,vv2,vv3

etk=0.0_r8
do i=1,nat
 vv1=velo(1,i)**2
 vv2=velo(2,i)**2
 vv3=velo(3,i)**2
 etk=etk+mass(i)*(vv1+vv2+vv3)
enddo
etk=etk*0.5_r8
return
end


! Velocity Verley integrator
subroutine MD_VVstep(switch,dt)
use fiso, only: r8
use parm, only: nat, grad,xyz
use MDdat, only: velo, mass
implicit none
integer i,j,switch
real(r8) a(3,nat)
real(r8) dt,dt2

! acceleration
!m=m*amu2au
a(1,1:nat)=-grad(1,1:nat)/mass(1:nat)
a(2,1:nat)=-grad(2,1:nat)/mass(1:nat)
a(3,1:nat)=-grad(3,1:nat)/mass(1:nat)


dt2=dt*dt

if(switch==1) then
! dt/2 step
do i =1,nat
 do j=1,3
  velo(j,i)=velo(j,i)+a(j,i)*0.5_r8*dt
  xyz(j,i)=xyz(j,i)+velo(j,i)*dt+a(j,i)*0.5_r8*dt2
 enddo
enddo
elseif(switch==2) then
do i =1,nat
 do j=1,3
  velo(j,i)=velo(j,i)+a(j,i)*0.5_r8*dt
!  xyz(j,i)=xyz(j,i)+velo(j,i)*dt
 enddo
enddo
endif

end subroutine


!***************
!* THERMOSTATS *
!***************
! Berendsen thermostat
subroutine thermostat_berendsen(dt,etk)
use MDdat, only: tau,velo,temp0
use parm, only: nat
use fiso, only: r8
!use constant, only: kb_au
implicit none
integer i,j
real(r8) dt,etk
real(r8) fscal,itemp,ratio,dtt


call getTemp(etk,itemp)

ratio=(Temp0/iTemp-1.0_r8)
dtt=dt/tau*2.0_r8
fscal=sqrt(1.0_r8+dtt*ratio)

! scale velocities
do i=1,nat
 do j=1,3
  velo(j,i)=velo(j,i)*fscal
 enddo
enddo
end subroutine

! simple velocity rescaling
subroutine thermostat_scale(etk)
use MDdat, only: velo,temp0
use parm, only: nat
use fiso, only: r8
implicit none
integer i,j
real(r8) etk
real(r8) fscal,itemp,ratio

! instant. temperature
call getTemp(etk,itemp)

ratio=(Temp0/iTemp)
fscal=sqrt(ratio)

! scale velocities
do i=1,nat
 do j=1,3
  velo(j,i)=velo(j,i)*fscal
 enddo
enddo

end subroutine


! Nose-Hoover thermostat following Frenkel&Smith 
subroutine thermostat_NH(dt,etk)
use MDdat, only: velo,temp0,c_nh,m_nh,v_nh,p_nh
use parm, only: nvar
use constant, only: kb_au
use fiso, only: r8
implicit none

real(r8) dt2,dt4,dt8,dt,etk,fscal,T

!T=temp0
T=temp0*kb_au

dt2=0.5_r8*dt
dt4=0.25_r8*dt
dt8=0.125_r8*dt

c_nh(2)=(m_nh(1)*v_nh(1)*v_nh(1)-T)/m_nh(2)
!c_nh(2)=(m_nh(1)*v_nh(1)*v_nh(1)-T)
v_nh(2)=v_nh(2)+c_nh(2)*dt4
v_nh(1)=v_nh(1)*exp(-v_nh(2)*dt8)

c_nh(1)=(2d0*etk-nvar*T)/m_nh(1)
v_nh(1)=v_nh(1)+c_nh(1)*dt4
v_nh(1)=v_nh(1)*exp(-v_nh(2)*dt8)
p_nh(1)=p_nh(1)+v_nh(1)*dt2
p_nh(2)=p_nh(2)+v_nh(2)*dt2

fscal=exp(-v_nh(1)*dt2)

! scale Ekin
etk=etk*fscal*fscal

! scale vecolities
velo=velo*fscal

v_nh(1)=v_nh(1)*exp(-v_nh(2)*dt8)
c_nh(1)=(2.0_r8*etk-nvar*T)/m_nh(1)
v_nh(1)=v_nh(1)+c_nh(1)*dt4
v_nh(1)=v_nh(1)*exp(-v_nh(2)*dt8)
c_nh(2)=(m_nh(1)*v_nh(1)*v_nh(1)-T)/m_nh(2)
v_nh(2)=v_nh(2)+c_nh(2)*dt4

return
end subroutine


subroutine thermostat_NH2(dt,etk)
! 2-chain NHC after Frenkel book
use MDdat, only: velo,temp0,m_nh,v_nh,g_nh,r_nh
use parm, only: nvar
use constant, only: kb_au
use fiso, only: r8
implicit none

real(r8) dt2,dt4,dt8,dt,etk
real(r8) kbt,ke_in,scalefac

ke_in=etk
kbT=temp0*kb_au
dt2=0.5_r8*dt
dt4=0.25_r8*dt
dt8=0.125_r8*dt

!      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kbt)/m_nh(2)
      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kbt)
      v_nh(2) = v_nh(2) + g_nh(2)*dt4
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
      g_nh(1) = (2.0_r8*KE_in - nvar*kbt)/m_nh(1)
      v_nh(1) = v_nh(1) + g_nh(1)*dt4
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
      r_nh(1) = r_nh(1) + v_nh(1)*dt2
      r_nh(2) = r_nh(2) + v_nh(2)*dt2

      scalefac = dexp(-v_nh(1)*dt2)
      KE_in = KE_in*scalefac*scalefac

      g_nh(1) = (2.0d0*KE_in - nvar*kbt)/m_nh(1)
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
      v_nh(1) = v_nh(1) + g_nh(1)*dt4
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kbt)/m_nh(2)
      v_nh(2) = v_nh(2) + g_nh(2)*dt4

velo=velo*scalefac

return
end subroutine


subroutine get_p_virial(temp,p)
! virial pressure
! P=1/V ( 1/3 sum_i (m_iv_i*v_i) + 1/3 sum_i (r_i * f_i))
use parm, only: nat, grad,xyz,molvol
use fiso, only: r8
use constant, only: kb_au
use MDdat, only: velo, mass
implicit none
integer i,j
real(r8), intent(out):: p
real(r8) :: s,rij,V,temp


s=0.0_r8
do i=1,nat
 do j=1,3
 s=1
 enddo
enddo
p=s/(3.0_r8*molvol)

end subroutine
