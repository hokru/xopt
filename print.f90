!**************************************************************************************
!    Copyright 2013, 2014 Holger Kruse                                                *
!                                                                                     *
!    This file is part of Xopt.                                                       *
!                                                                                     *
!    Xopt is free software: you can redistribute it and/or modify                     *
!    it under the terms of the GNU Lesser General Public License as published by      *
!    the Free Software Foundation, either version 3 of the License, or                *
!    (at your option) any later version.                                              *
!                                                                                     *
!    xopt is distributed in the hope that it will be useful,                        *
!    but WITHOUT ANY WARRANTY; without even the implied warranty of                   *
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    *
!    GNU Lesser General Public License for more details.                              *
!                                                                                     *
!    You should have received a copy of the GNU Lesser General Public License         *
!    along with Xopt.  If not, see <http://www.gnu.org/licenses/>.                    *
!                                                                                     *
!**************************************************************************************
subroutine head
character(80) host,dir
!*************************
!* An eXternal OPTimizer *
!* ***********************
call get_environment_variable('HOSTNAME',host)
call get_environment_variable('PWD',dir)
print*,'             '
print*,' -------------------------- '
print*,'|                         | '
print*,'|      \_/ _  _ |_        | '
print*,'|      / \(_)|_)|_        | '
print*,'|            |            | '
print*,' ------------------------- '
print*,'             '
print*,' An eXternal OPTimizer by'
print*,'  Holger Kruse (mail2holger@gmail.com)'
print*,'             '
write(*,'($,a)') ' Date: '
call timestamp()
write(*,'(2a)') ' Host:   ',trim(host)
write(*,'(2a)') ' working directory: ',trim(dir)
print*,'             '
include 'version.dat'
print*,'             '
end subroutine head


subroutine p_info
use popt
use logic
use parm
integer io
logical da
character(120) aa
io=6

write(io,'(2x,a)') ' Molecular data:'
write(io,'(3x,a,I4)') ' Atoms   = ',nat
!call composition

! general stufff
write(io,'(2x,a)') ' Convergence criteria:'
write(io,'(3x,a,ES10.3)') ' energy    = ',econv
write(io,'(3x,a,ES10.3)') ' gnorm     = ', gconv
write(io,'(3x,a,ES10.3)') ' max grad  = ',maxgrad
write(io,'(3x,a,ES10.3)') ' max displ = ',dconv
write(io,'(a)')
write(io,'(2x,a)') ' Setup: ' 
write(io,'(3x,a,F7.4)') ' max displacement = ', maxd
write(io,'(3x,a,I4)')   ' max iter         = ',maxiter
write(io,'(a)')

! selected methods

write(io,'(2x,a,$)')      ' selected program : '
if(tm) write(io,'(5x,a,L)')   '    --> Turbomole'
if(orca) write(io,'(5x,a,L)') '    --> ORCA'
if(mopac) write(io,'(5x,a,L)')'    --> mopac12'
if(gaus) write(io,'(5x,a,L)') '    --> Gaussian09'
if(openbabel) write(io,'(5x,a,L)') '    --> openbabel'

write(io,'(2x,a,I2,$)')   ' Hessian update   : ',iupdate
if(iupdate==1)   write(io,'(3x,a)') '    --> SR1-BFGS'
if(iupdate==2)   write(io,'(3x,a)') '    --> SR1'
if(iupdate==3)   write(io,'(3x,a)') '    --> BFGS'
if(iupdate==4)   write(io,'(3x,a)') '    --> SS-BFGS'
if(iupdate==5)   write(io,'(3x,a)') '    --> SS-SR1-BFGS'

write(io,'(2x,a,I2,$)')   ' program mode     : ', mode
if(mode==1)      write(io,'(3x,a)') '    --> optimization'
if(mode>1)       write(io,'(3x,a)') '    --> transition state search'

write(io,'(2x,a,I2,$)')   ' convergence mode : ',iconv
if(iconv==4) write(io,'(3x,a)') '    --> all criteria'
if(iconv==2) write(io,'(3x,a)') '    --> energy + gradient norm'
if(iconv==1) write(io,'(3x,a)') '    --> energy only'

write(io,'(2x,a,I2,$)')   ' optimizer mode   : ',iopt
if(iopt==11) write(io,'(3x,a)') '    --> cartesian + RFO'
if(iopt==12) write(io,'(3x,a)') '    --> cartesian + SI-RFO'
if(iopt==21) write(io,'(3x,a)') '    --> approx. normal coords + RFO'
if(iopt==22) write(io,'(3x,a)') '    --> approx. normal coords + SI-RFO'
write(io,'(a)')


if(ppot.or.d3.or.gcp) then 
write(io,'(2x,a)')   ' helper & options  : '
if(gcp) write(io,'(5x,2a)')   '    --> gcp ', trim(gcplevel)
if(d3) write(io,'(5x,2a)')   '    --> dftd3 ', trim(d3func)
if(ppot) write(io,'(5x,a)')   '    --> ppot '
if(mopac) then
 inquire(file='SETUP',exist=da)
 if(da) then
  open(88,file='SETUP')
  read(88,'(a)') aa
  write(*,*) 'SETUP: ',trim(aa)
 else
  stop ' no SETUP file found!!'
endif

endif
write(io,'(a)')
endif

end subroutine p_info



!subroutine pciter(iter,e,de,gnorm,gmax,dnorm,disp,lam)
subroutine pciter(iter,e,de,gnorm,gmax,dnorm,disp,lam,EE,GG,GM,DM)
!use popt
implicit none
real(8) e,de,disp,gmax,lam,dnorm,gnorm
integer iter,io
logical EE,GG,GM,DM

io=6
!write(io,'(2x,i4,F18.8,2x,6(ES10.3,2x),2L)') iter,e,de,gnorm,gmax,dnorm,disp,lam,GG,EE
!write(io,'(2x,i4,F18.8,2x,ES10.3,1x,L,2x,ES10.3,1x,L,2x,ES10.3,1x,L,2x,F10.3,2x,F10.3,1x,L,2x,ES10.3)') iter,e,de,EE,gnorm,GG,gmax,GM,dnorm,disp,DM,lam
write(io,'(2x,i4,F18.8,2x,ES10.3,1x,L,2x,F8.5,1x,L,2x,ES10.3,1x,L,2x,F10.3,2x,F9.4,1x,L,2x,ES10.3)') iter,e,de,EE,gnorm,GG,gmax,GM,dnorm,disp,DM,lam
end subroutine pciter


subroutine endopt(dE,gnorm,gmax,dnorm)
use popt
implicit none
real(8), intent(in) :: dE,gnorm,gmax,dnorm
print*,'CONVERGED !'
write(6,'(3x,a)') '               criteria   actual value'
write(6,'(3x,a,2(ES10.3,2x))') ' energy    = ',econv,dE
write(6,'(3x,a,2(ES10.3,2x))') ' gnorm     = ',gconv,gnorm
write(6,'(3x,a,2(ES10.3,2x))') ' max grad  = ',maxgrad,gmax
write(6,'(3x,a,2(ES10.3,2x))') ' max displ = ',dconv,dnorm
end subroutine

