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
subroutine eval_options (infile,outfile)
use logic
use popt
implicit none
integer i,maxarg
character(80), allocatable :: arg(:)
character(120) infile,outfile
character(80) ftmp,aa,narg
integer s2i
real(8) s2r


 call getarg(1,infile)
 infile=trim(infile)
 if(infile=='') stop 'no structure file!!'
 if(infile=='-h'.or.infile=='-help') stop 'help..'!call help

maxarg=iargc()
if(maxarg.gt.0) then

 allocate(arg(maxarg))
   do i=1,maxarg
     call getarg(i,arg(i))
   enddo

 do i=2,maxarg
  ftmp=arg(i)
  if(i/=maxarg)narg=trim(arg(i+1))
  if(index(ftmp,'-h ').ne.0) then
   print*,' help!..Help? HELP!!! ...h e l p ? !...  .'
  stop
  endif
!  if(index(ftmp,'-noprint ').ne.0) echo=.false.
  ! <<<<<         PROGS              >>>>>>
  if(index(ftmp,'-tm ').ne.0) tm=.true.
  if(index(ftmp,'-readonly ').ne.0) readonly=.true.
  if(index(ftmp,'-tmhuge').ne.0) tmhuge=.true.
  if(index(ftmp,'-g09 ').ne.0) gaus=.true.
  if(index(ftmp,'-mop ').ne.0) mopac=.true.
  if(index(ftmp,'-mopac ').ne.0) mopac=.true.
  if(index(ftmp,'-orca ').ne.0) orca=.true.
  if(index(ftmp,'-obabel ').ne.0) openbabel=.true.
  if(index(ftmp,'-relax ').ne.0) relax=.true.
  if(index(ftmp,'-driver ').ne.0) driver=.true.
  if(index(ftmp,'-ppot ').ne.0) then
   ppot=.true.
   tm=.true.
  endif
  if(index(ftmp,'-apbs ').ne.0) APBS=.true.
  if(index(ftmp,'-oniom ').ne.0) qmmm=.true.
  if(index(ftmp,'-amber ').ne.0) then
    amber=.true.
    maxd=0.01
    maxiter=1000
  endif
  if(index(ftmp,'-rhess ').ne.0) then
   readHess=.true.
  endif
  if(index(ftmp,'-restart ').ne.0) then
   restart=.true.
  endif
  if(index(ftmp,'-d3 ').ne.0) then
    d3=.true.
    d3func=trim(arg(i+1))
  endif
  if(index(ftmp,'-gcp ').ne.0) then
    gcp=.true.
    gcplevel=trim(arg(i+1))
  endif
  if(index(ftmp,'-numhess ').ne.0) hmodel='numerical'

  ! <<<<<         Minimizer              >>>>>>
  if(index(ftmp,'-cart ').ne.0) iopt=11
  if(index(ftmp,'-cart-si ').ne.0) iopt=12
  if(index(ftmp,'-anc ').ne.0) iopt=21
  if(index(ftmp,'-anc-si ').ne.0) iopt=22
  if(index(ftmp,'-anc-cg ').ne.0) iopt=23
  if(index(ftmp,'-cgrfo ').ne.0) iopt=24
  if(index(ftmp,'-int ').ne.0) iopt=31
  if(index(ftmp,'-sg1-bfgs ').ne.0) iupdate=1
  if(index(ftmp,'-sg1 ').ne.0) iupdate=2
  if(index(ftmp,'-bfgs ').ne.0) iupdate=3
  if(index(ftmp,'-ss-bfgs ').ne.0) iupdate=4
!  if(index(ftmp,'-ss-sg1-bfgs ').ne.0) iupdate=5

  ! <<<<<         Setup             >>>>>>
  if(index(ftmp,'-c ').ne.0) then ! max iterations
    maxiter=s2i(narg)
  endif
  if(index(ftmp,'-iconv ').ne.0) then ! max iterations
    iconv=s2i(narg)
  endif
  if(index(ftmp,'-dpl ').ne.0) then ! max displacement
    maxd=s2r(narg)
  endif
  if(index(ftmp,'-gconv ').ne.0) then ! gnorm threshold
    gconv=s2r(narg)
  endif
  if(index(ftmp,'-econv ').ne.0) then ! energy threshold
    econv=s2r(narg)
  endif
  ! <<<<<         Compound keywords           >>>>>>
  if(index(ftmp,'-tight ').ne.0) tightopt=.true.
  if(index(ftmp,'-loose ').ne.0) preopt=.true.



  ! <<<<<         special keywords           >>>>>>
  if(index(ftmp,'-amber-pb ').ne.0) then
   print*, ' writing  amber.in for MMPB optimizations'
   call GetAmberinput('pb')
    amber=.true.
    maxd=0.01
    maxiter=1000
  endif
  if(index(ftmp,'-amber-gb ').ne.0) then
   print*, ' writing  amber.in for MMGB optimizations'
   call GetAmberinput('gb')
    amber=.true.
    maxd=0.01
    maxiter=1000
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


end subroutine



