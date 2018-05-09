! hydrogen bonding analysis for a cluster
subroutine hbonds(nat,iat,xyz,bond,do_assign)
use fiso, only: r8,stdout
use internals, only: frag_hval,frag_hcast,frag_nh
use constant, only: au2ang
implicit none
integer nat,iat(nat),i,j,k
character(2) esym
real(r8) ang,anggrad,xyz(3,nat)
integer acc(6),bond(nat,nat)
real(r8) rthr,athr,a1,a2,rYX,rYH
real(r8) dbond,r
logical do_assign

! distance and angle criteria (simple!)
rthr=3.0/au2ang ! Ang
athr=35_r8 ! deg D-H--A = 145.0 to 215.0
a1=180_r8-athr
a2=180_r8+athr
! print*,'thresholds:'
! print'(a,F8.1)',     ' max  Hbond length [A]   :',rthr
! print'(a,2(F8.1,1x))',' Hbond angle range [deg] :', a1,a2
! print*,''

frag_nh=0
! list of active sites: N,O,F,P,S,Cl
acc(1)=7
acc(2)=8
acc(3)=9
acc(4)=15
acc(5)=16
acc(6)=17
i=0;j=0;k=0
do i=1,nat ! Y
 if(.not.any(acc==iat(i))) cycle
 do j=1,nat ! H
   if(i==j) cycle 
   if(bond(i,j)/=1) cycle
   if(iat(j)/=1) cycle 
   do k=1,nat  ! X
     if(k==i.or.k==j) cycle
     if(.not.any(acc==iat(k))) cycle
       r=dbond(xyz(:,j),xyz(:,k))
       rYX=dbond(xyz(:,i),xyz(:,k))
       rYH=dbond(xyz(:,i),xyz(:,j))
       call angle(xyz,i,j,k,ang,anggrad)
       if(r<=rthr.and.anggrad>=a1.and.anggrad<=a2) then
        frag_nh=frag_nh+1
        if(do_assign) then
          frag_hval(frag_nh)=r
          frag_hcast(frag_nh,1)=j
          frag_hcast(frag_nh,2)=k
          bond(j,k)=2 ! padding bond matrix
          bond(k,j)=2 ! padding bond matrix
        endif
       endif
   enddo
 enddo
enddo

if(do_assign) write(stdout,'(2x,a,I4)') ' H-bonds found ',frag_nh

end subroutine


subroutine connect_frag(nat,iat,xyz,bond)
! goal: connect closest heavy-atom pair between fragment
! for now connect first atom of every fragment
use fiso, only: r8,stdout
use internals, only: frag_bval, frag_bcast, frag_nb
use logic, only: debug
implicit none
integer i,j,k,l,nat,iat(nat)
real(r8) xyz(3,nat)
integer bond(nat,nat)
! fragment data
integer, parameter :: maxfrag=5000
integer :: ifrag(nat,maxfrag),idx  !frag info array
integer :: idn,idf  ! frag_atom_index, frag_index
integer :: flen(maxfrag)       ! # atoms per fragment
integer :: iidn(maxfrag)       ! # atom index per fragment
integer :: nfrag  ! number of fragments
logical assigned
character(255) fmt
integer fpairs,ii,jj,io
real(r8) rab, dbond
integer iheavy,jheavy

ifrag=0
flen=0
iidn=0
nfrag=0
! initial fragment = atomnr 1
! flen(1)=1 
! ifrag(flen(1),1)=1

do i=1,nat
assigned=.false.

 do idf=1,nfrag
  do k=1,flen(idf)
   idx=ifrag(k,idf)
    if(bond(i,idx)==1.and..not.assigned) then
      assigned=.true.
      flen(idf)=flen(idf)+1
      ifrag(flen(idf),idf)=i
    endif
  enddo
 enddo

if(.not.assigned) then
 nfrag=nfrag+1
 iidn(nfrag)=1
 idn=iidn(idf)
 if(nfrag>maxfrag) then
    write(stdout,'(a,I6)')' increase maxfrag =',maxfrag
    stop 
 endif
 flen(nfrag)=1
 ifrag(idn,nfrag)=i
endif

enddo


write(stdout,'(2x,a,I4)') ' fragments found ',nfrag
open(newunit=io,file='xopt.intcoord')
do i=1,nfrag
  write(io,'(2x,a,I4)')'Fragment: ',i
  write(io,'(2x,a,I4)')'    size: ',flen(i)
  write(io,'(2x,a)')   '   atoms: '
  write(fmt,'("(3x",I5,1x,"I5)")') flen(i)
  write(io,fmt) ifrag(1:flen(i),i)
enddo
close(io)

!write(stdout,'(2x,a,I3,a)') ' making ',fpairs,' inter-fragment bonds (excluding H-bonds) '


fpairs=(nfrag*(nfrag-1))/2
allocate(frag_bcast(fpairs,2),frag_bval(fpairs))
idx=0

do i=1,nfrag-1
  do j=i,nfrag
   if(i==j) cycle
    ii=ifrag(1,i)
    jj=ifrag(1,j)
! check if bound already
    do k=1,flen(i)
        do l=1,flen(j)
         if(bond(ifrag(k,i),ifrag(l,j))>0) goto 111
        enddo
    enddo
! find heaviest atom in fragment
    iheavy=ii
    do k=1,flen(i)
    if (iat(ifrag(k,i))>iheavy) iheavy=ifrag(k,i)
    enddo
    jheavy=jj
    do l=1,flen(i)
    if (iat(ifrag(l,i))>iheavy) jheavy=ifrag(l,i)
    enddo
    ! print*,iheavy,jheavy
!    if(bond(ii,jj)==2) cycle ! exclude already assigned H-bonds
    ! if(bond(ii,jj)>0) cycle ! exclude if already bonded
    idx=idx+1
    frag_bcast(idx,1)=ii
    frag_bcast(idx,2)=jj
    bond(ii,jj)=3 ! padding
    bond(jj,ii)=3 
    frag_bval(idx)=dbond(xyz(1,ifrag(1,i)),xyz(1,ifrag(1,j)))
    111 continue
  enddo
enddo
frag_nb=idx
write(stdout,'(2x,a,I10,a)') ' inter-fragment bonds: ',frag_nb
write(stdout,'(2x,a,I10,a)') ' fragment pairs      : ',fpairs

end subroutine
