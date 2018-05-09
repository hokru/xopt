! internal coordinate routines
! We need a set of primitives for the model Hessians (get_primitives)
! But for the actual coordinate transformation we need more (make_intcoords)

subroutine get_primitives(nat,iat,xyz)
! use constant, only: au2ang
use fiso, only: r8,stdout
use logic, only: debug,large
use internals
implicit none
integer, intent(in) :: nat,iat(nat)
real(r8), intent(in) :: xyz(3,nat)
integer, allocatable :: bond(:,:),cn(:)
integer :: kk,i,j,k,l
real(r8) :: rab,ang,anggrad,dbond,agrad,di360,word2mb
integer :: ii,jj,ll
integer :: it
character(255) :: file_out
character(2) :: esym
integer, allocatable :: bpair(:,:)
integer i1,i2,j1,j2,nbond

! faster connectivity handling?
integer :: inat(nat),ic  ! number of attached atoms
integer :: ipair12(nat,1000) ! id of attached atoms, max 8

allocate(bond(nat,nat),cn(nat))
! bond matrix
! contains: 1=normal bond
!           2=H-bond
!           3=fragment bond
call bondmatrix(nat,iat,xyz,bond,cn)

! count internals

!counting bonds
k=0
kk=0
do i=1,nat-1
 ic=1
 do j=i+1,nat
   k=k+1
   if(i==j) cycle
   if(bond(i,j)==1) then
    kk=kk+1
  endif
 enddo
enddo
int_nb=kk


! pad bond matrix with H-bonds
! and interfragment bonds
! so that additional angles+torsions will be generated
call hbonds(nat,iat,xyz,bond,.false.)
allocate(frag_hval(frag_nh),frag_hcast(frag_nh,2))
frag_hval=0; frag_hcast=0
call hbonds(nat,iat,xyz,bond,.true.)
call connect_frag(nat,iat,xyz,bond)

! 1-2 pairs
ipair12=0
inat=0
do i=1,nat
 do j=1,nat
  if(i==j) cycle
!  if(bond(i,j)==1) then
  if(bond(i,j)>0) then ! this includes H- and F-bonds
    inat(i)=inat(i)+1
    ipair12(i,inat(i))=j
!    if(inat(i)>8) then ! this does not make sense with fragments
!      call error('too many connected atoms!')
!    endif
  endif
 enddo
enddo

if(debug) then
    print*, '1-2 pair list'
    do i=1,nat
      print*, 'elem:',esym(iat(i)),' #n12:',inat(i)
      print*, '->', ipair12(i,1:inat(i))
    enddo
endif

!# for a large system get only the bonds and do approx. vdw bonding pairs
! test
!--------EXPERIMENTAL START------------
int_nLJ=0
if (large) then
  print*, 'experimental LARGE mode'
  frag_nh=0; frag_nh=0; int_na=0;int_nt=0
  call read_intcoord('xopt.control',nat,xyz)
  kk=0
  do i=1,nat-1
  do j=i+1,nat
    if(i==j) cycle
    if(bond(i,j)==1) then
      kk=kk+1
      rab=dbond(xyz(1,i),xyz(1,j))
      int_bcast(kk,1)=i
      int_bcast(kk,2)=j
      int_bval(kk)=rab
    elseif(bond(i,j)==0) then
      rab=dbond(xyz(1,i),xyz(1,j))
      if(rab>int_lj_cut) cycle
      int_nLJ=int_nLJ+1
    endif
  enddo
  enddo
  allocate(int_LJcast(int_nLJ,2),int_LJval(int_nLJ))
  kk=0
  do i=1,nat-1
  do j=i+1,nat
    if(i==j) cycle
    if(bond(i,j)==0) then
      rab=dbond(xyz(1,i),xyz(1,j))
      if(rab>int_lj_cut) cycle
      kk=kk+1
      int_LJcast(kk,1)=i
      int_LJcast(kk,2)=j
      int_LJval(kk)=rab
    endif
  enddo
  enddo
  goto 666
endif
!--------EXPERIMENTAL END------------


nbond=int_nb+frag_nh+frag_nb
if(nbond>10000) call exclaim('over 10 000 bonds! ')
allocate(bpair(nbond,2))
bpair=0
kk=0
do i=1,nat-1
do j=i+1,nat
  if(i==j) cycle
  if(bond(i,j)>0) then
   kk=kk+1
   bpair(kk,1)=i
   bpair(kk,2)=j
 endif
enddo
enddo
if(kk/=nbond) call error('no no no')

! count angles
kk=0
do i=1,nat
 do j=1,inat(i)-1
   do k=j+1,inat(i)
    kk=kk+1
   enddo
 enddo
enddo
int_na=kk

! count torsions
ic=0
do i=1,nat
 do j=1,inat(i)
  jj=ipair12(i,j)
  if(i==jj) cycle
   do k=1,inat(jj)
    kk=ipair12(jj,k)
    if(kk==jj.or.i==kk) cycle
     do l=1,inat(kk)
      ll=ipair12(kk,l)
      if(ll==kk.or.ll==jj.or.ll<=i) cycle
      ic=ic+1
     enddo
   enddo
 enddo
enddo
int_nt=ic

! read custom internals and allocate arrays
call read_intcoord('xopt.control',nat,xyz)


! BONDS
k=0
kk=0
do i=1,nat-1
 do j=i+1,nat
   k=k+1
   if(i==j) cycle
   if(bond(i,j)==1) then
    kk=kk+1
    rab=dbond(xyz(1,i),xyz(1,j))
    int_bcast(kk,1)=i
    int_bcast(kk,2)=j
    int_bval(kk)=rab
  endif
 enddo
enddo

! ANGLES
ic=0
do i=1,nat
 do j=1,inat(i)
   jj=ipair12(i,j)
   do k=1,inat(jj)
      kk=ipair12(jj,k)
      if(i>=kk) cycle
      ic=ic+1
      call  angle(xyz,i,jj,kk,ang,agrad)
      int_acast(ic,1)=i
      int_acast(ic,2)=jj
      int_acast(ic,3)=kk
      int_aval(ic)=ang
   enddo
 enddo
enddo

! TORSIONS
ic=0
do i=1,nat
 do j=1,inat(i)
  jj=ipair12(i,j)
  if(i==jj) cycle
   do k=1,inat(jj)
    kk=ipair12(jj,k)
    if(kk==jj.or.i==kk) cycle
     do l=1,inat(kk)
      ll=ipair12(kk,l)
      if(ll==kk.or.ll==jj.or.ll<=i) cycle
      ic=ic+1
      call dihed(xyz,i,jj,kk,ll,ang,agrad)
      int_tcast(ic,1)=i
      int_tcast(ic,2)=jj
      int_tcast(ic,3)=kk
      int_tcast(ic,4)=ll
      int_tval(ic)=ang
     enddo
   enddo
 enddo
enddo


deallocate(bond)

  if(debug) then
    call print_primitives(.true.,.true.)
  else
  ! call print_primitives(.true.,.true.)
    call print_primitives(.false.,.true.)
  endif
  666 continue
  nints=int_nb+int_na+int_nt+frag_nh+frag_nb+int_nLJ
  call message('(see xopt.intcoord)')
  write(stdout,'(2x,a,1x,I7)') ' # primitives :', nints
  write(stdout,'(2x,a,1x,I7)') ' # bonds      :', int_nb
  write(stdout,'(2x,a,1x,I7)') ' # angles     :', int_na
  write(stdout,'(2x,a,1x,I7)') ' # torsions   :', int_nt
  write(stdout,'(2x,a,1x,I7)') ' # LJ bonds   :', int_nLJ
  write(stdout,'(2x,a)')       ' -------------'
  write(stdout,'(2x,a,1x,I7)') ' # H-bonds    :', frag_nh
  write(stdout,'(2x,a,1x,I7)') ' # frag bonds :', frag_nb
  write(stdout,'(a)') '' 
  write(stdout,'(2x,a)') 'Memory for internals:' 
  call p_memory(dble(nat)*3.0*dble(nints),'B-matrix')
  call p_memory(dble(nints),'H(diag)')
  call p_memory(dble(nints)*dble(nints),'H(int)')

end subroutine




subroutine print_primitives(do_screen,do_file)
use fiso, only: r8,stdout
use internals
use parm, only: iat
use constant, only: au2ang
implicit none
character(2) esym
integer it,i,j,k,l,io
real(r8) anggrad,di360
logical do_screen, do_file

open(newunit=io,file='xopt.intcoord',position='append')

call message('internal primitives:')
do it=1,int_nb
  i=int_bcast(it,1)
  j=int_bcast(it,2)
  if(do_screen) write(stdout,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'bond    ',&
        i,esym(iat(i)),j,esym(iat(j)),int_bval(it)*au2ang
  if(do_file) write(io,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'bond    ',&
        i,esym(iat(i)),j,esym(iat(j)),int_bval(it)*au2ang
enddo

do it=1,int_nLJ
  i=int_LJcast(it,1)
  j=int_LJcast(it,2)
  if(do_screen) write(stdout,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'LJ    ',&
        i,esym(iat(i)),j,esym(iat(j)),int_LJval(it)*au2ang
  if(do_file) write(io,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'LJ    ',&
        i,esym(iat(i)),j,esym(iat(j)),int_LJval(it)*au2ang
enddo


do it=1,int_na
 i=int_acast(it,1)
 j=int_acast(it,2)
 k=int_acast(it,3)
 if(do_screen) write(stdout,'(2x,a,1x,3(I5,''['',a2,'']''),1x,F8.4)')'angle   ',&
       i,esym(iat(i)),j,esym(iat(j)),k,esym(iat(k)),anggrad(int_aval(it))
 if(do_file) write(io,'(2x,a,1x,3(I5,''['',a2,'']''),1x,F8.4)')'angle   ',&
       i,esym(iat(i)),j,esym(iat(j)),k,esym(iat(k)),anggrad(int_aval(it))
enddo

do it=1,int_nt
 i=int_tcast(it,1)
 j=int_tcast(it,2)
 k=int_tcast(it,3)
 l=int_tcast(it,4)
 if(do_screen) write(stdout,'(2x,a,1x,4(I5,''['',a2,'']''),1x,F8.4)')'torsion ',&
       i,esym(iat(i)),j,esym(iat(j)),k,esym(iat(k)),l,esym(iat(l)),di360(int_tval(it))
 if(do_file) write(io,'(2x,a,1x,4(I5,''['',a2,'']''),1x,F8.4)')'torsion ',&
       i,esym(iat(i)),j,esym(iat(j)),k,esym(iat(k)),l,esym(iat(l)),di360(int_tval(it))
enddo

do it=1,frag_nh
 i=frag_hcast(it,1)
 j=frag_hcast(it,2)
 if(do_screen) write(stdout,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'H-bond  ',&
       i,esym(iat(i)),j,esym(iat(j)),frag_hval(it)*au2ang
 if(do_file) write(io,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'H-bond  ',&
       i,esym(iat(i)),j,esym(iat(j)),frag_hval(it)*au2ang
enddo

do it=1,frag_nb
 i=frag_bcast(it,1)
 j=frag_bcast(it,2)
 if(do_screen) write(stdout,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'F-bond  ',&
       i,esym(iat(i)),j,esym(iat(j)),frag_bval(it)*au2ang
 if(do_file) write(io,'(2x,a,1x,2(I5,''['',a2,'']''),1x,F8.4)')'F-bond  ',&
       i,esym(iat(i)),j,esym(iat(j)),frag_bval(it)*au2ang
enddo

end subroutine

! fixed 0/360 switches in torsion differences
subroutine torsionfix(t0,t1,dt)
use fiso, only: r8
implicit none
logical qI,qIV,shift
real(r8) t0,t1,dt
! Check for the case that the torsion goes from the I to the IV quadrant and adjust accordingly
! we check for "greater/less equal" since we might want to reach 0 as target value

!  IV to I
qI=.false.
qIV=.false.
shift=.false.
 if(t1.ge.0.and.t1.le.90) qI=.true.
 if(t0.ge.270.and.t0.le.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
  t1=t1+360_r8
  dt=t1-t0
  goto 999
 endif

 qI=.false.
 qIV=.false.

! I to VI
qI=.false.
qIV=.false.
shift=.false.
 if(t0.ge.0.and.t0.le.90) qI=.true.
 if(t1.ge.270.and.t1.le.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
  t1=t1-360_r8
  dt=t1-t0
  goto 999
 endif

999 continue

 if(.not.shift)  dt=t1-t0

return
end


subroutine read_intcoord(infile,nat,xyz0)
! read user defined bonds/angles/torsions
! NEEDS TESTING
! use constant, only:au2ang
use internals
use fiso, only: r8,stdout
implicit none
integer, intent(in) :: nat
real(r8), intent(in) :: xyz0(3,nat)! coordinate that will serve as reference
character(*),intent(in) :: infile
character(80) :: aa,bb
logical :: da,fstr
integer :: ii,jj,kk,ll
integer :: s2i,io
real(r8) :: dih,dgrad,ang,dang
real(r8) :: bdist,dbond
integer :: iint,n_ints,nb,na,nt


da=.false.
inquire(file=infile,exist=da)

if(da) then
  write(stdout,*) 'Reading custom internals from '//trim(infile)
  open(newunit=io,file=infile,status='old')

  ! determine size
  nb=0
  na=0
  nt=0
  do
   read(io,'(a)',end=125) aa
   if(index(aa,'$prim').ne.0) then
    do
      read(io,'(a)',end=126) aa
      if(index(aa,'#').ne.0) cycle
      if(fstr(aa,'bond')) nb=nb+1
      if(fstr(aa,'ang'))   na=na+1
      if(fstr(aa,'dihed').or.fstr(aa,'torsion'))  nt=nt+1
     enddo
   endif
  enddo
  125 continue
  126 continue

  n_ints=nb+na+nt

  write(stdout,'(2x,a,1x,I4)') ' # custom primitives :', n_ints
  write(stdout,'(2x,a,1x,I4)') ' # bonds             :', nb
  write(stdout,'(2x,a,1x,I4)') ' # angles            :', na
  write(stdout,'(2x,a,1x,I4)') ' # torsions          :', nt
  write(stdout,'(a)')''

  int_nb=int_nb+nb
  int_na=int_na+na
  int_nt=int_nt+nt
  allocate(int_bval(int_nb),int_bcast(int_nb,2))
  allocate(int_aval(int_na),int_acast(int_na,3))
  allocate(int_tval(int_nt),int_tcast(int_nt,4))
  int_bval=0d0;int_aval=0d0;int_tval=0d0
  int_bcast=0;int_acast=0;int_tcast=0

  rewind(io)
  iint=0
  nb=0
  na=0
  nt=0
  do
   read(io,'(a)',end=124) aa
   if(index(aa,'$prim').ne.0) then
   write(stdout,'(a)')''
   write(stdout,'(a)')' user definied internal coordinates: '
   write(stdout,'(a)')''
    do
      read(io,'(a)',end=123) aa
      if(index(aa,'#').ne.0) cycle

      if(fstr(aa,'bond')) then
       iint=iint+1
       nb=nb+1
       call charXsplit(aa,bb,2)
        ii=s2i(bb)
       call charXsplit(aa,bb,3)
        jj=s2i(bb)
        bdist=dbond(xyz0(1,ii),xyz0(1,jj))
        int_bval(nb)=bdist
        int_bcast(nb,1)=ii
        int_bcast(nb,2)=jj
        write(stdout,'(2x,a,2I4,2x,a,F7.2)') 'bond [A]      ', ii,jj,' | value: ',bdist
      endif

      if(fstr(aa,'ang')) then
       iint=iint+1
       na=na+1
       call charXsplit(aa,bb,2)
        ii=s2i(bb)
       call charXsplit(aa,bb,3)
        jj=s2i(bb)
       call charXsplit(aa,bb,4)
        kk=s2i(bb)
       call angle(xyz0,ii,jj,kk,ang,dang)
        int_aval(na)=dang
        int_acast(na,1)=ii
        int_acast(na,2)=jj
        int_acast(na,3)=kk
       write(stdout,'(2x,a,3I4,2x,a,F7.2)') 'angle [deg]   ', ii,jj,kk,' | value: ',dang
       if( abs(dang-180_r8)<1e-2_r8) call error('(quasi)-linear angle found!')
      endif

      if(fstr(aa,'dihed').or.fstr(aa,'torsion')) then
       iint=iint+1
       nt=nt+1
       call charXsplit(aa,bb,2)
        ii=s2i(bb)
       call charXsplit(aa,bb,3)
        jj=s2i(bb)
       call charXsplit(aa,bb,4)
        kk=s2i(bb)
       call charXsplit(aa,bb,5)
        ll=s2i(bb)
       call dihed(xyz0,ii,jj,kk,ll,dih,dgrad)
       int_tval(nt)=dgrad
       int_tcast(nt,1)=ii
       int_tcast(nt,2)=jj
       int_tcast(nt,3)=kk
       int_tcast(nt,4)=ll
       write(stdout,'(2x,a,4I4,2x,a,F7.2)') 'torsion [deg] ',ii,jj,kk,ll,' | value: ',dgrad
      endif

      if(index(aa,'$').ne.0) then
       backspace(io)
       exit
      endif
    enddo
   endif

  enddo
  123 continue
  124 continue
  close(io)

  ! should not happen
  if(iint/=n_ints) stop 'I/O error in custom internals'

else
  allocate(int_bval(int_nb),int_bcast(int_nb,2))
  allocate(int_aval(int_na),int_acast(int_na,3))
  allocate(int_tval(int_nt),int_tcast(int_nt,4))
  int_bval=0d0;int_aval=0d0;int_tval=0d0
  int_bcast=0;int_acast=0;int_tcast=0
endif


end subroutine

! print all values about thresh for an input vector
subroutine print_threshold(n,vec,thresh)
use fiso, only: r8
implicit none
real(r8) vec(*), thresh
integer i,n

do i=1,n
  if(vec(i)>=thresh) print '(1x,F8.2)',vec(i)
enddo
end subroutine




subroutine get_Bmatrix(nat,xyz)
use internals
use fiso, only: r8
use parm, only: B
use logic, only: debug
implicit none
integer i,j,it,k,l,idx
integer ii,jj,kk,ll,io
real(r8) e(3),va(3),vb(3),vc(3),vd(3)
real(r8) eba(3),ebc(3),rba,rbc,rcd
real(r8) v12(3),v23(3),v34(3),v43(3),v32(3),v21(3),crosa(3),crosb(3)
real(r8) si,lab,x1,x2,x3,x4,co,co2,si2
integer nat
real(r8) xyz(3,nat)
real(r8) eig(nints)


b=0_r8
! bonds
idx=0
do it=1,int_nb
 idx=idx+1
 i=int_bcast(it,1)
 j=int_bcast(it,2)
 ii=3*i-2
 jj=3*j-2
 e(1:3)=xyz(1:3,i)-xyz(1:3,j)
 call veclen(e,lab)
 e=e/lab

 b(idx,ii+0)=e(1)
 b(idx,ii+1)=e(2)
 b(idx,ii+2)=e(3)
 b(idx,jj+0)=-e(1)
 b(idx,jj+1)=-e(2)
 b(idx,jj+2)=-e(3)
enddo

do it=1,int_nLJ
 idx=idx+1
 i=int_LJcast(it,1)
 j=int_LJcast(it,2)
 ii=3*i-2
 jj=3*j-2
 e(1:3)=xyz(1:3,i)-xyz(1:3,j)
 call veclen(e,lab)
 e=e/lab

 b(idx,ii+0)=e(1)
 b(idx,ii+1)=e(2)
 b(idx,ii+2)=e(3)
 b(idx,jj+0)=-e(1)
 b(idx,jj+1)=-e(2)
 b(idx,jj+2)=-e(3)
enddo

if(debug) print*,'bonds done'

! fragment coordinates
do it=1,frag_nh
 idx=idx+1
 i=frag_hcast(it,1)
 j=frag_hcast(it,2)
 ii=3*i-2
 jj=3*j-2
 e(1:3)=xyz(1:3,i)-xyz(1:3,j)
 call veclen(e,lab)
 e=e/lab

 b(idx,ii+0)=e(1)
 b(idx,ii+1)=e(2)
 b(idx,ii+2)=e(3)
 b(idx,jj+0)=-e(1)
 b(idx,jj+1)=-e(2)
 b(idx,jj+2)=-e(3)
enddo

do it=1,frag_nb
 idx=idx+1
 i=frag_bcast(it,1)
 j=frag_bcast(it,2)
 ii=3*i-2
 jj=3*j-2
 e(1:3)=xyz(1:3,i)-xyz(1:3,j)
 call veclen(e,lab)
 e=e/lab

 b(idx,ii+0)=e(1)
 b(idx,ii+1)=e(2)
 b(idx,ii+2)=e(3)
 b(idx,jj+0)=-e(1)
 b(idx,jj+1)=-e(2)
 b(idx,jj+2)=-e(3)
enddo



! angles
do it=1,int_na
 idx=idx+1
 i=int_acast(it,1)
 j=int_acast(it,2)
 k=int_acast(it,3)
 ii=i*3-2
 jj=j*3-2
 kk=k*3-2
  
  va(1:3)=xyz(1:3,i)
  vb(1:3)=xyz(1:3,j)
  vc(1:3)=xyz(1:3,k)
  call evec(va,vb,eba)
  call evec(vc,vb,ebc)
  call veclen2(va,vb,rba)
  call veclen2(vc,vb,rbc)
  co=dot_product(eba,ebc)
  si=sqrt(1.0_r8-co**2)
  e(1:3)=(co*eba(1:3)-ebc(1:3))/(rba*si)
  b(idx,ii+0)=e(1)
  b(idx,ii+1)=e(2)
  b(idx,ii+2)=e(3)
  e(1:3)=((rba-rbc*co)*eba(1:3)+(rbc-rba*co)*ebc(1:3))/(rba*rbc*si)
  b(idx,jj+0)=e(1)
  b(idx,jj+1)=e(2)
  b(idx,jj+2)=e(3)
  e(1:3)=(co*ebc(1:3)-eba(1:3))/(rbc*si)
  b(idx,kk+0)=e(1)
  b(idx,kk+1)=e(2)
  b(idx,kk+2)=e(3)

enddo

if(debug) print*,'angles done'
! torsions
do it=1,int_nt
 idx=idx+1
 i=int_tcast(it,1)
 j=int_tcast(it,2)
 k=int_tcast(it,3)
 l=int_tcast(it,4)
 ii=i*3-2
 jj=j*3-2
 kk=k*3-2
 ll=l*3-2

 va(1:3)=xyz(1:3,i)
 vb(1:3)=xyz(1:3,j)
 vc(1:3)=xyz(1:3,k)
 vd(1:3)=xyz(1:3,l)
 call evec(va,vb,v12)
 call evec(vb,vc,v23)
 call evec(vc,vd,v34)
 call evec(vd,vc,v43)
 call evec(vc,vb,v32)
 call evec(vb,va,v21)
 co=DOT_PRODUCT(v21,v23)
 co2=DOT_PRODUCT(v32,v34)
 si=sqrt(1.0_r8-co**2)
 si2=sqrt(1.0_r8-co2**2)
 call cross_prod(crosa,v12,v23)
 call cross_prod(crosb,v43,v32)
 call veclen2(va,vb,lab)
 call veclen2(vb,vc,rbc)
 call veclen2(vc,vd,rcd)
 x1=(rbc-lab*co)/(rbc*lab*si**2)
 x2=co2/(rbc*si2**2)
 x3=(rbc-rcd*co2)/(rbc*rcd*si2**2)
 x4=co/(rbc*si**2)

 e(1:3)=-crosa(1:3)/(lab*si**2)
  b(idx,ii+0)=e(1)
  b(idx,ii+1)=e(2)
  b(idx,ii+2)=e(3)
 e(1:3)=x1*crosa(1:3)+x2*crosb(1:3)
  b(idx,jj+0)=e(1)
  b(idx,jj+1)=e(2)
  b(idx,jj+2)=e(3)
 e(1:3)=x3*crosb(1:3)+x4*crosa(1:3)
  b(idx,kk+0)=e(1)
  b(idx,kk+1)=e(2)
  b(idx,kk+2)=e(3)
 e(1:3)=-crosb(1:3)/(rcd*si2**2)
  b(idx,ll+0)=e(1)
  b(idx,ll+1)=e(2)
  b(idx,ll+2)=e(3)

enddo
if(debug) print*,'dihedrals done'
if(idx/=nints) call error('B-matrix construction error')

! call screenMat(nints,nat*3,b,epsilon(1d0),.true.,'B-matrix')
if(debug) then
 open(newunit=io,file='xopt.Bmatrix')
 call printMat(io,nints,nat*3,b,'B-Matrix')
endif

return
end subroutine


subroutine get_G_matrix()
! calculates
! G = B*B^t and its pseudo-inverse
! returns Ginv
use fiso, only: r8,stdout
use internals, only: nints,Ginv,Gmat,Umat
use parm, only: nat3,B,nvar
use logic, only: debug, int_deloc
implicit none
integer i,nzero,k,first,j
real(r8) gvec(nints),gmat2(nints,nints)


write(stdout,'(a)', advance="no") 'making G matrix ...'
call matmult('n','t',nints,nints,nat3,b,b,gmat)
write(stdout,'(a)') ' done'
if(debug) call printMat(6,nints,nints,gmat,'G=B*B^t')

! diag
gmat2=gmat
call DiagSM(nints,gmat2,gvec)
! non-zero elements are our degress of freedom for deloc ints
 nzero=0
 do i=1,nints
  if(gvec(i)>1d-10) then
   first=i
   exit
  endif
 enddo
 do i=1,nints
  if(gvec(i)>1d-10) nzero=nzero+1
 enddo
 if(int_deloc) then
    nvar=nzero
    allocate(Umat(nints,nzero))
    k=0
    do i=1,nints
      do j=1,nzero
        umat(i,j)=gmat2(i,first+j)
      enddo
    enddo
    write(stdout,'(2x,a,I5)') 'degrees of freedom (deloc internals) : ',nzero
    !make new H(int,nonred)
    ! H(int,nonred)=U^T*H(int,diag)*U
 else
  nvar=nints
  write(stdout,'(a)', advance="no") 'making pseudo-inverse of G  ...'
  call mat_pinv(gmat,nints,nints,ginv,1d-10)
  write(*,'(a)') ' done'
  if(debug) call printMat(6,nints,nints,ginv,'G^-')
  write(stdout,'(2x,a,I5)') 'degrees of freedom (red internals) : ',nvar
 endif




end subroutine

! UNTESTED
subroutine Hcart2Hint
! H(int)=G^- * B * H(cart) * B^t * G^-
! H(int)=inhess(nints,nints)
use fiso, only:r8,stdout
use internals, only: nints,inthess,Ginv
use parm, only: nat3,B,chess
real(r8) :: scr1(nat3,nints),scr2(nat3,nints),scr3(nints,nints)

! Bt x G
call matmult('t','n',nat3,nints,nints,b,ginv,scr1)
! H(cart) x BtG
call matmult('n','n',nat3,nints,nat3,chess,scr1,scr2)
! B x HBtG
call matmult('n','n',nints,nints,nat3,b,scr2,scr3)
! G x BHBtG
call matmult('n','n',nints,nints,nints,ginv,scr3,inhess)
call printMat(stdout,nints,nints,inhess,"H(int)")
end subroutine


subroutine Hdiag2Cart()
use internals, only: hdiag,nints
use parm, only: B,nat3,chess
use fiso, only: r8
implicit none
real(r8) scr(nints,nat3)
integer i,j

scr=0.0
do j=1,nat3
 do i=1,nints
 scr(i,j)=hdiag(i)*B(i,j)
 enddo
enddo
call matmult('t','n',nat3,nat3,nints,b,scr,chess)

end subroutine

subroutine Hint2Cart(nvar,hint,chess)
!Hxyz=Bt*H*B
! note that the hint2xyz has essentially B^t !
! This is an unfortunate unconsistency...
! HERE B IS (nvar,nat3)
use fiso, only: r8
use parm, only: B,nat3
implicit none
integer :: nvar
real(r8) :: hint(nvar*(nvar+1)/2)
real(r8), intent(out) :: chess(nat3,nat3)
real(r8) :: h(nvar,nvar)
real(r8) :: scr(nvar,nat3)

call packM(nvar,h,hint,'unpack')
call matmult('n','n',nvar,nat3,nvar,hint,b,scr)
call matmult('t','n',nat3,nat3,nvar,b,scr,chess)
end subroutine

subroutine grad2int(gcart,gint)
! G^- x B x g(cart)
use fiso, only: r8
use internals, only: nints,Ginv
use parm, only: nat,B,chess,nat3
implicit none
real(r8),intent(in) :: gcart(3,nat)
real(r8),intent(out) :: gint(nints)
real(r8) :: gc(nat3),scr(nints)
integer i,j,k

k=0
do i=1,nat
  do j=1,3
    k=k+1
    gc(k)=gcart(j,i)
  enddo
enddo

call DGEMV('n',nints,nat3,1d0,B,nints,gc,1,1d0,scr,1)
call DGEMV('n',nints,nints,1d0,Ginv,nints,scr,1,1d0,gint,1)

end subroutine

subroutine proj_nonred(gint,hint)
! project g(int) and H(int) onto nonredundant space
! with projector P=GG- dim: (nints,nints)
! g(nonred)=P*g(int)
! H(nonred)=P*H(int)*P
! returns: new hint(packed) and gint
use internals,only: nints,Ginv,Gmat
use fiso, only: r8,stdout
use logic, only: debug
implicit none
real(r8),intent(inout) :: gint(nints)
real(r8) proj(nints,nints),gint2(nints)
real(r8) :: hint(nints*(nints+1)/2)
real(r8) :: h(nints,nints),scr(nints,nints)

!check programmer error

if(debug) write(stdout,'(a)', advance="no") 'making projector P(int) ...'
call matmult('n','n',nints,nints,nints,gmat,ginv,proj)
if(debug) write(stdout,'(a)') ' done'
if(debug) call printMat(6,nints,nints,proj,"P'(int)")

! project gradient g(int)
call dgemv('n',nints,nints,1d0,proj,nints,gint,1,1d0,gint2,1)
gint=gint2

! project H(int), unpack first
call packM(nints,h,hint,'unpack')
call matmult('n','n',nints,nints,nints,h,proj,scr)
call matmult('n','n',nints,nints,nints,proj,scr,h)
call packM(nints,h,hint,'pack')
end subroutine


subroutine get_Amatrix(Amat)
! A=B^t*G^-
! matrix for displacement backtransformation:
! dipl(cart)=A*displ(int)
use fiso, only: r8
use parm, only: B,nat3
use internals, only: Ginv,nints
implicit none
real(r8), intent(out) :: Amat(nat3,nints)

call matmult('t','n',nat3,nints,nints,b,ginv,amat)
end subroutine


subroutine dint2cart(dint,xyznew)
! dipl(cart)=A*displ(int)
! updates xyz coordinates
use internals, only: nints
use fiso, only: r8
use parm, only: nat,xyz0,nat3
implicit none
real(r8), intent(in) :: dint(nints)
real(r8), intent(out) :: xyznew(3,nat) 
real(r8) :: dcart(nat3),Amat(nat3,nints)
integer ii,i,j

call get_Amatrix(Amat)
call dgemv('n',nat3,nints,1d0,amat,nat3,dint,1,0d0,dcart,1)

ii=0
xyznew=0_r8
do i=1,nat
  do j=1,3
    ii=ii+1
    xyznew(j,i)=xyz0(j,i)+dcart(ii)
  enddo
enddo

end subroutine

subroutine make_intcoords()
use fiso, only: r8,stdout
use logic, only: debug
use internals
implicit none

! FOR NOW ONLY DELOCALS
call make_primitives()
allocate(Ginv(nints,nints),Gmat(nints,nints))
call get_G_matrix()


end subroutine


subroutine make_primitives()
use internals
use parm
use timing
implicit none
    call status1(timer)
    call message_head('* internal primitives *')
    call get_primitives(nat,iat,xyz)
    call status2(timer)

    call status1(timer)
    call message_head('* B-Matrix *')
    allocate(b(nints,nat3),hdiag(nints))
    call get_Bmatrix(nat,xyz)
    call status2(timer)
end subroutine
