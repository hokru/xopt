C hydrogen bond correction a la DH+ (M. Korth)
      subroutine hbsimple(echo,grad,n,at,xyz,hbscale,energy,g)
      implicit none 
      integer n,at(n)
      real*8 xyz(3,n),hbscale,g(3,n),energy
c gradient?
      logical grad
c printout?
      logical echo
  
      integer i,j,k,m,A,B,H,nhb,hbpar,typa,typb
      real*8 ,allocatable ::hbs(:,:)
      real*8 e,er,el,step,rab
c scale factors for N,O,F,P,S,Cl
c currently all set to the same value hbscale
      real*8 scalehb(6)       
c RAB^2 cutoff in Bohr^2
      real*8 :: thr
c short cutoff dftd3 radii in bohr scaled by shortcutscale
      real*8 r0ab(6,6),shortcut,cab
      real*8 :: shortcutscale=1.00

      scalehb(1:6)=hbscale

      thr=250.d0

      r0ab(1,1)=4.95580619390096    
      r0ab(2,1)=4.69521322049756    
      r0ab(2,2)=4.68973278168166    
      r0ab(3,1)=4.51361038303964    
      r0ab(3,2)=4.44293461883933    
      r0ab(3,3)=4.34561357746278    
      r0ab(4,1)=5.75705004325974    
      r0ab(4,2)=5.42861568913110    
      r0ab(4,3)=5.22773805284199    
      r0ab(4,4)=6.61725321392116    
      r0ab(5,1)=5.49040984445458    
      r0ab(5,2)=5.44335574472965    
      r0ab(5,3)=5.16462109512098    
      r0ab(5,4)=6.27011084754718    
      r0ab(5,5)=6.25631558644032    
      r0ab(6,1)=5.28216218036620    
      r0ab(6,2)=5.18862031695641    
      r0ab(6,3)=5.07977251245376    
      r0ab(6,4)=6.03124949908998    
      r0ab(6,5)=5.95698288505716    
      r0ab(6,6)=5.86684309279986    
      do i=1,6
         do j=1,i
            r0ab(j,i)=r0ab(i,j)
         enddo
      enddo
      r0ab=r0ab*shortcutscale

c count HBs for allocate
      nhb=0
      do i=1,n-1
         if(at(i).eq.7.or.at(i).eq.8.or.at(i).eq.9
     .  .or.at(i).eq.15.or.at(i).eq.16.or.at(i).eq.17)then
         do j=i+1,n
            if(at(j).eq.7.or.at(j).eq.8.or.at(j).eq.9
     .     .or.at(j).eq.15.or.at(j).eq.16.or.at(j).eq.17)then
            rab=((xyz(1,i)-xyz(1,j))**2
     .          +(xyz(2,i)-xyz(2,j))**2
     .          +(xyz(3,i)-xyz(3,j))**2)
            if(rab.lt.thr) then 
            do k=1,n
            if(at(k).eq.1)nhb=nhb+1
            enddo
            endif
            endif
         enddo
         endif
      enddo

      if(echo)then
         write(*,*) 
         write(*,*) 'HBSIMPLE correction switched on'
         write(*,*) '# H-bonds considered',nhb
      endif   

      allocate(hbs(3,nhb))

c assign HBs: 1=A, 2=B, 3=H
      nhb=0
      do i=1,n-1
         if(at(i).eq.7.or.at(i).eq.8.or.at(i).eq.9
     .  .or.at(i).eq.15.or.at(i).eq.16.or.at(i).eq.17)then
         do j=i+1,n
            if(at(j).eq.7.or.at(j).eq.8.or.at(j).eq.9
     .     .or.at(j).eq.15.or.at(j).eq.16.or.at(j).eq.17)then
            rab=((xyz(1,i)-xyz(1,j))**2
     .          +(xyz(2,i)-xyz(2,j))**2
     .          +(xyz(3,i)-xyz(3,j))**2)
            if(rab.lt.thr) then 
            do k=1,n
            if(at(k).eq.1)then
               nhb=nhb+1
               hbs(1,nhb)=i
               hbs(2,nhb)=j
               hbs(3,nhb)=k
            endif
            enddo
            endif
            endif
         enddo
         endif
      enddo

c compute energy
      if(echo.and.n.lt.100)write(*,*)'true H-bonds A-H-B, in kcal/mol' 

      energy=0
      do m=1,nhb
         A=hbs(1,m)
         B=hbs(2,m)
         H=hbs(3,m)
         typa=hbpar(at(A))
         typb=hbpar(at(B))
         shortcut=r0ab(typa,typb)
         cab=0.50d0*(scalehb(typa)+scalehb(typb))
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,e)
         if(echo.and.e.lt.-0.002.and.n.lt.100) 
     .   write(*,'(3I6,F12.5)')A,H,B,e*627.51
         energy=energy+e   
      enddo

      if(grad)then
c compute gradient (seems to be faster numerically than with
c longish maple code)
      g=0
      step=1.d-5
      do m=1,nhb
         A=hbs(1,m)
         B=hbs(2,m)
         H=hbs(3,m)
         typa=hbpar(at(A))
         typb=hbpar(at(B))
         cab=0.50d0*(scalehb(typa)+scalehb(typb))
         shortcut=r0ab(hbpar(at(A)),hbpar(at(B)))
         do j=1,3
         xyz(j,A)=xyz(j,A)+step
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,er)
         xyz(j,A)=xyz(j,A)-step*2.
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,el)
         xyz(j,A)=xyz(j,A)+step
         g(j,A)=g(j,A)+(er-el)/(2.0d0*step)
         enddo
         do j=1,3
         xyz(j,B)=xyz(j,B)+step
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,er)
         xyz(j,B)=xyz(j,B)-step*2.
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,el)
         xyz(j,B)=xyz(j,B)+step
         g(j,B)=g(j,B)+(er-el)/(2.0d0*step)
         enddo
         do j=1,3
         xyz(j,H)=xyz(j,H)+step
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,er)
         xyz(j,H)=xyz(j,H)-step*2.
         call eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,el)
         xyz(j,H)=xyz(j,H)+step
         g(j,H)=g(j,H)+(er-el)/(2.0d0*step)
         enddo
      enddo
      endif

      deallocate(hbs)
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine eabh(n,A,B,H,at,xyz,scalehb,shortcut,cab,energy)
      implicit none 
      integer A,B,H,n,at(n)
      real*8 xyz(3,n),scalehb(*),energy,shortcut,cab
  
      integer i,j,k,m,nhb
      real*8 rab,ang,xy,cosabh,d2ij,d2ik,d2jk,term
      real*8 aterm,dampm,damps,dampl,xm,ym,zm,rhm,rab2
      REAL*8 :: longcut=7.50388 !  5 angstroem
      REAL*8 :: alp=20

c AB distance
         rab2=   ((xyz(1,A)-xyz(1,B))**2
     .           +(xyz(2,A)-xyz(2,B))**2
     .           +(xyz(3,A)-xyz(3,B))**2)
         rab=sqrt(rab2)
c cos angle A-H-B
         D2IJ = (XYZ(1,A)-XYZ(1,H))**2+
     1          (XYZ(2,A)-XYZ(2,H))**2+
     2          (XYZ(3,A)-XYZ(3,H))**2
         D2JK = (XYZ(1,H)-XYZ(1,B))**2+
     1          (XYZ(2,H)-XYZ(2,B))**2+
     2          (XYZ(3,H)-XYZ(3,B))**2
         D2IK = (XYZ(1,A)-XYZ(1,B))**2+
     1          (XYZ(2,A)-XYZ(2,B))**2+
     2          (XYZ(3,A)-XYZ(3,B))**2
         XY = SQRT(D2IJ*D2JK+1.d-14)
         cosabh = 0.5D0 * (D2IJ+D2JK-D2IK) / XY 
c angle term
         aterm = 0.125d0*(cosabh-1.0d0)**4
c H-(AB) midpoint distance
         xm=0.5*(XYZ(1,A)+XYZ(1,B))
         ym=0.5*(XYZ(2,A)+XYZ(2,B))
         zm=0.5*(XYZ(3,A)+XYZ(3,B))
         rhm=sqrt((xyz(1,H)-xm)**2
     .           +(xyz(2,H)-ym)**2
     .           +(xyz(3,H)-zm)**2)
c H out of linear damping
         dampm=1.d0-1.d0/(1.d0+exp(-alp*(rhm/rab-1.d0)))
c short damping
         damps=1.d0/(1.d0+exp(-alp*(rab/shortcut-1.d0)))
c long  damping
         dampl=1.d0-1.d0/(1.d0+exp(-alp*(rab/longcut-1.d0)))
c all for this A-H-B
         energy=-dampl*damps*dampm*aterm*cab/rab2**2    
c        write(*,'(3I4,5F12.5)')
c    .   A,B,H,57*acos(cosabh),damps,dampl,dampm,term*627.5

         end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function hbpar(elem)
      integer elem
      if(elem.eq.7) hbpar=1
      if(elem.eq.8) hbpar=2
      if(elem.eq.9) hbpar=3
      if(elem.eq.15)hbpar=4
      if(elem.eq.16)hbpar=5
      if(elem.eq.17)hbpar=6
      end
