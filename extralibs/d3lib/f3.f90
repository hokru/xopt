!> @brief Hydrogen bond correction for semi-empirical methods
!>
!> @details Korth's Hydrogen correction subroutines for semi-empirical methods
!>          eg. PM6-DH+
!>
!> @author Jimmy Kromann
!> - Nov 2013
!>

module hp_hbcorr

    ! Parameters
    double precision, parameter :: hartree2kcal=627.509541d0 ! korth
    double precision, parameter :: bohr2angstroem=0.52917726d0 ! TODO Check with newest constants
    double precision, parameter :: pi=3.1415926535d0
    integer, parameter :: mxatm=2000

    ! Default F3 parameters
    double precision, parameter :: expo_x=2.d0
    double precision, parameter :: shortcut=4.535d0 !  2.4 angstroem
    double precision, parameter :: longcut=13.228d0 !  7.0 angstroem
    double precision, parameter :: distcut=19.842d0 ! 10.5 angstroem (zero longcut) 
    double precision, parameter :: covcut=2.268d0   !  1.2 angstroem (hb mid point)
    double precision, parameter :: distcut2=2.646d0 !  1.4 angstroem (zero covcut)

    ! Covalent radii in angstroems from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197
    ! values for metals decreased by 10 %, all multiplied by 4/3 both done as for D3
    double precision :: covrad(94)
    data covrad /&
    &0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,&
    &1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54,&
    &1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,&
    &1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39,&
    &1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26,&
    &1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57,&
    &1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,&
    &1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,&
    &1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58,&
    &1.52, 1.53, 1.54, 1.55 /

    ! Logical administration
    logical :: lmpchbond, virgin
    data virgin /.true./

    ! pm6-d3h+ parameters
    !data scale_n / -0.15d0/
    !data scale_o / -0.11d0/

    ! pm6-dh+ parameters
    ! data scale_nsp3 / -0.16d0/
    ! data scale_nsp2 / -0.16d0/
    ! data scale_osp3 / -0.12d0/
    ! data scale_osp2 / -0.12d0/

contains

!
!> @brief Distance between atom a and b
!>
!> @author Jimmy Kromann
!> - April 2013
!>
!> @param a
!> @param b
!>
double precision function distance(a, b, geo)
    implicit none
    double precision :: geo(3, mxatm)
    integer :: a, b
    distance = sqrt((geo(1,a) - geo(1,b))**2 + (geo(2,a) - geo(2,b))**2 + (geo(3,a) - geo(3,b))**2)
    distance = distance
end function distance
!
!> @brief
!>
!> @author Jimmy Kromann
!> - April 2013
!>
!>
double precision function angle(a, b, c, geo)
    implicit none
    integer, intent(in) :: a,b,c
    double precision :: geo(3,mxatm)
    double precision :: ab,cd,ac,de
    ab=(geo(1,a)-geo(1,b))**2+(geo(2,a)-geo(2,b))**2+(geo(3,a)-geo(3,b))**2
    cd=(geo(1,b)-geo(1,c))**2+(geo(2,b)-geo(2,c))**2+(geo(3,b)-geo(3,c))**2
    ac=(geo(1,a)-geo(1,c))**2+(geo(2,a)-geo(2,c))**2+(geo(3,a)-geo(3,c))**2
    de=sqrt(ab*cd)
    if(de.lt.1.d-20)then
    angle=0.d0
    return
    endif
    angle=0.5d0*(ab+cd-ac)/de
    angle=min(1.0D0,angle)
    angle=dmax1(-1.0D0,angle)
    angle=acos(angle)
end function angle
!
!> @brief
!>
!> @author Jimmy Kromann
!> - April 2013
!>
!>
double precision function torsion(i, j, k, l, geo)
    implicit none
    integer, intent(in) :: i, j, k, l
    double precision :: geo(3, mxatm)
    double precision :: xi1, xj1, xl1, yi1, yj1, yl1, zi1, zj1, zl1, dist, cosa
    double precision :: ddd, yxdist, xi2, xl2, yi2, yl2, costh, sinth, cosph, sinph, yj2, yi3,yl3
    xi1 = geo(1,i) - geo(1,k)
    xj1 = geo(1,j) - geo(1,k)
    xl1 = geo(1,l) - geo(1,k)
    yi1 = geo(2,i) - geo(2,k)
    yj1 = geo(2,j) - geo(2,k)
    yl1 = geo(2,l) - geo(2,k)
    zi1 = geo(3,i) - geo(3,k)
    zj1 = geo(3,j) - geo(3,k)
    zl1 = geo(3,l) - geo(3,k)
    ! ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
    dist = sqrt(xj1*xj1 + yj1*yj1 + zj1*zj1)
    cosa = zj1/dist
    cosa = min(1.0D0,cosa)
    cosa = dmax1(-1.0D0,cosa)
    ddd = 1.0D0 - cosa**2
    if (ddd <= 0.0D0) go to 10
    yxdist = dist*sqrt(ddd)
    if (yxdist > 1.0D-6) go to 20
    10 continue
    xi2 = xi1
    xl2 = xl1
    yi2 = yi1
    yl2 = yl1
    costh = cosa
    sinth = 0.D0
    go to 30
    20 continue
    cosph = yj1/yxdist
    sinph = xj1/yxdist
    xi2 = xi1*cosph - yi1*sinph
    xl2 = xl1*cosph - yl1*sinph
    yi2 = xi1*sinph + yi1*cosph
    yj2 = xj1*sinph + yj1*cosph
    yl2 = xl1*sinph + yl1*cosph
    ! ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
    costh = cosa
    sinth = yj2/dist
    30 continue
    yi3 = yi2*costh - zi1*sinth
    yl3 = yl2*costh - zl1*sinth
    torsion=dangle(xl2,yl3,xi2,yi3)
    ! 6.2831853  IS 2 * 3.1415926535 = 180 DEGREE
    if (torsion < 0.) torsion = pi*2.d0 + torsion
    if (torsion >= 6.28318530717959D0) torsion = 0.D0
    return
end function torsion
!
!> @brief
!>
!> @author Jimmy Kromann
!> - April 2013
!>
!>
double precision function dangle(a1,a2,b1,b2)
    implicit none
    double precision, INTENT(INOUT) :: a1,a2,b1,b2
    double precision :: zero, anorm, bnorm, sinth, costh, rcos
    zero = 1.0D-6
    if (abs(a1)>=zero .or. abs(a2)>=zero) then
    if (abs(b1)>=zero .or. abs(b2)>=zero) then
        anorm = 1.0D0/sqrt(a1**2 + a2**2)
        bnorm = 1.0D0/sqrt(b1**2 + b2**2)
        a1 = a1*anorm
        a2 = a2*anorm
        b1 = b1*bnorm
        b2 = b2*bnorm
        sinth = a1*b2 - a2*b1
        costh = a1*b1 + a2*b2
        costh = min(1.0D0,costh)
        costh = dmax1(-1.0D0,costh)
        rcos = acos(costh)
        if (abs(rcos) >= 4.0D-5) then
        if (sinth > 0.D0) rcos = 6.28318530717959D0 - rcos
        rcos = -rcos
        dangle = rcos
        return
        endif
    endif
    endif
    rcos = 0.0D0
    dangle = rcos
    return
end function dangle
!
!> @brief
!>
!> @author Jimmy Kromann
!> - April 2013
!>
!>
double precision function bonding(x, y, labels)
    implicit none
    integer, intent(in) :: x, y
    integer :: labels(mxatm)
    bonding = covrad(labels(x)) + covrad(labels(y))
    return
end function bonding


function cross_product(a,b)
    implicit none
    double precision, INTENT(IN) :: a(3),b(3)
    double precision cross_product(3)
    cross_product(1)=a(2)*b(3)-a(3)*b(2)
    cross_product(2)=a(3)*b(1)-a(1)*b(3)
    cross_product(3)=a(1)*b(2)-a(2)*b(1)
end function cross_product


double precision function vangle(a,b)
    implicit none
    double precision, INTENT(IN) :: a(3),b(3)
    vangle=acos(dot_product(a,b)/absolute(a)/absolute(b))
end function vangle


double precision function absolute(a)
    implicit none
    double precision, INTENT(IN) :: a(3)
    absolute=sqrt(a(1)**2+a(2)**2+a(3)**2)
end function absolute


!**********************************************************************


!> @brief
!>
!> @author Jimmy Kromann
!> - Nov 2013
!>
!> @param natoms, number of atoms
!> @param geo, geometry of the system
!> @param labels, the labels of the atoms
!>
subroutine hbond_energy(natoms, geo, labels, param, hb_energy)
    implicit none

    ! in
    integer :: mode, natoms
    double precision :: geo(3, mxatm)
    integer :: labels(mxatm)
    double precision :: param(2)

    ! out
    double precision :: hb_energy, gradient
    double precision :: hb_gradient(3, mxatm)

    mode = 1
    call hbond_f3(mode, natoms, geo, labels, param, hb_energy, hb_gradient)

end subroutine
!> @brief
!>
!> @author Jimmy Kromann
!> - Nov 2013
!>
!> @param natoms, number of atoms
!> @param geo, geometry of the system
!> @param labels, the labels of the atoms
!>
subroutine hbond_gradient(natoms, geo, labels, param, hb_energy, gradient)
    implicit none

    ! in
    integer :: mode, natoms
    double precision :: geo(3, mxatm)
    integer :: labels(mxatm)
    double precision :: param(2)

    ! out
    double precision :: hb_energy, gradient(3, mxatm)
    double precision :: hb_gradient(3, mxatm)


    !
    logical :: numgrad
    integer :: i, j

    numgrad = .False.

    mode = 2

    call hbond_f3(mode, natoms, geo, labels, param, hb_energy, hb_gradient)

    do j=1,natoms
        do i=1,3
            gradient(i,j) = gradient(i,j) + hb_gradient(i,j)
        enddo
    enddo

    ! Numerical gradient


end subroutine
!> @brief
!>
!> @author Jimmy Kromann
!> - Nov 2013
!>
!> @param natoms, number of atoms
!> @param geo, geometry of the system
!> @param labels, the labels of the atoms
!>
subroutine hbond_f3(mode, natoms, geo, labels, param, hb_energy, hb_gradient)
    implicit none

    ! in
    integer :: mode, natoms
    double precision :: geo(3, mxatm)
    integer :: labels(mxatm)
    double precision :: param(2)

    ! out
    double precision :: hb_energy, hb_gradient(3, mxatm)

    ! Administration
    logical :: debug, verbose

    ! HB system
    integer, allocatable :: hbs(:,:), bondlist(:,:)

    ! possible pairs
    integer :: npairs, heteroatoms
    double precision :: ab_dist, ha_dist, hb_dist

    ! specify pairs
    logical :: hbsa, hbsb
    integer, allocatable :: nbondsa(:), nbondsb(:)
    integer :: nbondsc
    double precision :: old_dist, xa_dist, xb_dist, xh_dist, xc_dist

    ! Angles
    double precision :: theta
    double precision :: cos_theta, cos_phi_a, cos_phi_b, cos_psi_a, cos_psi_b
    double precision :: cos_phi_a2, cos_phi_b2
    double precision :: cos_psi_a2, cos_psi_b2
    logical :: angle_check,torsion_check_a,torsion_check_b !MKnew2

    ! dE Angles
    double precision :: de_cos_theta
    double precision :: de_cos_phi_a, de_cos_phi_a2
    double precision :: de_cos_phi_b, de_cos_phi_b2

    double precision :: de_cos_psi_a, de_cos_psi_a2
    double precision :: de_cos_psi_b, de_cos_psi_b2

    ! Energy
    double precision :: scale_a, scale_b, scale_c
    double precision :: scale_n, scale_o
    double precision :: xy_dist, zh_dist, zh_inf
    double precision :: damping, damping_h, damping_l, damping_s
    double precision :: hb_correction

    ! Gradient
    logical :: lgradient
    double precision :: de_dr, de_dtheta, de_dphi_a, de_dphi_b, de_dpsi_a, de_dpsi_b, de_dpsi_a2, de_dpsi_b2
    double precision :: de_dh, de_dh_part, de_dl, de_dl_part, de_ds, de_ds_part, temp

    integer :: i, j, k, l

    if(virgin) then
        ! JCK Lazy way of updating the constants
        ! but makes it easiere to read what is happening
        covrad = 4.0d0/3.d0 * covrad / bohr2angstroem
        virgin = .False.
    endif

    ! Energy or Gradient calculation?
    select case(mode)
        case(1)
            lgradient = .False.
        case(2)
            lgradient = .True.
    end select

    ! Parameter
    scale_n = param(1)
    scale_o = param(2)

    ! Debug
    verbose = .False.
    debug = .False.

    heteroatoms=0
    do i=1, natoms
        if(labels(i)==7.or.labels(i)==8) heteroatoms = heteroatoms + 1
    enddo

    ! To avoid bugs
    if(heteroatoms==0) return

    allocate(hbs(heteroatoms*250, 10))
    allocate(bondlist(mxatm, 3))
    allocate(nbondsa(heteroatoms*250))
    allocate(nbondsb(heteroatoms*250))

!    print *, 'allo', heteroatoms*250

    npairs = 0
    hbs(:,:) = -1

    ! HBS Map
    ! 0:
    ! 1: Atom A
    ! 2:
    ! 3:
    ! 4:
    ! 5: Atom B
    ! 6:
    ! 7:
    ! 8:
    ! 9: Hydrogen atom
    !10: Status:
    !        -1 = nothing checked
    !      -666 = bad

    ! Find possible pairs
    do i=1,natoms

        if(labels(i).ne.7.and.labels(i).ne.8) cycle

        do j=i+1,natoms
            if(labels(j).ne.7.and.labels(j).ne.8) cycle

            ab_dist = distance(i, j, geo)

            if(ab_dist.gt.distcut) cycle

            do k=1,natoms
                if(labels(k)==1) then
                    ha_dist=distance(k, i, geo)
                    hb_dist=distance(k, j, geo)

                    if(ha_dist.gt.distcut2.and.hb_dist.gt.distcut2)cycle

                    npairs=npairs+1
                    hbs(npairs,1)=i
                    hbs(npairs,5)=j
                    hbs(npairs,9)=k

                endif
            enddo
        enddo
    enddo

!    print *, 'npairs', npairs

    ! Specify pairs
    do i=1,npairs
        nbondsa(i) = 0
        nbondsb(i) = 0
        bondlist(:,:) = 0

        ! Check how many atoms are bonded with
        ! the current a and b atoms for pair = 1
        do j=1,natoms
            xa_dist = distance(j, hbs(i,1), geo)
            xb_dist = distance(j, hbs(i,5), geo)

            if(xa_dist.lt.bonding(j, hbs(i,1), labels).and.hbs(i,1).ne.j)then
                nbondsa(i) = nbondsa(i) + 1
                bondlist(nbondsa(i),1) = j
            endif

            if(xb_dist.lt.bonding(j, hbs(i,5), labels).and.hbs(i,5).ne.j)then
                nbondsb(i) = nbondsb(i) + 1
                bondlist(nbondsb(i),2) = j
            endif

            if(nbondsa(i).gt.4.or.nbondsb(i).gt.4)then
                write(*,*) "Too many bonds for some atoms!"
                write(*,*) i, hbs(i,1), hbs(i,5), nbondsa(i), nbondsb(i)
                hbs(i,10) = -666
            endif

        enddo

        ! Check if the acceptor is attached to the donor.
        ! 1-4 and 1-3 interaction
        do j=1,nbondsa(i)
            if(bondlist(j,1).eq.hbs(i,5)) hbs(i,10) = -666
            do k=1,nbondsb(i)
                if(bondlist(j,1).eq.bondlist(k,2)) hbs(i,10) = -666
            enddo
        enddo

        ! Cycle if the specific HB is bad.
        if(hbs(i,10).eq.-666) cycle

        ! hbs1 part
        hbsa=.False.

        if(nbondsa(i).eq.3.or.nbondsa(i).eq.4)then

            old_dist=-1

            do k=1,nbondsa(i)
                xh_dist = distance(bondlist(k,1), hbs(i,9), geo)
                if(xh_dist.gt.old_dist)then
                    old_dist = xh_dist
                    hbs(i,2) = bondlist(k,1)
                endif
            enddo

            old_dist=-1

            do k=1,nbondsa(i)
                xh_dist = distance(bondlist(k,1), hbs(i,9), geo)
                if(xh_dist.gt.old_dist.and.bondlist(k,1).ne.hbs(i,2))then
                    old_dist = xh_dist
                    hbs(i,3) = bondlist(k,1)
                endif
            enddo

            old_dist=-1

            do k=1,nbondsa(i)
                xh_dist = distance(bondlist(k,1), hbs(i,9), geo)
                if(xh_dist.gt.old_dist.and.bondlist(k,1).ne.hbs(i,2).and.bondlist(k,1).ne.hbs(i,3))then
                    old_dist = xh_dist
                    hbs(i,4) = bondlist(k,1)
                endif
            enddo ! possible forth atom should be h

        elseif(nbondsa(i).eq.2)then

            old_dist=-1

            do k=1,nbondsa(i)
                xh_dist = distance(bondlist(k,1), hbs(i,9), geo)
                if(xh_dist.gt.old_dist)then
                    old_dist = xh_dist
                    hbs(i,2) = bondlist(k,1)
                endif
            enddo

            do k=1,nbondsa(i)
                if(bondlist(k,1).ne.hbs(i,2))then
                    hbs(i,3) = bondlist(k,1)
                endif
            enddo

            if(distance(hbs(i,1), hbs(i,9), geo).lt.bonding(hbs(i,1), hbs(i,9), labels))then !MK needed
                hbs(i,4)=hbs(i,9)
            else
                hbs(i,4)=hbs(i,1)
            endif

        elseif(nbondsa(i).eq.1)then

            hbs(i,2) = bondlist(1,1)

            ! need bonds on first bonded atom here
            nbondsc = 0

            do k=1,natoms
                xc_dist = distance(k, hbs(i,2), geo)
                if(xc_dist.lt.bonding(k, hbs(i,2), labels).and.hbs(i,2).ne.k)then
                    nbondsc = nbondsc + 1
                    bondlist(nbondsc,3) = k
                endif
            enddo

            old_dist=-1

            do k=1,nbondsc
                xh_dist = distance(bondlist(k,3), hbs(i,9), geo)
                if(xh_dist.gt.old_dist)then
                    old_dist = xh_dist
                    hbs(i,3) = bondlist(k,3)
                endif
            enddo

            if(distance(hbs(i,1), hbs(i,9), geo).lt.bonding(hbs(i,1), hbs(i,9), labels))then !MK needed
                hbs(i,4) = hbs(i,9)
            else
                hbs(i,4) = hbs(i,1)
            endif

        elseif(nbondsa(i).eq.0)then

            if(distance(hbs(i,1), hbs(i,9), geo).lt.bonding(hbs(i,1), hbs(i,9), labels))then !MK needed
                hbs(i,2) = hbs(i,9)
                hbs(i,3) = hbs(i,9)
                hbs(i,4) = hbs(i,9)
            else
                hbs(i,2) = hbs(i,1)
                hbs(i,3) = hbs(i,1)
                hbs(i,4) = hbs(i,1)
            endif

        else
            hbsa=.True.
        endif


        ! hbs2 part
        hbsb=.False.

        if(nbondsb(i).eq.3.or.nbondsb(i).eq.4)then

            old_dist=-1

            do k=1,nbondsb(i)
                xh_dist = distance(bondlist(k,2),hbs(i,9), geo)
                if(xh_dist.gt.old_dist)then
                    old_dist = xh_dist
                    hbs(i,6) = bondlist(k,2)
                endif
            enddo

            old_dist=-1

            do k=1,nbondsb(i)
                xh_dist = distance(bondlist(k,2),hbs(i,9), geo)
                if(xh_dist.gt.old_dist.and.bondlist(k,2).ne.hbs(i,6))then
                    old_dist = xh_dist
                    hbs(i,7) = bondlist(k,2)
                endif
            enddo

            old_dist=-1

            do k=1,nbondsb(i)
                xh_dist = distance(bondlist(k,2),hbs(i,9), geo)
                if(xh_dist.gt.old_dist.and.bondlist(k,2).ne.hbs(i,6).and.bondlist(k,2).ne.hbs(i,7))then
                old_dist = xh_dist
                hbs(i,8) = bondlist(k,2)
                endif
            enddo ! possible forth atom should be h

        elseif(nbondsb(i).eq.2)then

            old_dist=-1

            do k=1,nbondsb(i)
                xh_dist = distance(bondlist(k,2), hbs(i,9), geo)
                if(xh_dist.gt.old_dist)then
                    old_dist = xh_dist
                    hbs(i,6) = bondlist(k,2)
                endif
            enddo

            do k=1,nbondsb(i)
                if(bondlist(k,2).ne.hbs(i,6))then
                    hbs(i,7)=bondlist(k,2)
                endif
            enddo

            if(distance(hbs(i,5), hbs(i,9), geo).lt.bonding(hbs(i,5), hbs(i,9), labels))then !MK needed
                hbs(i,8)=hbs(i,9)
            else
                hbs(i,8)=hbs(i,5)
            endif

        elseif(nbondsb(i).eq.1)then

            hbs(i,6) = bondlist(1,2)

            ! need bonds on first bonded atom here
            nbondsc = 0

            do k=1,natoms
                xc_dist = distance(k, hbs(i,6), geo)
                if(xc_dist.lt.bonding(k, hbs(i,6), labels).and.hbs(i,6).ne.k)then
                    nbondsc = nbondsc + 1
                    bondlist(nbondsc,3) = k
                endif
            enddo

            old_dist=-1

            do k=1,nbondsc
                xh_dist = distance(bondlist(k,3), hbs(i,9), geo)
                if(xh_dist.gt.old_dist)then
                    old_dist = xh_dist
                    hbs(i,7) = bondlist(k,3)
                endif
            enddo

            if(distance(hbs(i,5), hbs(i,9), geo).lt.bonding(hbs(i,5), hbs(i,9), labels))then !MK needed
                hbs(i,8)=hbs(i,9)
            else
                hbs(i,8)=hbs(i,5)
            endif

        elseif(nbondsb(i).eq.0)then

            if(distance(hbs(i,5), hbs(i,9), geo).lt.bonding(hbs(i,5), hbs(i,9), labels))then !MK needed
                hbs(i,6) = hbs(i,9)
                hbs(i,7) = hbs(i,9)
                hbs(i,8) = hbs(i,9)
            else
                hbs(i,6) = hbs(i,5)
                hbs(i,7) = hbs(i,5)
                hbs(i,8) = hbs(i,5)
            endif

        else
            hbsb=.True.
        endif


        ! write pair data
        if(.False.)then
            write(*,*)"Possible pair found:"
            write(*,*)i,hbs(i,1),hbs(i,5),hbs(i,9)
            write(*,*)"atom nbonds a",hbs(i,1),nbondsa(i)
            do j=1,nbondsa(i)
                write(*,*)bondlist(j,1)
            enddo
            write(*,*)"atom nbonds b",hbs(i,5),nbondsb(i)
            do j=1,nbondsb(i)
                write(*,*)bondlist(j,2)
            enddo
            write(*,*)
        endif
    enddo


    ! check for problems
    if(hbsb.or.hbsa)then
        if(.true.)then
            write(*,*)
            write(*,*)"Found NO or unknown h-bond constellation." 
            write(*,*)"In the later case, please contact me to fix this ..."
        endif
    endif

    ! Write h-bond data
    if(verbose.and.debug.and..False.)then
        write(*,*)"Possible h-bonds found:"
        do i=1,npairs
            if(hbs(i,10).ne.-666) write(*,'(10I3)') i, hbs(i,4), hbs(i,3), hbs(i,2), &
                & hbs(i,1), hbs(i,9), hbs(i,5), hbs(i,6), hbs(i,7), hbs(i,8)
        enddo
        write(*,*)
        write(*,*)"Hbond corrections:"
    elseif(verbose)then
        write(*,*)"Hbonds:"
    endif

    ! Calculate internal coordinates

    ! theta = angle(a, b, hydrogen)
    ! phi_a = angle()
    ! phi_b = angle()
    ! psi_a = torsion()
    ! psi_b = torsion()

    ! HBS Map
    ! 0:
    ! 1: Atom i
    ! 2:
    ! 3:
    ! 4:
    ! 5: Atom j
    ! 6:
    ! 7:
    ! 8:
    ! 9: Hydrogen atom
    !10: Status:
    !        -1 = nothing checked
    !      -666 = bad

    hb_energy = 0.0d0
    hb_gradient = 0.0d0
    hb_correction = 0.0d0

    do i=1,npairs

        ! Check if valid HB
        if(hbs(i,10)==-666) cycle

        theta = angle(hbs(i,1), hbs(i,9), hbs(i,5), geo)
        cos_theta = -cos(theta) ! theta = PI - Theta

        ! If HB is at a 90 degree angle, skip,
        ! because that is not a physical HB.
        if(cos_theta.le.0) cycle

        ! Comment about ^2
        de_cos_theta = 2.0d0*cos_theta*sin(theta)

        !
        ! Atom A
        !
        call hb_angles(hbs(i,1), hbs(i,2), hbs(i,3), hbs(i,4), hbs(i,9), &
            & geo, labels, nbondsa(i), &
            & cos_phi_a, cos_psi_a, &
            & cos_phi_a2, cos_psi_a2, &
            & de_cos_phi_a, de_cos_psi_a, &
            & de_cos_phi_a2, de_cos_psi_a2, &
            & angle_check,torsion_check_a) !MKnew2

        ! Something failed when calculating the angles
        if(angle_check) cycle

        !
        ! Atom B
        !
        call hb_angles(hbs(i,5), hbs(i,6), hbs(i,7), hbs(i,8), hbs(i,9), &
            & geo, labels, nbondsb(i), &
            & cos_phi_b, cos_psi_b, &
            & cos_phi_b2, cos_psi_b2, &
            & de_cos_phi_b, de_cos_psi_b, &
            & de_cos_phi_b2, de_cos_psi_b2, &
            & angle_check,torsion_check_b) !MKnew2

        ! Something failed when calculating the angles
        if(angle_check) cycle

        !
        ! Parameters
        !
        select case (labels(hbs(i,1)))
            case (7)
                scale_a = scale_n
            case (8)
                scale_a = scale_o
        end select

        select case(labels(hbs(i,5)))
            case(7)
                scale_b = scale_n
            case(8)
                scale_b = scale_o
        end select

        scale_c = (scale_a + scale_b)/2.0d0

        !
        ! Distance
        !
        ha_dist = distance(hbs(i,9), hbs(i,1), geo) ! zh_inf == 1
        hb_dist = distance(hbs(i,9), hbs(i,5), geo) ! zh_inf != 1

        zh_dist = min(ha_dist, hb_dist)
        zh_inf = 1
        if(zh_dist.eq.hb_dist) zh_inf=2

        !
        ! Energy
        !
        damping_h = 1.d0 - 1.d0/(1.d0 + exp(-60.d0*(zh_dist/covcut-1.d0)))

        ! y-x damping
        xy_dist = distance(hbs(i,1), hbs(i,5), geo)

        ! short range
        damping_s = 1.d0/(1.d0 + exp(-100.d0*(xy_dist/shortcut-1.d0)))

        ! long range
        damping_l = 1.d0 - 1.d0/(1.d0+exp(-10.d0*(xy_dist/longcut-1.d0)))

        damping = damping_h*damping_s*damping_l


        hb_correction = scale_c/xy_dist**2.0d0*cos_theta**2*&
                        &cos_phi_a**2*cos_psi_a**2*&
                        &cos_phi_b**2*cos_psi_b**2*&
                        &damping

        hb_energy = hb_energy + hb_correction

        !
        ! Gradient
        !
        if(lgradient.and..False.) then

            ! Numerical gradient


        elseif(lgradient) then

            ! r
            de_dr = -2.0d0*scale_c/xy_dist**3*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            & *cos_phi_b**2*cos_psi_b**2*damping

            ! theta
            de_dtheta = scale_c/xy_dist**2*de_cos_theta*cos_phi_a**2*cos_psi_a**2 &
            & *cos_phi_b**2*cos_psi_b**2*damping

            ! phi
            de_dphi_a = scale_c/xy_dist**2*cos_theta**2*de_cos_phi_a*cos_psi_a**2 &
            & *cos_phi_b**2*cos_psi_b**2*damping

            de_dphi_b = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            & *de_cos_phi_b*cos_psi_b**2*damping

            ! psi
            de_dpsi_a = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*de_cos_psi_a &
            & *cos_phi_b**2*cos_psi_b**2*damping

            de_dpsi_b = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            & *cos_phi_b**2*de_cos_psi_b*damping

            !de_dpsi_a2 = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*de_cos_psi_a2 &
            !& *cos_phi_b**2*cos_psi_b**2*damping

            !de_dpsi_a2 = de_dpsi_a2 + scale_c/xy_dist**2*cos_theta**2*de_cos_phi_a2*cos_psi_a**2 &
            !& *cos_phi_b**2*cos_psi_b**2*damping

            !de_dpsi_b2 = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a2**2 &
            !& *cos_phi_b**2*de_cos_psi_b2*damping

            !de_dpsi_b2 = de_dpsi_b2 + scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            !& *de_cos_phi_b2*cos_psi_b**2*damping

            temp = exp(-60.d0*(zh_dist/covcut-1.d0))

            de_dh_part = 60.d0/covcut*temp*(1.0d0/(1.0d0+temp)**2)

            de_dh = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            & *cos_phi_b**2*cos_psi_b**2*de_dh_part*damping_s*damping_l

            temp = exp(-100.d0*(xy_dist/shortcut-1.d0))

            de_ds_part = -100.d0/shortcut*temp*(1.0d0/(1.0d0+temp)**2)

            de_ds = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            & *cos_phi_b**2*cos_psi_b**2*de_ds_part*damping_h*damping_l

            temp = exp(-10.d0*(xy_dist/longcut-1.d0))

            de_dl_part = 10.d0/longcut*temp*(1.0d0/(1.0d0+temp)**2)

            de_dl = scale_c/xy_dist**2*cos_theta**2*cos_phi_a**2*cos_psi_a**2 &
            & *cos_phi_b**2*cos_psi_b**2*de_dl_part*damping_h*damping_s

            de_dr = de_dr + de_ds + de_dl

            if(.False.)then
                write(*,*) "de_dr", de_dr
                write(*,*) "de_dtheta", de_dtheta
                write(*,*) "de_dphi_a", de_dphi_a
                write(*,*) "de_dphi_b", de_dphi_b
                write(*,*) "de_dpsi_a", de_dpsi_a
                write(*,*) "de_dpsi_b", de_dpsi_b
                !write(*,*) "de_dpsi_a2", de_dpsi_a2
                !write(*,*) "de_dpsi_b2", de_dpsi_b2
                write(*,*) "de_dh", de_dh
            endif

            ! Move from internal to cartesian coordinates
            ! and fill out hb_gradient

            call int2xyz(2, de_dr, hbs(i,1), hbs(i,5),0 ,0, geo, hb_gradient)

            call int2xyz(3, de_dtheta, hbs(i,1), hbs(i,9), hbs(i,5), 0, geo, hb_gradient)

            call int2xyz(3, de_dphi_a, hbs(i,2), hbs(i,1), hbs(i,9), 0, geo, hb_gradient)
            call int2xyz(3, de_dphi_b, hbs(i,6), hbs(i,5), hbs(i,9), 0, geo, hb_gradient)

            !MKnew2 start
            if(torsion_check_a)then
             call int2xyz(4, -de_dpsi_a, hbs(i,3), hbs(i,4), hbs(i,1), hbs(i,9), geo, hb_gradient) 
            else
             call int2xyz(4, -de_dpsi_a, hbs(i,3), hbs(i,2), hbs(i,1), hbs(i,9), geo, hb_gradient)
            endif
            !call int2xyz(4, de_dpsi_a2, hbs(i,3), hbs(i,2), hbs(i,1), hbs(i,4), geo, hb_gradient)
            if(torsion_check_b)then
             call int2xyz(4, -de_dpsi_b, hbs(i,7), hbs(i,8), hbs(i,5), hbs(i,9), geo, hb_gradient) 
            else
             call int2xyz(4, -de_dpsi_b, hbs(i,7), hbs(i,6), hbs(i,5), hbs(i,9), geo, hb_gradient)
            endif
            !call int2xyz(4, de_dpsi_b2, hbs(i,7), hbs(i,6), hbs(i,5), hbs(i,8), geo, hb_gradient)
            !MKnew2 stop

            if(zh_inf.eq.1)then
                call int2xyz(2, de_dh, hbs(i,9), hbs(i,1), 0, 0, geo, hb_gradient)
            else
                call int2xyz(2, de_dh, hbs(i,9), hbs(i,5), 0, 0, geo, hb_gradient)
            endif

        endif ! gradient

    enddo ! npairs

end subroutine



!> @brief Calculat Phi and Psi given a geometry
!>
!> @author Jimmy Kromann
!> - Nov 2013
!>
!> @param a, acceptor/donor atom index
!> @param b, 
!> @param c, 
!> @param h, Hydrogen index
!> @param geo, geometry of the system
!> @param labels, the labels of the atoms
!> @param nbonds,
!> @return cos_phi
!> @return cos_psi
subroutine hb_angles(a, b, c, d, h, &
                   & geo, labels, nbonds, &
                   & cos_phi, cos_psi, &
                   & cos_phi2, cos_psi2, &
                   & de_cos_phi, de_cos_psi, &
                   & de_cos_phi2, de_cos_psi2, &
                   & check,psi_check_set2) !MKnew2
    implicit none

    ! in
    integer :: a, b, c, d, h

    double precision :: geo(3, mxatm)
    integer :: labels(mxatm)
    integer :: nbonds

    ! out
    double precision :: cos_phi, cos_psi
    double precision :: cos_phi2, cos_psi2
    double precision :: de_cos_phi, de_cos_psi
    double precision :: de_cos_phi2, de_cos_psi2
    logical :: check

    !
    double precision :: phi, psi
    double precision :: phi_shift, phi_shift2

    double precision :: psi_check, psi_check_bac, psi_check_factor
    double precision :: psi_shift, psi_correct, psi_value, psi_value_2

    logical :: psi_check_set, psi_check_set2, planar !MK

    check = .False.

    !MK start
    ! R3X...H
    ! a X
    ! b R3
    ! c R2
    ! d R1
    ! h H
    !MK stop

    psi_check_set = .False.
    psi_check_set2 = .False.
    psi_check = 0.0d0

    ! Check for shifting angles
    if(labels(a)==8)then
        if(nbonds==1)then
            phi_shift = PI
            phi_shift2 = PI/(180.0/120.0)

            psi_shift = 0.0
            psi_check_set2 = .true. ! H...O=C-R1
        else
            phi_shift = PI/(180.0/109.48)
            phi_shift2 = phi_shift

            psi_shift = PI/(180.0/54.74)
        endif
    elseif(labels(a)==7)then
        if(nbonds==2)then
            phi_shift = PI/(180.0/120.0)
            phi_shift2 = phi_shift

            psi_shift = 0.0
        else
            !MK start
            planar=.false.
            ! check for double bonds
            if(a.ne.b.and.distance(a,b,geo).lt.bonding(a,b,labels)*3.0d0/4.0d0*0.97)planar=.true. !MKnew5
            if(a.ne.c.and.distance(a,c,geo).lt.bonding(a,c,labels)*3.0d0/4.0d0*0.97)planar=.true. !MKnew5
            if(a.ne.d.and.distance(a,d,geo).lt.bonding(a,d,labels)*3.0d0/4.0d0*0.97)planar=.true. !MKnew5
            if(.not.planar)then !NR3t
             phi_shift = PI/(180.0/109.48)
             phi_shift2 = phi_shift

             psi_shift = PI/(180.0/54.74)
             psi_check_set = .true.
            else !NR3p
             phi_shift = PI/(180.0/120.0)
             phi_shift2 = phi_shift

             psi_shift = PI/(180.0/90.0)
            endif
            !MK stop
        endif
    endif

    !MK now without interpolation between NR3t and NR3p 
    ! Torsion check is set to hard to false, because of alot of problems
    ! psi_check_set = .False.

    ! extrapolation between tetragonal and planar NR3 group
    if(psi_check_set)then

        !JCK psi_check = torsion(hbs(i,3), hbs(i,2), hbs(i,1), hbs(i,4), geo) !torsion2
        psi_check = torsion(c, b, a, d, geo) !torsion2

        !MKtest if(psi_check.le.-PI) psi_check = psi_check + 2.0 * PI
        !MKtest if(psi_check.gt.PI)  psi_check = psi_check - 2.0 * PI

        if(psi_check.lt.0)then
            psi_check =-PI - psi_check
        else
            psi_check = PI - psi_check
        endif

        psi_check_bac = psi_check ! save sign

        if(psi_check.lt.0) psi_check = -1.0 * psi_check

        if(psi_check.ne.0) then ! for grad
            psi_check_factor = psi_check_bac/psi_check
        else
            psi_check_factor = 0.0
        endif

        psi_check = psi_check*180.0/PI
        !MK start
        !psi_shift = psi_shift + PI/(180.0/((54.74-psi_check)/54.74*35.26))
        !phi_shift  = phi_shift - PI/(180.0/((54.74-psi_check)/54.74*19.48))
        !MK phi_shift  = phi_shift + PI/(180.0/((54.74-psi_check)/54.74*11.52))
        !MK stop
        phi_shift2 = phi_shift

        psi_check = psi_check_bac ! restore sign
    endif

    ! Calculate Phi
    phi = angle(b, a, h, geo)
    cos_phi = cos(phi_shift - phi)
    de_cos_phi = 2.0d0*cos_phi*sin(phi_shift-phi)

    !JCK Depends on torsion
    de_cos_phi2 = 0.0
    !MK start
    !if(psi_check_set)then
    !    de_cos_phi2 = -1.0/(54.74*19.48)*2*cos_phi*sin(phi_shift-phi)*psi_check_factor
    !    !de_cos_phi2 = 1.0/(54.74*11.52)*2*cos_phi*sin(phi_shift-phi)*psi_check_factor
    !endif
    !MK stop
    !JCK

    cos_phi2 = cos(phi_shift2-phi)

    if(cos_phi2.gt.cos_phi)then
        cos_phi = cos_phi2
        de_cos_phi = 2.0d0*cos_phi*sin(phi_shift2-phi) 
        !MK start
        !if(psi_check_set)then
            !de_cos_phi2 = -1.0/(54.74*19.48)*2*cos_phi*sin(phi_shift2-phi)*psi_check_factor
            !MK de_cos_phi2 = 1.0/(54.74*11.52)*2*cos_phi*sin(phi_shift2-phi)*psi_check_factor
        !endif
        !MK stop
    endif

    ! If the cos(phi) is less than 0, cycle
    if(cos_phi.le.0) check = .True.


    ! torsion angle
    ! psi_correct = torsion(hbs(i,3), hbs(i,2), hbs(i,1), hbs(i,9), geo)
    psi_correct = torsion(c, b, a, h, geo)
    if(psi_check_set2)psi_correct=torsion(c,d,a,h,geo) !MKnew 

    !MKtest if(.not.psi_check_set2)then !MKnew4
    !MKtest  if(psi_correct.le.-PI) psi_correct = psi_correct + 2.0 * PI
    !MKtest  if(psi_correct.gt.PI)  psi_correct = psi_correct - 2.0 * PI
    !MKtest endif !MKnew4

    if((.not.psi_check_set2).or.(abs(psi_correct*180.0/PI).gt.90))then
        if(psi_correct.lt.0)then
            psi_correct = -PI - psi_correct
        else
            psi_correct =  PI - psi_correct
        endif
    endif

    ! correction of NR3 torsion angle for through-bond case
    if(psi_check.lt.0.0d0)then ! negative torsion angle occupied by -NR3 r3

        psi_value = psi_shift - psi_correct
        !MKtest if(psi_value.le.-PI) psi_value = psi_value + 2.0 * PI
        !MKtest if(psi_value.gt.PI)  psi_value = psi_value - 2.0 * PI

        cos_psi = cos(psi_value)

        de_cos_psi = 2*cos_psi*sin(psi_value)
        !MK de_cos_psi2 = 1.0/(54.75*35.26)*2*cos_psi*sin(psi_value)*psi_check_factor
        de_cos_psi2 = 0.0d0

    elseif(psi_check.gt.0.0d0)then ! positive torsion angle occupied by -NR3 r3

        psi_value = -psi_shift - psi_correct

        !MKtest if(psi_value.le.-PI) psi_value = psi_value + 2.0 * PI
        !MKtest if(psi_value.gt.PI)  psi_value = psi_value - 2.0 * PI

        cos_psi = cos(psi_value)

        de_cos_psi = 2*cos_psi*sin(psi_value)
        !MK de_cos_psi2 = 1.0/(54.75*35.26)*2*cos_psi*sin(psi_value)*psi_check_factor
        de_cos_psi2 = 0.0d0

    else ! planar -NR3 or general case

        psi_value   =  psi_shift - psi_correct
        psi_value_2 = -psi_shift - psi_correct

        !MKtest if(psi_value.le.-PI) psi_value = psi_value + 2.0 * PI
        !MKtest if(psi_value.gt.PI)  psi_value = psi_value - 2.0 * PI

        !MKtest if(psi_value_2.le.-PI) psi_value_2 = psi_value_2 + 2.0 * PI
        !MKtest if(psi_value_2.gt.PI)  psi_value_2 = psi_value_2 - 2.0 * PI

        cos_psi = cos(psi_value)
        de_cos_psi = 2*cos_psi*sin(psi_value)

        de_cos_psi2 = 0.0d0
        !MK start
        !if(psi_check_set)then
        !    de_cos_psi2 = 1.0/(54.75*35.26)*2*cos_psi*sin(psi_value)*psi_check_factor
        !endif
        !MK stop
        cos_psi2 = cos(psi_value_2)

        if(cos_psi2.gt.cos_psi)then
            cos_psi = cos_psi2
            de_cos_psi = 2*cos_psi*sin(psi_value_2)
            !MK start
            !if(psi_check_set)then
            !    de_cos_psi = 1.0/(54.75*35.26)*2*cos_psi*sin(psi_value_2)*psi_check_factor
            !endif
            !MK stop
        endif

    endif

    !MK corrects for H...C=0 instead of C=0...H case
    if(distance(h, a, geo).gt.distance(h, b, geo).and.psi_check_set2)then
        cos_psi=0.d0
        de_cos_psi=0.0d0 !MKnew
        de_cos_psi2=0.0d0
    endif

    !MK corrects for cases with not enough atoms for defining psi
    !MKnewStart
    if(b.eq.a.or.b.eq.h)then 
        cos_phi=1.d0
        de_cos_phi=0.0d0
        de_cos_psi2=0.0d0
    endif
    if(c.eq.a.or.c.eq.h) then 
        cos_psi=1.d0
        de_cos_psi=0.d0
        de_cos_psi2=0.0d0
    endif
    if(d.eq.h)then
        cos_phi=1.d0
        de_cos_phi=0.0d0
        cos_psi=1.d0
        de_cos_psi=0.0d0
        de_cos_psi2=0.0d0
    endif
    !MKnewStop

    ! If cos(psi) is less than zero, skip
    if(cos_psi.le.0) check = .True.


end subroutine


!> @brief Transform gradient from internal coordinates to cartesian
!>
!> @author Jimmy Kromann
!> - Nov 2013
!>
!> @param a, acceptor/donor atom index
subroutine int2xyz(i, g, a, b, c, d, geo, gradient)
    implicit none

    integer, intent(in) :: i, a, b, c, d
    double precision, intent(in) :: g
    double precision, intent(in) :: geo(3,mxatm)
    double precision :: gradient(3, mxatm)

    integer :: n

    double precision :: ad1, ad2, d1d2
    double precision :: cd1(3), cd2(3)
    double precision :: cd1d2(3), cd2d3(3)
    double precision :: tmp1, tmp2, tmp3, tmp4
    double precision :: ctmp1(3), ctmp2(3)
    double precision :: grad(3,4)

    if(i.eq.2)then ! distance
        if(a.eq.b) return
        ad1=absolute(geo(:,a)-geo(:,b))
        cd1=geo(:,a)-geo(:,b)

        gradient(:,a) = gradient(:,a) + g/ad1*cd1(:)
        gradient(:,b) = gradient(:,b) - g/ad1*cd1(:)
    elseif(i.eq.3)then ! angle
        if(a.eq.b.or.a.eq.c.or.b.eq.c)return
        ad1=absolute(geo(:,a)-geo(:,b))
        ad2=absolute(geo(:,c)-geo(:,b))
        d1d2=dot_product(geo(:,a)-geo(:,b),geo(:,c)-geo(:,b))
        cd1=geo(:,a)-geo(:,b)
        cd2=geo(:,c)-geo(:,b)
        ctmp1(:)=g/(1.0-d1d2**2/ad1**2/ad2**2)**0.5*(-1.0d0)*(cd2(:)/ad1/ad2-cd1(:)*d1d2/ad1**3/ad2)
        ctmp2(:)=g/(1.0-d1d2**2/ad1**2/ad2**2)**0.5*(-1.0d0)*(cd1(:)/ad1/ad2-cd2(:)*d1d2/ad1/ad2**3)

        gradient(:,a) = gradient(:,a) + ctmp1(:)
        gradient(:,b) = gradient(:,b) - ctmp1(:)-ctmp2(:)
        gradient(:,c) = gradient(:,c) + ctmp2(:)
    elseif(i.eq.4)then ! torsion
        if(a.eq.b.or.a.eq.c.or.a.eq.d.or.b.eq.c.or.b.eq.d.or.c.eq.d) return
        tmp1=vangle(geo(:,b)-geo(:,a),geo(:,c)-geo(:,b))
        tmp2=vangle(geo(:,c)-geo(:,b),geo(:,d)-geo(:,c))

        if(tmp1.eq.0.or.tmp1.eq.PI) return
        if(tmp2.eq.0.or.tmp2.eq.PI) return

        ad2=absolute(geo(:,c)-geo(:,b))
        cd1d2=cross_product(geo(:,b)-geo(:,a),geo(:,c)-geo(:,b))
        cd2d3=cross_product(geo(:,c)-geo(:,b),geo(:,d)-geo(:,c))

        tmp1=ad2*dot_product(geo(:,b)-geo(:,a),cd2d3)
        tmp2=dot_product(cd1d2,cd2d3)
        tmp3=tmp1**2+tmp2**2

        ! DEBUG: u=b-a | v=c-b | w=d-c | t=a-c | s=b-d | cuv=c12 | cvw=c23
        grad(:,1)=(cross_product(geo(:,c)-geo(:,b),cd2d3(:))*tmp1-cd2d3(:)*ad2*tmp2)*g/tmp3
        grad(:,2)=(-dot_product(geo(:,b)-geo(:,a),cd2d3(:))/ad2*tmp2*(geo(:,c)-geo(:,b))+cd2d3(:)*ad2*tmp2 &
        & +cross_product(geo(:,b)-geo(:,a),geo(:,d)-geo(:,c))*ad2*tmp2 &
        & +cross_product(geo(:,a)-geo(:,c),cd2d3(:))*tmp1-cross_product(cd1d2(:),geo(:,d)-geo(:,c))*tmp1)*g/tmp3
        grad(:,3)=( dot_product(geo(:,b)-geo(:,a),cd2d3(:))/ad2*tmp2*(geo(:,c)-geo(:,b)) &
        & +cross_product(geo(:,b)-geo(:,a),geo(:,b)-geo(:,d))*ad2*tmp2 &
        & +cross_product(geo(:,b)-geo(:,a),cd2d3(:))*tmp1-cross_product(cd1d2(:),geo(:,b)-geo(:,d))*tmp1)*g/tmp3
        grad(:,4)=(cross_product(geo(:,c)-geo(:,b),cd1d2(:))*tmp1+cd1d2(:)*ad2*tmp2)*g/tmp3

        gradient(:,a) = gradient(:,a) + grad(:,1)
        gradient(:,b) = gradient(:,b) + grad(:,2)
        gradient(:,c) = gradient(:,c) + grad(:,3)
        gradient(:,d) = gradient(:,d) + grad(:,4)
    else
        write(*,*) "ERROR: strange gradient problem for -DH+ in int2xyz - please check this!" !MK
    endif

    return

end subroutine int2xyz




end module hp_hbcorr
