!MODULE DEFINITIONS
!==================================================================================================
!Module defintion
!By Runyu Zhang & Tennyson Jonathan 1/Sep/2022
!Contains sizes, diffs, old, logic, stream, mass, eqm, time
!Special notice:: The module name contains the file name, it means cannot directly paste to other 
!                 files and use. There are some difference between the constant numbers, types and 
!                 default value.
!==================================================================================================
module dipj0dvr_sizes
    save
    integer :: lbra0
    integer :: nbra0
    integer :: lket0
    integer :: nket0
    integer :: lbra1
    integer :: nbra1
    integer :: lket1
    integer :: nket1
    integer :: ntheta
    integer :: nr1
    integer :: nr2
    integer :: neval0
    integer :: neval1
end module dipj0dvr_sizes

module dipj0dvr_diffs
    save
    integer :: nqe
    integer :: nqo
    integer :: nr21
    double precision alphae
    double precision :: betae
    double precision :: alphao
    double precision :: betao
    logical zsame
end module dipj0dvr_diffs

module dipj0dvr_old
    save
    integer :: npta
    integer :: nptb
    integer :: nptc
    integer :: max2d    ! upper bound on size of intermediate 2d hamiltonian
    integer :: max3d    ! upper bound on size of full 3d hamiltonian
    integer :: npta1
    integer :: nptb1
    integer :: nptc1
    integer :: max2d1
    integer :: max3d1
    logical :: ztheta ! T let theta be first in the order of solution;
                                ! F let theta be last in the order of solution.
    logical :: zr2r1  ! T let r_2 come before r_1 in the order of solution;
                                ! F let r_1 come before r_2 in the order of solution. (only idia > -2).
    logical :: zthet1
    logical :: zr2r11
end module dipj0dvr_old

module dipj0dvr_logic
    save
    integer :: iptot
    integer :: idia
    logical :: zembed ! T z axis is along r2, = f z axis is along r1.
                    ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.
    logical :: zdone
end module dipj0dvr_logic

module dipj0dvr_stream
    save
    integer :: ibra0
    integer :: ibra1
    integer :: iket0
    integer :: iket1
    integer :: iwave0
    integer :: iwave1
    integer :: ivc0
    integer :: ivc1
    integer :: ione
    integer :: itwo
end module dipj0dvr_stream

module dipj0dvr_mass
    save
    double precision :: xmass(3)
    double precision :: g1
    double precision :: g2
end module dipj0dvr_mass

module dipj0dvr_eqm
    save
    double precision :: ex(3)
    double precision :: ez(3)
    double precision :: tmass
end module dipj0dvr_eqm

module dipj0dvr_time
    save
    double precision :: ouser
    double precision :: osys
    double precision :: ototal
end module dipj0dvr_time


program dipj0dvr
!
! This program calculates transition intensities between vibrations, by
! using the theory set out in Le Sueur, CR, Miller, S, Tennyson, J &
! Sutcliffe BT, Mol Phys (1992) 76,1147.  The vibrational wavefunctions
! have been obtained using DVR theory in the program DVR3D.
!

    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: navail, maxq, ir2, itheta, i, j
    double precision :: dipx, dipz, dx, dz
    parameter (navail=500000)
    parameter (maxq=500)
    character*80 title
    character c1,c2
    double precision, dimension(navail) :: array

! this subroutine reads all the controlling data for the job
    call rddata(title,array(1),array(1+maxq),array(1+maxq*2),array(1+maxq*3),array(1+maxq*4),array(1+maxq*5),maxq)
    do ir2=1,nr2
     array(nr1+ir2)=array(maxq+ir2)
    enddo
    do itheta=1,ntheta
     array(nr1+nr2+itheta)=array(maxq*2+itheta)
    enddo
    if (.not.zsame) then
     do ir2=1,nr21
        array(nr1+nr2+ntheta+ir2)=array(maxq*4+ir2)
     enddo
    endif
    call messge(title)         ! this subroutine writes the header for the job
    call gtmain(array,navail)  ! this subroutine allocates memory and then calls rest of the program
    do i=1,360
     call DIPD(dipx, 1.8242d0, 1.8242d0, cos(i*3.14159d0/180.0d0),1)
     call DIPD(dipz, 1.8242d0, 1.8242d0, cos(i*3.14159d0/180.0d0),0)
     write(20,*),i,sqrt(dipx*dipx+dipz*dipz)
    enddo
    do i=1,200
     call DIPD(dx,1.8242d0*i/100.0d0,1.8242d0*i/100.0d0,-.030527d0,1)
     call DIPD(dz,1.8242d0*i/100.0d0,1.8242d0*i/100.0d0,-.030527d0,0)
     write(21,*),i,sqrt(dx*dx+dz*dz)
    enddo
    do i=1,200
     j=201-i
     call DIPD(dx,1.8242d0*i/100.0d0,1.8242d0*j/100.0d0,-.030527d0,1)
     call DIPD(dz,1.8242d0*i/100.0d0,1.8242d0*j/100.0d0,-.030527d0,0)
     write(22,*),i,sqrt(dx*dx+dz*dz)
    enddo
end program dipj0dvr

!========================================================== rddata =====
! This subroutine reads in the data needed from the various places: a
! data file on stream 5, iwave0, ivc0, and if iptot=2, iwave1 and iv1.
subroutine rddata(title,r1,r2,theta,r11,r21,theta1,maxq)
    
    use dipj0dvr_logic
    use dipj0dvr_stream
    use dipj0dvr_mass
    use dipj0dvr_eqm
    use dipj0dvr_sizes
    use dipj0dvr_old
    use dipj0dvr_diffs
    namelist/prt/ibra0,iket0,iwave0,ivc0,ibra1,iket1,iwave1,ivc1,zsame,iptot,zdone,ione,itwo
    implicit none

    integer :: ipar, jrot, kmin, nlim, idia1, ipar1, ntheta1, nr11, jrot1, kmin1,nlim1,nthet1, ipt, maxq
    double precision :: tol, amtoau
    double precision, dimension(maxq) :: r1,r2,theta,r11,r21
    double precision :: theta1(maxq),xmass1(3)
    double precision r1e, r2e, xcose, r1e0, d1e0,w1e0,r2e0,d2e0,w2e0, g11, g21,r1e1,d1e1,w1e1,r2e1,d2e1,w2e1,rmass
    logical :: zmors1, zmors2, zemb1

    parameter (tol=1d-8)
    parameter (amtoau=1.6605402d-27/0.91093897d-30)
    character*80 title
  

    ibra0 = 51
    iket0 = 52
    iwave0 = 26
    ivc0 = 54
    ibra1 = 61
    iket1 = 62
    iwave1 = 63
    ivc1 = 64
    ione = 71
    itwo = 72
    iptot = 0

    zembed = .true.
    zsame = .true.
    zdone = .true.


    read(5,prt)
    read(5,"(a80)") title
    read(5,"(5i5)") lbra0,lket0,lbra1,lket1
    read(5,"(5i5)") nbra0,nket0,nbra1,nket1
    read(5,"(3f20.0)") r1e,r2e,xcose

    if (ibra0 == iket0 .or. ibra0 == iwave0 .or.                   &
    &   ibra0 == ibra1 .or. ibra0 == iket1 .or.                    &
    &   ibra0 == iwave1 .or. (.not.zdone .and.                       &
    &              (ibra0 == ivc0 .or. ibra0 == ivc1 .or.           &
    &               ibra0 == ione .or. ibra0 == itwo))) then
        write(6,*)' **** stream ibra0 clashes with other streams ****'
        stop
    elseif (iket0 == iwave0 .or. iket0 == ibra1 .or.               &
    &        iket0 == iket1 .or.  iket0 == iwave1 .or.              &
    &              (.not.zdone .and.                                    &
    &              (iket0 == ivc0 .or. iket0 == ivc1 .or.           &
    &               iket0 == ione .or. iket0 == itwo))) then
        write(6,*)' **** stream iket0 clashes with other streams ****'
        stop
    elseif (iwave0 == ibra1 .or. iwave0 == iket1 .or.              &
    &        iwave0 == iwave1 .or. (.not.zdone .and.                  &
    &              (iwave0 == ivc0 .or. iwave0 == ivc1 .or.         &
    &               iwave0 == ione .or. iwave0 == itwo))) then
        write(6,*)' **** stream iwave0 clashes with other streams ****'
        stop
    elseif (ibra1 == iket1 .or. ibra1 == iwave1 .or.               &
    &              (.not.zdone .and.                                    &
    &              (ibra1 == ivc0 .or. ibra1 == ivc1 .or.           &
    &               ibra1 == ione .or. ibra1 == itwo))) then
        write(6,*)' **** stream ibra1 clashes with other streams ****'
        stop
    elseif (iket1 == iwave1 .or. (.not.zdone .and.                   &
    &              (iket1 == ivc0 .or. iket1 == ivc1 .or.           &
    &               iket1 == ione .or. iket1 == itwo))) then
        write(6,*)' **** stream iket1 clashes with other streams ****'
        stop
    elseif (.not.zdone .and.                                           &
    &              (iwave1 == ivc0 .or. iwave1 == ivc1 .or.         &
    &               iwave1 == ione .or. iwave1 == itwo)) then
        write(6,*)' **** stream iwave1 clashes with other streams ****'
        stop
    elseif (.not.zdone .and.                                           &
    &              (ivc0 == ivc1 .or.                                 &
    &               ivc0 == ione .or. ivc0 == itwo)) then
        write(6,*)' **** stream ivc0 clashes with other streams ****'
        stop
    elseif (.not.zdone .and.                                           &
    &              (ivc1 == ione .or. ivc1 == itwo)) then
        write(6,*)' **** stream ivc1 clashes with other streams ****'
        stop
    elseif (.not.zdone .and.                                           &
    &               ione == itwo) then
        write(6,*)' **** stream ione clashes with other streams ****'
        stop
    endif

    if (lbra0 <= 0) lbra0=1
    if (lket0 <= 0) lket0=1
    if (lbra1 <= 0) lbra1=1
    if (lket1 <= 0) lket1=1
    if ((lbra0 > nbra0 .and. nbra0 /= -1) .or. nbra0 == 0) then
        lbra0=1
        nbra0=0
    endif

    if ((lket0 > nket0 .and. nket0 /= -1) .or. nket0 == 0) then
        lket0=1
        nket0=0
    endif

    if ((lbra1 > nbra1 .and. nbra1 /= -1) .or. nbra1 == 0) then
        lbra1=1
        nbra1=0
    endif

    if ((lket1 > nket1 .and. nket1 /= -1) .or. nket1 == 0) then
        lket1=1
        nket1=0
    endif

    if (nbra0 /= 0 .and. nket0 /= 0) then
        open (unit=ibra0,form='UNFORMATTED')
        open (unit=iket0,form='UNFORMATTED')
        open (unit=iwave0,form='UNFORMATTED')
        if (.not.zdone) open (unit=ivc0, form='UNFORMATTED')
    endif

    if (nbra1 /= 0 .and. nket0 /= 0) then
        open (unit=ibra1,form='UNFORMATTED')
        open (unit=iket1,form='UNFORMATTED')
        open (unit=iwave1,form='UNFORMATTED')
        if (.not.zdone) open (unit=ivc1, form='UNFORMATTED')
    endif

    if (.not.zdone) then
        open (unit=ione,form='UNFORMATTED')
        open (unit=itwo,form='UNFORMATTED')
    endif

    if (zdone) then
        write(*,*)  'DEBUG :: READING FROM IWAVE, LINE=185' ! Lorenzo Lodi 30-Oct-2009
        read(iwave0) idia,ipar,ntheta,nr1,nr2,jrot,kmin,neval0,nlim
        read(iwave0) zembed,zmors1,zmors2,xmass,g1,g2
        read(iwave0) r1e0,d1e0,w1e0,r2e0,d2e0,w2e0
    else
        write(*,*)  'DEBUG :: READING FROM IWAVE, LINE=190' ! Lorenzo Lodi 30-Oct-2009
        read(iwave0)idia,ipar,npta,nptb,nptc,max2d,max3d,neval0
        read(iwave0)zembed,zmors1,zmors2,ztheta,zr2r1,xmass,g1,g2
        if (ztheta) then
            ntheta=npta
            if (zr2r1) then
                nr2=nptb
                nr1=nptc
            else
                nr2=nptc
                nr1=nptb
            endif
        else
            ntheta=nptc
            if (zr2r1) then
                nr2=npta
                nr1=nptb
            else
                nr2=nptb
                nr1=npta
            endif
        endif
    endif

    if (idia /= 2 .or. (idia == 2 .and. iptot /= 2)) then
        iptot=0
        nbra1=0
        nket1=0
        neval1=0
    endif

    if (ntheta > maxq .or. nr1 > maxq .or. nr2 > maxq) then
        write(6,*)' **** parameter MAXQ is too small ****'
        write(6,*)'       resize it to at least ',max(ntheta,nr1,nr2)
        write(6,*)'       and recompile'
        stop
    endif

    if (nbra0 == -1) nbra0=neval0
    if (nket0 == -1) nket0=neval0
    write(*,*)  'DEBUG :: READING FROM IWAVE, LINE=225' ! Lorenzo Lodi 30-Oct-2009
    read(iwave0) ! Record skip to conform to new dvr3drjz format. Lorenzo Lodi 30-Oct-2009
    call getrow(r1,nr1,iwave0)
    call getrow(r2,nr2,iwave0)
    call getrow(theta,ntheta,iwave0)

    if (iptot == 2 .and. (nbra1 == -1 .or. nket1 == -1           &
    &    .or. nbra1 >= lbra1 .or. nket1 >= lket1)) then
        if (zdone) then
            read(iwave1)idia1,ipar1,ntheta1,nr11,nr21,jrot1,kmin1,neval1,   &
            &                                                             nlim1
            read(iwave1)zemb1,zmors1,zmors2,xmass1,g11,g21
            read(iwave1)r1e1,d1e1,w1e1,r2e1,d2e1,w2e1
        else
            read(iwave1)idia1,ipar1,npta1,nptb1,nptc1,max2d1,max3d1,neval1
            read(iwave1)zemb1,zmors1,zmors2,zthet1,zr2r11,xmass1,g11,g21
            if (ztheta) then
                nthet1=npta1
                if (zr2r1) then
                    nr21=nptb1
                    nr11=nptc1
                else
                    nr21=nptc1
                    nr11=nptb1
                endif
            else
                nthet1=nptc1
                if (zr2r1) then
                    nr21=npta1
                    nr11=nptb1
                else
                    nr21=nptb1
                    nr11=npta1
                endif
            endif
            if (nbra1 == -1) nbra1=neval1
            if (nket1 == -1) nket1=neval1
        endif

        call getrow(r11,nr11,iwave1)
        call getrow(r21,nr21,iwave1)
        call getrow(theta1,nthet1,iwave1)
        if (((zembed .and. .not.zemb1) .or. (.not.zembed .and. zemb1))     &
        &  .or. idia /= idia1 .or. ntheta /= nthet1 .or. nr1 /= nr11)   &
        &  then
            write(6,*)' **** The two runs are incompatible ****'
            stop
        elseif (zsame .and. nr2 /= nr21) then
            write(6,*)' **** The two runs are incompatible ****'
            stop
        elseif (abs(xmass(1)-xmass1(1))/xmass(1) > tol .or.            &
        &          abs(xmass(2)-xmass1(2))/xmass(2) > tol .or.            &
        &          abs(xmass(3)-xmass1(3))/xmass(3) > tol) then
            write(6,*)' **** The two runs use different masses ****'
            stop
        elseif (abs(g1-g11) > tol .or. abs(g2-g21) > tol) then
            write(6,*)' **** The two runs use different g''s ****'
            stop
        elseif (zsame) then
            do 10 ipt=1,nr1
                if (abs(r1(ipt)-r11(ipt))/r1(ipt) > tol) then
                    write(6,*)' **** The two runs use different quadrature points ****'
                stop
                endif
10        continue
        do 20 ipt=1,nr2
            if (abs(r2(ipt)-r21(ipt))/r2(ipt) > tol) then
              write(6,*)' **** The two runs use different quadrature points ****'
              stop
            endif
20          continue
            do 30 ipt=1,ntheta
            if (abs(theta(ipt)-theta1(ipt))/theta(ipt) > tol) then
              write(6,*)' **** The two runs use different quadrature points ****'
              stop
            endif
30          continue
        endif
    endif

! These data are extra for the special case of evens and odds being done
! with different basis functions in r2
    if (.not.zsame) then
        if (.not.zdone) then
            read(5,"(3f20.0)")r2e0,d2e0,w2e0
            read(5,"(3f20.0)")r2e1,d2e1,w2e1
        endif
        read(5,"(5i5)")nqe,nqo
        rmass=amtoau/(1d0/xmass(1)+g1*g1/xmass(2)+(1d0-g1)**2/xmass(3))
        betae=dsqrt(w2e0*rmass)
        betao=dsqrt(w2e1*rmass)
        alphae=dfloat(idint(4d0*d2e0/betae))
        alphao=dfloat(idint(4d0*d2e1/betao))
    endif

    tmass=xmass(1)+xmass(2)+xmass(3)
    call conver(r1e,r2e,xcose,ex,ez)

    return
end

!========================================================== message =====
! This subroutine writes a header message giving information about the
! program and the wavefunctions used.
subroutine messge(title)
    use dipj0dvr_logic
    use dipj0dvr_mass
    use dipj0dvr_eqm
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: i, maxq
    double precision :: x1, xsin, xcos, r1, r2, xl, xm, sumx, sumz, rmass
    character*80 title
    write(6,"(1x,78a1///5x,'Program DIPJ0DVR - version 1.1 (16 Oct 1992)'///1x,78a1///1x,a80///)")('*',i=1,78),('*',i=1,78),title
    if (g1 == 0d0) then
    ! bondlength-bondangle coordinates
        write(6,"(1x,'The DVR calculation was performed in ',a,' coordinates'//)")'bondlength-bondangle'
    elseif (g2 == 0d0) then
    ! scattering coordinates
        write(6,"(1x,'The DVR calculation was performed in ',a,' coordinates'//)")'scattering'
    else
    ! all other general coordinates
        write(6,"(1x,'The DVR calculation was performed in ',a,' coordinates'//)")'generalised'
    endif
    if (zembed) then
        write(6,"(a100, a100)")'The coordinate system used for the DVR calculation',& 
        & 'of the wavefunction had the z axis embedded along the r2 direction'
    else
        write(6,"(a100,a100)")'The coordinate system used for the DVR calculation',&
        &'of the wavefunction had the z axis embedded along the r1 direction'
    endif
    write(6,104)
    do 10 i=1,3
        write(6,"(1x,'Atom ',i1,' (mass ',f12.6,' amu):',5x,f13.6,7x,f13.6)")i,xmass(i),ex(i),ez(i)
10  continue
    if (iptot /= 2) then
        write(6,"(a100,i3,a10,i3,a20,i3,a10,i3,a20)") 'Transition intensities are to be calculated for transitions from the ',&
     &lket0,' to ',nket0,' initial states to the ',lbra0,' to ',nbra0,' final states'
    else
        write(6,"(a100)")'Transition intensities are to be calculated for transitions'
        if (nket0 /= 0 .and. nket1 /= 0) then
            write(6,"(1x,a20,i3,a10,i3,a20,i3,a10,i3)")'from even states ',lket0,' to ',nket0,' and odd states ',lket1,' to ',nket1
        elseif (nket0 == 0) then
            write(6,"(a20,i3,a10,i3)")'from odd states ',lket1,' to ',nket1
        elseif (nket1 == 0) then
            write(6,"(a20,i3,a10,i3)")'from even states ',lket0,' to ',nket0
        endif
        if (nbra0 /= 0 .and. nbra1 /= 0) then
            write(6,"(a10,i3,a10,i3,a20,i3,a10,i3)")'to even states ',lbra0,' to ',nbra0,' and odd states ',lbra1,' to ',nbra1
        elseif (nbra0 == 0) then
            write(6,"(a20,i3,a10,i3)")'to odd states ',lbra1,' to ',nbra1
        elseif (nbra1 == 0) then
            write(6,"(1x,'  to ',a4,' states ',i3,' to ',i3//)")'even',lbra0,nbra0
        endif
        if (.not.zsame) write(6,112)nqe,idint(alphae),betae,nqo,idint(alphao),betao
    endif
    return
104 format(1x,'The equilibrium values of r1, r2, and cos(theta)'/     &
     &       1x,'correspond to the following atom coordinates:'//       &
     &       42x,'x',20x,'z')
112 format(1x,'The even and odd blocks of the Hamiltonian used',        &
     &       ' different parameters for r_2.'/                            &
     &       1x,'An average will be taken of the transition',             &
     &       ' intensities obtained using'/                               &
     &       1x,'quadrature based on the two types of spherical',         &
     &       ' oscillators used in the'/                                  &
     &       1x,'DVR calculations.'/                                      &
     &       i3,' points are used for functions with alpha=',i3,          &
     &       ' and beta=',f12.6,', and '/                                 &
     &       i3,' points are used for functions with alpha=',i3,          &
     &       ' and beta=',f12.6//)
    end

!== ======================================================== conver =====
subroutine conver(r1,r2,xcos,x,z)
    
    use dipj0dvr_logic
    use dipj0dvr_mass
    use dipj0dvr_eqm
    implicit none
    integer :: x1, i
    double precision :: xsin, xcos, r2, r1, xl, xm, sumx, sumz
    parameter (x1=1d0)
    double precision x(3),z(3),zint
    do 10 i=1,3
        x(i)=0d0
        z(i)=0d0
10  continue
    xsin=dsqrt(x1-xcos*xcos)
    if (g1 == 0d0) then
    ! bondlength-bondangle coordinates
        if (zembed) then
            z(1)=r2
            z(2)=r1*xcos
            x(2)=r1*xsin
        else
            z(2)=r1
            z(1)=r2*xcos
            x(1)=r2*xsin
        endif
    elseif (g2 == 0d0) then
! scattering coordinates
        if (zembed) then
            z(1)=r2
            z(2)=(x1-g1)*r1*xcos
            x(2)=(x1-g1)*r1*xsin
            z(3)=-g1*r1*xcos
            x(3)=-g1*r1*xsin
        else
            z(2)=(x1-g1)*r1
            z(3)=-g1*r1
            z(1)=r2*xcos
            x(1)=r2*xsin
        endif
    else
! all other types of coordinates
        xl=(x1-g1)/(x1-g1*g2)
        xm=(x1-g2)/(x1-g1*g2)
        if (zembed) then
            z(1)=r2*xm
            z(2)=r1*xl*xcos
            x(2)=r1*xl*xsin
            zint=r2*(x1-xm)
            x(3)=-x(2)*g1/(x1-g1)
            z(3)=-(zint+z(2)*g1)/(x1-g1)
        else
            z(2)=r1*xl
            z(1)=r2*xm*xcos
            x(1)=-r2*xm*xsin
            z(3)=-(z(2)*g1+r2*(x1-xm)*xcos)/(x1-g1)
            x(3)=r2*(x1-xm)*xsin/(x1-g1)
        endif
    endif
    sumx=0d0
    sumz=0d0
    do 20 i=1,3
        sumx=sumx+xmass(i)*x(i)
        sumz=sumz+xmass(i)*z(i)
20  continue
        sumx=sumx/tmass
        sumz=sumz/tmass
! Move so that centre-of-mass is at the origin
    do 30 i=1,3
        x(i)=x(i)-sumx
        z(i)=z(i)-sumz
30  continue
    return
end
!== ======================================================== gtmain =====
! This routine allocates memory, & calls the real main program (inmain)
subroutine gtmain(array,navail)
    use dipj0dvr_logic
    use dipj0dvr_sizes
    use dipj0dvr_old
    use dipj0dvr_diffs
    implicit none
    integer :: ibegin, istart, i1, i2, i3, i4, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, &
    & i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, navail
    double precision ::  array(navail)
    ! check that enough memory has been allocated.
    ibegin=nr1+nr2+ntheta+1
    istart=ibegin
    if (.not.zsame) istart=istart+nr21
    ! size of wavefunction arrays
    ! phibra
    i1=istart+ntheta*nr1*max(nr2,nr21)
    ! phiket
    i2=i1+ntheta*nr1*max(nr2,nr21)
    ! size of dipole array
    ! dipx
    i3=i2+ntheta*nr1*nr2
    ! dipz
    i4=i3+ntheta*nr1*nr2
    if (.not.zsame .and.                                                &
    & (min(nbra0,nket1) /= 0 .or. min(nbra1,nket0) /= 0)) then
    ! dipx
        i3=i2+ntheta*nr1*nr2*nr21
    ! dipz
        i4=i3+ntheta*nr1*nr2
    !  extra arrays to produce dipx.  These can all be overwritten.
    ! transe
        j1=i4+nr2*nr2
    ! transo
        j2=j1+nr21*nr21
    ! basise
        j3=j2+nr2*(nqe+nqo)
    ! basiso
        j4=j3+nr21*(nqe+nqo)
    ! reo
        j5=j4+nqe+nqo
    ! q
        j6=j5+2*max(nr2,nr21,nqe,nqo)
    ! wt
        j7=j6+2*max(nr2,nr21,nqe,nqo)
    ! b
        j8=j7+max(nr2,nr21,nqe,nqo)
    ! c
        j9=j8+max(nr2,nr21,nqe,nqo)
    ! dnorme
        j10=j9+nr2
    ! dnormo
        j11=j10+nr21
    ! crunch
        j12=j11+nr2*nr21*(nqe+nqo)
    endif
    ! size of transition intensity arrays
    ! Tx
    i5=i4+(nbra0+nbra1-lbra0-lbra1+2)*(nket0+nket1-lket0-lket1+2)
    ! Tz
    i6=i5+(nbra0+nbra1-lbra0-lbra1+2)*(nket0+nket1-lket0-lket1+2)
! size of pointer arrays
! these can be overlaid over phiket onwards, as they are only needed to
! generate the wavefunctions which are then stored.
    if (.not.zdone) then
! j
    i7=i2+max(nptc,nptc1)
! iidum
    i8=i7+nptb*max(nptc,nptc1)
! k
    i9=i8+max(nptb,nptb1)*max(nptc,nptc1)
! c2d
    i10=i9+max(max2d,max2d1)
! ndim2d
    i11=i10+max(nptc,nptc1)
! c1d
    i12=i11+max(max2d,max2d1)*max(npta,npta1)
! c3d
    i13=i12+max(max3d,max3d1)
! cprod
    i14=i13+max(npta,npta1)*max(nptb,nptb1)*max(max3d,max3d1)
    else    !The following 2 lines were added by L. Lodi to avoid run-time error "used but not defined" on i14
        i14=0
    endif
    i15=max(i6,i14)
    if (.not.zsame) i15=max(i15,j12)
! size of energy array
! evals
    i16=i15+neval0+neval1
    write(6,*)' space needed is    ',i16
    write(6,*)' space available is ',navail
    if (i16 > navail) then
        write(6,*)' **** memory allocated is not large enough **** '
        write(6,*)'       change parameter NAVAIL to at least ',i16
        write(6,*)
        write(6,*)'       and then recompile and rerun'
        stop
    endif
    call inmain(array(istart),array(istart),array(i1),array(i2),       &
! That's......... phi,          phibra,       phiket,   dipx,
     &            array(i3),array(i4),array(i5),array(i2),array(i15),    &
! ............... dipz,     Tx,       Tz        j,        evals,
     &            array(1),array(nr1+1),array(nr1+nr2+1),array(ibegin))
! ............... r1,      r2,          theta,           r21
    return
end

!== ======================================================== inmain =====
! This is the real main program, which calls all the other subroutines
! involved in this calculation
subroutine inmain(phi,phibra,phiket,dipx,dipz,Tx,Tz,j,              &
     &            evals,r1,r2,theta,r21)

    use dipj0dvr_logic
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    use dipj0dvr_time
    implicit none

! note: the dimensions of c1d, c2d and c3d may not seem to be correct
! here but enough space has been left for them so it doesn't matter
    integer :: j
    double precision :: evals(neval0+neval1), phi, dipx, dipz, r1,r2, theta, tx, phibra, phiket, tz, r21
    call timer(' ')

! Read in wavefunctions
    if (max(nbra0,nket0) /= 0) then
        if (zdone) then
            call rdphi(phi,evals(1),0,nr2,neval0)
        else
            call oldphi(phi,j,evals(1),0,neval0)
        endif

        if (iptot == 2) then
            call timer( &
            & 'Reading in wavefunctions for even symmetry block took')
        else
            call timer('Reading in wavefunctions took')
        endif
    endif
    if (iptot == 2 .and. max(nbra1,nket1) /= 0) then

        if (zdone) then
            call rdphi(phi,evals(neval0+1),1,nr21,neval1)
        else
            call oldphi(phi,j,evals(neval0+1),1,neval1)
        endif

        call timer( &
        &  'Reading in wavefunctions for odd symmetry block took')
    endif
! calculate dipole (in Eckart coordinates) at each DVR point
    call getmu(dipx,dipz,r1,r2,theta,nr2)
    call timer('Calculating dipole at DVR points took')
! If even and odd symmetry blocks use different DVR points, need to get
! transformation matrices, and do integral over dipole & basis functions
!
! Tx is passed as the starting point of arrays needed to produce dipx
    if (.not.zsame .and. iptot == 2 .and..not. &
        & (max(nbra0,nket0) == 0 .or. max(nbra1,nket1) == 0)) then
        call diffmu(dipx,Tx,r1,theta,max(nqe,nqo,nr2,nr21))
    endif
! Calculate transition intensities
    call getT(phibra,phiket,dipx,dipz,evals,Tx,Tz,0,0,nr2)
    if (iptot == 2) then
        call timer( &
        &'Calculating transition intensities for even-even transitions took &
        &')
    else
        call timer('Calculating transition intensities took')
    endif
! If symmetry used other transition intensities must also be calculated
    if (iptot == 2 .and. zsame) then
        call getT(phibra,phiket,dipx,dipz,evals,Tx,Tz,0,1,nr2)
        call timer( &
        &'Calculating transition intensities for even-odd transitions took'  &
        &)
        call getT(phibra,phiket,dipx,dipz,evals,Tx,Tz,1,1,nr2)
        call timer( &
        &'Calculating transition intensities for odd-odd transitions took')
    elseif (iptot == 2 .and. max(nbra1,nket1) /= 0) then
! if even and odd blocks use different DVR points, calculating the
! transition intensities is slightly more complicated
        call diffT(phibra,phiket,dipx,evals,Tx,Tz)
        call timer( &
        &'Calculating transition intensities for even-odd transitions took' &
        &)
! if even and odd blocks use different DVR points, the dipole must be
! recalculated at the DVR points for the odd block
        call getmu(dipx,dipz,r1,r21,theta,nr21)
        call timer( &
        &'Calculating dipole at DVR points of odd Hamiltonian took')
        call getT(phibra,phiket,dipx,dipz,evals,Tx,Tz,1,1,nr21)
        call timer( &
        &'Calculating transition intensities for odd-odd transitions took')
    endif
! After all sections of the transition intensity matrix have been
! calculated, it can be written out.
    call writeT(Tx,Tz,evals)
    return
end

!== ======================================================== rdphi  =====
! This subroutine reads in phi(a,b,c) from the stream supplied by DVR3D
! (iwave) and then stores it in two places for the bra and the ket
subroutine rdphi(phi,evals,itime,mr2,neval)
    use dipj0dvr_stream
    use dipj0dvr_sizes
    implicit none
    integer :: itime, iwave, nbra, nket, jbra, jket, ibra, iket, meval, i, min, max, mr2, neval
    double precision :: autocm
    double precision :: phi(ntheta*nr1*mr2), evals(neval)
    autocm = 2.19474624D+05
    if (itime == 0) then
        iwave=iwave0
        neval=neval0
        nbra=nbra0
        nket=nket0
        jbra=lbra0
        jket=lket0
        ibra=ibra0
        iket=iket0
    else
        iwave=iwave1
        neval=neval1
        nbra=nbra1
        nket=nket1
        jbra=lbra1
        jket=lket1
        ibra=ibra1
        iket=iket1
    endif
    write(*,*)  'DEBUG :: READING FROM IWAVE, LINE=760', iwave ! Lorenzo Lodi 30-Oct-2009
    read(iwave0) ! Record skip to conform to new dvr3drjz format. Lorenzo Lodi 30-Oct-2009
    read(iwave0) ! Record skip to conform to new dvr3drjz format. Lorenzo Lodi 30-Oct-2009
    read(iwave) meval
    if (meval < neval) neval=meval
    call getrow(evals,neval,iwave)
    do 300 i=1,min(jbra,jket)-1
    call getrow(phi,ntheta*nr1*mr2,iwave)
300 continue
    do 310 i=min(jbra,jket),max(nbra,nket)
    call getrow(phi,ntheta*nr1*mr2,iwave)
    call outrow(phi,ntheta*nr1*mr2,ibra)
    call outrow(phi,ntheta*nr1*mr2,iket)
310 continue
    do 320 i=1,neval
    evals(i)=autocm*evals(i)
320 continue
    rewind ibra
    rewind iket
    return
end
!== ======================================================== oldphi =====
! This subroutine works out phi(a,b,c) from the 1D 2D and 3D
! coefficients, and then stores it in two places for the bra and the ket
subroutine oldphi(phi,array,evals,itime,neval)
    use dipj0dvr_logic
    use dipj0dvr_stream
    use dipj0dvr_mass
    use dipj0dvr_old
    use dipj0dvr_eqm
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: itime, ibra, iket, iout, ivc, mpta, mptb, mptc, n2d, n3d, jbra, nbra, jket, nket, &
        & istart, i1, i2, i3, i4, i5, i6, i7, i8, neval
    double precision :: phi, evals(neval)
    double precision :: array((2+2*nptb)*nptc+max2d+npta*max2d+max3d*(1+npta*nptb))
    logical :: zt, zr
    if (itime == 0) then
        ibra=ibra0
        iket=iket0
        iout=iwave0
        ivc=ivc0
        mpta=npta
        mptb=nptb
        mptc=nptc
        n2d=max2d
        n3d=max3d
        jbra=lbra0
        nbra=nbra0
        jket=lket0
        nket=nket0
        zt=ztheta
        zr=zr2r1
    else
        ibra=ibra1
        iket=iket1
        iout=iwave1
        ivc=ivc1
        mpta=npta1
        mptb=nptb1
        mptc=nptc1
        n2d=max2d1
        n3d=max3d1
        jbra=lbra1
        nbra=nbra1
        jket=lket1
        nket=nket1
        zt=zthet1
        zr=zr2r11
    endif
! set up arrays for working out phi
    istart=1
    i1=istart+mptc
    i2=i1+mptc*mptb
    i3=i2+mptc*mptb
    i4=i3+n2d
    i5=i4+mptc
    i6=i5+n2d*mpta
    i7=i6+n3d
    i8=i7+n3d*mpta*mptb
    call getphi(phi,array(istart),array(i1),array(i2),array(i3),      &
! That's........  phi,j,            iidum,    k,        c2d,
     &            array(i4),array(i5),array(i6),array(i7),evals,        &
! ..............  ndim2d,   c1d,      c3d,      cprod,
     &            n2d,n3d,mpta,mptb,mptc,neval,ibra,iket,iout,ivc,      &
     &            ione,itwo,jbra,nbra,jket,nket,zt,zr)
    return
end

!== ======================================================== getphi =====
subroutine getphi(phi,j,iidum,k,c2d,ndim2d,c1d,c3d,cprod, &
     & evals,n2d,n3d,npta,nptb,nptc,neval,ibra,iket,   &
     & iout,ivc,ivec1,ivec2,jbra,nbra,jket,nket,       &
     & ztheta,zr2r1)
    implicit none
    integer :: ivc, i, ndum, ibeta, igamma, ntot, nham2, ivec1, ivec2, ialpha, iout, jj, &
        & meval, jdum, kdum, kk, jbra, kbet, min, max, nbra, nket, ibra, iket, jket, n2d, &
        & n3d, npta, nptb, nptc, neval
    double precision :: utocm, autocm, temp
    logical :: ztheta, zr2r1
    integer :: j(nptc), iidum(nptb*nptc), k(nptb,nptc), ndim2d(nptc)
    double precision :: phi(npta*nptb*nptc), c2d(n2d),c1d(n2d,npta),&
     & c3d(n3d),evals(neval), cprod(n3d,npta,nptb)
    utocm=2.19474624D+05
! k(b,c) and j(c) mark the size of the 1D and 2D Hamiltonians
! respectively.  They are read in from a different stream (ivc)
    read(ivc)(IIDUM(i),i=1,NPTB*NPTC)
    read(ivc)(J(i),i=1,NPTC)
    NDUM=0
    DO 10 IBETA=1,NPTB
    DO 10 IGAMMA=1,NPTC
        NDUM=NDUM+1
        K(IBETA,IGAMMA)=IIDUM(NDUM)
10  CONTINUE
! Calculate the cumulative size of the 2D Hamiltonian for each value
! of gamma
    ntot=0
    DO 20 IGAMMA=1,NPTC
        NHAM2=0
        DO 30 IBETA=1,NPTB
            NHAM2=NHAM2+K(IBETA,IGAMMA)
30      CONTINUE
        NDIM2D(IGAMMA)=NHAM2
        ntot=ntot+nham2
20  CONTINUE
    if(n3d>ntot) n3d=ntot
! Write the 1D and 2D  Hamiltonian coefficients to different streams
! in order to facilitate the calculations
    rewind ivec1
    rewind ivec2
    DO 40 IGAMMA=1,NPTC
        IF (NDIM2D(IGAMMA) == 0) GOTO 40
        DO 50 IALPHA=1,NPTA
            CALL GETROW(C1D(1,ialpha),NDIM2D(IGAMMA),IOUT)
            if (j(igamma)/=0) &
             & CALL OUTROW(C1D(1,ialpha),NDIM2D(IGAMMA),IVEC1)
50      CONTINUE
40  CONTINUE
    read (iout)(IIDUM(i),i=1,NPTB*NPTC)
    DO 41 IGAMMA=1,NPTC
        DO 60 JJ=1,J(IGAMMA)
            CALL GETROW(C2D,NDIM2D(IGAMMA),IOUT)
            CALL OUTROW(C2D,NDIM2D(IGAMMA),IVEC2)
60      CONTINUE
41  CONTINUE
    read(iout)(J(i),i=1,NPTC)
    read(iout) meval
    if (meval < neval) neval=meval
    CALL GETROW(EVALS,NEVAL,IOUT)
    do 120 i=1,neval
        evals(i)=autocm*evals(i)
120 continue
! Now we are ready to start the calculation
    JDUM=0
    REWIND IVEC1
    REWIND IVEC2
    DO 150 IGAMMA=1,NPTC
        DO 160 JJ=1,J(IGAMMA)
            JDUM=JDUM+1
            IF (NDIM2D(IGAMMA) == 0) GOTO 170
            CALL GETROW(C2D,NDIM2D(IGAMMA),IVEC2)
170         DO 190 IALPHA=1,NPTA
                IF (NDIM2D(IGAMMA) == 0 .OR. JJ /= 1) GOTO 200
                CALL GETROW(C1D(1,ialpha),NDIM2D(IGAMMA),IVEC1)
200             KDUM=0
                DO 210 IBETA=1,NPTB
                TEMP=0D0
                    DO 220 KK=1,K(IBETA,IGAMMA)
                        KDUM=KDUM+1
                        temp=temp+C2D(KDUM)*C1D(KDUM,ialpha)
220                 CONTINUE
                cprod(jdum,ialpha,ibeta)=temp
210             CONTINUE
190         CONTINUE
160     CONTINUE
150 CONTINUE
    do 300 i=1,min(jbra,jket)-1
        CALL GETROW(C3D,n3D,iout)
300 continue
    DO 310 I=min(jbra,jket),max(nbra,nket)
        CALL GETROW(C3D,n3D,iout)
        do 320 ialpha=1,npta
            do 320 ibeta=1,nptb
            JDUM=0
                DO 330 IGAMMA=1,NPTC
                    temp=0d0
                    DO 340 JJ=1,J(IGAMMA)
                        JDUM=JDUM+1
                        temp=temp+cprod(jdum,ialpha,ibeta)*c3d(jdum)
340                 continue
                    if (ztheta) then
                        if (zr2r1) then
                            phi((igamma-1)*npta*nptb+(ialpha-1)*nptb+ibeta)=temp
                        elseif (.not. zr2r1) then
                            phi((igamma-1)*nptb*npta+(ibeta-1)*npta+ialpha)=temp
                        endif
                    elseif (.not. ztheta) then
                        if (zr2r1) then
                            phi((ialpha-1)*nptb*nptc+(ibeta-1)*nptc+igamma)=temp
                        else
                            phi((ialpha-1)*nptc*nptb+(igamma-1)*nptb+ibeta)=temp
                        endif
                    endif
330             continue
320         continue
        call outrow(phi,npta*nptb*nptc,ibra)
        call outrow(phi,npta*nptb*nptc,iket)
310  continue
    rewind ibra
    rewind iket
    return
    end
!== ======================================================== getmu  =====
! This subroutine calculates the dipoles along the x and z axes of the
! Eckart embedding coordinate system, and stores them
subroutine getmu(dipx,dipz,r1,r2,theta,mr2)
    use dipj0dvr_sizes
    implicit none
    integer :: ir2, ir1, itheta
    double precision :: dx, dz, mr2
    double precision :: dipx(ntheta,nr1,mr2),dipz(ntheta,nr1,mr2),r1(nr1), &
    & r2(nr2),theta(ntheta)
    do 10 ir2=1,mr2
    do 10 ir1=1,nr1
    do 10 itheta=1,ntheta
    call dipcal(r1(ir1),r2(ir2),theta(itheta),dx,dz)
    dipx(itheta,ir1,ir2)=dx
    dipz(itheta,ir1,ir2)=dz
10  continue
    return
end
!== ======================================================== dipcal =====
! This subroutine calls dipd and hence obtains the dipole components
! along the axes of the embedding used for the wavefunctions.  It then
! converts these to components along the axes of the Eckart embedding
! (see Le Sueur, CR,  Miller, S, Tennyson, J & Sutcliffe, BT; Mol Phys
! (1992) 76, 1147 for the theory behind this conversion)
subroutine dipcal(r1,r2,xcos,dipx,dipz)
    use dipj0dvr_mass
    use dipj0dvr_eqm
    implicit none
    integer :: i
    double precision :: r1, r2, xcos, dx, dz, top, bottom, angle, datam, datan, gcos, &
        & dcos, gsin, dsin, dipx, dipz
    double precision :: x(3),z(3)
    call conver(r1,r2,xcos,x,z)
    call dipd(dx,r1,r2,xcos,1)
    call dipd(dz,r1,r2,xcos,0)
    top=0d0
    bottom=0d0
    do 10 i=1,3
        top = top +xmass(i)*(ex(i)*z(i)-ez(i)*x(i))
        bottom=bottom+xmass(i)*(ex(i)*x(i)+ez(i)*z(i))
10  continue
    angle=datan(top/bottom)
    gcos=dcos(angle)
    gsin=dsin(angle)
    dipx=gcos*dx+gsin*dz
    dipz=gcos*dz-gsin*dx
    return
end
!== ======================================================== getT   =====
! This subroutine calculates the transition intensities from the
! wavefunctions and the dipoles
subroutine getT(phibra,phiket,dipx,dipz,evals,Tx,Tz,ibtime,iktime,mr2)
    use dipj0dvr_logic
    use dipj0dvr_stream
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: ibtime, iktime, ibra, iket, nbra, nket, jbra, jket, mket, min, max, mbra, &
        & lket, lbra, iabc, mr2
    double precision :: tix, tiz
    double precision :: phibra(ntheta*nr1*mr2),phiket(ntheta*nr1*mr2), &
        & ipx(ntheta*nr1*mr2), dipx(ntheta*nr1*mr2), dipz(ntheta*nr1*mr2), &
        & evals(neval0+neval1), &
        & Tx(nbra0+nbra1-lbra0-lbra1+2,nket0+nket1-lket0-lket1+2), &
        & Tz(nbra0+nbra1-lbra0-lbra1+2,nket0+nket1-lket0-lket1+2)
    if (ibtime == iktime .and. ibtime == 0) then
        ibra=ibra0
        iket=iket0
        nbra=nbra0
        nket=nket0
        jbra=lbra0
        jket=lket0
    elseif (ibtime == iktime) then
        ibra=ibra1
        iket=iket1
        nbra=nbra1
        nket=nket1
        jbra=lbra1
        jket=lket1
    else
        ibra=ibra0
        iket=iket1
        nbra=max(nbra0,nket0)
        nket=max(nbra1,nket1)
        jbra=min(lbra0,lket0)
        jket=min(lbra1,lket1)
    endif
    rewind iket
    mket=min(jbra,jket)
    if (ibtime /= iktime) mket=jket
    mbra=min(jbra,jket)
    if (ibtime /= iktime) mbra=jbra
    do 10 lket=mket,nket
        call getrow(phiket,ntheta*nr1*mr2,iket)
        if (lket < jket) goto 10
        if (ibtime /= iktime .and. lket > min(nbra1,nket1) .and. &
            & lket < max(lbra1,lket1)) goto 10

        rewind ibra
        do 20 lbra=mbra,nbra
            call getrow(phibra,ntheta*nr1*mr2,ibra)
            if (lbra < jbra) goto 20
            if (ibtime == iktime .and. lket < lbra .and. &
            & nbra >= lket .and. lket >= jbra .and. &
            & nket >= lbra .and. lbra >= jket) goto 20

            if (ibtime /= iktime .and..not.  &
            & ((nbra1 >= lket .and. lket >= lbra1 .and. &
            & nket0 >= lbra .and. lbra >= lket0) .or. &
            & (nket1 >= lket .and. lket>= lket1 .and. &
            & nbra0 >= lbra .and. lbra >= lbra0))) goto 20
            tix=0d0
            tiz=0d0
            if (idia == 2 .and. ibtime == iktime) then
                do 30 iabc=1,ntheta*nr1*mr2
                    tiz=tiz+phibra(iabc)*dipz(iabc)*phiket(iabc)
30              continue
            elseif (ibtime /= iktime) then
                do 40 iabc=1,ntheta*nr1*mr2
                    tix=tix+phibra(iabc)*dipx(iabc)*phiket(iabc)
40              continue
            else
                do 50 iabc=1,ntheta*nr1*mr2
                    tix=tix+phibra(iabc)*dipx(iabc)*phiket(iabc)
                    tiz=tiz+phibra(iabc)*dipz(iabc)*phiket(iabc)
50              continue
            endif

            if (ibtime == iktime) then
                if (lket <= nket .and. lket >= jket .and. &
                & lbra <= nbra .and. lbra >= jbra) then
                    Tx(ibtime*(nbra0-lbra0+1)+lbra-jbra+1, &
                    & iktime*(nket0-lket0+1)+lket-jket+1)=tix
                    Tz(ibtime*(nbra0-lbra0+1)+lbra-jbra+1, &
                    & iktime*(nket0-lket0+1)+lket-jket+1)=tiz
                endif
            if (lket <= nbra .and. lket >= jbra .and.  &
            & lbra <= nket .and. lbra >= jket) then
                Tx(ibtime*(nbra0-lbra0+1)+lket-jbra+1, &
                & iktime*(nket0-lket0+1)+lbra-jket+1)=tix
                Tz(ibtime*(nbra0-lbra0+1)+lket-jbra+1,  &
                & iktime*(nket0-lket0+1)+lbra-jket+1)=tiz
            endif
            else
                if (lket <= nket1 .and. lket >= lket1 .and. &
                & lbra <= nbra0 .and. lbra >= lbra0) then
                    Tx(lbra-lbra0+1,nket0+lket-lket0-lket1+2)=tix
                    Tz(lbra-lbra0+1,nket0+lket-lket0-lket1+2)=tiz
                endif
                if (lket <= nbra1 .and. lket >= lbra1 .and.  &
                & lbra <= nket0 .and. lbra >= lket0) then
                    Tx(nbra0+lket-lbra0-lbra1+2,lbra-lket0+1)=tix
                    Tz(nbra0+lket-lbra0-lbra1+2,lbra-lket0+1)=tiz
                endif
            endif
20      continue
10  continue
    return
end
!== ======================================================== diffmu =====
! This subroutine calculates the dipoles along the x and z axes of the
! Eckart embedding coordinate system, and stores them
    subroutine diffmu(dipx,start,r1,theta,maxq)
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: istart, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, maxq
    double precision :: dipx(ntheta,nr1,nr2,nr21),r1(nr1),theta(ntheta), &
        & start(nr2*nr2+nr21*nr21+(nr2+nr21+1)*(nqe+nqo)+ &
        & 6*maxq+nr2+nr21)
! set up arrays
    istart=1
! transe
    j1=istart+nr2*nr2
! transo
    j2=j1+nr21*nr21
! basise
    j3=j2+nr2*(nqe+nqo)
! basiso
    j4=j3+nr21*(nqe+nqo)
! reo
    j5=j4+nqe+nqo
! q
    j6=j5+2*max(nr2,nr21,nqe,nqo)
! wt
    j7=j6+2*max(nr2,nr21,nqe,nqo)
! b
    j8=j7+max(nr2,nr21,nqe,nqo)
! c
    j9=j8+max(nr2,nr21,nqe,nqo)
! dnorme
    j10=j9+nr2
! dnormo
    j11=j10+nr21
! crunch
    j12=j11+nr2*nr21*(nqe+nqo)
    call gettra(start(istart),start(j1),start(j5),start(j6), &
! That's......... transe,       transo,   q,        wt
    & start(j7),start(j8),start(j9),start(j10),maxq)
! ............... b,        c,        dnorme,    dnormo
    call timer('Calculating transformation matrices took')
    call getbas(start(j2),start(j3),start(j4),start(j5),start(j6), &
! That's......... basise,   basiso,   reo,      q,        wt,
    & start(j7),start(j8),start(j9),start(j10),maxq)
! ............... b,        c,        dnorme,   dnormo
    call timer('Calculating basis functions took')
    call domult(dipx,start(istart),start(j1),start(j2),start(j3),r1, &
! That's......... dipx,transe,       transo,   basise,   basiso,   r1,
    & start(j4),theta,start(j11))
! ............... r2,       theta,crunch
    call timer('Integrating over dipole and transforming back took')
    return
end
!== ======================================================== gettra =====
! this routine calculates the transformation matrices for the r2 points
subroutine gettra(transe,transo,q,wt,b,c,dnorme,dnormo,maxq)
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: ipt, ifn, maxq
    double precision :: csx, csa, tsx, tsa, prefac
    double precision :: toler, transe(nr2,nr2),transo(nr21,nr21), &
        & q(maxq),wt(maxq),b(maxq),c(maxq),dnorme(nr2),dnormo(nr21)
    toler=1d-8
    call glagpt(q,wt,b,c,alphae+0.5d0,csx,csa,tsx,nr2)
    call basis(transe,q,dnorme,alphae,nr2-1,nr2-1,nr2)
    tsa=1d0/(dnorme(1)*dnorme(1))
    if (abs((csx-tsx)/tsx) > toler .or. &
    & abs((csa-tsa)/tsa) > toler) then
        write(6,1000)'even',nr2,csx,csa,tsx,tsa
        stop
    endif
    do 10 ipt=1,nr2
        prefac=dsqrt(wt(ipt))
        do 20 ifn=1,nr2
            transe(ifn,ipt)=transe(ifn,ipt)*prefac*dnorme(ifn)
20      continue
10  continue
    call glagpt(q,wt,b,c,alphao+0.5d0,csx,csa,tsx,nr21)
    call basis(transo,q,dnormo,alphao,nr21-1,nr21-1,nr21)
        tsa=1d0/(dnormo(1)*dnormo(1))
    if (abs((csx-tsx)/tsx) > toler .or.  &
    & abs((csa-tsa)/tsa) > toler) then
        write(6,1000)'odd',nr21,csx,csa,tsx,tsa
        stop
    endif
    do 110 ipt=1,nr21
        prefac=dsqrt(wt(ipt))
        do 120 ifn=1,nr21
          transo(ifn,ipt)=transo(ifn,ipt)*prefac*dnormo(ifn)
120     continue
110 continue
    return
1000  format(1x,'Computed and analyti! sums of points and weights ',       &
     &       'do not agree'/                                               &
     &       1x,'for calculation of ',a,' transformation matrix using '    &
     &       ,i4,' points'/                                                &
     &       20x,'points',29x,'weights'/                                   &
     &       1x,'Computed',2x,e25.13,10x,e25.13/                           &
     &       1x,'Analytic',2x,e25.13,10x,e25.13)
      end
!== ======================================================== basis  =====
subroutine basis(rbasis,q,dnorm,A,idim,maxfn,npts)
    implicit none
    integer :: ipt, npt1, ifn, ifail1, ifail2, ifnm1, ifnm2, idim, maxfn, npts
    double precision :: a, alpha, amx, afnm1, bin, count, fn, rbasis(0:idim,npts),q(npts),dnorm(0:maxfn),s14aaf
    alpha=A+0.5d0
! Calculate normalisation constants, which are given by
!                 {(npts + alpha)               }
!  N_nalpha = sqrt{(   n + alpha) n! (npts - n)!}
!
! Note that the factor sqrt{(npts + alpha)!} comes from the quadrature
! weights. This transfer of factors is to avoid overflow problems.
!
    do 10 ipt=1,npts
        rbasis(0,ipt)=1d0
        rbasis(1,ipt)=alpha+1d0-q(ipt)
10  continue
    count=dfloat(npts)+alpha
    npt1=npts+1
    bin=1d0
    do 30 ifn=npts-1,maxfn+1,-1
        fn=dfloat(npts-ifn)
        bin=bin*count/fn
        count=count-1d0
30  continue
    ifail1=0
    ifail2=0
    do 40 ifn=maxfn,0,-1
        fn=dfloat(npts-ifn)
        bin=bin*count/fn
        dnorm(ifn)=dsqrt(bin*s14aaf(dfloat(ifn+1),ifail1) &
        & *s14aaf(fn+1d0,ifail2))
        count=count-1d0
40  continue
! Calculate the rest of the basis functions at each quadrature point,
! using arecursion formula
    amx=alpha+1d0
    do 50 ifn=2,maxfn
        fn=dfloat(ifn)
        amx=amx+2d0
        afnm1=alpha+fn-1d0
        ifnm1=ifn-1
        ifnm2=ifn-2
        do 50 ipt=1,npts
            rbasis(ifn,ipt)=((amx-q(ipt))*rbasis(ifnm1,ipt) &
            & -afnm1*rbasis(ifnm2,ipt))/fn
50  continue
    return
    end
!== ======================================================== s14aaf =====
    double precision FUNCTION S14AAF(xx,ifail)
!   Mock-up of NAG routine of the same names for
!   evaluating Gamma functions in the range (0,1)
    implicit none
    integer :: ifail 
    double precision :: xx, gammln
    s14aaf=exp(gammln(xx))
    ifail=0
    RETURN
    END
!== ======================================================== gammln =====
    double precision FUNCTION GAMMLN(XX)
!     Gamma function routine from Numerical Recipes p. 157
    implicit none
    integer :: j
    DOUBLE PRECISION :: COF(6), HALF,ONE,X,TMP,SER,xx, &
    & FPF(3)= (/0.5D0, 1.0D0, 5.5D0/), &
    stp(7) = (/76.18009173D0,-86.50532033D0,24.01409822D0,         &
     &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/)
    X=XX-ONE
    TMP=X+FPF
    TMP=(X+HALF)*LOG(TMP)-TMP
    SER=ONE
    DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11  CONTINUE
    GAMMLN=TMP+LOG(STP*SER)
    RETURN
    END
!== ======================================================== glagpt =====
subroutine glagpt(q,wt,b,c,alpha,csx,csa,tsx,npts)
! This subroutine calculates the quadrature points and weights for
! Gauss-Laguerre quadrature according to the routine given by Stroud and
! Secrest.  Roots are searched for by using a recursion formula
! x_k+2-x_k+1              1      { 1+2.55k   1.26kalpha }
! -----------  approx= ---------- { ------- + ---------- } for k > 0
! x_k+1-x_k            1+0.3alpha {   1.9k      1+3.5k   }
!
! results are improved on by using LGROOT which calls LGRECR.
!
! LGRECR also calculates L_n-1alpha and dL_nalpha/dx evaluated at x_i
! by using a recursion formula
!                             P_n=(x-b_n)P_n-1 - c_nP_n-2
! These are needed for the weights, given by
!
!                                  n! Gamma(n+alpha+1)
!                      wt_i = -----------------------------------
!                             dL_nalpha(x_i)/dx L_n-1alpha(x_i)
! (This is the formula given by S&S, but their program actually uses
! (n-1)! Gamma(n+alpha) on the top.  The formula given is a misprint.)
!
! However, in this program the weights are given as (n-1)!/(alpha+n) on
! the top instead & the extra factor of Gamma(n+alpha+1) is transferred
! to the normalisation constants.  This is in order to avoid overflow.
!
    implicit none
    double precision :: q(npts),wt(npts),b(npts),c(npts)
    integer :: ifail, ipt, npts
    double precision :: cc, alpha, fa, fn,eps = 1d-12, csx, csa, tsx, xt, step, pt, dpn, &
        & fi, r1, r2, ratio, pn1, xt2, pt2,s14aaf
! zero calculated sum of points
    csx=0d0
! zero calculated sum of weights
    csa=0d0
! evaluate numerator for weights
    ifail=0
    cc=s14aaf(dfloat(npts),ifail)/(dfloat(npts)+alpha)
    fa=alpha+1d0
    b(1)=fa
    c(1)=0d0
    fn=1d0
    do 10 ipt=2,npts
        fa=fa+2d0
        b(ipt)=fa
        c(ipt)=fn*(alpha+fn)
        fn=fn+1d0
10  continue
! evaluate exact sum of points
    tsx=fn*(alpha+fn)
! calculate roots.  First two are calculated using a different formula
! from that given by S&S.
    xt=(1d0+alpha)*(2d0+alpha)/(1d0+3d0*fn+2d0*alpha)
    step=3d0*(1d0+alpha)/(1d0+3d0*fn+alpha)
    call lgrecr(pt,dpn,pn1,xt,npts,alpha,b,c)
20  xt2=xt+step
    call lgrecr(pt2,dpn,pn1,xt2,npts,alpha,b,c)
    if (dsign(1d0,pt)*dsign(1d0,pt2) <= 0d0) then
        pt=pt2
        q(1)=0.5d0*(xt+xt2)
    else
        pt=pt2
        xt=xt2
        goto 20
    endif
    call lgroot(q(1),npts,alpha,dpn,pn1,b,c,eps)
    wt(1)=cc/dpn/pn1
    csx=csx+q(1)
    csa=csa+wt(1)
    xt=q(1)
25  xt2=xt+step
    call lgrecr(pt2,dpn,pn1,xt2,npts,alpha,b,c)
    if (dsign(1d0,pt)*dsign(1d0,pt2) <= 0d0) then
        pt=pt2
        q(2)=0.5d0*(xt+xt2)
    else
        pt=pt2
        xt=xt2
        goto 25
    endif
    call lgroot(q(2),npts,alpha,dpn,pn1,b,c,eps)
    wt(2)=cc/dpn/pn1
    csx=csx+q(2)
    csa=csa+wt(2)
    if (npts <= 2) return
    do 30 ipt=3,npts
        fi=dfloat(ipt-2)
        r1=(1d0+2.55d0*fi)/(1.9d0*fi)
        r2=1.26d0*fi*alpha/(1d0+3.5d0*fi)
        ratio=(r1+r2)/(1d0+0.3d0*alpha)
        q(ipt)=q(ipt-1)+ratio*(q(ipt-1)-q(ipt-2))
        call lgroot(q(ipt),npts,alpha,dpn,pn1,b,c,eps)
        wt(ipt)=cc/dpn/pn1
        csx=csx+q(ipt)
        csa=csa+wt(ipt)
30   continue
    return
    end
!== ======================================================== LGRECR =====
SUBROUTINE LGRECR(PN,DPN,PN1,X,NN,ALF,B,C)
!     USES RECURRENCE RELATIONS TO SET UP POLYNOMIALS
!     THIS ROUTINE IS DUE TO STROUD & SECREST
    implicit none
    integer :: j, nn
    double precision :: p1, p, x, alf, dp1, dp,q, dq,pn, dpn, pn1, &
        & B(NN), C(NN)
    P1 = 1.0D0
    P = X - ALF - 1.0D0
    DP1 = 0.0D0
    DP = 1.0D0
    DO 1 J=2,NN
        Q  = (X-B(J))* P-C(J)* P1
        DQ = (X-B(J))*DP-C(J)*DP1 + P
        P1 = P
        P  = Q
        DP1= DP
1       DP = DQ
    PN = P
    DPN= DP
    PN1= P1
    RETURN
    END
!== ======================================================== LGROOT =====
SUBROUTINE LGROOT(X,NN,ALF,DPN,PN1,B,C,EPS)
!     IMPROVES THE APPROXIMATE ROOT X; IN ADDITION OBTAINS
!          DPN = DERIVATIVE OF P(N) AT X
!          PN1 = VALUE OF P(N-1) AT X
!     THIS ROUTINE IS DUE TO STROUD & SECREST
    IMPLICIT none
    integer :: itmax, iter, nn
    double precision :: B(NN),C(NN), p, dpn, pn1, x, alf, d, eps
    itmax=10
    ITER=0
1   ITER=ITER+1
    CALL LGRECR(P,DPN,PN1,X,NN,ALF,B,C)
    D = P/DPN
    X = X-D
    IF (ABS(D/X) <= EPS) RETURN
    IF (ITER - ITMAX) 1,2,2
2   WRITE(6,"(a50,I4,a20,D26.15,a20,D26.15)")'WARNING: NOCONVERGENCE AFTER',ITER,' ITERATIONS',&
        &'CURRENT DIFFERENCE',D,' & ROOT',X
    RETURN
    END

!== ======================================================== getbas =====
! this subroutine finds the quadrature points to be used in evaluating
! the dipole
subroutine getbas(basise,basiso,reo,q,wt,b,c,dnorme,dnormo,maxq)
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: ipt, ifaile, ifailo, ifn, maxq
    double precision :: toler=1d-8, csxe, csae, tsxe, csxo, csao,tsxo, tsae, factor, prefac, tsao,dsqrt,s14aaf
    double precision ::  basise(nr2,nqe+nqo),basiso(nr21,nqe+nqo), &
            & reo(nqe+nqo),q(2*maxq),wt(2*maxq),b(maxq),c(maxq),   &
            & dnorme(nr2),dnormo(nr21)
    if (nqe > 0) call glagpt(q,wt,b,c,alphae+0.5d0,csxe,csae,tsxe,nqe)
    do 100 ipt=1,nqe
        reo(ipt)=dsqrt(q(ipt)/betae)
100 continue
    if (nqo > 0) call glagpt(q(1+nqe),wt(1+nqe),b,c,alphao+0.5d0,csxo,csao,tsxo,nqo)
    do 200 ipt=1,nqo
        reo(ipt+nqe)=dsqrt(q(ipt+nqe)/betao)
200 continue
    if (nqe > 0) then
        call basis(basise(1,1),q,dnorme,alphae,nr2-1,nr2-1,nqe)
        tsae=1d0/(dnorme(1)*dnorme(1))
        if (abs((csxe-tsxe)/tsxe) > toler .or.  &
        & abs((csae-tsae)/tsae) > toler) then
            write(6,1000)'even',nqe,csxe,csae,tsxe,tsae
            stop
        endif
        call basis(basiso(1,1),q,dnormo,alphao,nr21-1,nr21-1,nqe)
    endif
    ifaile=0
    ifailo=0
    if (nqe > 0)  &
     & factor=dsqrt(s14aaf(dfloat(nqe)+alphae+1.5d0,ifaile)/ &
     & s14aaf(dfloat(nqe)+alphao+1.5d0,ifailo))

        do 110 ipt=1,nqe
            prefac=dsqrt(wt(ipt))
            do 120 ifn=1,nr2
                basise(ifn,ipt)=basise(ifn,ipt)*prefac*dnorme(ifn)
120         continue
            do 130 ifn=1,nr21
                basiso(ifn,ipt)=basiso(ifn,ipt)*prefac*dnormo(ifn) &
                & *q(ipt)**((alphao-alphae)/2d0)  &
                & *factor
130         continue
110     continue

        if (nqo > 0) then
        call basis(basise(1,1+nqe),q(1+nqe),dnorme,alphae,nr2-1, &
        & nr2-1,nqo)
        call basis(basiso(1,1+nqe),q(1+nqe),dnormo,alphao,nr21-1, &
        & nr21-1,nqo)

        tsao=1d0/(dnormo(1)*dnormo(1))
        if (abs((csxo-tsxo)/tsxo) > toler .or.  &
        & abs((csao-tsao)/tsao) > toler) then
            write(6,1000)'odd',nqo,csxo,csao,tsxo,tsao
            stop
        endif
    endif
    ifaile=0
    ifailo=0
    if (nqo > 0) &
     &factor=dsqrt(s14aaf(dfloat(nqo)+alphao+1.5d0,ifailo)/ &
     &s14aaf(dfloat(nqo)+alphae+1.5d0,ifaile))
        do 210 ipt=1,nqo
            prefac=dsqrt(wt(ipt+nqe))
            do 220 ifn=1,nr2
                basise(ifn,ipt+nqe)=basise(ifn,ipt+nqe)*prefac*dnorme(ifn) &
                & *q(ipt+nqe)**((alphae-alphao)/2d0)       &
                & *factor
220         continue
            do 230 ifn=1,nr21
                basiso(ifn,ipt+nqe)=basiso(ifn,ipt+nqe)*prefac*dnormo(ifn)
230         continue
210     continue
    return
1000  format(1x,'Computed and analyti! sums of points and weights ',     &
     &       'do not agree'/                                             &
     &       1x,'for calculation of ',a,' transformation matrix using '  &
     &       ,i4,' points'/                                              &
     &       20x,'points',29x,'weights'/                                 &
     &       1x,'Computed',2x,e25.13,10x,e25.13/                         &
     &       1x,'Analytic',2x,e25.13,10x,e25.13)
    end
!== ======================================================== domult =====
! this routine calculates transition moments for the specific case when
! r2 is run last, and the bra and the ket use different functions for r2
subroutine domult(dipx,transe,transo,basise,basiso,r1,r2,theta,crunch)
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: ico, ic1, ipt, ic0, ia, ib
    double precision :: divide, temp, nc0, even, odd, nc1, dx, dz, &
        & dipx(ntheta,nr1,nr2,nr21),transe(nr2,nr2), &
        & transo(nr21,nr21),basise(nr2,nqe+nqo), &
        & basiso(nr21,nqe+nqo),r1(nr1),r2(nqe+nqo),theta(ntheta), &
        & crunch(nr2,nr21,nqe+nqo)
    divide=2d0
    if (nqe == 0 .or. nqo == 0) divide=1d0
        do 100 ic0=1,nr2
            do 100 ic1=1,nr21
                do 100 ipt=1,nqe+nqo
                    temp=0d0
                    do 110 nc0=1,nr2
                        even=transe(nc0,ic0)*basise(nc0,ipt)
                        odd=0d0
                        do 120 nc1=1,nr21
                            odd=odd+basiso(nc1,ipt)*transo(nc1,ic1)
120                     continue
                        temp=temp+even*odd
110                 continue
                    crunch(ic0,ic1,ipt)=temp/divide
100     continue
    do 10 ia=1,ntheta
        do 10 ib=1,nr1
            do 20 ic0=1,nr2
                do 20 ic1=1,nr21
                    dipx(ia,ib,ic0,ic1)=0d0
20          continue
            do 30 ipt=1,nqe+nqo
                call dipcal(r1(ib),r2(ipt),theta(ia),dx,dz)
                do 40 ic1=1,nr21
                    do 40 ic0=1,nr2
                        dipx(ia,ib,ic0,ic1)=dipx(ia,ib,ic0,ic1)+crunch(ic0,ic1,ipt)*dx
40              continue
30             continue
10  continue
    return
    end

!== ======================================================== diffT  =====
! this routine calculates transition moments for the specific case when
! r2 is run last, and the bra and the ket use different functions for r2
subroutine diffT(phibra,phiket,dipx,evals,Tx,Tz)
    use dipj0dvr_stream
    use dipj0dvr_sizes
    use dipj0dvr_diffs
    implicit none
    integer :: ico, ic1, iab, ic0
    double precision :: ne, no, tix, &
        & phibra(ntheta*nr1,nr2),phiket(ntheta*nr1,nr21), &
        & dipx(ntheta*nr1,nr2,nr21),evals(neval0+neval1), &
        & Tx(nbra0+nbra1-lbra0-lbra1+2,nket0+nket1-lket0-lket1+2), &
        & Tz(nbra0+nbra1-lbra0-lbra1+2,nket0+nket1-lket0-lket1+2)
    rewind ibra0
    do 200 ne=min(lbra0,lket0),max(nbra0,nket0)
        call getrow(phibra,ntheta*nr1*nr2,ibra0)
        if (ne > min(nbra0,nket0) .and. ne < max(lbra0,lket0)) goto 200
        rewind iket1
        do 210 no=min(lbra1,lket1),max(nket1,nbra1)
            call getrow(phiket,ntheta*nr1*nr21,iket1)
            if (.not.((ne <= nbra0 .and. ne >= lbra0 .and. &
            & no <= nket1 .and. no >= lket1) .or. &
            & (ne <= nket0 .and. ne >= lket0 .and. &
            & no <= nbra1 .and. no >= lbra1))) goto 210
            tix=0d0
            do 220 ic0=1,nr2
                do 220 ic1=1,nr21
                    do 230 iab=1,ntheta*nr1
                        tix=tix+phibra(iab,ic0)*dipx(iab,ic0,ic1)* phiket(iab,ic1)
230                 continue
220         continue
            if (ne <= nbra0 .and. ne >= lbra0 .and. &
            & no <= nket1 .and. no >= lket1) then
                Tx(ne-lbra0+1,nket0+no-lket0-lket1+2)=tix
                Tz(ne-lbra0+1,nket0+no-lket0-lket1+2)=0d0
            endif
            if (ne <= nket0 .and. ne >= lket0 .and. &
            & no <= nbra1 .and. no >= lbra1) then
                Tx(nbra0+no-lbra0-lbra1+2,ne-lket0+1)=tix
                Tz(nbra0+no-lbra0-lbra1+2,ne-lket0+1)=0d0
            endif
210     continue
200 continue
    return
    end

!== ======================================================== writeT =====
! This subroutine writes out the transition intensities in the same
! format as DIPOLE does.
subroutine writeT(Tx,Tz,evals)
    use dipj0dvr_logic
    use dipj0dvr_sizes
    implicit none
    integer :: ie1, i1, ie2, i2, istream, nsize
    double precision :: autocm, autode, detosec, gzz, e1, xe1, e2, xe2, &
        & dd, dd3, sx, tzd,txd, sxd, t, a, vec, &
        & Tx(nbra0+nbra1-lbra0-lbra1+2,nket0+nket1-lket0-lket1+2),  &
        & Tz(nbra0+nbra1-lbra0-lbra1+2,nket0+nket1-lket0-lket1+2),  &
        & evals(neval0+neval1)
    character c1, c2;
    AUTOCM= 2.19474624D+05
! AUTODE CONVERTS ATOMIC UNITS TO DEBYE
    AUTODE= 2.5417662D0
! DETOSEC CONVERTS FROM S(F-I) IN DEBYE**2 TO SECONDS-1
    DETOSEC= 3.136186D-07
    WRITE(6,"(/,/,/)")
    WRITE(6,201)
    WRITE(6,"(/,/,/)")
    GZZ= evals(1)
    WRITE(6,"(5X,'GROUND ZERO =',E16.8,' CM-1')") GZZ
    WRITE(6,"(/,/,/)")
    WRITE(6,"(a100, a100)")' IE1   IE2     KET ENERGY   BRA ENERGY    FREQUENCY  Z',&
      &' TRANSITION    X TRANSITION       DIPOLE       S(F-I)      A-COEFFICIENT'
!   CALCULATE TRANSITION MOMENTS, LINE STRENGTHS AND A-COEFFICIENTS
    DO 1 IE1=1,Nbra0+nbra1-lbra0-lbra1+2
        c1=' '
        if (iptot == 2 .and. ie1 <= nbra0-lbra0+1) c1='e'
        if (iptot == 2 .and. ie1 > nbra0-lbra0+1) c1='o'
        i1=ie1+lbra0-1
        e1=evals(i1)
        if (ie1 > nbra0-lbra0+1) then
            i1=ie1-nbra0+lbra0-1+lbra1-1
            e1=evals(i1+neval0)
        endif
        XE1= (E1 - GZz)
        DO 2 IE2=1,Nket0+nket1-lket0-lket1+2
            c2=' '
            if (iptot == 2 .and. ie2 <= nket0-lket0+1) c2='e'
            if (iptot == 2 .and. ie2 > nket0-lket0+1) c2='o'
            i2=ie2+lket0-1
            e2=evals(i2)
            if (ie2 > nket0-lket0+1) then
                i2=ie2-nket0+lket0-1+lket1-1
                e2=evals(i2+neval0)
            endif
            XE2= (E2 - GZz)
            DD= XE2 - XE1
            DD3= ABS(DD*DD*DD)
            sx=tz(ie1,ie2)**2
            sx=sx+tx(ie1,ie2)**2
            tzd=tz(ie1,ie2)*autode
            TXD= Tx(IE1,IE2)*AUTODE
            SXD= TZD*TZD + TXD*TXD
            t=sqrt(sxd)
            A= SXD*DD3*DETOSEC
            WRITE(6,"(2(I4,1x,a1),3(3X,F10.3),5(2X,E13.6))") I1,c1,I2,c2,XE1,XE2,DD,TZD,TXD,T,SXD,A
2       CONTINUE
        WRITE(6,"(//)")
1   CONTINUE
    return
!200 FORMAT(/,/,/)
201 FORMAT(1H1,5X,'*************************************************'   &
     &      ,/,/,5X,'PRINT OUT OF DIPOLE TRANSITION MOMENTS AND S(F-1)'   &
     &      ,/,/,9X,'FREQUENCIES IN WAVENUMBERS',                         &
     &         /,9X,'TRANSITION MOMENTS IN DEBYE (2.54174A.U.)',          &
     &         /,9X,'S(F-I) IN DEBYE**2',                                 &
     &         /,9X,'EINSTEIN A-COEFFICIENT IN SEC-1',/,/,                &
     &      5X,'*************************************************')
!204 FORMAT(5X,'GROUND ZERO =',E16.8,' CM-1')
!205 FORMAT(/,' IE1   IE2     KET ENERGY   BRA ENERGY    FREQUENCY  Z' &
!      & // ' TRANSITION    X TRANSITION       DIPOLE       S(F-I)      A-COEFFICIENT',/)
!206 FORMAT(2(I4,1x,a1),3(3X,F10.3),5(2X,E13.6))
!207 FORMAT(//)
    end

!== ======================================================== GETROW =====
SUBROUTINE GETROW(VEC,NSIZE,ISTREAM)
    IMPLICIT none
    integer :: nsize, istream
    double precision :: vec
    DIMENSION VEC(NSIZE)
    READ(ISTREAM) VEC
    RETURN
END
!== ======================================================== OUTROW =====
SUBROUTINE OUTROW(VEC,NSIZE,ISTREAM)
    IMPLICIT none
    integer :: nsize, istream
    double precision :: vec
    DIMENSION VEC(NSIZE)
    WRITE(ISTREAM) VEC
    RETURN
END
!== ======================================================== timer  =====
subroutine timer(text)
    use dipj0dvr_time
    implicit none
    character text*(*)
    real*4 etime,inewtime(2)
    double precision total, user, system
    total = 0.0
    user = 0.0
    system = 0.0
!     total=dble(etime(inewtime))
!     user=dble(inewtime(1))
!     system=dble(inewtime(2))
    if (ototal /= 0d0) write(6,1000)text, total-ototal
    ouser=user
    osys=system
    ototal=total
    return
1000  format(1x,a,/,10x,f20.3,' seconds of user time',  &
        & /,10x,f20.3,' seconds of system time', &
        & /,10x,f20.3,' seconds in total',/)
      end
!========================================================== FCT    =====
double precision FUNCTION FCT(R1,R2,R3,NPA,C,INDEX,RO)
    IMPLICIT none
    integer :: i, j, npa
    double precision :: r2, r3, r1, f, C(NPA),INDEX(3,NPA), Q(3),R(3),RO(3)
    R(1)= R1
    R(2)= R2
    R(3)= R3
    DO 3 I=1,3
3   Q(I)= R(I) - RO(I)
    FCT= 0.
        DO 1 I=1,NPA
            F=1.
            DO 2 J=1,3
            IF (INDEX(J,I)/=0) THEN
                F=F*Q(J)**INDEX(J,I)
            END IF
2           CONTINUE
            FCT=FCT+C(I)*F
1       CONTINUE
    RETURN
    END
