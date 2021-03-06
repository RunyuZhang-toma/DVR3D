!=================================================================
! /H3+Jacobi.f90
! Ediss, Esp, m, ga1

!==================================================================
!Section 001
!/F90 file
!/HCN/dipole3.f90
!/H3+Jacobi
!/H3+Radau
!/water
!/source/dipj0dvr.f90
!module template /source/dipolez.f90
!new common constant dvr3drjz.f90

    implicit None
    ! integer /stream/ /timing/ /dim/ /sym/ /potiential/ /pot/
    integer :: &
            ezero, &
            ibase1, ibase2, ibra, ibra0,ibra1, idia, ifit, iket, ione, ipar1, ipar2, ipot, iprop, iptot, &
                ires, iscr, itera, itra, itime0, itwo, iwave0,iwave1, ivc0, ivc1, &
            jk1, jk2, iket0,iket1, jrot1, jrot2 &
            lket0, lket1, lbra0, lbra1 &
            kmin1, kmin2, &
            max2d, max2d1, max3d, max3d1, mbass, mbass1, mbass2, mblock, &
            nbin, nblock, nbmax1, nbmax2, nbra0, nbra1, ncoord, neval0, neval1, neval2, &
            nket0, nket1, nn2, npnt, npta, npta1, nptb, nptb1, nptc, nptc1, npnt1, npnt2, npot, nqe, nqo, &
            nr1, nr2, nr21, nrade, nrado, ntheta,


    ! logical /logic/ /mass/
    logical :: &
            zbisc, &
            zdone, &
            zembed, zezmbed, &
            zmors1, zmors2, &
            znco1, znco2, &
            zpmin, zprint, zpseg, &
            zrme1, zrme2, zrme3, zr2r1, zr2r11, &
            zsame, zstart,  &
            ztheta, zthet1, ztra,  &
            zuvvis
             


    ! real /head/ /mass/ /potential/
    real(kind=dp) :: 
            alphao, alphae, &
            betao, , betae, &
            cv, &
            der, &
            ex(3), ez(3), &
            g1, g2, &
            ototal, osys, ouser, &
            PI, &
            title, tmass, &
            xm1(3), xm2(3), xmass(3), xmassr(3)
                      

    ! character
    character(len=8) :: title(9)

!-------------------------------------------------------------------

!    zmors1 = .true.
!    zprint = .false.
!    ztra = .true.
!    zmors2 = .true.
!    zpmin = .false.
!    zstart = .false.

!    ires = 0
!    nblock = 1000
!    iket = 11
!    ibra = 12
!    itera = 13
!    iscr = 24

!    cv = 100
!    ifit = 100
!    der = 100

