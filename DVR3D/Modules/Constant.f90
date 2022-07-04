Module Constant

!==================================================================
!Section 000
!Copyright(C), 2022, UCL
!FileName: Example.f90
!Author: Runyu Zhang
!Data: 04/07/2022
!Description: This is a test mudole file for dipole3.f90

!==================================================================
!Section 001

    implicit None
    ! integer /stream/ /timing/ /dim/ /sym/
    integer :: iket, ibra, itera, iscr, ires, mblock, nblock, &
               itime0 &
               ncoord, npnt, npnt1, npnt2, nrade, nrado, npot, nbin, &
               nbmax1, nbmax2, mbass, mbass1, mbass2, kmin1, kmin2,  &
               jk1, jk2, neval1, neval2, nn2, ibase1, ibase2, ipot   &
               idia, ipar1, ipar2, jrot1, jrot2                      &

    ! logical /logic/ /mass/
    logical :: zmors1, znco1, znco2, zprint, zpmin, ztra, zstart, zmors2 &
               zembed, zbisc

    ! real /head/ /mass/
    real(kind=dp) :: title &
                     xmass(3), xm1(3), xm2(3), g1, g2

    ! character
    character(len=8) :: title(9)

!===================================================================
!Section 002

    zmors1 = .true.
    zprint = .false.
    ztra = .true.
    zmors2 = .true.
    zpmin = .false.
    zstart = .false.

    ires = 0
    nblock = 1000
    iket = 11
    ibra = 12
    itera = 13
    iscr = 24

