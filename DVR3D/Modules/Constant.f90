Module Constant

!==================================================================
!Section 000
!Copyright(C), 2022, UCL
!FileName: Constant.f90
!Author: Runyu Zhang
!Data: 04/07/2022
!Description: This is a test mudole file for dipole3.f90

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


!==================================================================

!EXample module template
!MODULE DEFINITIONS--------------------------------------------

module input
!  definition of the control input parameters
  save
  !line 1
  !namelist prt
  logical :: zprint ! =T supplies extra print out for debugging purposes
  logical :: ztra   ! =T writes data for spectra to stream ITRA
  logical :: zstart ! =T initiates the output file for the data for SPECTRA
  integer :: iket   ! input stream from DVR3DRJZ/ROTLEV3/ROTLEV3B for the ket (unformated)
  integer :: ibra   ! input stream for the bra (unformmatted)
  integer :: itra   ! output stream to SPECTRA (if ZTRA = T ) (unformated)
  namelist /prt/ zprint, ztra, zstart, iket, ibra, itra    
  
  !line 2
  character(len=72) :: title
  
  !line 3
  integer :: npot   ! number of gauss legendre integration points
  integer :: nv2    ! number of ket eigenfunctions considered
  integer :: nv1    ! number of bra eigenfunctions considered
  integer :: ibase2 ! number of lowest ket eigenfunctions skipped
  integer :: ibase1 ! number of lowest bra eigenfunctions skipped
  
  !line 4
  real(8) :: ezero  ! the ground state of the system in cm-1
end module input


!==================================================================================================

module dvr3drjz_file
! All common data from file dvr3drjz.f90
    save

    implicit none

!logical variable
! /prt/ /outp/
logical :: zpham = .false.      ! T request printing of the hamiltonian matrix
logical :: zprad = .false.      ! T request printing of the radial matrix elements
logical :: zpvec = .false.      ! T request printing of the eigenvectors
logical :: zrot  = .true.       ! F do vibrational part of rotational calculation by looping over k
logical :: zladd = .true.       ! T NALF kept constant as k increases
                                ! F NALF decreases with k (=f has a bug), (only if zrot = .true.)
logical :: zembed= .true.       ! T z axis is along r2, = f z axis is along r1.
                                ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.
logical :: zmors1= .true.       ! T use morse oscillator-like functions for r_1 coordinate;
                                ! F use spherical oscillator functions.
logical :: zmors2= .true.       ! T use morse oscillator-like functions for r_2 coordinate;
                                ! F use spherical oscillator functions.
logical :: zpmin = .false.      ! T requests only minimal printing.
logical :: zvec  = .false.      ! T store the eigenvectors from all the parts of the calculation
                                ! (1d,2d and 3d) on stream iout2.
                                ! further information relating to this (arrays iv1 and iv2) is
                                ! stored on stream iout1.
logical :: zquad2= .true.       ! T use the dvr quadrature approximation for the integrals of
                                ! the r_2^{-2} matrix, and hence make its dvr transformation diagonal.
                                ! F evaluate the r_2^{-2} integrals fully and perform the
                                ! full dvr transformation on them.
                                ! Note that zquad2 = f is not implemented for zmors2 = t
                                ! or for ztheta = f.
logical :: zdiag = .true.       ! F do not do final diagonalisation, instead the final Hamiltonian
                                ! matrix is written on units IDIAG1 and IDIAG2. 
logical :: zlmat = .false.      ! T requests printing of the L-matrix.
logical :: zcut  = .false.      ! T final dimension selected using an energy cut-off given  by emax2.
                                ! F final dimension determined by nham3.
logical :: zall  = .false.      ! T requests no truncation of the intermediate solution.
logical :: zlin  = .false.      ! T NALF kept constant as k increases
                                ! F NALF decreases with k (=f has a bug), (only if zrot = .true.)
logical :: zp1d  = .false.      ! T requests printing of the results of 1d calculations.
logical :: zp2d  = .false.      ! T requests printing of the results of 2d calculations.
logical :: zr2r1 = .true.       ! T let r_2 come before r_1 in the order of solution;
                                ! F let r_1 come before r_2 in the order of solution. (only idia > -2).
logical :: ztheta= .true.       ! T let theta be first in the order of solution;
                                ! F let theta be last in the order of solution,
logical :: ztran = .false.      ! T perform the transformation of the solution coefficients
                                ! to the expression for the wavefunction amplitudes at the grid
                                ! points. store the data on stream iwave.
                                ! ztran = T automatically sets zvec = t for idia > -2.
logical :: ztwod = .false.      ! T perform 2D calculation only at specified grid point.
logical :: zbisc                ! T place the Z-axis along the bisector
logical :: zperp = .false.      ! T place the Z-axis perpendicular to the molecule place
logical :: zpfun = .false.      ! F store energy levels on stream ilev

logical :: zpseg = .false.
logical :: zs0   = .false.
logical :: zs1   = .false.
logical :: zx    = .false.

!Integer variable
! /prt/ /outp/
integer :: idiag1 = 20          ! the final Hamiltonian matrix is written on units IDIAG1 and IDIAG2.
integer :: idiag2 = 21          ! the final Hamiltonian matrix is written on units IDIAG1 and IDIAG2.
integer :: iout1  = 24          ! stream for arrays iv1 and iv2, which record the sizes of
                                ! the truncated vctors. used when zvec = t.
integer :: iout2  = 25          ! stream for the 1d, 2d and 3d vectors for use when zvec = t.
integer :: iwave  = 26          ! stores the wavefunction amplitudes at the grid points when
                                ! ztran = t.
integer :: ilev   = 14          ! stream for final eigenvalues (formatted).
integer :: ieigs1 = 7           ! stream for eigenvalues of the 1d solutions.
integer :: ieigs2 = 2           ! stream for eigenvalues of the 2d solutions.
integer :: ivecs1 = 3           ! stream for eigenvectors of the 1d solutions.
integer :: ivecs2 = 4           ! stream for eigenvectors of the 2d solutions.
integer :: ivint  = 17          ! a scratch stream used for storing intermediate vectors in
                                ! building the final hamiltonian.
integer :: iband  = 15          ! scratch file used for storing bands of the final hamiltonian.
integer :: intvec = 16          ! a scratch stream for intermediate storage of the 2d vectors.
integer :: itime0

namelist /prt/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,&
                zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zperp,&
                idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,&
                ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,&
                zpseg

!Double precision real
!Common /oupb/
real(kind=dp) :: xp0
real(kind=dp) :: xp1
real(kind=dp) :: xp2

end module dvr3drjz_file
!unmention variable
!integer :: ncoord