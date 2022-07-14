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

! /size/


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

! /size/
! contral parameter from the problem
integer :: npnt         ! max(npnt1,npnt2)
integer :: npnt1        ! number of (gauss-laguerre) dvr points in r1
integer :: npnt2        ! number of (gauss-laguerre) dvr points in r2
integer :: npnta        ! the number of dvr points in
                        ! the coordinate to be treated first in the dvr successive
                        ! diagonalisation-truncation procedure
integer :: npntb        ! the number of dvr points in the coordinate to come second
integer :: npntc        ! the number of dvr points in the coordinate to come last
integer :: nalf         ! number of (gauss-legendre) dvr points in theta
integer :: nalf2
integer :: nmax1        ! max order of r1 radial laguerre polynomial ( = npnt1-1)
integer :: nmax2        ! max order of r2 radial laguerre polynomial ( = npnt2-1)
integer :: maxleg       ! max order of angular legendre polynomial   ( = nalf -1)
integer :: idvr         ! number of unique dvr points
integer :: nlim1        ! nmax1+1*(nmax1+1+1)/2
integer :: nlim2        ! nmax2+1*(nmax2+1+1)/2
integer :: neval        ! number of eigenvalues which have to actually be supplied as output
integer :: ncoord       ! number of vibrational coordinates explicitly considered
                        ! ncoord = 2: atom-diatom problem with diatom rigid
                        ! ncoord=2: also need lmax,lpot,idia,kmin
                        ! ncoord = 3: full 3-d triatomic problem
                        ! ncoord=3: all paramters required
integer :: jrot         ! total angular momentum of the molecule
integer :: kmin         ! zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
                        ! kmin=2 loop over both sym & anti-sym (zbisc=t only)
                        ! zrot=f, kmin=fixed value of k
integer :: idia         ! 1 scattering coordinates heteronuclear diatomic
                        ! 2 scattering coordinates homonuclear diatomic
                        ! -1 radau  coordinates hetronuclear diatomic
                        ! -2 radau  coordinates homonuclear  diatomic
                        ! 0 radau   coordinates with the z axis perpendicular to the molecular plane.
integer :: ipar         ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
integer :: max2d        ! upper bound on size of intermediate 2d hamiltonian
integer :: max2d2       ! max2d for smaller block (zbisc=t only)
integer :: max3d        ! upper bound on size of full 3d hamiltonian
integer :: max3d2       ! max3d for smaller block (zbisc=t only)
integer :: ndima        ! set equal to npnta at the start - used for dimensioning
integer :: ndimb        ! set equal to npntb at the start - used for dimensioning
integer :: ndimc        ! set equal to npntc at the start - used for dimensioning
integer :: iq


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

! /size/
real(kind=dp) :: emax1
real(kind=dp) :: emax2

! /split1/
real(kind=dp) :: re1
real(kind=dp) :: diss1
real(kind=dp) :: we1
real(kind=dp) :: beta1
real(kind=dp) :: ur1
real(kind=dp) :: urr1
real(kind=dp) :: a1
real(kind=dp) :: iu1

! /split2/
real(kind=dp) :: re2
real(kind=dp) :: diss2
real(kind=dp) :: we2
real(kind=dp) :: beta2
real(kind=dp) :: ur2
real(kind=dp) :: urr2
real(kind=dp) :: a2
real(kind=dp) :: iu2

! /mass/
real(kind=dp) :: xmass(3)
real(kind=dp) :: xmassr(3)
real(kind=dp) :: g1
real(kind=dp) :: g2

end module dvr3drjz_file
!unmention variable
!integer :: ncoord