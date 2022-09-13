!Constants.f90
!==================================================================================================
!Copyright(C). 2022, University College London
!File name: model.f90
!Author: Runyu Zhang & Jonathan Tennyson
!Version: 1.1
!Data: 7th/july/2022
!Description: this Fortran90 file contains the module template for the source folder common 
!             constants divided by the common group
!Dependency: Folder source
!===================================================================================================

module sizes
implicit none
!constant count: integer 13.

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
    
end module sizes

!===================================================================================================

module diffs
implicit none
!constant count: integer 3, real 4, ligical 1.

    integer :: nqe
    integer :: nqo
    integer :: nr21

    real(kind=dp) :: alphae
    real(kind=dp) :: betae
    real(kind=dp) :: alphao
    real(kind=dp) :: betao

    logical zsame
end module diffs

!===================================================================================================

module logic
implicit none
!constant count: integer 6, ligical 22.

    integer :: iptot
    integer :: idia
    integer :: iket = 11      ! input stream for the ket.
    integer :: itra = 13      ! output stream to program spectrm (if ztra).
                              ! note that for all times other than the dipole assumes 
                              ! that we have accessed the permanent dataset or file which has the
                              ! data from previous runs and that we are writing to the end of that file.
                              ! ************************************************
                              ! **  for the sake of safety you are therefore  **
                              ! **  advised to keep one previous edition as   **
                              ! **  backup!                                   **
                              ! ************************************************

    integer :: itra0
    integer :: ilev = 14      ! stream for final eigenvalues (formatted).

    logical :: zembed = .true.     ! T z axis is along r2, = f z axis is along r1.
    logical :: zdone = .true.      ! T use morse oscillator-like functions for r_1 coordinate;
                                   ! F use spherical oscillator functions.
    logical :: zmors1 = .true.     ! T use morse oscillator-like functions for r_1 coordinate;
    logical :: znco1
    logical :: znco2
    logical :: zprint  = .false.   ! T supplies extra print out for debugging purposes.           
    logical :: zpmin               ! T supplies less  print out for large runs.
    logical :: ztra = .true.       ! T writes out the data needed for program spectra 
                                   ! to calculate simulated spectra.

    logical :: zstart = .false.    ! T if we are writing out for spectra for the first time.
    logical :: zmors2 = .true.     ! T use morse oscillator-like functions for r_2 coordinate;
                                   ! F use spherical oscillator functions.

    logical :: zrme1 = .true.      ! F program calculates reduced matrix elements for the dipole order.
                                   ! We follow defintion of Lamouroux et al.
                                   ! http://dx.doi.org/10.1016/j.jqsrt.2014.06.011

    logical :: zrme2 = .true.      ! F program calculates reduced matrix elements for the quadrupole order.
    logical :: zrme3 = .false.     ! F program calculates reduced matrix elements for the octupole order.
    logical :: zpseg = .false.
    logical :: zuvvis
    logical :: zout                ! this should be true if the sorted line strengths are to be written
    logical :: zsort               ! if false subroutine sortsp is skipped.
    logical :: zspe                ! if false the program stops after sortsp. units of ispe are atomic units.
    logical :: zpfun               ! calculates the partition function from energy levels supplied from
                                   ! DVR3DRJZ and ROTLEV3/3B.
                                   ! if zpfun false, the partition function
                                   ! is set to q read in below.

    logical :: zncor = .false.
    logical :: zfit = .false.
    logical :: zform = .true.

end module logic

!===================================================================================================

module base
implicit none
!constant count: integer 2.

    integer :: ibase1 ! number of lowest ket eigenfunctions skipped
    integer :: ibase2 ! number of lowest bra eigenfunctions skipped

end module base

!===================================================================================================

module stream
implicit none
!constant count: integer 17.
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
    integer :: iket = 11      ! input stream from DVR3DRJZ/ROTLEV3/ROTLEV3B for the ket (unformated)
    integer :: ibra = 12      ! input stream for the bra (unformmatted)
    integer :: itra = 13      ! output stream to program spectrm (if ztra).
                              ! note that for all times other than the dipole assumes 
                              ! that we have accessed the permanent dataset or file which has the
                              ! data from previous runs and that we are writing to the end of that file.
                              ! ************************************************
                              ! **  for the sake of safety you are therefore  **
                              ! **  advised to keep one previous edition as   **
                              ! **  backup!                                   **
                              ! ************************************************

    integer :: iscr = 24      ! scratch file used for restart runs
                              ! holds hamiltonian file used if always

    integer :: ires = 0       ! a) ires (0) restart parameter
                              ! ires = 0, normal run
                              ! ires = 1, restart run

    integer :: mblock
    integer :: nblock = 1000  ! b) nblock (1000) number of k --> k' blocks to be attempted

end module stream

!===================================================================================================

module mass
implicit none
!constant count: real 5, logical 2.

     real(kind=dp) :: xmass(3)
     real(kind=dp) :: xmassr(3)
     real(kind=dp) :: g1
     real(kind=dp) :: g2
     real(kind=dp) :: ezero

     logical :: zembed = .true.    ! T z axis is along r2, = f z axis is along r1.
                                   ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.

     logical :: zbisc              ! T place the Z-axis along the bisector

end module mass

!===================================================================================================

module eqm
implicit none
!constant count: real 3.

     real(kind=dp) :: ex(3)
     real(kind=dp) :: ez(3)
     real(kind=dp) :: tmass

end module eqm

!===================================================================================================

module old
implicit none
!constant count: integer 10， logical 4.
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

     logical :: ztheta = .true. ! T let theta be first in the order of solution;
                                ! F let theta be last in the order of solution.

     logical :: zr2r1 = .true.  ! T let r_2 come before r_1 in the order of solution;
                                ! F let r_1 come before r_2 in the order of solution. (only idia > -2).

     logical :: zthet1
     logical :: zr2r11

end module old

!===================================================================================================

module time
implicit none
!constant count: real 3.

     real(kind=dp) :: ouser
     real(kind=dp) :: osys
     real(kind=dp) :: ototal

end module time

!===================================================================================================

module head
implicit none
!constant count: real 1.

     real(kind=dp) :: title

end module head

!===================================================================================================

module timing
implicit none
!constant count: integer 1.

     integer itime0

end module timing

!===================================================================================================

module sym
implicit none
!constant count: integer 5.

     integer :: idia     ! 1 scattering coordinates heteronuclear diatomic
                         ! 2 scattering coordinates homonuclear diatomic
                         ! -1 radau  coordinates hetronuclear diatomic
                         ! -2 radau  coordinates homonuclear  diatomic
                         ! 0 radau   coordinates with the z axis perpendicular to the molecular plane.

     integer :: ipar1
     integer :: ipar2
     integer :: jrot1
     integer :: jrot2

end module sym

!===================================================================================================

module dim
implicit none
!constant count: integer 33.

     integer :: neval         ! number of eigenvalues which have to actually be supplied as output
     integer :: lpot
     integer :: ncoord        ! number of vibrational coordinates explicitly considered
                              ! ncoord = 2: atom-diatom problem with diatom rigid
                              ! ncoord=2: also need lmax,lpot,idia,kmin
                              ! ncoord = 3: full 3-d triatomic problem
                              ! ncoord=3: all paramters required

     integer :: npnt          ! max(npnt1,npnt2) number of gauss-associated legendre grid points requested
     integer :: npnt1         ! number of (gauss-laguerre) dvr points in r1
     integer :: npnt2         ! number of (gauss-laguerre) dvr points in r2
     integer :: nrade
     integer :: nrado
     integer :: npot          ! number of Gauss-Legendre integration points used
                              ! in i5 format

     integer :: nbin          ! largest binomial coef. required for angular integration(+1)
     integer :: nbmax1
     integer :: nbmax2
     integer :: mbass         ! maximum size of vibrational problem (excluding linear geom)
     integer :: mbass1
     integer :: mbass2
     integer :: kmin1
     integer :: kmin2
     integer :: jk1
     integer :: jk2
     integer :: neval1
     integer :: neval2
     integer :: nn2
     integer :: ibase1        ! number of lowest ket eigenfunctions skipped
     integer :: ibase2        ! number of lowest bra eigenfunctions skipped
     integer :: ipot
     integer :: nn2
     integer :: lmax
     integer :: npropin
     integer :: nprt
     integer :: jrot          ! total angular momentum of the molecule
     integer :: idia
     integer :: ipar          ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
     integer :: nv1           ! number of bra eigenfunctions considered
                              ! if this is input as zero, all available
                              ! ket eigenfunctions will be considered when computing transitions.
                              ! in i5 format

end module dim

!===================================================================================================

module outp
implicit none
!constant count: integer 30, real 2, logical 34.

     integer :: kvecpb = 9
     integer :: ivec1 = 27
     integer :: idiag = 2
     integer :: idiag1 = 20        ! the final Hamiltonian matrix is written on units IDIAG1 and IDIAG2.
     integer :: idiag2 = 21        ! the final Hamiltonian matrix is written on units IDIAG1 and IDIAG2.
     integer :: iout1 = 24         ! stream for arrays iv1 and iv2, which record the sizes of
                                   ! the truncated vctors. used when zvec = t.

     integer :: iout2 = 25         ! stream for the 1d, 2d and 3d vectors for use when zvec = t.
     integer :: iwave = 26         ! stores the wavefunction amplitudes at the grid points when
     integer :: ilev = 14          ! stream for final eigenvalues (formatted).
                                   ! holds input/output of eigenvalues used if zpfun = .true.

     integer :: ieigs1 = 7         ! stream for eigenvalues of the 1d solutions.
     integer :: ieigs2 = 2         ! stream for eigenvalues of the 2d solutions.
     integer :: ivecs1 = 3         ! stream for eigenvectors of the 1d solutions.
     integer :: ivecs2 = 4         ! stream for eigenvectors of the 2d solutions.
     integer :: ivint = 17         ! a scratch stream used for storing intermediate vectors in
                                   ! building the final hamiltonian.

     integer :: iband = 15         ! scratch file used for storing bands of the final hamiltonian.
     integer :: intvec = 16        ! a scratch stream for intermediate storage of the 2d vectors.
     integer :: jscr = 7
     integer :: jvec = 3           ! holds output first  set eigenvalues & vectors used if zvec=.true.
     integer :: jvec2 = 2          ! holds output second set                        zvec=.true.
     integer :: kvec = 8           ! holds output first  set transformed vectors used if ztran=.true.
     integer :: kvec2 = 9          ! holds output second set   used if ztran=.true.
     integer :: iscr = 24          ! scratch file used for restart runs
                                   ! holds hamiltonian file used if always

     integer :: ires = 0           ! = 0  normal run
                                   ! = 1  restart from first  call to dgrot
                                   ! = 2  restart from second call to dgrot
                                   ! = 3  restart from first  call to dgrot, one diagonalisation only
                                   ! = -1 perform both transformations
                                   ! = -2 perform second transformation only
                                   ! = -3 perform first  transformation only
                                   ! (restart after zdiag=.false. run, ivec=irf1 and irf2 required)

     integer :: irf1 = 21          ! irf1   restart file one  used if zdiag=.false.
     integer :: irf2 = 22          ! irf2   restart file two  used if always
     integer :: ivec = 26          ! holds input  eigenvalues & eigenvectors used if always
     integer :: ivec2              ! holds input  eigenvalues & eigenvectors used if nblk > 2
     integer :: nploti = 1
     integer :: nplotf = 0
     integer :: ithre = -8

     real(kind=dp) :: toler = 0.0D0     ! convergence tolerance for the iterative diagonaliser
                                        ! toler = 0.0 gives machine accuracy

     real(kind=dp) :: thresh = 0.1D0    ! threshold for printing a coefficient if zpvec=.true.

     logical :: zpham = .false.    ! T request printing of the hamiltonian matrix
     logical :: zplot = .false.
     logical :: zprad = .false.    ! T request printing of the radial matrix elements
     logical :: zpvec = .false.    ! T request printing of the eigenvectors
     logical :: zrot = .true.      ! F do vibrational part of rotational calculation by looping over k
     logical :: zladd = .true.     ! T NALF kept constant as k increases
                                   ! F NALF decreases with k (=f has a bug), (only if zrot = .true.)

     logical :: zembed = .true.    ! T z axis is along r2, = f z axis is along r1.
                                   ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.

     logical :: zmors2 = .true.    ! T use morse oscillator-like functions for r_2 coordinate;
                                   ! F use spherical oscillator functions.

     logical :: zs0 = .false.
     logical :: zx = .false.
     logical :: zs1 = .false.
     logical :: zpmin = .false.    ! T requests only minimal printing.
     logical :: zvec = .false.     ! T store the eigenvectors from all the parts of the calculation
                                   ! eigenvalues and eigenvectors written to disk file.
                                   ! (1d,2d and 3d) on stream iout2.
                                   ! further information relating to this (arrays iv1 and iv2) is
                                   ! stored on stream iout1.

     logical :: zquad2 = .true.    ! T use the dvr quadrature approximation for the integrals of
                                   ! the r_2^{-2} matrix, and hence make its dvr transformation diagonal.
                                   ! F evaluate the r_2^{-2} integrals fully and perform the
                                   ! full dvr transformation on them.
                                   ! Note that zquad2 = f is not implemented for zmors2 = t
                                   ! or for ztheta = f.

     logical :: zdiag = .true.     ! F do not do final diagonalisation, instead the final Hamiltonian
                                   ! matrix is written on units IDIAG1 and IDIAG2. 

     logical :: zlmat = .false.    ! T requests printing of the L-matrix.
     logical :: zcut = .false.     ! T final dimension selected using an energy cut-off given  by emax2.
                                   ! F final dimension determined by nham3.

     logical :: zall = .false.     ! T requests no truncation of the intermediate solution.
     logical :: zlin = .false.     ! T NALF kept constant as k increases
                                   ! F NALF decreases with k (=f has a bug), (only if zrot = .true.)

     logical :: zp1d = .false.     ! T requests printing of the results of 1d calculations.
     logical :: zp2d = .false.     ! T requests printing of the results of 2d calculations.
     logical :: zr2r1 = .true.     ! T let r_2 come before r_1 in the order of solution;
                                   ! F let r_1 come before r_2 in the order of solution. (only idia > -2).

     logical :: ztheta = .true.    ! T let theta be first in the order of solution;
                                   ! F let theta be last in the order of solution,

     logical :: ztran = .false.    ! T perform the transformation of the solution coefficients
                                   ! to the expression for the wavefunction amplitudes at the grid
                                   ! points. store the data on stream iwave.
                                   ! ztran = T automatically sets zvec = t for idia > -2.

     logical :: zmors1 = .true.    ! T use morse oscillator-like functions for r_1 coordinate;
                                   ! F use spherical oscillator functions.

     logical :: ztwod = .false.    ! T perform 2D calculation only at specified grid point.
     logical :: zbisc              ! T place the Z-axis along the bisector
     logical :: zperp = .false.    ! T place the Z-axis perpendicular to the molecule place
     logical :: zpseg = .false.
     logical :: zpfun = .false.    ! F store energy levels on stream ilev
     logical :: zembed = .true.    ! T z axis is along r2, = f z axis is along r1.
                                   ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.

     logical :: zptra = .false.    ! print the transformed vectors.
     logical :: zdcore = .false.   ! T for in core diagonalisation
     logical :: z1da = .false.
   

end module outp

!===================================================================================================

module size
implicit none
!constant count: integer 69, real 2.

     integer :: npnt1    ! number of (gauss-laguerre) dvr points in r1
     integer :: npnt2    ! number of (gauss-laguerre) dvr points in r2
     integer :: nktot    ! number of k values
     integer :: kpar
     integer :: iqpar
     integer :: maxblk_odd
     integer :: nalf     ! number of (gauss-legendre) dvr points in theta
     integer :: nmax1    ! max order of r1 radial laguerre polynomial ( = npnt1-1)
     integer :: nmax2    ! max order of r2 radial laguerre polynomial ( = npnt2-1)
     integer :: nmax     ! number of dvr points in each radial coordinate
     integer :: maxblk   ! size of vibrational radial problem (even basis)
     integer :: maxleg   ! max order of angular legendre polynomial   ( = nalf -1)
     integer :: nalf2
     integer :: idvr     ! number of unique dvr points
     integer :: npnt     ! max(npnt1,npnt2) number of gauss-associated legendre grid points requested
     integer :: nlim1    ! nmax1+1*(nmax1+1+1)/2
     integer :: nlim2    ! nmax2+1*(nmax2+1+1)/2
     integer :: neval    ! number of eigenvalues which have to actually be supplied as output
     integer :: ncoord   ! number of vibrational coordinates explicitly considered
                         ! ncoord = 2: atom-diatom problem with diatom rigid
                         ! ncoord=2: also need lmax,lpot,idia,kmin
                         ! ncoord = 3: full 3-d triatomic problem
                         ! ncoord=3: all paramters required

     integer :: jrot     ! total angular momentum of the molecule
     integer :: kmin     ! zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
                         ! kmin=2 loop over both sym & anti-sym (zbisc=t only)
                         ! zrot=f, kmin=fixed value of k

     integer :: idia     ! 1 scattering coordinates heteronuclear diatomic
                         ! 2 scattering coordinates homonuclear diatomic
                         ! -1 radau  coordinates hetronuclear diatomic
                         ! -2 radau  coordinates homonuclear  diatomic
                         ! 0 radau   coordinates with the z axis perpendicular to the molecular plane.

     integer :: ipar     ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
     integer :: max2d    ! upper bound on size of intermediate 2d hamiltonian
     integer :: max3d    ! upper bound on size of full 3d hamiltonian
     integer :: max2d2   ! max2d for smaller block (zbisc=t only)
     integer :: max3d2   ! max3d for smaller block (zbisc=t only)
     integer :: npnta    ! the number of dvr points in
                         ! the coordinate to be treated first in the dvr successive
                         ! diagonalisation-truncation procedure

     integer :: npntb    ! the number of dvr points in the coordinate to come second
     integer :: npntc    ! the number of dvr points in the coordinate to come last
     integer :: ndima    ! set equal to npnta at the start - used for dimensioning
     integer :: ndimb    ! set equal to npntb at the start - used for dimensioning
     integer :: ndimc    ! set equal to npntc at the start - used for dimensioning
     integer :: iq
     integer :: idia     ! 1 scattering coordinates heteronuclear diatomic
                         ! 2 scattering coordinates homonuclear diatomic
                         ! -1 radau  coordinates hetronuclear diatomic
                         ! -2 radau  coordinates homonuclear  diatomic
                         ! 0 radau   coordinates with the z axis perpendicular to the molecular plane.
     integer :: ipar     ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
     integer :: lmax
     integer :: npnt1    ! number of (gauss-laguerre) dvr points in r1
     integer :: npnt2    ! number of (gauss-laguerre) dvr points in r2
     integer :: jrot     ! total angular momentum of the molecule
     integer :: kmin     ! zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
                         ! kmin=2 loop over both sym & anti-sym (zbisc=t only)
                         ! zrot=f, kmin=fixed value of k
     integer :: neval    ! number of eigenvalues which have to actually be supplied as output
     integer :: jk
     integer :: ifile
     integer :: NBASS    ! maximum dimension of rotational secular problem
     integer :: MBASS    ! maximum size of vibrational problem (excluding linear geom)
     integer :: IBASS    ! actual dimension of rotational secular problem
     integer :: NEVAL    ! number of eigenvalues which have to actually be supplied as output
     integer :: IPAR     ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
     integer :: IDIA     ! 1 scattering coordinates heteronuclear diatomic
                         ! 2 scattering coordinates homonuclear diatomic
                         ! -1 radau  coordinates hetronuclear diatomic
                         ! -2 radau  coordinates homonuclear  diatomic
                         ! 0 radau   coordinates with the z axis perpendicular to the molecular plane.

     integer :: nlim
     integer :: jrot     ! total rotational angular momentum
     integer :: KMIN     ! zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
                         ! kmin=2 loop over both sym & anti-sym (zbisc=t only)
                         ! zrot=f, kmin=fixed value of k
                         
     integer :: NEVAL2   ! neval for f block when kmin>1.
     integer :: MEVAL    ! number of eigenvalues computed in the vibrational problem
     integer :: KEVAL    ! number of eigenvectors used for iterations (=neval+4)
     integer :: NVIB     ! number of vibrational eigenvalues used in rotational prob.
     integer :: NBLK     ! number of k values
     integer :: LOFF     ! space required for all the off-diagonal blocks
     integer :: LOFF0    ! space required for the largest off-diagonal block
     integer :: kbass
     integer :: npnt1    ! number of (gauss-laguerre) dvr points in r1
     integer :: npnt2    ! number of (gauss-laguerre) dvr points in r2
     integer :: npntt
     integer :: kmax
     integer :: ndvr     ! maximum dimension of theta dvr grid used in vibrational problem
     integer :: iang     ! maximum number of discrete angles retained in vib. problem
     integer :: mxblk2   ! size of vibrational radial problem (odd  basis)
     integer :: mbass0   ! maximum size of vibrational problem (excluding linear geom)

     real(kind=dp) :: emax1
     real(kind=dp) :: emax2
     
end module size

!===================================================================================================

module split1
implicit none
!constant count: integer 1, real 7.
     integer :: iu1

     real(kind=dp) :: re1
     real(kind=dp) :: diss1
     real(kind=dp) :: we1
     real(kind=dp) :: beta1
     real(kind=dp) :: ur1
     real(kind=dp) :: urr1
     real(kind=dp) :: a1

end module split1

!===================================================================================================

module split2
implicit none
!constant count: integer 1, real 7.
     integer :: iu2

     real(kind=dp) :: re2
     real(kind=dp) :: diss2
     real(kind=dp) :: we2
     real(kind=dp) :: beta2
     real(kind=dp) :: ur2
     real(kind=dp) :: urr2
     real(kind=dp) :: a2

end module split2

!===================================================================================================

module oupb
implicit none
!constant count: real 3.

     real(kind=dp) :: xp0
     real(kind=dp) :: xp1
     real(kind=dp) :: xp2

end module oupb

!===================================================================================================

module pot
implicit none
     parameter (mxprop=1000)
     integer :: iprop(mxprop)
     ! dimension  iprop(mxprop)
end module pot

!===================================================================================================

module pb
implicit none
!constant count: integer 25
     integer :: inda1(100)
     integer :: inda2(100)
     integer :: indb1(100)
     integer :: indb2(100)
     integer :: indk(100)
     
     integer :: iqa
     integer :: iqb
     integer :: isa
     integer :: isb
     integer :: ipa
     integer :: ipb
     integer :: kmina
     integer :: kminb
     integer :: nka
     integer :: nkb
     integer :: nbassa
     integer :: nbassb
     integer :: nskipka
     integer :: nskipkb
     integer :: mevala
     integer :: mevalb
     integer :: ibassa
     integer :: ibassb
     integer :: nviba
     integer :: nvibb

end module pb

!===================================================================================================