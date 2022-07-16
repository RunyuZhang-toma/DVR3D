Module Constant

!==================================================================
!Section 000
!Copyright(C), 2022, UCL
!FileName: Constant.f90
!Author: Runyu Zhang
!Data: 04/07/2022
!Description: This is a test mudole file for dipole3.f90


!==================================================================
! include file: dvr3drjz.f90 dipj0dvr.f90 dipole3_with_rme.f90
!==================================================================================================

module dvr3drjz_file
! All common data from file dvr3drjz.f90
    save

    implicit none
!/ from file dvr3drjz.f90
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

logical :: zdcore = .false.     ! T for in core diagonalisation
logical :: z1da = .false.
logical :: zptra = .false.      ! print the transformed vectors.
logical :: zplot = .false.

! /diffs/
logical :: zsame

! /logic/
logical :: zembed
logical :: zdone
logical :: zncor = .false.
logical :: znco1
logical :: znco2
logical :: zprint = .false.     ! T supplies extra print out for debugging purposes.
logical :: zpmin                ! T supplies less  print out for large runs.
logical :: ztra = .true.        ! T writes out the data needed for program spectra 
                                ! to calculate simulated spectra.
logical :: zstart = .false.     ! T if we are writing out for spectra for the first time.
logical :: zrme1 = .true.       ! F program calculates reduced matrix elements for the dipole order.
                                ! We follow defintion of Lamouroux et al.
                                ! http://dx.doi.org/10.1016/j.jqsrt.2014.06.011
logical :: zrme2 = .true.       ! F program calculates reduced matrix elements for the quadrupole order.
logical :: zrme3 = .false.      ! F program calculates reduced matrix elements for the octupole order.
logical :: zout                 ! this should be true if the sorted line strengths are to be written
                                ! to the lineprinter. 
                                ! zout is set to true automatically if zspe is false.
logical :: zsort                ! if false subroutine sortsp is skipped.
logical :: zspe                 ! if false the program stops after sortsp. units of ispe are atomic units.
logical :: zpfun                ! calculates the partition function from energy levels supplied from
                                ! DVR3DRJZ and ROTLEV3/3B.
                                ! if zpfun false, the partition function
                                ! is set to q read in below.
logical :: zembed
logical :: zfit = .false.
logical :: zform = .true.
logical :: 


!==================================================================================================
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
integer :: jscr = 7
integer :: jvec = 3
integer :: jvec2 = 2
integer :: kvec = 8
integer :: kvec2 = 9
integer :: irf1 = 21
integer :: irf2 = 22
integer :: ivec = 26            ! input  eigenvalues & eigenvectors 
integer :: ivec1 = 27           ! input  eigenvalues & eigenvectors 
integer :: ivec2                ! input  eigenvalues & eigenvectors
integer :: ivec3                ! input  eigenvalues & eigenvectors
integer :: kvecpb = 9
integer :: nploti = 1
integer :: nplotf = 0
integer :: ithre = -8
integer :: idiag = 2


! /size/
! contral parameter from the problem
integer :: npnt         ! max(npnt1,npnt2) number of gauss-associated legendre grid points requested
integer :: npnt1        ! number of (gauss-laguerre) dvr points in r1
integer :: npnt2        ! number of (gauss-laguerre) dvr points in r2
integer :: npnta        ! the number of dvr points in
                        ! the coordinate to be treated first in the dvr successive
                        ! diagonalisation-truncation procedure
integer :: npntb        ! the number of dvr points in the coordinate to come second
integer :: npntc        ! the number of dvr points in the coordinate to come last
integer :: npntt
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
integer :: lmax
integer :: jk
integer :: ifile
integer :: nbass        ! maximum dimension of rotational secular problem
integer :: mbass
integer :: ibass        ! actual dimension of rotational secular problem
integer :: nlim
integer :: kmin         ! kmin=1 for sym. rotational basis, =0 for anti-sym.
                        ! kmin>1 loop over both.
integer :: neval2       ! neval for f block when kmin>1.
integer :: meval        ! number of eigenvalues computed in the vibrational problem
integer :: keval        ! number of eigenvectors used for iterations (=neval+4)
integer :: nvib         ! number of vibrational eigenvalues used in rotational prob.
integer :: nblk         ! number of k values
integer :: loff         ! space required for all the off-diagonal blocks
integer :: loff0        ! space required for the largest off-diagonal block
integer :: kbass
integer :: nmax         ! number of dvr points in each radial coordinate
integer :: maxblk       ! size of vibrational radial problem (even basis)
integer :: mxblk2       ! size of vibrational radial problem (odd  basis)
integer :: maxblk_even
integer :: maxblk_odd
integer :: kmax
integer :: ndvr         ! maximum dimension of theta dvr grid used in vibrational problem
integer :: iang         ! maximum number of discrete angles retained in vib. problem
integer :: nktot        ! number of k values
integer :: kpar
integer :: iqpar
integer :: nr           ! number of dvr points in each radial coordinate

! /sizes/
integer :: lbra0
integer :: lbra1
integer :: nbra0
integer :: nbra1
integer :: lket0
integer :: lket1
integer :: nket0
integer :: nket1
integer :: ntheta
integer :: nr1
integer :: nr2
integer :: neval0
integer :: neval1

! /diffs/
integer :: nqe
integer :: nqo
integer :: nr21
! /ligic/
integer :: iptot
integer :: idia
! /stream/
integer :: ibra = 12    ! input stream for the bra
integer :: ibra0
integer :: ibra1
integer :: itra = 13    ! output stream to program spectrm (if ztra).
                        ! note that for all times other than the dipole assumes 
                        ! that we have accessed the permanent dataset or file which has the
                        ! data from previous runs and that we are writing to the end of that file.
                        ! ************************************************
                        ! **  for the sake of safety you are therefore  **
                        ! **  advised to keep one previous edition as   **
                        ! **  backup!                                   **
                        ! ************************************************
integer :: iscr = 24    ! scratch file used for restart runs
integer :: ires = 0
integer :: iket = 11    ! input stream for the ket.
integer :: iket0
integer :: iket1
integer :: iwave0
integer :: iwave1
integer :: ivc0
integer :: ivc1
integer :: ione
integer :: itwo
integer :: mblock
integer :: nblock = 1000
!        e) outzrme1 (14) = program calculates reduced matrix elements for the dipole order.
!                      We follow defintion of Lamouroux et al. http://dx.doi.org/10.1016/j.jqsrt.2014.06.011
!        f) outzrme2 (15) = program calculates reduced matrix elements for the quadrupole order.
!                      We follow defintion of Lamouroux et al. http://dx.doi.org/10.1016/j.jqsrt.2014.06.011
!        fg outzrme3 (16) = program calculates reduced matrix elements for the octupole order.
!                      We follow defintion of Lamouroux et al. http://dx.doi.org/10.1016/j.jqsrt.2014.06.011
!
!        a) ires (0) restart parameter
!           ires = 0, normal run
!           ires = 1, restart run
!        b) nblock (1000) number of k --> k' blocks to be attempted
!
!       l3) npot, nv1, nv2 (all in i5 format)
!        a) npot  = number of Gauss-Legendre integration points used
!        b) nv1   = number of ket eigenfunctions considered.
!                   if this is input as zero, all available
!                   ket eigenfunctions will be considered when
!                   computing transitions.
!        c) nv2   = as above for the bra.
! /dim/
integer :: nrade
integer :: nrado
integer :: npot         ! number of Gauss-Legendre integration points used
integer :: nbin         ! largest binomial coef. required for angular integration(+1)
integer :: nbmax1
integer :: nbmax2
integer :: mbass        ! maximum size of vibrational problem (excluding linear geom)
integer :: mbass0       ! maximum size of vibrational problem (including linear geom)
integer :: mbass1
integer :: mbass2
integer :: kmin1
integer :: kmin2
integer :: jk1
integer :: jk2
integer :: nn2
integer :: ibase1       ! number of lowest ket eigenfunctions skipped
integer :: ibase2       ! number of lowest bra eigenfunctions skipped
integer :: ipot
integer :: nv1          ! number of bra eigenfunctions considered
integer :: nv2          ! number of ket eigenfunctions considered
integer :: nprt
integer :: npropin
integer :: lpot
! /sym/
integer :: ipar1
integer :: ipar2
integer :: jrot1
integer :: jrot2

! /pb/
integer :: inda1 = 100
integer :: inda2 = 100
integer :: indb1 = 100
integer :: indb2 = 100
integer :: indk = 100
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
integer :: mevlab
integer :: ibassa
integer :: ibassb
integer :: nviba
integer :: nvibb

! /pot/
integer :: iprop        ! vector with the information on which properties
                        ! will be considered in the run




namelist /prt/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,&
                zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zperp,&
                idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,&
                ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,&
                zpseg

!     outp holds information which controls the amount of printed output
!     toler: convergence tolerance for the iterative diagonaliser
!            toler = 0.0 gives machine accuracy
!     zpham: print matrix hamil if zpham = .true.
!     zpvec: print eigenvectors if zpvec = .true.
!     thresh: threshold for printing a coefficient if zpvec=.true.
!     zptra: print the transformed vectors.
!     zdcore: = .true. for in core diagonalisation
!     zdiag:  = .false. do not diagonalise the Hamiltonian matrix.
!     ztran:  = .true. transform eigenvectors back to original basis.
!     zvec:   = .true. eigenvalues and eigenvectors written to disk file.
!     zpfun:  = .true.  eigenvalues concatenated on stream ILEV.
!     stream         holds                              used if
!      ilev    input/output of eigenvalues              zpfun=.true.
!      ivec    input  eigenvalues & eigenvectors        always
!      ivec2   input  eigenvalues & eigenvectors        nblk .gt. 2
!      jvec    output first  set eigenvalues & vectors  zvec=.true.
!      jvec2   output second set                        zvec=.true.
!      kvec    output first  set transformed vectors    ztran=.true.
!      kvec2   output second set                        ztran=.true.
!      iscr    hamiltonian file                          always
!      irf1    restart file one                         zdiag=.false.
!      irf2    restart file two                         always
!
!     ires = 0  normal run
!          = 1  restart from first  call to dgrot
!          = 2  restart from second call to dgrot
!          = 3  restart from first  call to dgrot, one diagonalisation only
!          = -1 perform both transformations
!          = -2 perform second transformation only
!          = -3 perform first  transformation only
! (restart after zdiag=.false. run, ivec=irf1 and irf2 required)


!==================================================================================================
!Double precision real
! /oupb/
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
real(kind=dp) :: ezero

! /time/
real(kind=dp) :: ouser
real(kind=dp) :: osys
real(kind=dp) :: ototal

! /diffs/
real(kind=dp) :: alphae
real(kind=dp) :: betae
real(kind=dp) :: alphao
real(kind=dp) :: betao

! /eqm/
real(kind=dp) :: ex(3)
real(kind=dp) :: ez(3)
real(kind=dp) :: tmass

! /outp/
real(kind=dp) :: toler = 0.0D0
real(kind=dp) :: thresh = 0.1D0 ! threshold for printing a coefficient if zpvec=.true.












! /head/
!===================================================
!real(kind=dp) :: title
!===================================================

















end module dvr3drjz_file
!unmention variable
!integer :: ncoord

!=====================================================================================================
module input
!  definition of the control input parameters
  save
  !line 1
  !namelist prt
  logical :: zprint = .false. ! =T supplies extra print out for debugging purposes
  logical :: ztra = .true.  ! =T writes data for spectra to stream ITRA
  logical :: zstart = .false. ! =T initiates the output file for the data for SPECTRA
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
module lists
!linked list definitions for the complicated header reading procedure
!this is just an exercise, feel free to get rid of this mess (but it works ;-) 
  save
  type rgridnode
     real(8) :: rgrid 
     type (rgridnode),pointer :: next
  end type rgridnode
  type eigsnode
     real(8) :: eigs  
     type (eigsnode),pointer :: next
  end type eigsnode
  
  type headnode
     integer :: n3d    ! total number of grid points 
     integer :: nr     ! number of radial grid points
     integer :: jrot   ! total rotational angular momentum
     integer :: nval   ! number of eigenvalues
     integer :: kpar   ! laboratory coordinates parity
     integer :: ifpar  ! nuclear permutation parity
     integer :: nk     ! total number of k blocks
     type(rgridnode),pointer :: rgridfirst,rgridlist 
     type(eigsnode),pointer :: eigsfirst,eigslist
     type (headnode),pointer :: next 
  end type headnode
  type(headnode),pointer :: headfirst,headlist 
end module lists

module workdata
  save
  integer :: nk1, nk2, ispar1, ispar2, kz1, kz2, nval, nval1, nval2, nth1, nth2
  integer :: jrot1,jrot2,iqpar1,iqpar2,ifpar1,ifpar2,kpar1,kpar2,n3d1,n3d2,nr1,nr2
  real(8), allocatable :: eigs1(:),eigs2(:)
  real(8), allocatable :: rgrid1(:), rgrid2(:), thgrid1(:), thgrid2(:)  
  
  ! common grid
  integer :: nr,nth
  real(8) :: x,x1,x2
  real(8), allocatable :: thgrid(:),tanth(:),rgrid(:),wt(:)
  real(8), allocatable :: waves1(:,:),waves2(:,:)
  real(8), allocatable :: dipcxa(:),dipcya(:),dipcxb(:),dipcyb(:)
  real(8), allocatable :: pol1(:,:),pol2(:,:)

  ! dipole storage
  real(8), allocatable :: dpba(:,:)!,dpbb(:,:)
  real(8) :: pi,sumpb
end module workdata