!==================================================================================================
!Module defintion
!By Runyu Zhang & Tennyson Jonathan 1/Sep/2022
!Contains size, outp, oupb, timing, split1, split2, mass
!Special notice:: The module name contains the file name, it means cannot directly paste to other 
!                 files and use. There are some difference between the constant numbers, types and 
!                 default value.
!==================================================================================================
module dvr3drjz_size
    integer :: npnt     ! max(npnt1,npnt2) number of gauss-associated legendre grid points requested
    integer :: npnt1    ! number of (gauss-laguerre) dvr points in r1
    integer :: npnt2    ! number of (gauss-laguerre) dvr points in r2
    integer :: nalf     ! number of (gauss-legendre) dvr points in theta
    integer :: nalf2
    integer :: nmax1    ! max order of r1 radial laguerre polynomial ( = npnt1-1)
    integer :: nmax2    ! max order of r2 radial laguerre polynomial ( = npnt2-1)
    integer :: maxleg   ! max order of angular legendre polynomial   ( = nalf -1)
    integer :: idvr     ! number of unique dvr points
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

    double precision :: emax1
    double precision :: emax2

end module dvr3drjz_size

module dvr3drjz_outp
    logical :: zpham = .false.  !  zpham[f] = t requests printing of the hamiltonian matrix.
    logical :: zprad = .false.  !  zprad[f] = t requests printing of the radial matrix elements.
    logical :: zpvec = .false.  !  requests printing of the eigenvectors.
    logical :: zrot = .true.    ! F do vibrational part of rotational calculation by looping over k
    logical :: zladd = .true.   ! T NALF kept constant as k increases
                                ! F NALF decreases with k (=f has a bug), (only if zrot = .true.)
    logical :: zembed = .true.  ! T z axis is along r2, = f z axis is along r1.
                                ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.

    logical :: zmors2 = .true.  ! T use morse oscillator-like functions for r_2 coordinate;
                                ! F use spherical oscillator functions.
    logical :: zs0 = .false.
    logical :: zx = .false.
    logical :: zs1 = .false.
    logical :: zpmin = .false.  !  zpmin[f] = t requests only minimal printing.
    logical :: zvec = .false.   ! T store the eigenvectors from all the parts of the calculation
                                ! eigenvalues and eigenvectors written to disk file.
                                ! (1d,2d and 3d) on stream iout2.
                                ! further information relating to this (arrays iv1 and iv2) is
                                ! stored on stream iout1.
    logical :: zquad2 = .true.
    logical :: zdiag = .true.
    logical :: zlmat = .false.  !  zlmat[f] = t requests printing of the L-matrix.
    logical :: zcut = .false.   !  zcut[f]  = t final dimension selected using an energy cut-off given
                                !             by emax2.
                                !            = f final dimension determined by nham3.
    logical :: zall = .false.
    logical :: zlin = .false.
    logical :: zp1d = .false.   !  zp1d [f] = t requests printing of the results of 1d calculations.
    logical :: zp2d = .false.   !  zp2d [f] = t requests printing of the results of 2d calculations.
    logical :: zr2r1 = .true.   ! T let r_2 come before r_1 in the order of solution;
                                ! F let r_1 come before r_2 in the order of solution. (only idia > -2).
    logical :: ztheta = .true.  ! T let theta be first in the order of solution;
                                ! F let theta be last in the order of solution,
    logical :: ztran = .false.  ! T perform the transformation of the solution coefficients
                                ! to the expression for the wavefunction amplitudes at the grid
                                ! points. store the data on stream iwave.
                                ! ztran = T automatically sets zvec = t for idia > -2.
    logical :: zmors1 = .true.  ! T use morse oscillator-like functions for r_1 coordinate;
                                ! F use spherical oscillator functions.
    logical :: ztwod = .false.  ! T perform 2D calculation only at specified grid point.
    logical :: zbisc            ! T place the Z-axis along the bisector
    logical :: zperp = .false.  ! T place the Z-axis perpendicular to the molecule place
    logical :: zpseg = .false.  
    logical :: zpfun = .false.  ! F store energy levels on stream ilev

    integer :: idiag1 = 20      ! the final Hamiltonian matrix is written on units IDIAG1 and IDIAG2.
    integer :: idiag2 = 21      ! the final Hamiltonian matrix is written on units IDIAG1 and IDIAG2.
    integer :: iout1 = 24       ! stream for arrays iv1 and iv2, which record the sizes of
                                ! the truncated vctors. used when zvec = t.
    integer :: iout2 = 25       ! stream for the 1d, 2d and 3d vectors for use when zvec = t.
    integer :: iwave = 26       ! stores the wavefunction amplitudes at the grid points when
    integer :: ilev =14         ! stream for final eigenvalues (formatted).
                                ! holds input/output of eigenvalues used if zpfun = .true.

    integer :: ieigs1 = 7       ! stream for eigenvalues of the 1d solutions.
    integer :: ieigs2 = 2
    integer :: ivecs1 = 3       ! stream for eigenvectors of the 1d solutions.
    integer :: ivecs2 = 4       ! stream for eigenvectors of the 2d solutions.
    integer :: ivint = 17       ! a scratch stream used for storing intermediate vectors in
                                ! building the final hamiltonian.

    integer :: iband = 15       ! scratch file used for storing bands of the final hamiltonian.
    integer :: intvec = 16      ! a scratch stream for intermediate storage of the 2d vectors.
end module dvr3drjz_outp

module dvr3drjz_oupb
    double precision :: xp0
    double precision :: xp1
    double precision :: xp2
end module dvr3drjz_oupb

module dvr3drjz_timing
    integer :: itime0
end module dvr3drjz_timing

module dvr3drjz_split1
    integer :: iu1
    double precision :: re1
    double precision :: diss1
    double precision :: we1
    double precision :: beta1
    double precision :: ur1
    double precision :: urr1
    double precision :: a1
end module dvr3drjz_split1

module dvr3drjz_split2
    integer :: iu2
    double precision :: re2
    double precision :: diss2
    double precision :: we2
    double precision :: beta2
    double precision :: ur2
    double precision :: urr2
    double precision :: a2
end module dvr3drjz_split2

module dvr3drjz_mass
    double precision :: g1, g2
    double precision :: xmass(3), xmassr(3)
end module dvr3drjz_mass

program DVR3DRJZ
    call dvr3d
    stop
end

!######################################################################
    subroutine dvr3d
!
!     program               d v r 3 d r j z
!                           ~~~~~~~~~~~~~~~
!     should be cited as:
!         *********** add here ********************
!     program to do ro-vibrational calculations on triatomic systems
!     using general length, length, angle coordinates and a choice of
!     embeddings, in a multidimensional dvr in choice of coordinate
!     orders.
!     there are various options:
!     (a) calculation at a frozen angle                       ztwod  = t
!     (b) atom - rigid diatom calculations                    ncoord = 2
!     (c) triatomic calculations in all 3 dimensions          ncoord = 3
!     additionally, for Radau coordinates there are options to:
!     (d) place the z-axis along the bisector                 zbisc  = t
!     (e) place the z-axis perpendicular to the molecule plane zperp = t
!     see:
!
!     j.r.henderson & j.tennyson, chem.phys.lett., 173, 133 (1990).
!     j.r.henderson, phd thesis, university of london (1990).
!     j.r.henderson, j.tennyson & b.t. sutcliffe, j.chem.phys. 98, 7191 (1993).
!     j. tennyson & b.t. sutcliffe, int. j. quantum chem. 42, 941 (1992).
!
!     use as follows:
!     comments on namelist parameters (& defaults) in block data
!     input data read in subroutines insize & setcon
!     the program needs the following subroutines:
!            subroutine potv(v,r1,r2,x) should return the potential v
!            in hartrees for the point x = cos(theta) and
!            bondlengths r1 & r2 in bohr.
!     2. dsyev to do in core diagonalisation (lapack f90 routine).
!     the programme works in **** atomic units ***** :
!     1. the subroutine potv should return the potential in
!        hartrees for distances in bohr.
!     2. all input in setcon is in hartree or bohr except
!     3. the nuclear masses are read in atomic mass units & converted.
!     4. the eigenvalues are printed in both hartree & wavenumbers.
!
!     Rewritten into fortran 95 by Max Kostin and Jonathan Tennyson
    use dvr3drjz_outp
    use dvr3drjz_oupb
    use dvr3drjz_timing
    namelist/prt/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,&
                    zpseg
    implicit double precision (a-h,o-y), logical (z)

    write(6, '(5x"Program DVR3DRJZ (version of March 2002)")') 

!     read in namelist input data (defaults in block data)
    read(5,prt)

    call SYSTEM_CLOCK(itime0,irate2,imax2)

!     read in control parameters of problem.
    call insize

!     now do the calculation
    call ccmain

    call SYSTEM_CLOCK(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,'(/i10," secs CPU time used"/)') itime
    stop
end

!======================================================================
subroutine insize

!     set up common /size/ & write control parameters of problem
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb
    implicit double precision (a-h,o-y), logical (z)

!     common /size/ stores control parameters for the problem
!     npnt1: number of (gauss-laguerre) dvr points in r1
!     npnt2: number of (gauss-laguerre) dvr points in r2
!     nmax1: max order of r1 radial laguerre polynomial ( = npnt1-1)
!     nmax2: max order of r2 radial laguerre polynomial ( = npnt2-1)
!     maxleg:max order of angular legendre polynomial   ( = nalf -1)
!     max2d :upper bound on size of intermediate 2d hamiltonian
!     max2d2:max2d for smaller block (zbisc=t only)
!     max3d :upper bound on size of full 3d hamiltonian
!     max3d2:max3d for smaller block (zbisc=t only)
!     nalf : number of (gauss-legendre) dvr points in theta
!     idvr : number of unique dvr points
!     jrot : jrot is total angular momentum of the molecule
!     kmin : zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
!                    kmin=2 loop over both sym & anti-sym (zbisc=t only)
!            zrot=f, kmin=fixed value of k
!     npnt : max(npnt1,npnt2)
!     idia : = 1 scattering coordinates heteronuclear diatomic
!            = 2 scattering coordinates homonuclear diatomic
!            =-1 radau      coordinates hetronuclear diatomic
!            =-2 radau      coordinates homonuclear  diatomic
!            = 0 radau      coordinates with the z axis perpendicular to
!                the molecular plane.
!     ipar : parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
!     nlim1: =nmax1+1*(nmax1+1+1)/2
!     nlim2: =nmax2+1*(nmax2+1+1)/2
!     npnta: the number of dvr points in
!            the coordinate to be treated first in the dvr successive
!            diagonalisation-truncation procedure
!     npntb: the number of dvr points in the coordinate to come second
!     npntc: the number of dvr points in the coordinate to come last
!     ndima: set equal to npnta at the start - used for dimensioning
!     ndimb: set equal to npntb at the start - used for dimensioning
!     ndimc: set equal to npntc at the start - used for dimensioning
!     neval: number of eigenvalues which have to actually be supplied
!            as output
!     ncoord: number of vibrational coordinates explicitly considered
!     if (ncoord /= 3) some of the above are dummies, see below.

    character(len=8) title(9)
!     read in control parameters of problem:

!     ncoord = 2: atom-diatom problem with diatom rigid
!     ncoord = 3: full 3-d triatomic problem
!       read(5,5) ncoord
!     5 format(15i5)
    read(5,'(15i5)') ncoord
!     other paramters: see comments above on /size/ for meaning
!     if ncoord=2: also need lmax,lpot,idia,kmin
!     if ncoord=3: all paramters required

    read(5,'(15i5)') npnt2,jrot,neval,nalf,max2d,max3d,idia,&
                kmin,npnt1,ipar,max3d2
      
    zbisc = .false.
    if (jrot == 0) then
        zembed = .true.
        kmin = 0
        zrot = .false.
    endif
    if (idia == -2) then
        ncoord=3
        npnt1=npnt2
        idvr=nalf
        zmors1=zmors2
        ztheta=.false.
        if (jrot /= 0) zbisc=.true.
!        Q symmetry
!        iq=0 -> q=0("+")
!        iq=1 -> q=1("-")
        if (zperp) then
            iq=kmin
            kmin=1
            if (jrot == 0) kmin=0
        endif
    else
        if (ztran) zvec=.true.
    endif

    if (ztwod) then
        ncoord = 3
        idia = 1
        nalf = 1
        nalf2= 1
        idvr = 1
        zrot=.false.
        ztheta=.false.
        max3d=max2d
        neval=min(max2d,neval)
        goto 887
    endif

!     the gauss-legendre integration points are only acually
!     calculated for the half-range:
    nalf2 = (nalf+1)/2
!     are we doing atom, atom-rigid diatom or the full problem?
    if (ncoord == 2) then
!       atom - rigid diatom case, set dummy /size/
        if (zr2r1) then
            write(6,'(/10x,"atom - rigid diatom vibrational analysis with:"/)')
            npnt1 = 1
            nmax1 = 0
        else
            write(6,'(/5x,"Fixed - r2 vibrational analysis with:"/)')
            npnt2 = 1
            nmax2 = 0
        endif
    else
!        full case: print data about extra radial basis
        ncoord  = 3
        write(6,'(/5x,"Full triatomic vibrational problem with"/ &
                /5x,i5,3x,"radial r1 dvr points used,")') npnt1
    endif
      
  887 continue
!     maximum order of polynomials;
    maxleg = nalf - 1
    nmax2  = npnt2- 1
    nmax1  = npnt1- 1

    if (ztheta) then
        npnta = nalf/max(1,idia)
        if (zr2r1) then
            npntb = npnt2
            npntc = npnt1
        else
            npntb = npnt1
            npntc = npnt2
        endif
    else
        if (zr2r1) then
            npnta = npnt2
            npntb = npnt1
        else
            npnta = npnt1
            npntb = npnt2
        endif
        npntc = nalf/max(1,idia)
    endif

    if (idia == -2) then
        if (zrot .and. (jrot+kmin)>1) then
            max2d  = npnt1*(npnt1+1)/2
            max2d2 = npnt1*(npnt1-1)/2
        else
            max2d = npnt1*(npnt1+1-(ipar * 2))/2
            max2d2 = 0
        endif
        if (.not. zall)  max3d=min(max3d,max2d*nalf)
        if (zall) max3d=max2d*nalf
        if (zrot) then
            if (max3d2 <= 0) max3d2=max3d
            max3d2=min(max3d,max2d2*nalf,max3d2)
        endif
    else
        if (zall) then
            max2d = npnta*npntb
            max3d = max2d*npntc
        else
            max2d=min(max2d,npnta*npntb)
            max3d=min(max3d,npnta*npntb*npntc)
        endif
    endif

    if (neval <= 0) neval = 10
    neval=min(max3d,neval)
    if (ztwod) write(6,'(/5x,i5,3x,"radial r1 dvr points used,")') npnt1
    if (ncoord == 3) &
            write(6, '(5x,i5,3x,"radial r",i1," dvr points used,",&
            /5x,i5,3x,"angular dvr points used, with",&
            /5x,i5,3x,"lowest eigenvectors chosen from",&
            /5x,i5,3x,"maximum dimension secular problem"/)') npnt2,2,nalf,neval,max3d
    if (ncoord == 2 .and. zr2r1)&
          write(6,'(5x,i5,3x,"radial r",i1," dvr points used,",&
            /5x,i5,3x,"angular dvr points used, with",&
            /5x,i5,3x,"lowest eigenvectors chosen from",&
            /5x,i5,3x,"maximum dimension secular problem"/)') npnt2,2,nalf,neval,max3d
    if(ncoord == 2 .and. .not. zr2r1)&
          write(6,'(5x,i5,3x,"radial r",i1," dvr points used,",&
            /5x,i5,3x,"angular dvr points used, with",&
            /5x,i5,3x,"lowest eigenvectors chosen from",&
            /5x,i5,3x,"maximum dimension secular problem"/)') npnt1,1,nalf,neval,max3d
    if(idia == 2 .and. zperp) then
        write(6,'(/5x,"STOP!!!  ZPERP should be .true. for IDIA=2")')
        stop
    endif      
    read(5,'(9a8)')   title
    write(6,'(5x,"Title: ",9a8/)') title
    if (ncoord == 3) then
        if (zmors1)  write(6,'(5x,"Morse oscillators used for r",i1," basis")') 1
        if (.not. zmors1) write(6,'(5x,"Spherical oscillators used for r",i1," basis")') 1
    endif
    if (ncoord == 2 .and. .not. zr2r1) then
        if (zmors1)  write(6,'(5x,"Morse oscillators used for r",i1," basis")') 1
        if (.not. zmors1) write(6,'(5x,"Morse oscillators used for r",i1," basis")') 1
    else
        if (zmors2) write(6,'(5x,"Morse oscillators used for r",i1," basis")') 2
        if (.not. zmors2) write(6,'(5x,"Morse oscillators used for r",i1," basis")') 2
        if (zquad2) then
            write(6,'(/5x,"Quadrature approximation used for the r2**(-2) terms"/)')
        else
            write(6,'(/5x,"Quadrature approximation abandoned for r2**(-2) terms"/)')
        endif
    endif
    if (zall) write(6,'(/5x,"All solutions from lower dimensions have been used")')
    if (ztheta) then
        if (zr2r1) write(6,'(5x,"Problem solved in the order: theta -> r2 -> r1")')
        if (.not. zr2r1) write(6,'(5x,"Problem solved in the order: theta -> r1 -> r2")')
    else
        if (.not.ztwod) then
            if (zr2r1) write(6,'(5x,"Problem solved in the order: r2 -> r1 -> theta")')
            if (.not. zr2r1) write(6,'(5x,"Problem solved in the order: r1 -> r2 -> theta")')
        else
            if (zr2r1) write(6,'(5x,"Problem solved in the order: r2 -> r1")')
            if (.not. zr2r1) write(6,'(5x,"Problem solved in the order: r1 -> r2")')
        endif
    endif
    if (zcut) then
        write(6,'(5x,"Final basis selected using energy cut-off")')
    else
        if (zrot .and. zbisc .and. (jrot+kmin)>1) then
            write(6,'(/5x,"Final basis comprises",i5," lowest functions",&
                   " for even parity hamiltonian"/&
                5x,"Final basis comprises",i5," lowest functions",&
                   " for odd  parity hamiltonian")') max3d,max3d2

        else
           if (.not. zcut) write(6,'(5x,"final basis comprises",i5," lowest functions")') max3d
        endif
    endif
    if (.not.ztwod) then
        if (zlmat) write(6,'(/5x,"Printing of L-matrix requeste")')
        if (.not.zlmat) write(6,'(/5x,"Printing of L-matrix not requested")')
    endif

    if (zp1d) write(6,'(5x,"Printing of 1d eigenvalues requested")')
    if (zp2d) write(6,'(5x,"Printing of 2d eigenvalues requested")')
    if (zpham) write(6,'(5x,"Printing of hamiltonian matrix requested")')
    if (.not.zpham) write(6,'(5x,"Printing of hamiltonian matrix not requested")')
    if (zprad) write(6,'(5x,"Printing of radial matrix elements requested")')
    if (.not.zprad) write(6,'(5x,"Printing of radial matrix elements not requested")')
    if (zpvec) write(6,'(5x,"Printing of eigenvectors requested"/)')
    if (.not.zpvec) write(6,'(5x,"Printing of eigenvectors not requested"/)')

    if (zvec) then
        write(6,'(5x,"Eigenvalues & vectors   written to stream IOUT2 =",i4)') iout2
        write(6,'(5x,"Restart information     written to stream IOUT1 =",i4)') iout1
        if (zpseg==.true.) then 
            open(unit=iout1, form='unformatted',recordtype='segmented')
            open(unit=iout2, form='unformatted',recordtype='segmented')
        else 
            open(unit=iout1, form='unformatted')
            open(unit=iout2, form='unformatted')
        end if 
    endif

    if (ztran) then
        write(6,'(5x,"Wavefunction amplitudes written to stream IWAVE =",i4)') iwave
        if (zpseg==.true.) then
            open(unit=iwave, form='unformatted',recordtype='segmented')
        else 
           open(unit=iwave, form='unformatted')
        end if 
      endif

      if (abs(jrot) > 1) zpfun=.false.
      if (zpfun) then
        open(unit=ilev,form='formatted')
        if (jrot == 0 .and. mod(ipar,2) == 0) then
!           header on file ilev
            write(ilev,'(9a8)') title
            write(6,'(/5x,"Eigenvalues      written to start of stream ilev =",i4)') ilev
        else
!           position file ilev
  200       read(ilev,*,end=210,err=210)
            goto 200
  210       continue
! ******  inclusion of the following card is machine dependent ******
!           backspace ilev
            write(6,'(/5x,"Eigenvalues      written at end   of stream ilev =",i4)') ilev
        endif
    endif

    if (idia > 0) write(6,'(/5x,"Calculation performed in scattering coordinates")')
        if (idia <= 0 .and.  .not.zperp) write(6,'(/5x,"Calculation performed in Radau coordinates")')
        if (idia <= 0 .and. zperp)&
            write(6,'(5x,"Calculation performed in Radau coordinates with Z axis",/5x,"perpendicular to the plane")')
        if (zperp) write(6,"(/5x,'Q symmetry =', i3)") iq
        if (ztwod) goto 886
        if (abs(idia) == 2) then
            write(6,'(/5x,"Diatomic assumed homonuclear")')
            if (ipar == 1) then
                write(6,'(5x,"Odd parity functions in basis set")')
            else if (ipar == 0) then
                write(6,'(5x,"Even parity functions in basis set")')
        else
            write(6,'(5x,"Illegal value of ipar for idia = +/-2: STOP")')
            stop
        endif
        if (idia == 2) then
            idvr=nalf2
            if (2*idvr /= nalf) goto 960
        endif
    else
        write(6,'(/5x,"Diatomic assumed hetronuclear")')
        idvr=nalf
        ipar=0
    endif

    if (jrot /= 0) then
        jrot=abs(jrot)
        if (zrot) then
            if (kmin /= 0 .and. .not. zbisc) kmin=1
            write(6,'(/5x,"***  vibrational part of rot-vib calculation  ***")')
            write(6,'(/5x,"J =",i3," rotational state")') jrot
            if (kmin == 1) then
               write(6,'(12x,"with symmetric |Jk> + |J-k> functions in basis")')
            else if (kmin == 0) then
               write(6,'(12x,"with anti-symmetric |Jk> - |J-k> functions in basis")')
            else
               kmin=2
               write(6,'(12x,"loop over symmetric & anti-symmetric |jk> functions")')
            endif
            if (zladd) write(6,'(5x,"Number of angular grid points to be kept constant with k")')
            if (.not. zladd) write(6,'(/5x,"Nunber of angular grid points to decrease with k")')
        else
            write(6,&
      & '(5x,"J =",i3,"  k =",i3,/5x,"***  option to neglect coriolis interactions  ***")') jrot,kmin
            if (abs(kmin) > abs(jrot)) then
               write(6,'(5x,"Error: k greater than J. STOP")')
               stop
            endif 
        endif
        if (zbisc) then
            zembed=.false.
            write(6,'(/5x,"z axis embedded along the biscetor of r1 and r2")')
            if (zlin) write(6,"(/5x,'Removal of functions with theta = 0 enforced')")
        else
            if (zembed) write(6,"(/5x,'z axis embeded along the r',i1,' coordinate')") 2
            if (.not. zembed) write(6,"(/5x,'z axis embeded along the r',i1,' coordinate')") 1
        endif
    else
!        case j = 0
        write(6,"(/5x,'J =',i3,' rotational state')") jrot
    endif

    if (zx) then
        write(6,"(5x,'the rotational x-term is active')")
        xp0=1.d0
    else
        write(6,"(5x,'the rotational x-term is ignored')")
        xp0=0.d0
    end if

    if (zs0) then
    write(6,"(5x,'diagonal NBO vibrational corrections are implemented:')")
    write(6,"(5x,'Those corrections are active only for the Radau hamiltonian')")
    write(6,"(5x,'(BODC Corrections follow procedure by D. Schwenke...)')")
    xp1=1.d0
    else
    xp1=0.d0
    write(6,"(5x,'diagonal NBO vibrational corrections are ignored')")
    end if

    if (zs1) then
    write(6,"(5x,'NBO vib-rotational corrections are implemented:')")
    write(6,"(5x,'Those corrections are active only for the Radau hamiltonian')")
    write(6,"(5x,'(see MNRAS 480, 2597–2608 (2018) doi:10.1093/mnras/sty1877)')")
    xp2=1.d0
    else
    xp2=0.d0
    write(6,"(5x,'NBO vib-rotational corrections are ignored')")
    end if

    write (6,"(/5x,'Routine DSYEV to do in core diagonalisation')")
  886 continue

!     check input parameters are consistent
    npnt = max(npnt1,npnt2)
    nmax = max(nmax1,nmax2)
!     dimension of square matrices stored in triangular form :
    nlim1 = npnt1 * (npnt1+1) / 2
    nlim2 = npnt2 * (npnt2+1) / 2

!     declare dvr sizes for dimensioning the arrays
    ndima=npnta
    if (idia <= -2 .and. .not. zrot) ndima=ndima-ipar
    ndimb=npntb
    ndimc=npntc
!     store parameters on disk files requested

    if (zvec) write (iout2) idia,ipar,npnta,npntb,npntc,max2d,max3d,neval

    if (zmors2 .and. .not. zquad2) goto 961
    if (.not. ztheta .and. .not. zquad2) goto 962
    if (idia <= -2 .and. .not. zquad2) goto 963

    return
960 write(6,"(//6x,'** nalf must be even when idia=2: stop **')")
    stop
961 write(6,972)
972 format(//6x,'** can''t have zquad2 = f with zmors2 = t : stop **',&
              /6x,'               (not yet implemented)               ')
    stop
962 write(6,973)
973 format(//6x,'** can''t have zquad2 = f with ztheta = f : stop **',&
              /6x,'               (not yet implemented)               ')
    stop
963 write(6,974)
974 format(//6x,'** can''t have zquad2 = f with idia = -2: stop **',&
              /6x,'               (not yet implemented)               ')
    stop
end

!##############################################################################
subroutine ccmain

!     subroutine ccmain is the 'real' main programme & contains
!     the calls to the various subroutines which set & solve the
!     intermediate and the final hamiltonians.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb
    use dvr3drjz_mass
    use dvr3drjz_split1
    use dvr3drjz_split2
    implicit double precision (a-h,o-y), logical (z)

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dnorm1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r1m2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dnorm2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r2m2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bass1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bass2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: y1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wt1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: y2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wt2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: b
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: c
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: hbl1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: hbl2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xalf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: walf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pleg
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: plegw
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xlmatr
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xk1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xk2 
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r2m2t
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: jxcos
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: jwalf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sjwalf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pjac
 
    x0 = 0.0d0
    x1 = 1.0D0
    x2 = 2.0D0
    x8 = 8.0d0
    xp = 0.5d0
    toler = 1.0d-8

    ALLOCATE(dnorm1(0:nmax1),r1m2(nlim1),dnorm2(0:nmax2),&
            r2m2(nlim2),bass1(0:nmax1,npnt1),bass2(0:nmax2,npnt2),&
            y1(npnt1),r1(npnt1),wt1(npnt1),y2(npnt2),r2(npnt2),&
            wt2(npnt2),b(npnt+1),c(npnt+1),&
            hbl1(nlim1),hbl2(nlim2),xalf(idvr),walf(idvr),&
            xlmatr(idvr,idvr),pleg (0:maxleg,idvr),&
            xk1(npnt1,npnt1),xk2(npnt2,npnt2),r2m2t(npnt2,npnt2))

    if (zperp) ALLOCATE(sjwalf(idvr),&
                jxcos(idvr),jwalf(idvr),pjac(0:maxleg,idvr),&
                plegw(0:maxleg,idvr))

!     open the streams......(most not needed for symmetrised radau case)
    if (idia > -2) then
        if (zpseg==.true.) then 
        open(unit=ieigs1,form='unformatted',recordtype='segmented')
        open(unit=ieigs2,form='unformatted',recordtype='segmented')
        open(unit=ivecs1,form='unformatted',recordtype='segmented')
        open(unit=ivecs2,form='unformatted',recordtype='segmented')
        open(unit=ivint, form='unformatted',recordtype='segmented')
        open(unit=iband, form='unformatted',recordtype='segmented')
    else 
        open(unit=ieigs1,form='unformatted')
        open(unit=ieigs2,form='unformatted')
        open(unit=ivecs1,form='unformatted')
        open(unit=ivecs2,form='unformatted')
        open(unit=ivint, form='unformatted')
        open(unit=iband, form='unformatted')
       end if 
    endif

    if (zpseg==.true.) then
        open(unit=intvec,form='unformatted',recordtype='segmented')
    else 
        open(unit=intvec,form='unformatted')
    end if 
!     read in masses & basis set parameters
    call setcon(fixcos)
!.....and save them if necessary
    if (zvec) write(iout2) zembed,zmors1,zmors2,ztheta,zr2r1,xmass,&
                            g1,g2
!     set up binomial and normalisation arrays
    call setfac(dnorm1,dnorm2,cc1,cc2)

!     set up points, weights & basis for numerical integration

    if (ncoord == 2) then
!     in less than 3-d cases fix the 3-d paramters
        if (zr2r1) then
            bass1(0,1) = x1
            wt1(1) = x1
            r1(1) = re1
            call lagpt(2,y2,r2,wt2,b,c,cc2,bass2,dnorm2,npnt2,nmax2,&
                       zmors2,re2,beta2,a2,iu2)
            hbl2(1) = x0
!           setup kinetic energy & inertia integrals over r1
            if (zmors1) then
               fke = beta1 * beta1 / (x8 * ur1)
               call keints(hbl1,fke,nlim1,nmax1,iu1)
            else
               fke = beta1 / (x2 * ur1)
               call keint2(hbl1,fke,r1m2,dnorm1,nlim1,nmax1,npnt1,a1)
            endif
        else
            bass2(0,1) = x1
            wt2(1) = x1
            r2(1) = re2
            call lagpt(1,y1,r1,wt1,b,c,cc1,bass1,dnorm1,npnt1,nmax1,&
                       zmors1,re1,we1,a1,iu1)
            hbl1(1) = x0
!           setup kinetic energy & inertia integrals over r2
            if (zmors2) then
               fke = beta2 * beta2 / (x8 * ur2)
               call keints(hbl2,fke,nlim2,nmax2,iu2)
            else
               fke = beta2 / (x2 * ur2)
               call keint2(hbl2,fke,r2m2,dnorm2,nlim2,nmax2,npnt2,a2)
            endif
        endif
    else
        call lagpt(1,y1,r1,wt1,b,c,cc1,bass1,dnorm1,npnt1,nmax1,zmors1,&
                    re1,beta1,a1,iu1)
        if (idia > -2) then
            call lagpt(2,y2,r2,wt2,b,c,cc2,bass2,dnorm2,npnt2,nmax2,&
                       zmors2,re2,beta2,a2,iu2)
!           setup kinetic energy & inertia integrals over r2
            if (zmors2) then
              fke = beta2 * beta2 / (x8 * ur2)
              call keints(hbl2,fke,nlim2,nmax2,iu2)
            else
              fke = beta2 / (x2 * ur2)
              call keint2(hbl2,fke,r2m2,dnorm2,nlim2,nmax2,npnt2,a2)
            endif
        endif
!        setup kinetic energy & inertia integrals over r1
        if (zmors1) then
            fke = beta1 * beta1 / (x8 * ur1)
            call keints(hbl1,fke,nlim1,nmax1,iu1)
        else
            fke = beta1 / (x2 * ur1)
            call keint2(hbl1,fke,r1m2,dnorm1,nlim1,nmax1,npnt1,a1)
        endif
    endif

!     write the quadrature points to disk for zvec = .true.
    if (zvec) then
        write (iout2) r1
        if (idia > -2) write (iout2) r2
    endif
    if (ncoord == 3) then
    endif

!     take square roots of the weights
    wt1=sqrt(wt1)
    if (idia > -2) wt2=sqrt(wt2)

!     set up the transformed kinetic energy integrals,  t'(hbl) t
!                                                       ~  ~~~  ~
    call k1k2(xk1,hbl1,bass1,wt1,npnt1,nmax1,nlim1)
    if (idia > -2)&
        call k1k2(xk2,hbl2,bass2,wt2,npnt2,nmax2,nlim2)

!     ...... and the inertia integrals for spherical oscillators
    if (.not.zmors2 .and. .not.ztwod .and. idia > -2)&
        call k1k2(r2m2t,r2m2,bass2,wt2,npnt2,nmax2,nlim2)

!     some of the j>0 stuff to get the loop over k correct
    if (zrot) then
        kd = 1-min(kmin,1)
        ku = jrot
!       if looping over sym & anti-sym, do one extra k=1 block
        if (kmin == 2) ku=ku+1
    else
        kd = kmin
        ku = kmin
    endif
    kkz12 = 0
    kkz0  = 0

!     for j > 0, store r**(-2) term for rotlev3
    if (ztran) then
        if (ku > kd) then
            if (zembed) then
                if (zquad2) then
                    write(iwave) (xp5/(r2(i)*r2(i)*urr2),i=1,npnt2)
                else
                    write(iwave) ((ur2*r2m2t(i,j)/urr2,i=1,npnt2),j=1,npnt2)
                endif
            else
                write(iwave) (xp5/(r1(i)*r1(i)*urr1),i=1,npnt1)
            endif
        else
            jdia=max(1,idia)
            jstart=kmin
            if (mod(jstart,jdia) /= ipar) jstart=jstart+1
            nang=(maxleg-jstart)/jdia+1
            mbass=idvr*npnt1*npnt2
            if (idia == -2) mbass=idvr*max2d
            write(iwave) mbass,jstart,nang,mbass
        endif
        write(iwave) r1
        if (idia > -2) write(iwave) r2
    endif

    if (idia == -2) then
        max2d1=max2d
        max3d1=max3d
        DEALLOCATE(xk2,r2,r2m2t)
    endif

    DEALLOCATE(dnorm1,r1m2,dnorm2,r2m2,bass1,bass2,&
                y1,wt1,y2,wt2,b,c,hbl1,hbl2)

!     -------------  start rotational loop here  -------------
    do 40 kk=kd,ku
    if (kk <= jrot) then
        kz=kk
    else
        kz=1
        kmin=0
        if(.not.zperp) ipar=mod(ipar+jrot,2)
    endif

!     first rewind the scratch files for a calculation with j>0
!     and, if needed, reposition iout2 after set up data.
    if (kk > kd) then
        if (idia > -2) then
           rewind ieigs1
           rewind ieigs2
           rewind ivecs1
           rewind ivecs2
           rewind ivint
           rewind iband
        endif
        rewind intvec
        if (zvec) then
          rewind iout1
          rewind iout2
          do 45 ii=1,4
            read(iout2)
   45     continue
        endif
      endif

    realkz = dble(kz)
!     tswalf is the exact sum of weights for gauss-jacobi integration
    tswalf= x2**(kz+kz+1)/dble(kz+1)
    do 30 ia=1,kz
        tswalf=tswalf*dble(ia)/dble(kz+ia+1)
   30 continue

    if (zladd .or. kz == kd) then
        nidvr = idvr
        nang  = nalf
        nang2 = nalf2
        lincr = kz
    else
        lincr = 0
        if (idia /= 2) then
            nidvr = idvr - kz
            nang  = nalf - kz
            nang2 = (nang+1)/2
        else if (ipar == 0) then
            if(mod(kz,2)==1) kkz0 = kkz0 + 2
            nidvr = idvr - kkz0/2
            nang  = nalf - kkz0
            nang2 = (nang+1)/2
        else
            if(mod(kz,2)==0 .and. kz>0) kkz12 = kkz12 + 2
            nidvr = idvr - kkz12/2
            nang  = nalf - kkz12
            nang2 = (nang+1)/2
        endif
    endif
    if (.not. zladd) then
        if (ztheta) then
            npnta=nidvr
        else
            npntc=nidvr
        endif
    endif

    if (ztwod) then
        xalf(1) = fixcos
        goto 333
    endif

    if(.not. zperp) then
        call jacobi(nang,nang2,xalf,walf,realkz,realkz,cswalf,tswalf)
        write(6,1000) nang,kz,(xalf(ii),walf(ii),ii=1,nang2)
  1000   format(//i8,' point Gauss-associated Legendre integration with',&
              ' k =',i3//5x,'integration points',11x,'weights',&
               //(f23.15,d25.12))
        write(6,1010) cswalf,tswalf
 1010   format(/4x,'Computed sum of weights',d22.15,&
              /4x,'Exact    sum of weights',d22.15//)
          if (abs((cswalf-tswalf)/tswalf) > toler) then
             write(6,"(//5x,'Points & weights in error, adjust algorithm'//)")
             stop
          endif
        call allpts(xalf,walf,nang,nang2)
        alf = x0
    else
        realj = DBLE(jrot)
        alf = dsqrt(0.5d0*(((realj)**2+realj)-realkz**2))
        bet = alf
        call gaujac(xalf,walf,2*nang2,alf,bet)
        write(6,1020) nang,alf,(xalf(ii),walf(ii),ii=1,nang2)
 1020   format(//i8,' point Gauss-Jacobi integration with',&
             ' alpha = beta =',f7.3//5x,'integration points',11x,'weights',&
              //(f23.15,d25.12))
    endif

    if (zvec) write(iout2) xalf

    if (.not.zperp) then
!     set up Legendre polynomials for the transformation matrix
!     for z in plane case
    call asleg(pleg,maxleg,xalf,nidvr,kz,lincr)
    else
!     set up Jacobi polynomials for the transformation matrix
!     for z perpendicular case
    call jac_basis(nidvr,maxleg,alf,bet,xalf,pleg)
        do  i=1,nidvr
            walf(i) = sqrt(walf(i))
        enddo
    endif

!     save polinomials for rotlev3 or rotlev3b, or rotlev3z
    if (ztran) then
        write(iwave) xalf
        write(iwave) kz,maxleg,nidvr,lincr
        write(iwave) ((pleg(i,j)*walf(j),i=0,maxleg),j=1,nidvr)
    endif
      
!     build the transformed angular momentum matrix xlmatr;
    ipar0=0
    if (idia == 2 .and. ipar == 1) ipar0=1
    call lmatrx(xlmatr,pleg,walf,kz,ipar0,nidvr,lincr,alf)

  333 continue

!     for ab2 molecules in radau coordinates, use separate main
!     driving routine
    if (idia <= -2) then
        if (zrot) then
            if (kz > kd .and. .not.zperp) ipar=mod(ipar+1,2)
            if (ipar > 0) then
               max2d=max2d2
               max3d=max3d2
            else
               max2d=max2d1
               max3d=max3d1
            endif
        endif
        if (zperp) then
           call mkmain(xk1,xlmatr,r1,xalf,kz)
        else
           call nfmain(xk1,xlmatr,r1,xalf,kz)
        endif    
    else
        call jhmain(xk1,xk2,xlmatr,r1,r2,xalf,r2m2t,kz)
    endif
   40 continue

    DEALLOCATE(xalf,walf,xlmatr,pleg,xk1,r1)
    if (zperp) DEALLOCATE(jxcos,jwalf,sjwalf,pjac,plegw)
    if (idia > -2) DEALLOCATE(r2,xk2,r2m2t)
    return
end

!#############################################################################
subroutine setcon(fixcos)

!     read in masses & set constants for radial basis sets          #007
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_split1
    use dvr3drjz_split2
    use dvr3drjz_oupb
    use dvr3drjz_mass
    implicit double precision (a-h,o-y), logical (z)

!     amtoau converts amu (proton masses) to au (electron masses).
    amtoau = 1.8228883d03
    x0 = 0.0d0
    xp5 = 0.5d0
    x1 = 1.0d0
    x4 = 4.0d0

!     read cos(theta) for fixed angle 2-d calculation

    read(5,'(3f20.0)') fixcos
    if (ztwod) write(6,1088) fixcos
1088 format(//5x,'two-d fixed angle vibrational problem with'&
             //5x,'***** fixed value of cos(theta) =',f6.2,' *****'/)

!     read masses of the atoms in atomic mass units
!     first vibrational mass...
    read(5,'(3f20.0)')     xmass
!     .... then rotational mass
    read(5,'(3f20.0)')     xmassr
!     Default rotational mass to vibration mass if it is not set
    if (xmassr(1)<=x0) xmassr = xmass

!     read cut off energies
!     read parameters defining energy cut offs for each block
    read(5,'(3f20.0)') emax1,emax2
    if (.not. zall) then
        if (zcut) then
            if (idia > -2) write(6,"(//5x,'Cut-off energies in wavenumbers:',2d16.8/)") emax1,emax2
            if (idia == -2) write(6,"(//5x,'Final cut-off energy in wavenumbers:',2d16.8/)") emax2
        else
            if (idia > -2) write(6,"(//5x,'First cut-off energy in wavenumbers:',1d16.8/)") emax1
        endif
    endif
!     set default value of g1 and g2
    if (idia >= 1) then
!        scattering coordinates
        g1 = xmass(2) / (xmass(2) + xmass(3))
        g2 = x0
    else
!        radau coordinates
        a = sqrt(xmass(3) / (xmass(1)+xmass(2)+xmass(3)))
        b = xmass(2) / (xmass(1)+xmass(2))
        g1 = x1 - a / (a+b-a*b)
        g2 = x1 - a / (x1-b+a*b)
    endif
!     ncoord = 3: read parameters for r1 radial basis (see below)
!     ncoord = 2: read fixed r1 bondlength re1. diss1 & we1 dummy
    read(5,'(3f20.0)')     re1,diss1,we1
!     read parameters for r2 radial basis function,
!     for morse oscillator functions use the following:
!     re2: equilibrium bondlength of r2 coordinate (in bohr)
!     diss2: dissociation energy of the r2 coordinate (in hartree)
!     we2: fundamental stretching vibration of r2 (in hartree)
!     for spherical oscillator functions use the following:
!     re2 : dummy
!     diss2: order of laguerre polynomials used (dimensionless)
!     we2: fundamental stretching vibration of r2 (in hartree)
!     all are treated as variationally optimisable parameters.
    read(5,'(3f20.0)')     re2,diss2,we2
    write(6,'(/5X,"Vibrational nuclear mass in AMU:",3F12.6)') xmass
    if (jrot/=0) write(6,'(5X,"Rotational  nuclear mass in AMU:",3F12.6/)') xmassr
!     compute the effective moments of inertia
    ur1 = amtoau/(g2*g2/xmass(1)+x1/xmass(2)+(x1-g2)**2/xmass(3))
    ur2 = amtoau/(x1/xmass(1)+g1*g1/xmass(2)+(x1-g1)**2/xmass(3))
    urr1 = amtoau/(g2*g2/xmassr(1)+x1/xmassr(2)+(x1-g2)**2/xmassr(3))
    urr2 = amtoau/(x1/xmassr(1)+g1*g1/xmassr(2)+(x1-g1)**2/xmassr(3))
 
    if (ncoord == 3) goto 20
    if (zr2r1) then
        write(6,"(/5x,'r1 fixed bondlength =',f8.4,' bohr',' & reduced mass =',d16.7,' a.u.'/)") re1,ur1

    else
        write(6,"(/5x,'r2 fixed bondlength =',f8.4,' bohr',' & reduced mass =',d16.7,' a.u.'/)") re2,ur2
    endif
    if (zr2r1) goto 30
   20 continue
    if (zmors1) then
        write(6,1020) 1,re1,diss1,we1
1020 format(/5x,'Morse function parameters for r',i1,' basis',&
             /5x,'r equilibrium =',f8.4,' bohr, dissociation energy',&
        d15.7,' hartree &  vibrational frequency =',d15.7,' hartree')
        beta1 = we1 * sqrt(xp5*ur1/diss1)
        a1 = x4 * diss1 / we1
        iu1 = int(a1+xp5)
        write(6,1030) ur1,beta1,a1,iu1
1030 format(/5x,'Constants used to construct morse oscillators:',&
            /5x,'reduced mass =',d16.7,' a.u., beta =',f8.4,&
                ' (1/bohr), a =',d16.7,' and u =',i5)
    else
        a1=diss1
        beta1 = sqrt(we1 * ur1)
        write(6,1039) a1,we1,ur1,beta1
1039 format(/5x,'Spherical oscillator parameters for r1 basis:',&
            /5x,'alpha =',f10.5,&
                ' &  vibrational frequency =',d15.7,' hartree',&
            //5x,'Constants used to construct spherical oscillators:',&
            /5x,'reduced mass =',d16.7,' a.u., beta =',f12.6,&
                ' bohr**-2')
    endif
    if (ncoord == 2 .and. .not. zr2r1) goto 40
    if (idia == -2) goto 40
   30 continue
    if (zmors2) then
        write(6,1020) 2,re2,diss2,we2
        beta2 = we2 * sqrt(xp5*ur2/diss2)
        a2 = x4 * diss2 / we2
        iu2 = int(a2+xp5)
        write(6,1030) ur2,beta2,a2,iu2
    else
        a2=diss2
        beta2 = sqrt(we2 * ur2)
        write(6,1040) a2,we2,ur2,beta2
1040 format(/5x,'Spherical oscillator parameters for r2 basis:',&
            /5x,'alpha =',f10.5,&
                ' &  vibrational frequency =',d15.7,' hartree',&
            //5x,'Constants used to construct spherical oscillators:',&
            /5x,'reduced mass =',d16.7,' a.u., beta =',f12.6,&
                ' bohr**-2')
    endif
   40 continue
    if (ztran) then
         if (zrot) then
            if (zembed) then
               if (zquad2) then
                  nlim = npnt2
               else
                  nlim = npnt2*npnt2
               endif
            else
               nlim = npnt1
            endif
        else
            nlim = 0
        endif
        zncor=.not.zrot
        if (idia > -2) then
            write(iwave) idia,ipar,idvr,npnt1,npnt2,jrot,kmin,neval,nlim
            write(iwave) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
            write(iwave) re1,diss1,we1,re2,diss2,we2
        else
            write(iwave) idia,ipar,idvr,npnt1,npnt1,jrot,kmin,neval,nlim
            write(iwave) zembed,zmors1,zmors1,xmass,g1,g2,zncor,zquad2
            write(iwave) re1,diss1,we1,re1,diss1,we1
        endif
    endif
    return
end

!#########################################################################
subroutine setfac(dnorm1,dnorm2,cc1,cc2)

!     setfac initialises binomial array:                            #021
!       binom(i+1,j+1) = i! / (j! * (i-j)!)
!     and pseudo-normalisation array:
!       dnorm(m) = sqrt((m-1)! * binom(npnt+iu,npnt-m))
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_split1
    use dvr3drjz_split2
    use dvr3drjz_oupb
    implicit double precision (a-h,o-y), logical (z)


    double precision, dimension(0:nmax1) :: dnorm1
    double precision, dimension(0:nmax2) :: dnorm2
    double precision, dimension(npnt+1) :: fact

    data xp5/0.5d0/,x1/1.0d0/
    fact(1) = x1
    count = x1
    do 10 i=1,npnt
    fact(i+1) = count * fact(i)
    count = count + x1
   10 continue
    if (ncoord == 2) then
        if (.not. zr2r1) then
          if (zmors1) then
            alf=dble(iu1)
          else
            alf=a1+xp5
          endif
          call norms(dnorm1,fact,cc1,alf,npnt1,nmax1)
        else
          if (zmors2) then
            alf=dble(iu2)
          else
            alf=a2+xp5
          endif
          call norms(dnorm2,fact,cc2,alf,npnt2,nmax2)
        endif
    else
        if (zmors1) then
          alf=dble(iu1)
        else
          alf=a1+xp5
        endif
        call norms(dnorm1,fact,cc1,alf,npnt1,nmax1)
        if (idia > -2) then
           if (zmors2) then
             alf=dble(iu2)
           else
             alf=a2+xp5
           endif
           call norms(dnorm2,fact,cc2,alf,npnt2,nmax2)
        endif
    endif

    return
end

!##############################################################################
subroutine norms(dnorm,fact,cc,alf,npnt,nmax)
!     set up factors for normalising the radial basis functions     #022

    implicit double precision (a-h,o-y)

    double precision, dimension(0:nmax) :: dnorm
    double precision, dimension(npnt+1) :: fact
    double precision, dimension(npnt+1) :: bin

    data x1/1.0d0/
    count = dble(npnt) + alf
!     cc is exact sum of weights for npnt Gauss-Laguerre integration
    cc = fact(npnt) / count
!     normalisation array for l(i,alf): first set up binomials
    npt1 = npnt + 1
    bin(npt1) = x1
    do 20 i=1,npnt
    n = npt1 - i
    bin(n) = bin(n+1) * count / dble(i)
    count = count - x1
   20 continue
    npt1 = npt1 + 1
    do 30 i=0,nmax
    dnorm(i) = sqrt(bin(i+1)*fact(i+1)*fact(npt1-i-1))
   30 continue
    return
end

!############################################################################
subroutine lagpt(ir,y,r,wt,b,c,cc,bass,dnorm,npnt,nmax,zmorse,&
                       re,beta,a,iu)

!     subroutine lagpt gets integration points and weights for      #015
!     npnt gauss laguerre integration and sets up basis
!     functions at the integration points.

    implicit double precision(a-h,o-y), logical (z)

    double precision, dimension(npnt+1) :: b
    double precision, dimension(npnt+1) :: c
    double precision, dimension(npnt) :: y
    double precision, dimension(npnt) :: r
    double precision, dimension(npnt) :: wt
    double precision, dimension(0:nmax) :: dnorm
    double precision, dimension(0:nmax,npnt) :: bass
         
    data x0,xp5,x1,x2/0.0d0,0.5d0,1.0d0,2.0d0/,toler/1.0d-8/

    if (zmorse) alf=dble(iu)
    if (.not. zmorse) alf = a + xp5
    alfm1=alf-x1
!     set up integration points and weights.
    call laguer(npnt,y,wt,alf,b,c,csx,csa,tsx,cc)
    tsa = x1 / (dnorm(0) * dnorm(0))
    write(6,1000) npnt,ir
 1000 format(/,i8,' point Gauss-Laguerre integration',&
             /,5x,'integration points',11x,'weights',9x,&
                  'corresponding r',i1,/)
    do 60 i=1,npnt
    if (zmorse) then
!         calculate potential at r = re+beta(**-1)*ln(a/y)
        r(i) = re + dlog(a/y(i)) / beta
    else
!         calculate potential at r = sqrt(y/beta)
        r(i) = sqrt(y(i)/beta)
    endif
    write(6," (f23.15,d25.12,f13.5)") y(i),wt(i),r(i)
    if (r(i) < x0) &
    write(6,"(5x,'***** warning: for integration point',i3,', r less than zero *****')") i


!     calculate unnormalised laguerre polynomials at y

!     polynomial of order 0
    bass(0,i) = x1
    if (nmax < 1) goto 70
!     polynomial of order 1
    amx = alf + x1 - y(i)
    bass(1,i) = amx
!     use recurrence relationships for polynomials of order > 2
!     n * l(n,alf) = (2*n+alf-1-x)*l(n-1,alf) - (n+alf-1)*l(n-2,alf)
    en = x1
    do 80 n=2,nmax
    en = en + x1
    amx = amx + x2
    bass(n,i) = (amx * bass(n-1,i) - (alfm1+en) * bass(n-2,i)) / en
   80 continue
   70 continue

    do 90 n2=0,nmax
!     normalise polynomials
    bass(n2,i) = bass(n2,i) * dnorm(n2)
   90 continue
   60 continue
!     check that the correct points & weights have been generated
    write(6,1020) csx,csa,tsx,tsa
1020 format(/4x,'Computed sum of points',d22.15,' & weights',d22.15,&
             /4x,'Exact    sum of points',d22.15,' & weights',d22.15)
    if (abs((csx-tsx)/tsx) > toler) goto 900
    if (abs((csa-tsa)/tsa) > toler) goto 900
    return
900 write(6,"(//5x,'points & weights in error, adjust algorithm',//)")
    stop
end

!##########################################################################
subroutine laguer(nn,x,a,alf,b,c,csx,csa,tsx,cc)

!     calculates points & weights for gauss-laguerre integration    #016
!     see:
!     "gaussian quadrature formulas" by a.h.stroud & d.secrest
!      1966, prentice-hall, p.32.
!     **** version to avoid overflows (j.t. 25/11/81) ****
!     calculates weights divided by gamma(nn+alf+1)
!     this is an initialsation entry

    implicit double precision (a-h,o-y)
    double precision, dimension(nn) :: x
    double precision, dimension(nn) :: a
    double precision, dimension(nn+1) :: b
    double precision, dimension(nn+1) :: c      

    eps = 1.0d-12
    x1 = 1.0d0
    csx=0.0d0
    csa=0.0d0
    fa=alf+1.0d0
!     cc = n!                      denominator for pseudo-weights: a
!     b(n) = (alf + 2n -1)             b & c for recurrence relation
!     c(n) = (n - 1) * ( alf + n - 1)
    b(1)=fa
    c(1)=0.0d0
    fn=1.0d0
    do 1 j=2,nn
        fa=fa+2.0d0
        b(j)=fa
        c(j)=fn*(alf+fn)
        fn=fn+1.0d0
    1 continue
    tsx=fn*(alf+fn)
    xt1=0.0d0
!     formulas for initial point & step chosen because they work!
    xt=(1.0d0+alf)*(2.0d0+alf)/(1.0d0+3.0d0*fn+2.0d0*alf)
    step=3.0d0*(1.0d0+alf)/(1.0d0+3.0d0*fn+alf)
    call lgrecr(pt,dpn,pn1,xt,nn,alf,b,c)

    do 7 i=1,nn
    if (i > 2) goto 4
!     smallest two zeros: found by "brute force" search
    2 xt2 = xt + step
    call lgrecr(pt2,dpn,pn1,xt2,nn,alf,b,c)
    if (dsign(x1,pt)*dsign(x1,pt2) > 0.0d0) goto 5
    pt = pt2
    xt = 0.5d0 * (xt + xt2)
    go to 6
    5 pt = pt2
    xt = xt2
    go to 2
!     all other zeros: found using formula of stroud & secrest
    4 fi = dble(i-2)
    r1 = (1.0d0+2.55d0*fi)/(1.9d0*fi)
    r2 = 1.26d0*fi*alf/(1.0d0+3.5d0*fi)
    ratio = (r1+r2)/(1.0d0+0.3d0*alf)
    xt = xt + ratio*(xt-xt2)

    6 call lgroot(xt,nn,alf,dpn,pn1,b,c,eps)
    xt2=xt1
    xt1=xt
    x(i) = xt
    a(i) = cc/dpn/pn1
    csx = csx + xt
    csa = csa + a(i)
    7 continue
    return
    end

!##############################################################################
    subroutine lgroot(x,nn,alf,dpn,pn1,b,c,eps)

!     improves the approximate root x; in addition obtains          #017
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x
!     this routine is due to stroud & secrest (see subroutine laguer)

    implicit double precision (a-h,o-y)

    double precision, dimension(nn+1) :: b
    double precision, dimension(nn+1) :: c     

    itmax = 10
    iter=0
    1 iter=iter+1
    call lgrecr(p,dpn,pn1,x,nn,alf,b,c)
    d = p/dpn
    x = x-d
    if (abs(d/x) <= eps) return
    if (iter < itmax) goto 1 
    write(6,100) iter,d,x
100 format(5x,'warning: noconvergence after',i4,' iterations',&
            /,5x,'current difference',d26.15,' & root',d26.15)
    return
end

!###########################################################################
subroutine lgrecr(pn,dpn,pn1,x,nn,alf,b,c)

!     uses recurrence relations to set up polynomials               #018
!     this routine is due to stroud & secrest (see subroutine laguer)

    implicit double precision (a-h,o-y)

    double precision, dimension(nn+1) :: b
    double precision, dimension(nn+1) :: c     

    p1 = 1.0d0
    p = x - alf - 1.0d0
    dp1 = 0.0d0
    dp = 1.0d0
    do 1 j=2,nn
        q  = (x-b(j))* p-c(j)* p1
        dq = (x-b(j))*dp-c(j)*dp1 + p
        p1 = p
        p  = q
        dp1= dp
        dp = dq
    1 continue
    pn = p
    dpn= dp
    pn1= p1
    return
end

!###########################################################################
subroutine keints(hbl,fke,nlim,nmax,iu)

!     keints calculates analytic kinetic energy integrals over r    #012
!     for morse oscillator-like functions
    use dvr3drjz_outp
    use dvr3drjz_oupb
    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(nlim) ::  hbl

    x0 = 0.0d0
    index = 0
    do 10 n2=1,nmax+1
        do 20 n1=1,n2
          index=index+1
          if (n1 == n2) then
!             special case:  n1 = n2
              hbl(index) = fke * dble(2*(n2-1)*(n2+iu)+iu+1)
          else if (n1+2 == n2) then
!                  special case:  n2 = n1 + 2
                   hbl(index) = - fke * sqrt(dble((iu+n1+1)*(iu+n1))*dble((n1+1)*n1))
          else
!           n1+1 = n2 or n2 > n1 + 2  all matrix elements are zero
            hbl(index) = x0
          endif
   20   continue
   10 continue
    if (.not. zprad) return
!     write kinetic energy integrals if requested
    write(6,"(//,5x,'radial kinetic energy matrix calculated analytically',/)")
    call symout(hbl,nmax+1)
    return
end

!#######################################################################
subroutine keint2(hbl,fke,rm2,dnorm,nlim,nmax,npnt,alf)

!     keint2 calculates analytic kinetic energy integrals over r2   #013
!     and moment of intertia integral for spherical oscillator functions
    use dvr3drjz_outp
    use dvr3drjz_oupb
    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(nlim) :: hbl
    double precision, dimension(nlim) :: rm2
    double precision, dimension(*) :: dnorm

    x0 = 0.0d0
    xp5 = 0.5d0
    x1 = 1.0d0

    gam = fke / (alf + xp5)
    do 10 n1=1,npnt
        gam = gam /(dble(n1)+alf+xp5)
   10 continue  
    fact = x1
    fn = x0
    sum = gam
    do 20 n1 = 1,nmax+1
    index = (n1 * (n1+1)) / 2
    do 30 n2 = n1,nmax+1
    rm2(index) = dnorm(n1) * dnorm(n2) * sum
    index = index + n2
   30 continue
    gam = (fn+alf+xp5) * gam
    fn = fn + x1
    fact = fn * fact
    sum = sum + gam / fact
   20 continue
    fact = alf * (alf + x1)
    do 40 index=1,nlim
    hbl(index) = - fact * rm2(index)
   40 continue
    index = 0
    fn = - x1
    do 50 n2=1,nmax+1
    fn = fn + x1
    index=index+n2
!     special case:  n1 = n2
    hbl(index) = hbl(index) + fke * (fn+fn+alf+1.5d0)
!     special case:  n2 = n1 + 1
    if (n2 > 1)&
        hbl(index-1) = hbl(index-1) + fke * sqrt(fn*(fn+alf+xp5))
   50 continue
    if (.not. zprad) return
!     write kinetic energy & inertia integrals if requested
    write(6,"(//5x,'radial kinetic energy matrix calculated analytically'/)")
    call symout(hbl,nmax+1)
    write(6,"(//5x,'moment of inertia matrix calculated analytically'/)")
    call symout(rm2,nmax+1)
    return
end

!##########################################################################
    subroutine k1k2(xk,hbl,bass,wt,npnt,nmax,nlim)

!     set up the transformed kinetic energy integrals,  t'(hbl) t
!                                                       ~  ~~~  ~
!     (note that the radial basis functions are already normalised)

    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(npnt,npnt) :: xk
    double precision, dimension(nlim) :: hbl
    double precision, dimension(0:nmax,npnt) :: bass
    double precision, dimension(npnt) :: wt

    xk = 0.0d0

    do 10 k=1,npnt
        do 20 kp=1,k
          wtkkp = wt(k)*wt(kp)
          do 30 m=0,nmax
            t = bass(m,k) * wtkkp
            do 40 mp=0,nmax
              in = max(m,mp) * (max(m,mp)+1)/2 + min(m,mp) + 1 
              xk(k,kp) = xk(k,kp) + (hbl(in) * t * bass(mp,kp))
   40       continue
   30     continue
        xk(kp,k)=xk(k,kp)
   20   continue
   10 continue
    return
end

!############################################################################
subroutine jacobi(nn,nn2,x,a,alf,bta,csa,tsa)

!     calculates zeros x(i) of the nn'th order jacobi polynomial
!     pn(alf,bta) for the segment (-1,1) & and corresponding weights
!     for gauss-jacobi integration. this routine uses a brute force
!     search for zeros due to GJ Harris (2001).
!     note that for our purposes, alf= bta= nu.


    implicit double precision(a-h,o-z)
    double precision, dimension(nn) :: x,a,b,c,xt
    x0 = 0.0d0
    x1 = 1.0d0
    x2 = 2.0d0
    x3 = 3.0d0
    x4 = 4.0d0
    eps = 1.0d-12

    fn= dble(nn)
    csa= x0
    c(1) = x0
    b(1) = x0
    do 10 i=2,nn
    xi= float(i)
    b(i) = x0
    c(i)= x4*(xi-x1)*(alf+xi-x1)*(bta+xi-x1)*(alf+bta+xi-x1)/&
            ((alf+bta+x2*xi-x1)*(alf+bta+x2*xi-x2)*&
            (alf+bta+x2*xi-x2)*(alf+bta+x2*xi-x3))
10 continue
    cc=tsa
    do 15 i=2,nn
    cc= cc*c(i)
15 continue

! step through the jacobi polynomial and 
! notes a first guess for the posititions of the zeros in the 
! array xt(ii). Zeros are found by looking for a change in sign.

    ii=0
    pm1=x1
    xxx=x1
30   continue
    call recur(p,dp,pn1,xxx,nn,alf,bta,b,c)

    if (pm1*p < x0) then
        pm1 = -pm1
        ii = ii +1
        xt(ii)=xxx
    endif

    if (ii == nn2) then
        do 40 i=1,nn2
        call recur(ptemp,dp,pn1,xt(i),nn,alf,bta,b,c)
40      continue
    else
        xxx=xxx-0.0001
        if (xxx > -0.002) goto 30
        write(6,*) "Incorrect number",ii-1," of zeros found in JACOBI"
        stop
    endif 

    do 20 i=1,nn2
    call root(xt(i),nn,alf,bta,dpn,pn1,b,c,eps)
    x(i)= xt(i)
    a(i)= cc/(dpn*pn1)
    csa= csa + a(i) + a(i)
20  continue
    if (2*nn2 /= nn) csa=csa-a(nn2)
    return
    end

!#####################################################################
    subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)

!     improves the approximate root x; in addition obtains
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x.

    implicit double precision(a-h,o-z)
    double precision, dimension(nn) :: b
    double precision, dimension(nn) :: c
   
    iter= 0
1   iter= iter + 1
    call recur(p,dp,pn1,x,nn,alf,bta,b,c)
    d = p/dp
    x = x - d
    if(abs(d) <= eps) goto 3 
    if(iter < 10) goto 1 
3   dpn= dp
    return
    end

!===================================================================
    subroutine recur(pn,dpn,pn1,x,nn,alf,bta,b,c)
    implicit double precision(a-h,o-z)
    double precision, dimension(nn) :: b
    double precision, dimension(nn) :: c
    
    x0 = 0.0d0
    x1 = 1.0d0
    x2 = 2.0d0

    p1= x1
    p= x + (alf-bta)/(alf + bta + x2)
    dp1= x0
    dp= x1
    do 1 j=2,nn
        q= (x - b(j))*p - c(j)*p1
        dq= (x - b(j))*dp + p - c(j)*dp1
        p1= p
        p= q
        dp1= dp
        dp= dq
1   continue
    pn= p
    dpn= dp
    pn1= p1
    return
    end

!################################################################
    subroutine allpts(xalf,walf,nang,nang2)

!     takes the points & weights generated by legpt for the half-range
!     and creates new arrays for the full-range (-1,+1).
    use dvr3drjz_size
    implicit double precision (a-h,o-y)
    double precision, dimension(idvr) :: xalf,walf

    scale = dble(max(1,idia))
    do 10 i=1,nang2
    walf(i) = sqrt(scale*walf(i))
    if (idia == 2) goto 10
    xalf(nang+1-i) = xalf(i)
    xalf(i)        =-xalf(i)
    walf(nang+1-i) = walf(i)
10 continue
  
    return
    end

!########################################################################
    subroutine asleg(pleg,lmax,x,nn2,kz,lincr)

!     calculate polynomials 1 to lmax at x = cos(theta) for m = 0 or 1,
!     using the routine of press et al, page 182,
!     for the polynomial part of associated legendre functions.
!     we have removed sin(theta)**2 for nu = 1.
!     this enables us to use jacobi integration with alf = bta = nu,
!     using routines derived from beidenharn and louck.

    implicit double precision (a-h,o-z)
    double precision, dimension(0:lmax,nn2) :: pleg
    double precision, dimension(nn2) :: x
    double precision, dimension(0:lmax) :: pnorm

    x1 = 1.0d0
    x2 = 2.0d0

    m = kz
    if (m<0) goto 999
    do 10 i=1,nn2
    pmm = x1
    fact = x1
    do 11 j=1,m
        pmm = -pmm * fact
        fact = fact + x2
   11 continue

    pleg(0,i) = pmm
    pmmp1= x(i)*(m+m+1)*pmm
    pleg(1,i)= pmmp1
    ll=1
    do 2 l= 2+m,lmax+lincr
    r2lm1 = dble(l+l-1)
    rlpmm1= dble(l+m-1)
    rlmm  = dble(l-m)
    pll= (x(i)*r2lm1*pmmp1 - rlpmm1*pmm)/rlmm
    pmm= pmmp1
    pmmp1= pll
    ll=ll+1
    pleg(ll,i)= pll
2     continue
10    continue

!     set up the normalisation constants
!     (pnorm)**2 = (2j + 1)/2   *   (j - k)! / (j + k)!
    jstart = m
    jj = -1
    do 13 j = jstart,lmax+lincr
    fact = x1
    do 12 i = j-m+1,j+m
        facti = dble(i)
        fact = fact * facti
   12 continue
    rj = dble(j)
    jj = jj + 1
    pnorm(jj) = (rj + rj + x1) / (fact + fact)
    pnorm(jj) = sqrt(pnorm(jj))
   13 continue
!     now normalise the polynomials
    do 14 i=1,nn2
        jj = -1
        do 15 j=jstart,lmax+lincr
            jj = jj + 1
            pleg(jj,i) = pleg(jj,i) * pnorm(jj)
15    continue
14 continue
    return
999  write(6,"(/,/,5x,'improper argument in subroutine asleg',/)")
    stop
    end

!###################################################################
    subroutine lmatrx(xlmatr,pleg,walf,kz,ipar0,nidvr,lincr,alf)

!     this subroutine sets up the lower triangle of the transformed
!     angular momentum matrix l(alpha,alpha')
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb
    implicit double precision (a-h,o-y), logical (z) 

    double precision, dimension(idvr,idvr) :: xlmatr
    double precision, dimension(0:maxleg,idvr) :: pleg
    double precision, dimension(idvr) :: walf

    x0 = 0.0d0

    if(zperp) then      
        jstart = 0
        iparam = jstart
    else
        jstart = kz
        iparam = lincr
    endif

    realkz = kz

    jdia=max(1,idia)
    jj0=-jdia
    if (.not.zperp .and. mod(jstart,jdia) /= ipar0) then
        jj0=jj0+1
        jstart=jstart+1
    endif
    do 10 k= 1,nidvr
        term =  walf(k)
        do 11 kp=k,nidvr
          sumj1=x0
          jj = jj0
          do 20 j=jstart,maxleg+iparam,jdia
            jj = jj + jdia
            if (zperp) then
               sumj1 = sumj1 + pleg(jj,k) * pleg(jj,kp) * &
                    dble((j+1+alf)*(j+alf))
            else
               sumj1 = sumj1 + pleg(jj,k) * pleg(jj,kp) * &
                    dble((j+1)*(j))
            endif
20     continue
        xlmatr(kp,k) = sumj1 * term * walf(kp)
        xlmatr(k,kp) = xlmatr(kp,k)
11   continue
10 continue   
    if (.not. zlmat) return
!     write xlmatr if requested
    write(6,"(5x,'L-matrix for kz =',i3,', ipar =',i2/)") kz,ipar0
    call sqout(xlmatr,nidvr)
    return
    end

!##########################################################################

    subroutine jhmain(xk1,xk2,xlmatr,r1,r2,xalf,r2m2t,kz)

!     this routine controls the dvr calculation in all cases except
!     symmetrised radau coordinates.
!     written by james henderson
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision(a-h,o-y),logical(z)

    double precision, dimension(npnt1) :: r1
    double precision, dimension(npnt2) :: r2
    double precision, dimension(idvr) :: xalf
    double precision, dimension(idvr,idvr) :: xlmatr
    double precision, dimension(npnt1,npnt1) :: xk1
    double precision, dimension(npnt2,npnt2) :: xk2
    double precision, dimension(npnt2,npnt2) :: r2m2t
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: eigs2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eig2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::ham3
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::cint
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::cintp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eigs2d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eval
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::evall
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vecs1l
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::vecs2l
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::vecs3l
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::phi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::psi 
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iv1
    INTEGER, ALLOCATABLE, DIMENSION(:) ::ndim2d
    INTEGER, ALLOCATABLE, DIMENSION(:) ::iv2
    INTEGER, ALLOCATABLE, DIMENSION(:,:) ::iv1l
    INTEGER, ALLOCATABLE, DIMENSION(:) ::iv2l
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ndim2l

    ALLOCATE(ham1(ndima,ndima),eigs1d(max2d),&
            eig1(ndima),iv1(ndimc,ndimb),vecs1d(max2d,ndima))

    term  = dble(jrot * jrot + jrot - (2 * kz * kz))
!     construct the one-dimensional hamiltonian matrix h1(npnta,npnta')
!     for each npntb and npntc, then solve by diagonalisation.
    icall = 0
    nsum = 0
    do 10 ione = 1,npntc
    nham2 = 0
    do 20 itwo = 1,npntb

    if (zr2r1) then
        call mkham1(ham1,xlmatr,ione,itwo,term,r1,r2,xalf,xk1,xk2,kz)
    else
        call mkham1(ham1,xlmatr,itwo,ione,term,r1,r2,xalf,xk1,xk2,kz)
    endif

!     diagonalise each block, saving the eigenvalues and vectors that
!     fall below the cut-off energy emax1.
    call diag(ham1,ndima,npnta,eig1)
    call cut1d(ham1,eig1,iv1(ione,itwo),eigs1d,vecs1d,nham2,icall)
   20 continue
    write(6,"(5x,' nham2 = ',i4)") nham2
    nsum = nsum + nham2
    if (ione == npntc) write(6,"(/5x,' sum = ',i5)") nsum
        if (nham2 > 0) then
!           dump the 1d eigenavlues & vectors to disk for each ione
            call outrow(eigs1d,nham2,ieigs1)
            do 2 ka=1,npnta
            if (zvec) call outrow(vecs1d(1,ka),nham2,iout2)
            call outrow(vecs1d(1,ka),nham2,ivecs1)
2    continue
    endif
10 continue

    if (zvec) write (iout1) iv1

    DEALLOCATE(ham1,eig1)
    ALLOCATE(ndim2d(ndimc),eigs2(max2d,ndimc),&
            ham2(max2d,max2d),eig2(max2d),iv2(ndimc),eigs2d(max3d))
      
!     now want to make the two-dimensional ham2(npntb,npntb',i,i'),
!     where i runs over the selected ham1(npnta,npnta') solutions.
    rewind ivecs1
    rewind ieigs1
    icall = 0
    low3d = 0
!     if zall then set zcut for convenience...
    if (zall) zcut = .true.
    do 31 ione = 1,npntc
    nham2 = 0
!     recall the size of ham2
    do 3 itwo = 1,npntb
    nham2 = nham2 + iv1(ione,itwo)
    3 continue
    ndim2d(ione) = nham2

    if ( nham2 > 0 ) then
        call mkham2(ham2,eigs1d,vecs1d,xk1,xk2,iv1,ione,nham2)
        if (.not. ztwod) then
           call diag(ham2,max2d,nham2,eig2)
        else
           call diag3d(ham2,nham2,eig2,kz)
           return
        endif

        if (zcut) then
          call cut2d(ham2,eig2,iv2(ione),nham2,low3d,icall)
        else
          do 4 ii = 1,nham2
            call outrow(ham2(1,ii),nham2,intvec)
            eigs2(ii,ione) = eig2(ii)
4     continue
          low3d = max3d
        endif

    endif

    if (ione == npntc .and. zcut)&
    & write(6,"(/i14,' eigenvalues selected below ',d20.10)")  low3d,emax2

31 continue

    if (.not. zcut) call choose(eigs2,ndim2d,ham2,iv2,low3d)

    if (zvec) write (iout1) iv2

!      call timer

    nham3 = low3d

!     save the required bits to disk if zdiag = .false.
!     first open the required disk file
    if (.not. zdiag) then
        if (zpseg==.true.) then 
        open(unit=idiag1,form='unformatted',recordtype='segmented')
        open(unit=idiag2,form='unformatted',recordtype='segmented')
    else 
        open(unit=idiag1,form='unformatted')
        open(unit=idiag2,form='unformatted')
    end if 
        write(idiag1) npnta,npntb,npntc,max2d,nham3,zcut
        call ioutro(iv1,ndimb*ndimc,idiag1)
        call ioutro(iv2,ndimc,idiag1)

!       need also the 2-d eigenvalues to build the final hamiltonian
        if (zcut) then
          rewind ieigs2
          do 55 i1 = 1,npntc
          iv = iv2(i1)
          if (iv > 0) call getrow(eig2,iv,ieigs2)
          if (iv > 0) call outrow(eig2,iv,idiag1)
   55     continue
        else
          call outrow(eigs2,max2d*ndimc,idiag1)
        endif

    endif
    DEALLOCATE(eigs1d,eig2)
    ALLOCATE(hband(max2d,nham3),&
               cint(ndima*ndimb,max2d),cintp(ndima*ndimb,max2d))

    call mkham3(nham3,ham2,cint,cintp,iv1,iv2,ndim2d,&
                  vecs1d,xk1,xk2,eigs2d,eigs2,hband,xlmatr,r1,r2,&
                  r2m2t,term)

    DEALLOCATE(ham2,cint,cintp,vecs1d,hband,iv1)
    ALLOCATE(ham3(nham3,nham3))

    call loadh(ham3,nham3,iv2,ndim2d,eigs2d,eigs2)

    DEALLOCATE(iv2,ndim2d,eigs2d,eigs2)
    ALLOCATE(eval(nham3))

    call diag3d(ham3,nham3,eval,kz)

    DEALLOCATE(ham3,eval)

!.....finally, can compute the actual wavefunction amplitude at
!.....the grid points if needed.
    if (ztran) then
        ALLOCATE(evall(neval),&
                iv1l(ndimc,ndimb),iv2l(ndimc),ndim2l(ndimc),&
                vecs1l(max2d,ndima),vecs2l(max2d,max2d),vecs3l(nham3),&
                phi(nham3,ndima,ndimb),psi(idvr,npnt1,npnt2))

        call trans(iv1l,iv2l,ndim2l,vecs1l,&
                vecs2l,vecs3l,phi,psi,evall,nham3)

        DEALLOCATE(iv1l,iv2l,ndim2l,vecs1l,vecs2l,vecs3l,phi,psi,evall)
    endif

    return
    end

!#########################################################################
    subroutine mkham1(ham1,xlmatr,i1,i2,term,r1,r2,xalf,xk1,xk2,kz)
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_split1
    use dvr3drjz_split2
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(ndima,ndima) :: ham1
    double precision, dimension(idvr,idvr) :: xlmatr
    double precision, dimension(npnt1) :: r1
    double precision, dimension(npnt2) :: r2
    double precision, dimension(idvr) :: xalf
    double precision, dimension(npnt1,npnt1) :: xk1
    double precision, dimension(npnt2,npnt2) :: xk2

    x0 = 0.0d0
    xp5 = 0.50d0
    x1 = 1.0d0

!     zero ham1
    ham1 = x0
!     zero rotational excitation term for j=0 cases
    wterm = x0
!.....theta first
    if (ztheta) then

        w1gama = xp5 / (r1(i1)*r1(i1)*ur1)
        w2beta = xp5 / (r2(i2)*r2(i2)*ur2)
        wsum = w1gama + w2beta

        if (jrot > 0) then
            if (zembed) then
!              have term * r2**(-2) term
               wterm = term * w2beta * ur2/urr2
            else
!              have term * r1**(-2) term
               wterm = term * w1gama * ur1/urr1
            endif
!     extra NBO term if vib mass /= rot mass
            if (kz > 0) then
               s1 = ur1/urr1-x1
               s2 = ur2/urr2-x1
               w3 = dble(kz*kz) * (s1*w1gama + s2*w2beta)*xp0
            endif
        endif
    endif

    do 10 k = 1,npnta
    if (ztheta) then
        call potv(v,r1(i1),r2(i2),xalf(k))
    else
        if (zr2r1)then
            call potv(v,r1(i2),r2(k),xalf(i1))
        else
            call potv(v,r1(k),r2(i1),xalf(i2))
        endif
    endif

    if (ztheta) then
        if (kz > 0) v = v + w3/(x1-xalf(k)**2)
        ham1(k,k) = v + wterm
        do 20 kp= 1,k
        ham1(k,kp) = ham1(k,kp) + xlmatr(k,kp)*wsum
   20    continue    

    else
        if (jrot > 0) then
            if (zembed) then
              if (zr2r1)then
                 wterm = (term * xp5) / (r2(k)*r2(k)*ur2)
              else
                 wterm = (term * xp5) / (r2(i1)*r2(i1)*ur2)
              endif
            else
              if (zr2r1)then
                 wterm = (term * xp5) / (r1(i2)*r1(i2)*ur1)
              else
                 wterm = (term * xp5) / (r1(k)*r1(k)*ur1)
              endif
            endif
        endif
        ham1(k,k) = v + wterm

        if (zr2r1) then
            do 30 kp= 1,k
            ham1(k,kp) = ham1(k,kp) + xk2(k,kp)
   30       continue
        else
            do 40 kp= 1,k
            ham1(k,kp) = ham1(k,kp) + xk1(k,kp)
   40       continue
        endif
    endif
   10 continue

    return
    end

!####################################################################
    subroutine mkham2(ham2,eigs1d,vecs1d,xk1,xk2,iv1,ione,nham2)
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)


    double precision, dimension(max2d,max2d) :: ham2
    double precision, dimension(npnt2,npnt2) :: xk2
    dimension iv1(ndimc,ndimb)
    double precision, dimension(max2d) :: eigs1d
    double precision, dimension(max2d,ndima) :: vecs1d
    double precision, dimension(npnt1,npnt1) :: xk1

!     zero ham2
    ham2 = 0.0d0

    do 10 i3= 1,npnta
        call getrow(vecs1d(1,i3),nham2,ivecs1)
        do 20 j = 1,nham2
           do 30 k = 1,j
             ham2(j,k) = ham2(j,k) + vecs1d(j,i3)*vecs1d(k,i3)
   30      continue
   20    continue
   10 continue
!     must now multiply by xk1 or xk2
    ivbsm = 0
    do 50 itwo = 1,npntb
        ivbpsm = 0
        ivb = iv1(ione,itwo)
        do 60 itwop = 1,itwo
          if (ztheta) then
            if (zr2r1) then
              xkterm = xk2(itwo,itwop)
            else
             xkterm = xk1(itwo,itwop)
            endif
          else
           if (zr2r1) then
           xkterm = xk1(itwo,itwop)
        else
           xkterm = xk2(itwo,itwop)
        endif
    endif
    ivbp = iv1(ione,itwop)
    do 70 j  = 1,ivb
        ind1 = ivbsm + j
        do 77 jp = 1,ivbp
          ind2 = ivbpsm + jp
          ham2(ind1,ind2) = ham2(ind1,ind2) * xkterm
   77   continue
   70 continue
    ivbpsm = ivbpsm + ivbp
   60 continue
    ivbsm = ivbsm + ivb
   50 continue

!     now add the 1-d eigenvalues along the diagonal
    call getrow(eigs1d,nham2,ieigs1)
    do 80 nn = 1,nham2
        ham2(nn,nn) = ham2(nn,nn) + eigs1d(nn)
   80 continue
    return
    end

!################################################################################
    subroutine mkham3(nham3,ham2,cint,cintp,iv1,iv2,ndim2d,&
                        vecs1d,xk1,xk2,eigs2d,eigs2,hband,xlmatr,&
                        r1,r2,r2m2t,term)

!     build the final 3-d hamiltonian matrix.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_split1
    use dvr3drjz_split2
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(ndima*ndimb,max2d) :: cint
    double precision, dimension(max2d,max2d) :: ham2
    dimension iv1(ndimc,ndimb),iv2(ndimc),ndim2d(npntc)
    double precision, dimension(max2d,ndima) :: vecs1d
    double precision, dimension(ndima*ndimb,max2d) :: cintp
    double precision, dimension(npnt1,npnt1) :: xk1
    double precision, dimension(nham3) :: eigs2d
    double precision, dimension(npnt2,npnt2) :: xk2
    double precision, dimension(max2d,nham3) :: hband
    double precision, dimension(max2d,ndimc) :: eigs2
    double precision, dimension(idvr,idvr) :: xlmatr
    double precision, dimension(npnt1) :: r1
    double precision, dimension(npnt2) :: r2
    double precision, dimension(nham3) :: work3
    double precision, dimension(npnt2,npnt2) :: r2m2t

    xp5 = 0.50d0

!     if zdiag = .false. want eigs2d now
    if (.not. zdiag) then
        rewind ieigs2
        neig = 0
        do 184 ione = 1,npntc
          iv = iv2(ione)
          if (.not. zcut) then
            do 19 ii = 1,iv
              neig = neig + 1
              eigs2d(neig) = eigs2(ii,ione)
   19       continue
          endif
          if (iv > 0 .and. zcut) call getrow(eigs2d,iv,ieigs2)
  184   continue
    endif

!     first do the intermediate transformation
    ncint = npnta*npntb
    ndimt = ndima*ndimb
    rewind ivecs2
    rewind ivecs1
    do 10 ione = 1,npntc
!       recall the size of the 2-d vectors
        nham2=ndim2d(ione)
        if (nham2 == 0) goto 10
        cint = 0.0d0
!     bring back the 1-d vectors for each i1
        do 23 kk = 1,npnta
           if (iv2(ione) > 0) then
             call getrow(vecs1d(1,kk),nham2,ivecs1)
           else
             read (ivecs1)
           endif
23    continue
        do 20 j = 1,iv2(ione)
!     bring back the 2-d vectors for each npntc
           call getrow(ham2(1,j),nham2,ivecs2)
           ind2 = 0
           do 30 k = 1,npnta
             ind1 = 0
             do 40 itwo = 1,npntb
               ind2 = ind2 + 1
               do 50 i = 1,iv1(ione,itwo)
                 ind1 = ind1 + 1
                 cint(ind2,j) = cint(ind2,j) + ham2(ind1,j)*vecs1d(ind1,k)
   50          continue
   40        continue
   30      continue
           if (.not. ztheta) then
             in2 = 0
             do 37 ia = 1,npnta
               do 47 ib = 1,npntb
                 if (zr2r1) then
                    ii1 = ib
                    ii2 = ia
                 else
                    ii1 = ia
                    ii2 = ib
                 endif
                 w1 = xp5 / (r1(ii1)*r1(ii1)*ur1)
                 w2 = xp5 / (r2(ii2)*r2(ii2)*ur2)
                 wsum = w1 + w2
                 in2 = in2 + 1
                 cint(in2,j) = cint(in2,j)*sqrt(wsum)
   47          continue
   37        continue
           endif
   20    continue
!     store cint on disk for each npntc
        if (iv2(ione) > 0) call outrow(cint,ndimt*iv2(ione),ivint)
   10 continue

!     now do the second part of the transformation
!      endfile ivint


    rewind ivint
    if (ztheta) then
        length = 0
        do 61 ione = 1,npntc
!        set the hband to zero
        hband = 0.0d0
        ivsm = 0
        nham2 = ndim2d(ione)
        if (nham2 == 0) goto 61
        iv = iv2(ione)
        if (iv > 0) call getrow(cint,ndimt*iv,ivint)
        rewind ivint
        ivpsm = 0
        do 51 ionep = 1,ione
           if (ione == ionep) then
             if (zquad2) then
               rm2t = 0.0d0
             else
               d1r2 = r2m2t(ione,ionep)
               d2r2 = xp5 / (r2(ione)*r2(ione)*ur2)
               rm2t = d1r2 - d2r2
             endif
           else
             if (zquad2) then
               rm2t = 0.0d0
             else
               rm2t = r2m2t(ione,ionep)
             endif
           endif
           if (zr2r1) then
             xkterm = xk1(ione,ionep)
           else
             xkterm = xk2(ione,ionep)
             if(jrot > 0 .and. .not. zquad2 .and. zembed) xkterm = xkterm + term*rm2t
           endif
           ivp = iv2(ionep)
           if (ivp > 0) call getrow(cintp,ndimt*ivp,ivint)
           do 41 j=1,iv
             ind1 = ivsm + j
             do 31 jp =1,ivp
               ind2 = ivpsm + jp
               ind3 = 0
               do 81 k = 1,ncint
                  ind3 = ind3 + 1
                  hband(ind1,ind2) = hband(ind1,ind2)&
                                     + cint(ind3,j)*cintp(ind3,jp)*xkterm
   81          continue
               ind3 = 0
               do 83 na = 1,npnta
                  do 82 k = 1,npntb
                    ind3 = ind3 + 1
                    ind4 = k
                    do 84 nap = 1,npnta
                       hband(ind1,ind2) = hband(ind1,ind2) + cint(ind3,j)*cintp(ind4,jp)&
                                           * xlmatr(na,nap) * rm2t
                       ind4 = ind4 + npntb
   84               continue
   82             continue
   83          continue
   31        continue
   41      continue
        ivpsm = ivpsm + ivp
   51    continue
        if (zdiag) then
           do 71 jj=1,ind2
             if (iv > 0) call outrow(hband(1,jj),iv,iband)
   71      continue
           ndim2d(ione)=ind2
        else
           do 93 isave = 1,iv
             length = length + 1
             do 94 jsave = 1,ind2
               work3(jsave)= hband(isave,jsave)
   94        continue
             work3(length) = work3(length) + eigs2d(length)
             call outrow(work3,length,idiag2)
   93      continue
        endif
   61 continue

    else
        length = 0
        do 66 ione = 1,npntc
!         set the hband to zero
        hband=0.0d0
        ivsm = 0
        if (ndim2d(ione) == 0) goto 66
        iv = iv2(ione)
        if (iv > 0) call getrow(cint,ndimt*iv,ivint)
        rewind ivint
        ivpsm = 0
        do 56 ionep = 1,ione
            xkterm = xlmatr(ione,ionep)
            ivp = iv2(ionep)
            if (ivp > 0) call getrow(cintp,ndimt*ivp,ivint)
            do 46 j=1,iv
              ind1 = ivsm + j
              do 36 jp =1,ivp
                ind2 = ivpsm + jp
                ind3 = 0
                do 86 k = 1,ncint
                  ind3 = ind3 + 1
                  hband(ind1,ind2) = hband(ind1,ind2)&
                               + cint(ind3,j)*cintp(ind3,jp)*xkterm
   86           continue
   36         continue
   46       continue
            ivpsm = ivpsm + ivp
   56   continue
        if (zdiag) then
            do 76 jj=1,ind2
              if (iv > 0) call outrow(hband(1,jj),iv,iband)
   76       continue
            ndim2d(ione)=ind2
        else
            do 95 isave = 1,iv
              length = length + 1
              do 96 jsave = 1,ind2
                work3(jsave)= hband(isave,jsave)
   96         continue
            work3(length) = work3(length) + eigs2d(length)
            call outrow(work3,length,idiag2)
   95       continue
        endif
   66  continue
    endif

    if (.not. zdiag) then
        write(6,"(/5x,'hamiltonian written to disk - not diagonalised')")
        write(6,"(/5x,'hamiltonian bands written to stream idiag =',i4/)") idiag2
        stop
    endif
    return
    end

    subroutine loadh(ham3,nham3,iv2,ndim2d,eigs2d,eigs2)

!     load the final 3-d hamiltonian matrix.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)


    dimension iv2(ndimc),ndim2d(npntc)
    double precision, dimension(nham3,nham3) :: ham3
    double precision, dimension(nham3) :: eigs2d
    double precision, dimension(max2d,ndimc) :: eigs2

    ham3 = 0.0d0
    rewind iband
    ivsm = 0
    do 62 ione = 1,npntc
        if (ndim2d(ione) == 0) goto 62
        iv = iv2(ione)
        do 72 jj=1,ndim2d(ione)
          if (iv > 0) call getrow(ham3(ivsm + 1,jj),iv,iband)
   72   continue
        ivsm = ivsm + iv
   62 continue

!     now put the 2-d eigensolutions along the diagonal

    rewind ieigs2
    nn = 0
    do 183 ione = 1,npntc
    iv = iv2(ione)
    if (.not. zcut) then
       do 9 ii = 1,iv
          eigs2d(ii) = eigs2(ii,ione)
    9  continue
       endif
       if (iv > 0 .and. zcut) call getrow(eigs2d,iv,ieigs2)
       do 185 i = 1,iv
         nn = nn + 1
         ham3(nn,nn) = ham3(nn,nn) + eigs2d(i)
  185  continue
  183 continue

    if (zpham) call wrtham(ham3,nham3)

    return
    end

!########################################################################
!      subroutine diag(ham,maxham,nham,eig)
!        ! Diag Divide And Conquer
!        ! new routine uses the faster dysevd lapack routine
!        ! but does require an extra N**2 workspace,
!        ! which for very large runs may be restrictive
!        implicit none
!        ! external vars
!        integer :: maxham,nham
!        double precision :: ham(maxham,nham)
!        double precision :: eig(nham)
!        ! internal vars
!        integer :: ifail,lwork,liwork
!        double precision :: work(1 + 6*nham + 2*nham**2)
!        integer :: iwork(3 + 5*nham)

!        lwork = 1 + 6*nham + 2*nham**2
!        liwork = 3 + 5*nham
!        ifail = 0
        
!        call dsyevd('V','L',nham,ham,maxham,eig,work,lwork,iwork,liwork,ifail)
        
!        if (ifail /= 0) write(6,100) ifail
!        return
!100     format(' diagonalisation has failed with, ifail=',i3)
!      end subroutine diag
      
    SUBROUTINE diag(ham,maxham,nham,eig)

!     diagonalise the appropriate hamiltonian matrice
    use dvr3drjz_outp
    implicit real*8 (a-h,o-y), logical (z)

    real*8, dimension(maxham,nham) :: ham
    real*8, dimension(nham) :: eig
    real*8, dimension(maxham*3) :: work
      
    ifail=0
    nnham=maxham*3
    call dsyev ('V','L',nham,ham,maxham,eig,work,nnham,ifail)

    if (ifail /= 0) &
    &write(6,"(' diagonalisation has failed with, ifail=',i3)") ifail
    return
    END SUBROUTINE diag
!########################################################################
    subroutine diag3d(ham3,nham3,eval,kz)
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)
    double precision, dimension(nham3,nham3) :: ham3
    double precision, dimension(nham3) :: eval
    double precision, dimension(nham3) :: evalcm

!     autocm converts atomic units (hartree) to cm-1.
    autocm = 2.19474624d+05
    x0 = 0.0d0
    if (zrot) then
        write(6,"(/5x,'Solutions with J =',i3,' k =',i3)") jrot,kz
        if (idia == -2) then
            if (ipar == 0)  write(6,"(/5x,'even parity solutions')")
            if (ipar == 1)  write(6,"(/5x,'odd parity solutions')")
        endif
    endif
!     we need this in special cases (nham3 < neval)
    meval=min(neval,nham3)
    if (zvec) write(iout2) meval

!     diagonalise  the hamiltonian. the vectors are to overwrite the
!     hamiltonian.
    call diag(ham3,nham3,nham3,eval)

!     print eigenvalues in atomic untis & wavenumbers

    if (.not.zpmin) then
        write(6,"(/5x,'lowest',i4,' eigenvalues in hartrees:',/)") meval
        write(6,"(5d24.12/)") (eval(i),i=1,meval)
    endif
    if (zpfun) then
        ip=jrot-kmin
        if (jrot /= 0 .and. ip /= 1) goto 10
        jdia=max(0,idia)
        jpar=min(jdia,ipar)
        isym=abs(min(0,idia))
        if (ipar == 1) isym=-isym
        write(ilev,"(6i6)") jrot,ip,jdia,jpar,isym,meval
        write(ilev,"(4d20.12)") (eval(i),i=1,meval)
    endif
   10 continue

!     save the eigenvalues if needed
    if (zvec) write (iout2) (eval(i),i=1,meval)
    do 20 i=1,meval
    evalcm(i) = eval(i) * autocm
   20 continue
    write(6,"(//5x,'Lowest',i4,' eigenvalues in wavenumbers:'/)") meval
    write(6,"(5d24.12/)") (evalcm(i),i=1,meval)

!     if requested print the eigenvectors
    if (zpvec) then
        write(6,"(//'  Eigenvectors'/)")
        do 40 i=1,meval
        write(6,"(1x,10f13.7)") (ham3(j,i),j=1,nham3)
   40    continue
    endif
!     write the final vectors to disk if required
    if (zvec) then
        do 60 l = 1,meval
        call outrow(ham3(1,l),nham3,iout2)
   60    continue
    endif
    if (jrot /= 0) return
    if (abs(idia) == 2 .and. ipar == 1) then
        ii=1
        ezero=x0
        read(5,'(f20.0)',end=55) ezero
   55    continue
    else
        ezero=evalcm(1)
        ii=2
    endif
    write(6,"(//'  Band origins in wavenumbers:'/)")
    write(6,"(5f15.6/)") (evalcm(i)-ezero,i=ii,meval)

    return
    end

!########################################################################
    subroutine choose(eigs2,ndim2d,ham2,iv2,low3d)

!     this routine chooses the max3d lowest eigenvalues from eigs2.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(max2d,ndimc) :: eigs2
    dimension iv2(ndimc)
    dimension ndim2d(ndimc)
    double precision, dimension(max2d,max2d) :: ham2

    autocm = 2.19474624d+05

    eigmin = eigs2(1,1)
    nhamsm = 0
    do 160 i=1,npntc
        eigmin = min(eigmin,eigs2(1,i))
        nhamsm = nhamsm + ndim2d(i)
        iv2(i) = 1
160 continue
    if (nhamsm < max3d) low3d = nhamsm
    ipt = 1
    do 200 n=1,low3d
210 if(iv2(ipt) <= ndim2d(ipt)) then
        eigvib = eigs2(iv2(ipt),ipt)
        jpt = ipt
    else
        ipt = ipt + 1
        goto 210
    endif
    do 220 j=ipt+1,npntc
    if(iv2(j) > ndim2d(j)) goto 220
    if(eigs2(iv2(j),j) >= eigvib) goto 220
    eigvib = eigs2(iv2(j),j)
    jpt = j
220 continue
!     keep the eigenvalue:
    iv2(jpt) = iv2(jpt) + 1
200 continue
!     store the number of eigenvalues selected for each alpha
    iv2(1) = iv2(1) - 1
    ivib = iv2(1)
    write(6,"(//5x,'Selection outcome for the final, contracted basis:'/)")
    write(6,"(5x,'npntc =',i3, ',', ' no. of eigenvectors =',i3 )")  1,iv2(1)
    do 230 i=2,npntc
       iv2(i) = iv2(i) - 1
       write(6,"(5x,'npntc =',i3, ',', ' no. of eigenvectors =',i3 )")  i,iv2(i)
       ivib = max(ivib,iv2(i))
  230 continue

    write(6,"(/i14,' eigenvalues selected from ',d20.10,' to',d20.10,' hartrees')") &
        &low3d,eigmin,eigvib
    write(6,"(/i14,' eigenvalues selected from ',d20.10,' to',d20.10,' cm-1')") &
        & low3d,eigmin*autocm,eigvib*autocm



    if (zp2d) then
        write(6,"(//5x,'2d eigenvalues in wavenumbers:'/)")
        do 32 jone = 1,npntc
        ivj = iv2(jone)
        write(6,"(5d24.12/)") (eigs2(i,jone)*autocm,i=1,ivj)
        write(6,"(/5x,'------------------------------'/)")
   32    continue
    endif      

!     save the vectors chosen
    rewind intvec
    do 31 ione = 1,npntc
    ivm = iv2(ione)
    nham2 = ndim2d(ione)
    do 7 i = 1,nham2
        call getrow(ham2(1,i),nham2,intvec)
    7 continue
    do 30 ind = 1,ivm
    if (zvec) then
        if (nham2 > 0) call outrow(ham2(1,ind),nham2,iout2)
    endif
        if (nham2 > 0) call outrow(ham2(1,ind),nham2,ivecs2)
   30 continue
   31 continue

    return
    end

!#####################################################################
    subroutine cut1d(ham1,eig1,ivn,eigs1d,vecs1d,nham2,icall)

!     this routine selects all the eigenvalues that are lower than the
!     the cut-off energy emax1, which is user-supplied in wavenumbers.
!     these eigenvalues & their corresponding vectors are then saved
!     to disk.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)


    double precision, dimension(ndima,ndima) :: ham1
    double precision, dimension(ndima) :: eig1
    double precision, dimension(max2d) :: eigs1d
    double precision, dimension(max2d,ndima) :: vecs1d

    autocm = 2.19474624d+05
!     change emax1 to hartree for the selection
    emaxau=emax1/autocm

    icall = icall + 1
    if (.not. zall) then
        ivn = 0
        do 10 n=1,npnta
        if (eig1(n) > emaxau) goto 20
        ivn = ivn + 1
10   continue
20   iv = ivn
    else
        ivn = npnta
        iv  = ivn
    endif

    if (zp1d) then
        if (icall == 1) &
            &write(6,"(//5x,'1d eigenvalues in wavenumbers:'/)")
        write(6,"(5d24.12/)") (eig1(i)*autocm,i=1,npnta)
        write(6,"(/5x,'------------------------------'/)")
    endif

!     save the vectors and eigenvalues (overwrite for each gamma).
    do 30 i=1,iv
        nham2 = nham2 + 1

    if (nham2 > max2d .and. .not. zall)  then
        write(6,"(//6x,'**** core exceeded: reduce cut-off emax1 ****')")
        stop
    endif

        eigs1d(nham2)   = eig1(i)
    do 31 j=1,npnta
        vecs1d(nham2,j) = ham1(j,i)
   31 continue
   30 continue

    return
    end

!####################################################################
    subroutine cut2d(ham2,eig2,ivm,nham2,low3d,icall)

!     this routine selects all the eigenvalues that are lower than the
!     the cut-off energy emax2, which is user-supplied in wavenumbers.
!     these eigenvalues & their corresponding vectors are then saved
!     to disk.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(max2d,max2d) :: ham2
    double precision, dimension(nham2) :: eig2

    autocm = 2.19474624d+05
!    change emax2 to hartree for the selection
    emaxau=emax2/autocm

    icall = icall + 1
    if (.not. zall) then
        ivm = 0
        do 10 n=1,nham2
        if(eig2(n) > emaxau) goto 20
        ivm = ivm + 1
   10    continue
   20    low3d = low3d + ivm
    else
        ivm = nham2
        low3d = low3d + ivm
    endif

    if (low3d > max3d .and. .not. zall)  then
        write(6,"(//6x,'**** core exceeded: reduce cut-off emax2 ****')")
        stop
    endif

    if (zp2d) then
        if (icall == 1) write(6,"(//5x,'2d eigenvalues in wavenumbers:'/)")
        write(6,"(5d24.12/)") (eig2(i)*autocm,i=1,nham2)
        write(6,"(/5x,'------------------------------'/)")
    endif

!     save the vectors and eigenvalues
    if (ivm > 0) call outrow(eig2,ivm,ieigs2)
    do 30 ind = 1,ivm
    if (zvec) call outrow(ham2(1,ind),nham2,iout2)
    call outrow(ham2(1,ind),nham2,ivecs2)
   30 continue
    return
    end

!########################################################################
    subroutine nfmain(hr,htheta,r,theta,kz)

!     this routine controls the dvr calculation in the case of
!     symmetrised radau coordinates.
!     written by nic fulton, feb 1993.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb
    implicit double precision(a-h,o-y),logical(z)


    double precision, dimension(npnt,npnt) :: hr
    double precision, dimension(nalf,nalf) :: htheta
    double precision, dimension(npnt) :: r
    double precision, dimension(nalf) :: theta
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::ham2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eig2,eigtmp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::vecs2d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eigs2d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::ham3
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eig3
    INTEGER, ALLOCATABLE, DIMENSION(:,:) ::iv2
    INTEGER, ALLOCATABLE, DIMENSION(:) ::nv2
    x4 = 4.0d0
    x8 = 8.0d0
    x16 = 1.6d1

    ALLOCATE(ham2(max2d,max2d),eig2(max2d),iv2(2,nalf),& 
            vecs2d(max2d,max3d),eigs2d(max3d),nv2(nalf))
    if (.not.zcut) ALLOCATE(eigtmp(nalf*max2d))

! Schwenke rotational corrections (times 1.1 !!!)

!      asw=0.98205515d-8*(1.10d0)*(0.9d0)
    asw=0.98205515d-8*(1.10d0)*(0.06d0)!(0.194d0)!*(0.06d0)!
!      bsw=7.3757478d-9*(1.10d0)*(0.9d0)
    bsw=7.3757478d-9*(1.10d0)*(0.06d0)!(0.194d0)!*(0.06d0)!
!      csw=1.7082113d-8*(1.10d0)*(-0.55d0)
    csw=1.7082113d-8*(1.10d0)*(-0.21d0)!(0.14d0)!*(-0.21d0)!
!       asw=1.080260665d-8
!       bsw=8.11332258d-9
!       csw=1.87903243d-8

    xtscw=(dble(jrot * jrot + jrot)*(asw+bsw)*0.5d0+dble(kz*kz)& 
        *(csw-(bsw+asw)*0.5d0))*xp2
    if (jrot /= 0) then
        term  = dble(jrot * jrot + jrot - kz * kz) / x8
        term2 = dble(jrot * jrot + jrot - 3 * kz * kz) / x4
        if (abs(kz) == 1) then
            term3 = dble(jrot * jrot + jrot) / x16
            if (kmin >= 1 .and. zrot) term3 = -term3
        endif
    endif

! no need for the 1d diagonalisations as there is no possiblity of
! truncation as the symmetry would be broken.
    write(6,"(5x,'starting the 2d calculations.')")
    call nftim('beginning of 2d loop')
    do 30 igamma = 1,nalf
  
        nv2(igamma) = (npnt*(npnt+1-ipar*2))/2
        nham2 = nv2(igamma)
        call blc2d1(theta(igamma),r,hr,ham2,nham2,&
                        term,term2,term3,xtscw,kz)

        call diag(ham2,nham2,nham2,eig2)

        if (.not. zcut) then
          call choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,&
                    vecs2d,nv2,eigtmp)

          if(igamma == nalf .and. .not. zpmin)&
              write(6,"(5x,'for gamma = ',i2,' selected ',i3,' energies.')") &
                &(itmp, iv2(2,itmp),itmp=1,nalf)
        else
          call cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)
          if (.not. zpmin) write(6,"(5x,'for gamma = ',i2,' selected ',i3,' energies.')") &
            & igamma, iv2(2,igamma)
        endif
30    continue
    call testiv(iv2,nbass)

    nham3 = iv2(1,nalf) + iv2(2,nalf) - 1
    if (.not. zpmin) write(6,"(25x,'total = ',i5)") nham3
    call nftim('end of 2d loop')
    write(6,"(5x,'building the 3d hamiltonian.')")

    if (.not.zcut) DEALLOCATE(eigtmp)
    ALLOCATE(ham3(nham3,nham3))

    call bloc3d(htheta,ham3,eigs2d,vecs2d,iv2,nv2,ham2,nham3,r)

    DEALLOCATE(ham2,eig2,eigs2d,nv2)
 
    call nftim('end of 3d ham building')
    write(6,"(5x,'diagonalising the 3d hamiltonian.')")
    ALLOCATE(eig3(nham3))

    call diag3d(ham3,nham3,eig3,kz)

    call nftim('end of diagonalising 3d')

    if (ztran) call transr(iv2,vecs2d,ham3,eig3,nham3,nbass)

    DEALLOCATE(iv2,vecs2d,ham3,eig3)
    return
    end

!    ***********************************************************************

    subroutine blc2d1(xcos,r,hr,ham2,nham2,term,term2,term3,xtscw,kz)
    use dvr3drjz_size
    use dvr3drjz_split1
    use dvr3drjz_oupb

    implicit double precision(a-h,o-y),logical(z)

    double precision, dimension(npnt,npnt) :: hr
    double precision, dimension(npnt,npnt) :: hrpb
    double precision, dimension(npnt) :: r
    double precision, dimension(nham2,nham2) :: ham2
    xp5 = 0.5d0
    x1 = 1.0d0
    factr2 = sqrt(xp5)

    ham2 = 0.0d0
    xsen2=1.d0-xcos**2

! Schwenke radial corrections (it looks like a mass correction)
! Kolya
    hrpb=hr*(1.d0-xp1*0.4150706813595d-3*1.25d0)
!      hrpb=hr
! I'm cancelling the rotational correction!!!!!!!!!!!!
!      xtscw=0.d0


!     allow for j > 0 case
    if (jrot /= 0) then
        fact =  term + term2 / (x1 - xcos)
        if (kz == 1) fact = fact + term3 * (x1 + xcos)/(x1-xcos)
        ia = 0
        do 15 ibeta=1,npnt
          do 25 ialpha=1,ibeta-ipar
            walpha = xp5 / (r(ialpha)*r(ialpha))
            wbeta = xp5 / (r(ibeta)*r(ibeta))
! Attenzione: ho inserito Schwenke rotational corrections qui!!!!!!
!(il segno meno perche' esse sono effettivamente negative, almeno spero)
            wsum = (walpha + wbeta) * fact / urr1 - xtscw

! NBO term in Radau case (wsum2)

        wsum2=kz*kz*(walpha+wbeta)*(1.d0/ur1-1.d0/urr1)*xp0/xsen2

       ham2(ialpha+ia,ialpha+ia)=ham2(ialpha+ia,ialpha+ia)+wsum-wsum2

 25       continue
          ia=ia+ibeta-ipar
          if (ipar == 0) ham2(ia,ia) = ham2(ia,ia) + wsum - wsum2
 15     continue
    endif
    q=1.0d0
    if(ipar == 1) q=-1.0d0
    iap=0
    do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
                if(ibeta == ibetap)&
                    ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)+hrpb(ialphp,ialpha)
                if(ibeta == ialphp)&
                    ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)+q*hrpb(ibetap,ialpha)
                if(ialpha == ibetap)&
                    ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)+q*hrpb(ialphp,ibeta)
                if(ialpha == ialphp)&
                    ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)+hrpb(ibetap,ibeta)
                if(ialpha == ialphp .and. ibeta == ibetap) then
                    call potv(v,r(ialpha),r(ibeta),xcos)
                    ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)+v
                endif
                if(ialpha == ibetap .and. ibeta == ialphp) then
                    call potv(v,r(ibeta),r(ialpha),xcos)
                    ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)+q*v
                endif
                if(ialphp == ibetap) ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)*factr2
                if(ialpha == ibeta) ham2(ialphp+iap,ialpha+ia)=&
                    ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
        ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
    return
end

!     ***********************************************************************

    subroutine choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d,nv2,eigtmp)

!     this routine chooses the max3d lowest eigenvalues from eigs2.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)


    double precision, dimension(nham2) :: eig2
    double precision, dimension(nham2,nham2) :: ham2
    dimension iv2(2,nalf)
    double precision, dimension(max3d) :: eigs2d
    double precision, dimension(nalf*max2d) :: eigtmp
    double precision, dimension(max2d,max3d) :: vecs2d
    dimension nv2(nalf)
    save itotal
    autocm = 2.19474624d+05

    if (igamma == 1) then
        itotal = 0
        rewind intvec
    endif

    do 10 i=1,nham2
        write(intvec) eig2(i)
        write(intvec) (ham2(j,i),j=1,nham2)
10    continue

    iprev=itotal
    inew=nham2
    do 20 ipos = itotal + nham2,1,-1
        if(iprev /= 0) then
          if(inew /= 0) then
            if(eig2(inew) < eigtmp(iprev)) then
              eigtmp(ipos) = eigtmp(iprev)
              iprev = iprev - 1
            else
              eigtmp(ipos) = eig2(inew)
              inew = inew - 1
            endif
          else
            goto 30
          endif
        else
          eigtmp(ipos) = eig2(inew)
          inew = inew - 1
        endif
20    continue
30    continue

    itotal = itotal + nham2

    if (igamma == nalf) then
        rewind intvec

        if(itotal <= max3d) then
          emax = eigtmp(itotal)
          write(6,100) itotal,eigtmp(1),emax
        else
          emax = eigtmp(max3d)
          write(6,100) max3d,eigtmp(1),emax
          write(6,110) max3d,eigtmp(1)*autocm,emax*autocm
        endif
        itotal = 1
        do 40 i=1,nalf
          iv2(1,i) = itotal
          ichose=0
          do 50 j=1,nv2(i)
            read(intvec) eig
            if(eig <= emax) then
              eigs2d(itotal) = eig
              read(intvec) (vecs2d(k,itotal),k=1,nv2(i))
              ichose = ichose + 1
              itotal = itotal + 1
            else
              read(intvec)
            endif
50        continue
          iv2(2,i) = ichose
40      continue
    endif

100   format(5x,'selecting ',i5,' energies between ',e14.6,&
            ' and ',e14.6,' in hartrees')
110   format(5x,'selecting ',i5,' energies between ',e14.6,&
            ' and ',e14.6,' in wavenumbers')
    return
    end

!     ***********************************************************************

    subroutine cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)

!     this routine selects all the eigenvalues that are lower than the
!     the cut-off energy emax1, which is user-supplied in wavenumbers.
!     these eigenvalues & their corresponding vectors are then saved
!     in the array vecs1d.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision(a-h,o-y),logical(z)

    double precision, dimension(nham2) :: eig2
    double precision, dimension(nham2,nham2) :: ham2
    dimension iv2(2,nalf)
    double precision, dimension(max3d) :: eigs2d
    double precision, dimension(max2d,max3d) :: vecs2d
    npos
    autocm = 2.19474624d+05
    if (igamma == 1) npos = 1

!     change emax2 to hartree for the selection
    emaxau=emax2/autocm

    nvec = 0
    do 10 ialpha = 1,nham2
        if (eig2(ialpha) < emaxau) then
          ntot = npos + nvec
          if (ntot > max3d) then
             write(6,140) max3d, emax2
             stop
          endif
          eigs2d(ntot) = eig2(ialpha)
          do 20 j = 1,nham2
            vecs2d(j,ntot) = ham2(j,ialpha)
20        continue
          nvec = nvec + 1
        endif
10    continue
    iv2(1,igamma) = npos
    iv2(2,igamma) = nvec
    npos = npos + nvec
    return
140   format('number of 2d eigenvalues greater than ',i4,&
            ' increase max3d or reduce emax2 which is ',e12.5,&
            ' cm-1.')
    end

!##################################################################
    subroutine testiv(iv,nbass)

!     selection vectors for the bisector embedding to ensure that
!     singular region of theta = 0 is not sampled when j > 0.
!     also calculate which angular grid points are redundant.
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision (a-h,o-y), logical (z)

    dimension iv(2,nalf)

    if (jrot > 0) then
!        first find the extent of the functions in low theta direction
        do 20 iend=nalf/2,nalf
        if (iv(2,iend) == 0) goto 30
   20    continue
        iend=nalf
!        have we saved functions beyond the point of zero amplitude ?
   30    ioff=0
        do 40 i=iend+1,nalf
!        if so remove them
        if (iv(2,i) > 0) then
            ioff=ioff+iv(2,i)
            iv(2,i)=0
            iv(1,i)=iv(1,i-1)
        endif
   40    continue

!         zlin = .true.

        if (ioff > 0) write(6,960) ioff,iv(1,iend)
960    format(/5x,'*** warning:',i4,' functions removed from theta =',&
                ' 0 region'/8x,' basis reset to nham3 =',i5)
        if (iv(2,nalf) > 0) then
            if (zlin) then
!!              if there are still theta=0 functions, remove them
               write(6,987) iv(2,nalf),iv(1,nalf)-1
               iv(2,nalf)=0
            else
!!              wavefunction has amplitude all the way to theta = 0: stop
               write(6,950)
               stop
            endif
        endif
950    format(/5x,'bisector embedding: ',&
                    'wavefunction has amplitude for theta = 0. stop.')
987    format(/5x,'*** warning: zlin = t forced the removal of',i4,&
               ' functions at theta = 0'&
                /8x,' basis reset to  nham3 =',i5)
    endif

    iang=0
    do 60 ii=1,nalf
        if (iv(2,ii) > 0) iang=iang+1
60 continue
    nbass=iang*max2d
    if (ztran) then
        write(iwave) iang,nbass
        write(iwave) (iv(2,ii),ii=1,nalf)
    endif
    return
    end

!***********************************************************************

    subroutine bloc3d(htheta,ham3,eigs2d,vecs2d,iv2,nv2,ham2,nham3,r)
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_oupb

    implicit double precision(a-h,o-y),logical(z)


    double precision, dimension(nalf,nalf) ::  htheta
    double precision, dimension(nham3,nham3) ::  ham3
    double precision, dimension(max3d) :: eigs2d
    double precision, dimension(max2d,max3d) :: vecs2d
    dimension iv2(2,nalf),nv2(nalf)
    double precision, dimension(max2d,max2d) :: ham2
    double precision, dimension(npnt) :: r

! zero ham3
    ham3 = 0.0d0

    ndim1g = 1
    do 10 igamma = 1,nalf
        if(iv2(2,igamma) > 0) then
          ndim2g = 1
          do 20 igammp = 1,igamma
            if(iv2(2,igammp) > 0) then
              call blc2d2(r,igamma,igammp,htheta,ham2,nv2(igamma))

              call vecmul(vecs2d(1,iv2(1,igamma)),nv2(igamma),&
                          iv2(2,igamma),vecs2d(1,iv2(1,igammp)),&
                          nv2(igammp),iv2(2,igammp),max2d,&
                          ham2,ham3(ndim1g,ndim2g),nham3)
              ndim2g = ndim2g + iv2(2,igammp)
            endif
20        continue
          ndim1g = ndim1g + iv2(2,igamma)
        endif
10    continue

    do 60 i = 1,nham3
        ham3(i,i) = ham3(i,i) + eigs2d(i)
60    continue
    if (zpham) call wrtham(ham3,nham3)

    return
    end

!    ***********************************************************************

    subroutine blc2d2(r,igamma,igammp,htheta,ham2,nham2)
    use dvr3drjz_size
    use dvr3drjz_split1
    use dvr3drjz_oupb

    implicit double precision(a-h,o-y),logical(z)


    double precision, dimension(nalf,nalf) ::  htheta
    double precision, dimension(npnt) :: r
    double precision, dimension(nham2,nham2) :: ham2
    xp5 = 0.5d0
    factr2 = sqrt(xp5)

    ham2 = 0.0d0

!      q=0.0d0
    q=1.0d0
    if(ipar == 1) q=-1.0d0
    iap=0
    do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
              if(ialpha == ialphp .and. ibeta == ibetap) then
                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
                wsum = walpha + wbeta

! Schwenke angular term in the program (it is added to the normal term)
! Kolya
                wsum=wsum-xp1*1.5468847900d-8*1.250d0

                ham2(ialphp+iap,ialpha+ia)=& 
                ham2(ialphp+iap,ialpha+ia)+htheta(igammp,igamma)*wsum
              endif
              if(ialpha == ibetap .and. ibeta == ialphp) then
                walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
                wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
                wsum = walpha + wbeta
                ham2(ialphp+iap,ialpha+ia)=& 
                ham2(ialphp+iap,ialpha+ia)+q*htheta(igammp,igamma)*wsum
              endif
              if(ialphp == ibetap) ham2(ialphp+iap,ialpha+ia)=& 
                ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha == ibeta) ham2(ialphp+iap,ialpha+ia)=& 
                ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
    return
    end

!     ***********************************************************************

    subroutine vecmul(veca,idima,jdima,vecb,idimb,jdimb,nvecln,hama,& 
                        hamb,ndim)

! this routine does hb = hb + va^T * ha * vb
! could be replaced by blas?

    implicit double precision(a-h,o-y),logical(z)
    double precision, dimension(nvecln,jdima) ::  veca
    double precision, dimension(nvecln,jdimb) :: vecb
    double precision, dimension(idima,idimb) :: hama
    double precision, dimension(ndim,jdimb) :: hamb

    do 10 ib = 1,idimb
        do 20 ia = 1,idima
          temp1 = hama(ia,ib)
          if (temp1 == 0.0d0) goto 20
          do 30 ja = 1,jdima
            temp2 = veca(ia,ja) * temp1
            do 40 jb = 1,jdimb
              hamb(ja,jb) = hamb(ja,jb) +  vecb(ib,jb) * temp2
40          continue
30        continue
20      continue
10    continue
    return
    end

!     ***********************************************************************

    subroutine transr(iv2,vecs2d,ham3,eig3,nham3,nbass)
    use dvr3drjz_size
    use dvr3drjz_size      
    use dvr3drjz_outp
    use dvr3drjz_split1
    use dvr3drjz_oupb
    use dvr3drjz_mass
    implicit double precision(a-h,o-y),logical(z)      


    dimension iv2(2,nalf)
    double precision, dimension(max2d,max3d) :: vecs2d
    double precision, dimension(nham3,nham3) :: ham3
    double precision, dimension(nham3) :: eig3
    double precision, dimension(nbass) :: wvfunc

    call nftim('start of wavefunction generation')
    write(iwave) neval
    call outrow(eig3,neval,iwave)

    do 40 level3=1,neval
        index=1
        wvfunc = 0.0d0
        index = 1    
    do 30 irs=1,(npnt*(npnt+1-2*ipar))/2
          do 300 igamma=1,nalf
            if (iv2(2,igamma) == 0) goto 300
            do 20 level2=0,iv2(2,igamma)-1
              ivec2=iv2(1,igamma)+level2
              wvfunc(index)=wvfunc(index)+&
              vecs2d(irs,ivec2)*ham3(ivec2,level3)
20          continue
            index=index+1
300       continue
30      continue
        write(iwave) wvfunc
40    continue

    call nftim('end of wavefunction generation')
    return
    end

!     ***********************************************************************

    subroutine trans(iv1l,iv2l,ndim2l,vecs1l,&
                       vecs2l,vecs3l,phi,psi,evall,nham3)

!     if ztran then this routine transforms the sets of 1d, 2d and 3d
!     coefficients to psi, the wavefunction amplitudes at the dvr points
    use dvr3drjz_size
    use dvr3drjz_outp
    use dvr3drjz_split1
    use dvr3drjz_split2
    use dvr3drjz_oupb
    use dvr3drjz_mass
    implicit double precision (a-h,o-y), logical (z)

    double precision, dimension(max2d,ndima) :: vecs1l
    double precision, dimension(max2d,max2d) :: vecs2l
    double precision, dimension(nham3) :: vecs3l
    dimension iv1l(ndimc,ndimb)
    dimension ndim2l(ndimc)
    double precision, dimension(idvr,npnt1,npnt2) :: psi
    dimension iv2l(ndimc)
    double precision, dimension(neval) :: evall
    double precision, dimension(nham3,ndima,ndimb) :: phi

    rewind iout1
    rewind iout2
    rewind ivecs1
    rewind ivecs2

!     skip header on iout2

    do 10 i=1,5
        read(iout2)
   10 continue

!     read in iv1 and iv2 from a seperate stream (iout1)
    call igetro(iv1l,ndimb*ndimc,iout1)
    call igetro(iv2l,ndimc,iout1)

    do 32 ione = 1,npntc
    nham2 = 0
!     recall the size of ham2 for each npntc
      do 3 itwo = 1,npntb
        nham2 = nham2 + iv1l(ione,itwo)
    3 continue
    ndim2l(ione) = nham2
   32 continue

!     now rewrite the vectors to different streams

    do 35 ione = 1,npntc
    nham2 = ndim2l(ione)
    if (nham2 > 0) then
        do 5 kk=1,npnta
        call getrow(vecs1l(1,kk),nham2,iout2)
        if (iv2l(ione) /= 0) call outrow(vecs1l(1,kk),nham2,ivecs1)
    5    continue
    endif
   35 continue
    rewind ivecs1

    do 31 ione = 1,npntc
    nham2 = ndim2l(ione)
    ivm = iv2l(ione)
    if (nham2 > 0) then
        do 30 ind = 1,ivm
        call getrow(vecs2l(1,ind),nham2,iout2)
        call outrow(vecs2l(1,ind),nham2,ivecs2)
   30    continue
    endif
   31 continue
    rewind ivecs2

    read (iout2) meval
    write(iwave) meval

!     call the eigenvalues into evall
    call getrow(evall,meval,iout2)
    call outrow(evall,meval,iwave)

!     now ready to start reading the 3-d vectors

!     1-d vectors:  ivecs1 (stream 26)
!     2-d vectors:  ivecs2 (stream 27)
!     3-d vectors:  iout2  (stream 25)

    ind2 = 0
    do 58 ic = 1,npntc
    nham2 = ndim2l(ic)
    do 56 j  = 1,iv2l(ic)
    ind2 = ind2 + 1
    if (nham2 > 0) call getrow(vecs2l(1,j),nham2,ivecs2)
    do 54 ia = 1,npnta
    if (nham2>0 .and. j==1)&
    call getrow(vecs1l(1,ia),nham2,ivecs1)
    ind1 = 0
    do 52 ib = 1,npntb
    sum1 = 0.0d0
    do 50 k  = 1,iv1l(ic,ib)
    ind1 = ind1 + 1
    sum1 = sum1 + vecs1l(ind1,ia)*vecs2l(ind1,j)
50 continue
    phi(ind2,ia,ib) = sum1
52 continue
54 continue
56 continue
58 continue

    do 959 ll=1,meval

    call getrow(vecs3l,nham3,iout2)

!     set psi to zero
    psi = 0.0d0

    if (ztheta) then

        if (.not. zr2r1) then

          ind2 = 0
          do 158 ic = 1,npnt2
          do 156 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 154 ia = 1,idvr
          do 152 ib = 1,npnt1
          psi(ia,ib,ic) = psi(ia,ib,ic) + vecs3l(ind2)*phi(ind2,ia,ib)
152     continue
154     continue
156     continue
158     continue

        else

          ind2 = 0
          do 258 ic = 1,npnt1
          do 256 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 254 ia = 1,idvr
          do 252 ib = 1,npnt2
          psi(ia,ic,ib) = psi(ia,ic,ib) + vecs3l(ind2)*phi(ind2,ia,ib)
252     continue
254     continue
256     continue
258     continue

        endif

    else

        if (.not. zr2r1) then

          ind2 = 0
          do 358 ic = 1,idvr
          do 356 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 354 ia = 1,npnt1
          do 352 ib = 1,npnt2
          psi(ic,ia,ib) = psi(ic,ia,ib) + vecs3l(ind2)*phi(ind2,ia,ib)
352     continue
354     continue
356     continue
358     continue

        else

          ind2 = 0
          do 458 ic = 1,idvr
          do 456 j  = 1,iv2l(ic)
          ind2 = ind2 + 1
          do 454 ia = 1,npnt2
          do 452 ib = 1,npnt1
          psi(ic,ib,ia) = psi(ic,ib,ia) + vecs3l(ind2)*phi(ind2,ia,ib)
452     continue
454     continue
456     continue
458     continue

        endif

    endif

!     now save the wavefunction for each level (to stream 28)
    call outrow(psi,idvr*npnt1*npnt2,iwave)

959 continue
    return
    end

!#############################################################################

    subroutine sqout(sqmat,ndim)
!     print lower triangle of square matrix
    double precision sqmat(ndim,ndim)
    do 30 i=1,ndim
    write(6,"(10f13.7)") (sqmat(i,j),j=1,i)
30 continue
    return
    end

!###########################################################################
    subroutine symout(symmat,ndim)

!     print out lower triangle of symmetric matrices                #008

    double precision symmat(1)
    ip=0
    3 llow=10*ip+1
    lup=min(llow+9,ndim)
    write(6,"(/,i11,9i13)") (i,i=llow,lup)
    ind0=llow*(llow+1)/2
    do 5 i=llow,ndim
    itop=min(lup,i)
    ktop=itop-llow+ind0
    write(6,"(i4,f12.7,9f13.7)") i,(symmat(j),j=ind0,ktop)
    ind0=ind0+i
    5 continue
    if(lup>=ndim) return
    ip=ip+1
    go to 3
    end

    subroutine wrtham(hamil,nham)
!     print hamiltonian matrix                                      #011
    double precision hamil(nham,nham)
    write(6,"(5x,'hamiltonian matrix'/)")
    do 30 i=1,nham
    write(6,"(10f13.7)") (hamil(i,j),j=1,i)
   30 continue
    return
    end

    subroutine getrow(row,nrow,iunit)
!     fast non-formatted read
    double precision row(nrow)
    read(iunit) row
    return
    end

    subroutine igetro(ivec,nsize,istream)
    dimension ivec(nsize)
    read(istream) ivec
    return
    end

    subroutine outrow(row,nrow,iunit)
!     fast non-formatted write                                      #025
    double precision row(nrow)
    write(iunit) row
    return
    end

    subroutine ioutro(ivec,nsize,istream)
    dimension ivec(nsize)
    write(istream) ivec
    return
    end
    subroutine nftim(text)
    use dvr3drjz_timing
    character text*(*)
    write(6,"(/)")
    write(6,*) 'Time at ',text,' is.........'
    call SYSTEM_CLOCK(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,' secs CPU time used'/) ")itime
    return
    end

!########################################################################
    subroutine mkmain(hr,htheta,r,theta,kz)

!     this routine controls the dvr calculation in the case of
!     symmetrised radau coordinates with z axes perpendicular to the 
!     molecular plane.
!     written by max kostin, 2001.
    use dvr3drjz_size
    use dvr3drjz_outp
    implicit double precision(a-h,o-y),logical(z)

    double precision, dimension(npnt,npnt) :: hr
    double precision, dimension(nalf,nalf) :: htheta
    double precision, dimension(npnt) :: r
    double precision, dimension(nalf) :: theta
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::ham2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eig2,eigtmp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::vecs2d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eigs2d
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::ham3
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::eig3
    INTEGER, ALLOCATABLE, DIMENSION(:,:) ::iv2
    INTEGER, ALLOCATABLE, DIMENSION(:) ::nv2

    ALLOCATE(ham2(max2d,max2d),eig2(max2d),iv2(2,nalf),&
               vecs2d(max2d,max3d),eigs2d(max3d),nv2(nalf))
    if (.not.zcut) ALLOCATE(eigtmp(nalf*max2d))

! no need for the 1d diagonalisations as there is no possiblity of
! truncation as the symmetry would be broken.
    write(6,"(5x,'starting the 2d calculations.')")
    call nftim('beginning of 2d loop')

    do 30 igamma = 1,nalf     
        nv2(igamma) = (npnt*(npnt+1-ipar*2))/2      
        nham2 = nv2(igamma)
        call z_blc2d1(theta(igamma),r,hr,ham2,nham2,&
                          term,term2,kz)

        call diag(ham2,nham2,nham2,eig2)

        if (.not. zcut) then
          call choosr(igamma,nham2,eig2,ham2,iv2,eigs2d,&
                    vecs2d,nv2,eigtmp)

          if(igamma == nalf .and. .not. zpmin)&
            write(6,"(5x,'for gamma = ',i2,' selected ',i3,' energies.')") &
            &(itmp, iv2(2,itmp),itmp=1,nalf)
        else
           call cut2dr(igamma,nham2,eig2,ham2,iv2,eigs2d,vecs2d)
          if (.not. zpmin) write(6,"(5x,'for gamma = ',i2,' selected ',i3,' energies.')")&
            & igamma, iv2(2,igamma)
        endif
30    continue

    call testiv(iv2,nbass)

    nham3 = iv2(1,nalf) + iv2(2,nalf) - 1
    if (.not. zpmin) write(6,"(25x,'total = ',i5)") nham3
    call nftim('end of 2d loop')
    write(6,"(5x,'building the 3d hamiltonian.')")

    if (.not.zcut) DEALLOCATE(eigtmp)
    ALLOCATE(ham3(nham3,nham3))

    call bloc3d(htheta,ham3,eigs2d,vecs2d,&
                iv2,nv2,ham2,nham3,r)

    DEALLOCATE(ham2,eig2,eigs2d,nv2)
 
    call nftim('end of 3d ham building')
    write(6,"(5x,'diagonalising the 3d hamiltonian.') ")
    ALLOCATE(eig3(nham3))

    call diag3d(ham3,nham3,eig3,kz)

    call nftim('end of diagonalising 3d')

    if (ztran) call transr(iv2,vecs2d,ham3,eig3,nham3,nbass)

    DEALLOCATE(iv2,vecs2d,ham3,eig3)
    return
    
    end
      
!########################################################################     


    subroutine z_blc2d1(xcos,r,hr,ham2,nham2,term,term2,kz)
    use dvr3drjz_size
    use dvr3drjz_split1
    implicit double precision(a-h,o-y),logical(z)

    double precision, dimension(npnt,npnt) :: hr
    double precision, dimension(npnt) :: r
    double precision, dimension(nham2,nham2) :: ham2
    data xp5/0.5d0/,x1/1.0d0/,x4/4.0d0/
    factr2 = sqrt(xp5)

    ham2 = 0.0d0
    realj = DBLE(jrot)
    realkz = DBLE(kz)

!     allow for j > 0 case
    if (jrot /= 0) then
! for q=0(+) term goes with -
! for q=1(-) term goes with +
        if(iq == 0)then
            term= -realj*(realj+x1)/x4
        else
            term= realj*(realj+x1)/x4
        endif
        fact = (realkz**2)/x4
        if (kz == 1) fact = fact + term * (xcos)/(x1-(xcos)**2)

!-----  Extra NBO term ------
        s=1.0d0-urr1/ur1
        term2=(realj**2+realj-realkz**2)/4.0d0
        fact2=s*term2*x1/(x1-(xcos)**2)
!----------------------------
         
        ia = 0
        do 15 ibeta=1,npnt
            do 25 ialpha=1,ibeta-ipar
!           walpha = xp5 / (r(ialpha)*r(ialpha)*ur1)
!           wbeta = xp5 / (r(ibeta)*r(ibeta)*ur1)
               walpha = xp5 / (r(ialpha)*r(ialpha)*urr1)
               wbeta = xp5 / (r(ibeta)*r(ibeta)*urr1)
               wsum = (walpha + wbeta) * fact
               wsum2= (walpha + wbeta) * fact2
               ham2(ialpha+ia,ialpha+ia) = ham2(ialpha+ia,ialpha+ia)+wsum+wsum2
25          continue
            ia=ia+ibeta-ipar
            if (ipar == 0)ham2(ia,ia) = ham2(ia,ia) + wsum + wsum2
 15     continue
    endif

    q=x1
    if(ipar == 1) q=-x1
    iap=0
    do 10 ibetap=1,npnt
        ia=0
        do 20 ibeta=1,npnt
          do 30 ialphp=1,ibetap-ipar
            do 40 ialpha=1,ibeta-ipar
              if(ibeta == ibetap)&
                 ham2(ialphp+iap,ialpha+ia)=&
                 ham2(ialphp+iap,ialpha+ia)+hr(ialphp,ialpha)
              if(ibeta == ialphp)&
                 ham2(ialphp+iap,ialpha+ia)=&
                 ham2(ialphp+iap,ialpha+ia)+q*hr(ibetap,ialpha)
              if(ialpha == ibetap)&
                 ham2(ialphp+iap,ialpha+ia)=&
                 ham2(ialphp+iap,ialpha+ia)+q*hr(ialphp,ibeta)
              if(ialpha == ialphp)&
                 ham2(ialphp+iap,ialpha+ia)=&
                 ham2(ialphp+iap,ialpha+ia)+hr(ibetap,ibeta)
              if(ialpha == ialphp .and. ibeta == ibetap) then
                call potv(v,r(ialpha),r(ibeta),xcos)
                ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)+v
              endif
              if(ialpha == ibetap .and. ibeta == ialphp) then
                call potv(v,r(ibeta),r(ialpha),xcos)
                ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)+q*v
              endif
              if(ialphp == ibetap) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
              if(ialpha == ibeta) ham2(ialphp+iap,ialpha+ia)=&
                ham2(ialphp+iap,ialpha+ia)*factr2
40          continue
30        continue
          ia=ia+ibeta-ipar
20      continue
        iap=iap+ibetap-ipar
10    continue
    return
    end

!#################################################################

    MODULE constants
    IMPLICIT NONE
    INTEGER, PARAMETER :: real_kind=SELECTED_REAL_KIND(8,40)
    END MODULE constants

    SUBROUTINE jac_basis(nn,nb,alf,bet,x,basis)
    USE constants
    IMPLICIT NONE
    INTEGER :: n,i,nb,nn
    REAL(KIND=real_kind) :: alf, bet
    REAL(KIND=real_kind) :: x(nn),basis(0:nb,nn)
    REAL(KIND=real_kind) :: bass(0:nb,nn),norm(0:nb)

    CALL norms2(norm,nb,alf,bet)
    CALL jacgt(x,bass,alf,bet,nn,nb)

    DO i=0,nb
        DO n=1,nn
            basis(i,n)=bass(i,n)*norm(i)
        END DO
    END DO

    RETURN
    END SUBROUTINE jac_basis

!##################################################################

    SUBROUTINE jacgt(x,bass,alf,bet,nn,nb)
    USE constants
    IMPLICIT NONE
    INTEGER :: n,I,nn,nb
    REAL(KIND=real_kind),EXTERNAL :: gammln
    REAL(KIND=real_kind) :: alf, bet,lmd, x0,x1,x2
    REAL(KIND=real_kind) :: x(nn),bass(0:nb,nn)
    REAL(KIND=real_kind) :: A1n(nb),A2n(nb),A3n(nb),A4n(nb)
    DATA x0,x1,x2/0.0d0,1.0d0,2.0d0/
    lmd=alf+bet+x1
    DO n=1,nb
        A1n(n)=x2*(n+x1)*(n+lmd)*(x2*n+lmd-x1)
        A2n(n)=(x2*n+lmd)*(alf*alf-bet*bet)
        A3n(n)=(x2*n+lmd-x1)*(x2*n+lmd)*(x2*n+lmd+x1)
        A4n(n)=x2*(n+alf)*(n+bet)*(x2*n+lmd+x1)
    END DO   
    bass=x0
    DO 60 I=1,nn
    bass(0,I)=x1
    IF(nb<1) GO TO 70
    bass(1,I)=(alf-bet+(lmd+x1)*x(I))/x2
    DO 80 n=2,nb
        bass(n,I)=((A2n(n-1)+A3n(n-1)*x(I))*bass(n-1,I)&
        &         -A4n(n-1)*bass(n-2,I))/A1n(n-1)
 80   CONTINUE
 70   CONTINUE 
 60   CONTINUE
    RETURN
    END SUBROUTINE jacgt

!##################################################################

    SUBROUTINE gaujac(x,w,n,alf,bet)
    USE constants
    IMPLICIT NONE
    INTEGER :: n
    REAL(KIND=real_kind):: alf,bet,x1,x2,x3
    REAL(KIND=real_kind):: w(n),x(n)
    REAL(KIND=real_kind),PARAMETER :: EPS=3.0D-14
    INTEGER,PARAMETER :: MAXIT=10
    INTEGER :: i,its,j
    REAL(KIND=real_kind)::alfbet,an,bn,r1,r2,r3
    REAL(KIND=real_kind)::c1,c2,c3,p1,p2,p3,pp,temp,z,z1
    REAL(KIND=real_kind),EXTERNAL :: gammln
    DATA x1,x2,x3/1.0d0,2.0d0,3.0d0/
    DO 13 i=1,n
        IF(i==1)THEN
            an=alf/DBLE(n)
            bn=bet/DBLE(n)
            r1=(x1+alf)*(2.78D0/(4.0D0+DBLE(n*n))+0.768D0*an/DBLE(n))
            r2=x1+1.48D0*an+0.96D0*bn+0.452D0*an*an+0.83D0*an*bn
            z =x1-r1/r2
        ELSE IF(i==2)THEN
            r1=(4.1D0+alf)/((x1+alf)*(x1+0.156D0*alf))
            r2=x1+0.06D0*(DBLE(n)-8.0D0)*(1.0D0+0.12D0*alf)/DBLE(n)
            r3=x1+0.012*bet*(x1+0.25D0*DABS(alf))/DBLE(n)
            z=z-(x1-z)*r1*r2*r3
        ELSE IF(i==3)THEN
            r1=(1.67D0+0.28D0*alf)/(x1+0.37D0*alf)
            r2=x1+0.22D0*(DBLE(n)-8.0D0)/DBLE(n)
            r3=x1+8.0D0*bet/((6.28D0+bet)*DBLE(n*n))
            z=z-(x(1)-z)*r1*r2*r3
        ELSE IF(i==n-1)THEN
            r1=(x1+0.235D0*bet)/(0.766D0+0.119D0*bet)
            r2=x1/(x1+0.639D0*(DBLE(n)-4.0D0)/&
            &       (x1+0.71D0*(DBLE(n)-4.0D0)))
            r3=x1/(x1+20.0D0*alf/((7.5D0+alf)*DBLE(n*n)))
            z=z+(z-x(n-3))*r1*r2*r3
        ELSE IF(i==n)THEN
            r1=(x1+0.37D0*bet)/(1.67D0+0.28D0*bet)
            r2=x1/(x1+0.22D0*DBLE(n-8)/DBLE(n))
            r3=x1/(x1+8.0D0*alf/((6.28D0+alf)*DBLE(n*n)))
            z=z+(z-x(n-2))*r1*r2*r3
        ELSE
            z=x3*x(i-1)-x3*x(i-2)+x(i-3)
        ENDIF
        alfbet=alf+bet

        DO 12 its=1,MAXIT
            temp=x2+alfbet
            p1=(alf-bet+temp*z)/x2
            p2=x1
            DO 11 j=2,n
               p3=p2
               p2=p1
               temp=x2*DBLE(j)+alfbet
               c1=x2*DBLE(j)*(DBLE(j)+alfbet)*(temp-x2)
               c2=(temp-x1)*(alf*alf-bet*bet+temp* &
              &               (temp-x2)*z)
               c3=x2*(DBLE(j-1)+alf)*(DBLE(j-1)+bet)*temp
               p1=(c2*p2-c3*p3)/c1
 11         CONTINUE 
             pp=(DBLE(n)*(alf-bet-temp*z)*p1+x2*(DBLE(n)+alf)* &
           &    (DBLE(n)+bet)*p2)/(temp*(x1-z*z))
            z1=z
            z=z1-p1/pp
            IF(ABS(z-z1)<=EPS) GOTO 1
 12         CONTINUE            
     1   x(i)=z 
        w(i)=DEXP(gammln(alf+DBLE(n))+gammln(bet+DBLE(n))    &
        &    -gammln(DBLE(n+1))-gammln(DBLE(n)+alfbet+x1))&
        &      *temp*x2**alfbet/(pp*p2)

 13    CONTINUE
    RETURN
    END SUBROUTINE gaujac

    FUNCTION gammln(XX)
    USE constants
    IMPLICIT NONE
    INTEGER :: j
    REAL(KIND=real_kind)::GAMMLN,XX
    REAL(KIND=real_kind)::SER,STP,TMP,X,COF(6)
    REAL(KIND=real_kind)::HALF,ONE,FPF
    DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
      &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
    DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
    X=XX-ONE
    TMP=X+FPF
    TMP=(X+HALF)*LOG(TMP)-TMP
    SER=ONE
    DO 11 J=1,6
    X=X+ONE
    SER=SER+COF(J)/X
!        SER=SER+COF(J)/(X+0.000001d0)
11  CONTINUE
    GAMMLN=TMP+LOG(STP*SER)
    RETURN
    END FUNCTION gammln

!####################################################################

    SUBROUTINE norms2(norm,nn,alf,bet)
    USE constants
    IMPLICIT NONE
    INTEGER :: n,nn
    REAL(KIND=real_kind) :: alf, bet,lmd,x1,x2,norm1
    REAL(KIND=real_kind) :: norm(0:nn)
    REAL(KIND=real_kind) :: a1,a2,a3,a4
    REAL(KIND=real_kind), EXTERNAL :: gammln
    DATA x1,x2/1.0d0,2.0d0/
    lmd=alf+bet+x1
    do n=0,nn
        a1=gammln(DBLE(n+1))
        a2=gammln(DBLE(n)+lmd)
        a3=gammln(DBLE(n)+alf+x1)
        a4=gammln(DBLE(n)+bet+x1)
        norm1=2**(-lmd)*(x2*DBLE(n)+lmd)*exp(a1+a2-a3-a4)
        norm(n)=SQRT(norm1)
    end do
    END SUBROUTINE norms2