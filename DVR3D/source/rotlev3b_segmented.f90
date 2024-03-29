!==================================================================================================
!Module defintion
!By Runyu Zhang & Tennyson Jonathan 1/Sep/2022
!Contains: size, outp, timing
!Special notice:: The module name contains the file name, it means cannot directly paste to other 
!                 files and use. There are some difference between the constant numbers, types and 
!                 default value.
!==================================================================================================
module rotlev3b_seg_size
    integer :: nbass   ! maximum dimension of rotational secular problem
    integer :: mbass     ! maximum size of vibrational problem (excluding linear geom)
    integer :: IBASS    ! actual dimension of rotational secular problem
    integer :: NEVAL    ! number of eigenvalues which have to actually be supplied as output
    integer :: IPAR     ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
    integer :: nmax     ! number of dvr points in each radial coordinate
    integer :: maxblk   ! size of vibrational radial problem (even basis)
    integer :: jrot     ! total angular momentum of the molecule
    integer :: kmin     ! zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
                         ! kmin=2 loop over both sym & anti-sym (zbisc=t only)
                         ! zrot=f, kmin=fixed value of k
    integer :: kmax
    integer :: meval    ! number of eigenvalues computed in the vibrational problem
    integer :: ndvr  ! maximum dimension of theta dvr grid used in vibrational problem
    integer :: iang  ! maximum number of discrete angles retained in vib. problem
    integer :: npnt   ! max(npnt1,npnt2) number of gauss-associated legendre grid points requested
    integer :: keval    ! number of eigenvectors used for iterations (=neval+4)
    integer :: nvib  ! number of vibrational eigenvalues used in rotational prob.
    integer :: mxblk2  ! size of vibrational radial problem (odd  basis)
    integer :: neval2   ! neval for f block when kmin>1.
    integer :: nblk    ! number of k values
    integer :: LOFF     ! space required for all the off-diagonal blocks
    integer :: LOFF0    ! space required for the largest off-diagonal block
    integer :: mbass0  ! maximum size of vibrational problem (excluding linear geom)
end module rotlev3b_seg_size

module rotlev3b_seg_outp

    integer :: ilev = 14     ! stream for final eigenvalues (formatted).
                                   ! holds input/output of eigenvalues used if zpfun = .true.
    integer :: ivec = 26     ! holds input  eigenvalues & eigenvectors used if always
    integer :: ivec2 = 4   ! holds input  eigenvalues & eigenvectors used if nblk > 2
    integer :: jvec = 3    ! holds output first  set eigenvalues & vectors used if zvec=.true.
    integer :: jvec2 = 2    ! holds output second set                        zvec=.true.
    integer :: kvec = 8    ! holds output first  set transformed vectors used if ztran=.true.
    integer :: kvec2 = 9    ! holds output second set   used if ztran=.true.
    integer :: iscr = 1    ! scratch file used for restart runs
                                   ! holds hamiltonian file used if always

    integer :: ires = 0     ! = 0  normal run
                                   ! = 1  restart from first  call to dgrot
                                   ! = 2  restart from second call to dgrot
                                   ! = 3  restart from first  call to dgrot, one diagonalisation only
                                   ! = -1 perform both transformations
                                   ! = -2 perform second transformation only
                                   ! = -3 perform first  transformation only
                                   ! (restart after zdiag=.false. run, ivec=irf1 and irf2 required)
    integer :: irf1 = 21     ! irf1   restart file one  used if zdiag=.false.
    integer :: irf2  = 22     ! irf2   restart file two  used if always
    logical :: zpham = .false.     ! T request printing of the hamiltonian matrix
    logical :: zpvec = .false.    ! T request printing of the eigenvectors
    logical :: zvec = .false.     ! T store the eigenvectors from all the parts of the calculation
                                   ! eigenvalues and eigenvectors written to disk file.
                                   ! (1d,2d and 3d) on stream iout2.
                                   ! further information relating to this (arrays iv1 and iv2) is
                                   ! stored on stream iout1.
    logical :: ztran = .false.     ! T perform the transformation of the solution coefficients
                                   ! to the expression for the wavefunction amplitudes at the grid
                                   ! points. store the data on stream iwave.
                                   ! ztran = T automatically sets zvec = t for idia > -2.
    logical :: zptra = .false.     ! print the transformed vectors.
    logical :: zpfun = .false.    ! F store energy levels on stream ilev
    logical :: zdcore = .false.    ! T for in core diagonalisation
    logical :: zdiag = .false.    ! F do not do final diagonalisation, instead the final Hamiltonian
                                   ! matrix is written on units IDIAG1 and IDIAG2. 
    double precision :: toler = 0.0d0    ! convergence tolerance for the iterative diagonaliser
                                        ! toler = 0.0 gives machine accuracy

    double precision :: thresh = 0.1d0    ! threshold for printing a coefficient if zpvec=.true.

end module rotlev3b_seg_outp

module rotlev3b_seg_timing
    integer :: itime0
end module rotlev3b_seg_timing

!     dummy main program                                           #001
    call rotlev3b
    stop
    end

!#######################################################################
    subroutine rotlev3b

!     program               R O T L E V 3 B
!
!     program to do rotational analysis using vibrational output of
!     dvr3drjz run with idia = -2 and zrot = .true.
!
!     see:
!     j.tennyson & b.t.sutcliffe, mol.phys. 58, 1067 (1986).
!     b.t.sutcliffe, j.tennyson & s.miller, computer. phys. comms.,
!     51, 73 (1988).
!     j.tennyson & b.t.sutcliffe, int.j.quantum chem., 42, 941 (1992).
!
!     use as follows:
!     comments on namelist parameters (& defaults) in block data
!     the program needs the following subroutines:
!     1. limited card input which is read in subroutine insize
!        and an input file on stream ivec (& ivec2) from triatom
!     2. f02fjf to do iterative diagonalisation (nag routine).
!     or dsyev to do in core diagonalisation (lapack f90 routine).
!     the program works in **** atomic units ***** :
!     Fortan90 version with dynamic arrays by Max Kostin & Jonathan Tennyson
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    use rotlev3b_seg_timing
    namelist/prt/ toler,thresh,zpham,zpvec,zvec,ztran,zptra,&
                    zpfun,ilev,ivec,ivec2,jvec,jvec2,kvec,kvec2,&
                    zdiag,zdcore,iscr,ires,irf1,irf2
    implicit double precision (a-h,o-y), logical (z)

    INTEGER :: count_0, count_rate, count_max, walltime , tstart, tend
      
    write(6,"('1',4x,'Program ROTLEV3B (version of March 2002):'/)")
!     read in namelist input data (defaults in block data)
    read(5,prt)

    call nftim('beginning')
      
    call SYSTEM_CLOCK(itime0,irate2,imax2)
!     read in control parameters of problem.
    call insize

!     first for select, then onto main program
    call select 

    call SYSTEM_CLOCK(itime0,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,' secs CPU time used'/)")itime
    stop
    end

!#######################################################################
 
!########################################################################
    subroutine insize

!     set up common /size/ & write control parameters of problem    #004
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)

!     common /size/ stores control parameters for the problem
!     nbass: maximum dimension of rotational secular problem
!     ibass: actual dimension of rotational secular problem
!     mbass: maximum size of vibrational problem (excluding linear geom)
!     mbass0: maximum size of vibrational problem (including linear geom)
!     maxblk: size of vibrational radial problem (even basis)
!     mxblk2: size of vibrational radial problem (odd  basis)
!     ndvr : maximum dimension of theta dvr grid used in vibrational problem
!     iang : maximum number of discrete angles retained in vib. problem
!     npnt : number of gauss-associated legendre grid points requested
!     jrot : total rotational angular momentum
!     ipar : parity of vibrational basis with lowest k:
!            ipar=0 for even & =1 for odd
!     kmin : kmin=1 for sym. rotational basis, =0 for anti-sym.
!            kmin>1 loop over both.
!     nmax : number of dvr points in each radial coordinate
!     meval: number of eigenvalues computed in the vibrational problem
!     nvib : number of vibrational eigenvalues used in rotational prob.
!     neval: number of eigenvalues which have to actually be computed
!     neval2: neval for f block when kmin>1.
!     keval: number of eigenvectors used for iterations (=neval+4)
!     nblk : number of k values
!     loff : space required for all the off-diagonal blocks
!     loff0: space required for the largest off-diagonal block

    character(len=8) title(9)
    data x0/0.0d0/

    open(unit=ivec,form='unformatted',recordtype='segmented')
    open(unit=irf2,form='unformatted')

    read(ivec) idia,ipar,ndvr,nmax,nmax,jrot,kmin0,meval,nlim
    if (.not. zdiag) then
        open(unit=irf1,form='unformatted')
        write(irf1) idia,ipar,ndvr,nmax,nmax,jrot,kmin0,meval,nlim
    endif

    read(5,"(12i5)")  nvib,neval,kmin,ibass,neval2,npnt
 
    maxblk=nmax*(nmax+1)/2
    mxblk2=nmax*(nmax-1)/2
    mbass0=ndvr*maxblk

    if (kmin /= kmin0 .and. kmin0 /= 2) goto 960
    if (ires < 0) ztran=.true.
    if (ires > 0) zdiag=.true.

!     compute size of rotational secular problem

    nvib=min(nvib,meval)
    if (kmin > 0) then
        nblk=jrot+1
    else
        nblk=jrot
        if (kmin0 == 2) ipar=1-ipar
    endif
    nbass=nblk*nvib
    if (ibass > 0) nbass=min(nbass,ibass)
    if (neval <= 0) neval = 10
    neval=min(neval,nbass)
    if (neval2 <= 0) neval2 = neval
    npnt=max(ndvr,npnt)

    write(6,1000) meval,mbass0,ndvr,nvib,npnt,neval,nbass
1000 format(/5x,'Rotational part of rot-vib calculation  with:',&
             /i9,3x,'lowest vibrational eigenvectors supplied from',&
             /i9,3x,'dimension vibration secular problem with',&
             /i9,3x,'angular dvr points.',&
             /i9,3x,'lowest vibrational eigenvectors actually used',&
             /i9,3x,'point gauss-associated legendre integration',&
             /i9,3x,'lowest rotational eigenvectors required for',&
             /i9,3x,'dimension rotation secular problem')
    if (ibass > 0) write(6,"(12x,'with basis selected by energy ordering')")
    read(5,"(9a8)")   title
    write(6,"(/5x,'Title:',9a8)") title
    if (.not.zdiag) then
        write(6,1012) iscr
1012    format(/5x,'Hamiltonian construction only requested',&
            /5x,'to be written to stream ISCR =',i4)
        zpham=.false.
    else
        if (zdcore) then
            write(6,1013)
1013 format(/5x,'Diagonalisation performed in core using lapack',&
                ' routine dsyev')
        else
            write(6,"(/5x,'Diagonalisation performed iteratively using f02fjf')")
            if (toler /= x0) write(6,"(5x,'Eigenvalue convergence tolerance, TOLER =',d12.3)") toler
            if (toler == x0) write(6,"(5x,'Eigenvalues converged to machine accuracy')")
        endif
    endif
    if (ires /= 0) write(6,"(/5x,'***** restart run, IRES =',i2,' *****')") ires
    if (ires < 0) write(6,"(/'      transformation only')")
    if (zpham) write(6,"(/5x,'Printing of hamiltonian matrix requested')")
    if (.not.zpham) write(6,"(/5x,'Printing of hamiltonian matrix not requested')")
    if (zpvec) write(6,1040) thresh
1040 format(5x,'Printing of eigenvector coefficients greater than',&
                ' thresh =',f5.2,' requested')
    if (.not.zpvec) write(6,"(5x,'Printing of eigenvectors not requested')")
    if (ztran) then
        zvec = .true.
        if (zptra) write(6,"(5x,'Printing of transformed vectors requested')")
        if (.not.zptra) write(6,"(5x,'Printing of transformed vectors not requested')")
    endif
    if (ires == 0) write(6,"(a100,i4)") 'DVR3DRJZ data  to be read         from stream IVEC =',ivec
    if (ires /= 0) write(6,"(a100,i4)") 'DVR3DRJZ data  to be read         from stream IRTF1 =',irf1
    if (nblk > 2) write(6,"(a100,i4)") 'and from stream IVEC2 =',ivec2
    if (zdiag) write(6,"(a100,i4)") 'Hamiltonian scratch file to be held on stream ISCR  =',iscr
    if (zpfun) write(6,"(a100,i4)") 'Eigenvalues    to be written to end of stream ILEV  =',ilev
    if (zvec) write(6,"(a100,i4)") 'Eigenvalues & vectors to be written to stream JVEC  =',jvec
    if (zvec .and. kmin > 1) write(6,"(a100,i4)") 'Second set to be written to stream JVEC2 =',jvec2
    if (ztran) write(6,"(a100,i4)") 'Transformed vectors to be written to stream KVEC  =',kvec
    if (ztran .and. kmin > 1) write(6,"(a100,i4)") kvec2
 1059 format( 5x,'Second set            to be written to stream ',&
             'KVEC2 =',i4)
    if (.not. zdiag .and. ires==0) write(6,"(a100,i4)") 'Restart data file one to be written to stream IRF1  =',irf1
    if (ires == 0) write(6,"(a100,i4)") 'Restart data file two to be written to stream IRF2  =',irf2
    if (ires /= 0) write(6,"(a100,i4)") 'Restart data file two to be read from  stream IRF2  =',irf2

    write(6,"(/5x,'Radau coordinates used')")

    if (ipar == 0) write(6,"('Diatomic assumed homonuclear even parity functions in basis set')")
    if (ipar == 1) write(6,"('Diatomic assumed homonuclear odd parity functions in basis set')")
    write(6,"(/5x,'J =',i3,' rotational state')") jrot
    if (kmin == 0) write(6,"(12x,'with anti-symmetric |Jk> - |J-k> functions in basis')")
    if (kmin == 1) write(6,"(12x,'with symmetric |Jk> + |J-k> functions in basis')")
    if (kmin > 1) write(6,"(12x,'loop over symmetric and anti-symmetric functions')")
    return

960 write(6,"(//5x,'**** kmin =',i3,' but kmin0 =',i3,' STOP ****',/)") kmin,kmin0
    stop
    end

!#######################################################################

    subroutine select
 
!     subroutine select determines which vibrational basis          #007
!     functions are to be used
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)

      
    DOUBLE PRECISION, DIMENSION(NVIB,NBLK) :: EVIB
    DIMENSION IV(NBLK),MVIB(NBLK),nkbas(NBLK),lmin(NBLK),lbasis(NBLK) 
    DOUBLE PRECISION, DIMENSION(3) :: xmass
    DOUBLE PRECISION, DIMENSION(nmax) :: r

    character(len=4) symm(2)
    data symm/'even','odd '/
!     read energies from file ivec,
    if (ires == 0) then
    if (zdiag) then
!        first skip matrix elements
        do 10 i=1,4
        read(ivec,end=900)
10    continue
    else
!       or copy them for a restart run
        read(ivec) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
        write(irf1) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
        read(ivec) re1,diss1,we1,re2,diss2,we2
        write(irf1) re1,diss1,we1,re2,diss2,we2
        read(ivec) dummy
        write(irf1) dummy
        read(ivec) r
        write(irf1) r
    endif
!     read in next set of vibrational vectors & set diagonal elements
    ipt=0
    mbass=0
    iang=0
    do 100 ioff=1,nblk
    read(ivec,end=900) 
    read(ivec,end=900) k2,maxleg,idvr,lincr
    read(ivec,end=900)
!     if k2=0 and we are doing an f parity calculation, read k=1 from
!     end of file
    if (k2 == 0 .and. kmin == 0) then
        call endiv(ivec,jrot+1)
        read(ivec,end=900) k2,maxleg,idvr,lincr
        if (abs(k2) /= 1) then
            write(6,950) ivec,k2
950       format(//' Last block on stream',i3,' has k =',i3,&
                    /' 1 expected: STOP')
            stop
        endif
        read(ivec,end=900)
        read(ivec,end=900) iang2,ibass2
        read(ivec,end=900)
        iang=max(iang,iang2)
        mbass=max(mbass,ibass2)
        read(ivec,end=900) meval2
        mvib(ioff)=min(nvib,ibass2,meval2)
        call getrow(evib(1,ioff),mvib(ioff),ivec)
        call reseti(ivec)
    else
        read(ivec,end=900) iang2,ibass2
        read(ivec,end=900)
        iang=max(iang,iang2)
        mbass=max(mbass,ibass2)
 
        read(ivec,end=900) meval2
        mvib(ioff)=min(nvib,ibass2,meval2)
        call getrow(evib(1,ioff),mvib(ioff),ivec)
        do 30 i=1,meval2
        read(ivec,end=900)
30   continue
    endif
    ipt=ipt+mvib(ioff)
    if (ztran) then
        lmin(ioff)=k2
        nang=maxleg+lincr-k2+1
        nrad=ibass2/iang2
        lbasis(ioff)=nang
        nkbas(ioff)=nrad*nang
    endif
100 continue
!     reposition file ivec
    rewind ivec
    do 110 i=1,3
    read(ivec)
110 continue
    if (ibass <= 0 .or. ibass >= ipt) then
        nbass=ipt
        ivib=nvib
    else
 
!        select the nbass lowest functions
 
        emin=evib(1,1)
        iv=1
        do 160 i=1,nblk
        emin=min(emin,evib(1,i))
160    continue
        ipt=1
        do 200 n=1,ibass
210     if (iv(ipt) <= mvib(ipt)) then
            evibr=evib(iv(ipt),ipt)
            jpt=ipt
        else
            ipt=ipt+1
            goto 210
        endif
        do 220 j=ipt+1,nblk
        if (iv(j) > mvib(j)) goto 220
        if (evib(iv(j),j) >= evibr) goto 220
        evibr=evib(iv(j),j)
        jpt=j
220    continue
!        keep the basis function
        iv(jpt)=iv(jpt)+1
200    continue
!        store the number of functions selected for each k
        mvib(1)=iv(1)-1
        ivib=mvib(1)
        do 230 i=2,nblk
        mvib(i)=iv(i)-1
        ivib=max(ivib,mvib(i))
230    continue
        write(6,1000) nbass,emin,evibr
 1000 format(/i10,' functions selected from e =',d20.10,' to',d20.10,&
             ' hartree')
    endif
!     store re-start info
    write(irf2) nbass,neval,ipar,jrot,kmin,nblk,mbass
    write(irf2) mvib,nkbas,lmin,lbasis

    else
    read(irf2) nbass,neval,ipar,jrot,kmin,nblk,mbass
    read(irf2) mvib,nkbas,lmin,lbasis

    ivib=mvib(1)
    do 232 i=2,nblk
    ivib=max(ivib,mvib(i))
232   continue
    endif
 
!     determine  and print basis set labels
 
    write(6,"(//5x,' basis functions selected')")
    if (kmin == 0) then
        kz=1
    else
        kz=0
    endif
    loff0=mvib(1)*mvib(2)
    loff=loff0
    ipu=mvib(1)
    write(6,"(5x,'k =',i3,', i runs from',i5,' to',i5,2x,a4)") kz,1,ipu,symm(mod(ipar+kz,2)+1)
    kz=kz+1
    ipd=ipu
    ipu=ipd+mvib(2)
    write(6,"(5x,'k =',i3,', i runs from',i5,' to',i5,2x,a4)") kz,ipd+1,ipu,symm(mod(ipar+kz,2)+1)
    do 310 ioff=3,nblk
    kz=kz+1
    leng=mvib(ioff-1)*mvib(ioff)
    loff=loff+leng
    loff0=max(loff0,leng)
    leng=mvib(ioff-2)*mvib(ioff)
    loff=loff+leng
    loff0=max(loff0,leng)
    ipd=ipu
    ipu=ipd+mvib(ioff)
    write(6,"(5x,'k =',i3,', i runs from',i5,' to',i5,2x,a4)") kz,ipd+1,ipu,symm(mod(ipar+kz,2)+1)
    if(ipu==nbass) exit
310 continue

    nblk = kz + 1

    nvib=ivib
    if (zdcore) then
        keval=nbass
        lwork=max(loff0,3*nbass)
    else
        keval=min(max(neval,neval2)+4,nbass)
        lwork=3*keval+max(keval*keval,nbass+nbass)
    endif
    call vrmain(mvib,nkbas,lmin,lbasis,idvr)

    return
!     error on unit ivec
900 rewind ivec
    i=1
920 read(ivec,end=930)
    i=i+1
    goto 920
930 write(6,"(/'   End of stream IVEC =',i3,' at record',i6)") ivec,i
    stop
    end

!#####################################################################
    subroutine vrmain(mvib,nkbas,lmin,lbasis,idvr)

!     subroutine vrmain is the 'real' main program & contains       #006
!     the calls to the various subroutines which set & solve hamil
!
    use rotlev3b_seg_outp
    use rotlev3b_seg_size

    implicit double precision (a-h,o-y), logical (z)


    DIMENSION MVIB(NBLK),nkbas(NBLK), lmin(NBLK),lbasis(NBLK)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG,eval
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vec
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: radmee,radmoo,radmeo,radmoe
    data x0/0.0d0/

    if (abs(ires) == 2) goto 100
    if (ires < 0) then
        allocate(eval(neval))
        goto 50
    endif
    if (ires == 0) then
        ALLOCATE(radmee(maxblk),radmoo(mxblk2),radmeo(mxblk2),radmoe(mxblk2))

!        2d symmetrised radial matrix elements from 1d componants

        call radint(radmee,radmoo,radmeo,radmoe)

!        set up the hamiltonian matrix (not for restart runs)
        call solrt(radmee,radmoo,radmeo,radmoe,mvib)

        write(6,"(/5x,'hamiltonian construction complete')")
!         call timer
	call nftim('end of hamiltonian construction')
    DEALLOCATE(radmee,radmoo,radmeo,radmoe)
    endif
 
    if (.not.zdiag) stop
    if (ires >= 0) then
        ezero=x0
        read(5,"(f20.0)",end=555) ezero
555    continue
 
!        load the hamiltonian matrix elements
 
        ibass=nbass
        if (zdcore) then
            keval=nbass
            lwork=max(loff0,3*nbass)
            allocate(vec(nbass,nbass),diag(lwork),eval(keval))
        else
            keval=min(ibass,neval+4)
            allocate(vec(nbass,keval),diag(nbass+loff),eval(keval))
        endif

         all loadh(diag,mvib,vec,1)
 
!        diagonalise the hamiltonian (twice if requested)
 
        if (kmin == 0) write(6,1000) jrot,ibass
 1000 format('1'/5x,'J =',i3,' rotational state,',i7,' basis functions'&
        /12x,'f parity, anti-symmetric |jk> - |j-k> functions in basis')
        if (kmin >= 1) write(6,1010) jrot,ibass
 1010 format('1'/5x,'J =',i3,' rotational state,',i7,' basis functions'&
        /12x,'e parity, symmetric |jk> + |j-k> functions in basis')
        if (ipar == 0) write(6,"(12x,'even parity radial functions in basis set')")
        if (ipar == 1) write(6,"(12x,'odd parity radial functions in basis set')")
 
        call dgrot(diag,mvib,eval,vec,1,ezero,lwork)
 
        write(6,"(/5x,'diagonalisation complete')")
!         call timer
	call nftim('end of diagonalisation')
    endif
 
!     transform the eigenvectors if requested
    if  (ztran .or. zdcore) deallocate(vec,diag)  
50 if (ztran) then 
        call dstore(mvib,1,nkbas,lmin,lbasis,eval,idvr)
        call nftim('end of transformation')
    endif
    if (kmin<=1 .or. abs(ires)==3) goto 200
 
!     diagonalise/transform a second time if kmin > 1
 
100  ibass=nbass-mvib(1)
    neval=min(neval2,ibass)
    if (.not. zdcore) keval=min(ibass,neval+4)
    if (      zdcore) keval=ibass
    if (abs(ires)==2) allocate(eval(keval))
    jvec=jvec2
    ipar=1-ipar
    kvec=kvec2


    if (jrot==1) then
!     J=1f: treat as a special case
        write(6,1001) jrot
1001 format('1'/5x,'J =',i3,' rotational state special case for'&
        /12x,'f parity, anti-symmetric |jk> - |j-k> functions in basis')
        if (ipar == 0) write(6,"(12x,'even parity radial functions in basis set')")
        if (ipar == 1) write(6,"(12x,'odd parity radial functions in basis set')")
        call dstore1(2,eval,idvr,ezero)
        goto 200
    endif

    if (ires >= 0) then
        write(6,1000) jrot,ibass
        if (ipar == 0) write(6,"(12x,'even parity radial functions in basis set')")
        if (ipar == 1) write(6,"(12x,'odd parity radial functions in basis set')")
 
!        re-load the hamiltonian matrix elements
!        (for a restart run, first position the file)
        if (ires==2 .or. ztran .or. zdcore) then
            if (zdcore) then
               lwork=max(loff0,3*nbass)
               allocate(vec(nbass,nbass),diag(lwork))
            else
               allocate(vec(nbass,keval),diag(nbass+loff))
            endif
        endif 
        call loadh(diag,mvib,vec,2)
        nblk=nblk-1
        call dgrot(diag,mvib(2),eval,vec,2,ezero,lwork)

        write(6,"(/5x,'diagonalisation complete')")
!        call timer
	call nftim('transformation')
    else
        nblk=nblk-1
    endif
    deallocate(diag,vec)
        
!     transform the eigenvectors if requested
    if (ztran) then 
            call dstore(mvib(2),2,nkbas(2),lmin(2),lbasis(2),eval,idvr)
            call nftim('end of transformation')
    endif
200   deallocate(eval)
    return
    end


!#######################################################################
    subroutine radint(radmee,radmoo,radmeo,radmoe)
 
!     subroutine radint calculates the two-dimensional radial basis
!     functions between two symmetrised orthogonal coordinates.
    use rotlev3b_seg_outp
    use rotlev3b_seg_size

    implicit double precision (a-h,o-y), logical (z)

    DOUBLE PRECISION, DIMENSION(nmax) :: rm2
    DOUBLE PRECISION, DIMENSION(maxblk) :: radmee
    DOUBLE PRECISION, DIMENSION(mxblk2) :: radmoo
    DOUBLE PRECISION, DIMENSION(mxblk2) :: radmeo
    DOUBLE PRECISION, DIMENSION(mxblk2) :: radmoe
 
!     first read 1d radial matrix elements from file
    read(ivec) rm2
 
!     then use this data to construct the radial matrices
    call mkrad(radmee,nmax,rm2,0,0,0)
    if (mxblk2 > 0) then
        call mkrad(radmoo,nmax,rm2,1,0,0)
        call mkrad(radmeo,nmax,rm2,1,1,0)
        call mkrad(radmoe,nmax,rm2,1,1,1)
    endif
 
    return
    end
    subroutine mkrad(radmat,nmax,rm2,iq,ip,is)
 
!     construct symmetrised matrix elements with the appropriate parity
 
    implicit double precision (a-h,o-y), logical (z)
    dimension radmat(*),rm2(*)
 
    symm  = dble(1-2*ip)
    sign  = dble(1-2*is)
 
    ipt=0
    do 10 i1=1,nmax
    do 20 i2=1,i1-iq
    ipt=ipt+1
    radmat(ipt) = sign * (rm2(i2) + symm*rm2(i1))
20 continue
10 continue
    return
    end

!#########################################################################
    subroutine angin1(angmat,xalf,walf,pleg1,pleg2,mn,k1,&
                    plega,plegb,iv1,iv2,nang1,nang2,angfac)
 
!     angin1 calculates the angular integral between blocks k and k+1
!     required for the coriolis coupling term for orthogonal (radau)
!     coordinates in the bisector embedding. this is done by first
!     evaluating the matrix element in an fbr using gaussian quadrature
!     and then transforming to the appropriate dvrs.
    use rotlev3b_seg_size

    implicit double precision (a-h,o-y), logical (z)
 

    DIMENSION iv1(ndvr),iv2(ndvr)
    DOUBLE PRECISION, DIMENSION(ndvr,ndvr) :: fbrmat,pleg1,pleg2
    DOUBLE PRECISION, DIMENSION(iang,iang) :: angmat
    DOUBLE PRECISION, DIMENSION(npnt) :: xalf,walf
    DOUBLE PRECISION, DIMENSION(nang1,npnt) :: plega
    DOUBLE PRECISION, DIMENSION(nang2,npnt) :: plegb
    data x0/0.0d0/,xp5/0.5d0/,x1/1.0d0/,x2/2.0d0/,toler/1.0d-8/
 
!     first: set up an npnt gauss-associated legendre quadrature scheme
    realk = dble(k1)
    npnt2 = (npnt+1)/2
!     tswalf is the ]xact sum of weights for gauss-jacobi integration
    tswalf= x2**(k1+k1+1)/dble(k1+1)
    do 10 i=1,k1
    tswalf=tswalf*dble(i)/dble(k1+i+1)
10 continue
    call jacobi(npnt,npnt2,xalf,walf,realk,realk,cswalf,tswalf)
    write(6,1000) npnt,k1,(xalf(i),walf(i),i=1,npnt2)
 1000 format(//i8,' point gauss-associated legendre integration with ',&
             'k =' ,i3//5x,'integration points',11x,'weights',&
              /(f23.15,d25.12))
      write(6,1010) cswalf,tswalf
1010 format(/5x,'computed sum of weights',d26.15,&
            /5x,'exact    sum of weights',d26.15)
    if (abs((cswalf-tswalf)/tswalf) > toler) then
        write(6,"(//5x,'points & weights in error, adjust algorithm'//)")
        stop
    endif
    do 20 i=1,npnt2
    xalf(npnt+1-i) = xalf(i)
    xalf(i)        =-xalf(i)
    walf(npnt+1-i) = walf(i)
20 continue
!     evaluate the polynomials at the quadrature points
    call asleg(plega,fbrmat,nang1-1,xalf,npnt,k1)
!     return if the matrix elements are not actually needed
    if (mn == 0) return
    call asleg(plegb,fbrmat,nang2-1,xalf,npnt,k1+1)
!     compute the fbr matrix elements
    fbrmat = x0
    do 30 i=1,npnt
    term = walf(i) * (x1 + xalf(i))
    do 40 j=1,nang1
    fact = term * plega(j,i)
    do 50 k=1,nang2
    fbrmat(j,k) = fbrmat(j,k) + fact * plegb(k,i)
50 continue
40 continue
30 continue 
!     final step: transform the fbr matrix elements to the dvr
 
    kkp1=k1*(k1+1)
    xkph=dble(k1)+xp5
    i1=0
    do 60 i= 1,nang1
    if (iv1(i) == 0) goto 60
    i1=i1+1
    i2=0
    do 65 ip=1,nang2
    if (iv2(ip) == 0) goto 65
    i2=i2+1
    sum1=x0
    sum2=x0
    jj=k1
    do 70 j =1,nang1
    if (j > 1) then
        jj=jj+1
        sum1 = sum1 + pleg1(j,i) * pleg2(j-1,ip)&
            * sqrt(dble(jj*(jj+1)-kkp1))
    endif
    do 75 jp=1,nang2
    sum2 = sum2 + pleg1(j,i) * pleg2(jp,ip) * fbrmat(j,jp)
75 continue
70 continue
    angmat(i1,i2) = (sum1 + xkph*sum2) * angfac
65 continue
60 continue
    return
    end

!########################################################################
    subroutine angin2(angmat,xalf,walf,pleg1,pleg2,k1,&
                    plega,plegb,iv1,iv2,nang1,nang2,angfac)
 
!     angin2 calculates the angular integral between blocks k and k+2
!     required for the coriolis coupling term for orthogonal (radau)
!     coordinates in the bisector embedding. this is done by first
!     evaluating the matrix element in an fbr using gaussian quadrature
!     and then transforming to the appropriate dvrs.
    use rotlev3b_seg_size

    implicit double precision (a-h,o-y), logical (z)


    DIMENSION iv1(ndvr),iv2(ndvr)
    DOUBLE PRECISION, DIMENSION(ndvr,ndvr) :: fbrmat,pleg1,pleg2
    DOUBLE PRECISION, DIMENSION(iang,iang) :: angmat
    DOUBLE PRECISION, DIMENSION(npnt) :: xalf,walf
    DOUBLE PRECISION, DIMENSION(nang1,npnt) :: plega
    DOUBLE PRECISION, DIMENSION(nang2,npnt) :: plegb 
 
    data x0/0.0d0/,x1/1.0d0/
!     evaluate the polynomials at the quadrature points
    call asleg(plegb,fbrmat,nang2-1,xalf,npnt,k1+2)
!     now compute the fbr matrix elements
    fbrmat = x0
    do 30 i=1,npnt
    term = walf(i) * (x1 + xalf(i))**2
    do 40 j=1,nang1
    fact = term * plega(j,i)
    do 50 k=1,nang2
    fbrmat(j,k) = fbrmat(j,k) + fact * plegb(k,i)
50 continue
40 continue
30 continue 
!     second step: transform the fbr matrix elements to the dvr
 
    i1=0
    do 60 i= 1,nang1
    if (iv1(i) == 0) goto 60
    i1=i1+1
    i2=0
    do 65 ip=1,nang2
    if (iv2(ip) == 0) goto 65
    i2=i2+1
    sum=x0
    do 70 j =1,nang1
    do 75 jp=1,nang2
    sum = sum + pleg1(j,i) * pleg2(jp,ip) * fbrmat(j,jp)
    75 continue
    70 continue
    angmat(i1,i2) = sum * angfac
    65 continue
    60 continue
    return
    end

!#################################################################
    subroutine gasleg(nn,nn2,x,a,alf,bta,b,c,csa,tsa)
 
!     calculates zeros x(i) of the nn'th order associated legendre
!     polynomial pn(alf,bta) for the segment (-1,1) & and corresponding
!     weights for gauss-jacobi integration. this routine, and those
!     which follow, are due to a.h.stroud and d. secrest, "gaussian
!     integration formulas", 1966, prentice hall, page 29.
!     note that for our purposes, alf= bta= nu.
 
    implicit double precision(a-h,o-z)
    DOUBLE PRECISION, DIMENSION(nn) :: x,a,b,c
    data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,eps/1.0d-12/,&
        x3/3.0d0/,x4/4.0d0/,x8/8.0d0/
    fn= dble(nn)
    csa= x0
    c(1) = 0.0d0
    b(1) = 0.0d0
    do 21 i=2,nn
    xi= dble(i)
    b(i) = 0.0d0
    c(i)= x4*(xi-x1)*(alf+xi-x1)*(bta+xi-x1)*(alf+bta+xi-x1)/&
            ((alf+bta+x2*xi-x1)*(alf+bta+x2*xi-x2)*&
            (alf+bta+x2*xi-x2)*(alf+bta+x2*xi-x3))
21    continue
    cc=tsa
    do 1 i=2,nn
    cc= cc*c(i)
 1    continue
    do 12 i=1,nn2
    if (i == 1) then
!        largest zero
        an= alf/fn
        bn= bta/fn
        r1= (x1 + alf)*(2.78d0/(x4 + fn*fn) +0.768*an/fn)
        r2= x1 + 1.48d0*an + 0.96d0*bn + 0.452*an*an + 0.83d0*an*bn
        xt= x1 - r1/r2
    else if (i == 2) then
!        second zero
        r1= (4.1d0 + alf)/((x1 + alf)*(x1 + 0.156*alf))
        r2= x1 + 0.06d0*(fn - x8)*(x1 + 0.12d0*alf)/fn
        r3= x1 + 0.012*bta*(x1 + abs(alf)/x4)/fn
        ratio= r1*r2*r3
        xt= xt - ratio*(x1 - xt)
    else if (i == 3) then
!        third zero
        r1= (1.67d0 + 0.28d0*alf)/(x1 + 0.37d0*alf)
        r2= x1 + 0.22d0*(fn - x8)/fn
        r3= x1 + x8*bta/((6.28d0 + bta)*fn*fn)
        ratio= r1*r2*r3
        xt= xt - ratio*(x(1) - xt)
    else
!        middle zeros
        xt= x3*x(i-1) - x3*x(i-2) + x(i-3)
    endif
 
    call root(xt,nn,alf,bta,dpn,pn1,b,c,eps)
    x(i)= xt
    a(i)= cc/(dpn*pn1)
    csa= csa + a(i)*x2
12  continue
    if (2*nn2 /= nn) csa = csa - a(nn2)
    return
    end


!############################################################################
    subroutine jacobi(nn,nn2,x,a,alf,bta,csa,tsa)

!     calculates zeros x(i) of the nn'th order jacobi polynomial
!     pn(alf,bta) for the segment (-1,1) & and corresponding weights
!     for gauss-jacobi integration. this routine uses a brute force
!     search for zeros due to GJ Harris (2001).
!     note that for our purposes, alf= bta= nu.


    implicit real*8(a-h,o-z)
    real*8, dimension(nn) :: x,a,b,c,xt
    data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/,& 
        eps/1.0d-12/,xstep/1.0d-6/
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
        xt(ii)=xxx-0.5*xstep
    endif

    if (ii == nn2) then
        do 40 i=1,nn2
        call recur(ptemp,dp,pn1,xt(i),nn,alf,bta,b,c)
40      continue
    else
        xxx=xxx-xstep
        if (xxx > -1.5*xstep) goto 30
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

!**********************************************************************c
!                                                **020
    subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)
 
!     improves the approximate root x; in addition obtains
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x.
 
    implicit double precision(a-h,o-z)
    DOUBLE PRECISION, DIMENSION(nn) :: b,c
    iter= 0
1     iter= iter + 1
    call recur(p,dp,pn1,x,nn,alf,bta,b,c)
    d = p/dp
    x = x - d
    if(abs(d) <= eps) goto 3
    if (iter < 10) goto 1
3     dpn= dp
    return
    end
!**********************************************************************c
 
    subroutine recur(pn,dpn,pn1,x,nn,alf,bta,b,c)
    implicit double precision(a-h,o-z)
    DOUBLE PRECISION, DIMENSION(nn) :: b,c
    data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
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
1     continue
    pn= p
    dpn= dp
    pn1= p1
    return
    end

    subroutine asleg(pleg,pnorm,lmax,x,nn2,m)
 
!     calculate polynomials 1 to lmax at x = cos(theta) for m = 0 or 1,
!     using the routine of press et al, numerical recipes, p. 182,
!     for the polynomial part of associated legendre functions.
!     a factor of sin(theta)**m has been removed from all functions;
!     this enables us to use jacobi integration with alf = bta = m,
!     using routines refived from beidenharn and louck.
 
    implicit double precision (a-h,o-z)
    dimension pleg(0:lmax,nn2),x(nn2),pnorm(0:lmax)
    data x1/1.0d0/,x2/2.0d0/

    if (m < 0) goto 999
    do 10 i=1,nn2

      !For high J and high k the value pmm can become too large over this loop
      !Therefore in these instances we need pmm to be very small to begin with
      !We later divide by this factor to achieve the originally required number
      !For so2 at J = 200 and k > 90 we find this necessary, so we add this conditional
	if(m<90)then
	  pmm = x1
	else
	  pmm = 1d-250
	end if

    fact = x1
    do 11 j=1,m
    pmm = -pmm * fact
    fact = fact + x2
11 continue

      !Even when we use the above trick to reduce the size of pmm, it can still be too large for subsequent loops
      !We therefore conditionally divide pmm by a large enough factor - we take this to be for k > 150 (for so2)
      !write(*,*) pmm
      if(m>150)then
	pmm = pmm/1.0d200
    end if
 
    pleg(0,i) = pmm
    pmmp1= x(i)*(m+m+1)*pmm
    pleg(1,i)= pmmp1

    ll=1
    do 2 l= 2+m,lmax+m
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
    jj = -1
    do 13 j = m,lmax+m
    fact = x1
    do 12 i = j-m+1,j+m
  
!After a certain point, square-rooting this factor doesn't reduce the number sufficiently to a coping level
!Therefore we take the fourth root instead. We set the cut-off point to be the same as that for pmm, k > 90
	if(m<90)then
	   facti = sqrt(dble(i))
	else
	   facti = dble(i)**(1.0d0/4.0d0)
	end if
      fact = fact * facti

12 continue
    rj = dble(j)
    jj = jj + 1

    if(m<90)then
	   pnorm(jj) = sqrt((rj + rj + x1) /2)/fact
	else
	   pnorm(jj) = (((rj + rj + x1) /2)**(1.0d0/4.0d0))/fact
    end if
!      pnorm(jj) = sqrt(dble(j+j+1) / (fact + fact))

13 continue
!     now normalise the polynomials
    do 14 i=1,nn2
	do 15 jj=0,lmax
	  pleg(jj,i) = pleg(jj,i) * pnorm(jj)
	   !This is where we multiply by the large factor initially divided from pmm at the beginning of its loop
	    if(m>=90)then
	      pleg(jj,i) = pleg(jj,i) * 1.0d250 * pnorm(jj)
	    end if
	    
	    !And for the k's which are greater than 150, we need to divide by the reducing factor again
	    if(m>150)then
	      pleg(jj,i) = pleg(jj,i) * 1.0d200
	    end if
15 continue
14 continue
      
    return
999   write(6,"(//5x,'improper argument in subroutine asleg'/)")
    stop
    end

    subroutine wrtho(diag,offdg,mvib,nbass,nblk)
!     print hamiltonian matrix                                      #009
    implicit double precision (a-h,o-y)
    dimension diag(nbass),offdg(*),mvib(nblk)
    write(6,1010) diag
 1010 format('1',5x,'hamiltonian matrix: diagonal elements',&
            /(10f13.8))
 
    ioff1=1
    ioff2=mvib(1)*mvib(2)
    num=1
    write(6,"(//5x,'off-diagonal block number',i3/(10f13.8))") num,(offdg(i),i=ioff1,ioff2)
    do 10 k=3,nblk
    ioff1=ioff2+1
    ioff2=ioff2+mvib(k-2)*mvib(k)
    num=num+1
    write(6,"(//5x,'off-diagonal block number',i3/(10f13.8))") num,(offdg(i),i=ioff1,ioff2)
    ioff1=ioff2+1
    ioff2=ioff2+mvib(k-1)*mvib(k)
    num=num+1
    write(6,"(//5x,'off-diagonal block number',i3/(10f13.8))") num,(offdg(i),i=ioff1,ioff2)
10 continue
    return
    end
    subroutine wrthi(hamil,nham)
!     print hamiltonian matrix in core version
    double precision hamil(nham,nham)
    write(6,"(5x,'hamiltonian matrix'/)")
    do 30 i=1,nham
    write(6,"(10f13.7)") (hamil(i,j),j=1,i)
30 continue
    return
    end

    subroutine loadh(diag,mvib,hamil,itime)
 
!     subroutine loadh loads the hamiltonian matrix from disk
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)

    DIMENSION MVIB(NBLK)
    DOUBLE PRECISION, dimension(*) :: diag
    DOUBLE PRECISION, dimension(ibass,ibass) :: hamil
    data x0/0.0d0/
    if (zdcore) then
        hamil=x0
        ist2=0
    endif
 
    ipt=1+ibass
    if (itime > 1) then
!     first the blocks involving k=1 from the end of iscr
        if (.not. zdcore) then
            leng=mvib(3)*mvib(2)
            if (leng > 0)  call getrow(diag(ipt),leng,iscr)
            ipt=ipt+leng
            if (nblk > 3) then
               leng=mvib(4)*mvib(2)
               if (leng > 0) call getrow(diag(ipt),leng,iscr)
               ipt=ipt+leng
            endif
!           then the diagonal elements
            call getrow(diag,ibass,iscr)
        else
            ist2=ist2+mvib(2)
            leng=mvib(3)*mvib(2)
            if (leng > 0) call getrow(diag,leng,iscr)
            ipt=0
            do 10 i1=1,mvib(2)
            do 20 i2=ist2+1,ist2+mvib(3)
            ipt=ipt+1
            hamil(i1,i2)=diag(ipt)
            hamil(i2,i1)=diag(ipt)
   20       continue
   10       continue
            if (nblk > 3) then
               ist3=ist2+mvib(3)
               leng=mvib(4)*mvib(2)
               if (leng > 0) call getrow(diag,leng,iscr)
               ipt=ipt+leng
               ipt=0
               do 30 i1=1,mvib(2)
               do 40 i2=ist3+1,ist3+mvib(4)
               ipt=ipt+1
               hamil(i1,i2)=diag(ipt)
               hamil(i2,i1)=diag(ipt)
   40          continue
   30          continue
            endif
!           then the diagonal elements
            call getrow(diag,ibass,iscr)
            do 60 i=1,ibass
            hamil(i,i)=diag(i)
   60       continue
        endif
    endif
!     reposition iscr (skipping k=0 and k=1 blocks if already read)
    rewind iscr
    if (itime > 1) then
        if (mvib(1)*mvib(2) > 0) read(iscr)
        if (mvib(1)*mvib(3) > 0) read(iscr)
        if (mvib(2)*mvib(3) > 0) read(iscr)
        if (nblk > 3) then
        if (mvib(2)*mvib(4) > 0) read(iscr)
        end if
        j0=4
    else
        j0=2
    endif
 
    if (.not. zdcore) then
        do 100 j=j0,nblk
        leng=mvib(j)*mvib(j-1)
        if (leng > 0) call getrow(diag(ipt),leng,iscr)
        ipt=ipt+leng
        if (j < nblk) then
            leng=mvib(j+1)*mvib(j-1)
            if (leng > 0) call getrow(diag(ipt),leng,iscr)
            ipt=ipt+leng
        endif
  100    continue
        if (itime <= 1) call getrow(diag,nbass,iscr)
 
!     print hamiltonian matrix- if requested
 
        if (zpham)& 
            call wrtho(diag,diag(ibass+1),mvib(itime),ibass,nblk+1-itime)
    else
        do 150 j=j0,nblk
        ist1=ist2
        ist2=ist2+mvib(j-1)
        leng=mvib(j)*mvib(j-1)
        if (leng > 0) call getrow(diag,leng,iscr)
        ipt=0
        do 130 i1=ist1+1,ist1+mvib(j-1)
        do 135 i2=ist2+1,ist2+mvib(j)
        ipt=ipt+1
        hamil(i1,i2)=diag(ipt)
        hamil(i2,i1)=diag(ipt)
  135    continue
  130    continue
        if (j < nblk) then
            ist3=ist2+mvib(j)
            leng=mvib(j+1)*mvib(j-1)
            if (leng > 0) call getrow(diag,leng,iscr)
            ipt=ipt+leng
            ipt=0
            do 140 i1=ist1+1,ist1+mvib(j-1)
            do 145 i2=ist3+1,ist3+mvib(j+1)
            ipt=ipt+1
            hamil(i1,i2)=diag(ipt)
            hamil(i2,i1)=diag(ipt)
  145       continue
  140       continue
        endif
  150    continue
        if (itime == 1) then
            call getrow(diag,ibass,iscr)
            do 160 i=1,ibass
            hamil(i,i)=diag(i)
  160      continue
        endif
!     PRINT HAMILTONIAN MATRIX- IF REQUESTED

        IF (ZPHAM) CALL WRTHI(hamil,ibass)

    endif
    return
    end

!#########################################################################
    subroutine solrt(radmee,radmoo,radmeo,radmoe,mvib)
 
!     subroutine solrt sets up non-zero parts of hamiltonian        #010
!     including the computation of the k dependent angular matrix
!     elements
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)
 
    DIMENSION MVIB(NBLK)
    DIMENSION iv1(ndvr),iv2(ndvr)
    DOUBLE PRECISION, DIMENSION(npnt) :: xalf,walf
    DOUBLE PRECISION, DIMENSION(maxblk) :: radmee
    DOUBLE PRECISION, DIMENSION(mxblk2) :: radmoo,radmeo,radmoe
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pleg1,pleg2,plega,plegb
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG,OFFDG
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COEF1,COEF2,angmat
 
    data x2/2.0d0/,x16/16.0d0/,sqrt2/1.4142135623731d0/

    ALLOCATE(COEF1(MBASS,NVIB),COEF2(MBASS,NVIB),angmat(iang,iang),&
            pleg1(ndvr,ndvr),pleg2(ndvr,ndvr),diag(nbass),offdg(loff0),&
            plega(ndvr,npnt),plegb(ndvr,npnt))

    open(unit=iscr,form='unformatted')
    rewind iscr
 
    jjp1 = jrot * (jrot+1)
    ip = ipar
 
!     position ket input file (and load first ket data)
!     for nblk == 2 also load bra file
  
    read(ivec)
    read(ivec)
    read(ivec) k1,maxleg,nang1
!     if k1=0 and we are doing an f parity calculation, read k=1 from
!     end of file
    if (k1 == 0 .and. kmin == 0) then
        read(ivec)
        call endiv(ivec,jrot+1)
 
        read(ivec) k1,maxleg,nang1
        call getrow(pleg1,nang1*ndvr,ivec)
        read(ivec) iang1,ibass1
        read(ivec) (iv1(i),i=1,nang1)
        if (mvib(1) > 0) then
            read(ivec)
            call getrow(diag(1),mvib(1),ivec)
            if (nblk <= 2)&
               call rdcoef(coef1,ibass1,mvib(1),mvib(1),ivec)
        endif
        call reseti(ivec)
    else
        call getrow(pleg1,nang1*ndvr,ivec)
        read(ivec) iang1,ibass1
        read(ivec) (iv1(i),i=1,nang1)
        read(ivec) meval2
        if (mvib(1) > 0) then
            call getrow(diag(1),mvib(1),ivec)
            if (nblk <= 2) then
               call rdcoef(coef1,ibass1,mvib(1),meval2,ivec)
            else
               call rdcoef(coef1,ibass1,0,meval2,ivec)
            endif
        else
               call rdcoef(coef1,ibass1,0,meval2+1,ivec)
        endif
    endif
    idpt=mvib(1)+1
 
    read(ivec)
    read(ivec) k2,maxleg,nang2
    call getrow(pleg2,nang2*ndvr,ivec)
    read(ivec) iang2,ibass2
    read(ivec) (iv2(i),i=1,nang2)
 
    read(ivec) meval2
    if (mvib(2) > 0) then
        call getrow(diag(idpt),mvib(2),ivec)
        idpt=idpt+mvib(2)
    else
        read(ivec)
    endif
    call rdcoef(coef2,ibass2,mvib(2),meval2,ivec)
!     position second vectors file if it is required
    if (nblk > 2) then
        open(unit=ivec2,form='unformatted',recordtype='segmented')
        rewind ivec2
        do 10 i=1,5
        read(ivec2)
   10    continue
    endif
 
!     start a new band of the hamiltonian matrix: loop over ket data
 
    do 300 ioff=2,nblk
!     read in next set of bra vibrational vectors & diagonal elements
    if (nblk > 2) then
        read(ivec2)
        read(ivec2) k1,maxleg,nang1
!     if k1=0 and we are doing an f parity calculation, read k=1 from
!     end of file
        if (k1 == 0 .and. kmin == 0) then
            read(ivec2)
            call endiv(ivec2,jrot+1)
 
            read(ivec2) k1,maxleg,nang1
            call getrow(pleg1,nang1*ndvr,ivec2)
            read(ivec2) iang1,ibass1
            read(ivec2) (iv1(i),i=1,nang1)
            read(ivec2)
            if (mvib(1) > 0) then
               call getrow(diag(1),mvib(1),ivec2)
               call rdcoef(coef1,ibass1,mvib(1),mvib(1),ivec2)
            endif
            call reseti(ivec2)
        else
            call getrow(pleg1,nang1*ndvr,ivec2)
            read(ivec2) iang1,ibass1
            read(ivec2) (iv1(i),i=1,nang1)
            read(ivec2) meval2
            read(ivec2)
 
            call rdcoef(coef1,ibass1,mvib(ioff-1),meval2,ivec2)
        endif
    endif
 
!     compute angular off diagonal elements (kp = k+1)
 
    mn = mvib(ioff-1)*mvib(ioff)
!     angular factor for the present off-diagonal block
      angfac = +sqrt(dble(jjp1-k1*(k1+1)))/x2
    if (k1 == 0) angfac = sqrt2 * angfac

    call angin1(angmat,xalf,walf,pleg1,pleg2,mn,k1,&
                plega,plegb,iv1,iv2,nang1,nang2,angfac)
    if (mn > 0) then
        if (ip == 0) then
            call solofd(mn,radmeo,angmat,nmax,1,1,0,offdg,iang,&
                        iang1,iang2,mvib(ioff-1),mvib(ioff),coef1,coef2,&
                        ibass1,ibass2)
        else
            call solofd(mn,radmoe,angmat,nmax,1,0,1,offdg,iang,&
                        iang1,iang2,mvib(ioff-1),mvib(ioff),coef1,coef2,&
                        ibass1,ibass2)
        endif
    endif
 
!     next block: kp = k+2 (skip if we are near end of matrix)
 
    if (ioff == nblk) goto 300
    read(ivec)
    read(ivec) k2,maxleg,nang2
    call getrow(pleg2,nang2*ndvr,ivec)
    read(ivec) iang2,ibass2
    read(ivec) (iv2(i),i=1,nang2)
    read(ivec) meval2
 
    if (mvib(ioff+1) > 0) then
        call getrow(diag(idpt),mvib(ioff+1),ivec)
        idpt=idpt+mvib(ioff+1)
    else
        read(ivec)
    endif
    call rdcoef(coef2,ibass2,mvib(ioff+1),meval2,ivec)
 
!     compute angular off diagonal elements (if needed)
 
    mn = mvib(ioff-1)*mvib(ioff+1)
      if (mn > 0) then
!        angular factor for the present off-diagonal block
        angfac = -sqrt(dble((jjp1-(k1+1)*(k1+2))*(jjp1-k1*(k1+1))))/x16
        if (k1 == 0) angfac = sqrt2 * angfac
        call angin2(angmat,xalf,walf,pleg1,pleg2,k1,&
                plega,plegb,iv1,iv2,nang1,nang2,angfac)
 
        if (ip == 0) then
            call solofd(mn,radmee,angmat,nmax,0,0,0,offdg,iang,&
                    iang1,iang2,mvib(ioff-1),mvib(ioff+1),coef1,coef2,&
                    ibass1,ibass2)
        else
            call solofd(mn,radmoo,angmat,nmax,1,0,0,offdg,iang,&
                    iang1,iang2,mvib(ioff-1),mvib(ioff+1),coef1,coef2,&
                    ibass1,ibass2)
        endif
    endif
    ip=1-ip
  300 continue
 
!     place diagonal elements at end of scratch file
 
    call outrow(diag,nbass,iscr)
    if (kmin < 2) goto 600
 
!     compute k=1 f blocks and place them at the end of the scratch file
 
    read(ivec)
    read(ivec) k1,maxleg,nang1
    call getrow(pleg1,nang1*ndvr,ivec)
    read(ivec) iang1,ibass1
    read(ivec) (iv1(i),i=1,nang1)
    if (mvib(2) <= 0 .or. nblk<=2) goto 500
    read(ivec)
    call getrow(diag(mvib(1)+1),mvib(2),ivec)
    call rdcoef(coef1,ibass1,mvib(2),mvib(2),ivec)
 
    call reseti(ivec)
!     k=2 data
    read(ivec)
    read(ivec) k2,maxleg,nang2
    call getrow(pleg2,nang2*ndvr,ivec)
    read(ivec) iang2,ibass2
    read(ivec) (iv2(i),i=1,nang2)
    read(ivec) meval2
    read(ivec)
    call rdcoef(coef2,ibass2,mvib(3),meval2,ivec)
 
!     compute angular off diagonal elements (kp = k+1)
 
    mn = mvib(2)*mvib(3)
    ip = 1-ipar
!     angular factor for the present off-diagonal block
    angfac = +sqrt(dble(jjp1-2))/x2
    call angin1(angmat,xalf,walf,pleg1,pleg2,mn,k1,&
                plega,plegb,iv1,iv2,nang1,nang2,angfac)
    if (mn > 0) then
        if (ip == 0) then
            call solofd(mn,radmeo,angmat,nmax,1,1,0,offdg,iang,&
                    iang1,iang2,mvib(2),mvib(3),coef1,coef2,&
                    ibass1,ibass2)
        else
            call solofd(mn,radmoe,angmat,nmax,1,0,1,offdg,iang,&
                    iang1,iang2,mvib(2),mvib(3),coef1,coef2,&
                    ibass1,ibass2)
        endif
    endif
 
!     next block: kp = k+2 (skip if we are near end of matrix)
 
    if (nblk <= 3) goto 500
    read(ivec)
    read(ivec) k2,maxleg,nang2
    call getrow(pleg2,nang2*ndvr,ivec)
    read(ivec) iang2,ibass2
    read(ivec) (iv2(i),i=1,nang2)
    read(ivec)
    read(ivec) meval2
 
    call rdcoef(coef2,ibass2,mvib(4),mvib(4),ivec)
 
!     compute angular off diagonal elements (if needed)
 
    mn = mvib(2)*mvib(4)
    if (mn > 0) then
!       angular factor for the present off-diagonal block
        angfac = -sqrt(dble((jjp1-6)*(jjp1-2)))/x16
        call angin2(angmat,xalf,walf,pleg1,pleg2,k1,&
                plega,plegb,iv1,iv2,nang1,nang2,angfac)
 
        if (ip == 0) then
            call solofd(mn,radmee,angmat,nmax,0,0,0,offdg,iang,&
                     iang1,iang2,mvib(2),mvib(4),coef1,coef2,&
                     ibass1,ibass2)
        else
            call solofd(mn,radmoo,angmat,nmax,1,0,0,offdg,iang,&
                     iang1,iang2,mvib(2),mvib(4),coef1,coef2,&
                     ibass1,ibass2)
        endif
    endif
 
!     place f block diagonal elements at end of scratch file
 
500 call outrow(diag(mvib(1)+1),nbass-mvib(1),iscr)
600 DEALLOCATE(pleg1,pleg2,plega,plegb,DIAG,OFFDG,angmat,COEF1,COEF2)
    return
    end
    subroutine rdcoef(coef,idim,mvib,mev,iv)
!     read first step vector array coef from unit iv skipping surplus vectors
    implicit double precision (a-h,o-y), logical (z)
    dimension coef(idim,max(mvib,1))
 
    do 10 i=1,mvib
    call getrow(coef(1,i),idim,iv)
10 continue
    do 20 i=mvib+1,mev
    read(iv)
20 continue
    return
    end
    subroutine endiv(iv,nskip)
!     move input unit iv to end to read k=1f block data
    do 10 j=1,nskip-1
    read(iv)
    read(iv)
    read(iv) mev
    read(iv)
    do 20 i=1,mev
    read(iv)
20 continue
    read(iv)
    read(iv)
    read(iv)
   10 continue
!     last time, don't read angular header record
    read(iv)
    read(iv)
    read(iv) mev
    read(iv)
    do 30 i=1,mev
    read(iv)
30 continue
    read(iv)
    return
    end
    subroutine reseti(iv)
!     reposition input unit iv after reading k=1f block data from end
    rewind iv
    do 10 j=1,5
    read(iv)
10 continue
    do 20 j=1,2
    read(iv)
    read(iv)
    read(iv)
    read(iv)
    read(iv)
    read(iv) mev
    read(iv)
    do 30 i=1,mev
    read(iv)
30 continue
20 continue
    return
     end

!########################################################################
    subroutine solofd(mn,radmat,angmat,nmax1,iq,iq1,iq2,offdg,iang11,& 
                iang1,iang2,mvib1,mvib2,coef1,coef2,ibass1,ibass2)

!     construct off-diagonal matrix elements from radial and angular
!     matrix elements and first step vectors.
!     New algorithm introduced. JT June 2012.
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)

    DOUBLE PRECISION, DIMENSION(mn) :: offdg
    DOUBLE PRECISION, DIMENSION(*) :: coef1,coef2
    DOUBLE PRECISION, DIMENSION(*) :: radmat
    DOUBLE PRECISION, DIMENSION(iang11,iang11) :: angmat
    DOUBLE PRECISION, DIMENSION(iang1,mvib2) :: pdg
    data x0/0.0d0/,x1/1.0d0/
    offdg = x0
    ir1=0
    ir2=0
    ir=0
    do 150 i1=1,nmax1
    do 200 i2=1,i1-iq
    i0=ir1*iang1+1
    ir1=ir1+1
    j0=ir2*iang2+1
    ir2=ir2+1
    ir=ir+1
!Ala'a Azzam modification 18 sep 2012

    pdg = x0
    call dgemm('N','N',iang1,mvib2,iang2,radmat(ir),angmat,iang11,&
                coef2(j0),ibass2,x1,pdg,iang1)

    call dgemm('T','N',mvib2,mvib1,iang1,x1,pdg,iang1,&
                coef1(i0),ibass1,x1,offdg,mvib2)
200 continue
    ir1=ir1+iq1
    ir2=ir2+iq2
150 continue
!     dump completed block to unit iscr
    call outrow(offdg,mn,iscr)
    return
    end
!##########################################################################
    subroutine dgrot(diag,mvib,eval,vec,k1,ezero,lwork)
 
!     subroutine diag solves the eigenvalue problem:                #012
!          hamil * vec = eval * vec
!     by using iterative nag routine f02fjf to do diagonalisations.
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)


    DOUBLE PRECISION, DIMENSION(NEVAL) :: EVAL
    DOUBLE PRECISION, DIMENSION(*) :: DIAG
    DOUBLE PRECISION, DIMENSION(IBASS,*) :: VEC
    DIMENSION MVIB(NBLK),IBIG(IBASS)

!          autocm converts atomic units (hartree) to cm-1.
    data autocm/2.19474624d+05/,x0/0.0d0/

    if (zdcore) then
        ifail=0
        call dsyev('V','U',ibass,vec,ibass,eval,diag,lwork,ifail)
        if (ifail /= 0) write(6,"(/5x,'LAPACK routine SSYEV returned INFO =',i5)  ") ifail 
    else
        LWORK=3*KEVAL+MAX(KEVAL*KEVAL,NBASS+NBASS)
        call dgiter(diag,mvib,eval,vec,lwork)
    endif
 
!     print eigenvalues in atomic units & wavenumbers
 
    write(6,"(//5x,'lowest',i4,' eigenvalues in Hartrees',/)") neval
    write(6,"(5d24.12/)") eval
    if (zpfun) then
        if (k1 == 1) then
            open(unit=ilev,form='formatted')
            rewind ilev
  200       read(ilev,*,end=210,err=210)
            goto 200
  210       continue
! ***** inclusion of the following card is machine dependent *****
            backspace ilev
        endif
        ip=1-kmin
        if (kmin > 1) ip=k1-1
        write(ilev,"(5i4,i7) ") jrot,ip,0,0,(2-4*ipar),neval
        write(ilev,"(4d20.12)") eval
    endif
    if (zvec) then
!        write eigenvalues, eigenvectors, etc to stream jvec
        kz=kmin
        if (kmin > 1) kz=2-k1
        open(unit=jvec,form='unformatted')
        rewind jvec
        write(jvec) jrot,kz,ipar,neval,ibass
        write(jvec) mvib
        call outrow(eval,neval,jvec)
        mend=0
        do 100 k=1,nblk
        mbeg=mend+1
        mend=mend+mvib(k)
        if(mvib(k)>0) write(jvec) ((vec(j,i),j=mbeg,mend),i=1,neval)
  100    continue
    endif
    do 60 i=1,neval
    eval(i) = eval(i) * autocm - ezero
   60 continue
    write(6,1010) neval,ezero
 1010 format(//5x,'lowest',i4,' eigenvalues in wavenumbers relative to',&
                  ' ezero =',d24.12/)
    write(6,"(5d24.12/)") eval

    if (zpvec) then
        if (thresh <= x0) then
!             print complete eigenvectors
              write(6,"(//'    eigenvectors',/)")
              do 70 i=1,neval
              write(6,"(/(1x,10f13.7))") (vec(j,i),j=1,ibass)
   70         continue
          else
!             print largest componants of the eigenvectors
              write(6,"(a100,f5.2)") 'eigenvector componants greater than thresh =',thresh
              do 80 i=1,neval
              vbig=x0
              ipt=0
              do 90 j=1,ibass
              vv=abs(vec(j,i))
              if (vv > thresh) then
                  ipt=ipt+1
                  ibig(ipt)=j
              endif
              if (ipt <= 0 .and. vv > vbig) then
                  vbig=vv
                  ibig(1)=j
              endif
   90         continue
              write(6,"('  vector',i3,5(i7,f14.10)/(11x,5(i7,f14.10)))") i,(ibig(j),vec(ibig(j),i),j=1,max(1,ipt))
 !1060         format('  vector',i3,5(i7,f14.10)/(11x,5(i7,f14.10)))
   80         continue
        endif
    endif
    return
    end

!####################################################################
    subroutine dgiter(diag,mvib,eval,vec,lwork)
 
!     subroutine dgiter solves the eigenvalue problem:
!          hamil * vec = eval * vec
!     by using iterative nag routine f02fjf to do diagonalisations.
    use rotlev3b_seg_outp

    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)

    double precision, external :: vecvec
    external matvec,f02fjz

    DOUBLE PRECISION, DIMENSION(KEVAL) :: EVAL
    DOUBLE PRECISION, DIMENSION(*) :: DIAG
    DOUBLE PRECISION, DIMENSION(IBASS,keval) :: VEC
    DIMENSION MVIB(NBLK)
    DOUBLE PRECISION, DIMENSION(lwork) :: WORK

    data x0/0.0d0/,x1/1.0d0/,emax/1.0d50/,noffd/1/
 
!     create some starting vectors by using the diagonal elements
    vec=x0
    esmall=-emax
    do 10 i=1,keval
    ebig=emax
    do 20 j=1,ibass
    if (diag(j) > ebig) goto 20
    if (diag(j) <= esmall) goto 20
    ind=j
    ebig=diag(j)
20 continue
    vec(ind,i)=x1
    esmall=ebig
10 continue
 
!     shift diagonals to ensure we get the lowest eigenvalules
 
    eshift=ebig
    do 30 i=1,ibass
    eshift=max(eshift,diag(i))
30 continue
    do 40 i=1,ibass
    diag(i)=diag(i)-eshift
40 continue
 
!     diagonalise  the hamiltonian

    ifail=1
    noits=10000
    ndiag=ibass+1
    call f02fjf(ibass,neval,keval,noits,toler,vecvec,matvec,f02fjz,&
                keval,vec,ibass,eval,work,lwork,diag,noffd,&
                mvib,ndiag,ifail)
    if (ifail /= 0) write(6,"(//5x,'f02fjf returned ifail =',i3)") ifail
 
    write(6,"(/5x,'convergence after noits =',i6,' iterations')") noits
!     shift eigenvalues back
    do 50 i=1,neval
    eval(i)=eval(i)+eshift
50 continue
    return
    end

!#######################################################################
    subroutine dstore(mvib,itra,nkbas,lmin,lbasis,energy,idvr)
 
!     dstore transforms the eigenvectors of the second variational
!     step into ones for the first step basis and stores the
!     results in a form suitable for program dipole3.
!     This version treats each k-block seperately.
      use rotlev3b_seg_outp
      use rotlev3b_seg_size
      implicit double precision(a-h,o-y), logical(z)

      DIMENSION MVIB(nblk),NKBAS(nblk),lmin(nblk),lbasis(nblk)
!      DOUBLE PRECISION, DIMENSION(3) :: XMASS
!      DOUBLE PRECISION, DIMENSION(mbass0,neval) :: d
!      DOUBLE PRECISION, DIMENSION(nvib*neval) :: c
      DOUBLE PRECISION, DIMENSION(KEVAL) :: ENERGY
!      DOUBLE PRECISION, DIMENSION(nmax) :: r
!      DIMENSION ivt(ndvr)
!      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: B

    DOUBLE PRECISION, allocatable :: XMASS (:)
    DOUBLE PRECISION, allocatable :: d (:,:)
    DOUBLE PRECISION, allocatable :: c (:,:)
    DOUBLE PRECISION, allocatable :: r (:)
    DIMENSION ivt(ndvr)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: B

    allocate(XMASS(3))
    allocate(d(mbass0,neval))
    allocate(c(nvib,neval))
    allocate(r(nmax))
      
    data x0/0.0d0/,x1/1.0d0/
 
    write(6,1000) ivec,jvec,kvec
1000 format('1'/5x,'eigenvector transformation:',&
            /5x,'input:   dvr3d data,          ivec =',i3,&
            /5x,'         rotlev3b data        jvec =',i3,&
            /5x,'output:  transformed vectors, kvec =',i3/)
!     read dvr3d header
    rewind ivec
    read(ivec) idia,ipar,idvr,npnt1,npnt2,jrot0,kmin0,mval
    read(ivec) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
    read(ivec) re1,diss1,we1,re2,diss2,we2
    read(ivec)
    read(ivec) r
!     read rotlev header
    rewind jvec
    read(jvec) jrot1,kmin1,ipar1,nval,ibass
    read(jvec) mvib
    nval=min(nval,neval)
!     check for compatability
    if (jrot1 /= abs(jrot0)) then
        write(6,900) jrot1,abs(jrot0)
900     format(/5x,'j levels mismatched',&
                 /5x,'jrot1 =',i3,'  jrot0 =',i3)
        stop
    endif
    if (kmin1 > kmin0) then
        write(6,"(/5x,'kmin1 and kmin0 incompatible')")
        stop
    endif
    if (itra == 2) ipar=1-ipar
    if (itra==1 .and. kmin0==2 .and. kmin1==0)  ipar=1-ipar
    if (ipar /= ipar1) then
        write(6,"(/5x,'parities mismatched'/5x,'ipar= ',i2,'  ipar1= ',i2)") ipar,ipar1
        stop
    endif
    if (zptra) write(6,"(/5x,'j =',i3,' kmin =',i2,'  parity =',i2,' nval =',i5/)") jrot1,kmin1,ipar,nval
!     size of final basis?
    nend=0
    nkmax=0
    do 20 k=1,nblk
    nend=nend+nkbas(k)
    nkmax=max(nkmax,nkbas(k))
   20 continue 
!     write header on new file
    open(unit=kvec,form='unformatted',recordtype='segmented')
    write(kvec) idia,ipar,idvr,npnt1,npnt2,jrot1,kmin1,nval
    write(kvec) zembed,zmors1,zmors2,xmass,g1,g2,zncor
    write(kvec) re1,diss1,we1,re2,diss2,we2
    write(kvec) nend,lmin,lbasis,nkbas
!     transfer the dvr points
    write(kvec) r
    write(kvec) nval
    call getrow(energy,nval,jvec)
    call outrow(energy,nval,kvec)
 
!     start transformation step
    allocate(b(nkmax,nvib))
    iiout=0
    do 50 k=1,nblk
!     transform dvr basis to fbr in theta for current k-block
    call tofbr(nkbas(k),mvib(k),ivt,b,itra,idvr)
    if (mvib(k) <= 0) goto 50
!     read in the untransformed vectors for current k-block
    call getrow(c,mvib(k)*nval,jvec)

!     the code below here performs a standard matrix multiply of the form
!     D = B.C
!     where the matrices are dimensioned d(nkbas(k),nval),
!     b(nkbas(k),mvib(k))  and c(mvib(k),nval). 

    d = x0
    call dgemm('N','N',nkbas(k),nval,mvib(k),x1,b,nkbas(k),c,mvib(k),&
                x0,d,nkbas(k))
    call outrow(d,nkbas(k)*nval,kvec)
        if (zptra) then
            write(6,"(/'Vectors for K-block',i3/)") k
            do 60 i=1,nval
            write(6,"(1X,10F13.7)") (d(j,i),j=1,nkbas(k))
   60       continue
        endif
   50 continue
    deallocate(b) 
    write(6,1050) jrot1,kmin1,ipar,nval,nend
 1050 format(//5x,'transformation completed successfully',&
            //5x,'for rotational state: jrot =',i3,', kmin =',i3,&
            ', ipar =',i3,'.  stored',i6,' vib-rot levels',&
            /5x,'nend =',i7)
    rewind jvec
    rewind kvec
    return
    end
!#######################################################################
    subroutine dstore1(itra,energy,idvr,ezero)
                        
!     `Transformation' step for J=1f special case      
!     RESULTS IN A FORM SUITABLE FOR program DIPOLE3.
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision(a-h,o-y), logical(z)

    DOUBLE PRECISION, DIMENSION(3) :: XMASS
    DOUBLE PRECISION, DIMENSION(mbass) :: d
    DOUBLE PRECISION, DIMENSION(KEVAL) :: ENERGY
    DOUBLE PRECISION, DIMENSION(nmax) :: r
    DOUBLE PRECISION, DIMENSION((ndvr+1)**2) :: pleg
    DIMENSION ivt(ndvr)
    data autocm/2.19474624d+05/
 
!     read dvr3d header
    rewind ivec
    read(ivec) idia,ipar,idvr0,npnt1,npnt2,jrot0,kmin0,mval
    read(ivec) zembed,zmors1,zmors2,xmass,g1,g2,zncor,zquad2
    read(ivec) re1,diss1,we1,re2,diss2,we2
    read(ivec)
    read(ivec) r
    kmin1=0
    if (itra == 2) ipar=1-ipar
    if (itra==1 .and. kmin0==2 .and. kmin1==0)  ipar=1-ipar
    if (ztran) then
!       write header on new file
        open(unit=kvec,form='unformatted',recordtype='segmented')
        write(kvec) idia,ipar,idvr,npnt1,npnt2,jrot0,kmin1,neval
        write(kvec) zembed,zmors1,zmors2,xmass,g1,g2,zncor
        write(kvec) re1,diss1,we1,re2,diss2,we2
    endif
!     transfer the dvr points
    nang=idvr
    max2d  = npnt1*(npnt1+1-2*ipar)/2
    mbass=idvr*max2d
    if (ztran) write(kvec) mbass,0,idvr,mbass
    if (ztran) write(kvec) r
!     if file contains e and f, we are doing an f parity calculation, 
!     read k=1 from  end of file
    if (kmin0 == 2) then
        read(ivec)
        read(ivec)
        read(ivec)
        call endiv(ivec,2)
        backspace ivec
    endif
    call getrow(pleg,idvr,ivec)
    if (ztran) call outrow(pleg,idvr,kvec)
    read(ivec)  k1,maxleg,nang1
    if (ztran) write(kvec) k1,maxleg,nang1
    call getrow(pleg,nang1*ndvr,ivec)
    if (ztran) call outrow(pleg,nang1*ndvr,kvec)
    read(ivec)  iang1,ibass1
    if (ztran) write(kvec) iang1,ibass1
    read(ivec)  (ivt(i),i=1,nang1)
    if (ztran) write(kvec) (ivt(i),i=1,nang1)
    read(ivec) meval0
    if (ztran) write(kvec) neval
    call getrow(energy,neval,ivec)
    if (ztran) call outrow(energy,neval,kvec)
 
!     print eigenvalues in atomic units & wavenumbers
 
    write(6,"(//5x,'lowest',i4,' eigenvalues in hartrees',/)") neval
    write(6,"(5d24.12/)") (energy(i),i=1,neval)
    if (zpfun) then
        if (itra == 1) then
            open(unit=ilev,form='formatted')
            rewind ilev
  200       read(ilev,*,end=210,err=210)
            goto 200
  210       continue
! ***** inclusion of the following card is machine dependent *****
            backspace ilev
        endif
        ip=0
        write(ilev,"(5i4,i7)") jrot,ip,0,0,(2-4*ipar),neval
        write(ilev,"(4d20.12)") (energy(i),i=1,neval)
    endif
    do 60 i=1,neval
    energy(i) = energy(i) * autocm - ezero
   60 continue
    write(6,1010) neval,ezero
 1010 format(//5x,'lowest',i4,' eigenvalues in wavenumbers relative to',&
                  ' ezero =',d24.12/)
    write(6,"(5d24.12/)") (energy(i),i=1,neval)
! 1020 format(5d24.12/)
    if (.not.ztran) return

    WRITE(6,2000) ivec,kvec
 2000 FORMAT(5X,'EIGENVECTOR DUMMY TRANSFORMATION J=1f case:',&
             /5X,'INPUT:  DVR3DRJZ data,         IVEC  =',I3,&
             /5X,'OUTPUT: transformed vectors,   KVEC  =',I3/)
    if (zptra) write(6,2010) jrot0,kmin1,ipar,neval
 2010 format(/5x,'j =',i3,' kmin =',i2,'  parity =',i2,' neval =',i5/&
              /'Vectors for K-block',i3/)

!     copy across DVR3DRJZ vectors
    do 50 k=1,neval
    call getrow(d,ibass1,ivec)
    call outrow(d,ibass1,kvec)
    if (zptra) write(6,"(1X,10F13.7)") (d(j),j=1,ibass1)
   50 continue
 
    write(6,1050) jrot0,kmin1,ipar,neval,ibass1
 1050 format(//5x,'Transfer of vectors completed successfully',&
            //5x,'for rotational state: jrot =',i3,', kmin =',i3,&
            ', ipar =',i3,'.  stored',i6,' vib-rot levels',&
            /5x,'ibass =',i7)
    close(unit=kvec)
    return
    end
!#######################################################################
    subroutine tofbr(nkbas,mvib,iv1,fbrvec,itra,idvr)
 
!     take the dvr vectors from unit ivec and transform them to fbr in
!     theta. also construct pointer arrays for dipole3.
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision(a-h,o-y), logical(z)


    DOUBLE PRECISION, DIMENSION(*) :: fbrvec
    DOUBLE PRECISION, DIMENSION(ndvr,idvr) :: pleg
    DIMENSION iv1(ndvr)
 
    read(ivec)
    read(ivec) kz,maxleg,idvr,lincr
!     if k=0 and we are doing an f parity calculation, read k=1 from
!     end of file
    if (kz == 0 .and. (kmin == 0 .or. itra == 2)) then
        read(ivec)
        call endiv(ivec,jrot+1)
 
        read(ivec) kz,maxleg,idvr,lincr
        call getrow(pleg,ndvr*idvr,ivec)
        read(ivec) iang1,ibass1
        read(ivec) (iv1(ii),ii=1,idvr)
        ivres=1
    else
        call getrow(pleg,ndvr*idvr,ivec)
        read(ivec) iang1,ibass1
        read(ivec) (iv1(ii),ii=1,idvr)
        ivres=0
    endif
    nrad=ibass1/iang1
    nang=maxleg+lincr-kz+1
    read(ivec) meval
    read(ivec)
!     read basis vectors and transform to associated legendres
    if (mvib > 0) &
        call jtran(fbrvec,mvib,pleg,idvr,nrad,nang,ibass1,iv1,nkbas)
    if (ivres == 0) then
!        skip vectors that are not needed
        do 106 l=mvib+1,meval
        read(ivec)
 106     continue
    else
        call reseti(ivec)
    endif
    return
    end

!#######################################################################
    subroutine jtran(coef,mvib,pleg,idvr,nrad,nang,ibass1,iv,nkbas)
    use rotlev3b_seg_outp
    use rotlev3b_seg_size
    implicit double precision (a-h,o-y), logical (z)

    DOUBLE PRECISION, DIMENSION(ndvr,idvr) :: pleg
    DOUBLE PRECISION, DIMENSION(iang,nrad) :: dvrvec
    DOUBLE PRECISION, DIMENSION(nkbas,mvib) :: coef
    DOUBLE PRECISION, DIMENSION(nrad) :: sumk
    DIMENSION iv(nang)
    data x0/0.0d0/

!     transform back to the original fbr-type basis in the
!     associated legendre functions

    do 10 l=1,mvib
!     first read in a new vector
    iang1=ibass1/nrad
    read(ivec)((dvrvec(i,j),i=1,iang1),j=1,nrad)
    do 20 j=1,nang
    sumk=x0
    kk=0
    do 40 k=1,idvr
    if (iv(k) == 0) goto 40
    kk=kk+1
    do 50 mn=1,nrad
    sumk(mn)=sumk(mn) + dvrvec(kk,mn) * pleg(j,k)
   50 continue
   40 continue
    ipt=j
    do 60 mn=1,nrad
    coef(ipt,l) = sumk(mn)
    ipt=ipt+nang
   60 continue
   20 continue
 
   10 continue
    return
    end

!######################################################################
    subroutine getrow(row,nrow,iunit)
!     fast non-formatted read                                       #015
    DOUBLE PRECISION, DIMENSION(NROW) :: ROW
    read(iunit) row
    return
    end

!######################################################################
    subroutine outrow(row,nrow,iunit)
!     fast non-formatted write                                      #016
    DOUBLE PRECISION, DIMENSION(NROW) :: ROW     
    write(iunit) row
    return
    end

!######################################################################
    FUNCTION VECVEC(IFLAG,N,Z,W,D1,D2,D3,D4)

!     VECVEC RETURNS THE DOT PRODUCT W.Z FOR F02FJF                 

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DOUBLE PRECISION, DIMENSION(N) :: W
    DOUBLE PRECISION, DIMENSION(N) :: Z

    VECVEC=DDOT(N,W,1,Z,1)
    RETURN
    END

!######################################################################
    subroutine matvec(iflag,n,z,w,hamil,noffd,mvib,ndiag)
 
!     matvec performs w = hamil * z for f02fjf                      #019
!     written to take advantage of the block structure
!     note:
!     hamil contains arrays diag & offdg relying on them being
!     adjacent in the dynamic store allocation
    use rotlev3b_seg_size
    implicit double precision (a-h,o-z)

    DOUBLE PRECISION, DIMENSION(IBASS) :: W,Z
    DOUBLE PRECISION, DIMENSION(*) :: HAMIL
    DIMENSION MVIB(NBLK)

    parameter (x1=1.0d0)
!     the diagonal diagonal contribution
    do 10 i=1,ibass
    w(i) = z(i) * hamil(i)
   10 continue
!     then the off-diagonal blocks
    ioff=ndiag
    i2=1
    do 20 k=1,nblk
    if (mvib(k) < 1) goto 20
!     (k,k-2) block
    if (k > 2 .and. mvib(k-2) > 0) then
        j1=i1-mvib(k-2)
        joff2=ioff-mvib(k-2)*mvib(k)
        call dgemv('N',mvib(k),&
                mvib(k-2),x1,hamil(joff2),mvib(k),z(j1),1,x1,w(i2),1)
    endif
!     (k,k-1) block
    if (k > 1 .and. mvib(k-1) > 0) call dgemv('N',mvib(k),&
                mvib(k-1),x1,hamil(ioff),mvib(k),z(i1),1,x1,w(i2),1)
    if (k>1.and.k<nblk) ioff=ioff2+mvib(k-1)*mvib(k+1)
    i1=i2
    i2=i2+mvib(k)
!     (k,k+1) block
    if (k < nblk .and. mvib(k+1) > 0) then
        call dgemv('T',mvib(k+1),&
                mvib(k),x1,hamil(ioff),mvib(k+1),z(i2),1,x1,w(i1),1)
        ioff2=ioff+mvib(k+1)*mvib(k)
    endif
!     (k,k+2) block
    if (k < nblk-1 .and. mvib(k+2) > 0) then
        j2=i2+mvib(k+1)
        call dgemv('T',mvib(k+2),&
                mvib(k),x1,hamil(ioff2),mvib(k+2),z(j2),1,x1,w(i1),1)
    endif
   20 continue
    return
    end

!###########################

    subroutine nftim(text)
    common/timing/itime0
    character text*(*)
    write(6,"(/)")
    write(6,*) 'Time at ',text,' is.........'
    call SYSTEM_CLOCK(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,' secs CPU time used'/)")itime
    return
    end
