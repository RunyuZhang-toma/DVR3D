!==================================================================================================
!Module defintion
!By Runyu Zhang & Tennyson Jonathan 1/Sep/2022
!Contains: size, timing, mass, dim
!Special notice:: The module name contains the file name, it means cannot directly paste to other 
!                 files and use. There are some difference between the constant numbers, types and 
!                 default value.
!==================================================================================================

module xpect3_logic
    logical :: zmors1   ! T use morse oscillator-like functions for r_1 coordinate;
    logical :: zncor = .false.
    logical :: zprint = .false.   ! T supplies extra print out for debugging purposes.  
    logical :: ztra = .true.    ! T writes out the data needed for program spectra 
                                   ! to calculate simulated spectra.
    logical :: zstart = .false.     ! T if we are writing out for spectra for the first time.
    logical :: zmors2       ! T use morse oscillator-like functions for r_2 coordinate;
                                   ! F use spherical oscillator functions.
    logical :: zfit = .false.
    logical :: zform = .true.
    integer :: iket = 11   ! input stream for the ket.
    integer :: itra = 12    ! output stream to program spectrm (if ztra).
                              ! note that for all times other than the dipole assumes 
                              ! that we have accessed the permanent dataset or file which has the
                              ! data from previous runs and that we are writing to the end of that file.
                              ! ************************************************
                              ! **  for the sake of safety you are therefore  **
                              ! **  advised to keep one previous edition as   **
                              ! **  backup!                                   **
                              ! ************************************************
    integer :: itra0 = 28
    integer :: ilev = 14   ! stream for final eigenvalues (formatted).

end module xpect3_logic

module xpect3_timing
    integer :: itime0
end module xpect3_timing

module xpect3_mass
    logical :: zembed   ! T z axis is along r2, = f z axis is along r1.
                                   ! only used if J > 0 ZBISC = in JHMAIN ie if zbisc=f and zperp=f.

    logical :: zbisc     ! T place the Z-axis along the bisector
    double precision :: xmass(3)
    double precision :: g1, g2
end module xpect3_mass

module xpect3_dim
    integer :: ncoord    ! number of vibrational coordinates explicitly considered
                              ! ncoord = 2: atom-diatom problem with diatom rigid
                              ! ncoord=2: also need lmax,lpot,idia,kmin
                              ! ncoord = 3: full 3-d triatomic problem
                              ! ncoord=3: all paramters required
    integer :: npnt          ! max(npnt1,npnt2) number of gauss-associated legendre grid points requested
    integer :: npnt1         ! number of (gauss-laguerre) dvr points in r1
    integer :: npnt2         ! number of (gauss-laguerre) dvr points in r2
    integer :: nrade
    integer :: nrado
    integer :: lpot
    integer :: npot    ! number of Gauss-Legendre integration points used
                              ! in i5 format
    integer :: nbmax1
    integer :: mbass  ! maximum size of vibrational problem (excluding linear geom)
    integer :: kmin1
    integer :: jk1
    integer :: neval
    integer :: nn2
    integer :: lmax
    integer :: npropin
    integer :: nprt
    integer :: jrot
    integer :: idia
    integer :: ipar ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
     integer :: nv1           ! number of bra eigenfunctions considered
                              ! if this is input as zero, all available
                              ! ket eigenfunctions will be considered when computing transitions.
                              ! in i5 format

end module xpect3_dim

!  Xpect takes expectation values of operators defined in prop.
!  Version considers diagonal case only. Nikolai Zobov (2002)
program leek
    call xpect3
    stop
end program leek
!**************************************************001
    subroutine xpect3
!
!     xpect3, and the various subroutines it calls, calculate the
!     expectation values of vib-rot eigenfunctions with respect
!     the properties defined by subroutine prop.
!     
    use xpect3_logic
    use xpect3_mass
    use xpect3_timing
    implicit double precision (a-h,o-y), logical (z)
    namelist/prt/ zprint,ztra,zstart,zfit,zform,&
                    iket,itra,itra0,ilev
    character(len=8) title(9)

    write(6,"(//5x,'program xpect3 (version of 19 March 2003):',/)")
!200   format(//5x,'program xpect3 (version of 19 March 2003):',/)
!
!     read in logical control parameters zprint,ztra,zstart,zfit,zform
!
!     zfit = .true. expectation values for potential fitting 
!                  (relative to the ground state)
!  *  zprint= .true. gives extra print out - useful for debugging.
!  *  ztra= .true. outputs the property data to strean itra
!  *  the default values of all three parameters are .false.
!     the variable "zform" controls whether the derivatives are outputted
!     formatted (true) or unformatted (false - default). (jhs)

    read(5,prt)
    read(5,"(9a8)") title

    write(6,"(/5x,9a8)") title

!     read control parameters
    call system_clock(itime0,irate2,imax2)
    call insize
!     now do the calculations

    call core

    call system_clock(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,'secs cpu time used'/)") itime     
    stop
end subroutine xpect3
      
!**************************************************002
      
!**************************************************003
subroutine insize
!
!     subroutine insize reads in the parameters which control the
!     overall size of the problem. some of these are passed from the
!     output streams of programs triatom and rotlev via the input
!     stream iket.
!     this means that many parameters are duplicated, giving the program
!     the opportunity to check that the bra and ket are compatible.
!     integration parameters for the two radial co-ordinates and the
!     angular functions must be set by the user, and are inputted on
!     stream 5 (five).
!
    use xpect3_logic
    use xpect3_mass
    use xpect3_timing
    implicit double precision (a-h,o-y),logical (z)
!
!     the following size parameters have these functions:
!
!     mbass: the total number of basis functions for the ket.
!     npnt1: number of gauss-laguerre dvr points for r1
!     npnt2: number of gauss-laguerre dvr points for r2
!     lmax : max number of angular basis functions in ket
!     lpot : max order of legendre polys in property expansion
!     npot : number of gauss-legendre integration points 
!     npnt : max(npnt1,npnt2)
!     idia : = 1 scattering coordinates hetronuclear diatomic
!            = 2 scattering coordinates homonuclear diatomic
!            = -1 radau coordinates hetronuclear diatomic
!            = -2 radau coordinates homonuclear diatomic
!     ipar : parity of basis - if |idia|=2: ipar=0 for even & =1 for odd
!          : taken from ipar1 and ipar2
!     kmin : kmin=1 for sym. rotational basis, =0 for anti-sym.
!          : for non-coriolis calculations, kmin= k.
!     neval: number of eigenvalues supplied from rotlev or triatom
!     ncoord: number of vibrational coordinates explicitly considered
!     if (ncoord /= 3) some of the above are dummies, see below.
!     nprop: number of properties to be considered
!     iprop: vector with the information on which properties
!            will be considered in the run
!
    parameter (mxprop=1000)
    common /pot/ iprop(mxprop)
!     save masses & g in case they are needed in the properties routine
!
!     read in control parameters of problem:
!     ncoord = 2: atom-diatom problem with diatom rigid
!     ncoord = 3: full 3-d triatomic problem
!
!     read in of data for the wavefunction
!
    if (zfit) write(6,"(/5x,'Property analysis for potential energy fitting')")
    if (.not. zfit) write(6,"(/5x,'Property analysis to obtain expectation values')")
    write(6,"(a100,i4)") 'input ket wavefunctions read from stream iket  =',iket
    if (ztra) then
        write(6,"(5x,a100,i4)") 'output transition data written to stream itra  =',itra
        if (zform) write(6,"(5x,'file written formatted')")
        if (.not.zform) write(6,"(5x,'file written unformatted')")
        if (.not. zfit) then
            if (zstart) write(6,"(5x,'file written from the start')")
            if (.not.zstart) write(6,"(5x,'data appended to the end of file')")
        endif
    endif

    open(unit=iket,form='unformatted')
    read(iket) idia,ipar,lmax,nmax11,nmax12,jrot,kmin1,neval
    read(iket) zembed,zmors1,zmors2,xmass,g1,g2,zncor
    read(iket) re1,diss1,we1,re2,diss2,we2

    ncoord=3
      if (idia > -2) then
        if (nmax11 == 1) ncoord=2
        zbisc=.false.
    else
        zbisc=.true.
        zembed=.true.
    endif

    if (zfit .and. max(jrot,ipar) > 0) write(6,1005) itra0
1005 format(  5x,'Ground state data read from       stream ',&
            'itra0 =',i4)
    write(6,1006) ilev
1006 format(  5x,'Ground state energy read from     stream ',&
            'ilev   ',i4)
      
!     set other parameters to maximum values

    npnt1= nmax11
    npnt2= nmax12
!
    npnt=max(npnt1,npnt2)
!
!     read in of integration parameters and states to be considered
!     zeros will give the default values
!
    read(5,108) lpot,npropin,nprt,nv1
    write(*,*)lpot,npropin,nprt,nv1,neval
108  format(5i5)

    if (nv1 > neval) then
        write(6,*)'NV1 > NEVAL! (NV1 MUST BE <= NEVAL)'
        write(6,*)'STOPPING!'
        stop
    else if (nv1 == 0) then
        nv1 = neval
    endif


    if (npropin == 0) then
        write(6,*)'NO PROPERTIES TO BE CONSIDERED! (NPROPIN = 0)'
        write(6,*)'STOPPING!'
        stop
    endif

if (npropin > mxprop) then
        write(6,9999) npropin,mxprop
9999    format(/,5x,'npropin =',i3,' requested, but mxprop =',i3,&
            /,5x,' ***** stop *****')
    endif
    if (nprt < npropin)then
        write(6,*)'WARNING NPROPIN SHOULD BE LESS THAN OR EQUAL TO NPRT!'
        write(6,*)'SETTING TO NPROPIN!'
        nprt=npropin
    endif

    if (lpot == 0) lpot = lmax+jrot+mod(lmax+jrot,2)
    if (jrot == 0) lpot = lmax

    iprop = 0
    read(5,"(20i5)") (iprop(i),i=1,npropin)

    if (iprop(1) <= 0) then
        do i=1,npropin
            iprop(i)=i 
        enddo
    else
        if (iprop(npropin)==0)then
            write(*,*)'WARNING NUMBER OF PROPERTIES CONSIDERED MUST BE EQUAL'
            write(*,'(A,I4)')' TO NPROPIN! SETTING TO 1, ...,',NPROPIN
        endif
    endif

    write(6,201) npropin,nprt,lpot,(iprop(i),i=1,npropin)
201   format(/,5x,'Number of properties to be considered   =',i5,&
            /,5x,'Number of properties computed in PROP   =',i5,&
            /,5x,'Number of angular DVR points            =',i5,&
            /,5x,'properties to be considered are numbered:',&
            /,5x,(25i5))

    if (zprint) write(6,*) '     Debug printing requested'
!
!     are we doing a non coriolis problem?
!     set up jk for the two cases and total basis set sizes
!
    if (jrot==0) kmin1= 1
    if (zncor) then
        jk1= 1
        jrot=abs(jrot)
    else
        jk1= jrot + kmin1
    endif
    if (idia > -2) then
        nrade=npnt1*npnt2
        nrado=nrade
        mbass=nrade*jk1*lmax
    else
        nrade=npnt1*(npnt1+1)/2
        nrado=npnt1*(npnt1-1)/2
        jt=jk1/2
        mbass=(nrade+nrado)*jt
        if (2*jt /= jk1) then
           if (ipar == 0) mbass=mbass+nrade
           if (ipar /= 0) mbass=mbass+nrado
        endif
        mbass=mbass*lmax
    endif


!     reset number of vib-rot functions to be considered

!     write out data from triatom/rotlev runs

    write(6,300) mbass,nmax11,nmax12,lmax,jrot,ipar,kmin1,neval,nv1
300   format(5x,'total number of basis functions   = ',i5,/&
            5x,'number of r1 radial dvr points    = ',i5,/&
            5x,'number of r2 radial dvr points    = ',i5,/&
            5x,'number of angular functions       = ',i5,/&
            5x,'total angular momentum            = ',i5,/&
            5x,'ipar, vibrational parity          = ',i5,/&
            5x,'kmin= 1-p, where j + p is parity, = ',i5,/&
            5x,'total vib-rot functions available = ',i5,/&
            5x,'vib-rot functions taken (0 = all) = ',i5,/)

!     space for angular integration

    npot = lpot
    if (jrot > 0) then
        ltop = lpot
        if (idia == 2 .and. mod(ltop,2) /= ipar) ltop=ltop+1
        npot=((ltop+2)/2)*2
        nn2= npot/2
    endif

    return
end subroutine insize
   
      
!**************************************************004
subroutine core
!
    use xpect3_logic
    use xpect3_dim
    implicit double precision (a-h,o-y), logical (z)

    integer, dimension (jk1) :: nbass1,lmin1

!     generate the subindex arrays needed for trans

    call genind(nbass1,lmin1,mbass,jk1,nbmax1,iket)
    call x3main(nbass1,lmin1)
         
    stop
end subroutine core

!**************************************************005
    subroutine genind(nbass,lmin,mbass,jk,nbmax,ivec)
!
!     this subroutine sets up the sub-index arrays needed for trans
!     and reads in the basis function labels for the vectors.
!     it calculates nbmax, the largest value of nbass, neeeded to
!     dimension the space needed for the d-coefficients.
!
    use xpect3_mass
    use xpect3_logic
    implicit double precision(a-h,o-y),logical (z)

    integer, dimension (jk) :: nbass,lmin,lbass

    if (zprint) write(6,"(//,'   j + kmin =',i3,'   mbass=',i6,/)") jk,mbass

!     read in the basis function labels
    read(ivec) mbass0,lmin,lbass,nbass
    if (mbass0>mbass) goto 999

!     generate the sub-index arrays and  find nbmax

    nbmax=nbass(1)
    do 2 k=2,jk
    nbmax=max(nbmax,nbass(k))
2     continue
    if (zprint) then
        write(6,"(//,5x,'indices generated by genind',/)")
        write(6,"(5x,'nbmax= ',i5,/,5x,'nbass, lmin and lbass follow',/)") nbmax
        write(6,*) (nbass(k), k=1,jk),(lmin(k), k=1,jk),&
                (lbass(k), k=1,jk)
    endif
    return
999   write(6,"(//,5x,a100,i3,a20,i5,a20,i5,a20)") 'basis function dimensions in error unit =',ivec,&
            &' mbass =',mbass,' expected, =',mbass0,' found'
    stop
    end

!**************************************************006
    subroutine x3main(nbass1,lmin1)

!     in this part, all data etc for the ket are labelled 1;
!     all data etc for the bra are labelled 2.
!
    use xpect3_timing
    use xpect3_dim
    use xpect3_mass
    use xpect3_logic
    implicit double precision (a-h,o-y), logical (z)
    parameter (mxprop=1000)
    common /pot/ iprop(mxprop)

    integer, dimension (jk1) :: nbass1,lmin1
    double precision, dimension(neval) :: eh
    double precision, dimension(npot) :: xd,wtd,xalf
    double precision, dimension(npropin) :: prin,trin
!     arrays for lagpt
    double precision, allocatable :: d0(:,:)
    double precision, dimension(npnt1) :: r1 
    double precision, dimension(npnt2) :: r2
!     arrays for trans
    double precision, dimension(nv1,npropin) :: props
    double precision, allocatable :: dc1(:)


!     array for dsrd
    integer, dimension (max(npot,lmax)) :: iv
    data x0/0.0d0/,xp5/0.5d0/

    allocate(d0(nrade*npot,npropin), dc1(neval*nbmax1/lmax*npot))

!***************************************************************************
!     read in radial dvr grid points 

    call getrow(r1,npnt1,iket)
    if (idia > -2)  call getrow(r2,npnt2,iket)
    if (jk1 <= 1) then !jrot==0) then
        !if (idia == -1) then
        read(iket) xalf
        read(iket) 
        read(iket)
        if (idia > -1) then
            read(iket) xalf
            do i=1,lmax
               xalf(2*lmax+1-i) = -xalf(i)
            end do
        else if (idia == -2) then
            read(iket)
            read(iket)
        endif
    endif

!     read in of energies for the ket

    read(iket) neval
    read(iket) eh

! J = 0
    if (ipar==0.and.jrot==0.and.zprint) then
        write(6,*) 'number of eigenvalues: ', neval
        write(6,"(//5x,'energies for the ket in cm-1',/)")
        write(6,"(5d24.12,/)") (eh(i)*219474.63,i=1,neval-1)

    end if

!     zero properties array
    props = x0

    if (jrot==0) then
!     J=0 case
        call dsrd(dc1,iket,mbass,nbass1(1),jk1,iv,ibass)
         
        if (idia==-2) then
            iadd=ipar
        else
            i0=npnt1
        endif
        ii=0
        do i2=1,npnt2
            if (idia==-2) then
               rr2=r1(i2)
               i0=i2-iadd
            else
               rr2=r2(i2)
            endif
            do i1=1,i0

! ######################################################################
! ###### IDIA = -2 RADAU HOMONUCLEAR ###################################
! ######################################################################
    if (idia==-2) then
            do 20 j=1,npot
                     ! Calling prop for each r1,r2,and alpha point, except 
                     ! symmetrically equilavent points.
            if(iv(j)==0) go to 20
            ii=ii+1
                     
            call prop(prin,r1(i1),rr2,xalf(j),npropin,nprt )
                          
            call prop(trin,rr2,r1(i1),xalf(j),npropin,nprt )
                     
            do n=1,npropin
                d0(ii,n)=xp5*(prin(n)+trin(n))
            end do
20       continue

! ######################################################################
! ###### IDIA = -1 RADAU HETERONUCLEAR #################################
! ###### IDIA = 1 JACOBI HETERONUCLEAR #################################
! ###### IDIA = 2 JACOBI HOMONUCLEAR   #################################
! ######################################################################
    else
        do 21 j=1,lmax
            if(iv(j)==0) go to 21
                ii=ii+1
                     
                call prop(prin,r1(i1),rr2,xalf(j),npropin,nprt)   
                if (idia == 2) then
                    call prop(trin,r1(i1),rr2,-xalf(j),npropin,nprt)
                    do n=1,npropin
                        d0(ii,n)=xp5*(prin(n)+trin(n))
                    end do
                else
                    do n=1,npropin
                        d0(ii,n)=prin(n)
                    end do
                endif
21                  continue 
            end if
        end do
        end do
                
                !write(32,*)'*'
                !do n=1,npropin
                !   write(32,*)(ii,n)
                !enddo


!       call to trans
    write(6,"(/,5x,'trans called after')")
    call system_clock(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,' secs cpu time used'/)") itime
    call trans(props,d0,dc1,ibass)
            
    else
!       J > 0 case
!       first generate angular grid and then call props

! Kolya
!  print *,'x3main: '
!  print *, d0(1,1),r1,r2,xd,wtd

    call lagpt(d0,r1,r2,xd,wtd)

    do 10 k1= 1,jk1
    if(nbass1(k1)==0) goto 10
    ndbass = nbass1(k1)*npot/lmax
    ibass = nbass1(k1)*neval
    kz = k1 -  kmin1

    call dsrd(dc1,iket,mbass,ndbass,jk1,iv,ibass)
    call rest(dc1,ndbass,jk1,xd,wtd,kz,lmin1)

!     call to trans
    write(6,"(/,5x,'trans called after')")
    call system_clock(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,' secs cpu time used'/)") itime
    ibasst = ibass/neval
    call trans(props,d0,dc1,ndbass)
10      continue
    end if

    deallocate(dc1)

!     end of property calculation

    write(6,"(/,5x,'outp called after')")
    call system_clock(itime2,irate2,imax2)
    itime=(itime2-itime0)/irate2
    write(6,"(/i10,' secs cpu time used'/)") itime

    call outpt2(props,eh)

    return

    stop
    end

!**************************************************011
    subroutine lagpt(d0,r1,r2,xd,wtd)

!     subroutine lagpt obtains values of the dipole at the radial
!     dvr points and angular integration points

    implicit double precision(a-h,o-y), logical (z)
    common/dim/ ncoord,npnt,npnt1,npnt2,nrade,nrado,&
                lpot,npot,nbmax1,mbass,kmin1,jk1,neval,&
                nn2,lmax,npropin,nprt,jrot,idia,ipar
 
    double precision, dimension(nrade*npot,npropin) :: d0
    double precision, dimension(npnt1) :: r1
    double precision, dimension(npnt2) :: r2
    double precision, dimension(npot) :: xd,wtd,b,c
    double precision, dimension(npropin) :: prin,trin
 
    data x0/0.0d0/,toler/1.0d-8/,xp5/0.5d0/,&
        x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/

!     set up points & weights for npot point angular integration

    xnu= x0
    alf= xnu
    bta= alf
    do i=2,npot
        xi= dble(i)
        b(i)= (alf+bta)*(bta-alf)/&
            ((alf+bta+x2*xi)*(alf+bta+x2*xi-x2))
        c(i)= x4*(xi-x1)*(alf+xi-x1)*(bta+xi-x1)*(alf+bta+xi-x1)/&
            ((alf+bta+x2*xi-x1)*(alf+bta+x2*xi-x2)*&
             (alf+bta+x2*xi-x2)*(alf+bta+x2*xi-x3))
    end do
  
    call jacobi(npot,nn2,xd,wtd,alf,bta,b,c,csa,tsa)

    write(6,1000) npot,0,(xd(i),wtd(i),i=1,nn2)
 1000 format(//,i8,' point gauss-associated legendre integration',&
            ' with k =',i2,&
            //,5x,'integration points',11x,'weights',&
            /,(f23.15,d25.12))
    write(6,1010) csa,tsa
 1010 format(/,5x,'computed sum of weights',d26.15,&
            /,5x,'exact    sum of weights',d26.15//)
    if (abs((csa-tsa)/tsa) > toler) goto 930
!     define other integration points
    do i=1,nn2
        j=i+nn2
        xd(j)=-xd(nn2+1-i)
        wtd(j)=wtd(nn2+1-i)
    end do
    do i = 1, 2*nn2
        xd(i) = -xd(i)
    end do

!     calculate properties at (r1,r2,cos\theta)

    i0=npnt1
      
    icall = 0
    ii = 0
    do i2=1,npnt2
        if (idia == -2) then
           rr2=r1(i2)
           i0=i2
        else
           rr2=r2(i2)
        end if
        do i1=1,i0
          do j=1,npot
            ii=ii+1
            icall = icall + 1
            call prop(prin,r1(i1),rr2,xd(j),npropin,nprt )
            if (idia == 2) then
               call prop(trin,r1(i1),rr2,-xd(j),npropin,nprt )
               do n=1,npropin
                 d0(ii,n)=xp5*(prin(n)+trin(n))
               end do          
            else if (idia == -2) then
               call prop(trin,rr2,r1(i1),xd(j),npropin,nprt )
               do n=1,npropin
                 d0(ii,n)=xp5*(prin(n)+trin(n))
               end do          
            else
               do n=1,npropin
                 d0(ii,n)=prin(n)
               end do
            endif
          end do
        end do
    end do

    return
930 write(6,"(/,5x,'gauss-legendre weights in error: adjust algorithm')")
 ! 940 format(/,5x,'gauss-legendre weights in error: adjust algorithm')
    stop
    end

!===========================================================================
    subroutine jacobi(nn,nn2,x,a,alf,bta,b,c,csa,tsa)

!     calculates zeros x(i) of the nn'th order jacobi polynomial
!     pn(alf,bta) for the segment (-1,1) & and corresponding weights
!     for gauss-jacobi integration. this routine uses a brute force
!     search for zeros due to GJ Harris (2001).
!     note that for our purposes, alf= bta= nu.


    implicit double precision(a-h,o-z)
    double precision, dimension(nn) :: x,a,b,c,xt
    data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x6/6.0d0/,eps/1.0d-12/
    fn= dble(nn)
    csa= x0
    beta= x1
    if (alf==1) beta = beta/x6
    cc= beta*x2**(alf+bta+x1)
!      tsa= cc/x2
    tsa= cc
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
        if (xxx > -0.0002) goto 30
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

!***********************************************019
    subroutine root(x,nn,alf,bta,dpn,pn1,b,c,eps)

!     improves the approximate root x; in addition obtains
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x.

    implicit double precision(a-h,o-z)
    double precision, dimension(nn) :: b,c
    iter= 0
1    iter= iter + 1
    call recur(p,dp,pn1,x,nn,alf,bta,b,c)
    d = p/dp
    x = x - d
    if (abs(d) > eps .and. iter < 10) goto 1
    dpn= dp
    return
    end
      
      
!**************************************************020
    subroutine recur(pn,dpn,pn1,x,nn,alf,bta,b,c)
    implicit double precision(a-h,o-z)
    double precision, dimension(nn) :: b,c
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

!*******************************************************020
    subroutine dsrd(d,ivec,mmbass,nbass,jk,iv,ibass)

!     subroutine to read d coefficients from dstore data
!     if read directly from dvr3d (i.e. j <= 0), then transform to an
!     fbr in theta
    use xpect3_logic
    use xpect3_mass
    use xpect3_dim
    implicit double precision (a-h,o-y), logical(z)
    double precision, dimension(neval,nbass) :: d
    double precision, allocatable :: temp(:),plegd(:)
    integer, dimension(max(npot,lmax))  :: iv
      
    if (jk>1) then
        nread = nbass/npot*lmax
        d = 0.0d0
        read(ivec) ((d(i,j),j=1,nread),i=1,neval)
!        orthonormality test
	sumd11 = 0.0d0
	do i=1,nread
        sumd11 = sumd11 + d(1,i)*d(1,i)
	end do
         
	write(6,*) 'dsrd: nread=', nread,'sumd11=', sumd11
         
    else
        rewind ivec
        do 30 i=1,6
            read(ivec)
30      continue
        if (.not. zbisc) read(ivec)
        if(idia>-2) then
           rewind ivec
           do i=1,7
             read(ivec)
           end do
           read(ivec) kz,maxleg,nidvr,lincr
        else
           read(ivec) kz,maxleg,nidvr
        end if
   
        leng=(maxleg+1)*nidvr
         ! New allocation below.
        allocate(plegd(leng))
        plegd = 0.0d0

        call getrow(plegd,leng,ivec)
         !call getrow(plegd,nbass,ivec)
         
        if (zbisc) then            
            read(ivec) iang,ibass
            read(ivec) (iv(ii),ii=1,nidvr)
        else
            iang=nidvr
            ibass=mmbass
            do i=1,nidvr
              iv(i)=1
            end do
        end if
        read(ivec)   
        read(ivec)
         
        allocate(temp(ibass))
        temp = 0.0d0
         
        do l=1, neval
            call getrow(temp,ibass,ivec)
            do i=1, ibass
               d(l,i)=temp(i)
            end do
        end do

    end if

    deallocate(plegd,temp)

    return
    end

!******************************************************020a
    subroutine rest(d,nbass,jk,xd,wtd,kz,lmin1)

!     rest of dsrd
!     subroutine to read d coefficients from dstore data
!     if read directly from dvr3d (i.e. j <= 0), then transform to an
!     fbr in theta
    use xpect3_logic
    use xpect3_mass
    use xpect3_dim
    implicit double precision (a-h,o-y), logical(z)

    integer, dimension (jk) :: lmin1
    double precision, dimension(neval,nbass) :: d,temp
    double precision, dimension(npot,npot) :: plegd
    double precision, dimension(npot) :: xd,wtd

    data x0/0.0d0/,x1/1.0d0/
     
    temp = d
    d = x0 

    nrad = nbass/npot
    jdia=max(1,idia)
    nang = lmax*jdia
    nang2 = lmax
    npot2 = npot/jdia
!     calculate the legendre function
    if (idia>-2.and.jk==0) kz=0
    call asleg(plegd,d,nang-1,xd,npot,kz)
    if (idia==2) then 
        do i=1,npot
          scale = sqrt(wtd(i))
          jj = lmin1(kz+1) - kz + 1
          do j=1,nang2
            plegd(i,j) =  scale*plegd(i,jj)
            jj=jj+jdia
          end do
        end do
    else
        do i=1,npot
          scale = sqrt(wtd(i))
          do j=1,npot
            plegd(i,j) = scale*plegd(i,j)
          end do
        end do
    end if
!        orthonormality test
    do k1=1,nang2
      do k2=1,nang2
        sum = x0
         do j=1,npot
           sum = sum + plegd(j,k1)*plegd(j,k2)
         end do
        write(6,*) 'rest: k1,k2 =',k1,k2,'sum plegd  = ', sum
      end do
    end do


!     evaluate wavefunction at angular integration points

    beta=x1
    d = x0
    ipt=1
    jpt=1

    do  mn=1,nrad
        call dgemm('n','t',neval,npot,nang2,beta,temp(1,jpt),neval,plegd,&
            npot,beta,d(1,ipt),neval)
        ipt = ipt + npot
        jpt = jpt + nang2
    end do 
  
    return
    end

!******************************************************021
    subroutine asleg(pleg,pnorm,lmax,x,nn2,m)

!     calculate polynomials 1 to lmax at x = cos(theta) for m = 0 or 1,
!     using the routine of press et al, numerical recipes, p. 182,
!     for the polynomial part of associated legendre functions.
!     a factor of sin(theta)**m has NOT been removed from all functions.

    implicit double precision (a-h,o-z)

    double precision, dimension(nn2,0:lmax) :: pleg
    double precision, dimension(nn2) :: x
    double precision, dimension(0:lmax) :: pnorm

    data x1/1.0d0/,x2/2.0d0/

    if (m < 0) goto 999
    do 10 i=1,nn2
    if (abs(x(i)) > x1) goto 999
    pmm = x1
    fact = x1
    somx2=sqrt((x1-x(i))*(x1+x(i)))
    do j=1,m
        pmm = -pmm * fact * somx2
        fact = fact + x2
    end do
    pleg(i,0)= pmm
    pmmp1= x(i)*dble(m+m+1)*pmm
    pleg(i,1)= pmmp1
    ll=1

    do l= 2+m,lmax+m
        r2lm1 = dble(l+l-1)
        rlpmm1= dble(l+m-1)
        rlmm  = dble(l-m)
        pll= (x(i)*r2lm1*pmmp1 - rlpmm1*pmm)/rlmm
        pmm= pmmp1
        pmmp1= pll
        ll=ll+1
        pleg(i,ll)= pll
    end do
10    continue

!     set up the normalisation constants
!     (pnorm)**2 = (2j + 1)/2   *   (j - k)! / (j + k)!
    jj = -1
    do 13 j = m,lmax+m
    fact = x1
      do i = j-m+1,j+m
        facti = dble(i)
        fact = fact * facti
      end do
    jj = jj + 1
    pnorm(jj) = sqrt(dble(j+j+1) / (fact + fact))
   13 continue
!     now normalise the polynomials
    do jj=0,lmax
        do i=1,nn2
          pleg(i,jj) = pleg(i,jj) * pnorm(jj)
        end do
    end do
    return
999   write(6,"(//5x,'improper argument in subroutine asleg'/)")
    stop
    end



!*********************************************************022

    subroutine trans(t,d0,dc1,ndbass)

!     subroutine trans is the main working routine in program xpect3
!     it carries out the necessary angular integrations using 3-j
!     symbols as defined by brink and satchler "angular momentum".
!     trans is called by x3main for each k-k' overlap integral.

!     it uses a fast matrix multiplication routine mxmb.
    use xpect3_logic
    use xpect3_mass
    use xpect3_dim
    implicit double precision (a-h,o-y), logical (z)
    double precision, dimension(nv1,npropin) :: t   

    double precision, dimension(npot*nrade,npropin) :: d0
    double precision, dimension(neval,ndbass) :: dc1
    double precision, dimension(neval) :: sumd
!.......................................................................


    nrad = ndbass/npot
    if (zprint) write(6,"('  trans: ndbass, nrad = ',2i6)") ndbass, nrad


!    take square of the wavefunction   
    dc1=dc1**2

!    orthonormality test
    sumdc111 = 0.0d0
    do i=1,ndbass
       sumdc111 = sumdc111 + dc1(1,i)
    end do

    sumd = 0.0d0

    do j=1,neval
        do i=1,ndbass
           sumd(j) = sumd(j) + dc1(j,i)
        enddo
    enddo
         
!     write(6,*) 'trans: sum dc1(1,i) 1-1,',ndbass,'=',sumdc111
!     do j=1,neval
!        write(6,*) 'trans: sum dc1(',j,'i) 1-1,',ndbass,'=',sumd(j)
!     enddo

    if (jrot==0) then
        do i=1, nv1!neval
          do n=1, npropin
            do j=1, ndbass
              t(i,n) = t(i,n) + d0(j,n)*dc1(i,j)
            end do
          end do
        end do 
    else

        n0 = npnt1
        i0 = 1
        if(nrad==nrade) i0 = 0
        do i=1, nv1!neval
          do n=1, npropin
            ic = 0
            ip = 0
            do n2 = 1,npnt2
              if(idia==-2) n0 = n2 - i0
              do n1 = 1,n0
                do j = 1,npot
                  ic = ic + 1
                  ip = ip + 1
                  t(i,n) = t(i,n) + d0(ip,n)*dc1(i,ic)
                end do
              end do
              if(i0==1) ip = ip + npot
            end do 
          end do
        end do

    end if

    return
    end subroutine trans


!**************************************************026a
    subroutine outpt2(props,eh)
!
!     outpt2 outputs the data at the end of a property calculation
!
    use xpect3_logic
    use xpect3_mass
    use xpect3_dim
    implicit double precision(a-h,o-y), logical(z)

    parameter (mxprop=1000)
    common /pot/ iprop(mxprop)
    double precision, dimension (neval) :: eh
    double precision, dimension (neval) :: ehtmp
    double precision, dimension (nv1,npropin) :: props,xprop

!      double precision, dimension (700,300000) :: props,xprop
    double precision, dimension (npropin) :: pground

    character(len=2) s1
!     autocm converts atomic units to wavenumbers
    data autocm/ 2.19474624d+05/,x0/0.0d0/,s1/'  '/

    if(zncor) then
        j1= jk1
    else
        j1= jk1-kmin1
    endif

    write(6,"(//5x,a20,i4,a20,i4,a20,i4,a20,i4,//)") 'j1=',jrot, '  kmin1=',kmin1,&
                & '  idia=',idia,'  ipar=', ipar


    ezero2=x0
    read(ilev,*,end=510,err=510)
    read(ilev,*,end=510,err=510)
    read(ilev,"(d20.12)",end=510,err=510) ezero2
510 continue   
    ezero = ezero2 * autocm
    write(6,"(5x,'ground zero in hartrees ',d16.8)") ezero2
    write(6,"(5x,'ground zero =',d16.8,' cm-1')") ezero

!
    if (ipar==0.and.kmin1==1) then
        s1 = 'a1' 
    else if (ipar==1.and.kmin1==0) then
        s1 = 'b2' 
    else if (ipar==1.and.kmin1==1) then
        s1 = 'a2' 
    else if (ipar==0.and.kmin1==0) then
        s1 = 'b1'
    end if

    if (ztra .and. zform) then
        open(unit=itra,form='formatted')
    else if (ztra) then
        open(unit=itra,form='unformatted')
    endif
    write(6,"(//,' J     ie       energy   property')")

    if (ipar==0.and.jrot==0.and.zfit) then
! Kolya
        do ie=2,nv1!neval
          do n=1,npropin
            xprop(ie,n) = props(ie,n) - props(1,n)
          end do
        end do
        do n=1,npropin
          write(itra0,*) props(1,n)
        end do       

        do 4 ie=1,nv1-1 !neval-1
        xeh=eh(ie+1)-ezero2
        xe1 = xeh*autocm
        if (ztra) write(itra,308) jrot,s1,ie,xe1,xeh
308     format(i2,1x,a2,i4,3x,f10.3,5(2x,d24.16)/&
            (17x,5(2x,d13.6)))
        do n=1, npropin
           write(12,*) xprop(ie+1,n)
        end do      
        write(6,308) jrot,s1,ie,xe1,(xprop(ie+1,n),n=1,npropin)
 4      continue

    else if (zfit) then

!   Read in ground state shifts
        if (zform) then
            open(unit=itra0,form='formatted')
            read(itra0,*) pground
        else
            open(unit=itra0,form='unformatted')
            read(itra0)
            read(itra0) pground
        endif
        close (unit=itra0)

        do ie=1,nv1 !neval
            do n=1,npropin
               xprop(ie,n) = props(ie,n) - pground(n)
            end do
        end do

        do ie=1,nv1 !neval
          xeh= eh(ie)-ezero2
          xe1 = xeh*autocm
           if (zform) then
              if (ztra) write(itra,308) jrot,s1,ie,xe1,xeh
              do n = 1,npropin
                if(ztra .and. xprop(ie,n)/=x0) write(itra,*) xprop(ie,n)
              end do
           else
              if (ztra) write(itra) jrot,s1,ie,xe1,xeh
              do n = 1,npropin
                if(ztra .and. xprop(ie,n)/=x0) write(itra) xprop(ie,n)
              end do
            endif
          write(6,308) jrot,s1,ie,xe1,(xprop(ie,n),n=1,npropin)
        end do

    else
        do ie=1,nv1 !neval
          xeh=eh(ie)-ezero2
          xe1 = xeh*autocm
          write(6,308) jrot,s1,ie,xe1,(props(ie,n),n=1,npropin)
        enddo

!     write out to stream itra?

        if (.not. ztra) return

          if (idia > 0) then
            isym=0
          else
            isym=abs(idia)
            if (ipar == 1) isym=-isym
          endif

        if (zform) then
!         is this the first write-out?
          if(.not.zstart) then
10          read(itra,*,end=90)
            goto 10
90         continue
! *****   inclusion of the following card is machine dependent *****
            backspace itra
           endif
!
        write(itra,"(i4,1x,i4,1x,i4,1x,i4)") jrot,nv1,ipar,kmin1!,neval!neval,ipar,kmin1
        write(itra,*) npropin,max(idia,0),isym,ezero
        write(itra,*) (iprop(n),n=1,npropin)
          
        ehtmp = eh*autocm

        call outro2(ehtmp,nv1,itra)

!          write(itra,*) eh*autocm
!701       format(i4,1x,i4,1x,i4,1x,i4)
          
        do n=1,npropin
            call outro2(props(1,n),nv1,itra)!neval,itra)
        end do

        else ! UNFORMATTED
!         is this the first write-out?
        if(.not.zstart) then
15        read(itra,end=95)
        goto 15
95        continue
! *****   inclusion of the following card is machine dependent *****
        backspace itra
        endif
!
    write(itra) jrot,neval,ipar,kmin1
    write(itra) npropin,max(idia,0),isym,ezero
    write(itra) (iprop(n),n=1,npropin)
    
    ehtmp = eh*autocm

    call outro2(ehtmp,nv1,itra)

    do n=1,npropin
        call outrow(props(1,n),nv1,itra)!neval,itra)
    end do
    endif
    endif
  
    close (unit=itra)
    return
    end
!                                                **027
    subroutine getrow(row,nrow,iunit)
!
    implicit double precision (a-h,o-y)
    double precision, dimension(nrow) :: row
    read(iunit) row
    return
    end
!                                                **028
    subroutine outrow(row,nrow,iunit)
!
    implicit double precision (a-h,o-y)
    double precision, dimension(nrow) :: row
    write(iunit) row
    return
    end
!                                                **028a
    subroutine outro2(row,nrow,iunit)
!
    implicit double precision (a-h,o-y)
    double precision, dimension(nrow) :: row
    write(iunit,*) row
    return
    end
!                                                **029
    subroutine dprint(dz0,dz1,nrad,lpot,jdia,npropin)

!     subroutine to print out the property integrals calculated
!     in subroutine lagpt for debugging purposes.
!
    use xpect3_mass
    implicit double precision (a-h,o-y), logical (z)
    double precision, dimension (nrad,npropin) :: dz0
    double precision, dimension (nrad,npropin,lpot) :: dz1
    write(6,"(//,5x,'p matrices calculated by subroutine lagpt',/)")
    do 1 n=1,npropin
    do 11 i=1,nrad
    write(6,*) dz0(i,n)
    write(6,*) (dz1(i,n,k),k=jdia,lpot,jdia)
11    continue
 1    continue
    return
    end

