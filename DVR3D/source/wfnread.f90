!==================================================================================================
!Module defintion
!By Runyu Zhang & Tennyson Jonathan 1/Sep/2022
!Contains: size
!Special notice:: The module name contains the file name, it means cannot directly paste to other 
!                 files and use. There are some difference between the constant numbers, types and 
!                 default value.
!==================================================================================================

module wfnread_size
    integer :: idia    ! 1 scattering coordinates heteronuclear diatomic
                         ! 2 scattering coordinates homonuclear diatomic
                         ! -1 radau  coordinates hetronuclear diatomic
                         ! -2 radau  coordinates homonuclear  diatomic
                         ! 0 radau   coordinates with the z axis perpendicular to the molecular plane.
    integer :: ipar   ! parity of basis - if idia=+/-2: ipar=0 for even & =1 for odd
    integer :: lmax
    integer :: npnt1    ! number of (gauss-laguerre) dvr points in r1
    integer :: npnt2    ! number of (gauss-laguerre) dvr points in r2
    integer :: jrot   ! total angular momentum of the molecule
    integer :: kmin   ! zrot=t, kmin=1 sym. rot. basis, =0 anti-sym.
                         ! kmin=2 loop over both sym & anti-sym (zbisc=t only)
                         ! zrot=f, kmin=fixed value of k
    integer :: neval   ! number of eigenvalues which have to actually be supplied as output
    integer :: jk
    integer :: ifile

end module wfnread_size

program reader
! This is a small utility to read the wavefunction files.
! At present it just echos the reads to standard output but 
! it can adapted for other functions such as plotting wavefunctions

    use wfnread_size

! file is assumed to be attached to unit 26
    ifile=26
    open(unit=ifile,form='unformatted')

! reader header record to determine file structure  
    read(ifile) idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval
    write(*,*) idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval
    If (jrot==0) kmin=1
    jk=jrot+kmin



    If (idia == -2 .and. jk > 1) then
        print *, "Radau, 8 or 9"
        call read_8or9_radau
    elseif (idia == -2) then
        print *, "Radau, 26"
        call read_26_radau
    elseif (jk > 1) then
        print *, "Jacobi, 8 or 9"
        call read_8or9_jacobi
    else
        print *, "Jacobi, 26"
        call read_26_jacobi
    endif

end program reader


!###################################################################
subroutine read_8or9_jacobi
    use wfnread_size
    implicit double precision (a-h,o-y), logical (z)

    integer, allocatable ::  nbass(:), lmin(:), lbass(:)
    double precision, allocatable:: r1(:), r2(:), e(:), temp(:,:)
    dimension xm(3)
    

    read(ifile) zembed,zmorse1,zmorse2,xm,g1,g2,zncor
    write(*,*) zembed,zmorse1,zmorse2,xm,g1,g2,zncor

    read(ifile) re2,diss2,we2,re2,diss2,we2
    write(*,*) re2,diss2,we2,re2,diss2,we2

    allocate( nbass(jk),lmin(jk),lbass(jk) )
    read(ifile) mbass0,lmin,lbass,nbass
    write(*,*) mbass0,lmin,lbass,nbass

    !read radial grid points
    allocate( r1(npnt1), r2(npnt2) )
    call getrow(r1,npnt1,ifile)
    write(*,*) r1
    write(*,*)
    if (idia > -2) then
        call getrow(r2,npnt2,ifile) !only if idia > -2
        write(*,*) r2
    endif

    !read in the energies, Hartrees
    read(ifile) neval
    write(*,*) neval 
    allocate( e(neval) )
    read(ifile) e
    write(*,*) e

    ne=neval

! the array temp consists of neval eigenvectors. each eigenvector
! consists of the values of the wavefunction as expressed on the DVR
! grid in the 3 coordinates, written as a 3d array. thus if the
! eigenvector is placed an array vec(NALF, NPNT1, NPNT2), NALF, NPNT1,
! NPNT2 will give the location on the DVR grid for that value of the
! wavefunction.

    print *, "nbass ", nbass(1)
    print *, "ne ", ne
    print *, "jk ", jk
    do jay1=0, jk-1
        print *, "k=", jay1+1-kmin, " block"
    allocate ( temp(ne,nbass(1)) )
    read(ifile) ((temp(i,j),j=1,nbass(1)),i=1,ne)
    
    !  read(ifile) temp
    write(*,*) temp
    deallocate( temp )
    enddo

end subroutine read_8or9_jacobi

!##################################################################

subroutine read_26_jacobi
    use wfnread_size
    implicit double precision (a-h,o-y), logical (z)

    integer :: mbass0,lmin,lbass,nbass
    !  double precision, allocatable:: r1(:), r2(:), e(:), d(:), plegd(:,:)
    double precision, allocatable:: e(:), d(:)
    dimension xm(3)

    read(ifile) zembed,zmorse1,zmorse2,xm,g1,g2,zncor
    write(*,*) zembed,zmorse1,zmorse2,xm,g1,g2,zncor

    read(ifile) re2,diss2,we2,re2,diss2,we2
    write(*,*) re2,diss2,we2,re2,diss2,we2

    do i=1, 6
        read(ifile)
    enddo

    !read in the energies, Hartrees  
    read(ifile) neval
    write(*,*) neval
    allocate( e(neval) )
    read(ifile) e
    write(*,*) e

    rewind ifile
    do i=1, 6
        read(ifile)
    enddo
    
    !  if (.not. zbisc)
    read(ifile)
 
    read(ifile) kz,maxleg,nidvr
    write(*,*) kz,maxleg,nidvr

    ! call getrow(plegd,leng,ifile)
    ! write(*,*) plegd
    read(ifile)
    read(ifile)
    read(ifile)

    iang=nidvr
    ibass=mbass0
    ne=neval
    print *, "ne ", ne
    print *, "ibass ", ibass
    print *, "nbass ", nbass
    print *, "iang ", iang
    

! D is consists of the values of the wavefunction as expressed on the
! the DVR grid in the 3 coordinates, writen as a 3d array d(NALF,
! NPNT1, NPNT2) NALF, NPNT1, NPNT2 will give the location on the DVR
! grid for that value of the wavefunction.

    allocate( d(iang*npnt1*npnt2) ) 
    do i=1, ne
        write(*,*) "eigenvector ", i
        call getrow(d,(iang*npnt1*npnt2),ifile)
    write(*,*) d
    enddo

end subroutine read_26_jacobi

!###################################################################

subroutine read_26_radau
    use wfnread_size
    implicit double precision (a-h,o-y), logical (z)

    integer, allocatable ::  nbass(:), lmin(:), lbass(:), iv(:)
    double precision, allocatable:: e(:), d(:)
    dimension xm(3)

    read(ifile) zembed,zmors1,zmors1,xmass,g1,g2,zncor,zquad2
    print *, zembed,zmors1,zmors1,xm,g1,g2,zncor,zquad2

    read(ifile) re1,diss1,we1,re1,diss1,we1
    print *, re1,diss1,we1,re1,diss1,we1

    do i=1,5
        read(ifile)
    enddo
    
    read(ifile) iang2,ibass2
    print *, "iang2,ibass2 ",  iang2,ibass2
    read(ifile)

    read(ifile) meval2
    print *, "meval2 ",  meval2

!get the enrgy levels
    allocate( e(meval2) )
    call getrow(e,meval2,ifile)
    print *, "energy levels in Hartrees"
    print *, e


! D is consists of the values of the wavefunction as expressed on the
! the DVR grid in the 3 coordinates, writen as a 3d array d(NALF,
! NPNT1, NPNT2) NALF, NPNT1, NPNT2 will give the location on the DVR
! grid for that value of the wavefunction.

!get the eigenvectors
    allocate( d(ibass2) )
    do i=1, meval2
        call getrow(d,ibass2,ifile)
        print *, "eigenvector ", i
        print *, d
    enddo

end subroutine read_26_radau

!###################################################################

subroutine read_8or9_radau
    use wfnread_size
    implicit double precision (a-h,o-y), logical (z)

    integer, allocatable ::  nbass(:), lmin(:), lbass(:)
    double precision, allocatable:: e(:), d(:)
    dimension xm(3)

    read(ifile) zembed,zmorse1,zmorse2,xm,g1,g2,zncor
    write(*,*) zembed,zmorse1,zmorse2,xm,g1,g2,zncor

    read(ifile) re2,diss2,we2,re2,diss2,we2
    write(*,*) re2,diss2,we2,re2,diss2,we2

    allocate( nbass(jk),lmin(jk),lbass(jk) )
    read(ifile) mbass0,lmin,lbass,nbass
    write(*,*) mbass0,lmin,lbass,nbass

    read(ifile)

    !read in the energies, Hartrees
    read(ifile) neval
    write(*,*) neval 
    allocate( e(neval) )
    read(ifile) e
    write(*,*) e

! the array d consists of neval eigenvectors. Thus d(1:nbass(k))
! represents the 1st eigenvector; d( nbass(k)+1:nbass(k) ) the second,
! etc. each eigenvector consists of the values of the wavefunction as
! expressed on the DVR grid in the 3 coordinates, written as a 3D
! array. thus if the eigenvector is placed an array vec(NALF, NPNT1,
! NPNT2), NALF, NPNT1, NPNT2 will give the location on the DVR grid for
! that value of the wavefunction.

  !read in wavefunctions, each k block is read in separtly
    do k=1, jk
        allocate( d(nbass(k)*neval) )
        print *, "K block ", k-1
        read(ifile) d
        print *, d
        deallocate( d )
    enddo
end subroutine read_8or9_radau

!#####################################################################

subroutine getrow(row,nrow,iunit)
    implicit double precision (a-h,o-y)
    double precision, dimension(nrow) :: row      
    read(iunit) row
    return
end subroutine getrow

