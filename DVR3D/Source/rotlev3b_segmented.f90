!     dummy main program                                           #001
      call rotlev3b
      stop
      end

!#######################################################################
      subroutine rotlev3b
      implicit double precision (a-h,o-y), logical (z)
      common /size/ nbass,mbass,ibass,neval,ipar,nmax,maxblk,jrot,&
                    kmin,kmax,meval,ndvr,iang,npnt,keval,nvib,mxblk2,neval2,&
                    nblk,loff,loff0,mbass0
      common /outp/ toler,thresh,zpham,zpvec,zvec,ztran,zptra,&
                    zpfun,ilev,ivec,ivec2,jvec,jvec2,kvec,kvec2,&
                    zdiag,zdcore,iscr,ires,irf1,irf2,zpseg
      namelist/prt/ toler,thresh,zpham,zpvec,zvec,ztran,zptra,&
                    zpfun,ilev,ivec,ivec2,jvec,jvec2,kvec,kvec2,&
                    zdiag,zdcore,iscr,ires,irf1,irf2,zpseg
      common/timing/itime0

     INTEGER :: count_0, count_rate, count_max, walltime , tstart, tend
      
      write(6,1000)
 1000 format('1',4x,'Program ROTLEV3B (version of March 2002):'/)
!     read in namelist input data (defaults in block data)
      read(5,prt)

      call nftim('beginning')
      
      call SYSTEM_CLOCK(itime0,irate2,imax2)
!     read in control parameters of problem.
      call insize

!     first for select, then onto main program
      call select 

      call SYSTEM_CLOCK(itime2,irate2,imax2) ! MP??? before was itime0
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,' secs CPU time used'/)

      stop
      end

!#######################################################################
      block data
      implicit double precision (a-h,o-y), logical (z)
      common /outp/ toler,thresh,zpham,zpvec,zvec,ztran,zptra,&
                    zpfun,ilev,ivec,ivec2,jvec,jvec2,kvec,kvec2,&
                    zdiag,zdcore,iscr,ires,irf1,irf2, zpseg
      data toler/0.0d0/,thresh/0.1d0/,zpham/.false./,zpvec/.false./,&
           ivec/26/,zvec/.false./,jvec/3/,jvec2/2/,iscr/1/,ires/0/,&
           ivec2/4/,zpfun/.false./,ilev/14/,kvec/8/,kvec2/9/,&
           zdiag/.true./,ztran/.false./,zptra/.false./,zdcore/.false./,&
           irf1/21/,irf2/22/,zpseg/.true./
      end
!########################################################################
