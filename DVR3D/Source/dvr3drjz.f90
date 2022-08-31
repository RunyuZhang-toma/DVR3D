      program DVR3DRJZ
      use dvrmod
      call dvr3d
      stop
      end

!######################################################################
      subroutine dvr3d
      implicit double precision (a-h,o-y), logical (z)
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,&
                    zpseg
      common /oupb/   xp0,xp1,xp2
      namelist/prt/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,&
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,&
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zperp,&
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,&
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,&
                    zpseg
      common/timing/itime0

      write(6,1000)
 1000 format(5x,'Program DVR3DRJZ (version of March 2002)')

!     read in namelist input data (defaults in block data)
      read(5,prt)

      call SYSTEM_CLOCK(itime0,irate2,imax2)

!     read in control parameters of problem.
      call insize

!     now do the calculation
      call ccmain

      call SYSTEM_CLOCK(itime2,irate2,imax2)
      itime=(itime2-itime0)/irate2
      write(6,1)itime
 1    format(/i10,' secs CPU time used'/)
      stop
      end

!##############################################################################
      block data 
      implicit logical (z)
      common /outp/ zpham,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,& 
                    zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,& 
                    zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,& 
                    idiag1,idiag2,iout1,iout2,iwave,zpfun,ilev,& 
                    ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,&
                    zpseg
      data zpham/.false./,zprad/.false./,zpvec/.false./,zrot/.true./,& 
           zladd/.true./,zembed/.true./,zmors2/.true./,& 
           zpmin/.false./,zvec/.false./,zquad2/.true./,zcut/.false./,& 
           zdiag/.true./,zlmat/.false./,zall/.false./,& 
           zp1d/.false./,zp2d/.false./,zr2r1/.true./,ztheta/.true./,& 
           zmors1/.true./,ztran/.false./,ztwod/.false./,zperp/.false./,& 
            zx/.false./,zs0/.false./,zs1/.false./,zpseg/.false./,& 
           ieigs1/7/,ivecs1/3/,ieigs2/2/,ivecs2/4/,ivint/17/,& 
           iband/15/,intvec/16/,idiag1/20/,idiag2/21/,iout1/24/,& 
           iout2/25/,iwave/26/,zlin/.false./,zpfun/.false./,ilev/14/
      end

!############################################################################
      

