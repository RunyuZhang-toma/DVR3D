      program shallot
      call dipole3b
      stop
      end

      subroutine dipole3b
      implicit double precision (a-h,o-y), logical (z)
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,&
      zmors2,zrme1,zrme2,zrme3
      namelist/prt/ zprint, zpmin, ztra, zstart,zrme1,zrme2,zrme3,&
                    iket, ibra, itra, iscr, ires, nblock
      common /head/ title
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      common/timing/itime0
      character(len=8) title(9)

      write(6,200)
200   format(//,5x,'Program DIPOLE3 (version of January 2020):',/)
      call SYSTEM_CLOCK(itime0,irate2,imax2)
      read(5,prt)
      read(5,100) title
100   format(9a8)
      write(6,202) title
202   format(5x,9a8)
      call insize
      call main
      stop
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      block data
      implicit double precision (a-h,o-y), logical (z)
      common /logic/ zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2,zrme1,zrme2,zrme3,&
                     zpseg
      common /stream/ iket, ibra, itra, iscr, ires, mblock, nblock
      data zmors1/.true./,zprint/.false./,ztra/.true./,zrme1/.false./,zrme2/.false./&
           zrme3/.false./,zmors2/.true./, zpmin /.false./, ires/0/, nblock/1000/,&
           zstart/.false./, iket/11/, ibra/12/, itra/13/,iscr/24/,zpseg/.true./
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
