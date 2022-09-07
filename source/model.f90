module sizes
implicit none
    integer lbra0,nbra0,lket0,nket0,lbra1,nbra1,lket1,nket1,ntheta,nr1,nr2,neval0,neval1
    
end module sizes

module diffs
implicit none
    integer nqe,nqo,nr21
    double precision alphae,betae,alphao,betao
    logical zsame
end module diffs


module logic
implicit none
    logical zembed,zdone,zmors1,znco1,znco2,zprint,zpmin,ztra,zstart,zmors2,zrme1,zrme2,zrme3,zpseg,zuvvis,zout, zsort, zspe, zpfun,zncor,zfit,zform

    integer iptot,idia,iket,itra,itra0,ilev
end module logic

module base
implicit none
    integer ibase1,ibase2
end module base

module stream
implicit none
    integer ibra0,ibra1,iket0,iket1,iwave0,iwave1,ivc0,ivc1,ione,itwo,iket, ibra, itra, iscr, ires, mblock, nblock

end module stream

module mass
implicit none
     dimension xmass(3),xmassr(3)
     double precision g1,g2,ezero
     logical zembed,zbisc

end module mass

module eqm
implicit none
     dimension ex(3),ez(3)
     double precision tmass
end module eqm

module old
implicit none
     logic ztheta,zr2r1, zthet1,zr2r11
     integer npta, nptb, nptc, max2d, max3d,npta1,nptb1,nptc1,max2d1,max3d1
end module old

module time
implicit none
     double precision ouser,osys,ototal
end module time

module head
implicit none
     double precision title
end module head

module timing
implicit none
     integer itime0
end module timing

module sym
implicit none
     integer idia,ipar1,ipar2,jrot1,jrot2
end module sym

module dim
implicit none
     
     integer neval,lpot,ncoord,npnt,npnt1,npnt2,nrade,nrado,npot,nbin,nbmax1,nbmax2,mbass1,mbass2,mbass,kmin1,kmin2,jk1,jk2,neval1,neval2,nn2,ibase1,ibase2,ipot,nn2,lmax,npropin,nprt,jrot,idia,ipar,nv1
end module dim

module outp
implicit none
     logic zpham,zplot,zprad,zpvec,zrot,zladd,zembed,zmors2,zs0,zx,zs1,zpmin,zvec,zquad2,zdiag,zlmat,zcut,zall,zlin,zp1d,zp2d,zr2r1,ztheta,ztran,zmors1,ztwod,zbisc,zperp,zpseg,zpfun,zembed,zptra,zdcore,z1da
   
     integer kvecpb,ivec1,idiag,idiag1,idiag2,iout1,iout2,iwave,ilev,ieigs1,ivecs1,ieigs2,ivecs2,ivint,iband,intvec,jscr,jvec,jvec2,kvec,kvec2,iscr,ires,irf1,irf2,ivec,ivec2,nploti,nplotf,ithre 
     double precision toler,thresh
end module outp

module size
implicit none
     integer npnt1,npnt2,nktot,kpar,iqpar,maxblk_odd,nalf,nmax1,nmax2,nmax,maxblk,maxleg,nalf2,idvr,npnt,nlim1,nlim2,neval,ncoord,jrot,kmin,idia,ipar,max2d,max3d,max2d2,max3d2,npnta,npntb,npntc,ndima,ndimb,ndimc,iq,idia,ipar,lmax,npnt1,npnt2,jrot,kmin,neval,jk,ifile,NBASS,MBASS,IBASS,NEVAL,IPAR,IDIA,nlim,jrot,KMIN,NEVAL2,MEVAL,KEVAL,NVIB,NBLK,LOFF,LOFF0,kbass,npnt1,npnt2,npntt,kmax,ndvr,iang,mxblk2,mbass0

     double precision emax1,emax2
     
end module size

module split1
implicit none
     double precision re1,diss1,we1,beta1,ur1,urr1,a1
     integer iu1
end module split1

module split2
implicit none
     double precision re2,diss2,we2,beta2,ur2,urr2,a2
     integer iu2
end module split2

module oupb
implicit none
     double precision  xp0,xp1,xp2
end module oupb

module pot
implicit none
     parameter (mxprop=1000)

     dimension  iprop(mxprop)
end module pot

module pb
implicit none
     dimension inda1(100),inda2(100),indb1(100),indb2(100),indk(100)
     integer  iqa,iqb,isa,isb,ipa,ipb,kmina,kminb,nka,nkb,nbassa,nbassb,nskipka,nskipkb,mevala,mevalb,ibassa,ibassb,nviba,nvibb
end module pb
