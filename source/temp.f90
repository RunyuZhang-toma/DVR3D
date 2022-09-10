!model.f90
!====================================================================================================
!Copyright(C). 2022, University College London
!File name: model.f90
!Author: Runyu Zhang & Jonathan Tennyson
!Version: 1.1
!Data: 7th/Sep/2022
!Description: this Fortran90 file contains the module for the source folder common constants
!             divided by the common group
!Dependency: Folder source
!===================================================================================================
!include 'common.f90'

module dipj0dvrtemp
implicit none
integer :: idia
integer :: ipar
integer :: npnt1
integer :: npnt2
integer :: jrot
integer :: kmin
integer :: neval
integer :: irate2
integer :: imax2
integer :: itime2
integer :: itime
end module dipj0dvrtemp


module rddatatemp
implicit none
data ibra0,iket0,iwave0,ivc0,ibra1,iket1,iwave1,ivc1,ione,itwo &
     &    /51,   52,   26,    54,  61,   62,   63,    64,  71,  72/
data zembed,zsame/.true.,.true./
data zdone/.true./
data iptot/0/
end module rddatatemp


module inmaintemp
implicit none
!constant count: integer 13.
    integer :: j
    double precision :: phi
    double precision :: phibra
    double precision :: phiket
    double precision :: dipx
    double precision :: dipz
    double precision :: Tx
    double precision :: Tz
    double precision :: evals
    double precision :: r1
    double precision :: r2
    double precision :: theta
    double precision :: r21
end module inmaintemp

!===================================================================================================

module rdphitemp
implicit none
!constant count: integer 3, real 4, ligical 1.

    integer :: itime
    integer :: mr2
    integer :: neval

    double precision :: phi
    double precision :: evals
    double precision :: autocm
    data autocm/2.19474624D+05/
end module rdphitemp

!===================================================================================================

module oldphitemp
implicit none
!constant count: integer 3, real 4, ligical 1.
    integer :: itime
    integer :: neval
    double precision :: phi
    double precision :: evals
    double precision :: array
    data autocm/2.19474624D+05/
end module oldphitemp

!===================================================================================================

module getphitemp
implicit none
!constant count: integer 3, real 4, ligical 1.

    integer :: j
    integer :: iidum
    integer :: k
    integer :: ndim2d
    integer :: n2d
    integer :: n3d
    integer :: npta
    integer :: nptb
    integer :: nptc
    integer :: neval
    integer :: ibra
    integer :: iket
    integer :: iout
    integer :: ivc
    integer :: ivec1
    integer :: ivec2
    integer :: jbra
    integer :: nbra
    integer :: jket
    integer :: nket
    logical :: ztheta
    logical :: zr2r1
    double precision :: phi
    double precision :: c2d
    double precision :: c1d
    double precision :: c3d
    double precision :: cprod
    double precision :: evals
    double precision :: autocm
end module getphitemp

!===================================================================================================

module getmutemp
implicit none
!constant count: integer 13.
    integer :: mr2
    double precision :: dipx
    double precision :: dipz
    double precision :: evals
    double precision :: r1
    double precision :: r2
    double precision :: theta
end module getmutemp

!===================================================================================================

module dipcaltemp
implicit none
!constant count: integer 13.
    integer :: mr2
    double precision :: dipx
    double precision :: dipz
    double precision :: xcos
    double precision :: r1
    double precision :: r2
end module dipcaltemp

!===================================================================================================

module getTtemp
implicit none
!constant count: integer 13.
    integer :: ibtime
    integer :: iktime
    integer :: mr2
    double precision :: phibra
    double precision :: phiket
    double precision :: dipx
    double precision :: dipz
    double precision :: Tx
    double precision :: Tz
    double precision :: evals
end module getTtemp


!===================================================================================================

module diffmutemp
implicit none
!constant count: integer 13.
    integer :: maxq
    double precision :: dipx
    double precision :: start
    double precision :: theta
    double precision :: r1
end module diffmutemp

!===================================================================================================

module gettratemp
implicit none
!constant count: integer 13.
    integer :: maxq
    double precision :: transe
    double precision :: transo
    double precision :: q
    double precision :: wt
    double precision :: b
    double precision :: c
    double precision :: dnorme
    double precision :: dnormo
    double precision :: toler
end module gettratemp

!===================================================================================================

module basistemp
implicit none
!constant count: integer 13.
    integer :: idim
    integer :: maxfn
    integer :: npts
    double precision :: rbasis
    double precision :: q
    double precision :: dnorm
    double precision :: A
    double precision :: alpha
end module basistemp

!===================================================================================================

module glagpttemp
implicit none
!constant count: integer 13.
    integer :: npts
    double precision :: csx
    double precision :: alpha
    double precision :: q
    double precision :: wt
    double precision :: b
    double precision :: c
    double precision :: csa
    double precision :: tsx
    double precision :: eps
    data eps/1d-12/
end module glagpttemp

!===================================================================================================

module LGRECRtemp
implicit none
!constant count: integer 13.
    integer :: NN
    integer :: itmax
    double precision :: PN
    double precision :: DPN
    double precision :: PN1
    double precision :: X
    double precision :: B
    double precision :: C
    double precision :: ALF
    double precision :: EPS
end module LGRECRtemp

!===================================================================================================

module getbastemp
implicit none
!constant count: integer 13.
    integer :: maxq
    double precision :: basise
    double precision :: basiso
    double precision :: reo
    double precision :: q
    double precision :: wt
    double precision :: b
    double precision :: c
    double precision :: dnorme
    double precision :: dnormo
    double precision :: toler
end module getbastemp

!===================================================================================================

module domulttemp
implicit none
!constant count: integer 13.
    double precision :: basise
    double precision :: basiso
    double precision :: r1
    double precision :: r2
    double precision :: theta
    double precision :: crunch
    double precision :: transo
    double precision :: transe
    double precision :: dipx
end module domulttemp
   

!===================================================================================================

module diffTtemp
implicit none
!constant count: integer 13.
    double precision :: Tz
    double precision :: Tx
    double precision :: evals
    double precision :: phiket
    double precision :: phibra
    double precision :: dipx
    double precision :: AUTOCM
    double precision :: AUTODE
    double precision :: DETOSEC
end module diffTtemp

!===================================================================================================

module insizetemp
implicit none
!constant count: integer 13.
    double precision :: toler
    data toler/1.0d-3/
end module insizetemp

!===================================================================================================

module genindtemp
implicit none
!constant count: integer 1, real 7.
     integer :: nbass
     integer :: mbass
     integer :: jk
     integer :: nbmax
     integer :: ivec

end module genindtemp

!===================================================================================================

module setfactemp
implicit none
!constant count: real 3.
     double precision :: binom
     double precision :: x1
     integer :: nbin
     data x1/1.0d0/
end module setfactemp

!===================================================================================================

module lagpttemp
implicit none
     integer :: nu
     double precision :: x0
     double precision :: toler
     double precision :: x1
     double precision :: x2
     double precision :: x3
     double precision :: x4
    data x0/0.0d0/,toler/1.0d-8/,&
        x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/
end module lagpttemp

!===================================================================================================

module aslegtemp
implicit none
!constant count: integer 25
     integer :: lmax
     integer :: m
     integer :: ipot
     double precision :: x
     double precision :: x1
     double precision :: x2
     data x1/1.0d0/,x2/2.0d0/
end module aslegtemp

!===================================================================================================

module jacobitemp
!constant count: integer 25
     integer :: nn
     double precision :: x1
     double precision :: x2
     double precision :: alf
     double precision :: bta
     double precision :: csa
     double precision :: tsa
     double precision :: x0
     double precision :: x4
     double precision :: x6
     double precision :: x8
     double precision :: eps
     data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x3/3.0d0/,x4/4.0d0/,x6/6.0d0/,&
          x8/8.0d0/,eps/1.0d-12/
end module jacobitemp

!===================================================================================================

module recurtemp
!constant count: integer 25
     integer :: nn
     double precision :: x0
     double precision :: x1
     double precision :: x2
     double precision :: alf
     double precision :: bta
     double precision :: tsa
     double precision :: dpn
     double precision :: pn1
     double precision :: pn
     data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
end module recurtemp

!===================================================================================================

module dsrdtemp
!constant count: integer 25
     integer :: ivec
     integer :: mmbass
     integer :: nbass
     integer :: ne
     integer :: jk
     integer :: ipar
     integer :: ibase
     integer :: kk
     integer :: nu
     integer :: jay_ipar
     integer :: iz
     double precision :: x0
     double precision :: x1
     double precision :: x2
     double precision :: d
     double precision :: temp
     double precision :: kneed
     double precision :: kbeg
     double precision :: xd
     data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
end module dsrdtemp

!===================================================================================================

module jtrantemp
!constant count: integer 25
     integer :: nrad
     integer :: mvib
     integer :: maxleg
     integer :: idvr
     integer :: kz
     integer :: ipar
     integer :: ivec
     integer :: iv
     integer :: iang
     integer :: ibass
     integer :: ibase
     integer :: nu
     integer :: jay_ipar

     double precision :: x0
     double precision :: pleg
     double precision :: dvrvec
     double precision :: temp
     data x0/0.0d0/
end module jtrantemp

!===================================================================================================

module transtemp
!constant count: integer 25
     integer :: k2
     integer :: k1
     integer :: nu
     integer :: ipar
     integer :: NCPUS

     double precision :: x0
     double precision :: t
     double precision :: dipol
     double precision :: binom
     double precision :: dc1
     double precision :: dc2
     double precision :: xfac
     double precision :: order

     data x0/0.0d0/
end module transtemp

!===================================================================================================

module specttemp
!constant count: integer 25
     double precision :: detosec
     double precision :: autode
     double precision :: x0
     double precision :: autocm

     data autocm/2.19474624d+05/,x0/0.0d0/,&
           autode/2.5417662d0/,&
           detosec/3.136186d-07/
end module specttemp

!===================================================================================================

module DIPDtemp
!constant count: integer 25
     integer :: NU

     double precision :: DIPC
     double precision :: RME
     double precision :: R1
     double precision :: R2
     double precision :: XCOS
     double precision :: X1
     double precision :: X0
     double precision :: TINY
     double precision :: X2
     double precision :: PI

     DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/,PI/3.1415927D0/
end module DIPDtemp

!===================================================================================================

module datatemp
!constant count: integer 25
     data zmors1/.true./, zprint/.false./, ztra/.true./,&
           zmors2/.true./, zpmin /.false./, ires/0/, nblock/1000/,&
           zstart/.false./, iket/11/, ibra/12/, itra/13/, iscr/24/,&
           zuvvis/.false./,zpseg/.false./
end module datatemp

!===================================================================================================

module dattemp
!constant count: integer 25
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
end module dattemp

!===================================================================================================

module ccmaintemp
implicit none
     integer :: nu
     double precision :: x0
     double precision :: toler
     double precision :: x1
     double precision :: x2
     double precision :: x8
     double precision :: xp5
     double precision :: sqrt2
     data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,x8/8.0d0/,xp5/0.5d0/,&
           toler/1.0d-8/,sqrt2/1.4142135623731d0/
end module ccmaintemp

!===================================================================================================

module setcontemp
implicit none
     double precision :: x0
     double precision :: x1
     double precision :: x4
     double precision :: amtoau
     double precision :: xp5
     double precision :: fixcos
     data amtoau/1.8228883d03/
     data x0,xp5,x1,x4/0.0d0,0.5d0,1.0d0,4.0d0/
end module setcontemp

!===================================================================================================

module nfmaintemp
implicit none
     double precision :: x16
     double precision :: x8
     double precision :: x4
  
     data x4/4.0d0/,x8/8.0d0/,x16/1.6d1/
end module nfmaintemp

!===================================================================================================

module data1temp
implicit none
     data toler/0.0d0/,thresh/0.1d0/,zpham/.false./,zpvec/.false./,&
           ivec/26/,zvec/.false./,jvec/3/,jvec2/2/,iscr/1/,ires/0/,&
           ivec2/4/,zpfun/.false./,ilev/14/,kvec/8/,kvec2/9/,&
           zdiag/.true./,ztran/.false./,zptra/.false./,zdcore/.false./,&
           irf1/21/,irf2/22/
end module data1temp

!===================================================================================================


module data2temp
implicit none
      data zmors1/.true./, zprint/.false./,ztra/.true./,zrme1/.true./,zrme2/.true./&
           zrme3/.false./,zmors2/.true./, zpmin /.false./, ires/0/, nblock/1000/,&
           zstart/.false./, iket/11/, ibra/12/, itra/13/, iscr/24/,zpseg/.false./
end module data1temp

!===================================================================================================

module dmaintemp
implicit none
      double precision :: x0
      double precision :: x1
      double precision :: x2
      double precision :: dwl
      double precision :: hc
      double precision :: bk
      double precision :: autocm
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,dwl/0.0d0/
      data hc/ 1.9864476d-16/, &
           bk/ 1.3806581d-16/, &
           autocm/ 2.19474624d+05/
end module dmaintemp

!===================================================================================================

module threejtemp
implicit none
      double precision :: zero
      double precision :: one
      data zero,one/0.0d0,1.0d0/
end module threejtemp

!===================================================================================================

module gaujactemp
implicit none
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
      &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
end module gaujactemp

!===================================================================================================

module laguertemp
implicit none
      double precision :: eps
      double precision :: x1

      data eps/1.0d-12/,x1/1.0d0/
end module laguertemp

!===================================================================================================
module itmaxtemp
implicit none
      integer :: itmax

      data itmax/10/
end module itmaxtemp

!===================================================================================================

module dmaintemp
implicit none
      double precision :: x0
      double precision :: x1
      double precision :: x2
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/
end module dmaintemp

!===================================================================================================
module mkham1temp
implicit none
      double precision :: x0
      double precision :: x1
      double precision :: xp5
      data x0/0.0d0/,xp5/0.50d0/,x1/1.0d0/
end module mkham1temp

!===================================================================================================

module diag3dtemp
implicit none
      double precision :: autocm
      double precision :: x0
      
      data autocm/2.19474624d+05/
      data x0/0.0d0/
end module diag3dtemp

!===================================================================================================

module SOLRTtemp
implicit none
      double precision :: sqrt2
      double precision :: x0
      double precision :: x16
      data x2/2.0d0/,x16/16.0d0/,sqrt2/1.4142135623731d0/
end module SOLRTtemp

!===================================================================================================



module DGROTtemp
implicit none
      double precision :: AUTOCM
      double precision :: x0
      
      DATA AUTOCM/2.19474624D+05/,x0/0.0d0/
end module DGROTtemp

!===================================================================================================

module gaslegtemp
implicit none
      double precision :: eps
      double precision :: x0
      double precision :: x1
      double precision :: x2
      double precision :: x3
      double precision :: x4
      double precision :: x8
      double precision :: xstep
      
      data x0/0.0d0/,x1/1.0d0/,x2/2.0d0/,eps/1.0d-12/,&
           x3/3.0d0/,x4/4.0d0/,x8/8.0d0/,xstep/1.0d-6/
end module gaslegtemp

!===================================================================================================

module DGitertemp
implicit none
      double precision :: X1
      double precision :: X0
      double precision :: EMAX
      integer :: noffd
      DATA X0/0.0D0/,X1/1.0D0/,EMAX/1.0D50/,noffd/1/
end module DGitertemp

!===================================================================================================


module timertemp
implicit none
      double precision :: total
      double precision :: user
      double precision :: system
  
      data total/0.0/,user/0.0/,system/0.0/
end module timertemp

!===================================================================================================