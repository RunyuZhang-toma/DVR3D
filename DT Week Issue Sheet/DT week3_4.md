# DT Week 3/4 issue sheet

1. Constants module complete.

* Some constant have no comments to indicate what do they mean.
* Some constants have no default value.
* Some subroutines defined itself module. Should we keep if or add to our constant module?

* In the source/spectra.f90 file the following two parts, there are two namelists but I cannot find the constants declaration.

* Could we use namelist directly as the the variable declaration or did I miss some code in other file?

``` fortran
! file spectra.f90
! 163
namelist/prt/ zout, zsort, zspe, zpfun, itra, ilev, ispe, item, &
wsmax,wsmin, emin, emax, jmax, smin, gz,zpseg

! 640
namelist /spe/ emin1,emax1,jmax,zplot,zemit,iplot,zfreq,zeinst, &
emin2,emax2,zprof,idat,zene,tinte,zlist,ilist, &
zdop,prthr,npoints,xmolm
    
```

2. Start to analysis the segmented code section.

3. Still waiting for the large test data.

4. Start to plan the final report with latex.

   1. Title page

   2. Abstract

   3. Contents

      * Introduction(current version)

      * Problems and shortages

      * Expectations

      * Environments

      * Project file structure

      * Version control

      * What I have done

      * Test and debug

      * Performance

      * Evaluation

   4. Table of contents

   5. illustrative material