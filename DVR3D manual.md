# DVR3D full manual

### Original project git repo: https://github.com/ExoMol/dvr3d

### Current project git repo: https://github.com/RunyuZhang-toma/DVR3D



### Dependency Library:  

1. intel ifortran(standalone version): https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.4vpj2f
2. intel math kernal Library(standalone version): https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html
3. makefile
4. Ubuntu Linux LTS 20.04 x64
5. gcc
6. VMware



### Development Platform

* Intel 10980xe
* Asus R6EE motherboard
* Crucial 16G*8 ddr4 3600
* Samsung PM983zet 960G
* Microsoft windows11 professional edition



#### Virtual Machine:

* CPU 8 cores
* Ram 8 GB
* ROM 50 GB



### Commands Requirements

``` bash
source ~/.bashrc
source /opt/intel/oneapi/compiler/latest/env/vars.sh # default location
source /home/toma/intel/oneapi/setvars.sh # default location

# compile and run in the sample folder such as ~/DVR3D/HCN
make dvr.out # in the project location
sh run.sh # in the HCN file
ll -htr
tail -n20 result.HCN.
tail -n20 result.HCN.*
```



### Code format rules

``` fortran
!For program file head
!=========================================================
!Copyright(C), 2022, UCL
!FileName: filename!Author: author
!Version: version***
!Data: Last modify data
!Description: main 
!			 function
!             interface
!             dependency!Other:!Function List: 1
!History Version info:
!=========================================================
```

```fortran
!========================================================= !Copyright(C). 2022, University College London
!File name: model.f90
!Author: Runyu Zhang & Jonathan Tennyson
!Version: 1.1
!Data: 7th/july/2022
!Description: this Fortran90 file contains the module !template for the source folder common
! constants divided by the common group
!Dependency: Folder source !=========================================================
```

``` fortran
!Comment for module definition

!========================================================
!By Runyu Zhang & Tennyson Jonathan 1/Sep/2022
!Contains size, outp, oupb, timing, split1, split2, mass 
!Special notice:: The module name contains the file name, 
!it means cannot directly paste to other
!files and use. There are some difference
!between the constant numbers , types and
! default value. !=========================================================
```

### 

### Updates

* Relation Operators Update
* Declaration Standard
* Constants Documents
* Constants Modulization
* Code Format Standardization
* Read/Write update
* Comment
