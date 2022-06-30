# DVR3D full manual

### Original project git repo: https://github.com/ExoMol/dvr3d

### Current project git repo: https://github.com/RunyuZhang-toma/DVR3D

### Dependency Library and :  

1. intel ifortran(standalone version): https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.4vpj2f
2. intel math kernal Library(standalone version): https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html
3. makefile
4. Ubuntu Linux LTS 20.04 x64
5. gcc
6. VMware

### Command Requirements

``` bash
source ~/.bashrc
source /opt/intel/oneapi/compiler/latest/env/vars.sh # default location
source /home/toma/intel/oneapi/setvars.sh # default location
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
! For subroutine / function / main function
!=========================================================
!Section ***
	subroutine name
	implicit
	common
	Data
	integer
	real
	allocate
	...
	end
!=========================================================
```

``` fortran
!Comment Section

!---------------------------------------------------------
!<Section number> <Section title>
!---------------------------------------------------------

!(same indent as code) <Comment>
```

### Current changed format

1. Regular indent 6.
2. Comment notation ! in the front line, Comment align to the code.
3. Every function, module, subroutine and main function have one blank line with before and after section.
4. if statement in one line mode should have a blank line after it.
5. The maximum line length permitted is 80 for old editor and terminals support.
6. Use >, >=, ==, <, <=, /= instead of .gt., .ge., .eq., .lt., .le., .ne. in logical comparisons. The new syntax, being closer to standard mathematical notation, should be clearer.
7. !(incomplete)Separate the information to be output from the formatting information on how to output it on I/O statements.  That is don't put text inside the brackets of the I/O statement.
