# =========================================================================
# NAME
#
#     gpse Makefile
#
# SYNOPSIS
#
# DESCRIPTION
#
#     This Makefile ...
#
# OPTIONS
#
# SEE ALSO
#
# KNOWN BUGS
#
# ADDITIONAL NOTES
#
# HISTORY
#
# AUTHOR(S)
#
#     Marty Kandes
#     Computational Science Research Center
#     San Diego State University
#     5500 Campanile Drive
#     San Diego, California 92182
#
# COPYRIGHT
#     
#     Copyright (c) 2014 Martin Charles Kandes
#
# LAST UPDATED
#
#     Wednesday, October 29th, 2014
#
# -------------------------------------------------------------------------

#     Specify the SHELL that will interpret this Makefile. This line avoids
#     issues on systems where the SHELL variable might be inherited from
#     the environment. 

SHELL    := /bin/bash

#     Get hostname of the machine gpse will be compiled and run on.

HNAME    := $(shell hostname)

#     Provide directory paths required for compilation.

DIR_EXE  := ./
DIR_SRC  := ./source/
DIR_OBJ  := ./objects/
DIR_MOD  := ./modules/

#     Specify general user-defined compilation options.

MPIF90   := mpif90
COMPILER := gfortran
DEBUG    := OFF
OPTIMIZE := ON
OPENMP   := ON

#     Specify compiler-specific compilation options.
#
#     GNU Fortran (GCC) 4.8.2 ( ivymike )
#     GNU Fortran (GCC) 4.1.2 20080704 (Red Hat 4.1.2-48) ( dugong.sdsu.edu )

STD_GNU  := -fimplicit-none -ffree-form -ffree-line-length-none -fall-intrinsics
INT_GNU  :=
REAL_GNU := -fdefault-real-8
OPT_GNU  := -O3 -ffast-math -funroll-loops
OMP_GNU  := -fopenmp
AUTO_GNU :=
CHK_GNU  := -fcheck=all
DBG_GNU  := -fmodule-private -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan
WRN_GNU  := -Wall -fmax-errors=0 -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized

#     ifort (IFORT) 11.1 20100414 ( dugong.sdsu.edu )

STD_IFORT  := -stand f03 -implicitnone -free -m64 -static -static-intel -static-libgcc -mkl=parallel
INT_IFORT  := -integer-size 64
REAL_IFORT := -real-size 64
OPT_IFORT  := -archSSE4.1 -fast -fp-model fast=2 unroll-agressive -vec -opt-report 3 -opt-report-file=gpse.opt -vec-report5
OPM_IFORT  := -openmp -openmp-report2
AUTO_IFORT := -parallel -par-num-threads=4 -par-schedule=auto -par-report3
CHK_IFORT  := -check all
DBG_IFORT  := -debug all
WRN_IFORT  := -warn all

#     pgf90 10.9-0 64-bit target on x86-64 Linux -tp istanbul-64 ( dugong.sdsu.edu )

STD_PGI    := -Mfree -m64 -mcmodel=medium -Minfo -Mneginfo -Mpreprocess -Dpgf95
INT_PGI    := -i4
REAL_PGI   := -r8
OPT_PGI    := -fast#-Mipa=fast,inline -Munroll
OMP_PGI    := -mp=allcores
AUTO_PGI   := -Mconcur=allcores
CHK_PGI    := -C
DBG_PGI    :=
WRN_PGI    :=

#     Set compiler-specific options.

ifeq ($(COMPILER),gfortran)

   STD  := $(STD_GNU)
   INT  := $(INT_GNU)
   REAL := $(REAL_GNU)
   OPT  := $(OPT_GNU)
   OMP  := $(OMP_GNU)
   CHK  := $(CHK_GNU)
   DBG  := $(DBG_GNU)
   WRN  := $(WRN_GNU)

endif

ifeq ($(COMPILER),ifort)

   STD  := $(STD_IFORT)
   INT  := $(INT_IFORT)
   REAL := $(REAL_IFORT)
   OPT  := $(OPT_IFORT)
   OMP  := $(OMP_IFORT)
   AUTO := $(AUTO_IFORT)
   CHK  := $(CHK_IFORT)
   DBG  := $(DBG_IFORT)
   WRN  := $(WRN_IFORT)

endif

ifeq ($(COMPILER),pgf90)

   STD  := $(STD_PGI)
   INT  := $(INT_PGI)
   REAL := $(REAL_PGI)
   OPT  := $(OPT_PGI)
   OMP  := $(OMP_PGI)
   AUTO := $(AUTO_PGI)
   CHK  := $(CHK_PGI)
   DBG  := $(DBG_PGI)
   WRN  := $(WRN_PGI)

endif

FFLAGS := $(STD) $(INT) $(REAL)

ifeq ($(DEBUG),ON)

   FFLAGS += -O0 -C -g $(CHK) $(DBG) $(WRN)

endif

ifeq ($(OPTIMIZE),ON)

   FFLAGS += $(OPT)

endif

ifeq ($(OPENMP),ON)

   FFLAGS += $(OMP)

else

   FFLAGS += $(AUTO)

endif

all : gpse.x

gpse.x : gpse.o evua.o grid.o io.o math.o psi.o vex.o 
	$(MPIF90) $(FFLAGS) -o $@ $^

gpse.o : gpse.f90 evua.mod grid.mod io.mod math.mod psi.mod vex.mod 
	$(MPIF90) $(FFLAGS) -c $<

evua.o evua.mod : evua.f90 grid.mod math.mod
	$(MPIF90) $(FFLAGS) -c $<

grid.o grid.mod : grid.f90
	$(MPIF90) $(FFLAGS) -c $<

io.o io.mod : io.f90 grid.mod 
	$(MPIF90) $(FFLAGS) -c $<

math.o math.mod : math.f90
	$(MPIF90) $(FFLAGS) -c $<

psi.o psi.mod : psi.f90 grid.mod math.mod
	$(MPIF90) $(FFLAGS) -c $<

vex.o vex.mod : vex.f90
	$(MPIF90) $(FFLAGS) -c $<

.PHONY : printfflags
.NOTPARALLEL : printfflags
printfflags :
	@echo $(FFLAGS)

.PHONY : clean
clean : 
	rm -f *.x *.o *.mod

# =========================================================================
