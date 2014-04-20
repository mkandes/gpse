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
#     Sunday, March 30th, 2014
#
# -------------------------------------------------------------------------

#     Specify the SHELL that interprets this Makefile. 

SHELL    = /bin/bash

#     Get machine architecture information

ARCH     = $(shell uname -m)

#     Provide directory and file location information.

DIR_EXE  = ./
DIR_SRC  = ./source/
DIR_OBJ  = ./objects/
DIR_MOD  = ./modules/

#     Specify general user-defined compilation options here.

COMPILER = gfortran
INTEGER  = DOUBLE
REAL     = DOUBLE
DEBUG    = NO
OPTIMIZE = NO
OPENMP   = ON
MPI      = OFF

#     Specify compiler-specific compilation options here.
#
#     gfortran:

STD_GNU  = -fimplicit-none -ffree-form -ffree-line-length-none -fall-intrinsics -fdefault-real-8
INT_GNU  =
RE_GNU   = 
DBG_GNU  = -fmodule-private -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan
CHK_GNU  = -fcheck=all
WRN_GNU  = -Wall -fmax-errors=0 -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized
OPT_GNU  = -O3
OMP_GNU  = -fopenmp

#     ifort: 

STD_INTC = 
INT_INTC = 
RE_INTC  = 
DBG_INTC = -debug all
CHK_INTC = -check all
WRN_INTC = -warn all
OPT_INTC = -O3
OMP_INTC = -openmp
MPI_INTC =

#     Set compiler-specific options.

ifeq ($(COMPILER),gfortran)

   FC    = gfortran
   STD   = $(STD_GNU)
   INT   = $(INT_GNU)
   RE    = $(RE_GNU)
   DBG   = $(DBG_GNU)
   CHK   = $(CHK_GNU)
   WRN   = $(WRN_GNU)
   OPT   = $(OPT_GNU)
   OMP   = $(OMP_GNU)

endif

FC = mpif90
FFLAGS = $(STD) $(OPT) $(OMP)

all : gpse.x

gpse.x : gpse.o grid.o io.o math.o psi.o vex.o lib_vtk_io.o lib_base64.o ir_precision.o
	$(FC) $(FFLAGS) -o $@ $^

gpse.o : gpse.f90 grid.mod io.mod math.mod psi.mod vex.mod lib_vtk_io.mod lib_base64.mod ir_precision.mod
	$(FC) $(FFLAGS) -c -cpp $<

grid.o grid.mod : grid.f90
	$(FC) $(FFLAGS) -c -cpp $<

io.o io.mod : io.f90 grid.mod lib_vtk_io.mod
	$(FC) $(FFLAGS) -c -cpp $<

math.o math.mod : math.f90
	$(FC) $(FFLAGS) -c -cpp $<

psi.o psi.mod : psi.f90 math.mod
	$(FC) $(FFLAGS) -c -cpp $<

vex.o vex.mod : vex.f90
	$(FC) $(FFLAGS) -c -cpp $<

lib_vtk_io.o lib_vtk_io.mod : lib_vtk_io.f90 lib_base64.mod ir_precision.mod
	$(FC) $(FFLAGS) -c -cpp $<

lib_base64.o lib_base64.mod : lib_base64.f90 ir_precision.mod
	$(FC) $(FFLAGS) -c -cpp $<

ir_precision.o ir_precision.mod : ir_precision.f90
	$(FC) $(FFLAGS) -c -cpp $<

clean : 
	rm -f *.x *.o *.mod

# =========================================================================
