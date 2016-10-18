# =========================================================================
# Makefile : GPSE
#
# TESTED
#
#    GNU Make 3.81
#    GNU Fortran (Homebrew gcc 5.3.0) 5.3.0
#
#    GNU Make 3.81
#    GNU Fortran (GCC) 4.4.7 20120313 (Red Hat 4.4.7-17)
#
# LAST TESTED
#
#    GNU Make 3.81
#    GNU Fortran (Homebrew gcc 5.3.0) 5.3.0
#    Monday, October 17th, 2016
#
# -------------------------------------------------------------------------
#
#    Specify the SHELL that will interpret this Makefile. This line avoids
#    issues on systems where the SHELL variable might be inherited from
#    the environment. 

SHELL := /bin/bash

#    Get hostname of the machine the source will be compiled and run on.

HOSTNAME := $(shell hostname)

#    Specify general user-defined compilation options.

MPIF90       := mpif90
COMPILER     := gfortran
DEBUG        := OFF
OPTIMIZATION := ON
OPENMP       := OFF

# Set compiler-specific options.

ifeq ($(COMPILER),gfortran)

   STANDARD_OPTIONS     := -fimplicit-none -fmodule-private \
                           -ffree-form -ffree-line-length-none -std=gnu
   INTEGER_OPTIONS      := # -fdefault-integer-8 ! need to come up with 
                           # better way to provide overall support for 
                           # INT64; perhaps rely on MPI-3 mpi_08 module?
   REAL_OPTIONS         := -fdefault-real-8
   OPTIMIZATION_OPTIONS := -O2 -mtune=native
   OPENMP_OPTIONS       := -fopenmp
   CHECK_OPTIONS        := -fcheck=all
   DEBUG_OPTIONS        := -ffpe-trap=invalid,overflow -fbacktrace \
                           -fdump-core -finit-real=nan
   WARNING_OPTIONS      := -Wall -fmax-errors=0 -Wno-array-temporaries \
                           -Warray-bounds -Wcharacter-truncation \
                           -Wline-truncation -Wconversion-extra \
                           -Wimplicit-interface -Wimplicit-procedure \
                           -Wunderflow -Wextra -Wuninitialized

endif

COMPILER_OPTIONS := $(STANDARD_OPTIONS) $(INTEGER_OPTIONS) \
                    $(REAL_OPTIONS)

ifeq ($(DEBUG),ON)

   COMPILER_OPTIONS += -O0 -C -g $(CHECK_OPTIONS) $(DEBUG_OPTIONS) \
                       $(WARNING_OPTIONS)

endif

ifeq ($(OPTIMIZATION),ON)

   COMPILER_OPTIONS += $(OPTIMIZATION_OPTIONS)

endif

ifeq ($(OPENMP),ON)

   COMPILER_OPTIONS += $(OPENMP_OPTIONS)

endif

SOURCE_DIR := source
BUILD_DIR  := build
TARGET     := gpse.x 
SOURCES    := $(shell find $(SOURCE_DIR) -type f -name *.f90)
SOURCES    := $(filter-out $(SOURCE_DIR)/cli.f90, $(SOURCES))
OBJECTS    := $(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.f90=.o))
SOURCES    := $(filter-out $(SOURCE_DIR)/gpse.f90, $(SOURCES))
MODULES    := $(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.f90=.mod))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(MPIF90) $(COMPILER_OPTIONS) -o $@ $^

$(BUILD_DIR)/gpse.o: $(SOURCE_DIR)/gpse.f90 $(MODULES)
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/gpse.f90 -o $(BUILD_DIR)/gpse.o

$(BUILD_DIR)/evua.o $(BUILD_DIR)/evua.mod: $(SOURCE_DIR)/evua.f90 $(BUILD_DIR)/math.mod
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/evua.f90 -o $(BUILD_DIR)/evua.o

$(BUILD_DIR)/grid.o $(BUILD_DIR)/grid.mod: $(SOURCE_DIR)/grid.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/grid.f90 -o $(BUILD_DIR)/grid.o

$(BUILD_DIR)/grk4.o $(BUILD_DIR)/grk4.mod: $(SOURCE_DIR)/grk4.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/grk4.f90 -o $(BUILD_DIR)/grk4.o

$(BUILD_DIR)/io.o $(BUILD_DIR)/io.mod: $(SOURCE_DIR)/io.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/io.f90 -o $(BUILD_DIR)/io.o

$(BUILD_DIR)/math.o $(BUILD_DIR)/math.mod: $(SOURCE_DIR)/math.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/math.f90 -o $(BUILD_DIR)/math.o

$(BUILD_DIR)/psi.o $(BUILD_DIR)/psi.mod: $(SOURCE_DIR)/psi.f90 $(BUILD_DIR)/math.mod
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/psi.f90 -o $(BUILD_DIR)/psi.o

$(BUILD_DIR)/rot.o $(BUILD_DIR)/rot.mod: $(SOURCE_DIR)/rot.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/rot.f90 -o $(BUILD_DIR)/rot.o

$(BUILD_DIR)/vex.o $(BUILD_DIR)/vex.mod: $(SOURCE_DIR)/vex.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIF90) -J$(BUILD_DIR) $(COMPILER_OPTIONS) -c $(SOURCE_DIR)/vex.f90 -o $(BUILD_DIR)/vex.o

.PHONY: clean
clean:
	@echo " Cleaning...";
	$(RM) -r $(BUILD_DIR) $(TARGET) *.output *.vtk

# =========================================================================
