FC=mpif90
FFLAGS=-ffree-form -ffree-line-length-none -fimplicit-none -O3 -pedantic-errors -Wall

all : gpse.x

gpse.x : gpse.o grid.o math.o psi.o
	$(FC) -o gpse.x $(FFLAGS) gpse.o grid.o math.o psi.o

gpse.o : gpse.f90 grid.mod math.mod psi.mod
	$(FC) -c $(FFLAGS) gpse.f90

psi.o psi.mod : psi.f90 grid.mod math.mod
	$(FC) -c $(FFLAGS) psi.f90

grid.o grid.mod : grid.f90
	$(FC) -c $(FFLAGS) grid.f90

math.o math.mod : math.f90
	$(FC) -c $(FFLAGS) math.f90

clean : 
	rm -f *.x *.o *.mod
