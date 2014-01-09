FC=gfortran
FFLAGS=-ffree-form -ffree-line-length-none -fimplicit-none -O3 -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -pedantic-errors -Wall

all : gpse.x

gpse.x : gpse.o math.o
	$(FC) -o gpse.x $(FFLAGS) gpse.o math.o

gpse.o : gpse.f90 math.mod
	$(FC) -c $(FFLAGS) gpse.f90

math.o math.mod : math.f90
	$(FC) -c $(FFLAGS) math.f90

clean : 
	rm -f *.x *.o *.mod
