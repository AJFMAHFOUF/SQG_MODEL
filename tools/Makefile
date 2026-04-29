SOURCE	= compute_legendre_polynomials.f90 factrl.f90 gammln.f90 gauleg.f90 plgndr.f90

FC = gfortran
FFLAGS = -fdefault-real-8 

OBJECTS	=	$(SOURCE:.f90=.o)

.SUFFIXES:
.SUFFIXES:	.o .f90

.f90.o :
	$(FC)  $(FFLAGS) -c $*.f90

main:	$(OBJECTS)
	$(FC) -o main.exe ${OBJECTS} 

clean:
	 \rm -f $(OBJECTS) *~

fort:
	 \rm -f  fort.* *~
