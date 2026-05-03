SOURCE	= fft99_modified.f90 parameters.f90 compute_pvorticity_tendency.f90  init.f90 \
legt_i.f90 save_output.f90 convert_vor2uv.f90 convert_uv2pvor.f90 fft_d.f90 \
j_index2.f90 main_sqg.f90 d_legpol.f90 fft_i.f90 legt_d.f90 numerical_diffusion.f90 inverse_pvorticity.f90 compute_pvorticity_dissipation.f90 \
model.f90 compute_balance_equation.f90 convert_psi2uv.f90

FC = gfortran
FFLAGS = -fdefault-real-8 -O3 
#FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -frecursive
#FFLAGS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

OBJECTS	=	$(SOURCE:.f90=.o)

.SUFFIXES:
.SUFFIXES:	.o .f90

.f90.o :
	$(FC)  $(FFLAGS) -c $*.f90

main:	$(OBJECTS)
	$(FC) -o main ${OBJECTS} 

clean:
	 \rm -f $(OBJECTS) *~

fort:
	 \rm -f  fort.* *~
