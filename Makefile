.SUFFIXES :
.SUFFIXES : .f90 .o

FC = gfortran
FFLAGS = -O2 -I/usr/local/include -Wall -g -fcheck=array-temps,bounds,do,mem,pointer,recursion -frecursive
LD = f90
LDFLAGS = -L/usr/local/lib
LDLIBS = -lfftw3 -lstdlib

SRC = ${wildcard *.f90} 
OBJ = ${SRC:.f90=.o}
TARGET=adv

$(TARGET) : $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LDLIBS) -o $@

analysis_module.o : grid_module.o
cascade_interpolation_module.o : search_module.o sort_module.o spline_interpolate_module.o
time_module.o: planet_module.o
lsqr.o : lsqrblas.o
glatwgt_module.o : math_module.o
legendre_transform_module.o : glatwgt_module.o alf_module.o fft_module.o
init_module.o : math_module.o planet_module.o sphere_module.o legendre_transform_module.o time_module.o
upstream_module.o : sphere_module.o grid_module.o time_module.o interpolate_module.o
upstream3d_module.o : sphere_module.o grid_module.o time_module.o interpolate3d_module.o
interpolate_module.o : math_module.o grid_module.o sphere_module.o bicubic_module.o polint_module.o
interpolate3d_module.o : math_module.o grid_module.o sphere_module.o tricubic_module.o
grid_module.o : legendre_transform_module.o init_module.o uv_module.o uv_hadley_module.o
euler_module.o : planet_module.o grid_module.o time_module.o legendre_transform_module.o uv_module.o
semilag_module.o : grid_module.o time_module.o legendre_transform_module.o upstream3d_module.o interpolate3d_module.o
forward_semilag_module.o : grid_module.o time_module.o legendre_transform_module.o upstream_forward_module.o interpolate3d_module.o
nisl_module.o : grid_module.o time_module.o legendre_transform_module.o upstream_module.o sphere_module.o interpolate1d_module.o interpolate3d_module.o upstream3d_module.o
direction_module.o : grid_module.o time_module.o legendre_transform_module.o upstream_module.o sphere_module.o interpolate1d_module.o interpolate3d_module.o upstream3d_module.o
main.o : grid_module.o time_module.o semilag_module.o analysis_module.o new_diagram_module.o nisl_module.o
new_diagram_module.o : grid_module.o time_module.o legendre_transform_module.o upstream_module.o sphere_module.o interpolate1d_module.o interpolate3d_module.o upstream3d_module.o interpolate_module.o

clean :
	rm -f *.o *.mod $(TARGET) *.dat $(TARGET).log *.txt

.f90.o :
	$(FC) $(FFLAGS) $< -c
