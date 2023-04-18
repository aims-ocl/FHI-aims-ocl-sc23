# for ifort
 FC = ifort
FFLAGS = -O3
#  FFLAGS = -O2 -warn unused -warn interfaces -warn all   -g -fpstkchk -check bounds -check uninit -traceback -fpe
 F90FLAGS = $(FFLAGS)
 ARCH = arch_generic.f90

# p4, p5
#FC = xlf90
#FFLAGS = -qfixed -O3
#F90FLAGS = -O3
#ARCH = arch_alpha.f90

# compaq
#FC = f95
#FFLAGS = -O5 -check nopower
#F90FLAGS = $(FFLAGS)
#ARCH = arch_alpha.f90

MODULES = vector.o \
          vector_array.o \
          cluster.o \
          geometry.o \
          geometry_list.o \
          control.o \
	  basin_hopping.o \
	  dft_file.o \
	  conjugate_gradient.o \
          lennard_jones.o \
          map_to_one_dimension.o \
          spline.o \
	  pair_pot.o \
	  cartesian_ylm.o

SUBROUTINES = read_control.o \
	      read_geo.o \
              parse_control.o \
	      compare_distances.o \
	      compare_angles.o \
              check_cavity.o \
	      center_structure.o \
	      d1mach.o \
	      check_and_add_list_of_minima.o \
              check_and_add_list_of_minima_v2.o \
	      write_molden_file.o \
	      check_hard_sphere.o \
	      write_status.o \
              get_max_force.o \
	      get_line_force.o \
              get_forces_and_energy.o \
              zufall_sub.o \
              check_compound.o \
	      dgtsv.o \
              xerbla.o \
	      get_pair_energy.o	\
              get_highest_pair_energy.o \
	      get_index_of_minimum.o \
	      binomial_coefficient.o \
              factorial.o \
	      rotate.o \
	      diff_distances.o \
	      nn_optimize.o \
	      pre_relax_hard_sphere.o \
	      check_single_dissociated.o

# Begin make

all: effernir.250209.x

alpha: effernir.250209.alpha.x

effernir.250209.x: $(MODULES) main.o $(SUBROUTINES) 
	$(FC) $(F90FLAGS) -o $@ main.o $(SUBROUTINES) $(MODULES) $(ARCH) $(LIBRARIES) ; \
	cp effernir.250209.x ~/codes/fhi-aims/bin 

effernir.250209.alpha.x: $(MODULES) main.o $(SUBROUTINES) 
	$(FC) $(F90FLAGS) -o $@ main.o $(SUBROUTINES) $(MODULES) $(ARCH) $(LIBRARIES)

# Modules

vector_class.mod: vector.o
vector_array_class.mod: vector_array.o
cluster.mod: cluster.o
basin_hopping.mod: basin_hopping.o
control.mod: control.o
geometry.mod: geometry.o
geometry_list.mod: geometry_list.o
dft_file.mod: dft_file.o
arch_specific.mod: arch_specific.o
conjugate_gradient.mod: conjugate_gradient.o
lennard_jones.mod: lennard_jones.o
map_to_one_dimension.mod: map_to_one_dimension.o
spline.mod: spline.o
pair_potential.mod: pair_pot.o
cartesian_ylm.mod: cartesian_ylm.o

vector.o: vector.f90
	rm vector.mod ; \
        $(FC) $(F90FLAGS) -c vector.f90

vector_array.o: vector_array.f90 vector_class.mod
	rm vector_array.mod ; \
	$(FC) $(F90FLAGS) -c vector_array.f90

cluster.o: cluster.f90 vector_array_class.mod control.mod
	rm cluster.mod ; \
	$(FC) $(F90FLAGS) -c cluster.f90

basin_hopping.o: basin_hopping.f90 vector_array_class.mod cluster.mod control.mod \
	geometry.mod geometry_list.mod dft_file.mod arch_specific.mod conjugate_gradient.mod zufall_sub.o \
	pair_potential.mod
	rm basin_hopping.mod ; \
	$(FC) $(F90FLAGS) -c basin_hopping.f90

control.o: control.f90 
	rm control.mod ; \
	$(FC) $(F90FLAGS) -c control.f90

geometry.o: geometry.f90 vector_array_class.mod cluster.mod cartesian_ylm.mod
	rm geometry.mod ; \
	$(FC) $(F90FLAGS) -c geometry.f90

geometry_list.o: geometry_list.f90 geometry_class.mod cluster.mod control.mod
	rm geometry_list.mod ; \
	$(FC) $(F90FLAGS) -c geometry_list.f90

dft_file.o: dft_file.f90 cluster.mod
	rm dft_file.mod ; \
	$(FC) $(F90FLAGS) -c dft_file.f90

arch_specific.o: $(ARCH)
	rm arch_specific.mod ; \
        $(FC) -o $@ $(F90FLAGS) $(CBFLAGS) -c $(ARCH)

conjugate_gradient.o: conjugate_gradient.f90 vector_array_class.mod map_to_one_dimension.mod cluster.mod lennard_jones.mod
	rm conjugate_gradient.mod ; \
	$(FC) $(F90FLAGS) -c conjugate_gradient.f90

lennard_jones.o: lennard_jones.f90 cluster.mod vector_array_class.mod
	rm lennard_jones.mod ; \
	$(FC) $(F90FLAGS) -c lennard_jones.f90

map_to_one_dimension.o: map_to_one_dimension.f90 vector_array_class.mod cluster.mod dft_file.mod
	rm map_to_one_dimension.mod ; \
	$(FC) $(F90FLAGS) -c map_to_one_dimension.f90

spline.o: spline.f 
	rm spline.mod ; \
	$(FC) $(FFLAGS) -c spline.f

pair_pot.o: pair_pot.f90 cluster.mod spline.mod
	rm pair_pot.mod ; \
	$(FC) $(F90FLAGS) -c pair_pot.f90

cartesian_ylm.o: cartesian_ylm.f90 cluster.mod
	rm cartesian_ylm.mod ; \
	$(FC) $(F90FLAGS) -c cartesian_ylm.f90

# Subroutines, in order of appearance

main.o: main.f90 constants.f90 vector_array_class.mod cluster.mod control.mod dft_file.mod basin_hopping.mod pair_potential.mod \
	cartesian_ylm.mod
	$(FC) $(F90FLAGS) -c main.f90 

read_control.o: read_control.f90 cluster.mod control.mod basin_hopping.mod geometry_class.mod pair_potential.mod
	$(FC) $(F90FLAGS) -c read_control.f90

read_geo.o: read_geo.f90 constants.f90 cluster.mod
	$(FC) $(F90FLAGS) -c read_geo.f90

parse_control.o: parse_control.f90 cluster.mod pair_potential.mod
	$(FC) $(F90FLAGS) -c parse_control.f90

write_molden_file.o: write_molden_file.f90
	$(FC) $(F90FLAGS) -c write_molden_file.f90

compare_distances.o: compare_distances.f90	
	$(FC) $(F90FLAGS) -c compare_distances.f90

compare_angles.o: compare_angles.f90 cluster.mod control.mod geometry.mod	
	$(FC) $(F90FLAGS) -c compare_angles.f90

check_and_add_list_of_minima.o: check_and_add_list_of_minima.f90	
	$(FC) $(F90FLAGS) -c check_and_add_list_of_minima.f90

check_and_add_list_of_minima_v2.o: check_and_add_list_of_minima_v2.f90	
	$(FC) $(F90FLAGS) -c check_and_add_list_of_minima_v2.f90

check_cavity.o: check_cavity.f90	
	$(FC) $(F90FLAGS) -c check_cavity.f90

center_structure.o: center_structure.f90	
	$(FC) $(F90FLAGS) -c center_structure.f90

d1mach.o: d1mach.f	
	$(FC) $(FFLAGS) -c d1mach.f

check_hard_sphere.o: check_hard_sphere.f90 vector_array_class.mod control.mod
	$(FC) $(F90FLAGS) -c check_hard_sphere.f90

write_status.o: write_status.f90
	$(FC) $(F90FLAGS) -c write_status.f90

get_forces_and_energy.o: get_forces_and_energy.f90 control.mod vector_array_class.mod lennard_jones.mod
	$(FC) $(F90FLAGS) -c get_forces_and_energy.f90

get_max_force.o: get_max_force.f90  vector_array_class.mod
	$(FC) $(F90FLAGS) -c get_max_force.f90

get_line_force.o: get_line_force.f90  vector_array_class.mod
	$(FC) $(F90FLAGS) -c get_line_force.f90

zufall_sub.o: zufall_sub.f 
	$(FC) $(FFLAGS) -c zufall_sub.f

check_compound.o: check_compound.f90 vector_array_class.mod cluster.mod
	$(FC) $(F90FLAGS) -c check_compound.f90 

dgtsv.o: dgtsv.f 
	$(FC) $(FFLAGS)	-c dgtsv.f

xerbla.o: xerbla.f
	$(FC) $(FFLAGS)	-c xerbla.f

get_pair_energy.o: get_pair_energy.f90 cluster.mod pair_potential.mod lennard_jones.mod vector_array_class.mod
	$(FC) $(F90FLAGS) -c get_pair_energy.f90

get_highest_pair_energy.o: get_highest_pair_energy.f90 cluster.mod vector_array_class.mod
	$(FC) $(F90FLAGS) -c get_highest_pair_energy.f90

get_index_of_minimum.o: get_index_of_minimum.f90 geometry_class.mod geometry_list.mod cluster.mod
	$(FC) $(F90FLAGS) -c get_index_of_minimum.f90

binomial_coefficient.o: binomial_coefficient.f90 
	$(FC) $(F90FLAGS) -c binomial_coefficient.f90

factorial.o: factorial.f90
	$(FC) $(F90FLAGS) -c factorial.f90

rotate.o: rotate.f90 vector_class.mod
	$(FC) $(F90FLAGS) -c rotate.f90

diff_distances.o: diff_distances.f90 control.mod cluster.mod geometry_class.mod
	$(FC) $(F90FLAGS) -c diff_distances.f90

nn_optimize.o: nn_optimize.f90 vector_array_class.mod cluster.mod
	$(FC) $(F90FLAGS) -c nn_optimize.f90

pre_relax_hard_sphere.o: pre_relax_hard_sphere.f90 vector_array_class.mod control.mod
	$(FC) $(F90FLAGS) -c pre_relax_hard_sphere.f90 

check_single_dissociated.o: check_single_dissociated.f90 vector_array_class.mod
	$(FC) $(F90FLAGS) -c check_single_dissociated.f90

clean:
	rm -f *.o *.mod ; 
