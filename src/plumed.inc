# This is included in the Makefile for the metadynamics
# plugin PLUMED compilation

PLUMED_OBJECTS_U = \
	biasexchange.o \
	dafed.o \
	hills.o \
	metadyn.o \
	ptmetad.o \
	read_restraint.o \
	restraint.o \
	restraint_adf.o \
	restraint_alfabeta.o \
	restraint_alpharmsd.o \
	restraint_angle.o \
	restraint_antibetarmsd.o \
	restraint_cmap.o \
	restraint_coord.o \
	restraint_dihcor.o \
	restraint_dipole.o \
	restraint_dist.o \
	restraint_elstpot.o \
	restraint_energy.o \
	restraint_func.o \
	restraint_hbonds.o \
	restraint_helix.o \
	restraint_histogram.o \
	restraint_mindist.o \
	restraint_parabetarmsd.o \
	restraint_pca.o \
	restraint_poly.o \
	restraint_position.o \
	restraint_puckering.o \
	restraint_rdf.o \
	restraint_rgyr.o \
	restraint_spath.o \
	restraint_sprint.o \
	restraint_torsion.o \
	restraint_waterbridge.o \
	restraint_zpath.o \
	testderivatives.o

PLUMED_OBJECTS:=$(addprefix $(OBJDIR)/,$(PLUMED_OBJECTS_U))

$(OBJDIR)/%.o: external/PLUMED-1.3.0/%.c
	$(CC) $(CFLAGS) -DFHIAIMS -c $< -o $@
