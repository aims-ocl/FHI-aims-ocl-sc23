# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.
SOURCES := $(filter-out %tests.f90,$(subst F90,f90,$(notdir $(wildcard $(SRCDIR)/*.?90))))

VPATH = $(SRCDIR)

all:

dist: $(addprefix src/,$(SOURCES)) version.txt

src/%.f90: %.f90 | src
	cp $^ $@

src/%.f90: %.F90 | src
	gfortran -E -P $^ -o $@ -DWITH_MPI -DWITH_SCALAPACK -DWITH_MPIFH -DWITH_ELSI

version.txt:
	git -C $(SRCDIR)/.. rev-parse HEAD >$@

src:
	mkdir $@

clean:
	rm -f {serial,scalapack_mpi}/*.f90
