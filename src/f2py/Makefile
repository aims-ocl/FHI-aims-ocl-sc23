##############################################################################
# Read & postprocess configuration from ../$(CFG).cfg (../default.cfg)

all:

CFG ?= default
include ../$(CFG).cfg
-include $(CFG).cfg

$(eval $(shell egrep '^ *AIMSVERSION' ../Makefile.backend))

ifeq ($(BUILDDIR),)
  F2PY_BUILDDIR ?= ..
else
  F2PY_BUILDDIR ?= ../$(BUILDDIR)
endif
ifeq ($(OBJDIR),)
  F2PY_OBJDIR ?= $(F2PY_BUILDDIR)
else
  F2PY_OBJDIR ?= ../$(OBJDIR)
endif
ifeq ($(MODDIR),)
  F2PY_MODDIR ?= $(F2PY_BUILDDIR)
  MODDIR = $(F2PY_MODDIR)
else
  F2PY_MODDIR ?= ../$(MODDIR)
endif
AIMS_SO = $(PWD)/$(F2PY_OBJDIR)/libaims-$(AIMSVERSION).so
FFLAGS+=-I$(F2PY_MODDIR)
F90 ?= $(FC)

ifeq ($(strip $(FC)),)
  $(info Variable FC is not set.  Please make sure that you)
  $(info have set the configuration variables in Makefile.)
  $(error Variable FC needs to be set)
endif

ifeq ($(F2PY_FCOMP),)
  ifeq ($(FC),ifort)
    F2PY_FCOMP=intelem
  else ifeq ($(FC),gfortran)
    F2PY_FCOMP=gnu95
  else
    $(error Variable F2PY_FCOMP not set)
  endif
endif

##############################################################################
# Module _ylm

ylm.pyf: ../external/ylm.f ../external/ylm_real.f ../external/drc3jj.f \
         ../external/drc3jm.f ../triple_Y.f90
_ylm.so:$(AIMS_SO)
all: _ylm.so

##############################################################################
# Module _triple_y.so

triple_y.pyf: ../external/drc3jj.f ../external/drc3jm.f ../triple_Y.f90
_triple_y.so:$(AIMS_SO)
all: _triple_y.so

##############################################################################
# Module _bspline

bspline.pyf: ../bspline.f90 val_vec_bspline.f90

_bspline.so:$(AIMS_SO) val_vec_bspline.o
all: _bspline.so

##############################################################################
# Module _wigner_rot

wigner_rot.pyf: ../wigner_small_d.f90
_wigner_rot.so:$(AIMS_SO)
all: _wigner_rot.so

##############################################################################
# check_auxmat.x

check_periodic_auxmat.x: init_aims.o
check_auxmat.x: init_aims.o
all: check_auxmat.x

##############################################################################
# others ...

test_basis.x: init_aims.o
test_eigen.x: init_aims.o
test_lvlpairs.x: init_aims.o

##############################################################################
# General compilation rules

# %.mod:%.o
# 	touch $@

%.o:%.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o:%.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

%.x:%.o $(AIMS_SO)
	$(MPIFC) $(F90FLAGS) $(LDFLAGS) -o $@ $^ $(LAPACKBLAS)

gen_%.pyf:../%.f90
	-test -f $@ && rm $@
	f2py -h $@ -m _gen_$* $^
	sed -i '/character.*allocatable/d' $@  # char...*10 a(:) breaks f2py
	sed -i '/\w\{45\}/d' $@                # long names break f2py
	sed -i '/^ \+type(/d' $@               # types are not really supported

%.pyf:
	test -n "$^" || { echo "No definition for $@ given."; false; }
	-test -f $@ && rm $@
	f2py -h $@ -m _$* $^
	sed -i '/character.*allocatable/d' $@  # char...*10 a(:) breaks f2py
	sed -i '/\w\{45\}/d' $@                # long names break f2py
	sed -i '/^ \+type(/d' $@               # types are not really supported

_%.so:%.pyf dummy.f90 $(AIMS_SO)
	f2py --link-lapack --fcompiler=$(F2PY_FCOMP) \
             --f90flags="-I$(F2PY_MODDIR)" \
             -c $< $(filter-out %.mod %.pyf,$^)


##############################################################################
# Auxiliary targets


clean:
	rm -f *.o *.so *.pyc *.mod *.x TAGS *.pyf *~

TAGS:../TAGS

../TAGS: FORCE
	cd ..; $(MAKE) TAGS

lib:
	cd ..; $(MAKE) lib

.SECONDARY:%.o
.PHONY: FORCE
