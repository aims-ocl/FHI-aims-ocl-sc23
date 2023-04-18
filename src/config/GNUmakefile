#  Emacs, this is: -*- mode: makefile -*-

#************************************************************************
#
# Read configuration
#
#************************************************************************

# Default configuration is expected in make.cfg.
CFG ?= default
include $(CFG).cfg
override_backend_all: default  # Must come before $(CFG)


#************************************************************************
#
# Include backend
#
#************************************************************************

include Makefile.backend

#************************************************************************
#
# Include main target Makefiles
#
#************************************************************************

Makefile.tmp:Makefile
	sed '/^ *include Makefile.backend/ s/^/# /'  <$< >$@
Makefile.%.tmp:Makefile.%
	sed '/^ *include Makefile.backend/ s/^/# /'  <$< >$@
include Makefile.tmp
include Makefile.shm.tmp
include Makefile.meta.tmp
include Makefile.cuba.tmp

#************************************************************************
#
# Default binary (all configured features)
#
#************************************************************************

# Name of the default binary:
ifeq ($(strip $(BINNAME)),)
  BINARY = $(BINDIR)/aims.$(AIMSVERSION).x
else
  BINARY = $(BINDIR)/aims.$(AIMSVERSION).$(BINNAME).x
endif
default: $(BINARY)

# BUILD accumulates objects to be compiled and linked in.
# LINK contains additional pre-built libraries for the linker.
BUILD = $(OBJDIR)/write_version_stamp.o $(MODULES) $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) $(ANACON) $(METIS_STUB)
LINK = 

# If new object files or libraries need to be added, they should go somewhere
# here.  Be aware that the order matters.
ifeq ($(strip $(USE_MPI)),yes)
  LINK_FC = $(MPIFC)
else
  LINK_FC = $(FC)
  BUILD += $(MPI_STUBS)
endif
ifneq ($(strip $(SHM)),)
  # shm.o should also go into $(OBJDIR):
  BUILD += $(patsubst shm.o,$(OBJDIR)/shm.o,$(SHM))
else
  BUILD += $(SHM_STUBS)
endif
ifneq ($(strip $(CUBA)),)
  LINK += $(CUBA)
else
  BUILD += $(CUBA_STUBS)
endif
ifneq ($(strip $(META)),)
  BUILD += $(PLUMED_OBJECTS)
else
  BUILD += $(PLUMED_STUB)
endif
ifneq ($(strip $(SCALAPACK)),)
  LINK += $(SCALAPACK)
else
  BUILD += $(SCALAPACK_STUBS)
endif
ifeq ($(strip $(PARSER)),yes)
  BUILD += $(PARSER_OBJECTS)
else
  BUILD += $(PARSER_STUBS)
endif
# Trouble is, CC is set to "cc" by default.
# In most cases, this is just fine, 
ifeq ($(origin CC),default)
  BUILD += $(C_FILES_STUBS)
else
  BUILD += $(C_FILES)
endif

# LAPACKBLAS last because things might depend on it.
LINK += $(LAPACKBLAS)

$(BINARY): $(BUILD)
	$(LINK_FC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ $^ $(LINK)
	@echo "Built $@."

%.x: $(BUILDDIR)/%.o $(BUILDDIR)/libaims-$(AIMSVERSION).so
	$(LINK_FC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ $^ $(LINK)
	@echo "Built $@."
