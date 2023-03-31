test_DIR = $(TOP_DIR)/test
VPATH := $(VPATH):$(test_DIR)

# Define common codes used for the tests
TEST_COMMON = tst_ncdf_utils.o

TEST_EXE_  = tst_ncdf_3 tst_ncdf_4 tst_ncdf_att
TEST_EXE_ += tst_ncdf_crt

TEST ?= $(TEST_EXE_)

# Create the objects
TEST_OBJECTS = $(addsuffix .o, $(TEST))


ifeq ($(STATIC),1)
  TEST_LIB    = $(NCDF_LIB_STATIC)
else
  ifeq ($(SHARED),1)
    TEST_LIB  = $(NCDF_LIB_SHARED)
  else
    $(error No library is build, the tests cannot be performed.)
  endif
endif


# Create correct dependency schemes without adding too
# many flags in the direct compilation targets.
# Ensure the common utility is compiled before any test
$(TEST_OBJECTS): $(TEST_COMMON)
# Ensure the executables depend on the utility and the library
$(TEST): $(TEST_COMMON) $(TEST_LIB)
# There are *.inc files in the src directory that is
# needed
$(TEST): INCLUDES += -I$(src_DIR)


# Add source includes

.PHONY: test
.NOTPARALLEL: test
test: $(TEST_LIB) $(TEST)


tst_ncdf_3: tst_ncdf_3.o
	$(FC) -o $@ $(FFLAGS) $(INCLUDES) $^ $(LIBS) $(LDFLAGS)
	./$@
ifeq ($(MPI),1)
	mpirun -np 2 ./$@
endif


tst_ncdf_4: tst_ncdf_4.o
	$(FC) -o $@ $(FFLAGS) $(INCLUDES) $^ $(LIBS) $(LDFLAGS)
	./$@
ifeq ($(MPI),1)
	mpirun -np 2 ./$@
endif


tst_ncdf_att: tst_ncdf_att.o
	$(FC) -o $@ $(FFLAGS) $(INCLUDES) $^ $(LIBS) $(LDFLAGS)
	./$@
ifeq ($(MPI),1)
	mpirun -np 2 ./$@
endif


# This should not run in parallel, however, it could be performed.
tst_ncdf_crt: tst_ncdf_crt.o
	$(FC) -o $@ $(FFLAGS) $(INCLUDES) $^ $(LIBS) $(LDFLAGS)
	./$@


.PHONY: clean-test
clean-test:
	-rm -f $(TEST) $(LIB_TEST) $(TEST_COMMON) *.o *.mod *.nc

clean: clean-test


# Local Variables:
#  mode: makefile-gmake
# End:
