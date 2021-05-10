
SHELL := /bin/bash
COMPILE_DIR = $(shell pwd)

# Check if compile directory exists, returns error code and exit otherwise
EXISTS = $(firstword $(wildcard $(COMPILE_DIR)) )
$(info $$COMPILE_DIR is [${COMPILE_DIR}] )
ifeq ($(EXISTS), )
      $(error ${COMPILE_DIR} directory does not exist)
endif

# Basic options
#   - COMPUTER..............: current choice is DESKTOP, QUARTZ, LASSEN for preset configurations of BLAS and LAPACK libraries
#   - COMPILER .............: type of compiler: GFORTRAN, IFORT, CRAY or IBM for preset compiler options
#   - FORTRAN_MPI ..........: actual name of the Fortran compiler. Standard choices are gfortran, gfortran, ifort, xlf-gpu.
#                             Attention:
#                               - On Lassen, it seems necessary to use a compiler with suffix -gpu
#                               - In order to use HDF5, one should most likely use FORTRAN_MPI = h5pfc
COMPUTER    = DESKTOP
COMPILER    = IFORT
FORTRAN_MPI = ifort

# Build options (TRUE or FALSE)
#   - DEBUG ................: activates debugging options for each compiler
#   - VERBOSE ..............: enables verbose compilation
DEBUG   = FALSE
VERBOSE = FALSE

#  HPC options (integer numbers, 0: inactive, 1: active)
#   - USE_ESSL .............: uses IBM ESSL library for LAPACK and BLAS (may require SWITCH_DIAG=2)
#   - USE_OPENMP ...........: enables OpenMP multithreading
#   - USE_MPI ..............: enables MPI distributed parallelism
#   - USE_LAPACK ...........: links to Lapack library for advanced linear algebra
#   - USE_BLAS .............: links to BLAS library for basic linear algebra operations
#   - USE_MKL ..............: links to Intel MKL implementation of BLAS and LAPACK
#   - USE_MAGMA ............: links to MAGMA for a GPU implementation of BLAS and LAPACK
#   - USE_MANYCORES ........: uses several MPI tasks per HFB configuration, number is M_GRID x N_GRID
#   - USE_SCALAPACK ........: uses ScaLapack for HFB matrix diagonalization in ISIMPY=0 case (may require SWITCH_DIAG=4)
USE_ESSL      = 0
USE_OPENMP    = 1
USE_MPI       = 0
USE_LAPACK    = 0
USE_BLAS      = 0
USE_MKL       = 1
USE_MAGMA     = 0
USE_SCALAPACK = 0
USE_MANYCORES = 0
ifeq ($(USE_MANYCORES),1)
      M_GRID  = 2
      N_GRID  = 1
else
      M_GRID  =
      N_GRID  =
endif

#======================================================================#
#  LIBRARIES
#======================================================================#

# OpenMP
OPENMP_LIBDIR  =
OPENMP_LIBRARY = gomp

# Lapack
LAPACK_INCLUDE =
LAPACK_LIBDIR  = $(HOME)/local
LAPACK_LIBRARY = lapack

# Blas
BLAS_INCLUDE   =
BLAS_LIBDIR    = $(HOME)/local
BLAS_LIBRARY   = blas

# MKL
MKL_VERSION    = mkl-2020.0
MKL_INCLUDE    = $(MKLROOT)/include
MKL_LIBDIR     = $(MKLROOT)/lib/intel64
MKL_LIBRARY    = mkl_intel_lp64 mkl_intel_thread mkl_core mkl_avx2

# Magma
MAGMA_INCLUDE  =
MAGMA_LIBDIR   = $(MAGMA_DIR)/lib
MAGMA_LIBRARY  = magma

# ESSL
ESSL_INCLUDE   = $(ESSLHEADERDIR)
ESSL_LIBDIR    = $(ESSLLIBDIR64)
ESSL_LIBRARY   = esslsmp lapackforessl

# ScaLapack
SCALAPACK_INCLUDE = $(HOME)/local
SCALAPACK_LIBDIR  = $(HOME)/local
SCALAPACK_LIBRARY = scalapack

# Preset configurations
ifeq ($(COMPUTER),DESKTOP)
      USE_LAPACK = 0
      USE_BLAS   = 0
      USE_MKL    = 1
      USE_MAGMA  = 1 
endif
ifeq ($(COMPUTER),QUARTZ)
      USE_LAPACK = 0
      USE_BLAS   = 0
      USE_MKL    = 1
      USE_MAGMA  = 0

      OPENMP_LIBRARY = iomp5
      MKL_VERSION    = mkl-2020.0
      MKL_INCLUDE    = /usr/tce/packages/mkl/$(MKL_VERSION)/include
      MKL_LIBDIR     = /usr/tce/packages/mkl/$(MKL_VERSION)/lib
      MKL_LIBRARY    = mkl_intel_lp64 mkl_intel_thread mkl_core mkl_avx2 iomp5
endif
ifeq ($(COMPUTER),LASSEN)
      USE_LAPACK = 0
      USE_BLAS   = 0
      USE_MKL    = 0
      USE_MAGMA  = 0
      USE_ESSL   = 1

      OPENMP_LIBRARY = xlsmp
      MAGMA_LIBDIR   = $(MAGMA_DIR)/lib
      MAGMA_LIBRARY  = magma
      ESSL_INCLUDE   = $(ESSLHEADERDIR)
      ESSL_LIBDIR    = $(ESSLLIBDIR64)
      ESSL_LIBRARY   = esslsmp lapackforessl
endif

#======================================================================#
#  Nothing beyond this line should be changed, in principle            #
#======================================================================#

# Define the preprocessor variables for each platform
PREPROC_VAR = -DUSE_OPENMP=$(USE_OPENMP) \
              -DUSE_MPI=$(USE_MPI) \
              -DUSE_MANYCORES=$(USE_MANYCORES) \
              -DUSE_HDF5=$(USE_HDF5) \
              -DSWITCH_QUAD=$(SWITCH_QUAD) \
              -DSWITCH_PORT=$(SWITCH_PORT) \
              -DSWITCH_DIAG=$(SWITCH_DIAG) \
              -DSWITCH_VECT=$(SWITCH_VECT) \
              -DUSE_SCALAPACK=$(USE_SCALAPACK) \
              -DM_GRID=$(M_GRID) \
              -DN_GRID=$(N_GRID) \
              -DUSE_INERTIA=$(USE_INERTIA)\
              -DUSE_COULEX=$(USE_COULEX) \
              -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
              -DGOGNY_HYPER=$(GOGNY_HYPER) \
              -DUSE_ESSL=$(USE_ESSL) \
              -DUSE_LOCALI=$(USE_LOCALI) \
              -DHFODD_USE_CPP=$(USE_CPP)

# Defining compiler options for: INTEL FORTRAN COMPILER (ifort)
ifeq ($(COMPILER),IFORT)
      PRECISION    = -i4 -r8
      FORMAT_F77   = -fixed -80
      FORMAT_F90   = -free -extend_source
      FCFLAGS      = -cpp -fPIC $(PREPROC_VAR)
      LDFLAGS      = -fPIC
      LDLIBS       =
      ifeq ($(DEBUG),FALSE)
            FCFLAGS += -O3
      else
            FCFLAGS += -g -check all -warn nodeclarations -warn nounused
      endif
      ifeq ($(STATIC_LIB),TRUE)
            LDFLAGS += -static
      endif
      ifeq ($(USE_OPENMP),1)
            FCFLAGS += -qopenmp
            LDFLAGS += -qopenmp
	endif
      ifeq ($(LARGE_MEMORY),TRUE)
            FCFLAGS += -mcmodel=medium -shared-intel
      endif
      ifeq ($(USE_CPP),1)
            LDLIBS += -lstdc++
      endif
endif

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER),GFORTRAN)
      PRECISION    = -fdefault-real-8 -fdefault-double-8
      FORMAT_F77   = -ffixed-form
      FORMAT_F90   = -ffree-form -ffree-line-length-none
      FCFLAGS      = -cpp -fpic $(PREPROC_VAR)
      LDFLAGS      = -fpic
      LDLIBS       =
      ifeq ($(DEBUG),FALSE)
            FCFLAGS += -O3 -march=native -mavx2 -mfma -ffp-contract=fast
            #FCFLAGS += -O1
      else
            FCFLAGS += -g -O0 -Wall \
                       -Warray-bounds -Wunderflow -Warray-temporaries \
                       -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check \
                       -Wconversion -Wuninitialized -pedantic \
                       -finit-real=nan \
                       -ftrapv
      endif
      ifeq ($(STATIC_LIB),TRUE)
            LDFLAGS += -static
      endif
      ifeq ($(USE_OPENMP),1)
            FCFLAGS += -fopenmp
      endif
      ifeq ($(LARGE_MEMORY),TRUE)
            FCFLAGS += -mcmodel=medium
      endif
      ifeq ($(USE_CPP),1)
            LDLIBS += -lstdc++
      endif
endif

# Defining compiler options for: CRAY FORTRAN COMPILER
ifeq ($(COMPILER),CRAY)
      PRECISION  = -s real64 -d p
      FORMAT_F77 = -f fixed
      FORMAT_F90 = -f free
      FCFLAGS    = -e Z -h Pic $(PREPROC_VAR)
      LDFLAGS    = -h Pic
      LDLIBS     =
      ifeq ($(DEBUG),FALSE)
            FCFLAGS += -O3
      else
            FCFLAGS += -e c -e D
      endif
      ifeq ($(STATIC_LIB),TRUE)
            LDFLAGS += --staticlink
      endif
      ifeq ($(USE_OPENMP),0)
            FCFLAGS += -h noomp
      endif
      ifeq ($(USE_CPP),1)
            LDLIBS += -lstdc++
      endif
endif

# Defining compiler options for: IBM FORTRAN COMPILER (xlf90)
ifeq ($(COMPILER),IBM)
      PRECISION    = -qautodbl=dbl4
      FORMAT_F77   = -qfixed -qsuffix=cpp=f
      FORMAT_F90   = -qfree=f90 -qsuffix=cpp=f90
      FCFLAGS      = -qnoescape -qpic $(foreach opt,$(PREPROC_VAR),'-WF,$(opt)')
      LDFLAGS      = -qpic
      LDLIBS       =
      ifeq ($(DEBUG),FALSE)
            FCFLAGS += -qstrict -O2 -qarch=pwr9 -qtune=pwr9
      else
            FCFLAGS += -g9 -C -qflttrap -qcheck=all
      endif
      ifeq ($(USE_OPENMP),1)
            FCFLAGS += -qsmp=omp
      endif
      ifeq ($(USE_CPP),1)
            LDLIBS += -lstdc++
      endif
endif

#======================================================================#
#  VARIABLES FOR COMPILER OPTIONS AND LIBRARIES
#======================================================================#

# Name of the main executable
TEST_DENSHF_EXE = test_denshf

# Add libraries
ifeq ($(USE_ESSL),1)
      CXXFLAGS += $(addprefix -I,$(ESSL_INCLUDE))
      FCFLAGS  += $(addprefix -I,$(ESSL_INCLUDE))
      LDFLAGS  += $(addprefix -L,$(ESSL_LIBDIR)) $(addprefix -R,$(ESSL_LIBDIR))
      LDLIBS   += $(addprefix -l,$(ESSL_LIBRARY))
endif
ifeq ($(USE_OPENMP),1)
      LDFLAGS  += $(addprefix -L,$(OPENMP_LIBDIR))
      LDLIBS   += $(addprefix -l,$(OPENMP_LIBRARY))
endif
ifeq ($(USE_LAPACK),1)
      CXXFLAGS += $(addprefix -I,$(LAPACK_INCLUDE))
      FCFLAGS  += $(addprefix -I,$(LAPACK_INCLUDE))
      LDFLAGS  += $(addprefix -L,$(LAPACK_LIBDIR))
      LDLIBS   += $(addprefix -l,$(LAPACK_LIBRARY))
endif
ifeq ($(USE_BLAS),1)
      CXXFLAGS += $(addprefix -I,$(BLAS_INCLUDE))
      FCFLAGS  += $(addprefix -I,$(BLAS_INCLUDE))
      LDFLAGS  += $(addprefix -L,$(BLAS_LIBDIR))
      LDLIBS   += $(addprefix -l,$(BLAS_LIBRARY))
endif
ifeq ($(USE_MKL),1)
      CXXFLAGS += $(addprefix -I,$(MKL_INCLUDE))
      FCFLAGS  += $(addprefix -I,$(MKL_INCLUDE))
      LDFLAGS  += $(addprefix -L,$(MKL_LIBDIR))
      LDLIBS   += $(addprefix -l,$(MKL_LIBRARY))
endif
ifeq ($(USE_MAGMA),1)
      CXXFLAGS += $(addprefix -I,$(MAGMA_INCLUDE))
      FCFLAGS  += $(addprefix -I,$(MAGMA_INCLUDE))
      LDFLAGS  += $(addprefix -L,$(MAGMA_LIBDIR))
      LDLIBS   += $(addprefix -l,$(MAGMA_LIBRARY))
endif
ifeq ($(USE_SCALAPACK),1)
      CXXFLAGS += $(addprefix -I,$(SCALAPACK_INCLUDE))
      FCFLAGS  += $(addprefix -I,$(SCALAPACK_INCLUDE))
      LDFLAGS  += $(addprefix -L,$(SCALAPACK_LIBDIR))
      LDLIBS   += $(addprefix -l,$(SCALAPACK_LIBRARY))
endif

#======================================================================#
#  VARIABLE EXPORTS
#======================================================================#

# Export use options
export USE_MPI
export USE_OPENMP
export USE_ESSL
export USE_LAPACK
export USE_BLAS
export USE_MKL
export USE_MAGMA
export USE_SCALAPACK
# Export main compiler options
export COMPUTER
export FORTRAN_MPI
export COMPILER
export PRECISION
export FORMAT_F90
export FCFLAGS_EXT  ::= $(FCFLAGS)
export LDFLAGS_EXT  ::= $(LDFLAGS)
export LDLIBS_EXT   ::= $(LDLIBS)
export LINEAR_ALGEBRA
export VERBOSE

#======================================================================#
#  BUILDING
#======================================================================#

build = $(TEST_DENSHF_EXE)

# Build test_denshf
all: $(build)

# Cleaning
clean ::
	$(MAKE) clean -C doc
	-rm -f *.o *.i *.lst *.oo *.ipo *.ipa *.mod *_PREPROCESSED*

# Deep cleaning
mrproper ::
	$(MAKE) clean -C doc
	-rm -f *.o *.i *.lst *.oo *.ipo *.ipa *.mod *_PREPROCESSED* $(TEST_DENSHF_EXE)

#======================================================================#
#  DEPENDENCIES
#======================================================================#

# Build required dependencies for each source file
HFODD_TYPES_DEPS =
HFODD_SIZES_DEPS =
BASIS_DEPS = hfodd_sizes.o hfodd_types.o
TEST_DENSHF_DEPS = basis.o hfodd_sizes.o hfodd_types.o
TEST_DENSHF_ALL_DEPS = basis.o hfodd_sizes.o hfodd_types.o

#======================================================================#
#  RULES
#======================================================================#

hfodd_types.o : hfodd_types.f90 $(HFODD_TYPES_DEPS)
hfodd_sizes.o : hfodd_sizes.f90 $(HFODD_SIZES_DEPS)
basis.o       : basis.f90 $(BASIS_DEPS)
test_denshf.o : test_denshf.f90 $(TEST_DENSHF_DEPS)

%.o: %.f90
ifeq ($(VERBOSE),1)
	$(FORTRAN_MPI) -c $< $(FORMAT_F90) $(PRECISION) $(FCFLAGS)
else
	@echo "Compiling \"$@\"..."
	@$(FORTRAN_MPI) -c $< $(FORMAT_F90) $(PRECISION) $(FCFLAGS)
endif

# executable
$(TEST_DENSHF_EXE): $(TEST_DENSHF_ALL_DEPS) test_denshf.o
ifeq ($(VERBOSE),1)
	$(FORTRAN_MPI) $(LDFLAGS) $(FORMAT_F90) $(PRECISION) -o $@ $^ $(LDLIBS)
else
	@echo "Linking \"$@\"..."
	@$(FORTRAN_MPI) $(LDFLAGS) $(FORMAT_F90) $(PRECISION) -o $@ $^ $(LDLIBS)
endif

