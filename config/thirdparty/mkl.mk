# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile settings for MKL.
# ------------------------------------------------------------------

MKL_PATH      ?= /opt/intel/mkl
MKL_LIB        = $(MKL_PATH)/lib/intel64
MKL_CXXFLAGS   = -DHAS_MKL
MKL_LIBFLAGS   = -lrt -L$(MKL_PATH) \
	-lmkl_core -liomp5 -lpthread -lmkl_intel_thread

BLAS_LIBFLAGS   = $(MKL_LIBFLAGS) -lmkl_intel_lp64
LAPACK_CXXFLAGS = -DHAS_BLAS -D HAS_LAPACK=1 -DMKL_LAPACK
LAPACK_CCFLAGS  = -DHAS_BLAS -D HAS_LAPACK=1 -DMKL_LAPACK
LAPACK_LIBFLAGS = $(MKL_LIBFLAGS) -lmkl_intel_lp64
