# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile settings for Lapack from Netlib.
# ------------------------------------------------------------------

LAPACK_PATH    ?= /usr/local/lib
LAPACK_LIB      = $(LAPACK_PATH)

BLAS_CXXFLAGS   = -DHAS_BLAS
BLAS_CCFLAGS    = -DHAS_BLAS
BLAS_LIBFLAGS   = -L$(LAPACK_LIB) -lrefblas
LAPACK_CXXFLAGS = -DHAS_BLAS -DHAS_LAPACK=LAPACK
LAPACK_CCFLAGS  = -DHAS_BLAS -DHAS_LAPACK=LAPACK
LAPACK_LIBFLAGS = -L$(LAPACK_LIB)  -llapack -lrefblas #-fopenmp #$(OMP_LIBFLAGS)
#LAPACK_LIBFLAGS = $(patsubst %, $(LAPACK_LIB)/lib%.a, lapack refblas)# -fopenmp
