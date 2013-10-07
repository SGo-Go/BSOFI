# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile settings for Atlas.
# ------------------------------------------------------------------

ATLAS_PATH    ?= /usr/local
ATLAS_LIB      = $(ATLAS_PATH)/lib
ATLAS_CXXFLAGS = -DHAS_ATLAS
ATLAS_LIBFLAGS = -L$(ATLAS_LIB) -lf77refblas -latlas #-lcblas
#ATLAS_LIBFLAGS = $(patsubst %, $(ATLAS_LIB)/lib%.a, atlas f77refblas)

#BLAS_CXXFLAGS   = -DHAS_BLAS
BLAS_LIBFLAGS   = -L$(ATLAS_LIB) -lf77blas -lcblas -latlas
LAPACK_CXXFLAGS = -DHAS_BLAS -D HAS_LAPACK -DHAS_ATLAS
LAPACK_CCFLAGS  = -DHAS_BLAS -D HAS_LAPACK -DHAS_ATLAS
LAPACK_LIBFLAGS:= -L$(ATLAS_LIB)  -llapack -lf77blas -lcblas -latlas
#LAPACK_LIBFLAGS = $(ATLAS_LIB)/liblapack.a
