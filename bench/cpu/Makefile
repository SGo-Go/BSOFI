# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile of benchmark for pure CPU implementation 
# ------------------------------------------------------------------

include ../../config.mk

THIS_DIR = $(TOP_DIR)/bench/cpu
TARGET   = bench

SOURCES   =  cpu_bench.cpp
CXXFLAGS += -DUSE_PROF

INCLUDES += -I$(TOP_DIR) -I$(TOP_DIR)/bench/include -I$(TOP_DIR)/include
CXXFLAGS += $(LAPACK_CXXFLAGS) #$(CUDA_CXXFLAGS) $(MAGMA_CXXFLAGS)
CFLAGS   += $(LAPACK_CXXFLAGS) 

LIBS     += $(BSOFI_LIB) $(LAPACK_LIBFLAGS) # -limf -lirc $(CXXLIBS) \
# $(CUDA_LIBFLAGS) $(MAGMA_LIBFLAGS)

HEADERS   = setup.h $(wildcard $(TOP_DIR)/bench/include/*.h)
HEADERSXX = $(HEADERS) #$(wildcard $(TOP_DIR)/include/*.hpp)

include $(TOP_DIR)/config/makefile-general-rules.mk
