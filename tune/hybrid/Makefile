# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile of benchmark for hybrid CPU+GPU implementation 
# ------------------------------------------------------------------

include ../../config.mk

THIS_DIR = $(TOP_DIR)/tune/hybrid
TARGET   = tune

SOURCES   =  tune_hybrid.c
CXXFLAGS += #-DUSE_PROF

INCLUDES += -I$(TOP_DIR) -I$(TOP_DIR)/bench/include -I$(TOP_DIR)/include
CXXFLAGS += $(LAPACK_CXXFLAGS) $(CUDA_CXXFLAGS)
CFLAGS   += $(LAPACK_CXXFLAGS) $(CUDA_CFLAGS) #$(CUDA_CXXFLAGS)

LIBS     += $(BSOFI_LIB) $(LAPACK_LIBFLAGS) \
	$(CUDA_LIBFLAGS) $(MAGMA_LIBFLAGS) #$(CXXLIBS)

HEADERS   = setup.h $(wildcard $(TOP_DIR)/tune/include/*.h) $(wildcard $(TOP_DIR)/include/bench/*.h)
HEADERSXX = $(HEADERS) #$(wildcard $(TOP_DIR)/include/*.hpp)

include $(TOP_DIR)/config/makefile-general-rules.mk
