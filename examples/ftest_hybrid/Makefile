# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile that shows how to use BSOF/I CPU+GPU codes with Fortran
# ------------------------------------------------------------------

TOP_DIR = ../..
include ../../config.mk

THIS_DIR = $(TOP_DIR)/examples/ftest_hybrid
TARGET   = test

SOURCES   = ftest_hybrid.f90

LINKER    = $(FC)

INCLUDES += -I$(TOP_DIR) -I$(TOP_DIR)/include

CXXFLAGS += -DUSE_PROF $(LAPACK_CXXFLAGS) \
	$(CUDA_CXXFLAGS) $(MAGMA_CXXFLAGS) # $(DBG_CXXFLAGS) 

CFLAGS   += $(LAPACK_CXXFLAGS) \
	$(CUDA_CXXFLAGS) $(MAGMA_CXXFLAGS) # $(DBG_CXXFLAGS) 

LIBS     +=  $(BSOFI_LIB) \
	$(LAPACK_LIBFLAGS) $(CXXLIBS) $(CUDA_LIBFLAGS) $(MAGMA_LIBFLAGS)

HEADERS   = bench_bso.hpp gemm.hpp setup.hpp #$(wildcard $(TOP_DIR)/include/*.h)
HEADERSXX = $(HEADERS) $(wildcard $(TOP_DIR)/include/*.hpp)

#include ../config/makefile-general-rules.mk

OBJ     ?= $(patsubst %.f90, %.o, \
	$(patsubst %.f, %.o, \
	$(patsubst %.c, %.o, $(SOURCES:.cpp=.o))))
LIBS    += $(if $(filter %.f, $(SOURCES)), $(FLIBS),)
LINKER  ?= $(CXX)
.PHONY  : clean install all #$(TARGET) 
ifneq   "$(TARGET)" ""
all: $(TARGET)
else
all:;@echo 'Nothing to make!'
endif

$(TARGET): $(OBJ)
	$(LINKER) $(^) $(LIBS) -o $@

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp $(HEADERSXX)
	$(CXX) $(CXXFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

#$(OBJ:.o=$(THIS_DIR)/.o)
# $(addprefix $(THIS_DIR)/, $(OBJ))

clean:
	$(RM) $(OBJ) 
	$(RM) $(addprefix $(THIS_DIR)/, *~ $(TARGET))

install: $(TARGET)
	$(MKDIR) $(INSTALL_DIR)
	$(CP) $(TARGET) $(INSTALL_DIR)/$(TARGET)
