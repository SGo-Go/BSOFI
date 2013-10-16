# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Makefile for library.
# ------------------------------------------------------------------

TOP_DIR = ..
include ../config.mk

THIS_DIR = $(TOP_DIR)/src
TARGET   = $(BSOFI_LIB)

SOURCES   = bsofi_cpu.c bsofi_callback.c bsofi_lapack.c 

#ifneq "$(CUDA_PATH)" ""
SOURCES  += bsofi_hybrid.c
#endif

ifneq "$(MAGMA_PATH)" ""
SOURCES  += bsofi_gpu.c 
endif

INCLUDES += -I$(TOP_DIR) -I$(TOP_DIR)/include -I$(TOP_DIR)/config/include
CFLAGS   += $(LAPACK_CXXFLAGS) \
	$(CUDA_CXXFLAGS) $(MAGMA_CXXFLAGS) # $(DBG_CXXFLAGS) 

LIBS     +=  -limf -lirc  $(LAPACK_LIBFLAGS) $(CXXLIBS) $(CUDA_LIBFLAGS) $(MAGMA_LIBFLAGS)

HEADERS   = $(wildcard $(TOP_DIR)/include/*.h)
HEADERSXX = $(HEADERS) $(wildcard $(TOP_DIR)/include/*.hpp)

#include ../config/makefile-general-rules.mk

OBJ     ?= $(patsubst %.f90, %.o, \
	$(patsubst %.f, %.o, \
	$(patsubst %.c, %.o, $(SOURCES:.cpp=.o))))
LIBS    += $(if $(filter %.f, $(SOURCES)), $(FLIBS),)
LINKER  ?= $(CXX)
.PHONY  : clean install all lib #$(TARGET) 
ifneq   "$(TARGET)" ""
all: $(TARGET)
else
all:;@echo 'Nothing to make!'
endif

$(TARGET): $(OBJ)
	$(ARCH) $(ARFLAG) $(TARGET) $(OBJ)
	$(RANLIB) $(TARGET)

lib: 	$(TARGET)

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
	$(RM) $(TARGET) 
	$(RM) $(OBJ) 
	echo "$(MAGMA_PATH) $(MAGMA_CFLAGS)"

install: $(TARGET)
	$(MKDIR) $(INSTALL_DIR)
	$(CP) $(TARGET) $(INSTALL_DIR)/$(TARGET)

