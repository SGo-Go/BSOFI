# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Script with makefile default settings.
# ------------------------------------------------------------------

BSOFI_LIB := $(TOP_DIR)/bsofi.la

INSTALL_DIR ?= $(TOP_DIR)/bin

CC   ?= icc 
CXX  ?= icpc
FC   ?= ifort
NVCC ?= nvcc

MAKE   = make
ARCH   = ar cr
RANLIB = ranlib
CP     = cp 
CD     = cd 
MKDIR  = mkdir -p 
RM     = rm -f 


FLIBS          = $(if $(filter $(FC),  gfortran), -lgfortran, -lg2c)
CLIBS          = -lc
CXXLIBS        = -lstdc++ 

OMP_CXXFLAGS   = $(if $(filter $(CXX), g++), -fopenmp, -openmp)
OMP_CFLAGS     = $(if $(filter $(CC),  gcc g++), -fopenmp, -openmp)
OMP_FFLAGS     = $(if $(filter $(FC),  gfortran g77), -fopenmp, -openmp)
#-fopenmp
OMP_LIBFLAGS   = $(if $(filter $(CXX), g++), -lrt -lgomp, \
	-openmp)

DBG_CXXFLAGS   = -g -DDEBUG
DBG_CFLAGS     = -g -DDEBUG
DBG_FFLAGS     = -g -DDEBUG


VERBOSE_CXXFLAGS = $(if $(filter $(CXX), g++), , \
	-Wall) #-ftree-vectorizer-verbose=2/-vec-report3
VERBOSE_CFLAGS   = -Wall #\
	$(if $(filter  $(CC), g++ gcc), -ftree-vectorizer-verbose=2, -vec-report2)

VEC_CXXFLAGS   = $(if $(filter $(CXX), icc icpc), -restrict,)
VEC_CFLAGS     = $(if $(filter  $(CC), icc icpc), -restrict, \
	$(if $(filter $(CXX), gcc), -ftree-vectorize -ffast-math -msse3 -mssse3 -fbuiltin))
VEC_FFLAGS     = #-restrict

OPT_CXXFLAGS   := -O3 $(VEC_CXXFLAGS)
OPT_CFLAGS     := -O3 $(VEC_CFLAGS)
OPT_FFLAGS     := -O3 $(VEC_FFLAGS)
OPT_CUFLAGS    := -O3 -Xptxas -v -m 64 -arch compute_20
# CU_FLAGS = -g -Xptxas -v -m 64 -arch sm_20

CUDA_PATH     ?= /usr/local/cuda
ifneq "$(CUDA_PATH)" ""
CUDA_CFLAGS    = -I$(CUDA_PATH) -DHAS_CUBLAS
CUDA_CXXFLAGS  = -I$(CUDA_PATH)/include -DHAS_CUBLAS
CUDA_LIBFLAGS  = -L$(CUDA_PATH)/lib64 -lcublas -lcudart 
endif 

BOOST_PATH    ?= /usr
ifneq "$(BOOST_PATH)" ""
BOOST_INCLUDE  = $(BOOST_PATH)/include
BOOST_LIB      = $(BOOST_PATH)/lib
BOOST_CXXFLAGS = -I$(BOOST_INCLUDE) -DHAS_BOOST
BOOST_LIBFLAGS = $(BOOST_LIB)/libboost_program_options.a
endif 

LAPACK ?= mkl
include $(TOP_DIR)/config/thirdparty/$(LAPACK).mk

MAGMA_PATH    ?= /usr/local
ifneq "$(MAGMA_PATH)" ""
MAGMA_INCLUDE  = $(MAGMA_PATH)/include
MAGMA_LIB      = $(MAGMA_PATH)/lib
MAGMA_CFLAGS   = -I$(MAGMA_INCLUDE) -DHAS_MAGMA
MAGMA_CXXFLAGS = -I$(MAGMA_INCLUDE) -DHAS_MAGMA
MAGMA_LIBFLAGS = -L$(MAGMA_LIB) -lmagma -lmagmablas
endif 

QUEST_PATH    ?= /usr/local/QUEST
QUEST_LIBFLAGS = $(QUEST_PATH)/libdqmc.a
QUEST_CXXFLAGS = -I$(QUEST_PATH)/SRC -DHAS_QUEST

FFLAGS    = $(OPT_FFLAGS) $(OMP_FFLAGS)
CFLAGS    = $(VERBOSE_CFLAGS) $(OPT_CFLAGS) $(OMP_CFLAGS) $(INCLUDES)
CXXFLAGS  = $(VERBOSE_CXXFLAGS) $(OPT_CXXFLAGS) $(OMP_CXXFLAGS) $(INCLUDES)
CUFLAGS   = $(OPT_CUFLAGS)
LIBS      = $(OMP_LIBFLAGS)

INCLUDES := -I. 
THIS_DIR := .

SOURCES  :=
