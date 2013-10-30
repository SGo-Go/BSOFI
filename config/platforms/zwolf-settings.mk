# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Auto-selection of configuration.

CUDA_PATH      = /usr/local/cuda
BOOST_PATH     = /home/sergiy/opt/boost_1_52_0
MKL_PATH       = /opt/intel/mkl/lib/intel64
MAGMA_PATH     = /home/sergiy/opt/GPU/magma-1.3.0
#/home/aetdom/magma-1.2.1
QUEST_PATH     = /home/sergiy/opt/QUEST/last


ARPACK_PATH    = 

CC   = icc -DBSOFI_PLATFORM=1
CXX  = icc -DBSOFI_PLATFORM=1 #icpc
FC   = gfortran
