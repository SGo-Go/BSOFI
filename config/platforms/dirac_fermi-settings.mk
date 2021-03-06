# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Default configuration for Dirac at NERSC.

CUDA_PATH      = /usr/common/usg/cuda/5.5
MKL_PATH       = /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117
BOOST_PATH     = /usr/common/usg/boost/1.50
MAGMA_PATH     = #/usr/common/usg/magma/1.1.0
QUEST_PATH     = 
ATLAS_PATH     = /usr/common/usg/ATLAS/3.9.15
LAPACK         = mkl

CC   = icc -DBSOFI_PLATFORM=10 #pgcc
CXX  = icc -DBSOFI_PLATFORM=10 #pgcc
FC   = #gfortran
