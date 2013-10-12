# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Auto-selection of configuration.

CUDA_PATH      = /usr/common/usg/cuda/5.0
BOOST_PATH     = #/usr/opt
MKL_PATH       = /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117
#/usr/common/usg/intel/mkl/10.2.5.035/lib/em64t
MAGMA_PATH     = #/home/sergiy/opt/GPU/magma-1.3.0
QUEST_PATH     = 
ATLAS_PATH     = /usr/common/usg/ATLAS/3.9.15
LAPACK = atlas

CC   = gcc
CXX  = gcc #icpc
FC   = gfortran
