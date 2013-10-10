# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Auto-selection of configuration.

CUDA_PATH      = #/usr/local/cuda
BOOST_PATH     = /usr#/home/sergiy/opt/boost_1_52_0
MKL_PATH       = #/opt/intel/mkl/lib/intel64
MAGMA_PATH     = #/usr/local#/home/aetdom/magma-1.2.1
ATLAS_PATH     = /home/sergiy/ucd/la/ATLAS3.10.0/BLDdir
LAPACK_PATH    = /home/sergiy/ucd/la/lapack-3.4.2

LAPACK         = atlas

CC   = gcc
CXX  = g++
FC   = gfortran
LINK = gcc #g++

