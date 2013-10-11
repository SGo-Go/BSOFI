/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Convenient wraps for CUBLAS routines.
 *   To be able to compile our codes on platforms without CUDA,
 *   we substitute CUBLAS calls by LAPACK calls if needed.
 */

#ifndef __BSOFI_CUBLAS_H__
#define __BSOFI_CUBLAS_H__

#include <config.h>

#ifdef HAS_CUBLAS
#  include <cuda_runtime.h>
#  include <cublas.h>

/* @TODO Fix issue with `if' in macro */

/* #  define HANDLE_CUBLAS_ERROR(__code, __message) __code */
#  define HANDLE_CUBLAS_ERROR(__code, __message)			\
  if (__code != CUBLAS_STATUS_SUCCESS) {				\
    DBGERROR("CUBLAS error in " __message);				\
    cudaFree(dwork); cublasShutdown();/* cublasDestroy(handle); */	\
    *info = -1;								\
    return -1;								\
  }


#  define cublasXlaset(__t, __m, __n, __A, __lda, __B, __ldb)		\
  HANDLE_CUBLAS_ERROR( cublasSetMatrix					\
		       (__m, __n, sizeof(scalar_t),			\
			__A, __lda, __B, __ldb), "cublasSetMatrix")

#  define cublasXlaget(__t, __m, __n, __A, __lda, __B, __ldb)		\
  HANDLE_CUBLAS_ERROR( cublasGetMatrix					\
		       (__m, __n, sizeof(scalar_t),			\
			__A, __lda, __B, __ldb), "cublasGetMatrix")

#  ifdef __SINGLE_PREC__
#    define cublasXgemm  cublasDgemm
#  else
#    define cublasXgemm  cublasSgemm
#  endif

#else
#  include <third_party/lapack.hpp>

#  define cublasXgemm  lapackXgemm
#  define cublasXlaget lapackXlaget
#  define cublasXlaset lapackXlaset
#  define cudaDeviceSynchronize() 
#endif

#endif
