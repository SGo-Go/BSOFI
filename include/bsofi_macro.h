/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Common macro related to implementations of BSOF/I subroutines
 */

#ifndef __BSOFI_MACRO_H__
#define __BSOFI_MACRO_H__

#include <config.h>
#include "third_party/cublas.h"

#define GPU_ALIGN 64 /* 32 */
#define GET_LD_DEV(__len) \
  ((((__len) + (GPU_ALIGN)-1)/(GPU_ALIGN))*(GPU_ALIGN))

#define IDX_OFFSET(__ld,__i, __j) ((__j)*(__ld) + __i)
#define BLK_OFFSET(__ld,__i, __j) (n*((__j)*(__ld) + __i))

#ifndef A_BLK
#define A_IDX_OFFSET(__i, __j) ((__j)*lda + __i)
#define A_BLK_OFFSET(__i, __j) (n*((__j)*lda + __i))

#define A_IDX(__i, __j) A[A_IDX_OFFSET(__i, __j)]
#define A_BLK(__i, __j) A + A_BLK_OFFSET(__i, __j)
#endif

#include <stdio.h>
#define CHECK_PTR(__ptr)						\
  if(!(__ptr)) {							\
    /* DBGERROR("Invalid pointed" #__ptr); */				\
    exit(-1);								\
  }


typedef struct _bsofi_profile_t {
  double cpu;
  double gpu;
  double memcpu;
  double memset;
  double memget;
  double sync;
} bsofi_profile_t;

#ifdef USE_PROF

#include <timer.h>
#include <stdio.h>

#  define BENCH_CUMMULATIVE(__total_time, __code_do) {	\
    tim1 = getwalltime();				\
    {__code_do;}					\
    __total_time += elapsed(getwalltime(), tim1); }
#  define RESET_BSOFI_PROFILE(__counter) \
  __counter = (const struct _bsofi_profile_t){0}

#else
#  define BENCH_CUMMULATIVE(__total_time, __code_do) {__code_do;}
#  define RESET_BSOFI_PROFILE(__counter) 
#endif

/************************************************************
 *  CUDA initialization/finalization macro
 ************************************************************/

#ifdef HAS_CUBLAS
/* #  define CHECK_CUMALLOC(__code) __code */
#  define CHECK_CUMALLOC(__code)					\
  if (cudaSuccess != (__code)) {					\
    DBGERROR("CUDA: GPU device memory allocation failed");		\
    /* cudaFree(dwork); */						\
    cublasShutdown();/* cublasDestroy(handle); */			\
    return -1;								\
  }

#  define CUBLAS_INIT(__handle)					\
  if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {			\
    DBGERROR("CUBLAS initialization failed");			\
    exit(-1);							\
  }

#  define CUBLAS_FINALIZE(__handle)		\
  cublasShutdown() /* cublasDestroy(handle); */

#else
#  define CUBLAS_INIT(__handle)
#  define CUBLAS_FINALIZE(__handle)
#  define CHECK_CUMALLOC(__code) __code
#endif

/************************************************************
 *  Workspace allocation/deallocation macro for CUBLAS codes
 ************************************************************/
#ifdef HAS_CUBLAS
#  include <cuda_runtime.h>
#  include <cublas.h>

#  define CUBLAS_HOSTALLOC(__ptr, __type, __size)			\
  if ( cudaSuccess !=							\
       cudaMallocHost( (void**) &__ptr, (__size)*sizeof(__type) )) {	\
    DBGERROR("CUDA pinned malloc failed for: %s", #__ptr );	\
    exit(-1);							\
  }

#  define CUBLAS_HOSTFREE(__ptr)			\
  cudaFreeHost( ptr );

#  define CUBLAS_DEVALLOC(__ptr, __type, __size)			\
  if ( cudaSuccess !=							\
       cudaMalloc( (void**) &__ptr, (__size)*sizeof(__type) )) {	\
    DBGERROR("CUDA device malloc failed for: %s", #__ptr );		\
    exit(-1);								\
  }

#  define CUBLAS_DEVFREE(__ptr)			\
  cudaFree( __ptr );

#else
/* #  define CUBLAS_HOSTALLOC(__ptr, __type, __size) LAPACK_MALLOC( __ptr, __type, __size ) */
/* #  define CUBLAS_HOSTFREE (__ptr)             LAPACK_FREE(__ptr) */
/* #  define CUBLAS_DEVALLOC (__ptr, __type, __size) LAPACK_MALLOC( __ptr, __type, __size ) */
/* #  define CUBLAS_DEVFREE  (__ptr)             LAPACK_FREE(__ptr) */
#  define CUBLAS_HOSTALLOC LAPACK_MALLOC
#  define CUBLAS_HOSTFREE  LAPACK_FREE
#  define CUBLAS_DEVALLOC  LAPACK_MALLOC
#  define CUBLAS_DEVFREE   LAPACK_FREE
#endif 

/************************************************************
 *  Workspace allocation/deallocation macro for LAPACK codes
 ************************************************************/
#include <stdlib.h>
#define LAPACK_MALLOC(__ptr, __type, __size)			\
  if ( 0 == (__ptr = malloc((__size)*sizeof(__type)))) {	\
    DBGERROR("malloc failed for: %s", #__ptr );			\
    exit(-1);							\
  }

#define LAPACK_FREE(__ptr)			\
  free(__ptr);

#endif
