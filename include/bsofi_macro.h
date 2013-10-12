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

#ifdef HAS_CUBLAS
/* @TODO Fix issue with `if' in macro */
/* #  define CHECK_CUMALLOC(__code) __code */
#  define CHECK_CUMALLOC(__code)					\
  if (cudaSuccess != (__code)) {					\
    DBGERROR("CUDA: GPU device memory allocation failed");		\
    /* cudaFree(dwork); */						\
    cublasShutdown();/* cublasDestroy(handle); */			\
    return -1;								\
  }

#  define CUBLAS_INIT()						\
  if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {			\
    DBGERROR("CUBLAS initialization failed");			\
    exit(-1);							\
  }

#  define CUBLAS_FINALIZE()			\
  cublasShutdown() /* cublasDestroy(handle); */

#else
#  define CUBLAS_INIT()
#  define CUBLAS_FINALIZE()
#  define CHECK_CUMALLOC(__code) __code
#endif

#endif