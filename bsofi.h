/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *  Host header file
 */

#ifndef __BSOFI_H__
#define __BSOFI_H__

#include "bsofi_flops.h"

#ifdef __cplusplus
extern "C" {
#endif
  
#ifdef HAS_CUDA
#include "bsofi_hybrid.h"
#endif
#include "bsofi_cpu.h"

#ifdef HAS_MAGMA
#include "bsofi_gpu.h"
#endif

#ifdef __cplusplus
}
#endif

#undef  lapackXbsoi
#define lapackXbsoi lapackXbsoi_ormqr //orgqr //

#ifdef USE_FORTRAN
/* Callback names definitions for use of switches 
   to the best of CPU and CPU+GPU codes in Fortran */
/* #include <cuda_runtime.h> */

typedef struct __bsofiHandle_t {
#ifdef HAS_CUDA
    cublasHandle_t handle;
#else
    int handle;
#endif
    scalar_t *dwork;
    int ldwork;
    int switchBsoftri;
    int switchBsoi;
  } bsofiHandle_t;

#  define cb_hybridXinitH       FROUTINE(hybridxinith, HYBRIDXINITH)
#  define cb_hybridXfinalizeH   FROUTINE(hybridxfinalizeh, HYBRIDXFINALIZEH)
#  define cb_hybridXbsofiLWorkH FROUTINE(hybridxbsofilworkh, HYBRIDXBSOFILWORKH)
#  define cb_hybridXbsoftriH    FROUTINE(hybridxbsoftrih, HYBRIDXBSOFTRIH)
#  define cb_hybridXbsoiH       FROUTINE(hybridxbsoih, HYBRIDXBSOIH)

int cb_hybridXinitH(/* bsofiHandle_t* handle*/int* _handle, int* n, int* L, int* lda);
int cb_hybridXbsoftriH(/* bsofiHandle_t* handle*/int* _handle, 
		       int* n, int* L, scalar_t A[], int* lda, 
		       scalar_t tau[], scalar_t work[], int* lwork, int* info);
int cb_hybridXbsoiH(/* bsofiHandle_t* handle*/int* _handle, 
		    int* n, int* L, scalar_t A[], int* lda, 
		    scalar_t tau[], scalar_t work[], int* lwork, int* info);
int cb_hybridXfinalizeH(/* bsofiHandle_t* handle*/int* _handle);
int cb_hybridXbsofiLWorkH(/* bsofiHandle_t* handle*/int* _handle, int* n, int* L, int* lda);

#endif 

#endif
