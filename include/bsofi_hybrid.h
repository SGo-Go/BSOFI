/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Header for CPU+GPU implementation of routines based on LAPACK
 *   and CuBLAS subroutines calls
 */

#ifndef __BSOFI_HYBRID_H__
#define __BSOFI_HYBRID_H__

/* Settings */
#include <config.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include "bsofi_macro.h"

/* Names with plain C mangling */
#ifdef __cplusplus
extern "C" {
#endif

  int hybridXbsofiLWork(int n, int L, int lda, int lworkHostDevice);

  /************************************************************
   * Hybrid solvers
   ************************************************************/
#include "bsofi_cpu.h"

#define hybridXbsof  lapackXbsof
#define hybridXbstri lapackXbstri

  int hybridXbsoftri(cublasHandle_t handle, int n, int L, 
		     scalar_t A[], int lda, scalar_t _tau[], 
		     scalar_t work[], int lwork, scalar_t dwork[], int ldwork,
		     /* int Lswitch,  */
		     int* info);

  int hybridXbsoi(cublasHandle_t handle, int n, int L, 
		  scalar_t A[], int lda, scalar_t _tau[], 
		  scalar_t work[], int lwork, scalar_t dwork[], int ldwork, 
		  /* int Lswitch,  */
		  int* info);

#ifdef USE_FORTRAN
  /* Callback names definitions for Fortran */
#  define cb_hybridXbsofiLWork  FROUTINE(hybridxbsofilwork, HYBRIDXBSOFILWORK)
#  define cb_hybridXbsoftri     FROUTINE(hybridxbsoftri, HYBRIDXBSOFTRI)
#  define cb_hybridXbsoi        FROUTINE(hybridxbsoi, HYBRIDXBSOI)

  /* Callback functions declaration for Fortran calls */
  int cb_hybridXbsofiLWork(int* n, int* L, int* lda, int* lworkHostDevice);

  int cb_hybridXbsoftri(/* cublasHandle_t handle, */ int* n, int* L, 
			scalar_t A[], int* lda, scalar_t tau[], 
			scalar_t work[], int* lwork, scalar_t dwork[], int* ldwork,
			/* int Lswitch,  */
			int* info);

  int cb_hybridXbsoi(/* cublasHandle_t handle,  */ int* n, int* L, 
		     scalar_t A[], int* lda, scalar_t tau[], 
		     scalar_t work[], int* lwork, scalar_t dwork[], int* ldwork, 
		     /* int Lswitch,  */
		     int* info);


  /* Extra-callback names definitions for Fortran */
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

#ifdef USE_PROF
  extern bsofi_profile_t profHybridBSOI;
  extern bsofi_profile_t profHybridBSOFTRI;
#endif

#ifdef __cplusplus
    }
#endif

#endif
