/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Header for GPU implementation of routines based on Magma
 *   and CuBLAS subroutines calls
 */

#ifndef __BSOFI_GPU_H__
#define __BSOFI_GPU_H__

/* Settings */
#include <config.h>
#define gpuXbstri gpuXbstri_rseq //rbatched
#define gpuXbsoi  gpuXbsoi_batched // orgqr ormqr 

/* Callback names definitions for Fortran */
#define cb_gpuXbsofiLWork FROUTINE(gpuxbsofilwork, GPUXBSOFILWORK)
#define cb_gpuXbsof       FROUTINE(gpuxbsof, GPUXBSOF)
#define cb_gpuXbstri      FROUTINE(gpuxbstri, GPUXBSTRI)
#define cb_gpuXbsoi       FROUTINE(gpuxbsoi, GPUXBSOI)

/* Names with plain C mangling */
#ifdef __cplusplus
extern "C" {
#endif
  /* int gpuXbsoi(int n, int L, scalar_t A[], int lda,  */
  /* 	     scalar_t _tau[], scalar_t *dT,  */
  /* 	     scalar_t work[], int lwork, int* info); */

  int gpuXbsof(int n, int L, 
	       scalar_t dA[], int ldda, scalar_t _tau[], 
	       scalar_t work[], int lwork, scalar_t* dT, int* info);

  int gpuXbstri_rseq(int n, int L, scalar_t dA[], int ldda, 
		     scalar_t dwork[], int ldwork, 
		     int* info);

  int gpuXbstri_rbatched(int n, int L, scalar_t dA[], int ldda, 
			 scalar_t dwork[], int ldwork, 
			 int* info);

  int gpuXbsoi_ormqr(int n, int L, 
		     scalar_t dA[], int ldda, scalar_t tau[], 
		     scalar_t work[], int lwork, 
		     scalar_t *dT, scalar_t dwork[], int ldwork, 
		     int* info);

  int gpuXbsoi_orgqr(int n, int L, 
		     scalar_t dA[], int ldda, scalar_t _tau[], 
		     scalar_t work[], int lwork, 
		     scalar_t *_dT, scalar_t dwork[], int ldwork, 
		     int* info);

  int gpuXbsoi_batched(int n, int L, 
		       scalar_t dA[], int ldda, scalar_t _tau[], 
		       scalar_t work[], int lwork, 
		       scalar_t *_dT, scalar_t dwork[], int ldwork, 
		       int* info);

  /* Callback functions declaration for Fortran calls */
  int cb_gpuXbsofiLWork(int* n, int* L, int* lda);
  
  int cb_gpuXbsof(int* n, int* L, 
		  scalar_t A[], int* lda, scalar_t _tau[], 
		  scalar_t work[], int* lwork, int* info);
  
  int cb_gpuXbsoi_orgqr(int* n, int* L, 
			scalar_t A[], int* lda, scalar_t _tau[], 
			scalar_t work[], int* lwork, int* info);

  int cb_gpuXbsoi_ormqr(int* n, int* L, 
		     scalar_t A[], int* lda, scalar_t _tau[], 
		     scalar_t work[], int* lwork, int* info);

  int cb_gpuXbstri(int* n, int* L, 
		   scalar_t A[], int* lda, 
		   // scalar_t work[], int lwork, 
		int* info);

#ifdef __cplusplus
    }
#endif


/************************************************************
 * Auxiliary GPU routines
 ************************************************************/
int get_QR_T_device_mem(int m, int n, int k);

void gpuXcpyzero_gpu(char uplo, int m, int n, 
		     scalar_t *A, int ldda, scalar_t *B, int lddb);

void gpuXzero_gpu(char uplo, int m, int n, 
		  scalar_t *dA, int ldda);

#include "magma.h"
#define gpuXgemm_gpu   gpuDgemm
#define gpuXtrmm_gpu   gpuDtrmm
#define gpuXtrsm_gpu   magmablas_dtrsm
#define gpuXlacpy_gpu  magmablas_dlacpy
#define gpuXtrtri_gpu  magma_dtrtri_gpu
#define gpuXgeqrf_gpu  magma_dgeqrf_gpu
#define gpuXormqr_gpu  magma_dormqr_gpu
#define gpuXorgqr_gpu  magma_dorgqr_gpu

#define gpuXgemm_batched cublasDgemmBatched
#define gpuXtrsm_batched cublasDtrsmBatched
#define gpuXtrmm_batched cublasDgemmBatched

#endif
