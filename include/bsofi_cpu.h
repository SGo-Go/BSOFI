/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Header for CPU implementation of routines based on LAPACK
 *   subroutines calls
 */

#ifndef __BSOFI_CPU_H__
#define __BSOFI_CPU_H__

/* Settings */
#include <config.h>
#define lapackXbsoi lapackXbsoi_orgqr

/* Names with plain C mangling */
#ifdef __cplusplus
extern "C" {
#endif

  int lapackXbsofiLWork(int n, int L, int lda);

  int lapackXbsof(int n, int L, 
		  scalar_t A[], int lda, scalar_t _tau[], 
		  scalar_t work[], int lwork, int* info);

  int lapackXbsoi_orgqr(int n, int L, 
			scalar_t A[], int lda, scalar_t _tau[], 
			scalar_t work[], int lwork, int* info);

  int lapackXbsoi_ormqr(int n, int L, 
			scalar_t A[], int lda, scalar_t _tau[], 
			scalar_t work[], int lwork, int* info);

  int lapackXbstri(int n, int L, 
		   scalar_t A[], int lda, 
		   // scalar_t work[], int lwork, 
		   int* info);



#ifdef USE_FORTRAN

/* Callback names definitions for Fortran */
#  define cb_cpuXbsofiLWork FROUTINE(lapackxbsofilwork, LAPACKXBSOFILWORK)
#  define cb_cpuXbsof       FROUTINE(lapackxbsof, LAPACKXBSOF)
#  define cb_cpuXbstri      FROUTINE(lapackxbstri, LAPACKXBSTRI)
#  define cb_cpuXbsoi_ormqr FROUTINE(lapackxbsoi_ormqr, LAPACKXBSOI_ORMQR)
#  define cb_cpuXbsoi_orgqr FROUTINE(lapackxbsoi_orgqr, LAPACKXBSOI_ORGQR)

  /* Callback functions declaration for Fortran calls */
  int cb_cpuXbsofiLWork(int* n, int* L, int* lda);

  int cb_cpuXbsof(int* n, int* L, 
		  scalar_t A[], int* lda, scalar_t _tau[], 
		  scalar_t work[], int* lwork, int* info);

  int cb_cpuXbsoi_orgqr(int* n, int* L, 
			scalar_t A[], int* lda, scalar_t _tau[], 
			scalar_t work[], int* lwork, int* info);

  int cb_cpuXbsoi_ormqr(int* n, int* L, 
			scalar_t A[], int* lda, scalar_t _tau[], 
			scalar_t work[], int* lwork, int* info);

  int cb_cpuXbstri(int* n, int* L, 
		   scalar_t A[], int* lda, 
		   // scalar_t work[], int lwork, 
		   int* info);
#endif

#ifdef __cplusplus
}
#endif
#endif
