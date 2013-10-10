/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Convenient wraps for LAPACK routines
 */

#ifndef __BSOFI_LAPACK_H__
#define __BSOFI_LAPACK_H__

void lapackXgemm(const char transa, const char transb, 
		 int m, int n, int k, 
		 scalar_t alpha, scalar_t a[], int lda, 
		 scalar_t b[], int ldb, scalar_t beta, 
		 scalar_t c[], int ldc);
void lapackXtrmm(const char side, const char uplo, 
		 const char transa, const char diag,
		 int m, int n, 
		 scalar_t alpha, scalar_t A[], int lda,
		 scalar_t B[], int ldb);
int lapackXgeqp3(int m, int n, scalar_t A[], int lda, 
		 int jpvt[], scalar_t tau[], 
		 scalar_t work[], int lwork, int* info);
int lapackXgeqrf(int m, int n, 
		 scalar_t A[], int lda, scalar_t tau[], 
		 scalar_t work[], int lwork, int* info);
int lapackXormqr(const char side, const char trans,
		 int m, int n, int k,
		 scalar_t a[], int lda,
		 scalar_t tau[], scalar_t c[], int ldc,
		 scalar_t work[], int lwork, int *info);
int lapackXorgqr(int m, int n, int k,
		 scalar_t a[], int lda, scalar_t tau[],
		 scalar_t work[], int lwork, int *info);
void lapackXtrsm(const char side, const char uplo, 
		 const char transa, const char diag,
		 int m, int n, 
		 scalar_t alpha, scalar_t A[], int lda,
		 scalar_t B[], int ldb);
int lapackXtrtri(const char uplo, const char diag, int n, 
		 scalar_t A[], int lda, int *info);
void lapackXlacpy(const char uplo, int m, int n,
		  scalar_t A[], int lda, scalar_t B[], int ldb);
int lapackXgetrf(int m, int n, scalar_t A[], int lda,
		 int ipiv[], int *info);
int lapackXgetri(int n, scalar_t A[], int lda, int ipiv[],
		 scalar_t work[], int lwork, int *info);
scalar_t lapackXnrm2(int n, scalar_t x[], int incx);
scalar_t lapackXlange(char norm, int m, int n, 
		      scalar_t A[], int lda, scalar_t work[]);
void lapackXcpyzero(char uplo, int m, int n,
		    scalar_t A[], int lda, scalar_t B[], int ldb);
void lapackXzero(char uplo, int m, int n, scalar_t A[], int lda);

#define lapackXlaset lapackXlacpy
#define lapackXlaget lapackXlacpy

#endif
