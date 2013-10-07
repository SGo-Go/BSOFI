/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Convenient wraps for LAPACK routines
 */

#ifndef __LAPACK_HPP__
#define __LAPACK_HPP__

//#include <config.h>
#ifdef __LAPACK_HPP__
#include "third_party/lapack.h"
#else
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#ifdef __SINGLE_PREC__
#  define LAROUTINE(f77name, F77NAME) FROUTINE(CONCAT(s, f77name), CONCAT(S, F77NAME))
#else
#  define LAROUTINE(f77name, F77NAME) FROUTINE(CONCAT(d, f77name), CONCAT(D, F77NAME))
#endif

#define ROUTINE_REPLACE(OLDNAME, NEWNAME, RET_TYPE, ...)	\
  extern "C" {							\
    RET_TYPE OLDNAME(__VA_ARGS__);				\
    RET_TYPE NEWNAME(__VA_ARGS__);				\
  }								\
  RET_TYPE NEWNAME(__VA_ARGS__)			

#define LAPACK_ROUTINE_REPLACE(NAME, RET_TYPE, ...)	\
  extern "C" {						\
    RET_TYPE PRROUTINE(NAME)(__VA_ARGS__);		\
    RET_TYPE LAROUTINE(NAME)(__VA_ARGS__);		\
  }							\
  RET_TYPE PRROUTINE(NAME)(__VA_ARGS__)			

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#define xGEMM  LAROUTINE(gemm,  GEMM)
#define xTRMM  LAROUTINE(trmm,  TRMM)

#define xGEQRF LAROUTINE(geqrf, GEQRF)
#define xORMQR LAROUTINE(ormqr, ORMQR)
#define xORGQR LAROUTINE(orgqr, ORGQR)
#define xGEQP3 LAROUTINE(geqp3, GEQP3)

#define xTRSM  LAROUTINE(trsm,  TRSM)
#define xLACPY LAROUTINE(lacpy, LACPY)
#define xGETRF LAROUTINE(getrf, GETRF)
#define xGETRI LAROUTINE(getri, GETRI)
#define xTRTRI LAROUTINE(trtri, TRTRI)

#define xNRM2  LAROUTINE(nrm2,  NRM2)
#define xLANGE LAROUTINE(lange, LANGE)

#ifdef __cplusplus
extern "C" {
#endif
  scalar_t xNRM2(int *n, scalar_t x[], int *incx);

  scalar_t xLANGE(const char *norm, int *m, int *n, 
		  scalar_t A[], int *lda, scalar_t work[]);

  void xGEQP3(const int* m, const int* n, 
	      scalar_t A[], const int* lda, 
	      int jpvt[], scalar_t tau[], 
	      scalar_t work[], const int* lwork, 
	      int* info);

  void xGEMM(const char *transa, const char *transb, 
	     const int *m, const int *n, const int *k, 
	     const scalar_t *alpha, scalar_t a[], const int *lda, 
	     scalar_t b[], const int *ldb, const scalar_t *beta, 
	     scalar_t c[], const int *ldc);

  void xTRMM(const char* side, const char* uplo, 
	     const char* transa, const char* diag,
	     const int *m, const int *n,
	     const scalar_t *alpha, scalar_t A[], const int *lda, 
	     scalar_t B[], const int *ldb);

  void xGEQRF(const int* m, const int* n, 
	      scalar_t A[], const int* lda, 
	      scalar_t tau[], 
	      scalar_t work[], const int* lwork, 
	      int* info);

  void xORMQR(const char* side, const char* trans,
	      const int *m, const int *n, const int *k,
	      scalar_t a[], const int *lda,
	      scalar_t tau[], scalar_t c[], const int *ldc,
	      scalar_t work[], const int *lwork, int *info);

  void xORGQR(const int *m, const int *n, const int *k,
	      scalar_t a[], const int *lda, scalar_t tau[],
	      scalar_t work[], const int *lwork, int *info);

  void xTRSM(const char* side, const char* uplo, 
	     const char* transa, const char* diag,
	     const int *m, const int *n,
	     const scalar_t *alpha,
	     scalar_t A[], const int *lda,
	     scalar_t B[], const int *ldb);

  void xTRTRI(const char* uplo, const char* diag, 
	      int* n, scalar_t A[], int* lda, int *info);

  void xLACPY(const char* uplo,
	      const int *m, const int *n,
	      scalar_t A[], const int *lda,
	      scalar_t B[], const int *ldb);

  void xGETRF(const int *m, const int *n,
	      scalar_t A[], const int *lda,
	      int ipiv[], int *info);

  void xGETRI(const int *n,
	      scalar_t A[], const int *lda,
	      int ipiv[],
	      scalar_t work[], const int *lwork,
	      int *info);

#ifdef __cplusplus
}
#endif

inline
void lapackXgemm(const char transa, const char transb, 
		 int m, int n, int k, 
		 scalar_t alpha, scalar_t a[], int lda, 
		 scalar_t b[], int ldb, scalar_t beta, 
		 scalar_t c[], int ldc) {
  const char strTransa[2] = {transa, 0};
  const char strTransb[2] = {transb, 0};
  xGEMM(strTransa, strTransb, &m, &n, &k, &alpha, a, &lda, 
	b, &ldb, &beta, c, &ldc);
}


inline 
void lapackXtrmm(const char side, const char uplo, 
		 const char transa, const char diag,
		 int m, int n, 
		 scalar_t alpha, scalar_t A[], int lda,
		 scalar_t B[], int ldb) {
  const char strSide[2]   = {side, 0};
  const char strUplo[2]   = {uplo, 0};
  const char strTransa[2] = {transa, 0};
  const char strDiag[2]   = {diag, 0};
  xTRMM(strSide, strUplo, strTransa, strDiag,
	&m, &n, &alpha, A, &lda, B, &ldb);
}



inline 
int lapackXgeqp3(int m, int n, scalar_t A[], int lda, 
		 int jpvt[], scalar_t tau[], 
		 scalar_t work[], int lwork, int* info) {
  //int info;
  xGEQP3(&m, &n, A, &lda, jpvt, tau, work, &lwork, info);
  return *info;
}

inline 
int lapackXgeqrf(int m, int n, 
		 scalar_t A[], int lda, scalar_t tau[], 
		 scalar_t work[], int lwork, int* info) {
  xGEQRF(&m, &n, A, &lda, tau, work, &lwork, info);
  return *info;
}

inline 
int lapackXormqr(const char side, const char trans,
		 int m, int n, int k,
		 scalar_t a[], int lda,
		 scalar_t tau[], scalar_t c[], int ldc,
		 scalar_t work[], int lwork, int *info) {
  const char strSide[2]   = {side, 0};
  const char strTrans[2]  = {trans, 0};
  xORMQR(strSide, strTrans, &m, &n, &k, a, &lda, 
	 tau, c, &ldc, work, &lwork, info);
  return *info;
}

inline 
int lapackXorgqr(int m, int n, int k,
		 scalar_t a[], int lda, scalar_t tau[],
		 scalar_t work[], int lwork, int *info) {
  xORGQR(&m, &n, &k, a, &lda, tau, work, &lwork, info);
  return *info;
}

inline 
void lapackXtrsm(const char side, const char uplo, 
		 const char transa, const char diag,
		 int m, int n, 
		 scalar_t alpha, scalar_t A[], int lda,
		 scalar_t B[], int ldb) {
  const char strSide[2]   = {side, 0};
  const char strUplo[2]   = {uplo, 0};
  const char strTransa[2] = {transa, 0};
  const char strDiag[2]   = {diag, 0};
  xTRSM(strSide, strUplo, strTransa, strDiag,
	&m, &n, &alpha, A, &lda, B, &ldb);
}

inline 
int lapackXtrtri(const char uplo, const char diag, int n, 
		  scalar_t A[], int lda, int *info) {
  const char strUplo[2]   = {uplo, 0};
  const char strDiag[2]   = {diag, 0};
  xTRTRI(strUplo, strDiag, &n, A, &lda, info);
  return *info;
}

inline 
void lapackXlacpy(const char uplo, int m, int n,
		  scalar_t A[], int lda, scalar_t B[], int ldb) {
  const char strUplo[2]   = {uplo, 0};
  xLACPY(strUplo, &m, &n, A, &lda, B, &ldb);
}

inline 
int lapackXgetrf(int m, int n, scalar_t A[], int lda,
		 int ipiv[], int *info) {
  xGETRF(&m, &n, A, &lda, ipiv, info);
  return *info;
}

inline 
int lapackXgetri(int n, scalar_t A[], int lda, int ipiv[],
		 scalar_t work[], int lwork, int *info) {
  xGETRI(&n, A, &lda, ipiv, work, &lwork, info);
  return *info;
}

inline 
scalar_t lapackXnrm2(int n, scalar_t x[], int incx) {
  return xNRM2(&n, x, &incx);
}

inline 
scalar_t lapackXlange(char norm, int m, int n, 
		      scalar_t A[], int lda, scalar_t work[]) {
  const char strNorm[2]   = {norm, 0};
  return xLANGE(strNorm, &m, &n, A, &lda, work);
}

#undef  A_IDX
#undef  B_IDX
#define A_IDX(__i, __j)   A[ lda*(__j) + (__i)]
#define B_IDX(__i, __j)   B[ ldb*(__j) + (__i)]

//#include <string.h>

inline 
void lapackXcpyzero(char uplo, int m, int n,
		    scalar_t A[], int lda, scalar_t B[], int ldb) {
  // const char strUplo[2]   = {'A', 0};
  // xLACPY(strUplo, &m, &n, A, &lda, B, &ldb);
  int j, i; 
  switch(uplo) {
  case 'L':
    for (j = 0; j < n; j++)
      for (i = j + 1; i < m; i++)
	{B_IDX(i, j) = A_IDX(i, j); A_IDX(i, j) = 0;}
    break;
  case 'U':
    for (j = 0; j < n; j++)
      for (i = 0; i <= j; i++)
	{B_IDX(i, j) = A_IDX(i, j); A_IDX(i, j) = 0;}
    break;
  default:
    for (j = 0; j < n; j++)
      for (i = 0; i < m; i++) 
	{B_IDX(i, j) = A_IDX(i, j); A_IDX(i, j) = 0;}
  }
}

inline 
void lapackXzero(char uplo, int m, int n, scalar_t A[], int lda) {
  int j, i; 
  switch(uplo) {
  case 'L':
    for (j = 0; j < n; j++)
      for (i = j; i < m; i++)
	A_IDX(i, j) = 0;
    break;
  case 'U':
    for (j = 0; j < n; j++)
      for (i = 0; i <= j; i++)
	A_IDX(i, j) = 0;
    break;
  default:
    for (j = 0; j < n; j++)
      for (i = 0; i < m; i++) {
	A_IDX(i, j) = 0;
      }
  }
}
#undef A_IDX
#undef B_IDX

#endif
#endif
