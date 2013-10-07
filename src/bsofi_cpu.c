/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 *  Univ. of California, Davis
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Pure CPU implementation of routines based on LAPACK 
 *   subroutines calls)
 */

#include "bsofi.h"
#include <config.h>
#include <third_party/lapack.hpp>

/* ===================================================================
   BSOF on CPU
   =================================================================== */

#define A_IDX_OFFSET(__i, __j) ((__j)*lda + __i)
#define A_BLK_OFFSET(__i, __j) (n*((__j)*lda + __i))

#define A_IDX(__i, __j) A[A_IDX_OFFSET(__i, __j)]
#define A_BLK(__i, __j) A + A_BLK_OFFSET(__i, __j)

#ifdef __cplusplus
//#include "profile.hpp"
#endif

int lapackXbsofiLWork(int n, int L, int lda)
{
  int lworkTotal = 0;
  scalar_t *A = 0; scalar_t *tau = 0; 
  scalar_t work[1]; int lwork = -1; int info;
  lapackXbsof(n, L, A, lda, tau, work, lwork, &info);
  lworkTotal = work[0];
  lwork = -1; 
  lapackXbsoi(n, L, A, lda, tau, work, lwork, &info);
  if(lworkTotal < work[0]) lworkTotal = work[0];
  return lworkTotal;
}

/**
 * Purpose:
 *   {\bf lapackXbsof} computes a structured orthogonal factorization 
 *   of a $N$-by-$N$ $p$-cyclic matrix $A$ where $p = L$, $N = n*L$. 
 *   The factorization has the form
 *   $A = Q^{L-1} * ... Q^1 * R$
 *   R is upper triangular, $Q^k$ are orthogonal.
 *   This is the Level 3 BLAS version of the algorithm.
 * ********************************************************************
 * Arguments:
 * @arg n (input) INTEGER
 *   The block size of the matrix $A$ ($n > 0$).
 * @arg L (input) INTEGER
 *   The number of row/column blocks of the matrix $A$ ($L > 4$).
 * @arg A (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
 *   On entry, the $N$-by-$N$ matrix to be factored.
 *   On exit, the Householder reflectors of factors $Q^k$ 
 *   and factor $R$ from the factorization.
 *   The unit scaling coefficients of Householder reflectors are not stored.
 * @arg lda (input) INTEGER
 *   The leading dimension of the array A. $lda >= max(1,n*L)$.
 * @arg tau (output) DOUBLE_PRECISION array, dimension ($n*L$)
 *   Scaling coefficients of the elementary Householder reflectors.
 * @arg work (workspace/output) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))
 *   On exit, if $info = 0$, WORK(1) returns the optimal LWORK.
 * @arg lwork (input) INTEGER
 *   The dimension of the array WORK. 
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 * ********************************************************************
 * Further Details:
 *   The matrices $Q^k$ are represented as extensions of matrices $Q^{(k)}$
 *     $Q^k = I_{n(k-1)}\oplus Q^{(k)} \oplus I_{n(L-k-1)}$.
 *   whereas $Q^{(k)}$ are represented as products of elementary reflectors
 *     $Q^{(k)} = H^{(k)}(1) H^{(k)}(2) . . . H^{(k)}(n)$.
 *   Each $H^{(k)}(i)$ has the form
 *     $H^{(k)}(i) = I - \tau * v * v'$
 *   where $\tau$ is a real scalar, and $v$ is a real vector with
 *   $v(1:i-1) = 0$ and $v(i) = 1$; $v(i+1:m)$ is stored on exit 
 *   in A(k*n + i + 1:(k+2)*n,k*n+i) and $\tau$ in tau(k*n + i).
 */
int lapackXbsof(int n, int L, 
		scalar_t A[], int lda, scalar_t _tau[], 
		scalar_t work[], int lwork, int* info)
{
  int i;
  scalar_t *Rii, *Riip1, *Rin, *tau = _tau;
  scalar_t max_work0; 

  if(lwork == -1) {
    lapackXgeqrf(2*n, 2*n, A, lda, tau, work, lwork, info);
    max_work0 = work[0];
    lapackXormqr('L', 'T', 2*n, n, n,
		 A, lda, tau, A, lda, work, lwork, info);
    if(max_work0 > work[0]) work[0] = max_work0;
    return *info;
  }

  for(i = 0; i < L-2; i++) {
    Rii   = A_BLK(i, i); 
    Riip1 = A_BLK(i, i+1); 
    Rin   = A_BLK(i, L-1);
    lapackXgeqrf(2*n, n, Rii, lda, tau, work, lwork, info);
    lapackXormqr('L', 'T', 2*n, n, n,
    		 Rii, lda, tau, Riip1, lda, work, lwork, info);
    lapackXormqr('L', 'T', 2*n, n, n,
     		 Rii, lda, tau, Rin, lda, work, lwork, info);
    tau +=n;
  }

  Rii   = A_BLK(L-2,L-2); 
  lapackXgeqrf(2*n, 2*n, Rii, lda, tau, work, lwork, info);
  return *info;
}

/**
 * Purpose:
 *   {\bf lapackXbstri} computes inverse of factor $R$ for a factorized 
 *   $N$-by-$N$ $p$-cyclic matrix $A$ where $p = L$, $N = n*L$. 
 * ********************************************************************
 * Arguments:
 * @arg n (input) INTEGER
 *   The block size of the matrix $A$ ($n > 0$).
 * @arg L (input) INTEGER
 *   The number of row/column blocks of the matrix $A$ ($L > 4$).
 * @arg A (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
 *   On entry, matrix $A$ factorized by 
 *   structured orthogonal factorization method. Namely, 
 *   the Householder reflectors of factors $Q^k$ 
 *   and factor $R$ from the factorization.
 *   On exit, the Householder reflectors of factors $Q^k$ 
 *   and inverse of factor $R$.
 *   The unit scaling coefficients of Householder reflectors are not stored.
 * @arg lda (input) INTEGER
 *   The leading dimension of the array A. $lda >= max(1,n*L)$.
 * @arg tau (output) DOUBLE_PRECISION array, dimension ($n*L$)
 *   Scaling coefficients of the elementary Householder reflectors.
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 */
int lapackXbstri(int n, int L, 
		 scalar_t A[], int lda, 
		 //scalar_t work[], int lwork, 
		 int* info)
{
  int N = n*L, i = L - 3;
  //scalar_t *W = A_BLK(L-1, 0);

  if(L <= 2)
    lapackXtrtri('U', 'N', N, A_BLK(0, 0), lda, info);
  else {
    lapackXtrtri('U', 'N', 3*n, A_BLK(L-3, L-3), lda, info);

    if(L > 3) {
      lapackXtrmm('R', 'U', 'N', 'N', n*(L-3), n,  1, A_BLK(L-1, L-1), lda, A_BLK(0, L-1), lda);
      for(i = L - 4; i >= 0; i--) {
	lapackXtrtri('U', 'N', n, A_BLK(i, i), lda, info);
	lapackXtrmm('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, A_BLK(i, L-1), lda); 

	lapackXtrmm('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, A_BLK(i, i+1), lda);
	lapackXgemm('N', 'N', n, n*(L-i-2),n, 1, A_BLK(i, i+1), lda, A_BLK(i+1,i+2), lda, 1, A_BLK(i, i+2), lda);
	lapackXtrmm('R', 'U', 'N', 'N', n, n, 1, A_BLK(i+1, i+1), lda, A_BLK(i, i+1), lda); 

	// lapackXcpyzero('A', n, n, A_BLK(i, i+1), lda, W, lda);
	// lapackXtrmm('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, W, lda);
	// lapackXgemm('N', 'N', n, n*(L-i-1),n, 1, W, lda, A_BLK(i+1,i+1), lda, 1, A_BLK(i, i+1), lda);
      }
      // lapackXzero('A', n, n, W, lda);
    }
  }
  return *info;
}

/**
 * Purpose:
 *   {\bf lapackXbsoi_ormqr} computes inverse of a factorized 
 *   $N$-by-$N$ $p$-cyclic matrix $A$ with pre-inverted factor $R$
 *   where $p = L$, $N = n*L$. 
 *   Uses xORMQR subroutine for blocked Householder updates
 * ********************************************************************
 * Arguments:
 * @arg n (input) INTEGER
 *   The block size of the matrix $A$ ($n > 0$).
 * @arg L (input) INTEGER
 *   The number of row/column blocks of the matrix $A$ ($L > 4$).
 * @arg A (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
 *   On entry, matrix $A$ factorized by 
 *   structured orthogonal factorization method with inverted factor $R$. 
 *   Namely, 
 *   the Householder reflectors of factors $Q^k$ 
 *   and factor $R^{-1}$.
 *   The unit scaling coefficients of Householder reflectors are not stored.
 *   On exit, inverse of matrix $A$.
 * @arg lda (input) INTEGER
 *   The leading dimension of the array A. $lda >= max(1,n*L)$.
 * @arg tau (output) DOUBLE_PRECISION array, dimension ($n*L$)
 *   Scaling coefficients of the elementary Householder reflectors.
 * @arg work (workspace/output) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))
 *   On exit, if $info = 0$, WORK(1) returns the optimal LWORK.
 * @arg lwork (input) INTEGER
 *   The dimension of the array WORK. 
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 */
int lapackXbsoi_ormqr(int n, int L, 
		      scalar_t A[], int lda, scalar_t _tau[], 
		      scalar_t work[], int lwork, int* info)
{
  int N = n*L, i, ldq = 2*n;
  scalar_t *R1i=A, *Qii, *tau = _tau + (L-2)*n;
  int lworkXormqr = lwork - 4*n*n;
  scalar_t *Q = work, *workXormqr = work + 4*n*n;  

  if(lwork == -1) {
    lapackXormqr('R', 'T', N, 2*n, 2*n,
		 Q, ldq, tau, R1i, lda, work, -1, info);
    work[0] += 4*n*n;
    return *info;
  }

  R1i = A_BLK(0,L-2); 
  Qii = A_BLK(L-2,L-2); 

  //lapackXlacpy('A', 2*n, 2*n, Qii, lda, Q, ldq);
  lapackXcpyzero('L', 2*n, 2*n, Qii, lda, Q, ldq);
  // lapackXorgqr  (2*n, 2*n, 2*n, Q, ldq, tau, work, lwork, info);

  lapackXormqr('R', 'T', N, 2*n, 2*n,
   	       Q, ldq, tau, R1i, lda, workXormqr, lworkXormqr, info);

  for(i = L-3; i >= 0; i--) {
    R1i = A_BLK(0,i); 
    Qii = A_BLK(i,i); 
    tau -=n;

    lapackXcpyzero('L', 2*n, n, Qii, lda, Q, ldq);
    lapackXormqr('R', 'T', N, 2*n, n,
     		 Q, ldq, tau, R1i, lda, workXormqr, lworkXormqr, info);
  }

  return *info;
}

/**
 * Purpose:
 *   {\bf lapackXbsoi_ormqr} computes inverse of a factorized 
 *   $N$-by-$N$ $p$-cyclic matrix $A$ with pre-inverted factor $R$
 *   where $p = L$, $N = n*L$. 
 *   Uses expricit representation of $Q^{(k)}$ 
 *   for blocked Householder updates
 *********************************************************************
 * Arguments:
 * @arg n (input) INTEGER
 *   The block size of the matrix $A$ ($n > 0$).
 * @arg L (input) INTEGER
 *   The number of row/column blocks of the matrix $A$ ($L > 4$).
 * @arg A (input/output) DOUBLE_PRECISION array, dimension (LDA,N)
 *   On entry, matrix $A$ factorized by 
 *   structured orthogonal factorization method with inverted factor $R$. 
 *   Namely, 
 *   the Householder reflectors of factors $Q^k$ 
 *   and factor $R^{-1}$.
 *   The unit scaling coefficients of Householder reflectors are not stored.
 *   On exit, inverse of matrix $A$.
 * @arg lda (input) INTEGER
 *   The leading dimension of the array A. $lda >= max(1,n*L)$.
 * @arg tau (output) DOUBLE_PRECISION array, dimension ($n*L$)
 *   Scaling coefficients of the elementary Householder reflectors.
 * @arg work (workspace/output) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))
 *   On exit, if $info = 0$, WORK(1) returns the optimal LWORK.
 * @arg lwork (input) INTEGER
 *   The dimension of the array WORK. 
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 */
int lapackXbsoi_orgqr(int n, int L, 
		scalar_t A[], int lda, scalar_t _tau[], 
		scalar_t work[], int lwork, int* info)
{
  int N = n*L, i, ldq = 2*n;
  scalar_t /* *R1i, *Qii, */ *tau = _tau + (L-2)*n;
  int lworkXorgqr = lwork; scalar_t *workXorgqr = work;

  scalar_t *Q = workXorgqr;
  lworkXorgqr -= 2*n*ldq; workXorgqr += 2*n*ldq;

  scalar_t *hX1icp = workXorgqr;
  lworkXorgqr -= 2*n*lda; workXorgqr += 2*n*lda;

  if(lwork == -1) {
    lapackXorgqr  (2*n, 2*n, 2*n, Q, ldq, tau, work, -1, info);
    work[0] += 2*n*(ldq + lda);
    return *info;
  }

  /* R1i = A_BLK(0,L-2);  */
  /* Qii = A_BLK(L-2,L-2);  */

  //lapackXlacpy('A', 2*n, 2*n, Qii, lda, Q, ldq);
  lapackXcpyzero('L', 2*n, 2*n, A_BLK(L-2,L-2), lda, Q, ldq);
  lapackXorgqr  (2*n, 2*n, 2*n, Q, ldq, tau, workXorgqr, lworkXorgqr, info);
  //lapackXgemm   ('N', 'T', N, 2*n,2*n, 1, R1i, lda, Q, ldq, 1, R1i, lda);
  lapackXlacpy  ('A', N, 2*n, A_BLK(0,L-2), lda, hX1icp, lda);
  lapackXgemm   ('N', 'T', N, 2*n, 2*n, 1, hX1icp, lda, Q, ldq, 0, A_BLK(0,L-2), lda);

  for(i = L-3; i >= 0; i--) {
    /* R1i = A_BLK(0,i);  */
    /* Qii = A_BLK(i,i);  */
    tau -=n;

    // lapackXcpyzero('L', 2*n, n, Qii, lda, Q, ldq);
    // lapackXorgqr  (2*n, 2*n, n, Q, ldq, tau, workXorgqr, lworkXorgqr, info);
    // lapackXgemm   ('N', 'T', N, 2*n,2*n, 1, R1i, lda, Q, ldq, 1, R1i, lda);

    lapackXcpyzero('L', 2*n, n, A_BLK(i,i), lda, Q, ldq);
    lapackXorgqr  (2*n, 2*n, n, Q, ldq, tau, workXorgqr, lworkXorgqr, info);

    lapackXlacpy  ('A', (i+2)*n, n, A_BLK(0, i), lda, hX1icp, lda);
    lapackXlacpy  ('A', N, n, A_BLK(0, i + 1), lda, hX1icp + A_BLK_OFFSET(0, 1), lda);
    lapackXgemm   ('N', 'T', (i+1)*n, 2*n, 2*n, 1, hX1icp, lda, Q, ldq, 0, A_BLK(0, i), lda);
    lapackXgemm   ('N', 'T', N - (i+1)*n, 2*n, n, 1, 
		   hX1icp + A_BLK_OFFSET(i+1, 1), lda, Q + n*ldq, ldq, 0, A_BLK(i+1, i), lda);
  }

  return *info;
}
