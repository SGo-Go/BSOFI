/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   CPU+GPU implementation of routines based on LAPACK 
 *   and CuBLAS subroutines calls)
 */

#include "bsofi.h"
#include <config.h>
#include <third_party/lapack.h>
#include "bsofi_macro.h"
#include <third_party/cublas.h>
#include "bsofi_params.h"

/* ===================================================================
   BSOF on CPU
   =================================================================== */

/**
 * Purpose:
 *   {\bf hybridXbsofiLWork} estimates workspace required to invert
 *   a $N$-by-$N$ $p$-cyclic matrix $A$ where $p = L$, $N = n*L$
 *   by structured orthogonal invertion method (BSOF/I subroutines). 
 * ********************************************************************
 * Arguments:
 * @arg handle (input) cublasHandle_t
 *   The handler for threading in CuBLAS
 * @arg n (input) INTEGER
 *   The block size of the matrix $A$ ($n > 0$).
 * @arg L (input) INTEGER
 *   The number of row/column blocks of the matrix $A$ ($L > 4$).
 * @arg lda (input) INTEGER
 *   The leading dimension of the array A. $lda >= max(1,n*L)$.
 * @arg lworkHostDevice (input) INTEGER
 *   CPU/GPU indicator. 
 *   = -1 -- CPU
 *   = -2 -- GPU
 */
int hybridXbsofiLWork(int n, int L, int lda, int lworkHostDevice)
{
  cublasHandle_t handle = 0;
  int lworkTotal = 0;
  scalar_t *A = 0; scalar_t *tau = 0; 
  scalar_t dwork[1]; int ldwork = -1; 
  scalar_t work[1]; int info;
  hybridXbsoftri(handle, n, L, A, lda, tau, 
		 work, lworkHostDevice, dwork, ldwork, &info);
  lworkTotal = work[0]; 
  hybridXbsoi   (handle, n, L, A, lda, tau, 
		 work, lworkHostDevice, dwork, ldwork, &info);
  if(lworkTotal < work[0]) {lworkTotal = work[0];}
  return lworkTotal;
}

/**
 * Purpose:
 *   {\bf hybridXbsoftri} computes a structured orthogonal factorization 
 *   of a $N$-by-$N$ $p$-cyclic matrix $A$ where $p = L$, $N = n*L$
 *   and inverts factor $R$. 
 *   The factorization has the form
 *   $A = Q^{L-1} * ... Q^1 * R$
 *   $R$ is an upper triangular, $Q^k$ are orthogonal.
 *   This is the Level 3 BLAS version of the algorithm.
 * ********************************************************************
 * Arguments:
 * @arg handle (input) cublasHandle_t
 *   The handler for threading in CuBLAS
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
 * @arg dwork (workspace/output) DOUBLE_PRECISION array on the GPU, dimension (MAX(1,LDWORK))
 *   Work array.
 * @arg ldwork (input) INTEGER
 *   The dimension of the array DWORK. 
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
int hybridXbsoftri(cublasHandle_t handle, int n, int L, 
		   scalar_t A[], int lda, scalar_t _tau[], 
		   scalar_t work[], int lwork, scalar_t dwork[], int ldwork,
		   int* info)
{
  int i, l_switch, l_switch_prev;
  scalar_t *tau = _tau;
  scalar_t max_work0; 
  int lddx = GET_LD_DEV(n), lddw = GET_LD_DEV(n);

  scalar_t *dX;
  scalar_t *dXarg;
  scalar_t *dXprod;
  scalar_t *dW;

  /* Get maximum value of l_switch for inversion */
  double kappa_R = KAPPA_R(n);
  int ci         = C_I(n);

  l_switch = get_li(0, L, kappa_R, ci);

  dXarg  = dwork                   - (L - l_switch - 1)*n*lddx;
  dXprod = dwork + lddx*n*l_switch - (L - l_switch - 1)*n*lddx;
  dW     = dwork + 2*lddx*n*l_switch;

  if(lwork < 0) {
    if(lwork == -1) {
      /* Workspace for factorization */
      lapackXgeqrf(2*n, 2*n, A, lda, tau, work, lwork, info);
      max_work0 = work[0];
      lapackXormqr('L', 'T', 2*n, n, n,
		   A, lda, tau, A, lda, work, lwork, info);
      if(max_work0 > work[0]) work[0] = max_work0;
    } else if (lwork == -2) {
      /* Workspace for triangular inversion */
      work[0] = 2*lddx*n*l_switch + lddw*n;
    }
    return *info;
  }

  /************************************************************
   * Factorization code
   ************************************************************/
  for(i = 0; i < L-2; i++) {
    lapackXgeqrf(2*n, n, A_BLK(i, i), lda, tau, work, lwork, info);
    lapackXormqr('L', 'T', 2*n, n, n,
    		 A_BLK(i, i), lda, tau, A_BLK(i, i+1), lda, work, lwork, info);
    lapackXormqr('L', 'T', 2*n, n, n,
     		 A_BLK(i, i), lda, tau, A_BLK(i, L-1), lda, work, lwork, info);
    tau +=n;
  }

  lapackXgeqrf(2*n, 2*n, A_BLK(L-2,L-2), lda, tau, work, lwork, info);

  /************************************************************
   * Inversion code
   ************************************************************/
  if(L <= 2) { 
    lapackXtrtri('U', 'N', n*L, A_BLK(0, 0), lda, info);
  } else {
    i = L - 3;
    lapackXtrtri('U', 'N', 3*n, A_BLK(L-3, L-3), lda, info);

    if(L > 3) {
      lapackXtrmm ('R', 'U', 'N', 'N', n*(L-3), n,  1, 
		   A_BLK(L-1, L-1), lda, A_BLK(0, L-1), lda);

      l_switch_prev = 0;
      for(i = L - 4; i >= 0; i--) {
	l_switch = get_li(i, L, kappa_R, ci);
	if(l_switch > l_switch_prev) {
	  cublasXlaset ('A', n, n*(l_switch - l_switch_prev), 
			A_BLK(i+1, L-1 - l_switch), lda, 
			dXarg + BLK_OFFSET(lddx, 0, L-1 - l_switch),  lddx);
	  l_switch_prev = l_switch;
	}

	lapackXtrtri('U', 'N', n, A_BLK(i, i), lda, info);
	lapackXtrmm ('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, A_BLK(i, i+1), lda);
	if (l_switch) 
	  cublasXlaset('A', n, n, A_BLK(i, i+1), lda, dW, lddw);

	lapackXtrmm ('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, A_BLK(i, L-1), lda); 

	if (l_switch) 
	  cublasXgemm ('N', 'N', n, n*l_switch, n, 1, dW, lddw, 
		       dXarg  + BLK_OFFSET(lddx, 0, L-1 - l_switch), lddx, 0, 
		       dXprod + BLK_OFFSET(lddx, 0, L-1 - l_switch), lddx);

	lapackXgemm ('N', 'N', n, n*(L-i-3 - l_switch), n, 1, 
		     A_BLK(i, i+1), lda, A_BLK(i+1,i+2), lda, 0, A_BLK(i, i+2), lda);

	lapackXgemm ('N', 'N', n, n, n, 1, A_BLK(i, i+1), lda, A_BLK(i+1,L-1), lda, 1, A_BLK(i, L-1), lda);
	lapackXtrmm ('R', 'U', 'N', 'N', n, n, 1, A_BLK(i+1, i+1), lda, A_BLK(i, i+1), lda); 

	if (l_switch) 
	  cublasXlaget('A', n, n*l_switch,
		       dXprod + BLK_OFFSET(lddx, 0, L-1 - l_switch), lddx,
		       A_BLK(i, L-1 - l_switch), lda);

	dX = dXarg; dXarg = dXprod; dXprod = dX;
      }
    }
  }
  return *info;
}

/**
 * Purpose:
 *   {\bf hybridXbsoi} computes inverse of a factorized 
 *   $N$-by-$N$ $p$-cyclic matrix $A$ with pre-inverted factor $R$
 *   where $p = L$, $N = n*L$. 
 *   Uses expricit representation of $Q^{(k)}$ 
 *   for blocked Householder updates.
 *********************************************************************
 * Arguments:
 * @arg handle (input) cublasHandle_t
 *   The handler for threading in CuBLAS
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
 *   On exit, if $info = 0$, $lwork=-1$, WORK(1) returns the optimal LWORK;
 *   if $info = 0$, $lwork=-1$, WORK(1) returns the optimal LWORK.
 * @arg lwork (input) INTEGER
 *   The dimension of the array WORK. 
 * @arg dwork (workspace/output) DOUBLE_PRECISION array on the GPU, dimension (MAX(1,LDWORK))
 *   Work array.
 * @arg ldwork (input) INTEGER
 *   The dimension of the array DWORK. 
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 */
int hybridXbsoi(cublasHandle_t handle, int n, int L, 
		scalar_t A[], int lda, scalar_t _tau[], 
		scalar_t work[], int lwork, scalar_t dwork[], int ldwork, 
		// int Lswitch, 
		int* info)
{
  int n_corr = (int)(n*pow(L>>1,1./3));
  double kappa_Q = KAPPA_Q(n_corr);
  int ck1 = C_K1(n);
  int ck2 = C_K1(n);

  int k;
  int lk_old, lk = 0;

  /* cublasStatus_t stat = 0; */
  int ldq = 2*n;
  scalar_t *tau = _tau + (L-2)*n;
  int lworkXorgqr = lwork; scalar_t *workXorgqr = work;

  scalar_t *Q  = workXorgqr;
  lworkXorgqr -= 2*n*ldq; 
  workXorgqr  += 2*n*ldq;

  scalar_t *hWtop  = workXorgqr;
  int lhtop = L*n;

  lworkXorgqr -= 2*n*lhtop; 
  workXorgqr  += 2*n*lhtop; 

  int lddq = GET_LD_DEV(ldq), ldbottom = GET_LD_DEV(L*n);
  scalar_t *dQ       = dwork;
  scalar_t *hWbottom = dwork + 2*n*lddq;
  scalar_t *hXbottom = dwork + 2*n*lddq + 2*n*ldbottom;
  scalar_t *houtAbottom;

  if(lwork == -1) {
    lapackXorgqr  (2*n, 2*n, 2*n, Q, ldq, tau, work, -1, info);
    work[0] += 2*n*ldq + 2*n*lhtop;
    return *info;
  } else if (lwork == -2) {
    work[0] = 2*n*lddq + 4*n*ldbottom;
    return *info;
  }

  /* Assume that in the beginning 
     GPU does not have any data about matrix X */
  lk = 0;
  /* Computation cycle over blocks of reflectors */
  for(k = L-2; k >= 0; k--) {
    tau -=n;
    /* Switch multiplication result buffer and multiplication input buffer on GPU */
    houtAbottom  = hWbottom; hWbottom = hXbottom; hXbottom  = houtAbottom;

    /* First iteration is special in the sence 
       that we have 2*n Householder reflectors */
    lapackXcpyzero('L', (lk ? 1 : 2)*n, n, A_BLK(k,k), lda, Q, ldq);
    /* Reconstruct Q_k from Householder reflectors */
    lapackXorgqr  (2*n, (lk ? 1 : 2)*n, n, 
		   Q, ldq, tau, workXorgqr, lworkXorgqr, info);

    /* Copy column blocks from GPU of CPU if needed */
    lk_old = lk; lk = get_lk(k, L, kappa_Q, ck1, ck2);
    if(lk > lk_old)
      cublasXlaset  ('A', (lk - lk_old)*n, n, 
		     A + BLK_OFFSET(lda, L - lk, k+1), lda, 
		     hWbottom + BLK_OFFSET(ldbottom, L - lk, 0), ldbottom);

    /* DBGERROR("%d %d %d %d", lhtop, ldbottom, lk_old, lk); */

    cudaDeviceSynchronize();
    
    /* Multiplication on GPU + CPU */
    if(k + lk >= L) { /* if k + l_k > p */
      if(ldbottom > 0) {
	cublasXlaset  ('A', 2*n, 2*n, 
		       Q, ldq, 
		       dQ,  lddq);
	cublasXlaset  ('A', (k + 1 - L + lk)*n, n, 
		       A + BLK_OFFSET(lda, L - lk, k),   lda, 
		       hWbottom + BLK_OFFSET(ldbottom, L - lk, 1), ldbottom);
	/* DBGERROR("1"); */
	cublasXgemm   ('N', 'T', lk*n, 2*n, n, 1,
		       hWbottom + BLK_OFFSET(ldbottom, L - lk, 0), ldbottom, 
		       dQ + BLK_OFFSET(lddq, 0, 1), lddq, 0, 
		       hXbottom + BLK_OFFSET(ldbottom, L - lk, 0), ldbottom);
	/* DBGERROR("2 %d", (k + 1 - L + lk)); */
	cublasXgemm   ('N', 'T', (k + 1 - L + lk)*n, 2*n, n, 1, 
		       hWbottom + BLK_OFFSET(ldbottom, L - lk, 1), ldbottom,
		       dQ, lddq, 1, 
		       hXbottom + BLK_OFFSET(ldbottom, L - lk, 0), ldbottom);
      }
      if(lhtop > 0) {
	lapackXlacpy  ('A',      (L - lk)*n, 2*n, 
		       A + BLK_OFFSET(lda, 0, k), lda, 
		       hWtop, lhtop);
	/* DBGERROR("3"); */
	lapackXgemm   ('N', 'T', (L - lk)*n, 2*n, 2*n, 1, 
		       hWtop, lhtop, 
		       Q, ldq, 0, 
		       A + BLK_OFFSET(lda, 0, k), lda);
      }

    } else { /* if k + l_k <= p */
      if(ldbottom > 0) {
	cublasXlaset  ('A', 2*n, n, 
		       Q  + BLK_OFFSET(ldq,  0, 1), ldq, 
		       dQ + BLK_OFFSET(lddq, 0, 1), lddq);
	/* DBGERROR("4"); */
	cublasXgemm   ('N', 'T', lk*n, 2*n, n, 1,
		       hWbottom + BLK_OFFSET(ldbottom, L - lk, 0), ldbottom, 
		       dQ + BLK_OFFSET(lddq, 0, 1), lddq, 0, 
		       hXbottom + BLK_OFFSET(ldbottom, L - lk, 0), ldbottom);
      }
      if(lhtop > 0) {
	lapackXlacpy  ('A', (k+1)*n,      2*n, 
		       A + BLK_OFFSET(lda, 0, k), lda, 
		       hWtop, lhtop);
	/* DBGERROR("5"); */
	lapackXgemm   ('N', 'T', (k+1)*n, 2*n, 2*n, 1, 
		       hWtop, lhtop, 
		       Q, ldq, 0, 
		       A + BLK_OFFSET(lda, 0, k), lda);
	if(k + 1 < L - lk) {
	  lapackXlacpy  ('A', (L - lk - k - 1)*n, n, 
			 A + BLK_OFFSET(lda, k+1, k+1), lda, 
			 hWtop + BLK_OFFSET(lhtop, k+1, 1), lhtop);
	  /* DBGERROR("6"); */
	  lapackXgemm   ('N', 'T', (L - lk - k - 1)*n, 2*n, n, 1,
			 hWtop + BLK_OFFSET(lhtop, k+1, 1), lhtop, 
			 Q + BLK_OFFSET(ldq, 0, 1), ldq, 0, 
			 A + BLK_OFFSET(lda, k+1, k), lda);
	}
      }
    }

    if(ldbottom > 0) {
      cublasXlaget  ('A', lk*n, n, 
		     hXbottom + BLK_OFFSET(ldbottom, L - lk, 1), ldbottom, 
		     A + BLK_OFFSET(lda, L - lk, k + 1), lda);
    }
  }
  if(ldbottom > 0) {
    cublasXlaget  ('A', lk*n, n, 
		   hXbottom + BLK_OFFSET(ldbottom, 0, 0), ldbottom, 
		   A + BLK_OFFSET(lda, L - lk, 0), lda);
  }
  cudaDeviceSynchronize();
  return *info;
}
