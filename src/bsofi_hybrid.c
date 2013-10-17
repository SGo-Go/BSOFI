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
 * Attention: 
 *   This implementation of get_lk assumes static load balancing 
 *   between host and device 
 */
int get_lk(int k, int L, double kappa_Q, int ck1) /* const  */
{
  double lk;
  /* if(k*(KAPPA_Q + 2) <= KAPPA_Q*L - 2) */
  /*   lk = (L + k + 2.)/(1. + KAPPA_Q); */
  /* else */
  /*   lk = (L + 0.5*KAPPA_Q*(L - k) + 1.)/(1. + KAPPA_Q); */
  lk = (kappa_Q*L - 2.)/(kappa_Q + 2.) - ck1;
  if(lk < 0) lk = 0; 
  else if (lk > L) lk = L;
  return (int)(lk + 0.5);
}


int get_li(int i, int L, double kappa_R) /* const  */
{
  double li = (L - i + 1.)/(1. + kappa_R);
  if(li < 0) li = 0; 
  else if (li > L - i - 3) li = L - i - 3;
  return (int)(li + 0.5);
}


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
  scalar_t work[1]; /* int lwork = -1; */ int info;
  /* hybridXbsof(n, L, A, lda, tau, work, lworkHostDevice, &info); */
  hybridXbsoftri(handle, n, L, A, lda, tau, 
		 work, lworkHostDevice, dwork, ldwork, /* Fswitch, */ &info);
  lworkTotal = work[0]; 
  hybridXbsoi   (handle, n, L, A, lda, tau, 
		 work, lworkHostDevice, dwork, ldwork, /* Lswitch, */ &info);
  if(lworkTotal < work[0]) {lworkTotal = work[0];}
  return lworkTotal;
}

#ifdef USE_PROF
timeval_t tim1;
bsofi_profile_t profHybridBSOI, profHybridBSOFTRI;
#endif

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
		   //int Fswitch, 
		   int* info)
{
  RESET_BSOFI_PROFILE(profHybridBSOFTRI);

  int i, l_switch, l_switch_prev;
  scalar_t /* *Rii, *Riip1, *Rin, */ *tau = _tau;
  scalar_t max_work0; 
  int lddx = GET_LD_DEV(n), lddw = GET_LD_DEV(n);

  scalar_t *dX;
  scalar_t *dXarg;
  scalar_t *dXprod;
  scalar_t *dW;

  /* Get maximum value of l_switch for inversion */
  double kappa_R = GET_KAPPA_R(n);
  l_switch = get_li(0, L, kappa_R); /* L*Fswitch/1000; if(l_switch > L-3) l_switch = L-3; */

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
      //work[0] = 2*lddx*n*(L-3) + lddw*n;
    } else if (lwork == -2) {
      /* Workspace for triangular inversion */
      work[0] = 2*lddx*n*l_switch + lddw*n;
    }
    return *info;
  }

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
    BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrtri('U', 'N', n*L, A_BLK(0, 0), lda, info));
  } else {
    i = L - 3;
    BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrtri('U', 'N', 3*n, A_BLK(L-3, L-3), lda, info));

    if(L > 3) {
      BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrmm ('R', 'U', 'N', 'N', n*(L-3), n,  1, 
							       A_BLK(L-1, L-1), lda, A_BLK(0, L-1), lda));

      l_switch_prev = 0;
      for(i = L - 4; i >= 0; i--) {
	l_switch = get_li(i, L, kappa_R); /* l_switch = ((L - i)*Fswitch)/1000; if(l_switch > L-i-3) l_switch = L-i-3; */
	if(l_switch > l_switch_prev) {
	  BENCH_CUMMULATIVE(profHybridBSOFTRI.memset, cublasXlaset ('A', n, n*(l_switch - l_switch_prev), 
								    A_BLK(i+1, L-1 - l_switch), lda, 
								    dXarg + BLK_OFFSET(lddx, 0, L-1 - l_switch),  lddx));
	  l_switch_prev = l_switch;
	}

	BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrtri('U', 'N', n, A_BLK(i, i), lda, info));
	BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrmm ('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, A_BLK(i, i+1), lda));
	if (l_switch) 
	  BENCH_CUMMULATIVE(profHybridBSOFTRI.memset, cublasXlaset('A', n, n, A_BLK(i, i+1), lda, dW, lddw));

	BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrmm ('L', 'U', 'N', 'N', n, n, -1, A_BLK(i, i), lda, A_BLK(i, L-1), lda)); 

	if (l_switch) 
	  BENCH_CUMMULATIVE(profHybridBSOFTRI.gpu,    cublasXgemm ('N', 'N', n, n*l_switch, n, 1, dW, lddw, 
								   dXarg  + BLK_OFFSET(lddx, 0, L-1 - l_switch), lddx, 0, 
								   dXprod + BLK_OFFSET(lddx, 0, L-1 - l_switch), lddx));

	BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXgemm ('N', 'N', n, n*(L-i-3 - l_switch), n, 1, 
								 A_BLK(i, i+1), lda, A_BLK(i+1,i+2), lda, 0, A_BLK(i, i+2), lda));

	BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXgemm ('N', 'N', n, n, n, 1, A_BLK(i, i+1), lda, A_BLK(i+1,L-1), lda, 1, A_BLK(i, L-1), lda));
	BENCH_CUMMULATIVE(profHybridBSOFTRI.cpu,    lapackXtrmm ('R', 'U', 'N', 'N', n, n, 1, A_BLK(i+1, i+1), lda, A_BLK(i, i+1), lda)); 

	if (l_switch) 
	  BENCH_CUMMULATIVE(profHybridBSOFTRI.memget, cublasXlaget('A', n, n*l_switch,
								   dXprod + BLK_OFFSET(lddx, 0, L-1 - l_switch), lddx,
								   A_BLK(i, L-1 - l_switch), lda));

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
  RESET_BSOFI_PROFILE(profHybridBSOI);

  double kappa_Q = GET_KAPPA_Q(n);
  int ck1 = GET_Ck1(n);
  int k_switch = get_lk(0,L, kappa_Q, ck1);
  /* printf("\n\n%d\n\n", k_switch); */

  /* cublasStatus_t stat = 0; */
  int i, ldq = 2*n;
  scalar_t *tau = _tau + (L-2)*n;
  int lworkXorgqr = lwork; scalar_t *workXorgqr = work;

  scalar_t *Q  = workXorgqr;
  lworkXorgqr -= 2*n*ldq; 
  workXorgqr  += 2*n*ldq;

  scalar_t *hWtop  = workXorgqr;
  int lhwtop = (k_switch)*n;

  lworkXorgqr -= 2*n*lhwtop; 
  workXorgqr  += 2*n*lhwtop; 

  int lddq = GET_LD_DEV(ldq), lhwbottom = GET_LD_DEV((L-k_switch)*n);
  scalar_t *dQ       = dwork;
  scalar_t *hWbottom = dwork + 2*n*lddq;
  scalar_t *hRbottom = dwork + 2*n*lddq + 2*n*lhwbottom;
  scalar_t *houtAbottom;

  /* printf("\n%d %d\n", lhwtop,lhwbottom); */

  if(lwork == -1) {
    lapackXorgqr  (2*n, 2*n, 2*n, Q, ldq, tau, work, -1, info);
    work[0] += 2*n*ldq + 2*n*lhwtop;
    return *info;
  } else if (lwork == -2) {
    work[0] = 2*n*lddq + 4*n*lhwbottom;
    return *info;
  }

  BENCH_CUMMULATIVE(profHybridBSOI.memcpu, lapackXcpyzero('L', 2*n, 2*n, A_BLK(L-2,L-2), lda, Q, ldq));
  BENCH_CUMMULATIVE(profHybridBSOI.memset, cublasXlaset  ('A', n*(L-k_switch), 2*n, A_BLK(k_switch,L-2), lda, hWbottom, lhwbottom));

  BENCH_CUMMULATIVE(profHybridBSOI.cpu,    lapackXorgqr  (2*n, 2*n, 2*n, Q, ldq, tau, workXorgqr, lworkXorgqr, info));
  BENCH_CUMMULATIVE(profHybridBSOI.memset, cublasXlaset  ('A', 2*n, 2*n, Q, ldq, dQ,  lddq));

  BENCH_CUMMULATIVE(profHybridBSOI.gpu,    cublasXgemm   ('N', 'T', n*(L-k_switch), 2*n, 2*n, 1, hWbottom, lhwbottom, dQ, lddq, 0, hRbottom, lhwbottom));
  BENCH_CUMMULATIVE(profHybridBSOI.memget, cublasXlaget  ('A', n*(L-k_switch), n, hRbottom + BLK_OFFSET(lhwbottom, 0, 1), lhwbottom, A_BLK(k_switch,L-1), lda));

  if(lhwtop > 0) {
    BENCH_CUMMULATIVE(profHybridBSOI.memcpu, lapackXlacpy  ('A', n*k_switch, 2*n, A_BLK(0,L-2), lda, hWtop,  lhwtop));
    BENCH_CUMMULATIVE(profHybridBSOI.cpu,    lapackXgemm   ('N', 'T', n*(k_switch),   2*n, 2*n, 1, hWtop,  lhwtop,  Q, ldq, 0, A_BLK(0,L-2), lda));
}

  for(i = L-3; i >= 0; i--) {
    tau -=n;

    BENCH_CUMMULATIVE(profHybridBSOI.memcpu, lapackXcpyzero('L', 2*n, n, A_BLK(i,i), lda, Q, ldq));
    BENCH_CUMMULATIVE(profHybridBSOI.cpu,    lapackXorgqr  (2*n, 2*n, n, Q, ldq, tau, workXorgqr, lworkXorgqr, info));

    houtAbottom  = hWbottom; hWbottom = hRbottom; hRbottom  = houtAbottom;

    BENCH_CUMMULATIVE(profHybridBSOI.sync,   cudaDeviceSynchronize());
    if(i + 1 > k_switch) {
      BENCH_CUMMULATIVE(profHybridBSOI.memset, cublasXlaset  ('A', 2*n, 2*n, Q, ldq, dQ,  lddq));
      BENCH_CUMMULATIVE(profHybridBSOI.memset, cublasXlaset  ('A', (i+1 - k_switch)*n, n, A_BLK(k_switch, i),   lda, 
							      hWbottom + BLK_OFFSET(lhwbottom, 0, 1), lhwbottom));
      BENCH_CUMMULATIVE(profHybridBSOI.gpu,    cublasXgemm   ('N', 'T', (L - k_switch)*n, 2*n, n, 1,
							      hWbottom + BLK_OFFSET(lhwbottom, 0, 0), lhwbottom, 
							      dQ + BLK_OFFSET(lddq, 0, 1), lddq, 0, hRbottom, lhwbottom));
      BENCH_CUMMULATIVE(profHybridBSOI.gpu,    cublasXgemm   ('N', 'T', (i+1 - k_switch)*n, 2*n, n, 1, 
							      hWbottom + BLK_OFFSET(lhwbottom, 0, 1), lhwbottom,
							      dQ, lddq, 1, hRbottom, lhwbottom));

      if(lhwtop > 0) {
	BENCH_CUMMULATIVE(profHybridBSOI.memcpu, lapackXlacpy  ('A',      k_switch*n, 2*n, A_BLK(0, i), lda, hWtop, lhwtop));
	BENCH_CUMMULATIVE(profHybridBSOI.cpu,    lapackXgemm   ('N', 'T', k_switch*n, 2*n, 2*n, 1, hWtop, lhwtop, Q, ldq, 0, A_BLK(0, i), lda));
      }

    } else {
      BENCH_CUMMULATIVE(profHybridBSOI.memset, cublasXlaset  ('A', 2*n, n, Q + BLK_OFFSET(ldq, 0, 1), ldq, dQ + BLK_OFFSET(lddq, 0, 1),  lddq));
      BENCH_CUMMULATIVE(profHybridBSOI.gpu,    cublasXgemm   ('N', 'T', (L - k_switch)*n, 2*n, n, 1,
							      hWbottom + BLK_OFFSET(lhwbottom, 0, 0), lhwbottom, 
							      dQ + BLK_OFFSET(lddq, 0, 1), lddq, 0, hRbottom, lhwbottom));

      if(lhwtop > 0) {
	BENCH_CUMMULATIVE(profHybridBSOI.memcpu, lapackXlacpy  ('A', (i+1)*n,      2*n, A_BLK(0, i), lda, hWtop, lhwtop));
	BENCH_CUMMULATIVE(profHybridBSOI.cpu,    lapackXgemm   ('N', 'T', (i+1)*n, 2*n, 2*n, 1, hWtop, lhwtop, Q, ldq, 0, A_BLK(0, i), lda));
      }
      if(i + 1 < k_switch) {
	if(lhwtop > 0) {
	  BENCH_CUMMULATIVE(profHybridBSOI.memcpu, lapackXlacpy  ('A', (k_switch - (i+1))*n, n, A_BLK(i+1, i+1), lda, 
								  hWtop + BLK_OFFSET(lhwtop, i+1, 1), lhwtop));
	  BENCH_CUMMULATIVE(profHybridBSOI.cpu,    lapackXgemm   ('N', 'T', (k_switch - (i+1))*n, 2*n, n, 1,
								  hWtop + BLK_OFFSET(lhwtop,i+1, 1), lhwtop, Q + n*ldq, ldq, 0, A_BLK(i+1, i), lda));
	}
      }

    }

    BENCH_CUMMULATIVE(profHybridBSOI.memget, cublasXlaget  ('A', n*(L-k_switch), n, hRbottom + BLK_OFFSET(lhwbottom, 0, 1), lhwbottom, 
							    A_BLK(k_switch, i + 1), lda));
  }
  BENCH_CUMMULATIVE(profHybridBSOI.memget, cublasXlaget  ('A', n*(L-k_switch), n, hRbottom + BLK_OFFSET(lhwbottom, 0, 0), lhwbottom, 
							  A_BLK(k_switch, 0), lda));

  BENCH_CUMMULATIVE(profHybridBSOI.sync,   cudaDeviceSynchronize());
  return *info;
}
