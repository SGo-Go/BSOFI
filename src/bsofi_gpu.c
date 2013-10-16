/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   GPU implementation of routines based on Magma
 *   and CuBLAS subroutines calls
 */

#include <config.h>

#include <magma.h>
#include "bsofi.h"

#include "bsofi_macro.h"
#define gpuDgemm cublasDgemm
#define gpuDtrmm cublasDtrmm

#define dA_IDX_OFFSET(__i, __j) ((__j)*ldda + __i)
#define dA_BLK_OFFSET(__i, __j) (n*((__j)*ldda + __i))

#define dA_IDX(__i, __j) dA[dA_IDX_OFFSET(__i, __j)]
#define dA_BLK(__i, __j) dA + dA_BLK_OFFSET(__i, __j)

/**
 * Purpose:
 *   {\bf gpuXbsof} computes a structured orthogonal factorization 
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
 * @arg dT (workspace/output) DOUBLE_PRECISION array on the GPU
 *   It contains $L$ blocks.
 *   Each block starts with subblock that store the triangular $T$ matrices, 
 *   followed by the subblock of the diagonal inverses for the $R$ matrix. 
 *   The rest of the array is used as workspace.
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
int gpuXbsof(int n, int L, 
	     scalar_t dA[], int ldda, scalar_t _tau[], 
	     scalar_t work[], int lwork, scalar_t* _dT, int* info)
{
  int N = n*L, i;
  scalar_t *Rii, *Riip1, *Rin;
  scalar_t *tau = _tau, *dT = _dT;

  magma_int_t nb            = magma_get_dgeqrf_nb(2*n);
  magma_int_t dt_size       = get_QR_T_device_mem(2*n, n, nb);
  magma_int_t dt_size_total = (L-2)*dt_size + get_QR_T_device_mem(2*n, 2*n, nb);

  if(lwork == -1) {
    //lapackXgeqrf(2*n, 2*n, dA, ldda, tau, work, -1, info);
    //lwork = 0;work[0];
    gpuXormqr_gpu('L', 'T', 2*n, n, n,
		     dA, ldda, tau, dA, ldda, work, -1, dT, nb, info);
    if(lwork > work[0]) work[0] = lwork;
    if(lwork > work[0]) work[0] = lwork;
    return *info;
  } else if (lwork == -2)  {
    work[0] = dt_size_total;
    return 0;
  }

  for(i = 0; i < L-2; i++) {
    Rii   = dA_BLK(i,i); 
    Riip1 = dA_BLK(i, i+1); 
    Rin   = dA_BLK(i,L-1);

    gpuXgeqrf_gpu(2*n, n, Rii, ldda, tau, dT, info);
    /* magma_dormqr2_gpu('L', 'T', 2*n, n, n, */
    /* 		      Rii, ldda, tau, Riip1, ldda, dT, dt_size_total, info); */
    /* printf("%d\n",1); */
    /* magma_dormqr2_gpu('L', 'T', 2*n, n, n, */
    /* 		      Rii, ldda, tau, Rin, ldda, dT, dt_size_total, info); */

    gpuXormqr_gpu('L', 'T', 2*n, n, n,
    		Rii, ldda, tau, Riip1, ldda, work, lwork, dT, nb, info);
    gpuXormqr_gpu('L', 'T', 2*n, n, n,
    		  Rii, ldda, tau, Rin, ldda, work, lwork, dT, nb, info);
    tau +=n;
    dT += dt_size; 
    /* dt_size_total -= dt_size; */
  }

  Rii   = dA_BLK(L-2,L-2); 
  gpuXgeqrf_gpu(2*n, 2*n, Rii, ldda, tau, dT, info);
  return *info;
}

/************************************************************
 * xBSTRI routine implementations
 ************************************************************/

/**
 * Purpose:
 *   {\bf gpuXbstri_rseq} computes inverse of factor $R$ for a factorized 
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
 * @arg dwork (workspace/output) DOUBLE_PRECISION array on the GPU, dimension (MAX(1,LWORK))
 *   Work array.
 * @arg ldwork (input) INTEGER
 *   The dimension of the array DWORK. 
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 */
int gpuXbstri_rseq(int n, int L, 
		  scalar_t dA[], int ldda, 
		  scalar_t dwork[], int ldwork, 
		  int* info)
{
  int N = n*L, i = L - 3;

  if(L <= 2)
    gpuXtrtri_gpu('U', 'N', N, dA_BLK(0, 0), ldda, info);
  else {
    gpuXtrtri_gpu('U', 'N', 3*n, dA_BLK(L-3, L-3), ldda, info);

    if(L > 3) {
      gpuXtrmm_gpu('R', 'U', 'N', 'N', n*(L-3), n,  1, dA_BLK(L-1, L-1), ldda, dA_BLK(0, L-1), ldda);
      for(i = L - 4; i >= 0; i--) {
	gpuXtrtri_gpu('U', 'N', n, dA_BLK(i, i), ldda, info);
	gpuXtrmm_gpu ('L', 'U', 'N', 'N', n, n, -1, dA_BLK(i, i), ldda, dA_BLK(i, L-1), ldda); 

	gpuXtrmm_gpu ('L', 'U', 'N', 'N', n, n, -1, dA_BLK(i, i), ldda, dA_BLK(i, i+1), ldda);
	gpuXgemm_gpu ('N', 'N', n, n*(L-i-2), n, 1, dA_BLK(i, i+1), ldda, 
		      dA_BLK(i+1,i+2), ldda, 1, dA_BLK(i, i+2), ldda);
	gpuXtrmm_gpu ('R', 'U', 'N', 'N', n, n, 1, dA_BLK(i+1, i+1), ldda, dA_BLK(i, i+1), ldda); 
      }
    }
  }
  return *info;
}

int gpuXbstri_rbatched(int n, int L, 
		       scalar_t dA[], int ldda, 
		       scalar_t dwork[], int ldwork, 
		       int* info)
{
  int N = n*L, i;

  if(L <= 2)
    gpuXtrtri_gpu('U', 'N', N, dA, ldda, info);
  else {
    gpuXtrtri_gpu('U', 'N', 3*n, dA_BLK(L-3, L-3), ldda, info);

    if(L > 3) {
      gpuXtrmm_gpu('R', 'U', 'N', 'N', n*(L-3), n,  1, dA_BLK(L-1, L-1), ldda, dA_BLK(0, L-1), ldda);

      for(i = 0; i < L - 3; i++) {
	gpuXtrtri_gpu('U', 'N', n, dA_BLK(i, i), ldda, info);

	gpuXtrmm_gpu ('L', 'U', 'N', 'N', n, n, -1, dA_BLK(i, i), ldda, dA_BLK(i, i+1), ldda);
	gpuXtrmm_gpu ('L', 'U', 'N', 'N', n, n, -1, dA_BLK(i, i), ldda, dA_BLK(i, L-1), ldda); 
      }

      for(i = L - 4; i >= 0; i--) {
	gpuXgemm_gpu ('N', 'N', n, n*(L-i),n, 1, dA_BLK(i, i+1), ldda, 
		      dA_BLK(i+1,i+2), ldda, 1, dA_BLK(i, i+2), ldda);
	gpuXtrmm_gpu ('R', 'U', 'N', 'N', n, n, 1, dA_BLK(i+1, i+1), ldda, dA_BLK(i, i+1), ldda); 
      }
    }
  }
  return *info;
}

/************************************************************
 * xBSOI routine implementations
 ************************************************************/

int gpuXbsoi_ormqr(int n, int L, 
		   scalar_t dA[], int ldda, scalar_t _tau[], 
		   scalar_t work[], int lwork, 
		   scalar_t *_dT, scalar_t dwork[], int ldwork, 
		   int* info)
{
  magma_int_t nb      = magma_get_dgeqrf_nb(2*n);
  magma_int_t dt_size = get_QR_T_device_mem(2*n, n, nb);
  magma_int_t dt_size_total = (L-2)*dt_size + get_QR_T_device_mem(2*n, 2*n, nb);

  int N = n*L, i, lddq = GET_LD_DEV(2*n);
  scalar_t *R1i = dA, *Qii, 
    *tau = _tau + (L-2)*n, *dT = _dT + (L-2)*dt_size;
  scalar_t *dQ = dwork;
  

  if(lwork == -1) {
    *info = 0;
    //lapackXgeqrf(2*n, 2*n, dA, ldda, tau, work, -1, info);
    //lwork = 0;work[0];
    gpuXormqr_gpu('R', 'T', N, 2*n, 2*n,
		  dQ, lddq, tau, R1i, ldda, work, -1, dT, nb, info);
    //if(lwork > work[0]) work[0] = lwork;
    // work[0] *=10;
    //DBG("!!!" << lwork << " " << work[0])
    return *info;
  } else if (lwork == -2 || ldwork == -1) {
    work[0] = 2*n*lddq;
    return 0;
  }

  R1i = dA_BLK(0,L-2); 
  Qii = dA_BLK(L-2,L-2); 

  gpuXcpyzero_gpu('L', 2*n, 2*n, Qii, ldda, dQ, lddq);
  gpuXormqr_gpu('R', 'T', N, 2*n, 2*n,
   		dQ, lddq, tau, R1i, ldda, work, lwork, dT, nb, info);
  //DBG("!!!" << lwork)
  for(i = L-3; i >= 0; i--) {
    R1i  = dA_BLK(0,i); 
    Qii  = dA_BLK(i,i); 
    tau -= n; dT -= dt_size;

    gpuXcpyzero_gpu('L', 2*n, n, Qii, ldda, dQ, lddq);
    gpuXormqr_gpu('R', 'T', N, 2*n, n,
     		  dQ, lddq, tau, R1i, ldda, work, lwork, dT, nb, info);
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
 *   On exit, if $info = 0$, $lwork=-1$, WORK(1) returns the optimal LWORK;
 *   if $info = 0$, $lwork=-1$, WORK(1) returns the optimal LWORK.
 * @arg lwork (input) INTEGER
 *   The dimension of the array WORK. 
 * @arg dT (workspace/output) DOUBLE_PRECISION array on the GPU
 *   It contains $L$ blocks.
 *   Each block starts with subblock that store the triangular $T$ matrices, 
 *   followed by the subblock of the diagonal inverses for the $R$ matrix. 
 *   The rest of the array is used as workspace.
 * @arg dwork (workspace/output) DOUBLE_PRECISION array on the GPU, dimension (MAX(1,LDWORK))
 *   Work array.
 * @arg ldwork (input) INTEGER
 *   The dimension of the array DWORK. 
 * @arg info (output) INTEGER
 *   $= 0$: successful exit
 *   $< 0$: if $info = -i$, the $i$th argument had an illegal value.
 */
int gpuXbsoi_orgqr(int n, int L, 
		   scalar_t dA[], int ldda, scalar_t _tau[], 
		   scalar_t work[], int lwork, 
		   scalar_t *_dT, scalar_t dwork[], int ldwork, 
		   int* info)
{
  magma_int_t nb      = magma_get_dgeqrf_nb(2*n);
  magma_int_t dt_size = get_QR_T_device_mem(2*n, n, nb);
  magma_int_t dt_size_total = (L-2)*dt_size + get_QR_T_device_mem(2*n, 2*n, nb);

  int N = n*L, i, lddq = GET_LD_DEV(2*n);
  scalar_t *X1i = dA, *Qii, 
    *tau = _tau + (L-2)*n, *dT = _dT + (L-2)*dt_size;
  scalar_t *dQ = dwork;
  scalar_t *dX1icp = dwork + 2*n*lddq;

  if(lwork == -1) {
    *info = 0;
    //gpuXorgqr_gpu(2*n, 2*n, 2*n, dQ, lddq, tau, dT, nb, info);
    work[0] = 0;
    return *info;
  } else if (lwork == -2 || ldwork == -1) {
    work[0] = 2*n*lddq + 2*n*ldda;
    return 0;
  }

  X1i = dA_BLK(0,L-2); 
  Qii = dA_BLK(L-2,L-2); 

  gpuXcpyzero_gpu('L', 2*n, 2*n, Qii, ldda, dQ, lddq);
  gpuXorgqr_gpu  (2*n, 2*n, 2*n, dQ, lddq, tau, dT, nb, info);
  gpuXlacpy_gpu  ('A', N, 2*n, X1i, ldda, dX1icp, ldda);
  gpuXgemm_gpu   ('N', 'T', N, 2*n, 2*n, 1, dX1icp, ldda, dQ, lddq, 0, X1i, ldda);

  int N_panel = N - n;
  for(i = L-3; i >= 0; i--) {
    X1i  = dA_BLK(0,i); 
    Qii  = dA_BLK(i,i); 
    tau -= n; dT -= dt_size; N_panel -=n;


    gpuXcpyzero_gpu('L', 2*n, n, Qii, ldda, dQ, lddq);
    gpuXorgqr_gpu  (2*n, 2*n, n, dQ, lddq, tau, dT, nb, info);
    gpuXlacpy_gpu  ('A', N_panel + n, n, X1i, ldda, dX1icp, ldda);
    gpuXlacpy_gpu  ('A', N, n, X1i + dA_BLK_OFFSET(0, 1), ldda, dX1icp + dA_BLK_OFFSET(0, 1), ldda);
    gpuXgemm_gpu   ('N', 'T', N, 2*n, 2*n, 1, dX1icp, ldda, dQ, lddq, 0, X1i, ldda); //TODO << find bugs (don't work)
  }

  return *info;
}

int gpuXbsoi_batched(int n, int L, 
		     scalar_t dA[], int ldda, scalar_t _tau[], 
		     scalar_t work[], int lwork, 
		     scalar_t *_dT, scalar_t dwork[], int ldwork, 
		     int* info)
{
  magma_int_t nb      = magma_get_dgeqrf_nb(2*n);
  magma_int_t dt_size = get_QR_T_device_mem(2*n, n, nb);
  magma_int_t dt_size_total = (L-2)*dt_size + get_QR_T_device_mem(2*n, 2*n, nb);

  int N = n*L, i, lddq = GET_LD_DEV(2*n);
  scalar_t *X1i = dA, *Qii, 
    *tau = _tau + (L-2)*n, *dT = _dT + (L-2)*dt_size;
  scalar_t *dQ = dwork;
  scalar_t *dW = dwork + 2*n*lddq;

  if(lwork == -1) {
    *info = 0;
    //gpuXorgqr_gpu(2*n, 2*n, 2*n, dQ, lddq, tau, dT, nb, info);
    work[0] = 0;
    return *info;
  } if (lwork == -2)  {
    work[0] = 2*n*lddq + 2*n*ldda;
    return 0;
  }

  gpuXcpyzero_gpu('L', 2*n, 2*n, dA_BLK(L-2,L-2), ldda, dQ, lddq);
  gpuXorgqr_gpu  (2*n, 2*n, 2*n, dQ, lddq, tau, dT, nb, info);
  gpuXlacpy_gpu  ('A', N, 2*n, dA_BLK(0,L-2), ldda, dW, ldda);
  gpuXgemm_gpu   ('N', 'T', N, 2*n, 2*n, 1, dW, ldda, dQ, lddq, 0, dA_BLK(0,L-2), ldda);

  for(i = L-3; i >= 0; i--) {
    // X1i  = dA_BLK(0,i); 
    // Qii  = dA_BLK(i,i); 
    tau -= n; dT -= dt_size;

    gpuXcpyzero_gpu('L', 2*n, n, dA_BLK(i,i), ldda, dQ, lddq);
    gpuXorgqr_gpu  (2*n, 2*n, n, dQ, lddq, tau, dT, nb, info);

    gpuXlacpy_gpu  ('A', (i+2)*n, n, dA_BLK(0, i), ldda, dW, ldda);
    gpuXlacpy_gpu  ('A', N, n, dA_BLK(0, i + 1), ldda, dW + dA_BLK_OFFSET(0, 1), ldda);
    gpuXgemm_gpu   ('N', 'T', (i+1)*n, 2*n, 2*n, 1, dW, ldda, dQ, lddq, 0, dA_BLK(0, i), ldda);
    gpuXgemm_gpu   ('N', 'T', N - (i+1)*n, 2*n, n, 1, 
      		    dW + dA_BLK_OFFSET(i+1, 1), ldda, dQ+ n*lddq, lddq, 0, dA_BLK(i+1, i), ldda);
  }
  return *info;
}

/************************************************************
 * Auxiliary routines
 ************************************************************/
inline
magma_int_t get_QR_T_device_mem(int m, int n, int nb)
{
  magma_int_t mn = m > n ? m : n;
  magma_int_t dt_size = (2*mn + GET_LD_DEV(n))*nb;
  return dt_size;
}

inline
void gpuXcpyzero_gpu(char uplo, int m, int n, 
		     scalar_t *dA, int ldda, scalar_t *dB, int lddb)
{
  switch(uplo) {
  case 'L':
    gpuXlacpy_gpu('A', m, n, dA, ldda, dB, lddb); break;
  default:
    gpuXlacpy_gpu(uplo, m, n, dA, ldda, dB, lddb);
  }
  gpuXzero_gpu(uplo, m, n, dA, ldda);
}

inline
void gpuXzero_gpu(char uplo, int m, int n, 
		  scalar_t *dA, int ldda)
{
  cudaError_t info;
  int j, i; 
  switch(uplo) {
  case 'L':
    for (j = 0; j < n-1; j++)
      info = cudaMemset(dA + dA_IDX_OFFSET(j+1,j), 0x00, (m - j - 1)*sizeof(double));
    break;
  case 'U':
    for (j = 1; j < n; j++)
      info = cudaMemset(dA + dA_IDX_OFFSET(0,j), 0x00, (j - 1)*sizeof(double));
    break;
  default:
    for (j = 0; j < n; j++) {
      info = cudaMemset(dA + dA_IDX_OFFSET(0,j), 0x00, m*sizeof(double));
    }
  }
  
}
