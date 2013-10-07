/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *  Example of use BSOF/I routines based on Magma calls, in C/C++
 */

#include "bsofi.h"
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <timer.h>

#include "bsofi_macro.h"

scalar_t lapackXnrm2(int n, scalar_t x[], int incx);
scalar_t lapackXlange(char norm, int m, int n, 
		      scalar_t A[], int lda, scalar_t work[]);

int lapackXgetrf(int m, int n, scalar_t A[], int lda,
		 int ipiv[], int *info);
int lapackXgetri(int n, scalar_t A[], int lda, int ipiv[],
		 scalar_t work[], int lwork, int *info);
void lapackXgemm(const char transa, const char transb, 
		 int m, int n, int k, 
		 scalar_t alpha, scalar_t a[], int lda, 
		 scalar_t b[], int ldb, scalar_t beta, 
		 scalar_t c[], int ldc);
void lapackXlacpy(const char uplo, int m, int n,
		  scalar_t A[], int lda, scalar_t B[], int ldb);

void Init_DQMC_matrix(int n, int nb, scalar_t A[])
{
  int nnb=n*nb;
  int i, j, k;

  for(j=0; j < nnb; j++)
    for(i=0; i < nnb; i++)
      A[j*nnb + i] = (i == j) ? 1. : 0.;

  srand((unsigned int)time(NULL));
  for(k=0; k < nb; k++)
    for(j=k*n; j < (k+1)*n; j++)
      for(i=((k+1)%nb)*n; i < ((k+1)%nb + 1)*n; i++) {
	A[j*nnb + i] = (scalar_t)rand()/(scalar_t)RAND_MAX;
      }
}

#define n    1024
#define nb   10

int main(void)
{
  int nnb = n*nb;
  scalar_t *A, *dA = 0;
  scalar_t errGETRI, errBSOFImagma;
  timeval_t tim1, tim2;
  double t_inv;

  /*************************************************************
   *
   *************************************************************/
  scalar_t *tauBsofi = 0, s;
  int info, i;
  scalar_t *hW = 0, *dW = 0, *dT = 0;
  int lhw = 1, ldw = 1, ldT = 1, 
    lda = nnb, ldda = GET_LD_DEV(nnb);

  scalar_t tmp[1];

  /*************************************************************
   * Init devices and space
   *************************************************************/
  if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {
    fprintf(stderr, "ERROR: CUBLAS initialization failed\n");
    exit(-1);
  }

  gpuXbsof(n, nb, dA, ldda, tauBsofi, tmp, -1, dT, &info);
  if (lhw < tmp[0]) lhw = tmp[0];
  gpuXbsoi(n, nb, dA, ldda, tauBsofi, tmp, -1, dT, dW, 0, &info);
  if (lhw < tmp[0]) lhw = tmp[0];
  
  /* ldw = (2*nnb + ldda)*(magma_get_dgeqrf_nb(nnb)); */
  gpuXbsof(n, nb, dA, ldda, tauBsofi, tmp, -2, dT, &info);
  if (ldT < tmp[0]) ldT = tmp[0];

  gpuXbsoi(n, nb, dA, ldda, tauBsofi, tmp, -2, dT, dW, -1, &info);
  if (ldw < tmp[0]) ldw = tmp[0];

  CHECK_PTR(A  = malloc(nnb*lda*sizeof(scalar_t)));
  /* CHECK_PTR(Ai = malloc(nnb*ldda*sizeof(scalar_t))); */

  CHECK_PTR(hW       = malloc(lhw*sizeof(scalar_t)));
  CHECK_PTR(tauBsofi = malloc(nnb*sizeof(scalar_t)));


#define CU_CHECK_MALLOC(__code)						\
  if (cudaSuccess != __code) {						\
    fprintf( stderr, "CUDA: GPU device memory allocation failed\n");	\
    exit(-1);								\
  }

  CU_CHECK_MALLOC(cudaMalloc( (void**) &dA, nnb*ldda*sizeof(scalar_t) ));
  CU_CHECK_MALLOC(cudaMalloc( (void**) &dW, ldw*sizeof(scalar_t) ));
  CU_CHECK_MALLOC(cudaMalloc( (void**) &dT, ldT*sizeof(scalar_t) ));

  /*************************************************************
   * Init matrix
   *************************************************************/
  Init_DQMC_matrix(n, nb, A);

  /*************************************************************
   * MAGMA xBSOF/I inversion  
   *************************************************************/
#define BSOFI_PREPAR(__code_init)		\
  __code_init

#define BSOFI_TIMING(__code_do)				\
  tim1 = getwalltime();					\
  __code_do;						\
  tim2 = getwalltime();					\
  t_inv = elapsed(tim2, tim1);				\
  printf("Measured time : %e (seconds)\n", t_inv);

  /*************************************************************
   * GPU xBSOF/I inversion  
   *************************************************************/
  printf("\n===================== TEST =================================\n");
  printf("Parameters: n x L = %4d x %2d\n", n, nb);

  printf("Memory requirements: \n"					\
	 "  - matrix [A] = %d x %d (%4.2f Gb)\n"			\
	 "  - workspace CPU: %d (%4.2f Mb)\n"				\
	 "  - workspace GPU: %d+%d (%4.2f Mb)\n"				\
	 , nnb, nnb, (double)nnb*nnb/(double)(1<<(30-4)),
	 lhw, (double)lhw/(double)(1<<(20-4)),
	 ldw, ldT, (double)(ldw+ldT)/(double)(1<<(20-4)));

  magma_dsetmatrix( nnb, nnb, A, lda, dA, ldda);	\

  /*************************************************************
   * Factorize A and invert R
   *************************************************************/
  printf("BSOF  >> ");
  BSOFI_TIMING(gpuXbsof(n, nb, dA, ldda, tauBsofi, hW, lhw, dT, &info));

  /*************************************************************
   * Invert R
   *************************************************************/
  printf("BSTRI >> ");
  BSOFI_TIMING(gpuXbstri (n, nb, dA, ldda, dW, ldw, &info));

  /*************************************************************
   * Apply Q^T to R^{-1}
   *************************************************************/
  printf("BSOI  >> ");
  BSOFI_TIMING(gpuXbsoi(n, nb, dA, ldda, tauBsofi, hW, lhw, dT, dW, ldw, &info));

  /*************************************************************
   * Free space and close devices
   *************************************************************/
  /* magma_dgetmatrix( nnb, nnb, dA, ldda, A, lda ); */

  free(A);

  free(tauBsofi);
  free(hW);

  cudaFree (dW);
  cudaFree (dT);
  cudaFree (dA);
  cublasShutdown(); /* cublasDestroy(handle); */

  return 0;
}
