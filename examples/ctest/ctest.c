/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *  Example of use in C/C++
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

#define n   1728  /* 1024 1024 1728 512 */
#define nb  10    /* 45 45 16 26 90 */

int main(void)
{
  int nnb = n*nb;
  scalar_t *A, *Ai;
  /* scalar_t errGETRI, errBSOFIlapack, errBSOFImagma, errBSOFIhybrid; */
  timeval_t tim1, tim2;
  double t_inv;
  /* double sgnBSOFIlapack = 0; */

  /*************************************************************
   *
   *************************************************************/
  cublasHandle_t handle = 0;

  scalar_t *tauBsofi/* , s */;
  int info/* , i */;
  scalar_t *hW1 = 0; 
  scalar_t *dW = 0;
  int lw, lhw, ldw;
  /* int Lswitch, Lswitch_min = 4, Lswitch_max = nb/2; */

  /*************************************************************
   * Init devices and space
   *************************************************************/

  CUBLAS_INIT();

  lw       = lapackXbsofiLWork(n, nb, nnb);
  lhw      = hybridXbsofiLWork(n, nb, nnb, /* Lswitch_max, */ -1);
  lhw      = lhw > lw ? lhw : lw;
  ldw      = hybridXbsofiLWork(n, nb, nnb, /* Lswitch_min, */ -2);

  CHECK_PTR(A  = malloc(nnb*nnb*sizeof(scalar_t)));
  CHECK_PTR(Ai = malloc(nnb*nnb*sizeof(scalar_t)));

  CHECK_PTR(hW1      = malloc(lhw*sizeof(scalar_t)));
  CHECK_PTR(tauBsofi = malloc(nnb*sizeof(scalar_t)));

#ifdef HAS_CUBLAS
  CHECK_CUMALLOC(cudaMalloc( (void**) &dW, ldw*sizeof(scalar_t) ));
#else
  CHECK_PTR(dW = malloc(ldw*sizeof(scalar_t)));
#endif
  /*************************************************************
   * Init matrix
   *************************************************************/
  Init_DQMC_matrix(n, nb, A);

  /*************************************************************
   * LAPACK xBSOF/I inversion  
   *************************************************************/
#define BSOFI_PREPAR(__code_init)		\
  __code_init

#define BSOFI_TIMING(__code_do)				\
  lapackXlacpy('A', nnb, nnb, A, nnb, Ai, nnb);		\
  tim1 = getwalltime();					\
  __code_do;						\
  tim2 = getwalltime();					\
  t_inv = elapsed(tim2, tim1);				\
  printf("Measured time : %e (seconds)\n", t_inv);

  /*************************************************************
   * CPU/GPU xBSOF/I inversion  
   *************************************************************/
  printf("Memory requirements: \n"					\
	 "  - matrix [A] = %d x %d (%4.2f Gb)\n"			\
	 "  - workspace CPU: %d (%4.2f Mb)\n"				\
	 "  - workspace GPU: %d (%4.2f Mb)\n"				\
	 , nnb, nnb, (double)nnb*nnb/(double)(1<<(30-4)),
	 lhw, (double)lhw/(double)(1<<(20-4)),
	 ldw, (double)ldw/(double)(1<<(20-4)));
  BSOFI_PREPAR();

  /*************************************************************
   * Factorize A and invert R
   *************************************************************/

  printf("\n===================== BSOFTRI ==============================\n");
  printf("Pure CPU(h) >> ");
  BSOFI_TIMING(lapackXbsof (n, nb, Ai, nnb, tauBsofi, hW1, lw, &info);
	       lapackXbstri(n, nb, Ai, nnb, &info));

  printf("Pure CPU    >> ");
  BSOFI_TIMING(lapackXbsof (n, nb, Ai, nnb, tauBsofi, hW1, lw, &info);
	       lapackXbstri(n, nb, Ai, nnb, &info));
  scalar_t t_cpu = t_inv;

  printf("GPU+CPU(h)  >> ");
  BSOFI_TIMING(hybridXbsoftri (handle, n, nb, Ai, nnb, tauBsofi, 
			       hW1, lhw, dW, ldw, /* L_frac,  */ &info));
  double t_opt = t_inv;
  printf("GPU+CPU     >> ");
  BSOFI_TIMING(hybridXbsoftri (handle, n, nb, Ai, nnb, tauBsofi, 
			       hW1, lhw, dW, ldw, /* L_frac,  */ &info));
  t_opt = t_opt < t_inv? t_opt : t_inv;

  printf("================== PERF IMPROVEMENT ========================\n");
  printf("n x L = %4d x %2d  T = %e  R = %4.2f%%\n", n, nb, t_opt, 100*(t_cpu - t_opt)/t_cpu);

  /*************************************************************
   * Apply Q^T to R^{-1}
   *************************************************************/
  lapackXlacpy('A', nnb, nnb, A, nnb, Ai, nnb);

  printf("\n=====================  BSOI  ===============================\n");
  printf("Pure CPU(h) >> ");
  BSOFI_TIMING(lapackXbsoi (n, nb, Ai, nnb, tauBsofi, hW1, lw, &info));
  t_cpu = t_inv;
  printf("Pure CPU    >> ");
  BSOFI_TIMING(lapackXbsoi (n, nb, Ai, nnb, tauBsofi, hW1, lw, &info));
  t_cpu = t_inv;

  printf("GPU+CPU(h)  >> ");
  BSOFI_TIMING(hybridXbsoi (handle, n, nb, Ai, nnb, tauBsofi, hW1, lhw, dW, ldw, /* Lswitch_opt, */ &info));
  t_opt = t_inv;
  printf("GPU+CPU     >> ");
  BSOFI_TIMING(hybridXbsoi (handle, n, nb, Ai, nnb, tauBsofi, hW1, lhw, dW, ldw, /* Lswitch_opt, */ &info));
  t_opt = t_opt < t_inv? t_opt : t_inv;

  printf("================== PERF IMPROVEMENT ========================\n");
  printf("n x L = %4d x %2d  T = %e  R = %4.2f%%\n", n, nb, t_opt, 100*(t_cpu - t_opt)/t_cpu);


  /*************************************************************
   * Free space and close devices
   *************************************************************/
  free(A);
  free(Ai);

  free(tauBsofi);
  free(hW1);

#ifdef HAS_CUBLAS
  cudaFree (dW);
#else
  free (dW);
#endif

  CUBLAS_FINALIZE();

  return 0;
}
