/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Benchmark for CPU+GPU implementation of BSOF/I routines
 */

#include "setup.h"
#include <config.h>

#include <bsofi.h>
#include <third_party/lapack.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <timer.h>

#include <bench/bsofi_random.h>
#include <bench/bench_macro.h>
#include <bench/flops.h>

#include "options.h"

#ifndef MATR_INIT
#  define MATR_INIT qmc_matr_init
#endif

int process(int threads, int tests, int n_process, int L_first)
{
  timeval_t start, end;
  double cpu_perf, cpu_time, cpu_flops;

 
  cublasHandle_t handle = 0;
  int info;
#ifdef BENCH_CPU_BSOFI
  scalar_t tmp[1]; 
#endif

  int L_last, L_step;
  int n, L;
  int N, lda, ldda;
  int lwork, ldwork;
  
  scalar_t *work = 0, *dwork = 0;
  scalar_t *A = 0, *tau = 0;
  scalar_t *hA = 0;

#ifdef BENCH_GEMM
  scalar_t *dX, *dW, *dQ;
  int lddq, ldbottom;
#endif

  /************************************************************
   *             Define matrix sizes for benchmarking
   ************************************************************/
  L_first = L_first > 4 ? L_first : 4;
  L_last  = L_first + tests;
  L_step  = 1;
  n       = n_process > 0 ? n_process : 64;

  /*************************************************************
   *                        Init device
   *************************************************************/
#ifdef HAS_CUBLAS
  if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {
    DBGERROR("CUBLAS initialization failed");
    exit(-1);
  }
#endif

  /************************************************************
   *             Query for Lapack workspace size
   ************************************************************/
  L     = L_last;
  N     = n*L;
  lda   = N;
  ldda  = GET_LD_DEV(N);

  lwork    = hybridXbsofiLWork(n, L,  lda, /* Lswitch_max, */ -1);
  ldwork   = hybridXbsofiLWork(n, L, ldda, /* Lswitch_min, */ -2);

#ifdef BENCH_CPU_BSOFI
  lapackXbsof(n, L, A, lda, tau, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
  lapackXbsoi_ormqr(n, L, A, lda, tau, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
  lapackXbsoi_orgqr(n, L, A, lda, tau, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
#endif

  /************************************************************
   *                   Allocate memory 
   ************************************************************/
  DBGERROR("Mem consumption: [A]=%.2fMb work=%.2fMb dwork=%.2fMb", 
	   (N*lda)*(sizeof(scalar_t)/1048576.), 
	   lwork  *(sizeof(scalar_t)/1048576.),
	   ldwork *(sizeof(scalar_t)/1048576.));

  LAPACK_MALLOC(   tau,   scalar_t,    N      );
  LAPACK_MALLOC(  work,   scalar_t,    lwork  );

  LAPACK_MALLOC(     A,   scalar_t,    N * lda);

  CUBLAS_DEVALLOC ( dwork,   scalar_t,    ldwork  );

#ifdef BENCH_CPU_BSOFI
  DBGERROR("Mem consumption: [hA]=%.2fMb", 
	   (N*lda)*(sizeof(scalar_t)/1048576.));

  CUBLAS_HOSTALLOC(    hA,   scalar_t,    N * lda );
#else
  hA = A;
#endif  

#ifdef BENCH_GEMM
  lddq = GET_LD_DEV(2*n); ldbottom = GET_LD_DEV(L*n);

  dQ = dwork;
  dW = dwork + 2*n*lddq;
  dX = dwork + 2*n*lddq + 2*n*ldbottom;
#endif

  /************************************************************
   *              Title string in tabular output
   ************************************************************/
  MESSAGE(FORMAT_TITLE_SIZE"      "      );
#ifdef BENCH_GEMM
  MESSAGE(FORMAT_TITLE_PARAM(GEMM_GPU)"    ");
#endif
  MESSAGE(FORMAT_TITLE_PARAM(BSOFTRI_HYB)"    ");
  MESSAGE(FORMAT_TITLE_PARAM(BSOI_HYB)"     ");
#ifdef BENCH_CPU_BSOFI
  MESSAGE(FORMAT_TITLE_PARAM(BSOF)"     ");
  MESSAGE(FORMAT_TITLE_PARAM(BSTRI)"    ");
  MESSAGE(FORMAT_TITLE_PARAM(BSOI)"     ");
#endif
  MESSAGE("\n");

  /************************************************************
   *  Iterations over deffierent sizes of p-cyclic matrices
   ************************************************************/
  for( L=L_first; L < L_last; L+=L_step ) {
    N = n*L;
    lda = N; ldda = GET_LD_DEV(N);

    /********************** Init ********************************/
    MATR_INIT(n, L, A, lda);
    MESSAGE( FORMAT_SIZE "  ", n, L);

    /********************* bench xGEMM ************************/
#ifdef BENCH_GEMM
    cudaDeviceSynchronize();
    BENCH((FLOPS_DGEMM((double)L*n, (double)2*n, (double)n)),
	  cublasXgemm   ('N', 'T', L*n, 2*n, n, 1, 
			 dW, ldbottom,
			 dQ, lddq, 0, 
			 dX, ldbottom);
	  cudaDeviceSynchronize());
#endif

    /********************* CPU bench ****************************/
#ifdef BENCH_CPU_BSOFI
    lapackXlacpy('A', N, N, A, lda, hA, lda);
    BENCH((FLOPS_BSOF((double)n, (double)L)),
	  lapackXbsof (n, L, hA, lda, tau, work, lwork, &info));

    BENCH((FLOPS_BSTRI((double)n, (double)L)),
    	  lapackXbstri(n, L, hA, lda, &info));

#  ifdef BENCH_BSOI_AUTO
    if(n < 64 || n > 128) {
#  endif
#  if defined(BENCH_BSOI_AUTO) || defined(BENCH_BSOI_ORMQR)
    BENCH((FLOPS_BSOI((double)n, (double)L)),
	  lapackXbsoi_ormqr(n, L, A, lda, tau, work, lwork, &info));
#  endif
#  ifdef BENCH_BSOI_AUTO
    }else {
#  endif
#  if defined(BENCH_BSOI_AUTO) || defined(BENCH_BSOI_ORGQR)
      BENCH((FLOPS_BSOI((double)n, (double)L)),
    	    lapackXbsoi_orgqr(n, L, A, lda, tau, work, lwork, &info));
#  endif
#  ifdef BENCH_BSOI_AUTO
    }
#  endif
#endif

    /********************* CPU+GPU bench ************************/
#ifdef BENCH_CPU_BSOFI
    lapackXlacpy('A', N, N, A, lda, hA, lda);
#endif  
    BENCH((FLOPS_BSOFTRI((double)n, (double)L)),
	  hybridXbsoftri (handle, n, L, hA, lda, tau,
			  work, lwork, dwork, ldwork,/* L_frac,  */ &info));

    BENCH((FLOPS_BSOI((double)n, (double)L)),
	  hybridXbsoi (handle, n, L, A, lda, tau, 
		       work, lwork, dwork, ldwork,/* L_frac,  */ &info));

    MESSAGE( "\n" );
    MSGFLUSH();
  }

  /************************************************************
   *                    Memory clean up
   ************************************************************/
  LAPACK_FREE(   tau   );
  LAPACK_FREE(   work  );

  LAPACK_FREE(    A    );

#ifdef BENCH_CPU_BSOFI
  CUBLAS_HOSTFREE(    hA   );
#endif  
  CUBLAS_DEVFREE (  dwork  );

  /*************************************************************
   *                  Shutdown device
   *************************************************************/
#ifdef HAS_CUBLAS
  cublasShutdown(); /* cublasDestroy(handle); */
#endif
  return 0;
}
