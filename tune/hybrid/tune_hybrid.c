/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Tuning CPU+GPU implementation of BSOF/I routines
 */

#include "setup.h"
#include <config.h>

#include <bsofi.h>
#include <third_party/lapack.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <timer.h>

#include <bench/flops.h>
#include <bench/bench_macro.h>
#include <bench/bsofi_random.h>
#include "options.h"


/************************************************************
 *                  Benchmarking macro
 ************************************************************/
#undef BENCH_VARS_DEF
#define BENCH_VARS_DEF				\
  timeval_t start, end;				\
  double cpu_perf, cpu_time, cpu_flops;		\
  int i_bench

/* #define BENCH_CPU_VARS_SET				\ */
/*   timeval_t start, end;					\ */
/*   double cpu_perf, cpu_time, cpu_flops */

#define BENCH_TUNER(__repeats, __flops,			\
		    __code_init, __code_do) {		\
    hA = hwork;						\
    hB = hwork + 2*n*lda;				\
    hC = hwork + 4*n*lda;				\
    dA = dwork;						\
    dB = dwork + 2*n*ldda;				\
    dC = dwork + 4*n*ldda;				\
    __code_init;					\
    cpu_flops = (__flops)/ 1e6;				\
    {__code_do;}					\
    info = 0;						\
    start = getwalltime();				\
    for(i_bench = 0; i < (__repeats); i_bench++) {	\
      hA +=n; hB +=n; hC +=n; dA +=n; dB +=n; dC +=n;	\
      __code_do;					\
    }							\
    end   = getwalltime();				\
    if (info != 0)					\
      DBGERROR("returned error %d", (int) info);	\
    cpu_time = elapsed(end, start)/(__repeats);		\
    cpu_perf = cpu_flops / cpu_time;			\
    printf( TUNE_FORMAT_PARAM,				\
	    cpu_perf, cpu_flops, cpu_time);		\
  }

#define BENCH_CPU(__repeats, __flops, __code)			\
  BENCH_TUNER(__repeats, __flops,				\
	      memcpy(hwork, work, lwork*sizeof(scalar_t)),	\
	      __code)

#define BENCH_GPU(__repeats, __flops, __code)			\
  BENCH_TUNER(__repeats, __flops,					\
	      {								\
		cublasXlaset ('A',   N, 2*n, A, lda, dA,  ldda);	\
		cublasXlaset ('A',   N, 2*n, B, lda, dB,  ldda);	\
		cublasXlaset ('A',   N, 2*n, C, lda, dC,  ldda);	\
		cublasXlaset ('A', 2*n, 2*n, Q, ldq, dQ,  lddq);	\
	      }, {__code; cudaDeviceSynchronize();})


int process(int threads, int tests, int n_last, int L)
{
  BENCH_VARS_DEF;
 
  cublasHandle_t handle = 0;
  int i, info/* , bench_repeats */;

  int n_first, n_step;
  int n;
  int N; 
  int lda, ldda, ldq, lddq;
  int lwork, lhwork, ldwork;
  
  scalar_t *work = 0;
  scalar_t *A = 0, *B = 0, *C = 0, *Q = 0;
  scalar_t *hwork = 0, *dwork = 0;
  scalar_t *hA = 0, *hB = 0, *hC = 0, *hQ = 0;
  scalar_t *dA = 0, *dB = 0, *dC = 0, *dQ = 0;
  scalar_t *tau = 0;

  /************************************************************
   *             Define matrix sizes for benchmarking
   ************************************************************/
  n_step   = 128;
  n_last   = n_last > 0 ? n_last : 1024;
  n_first  = n_last - tests - 1;

  /*************************************************************
   *                        Init device
   *************************************************************/
  CUBLAS_INIT(handle);
  srand((unsigned int)time(NULL));

  /************************************************************
   *             Query for Lapack workspace size
   ************************************************************/
  n = n_last;
  /* L     = 7; */
  N     = MAXN;
  lda   = N;
  ldda  = GET_LD_DEV(N);
  ldq   = 2*n;
  lddq  = GET_LD_DEV(ldq);

  lwork    = 2*n*ldq + 6*n*lda;
  lhwork   = lwork;
  ldwork   = 2*n*lddq + 6*n*ldda;

  /************************************************************
   *                   Allocate memory 
   ************************************************************/
  LAPACK_MALLOC(   tau,   scalar_t,    N        );
  LAPACK_MALLOC(  work,   scalar_t,    lwork    );

  CUBLAS_HOSTALLOC( hwork,   scalar_t,    lhwork );
  CUBLAS_DEVALLOC ( dwork,   scalar_t,    ldwork );

  A = work;
  B = work + 2*n*lda;
  C = work + 4*n*lda;
  Q = work + 6*n*lda;

  hA = hwork; 
  hB = hwork + 2*n*lda; 
  hC = hwork + 4*n*lda; 
  hQ = hwork + 6*n*lda;

  dA = dwork; 
  dB = dwork + 2*n*ldda; 
  dC = dwork + 4*n*ldda; 
  dQ = dwork + 6*n*ldda;

  /************************************************************
   *              Title string in tabular output
   ************************************************************/
  printf(TUNE_FORMAT_TITLE_SIZE);
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_GEMM)"     ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_GEMM2)"    ");
  printf(TUNE_FORMAT_TITLE_PARAM(GPU_GEMM)"     ");
  printf(TUNE_FORMAT_TITLE_PARAM(GPU_GEMM2)"    ");

  printf(TUNE_FORMAT_TITLE_PARAM(CPU_TRMML)"    ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_TRMMR)"    ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_TRTRI)"    ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_TRTRI2)"   ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_GEQRF)"    ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_ORGQR)"    ");
  printf(TUNE_FORMAT_TITLE_PARAM(CPU_ORMQR)"    ");
  printf("\n");

  /************************************************************
   *  Iterations over deffierent sizes of p-cyclic matrices
   ************************************************************/
  for( n = n_last; n > n_first; n-=n_step ) {
    /* N = n*L; */
    /* lda = N; ldda = GET_LD_DEV(N); */

    /********************** Init ********************************/
    for(i = 0; i < lwork; i++)
      work[i] = (scalar_t)rand()/(scalar_t)RAND_MAX;

    L = N/n - 1;
    printf(TUNE_FORMAT_SIZE, n, L);

    /********************* CPU+GPU bench ************************/
    memcpy(hwork, work, lwork*sizeof(scalar_t));
    /* ("   "(CPU_GEMM)"     "); */
    BENCH_CPU(L, (FLOPS_DGEMM((double)n, (double)n, (double)n)),
	      lapackXgemm ('N', 'N', n, n, n, 1, 
			   hA, lda, hB, lda, 0, hC, lda));

    /* ("   "(CPU_GEMM2)"    "); */
    BENCH_CPU(L, (FLOPS_DGEMM((double)n, (double)2*n, (double)n)),
	      lapackXgemm('N', 'T', n, 2*n, n, 1, 
			  hA, lda, hQ + n*ldq, ldq, 1, hC, lda));


    /* ("   "(GPU_GEMM)"     "); */
    BENCH_GPU(L, (FLOPS_DGEMM((double)n, (double)n, (double)n)),
	      cublasXgemm ('N', 'N', n, n, n, 1, 
			   dA, ldda, dB, ldda, 0, dC, ldda);
	      cublasXlaget('A', n, n, dC, ldda, hC, lda));

    /* ("   "(GPU_GEMM2)"    "); */
    BENCH_GPU(L, (FLOPS_DGEMM((double)n, (double)2*n, (double)n)),
    	      cublasXgemm('N', 'T', n, 2*n, n, 1,
    			  dA, ldda, dQ + n*lddq, lddq, 1, dC, ldda);
    	      cublasXlaget('A', n, n, dC + n*ldda, ldda, hA + n*lda, lda));

    /* ("   "(CPU_TRMML)"    "); */
    BENCH_CPU(L, (FLOPS_DTRMM(MagmaLeft, (double)n, (double)n)),
	      lapackXtrmm ('L', 'U', 'N', 'N', n, n, -1, hA, lda, hB, lda));

    /* ("   "(CPU_TRMMR)"    "); */
    BENCH_CPU(L, (FLOPS_DTRMM(MagmaRight, (double)n, (double)n)),
	      lapackXtrmm ('R', 'U', 'N', 'N', n, n, 1, hC, lda, hB, lda)); 

    /* ("   "(CPU_TRTRI)"    "); */
    BENCH_CPU(L, (FLOPS_DTRTRI((double)n)),
	      lapackXtrtri('U', 'N', n, hA, lda, &info));

    /* ("   "(CPU_TRTRI2)"   "); */
    BENCH_CPU(L, (FLOPS_DTRTRI((double)2*n)),
	      lapackXtrtri('U', 'N', 2*n, hA, lda, &info));

    /* ("   "(CPU_GEQRF)"    "); */
    BENCH_CPU(L, (FLOPS_DGEQRF((double)2*n, (double)n)),
	      lapackXgeqrf(2*n, n, hA, lda, tau, hQ, 2*n*ldq, &info));

    /* ("   "(CPU_ORGQR)"    "); */
    BENCH_CPU(L, (FLOPS_DORGQR((double)2*n, (double)2*n, (double)n)),
	      lapackXorgqr  (2*n, 2*n, n, hA, ldq, tau, hQ, 2*n*ldq, &info));

    /* ("   "(CPU_ORMQR)"    "); */
    BENCH_CPU(L, (FLOPS_DORMQR((double)2*n, (double)n, (double)n, MagmaLeft)),
	      lapackXormqr('L', 'T', 2*n, n, n,
			   hA, lda, tau, hB, lda, hQ, 2*n*ldq, &info));


    printf( "\n" );
    fflush(stdout);
  }

  /************************************************************
   *                    Memory clean up
   ************************************************************/
  LAPACK_FREE(   work  );
  LAPACK_FREE(   tau   );

  CUBLAS_HOSTFREE(  hwork  );
  CUBLAS_DEVFREE (  dwork  );

  /*************************************************************
   *                  Shutdown device
   *************************************************************/
  CUBLAS_FINALIZE(handle);
  return 0;
}
