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
#include <third_party/cublas.h>

#include <stdlib.h>
#include <time.h>
#include <timer.h>

#include <bench/flops.h>
//#include <bench/bench_macro.h>
#include <bench/bsofi_random.h>
#include "options.h"


/************************************************************
 *                  Benchmarking macro
 ************************************************************/
#undef BENCH_VARS_DEF
#define BENCH_VARS_DEF				\
  timeval_t start, end;				\
  double cpu_perf, cpu_time, cpu_flops

/* #define BENCH_CPU_VARS_SET				\ */
/*   timeval_t start, end;					\ */
/*   double cpu_perf, cpu_time, cpu_flops */

#define BENCH_TUNER(__flops,				\
		    __code_init, __code_do) {		\
    {__code_init;}					\
    cpu_flops = (__flops)/ 1e6;				\
    info = 0;						\
    start = getwalltime();				\
    {__code_do;}					\
    end   = getwalltime();				\
    if (info != 0)					\
      DBGERROR("returned error %d", (int) info);	\
    cpu_time = elapsed(end, start);			\
    cpu_perf = cpu_flops / cpu_time / 1e3;		\
    MESSAGE( TUNE_FORMAT_GEMM_PARAM, L, cpu_time);	\
  }

#define BENCH_CPU(__flops, __code)		\
  BENCH_TUNER(__flops, , __code)

#define BENCH_GPU(__flops, __code)			\
  BENCH_TUNER(__flops,					\
	      {}, {					\
		cudaDeviceSynchronize();		\
		__code;					\
		cudaDeviceSynchronize();})

int process(int n_first, int n_last, int n_step, int L)
{
  BENCH_VARS_DEF;
 
  cublasHandle_t handle = 0;
  int i, info/* , bench_repeats */;

  size_t L_last, L_first, L_step;
  /* size_t n_first, n_step; */
  size_t n;
  size_t N;
  size_t lda, ldda, ldq, lddq;
  size_t lhwork, ldwork;
  
  scalar_t *hwork = 0, *dwork = 0;
  scalar_t *hA = 0/* , *hB = 0 */, *hC = 0, *hQ = 0;
  scalar_t *dA = 0/* , *dB = 0 */, *dC = 0, *dQ = 0;

  /************************************************************
   *             Define matrix sizes for benchmarking
   ************************************************************/
  /* n_step   = 128; */
  if(n_first < NMIN) n_first = NMIN;
  if(n_last  < NMIN) n_last  = NMAX;
  if(n_step  < N_MIN_STEP) n_first = N_MIN_STEP;

  /*************************************************************
   *                        Init device
   *************************************************************/
  CUBLAS_INIT(handle);
  srand((unsigned int)time(NULL));

  /************************************************************
   *             Query for Lapack workspace size
   ************************************************************/
  lhwork   = 2*NMIN*(2*NMIN) + 3*NMIN*(LMAX*NMIN);
  ldwork   = 2*NMIN*GET_LD_DEV(2*NMIN) + 3*NMIN*GET_LD_DEV(LMAX*NMIN);

  /************************************************************
   *                   Allocate memory 
   ************************************************************/

  LAPACK_MALLOC(  hwork,   scalar_t,    lhwork   );
  /* CUBLAS_HOSTALLOC( hwork,   scalar_t,    lhwork ); */
  CUBLAS_DEVALLOC ( dwork,   scalar_t,    ldwork );

  /************************************************************
   *              Title string in tabular output
   ************************************************************/

  /************************************************************
   *  Iterations over deffierent sizes of p-cyclic matrices
   ************************************************************/
  for( n = n_last; n > n_first; n-=n_step ) {
    L_last  = (int)((double)(NMAX*NMAX*NMAX)*LMIN/(double)(n*n*n));
    L_first = 2;
    if((L_last - L_first)/64 > 0) L_step  = L_last/64;
    else              L_step  = 1;

    L = L_last;
    lda  = n*L;
    ldq  = 2*n;
    ldda = GET_LD_DEV(lda);
    lddq = GET_LD_DEV(ldq);

    hA = hwork; 
    hC = hwork +   n*lda; 
    hQ = hwork + 3*n*lda;

    dA = dwork; 
    dC = dwork +   n*ldda; 
    dQ = dwork + 3*n*ldda;

    /********************** Init ********************************/
    for(i = 0; i < lhwork; i++)
      hwork[i] = (scalar_t)rand()/(scalar_t)RAND_MAX;
    
    MESSAGE(TUNE_FORMAT_GEMM_SIZE, n);
    for(L = L_last; L >= L_first; L -= L_step) {
      N = n*L;
      /* MESSAGE("\n%d %d %d\n", L, lda, N);     fflush(stdout); */
      /********************* CPU+GPU bench ************************/
      BENCH_CPU((FLOPS_DGEMM((double)N, (double)2*n, (double)n)),
		lapackXgemm('N', 'T', N, 2*n, n, 1, 
			    hA, lda,  hQ, ldq,  1, hC, lda));

      cublasXlaset ('A',   N, n, hA, lda, dA,  ldda);
      cublasXlaset ('A', 2*n, n, hQ, ldq, dQ,  lddq);
      BENCH_GPU((FLOPS_DGEMM((double)N, (double)2*n, (double)n)),
		cublasXgemm('N', 'T', N, 2*n, n, 1,
			    dA, ldda, dQ, lddq, 1, dC, ldda);
		cudaDeviceSynchronize();
		cublasXlaget('A', N, n, dC + n*ldda, ldda, hC + n*lda, lda));
    }
    MESSAGE( "\n" );
    fflush(stdout);
  }

  /************************************************************
   *                    Memory clean up
   ************************************************************/
  LAPACK_FREE(   hwork  );
  /* CUBLAS_HOSTFREE(  hwork  ); */
  CUBLAS_DEVFREE (  dwork  );

  /*************************************************************
   *                  Shutdown device
   *************************************************************/
  CUBLAS_FINALIZE(handle);
  return 0;
}
