/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Benchmark for pure CPU implementation of BSOF/I routines
 */

#include "setup.h"
#include <config.h>

#include <bsofi.h>
#include <third_party/lapack.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <timer.h>

#include "bench_macro.h"
#include "qmc.h"
#include "options.h"
#include <bench/flops.h>

#ifndef MATR_INIT
#  define MATR_INIT qmc_matr_init
#endif

int process(int threads, int tests, int n_process, int L_first)
{
  timeval_t start, end;
  double cpu_perf, cpu_time, cpu_flops;

  int info;
  scalar_t tmp[1]; 

  int L_last, L_step;
  int n, L;
  int N, lda;
  int lwork;
  
  scalar_t *A = 0, *tau = 0;
  scalar_t *hA = 0;
  scalar_t *work = 0;

#ifdef BENCH_DENSE_CODE
#  ifdef BENCH_GETRI
  int *ipiv = 0;
#  endif
#  ifdef BENCH_TRSM
  int i;
  scalar_t *hR = 0;
#  endif
#  if defined(BENCH_TRSM) || defined(BENCH_TRTRI)
  scalar_t *hQ = 0;
#  endif
#endif

  /************************************************************
   *             Define matrix sizes for benchmarking
   ************************************************************/
  L_first = L_first > 4 ? L_first : 4;
  L_last  = L_first + tests;
  L_step  = 1;
  n       = n_process > 0 ? n_process : 64;

  /************************************************************
   *             Query for Lapack workspace size
   ************************************************************/
  L     = L_last;
  N     = n*L;
  lda   = N;

  lapackXbsof(n, L, A, lda, tau, tmp, -1, &info);
  lwork = tmp[0];
  lapackXbsoi_ormqr(n, L, A, lda, tau, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
  lapackXbsoi_orgqr(n, L, A, lda, tau, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
#ifdef BENCH_DENSE_CODE
#  ifdef BENCH_GETRI
  lapackXgetri(N, hA, lda, ipiv, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
#  endif
#  if defined(BENCH_TRSM) || defined(BENCH_TRTRI)
  lapackXgeqrf(N, N, hA, lda, tau, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
  lapackXormqr('R', 'T', N, N, N, hQ, lda, tau, hA, lda, tmp, -1, &info);
  if (lwork < tmp[0]) lwork = tmp[0];
#  endif
#endif

  /************************************************************
   *                   Allocate memory 
   ************************************************************/
  LAPACK_MALLOC(     A,   scalar_t,    N * lda);
  LAPACK_MALLOC(   tau,   scalar_t,    N      );
  LAPACK_MALLOC(    hA,   scalar_t,    N * lda);
  LAPACK_MALLOC(  work,   scalar_t,    lwork  );
#if defined(BENCH_DENSE_CODE) && (defined(BENCH_TRSM) || defined(BENCH_TRTRI))
  LAPACK_MALLOC(    hQ,   scalar_t,    N * lda);
#endif
#if defined(BENCH_GETRI) && defined(BENCH_DENSE_CODE)
  LAPACK_MALLOC(   ipiv,  int,         N      );
#endif

  /************************************************************
   *              Title string in tabular output
   ************************************************************/
  printf(FORMAT_TITLE_SIZE"      "      );
#ifdef BENCH_DENSE_CODE
#  ifdef BENCH_GETRI
  printf(FORMAT_TITLE_PARAM(GETRF)"    ");
  printf(FORMAT_TITLE_PARAM(GETRI)"    ");
#  endif
#  if defined(BENCH_TRSM) || defined(BENCH_TRTRI)
  printf(FORMAT_TITLE_PARAM(GEQRF)"    ");
  printf(FORMAT_TITLE_PARAM(LACPY)"    ");
#    ifndef BENCH_TRSM
  printf(FORMAT_TITLE_PARAM(TRTRI)"    ");
#    else
  printf(FORMAT_TITLE_PARAM(TRSM)"     ");
#    endif
  printf(FORMAT_TITLE_PARAM(ORMQR)"    ");
#  endif
#endif
  printf(FORMAT_TITLE_PARAM(BSOF)"     ");
  printf(FORMAT_TITLE_PARAM(BSTRI)"    ");
  printf(FORMAT_TITLE_PARAM(BSOI)"     ");
  printf("\n");

  /************************************************************
   *  Iterations over deffierent sizes of p-cyclic matrices
   ************************************************************/
  for( L=L_first; L < L_last; L+=L_step ) {
    N = n*L;
    lda   = N; 

    MATR_INIT(n, L, A, lda);
    printf( FORMAT_SIZE "  ", n, L);

#ifdef BENCH_DENSE_CODE
    /************************************************************/
#  ifdef BENCH_GETRI /* BENCH_GETRI */
#      warning Use LUP
    lapackXlacpy('A', N, N, A, lda, hA, lda);
    BENCH((FLOPS_DGETRF((double)N, (double)N)),
    	  lapackXgetrf(N, N, hA, lda, ipiv, &info));

    BENCH((FLOPS_DGETRI((double)N)),
    	  lapackXgetri(N, hA, lda, ipiv, work, lwork, &info));
#  endif
    /************************************************************/
#  if defined(BENCH_TRSM) || defined(BENCH_TRTRI)
    lapackXlacpy('A', N, N, A, lda, hA, lda);
    BENCH((FLOPS_DGEQRF((double)N, (double)N)),
    	  lapackXgeqrf(N, N, hA, lda, tau, work, lwork, &info));

#    ifndef BENCH_TRSM /* BENCH_TRTRI */
#      warning Use QR with TRTRI
    BENCH(N*(N-1)/2,
	  lapackXcpyzero('L', N, N, hA, lda, hQ, lda));

    BENCH((FLOPS_DTRTRI((double)N)),
    	  lapackXtrtri('U', 'N', N, hA, lda, &info));

    BENCH((FLOPS_DORMQR((double)N, (double)N, (double)N, 'R')),
	  lapackXormqr('R', 'T', N, N, N, hQ, lda, tau, hA, lda, work, lwork, &info));

#    else  /* BENCH_TRSM */
#      warning Use QR with TRSM 
    hR = hQ; 
    BENCH(N*N, 
	  lapackXzero('A', N, N, hR, lda);
	  for(i = 0; i < N; i++) *(hR + i*(lda+1)) = 1);

    BENCH((FLOPS_DTRSM('L', (double)N, (double)N)),
    	  lapackXtrsm ('L', 'U', 'N', 'N', N, N, 1, hA, lda, hR, lda));

    BENCH((FLOPS_DORMQR((double)N, (double)N, (double)N, 'R')),
    	  lapackXormqr ('R', 'T', N, N, N, hA, lda, tau, hR, lda, work, lwork, &info));
#    endif
#  endif
#endif

    /************************************************************/

    lapackXlacpy('A', N, N, A, lda, hA, lda);
    BENCH((FLOPS_BSOF((double)n, (double)L)),
	  lapackXbsof (n, L, hA, lda, tau, work, lwork, &info));

    BENCH((FLOPS_BSTRI((double)n, (double)L)),
    	  lapackXbstri(n, L, hA, lda, &info));

#ifdef BENCH_BSOI_AUTO
    if(n < 64 || n > 128) {
#endif
#if defined(BENCH_BSOI_AUTO) || defined(BENCH_BSOI_ORMQR)
    BENCH((FLOPS_BSOI((double)n, (double)L)),
	  lapackXbsoi_ormqr(n, L, A, lda, tau, work, lwork, &info));
#endif
#ifdef BENCH_BSOI_AUTO
    }else {
#endif
#if defined(BENCH_BSOI_AUTO) || defined(BENCH_BSOI_ORGQR)
      BENCH((FLOPS_BSOI((double)n, (double)L)),
    	    lapackXbsoi_orgqr(n, L, A, lda, tau, work, lwork, &info));
#endif
#ifdef BENCH_BSOI_AUTO
    }
#endif

    printf( "\n" );
    fflush(stdout);
  }

  /************************************************************
   *                    Memory clean up
   ************************************************************/
  LAPACK_FREE(    A    );
  LAPACK_FREE(   tau   );
  LAPACK_FREE(    hA   );
#if defined(BENCH_DENSE_CODE) && (defined(BENCH_TRSM) || defined(BENCH_TRTRI))
  LAPACK_FREE(    hQ   );
#endif
  LAPACK_FREE(   work  );
#if defined(BENCH_GETRI) && defined(BENCH_DENSE_CODE)
  LAPACK_FREE(   ipiv  );
#endif

  /* Shutdown */
  return 0;
}
