/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Some macro to make easier writting benchmarks.
 */

#ifndef TESTINGS_H
#define TESTINGS_H

/************************************************************
 *                 General purpose macro 
 ************************************************************/

#ifndef min
#define min(a,b)  (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif

/************************************************************
 *                 Output formatting 
 ************************************************************/

#ifdef MAKE_BENCH_REPORT
#  define FORMAT_TITLE_SIZE  " n   L "
#  define FORMAT_TITLE_PARAM(__parname) #__parname "_G " #__parname "_I " #__parname "_T "
#  define FORMAT_SIZE  "%3d %02d "
#  define FORMAT_PARAM "%4e %4e %4e "
#else
#  define FORMAT_TITLE_SIZE  " [nxL] "
#  define FORMAT_TITLE_PARAM(__parname) #__parname ", GFlops (MFlops : ms)"
#  define FORMAT_SIZE  "%3dx%02d   "
#  define FORMAT_PARAM "%6.2f(%10.2f:%9.2f)   "
#endif

/************************************************************
 *                  Benchmarking macro
 ************************************************************/
#define BENCH(__FLOPS, __CODE) {			\
    cpu_flops = (__FLOPS)/ 1e6;				\
    start = getwalltime();				\
    {__CODE;}						\
    end   = getwalltime();				\
    if (info != 0)					\
      DBGPRINTF("returned error %d\n", (int) info);	\
    cpu_time = elapsed(end, start);			\
    cpu_perf = cpu_flops / cpu_time;			\
    printf( FORMAT_PARAM,				\
	    cpu_perf, cpu_flops, cpu_time);		\
  }

/************************************************************
 *  Workspace allocation/deallocation macro for LAPACK codes
 ************************************************************/
#define LAPACK_MALLOC( ptr, type, size )			\
  if ( 0 == (ptr = malloc((size)*sizeof(type)))) {		\
    DBGERROR("malloc failed for: %s\n", #ptr );			\
    exit(-1);							\
  }

#define LAPACK_FREE(ptr)			\
  free(ptr);

/************************************************************
 *  Workspace allocation/deallocation macro for CUBLAS codes
 ************************************************************/
#ifdef HAS_CUBLAS
#  define CUBLAS_HOSTALLOC( ptr, type, size )				\
  if ( cudaSuccess !=							\
       cudaMallocHost( (void**) &ptr, (size)*sizeof(type) )) {		\
    DBGERROR("CUDA pinned malloc failed for: %s\n", #ptr );	\
    exit(-1);								\
  }

#  define CUBLAS_HOSTFREE(ptr)			\
  cudaFreeHost( ptr );

#  define CUBLAS_DEVALLOC( ptr, type, size )				\
  if ( cudaSuccess !=							\
       cudaMalloc( (void**) &ptr, (size)*sizeof(type) )) {		\
    DBGERROR("CUDA device malloc failed for: %s\n", #ptr );		\
    exit(-1);								\
  }

#  define CUBLAS_DEVFREE(ptr)			\
  cudaFree( ptr );

#else
/* #  define CUBLAS_HOSTALLOC(__ptr, __type, __size) LAPACK_MALLOC( __ptr, __type, __size ) */
/* #  define CUBLAS_HOSTFREE (__ptr)             LAPACK_FREE(__ptr) */
/* #  define CUBLAS_DEVALLOC (__ptr, __type, __size) LAPACK_MALLOC( __ptr, __type, __size ) */
/* #  define CUBLAS_DEVFREE  (__ptr)             LAPACK_FREE(__ptr) */
#  define CUBLAS_HOSTALLOC LAPACK_MALLOC
#  define CUBLAS_HOSTFREE  LAPACK_FREE
#  define CUBLAS_DEVALLOC  LAPACK_MALLOC
#  define CUBLAS_DEVFREE   LAPACK_FREE
#endif 

#endif /* TESTINGS_H */
