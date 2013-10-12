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
#define BENCH_VARS_DEF					\
  timeval_t start, end;					\
  double cpu_perf, cpu_time, cpu_flops

#define BENCH(__FLOPS, __CODE) {			\
    cpu_flops = (__FLOPS)/ 1e6;				\
    start = getwalltime();				\
    {__CODE;}						\
    end   = getwalltime();				\
    if (info != 0)					\
      DBGERROR("returned error %d", (int) info);	\
    cpu_time = elapsed(end, start);			\
    cpu_perf = cpu_flops / cpu_time;			\
    printf( FORMAT_PARAM,				\
	    cpu_perf, cpu_flops, cpu_time);		\
  }

#endif /* TESTINGS_H */
