/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Setup file for hybrid CPU+GPU benchmarks
 */

#ifndef __SETUP_HPP__
#define __SETUP_HPP__

/** Maximal value of OP = L*n^3 in tuning benchmarks */
#define  NMIN  128
#define  LMAX 1024
#define  NMAX 1024
#define  LMIN    2
#define  N_MIN_STEP 32

/** Output control: report or console (console is default) */
#define  MAKE_BENCH_REPORT

#ifdef MAKE_BENCH_REPORT
#  define TUNE_FORMAT_GEMM_SIZE  "%3d "
#  define TUNE_FORMAT_GEMM_PARAM "%d %4e "
#else
#  define TUNE_FORMAT_GEMM_SIZE  "%4d "
#  define TUNE_FORMAT_GEMM_PARAM "%4d %6.2f "
#endif

#include <config.h>
#include <bsofi_macro.h>

#endif
