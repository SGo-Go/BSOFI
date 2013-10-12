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

/** Maximal value of N = L*n in tuning benchmarks */
#define  MAXN (5*1024)

/** Output control: report or console (console is default) */
#define  MAKE_BENCH_REPORT

#ifdef MAKE_BENCH_REPORT
#  define TUNE_FORMAT_TITLE_SIZE  " n   L "
#  define TUNE_FORMAT_TITLE_PARAM(__parname) #__parname "_G " #__parname "_I " #__parname "_T "
#  define TUNE_FORMAT_SIZE  "%3d %02d "
#  define TUNE_FORMAT_PARAM "%4e %4e %4e "
#else
#  define TUNE_FORMAT_TITLE_SIZE  " n(reps) "
#  define TUNE_FORMAT_TITLE_PARAM(__parname) #__parname ", GFlops (MFlops : ms)"
#  define TUNE_FORMAT_SIZE  "%4d(%2d) "
#  define TUNE_FORMAT_PARAM "%6.2f(%10.2f:%9.2f)   "
#endif

#include <config.h>
#include <bsofi_macro.h>

#endif
