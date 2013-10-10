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

/** Choice of random p-cyclic matrices type */
#define MATR_INIT qmc_matr_init
/* #define MATR_INIT rand_matr_init */

/** Output control: report or console (console is default) */
#define  MAKE_BENCH_REPORT

/** Control if CPU routines benchmarking */
/* #define  BENCH_CPU_BSOFI */
/* #undef   BENCH_CPU_BSOFI   */

/** Choice of BSOI routine implementation (Attention: Use only one!) */
#define  BENCH_BSOI_ORMQR
/* #define  BENCH_BSOI_ORGQR */
/* #define  BENCH_BSOI_AUTO */

#include <config.h>
#include <bsofi_macro.h>
#endif
