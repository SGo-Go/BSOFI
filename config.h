/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenko@ucdavis.edu
 **********************************************************************
 * Description:
 *  Configuration file
 */

#ifndef __CONFIG_H__
#define __CONFIG_H__
#pragma once

/* #undef  __SINGLE_PREC__ */
#define USE_FORTRAN
#define F77MANGLING _

/**************************************************
 *      Definitions for portability               *
 **************************************************/

#define FROUTINE(f77name, F77NAME) CONCAT(f77name, F77MANGLING)
#define FCOMMONBLOCK(NAME) CONCAT(NAME, F77MANGLING)

#if defined(__C99__)						\
  || ( defined(__cplusplus) && (__cplusplus > 199711) )		\
  || ( defined(__INTEL_COMPILER) && __INTEL_COMPILER > 1000 )
#  define USERDEF_PRAGMA(directive) _Pragma(#directive)
#elif defined(_MSC_VER)
#  define USERDEF_PRAGMA(directive) __pragma(directive)
#else 
#  define USERDEF_PRAGMA(directive)
/* #warning User defined pragmas are not used */
#endif

#if defined __INTEL_COMPILER
#  define FORCE_INLINE inline
#  define RESTRICT __restrict__
#  define FORCE_VECT_LOOP  USERDEF_PRAGMA(ivdep)/*simd statement*/
#  define FORCE_VECT_LOOP2 USERDEF_PRAGMA(vector always) 
#  define ASSUME_ALIGNED(P) P

#elif defined __GNUC__
#  define FORCE_INLINE inline /* attribute(always_inline) */
#  define RESTRICT __restrict__
#  define FORCE_VECT_LOOP  USERDEF_PRAGMA(vector always)
#  define FORCE_VECT_LOOP2 USERDEF_PRAGMA(vector always)
#  if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)
#    define ASSUME_ALIGNED(P) __builtin_assume_aligned (P, USE_ALIGN)
#  else
#    define ASSUME_ALIGNED(P) P
#  endif

#elif defined _MSC_VER
#  define FORCE_INLINE inline /* __forceinline */
#  define RESTRICT __restrict
#  define FORCE_VECT_LOOP  USERDEF_PRAGMA(ivdep)
#  define FORCE_VECT_LOOP2 USERDEF_PRAGMA(vector always)
#  define ASSUME_ALIGNED(P) P
#endif

#ifdef _MSC_VER
#  define WARNING_PRAGMA(__text) __pragma(message #__text)
#else
#endif

#ifdef HAS_MKL
#  include <mkl_service.h>
#  define SET_NUM_BLAS_THREADS(N) mkl_set_num_threads(N)
#else
#  define SET_NUM_BLAS_THREADS(N) 
#endif

/**************************************************
 *   Macro for processor control                  *
 **************************************************/

#ifdef  __GNUC__
#define DISABLE_SSE_EXCEPTIONS()  {		\
    int aux;					\
    asm(					\
	"stmxcsr   %[aux]           \n\t"	\
	"orl       $32832, %[aux]   \n\t"	\
	"ldmxcsr   %[aux]           \n\t"	\
	: : [aux] "m" (aux));			\
  }
#else
#define DISABLE_SSE_EXCEPTIONS()
#endif

/**************************************************
 *   Definitions for profiling & benchmarking     *
 **************************************************/

#ifdef USE_PROF

#  define PRROUTINE(f77name, F77NAME)			\
  FROUTINE(CONCAT(my, f77name), CONCAT(MY, F77NAME))

#  define PROF(COUNTER, FLOPS, CODE) {		\
    if(!tid) {					\
      PROFILE_BEGIN(fcnt ## COUNTER);}		\
    CODE; if(!tid) {				\
      PROFILE_END(fcnt ## COUNTER);		\
      PROFILE_INC_FLOPS_CNT(COUNTER, FLOPS);}	\
  }
#else
#  define PROF(COUNTER, FLOPS, CODE) CODE
#endif

/**************************************************
 *          Macro for debugging                   *
 **************************************************/

#define DBGPRINTF(__fmt, ...)			\
  fprintf(stderr, "[file %s, line %d]: " __fmt,	\
	  __FILE__, __LINE__, __VA_ARGS__)

/* Valid for GCC, but doestn't work with ICC */
/* #define DBGERROR(__fmt, ...) DBGPRINTF("Error>> " __fmt "\n", __VA_ARGS__) */

#define DBGERROR(...)							\
  fprintf(stderr, "Error[%s, l%d]: "					\
	  __ARGS_FIRST__(__VA_ARGS__) "\n",				\
	  __FILE__, __LINE__ __ARGS_REST__(__VA_ARGS__))

/**************************************************
 *      Precision related definitions             *
 **************************************************/

#ifdef __SINGLE_PREC__
typedef float scalar_t;
#else
typedef double scalar_t;
#endif

/**************************************************
 *      Extra macro                               *
 **************************************************/

/* Concatenation of names */
#define CONCAT(STR1, STR2) CONCAT_NONMACRO(STR1, STR2)
#define CONCAT_NONMACRO(STR1, STR2) STR1 ## STR2

/* Expantion to the first argument in __VA_ARGS__
 * If there's only one argument, expands to nothing.  
 * If there is more than one argument, expands to a comma 
 * followed by everything but the first argument.  
 *
 * This trick works equally well with GCC and ICC
 * New macro: __VA_ARGS_FIRST__ and __VA_ARGS_REST__
 * Attantion: we implement support up to 9 arguments!
 * Based on:
 * http://stackoverflow.com/questions/5588855/standard-alternative-to-gccs-va-args-trick
 */
#define __ARGS_FIRST__(...) __ARGS_GET_FIRST__(__VA_ARGS__, __throwaway)
#define __ARGS_GET_FIRST__(__first, ...) __first

#define __ARGS_REST__(...) __ARGS_GET_REST__(__ARGS_NUM__(__VA_ARGS__), __VA_ARGS__)
#define __ARGS_GET_REST__(__qty, ...) __ARGS_PUT_NAMED__(__qty, __VA_ARGS__)
#define __ARGS_PUT_NAMED__(__qty, ...) __ARGS_PUT_##__qty(__VA_ARGS__)
#define __ARGS_PUT_ONE(__first)
#define __ARGS_PUT_TWOORMORE(first, ...) , __VA_ARGS__
#define __ARGS_NUM__(...)						\
  __ARGS_SELECT_10TH__(__VA_ARGS__, TWOORMORE, TWOORMORE, TWOORMORE, TWOORMORE, \
		       TWOORMORE, TWOORMORE, TWOORMORE, TWOORMORE, ONE, throwaway)
#define __ARGS_SELECT_10TH__(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, ...) a10

#endif
