/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Convenient wraps for Magma routines
 */

#ifndef __MAGMA_HPP__
#define __MAGMA_HPP__

//#include <config.h>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#if defined(HAS_MAGMA)
#  include <magma.h>
#ifdef __SINGLE_PREC__
#  define MAGMA_ROUTINE(f77name, F77name, F77NAME, SUFFIX) CONCAT(CONCAT(magma_s, f77name), SUFFIX)
#else
#  define MAGMA_ROUTINE(f77name, F77name, F77NAME, SUFFIX) CONCAT(CONCAT(magma_d, f77name), SUFFIX)
#endif

#elif defined(HAS_LAPACK)
#  include <third_party/lapack.hpp>
#  define MAGMA_ROUTINE(f77name, F77name, F77NAME, SUFFIX) CONCAT(lapack_x, F77name)

#else 
#  error "Platform lacks for LAPACK and/or Magma installation"
#endif


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#define magma_xGemm  MAGMA_ROUTINE(gemm,  Gemm,  GEMM, )

#define magma_xGeqrf MAGMA_ROUTINE(geqrf, Geqrf, GEQRF, )
#define magma_xOrmqr MAGMA_ROUTINE(ormqr, Ormqr, ORMQR, )
#define magma_xOrgqr MAGMA_ROUTINE(orgqr, Orgqr, ORGQR, )
#define magma_xGeqp3 MAGMA_ROUTINE(geqp3, Geqp3, GEQP3, )

#define magma_xTrsm  MAGMA_ROUTINE(trsm,  Trsm,  TRSM, )
#define magma_xLacpy MAGMA_ROUTINE(lacpy, Lacpy, LACPY, )
#define magma_xGetrf MAGMA_ROUTINE(getrf, Getrf, GETRF, )
#define magma_xGetri MAGMA_ROUTINE(getri, Getri, GETRI, _gpu)

#endif
