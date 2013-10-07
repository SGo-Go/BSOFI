/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Macro for estimation of Flops number required by BSOF/I codes
 */

#ifndef __BSOFI_FLOFS_H__
#define __BSOFI_FLOFS_H__

#include <config.h>
#define MagmaLeft  'L'
#define MagmaRight 'R'
#include <bench/flops.h>

#define FLOPS_BSOF(__n, __L) ((((__L)-2)*((2*FLOPS_DORMQR(2*(__n), (__n), (__n), 'L')) \
					  + (FLOPS_DGEQRF(2*(__n), (__n)))) \
			       + (FLOPS_DGEQRF(2*(__n), 2*(__n)))))

#define FLOPS_BSTRI(__n, __L) (FLOPS_DTRTRI(3*(__n)) + FLOPS_DTRMM('R', (__n)*(__L-3), __n) \
			       + ((__L)-3)*(FLOPS_DTRTRI(__n) + 2*FLOPS_DTRMM('L', __n, __n)) \
			       +(((__L)-3)*((__L)+2)/2)*(FLOPS_DGEMM(__n, __n, __n)))

#define FLOPS_BSOI(__n, __L) (((__L)-2)*(FLOPS_DORMQR((__L)*(__n), 2*(__n), (__n), 'R')) \
			      + (FLOPS_DORMQR((__L)*(__n), 2*(__n), 2*(__n), 'R')))

#define FLOPS_BSOFTRI(__n, __L) (FLOPS_BSOF(__n, __L) + FLOPS_BSTRI(__n, __L))

#endif
