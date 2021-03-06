/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *   Parameters for efficient workload distribution between host and 
 *   device (see tech. report)
 */

#ifndef __BSOFI_PARAMS_H__
#define __BSOFI_PARAMS_H__

/* #define BSOFI_PLATFORM DIRAC_FERMI */

#ifdef BSOFI_PLATFORM
#  include <bsofi_params_all.h>
#else
#  include <bsofi_params_default.h>
#endif

int get_lk(int k, int p, double kappa_Q, int ck1, int ck2);
int get_li(int i, int p, double kappa_R, int ci);
int get_lj(int j, int p, double kappa_C, int cj);

#endif
