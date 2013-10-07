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


/**************************************************
 *      Control of workload distribution          *
 **************************************************/
#define KAPPA_C 1.
#define KAPPA_R 1.
#define KAPPA_Q 1.

/* const int KAPPA_C = 1; */
/* const int KAPPA_R = 1; */
/* const int KAPPA_Q = 1; */

/**************************************************
 *      Control of switching between CPU          *
 *              and hybrid codes                  *
 **************************************************/

#define BSOI_SWITCH_SIZE    8192
#define BSOI_SWITCH_BLK     256
#define BSOFTRI_SWITCH_SIZE 10240
#define BSOFTRI_SWITCH_BLK  512

#endif
