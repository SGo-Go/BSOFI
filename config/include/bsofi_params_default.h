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

#ifndef __BSOFI_PARAMS_DEFAULT_H__
#define __BSOFI_PARAMS_DEFAULT_H__

#include <config.h>

/**************************************************
 *      Control of workload distribution          *
 **************************************************/
#define KAPPA_R_A  2.04316961e-01  
#define KAPPA_R_B -1.35555688e+03
#define KAPPA_R_C -2.08428702e+03
#define KAPPA_R_D  4.33011686e+02
#define KAPPA_Q_A  2.05985350e-01  
#define KAPPA_Q_B -2.83454838e+02  
#define KAPPA_Q_C -6.66607109e+02   
#define KAPPA_Q_D  2.60855477e+02

#  define KAPPA_R(__n) KAPPA_R_A
#  define KAPPA_Q(__n) KAPPA_Q_A

#  define C_K1(__n) 0

/**************************************************
 *      Control of switching between CPU          *
 *              and hybrid codes                  *
 **************************************************/

#endif
