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

#ifndef __BSOFI_PARAMS_DIRAC_FERMI_H__
#define __BSOFI_PARAMS_DIRAC_FERMI_H__

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

#ifdef USE_PMODEL
/* #  define GET_KAPPA(__n, __kappa)					\ */
/*   CONCAT(__kappa, A) + CONCAT(__kappa, B)/(CONCAT(__kappa, C)*(__n)/100.+CONCAT(__kappa, D)) */

#  define GET_KAPPA_R(__n)				\
  KAPPA_R_A + KAPPA_R_B/(KAPPA_R_C*(__n)/100.+KAPPA_R_D)
#  define GET_KAPPA_Q(__n)				\
  KAPPA_Q_A + KAPPA_Q_B/(KAPPA_Q_C*(__n)/100.+KAPPA_Q_D)

#  define GET_Ck1(__n) (((__n) > 512) ? 2 : 3)

#else
/* #  define GET_KAPPA(__n, __kappa) CONCAT(__kappa, A) */

#define GET_KAPPA_R(__n) KAPPA_R_A
#define GET_KAPPA_Q(__n) KAPPA_Q_A

#  define GET_Ck1(__n) 0
#endif

/* #define GET_KAPPA_R(__n) (GET_KAPPA(__n, KAPPA_R_)) */
/* #define GET_KAPPA_Q(__n) (GET_KAPPA(__n, KAPPA_Q_)) */

/**************************************************
 *      Control of switching between CPU          *
 *              and hybrid codes                  *
 **************************************************/
#define BSOI_SWITCH_SIZE    8192
#define BSOI_SWITCH_BLK     256
#define BSOFTRI_SWITCH_SIZE 10240
#define BSOFTRI_SWITCH_BLK  512

#endif
