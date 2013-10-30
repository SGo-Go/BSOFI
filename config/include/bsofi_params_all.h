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

#ifndef __BSOFI_PARAMS_ALL_H__
#define __BSOFI_PARAMS_ALL_H__

#define ZWOLF        1
#define ZWOLF06      2
#define DIRAC_FERMI 10
#define DIRAC_TESLA 11

/************************************************************/
#if BSOFI_PLATFORM == ZWOLF
#warning Optimized for ZWolf

#define KAPPA_R(__n) ((__n) < 155 ? (6.8024e-01 + 1.2383e+01/((__n) - (6.7813e+01))) : 8.2227e-01)
#define KAPPA_Q(__n) ((__n) < 181 ? (7.1411e-01 + 1.3100e+01/((__n) - (7.7861e+01))) : 8.4112e-01)
#define L_F(__n) ((__n) > 1266 ? 2 : ((__n) > 193 ? 3 : 4))
#define C_I(__n) 1
#define C_J(__n) 0
#define C_K1(__n) ((__n) > 1601 ? 3 : ((__n) > 1106 ? 4 : ((__n) > 626 ? 5 : ((__n) > 315 ? 6 : 7))))
#define C_K2(__n) ((__n) > 918 ? 2 : ((__n) > 231 ? 3 : 4))

#endif
/************************************************************/
#if BSOFI_PLATFORM == ZWOLF06
#warning Optimized for ZWolf (1 slot)

#define KAPPA_R(__n) ((__n) < 1863 ? (4.2121e-01 + 2.7216e+01/((__n) - (4.1433e+01))) : 4.3615e-01)
#define KAPPA_Q(__n) ((__n) < 1861 ? (4.1808e-01 + 1.7680e+01/((__n) - (4.6115e+01))) : 4.2783e-01)
#define L_F(__n) ((__n) > 143 ? 3 : ((__n) > 77 ? 4 : 5))
#define C_I(__n) ((__n) > 341 ? 0 : 1)
#define C_J(__n) ((__n) > 91 ? 0 : 1)
#define C_K1(__n) ((__n) > 1345 ? 2 : ((__n) > 905 ? 3 : ((__n) > 445 ? 4 : ((__n) > 153 ? 5 : ((__n) > 116 ? 6 : 7)))))
#define C_K2(__n) ((__n) > 1093 ? 1 : ((__n) > 181 ? 2 : ((__n) > 100 ? 3 : 4)))

#endif
/************************************************************/
#if BSOFI_PLATFORM == DIRAC_TESLA
#warning Optimized for Dirac (Tesla GPU)

#define KAPPA_R(__n) ((__n) < 678 ? (-3.8656e+00 + 2.1618e+04/((__n) - (-3.7944e+03))) : 9.6794e-01)
#define KAPPA_Q(__n) ((__n) < 619 ? (6.3164e-01 + 2.4263e+02/((__n) - (-1.1800e+02))) : 9.6085e-01)
#define L_F(__n) ((__n) > 490 ? 3 : 4)
#define C_I(__n) ((__n) > 278 ? 0 : 1)
#define C_J(__n) ((__n) > 1774 ? -1 : 0)
#define C_K1(__n) ((__n) > 1518 ? 2 : ((__n) > 966 ? 3 : ((__n) > 426 ? 4 : ((__n) > 302 ? 5 : ((__n) > 222 ? 6 : ((__n) > 175 ? 7 : ((__n) > 141 ? 8 : 9)))))))
#define C_K2(__n) ((__n) > 1242 ? 1 : ((__n) > 358 ? 2 : ((__n) > 198 ? 3 : 4)))

#endif
/************************************************************/
#if BSOFI_PLATFORM == DIRAC_FERMI
#warning Optimized for Dirac (Fermi GPU)

#define KAPPA_R(__n) ((__n) < 1185 ? (1.8926e-01 + 7.1841e+01/((__n) - (1.4126e+01))) : 2.5062e-01)
#define KAPPA_Q(__n) ((__n) < 800 ? (1.8603e-01 + 4.4045e+01/((__n) - (4.4763e+01))) : 2.4435e-01)
#define L_F(__n) ((__n) > 646 ? 3 : 4)
#define C_I(__n) ((__n) > 274 ? 0 : 1)
#define C_J(__n) ((__n) > 1446 ? -1 : ((__n) > 141 ? 0 : 1))
#define C_K1(__n) ((__n) > 1770 ? 2 : ((__n) > 1010 ? 3 : ((__n) > 462 ? 4 : ((__n) > 330 ? 5 : ((__n) > 256 ? 6 : ((__n) > 197 ? 7 : 8))))))
#define C_K2(__n) ((__n) > 1262 ? 1 : ((__n) > 374 ? 2 : ((__n) > 206 ? 3 : ((__n) > 135 ? 4 : 5))))

#endif

#endif
