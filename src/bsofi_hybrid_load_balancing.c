#include "bsofi_params.h"
/* #include "stdio.h" */

#define max(a,b) ((a)<(b)? (b):(a))

/** 
 */

#ifndef BSOFI_STATIC_LOAD
int get_lk(int k, int p, double kappa_Q, int ck1, int ck2) /* const  */
{
  k++;
  int lk2, lk1 = (int)((p + k + 2 + ck1)/(1. + kappa_Q) + 0.5);
  /* DBGERROR("lk1 = %d", lk1); */
  if(lk1 > p) lk1 = p;  
  if (lk1 + k <= p) return (lk1);
  else {
    lk2 = (int)((p + 0.5*kappa_Q*(p - k) + 1 + ck2)/(1. + kappa_Q) + 0.5);
    /* DBGERROR("lk2 = %d", lk1); */
    if(lk2 > p) lk2 = p;
    if (lk2 + k > p) return (lk2);
    else             return (lk1);
  }
}

int get_li(int i, int p, double kappa_R, int ci) /* const  */
{
  i++;
  int li = (int)((p - i + 1. + ci)/(1. + kappa_R) + 0.5);
  ci = 0; 
  if (li > p - i - 3) return (p - i - 3);
  else        return (li);
}

int get_lj(int j, int p, double kappa_C, int cj) /* const  */
{
  j++;
  int lj = (int)((j + 1. + cj)/(1. + kappa_C) + 0.5);
  cj = 0;
  if (lj > j - 2) return j - 2;
  else        return lj;
}

#else /* Implementation for staic balancing */

int get_lk(int k, int L, double kappa_Q, int ck1, int ck2) /* const  */
{
  double lk;
  /* if(k*(KAPPA_Q + 2) <= KAPPA_Q*L - 2) */
  /*   lk = (L + k + 2.)/(1. + KAPPA_Q); */
  /* else */
  /*   lk = (L + 0.5*KAPPA_Q*(L - k) + 1.)/(1. + KAPPA_Q); */
  lk = (kappa_Q*L - 2.)/(kappa_Q + 2.) - ck1;
  if(lk < 0) lk = 0; 
  else if (lk > L) lk = L;
  return (int)(lk + 0.5);
}

int get_li(int i, int L, double kappa_R, int ci) /* const  */
{
  double li = (L - i + 1.)/(1. + kappa_R);
  if(li < 0) li = 0; 
  else if (li > L - i - 3) li = L - i - 3;
  return (int)(li + 0.5);
}

#endif
