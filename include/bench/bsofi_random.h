/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *  Random p-cyclic matrices generation/output
 */

#ifndef __QMC_H__
#define __QMC_H__

#include <stdlib.h>
#include <time.h>
#include <string.h>

#undef  A_IDX
#undef  hA_IDX
#define A_IDX(i, j)   A[ lda*(j) + i]
#define hA_IDX(i, j) hA[ldha*(j) + i]

inline
void qmc_matr_init(int n, int L, scalar_t A[], int lda)
{
  int i_off, j, k;
  int ibeg_a; 
  int ibeg_b; 
  int jbeg, jend; 

  memset(A, '\0', sizeof(scalar_t)*n*L*lda);
  srand((unsigned int)time(NULL));

  for (k = 0; k < L; k++) {
    ibeg_a = k*n; 
    ibeg_b = ((k+1)%L)*n; 
    jbeg  = k*n; jend = (k+1)*n; 
    for (j = jbeg; j < jend; j++)
      for (i_off = 0; i_off < n; i_off++) {
	A_IDX(ibeg_a + i_off, j) = (scalar_t)rand()/(scalar_t)RAND_MAX;
	A_IDX(ibeg_b + i_off, j) = (scalar_t)rand()/(scalar_t)RAND_MAX;
      }
  }
}

inline
void rand_matr_init(int m, int n, scalar_t A[], int lda)
{
  int i, j; 
  
  srand((unsigned int)time(NULL));
  
  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
      A_IDX(i, j) = (scalar_t)rand()/(scalar_t)RAND_MAX;
}

inline
void zero_matr_init(int m, int n, scalar_t A[], int lda)
{
  int j, i; 
  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
      A_IDX(i, j) = 0;
}

#ifdef __cplusplus
inline
void chess_column_init(int n, int L, scalar_t A[], int lda, bool first_zero = true)
{
  int i;
  bool k = first_zero; 
  scalar_t *M = A;
  
  if(n*L > lda ) return;
  srand((unsigned int)time(NULL));

  for (i = 0; i < L; i++) {
    if(k)rand_matr_init(n, n, M, lda);
    else zero_matr_init(n, n, M, lda);
    k = !k; M += n;
  }
}

inline
void chess_matr_init(int n, int L, scalar_t A[], int lda, bool first_zero = true)
{
  int j;
  bool k = first_zero; 
  scalar_t *M = A;
  
  srand((unsigned int)time(NULL));

  for (j = 0; j < L; j++) {
    chess_column_init(n, L, M, lda, k);
      k = !k; M += lda*L;
  }
}
#endif 

//#include<iostream>
#include<stdio.h>
inline
void print_matrix( int m, int n, scalar_t* A, int lda)
{
  int j, i; 
  printf("\n");
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      printf("%5.2f ", A_IDX(i, j));
    printf("\n");
  }
}

inline
void print_diag( int m, int n, scalar_t* A, int lda)
{
  int i; 
  int mn = m<n?m:n;
  printf("[");
  for (i = 0; i < mn; i++) {
    printf("%5.2f ", A_IDX(i, i));
  }
  printf("]");
}


#undef A_IDX
#undef hA_IDX

#endif
