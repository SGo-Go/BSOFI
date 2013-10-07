/**********************************************************************
 * BSOF/I - Block structured factorization and inversion codes 
 *  for CPU+GPU platforms
 * Copyright (c) 2013, Sergiy Gogolenko
 * e-mail: sgogolenk@ucdavis.edu
 **********************************************************************
 * Description:
 *  Wraps for use in Fortran
 */

#include "config.h"
#include "bsofi.h"
#include "bsofi_params.h"
#include <string.h>


#ifndef USE_FORTRAN
#  warning "Variable USE_FORTRAN is not set to compile " ## __FILE__
#else 
/************************************************************
 * Callbacks for CPU BSOF/I routines
 ************************************************************/

int cb_cpuXbsofiLWork(int* n, int* L, int* lda)
{
  return lapackXbsofiLWork(*n, *L, *lda);
}
 
int cb_cpuXbsof(int* n, int* L, 
		scalar_t A[], int* lda, scalar_t tau[], 
		scalar_t work[], int* lwork, int* info)
{
  return lapackXbsof(*n, *L, A, *lda, tau, work, *lwork, info);
}

int cb_cpuXbsoi_orgqr(int* n, int* L, 
		      scalar_t A[], int* lda, scalar_t tau[], 
		      scalar_t work[], int* lwork, int* info)
{
  return lapackXbsoi_orgqr(*n, *L, A, *lda, tau, work, *lwork, info);
}

int cb_cpuXbsoi_ormqr(int* n, int* L, 
		      scalar_t A[], int* lda, scalar_t tau[], 
		      scalar_t work[], int* lwork, int* info)
{
  return lapackXbsoi_ormqr(*n, *L, A, *lda, tau, work, *lwork, info);
}

int cb_cpuXbstri(int* n, int* L, 
		 scalar_t A[], int* lda, 
		 // scalar_t work[], int lwork, 
		 int* info)
{
  return lapackXbstri(*n, *L, A, *lda, info);
}


/************************************************************
 * Callbacks for CPU+GPU BSOF/I routines
 ************************************************************/
#ifdef HAS_CUDA
int cb_hybridXbsofiLWork(int* n, int* L, int* lda, int* lworkHostDevice)
{
  return hybridXbsofiLWork(*n, *L, *lda, *lworkHostDevice);
}

int cb_hybridXbsoftri(/* cublasHandle_t handle, */ int* n, int* L, 
		      scalar_t A[], int* lda, scalar_t tau[], 
		      scalar_t work[], int* lwork, scalar_t dwork[], int* ldwork,
		      /* int Lswitch,  */
		      int* info)
{
  cublasHandle_t handle = 0;
  return hybridXbsoftri(handle, *n, *L, A, *lda, tau, work, *lwork, dwork, *ldwork, info);
}

int cb_hybridXbsoi(/* cublasHandle_t handle,  */ int* n, int* L, 
		   scalar_t A[], int* lda, scalar_t tau[], 
		   scalar_t work[], int* lwork, scalar_t dwork[], int* ldwork, 
		   /* int Lswitch,  */
		   int* info)
{
  cublasHandle_t handle = 0;
  return hybridXbsoi(handle, *n, *L, A, *lda, tau, work, *lwork, dwork, *ldwork, info);
}
#endif

/************************************************************
 * Simplified interface for CPU+GPU BSOF/I routines 
 * with Fortran
 ************************************************************/

#define UNPACK_HANDLE							\
  /* handle = malloc(sizeof(bsofiHandle_t)); */			\
  bsofiHandle_t handle; memcpy(&handle, _handle, sizeof(handle))

#define STORE_HANDLE				\
  memcpy(_handle, &handle, sizeof(handle))

int cb_hybridXinitH(/* bsofiHandle_t* handle */int* _handle, int* n, int* L, int* lda)
{
#ifdef HAS_CUDA
  UNPACK_HANDLE;
  handle.switchBsoi    = (*n >= BSOI_SWITCH_BLK    && (*n) * (*L) >= BSOI_SWITCH_SIZE);
  handle.switchBsoftri = (*n >= BSOFTRI_SWITCH_BLK && (*n) * (*L) >= BSOFTRI_SWITCH_SIZE);

  if (handle.switchBsoi || handle.switchBsoftri) {
    if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {
      /* fprintf(stderr, "ERROR: CUBLAS initialization failed\n"); */
      return -1;
    }

    handle.handle = 0;
    handle.ldwork = hybridXbsofiLWork(*n, *L, *lda, /* Lswitch_min, */ -2);
    /* printf("\n!!!!!cb_hybridXinitH = %d\n", handle.ldwork); */

    handle.dwork  = 0; 
    if (cudaSuccess != cudaMalloc( (void**) &handle.dwork, handle.ldwork*sizeof(scalar_t) )) {
      /* fprintf( stderr, "CUDA: GPU device memory allocation failed\n"); */
      return -1;
    }
  }

  STORE_HANDLE;
#endif
  return 0;
}

int cb_hybridXbsoftriH(/* bsofiHandle_t* handle, */int* _handle,  
		       int* n, int* L, scalar_t A[], int* lda, 
		       scalar_t tau[], scalar_t work[], int* lwork, int* info)
{
#ifdef HAS_CUDA
  UNPACK_HANDLE;
  if (handle.switchBsoftri) 
    return hybridXbsoftri(handle.handle, *n, *L, A, *lda, tau, work, *lwork, handle.dwork, handle.ldwork, info);
  else 
#endif
    {
      lapackXbsof(*n, *L, A, *lda, tau, work, *lwork, info);
      return lapackXbstri(*n, *L, A, *lda, info);
    }
}

int cb_hybridXbsoiH(/* bsofiHandle_t* handle, */int* _handle, 
		    int* n, int* L, scalar_t A[], int* lda, 
		    scalar_t tau[], scalar_t work[], int* lwork, int* info)
{
#ifdef HAS_CUDA
  UNPACK_HANDLE;
  if (handle.switchBsoi) 
    return hybridXbsoi(handle.handle, *n, *L, A, *lda, tau, work, *lwork, handle.dwork, handle.ldwork, info);
  else
#endif
    return lapackXbsoi_ormqr(*n, *L, A, *lda, tau, work, *lwork, info);
}

int cb_hybridXfinalizeH(/* bsofiHandle_t* handle*/int* _handle)
{
#ifdef HAS_CUDA
  UNPACK_HANDLE;
  cudaFree (handle.dwork);
  cublasShutdown(); /* cublasDestroy(handle); */
  /* free(handle); */
#endif
  return 0;
}

int cb_hybridXbsofiLWorkH(/* bsofiHandle_t* handle, */int* _handle, int* n, int* L, int* lda)
{
#ifdef HAS_CUDA
  /* int lworkHostDevice = -1; */
  UNPACK_HANDLE;
  if (handle.switchBsoi || handle.switchBsoftri)
    return hybridXbsofiLWork(*n, *L, *lda, -1/* *lworkHostDevice */);
  else
#endif
    return lapackXbsofiLWork(*n, *L, *lda);
}

/************************************************************
 * Callbacks for *nix timers
 ************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
# define fgetwalltime FROUTINE(fgetwalltime, FGETWALLTIME)
# define felapsed     FROUTINE(felapsed, FELAPSED)
  void   fgetwalltime(double tv[2]);
  double felapsed(double t0[2], double t1[2]);
#ifdef __cplusplus
}
#endif

/* #include <timer/sys_timer.h> */
/* void fgetwalltime(double tv[2]) { */
/*   struct timeval t = getwalltime(); */
/*   tv[0] = t.tv_sec; */
/*   tv[1] = t.tv_usec; */
/* } */

#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
void fgetwalltime(double tv[2]) {
  struct timeval t;
  gettimeofday(&t,NULL);
  tv[0] = t.tv_sec;
  tv[1] = t.tv_usec;
}

double felapsed(double t1[2], double t0[2])
{
  return ((double)(t1[0] - t0[0]) + 
	  ((double)(t1[1] - t0[1])/1000000.0));
}

#endif
