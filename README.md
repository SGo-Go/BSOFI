==================
BSOF/I README FILE
==================

Structured Orthogonal Inversion Codes for Block p-Cyclic Matrices

=====

This is a set of pure C99 codes for 
block structured factorization and inversion (BSOF/I)
of block p-cyclic matrices on CPU+GPU platforms.

* To COMPILE BSOF/I codes, create the settings.make file to indicate where 
  CUDA, CPU BLAS, MAGMA, and LAPACK are installed on your system. Example
  are given in zwolf.make, showing how to link to 
  MKL, CuBLAS and MAGMA on zWolf platform in UC Davis.  
   
* To TUNE CPU+GPU codes, you can modify the parameters 
  inserted in file 'include/bsofi_params.h'. The default 
  values are tuned for zWolf. 
  Parameters KAPPA_* control workload distribution between host and device.
  Parameters *_SWITCH_* specify when to switch between CPU and hybrid codes. 

* To TEST BSOF/I, compile codes and run 'test'.  

  Note that on some systems the experimental code must be linked against 
  sequential MKL.  If you notice excessively poor performance, change 
  make.inc to link against sequential MKL.

* To  CALL BSOF/I subroutines in own codes, see examples 
  'ftest_hybrid.f90', 'ftest_cpu.f90' (Fortran) and 'ctest.c' (C/C++).  
  The general scheme for Fortran is shown below

  integer  :: nnb, lw
  integer  ::handler(10)

  ! Initialize device and allocate GPU memory
  call hybridXinitH(n, nb, nnb, handler)

  ! Create working space on CPU
  lw = hybridXbsofiLWorkH(n, nb, nnb, -1)
  allocate(W1(lw))

  ! Inversion
  call hybridXbsoftriH(handler, n, nb, A, nnb, tauBsofi, W1, lw, info)
  call hybridXbsoiH(handler, n, nb, A, nnb, tauBsofi, W1, lw, info)

  ! Shutdown device and deallocate GPU memory
  call hybridXfinalizeH(handler)

For more INFORMATION, please refer to the technical report

