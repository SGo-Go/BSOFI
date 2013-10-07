! ------------------------------------------------------------------
!  BSOF/I - Block structured factorization and inversion codes 
!  for CPU+GPU platforms
!  Copyright (c) 2013, Sergiy Gogolenko
!  e-mail: sgogolenk@ucdavis.edu
! ------------------------------------------------------------------
!  Description:
!   Example of use with Fortran (pure CPU codes based on LAPACK 
!   subroutines calls)

subroutine Init_DQMC_matrix(n, nb, A)
  integer,  parameter :: wp = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1
  integer,  parameter :: STDOUT = 6        ! standard output

  integer,  intent(in)  :: n, nb
  real(wp), intent(out) :: A(n*nb,n*nb)

  ! ... Local var ...
  integer  :: nnb, i, j, k

  integer                :: count, nseeds
  integer, allocatable  :: seeds(:)

  ! ... Executable ...
  nnb   = n*nb

  ! making A
  A(1:nnb,1:nnb) = ZERO
  do i = 1, nnb
     A(i,i) = ONE
  end do

  ! Init non-zero values of A (generate random)
  call random_seed
  call random_seed(SIZE = nseeds)
  allocate(seeds(nseeds))
  call system_clock(count)
  seeds = count
  call random_seed(PUT  = seeds)
  deallocate(seeds)

  do k = 0, nb - 1
     do j = k*n + 1, (k+1)*n
        do i = mod(k+1,nb)*n + 1, (mod(k+1,nb)+1)*n
           call random_number(A(i,j))
        end do
     end do
  end do
end subroutine Init_DQMC_matrix

function Get_Error(n, nb, A, Ai) result(err)
  implicit none
  integer,  parameter :: wp = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1

  !double precision :: dnrm2
  real(wp) :: dnrm2
  external    dnrm2

  integer,  intent(in)  :: n, nb
  real(wp), intent(in)  :: A(n*nb,n*nb)
  real(wp), intent(in)  :: Ai(n*nb,n*nb)
  real(wp)              :: err

  real(wp), parameter   :: alpha = 1.0
  real(wp), parameter   :: beta  = 0.0

  ! ... Local var ...
  integer  :: nnb
  integer  :: i, j
  real(wp), pointer     :: AAi(:,:)

  ! ... Executable ...
  nnb   = n*nb

  allocate(AAi(nnb,nnb))

  call dgemm('N', 'N', nnb, nnb, nnb, alpha, A, nnb, Ai, nnb, beta, AAi, nnb);
  ! if (nnb.le.12) then
  !    write(STDOUT,*) "Product M*M^{-1} = "
  !    do, i=1,nnb
  !       write(STDOUT,"(100f6.2)") ( AAi(i,j), j=1,nnb )
  !    enddo
  ! end if

  do i = 1, nnb
     AAi(i,i) = AAi(i,i) - ONE
  enddo

  ! write(STDOUT,*)  "Abs error:", dnrm2(nnb*nnb,AAi,1)
  ! write(STDOUT,*)  "norm(A^{-1}) = ", dnrm2(nnb*nnb,Ai,1)
  ! write(STDOUT,*)  "norm(A)      = ", dnrm2(nnb*nnb,A,1)

  err = dnrm2(nnb*nnb,AAi,1)/(dnrm2(nnb*nnb,A,1)*dnrm2(nnb*nnb,Ai,1))

  deallocate(AAi)
end function Get_Error

subroutine InvertA_BSOFI(n, nb, A, s)
  integer,  parameter :: wp = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1
  integer,  parameter :: STDOUT = 6        ! standard output

  integer,  intent(in)  :: n, nb
  !real(wp), pointer, intent(out) :: A(:,:)
  real(wp), intent(out) :: A(n*nb,n*nb)
  real(wp), intent(out) :: s!sgn

  integer  :: lapackXbsofiLWork
  external    lapackXbsofiLWork

  ! ... Local var ...
  real(wp),   allocatable :: tauBsofi(:)
  real(wp),   allocatable :: W1(:)

  integer  :: nnb, lw
  integer  :: info 
  !real(wp) :: s

  ! ... Executable ...
  nnb   = n*nb

  ! Create working space 
  lw = lapackXbsofiLWork(n, nb, nnb)
  allocate(W1(lw))
  allocate(tauBsofi(nnb))

  ! Inversion
  call lapackXbsof(n, nb, A, nnb, tauBsofi, W1, lw, info)
  s = ONE
  do i = 1, nnb
     if (tauBsofi(i) .eq. ZERO) s = -s
     if (A(i,i) .lt. ZERO) s = -s
  enddo
  call lapackXbstri(n, nb, A, nnb, info)
  call lapackXbsoi_ormqr(n, nb, A, nnb, tauBsofi, W1, lw, info)

  !sgn = s

  deallocate(tauBsofi)
  deallocate(W1)
end subroutine InvertA_BSOFI

subroutine InvertA_GETRI(n, nb, A, s)
  integer,  parameter :: wp = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1
  integer,  parameter :: STDOUT = 6        ! standard output

  integer,           intent(in)  :: n, nb
  real(wp), intent(out) :: A(n*nb,n*nb)
  real(wp), intent(out) :: s!sgn

  integer  :: lapackXbsofiLWork
  external    lapackXbsofiLWork

  ! ... Local var ...
  real(wp),   allocatable :: W1(:)
  integer,    allocatable :: IW(:)

  integer  :: nnb, lw
  integer  :: info 
  !real(wp) :: s
  real(wp) :: query(1)

  ! ... Executable ...
  nnb   = n*nb

  ! Create working space 
  call dgetri(nnb, A, nnb, IW, query, -1, info)
  lw = nint(query(1))
  allocate(IW(nnb))
  allocate(W1(lw))

  ! Inversion
  call dgetrf(nnb, nnb, A, nnb, IW, info)
  s = ONE
  do i = 1, nnb
     if (IW(i).ne. i)  s = -s
     if (A(i,i) .lt. ZERO) s = -s
  enddo
  call dgetri(nnb, A, nnb, IW, W1, lw, info)

  !sgn = s

  deallocate(IW)
  deallocate(W1)
end subroutine InvertA_GETRI

!--------------------------------------------------------------------!

program TestInverse
  implicit none

  integer,  parameter :: wp = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1
  integer,  parameter :: STDOUT = 6        ! standard output

  real(wp) :: Get_Error
  external    Get_Error

  double precision  :: felapsed
  external             felapsed

  integer, parameter    :: n  = 216!128!512!
  integer, parameter    :: nb = 10
  integer, parameter    :: nnb = n*nb

  !real(wp), pointer     :: A(:,:)
  !real(wp), pointer     :: Ai(:,:)
  real(wp), allocatable :: A(:,:)
  real(wp), allocatable :: Ai(:,:)
  !real(wp), dimension(nnb,nnb) :: A
  !real(wp), dimension(nnb,nnb) :: Ai
  !real(wp)              :: A(n*nb,n*nb)
  !real(wp)              :: Ai(n*nb,n*nb)
  real(wp)              :: errGETRI, errBSOFI, sgnGETRI, sgnBSOFI
  integer               :: i,j
  real                  :: t1, t2
  real(wp)              :: t_inv, t_invBSOFI, t_invGETRI

  double precision      :: tim1(2), tim2(2)

  allocate(A(nnb,nnb))
  allocate(Ai(nnb,nnb))

  call Init_DQMC_matrix(n, nb, A)

  Ai(1:nnb,1:nnb) = A(1:nnb,1:nnb)
  call fgetwalltime(tim1) 
  call cpu_time(t1)
  call InvertA_BSOFI(n, nb, Ai, sgnBSOFI)
  call fgetwalltime(tim2) 
  call cpu_time(t2)
  t_invBSOFI = felapsed(tim2, tim1) 
  t_inv = t2 - t1
  errBSOFI = Get_Error(n, nb, A, Ai)

  Ai(1:nnb,1:nnb) = A(1:nnb,1:nnb)
  call fgetwalltime(tim1)
  call InvertA_GETRI(n, nb, Ai, sgnGETRI)
  call fgetwalltime(tim2)
  t_invGETRI = felapsed(tim2, tim1)
  errGETRI = Get_Error(n, nb, A, Ai)

  if (nnb.le.12) then
     write(STDOUT,*) "Matrix M"
     do, i=1,nnb
        write(STDOUT,"(100f6.2)") ( A(i,j), j=1,nnb )
        !print *, ( A(i,j), j=1,nnb )
     enddo

     write(STDOUT,*) "Inversion of M"
     do, i=1,nnb
        write(STDOUT,"(100f6.2)") ( Ai(i,j), j=1,nnb )
     enddo
  end if

  write(STDOUT,"(a,f14.10)") " Determinant sign (BSOFI):", sgnBSOFI
  write(STDOUT,"(a,f14.10)") " Determinant sign (GETRI):", sgnGETRI
  write(STDOUT,"(a,e12.4)") " Relative error (BSOFI):", errBSOFI
  write(STDOUT,"(a,e12.4)") " Relative error (GETRI):", errGETRI
  write(STDOUT,"(a,e12.4,a,e12.4)") " Invert M time  (BSOFI):", t_invBSOFI, " (seconds) CPU_TIME=", t_inv/12
  write(STDOUT,"(a,e12.4,a)") " Invert M time  (GETRI):", t_invGETRI, " (seconds)"

  deallocate(A)
  deallocate(Ai)
end program TestInverse
