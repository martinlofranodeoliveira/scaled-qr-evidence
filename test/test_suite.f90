program test_suite
  use, intrinsic :: iso_fortran_env, only: real64, int64, output_unit, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use scaled_qr_evidence
  implicit none

  interface
    subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
      import :: real64
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      real(real64), intent(inout) :: a(lda, *), b(ldb, *)
      real(real64), intent(out) :: s(*), work(*)
      real(real64), intent(in) :: rcond
      integer, intent(out) :: rank, iwork(*), info
    end subroutine dgelsd
  end interface

  call test_well_conditioned_and_lapack_oracle()
  call test_badly_scaled_columns()
  call test_rank_deficient_and_zero_columns()
  call test_invalid_inputs_and_nan()
  call test_budget_limits()

  write(output_unit, '(A)') "PASS all scaled-qr-evidence Fortran unit and LAPACK DGELSD oracle tests"

contains

  subroutine solve_lapack_dgelsd(m, n, nrhs, a_in, b_in, x_out, info_out)
    integer, intent(in) :: m, n, nrhs
    real(real64), intent(in) :: a_in(m, n), b_in(m, nrhs)
    real(real64), intent(out) :: x_out(n, nrhs)
    integer, intent(out) :: info_out

    real(real64) :: a_copy(m, n), b_copy(max(m, n), nrhs), s(min(m, n))
    real(real64) :: work_query(1)
    real(real64), allocatable :: work(:)
    integer :: iwork(2048), lwork, rank
    real(real64), parameter :: rcond = 1.0e-12_real64

    a_copy = a_in
    b_copy = 0.0_real64
    b_copy(1:m, 1:nrhs) = b_in

    lwork = -1
    call dgelsd(m, n, nrhs, a_copy, m, b_copy, max(m, n), s, rcond, rank, work_query, lwork, iwork, info_out)
    lwork = max(1024, int(work_query(1)) + 64)
    allocate(work(lwork))
    call dgelsd(m, n, nrhs, a_copy, m, b_copy, max(m, n), s, rcond, rank, work, lwork, iwork, info_out)
    deallocate(work)

    x_out(1:n, 1:nrhs) = b_copy(1:n, 1:nrhs)
  end subroutine solve_lapack_dgelsd

  subroutine test_well_conditioned_and_lapack_oracle()
    integer, parameter :: m = 12, n = 4, nrhs = 3
    real(real64) :: a(m, n), b(m, nrhs), x_lapack(n, nrhs)
    real(real64) :: rel_recon_err, orth_err, rel_sol_diff
    logical :: recon_ok
    integer :: i, j, rhs, info
    type(sqe_budget_t) :: budget
    type(sqe_result_t) :: res

    do j = 1, n
      do i = 1, m
        a(i, j) = cos(real(i * j, real64) * 0.37_real64) + real(i, real64) * 0.05_real64**(j - 1)
        if (i == j) a(i, j) = a(i, j) + 3.0_real64
      end do
    end do

    do rhs = 1, nrhs
      do i = 1, m
        b(i, rhs) = sin(real(i + 2 * rhs, real64) * 0.51_real64) + real(rhs, real64)
      end do
    end do

    call solve_scaled_qr(a, b, m, n, nrhs, 30, budget, res)
    if (res%status /= STATUS_COMPUTED) then
      write(error_unit, '(A,A)') "Expected COMPUTED on well-conditioned fixture, got: ", res%error_code
      stop 1
    end if

    call verify_qr_factorization(a, res, rel_recon_err, orth_err, recon_ok)
    if (.not. recon_ok) then
      write(error_unit, '(A,2ES14.6)') "QR reconstruction/orthogonality exceeded 5e-12: ", rel_recon_err, orth_err
      stop 1
    end if

    call solve_lapack_dgelsd(m, n, nrhs, a, b, x_lapack, info)
    if (info /= 0) then
      write(error_unit, '(A,I0)') "LAPACK DGELSD failed with info=", info
      stop 1
    end if

    do rhs = 1, nrhs
      rel_sol_diff = scaled_vec_norm2(res%x(1:n, rhs) - x_lapack(1:n, rhs)) / &
                     scaled_vec_norm2(x_lapack(1:n, rhs))
      if (rel_sol_diff > 1.0e-8_real64) then
        write(error_unit, '(A,ES14.6)') "Relative difference vs DGELSD exceeded 1e-8: ", rel_sol_diff
        stop 1
      end if
      if (res%normal_eq_rel_residual(rhs) > 5.0e-11_real64) then
        write(error_unit, '(A,ES14.6)') "Normalized first-order residual exceeded 5e-11: ", &
          res%normal_eq_rel_residual(rhs)
        stop 1
      end if
    end do
  end subroutine test_well_conditioned_and_lapack_oracle

  subroutine test_badly_scaled_columns()
    integer, parameter :: m = 8, n = 3, nrhs = 2
    real(real64) :: a(m, n), b(m, nrhs), x_lapack(n, nrhs), x_true(n, nrhs)
    real(real64) :: scales(n), rel_recon_err, orth_err, rel_diff
    logical :: recon_ok
    integer :: i, j, rhs, info
    type(sqe_budget_t) :: budget
    type(sqe_result_t) :: res

    ! Part A: Scaled columns within condition number <= 10000 compared with DGELSD
    scales = [1.0e-2_real64, 1.0_real64, 1.0e2_real64]
    do j = 1, n
      do i = 1, m
        a(i, j) = scales(j) * (1.0_real64 / real(i + j, real64) + merge(2.0_real64, 0.0_real64, i == j))
      end do
    end do
    do rhs = 1, nrhs
      do i = 1, m
        b(i, rhs) = real(i * rhs, real64) * 0.25_real64
      end do
    end do

    call solve_scaled_qr(a, b, m, n, nrhs, 30, budget, res)
    if (res%status /= STATUS_COMPUTED) then
      write(error_unit, '(A)') "Failed moderately scaled column test"
      stop 1
    end if

    call verify_qr_factorization(a, res, rel_recon_err, orth_err, recon_ok)
    if (.not. recon_ok) then
      write(error_unit, '(A)') "Failed QR reconstruction on moderately scaled columns"
      stop 1
    end if

    call solve_lapack_dgelsd(m, n, nrhs, a, b, x_lapack, info)
    do rhs = 1, nrhs
      rel_diff = scaled_vec_norm2(res%x(1:n, rhs) - x_lapack(1:n, rhs)) / &
                 scaled_vec_norm2(x_lapack(1:n, rhs))
      if (rel_diff > 1.0e-8_real64) then
        write(error_unit, '(A,ES14.6)') "Scaled DGELSD comparison failed: ", rel_diff
        stop 1
      end if
    end do

    ! Part B: Extreme column scaling (1e-20 to 1e20) with exact synthetic solution
    scales = [1.0e-20_real64, 1.0_real64, 1.0e20_real64]
    do j = 1, n
      do i = 1, m
        a(i, j) = scales(j) * (1.0_real64 / real(i + j, real64) + merge(2.0_real64, 0.0_real64, i == j))
      end do
    end do
    x_true(:, 1) = [2.0e20_real64, -3.5_real64, 1.25e-20_real64]
    x_true(:, 2) = [-1.5e20_real64, 2.25_real64, -0.5e-20_real64]
    b = matmul(a, x_true)

    call solve_scaled_qr(a, b, m, n, nrhs, 30, budget, res)
    if (res%status /= STATUS_COMPUTED) then
      write(error_unit, '(A)') "Failed extreme column scaling test"
      stop 1
    end if
    call verify_qr_factorization(a, res, rel_recon_err, orth_err, recon_ok)
    if (.not. recon_ok) then
      write(error_unit, '(A)') "Failed QR reconstruction on extreme column scaling"
      stop 1
    end if
    do rhs = 1, nrhs
      rel_diff = scaled_vec_norm2(res%x(1:n, rhs) - x_true(1:n, rhs)) / &
                 scaled_vec_norm2(x_true(1:n, rhs))
      if (rel_diff > 1.0e-8_real64) then
        write(error_unit, '(A,ES14.6)') "Extreme scaling exact recovery failed: ", rel_diff
        stop 1
      end if
    end do
  end subroutine test_badly_scaled_columns

  subroutine test_rank_deficient_and_zero_columns()
    integer, parameter :: m = 5, n = 3, nrhs = 1
    real(real64) :: a(m, n), b(m, nrhs)
    type(sqe_budget_t) :: budget
    type(sqe_result_t) :: res
    integer :: i

    ! Zero column
    a = 1.0_real64
    a(:, 2) = 0.0_real64
    b = 1.0_real64
    call solve_scaled_qr(a, b, m, n, nrhs, 20, budget, res)
    if (res%status /= STATUS_NUMERIC_UNRESOLVED) then
      write(error_unit, '(A)') "Zero column did not return NUMERIC_UNRESOLVED"
      stop 1
    end if

    ! Repeated column
    do i = 1, m
      a(i, 1) = real(i, real64)
      a(i, 2) = real(i * 2, real64)
      a(i, 3) = real(i * i + 1, real64)
    end do
    call solve_scaled_qr(a, b, m, n, nrhs, 20, budget, res)
    if (res%status /= STATUS_NUMERIC_UNRESOLVED) then
      write(error_unit, '(A)') "Repeated column did not return NUMERIC_UNRESOLVED"
      stop 1
    end if
  end subroutine test_rank_deficient_and_zero_columns

  subroutine test_invalid_inputs_and_nan()
    real(real64) :: a(4, 2), b(4, 1)
    type(sqe_budget_t) :: budget
    type(sqe_result_t) :: res

    a = reshape([1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                 0.0_real64, 1.0_real64, 1.0_real64, 0.0_real64], [4, 2])
    b(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]

    ! NaN entry
    a(2, 1) = ieee_value(0.0_real64, ieee_quiet_nan)
    call solve_scaled_qr(a, b, 4, 2, 1, 20, budget, res)
    if (res%status /= STATUS_INVALID_INPUT) then
      write(error_unit, '(A)') "NaN entry did not return INVALID_INPUT"
      stop 1
    end if

    ! Inf entry
    a(2, 1) = ieee_value(0.0_real64, ieee_positive_inf)
    call solve_scaled_qr(a, b, 4, 2, 1, 20, budget, res)
    if (res%status /= STATUS_INVALID_INPUT) then
      write(error_unit, '(A)') "Inf entry did not return INVALID_INPUT"
      stop 1
    end if

    ! Out-of-bounds magnitude (< 2^-100)
    a(2, 1) = 2.0_real64**(-105)
    call solve_scaled_qr(a, b, 4, 2, 1, 20, budget, res)
    if (res%status /= STATUS_INVALID_INPUT) then
      write(error_unit, '(A)') "Sub-threshold magnitude did not return INVALID_INPUT"
      stop 1
    end if

    ! Out-of-bounds threshold_k
    a(2, 1) = 1.0_real64
    call solve_scaled_qr(a, b, 4, 2, 1, 5, budget, res)
    if (res%status /= STATUS_INVALID_INPUT) then
      write(error_unit, '(A)') "Invalid threshold_k did not return INVALID_INPUT"
      stop 1
    end if
  end subroutine test_invalid_inputs_and_nan

  subroutine test_budget_limits()
    real(real64) :: a(6, 3), b(6, 1)
    type(sqe_budget_t) :: tiny_updates, tiny_scratch
    type(sqe_result_t) :: res
    integer :: i, j

    do j = 1, 3
      do i = 1, 6
        a(i, j) = real(i + j * j, real64)
      end do
    end do
    b(:, 1) = 1.0_real64

    tiny_updates%max_updates = 10_int64
    call solve_scaled_qr(a, b, 6, 3, 1, 20, tiny_updates, res)
    if (res%status /= STATUS_UNKNOWN_BUDGET) then
      write(error_unit, '(A)') "Tiny update budget did not return UNKNOWN_BUDGET"
      stop 1
    end if

    tiny_scratch%max_scratch_bytes = 64_int64
    call solve_scaled_qr(a, b, 6, 3, 1, 20, tiny_scratch, res)
    if (res%status /= STATUS_UNKNOWN_BUDGET) then
      write(error_unit, '(A)') "Tiny scratch budget did not return UNKNOWN_BUDGET"
      stop 1
    end if
  end subroutine test_budget_limits

end program test_suite
