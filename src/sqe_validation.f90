module sqe_validation
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use sqe_types
  implicit none
  private

  public :: validate_least_squares_inputs

contains

  subroutine validate_least_squares_inputs(a, b, m, n, nrhs, threshold_k, budget, &
                                           ok, status, status_name, error_code)
    real(real64), intent(in) :: a(:,:), b(:,:)
    integer, intent(in) :: m, n, nrhs, threshold_k
    type(sqe_budget_t), intent(in) :: budget
    logical, intent(out) :: ok
    integer, intent(out) :: status
    character(len=24), intent(out) :: status_name
    character(len=48), intent(out) :: error_code

    real(real64), parameter :: MIN_MAG = 2.0_real64**(-100)
    real(real64), parameter :: MAX_MAG = 2.0_real64**(100)
    real(real64) :: mag
    integer :: i, j

    ok = .false.
    status = STATUS_INVALID_INPUT
    status_name = "INVALID_INPUT"
    error_code = ""

    if (budget%max_updates <= 0_int64 .or. budget%max_scratch_bytes <= 0_int64 .or. &
        budget%max_output_bytes <= 0_int64) then
      status = STATUS_UNKNOWN_BUDGET
      status_name = "UNKNOWN_BUDGET"
      error_code = "nonpositive_budget"
      return
    end if

    if (n < 1 .or. n > SQE_MAX_N) then
      error_code = "invalid_column_count"
      return
    end if

    if (m < n .or. m > SQE_MAX_M) then
      error_code = "invalid_row_count"
      return
    end if

    if (nrhs < 1 .or. nrhs > SQE_MAX_NRHS) then
      error_code = "invalid_rhs_count"
      return
    end if

    if (threshold_k < SQE_MIN_K .or. threshold_k > SQE_MAX_K) then
      error_code = "invalid_threshold_k"
      return
    end if

    if (size(a, 1) < m .or. size(a, 2) < n .or. size(b, 1) < m .or. size(b, 2) < nrhs) then
      error_code = "array_extent_mismatch"
      return
    end if

    do j = 1, n
      do i = 1, m
        if (.not. ieee_is_finite(a(i, j))) then
          error_code = "nonfinite_matrix_entry"
          return
        end if
        mag = abs(a(i, j))
        if (mag > 0.0_real64) then
          if (mag < MIN_MAG .or. mag > MAX_MAG) then
            error_code = "matrix_magnitude_out_of_bounds"
            return
          end if
        end if
      end do
    end do

    do j = 1, nrhs
      do i = 1, m
        if (.not. ieee_is_finite(b(i, j))) then
          error_code = "nonfinite_rhs_entry"
          return
        end if
        mag = abs(b(i, j))
        if (mag > 0.0_real64) then
          if (mag < MIN_MAG .or. mag > MAX_MAG) then
            error_code = "rhs_magnitude_out_of_bounds"
            return
          end if
        end if
      end do
    end do

    ok = .true.
    status = STATUS_COMPUTED
    status_name = "COMPUTED"
    error_code = "ok"
  end subroutine validate_least_squares_inputs

end module sqe_validation
