module sqe_solver
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use sqe_types
  use sqe_validation
  implicit none
  private

  public :: solve_scaled_qr
  public :: scaled_vec_norm2

contains

  pure function scaled_vec_norm2(v) result(nrm)
    real(real64), intent(in) :: v(:)
    real(real64) :: nrm
    real(real64) :: scale, ssq, absxi
    integer :: i

    scale = 0.0_real64
    ssq = 1.0_real64
    do i = 1, size(v)
      if (v(i) /= 0.0_real64) then
        absxi = abs(v(i))
        if (scale < absxi) then
          ssq = 1.0_real64 + ssq * (scale / absxi)**2
          scale = absxi
        else
          ssq = ssq + (absxi / scale)**2
        end if
      end if
    end do
    nrm = scale * sqrt(ssq)
  end function scaled_vec_norm2

  subroutine solve_upper_triangular(r, n, rhs, sol)
    real(real64), intent(in) :: r(SQE_MAX_N, SQE_MAX_N)
    integer, intent(in) :: n
    real(real64), intent(in) :: rhs(SQE_MAX_N)
    real(real64), intent(out) :: sol(SQE_MAX_N)
    integer :: i, j
    real(real64) :: s

    sol = 0.0_real64
    do i = n, 1, -1
      s = rhs(i)
      do j = i + 1, n
        s = s - r(i, j) * sol(j)
      end do
      sol(i) = s / r(i, i)
    end do
  end subroutine solve_upper_triangular

  subroutine solve_lower_triangular_transpose(r, n, rhs, sol)
    real(real64), intent(in) :: r(SQE_MAX_N, SQE_MAX_N)
    integer, intent(in) :: n
    real(real64), intent(in) :: rhs(SQE_MAX_N)
    real(real64), intent(out) :: sol(SQE_MAX_N)
    integer :: i, j
    real(real64) :: s

    sol = 0.0_real64
    do i = 1, n
      s = rhs(i)
      do j = 1, i - 1
        s = s - r(j, i) * sol(j)
      end do
      sol(i) = s / r(i, i)
    end do
  end subroutine solve_lower_triangular_transpose

  subroutine solve_scaled_qr(a, b, m, n, nrhs, threshold_k, budget, res)
    real(real64), intent(in) :: a(:,:), b(:,:)
    integer, intent(in) :: m, n, nrhs, threshold_k
    type(sqe_budget_t), intent(in) :: budget
    type(sqe_result_t), intent(out) :: res

    logical :: valid
    integer :: st
    character(len=24) :: st_name
    character(len=48) :: err_code

    real(real64) :: a_scaled(SQE_MAX_M, SQE_MAX_N)
    real(real64) :: a_work(SQE_MAX_M, SQE_MAX_N)
    real(real64) :: tau(SQE_MAX_N)
    real(real64) :: v(SQE_MAX_M)
    real(real64) :: col_nrm, max_nrm, beta, x1, v_norm, dot_val, thresh
    real(real64) :: a_frob, r_norm, x_norm, atr_norm, denom, eta, best_g_norm, cand_g_norm
    real(real64) :: qtb(SQE_MAX_M), y_cur(SQE_MAX_N), y_cand(SQE_MAX_N)
    real(real64) :: r_vec(SQE_MAX_M), g_vec(SQE_MAX_N), z_vec(SQE_MAX_N), dy_vec(SQE_MAX_N)
    real(real64) :: atr_vec(SQE_MAX_N), tmp_col(SQE_MAX_M), tmp_scale
    integer :: i, j, k, p, rhs_idx, iter, tmp_piv
    integer(int64) :: updates, scratch_bytes

    res%m = m
    res%n = n
    res%nrhs = nrhs
    res%threshold_k = threshold_k

    call validate_least_squares_inputs(a, b, m, n, nrhs, threshold_k, budget, &
                                       valid, st, st_name, err_code)
    if (.not. valid) then
      res%status = st
      res%status_name = st_name
      res%error_code = err_code
      return
    end if

    ! Account scratch storage in bytes (8 bytes per real64 element)
    scratch_bytes = int(m * n * 2 + m * m + n * n + m * 4 + n * 8, int64) * 8_int64
    res%scratch_bytes_used = scratch_bytes
    if (scratch_bytes > budget%max_scratch_bytes) then
      res%status = STATUS_UNKNOWN_BUDGET
      res%status_name = "UNKNOWN_BUDGET"
      res%error_code = "scratch_storage_budget_exceeded"
      return
    end if

    updates = 0_int64
    a_scaled = 0.0_real64
    a_work = 0.0_real64

    ! Step 1: Column scaling (equilibration)
    do j = 1, n
      res%pivots(j) = j
      updates = updates + int(m, int64)
      if (updates > budget%max_updates) then
        res%status = STATUS_UNKNOWN_BUDGET
        res%status_name = "UNKNOWN_BUDGET"
        res%error_code = "updates_budget_exceeded"
        res%updates_used = updates
        return
      end if
      col_nrm = scaled_vec_norm2(a(1:m, j))
      if (col_nrm <= 0.0_real64 .or. .not. ieee_is_finite(col_nrm)) then
        res%status = STATUS_NUMERIC_UNRESOLVED
        res%status_name = "NUMERIC_UNRESOLVED"
        res%error_code = "zero_or_invalid_column_norm"
        res%updates_used = updates
        return
      end if
      res%col_scales(j) = 1.0_real64 / col_nrm
      do i = 1, m
        a_scaled(i, j) = a(i, j) * res%col_scales(j)
        a_work(i, j) = a_scaled(i, j)
      end do
      updates = updates + int(m, int64)
    end do

    ! Step 2: Column-pivoted Householder QR factorization
    tau = 0.0_real64
    thresh = 2.0_real64**(-threshold_k)

    do j = 1, n
      ! Find pivot column in j..n with maximum residual 2-norm
      p = j
      max_nrm = -1.0_real64
      do k = j, n
        updates = updates + int(m - j + 1, int64)
        col_nrm = scaled_vec_norm2(a_work(j:m, k))
        if (col_nrm > max_nrm) then
          max_nrm = col_nrm
          p = k
        end if
      end do

      if (updates > budget%max_updates) then
        res%status = STATUS_UNKNOWN_BUDGET
        res%status_name = "UNKNOWN_BUDGET"
        res%error_code = "updates_budget_exceeded"
        res%updates_used = updates
        return
      end if

      if (p /= j) then
        tmp_col(1:m) = a_work(1:m, j)
        a_work(1:m, j) = a_work(1:m, p)
        a_work(1:m, p) = tmp_col(1:m)

        tmp_col(1:m) = a_scaled(1:m, j)
        a_scaled(1:m, j) = a_scaled(1:m, p)
        a_scaled(1:m, p) = tmp_col(1:m)

        tmp_piv = res%pivots(j)
        res%pivots(j) = res%pivots(p)
        res%pivots(p) = tmp_piv

        tmp_scale = res%col_scales(j)
        res%col_scales(j) = res%col_scales(p)
        res%col_scales(p) = tmp_scale
      end if

      ! Compute Householder vector v(j:m) and tau(j)
      col_nrm = scaled_vec_norm2(a_work(j:m, j))
      if (col_nrm == 0.0_real64) then
        res%status = STATUS_NUMERIC_UNRESOLVED
        res%status_name = "NUMERIC_UNRESOLVED"
        res%error_code = "zero_pivot_encountered"
        res%updates_used = updates
        return
      end if

      x1 = a_work(j, j)
      if (x1 >= 0.0_real64) then
        beta = -col_nrm
      else
        beta = col_nrm
      end if

      v = 0.0_real64
      v(j) = x1 - beta
      if (j < m) then
        v(j+1:m) = a_work(j+1:m, j)
      end if
      v_norm = scaled_vec_norm2(v(j:m))
      if (v_norm > 0.0_real64) then
        v(j:m) = v(j:m) / v_norm
        tau(j) = 2.0_real64
      else
        tau(j) = 0.0_real64
      end if

      a_work(j, j) = beta
      if (j < m) then
        a_work(j+1:m, j) = v(j+1:m)
      end if

      ! Rank screening check
      if (j == 1) then
        if (abs(a_work(1, 1)) <= 0.0_real64) then
          res%status = STATUS_NUMERIC_UNRESOLVED
          res%status_name = "NUMERIC_UNRESOLVED"
          res%error_code = "rank_screen_zero_leading_pivot"
          res%updates_used = updates
          return
        end if
      else
        if (abs(a_work(j, j)) <= thresh * abs(a_work(1, 1))) then
          res%status = STATUS_NUMERIC_UNRESOLVED
          res%status_name = "NUMERIC_UNRESOLVED"
          res%error_code = "rank_screen_pivot_below_threshold"
          res%updates_used = updates
          return
        end if
      end if

      ! Apply reflector to trailing columns j+1..n
      if (tau(j) > 0.0_real64) then
        do k = j + 1, n
          dot_val = v(j) * a_work(j, k)
          do i = j + 1, m
            dot_val = dot_val + v(i) * a_work(i, k)
          end do
          a_work(j, k) = a_work(j, k) - tau(j) * v(j) * dot_val
          do i = j + 1, m
            a_work(i, k) = a_work(i, k) - tau(j) * v(i) * dot_val
          end do
          updates = updates + int(2 * (m - j + 1), int64)
        end do
      end if

      if (updates > budget%max_updates) then
        res%status = STATUS_UNKNOWN_BUDGET
        res%status_name = "UNKNOWN_BUDGET"
        res%error_code = "updates_budget_exceeded"
        res%updates_used = updates
        return
      end if
    end do

    ! Extract R(1:n, 1:n)
    res%r_upper = 0.0_real64
    do j = 1, n
      do i = 1, j
        res%r_upper(i, j) = a_work(i, j)
      end do
    end do

    ! Construct explicit orthogonal Q(1:m, 1:m) by backward accumulation
    res%q_full = 0.0_real64
    do i = 1, m
      res%q_full(i, i) = 1.0_real64
    end do

    do j = n, 1, -1
      if (tau(j) > 0.0_real64) then
        v = 0.0_real64
        v(j+1:m) = a_work(j+1:m, j)
        ! Recover v(j) from unit norm property: v(j)^2 + sum_{i>j} v(i)^2 = 1
        dot_val = 0.0_real64
        do i = j + 1, m
          dot_val = dot_val + v(i) * v(i)
        end do
        if (dot_val < 1.0_real64) then
          v(j) = sqrt(max(0.0_real64, 1.0_real64 - dot_val))
          ! Note: sign of v(j) does not affect v*v^T if we stored v(j+1:m) with the same sign as v(j).
          ! Wait: v(j) had sign of (x1 - beta). Since beta has opposite sign of x1 (or negative if x1>=0),
          ! x1 - beta is ALWAYS >= 0 when x1 >= 0, and ALWAYS < 0 when x1 < 0!
          ! Even simpler: since beta = -sign(x1)*||x||, x1 - beta has the same sign as -beta!
          if (res%r_upper(j, j) > 0.0_real64) then
            v(j) = -v(j)
          end if
        else
          v(j) = 0.0_real64
        end if

        do k = j, m
          dot_val = 0.0_real64
          do i = j, m
            dot_val = dot_val + v(i) * res%q_full(i, k)
          end do
          do i = j, m
            res%q_full(i, k) = res%q_full(i, k) - tau(j) * v(i) * dot_val
          end do
          updates = updates + int(2 * (m - j + 1), int64)
        end do
      end if
    end do

    if (updates > budget%max_updates) then
      res%status = STATUS_UNKNOWN_BUDGET
      res%status_name = "UNKNOWN_BUDGET"
      res%error_code = "updates_budget_exceeded"
      res%updates_used = updates
      return
    end if

    ! Compute Frobenius norm of original A
    a_frob = 0.0_real64
    do j = 1, n
      col_nrm = scaled_vec_norm2(a(1:m, j))
      a_frob = hypot(a_frob, col_nrm)
    end do

    ! Step 3: Solve for each RHS with up to 3 iterative refinement corrections
    do rhs_idx = 1, nrhs
      ! Compute Q^T * b
      qtb = 0.0_real64
      do i = 1, n
        dot_val = 0.0_real64
        do k = 1, m
          dot_val = dot_val + res%q_full(k, i) * b(k, rhs_idx)
        end do
        qtb(i) = dot_val
      end do
      updates = updates + int(m * n, int64)

      call solve_upper_triangular(res%r_upper, n, qtb, y_cur)
      updates = updates + int(n * n, int64)

      ! Compute initial normal-equation gradient g = A_scaled^T * (b - A_scaled * y_cur)
      r_vec(1:m) = b(1:m, rhs_idx) - matmul(a_scaled(1:m, 1:n), y_cur(1:n))
      g_vec(1:n) = matmul(transpose(a_scaled(1:m, 1:n)), r_vec(1:m))
      best_g_norm = scaled_vec_norm2(g_vec(1:n))
      updates = updates + int(2 * m * n, int64)

      do iter = 1, 3
        if (best_g_norm <= 0.0_real64) exit
        call solve_lower_triangular_transpose(res%r_upper, n, g_vec, z_vec)
        call solve_upper_triangular(res%r_upper, n, z_vec, dy_vec)
        y_cand(1:n) = y_cur(1:n) + dy_vec(1:n)

        r_vec(1:m) = b(1:m, rhs_idx) - matmul(a_scaled(1:m, 1:n), y_cand(1:n))
        g_vec(1:n) = matmul(transpose(a_scaled(1:m, 1:n)), r_vec(1:m))
        cand_g_norm = scaled_vec_norm2(g_vec(1:n))
        updates = updates + int(2 * m * n + n * n, int64)

        if (cand_g_norm < best_g_norm) then
          y_cur(1:n) = y_cand(1:n)
          best_g_norm = cand_g_norm
        else
          exit
        end if
      end do

      if (updates > budget%max_updates) then
        res%status = STATUS_UNKNOWN_BUDGET
        res%status_name = "UNKNOWN_BUDGET"
        res%error_code = "updates_budget_exceeded"
        res%updates_used = updates
        return
      end if

      ! Unscale and unpivot into x(:, rhs_idx)
      do j = 1, n
        p = res%pivots(j)
        res%x(p, rhs_idx) = y_cur(j) * res%col_scales(j)
        if (.not. ieee_is_finite(res%x(p, rhs_idx))) then
          res%status = STATUS_NUMERIC_UNRESOLVED
          res%status_name = "NUMERIC_UNRESOLVED"
          res%error_code = "nonfinite_solution_component"
          res%updates_used = updates
          return
        end if
      end do

      ! Independently accumulate unscaled residual r = b - A * x and A^T * r
      do i = 1, m
        dot_val = b(i, rhs_idx)
        do j = 1, n
          dot_val = dot_val - a(i, j) * res%x(j, rhs_idx)
        end do
        r_vec(i) = dot_val
      end do
      r_norm = scaled_vec_norm2(r_vec(1:m))
      res%residual_norm(rhs_idx) = r_norm

      do j = 1, n
        dot_val = 0.0_real64
        do i = 1, m
          dot_val = dot_val + a(i, j) * r_vec(i)
        end do
        atr_vec(j) = dot_val
      end do
      atr_norm = scaled_vec_norm2(atr_vec(1:n))
      x_norm = scaled_vec_norm2(res%x(1:n, rhs_idx))
      updates = updates + int(2 * m * n, int64)

      denom = a_frob * r_norm + (a_frob * a_frob) * x_norm
      if (denom > 0.0_real64) then
        eta = atr_norm / denom
      else
        eta = 0.0_real64
      end if
      res%normal_eq_rel_residual(rhs_idx) = eta

      if (.not. ieee_is_finite(eta) .or. eta > 5.0e-11_real64) then
        res%status = STATUS_NUMERIC_UNRESOLVED
        res%status_name = "NUMERIC_UNRESOLVED"
        res%error_code = "first_order_residual_exceeded"
        res%updates_used = updates
        return
      end if
    end do

    res%updates_used = updates
    res%status = STATUS_COMPUTED
    res%status_name = "COMPUTED"
    res%error_code = "ok"
  end subroutine solve_scaled_qr

end module sqe_solver
