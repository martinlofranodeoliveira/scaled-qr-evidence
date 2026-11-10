module sqe_verifier
  use, intrinsic :: iso_fortran_env, only: real64
  use sqe_types
  use sqe_solver, only: scaled_vec_norm2
  implicit none
  private

  public :: verify_qr_factorization

contains

  subroutine verify_qr_factorization(a, res, rel_recon_err, orth_err, passed)
    real(real64), intent(in) :: a(:,:)
    type(sqe_result_t), intent(in) :: res
    real(real64), intent(out) :: rel_recon_err
    real(real64), intent(out) :: orth_err
    logical, intent(out) :: passed

    integer :: m, n, i, j, k, p
    real(real64) :: a_sp(SQE_MAX_M, SQE_MAX_N)
    real(real64) :: qr(SQE_MAX_M, SQE_MAX_N)
    real(real64) :: qtq_minus_i(SQE_MAX_M, SQE_MAX_M)
    real(real64) :: diff_col(SQE_MAX_M)
    real(real64) :: a_sp_frob, recon_frob, s

    m = res%m
    n = res%n
    rel_recon_err = 1.0_real64
    orth_err = 1.0_real64
    passed = .false.

    if (res%status /= STATUS_COMPUTED) return

    ! Reconstruct A_scaled * P
    a_sp_frob = 0.0_real64
    do j = 1, n
      p = res%pivots(j)
      do i = 1, m
        a_sp(i, j) = a(i, p) * res%col_scales(j)
      end do
      a_sp_frob = hypot(a_sp_frob, scaled_vec_norm2(a_sp(1:m, j)))
    end do

    ! Compute Q(1:m, 1:n) * R(1:n, 1:n)
    recon_frob = 0.0_real64
    do j = 1, n
      do i = 1, m
        s = 0.0_real64
        do k = 1, j
          s = s + res%q_full(i, k) * res%r_upper(k, j)
        end do
        qr(i, j) = s
        diff_col(i) = a_sp(i, j) - qr(i, j)
      end do
      recon_frob = hypot(recon_frob, scaled_vec_norm2(diff_col(1:m)))
    end do

    if (a_sp_frob > 0.0_real64) then
      rel_recon_err = recon_frob / a_sp_frob
    else
      rel_recon_err = recon_frob
    end if

    ! Compute ||Q^T Q - I_m||_F
    orth_err = 0.0_real64
    do j = 1, m
      do i = 1, m
        s = 0.0_real64
        do k = 1, m
          s = s + res%q_full(k, i) * res%q_full(k, j)
        end do
        if (i == j) s = s - 1.0_real64
        qtq_minus_i(i, j) = s
      end do
      orth_err = hypot(orth_err, scaled_vec_norm2(qtq_minus_i(1:m, j)))
    end do

    passed = (rel_recon_err <= 5.0e-12_real64) .and. (orth_err <= 5.0e-12_real64)
  end subroutine verify_qr_factorization

end module sqe_verifier
