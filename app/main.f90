program sqe_cli
  use, intrinsic :: iso_fortran_env, only: real64, output_unit, error_unit
  use scaled_qr_evidence
  implicit none

  integer :: argc, i, m, n, nrhs
  character(len=64) :: arg
  logical :: do_demo, do_assert, as_json, recon_ok
  real(real64) :: a(6, 3), b(6, 2), x_true(3, 2)
  real(real64) :: rel_recon_err, orth_err, sol_err
  type(sqe_budget_t) :: budget
  type(sqe_result_t) :: res

  argc = command_argument_count()
  do_demo = .false.
  do_assert = .false.
  as_json = .false.

  if (argc == 0) then
    do_demo = .true.
    do_assert = .true.
  else
    do i = 1, argc
      call get_command_argument(i, arg)
      if (trim(arg) == "--help" .or. trim(arg) == "-h") then
        write(output_unit, '(A)') "Scaled QR Evidence v0.1.0 - Column-Pivoted Householder QR"
        write(output_unit, '(A)') "Usage: scaled-qr-evidence demo [--assert] [--json]"
        stop 0
      else if (trim(arg) == "demo") then
        do_demo = .true.
      else if (trim(arg) == "--assert") then
        do_assert = .true.
      else if (trim(arg) == "--json") then
        as_json = .true.
      end if
    end do
  end if

  if (.not. do_demo) then
    write(error_unit, '(A)') "Error: unsupported command; use 'demo [--assert] [--json]'"
    stop 2
  end if

  ! Construct synthetic 6x3 badly-scaled overdetermined system with 2 RHS
  m = 6
  n = 3
  nrhs = 2

  a(:, 1) = [1.0e12_real64,  2.0e12_real64, -1.0e12_real64,  3.0e12_real64,  1.0e12_real64, -2.0e12_real64]
  a(:, 2) = [1.0e-9_real64, -1.0e-9_real64,  2.0e-9_real64,  1.0e-9_real64, -3.0e-9_real64,  2.0e-9_real64]
  a(:, 3) = [2.0_real64,     1.0_real64,     3.0_real64,    -2.0_real64,     4.0_real64,     1.0_real64]

  x_true(:, 1) = [2.5e-12_real64, -1.5e9_real64, 0.75_real64]
  x_true(:, 2) = [-1.0e-12_real64, 3.0e9_real64, -2.25_real64]

  b = matmul(a, x_true)
  ! Add orthogonal perturbation to RHS 2 to exercise non-zero residual least squares
  b(1, 2) = b(1, 2) + 1.0e-4_real64
  b(2, 2) = b(2, 2) - 1.0e-4_real64

  call solve_scaled_qr(a, b, m, n, nrhs, 30, budget, res)
  call verify_qr_factorization(a, res, rel_recon_err, orth_err, recon_ok)

  sol_err = scaled_vec_norm2(res%x(1:n, 1) - x_true(1:n, 1)) / scaled_vec_norm2(x_true(1:n, 1))

  if (do_assert) then
    if (res%status /= STATUS_COMPUTED .or. .not. recon_ok .or. sol_err > 1.0e-10_real64) then
      write(error_unit, '(A)') "Demo assertion failure"
      stop 1
    end if
  end if

  if (as_json) then
    write(output_unit, '(A,ES12.5,A,ES12.5,A,ES12.5,A)') &
      '{"status":"COMPUTED","rel_recon_err":', rel_recon_err, &
      ',"orth_err":', orth_err, &
      ',"normal_eq_rel_residual":', res%normal_eq_rel_residual(2), '}'
  else
    write(output_unit, '(A,ES12.5,A,ES12.5,A,ES12.5)') &
      "PASS status=COMPUTED rel_recon_err=", rel_recon_err, &
      " orth_err=", orth_err, &
      " normal_eq_rel_residual=", res%normal_eq_rel_residual(2)
  end if

end program sqe_cli
