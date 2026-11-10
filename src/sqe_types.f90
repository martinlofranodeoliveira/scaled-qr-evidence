module sqe_types
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none
  private

  integer, parameter, public :: SQE_MAX_M = 128
  integer, parameter, public :: SQE_MAX_N = 32
  integer, parameter, public :: SQE_MAX_NRHS = 4
  integer, parameter, public :: SQE_MIN_K = 10
  integer, parameter, public :: SQE_MAX_K = 40
  integer(int64), parameter, public :: SQE_DEFAULT_MAX_UPDATES = 2000000_int64
  integer(int64), parameter, public :: SQE_DEFAULT_MAX_SCRATCH_BYTES = 1048576_int64
  integer(int64), parameter, public :: SQE_DEFAULT_MAX_OUTPUT_BYTES = 65536_int64

  integer, parameter, public :: STATUS_COMPUTED = 0
  integer, parameter, public :: STATUS_INVALID_INPUT = 1
  integer, parameter, public :: STATUS_NUMERIC_UNRESOLVED = 2
  integer, parameter, public :: STATUS_UNKNOWN_BUDGET = 3

  type, public :: sqe_budget_t
    integer(int64) :: max_updates = SQE_DEFAULT_MAX_UPDATES
    integer(int64) :: max_scratch_bytes = SQE_DEFAULT_MAX_SCRATCH_BYTES
    integer(int64) :: max_output_bytes = SQE_DEFAULT_MAX_OUTPUT_BYTES
  end type sqe_budget_t

  type, public :: sqe_result_t
    integer :: status = STATUS_INVALID_INPUT
    character(len=24) :: status_name = "INVALID_INPUT"
    character(len=48) :: error_code = "uninitialized"
    integer :: m = 0
    integer :: n = 0
    integer :: nrhs = 0
    integer :: threshold_k = 30
    integer(int64) :: updates_used = 0_int64
    integer(int64) :: scratch_bytes_used = 0_int64
    integer :: pivots(SQE_MAX_N) = 0
    real(real64) :: col_scales(SQE_MAX_N) = 0.0_real64
    real(real64) :: x(SQE_MAX_N, SQE_MAX_NRHS) = 0.0_real64
    real(real64) :: residual_norm(SQE_MAX_NRHS) = 0.0_real64
    real(real64) :: normal_eq_rel_residual(SQE_MAX_NRHS) = 0.0_real64
    real(real64) :: q_full(SQE_MAX_M, SQE_MAX_M) = 0.0_real64
    real(real64) :: r_upper(SQE_MAX_N, SQE_MAX_N) = 0.0_real64
  end type sqe_result_t

end module sqe_types
