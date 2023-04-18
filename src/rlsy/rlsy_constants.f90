module rlsy_constants
!! Sensible place to store choices of precision.
! I use this instead of iso_fortran_env to be consistent with ELSI.
use, intrinsic :: ISO_C_BINDING, only: c_float,c_double,c_int32_t,c_int64_t
use, intrinsic :: iso_fortran_env, only: REAL128
implicit none

!> precision definitions, consistent with ELSI
integer, parameter :: r4 = c_float
integer, parameter :: r8 = c_double
integer, parameter :: r16 = REAL128
integer, parameter :: i4 = c_int32_t
integer, parameter :: i8 = c_int64_t

!> imaginary i
complex(r8), parameter :: rl_imag = (0.0_r8,1.0_r8)
!> pi
real(r8),parameter :: rl_pi=3.141592653589793_r8
!> 2*pi
real(r8),parameter :: rl_twopi=6.283185307179586_r8

! Some default tolerances

!> Tolerance for realspace distances to be 0
real(r8), parameter :: rl_tol=1E-5_r8
!> Tolerance for realspace squared distances to be 0
real(r8), parameter :: rl_sqtol=rl_tol**2
!> Tolerance for reciprocal distances to be 0
real(r8), parameter :: rl_rectol=1E-6_r8
!> Tolerance for reciprocal squared distances
real(r8), parameter :: rl_sqrectol=rl_rectol**2
!> Tolerance for angles in degrees to be 0
real(r8), parameter :: rl_degreetol=1E-4_r8
!> Tolerance for angles in radians to be 0
real(r8), parameter :: rl_radiantol=rl_degreetol*180.0_r8/rl_pi
!> Tolerance for phonon frequencies to be 0.
real(r8), parameter :: rl_freqtol=rl_tol*1E-4_r8
!> Tolerance for phonon group velocities to be zero
!real(r8), parameter :: rl_phonongroupveltol=rl_tol*1E-5_r8
!> Tolerance for temperatures to be 0, in K
real(r8), parameter :: rl_temperaturetol=1E-3_r8
!> large number
real(r8), parameter :: rl_huge=huge(1.0_r8)
!> small number
real(r8), parameter :: rl_tiny=tiny(1.0_r8)
!> large integer
real(r8), parameter :: rl_hugeint=huge(1)

!> Variable that holds exit status, for catching exceptions
integer, save :: rl_status=0
!> IO unit for printing to stdout/something else.
integer, save :: rl_iou=6

! Exit codes for when things go really wrong.

!> Something has an unexpected dimension
integer, parameter :: rl_exitcode_baddim=1
!> BLAS or LAPACK fails
integer, parameter :: rl_exitcode_blaslapack=2
!> Unphysical value, i.e. negative temperature or something like that.
integer, parameter :: rl_exitcode_physical=3
!> Bad symmetry
integer, parameter :: rl_exitcode_symmetry=4
!> Something off with the arguments sent to the routine
integer, parameter :: rl_exitcode_param=5
!> IO error
integer, parameter :: rl_exitcode_io=6
!> MPI problem
integer, parameter :: rl_exitcode_mpi=7
!> Memory problem
integer, parameter :: rl_exitcode_memory=8
end module rlsy_constants
