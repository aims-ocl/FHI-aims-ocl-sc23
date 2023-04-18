!! Copyright (C) 2015 T. Markovich, M. Forsythe
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the GNU Lesser General Public
!! License as published by the Free Software Foundation; either
!! version 3.0 of the License, or (at your option) any later version.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library.

module mbdvdw_module
  use mbdvdw_interface_module
  use quadrature_grid_module, only: &
      casimir_omega, casimir_omega_weight, generate_grid
  use mbd_interface, only: print_log
  use localorb_io, only: default_unit, use_unit

  implicit none
  private
  save

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: VdW Radii at various levels of theory and their derivatives.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:),     allocatable :: R_vdw_free

  real(dp), dimension(:),     allocatable :: R_TS_VdW
  real(dp), dimension(:,:,:), allocatable :: dR_TS_VdWdR
  real(dp), dimension(:,:,:), allocatable :: dR_TS_VdWdh
  real(dp), dimension(:, :),  allocatable :: dR_TS_VdWdV

  real(dp), dimension(:),     allocatable :: R_MBD_VdW
  real(dp), dimension(:,:,:), allocatable :: dR_MBD_VdWdR
  real(dp), dimension(:,:,:), allocatable :: dR_MBD_VdWdh
  real(dp), dimension(:, :),  allocatable :: dR_MBD_VdWdV

  real(dp), dimension(:),     allocatable, public :: R_MBD_VdW_sl
  real(dp), dimension(:,:,:), allocatable, public :: dR_MBD_VdWdR_sl
  real(dp), dimension(:,:,:), allocatable, public :: dR_MBD_VdWdh_sl
  real(dp), dimension(:,:),   allocatable :: dR_MBD_VdWdV_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Screened oscillator frequency
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:),     allocatable :: C6_free

  real(dp), dimension(:),     allocatable, public :: omega_scs
  real(dp), dimension(:,:,:), allocatable :: domegadR
  real(dp), dimension(:,:,:), allocatable, public :: domegadh
  real(dp), dimension(:, :),  allocatable :: domegadV

  real(dp), dimension(:),     allocatable, public :: omega_scs_sl
  real(dp), dimension(:,:,:), allocatable, public :: domegadR_sl
  real(dp), dimension(:,:,:), allocatable, public :: domegadh_sl
  real(dp), dimension(:,:),   allocatable, public :: domegadV_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Polarizabilities at various levels of theory and their derivatives
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp), dimension(:),     allocatable :: alpha_free

  real(dp), dimension(:),     allocatable :: alpha_ts
  real(dp), dimension(:,:,:), allocatable :: dalpha_tsdR
  real(dp), dimension(:,:,:), allocatable :: dalpha_tsdh
  real(dp), dimension(:, :),  allocatable :: dalpha_tsdV

  real(dp), dimension(:),     allocatable, public :: alpha_0
  real(dp), dimension(:,:,:), allocatable, public :: dalpha_0dR
  real(dp), dimension(:,:,:), allocatable, public :: dalpha_0dh
  real(dp), dimension(:, :),  allocatable :: dalpha_0dV

  real(dp), dimension(:),     allocatable :: alpha_0_sl
  real(dp), dimension(:,:,:), allocatable :: dalpha_0dR_sl
  real(dp), dimension(:,:,:), allocatable :: dalpha_0dh_sl
  real(dp), dimension(:,:),   allocatable :: dalpha_0dV_sl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: VdW Radii at various levels of theory and their derivatives.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp), dimension(:),     allocatable :: sigma
  real(dp), dimension(:,:,:), allocatable :: dsigmadR
  real(dp), dimension(:,:,:), allocatable :: dsigmadh
  real(dp), dimension(:, :),  allocatable :: dsigmadV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Working variables
  ! molecular_polarizability: This is the 3x3 tensor that's been contracted
  !                                      to be the polarizability in each cartesian direction
  ! A_matrix: This is the 3Nx3N frequency dependent matrix that we'll
  !                populate by (A + T)^-1. In the nomenclature of the MBD forces paper
  !                we'll understand this paper as having a bar over top
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:,:),     allocatable :: A_matrix
  real(dp), dimension(:,:,:),   allocatable :: dA_matrixdR
  real(dp), dimension(:,:,:),   allocatable :: dA_matrixdh
  real(dp), dimension(:,:,:),   allocatable :: dA_matrixdV

  real(dp), dimension(:, :),    allocatable :: Cpq
  real(dp), dimension(:,:,:),   allocatable :: dCpqdR
  real(dp), dimension(:,:,:),   allocatable :: dCpqdh
  real(dp), dimension(:,:,:),   allocatable :: dCpqdV

  complex(dp), dimension(:,:,:),    allocatable, public :: Cpq_c
  complex(dp), dimension(:,:,:),   allocatable, public :: dCpqdR_c
  complex(dp), dimension(:,:,:),   allocatable, public :: dCpqdh_c
  complex(dp), dimension(:,:,:),   allocatable :: dCpqdV_c

  real(dp) :: interacting_energy
  real(dp) :: non_interacting_energy
  real(dp), public  :: EmbdvdW

  real(dp), dimension(:,:), allocatable :: int_FmbdVdW
  real(dp), dimension(:,:), allocatable :: nonint_FmbdVdW
  real(dp), dimension(:,:), allocatable, public  :: FmbdVdW

  real(dp), dimension(:,:), allocatable :: int_HmbdVdW
  real(dp), dimension(:,:), allocatable :: nonint_HmbdVdW
  real(dp), dimension(:,:), allocatable, public  :: HmbdVdW

  real(dp), dimension(:),   allocatable :: int_Uprefactor
  real(dp), dimension(:),   allocatable :: nonint_Uprefactor
  real(dp), dimension(:),   allocatable :: Uprefactor
  real(dp), dimension(:),   allocatable, public  :: UmbdvdW

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Variable Declarations: Simulation parameters
  ! beta: Damping factor associated with the functional used. Only defined for
  !         PBE, beta=0.83
  !         PBE0, beta=0.85
  !         HSE, beta=0.85
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp)               :: beta
  real(dp)               :: dip_cutoff
  real(dp)               :: supercell_cutoff
  logical, public                 :: do_forces
  logical, public                 :: mbd_conv_elec
  logical, public                 :: mbd_first_step
  real(dp)                        :: omega_to_print

  REAL(DP) :: h_(3,3)                                 !HK: PW and CP use difference var for cell (at) and (h), so this is wrapper between the two
  REAL(DP), public :: ainv_(3,3)                              !HK: PW and CP use difference var for cell (at) and (h), so this is wrapper between the two


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Super Lattice Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, public                                 :: nat_sl
  real(dp), dimension(3, 3)              :: h_sl
  real(dp), dimension(3, 3)              :: ainv_sl
  real(dp), dimension(:, :),  allocatable :: tau_sl, tau_sl_s
  real(dp), dimension(:,:),   allocatable          :: tau
  real(dp), dimension(:,:),   allocatable :: tau_s
  integer                                          :: sl_i, sl_j, sl_k
  real(dp), dimension(3), public                           :: sl_mult
  integer, dimension(:), allocatable :: orig_idx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Kpoint Grid Variables and Ewald variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp), dimension(:, :), allocatable :: k_grid
  real(dp), public :: dadh(3, 3), ruc(3, 3), Rc, Gc, ewald_cutoff, nu
  real(dp) :: dRcdh(3, 3), dGcdh(3, 3), dnudh(3,3)
  integer, public :: nk
  complex(dp), public :: II = cmplx(0, 1.0_DP, kind=DP)
  logical :: do_ewald, do_recip, do_cmplx, MBD_PER_ERROR

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! assorted constants
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp) :: SQRTPI = 1.772453850905516027

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parallel variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type pairs
    integer p
    integer q
    integer CPU
  end type pairs

  integer                           :: num_pairs
  type(pairs), dimension(:), allocatable     :: unique_pairs
  type(pairs), dimension(:), allocatable     :: pairs_scs
  integer, dimension(:), allocatable         :: f_cpu_id, h_cpu_id, v_cpu_id, k_cpu_id
  integer                                    :: max_proc_forces, max_proc_sc, max_proc_h, max_proc_k
  integer                                    :: max_proc_pairs

  integer :: npts

  integer, dimension(:), allocatable :: n_pairs, n_comps_f, n_comps_v, n_comps_h, n_comps_k


  character(len=80), parameter :: neg_eigvals_msg(9) = [character(len=80) :: &
    "Number of negative eigenvalues larger than threshold, aborting.", &
    "The MBD energy is obtained by diagonalizing an MBD Hamiltonian. When", &
    "a system is overpolarized, the digonalization produces negative", &
    "eigenvalues, which are unphysical. This happens often with ionic and", &
    "metallic systems, for which the MBD@rsSCS parametrization doesn't work well.", &
    "When the number of negative eigenvalues is smaller than a certain threshold", &
    "(default 3), they are set to zero. When the threshold is exceeded, the", &
    "calculation aborts. For metallic systems, a potential solution is to use", &
    "the MBD@surf method by manually specifying the vdW parameters." &
  ]

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Method Prototypes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  public :: mbdvdw_initialize, mbdvdw_calculate, mbdvdw_cleanup, mbdvdw_finalize, inv3x3_mat
  ! Public methods only for debugging
  public :: compute_recucell, mbdvdw_b, mbdvdw_dbdh, mbdvdw_dbdr, mbdvdw_c, mbdvdw_dcdh, mbdvdw_dcdr, mbdvdw_cellvol, &
            mbdvdw_compute_drdh, mbdvdw_compute_dRdh_only, mbdvdw_compute_dRdR, mbdvdw_compute_dRdR_only, &
            mbdvdw_compute_tdip, mbdvdw_compute_TLR, mbdvdw_compute_TLR_complex, mbdvdw_TLR_ewald_realspace, &
            mbdvdw_TGG_complex, mbdvdw_TGG, mbdvdw_tdip_ewald, mbdvdw_ewaldsum, mbdvdw_ewald_params, mbdvdw_dewald_params_dh, &
            mbdvdw_dcellvoldh, mbdvdw_construct_sl, mbdvdw_compute_dRdR_ewald, mbdvdw_drec_vecdh, mbdvdw_drec_vec_norm_dh, &
            mbdvdw_compute_TSR, mbdvdw_effqts

  CONTAINS

  subroutine inv3x3_mat(A,T)
    IMPLICIT NONE
    real(DP), intent(in) :: A(3,3)
    real(DP) :: D
    real(DP), intent(out) :: T(3,3)
        D = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)- &
        A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
        A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
        T(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/D
        T(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/D
        T(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/D
        T(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/D
        T(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/D
        T(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/D
        T(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/D
        T(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/D
        T(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/D
      return
    end subroutine inv3x3_mat

  function determinant3x3_mat(A) result(det)
    implicit none
    real(dp), intent(in) :: A(3, 3)
    real(dp) :: det

    det = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)- &
          A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
          A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
    end function

  subroutine mbdvdw_zeroout()
    R_TS_VdW    = 0.0_DP
    R_MBD_VdW    = 0.0_DP
    alpha_ts = 0.0_DP
    alpha_0 = 0.0_DP
    omega_scs = 0.0_DP
    sigma = 0.0_DP

    if(vdw_self_consistent) then
      dR_TS_VdWdV = 0.0_DP
      dR_MBD_VdWdv = 0.0_DP
      domegadV = 0.0_DP
      dalpha_tsdV = 0.0_DP
      dalpha_0dV = 0.0_DP
      dsigmadV = 0.0_DP
      end if

    if(do_forces) then
      dR_TS_VdWdR = 0.0_DP
      dR_TS_VdWdh = 0.0_DP
      dR_MBD_VdWdR = 0.0_DP
      dR_MBD_VdWdh = 0.0_DP
      domegadR = 0.0_DP
      domegadh = 0.0_DP
      dalpha_tsdR = 0.0_DP
      dalpha_tsdh = 0.0_DP
      dalpha_0dR = 0.0_DP
      dalpha_0dh = 0.0_DP
      dsigmadR = 0.0_DP
      dsigmadh = 0.0_DP
    end if

    nonint_Uprefactor = 0.0_DP
    int_Uprefactor = 0.0_DP
    Uprefactor = 0.0_DP

    nonint_FmbdVdW = 0.0_DP
    int_FmbdVdW = 0.0_DP
    FmbdVdW = 0.0_DP

    nonint_HmbdVdW = 0.0_DP
    int_HmbdVdW = 0.0_DP
    HmbdVdW = 0.0_DP

    tau_s = 0.0_DP
    n_pairs = 0

    end subroutine mbdvdw_zeroout

  subroutine mbdvdw_get_is(i_f, s, i)
    implicit none
    integer, intent(out) :: s, i
    integer, intent(in) :: i_f
    i = modulo(i_f-1, 3) + 1
    s = ceiling(i_f/dble(3))
    end subroutine mbdvdw_get_is

  subroutine mbdvdw_get_isk(i_f, s, i, k)
    implicit none
    integer, intent(out) :: s, i, k
    integer, intent(in) :: i_f
    k = modulo(i_f-1, nk) + 1
    s = ceiling(dble(i_f)/dble(3*nk))
    i = ceiling(dble(i_f)/dble(nk))
    i = modulo(i - 1, 3) + 1
    end subroutine mbdvdw_get_isk

  subroutine mbdvdw_get_ks(i_f, s, k)
    implicit none
    integer, intent(out) :: s, k
    integer, intent(in) :: i_f
    k = modulo(i_f-1, nk) + 1
    s = ceiling(dble(i_f)/nk)
    end subroutine mbdvdw_get_ks

  subroutine mult_TNN(N, r_mat, c_mat, result)
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in) :: r_mat(N, N), c_mat(N, N)
    real(dp), intent(out) :: result(N, N)
    real(dp) :: temp1(N, N)
    temp1 = 0.0_DP
    call dgemm('N', 'N', N, N, N, 1.0_DP, c_mat, N, r_mat, N, 0.0_DP, temp1, N)
    call dgemm('T', 'N', N, N, N, 1.0_DP, r_mat, N, temp1, N, 0.0_DP, result, N)
    end subroutine mult_TNN

  subroutine mult_CNN(N, r_mat, c_mat, result)
    implicit none
    integer, intent(in) :: N
    complex(dp), intent(in) :: r_mat(N, N), c_mat(N, N)
    complex(dp), intent(out) :: result(N, N)
    complex(dp) :: temp1(N, N)
    temp1 = 0.0_DP
    call zgemm('N', 'N', N, N, N, (1.0_DP, 0.0_DP), c_mat, N, r_mat, N, (0.0_DP, 0.0_DP), temp1, N)
    call zgemm('C', 'N', N, N, N, (1.0_DP, 0.0_DP), r_mat, N, temp1, N, (0.0_DP, 0.0_DP), result, N)
    end subroutine mult_CNN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These three functions were taken from Jan Herman's code

  function make_g_grid(n1, n2, n3) result(g_grid)
      integer, intent(in) :: n1, n2, n3
      real(dp) :: g_grid(n1*n2*n3, 3)

      integer :: g_kpt(3), i_kpt, kpt_range(3)
      real(dp) :: g_kpt_shifted(3)

      g_kpt = (/ 0, 0, -1 /)
      kpt_range = (/ n1, n2, n3 /)
      do i_kpt = 1, n1*n2*n3
          call shift_cell (g_kpt, (/ 0, 0, 0 /), kpt_range-1)
          g_kpt_shifted = dble(g_kpt)+mbd_vdw_kgrid_shift
          where (2*g_kpt_shifted > kpt_range)
              g_kpt_shifted = g_kpt_shifted-dble(kpt_range)
          end where
          g_grid(i_kpt, :) = g_kpt_shifted/kpt_range
      end do
    end function make_g_grid


  function make_k_grid(n1, n2, n3, uc) result(k_grid)
      real(dp), intent(in) :: uc(3,3)
      integer, intent(in) :: n1, n2, n3
      real(dp) :: g_grid(n1*n2*n3, 3)
      real(dp) :: k_grid(n1*n2*n3, 3)
      integer :: i_kpt
      real(dp) :: ruc(3, 3)

      g_grid = make_g_grid(n1, n2, n3)
      ruc = compute_recucell(uc)
      do i_kpt = 1, n1*n2*n3
          k_grid(i_kpt, :) = matmul(g_grid(i_kpt, :), ruc)
      end do
    end function make_k_grid

  subroutine shift_cell(ijk, first_cell, last_cell)
      integer, intent(inout) :: ijk(3)
      integer, intent(in) :: first_cell(3), last_cell(3)

      integer :: i_dim, i

      do i_dim = 3, 1, -1
          i = ijk(i_dim)+1
          if (i <= last_cell(i_dim)) then
              ijk(i_dim) = i
              return
          else
              ijk(i_dim) = first_cell(i_dim)
          end if
      end do
    end subroutine

! The above three functions were taken from Jan Herman's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mbdvdw_compute_dRdR(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dRmatdR, dTdipdR)
    implicit none
    integer, intent(in)   :: s, i, p, q
    real(dp), intent(in)  :: Rpq(3), Rpq_norm
    real(dp), intent(out) :: dRpqdR(3), dRpq_normdr, dRmatdR(3,3), dTdipdR(3,3)
    integer :: i_idx, j_idx

    call mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
    do i_idx=1,3,1
      do j_idx=1,3,1
        dRmatdR(i_idx, j_idx) = (Rpq(i_idx)*dRpqdR(j_idx) + dRpqdR(i_idx)*Rpq(j_idx))/Rpq_norm**5.0_DP&
                        - 5.0_DP*Rpq(i_idx)*Rpq(j_idx)*dRpq_normdR/Rpq_norm**6.0_DP
        dTdipdR(i_idx, j_idx) = -3.0_DP*dRmatdR(i_idx, j_idx)
        end do
      dTdipdR(i_idx, i_idx) = dTdipdR(i_idx, i_idx) - 3.0_DP*dRpq_normdR/Rpq_norm**4.0_DP
      end do
    end subroutine mbdvdw_compute_dRdR

  subroutine mbdvdw_compute_dRdR_ewald(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dTdipdR)
    implicit none
    integer, intent(in)   :: s, i, p, q
    real(dp), intent(in)  :: Rpq(3), Rpq_norm
    real(dp), intent(out) :: dRpqdR(3), dRpq_normdr, dTdipdR(3,3)
    integer :: i_idx, j_idx
    real(dp) :: b, c, dbdr, dcdr, top(3,3), dtop(3,3), r2, r3, r4, r5, r6

    call mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
    call mbdvdw_b(Rpq_norm, b)
    call mbdvdw_dbdr(Rpq_norm, dRpq_normdR, dbdr)
    call mbdvdw_c(Rpq_norm, c)
    call mbdvdw_dcdr(Rpq_norm, dRpq_normdR, dcdr)

    r2 = rpq_norm*rpq_norm
    r3 = r2*rpq_norm
    r4 = r2*r2
    r5 = r4*rpq_norm
    r6 = r4*r2

    do i_idx=1,3,1
      do j_idx=1,3,1
        top(i_idx, j_idx) = -(Rpq(i_idx)*Rpq(j_idx)*C)
        dtop(i_idx, j_idx) = -(dRpqdr(i_idx)*Rpq(j_idx)*C + Rpq(i_idx)*dRpqdr(j_idx)*C + Rpq(i_idx)*Rpq(j_idx)*dcdr)
        end do
        top(i_idx, i_idx) = top(i_idx, i_idx) + Rpq_norm*Rpq_norm*B
        dtop(i_idx, i_idx) = dtop(i_idx, i_idx) + (2.0_dp*Rpq_norm*B*dRpq_normdr + Rpq_norm*Rpq_norm*dBdr)
      end do
    dTdipdr = dtop/R5 - 5.0_dp*top*dRpq_normdr/R6
    end subroutine mbdvdw_compute_dRdR_ewald

  subroutine mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
    implicit none
    integer, intent(in)   :: s, i, p, q
    real(dp), intent(in)  :: Rpq(3), Rpq_norm
    real(dp), intent(out) :: dRpqdR(3), dRpq_normdr

    dRpqdR = 0.0_DP
    dRpq_normdR = 0.0_DP
    if(orig_idx(q).eq.s) then
      dRpqdR(i) = 1.0_DP
      dRpq_normdR = Rpq(i)/(Rpq_norm)
      else if(orig_idx(p).eq.s) then
      dRpqdR(i) = -1.0_DP
      dRpq_normdR = -Rpq(i)/(Rpq_norm)
      end if
    if(orig_idx(p).eq.orig_idx(q)) then
      dRpqdR = 0.0_dp
      dRpq_normdR = 0.0_dp
      end if
    end subroutine mbdvdw_compute_dRdR_only

  subroutine mbdvdw_compute_dRdh(s, i, Rpq, Rpq_norm, Spq, dRpqdh, dRpq_normdh, dRmatdh, dTdipdh)
    implicit none
    integer, intent(in)   :: s, i
    real(dp), intent(in)  :: Rpq(3), Rpq_norm, Spq(3)
    real(dp), intent(out) :: dRpqdh(3), dRpq_normdh, dRmatdh(3,3), dTdipdh(3,3)
    integer :: i_idx, j_idx

    call mbdvdw_compute_dRdh_only(s, i, Rpq, Rpq_norm, Spq, dRpqdh, dRpq_normdh)

    do i_idx=1,3,1
      do j_idx=1,3,1
        dRmatdH(i_idx, j_idx) = (Rpq(i_idx)*dRpqdH(j_idx) + dRpqdH(i_idx)*Rpq(j_idx))/Rpq_norm**5.0_DP&
                        - 5.0_DP*Rpq(i_idx)*Rpq(j_idx)*dRpq_normdH/Rpq_norm**6.0_DP
        dTdipdH(i_idx, j_idx) = -3.0_DP*dRmatdH(i_idx, j_idx)
        end do
      dTdipdH(i_idx, i_idx) = dTdipdH(i_idx, i_idx) - 3.0_DP*dRpq_normdH/Rpq_norm**4.0_DP
      end do
    return
    end subroutine mbdvdw_compute_dRdh

  subroutine mbdvdw_compute_dRdh_ewald(s, i, Rpq, Rpq_norm, Spq, dadh, dRpqdh, dRpq_normdh, dTdipdh)
    implicit none
    integer, intent(in)   :: s, i
    real(dp), intent(in)  :: Rpq(3), Rpq_norm, Spq(3), dadh
    real(dp), intent(out) :: dRpqdh(3), dRpq_normdh, dTdipdh(3,3)
    integer :: i_idx, j_idx
    real(dp) :: top(3,3), dtop(3,3), r2, r3, r4, r5, r6, b, c, dbdh, dcdh

    call mbdvdw_compute_dRdh_only(i, s, Rpq, Rpq_norm, Spq, dRpqdh, dRpq_normdh)
    call mbdvdw_b(Rpq_norm, b)
    call mbdvdw_dbdh(Rpq_norm, dRpq_normdh, dadh, dbdh)
    call mbdvdw_c(Rpq_norm, c)
    call mbdvdw_dcdh(Rpq_norm, dRpq_normdh, dadh, dcdh)

    r2 = rpq_norm*rpq_norm
    r3 = r2*rpq_norm
    r4 = r2*r2
    r5 = r4*rpq_norm
    r6 = r4*r2

    do i_idx=1,3,1
      do j_idx=1,3,1
        top(i_idx, j_idx) = -(Rpq(i_idx)*Rpq(j_idx)*C)
        dtop(i_idx, j_idx) = -(dRpqdh(i_idx)*Rpq(j_idx)*C + Rpq(i_idx)*dRpqdh(j_idx)*C + Rpq(i_idx)*Rpq(j_idx)*dcdh)
        end do
        top(i_idx, i_idx) = top(i_idx, i_idx) + Rpq_norm*Rpq_norm*B
        dtop(i_idx, i_idx) = dtop(i_idx, i_idx) + (2.0_dp*Rpq_norm*B*dRpq_normdh + Rpq_norm*Rpq_norm*dBdh)
      end do
    dTdipdh = dtop/R5 - 5.0_dp*top*dRpq_normdh/R6
    end subroutine mbdvdw_compute_dRdh_ewald

  subroutine mbdvdw_compute_dRdh_only(s, i, Rpq, Rpq_norm, Spq, dRpqdh, dRpq_normdh)
    implicit none
    integer, intent(in)   :: s, i
    real(dp), intent(in)  :: Rpq(3), Rpq_norm, Spq(3)
    real(dp), intent(out) :: dRpqdh(3), dRpq_normdh

    dRpqdh = 0.0_DP
    dRpq_normdh = 0.0_DP
    dRpq_normdh = -Rpq(s)*Spq(i)*dble(sl_mult(i))/(rpq_norm)
    dRpqdh(s) = -Spq(i)*dble(sl_mult(i))

    return
    end subroutine mbdvdw_compute_dRdh_only

  subroutine mbdvdw_compute_tdip(Rpq, R3, R5, tdip)
    real(dp), intent(in) :: Rpq(3), R3, R5
    real(dp), intent(out) :: tdip(3,3)
    integer i, j

    do i=1,3,1
      do j=1,3,1
        Tdip(i, j) = -3.0_dp*Rpq(i)*Rpq(j)/(R5)
        end do
      Tdip(i, i) = Tdip(i, i) + 1.0_DP/(R3)
      end do
    end subroutine mbdvdw_compute_tdip

  subroutine mbdvdw_b(R, b)
    real(dp), intent(in) :: R
    real(dp), intent(out) :: b
    real(dp) :: R2, t1, t2, t3, a
    a = ewald_cutoff
    R2 = R*R
    t1 = erfc(a*r)
    t2 = 2.0_dp*a*R/SQRTPI
    t3 = exp(-(a*R)**2.0_dp)
    b = t1 + t2*t3
    end subroutine mbdvdw_b

  subroutine mbdvdw_dbdR(R, dRdR, dbdR)
    real(dp), intent(in) :: R, dRdR
    real(dp), intent(out) :: dbdR
    real(dp) :: R2, t1, t2, t3, a
    real(dp) :: d1, d2, daRadR, pre
    R2 = R*R
    a = ewald_cutoff

    t1 = erfc(-a*R)
    t2 = (2.0_dp*a*R/SQRTPI)
    t3 = exp(-a*a*R2)

    daRadR = a*dRdR
    pre = (2.0_dp/SQRTPI)

    d1 = -pre*daRadR*t3
    d2 = pre*(daRadR*t3 - 2.0_dp*(a*R)**2.0_dp*daRadR*t3)

    dbdR = d1 + d2
    end subroutine mbdvdw_dbdR

  subroutine mbdvdw_dbdh(R, dRdh, dadh, dbdh)
    real(dp), intent(in) :: R, dRdh
    real(dp), intent(out) :: dbdh
    real(dp) :: R2, t1, t2, t3, dadh, a
    real(dp) :: d1, d2, daRadR, pre
    R2 = R*R
    a = ewald_cutoff

    t1 = erfc(-a*R)
    t2 = (2.0_dp*a*R/SQRTPI)
    t3 = exp(-a*a*R2)

    daRadR = dadh*R + a*dRdh
    pre = (2.0_dp/SQRTPI)

    d1 = -pre*daRadR*t3
    d2 = pre*(daRadR*t3 - 2.0_dp*(a*R)**2.0_dp*daRadR*t3)

    dbdh = d1 + d2
    end subroutine mbdvdw_dbdh

  subroutine mbdvdw_c(R, c)
    real(dp), intent(in) :: R
    real(dp), intent(out) :: c
    real(dp) :: t1, t2, t3, t4, a
    a = ewald_cutoff
    t1 = 3.0_dp*erfc(a*R)
    t2 = 2.0_dp*a*R/SQRTPI
    t3 = 3.0_dp + 2.0_dp*(a*R)**2.0_dp
    t4 = exp(-a*a*R*R)
    c = (t1+t2*t3*t4)
    end subroutine mbdvdw_c

  subroutine mbdvdw_dcdR(R, dRdR, dcdR)
    real(dp), intent(in) :: R, dRdR
    real(dp), intent(out) :: dcdR
    real(dp) :: t1, t2, t3, t4, a
    real(dp) :: dt1, dt2, dt3, dt4, pre, daRadR

    a = ewald_cutoff
    t1 = 3.0_dp*erfc(a*R)
    t2 = 2.0_dp*a*R/SQRTPI
    t3 = 3.0_dp + 2.0_dp*(a*R)**2.0_dp
    t4 = exp(-a*a*R*R)

    pre = 2.0_dp/SQRTPI
    daRadR = a*dRdR

    dt1 = -3.0_dp*pre*daRadR*t4
    dt2 = pre*daRadR
    dt3 = 4.0_dp*a*R*daRadR
    dt4 = (-2.0_dp)*a*R*daRadR*t4

    dcdR = dt1 + dt2*t3*t4+ t2*dt3*t4 + t2*t3*dt4
    end subroutine mbdvdw_dcdR

  subroutine mbdvdw_dcdh(R, dRdh, dadh, dcdh)
    real(dp), intent(in) :: R, dRdh, dadh
    real(dp), intent(out) :: dcdh
    real(dp) :: t1, t2, t3, t4, a
    real(dp) :: dt1, dt2, dt3, dt4, pre, daRadR

    a = ewald_cutoff
    t1 = 3.0_dp*erfc(a*R)
    t2 = 2.0_dp*a*R/SQRTPI
    t3 = 3.0_dp + 2.0_dp*(a*R)**2.0_dp
    t4 = exp(-a*a*R*R)

    pre = 2.0_dp/SQRTPI
    daRadR = dadh*R + a*dRdh

    dt1 = -3.0_dp*pre*daRadR*t4
    dt2 = pre*daRadR
    dt3 = 4.0_dp*a*R*daRadR
    dt4 = (-2.0_dp)*a*R*daRadR*t4

    dcdh = dt1 + dt2*t3*t4+ t2*dt3*t4 + t2*t3*dt4
    end subroutine mbdvdw_dcdh

  subroutine mbdvdw_kmat(k, k_norm, kmat)
    implicit none
    real(dp), intent(in) :: k(3), k_norm
    real(dp), intent(out) :: kmat(3, 3)
    integer :: i, j

    do i = 1, 3, 1
      do j = 1, 3, 1
        kmat(i, j) = k(i) * k(j)
        end do
      end do
    kmat = kmat/(k_norm*k_norm)
    return
    end subroutine mbdvdw_kmat

  subroutine mbdvdw_dkmatdh(k, dk, k_norm, dk_norm, dkmat)
    implicit none
    real(dp), intent(in) :: k(3), dk(3), k_norm, dk_norm
    real(dp), intent(out) :: dkmat(3, 3)
    real(dp) :: k2, k3
    integer :: i, j
    k2 = k_norm*k_norm
    k3 = k2*k_norm
    do i = 1, 3, 1
      do j = 1, 3, 1
        dkmat(i, j) = (dk(i)*k(j) + k(i)*dk(j))/k2
        dkmat(i, j) = dkmat(i, j) - 2.0_dp*k(i)*k(j)*(dk_norm)/k3
        end do
      end do
    return
    end subroutine mbdvdw_dkmatdh

  function compute_recucell(h_in) result(rec_ucell)
    implicit none
    real(dp), intent(in) :: h_in(3, 3)
    real(dp) :: rec_ucell(3, 3), temp(3, 3), ht(3,3)

    ht = transpose(h_in)
    call inv3x3_mat(ht, temp)
    rec_ucell = 2.0_DP*pi*transpose(temp)
    end function

  subroutine mbdvdw_drec_vecdh(rec_vec, a_in, gamma, sigma, drec_vec)
    implicit none
    real(dp), intent(in) :: rec_vec(3), a_in(3, 3)
    integer, intent(in) :: gamma, sigma
    real(dp), intent(out) :: drec_vec(3)
    integer :: i
    do i = 1, 3, 1
      drec_vec(i) = 1.0_DP*rec_vec(gamma)*a_in(sigma, i)/(2.0_DP*pi)*dble(sl_mult(sigma))
    end do
    end subroutine mbdvdw_drec_vecdh

  subroutine mbdvdw_drec_vec_norm_dh(rec_vec, rec_vec_norm, a_in, gamma, sigma, drec_vec_norm)
    implicit none
    real(dp), intent(in) :: rec_vec(3), rec_vec_norm, a_in(3, 3)
    integer, intent(in) :: gamma, sigma
    real(dp), intent(out) :: drec_vec_norm
    drec_vec_norm = (rec_vec(1)*a_in(sigma, 1) + rec_vec(2)*a_in(sigma, 2) + rec_vec(3)*a_in(sigma, 3))
    drec_vec_norm = rec_vec(gamma)*drec_vec_norm/(rec_vec_norm*2.0_dp*pi)*sl_mult(sigma)
    end subroutine mbdvdw_drec_vec_norm_dh

  subroutine mbdvdw_cellvol(h_in, nu)
    implicit none
    real(dp), intent(in) :: h_in(3,3)
    real(dp), intent(out) :: nu
    nu = abs(determinant3x3_mat(h_in))
    end subroutine mbdvdw_cellvol

  subroutine mbdvdw_dcellvoldh(ainv_in, gamma, sigma, nu, dnu)
    implicit none
    real(dp), intent(in) :: ainv_in(3,3)
    integer, intent(in) :: gamma, sigma
    real(dp), intent(in) :: nu
    real(dp), intent(out) :: dnu
    dnu = -nu*ainv_in(sigma, gamma)
    end subroutine mbdvdw_dcellvoldh

  subroutine mbdvdw_ewald_params(nu)
    implicit none
    real(dp), intent(in) :: nu
    real(dp) :: nu13

    nu13 = nu**(1.0_dp/3.0_dp)
    ewald_cutoff = 2.5_DP/nu13
    Rc = 2.4_dp*nu13
    Gc = 25.0_dp/nu13
    end subroutine mbdvdw_ewald_params

  subroutine mbdvdw_dewald_params_dh(ainv_in, nu)
    implicit none
    real(dp), intent(in) :: nu
    real(dp), intent(in) :: ainv_in(3,3)
    real(dp) :: nu43, nu23, dnu
    integer :: i, s

    nu43 = nu**(4.0_dp/3.0_dp)
    nu23 = nu**(2.0_dp/3.0_dp)

    do s = 1, 3, 1
      do i = 1, 3, 1
        call mbdvdw_dcellvoldh(ainv_in, s, i, nu, dnu)
        dadh(s, i) = -2.5_dp/(3.0_dp*nu43)*dnu
        dRcdh(s, i) = 2.4_dp/(3.0_dp*nu23)*dnu
        dGcdh(s, i) = -25.0_dp/(3.0_dp*nu43)*dnu
        end do
      end do
    end subroutine mbdvdw_dewald_params_dh

  function mbdvdw_circumscribe(uc, radius) result(sc)
    ! Taken from Jan Herman's code
    real(dp), intent(in) :: uc(3, 3), radius
    integer :: sc(3)
    real(dp) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = compute_recucell(uc)
    forall (i = 1:3) layer_sep(i) = sum(uc(i, :)*ruc(i, :)/dsqrt(sum(ruc(i, :)**2)))
    sc = ceiling(radius/layer_sep+0.5_dp)
    if(mbd_vdw_vacuum(1)) sc(1) = 0
    if(mbd_vdw_vacuum(2)) sc(2) = 0
    if(mbd_vdw_vacuum(3)) sc(3) = 0
    end function

  subroutine mbdvdw_cleanup_postscs()

    if(allocated(A_matrix))        deallocate(A_matrix)
    if(allocated(dA_matrixdR))     deallocate(dA_matrixdR)
    if(allocated(dA_matrixdh))     deallocate(dA_matrixdh)
    if(allocated(dA_matrixdV))     deallocate(dA_matrixdV)

    end subroutine mbdvdw_cleanup_postscs

  subroutine mbdvdw_cleanup()

    if(allocated(Cpq))              deallocate(Cpq)
    if(allocated(dCpqdR))           deallocate(dCpqdR)
    if(allocated(dCpqdh))           deallocate(dCpqdh)
    if(allocated(dCpqdV))           deallocate(dCpqdV)

    if(allocated(Cpq_c))              deallocate(Cpq_c)
    if(allocated(dCpqdR_c))           deallocate(dCpqdR_c)
    if(allocated(dCpqdh_c))           deallocate(dCpqdh_c)
    if(allocated(dCpqdV_c))           deallocate(dCpqdV_c)

    if(allocated(n_comps_k))        deallocate(n_comps_k)
    if(allocated(k_cpu_id))         deallocate(k_cpu_id)
    if(allocated(k_grid))           deallocate(k_grid)
    end subroutine mbdvdw_cleanup

  subroutine mbdvdw_finalize()
    call mbdvdw_cleanup_postscs()
    call mbdvdw_cleanup()

    if(allocated(tau))              deallocate(tau)
    if(allocated(tau_s))            deallocate(tau_s)
    if(allocated(tau_sl))           deallocate(tau_sl)
    if(allocated(tau_sl_s))         deallocate(tau_sl_s)
    if(allocated(n_pairs))          deallocate(n_pairs)

    if(allocated(R_TS_VdW))        deallocate(R_TS_VdW)
    if(allocated(dR_TS_VdWdR))     deallocate(dR_TS_VdWdR)
    if(allocated(dR_TS_VdWdh))     deallocate(dR_TS_VdWdh)
    if(allocated(dR_TS_VdWdV))     deallocate(dR_TS_VdWdV)

    if(allocated(alpha_ts))        deallocate(alpha_ts)
    if(allocated(dalpha_tsdR))     deallocate(dalpha_tsdR)
    if(allocated(dalpha_tsdh))     deallocate(dalpha_tsdh)
    if(allocated(dalpha_tsdV))     deallocate(dalpha_tsdV)

    if(allocated(R_MBD_VdW))       deallocate(R_MBD_VdW)
    if(allocated(dR_MBD_VdWdR))    deallocate(dR_MBD_VdWdR)
    if(allocated(dR_MBD_VdWdh))    deallocate(dR_MBD_VdWdh)
    if(allocated(dR_MBD_VdWdV))    deallocate(dR_MBD_VdWdV)

    if(allocated(domegadR))        deallocate(domegadR)
    if(allocated(domegadh))        deallocate(domegadh)
    if(allocated(domegadV))        deallocate(domegadV)

    if(allocated(alpha_0))         deallocate(alpha_0)
    if(allocated(dalpha_0dR))      deallocate(dalpha_0dR)
    if(allocated(dalpha_0dh))      deallocate(dalpha_0dh)
    if(allocated(dalpha_0dV))      deallocate(dalpha_0dV)

    if(allocated(omega_scs))       deallocate(omega_scs)
    if(allocated(sigma))           deallocate(sigma)
    if(allocated(dsigmadR))        deallocate(dsigmadR)
    if(allocated(dsigmadh))        deallocate(dsigmadh)
    if(allocated(dsigmadV))        deallocate(dsigmadV)

    if(allocated(alpha_0_sl))       deallocate(alpha_0_sl)
    if(allocated(dalpha_0dR_sl))    deallocate(dalpha_0dR_sl)
    if(allocated(dalpha_0dh_sl))    deallocate(dalpha_0dh_sl)

    if(allocated(omega_scs_sl))     deallocate(omega_scs_sl)
    if(allocated(domegadR_sl))      deallocate(domegadR_sl)
    if(allocated(domegadh_sl))      deallocate(domegadh_sl)

    if(allocated(R_MBD_VdW_sl))     deallocate(R_MBD_VdW_sl)
    if(allocated(dR_MBD_VdWdR_sl))  deallocate(dR_MBD_VdWdR_sl)
    if(allocated(dR_MBD_VdWdh_sl))  deallocate(dR_MBD_VdWdh_sl)

    if(allocated(alpha_free))       deallocate(alpha_free)
    if(allocated(C6_free))          deallocate(C6_free)
    if(allocated(R_vdw_free))       deallocate(R_vdw_free)

    if(allocated(nonint_FmbdVdW))   deallocate(nonint_FmbdVdW)
    if(allocated(int_FmbdVdW))      deallocate(int_FmbdVdW)
    if(allocated(FmbdVdW))          deallocate(FmbdVdW)

    if(allocated(nonint_HmbdVdW))   deallocate(nonint_HmbdVdW)
    if(allocated(int_HmbdVdW))      deallocate(int_HmbdVdW)
    if(allocated(HmbdVdW))          deallocate(HmbdVdW)

    if(allocated(nonint_Uprefactor))   deallocate(nonint_Uprefactor)
    if(allocated(int_Uprefactor))      deallocate(int_Uprefactor)
    if(allocated(Uprefactor))          deallocate(Uprefactor)
    if(allocated(UmbdVdW))             deallocate(UmbdvdW)

    if(allocated(unique_pairs))        deallocate(unique_pairs)
    if(allocated(pairs_scs))       deallocate(pairs_scs)

    if(allocated(n_comps_f))        deallocate(n_comps_f)
    if(allocated(n_comps_h))        deallocate(n_comps_h)
    if(allocated(n_comps_v))        deallocate(n_comps_v)

    if(allocated(f_cpu_id))        deallocate(f_cpu_id)
    if(allocated(h_cpu_id))        deallocate(h_cpu_id)
    if(allocated(v_cpu_id))        deallocate(v_cpu_id)
    end subroutine mbdvdw_finalize

  subroutine mbdvdw_initialize(h_in)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This method will initialize the MBD calc by initalizing the associated constants
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer :: i_atom, Ndim
    REAL(DP), INTENT(IN), OPTIONAL :: h_in(3,3)

    CALL start_clock('mbd_init')

    call hirshfeld_initialize()

    dip_cutoff = 180.0_DP
    mbd_first_step = .true.
    mbd_conv_elec = .false.
    supercell_cutoff = mbd_vdw_supercell/0.5291772490_DP
    do_forces = mbd_vdw_forces

    if(mbd_vdw_verbosity.eq.-2) mbd_vdw_verbosity = iverbosity

    do_ewald = mbd_vdw_ewald.and.(.not.mbd_vdw_isolated)
    if(mbd_vdw_low_dim) do_ewald = .false.
    if(mbd_vdw_vacuum(1)) do_ewald = .false.
    if(mbd_vdw_vacuum(2)) do_ewald = .false.
    if(mbd_vdw_vacuum(3)) do_ewald = .false.

    nk = mbd_vdw_kgrid(1)*mbd_vdw_kgrid(2)*mbd_vdw_kgrid(3)
    MBD_PER_ERROR = .false.
    if(.not.mbd_vdw_isolated) then
      if(nk.eq.0) then
        if(supercell_cutoff.eq.0 .and. me_image == root_image) then
          ! Case where neither supercell nor kgrid are set
          write(use_unit, *) "WARNING!! NEITHER A K-GRID NOR A SUPERCELL WERE SET. PLEASE SET ONE"
          MBD_PER_ERROR = .true.
          else
          ! Case where supercell is set but kgrid is not
          do_recip = .false.
          end if
        else
        ! Case where supercell is not set but kgrid is set
        do_recip = .true.
        supercell_cutoff = 1.0_DP
        end if
      end if
    do_cmplx = do_recip

    if(.not.mbd_vdw_isolated) then
        sl_i = ceiling(supercell_cutoff / dsqrt(h_in(1,1)**2.0_DP + h_in(2,1)**2.0_DP + h_in(3,1)**2.0_DP ))
        sl_j= ceiling(supercell_cutoff / dsqrt(h_in(1,2)**2.0_DP + h_in(2,2)**2.0_DP + h_in(3,2)**2.0_DP ))
        sl_k = ceiling(supercell_cutoff / dsqrt(h_in(1,3)**2.0_DP + h_in(2,3)**2.0_DP + h_in(3,3)**2.0_DP ))
      else
        sl_i = 1
        sl_j = 1
        sl_k = 1
      end if
    if(mbd_vdw_vacuum(1)) sl_i = 1
    if(mbd_vdw_vacuum(2)) sl_j = 1
    if(mbd_vdw_vacuum(3)) sl_k = 1

    nat_sl = nat*sl_i*sl_j*sl_k

    if(.not.do_recip) then
      if(mbd_vdw_verbosity.ge.0 .and. me_image == root_image) then
        write(use_unit, *) 'Using a supercell cutoff of', supercell_cutoff
        write(use_unit, *) 'Constructed a superlattice of size ', sl_i, sl_j, sl_k
        write(use_unit, *) 'Number of atoms in the superlattice is', nat_sl
        end if
      end if

    sl_mult = dble(1.0_DP)

    if(.not.allocated(R_vdw_free))      allocate(R_vdw_free(nat))
    if(.not.allocated(R_TS_VdW))        allocate(R_TS_VdW(nat));            R_TS_VdW    = 0.0_DP
    if(.not.allocated(R_MBD_VdW))       allocate(R_MBD_VdW(nat));            R_MBD_VdW    = 0.0_DP
    if(.not.allocated(C6_free))         allocate(C6_free(nat));              C6_free = 0.0_DP
    if(.not.allocated(tau))             allocate(tau(3, nat))
    if(.not.allocated(alpha_free))      allocate(alpha_free(nat))
    if(.not.allocated(alpha_ts))        allocate(alpha_ts(nat));             alpha_ts = 0.0_DP
    if(.not.allocated(alpha_0))         allocate(alpha_0(nat));              alpha_0 = 0.0_DP
    if(.not.allocated(omega_scs))       allocate(omega_scs(nat));            omega_scs = 0.0_DP
    if(.not.allocated(sigma))           allocate(sigma(nat));                sigma = 0.0_DP
    if(.not.allocated(dR_TS_VdWdV))     allocate(dR_TS_VdWdV(nat, nat)); dR_TS_VdWdV = 0.0_DP
    if(.not.allocated(dR_MBD_VdWdV))    allocate(dR_MBD_VdWdV(nat, nat));dR_MBD_VdWdv = 0.0_DP
    if(.not.allocated(domegadV))        allocate(domegadV(nat, nat));    domegadV = 0.0_DP
    if(.not.allocated(dalpha_tsdV))     allocate(dalpha_tsdV(nat, nat)); dalpha_tsdV = 0.0_DP
    if(.not.allocated(dalpha_0dV))      allocate(dalpha_0dV(nat, nat));  dalpha_0dV = 0.0_DP
    if(.not.allocated(dsigmadV))        allocate(dsigmadV(nat, nat));    dsigmadV = 0.0_DP
    if(.not.allocated(dR_TS_VdWdR))     allocate(dR_TS_VdWdR(nat, nat, 3)); dR_TS_VdWdR = 0.0_DP
    if(.not.allocated(dR_TS_VdWdh))     allocate(dR_TS_VdWdh(nat, 3  , 3)); dR_TS_VdWdh = 0.0_DP
    if(.not.allocated(dR_MBD_VdWdR))    allocate(dR_MBD_VdWdR(nat, nat, 3)); dR_MBD_VdWdR = 0.0_DP
    if(.not.allocated(dR_MBD_VdWdh))    allocate(dR_MBD_VdWdh(nat, 3  , 3)); dR_MBD_VdWdh = 0.0_DP
    if(.not.allocated(domegadR))        allocate(domegadR(nat, nat, 3));     domegadR = 0.0_DP
    if(.not.allocated(domegadh))        allocate(domegadh(nat, 3, 3));       domegadh = 0.0_DP
    if(.not.allocated(dalpha_tsdR))     allocate(dalpha_tsdR(nat, nat, 3));  dalpha_tsdR = 0.0_DP
    if(.not.allocated(dalpha_tsdh))     allocate(dalpha_tsdh(nat, 3  , 3));  dalpha_tsdh = 0.0_DP
    if(.not.allocated(dalpha_0dR))      allocate(dalpha_0dR(nat, nat, 3));   dalpha_0dR = 0.0_DP
    if(.not.allocated(dalpha_0dh))      allocate(dalpha_0dh(nat, 3, 3));     dalpha_0dh = 0.0_DP
    if(.not.allocated(dsigmadR))        allocate(dsigmadR(nat, nat, 3));     dsigmadR = 0.0_DP
    if(.not.allocated(dsigmadh))        allocate(dsigmadh(nat, 3  , 3));     dsigmadh = 0.0_DP

    if(.not.allocated(nonint_Uprefactor))  allocate(nonint_Uprefactor(nat));       nonint_Uprefactor = 0.0_DP
    if(.not.allocated(int_Uprefactor))     allocate(int_Uprefactor(nat));          int_Uprefactor = 0.0_DP
    if(.not.allocated(Uprefactor))         allocate(Uprefactor(nat));              Uprefactor = 0.0_DP

    if(.not.allocated(nonint_FmbdVdW))  allocate(nonint_FmbdVdW(nat, 3));    nonint_FmbdVdW = 0.0_DP
    if(.not.allocated(int_FmbdVdW))     allocate(int_FmbdVdW(nat, 3));       int_FmbdVdW = 0.0_DP
    if(.not.allocated(FmbdVdW))         allocate(FmbdVdW(nat, 3));           FmbdVdW = 0.0_DP

    if(.not.allocated(nonint_HmbdVdW))  allocate(nonint_HmbdVdW(3, 3));      nonint_HmbdVdW = 0.0_DP
    if(.not.allocated(int_HmbdVdW))     allocate(int_HmbdVdW(3, 3));         int_HmbdVdW = 0.0_DP
    if(.not.allocated(HmbdVdW))         allocate(HmbdVdW(3, 3));             HmbdVdW = 0.0_DP

    if(.not.allocated(tau_s))           allocate(tau_s(3, nat));            tau_s = 0.0_DP
    if(.not.allocated(n_pairs))         allocate(n_pairs(nproc_image));     n_pairs = 0

    if(.not.allocated(n_comps_f))       allocate(n_comps_f(nproc_image)); n_comps_f = 0
    if(.not.allocated(n_comps_h))       allocate(n_comps_h(nproc_image)); n_comps_h = 0
    if(.not.allocated(n_comps_v))       allocate(n_comps_v(nproc_image)); n_comps_v = 0

    if(.not.allocated(f_cpu_id))        allocate(f_cpu_id(3*nat)); f_cpu_id = 0
    if(.not.allocated(h_cpu_id))        allocate(h_cpu_id(9)); h_cpu_id = 0
    if(.not.allocated(v_cpu_id))        allocate(v_cpu_id(nat)); v_cpu_id = 0
    Ndim=dfftp%nnr
    if(.not.allocated(UmbdvdW)) ALLOCATE(UmbdvdW(Ndim)); UmbdvdW=0.0_DP
    call generate_grid(mbd_vdw_n_quad_pts, npts)

    do i_atom=1, nat, 1
      call GetVdWParam(atm(ityp(i_atom)), C6_free(i_atom), alpha_free(i_atom), R_vdw_free(i_atom))
    end do

    beta = mbd_vdw_beta

    CALL stop_clock('mbd_init')
    end subroutine mbdvdw_initialize

  subroutine mbdvdw_pbc(tauin, hin, ainvin, natin)
    implicit none
    integer,  intent(in):: natin
    integer              :: i_atom
    REAL(DP), intent(inout) :: tauin(3, natin)
    REAL(DP), intent(in) :: hin(3, 3)
    real(dp), intent(in) :: ainvin(3, 3)
    real(dp) :: tau_s(3)

    do i_atom = 1, natin, 1
      tau_s(1)=ainvin(1,1)*tauin(1,i_atom)+ainvin(1,2)*tauin(2,i_atom)+ainvin(1,3)*tauin(3,i_atom)   ! s = h_^-1 r
      tau_s(2)=ainvin(2,1)*tauin(1,i_atom)+ainvin(2,2)*tauin(2,i_atom)+ainvin(2,3)*tauin(3,i_atom)   ! s = h_^-1 r
      tau_s(3)=ainvin(3,1)*tauin(1,i_atom)+ainvin(3,2)*tauin(2,i_atom)+ainvin(3,3)*tauin(3,i_atom)   ! s = h_^-1 r

      ! tau_s(1)=tau_s(1)-floor(tau_s(1))   ! impose PBC on s in range: [0,1)
      ! tau_s(2)=tau_s(2)-floor(tau_s(2))   ! impose PBC on s in range: [0,1)
      ! tau_s(3)=tau_s(3)-floor(tau_s(3))   ! impose PBC on s in range: [0,1)

      tauin(1,i_atom)=hin(1,1)*tau_s(1)+hin(1,2)*tau_s(2)+hin(1,3)*tau_s(3)   ! r = h_ s
      tauin(2,i_atom)=hin(2,1)*tau_s(1)+hin(2,2)*tau_s(2)+hin(2,3)*tau_s(3)   ! r = h_ s
      tauin(3,i_atom)=hin(3,1)*tau_s(1)+hin(3,2)*tau_s(2)+hin(3,3)*tau_s(3)   ! r = h_ s
    end do

    end subroutine mbdvdw_pbc

  subroutine mbdvdw_effqts(omega)
    implicit none
    real(dp), intent(in) :: omega
    real(dp) :: omega_free, gamma, xi, pade_approx, lambda
    integer :: i_atom, s, i, ias

    call start_clock('mbd_effqts')

    if(vdw_self_consistent) then
      dsigmadV = 0.0_DP
      dR_TS_VdWdV = 0.0_DP
      dalpha_tsdV = 0.0_DP
      end if

    if(do_forces) then
      dsigmadR = 0.0_DP
      dR_TS_VdWdR = 0.0_DP
      dalpha_tsdR = 0.0_DP
      dsigmadH = 0.0_DP
      dR_TS_VdWdH = 0.0_DP
      dalpha_tsdH = 0.0_DP
      end if

    do i_atom=1, nat, 1
      ias = ityp(i_atom)
      omega_free = ((4.0_DP/3.0_DP)*C6_free(i_atom)/(alpha_free(i_atom)**2.0_DP))
      ! Pade Approx
      pade_approx = 1.0_DP/(1.0_DP + (omega/omega_free)**2.0_DP )
      ! Computes sigma
      gamma = (1.0_DP/3.0_DP)*dsqrt(2.0_DP/pi)*pade_approx*(alpha_free(i_atom)/vfree(ias))
      gamma = gamma**(1.0_DP/3.0_DP)
      sigma(i_atom) = gamma*VefftsvdW(i_atom)**(1.0_DP/3.0_DP)

      ! Computes R_TS
      xi = R_vdw_free(i_atom)/(vfree(ias))**(1.0_DP/3.0_DP)
      R_TS_VdW(i_atom) = xi*VefftsvdW(i_atom)**(1.0_DP/3.0_DP)

      ! Computes alpha_ts
      lambda = pade_approx*alpha_free(i_atom)/vfree(ias)
      alpha_ts(i_atom) = lambda*VefftsvdW(i_atom)

      if(vdw_self_consistent) then
        do s = 1, nat, 1
          if(i_atom.eq.s) then
            dsigmadV(i_atom, s) = gamma/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))
            dR_TS_VdWdV(i_atom, s) = xi/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))
            dalpha_tsdV(i_atom, s) = lambda
            end if
          end do ! s loop
        end if ! Test for vdw_self_consistent
      if(do_forces) then
        ! Loop for dR
        do s = 1, nat, 1
          do i = 1, 3, 1
            dsigmadR(i_atom, s, i) = gamma/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdR(i_atom, s, i)
            dR_TS_VdWdR(i_atom, s, i) = xi/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdR(i_atom, s, i)
            dalpha_tsdR(i_atom, s, i) = lambda*dveffdR(i_atom, s, i)
            end do ! s loop
          end do ! i loop
        ! Loop for dH
        if(.not.mbd_vdw_isolated) then
          do s = 1, 3, 1
            do i = 1, 3, 1
              dsigmadh(i_atom, s, i) = gamma/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdh(i_atom, s, i)
              dR_TS_VdWdh(i_atom, s, i) = xi/(3.0_DP*VefftsvdW(i_atom)**(2.0_DP/3.0_DP))*dveffdh(i_atom, s, i)
              dalpha_tsdh(i_atom, s, i) = lambda*dveffdh(i_atom, s, i)
              end do ! S loop
            end do ! i loop
          end if ! Test to see if we need to compute stresses
        end if ! Test for doing forces
      end do ! i_atom loop
    call stop_clock('mbd_effqts')

    end subroutine mbdvdw_effqts

  subroutine mbdvdw_para_init()
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initializes the parallel part of the mbd computation ! TODO: rewrite
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    integer :: p, q, counter
    integer :: n_pairs_per_proc, cpu

    num_pairs = (nat*nat - nat)/2 + nat
    n_pairs = 0
    if(.not.allocated(pairs_scs))        allocate(pairs_scs(num_pairs))
    n_pairs_per_proc = floor(dble(num_pairs)/nproc_image)
    cpu = 0
    n_pairs = 0
    counter = 1

    !HK: TODO (comment 1)
    do p = 1, nat, 1
      do q = p, nat, 1
        pairs_scs(counter)%p = p
        pairs_scs(counter)%q = q
        counter = counter + 1
      end do
    end do

    do counter = 0, num_pairs-1, 1
      n_pairs(modulo(counter, nproc_image)+1) = n_pairs(modulo(counter, nproc_image)+1) + 1
    end do

    !HK: TODO (comment 1)
    p = 1
    do counter = 1, num_pairs, 1
      pairs_scs(counter)%cpu = cpu
      if((counter.lt.num_pairs)) then
        if(p.eq.n_pairs(cpu+1)) then
          cpu = cpu + 1
          p=0
        end if
      end if
      p=p+1
    end do

    if(num_pairs.ge.nproc_image) then
      max_proc_pairs = nproc_image
    else
      max_proc_pairs = num_pairs
    end if

    RETURN
    end subroutine mbdvdw_para_init

  subroutine mbdvdw_para_init_sl()
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initializes the parallel part of the mbd computation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer :: p, q, counter
    integer :: n_pairs_per_proc, cpu

    num_pairs = (nat_sl*nat_sl-nat_sl)/2 + nat_sl
    n_pairs = 0
    if(.not.allocated(unique_pairs))        allocate(unique_pairs(num_pairs))
    n_pairs_per_proc = floor(dble(num_pairs)/nproc_image)
    cpu = 0
    n_pairs = 0
    counter = 1

    !HK: TODO (comment 1)
    do p = 1, nat_sl, 1
      do q = p, nat_sl, 1
        unique_pairs(counter)%p = p
        unique_pairs(counter)%q = q
        counter = counter + 1
      end do
    end do

    do counter = 0, num_pairs-1, 1
      n_pairs(modulo(counter, nproc_image)+1) = n_pairs(modulo(counter, nproc_image)+1) + 1
    end do

    !HK: TODO (comment 1)
    p = 1
    do counter = 1, num_pairs, 1
      unique_pairs(counter)%cpu = cpu
      if((counter.lt.num_pairs)) then
        if(p.eq.n_pairs(cpu+1)) then
          cpu = cpu + 1
          p=0
          end if
        end if
        p=p+1
      end do

    if(num_pairs.ge.nproc_image) then
      max_proc_pairs = nproc_image
      else
      max_proc_pairs = num_pairs
      end if

    RETURN

    end subroutine mbdvdw_para_init_sl

  subroutine mbdvdw_para_init_forces(n, cpu_id, ncomps, max_proc)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initializes the parallel part of the mbd computation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: cpu_id(n), ncomps(nproc_image), max_proc
    integer :: p, s, counter
    integer :: cpu

    cpu = 0
    ncomps = 0

    counter = 1
    p = 1

    do counter = 0, n-1, 1
      ncomps(modulo(counter, nproc_image)+1) = ncomps(modulo(counter, nproc_image)+1) + 1
    end do

    do s = 1, n, 1
      cpu_id(s) = cpu
      if((s.lt.n).and.(p.eq.ncomps(cpu+1))) then
        cpu = cpu + 1
        p = 0
        end if
      p = p + 1
      end do

    if(n.ge.nproc_image) then
      max_proc = nproc_image
      else
      max_proc = n
      end if

    end subroutine mbdvdw_para_init_forces

  subroutine mbdvdw_noninteracting_energy(energy, forcedR, forcedh, forcedV)
    implicit none
    real(dp), intent(out) :: energy
    real(dp), dimension(nat, 3), intent(out) :: forcedR
    real(dp), dimension(3,3), intent(out) :: forcedh
    real(dp), dimension(nat), intent(out) :: forcedV
    integer :: p
    ! Local variable

    call start_clock('mbd_nonint_energy')
    energy = 0.0_DP
    forcedR = 0.0_DP
    forcedh = 0.0_DP
    forcedV = 0.0_DP

    !$omp parallel do reduction(+:energy,forcedR,forcedh,forcedV)
    do p=1, nat, 1
      energy = energy + omega_scs(p)
      if(do_forces) then
        forcedR = forcedR + domegadR(p, :, :)
        forcedh = forcedh + domegadh(p, :, :)
        end if
      if(vdw_self_consistent) forcedV(:) = forcedV(:) + domegadV(p, :)
    end do
    !$omp end parallel do

    if(vdw_self_consistent) call mp_sum(forcedV, intra_image_comm)
    if(mbd_vdw_forces) call mp_sum(forcedR, intra_image_comm)
    if(mbd_vdw_forces) call mp_sum(forcedh, intra_image_comm)

    call stop_clock('mbd_nonint_energy')
    return
    end subroutine mbdvdw_noninteracting_energy

  subroutine mbdvdw_interacting_energy(energy, forcedR, forcedh, forcedV)
    implicit none
    real(dp), intent(out) :: energy
    real(dp), dimension(nat, 3), intent(out) :: forcedR
    real(dp), dimension(3,3), intent(out) :: forcedh
    real(dp), dimension(nat), intent(out) :: forcedV
    real(dp), dimension(3*nat_sl, 3*nat_sl) :: temp
    integer :: num_negative, i_atom, s, i, j, i_f, cnt_v, cnt_h, cnt_f

    ! lapack work variables
    integer :: LWORK, errorflag
    real(dp) :: WORK((3*nat_sl)*(3+(3*nat_sl)/2)), eigenvalues(3*nat_sl), eval_inv(3*nat_sl)

    call start_clock('mbd_int_energy')
    eigenvalues = 0.0_DP
    forcedR = 0.0_DP
    forcedh = 0.0_DP
    energy = 0.0_DP
    num_negative = 0
    forcedV = 0.0_DP

    errorflag=0
    LWORK=3*nat_sl*(3+(3*nat_sl)/2)
    call DSYEV('V', 'U', 3*nat_sl, Cpq, 3*nat_sl, eigenvalues, WORK, LWORK, errorflag)

    if(errorflag.eq.0) then
      do i_atom=1, 3*nat_sl, 1
        if(eigenvalues(i_atom).ge.0.0_DP) then
          eigenvalues(i_atom) = dsqrt(eigenvalues(i_atom))
          eval_inv(i_atom) = 1.0_DP/(2.0_DP*eigenvalues(i_atom))
          energy = energy + eigenvalues(i_atom)
        else
          num_negative = num_negative + 1
        end if
      end do

      if(num_negative.ge.1 .and. me_image == root_image) then
        write(use_unit, '(3X," WARNING: Found ", I3, " Negative Eigenvalues.")') num_negative
        end if
      if (num_negative > mbd_max_negative_eigvals) then
        if (me_image == root_image) write (use_unit, '(A)') neg_eigvals_msg
        energy = 0.0_DP
        forcedR = 0.0_DP
        forcedh = 0.0_DP
        forcedV = 0.0_DP
        return
      end if
    end if
    energy = energy*nat/nat_sl

    if(vdw_self_consistent) then
      cnt_v = 1
      do i_f = 1, nat, 1
        if(v_cpu_id(i_f).eq.me_image) then
          call mult_TNN(3*nat_sl, Cpq, dCpqdV(:, :, cnt_v), temp)
          do j=1,3*nat_sl,1
            if(eigenvalues(j).ge.0.0_DP) forcedV(i_f) = forcedV(i_f) + eval_inv(j)*temp(j,j)
            end do
          cnt_v = cnt_v + 1
          end if ! test if this mpi_rank should do this calc
        end do ! i_f loop for dv
      forcedV = forcedV*nat/nat_sl
      end if ! do forces test

    if(do_forces) then
      ! dR
      cnt_f = 1
      do i_f = 1, 3*nat, 1
        if(f_cpu_id(i_f).eq.me_image) then
          call mbdvdw_get_is(i_f, s, i)
          call mult_TNN(3*nat_sl, Cpq, dCpqdR(:, :, cnt_f), temp)
          do j=1,3*nat_sl,1
            if(eigenvalues(j).ge.0.0_DP) forcedR(s, i) = forcedR(s, i) + eval_inv(j)*temp(j,j)
            end do
          cnt_f = cnt_f + 1
          end if ! tests if this mpi_rank should do this force calc
        end do ! i_f loop for dR
      forcedR = forcedR*nat/nat_sl
      ! dh
      if(.not.mbd_vdw_isolated) then
        cnt_h = 1
        do i_f = 1, 9, 1
          if(h_cpu_id(i_f).eq.me_image) then
            call mbdvdw_get_is(i_f, s, i)
            call mult_TNN(3*nat_sl, Cpq, dCpqdh(:, :, cnt_h), temp)
            do j=1,3*nat_sl,1
              if(eigenvalues(j).ge.0.0_DP) forcedh(s, i) = forcedh(s, i) + eval_inv(j)*temp(j,j)
              end do
            cnt_h = cnt_h + 1
            end if  ! tests if this mpi_rank should do this force calc
          end do ! i_f loop for dh
        end if ! test to check if we should compute stresses
      forcedh = forcedh*nat/nat_sl
      end if ! do forces test

    if(vdw_self_consistent) call mp_sum(forcedV, intra_image_comm)
    if(do_forces) call mp_sum(forcedR, intra_image_comm)
    if(do_forces.and.(.not.mbd_vdw_isolated)) call mp_sum(forcedh, intra_image_comm)

    call stop_clock('mbd_int_energy')
    return
    end subroutine mbdvdw_interacting_energy

  subroutine mbdvdw_interacting_energy_complex(energy, forcedR, forcedh, forcedV)
    implicit none
    real(dp), intent(out) :: energy
    real(dp), dimension(nat, 3), intent(out) :: forcedR
    real(dp), dimension(3,3), intent(out) :: forcedh
    real(dp), dimension(nat), intent(out) :: forcedV
    complex(dp), dimension(3*nat_sl, 3*nat_sl) :: temp
    integer :: num_negative, i_atom, s, i, j, i_f, cnt_v, cnt_h, cnt_f, k, nmax, ik

    ! lapack work variables
    integer :: lwork, errorflag, do_k(nk), did_k(nk)
    real(dp) :: eigenvalues(3*nat_sl, nk), eval_inv(3*nat_sl, nk), remainder, energy_arr(nk)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    complex(dp) :: lwork_cmplx

    call start_clock('mbd_int_energy')
    eigenvalues = 0.0_DP
    forcedR = 0.0_DP
    forcedh = 0.0_DP
    energy = 0.0_DP
    num_negative = 0
    forcedV = 0.0_DP
    energy_arr = 0.0_dp
    do_k = 0
    did_k = 0

    nmax = nk
    if(vdw_self_consistent) then
      do i_f = 1, nat*nk, 1
        call mbdvdw_get_ks(i_f, s, k)
        do_k(k) = 1
        end do
      end if
    if(do_forces) then
      do i_f = 1, 3*nat*nk, 1
        call mbdvdw_get_isk(i_f, s, i, k)
        do_k(k) = 1
        end do
      if(.not.mbd_vdw_isolated) then
        do i_f = 1, 9*nk, 1
          call mbdvdw_get_isk(i_f, s, i, k)
          do_k(k) = 1
          end do
        end if
      end if

    if((.not.do_forces).and.(.not.vdw_self_consistent)) then
      do i_f = 1, nk, 1
        if(k_cpu_id(i_f).eq.me_image) then
          do_k(i_f) = 1
          end if
        end do
      end if

    errorflag=0
    if(.not.allocated(rwork)) allocate(rwork(max(1, 3*(3*nat_sl)-2)))
    call ZHEEV('V', 'U', 3*nat_sl, Cpq_c(:,:,1), 3*nat_sl, eigenvalues(:, 1), lwork_cmplx, -1, rwork, errorflag)
    lwork = nint(dble(lwork_cmplx))
    if(.not.allocated(work)) allocate (work(lwork))

    eigenvalues = 0.0_dp
    num_negative = 0
    do ik = 1, nk, 1
      if(do_k(ik).eq.1) then
        work = 0.0_dp
        rwork = 0.0_dp
        call ZHEEV('V', 'U', 3*nat_sl, Cpq_c(:,:,ik), 3*nat_sl, eigenvalues(:,ik), work, lwork, rwork, errorflag)
        do i_atom = 1, 3*nat_sl, 1
          if(eigenvalues(i_atom, ik).gt.0.0_DP) then
            eigenvalues(i_atom, ik) = dsqrt(eigenvalues(i_atom, ik))
            eval_inv(i_atom, ik) = 1.0_DP/(2.0_DP*eigenvalues(i_atom, ik))
            else
            num_negative = num_negative + 1
            end if
          end do
        end if
      end do
    if(num_negative.ge.1 .and. me_image == root_image) then
      write(use_unit, '(3X," WARNING: Found ", I3, " Negative Eigenvalues.")') num_negative
      end if
      if (num_negative > mbd_max_negative_eigvals) then
        if (me_image == root_image) write (use_unit, '(A)') neg_eigvals_msg
        energy = 0.0_DP
        forcedR = 0.0_DP
        forcedh = 0.0_DP
        forcedV = 0.0_DP
        return
      end if
    do i_f = 1, nk, 1
      if(do_k(i_f).eq.1) then
        did_k(i_f) = 1
        do i_atom=1, 3*nat_sl, 1
          if(eigenvalues(i_atom, i_f).gt.0.0_DP) then
            energy_arr(i_f) = energy_arr(i_f) + eigenvalues(i_atom, i_f)
            end if
          end do
        end if
      end do

    if(allocated(work)) deallocate(work)
    if(allocated(rwork)) deallocate(rwork)

    remainder = 0.0_dp

    if(vdw_self_consistent) then
      cnt_v = 1
      do i_f = 1, nat*nk, 1
        if(v_cpu_id(i_f).eq.me_image) then
          call mbdvdw_get_ks(i_f, s, k)
          call mult_CNN(3*nat_sl, Cpq_c(:,:,k), dCpqdV_c(:, :, cnt_v), temp)
          do j=1,3*nat_sl,1
            if(eigenvalues(j, k).gt.0.0_DP) forcedV(s) = forcedV(s) + eval_inv(j, k)*real(temp(j,j))
            end do
          cnt_v = cnt_v + 1
          end if ! statement to test if this self consistent component should be computed
        end do ! Loop over self consistent components
      forcedV = forcedV*nat/nat_sl/nk
      end if ! Test to see if self consistent forces should be done

    if(do_forces) then
      ! dR
      cnt_f = 1
      do i_f = 1, 3*nat*nk, 1
        if(f_cpu_id(i_f).eq.me_image) then
          call mbdvdw_get_isk(i_f, s, i, k)
          call mult_CNN(3*nat_sl, Cpq_c(:,:,k), dCpqdR_C(:, :, cnt_f), temp)
          do j=1,3*nat_sl,1
            if(eigenvalues(j, k).gt.0.0_DP) then
              forcedR(s, i) = forcedR(s, i) + eval_inv(j, k)*real(temp(j,j))
              remainder = remainder + imag(temp(j,j))
              end if
            end do
          cnt_f = cnt_f + 1
          end if ! statement to test if this forces should be computed
        end do
      forcedR = forcedR*nat/nat_sl/nk
      ! dH
      if(.not.mbd_vdw_isolated) then
        cnt_h = 1
        do i_f = 1, 9*nk, 1
          if(h_cpu_id(i_f).eq.me_image) then
            call mbdvdw_get_isk(i_f, s, i, k)
            call mult_CNN(3*nat_sl, Cpq_c(:,:,k), dCpqdh_C(:, :, cnt_h), temp)
            do j=1,3*nat_sl,1
              if(eigenvalues(j, k).gt.0.0_DP) then
                forcedh(s, i) = forcedh(s, i) + eval_inv(j, k)*real(temp(j,j))
                remainder = remainder + imag(temp(j,j))
                end if
              end do
            cnt_h = cnt_h + 1
            end if ! statement to test if this forces should be computed
          end do
        end if
      forcedh = forcedh*nat/nat_sl/nk
      end if ! Test to see if forces and stresses should be calculated

    call mp_sum(energy_arr, intra_image_comm)
    call mp_sum(did_k, intra_image_comm)

    do i_f = 1, nk, 1
      energy = energy + energy_arr(i_f)/did_k(i_f)
      end do
    energy = energy*nat/nat_sl
    energy = energy/nk

    if(vdw_self_consistent) call mp_sum(forcedV, intra_image_comm)
    if(do_forces) call mp_sum(forcedR, intra_image_comm)
    if(do_forces.and.(.not.mbd_vdw_isolated)) call mp_sum(forcedh, intra_image_comm)

    call stop_clock('mbd_int_energy')
    return
    end subroutine mbdvdw_interacting_energy_complex


  subroutine mbdvdw_calculate_screened_pol(alpha_isotropic, dalpha_isotropicdR, dalpha_isotropicdh, dalpha_isotropicdV)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This method simply performs the trace of the A matrix found in the SCS
    ! procedure. The equation is
    ! $\bar{\alpha}_p(i \omega) = \frac{1}{3} Tr \left[ \sum_q \bar{A}_{pq}(i \omega) \right]$
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    ! variables to be returned
    real(dp), dimension(nat), intent(out) :: alpha_isotropic
    real(dp), dimension(nat, nat, 3), intent(out) :: dalpha_isotropicdR
    real(dp), dimension(nat, 3  , 3), intent(out) :: dalpha_isotropicdh
    real(dp), dimension(nat, nat), intent(out) :: dalpha_isotropicdV
    ! local variables
    integer :: p_atom, q_atom, i_index, j_index, s, i, i_idx, cnt_f, cnt_h, cnt_v, i_f
    call start_clock('mbd_calculate_screened_pol')

    alpha_isotropic = 0.0_DP
    if(vdw_self_consistent) dalpha_isotropicdV = 0.0_DP
    if(do_forces) dalpha_isotropicdR = 0.0_DP
    if(do_forces) dalpha_isotropicdh = 0.0_DP

    ! This set of loops goes over all atoms, and traces along the rows. We'll create a 3x3 matrix that takes
    ! care to preserve the directional information. Then, we'll compute the eigenvalues of the 3x3
    ! matrix, sum them, and throw them into the isotropic alpha.

    !HK: TODO omp-lize progress line
    do p_atom=1, nat, 1
      do q_atom=1, nat, 1
        do i_idx=1,3,1
          cnt_v = 1
          cnt_h = 1
          cnt_f = 1
          i_index = (3*p_atom - 3)
          j_index = (3*q_atom - 3)

          if(vdw_self_consistent) then
            do i_f = 1, nat, 1
              if(v_cpu_id(i_f).eq.me_image) then
                dalpha_isotropicdV(p_atom, i_f) = dalpha_isotropicdV(p_atom, i_f)+dA_matrixdV(i_index+i_idx, j_index+i_idx, cnt_v)
                cnt_v = cnt_v + 1
                end if ! checks if this MPI rank computes this component
              end do ! Loop over all possible components
            end if ! if to check if self consistent derivatives need to be computed

          if(do_forces) then
            ! dR
            do i_f = 1, 3*nat, 1
              if(f_cpu_id(i_f).eq.me_image) then
                call mbdvdw_get_is(i_f, s, i)
                dalpha_isotropicdR(p_atom, s, i)=dalpha_isotropicdR(p_atom, s, i)+dA_matrixdR(i_index+i_idx, j_index+i_idx, cnt_f)
                cnt_f = cnt_f + 1
                end if ! checks if this MPI rank computes this component
              end do ! Loop over all possible components
            ! dH
            if(.not.mbd_vdw_isolated) then
              do i_f = 1, 9, 1
                if(h_cpu_id(i_f).eq.me_image) then
                  call mbdvdw_get_is(i_f, s, i)
                  dalpha_isotropicdh(p_atom, s, i)=dalpha_isotropicdh(p_atom, s, i)+dA_matrixdh(i_index+i_idx, j_index+i_idx, cnt_h)
                  cnt_h = cnt_h + 1
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              end if
            end if ! if to check if forces and stresses need to be computed

          end do
        do i_idx=1,3,1
          i_index = (3*p_atom - 3)
          j_index = (3*q_atom - 3)
          alpha_isotropic(p_atom)=alpha_isotropic(p_atom)+A_matrix(i_index+i_idx, j_index+i_idx)
          end do
        end do
      end do

    alpha_isotropic = alpha_isotropic/3.0_DP
    if(do_forces) dalpha_isotropicdR = dalpha_isotropicdR/3.0_DP
    if(do_forces) dalpha_isotropicdh = dalpha_isotropicdh/3.0_DP
    if(vdw_self_consistent) dalpha_isotropicdV = dalpha_isotropicdV/3.0_DP

    call stop_clock('mbd_calculate_screened_pol')
    return
    end subroutine mbdvdw_calculate_screened_pol

  subroutine mbdvdw_compute_TSR(p, q, RPQ, Spq_lat, TSR, dTSRdR, dTSRdh, dTSRdV)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This method is going to compute TSR and its relevant derivatives.
    ! The formula in latex for TSR is:
    ! T^{ij}_{SR} = \left[  1  - f^{TS}_{damp} \right] T^{ij}
    ! with
    ! T^{ij} = T_{dip}^{ij} \left[  erf\left[ \frac{R}{\sigma} \right] - \frac{2 R}{\sqrt{\pi} \sigma} e^{- \frac{R^2}{\sigma^2} } \right] + \frac{4}{\sqrt{\pi}} \frac{R^iR^j}{\sigma^3 R^2} e^{-\frac{R^2}{\sigma^2}}
    ! and
    ! T^{ij}_{dip} = \frac{-3 R^i R^j + R^2 \delta_{ij}}{ R^5 }
    !
    ! and \partial T_{SR} is given by
    ! \partial \mathbf{T}^{ij}_{\rm SR} &=& \mathbf{T}^{ij} \; \partial f\left(Z^{\rm TS}\right) + \left[1- f\left(Z^{\rm TS}\right) \right] \partial \mathbf{T}^{ij},
    ! \partial \mathbf{T}^{ij} &=& -3 \left[{\rm erf}\left[\zeta\right]   -\frac{h_(\zeta)}{2\zeta}   \right]  \partial \mathbf{T}_{\rm dip}^{ij}  \\
    ! &&  - \frac{1}{3} \zeta\, h_(\zeta)  \left[   \partial \mathbf{T}_{\rm dip}^{ij} - \frac{\delta_{ij}}{R^4} \partial R \right]    \nonumber \\
    ! && +  \left[ \mathbf{T}_{\rm dip}^{ij} + \frac{R^iR^j}{ R^5} \Big[ 3  - 2\zeta^2  \Big]\right]  h_(\zeta)\,  \partial \zeta,   \nonumber
    !
    ! Yeah... that's super ugly.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer, intent(in) :: p, q
    ! Variables needed for Tsr
    real(dp), dimension(3, 3), intent(out):: TSR
    real(dp) :: Rpq_norm, Spq, Sigma_pq, R_vdw_pq, Z, fermi_fn, zeta, U, W, gaussian, dfn_pre
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    real(dp), dimension(3,3) :: Tdip, Rmat, TGG

    ! Variables needed for dTsrdR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTSRdR
    real(dp) :: dRpq_normdR, dsigma_pqdR, dZetadR, dZdR, dSpqdR, dFermi_fndR, dUdR, dWdR
    real(dp), dimension(3) :: dRpqdR
    real(dp), dimension(3, 3) :: dRmatdR, dTdipdR

    ! Variables needed for dTsrdh
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTSRdh
    real(dp) :: dRpq_normdH, dsigma_pqdh, dZetadh, dZdh, dSpqdh, dFermi_fndh, dUdh, dWdh
    real(dp), dimension(3) :: dRpqdh
    real(dp), dimension(3, 3) :: dRmatdh, dTdipdh

    ! Variables needed for dTsrdV
    real(dp), dimension(3, 3, nat), intent(out) :: dTSRdV
    real(dp) :: dsigma_pqdV, dZetadV, dZdV, dSpqdV, dFermi_fndV, dUdV, dWdV

    ! Loop variables
    integer :: i, j, s

    real(dp) :: derf

    ! Computes the cartesian distance from the vector distance
    Rpq_norm = dsqrt(Rpq(1) **2.0_DP + Rpq(2)**2.0_DP + Rpq(3)**2.0_DP)

    ! Computes the effective correlation length of the interaction potential
    ! defined from the widths of the QHO Gaussians
    Sigma_pq = dsqrt(sigma(p)**2.0_DP + sigma(q)**2.0_DP)
    ! Computes the damping radius
    R_VdW_pq = R_TS_VdW(p) + R_TS_VdW(q)
    Spq = beta*R_VdW_pq
    Z = 6.0_DP*(Rpq_norm/Spq - 1.0_DP)
    zeta = Rpq_norm/Sigma_pq
    fermi_fn = 1.0_DP
    dfn_pre = 0.0_dp

    ! computes the fermi damping function. The latex for this is
    ! f_{damp}(R_{pq}) = \frac{1}{ 1 + exp( - Z(R_{pq}) ) }
    ! where Z = 6 \left( \frac{R_{pq} }{ S_{pq}} - 1 \right)
    ! and S_{pq} = \beta \left(  R_{p, VdW} + R_{q, VdW} \right)
    if(Z.le.35.0_DP) then
      fermi_fn = 1.0_DP/(1.0_DP+exp( -Z))
      dfn_pre = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP
      ! Computes the factors for U
      ! U = {\rm erf}\left[\zeta\right] -  \frac{2}{\sqrt{\pi}}\zeta \exp\left[-\zeta^2\right]
      zeta = Rpq_norm/Sigma_pq
      end if

    if(zeta.ge.6.0_DP) then
      U = 1.0_DP
      W = 0.0_DP
      gaussian = 0.0_dp
      else
      gaussian = dexp(-zeta*zeta)
      U = derf(zeta)-(2.0_DP*zeta)/SQRTPI*gaussian
      ! Computes the first half of the factor for W before we multiply by the R^i R^j tensor
      ! \mathbf{W}^{ij} &\equiv&   \left(\frac{R^i R^j}{R^5}\right) \, \frac{4}{\sqrt{\pi}}  \zeta^3  \exp\left[-\zeta^2\right]
      W = 4.0_DP*zeta**3.0_DP/SQRTPI*gaussian
      end if

    ! Loops over the cartesian coordinates to compute the R^i*R^j
    ! matrix that goes in to constructing the dipole matrix
    do i=1,3,1
      do j=1,3,1
        Rmat(i, j) = Rpq(i)*Rpq(j)/(Rpq_norm**5.0_DP)
        Tdip(i, j) = -3.0_DP*Rmat(i, j)
        end do
      ! This just applies the kronheker delta center for the dipole tensor. Recall
      ! that T^{ij}_{dip} = \frac{-3 R^i R^j + R^2 \delta_{ij}}{ R^5 }
      Tdip(i, i) = Tdip(i, i) + 1.0_DP/(Rpq_norm**3.0_DP)
      end do

    ! Computes the short range dipole coupling quantity using the fact that
    ! T = T_{dip}\left[ U \right] + W
    TGG = (Tdip*U + Rmat*W)
    TSR = (1.0_DP - fermi_fn)*TGG

    if(vdw_self_consistent) then
      do s = 1, nat, 1
        dsigma_pqdV = (sigma(p)*dsigmadV(p, s) + sigma(q)*dsigmadV(q, s))/Sigma_pq
        dSpqdV = beta*(dR_TS_VdWdV(p, s) + dR_TS_VdWdV(q, s))
        dZetadV = -Rpq_norm*dsigma_pqdV/Sigma_pq**2.0_DP
        dZdV = 6.0_DP*( -Rpq_norm*dSpqdV/Spq**2.0_DP)
        dFermi_fndV = dfn_pre*dZdV

        if(zeta.ge.6.0_DP) then
          dUdV = 0.0_DP
          dWdV = 0.0_DP
          else
          dUdV = dZetadV*zeta*zeta*4.0_DP*gaussian/SQRTPI
          dWdV = 4.0_DP/SQRTPI*zeta*zeta*gaussian*dzetadV*(3.0_DP - 2.0_DP*zeta*zeta)
          end if
        dTSRdV(:, :, s) = -TGG*dFermi_fndV + (1.0_DP - fermi_fn)*(Tdip*dUdV+Rmat*dWdV)
        end do ! Loop over all possible components
      end if ! if to check if self consistent derivatives need to be computed

    if(do_forces) then
      ! dR
      do s = 1, nat, 1
        do i = 1, 3, 1
          call mbdvdw_compute_dRdR(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dRmatdR, dTdipdR)
          dsigma_pqdR = (sigma(p)*dsigmadR(p, s, i) + sigma(q)*dsigmadR(q, s, i))/Sigma_pq
          dZetadR = dRpq_normdR/Sigma_pq - Rpq_norm*dsigma_pqdR/Sigma_pq**2.0_DP
          dSpqdR = beta*(dR_TS_VdWdR(p, s, i) + dR_TS_VdWdR(q, s, i))
          dZdR = 6.0_DP*( dRpq_normdR/Spq - Rpq_norm*dSpqdR/Spq**2.0_DP)
          dFermi_fndR = dfn_pre*dZdR
          if(zeta.ge.6.0_DP) then
            dUdR = 0.0_DP
            dWdR = 0.0_DP
            else
            dUdR = dZetadR*zeta*zeta*4.0_DP*gaussian/SQRTPI
            dWdR = 4.0_DP/SQRTPI*zeta*zeta*gaussian*dzetadR*(3.0_DP - 2.0_DP*zeta*zeta)
            end if
          dTSRdr(:, :, s, i) = -TGG*dFermi_fndR + (1.0_DP - fermi_fn)*(dTdipdR*U + Tdip*dUdR + W*dRmatdR + Rmat*dWdR)
          end do ! loop over cartesian components
        end do ! Loop over all atoms
      ! dH
      if(.not.mbd_vdw_isolated) then
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdh(s, i, Rpq, Rpq_norm, Spq_lat, dRpqdh, dRpq_normdh, dRmatdh, dTdipdh)
            dsigma_pqdH = (sigma(p)*dsigmadH(p, s, i) + sigma(q)*dsigmadH(q, s, i))/Sigma_pq
            dZetadH = dRpq_normdH/Sigma_pq - Rpq_norm*dsigma_pqdH/Sigma_pq**2.0_DP
            dSpqdH = beta*(dR_ts_VdWdH(p, s, i) + dR_ts_VdWdH(q, s, i))
            dZdH = 6.0_DP*( dRpq_normdH/Spq - Rpq_norm*dSpqdH/Spq**2.0_DP)
            dFermi_fndH = dfn_pre*dZdH
            if(zeta.ge.6.0_DP) then
              dUdH = 0.0_DP
              dWdH = 0.0_DP
              else
              dUdH = dZetadH*zeta*zeta*4.0_DP*gaussian/SQRTPI
              dWdH = 4.0_DP/SQRTPI*zeta*zeta*gaussian*dzetadH*(3.0_DP - 2.0_DP*zeta*zeta)
              end if
            dTSRdH(:, :, s, i) = -TGG*dFermi_fndH + (1.0_DP - fermi_fn)*(dTdipdH*U + Tdip*dUdH + W*dRmatdH + Rmat*dWdH)
            end do ! loop over all cell vector components
          end do ! Loop over all cell vectors
        end if
      end if ! if to check if forces and stresses need to be computed

    end subroutine mbdvdw_compute_TSR

  subroutine mbdvdw_TGG(fn, p, q, n, h_in, ainv_in, tau_in, T, dTdR, dTdh, dTdV, deb)
    integer, intent(in) :: p, q, n, fn
    real(dp), intent(in) :: h_in(3,3), tau_in(3,n), ainv_in(3,3)
    integer, intent(in), optional :: deb
    real(dp), dimension(3, 3), intent(out):: T
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTdR
    real(dp), dimension(3, 3, nat), intent(out) :: dTdV
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTdh
    real(dp), dimension(3) :: Rpq, Rpq_lat
    real(dp), dimension(3) :: Spq_lat, Spq
    real(dp), dimension(3, 3) :: T_temp
    real(dp), dimension(3, 3, nat, 3) :: dTdR_temp
    real(dp), dimension(3, 3, nat)  :: dTdV_temp
    real(dp), dimension(3, 3, 3, 3) :: dTdh_temp
    real(dp) :: rnorm, Rc
    integer :: debug, sc(3)
    integer :: n1, n2, n3

    debug = 0
    if(present(deb)) debug = 1

    T = 0.0_DP
    dTdh = 0.0_DP
    dTdR = 0.0_DP
    dTdV = 0.0_DP

    ! Computes the vector distance between two atoms
    Rpq = tau_in(:, p) - tau_in(:, q)
    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)

    ! Spq(1) = Spq(1) - anint(Spq(1))
    ! Spq(2) = Spq(2) - anint(Spq(2))
    ! Spq(3) = Spq(3) - anint(Spq(3))

    if(fn.eq.1) Rc = 20.0_dp ! 20 \AA in Bohr
    if(fn.eq.2) Rc = (2.0_dp/mbd_vdw_econv_thr)**(1.0_dp/3.0_dp)

    Rpq = tau_in(:, p) - tau_in(:, q)
    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)
    sc = mbdvdw_circumscribe(h_in, Rc)
    do n1 = -sc(1), sc(1), 1
      do n2 = -sc(2), sc(2), 1
        do n3 = -sc(3), sc(3), 1
          if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0)) then
            if(p.eq.q) cycle
            end if
          Spq_lat(1) = Spq(1) + DBLE(n1)
          Spq_lat(2) = Spq(2) + DBLE(n2)
          Spq_lat(3) = Spq(3) + DBLE(n3)

          Rpq_lat(1) = h_in(1,1)*Spq_lat(1) + h_in(1,2)*Spq_lat(2) + h_in(1,3)*Spq_lat(3)
          Rpq_lat(2) = h_in(2,1)*Spq_lat(1) + h_in(2,2)*Spq_lat(2) + h_in(2,3)*Spq_lat(3)
          Rpq_lat(3) = h_in(3,1)*Spq_lat(1) + h_in(3,2)*Spq_lat(2) + h_in(3,3)*Spq_lat(3)
          rnorm = dsqrt(Rpq_lat(1)**2.0_dp + Rpq_lat(2)**2.0_dp + Rpq_lat(3)**2.0_dp)
          if(rnorm.ge.Rc) cycle
          if(fn.eq.1) call mbdvdw_compute_TSR(p, q, Rpq_lat, Spq_lat, T_temp, dTdR_temp, dTdh_temp, dTdV_temp)
          if(fn.eq.2) call mbdvdw_compute_TLR(p, q, Rpq_lat, Spq_lat, T_temp, dTdR_temp, dTdh_temp, dTdV_temp)
          T = T + T_temp
          if(do_forces) dTdR = dTdR + dTdR_temp
          if(do_forces) dTdh = dTdh + dTdh_temp
          if(vdw_self_consistent) dTdV = dTdV + dTdV_temp
          end do ! loop over n3
        end do ! loop over n2
      end do !loop over n1
    end subroutine mbdvdw_TGG

  subroutine mbdvdw_compute_TLR_complex(p, q, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh, dTLRdV)
    implicit none
    ! IO Variables
    integer, intent(in) :: p, q
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    complex(dp), dimension(3, 3, nk), intent(out):: TLR
    complex(dp), dimension(3, 3, nat, 3, nk), intent(out) :: dTLRdr
    complex(dp), dimension(3, 3, 3, 3, nk), intent(out) :: dTLRdh
    complex(dp), dimension(3, 3, nat, nk), intent(out) :: dTLRdV

    real(dp) :: Rpq_norm, dot_prod, dRpq_normdh, dRpq_normdr, dkdh(3)
    real(dp) :: tlr_pre(3,3), dtlrdr_pre(3,3,nat,3), dtlrdh_pre(3,3,3,3), dtlrdv_pre(3,3,nat)
    real(dp), dimension(3) :: dRpqdR, dRpqdh, kvec
    complex(dp) :: expfactor, dexp_factor

    ! Loop variables
    integer :: i, s, ik

    Rpq_norm = 1.d0
    call mbdvdw_compute_TLR(p, q, Rpq, Spq_lat, tlr_pre, dtlrdr_pre, dtlrdh_pre, dtlrdv_pre)

    do ik = 1, nk, 1
      kvec = k_grid(ik, :)
      dot_prod = dot_product(kvec, Rpq)
      expfactor = exp(II*dot_prod)
      TLR(:, :, ik) = TLR_pre*expfactor

      if(vdw_self_consistent) then
        do s = 1, nat, 1
          dTLRdV(:, :, s, ik) = dtlrdv_pre(:,:,s)*exp(II*dot_prod)
          end do ! Loop over all possible components
        end if ! if to check if self consistent derivatives need to be computed

      if(do_forces) then
        ! dR
        do s = 1, nat, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
            dexp_factor = II*(dot_product(dRpqdR, kvec) + 0.0_dp)*expfactor
            dTLRdr(:, :, s, i, ik) = dtlrdr_pre(:,:,s,i)*expfactor &
                                    + TLR_pre*dexp_factor
            end do ! loop over cartesian components
          end do ! Loop over all atoms
        ! dH
        if(.not.mbd_vdw_isolated) then
          do s = 1, 3, 1
            do i = 1, 3, 1
              call mbdvdw_compute_dRdh_only(s, i, Rpq, Rpq_norm, Spq_lat, dRpqdh, dRpq_normdh)
              call mbdvdw_drec_vecdh(kvec, ruc, s, i, dkdh)
              dTLRdh(:, :, s, i, ik) = dtlrdh_pre(:,:,s,i)*expfactor &
                                      + II*(dot_product(dRpqdh, kvec) + dot_product(Rpq, dkdh))*expfactor*TLR_pre
              end do ! loop over all cell vector components
            end do ! Loop over all cell vectors
          end if
        end if ! if to check if forces and stresses need to be computed
      end do ! loop over ik

    end subroutine mbdvdw_compute_TLR_complex

  subroutine mbdvdw_tdip_ewald(R, Rnorm, tdip)
    implicit none
    real(dp), intent(in) :: R(3), Rnorm
    real(dp), intent(out) :: tdip(:,:)
    real(dp) :: b, c, R2, R3, R5
    integer :: i, j
    call mbdvdw_b(Rnorm, b)
    call mbdvdw_c(Rnorm, c)
    R2 = Rnorm*Rnorm
    R3 = R2*Rnorm
    R5 = R3*R2
    do i = 1, 3, 1
      do j = 1, 3, 1
        Tdip(i,j) = -C*R(i)*R(j)/R5
        end do
        Tdip(i,i) = Tdip(i,i) + b/R3
      end do

    end subroutine mbdvdw_tdip_ewald

  subroutine mbdvdw_TLR_ewald_realspace_real(p, q, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh, dTLRdV)
    implicit none
    ! IO Variables
    integer, intent(in) :: p, q
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    real(dp), dimension(3, 3), intent(out):: TLR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTLRdr
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTLRdh
    real(dp), dimension(3, 3, nat), intent(out) :: dTLRdV
    real(dp), dimension(3,3) :: Tdip, TdipEwald, dTdipdR, dRmatdR, dTdipdh, dRmatdh, dTdip_ewald_dh
    real(dp) :: Z, fermi_fn, R2, R3, R5, dRpq_normdR, dZdR, dSpqdR, dFermi_fndR, dRpq_normdh, dZdh, dSpqdh, dFermi_fndh, &
                dFermi_fndV, Rpq_norm, Spq, R_vdw_pq, dSpqdV, dtdip_ewald_dr(3,3), dzdv, Spq2, df_pre
    real(dp), dimension(3) :: dRpqdR, dRpqdh
    ! Loop variables
    integer :: i, s

    Rpq_norm = dsqrt(Rpq(1)**2.0_DP + Rpq(2)**2.0_DP + Rpq(3)**2.0_DP)

    R2 = Rpq_norm*Rpq_norm
    R3 = R2*Rpq_norm
    R5 = R3*R2
    call mbdvdw_compute_tdip(Rpq, R3, R5, Tdip)
    call mbdvdw_tdip_ewald(Rpq, Rpq_norm, TdipEwald)
    R_VdW_pq = R_MBD_VdW_sl(p) + R_MBD_VdW_sl(q)     ! Computes the damping radius: note the use of screened effective vdW radii
    Spq = beta*R_VdW_pq
    Spq2 = Spq*Spq
    Z = 6.0_DP*(Rpq_norm/Spq - 1.0_DP)
    fermi_fn = 1.0_DP
    df_pre = 0.0_dp
    if(Z.le.35.0_DP) fermi_fn = 1.0_DP/(1.0_DP+exp( -Z))
    if(Z.le.35.0_dp) df_pre = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP
    TLR = TdipEwald + (fermi_fn - 1.0_dp)*Tdip

    if(vdw_self_consistent) then
      do s = 1, nat, 1
        dSpqdV = beta*(dR_MBD_VdWdV_sl(p, s) + dR_MBD_VdWdV_sl(q, s))
        dZdV = 6.0_DP*( -Rpq_norm*dSpqdV/Spq**2.0_DP)
        dFermi_fndV = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdV
        dTLRdV(:, :, s) = Tdip*dFermi_fndV
        end do ! Loop over all possible components
      end if ! if to check if self consistent derivatives need to be computed

    if(do_forces) then
      ! dR
      do s = 1, nat, 1
        do i = 1, 3, 1
          call mbdvdw_compute_dRdR_ewald(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dTdip_ewald_dR)
          call mbdvdw_compute_dRdR(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dRmatdR, dTdipdR)
          dSpqdR = beta*(dR_MBD_VdWdR_sl(p, s, i) + dR_MBD_VdWdR_sl(q, s, i))
          dZdR = 6.0_DP*(dRpq_normdR/Spq - Rpq_norm*dSpqdR/Spq2)
          dFermi_fndR = df_pre*dZdR
          dTLRdr(:, :, s, i) = dTdip_ewald_dR + dFermi_fndR*Tdip + (fermi_fn - 1.0_dp)*dTdipdR
          end do ! loop over cartesian components
        end do ! Loop over all atoms
      ! dH
      if(.not.mbd_vdw_isolated) then
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdh_ewald(i, s, Rpq, Rpq_norm, Spq_lat, dadh(s, i), &
                                            dRpqdh, dRpq_normdh, dTdip_ewald_dh)
            call mbdvdw_compute_dRdh(s, i, Rpq, Rpq_norm, Spq_lat, &
                                            dRpqdh, dRpq_normdh, dRmatdh, dTdipdh)
            dSpqdh = beta*(dR_MBD_VdWdh_sl(p, s, i) + dR_MBD_VdWdh_sl(q, s, i))
            dZdh = 6.0_DP*(dRpq_normdh/Spq - Rpq_norm*dSpqdh/Spq2)
            dFermi_fndh = df_pre*dZdh
            dTLRdh(:, :, s, i) = dTdip_ewald_dh + dFermi_fndh*Tdip + (fermi_fn - 1.0_dp)*dTdipdh
            end do ! loop over all cell vector components
          end do ! Loop over all cell vectors
        end if
      end if ! if to check if forces and stresses need to be computed


    end subroutine mbdvdw_TLR_ewald_realspace_real

  subroutine mbdvdw_TLR_ewald_realspace(p, q, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh, dTLRdV)
    implicit none
    ! IO Variables
    integer, intent(in) :: p, q
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    complex(dp), dimension(3, 3, nk), intent(out):: TLR
    complex(dp), dimension(3, 3, nat, 3, nk), intent(out) :: dTLRdr
    complex(dp), dimension(3, 3, 3, 3, nk), intent(out) :: dTLRdh
    complex(dp), dimension(3, 3, nat, nk), intent(out) :: dTLRdV
    real(dp), dimension(3,3) :: Tdip, TdipEwald, dTdipdR, dTdip_ewald_dR, dRmatdR, dTdip_ewald_dh, dTdipdh, dRmatdh
    real(dp) :: Z, fermi_fn, R2, R3, R5, dRpq_normdR, dZdR, dSpqdR, dFermi_fndR, dRpq_normdh, dZdh, dSpqdh, dFermi_fndh, &
                dFermi_fndV, Rpq_norm, Spq, dot_prod, R_vdw_pq, dSpqdV, dZdV
    real(dp) :: tlr_pre(3,3), dtlrdr_pre(3,3,nat,3), dtlrdh_pre(3,3,3,3), dtlrdv_pre(3,3,nat)
    real(dp), dimension(3) :: dRpqdR, dRpqdh, kvec, dkdh
    complex(dp) :: expfactor, dexp_factor
    ! Loop variables
    integer :: i, s, ik

    Rpq_norm = dsqrt(Rpq(1)**2.0_DP + Rpq(2)**2.0_DP + Rpq(3)**2.0_DP)

    R2 = Rpq_norm*Rpq_norm
    R3 = R2*Rpq_norm
    R5 = R3*R2
    call mbdvdw_compute_tdip(Rpq, R3, R5, Tdip)
    call mbdvdw_tdip_ewald(Rpq, Rpq_norm, TdipEwald)
    R_VdW_pq = R_MBD_VdW_sl(p) + R_MBD_VdW_sl(q)     ! Computes the damping radius: note the use of screened effective vdW radii
    Spq = beta*R_VdW_pq
    Z = 6.0_DP*(Rpq_norm/Spq - 1.0_DP)
    fermi_fn = 1.0_DP
    if(Z.le.35.0_DP) fermi_fn = 1.0_DP/(1.0_DP+exp( -Z))
    TLR_pre = TdipEwald + (fermi_fn - 1.0_dp)*Tdip

    if(vdw_self_consistent) then
      do s = 1, nat, 1
        dSpqdV = beta*(dR_MBD_VdWdV_sl(p, s) + dR_MBD_VdWdV_sl(q, s))
        dZdV = 6.0_DP*( -Rpq_norm*dSpqdV/Spq**2.0_DP)
        dFermi_fndV = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdV
        dTLRdV_pre(:, :, s) = Tdip*dFermi_fndV
        end do ! Loop over all possible components
      end if ! if to check if self consistent derivatives need to be computed

    if(do_forces) then
      ! dR
      do s = 1, nat, 1
        do i = 1, 3, 1
          call mbdvdw_compute_dRdR_ewald(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dTdip_ewald_dR)
          call mbdvdw_compute_dRdR(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dRmatdR, dTdipdR)
          dSpqdR = beta*(dR_MBD_VdWdR_sl(p, s, i) + dR_MBD_VdWdR_sl(q, s, i))
          dZdR = 6.0_DP*(dRpq_normdR/Spq - Rpq_norm*dSpqdR/Spq**2.0_DP)
          dFermi_fndR = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdR
          dTLRdr_pre(:, :, s, i) = dTdip_ewald_dR + dFermi_fndR*Tdip + (fermi_fn - 1.0_dp)*dTdipdR
          end do ! loop over cartesian components
        end do ! Loop over all atoms
      ! dH
      if(.not.mbd_vdw_isolated) then
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdh_ewald(i, s, Rpq, Rpq_norm, Spq_lat, dadh(s, i), &
                                            dRpqdh, dRpq_normdh, dTdip_ewald_dh)
            call mbdvdw_compute_dRdh(s, i, Rpq, Rpq_norm, Spq_lat, &
                                            dRpqdh, dRpq_normdh, dRmatdh, dTdipdh)
            dSpqdh = beta*(dR_MBD_VdWdh_sl(p, s, i) + dR_MBD_VdWdh_sl(q, s, i))
            dZdh = 6.0_DP*(dRpq_normdh/Spq - Rpq_norm*dSpqdh/Spq**2.0_DP)
            dFermi_fndh = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdh
            dTLRdh_pre(:, :, s, i) = dTdip_ewald_dh + dFermi_fndh*Tdip + (fermi_fn - 1.0_dp)*dTdipdh
            end do ! loop over all cell vector components
          end do ! Loop over all cell vectors
        end if
      end if ! if to check if forces and stresses need to be computed

    do ik = 1, nk, 1
      kvec = k_grid(ik, :)
      dot_prod = dot_product(kvec, Rpq)
      expfactor = exp(II*dot_prod)
      TLR(:, :, ik) = tlr_pre*expfactor
      if(vdw_self_consistent) then
        do s = 1, nat, 1

          dTLRdV(:, :, s, ik) = dtlrdv_pre(:,:,s)*exp(II*dot_prod)
          end do ! Loop over all possible components
        end if ! if to check if self consistent derivatives need to be computed

      if(do_forces) then
        ! dR
        do s = 1, nat, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
            dexp_factor = II*(dot_product(dRpqdR, kvec) + 0.0_dp)*expfactor
            dTLRdr(:, :, s, i, ik) = dtlrdr_pre(:,:,s,i)*expfactor &
                                   + II*dot_product(dRpqdR,kvec)*TLR(:,:,ik)
            end do ! loop over cartesian components
          end do ! Loop over all atoms
        ! dH
        if(.not.mbd_vdw_isolated) then
          do s = 1, 3, 1
            do i = 1, 3, 1
              call mbdvdw_compute_dRdh_only(s, i, Rpq, Rpq_norm, Spq_lat, dRpqdh, dRpq_normdh)
              call mbdvdw_drec_vecdh(kvec, ruc, s, i, dkdh)
              dTLRdh(:, :, s, i, ik) = dtlrdh_pre(:,:,s,i)*expfactor &
                                       + II*(dot_product(dRpqdh,kvec) + dot_product(Rpq, dkdh))*TLR(:,:,ik)
              end do ! loop over all cell vector components
            end do ! Loop over all cell vectors
          end if
        end if ! if to check if forces and stresses need to be computed
      end do ! loop over ik
    end subroutine mbdvdw_TLR_ewald_realspace

  subroutine mbdvdw_TLR_ewald_recipspace(p, q, G, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh)
    integer, intent(in) :: p, q
    real(dp), intent(in) :: G(3), Rpq(3), Spq_lat(3)
    complex(dp), dimension(3, 3, nk), intent(out):: TLR
    complex(dp), dimension(3, 3, nat, 3, nk), intent(out) :: dTLRdr
    complex(dp), dimension(3, 3, 3, 3, nk), intent(out) :: dTLRdh

    real(dp) :: a2, k(3), k_norm, k2, kmat(3,3), Rpq_norm, dRpqdR(3), dRpq_normdR, dRpqdh(3), dRpq_normdh, &
                a, dk_norm, dkmat(3,3), dGdh(3), dkdh(3)
    complex(dp) :: exp_factor, exp2_factor, dexp_factor, pre_factor, dexp2_factor
    integer :: ik, s, i

    a = ewald_cutoff
    exp_factor = exp(-ii*dot_product(G, Rpq))
    a2 = 4.0_dp*a*a
    Rpq_norm = 1.d0

    do ik = 1, nk, 1
      k = G + k_grid(ik, :)
      k2 = k(1)**2.0_dp + k(2)**2.0_dp + k(3)**2.0_dp
      k_norm = dsqrt(k2)
      if (k_norm.le.1e-15) then
          TLR(:, :, ik) = 0.0_DP
          if (do_forces) then
              dTLRdr(:, :, :, :, ik) = 0.0_DP
              if (.not. mbd_vdw_isolated) dTLRdh(:, :, :, :, ik) = 0.0_DP
          end if
          cycle
      end if
      exp2_factor = exp(-k2/a2)
      call mbdvdw_kmat(k, k_norm, kmat)
      pre_factor = exp_factor*exp2_factor
      TLR(:,:,ik) = pre_factor*kmat

      if(do_forces) then
        ! dR
        do s = 1, nat, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
            dexp_factor = -II*(dot_product(dRpqdR, G) + 0.0_dp)*exp_factor
            dTLRdr(:, :, s, i, ik) = kmat*exp2_factor*dexp_factor
            end do ! loop over cartesian components
          end do ! Loop over all atoms
        ! dH
        if(.not.mbd_vdw_isolated) then
          do s = 1, 3, 1
            do i = 1, 3, 1
              call mbdvdw_drec_vecdh(G, ruc, s, i, dGdh)
              call mbdvdw_compute_dRdh_only(s, i, Rpq, Rpq_norm, Spq_lat, dRpqdh, dRpq_normdh)

              call mbdvdw_drec_vec_norm_dh(k, k_norm, ruc, s, i, dk_norm)
              call mbdvdw_drec_vecdh(k, ruc, s, i, dkdh)
              call mbdvdw_dkmatdh(k, dkdh, k_norm, dk_norm, dkmat)
              dexp2_factor = k_norm*(k_norm*dadh(s,i) - dk_norm*a)/(2.0_dp*a*a*a)*exp2_factor
              dexp_factor = -II*(dot_product(dGdh, Rpq) + dot_product(G, dRpqdh))*exp_factor
              dTLRdh(:, :, s, i, ik) = dkmat*exp2_factor*exp_factor &
                                      +kmat*dexp2_factor*exp_factor &
                                      +kmat*exp2_factor*dexp_factor

              end do ! loop over all cell vector components
            end do ! Loop over all cell vectors
          end if
        end if ! if to check if forces and stresses need to be computed
      end do
    end subroutine mbdvdw_TLR_ewald_recipspace

  subroutine mbdvdw_TLR_ewald_recipspace_real(p, q, G, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh)
    integer, intent(in) :: p, q
    real(dp), intent(in) :: G(3), Rpq(3), Spq_lat(3)
    real(dp), dimension(3, 3), intent(out):: TLR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTLRdr
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTLRdh

    real(dp) :: a2, k(3), k_norm, k2, kmat(3,3), Rpq_norm, dRpqdR(3), dRpq_normdR, dRpqdh(3), dRpq_normdh, &
                a, dk_norm, dkmat(3,3), dGdh(3), a3, dot_prod
    real(dp) :: exp_factor, exp2_factor, dexp_factor, pre_factor, dexp2_factor, dexp
    integer :: s, i

    a = ewald_cutoff
    dot_prod = dot_product(G, Rpq)
    exp_factor = cos(-dot_prod)! exp(-ii*dot_product(G, Rpq))
    dexp = sin(-dot_prod)
    a2 = 4.0_dp*a*a
    a3 = 2.0_dp*a*a*a
    Rpq_norm = 1.d0

    k = G
    k2 = k(1)**2.0_dp + k(2)**2.0_dp + k(3)**2.0_dp
    k_norm = sqrt(k2)
    exp2_factor = exp(-k2/a2)
    call mbdvdw_kmat(k, k_norm, kmat)
    pre_factor = exp_factor*exp2_factor
    TLR = pre_factor*kmat

    if(do_forces) then
      ! dR
      do s = 1, nat, 1
        do i = 1, 3, 1
          call mbdvdw_compute_dRdR_only(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR)
          dexp_factor = (dot_product(dRpqdR, G))*dexp
          dTLRdr(:, :, s, i) = kmat*exp2_factor*dexp_factor
          end do ! loop over cartesian components
        end do ! Loop over all atoms
      ! dH
      if(.not.mbd_vdw_isolated) then
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdh_only(s, i, Rpq, Rpq_norm, Spq_lat, dRpqdh, dRpq_normdh)

            call mbdvdw_drec_vecdh(G, ruc, s, i, dGdh)
            call mbdvdw_drec_vec_norm_dh(k, k_norm, ruc, s, i, dk_norm)
            call mbdvdw_dkmatdh(G, dGdh, k_norm, dk_norm, dkmat)

            dexp2_factor = k_norm*(k_norm*dadh(s,i) - dk_norm*a)/(a3)*exp2_factor
            dexp_factor = (dot_product(dGdh, Rpq) + dot_product(G, dRpqdh))*dexp
            dTLRdh(:, :, s, i) = dkmat*exp2_factor*exp_factor &
                                    +kmat*dexp2_factor*exp_factor &
                                    +kmat*exp2_factor*dexp_factor

            end do ! loop over all cell vector components
          end do ! Loop over all cell vectors
        end if
      end if ! if to check if forces and stresses need to be computed
    end subroutine mbdvdw_TLR_ewald_recipspace_real

  subroutine mbdvdw_TGG_complex(fn, p, q, n, h_in, ainv_in, tau_in, T, dTdR, dTdh, dTdV)
    integer, intent(in) :: p, q, n, fn
    real(dp), intent(in) :: h_in(3,3), tau_in(3,n), ainv_in(3,3)
    complex(dp), dimension(3, 3, nk), intent(out):: T
    complex(dp), dimension(3, 3, nat, 3, nk), intent(out) :: dTdR
    complex(dp), dimension(3, 3, nat, nk), intent(out) :: dTdV
    complex(dp), dimension(3, 3, 3, 3, nk), intent(out) :: dTdh
    real(dp), dimension(3) :: Rpq, Rpq_lat
    real(dp), dimension(3) :: Spq_lat, Spq
    complex(dp), dimension(3, 3, nk) :: T_temp
    complex(dp), dimension(3, 3, nat, 3, nk) :: dTdR_temp
    complex(dp), dimension(3, 3, nat, nk)  :: dTdV_temp
    complex(dp), dimension(3, 3, 3, 3, nk) :: dTdh_temp
    real(dp) :: max_value, rnorm
    integer :: n_period
    logical :: converged
    integer :: n1, n2, n3, i_idx, j_idx, k_idx
    integer :: conv  ! auxiliary convergence counter

    T = cmplx(0.0_dp, 0.0_dp, kind=DP)
    dTdh = cmplx(0.0_dp, 0.0_dp, kind=DP)
    dTdR = cmplx(0.0_dp, 0.0_dp, kind=DP)
    dTdV = cmplx(0.0_dp, 0.0_dp, kind=DP)

    Rpq = tau_in(:, p) - tau_in(:, q)

    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)

    ! Spq(1) = Spq(1) - anint(Spq(1))
    ! Spq(2) = Spq(2) - anint(Spq(2))
    ! Spq(3) = Spq(3) - anint(Spq(3))

    ! Computes the vector distance between two atoms
    n_period = 0
    max_value = 0.0_DP
    converged = .false.
    do while(.not.converged)
      max_value = 0.0_DP
      conv=0
      do n1 = -n_period, n_period
        do n2 = -n_period, n_period
          do n3 = -n_period, n_period
            if((abs(n1).eq.n_period).or.(abs(n2).eq.n_period).or.(abs(n3).eq.n_period)) then
              i_idx = n1
              j_idx = n2
              k_idx = n3
              if(mbd_vdw_vacuum(1)) then
                if(n1.ne.0) cycle
                i_idx = 0
                end if
              if(mbd_vdw_vacuum(2)) then
                if(n2.ne.0) cycle
                j_idx = 0
                end if
              if(mbd_vdw_vacuum(3)) then
                if(n3.ne.0) cycle
                k_idx = 0
                end if

              Spq_lat(1) = Spq(1) + DBLE(i_idx)
              Spq_lat(2) = Spq(2) + DBLE(j_idx)
              Spq_lat(3) = Spq(3) + DBLE(k_idx)

              Rpq_lat(1) = h_in(1,1)*Spq_lat(1) + h_in(1,2)*Spq_lat(2) + h_in(1,3)*Spq_lat(3)
              Rpq_lat(2) = h_in(2,1)*Spq_lat(1) + h_in(2,2)*Spq_lat(2) + h_in(2,3)*Spq_lat(3)
              Rpq_lat(3) = h_in(3,1)*Spq_lat(1) + h_in(3,2)*Spq_lat(2) + h_in(3,3)*Spq_lat(3)

              rnorm = sqrt(Rpq_lat(1)**2.0_dp + Rpq_lat(2)**2.0_dp + Rpq_lat(3)**2.0_dp)
              if(rnorm.le.1e-12) then
                conv = conv + 1
                cycle
                end if
              if(fn.eq.2) call mbdvdw_compute_TLR_complex(p, q, Rpq_lat, Spq_lat, T_temp, dTdR_temp, dTdh_temp, dTdV_temp)
              if(abs(maxval(abs(T_temp))).ge.mbd_vdw_econv_thr) conv = conv + 1
              T = T + T_temp
              if(vdw_self_consistent) dTdV = dTdV + dTdV_temp
              if(do_forces) dTdR = dTdR + dTdR_temp
              if(do_forces) dTdh = dTdh + dTdh_temp
              end if
            end do
          end do
        end do
      if((conv.eq.0).or.mbd_vdw_isolated) converged = .true.
      n_period = n_period + 1
      end do

    end subroutine mbdvdw_TGG_complex

  subroutine mbdvdw_ewaldsum_real(p, q, nelem, h_in, ainv_in, tau_in, TLR, dTLRdR, dTLRdh, dTLRdV)
    integer, intent(in) :: p, q, nelem
    real(dp), intent(in) :: h_in(3,3), ainv_in(3,3), tau_in(3, nelem)
    real(dp), intent(out) :: TLR(3,3), dTLRdR(3,3,nat,3), dTLRdh(3,3,3,3), dTLRdV(3,3,nat)
    real(dp) :: TLR_tmp(3,3), dTLRdR_tmp(3,3,nat,3), dTLRdh_tmp(3,3,3,3), dTLRdV_tmp(3,3,nat)
    real(dp) :: Rpq(3), Spq(3), G(3), Spq_lat(3), Rpq_lat(3), nu2, surface_term, dSurface_term
    real(dp) :: a3, knorm, kvec(3), rnorm, se, dSe
    integer :: sc(3), n1, n2, n3, n(3), s, i, j

    Rpq = tau_in(:, p) - tau_in(:, q)
    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)

    Rpq(1) = h_in(1,1)*Spq(1) + h_in(1,2)*Spq(2) + h_in(1,3)*Spq(3)
    Rpq(2) = h_in(2,1)*Spq(1) + h_in(2,2)*Spq(2) + h_in(2,3)*Spq(3)
    Rpq(3) = h_in(3,1)*Spq(1) + h_in(3,2)*Spq(2) + h_in(3,3)*Spq(3)

    TLR = 0.0_dp
    dTLRdh = 0.0_dp
    dTLRdR = 0.0_dp
    dTLRdV = 0.0_dp

    ! Reciprocal space bits
    sc = mbdvdw_circumscribe(ruc, Gc)

    do n1 = -sc(1), sc(1), 1
      do n2 = -sc(2), sc(2), 1
        do n3 = -sc(3), sc(3), 1
          if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0)) cycle
          n = (/ n1, n2, n3 /)
          G = matmul(dble(n), ruc)
          if (sum(G**2) > Gc**2) cycle
          call mbdvdw_TLR_ewald_recipspace_real(p, q, G, Rpq, Spq, TLR_tmp, dTLRdR_tmp, dTLRdh_tmp)
          TLR = TLR + TLR_tmp
          if(do_forces) dTLRdR = dTLRdR + dTLRdR_tmp
          if(do_forces) dTLRdh = dTLRdh + dTLRdh_tmp
          end do ! n3
        end do ! n2
      end do ! n1
    nu2 = nu*nu
    if(do_forces) then
      dTLRdR = 4.0_dp*pi*dTLRdR/nu
      do s = 1, 3, 1
        do i = 1, 3, 1
          dTLRdh(:,:,s,i) = 4.0_dp*pi*(dTLRdh(:,:,s,i)/nu - TLR(:,:)*dnudh(s,i)/(nu2))
          end do
        end do
      end if
    TLR = 4.0_dp*pi*TLR/nu

    ! Self Energy. Should have \delta_{pq} \delta_{ij}
    if(p.eq.q) then
      a3 = ewald_cutoff*ewald_cutoff*ewald_cutoff
      se = -4.0_dp*a3/(3.0_dp*SQRTPI)
      do i = 1, 3, 1
        TLR(i, i) = TLR(i, i) + se
        if(do_forces) then
          do s = 1, 3, 1
            do j = 1, 3, 1
              dSe = -4.0_dp*ewald_cutoff*ewald_cutoff*dadh(s, j)/(SQRTPI)
              dTLRdh(i, i, s, j) = dTLRdh(i, i, s, j) + dSe
              end do ! loop over unit cell vector components
            end do ! loop over unit cell vectors
          end if ! If statement to check it we need to do the forces
        end do ! Loop to apply the \delta_{ij}
      end if ! If statement to apply the \delta_{pq}

    ! Surface term. Should have \delta(k=0)n \delta_{ij}
    kvec = k_grid(1, :)
    knorm = dsqrt(kvec(1)**2.0_dp + kvec(2)**2.0_dp + kvec(3)**2.0_dp)
    if(knorm.le.1e-12) then
      surface_term = 4.0_dp*pi/(3.0_dp*nu)
      do i = 1, 3, 1
        tlr(i, i) = tlr(i, i) + surface_term
        if(do_forces) then
          do s = 1, 3, 1
            do j = 1, 3, 1
              dSurface_term = -4.0_dp*pi/(3.0_dp*nu*nu)*dnudh(s, j)
              dTlrdh(i, i, s, j) = dTlrdh(i, i, s, j) + dSurface_term
              end do ! loop over all unit cell vector components
            end do ! loop over all unit cell vectors
          end if ! if statement to check if we should do forces
        end do ! loop to apply the \delta_{ij}
      end if ! if to see if k=0


    ! Realspace bits
    Rpq = tau_in(:, p) - tau_in(:, q)
    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)
    sc = mbdvdw_circumscribe(h_in, Rc)
    do n1 = -sc(1), sc(1), 1
      do n2 = -sc(2), sc(2), 1
        do n3 = -sc(3), sc(3), 1
          if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0)) then
            if(p.eq.q) cycle
            end if
          Spq_lat(1) = Spq(1) + DBLE(n1)
          Spq_lat(2) = Spq(2) + DBLE(n2)
          Spq_lat(3) = Spq(3) + DBLE(n3)

          Rpq_lat(1) = h_in(1,1)*Spq_lat(1) + h_in(1,2)*Spq_lat(2) + h_in(1,3)*Spq_lat(3)
          Rpq_lat(2) = h_in(2,1)*Spq_lat(1) + h_in(2,2)*Spq_lat(2) + h_in(2,3)*Spq_lat(3)
          Rpq_lat(3) = h_in(3,1)*Spq_lat(1) + h_in(3,2)*Spq_lat(2) + h_in(3,3)*Spq_lat(3)
          rnorm = dsqrt(Rpq_lat(1)**2.0_dp + Rpq_lat(2)**2.0_dp + Rpq_lat(3)**2.0_dp)
          if(rnorm.ge.Rc) cycle
          call mbdvdw_TLR_ewald_realspace_real(p, q, Rpq_lat, Spq_lat, TLR_tmp, dTLRdR_tmp, dTLRdh_tmp, dTLRdV_tmp)
          TLR = TLR + TLR_tmp
          if(do_forces) dTLRdR = dTLRdR + dTLRdR_tmp
          if(do_forces) dTLRdh = dTLRdh + dTLRdh_tmp
          if(vdw_self_consistent) dTLRdV = dTLRdV + dTLRdV_tmp
          end do ! loop over n3
        end do ! loop over n2
      end do !loop over n1

    end subroutine mbdvdw_ewaldsum_real

  subroutine mbdvdw_ewaldsum(p, q, nelem, h_in, ainv_in, tau_in, TLR, dTLRdR, dTLRdh, dTLRdV)
    integer, intent(in) :: p, q, nelem
    real(dp), intent(in) :: h_in(3,3), ainv_in(3,3), tau_in(3, nelem)
    complex(dp), intent(out) :: TLR(3,3,nk), dTLRdR(3,3,nat,3,nk), dTLRdh(3,3,3,3,nk), dTLRdV(3,3,nat,nk)
    complex(dp) :: TLR_tmp(3,3,nk), dTLRdR_tmp(3,3,nat,3,nk), dTLRdh_tmp(3,3,3,3,nk), dTLRdV_tmp(3,3,nat,nk)
    real(dp) :: Rpq(3), Spq(3), G(3), Spq_lat(3), Rpq_lat(3), nu2, surface_term, dSurface_term
    real(dp) :: a3, knorm, kvec(3), rnorm, se, dSe
    integer :: sc(3), n1, n2, n3, s, i, ik, j, n(3)

    Rpq = tau_in(:, p) - tau_in(:, q)
    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)

    Rpq(1) = h_in(1,1)*Spq(1) + h_in(1,2)*Spq(2) + h_in(1,3)*Spq(3)
    Rpq(2) = h_in(2,1)*Spq(1) + h_in(2,2)*Spq(2) + h_in(2,3)*Spq(3)
    Rpq(3) = h_in(3,1)*Spq(1) + h_in(3,2)*Spq(2) + h_in(3,3)*Spq(3)

    TLR = 0.0_dp
    dTLRdh = 0.0_dp
    dTLRdR = 0.0_dp
    dTLRdV = 0.0_dp

    dnudh = 0.0_dp
    ruc = compute_recucell(h_in)
    call mbdvdw_cellvol(h_in, nu)
    do s = 1, 3, 1
      do i = 1, 3, 1
        call mbdvdw_dcellvoldh(ainv_in, s, i, nu, dnudh(s, i))
        end do
      end do
    call mbdvdw_ewald_params(nu)
    call mbdvdw_dewald_params_dh(ainv_in, nu)

    ! Reciprocal space bits
    sc = mbdvdw_circumscribe(ruc, Gc)
    do n1 = -sc(1), sc(1), 1
      do n2 = -sc(2), sc(2), 1
        do n3 = -sc(3), sc(3), 1
          n = (/ n1, n2, n3 /)
          G = matmul(dble(n), ruc)
          if (sum(G**2) > Gc**2) cycle
          call mbdvdw_TLR_ewald_recipspace(p, q, G, Rpq, Spq, TLR_tmp, dTLRdR_tmp, dTLRdh_tmp)
          TLR = TLR + TLR_tmp
          if(do_forces) dTLRdR = dTLRdR + dTLRdR_tmp
          if(do_forces) dTLRdh = dTLRdh + dTLRdh_tmp
          end do ! n3
        end do ! n2
      end do ! n1
    nu2 = nu*nu
    if(do_forces) then
      dTLRdR = 4.0_dp*pi*dTLRdR/nu
      do ik = 1, nk, 1
        do s = 1, 3, 1
          do i = 1, 3, 1
            dTLRdh(:,:,s,i,ik) = 4.0_dp*pi*(dTLRdh(:,:,s,i,ik)/nu - TLR(:,:,ik)*dnudh(s,i)/(nu2))
            end do
          end do
        end do
      end if
    TLR = 4.0_dp*pi*TLR/nu

    ! Self Energy. Should have \delta_{pq} \delta_{ij}
    if(p.eq.q) then
      a3 = ewald_cutoff*ewald_cutoff*ewald_cutoff
      se = -4.0_dp*a3/(3.0_dp*SQRTPI)
      do ik = 1, nk, 1
        do i = 1, 3, 1
          TLR(i, i, ik) = TLR(i, i, ik) + se
          if(do_forces) then
            do s = 1, 3, 1
              do j = 1, 3, 1
                dSe = -4.0_dp*ewald_cutoff*ewald_cutoff*dadh(s, j)/(SQRTPI)
                dTLRdh(i, i, s, j, ik) = dTLRdh(i, i, s, j, ik) + dSe
                end do ! loop over unit cell vector components
              end do ! loop over unit cell vectors
            end if ! If statement to check it we need to do the forces
          end do ! Loop to apply the \delta_{ij}
        end do ! loop over ik
      end if ! If statement to apply the \delta_{pq}

    ! Surface term. Should have \delta(k=0)n \delta_{ij}
    do ik = 1, nk, 1
      kvec = k_grid(ik, :)
      knorm = dsqrt(kvec(1)**2.0_dp + kvec(2)**2.0_dp + kvec(3)**2.0_dp)
      if(knorm.le.1e-12) then
        surface_term = 4.0_dp*pi/(3.0_dp*nu)
        do i = 1, 3, 1
          tlr(i, i, ik) = tlr(i, i, ik) + surface_term
          if(do_forces) then
            do s = 1, 3, 1
              do j = 1, 3, 1
                dSurface_term = -4.0_dp*pi/(3.0_dp*nu*nu)*dnudh(s, j)
                dTlrdh(i, i, s, j, ik) = dTlrdh(i, i, s, j, ik) + dSurface_term
                end do ! loop over all unit cell vector components
              end do ! loop over all unit cell vectors
            end if ! if statement to check if we should do forces
          end do ! loop to apply the \delta_{ij}
        end if ! if to see if k=0
      end do ! loop over ik


    ! Realspace bits
    Rpq = tau_in(:, p) - tau_in(:, q)
    Spq(1) = ainv_in(1,1)*rpq(1) + ainv_in(1,2)*rpq(2) + ainv_in(1,3)*rpq(3)
    Spq(2) = ainv_in(2,1)*rpq(1) + ainv_in(2,2)*rpq(2) + ainv_in(2,3)*rpq(3)
    Spq(3) = ainv_in(3,1)*rpq(1) + ainv_in(3,2)*rpq(2) + ainv_in(3,3)*rpq(3)
    sc = mbdvdw_circumscribe(h_in, Rc)
    do n1 = -sc(1), sc(1), 1
      do n2 = -sc(2), sc(2), 1
        do n3 = -sc(3), sc(3), 1
          if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0)) then
            if(p.eq.q) cycle
            end if
          Spq_lat(1) = Spq(1) + DBLE(n1)
          Spq_lat(2) = Spq(2) + DBLE(n2)
          Spq_lat(3) = Spq(3) + DBLE(n3)

          Rpq_lat(1) = h_in(1,1)*Spq_lat(1) + h_in(1,2)*Spq_lat(2) + h_in(1,3)*Spq_lat(3)
          Rpq_lat(2) = h_in(2,1)*Spq_lat(1) + h_in(2,2)*Spq_lat(2) + h_in(2,3)*Spq_lat(3)
          Rpq_lat(3) = h_in(3,1)*Spq_lat(1) + h_in(3,2)*Spq_lat(2) + h_in(3,3)*Spq_lat(3)
          rnorm = dsqrt(Rpq_lat(1)**2.0_dp + Rpq_lat(2)**2.0_dp + Rpq_lat(3)**2.0_dp)
          if(rnorm.ge.Rc) cycle
          n = (/ n1, n2, n3 /)
          call mbdvdw_TLR_ewald_realspace(p, q, Rpq_lat, Spq_lat, TLR_tmp, dTLRdR_tmp, dTLRdh_tmp, dTLRdV_tmp)
          TLR = TLR + TLR_tmp
          if(do_forces) dTLRdR = dTLRdR + dTLRdR_tmp
          if(do_forces) dTLRdh = dTLRdh + dTLRdh_tmp
          if(vdw_self_consistent) dTLRdV = dTLRdV + dTLRdV_tmp
          end do ! loop over n3
        end do ! loop over n2
      end do !loop over n1
    end subroutine mbdvdw_ewaldsum

  subroutine mbdvdw_compute_TLR(p, q, Rpq, Spq_lat, TLR, dTLRdR, dTLRdh, dTLRdV)
    implicit none

    ! IO Variables
    integer, intent(in) :: p, q
    real(dp), dimension(3), intent(in) :: Rpq, Spq_lat
    real(dp), dimension(3, 3), intent(out):: TLR
    real(dp), dimension(3, 3, nat, 3), intent(out) :: dTLRdr
    real(dp), dimension(3, 3, 3, 3), intent(out) :: dTLRdh
    real(dp), dimension(3, 3, nat), intent(out) :: dTLRdV

    real(dp) :: Rpq_norm, Spq
    real(dp) :: R_vdw_pq
    real(dp), dimension(3,3) :: Tdip, dTdipdR, dRmatdR, dTdipdh, dRmatdh
    real(dp) :: Z, fermi_fn, R2, R3, R5, dRpq_normdR, dZdR, dSpqdR, dFermi_fndR, dRpq_normdh, dZdh, dSpqdh, dFermi_fndh, &
                dZdV, dSpqdV, dFermi_fndV
    real(dp), dimension(3) :: dRpqdR, dRpqdh

    ! Loop variables
    integer :: i, s

    ! Computes the cartesian distance from the vector distance
    Rpq_norm = dsqrt(Rpq(1)**2.0_DP + Rpq(2)**2.0_DP + Rpq(3)**2.0_DP)
    if(Rpq_norm.le.1e-12) print*, p, q, Rpq
    R2 = Rpq_norm*Rpq_norm
    R3 = R2*Rpq_norm
    R5 = R3*R2
    call mbdvdw_compute_tdip(Rpq, R3, R5, Tdip)
    R_VdW_pq = R_MBD_VdW_sl(p) + R_MBD_VdW_sl(q)     ! Computes the damping radius: note the use of screened effective vdW radii
    Spq = beta*R_VdW_pq
    Z = 6.0_DP*(Rpq_norm/Spq - 1.0_DP)
    fermi_fn = 1.0_DP
    if(Z.le.35.0_DP) fermi_fn = 1.0_DP/(1.0_DP+exp( -Z))
    TLR = fermi_fn*Tdip

    if(vdw_self_consistent) then
      do s = 1, nat, 1
        dSpqdV = beta*(dR_MBD_VdWdV_sl(p, s) + dR_MBD_VdWdV_sl(q, s))
        dZdV = 6.0_DP*( -Rpq_norm*dSpqdV/Spq**2.0_DP)
        dFermi_fndV = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdV
        dTLRdV(:, :, s) = Tdip*dFermi_fndV
        end do ! Loop over all possible components
      end if ! if to check if self consistent derivatives need to be computed

    if(do_forces) then
      ! dR
      do s = 1, nat, 1
        do i = 1, 3, 1
          call mbdvdw_compute_dRdR(s, i, p, q, Rpq, Rpq_norm, dRpqdR, dRpq_normdR, dRmatdR, dTdipdR)
          dSpqdR = beta*(dR_MBD_VdWdR_sl(p, s, i) + dR_MBD_VdWdR_sl(q, s, i))
          dZdR = 6.0_DP*(dRpq_normdR/Spq - Rpq_norm*dSpqdR/Spq**2.0_DP)
          dFermi_fndR = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdR
          dTLRdr(:, :, s, i) = Tdip*dFermi_fndR + fermi_fn*dTdipdR
          end do ! loop over cartesian components
        end do ! Loop over all atoms
      ! dH
      if(.not.mbd_vdw_isolated) then
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_compute_dRdh(s, i, Rpq, Rpq_norm, Spq_lat, dRpqdh, dRpq_normdh, dRmatdh, dTdipdh)
            dSpqdh = beta*(dR_MBD_VdWdh_sl(p, s, i) + dR_MBD_VdWdh_sl(q, s, i))
            dZdh = 6.0_DP*(dRpq_normdh/Spq - Rpq_norm*dSpqdh/Spq**2.0_DP)
            dFermi_fndh = exp(-Z)/(1.0_DP + exp(-Z))**2.0_DP*dZdh
            dTLRdh(:, :, s, i) = Tdip*dFermi_fndh + fermi_fn*dTdipdh
            end do ! loop over all cell vector components
          end do ! Loop over all cell vectors
        end if
      end if ! if to check if forces and stresses need to be computed
    end subroutine mbdvdw_compute_TLR

  subroutine mbdvdw_send_recv_forces(my_dT, collected_dT, max_proc, npairs, splits, foff, off)
    implicit none
    real(dp), intent(in) :: my_dT(:,:,:,:)
    real(dp), intent(inout) :: collected_dT(:,:,:,:)
    integer, intent(in) :: npairs(:), splits(:), off(:), foff(:), max_proc
    real(dp), allocatable :: temp_dt(:,:,:,:)
    integer :: ind_p, ind_f, i, j, ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    do ind_p = 1, npairs(me_image+1), 1
      do ind_f = 1, splits(me_image+1), 1
        collected_dT(off(me_image+1)+ind_p,:,:,ind_f) = &
                      my_dT(ind_p,:,:,foff(me_image+1)+ind_f)
        end do
      end do

    do i = 0, max_proc-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs

        if((j.eq.me_image).and.(i.ne.me_image)) then
          if(.not.allocated(temp_dT))    allocate(temp_dT(npairs(j+1), 3, 3, splits(i+1)))
          do ind_f = 1, splits(i+1), 1
            temp_dT(:, :, :, ind_f) = my_dT(:, :, :, foff(i+1) + ind_f)
            end do
          call mpi_send(temp_dT,&                               ! Send buffer
                  npairs(j+1)*3*3*splits(i+1),& ! Count of communicated data
                  mpi_double_precision,&                   ! Data type
                  i,&                                      ! Destination
                  j,&                                      ! Tag (who it was sent from)
                  intra_image_comm,&                       ! Communicator to send with
                  ierr)                                    ! Error flag
          if(allocated(temp_dT))  deallocate(temp_dT)

        else

          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_dT))    allocate(temp_dT(npairs(j+1), 3, 3, splits(i+1)))
            call mpi_recv(temp_dT,&                              ! Recv buffer
                    npairs(j+1)*3*3*splits(i+1),& ! Count of communicated data
                    mpi_double_precision,&                  ! Data type
                    j,&                                     ! Source (who it was sent from)
                    j,&                                     ! Tag (also who it was sent from)
                    intra_image_comm,&                      ! Communicator to send with
                    rstatus,&                               ! Status flag
                    ierr)                                   ! Error flag
            do ind_f = 1, npairs(j+1), 1
              collected_dT(ind_f+off(j+1),:,:,:) = temp_dT(ind_f,:,:,:)
              end do
            if(allocated(temp_dT))  deallocate(temp_dT)
            end if
          end if
        end do
      end do
    if(allocated(temp_dT)) deallocate(temp_dT)

    end subroutine mbdvdw_send_recv_forces

    subroutine mbdvdw_send_recv_forces_complex(my_dT, collected_dT, max_proc, npairs, splits, foff, off)
      implicit none
      complex(dp), intent(in) :: my_dT(:,:,:,:)
      complex(dp), intent(inout) :: collected_dT(:,:,:,:)
      integer, intent(in) :: npairs(:), splits(:), off(:), foff(:), max_proc
      complex(dp), allocatable :: temp_dt(:,:,:,:)
      integer :: ind_p, ind_f, i, j, ierr
      integer, dimension(MPI_STATUS_SIZE) :: rstatus

      do ind_p = 1, npairs(me_image+1), 1
        do ind_f = 1, splits(me_image+1), 1
          collected_dT(off(me_image+1)+ind_p,:,:,ind_f) = &
                        my_dT(ind_p,:,:,foff(me_image+1)+ind_f)
          end do
        end do

      do i = 0, max_proc-1, 1                      ! Outer loop is the loop for recvs
        do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs

          if((j.eq.me_image).and.(i.ne.me_image)) then
            if(.not.allocated(temp_dT))    allocate(temp_dT(npairs(j+1), 3, 3, splits(i+1)))
            do ind_f = 1, splits(i+1), 1
              temp_dT(:, :, :, ind_f) = my_dT(:, :, :, foff(i+1) + ind_f)
              end do
            call mpi_send(temp_dT,&                               ! Send buffer
                    npairs(j+1)*3*3*splits(i+1),& ! Count of communicated data
                    mpi_complex16,&                          ! Data type
                    i,&                                      ! Destination
                    j,&                                      ! Tag (who it was sent from)
                    intra_image_comm,&                       ! Communicator to send with
                    ierr)                                    ! Error flag
            if(allocated(temp_dT))  deallocate(temp_dT)

          else

            if((i.eq.me_image).and.(j.ne.me_image)) then
              if(.not.allocated(temp_dT))    allocate(temp_dT(npairs(j+1), 3, 3, splits(i+1)))
              call mpi_recv(temp_dT,&                              ! Recv buffer
                      npairs(j+1)*3*3*splits(i+1),& ! Count of communicated data
                      mpi_complex16,&                  ! Data type
                      j,&                                     ! Source (who it was sent from)
                      j,&                                     ! Tag (also who it was sent from)
                      intra_image_comm,&                      ! Communicator to send with
                      rstatus,&                               ! Status flag
                      ierr)                                   ! Error flag
              do ind_f = 1, npairs(j+1), 1
                collected_dT(ind_f+off(j+1),:,:,:) = temp_dT(ind_f,:,:,:)
                end do
              if(allocated(temp_dT))  deallocate(temp_dT)
              end if
            end if
          end do
        end do
      if(allocated(temp_dT)) deallocate(temp_dT)

    end subroutine mbdvdw_send_recv_forces_complex

  subroutine mbdvdw_SCS()
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This method solves the self consistent screen equation to solve for \bar{A}
    ! It does by solving the equation
    ! \bar{A}(i\omega) = \left(  A^{-1} + TSR \right)^{-1}
    ! Because A is diagonal, the inverse is trivial, and will be dealt with in the inner loops.
    ! TSR will be generated for each pair of atoms
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    ! Local variables
    real(dp), dimension(3, 3) :: TSR
    real(dp), dimension(3*nat, 3*nat) :: temp, temp1
    real(dp), dimension(:,:,:,:), allocatable:: dTSRdR
    real(dp), dimension(:,:,:), allocatable :: dTSRdv, TsrdV, TsrdR, Tsrdh
    real(dp), dimension(3,3,3,3) :: dTSRdh
    real(dp), dimension(:,:,:),   allocatable :: my_tsr, collected_tsr, temp_tsr
    real(dp), dimension(:,:,:,:), allocatable :: my_dtsrdr, collected_dtsrdr
    real(dp), dimension(:,:,:,:), allocatable :: my_dtsrdh, collected_dtsrdh
    real(dp), dimension(:,:,:,:), allocatable :: my_dtsrdV, collected_dtsrdV

    ! Looping varlables
    integer :: p, q, i_idx, j_idx, i_index, j_index, i, s, i_f, i_index_t, j_index_t

    ! lapack variables
    integer :: errorflag
    integer, dimension(3*nat) :: IPIV
    real(dp), dimension(3*nat) :: WORK

    !MPI Vars
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    integer, dimension(nproc_image) :: offsets, force_offsets, sc_offsets, force_offsets_h

    ! Looping varlables
    integer :: counter, counter_amat, cnt_v, cnt_h, cnt_f
    integer :: cpuid, j, my_num_pairs, me
    real(dp) :: divisor

    if(.not.allocated(dTSRdr)) allocate(dTSRdr(3, 3, nat, 3)); dTSRdR = 0.0_DP
    if(.not.allocated(dTSRdV)) allocate(dTSRdV(3, 3, nat)); dTSRdV = 0.0_DP
    call start_clock('mbd_scs')

    call start_clock('mbd_pinit')
    call mbdvdw_para_init()
    call stop_clock('mbd_pinit')

    me = me_image+1

    my_num_pairs = n_pairs(me_image+1)
    offsets(1) = 0
    force_offsets(1) = 0
    sc_offsets(1) = 0
    force_offsets_h(1) = 0

    do counter = 2, nproc_image, 1
      offsets(counter) = offsets(counter-1) + n_pairs(counter-1)
      force_offsets(counter) = force_offsets(counter - 1) + n_comps_f(counter - 1)
      force_offsets_h(counter) = force_offsets_h(counter - 1) + n_comps_h(counter - 1)
      sc_offsets(counter) = sc_offsets(counter - 1) + n_comps_v(counter - 1)
      end do

    if(.not.allocated(my_tsr))           allocate(my_tsr(my_num_pairs, 3, 3)); my_tsr = 0.0_DP
    if(.not.allocated(collected_tsr))    allocate(collected_tsr(num_pairs, 3, 3)); collected_tsr = 0.0_DP

    if(vdw_self_consistent) then
      if(.not.allocated(my_dtsrdV))           allocate(my_dtsrdV(my_num_pairs, 3, 3, nat)); my_dtsrdV = 0.0_DP
      if(.not.allocated(collected_dtsrdV))    allocate(collected_dtsrdV(num_pairs, 3, 3, n_comps_v(me)));
      collected_dtsrdV = 0.0_DP
      end if

    if(do_forces) then
      if(.not.mbd_vdw_isolated) then
        if(.not.allocated(my_dtsrdh))           allocate(my_dtsrdh(my_num_pairs, 3, 3, 9)); my_dtsrdh = 0.0_DP
        if(.not.allocated(collected_dtsrdh))    allocate(collected_dtsrdh(num_pairs, 3, 3, n_comps_h(me)));
        collected_dtsrdh = 0.0_DP
        end if

      if(.not.allocated(my_dtsrdr))           allocate(my_dtsrdr(my_num_pairs, 3, 3, nat*3)); my_dtsrdr = 0.0_DP
      if(.not.allocated(collected_dtsrdr))    allocate(collected_dtsrdr(num_pairs, 3, 3, n_comps_f(me)));
      collected_dtsrdr = 0.0_DP
      end if

    call start_clock('mbd_lut_tsr')
    counter = 1
    do i = 1, num_pairs, 1
      p = pairs_scs(i)%p
      q = pairs_scs(i)%q
      cpuid = pairs_scs(i)%cpu

      if((cpuid+1).eq.me) then
        sl_mult = 1.0_DP
        call mbdvdw_TGG(1, p, q, nat, h_, ainv_, tau, Tsr, dTsrdR, dTsrdh, dTsrdV)

        my_tsr(counter,:,:) = tsr

        if(vdw_self_consistent) my_dtsrdv(counter,:,:,:) = dtsrdv

        if(do_forces) then
          do s = 1, nat, 1
            do i_f = 1, 3, 1
              my_dtsrdr(counter,:,:,3*(s-1)+i_f) = dtsrdr(:,:,s,i_f)
              end do
            end do
          if(.not.mbd_vdw_isolated) then
            do s = 1, 3, 1
              do i_f = 1, 3, 1
                my_dtsrdh(counter,:,:,3*(s-1)+i_f) = dtsrdh(:,:,s,i_f)
                end do
              end do
            end if
          end if
        counter = counter + 1
        end if
      end do
    call stop_clock('mbd_lut_tsr')

    !!!!!!!!!!!!!!!!!!
    ! Send and Receives for tsr
    !!!!!!!!!!!!!!!!!!
    do counter = 1, n_pairs(me), 1
      collected_tsr(counter + offsets(me_image+1), :, :) = my_tsr(counter, :, :)
      end do

    do i = 0, nproc_image-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          call mpi_send(my_tsr,&                  ! Send buffer
                  my_num_pairs*3*3,&        ! Count of communicated data
                  mpi_double_precision,&    ! Data type
                  i,&                       ! Destination
                  j,&                       ! Tag (who it was sent from)
                  intra_image_comm,&        ! Communicator to send with
                  ierr)                     ! Error flag
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_tsr))    allocate(temp_tsr(n_pairs(j+1), 3, 3))
            call mpi_recv(temp_tsr,&                    ! Recv buffer
                    n_pairs(j+1)*3*3,&            ! Count of communicated data
                    mpi_double_precision,&        ! Data type
                    j,&                           ! Source (who it was sent from)
                    j,&                           ! Tag (also who it was sent from)
                    intra_image_comm,&            ! Communicator to send with
                    rstatus,&                     ! Status flag
                    ierr)                         ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_tsr(counter + offsets(j+1), :, :) = temp_tsr(counter, :, :)
              end do
            if(allocated(temp_tsr)) deallocate(temp_tsr)
            end if
          end if
        end do
      end do

    if(allocated(my_tsr)) deallocate(my_tsr)

    if(do_forces) call mbdvdw_send_recv_forces(my_dtsrdr, collected_dtsrdr, max_proc_forces, &
                                               n_pairs, n_comps_f, force_offsets, offsets)
    if(.not.mbd_vdw_isolated) then
      if(do_forces) call mbdvdw_send_recv_forces(my_dtsrdh, collected_dtsrdh, max_proc_h, &
                                                 n_pairs, n_comps_h, force_offsets_h, offsets)
      end if
    if(vdw_self_consistent) call mbdvdw_send_recv_forces(my_dtsrdv, collected_dtsrdv, max_proc_sc,&
                                               n_pairs, n_comps_v, sc_offsets, offsets)

    call start_clock('mbd_a_mat')
    counter_amat = 1
    if(.not.allocated(A_matrix))        allocate(A_matrix(3*nat, 3*nat));                           A_matrix = 0.0_DP
    if(vdw_self_consistent) then
      if(.not.allocated(dA_matrixdV))     allocate(dA_matrixdV(3*nat, 3*nat, n_comps_v(me)));       dA_matrixdV = 0.0_DP
      if(.not.allocated(TsrdV))          allocate(TSRdV(3, 3, n_comps_v(me))); TsrdV = 0.0_DP
      end if

    if(do_forces) then
      if(.not.allocated(dA_matrixdR))     allocate(dA_matrixdR(3*nat, 3*nat, n_comps_f(me)));     dA_matrixdR = 0.0_DP
      if(.not.allocated(Tsrdr))           allocate(TSRdr(3, 3, n_comps_f(me))); Tsrdr = 0.0_DP
      if(.not.mbd_vdw_isolated) then
        if(.not.allocated(dA_matrixdh))     allocate(dA_matrixdh(3*nat, 3*nat, n_comps_h(me)));     dA_matrixdh = 0.0_DP
        if(.not.allocated(Tsrdh))           allocate(TSRdh(3, 3, n_comps_h(me))); Tsrdh = 0.0_DP
        end if
      end if

    ! only loops over the upper triangle
    counter = 1
    do p = 1, nat, 1
      do q = p, nat, 1
        TSR = collected_tsr(counter, :, :)
        if(vdw_self_consistent) TSRdv(:,:,:) = collected_dtsrdv(counter,:,:,:)
        if(do_forces) then
            TSRdr(:,:,:) = collected_dtsrdr(counter, :, :, :)
            if(.not.mbd_vdw_isolated) TSRdh(:,:,:) = collected_dtsrdh(counter, :, :, :)
        end if
        counter = counter + 1
        do i_idx = 1, 3, 1
          do j_idx = i_idx, 3, 1
            i_index = (3*p - 3 + i_idx)
            j_index = (3*q - 3 + j_idx)
            i_index_T = (3*p - 3 + j_idx)
            j_index_T = (3*q - 3 + i_idx)
            if(p.eq.q) then
              if(i_idx.eq.j_idx) then
                A_matrix(i_index, j_index) = 1.0_DP/alpha_ts(p)
                divisor = 1.0_DP/alpha_ts(p)**2.0_DP
                if(vdw_self_consistent) then
                  cnt_v = 1
                  do i_f = 1, nat, 1
                    if(v_cpu_id(i_f).eq.me_image) then
                      s = i_f
                      dA_matrixdV(i_index, j_index, cnt_v) = -dalpha_tsdV(p, i_f)*divisor
                      cnt_v = cnt_v + 1
                      end if ! checks if this MPI rank computes this component
                    end do ! Loop over all possible components
                  end if ! if to check if self consistent derivatives need to be computed
                if(do_forces) then
                  ! dR
                  cnt_f = 1
                  do i_f = 1, 3*nat, 1
                    if(f_cpu_id(i_f).eq.me_image) then
                      call mbdvdw_get_is(i_f, s, i)
                      dA_matrixdR(i_index, j_index, cnt_f) = -dalpha_tsdR(p, s, i)*divisor
                      cnt_f = cnt_f + 1
                      end if ! checks if this MPI rank computes this component
                    end do ! Loop over all possible components
                  ! dH
                  if(.not.mbd_vdw_isolated) then
                    cnt_h = 1
                    do i_f = 1, 9, 1
                      if(h_cpu_id(i_f).eq.me_image) then
                        call mbdvdw_get_is(i_f, s, i)
                        dA_matrixdh(i_index, j_index, cnt_h) = -dalpha_tsdh(p, s, i)*divisor
                        cnt_h = cnt_h + 1
                        end if ! checks if this MPI rank computes this component
                      end do ! Loop over all possible components
                    end if
                  end if ! if to check if forces and stresses need to be computed
                end if ! i_idx = j_idx if statement
              end if ! p = q iff statement
            A_matrix(i_index, j_index) = A_matrix(i_index, j_index) + TSR(i_idx, j_idx)
            A_matrix(j_index, i_index) = A_matrix(i_index, j_index)
            A_matrix(i_index_T, j_index_T) = A_matrix(i_index, j_index)
            A_matrix(j_index_T, i_index_T) = A_matrix(j_index, i_index)
            if(vdw_self_consistent) then
              cnt_v = 1
              do i_f = 1, nat, 1
                if(v_cpu_id(i_f).eq.me_image) then
                  dA_matrixdV(i_index,j_index,cnt_v)=dA_matrixdV(i_index,j_index,cnt_v)+TSRdV(i_idx,j_idx,cnt_v)
                  dA_matrixdV(j_index,i_index,cnt_v)=dA_matrixdV(i_index,j_index,cnt_v)
                  dA_matrixdV(i_index_T,j_index_T,cnt_v)=dA_matrixdV(i_index,j_index,cnt_v)
                  dA_matrixdV(j_index_T,i_index_T,cnt_v)=dA_matrixdV(j_index,i_index,cnt_v)
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              end if ! if to check if self consistent derivatives need to be computed

            if(do_forces) then
              ! dR
              cnt_f = 1
              do i_f = 1, 3*nat, 1
                if(f_cpu_id(i_f).eq.me_image) then
                  call mbdvdw_get_is(i_f, s, i)
                  dA_matrixdR(i_index,j_index,cnt_f)=dA_matrixdR(i_index,j_index,cnt_f)+TSRdR(i_idx,j_idx,cnt_f)
                  dA_matrixdR(j_index,i_index,cnt_f)=dA_matrixdR(i_index, j_index, cnt_f)
                  dA_matrixdR(i_index_T,j_index_T,cnt_f)=dA_matrixdR(i_index,j_index,cnt_f)
                  dA_matrixdR(j_index_T,i_index_T,cnt_f)=dA_matrixdR(j_index,i_index,cnt_f)
                  cnt_f = cnt_f + 1
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              ! dH
              if(.not.mbd_vdw_isolated) then
                cnt_h = 1
                do i_f = 1, 9, 1
                  if(h_cpu_id(i_f).eq.me_image) then
                    call mbdvdw_get_is(i_f, s, i)
                    dA_matrixdh(i_index,j_index,cnt_h)=dA_matrixdh(i_index, j_index,cnt_h)+TSRdh(i_idx,j_idx,cnt_h)
                    dA_matrixdh(j_index,i_index,cnt_h)=dA_matrixdh(i_index, j_index,cnt_h) !fill in symmetric elements
                    dA_matrixdh(i_index_T,j_index_T,cnt_h)=dA_matrixdh(i_index,j_index,cnt_h)
                    dA_matrixdh(j_index_T,i_index_T,cnt_h)=dA_matrixdh(j_index,i_index,cnt_h)
                    cnt_h = cnt_h + 1
                    end if ! checks if this MPI rank computes this component
                  end do ! Loop over all possible components
                end if
              end if ! if to check if forces and stresses need to be computed
            end do ! j_idx loop
          end do ! i_idx loop
        end do ! q loop
      end do ! q loop
    call stop_clock('mbd_a_mat')
    ! Performs the LU Decomposition of the A_matrix in preparation for inverting it
    call DGETRF(3*nat, 3*nat, A_matrix, 3*nat, IPIV, errorflag)
      if(errorflag.ne.0 .and. me_image == root_image) then
        write(use_unit, '(3X, "ERROR IN LU DECOMPOSITION")')
        end if
     call DGETRI(3*nat, A_matrix, 3*nat, IPIV, WORK,3*nat,errorflag )
      if(errorflag.ne.0 .and. me_image == root_image) then
        write(use_unit, '(3X, "ERROR IN INVERSION")')
        end if
    ! At this point the A_matrix is now screened.  We complete the derivative of A by again
    ! applying the formula for the derivative of an inverse: \partial\bar{A} = -\bar{A} [\partial inv(\bar{A}) ] \bar{A}
    call start_clock('mbd_a_force')

    if(vdw_self_consistent) then
      cnt_v = 1
      do i_f = 1, nat, 1
        if(v_cpu_id(i_f).eq.me_image) then
          call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                dA_matrixdV(:, :, cnt_v), 3*nat,&
                A_matrix, 3*nat, 0.0_DP, temp1, 3*nat)

          call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                A_matrix, 3*nat,&
                temp1, 3*nat, 0.0_DP, temp, 3*nat)
          dA_matrixdV(:, :, cnt_v) = -temp
          cnt_v = cnt_v + 1
          end if ! checks if this MPI rank computes this component
        end do ! Loop over all possible components
      end if ! if to check if self consistent derivatives need to be computed

    if(do_forces) then
      ! dR
      cnt_f = 1
      do i_f = 1, 3*nat, 1
        if(f_cpu_id(i_f).eq.me_image) then
          call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                dA_matrixdR(:, :, cnt_f), 3*nat,&
                A_matrix, 3*nat, 0.0_DP, temp1, 3*nat)

          call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                A_matrix, 3*nat,&
                temp1, 3*nat, 0.0_DP, temp, 3*nat)
          dA_matrixdR(:, :, cnt_f) = -temp
          cnt_f = cnt_f + 1
          end if ! checks if this MPI rank computes this component
        end do ! Loop over all possible components
      ! dH
      if(.not.mbd_vdw_isolated) then
        cnt_h = 1
        do i_f = 1, 9, 1
          if(h_cpu_id(i_f).eq.me_image) then
            call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                  dA_matrixdh(:, :, cnt_h), 3*nat,&
                  A_matrix, 3*nat, 0.0_DP, temp1, 3*nat)

            call dgemm('N', 'N', 3*nat, 3*nat, 3*nat, 1.0_DP, &
                  A_matrix, 3*nat,&
                  temp1, 3*nat, 0.0_DP, temp, 3*nat)
            dA_matrixdh(:, :, cnt_h) = -temp
            cnt_h = cnt_h + 1
            end if ! checks if this MPI rank computes this component
          end do ! Loop over all possible components
        end if
      end if ! if to check if forces and stresses need to be computed
    call stop_clock('mbd_a_force')
    call stop_clock('mbd_scs')

    return
    end subroutine mbdvdw_SCS

  subroutine mbdvdw_construct_hamiltonian_complex()
    implicit none
    complex(dp), dimension(3, 3, nk) :: TLR
    complex(dp), dimension(:, :, :, :, :), allocatable :: dTLRdr
    complex(dp), dimension(3, 3, 3, 3, nk) :: dTLRdh
    complex(dp), dimension(:,:,:,:), allocatable :: dTLRdV
    complex(dp), dimension(:, :, :), allocatable    :: Tlrdh, Tlrdr, TlrdV
    complex(dp), dimension(:,:,:,:), allocatable :: my_tlr, collected_tlr, temp_tlr
    complex(dp), dimension(:,:,:,:), allocatable :: my_dtlrdr, collected_dtlrdr
    complex(dp), dimension(:,:,:,:), allocatable :: my_dtlrdh, collected_dtlrdh
    complex(dp), dimension(:,:,:,:),   allocatable :: my_dtlrdV, collected_dtlrdV
    complex(dp) :: prefactor, pre1, pre2, der, aSqrt, omMult, a, self_energy

    !MPI Vars
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    integer, dimension(nproc_image) :: offsets, force_offsets, force_offsets_h, sc_offsets

    ! Looping varlables
    integer :: p, q, i_idx, j_idx, i_index, j_index, i, s, counter, i_f, cnt_v, cnt_f, cnt_h
    integer :: i_index_T, j_index_T
    integer :: cpuid, j, my_num_pairs, cnt, i_k

    a = ewald_cutoff
    self_energy = 4.0_dp*a**3.0_dp/(3.0_dp*SQRTPI)

    ! Supercell variables
    call start_clock('mbd_construct_hamiltonian')
    if(do_recip) supercell_cutoff = 1.0_DP
    call mbdvdw_construct_sl()

    ruc = compute_recucell(h_sl)
    if(do_ewald) then
      call mbdvdw_cellvol(h_, nu)
      call mbdvdw_ewald_params(nu)
      if(do_forces) then
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_dcellvoldh(ainv_, s, i, nu, dnudh(s, i))
            end do
          end do
        call mbdvdw_dewald_params_dh(ainv_, nu)
        end if
      end if

    call start_clock('mbd_pinit_sl')
    call mbdvdw_para_init_sl()
    call stop_clock('mbd_pinit_sl')

    call start_clock('mbd_lut_tlr')
    my_num_pairs = n_pairs(me_image+1)
    if(.not.allocated(dTLRdr)) allocate(dTLRdr(3, 3, nat, 3, nk)); dTLRdR = 0.0_DP
    if(.not.allocated(dTLRdV)) allocate(dTLRdV(3, 3, nat, nk)); dTLRdV = 0.0_DP

    offsets(1) = 0
    force_offsets(1) = 0
    force_offsets_h(1) = 0
    sc_offsets(1) = 0
    do counter = 2, nproc_image, 1
      offsets(counter) = offsets(counter-1) + n_pairs(counter-1)
      force_offsets(counter) = force_offsets(counter - 1) + n_comps_f(counter - 1)
      force_offsets_h(counter) = force_offsets_h(counter - 1) + n_comps_h(counter - 1)
      sc_offsets(counter) = sc_offsets(counter - 1) + n_comps_v(counter - 1)
      end do

    if(.not.allocated(my_tlr))           allocate(my_tlr(n_pairs(me_image+1), 3, 3, nk)); my_tlr = 0.0_DP
    if(.not.allocated(collected_tlr))    allocate(collected_tlr(num_pairs, 3, 3, nk)); collected_tlr = 0.0_DP

    if(do_forces) then
      if(.not.mbd_vdw_isolated) then
        if(.not.allocated(my_dtlrdh))           allocate(my_dtlrdh(n_pairs(me_image+1), 3, 3, 9*nk)); my_dtlrdh = 0.0_DP
        if(.not.allocated(collected_dtlrdh))    allocate(collected_dtlrdh(num_pairs, 3, 3, n_comps_h(me_image+1)));
        collected_dtlrdh = 0.0_DP
        end if

      if(.not.allocated(my_dtlrdr))           allocate(my_dtlrdr(n_pairs(me_image+1), 3, 3, nk*nat*3)); my_dtlrdr = 0.0_DP
      if(.not.allocated(collected_dtlrdr))    allocate(collected_dtlrdr(num_pairs, 3, 3, n_comps_f(me_image+1)));
      collected_dtlrdr = 0.0_DP
      end if

    if(vdw_self_consistent) then
      if(.not.allocated(my_dtlrdv))           allocate(my_dtlrdv(n_pairs(me_image+1), 3, 3, nk*nat)); my_dtlrdv = 0.0_DP
      if(.not.allocated(collected_dtlrdv))    allocate(collected_dtlrdv(num_pairs, 3, 3, n_comps_v(me_image+1)));
      collected_dtlrdv = 0.0_DP
      end if

    counter = 1
    do i = 1, num_pairs, 1
      p = unique_pairs(i)%p
      q = unique_pairs(i)%q
      cpuid = unique_pairs(i)%cpu

      if(cpuid.eq.me_image) then
        if(do_ewald) then
          call mbdvdw_ewaldsum(p, q, nat_sl, h_sl, ainv_sl, tau_sl, TLR, dTLRdR, dTLRdh, dTLRdV)
          else
          call mbdvdw_TGG_complex(2, p, q, nat_sl, h_sl, ainv_sl, tau_sl, TLR, dTLRdR, dTLRdh, dTLRdV)
          end if

        my_tlr(counter, :, :, :) = tlr
        if(vdw_self_consistent) then
          do s = 1, nat, 1
            do i_k = 1, nk, 1
              cnt = nk*(s-1) + i_k
              my_dtlrdV(counter,:,:,cnt) = dTlrdV(:,:,s,i_k)
              end do
            end do
          end if
        if(do_forces) then
          do s = 1, nat, 1
            do i_f = 1, 3, 1
              do i_k = 1, nk, 1
                cnt = i_k + 3*nk*(s-1) + nk*(i_f-1)
                my_dtlrdr(counter,:,:,cnt) = dtlrdr(:,:,s,i_f,i_k)
                end do
              end do
            end do
          if(.not.mbd_vdw_isolated) then
            do s = 1, 3, 1
              do i_f = 1, 3, 1
                do i_k = 1, nk, 1
                  cnt = i_k + 3*nk*(s-1) + nk*(i_f-1)
                  my_dtlrdh(counter,:,:,cnt) = dtlrdh(:,:,s,i_f,i_k)
                  end do
                end do
              end do
            end if
          end if
        counter = counter + 1
        end if
      end do
    call stop_clock('mbd_lut_tlr')

    call start_clock('mbd_cpq_comm')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Send and Receives for tlr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do counter = 1, n_pairs(me_image+1), 1
      collected_tlr(counter + offsets(me_image+1), :, :, :) = my_tlr(counter, :, :, :)
      end do

    do i = 0, nproc_image-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          call mpi_send(my_tlr,&                  ! Send buffer
                  my_num_pairs*3*3*nk,&        ! Count of communicated data
                  mpi_complex16,&    ! Data type
                  i,&                       ! Destination
                  j,&                       ! Tag (who it was sent from)
                  intra_image_comm,&        ! Communicator to send with
                  ierr)                     ! Error flag
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_tlr))    allocate(temp_tlr(n_pairs(j+1), 3, 3, nk))
            call mpi_recv(temp_tlr,&                    ! Recv buffer
                    n_pairs(j+1)*3*3*nk,&            ! Count of communicated data
                    mpi_complex16,&        ! Data type
                    j,&                           ! Source (who it was sent from)
                    j,&                           ! Tag (also who it was sent from)
                    intra_image_comm,&            ! Communicator to send with
                    rstatus,&                     ! Status flag
                    ierr)                         ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_tlr(counter + offsets(j+1), :, :, :) = temp_tlr(counter, :, :, :)
              end do
            if(allocated(temp_tlr)) deallocate(temp_tlr)
            end if
          end if
        end do
      end do
    if(allocated(my_tlr)) deallocate(my_tlr)

    if(do_forces) call mbdvdw_send_recv_forces_complex(my_dtlrdr, collected_dtlrdr, max_proc_forces, &
                                               n_pairs, n_comps_f, force_offsets, offsets)
    if(do_forces.and.(.not.mbd_vdw_isolated)) call mbdvdw_send_recv_forces_complex(my_dtlrdh, collected_dtlrdh, &
                                              max_proc_h, n_pairs, n_comps_h, force_offsets_h, offsets)
    if(vdw_self_consistent) call mbdvdw_send_recv_forces_complex(my_dtlrdv, collected_dtlrdv, max_proc_sc,&
                                               n_pairs, n_comps_v, sc_offsets, offsets)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Putting it all together!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call stop_clock('mbd_cpq_comm')
      ! only loops over the upper triangle
      if(.not.allocated(Cpq_c))             allocate(Cpq_c(3*nat_sl, 3*nat_sl, nk));              Cpq_c = 0.0_DP
      call start_clock('mbd_build_c')
      if(vdw_self_consistent) then
        if(.not.allocated(TLRdV))           allocate(TLRdV(3, 3, n_comps_v(me_image+1))); TLRdV = 0.0_DP
        if(.not.allocated(dCpqdV_c))          allocate(dCpqdV_c(3*nat_sl, 3*nat_sl, n_comps_v(me_image+1)));   dCpqdV_c = 0.0_DP
        end if
      if(do_forces) then
        if(.not.allocated(dCpqdR_c))          allocate(dCpqdR_c(3*nat_sl, 3*nat_sl, n_comps_f(me_image+1)));   dCpqdR_c = 0.0_DP
        if(.not.allocated(Tlrdr))           allocate(TLRdr(3, 3, n_comps_f(me_image+1))); TlrdR = 0.0_DP
        if(.not.mbd_vdw_isolated) then
          if(.not.allocated(dCpqdh_c))          allocate(dCpqdh_c(3*nat_sl, 3*nat_sl, n_comps_h(me_image+1)));   dCpqdh_c = 0.0_DP
          if(.not.allocated(Tlrdh))           allocate(TLRdh(3, 3, n_comps_h(me_image+1))); TlrdH = 0.0_DP
          end if
        end if
    counter = 1
    do p = 1, nat_sl, 1
      do q = p, nat_sl, 1
        aSqrt = dsqrt(alpha_0_sl(p)*alpha_0_sl(q))
        omMult = omega_scs_sl(p)*omega_scs_sl(q)
        prefactor = omMult*aSqrt
        TLR = collected_tlr(counter, :, :, :)
        if(do_forces) TLRdr = collected_dtlrdr(counter, :, :, :)
        if(do_forces.and.(.not.mbd_vdw_isolated)) TLRdh = collected_dtlrdh(counter, :, :, :)
        if(vdw_self_consistent) TLRdV = collected_dtlrdV(counter, :, :, :)
        counter = counter + 1
        do i_idx = 1, 3, 1
          do j_idx = i_idx, 3, 1
            i_index = (3*p - 3 + i_idx)
            j_index = (3*q - 3 + j_idx)
            i_index_T = (3*p - 3 + j_idx)
            j_index_T = (3*q - 3 + i_idx)
            if(p.eq.q) then
              if(i_idx.eq.j_idx) then
                do i_k = 1, nk, 1
                  Cpq_c(i_index, j_index, i_k) = omega_scs_sl(p)**2.0_DP
                  end do ! i_k loop
                if(vdw_self_consistent) then
                  cnt_v = 1
                  do i_f = 1, nat*nk, 1
                    if(v_cpu_id(i_f).eq.me_image) then
                      call mbdvdw_get_ks(i_f, s, i_k)
                      dCpqdV_c(i_index, j_index, cnt_v) = 2.0_DP*omega_scs_sl(p)*domegadV_sl(p, s)
                      cnt_v = cnt_v + 1
                      end if ! checks if this MPI rank computes this component
                    end do ! Loop over all possible components
                  end if ! if to check if self consistent derivatives need to be computed
                if(do_forces) then
                  ! dR
                  cnt_f = 1
                  do i_f = 1, 3*nat*nk, 1
                    if(f_cpu_id(i_f).eq.me_image) then
                      call mbdvdw_get_isk(i_f, s, i, i_k)
                      dCpqdR_c(i_index,j_index,cnt_f)=2.0_DP*omega_scs_sl(p)*domegadR_sl(p,s,i)
                      cnt_f = cnt_f + 1
                      end if ! checks if this MPI rank computes this component
                    end do ! Loop over all possible components
                  ! dH
                  if(.not.mbd_vdw_isolated) then
                    cnt_h = 1
                    do i_f = 1, 9*nk, 1
                      if(h_cpu_id(i_f).eq.me_image) then
                        call mbdvdw_get_isk(i_f, s, i, i_k)
                        dCpqdh_c(i_index,j_index,cnt_h)=2.0_DP*omega_scs_sl(p)*domegadh_sl(p,s,i)
                        cnt_h = cnt_h + 1
                        end if ! checks if this MPI rank computes this component
                      end do ! Loop over all possible components
                    end if
                  end if ! if to check if forces and stresses need to be computed
                end if ! i_idx = j_idx if statement
              end if ! p = q iff statement
            do i_k = 1, nk, 1
              Cpq_c(i_index, j_index, i_k) = Cpq_c(i_index, j_index, i_k) +  prefactor*TLR(i_idx, j_idx, i_k)
              Cpq_c(j_index, i_index, i_k) = conjg(Cpq_c(i_index, j_index, i_k))
              Cpq_c(i_index_T, j_index_T, i_k) = Cpq_c(i_index, j_index, i_k)
              Cpq_c(j_index_T, i_index_T, i_k) = Cpq_c(j_index, i_index, i_k)
              end do ! i_k loop
            if(vdw_self_consistent) then
              cnt_v = 1
              do i_f = 1, nat*nk, 1
                if(v_cpu_id(i_f).eq.me_image) then
                  call mbdvdw_get_ks(i_f, s, i_k)
                  pre1=(omega_scs_sl(p)*domegadV_sl(q,s)+omega_scs_sl(q)*domegadV_sl(p,s))*asqrt
                  pre2=omMult*(alpha_0_sl(q)*dalpha_0dV_sl(p,s)+alpha_0_sl(p)*dalpha_0dV_sl(q,s))
                  pre2=pre2/(2.0_DP*asqrt)
                  der = (pre1+pre2)*TLR(i_idx,j_idx,i_k)+prefactor*TLRdV(i_idx,j_idx,cnt_v)
                  dCpqdV_c(i_index,j_index,cnt_v)=dCpqdV_c(i_index,j_index,cnt_v)+der
                  dCpqdV_c(j_index,i_index,cnt_v)=conjg(dCpqdV_c(i_index,j_index,cnt_v))
                  dCpqdV_c(i_index_T,j_index_T,cnt_v) = dCpqdV_c(i_index,j_index,cnt_v)
                  dCpqdV_c(j_index_T,i_index_T,cnt_v) = dCpqdV_c(j_index,i_index,cnt_v)
                  cnt_v = cnt_v + 1
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              end if ! if to check if self consistent derivatives need to be computed

            if(do_forces) then
              ! dR
              cnt_f = 1
              do i_f = 1, 3*nat*nk, 1
                if(f_cpu_id(i_f).eq.me_image) then
                  call mbdvdw_get_isk(i_f, s, i, i_k)
                  pre1=omega_scs_sl(p)*domegadR_sl(q,s,i)+omega_scs_sl(q)*domegadR_sl(p,s,i)
                  pre1=pre1*asqrt
                  pre2=omMult
                  pre2=pre2*(alpha_0_sl(q)*dalpha_0dR_sl(p,s,i)+alpha_0_sl(p)*dalpha_0dR_sl(q,s,i))
                  pre2=pre2/(2.0_DP*asqrt)
                  der = (pre1+pre2)*TLR(i_idx,j_idx,i_k)+prefactor*TLRdR(i_idx,j_idx,cnt_f)
                  dCpqdR_c(i_index,j_index,cnt_f)=dCpqdR_c(i_index,j_index,cnt_f)+der
                  dCpqdR_c(j_index,i_index,cnt_f)=conjg(dCpqdR_c(i_index,j_index,cnt_f))
                  dCpqdR_c(i_index_T,j_index_T,cnt_f) = dCpqdR_c(i_index,j_index,cnt_f)
                  dCpqdR_c(j_index_T,i_index_T,cnt_f) = dCpqdR_c(j_index,i_index,cnt_f)
                  cnt_f = cnt_f + 1
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              ! dH
              if(.not.mbd_vdw_isolated) then
                cnt_h = 1
                do i_f = 1, 9*nk, 1
                  if(h_cpu_id(i_f).eq.me_image) then
                    call mbdvdw_get_isk(i_f, s, i, i_k)
                    pre1=omega_scs_sl(p)*domegadh_sl(q,s,i)+omega_scs_sl(q)*domegadh_sl(p,s,i)
                    pre1=pre1*asqrt
                    pre2=omMult
                    pre2=pre2*(alpha_0_sl(q)*dalpha_0dh_sl(p,s,i)+alpha_0_sl(p)*dalpha_0dh_sl(q,s,i))
                    pre2=pre2/(2.0_DP*asqrt)
                    der = (pre1+pre2)*TLR(i_idx,j_idx,i_k)+prefactor*TLRdh(i_idx,j_idx,cnt_h)
                    dCpqdh_c(i_index,j_index,cnt_h)=dCpqdh_c(i_index,j_index,cnt_h)+der
                    dCpqdh_c(j_index,i_index,cnt_h)=conjg(dCpqdh_c(i_index,j_index,cnt_h))
                    dCpqdh_c(i_index_T,j_index_T,cnt_h) = dCpqdh_c(i_index,j_index,cnt_h)
                    dCpqdh_c(j_index_T,i_index_T,cnt_h) = dCpqdh_c(j_index,i_index,cnt_h)
                    cnt_h = cnt_h + 1
                    end if ! checks if this MPI rank computes this component
                  end do ! Loop over all possible components
                end if
              end if ! if to check if forces and stresses need to be computed
            end do ! j_idx loop
          end do ! i_idx loop
        end do ! q loop
      end do ! q loop
    call stop_clock('mbd_build_c')
    if(allocated(collected_tlr)) deallocate(collected_tlr)
    if(allocated(collected_dtlrdr)) deallocate(collected_dtlrdr)
    if(allocated(collected_dtlrdh)) deallocate(collected_dtlrdh)
    if(allocated(collected_dtlrdV)) deallocate(collected_dtlrdv)

    call stop_clock('mbd_construct_hamiltonian')
    return
    end subroutine mbdvdw_construct_hamiltonian_complex

  subroutine mbdvdw_construct_sl()
    implicit none
    real(dp) :: phi

    ! Looping varlables
    integer :: p, q, i_idx, j_idx, k_idx, s, i

    if(.not.do_recip) then
      if(mbd_vdw_verbosity.ge.0 .and. me_image == root_image) then
        write(use_unit, *) 'using a superlattice of size ', sl_i, sl_j, sl_k
        write(use_unit, *) 'Number of atoms in the superlattice is', nat_sl
        end if
      end if

    if(.not.allocated(alpha_0_sl))      allocate(alpha_0_sl(nat_sl));                   alpha_0_sl = 0.0_DP
    if(.not.allocated(omega_scs_sl))    allocate(omega_scs_sl(nat_sl));                 omega_scs_sl = 0.0_DP
    if(.not.allocated(R_MBD_VdW_sl))    allocate(R_MBD_VdW_sl(nat_sl));                 R_MBD_VdW = 0.0_DP

    if(vdw_self_consistent) then
      if(.not.allocated(dalpha_0dV_sl))   allocate(dalpha_0dV_sl(nat_sl, nat));        dalpha_0dV_sl = 0.0_DP
      if(.not.allocated(domegadV_sl))     allocate(domegadV_sl(nat_sl, nat));          domegadV_sl = 0.0_DP
      if(.not.allocated(dR_MBD_VdWdV_sl)) allocate(dR_MBD_VdWdV_sl(nat_sl, nat));      dR_MBD_VdWdV_sl = 0.0_DP
      end if

    if(do_forces) then
      if(.not.allocated(dalpha_0dR_sl))   allocate(dalpha_0dR_sl(nat_sl, nat, 3));        dalpha_0dR_sl = 0.0_DP
      if(.not.allocated(domegadR_sl))     allocate(domegadR_sl(nat_sl, nat, 3));          domegadR_sl = 0.0_DP
      if(.not.allocated(dR_MBD_VdWdR_sl)) allocate(dR_MBD_VdWdR_sl(nat_sl, nat, 3));      dR_MBD_VdWdR_sl = 0.0_DP
      if(.not.mbd_vdw_isolated) then
        if(.not.allocated(dR_MBD_VdWdh_sl)) allocate(dR_MBD_VdWdh_sl(nat_sl, 3, 3));        dR_MBD_VdWdh_sl = 0.0_DP
        if(.not.allocated(domegadh_sl))     allocate(domegadh_sl(nat_sl, 3, 3));            domegadh_sl = 0.0_DP
        if(.not.allocated(dalpha_0dh_sl))   allocate(dalpha_0dh_sl(nat_sl, 3, 3));          dalpha_0dh_sl = 0.0_DP
        end if
    end if
    if(.not.allocated(tau_sl))          allocate(tau_sl(3, nat_sl)); tau_sl = 0.0_DP
    if(.not.allocated(tau_sl_s))          allocate(tau_sl_s(3, nat_sl)); tau_sl_s = 0.0_DP
    h_sl = 0.0_DP
    ainv_sl = 0.0_DP

    h_sl(:, 1) = h_(:, 1)*DBLE(sl_i)
    h_sl(:, 2) = h_(:, 2)*DBLE(sl_j)
    h_sl(:, 3) = h_(:, 3)*DBLE(sl_k)

    call inv3x3_mat(h_sl, ainv_sl)

    sl_mult(1) = DBLE(sl_i)
    sl_mult(2) = DBLE(sl_j)
    sl_mult(3) = DBLE(sl_k)

    if(allocated(orig_idx)) deallocate(orig_idx)
    allocate(orig_idx(nat_sl))

    q = 1
    do i_idx = 0, sl_i-1, 1
      do j_idx = 0, sl_j-1, 1
        do k_idx = 0, sl_k-1, 1
          do p=1, nat, 1
            orig_idx(q) = p
            phi = R_vdw_free(p)/(alpha_free(p))**(1.0_DP/3.0_DP)

            R_MBD_VdW_sl(q)=phi*alpha_0(p)**(1.0_DP/3.0_DP)
            alpha_0_sl(q) = alpha_0(p)
            omega_scs_sl(q) = omega_scs(p)

            if(vdw_self_consistent) then
              dalpha_0dV_sl(q, :) = dalpha_0dV(p, :)
              domegadV_sl(q, :) = domegadV(p, :)
              dR_MBD_VdWdV_sl(q,:) = phi/(3.0_DP*alpha_0(p)**(2.0_DP/3.0_DP))*dalpha_0dV(p, :)
              end if

            if(do_forces) then
              dalpha_0dR_sl(q, :, :) = dalpha_0dR(p, :, :)
              domegadR_sl(q, :, :) = domegadR(p, :, :)
              dR_MBD_VdWdR_sl(q,:,:) = phi/(3.0_DP*alpha_0(p)**(2.0_DP/3.0_DP))*dalpha_0dR(p, :, :)
              if(.not.mbd_vdw_isolated) then
                dalpha_0dh_sl(q, :, :) = dalpha_0dh(p, :, :)
                domegadh_sl(q, :, :) = domegadh(p, :, :)
                dR_MBD_VdWdh_sl(q,:,:) = phi/(3.0_DP*alpha_0(p)**(2.0_DP/3.0_DP))*dalpha_0dh(p, :, :)
                end if
              end if

            tau_sl(:, q) = tau(:, p) + i_idx*h_(:, 1) + j_idx*h_(:, 2) + k_idx*h_(:, 3)
            q = q + 1
            end do
          end do
        end do
      end do

      if(do_ewald) then
        ruc = compute_recucell(h_sl)
        call mbdvdw_cellvol(h_sl, nu)
        do s = 1, 3, 1
          do i = 1, 3, 1
            call mbdvdw_dcellvoldh(ainv_, s, i, nu, dnudh(s, i))
            end do
          end do
        call mbdvdw_ewald_params(nu)
        call mbdvdw_dewald_params_dh(ainv_, nu)

        end if

    end subroutine mbdvdw_construct_sl

  subroutine mbdvdw_construct_hamiltonian_real()
    implicit none
    real(dp), dimension(3, 3) :: TLR
    real(dp), dimension(:, :, :, :), allocatable :: dTLRdr
    real(dp), dimension(3, 3, 3  , 3) :: dTLRdh
    real(dp), dimension(:, :, :), allocatable    :: dTLRdV, Tlrdh, Tlrdr, TlrdV
    real(dp), dimension(:,:,:), allocatable :: my_tlr, collected_tlr, temp_tlr
    real(dp), dimension(:,:,:,:), allocatable :: my_dtlrdr, collected_dtlrdr
    real(dp), dimension(:,:,:,:), allocatable :: my_dtlrdh, collected_dtlrdh
    real(dp), dimension(:,:,:,:),   allocatable :: my_dtlrdV, collected_dtlrdV
    real(dp) :: prefactor, pre1, pre2, der, aSqrt, omMult

    !MPI Vars
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    integer, dimension(nproc_image) :: offsets, force_offsets, force_offsets_h, sc_offsets

    ! Looping varlables
    integer :: p, q, i_idx, j_idx, i_index, j_index, i, s, counter, i_f, cnt_v, cnt_f, cnt_h
    integer :: cpuid, j, my_num_pairs, i_index_T, j_index_T

    ! Supercell variables
    call start_clock('mbd_construct_hamiltonian')
    call mbdvdw_construct_sl()

    call start_clock('mbd_pinit_sl')
    call mbdvdw_para_init_sl()
    call stop_clock('mbd_pinit_sl')

    call start_clock('mbd_lut_tlr')
    my_num_pairs = n_pairs(me_image+1)
    if(.not.allocated(dTLRdr)) allocate(dTLRdr(3, 3, nat, 3)); dTLRdR = 0.0_DP
    if(.not.allocated(dTLRdV)) allocate(dTLRdV(3, 3, nat)); dTLRdV = 0.0_DP

    offsets(1) = 0
    force_offsets(1) = 0
    force_offsets_h(1) = 0
    sc_offsets(1) = 0
    do counter = 2, nproc_image, 1
      offsets(counter) = offsets(counter-1) + n_pairs(counter-1)
      force_offsets(counter) = force_offsets(counter - 1) + n_comps_f(counter - 1)
      force_offsets_h(counter) = force_offsets_h(counter - 1) + n_comps_h(counter - 1)
      sc_offsets(counter) = sc_offsets(counter - 1) + n_comps_v(counter - 1)
      end do

    if(.not.allocated(my_tlr))           allocate(my_tlr(n_pairs(me_image+1), 3, 3)); my_tlr = 0.0_DP
    if(.not.allocated(collected_tlr))    allocate(collected_tlr(num_pairs, 3, 3)); collected_tlr = 0.0_DP

    if(do_forces) then
      if(.not.mbd_vdw_isolated) then
        if(.not.allocated(my_dtlrdh))           allocate(my_dtlrdh(n_pairs(me_image+1), 3, 3, 9)); my_dtlrdh = 0.0_DP
        if(.not.allocated(collected_dtlrdh))    allocate(collected_dtlrdh(num_pairs, 3, 3, n_comps_h(me_image+1)));
        collected_dtlrdh = 0.0_DP
        end if

      if(.not.allocated(my_dtlrdr))           allocate(my_dtlrdr(n_pairs(me_image+1), 3, 3, nat*3)); my_dtlrdr = 0.0_DP
      if(.not.allocated(collected_dtlrdr))    allocate(collected_dtlrdr(num_pairs, 3, 3, n_comps_f(me_image+1)));
      collected_dtlrdr = 0.0_DP
      end if

    if(vdw_self_consistent) then
      if(.not.allocated(my_dtlrdv))           allocate(my_dtlrdv(n_pairs(me_image+1), 3, 3, nat)); my_dtlrdv = 0.0_DP
      if(.not.allocated(collected_dtlrdv))    allocate(collected_dtlrdv(num_pairs, 3, 3, n_comps_v(me_image+1)));
      collected_dtlrdv = 0.0_DP
      end if

    counter = 1
    do i = 1, num_pairs, 1
      p = unique_pairs(i)%p
      q = unique_pairs(i)%q
      cpuid = unique_pairs(i)%cpu

      if(cpuid.eq.me_image) then
        if(do_ewald) then
          call mbdvdw_ewaldsum_real(p, q, nat_sl, h_sl, ainv_sl, tau_sl, TLR, dTLRdR, dTLRdh, dTLRdV)
          else
          call mbdvdw_TGG(2, p, q, nat_sl, h_sl, ainv_sl, tau_sl, TLR, dTLRdR, dTLRdh, dTLRdV)
          end if

        my_tlr(counter, :, :) = tlr
        if(vdw_self_consistent) my_dtlrdV(counter, :, :, :) = dTlrdV(:,:,:)
        if(do_forces) then
          do s = 1, nat, 1
            do i_f = 1, 3, 1
              my_dtlrdr(counter,:,:,3*(s-1)+i_f) = dtlrdr(:,:,s,i_f)
              end do
            end do
          if(.not.mbd_vdw_isolated) then
            do s = 1, 3, 1
              do i_f = 1, 3, 1
                my_dtlrdh(counter,:,:,3*(s-1)+i_f) = dtlrdh(:,:,s,i_f)
                end do
              end do
            end if
          end if
        counter = counter + 1
        end if
      end do
    call stop_clock('mbd_lut_tlr')

    call start_clock('mbd_cpq_comm')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Send and Receives for tlr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do counter = 1, n_pairs(me_image+1), 1
      collected_tlr(counter + offsets(me_image+1), :, :) = my_tlr(counter, :, :)
      end do

    do i = 0, nproc_image-1, 1                      ! Outer loop is the loop for recvs
      do j = 0, max_proc_pairs-1, 1                   ! Loop for recvs
        if((j.eq.me_image).and.(i.ne.me_image)) then
          call mpi_send(my_tlr,&                  ! Send buffer
                  my_num_pairs*3*3,&        ! Count of communicated data
                  mpi_double_precision,&    ! Data type
                  i,&                       ! Destination
                  j,&                       ! Tag (who it was sent from)
                  intra_image_comm,&        ! Communicator to send with
                  ierr)                     ! Error flag
        else
          if((i.eq.me_image).and.(j.ne.me_image)) then
            if(.not.allocated(temp_tlr))    allocate(temp_tlr(n_pairs(j+1), 3, 3))
            call mpi_recv(temp_tlr,&                    ! Recv buffer
                    n_pairs(j+1)*3*3,&            ! Count of communicated data
                    mpi_double_precision,&        ! Data type
                    j,&                           ! Source (who it was sent from)
                    j,&                           ! Tag (also who it was sent from)
                    intra_image_comm,&            ! Communicator to send with
                    rstatus,&                     ! Status flag
                    ierr)                         ! Error flag
            do counter = 1, n_pairs(j+1), 1
              collected_tlr(counter + offsets(j+1), :, :) = temp_tlr(counter, :, :)
              end do
            if(allocated(temp_tlr)) deallocate(temp_tlr)
            end if
          end if
        end do
      end do
    if(allocated(my_tlr)) deallocate(my_tlr)

    if(do_forces) call mbdvdw_send_recv_forces(my_dtlrdr, collected_dtlrdr, max_proc_forces, &
                                               n_pairs, n_comps_f, force_offsets, offsets)
    if(.not.mbd_vdw_isolated) then
      if(do_forces) call mbdvdw_send_recv_forces(my_dtlrdh, collected_dtlrdh, max_proc_h, &
                                                 n_pairs, n_comps_h, force_offsets_h, offsets)
      end if
    if(vdw_self_consistent) call mbdvdw_send_recv_forces(my_dtlrdv, collected_dtlrdv, max_proc_sc,&
                                               n_pairs, n_comps_v, sc_offsets, offsets)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Putting it all together!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call stop_clock('mbd_cpq_comm')
      ! only loops over the upper triangle
      if(.not.allocated(Cpq))             allocate(Cpq(3*nat_sl, 3*nat_sl));              Cpq = 0.0_DP
      call start_clock('mbd_build_c')
      if(vdw_self_consistent) then
        if(.not.allocated(TLRdV))           allocate(TLRdV(3, 3, n_comps_v(me_image+1))); TLRdV = 0.0_DP
        if(.not.allocated(dCpqdV))          allocate(dCpqdV(3*nat_sl, 3*nat_sl, n_comps_v(me_image+1)));   dCpqdV = 0.0_DP
        end if
      if(do_forces) then
        if(.not.allocated(dTLRdr))          allocate(dTLRdr(3, 3, n_comps_f(me_image+1), 3))
        if(.not.allocated(dCpqdR))          allocate(dCpqdR(3*nat_sl, 3*nat_sl, n_comps_f(me_image+1)*3));   dCpqdR = 0.0_DP
        if(.not.allocated(Tlrdr))           allocate(TLRdr(3, 3, n_comps_f(me_image+1))); TlrdR = 0.0_DP
        if(.not.mbd_vdw_isolated) then
          if(.not.allocated(dCpqdh))          allocate(dCpqdh(3*nat_sl, 3*nat_sl, n_comps_h(me_image+1)*3));   dCpqdh = 0.0_DP
          if(.not.allocated(Tlrdh))           allocate(TLRdh(3, 3, n_comps_h(me_image+1))); TlrdH = 0.0_DP
          end if
        end if
    counter = 1
    do p = 1, nat_sl, 1
      do q = p, nat_sl, 1
        aSqrt = dsqrt(alpha_0_sl(p)*alpha_0_sl(q))
        omMult = omega_scs_sl(p)*omega_scs_sl(q)
        prefactor = omMult*aSqrt
        TLR = collected_tlr(counter, :, :)
        if(do_forces) TLRdr = collected_dtlrdr(counter, :, :, :)
        if(do_forces.and.(.not.mbd_vdw_isolated)) TLRdh = collected_dtlrdh(counter, :, :, :)
        if(vdw_self_consistent) TLRdV = collected_dtlrdV(counter, :, :, :)
        counter = counter + 1
        do i_idx = 1, 3, 1
          do j_idx = i_idx, 3, 1
            i_index = (3*p - 3 + i_idx)
            j_index = (3*q - 3 + j_idx)
            i_index_T = (3*p - 3 + j_idx)
            j_index_T = (3*q - 3 + i_idx)
            if(p.eq.q) then
              if(i_idx.eq.j_idx) then
                Cpq(i_index, j_index) = omega_scs_sl(p)**2.0_DP
                if(vdw_self_consistent) then
                  cnt_v = 1
                  do i_f = 1, nat, 1
                    if(v_cpu_id(i_f).eq.me_image) then
                      s = i_f
                      dCpqdV(i_index, j_index, cnt_v) = 2.0_DP*omega_scs_sl(p)*domegadV_sl(p, s)
                      cnt_v = cnt_v + 1
                      end if ! checks if this MPI rank computes this component
                    end do ! Loop over all possible components
                  end if ! if to check if self consistent derivatives need to be computed
                if(do_forces) then
                  ! dR
                  cnt_f = 1
                  do i_f = 1, 3*nat, 1
                    if(f_cpu_id(i_f).eq.me_image) then
                      call mbdvdw_get_is(i_f, s, i)
                      dCpqdR(i_index,j_index,cnt_f)=2.0_DP*omega_scs_sl(p)*domegadR_sl(p,s,i)
                      cnt_f = cnt_f + 1
                      end if ! checks if this MPI rank computes this component
                    end do ! Loop over all possible components
                  ! dH
                  if(.not.mbd_vdw_isolated) then
                    cnt_h = 1
                    do i_f = 1, 9, 1
                      if(h_cpu_id(i_f).eq.me_image) then
                        call mbdvdw_get_is(i_f, s, i)
                        dCpqdh(i_index,j_index,cnt_h)=2.0_DP*omega_scs_sl(p)*domegadh_sl(p,s,i)
                        cnt_h = cnt_h + 1
                        end if ! checks if this MPI rank computes this component
                      end do ! Loop over all possible components
                    end if
                  end if ! if to check if forces and stresses need to be computed
                end if ! i_idx = j_idx if statement
              end if ! p = q iff statement
            Cpq(i_index, j_index) = Cpq(i_index, j_index) +  prefactor*TLR(i_idx, j_idx)
            Cpq(j_index, i_index) = (Cpq(i_index, j_index))
            Cpq(i_index_T, j_index_T) = Cpq(i_index, j_index)
            Cpq(j_index_T, i_index_T) = Cpq(j_index, i_index)
            if(vdw_self_consistent) then
              cnt_v = 1
              do i_f = 1, nat, 1
                if(v_cpu_id(i_f).eq.me_image) then
                  s = i_f
                  pre1=(omega_scs_sl(p)*domegadV_sl(q,s)+omega_scs_sl(q)*domegadV_sl(p,s))*asqrt
                  pre2=omMult*(alpha_0_sl(q)*dalpha_0dV_sl(p,s)+alpha_0_sl(p)*dalpha_0dV_sl(q,s))
                  pre2=pre2/(2.0_DP*asqrt)
                  der = (pre1+pre2)*TLR(i_idx,j_idx)+prefactor*TLRdV(i_idx,j_idx,cnt_v)
                  dCpqdV(i_index,j_index,cnt_v)=dCpqdV(i_index,j_index,cnt_v)+der
                  dCpqdV(j_index,i_index,cnt_v)=(dCpqdV(i_index,j_index,cnt_v))
                  dCpqdV(i_index_T,j_index_T,cnt_v) = dCpqdV(i_index,j_index,cnt_v)
                  dCpqdV(j_index_T,i_index_T,cnt_v) = dCpqdV(j_index,i_index,cnt_v)
                  cnt_v = cnt_v + 1
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              end if ! if to check if self consistent derivatives need to be computed

            if(do_forces) then
              ! dR
              cnt_f = 1
              do i_f = 1, 3*nat, 1
                if(f_cpu_id(i_f).eq.me_image) then
                  call mbdvdw_get_is(i_f, s, i)
                  pre1=omega_scs_sl(p)*domegadR_sl(q,s,i)+omega_scs_sl(q)*domegadR_sl(p,s,i)
                  pre1=pre1*asqrt
                  pre2=omMult
                  pre2=pre2*(alpha_0_sl(q)*dalpha_0dR_sl(p,s,i)+alpha_0_sl(p)*dalpha_0dR_sl(q,s,i))
                  pre2=pre2/(2.0_DP*asqrt)
                  der = (pre1+pre2)*TLR(i_idx,j_idx)+prefactor*TLRdR(i_idx,j_idx,cnt_f)
                  dCpqdR(i_index,j_index,cnt_f)=dCpqdR(i_index,j_index,cnt_f)+der
                  dCpqdR(j_index,i_index,cnt_f)=(dCpqdR(i_index,j_index,cnt_f))
                  dCpqdR(i_index_T,j_index_T,cnt_f) = dCpqdR(i_index,j_index,cnt_f)
                  dCpqdR(j_index_T,i_index_T,cnt_f) = dCpqdR(j_index,i_index,cnt_f)
                  cnt_f = cnt_f + 1
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              ! dH
              if(.not.mbd_vdw_isolated) then
                cnt_h = 1
                do i_f = 1, 9, 1
                  if(h_cpu_id(i_f).eq.me_image) then
                    call mbdvdw_get_is(i_f, s, i)
                    pre1=omega_scs_sl(p)*domegadh_sl(q,s,i)+omega_scs_sl(q)*domegadh_sl(p,s,i)
                    pre1=pre1*asqrt
                    pre2=omMult
                    pre2=pre2*(alpha_0_sl(q)*dalpha_0dh_sl(p,s,i)+alpha_0_sl(p)*dalpha_0dh_sl(q,s,i))
                    pre2=pre2/(2.0_DP*asqrt)
                    der = (pre1+pre2)*TLR(i_idx,j_idx)+prefactor*TLRdh(i_idx,j_idx,cnt_h)
                    dCpqdh(i_index,j_index,cnt_h)=dCpqdh(i_index,j_index,cnt_h)+der
                    dCpqdh(j_index,i_index,cnt_h)=(dCpqdh(i_index,j_index,cnt_h))
                    dCpqdh(i_index_T,j_index_T,cnt_h) = dCpqdh(i_index,j_index,cnt_h)
                    dCpqdh(j_index_T,i_index_T,cnt_h) = dCpqdh(j_index,i_index,cnt_h)
                    cnt_h = cnt_h + 1
                    end if ! checks if this MPI rank computes this component
                  end do ! Loop over all possible components
                end if
              end if ! if to check if forces and stresses need to be computed
            end do ! j_idx loop
          end do ! i_idx loop
        end do ! q loop
      end do ! q loop

    call stop_clock('mbd_build_c')
    if(allocated(collected_tlr)) deallocate(collected_tlr)
    if(allocated(collected_dtlrdr)) deallocate(collected_dtlrdr)
    if(allocated(collected_dtlrdh)) deallocate(collected_dtlrdh)
    if(allocated(collected_dtlrdV)) deallocate(collected_dtlrdv)

    call stop_clock('mbd_construct_hamiltonian')
    return
    end subroutine mbdvdw_construct_hamiltonian_real

  subroutine mbdvdw_calculate(tauin, rhor, at_x_alat, bg_trans_by_alat)
    implicit none
    ! IO Variables
    ! Atomic density
    real(dp), intent(in) :: rhor(:,:)
    ! Coordinates
    real(dp) :: tauin(3, nat)
    REAL(DP), INTENT(IN), OPTIONAL :: at_x_alat(3,3)
    REAL(DP), INTENT(IN), OPTIONAL :: bg_trans_by_alat(3,3)

    ! Local variables
    real(dp) :: toAng
    integer :: i_freq, p, s, i, i_idx, i_f
    real(dp), dimension(nat) :: alpha_temp
    real(dp), dimension(nat, nat, 3) :: dalpha_tempdR
    real(dp), dimension(nat, 3  , 3) :: dalpha_tempdh
    real(dp), dimension(nat, nat   ) :: dalpha_tempdV
    real(dp), dimension(3,3) :: Htmp
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    logical :: alpha_zero_failed

    call start_clock('mbd')

    if(vdw_debug) mbd_vdw_verbosity = -1

    if (present(at_x_alat)) then
      ! PW case
      h_(:,:) = at_x_alat(:,:)
      call inv3x3_mat(h_, ainv_)
      else
      ! CP case
      h_(:,:) = h_(:,:)
      ainv_(:,:) = ainv_(:,:)
      end if
    call inv3x3_mat(h_, ainv_)

    ! This short circuits the MBD calls if the run is periodic and neither a supercell nor a kpoint grid are specified
    if(MBD_PER_ERROR) goto 10
    EmbdvdW = 0.0_DP
    FmbdVdW = 0.0_DP
    HmbdVdW = 0.0_DP
    Umbdvdw = 0.0_DP
    alpha_zero_failed = .false.
    tau = tauin
    toAng = 0.52917721092_DP

    do_forces = mbd_vdw_forces
    if(mbd_first_step) goto 10
    if((.not.vdw_self_consistent).and.(.not.mbd_conv_elec)) goto 10
    if(vdw_self_consistent.and.mbd_vdw_forces) then
      if(mbd_conv_elec) then
        do_forces = .true.
        else
        do_forces = .false.
        end if
      end if
    if((do_forces.and.mbd_conv_elec).and.(.not.vdw_self_consistent)) then
        !mbd_conv_elec = .false. !! Having this flag false is surpressing output (it gets aliased to mbd_scf_converged)
        !see the file get_total_energy.f90.
        !commenting out the above line does not seem to create any problem, and restores output.  Josh Berryman
 
        mbd_first_step = .true.
      end if

    if(vdw_debug) then
      mbd_conv_elec = .true.
      mbd_first_step = .false.
      end if

    call start_clock('mbd_veff')
    if(.not.vdw_self_consistent) vdw_lscreen = .true.
    call hirshfeld_calculate(tauin, rhor, h_, ainv_, .true.)
    if (do_forces) then
        dveffdr = -1.0_DP*dveffdr
        dveffdh = -1.0_DP*dveffdh
    end if

    call stop_clock('mbd_veff')

    if(do_recip) then
      if(.not.allocated(k_grid)) allocate(k_grid(nk, 3))
      k_grid = make_k_grid(mbd_vdw_kgrid(1), mbd_vdw_kgrid(2), mbd_vdw_kgrid(3), h_)
    else
      if(.not.allocated(k_grid)) allocate(k_grid(1, 3)); k_grid = 0.0_dp
      end if

    if(mbd_vdw_verbosity.ge.0 .and. me_image == root_image) then
      write(use_unit, *) "mbd_vdw_vacuum = ", mbd_vdw_vacuum
      write(use_unit, *) "mbd_vdw_isolated = ", mbd_vdw_isolated
      write(use_unit, *) "do_forces = ", do_forces
      write(use_unit, *) "do_complex = ", do_cmplx
      write(use_unit, *) "do_recip = ", do_recip
      write(use_unit, *) "do_ewald = ", do_ewald, mbd_vdw_ewald
      write(use_unit, *) "supercell cutoff = ", supercell_cutoff
      write(use_unit, *) "maxval cutoff = ", mbd_vdw_econv_thr
      write(use_unit, *) "nk = ", nk
      write(use_unit, *) "kpoint size = ", mbd_vdw_kgrid(1), mbd_vdw_kgrid(2), mbd_vdw_kgrid(3)
      write(use_unit, *) "h_", h_
      end if

    call start_clock('mbd_firstpinit')
    call mbdvdw_para_init()
    call stop_clock('mbd_firstpinit')

    ! Sets up the grid and computes the derivatives with respect to volume

    call mbdvdw_zeroout()
    call start_clock('mbd_fparainit')
    call mbdvdw_para_init_forces(3*nat, f_cpu_id, n_comps_f, max_proc_forces)
    call mbdvdw_para_init_forces(nat, v_cpu_id, n_comps_v, max_proc_sc)
    call mbdvdw_para_init_forces(9,   h_cpu_id, n_comps_h, max_proc_h)
    call stop_clock('mbd_fparainit')
    call start_clock('mbd_pbc')
    call mbdvdw_pbc(tau, h_, ainv_, nat)
    call stop_clock('mbd_pbc')

    sl_mult = (/1.0_dp, 1.0_dp, 1.0_dp/)
    if(allocated(orig_idx)) deallocate(orig_idx)
    allocate(orig_idx(nat))
    do p = 1, nat, 1
      orig_idx(p) = p
      end do

    call print_log('Total dynamic anisotropic polarizability, \alpha_{ij}(iu), (a.u.)')
    call print_log('-------------------------------------------------------------------')
    call print_log('     u            xx         yy         zz         yz         xz         xy')
    call start_clock('mbd_loop')
    do i_freq = 1,npts,1
      if(mbd_vdw_verbosity.ge.0 .and. me_image == root_image) then
          write(use_unit, '(3X, "omega = ", F14.6)') casimir_omega(i_freq)
      end if
      omega_to_print = casimir_omega(i_freq)
      alpha_temp = 0.0_DP
      dalpha_tempdR = 0.0_DP
      dalpha_tempdh = 0.0_DP
      dalpha_tempdV = 0.0_DP
      ! Compute $R_{VdW}$
      call mbdvdw_effqts(casimir_omega(i_freq))
      ! Performs the self consistent screening
      call mbdvdw_SCS()
        ! Computes the contracted alphas
        if (me_image == root_image) then
            write (default_unit, '(2x,"| ",f10.6,6(f11.3))') &
                casimir_omega(i_freq), sum(A_matrix(1::3, 1::3)), &
                sum(A_matrix(2::3, 2::3)), sum(A_matrix(3::3, 3::3)), &
                sum(A_matrix(2::3, 3::3)), sum(A_matrix(1::3, 3::3)), &
                sum(A_matrix(1::3, 2::3))
        end if
        if(i_freq.eq.1) then
          call mbdvdw_calculate_screened_pol(alpha_0, dalpha_0dR, dalpha_0dh, dalpha_0dV)
          if(me_image.eq.root_image) then
            do i_idx=1,nat,1
              if(alpha_0(i_idx).le.0.0 .and. me_image == root_image) then
                write(use_unit, '(3X, "Negative alpha! Aborting this iteration")')
                EmbdvdW = 0.0_DP
                Fmbdvdw = 0.0_DP
                Hmbdvdw = 0.0_DP
                alpha_zero_failed = .true.
                GOTO 20
                end if
              end do
              20 CONTINUE
              do i=1, nproc_image-1,1
                call mpi_send(alpha_zero_failed, 1, mpi_logical,i, 0, intra_image_comm,ierr)
                end do
              if(alpha_zero_failed) GOTO 10
            end if
          else
          alpha_temp = 0.0_DP
          call mbdvdw_calculate_screened_pol(alpha_temp, dalpha_tempdR, dalpha_tempdh, dalpha_tempdV)
          ! Performs the casimir polder integral
          do p=1, nat, 1
            omega_scs(p)=omega_scs(p)+(4.0_DP/pi)*casimir_omega_weight(i_freq)*&
                                    alpha_temp(p)**2.0_DP/(alpha_0(p)**2.0_DP)
            if(vdw_self_consistent) then
              do i_f = 1, nat, 1
                if(v_cpu_id(i_f).eq.me_image) then
                  domegadV(p, i_f) = domegadV(p, i_f) + casimir_omega_weight(i_freq)*(8.0_DP/pi)*&
                                      (alpha_temp(p)*dalpha_tempdV(p, i_f) - 1.0_DP/alpha_0(p)*&
                                      dalpha_0dV(p, i_f)*alpha_temp(p)**2.0_DP)/(alpha_0(p)**2.0)
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              end if ! if to check if self consistent derivatives need to be computed

            if(do_forces) then
              ! dR
              do i_f = 1, 3*nat, 1
                if(f_cpu_id(i_f).eq.me_image) then
                  call mbdvdw_get_is(i_f, s, i)
                  domegadR(p, s, i) =  domegadR(p, s, i) + casimir_omega_weight(i_freq)*(8.0_DP/pi)*&
                          (alpha_temp(p)*dalpha_tempdR(p, s, i) - 1.0_DP/alpha_0(p)*&
                            dalpha_0dR(p, s, i)*alpha_temp(p)**2.0_DP)/(alpha_0(p)**2.0)
                  end if ! checks if this MPI rank computes this component
                end do ! Loop over all possible components
              ! dH
              if(.not.mbd_vdw_isolated) then
                do i_f = 1, 9, 1
                  if(h_cpu_id(i_f).eq.me_image) then
                    call mbdvdw_get_is(i_f, s, i)
                    domegadh(p, s, i) =  domegadh(p, s, i) + casimir_omega_weight(i_freq)*(8.0_DP/pi)*&
                            (alpha_temp(p)*dalpha_tempdh(p, s, i) - 1.0_DP/alpha_0(p)*&
                              dalpha_0dh(p, s, i)*alpha_temp(p)**2.0_DP)/(alpha_0(p)**2.0)
                    end if ! checks if this MPI rank computes this component
                  end do ! Loop over all possible components
                end if
              end if ! if to check if forces and stresses need to be computed
            end do
          end if
      if((me_image.ne.root_image).and.(i_freq.eq.1)) then
        call mpi_recv(alpha_zero_failed,& ! Buffer
                1,& ! Count
                mpi_logical,& ! Type
                0,& ! Source
                0,& ! Tag
                intra_image_comm,& ! Communicator
                rstatus,& ! Status var
                ierr) ! Error flag
        if(alpha_zero_failed) then
          EmbdvdW = 0.0_DP
          Fmbdvdw = 0.0_DP
          Hmbdvdw = 0.0_DP
          GOTO 10
          end if
        end if
      end do
      if (me_image == root_image) then
          call stop_clock('mbd_loop')
          call print_log('--------------------------------------------------------------')
          call print_log('Partitioned atomic C6 coefficients and polarizabilities (a.u.)')
          call print_log('--------------------------------------------------------------')
          do i_idx = 1, nat
              write (default_unit, '(2x,A,I4,2x,A,f12.6,6x,f12.6)') &
                  "| ATOM", i_idx, trim(atm(ityp(i_idx))), &
                  3.d0/4*alpha_0(i_idx)**2*omega_scs(i_idx), alpha_0(i_idx)
          end do
          call print_log('--------------------------------------------------------------')
          call print_log('Computing MBD@rsSCS energy...')
      end if
      ! Compute the noninteracting energy.
      ! This is done by simply summing over the effective frequencies.
      ! \sum_p \bar{\omega}_p
      ! Where \bar{\omega}_p is given by
      ! \bar{\omega}_p = \frac{4}{\pi} \frac{ \int_0^{\infty} \left[ \alpha_p(i \omega) \right]^2 d\omega  }{ \left[ \alpha_p^0 \right]^2 }.
    call mbdvdw_noninteracting_energy(non_interacting_energy, nonint_FmbdVdW, nonint_HmbdVdW, nonint_Uprefactor)
    call mbdvdw_cleanup_postscs()

    if(vdw_self_consistent) then
      call mp_sum(dalpha_0dV, intra_image_comm)
      call mp_sum(domegadV, intra_image_comm)
      end if

    if(do_forces) then
      call mp_sum(dalpha_0dR, intra_image_comm)
      call mp_sum(domegadR, intra_image_comm)
      call mp_sum(dalpha_0dh, intra_image_comm)
      call mp_sum(domegadh, intra_image_comm)
      end if
    ! This method constructs C so that we can find its eigenvalues and sum them to find the interacting
    ! parts of the energy. This is done by constructing
    ! \mathcal{C}^{MBD}_{pq} = \delta_{pq} \bar{\omega}^2_p + \left(1 - \delta_pq \right) \bar{\omega}_p \bar{\omega}_q \sqrt{\bar{\alpha}_p^0 \bar{\alpha}_q^0} T_{pq}
    if(do_cmplx) then
      if(allocated(f_cpu_id)) deallocate(f_cpu_id)
      if(allocated(h_cpu_id)) deallocate(h_cpu_id)
      if(allocated(v_cpu_id)) deallocate(v_cpu_id)
      if(.not.allocated(f_cpu_id))        allocate(f_cpu_id(3*nat*nk)); f_cpu_id = 0
      if(.not.allocated(h_cpu_id))        allocate(h_cpu_id(9*nk)); h_cpu_id = 0
      if(.not.allocated(v_cpu_id))        allocate(v_cpu_id(nat*nk)); v_cpu_id = 0
      if(.not.allocated(k_cpu_id))        allocate(k_cpu_id(nk)); k_cpu_id = 0
      if(.not.allocated(n_comps_k))       allocate(n_comps_k(nproc_image)); n_comps_k = 0
      call mbdvdw_para_init_forces(nk, k_cpu_id, n_comps_k, max_proc_k)
      call mbdvdw_para_init_forces(3*nat*nk, f_cpu_id, n_comps_f, max_proc_forces)
      call mbdvdw_para_init_forces(nat*nk, v_cpu_id, n_comps_v, max_proc_sc)
      call mbdvdw_para_init_forces(9*nk,   h_cpu_id, n_comps_h, max_proc_h)

      call mbdvdw_construct_hamiltonian_complex()
      else
      call mbdvdw_construct_hamiltonian_real()
      end if

    ! Computes the interacting energy by diagonalizing the MBD hamiltonian
    if(do_cmplx) then
      call mbdvdw_interacting_energy_complex(interacting_energy, int_FmbdVdW, int_HmbdVdW, int_Uprefactor)
      else
      call mbdvdw_interacting_energy(interacting_energy, int_FmbdVdW, int_HmbdVdW, int_Uprefactor)
      end if
    if (interacting_energy == 0.d0) then
        alpha_zero_failed = .true.
        goto 10
    end if
    if(vdw_self_consistent) then
      Uprefactor = 0.5_DP*int_Uprefactor - 1.5_DP*nonint_Uprefactor
      CALL hirshfeld_wfforce(Uprefactor, UmbdvdW) !HK/RAD/TM: turn on/off self consistency
      end if
    if(me_image.eq.root_image) then
      EmbdvdW = 0.5_DP*interacting_energy - 1.5_DP*non_interacting_energy
      if(mbd_vdw_verbosity.ge.0) then
        write(use_unit, *) "nonint_energy = ", non_interacting_energy
        write(use_unit, *) "int_energy = ", interacting_energy
        write(use_unit, *) "mbdvdw energy = ", Embdvdw
        end if

      if(do_forces) then
        FmbdVdW = 0.5_DP*int_FmbdVdW - 1.5_DP*nonint_FmbdVdW
        HmbdVdW = 0.5_DP*int_HmbdVdW - 1.5_DP*nonint_HmbdVdW
        HmbdVdW = -1.0_DP*HmbdvdW
        if(mbd_vdw_isolated) HmbdvdW = 0.0_dp
        if(mbd_vdw_verbosity.ge.0) then
          write(use_unit, '(3X, "FmbdVdW")')
          write(use_unit, '(3X, "DX       DY       DZ")')
          do s = 1, nat, 1
            write(use_unit, '(3X, F14.10, F14.10, F14.10)') FmbdVdW(s, 1), FmbdVdW(s, 2), FmbdVdW(s, 3)
            end do

          Htmp = Hmbdvdw
          write(use_unit, '(3X, "HmbdVdW")')
          write(use_unit, '(3X, "DX       DY       DZ")')
          do s = 1, 3, 1
            write(use_unit, '(3X, F14.10, F14.10, F14.10)') Htmp(s, 1), Htmp(s, 2), Htmp(s, 3)
            end do

          Htmp = -2.0_DP*MATMUL( HmbdvdW, transpose(h_) )/nu
          write(use_unit, '(3X, "MBD Stresses [Ry]")')
          write(use_unit, '(3X, "DX       DY       DZ")')
          do s = 1, 3, 1
            write(use_unit, '(3X, F14.10, F14.10, F14.10)') Htmp(s, 1), Htmp(s, 2), Htmp(s, 3)
            end do

          end if
        end if
      end if

    10 CONTINUE

    if(alpha_zero_failed) then
      Embdvdw = 0.0_DP
      Fmbdvdw = 0.0_DP
      Hmbdvdw = 0.0_DP
      Umbdvdw = 0.0_DP
      end if

    call MPI_Bcast(EmbdvdW, & ! buffer
             1,       & ! count
             mpi_double_precision, & ! MPI_Datatype datatype
             root_image, & !root
             intra_image_comm, & ! Communicator
             ierr)

    call MPI_Bcast(FmbdvdW, & ! buffer
             3*nat,   & ! count
             mpi_double_precision, & ! MPI_Datatype datatype
             root_image, & !root
             intra_image_comm, & ! Communicator
             ierr)

    call MPI_Bcast(HmbdvdW, & ! buffer
             3*3,     & ! count
             mpi_double_precision, & ! MPI_Datatype datatype
             root_image, & !root
             intra_image_comm, & ! Communicator
             ierr)

    call mbdvdw_cleanup()
    call hirshfeld_cleanup()
    call stop_clock('mbd')
    call print_log('--------------------------------------------------------------')

    end subroutine mbdvdw_calculate

end module mbdvdw_module
