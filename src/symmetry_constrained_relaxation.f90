module symmetry_constrained_relaxation
  use dimensions  ! contains all input needed
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use fparser, only: initf, parsef, evalf
  use matrix_inversion, only: invmat_real_lapack  ! MOL for symm-constr. relaxation
  implicit none
  real*8, dimension(:), allocatable :: SCR_params_val  ! parameter values for parsing
  real*8, dimension(:,:), allocatable :: SCR_A_coords
  real*8, dimension(:,:), allocatable :: SCR_A_coords_inv
  real*8, dimension(:), allocatable   :: SCR_B_coords    ! vector for lattice traform
  real*8, dimension(:), allocatable   :: SCR_B_coords_original ! vector for lattice trafom (for writing files with same symmetry)
  real*8, dimension(:,:), allocatable :: SCR_A_lv
  real*8, dimension(:,:), allocatable :: SCR_A_lv_inv
  real*8, dimension(:), allocatable   :: SCR_B_lv        ! vector for coords traform
  real*8, dimension(:,:), allocatable :: SCR_cart2frac_jacobian
  real*8, dimension(:,:), allocatable :: SCR_force_jacobian_inv
  real*8, dimension(:,:), allocatable :: SCR_force_jacobian_sym
  real*8, dimension(:,:), allocatable :: SCR_force_jacobian

contains

  subroutine allocate_SCR()
    call aims_allocate(SCR_params_val, SCR_n_params, "SCR_params_val")
    call aims_allocate(SCR_A_coords,3*n_atoms,SCR_n_params_coords, "SCR_A_coords")
    call aims_allocate(SCR_A_coords_inv,SCR_n_params_coords,3*n_atoms, "SCR_A_coords_inv")
    call aims_allocate(SCR_A_lv,3*n_periodic,SCR_n_params_lv,"SCR_A_lv")
    call aims_allocate(SCR_A_lv_inv,SCR_n_params_lv,3*n_periodic,"SCR_A_lv_inv")
    call aims_allocate(SCR_B_coords,3*n_atoms,"SCR_B_coords")
    call aims_allocate(SCR_B_coords_original,3*n_atoms,"SCR_B_coords_original")
    call aims_allocate(SCR_B_lv,3*n_periodic,"SCR_B_lv")
    call aims_allocate(SCR_cart2frac_jacobian, 3*n_atoms, 3*n_atoms, "SCR_cart2frac_jacobian")
    call aims_allocate(SCR_force_jacobian, 3*n_atoms, SCR_n_params_coords, "SCR_force_jacobian")
    call aims_allocate(SCR_force_jacobian_sym, SCR_n_params_coords, SCR_n_params_coords, "SCR_force_jacobian_sym")
    call aims_allocate(SCR_force_jacobian_inv, SCR_n_params_coords, 3*n_atoms, "SCR_force_jacobian_inv")
  end subroutine allocate_SCR

  subroutine deallocate_SCR()
    call aims_deallocate(SCR_params_val, "SCR_params_val")
    call aims_deallocate(SCR_A_coords,"SCR_A_coords")
    call aims_deallocate(SCR_A_coords_inv, "SCR_A_coords_inv")
    call aims_deallocate(SCR_A_lv,"SCR_A_lv")
    call aims_deallocate(SCR_A_lv_inv,"SCR_A_lv_inv")
    call aims_deallocate(SCR_B_coords,"SCR_B_coords")
    call aims_deallocate(SCR_B_coords_original,"SCR_B_coords_original")
    call aims_deallocate(SCR_B_lv,"SCR_B_lv")
    call aims_deallocate(SCR_cart2frac_jacobian, "SCR_cart2frac_jacobian")
    call aims_deallocate(SCR_force_jacobian_inv, "SCR_force_jacobian_inv")
    call aims_deallocate(SCR_force_jacobian_sym, "SCR_force_jacobian_sym")
    call aims_deallocate(SCR_force_jacobian, "SCR_force_jacobian")
  end subroutine deallocate_SCR


  subroutine SCR_create_coords_matrices()
    real*8, dimension(SCR_n_params_coords,SCR_n_params_coords)    :: SCR_A_coords_sym
    real*8  :: SCR_res
    integer :: i, j, k, l

    call initf(3*n_atoms)
    k = 0
    do i=1,n_atoms
      do j=1,3
        k = k + 1
        call parsef (k, SCR_coord_str(j,i), SCR_params)
      end do
    end do
    !translational part
    SCR_params_val = 0.
    k = 0
    do i=1,n_atoms
      do j=1,3
        k = k + 1
        SCR_res = evalf(k, SCR_params_val)
        SCR_B_coords(k) = SCR_res
        SCR_B_coords_original(k) = SCR_res
      end do
    end do
    ! variable parts
    do l=SCR_n_params_lv+1,SCR_n_params
      SCR_params_val = 0.
      SCR_params_val(l) = 1.
      k = 0
      do i=1,n_atoms
        do j=1,3
          k = k + 1
          SCR_res = evalf(k, SCR_params_val)
          SCR_A_coords(k,l-SCR_n_params_lv) = SCR_res - SCR_B_coords(k)
        end do
      end do
    end do

    ! Invert matrix
    if (SCR_n_params_coords > 0) then
      SCR_A_coords_sym = matmul( transpose(SCR_A_coords), SCR_A_coords ) ! ^(-1) * transpose(SCR_A_coords)
      call invmat_real_lapack(SCR_A_coords_sym)
      SCR_A_coords_inv = matmul( SCR_A_coords_sym, transpose(SCR_A_coords) )
    else
      SCR_A_coords = 0
      SCR_A_coords_inv = 0
    endif
  end subroutine SCR_create_coords_matrices

  subroutine SCR_create_lv_matrices()
    use constants, only: bohr
    real*8, dimension(SCR_n_params_lv,SCR_n_params_lv)           :: SCR_A_lv_sym
    real*8  :: SCR_res
    integer :: i, j, k, l

    ! use fparser to parse functions from strings
    call initf(3*n_periodic)
    k = 0
    do i=1,n_periodic
      do j=1,3
        k = k + 1
        call parsef (k, SCR_lv_str(j,i), SCR_params)
      end do
    end do
    !translational part
    SCR_params_val = 0.
    k = 0
    do i=1,n_periodic
      do j=1,3
        k = k + 1
        SCR_res = evalf(k, SCR_params_val)
        SCR_B_lv(k) = SCR_res / bohr
      end do
    end do
    if ((SCR_n_params_lv > 0) .and. (relax_unit_cell .ne. 0)) then
      ! variable parts
      do l=1,SCR_n_params_lv
        SCR_params_val = 0.
        SCR_params_val(l) = 1.
        k = 0
          do i=1,n_periodic
            do j=1,3
              k = k + 1
              SCR_res = evalf(k, SCR_params_val)
              SCR_A_lv(k,l) = SCR_res - SCR_B_lv(k)
            end do
          end do
      end do
    end if

    ! Invert matrix
    if ((SCR_n_params_lv > 0) .and. (relax_unit_cell .ne. 0)) then
        SCR_A_lv_sym     = matmul( transpose(SCR_A_lv), SCR_A_lv ) ! ^(-1) * transpose(SCR_A_lv)
        call invmat_real_lapack(SCR_A_lv_sym)
        SCR_A_lv_inv     = matmul( SCR_A_lv_sym, transpose(SCR_A_lv) )
    else
        SCR_A_lv = 0
        SCR_A_lv_inv = 0
    end if
  end subroutine SCR_create_lv_matrices

  subroutine SCR_generate_force_Jacobian(lv)
      real*8, dimension(3,n_periodic), intent(in) :: lv
      integer :: i_atom

      SCR_force_jacobian = 0.d0
      SCR_cart2frac_jacobian = 0.d0
      SCR_force_jacobian_sym = 0.d0
      SCR_force_jacobian_inv = 0.d0

      if (SCR_n_params > 0) then
          do i_atom = 0, 3*(n_atoms-1),3
              SCR_cart2frac_jacobian(i_atom+1:i_atom+3, i_atom+1:i_atom+3) = lv
          end do

          ! Get the full Jacobian
          SCR_force_jacobian = matmul(SCR_cart2frac_jacobian, SCR_A_coords)
          SCR_force_jacobian_sym = matmul( transpose(SCR_force_jacobian), SCR_force_jacobian ) ! ^(-1) * transpose(SCR_A_coords)

          call invmat_real_lapack(SCR_force_jacobian_sym)
          SCR_force_jacobian_inv = matmul( SCR_force_jacobian_sym, transpose(SCR_force_jacobian) )
      end if

  end subroutine

  subroutine SCR_transform(coords, coords_sr)
    ! coords need to be fractional coordinates
    real*8, dimension(3,n_atoms), intent(in)  :: coords
    real*8, dimension(3,n_atoms), intent(out) :: coords_sr
    real*8, dimension(3*n_atoms)           :: coords_1d
    real*8, dimension(SCR_n_params_coords) :: coords_sr_1d
    integer :: i
    coords_1d   = pack(coords, mask=.true.)
    coords_1d   = coords_1d - SCR_B_coords
    coords_sr_1d = matmul(SCR_A_coords_inv, coords_1d)
    coords_sr    = 0.0d0
    coords_sr    = reshape(coords_sr_1d, (/3,n_atoms/), (/0.d0/))
  end subroutine SCR_transform

  subroutine SCR_backtransform(coords_sr, coords, to_file)
    logical, intent(in) :: to_file
    real*8, dimension(3,n_atoms), intent(in) :: coords_sr
    real*8, dimension(3,n_atoms), intent(out) :: coords
    real*8, dimension(3*n_atoms)            :: coords_1d, length
    real*8, dimension(SCR_n_params_coords)  :: coords_sr_1d
    integer :: i

    coords_1d   = pack(coords_sr, mask=.true.)
    coords_sr_1d = coords_1d(1:SCR_n_params_coords)
    coords_1d      = matmul( SCR_A_coords, coords_sr_1d )
    if (to_file) then
      coords_1d      = coords_1d + SCR_B_coords_original
    else
      coords_1d      = coords_1d + SCR_B_coords
      ! map into center cell
      length = coords_1d - 1D-8
      SCR_B_coords = SCR_B_coords - nint(length)
      coords_1d = length - nint(length) + 1D-8
    end if
    coords  = 0.d0
    coords  = reshape(coords_1d, (/3,n_atoms/), (/0.d0/))
  end subroutine SCR_backtransform

  subroutine SCR_symmetrize(coords)
    real*8, dimension(3,n_atoms), intent(inout) :: coords
    real*8, dimension(3,n_atoms)  :: coords_sr
    call SCR_transform(coords, coords_sr)
    call SCR_backtransform(coords_sr, coords, .false.)
  end subroutine SCR_symmetrize

  subroutine SCR_transform_lv(lv, lv_sr)
    ! coords need to be fractional coordinates
    real*8, dimension(3,n_periodic), intent(in)  :: lv
    real*8, dimension(3,n_periodic), intent(out) :: lv_sr
    real*8, dimension(3*n_periodic)    :: lv_1d
    real*8, dimension(SCR_n_params_lv) :: lv_sr_1d
    integer :: i

    lv_1d   = pack(lv, mask=.true.)
    lv_1d   = lv_1d - SCR_B_lv
    lv_sr_1d = matmul( SCR_A_lv_inv, lv_1d )
    lv_sr    = 0.0d0
    lv_sr    = reshape(lv_sr_1d, (/3,n_periodic/), (/0.d0/))
  end subroutine SCR_transform_lv

  subroutine SCR_backtransform_lv(lv_sr, lv)
    real*8, dimension(3,n_periodic), intent(in)  :: lv_sr
    real*8, dimension(3,n_periodic), intent(out) :: lv
    real*8, dimension(3*n_periodic)    :: lv_1d
    real*8, dimension(SCR_n_params_lv) :: lv_sr_1d

    lv_1d   = pack(lv_sr, mask=.true.)
    lv_sr_1d = lv_1d(1:SCR_n_params_lv)
    lv_1d      = matmul( SCR_A_lv, lv_sr_1d )
    lv_1d      = lv_1d + SCR_B_lv
    lv  = 0.d0
    lv  = reshape(lv_1d, (/3,n_periodic/), (/0.d0/))
  end subroutine SCR_backtransform_lv

  subroutine SCR_symmetrize_lv(lv)
    real*8, dimension(3,n_periodic), intent(inout) :: lv
    real*8, dimension(3,n_periodic)  :: lv_sr
    call SCR_transform_lv(lv, lv_sr)
    call SCR_backtransform_lv(lv_sr, lv)
  end subroutine SCR_symmetrize_lv

  ! Transform forces on fractional atomic positions to forces on atom parameters
  subroutine SCR_transform_forces(forces, forces_sr)
    real*8, dimension(3,n_atoms), intent(in)   :: forces
    real*8, dimension(3*n_atoms)               :: forces_1d
    real*8, dimension(SCR_n_params_coords)     :: forces_sr_1d
    real*8, dimension(3,n_atoms), intent(out)  :: forces_sr
    integer :: i

    forces_1d    = pack(forces, mask=.true.)
    forces_sr_1d = matmul( transpose(SCR_force_jacobian), forces_1d)
    forces_sr    = 0.d0
    forces_sr    = reshape(forces_sr_1d, (/3,n_atoms/), (/0.d0/))
  end subroutine SCR_transform_forces

  !Backtransform forces on atom parameters to forces on fractional atomic positions
  subroutine SCR_backtransform_forces(forces_sr, forces)
    real*8, dimension(3,n_atoms), intent(in)   :: forces_sr
    real*8, dimension(3*n_atoms)               :: forces_1d
    real*8, dimension(SCR_n_params_coords)     :: forces_sr_1d
    real*8, dimension(3,n_atoms), intent(out)  :: forces
    integer :: i

    forces_1d   = pack(forces_sr, mask=.true.)
    forces_sr_1d = forces_1d(1:SCR_n_params_coords)
    forces_1d    = matmul( transpose(SCR_force_jacobian_inv), forces_sr_1d)
    forces       = 0.d0
    forces       = reshape(forces_1d, (/3,n_atoms/), (/0.d0/))
  end subroutine SCR_backtransform_forces

  ! Symmetrize atomic forces
  subroutine SCR_symmetrize_forces(forces)
    real*8, dimension(3,n_atoms), intent(inout) :: forces
    real*8, dimension(3,n_atoms) :: forces_sr
    call SCR_transform_forces(forces, forces_sr)
    call SCR_backtransform_forces(forces_sr, forces)
  end subroutine SCR_symmetrize_forces

  ! Transform forces on lattice vectors to forces on lattice vector parameters
  subroutine SCR_transform_forces_lv(forces_lv, forces_lv_sr)
    real*8, dimension(3,n_periodic), intent(in)   :: forces_lv
    real*8, dimension(3*n_periodic)               :: forces_lv_1d
    real*8, dimension(SCR_n_params_lv)     :: forces_lv_sr_1d
    real*8, dimension(3,n_periodic), intent(out)  :: forces_lv_sr

    forces_lv_1d    = pack(forces_lv, mask=.true.)
    forces_lv_sr_1d = matmul( transpose(SCR_A_lv), forces_lv_1d)
    forces_lv_sr    = 0.d0
    forces_lv_sr    = reshape(forces_lv_sr_1d, (/3,n_periodic/), (/0.d0/))
  end subroutine SCR_transform_forces_lv

  ! Backtransform forces on lattice vector parameters to forces on lattice vectors
  subroutine SCR_backtransform_forces_lv(forces_lv_sr, forces_lv)
    real*8, dimension(3,n_periodic), intent(in)   :: forces_lv_sr
    real*8, dimension(3*n_periodic)               :: forces_lv_1d
    real*8, dimension(SCR_n_params_lv)     :: forces_lv_sr_1d
    real*8, dimension(3,n_periodic), intent(out)  :: forces_lv

    forces_lv_1d   = pack(forces_lv_sr, mask=.true.)
    forces_lv_sr_1d = forces_lv_1d(1:SCR_n_params_lv)
    forces_lv_1d    = matmul( transpose(SCR_A_lv_inv), forces_lv_sr_1d)
    forces_lv       = 0.d0
    forces_lv       = reshape(forces_lv_1d, (/3,n_periodic/), (/0.d0/))
  end subroutine SCR_backtransform_forces_lv

  ! Symmetrize forces on lattice vectors
  subroutine SCR_symmetrize_forces_lv(forces_lv)
    real*8, dimension(3,n_periodic), intent(inout) :: forces_lv
    real*8, dimension(3,n_periodic) :: forces_lv_sr
    call SCR_transform_forces_lv(forces_lv, forces_lv_sr)
    call SCR_backtransform_forces_lv(forces_lv_sr, forces_lv)
  end subroutine SCR_symmetrize_forces_lv

  subroutine SCR_reduce_hessian(full_hessian, hessian_sr, lv)
    use bravais, only: get_cell_volume
    use distributed_hessian, only: nrow_hess, ncol_hess
    real*8, dimension(nrow_hess, ncol_hess), intent(in) :: full_hessian
    real*8, dimension(3, 3), intent(in) :: lv
    real*8, dimension(3*n_atoms,3*n_atoms) :: hessian_coords
    real*8, dimension(3*n_periodic,3*n_periodic) :: hessian_lv
    real*8, dimension(SCR_n_params_coords,3*n_atoms) :: hessian_trans_coords
    real*8, dimension(SCR_n_params_lv,3*n_periodic) :: hessian_trans_lv
    real*8, dimension(SCR_n_params_coords,SCR_n_params_coords) :: hessian_sr_coords
    real*8, dimension(SCR_n_params_lv,SCR_n_params_lv) :: hessian_sr_lv
    real*8 :: SCR_scaling

    real*8, dimension(SCR_n_params,SCR_n_params), intent(out) :: hessian_sr
    ! Separate the atomic and lattice components of the Hessian into separate matrices
    ! Cross coupling is ignored by the symmetry constraints
    ! Individually transform atomic and lattice parts of the Hessian and store in hessian_sr
    hessian_sr = 0.d0
    hessian_sr_lv = 0.d0
    hessian_sr_coords = 0.d0
    if (SCR_n_params_coords > 0) then
      ! TARP: Transform Hessian to parameter Space
      ! Get the Jacobian for Cartesian -> Fractional * Fractional -> parameter
      hessian_coords = full_hessian(1:3*n_atoms, 1:3*n_atoms)

      call get_cell_volume(lv, SCR_scaling)
      SCR_scaling = SCR_scaling**(1./3.)

      call SCR_generate_force_Jacobian(lv)

      hessian_trans_coords = matmul( transpose(SCR_force_jacobian), hessian_coords)
      hessian_sr_coords = matmul(hessian_trans_coords, SCR_force_jacobian)
      hessian_sr(1:SCR_n_params_coords, 1:SCR_n_params_coords) = hessian_sr_coords / SCR_scaling**2.0
    end if

    if (SCR_n_params_lv > 0) then
      hessian_lv = full_hessian(3*n_atoms+1:, 3*n_atoms+1:)
      hessian_trans_lv = matmul(transpose(SCR_A_lv), hessian_lv)
      hessian_sr_lv = matmul(hessian_trans_lv, SCR_A_lv)
      hessian_sr(SCR_n_params_coords+1:, SCR_n_params_coords+1:) = hessian_sr_lv(:,:)
    end if

  end subroutine SCR_reduce_hessian

end module symmetry_constrained_relaxation
