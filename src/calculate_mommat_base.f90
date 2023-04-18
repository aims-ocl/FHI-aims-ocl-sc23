!****h* FHI-aims/calculate_mommat
!  NAME
!    calculate_mommat 
!  SYNOPSIS

module calculate_mommat_base
!  PURPOSE
!  This module contains all routines related to the evaluation of the Momentum-
!  matarix-element as postprocessing
!
!  AUTHOR
!  HISTORY
!    Development version, FHI-aims (2010).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE
  use localorb_io
  use mpi_tasks
  implicit none
  real*8, dimension(:,:),allocatable :: work_ovl_mom
  real*8, dimension(:), allocatable :: mommat_full_oned_up
  real*8, dimension(:), allocatable :: mommat_full_oned_low
  real*8, dimension(:), allocatable :: mommat_full_oned_two_up
  real*8, dimension(:), allocatable :: mommat_full_oned_two_low
  real*8, dimension(:), allocatable :: mommat_full_w_up
  real*8, dimension(:), allocatable :: mommat_full_w_low
  complex*16, dimension(:), allocatable :: mommat_full_w_complex_up
  complex*16, dimension(:), allocatable :: mommat_full_w_complex_low
  complex*16, dimension(:), allocatable :: moment_one
  complex*16, dimension(:), allocatable :: moment_two
  complex*16, dimension(:), allocatable :: moment_three

 contains

  !    o mommat_full_oned_up/low upper/lower triangle of Momentummatrix (one 
  !      component x,y or z) for all cells in atom centered basis
  !    o mommat_full_oned_up/low upper/lower triangle of Momentummatrix (on 
  !      component x,y or z) for all cells in atom centered basis, use to store 
  !      different component
  !    o mommat_full_w_up/low -- upper/lower triangle of Momentummatrix (one 
  !      component x,y or z) fourier transformed in atom centered basis, real
  !    o mommat_full_w_up/low -- upper/lower triangle of Momentummatrix (one 
  !      component x,y or z) fourier transformed in atom centered basis, complex
  !    o moment_one/two/thre -- Momentummatrix (one component x,y or z) 
  !      fourier transformed (single k-point) in KS-basis


! Allocate functions









subroutine allocate_mommat
!  PURPOSE
!    allocation of momentum matrix
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  integer :: info
    if (.not.allocated(mommat_full_oned_up)) then
      call aims_allocate( mommat_full_oned_up, n_hamiltonian_matrix_size,   "mommat_full_oned_up" )
    end if
    if (.not.allocated(mommat_full_oned_low)) then
      call aims_allocate( mommat_full_oned_low, n_hamiltonian_matrix_size, "mommat_full_oned_low" )
    end if
end subroutine allocate_mommat

subroutine allocate_mommat_two
!  PURPOSE
!    allocation of momentum matrix (mixed components)
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  integer :: info
  if (.not.allocated(mommat_full_oned_two_up)) then
     call aims_allocate( mommat_full_oned_two_up, n_hamiltonian_matrix_size,   "mommat_full_oned_two_up" )
  end if
  if (.not.allocated(mommat_full_oned_two_low)) then
     call aims_allocate( mommat_full_oned_two_low, n_hamiltonian_matrix_size, "mommat_full_oned_two_low" )
  end if
end subroutine allocate_mommat_two

subroutine allocate_mommat_k
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in atom centered basis and in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  integer :: info
  if(packed_matrix_format==PM_none)then
     call aims_allocate( work_ovl_mom, n_centers_basis_I, n_centers_basis_I,               "work_ovl_mom" )
  else
     ! ONLY dummy allocation - no later code should ever touch these arrays
     ! This is the normal case for pretty much every run, non-packed storage for
     ! periodic systems should never be used.
     call aims_allocate( work_ovl_mom,1, 1,                                                "work_ovl_mom" )
  end if
  if (real_eigenvectors) then
      if (.not.allocated( mommat_full_w_up)) then
        call aims_allocate( mommat_full_w_up , (n_basis+1)*n_basis/2,                  "mommat_full_w_up" )
      end if
      if (.not.allocated( mommat_full_w_low)) then
        call aims_allocate( mommat_full_w_low, (n_basis+1)*n_basis/2,                 "mommat_full_w_low" )
      end if
     !ALLOCATE WITH 1 ELEMENT
      if (.not.allocated( mommat_full_w_complex_up)) then
        call aims_allocate( mommat_full_w_complex_up, 1,                       "mommat_full_w_complex_up" )
      end if
      if (.not.allocated( mommat_full_w_complex_low)) then
        call aims_allocate( mommat_full_w_complex_low, 1,                     "mommat_full_w_complex_low" )
      end if
  else
      if (.not.allocated( mommat_full_w_complex_up)) then
        call aims_allocate( mommat_full_w_complex_up, (n_basis+1)*n_basis/2,   "mommat_full_w_complex_up" )
      end if
      if (.not.allocated( mommat_full_w_complex_low)) then
        call aims_allocate( mommat_full_w_complex_low, (n_basis+1)*n_basis/2, "mommat_full_w_complex_low" )
      end if
      !ALLOCATE WITH 1 ELEMENT
      if (.not.allocated( mommat_full_w_up)) then
        call aims_allocate( mommat_full_w_up, 1,                                       "mommat_full_w_up" )
      end if
      if (.not.allocated( mommat_full_w_low)) then
        call aims_allocate( mommat_full_w_low, 1,                                     "mommat_full_w_low" )
      end if
  endif
end subroutine allocate_mommat_k

subroutine allocate_moment_one_p1_SOC(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    n_state_min_in -- Minimum state considered
!    n_state_max_in -- Maximum state considered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer :: info
  if (.not.allocated( moment_one)) then
    call aims_allocate( moment_one, &
           (((n_state_max_in-n_state_min_in+1)+1)* (n_state_max_in-n_state_min_in+1)/2)*2,  "moment_one" )
  end if
end subroutine allocate_moment_one_p1_SOC

subroutine allocate_moment_two_p1_SOC(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    n_state_min_in -- Minimum state considered
!    n_state_max_in -- Maximum state considered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer :: info
    if (.not.allocated( moment_two)) then
      call aims_allocate( moment_two, (((n_state_max_in-n_state_min_in+1)+1)*&
        (n_state_max_in-n_state_min_in+1)/2)*2,                                              "moment_two" )
    end if
end subroutine allocate_moment_two_p1_SOC

subroutine allocate_moment_three_p1_SOC(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    n_state_min_in -- Minimum state considered
!    n_state_max_in -- Maximum state considered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer :: info
    if (.not.allocated( moment_three)) then
      call aims_allocate( moment_three, (((n_state_max_in-n_state_min_in+1)+1)*&
        (n_state_max_in-n_state_min_in+1)/2)*2,                                            "moment_three" )
    end if
end subroutine allocate_moment_three_p1_SOC


subroutine allocate_moment_one(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    n_state_min_in -- Minimum state considered
!    n_state_max_in -- Maximum state considered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer :: info
  if (.not.allocated( moment_one)) then
    call aims_allocate( moment_one, (((n_state_max_in-n_state_min_in+1)+1)*&
        (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2),                        "moment_one" )
  end if
end subroutine allocate_moment_one

subroutine allocate_moment_two(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    n_state_min_in -- Minimum state considered
!    n_state_max_in -- Maximum state considered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer :: info
    if (.not.allocated( moment_two)) then
      call aims_allocate( moment_two, (((n_state_max_in-n_state_min_in+1)+1)*&
        (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2),                        "moment_two" )
    end if
end subroutine allocate_moment_two

subroutine allocate_moment_three(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of momentum matrix, integrated, Fourier transformed, 
!    in KS-basis
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use aims_memory_tracking, only : aims_allocate
!  ARGUMENTS
!  INPUTS
!    n_state_min_in -- Minimum state considered
!    n_state_max_in -- Maximum state considered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  integer :: info
    if (.not.allocated( moment_three)) then
      call aims_allocate( moment_three, (((n_state_max_in-n_state_min_in+1)+1)*&
        (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2),                      "moment_three" )
    end if
end subroutine allocate_moment_three
subroutine clean_mommat
!  PURPOSE
!    deallocation of arrays
!  USES
   use runtime_choices
   use aims_memory_tracking, only : aims_deallocate
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  !write(use_unit,*) "TZ: state 0 "
  if (allocated(work_ovl_mom))  call aims_deallocate( work_ovl_mom,            "work_ovl_mom" )
  if (allocated(moment_one))    call aims_deallocate( moment_one,                "moment_one" )
  !write(use_unit,*) "TZ: state 0.1"
  if (allocated(moment_two))    call aims_deallocate( moment_two,                "moment_two" )
  if (allocated(moment_three))  call aims_deallocate( moment_three,            "moment_three" )
  !write(use_unit,*) "TZ:state 1 "
  if (allocated(mommat_full_w_complex_up))  &
           call aims_deallocate( mommat_full_w_complex_up,            "mommat_full_w_complex" )
  if (allocated(mommat_full_w_up)) call aims_deallocate( mommat_full_w_up, "mommat_full_w_up" )
  !write(use_unit,*) "TZ:state 2 "
  if (allocated(mommat_full_w_complex_low))  &
           call aims_deallocate( mommat_full_w_complex_low,       "mommat_full_w_complex_low" )
  if (allocated(mommat_full_w_low)) &
           call aims_deallocate( mommat_full_w_low,                       "mommat_full_w_low" )
end subroutine clean_mommat

subroutine clean_mommat_final
!  PURPOSE
!    deallocation of big arrays
!  USES
   use runtime_choices
   use aims_memory_tracking, only : aims_deallocate
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(mommat_full_oned_up))  call aims_deallocate( mommat_full_oned_up,   "mommat_full_oned_up" )
  if (allocated(mommat_full_oned_low)) call aims_deallocate( mommat_full_oned_low, "mommat_full_oned_low" )
  if (allocated(mommat_full_oned_two_up))  &
             call aims_deallocate( mommat_full_oned_two_up,                     "mommat_full_oned_two_up" )
  if (allocated(mommat_full_oned_two_low))  &
             call aims_deallocate( mommat_full_oned_two_low,                   "mommat_full_oned_two_low" )

end subroutine clean_mommat_final

subroutine get_state_minmax_K(KS_eigen, n_state_min_out, n_state_max_out)

  !  PURPOSE
  !   Get the highest and lowest states that are smaller/bigger than E_max/E_min 
  !
  ! USES
    use constants, only: pi
    use pbc_lists
    use runtime_choices
    use aims_memory_tracking, only : aims_allocate, aims_deallocate
    use basis
    use dimensions, only : n_states, n_spin, n_k_points
    implicit none
    real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) ::     KS_eigen
    integer, INTENT(OUT) :: n_state_min_out
    integer, INTENT(OUT) :: n_state_max_out
    real*8 :: homo_level, lumo_level 

    real*8 :: Emin_ha
    real*8 :: Emax_ha
    integer :: i, i_k_point, i_spin
    integer, allocatable :: n_state_min_k(:,:)
    integer, allocatable :: n_state_max_k(:,:)
    call aims_allocate( n_state_min_k, n_k_points,n_spin, "n_state_min_k" )
    call aims_allocate( n_state_max_k, n_k_points,n_spin, "n_state_max_k" )
    n_state_min_k = 1 
    n_state_max_k = 0 

    Emin_ha=Emin/hartree
    Emax_ha=Emax/hartree
    do i_k_point = 1, n_k_points, 1
      do i_spin = 1, n_spin, 1
        do i = 1, n_states, 1
          if (KS_eigen(i,i_spin,i_k_point)<=Emin_ha) then
            n_state_min_k(i_k_point,i_spin)=n_state_min_k(i_k_point,i_spin)+1
          endif
          if (KS_eigen(i,i_spin,i_k_point)<=Emax_ha) then
            n_state_max_k(i_k_point,i_spin)=n_state_max_k(i_k_point,i_spin)+1
          endif
        enddo
      enddo
      n_state_min_out = min(minval(n_state_min_k(:,1)), &
           minval(n_state_min_k(:,n_spin)))
      n_state_max_out = max(maxval(n_state_max_k(:,1)), &
           maxval(n_state_max_k(:,n_spin)))
    enddo ! k-point loop
    if (n_state_min_out>=n_states) then
       n_state_min_out=n_states-1
    endif
    if (n_state_max_out>n_states) then
       n_state_max_out=n_states
    endif
    if (n_state_max_out==0) then
       n_state_max_out=2
    endif
    if (n_state_min_out>=n_state_max_out) then
       n_state_min_out=n_state_max_out-1
    endif
    call aims_deallocate( n_state_min_k, "n_state_min_k" )
    call aims_deallocate( n_state_max_k, "n_state_max_k" )
end subroutine get_state_minmax_K

subroutine calculate_mommat_p0( partition_tab_std, basis_l_max, &
                       mommat_full_oned_up, mommat_full_oned_low, coord)
!  PURPOSE
!  The subroutine integrates the momentum matrix
!  using a fixed basis set (no adaptive modifications).
!  
!  Adapted from integrate_olvp_matrix and integrate_hamitonian_matrix
!
!  USES
  use dimensions,            only : n_full_points, n_species, n_centers_integrals, &
                                    n_max_batch_size, n_max_compute_fns_ham, &
                                    n_centers_basis_I, n_basis_fns, l_wave_max, &
                                    n_hamiltonian_matrix_size, n_max_compute_atoms, &
                                    n_max_compute_ham, n_my_batches, n_periodic, &
                                    n_spin, n_centers
  use basis,                 only : basis_wave_ordered, basis_deriv_ordered
  use runtime_choices,       only : prune_basis_once, use_local_index
  use grids,                 only : batch_of_points, batches
  use pbc_lists,             only : centers_basis_integrals, inv_centers_basis_integrals 
  use synchronize_mpi_basic, only : sync_vector, sync_real_number
  use load_balancing,        only : batch_perm, use_batch_permutation, &
                                    get_batch_weights, set_batch_weights
  use species_data,          only : species_name
  use physics,               only : deallocate_vdw_splines
  use aims_memory_tracking , only : aims_allocate, aims_deallocate
  implicit none
!  ARGUMENTS
  real*8, target, dimension(n_full_points), intent(in)  :: partition_tab_std
  integer,        dimension(n_species),     intent(in)  :: basis_l_max
  ! when this routine is called, dipole_mat has either the dimension
  ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
  ! so we declare it here as a 1D assumed size array
  real*8,         dimension(*),             intent(out) :: mommat_full_oned_up
  real*8,         dimension(*),             intent(out) :: mommat_full_oned_low
  integer,                                  intent(in)  :: coord
!  INPUTS
!  o partition_tab_std -- values of partition function
!  o basis_l_max -- maximum l component of basis functions
!  OUTPUT
!  o overlap_matrix -- overlap matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  !  local variables
  real*8,  dimension(:,:,:), allocatable :: work_ham
  real*8,  dimension(:,:),   allocatable :: work_ovl

  integer                                :: l_ylm_max
  integer, dimension(:,:),   allocatable :: index_lm
  real*8,  dimension(:,:),   allocatable :: ylm_tab

  real*8,  dimension(:,:),   allocatable :: dylm_dtheta_tab
  real*8,  dimension(:,:),   allocatable :: scaled_dylm_dphi_tab
  real*8,  dimension(:),     allocatable :: dist_tab_full
  real*8,  dimension(:,:),   allocatable :: dir_tab_full_norm
  real*8,  dimension(:),     allocatable :: i_r_full

  real*8,  dimension(:),     allocatable :: radial_wave
  real*8,  dimension(:),     allocatable :: radial_wave_deriv
  real*8,  dimension(:,:),   allocatable :: wave
  real*8,  dimension(:,:,:), allocatable :: gradient_basis_wave
  real*8 coord_current(3)
  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 i_r(n_centers_integrals)
  real*8 dir_tab(3, n_centers_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_integrals)

  !    Auxiliary Hamiltonian matrix, to sum up contributions from only a 
  !    single integration shell. The hope is that such a separate treatment 
  !    will allow to minimize numerical noise introduced through ZORA
  real*8, dimension(:), allocatable :: mommat_shell
  !     optimal accounting for matrix multiplications: only use points with 
  !     nonzero components
  integer :: n_points

  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  !     for pruning of atoms, radial functions, and basis functions, to only 
  !     the relevant ones ...
  integer :: n_compute_a, n_compute_c
  integer :: i_basis(n_centers_basis_I)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_max_compute_atoms)
  integer :: spline_array_end(n_max_compute_atoms)

! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms)

      ! indices for basis functions that are nonzero at current point
!
! VB: Warning - n_max_compute_fns_ham is outrightly wrong for the density
!     update!!
!
!     This should MUCH rather be counted in prune_basis for each batch and then
!     passed here, that would be the appropriate maximum!
!
  integer :: rad_index(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_ham)
  integer :: l_index(n_max_compute_fns_ham)
  integer :: l_count(n_max_compute_fns_ham)
  integer :: fn_atom(n_max_compute_fns_ham)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute
  integer :: zero_index_point(n_max_compute_ham)

  ! active atoms in current batch
  integer :: n_batch_centers
  integer :: batch_center(n_centers_integrals)

  ! counters
  integer i_index, i_l, i_m

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_my_batch

  character*200 :: info_str
  integer :: info

  ! Load balancing stuff
  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches 
                                                     ! actually used

  integer dim_mommat ! actual dimension of momentum matrices

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  !  begin work
  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating real-space momentum matrix: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating real-space momentum matrix: batch-based integration."
  endif
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  ! begin with general allocations
  l_ylm_max = l_wave_max
  call aims_allocate( ylm_tab, (l_ylm_max+1)**2, n_max_compute_atoms,                           "ylm_tab" )
  call aims_allocate( index_lm, -l_ylm_max,l_ylm_max, 0,l_ylm_max,                             "index_lm" )
  call aims_allocate( radial_wave, n_max_compute_fns_ham,                                   "radial_wave" )
  call aims_allocate( radial_wave_deriv, n_max_compute_fns_ham,                       "radial_wave_deriv" )
  call aims_allocate( wave, n_max_compute_ham, n_max_batch_size,                                   "wave" )
  call aims_allocate( gradient_basis_wave, n_max_compute_ham,3,n_max_batch_size,    "gradient_basis_wave" )
  call aims_allocate( dylm_dtheta_tab, (l_ylm_max+1)**2, n_max_compute_atoms,           "dylm_dtheta_tab" )
  call aims_allocate( scaled_dylm_dphi_tab, (l_ylm_max+1)**2, n_max_compute_atoms, "scaled_dylm_dphi_tab" )
  call aims_allocate( dist_tab_full, n_centers_integrals,                           "n_centers_integrals" )
  call aims_allocate( dir_tab_full_norm, 3,n_centers_integrals,                       "dir_tab_full_norm" )
  call aims_allocate( i_r_full, n_centers_integrals,                                           "i_r_full" )
  call aims_allocate( mommat_shell, n_max_compute_ham*n_max_compute_ham,                   "mommat_shell" )

  call aims_allocate( work_ham, 1, 1, n_spin,                                                  "work_ham" )
  call aims_allocate( work_ovl, 1, 1,                                                          "work_ovl" )
  !-----------------------------------------------------------------------------
  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points 
  ! (for load balancing) or to standard batches / arrays (no load balancing)
  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then
    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    call aims_allocate( ins_idx, batch_perm(n_bp)%n_basis_local,                                "ins_idx" )

    dim_mommat = batch_perm(n_bp)%n_local_matrix_size
  else
    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

    dim_mommat = n_hamiltonian_matrix_size
  endif

  if(get_batch_weights) call aims_allocate( batch_times, n_my_batches_work,                 "batch_times" )
  !-----------------------------------------------------------------------------

  ! initialize
  mommat_full_oned_up(1:dim_mommat)  = 0.d0
  mommat_full_oned_low(1:dim_mommat) = 0.d0
  i_basis_fns_inv = 0

  ! initialize index_lm
  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  i_full_points = 0
  i_full_points_2 = 0

  call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_my_batch = 1, n_my_batches_work, 1
    if (get_batch_weights) time_start = mpi_wtime()

    n_compute_c = 0
    n_compute_a = 0
    i_basis = 0
    i_point = 0

    gradient_basis_wave = 0.0d0
    partition           = 0.0d0
    mommat_shell        = 0.0d0

    ! loop over one batch
    do i_index = 1, batches_work(i_my_batch)%size, 1
      i_full_points_2 = i_full_points_2 + 1
      if (partition_tab(i_full_points_2).gt.0.d0) then
        i_point = i_point+1

        ! get current integration point coordinate
        coord_current(:) = batches_work(i_my_batch) % points(i_index) &
                           % coords(:)
        if (n_periodic > 0)then
          call map_to_center_cell(coord_current(1:3) )
        end if

        ! compute atom-centered coordinates of current integration point,
        ! as viewed from all atoms
        call tab_atom_centered_coords_p0( coord_current, &
             dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
             n_centers_integrals, centers_basis_integrals )

        ! determine which basis functions are relevant at current 
        ! integration point, and tabulate their indices

        ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) 
        ! are actually needed
        if (.not.prune_basis_once) then
          call prune_basis_p2( dist_tab_sq(1,i_point), n_compute_c, &
               i_basis, n_centers_basis_I, n_centers_integrals, &
               inv_centers_basis_integrals )
        end if
      end if
    end do ! end loop over one batch

    if (prune_basis_once) then
      n_compute_c = batches_work(i_my_batch)%batch_n_compute
      i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
    end if
    ! from list of n_compute active basis functions in batch, collect all 
    ! atoms that are ever needed in batch.
    call collect_batch_centers_p2 &
         ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, &
           inv_centers_basis_integrals, n_batch_centers, batch_center )

    n_points = i_point

    ! Perform actual integration if more than 0 basis functions
    ! are actually relevant on the present angular shell ...
    if (n_compute_c.gt.0) then
      i_point = 0
      ! loop over one batch of integration points
      do i_index = 1, batches_work(i_my_batch)%size, 1
        ! Increment the (global) counter for the grid, to access storage 
        ! arrays
        i_full_points = i_full_points + 1
        if (partition_tab(i_full_points).gt.0.d0) then
          i_point = i_point+1

          ! for all integrations
          n_compute_atoms = 0
          n_compute_fns = 0
          partition(i_point) = partition_tab(i_full_points)

          ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) 
          ! if needed). Are stored in a compact spline array that can be 
          ! accessed by spline_vector_waves,without any copying and 
          ! without doing any unnecessary operations. The price is that 
          ! the interface is no longer explicit in terms of physical 
          ! objects. See shrink_fixed_basis() for details regarding the 
          ! reorganized spline arrays.
          call tab_global_geometry_p0 &
               ( dist_tab_sq(1,i_point), &
                dir_tab(1,1,i_point), &
                dist_tab_full, &
                i_r_full, &
                dir_tab_full_norm, &
                n_centers_integrals,  centers_basis_integrals )

          call prune_radial_basis_p2 &
               ( n_max_compute_atoms, n_max_compute_fns_ham, &
                 dist_tab_sq(1,i_point), dist_tab(1,i_point), &
                 dir_tab(1,1,i_point), n_compute_atoms, atom_index, &
                 atom_index_inv, n_compute_fns, i_basis_fns, &
                 i_basis_fns_inv, i_atom_fns, spline_array_start, &
                 spline_array_end, n_centers_integrals, &
                 centers_basis_integrals, n_compute_c, i_basis, &
                 n_batch_centers, batch_center, &
                 one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                 fn_atom, n_zero_compute, zero_index_point )

          ! Tabulate distances, unit vectors, and inverse logarithmic 
          ! grid units for all atoms which are actually relevant
          call tab_local_geometry_p2 &
               ( n_compute_atoms, atom_index, &
                 dist_tab(1,i_point), i_r )

          ! compute trigonometric functions of spherical coordinate
          ! angles of current integration point, viewed from all atoms
          call tab_trigonom_p0 &
               ( n_compute_atoms, dir_tab(1,1,i_point),  &
                 trigonom_tab )

          ! tabulate distance and Ylm's w.r.t. other atoms
          call tab_wave_ylm_p0 &
               ( n_compute_atoms, atom_index,  &
                 trigonom_tab, basis_l_max,  &
                 l_ylm_max, &
                 ylm_tab )

          ! Now evaluate radial functions
          ! from the previously stored compressed spline arrays
          call evaluate_radial_functions_p0  &
               ( spline_array_start, spline_array_end,  &
                 n_compute_atoms, n_compute_fns,   &
                 dist_tab(1,i_point), i_r,  &
                 atom_index, i_basis_fns_inv,  &
                 basis_wave_ordered, radial_wave,  &
                 .false. , n_compute_c, n_max_compute_fns_ham )

          ! tabulate total wave function value for each basis function
          call evaluate_waves_p2  &
               ( n_compute_c, n_compute_atoms, n_compute_fns, &
                 l_ylm_max, ylm_tab, one_over_dist_tab,   &
                 radial_wave(1), wave(1,i_point), &
                 rad_index, wave_index, l_index, l_count, fn_atom, &
                 n_zero_compute, zero_index_point )

          call tab_gradient_ylm_p0  &
               ( trigonom_tab(1,1), basis_l_max,   &
                 l_ylm_max, n_compute_atoms, atom_index,  &
                 ylm_tab(1,1),   &
                 dylm_dtheta_tab(1,1),   &
                 scaled_dylm_dphi_tab(1,1)  )
          call evaluate_radial_functions_p0  &
               ( spline_array_start, spline_array_end,  &
                 n_compute_atoms, n_compute_fns,   &
                 dist_tab(1,i_point), i_r,  &
                 atom_index, i_basis_fns_inv,  &
                 basis_deriv_ordered,   &
                 radial_wave_deriv(1), .true.,  &
                 n_compute_c, n_max_compute_fns_ham )

          ! and finally, assemble the actual gradients
          call evaluate_wave_gradient_p2  &
               ( n_compute_c, n_compute_atoms, n_compute_fns, &
                 one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1), &
                 l_ylm_max, ylm_tab,  &
                 dylm_dtheta_tab,  &
                 scaled_dylm_dphi_tab,  &
                 radial_wave,  &
                 radial_wave_deriv,  &
                 gradient_basis_wave(1:n_compute_c,1:3,i_point),  &
                 rad_index, wave_index, l_index, l_count, fn_atom, &
                 n_zero_compute, zero_index_point )
! new version: gradient_basis_wave(1,1,i_point) - probably reason to crash
! solll 1:n_compute, 1:3, i_point sein!!!

          ! Reset i_basis_fns_inv
          i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
        end if ! partition_tab(i_full_points).gt.0.d0
      end do ! i_index
      ! end loop over one batch

      call evaluate_mommat( & 
           n_points, partition, n_compute_c, &
           gradient_basis_wave, &
           n_max_compute_ham, wave, &
           mommat_shell, coord )

      if (use_batch_permutation > 0) then
        ! If use_batch_permutation > 0 is set, the local hamiltonian is
        ! always stored in full form for the local basis functions

        ! Get position of basis functions of current batch within local 
        ! hamiltonian
        do i=1,n_compute_c
          ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
        end do

        ! Insert hamiltonian_shell of current batch
        do i=1,n_compute_c
          i_off = (ins_idx(i)*(ins_idx(i)-1))/2
          do j=1,i ! n_compute_c
            mommat_full_oned_up(ins_idx(j)+i_off) = &
                 mommat_full_oned_up(ins_idx(j)+i_off) &
                 + mommat_shell(j+(i-1)*n_compute_c)
            mommat_full_oned_low(ins_idx(j)+i_off) = &
                 mommat_full_oned_low(ins_idx(j)+i_off)&
                 + mommat_shell(i+(j-1)*n_compute_c)
          end do
        end do
      else
        !keep upper and lower Triangle
        call update_full_matrix_mommat &
             ( n_compute_c, n_compute_c,  &
               i_basis(1), &
               mommat_shell,  &
               mommat_full_oned_up(1) )
        call transpose_moment_shell(mommat_shell,n_compute_c)
        call update_full_matrix_mommat &
             ( n_compute_c, n_compute_c,  &
               i_basis(1), &
               mommat_shell, &
               mommat_full_oned_low(1) )
      endif
    else
      i_full_points = i_full_points + batches_work(i_my_batch)%size
    end if ! n_compute.gt.0

    if (get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start
  end do ! end loop over batches

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for '//&
                                   'integration: real work ', &
                                   time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)

  !     synchronise the hamiltonian
  if(.not. use_local_index) call sync_vector(mommat_full_oned_up, &
                                             n_hamiltonian_matrix_size )
  if(.not. use_local_index) call sync_vector(mommat_full_oned_low, &
                                             n_hamiltonian_matrix_size )

  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)

  if (allocated(ylm_tab)) then
     call aims_deallocate( ylm_tab,                           "ylm_tab" )
  end if
  if (allocated(index_lm)) then
     call aims_deallocate( index_lm,                         "index_lm" )
  end if

  if (allocated(mommat_shell)) then
     call aims_deallocate( mommat_shell,                 "mommat_shell" )
  end if

  if (allocated(ins_idx)) then
     call aims_deallocate( ins_idx,                           "ins_idx" )
  endif

  if (allocated(batch_times)) then
     call aims_deallocate( batch_times,                   "batch_times" )
  endif

  if (allocated(radial_wave)) then
     call aims_deallocate( radial_wave,                   "radial_wave" )
  endif

  if (allocated(radial_wave_deriv)) then
     call aims_deallocate( radial_wave_deriv,       "radial_wave_deriv" )
  endif

  if (allocated(wave)) then
     call aims_deallocate( wave,                                 "wave" )
  endif

  if (allocated(gradient_basis_wave)) then
     call aims_deallocate( gradient_basis_wave,   "gradient_basis_wave" )
  endif

  if (allocated(dylm_dtheta_tab)) then
     call aims_deallocate( dylm_dtheta_tab,           "dylm_dtheta_tab" )
  endif

  if (allocated(scaled_dylm_dphi_tab)) then
     call aims_deallocate( scaled_dylm_dphi_tab, "scaled_dylm_dphi_tab" )
  endif

  if (allocated(dist_tab_full)) then
     call aims_deallocate( dist_tab_full,               "dist_tab_full" )
  endif

  if (allocated(dir_tab_full_norm)) then
     call aims_deallocate( dir_tab_full_norm,       "dir_tab_full_norm" )
  endif

  if (allocated(i_r_full)) then
     call aims_deallocate( i_r_full,                         "i_r_full" )
  endif

  if (allocated(work_ham)) then
     call aims_deallocate( work_ham,                          "work_ham" )
  end if

  if (allocated(work_ovl)) then
     call aims_deallocate( work_ovl,                          "work_ovl" )
  end if 
end subroutine calculate_mommat_p0

subroutine calculate_dipmat_p0( partition_tab_std, basis_l_max, &
                       dipmat_full_oned_up, dipmat_full_oned_low, coord)

!  PURPOSE
!  The subroutine integrates the dipole matrix
!  using a fixed basis set (no adaptive modifications).
!  
!  Adapted from integrate_olvp_matrix and integrate_hamiltonian_matrix
!
!  USES
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use plus_u
  use mpi_utilities
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data, only: species_name
  use load_balancing
  use pbc_lists, only: centers_basis_integrals, inv_centers_basis_integrals
  use physics, only: deallocate_vdw_splines
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points) :: partition_tab_std
  integer basis_l_max (n_species)
  integer coord
  ! when this routine is called, dipole_mat has either the dimension
  ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
  ! so we declare it here as a 1D assumed size array
  real*8 :: dipmat_full_oned_up(*)
  real*8 :: dipmat_full_oned_low(*)

!  INPUTS
!  o partition_tab_std -- values of partition function
!  o basis_l_max -- maximum l component of basis functions
!
!  OUTPUT
!  o overlap_matrix -- overlap matrix
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  real*8 coord_current_new(3)
  real*8, dimension(:,:),allocatable :: coordinates
  !  local variables
  real*8, dimension(:,:,:),allocatable :: work_ham
  real*8, dimension(:,:),allocatable :: work_ovl

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:), allocatable :: dist_tab_full
  real*8, dimension(:,:), allocatable :: dir_tab_full_norm
  real*8, dimension(:), allocatable :: i_r_full

  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:,:)  ,allocatable:: wave
  real*8 coord_current(3)
  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 i_r(n_centers_integrals)
  real*8 dir_tab(3, n_centers_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_integrals)


  !    Auxiliary Hamiltonian matrix, to sum up contributions from only a 
  !    single integration shell. The hope is that such a separate treatment 
  !    will allow to minimize numerical noise introduced through ZORA
  real*8, dimension(:,:), allocatable :: dipolemat_shell
  !     optimal accounting for matrix multiplications: only use points with 
  !     nonzero components
  integer :: n_points

  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  !     for pruning of atoms, radial functions, and basis functions, to only 
  !     the relevant ones ...

  integer :: n_compute_a, n_compute_c
  integer :: i_basis(n_centers_basis_I)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_max_compute_atoms)
  integer :: spline_array_end(n_max_compute_atoms)

! VB - renewed index infrastructure starts here

      real*8 one_over_dist_tab(n_max_compute_atoms)

      ! indices for basis functions that are nonzero at current point
!
! VB: Warning - n_max_compute_fns_ham is outrightly wrong for the density
!     update!!
!
!     This should MUCH rather be counted in prune_basis for each batch and then
!     passed here, that would be the appropriate maximum!
!

      integer :: rad_index(n_max_compute_atoms)
      integer :: wave_index(n_max_compute_fns_ham)
      integer :: l_index(n_max_compute_fns_ham)
      integer :: l_count(n_max_compute_fns_ham)
      integer :: fn_atom(n_max_compute_fns_ham)

      ! indices for known zero basis functions at current point
      integer :: n_zero_compute
      integer :: zero_index_point(n_max_compute_ham)

      ! active atoms in current batch
      integer :: n_batch_centers
      integer :: batch_center(n_centers_integrals)

  !     for splitting of angular shells into "octants"


  !  counters




  integer i_index, i_l, i_m

  integer  i_center_L




  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_my_batch

  character*200 :: info_str
  integer :: info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches 
                                                     ! actually used

  integer dim_overlap_matrix ! actual dimension of overlap matrix

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  !  begin work

  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A,A)')"Integrating Momentum Matrix with loadbalancing"
  else
    write(info_str,'(2X,A,A)')"Integrating momentum matrix."
  endif
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  !     begin with general allocations
  l_ylm_max = l_wave_max
  call aims_allocate( coordinates, 3, n_max_batch_size,                         "coordinates" )
  call aims_allocate( ylm_tab, (l_ylm_max+1)**2, n_max_compute_atoms,               "ylm_tab" )
  call aims_allocate( index_lm, -l_ylm_max,l_ylm_max, 0,l_ylm_max,                 "index_lm" )
  call aims_allocate( radial_wave, n_max_compute_fns_ham,                       "radial_wave" )
  call aims_allocate( wave, n_max_compute_ham, n_max_batch_size,                       "wave" )
  call aims_allocate( dist_tab_full, n_centers_integrals,                     "dist_tab_full" )
  call aims_allocate( dir_tab_full_norm, 3,n_centers_integrals,           "dir_tab_full_norm" )
  call aims_allocate( i_r_full, n_centers_integrals,                               "i_r_full" )
  call aims_allocate( dipolemat_shell, n_max_compute_ham,n_max_compute_ham, "dipolemat_shell" )
  call aims_allocate( work_ham, 1, 1, n_spin,                                      "work_ham" )
  call aims_allocate( work_ovl, 1, 1,                                              "work_ovl" )
  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points 
  ! (for load balancing) or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    call aims_allocate( ins_idx, batch_perm(n_bp)%n_basis_local,                    "ins_idx" )

    dim_overlap_matrix = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

    dim_overlap_matrix = n_hamiltonian_matrix_size

  endif

  if(get_batch_weights) call aims_allocate( batch_times, n_my_batches_work,     "batch_times" )

  !-----------------------------------------------------------------------------

  ! initialize

  dipmat_full_oned_up(1:dim_overlap_matrix) = 0.d0
  dipmat_full_oned_low(1:dim_overlap_matrix) = 0.d0
  i_basis_fns_inv = 0
  ! initialize index_lm

  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  i_full_points = 0
  i_full_points_2 = 0

  call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_my_batch = 1, n_my_batches_work, 1

        if(get_batch_weights) time_start = mpi_wtime()

        n_compute_c = 0
        n_compute_a = 0
        i_basis = 0

        i_point = 0

        ! loop over one batch
        do i_index = 1, batches_work(i_my_batch)%size, 1

           i_full_points_2 = i_full_points_2 + 1

           if (partition_tab(i_full_points_2).gt.0.d0) then

              i_point = i_point+1

              !     get current integration point coordinate
              coord_current(:) = batches_work(i_my_batch) % points(i_index) &
                                 % coords(:)

              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if

              !     compute atom-centered coordinates of current integration point,
              !     as viewed from all atoms
              call tab_atom_centered_coords_p0( coord_current, &
                   dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
                               n_centers_integrals, centers_basis_integrals )

              !    determine which basis functions are relevant at current 
              !    integration point, and tabulate their indices

              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) 
              ! are actually needed
              if (.not.prune_basis_once) then
                 call prune_basis_p2( dist_tab_sq(1,i_point), n_compute_c, &
                      i_basis, n_centers_basis_I, n_centers_integrals, &
                      inv_centers_basis_integrals )
              end if
           end if

        enddo ! end loop over one batch

        if (prune_basis_once) then
           n_compute_c = batches_work(i_my_batch)%batch_n_compute
           i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
        end if
        ! from list of n_compute active basis functions in batch, collect all 
        ! atoms that are ever needed in batch.
        call collect_batch_centers_p2 &
        ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, &
          inv_centers_basis_integrals, n_batch_centers, batch_center &
        )

        n_points = i_point

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_c.gt.0) then

           i_point = 0

           ! loop over one batch of integration points
           do i_index = 1, batches_work(i_my_batch)%size, 1

              ! Increment the (global) counter for the grid, to access storage 
              ! arrays
              i_full_points = i_full_points + 1

              if (partition_tab(i_full_points).gt.0.d0) then

                 i_point = i_point+1

                 ! for all integrations
                 partition(i_point) = partition_tab(i_full_points)

                 n_compute_atoms = 0
                 n_compute_fns = 0

                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) 
                 ! if needed). Are stored in a compact spline array that can be 
                 ! accessed by spline_vector_waves,without any copying and 
                 ! without doing any unnecessary operations. The price is that 
                 ! the interface is no longer explicit in terms of physical 
                 ! objects. See shrink_fixed_basis() for details regarding the 
                 ! reorganized spline arrays.

                 call tab_global_geometry_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dir_tab(1,1,i_point), &
                      dist_tab_full, &
                      i_r_full, &
                      dir_tab_full_norm, &
                      n_centers_integrals,  centers_basis_integrals)

                 call prune_radial_basis_p2 &
                   ( n_max_compute_atoms, n_max_compute_fns_ham, &
                     dist_tab_sq(1,i_point), dist_tab(1,i_point), &
                     dir_tab(1,1,i_point), n_compute_atoms, atom_index, &
                     atom_index_inv, n_compute_fns, i_basis_fns, &
                     i_basis_fns_inv, i_atom_fns, spline_array_start, &
                     spline_array_end, n_centers_integrals, &
                     centers_basis_integrals, n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                     fn_atom, n_zero_compute, zero_index_point &
                    )

                 ! Tabulate distances, unit vectors, and inverse logarithmic 
                 ! grid units for all atoms which are actually relevant
                 call tab_local_geometry_p2 &
                      ( n_compute_atoms, atom_index, &
                        dist_tab(1,i_point),  &
                      i_r )

                 ! compute trigonometric functions of spherical coordinate
                 ! angles of current integration point, viewed from all atoms
                 call tab_trigonom_p0 &
                      ( n_compute_atoms, dir_tab(1,1,i_point),  &
                      trigonom_tab )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm_p0 &
                      ( n_compute_atoms, atom_index,  &
                      trigonom_tab, basis_l_max,  &
                      l_ylm_max, &
                      ylm_tab )

                 ! Now evaluate radial functions
                 ! from the previously stored compressed spline arrays
                 call evaluate_radial_functions_p0  &
                      ( spline_array_start, spline_array_end,  &
                      n_compute_atoms, n_compute_fns,   &
                      dist_tab(1,i_point), i_r,  &
                      atom_index, i_basis_fns_inv,  &
                      basis_wave_ordered, radial_wave,  &
                      .false. , n_compute_c, n_max_compute_fns_ham )

                 ! tabulate total wave function value for each basis function
                 call evaluate_waves_p2  &
                   ( n_compute_c, n_compute_atoms, n_compute_fns, &
                     l_ylm_max, ylm_tab, one_over_dist_tab,   &
                     radial_wave(1), wave(1,i_point), &
                     rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point &
                   )


                 ! Reset i_basis_fns_inv
                 i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
		 coord_current_new(:) = batches_work(i_my_batch) &
                                        % points(i_index)  % coords(:)
		 if(n_periodic > 0)then
		   call map_to_center_cell(coord_current_new(1:3) )
		 end if
		 coordinates(1,i_point)=coord_current_new(1)
		 coordinates(2,i_point)=coord_current_new(2)
		 coordinates(3,i_point)=coord_current_new(3)
               end if

              ! end loop over one batch
           enddo
	      call evaluate_dipmat( & 
		n_points, partition(1), n_compute_c, &
		coordinates, &
		n_max_compute_ham, wave(1,1), &
		dipolemat_shell,coord)
           if(use_batch_permutation > 0) then
              ! If use_batch_permutation > 0 is set, the local hamiltonian is
              ! always stored in full form for the local basis functions

              ! Get position of basis functions of current batch within local 
              ! hamiltonian
              do i=1,n_compute_c
                 ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
              enddo

              ! Insert hamiltonian_shell of current batch
              do i=1,n_compute_c
                 i_off = (ins_idx(i)*(ins_idx(i)-1))/2
                 do j=1,i ! n_compute_c
                    dipmat_full_oned_up(ins_idx(j)+i_off) = &
                                     dipmat_full_oned_up(ins_idx(j)+i_off) &
                                   + dipolemat_shell(j,(i-1)*n_compute_c)
                    dipmat_full_oned_low(ins_idx(j)+i_off) = &
                                     dipmat_full_oned_low(ins_idx(j)+i_off)&
                                   + dipolemat_shell(i,(j-1)*n_compute_c)
                 enddo
              enddo
           else
              !keep upper and lower Triangle
	      call update_full_matrix_mommat &
		      ( n_compute_c, n_compute_c,  &
		      i_basis, &
		      dipolemat_shell,  &
		      dipmat_full_oned_up)
              call transpose_moment_shell(dipolemat_shell,n_compute_c)
	      call update_full_matrix_mommat &
		      ( n_compute_c, n_compute_c,  &
		      i_basis, &
		      dipolemat_shell, &
		      dipmat_full_oned_low)
           !enddo
           endif
        else
           i_full_points = i_full_points + batches_work(i_my_batch)%size
        end if ! n_compute.gt.0

        if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  enddo ! end loop over batches

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for '//&
                                   'integration: real work ', &
                                   time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)

  !     synchronise the hamiltonian
  if(.not. use_local_index) call sync_vector(dipmat_full_oned_up, &
                                             n_hamiltonian_matrix_size )
  if(.not. use_local_index) call sync_vector(dipmat_full_oned_low, &
                                             n_hamiltonian_matrix_size )

  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)


  if (allocated(ylm_tab)) then
     call aims_deallocate( ylm_tab,                     "ylm_tab" )
  end if
  if (allocated(index_lm)) then
     call aims_deallocate( index_lm,                   "index_lm" )
  end if

  if (allocated(dipolemat_shell)) then
     call aims_deallocate( dipolemat_shell,     "dipolemat_shell" )
  end if

  if (allocated(coordinates)) then
     call aims_deallocate( coordinates,             "coordinates" )
  end if

  if (allocated(ins_idx)) then
     call aims_deallocate( ins_idx,                     "ins_idx" )
  endif

  if (allocated(batch_times)) then
     call aims_deallocate( batch_times,             "batch_times" )
  endif

  if (allocated(radial_wave)) then
     call aims_deallocate( radial_wave,             "radial_wave" )
  endif

  if (allocated(wave)) then
     call aims_deallocate( wave,                           "wave" )
  endif


  if (allocated(dist_tab_full)) then
     call aims_deallocate( dist_tab_full,         "dist_tab_full" )
  endif

  if (allocated(dir_tab_full_norm)) then
     call aims_deallocate( dir_tab_full_norm, "dir_tab_full_norm" )
  endif

  if (allocated(i_r_full)) then
     call aims_deallocate( i_r_full,                   "i_r_full" )
  endif

  if (allocated(work_ham)) then
     call aims_deallocate( work_ham,                          "work_ham" )
  end if

  if (allocated(work_ovl)) then
     call aims_deallocate( work_ovl,                          "work_ovl" )
  end if 

end subroutine calculate_dipmat_p0

subroutine evaluate_mommat( & 
     n_points, partition, n_compute, gradient_basis_wave, &
     n_basis_list, wave, mommat_shell_out, i_coord )
! PURPOSE
!   Evaluates Momentummatrix for one grid batch,
! USES
  implicit none
! ARGUMENTS
  integer, intent(in)  :: n_basis_list
  integer, intent(in)  :: n_points  
  integer, intent(in)  :: n_compute 
  integer, intent(in)  :: i_coord
  real*8,  intent(in)  :: partition(n_points)
  real*8,  intent(in)  :: wave(n_basis_list, n_points)
  real*8,  intent(in)  :: gradient_basis_wave(n_basis_list, 3, n_points) 
  real*8,  intent(out) :: mommat_shell_out(n_compute, n_compute)
! INPUTS
! o n_points -- number of points in this batch
! o n_compute -- number of non zero basis functions in this batch
! o gradient_basis_wave -- the gradient of the basis function (x,y or z)
! o wave -- basis wave
! o i_coord -- component to calculate
! OUTPUT
! o mommat_shell
! AUTHOR
!   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
! SEE ALSO
!   Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!   Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!   "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!   Computer Physics Communications (2008), submitted.
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
! HISTORY
!   Release version, FHI-aims (2008).
! SOURCE
! begin work
  real*8 wave_compute(n_compute, n_points)
  real*8 grad_wave_compute(n_compute, n_points) 
! counters
  integer :: i_point

  mommat_shell_out = 0.0
! Condense basis functions to only those that are used at present
! integration points
  do i_point = 1, n_points, 1
    wave_compute(1:n_compute, i_point) =  &
         sqrt(partition(i_point))* &
         wave(1:n_compute, i_point)
    grad_wave_compute(1:n_compute, i_point) =  &
         sqrt(partition(i_point))*gradient_basis_wave(1:n_compute,i_coord,&
         i_point)   
  end do

  call dgemm('N', 'T', n_compute, n_compute,  &
       n_points, 1.0d0,  &
       wave_compute, n_compute, grad_wave_compute,  &
       n_compute, 0.0d0, mommat_shell_out, n_compute )
end subroutine evaluate_mommat

subroutine evaluate_dipmat( & 
     n_points, partition, n_compute, coords, &
     n_basis_list, wave, dipmat_shell_out,i_coord)

!  PURPOSE
!  Evaluates Dipolematrix for one grid batch,
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer, Intent(IN) :: n_basis_list
  integer, Intent(IN) :: n_points  
  integer, Intent(IN) :: n_compute 
  real*8, Intent(IN)  :: coords(1:3,n_points)
  integer, Intent(IN) :: i_coord
  real*8, Intent(IN)  :: partition(n_points)
  real*8, Intent(IN)  :: wave(n_basis_list, n_points)
  real*8, Intent(OUT) :: dipmat_shell_out( n_compute, n_compute)

! INPUTS
! o n_points -- number of points in this batch
! o n_compute -- number of non zero basis functions in this batch
! o gradient_basis_wave -- the gradient of the basis function (x,y or z)
! o wave -- basis wave
! o i_coord -- component to calculate
!  OUTPUT
! o dipmat_shell
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     begin work
      real*8 wave_compute(n_compute,n_points)
      real*8 grad_wave_compute(n_compute,n_points) 
!     counters
      integer :: i_point

      dipmat_shell_out = 0.0
!     Condense basis functions to only those that are used at present
!     integration points
      do i_point = 1, n_points, 1
         wave_compute(1:n_compute, i_point) =  &
              sqrt(partition(i_point))* &
                wave(1:n_compute, i_point) *  coords(i_coord,i_point)
         grad_wave_compute(1:n_compute, i_point) =  &
              sqrt(partition(i_point))*wave(1:n_compute, i_point)
                   
      enddo

        call dgemm('N', 'T', n_compute, n_compute,  &
            n_points, 1.0d0,  &
             wave_compute, n_compute, grad_wave_compute,  &
             n_compute, 0.0d0, dipmat_shell_out, n_compute )

end subroutine evaluate_dipmat


subroutine update_full_matrix_mommat &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell, matrix &
     )

!  PURPOSE
!  Subroutine update_full_matrix adds a part of the integrals in a
!  matrix (overlap or Hamiltonian matrix)
!  (only for the n_compute basis functions that are nonzero at the
!  current integration shell) to the full matrix (which contains
!  all n_basis basis functions). The link between i_compute = 1 ... n_compute
!  and i_basis = 1 ... n_basis is provided by the index array 
!  i_basis(i_compute).
!
!  USES

  use dimensions
  use runtime_choices
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  real*8 matrix_shell(n_compute_c,n_compute_a)
  real*8 matrix( n_hamiltonian_matrix_size)

!  INPUTS
!   o n_compute_c == n_compute_a -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
!   o matrix_shell -- hamiltonian / overlap_matrix of relevant basis functions
!
!  OUTPUT
!   o matrix -- date from matrix_shell is added here
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
!







  !  local variables

  !     auxiliary matrices for Level 3 Blas matrix multiplications

  !     counters

  integer :: i_compute_1
  integer :: i_compute_2
  integer :: i_offset
  integer :: i_index_real

  integer :: i_cell_1

  integer :: i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)

!  integer :: direct_citing(n_compute_a,n_compute_a )
!  real*8 :: data(n_cells_in_hamiltonian, n_basis)
!  real*8,dimension(:,:), allocatable :: data
  integer:: i_cell_old

!NEC_CB
  integer::help1(n_compute_a)
  integer::help2(n_compute_a)

  !     begin work

  !      now add the aux. ham. matrix to the actual ham. matrix
  !      this requires translating between the actually computed matrix elements 
  !      and the full ham. matrix ...

  ! When basis is smaller than n_basis

  !      write(use_unit,*) n_compute_c,n_compute_a

  select case(packed_matrix_format)

  case(PM_none)

     i_index_real = 0
     do i_compute_2 = 1, n_compute_a, 1

        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2


        do i_compute_1 = 1,i_compute_2,1


           i_index_real = i_offset + i_basis(i_compute_1) 

           matrix(i_index_real) = matrix(i_index_real) &
                + matrix_shell(i_compute_1, i_compute_2)


        enddo
     enddo


  case(PM_index) !-------------------------------------------------------------


     if(n_periodic == 0)then

        do i_compute_1 = 1, n_compute_a, 1

           
!           write(use_unit,*) 'i_compute_2', i_compute_2
           
           i_start =  index_hamiltonian(1,1, i_basis(i_compute_1))
           i_end   =  index_hamiltonian(2,1, i_basis(i_compute_1))


!           write(use_unit,*) '-'

           do i_compute_2 = 1,i_compute_1,1
              
              i_basis_2 = i_basis(i_compute_2)
              
              place: do i_place = i_start, i_end, 1
                 
                 if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
!NEC_CB             matrix(i_place) = matrix(i_place) + 
!                   matrix_shell(i_compute_2, i_compute_1)
                    if (i_compute_2.le.i_compute_1) then
                       matrix(i_place) = matrix(i_place) + &
                                         matrix_shell(i_compute_2, i_compute_1)
                    else
                       matrix(i_place) = matrix(i_place) + &
                                         matrix_shell(i_compute_2, i_compute_1)
                    end if

                    i_index_real = i_place
                    exit place 

                 
                 else if(column_index_hamiltonian( i_place) > i_basis_2)then
                    i_index_real = i_place
                    exit place  

                 end if
              end do place
              i_start = i_index_real
           end do
        end do
        


           
        
     else ! Periodic case----------------



        ! Unfortunately the periodic systems can not use the searching routine 
        ! used now in the clusters. This is because the peridic systems have 
        ! extra packing for supercell information.

!NEC_CB Do not read multi-indirectly addressed indices multiple times
        do i_compute_1 = 1, n_compute_a, 1
          help1(i_compute_1)=Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1)=&
                          center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do

        do i_compute_1 = 1, n_compute_a, 1

           i_basis_1 = help1(i_compute_1) !Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = help2(i_compute_1)
                         !center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))


           offset_end = -1
           offset = -1

           do i_cell_1 = 1, n_cells 


              i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

              offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
              offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

           end do



           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = help1(i_compute_2)
                                          !Cbasis_to_basis(i_basis(i_compute_2))



              if(i_basis_2 <= i_basis_1)then

                 i_cell    =  help2(i_compute_2)
                         !center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))


                 place_2: do i_place = offset(i_cell), offset_end(i_cell),1 

!                    write(use_unit,*) offset(i_cell), offset_end(i_cell)
                 
                    if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
!NEC_CB                       matrix(i_place) = matrix(i_place) +
!                                         matrix_shell(i_compute_2, i_compute_1)
                       if (i_compute_2.le.i_compute_1) then
                         matrix(i_place) = matrix(i_place) + &
                                          matrix_shell(i_compute_2, i_compute_1)
                       else
                         matrix(i_place) = matrix(i_place) + &
                                          matrix_shell(i_compute_2, i_compute_1)
                       end if

                       exit place_2
                 
                    else if(column_index_hamiltonian( i_place) > i_basis_2)then
                       exit place_2  


                    end if
                 end do place_2

                 offset(i_cell) = i_place

              end if
           end do
        end do

     end if ! n_periodic == 0 - else

  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select
end subroutine update_full_matrix_mommat

subroutine transpose_moment_shell( hamiltonian_shell, n_compute )

!  PURPOSE
!
!   Makes traspose of the table.
!   This is needed, because the dimensions of the gradient tables are different 
!   in subroutines (n_compute)  than in calling main subroutine (n_max_compute).
!
   implicit none 
!  ARGUMENTS
 
    real*8  :: hamiltonian_shell( n_compute,n_compute)
    integer :: n_compute

!  INPUTS
!    o hamiltonian_shell -- table before transpose
!    o  n_compute -- dimension of the table
!  OUTPUT
!    o hamiltonian_shell -- table after transpose
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    hamiltonian_shell = transpose(hamiltonian_shell)

end subroutine transpose_moment_shell

subroutine calc_moment_p0(moment, mommat_full_w_up, &
                          mommat_full_w_low, mommat_full_w_complex_up,&
                          mommat_full_w_complex_low, KS_vec, KS_vec_complex,&
                          k_point, i_coord, n_state_min_in, n_state_max_in)

  !  PURPOSE
  !   Calculete Momentum-Matrix elements by summing up of KS EV coefficents 
  !   (basis transformation)
  !
  ! USES
  use runtime_choices, only : real_eigenvectors
  use dimensions,      only : n_basis, n_states, n_spin
  implicit none
  !  ARGUMENTS
  integer, INTENT(IN):: n_state_min_in
  integer, INTENT(IN):: n_state_max_in
  complex*16, INTENT(IN)::     KS_vec_complex(n_basis, n_states, n_spin)
  real*8, INTENT(IN)::     KS_vec(n_basis, n_states, n_spin)
  real*8, INTENT(IN)::     mommat_full_w_up((n_basis+1)*n_basis/2)
  real*8, INTENT(IN)::     mommat_full_w_low((n_basis+1)*n_basis/2)
  complex*16, INTENT(IN):: mommat_full_w_complex_up((n_basis+1)*n_basis/2)
  complex*16, INTENT(IN):: mommat_full_w_complex_low((n_basis+1)*n_basis/2)
  integer, INTENT(IN):: k_point
  integer, INTENT(IN):: i_coord
  complex*16, INTENT(OUT) ::  moment((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  !  INPUTS
  !    o KS_eigenvector_complex -- Eigencoefficents at k_point (complex)
  !    o KS_eigenvector -- Eigencoefficents at k_point (real)
  !    o mommat_w_up/low  -- Momentummatrix (upper/lower triangle)
  !                             in atomic basis (real)
  !    o mommat_w_complex_up/low -- Momentummatrix (upper/lower triangle) 
  !                                    in atomic basis (complex)
  !    o k_point -- k-point wanted to be calculated 
  !  OUTPUT
  !    o dipelement -- Momentummatrix in KS-basis
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE
  integer:: i_spin,j_spin,i_basis_2, i_basis_1, n_state, m_state
  integer::  num, num_base
  real*8 :: reconstmommat(n_basis, n_basis), dummy(n_basis)
  complex*16 ::  reconstmommatcpl(n_basis, n_basis), dummycpl(n_basis)
  logical :: newflag
  real*8, external :: ddot
  complex*16, external :: zdotc

! I AM LEAVING IT FALSE (OLD WAY KEPT FOR NOW) UNTIL IT IS THOROUGHLY TESTED
  newflag=.true.
! MR OPTMIZE: build a new matrix that is actually "the" mommat matrix:
  num_base = 0
  if(real_eigenvectors) then
     do i_basis_2 = 1, n_basis
        do i_basis_1 = 1, i_basis_2
           num_base=num_base+1
           reconstmommat(i_basis_1, i_basis_2)=mommat_full_w_up(num_base)
           reconstmommat(i_basis_2, i_basis_1)=mommat_full_w_low(num_base)
        enddo
     enddo  
  else
     do i_basis_2 = 1, n_basis
        do i_basis_1 = 1, i_basis_2
           num_base=num_base+1
           reconstmommatcpl(i_basis_1, i_basis_2)=mommat_full_w_complex_up(num_base)
           reconstmommatcpl(i_basis_2, i_basis_1)=conjg(mommat_full_w_complex_low(num_base))
        enddo
     enddo 
  endif

  num = 0
  moment = (0.0,0.0)
  do i_spin = 1, n_spin
     do j_spin = i_spin, n_spin
       do n_state = n_state_min_in, n_state_max_in
         do m_state = n_state, n_state_max_in
           num = num + 1
           moment(num) = (0.0,0.0)
           ! MR OPT: NOW DO KS^T M KS on the basis -- first M*KS then dot prod
           if(newflag)then
             if(real_eigenvectors)then
               call DGEMV('N',n_basis,n_basis,1.d0, reconstmommat(:,:), n_basis, & 
                    KS_vec(:, m_state, j_spin), 1, 0.d0 ,dummy(:), 1)
!               moment(num)=dot_product(KS_vec(:, n_state, i_spin), dummy(:))! 
               moment(num)=ddot(n_basis, KS_vec(:, n_state, i_spin), 1, dummy(:), 1)
             else
               call ZGEMV('N',n_basis,n_basis,(1d0,0.d0), reconstmommatcpl(:,:), n_basis, & 
                    KS_vec_complex(:, m_state, j_spin), 1, (0.d0,0.d0) ,dummycpl(:), 1)
                    moment(num)=zdotc(n_basis, KS_vec_complex(:, n_state, i_spin), 1, dummycpl(:), 1)
             endif
             ! MR: END OPTIMIZATION
           else
             num_base = 0
             do i_basis_2 = 1, n_basis
               do i_basis_1 = 1, i_basis_2
                 num_base = num_base + 1
                 if(real_eigenvectors)then
                   if (i_basis_1==i_basis_2)then
                     moment(num) = moment(num) &
                          +  KS_vec(i_basis_1,n_state,i_spin) &
                          *mommat_full_w_up(num_base)*&
                          KS_vec(i_basis_2,m_state,j_spin)  
                   else
                     moment(num) = moment(num) &
                          + KS_vec(i_basis_1,n_state,i_spin) &
                          *(mommat_full_w_up(num_base))*&
                          KS_vec(i_basis_2,m_state,j_spin) &
                          + KS_vec(i_basis_2,n_state,i_spin)&
                          *mommat_full_w_low(num_base)*&
                          KS_vec(i_basis_1,m_state,j_spin)
                   endif
                 else
                   if (i_basis_1==i_basis_2)then
                     moment(num) = moment(num)&
                          + conjg(KS_vec_complex(i_basis_1,n_state,i_spin)) &
                          *(mommat_full_w_complex_up(num_base))*&
                          (KS_vec_complex(i_basis_2,m_state,j_spin))
                   else
                     moment(num) = moment(num)&
                          + conjg(KS_vec_complex(i_basis_1,n_state,i_spin)) &
                          *(mommat_full_w_complex_up(num_base))*&
                          (KS_vec_complex(i_basis_2,m_state,j_spin)) &
                          + ((conjg(KS_vec_complex(i_basis_2,n_state,i_spin))&
                          *conjg(mommat_full_w_complex_low(num_base))*&
                          (KS_vec_complex(i_basis_1,m_state,j_spin))))
                   endif
                 end if ! realeig
               end do ! i_basis_1
             end do ! i_basis_2
           endif ! newflag
         end do ! mstate
      end do ! nstate
    enddo ! jspin
  enddo ! ispin
end subroutine calc_moment_p0

subroutine calc_moment_p0_SOC(moment, mommat_full_w_complex_up, mommat_full_w_complex_low, &
                           KS_vec_SOC, k_point, i_coord, n_state_min_in, n_state_max_in)

  !  PURPOSE
  !   Calculete Momentum-Matrix elements by summing up of KS EV coefficents 
  !   (basis transformation) for SOC-perturbed eigenvectors
  !
  ! USES
  use dimensions,           only : n_basis
  use dimensions_soc,       only : n_basis_soc, n_saved_states_soc
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  implicit none
  !  ARGUMENTS
  integer,                                                   INTENT(IN) :: n_state_min_in
  integer,                                                   INTENT(IN) :: n_state_max_in
  complex*16, dimension(n_basis_soc, n_saved_states_soc, 1), INTENT(IN) :: KS_vec_SOC
  complex*16, dimension((n_basis+1)*n_basis/2),              INTENT(IN) :: mommat_full_w_complex_up
  complex*16, dimension((n_basis+1)*n_basis/2),              INTENT(IN) :: mommat_full_w_complex_low
  integer,                                                   INTENT(IN) :: k_point
  integer,                                                   INTENT(IN) :: i_coord
  complex*16, INTENT(OUT) ::  moment((((n_state_max_in-n_state_min_in+1)+1) &
                  *(n_state_max_in-n_state_min_in+1)/2)*2)
  !  INPUTS
  !    o KS_vec_SOC  -- Eigenvector at k_points by SOC version (always complex
  !                        in SOC)
  !    o mommat_w_up/low  -- Momentummatrix (upper/lower triangle)
  !                             in atomic basis (real)
  !    o mommat_w_complex_up/low -- Momentummatrix (upper/lower triangle) 
  !                                    in atomic basis (complex)
  !    o k_point -- k-point wanted to be calculated 
  !  OUTPUT
  !    o dipelement -- Momentummatrix in KS-basis
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !    the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Created in 20XX
  !  SOURCE
  integer    :: i_spin,j_spin,i_basis_2, i_basis_1, n_state, m_state, basis_offset
  integer    ::  num, num_base
  complex*16, dimension(:,:), allocatable :: reconstmommatcpl
  complex*16, dimension(:),   allocatable :: dummycpl
  logical :: newflag
  real*8, external     :: ddot
  complex*16, external :: zdotc
 
  call aims_allocate( reconstmommatcpl, n_basis, n_basis, "reconstmommatcpl" )
  call aims_allocate(         dummycpl, n_basis,                  "dummycpl" )

! I AM LEAVING IT FALSE (OLD WAY KEPT FOR NOW) UNTIL IT IS THOROUGHLY TESTED
  newflag = .false.
! MR OPTMIZE: build a new matrix that is actually "the" mommat matrix:
  num_base = 0
! For the reason, that the SOC engienvector is always complex
  do i_basis_2 = 1, n_basis
    do i_basis_1 = 1, i_basis_2
      num_base=num_base+1
      reconstmommatcpl(i_basis_1, i_basis_2)=mommat_full_w_complex_up(num_base)
      reconstmommatcpl(i_basis_2, i_basis_1)=conjg(mommat_full_w_complex_low(num_base))
    enddo
  enddo 

  num = 0
  moment = (0.0,0.0)
  do i_spin = 1, 2
    basis_offset = (i_spin-1)*n_basis
    do n_state = n_state_min_in, n_state_max_in
      do m_state = n_state, n_state_max_in
        num = num + 1
        !moment(num) = (0.0,0.0)
        ! MR OPT: NOW DO KS^T M KS on the basis -- first M*KS then dot prod
!        if(newflag) then
!          call ZGEMV('N',n_basis,n_basis,(1d0,0.d0), reconstmommatcpl(:,:), n_basis, & 
!                     KS_vec_SOC(:, m_state,i_spin), 1, (0.d0,0.d0) ,dummycpl(:), 1)
!                     moment(num)=zdotc(n_basis, KS_vec_SOC(:, n_state,i_spin), 1, dummycpl(:), 1)
!        ! MR: END OPTIMIZATION
!        else
          num_base = 0
          do i_basis_2 = 1, n_basis
            do i_basis_1 = 1, i_basis_2
              num_base = num_base + 1
              if (i_basis_1==i_basis_2)then
                moment(num) = moment(num)&
                     + conjg(KS_vec_SOC(i_basis_1 + basis_offset, n_state,1 )) &
                     *(mommat_full_w_complex_up(num_base))  &
                     *(KS_vec_SOC(i_basis_2 + basis_offset ,m_state, 1))
              else
                moment(num) = moment(num)&
                     + conjg(KS_vec_SOC( i_basis_1 + basis_offset, n_state, 1)) &
                     *(mommat_full_w_complex_up(num_base)) &
                     *(KS_vec_SOC(i_basis_2 + basis_offset ,m_state, 1)) &
                     + ((conjg(KS_vec_SOC(i_basis_2 + basis_offset ,n_state, 1))&
                     *conjg(mommat_full_w_complex_low(num_base)) &
                     *(KS_vec_SOC(i_basis_1 + basis_offset, m_state, 1))))
              endif
            end do ! i_basis_1
          end do ! i_basis_2
!        endif ! newflag
      end do ! mstate
    end do ! nstate
  enddo ! ispin

  call aims_deallocate( reconstmommatcpl, "reconstmommatcpl" )
  call aims_deallocate(         dummycpl,         "dummycpl" )

end subroutine calc_moment_p0_SOC

subroutine calc_dielectric_function(dipelementxi,dipelementxj,die_el,omegapl,&
                     chemical_potential,KS_eigen,k_weight,occ_numbers,&
                     use_absorption_in,widthone_in,widthtwo_in,n_state_min_in,&
                     n_state_max_in)

  !  PURPOSE
  !   Sum up Momentummatrix elements to dielectric function (Fermis Golden rule)
  !   or absoption
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  complex*16, INTENT(IN) :: dipelementxi((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2)&
                                          *((n_spin*(n_spin+1))/2))
  complex*16, INTENT(IN):: dipelementxj ((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2)&
                                          *((n_spin*(n_spin+1))/2))
  real*8, INTENT(INOUT):: die_el(n_omega)
  real*8, INTENT(INOUT):: omegapl
  real*8, INTENT(IN):: chemical_potential
  real*8, INTENT(IN):: KS_eigen(n_states, n_spin)
  real*8, INTENT(IN):: k_weight
  real*8, INTENT(IN):: occ_numbers(n_states,n_spin)
  logical, INTENT(IN):: use_absorption_in
  real*8, INTENT(IN):: widthone_in
  real*8, INTENT(IN):: widthtwo_in
! INPUTS
! o dipelementxi -- i component of momentummatrix
! o dipelementxj -- j component of momentummatrix
! o die_el -- dielectric function
! o omegapl -- plasma frequency
! o chemical_potential -- \epsilon_F
! o KS_eigen -- eigenvalues at current k_point
! o k_weight -- weight of current k_point
! o occ_numbers -- occupation numbers
! o widthone_in -- broadening of Gaussian/Lorentzian
! o widthtwo_in -- Smearing in Fermi functions
! o n_state_min_in -- Minimum state concidered
! o n_state_max_in -- Maximum state concidered
!  OUTPUT
! o die_el -- dielectric function
! o omegapl -- plasma frequency
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  real*8:: dipmult((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  real*8:: ohm((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  real*8:: dfermi((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  real*8:: omega
  real*8:: fermione
  real*8:: fermitwo
  real*8:: scaling
  real*8:: width
  real*8:: widthf
  real*8:: norm_lorentz
  real*8:: norm_gauss
  integer:: n_homo

  integer:: n_state
  integer:: m_state
  integer:: num
  integer:: omegaind
  integer:: i_state
  integer:: i_spin , j_spin

  width=widthone_in/hartree   ! In hartree now
  widthf=widthtwo_in/hartree
  norm_lorentz=(1.0/(width*pi))
  norm_gauss=sqrt(1.0/(2.0*pi))*(1.0/width)
  dipmult=0
  dfermi=0
  ohm=0
  num=0
  do i_spin = 1, n_spin
    do j_spin = i_spin, n_spin
      do n_state = n_state_min_in, n_state_max_in
	fermitwo=(1.0/(1.0+exp((KS_eigen(n_state,i_spin)-chemical_potential)/&
                                                                     widthf)))
	do m_state = n_state, n_state_max_in
	  num = num + 1
	  ohm(num)=KS_eigen(m_state,j_spin)-KS_eigen(n_state,i_spin)
	  fermione=(1.0/(1.0+exp((KS_eigen(m_state,1)-chemical_potential)&
			    /widthf)))
	  dipmult(num)=abs(dipelementxi(num)*dipelementxj(num))*k_weight*4.*pi**2/cell_volume
	  dfermi(num)=(fermitwo-fermione)  
	  if (n_state==m_state) then
		omegapl=omegapl+norm_gauss*exp(-0.5*&
		((KS_eigen(n_state,i_spin)-chemical_potential)**2/&
		((width)**2)))*dipmult(num)
	  endif 
	enddo
      enddo
    enddo
  enddo

  if (use_absorption_in) then !calculate absorption
    n_homo = 0
    do i_spin = 1, n_spin, 1  !Where is HOMO
	do i_state = 1, n_states
	  if(occ_numbers(i_state,i_spin) .gt. 1.e-12) then
	  n_homo = i_state
	  endif
      enddo
    enddo 

    do omegaind = 0, n_omega-1
        omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
        num=0
        do i_spin = 1, n_spin
          do j_spin = i_spin, n_spin
	      do n_state = n_state_min_in, n_state_max_in
		do m_state = n_state, n_state_max_in
		    num = num + 1
		    if ((n_state<=n_homo).and.(m_state>n_homo))then
			!Lorentz
			!scaling=norm_lorentz*(1.0/(1+((ohm(num)-omega)/ &
			!        (width))**2)) 
		      if (use_gauss) then
			!Gaussian
			scaling=norm_gauss*exp(-0.5*&
				(((ohm(num)-omega)**2)/((width)**2))) 
		      else
			!Lorentz
			scaling=norm_lorentz*(1.0/(1+((ohm(num)-omega)/ &
				(width))**2)) 
		      endif
		      die_el(omegaind+1)=die_el(omegaind+1)+dipmult(num)&
					    *scaling/(pi*(omega))
		    endif
		enddo
	      enddo
	  enddo
	enddo
    enddo
  else !calculate dielectric function
    do omegaind = 0, n_omega-1
        omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
        num=0
        do i_spin = 1, n_spin
          do j_spin = i_spin, n_spin
	    do n_state = n_state_min_in, n_state_max_in
	      do m_state = n_state, n_state_max_in
		  num = num + 1
		  scaling=norm_gauss*exp(-0.5*&
			      (((ohm(num)-omega)**2)/((width)**2))) 
		  die_el(omegaind+1)=die_el(omegaind+1)+dipmult(num)*&
			    dfermi(num)*scaling/((omega**2))
	      enddo
	    enddo
	  enddo
	enddo
    enddo
  endif
end subroutine calc_dielectric_function


subroutine calc_dielectric_function_p0(dipelementxi,dipelementxj,re_die_el,&
                     im_die_el, omegapl,chemical_potential,KS_eigen,&
                     k_weight,occ_numbers,&
                     use_absorption_in,widthone_in,n_state_min_in,&
                     n_state_max_in)

  !  PURPOSE
  !   Sum up Momentummatrix elements to dielectric function (Fermis Golden rule)
  !   This Version gets Imaginary and Real part directly and 
  !   calculates intra-band contribution.
  !   Calculation for absoption is unchanged
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  complex*16, INTENT(IN) :: dipelementxi((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2)&
                                          *((n_spin*(n_spin+1))/2))
  complex*16, INTENT(IN):: dipelementxj ((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2)&
                                          *((n_spin*(n_spin+1))/2))
  real*8, INTENT(INOUT):: re_die_el(n_omega)
  real*8, INTENT(INOUT):: im_die_el(n_omega)
  real*8, INTENT(INOUT):: omegapl
  real*8, INTENT(IN):: chemical_potential
  real*8, INTENT(IN):: KS_eigen(n_states, n_spin)
  real*8, INTENT(IN):: k_weight
  real*8, INTENT(IN):: occ_numbers(n_states,n_spin)
  logical, INTENT(IN):: use_absorption_in
  real*8, INTENT(IN):: widthone_in
! INPUTS
! o dipelementxi -- i component of momentummatrix
! o dipelementxj -- j component of momentummatrix
! o die_el -- dielectric function
! o omegapl -- plasma frequency (actually omegapl**2) 
! o chemical_potential -- \epsilon_F
! o KS_eigen -- eigenvalues at current k_point
! o k_weight -- weight of current k_point
! o occ_numbers -- occupation numbers
! o widthone_in -- broadening of Gaussian/Lorentzian
! o widthtwo_in -- Smearing in Fermi functions
! o n_state_min_in -- Minimum state concidered
! o n_state_max_in -- Maximum state concidered
!  OUTPUT
! o die_el -- dielectric function
! o omegapl -- plasma frequency
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  complex*16:: dipmult((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  real*8:: ohm((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  real*8:: dfermi((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2))
  real*8:: omega
  real*8:: scaling
  real*8:: width
  real*8:: norm_lorentz
  real*8:: norm_gauss
  complex*16:: sigma
  complex*16:: tau
  integer:: n_homo

  integer:: n_state
  integer:: m_state
  integer:: num
  integer:: omegaind
  integer:: i_state
  integer:: i_spin , j_spin
  integer:: occmax
  
  if (n_spin.gt.1) then
    occmax = 1
  else
    occmax = 2
  endif
  !occmax = 1 

  width=widthone_in/hartree   ! In hartree now
  norm_lorentz=(1.0/(width*pi))
  norm_gauss=sqrt(1.0/(2.0*pi))*(1.0/width)
  dipmult=0
  dfermi=0
  ohm=0
  num=0
  do i_spin = 1, n_spin
    do j_spin = i_spin, n_spin
      do n_state = n_state_min_in, n_state_max_in
	do m_state = n_state, n_state_max_in
	  num = num + 1
	  ohm(num)=KS_eigen(m_state,j_spin)-KS_eigen(n_state,i_spin)
	  dipmult(num)=dipelementxi(num)* conjg(dipelementxj(num))
	  if (n_state==m_state) then
		omegapl=omegapl+(exp(-0.5*norm_gauss*&
		((KS_eigen(n_state,i_spin)-chemical_potential)**2/&
		((width)**2)))*abs(dble (dipelementxi(num))**2 +&
		aimag(dipelementxi(num)**2)))*&
		occmax*4.0*pi/cell_volume
	  endif 
	enddo
      enddo
    enddo
  enddo
  if (use_absorption_in) then !calculate absorption
    n_homo = 0
    do i_spin = 1, n_spin, 1  !Where is HOMO
	do i_state = 1, n_states
	  if(occ_numbers(i_state,i_spin) .gt. 1.e-12) then
	  n_homo = i_state
	  endif
      enddo
    enddo 

    do omegaind = 0, n_omega-1
        omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
        num=0
        do i_spin = 1, n_spin
          do j_spin = i_spin, n_spin
	      do n_state = n_state_min_in, n_state_max_in
		do m_state = n_state, n_state_max_in
		    num = num + 1
		    if ((n_state<=n_homo).and.(m_state>n_homo))then
			!Lorentz
			!scaling=norm_lorentz*(1.0/(1+((ohm(num)-omega)/ &
			!        (width))**2)) 
		      if (use_gauss) then
			!Gaussian
			scaling=norm_gauss*exp(-0.5*&
				(((ohm(num)-omega)**2)/((width)**2))) 
		      else
			!Lorentz
			scaling=norm_lorentz*(1.0/(1+((ohm(num)-omega)/ &
				(width))**2)) 
		      endif
		      re_die_el(omegaind+1)=re_die_el(omegaind+1)+dipmult(num)&
					    *occmax*k_weight*4.*(pi**2)*scaling/&
                                            (omega*cell_volume)
		    endif
		enddo
	      enddo
	  enddo
	enddo
    enddo
  else !calculate dielectric function
    tau = dcmplx(0d0,width)
    do omegaind = 0, n_omega-1
        omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
        num=0
        do i_spin = 1, n_spin
          do j_spin = i_spin, n_spin
	    do n_state = n_state_min_in, n_state_max_in
	      do m_state = n_state, n_state_max_in
		num = num + 1
                if (KS_eigen(n_state,j_spin) <= chemical_potential .and. &
                    KS_eigen(m_state,i_spin) >= chemical_potential) then

                     sigma = (1.d0 / (ohm(num) + width)) * &
                          ((occmax * dipmult(num)) /(omega-ohm(num)+tau)+&
		       conjg(occmax * dipmult(num))/(omega+ohm(num)+tau))
		     sigma = dcmplx(0d0,1d0)* (k_weight / cell_volume) * sigma 
		     
		     re_die_el(omegaind+1) = re_die_el(omegaind+1) +  &
                                         (4.0*pi) * aimag(sigma / (omega+tau))
		     im_die_el(omegaind+1) = im_die_el(omegaind+1) +  &
                                         (4.0*pi) * dble(sigma / (omega+tau))
                endif
	      enddo
	    enddo
	  enddo
	enddo
    enddo
  endif
end subroutine calc_dielectric_function_p0


subroutine calc_dielectric_function_p1_SOC(dipelementxi,dipelementxj,re_die_el,&
                     im_die_el, omegapl,chemical_potential_SOC,KS_eigen_SOC,&
                     k_weight,occ_numbers_SOC,&
                     use_absorption_in,widthone_in,n_state_min_in,&
                     n_state_max_in)

  !  PURPOSE
  !   Sum up Momentummatrix elements to dielectric function (Fermis Golden rule)
  !   This Version gets Imaginary and Real part directly and 
  !   calculates intra-band contribution.
  !   Calculation for absoption is unchanged
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  complex*16, INTENT(IN) :: dipelementxi((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2)*2)
  complex*16, INTENT(IN):: dipelementxj ((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2)*2)
  real*8, INTENT(INOUT):: re_die_el(n_omega)
  real*8, INTENT(INOUT):: im_die_el(n_omega)
  real*8, INTENT(INOUT):: omegapl
  real*8, INTENT(IN):: chemical_potential_SOC
  real*8, INTENT(IN):: KS_eigen_SOC(2*n_states)
  real*8, INTENT(IN):: k_weight
  real*8, INTENT(IN):: occ_numbers_SOC(2*n_states)
  logical, INTENT(IN):: use_absorption_in
  real*8, INTENT(IN):: widthone_in
! INPUTS
! o dipelementxi -- i component of momentummatrix
! o dipelementxj -- j component of momentummatrix
! o die_el -- dielectric function
! o omegapl -- plasma frequency (actually omegapl**2) 
! o chemical_potential -- \epsilon_F
! o KS_eigen -- eigenvalues at current k_point
! o k_weight -- weight of current k_point
! o occ_numbers -- occupation numbers
! o widthone_in -- broadening of Gaussian/Lorentzian
! o widthtwo_in -- Smearing in Fermi functions
! o n_state_min_in -- Minimum state concidered
! o n_state_max_in -- Maximum state concidered
!  OUTPUT
! o die_el -- dielectric function
! o omegapl -- plasma frequency
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  complex*16:: dipmult((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*2)
  real*8:: ohm((((n_state_max_in-n_state_min_in+1)+1)&
                  *(n_state_max_in-n_state_min_in+1)/2)*2)
  real*8:: dfermi((((n_state_max_in-n_state_min_in+1)+1)&
                *(n_state_max_in-n_state_min_in+1)/2)*2)
  real*8:: omega
  real*8:: scaling
  real*8:: width
  real*8:: norm_lorentz
  real*8:: norm_gauss
  complex*16:: sigma
  complex*16:: tau
  integer:: n_homo

  integer:: n_state
  integer:: m_state
  integer:: num
  integer:: omegaind
  integer:: i_state
  integer:: i_spin , j_spin
  integer:: occmax
  
 ! as using the spin orbit coupling, the occmax define to 1
  occmax = 1

  width=widthone_in/hartree   ! In hartree now
  norm_lorentz=(1.0/(width*pi))
  norm_gauss=sqrt(1.0/(2.0*pi))*(1.0/width)
  dipmult=0
  dfermi=0
  ohm=0
  num=0
  
  !TZ: define the momentum matrix dipmult(num) and eigenvalue number ohm(num)
  do i_spin = 1, 2
   do n_state = n_state_min_in, n_state_max_in
	   do m_state = n_state, n_state_max_in
	      num = num + 1
	      ohm(num)=KS_eigen_SOC(m_state)-KS_eigen_SOC(n_state)
	      dipmult(num)=dipelementxi(num)* conjg(dipelementxj(num))
	      if (n_state==m_state) then
		      omegapl=omegapl+(exp(-0.5*norm_gauss*&
		      ((KS_eigen_SOC(n_state)-chemical_potential_SOC)**2/&
		      ((width)**2)))*abs(dble (dipelementxi(num))**2 +&
		      aimag(dipelementxi(num)**2)))*&
		      occmax*4.0*pi/cell_volume
	      endif 
	   enddo
   enddo
  enddo


  if (use_absorption_in) then !calculate absorption
    n_homo = 0
	 ! find the homo states 
    do i_state = 1, 2*n_states
	   if(occ_numbers_SOC(i_state) .gt. 1.e-12) then
	      n_homo = i_state
	   endif
    enddo
 

    do omegaind = 0, n_omega-1
      omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
      num=0
      do i_spin = 1, 2
	      do n_state = n_state_min_in, n_state_max_in
		      do m_state = n_state, n_state_max_in
		         num = num + 1
		         if ((n_state<=n_homo).and.(m_state>n_homo))then
			!Lorentz scaling definition
			      !   scaling=norm_lorentz*(1.0/(1+((ohm(num)-omega)/ &
			      !   (width))**2)) 
		            if (use_gauss) then
			!Gaussian
			            scaling=norm_gauss*exp(-0.5*&
				            (((ohm(num)-omega)**2)/((width)**2))) 
		            else
			!Lorentz
			            scaling=norm_lorentz*(1.0/(1+((ohm(num)-omega)/ &
				            (width))**2)) 
		            endif
		            re_die_el(omegaind+1)=re_die_el(omegaind+1)+dipmult(num)&
					      *occmax*k_weight*2.*(pi**2)*scaling/&
                     (omega*cell_volume)
		         endif
		      enddo
	      enddo
	   enddo
    enddo
   else !calculate dielectric function
    tau = dcmplx(0d0,width)
    do omegaind = 0, n_omega-1
        omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
        num=0
        do i_spin = 1, 2
	      do n_state = n_state_min_in, n_state_max_in
	         do m_state = n_state, n_state_max_in
		         num = num + 1
               if (KS_eigen_SOC(n_state) <= chemical_potential_SOC .and. &
                    KS_eigen_SOC(m_state) >= chemical_potential_SOC) then

                  sigma = (1.d0 / (ohm(num) + width))  &
                        * ((occmax * dipmult(num)) /(omega-ohm(num)+tau) &
		                  + conjg(occmax * dipmult(num))/(omega+ohm(num)+tau))
		            sigma = dcmplx(0d0,1d0)* (k_weight / cell_volume) * sigma 
		     
		            re_die_el(omegaind+1) = re_die_el(omegaind+1) +  &
                                         (4.0*pi) * aimag(sigma / (omega+tau))
		            im_die_el(omegaind+1) = im_die_el(omegaind+1) +  &
                                         (4.0*pi) * dble(sigma / (omega+tau))
               endif
	         enddo
	      enddo
      enddo
   enddo
  endif
end subroutine calc_dielectric_function_p1_SOC
!
















!!!!!!!!!!!  Output functions  !!!!!!!!!!!!!!!!!!





subroutine out_moment(moment_one_in,moment_two_in,moment_three_in, &
                         KS_Eigenvalue, k_point, n_state_min_in, n_state_max_in)

  !  PURPOSE
  !   Write Momentummatrixelements (x,y,z) to file for one k-point
  !
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  integer, INTENT(IN) :: n_state_min_in
  integer, INTENT(IN) :: n_state_max_in
  complex*16, dimension((((n_state_max_in-n_state_min_in+1)+1)*&
       (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)),&
       INTENT(IN) :: moment_one_in
  complex*16, dimension((((n_state_max_in-n_state_min_in+1)+1)*&
       (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)),&
       INTENT(IN) :: moment_two_in
  complex*16, dimension((((n_state_max_in-n_state_min_in+1)+1)*&
       (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)),&
       INTENT(IN) :: moment_three_in
  integer, INTENT(IN) :: k_point
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigenvalue
  ! Write dipelement to file
  !  INPUTS
  !    o dipelement_i -- n_states+1 x n_states/2 Momentummatrix
  !    o k_point -- k-point wanted to be calculated 
  !    o  KS_Eigenvalue
  !    o  n_state_min_in/n_state_max_in Maximum/Minimum state concidered
  !  OUTPUT
  !    o file dipelement_k_point.dat
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: n_state, m_state, i_spin, j_spin, num
  CHARACTER(len=5) :: value
  CHARACTER(len=1) :: valuek1
  CHARACTER(len=2) :: valuek2
  CHARACTER(len=3) :: valuek3
  CHARACTER(len=4) :: valuek4
  CHARACTER(len=5) :: valuek5
  CHARACTER(len=55) :: fmt
  CHARACTER(len=25) :: name
  CHARACTER(len=25) :: name2
  real*8 :: PP

  write(value,'(I3)') 7
  fmt = '(I5, I5, F10.4, I5, I5, F10.4, ' // value // 'ES14.4)'
  if (k_point<=9) then
     write(valuek1,'(I1)') k_point
     name = 'element_k_'//trim(valuek1)//'.dat'
     name2 = 'element_k_'//trim(valuek1)
  elseif (k_point<=99) then
     write(valuek2,'(I2)') k_point
     name = 'element_k_'//trim(valuek2)//'.dat'
     name2 = 'element_k_'//trim(valuek2)
  elseif (k_point<=999) then
     write(valuek3,'(I3)') k_point
     name = 'element_k_'//trim(valuek3)//'.dat'
     name2 = 'element_k_'//trim(valuek3)
  elseif (k_point<=9999) then
     write(valuek4,'(I4)') k_point
     name = 'element_k_'//trim(valuek4)//'.dat'
     name2 = 'element_k_'//trim(valuek4)
  elseif (k_point<=99999) then
     write(valuek5,'(I5)') k_point
     name = 'element_k_'//trim(valuek5)//'.dat'
     name2 = 'element_k_'//trim(valuek5)
  endif

  open(unit=9, file=name,ACTION='WRITE')
     WRITE (9,120) k_point, k_point_list(k_point,1:3)
     120 FORMAT ('# K-point:',I5,' at', 3F10.6)
     if (flag_out_dipmat .or. flag_out_dipmat_k_k) then
        WRITE (9,121)
        121 FORMAT ('# KS state i, spin i, KS energy i [eV], KS state j, &
                 & spin j, KS energy j &
                 & [eV], Re(x_ij), Im(x_ij), Re(y_ij), Im(y_ij), &
                 & Re(z_ij), Im(z_ij), (x_ij^2+y_ij^2+z_ij^2)/3')
     else
         WRITE (9,122)
         122 FORMAT ('# KS state i, spin i, KS energy i [eV], KS state j, &
                 & spin j, KS energy j &
                 & [eV], Re(p^x_ij), Im(p^x_ij), Re(p^y_ij), Im(p^y_ij), &
                 & Re(p^z_ij), Im(p^z_ij), (p^x_ij^2+p^y_ij^2+p^z_ij^2)/3')
     endif
     num = 0
     do i_spin = 1, n_spin
       do j_spin = i_spin, n_spin
	  do n_state = n_state_min_in,n_state_max_in
	    do m_state = n_state, n_state_max_in
		num = num + 1
		PP = dble((moment_one_in(num)*conjg(moment_one_in(num))+&
			  moment_two_in(num)*conjg(moment_two_in(num))+&
			  moment_three_in(num)*conjg(moment_three_in(num)))) 
		WRITE (9,fmt) n_state, i_spin,&
                              KS_eigenvalue(n_state, i_spin, k_point)*hartree, &
	                      m_state, j_spin, &
                              KS_eigenvalue(m_state, j_spin, k_point)*hartree, &
			   Real(moment_one_in(num)), aImag(moment_one_in(num)),&
			   Real(moment_two_in(num)), aImag(moment_two_in(num)),&
		       Real(moment_three_in(num)), aImag(moment_three_in(num)),&
			      (PP)/3.
	    end do
	  end do
       enddo
     enddo
  close(unit=9)

end subroutine out_moment

subroutine out_dielectric_function(die_el,ep1_in,ep2_in)


  !  PURPOSE
  !   Write Dielectric function \epsilon_ep1_ep2 to file
  !
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  real*8, INTENT(IN) :: die_el(n_omega)
  CHARACTER(len=1), INTENT(IN):: ep1_in
  CHARACTER(len=1), INTENT(IN):: ep2_in
  ! Write dipelement to file
  !  INPUTS
  !    o die_el -- n_omega Dielectric function
  !  OUTPUT
  !    o file
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: num
  integer:: omega
  CHARACTER(len=55) :: fmt
  CHARACTER(len=25) :: name


  fmt = '(2ES17.4)'
  name = "dielectric_"//trim(ep1_in)//"_"//trim(ep2_in)//".out"
  num=0
  open(unit=8, file=name,ACTION='WRITE')
  WRITE (8,122) trim(ep1_in), trim(ep2_in)
  122 FORMAT ('# \omega [eV], Im(\epsilon_{',A1,'_',A1,'})')
     do omega=0,n_omega-1
         num=num+1
         if (abs((die_el(num)))>1.0E-30)then
         WRITE (8,fmt) omega_min+omega*((omega_max-omega_min)/n_omega), &
                       die_el(num)
         else
         WRITE (8,fmt) omega_min+omega*((omega_max-omega_min)/n_omega), 0.0
         endif
     end do
  close(unit=8)

end subroutine out_dielectric_function






subroutine out_dielectric_function_p0(re_die_el,im_die_el,omegapl,ep1_in,ep2_in,&
                                      widthone_in)


  !  PURPOSE
  !   Write Dielectric function \epsilon_ep1_ep2 to file, Intra and Inter-band
  !   contributions, Real and Imaginary part.
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  real*8, INTENT(IN) :: re_die_el(n_omega)
  real*8, INTENT(IN) :: im_die_el(n_omega)
  real*8, INTENT(IN) :: omegapl
  CHARACTER(len=1), INTENT(IN):: ep1_in
  CHARACTER(len=1), INTENT(IN):: ep2_in
  real*8, INTENT(IN):: widthone_in
  ! Write dipelement to file
  !  INPUTS
  !    o die_el -- n_omega Dielectric function
  !  OUTPUT
  !    o file
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE
  
  complex*16 :: sigma
  complex*16 :: tau  
  real*8 :: re_intra_die_el
  real*8 :: im_intra_die_el
  integer:: num
  integer:: omegaind
  real*8 :: omega
  real*8 :: width
  CHARACTER(len=55) :: fmt
  CHARACTER(len=25) :: name


  fmt = '(7ES17.4)'
  name = "dielectric_"//trim(ep1_in)//"_"//trim(ep2_in)//".out"
  num=0
  open(unit=8, file=name,ACTION='WRITE')
  WRITE (8,522) trim(ep1_in), trim(ep2_in),trim(ep1_in), trim(ep2_in),&
                &trim(ep1_in), trim(ep2_in),trim(ep1_in), trim(ep2_in)
  522 FORMAT ('# \omega [eV],  Re(\epsilon_{',A1,'_',A1,'}^{inter}), '&
              &'Im(\epsilon_{',A1,'_',A1,'}^{inter}),  '&
              &'Re(\epsilon_{',A1,'_',A1,'}^{intra}), '&
              &'Im(\epsilon_{',A1,'_',A1,'}^{intra}), \epsilon^{inter}+'&
              &'\epsilon^{intra}')
  WRITE (8,423) sqrt(abs(omegapl))*hartree
  423 FORMAT ('# \omega_pl [eV]:  ',ES17.4)
  
     width=widthone_in/hartree
     tau = dcmplx(0d0,width)
     do omegaind=0,n_omega-1

         omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree
         sigma = (dcmplx(((sqrt(abs(omegapl))**2  ) / (4d0 * pi) ),0d0))/ &
                 & (dcmplx(width,0d0) - dcmplx(0d0,omega))
         re_intra_die_el= &
                 & (4.0*pi) * aimag(sigma / (dcmplx(omega,0d0)+tau))
         im_intra_die_el= &
                 & (4.0*pi) * dble(sigma / (dcmplx(omega,0d0)+tau))
         num=num+1
         if (trim(ep1_in) .ne. trim(ep2_in))then
           WRITE (8,fmt) omega_min+omegaind*((omega_max-omega_min)/n_omega), &
                       -re_die_el(num), im_die_el(num), -re_intra_die_el, &
                       im_intra_die_el,&
                       -re_die_el(num)-re_intra_die_el, &
                       im_die_el(num)+im_intra_die_el
         else
           WRITE (8,fmt) omega_min+omegaind*((omega_max-omega_min)/n_omega), &
                       1.0-re_die_el(num), im_die_el(num), &
                       -re_intra_die_el, im_intra_die_el,&
                       1.0-re_die_el(num)-re_intra_die_el, &
                       im_die_el(num)+im_intra_die_el           
         endif
     end do
  close(unit=8)

end subroutine out_dielectric_function_p0
subroutine out_absorption(die_el,ep1_in)


  !  PURPOSE
  !   Write Absorption \alpha_ep1 to file
  !
  ! USES

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  real*8, INTENT(IN) :: die_el(n_omega)
  CHARACTER(len=1), INTENT(IN):: ep1_in
  ! Write dipelement to file
  !  INPUTS
  !    o die_el -- n_omega Dielectric function/Absorptiom
  !  OUTPUT
  !    o file
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: num
  integer:: omega
  CHARACTER(len=55) :: fmt
  CHARACTER(len=16) :: name


  fmt = '(2ES17.4)'
  name = "absorption_"//trim(ep1_in)//".out"
  num=0
  open(unit=8, file=name,ACTION='WRITE')
  WRITE (8,122) trim(ep1_in)
  122 FORMAT ('# \omega [eV], \alpha_{',A1,'}')
     do omega=0,n_omega-1
         num=num+1
         if (abs((die_el(num)))>1.0E-30)then
         WRITE (8,fmt) omega_min+omega*((omega_max-omega_min)/n_omega), &
                       die_el(num)
         else
         WRITE (8,fmt) omega_min+omega*((omega_max-omega_min)/n_omega), 0.0
         endif
     end do
  close(unit=8)

end subroutine out_absorption

subroutine calc_greenwood(dipelementxi,dipelementxj,dielectric_function,seebeck,abtew_cond, & 
                           abtew_seebeck, fermideriv, omegapl, chemical_potential,KS_eigen,k_weight, &
                             widthone_in, widthtwo_in,n_state_min_in, n_state_max_in)

  !  PURPOSE
  !   Sum up Momentummatrix elements to dielectric function, 
  !   optical conductivity, Seebeck coefficient (Fermis Golden rule)
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  real*8:: omega,omegaev,energie
  real*8:: k_weight
  real*8:: widthone_in
  real*8:: widthtwo_in


  real*8, INTENT(INOUT):: dielectric_function(n_omega)
  real*8, INTENT(INOUT):: seebeck(n_omega)
  real*8, INTENT(INOUT):: abtew_cond (n_greenenergy)
  real*8, INTENT(INOUT):: abtew_seebeck (n_greenenergy)
  real*8, INTENT(INOUT):: omegapl
  real*8, INTENT(IN):: fermideriv (n_greenenergy)
  real*8, INTENT(IN):: chemical_potential


  real*8:: Emin_ha, Emax_ha
  real*8:: norm_gauss
  integer n_state_min_in
  integer n_state_max_in
  complex*16, INTENT(IN) :: dipelementxi((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2))
  complex*16, INTENT(IN):: dipelementxj ((((n_state_max_in-n_state_min_in+1)+1)&
                                          *(n_state_max_in-n_state_min_in+1)/2))
  real*8, INTENT(IN):: KS_eigen(n_states, n_spin)
!  complex*16:: dipelementxi ((n_states+1)*n_states/2)
!  complex*16:: dipelementxj ((n_states+1)*n_states/2)
!  real*8:: dipmultfeld ((n_states+1)*n_states/2)

  real*8:: dipmult((((n_state_max_in-n_state_min_in+1)+1)&
                    *(n_state_max_in-n_state_min_in+1)/2))
  real*8:: ohm((((n_state_max_in-n_state_min_in+1)+1)&
                    *(n_state_max_in-n_state_min_in+1)/2))
  real*8:: dfermi((((n_state_max_in-n_state_min_in+1)+1)&
                    *(n_state_max_in-n_state_min_in+1)/2))


  real*8:: scaling, scaling_one, scaling_two
  real*8:: width,widthf
  real*8:: fermione
  real*8:: fermitwo
  
  integer:: n_state,m_state,num,omegaind, energieind

!  CHARACTER(len=85) :: fmt  ! only for temporary output!!!!
!  CHARACTER(len=25) :: name
    width=widthone_in/hartree
    widthf=widthtwo_in/hartree
    Emin_ha=Emin/hartree
    Emax_ha=Emax/hartree
    norm_gauss=sqrt(1.0/(2.0*pi))*(1.0/width)

  dipmult=0.0d0
  dfermi=0.0d0
  ohm=0.0d0
  num=0

  do n_state = n_state_min_in, n_state_max_in
    fermitwo=(1.0/(1.0+exp((KS_eigen(n_state,1)-chemical_potential)/widthf)))
    do m_state = n_state, n_state_max_in
      num = num + 1
      ohm(num)=KS_eigen(m_state,1)-KS_eigen(n_state,1)
      fermione=(1.0/(1.0+exp((KS_eigen(m_state,1)-chemical_potential)&
                         /widthf)))
      dipmult(num)=abs(dipelementxi(num)*dipelementxj(num))
!      dipmult(num)=abs(dipelementxi(num)*dipelementxj(num))*k_weight/pi
      dfermi(num)=(fermitwo-fermione)  
      if (n_state==m_state) then
             omegapl=omegapl+norm_gauss*exp(-0.5*((KS_eigen(n_state,1) - &
                chemical_potential)**2/((width)**2)))*dipmult(num)*k_weight/pi ! k-weight and pi added...
      endif 
    enddo
  enddo

! Begin of block that calculates the omega-->0 ("direct current") value of the el.transport coeffs. 
! triggered by "flag_out_dclimit" (see Abtew et al. PRB 76 (2007)) 
if (flag_out_dclimit) then

! region-restriction wrt Fermi-deriv can be set up already here. 
! relevant_energieind_min=   ONLY two numbers need to be know here to directly 
! enter the following do loop at these values and not scanning the entire range. 

   do energieind=0, n_greenenergy-1
      energie=Emin_ha+energieind*((Emax_ha - Emin_ha)/n_greenenergy)
      num=0
    if(fermideriv(energieind+1)<1E-30) then
      abtew_cond(energieind+1)=0.0d0
    else
      do n_state = n_state_min_in, n_state_max_in
			scaling_one=norm_gauss*exp(-0.5*((( KS_eigen(n_state,1)-energie  )**2)/((width)**2)))  
         do m_state = n_state, n_state_max_in
            num = num + 1
		  if (n_state/=m_state) then
      		scaling_two = norm_gauss &
               * exp(-0.5*((( KS_eigen(m_state,1)-energie  )**2)/((width)**2)))
            abtew_cond(energieind+1) = abtew_cond(energieind+1) + k_weight &
               *dipmult(num)*scaling_one*scaling_two*fermideriv(energieind+1)
            abtew_seebeck(energieind+1) = abtew_seebeck(energieind+1) &
               + k_weight*dipmult(num) * scaling_one*scaling_two &
               * fermideriv(energieind+1) &
               * ((KS_eigen(n_state,1) + KS_eigen(m_state,1))/2 &
                   - chemical_potential)
		  endif ! n!=m condition
	    enddo ! end m-state
	  enddo ! end n-state
    endif ! end of fermidervi>1e-30 condition
    enddo ! end of energy loop
endif ! end of dc-limit block



    do omegaind = 0, n_omega-1
        omegaev=omega_min+omegaind*((omega_max-omega_min)/n_omega)
        num=0
        omega = omegaev/hartree
	  do n_state = n_state_min_in, n_state_max_in
	    do m_state = n_state, n_state_max_in
	      num = num + 1

           if(n_state/=m_state) then 
              scaling=norm_gauss*exp(-0.5*(((ohm(num)-omega)**2)/((width)**2)))
              dielectric_function(omegaind+1) = dielectric_function(omegaind+1)+k_weight*dipmult(num)*dfermi(num)*scaling     
!/(pi*(omega**2)) has been outcommted part from björn's original
              seebeck(omegaind+1)=seebeck(omegaind+1)+k_weight*dipmult(num)*dfermi(num)*scaling * &
                                  ((KS_eigen(n_state,1)+KS_eigen(m_state,1))/2-chemical_potential) 
! last seebeck term can be handled as ohm or dfermi as well!!!
          endif
	    enddo
	  enddo
    enddo


end subroutine calc_greenwood

subroutine out_greenwood(die_el,seebeck,abtew, abtewseebeck, fermideriv, ep1_in,ep2_in, zellgr, chempot)
  !  PURPOSE
  !   Write Dielectric function \epsilon_ep1_ep2 to file
  !
  ! USES

  use constants, only: pi
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none

  !  ARGUMENTS
  real*8, INTENT(IN) :: die_el(n_omega)
  real*8, INTENT(IN) :: seebeck(n_omega)
  real*8, INTENT(IN) :: abtew(n_greenenergy)
  real*8, INTENT(IN) :: abtewseebeck(n_greenenergy)
  real*8, INTENT(IN) :: fermideriv(n_greenenergy)

  CHARACTER(len=1), INTENT(IN):: ep1_in
  CHARACTER(len=1), INTENT(IN):: ep2_in
  real*8:: zellgr
  real*8:: chempot

  ! Write dipelement to file
  !  INPUTS
  !    o die_el -- n_omega Dielectric function
  !  OUTPUT
  !    o file
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: num
  integer:: omega
  real*8:: omega_value, energie_value, sigma, abtewvalue
  real*8:: abtewseebeckvalue
  real*8:: prefactorL11
  real*8:: prefactorL12

  CHARACTER(len=55) :: fmt1
  CHARACTER(len=55) :: fmt2
  CHARACTER(len=55) :: fmt3
  CHARACTER(len=45) :: nameOmega
  CHARACTER(len=55) :: nameDClimit
  CHARACTER(len=10) :: fermi
  ! file indices should be defined here already

  write(fermi, '(f7.4)')(chempot*hartree)


  nameOmega = "dielectric_sparse_matrix_"//trim(ep1_in)//"_"//trim(ep2_in)//"_"//trim(fermi)//".out"
  nameDClimit= "DC_limit_transport_sparse_matrix_"//trim(ep1_in)//"_"//trim(ep2_in)//"_"//trim(fermi)//".out"    


! equals (2pi * e*e * hbar )/(1 * CellVol * mass * 100(cm) ) - 
! EXPLICIT prefactorL11=(2*pi * ((-1.602176565E-19)**2) * (4.135667E-15/(2*pi)) * 9.109382E-31) / 
! (100 * 1*(zellgr * bohr**3)*1E-30*(9.109382E-31)**2) ! is (Coulumb**2 (eV*s) * kg )/(meter**3 * kg**2 ), 100 is for m->cm
 prefactorL11= 1.16540652466953E6 * (1/(1*(zellgr * bohr**3))) 

! EXPLICIT prefactorL12= (2*pi * -1.602176565E-19 * (4.135667E-15/(2*pi)) (hbar in eV*s - as forthcoming omega-division is in eV as well) 
! *4.35974434E-18 (Ha-in-Joule) * 9.109382E-31 )/(1 * (zellgr * bohr**3) * 1E-30*(9.109382E-31)**2))
 prefactorL12 = -3.17123256619787E9 *1E6 * 0.01 * (1/(1*(zellgr * bohr**3)))  

! pref2 equals (2pi * e * hbar * Hartree-to-Joule ) / (1 * CellVol * mass  ), factor 1e6 to give muV, factor 0.01 to 
! correct the Ohm*cm unit of sigma that enters S (not L12) , negative sign is for single elementary charge


if(flag_out_dclimit) then

  fmt1 = '(5ES17.4)'
  open(unit=8, file=nameOmega,ACTION='WRITE')
  WRITE (8,122) trim(ep1_in), trim(ep2_in)
  122 FORMAT ('# omega (eV), Im(\epsilon_{',A1,'_',A1,&
     '})   L11-Conductivity (Ohm-1*cm-1)       L12      (L12/(T*L11))', &
     '-Seebeck(muV/K)')
  num=0
     do omega=0,n_omega-1
         num=num+1
         omega_value=omega_min+omega*((omega_max-omega_min)/n_omega)
         energie_value=Emin+omega*((Emax-Emin)/n_omega)
         sigma=prefactorL11*die_el(num)/omega_value
             if ( sigma < 1.0E-30)then
               sigma=0.0d0
             endif
         WRITE (8,fmt1) omega_value, die_el(num)/(pi*((omega_value/hartree)*(omega_value/hartree ))), sigma, &
               (prefactorL12*seebeck(num)/omega_value) ,(prefactorL12*seebeck(num)/omega_value)/ & 
               (11604.519*widthtwo* prefactorL11*die_el(num)/omega_value)

!        WRITE (8,fmt1) omega_value, die_el(num)/(pi*(omega_value**2)), sigma, (prefactorL12*seebeck(num)/omega_value), &
!           (prefactorL12*seebeck(num)/omega_value)/(11604.519*widthtwo* prefactorL11*die_el(num)/omega_value) 
!            this output differs in hartree or eV for dielectric function

     end do
  close(unit=8)

  fmt2 = '(4ES17.4)'
  open(unit=9, file=nameDClimit,ACTION='WRITE')
  WRITE (9,125)
  125 FORMAT ('# Energy (eV)     AbtewSigmaValue   AbtewSeebeckValue   Fermideriv')
  num=0
     do omega=0,n_greenenergy-1
         num=num+1
         energie_value=Emin+omega*((Emax-Emin)/n_greenenergy)
         abtewvalue=abtew(num)
         abtewseebeckvalue=abtewseebeck(num)
             if (abtewvalue<1.0E-30)then
               abtewvalue=0.0d0
             endif   
             if (abs(abtewseebeckvalue) < 1.0E-30)then
               abtewseebeckvalue=0.0d0
             endif   
         WRITE (9,fmt2) energie_value, abtewvalue, abtewseebeckvalue, fermideriv(num)
     end do
  close(unit=9)


else ! case where no d.c.-limit has been asked for

  fmt3 = '(5ES17.4)'
  open(unit=8, file=nameOmega,ACTION='WRITE')
  WRITE (8,124) trim(ep1_in), trim(ep2_in)
  124 FORMAT ('# omega (eV), Im(\epsilon_{',A1,'_',A1,&
     '})   L11-Conductivity (Ohm-1*cm-1)       L12      (L12/(T*L11))-', &
     'Seebeck(muV/K)')
  num=0
     do omega=0,n_omega-1
         num=num+1
         omega_value=omega_min+omega*((omega_max-omega_min)/n_omega)
         sigma=prefactorL11*die_el(num)/omega_value
          if ( sigma<1.0E-30)then
            sigma=0.0d0
          endif
         WRITE (8,fmt3) omega_value, die_el(num)/(pi*((omega_value/hartree)*(omega_value/hartree ))), sigma,&
                (prefactorL12*seebeck(num)/omega_value) ,(prefactorL12*seebeck(num)/omega_value)/ &
                (11604.519*widthtwo* prefactorL11*die_el(num)/omega_value)
!        WRITE (8,fmt3) omega_value, die_el(num)/(pi*(omega_value**2)), sigma, &
!            (prefactorL12*seebeck(num)/omega_value) ,(prefactorL12*seebeck(num)/omega_value)/(11604.519*widthtwo* &
!             prefactorL11*die_el(num)/omega_value) ! differs in hartree or eV for dielectric function
     end do
  close(unit=8)


endif


end subroutine out_greenwood

function ergebnis(abtew, abtewseebeck, zellgr, switch)
!  use pbc_lists
!  use geometry
!  use dimensions
  use runtime_choices
!  use basis
  implicit none

  real*8:: abtew(n_greenenergy)
  real*8:: abtewseebeck(n_greenenergy)
  real*8:: summe_cond, summe_seeb
  real*8:: Emin_ha, Emax_ha, ergebnis, prefactorAbtewCond, prefactorAbtewSeeb, zellgr
  integer :: i, maxcounter, switch


    Emin_ha=Emin/hartree
    Emax_ha=Emax/hartree
    summe_cond=0.0d0
    summe_seeb=0.0d0   
    maxcounter=n_greenenergy/3 ! for more detailed simpsons rule; modulo abfrage in read_control. 

select case (switch)
  case(1) ! Conducivity Integral

!  prefactorAbtew=(2*pi*((-1.602176565E-19)**2) * (1.519829877E-16/(2*pi)) * 9.109382E-31)/&
! (100 * 1 * (zellgr * bohr**3)*1E-30 * (9.109382E-31)**2  ) 
! (e**2 * hbar * m) / (V * m**2) is in units (Coulumb**2 (Ha*s) * kg )/(meter**3 * kg**2 ) 
! upper kg has to be added to due to absence of electronic mass in the momentum matrix-elements, 
!hbar-Hatree is compatible with 1/Ha unit of the summation

    prefactorAbtewCond = 42827.9032811 * (1/(zellgr * bohr**3))   ! is the conductivity prefactor. 

   do i=1, maxcounter
     summe_cond= summe_cond + abtew(3*i) + 3*abtew(3*i-1) + 3*abtew(3*i-2) + abtew(3*i-3) ! gives 1/Ha unit
   enddo
   ergebnis=prefactorAbtewCond * 0.375d0*((Emax_ha-Emin_ha)/n_greenenergy) * summe_cond ! Simpson Rule. has Hartree dimension!

  case(2) ! Seebeck Integral

!  prefactorAbtew=(2*pi*((-1.602176565E-19)**2) * (1.519829877E-16/(2*pi)) * 9.109382E-31)  &
! /  (100 * 1 * (zellgr * bohr**3)*1E-30 * (9.109382E-31)**2  ) 
! (e**2 * hbar * m) / (V * m**2) is (Coulumb**2 (Ha*s) * kg )/(meter**3 * kg**2 ) upper kg has 
! to be added to due to absence of electronic mass in the momentum matrix-elements
! prefactorAbtew-seebeck= (2pi * (-1.602e-19)(A*s) * (6.62606957E-34/2pi)(V*A*s*s)  * (9.109382E-31)kg )&
! / (1 * (zellgr*bohr**3)*1E-30 * (9.109382E-31)**2(kg**2)   ) 

    prefactorAbtewSeeb = -1.165406543E8 * (1/(zellgr * bohr**3))

 
   do i=1, maxcounter
     summe_seeb= summe_seeb + abtewseebeck(3*i) + 3*abtewseebeck(3*i-1) + 3*abtewseebeck(3*i-2) + abtewseebeck(3*i-3) 
   enddo
                     
    ergebnis=prefactorAbtewSeeb * 0.375d0*((Emax_ha-Emin_ha)/n_greenenergy) * summe_seeb
endselect
end function

subroutine calc_fermideriv(fermideriv, chemical_potential)

  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none
  integer :: i
  real*8:: Emin_ha, Emax_ha, energie, chemical_potential, widthf
  real*8:: fermideriv(n_greenenergy)

   widthf=widthtwo/hartree
   Emin_ha=Emin/hartree
   Emax_ha=Emax/hartree

  do i=1, n_greenenergy

   energie=Emin_ha+i*((Emax_ha - Emin_ha)/n_greenenergy)
   fermideriv(i)=(exp((energie+chemical_potential)/widthf))/(widthf* &
                 ((exp(chemical_potential/widthf)+exp(energie/widthf))**2))

  enddo 

end subroutine calc_fermideriv


end module calculate_mommat_base
