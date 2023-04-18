!****h* FHI-aims/hartree_potential_storage
!  NAME
!    hartree_potential_storage - provides permanent storage for hartree potential calculations
!  SYNOPSIS
module hartree_potential_storage
  !  PURPOSE
  !
  !    Provides permanent storage for rho_multipole etc as well as allocation/deallocation
  !    routines and a routine for retrieving the rho_multipole spline coefficients
  !
  !    Maybe later on update_hartree_potential_p1/sum_up_whole_potential_p1 should also go
  !    to this module (renaming it to hartree_potential or so)
  !
  !  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use spline
  use localorb_io
  use pbc_lists, only : centers_hartree_potential, center_to_atom
  use synchronize_mpi, only : sync_find_max
  use mpi_utilities
  use, INTRINSIC ::ISO_C_BINDING

  implicit none

  !  AUTHOR
  !    FHI-aims team.
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  INPUTS
  !    none
  !  OUTPUT
  !    none
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
  !******

  private ! all module data is private, especially all data imported by use statements above!

  ! Number of atoms stored in rho_multipole.
  ! This is different on every task (unless storage of all atoms is forced)
  integer, public :: n_rho_multipole_atoms

  ! List of atoms stored in rho_multipole:
  integer, allocatable, public :: i_rho_multipole_atoms(:)

  ! For every atom: index in rho_multipole or 0 if not stored in rho_multipole
  integer, allocatable, public :: rho_multipole_index(:)

  ! rho_multipole
  real*8, allocatable, public :: rho_multipole(:,:,:)
  real*8, pointer, public :: rho_multipole_shmem(:,:,:)
  type(C_PTR), public :: rho_multipole_shmem_ptr
  integer, public :: rho_multipole_shmem_win
  logical, public :: use_rho_multipole_shmem = .true.
  integer, public :: rho_multipole_tile_size = 2048

  real*8, allocatable, public :: rho_multipole_supercell(:,:,:)

  real*8, allocatable, public :: delta_v_hartree_part_at_zero_supercell(:)
  real*8, allocatable, public :: delta_v_hartree_deriv_l0_at_zero_supercell(:,:)


  ! multipole moments on the original "radial" integration grid
  ! - i.e., before any spline errors could be introduced by splining to 
  !   to the denser logarithmic grid.
  real*8, dimension(:,:),     allocatable, public :: original_multipole_moments

  ! variables for compensating density to offset residual multipole (charge only) components
  ! if requested
  real*8, dimension(:), allocatable, public :: compensation_norm
  real*8, dimension(:), allocatable, public :: compensation_radius
  real*8, public :: total_compensated_charge

  ! public routines
  public :: initialize_hartree_potential_storage
  public :: reset_hartree_potential_storage
  public :: cleanup_hartree_potential_storage
  public :: get_rho_multipole_spl
  public :: get_rho_multipole_supercell_spl
  public :: get_multipole_moments_on_original_grid

  public :: compensating_density

contains

  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/initialize_hartree_potential_storage
  !  NAME
  !    initialize_hartree_potential_storage
  !  SYNOPSIS

  subroutine initialize_hartree_potential_storage()

    !  PURPOSE
    !    Prepare module data of hartree_potential_storage.
    !  USES

    use mpi_tasks, only : check_allocation
    use physics, only : partition_tab

    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


    integer istat, mpierr, win_tmp, winsize, disp_unit, tmp_int
    integer rho_multipole_shape(3)
    integer(kind = MPI_ADDRESS_KIND) size_bytes, size_bytes_tmp

    if (myid .eq. 0) print*, "Enter initialize_hartree_potential_storage"
    call localorb_info('  Initialize hartree_potential_storage')
    ! Allocate and set atom lists

    allocate(i_rho_multipole_atoms(n_atoms),stat=istat)
    call check_allocation(istat, 'i_rho_multipole_atoms')

    allocate(rho_multipole_index(n_atoms),stat=istat)
    call check_allocation(istat, 'rho_multipole_index')

    n_rho_multipole_atoms = 0

    call set_hartree_potential_atom_lists(n_my_batches, batches, partition_tab)

    ! Finally allocate rho_multipole to the needed size

    if (use_rho_multipole_shmem) then
      allocate(rho_multipole(1, 1, 1),stat=istat)
    else
      allocate(rho_multipole((l_pot_max+1)**2, n_max_radial+2, n_rho_multipole_atoms),stat=istat)
    endif
    call check_allocation(istat, 'rho_multipole')
    rho_multipole = 0

    ! ----------------------------------------------------------------

    if (use_rho_multipole_shmem) then

      rho_multipole_shape(1) = (l_pot_max+1)**2
      rho_multipole_shape(2) = n_max_radial+2
      rho_multipole_shape(3) = n_rho_multipole_atoms
      size_bytes = rho_multipole_shape(1) * rho_multipole_shape(2) * rho_multipole_shape(3)
      size_bytes_tmp = size_bytes

      call MPI_Win_allocate_shared(size_bytes_tmp, 8, MPI_INFO_NULL, shm_comm, rho_multipole_shmem_ptr, rho_multipole_shmem_win, mpierr)
      call check_allocation(mpierr, 'MPI_Win_allocate_shared:rho_multipole_shmem_ptr')

      win_tmp = rho_multipole_shmem_win
      call MPI_Win_shared_query(win_tmp, 0, winsize, disp_unit, rho_multipole_shmem_ptr, mpierr)
      if (winsize /= size_bytes) then
        tmp_int = size_bytes
        print*, "disp_unit:", disp_unit, "winsize:", winsize, ", size_bytes:", size_bytes
        call check_allocation(winsize-tmp_int, 'MPI_Win_shared_query:winsize-size_bytes')
      end if
      call check_allocation(mpierr, 'MPI_Win_shared_query:rho_multipole_shmem_ptr')

      call C_F_pointer(rho_multipole_shmem_ptr, rho_multipole_shmem, rho_multipole_shape)

      if (shm_rank .eq. 0) then
        call MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, rho_multipole_shmem_win, mpierr)
        call MPI_Win_sync(rho_multipole_shmem_win, mpierr)
        rho_multipole_shmem = 0
        call MPI_Win_unlock(0, rho_multipole_shmem_win, mpierr)
      end if

    endif

    ! ----------------------------------------------------------------

    if (packed_matrix_format /= PM_none .and. n_periodic > 0.and.use_DFPT_phonon.or.use_friction) then
     allocate(rho_multipole_supercell((l_pot_max+1)**2, n_max_radial+2, n_centers_in_sc_DFPT),stat=istat)
     call check_allocation(istat, 'rho_multipole_supercell')
     rho_multipole_supercell = 0.0d0
     
     allocate(delta_v_hartree_part_at_zero_supercell(n_centers_in_sc_DFPT),stat=istat)
     call check_allocation(istat, 'delta_v_hartree_part_at_zero  ')
     delta_v_hartree_part_at_zero_supercell = 0.0d0

     allocate(delta_v_hartree_deriv_l0_at_zero_supercell(3,n_centers_in_sc_DFPT),stat=istat)
     call check_allocation(istat, 'delta_v_hartree_deriv_l0_at_zero_supercell')
     delta_v_hartree_deriv_l0_at_zero_supercell = 0.0d0
    endif

    if (.not.allocated(original_multipole_moments)) then
      allocate( original_multipole_moments( ( l_pot_max + 1)**2, n_atoms),stat=istat)
      call check_allocation(istat, 'original_multipole_moments    ')
    end if

    if (compensate_multipole_errors) then
      if (.not.allocated(compensation_norm)) then
        allocate( compensation_norm( n_atoms),stat=istat)
        call check_allocation(istat, 'compensation_norm    ')
      end if

      if (.not.allocated(compensation_radius)) then
        allocate( compensation_radius( n_atoms),stat=istat)
        call check_allocation(istat, 'compensation_radius    ')
      end if
    end if

    if (myid .eq. 0) print*, "Return initialize_hartree_potential_storage"

  end subroutine initialize_hartree_potential_storage
  !******
  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/reset_hartree_potential_storage
  !  NAME
  !    reset_hartree_potential_storage
  !  SYNOPSIS

  subroutine reset_hartree_potential_storage(n_my_batches_work, batches_work, partition_tab_work)

    !  PURPOSE
    !    Reinitializes the data of hartree_potential_storage if load balancing is in effect.
    !  USES

    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(in) :: n_my_batches_work
    type (batch_of_points), intent(in) :: batches_work(n_my_batches_work)
    real*8, intent(in) :: partition_tab_work(*)

    !  INPUTS
    !    n_my_batches_work -- number of batches for current distribution
    !    batches_work -- batches for current distribution
    !    partition_tab_work -- partition_tab_work for current distribution
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


    integer istat, mpierr, win_tmp, winsize, disp_unit, tmp_int
    integer rho_multipole_shape(3)
    integer(kind = MPI_ADDRESS_KIND) size_bytes, size_bytes_tmp


    call localorb_info('  Reset hartree_potential_storage for load balanced distribution.')

    ! If use_distributed_spline_storage is not set, there is nothing to do.
    ! Especially we leave rho_multipole as is if it should be needed for postprocessing ...

    if(.not.use_distributed_spline_storage) return

    ! Set new atom lists

    call set_hartree_potential_atom_lists(n_my_batches_work, batches_work, partition_tab_work)

    ! deallocate/reallocate rho_multipole
    ! RJ: Is there any reason to keep the old values (for postprocessing or so) ???

    deallocate(rho_multipole)
    ! call MPI_Win_free(rho_multipole_shmem_win, istat)
    ! call check_allocation(istat, 'MPI_Win_free:rho_multipole_shmem_win')

    if (use_rho_multipole_shmem) then
      allocate(rho_multipole(1, 1, 1),stat=istat)
    else
      allocate(rho_multipole((l_pot_max+1)**2, n_max_radial+2, n_rho_multipole_atoms),stat=istat)
    endif
    call check_allocation(istat, 'rho_multipole')
    rho_multipole = 0

    ! ----------------------------------------------------------------

    if (use_rho_multipole_shmem) then
      rho_multipole_shape(1) = (l_pot_max+1)**2
      rho_multipole_shape(2) = n_max_radial+2
      rho_multipole_shape(3) = n_rho_multipole_atoms
      size_bytes = rho_multipole_shape(1) * rho_multipole_shape(2) * rho_multipole_shape(3)
      size_bytes_tmp = size_bytes

      call MPI_Win_allocate_shared(size_bytes_tmp, 8, MPI_INFO_NULL, shm_comm, rho_multipole_shmem_ptr, rho_multipole_shmem_win, mpierr)
      call check_allocation(mpierr, 'MPI_Win_allocate_shared:rho_multipole_shmem_ptr')

      win_tmp = rho_multipole_shmem_win
      call MPI_Win_shared_query(win_tmp, 0, winsize, disp_unit, rho_multipole_shmem_ptr, mpierr)
      if (winsize /= size_bytes) then
        tmp_int = size_bytes
        print*, "disp_unit:", disp_unit, "winsize:", winsize, ", size_bytes:", size_bytes
        call check_allocation(winsize-tmp_int, 'MPI_Win_shared_query:winsize-size_bytes')
      end if
      call check_allocation(mpierr, 'MPI_Win_shared_query:rho_multipole_shmem_ptr')

      call C_F_pointer(rho_multipole_shmem_ptr, rho_multipole_shmem, rho_multipole_shape)

      if (shm_rank .eq. 0) then
        call MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, rho_multipole_shmem_win, mpierr)
        call MPI_Win_sync(rho_multipole_shmem_win, mpierr)
        rho_multipole_shmem = 0
        call MPI_Win_unlock(0, rho_multipole_shmem_win, mpierr)
      end if
    endif
    ! ----------------------------------------------------------------

    if (packed_matrix_format /= PM_none .and. n_periodic > 0.and.use_DFPT_phonon.or.use_friction) then
     allocate(rho_multipole_supercell((l_pot_max+1)**2, n_max_radial+2, n_centers_in_sc_DFPT),stat=istat)
     call check_allocation(istat, 'rho_multipole_supercell')
     rho_multipole_supercell = 0.0d0
     allocate(delta_v_hartree_part_at_zero_supercell(n_centers_in_sc_DFPT),stat=istat)
     call check_allocation(istat, 'delta_v_hartree_part_at_zero  ')
     delta_v_hartree_part_at_zero_supercell = 0.0d0

     allocate(delta_v_hartree_deriv_l0_at_zero_supercell(3,n_centers_in_sc_DFPT),stat=istat)
     call check_allocation(istat, 'delta_v_hartree_deriv_l0_at_zero_supercell')
     delta_v_hartree_deriv_l0_at_zero_supercell = 0.0d0

    endif

  end subroutine reset_hartree_potential_storage
  !******
  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/set_hartree_potential_atom_lists
  !  NAME
  !    set_hartree_potential_atom_lists
  !  SYNOPSIS

  subroutine set_hartree_potential_atom_lists(n_my_batches_work, batches_work, partition_tab_work)

    !  PURPOSE
    !    Internal routine for finding the atoms involved in Hartree potential
    !    on current task and setting the related lists
    !  USES

    use mpi_tasks, only: myid, n_tasks
    implicit none

    !  ARGUMENTS

    integer, intent(in) :: n_my_batches_work
    type (batch_of_points), intent(in) :: batches_work(n_my_batches_work)
    real*8, intent(in) :: partition_tab_work(*)

    !  INPUTS
    !    n_my_batches_work -- number of batches for current distribution
    !    batches_work -- batches for current distribution
    !    partition_tab_work -- partition_tab_work for current distribution
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


    logical :: need_atom(n_atoms)

    real*8 dist_min_to_center
    real*8 coord_current(3)
    real*8 dist_tab_sq
    real*8 dir_tab(3)
    real*8 max_radius_sq

    integer i_center, current_center, current_atom, i_full_points, i_batch, i_atom, i_index, istat
    integer max_rho_multipole_atoms

    character*150 :: info_str

    ! Flag all atoms which must be included in rho_multipole
    ! If use_distributed_spline_storage is set to .false., we are including all atoms
    ! which may be necessary for postprocessing or output.

    ! Maybe the name "use_distributed_spline_storage" is not appropriate any more and should be changed.
    ! This flag is also used for the Kerker splines so it might be better to have 2 different flags here.

    if(use_distributed_spline_storage) then

      ! Check for every atom if one of the centers related with this atom is more near than
      ! the outer radius in r_radial - in this case the atom is potentially needed
      ! - for settung up rho_multipole
      ! - for evaluating rho_multipole splines in sum_up_whole_potential

      ! VB: This logic may be flawed, at least for periodic systems.
      !
      !     - The role of r_radial here is not clear to me 
      !     - In periodic systems, the reliance on the list of atomic images
      !       in centers_hartree_potential is questionable, this list is limited
      !       and not checked with respect to r_radial
      !     - I would have thought that using map_to_center_cell on each grid point
      !       would be needed, as the coordinates of each point reduced to the zeroth 
      !       unit cellmatter for each batch. Must check in sum_up_whole_potential that 
      !       this is really true, though.

      need_atom(:) = .false. ! By default none is needed


      do i_center = 1, n_centers_hartree_potential, 1

        current_center = centers_hartree_potential(i_center)
        current_atom   = center_to_atom(current_center)
        max_radius_sq  = r_radial(n_radial(species(current_atom)),species(current_atom))**2

        ! For every center in centers_hartree_potential:
        ! Get minimum distance to any integration point owned by this task

        dist_min_to_center = 1.d300
        i_full_points = 0

        do i_batch = 1, n_my_batches_work

          ! loop over one batch
          do i_index = 1, batches_work(i_batch)%size, 1

            i_full_points = i_full_points + 1

            if (partition_tab_work(i_full_points).gt.0.d0) then

              ! get current integration point coordinate
              coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)

              call tab_single_atom_centered_coords_p0 &
                       ( current_center, &
                       coord_current,  &
                       dist_tab_sq,  &
                       dir_tab )

              if(dist_tab_sq<dist_min_to_center) dist_min_to_center = dist_tab_sq
            endif
          enddo
        enddo

        ! Decide if we need current_atom.
        ! To be safe against rounding errors, we use max_radius_sq*1.001
        if(dist_min_to_center < max_radius_sq*1.001) need_atom(current_atom) = .true.

      enddo

      ! There are some load distribution loops assigning atoms in a round robin way across tasks.
      ! Include these atoms also - this should be avoidable by a more intelligent way of
      ! assigning these atoms to tasks, but the few extra atoms shouldn't cost too much

      do i_atom = 1, n_atoms
        if (mod(i_atom-1,n_tasks) == myid) need_atom(i_atom) = .true.
      enddo

    else

      ! if use_distributed_spline_storage is not set, just include all atoms in rho_multipole
      need_atom(:) = .true.

    endif

    ! If n_rho_multipole_atoms > 0 (i.e. when called from reset_hartree_potential_storage),
    ! we must also include the old atoms since otherways update_hartree_potential_p1 stops working.
    ! Maybe this can be optimized later by changing update_hartree_potential_p1, but the additional
    ! space needed for that shouldn't be too serious.

    do i_atom = 1, n_rho_multipole_atoms
      need_atom(i_rho_multipole_atoms(i_atom)) = .true.
    enddo

    ! Get the atoms for which we need rho_multipole

    n_rho_multipole_atoms = 0
    do i_atom = 1, n_atoms
      if(need_atom(i_atom)) then
        n_rho_multipole_atoms = n_rho_multipole_atoms+1
        i_rho_multipole_atoms(n_rho_multipole_atoms) = i_atom
      endif
    enddo

    call sync_find_max(n_rho_multipole_atoms, max_rho_multipole_atoms)
    write(info_str,'(2x,a,i12)') 'Max. number of atoms included in rho_multipole: ',max_rho_multipole_atoms
    call localorb_info(info_str)

    rho_multipole_index(:) = 0

    do i_atom = 1, n_rho_multipole_atoms
      rho_multipole_index(i_rho_multipole_atoms(i_atom)) = i_atom
    enddo

  end subroutine set_hartree_potential_atom_lists
  !******
  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/cleanup_hartree_potential_storage
  !  NAME
  !    cleanup_hartree_potential_storage
  !  SYNOPSIS

  subroutine cleanup_hartree_potential_storage

    !  PURPOSE
    !    Deallocate module arrays.
    !  USES

    implicit none

    !  ARGUMENTS

    ! none

    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


    if (allocated(i_rho_multipole_atoms))    deallocate(i_rho_multipole_atoms)
    if (allocated(rho_multipole_index))      deallocate(rho_multipole_index)
    if (allocated(rho_multipole))            deallocate(rho_multipole)
    if (allocated(rho_multipole_supercell))  deallocate(rho_multipole_supercell)
    if (allocated(delta_v_hartree_part_at_zero_supercell)) deallocate(delta_v_hartree_part_at_zero_supercell)
    if (allocated(delta_v_hartree_deriv_l0_at_zero_supercell)) & 
                                             deallocate(delta_v_hartree_deriv_l0_at_zero_supercell)

    if (allocated(original_multipole_moments)) deallocate(original_multipole_moments)

    if (allocated(compensation_radius)) deallocate(compensation_radius)
    if (allocated(compensation_norm))   deallocate(compensation_norm)

    n_rho_multipole_atoms = 0 ! important for set_hartree_potential_atom_lists

  end subroutine cleanup_hartree_potential_storage
  !******
  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/get_rho_multipole_spl
  !  NAME
  !    get_rho_multipole_spl
  !  SYNOPSIS

  subroutine get_rho_multipole_spl(rho_multipole_spl, spl_atom)

    !  PURPOSE
    !    Delivers the spline coefficients for rho_multipole
    !  USES

    use mpi_tasks, only: myid, aims_stop
    implicit none

    !  ARGUMENTS

    real*8, intent(out) :: rho_multipole_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2)
    integer, intent(in) :: spl_atom

    !  INPUTS
    !    o spl_atom -- the atom for which the spline coefficients are needed
    !  OUTPUTS
    !    o rho_multipole_spl -- the spline coefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


    integer i_radial, l_h_dim, n_rad, j, i_atom_index

    real*8 :: i_r_outer
    real*8 :: delta, delta_2, delta_3

    ! Get index of spl_atom in rho_multipole and check if rho_multipole for this atom is stored

    i_atom_index = rho_multipole_index(spl_atom)
    if(i_atom_index<=0) then
      print '(2(a,i5))','ID ',myid,' INTERNAL ERROR get_rho_multipole_spl - need atom! Atom: ',spl_atom
      call aims_stop
    endif

    l_h_dim = (l_hartree(species(spl_atom))+1)**2
    n_rad   = n_radial(species(spl_atom))

    ! First, simply produce spline coefficients from the charge density as
    ! given on the radial integration grid.

    if (.not. use_rho_multipole_shmem) then
      rho_multipole_spl(1:l_h_dim,1,1:n_rad+2) = rho_multipole(1:l_h_dim,1:n_rad+2,i_atom_index)
    else
      rho_multipole_spl(1:l_h_dim,1,1:n_rad+2) = rho_multipole_shmem(1:l_h_dim,1:n_rad+2,i_atom_index)
    endif

    ! Spline interpolation
    call cubic_spline_v2(rho_multipole_spl, (l_pot_max+1)**2, n_max_spline, n_max_radial+2, &
                         n_rad+2, l_h_dim)


    ! Splines are now tabulated up to r = infinity in principle.
    ! NOW, "doctor" all splines in the far field:
    ! We know that no density must occur outside the radius of the free atom,
    ! because the partition table is zero there.
    ! Therefore, extrapolate from last radial shell inside free atom radius to 
    ! become zero at multipole_radius_free

    ! find outermost radial grid point that is possibly non-zero
    i_radial = n_rad
    do while ( ( r_radial(i_radial,species(spl_atom)) .ge. &
         multipole_radius_free(species(spl_atom)) ) &
         .and.(i_radial.gt.1) )

       rho_multipole_spl(1:l_h_dim,:,i_radial+1) = 0.d0

       i_radial = i_radial - 1
    enddo

    ! Outermost atom radius in units of the radial integration grid
    i_r_outer = invert_radial_grid &
         ( multipole_radius_free(species(spl_atom)), &
         n_rad, &
         scale_radial(species(spl_atom)) )

    delta = dble(i_r_outer - i_radial)
    delta_2 = delta*delta
    delta_3 = delta_2*delta

    ! This is an ugly hack because the element i_radial+1 in rho_multipole_spl
    ! now corresponds to radial shell r_radial(i_radial). What a mess.
    i_radial = i_radial + 1

    ! doctor the spline coefficients at the outermost finite value
    ! i_radial to go smoothly to zero at multipole_radius_free

    do j = 1, l_h_dim

      rho_multipole_spl( j, 3, i_radial) = &
               - 3.d0 / delta_2 * rho_multipole_spl( j, 1, i_radial) &
               - 2.d0 / delta   * rho_multipole_spl( j, 2, i_radial)

      rho_multipole_spl( j, 4, i_radial) = &
                 2.d0 / delta_3 * rho_multipole_spl( j, 1, i_radial) &
               + 1.d0 / delta_2 * rho_multipole_spl( j, 2, i_radial)

    enddo

  end subroutine get_rho_multipole_spl

  !******
  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/get_rho_multipole_supercell_spl
  !  NAME
  !    get_rho_multipole_supercell_spl
  !  SYNOPSIS

  subroutine get_rho_multipole_supercell_spl(rho_multipole_supercell_spl, spl_center, spl_atom)

    !  PURPOSE
    !    Delivers the spline coefficients for rho_multipole_supercell
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(out) :: rho_multipole_supercell_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2)
    integer, intent(in) :: spl_center
    integer, intent(in) :: spl_atom

    !  INPUTS
    !    o spl_center -- the atom for which the spline coefficients are needed
    !  OUTPUTS
    !    o rho_multipole_spl -- the spline coefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


    integer i_radial, l_h_dim, n_rad, j

    real*8 :: i_r_outer
    real*8 :: delta, delta_2, delta_3


    l_h_dim = (l_hartree(species(spl_atom))+1)**2
    n_rad   = n_radial(species(spl_atom))

    ! First, simply produce spline coefficients from the charge density as
    ! given on the radial integration grid.

    rho_multipole_supercell_spl(1:l_h_dim,1,1:n_rad+2) =  & 
             rho_multipole_supercell(1:l_h_dim,1:n_rad+2,spl_center)

    ! Spline interpolation
    call cubic_spline_v2(rho_multipole_supercell_spl, (l_pot_max+1)**2, n_max_spline, n_max_radial+2, &
                         n_rad+2, l_h_dim)


    ! Splines are now tabulated up to r = infinity in principle.
    ! NOW, "doctor" all splines in the far field:
    ! We know that no density must occur outside the radius of the free atom,
    ! because the partition table is zero there.
    ! Therefore, extrapolate from last radial shell inside free atom radius to 
    ! become zero at multipole_radius_free

    ! find outermost radial grid point that is possibly non-zero
    i_radial = n_rad
    do while ( ( r_radial(i_radial,species(spl_atom)) .ge. &
         multipole_radius_free(species(spl_atom)) ) &
         .and.(i_radial.gt.1) )

       rho_multipole_supercell_spl(1:l_h_dim,:,i_radial+1) = 0.d0

       i_radial = i_radial - 1
    enddo

    ! Outermost atom radius in units of the radial integration grid
    i_r_outer = invert_radial_grid &
         ( multipole_radius_free(species(spl_atom)), &
         n_rad, &
         scale_radial(species(spl_atom)) )

    delta = dble(i_r_outer - i_radial)
    delta_2 = delta*delta
    delta_3 = delta_2*delta

    ! This is an ugly hack because the element i_radial+1 in rho_multipole_spl
    ! now corresponds to radial shell r_radial(i_radial). What a mess.
    i_radial = i_radial + 1

    ! doctor the spline coefficients at the outermost finite value
    ! i_radial to go smoothly to zero at multipole_radius_free

    do j = 1, l_h_dim

      rho_multipole_supercell_spl( j, 3, i_radial) = &
               - 3.d0 / delta_2 * rho_multipole_supercell_spl( j, 1, i_radial) &
               - 2.d0 / delta   * rho_multipole_supercell_spl( j, 2, i_radial)

      rho_multipole_supercell_spl( j, 4, i_radial) = &
                 2.d0 / delta_3 * rho_multipole_supercell_spl( j, 1, i_radial) &
               + 1.d0 / delta_2 * rho_multipole_supercell_spl( j, 2, i_radial)

    enddo

  end subroutine get_rho_multipole_supercell_spl



  !******
  !------------------------------------------------------------------------------
  !****s* hartree_potential_storage/get_multipole_moments_on_original_grid
  !  NAME
  !    get_multipole_moments_on_original_grid
  !  SYNOPSIS

  subroutine get_multipole_moments_on_original_grid(i_atom)

    !  PURPOSE
    !    Evaluates the multipole components of the density on the original
    !    grid that they are computed on - the "radial" grid. In principle,
    !    these components should be exact, since the density is exactly
    !    normalized on that grid.
    !
    !    In practice, they are not exact. It is important to find
    !    out how to compensate them properly.
    !  USES

    use mpi_tasks, only: myid, aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(in) :: i_atom

    !  INPUTS
    !    o i_atom -- the atom for which the moment is evaluated here.
    !  OUTPUTS
    !    o No output. We just set the values of original_multipole_moments .
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2012).
    !  SOURCE

    integer :: i_atom_index
    integer :: i_l, i_m, index_lm
    integer :: i_radial

    i_atom_index = rho_multipole_index(i_atom)
    if(i_atom_index<=0) then
      print '(2(a,i5))','ID ',myid, & 
        ' INTERNAL ERROR in get_multipole_moments_on_original_grid - need atom! Atom: ', &
        i_atom
      call aims_stop
    endif

    index_lm = 0

    if (.not. use_rho_multipole_shmem) then
      do i_l = 0, l_hartree(species(i_atom)), 1
         do i_m = -i_l, i_l, 1
            index_lm = index_lm+1

            ! The radial shells on which the density is normalized start at
            ! i_radial = 2 and end at i_radial = n_radial(species(i_atom))+1
            ! i_radial = 1 is r=0 and i_radial = n_radial(species(i_atom))+2 
            ! is r=infinity
            do i_radial = 2, n_radial(species(i_atom))+1, 1

               original_multipole_moments (index_lm, i_atom) = &
                 original_multipole_moments (index_lm, i_atom) + &
                 rho_multipole(index_lm,i_radial,i_atom_index) * w_radial(i_radial-1,species(i_atom)) * &
                 r_radial(i_radial-1,species(i_atom))**2 * &
                 r_radial(i_radial-1,species(i_atom))**(i_l) / (2.d0*dble(i_l)+1.d0)

                !  rho_multipole_shmem(index_lm,i_radial,i_atom_index) * w_radial(i_radial-1,species(i_atom)) * &
            enddo

         enddo
      enddo
    else
      do i_l = 0, l_hartree(species(i_atom)), 1
         do i_m = -i_l, i_l, 1
            index_lm = index_lm+1

            ! The radial shells on which the density is normalized start at
            ! i_radial = 2 and end at i_radial = n_radial(species(i_atom))+1
            ! i_radial = 1 is r=0 and i_radial = n_radial(species(i_atom))+2 
            ! is r=infinity
            do i_radial = 2, n_radial(species(i_atom))+1, 1

               original_multipole_moments (index_lm, i_atom) = &
                 original_multipole_moments (index_lm, i_atom) + &
                 rho_multipole_shmem(index_lm,i_radial,i_atom_index) * w_radial(i_radial-1,species(i_atom)) * &
                 r_radial(i_radial-1,species(i_atom))**2 * &
                 r_radial(i_radial-1,species(i_atom))**(i_l) / (2.d0*dble(i_l)+1.d0)

            enddo

         enddo
      enddo
    endif ! if (.not. use_rho_multipole_shmem) then

  end subroutine get_multipole_moments_on_original_grid
  !******
  !------------------------------------------------------------------------------
  !****f* hartree_potential_storage/compensating_density
  !  NAME
  !    compensating_density
  !  SYNOPSIS

  real*8 function compensating_density (radius, r_outer, l)

    !  PURPOSE
    !
    !  A compensating radial density of the form
    !
    !  P_l(r) = r^l * [ 2*(r/r_outer)^3 - 3*(r/r_outer)^2 + 1 ]
    !
    !  i.e., a third-order polynomial with zero derivative at 
    !  r = 0 and at r = r_outer, times an angular momentum prefactor
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: radius
    real*8, intent(IN) :: r_outer
    integer, intent(IN) :: l

    !  INPUTS
    !    o radius  - the 1-d coordinate at which the compensating density 
    !               is evaluated
    !    o r_outer - outermost radial coordinate, at which the compensating density 
    !                becomes zero.
    !    o l       - the angular momentum for which the compensating density is evaluated. 
    !  OUTPUTS
    !    o The function itself.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2012).
    !  SOURCE

    real*8 rl, rfrac, rfrac2, rfrac3
    integer :: i_l

    if (radius.ge.r_outer) then
      compensating_density = 0.d0
    else
      rl = 1.d0
      do i_l = 1, l, 1
        rl = rl * radius
      enddo

      rfrac = radius/r_outer
      rfrac2 = rfrac*rfrac
      rfrac3 = rfrac2*rfrac

      compensating_density = rl * ( 2.d0 * rfrac3 - 3.d0 * rfrac2 + 1.d0 )
    end if

    return
  end function compensating_density
  !******

end module hartree_potential_storage
!******
