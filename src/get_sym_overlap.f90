!****s* FHI-aims/get_sym_overlap
!  NAME
!    get_sym_overlap
!  SYNOPSIS

subroutine get_sym_overlap  &
     ()

!  PURPOSE
!
!  Wrapper function for symmetry operation
!
!  USES
  use spglib_symmetry,only: destroy_symmetry, write_symm_info, out_symm_mats,&
                            spg_rotations, destroy_symmats, num_symmetry
  ! use pbc_lists, only: n_centers,coords_center, n_full_points,center_to_atom
  use dimensions, only: l_wave_max, use_hf_kspace, &
                        n_atoms, n_k_points, n_basis,n_k_points_task
                        ! n_max_radial,n_max_angular, n_hamiltonian_matrix_size,n_my_batches
  use constants, only:pi
  use sym_base, only: symmattoeuler,get_TVlmm, map_atoms, reset_map, inv,&
                      construct_C_full,construct_Delta,apply_T, Delta,&
                      get_density_rotation, Delta_matrix_full,rotations_at_task,&
                      n_rotations_at_task, map_atom, map_atoms_all, &
		      destroy_symmetry_arrays
                      ! map_centers_syms, get_inv_sym
                      ! rotations_reciprocal, translations_reciprocal, &
                      ! count_points,map_real_space, map_real,local_partition_tab,&
                      ! prepare_local_partition_tab,partition_tab_sym,&
                      ! hartree_partition_tab_sym,clean_local_partition_tab,&
                      ! partition_tab_sym_batches, hartree_partition_tab_sym_batches,&
                      ! map_grid
  use mpi_tasks,only: check_allocation, aims_stop, myid, n_tasks, mpi_comm_global
  ! use grids,only: batches,n_angular,n_radial,r_angular,r_radial
  ! use physics,only: rho, partition_tab,hartree_partition_tab, KS_eigenvector,KS_eigenvector_complex
  use runtime_choices,only: use_spg_full_Delta !packed_matrix_format, PM_none
  use geometry,only:lattice_vector
  implicit none

!  ARGUMENTS

!  INPUTS

!  OUTPUT
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
  real*8  :: alpha, beta, gamma, det
  integer :: is, p, i_k, i_k_point, i_sym, num_rots
  integer :: info
  real*8  :: temp(3,3), temp_inv(3,3), ang(3)
  complex*16, dimension(:,:,:), allocatable :: C_full
  complex*16, dimension(:,:,:), allocatable :: TVlmm

! real space grid
!   integer :: n_points_in_grid
!   integer :: i_my_batch,i_index,i_point
!   real*8 :: partition_tab_local, hartree_partition_tab_local
!   real*8 :: coord_current(3)
!   integer :: current_atom,current_radial,current_angular, num_inv_symmetry
!   real*8, dimension(:,:,:), allocatable :: rho_transfrom
!   integer :: i_radial,i_atom,i_angular,i_grid

  ! Set all indexing arrays to comply with FHI-aims conventions
  ! e.g. mapping between reduced (distributed) and full k-point se
  
  call destroy_symmetry_arrays()

  call reset_map()
  ! Map symmetry transformed atoms to 1st unitcell and caluclate shift for k-phase 
  if (use_hf_kspace)then
      call map_atoms_all()
      num_rots = num_symmetry
  else
      call map_atoms_all()
      call map_atoms()
      num_rots = n_rotations_at_task
  endif
  ! Transformation matrix from basis of complex Y_lm to real Y_lm
  if (.not.allocated(C_full))then
    allocate ( C_full(0:l_wave_max,-l_wave_max:l_wave_max,&
                     -l_wave_max:l_wave_max) ,stat=info)
      call check_allocation(info, 'C_full                       ')
  endif
  call construct_C_full(l_wave_max, C_full)
  ! Apply FHI-aims Y_lm sign convention
  call apply_T(l_wave_max,C_full)
  ! Wigner D-Matrix - Rotation matrix in basis of complex spherical harmonics
  if (.not.allocated(TVlmm))then
    allocate ( TVlmm(0:l_wave_max,-l_wave_max:l_wave_max,&
                     -l_wave_max:l_wave_max) ,stat=info)
      call check_allocation(info, 'TVlmm                       ')
  endif
  ! Delta Matrix - Rotation matrix in basis of real spherical harmonics
  if (.not.allocated(Delta))then
    allocate ( Delta(0:l_wave_max,-l_wave_max:l_wave_max,&
                     -l_wave_max:l_wave_max,num_rots) ,stat=info)
      call check_allocation(info, 'Delta                       ')
  endif
  do is = 1, num_rots, 1
    ! We have to transform the rotation matrix from fractional coordinates 
    ! (coordinates of the lattice vectors) to cartesian coordinates to 
    ! determine the Euler angles
    if (use_hf_kspace)then
        i_sym = is
    else
        i_sym = rotations_at_task(is)
    endif
    temp = (matmul((lattice_vector),matmul((spg_rotations(1:3,1:3,i_sym)),(inv(lattice_vector)))))
    temp_inv = inv(temp)
    !debug write(use_unit,'(A,I5)') '# Rotation:', rotations_at_task(is)
    !debug write(use_unit,'(3F11.4)') temp_inv(1:3,1)
    !debug write(use_unit,'(3F11.4)') temp_inv(1:3,2)
    !debug write(use_unit,'(3F11.4)') temp_inv(1:3,3)
    !debug write(use_unit,'(A)') '--------------------------'
    !debug write(use_unit,(3F11.4)) spg_shift(1:3,rotations_at_task(is))
    !debug write(use_unit,'(A)') '-----------------------------------'
    
    ! Determine if rotation proper/improper
    det = temp(1,2) * temp(2,3) * temp(3,1) &
        - temp(1,3) * temp(2,2) * temp(3,1) &
        + temp(1,3) * temp(2,1) * temp(3,2) &
        - temp(1,1) * temp(2,3) * temp(3,2) &
        + temp(1,1) * temp(2,2) * temp(3,3) &
        - temp(1,2) * temp(2,1) * temp(3,3) 


    if (det > 0d0) then
         p = 1
    else
         p = - 1
         temp_inv = -temp_inv
    endif


    ! Get the Euler angles
    call symmattoeuler(temp_inv,alpha, beta, gamma)
    !debug write(use_unit,'(A,1F11.4,1F11.4,1F11.4,I4)') 'euler angles',alpha, beta, gamma, p
    ! Ylm transformation matrix for all rotations 
    call get_TVlmm(TVlmm,alpha, beta, gamma,p)
    call construct_Delta(l_wave_max,TVlmm,C_full,Delta(:,:,:,is))
    !debug write(use_unit,'(I5,I5)') myid, n_rotations_at_task
    !debug do is = 1, n_rotations_at_task, 1
    !debug   write(use_unit,'(I5,I5,I5,I5)') myid, is, rotations_at_task(is), rotations_at_task_map(rotations_at_task(is))
    !debug enddo


  enddo

  if(use_spg_full_Delta)then
    ! Full Delta - Rotation matrix in basis of real spherical harmonics
    if (.not.allocated(Delta_matrix_full))then
      allocate ( Delta_matrix_full(n_basis,n_basis,num_rots,n_k_points_task &
                      ) ,stat=info)
        call check_allocation(info, 'Delta_matrix_full           ')
    endif
    ! Full rotation matrix for complex spherical harmonics including phase factor
    ! Matrix is constructed for all rotations and k-points at task
    do is = 1, num_rots, 1
      if (use_hf_kspace)then
	  i_sym = is
      else
	  i_sym = rotations_at_task(is)
      endif
      i_k = 0
      do i_k_point = 1, n_k_points, 1
        if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
            i_k = i_k + 1
            call get_density_rotation(i_sym, i_k_point, &
                map_atom(:,is), Delta(:,:,:,is), Delta_matrix_full(:,:,is,i_k))
        end if !myid == mod(k...
      end do ! i_k_point
    enddo

    if(allocated(Delta)) deallocate(Delta)
    if(allocated(map_atom)) deallocate(map_atom)
    !call destroy_symmats()
  endif

  if(allocated(C_full)) deallocate(C_full)
  if(allocated(TVlmm)) deallocate(TVlmm)

!  All the following is for real space grid reduction and reconstruction

!  Invert rotations (including mirror symmetry)
!   call get_inv_sym(num_symmetry, spg_rotations, spg_shift, 0, num_inv_symmetry)
!  Find irreducible centers of periodic replica
!   call map_centers_syms(n_centers,coords_center,&
!                       center_to_atom,num_inv_symmetry, translations_reciprocal,&
!                       rotations_reciprocal)
!       
!   call count_points(n_points_in_grid)
!   call map_real_space(num_symmetry,spg_rotations,spg_shift,n_points_in_grid)

!  call prepare_local_partition_tab()


 !  if (.not.allocated(partition_tab_sym_batches))then
!     allocate ( partition_tab_sym_batches(n_full_points) ,stat=info)
!       call check_allocation(info, 'partition_tab_sym_batches                ')
!   endif
!   if (.not.allocated(hartree_partition_tab_sym_batches))then
!     allocate ( hartree_partition_tab_sym_batches(n_full_points) ,stat=info)
!       call check_allocation(info, 'hartree_partition_tab_sym_batches        ')
 !  endif
! 
 ! if (.not.allocated(rho_transfrom))then
!     allocate ( rho_transfrom(n_atoms,n_max_radial,n_max_angular) ,stat=info)
!       call check_allocation(info, 'rho_transfrom               ')
!   endif

!   rho_transfrom = 0d0

 !  i_point = 0

!   do i_my_batch = 1, n_my_batches, 1

 !     do i_index = 1, batches(i_my_batch)%size, 1

 !        i_point = i_point + 1

 !        current_atom    = batches(i_my_batch) % points(i_index) % index_atom
  !       current_radial  = batches(i_my_batch) % points(i_index) % index_radial
  !       current_angular = batches(i_my_batch) % points(i_index) % index_angular

   !      rho_transfrom(current_atom,current_radial,current_angular) = rho(1,i_point)

!      enddo
!   enddo


!   i_point = 0
!   do i_my_batch = 1, n_my_batches, 1

!      do i_index = 1, batches(i_my_batch)%size, 1

!         i_point = i_point + 1

!         coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
 !        current_atom    = batches(i_my_batch) % points(i_index) % index_atom
 !        current_radial  = batches(i_my_batch) % points(i_index) % index_radial
 !        current_angular = batches(i_my_batch) % points(i_index) % index_angular

 !        call local_partition_tab( partition_tab_local, hartree_partition_tab_local, &
 !                                coord_current, current_atom, &
 !                                current_radial, current_angular )
 !        partition_tab_sym_batches(i_point) = partition_tab_sym(current_atom,current_radial,current_angular)
 !        hartree_partition_tab_sym_batches(i_point) = hartree_partition_tab_sym(current_atom,current_radial,current_angular)
        !write(use_unit,*) i_my_batch, i_index, i_point, rho(1,i_point), rho_transfrom(current_atom,current_radial,current_angular)
        !write(use_unit,*) myid,partition_tab_sym(current_atom,current_radial,current_angular), partition_tab(i_point)
        !write(use_unit,*) myid,hartree_partition_tab_sym(current_atom,current_radial,current_angular),hartree_partition_tab(i_point)
 !     enddo
 !  enddo

 !  i_grid = 0
 !  do i_atom = 1, n_atoms, 1
 !    do i_radial = 1, n_radial(species(i_atom)), 1
  !      do i_angular = 1, n_angular( i_radial,species(i_atom) )
 !          i_grid = i_grid + 1
!           write(use_unit,*) i_atom,i_radial,i_angular,map_grid(1:3,i_grid)
!           write(use_unit,*) partition_tab_sym(i_atom,i_radial,i_angular), rho_transfrom(i_atom,i_radial,i_angular), rho_transfrom(map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid))

 !       end do
 !     end do
 !  end do
 !  call clean_local_partition_tab()

  call destroy_symmetry()
end subroutine get_sym_overlap
