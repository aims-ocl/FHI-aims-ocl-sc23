!****s* FHI-aims/get_n_compute_maxes_p1
!  NAME
!    get_n_compute_maxes_p1
!  SYNOPSIS 
subroutine get_n_compute_maxes_p1( partition_tab )
!  PURPOSE
!    Computes the threadwise maximal and average number of non-zero basis functions.
!    Selects additinally the method for density update if 
!    autoselect_density_method is set.
!  USES
  use dimensions,            only : n_full_points, n_centers_integrals, n_centers, n_max_batch_size, &
                                    n_basis_fns, n_avg_compute_dens, n_centers_basis_I, n_centers_basis_T, &
                                    n_int_points_total, n_max_compute_ang, n_max_compute_atoms, &
                                    n_max_compute_dens, n_max_compute_fns_dens, n_max_compute_fns_ham, &
                                    n_max_compute_ham, n_my_batches, n_periodic, use_gga
  use grids,                 only : batches
  use runtime_choices,       only : autoselect_density_method, got_n_compute_maxes, keep_restart_info, &
                                    lda_update_switch, measure_forces, packed_matrix_format, PM_index, &
                                    prune_basis_once, restart_file_exists, restart_read, restart_write, &
                                    split_updates, use_density_matrix, use_load_balancing, &
                                    use_metis_batch_distribution, use_scalapack, out_batch_statistics,&
                                    flag_rel, REL_x2c, REL_4c_dks, &
                                    use_rlsy_symmetry
  use localorb_io,           only : output_priority, use_unit, OL_norm, OL_low, localorb_info, &
                                    localorb_allinfo
  use mpi_tasks,             only : myid
  use pbc_lists,             only : centers_basis_integrals, inv_centers_basis_integrals
  use synchronize_mpi_basic, only : sync_find_max, bcast_logical
  use synchronize_mpi,       only : sync_n_avg_compute
  use physics,               only : n_electrons
  use timing,                only : relaxation_step_number, time_forces_orbital_update, time_density
  use rlsy_interface,        only : rlsy_h
  use rel_x2c_mod

  implicit none
!  ARGUMENTS
  real*8, dimension(n_full_points), intent(in) :: partition_tab
! the following parameters are needed for the automatic selection of charge density/force update
  real*8, parameter :: gga_prefac1 = 1.1
  real*8, parameter :: gga_prefac2 = 0.82
  real*8, parameter :: gga_prefac3 = 0.62
  real*8, parameter :: lda_prefac1 = 1.55
  real*8, parameter :: lda_prefac2 = 1.1


!  INPUTS
!    o partition_tab     -- the partition tab
!  OUTPUT
!    n_max_compute_ham/dens, n_max_compute_fns_ham/dens, n_max_compute_atoms and
!    n_avg_compute_dens are set on exit
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
! SOURCE


 ! real*8, dimension( n_full_points_hamiltonian_integrals) :: hamiltonian_partition_tab
  
  ! locals
  integer :: i_center, i_center_L
  integer :: i_my_batch, i_index
  integer :: i_point

  integer :: n_compute_c, n_compute_a, n_compute_small
  integer :: i_full_points_A, i_full_points_C, i_full_points_2C
  integer :: i_full_points, i_full_points_2, i_full_points_3

  real*8 :: coord_current(3)
  real*8 :: dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 :: dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 :: dir_tab(3,n_centers_integrals, n_max_batch_size)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

  integer :: n_centers_max, n_points
  integer :: n_max

  integer,dimension(:),allocatable :: i_basis, i_basis_small

  character*200 :: info_str
  integer :: mpierr

  logical, dimension(:), allocatable :: my_functions

  character*40 :: filename
  character*8  :: myid_str
  integer      :: n_min_compute_fns_ham_batch, n_max_compute_fns_ham_batch
  integer      :: n_min_compute_atoms_batch, n_max_compute_atoms_batch
  real*8       :: min_partition_tab_batch, max_partition_tab_batch

  integer      :: i_irr_batch, i_irr_index

  if (got_n_compute_maxes) return

  write (info_str,'(2X,A)') "Obtaining max. number of non-zero basis functions in each batch (get_n_compute_maxes)."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)


  n_centers_max = MAX(n_centers_basis_I, n_centers_basis_T)

  allocate(i_basis(n_centers_max))
  if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) allocate(i_basis_small(n_centers_basis_I_small))

  n_max_compute_ham = 0
  n_max_compute_fns_ham = 0
  n_max_compute_atoms = 0

  n_max_compute_dens = 0
  n_max_compute_fns_dens = 0

  n_avg_compute_dens = 0

  n_max_compute_ang = 0

  i_full_points_C = 0
  i_full_points_2C = 0
  i_full_points_A = 0

  i_basis_fns_inv = 0

  if (use_metis_batch_distribution) then
     allocate(my_functions(n_centers_basis_I))
     my_functions = .false.
  end if

  if (out_batch_statistics) then
    write (info_str,'(2X,A)') "Also outputting batch statistics to disk (one file for each MPI task)."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    write(myid_str,'(I8)')   myid
    write(filename,'(A,A)') "batch_statistics_task_", adjustl(myid_str)
    open(66, file=trim(filename) // ".dat")
 
    write(66, '(A)') &
         "batch_num   batch_size   min(radial_fns)   max(radial_fns)   &
         &max(basis_fns)   min(atoms)   max(atoms)   &
         &min(partition)   max(partition)" 
  end if

  do i_my_batch = 1, n_my_batches, 1

        n_compute_c = 0; n_compute_small = 0
        n_compute_a = 0
        i_basis = 0; if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) i_basis_small=0
     
        i_point = 0

        ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points_2C = i_full_points_2C + 1
           
           if (partition_tab(i_full_points_2C).gt.0.d0) then

              i_point = i_point+1

              ! get current integration point coordinate
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
              
              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if
              
              ! compute atom-centered coordinates of current integration point,
              ! as viewed from all atoms
              call tab_atom_centered_coords_p0 ( coord_current,  &
                   dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
                   n_centers_integrals, centers_basis_integrals )
              
              ! determine which basis functions are relevant at current integration point,
              ! and tabulate their indices
              
              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
              call prune_basis_p0X ( dist_tab_sq(1,i_point), n_compute_a, n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, i_point  )
              if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                call prune_basis_p2_rel ( dist_tab_sq(1,i_point), n_compute_small, i_basis_small,&
                   n_centers_basis_I_small, n_centers_integrals, inv_centers_basis_integrals )
              endif
!NEC_CB Using version which saves massively indirect addressed arrays
              
           end if
        enddo  ! end loop over one part of the angular division              

        n_points = i_point
        
        if (prune_basis_once) then
           batches(i_my_batch)%batch_n_compute = n_compute_a
           allocate(batches(i_my_batch)%batch_i_basis(n_compute_a))
           batches(i_my_batch)%batch_i_basis = i_basis(1:n_compute_a)
        end if

        if (use_metis_batch_distribution) then
           my_functions(i_basis(1:n_compute_a)) = .true.
        end if

        n_max_compute_ham = MAX(n_compute_c, n_compute_small, n_max_compute_ham)

        n_max_compute_ang = MAX(n_max_compute_ang, n_points)
        
        n_max_compute_dens = MAX(n_compute_c, n_max_compute_dens)
        n_avg_compute_dens = n_avg_compute_dens + dble(n_compute_c*n_points)

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_a.gt.0) then

           i_point = 0

           n_min_compute_fns_ham_batch = 0
           n_min_compute_atoms_batch   = 0
           n_max_compute_fns_ham_batch = 0
           n_max_compute_atoms_batch   = 0
           min_partition_tab_batch     = 0.0d0
           max_partition_tab_batch     = 0.0d0

           ! loop over one division of the angular grid
           do i_index = 1, batches(i_my_batch)%size, 1

              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points_C = i_full_points_C + 1
              
              if (partition_tab(i_full_points_C).gt.0.d0) then
                 
                 i_point = i_point+1

                 n_compute_atoms = 0
                 n_compute_fns = 0
!                 i_basis_fns_inv = 0

                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                 ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                 ! without any copying and without doing any unnecessary operations. 
                 ! The price is that the interface is no longer explicit in terms of physical 
                 ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                 call prune_radial_basis_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dist_tab(1,i_point), &
                      dir_tab(1,1,i_point), &
                      n_compute_atoms, atom_index, atom_index_inv, &
                      n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                      i_atom_fns, spline_array_start, spline_array_end, &
                      n_centers_integrals, centers_basis_integrals)
                 
                 ! Max values across all batches on this MPI task
                 n_max_compute_fns_ham = MAX(n_compute_fns, n_max_compute_fns_ham)
                 n_max_compute_fns_dens = MAX(n_compute_fns, n_max_compute_fns_dens)

                 n_max_compute_atoms = MAX(n_compute_atoms, n_max_compute_atoms)

                 ! Min/max values for this batch only (used when outputting batch statistics to disk)
                 n_max_compute_fns_ham_batch   = MAX(n_compute_fns, n_max_compute_fns_ham_batch)
                 n_max_compute_atoms_batch     = MAX(n_compute_atoms, n_max_compute_atoms_batch)
                 max_partition_tab_batch       = MAX(partition_tab(i_full_points_C), max_partition_tab_batch)

                 if (i_index .eq. 1) then
                   n_min_compute_fns_ham_batch = n_compute_fns
                   n_min_compute_atoms_batch   = n_compute_atoms
                   min_partition_tab_batch     = partition_tab(i_full_points_C)
                 else
                   n_min_compute_fns_ham_batch = MIN(n_compute_fns, n_min_compute_fns_ham_batch)
                   n_min_compute_atoms_batch   = MIN(n_compute_atoms, n_min_compute_atoms_batch)
                   min_partition_tab_batch     = MIN(partition_tab(i_full_points_C), min_partition_tab_batch)
                 end if
              end if

           end do ! end loop over a batch

           if (out_batch_statistics) &
                write(66, '(I9, 3X, I10, 3X, I15, 3X, I15, 3X, I14, 3X, I10, 3X, I10, 3X, E14.2, 3X, E14.2)') &
                     i_my_batch, batches(i_my_batch)%size, &
                     n_min_compute_fns_ham_batch, n_max_compute_fns_ham_batch, n_compute_c, &
                     n_min_compute_atoms_batch, n_max_compute_atoms_batch, &
                     min_partition_tab_batch, max_partition_tab_batch
        else
          ! must increment grid counter in this case too, else we get an inconsistent
          ! partition_tab and all sorts of trouble

           i_full_points_C = i_full_points_C + batches(i_my_batch)%size

           if (out_batch_statistics) &
                write(66, '(I9, 3X, I10, 3X, I15, 3X, I15, 3X, I14, 3X, I10, 3X, I10, 3X, E10.2, 3X, E10.2)') &
                     i_my_batch, batches(i_my_batch)%size, &
                     0, 0, n_compute_c, &
                     0, 0, &
                     0.0d0, 0.0d0
        end if ! end if (n_compute.gt.0)
        
     ! end if ! end distribution of tasks over threads
     
  end do ! end loop over batches

  if (out_batch_statistics) then
    close(66)
  end if

  ! YY symmetry
  if (use_rlsy_symmetry) then

  !n_max_compute_ham = 0
  !n_max_compute_fns_ham = 0
  !n_max_compute_atoms = 0

  !n_max_compute_dens = 0
  !n_max_compute_fns_dens = 0

  n_avg_compute_dens = 0

  !n_max_compute_ang = 0

  i_full_points_C = 0
  i_full_points_2C = 0
  i_full_points_A = 0

  i_basis_fns_inv = 0

  if (out_batch_statistics) then
    write (info_str,'(2X,A)') "Also outputting batch statistics to disk (one file for each MPI task)."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    write(myid_str,'(I8)')   myid
    write(filename,'(A,A)') "batch_statistics_symm_task_", adjustl(myid_str)
    open(66, file=trim(filename) // ".dat")
 
    write(66, '(A)') &
         "batch_num   batch_size   min(radial_fns)   max(radial_fns)   &
         &max(basis_fns)   min(atoms)   max(atoms)   &
         &min(partition)   max(partition)" 
  end if

  do i_irr_batch = 1, rlsy_h%grid%n_irr_batch, 1

        n_compute_c = 0; n_compute_small = 0
        n_compute_a = 0
        i_basis = 0; if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) i_basis_small=0
     
        i_point = 0

        ! loop over one batch
        do i_irr_index = 1, rlsy_h%grid%irr_batch(i_irr_batch)%n_point, 1

           i_full_points_2C = i_full_points_2C + 1
           
           !if (partition_tab(i_full_points_2C).gt.0.d0) then

           i_point = i_point+1

           ! get current integration point coordinate
           !coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
           coord_current(:) = rlsy_h%grid%irr_batch(i_irr_batch)%folded_coordinate(:,i_point) 
           
           if(n_periodic > 0)then
              call map_to_center_cell(coord_current(1:3) )
           end if
           
           ! compute atom-centered coordinates of current integration point,
           ! as viewed from all atoms
           call tab_atom_centered_coords_p0 ( coord_current,  &
                dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
                n_centers_integrals, centers_basis_integrals )
           
           ! determine which basis functions are relevant at current integration point,
           ! and tabulate their indices
           
           ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
           call prune_basis_p0X ( dist_tab_sq(1,i_point), n_compute_a, n_compute_c, i_basis,  &
                n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, i_point  )
           !call prune_basis_p0 ( dist_tab_sq(1,i_point), n_compute_a, n_compute_c, i_basis,  &
           !     n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals )
           if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
             call prune_basis_p2_rel ( dist_tab_sq(1,i_point), n_compute_small, i_basis_small,&
                n_centers_basis_I_small, n_centers_integrals, inv_centers_basis_integrals )
           endif
!NEC_CB Usiversion which saves massively indirect addressed arrays
           
           !end if
        enddo  ! end loop over one part of the angular division              

        n_points = i_point
        
        !if (prune_basis_once) then
        !   batches(i_my_batch)%batch_n_compute = n_compute_a
        !   allocate(batches(i_my_batch)%batch_i_basis(n_compute_a))
        !   batches(i_my_batch)%batch_i_basis = i_basis(1:n_compute_a)
        !end if

        if (use_metis_batch_distribution) then
           my_functions(i_basis(1:n_compute_a)) = .true.
        end if

!        write(*,*) "YY get n_max_compute_ham n_compute_a", n_compute_a

        n_max_compute_ham = MAX(n_compute_a, n_compute_c, n_compute_small, n_max_compute_ham) 

        n_max_compute_ang = MAX(n_max_compute_ang, n_points) 
        
        n_max_compute_dens = MAX(n_compute_c, n_max_compute_dens) 
        n_avg_compute_dens = n_avg_compute_dens + dble(n_compute_c*n_points)

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_a.gt.0) then

           i_point = 0

           n_min_compute_fns_ham_batch = 0
           n_min_compute_atoms_batch   = 0
           n_max_compute_fns_ham_batch = 0
           n_max_compute_atoms_batch   = 0
           min_partition_tab_batch     = 0.0d0
           max_partition_tab_batch     = 0.0d0

           ! loop over one division of the angular grid
           do i_irr_index = 1, rlsy_h%grid%irr_batch(i_irr_batch)%n_point, 1

              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points_C = i_full_points_C + 1
              
              !if (partition_tab(i_full_points_C).gt.0.d0) then
              
              i_point = i_point+1

              n_compute_atoms = 0
              n_compute_fns = 0
!              i_basis_fns_inv = 0

              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
              ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
              ! without any copying and without doing any unnecessary operations. 
              ! The price is that the interface is no longer explicit in terms of physical 
              ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

              call prune_radial_basis_p0 &
                   ( dist_tab_sq(1,i_point), &
                   dist_tab(1,i_point), &
                   dir_tab(1,1,i_point), &
                   n_compute_atoms, atom_index, atom_index_inv, &
                   n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                   i_atom_fns, spline_array_start, spline_array_end, &
                   n_centers_integrals, centers_basis_integrals)
              
              ! Max values across all batches on this MPI task
              n_max_compute_fns_ham = MAX(n_compute_fns, n_max_compute_fns_ham)
              n_max_compute_fns_dens = MAX(n_compute_fns, n_max_compute_fns_dens)

              n_max_compute_atoms = MAX(n_compute_atoms, n_max_compute_atoms)

              ! Min/max values for this batch only (used when outputting batch statistics to disk)
              n_max_compute_fns_ham_batch   = MAX(n_compute_fns, n_max_compute_fns_ham_batch)
              n_max_compute_atoms_batch     = MAX(n_compute_atoms, n_max_compute_atoms_batch) 
              max_partition_tab_batch       = MAX(partition_tab(i_full_points_C), max_partition_tab_batch) 

              if (i_irr_index .eq. 1) then
                n_min_compute_fns_ham_batch = n_compute_fns
                n_min_compute_atoms_batch   = n_compute_atoms
                min_partition_tab_batch     = partition_tab(i_full_points_C)
              else
                n_min_compute_fns_ham_batch = MIN(n_compute_fns, n_min_compute_fns_ham_batch)
                n_min_compute_atoms_batch   = MIN(n_compute_atoms, n_min_compute_atoms_batch)
                min_partition_tab_batch     = MIN(partition_tab(i_full_points_C), min_partition_tab_batch)
              end if
              !end if

           end do ! end loop over a batch

           if (out_batch_statistics) &
                write(66, '(I9, 3X, I10, 3X, I15, 3X, I15, 3X, I14, 3X, I10, 3X, I10, 3X, E14.2, 3X, E14.2)') &
                     i_irr_batch, rlsy_h%grid%irr_batch(i_irr_batch)%n_point, &
                     n_min_compute_fns_ham_batch, n_max_compute_fns_ham_batch, n_compute_c, &
                     n_min_compute_atoms_batch, n_max_compute_atoms_batch, &
                     min_partition_tab_batch, max_partition_tab_batch
        else
          ! must increment grid counter in this case too, else we get an inconsistent
          ! partition_tab and all sorts of trouble

           i_full_points_C = i_full_points_C + batches(i_my_batch)%size

           if (out_batch_statistics) &
                write(66, '(I9, 3X, I10, 3X, I15, 3X, I15, 3X, I14, 3X, I10, 3X, I10, 3X, E10.2, 3X, E10.2)') &
                     i_irr_batch, rlsy_h%grid%irr_batch(i_irr_batch)%n_point, &
                     0, 0, n_compute_c, &
                     0, 0, &
                     0.0d0, 0.0d0
        end if ! end if (n_compute.gt.0)
        
     ! end if ! end distribution of tasks over threads
     
  end do ! end loop over batches

  if (out_batch_statistics) then
    close(66)
  end if

  end if ! use_rlsy_symmetry
  ! YY symmetry

  if(use_load_balancing) then
     ! When load balancing is used, we need the global max numbers since batches are permuted
     n_max = n_max_compute_ham;      call sync_find_max(n_max, n_max_compute_ham)
     n_max = n_max_compute_fns_ham;  call sync_find_max(n_max, n_max_compute_fns_ham)
     n_max = n_max_compute_atoms;    call sync_find_max(n_max, n_max_compute_atoms)
     n_max = n_max_compute_dens;     call sync_find_max(n_max, n_max_compute_dens)
     n_max = n_max_compute_fns_dens; call sync_find_max(n_max, n_max_compute_fns_dens)
     n_max = n_max_compute_ang;      call sync_find_max(n_max, n_max_compute_ang)
     write(info_str,'(2X,A,I8)') "| Maximal number of non-zero basis functions: ",n_max_compute_ham
     call localorb_info(info_str)
  endif

  if (use_metis_batch_distribution .and. output_priority.le.1) then 
     write(info_str, '(2X,A,I8,A,A,I8,A,I3)') &
     & "| Using ", COUNT(my_functions), " distinct basis functions ", &
     & "out of ", n_centers_basis_I, " in task ", myid
     call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  end if

  call sync_n_avg_compute( n_avg_compute_dens )
  n_avg_compute_dens = n_avg_compute_dens / dble(n_int_points_total)
  if(.not.use_load_balancing .and. output_priority .le. OL_low) then
     write(info_str,'(2X,A,I8,A,I5)') &
     & "| Maximal number of non-zero basis functions: ", &
     & n_max_compute_ham, " in task ", myid
     call localorb_allinfo(info_str, use_unit, '(A)', OL_low)
  end if

! autoselection of density update method
  if (autoselect_density_method .and. (relaxation_step_number == 0)) then
     write(info_str,'(2X,A)') "Selecting the method for density update."
     call localorb_info(info_str,use_unit,'(A)')
     if (use_gga) then
        if (n_electrons > 2 * gga_prefac1 * n_avg_compute_dens) then
           use_density_matrix = .true.
        elseif (n_electrons > 2 * gga_prefac2 * n_avg_compute_dens) then
           use_density_matrix = .false.
           split_updates= .true.
        elseif (n_electrons > 2 * gga_prefac3 * n_avg_compute_dens) then
           use_density_matrix = .false.
           measure_forces = .true. !try both force updates and select the faster one
           split_updates = .false.
        else
           use_density_matrix = .false.
        end if
     else
        if (n_electrons > 2 * lda_prefac1 * n_avg_compute_dens) then
           use_density_matrix = .true.
        elseif (n_electrons > 2 * lda_prefac2 * n_avg_compute_dens .and. use_scalapack) then
        lda_update_switch = .true. ! density update via density matrix; if forces needed: updates via orbitals
           use_density_matrix = .true.
        else
           use_density_matrix = .false.
        end if
     end if

     if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) use_density_matrix = .true.

     if (use_density_matrix .and. .not. lda_update_switch) then
        write(info_str,'(2X,A)') "Density matrix based charge density update selected."
        call localorb_info(info_str,use_unit,'(A)')
     elseif (split_updates) then
        write(info_str,'(2X,A,A)') &
            "split updates selected: force will be updated via density matrix ", &
            "and electron density via Kohn Sham orbitals."
        call localorb_info(info_str,use_unit,'(A)')
     elseif (measure_forces) then
        if (myid.eq.0) then 
           write(use_unit,'(2X,A)')  "Electron density will be updated via Kohn Sham orbitals; " 
           write(use_unit,'(2X,A)')  "first relaxation step: forces will be updated via orbitals;"
           write(use_unit,'(2X,A)')  "second relaxation step: forces will be updated via density matrix;"
           write(use_unit,'(2X,A)')  "finally, the fastest method will be selected"
        end if
     elseif (lda_update_switch) then
        write(info_str,'(2X,A,A)') &
            "If forces are needed both density and forces will be updated via orbitals, ", &
            "otherwise density will be updated via density matrix"
        call localorb_info(info_str,use_unit,'(A)')
     else
        write(info_str,'(2X,A)') "Loop over occupied states selected for charge density update."
        call localorb_info(info_str,use_unit,'(A)')
     end if
     ! the following  if-statement is copied from read_control.f90 as the final decision for
     ! use_density_matrix is made here
     if ( (use_density_matrix) .and. (use_scalapack) .and. .not. lda_update_switch ) then
        packed_matrix_format = PM_index
        if (myid.eq.0) then
           write(use_unit,'(2X,A)') 'Scalapack and density-matrix based density update requested:'
           write(use_unit,'(2X,A)') 'Switching to packed matrix format "index" by default.'
        end if

        call bcast_logical( keep_restart_info, 0)
        call bcast_logical( restart_read, 0)
        call bcast_logical( restart_write, 0)
        call bcast_logical(restart_file_exists,0)
     end if
     if ((measure_forces .or. lda_update_switch) .and. use_scalapack) then
        packed_matrix_format = PM_index
        if (myid.eq.0) then
           write(use_unit,'(2X,A)') 'Scalapack and automatic selection of force update requested:'
           write(use_unit,'(2X,A)') 'Switching to packed matrix format "index" by default.'
        end if
     end if
     if (split_updates .and. use_scalapack) then
        packed_matrix_format = PM_index
        if (myid.eq.0) then 
           write(use_unit,'(2X,A)') 'Scalapack and split updates requested:'
           write(use_unit,'(2X,A)') 'Switching to packed matrix format "index" by default.'
        end if
     end if


     if (.not. measure_forces ) then ! the update method should only be chosen once
        autoselect_density_method =.false.
     end if
  elseif (autoselect_density_method .and. (relaxation_step_number == 1)) then
     time_forces_orbital_update=time_density
     split_updates=.true.
  elseif (autoselect_density_method .and. (relaxation_step_number == 2)) then
     if (myid.eq.0) then
        write(use_unit,'(2X,A,F13.2,A)') 'Time for orbital based force update:',time_forces_orbital_update,'s'
        write(use_unit,'(2X,A,F13.2,A)') 'Time for density based force update:',time_density,'s'
     end if
     if (time_density - time_forces_orbital_update > 0 ) then
        split_updates=.false.
        if (myid.eq.0) then
           write(use_unit,'(2X,A)') 'From now on: orbital based force update selected'
        end if
     else
        split_updates=.true.
        if (myid.eq.0) then
           write(use_unit,'(2X,A)') 'From now on: density matrix based force update selected'
        end if
     end if

     autoselect_density_method = .false. ! autoselection should only be called once during a calculation
  end if

  if( allocated(i_basis))then
     deallocate(i_basis)
  end if
  if( allocated(i_basis_small))then
     deallocate(i_basis_small)
  end if
  if (allocated(my_functions)) then
     deallocate(my_functions)
  end if
  
  got_n_compute_maxes = .true.
  
end subroutine get_n_compute_maxes_p1
!******
