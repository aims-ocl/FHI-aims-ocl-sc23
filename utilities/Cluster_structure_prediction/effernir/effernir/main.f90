program effernir

  use vector_array_class
  use basin_hopping, bh_allocate => allocate, bh_deallocate => deallocate, bh_optimize => optimize, &
       bh_rebuild_linked_list => rebuild_linked_list, bh_analyze_moves => analyze_moves
  use cluster
  use control
  use dft_file
  use arch_specific
  use geometry_list
  use geometry_class, g_allocate => allocate, g_deallocate => deallocate
  use pair_potential, pp_deallocate => deallocate
  use cartesian_ylm
  
  implicit none
  
  ! local variables
  integer :: status
  type (vector_array) :: final_coords
  real*8, dimension(:), allocatable :: final_spin_per_atom
  
  call parse_control()
  call cluster_allocate()

  call read_control()
  call read_geo("geometry.in.basic")
  
  call initialize_cartesian_ylm()

  if (potential_flag .eq. 'external') then
     call initialize_dft_file()
  end if

  call allocate_array(final_coords, n_atoms)
  if (spin_polarized) then
     allocate(final_spin_per_atom(n_atoms))
  end if
  
  select case(simulation_flag)
     
  case ('basin-hopping')

     call bh_allocate()
     call bh_optimize(final_coords, final_spin_per_atom, status)
     if (status .eq. 3) then
        call bh_analyze_moves()
     end if
     call bh_deallocate()

  case ('bh-analysis') 

     call bh_allocate()
     ! create linked_list from the dft_data.dat-file to be able to play around with
     ! norm and parameters to distinguish structures
     call bh_rebuild_linked_list()
     ! run once again through the bh_log-file with the linked list of geometries in mind
     ! to identify the moves that brought the system into another basin
     call bh_analyze_moves()
     call bh_deallocate()

  case ('sample-hop-matrix')

     call bh_allocate()
     ! create linked_list from the dft_data.dat-file to be able to play around with
     ! norm and parameters to distinguish structures
     call bh_rebuild_linked_list()
     call sample_hop_matrix(final_coords, final_spin_per_atom, status)

  case default
     
     write(6,*) simulation_flag, "?. What do you want from me??"
     write(6,*) "* Aborting."
     stop
     
  end select
  
  if (potential_flag .eq. 'external') then
     call write_status(status)
  
     if ((status .eq. 2) .or. (status .eq. 3)) then
        ! tell the wrapper to stop optimization
        open (20, FILE="finished.tmp")
        write (20, *) "Many greetings! Check the wrapper script please. File should have been deleted."
        write (20, *) "R.Gehrke"
        close(20)
     else

        ! new energy and forces needed
        call write_geometry_file(final_coords, ini_charge_per_atom, final_spin_per_atom)

     end if
     
     call cleanup_dft_file()

  end if
  
  if (allocated(final_spin_per_atom)) then
     deallocate(final_spin_per_atom)
  end if

  if (use_pair_pot) then
     call pp_deallocate()
  end if
  call cleanup_cartesian_ylm()

end program effernir
