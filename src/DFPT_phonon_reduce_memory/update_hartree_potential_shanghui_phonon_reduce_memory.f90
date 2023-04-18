!****s* FHI-aims/update_hartree_potential_shanghui_phonon_reduce_memory
!  NAME
!    update_hartree_potential_shanghui_phonon_reduce_memory
!  SYNOPSIS

subroutine update_hartree_potential_shanghui_phonon_reduce_memory &
     ( partition_tab, delta_rho, &
     delta_v_hartree_part_at_zero, &
     delta_v_hartree_deriv_l0_at_zero, &
     multipole_moments, multipole_radius_sq, &
     l_hartree_max_far_distance,       &
     outer_potential_radius )

!  PURPOSE
!  Input delta_rho, get rho_multipole in the hartree_potential_storage moldue.
!  Which can be used in the next subrution sum_up_whole_potential.f90 to get
!  delta_potential. 
!
!  The subroutine lso calls integrate_hartree_log_grid to calculate Hartree potential
!  in the form of multipole expansion. -----> I do not think
!  integrate_hartree_log_grid can calculate Hartree potential. it is only give
!  some point at ZERO.  
!  
!  shanghui 2012.06.28
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use hartree_potential_storage
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points)                 :: partition_tab  

!----------------shanghui changed here from p1 version--------
  !real*8, dimension(n_full_points)                 :: free_rho_superpos
  !real*8, dimension(n_spin, n_full_points)         :: rho
  real*8, dimension(n_full_points)                 :: delta_rho
!----------------shanghui end changed here--------------------


  real*8, dimension(n_atoms)                       :: delta_v_hartree_part_at_zero
  real*8, dimension(3,n_atoms)                     :: delta_v_hartree_deriv_l0_at_zero
  real*8, dimension( ( l_pot_max + 1)**2, n_atoms) :: multipole_moments
  real*8, dimension(n_atoms)                       :: multipole_radius_sq
  integer, dimension(n_atoms)                      :: l_hartree_max_far_distance
  real*8, dimension(0:l_pot_max, n_atoms)          :: outer_potential_radius


!  INPUTS
!   o  partition_tab -- values of partition function
!   o  delta_rho -- electron density
!
!  OUTPUT
!   o  delta_v_hartree_part_at_zero --  multipole expansion of the Hartree potential at origin of atoms
!   o  delta_v_hartree_deriv_l0_at_zero -- derivates of  multipole expansion of the Hartree potential at origin of atoms
!   o  multipole_moments -- multipole moments
!   o  multipole_radius_sq -- (radius of multipole moments)**2
!   o  l_hartree_max_far_distance -- maximum l for far distance Hartree potential (periodic systems)
!   o  outer_potential_radius -- outer radius of multipole expansion of Hartree potential, which needs spline. 
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





!  rho_multipole is the angular integral over
!  Y_lm(Omega) * partition_fn(at,r,Omega) * delta_rho(r,Omega),
!  Delley 1990, Eq. (11) 
!  i.e. rho_multipole = rho_multipole(r) . For convenience


  real*8 dir_tab(3)
  real*8 temp_rho_new
  real*8 multipole_radius
  real*8 ylm_tab((l_pot_max+1)**2)
  real*8, allocatable :: current_rho_multipole_spl(:,:,:), current_rho_multipole(:,:)

  integer l_h_dim(n_atoms)
  real*8 :: drho

!  counters

  integer i_atom
  integer i_index
  integer current_atom, current_radial, current_angular
  integer i_full_points
  integer i_spin
  integer i_my_batch
  integer i_atom_index

  ! for shmem
  real*8, allocatable :: delta_v_hartree_part_spl(:,:,:)
  integer num_shm_procs, my_shm_id, n_bytes

  character*100 :: info_str



!  begin work

  write(info_str,'(2X,A)') " "
  call localorb_info(info_str, use_unit, '(A)', OL_norm )
  write(info_str,'(2X,A,A)') "Evaluating partitioned Hartree potential by multipole expansion."
  call localorb_info(info_str, use_unit, '(A)', OL_norm )


  allocate(current_rho_multipole((l_pot_max+1)**2, n_max_radial+2),stat=i_index)
  call check_allocation(i_index, 'current_rho_multipole')

  allocate(current_rho_multipole_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2),stat=i_index)
  call check_allocation(i_index, 'current_rho_multipole_spl')


  do i_atom = 1, n_atoms, 1
    l_h_dim(i_atom) = (l_hartree(species(i_atom))+1)**2
  end do

  rho_multipole = 0.0d0
  i_full_points = 0

  do i_my_batch = 1, n_my_batches, 1
     
    do i_index = 1, batches(i_my_batch)%size, 1
           
      current_atom = batches(i_my_batch) % points(i_index) % index_atom
      current_radial = batches(i_my_batch) % points(i_index) % index_radial
      current_angular = batches(i_my_batch) % points(i_index) % index_angular

      ! run for the first time through the integration grid to tabulate 
      ! integral_1_zero_r and thus calculate integral_zero_infinity
      !
      ! meanwhile tabulate first greens-function solution of partitioned potential
      ! and difference charge density
      
      i_full_points = i_full_points + 1

      ! execute only if partition_tab.gt.0 here, i.e. if the integration point
      ! makes sense
      !!!!! This is NOT the usual integration partition tab
      !!!!! This is the HARTREE partition tab. 
      !!!!! Thus this criterion is not the same as in all other integrations.
      !!!!! We must fix that one day.
      if (partition_tab(i_full_points).gt.0.d0) then

        ! Check if the atom has been correctly dedected as mine
        i_atom_index = rho_multipole_index(current_atom)
        if(i_atom_index<=0) then
          print '(2(a,i5))','ID ',myid,' INTERNAL ERROR update_hartree_potential_p1 - need atom! Atom: ',current_atom
          call aims_stop
        endif

        ! compute atom-centered coordinates of current integration point,
        ! BUT HERE ONLY FOR PRESENT ATOM!
        dir_tab(:) = r_angular( : , current_angular, current_radial, species(current_atom))

        ylm_tab (1:l_h_dim(current_atom)) = &
             local_ylm_tab(1:l_h_dim(current_atom),current_angular, &
             lebedev_grid_index(current_radial,species(current_atom)))

        ! calculate contribution to the angular parts of the two integrals which
        ! are equal 

        ! notice that the radial and angular integration weights are already part of 
        ! partition_tab; we need not multiply them in separately.

        ! consider only difference density for hartree potential

        temp_rho_new = 0.d0
!        do i_spin = 1, n_spin, 1
!           temp_rho_new = temp_rho_new + rho(i_spin, i_full_points)
!        enddo

         temp_rho_new = delta_rho(i_full_points) 
            
!!!        rho_multipole(1, current_radial+1, i_atom_index) = &
!!!               temp_rho_new/ylm_tab(1) 
        ! implied loop over all (l,m) from (0,0) to (l_hartree,l_hartree)
        rho_multipole(1:l_h_dim(current_atom), current_radial+1, i_atom_index) = &
             rho_multipole(1:l_h_dim(current_atom), current_radial+1, i_atom_index) + &
             ylm_tab(1:l_h_dim(current_atom)) * partition_tab(i_full_points) * temp_rho_new 
         !if(current_radial.eq.1.and.i_atom_index.eq.1) then
         !write(use_unit,*) '----------------------------------'
         !write(use_unit,*) i_full_points,delta_rho(i_full_points)
         !write(use_unit,*) rho_multipole(1, current_radial+1, i_atom_index)
         !endif

      end if
    end do !    end loop over a batch
  end do !     end loop over batches


  ! Add all rho_multipole contributions

  do i_atom = 1, n_atoms

    i_atom_index = rho_multipole_index(i_atom)
    if(i_atom_index > 0) then
      current_rho_multipole(:,:) = rho_multipole(:,:,i_atom_index)
    else
      current_rho_multipole(:,:) = 0.
    endif

    call sync_vector(current_rho_multipole, ((l_pot_max+1)**2)*(n_max_radial+2))

    if(i_atom_index > 0) then
      rho_multipole(:,:,i_atom_index) = current_rho_multipole(:,:)
    endif

  enddo

  ! Set boundaries on rho_multipole

  do i_atom_index = 1, n_rho_multipole_atoms

    i_atom = i_rho_multipole_atoms(i_atom_index)

    ! At Infinity (n_radial+2)
    ! enforce zero charge at infinity explicitly by a trick:
    ! (used in splines later on)
    rho_multipole(1:l_h_dim(i_atom), n_radial(species(i_atom))+2, i_atom_index) = 0.d0

    ! At zero
    drho = (rho_multipole(1,2,i_atom_index) - rho_multipole(1,3,i_atom_index)) &
    &     /    (r_radial(1,species(i_atom)) - r_radial(2,species(i_atom)))
    if (legacy_monopole_extrapolation) then
       ! Backwards compatible but inaccurate choice of sign.
       ! Please note that the influence of this boundary rapidly decreases
       ! with increasing grid density.
       rho_multipole(1,1,i_atom_index) = rho_multipole(1,2,i_atom_index) &
       &                                 + drho * r_radial(1,species(i_atom))
    else
       rho_multipole(1,1,i_atom_index) = rho_multipole(1,2,i_atom_index) &
       &                                 - drho * r_radial(1,species(i_atom))
    end if
    rho_multipole(2:l_h_dim(i_atom),1,i_atom_index) = 0.d0

  enddo

  ! Calculate output variables
  ! the variables must be set to 0 since they are sync'd at the end
  delta_v_hartree_part_at_zero(:) = 0
  delta_v_hartree_deriv_l0_at_zero(:,:) = 0
  multipole_moments(:,:) = 0
  multipole_radius_sq(:) = 0
  l_hartree_max_far_distance(:) = 0
  outer_potential_radius(:,:) = 0

  do i_atom = 1, n_atoms

    if (mod(i_atom-1,n_tasks) == myid) then

      call get_rho_multipole_spl(current_rho_multipole_spl, i_atom)


      if (force_hartree_log_grid) then
        !write(use_unit,*) 'n_radius,in update_hartree_potential.f90',n_radial
        !write(use_unit,*) rho_multipole(1, :, 1)
        call integrate_hartree_log_grid &
                (i_atom, current_rho_multipole_spl, &
                delta_v_hartree_part_at_zero(i_atom), &
                delta_v_hartree_deriv_l0_at_zero(1:3,i_atom), &
                multipole_moments(1,i_atom), &
                multipole_radius, &
                l_hartree_max_far_distance(i_atom), &
                outer_potential_radius (0:l_pot_max, i_atom) )
        !stop 
        !RJ: for safety only
        multipole_radius = min(multipole_radius,r_radial(n_radial(species(i_atom)),species(i_atom)))
        multipole_radius_sq(i_atom) = multipole_radius**2           

      else
           
        write(use_unit,*) 'Only log integrals are supposted'

      end if

!     end distribution over threads
    end if

!     end loop over atoms 
  end do

    !------shanghui test multipole_moments
    !write(use_unit,*) 'MM',multipole_moments(1,1)
  ! synchronize results

  call sync_vector(delta_v_hartree_part_at_zero,n_atoms)
  call sync_vector(delta_v_hartree_deriv_l0_at_zero,3*n_atoms)
  call sync_vector(multipole_moments,(l_pot_max+1)**2*n_atoms)
  call sync_vector(multipole_radius_sq,n_atoms)
  call sync_integer_vector(l_hartree_max_far_distance,n_atoms)
  call sync_vector(outer_potential_radius,(l_pot_max+1)*n_atoms)

    !write(use_unit,*) 'after sync,MM',multipole_moments(1,1)

  if(communication_type.eq.shmem_comm) then

    if(use_distributed_spline_storage) then
      ! We need complete rho_multipole array (well, not really complete ...)
      print *,'use_distributed_spline_storage and shmem communication can not be used together'
      call aims_stop
    endif

    call aims_shm_n_procs(num_shm_procs)
    call aims_shm_myid(my_shm_id)

    allocate(delta_v_hartree_part_spl((l_pot_max+1)**2, 2, n_hartree_grid))

    ! Calculate spline coefficients and put them into shared memory

    n_bytes = (l_pot_max+1)**2 * 2 * n_hartree_grid * 8
    do i_atom = 1,n_atoms
      if( mod(i_atom-1,num_shm_procs) == my_shm_id) then
        call get_rho_multipole_spl(current_rho_multipole_spl, i_atom)
        call integrate_delta_v_hartree( &
             current_rho_multipole_spl, &
             delta_v_hartree_part_spl, 2, i_atom)
        call aims_shm_put(delta_v_hartree_part_spl, (i_atom-1)*n_bytes, n_bytes)
      endif
    enddo

    deallocate(delta_v_hartree_part_spl)

  endif



  ! Now that we have all multipole components (after synchronisation of all MPI threads),
  ! output a checksum of all extrapolated point charges in the Hartree potential
  write(info_str,'(2X,A)') "| Analytical far-field extrapolation by fixed multipoles:"
  call localorb_info(info_str, use_unit,'(A)', OL_norm )
  write(info_str,'(2X,A,E14.6)') "| Hartree multipole sum: apparent total charge = ", &
       - sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4)
  call localorb_info(info_str, use_unit,'(A)', OL_norm )



  deallocate(current_rho_multipole)
  deallocate(current_rho_multipole_spl)

end subroutine  update_hartree_potential_shanghui_phonon_reduce_memory
!******
