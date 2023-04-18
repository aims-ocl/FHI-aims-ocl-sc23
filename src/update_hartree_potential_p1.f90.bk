!****s* FHI-aims/update_hartree_potential_p1
!  NAME
!    update_hartree_potential_p1
!  SYNOPSIS

subroutine update_hartree_potential_p1 &
     ( partition_tab, free_rho_superpos, &
     rho, rho_multipole_old, &
     delta_v_hartree_part_at_zero, &
     delta_v_hartree_deriv_l0_at_zero, &
     multipole_moments, multipole_radius_sq, &
     l_hartree_max_far_distance, rho_radial_multipole_old, &
     outer_potential_radius, AS_stress_on, subtract_free_rho)

!  PURPOSE
!  computes the partitioned density (minus the density of free atoms) to
!  calculate the partitioned hartree potential. Calculates multipole expansion
!  from partitioned density.
!
!  The subroutine also calls integrate_hartree_log_grid to calculate Hartree potential
!  in the form of multipole expansion.
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
  use analytical_stress
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points)                 :: partition_tab  
  real*8, dimension(n_full_points)                 :: free_rho_superpos
  real*8, dimension(n_spin, n_full_points)         :: rho
  real*8, dimension(n_full_points)                 :: rho_multipole_old

  real*8, dimension(n_atoms)                       :: delta_v_hartree_part_at_zero
  real*8, dimension(3,n_atoms)                     :: delta_v_hartree_deriv_l0_at_zero
  real*8, dimension( ( l_pot_max + 1)**2, n_atoms) :: multipole_moments
  real*8, dimension(n_atoms)                       :: multipole_radius_sq
  integer, dimension(n_atoms)                      :: l_hartree_max_far_distance
  real*8, dimension((l_pot_max+1)**2, n_max_radial+2, n_atoms) :: rho_radial_multipole_old
  real*8, dimension(0:l_pot_max, n_atoms)          :: outer_potential_radius
  logical                                          :: AS_stress_on
  logical, intent(in)                              :: subtract_free_rho

!  INPUTS
!   o  partition_tab -- values of partition function
!   o  free_rho_superpos -- superposition of free atom electron densities
!   o  rho -- electron density
!   o  rho_multipole_old -- multipole expansion of electron density from previous iteration
!   o  rho_radial_multipole_old -- electron density from previous iteration
!   o  subtract_free_rho -- whether to subtract the superposition of free atom densitites from rho.
!                           If false, the multipole expansion is performed directly on rho and
!                           not delta_rho.
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

  real*8 :: multipole_moment_error

  real*8 :: leftover_charge_per_atom

!  counters

  integer i_atom
  integer i_index
  integer current_atom, current_radial, current_angular
  integer i_full_points
  integer i_spin
  integer i_my_batch
  integer i_atom_index

  integer i_l, i_m

  ! for shmem
  real*8, allocatable :: delta_v_hartree_part_spl(:,:,:)
  integer num_shm_procs, my_shm_id, n_bytes

  character*100 :: info_str



!  begin work

  write(info_str,'(2X,A)') " "
  call localorb_info( info_str, use_unit,'(A)', OL_norm )
  write(info_str,'(2X,A,A)') "Evaluating partitioned Hartree potential by multipole expansion."
  call localorb_info( info_str, use_unit,'(A)', OL_norm )


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

        ! Check if the atom has been correctly detected as mine
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
        do i_spin = 1, n_spin, 1
           temp_rho_new = temp_rho_new + rho(i_spin, i_full_points)
        enddo


        if (flag_delta_rho_in_multipole) &
             & temp_rho_new = temp_rho_new - rho_multipole_old(i_full_points)
        if (subtract_free_rho) &
             & temp_rho_new = temp_rho_new - pi4_inv * free_rho_superpos(i_full_points)
  
        ! implied loop over all (l,m) from (0,0) to (l_hartree,l_hartree)
        rho_multipole(1:l_h_dim(current_atom), current_radial+1, i_atom_index) = &
             rho_multipole(1:l_h_dim(current_atom), current_radial+1, i_atom_index) + &
             ylm_tab(1:l_h_dim(current_atom)) * partition_tab(i_full_points) * temp_rho_new 

      end if
    end do !    end loop over a batch
  end do !     end loop over batches

  if(flag_delta_rho_in_multipole)then
    do i_atom_index = 1, n_rho_multipole_atoms
      i_atom = i_rho_multipole_atoms(i_atom_index)
      rho_multipole(:,:,i_atom_index) = rho_multipole(:,:,i_atom_index) + rho_radial_multipole_old(:,:,i_atom)
      rho_radial_multipole_old(:,:,i_atom) = rho_radial_multipole_old(:,:,i_atom) * (1.d0-multipole_feed_back_parameter) &
           + multipole_feed_back_parameter * rho_multipole(:,:,i_atom_index)
    enddo
  end if
! SR: to save these and use them in other routines
!   if (solvent_method.eq.SOLVENT_MPB) then
!     rho_multipole_save_for_lpb(:, :, :)&
!       = rho_multipole(:, 2:, :)
!   end if

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

  original_multipole_moments(:,:) = 0.d0
  if (compensate_multipole_errors) then
     compensation_norm(:)   = 0.d0
     compensation_radius(:) = 0.d0
     total_compensated_charge = 0.d0
  end if

  do i_atom = 1, n_atoms

    if (mod(i_atom-1,n_tasks) == myid) then

      if (compensate_multipole_errors) then
         call get_multipole_moments_on_original_grid(i_atom)
      end if

      call get_rho_multipole_spl(current_rho_multipole_spl, i_atom)
      !! CC-DEBUG
      !! if ( AS_stress_on .and. AS_flag_read_mm_moments ) then
      !!   if ( .not.allocated(AS_mm_splines_save) ) then
      !!     allocate(AS_mm_splines_save(1:(l_pot_max+1)**2 ,  1:n_max_spline, 1:n_max_radial+2 ))
      !!     AS_mm_splines_save = current_rho_multipole_spl
      !!   end if
      !!   call AS_read_mm_splines(current_rho_multipole_spl, (l_pot_max+1)**2 ,  n_max_spline, n_max_radial+2 , i_atom )
      !! end if


      if (force_hartree_log_grid) then

        call integrate_hartree_log_grid &
                (i_atom, current_rho_multipole_spl, &
                delta_v_hartree_part_at_zero(i_atom), &
                delta_v_hartree_deriv_l0_at_zero(1:3,i_atom), &
                multipole_moments(1,i_atom), &
                multipole_radius, &
                l_hartree_max_far_distance(i_atom), &
                outer_potential_radius (0:l_pot_max, i_atom) )
    
        !RJ: for safety only
        !VB: I can not see that the following modification would make sense,
        !    hence commented out again. Time will tell whether there is a problem,
        !    but if I remember correctly, there was no specific problem that warranted
        !    this change. I see a potential way how it could do harm though.
        !
        !    multipole_radius = min(multipole_radius,r_radial(n_radial(species(i_atom)),species(i_atom)))

        multipole_radius_sq(i_atom) = multipole_radius**2           

      else
           
        write(use_unit,*) 'Only log integrals are supported'
        stop

      end if

!     end distribution over threads
    end if

!     end loop over atoms 
  end do

  ! synchronize results

  call sync_vector(delta_v_hartree_part_at_zero,n_atoms)
  call sync_vector(delta_v_hartree_deriv_l0_at_zero,3*n_atoms)
  call sync_vector(multipole_moments,(l_pot_max+1)**2*n_atoms)
  call sync_vector(multipole_radius_sq,n_atoms)
  call sync_integer_vector(l_hartree_max_far_distance,n_atoms)
  call sync_vector(outer_potential_radius,(l_pot_max+1)*n_atoms)

  if (compensate_multipole_errors) then
     call sync_vector(original_multipole_moments,(l_pot_max+1)**2*n_atoms)

     call sync_vector(compensation_norm,n_atoms)
     call sync_vector(compensation_radius,n_atoms)
     call sync_real_number(total_compensated_charge)

     write(info_str,'(2X,A,E14.6)') "| Original multipole sum: apparent total charge = ", &
       - sum( original_multipole_moments(1,1:n_atoms)) * sqrt(pi4)
     call localorb_info(info_str,use_unit,'(A)', OL_norm )

     write(info_str,'(2X,A,E14.6)') "| Sum of charges compensated after spline to logarithmic grids = ", &
       - total_compensated_charge * sqrt(pi4)
     call localorb_info(info_str,use_unit,'(A)', OL_norm )
  end if

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
  call localorb_info(info_str,use_unit,'(A)', OL_norm )
  write(info_str,'(2X,A,E14.6)') "| Hartree multipole sum: apparent total charge = ", &
       - sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4)
  call localorb_info(info_str,use_unit,'(A)', OL_norm )

  ! provide short form of apparent total charge in MD_light output
  if (output_level .eq. 'MD_light') then
     write(info_str,'(E10.2)') &
       - sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4)
     call localorb_info(info_str,use_unit,'(A,$)',OL_high)
  end if

  ! The following tag is legacy only and should be removed in a future version.
  ! As implemented, it will introduce discontinuities in the density.
  if (distribute_leftover_charge) then
      write(info_str,*) 'Apparent charge will be equally distributed among all atoms.'
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(1X,A)') '*** This is a legacy tag only. It can be used to give a defined dipole'
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(1X,A)') '*** correction in some slab calculations, but can introduce significant'
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(1X,A)') '*** total energy errors, including discontituities.'
      call localorb_info(info_str,use_unit,'(A)')
      leftover_charge_per_atom =  sum( multipole_moments(1,1:n_atoms)) /n_atoms
      do i_atom=1,n_atoms,1
          multipole_moments(1,i_atom) =  multipole_moments(1,i_atom) -leftover_charge_per_atom
      enddo
  endif
   
  if(out_zero_multipoles)then
     write(info_str,'(2X,A)') " "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A)') "| Atom-resolved charges in the electrostatic potential : "
     call localorb_info(info_str,use_unit,'(A)')
     do i_atom = 1, n_atoms, 1
        write(info_str,'(2X,A,1X,I8,A,F15.8)') '| Atom ',i_atom,' : ',-multipole_moments(1,i_atom)*sqrt(pi4)
        call localorb_info(info_str,use_unit,'(A)')
     end do
     write(info_str,'(2X,A)') " "
     call localorb_info(info_str,use_unit,'(A)')
  end if

  ! Write out all multipole components of the Hartree potential.
  ! This is mainly useful for debugging purposes, and in general is an output
  ! option we should have had ages ago.
  if(out_hartree_multipoles)then
     write(info_str,'(2X,A)') "| "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A)') "| Multipole moments used in the electrostatic potential : "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A)') "| What we write is -moment*sqrt(pi/4) , which gives charges in e etc."
     call localorb_info(info_str,use_unit,'(A)')
     do i_atom = 1, n_atoms, 1
        write(info_str,'(2X,A)') "| "
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,1X,I8,A)') '| Atom ',i_atom, ":   l   m          moment"
        call localorb_info(info_str,use_unit,'(A)')
        i_index = 0
        do i_l = 0, l_hartree(species(i_atom)), 1
           do i_m = -i_l, i_l, 1
              i_index = i_index+1
              write(info_str,'(2X,A,16X,I3,1X,I3,1X,F15.8)') '| ',i_l,i_m,-multipole_moments(i_index,i_atom)*sqrt(pi4)
              call localorb_info(info_str,use_unit,'(A)')
           end do
        end do
     end do
     write(info_str,'(2X,A)') " "
     call localorb_info(info_str,use_unit,'(A)')
  end if



  deallocate(current_rho_multipole)
  deallocate(current_rho_multipole_spl)

end subroutine update_hartree_potential_p1
!******
