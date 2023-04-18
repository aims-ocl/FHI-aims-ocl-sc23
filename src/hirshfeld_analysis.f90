!****s* FHI-aims/hirshfeld_analysis
!  NAME
!   hirshfeld_analysis
!  SYNOPSIS

subroutine hirshfeld_analysis ( )

!  PURPOSE
! Subroutine hirshfeld_analysis
!
! computes a Hirshfeld-type charge-fragment analysis as defined in
!
! Theor. Chim. Acta (Berl.) 44, 129-138 (1977)
!
! For a completely unique prescription, consider also "Hirshfeld-I",
! where the partial charges of the fragments used in the weighting function
! are iterated until they are equal to the partial charges that emerge from
! the Hirshfeld analysis.
!
! The output written by the present routine are partial charges, dipole moments, 
! and quadrupole moments on individual atoms.
! 
!
! AT -- includes Hirshfeld volume partitioning
!
!  USES

  use dimensions
  use runtime_choices
  use constants
  use grids
  use geometry
  use pbc_lists
  use free_atoms
  use physics
  use spline
  use vdw_correction
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use species_data
  use python_interface, only: run_python_hook, python_hooks

  implicit none

!  ARGUMENTS
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

  ! local variables

  character*100 :: info_str
  character*100 :: tmp_str

  real*8, dimension(:), allocatable   :: free_int 
  real*8, dimension(:), allocatable   :: hirsh_int
  real*8, dimension(:), allocatable   :: hirshfeld_spin_moment
  real*8, dimension(:,:), allocatable :: hirshfeld_dipole
  real*8, dimension(:,:), allocatable :: hirshfeld_quadrupole
  real*8, dimension(:,:,:), allocatable :: reference_rho_spl
  real*8, dimension(3) :: coord_current

  real*8, dimension(n_species) :: r_grid_min_sq

  real*8 dist_tab_sq(n_centers_basis_integrals)
  real*8 dist_tab(n_centers_basis_integrals)
  real*8 dir_tab(3,n_centers_basis_integrals)
  real*8 dir_tab_norm(3,n_centers_basis_integrals)
  real*8 i_r(n_centers_basis_integrals)


  real*8 :: aux_dens
  real*8 :: full_dens
  real*8 :: partition_norm
  real*8 :: hirshfeld_partition

  real*8 :: current_i_r
  real*8 :: temp_rho_new
  real*8 deformation_density(n_spin)

  real*8 :: dipole_moment

  logical :: point_on_atom

  ! counters

  integer :: i_atom

  integer :: i_my_batch
  integer :: i_full_points
  integer :: i_index

  integer :: i_center_L
  integer :: i_center

  integer :: i_spin
  integer :: i_coord
  integer :: i_coord_2

  integer :: i_sec_moment

  integer :: output_priority_old

  ! shorthand for identification of current grid point in batches
  integer :: current_atom, current_radial, current_angular
  

  ! begin work

  ! If output was requested under all circumstances, must reduce the
  ! output priority here. Must reset at end of subroutine!
  if (out_hirshfeld_always) then
     output_priority_old = output_priority
     if (output_level.eq.'MD_light') output_priority = 1
  end if

  write(info_str,'(2X,A)') "Performing Hirshfeld analysis of fragment charges and moments."
  call localorb_info( info_str, use_unit,'(A)',OL_norm)

  ! Error trap - cannot perform consistent Hirshfeld analysis for certain functionals yet!
  !
  ! PBE0, HSE, B3LYP (1,7,10) are enabled although the atomic reference densities are LDA/GGA
  ! Hartree-Fock not yet enabled.
  !
  if ( .not. ( (flag_xc.eq.1) .or. (flag_xc.eq.3) .or. (flag_xc.eq.6) .or. &
               (flag_xc.eq.7) .or. (flag_xc.eq.8) .or. (flag_xc.eq.9) .or. & 
               (flag_xc.eq.10) .or. (flag_xc.eq.11) .or. (flag_xc.eq.12) .or. &
               (flag_xc.eq.5)  .or. (flag_xc.eq.14) .or. (flag_xc.eq.4) &  !SAG 
               .or. (flag_xc == 17) &
            ) & 
     ) then

    write(info_str,'(1X,A)') & 
    "* Error: Cannot perform a consistent Hirshfeld analysis for your chosen XC functional."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    "* Reason: We do not calculate consistent atomic reference densities for this functional."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    "* If you think that this error message does not apply and wish to compute "
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    "* Hirshfeld charges anyway, please disable it in the source code and recompile."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    "* Your other option is to actually code the calculation of consistent "
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    "* reference atom densities for this functional. That would be much appreciated!"
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    " "
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') & 
    "* Leaving Hirshfeld charge determination without doing anything."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    !With the logical operator below set to false we guarantee a vdW energy and a vdW potential set to 0
    hirsh_analysis =  .false.
    return
  end if 

! Warning for B3LYP using BLYP densities
if (flag_xc.eq.10 .or. flag_xc == 17) then
    write(info_str,'(1X,A)') &
    "* Warning: The Hirshfeld analysis for your chosen XC functional is not consistent."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') &
    "* Reason: We do not calculate consistent atomic reference densities for this functional."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    select case (flag_xc)
    case (10)
        write(info_str,'(1X,A)') &
            "* B3LYP is calculated with BLYP densities. "
    case (17)
        write(info_str,'(1X,A)') &
            "* PBEsol is calculated with PBE densities. "
    end select
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(1X,A)') &
    "* Errors might be bigger, please check."
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
endif


  ! allocate needed arrays

  if (.not. allocated(hirshfeld_charge)) allocate ( hirshfeld_charge(n_atoms) )
  allocate ( free_int(n_atoms) )
  allocate ( hirsh_int(n_atoms) )
  allocate ( hirshfeld_dipole(3,n_atoms) )
  allocate ( hirshfeld_quadrupole(6,n_atoms) )
  if (spin_treatment .eq. 1) then
    allocate ( hirshfeld_spin_moment(n_atoms) )
    hirshfeld_spin_moment = 0.d0
  else ! dummy allocation to satisfy -check pointers
    allocate ( hirshfeld_spin_moment(1) )
    hirshfeld_spin_moment = 0.d0
  endif

  hirshfeld_charge     = 0.d0
  free_int             = 0.d0
  hirsh_int            = 0.d0
  hirshfeld_dipole     = 0.d0
  hirshfeld_quadrupole = 0.d0
  hirsh_analysis =  .true.
  ! Initialize charges for the Hirshfeld scheme. In preparation for a possible
  ! self-consistent Hirshfeld-I scheme later, we use temporary arrays for the
  ! "stockholder" atomic fragment densities and do not use evaluate_partition ().

  ! reference_rho_spl is not distributed across CPU's but should (and can) be distributed
  ! in exactly the same fashion as delta_v_hartree_part_spl in update_hartree_potential_p1, if needed.

  if(.not. use_distributed_spline_storage) then
     allocate (reference_rho_spl(n_max_spline,n_max_grid, n_atoms))
     do i_atom = 1, n_atoms, 1
        
        reference_rho_spl(:,:,i_atom) = free_rho_spl(:,:,species(i_atom))
        
     end do
  end if
  
  ! begin loop over batches of integration points
  i_full_points = 0

  ! need this for efficient evaluation of whether or not a point is "on an atom"
  r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)

  do i_my_batch = 1, n_my_batches, 1
     
        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_full_points = i_full_points + 1

           coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

           call tab_atom_centered_coords_p0 &
           ( coord_current, &
             dist_tab_sq, &
             dir_tab, &
             n_centers_basis_integrals, centers_basis_integrals )
           
           point_on_atom = .false.
           do i_center = 1, n_centers_basis_integrals, 1
              ! need to do exactly the same check as what is being done in the initialization of the partition tab:
              ! a point is defined to be "on an atom" when it is inside the innermost logarithmic grid shell of that atom
              ! - at this point it is impossible to spline the density of the logarithmic grid. 
              if ( dist_tab_sq(i_center).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
                 point_on_atom = .true.
                 exit ! exit do loop
              end if
           end do           

           ! avoid any evaluations on grid points that lie on another atom's nucleus
           if (.not.point_on_atom) then

              current_atom = batches(i_my_batch) % points(i_index) % index_atom
              current_radial = batches(i_my_batch) % points(i_index) % index_radial
              current_angular = batches(i_my_batch) % points(i_index) % index_angular

              if (empty(current_atom)) cycle

!             tabulate current integration point as it appears from spherical
!             coordinates centered at each atom 
              call tab_global_geometry_p0 &
             ( dist_tab_sq, &
               dir_tab, &
               dist_tab, &
               i_r, &
               dir_tab_norm, &
               n_centers_basis_integrals, centers_basis_integrals )

             ! determine the appropriate Hirshfeld partition weight explicitly.

             partition_norm = 0.d0

             if(.not. use_distributed_spline_storage) then
                do i_center_L = 1, n_centers_basis_integrals, 1
                  
                   i_center = centers_basis_integrals(i_center_L)

                   if (empty(center_to_atom(i_center))) cycle

                   aux_dens = val_spline & 
                        ( i_r(i_center_L), reference_rho_spl(1,1,center_to_atom(i_center)), & 
                        n_grid(species_center(i_center)) )

                   partition_norm = partition_norm + aux_dens

                enddo


                ! present atom
                ! recompute i_r because I am not sure that the periodic lists retain the order 
                ! that the atoms in the initial unit cell are the first n_atoms atoms in each list.

                current_i_r = invert_log_grid & 
                     ( r_radial(current_radial, species(current_atom)), &
                     r_grid_min(species(current_atom)), r_grid_inc(species(current_atom)) )
                
                aux_dens = val_spline & 
                     ( current_i_r, reference_rho_spl(1,1,current_atom), & 
                     n_grid(species(current_atom)) )


             else
                ! Comment, VB: 
                ! This is a special case where the general reference_rho_spl
                ! (a per-atom quantity, introduced to later allow more sophisticated
                ! Hirshfeld-type analyses, e.g., Hirshfeld_I ,etc - where one needs per-atom
                ! partitioning densities) cannot be simply used, since teh spline array would
                ! get too large for large structures (100-1000 atoms) on the BlueGene (512 MB per thread, max.). 
                ! Thus, we only use the splined free-atom density, which is a per-species quantity, not per atom.
                !
                ! If anyone wants to make this work with reference_rho_spl per atom, that array would have to 
                ! be distributed across CPU's, and accessed in a very similar way as the splines in
                ! sum_up_whole_potential () - see that subroutine for details. Not hard to do, we know exactly how
                ! to do it, but someone would have to take up the rewrite & testing effort ...

                do i_center_L = 1, n_centers_basis_integrals, 1
                  
                   i_center = centers_basis_integrals(i_center_L)

                   if (empty(center_to_atom(i_center))) cycle

                   aux_dens = val_spline & 
                        ( i_r(i_center_L),  free_rho_spl(1,1,species_center(i_center)), & 
                        n_grid(species_center(i_center)) )

                   partition_norm = partition_norm + aux_dens

                enddo

                ! present atom
                ! recompute i_r because I am not sure that the periodic lists retain the order 
                ! that the atoms in the initial unit cell are the first n_atoms atoms in each list.
                
                current_i_r = invert_log_grid & 
                     ( r_radial(current_radial, species(current_atom)), &
                     r_grid_min(species(current_atom)), r_grid_inc(species(current_atom)) )
                
                aux_dens = val_spline & 
                     ( current_i_r,  free_rho_spl(1,1,species(current_atom)), &
                     n_grid(species(current_atom)) )

             end if


             ! "naked" Hirshfeld weight at this point
             hirshfeld_partition = aux_dens / partition_norm
  
             ! ... and "dressed" with additional integration weights
             hirshfeld_partition = & 
             hirshfeld_partition * w_radial(current_radial, species(current_atom)) & 
             * (r_radial(current_radial, species(current_atom))**2.d0) & 
             * w_angular(current_angular, current_radial, species(current_atom)) & 
             * 4.d0 * pi

             ! Hirshfeld partition weight determined. Can now finally move on to actually use it. 
  
             ! temp_rho_new will be the difference between the actual charge density and the
             ! Hirshfeld superposition of atomic fragment reference densities
             temp_rho_new = 0.d0
             do i_spin = 1, n_spin, 1
               temp_rho_new = temp_rho_new + rho(i_spin, i_full_points)
             enddo
             full_dens = temp_rho_new

             temp_rho_new = temp_rho_new - pi4_inv * partition_norm

             ! Integrate Hirshfeld partial charge
             hirshfeld_charge(current_atom) = hirshfeld_charge(current_atom) &
               - hirshfeld_partition * temp_rho_new

           
             ! AT -- Hirshfeld volume partitioning
             free_int(current_atom) = free_int(current_atom) + &
                       r_radial(current_radial, species(current_atom))**3.d0 * & 
                       aux_dens &
                       *  w_radial(current_radial, species(current_atom)) &
                       * (r_radial(current_radial, species(current_atom))**2.d0) &
                       * w_angular(current_angular, current_radial,species(current_atom)) 

             hirsh_int(current_atom) = hirsh_int(current_atom) + &
                      r_radial(current_radial,species(current_atom))**3.d0 * &
                      hirshfeld_partition * full_dens 


             ! fragment dipole moment components
             do i_coord = 1,3,1
               hirshfeld_dipole(i_coord, current_atom) = hirshfeld_dipole(i_coord, current_atom) &
                 - hirshfeld_partition * temp_rho_new * r_radial(current_radial, species(current_atom))   &
                   * r_angular(i_coord, current_angular, current_radial, species(current_atom))
                 
             enddo
!test
!               if (current_atom.eq.2) then
!                 write(use_unit,'(3(2X,F15.8))') hirshfeld_dipole(:,current_atom)
!               end if
!test end


             ! fragment charge second moments
             i_sec_moment = 0
             do i_coord = 1,3,1
               do i_coord_2 = 1, i_coord, 1
                 i_sec_moment = i_sec_moment + 1
                 hirshfeld_quadrupole(i_sec_moment, current_atom) = hirshfeld_quadrupole(i_sec_moment, current_atom) &
                   - hirshfeld_partition * temp_rho_new * (r_radial(current_radial, species(current_atom))**2.d0)    &
                     * r_angular(i_coord, current_angular, current_radial, species(current_atom))                    &
                     * r_angular(i_coord_2, current_angular, current_radial, species(current_atom))
               enddo      
             enddo

! for spinpolarised calculations: compute spin-moment per atom  
         if (spin_treatment .eq. 1) then
        deformation_density = 0.d0
        do i_spin = 1, n_spin, 1
            deformation_density(i_spin) = rho(i_spin, i_full_points) - free_rho_spl(1,1,species(current_atom))/2.0
        enddo   
        hirshfeld_spin_moment(current_atom) = hirshfeld_spin_moment(current_atom) + hirshfeld_partition &
            * ( deformation_density(1) - deformation_density(2) )
         endif
    
           ! end if point_on_atom
           end if 

        !    end loop over a batch
        end do
     !     end distribution over threads
     !end if
  !     end loop over batches
  end do

  ! Synchronize integrals from different threads.
  call sync_hirshfeld ( hirshfeld_charge, hirshfeld_spin_moment, hirshfeld_dipole, hirshfeld_quadrupole, &
                        free_int, hirsh_int )


  if (use_vdw_correction_hirshfeld.or.use_mbd_old &
      .or. use_mbd_dev .or. use_libmbd) then
     do i_atom = 1, n_atoms, 1
        hirshfeld_volume(i_atom) = hirsh_int(i_atom) / free_int(i_atom)
     enddo
  end if
  !for the self-consistent vdW we need the free volumes and the hirshfeld weights, the three vectors are allocated in the vdw_correction.f90 routine
  if (use_vdw_correction_hirshfeld_sc .or. use_mbd_std) then
     do i_atom = 1, n_atoms, 1
        hirshfeld_volume(i_atom) = hirsh_int(i_atom) / free_int(i_atom)
        freeintegral(i_atom) = free_int(i_atom)
        hirshfeldw(i_atom) = hirsh_int(i_atom)
     enddo
  end if

  write(info_str,'(2X,A)') "----------------------------------------------------------------------"
  call localorb_info( info_str, use_unit,'(A)',OL_norm)
 
  do i_atom = 1, n_atoms, 1
  if(.not.species_pseudoized(species(i_atom))) then
    write(info_str,'(2X,A,I5,A,A9)') "| Atom ", i_atom, ": ", & 
      species_name(species(i_atom))
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8)') "|   Hirshfeld charge        : ", &
      hirshfeld_charge(i_atom)
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    if (use_pimd_wrapper .and. ipi_hirshfeld) then
       write(tmp_str,'(1X,A,I5,A,A9)') "| Atom ", i_atom, ": ", &
         species_name(species(i_atom))
       comm_string=trim(comm_string) // trim(tmp_str)
       write(tmp_str,'(1X, A,F15.8,A)') " Hirshfeld charge : ", &
         hirshfeld_charge(i_atom), " | "
       comm_string= trim(comm_string) // trim(tmp_str)
    endif

    write(info_str,'(2X,A,F15.8)') "|   Free atom volume        : ", &
      free_int(i_atom)
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8)') "|   Hirshfeld volume        : ", &
      hirsh_int(i_atom)
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    if ( spin_treatment .eq. 1) then
        write(info_str,'(2X,A,F15.8)') "|   Hirshfeld spin moment   : ", &
            hirshfeld_spin_moment(i_atom)
        call localorb_info( info_str, use_unit,'(A)',OL_norm)
    endif

    write(info_str,'(2X,A,3(F15.8,2X))') "|   Hirshfeld dipole vector : ", &
      ( hirshfeld_dipole(i_coord, i_atom) * bohr , i_coord = 1,3,1 )
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    dipole_moment = 0.d0
    do i_coord = 1,3,1
      dipole_moment = dipole_moment + hirshfeld_dipole(i_coord, i_atom)**2.d0
    enddo
    dipole_moment = sqrt(dipole_moment) * bohr
    write(info_str,'(2X,A,F15.8)') "|   Hirshfeld dipole moment : ", &
      dipole_moment
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,3(F15.8,2X))') "|   Hirshfeld second moments: ", &
      hirshfeld_quadrupole(1, i_atom) * bohr**2.d0, hirshfeld_quadrupole(2, i_atom) * bohr**2.d0, &
      hirshfeld_quadrupole(4, i_atom) * bohr**2.d0
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,24X,3(F15.8,2X))') "|     ", &
      hirshfeld_quadrupole(2, i_atom) * bohr**2.d0, hirshfeld_quadrupole(3, i_atom) * bohr**2.d0, & 
      hirshfeld_quadrupole(5, i_atom) * bohr**2.d0
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,24X,3(F15.8,2X))') "|     ", &
      hirshfeld_quadrupole(4, i_atom) * bohr**2.d0, hirshfeld_quadrupole(5, i_atom) * bohr**2.d0, & 
      hirshfeld_quadrupole(6, i_atom) * bohr**2.d0
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A)') "----------------------------------------------------------------------"
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

  endif
  enddo

!endif 
  call localorb_info ('')


  ! finally deallocate everything, well, whatever was allocated initially.
  if (allocated(free_int))              deallocate ( free_int )
  if (allocated(hirsh_int))             deallocate ( hirsh_int )
  if (allocated(hirshfeld_spin_moment)) deallocate ( hirshfeld_spin_moment )
  if (allocated(hirshfeld_dipole))      deallocate ( hirshfeld_dipole )
  if (allocated(hirshfeld_quadrupole))  deallocate ( hirshfeld_quadrupole )
  if (allocated(free_int))              deallocate ( free_int )
  if (allocated(hirsh_int))             deallocate ( hirsh_int )
  if (allocated(reference_rho_spl))     deallocate ( reference_rho_spl )

  ! If output was requested under all circumstances, must reduce the
  ! output priority for this subroutine. Must reset at end of subroutine!
  if (out_hirshfeld_always) then
     if (output_level.eq.'MD_light') output_priority = output_priority_old
  end if

  if (python_hooks%post_hirshfeld%registered) then
      call run_python_hook(python_hooks%post_hirshfeld)
  end if

end subroutine hirshfeld_analysis
!****** 
