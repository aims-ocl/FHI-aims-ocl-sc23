!****s* FHI-aims/output_molecular_dynamics
!  NAME
!    output_molecular_dynamics
!  SYNOPSIS
subroutine output_molecular_dynamics(free_energy)
!  PURPOSE
!    output routine for molecular dynamics
!    note that this routine outputs all the information for the LAST set of coordinates (not the freshly predicted ones!)
!    this is due to the fact that the generalized leap-frog algorithms don't determine the velocity at a given point until 
!    they actually predict the NEXT set of coordinates. 
!  USES
  use constants, only: boltzmann_kB, MD_KE_factor, bohr, hartree
  use molecular_dynamics, only: v_half, r_last, v_last, BDP_conint, BDP_factor,&
      langham_last, MD_Epot_last, MD_H0, MD_high_order_i, &
      MD_max_high_order_index, n_atoms_MD, s_dot_NP_last, s_NP_last, tsystem, &
      tsystem_last, v_last, counter, calculate_kinetic_energy, &
      calculate_kinetic_energy_nose
  use dimensions, only: compute_heat_flux, n_atoms, MD_QH_init_segments, &
      n_periodic, use_harmonic_pot_only, use_MD_QH_init, use_reffree_AS, &
      use_thermodynamic_integration
  use species_data, only: species_pseudoized, species_name, species_m
  use geometry, only: species, coords
  use runtime_choices, only: MD_ensemble, MD_g_DOF, MD_Q_NP, MD_temperature, &
      MD_time, output_level, MD_QH_first_atom, MD_QH_last_atom
  use localorb_io, only: localorb_info, use_unit
  use timing, only: MD_stepcount
  use thermodynamic_integration, only: output_tdi_statistics
  use heat_flux, only: assemble_heat_flux, print_atomic_stress_and_heat_flux
  implicit none
!  ARGUMENTS
  real*8 :: free_energy
!  INPUTS
!    free_energy - energy from scf solver 
!          All other inputs are scavenged from the various modules in which they reside
!  OUTPUT
!    none
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
!
  character*120 :: info_str
  character*5 :: andersen_update
  real*8 :: kinetic_energy, temperature, nose_Hamiltonian, nose_poincare_hamiltonian
  integer :: i_atom, i_coord, i_coord2, i_step,i_segment
  real*8, dimension(3,n_atoms) :: coords_temp, v_out, r_out
  real*8, dimension(3) :: length
  logical :: full_MD_output

  write(info_str,'(A)') "------------------------------------------------------------"
  call localorb_info(info_str)
  
  ! check at which point in the MD step we are, if it ain't the end point, do a separate output of intermediate quantities
  if ((MD_ensemble.eq.'NVE_4th_order').and.(MD_high_order_i.ne.1)) then  ! on 1 we know all of the previous coords, do an output then!
     full_MD_output = .false.
     if (MD_high_order_i.eq.0) then
        i_step = MD_max_high_order_index-1
     else
        i_step = MD_high_order_i
     end if
     write(info_str,'(2X,A,I1,A,I1,A)') "Intermediate step number ",i_step," of ",MD_max_high_order_index-1,& 
                                        " in high order integration of MD step"
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A)') "| This output does NOT correspond to a full time step and should NOT be used for data analysis!"
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A)') "| The current coordinates and velocities are given only to keep you informed of what is going on."
     call localorb_info(info_str,use_unit,'(A)')
     r_out(:,:) = coords(:,:)
     v_out(:,:) = v_half(:,:)
  else
     full_MD_output = .true.
  end if

  ! only do this if we are at the end of one full MD step
  if (use_thermodynamic_integration) then
    if ( use_harmonic_pot_only ) then
      write(info_str,'(A)') "  Harmonic Potential: "
      call localorb_info(info_str)
      write(info_str,'(A)') "  All potential energy values below have been calculated on the basis of the harmonic potential alone."
      call localorb_info(info_str)
      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info(info_str)
    else
      write(info_str,'(A)') "  Thermodynamic integration:"
      call localorb_info(info_str)
      write(info_str,'(A)') "  All values given in the MD section refer to the hybrid system"
      call localorb_info(info_str)
      write(info_str,'(A)') "------------------------------------------------------------"
      call localorb_info(info_str)
    end if
  end if
  if (full_MD_output) then
     if (MD_stepcount.eq.0) then
        write(info_str,'(2X,A)') 'Initial conditions for Born-Oppenheimer Molecular Dynamics:'
     else 
        write(info_str,'(2X,A)') 'Advancing structure using Born-Oppenheimer Molecular Dynamics:'
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)') 'Complete information for previous time-step:'
     end if
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A,I8)')                "| Time step number          : ", MD_stepcount-1
     call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(2X,A,1X,E30.15,A)')      "| Simulation time           : ", tsystem_last," ps"
     call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(2X,A,1X,E30.15,A)')      "| Electronic free energy    : ", MD_Epot_last * hartree, " eV"
     call localorb_info(info_str,use_unit,'(A)')     
     if ((MD_ensemble.eq."NVE"           ) .or. &
         (MD_ensemble.eq."NVE_4th_order" ) .or. &
         (MD_ensemble.eq."NVE_damped"    ) .or. &
         (MD_ensemble.eq."NVT_andersen"  ) .or. &
         (MD_ensemble.eq."NVT_parrinello") .or. &
         (MD_ensemble.eq."GLE_thermostat") .or. &
         (MD_ensemble.eq."NVT_berendsen" )) then
        call calculate_kinetic_energy(v_last,kinetic_energy)
        if ((use_MD_QH_init).and.(MD_QH_init_segments .gt. 1)) then
          write (info_str,'(2X,A)')               "| Temperature (nuclei)      : "
          call localorb_info(info_str,use_unit,'(A)')

          temperature = 2d0*kinetic_energy/(dble(3*n_atoms_MD)*boltzmann_kB)
          write (info_str,'(2X,A,1X,E30.15,A)')   "|  - Overall                : ", temperature, " K"
          call localorb_info(info_str,use_unit,'(A)')
          
          do i_segment=1,MD_QH_init_segments
            call calculate_kinetic_energy(v_last,kinetic_energy,MD_QH_first_atom(i_segment),MD_QH_last_atom(i_segment))
            temperature = 2d0*kinetic_energy/(dble(3*(MD_QH_last_atom(i_segment)-MD_QH_first_atom(i_segment)+1))*boltzmann_kB)
            write (info_str,'(2X,A,1X,I4,A,1X,E30.15,A)')   "|  - Segment         ",i_segment,    "  : ", temperature, " K"
            call localorb_info(info_str,use_unit,'(A)')
          end do
          call calculate_kinetic_energy(v_last,kinetic_energy)
        else
          temperature = 2d0*kinetic_energy/(dble(3*n_atoms_MD)*boltzmann_kB)
          write (info_str,'(2X,A,1X,E30.15,A)')   "| Temperature (nuclei)      : ", temperature, " K"
          call localorb_info(info_str,use_unit,'(A)')
        end if
        write (info_str,'(2X,A,1X,E30.15,A)')   "| Nuclear kinetic energy    : ", kinetic_energy * hartree, " eV"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(2X,A,1X,E30.15,A)')   "| Total energy (el.+nuc.)   : ", (kinetic_energy+MD_Epot_last) * hartree, " eV"
        call localorb_info(info_str,use_unit,'(A)')
	! print conserved quantity in Bussi-Donadio-Parrinello Thermostat
        if (MD_ensemble.eq.'NVT_parrinello') then
		write(info_str,'(2X,A,1X,E30.15,A)') "| BDP pseudo-Hamiltonian    : ", (kinetic_energy+MD_Epot_last+BDP_conint)*hartree, " eV"
           call localorb_info(info_str,use_unit,'(A)')
		write(info_str,'(2X,A,1X,E30.15)')   "| BDP vel. scaling factor   : ", BDP_factor 
           call localorb_info(info_str,use_unit,'(A)')
        end if
        if (MD_ensemble.eq.'GLE_thermostat') then
		write(info_str,'(2X,A,1X,E30.15,A)') "| GLE pseudo hamiltonian    : ", (kinetic_energy+MD_Epot_last+langham_last)*hartree, " eV"
           call localorb_info(info_str,use_unit,'(A)')
        end if        
     else if ((MD_ensemble.eq."NVT_nose-poincare")  .or. &
              (MD_ensemble.eq."NVT_nose-hoover"  )) then
        call calculate_kinetic_energy_nose(v_last,1d0,kinetic_energy)
        temperature = 2d0*kinetic_energy/(dble(3*n_atoms_MD)*boltzmann_kB)
        write (info_str,'(2X,A,1X,E30.15,A)')   "| Temperature (nuclei)      : ", temperature, " K"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(2X,A,1X,E30.15,A)')   "| Nuclear kinetic energy    : ", kinetic_energy * hartree, " eV"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(2X,A,1X,E30.15,A)')   "| Total energy (el.+nuc.)   : ", (kinetic_energy+MD_Epot_last) * hartree, " eV"
        call localorb_info(info_str,use_unit,'(A)')
        if (MD_ensemble.eq.'NVT_nose-hoover') then
           write (info_str,'(2X,A,1X,E30.15)')  "| Thermostat value ln(s)    : ", s_NP_last
        else 
           write (info_str,'(2X,A,1X,E30.15)')  "| Thermostat value    s     : ", s_NP_last
        end if
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A,1X,E30.15,A)')    "| Thermostat momentum pi    : ", s_dot_NP_last*hartree*MD_KE_factor," eV ps"
        call localorb_info(info_str,use_unit,'(A)')
        ! calculate value of Nosé-Hoover Hamiltonian, this should be conserved
        if (MD_ensemble.eq.'NVT_nose-hoover') then
           nose_Hamiltonian = (kinetic_energy + &
                MD_Epot_last   + &
                MD_KE_factor*s_dot_NP_last*s_dot_NP_last/(2d0*MD_Q_NP)+ &
                MD_temperature*dble(MD_g_DOF)*boltzmann_kB*s_NP_last)*hartree
           write(info_str,'(2X,A,1X,E30.15,A)') "| Nose-Hoover Hamiltonian   : ", nose_Hamiltonian, " eV"
           call localorb_info(info_str,use_unit,'(A)')
        end if
        if (MD_ensemble.eq.'NVT_nose-poincare') then         
           nose_poincare_Hamiltonian = s_NP_last*(kinetic_energy/(s_NP_last*s_NP_last) + &
                MD_Epot_last   + &
                MD_KE_factor*s_dot_NP_last*s_dot_NP_last/(2d0*MD_Q_NP)+ &
                MD_temperature*dble(MD_g_DOF)*boltzmann_kB*log(s_NP_last)-MD_H0)*hartree
           write(info_str,'(2X,A,1X,E30.15,A)')"| Nose-Poincare Hamiltonian : ", nose_poincare_Hamiltonian, " eV"
           call localorb_info(info_str,use_unit,'(A)')
        end if
     end if
     r_out(:,:) = r_last(:,:)
     v_out(:,:) = v_last(:,:)
  end if
  write(info_str,'(A)') "------------------------------------------------------------"
  call localorb_info(info_str)

  ! Output for the thermodynamic integration
  if (use_thermodynamic_integration .or. use_reffree_AS) then
    call output_TDI_statistics(kinetic_energy)
  end if
  
  !! OLD BACKFOLDING
  !if (n_periodic.gt.0) then
  !   ! move atoms in back in original unit cell ...
  !   ! use same technology as map_to_center_cell - with slightly shifted boundaries.
  !   ! moreover, this is not time critical, so if statements don't really matter
  !   coords_temp(:,:) = 0d0
  !   do i_atom = 1, n_atoms
  !      length(:) = 0d0
  !      do i_coord = 1, 3
  !         ! Expand coord in lattice vectors
  !         do i_coord2 = 1, 3
  !            length(i_coord) = length(i_coord) + &
  !                 map_to_center_cell_matrix(i_coord,i_coord2)*r_out(i_coord2,i_atom)
  !            !CC: The original & wrong line was:
  !            !length(i_coord) = length(i_coord) + &
  !	       !    map_to_center_cell_matrix(i_coord2,i_coord)*r_out(i_coord,i_atom)
  !         end do
  !         ! transform the lattice vectors to central Wigner Seitz cell, 
  !         ! first to [-1,1]
  !         length(i_coord) = length(i_coord)-dble(int(length(i_coord)))
  !         ! then to [0,1]
  !         if (length(i_coord).lt.0d0) length(i_coord) = length(i_coord) + 1d0
  !      end do
  !      ! and then back to lattice coordinates
  !      do i_coord = 1, 3
  !         do i_coord2 = 1, 3
  !            coords_temp(i_coord,i_atom) = coords_temp(i_coord,i_atom) + &
  !                 lattice_vector(i_coord,i_coord2)*length(i_coord2)
  !         end do
  !      end do
  !   end do
  !else
  !   coords_temp(:,:) = r_out(:,:) 
  !end if

  !CC New Backfolding:
  if (n_periodic.gt.0) then
     ! move atoms in back in original unit cell ...
     ! use same technology as map_to_center_cell - with slightly shifted boundaries.
     ! moreover, this is not time critical, so if statements don't really matter
     coords_temp(:,:) = r_out(:,:) 
     do i_atom = 1, n_atoms
       call map_to_first_octant(coords_temp(:,i_atom))
     end do
  else
     coords_temp(:,:) = r_out(:,:) 
  end if

  
  if (.not.full_MD_output) then
     write(info_str,'(2X,A)') 'Intermediate atomic structure used for the force and energy evaluation: '     
  else
     if (tsystem.lt.MD_time) then 
       write(info_str,'(2X,A)') 'Atomic structure (and velocities) as used in the preceding time step: '
     else
       write(info_str,'(2X,A)') 'Final atomic structure (and velocities) as used in the preceding time step: '
     end if
  end if
  call localorb_info(info_str,use_unit,'(A)')
  if (MD_ensemble.eq.'NVT_andersen') then
     write(info_str,'(22X,A,13X,A,13X,A,7X,A,8X,A)') "x [A]","y [A]","z [A]", "Atom", 'stochastic thermostat update'
  else
     write(info_str,'(22X,A,13X,A,13X,A,7X,A)') "x [A]","y [A]","z [A]", "Atom"
  end if
  call localorb_info(info_str,use_unit,'(A)')
  do i_atom = 1, n_atoms
     if(species_pseudoized(species(i_atom))) then
     write(info_str,'(A,3(2X,F16.8),2X,A)') '  pseudocore  ',& 
                                            (coords_temp(i_coord,i_atom)*bohr, i_coord = 1, 3),species_name(species(i_atom))

     else 
     write(info_str,'(A,3(2X,F16.8),2X,A)') '   atom       ',& 
                                            (coords_temp(i_coord,i_atom)*bohr, i_coord = 1, 3),species_name(species(i_atom))
     endif
     call localorb_info(info_str,use_unit,'(A)')
     if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
        write(info_str,'(A,3(2X,F18.8))'     ) '     velocity ',& 
                                               (v_out(i_coord,i_atom)*bohr/species_m(species(i_atom)), i_coord = 1, 3)
     else if (MD_ensemble.eq.'NVT_andersen') then
        if (counter(i_atom).eq.1) then             !Benedikt, set a scratch counter for velocity updates
           andersen_update = 'yes'
        else
           andersen_update = 'no'
        end if
        write(info_str,'(A,3(2X,F18.8),7X,A)'     ) '     velocity ',(v_out(i_coord,i_atom)*bohr, i_coord = 1, 3),andersen_update
!        call localorb_info(info_str,use_unit,'(A)')
!Benedikt, check velocities        write(info_str,'(A,F18.8,2(2X,F18.8),7X,A)'     ) '     velocity_c ',(v_comparison(i_coord,i_atom)*bohr, i_coord = 1, 3)
     else
        write(info_str,'(A,3(2X,F18.8))'     ) '     velocity ',(v_out(i_coord,i_atom)*bohr, i_coord = 1, 3)
     end if
     call localorb_info(info_str,use_unit,'(A)')
  end do
  write(info_str,'(A)') "------------------------------------------------------------"
  call localorb_info(info_str)     
  
  ! CC: Print updated MD structure as well, in spite of the fact that the associated
  !     velocities are not yet known.
  if((output_level .ne. 'MD_light').and.(tsystem.lt.MD_time)) then
    write(info_str,'(2X,A)') 'Structure to be used in the next time step: '
    call localorb_info(info_str)     

    write(info_str,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') "|","Atom ","x [A]","y [A]","z [A]"
    call localorb_info(info_str,use_unit,'(A)')

    !CC New Backfolding:
    if (n_periodic.gt.0) then
       ! move atoms in back in original unit cell ...
       ! use same technology as map_to_center_cell - with slightly shifted boundaries.
       ! moreover, this is not time critical, so if statements don't really matter
       coords_temp(:,:) = coords(:,:) 
       do i_atom = 1, n_atoms
         call map_to_first_octant(coords_temp(:,i_atom))
       end do
    else
       coords_temp(:,:) = coords(:,:) 
    end if

    do i_atom = 1, n_atoms
      write(info_str,'(2X,A1,I5,A,A2,3(2X,F15.6))')    & 
        & "|",i_atom, ": Species ", species_name(species(i_atom)), (coords_temp(i_coord,i_atom)*bohr, i_coord=1,3,1)
      call localorb_info(info_str,use_unit,'(A)')
    end do

    write(info_str,'(A)') "------------------------------------------------------------"
    call localorb_info(info_str)     

  end if

  ! CC: Compute and output the heat flux
  if (compute_heat_flux) then
    call assemble_heat_flux(v_out)
    call print_atomic_stress_and_heat_flux()
  end if
  
end subroutine output_molecular_dynamics
!******
