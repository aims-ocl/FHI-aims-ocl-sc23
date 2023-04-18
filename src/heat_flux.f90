MODULE heat_flux
  IMPLICIT NONE

  save
  
  real*8,dimension(:,:,:),allocatable   :: HF_stress_per_atom_MP_CO
  real*8,dimension(:,:,:),allocatable   :: HF_stress_per_atom_MP_AT
  real*8,dimension(:,:,:),allocatable   :: HF_stress_per_atom_MP_EL
  real*8,dimension(:,:),allocatable     :: HF_stress_per_atom_PU
  real*8,dimension(:,:,:),allocatable   :: HF_stress_per_atom_vdw
  real*8,dimension(:,:,:),allocatable   :: HF_stress_per_atom_vdw_change
  real*8,dimension(:,:),allocatable     :: HF_rho_multipole_per_atom
  real*8,dimension(:,:),allocatable     :: HF_rho_free_per_atom
  real*8,dimension(:,:,:,:),allocatable :: HF_dde_rho
  real*8,dimension(:,:,:),allocatable   :: HF_stress_per_atom
  real*8,dimension(1:3)                 :: J


  contains

  ! CC: Assemble stress per atom and heat flux
  subroutine assemble_heat_flux(velocities)
   use dimensions
   use analytical_stress
   use runtime_choices
   use geometry, only: cell_volume
   use localorb_io, only: localorb_info, use_unit
   use applicable_citations, only: cite_reference
   IMPLICIT NONE
   real*8,dimension(1:3,1:n_atoms),intent(in),optional :: velocities
   integer                                             :: i_atom,i_coord
   integer                                             :: l_index(9),m_index(9),map_pu(9)
   character*400                                       :: info_str
   real*8,dimension(1:3,1:3)                           :: stress
   real*8,dimension(1:3,1:n_atoms)                     :: J_per_atom

   ! Map lm-index 
   !  1 4 5
   !  7 2 6
   !  8 9 3
   l_index    = (/ 1, 2, 3, 1, 1, 2, 2, 3, 3 /)
   m_index    = (/ 1, 2, 3, 2, 3, 3, 1, 1, 2 /)
   ! Map hessian-index 
   !  1 2 3
   !  7 4 5
   !  8 9 6
   map_pu     = (/ 1, 4, 6, 2, 3, 5, 7, 8, 9 /)

   ! (A) Compute stress per atom from individual stress components
   if (.not.allocated(HF_stress_per_atom)) allocate(HF_stress_per_atom(1:3,1:3,1:n_atoms))
   HF_stress_per_atom(:,:,:) = 0.0d0
   stress(:,:)               = 0.0d0
   do i_coord=1,AS_components,1
     HF_stress_per_atom(l_index(i_coord),m_index(i_coord),1:n_atoms) = &
       HF_stress_per_atom_PU(map_pu(i_coord),1:n_atoms)                      + &
       HF_stress_per_atom_MP_AT(l_index(i_coord),m_index(i_coord),1:n_atoms) + &
       HF_stress_per_atom_MP_EL(l_index(i_coord),m_index(i_coord),1:n_atoms) + &
       HF_stress_per_atom_MP_CO(l_index(i_coord),m_index(i_coord),1:n_atoms) 
     if (use_vdw_correction_hirshfeld ) then
       HF_stress_per_atom(l_index(i_coord),m_index(i_coord),1:n_atoms) = &
         HF_stress_per_atom(l_index(i_coord),m_index(i_coord),1:n_atoms) + &
         HF_stress_per_atom_vdw(l_index(i_coord),m_index(i_coord),1:n_atoms)
     end if
     stress(l_index(i_coord),m_index(i_coord)) = sum(HF_stress_per_atom(l_index(i_coord),m_index(i_coord),1:n_atoms))
     ! Check for concistency:
     if ( ( abs( ( stress(l_index(i_coord),m_index(i_coord))/cell_volume & 
               - analytical_stress_tensor(l_index(i_coord),m_index(i_coord)) ) &
               / analytical_stress_tensor(l_index(i_coord),m_index(i_coord)) ) .gt. 1.0d-6 ) .and. &
        ( abs(stress(l_index(i_coord),m_index(i_coord))/cell_volume) .gt. 1.0d-8 ) ) then
       write(info_str,'(2X,A,A,2I2)') "*** WARNING: Stress used for Heat Flux", &
            " is NOT consistent with the Analytical Stress in Component: ", l_index(i_coord),m_index(i_coord)  
       call localorb_info(info_str,use_unit,'(2X,A)')
       write(info_str,'(2X,A,1E30.16)') "     Relative Deviance (%) : ", &
            abs( ( stress(l_index(i_coord),m_index(i_coord))/cell_volume &
               - analytical_stress_tensor(l_index(i_coord),m_index(i_coord)) ) &
               / analytical_stress_tensor(l_index(i_coord),m_index(i_coord)) ) *100.0d0
       call localorb_info(info_str,use_unit,'(2X,A)')
       write(info_str,'(2X,A,1E30.16,A)') "     Analytical Stress : ", &
            analytical_stress_tensor(l_index(i_coord),m_index(i_coord))*Hartree/bohr**3.0d0," eV/A**3"
       call localorb_info(info_str,use_unit,'(2X,A)')
       write(info_str,'(2X,A,1E30.16,A)') "     Heat Flux  Stress : ", &
            stress(l_index(i_coord),m_index(i_coord))/cell_volume*Hartree/bohr**3.0d0," eV/A**3"
       call localorb_info(info_str,use_unit,'(2X,A)')
     end if
   end do
   ! Symmetrize ...
   if (AS_components .eq. 9) then
     HF_stress_per_atom(1,2,1:n_atoms) = 0.5d0*(HF_stress_per_atom(1,2,1:n_atoms) &
                                              + HF_stress_per_atom(2,1,1:n_atoms))
     HF_stress_per_atom(1,3,1:n_atoms) = 0.5d0*(HF_stress_per_atom(1,3,1:n_atoms) &
                                              + HF_stress_per_atom(3,1,1:n_atoms))
     HF_stress_per_atom(2,3,1:n_atoms) = 0.5d0*(HF_stress_per_atom(2,3,1:n_atoms) &
                                             + HF_stress_per_atom(3,2,1:n_atoms))
   end if
   HF_stress_per_atom(2,1,1:n_atoms) = HF_stress_per_atom(1,2,1:n_atoms)
   HF_stress_per_atom(3,1,1:n_atoms) = HF_stress_per_atom(1,3,1:n_atoms)
   HF_stress_per_atom(3,2,1:n_atoms) = HF_stress_per_atom(2,3,1:n_atoms)

   ! (B) Get diffusive heat flux from stress*v: 
   J(:) = 0.0d0
   J_per_atom(:,:) = 0.0d0
   if (present(velocities)) then
     do i_atom=1,n_atoms   
        J_per_atom(:,i_atom) = J_per_atom(:,i_atom) + &
          MATMUL( HF_stress_per_atom(:,:,i_atom),velocities(:,i_atom) )
     end do
     do i_coord=1,3,1
       J(i_coord)      = sum(J_per_atom(i_coord,1:n_atoms))
     end do
   end if

   ! (C) Deallocate to be ready for next step
   if(allocated( HF_rho_multipole_per_atom )) deallocate( HF_rho_multipole_per_atom ) 
   if(allocated( HF_rho_free_per_atom      )) deallocate( HF_rho_free_per_atom      )
   if(allocated( HF_dde_rho                )) deallocate( HF_dde_rho                )

   ! (D) Register citation:
   call cite_reference("HeatFluxGK")

  end subroutine assemble_heat_flux


  ! FlK: print atomic stress to stdout and ipi comm_string
  subroutine print_atomic_stress()
    use dimensions
    use constants
    use runtime_choices
    use localorb_io, only: localorb_info, use_unit, comm_string
    IMPLICIT NONE
    integer                             :: i_atom, i_coord
    character*400                       :: tmp_str

    write(tmp_str,'(2X,A)') "- Per atom stress (eV) used for heat flux calculation:"
    call localorb_info(tmp_str, use_unit, '(A)')

    if (use_pimd_wrapper .and. ipi_atomic_stress) then
      comm_string= trim(comm_string) // new_line('A') // trim(tmp_str) // new_line('A')
    end if

    write(tmp_str,'(2X,A)') "    Atom   | Stress components (1,1), (2,2), (3,3), (1,2), (1,3), (2,3)"
    call localorb_info(tmp_str, use_unit, '(A)')

    if (use_pimd_wrapper .and. ipi_atomic_stress) then
      comm_string= trim(comm_string) // trim(tmp_str) // new_line('A')
    end if

    write(tmp_str,'(2X,A)') "  -------------------------------------------------------------------"
    call localorb_info(tmp_str, use_unit, '(A)')
    if (use_pimd_wrapper .and. ipi_atomic_stress) then
      comm_string= trim(comm_string) // trim(tmp_str) // new_line('A')
    end if

    !Output:
    do i_atom=1,n_atoms   
      write(tmp_str,'(4X,I8,A,6E20.10)') i_atom,' | ',  &
                                        ( HF_stress_per_atom(1 ,1 ,i_atom) * Hartree  ), &
                                        ( HF_stress_per_atom(2 ,2 ,i_atom) * Hartree  ), &
                                        ( HF_stress_per_atom(3 ,3 ,i_atom) * Hartree  ), &
                                        ( HF_stress_per_atom(1 ,2 ,i_atom) * Hartree  ), &
                                        ( HF_stress_per_atom(1 ,3 ,i_atom) * Hartree  ), &
                                        ( HF_stress_per_atom(2 ,3 ,i_atom) * Hartree  ) 
      call localorb_info(tmp_str, use_unit, '(A)')

      if (use_pimd_wrapper .and. ipi_atomic_stress) then
        comm_string= trim(comm_string) // trim(tmp_str) // new_line('A')
      end if
    end do
    write(tmp_str,'(2X,A)') "  -------------------------------------------------------------------"
    call localorb_info(tmp_str, use_unit, '(A)')

    if (use_pimd_wrapper .and. ipi_atomic_stress) then
      comm_string= trim(comm_string) // trim(tmp_str) // new_line('A')
    end if

  end subroutine


  ! print the heat flux to stdout
  subroutine print_heat_flux()
    use dimensions
    use constants
    use localorb_io, only: localorb_info, use_unit
    IMPLICIT NONE
    integer                             :: i_coord
    character*400                       :: tmp_str

    write(tmp_str,'(4X,A,3E20.10)') "   Cartesian Components of the Conductive Heat Flux (eV*AA/ps): ", &
                                                           ( J(i_coord) * Hartree *bohr, i_coord=1,3,1 )
    call localorb_info(tmp_str, use_unit, '(A)')
  end subroutine


  ! Output atomic stress and heat flux to stdout
  subroutine print_atomic_stress_and_heat_flux()
    IMPLICIT NONE

    call print_atomic_stress()
    call print_heat_flux()

  end subroutine

  ! CC: Legacy debug routine, now replaced by
  !     assemble_heat_flux and print_heat_flux
  subroutine compute_and_print_heat_flux_old(HF_coords,velocities, Epot, Ekin)
   use analytical_stress
   use dimensions
   use physics
   use runtime_choices
   use synchronize_mpi
   use energy_density
   use molecular_dynamics
   use geometry
   use localorb_io, only: localorb_info
   IMPLICIT NONE
   real*8,dimension(1:3,1:n_atoms),intent(in) :: HF_coords 
   real*8,dimension(1:3,1:n_atoms),intent(in) :: velocities
   integer                           :: i_atom,i_coord,i_comp
   integer                           :: l_index(9),m_index(9),map_pu(9)
   character*400                     :: info_str
   real*8,dimension(1:3,1:3)           :: sum_pu,sum_MP_AT,sum_MP_EL,sum_MP_CO,stress
   ! real*8,dimension(1:3,1:3,1:n_atoms) :: HF_stress_per_atom
   real*8,dimension(1:3,1:n_atoms)     :: J_per_atom_1,     J_per_atom_2
   real*8,dimension(1:3,1:n_atoms)     :: J_per_atom_1_symm,J_per_atom_2_symm
   real*8,dimension(1:3)               :: J_1,J_2,J_1_symm,J_2_symm
   real*8,dimension(1:3,1:n_atoms)     :: J_conv_per_atom_0, J_conv_per_atom_t
   real*8,dimension(1:3)               :: J_conv_0, J_conv_t
   real*8,dimension(1:3,1:n_atoms)     :: HF_moment_per_atom_0, HF_moment_per_atom_t
   real*8,dimension(1:3)               :: HF_moment_0, HF_moment_t
   real*8,dimension(1:3,1:n_atoms)     :: HF_moment_per_atom_pbc_0, HF_moment_per_atom_pbc_t
   real*8,dimension(1:3)               :: HF_moment_pbc_0, HF_moment_pbc_t
   real*8                              :: HF_Ekin, HF_totEkin, HF_totE
   real*8                              :: Epot,Ekin
   real*8,dimension(1:3)               :: PBC_shift
   real*8,dimension(1:3,1:n_atoms)     :: HF_coords_pbc

   ! Map lm-index 
   !  1 4 5
   !  7 2 6
   !  8 9 3
   l_index    = (/ 1, 2, 3, 1, 1, 2, 2, 3, 3 /)
   m_index    = (/ 1, 2, 3, 2, 3, 3, 1, 1, 2 /)
   ! Map hessian-index 
   !  1 2 3
   !  7 4 5
   !  8 9 6
   map_pu     = (/ 1, 4, 6, 2, 3, 5, 7, 8, 9 /)

   sum_pu(:,:)=0.0d0
   sum_MP_AT(:,:)=0.0d0
   sum_MP_EL(:,:)=0.0d0
   sum_MP_CO(:,:)=0.0d0

   ! (A) Compute stress per atom from individual stress components
   if (.not.allocated(HF_stress_per_atom)) allocate(HF_stress_per_atom(1:3,1:3,1:n_atoms))
   HF_stress_per_atom(:,:,:) = 0.0d0
   do i_coord=1,AS_components,1
     HF_stress_per_atom(l_index(i_coord),m_index(i_coord),1:n_atoms) = &
       HF_stress_per_atom_PU(map_pu(i_coord),1:n_atoms)                      + &
       HF_stress_per_atom_MP_AT(l_index(i_coord),m_index(i_coord),1:n_atoms) + &
       HF_stress_per_atom_MP_EL(l_index(i_coord),m_index(i_coord),1:n_atoms) + &
       HF_stress_per_atom_MP_CO(l_index(i_coord),m_index(i_coord),1:n_atoms) 
       stress(l_index(i_coord),m_index(i_coord)) = sum(HF_stress_per_atom(l_index(i_coord),m_index(i_coord),1:n_atoms))
       ! Check for concistency:
       if ( abs(stress(l_index(i_coord),m_index(i_coord))/cell_volume - analytical_stress_tensor(l_index(i_coord),m_index(i_coord))) .gt. 1d-14 ) then
         write(info_str,'(2X,A,2I2,4E20.10)') " *** ERROR HF: Stress inconsistent! Diff: ", &
           l_index(i_coord),m_index(i_coord), (stress(l_index(i_coord),m_index(i_coord))/cell_volume - analytical_stress_tensor(l_index(i_coord),m_index(i_coord)))*Hartree/bohr**3.0d0, &
           stress(l_index(i_coord),m_index(i_coord))/cell_volume*Hartree/bohr**3.0d0, &
           stress(m_index(i_coord),l_index(i_coord))/cell_volume*Hartree/bohr**3.0d0, &
           analytical_stress_tensor(l_index(i_coord),m_index(i_coord))*Hartree/bohr**3.0d0
         call localorb_info(info_str,use_unit,'(2X,A)')
       end if
   end do

   do i_atom=1,n_atoms   
     write(info_str,'(2X,A,I4,6E20.10)') " PU   : ", i_atom, HF_stress_per_atom_PU( map_pu(1),i_atom),  HF_stress_per_atom_PU( map_pu(2),i_atom),  &
                                                             HF_stress_per_atom_PU( map_pu(3),i_atom),  HF_stress_per_atom_PU( map_pu(4),i_atom),  &
                                                             HF_stress_per_atom_PU( map_pu(5),i_atom),  HF_stress_per_atom_PU( map_pu(6),i_atom)
     call localorb_info(info_str,use_unit,'(2X,A)')
   end do
   do i_atom=1,n_atoms   
     write(info_str,'(2X,A,I4,6E20.10)') " MP_AT: ", i_atom, HF_stress_per_atom_MP_AT( l_index(1),m_index(1),i_atom),  HF_stress_per_atom_MP_AT( l_index(2),m_index(2),i_atom),  &
                                                             HF_stress_per_atom_MP_AT( l_index(3),m_index(3),i_atom),  HF_stress_per_atom_MP_AT( l_index(4),m_index(4),i_atom),  &
                                                             HF_stress_per_atom_MP_AT( l_index(5),m_index(5),i_atom),  HF_stress_per_atom_MP_AT( l_index(6),m_index(6),i_atom)
     call localorb_info(info_str,use_unit,'(2X,A)')
   end do
   do i_atom=1,n_atoms   
     write(info_str,'(2X,A,I4,6E20.10)') " MP_EL: ", i_atom, HF_stress_per_atom_MP_EL( l_index(1),m_index(1),i_atom),  HF_stress_per_atom_MP_EL( l_index(2),m_index(2),i_atom),  &
                                                             HF_stress_per_atom_MP_EL( l_index(3),m_index(3),i_atom),  HF_stress_per_atom_MP_EL( l_index(4),m_index(4),i_atom),  &
                                                             HF_stress_per_atom_MP_EL( l_index(5),m_index(5),i_atom),  HF_stress_per_atom_MP_EL( l_index(6),m_index(6),i_atom)
     call localorb_info(info_str,use_unit,'(2X,A)')
   end do
   write(info_str,'(2X,A,6E20.10)') " SUM PU   : ", sum(HF_stress_per_atom_PU( map_pu(1),1:n_atoms)),  sum(HF_stress_per_atom_PU( map_pu(2),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_PU( map_pu(3),1:n_atoms)),  sum(HF_stress_per_atom_PU( map_pu(4),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_PU( map_pu(5),1:n_atoms)),  sum(HF_stress_per_atom_PU( map_pu(6),1:n_atoms))
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,6E20.10)') " SUM MP_AT: ", sum(HF_stress_per_atom_MP_AT( l_index(1),m_index(1),1:n_atoms)),  sum(HF_stress_per_atom_MP_AT( l_index(2),m_index(2),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_MP_AT( l_index(3),m_index(3),1:n_atoms)),  sum(HF_stress_per_atom_MP_AT( l_index(4),m_index(4),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_MP_AT( l_index(5),m_index(5),1:n_atoms)),  sum(HF_stress_per_atom_MP_AT( l_index(6),m_index(6),1:n_atoms))
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,6E20.10)') " SUM MP_EL: ", sum(HF_stress_per_atom_MP_EL( l_index(1),m_index(1),1:n_atoms)),  sum(HF_stress_per_atom_MP_EL( l_index(2),m_index(2),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_MP_EL( l_index(3),m_index(3),1:n_atoms)),  sum(HF_stress_per_atom_MP_EL( l_index(4),m_index(4),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_MP_EL( l_index(5),m_index(5),1:n_atoms)),  sum(HF_stress_per_atom_MP_EL( l_index(6),m_index(6),1:n_atoms))
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,6E20.10)') " SUM MP_CO: ", sum(HF_stress_per_atom_MP_CO( l_index(1),m_index(1),1:n_atoms)),  sum(HF_stress_per_atom_MP_CO( l_index(2),m_index(2),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_MP_CO( l_index(3),m_index(3),1:n_atoms)),  sum(HF_stress_per_atom_MP_CO( l_index(4),m_index(4),1:n_atoms)),  &
                                                    sum(HF_stress_per_atom_MP_CO( l_index(5),m_index(5),1:n_atoms)),  sum(HF_stress_per_atom_MP_CO( l_index(6),m_index(6),1:n_atoms))
   call localorb_info(info_str,use_unit,'(2X,A)')

   ! Symmetrize ...
   if (AS_components .ne. 9) then
     HF_stress_per_atom(2,1,1:n_atoms) = HF_stress_per_atom(1,2,1:n_atoms)
     HF_stress_per_atom(3,1,1:n_atoms) = HF_stress_per_atom(1,3,1:n_atoms)
     HF_stress_per_atom(3,2,1:n_atoms) = HF_stress_per_atom(2,3,1:n_atoms)
   end if

   ! (B) Get diffusive heat flux from stress*v
   J_per_atom_1(:,:) = 0.0d0
   J_per_atom_2(:,:) = 0.0d0
   J_per_atom_1_symm(:,:) = 0.0d0
   J_per_atom_2_symm(:,:) = 0.0d0
   ! Fancy matmul preferable, but not critical...
   do i_comp=1,3,1
     do i_coord=1,3
       do i_atom=1,n_atoms   
         J_per_atom_1(i_comp,i_atom) = J_per_atom_1(i_comp,i_atom) + ( HF_stress_per_atom(i_comp ,i_coord,i_atom)*velocities(i_coord,i_atom) )
         J_per_atom_2(i_comp,i_atom) = J_per_atom_2(i_comp,i_atom) + ( HF_stress_per_atom(i_coord,i_comp, i_atom)*velocities(i_coord,i_atom) )
         !! write(info_str,'(2X,3I4,1E20.10)') i_comp,i_coord,i_atom,HF_stress_per_atom(i_comp ,i_coord,i_atom)
         !! call localorb_info(info_str,use_unit,'(2X,A)')
         if ( i_comp .le. i_coord) then
           J_per_atom_1_symm(i_comp,i_atom) = J_per_atom_1_symm(i_comp,i_atom) + ( HF_stress_per_atom(i_comp ,i_coord,i_atom)*velocities(i_coord,i_atom) )
           J_per_atom_2_symm(i_comp,i_atom) = J_per_atom_2_symm(i_comp,i_atom) + ( HF_stress_per_atom(i_coord,i_comp, i_atom)*velocities(i_coord,i_atom) )
         else
           J_per_atom_1_symm(i_comp,i_atom) = J_per_atom_1_symm(i_comp,i_atom) + ( HF_stress_per_atom(i_coord,i_comp ,i_atom)*velocities(i_coord,i_atom) )
           J_per_atom_2_symm(i_comp,i_atom) = J_per_atom_2_symm(i_comp,i_atom) + ( HF_stress_per_atom(i_comp, i_coord,i_atom)*velocities(i_coord,i_atom) )
         end if
       end do
     end do
     J_1(i_comp)      = sum(J_per_atom_1(i_comp,1:n_atoms))
     J_2(i_comp)      = sum(J_per_atom_2(i_comp,1:n_atoms))
     J_1_symm(i_comp) = sum(J_per_atom_1_symm(i_comp,1:n_atoms))
     J_2_symm(i_comp) = sum(J_per_atom_2_symm(i_comp,1:n_atoms))
   end do

   !! ! (C) Get convective heat flux 
   !! ! C.1 Total MD energy 
   !! call calculate_kinetic_energy(velocities,HF_totEkin)
   !! HF_totE =sum(HF_energy_per_atom(1:n_atoms)) + HF_totEkin

   !! ! C.2 Get average/eq values and respective differences 
   !! if (.not.allocated(HF_Eq_energy_per_atom)) then
   !!   allocate(HF_Eq_energy_per_atom(1:n_atoms))
   !!   !! NOT NECESSARY call HF_read_energy_per_atom()
   !!   HF_Eq_energy_per_atom = 0.0d0
   !!   HF_sum_Eq_energy = sum(HF_Eq_energy_per_atom(1:n_atoms))
   !!   ! For averaging
   !!   allocate(HF_Avg_delta_energy_per_atom_0(1:n_atoms))
   !!   allocate(HF_Avg_delta_energy_per_atom_t(1:n_atoms))
   !!   HF_Avg_delta_energy_per_atom_0(:) = 0.0d0
   !!   HF_Avg_delta_energy_per_atom_t(:) = 0.0d0
   !!   HF_step = 0
   !!   HF_Delta_totE_Eq_E_0 = (HF_totE - HF_sum_Eq_energy) / dble(n_atoms)
   !! end if
   !! HF_Delta_totE_Eq_E_t = (HF_totE - HF_sum_Eq_energy) / dble(n_atoms)


   !! ! Generate PBC cleaned coords for HF Moment
   !! if (.not.allocated(HF_old_coords)) then
   !!   allocate(HF_old_coords(1:3,1:n_atoms))
   !!   HF_old_coords = HF_coords
   !!   do i_atom=1,n_atoms
   !!     call map_to_center_cell(HF_old_coords(:,i_atom))
   !!   end do
   !! end if
   !! !!Make motion continuous in PBC
   !! do i_atom=1,n_atoms
   !!   !Difference wrt to old positions
   !!   PBC_shift(1:3) = HF_coords(1:3,i_atom) -  HF_old_coords(1:3,i_atom)
   !!   ! Map to center cell gives PBC cleaned distance
   !!   call map_to_center_cell(PBC_shift(:))
   !!   ! Shift to continuous motion 
   !!   HF_coords_pbc(1:3,i_atom) = HF_old_coords(1:3,i_atom) + PBC_shift(1:3)
   !!   ! Remember old unfolded position for next run
   !!   HF_old_coords(1:3,i_atom) = HF_coords_pbc(1:3,i_atom)
   !! end do
   !! 
   !! ! C.3. Convective heat flux / Helfand Moment
   !! HF_step = HF_step + 1
   !! do i_atom=1,n_atoms
   !!   ! get kinetic_energy
   !!   call calculate_kinetic_energy(velocities,HF_Ekin,i_atom,i_atom)
   !!   J_conv_per_atom_0(1:3,i_atom)        = ( HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_0 ) * velocities(1:3,i_atom)
   !!   J_conv_per_atom_t(1:3,i_atom)        = ( HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_t ) * velocities(1:3,i_atom)
   !!   HF_moment_per_atom_0(1:3,i_atom)     = ( HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_0 ) * HF_coords(1:3,i_atom)
   !!   HF_moment_per_atom_t(1:3,i_atom)     = ( HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_t ) * HF_coords(1:3,i_atom)
   !!   HF_moment_per_atom_pbc_0(1:3,i_atom) = ( HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_0 ) * HF_coords_pbc(1:3,i_atom)
   !!   HF_moment_per_atom_pbc_t(1:3,i_atom) = ( HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_t ) * HF_coords_pbc(1:3,i_atom)
   !!   ! Avg for debug
   !!   HF_Avg_delta_energy_per_atom_0(i_atom) = HF_Avg_delta_energy_per_atom_0(i_atom) + (HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_0 )
   !!   HF_Avg_delta_energy_per_atom_t(i_atom) = HF_Avg_delta_energy_per_atom_t(i_atom) + (HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_t )
   !!   write(info_str,'(A,I4,4EN30.20)') "*** E STATS per atom",i_atom, &
   !!     (HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_0 )*Hartree, &
   !!     (HF_Ekin + HF_energy_per_atom(i_atom) - HF_Eq_energy_per_atom(i_atom) - HF_Delta_totE_Eq_E_t )*Hartree, &
   !!     HF_Avg_delta_energy_per_atom_0(i_atom)*Hartree/dble(HF_step), &
   !!     HF_Avg_delta_energy_per_atom_t(i_atom)*Hartree/dble(HF_step)
   !!   call localorb_info(info_str,use_unit,'(2X,A)')
   !! end do
   !! ! Full stat
   !! write(info_str,'(A,4EN30.20)') "*** E STATS ALL atoms", &
   !!   (HF_totE - HF_sum_Eq_energy - ( HF_Delta_totE_Eq_E_0 *dble(n_atoms)) )*Hartree, &
   !!   (HF_totE - HF_sum_Eq_energy - ( HF_Delta_totE_Eq_E_t *dble(n_atoms)) )*Hartree, &
   !!   sum(HF_Avg_delta_energy_per_atom_0(1:n_atoms))/dble(n_atoms)*Hartree/dble(HF_step), &
   !!   sum(HF_Avg_delta_energy_per_atom_t(1:n_atoms))/dble(n_atoms)*Hartree/dble(HF_step) 
   !! call localorb_info(info_str,use_unit,'(2X,A)')
   !! ! Sum everything up
   !! do i_coord=1,3
   !!   J_conv_0(i_coord)        = sum( J_conv_per_atom_0(i_coord,1:n_atoms) ) 
   !!   J_conv_t(i_coord)        = sum( J_conv_per_atom_t(i_coord,1:n_atoms) ) 
   !!   HF_moment_0(i_coord)     = sum( HF_moment_per_atom_0(i_coord,1:n_atoms) )
   !!   HF_moment_t(i_coord)     = sum( HF_moment_per_atom_t(i_coord,1:n_atoms) )
   !!   HF_moment_pbc_0(i_coord) = sum( HF_moment_per_atom_pbc_0(i_coord,1:n_atoms) )
   !!   HF_moment_pbc_t(i_coord) = sum( HF_moment_per_atom_pbc_t(i_coord,1:n_atoms) )
   !! end do
   !! ! Safety check:
   !! if ( abs(HF_totE-Ekin-Epot) .gt. 1d-12 ) then
   !!   write(info_str,'(A,7EN30.20)') "*** ENERGY ERROR",(sum(HF_energy_per_atom(1:n_atoms))-Epot)*Hartree,(HF_totEkin-Ekin)*hartree, (HF_totE-Epot-Ekin)*hartree
   !!   call localorb_info(info_str,use_unit,'(2X,A)')
   !!   !! call aims_stop
   !! end if


   !Output:
   do i_atom=1,n_atoms   
     write(info_str,'(2X,A,I8,3E20.10)') " J1  : ", i_atom,  ( J_per_atom_1(i_coord,i_atom) * Hartree * bohr, i_coord=1,3,1 )
     call localorb_info(info_str,use_unit,'(2X,A)')
     write(info_str,'(2X,A,I8,3E20.10)') " J2  : ", i_atom,  ( J_per_atom_2(i_coord,i_atom) * Hartree * bohr, i_coord=1,3,1 )
     call localorb_info(info_str,use_unit,'(2X,A)')
     write(info_str,'(2X,A,I8,3E20.10)') " J1S : ", i_atom,  ( J_per_atom_1_symm(i_coord,i_atom) * Hartree * bohr, i_coord=1,3,1 )
     call localorb_info(info_str,use_unit,'(2X,A)')
     write(info_str,'(2X,A,I8,3E20.10)') " J2S : ", i_atom,  ( J_per_atom_2_symm(i_coord,i_atom) * Hartree * bohr, i_coord=1,3,1 )
     call localorb_info(info_str,use_unit,'(2X,A)')
     write(info_str,'(2X,A,I8,9E20.10)') " HF_stress_per_atom : ", i_atom,  ( HF_stress_per_atom(1 ,1 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(2 ,2 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(3 ,3 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(1 ,2 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(1 ,3 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(2 ,3 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(2 ,1 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(3 ,1 ,i_atom) * Hartree  ), &
                                                                            ( HF_stress_per_atom(3 ,2 ,i_atom) * Hartree  )
     call localorb_info(info_str,use_unit,'(2X,A)')
   end do
   write(info_str,'(2X,A,3E20.10)') " J_1_all  : ", ( J_1(i_coord) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_2_all  : ", ( J_2(i_coord) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_1_allS : ", ( J_1_symm(i_coord) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_2_allS : ", ( J_2_symm(i_coord) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_0    : ", ( J_conv_0(i_coord) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_t    : ", ( J_conv_t(i_coord) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   ! Sums
   write(info_str,'(2X,A,3E20.10)') " J_c_0 + J_1  : ", ( (J_conv_0(i_coord) + J_1(i_coord)      ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_0 + J_2  : ", ( (J_conv_0(i_coord) + J_2(i_coord)      ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_0 + J_1_S: ", ( (J_conv_0(i_coord) + J_1_symm(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_0 + J_2_S: ", ( (J_conv_0(i_coord) + J_2_symm(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_t + J_1  : ", ( (J_conv_t(i_coord) + J_1(i_coord)      ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_t + J_2  : ", ( (J_conv_t(i_coord) + J_2(i_coord)      ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_t + J_1_S: ", ( (J_conv_t(i_coord) + J_1_symm(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " J_c_t + J_2_S: ", ( (J_conv_t(i_coord) + J_2_symm(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " HfMom_0      : ", ( (HF_moment_0(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " HfMom_t      : ", ( (HF_moment_t(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " HfMom_pbc_0  : ", ( (HF_moment_pbc_0(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')
   write(info_str,'(2X,A,3E20.10)') " HfMom_pbc_t  : ", ( (HF_moment_pbc_t(i_coord) ) * Hartree *bohr, i_coord=1,3,1 )
   call localorb_info(info_str,use_unit,'(2X,A)')

   ! Deallocate to be ready for next step
   if(allocated( HF_rho_multipole_per_atom )) deallocate( HF_rho_multipole_per_atom ) 
   if(allocated( HF_rho_free_per_atom      )) deallocate( HF_rho_free_per_atom      )
   if(allocated( HF_dde_rho                )) deallocate( HF_dde_rho                )

  end subroutine compute_and_print_heat_flux_old



END MODULE heat_flux
