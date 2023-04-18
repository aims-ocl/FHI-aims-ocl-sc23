!****s* FHI-aims/update_hartree_potential_shanghui_p1
!  NAME
!    update_hartree_potential_shanghui_p1
!  SYNOPSIS

subroutine update_hartree_potential_shanghui_p1 &
     ( i_center_input, partition_tab, delta_rho)

!  PURPOSE
!  Input delta_rho(1:n_centers), get rho_multipole_supercell in the hartree_potential_storage moldue.
!  Which can be used in the next subrution sum_up_whole_potential.f90 to get
!  delta_potential. 
!
!  
!  shanghui 2015.01
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
  use hartree_potential_storage , only : rho_multipole_supercell 
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer, intent(in)    :: i_center_input
  real*8, dimension(n_full_points),intent(in)   :: partition_tab  
  real*8, dimension(n_centers_in_sc_DFPT,n_full_points), intent(in) :: delta_rho

!  INPUTS
!   o  partition_tab -- values of partition function
!   o  delta_rho -- electron density
!
!  OUTPUT
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
!  real*8 multipole_radius
  real*8 ylm_tab((l_pot_max+1)**2)
  real*8, allocatable :: current_rho_multipole(:,:) 

  integer l_h_dim(n_atoms)
  real*8 :: drho

   


!  counters

  integer i_cell_input,i_cell_multipole,i_cell_delta_rho
  integer i_cell_multipole_in_hamiltonian
  integer i_center_multipole, i_center_delta_rho
  integer i_center_multipole_in_hamiltonian 
  integer i_atom_input, i_atom_multipole
  integer i_index
  integer current_atom, current_radial, current_angular
  integer i_full_points
  integer i_spin
  integer i_my_batch


  character*100 :: info_str



!  begin work

  write(info_str,'(2X,A)') " "
  call localorb_info(info_str, use_unit, '(A)', OL_norm )
  write(info_str,'(2X,A,A)') "Evaluating partitioned Hartree potential by multipole expansion."
  call localorb_info(info_str, use_unit, '(A)', OL_norm )


  allocate(current_rho_multipole((l_pot_max+1)**2, n_max_radial+2),stat=i_index)
  call check_allocation(i_index, 'current_rho_multipole')



  do i_atom_multipole = 1, n_atoms, 1
    l_h_dim(i_atom_multipole) = (l_hartree(species(i_atom_multipole))+1)**2
  end do

!--------------------------(1) near field rho_multopole------------
  rho_multipole_supercell = 0.0d0
   i_cell_input = center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center_input) 
   i_atom_input = center_in_sc_DFPT_to_atom(i_center_input)  

 do i_cell_multipole = 1, n_cells_in_sc_DFPT
    i_cell_delta_rho = cell_diff_sc_DFPT(i_cell_input,i_cell_multipole)

  !========== i_center_delta_rho : index used in delta_rho, refer to i_center_input under translation invariable
    i_center_delta_rho = cell_and_atom_to_center_sc_DFPT(i_cell_delta_rho, i_atom_input)  
  !====================================================================================   
  
  i_full_points = 0
  do i_my_batch = 1, n_my_batches, 1
     
    do i_index = 1, batches(i_my_batch)%size, 1
           
      current_atom = batches(i_my_batch) % points(i_index) % index_atom
  !========== i_center_multipole : index used in rho_multipole_supercell=========================   
      i_center_multipole = cell_and_atom_to_center_sc_DFPT(i_cell_multipole,current_atom)
  !====================================================================================   
      current_radial = batches(i_my_batch) % points(i_index) % index_radial
      current_angular = batches(i_my_batch) % points(i_index) % index_angular
 
      i_full_points = i_full_points + 1
 
      ! execute only if partition_tab.gt.0 here, i.e. if the integration point
      ! makes sense
      !!!!! This is NOT the usual integration partition tab
      !!!!! This is the HARTREE partition tab. 
      !!!!! Thus this criterion is not the same as in all other integrations.
      !!!!! We must fix that one day.
      if (partition_tab(i_full_points).gt.0.d0) then
 
 
        ! compute atom-centered coordinates of current integration point,
        ! BUT HERE ONLY FOR PRESENT ATOM!
        dir_tab(:) = r_angular( : , current_angular, current_radial, species(current_atom))
 
        ylm_tab (1:l_h_dim(current_atom)) = &
             local_ylm_tab(1:l_h_dim(current_atom),current_angular, &
             lebedev_grid_index(current_radial,species(current_atom)))
            
        ! implied loop over all (l,m) from (0,0) to (l_hartree,l_hartree)
        rho_multipole_supercell(1:l_h_dim(current_atom), current_radial+1, i_center_multipole) = &
        rho_multipole_supercell(1:l_h_dim(current_atom), current_radial+1, i_center_multipole) + &
             ylm_tab(1:l_h_dim(current_atom)) * partition_tab(i_full_points) * & 
                               delta_rho(i_center_delta_rho,i_full_points) 
 
      end if
    end do !    end loop over a batch
  end do !      end loop over batches
 end do   !      end loop over i_cell_multipole
 
  ! Add all rho_multipole contributions
 

  do i_center_multipole = 1, n_centers_in_sc_DFPT 

    current_rho_multipole(:,:) = rho_multipole_supercell(:,:,i_center_multipole)
    call sync_vector(current_rho_multipole, ((l_pot_max+1)**2)*(n_max_radial+2))
    rho_multipole_supercell(:,:,i_center_multipole) = current_rho_multipole(:,:)
 
  enddo
 
  ! Set boundaries on rho_multipole
 
  
  do i_center_multipole = 1, n_centers_in_sc_DFPT 

     i_atom_multipole = center_in_sc_DFPT_to_atom(i_center_multipole) 
    ! At Infinity (n_radial+2)
    ! enforce zero charge at infinity explicitly by a trick:
    ! (used in splines later on)
    rho_multipole_supercell(1:l_h_dim(i_atom_multipole), n_radial(species(i_atom_multipole))+2, & 
                            i_center_multipole) = 0.0d0
 
    ! At zero
    drho = (rho_multipole_supercell(1,2,i_center_multipole)  & 
           -rho_multipole_supercell(1,3,i_center_multipole)) &
    &     /    (r_radial(1,species(i_atom_multipole)) - r_radial(2,species(i_atom_multipole)))
 
    if (legacy_monopole_extrapolation) then
       ! Backwards compatible but inaccurate choice of sign.
       ! Please note that the influence of this boundary rapidly decreases
       ! with increasing grid density.
       rho_multipole_supercell(1,1,i_center_multipole) = rho_multipole_supercell(1,2,i_center_multipole) &
       &                                 + drho * r_radial(1,species(i_atom_multipole))
    else
       rho_multipole_supercell(1,1,i_center_multipole) = rho_multipole_supercell(1,2,i_center_multipole) &
       &                                 - drho * r_radial(1,species(i_atom_multipole))
    end if
    rho_multipole_supercell(2:l_h_dim(i_atom_multipole),1,i_center_multipole) = 0.0d0
  enddo
 


  deallocate(current_rho_multipole)


end subroutine  update_hartree_potential_shanghui_p1
!******
