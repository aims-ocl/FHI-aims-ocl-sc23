module KH_core_states


    use dimensions
    use runtime_choices
    use grids
    use geometry
    use basis
    use mpi_utilities
    use synchronize_mpi
    use localorb_io
    use constants
    use species_data
    use pbc_lists
    use spline, only: val_spline

    implicit none

    real*8, allocatable, private, dimension(:,:,:,:) :: Core_KH_eigenvectors
    real*8, allocatable, private, dimension(:,:,:)   :: Core_KH_eigenvalues
    real*8, allocatable, private, dimension(:,:,:)   :: Core_free_eigenvalues
 
    integer, allocatable, private, dimension(:) :: basis_to_core_basis

    integer, private:: n_max_KH_core_basis_species

    integer,private:: n_total_states
    integer,private:: n_KH_valence_states

    integer,private:: n_max_KH_core_states_species

    real*8, private, dimension(:,:), allocatable :: v_psi

    integer:: n_KH_core_states





contains

!--------------------------------------------------------------
 subroutine initialize_KH_core_states

    use lapack_wrapper
    integer:: atoms_in_species(n_species)
    integer:: i_atom, i_species, i_basis



    if(.not. flag_KH_core_states)   return

    n_total_states = n_states


    n_KH_core_states = 0
    do i_atom = 1, n_atoms
       n_KH_core_states = n_KH_core_states + n_KH_core_states_species(species(i_atom))
    end do

    if(n_KH_core_states == 0)then
       flag_KH_core_states = .false.
       return
    end if


    if(n_KH_core_states > 0)then
       
       n_KH_valence_states = n_total_states -  n_KH_core_states

       n_KH_core_basis_species = 0
       atoms_in_species = 0

       do i_basis = 1, n_basis

          n_KH_core_basis_species(species(Cbasis_to_atom(i_basis))) = n_KH_core_basis_species(species(Cbasis_to_atom(i_basis))) + 1

       end do
 
       do i_atom = 1, n_atoms
          atoms_in_species (species(i_atom))    =  atoms_in_species(species(i_atom)) + 1
       end do

      
       do i_species = 1, n_species
          if(atoms_in_species(i_species) >0)then
             n_KH_core_basis_species(i_species) = n_KH_core_basis_species(i_species)/ atoms_in_species(i_species)
          end if
       end do



       n_max_KH_core_states_species = 0
       n_max_KH_core_basis_species = 0
       do i_species = 1, n_species
          n_max_KH_core_states_species = max(n_max_KH_core_states_species, n_KH_core_states_species(i_species))
          n_max_KH_core_basis_species  = max(n_max_KH_core_basis_species,  n_KH_core_basis_species(i_species))
       end do
       
       allocate(Core_KH_eigenvectors(n_max_KH_core_basis_species, n_max_KH_core_states_species, n_spin,n_atoms))
       allocate(Core_KH_eigenvalues( n_max_KH_core_states_species, n_spin, n_atoms))
       allocate(Core_free_eigenvalues( n_max_KH_core_states_species, n_spin, n_atoms))
       
       Core_KH_eigenvectors = 0.d0
       Core_KH_eigenvalues = 0.d0
       

       
       if(myid == 0)then
          
          write(use_unit,*) 
          write(use_unit,'(A)') '  KH core states information:'

          write(use_unit,'(A,I5)') '  | Total number of states:                ', n_total_states
          write(use_unit,'(A,I5)') '  | Total number of KH core states:        ', n_KH_core_states
          write(use_unit,'(A,I5)') '  | Total number of valence states:        ', n_KH_valence_states
          write(use_unit,'(A,I5)') '  | Max KH core states in species:         ', n_max_KH_core_states_species

          do i_species = 1, n_species
             write(use_unit,'(A,A,A,I5)')    '  | KH core states of ',trim(species_name(i_species)),':                  ', &
                  n_KH_core_states_species(i_species)
          end do

          write(use_unit,'(A,I5)') '  | Total number of basis functions:       ', n_basis 
          write(use_unit,'(A,I5)') '  | Max KH core basis functions in species:',n_max_KH_core_basis_species

          do i_species = 1, n_species
             write(use_unit,'(A,A,A,I5)')    '  | KH core basis of ',trim(species_name(i_species)),':                   ', &
                  n_KH_core_basis_species(i_species)
          end do

          write(use_unit,*) 
       end if

    end if

    call  allocate_work_space(max(n_basis,n_max_KH_core_basis_species))
    allocate( v_psi(n_hamiltonian_matrix_size, n_spin))
    allocate(basis_to_core_basis(n_basis))


    

  end subroutine initialize_KH_core_states
!--------------------------------------------------------------

  subroutine deallocate_KH_core_states

    if(allocated( Core_KH_eigenvectors ))    deallocate(Core_KH_eigenvectors)
    if(allocated( Core_KH_eigenvalues  ))    deallocate(Core_KH_eigenvalues)
    if(allocated( Core_free_eigenvalues))    deallocate(Core_free_eigenvalues)
    if(allocated( v_psi                ))    deallocate( v_psi )
    if(allocated( basis_to_core_basis  ))    deallocate(basis_to_core_basis)


  end subroutine deallocate_KH_core_states
!--------------------------------------------------------------


  subroutine evaluate_KH_core_states &
       (   overlap_matrix, rho, rho_gradient, partition_tab, hartree_potential, &
           kinetic_density )

    use pbc_lists
    use lapack_wrapper
    use free_atoms


    real*8:: V_monopole_log_spl(4, n_max_grid, n_species )
    real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
    real*8, dimension(n_full_points) ::            partition_tab
    real*8, dimension(n_spin, n_full_points) ::    rho
    real*8, dimension(n_full_points) ::            hartree_potential
    real*8, dimension(n_spin, n_full_points) ::    kinetic_density

    real*8 overlap_matrix( n_hamiltonian_matrix_size)

    integer,dimension(:),allocatable:: core_basis, core_basis_inv, loc_min

    integer:: n_core_b

    real*8, dimension(:),allocatable:: core_ovl
    real*8, dimension(:,:),allocatable:: core_ham, relativistic_terms

    integer:: i_basis_1, i_basis_2, i_b_1, i_b_2,i_b_3, i_atom, i_size, i_index, i_state

    integer, dimension(:,:),allocatable:: b_index
    real*8, dimension(:),allocatable:: core_ovl_old
    real*8, dimension(:,:),allocatable:: core_ham_old
    real*8:: eigen_old
    integer:: i_grid, i_l
    integer:: temp(1000)

    real*8,allocatable,dimension(:)::eigenvalues,eigenvalues_free, eigenvalues_core

    integer,allocatable,dimension(:):: eigenvalues_place, eigenvalues_place_temp
        

    integer:: n_KH_core_basis_prune, i_term_max
    real*8:: time1, time2, eigen, norm
    integer:: i_cell, i_eigen, i_spin


      real*8:: eigenvectors(n_max_KH_core_basis_species, n_max_KH_core_states_species)

      

    if(.not. flag_KH_core_states)   return
    if(myid==0) write(use_unit,*) 'Evaluate KH core states'




    call cpu_time(time1)


    allocate(core_ovl(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2))
    allocate(core_ham(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2,n_spin))
    allocate(b_index(n_max_KH_core_basis_species, n_max_KH_core_basis_species))

    allocate(core_basis(n_max_KH_core_basis_species))
    allocate(core_basis_inv(n_basis))
    allocate(loc_min(n_max_KH_core_basis_species))

    allocate(relativistic_terms(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2,n_spin))

    allocate(core_ovl_old(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2))
    allocate(core_ham_old(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2,n_spin))

    call  allocate_work_space(max(n_basis,n_max_KH_core_basis_species))

    allocate(eigenvalues( n_max_KH_core_states_species))
    allocate(eigenvalues_free( n_max_KH_core_states_species))
    allocate(eigenvalues_core( n_max_KH_core_basis_species ))


    allocate(eigenvalues_place_temp( n_max_KH_core_basis_species))
    allocate(eigenvalues_place( n_max_KH_core_basis_species))


    call get_KH_free_atoms_potential( V_monopole_log_spl )

    call integrate_v_psi ( hartree_potential, rho, rho_gradient, kinetic_density, &
                           partition_tab, l_shell_max,  v_psi )


 
    basis_to_core_basis = 0    




    do i_atom = 1, n_atoms

       if(n_KH_core_states_species(species(i_atom))> 0)then
          
       
       eigenvalues_core = 0.d0
       core_basis_inv = 0.d0
       i_index = 0


       do i_basis_1 = 1, n_basis,1
          if( Cbasis_to_atom(i_basis_1) == i_atom)then
             
             i_index = i_index + 1
             core_basis(i_index) = i_basis_1
             core_basis_inv(i_basis_1) = i_index

             basis_to_core_basis(i_basis_1) = i_index

             
          end if
       end do
       
       n_core_b = i_index


       i_index = 0
       b_index = 0.d0
       do i_b_2 = 1, n_core_b,1
          do i_b_1 = 1, i_b_2,1
             
             i_basis_1 = core_basis(i_b_1)
             i_basis_2 = core_basis(i_b_2)
             
             i_index = i_index + 1
             b_index(i_b_1, i_b_2) =  i_index
             
          end do
       end do
       
       

       core_ovl = 0.d0
       core_ham = 0.d0


       select case(packed_matrix_format)

       case(PM_index)

          do i_cell = 1,n_cells_in_hamiltonian-1
             do i_b_2 = 1, n_core_b

             i_basis_2 = core_basis(i_b_2) 
             
             if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

                i_index = index_hamiltonian(1,i_cell, i_basis_2)-1
                
                do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
                   
                   i_index = i_index + 1
                   i_b_1 = core_basis_inv( column_index_hamiltonian(i_index))

                   if(i_b_1 > 0)then
                   

                      core_ham( b_index(i_b_1, i_b_2),1:n_spin) = core_ham( b_index(i_b_1, i_b_2),1:n_spin) & 
                                                                  + v_psi(i_index,1:n_spin)
                      core_ovl( b_index(i_b_1, i_b_2)) = core_ovl( b_index(i_b_1, i_b_2)) + overlap_matrix(i_index)

                   end if
                   
                end do
             end if
          end do
       end do

    case(PM_none)


       i_index = 0
       do i_basis_2 = 1, n_centers_basis_T
          do i_basis_1 = 1, i_basis_2

             i_index = i_index +1


             i_b_1 = core_basis_inv( Cbasis_to_basis( i_basis_1))
             i_b_2 = core_basis_inv(Cbasis_to_basis( i_basis_2))

             if(i_b_1 > 0 .and. i_b_2 > 0)then

                core_ham( b_index(i_b_1, i_b_2),1:n_spin) = core_ham( b_index(i_b_1, i_b_2),1:n_spin) + v_psi(i_index,1:n_spin)
                core_ovl( b_index(i_b_1, i_b_2)) = core_ovl( b_index(i_b_1, i_b_2)) + overlap_matrix(i_index)
                
             end if
          end do
       end do

    end select


    i_b_2 = 0

    do i_b_1 = 1, n_atomic(species(i_atom)), 1

       do i_l = -atomic_l(species(i_atom), i_b_1), atomic_l(species(i_atom), i_b_1),1

          i_b_2 = i_b_2 + 1
          eigenvalues_core(i_b_2) = free_wave_eigenval(species(i_atom),i_b_1 )
          eigenvalues_place_temp(i_b_2) = i_b_1

       end do
    end do
           



      do i_b_1 = 1, n_KH_core_states_species(species(i_atom))
         loc_min = minloc( eigenvalues_core)
          i_b_2 = loc_min(1)

          eigenvalues_free(i_b_1) =  eigenvalues_core(i_b_2)
          eigenvalues_place(i_b_1) = eigenvalues_place_temp(i_b_2)
          
          eigenvalues_core(i_b_2) = 100.d0
          
       end do


       core_ham_old = core_ham
       core_ovl_old = core_ovl



       do i_eigen = 1, n_KH_core_states_species(species(i_atom))


             do i_spin = 1, n_spin


                   call  integrate_KH_kinetic_term(core_basis,  n_KH_core_basis_species(species(i_atom)), & 
                         relativistic_terms(1,i_spin), i_atom, & 
                         free_wave_eigenval(species(i_atom),eigenvalues_place(i_eigen)),  & 
                         V_monopole_log_spl, b_index, core_ham_old )

                Core_free_eigenvalues( i_eigen, i_spin, i_atom) =  free_wave_eigenval(species(i_atom),eigenvalues_place(i_eigen))
                
                core_ham(:,i_spin)  =  core_ham_old(:,i_spin)  +  relativistic_terms(:,i_spin)

             end do


             eigen = Core_KH_eigenvalues( 1, 1, i_atom)

             do i_spin = 1, n_spin

                core_ovl = core_ovl_old

                call  real_lapack_solver (n_core_b, n_KH_core_states_species(species(i_atom)), & 
                     core_ovl, core_ham(1,i_spin), 1d-12, &
                     .false., eigenvalues(1),  eigenvectors(1:n_core_b, 1:n_KH_core_states_species(species(i_atom))))


                Core_KH_eigenvalues( i_eigen, i_spin, i_atom) = eigenvalues(i_eigen)
                Core_KH_eigenvectors(1:n_core_b, i_eigen, i_spin, i_atom)  =  eigenvectors(1:n_core_b, i_eigen)

             end do





          if(  use_small_component )then

             do i_spin = 1, n_spin
             

                call  normalize(core_basis,  n_KH_core_basis_species(species(i_atom)), relativistic_terms(1,i_spin), &
                     i_atom,  free_wave_eigenval(species(i_atom),eigenvalues_place(i_eigen)),  V_monopole_log_spl, b_index )
                

                norm  = 0.d0
             
                do i_b_2 = 1, n_core_b
                   do i_b_1 = 1, i_b_2       

                      norm = norm +  Core_KH_eigenvectors(i_b_1, i_eigen, i_spin,i_atom) & 
                           * Core_KH_eigenvectors(i_b_2, i_eigen, i_spin,i_atom) &
                           *  relativistic_terms( b_index(i_b_1, i_b_2),i_spin)

                   end do
                end do
                do i_b_2 = 1, n_core_b
                   do i_b_1 = 1, i_b_2-1
                   
                      norm = norm +  Core_KH_eigenvectors(i_b_1, i_eigen, i_spin,i_atom) & 
                           * Core_KH_eigenvectors(i_b_2, i_eigen, i_spin,i_atom) &
                           *  relativistic_terms( b_index(i_b_1, i_b_2),i_spin)
                   
                   end do
                end do


                do i_b_1 = 1, n_core_b
                   Core_KH_eigenvectors(i_b_1, i_eigen, i_spin,i_atom) = & 
                     Core_KH_eigenvectors(i_b_1, i_eigen, i_spin,i_atom)/sqrt(norm + 1)
                end do

             end do
          end if
       end do

    end if ! (n_KH_core_states_species(species(i_atom))> 0)
 end do ! i_atom


    deallocate(core_ovl)
    deallocate(core_ham)
    deallocate(b_index)
    deallocate(core_basis)
    deallocate(core_basis_inv)
    deallocate(relativistic_terms)
    deallocate(eigenvalues)
    deallocate(eigenvalues_free)




     end subroutine evaluate_KH_core_states

!-------------------------------------------------------------------


    subroutine integrate_KH_kinetic_term( basis_list, n_basis_list, integrals, i_atom, eigen, & 
                                          V_monopole_log_spl, b_index, core_ham_old )
      
      use spline
      use free_atoms

      integer:: i_basis_1, i_basis_2, i_function_1,i_function_2, i_atom, i_grid, i_basis_L, i_basis_L2, i_state, l_1,l_2, m_1,m_2
      real*8::  V_radial_deriv, eige, V_radial
      integer:: n_basis_list
      integer:: basis_list(n_basis_list)
      real*8:: eigen, eigen_old
      real*8:: alpha
      integer:: temp(n_basis_list)
      real*8:: V_monopole_log_spl(4, n_max_grid, n_species )
      integer:: b_index(n_max_KH_core_basis_species, n_max_KH_core_basis_species)
      real*8:: integrals(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2)
      real*8:: core_ham_old(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2)

      real*8::  i_r, wave_1, wave_2, wave_deriv, wave_kin
      integer:: i_radial

      eige = 1000.d0
      integrals = 0.d0

     ! Coefficient for the integrations
     alpha = log(r_grid_inc(species(i_atom)))

     eigen_old = 0.d0


      
      do i_basis_L2 = 1, n_basis_list,1
         do i_basis_L = 1, i_basis_L2,1


         i_basis_1    = basis_list(i_basis_L)
         i_basis_2    = basis_list(i_basis_L2)

         i_function_1 = basis_fn(i_basis_1)
         i_function_2 = basis_fn(i_basis_2)

         i_atom     = Cbasis_to_atom(i_basis_1)
         l_1 = basis_l(i_basis_1)
         l_2 = basis_l(i_basis_2)

         m_1 = basis_m(i_basis_1)
         m_2 = basis_m(i_basis_2)


         if(l_1 == l_2 .and. m_1 == m_2)then


         do i_radial = 1, n_radial(species(i_atom)), 1

   
            i_r  = invert_log_grid(  r_radial(i_radial, species(i_atom)), &
                 r_grid_min(species(i_atom)), r_grid_inc(species(i_atom)))
 
            
            wave_2 = val_spline( i_r,  basis_wave_spl(1,1,i_function_2), n_grid(species(i_atom)))
            wave_1 = val_spline( i_r,  basis_wave_spl(1,1,i_function_1), n_grid(species(i_atom)))

            wave_kin = val_spline( i_r,  basis_kinetic_spl(1,1,i_function_1), n_grid(species(i_atom)))

            wave_deriv = val_spline( i_r,  basis_deriv_spl(1,1,i_function_1), n_grid(species(i_atom)))


            V_radial = val_spline( i_r,  V_monopole_log_spl(1,1,species(i_atom)) , n_grid(species(i_atom))) &
                 - species_z(species(i_atom))/  r_radial(i_radial, species(i_atom))

            V_radial_deriv =  val_spline_deriv( i_r, V_monopole_log_spl(1,1, species(i_atom)), n_grid( species(i_atom))) &
!
                 / (log(r_grid_inc( species(i_atom))) &
                 * r_radial(i_radial, species(i_atom))) &
                + species_z(species(i_atom))/  (r_radial(i_radial, species(i_atom))**2)



            integrals(b_index(i_basis_L,i_basis_L2)) = integrals(b_index(i_basis_L,i_basis_L2)) &
                 + wave_2  & 
                 * w_radial(i_radial, species(i_atom)) * ( &
!
!
!
              wave_kin  * &
              ( 2 * light_speed_sq &
              / (  2 * light_speed_sq + eigen  -  V_radial)) &
              !
             - (light_speed_sq &
              / (  2*light_speed_sq + eigen -   V_radial)**2) &
!
!
              * V_radial_deriv &
               * (   wave_deriv  & 
               -    wave_1  &
               /  r_radial(i_radial, species(i_atom))))
!
!
!                 +  V_radial * wave_1* wave_2 * w_radial(i_radial, species(i_atom)) ! * r_radial(i_radial, species(i_atom))**2


         end do
         
!         integrals(b_index(i_basis_L,i_basis_L2)) = integrals(b_index(i_basis_L,i_basis_L2)) + core_ham_old(b_index(i_basis_L,i_basis_L2))


      else




      end if
   end do
   end do



 end subroutine integrate_KH_kinetic_term
!------------------------------------------------------


!-------------------------------------------------------------------------------------------------------------------------

    subroutine normalize(basis_list, n_basis_list, integrals, i_atom, eigen, V_monopole_log_spl, b_index )
      
      use spline
      use free_atoms

      integer:: i_basis_1, i_basis_2, i_function_1,i_function_2, i_atom,  i_basis_L, i_basis_L2,  l_1,l_2, m_1,m_2
      real*8::  V_radial_deriv, eige, V_radial
      integer:: n_basis_list
      integer:: basis_list(n_basis_list)
      real*8:: eigen, eigen_old
      real*8:: alpha
      integer:: temp(n_basis_list)
      real*8:: V_monopole_log_spl(4, n_max_grid, n_species )
      integer:: b_index(n_max_KH_core_basis_species, n_max_KH_core_basis_species)
      real*8:: integrals(n_max_KH_core_basis_species*(n_max_KH_core_basis_species+1)/2)

      real*8:: i_r, wave_1, wave_2, wave_deriv, wave_kin
      integer:: i_radial


      eige = 1000.d0
      integrals = 0.d0


     ! Coefficient for the integrations
     alpha = log(r_grid_inc(species(i_atom)))

     eigen_old = 0.d0




      
      do i_basis_L2 = 1, n_basis_list,1
         do i_basis_L = 1, i_basis_L2,1

         i_basis_1    = basis_list(i_basis_L)
         i_basis_2    = basis_list(i_basis_L2)

         i_function_1 = basis_fn(i_basis_1)
         i_function_2 = basis_fn(i_basis_2)

         i_atom     = Cbasis_to_atom(i_basis_1)
         l_1 = basis_l(i_basis_1)
         l_2 = basis_l(i_basis_2)

         m_1 = basis_m(i_basis_1)
         m_2 = basis_m(i_basis_2)


         if(l_1 == l_2 .and. m_1 == m_2)then


         do i_radial = 1, n_radial(species(i_atom)), 1

   
            i_r  = invert_log_grid(  r_radial(i_radial, species(i_atom)), &
                 r_grid_min(species(i_atom)), r_grid_inc(species(i_atom)))
 
            
            wave_2 = val_spline( i_r,  basis_wave_spl(1,1,i_function_2), n_grid(species(i_atom)))
            wave_1 = val_spline( i_r,  basis_wave_spl(1,1,i_function_1), n_grid(species(i_atom)))

            wave_kin = val_spline( i_r,  basis_kinetic_spl(1,1,i_function_1), n_grid(species(i_atom)))

            wave_deriv = val_spline( i_r,  basis_deriv_spl(1,1,i_function_1), n_grid(species(i_atom)))


            V_radial = val_spline( i_r,  V_monopole_log_spl(1,1,species(i_atom)) , n_grid(species(i_atom))) &
                 - species_z(species(i_atom))/  r_radial(i_radial, species(i_atom))

            V_radial_deriv =  val_spline_deriv( i_r, V_monopole_log_spl(1,1, species(i_atom)), n_grid( species(i_atom))) &
!
                 / (log(r_grid_inc( species(i_atom))) &
                 * r_radial(i_radial, species(i_atom))) &
                + species_z(species(i_atom))/  (r_radial(i_radial, species(i_atom))**2)



           integrals(b_index(i_basis_L,i_basis_L2)) = integrals(b_index(i_basis_L,i_basis_L2)) &
!
                 + wave_2  & 
                 * w_radial(i_radial, species(i_atom)) * ( &
!
!
!
              wave_kin  * &
              ( 2 * light_speed_sq &
              / (  2 * light_speed_sq + eigen  -  V_radial)**2) &
              !
             - (2*light_speed_sq &
              / (  2*light_speed_sq + eigen -   V_radial)**3) &
!
!
              * V_radial_deriv &
               * (   wave_deriv  & 
               -    wave_1  &
               /  r_radial(i_radial, species(i_atom))))
!
!



         end do
      end if
   end do
   end do


 end subroutine normalize
 





!----------------------------------------------------------------------------------------------------------------
  subroutine add_KH_states(KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex)

    real*8,     dimension(n_basis, n_states,n_spin,n_k_points_task) ::  KS_eigenvector
    complex*16, dimension(n_basis, n_states,n_spin,n_k_points_task) ::  KS_eigenvector_complex
    real*8,     dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue

    integer:: i_atom, i_state, i_s, i_b, i_basis
    integer::  i_k,i_k_point, i_spin
    real*8:: time1



    if(.not. flag_KH_core_states)   return


       call cpu_time(time1)


       do i_spin = 1, n_spin
          do i_k_point = 1, n_k_points,1
             i_s = 0
             do i_state = 1, n_max_KH_core_states_species
                do i_atom = 1, n_atoms

                   if(i_state <= n_KH_core_states_species(species(i_atom)))then
                      i_s = i_s + 1


                      KS_eigenvalue(i_s,i_spin,i_k_point) = Core_KH_eigenvalues( i_state, i_spin, i_atom)
                   


                   end if
                end do
             end do
          end do
       end do


       if(real_eigenvectors)then

          if(n_periodic == 0)then
             do i_spin = 1, n_spin

                i_s = 0
                do i_state = 1, n_max_KH_core_states_species
                   do i_atom = 1, n_atoms
                      if(i_state <= n_KH_core_states_species(species(i_atom)))then 
                   
                         i_b = 0
                         i_s = i_s + 1
                         do i_basis = 1, n_basis
                      
                            if(Cbasis_to_atom(i_basis) == i_atom)then
                               i_b = i_b + 1


                               KS_eigenvector(i_basis,i_s,i_spin,1) = Core_KH_eigenvectors( i_b, i_state, i_spin, i_atom)
                            else
                               KS_eigenvector(i_basis,i_s,i_spin,1) = 0.d0
                            end if
                         end do
                      end if
                   end do
                end do
             end do

          else

             do i_spin = 1, n_spin
                i_k = 0
                do i_k_point = 1, n_k_points,1
                   
                   if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                       
                      i_k = i_k + 1


                      i_s = 0
                      do i_state = 1, n_max_KH_core_states_species
                         do i_atom = 1, n_atoms
                            if(i_state <= n_KH_core_states_species(species(i_atom)))then 
                   
                               i_b = 0
                               i_s = i_s + 1
                               do i_basis = 1, n_basis
                      
                                  if(Cbasis_to_atom(i_basis) == i_atom)then
                                     i_b = i_b + 1


                                     KS_eigenvector(i_basis,i_s,i_spin,i_k) = Core_KH_eigenvectors( i_b, i_state, i_spin, i_atom)
                                  else
                                     KS_eigenvector(i_basis,i_s,i_spin,i_k) = 0.d0
                                  end if
                               end do
                            end if
                         end do
                      end do
                   end if
                end do
             end do
          end if

       else ! COMPLEX eigenvectors

          do i_spin = 1, n_spin
             i_k = 0
             do i_k_point = 1, n_k_points,1
                    
                if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                
                   i_k = i_k + 1



                   i_s = 0
                   do i_state = 1, n_max_KH_core_states_species
                      do i_atom = 1, n_atoms
                         if(i_state <= n_KH_core_states_species(species(i_atom)))then 

                            i_b = 0
                            i_s = i_s + 1
                            do i_basis = 1, n_basis
                            
                               if(Cbasis_to_atom(i_basis) == i_atom)then
                                  i_b = i_b + 1


                                  KS_eigenvector_complex(i_basis,i_s,i_spin,i_k) = & 
                                    Core_KH_eigenvectors( i_b, i_state, i_spin, i_atom)
                               else
                                  KS_eigenvector_complex(i_basis,i_s,i_spin,i_k) = (0.d0,0.d0)
                               end if
                            end do
                         end if
                      end do
                   end do
                end if
             end do
          end do
       end if





  end subroutine add_KH_states
!----------------------------------------------------------------------------------------------------------------
                      
  subroutine get_KH_free_atoms_potential(  V_monopole_log_spl )
    
    use spline
    use free_atoms


    real*8:: V_monopole_log_spl(4, n_max_grid, n_species )
    integer:: i_species, i_grid
    real*8:: V_monopole_log(n_max_grid)



    do i_species = 1, n_species
            
       do i_grid = 1,   n_grid(i_species)
          V_monopole_log(i_grid) = free_potential_spl(1,i_grid,i_species) + species_z(i_species)/ r_grid(i_grid, i_species)     
       end do
            
       call cubic_spline (   V_monopole_log(1), n_grid(i_species), &
            V_monopole_log_spl(1,1, i_species))
       
    end do



  end subroutine get_KH_free_atoms_potential

!----------------------------------------------------------------------------------------------------------------

    subroutine small_component(n_compute_c, i_basis, n_atom_list, dist_tab_sq, atom_index_inv, &
         atom_index, gradient_basis_wave, rho, n_basis_list,i_spin)

      use free_atoms

      integer:: n_points, n_compute_c, n_atom_list, i_center, n_basis_list
      integer:: i_basis(n_compute_c), atom_index_inv(n_centers), atom_index(n_atom_list)
      integer:: i_basis_L, i_basis_1, i_atom, i_point, i_state, i_core_basis, i_spin,i_center_L
      real*8:: dist_tab_sq(n_centers_integrals)
      real*8:: rr_radial, i_r, V_radial, eigen, vec
      real*8, dimension(n_basis_list, 3) :: gradient_basis_wave
      real*8,dimension(:,:,:), allocatable:: orb
      real*8 :: rho
      real*8:: V_monopole_log_spl(4, n_max_grid, n_species )


     if(.not. use_small_component) return

     call get_KH_free_atoms_potential(  V_monopole_log_spl )

     allocate(orb(n_max_KH_core_states_species, n_atom_list, 3))
      orb = 0.d0
      
      



      do i_basis_L = 1, n_compute_c,1


         i_basis_1  = i_basis(i_basis_L)

         i_core_basis = basis_to_core_basis(Cbasis_to_basis(i_basis_1))
         
         i_center   =  Cbasis_to_center(i_basis_1)
         i_atom     =  Cbasis_to_atom(i_basis_1)         
         i_center_L = atom_index_inv(i_center)
         
         if(i_center_L > 0)then
                        
               
               rr_radial = sqrt( dist_tab_sq(atom_index_inv(i_center)))

               i_r  = invert_log_grid(  rr_radial, &
                    r_grid_min(species(i_atom)), r_grid_inc(species(i_atom)))

               ! V_radial = val_spline( i_r, free_potential_spl(1,1,species(i_atom)), n_grid(species(i_atom))) 

               V_radial = val_spline( i_r,  V_monopole_log_spl(1,1,species(i_atom)) , n_grid(species(i_atom))) &
                    - species_z(species(i_atom))/  rr_radial
            
               do i_state = 1,  n_KH_core_states_species(species(i_atom))


                  vec =  Core_KH_eigenvectors(i_core_basis, i_state, i_spin, i_atom )
                  
                  eigen = Core_free_eigenvalues(i_state, i_spin, i_atom)


                  orb(i_state, i_center_L,1:3) =  orb(i_state, i_center_L, 1:3) + &
                       light_speed/( 2 * light_speed_sq + eigen - V_radial) * vec * gradient_basis_wave(i_basis_L, 1:3)

                  
               end do
         end if
      end do
      
      rho = 0.d0

      do i_center_L = 1,n_atom_list
         do i_state = 1, n_KH_core_states_species(species_center(atom_index(i_center_L)))
               

               if(n_spin == 1)then
                  rho = rho + 2*dot_product(orb(i_state, i_center_L,1:3), orb(i_state, i_center_L,1:3))
               else
                  rho = rho +  dot_product(orb(i_state, i_center_L,1:3), orb(i_state, i_center_L,1:3))
               end if

               
            end do
      end do

      deallocate(orb)


 end subroutine small_component
 










end module KH_core_states
