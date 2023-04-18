!---------------------------------------------------------------------
!******
!****s*prepare_embedding_grids
!  NAME
!   prepare_embedding_grids
!  SYNOPSIS
     subroutine prepare_embedding_grids( )
!      PURPOSE
!      this routine checks which of the pseudocores is too close to the qm region, 
!      and therefore should get a grid
!      IMPORTANT: this routine should be called also in reintialize_scf()
!
!
!
!      USES
       use dimensions
       use geometry
       use free_atoms
       use species_data
       use pseudodata
       use basis, only: atom_radius_sq

!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

       integer :: dummy, i, j, i_pp_atom, i_atom
       real*8, dimension(:,:) :: coords_old(3,n_real_atoms)
       integer, dimension(:) :: species_old(n_real_atoms)
       logical, dimension(:) :: empty_old(n_real_atoms)
       real*8 :: initial_rho_spl_old(n_max_spline, n_max_grid, &
                n_ini_type, n_spin)
       integer, dimension(:) :: atom_type_old(n_real_atoms)
       real*8 :: initial_drho_dr_spl_old (n_max_spline, n_max_grid, &
                  n_ini_type, n_spin)

    real*8 :: dis
    integer :: i_multipole
    integer, dimension(n_multipoles) :: offending_multipole


       character*130 :: info_str

          n_atoms = n_real_atoms          

! preparation for pseudocores
          if (allocated(pp_in_qm)) then
            deallocate(pp_in_qm)
          end if
          allocate(pp_in_qm(n_pp_atoms))


  ! count the number of pseudocores which should better be equipped with some grids

          call get_pp_distance2atom()

  ! and now we do the same for multipoles if any present

          if (n_multipoles.gt.0) then
               call min_multipole_dist( lattice_vector, coords, species, atom_radius_sq, multipole_coords, & 
                   &                        empty_coords, dis, i_atom, i_multipole, offending_multipole)
          end if


   !      write(use_unit,*) 'offending_multipole', offending_multipole
   !      write(use_unit,*) 'maxval(offending_multipole)', maxval(offending_multipole)

          n_mp_in_qm = 0

          do i_multipole = 1, n_multipoles
              if(offending_multipole(i_multipole).eq.0) exit
              n_mp_in_qm = n_mp_in_qm + 1
          enddo
 
         write(use_unit,*) 'n_mp_in_qm', n_mp_in_qm
         

          if ((n_pp_in_qm.gt.0).or.(n_mp_in_qm.gt.0)) then
  
             n_atoms = n_real_atoms + n_pp_in_qm

             if(n_pp_in_qm.gt.0) then
               write(info_str,'(2X, A, I3,2x, 2A)') 'A number of ',n_pp_in_qm,'pseudocores are too close to the ',&
                    'QM-region, and have been added to n_atoms. '
               call localorb_info(info_str) 
               write(info_str,'(2X,A, I4)') 'Updated number of atoms:', n_atoms
               call localorb_info(info_str) 

               write(info_str,*) ''
               call localorb_info(info_str)
             endif

             n_atoms = n_real_atoms + n_pp_in_qm + n_mp_in_qm
             
             if(n_mp_in_qm.gt.0) then
               write(info_str,'(2X, A, I3,2x, 2A)') 'A number of ',n_mp_in_qm,'monopoles are too close to the ',&
                    'QM-region. A ghost atom of species "Emptium" is placed on top of those monopoles.'
               call localorb_info(info_str) 
               write(info_str,'(2X,A, I4)') 'Updated number of atoms:', n_atoms
               call localorb_info(info_str) 

               write(info_str,*) ''
               call localorb_info(info_str)

             endif

!   re-allocate species and add pp_species accordingly

!   for all the pseudoions touching the qm zone,
!   we need to re-allocate following arrays:
!     - species(n_atoms)
!     - coords(3,n_atoms)
!     - empty(n_atoms)

             if(allocated(pp_atom2atom)) then
               deallocate(pp_atom2atom)
             end if
             allocate(pp_atom2atom(n_atoms))             
             pp_atom2atom = -1 

             species_old = species

             if(allocated(species)) then
               deallocate(species)
             end if
             allocate(species(n_atoms))
          
             do i = 1, n_real_atoms
                species(i) = species_old(i) 
             end do

             i=0
             do i_pp_atom = 1, n_pp_atoms
                if (pp_in_qm(i_pp_atom)) then
                i = i + 1
                species(n_real_atoms + i) = pp_species2species(pp_species(i_pp_atom))
                endif
             enddo


!   re-allocate coords and add pp_coords accordingly


             coords_old = coords

             if(allocated(coords)) then
               deallocate(coords)
             end if
             allocate(coords(3,n_atoms))

             do i = 1, n_real_atoms
                coords(1:3,i) = coords_old(1:3,i) 
             end do

             i=0
             do i_pp_atom = 1, n_pp_atoms
                if (pp_in_qm(i_pp_atom)) then
                i = i + 1
                coords(1:3,n_real_atoms + i) = pp_coords(1:3,i_pp_atom)
                endif
             enddo



             empty_old = empty

             if(allocated(empty)) then
               deallocate(empty)
             end if
             allocate(empty(n_atoms))

             do i = 1, n_real_atoms
                empty(i) = empty_old(i) 
             end do


             do i = 1, n_pp_in_qm
                empty(n_real_atoms + i) = .false.
             enddo


! for the spin polarized case we have to reallocate some more arrays
             if(use_initial_rho) then
 
                atom_type_old = atom_type

                if(allocated(atom_type)) then
                  deallocate(atom_type)
                end if

                allocate(atom_type(n_atoms))

                do i = 1, n_real_atoms
                   atom_type(i) = atom_type_old(i)
                enddo


                j = 0
                do i_pp_atom = 1, n_pp_atoms
                   if (pp_in_qm(i_pp_atom)) then
                      j = j + 1
                      atom_type(n_real_atoms + j) = n_ini_type + pp_species(i_pp_atom)
                   end if
                enddo


                initial_rho_spl_old = initial_rho_spl

                if(allocated(initial_rho_spl)) then
                  deallocate(initial_rho_spl)
                end if

                allocate (initial_rho_spl(n_max_spline, n_max_grid, &
                n_ini_type + n_pp_species, n_spin))
                
                do i= 1, n_ini_type 
                   initial_rho_spl(:,:,i,:) = initial_rho_spl_old(:,:,i,:) 
                enddo 

                do i = 1,n_pp_species
                   initial_rho_spl(:,:,n_ini_type+i,:) = 0.d0 
                enddo

                if (use_density_gradient) then
                    initial_drho_dr_spl_old = initial_drho_dr_spl


                   if(allocated(initial_drho_dr_spl)) then
                     deallocate(initial_drho_dr_spl)
                   end if

                   allocate (initial_drho_dr_spl(n_max_spline, n_max_grid, &
                        n_ini_type + n_pp_species, n_spin))
 

                   do i= 1, n_ini_type 
                     initial_drho_dr_spl(:,:,i,:) = initial_drho_dr_spl_old(:,:,i,:) 
                   enddo 

                   do i = 1,n_pp_species
                     initial_drho_dr_spl(:,:,n_ini_type+i,:) = 0.d0 
                   enddo


                end if



             end if ! use_initial_rho


             i_pp_atom = 0
             do i_atom = 1, n_atoms
                if (species_pseudoized(species(i_atom))) then
                  i_pp_atom = i_pp_atom + 1
                  if(pp_in_qm(i_pp_atom)) then
                    pp_atom2atom(i_atom) = i_pp_atom 
                  endif
                endif
             enddo


          elseif(n_pp_in_qm.eq.0) then

!             use_embedding_pp = .false.

          endif !(n_pp_in_qm.gt.0)


 
     end subroutine prepare_embedding_grids

