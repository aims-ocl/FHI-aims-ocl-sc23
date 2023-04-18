!****h* FHI-aims/pseudodata
!  NAME
!     pseudodata
!  SYNOPSIS

      module pseudodata

!  PURPOSE
!
! 
!
!  Subroutines:
!  * allocate_pseudoarrays
!  * get_pseudodata
!  * spline_and_sort_pseudodata
!  * read_pseudopotdata
!  * divide_pseudopot 
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2012).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  SOURCE


!      use

      use mpi_tasks
      use localorb_io

      implicit none

!  global variable declarations

      real*8, dimension(:), allocatable :: pp_charge
      integer, dimension(:), allocatable :: pp_local_component

      real*8, dimension(:,:), allocatable :: pseudo_wave
      real*8, dimension(:,:), allocatable :: pseudo_pot
      real*8, dimension(:,:), allocatable :: pseudo_grid
      real*8, dimension(:,:), allocatable :: pseudo_chi    
!  pseudo_chi = pseudo_wave*nonlocal_pot/(pseudo_grid*||pseudo_wave*nonlocal_pot||^0.5)

      real*8, dimension(:,:,:), allocatable :: pseudo_wave_spl
      real*8, dimension(:,:,:), allocatable :: pseudo_chi_spl
!      real*8, dimension(:,:,:), allocatable :: nonlocal_pseudopot_spl
      real*8, dimension(:,:,:), allocatable :: local_pseudopot_spl
      integer, dimension(:), allocatable :: pp_basis_atom
      integer, dimension(:), allocatable :: pp_basis_species 
      integer, dimension(:), allocatable :: pp_basis_l
      integer, dimension(:), allocatable :: pp_basis_m
      integer, dimension(:), allocatable :: pp_basis_fn
      integer, dimension(:), allocatable :: pp_basisfn_l
      integer, dimension(:), allocatable :: pp_basisfn_species

      real*8, dimension(:), allocatable :: pp_r_grid_min
      real*8, dimension(:), allocatable :: pp_r_grid_inc

      real*8, dimension(:), allocatable :: pp_outer_radius
  
      real*8, dimension(:), allocatable :: pp_psi_chi_overlap, E_l_KB
      real*8, dimension(:,:), allocatable :: basiswave_pp_overlap
      integer, dimension(:,:), allocatable :: nonzero_overlap_entries
      integer, dimension(:), allocatable :: n_max_basis_overlap
      
! if needed, the derivative of basiswave_pp_overlap with respect
! to the position of the pseudocores
      real*8, dimension(:,:,:), allocatable :: d_basiswave_pp_overlap      

      real*8, dimension(:,:), allocatable :: pp_basiswave_overlap      
      real*8, dimension(:,:,:), allocatable :: d_pp_basiswave_overlap      


      integer, dimension(:,:), allocatable :: Lsp2n_pp_basis_fnLsp
      integer, dimension(:,:,:), allocatable :: Lsp2_pp_basis_fn
!      integer, dimension(:,:,:), allocatable :: pp_fnL_to_col


      real*8, dimension(:,:), allocatable :: nonlocal_pseudopot
      real*8, dimension(:,:), allocatable :: local_pseudopot 
      real*8, dimension(:,:,:), allocatable :: partial_core_dens_spl
      real*8, dimension(:,:,:), allocatable :: partial_core_dens_deriv_spl


      real*8, dimension(:), allocatable :: whole_local_pseudpot_on_intgrid
      real*8, dimension(:), allocatable :: partial_core_rho
      real*8, dimension(:,:), allocatable :: partial_core_rho_grad
      real*8, dimension(:), allocatable :: pp_norm
      real*8, dimension(:), allocatable :: localpot_outer_radius
      real*8, dimension(:,:), allocatable :: nonlocal_matrix
      real*8 :: en_nonlocal

      logical, dimension(:), allocatable :: pp_in_qm
      integer, dimension(:), allocatable :: pp_atom2atom

      integer, dimension(:,:), allocatable :: pp_function_index

      real*8, dimension(:,:,:,:), allocatable :: pseudo_fctn_and_pot  

      real*8, dimension(:,:), allocatable :: partial_core_dens
      real*8, dimension(:,:), allocatable :: partial_core_dens_deriv  

      integer, dimension(:), allocatable :: n_pp_fns      !number of different radial parts of pseudowavefunctions of pp_species
      integer, dimension(:), allocatable :: n_points_pp_fn    !number of gridpoints the pseudowavefunctions are defined on  

      integer, dimension(:), allocatable :: species2pp_species
      integer, dimension(:), allocatable :: pp_species2species

      integer :: max_n_pp_basis_fnLsp, max_pp_basis_l

      real*8, dimension(:,:), allocatable :: pp_nonlocal_forces
      real*8, dimension(:,:), allocatable :: pp_nonlocal_forces_2
      real*8, dimension(:,:), allocatable :: pp_nlcc_forces

      real*8, dimension(:,:,:,:), allocatable :: pp_nonlocal_forces_matrix

      real*8 :: en_xc_nlcc


      integer :: n_pp_basis

      integer :: n_pseudofn           !number of pseudowavefunction of certain species
      integer :: n_points_pseudofn    !number of gridpoints the pseudopotential are defined on  

      integer :: n_pp_species_counter_read_in

      contains
!---------------------------------------------------------------------
!******
!****s* pseudodata/prepare_pseudodata
!  NAME
!   prepare_pseudodata
!  SYNOPSIS
     subroutine prepare_pseudocoregrids( )
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
       use pbc_lists, only: check_multipole_distances

!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

       integer :: dummy, i, j, i_pp_atom, i_atom
       real*8, dimension(:,:) :: coords_old(3,n_real_atoms)
       integer, dimension(:) :: coord_basis_atoms_old(n_real_atoms)
       integer, dimension(:) :: species_old(n_real_atoms)
       logical, dimension(:) :: empty_old(n_real_atoms)
       real*8 :: initial_rho_spl_old(n_max_spline, n_max_grid, &
                n_ini_type, n_spin)
       integer, dimension(:) :: atom_type_old(n_real_atoms)
       real*8 :: initial_drho_dr_spl_old (n_max_spline, n_max_grid, &
                  n_ini_type, n_spin)

       character*130 :: info_str


!       if(myid==0) then

          n_atoms = n_real_atoms          

! preparation for pseudocores
          if (allocated(pp_in_qm)) then
            deallocate(pp_in_qm)
          end if
          allocate(pp_in_qm(n_pp_atoms))
          pp_in_qm = .false.

!          if(use_qmmm) then
            ! we need to treat all pseudocores as atoms, otherwise 
            ! forces are not printed correctly for chemshell 
            n_pp_in_qm = n_pp_atoms
            
            pp_in_qm = .true.

!          else
    
!             call get_pp_distance2atom()

!          end if 

     

 
          if (n_pp_in_qm.gt.0) then
  
             n_atoms = n_real_atoms + n_pp_in_qm

             write(info_str,'(2X, A, I3,2x, 2A)') 'A number of ',n_pp_in_qm,'pseudocores are too close to the ',&
                    'QM-region, and have been added to n_atoms. '
             call localorb_info(info_str) 
             write(info_str,'(2X,A, I4)') 'Updated number of atoms:', n_atoms
             call localorb_info(info_str) 

             write(info_str,*) ''
             call localorb_info(info_str)


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

!--------
             coord_basis_atoms_old = coord_basis_atoms

             if(allocated(coord_basis_atoms)) then
               deallocate(coord_basis_atoms)
             end if
             allocate(coord_basis_atoms(n_atoms))

             do i = 1, n_real_atoms
                coord_basis_atoms(i) = coord_basis_atoms_old(i) 
             end do

             i=0
             do i_pp_atom = 1, n_pp_atoms
                if (pp_in_qm(i_pp_atom)) then
                i = i + 1
                ! 0 means cartesian coordinates
                coord_basis_atoms(n_real_atoms + i) = 0
                endif
             enddo



! ------

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
!             if(use_initial_rho) then
! 
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


!                initial_rho_spl_old = initial_rho_spl

!                if(allocated(initial_rho_spl)) then
!                  deallocate(initial_rho_spl)
!                end if

!                allocate (initial_rho_spl(n_max_spline, n_max_grid, &
!                n_ini_type + n_pp_species, n_spin))
                
!                do i= 1, n_ini_type 
!                   initial_rho_spl(:,:,i,:) = initial_rho_spl_old(:,:,i,:) 
!                enddo 

!                do i = 1,n_pp_species
!                   initial_rho_spl(:,:,n_ini_type+i,:) = 0.d0 
!                enddo

!                if (use_density_gradient) then
!                    initial_drho_dr_spl_old = initial_drho_dr_spl


!                   if(allocated(initial_drho_dr_spl)) then
!                     deallocate(initial_drho_dr_spl)
!                   end if

!                   allocate (initial_drho_dr_spl(n_max_spline, n_max_grid, &
!                        n_ini_type + n_pp_species, n_spin))
 

!                   do i= 1, n_ini_type 
!                     initial_drho_dr_spl(:,:,i,:) = initial_drho_dr_spl_old(:,:,i,:) 
!                   enddo 

!                   do i = 1,n_pp_species
!                     initial_drho_dr_spl(:,:,n_ini_type+i,:) = 0.d0 
!                   enddo


!                end if



!             end if ! use_initial_rho


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

!       endif! (myid==0)


     end subroutine prepare_pseudocoregrids


!******

!---------------------------------------------------------------------
!******
!****s* pseudodata/prepare_pseudodata
!  NAME
!   prepare_pseudodata
!  SYNOPSIS
     subroutine prepare_more_pseudocoregrids( )
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

!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

       integer :: dummy, i, j, i_pp_atom, i_atom
       real*8 :: initial_rho_spl_old(n_max_spline, n_max_grid, &
                n_ini_type, n_spin)
       real*8 :: initial_drho_dr_spl_old (n_max_spline, n_max_grid, &
                  n_ini_type, n_spin)

       character*130 :: info_str


!       if(myid==0) then

          n_atoms = n_real_atoms          

          if (n_pp_in_qm.gt.0) then
  
             n_atoms = n_real_atoms + n_pp_in_qm


! for the spin polarized case we have to reallocate some more arrays
             if(use_initial_rho) then
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



          endif !(n_pp_in_qm.gt.0)

!       endif! (myid==0)

 
     end subroutine prepare_more_pseudocoregrids


!---------------------------------------------------------------------
!******
!****s* pseudodata/prepare_pseudodata
!  NAME
!   prepare_pseudodata
!  SYNOPSIS
     subroutine reprepare_pseudocoregrids( )
!      PURPOSE
!      this routine checks which of the pseudocores is too close to the qm region, 
!      and therefore should get a grid
!
!
!      USES
       use dimensions
       use geometry
       use species_data
       use physics

!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

       integer :: i_pp_atom, i_atom

       if (use_embedding_pp) then

       ! first, update pp_coords
       do i_atom = 1, n_atoms
        if(species_pseudoized(species(i_atom))) then
          i_pp_atom = pp_atom2atom(i_atom)
          pp_coords(:,i_pp_atom) =   coords(:,i_atom)     
        endif
       enddo
 
       ! set forces to zero:

       pp_nonlocal_forces = 0.d0

       pp_nonlocal_forces_2 = 0.d0
       pseudocore_forces = 0.d0

       if(allocated(pp_nlcc_forces)) then
          pp_nlcc_forces = 0.d0 
       endif
       if(allocated(nlcc_forces)) then
         nlcc_forces = 0.d0
       endif

       call prepare_pseudocoregrids( )
       end if

 
     end subroutine reprepare_pseudocoregrids


!******

!---------------------------------------------------------------------
!******
!****s* pseudodata/get_pseudodata
!  NAME
!   get_pseudodata
!  SYNOPSIS
     subroutine get_pseudodata( )
!      PURPOSE
!      this is just a wrapper around the allocate and spline_and_sort pseudodata
!
!      USES
       use dimensions
!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

       integer :: dummy

       call allocate_pseudoarrays()
       call spline_and_sort_pseudodata()
       if (use_nonlinear_core) then
            call check_partial_core_norm()
       endif

     end subroutine get_pseudodata


!******

!---------------------------------------------------------------------
!****s* pseudodata/allocate_pseudoarrays
!  NAME
!    allocate_pseudoarrays
!  SYNOPSIS

     subroutine allocate_pseudospecies_arrays( )

!      PURPOSE:
!       all arrays with dimension n_pp_species are allocated here

       use dimensions

      
       n_pp_species_counter_read_in = 0

       allocate(species_nonlinear_core(n_species))
       species_nonlinear_core  = .false.      

       allocate( species2pp_species(n_species))
       species2pp_species = 0
       allocate( pp_species2species(n_pp_species))
       pp_species2species = 0
       allocate( n_pp_fns(n_pp_species))
       n_pp_fns = 0
       allocate( n_points_pp_fn(n_pp_species))
       n_points_pp_fn = 0

       allocate(pp_path_string(n_species))
       pp_path_string = ''
      
!      allocate (pp_atoms_in_structure(n_pp_species))

       allocate (pp_charge(n_pp_species))
       pp_charge = 0d0       
       allocate (pp_local_component(n_pp_species))
       pp_local_component = 0d0

       allocate ( pp_r_grid_min(n_pp_species) )
       allocate ( pp_r_grid_inc(n_pp_species) )

       pp_r_grid_min = 0d0
       pp_r_grid_inc = 0d0         

       if(use_nonlinear_core) then
             allocate(partial_core_dens(n_max_points_pp_fn, n_pp_species))
             allocate(partial_core_dens_deriv(n_max_points_pp_fn, n_pp_species))
             partial_core_dens = 0.0d0
             partial_core_dens_deriv = 0.0d0
       endif

       allocate(pseudo_fctn_and_pot(n_pp_species,0:n_max_pp_fn-1, 3,n_max_points_pp_fn))
       pseudo_fctn_and_pot = 0d0

       allocate(pp_function_index(n_pp_species,n_max_pp_fn))


     end subroutine allocate_pseudospecies_arrays

!---------------------------------------------------------------------
!****s* pseudodata/allocate_pseudoarrays
!  NAME
!    allocate_pseudoarrays
!  SYNOPSIS

     subroutine allocate_pseudoarrays( )

!      PURPOSE:
!       all arrays concerning pseudopotential data are allocated here, besides the
!       compact read-in-array pseudo_fctn_and_pot - see species_data.f90



       use dimensions
       use runtime_choices
       use geometry
       use basis
       use species_data
       use mpi_tasks
       use physics
!       AUTHOR
!         FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!       SEE ALSO
!         Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!         Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!         "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!         Computer Physics Communications (2008), submitted.
!       COPYRIGHT
!        Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!        e.V. Please note that any use of the "FHI-aims-Software" is subject to
!        the terms and conditions of the respective license agreement."
!       HISTORY
!         Release version, FHI-aims (20012).
!       SOURCE
       implicit none
       !  local variables
       integer :: i, i_pp_basis_fns, i_pp_basis, i_l, i_m
       integer :: info, j, i_pp_species, i_pp_atom


         allocate ( pseudo_wave(n_max_points_pp_fn, n_pp_basis_fns) )
         allocate ( pseudo_pot(n_max_points_pp_fn, n_pp_basis_fns) )
         allocate ( pseudo_grid(n_max_points_pp_fn, n_pp_basis_fns) )
         allocate ( pseudo_chi(n_max_points_pp_fn, n_pp_basis_fns) )
         pseudo_wave = 0.d0 
         pseudo_pot = 0.d0  
         pseudo_grid = 0.d0  
         pseudo_chi = 0.d0  

         allocate ( pp_basisfn_l(n_pp_basis_fns))
         allocate ( pp_basisfn_species(n_pp_basis_fns))
         allocate ( pp_outer_radius(n_pp_basis_fns))
         pp_basisfn_l = 0
         pp_basisfn_species = 0
         pp_outer_radius = 0

         i_pp_basis=0
         do i_pp_atom = 1, n_pp_atoms
             i_pp_species = pp_species(i_pp_atom)
             do i_l = 0, n_pp_fns(i_pp_species)-1
                 do i_m = -i_l, i_l
                     i_pp_basis = i_pp_basis + 1
                 end do
             end do
         enddo
         n_pp_basis = i_pp_basis



         allocate ( pp_basis_atom(n_pp_basis))
         allocate ( pp_basis_species(n_pp_basis) )
         allocate ( pp_basis_l(n_pp_basis) )
         allocate ( pp_basis_m(n_pp_basis) )
         allocate ( pp_basis_fn(n_pp_basis) )

         pp_basis_atom = 0
         pp_basis_species = 0  
         pp_basis_l = 0  
         pp_basis_fn = 0 


         if (.not.allocated(pseudo_wave_spl)) then
             allocate ( pseudo_wave_spl(n_max_spline,n_max_points_pp_fn,n_pp_basis_fns),stat=info )
!           call check_allocation(info, 'pseudo_wave_spl                ')
         end if
         if (.not.allocated(pseudo_chi_spl)) then
             allocate ( pseudo_chi_spl(n_max_spline,n_max_points_pp_fn,n_pp_basis_fns),stat=info )
!           call check_allocation(info, 'pseudo_wave_spl                ')
         end if
         if (.not.allocated(local_pseudopot_spl)) then
             allocate ( local_pseudopot_spl(n_max_spline,n_max_points_pp_fn,n_pp_species),stat=info )
!           call check_allocation(info, 'pseudo_pot_spl                ')
         end if
       

         local_pseudopot_spl = 0.d0 
         pseudo_wave_spl = 0.d0
         pseudo_chi_spl = 0.d0

         if (.not.allocated(pp_psi_chi_overlap)) then
             allocate(pp_psi_chi_overlap(n_pp_basis_fns))
         endif
         pp_psi_chi_overlap = 0.d0

         if (.not.allocated(E_l_KB)) then
             allocate(E_l_KB(n_pp_basis_fns))
         endif
         E_l_KB = 0.d0

         if (.not.allocated(basiswave_pp_overlap)) then
             allocate(basiswave_pp_overlap(n_basis, n_pp_basis))
!         allocate(basiswave_pp_overlap(n_pp_basis, n_basis))
         endif
         basiswave_pp_overlap = 0.d0

! DB 120112 remove after debugging
         if (.not.allocated(pp_basiswave_overlap)) then
             allocate(pp_basiswave_overlap(n_pp_basis, n_basis))
         endif
         pp_basiswave_overlap = 0.d0

         if (.not.allocated(d_pp_basiswave_overlap)) then
             allocate(d_pp_basiswave_overlap(3, n_pp_basis, n_basis))
         endif
         d_pp_basiswave_overlap = 0.d0

         if(use_forces) then
            if (.not.allocated(d_basiswave_pp_overlap)) then
!            allocate(d_basiswave_pp_overlap(3, n_pp_basis, n_basis))
                allocate(d_basiswave_pp_overlap(3, n_basis, n_pp_basis))
            endif
            d_basiswave_pp_overlap = 0.d0

            if (.not.allocated(pp_nonlocal_forces)) then
                allocate(pp_nonlocal_forces(3,n_pp_atoms))
            endif
            pp_nonlocal_forces=0.d0

            if (.not.allocated(pp_nonlocal_forces_2)) then
                allocate(pp_nonlocal_forces_2(3,n_real_atoms))
            endif
            pp_nonlocal_forces_2=0.d0

            if (.not.allocated(pp_nlcc_forces)) then
                allocate(pp_nlcc_forces(3,n_pp_atoms))
            endif
            pp_nlcc_forces=0.d0


         endif

         if (.not.allocated(nonlocal_pseudopot)) then
             allocate(nonlocal_pseudopot(n_max_points_pp_fn, n_pp_basis_fns))
         endif
         if (.not.allocated(local_pseudopot)) then
             allocate(local_pseudopot(n_max_points_pp_fn, n_pp_species))
         endif

         nonlocal_pseudopot = 0.0d0
         local_pseudopot = 0.0d0

         if (.not.allocated(pp_norm)) then
             allocate ( pp_norm(n_pp_basis_fns))
         endif

         if (.not.allocated(nonlocal_matrix)) then
             allocate ( nonlocal_matrix(n_basis, n_basis))
         endif

         if (.not.allocated(localpot_outer_radius)) then
             allocate ( localpot_outer_radius(n_pp_species))
         endif        

!	 allocate ( pp_in_qm(n_pp_atoms))

         if(use_nonlinear_core) then
            if (.not.allocated(partial_core_dens_spl)) then
                allocate ( partial_core_dens_spl(n_max_spline,n_max_points_pp_fn,n_pp_species),stat=info )
            end if

            if (.not.allocated(partial_core_dens_deriv_spl)) then
                allocate ( partial_core_dens_deriv_spl(n_max_spline,n_max_points_pp_fn,n_pp_species),stat=info )
            end if

         endif

        end subroutine allocate_pseudoarrays
!******

!----------------------------------------------------------------------------
!****s* pseudodata/spline_and_sort_pseudodata
!  NAME
!    spline_and_sort_pseudodata
!  SYNOPSIS

      subroutine spline_and_sort_pseudodata()

      !  PURPOSE
      !
      !  Splines and sorts pseudo wave functions and the pseudopotentials    
      !  Input is the compressed array pseudo_fctn_and_pot
      !  Arrays pseudo_wave & pseudo_wave_spl are parsed in and allocated
      !    and finally loaded with data.    


      ! USES      
      use species_data
      use dimensions
      use spline
      use geometry
!       use basis

      implicit none
 

      !  INPUTS
      !    o pseudo_fctn_and_pot -- array of whole pseudo input data
      !
      !  OUTPUTS
      !    o pseudo_wave_spl(n_max_spline,n_max_points_pp_fn,n_pp_basis_fns) -- array of splined pseudo wave fctn u = r*f
      !    o pseudo_pot_spl(n_max_spline,n_max_points_pp_fn,n_pp_basis_fns)  -- array of splined pseudo pot fctn
      !    o pp_basis_species(n_pp_basis) -- species i_pp_basis belongs to
      !    o pp_basis_l(n_pp_basis) -- angular momentum of i_pp_basis
      !    o pp_basis_m(n_pp_basis) -- z-component of the angular momentum
      !
      !    o pp_basisfn_l(n_pp_basis_fns) -- angular momentum of the distinctive radial parts i_pp_basis_fns
      !    o pp_basisfn_species(n_pp_basis_fns) -- species i_pp_basis_fns belongs to
      !    o pp_r_grid_min(n_pp_species) -- rmin 
      !    o pp_r_grid_inc(n_pp_species) -- incrementation of loggrid
      !
      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2011).
      !  SOURCE


!     local variables
      integer, dimension(:,:,:), allocatable :: Lsp2_pp_basis_sp


!     counter
      integer :: i_pp_species, i_pp_basis, i_pp_atom, i, i_fnLsp
      integer :: i_l, i_m, N, i_pp_basis_fns, n_pp_fnLsp, L
      integer :: dummy
      character(*), parameter :: func = 'spline_and_sort_pseudodata'

! --- sort & read data ----


      i_pp_basis_fns = 0 
      i_pp_species = 1
      do i_pp_species = 1, n_pp_species
          do i_l = 0,n_pp_fns(i_pp_species)-1
             i_pp_basis_fns = i_pp_basis_fns + 1
             pp_basisfn_l(i_pp_basis_fns) = i_l 
             pp_basisfn_species(i_pp_basis_fns) = i_pp_species
             pp_function_index(i_pp_species, i_l + 1 ) = i_pp_basis_fns
          enddo
      enddo


      max_pp_basis_l = maxval(pp_basisfn_l)

      i_pp_basis = 0 

      
      do i_pp_atom = 1,n_pp_atoms
        do i_l = 0, n_pp_fns(pp_species(i_pp_atom)) -1
           do i_m = -i_l, i_l 
           i_pp_basis = i_pp_basis + 1
           pp_basis_atom(i_pp_basis) = i_pp_atom
           pp_basis_species(i_pp_basis) = pp_species(i_pp_atom)
           pp_basis_l(i_pp_basis) = i_l
           pp_basis_m(i_pp_basis) = i_m
           pp_basis_fn(i_pp_basis) = pp_function_index(pp_species(i_pp_atom), i_l+1)
           end do
        end do
      end do


      do i_pp_species = 1, n_pp_species
        do i_l = 0, max_pp_basis_l
           n_pp_fnLsp = count(pp_basisfn_l == i_l .and. pp_basisfn_species == i_pp_species)
           max_n_pp_basis_fnLsp = max(max_n_pp_basis_fnLsp, n_pp_fnLsp)
        end do
      end do

      allocate (Lsp2n_pp_basis_fnLsp(0:max_pp_basis_l, n_pp_species))
      allocate (Lsp2_pp_basis_fn(max_n_pp_basis_fnLsp,0:max_pp_basis_l, n_pp_species))


      Lsp2n_pp_basis_fnLsp = 0
      Lsp2_pp_basis_fn = 0
      do i_pp_basis_fns = 1, n_pp_basis_fns
        i_pp_species = pp_basisfn_species(i_pp_basis_fns)
        L = pp_basisfn_l(i_pp_basis_fns)
        Lsp2n_pp_basis_fnLsp(L, i_pp_species) = Lsp2n_pp_basis_fnLsp(L, i_pp_species) + 1
        i_fnLsp = Lsp2n_pp_basis_fnLsp(L, i_pp_species)
        if (i_fnLsp > max_n_pp_basis_fnLsp) then
           call aims_stop('max_n_bas_fnLsp too small', func)
        else if (any(Lsp2_pp_basis_fn == i_pp_basis_fns)) then
           call aims_stop('Double usage of i_pp_basis_fn', func)
        end if
        Lsp2_pp_basis_fn(i_fnLsp, L, i_pp_species) = i_pp_basis_fns
      end do
      if (sum(Lsp2n_pp_basis_fnLsp) /= n_pp_basis_fns) then
         call aims_stop('sum(Lsp2n_pp_basis_fnLsp) mismatch', func)
      end if


! indexing is done, now:
!          putting data on arrays

      do i_pp_basis_fns = 1,n_pp_basis_fns
         i_pp_species = pp_basisfn_species(i_pp_basis_fns)
         do i = 1, n_points_pp_fn(i_pp_species)
            i_l = pp_basisfn_l(i_pp_basis_fns) 
            pseudo_wave(i, i_pp_basis_fns) = pseudo_fctn_and_pot(i_pp_species,i_l,2,i)
            pseudo_pot(i, i_pp_basis_fns) = pseudo_fctn_and_pot(i_pp_species,i_l,3,i)
            pseudo_grid(i, i_pp_basis_fns) = pseudo_fctn_and_pot(i_pp_species,i_l,1,i)
         end do
      end do


! folllowing data is needed in put_pot_on_integration_grid
      do i_pp_basis_fns = 1,n_pp_basis_fns
          pp_r_grid_min(pp_basisfn_species(i_pp_basis_fns)) = pseudo_grid(1, i_pp_basis_fns)
          pp_r_grid_inc(pp_basisfn_species(i_pp_basis_fns)) = &
             pseudo_grid(2, i_pp_basis_fns)/pseudo_grid(1, i_pp_basis_fns)
      enddo
 

      call divide_pseudopot()


      do i_pp_basis_fns = 1, n_pp_basis_fns
         if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
          pseudo_chi(:,i_pp_basis_fns) = & 
            pseudo_wave(:,i_pp_basis_fns) * nonlocal_pseudopot(:,i_pp_basis_fns)
         end if
      end do

!  \chi = u * \delta V

      call get_pp_norm()  



! changed at dd/mm/yy 12/04/13 
      do i_pp_basis_fns = 1, n_pp_basis_fns


        if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
          pseudo_chi(:,i_pp_basis_fns) = & 
             pseudo_chi(:,i_pp_basis_fns) / (pp_norm(i_pp_basis_fns))**(1./2.) 
        end if
      end do


! \chi = u*\delta V / sqrt{\int dr (u*\deltaV)^2}

! --- spline data
      do i_pp_species = 1, n_pp_species
         call cubic_spline ( local_pseudopot(1,i_pp_species), &
                             n_points_pp_fn(i_pp_species), &
                             local_pseudopot_spl(1,1,i_pp_species) &
                            )
      enddo


      do i_pp_species = 1, n_pp_species
         localpot_outer_radius(i_pp_species) = & 
             maxval(pseudo_fctn_and_pot(i_pp_species,:,1,:))
      end do

! ----------------------
      do i_pp_basis_fns = 1, n_pp_basis_fns
         call cubic_spline ( pseudo_wave(1,i_pp_basis_fns), &
                             n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns)), &
                             pseudo_wave_spl(1,1,i_pp_basis_fns) &
                            )

         call cubic_spline ( pseudo_chi(1,i_pp_basis_fns), &
                             n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns)), &
                             pseudo_chi_spl(1,1,i_pp_basis_fns) &
                            )

      end do
     


      do i_pp_species = 1, n_pp_species

         if(species_nonlinear_core(pp_species2species(i_pp_species))) then

           
           call cubic_spline ( partial_core_dens(1,i_pp_species), &
                             n_points_pp_fn(i_pp_species), &
                             partial_core_dens_spl(1,1,i_pp_species) &
                            )

           call cubic_spline ( partial_core_dens_deriv(1,i_pp_species), &
                             n_points_pp_fn(i_pp_species), &
                             partial_core_dens_deriv_spl(1,1,i_pp_species) &
                            )


         endif
      end do


      ! divide by two when spin collinear (0.5 for every spin channel)
      if(use_nonlinear_core .and. (n_spin .eq. 2)) then
          partial_core_dens_spl = 0.5*partial_core_dens_spl
          partial_core_dens_deriv_spl = 0.5*partial_core_dens_deriv_spl
      endif




!      if(allocated(partial_core_dens)) then
!         deallocate(partial_core_dens)
!      endif

!      if(allocated(partial_core_dens_deriv)) then
!         deallocate(partial_core_dens_deriv)
!      endif

! --- get outer radius -

      do i_pp_basis_fns = 1, n_pp_basis_fns
         if(pp_basisfn_l(i_pp_basis_fns) .eq. pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
            pp_outer_radius(i_pp_basis_fns) = 1.
            cycle
         endif 
         N = UBOUND(pseudo_wave,1)                                                                           
         do i= 0,N
         if (pseudo_chi(N-i,i_pp_basis_fns).ne. 0.d0) then
             exit
         end if
         end do
         pp_outer_radius(i_pp_basis_fns) = pseudo_grid(N-i + 1,i_pp_basis_fns)
      end do
      deallocate ( local_pseudopot)
      deallocate ( nonlocal_pseudopot)
!      deallocate ( pseudo_wave)      
      deallocate ( pseudo_pot)
      deallocate ( pseudo_fctn_and_pot)  !check whether those arrays are listed in some cleanup routine

    end subroutine spline_and_sort_pseudodata

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!****s* pseudopot/read_pseudopotdata
!  NAME
!    read_pseudopotdata
!  SYNOPSIS

      subroutine read_pseudopotdata(i_species)

      !  PURPOSE
      !
      !  reads in the pseudowavefunction from file *species*.cpi
   
    
      ! USES
      use species_data
      use dimensions

      implicit none

     !  ARGUMENTS
      
      integer,intent(in) :: i_species
      integer :: l

      !  INPUTS
      !    o i_species -- number of species (not pp_species)
      !    o l -- angular momentum of the wanted radial function
      !  OUTPUTS
      !    o n_points_pseudofn -- number of grid points the radial function is defined on
      !    o pseudo_fctn_and_pot(i_pp_species,l,3,n_points_pseudofn) -- (1,n_points_pseudofn) gives the position r on loggrid
      !                                         (2,n_points_pseudofn) gives the radial function u_l(r) on loggrid 
      !                                         (3,n_points_pseudofn) gives the value of potential V_l^{KB}
      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2012).
      !  SOURCE


!     local variables
      integer :: i, i_pp_species
      character*80 :: dummy
      real*8 :: r_dummy
      integer :: i_dummy

      i_pp_species = species2pp_species(i_species)


      open(146,file=trim(pp_path_string(i_species)))   


      do i = 1,11
         read(146,*) dummy
      end do

      ! skip lines in order to allow for correct l-angular function can be read in

      do l = 0,n_pp_fns(i_pp_species)-1, 1

        read(146,*) dummy
        do i = 1, n_points_pp_fn(i_pp_species)  
             read(146,*) i_dummy, pseudo_fctn_and_pot(i_pp_species,l,1,i), &
                         pseudo_fctn_and_pot(i_pp_species,l,2,i), &
                         pseudo_fctn_and_pot(i_pp_species,l,3,i)
        end do

      enddo


      if(species_nonlinear_core(i_species)) then

        do i = 1, n_points_pp_fn(i_pp_species) 
          read(146,*) r_dummy, partial_core_dens(i, i_pp_species), &
                      partial_core_dens_deriv(i, i_pp_species)
        end do


        ! scaling of 0.5 according to Porezag, Pederson, Liu
!        partial_core_dens(:, i_pp_species) = 0.5*partial_core_dens(:, i_pp_species)
!        partial_core_dens_deriv(:, i_pp_species) = 0.5*partial_core_dens_deriv(:, i_pp_species)


      endif

      close(146)

      end subroutine read_pseudopotdata


!******

!-----------------------------------------------------------------------------------
!******
!****s* pseudodata/divide_pseudopot
!  NAME
!   divide_pseudopot
!  SYNOPSIS
     subroutine divide_pseudopot( )
!      PURPOSE
!      divides pseudo_pot into local and non-local part
!
!      USES
       use dimensions
       use geometry

!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

!      counter
       integer :: i_pp_basis_fns, i_pp_species, i, i_pp_atom
       real*8 :: dummy(n_max_points_pp_fn, n_pp_species)

! ----local pseudopot -----------

       do i_pp_basis_fns = 1, n_pp_basis_fns
           i_pp_species = pp_basisfn_species(i_pp_basis_fns)
           if (pp_basisfn_l(i_pp_basis_fns).eq.pp_local_component(i_pp_species)) then
               do i = 1, ubound(pseudo_pot(:,i_pp_basis_fns),1)
               dummy(i,i_pp_species) = pseudo_pot(i,i_pp_basis_fns )
               enddo 
           endif
       enddo

       do i_pp_species = 1, n_pp_species
             local_pseudopot(:,i_pp_species) = dummy(:, i_pp_species)
       enddo


! --- nonlocal pseudpot

       do i_pp_basis_fns = 1, n_pp_basis_fns
           i_pp_species = pp_basisfn_species(i_pp_species)
           if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(i_pp_species)) then
           do i= 1,n_max_points_pp_fn
               nonlocal_pseudopot(i, i_pp_basis_fns) = pseudo_pot(i, i_pp_basis_fns) &
               - local_pseudopot(i,pp_basisfn_species(i_pp_basis_fns)) 
           enddo
           endif
       enddo



     end subroutine divide_pseudopot

!----------------------------------------------------------------------------
!****s* pseudodata/put_localpot_on_integration_grid
!  NAME
!   put_localpot_on_integration_grid
!  SYNOPSIS
     subroutine put_localpot_on_integration_grid( )
!      PURPOSE
!
!      KEEP IN MIND
!      only call this routine after the global integration grid has been set up,
!      so after partition_grid() has been called 
!
!
!      OUTPUT
!      o whole_local_pseudpot_on_intgrid : self-explanatory
!
!      USES
       use dimensions
       use grids
       use geometry
       use spline
 
!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

       real*8, dimension(:,:), allocatable :: points_on_grid

       real*8 :: relvec(3), dist, potential, ilog

       character*130 :: info_str
! counter
       integer :: i_full_points, i_my_batch, i_index
       integer :: i_pp_atom, i_pp_species


       allocate(points_on_grid(3, n_full_points))
  
! first, ensure that arrays are not allocated already (in case of relaxation for example)

       if(allocated(whole_local_pseudpot_on_intgrid)) then 
          deallocate(whole_local_pseudpot_on_intgrid)
       endif 

       allocate( whole_local_pseudpot_on_intgrid(n_full_points))
        whole_local_pseudpot_on_intgrid = 0d0

       i_full_points = 0
       do i_my_batch = 1, n_my_batches, 1
          do i_index = 1, batches(i_my_batch)%size, 1 
             i_full_points = i_full_points + 1
             points_on_grid(:, i_full_points) = batches(i_my_batch)%points(i_index)%coords

             do i_pp_atom = 1, n_pp_atoms
                i_pp_species = pp_species(i_pp_atom)

                relvec = points_on_grid(:, i_full_points) - pp_coords(:, i_pp_atom)
                dist = sqrt(sum(relvec**2))
                ilog = invert_log_grid(dist, &
                &                    pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))


                if(dist.lt.(1.d-8)) then
                   potential = 0.d0
                else if (dist > localpot_outer_radius(i_pp_species)) then
                   potential = -pp_charge(pp_species(i_pp_atom))/dist
                else
                   potential = val_spline ( ilog, local_pseudopot_spl(:,:,i_pp_species), &
                                       & n_points_pp_fn(i_pp_species) )
                endif
       

                whole_local_pseudpot_on_intgrid(i_full_points) = &
                    & whole_local_pseudpot_on_intgrid(i_full_points) + potential


             end do

          end do
       end do

       write(info_str,'(A)') ''
       call localorb_info(info_str)

       write(info_str,'(2X,A)') '*** local pseudopotential part successfully written on global intgrid ***'
       call localorb_info(info_str) 

       write(info_str,'(A)') ''
       call localorb_info(info_str)

    !   whole_local_pseudpot_on_intgrid = 0d0


     deallocate( points_on_grid)

     end subroutine put_localpot_on_integration_grid



!******************************************
!****s* pseudodata/get_pp_norm
!  NAME
!    get_pp_norm
!  SYNOPSIS


  subroutine get_pp_norm()

    !  PURPOSE
    !
    !    Calculates the integral I = int_0^inf dr (u_l*delta V_l)^2.
    !
    !  USES
!    use dimensions
    use dimensions           ! n_points_pp_fn
    use sbt_overlap
    use runtime_choices
    use spline, only: cubic_spline
    implicit none

    !  ARGUMENTS


    !  INPUTS
    !
    !  OUTPUTS
    !    o pp_norm(n_pp_basis_fns) :: integral I
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2012).
    !  SOURCE


!     local variables

    integer :: N
    real*8 :: dr, dlnr, r_i, lnr0, lnrange

    real*8, dimension(:), allocatable :: integrand
    real*8, dimension(:,:), allocatable :: integrand_spl
    real*8, dimension(:), allocatable :: fctn

    ! counter

    integer :: i_pp_basis_fns, i, j


    N = sbtgrid_N      ! N = sbtgrid_N from runtime_choice.f90  
    lnr0 = sbtgrid_lnr0
    lnrange = sbtgrid_lnrange

    dlnr = lnrange / N

    pp_norm(:) = 0d0



    do i_pp_basis_fns = 1, n_pp_basis_fns

        allocate(fctn(N))
        if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
        allocate(integrand(n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns))))
        allocate(integrand_spl(4,n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns))))
        integrand(:) = (pseudo_wave(:, i_pp_basis_fns)*nonlocal_pseudopot(:,i_pp_basis_fns))**2 &
                        *pseudo_grid(:,i_pp_basis_fns)

! pseudo_grid is multiplied here because sbt_import_spline scales with 1/r  

        call cubic_spline ( integrand(1), &
                             n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns)), &
                             integrand_spl(1,1) &
                            )

        call sbt_import_spline(N, fctn, lnr0, lnrange, &
                               pp_basisfn_l(i_pp_basis_fns),n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns)), &
                               integrand_spl(:,:), pseudo_grid(1,i_pp_basis_fns), &
                               (pseudo_grid(2,i_pp_basis_fns)/ pseudo_grid(1,i_pp_basis_fns)))

        do i = 1,N
           dr = exp(lnr0 + i*dlnr) - exp(lnr0 + (i-1)*dlnr)
           pp_norm(i_pp_basis_fns) = pp_norm(i_pp_basis_fns) + fctn(i)*dr

        end do


        deallocate(integrand)
        deallocate(integrand_spl)


        else
        pp_norm(i_pp_basis_fns) = 0d0

        end if
        deallocate(fctn)


     end do 





  end subroutine get_pp_norm
!******************************************
!******************************************
!****s* pseudodata/check_partial_core_norm
!  NAME
!    check_partial_core_norm
!  SYNOPSIS


  subroutine check_partial_core_norm()

    !  PURPOSE
    !
    !  USES
    use dimensions           ! n_points_pp_fn
    use species_data

    implicit none


    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


!     local variables

    integer :: N
    real*8 :: dr
    real*8 :: partial_core_norm

    character*130 :: info_str

    real*8, dimension(:), allocatable :: integrand

    ! counter

    integer :: i_pp_species, i

    partial_core_norm = 0d0



      allocate(integrand(n_points_pp_fn(1)))

      integrand(:) = partial_core_dens(:,1)

      do i = 2,n_points_pp_fn(1)-1

         dr = (pseudo_grid(i+1,1) - pseudo_grid(i-1,1))/2
         partial_core_norm = partial_core_norm + integrand(i)*dr*(pseudo_grid(i,1))**2

      end do

      deallocate(integrand)

      write(info_str,'(A,A,F6.3)') '  Appling nonlinear core correction. ',&
                                  'Partial core density:', partial_core_norm
      call localorb_info(info_str) 

      write(info_str,*)
      call localorb_info(info_str)

  end subroutine check_partial_core_norm
!******************************************

!******************************************
!****s* pseudodata/get_nonlocal_pot
!  NAME
!    get_nonlocal_pot
!  SYNOPSIS


  subroutine get_nonlocal_pot()

    !  PURPOSE
    !
    !    Calculates the coefficients E_l^{KB} of the nonlocal potential 
    !    and concatenate with basiswave_pp_overlap to nonlocal_matrix
    !     
    !  USES
    use dimensions           ! n_points_pp_fn
    use sbt_overlap
    use runtime_choices
    use spline, only: cubic_spline
    implicit none


!DB: 02/16/12 : E_l^{KB} calculated here perfectly match those from pswatch
!DB: 06/04/12 : We found an inconsistance between the documentaion (M.Fuchs paper)
!               and the actual output of FHI98PP.
!               After a complete check of FHI98PP and this pseudopot infrastructure
!               here, we are sure about that.

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2012).
    !  SOURCE


!     local variables
!    real*8, dimension(n_pp_basis_fns) :: E_l_KB
    real*8, dimension(:), allocatable :: integrand
    real*8, dimension(:,:), allocatable :: integrand_spl
    real*8, dimension(:), allocatable :: fctn
    integer :: N
    real*8 :: dr, dlnr, r_i, lnr0, lnrange, dummy

    ! counter

    integer :: i_pp_basis_fns, i_basis_1, i_basis_2, i_pp_basis, i, k, l


    N = sbtgrid_N      ! N = sbtgrid_N from runtime_choice.f90  
    lnr0 = sbtgrid_lnr0
    lnrange = sbtgrid_lnrange

    dlnr = lnrange / N

    pp_psi_chi_overlap(:) = 0d0

    do i_pp_basis_fns = 1, n_pp_basis_fns


        allocate(fctn(N))
        if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
        allocate(integrand(n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns))))
        allocate(integrand_spl(4,n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns))))

!        integrand(:) = pseudo_chi(:,i_pp_basis_fns)*pseudo_wave(:,i_pp_basis_fns) &
!                      *pseudo_grid(:,i_pp_basis_fns)

!! dd/mm/yy 29/01/14: to make it more understandable we do the necessary scaling with pp_norm here.
!                     Units are messed up (this is the bespoken ambiguity in Fuchs' paper)
!                     (according to ppcheck.f in FHi98PP we still have multiply with sqrt(pp_norm))
        integrand(:) = pseudo_chi(:,i_pp_basis_fns)*pseudo_wave(:,i_pp_basis_fns) &
                      *pseudo_grid(:,i_pp_basis_fns)*(pp_norm(i_pp_basis_fns))**(1./2.)


! this is already scaled with r to compensate scaling of sbt_import_spline

         call cubic_spline ( integrand(1), &
                             n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns)), &
                             integrand_spl(1,1) &
                            )

        call sbt_import_spline(N, fctn, lnr0, lnrange, &
                               pp_basisfn_l(i_pp_basis_fns),n_points_pp_fn(pp_basisfn_species(i_pp_basis_fns)), &
                               integrand_spl(:,:), pseudo_grid(1,i_pp_basis_fns), &
                               (pseudo_grid(2,i_pp_basis_fns)/ pseudo_grid(1,i_pp_basis_fns)))

        do i = 1,N
           dr = exp(lnr0 + i*dlnr) - exp(lnr0 + (i-1)*dlnr)
           r_i = exp(lnr0) +dr
           pp_psi_chi_overlap(i_pp_basis_fns) = pp_psi_chi_overlap(i_pp_basis_fns) + fctn(i)*dr

        end do

        deallocate(integrand)
        deallocate(integrand_spl)


        else
        pp_psi_chi_overlap(i_pp_basis_fns) = 0d0

        end if
        deallocate(fctn)


    end do 


    do i_pp_basis_fns =1, n_pp_basis_fns
!   only if the angular momentum of i_pp_basis_fns does not correspond to the 
!   local component of the species i_pp_basis_fns referes to
        if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
        E_l_KB(i_pp_basis_fns) = pp_norm(i_pp_basis_fns)/pp_psi_chi_overlap(i_pp_basis_fns)




        else
        E_l_KB(i_pp_basis_fns) = 0d0
        end if


    end do

! all these numbers match those from ppcheck.f in FHI98PP
!        if(myid==0) then
!        do i = 1, n_pp_basis_fns
!          write(use_unit,*)  'angular momentum :  ', i-1
!          write(use_unit,*)  'pp_psi_chi_overlap:  ', pp_psi_chi_overlap(i)
!          write(use_unit,*)  'E_l_KB:  ',E_l_KB(i)
!          write(use_unit,*)  'pp_norm:  ',pp_norm(i)
!        enddo
!        endif




!TODO:DB 10.11.: for huge systems saving nonlocal_matrix can become a bottle neck .. 
!                we should think about something to solve that, one day


     nonlocal_matrix = 0.d0

     do i_pp_basis = 1, n_pp_basis
        i_pp_basis_fns = pp_basis_fn(i_pp_basis)
        dummy = 0
        if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then

     do k = 1, n_max_basis_overlap(i_pp_basis)

           i_basis_1 = nonzero_overlap_entries(k,i_pp_basis)

        do l = 1, n_max_basis_overlap(i_pp_basis)

           i_basis_2 = nonzero_overlap_entries(l,i_pp_basis)

!     do i_basis_1 = 1, n_basis
!     do i_basis_2 = 1, n_basis

!     do i_pp_basis = 1, n_pp_basis
!        i_pp_basis_fns = pp_basis_fn(i_pp_basis)
!        dummy = 0
!        if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then

!DB changes on 10/04/13

         dummy =  basiswave_pp_overlap(i_basis_1,i_pp_basis)*E_l_KB(i_pp_basis_fns)*basiswave_pp_overlap(i_basis_2,i_pp_basis)
!   although this follows exactly M.Fuchs instructions, this is actually wrong. One can verify that by comparing everything step 
!   by step with the intermediate and final results of "pswatch" (FHI98PP)

!        else 
!         dummy = 0.d0
!        endif

        nonlocal_matrix(i_basis_1,i_basis_2) = nonlocal_matrix(i_basis_1,i_basis_2) + dummy  


     end do 
     end do
     endif
     end do


!     nonlocal_matrix = 0.d0

  end subroutine get_nonlocal_pot


!******************************************
!****s* pseudodata/get_pp_nonlocal_force_matrix
!  NAME
!    get_pp_nonlocal_force_matrix
!  SYNOPSIS
!

  subroutine  get_pp_nonlocal_force_matrix()

    !  PURPOSE
    !
    !    Calculates the nonlocal force matrix
    !    d_ 
    !    and concatenate with basiswave_pp_overlap to nonlocal_matrix
    !     
    !  USES
    use dimensions           ! n_points_pp_fn
    use runtime_choices

    implicit none

    integer :: k, l, i_pp_basis, i_pp_basis_fns, i_coord, i_pp_atom

    if(.not.allocated(pp_nonlocal_forces_matrix)) then
    allocate(pp_nonlocal_forces_matrix(3,n_basis,n_basis,n_pp_atoms))
    endif
    pp_nonlocal_forces_matrix = 0.d0


!     do i_pp_atom = 1, n_pp_atoms
     do k = 1, n_basis
        do l = 1, n_basis
        do i_pp_basis = 1, n_pp_basis
           i_pp_basis_fns = pp_basis_fn(i_pp_basis)
           if (pp_basisfn_l(i_pp_basis_fns).ne.pp_local_component(pp_basisfn_species(i_pp_basis_fns))) then
           i_pp_atom = pp_basis_atom(i_pp_basis)

              do i_coord = 1,3,1

!              pp_nonlocal_forces_matrix(i_coord,k,l,i_pp_atom) = pp_nonlocal_forces_matrix(i_coord,k,l,i_pp_atom)  + & 
!                   (-d_basiswave_pp_overlap(i_coord,k, i_pp_basis)/pp_psi_chi_overlap(i_pp_basis_fns)*&
!                   pp_basiswave_overlap(i_pp_basis,l) + &
!                   basiswave_pp_overlap(k, i_pp_basis)/pp_psi_chi_overlap(i_pp_basis_fns)*&
!                   d_pp_basiswave_overlap(i_coord,i_pp_basis,l))

! 29.1.14 wieder konsistent mit der energie machen
              pp_nonlocal_forces_matrix(i_coord,k,l,i_pp_atom) = pp_nonlocal_forces_matrix(i_coord,k,l,i_pp_atom)  + & 
                   (-d_basiswave_pp_overlap(i_coord,k, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
                   pp_basiswave_overlap(i_pp_basis,l) + &
                   basiswave_pp_overlap(k, i_pp_basis)*E_l_KB(i_pp_basis_fns)*&
                   d_pp_basiswave_overlap(i_coord,i_pp_basis,l))

              enddo

           endif
        enddo
        enddo
     enddo
!     enddo
  end subroutine  get_pp_nonlocal_force_matrix


!******************************************
!****s* pseudodata/get_pp_distances2atom
!  NAME
!    get_pp_distance2atom
!  SYNOPSIS


  subroutine get_pp_distance2atom()

    !  PURPOSE
    !
    !    Measures the distance from the QM region
    !     and writing result to logical array
    !    We need this array for setting up the extended integration grid
    !     for the qm-embedding. Same needs to be done for multipoles 
    !
    !    NOTE: only call this routine after get_grids have been called
    !      
    !     
    !  USES
    use dimensions
    use geometry          
    use grids
    use basis, only: outer_radius
  
    implicit none


    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2012).
    !  SOURCE


!     local variables

    real*8 :: distance, min_distance
    real*8 :: radius, qm_border
    integer :: closest_atom    

   ! counter

    integer :: i_pp_atom, i_atom

     n_pp_in_qm = 0

     closest_atom = 0

     qm_border = max(maxval(outer_radius),maxval(r_radial))

     do i_pp_atom = 1, n_pp_atoms
        min_distance = 300000.    
	do i_atom = 1, n_real_atoms

	   distance = (coords(1,i_atom) - pp_coords(1,i_pp_atom))**2 + &
                 & (coords(2,i_atom) - pp_coords(2,i_pp_atom))**2 + &
                 & (coords(3,i_atom) - pp_coords(3,i_pp_atom))**2

           if (distance.lt.min_distance) then
	       closest_atom = i_atom
               min_distance = distance
           endif

	enddo 

! DB: 04/26/12 definition if radius is to be optimized
!           so far we are on the safe but also slow side
        radius = maxval(pp_outer_radius) + max(maxval(outer_radius),maxval(r_radial))
        radius = radius**2
        
	if (min_distance.lt.radius) then
            pp_in_qm(i_pp_atom) = .true.
            n_pp_in_qm = 1 + n_pp_in_qm
	else
            pp_in_qm(i_pp_atom) = .false.                 
	endif
     enddo

  

  end subroutine get_pp_distance2atom
!----------------------------------------------------------------------------
!******************************************************
!****s* pseudodata/cleanup_pseudodata
!  NAME
!   cleanup_pseudodata
!  SYNOPSIS

!
!   This routine deallocates all pseudodata arrays except 
!       o whole_local_pseudopot_on_intgrid
!       o nonlocal_matrix 
!
!   This routine should only be called after those two arrays
!      have been computed, so after put_localpot_on_integration_grid 
!      has been called --> right before entering the scf_cycle

     subroutine cleanup_pseudodata( )
!      PURPOSE
!
!      USES
       use dimensions
       use geometry

!  AUTHOR
!      FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SOURCE
       implicit none

!      counter


! data arrays in dimensions.f90
       deallocate(n_pp_fns)
       deallocate(n_points_pp_fn)
       deallocate(pp_function_index)

!! data arrays in species_data.f90
!       deallocate(pseudo_fctn_and_pot)

! data arrays in geometry.f90
       if (allocated(pp_coords)) then
          deallocate(pp_coords)
       end if
       if (allocated(pp_species)) then
          deallocate(pp_species)
       end if
! data arrays in pseudodata.f90
       if (allocated(pseudo_wave)) then
          deallocate(pseudo_wave)
       end if
       if (allocated(pseudo_grid)) then
          deallocate(pseudo_grid)
       end if
       if (allocated(pseudo_chi)) then
          deallocate(pseudo_chi)
       end if
       
       if (allocated(pp_basisfn_l)) then
          deallocate(pp_basisfn_l)
       end if
       if (allocated(pp_basisfn_species)) then
          deallocate(pp_basisfn_species)
       end if
       if (allocated(pp_outer_radius)) then
          deallocate(pp_outer_radius)
       end if

       if (allocated(pp_basis_atom)) then
          deallocate(pp_basis_atom)
       end if 
       if (allocated(pp_basis_species)) then
          deallocate(pp_basis_species)
       end if
       if (allocated(pp_basis_l)) then
          deallocate(pp_basis_l)
       end if 
       if (allocated(pp_basis_m)) then
          deallocate(pp_basis_m)
       end if
       if (allocated(pp_basis_fn)) then
          deallocate(pp_basis_fn)
       end if
       if (allocated(pp_r_grid_min)) then
          deallocate(pp_r_grid_min)
       end if 
       !deallocate(pp_r_grid_inc)

       !deallocate(pseudo_wave_spl)
       !deallocate(pseudo_chi_spl)
       !deallocate(local_pseudopot_spl)

       !deallocate(pp_psi_chi_overlap)
       !deallocate(basiswave_pp_overlap)
       !deallocate(E_l_KB)

       if(allocated(d_basiswave_pp_overlap)) then
         deallocate(d_basiswave_pp_overlap)
       endif
       if(allocated(pp_nonlocal_forces)) then
         deallocate(pp_nonlocal_forces)
       end if
       if(allocated(pp_nonlocal_forces_2)) then
         deallocate(pp_nonlocal_forces_2)
       end if
       if(allocated(pp_nlcc_forces)) then
         deallocate(pp_nlcc_forces)
       end if
       if(allocated(pp_nonlocal_forces_matrix)) then
         deallocate(pp_nonlocal_forces_matrix)
       end if
       
!!! Arrays below are deallocated in spline_and_sort_pseudodata as well
       if (allocated(nonlocal_pseudopot)) then
          deallocate(nonlocal_pseudopot)
       end if
       if (allocated(local_pseudopot)) then 
           deallocate(local_pseudopot)
       end if
       if (allocated(pseudo_pot)) then 
           deallocate(pseudo_pot)
       end if
       if (allocated(pseudo_fctn_and_pot)) then
          deallocate(pseudo_fctn_and_pot)
       end if
!!!          
       if (allocated(pp_norm)) then 
          deallocate(pp_norm)
       end if
       if (allocated(nonlocal_matrix)) then
          deallocate(nonlocal_matrix)
       end if
       if (allocated(localpot_outer_radius)) then
          deallocate(localpot_outer_radius)
       end if
       if (allocated(Lsp2n_pp_basis_fnLsp)) then
          deallocate(Lsp2n_pp_basis_fnLsp)
       end if
       if (allocated(Lsp2_pp_basis_fn)) then
          deallocate(Lsp2_pp_basis_fn)
       end if
       if (allocated(partial_core_dens_spl)) then
             deallocate ( partial_core_dens_spl)
       end if

       if (allocated(partial_core_dens_deriv_spl)) then
             deallocate ( partial_core_dens_deriv_spl)
       end if
!data arrays in allocate_pseudospecies_arrays

      deallocate(species_nonlinear_core)
      deallocate(species2pp_species)
      deallocate(pp_species2species)
      deallocate(pp_path_string)
      deallocate(pp_charge)
      deallocate(pp_local_component)

      if (allocated(pp_r_grid_min)) then
         deallocate(pp_r_grid_min)
      end if
      if (allocated(pp_r_grid_inc)) then
         deallocate(pp_r_grid_inc)
      end if
      if(allocated(partial_core_dens)) then
         deallocate(partial_core_dens)
      end if
      
      if(allocated(partial_core_dens_deriv)) then
         deallocate(partial_core_dens_deriv)
      end if

     end subroutine cleanup_pseudodata


end module pseudodata
