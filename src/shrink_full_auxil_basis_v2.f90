!****s* FHI-aims/shrink_full_auxil_basis_v2
!  NAME
!   shrink_full_auxil_basis_v2
!  SYNOPSIS

 subroutine shrink_full_auxil_basis_v2 &
        ( )

!  PURPOSE
!   Collects all the possible product functions below some maximal angular
!   momentum max_l_prodbas.
!   orthogonalizes extra basis functions to atomic ones, for each species.
!   throws out all basis functions which are substantially identical to another one
!   gathers basis functions in a compressed array for numerics
!   implicitly sort all basis functions into blocks of identical l per species
!   labels all basis functions to remember their origin

!  USES
      use dimensions
      use runtime_choices
      use localorb_io
      use grids
      use geometry
      use basis
      use species_data
      use spline
      use prodbas
      use hartree_fock
      use mpi_tasks
      use constants
      use pbc_lists
      use synchronize_mpi
      use basbas_fn_coulomb
      implicit none

!  INPUTS
!    none
!  OUTPUTS
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


!  FIXME: For optimum integrations, the present subroutine should also interpolate all
!         wave functions as splines onto the radial integration grid, and _then_ perform
!         all normalisation integrals on this grid. THEN, we would have a numerically
!         well-defined Hamiltonian matrix later on.
!         The present case will basbasuce exact orthonormality on the log grid, but only
!         approximate orthonormality on the later radial integration grid.



!   max_n_prodbas is the highest n index of the orginal basis used to construct the
!   product basis

!     output

!     all output is found in module basis

!  local variables

!     i_first_fn : First basis function for current atom and current angular momentum shell

      integer i_first_fn
      integer i_first_species_fn
      real*8 norm

      integer n_reject

      integer fn_l_high
      integer fn_l_low

!     Arrays for intermediate storage/handling of radial functions before spline
!     basis_wave is compacted array of basis functions
!     basis_deriv is compacted array of basis function derivatives (if needed)
!     basis_kinetic : directly calculated kinetic energy term for all basis functions
!                which are not free-atom like basis functions. [Treat free-atom like basis functions
!                separately so we can subtract out the free-atom basis-defining potential later.]

      real*8, dimension(:,:), allocatable :: basbas_wave
      integer, dimension(:), allocatable :: species_basbas_fns
      integer, dimension(:), allocatable :: fn_l_basbas

!  counters

      integer i_species, i_species_1
      integer i_basis, i_function, i_atom, i_atom_1, i_atom_ind, i_shell, i_index
      integer i_m, i_l, i_grid, i_type
      integer i_n, i_n_1,i_l_1
      integer i_basbas_fn
      integer i_prev
!SVL
      integer i_basbas, i_basis_1, i_task, i_basbas_1, i_cell_1, i_cell_2, i_cell_3, i_cell, i_cell_tmp
      integer :: i_x, i_y, i_z, i_basbas_2, i_periodic
      real*8 :: coords_cell(3), outer_radius_pb, dist, rad_sum
      logical :: found

!  functions

      real*8 int_log_mesh
      real*8 int_log_coulomb_metric
      real*8 get_inner_max

!  debug

      integer species_atom(n_species)
      integer n_func(n_species)
      integer fn_species(n_max_basis_fns*n_max_basis_fns*10)
      character*8 fn_type(n_max_basis_fns)
      integer fn_l(n_max_basis_fns,n_species)
      integer fn_n(n_max_basis_fns,n_species)
      integer fn_index(n_max_basis_fns,n_species)

      logical flag_species_incl
      integer prodbas_metric
      logical :: use_hse_localization

      integer success

      ! For output using localorb_io
      character*200 :: info_str

      character*16 output_name
      real*8 scalar_basbas(n_max_basis_fns)

      character(*), parameter :: func = 'shrink_full_auxil_basis_v2'

!  begin work

      call localorb_info (' ',use_unit,'(2X,A)', OL_norm  )
      call localorb_info ('--------------------------------------------',use_unit,'(2X,A)', OL_norm  )

      write(info_str,'(A)') "Constructing auxiliary basis (full product) ..."
      call localorb_info (info_str,use_unit,'(2X,A)', OL_high  )
      call localorb_info (' ',use_unit,'(2X,A)', OL_norm  )
      use_hse_localization = (use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse &
          .and. .not. use_dftpt2_and_hse)

      i_task = myid + 1

!     allocations - local variables for basis function storage
      allocate ( basbas_wave (n_max_grid, &
                 n_max_basis_fns*10) )
      basbas_wave = 0.d0

      allocate ( species_basbas_fns ( &
                 n_max_basis_fns*10) )
      allocate ( fn_l_basbas ( &
                 n_max_basis_fns*10) )

!     determine basis _functions_ first

      species_atom = 0
      n_func = 0
      fn_index = 0
      i_function = 0

! how to orthonormalize the auxiliary basis: = 1, Coulomb metric, =0, Normal metric
      if (RI_type == RI_SVS) then
         prodbas_metric = 0   ! Normal metric 
      else
         prodbas_metric = 1   ! Coulomb metric
      endif

!  extract the information about the radial function from the original basis
      do i_species = 1, n_species
        flag_species_incl = .false.
        do i_basis = 1, n_ext

          i_atom = ext_atom (i_basis)
          i_l    = ext_l (i_basis)
          i_m    = ext_m (i_basis)
          i_n    = extfn_n (ext_fn(i_basis))
          i_species_1 = species(i_atom)

          if(i_species_1.eq.i_species) then

            flag_species_incl = .true.

            if (species_atom(i_species).eq.0) then
               species_atom (i_species) = i_atom
            endif

!   only count the first atom belonging to this particular species
            if (i_atom.eq.species_atom(i_species) .and. &
                 i_m .eq. -i_l)  then

               n_func(i_species) = n_func(i_species) + 1
               fn_l(n_func(i_species),i_species) = i_l
               fn_n(n_func(i_species),i_species) = i_n

               i_function = i_function +1
               fn_index (n_func(i_species),i_species) = i_function

            endif
           endif
           
        enddo
        if(.not.flag_species_incl) then
            i_function = i_function + n_ext_fn_species(i_species)
        endif
      enddo

        i_basbas_fn = 0
        n_basbas_supercell = 0
        n_basbas = 0
        n_reject = 0
        do i_species = 1, n_species, 1

           i_first_species_fn = i_basbas_fn + 1
           i_basis =0

          if( max_l_prodbas(i_species).gt.2*ext_l_shell_max(i_species) ) &
                        then
              max_l_prodbas(i_species)= 2*ext_l_shell_max(i_species)
          endif

          do i_l = 0, max_l_prodbas(i_species), 1
            i_first_fn = i_basbas_fn + 1

! include the original "single" basis
            if (flag_basis.eq.FLAG_BASIS_INCL_ONEPARTICLE) then
              do i_function = 1, n_func(i_species), 1

                 i_l_1 = fn_l (i_function, i_species)

                 if(i_l_1 .eq. i_l) then
                   i_basbas_fn = i_basbas_fn +1
                   species_basbas_fns (i_basbas_fn) = i_species
                   fn_l_basbas(i_basbas_fn) = i_l
                   do i_grid = 1, n_grid(i_species)
                       basbas_wave(i_grid,i_basbas_fn) &
                        = basis_wave_spl( 1, i_grid, &
                          fn_index(i_function, i_species) )
                   enddo
                   fn_species(i_basbas_fn)=i_species
                   i_basis = i_basis + 2*i_l_1 +1

                 endif
                 do i_atom_1 = 1, n_atoms
                    if(species(i_atom_1).eq.i_species)then


                       do i_cell = 1, n_cells_pairs
                          found = .false.
                          coords_cell(:) = coords_center(:,i_atom_1)
                          do i_periodic = 1, n_periodic
                             coords_cell(:) = coords_cell(:)+&
                                  cell_index_pairs(i_cell,i_periodic)*&
                                  lattice_vector(:,i_periodic)
                          enddo
                          do i_atom = 1, n_atoms
                             rad_sum = outer_radius(i_function)+&
                                  atom_radius(species_center(i_atom))
                             dist = (coords_cell(1) - coords_center(1,i_atom)) * &
                                  (coords_cell(1) - coords_center(1,i_atom)) + &
                                  (coords_cell(2) - coords_center(2,i_atom)) * &
                                  (coords_cell(2) - coords_center(2,i_atom)) + &
                                  (coords_cell(3) - coords_center(3,i_atom)) * &
                                  (coords_cell(3) - coords_center(3,i_atom))
                             if(sqrt(dist).lt.rad_sum)then
                                found = .true.
                                exit
                             endif
                          enddo
                          if(found)&
                               n_basbas_supercell = n_basbas_supercell + 2*i_l_1 +1
                               
                       enddo
                    endif
                 enddo
              enddo
            endif


!   construct product basis
            do i_function = 1, n_func(i_species), 1
             i_n = fn_n(i_function, i_species)
             do i_prev = 1, i_function, 1
              i_n_1 = fn_n(i_prev, i_species)

!              if ( i_n .le. max_n_prodbas(i_species) .and. &
!                   i_n_1 .le. max_n_prodbas(i_species) ) then

                fn_l_high = &
                  fn_l(i_function, i_species) + &
                  fn_l(i_prev, i_species)

                fn_l_low = &
                 abs( fn_l(i_function, i_species) - &
                  fn_l(i_prev, i_species))

                if (i_l.ge.fn_l_low .and. &
                          i_l.le.fn_l_high) then

                   i_basbas_fn =  i_basbas_fn +1
                   species_basbas_fns (i_basbas_fn) = i_species
                   fn_l_basbas(i_basbas_fn) = i_l

                   do i_grid = 1, n_grid(i_species)
                       basbas_wave(i_grid,i_basbas_fn) &
                        = ext_wave_spl( 1, i_grid, &
                          fn_index(i_function, i_species) ) * &
                          ext_wave_spl( 1, i_grid, &
                          fn_index(i_prev, i_species) ) &
                          /r_grid(i_grid, i_species)
                   enddo

!      first normalize the individual wave function
                   if(prodbas_metric.eq.0) then
                      norm = int_log_mesh ( &
                           basbas_wave(1, i_basbas_fn), &
                           basbas_wave(1, i_basbas_fn), &
                           n_grid (i_species), &
                           r_grid (1, i_species) )
                   elseif (prodbas_metric.eq.1) then
                      norm = int_log_coulomb_metric ( i_l, &
                           basbas_wave(1, i_basbas_fn), &
                           basbas_wave(1, i_basbas_fn), &
                           n_grid (i_species), &
                           r_grid (1, i_species) )
                   else
                      write(use_unit,'(2X,2A)') &
                         "* Error: No such metric for on-site auxiliary",&
                          " basis orthormalization exists."

                    stop
                 endif

                   do i_grid = 1, n_grid(i_species)
                    basbas_wave (i_grid,i_basbas_fn) = &
                    basbas_wave (i_grid,i_basbas_fn)/sqrt(norm)
                   enddo

!     then orthogonalize it to the previous ones
                   call orthonormalize_prodbas_fn &
                        ( prodbas_metric, i_l, i_basbas_fn, i_first_fn, &
                        n_max_basis_fns*10, basbas_wave, n_grid(i_species), &
                        r_grid(1,i_species), prodbas_acc(i_species), &
                        success &
                        )

                if (success .eq. 1) then

!                 if(myid.eq.0) then
!                  write(use_unit,'(2X,A,1X,I2,2X,A,1X,I2,2X,A,2I4,2X,A)' ) &
!                   "| Species", i_species, "l_shell",i_l , &
!                     ":   pair state", &
!                    i_function, i_prev, "   accepted."
!                 endif

                  i_basis = i_basis + 2*i_l + 1
                 do i_atom_1 = 1, n_atoms
                    if(species(i_atom_1).eq.i_species)then
                       do i_cell = 1, n_cells_pairs
                          found = .false.
                          coords_cell(:) = coords_center(:,i_atom_1)
                          do i_periodic = 1, n_periodic
                             coords_cell(:) = coords_cell(:)+&
                                  cell_index_pairs(i_cell,i_periodic)*&
                                  lattice_vector(:,i_periodic)
                          enddo
                          do i_atom = 1, n_atoms
                             rad_sum = atom_radius(i_species)+&
                                  atom_radius(species_center(i_atom))
                             dist = (coords_cell(1) - coords_center(1,i_atom)) * &
                                  (coords_cell(1) - coords_center(1,i_atom)) + &
                                  (coords_cell(2) - coords_center(2,i_atom)) * &
                                  (coords_cell(2) - coords_center(2,i_atom)) + &
                                  (coords_cell(3) - coords_center(3,i_atom)) * &
                                  (coords_cell(3) - coords_center(3,i_atom))
                             if(sqrt(dist).lt.rad_sum)then
                                found = .true.
                                exit
                             endif
                          enddo
                          if(found)&
                               n_basbas_supercell = n_basbas_supercell + 2*i_l +1

                       enddo
                    endif
                 enddo

                elseif (success .eq. -1) then
!                 if(myid.eq.0) then
!                  write(use_unit,'(2X,A,1X,I2,2X,A,1X,I2,2X,A,2I4,2X,A)' ) &
!                  "| Species", i_species, "l_shell", i_l, &
!                     ":   pair state", &
!                   i_function, i_prev, "  enough product function"
!                  write(use_unit,'(2X, A)') &
!                    "already for this l channel, stop here ."
!                 endif

                  n_reject = n_reject +1
                else
!                 if(myid.eq.0) then
!                  write(use_unit,'(2X,A,1X,I2,2X,A,1X,I2,2X, A,2I4,2X,A)' ) &
!                  "| Species", i_species, "l_shell",i_l, &
!                      ":   pair state", i_function, &
!                   i_prev, "   rejected: Linear dependence."
!                 endif

                   i_basbas_fn = i_basbas_fn - 1

                   n_reject = n_reject +1
                endif

                 fn_species(i_basbas_fn)=i_species
!  end of if i_l
            endif

!      end of if i_n
!            endif
!      end of do i_prev
            enddo
!      end of do i_function
           enddo

!   Localize basbas_fn within one l-channel
           call localize_basbas_fn_coulomb(i_l, i_first_fn, i_basbas_fn, &
           &      basbas_wave, n_grid(i_species), r_grid(:, i_species), &
           &      Adams_Moulton_integrator, use_hse_localization)

!          write(use_unit,*) "i_l", i_l , l_shell_max(i_species),
!     +          n_max_basis_fns
!       end of do i_l
          enddo

          ! Sort basbas_fn within one species
          if (.not. (use_hse .and. hse_omega_hf /= 0.d0).or.use_gw_and_hse &
              .or. use_dftpt2_and_hse) then
             call sort_species_basbas_fn(i_basbas_fn-i_first_species_fn+1, &
             &       n_grid(i_species), r_grid(:, i_species), wave_threshold, &
             &       basbas_wave(:, i_first_species_fn:i_basbas_fn), &
             &       fn_l_basbas(i_first_species_fn:i_basbas_fn))
          end if
          ! No need to update species_basbas_fns as this is done only within
          ! one species.


        n_basbas = n_basbas + i_basis * &
             atoms_in_structure (i_species)


!       end of do i_species
      enddo

      n_basbas_fns = i_basbas_fn

      call get_bas_dimensions(n_basbas_fns, fn_l_basbas(1:n_basbas_fns), &
      &                       species_basbas_fns(1:n_basbas_fns), &
      &                       max_basbas_L, max_n_basbas_fnLsp, n_basbas)

      call allocate_basbas()

      basbasfn_l = fn_l_basbas(1:n_basbas_fns)
      basbasfn_species = species_basbas_fns(1:n_basbas_fns)

      call generate_basbasfn(basbas_wave, basbasfn_l, basbasfn_species)

      call generate_Lsp_indexing(n_basbas_fns, basbasfn_l, basbasfn_species, &
      &                          max_basbas_L, max_n_basbas_fnLsp, n_basbas, &
      &                          Lsp2n_basbas_fnLsp, Lsp2basbas_fn, Lsp2basbas_sp, &
      &                          atom2basbas_off, sp2n_basbas_sp)
      max_n_basbas_sp = maxval(sp2n_basbas_sp)

      if(myid.eq.0) then
       write(use_unit, '(2X, A,1X, A, I8, 1X, A, I8, 1X, A)') &
        "| Shrink_full_auxil_basis :", "there are ", n_basbas_fns, &
          " radial auxiliary wave functions"
        write(use_unit, '(27X,A,I8,A)') &
         " accepted and", n_reject, " rejected."
       write(use_unit, '(2X,A,I20,A)') &
          "| Shrink_full_auxil_basis : there are totally", n_basbas_supercell, &
          " partial auxiliary wave functions."
      endif

! SVL For periodic systems, include auxiliary basis functions on 
!     all atoms, including atoms that touch the original unit cell

      i_basis = 0
      do i_atom_ind = 1, n_centers_ele_summation !SVL modified natoms -> n_centers_ele_summation
         i_atom = centers_ele_summation (i_atom_ind)
          i_species = species_center (i_atom) !SVL modified
         do i_basbas_fn = 1, n_basbas_fns
            i_species_1 = basbasfn_species ( i_basbas_fn)

            if(i_species .eq. i_species_1 ) then
               found = .false.
               do i_atom_1 = 1, n_atoms
                  rad_sum = atom_radius(i_species)+&
                       atom_radius(species_center(i_atom_1))
                  dist = (coords_center(1,i_atom) - coords_center(1,i_atom_1)) * &
                       (coords_center(1,i_atom) - coords_center(1,i_atom_1)) + &
                       (coords_center(2,i_atom) - coords_center(2,i_atom_1)) * &
                       (coords_center(2,i_atom) - coords_center(2,i_atom_1)) + &
                       (coords_center(3,i_atom) - coords_center(3,i_atom_1)) * &
                       (coords_center(3,i_atom) - coords_center(3,i_atom_1))
                  if(sqrt(dist).lt.rad_sum)then
                     found = .true.
                     exit
                  endif
               enddo
               if(found)then
                  i_l = basbasfn_l (i_basbas_fn)
                  
                  do i_m = -i_l, i_l, 1
                     
                     i_basis = i_basis +1
                     basbas_atom ( i_basis ) = i_atom
                     basbas_l ( i_basis ) = i_l
                     basbas_m ( i_basis ) = i_m
                     basbas_fn( i_basis ) = i_basbas_fn
!  end of do i_m
                  enddo
               endif
!  end of if i_species
            endif
!  end of do i_basbas_fn
           enddo
!  end of do i_atom
       enddo

       if (i_basis .ne. n_basbas_supercell) then

         write(use_unit, '(2X, A, A, 1X, A, I5, 1X, A, I5)' ) &
                     "Shrink_basbas_basis :", &
                     " the product basis number is not consistent. ", &
                     " there should be ", n_basbas_supercell,  "but only ", &
                      i_basis
       endif


!      Product basis distributed over different threads for periodic boundary conditions.
!      For the non-periodic case, this step may nor strictly be needed, done only to keep
!      a single line of code later.
!
!      Only call this if the actual distribution did change from the preceding step.
       if ( (.not.prodbas_tasks_pbc_distributed) .or. & 
            (previous_n_basbas_supercell.ne.n_basbas_supercell) ) then

            write(info_str,'(A)') & 
               "prodbas_tasks_distribution_pbc not up to date - distributing tasks."
            call localorb_info (info_str,use_unit,'(2X,A)', OL_high  )

            call prodbas_tasks_distribution_pbc ()

            prodbas_tasks_pbc_distributed = .true.
            previous_n_basbas_supercell = n_basbas_supercell
       end if

       ! Also distribute separately the prod basis from the first unit cell
       ! for parallel overlap and coulomb matrix calculations
       ! Only do this if the actual distribution has changed, to avoid creating
       ! needless MPI communicators
       if ( (.not.prodbas_tasks_distributed) .or. (previous_n_basbas.ne.n_basbas) .or. &
            (previous_prodbas_nb.ne.prodbas_nb) .or. (previous_n_states.ne.n_states) ) then
          write(info_str,'(A)') & 
            "Parallel distribution of product basis not up to date - distributing tasks."
          call localorb_info (info_str,use_unit,'(2X,A)', OL_high  )

          call prodbas_tasks_distribution ()

          prodbas_tasks_distributed = .true.
          previous_n_basbas = n_basbas
          previous_prodbas_nb = prodbas_nb
          previous_n_states = n_states
       else
          write(info_str,'(A)') & 
            "Parallel distribution of product basis is up to date."
          call localorb_info (info_str,use_unit,'(2X,A)', OL_high  )
       end if

       if(allocated(basbas_wave)) deallocate(basbas_wave)
       if(allocated(fn_l_basbas)) deallocate(fn_l_basbas)
       if(allocated(species_basbas_fns)) deallocate(species_basbas_fns)

 end subroutine shrink_full_auxil_basis_v2
!******
