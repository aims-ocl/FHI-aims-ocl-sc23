!  FIXME: For optimum integrations, the present subroutine should also interpolate all
!         wave functions as splines onto the radial integration grid, and _then_ perform
!         all normalisation integrals on this grid. THEN, we would have a numerically
!         well-defined Hamiltonian matrix later on.
!         The present case will produce exact orthonormality on the log grid, but only
!         approximate orthonormality on the later radial integration grid.
!****s* FHI-aims/shrink_fixed_basis
!  NAME
!   shrink_fixed_basis
!  SYNOPSIS
 subroutine shrink_fixed_basis (        )
!  PURPOSE
!    THIS IS HERE FOR HISTORICAL REASONS ONLY
!    IN PRACTICE SUPERCEDED BY shrink_fixed_basis_phi_thresh
!    Subroutine shrink_fixed_basis:
!    * orthogonalizes extra basis functions to atomic ones, for each species.
!    * throws out all basis functions which are substantially identical to another one
!    * gathers basis functions in a compressed array for numerics
!    * implicitly sort all basis functions into blocks of identical l per species
!    * labels all basis functions to remember their origin
!  USES
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use species_data
      use spline
      use mpi_tasks
      use localorb_io
      use pbc_lists
      use free_atoms
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

      implicit none

!  imported variables

!     input
!     most "input" is found in module species_data

!     output

!     all output is found in module basis

!  local variables

!     n_max : aux variables to determine loop counter maxima
!     i_first_fn : First basis function for current atom and current angular momentum shell
!     i_reject_wave : If the innermost maximum of a given wave function lies inside of
!                  r_radial(i_reject_wave), then that wave function will not be used.

      integer n_max
      integer i_first_fn
      character l_shell_str
      integer i_reject_wave
      real*8 r_inner_max

      integer function_index( n_species, n_basis_types, &
                              n_max_ind_fns )

!     Arrays for intermediate storage/handling of radial functions before spline
!     basis_wave is compacted array of basis functions
!     basis_deriv is compacted array of basis function derivatives (if needed)
!     basis_kinetic : directly calculated kinetic energy term for all basis functions
!                which are not free-atom like basis functions. [Treat free-atom like basis functions
!                separately so we can subtract out the free-atom basis-defining potential later.]

      real*8, dimension(:,:), allocatable :: basis_wave
      real*8, dimension(:,:), allocatable :: basis_deriv
      real*8, dimension(:,:), allocatable :: basis_kinetic
      real*8, dimension(:,:), allocatable :: basis_kinetic_scaled_zora

!     Arrays for a possible reindexing of all spline arrays
      integer :: i_species_index(n_max_basis_fns)
      integer :: perm_basis_fns(n_max_basis_fns)
      integer :: perm_basis_fns_inv(n_max_basis_fns)
      real*8 :: temp_radius(n_max_basis_fns)

      integer :: i_function_2, n_function, i_offset_spl, i_basis_2
      integer :: i_offset_spl_density

      integer :: i_offset
      integer :: species_offset

      integer :: outer_point

!  counters

      integer i_species, i_basis, i_function, i_atom, i_shell
      integer i_m, i_l, i_grid, i_type
      integer i_ionic, i_conf, i_hydro, i_gaussian
      integer i_atomic
      integer :: i_spline

!  functions

      real*8 int_log_mesh
      character l_to_str
      real*8 get_inner_max

!  debug

      integer fn_species(n_max_basis_fns)
      real*8 fn_eigenval(n_max_basis_fns)
      character*8 fn_type(n_max_basis_fns)
      integer fn_l(n_max_basis_fns)
      integer fn_n(n_max_basis_fns)
      character*16 output_name
      real*8 scalar_prod(n_max_basis_fns)
      integer i_prev

      real*8 write_radius, write_function, write_i_r

! Paula: these are for ordering the basis functions, so that the memory use in
! packing would be smaller.
      integer,    allocatable, dimension(:):: perm_basis, perm_basis_inv
      integer,    allocatable, dimension(:):: work
      character*8,allocatable, dimension(:):: work_c
      real*8,     allocatable, dimension(:):: radius

      real*8:: V_radial_deriv, basis_radial_deriv

!  begin work


      call localorb_info( &
           "Assembling full basis from fixed parts.", &
           use_unit,'(2X,A)')

!     allocations - local variables for basis function storage
!     basis_deriv is strictly _not_ needed unless use_basis_gradients is true.
!                 Allocate it here anyway (perhaps unnecessarily) so that
!                 orthonormalize_basis_fn can be called properly
!     FIXME - If we ever need to save memory, this would be a place.

      allocate ( basis_wave (n_max_grid, n_max_basis_fns) )
      allocate ( basis_deriv (n_max_grid, n_max_basis_fns) )
      allocate ( basis_kinetic (n_max_grid, n_max_basis_fns) )
      if (flag_rel .eq. REL_atomic_zora) then
        allocate( basis_kinetic_scaled_zora(n_max_grid,n_max_basis_fns))
      end if

!     determine basis _functions_ first

!     i_function counts distinct radial functions (at, type, n, l)
      i_function = 0
!     n_basis counts total number of accepted basis functions
      i_basis = 0
      n_basis = 0

      fn_eigenval = 0.d0

      do i_species = 1, n_species, 1

!       set i_reject_wave
        i_reject_wave = innermost_max(i_species)

!       FIXME: Patch until atomic wave functions are treated equivalent to
!              ionic, confined, hydrogenic
        i_basis = 0

!       treat each l shell as a block
        do i_l = 0, l_shell_max(i_species), 1

          i_type = 0

!         memorize index of first basis function for given atom and angular momentum shell
          i_first_fn = i_function+1

          if (include_min_basis(i_species)) then
!           treat atomic-like wave functions first
            i_type = i_type + 1
            do i_atomic = 1, n_atomic(i_species), 1

              if (atomic_l(i_species,i_atomic).eq.i_l) then

                i_function = i_function+1
                i_basis = i_basis + (2*i_l+1)

                fn_species(i_function) = i_species
                fn_type(i_function) = "atomic"
                fn_l(i_function) = i_l
                fn_n(i_function) = atomic_n(i_species,i_atomic)
                if (core_n_max(i_l,i_species) .ge. &
                    atomic_n(i_species,i_atomic)) then
                  ! only for our own relativity, store the correct eigenvalue -
                  ! for all but core functions, this has to be zero.
                  fn_eigenval(i_function) = &
                    atomic_eigenval(i_species,i_atomic)
                end if

!               store current basis function in basis_wave
!               also store kinetic energy term of basis function
                do i_grid = 1, n_grid(i_species),1

                  basis_wave(i_grid, i_function) = &
                    atomic_wave(i_grid, i_species, i_atomic)

                  basis_kinetic(i_grid, i_function) = &
                    atomic_kinetic(i_grid, i_species, i_atomic)

                enddo

                if (use_basis_gradients) then
                  do i_grid = 1, n_grid(i_species),1

                    basis_deriv(i_grid, i_function) = &
                      atomic_wave_deriv(i_grid, i_species, i_atomic)

                  enddo
                end if

!               verify that chosen integration grid r_radial is suitable for given function
                r_inner_max = &
                get_inner_max ( basis_wave(1,i_function), &
                                r_grid(1,i_species), &
                                n_grid(i_species) )

!               if a minimal basis function does not meet the criterion, stop the calculation
!               entirely - the grids must be chosen differently.
                if (r_inner_max.lt.r_radial(i_reject_wave,i_species)) &
                  then
                  l_shell_str = l_to_str(i_l)

                  if (myid.eq.0) then
                     write(use_unit,'(1X,A,A)') &
                          "* Warning - the innermost maximum ", &
                          "of the minimal basis function"
                     write(use_unit,'(1X,A,I2,A,A,I2)') &
                          "* ", atomic_n(i_species,i_atomic), &
                          l_shell_str, &
                          " of species ",i_species
                     write(use_unit,'(1X,A,I1,A)') &
                          "* has its innermost maximum inside the ", &
                          i_reject_wave, "th radial integration shell."
                     write(use_unit,'(1X,A,A)') &
                          "* Adjust your integration grid to be ", &
                          "accurate enough before continuing."
                  end if

                  stop
               end if

               function_index(i_species, i_type, i_atomic) = &
                    i_function

            end if

         enddo

      end if

!     treat confined basis functions next
      if (n_conf(i_species).gt.0) then
         i_type = i_type+1

         do i_conf = 1, n_conf(i_species), 1

            if (conf_l(i_species,i_conf).eq.i_l) then

               i_function = i_function+1
               i_shell = conf_n(i_species, i_conf)

               fn_species(i_function) = i_species
               fn_type(i_function) = "confined"
               fn_l(i_function) = i_l
               fn_n(i_function) = i_shell

                if (core_n_max(i_l,i_species) .ge. &
                    i_shell) then
                  ! only for our own relativity, store the correct eigenvalue -
                  ! for all but core functions, this has to be zero.
                  fn_eigenval(i_function) = &
                    confined_eigenval(i_species,i_conf)
                end if


!     copy current basis function to basis_wave
               do i_grid = 1, n_grid(i_species), 1

                  basis_wave(i_grid, i_function) = &
                       confined_wave(i_grid, i_species, i_conf)

                  basis_kinetic(i_grid, i_function) = &
                       confined_kinetic(i_grid, i_species, i_conf)

               enddo

               if (use_basis_gradients) then
                  do i_grid = 1, n_grid(i_species), 1
                     basis_deriv(i_grid, i_function) = &
                          confined_wave_deriv(i_grid, i_species, i_conf)
                  enddo
                end if

!     verify that chosen integration grid r_radial is suitable for given function
                r_inner_max = &
                     get_inner_max ( basis_wave(1,i_function), &
                     r_grid(1,i_species), &
                     n_grid(i_species) )

                if (r_inner_max.ge.r_radial(i_reject_wave,i_species)) &
                     then
!     now enforce orthonormality of current basis function against all
!     previous basis functions of the same atom and angular momentum quantum number

                   call orthonormalize_basis_fn &
                        ( i_function, i_first_fn, basis_wave, &
                        n_grid(i_species), &
                        r_grid(1,i_species), basis_acc(i_species), &
                        basis_kinetic, basis_deriv, &
                        function_index(i_species, i_type, i_conf), .false. &
                        )

                else
                   function_index(i_species,i_type,i_conf) = -1
                end if

                l_shell_str = l_to_str(i_l)

                if ( function_index(i_species,i_type,i_conf).gt.0 &
                     ) then

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Confined orbital ", i_shell, l_shell_str, &
                           "accepted."
                   end if

                  i_basis = i_basis + (2*i_l+1)

                else if &
                  ( function_index(i_species,i_type,i_conf).eq.0 &
                  ) then

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Confined orbital ", i_shell, l_shell_str, &
                           "rejected: Linear dependence."
                   end if

                   i_function = i_function - 1

                else

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Confined orbital ", i_shell, l_shell_str, &
                           "rejected: Insufficient integration grid."
                   end if

                   i_function = i_function - 1

                end if

              end if

            enddo

!         end confined function
          end if

!         treat ionic basis functions next
          if (n_ionic(i_species).gt.0) then
            i_type = i_type+1

            do i_ionic = 1, n_ionic(i_species), 1

              if (ionic_l(i_species, i_ionic).eq.i_l) then

                i_function = i_function+1
                i_shell = ionic_n(i_species, i_ionic)

                fn_species(i_function) = i_species
                fn_type(i_function) = "ionic"
                fn_l(i_function) = i_l
                fn_n(i_function) = i_shell

                if (core_n_max(i_l,i_species) .ge. &
                    ionic_n(i_species,i_ionic)) then
                  ! only for our own relativity, store the correct eigenvalue -
                  ! for all but core functions, this has to be zero.
                  fn_eigenval(i_function) = &
                    ionic_eigenval(i_species,i_ionic)
                end if


!               copy current basis function to basis_wave
                do i_grid = 1, n_grid(i_species), 1

                  basis_wave(i_grid, i_function) = &
                    ionic_wave(i_grid, i_species, i_ionic)

                  basis_kinetic(i_grid, i_function) = &
                    ionic_kinetic(i_grid, i_species, i_ionic)

                enddo

                if (use_basis_gradients) then
                  do i_grid = 1, n_grid(i_species), 1
                    basis_deriv(i_grid, i_function) = &
                      ionic_wave_deriv(i_grid, i_species, i_ionic)
                  enddo
                end if

!               verify that chosen integration grid r_radial is suitable for given function
                r_inner_max = &
                get_inner_max ( basis_wave(1,i_function), &
                            r_grid(1,i_species), &
                            n_grid(i_species) )

                if (r_inner_max.ge.r_radial(i_reject_wave,i_species)) &
                  then
!                 now enforce orthonormality of current basis function against all
!                 previous basis functions of the same atom and angular momentum quantum number

                  call orthonormalize_basis_fn &
                  ( i_function, i_first_fn, basis_wave, &
                    n_grid(i_species), &
                    r_grid(1,i_species), basis_acc(i_species), &
                    basis_kinetic, basis_deriv, &
                    function_index(i_species, i_type, i_ionic), .false. &
                  )

                else
                  function_index(i_species,i_type,i_ionic) = -1
                end if

                l_shell_str = l_to_str(i_l)

                if ( function_index(i_species,i_type,i_ionic).gt.0 &
                   ) then

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Ionic orbital ", i_shell, l_shell_str, &
                           "accepted."
                   end if

                   i_basis = i_basis + (2*i_l+1)

                else if &
                        ( function_index(i_species,i_type,i_ionic).eq.0 &
                        ) then

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Ionic orbital ", i_shell, l_shell_str, &
                           "rejected: Linear dependence."
                   end if

                   i_function = i_function - 1

                else

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                          "| Species ", i_species, &
                          ": Ionic orbital ", i_shell, l_shell_str, &
                          "rejected: Insufficient integration grid."
                  end if

                  i_function = i_function - 1

               end if

            end if

         enddo

!     end ionic function
      end if

!     treat hydrogenic basis functions next
      if (n_hydro(i_species).gt.0) then
         i_type = i_type+1

         do i_hydro = 1, n_hydro(i_species), 1

            if (hydro_l(i_species, i_hydro).eq.i_l) then

               i_function = i_function+1
               i_shell = hydro_n(i_species, i_hydro)

               fn_species(i_function) = i_species
               fn_type(i_function) = "hydro"
               fn_l(i_function) = i_l
               fn_n(i_function) = i_shell

!     copy current basis function to basis_wave
               do i_grid = 1, n_grid(i_species), 1

                  basis_wave(i_grid, i_function) = &
                       hydro_wave(i_grid, i_species, i_hydro)

                  basis_kinetic(i_grid, i_function) = &
                       hydro_kinetic(i_grid, i_species, i_hydro)

               enddo

               if (use_basis_gradients) then
                  do i_grid = 1, n_grid(i_species), 1
                     basis_deriv(i_grid, i_function) = &
                          hydro_wave_deriv(i_grid, i_species, i_hydro)
                  enddo
               end if

!     verify that chosen integration grid r_radial is suitable for given function
               r_inner_max = &
                    get_inner_max ( basis_wave(1,i_function), &
                    r_grid(1,i_species), &
                    n_grid(i_species) )

               if (r_inner_max.ge.r_radial(i_reject_wave,i_species)) &
                    then
!     now enforce orthonormality of current basis function against all
!     previous basis functions of the same atom and angular momentum quantum number

                  call orthonormalize_basis_fn &
                       ( i_function, i_first_fn, basis_wave, &
                       n_grid(i_species), &
                       r_grid(1,i_species), basis_acc(i_species), &
                       basis_kinetic, basis_deriv, &
                       function_index(i_species, i_type, i_hydro), .false. &
                       )

               else
                  function_index(i_species,i_type,i_hydro) = -1
               end if

               l_shell_str = l_to_str(i_l)

               if ( function_index(i_species,i_type,i_hydro).gt.0 &
                    ) then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                          "| Species ", i_species, &
                          ": Hydrogenic orbital ", i_shell, l_shell_str, &
                          "accepted."
                  end if

                  i_basis = i_basis + (2*i_l+1)

               else if &
                       ( function_index(i_species,i_type,i_hydro).eq.0 &
                       ) then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                          "| Species ", i_species, &
                          ": Hydrogenic orbital ", i_shell, l_shell_str, &
                          "rejected: Linear dependence."
                  end if

                  i_function = i_function - 1

                else

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Hydrogenic orbital ", i_shell, &
                           l_shell_str, &
                           "rejected: Insufficient integration grid."
                   end if

                  i_function = i_function - 1

                end if

              end if

            enddo

!         end hydrogenic function
          end if

!         treat Gaussian basis functions next
          if (n_gaussian(i_species).gt.0) then
            i_type = i_type+1

            do i_gaussian = 1, n_gaussian(i_species), 1

              if (gaussian_l(i_species,i_gaussian).eq.i_l) then

                i_function = i_function+1
                i_shell = gaussian_n(i_species, i_gaussian)

                fn_species(i_function) = i_species
                fn_type(i_function) = "gaussian"
                fn_l(i_function) = i_l
                fn_n(i_function) = i_shell

!               copy current basis function to basis_wave
                do i_grid = 1, n_grid(i_species), 1

                  basis_wave(i_grid, i_function) = &
                    gaussian_wave(i_grid, i_species, i_gaussian)

                  basis_kinetic(i_grid, i_function) = &
                    gaussian_kinetic(i_grid, i_species, i_gaussian)

                enddo

                if (use_basis_gradients) then
                  do i_grid = 1, n_grid(i_species), 1
                    basis_deriv(i_grid, i_function) = &
                      gaussian_wave_deriv(i_grid, i_species, i_gaussian)
                  enddo
                end if

!               verify that chosen integration grid r_radial is suitable for given function
                r_inner_max = &
                get_inner_max ( basis_wave(1,i_function), &
                            r_grid(1,i_species), &
                            n_grid(i_species) )

                if (r_inner_max.ge.r_radial(i_reject_wave,i_species)) &
                  then
!                 now enforce orthonormality of current basis function against all
!                 previous basis functions of the same atom and angular momentum quantum number

                  call orthonormalize_basis_fn &
                  ( i_function, i_first_fn, basis_wave, &
                    n_grid(i_species), &
                    r_grid(1,i_species), basis_acc(i_species), &
                    basis_kinetic, basis_deriv, &
                    function_index(i_species, i_type, i_gaussian), .false. &
                  )

                else
                  function_index(i_species,i_type,i_gaussian) = -1
                end if

                l_shell_str = l_to_str(i_l)

                if ( function_index(i_species,i_type,i_gaussian).gt.0 &
                   ) then

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Gaussian orbital ", i_shell, l_shell_str, &
                           "accepted."
                   end if

                  i_basis = i_basis + (2*i_l+1)

                else if &
                  ( function_index(i_species,i_type,i_gaussian).eq.0 &
                  ) then

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Gaussian orbital ", i_shell, l_shell_str, &
                           "rejected: Linear dependence."
                   end if

                  i_function = i_function - 1

                else

                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,1X,I4,A,I3,1X,A,1X,A)') &
                           "| Species ", i_species, &
                           ": Gaussian orbital ", i_shell, l_shell_str, &
                           "rejected: Insufficient integration grid."
                   end if

                  i_function = i_function - 1

                end if

              end if

            enddo

!         end Gaussian function
          end if

        enddo

        n_basis = n_basis + i_basis * atoms_in_structure(i_species)

      enddo

!     now we know the actual basis size, and hence all necessary array dimensions
      n_basis_fns = i_function

!     verify n_states against n_basis
      if (n_states.gt.n_basis) then

         if (myid.eq.0) then
            write(use_unit,*)
            write(use_unit,*) "* You requested ", n_states, &
                 " Kohn-Sham states in the calculation,"
            write(use_unit,*) "* but there are only ", &
                 n_basis, &
                 " available basis states."
            write(use_unit,*) "* Reducing total number of ",&
                 " Kohn-Sham states to ", &
                 n_basis, "."
            write(use_unit,*)
         end if

        n_states = n_basis
      end if

! VB: Deal with specialities for kinetic energy here. If atomic_zora needed, must amend
!     basis_kinetic by various extra terms.
!

      if (flag_rel .eq. REL_atomic_zora)then

         do i_function = 1, n_basis_fns, 1
            do i_grid = 1, n_grid(fn_species(i_function)),1

              ! first, obtain radial derivative of atomic potential
              V_radial_deriv =  val_spline_deriv( dble(i_grid), &
              free_potential_spl(1,1, fn_species(i_function)), &
              n_grid(     fn_species(i_function)          ) ) &
              / (log(r_grid_inc( fn_species(i_function))) &
                    * r_grid(i_grid, fn_species(i_function)))


            ! ... and radial derivative of this radial function
            basis_radial_deriv = basis_deriv(i_grid, i_function)




            ! first, use unmodified kinetic energy part to create renormalization
            ! expression for "scaled ZORA"
             basis_kinetic_scaled_zora(i_grid,i_function) = &
               basis_kinetic(i_grid,i_function) &
               * 2 * light_speed_sq &
               / (  2 * light_speed_sq &
               - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
!
              -  2*light_speed_sq &
!
                * V_radial_deriv &
                * (   basis_radial_deriv &
                -    basis_wave(i_grid,i_function) &
                 / r_grid(i_grid, fn_species(i_function))) &
!
               / (  2*light_speed_sq &
               - free_potential_spl(1,i_grid,fn_species(i_function)))**3



            ! now, modify the kinetic energy part itself for on-site ("atomic") ZORA corrections
            basis_kinetic(i_grid,i_function) = &
               basis_kinetic(i_grid,i_function) &
               * 2 * light_speed_sq &
               / (  2 * light_speed_sq &
               - free_potential_spl(1,i_grid,fn_species(i_function))) &
!
               - light_speed_sq &
               / (  2*light_speed_sq &
               - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
!
                * V_radial_deriv &
                * (   basis_radial_deriv &
                -    basis_wave(i_grid,i_function) &
                 / r_grid(i_grid, fn_species(i_function)))

           end do
         end do

      end if ! REL_atomic_zora

      if (flag_rel .eq. REL_own) then

         do i_function = 1, n_basis_fns, 1

            ! amend kinetic energy by core eigenvalue (all other eigenvales are set to zero above)

            ! Notice that this treatment is still outrightly wrong. Really, we should be
            ! setting the correct kinetic energy term already when the basis functions are
            ! generated. Then, we could simply use the orthonormalizaton above to project out all core
            ! contributions correctly. Now, confined and ionic core radial functions will
            ! get the wrong kinetic energy ... as will any accidentally overlapping core contributions from
            ! the original non-orthonormalized basis functions.

            do i_grid = 1, n_grid(fn_species(i_function)),1

              ! first, obtain radial derivative of atomic potential
              V_radial_deriv =  val_spline_deriv( dble(i_grid), &
              free_potential_spl(1,1, fn_species(i_function)), &
              n_grid(     fn_species(i_function)          ) ) &
              / (log(r_grid_inc( fn_species(i_function))) &
                    * r_grid(i_grid, fn_species(i_function)))


            ! ... and radial derivative of this radial function
            basis_radial_deriv = basis_deriv(i_grid, i_function)

            ! now, modify the kinetic energy part itself for on-site ("atomic") ZORA corrections
            basis_kinetic(i_grid,i_function) = &
               basis_kinetic(i_grid,i_function) &
               * 2 * light_speed_sq &
               / (  2 * light_speed_sq + fn_eigenval(i_function) &
               - free_potential_spl(1,i_grid,fn_species(i_function))) &
!
               - light_speed_sq &
               / (  2*light_speed_sq + fn_eigenval(i_function) &
               - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
!
                * V_radial_deriv &
                * (   basis_radial_deriv &
                -    basis_wave(i_grid,i_function) &
                 / r_grid(i_grid, fn_species(i_function)))

            end do

         end do

      end if ! REL_own

      if (out_basis) then
!       Output accepted basis functions for consistency
        do i_function = 1, n_basis_fns, 1
          l_shell_str = l_to_str(fn_l(i_function))
          if (i_function.lt.10) then
            if (fn_n(i_function).lt.10) then
              write(output_name,'(A2,A1,I1,A1,I1,A1,I1,A1,A1,A4)') &
              fn_type(i_function),"_",fn_species(i_function),"_", &
              i_function, "_", &
              fn_n(i_function),"_",l_shell_str,".dat"
            else if  (fn_n(i_function).lt.100) then
              write(output_name,'(A2,A1,I1,A1,I1,A1,I2,A1,A1,A4)') &
              fn_type(i_function),"_",fn_species(i_function),"_", &
              i_function, "_", &
              fn_n(i_function),"_",l_shell_str,".dat"
            end if
          else if (i_function.lt.100) then
            if (fn_n(i_function).lt.10) then
              write(output_name,'(A2,A1,I1,A1,I2,A1,I1,A1,A1,A4)') &
              fn_type(i_function),"_",fn_species(i_function),"_", &
              i_function, "_", &
              fn_n(i_function),"_",l_shell_str,".dat"
            else if  (fn_n(i_function).lt.100) then
              write(output_name,'(A2,A1,I1,A1,I2,A1,I2,A1,A1,A4)') &
              fn_type(i_function),"_",fn_species(i_function),"_", &
              i_function, "_", &
              fn_n(i_function),"_",l_shell_str,".dat"
            end if
          end if

          if (myid.eq.0) then
             open(50, file=output_name)
             write(50,*) "# ",n_grid(fn_species(i_function))
             do i_grid = 1, n_grid(fn_species(i_function)), 1
                write(50,*) &
                     r_grid(i_grid, fn_species(i_function)), &
                     basis_wave(i_grid, i_function)
             enddo
             close(50)
          end if
!         output kinetic energy term also
          write(output_name,'(A4,A2,A1,I1,A1,I1,A1,A1,A4)') &
            "kin_", &
            fn_type(i_function),"_",fn_species(i_function),"_", &
            fn_n(i_function),"_",l_shell_str,".dat"
          open(50, file=output_name)
          write(50,*) "# ",n_grid(fn_species(i_function))
          do i_grid = 1, n_grid(fn_species(i_function)), 1
            write(50,*) &
              r_grid(i_grid, fn_species(i_function)), &
              basis_kinetic(i_grid, i_function)
          enddo
          close(50)
        enddo
      end if

!!!! VB: BEFORE this point, anything to do with basis function manipulation is handled.
!!!!     AFTER this point, only splining, and manipulation and reorganizing of splines.

!     can now allocate the actual basis storage arrays from module basis .
      call allocate_basis ( )

!     and can store the splined versions of each basis function.
!
!     Note. We also "doctor" the splines, so that they are strictly zero outside
!     a given radius, and so that they provide a smooth non-oscillating transition to zero.
!     This is not identical with outer_radius further below, which is detemined by a threshold
!     parameter, rather than a "strictly zero" criterion.
!
      do i_function = 1, n_basis_fns, 1

        ! determine first "strictly zero" point on logarithmic grid
        outer_point = n_grid(fn_species(i_function))
        do while ( ( outer_point.ge.1 ) &
          .and.(basis_wave(outer_point,i_function).eq.0.d0) )
          outer_point = outer_point - 1
        enddo
        if (outer_point.lt.n_grid(fn_species(i_function))) then
          outer_point = outer_point+1
        end if

        ! create radial function spline
        basis_wave_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline &
        ( basis_wave(1,i_function), outer_point, &
          basis_wave_spl(1,1,i_function) )
        ! fudge outermost segment of spline so that it goes to zero smoothly as
        ! a parabola (but with a discontinuous forst derivative at the first non-zero
        ! logarithmic grid point)
        basis_wave_spl(2,outer_point-1,i_function) = &
              -2.d0 * basis_wave_spl(1,outer_point-1,i_function)
        basis_wave_spl(3,outer_point-1,i_function) = &
              basis_wave_spl(1,outer_point-1,i_function)
        basis_wave_spl(4,outer_point-1,i_function) = 0.d0

        ! create radial kinetic energy part spline
        basis_kinetic_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline &
        ( basis_kinetic(1,i_function), outer_point, &
          basis_kinetic_spl(1,1,i_function) )
        ! fudge outermost segment of spline so that it goes to zero smoothly as
        ! a parabola (but with a discontinuous forst derivative at the first non-zero
        ! logarithmic grid point)
        basis_kinetic_spl(2,outer_point-1,i_function) = &
              -2.d0 * basis_kinetic_spl(1,outer_point-1,i_function)
        basis_kinetic_spl(3,outer_point-1,i_function) = &
              basis_kinetic_spl(1,outer_point-1,i_function)
        basis_kinetic_spl(4,outer_point-1,i_function) = 0.d0

        if (use_basis_gradients) then
          ! create radial derivative spline
          basis_deriv_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
          call cubic_spline &
          ( basis_deriv(1,i_function), outer_point, &
            basis_deriv_spl(1,1,i_function) )
          ! fudge outermost segment of spline so that it goes to zero smoothly as
          ! a parabola (but with a discontinuous forst derivative at the first non-zero
          ! logarithmic grid point)
          basis_deriv_spl(2,outer_point-1,i_function) = &
              -2.d0 * basis_deriv_spl(1,outer_point-1,i_function)
          basis_deriv_spl(3,outer_point-1,i_function) = &
              basis_deriv_spl(1,outer_point-1,i_function)
          basis_deriv_spl(4,outer_point-1,i_function) = 0.d0
        end if

        if (flag_rel .eq. REL_atomic_zora) then
          ! create spline of kinetic energy part needed in scaled zora
          basis_kinetic_scaled_zora_spl &
          (1:n_max_spline,1:n_max_grid,i_function) = 0.d0
          call cubic_spline &
          ( basis_kinetic_scaled_zora(1,i_function), &
            outer_point, &
            basis_kinetic_scaled_zora_spl(1,1,i_function) )
          ! fudge outermost segment of spline so that it goes to zero smoothly as
          ! a parabola (but with a discontinuous forst derivative at the first non-zero
          ! logarithmic grid point)
          basis_kinetic_scaled_zora_spl(2,outer_point-1,i_function) = &
              -2.d0 * basis_kinetic_scaled_zora_spl &
                     (1,outer_point-1,i_function)
          basis_kinetic_scaled_zora_spl(3,outer_point-1,i_function) = &
              basis_kinetic_scaled_zora_spl(1,outer_point-1,i_function)
          basis_kinetic_scaled_zora_spl(4,outer_point-1,i_function)=0.d0
        end if

      enddo ! end loop over basis fn splines

! VB: Test output option for radial & kinetic energy function
!     on a dense logarithmic grid for verification of the
!     behavior near the cutoff.
!
!test
!      open(51,file="write_fn.dat")
!      open(53,file="write_kin.dat")
!
!      write_radius = 0.01 / bohr
!      do while (write_radius.le.(5.5/bohr))
!        write_i_r = invert_log_grid
!     +    (write_radius,r_grid_min(1),r_grid_inc(1))
!
!        write_function = val_spline
!     +    (write_i_r, basis_wave_spl(1,1,4),n_grid(1))
!        write(51,*) write_radius*bohr, write_function
!
!        write_function = val_spline
!     +    (write_i_r, basis_kinetic_spl(1,1,4),n_grid(1))
!        write(53,*) write_radius*bohr, write_function
!
!
!        write_radius = write_radius * 1.002
!      enddo
!
!      close(51)
!      close(53)
!test end

      ! If requested, do a separate verification of the basis function behavior near the cutoff.
      if (force_smooth_cutoff) then
        do i_function = 1, n_basis_fns, 1

          ! determine first "strictly zero" point on logarithmic grid
          outer_point = n_grid(fn_species(i_function))
          do while ( ( outer_point.ge.1 ) &
            .and.(basis_wave(outer_point,i_function).eq.0.d0) )
            outer_point = outer_point - 1
          enddo
          if (outer_point.lt.n_grid(fn_species(i_function))) then
            outer_point = outer_point+1
          end if

          ! radial function first
          if ( ( basis_wave(outer_point-2,i_function) &
                 .gt.smooth_cutoff_limit) .or. &
               ( basis_wave(outer_point-1,i_function) &
                 .gt.smooth_cutoff_limit) ) then

            write(use_unit,'(1X,A)') &
            "* After radial function orthonormalization / packing:"
            write(use_unit,'(1X,A,A,I5,A)') &
            "* Warning: You requested strict checking of the radial ", &
            "function cutoff, and function # ", i_function, &
            " exceeds the threshold at its two outermost points."
            stop
          end if

          ! kinetic function next
          if ( ( basis_kinetic(outer_point-2,i_function) &
                 .gt.smooth_cutoff_limit) .or. &
               ( basis_kinetic(outer_point-1,i_function) &
                 .gt.smooth_cutoff_limit) ) then

            write(use_unit,'(1X,A)') &
            "* After radial function orthonormalization / packing:"
            write(use_unit,'(1X,A,A,I5,A)') &
            "* Warning: You requested strict checking of the radial ", &
            "function cutoff, and kinetic # ", i_function, &
            " exceeds the threshold at its two outermost points."
            stop
          end if

          ! radial derivative next
          if (use_basis_gradients) then
          if ( ( basis_deriv(outer_point-2,i_function) &
                 .gt.smooth_cutoff_limit) .or. &
               ( basis_deriv(outer_point-1,i_function) &
                 .gt.smooth_cutoff_limit) ) then

            write(use_unit,'(1X,A)') &
            "* After radial function orthonormalization / packing:"
            write(use_unit,'(1X,A,A,I5,A)') &
            "* Warning: You requested strict checking of the radial ", &
            "function cutoff, and deriv # ", i_function, &
            " exceeds the threshold at its two outermost points."
            stop
          end if
          end if

        enddo
      end if
      ! end cutoff verification.

!test Check orthonormalisation of all wave functions
!      write(use_unit,*)
!      write (use_unit,*) " Orthonormality of different basis functions:"
!      do i_function = 1, n_basis_fns, 1
!        do i_prev = 1, n_basis_fns, 1
!          if ( (fn_species(i_function).eq.fn_species(i_prev)).and.
!     +         (fn_l(i_function)   .eq.fn_l(i_prev)   )      ) then
!            scalar_prod(i_prev) =
!     +        int_log_mesh (
!     +          basis_wave(1,i_function), basis_wave(1,i_prev),
!     +          n_grid(fn_species(i_function)),
!     +          r_grid(1,fn_species(i_function))
!     +                     )
!            write(use_unit,*) " |",i_function, i_prev, scalar_prod(i_prev)
!          end if
!        enddo
!      enddo
!test end

!  Now index all basis functions across all atoms!
!     i_basis counts basis functions (at, type, n, l, m)
      i_basis = 0

      do i_atom = 1, n_atoms, 1

        do i_l = 0, l_shell_max(species(i_atom)), 1

          i_type = 0

          if (include_min_basis(species(i_atom))) then
!           treat atomic-like wave functions first
            i_type = i_type + 1

            do i_atomic = 1, n_atomic(species(i_atom)), 1

              if (atomic_l(species(i_atom), i_atomic).eq.i_l) then

!               label atomic basis functions
                do i_m = -i_l, i_l, 1

                  i_basis = i_basis+1
                  basis_atom(i_basis) = i_atom
                  basis_l(i_basis) = i_l
                  basis_m(i_basis) = i_m
                  basis_fn(i_basis) = &
                  function_index(species(i_atom), i_type, i_atomic)

                enddo

              end if

            enddo

          end if

!         treat confined basis functions next
          if (n_conf(species(i_atom)).gt.0) then
            i_type = i_type+1

            do i_conf = 1, n_conf(species(i_atom)), 1

              if (conf_l(species(i_atom),i_conf).eq.i_l) then

                if (function_index(species(i_atom), &
                                 i_type, i_conf).gt.0) then

!                 label basis function
                  do i_m = -i_l, i_l, 1
                    i_basis = i_basis+1
                    basis_atom(i_basis) = i_atom
                    basis_l(i_basis) = i_l
                    basis_m(i_basis) = i_m
                    basis_fn(i_basis) = &
                      function_index(species(i_atom),i_type,i_conf)
                  enddo

                end if

              end if

            enddo

          end if

!         treat ionic basis functions next
          if (n_ionic(species(i_atom)).gt.0) then
            i_type = i_type+1

            do i_ionic = 1, n_ionic(species(i_atom)), 1

              if (ionic_l(species(i_atom),i_ionic).eq.i_l) then

                if (function_index(species(i_atom), &
                                 i_type, i_ionic).gt.0) then

!                 label basis function
                  do i_m = -i_l, i_l, 1
                    i_basis = i_basis+1
                    basis_atom(i_basis) = i_atom
                    basis_l(i_basis) = i_l
                    basis_m(i_basis) = i_m
                    basis_fn(i_basis) = &
                      function_index(species(i_atom),i_type,i_ionic)
                  enddo

                end if

              end if

            enddo

          end if

!         treat hydrogenic basis functions next
          if (n_hydro(species(i_atom)).gt.0) then
            i_type = i_type+1

            do i_hydro = 1, n_hydro(species(i_atom)), 1

              if (hydro_l(species(i_atom),i_hydro).eq.i_l) then

                if (function_index(species(i_atom), &
                                 i_type, i_hydro).gt.0) then

!                 label basis function
                  do i_m = -i_l, i_l, 1
                    i_basis = i_basis+1
                    basis_atom(i_basis) = i_atom
                    basis_l(i_basis) = i_l
                    basis_m(i_basis) = i_m
                    basis_fn(i_basis) = &
                      function_index(species(i_atom),i_type,i_hydro)
                  enddo

                end if

              end if

            enddo

          end if

!         treat Gaussian basis functions next
          if (n_gaussian(species(i_atom)).gt.0) then
            i_type = i_type+1

            do i_gaussian = 1, n_gaussian(species(i_atom)), 1

              if (gaussian_l(species(i_atom),i_gaussian).eq.i_l) then

                if (function_index(species(i_atom), &
                                 i_type, i_gaussian).gt.0) then

!                 label basis function
                  do i_m = -i_l, i_l, 1
                    i_basis = i_basis+1
                    basis_atom(i_basis) = i_atom
                    basis_l(i_basis) = i_l
                    basis_m(i_basis) = i_m
                    basis_fn(i_basis) = &
                      function_index(species(i_atom),i_type,i_gaussian)
                  enddo

                end if

              end if

            enddo

          end if

        enddo

      enddo

!  All basis functions are stored and indexed. Now verify consistency.

!     FIXME: This can be removed after testing
      if (n_basis.ne.i_basis) then

         if (myid.eq.0) then
            write(use_unit,*) &
                 "* Inconsistent no. of basis fns in ",&
                 "shrink_fixed_basis."
            write(use_unit,*) &
                 "* n_basis = ", n_basis, ", i_basis = ", i_basis
         end if

         stop
      end if




!     determine the outer radius of each basis function u(r)
!     [i.e. the radius outside of which all integrations may be skipped
!      because u(r) is practically zero]
!     MUST CHECK BOTH ACTUAL FN AND SECOND DERIVATIVE
      do i_function = 1, n_basis_fns, 1
        i_grid = n_grid(fn_species(i_function))
        do while &
        ( (abs(basis_wave(i_grid,i_function)).le.wave_threshold) .and. &
          (abs(basis_kinetic(i_grid,i_function)).le.wave_threshold) &
          .and.(i_grid.gt.1) &
        )
          i_grid = i_grid-1
        enddo
        if (i_grid.le.1) then

           if (myid.eq.0) then
              write(use_unit,'(1X,A,A)') &
                   "* Warning - a basis function is ",&
                   "lower that the requested ", &
                   "threshold value for integrations everywhere."
              write(use_unit,'(1X,A,A)') &
                   "Species : ", species_name(fn_species(i_function))
              write(use_unit,'(1X,A,A)') &
                   "Type    : ", fn_type(i_function)
              write(use_unit,'(1X,A,I3,A,I3)') &
                   "(n,l)   : ", fn_n(i_function), ",", fn_l(i_function)
           end if

           stop
        end if
        if (i_grid.eq.n_grid(fn_species(i_function))) then
          outer_radius(i_function) = &
            r_grid(i_grid,fn_species(i_function))
        else
          outer_radius(i_function) = &
            r_grid(i_grid,fn_species(i_function))
        end if
!test
!          write(use_unit,'(1X,A,A)')
!     +    "Species : ", species_name(fn_species(i_function))
!          write(use_unit,'(1X,A,A)')
!     +    "Type    : ", fn_type(i_function)
!          write(use_unit,'(1X,A,I3,A,I3)')
!     +    "(n,l)   : ", fn_n(i_function), ",", fn_l(i_function)
!        write(use_unit,*) "Function ", i_function, ", outer_radius : ",
!     +    outer_radius(i_function)
!test end
      enddo
      basisfn_l = fn_l(1:n_basis_fns)
      basisfn_species = fn_species(1:n_basis_fns)
      basisfn_type = fn_type(1:n_basis_fns)
      basisfn_n = fn_n(1:n_basis_fns)

      if (use_basis_gradients) then
!       also check first derivative
        do i_function = 1, n_basis_fns, 1
          i_grid = n_grid(fn_species(i_function))
          do while &
          ( (abs(basis_deriv(i_grid,i_function)).le.wave_threshold) &
            .and. (i_grid.gt.1) &
          )
            i_grid = i_grid-1
          enddo
          if (i_grid.eq.n_grid(fn_species(i_function))) then
            if ( r_grid(i_grid,fn_species(i_function)) &
                 .gt.outer_radius(i_function) ) then
              outer_radius(i_function) = &
                r_grid(i_grid,fn_species(i_function))
            end if
          else
            if ( r_grid(i_grid+1,fn_species(i_function)) &
                 .gt.outer_radius(i_function) ) then
              outer_radius(i_function) = &
                r_grid(i_grid,fn_species(i_function))
            end if
          end if
!test
!          write(use_unit,'(1X,A,A)')
!     +    "Species : ", species_name(fn_species(i_function))
!          write(use_unit,'(1X,A,A)')
!     +    "Type    : ", fn_type(i_function)
!          write(use_unit,'(1X,A,I3,A,I3)')
!     +    "(n,l)   : ", fn_n(i_function), ",", fn_l(i_function)
!        write(use_unit,*) "Function ", i_function, ", outer_radius : ",
!     +    outer_radius(i_function)
!test end
        enddo
      end if

!------------ end order of basis functions----------






        atom_radius_sq = 0.d0

        do i_function = 1, n_basis_fns, 1

          outer_radius_sq(i_function) = outer_radius(i_function)**2
          atom_radius_sq(fn_species(i_function)) = &
            max (atom_radius_sq(fn_species(i_function)), &
                 outer_radius_sq(i_function) )

        enddo
        ! Make sure that multipole_radius_free >= atom_radius
        multipole_radius_free_sq = max(multipole_radius_free_sq, atom_radius_sq)
        atom_radius = sqrt(atom_radius_sq)

        basis_fn_atom = .false.

        do i_basis_2 = 1, n_basis, 1
          basis_fn_atom(basis_fn(i_basis_2),basis_atom(i_basis_2)) = &
              .true.
        enddo

        i_offset_spl = 0

        species_offset = 0

        do i_species = 1, n_species, 1

          basis_fn_start_spl(i_species) = i_offset_spl + 1

          i_function_2 = 0
          i_species_index = 0

          do i_basis_2 = 1, n_basis_fns, 1

            if (fn_species(i_basis_2).eq.i_species) then

               i_function_2 = i_function_2 + 1
               i_species_index(i_function_2) = i_basis_2

            end if

          enddo

          n_function = i_function_2
          n_basis_fn_species(i_species) = n_function

!test
!          write(use_unit,*) "n_function = ", n_function
!test end

          do i_basis_2 = 1, n_function, 1
             temp_radius(i_basis_2) = &
                 outer_radius_sq(i_species_index(i_basis_2))
          enddo

!test
!          write(use_unit,*) temp_radius
!test end

          call insertionsort( temp_radius, n_function, perm_basis_fns, &
              perm_basis_fns_inv )

          ! Index array linking the actual order of the spline array to
          ! the original order of radial functions
          ! for Hamiltonian evaluation
          do i_spline = 1, n_function, 1
            i_offset = i_spline + species_offset

            ! store index of spline function as a function of the radial fn index
            perm_basis_fns_spl(i_offset) = &
              species_offset + &
              perm_basis_fns_inv(i_spline)

            ! store the inverse, i.e. the index of the radial function as a
            ! function of the index used in the array of splined radial functions
            i_radial_fn(i_offset) = species_offset + &
              perm_basis_fns(i_spline)

          enddo

          ! store species offset and increment it for the next species
          spline_offset(i_species) = species_offset
          species_offset = species_offset + n_function

          ! and store all basis function splines in arrays in order of
          ! increasing outer_radius per species

          ! in principle, this loop also runs over _spline_ functions, not over
          ! radial functions in their original order ...
          do i_basis_2 = 1, n_function, 1

             i_offset_spl = i_offset_spl + 1

             basis_wave_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
                  basis_wave_spl(1:4,1:n_grid(i_species), &
                  i_species_index(perm_basis_fns_inv(i_basis_2)))

             basis_kinetic_ordered(i_offset_spl,1:4,1:n_grid(i_species)) &
                  = basis_kinetic_spl(1:4,1:n_grid(i_species), &
                  i_species_index(perm_basis_fns_inv(i_basis_2)))

             if(flag_rel==REL_atomic_zora)then
             basis_kinetic_scaled_zora_ordered &
                     (i_offset_spl,1:4,1:n_grid(i_species)) &
                   = basis_kinetic_scaled_zora_spl &
                     (1:4,1:n_grid(i_species), &
                   i_species_index(perm_basis_fns_inv(i_basis_2)))
             end if

             if (use_basis_gradients) then

              basis_deriv_ordered(i_offset_spl,1:4,1:n_grid(i_species)) &
                     = basis_deriv_spl(1:4,1:n_grid(i_species), &
                     i_species_index(perm_basis_fns_inv(i_basis_2)))

             end if

          enddo

       enddo

!test
!          write(use_unit,*) perm_basis_fns_spl
!          write(use_unit,*) i_radial_fn
!test end



      if (myid.eq.0) then
         write(use_unit,*)
         write(use_unit,'(2X,A)') &
              "Basis size parameters after reduction:"
         write (use_unit,'(2X,A,I8)') &
              "| Total number of radial functions: ", &
              n_basis_fns
         write (use_unit,'(2X,A,I8)') &
              "| Total number of basis functions : ", &
              n_basis
         write(use_unit,*)
      end if

!     can now deallocate the temporary basis variables

      if ( allocated(basis_wave) ) then
        deallocate ( basis_wave )
      end if
      if ( allocated(basis_deriv) ) then
        deallocate ( basis_deriv )
      end if
      if ( allocated(basis_kinetic) ) then
        deallocate ( basis_kinetic )
      end if
      if(allocated ( basis_kinetic_scaled_zora))then
         deallocate ( basis_kinetic_scaled_zora)
      end if

!test
!      stop
!test end

      return
 end subroutine shrink_fixed_basis
!******
