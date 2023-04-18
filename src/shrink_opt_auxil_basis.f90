!****s* FHI-aims/shrink_opt_auxil_basis
!  NAME
!   shrink_opt_auxil_basis
!  SYNOPSIS

      subroutine shrink_opt_auxil_basis &
        ( )

!  PURPOSE
!  Collects the minimal product basis plus optimized gaussian basis. This size
!  of this basis set should be much smaller than the "full" product basis.
!  Orthogonalizes extra basis functions to atomic ones, for each species.
!  throws out all basis functions which are substantially identical to another one
!  gathers basis functions in a compressed array for numerics
!  implicitly sort all basis functions into blocks of identical l per species
!  labels all basis functions to remember their origin

!  USES
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use species_data
      use spline
      use prodbas
      use basbas_fn_coulomb, only: localize_basbas_fn_coulomb, &
          sort_species_basbas_fn
      use constants
      use mpi_utilities
      use localorb_io, only: use_unit

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


!   input
!   most "input" is found in module species_data

!   max_n_prodbas is the highest n index of the orginal basis used to construct the
!   product basis

!     output

!     all output is found in module basis

!  local variables

!     i_first_fn : First basis function for current atom and current angular momentum shell
!     i_reject_wave : If the innermost maximum of a given wave function lies inside of
!                  r_radial(i_reject_wave), then that wave function will not be used.

      integer i_first_fn, i_first_species_fn
      character l_shell_str
      integer i_reject_wave
      real*8 r_inner_max
      real*8 norm
      real*8 alpha
      real*8 z_eff

      integer n_reject

      integer fn_l_high
      integer fn_l_low
      integer n_max_basbas_fns
      integer num_aux_gaussian


!     Arrays for intermediate storage/handling of radial functions before spline
!     basis_wave is compacted array of basis functions
!     basis_deriv is compacted array of basis function derivatives (if needed)
!     basis_kinetic : directly calculated kinetic energy term for all basis functions
!                which are not free-atom like basis functions. [Treat free-atom like basis functions
!                separately so we can subtract out the free-atom basis-defining potential later.]

      real*8, dimension(:,:), allocatable :: basbas_wave
      integer, dimension(:), allocatable :: fn_l_basbas
      integer, dimension(:), allocatable :: fn_species_basbas(:)

!  counters

      integer i_species, i_species_1
      integer i_basis, i_function, i_atom, i_shell
      integer i_m, i_l, i_grid, i_type
      integer i_n, i_n_1,i_l_1
      integer i_basbas_fn
      integer i_prev

!  functions

      real*8 int_log_mesh
      real*8 int_log_coulomb_metric
      real*8 get_inner_max
      character l_to_str

!  debug

      character*8 fn_type(n_max_basis_fns)
      integer fn_l(n_max_basis_fns,n_species)
      integer fn_n(n_max_basis_fns,n_species)

      integer success
      logical :: use_hse_localization

      character*16 output_name
      integer prodbas_metric

      character(*), parameter :: func = 'shrink_opt_auxil_basis'

!  begin work

      use_hse_localization = (use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse &
          .and. .not. use_dftpt2_and_hse)

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"--------------------------------------------"
        if(include_min_basis(1)) then
           write(use_unit,'(2X,2A)') &
           "Constructing auxiliary basis (minimal product + ", &
           "optimized guassian) ..."
        else
           write(use_unit,'(2X,2A)') &
           "Constructing auxiliary basis set from specified ", &
           "gaussian functions ..."
        endif
        write(use_unit,*)
      endif

! how to orthonormalize the auxiliary basis: = 1, Coulomb metric, =0,
! Normal metric
      prodbas_metric = 0
! Normal metric is used for the SVS version, and Coulomb metric is used
! for V version.
!      if(RI_type == RI_SVS) then
!       prodbas_metric = 0
!      else
!       prodbas_metric = 1
!      endif

!     allocations - local variables for basis function storage

      num_aux_gaussian = sum(n_aux_gaussian(:))
      n_max_basbas_fns=max(n_max_basis_fns*n_max_basis_fns,num_aux_gaussian,10)
      allocate ( basbas_wave(n_max_grid, n_max_basbas_fns) )

      allocate ( fn_l_basbas(n_max_basbas_fns) )
      allocate ( fn_species_basbas(n_max_basbas_fns) )


      i_basbas_fn = 0
      n_basbas = 0
      n_reject = 0
      do i_species = 1, n_species, 1
         i_first_species_fn = i_basbas_fn + 1

         if (n_aux_gaussian(i_species).gt.0) then
              call get_aux_gaussian_basis_fns &
               ( i_species &
               )
         endif

!       set i_reject_wave
         i_reject_wave = innermost_max(i_species)

         i_basis =0

         do i_l = 0, max_l_prodbas(i_species), 1
            i_first_fn = i_basbas_fn + 1

         if (include_min_basis(i_species)) then
!   construct product basis
           do i_function = 1, n_atomic(i_species), 1
             do i_prev = 1, i_function, 1

                fn_l_high = &
                  atomic_l(i_species,i_function) + &
                  atomic_l(i_species,i_function)

                fn_l_low = &
                 abs( atomic_l(i_species,i_function) - &
                  atomic_l(i_species,i_prev))

                if (i_l.ge.fn_l_low .and. &
                          i_l.le.fn_l_high) then

                   i_basbas_fn =  i_basbas_fn +1
                   fn_l_basbas(i_basbas_fn) = i_l

                   do i_grid = 1, n_grid(i_species)
                       basbas_wave(i_grid,i_basbas_fn) &
                        = atomic_wave(i_grid, i_species, i_function) * &
                          atomic_wave(i_grid, i_species, i_prev ) &
                          /r_grid(i_grid, i_species)
                   enddo

!      first normalize the individual wave function
                   if(prodbas_metric .eq. 0) then
                      norm = int_log_mesh ( &
                         basbas_wave(1, i_basbas_fn), &
                         basbas_wave(1, i_basbas_fn), &
                         n_grid (i_species), &
                         r_grid (1, i_species) )
                   elseif (prodbas_metric .eq. 1) then
                       norm = int_log_coulomb_metric ( i_l, &
                        basbas_wave(1, i_basbas_fn), &
                        basbas_wave(1, i_basbas_fn), &
                        n_grid (i_species), &
                        r_grid (1, i_species) )
                   else
                       write(use_unit,'(2X,2A)') &
                         "* Error: No such metric for on-site auxiliary",&
                         " basis orthormalizaton exists."
                   endif

                   do i_grid = 1, n_grid(i_species)
                    basbas_wave (i_grid,i_basbas_fn) = &
                    basbas_wave (i_grid,i_basbas_fn)/sqrt(norm)
                   enddo

!     then orthogonalize it to the previous ones
                  call orthonormalize_prodbas_fn &
                  ( prodbas_metric, i_l, i_basbas_fn, i_first_fn, & 
                    n_max_basbas_fns, basbas_wave, n_grid(i_species), &
                    r_grid(1,i_species), prodbas_acc(i_species), &
                    success &
                  )

                if (success .eq. 1) then

                 if(myid.eq.0) then
                   write(use_unit,'(2X,A,1X,I4,2X,A,1X,I2,2X,A,2I4,2X,A)' ) &
                   "| Species", i_species, "l_shell",i_l , &
                     ":   pair state", &
                     i_function, i_prev, "   accepted."
                  endif

                  i_basis = i_basis + 2*i_l + 1

                elseif (success .eq. -1) then
                 if(myid.eq.0) then
                  write(use_unit,'(2X,A,1X,I4,2X,A,1X,I2,2X,A,2I4,2X,A)' ) &
                  "| Species", i_species, "l_shell", i_l, &
                     ":   pair state", &
                   i_function, i_prev, "  enough product function"
                  write(use_unit,'(2X, A)') &
                    "already for this l channel, stop here ."
                 endif

                  i_basbas_fn = i_basbas_fn - 1

                  n_reject = n_reject +1
                else
                 if(myid.eq.0) then
                  write(use_unit,'(2X,A,1X,I4,2X,A,1X,I2,2X, A,2I4,2X,A)' ) &
                  "| Species", i_species, "l_shell",i_l, &
                      ":   pair state", i_function, &
                   i_prev, "   rejected: Linear dependence."
                 endif

                   i_basbas_fn = i_basbas_fn - 1

                   n_reject = n_reject +1
                endif

                 fn_species_basbas(i_basbas_fn)=i_species

!  end of if i_l
             endif

!      end of do i_prev
             enddo
!      end of do i_function
            enddo
!  end of if include_min_basis
           endif


! include the auxiliary gaussian basis
          if (n_aux_gaussian(i_species).gt.0) then

             do i_function = 1, n_aux_gaussian(i_species), 1

                 i_l_1 = aux_gaussian_l (i_species,i_function)
                 if(i_l_1 .eq. i_l) then

                   i_basbas_fn = i_basbas_fn +1
                   fn_l_basbas(i_basbas_fn) = i_l

                   do i_grid = 1, n_grid(i_species)
                       basbas_wave(i_grid,i_basbas_fn) &
                        = aux_gaussian_wave &
                          ( i_grid, i_species, i_function)
                   enddo

!      verify that chosen integration grid r_radial is suitable for given function
                  r_inner_max = &
                    get_inner_max ( basbas_wave(1,i_basbas_fn), &
                            r_grid(1,i_species), &
                            n_grid(i_species) )

                  if (r_inner_max.ge.r_radial(i_reject_wave,i_species)) &
                    then
!      now enforce orthonormality of current basis function against all
!      previous basis functions of the same atom and angular momentum quantum number

                    call orthonormalize_prodbas_fn &
                    ( prodbas_metric, i_l, i_basbas_fn, i_first_fn, &
                      n_max_basbas_fns,basbas_wave, n_grid(i_species), &
                      r_grid(1,i_species), prodbas_acc(i_species), &
                      success &
                    )

                   else
                      success = -1
                   endif

                   l_shell_str = l_to_str (i_l)
                   alpha = aux_gaussian_alpha(i_species, i_function, 1)

                   if ( success .gt.0 ) then

                     if (myid.eq.0) then
                       write(use_unit,'(2X,A,1X,I4,A,1X,A,1X,A,1X,f10.2,4X,A)') &
                           "| Species", i_species, &
                           ": Gaussian orbital ", l_shell_str, &
                             " alpha = ", alpha, &
                           " accepted."
                     end if

                     i_basis = i_basis + (2*i_l+1)

                   else if ( success.eq.0 ) then
                     if (myid.eq.0) then
                       write(use_unit,'(2X,A,1X,I4,A,1X,A,1X,A,1X,f10.2,4X,A)') &
                           "| Species", i_species, &
                           ": Gaussian orbital ", l_shell_str, &
                             " alpha = ", alpha, &
                           " rejected: Linear dependence"
                     end if

                     i_basbas_fn = i_basbas_fn - 1
                     n_reject = n_reject +1

                   else
                     if (myid.eq.0) then
                       write(use_unit,'(2X,A,1X,I4,A,1X,A,1X,A,1X,f10.2,4X,A)') &
                           "| Species", i_species, &
                           ": Gaussian orbital ", l_shell_str, &
                             " alpha = ", alpha , &
                           " rejected: Insufficient integration grid."
                     end if

                     i_basbas_fn = i_basbas_fn - 1
                     n_reject = n_reject +1
                   endif

                   fn_species_basbas(i_basbas_fn)=i_species

                 endif
              enddo
!  endif guassian function
            endif

            ! Localize basbas_fn within one l-channel
            call localize_basbas_fn_coulomb(i_l, i_first_fn, i_basbas_fn, &
            &      basbas_wave, n_grid(i_species), r_grid(:, i_species), &
            &      Adams_Moulton_integrator, use_hse_localization)

!          write(use_unit,*) "i_l", i_l , l_shell_max(i_species),
!     +          n_max_basis_fns
!       end of do i_l
          enddo

          ! Sort basbas_fn within one species
          if (.not. (use_hse .and. hse_omega_hf /= 0.d0) .or. use_gw_and_hse &
              .or. use_dftpt2_and_hse) then
             call sort_species_basbas_fn(i_basbas_fn-i_first_species_fn+1, &
             &       n_grid(i_species), r_grid(:, i_species), wave_threshold, &
             &       basbas_wave(:, i_first_species_fn:i_basbas_fn), &
             &       fn_l_basbas(i_first_species_fn:i_basbas_fn))
          end if
          ! No need to update fn_species_basbas as this is done only within
          ! one species.

        n_basbas = n_basbas + i_basis * &
                atoms_in_structure (i_species)
!       end of do i_species
      enddo


! SVL
      n_basbas_supercell = n_basbas

      n_basbas_fns = i_basbas_fn

      ! Well, n_basbas is already set.  But resetting doesn't hurt.
      call get_bas_dimensions(n_basbas_fns, fn_l_basbas(1:n_basbas_fns), &
      &                       fn_species_basbas(1:n_basbas_fns), &
      &                       max_basbas_L, max_n_basbas_fnLsp, n_basbas)

      call allocate_basbas()

      basbasfn_l = fn_l_basbas(1:n_basbas_fns)
      basbasfn_species = fn_species_basbas(1:n_basbas_fns)

      call generate_basbasfn(basbas_wave, basbasfn_l, basbasfn_species)

      if(myid.eq.0) then
         write(use_unit, '(2X, A,1X, A, I4, 1X, A, I4, 1X, A)') &
           "| Shrink_opt_auxil_basis :", "there are ", n_basbas_fns, &
           "  radial auxiliary wave functions accepted "
         write(use_unit, '(27X,A,I4,A)') &
            " and", n_reject, " rejected."
         write(use_unit, '(2X,A,I5,A)') &
            "| Shrink_opt_auxil_basis : there are", n_basbas, &
            " partial auxiliary wave functions. "
      endif

      call generate_full_bas(n_basbas_fns, basbasfn_l, basbasfn_species, &
      &                      max_basbas_L, max_n_basbas_fnLsp, n_basbas, &
      &                      basbas_atom, basbas_l, basbas_m, basbas_fn, &
      &                      Lsp2n_basbas_fnLsp, Lsp2basbas_fn, Lsp2basbas_sp, &
      &                      atom2basbas_off, sp2n_basbas_sp)
      max_n_basbas_sp = maxval(sp2n_basbas_sp)

!    product basis distributed over different threads
       call prodbas_tasks_distribution ()


       deallocate(basbas_wave)
       deallocate(fn_l_basbas)
       deallocate(fn_species_basbas)
       end subroutine shrink_opt_auxil_basis
!******
