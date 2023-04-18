!****s* FHI-aims/shrink_full_auxil_basis
!  NAME
!   shrink_full_auxil_basis
!  SYNOPSIS

      subroutine shrink_full_auxil_basis_aux &
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
      use grids
      use geometry
      use basis
      use species_data
      use spline
      use prodbas
      use basbas_fn_coulomb, only: localize_basbas_fn_coulomb, &
          sort_species_basbas_fn
      use mpi_tasks
      use constants
      use localorb_io
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

      integer i_first_fn, i_first_species_fn
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
      integer i_basis, i_function, i_atom
      integer i_m, i_l, i_grid
      integer i_n, i_n_1,i_l_1
      integer i_basbas_fn
      integer i_prev

!  functions

      real*8 int_log_mesh
      real*8 int_log_coulomb_metric

!  debug

      integer species_atom(n_species)
      integer n_func(n_species)
      integer fn_l(n_max_basis_fns,n_species)
      integer fn_n(n_max_basis_fns,n_species)
      integer fn_index(n_max_basis_fns,n_species)

      logical flag_species_incl
      integer prodbas_metric
      logical :: use_hse_localization

      integer success

      character*150 :: info_str
      character(*), parameter :: func = 'shrink_full_auxil_basis'

!  begin work

      use_hse_localization = ((use_hse .and. hse_omega_hf /= 0.d0 &
          .and. .not. use_gw_and_hse) &
          .or. lrc_pt2_started)

      if(myid.eq.0) then
        write (use_unit,*)
        write (use_unit,*)"--------------------------------------------"
        write (use_unit,'(2X,A)') &
          "Constructing auxiliary basis (full product) ..."
        write (use_unit,*)
      endif

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

        fn_l_basbas = 0
        i_basbas_fn = 0
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
                   i_basis = i_basis + 2*i_l_1 +1

                 endif
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
!!
!! Testing 
!                    write(use_unit,*) " "
!                    write(use_unit,*) "Aux basis:", i_basbas_fn
!                    do i_grid = 1, n_grid(i_species)
!                        write(use_unit,*) i_grid, basbas_wave (i_grid,i_basbas_fn)
!                    enddo       
!!
!!

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

!  end of if i_l
            endif

!      end of if i_n
!            endif
!      end of do i_prev
            enddo
!      end of do i_function
           enddo

           ! Localize basbas_fn within one l-channel.
           call localize_basbas_fn_coulomb(i_l, i_first_fn, i_basbas_fn, &
           &      basbas_wave, n_grid(i_species), r_grid(:, i_species), &
           &      Adams_Moulton_integrator, use_hse_localization)

!          write(use_unit,*) "i_l", i_l , l_shell_max(i_species),
!     +          n_max_basis_fns
!       end of do i_l
          enddo

          ! Sort basbas_fn within one species
          if (.not. (use_hse .and. hse_omega_hf /= 0.d0).or. use_gw_and_hse &
              .or. use_dftpt2_and_hse) then
             ! Do not do this for HSE as we would need the field extent,
             ! which is comparably expensive.
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

! SVL
      n_basbas_supercell = n_basbas

      n_basbas_fns = i_basbas_fn

      ! Well, n_basbas is already set.  But resetting doesn't hurt.
      call get_bas_dimensions(n_basbas_fns, fn_l_basbas(1:n_basbas_fns), &
      &                       species_basbas_fns(1:n_basbas_fns), &
      &                       max_basbas_L, max_n_basbas_fnLsp, n_basbas)

      call allocate_basbas()

      basbasfn_l = fn_l_basbas(1:n_basbas_fns)
      basbasfn_species = species_basbas_fns(1:n_basbas_fns)

      call generate_basbasfn(basbas_wave, basbasfn_l, basbasfn_species)

      if(myid.eq.0) then
       write(use_unit, '(2X, A,1X, A, I8, 1X, A, I8, 1X, A)') &
        "| Shrink_full_auxil_basis :", "there are ", n_basbas_fns, &
          " radial auxiliary wave functions"
        write(use_unit, '(27X,A,I8,A)') &
         " accepted and", n_reject, " rejected."
       write(use_unit, '(2X,A,I20,A)') &
          "| Shrink_full_auxil_basis : there are totally", n_basbas, &
          " partial auxiliary wave functions."
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
       deallocate(species_basbas_fns)
       end subroutine shrink_full_auxil_basis_aux
!******
