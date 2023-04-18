!****h* FHI-aims/free_atoms
!  NAME
!    free_atoms
!  SYNOPSIS
 
     module free_atoms

!  PURPOSE
!    contains all calculated free-atom data, like densities, potentials, wave functions
!
!  Subroutines:
!  * allocate_free_atoms
!  * cleanup_free_atoms
!
!     Free atom densities, potential, wave functions etc:
!
!     free_wave: tabulated radial wave functions of all free species
!     free_wave_spl: splined radial wave functions of all free species
!     free_wave_eigenval : tabulated eigenvalues for all free atomic wave functions.
!     free_potential: Atomic(!) effective potential, tabulated on n_grid radial
!                grid points r_grid - for each species
!     free_potential_spl: Cubic splined version of free_potential
!                Only really needed for core_level_shift
!     free_pot_es : Electrostatic part of the free atom potential
!     free_pot_es_at_zero :  Electrostatic part of the free atom potential at the position of nucleus "zero-limit".
!     free_pot_es_spl : Cubic splined version of free_pot_es, on logarithmic grid
!     free_rho: spherical atomic charge density, tabulated on n_grid
!                radial grid points r_grid - for each species
!     free_rho_spl : cubic splined version of free_rho, on logarithmic grid.
!     free_d_rho_dr_spl : splined version of the radial derivative of the free atom density
!     partition_rho_spl : cubic splined free-atom density but for use in the partition tables only
!     partition_d_rho_dr_spl : splined version of the corresponding radial derivative
!     atom_type : lookup table that maps atom onto initialization type
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
!  USES
      implicit none
!  SOURCE
      real*8, dimension(:,:,:),   allocatable :: free_wave
      real*8, dimension(:,:,:),   allocatable :: free_wave_small
      real*8, dimension(:,:,:,:), allocatable :: free_wave_spl
      real*8, dimension(:,:,:),   allocatable :: free_wave_deriv
      real*8, dimension(:,:,:),   allocatable :: free_wave_small_deriv
      real*8, dimension(:,:,:,:), allocatable :: free_wave_deriv_spl
      real*8, dimension(:,:,:),   allocatable :: free_kinetic
      real*8, dimension(:,:,:,:), allocatable :: free_kinetic_spl
      real*8, dimension(:,:),     allocatable :: free_wave_eigenval
      real*8, dimension(:,:),     allocatable :: free_potential
      real*8, dimension(:,:),     allocatable :: free_pot_es
      real*8, dimension(:),       allocatable :: free_pot_es_at_zero
      real*8, dimension(:,:),     allocatable :: free_rho
      real*8, dimension(:,:,:),   allocatable :: initial_rho
      real*8, dimension(:,:,:),   allocatable :: initial_drho_dr
      real*8, dimension(:,:),     allocatable :: initial_pot_es
      real*8, dimension(:,:,:,:), allocatable :: initial_rho_spl
      real*8, dimension(:,:,:,:), allocatable :: initial_drho_dr_spl
      real*8, dimension(:,:,:),   allocatable :: initial_pot_es_spl
      real*8, dimension(:,:,:),   allocatable :: free_potential_spl
      real*8, dimension(:,:,:),   allocatable :: free_pot_es_spl
      real*8, dimension(:,:,:),   allocatable :: free_rho_spl
      real*8, dimension(:,:,:),   allocatable :: free_drho_dr_spl
      real*8, dimension(:,:,:),   allocatable :: partition_rho_spl
      real*8, dimension(:,:,:),   allocatable :: hartree_partition_rho_spl
      real*8, dimension(:,:,:),   allocatable :: hartree_partition_drho_dr_spl

      ! If we charge-compensate the free-atom superposition for the Hartree potential on the 3D integration grid,
      ! we must also renormalize the appropriate free-atom splines for the Hartree potential
      ! This change is only relevant if the compensate_multipole_errors flag is used.
      real*8, dimension(:,:,:),   allocatable :: renormalized_free_rho_spl
      real*8, dimension(:,:,:),   allocatable :: renormalized_free_drho_dr_spl


      ! For charged periodic systems, a uniform charge background will feel this average potential.
      real*8 :: average_free_es_pot

      ! In case the average potential is computed based on the atoms, here are the
      ! potential and the volume within which the potential is non-zero
      real*8, dimension(:), allocatable :: free_atom_average_es_pot
      real*8, dimension(:), allocatable :: free_atom_average_es_vol

      contains
!******

!****s* free_atoms/allocate_free_atoms
!  NAME
!    allocate_free_atoms
!  SYNOPSIS
        subroutine allocate_free_atoms( )
!  PURPOSE
!    allocation of free atom data
!  USES
      use dimensions,      only : n_ini_type, n_max_grid, n_max_ind_fns, n_max_spline, &
                                  n_species, n_spin, use_basis_gradients, use_density_gradient, &
                                  use_initial_rho, use_partition_deriv, use_prodbas
      use runtime_choices, only : analytic_potential_average, flag_hartree_partition_type, &
                                  multipole_interpolation_style, partition_type, &
                                  flag_rel, REL_x2c, REL_4c_dks
      use rel_x2c_mod,     only : free_den_diff, free_den_diff_spl, free_drho_dr_diff_spl
! INPUT
! none
! OUTPUT
! none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        allocate(free_wave(n_max_grid,n_species,n_max_ind_fns))
        allocate( free_wave_spl &
          (n_max_spline,n_max_grid,n_max_ind_fns,n_species) )
        allocate(free_wave_eigenval(n_species,n_max_ind_fns))
        allocate(free_kinetic(n_max_grid,n_species,n_max_ind_fns))
        allocate(free_kinetic_spl &
          (n_max_spline,n_max_grid,n_max_ind_fns,n_species) )

        allocate (free_potential(n_max_grid, n_species))
        allocate (free_potential_spl(n_max_spline,n_max_grid,n_species))
        allocate (free_pot_es(n_max_grid, n_species))
        allocate ( free_pot_es_at_zero(n_species))
        allocate (free_pot_es_spl(n_max_spline,n_max_grid, n_species))
        allocate (free_rho(n_max_grid, n_species))
        allocate (free_rho_spl(n_max_spline,n_max_grid, n_species))
        allocate (renormalized_free_rho_spl(n_max_spline,n_max_grid, n_species))

        if (flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) then
           allocate(free_wave_small(n_max_grid,n_species,n_max_ind_fns))
           allocate (free_den_diff(n_max_grid, n_species))
           allocate (free_den_diff_spl(n_max_spline,n_max_grid, n_species))
        end if

        if (use_initial_rho) then
           allocate (initial_rho(n_max_grid, n_ini_type, n_spin))
           allocate (initial_rho_spl(n_max_spline, n_max_grid, &
                n_ini_type, n_spin))
           allocate (initial_pot_es(n_max_grid, n_ini_type))
           allocate (initial_pot_es_spl(n_max_spline, n_max_grid, &
                n_ini_type))
        end if

        if (use_basis_gradients) then
           allocate(free_wave_deriv(n_max_grid,n_species,n_max_ind_fns))
           allocate( free_wave_deriv_spl &
                (n_max_spline,n_max_grid,n_max_ind_fns,n_species) )
        end if       

        if (flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) then
           if(.not.allocated(free_wave_deriv))then
             allocate(free_wave_deriv(n_max_grid,n_species,n_max_ind_fns))
           endif
           if(.not.allocated(free_wave_small_deriv))then
             allocate(free_wave_small_deriv(n_max_grid,n_species,n_max_ind_fns))
           endif
        end if
 
        if (use_density_gradient) then
          allocate (free_drho_dr_spl(n_max_spline,n_max_grid,n_species))
          allocate (renormalized_free_drho_dr_spl(n_max_spline,n_max_grid,n_species))
          if (flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) then
             allocate (free_drho_dr_diff_spl(n_max_spline,n_max_grid,n_species))
          endif
          if (use_initial_rho) then
             allocate (initial_drho_dr(n_max_grid, n_ini_type, n_spin))
             allocate (initial_drho_dr_spl(n_max_spline, n_max_grid, &
                  n_ini_type, n_spin))
          end if
       end if

        !---------shanghui add for free_drho_dr_spl for gradient_partition----
        if(use_partition_deriv) then
        allocate (free_drho_dr_spl(n_max_spline,n_max_grid,n_species))
        endif
        !---------shanghui end add for free_drho_dr_spl for gradient_partition----

        ! separate allocation of arrays for densities used in partition tabs ...
        ! not much memory, but if we ever make stratmann the default, these
        ! allocations should not be performed any more, on the BG they contribute
        ! to the per-thread memory overhead

        allocate (hartree_partition_rho_spl(n_max_spline,n_max_grid, n_species))
        if (multipole_interpolation_style.eq.1) then
            ! currently only a corner case, may be needed if we ever implement the
            ! "grid derivative" force terms ... but then for the partition_tab,
            ! not for the hartree_partition_tab ...
            allocate (hartree_partition_drho_dr_spl(n_max_spline,n_max_grid,n_species))
        else
            ! dummy allocation, array is not used
            allocate (hartree_partition_drho_dr_spl(1,1,1))
        end if

        if ( (flag_hartree_partition_type.ne.partition_type) .or. use_prodbas ) then
          ! for most cases, this extra allocation is superfluous anyway and should be 
          ! tightened drastically if ever another memory problem crops up,
          ! for example on the BG
          allocate (partition_rho_spl(n_max_spline,n_max_grid, n_species))
        end if

        if (analytic_potential_average) then
           allocate(free_atom_average_es_pot(n_species))
           allocate(free_atom_average_es_vol(n_species))
        end if

        end subroutine allocate_free_atoms
!******

!****s* free_atoms/cleanup_free_atoms
!  NAME
!    cleanup_free_atoms
!  SYNOPSIS
        subroutine cleanup_free_atoms( )
!  PURPOSE
!    deallocate free atom data
!  USES
      use rel_x2c_mod,     only : free_den_diff, free_den_diff_spl, free_drho_dr_diff_spl
      implicit none
! INPUT
! none
! OUTPUT
! none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        if (allocated(free_wave)) then
          deallocate(free_wave)
        end if
        if (allocated(free_wave_small)) then
          deallocate(free_wave_small)
        end if
        if (allocated(free_wave_spl)) then
           deallocate(free_wave_spl)
        end if
        if (allocated(free_wave_eigenval)) then
          deallocate(free_wave_eigenval)
        end if
        if (allocated(free_wave_deriv)) then
           deallocate(free_wave_deriv)
        end if
        if (allocated(free_wave_small_deriv)) then
           deallocate(free_wave_small_deriv)
        end if
        if (allocated(free_wave_deriv_spl)) then
           deallocate(free_wave_deriv_spl)
        end if
        if (allocated(free_kinetic)) then
          deallocate(free_kinetic)
        end if
        if (allocated(free_kinetic_spl)) then
          deallocate(free_kinetic_spl)
        end if
          if (allocated(free_potential)) then
            deallocate( free_potential )
          end if
          if (allocated(free_potential_spl)) then
            deallocate( free_potential_spl )
          end if
          if (allocated(free_pot_es)) then
            deallocate( free_pot_es )
          end if
          if(allocated( free_pot_es_at_zero))then
             deallocate(free_pot_es_at_zero)
          end if
          if (allocated(free_pot_es_spl)) then
            deallocate( free_pot_es_spl )
          end if
          if (allocated(free_rho)) then
            deallocate( free_rho )
          end if
          if (allocated(free_rho_spl)) then
            deallocate( free_rho_spl )
          end if
          if (allocated(free_den_diff)) then
            deallocate( free_den_diff )
          end if
          if (allocated(free_den_diff_spl)) then
            deallocate( free_den_diff_spl )
          end if
          if (allocated(renormalized_free_rho_spl)) then
            deallocate( renormalized_free_rho_spl )
          end if
          if (allocated(partition_rho_spl)) then
            deallocate( partition_rho_spl )
          end if
          if (allocated(hartree_partition_rho_spl)) then
            deallocate( hartree_partition_rho_spl )
          end if
          if (allocated(free_drho_dr_spl)) then
            deallocate (free_drho_dr_spl)
          end if
          if (allocated(free_drho_dr_diff_spl)) then
            deallocate (free_drho_dr_diff_spl)
          end if
          if (allocated(renormalized_free_drho_dr_spl)) then
            deallocate (renormalized_free_drho_dr_spl)
          end if
          if (allocated(hartree_partition_drho_dr_spl)) then
            deallocate (hartree_partition_drho_dr_spl)
          end if
          if (allocated(initial_rho_spl)) then
             deallocate( initial_rho_spl )
          end if
          if (allocated(initial_rho)) then
            deallocate( initial_rho )
          end if
          if (allocated(initial_drho_dr_spl)) then
             deallocate( initial_drho_dr_spl )
          end if
          if (allocated(initial_drho_dr)) then
             deallocate( initial_drho_dr )
          end if
          if (allocated(initial_pot_es)) then
             deallocate( initial_pot_es )
          end if
          if (allocated(initial_pot_es_spl)) then
             deallocate( initial_pot_es_spl )
          end if
          if (allocated(free_atom_average_es_pot)) then
             deallocate(free_atom_average_es_pot)
          end if
          if (allocated(free_atom_average_es_vol)) then
             deallocate(free_atom_average_es_vol)
          end if

        end subroutine cleanup_free_atoms
!******

!****s* free_atoms/initialize_initial_rho
!  NAME
!    initialize_initial_rho
!  SYNOPSIS
      subroutine initialize_initial_rho()
!  PURPOSE
!    get initial density for spin-polarized atoms 
!  USES
      use dimensions,  only : use_density_gradient
      use localorb_io, only : use_unit, localorb_info
      implicit none
! INPUT
! none
! OUTPUT
! none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

       ! local variables
       character*100 :: info_str

       initial_rho     = 0.d0
       initial_rho_spl = 0.d0
       if (use_density_gradient) then
          initial_drho_dr     = 0.d0
          initial_drho_dr_spl = 0.d0
       end if

       write (info_str, '(2X,A)') "Spin-polarized or charged system:"
       call localorb_info(info_str, use_unit,'(A)')
       write (info_str, '(2X,A,A)') "Charge density " , &
            "initialized according to selected moments and charges."
       call localorb_info(info_str, use_unit,'(A)')

       call get_free_atoms_polarized()

      end subroutine initialize_initial_rho
!******

!****s* free_atoms/get_atomic_occ_numbers
!  NAME
!    get_atomic_occ_numbers
!  SYNOPSIS
      subroutine get_atomic_occ_numbers (n_electrons, moment, &
           kind_of_initial, current_species, atomic_occ_numbers, write_out)
!  PURPOSE
!    calculate occupation numbers for free atoms
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  USES
      use dimensions,  only : n_wave_max, l_wave_max, n_spin
      use localorb_io, only : use_unit, localorb_info
      implicit none
!  ARGUMENTS
      real*8 :: n_electrons
      real*8, dimension(n_spin) :: n_electrons_polarized
      real*8 :: moment
      integer :: kind_of_initial
      integer :: current_species
      logical, optional :: write_out
      real*8, dimension(1:n_wave_max, 0:l_wave_max, n_spin) :: &
           atomic_occ_numbers
!  INPUTS
!    o n_electrons -- number of electrons
!    o number of electrons in each spin channel
!    o moment -- specified spin moment to result from the occupation numbers
!    o kind_of_initial -- initialization_type
!    o current_species -- species number of the current atom to be initialized. 
!                         This is needed here only for verification purposes
!                         (we check whether the needed atomic shells were
!                         correctly set up when the input was read).
!    o write_out -- whether or not to print output
!  OUTPUTS
!    o atomic_occ_numbers -- final initialized atomic occupation numbers 
!  SOURCE
      ! local variables
      integer :: spin_multiplier
      character*100 :: info_str
      logical :: local_write_out
      integer :: total_electrons

      ! counter
      integer :: i_spin
      integer :: i_l
      integer :: i_n
      integer :: start_n
      integer :: start_l

      if (present(write_out)) then
         local_write_out = write_out
      else
         local_write_out = .false.
      end if

      ! Store n_electrons here
      total_electrons = n_electrons

      atomic_occ_numbers = 0.d0

      ! choose orbitals in the order of the "step"-scheme (?)
      ! which gives roughly the right order of orbitals
      ! if you neglect the exceptions
      ! each orbital is then occupied according to hund's rule

      if (local_write_out) then
         call localorb_info('', use_unit)
         write(info_str,'(2X,A,F8.4)') "# electrons", n_electrons
         call localorb_info(info_str, use_unit,'(A)')
      end if

      spin_multiplier = mod(n_spin, 2) + 1

      if (kind_of_initial .eq. 1) then
         if (local_write_out) then
            call localorb_info( &
                 "apply hund's rule on atomic states.", &
                 use_unit,'(2X,A)')
         else if (kind_of_initial .eq. 2) then
            write(info_str,'(2X,A,F8.4)') &
                 "desired atomic moment for this type:", moment
            call localorb_info(info_str)
         end if
      end if

      if (local_write_out) then
         write(info_str,'(2X,A)') &
              "List of occupation numbers for free atom basis: "
         call localorb_info(info_str,use_unit,'(A)')
      end if
      if (kind_of_initial .eq. 1) then
         if (local_write_out) then
            if (n_spin .eq. 2) then
               write(info_str,'(5X,A,4X,A,3X,A,1X,A)') &
                    "n", "l", "spin up", "spin down"
            else
               write(info_str,'(5X,A,4X,A,4X,A)') &
                    "n", "l", "occ"
            end if
            call localorb_info(info_str)
         end if

         ! just apply hund's rule so that atomic moment is fixed

         ! The check_occupation_dimensions is needed because the infrastructure for
         ! this routine relies on the initial specification of valence shells
         ! for the neutral atom. There are, however, rare corner cases (unphysically
         ! charged negative ions) that could exceed the array allocation and initialization
         ! for the neutral atom. This can be circumvented by adjusting control.in, which
         ! is now documented below (see check_occupation_dimensions subroutine header). 
         ! It's a hack; if anyone really needs this in production
         ! (not just by accident after specifying an unphysical negative ion),
         ! please consider contributing the right infrastructure. It's not hard; 
         ! just a day of coding and testing, I believe.
         ! 
         ! The present version just prevents a run with incorrect dimensions, which is
         ! rather important from a user perspective.

         start_n = 1
         start_l = 0
         do while (n_electrons .gt. 0)

            i_n = start_n
            do i_l = start_l, 0, -1

               if (n_electrons .ge. 2.d0 * (2*i_l + 1)) then
               ! enough electrons to fill the complete subshell
                  do i_spin = 1, n_spin, 1
                     ! check dimension for corner cases first, fill then..
                     call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          spin_multiplier * (2*i_l + 1)
                  end do
               else
               ! partially filled subshell, so apply hund's rule
                  if (n_spin .eq. 1) then
                  ! no spin-polarization, so just put remaining electrons
                  ! in the last shell
                     i_spin = 1
                     ! check dimension for corner cases first, fill then..
                     call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          n_electrons
                  else
                  ! hund's rule
                     if (n_electrons .le. (2*i_l + 1)) then
                        i_spin = 1
                        ! check dimension for corner cases first, fill then..
                        call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons
                     else
                        i_spin = 1
                        ! check dimension for corner cases first, fill then..
                        call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             (2*i_l + 1)
                        i_spin = 2
                        ! check dimension for corner cases first, fill then..
                        call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons - (2*i_l + 1)
                     end if
                  end if
               end if

               if (local_write_out) then
                  if (n_spin .eq. 2) then
                     write (info_str,'(2X,I4,1X,I4,1X,F8.4,F8.4)') &
                          i_n, i_l, &
                          atomic_occ_numbers(i_n, i_l, 1), &
                          atomic_occ_numbers(i_n, i_l, 2)
                  else
                     write (info_str,'(2X,I4,1X,I4,1X,F8.4)') i_n, i_l, &
                          atomic_occ_numbers(i_n, i_l, 1)
                  end if
                  call localorb_info( info_str )
               end if

               n_electrons = n_electrons - 2.d0 * (2*i_l + 1)
               if (n_electrons .le. 0) then
                  exit
               end if

               i_n = i_n+1
            end do

            start_n = start_n + 1

            if (n_electrons .gt. 0) then

               i_n = start_n
               do i_l = start_l, 0, -1

                  if (n_electrons .ge. 2.d0 * (2*i_l + 1)) then
                  ! enough electrons to fill the complete subshell
                     do i_spin = 1, n_spin, 1
                        ! check dimension for corner cases first, fill then..
                        call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             spin_multiplier * (2*i_l + 1)
                     end do
                  else
                  ! partially filled subshell, so apply hund's rule
                     if (n_spin .eq. 1) then
                  ! no spin-polarization, so just put remaining electrons
                  ! in the last shell
                        i_spin = 1
                        ! check dimension for corner cases first, fill then..
                        call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons
                     else
                  ! hund's rule
                        if (n_electrons .le. (2*i_l + 1)) then
                           i_spin = 1
                           ! check dimension for corner cases first, fill then..
                           call check_occupation_dimensions & 
                                ( i_n, i_l, i_spin, current_species)
                           atomic_occ_numbers(i_n, i_l, i_spin) = &
                                n_electrons
                        else
                           i_spin = 1
                           call check_occupation_dimensions & 
                                ( i_n, i_l, i_spin, current_species)
                           atomic_occ_numbers(i_n, i_l, i_spin) = &
                                (2*i_l + 1)
                           i_spin = 2
                           call check_occupation_dimensions & 
                                ( i_n, i_l, i_spin, current_species)
                           atomic_occ_numbers(i_n, i_l, i_spin) = &
                                n_electrons - (2*i_l + 1)
                        end if
                     end if
                  end if

                  if (local_write_out) then
                     if (n_spin .eq. 2) then
                        write (info_str,'(2X,I4,1X,I4,1X,F8.4,F8.4)') &
                             i_n, i_l, &
                             atomic_occ_numbers(i_n, i_l, 1), &
                             atomic_occ_numbers(i_n, i_l, 2)
                     else
                        write (info_str,'(2X,I4,1X,I4,1X,F8.4)') &
                             i_n, i_l, atomic_occ_numbers(i_n, i_l, 1)
                     end if
                     call localorb_info( info_str )
                  end if

                  n_electrons = n_electrons - 2.d0 * (2*i_l + 1)
                  if (n_electrons .le. 0) then
                     exit
                  end if

                  i_n = i_n+1
               end do
               start_l = start_l + 1
            end if
         end do

         ! For some atoms we introduce special cases. Hund's rules are all nice and well, but
         ! they are ambiguous for some atoms. For example, Hund's rules with the order of shells above
         ! for Cu say 4s^2 3d^9, but the actual occupation (in theory and experiment) is
         ! 4s^1 3d^10
         if (total_electrons.eq.29) then ! This is a Cu-like atom: Ni^-, Cu(neutral), Zn^-, Ga^2- etc.
            write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
            call localorb_info( info_str )
            write(info_str,'(2X,A)') "* Found a Cu-like atom (29 electrons) - enforcing 4s^1 3d^10 ."
            call localorb_info( info_str )
            if (n_spin.eq.1) then
              atomic_occ_numbers(4, 0, 1) =  1.0
              atomic_occ_numbers(3, 2, 1) = 10.0
            else
              atomic_occ_numbers(4, 0, 2) =  1.d-8
              atomic_occ_numbers(3, 2, 2) =  5.0 - 1.d-8
            end if
         else if (total_electrons.eq.47) then ! This is a Ag-like atom
            write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
            call localorb_info( info_str )
            write(info_str,'(2X,A)') "* Found a Ag-like atom (47 electrons) - enforcing 5s^1 4d^10 ."
            call localorb_info( info_str )
            if (n_spin.eq.1) then
              atomic_occ_numbers(5, 0, 1) =  1.0
              atomic_occ_numbers(4, 2, 1) = 10.0
            else
              atomic_occ_numbers(5, 0, 2) =  1.d-8
              atomic_occ_numbers(4, 2, 2) =  5.0 - 1.d-8
            end if
         else if (total_electrons.eq.58) then ! this is a Ce-like atom
            write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
            call localorb_info( info_str )
            write(info_str,'(2X,A)') "* Found a Ce-like atom (58 electrons) - enforcing 6s^2 5d^1 4f^1 (triplet)."
            call localorb_info( info_str )
            if (n_spin.eq.1) then
              atomic_occ_numbers(6, 0, 1) =  2.0
              atomic_occ_numbers(5, 2, 1) =  1.0
              atomic_occ_numbers(4, 3, 1) =  1.0
            else
              atomic_occ_numbers(6, 0, 1) =  1.
              atomic_occ_numbers(5, 2, 1) =  1.
              atomic_occ_numbers(4, 3, 1) =  1.
            end if
         else if (total_electrons.eq.64) then ! this is a Gd-like atom
            write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
            call localorb_info( info_str )
            write(info_str,'(2X,A)') "* Found a Gd-like atom (64 electrons) - enforcing 6s^2 5d^1 4f^7 (nonet)."
            call localorb_info( info_str )
            if (n_spin.eq.1) then
              atomic_occ_numbers(6, 0, 1) =  2.0
              atomic_occ_numbers(5, 2, 1) =  1.0
              atomic_occ_numbers(4, 3, 1) =  7.0
            else
              atomic_occ_numbers(6, 0, 1) =  1.
              atomic_occ_numbers(5, 2, 1) =  1.
              atomic_occ_numbers(4, 3, 1) =  7.
              atomic_occ_numbers(6, 0, 2) =  1.
              atomic_occ_numbers(5, 2, 2) =  0.
              atomic_occ_numbers(4, 3, 2) =  0.
            end if
         else if (total_electrons.eq.79) then ! This is a Au-like atom
            write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
            call localorb_info( info_str )
            write(info_str,'(2X,A)') "* Found a Au-like atom (79 electrons) - enforcing 6s^1 5d^10 ."
            call localorb_info( info_str )
            if (n_spin.eq.1) then
              atomic_occ_numbers(6, 0, 1) =  1.0
              atomic_occ_numbers(5, 2, 1) = 10.0
            else
              atomic_occ_numbers(6, 0, 2) =  1.d-8
              atomic_occ_numbers(5, 2, 2) =  5.0 - 1.d-8
            end if
         end if

      else if (kind_of_initial .eq. 2) then

         ! atomic moment is chosen

         ! evaluate n_electrons_polarized according
         ! to selected initial moment

         n_electrons_polarized(1) = (n_electrons + moment) / 2.d0
         n_electrons_polarized(2) = (n_electrons - moment) / 2.d0

         if (local_write_out) then
            write(info_str,'(2X,A,F8.4,F8.4)') "# electrons polarized", &
                 n_electrons_polarized(1), n_electrons_polarized(2)
            call localorb_info( info_str )
            write(info_str,'(5X,A,4X,A,3X,A,1X,A)') &
                 "n", "l", "spin up", "spin down"
            call localorb_info( info_str )
         end if
         do i_spin = 1, n_spin, 1

            start_n = 1
            start_l = 0

            do while (n_electrons_polarized(i_spin) .gt. 0)

               i_n = start_n
               do i_l = start_l, 0, -1

                  if (n_electrons_polarized(i_spin) .ge. (2*i_l + 1)) &
                       then
                     ! enough electrons to fill the complete subshell
                     call check_occupation_dimensions & 
                     ( i_n, i_l, i_spin, current_species)
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          (2*i_l + 1)
                  else
                     call check_occupation_dimensions & 
                          ( i_n, i_l, i_spin, current_species)
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          n_electrons_polarized(i_spin)
                  end if

!                  write(use_unit,*) i_spin, n_electrons_polarized(i_spin),
!     +                 moment
                  if (((i_spin .eq. 2) .or. &
                       ((i_spin .eq. 1) .and. &
                       (n_electrons_polarized(i_spin) .le. moment))) &
                       .and. (local_write_out)) then
                     write (info_str,'(2X,I4,1X,I4,1X,F8.4,1X,F8.4)') &
                          i_n, i_l, &
                          atomic_occ_numbers(i_n, i_l, 1), &
                          atomic_occ_numbers(i_n, i_l, 2)
                     call localorb_info( info_str )
                  end if

                  n_electrons_polarized(i_spin) = &
                       n_electrons_polarized(i_spin) &
                       - (2*i_l + 1)
                  if (n_electrons_polarized(i_spin) .le. 0) then
                     exit
                  end if

                  i_n = i_n+1
               end do

               start_n = start_n + 1

               if (n_electrons_polarized(i_spin) .gt. 0) then

                  i_n = start_n
                  do i_l = start_l, 0, -1

                     if (n_electrons_polarized(i_spin) .ge. (2*i_l + 1)) &
                          then
                     ! enough electrons to fill the complete subshell
                        call check_occupation_dimensions & 
                             ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             (2*i_l + 1)
                     else
                        call check_occupation_dimensions & 
                             ( i_n, i_l, i_spin, current_species)
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons_polarized(i_spin)
                     end if

                     if (((i_spin .eq. 2) .or. &
                          ((i_spin .eq. 1) .and. &
                          (n_electrons_polarized(i_spin) .le. moment))) &
                          .and. (local_write_out)) then
                        write (info_str,'(2X,I4,1X,I4,1X,F8.4,1X,F8.4)') &
                             i_n, i_l, &
                             atomic_occ_numbers(i_n, i_l, 1), &
                             atomic_occ_numbers(i_n, i_l, 2)
                        call localorb_info( info_str )
                     end if

                     n_electrons_polarized(i_spin) = &
                          n_electrons_polarized(i_spin) - (2*i_l + 1)

                     if (n_electrons_polarized(i_spin) .le. 0) then
                        exit
                     end if

                     i_n = i_n+1
                  end do
                  start_l = start_l + 1
               end if
            end do
         end do

         ! For some atoms we introduce special cases. Hund's rules are all nice and well, but
         ! they are ambiguous for some atoms. For example, Hund's rules with the order of shells above
         ! for Cu say 4s^2 3d^9, but the actual occupation (in theory and experiment) is
         ! 4s^1 3d^10
         if (total_electrons.eq.29) then ! This is a Cu-like atom: Ni^-, Cu(neutral), Zn^-, Ga^2- etc.
            if (moment.eq.1.0) then      
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Cu-like atom (29 electrons) - enforcing 4s^1 3d^10 ."
              call localorb_info( info_str )
              if (n_spin.eq.2) then
                atomic_occ_numbers(4, 0, 2) =  1.d-8
                atomic_occ_numbers(3, 2, 2) =  5.0 - 1.d-8
              end if
            end if
         else if (total_electrons.eq.47) then ! This is a Ag-like atom.
            if (moment.eq.1.0) then      
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Ag-like atom (47 electrons) - enforcing 5s^1 4d^10 ."
              call localorb_info( info_str )
              if (n_spin.eq.2) then
                atomic_occ_numbers(5, 0, 2) =  1.d-8
                atomic_occ_numbers(4, 2, 2) =  5.0 - 1.d-8
              end if
            end if
         else if (total_electrons.eq.58) then
            if (moment.eq.2.) then
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Ce-like atom (58 electrons) - enforcing 6s^2 5d^1 4f^1 (triplet)."
              call localorb_info( info_str )
              if (n_spin.eq.2) then
                atomic_occ_numbers(6, 0, 1) =  1.
                atomic_occ_numbers(5, 2, 1) =  1.
                atomic_occ_numbers(4, 3, 1) =  1.
              end if
            else if (moment.eq.0.) then
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Ce-like atom (58 electrons) - enforcing 6s^2 5d^1 4f^1 (unpolarized!)."
              call localorb_info( info_str )
              if (n_spin.eq.2) then
                atomic_occ_numbers(5, 2, 1) =  0.5
                atomic_occ_numbers(4, 3, 1) =  0.5
                atomic_occ_numbers(5, 2, 1) =  0.5
                atomic_occ_numbers(4, 3, 1) =  0.5
              end if
            end if
         else if (total_electrons.eq.64) then
            if (moment.eq.8.) then
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Gd-like atom (64 electrons) - enforcing 6s^2 5d^1 4f^7 (nonet)."
              call localorb_info( info_str )
              if (n_spin.eq.2) then
                atomic_occ_numbers(6, 0, 1) =  1.
                atomic_occ_numbers(5, 2, 1) =  1.
                atomic_occ_numbers(4, 3, 1) =  7.
                atomic_occ_numbers(6, 0, 2) =  1.
                atomic_occ_numbers(5, 2, 2) =  0.
                atomic_occ_numbers(4, 3, 2) =  0.
              end if
            else
              write(info_str,'(1X,A)') "* WARNING: Found a Gd-like atom (64 electrons) with explicit moment for initialization."
              call localorb_info( info_str )
              write(info_str,'(1X,A)') "* WARNING: I do not know the correct atomic occupation of shells - it may be wrong!"
              call localorb_info( info_str )
            end if
         else if (total_electrons.eq.79) then ! This is a Au-like atom.
            if (moment.eq.1.0) then      
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Au-like atom (79 electrons) - enforcing 6s^1 5d^10 ."
              call localorb_info( info_str )
              if (n_spin.eq.2) then
                atomic_occ_numbers(6, 0, 2) =  1.d-8
                atomic_occ_numbers(5, 2, 2) =  5.0 - 1.d-8
              end if
            end if
         end if

      else
         write(use_unit,*) "chosen kind of initialization not implemented"
         write(use_unit,*) "* Abort."
         stop
      end if
!      call localorb_info('',use_unit)
    end subroutine get_atomic_occ_numbers
!******

!****s* free_atoms/check_occupation_dimensions
!  NAME
!    check_occupation_dimensions
!  SYNOPSIS
    subroutine check_occupation_dimensions & 
    ( current_n, current_l, current_spin, current_species)
!  PURPOSE
!
!   This subroutine is a hack that checks whether a given atomic shell n, l, is allowed to be 
!   occupied in a certain configuration.
!   This is purely a problem of how to initialize the array dimensions and content of
!   the code correctly. We fix it here only to mitigate any potential damage. A real
!   repair would be possible by adjusting the parse_geometry, parse_control and read_species_data
!   routines to know about additional shells that are not part of the neutral atom.
!
!   For example, if you specify that you would like to have a 
!   Cl atom with -2 as an initial charge, the Cl^-2 ion will have to access the 4s shell.
!   This is unphysical but legal from a technical point of view. However, the Cl atom that the
!   code knows has no 4s shell ... so the dimensions of all our arrays are not prepared to 
!   do this.
!
!   We could fix this, but a proper fix would be quite intricate (I think ultimately we would
!   need to create a new input item that allows the user to specify what exact occupation they
!   mean if the specify an unphysical atom.
!
!   This effort is not worth it, since this is a rare occurrence and there is a simple hack
!   to allow the right shell anyway: Include the desired shell in the 'valence' block of
!   the species definition. For Cl this would be (e.g.):
!
!#     valence basis states
!    valence      4  s   0.00001
!    valence      3  p   4.99999
! 
!   instead of the usual
!
!#     valence basis states
!    valence      3  s   2.000
!    valence      3  p   5.000
!
!  After this you could initialize with a Cl^2- electron since the needed shell is known and was
!  properly set up.
!
!  So if the code encounters such a case, it will simply stop and tell the user what to do.
!
!  If anyone needs this more frequently, please consider putting in the effort to create a 
!  proper way to allow for the exact functionality that you need. 
!
!  USES
      use dimensions,   only : n_spin
      use species_data, only : l_shell_max, valence_n_max
      use localorb_io,  only : use_unit, OL_high, localorb_info
      use mpi_tasks,    only : aims_stop_coll
      implicit none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!     Release version, FHI-aims (2013).
!     Created by VB.
!  ARGUMENTS
      integer :: current_n, current_l, current_spin, current_species
!  INPUTS
!    o current_n - main quantum number of current shell
!    o current_l - angular momentum quantum number of current shell
!    o current_spin - spin channel of current shell
!    o current_species - species number of current atom
!  SOURCE
      logical :: illegal_shell = .false.
      character*120 :: info_str

      if ( (current_n.gt.valence_n_max(current_l,current_species) ) & 
           .or. (current_l.gt.l_shell_max(current_species)) & 
           .or. (current_spin.gt.n_spin) )  then
         illegal_shell = .true.
      end if

      if (illegal_shell) then
         ! Write descriptive output error message and stop the code.

         write (info_str, '(1X,A)') &
           '* Error while trying to find the correct atomic occupation numbers to create initial density.'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A,A)') &
           '* The following shell needs to be occupied but is not ', &
              'allowed by the "valence" tag in control.in: '
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A,I8)') &
           '*   Species number                    : ', current_species
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A,I8)') &
           '*   Shell quantum number n            : ', current_n
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A,I8)') &
           '*   Angular momentum quantum number l : ', current_l
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A,I8)') &
           '*   Spin channel                      : ', current_spin
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* This may be due to a negative ion which would occupy a shell that does not exist '
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* in a free atom. '
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* One example is a Cl^2- ion which needs the 4s shell, although the Cl atom  '
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* only has shells occupied up to n=3 . '
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* A simple remedy is to change the "valence" tag in the species definition in control.in.'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* Example (for Cl): '
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '#     valence basis states'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '    valence      4  s   0.00001'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '    valence      3  p   4.99999'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* This would create the necessary 4 s shell for the Cl^2- ion (even if unphysical).'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* Of course, please be sure to check whether you actually want to do this calculation.'
         call localorb_info(info_str,use_unit,'(A)',OL_high)
         write (info_str, '(1X,A)') &
           '* '
         call localorb_info(info_str,use_unit,'(A)',OL_high)

         ! This problem is input related, so all processors should independently
         ! arrive here. The error message above must also be written, which is
         ! what the aims_stop_coll call hopefully ensures.
         call aims_stop_coll & 
         ('Error when searching for occupation numbers of initialized atoms in free_atoms.f90.', & 
          'check_occupation_dimensions' )

      end if

    end subroutine check_occupation_dimensions

  end module free_atoms
