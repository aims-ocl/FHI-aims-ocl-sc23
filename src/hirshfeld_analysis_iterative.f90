!  NAME
!    get_occ_numbers
!  SYNOPSIS
      subroutine get_occ_numbers (n_electrons, n_spin_1, moment, &
           kind_of_initial, atomic_occ_numbers, write_out)
!  PURPOSE
!    calculate occupation numbers for free atoms
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
      use dimensions
      use grids
      use constants
      use runtime_choices
      use spline
      use species_data
      use localorb_io
      use mpi_tasks, only: aims_stop, check_allocation
      implicit none

!  ARGUMENTS
      real*8 :: n_electrons
      integer :: n_spin_1
      real*8, dimension(n_spin_1) :: n_electrons_polarized
      real*8 :: moment
      integer :: kind_of_initial
      logical :: write_out
      real*8, dimension(1:n_wave_max, 0:l_wave_max, n_spin_1) :: &
           atomic_occ_numbers
!  INPUTS
!    o n_electrons -- number of electrons
!    o number of electrons in each spin channel
!    o moment -- ???????
!    o kind_of_initial -- initialization_type
!  OUTPUTS
!    o atomic_occ_numbers -- final initialized atomic occupation numbers 
!  SOURCE
      ! local variables
      integer :: spin_multiplier
      character*100 :: info_str
      integer :: total_electrons

      ! counter
      integer :: i_spin
      integer :: i_l
      integer :: i_n
      integer :: start_n
      integer :: start_l
      character(*), parameter :: func = 'hirshfeld_analysis_iterative:get_occ_numbers'

      ! Store n_electrons here
      total_electrons = nint(n_electrons)
      if (abs(total_electrons - n_electrons) > 1d-10) then
         call aims_stop('fractional number of electrons', func)
      end if

      atomic_occ_numbers = 0.d0

      ! choose orbitals in the order of the "step"-scheme (?)
      ! which gives roughly the right order of orbitals
      ! if you neglect the exceptions
      ! each orbital is then occupied according to hund's rule

      if (write_out) then
         call localorb_info('',use_unit)
         write(info_str,'(2X,A,F8.4)') "# electrons", n_electrons
         call localorb_info(info_str,use_unit,'(A)')
      end if

      spin_multiplier = mod(n_spin_1, 2) + 1

      if (kind_of_initial .eq. 1) then
         if (write_out) then
            call localorb_info( &
                 "apply hund's rule on atomic states.", &
                 6,'(2X,A)')
         else if (kind_of_initial .eq. 2) then
            write(info_str,'(2X,A,F8.4)') &
                 "desired atomic moment for this type:", moment
            call localorb_info(info_str)
         end if
      end if

      if (write_out) then
         write(info_str,'(2X,A)') &
              "List of occupation numbers for free atom basis: "
         call localorb_info(info_str,use_unit,'(A)')
      end if
      if (kind_of_initial .eq. 1) then
         if (write_out) then
            if (n_spin_1 .eq. 2) then
               write(info_str,'(5X,A,4X,A,3X,A,1X,A)') &
                    "n", "l", "spin up", "spin down"
            else
               write(info_str,'(5X,A,4X,A,4X,A)') &
                    "n", "l", "occ"
            end if
            call localorb_info(info_str)
         end if

         ! just apply hund's rule so that atomic moment is fixed

         start_n = 1
         start_l = 0
         do while (n_electrons .gt. 0)

            i_n = start_n
            do i_l = start_l, 0, -1

               if (n_electrons .ge. 2.d0 * (2*i_l + 1)) then
               ! enough electrons to fill the complete subshell
                  do i_spin = 1, n_spin_1, 1
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          spin_multiplier * (2*i_l + 1)
                  end do
               else
               ! partially filled subshell, so apply hund's rule
                  if (n_spin_1 .eq. 1) then
                  ! no spin-polarization, so just put remaining electrons
                  ! in the last shell
                     i_spin = 1
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          n_electrons
                  else
                  ! hund's rule
                     if (n_electrons .le. (2*i_l + 1)) then
                        i_spin = 1
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons
                     else
                        i_spin = 1
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             (2*i_l + 1)
                        i_spin = 2
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons - (2*i_l + 1)
                     end if
                  end if
               end if

               if (write_out) then
                  if (n_spin_1 .eq. 2) then
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
                     do i_spin = 1, n_spin_1, 1
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             spin_multiplier * (2*i_l + 1)
                     end do
                  else
                  ! partially filled subshell, so apply hund's rule
                     if (n_spin_1 .eq. 1) then
                  ! no spin-polarization, so just put remaining electrons
                  ! in the last shell
                        i_spin = 1
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons
                     else
                  ! hund's rule
                        if (n_electrons .le. (2*i_l + 1)) then
                           i_spin = 1
                           atomic_occ_numbers(i_n, i_l, i_spin) = &
                                n_electrons
                        else
                           i_spin = 1
                           atomic_occ_numbers(i_n, i_l, i_spin) = &
                                (2*i_l + 1)
                           i_spin = 2
                           atomic_occ_numbers(i_n, i_l, i_spin) = &
                                n_electrons - (2*i_l + 1)
                        end if
                     end if
                  end if

                  if (write_out) then
                     if (n_spin_1 .eq. 2) then
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
            if (n_spin_1.eq.1) then
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
            if (n_spin_1.eq.1) then
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
            if (n_spin_1.eq.1) then
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
            if (n_spin_1.eq.1) then
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
            if (n_spin_1.eq.1) then
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

         if (write_out) then
            write(info_str,'(2X,A,F8.4,F8.4)') "# electrons polarized", &
                 n_electrons_polarized(1), n_electrons_polarized(2)
            call localorb_info( info_str )
            write(info_str,'(5X,A,4X,A,3X,A,1X,A)') &
                 "n", "l", "spin up", "spin down"
            call localorb_info( info_str )
         end if
         do i_spin = 1, n_spin_1, 1

            start_n = 1
            start_l = 0

            do while (n_electrons_polarized(i_spin) .gt. 0)

               i_n = start_n
               do i_l = start_l, 0, -1

                  if (n_electrons_polarized(i_spin) .ge. (2*i_l + 1)) &
                       then
                     ! enough electrons to fill the complete subshell
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          (2*i_l + 1)
                  else
                     atomic_occ_numbers(i_n, i_l, i_spin) = &
                          n_electrons_polarized(i_spin)
                  end if

!                  write(use_unit,*) i_spin, n_electrons_polarized(i_spin),
!     +                 moment
                  if (((i_spin .eq. 2) .or. &
                       ((i_spin .eq. 1) .and. &
                       (n_electrons_polarized(i_spin) .le. moment))) &
                       .and. (write_out)) then
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
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             (2*i_l + 1)
                     else
                        atomic_occ_numbers(i_n, i_l, i_spin) = &
                             n_electrons_polarized(i_spin)
                     end if

                     if (((i_spin .eq. 2) .or. &
                          ((i_spin .eq. 1) .and. &
                          (n_electrons_polarized(i_spin) .le. moment))) &
                          .and. (write_out)) then
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
              if (n_spin_1.eq.2) then
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
              if (n_spin_1.eq.2) then
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
              if (n_spin_1.eq.2) then
                atomic_occ_numbers(6, 0, 1) =  1.
                atomic_occ_numbers(5, 2, 1) =  1.
                atomic_occ_numbers(4, 3, 1) =  1.
              end if
            else if (moment.eq.0.) then
              write(info_str,'(2X,A)') "* Attention: Adjustment of initial occupation."
              call localorb_info( info_str )
              write(info_str,'(2X,A)') "* Found a Ce-like atom (58 electrons) - enforcing 6s^2 5d^1 4f^1 (unpolarized!)."
              call localorb_info( info_str )
              if (n_spin_1.eq.2) then
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
              if (n_spin_1.eq.2) then
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
              if (n_spin_1.eq.2) then
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
    end subroutine get_occ_numbers



!****s* FHI-aims/get_ion_hirsh
!  NAME
!   get_ion_hirsh
!  SYNOPSIS

subroutine get_ion_hirsh &
(i_atom, charge_ion )

!  PURPOSE
!  Subroutine get_ion_hirsh provides splined ionic density for each species
!
!  this used to be a wrapper around Martin Fuchs' sr_atom.f
!  It is now a procedure that calls atomic solvers, but it still bears the 
!  marks of its original purpose.
!
!  We use a non-relativistic version for now. For a proper relativistic treatment, 
!  must have a radial relativistic Hamiltonian which can be extended directly to 
!  the full (non-radial-symmetric) potential. This is not the case at present.
!
!  USES

      use dimensions
      use runtime_choices
      use grids
      use species_data
      use free_atoms
      use spline
      use localorb_io
      use constants
      use synchronize_mpi, only: sync_atomic_solver_outputs
      use geometry,        only: species
      use psi_at_nucleus_mod, only: psi_at_nucleus_atomic
      use mpi_tasks, only: aims_stop, check_allocation
      implicit none

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
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
!
!  imported variables
   integer :: i_atom
   real*8 :: charge_ion

!  local variables

      real*8 drho_dr(n_max_grid, n_species)
      real*8 d2rho_dr2(n_max_grid, n_species)

!  variables for sratom - some of these should be set in control.in

!  sc_max_iter : Maximum number of self-consistency iterations for free atom
!  i_rel    : Masked version of flag_rel; i_rel = 1 <-> non-rel., i_rel=2 <-> Koelling-Harmon
!  i_grad   : if i_grad=1, calculate free-atom density derivatives (1st and 2nd)
!  i_exch   : Determines the exchange-correlation handling for sr_atom.
!             * i_exch > 0 : unpolarised, i_exch < 0 : spin-polarised version
!             * many implemented types of behavior / functionals, see fhipp for complete list
!             * relevant here: iexch = +-6 is PBE, iexch = +-3 is LDA
!  t_core   : Here, always true - we want to calculate core wavefunctions!
!  nuclear_charge : well, the nuclear charge ...
!  n_core   : number of frozen core shells in a pseudopot'l calculation?? This is not needed here!
!  n_valence: number of "valence" (n,l) shells in the calculation
!  m_shell  : This value seems to label the nature of each state - spin up (m=1) or spin down (m=2)
!             for each shell
!  occ_shell : occupation number of each shell
!  eigenval  : energy eigenvalues for each shell
!  i_grid_max: integer index of maximum grid point in integration for each shell
!  d_eig_sum : FIXME This is an apparently unnecessary variable which monitors the
!              change of the eigenvalue sum between successive scf iterations
!  r_peak    : FIXME: Unused. Outermost peak position of each individual wave function.
!  ini_pot :   initial radial potential for each individual atom. 
!  atom_pot :  s.-c. atomic potential FIXME: THIS HAS ONE INDEX TOO MANY - WE DO NOT NEED THE 
!              max_shells INDEX FOR ANYTHING, IT IS A HARTREE_FOCK RELIC!
!  atom_pot_es: electrostatic portion of s.-c. atomic potential FIXME: THIS HAS ONE INDEX TOO MANY 
!              - WE DO NOT NEED THE max_shells INDEX FOR ANYTHING, IT IS A HARTREE_FOCK RELIC!
!  density, density_deriv, density_2nd_deriv : radial atomic density and its derivatives.
!              FIXME: SAME AS atom_pot - THE max_shells INDEX IS !OMPLETELY USELESS!
!  core_density, core_dd, core_d2d : UNUSED - core density, density derivative, and second 
!              derivative, for frozen core calculations (never done here)
!  wave :      s.-c. atomic wave functions     
!  t_in_pot : if true, sratom() assumes an input non-selfconsistent potential - always false for now.
!  t_sic :    if true, sratom() does self-interaction correction - this is never the case. 
!  e_sic :    array of self-energies for SI! - not used
!  v_orb :    orbital-depedent potential(?) for sic - not used
!  t_kli :    kli is not used 
!  sv_kli :   kli is not used
!  svm   :    kli is not used
!  t_mgga :   meta-gga used? currently never.
!  wave_deriv : radial derivative of each wave function - might be useful after all.

      integer sc_max_iter
      integer i_rel
      integer i_grad
      integer i_exch
      logical t_core
      logical found
      integer n_max_ind_fns1      
      real*8 cut_free_atom_fixed

      real*8 nuclear_charge
      real*8 ::  dr_coef, alpha
      real*8 ::  residual_charge

      integer n_core
      integer n_valence
      integer i_fn
      integer i_l
      integer i_shell
      integer i_charge
      integer n
      integer l
      integer temp
      integer max_shell
      integer,save :: n_spin_1 = 2

      integer, dimension(:), allocatable :: n_shell
      integer, dimension(:), allocatable :: l_shell
      integer, dimension(:), allocatable :: n_shell1
      integer, dimension(:), allocatable :: l_shell1
      integer, dimension(:), allocatable :: m_shell
      integer, dimension(:), allocatable :: perm
      integer, dimension(:), allocatable :: inv_perm
      real*8, dimension(:), allocatable :: occ_shell
      real*8, dimension(:), allocatable :: occ_shell1
      real*8, dimension(:), allocatable :: eigenval
      integer, dimension(:), allocatable :: i_grid_max
      logical, dimension(:), allocatable :: core_type
      real*8, dimension(1:n_wave_max, 0:l_wave_max, 2) :: &
           atomic_occ

      real*8 d_eig_sum

      real*8, dimension(:), allocatable :: r_peak

      real*8, dimension(:), allocatable :: ini_pot
      real*8, dimension(:,:), allocatable :: atom_pot
      real*8, dimension(:,:), allocatable :: atom_pot_es
      real*8, dimension(:,:), allocatable :: atom_pot_xc
      real*8, dimension(:,:), allocatable :: density
      real*8, dimension(:,:), allocatable :: density_deriv
      real*8, dimension(:,:), allocatable :: density_2nd_deriv
      real*8, dimension(:,:), allocatable :: core_density
      real*8, dimension(:,:), allocatable :: core_dd
      real*8, dimension(:,:), allocatable :: core_d2d

      real*8, dimension(:,:), allocatable :: wave
      real*8, dimension(:,:), allocatable :: wave_deriv

      logical t_in_pot 
      logical t_sic
      real*8, dimension(:), allocatable :: e_sic
      real*8, dimension(:,:), allocatable :: v_orb

      logical t_kli
      real*8, dimension(:), allocatable :: sv_kli
      real*8, dimension(:), allocatable :: svm

      logical t_mgga

!  other aux variables

      integer n_max

      character*100 :: info_str

      real*8, dimension(:,:), allocatable :: kinetic

!  counters
      integer i_species
      integer i_grid
      integer i_channels
      integer i_function, info

!  begin work

      cut_free_atom_fixed = 10.d0
      i_species = species(i_atom)

!      write (info_str,'(2X,A,A)') &
!        "Creating wave function, potential, and density ", &
!        "for ions in Hirshfeld analysis."
!      call localorb_info('',use_unit)
!      call localorb_info(info_str,use_unit,'(A)')

!      write(use_unit,*) 
!      write(use_unit,*) " !reating wave function, potential, and density ",
!     +  "Fo free atoms."

!  allocations - needed because some of these arrays may collide with the
!  stack size on small machines if they were automatic arrays

      n_max_ind_fns1 = n_max_ind_fns + 3 

      allocate(n_shell(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'n_shell                       ') 

      allocate(n_shell1(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'n_shell1                       ')

      allocate(l_shell(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'l_shell                       ') 

      allocate(l_shell1(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'l_shell1                       ')

      allocate(m_shell(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'm_shell                       ') 

      allocate(perm(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'perm                       ')

      allocate(inv_perm(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'inv_perm                       ')

      allocate(occ_shell(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'occ_shell                     ') 

      allocate(occ_shell1(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'occ_shell1                     ')

      allocate(eigenval(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'eigenval                      ') 

      allocate(i_grid_max(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'i_grid_max                    ') 

      allocate(r_peak(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'r_peak                        ') 

      allocate(core_type(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'core_type                     ') 

      allocate(ini_pot(n_max_grid),stat=info)
      call check_allocation(info, 'ini_pot                       ') 

      allocate(atom_pot(n_max_grid, n_spin_1),stat=info)
      call check_allocation(info, 'atom_pot                      ') 

      allocate(atom_pot_es(n_max_grid, n_spin_1),stat=info)
      call check_allocation(info, 'atom_pot_es                   ') 

      allocate(atom_pot_xc(n_max_grid, n_spin_1),stat=info)
      call check_allocation(info, 'atom_pot_xc                   ') 

      allocate(density(n_max_grid, n_spin_1),stat=info)
      call check_allocation(info, 'density                       ') 

      allocate(density_deriv(n_max_grid, n_spin_1),stat=info)
      call check_allocation(info, 'density_deriv                 ') 

      allocate(density_2nd_deriv(n_max_grid, n_spin_1),stat=info)
      call check_allocation(info, 'density_2nd_deriv             ') 

      allocate(core_density(n_max_grid,n_spin_1),stat=info)
      call check_allocation(info, 'core_density                  ') 

      allocate(core_dd(n_max_grid,n_spin_1),stat=info)
      call check_allocation(info, 'core_dd                       ') 

      allocate(core_d2d(n_max_grid,n_spin_1),stat=info)
      call check_allocation(info, 'core_d2d                      ') 
                                                                                       
      allocate(wave(n_max_grid, n_max_ind_fns1),stat=info)
      call check_allocation(info, 'wave                          ') 

      allocate(e_sic(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'e_sic                         ') 

      allocate(v_orb(n_max_grid,n_max_ind_fns1),stat=info)
      call check_allocation(info, 'v_orb                         ') 

      allocate(sv_kli(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'sv_kli                        ') 

      allocate(svm(n_max_ind_fns1),stat=info)
      call check_allocation(info, 'svm                           ') 

      allocate(wave_deriv(n_max_grid, n_max_ind_fns1),stat=info)
      call check_allocation(info, 'wave_deriv                    ') 

      allocate(kinetic(n_max_grid, n_max_ind_fns1), stat=info)
      call check_allocation(info, 'kinetic                       ') 

!  input parameters for sr_atom(), from control.in 

      ! blanket setting for X!, if basis fn X! is the same as eventual self-consistent XC
      i_exch = flag_xc

!  for Hartree-Fock calculation --Ren
      if(flag_xc .eq. 0) then
        i_exch = 8 
      elseif(flag_xc .eq. 1) then
!  Hybrid-PBE0 calculation, using PBE instead as the atomic solver.
        i_exch = 6
      endif
      if(flag_xc .eq. 10) then
!  Hybrid-B3LYP calculation, using BLYP instead as the atomic solver.
        i_exch = 9
      endif
! HSE calculation or lc-wpbeh, using PBE
      if (flag_xc.eq.7.or.flag_xc.eq.23) then
        i_exch = 6
      endif
! PBEsol0 hybrid, using PBE
      if (flag_xc.eq.13) then
        i_exch = 6
      endif
!revPBE  
      if (flag_xc.eq.12) then
        i_exch = 15
      endif
!revPBE_vdw.  SAG
      if (flag_xc.eq.14) then
         i_exch = 15
      endif

!RPBE
      if (flag_xc.eq.11) then
        i_exch = 14
      endif

! PBEsol calculation, using PBE instead as the atomic solver
      if (flag_xc.eq.17) then
         i_exch = 6
      end if
!PBE_vdw   SAG
      if (flag_xc.eq.5) then
         i_exch = 6
      end if

! FIXME : A case should be added in external/vexcor.f in order to
! calculate the AM05 energies and potentials (using subroutine
! external/am05.f) and use that as atomic solver!
! AM05 calculation, using PW-LDA instead as the atomic solver
      if (flag_xc.eq.20) then
         i_exch = 8 
      end if
! VWN-LDA approximation, using VWN as initial guess.
! Unfortunately it is number 7 in vexcor.f
! I hope there will be no mixing with HSE
      if (flag_xc.eq.15) then
        i_exch = 7
      endif
! LDA's, no gradient derivatives needed   
      if ((i_exch.eq.3).or.(i_exch.eq.8).or.(i_exch.eq.7)) then
        i_grad = 0
!     PBE and its relatives or BLYP, should calculate the derivatives for free-atom
      else if ((i_exch.eq.6).or.(i_exch.eq.9).or.(i_exch.eq.14).or.(i_exch.eq.15)) then
        i_grad = 1
      end if

      if (use_density_gradient) then
        i_grad = 1
      end if

      if (flag_rel.eq.0) then
        i_rel=2
      else if ((flag_rel.eq.REL_zora).or.(flag_rel.eq.REL_atomic_zora)) &
            then
!TJ Changed to ZORA
        i_rel=9  
      else if (flag_rel.eq.REL_own.or.flag_rel==REL_KOLNING_HARMON)then
        i_rel = 10
      end if

!  more input parameters for sr_atom()
!  FIXME: Some of these should be adjustable in control.in !

      sc_max_iter = 100
      t_core = .true.
      n_core = 0
      t_in_pot = .false.
      t_sic = .false.
      t_mgga = .false.

!     We use KLI for free-atom densities
      t_kli = .true.
      i_exch = 12

!      t_kli = .false.

   if (n_spin_1 .eq. 2) i_exch=-i_exch


!  initialize grids
!  initialize also the core density parts, currently unused.
!  FIXME: This should really be done f90 style!

        do i_grid = 1, n_max_grid, 1
          free_potential(i_grid, i_species)   = 0.0d0          
          free_pot_es   (i_grid, i_species)   = 0.0d0          
          free_rho(i_grid, i_species) = 0.0d0          
          drho_dr(i_grid, i_species) = 0.0d0          
          d2rho_dr2(i_grid, i_species) = 0.0d0          
          do i_channels = 1,n_spin_1, 1
            core_density(i_grid,1) = 0.d0
            core_dd(i_grid,1) = 0.d0
            core_d2d(i_grid,1) = 0.d0
            density(i_grid,1) = 0.d0 
          enddo
!          do i_function = 1, n_max_ind_fns1, 1
!            free_wave(i_grid,i_species,i_function) = 0.d0
!          enddo 
        enddo

        write(info_str,'(2X,A,A,A,F15.8)') "Species: ",&
         species_name(i_species),"Charge:", charge_ion
        call localorb_info('',use_unit)
        call localorb_info(info_str,use_unit,'(A)')

        nuclear_charge = species_z(i_species)

!  filled states only
        n_valence = n_atomic(i_species)
        


        write(info_str,'(2X,A,A,A,I4,I4)') "Species,n_max_ind_fns1: ",&
         species_name(i_species),"N_valence:", n_valence, n_max_ind_fns1 
        call localorb_info('',use_unit)
        call localorb_info(info_str,use_unit,'(A)')

      call get_occ_numbers(nuclear_charge-charge_ion,n_spin_1,0.d0,1,atomic_occ,.true.)
!       call get_atomic_occ_numbers(nuclear_charge-charge_ion,0.d0,1,atomic_occ,.true.)
    
       n=1 
       l=0
       do i_function = 1, n_valence+3, 1
            n_shell (i_function) = n 
            l_shell (i_function) = l
            m_shell (i_function) = 1
            if (i_function.le.n_valence+3) then
               occ_shell(i_function) = atomic_occ(n_shell(i_function),l_shell(i_function),1)
            else 
              occ_shell(i_function)=0.0
            endif

            if (i_function.le.n_valence) then
              core_type(i_function) = core_fn(i_species,i_function)
            else 
              core_type(i_function) = .false.
            end if
            l=l+1
            if (l.gt.n-1) then
             n=n+1
             l=0
            endif
        enddo

!     Sort the states moving 0--occupied ones down 
      do i_fn = n_valence+3, 2, -1
       if ( (occ_shell(i_fn) .gt. 0.d0) .and. (occ_shell(i_fn-1) .eq. 0.d0) ) then
        n_shell1(1)=n_shell(i_fn)
        l_shell1(1)=l_shell(i_fn)
        occ_shell1(1)=occ_shell(i_fn)
        n_shell(i_fn)=n_shell(i_fn-1)
        l_shell(i_fn)=l_shell(i_fn-1)
        occ_shell(i_fn)=occ_shell(i_fn-1)
        n_shell(i_fn-1)=n_shell1(1)
        l_shell(i_fn-1)=l_shell1(1)
        occ_shell(i_fn-1)=occ_shell1(1)
        max_shell = i_fn-1
       end if 
      enddo


!     output free ion data
      call localorb_info('',use_unit)

      write(info_str,'(2X,A)') &
        "List of free ionic orbitals and eigenvalues: "
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(4X,A,4X,A,3X,A,6X,A,4X,A)') &
        "n", "l", "occ", "energy [Ha]", &
        "energy [eV]"
      call localorb_info(info_str,use_unit,'(A)')

      do i_fn = 1, n_valence, 1

        write(info_str,'(2X,I3,2X,I3,2X,F5.2,F15.6,F15.4)') &
          n_shell(i_fn), l_shell(i_fn), occ_shell(i_fn), &
          eigenval(i_fn), &
          eigenval(i_fn)*hartree
        call localorb_info(info_str,use_unit,'(A)')

      enddo
      call localorb_info('',use_unit)

! end output



        if (i_exch.lt.0) then
          write(info_str,'(1X,A,A)') & 
            "* m_shell not set properly for ", &
            "spin-polarized case!"
          call localorb_info(info_str)
          stop
        end if


        call atomini &
        (n_valence, n_grid(i_species), nuclear_charge, n_shell, &
         occ_shell, r_grid(1,i_species), ini_pot, eigenval &
        )

!       FIXME: This would be a superfluous copy, if atom_pot were not grossly overdimensioned
        do i_grid = 1, n_grid(i_species), 1
          atom_pot(i_grid, 1) = ini_pot(i_grid)
          atom_pot_xc(i_grid, 1) = 0.d0
          atom_pot_es(i_grid, 1) = 0.d0
        enddo

        ! Here we run the chosen atomic solver to generate candidate basis functions
        ! Should you wish to add in another solver, the values extracted from these 
        ! solvers are:
        !      eigenval, wave, atom_pot, atom_pot_es, density, density_deriv, density_2nd_deriv
        ! defined on the radial grid specified by r_grid (note that the l_channel
        ! index in these variables is legacy code that is not needed)
        if (atomic_solver .eq. ATOMIC_SOLVER_SRATOM) then
          call sratom &
          ( sc_max_iter, i_rel, i_grad, i_exch, t_core, nuclear_charge, &
            n_core, n_valence, n_shell, l_shell, m_shell, occ_shell, &
            eigenval, i_grid_max, d_eig_sum, &
            n_max_grid, n_max_ind_fns1, n_spin_1, &
            n_grid(i_species), r_grid(1,i_species),  &
!            .false.,cutoff_type(i_species), &
            .true., cutoff_type(i_species), &
            cut_free_atom_fixed,&
            scale_cutoff(i_species), w_cutoff(i_species), &
            r_peak, atom_pot, atom_pot_es,&
            atom_pot_xc, density, density_deriv, &
            density_2nd_deriv, core_density, core_dd, core_d2d, &
            wave, t_in_pot, t_sic, e_sic, v_orb, t_kli, sv_kli,&
            svm, t_mgga, wave_deriv, core_type &
          )
        else if (atomic_solver .eq. ATOMIC_SOLVER_ATOM_SPHERE) then
          ! atom_sphere_wrapper's argument list has been deliberately constructed to be an example of a generic 
          ! interface between aims and a radial atomic solver
          ! No modules are used within this wrapper, other than for system-level tasks like aims_stop.
          ! What you see here is an exhaustive list of what properties of the atom are needed from aims, what variables
          ! need to be exported back to aims to construct the basis sets, and nothing more (e.g. all work matrices reside 
          ! within the wrapper.) 

          call atom_sphere_wrapper &
          ( n_grid(i_species), r_grid(1,i_species), species_name(i_species), flag_xc, hybrid_coeff, hse_omega, &
            nuclear_charge, n_core, n_valence, n_max_ind_fns1, n_shell, l_shell, m_shell, occ_shell, n_max_grid, & 
            n_channels, n_spin_1, &
            .true., cutoff_type(i_species), cut_free_atom_fixed, scale_cutoff(i_species), &
                 w_cutoff(i_species), & ! Cutoff potential
            eigenval, wave, wave_deriv, density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
            kinetic, psi_at_nucleus_atomic(:n_valence,i_species)) ! Output
          call sync_atomic_solver_outputs( n_max_ind_fns1, n_max_grid, n_channels, n_spin_1, &
               eigenval, wave, wave_deriv, kinetic, &
               density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es )
        else
          call aims_stop('The atomic solver specified by the variable atomic_solver was not recognized, exiting.','func')  
        end if

!     Sort the states according to their energies
!      call insertionsort(eigenval, n_valence, perm, inv_perm)
!      do i_fn = 1, n_valence, 1
!        n_shell1(perm(i_fn)) = n_shell(i_fn)
!        l_shell1(perm(i_fn)) = l_shell(i_fn)
!        occ_shell1(perm(i_fn)) = occ_shell(i_fn)
!      enddo
!      n_shell(:) = n_shell1(:)
!      l_shell(:) = l_shell1(:)
!      occ_shell(:) = occ_shell1(:)

!        write(info_str,'(2X,A)') "after sort ..."
!        call localorb_info('',use_unit)
!        call localorb_info(info_str,use_unit,'(A)')



!     Distribute charge over atomic states with highest energy 
!      residual_charge = charge_ion

!     First case -- positive charge
!      if (residual_charge.gt.0.0) then
!         occ_shell(n_valence) = occ_shell(n_valence) - residual_charge
!         if (occ_shell(n_valence) .lt. 0.d0) occ_shell(n_valence) = 0.0 

!        found=.false.
!        do i_fn = 1, n_valence, 1
!          if ((occ_shell(i_fn).eq.0.0) .and. (.not. found)) then
!           i_charge=i_fn-1
!           found = .true.
!          endif 
!        enddo 
!        do while (residual_charge .gt. 0.0) 
!          occ_shell(i_charge) = occ_shell(i_charge) - residual_charge
!          if (occ_shell(i_charge).ge.0.0) then
!            residual_charge = 0.0
!          else
!            residual_charge = 0.0     ! Forbid desocupying lower-lying electron levels
!            Uncomment below two lines for allowing occupations of lower-lying
!            electron levels (and comment the "residual_charge = 0.0" line)
!            residual_charge = -occ_shell(i_charge)
!            occ_shell(i_charge) = 0.0
!            i_charge = i_charge - 1
!          end if
!        enddo         
!      endif


!     Second case -- negative charge
!      if (residual_charge.lt.0.0) then
!        do i_fn = 1, n_valence, 1
!          if (occ_shell(i_fn).lt.((1+2*l_shell(i_fn))*2.0) .and. residual_charge.lt.0.0) then
!            occ_shell(i_fn)=occ_shell(i_fn)-residual_charge

!            residual_charge=0.0
!            if (occ_shell(i_fn).gt.((1+2*l_shell(i_fn))*2.0)) then 
!               occ_shell(i_fn)=(1+2*l_shell(i_fn))*2.0
!            endif

!            residual_charge=((1+2*l_shell(i_fn))*2.0)-occ_shell(i_fn)
!            if (occ_shell(i_fn).gt.((1+2*l_shell(i_fn))*2.0)) then
!              occ_shell(i_fn)=(1+2*l_shell(i_fn))*2.0
!            endif
!          endif
!        enddo
!      endif

!     output free ion data
      call localorb_info('',use_unit)

      write(info_str,'(2X,A)') &
        "List of free ionic orbitals and eigenvalues: "
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(4X,A,4X,A,3X,A,6X,A,4X,A)') &
        "n", "l", "occ", "energy [Ha]", &
        "energy [eV]"
      call localorb_info(info_str,use_unit,'(A)')

      do i_fn = 1, n_valence, 1

        write(info_str,'(2X,I3,2X,I3,2X,F5.2,F15.6,F15.4)') &
          n_shell(i_fn), l_shell(i_fn), occ_shell(i_fn), &
          eigenval(i_fn), &
          eigenval(i_fn)*hartree
        call localorb_info(info_str,use_unit,'(A)')

      enddo
      call localorb_info('',use_unit)

! end output



        sc_max_iter = 100
        do i_grid = 1, n_max_grid, 1
          free_potential(i_grid, i_species)   = 0.0d0
          free_pot_es   (i_grid, i_species)   = 0.0d0
          free_rho(i_grid, i_species) = 0.0d0
          drho_dr(i_grid, i_species) = 0.0d0
          d2rho_dr2(i_grid, i_species) = 0.0d0
          do i_channels = 1,n_spin_1, 1
            core_density(i_grid,1) = 0.d0
            core_dd(i_grid,1) = 0.d0
            core_d2d(i_grid,1) = 0.d0
            density(i_grid,1) = 0.d0
          enddo
!          do i_function = 1, n_max_ind_fns1, 1
!            free_wave(i_grid,i_species,i_function) = 0.d0
!          enddo
        enddo



        call atomini &
        (n_valence, n_grid(i_species), nuclear_charge, n_shell, &
         occ_shell, r_grid(1,i_species), ini_pot, eigenval &
        )

!       FIXME: This would be a superfluous copy, if atom_pot were not grossly
!       overdimensioned
        do i_grid = 1, n_grid(i_species), 1
          atom_pot(i_grid, 1) = ini_pot(i_grid)
          atom_pot_xc(i_grid, 1) = 0.d0
          atom_pot_es(i_grid, 1) = 0.d0
        enddo


        if (atomic_solver .eq. ATOMIC_SOLVER_SRATOM) then
          call sratom &
          ( sc_max_iter, i_rel, i_grad, i_exch, t_core, nuclear_charge, &
            n_core, n_valence, n_shell, l_shell, m_shell, occ_shell, &
            eigenval, i_grid_max, d_eig_sum, &
            n_max_grid, n_max_ind_fns, n_spin_1, &
            n_grid(i_species), r_grid(1,i_species),  &
!            .false.,cutoff_type(i_species), &
            .true., cutoff_type(i_species), &
            cut_free_atom_fixed,&
            scale_cutoff(i_species), w_cutoff(i_species), &
            r_peak, atom_pot, atom_pot_es,&
            atom_pot_xc, density, density_deriv, &
            density_2nd_deriv, core_density, core_dd, core_d2d, &
            wave, t_in_pot, t_sic, e_sic, v_orb, t_kli, sv_kli,&
            svm, t_mgga, wave_deriv, core_type &
          )
        else if (atomic_solver .eq. ATOMIC_SOLVER_ATOM_SPHERE) then
          write(info_str, '(2X,A)') "You have chosen the atom_sphere atomic solver combined with Hirshfeld analysis."
          call localorb_info( info_str )
          write(info_str, '(2X,A)') "This functionality has not been tested, so we're stopping the code for now."
          call localorb_info( info_str )
          write(info_str, '(2X,A)') "Best guess:  it will work for LDA/GGA functionals, but will require some code &
               &changes for hybrids."
          call localorb_info( info_str )
          call aims_stop("Feel free to comment out this stop, but you're on your own.")
          call atom_sphere_wrapper &
          ( n_grid(i_species), r_grid(1,i_species), species_name(i_species), flag_xc, hybrid_coeff, hse_omega, &
            nuclear_charge, n_core, n_valence, n_max_ind_fns, n_shell, l_shell, m_shell, occ_shell, n_max_grid, & 
            n_channels, n_spin_1, &
            .true., cutoff_type(i_species), cut_free_atom_fixed, scale_cutoff(i_species), &
                 w_cutoff(i_species), & ! Cutoff potential
            eigenval, wave, wave_deriv, density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
            kinetic, psi_at_nucleus_atomic(:n_valence,i_species)) ! Output
          call sync_atomic_solver_outputs( n_max_ind_fns, n_max_grid, n_channels, n_spin_1, &
               eigenval, wave, wave_deriv, kinetic, &
               density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es )
        else
          call aims_stop('The atomic solver specified by the variable atomic_solver was not recognized, exiting.','func')  
        end if

!     output free ion data
      call localorb_info('',use_unit)

      write(info_str,'(2X,A)') &
        "List of free ionic orbitals and eigenvalues: "
      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(4X,A,4X,A,3X,A,6X,A,4X,A)') &
        "n", "l", "occ", "energy [Ha]", &
        "energy [eV]"
      call localorb_info(info_str,use_unit,'(A)')

      do i_fn = 1, n_valence, 1

        write(info_str,'(2X,I3,2X,I3,2X,F5.2,F15.6,F15.4)') &
          n_shell(i_fn), l_shell(i_fn), occ_shell(i_fn), &
          eigenval(i_fn), &
          eigenval(i_fn)*hartree
        call localorb_info(info_str,use_unit,'(A)')

      enddo
      call localorb_info('',use_unit)

! end output





!       atom_pot and density have an extra index which used to be shell-dependent (Hatree-Fock?),
!       then indicated spin up / down (Martin Fuchs), and is unused here. 

        do i_grid = 1, n_grid(i_species), 1
          free_pot_es(i_grid, i_species) = atom_pot_es(i_grid, 1)
          free_rho(i_grid, i_species) = density(i_grid, 1)
        enddo

!       spline rho
        call cubic_spline &
        ( free_rho(1,i_species), n_grid(i_species),  &
          free_rho_spl(1,1,i_species) )
        
!     deallocate local variables
      deallocate(wave_deriv)
      deallocate(svm)
      deallocate(sv_kli)
      deallocate(v_orb)
      deallocate(e_sic)
      deallocate(wave)
      deallocate(core_d2d)
      deallocate(core_dd)
      deallocate(core_density)
      deallocate(density_2nd_deriv)
      deallocate(density_deriv)
      deallocate(density)
      deallocate(atom_pot_xc)
      deallocate(atom_pot)
      deallocate(ini_pot)
      deallocate(core_type)
      deallocate(r_peak)
      deallocate(i_grid_max)
      deallocate(eigenval)
      deallocate(occ_shell)
      deallocate(occ_shell1)
      deallocate(perm)
      deallocate(inv_perm)
      deallocate(m_shell)
      deallocate(l_shell)
      deallocate(l_shell1)
      deallocate(n_shell)
      deallocate(n_shell1)
      return
end subroutine



!******
!****s* FHI-aims/hirshfeld_analysis_iterative
!  NAME
!   hirshfeld_analysis_iterative
!  SYNOPSIS

subroutine hirshfeld_analysis_iterative ( )

!  PURPOSE
! Subroutine hirshfeld_analysis
!
! computes a Hirshfeld-type charge-fragment analysis as defined in
!
! Theor. Chim. Acta (Berl.) 44, 129-138 (1977)
!
! For a completely unique prescription, consider also "Hirshfeld-I",
! where the partial charges of the fragments used in the weighting function
! are iterated until they are equal to the partial charges that emerge from
! the Hirshfeld analysis.
!
! The output written by the present routine are partial charges, dipole moments, 
! and quadrupole moments on individual atoms.
! 
!
! AT -- includes Hirshfeld volume partitioning
!
!  USES

  use dimensions
  use runtime_choices
  use constants
  use grids
  use geometry
  use pbc_lists
  use free_atoms
  use physics
  use spline
  use vdw_correction
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use species_data, only: species_name

  implicit none

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
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


  ! local variables

  character*100 :: info_str
  logical converged 

  real*8, dimension(:), allocatable   :: hirshfeld_charge_iter
  real*8, dimension(:), allocatable   :: prev_hirshfeld_charge
  real*8, dimension(:), allocatable   :: free_int 
  real*8, dimension(:), allocatable   :: hirsh_int
  real*8, dimension(:), allocatable   :: hirshfeld_spin_moment
  real*8, dimension(:,:), allocatable :: hirshfeld_dipole
  real*8, dimension(:,:), allocatable :: hirshfeld_quadrupole

  real*8, dimension(:,:,:), allocatable :: reference_rho_spl
  real*8, dimension(:,:,:), allocatable :: reference_rho_spl1

  real*8, dimension(:,:,:,:), allocatable :: rho_spl_ions
  logical, dimension(:,:), allocatable :: rho_spl_ions_exists

  real*8, dimension(3) :: coord_current

  real*8, dimension(n_species) :: r_grid_min_sq

  real*8 dist_tab_sq(n_centers_basis_integrals)
  real*8 dist_tab(n_centers_basis_integrals)
  real*8 dir_tab(3,n_centers_basis_integrals)
  real*8 dir_tab_norm(3,n_centers_basis_integrals)
  real*8 i_r(n_centers_basis_integrals)


  real*8 :: aux_dens
  real*8 :: aux_dens1
  real*8 :: full_dens
  real*8 :: partition_norm
  real*8 :: free_partition_norm
  real*8 :: hirshfeld_partition

  real*8 :: current_i_r
  real*8 :: temp_rho_new
  real*8 deformation_density(n_spin)

  real*8 :: dipole_moment

  real*8 charge_l (n_atoms)
  real*8 charge_u (n_atoms)

  logical :: point_on_atom

  ! counters

  integer :: i_atom

  integer :: i_my_batch
  integer :: i_full_points
  integer :: i_index

  integer :: i_center_L
  integer :: i_center

  integer :: i_spin
  integer :: i_coord
  integer :: i_coord_2

  integer :: i_sec_moment

  integer :: i_iter

  integer :: output_priority_old

  ! shorthand for identification of current grid point in batches
  integer :: current_atom, current_radial, current_angular


  ! begin work

  ! If output was requested under all circumstances, must reduce the
  ! output priority for this subroutine. Must reset at end of subroutine!
  if (out_hirshfeld_always) then
     output_priority_old = output_priority
     if (output_level.eq.'MD_light') output_priority = 1
  end if

  write(info_str,'(2X,A)') "Performing Hirshfeld-I analysis of fragment charges and moments."
  call localorb_info( info_str, use_unit,'(A)',OL_norm)

  ! allocate needed arrays

  allocate ( hirshfeld_charge_iter(n_atoms) )
  allocate ( prev_hirshfeld_charge(n_atoms) )
  allocate ( free_int(n_atoms) )
  allocate ( hirsh_int(n_atoms) )
  allocate ( hirshfeld_dipole(3,n_atoms) )
  allocate ( hirshfeld_quadrupole(6,n_atoms) )
  if (spin_treatment .eq. 1) then
	allocate ( hirshfeld_spin_moment(n_atoms) )
	hirshfeld_spin_moment = 0.d0
  endif

  hirshfeld_charge_iter     = 0.d0
  prev_hirshfeld_charge = 0.d0
  converged = .false. 
  i_iter = 0

  r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)

  ! Initialize charges for the Hirshfeld scheme. In preparation for a possible
  ! self-consistent Hirshfeld-I scheme later, we use temporary arrays for the
  ! "stockholder" atomic fragment densities and do not use evaluate_partition ().

  ! reference_rho_spl is not distributed across CPU's but should (and can) be distributed
  ! in exactly the same fashion as delta_v_hartree_part_spl in update_hartree_potential_p1, if needed.

  if(.not. use_distributed_spline_storage) then
    allocate (reference_rho_spl(n_max_spline,n_max_grid, n_atoms))
    allocate (reference_rho_spl1(n_max_spline,n_max_grid, n_atoms))
    allocate (rho_spl_ions(n_max_spline,n_max_grid,n_species,9)) 
    allocate (rho_spl_ions_exists(n_species,9))
  end if
   

  rho_spl_ions_exists(:,:) = .false.

  do while (converged .eqv. .false.)

  converged = .false.
  i_iter = i_iter + 1

  free_int             = 0.d0
  hirsh_int            = 0.d0
  hirshfeld_dipole     = 0.d0
  hirshfeld_quadrupole = 0.d0
  hirshfeld_charge_iter     = 0.d0

  ! reference_rho_spl is not distributed across CPU's but should (and can) be
  ! distributed
  ! in exactly the same fashion as delta_v_hartree_part_spl in
  ! update_hartree_potential_p1, if needed.
  if(.not. use_distributed_spline_storage) then
    do i_atom = 1, n_atoms, 1

      charge_l(i_atom) = real(floor(prev_hirshfeld_charge(i_atom)))
      charge_u(i_atom) = charge_l(i_atom) + 1.d0    
!      print *, species_name(species(i_atom)), prev_hirshfeld_charge(i_atom), charge_l(i_atom), charge_u(i_atom)
!      call get_ion_hirsh(i_atom, prev_hirshfeld_charge(i_atom))
      if (.not. rho_spl_ions_exists(species(i_atom),int(charge_l(i_atom))+4)) then
        call get_ion_hirsh(i_atom, charge_l(i_atom))
        rho_spl_ions(:,:,species(i_atom),int(charge_l(i_atom)+4)) = free_rho_spl(:,:,species(i_atom))
        rho_spl_ions_exists(species(i_atom),int(charge_l(i_atom))+4) = .true. 
      end if
      reference_rho_spl(:,:,i_atom) = rho_spl_ions(:,:,species(i_atom),int(charge_l(i_atom))+4)
!      call get_ion_hirsh(i_atom, prev_hirshfeld_charge(i_atom))
      if (.not. rho_spl_ions_exists(species(i_atom),int(charge_u(i_atom))+4)) then
        call get_ion_hirsh(i_atom, charge_u(i_atom))
        rho_spl_ions(:,:,species(i_atom),int(charge_u(i_atom))+4) = free_rho_spl(:,:,species(i_atom))
        rho_spl_ions_exists(species(i_atom),int(charge_u(i_atom))+4) = .true.
      end if
      reference_rho_spl1(:,:,i_atom) = rho_spl_ions(:,:,species(i_atom),int(charge_u(i_atom))+4)
    end do
  end if
  
  ! begin loop over batches of integration points
  i_full_points = 0

  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_full_points = i_full_points + 1

           coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

           call tab_atom_centered_coords_p0 &
           ( coord_current, &
             dist_tab_sq, &
             dir_tab, &
             n_centers_basis_integrals, centers_basis_integrals )
           
           point_on_atom = .false.
           do i_center = 1, n_centers_basis_integrals, 1
              if ( dist_tab_sq(i_center).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
                 point_on_atom = .true.
                 exit ! exit do loop
              end if
           end do

           ! avoid any evaluations on grid points that lie on another atom's nucleus
           if (.not.point_on_atom) then

              current_atom = batches(i_my_batch) % points(i_index) % index_atom
              current_radial = batches(i_my_batch) % points(i_index) % index_radial
              current_angular = batches(i_my_batch) % points(i_index) % index_angular

!test
!                write(use_unit,'(4X,6I10)') i_batch, i_index, i_full_points, current_atom, & 
!                                     current_radial, current_angular
!                write(use_unit,'(4X,F15.8)') r_radial(current_radial, species(current_atom))
!                write(use_unit,'(4X,3F15.8)') &
!                (r_angular(i_coord, current_angular, current_radial, species(current_atom)), i_coord=1,3,1)
! test end

!             tabulate current integration point as it appears from spherical
!             coordinates centered at each atom
 
              call tab_global_geometry_p0 &
             ( dist_tab_sq, &
               dir_tab, &
               dist_tab, &
               i_r, &
               dir_tab_norm, &
               n_centers_basis_integrals, centers_basis_integrals )

             ! determine the appropriate Hirshfeld partition weight explicitly.

             partition_norm = 0.d0
             free_partition_norm = 0.d0

             if(.not. use_distributed_spline_storage) then
                do i_center_L = 1, n_centers_basis_integrals, 1
                  
                   i_center = centers_basis_integrals(i_center_L)

                   aux_dens = val_spline & 
                        ( i_r(i_center_L), reference_rho_spl(1,1,center_to_atom(i_center)), & 
                        n_grid(species_center(i_center)) )
                   aux_dens1 = val_spline &
                        ( i_r(i_center_L),reference_rho_spl1(1,1,center_to_atom(i_center)), &
                        n_grid(species_center(i_center)) )

                   ! Fukui-function based finite difference approx. for ionic
                   ! density
                   aux_dens = (charge_u(center_to_atom(i_center))- &
                              prev_hirshfeld_charge(center_to_atom(i_center)))*aux_dens + & 
                              (prev_hirshfeld_charge(center_to_atom(i_center))-&
                              charge_l(center_to_atom(i_center)))*aux_dens1
                
 
                   partition_norm = partition_norm + aux_dens

                enddo


                ! present atom
                ! recompute i_r because I am not sure that the periodic lists retain the order 
                ! that the atoms in the initial unit cell are the first n_atoms atoms in each list.

                current_i_r = invert_log_grid & 
                     ( r_radial(current_radial, species(current_atom)), &
                     r_grid_min(species(current_atom)), r_grid_inc(species(current_atom)) )
                
                aux_dens = val_spline & 
                     ( current_i_r, reference_rho_spl(1,1,current_atom), & 
                     n_grid(species(current_atom)) )
                aux_dens1 = val_spline &
                     ( current_i_r, reference_rho_spl1(1,1,current_atom), &
                     n_grid(species(current_atom)) )

                ! Fukui-function based finite difference approx. for ionic
                ! density
                aux_dens = (charge_u(current_atom)-prev_hirshfeld_charge(current_atom))*aux_dens + &
                           (prev_hirshfeld_charge(current_atom)-charge_l(current_atom))*aux_dens1



             else
                ! Comment, VB: 
                ! This is a special case where the general reference_rho_spl
                ! (a per-atom quantity, introduced to later allow more sophisticated
                ! Hirshfeld-type analyses, e.g., Hirshfeld_I ,etc - where one needs per-atom
                ! partitioning densities) cannot be simply used, since teh spline array would
                ! get too large for large structures (100-1000 atoms) on the BlueGene (512 MB per thread, max.). 
                ! Thus, we only use the splined free-atom density, which is a per-species quantity, not per atom.
                !
                ! If anyone wants to make this work with reference_rho_spl per atom, that array would have to 
                ! be distributed across CPU's, and accessed in a very similar way as the splines in
                ! sum_up_whole_potential () - see that subroutine for details. Not hard to do, we know exactly how
                ! to do it, but someone would have to take up the rewrite & testing effort ...

                do i_center_L = 1, n_centers_basis_integrals, 1
                  
                   i_center = centers_basis_integrals(i_center_L)

                   aux_dens = val_spline & 
                        ( i_r(i_center_L),  free_rho_spl(1,1,species_center(i_center)), & 
                        n_grid(species_center(i_center)) )

                   partition_norm = partition_norm + aux_dens

                enddo

                ! present atom
                ! recompute i_r because I am not sure that the periodic lists retain the order 
                ! that the atoms in the initial unit cell are the first n_atoms atoms in each list.
                
                current_i_r = invert_log_grid & 
                     ( r_radial(current_radial, species(current_atom)), &
                     r_grid_min(species(current_atom)), r_grid_inc(species(current_atom)) )
                
                aux_dens = val_spline & 
                     ( current_i_r,  free_rho_spl(1,1,species(current_atom)), &
                     n_grid(species(current_atom)) )

             end if


             ! "naked" Hirshfeld weight at this point
             hirshfeld_partition = aux_dens / partition_norm
  
             ! ... and "dressed" with additional integration weights
             hirshfeld_partition = & 
             hirshfeld_partition * w_radial(current_radial, species(current_atom)) & 
             * (r_radial(current_radial, species(current_atom))**2.d0) & 
             * w_angular(current_angular, current_radial, species(current_atom)) & 
             * 4.d0 * pi

             ! Hirshfeld partition weight determined. Can now finally move on to actually use it. 

             ! temp_rho_new will be the difference between the actual charge density and the
             ! Hirshfeld superposition of atomic fragment reference densities
             temp_rho_new = 0.d0
             do i_spin = 1, n_spin, 1
               temp_rho_new = temp_rho_new + rho(i_spin, i_full_points)
             enddo
             full_dens = temp_rho_new
  
             temp_rho_new = temp_rho_new - pi4_inv * partition_norm

             ! Integrate Hirshfeld partial charge
             hirshfeld_charge_iter(current_atom) = hirshfeld_charge_iter(current_atom) &
               - hirshfeld_partition * temp_rho_new

           
             ! AT -- Hirshfeld volume partitioning
             free_int(current_atom) = free_int(current_atom) + &
                       r_radial(current_radial, species(current_atom))**3.d0 * & 
                       aux_dens &
                       *  w_radial(current_radial, species(current_atom)) &
                       * (r_radial(current_radial, species(current_atom))**2.d0) &
                       * w_angular(current_angular, current_radial,species(current_atom)) 

             hirsh_int(current_atom) = hirsh_int(current_atom) + &
                      r_radial(current_radial,species(current_atom))**3.d0 * &
                      hirshfeld_partition * full_dens 


             ! fragment dipole moment components
             do i_coord = 1,3,1
               hirshfeld_dipole(i_coord, current_atom) = hirshfeld_dipole(i_coord, current_atom) &
                 - hirshfeld_partition * temp_rho_new * r_radial(current_radial, species(current_atom))   &
                   * r_angular(i_coord, current_angular, current_radial, species(current_atom))
                 
             enddo
!test
!               if (current_atom.eq.2) then
!                 write(use_unit,'(3(2X,F15.8))') hirshfeld_dipole(:,current_atom)
!               end if
!test end


             ! fragment charge second moments
             i_sec_moment = 0
             do i_coord = 1,3,1
               do i_coord_2 = 1, i_coord, 1
                 i_sec_moment = i_sec_moment + 1
                 hirshfeld_quadrupole(i_sec_moment, current_atom) = hirshfeld_quadrupole(i_sec_moment, current_atom) &
                   - hirshfeld_partition * temp_rho_new * (r_radial(current_radial, species(current_atom))**2.d0)    &
                     * r_angular(i_coord, current_angular, current_radial, species(current_atom))                    &
                     * r_angular(i_coord_2, current_angular, current_radial, species(current_atom))
               enddo      
             enddo

! for spinpolarised calculations: compute spin-moment per atom	
	     if (spin_treatment .eq. 1) then
		deformation_density = 0.d0
		do i_spin = 1, n_spin, 1
			deformation_density(i_spin) = rho(i_spin, i_full_points) - free_rho_spl(1,1,species(current_atom))/2.0
		enddo	
		hirshfeld_spin_moment(current_atom) = hirshfeld_spin_moment(current_atom) + hirshfeld_partition &
			* ( deformation_density(1) - deformation_density(2) )
	     endif
	
           ! end if point_on_atom
           end if 

        !    end loop over a batch
        end do
  !     end loop over batches
  end do

  ! Synchronize integrals from different threads.
  call sync_hirshfeld ( hirshfeld_charge_iter, hirshfeld_spin_moment, hirshfeld_dipole, hirshfeld_quadrupole, &
                        free_int, hirsh_int )


  !if (use_vdw_correction_hirshfeld) then
  !  do i_atom = 1, n_atoms, 1
  !    hirshfeld_volume(i_atom) = hirsh_int(i_atom) / free_int(i_atom)
  !  enddo
  !end if

  ! Check whether Hirshfeld charges are converged in the iterative procedure

    converged = .true.
    do i_atom= 1, n_atoms, 1
      if (abs(hirshfeld_charge_iter(i_atom)) .gt. 0.0001) then
        converged = .false. 
      end if
    end do

   prev_hirshfeld_charge(:) = prev_hirshfeld_charge(:) + hirshfeld_charge_iter(:)

  ! Write output in appropriate form

  write(info_str,'(2X,A)') "----------------------------------------------------------------------"
  call localorb_info( info_str, use_unit,'(A)',OL_norm)

  write(info_str,'(2X,A,I5)') "| Hirshfeld-I iterations ", i_iter
  call localorb_info( info_str, use_unit,'(A)',OL_norm)

  do i_atom = 1, n_atoms, 1

    write(info_str,'(2X,A,I5,A,A2)') "| Atom ", i_atom, ": ", &
      species_name(species(i_atom))
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8)') "|   Hirshfeld charge        : ", &
      prev_hirshfeld_charge(i_atom)
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8)') "|   Reference ionic volume  : ", &
      free_int(i_atom)
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,F15.8)') "|   Hirshfeld volume        : ", &
      hirsh_int(i_atom)
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    if ( spin_treatment .eq. 1) then
        write(info_str,'(2X,A,F15.8)') "|   Hirshfeld spin moment   : ", &
                hirshfeld_spin_moment(i_atom)
        call localorb_info( info_str, use_unit,'(A)',OL_norm)
    endif

    write(info_str,'(2X,A,3(F15.8,2X))') "|   Hirshfeld dipole vector : ", &
      ( hirshfeld_dipole(i_coord, i_atom) * bohr , i_coord = 1,3,1 )
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    dipole_moment = 0.d0
    do i_coord = 1,3,1
      dipole_moment = dipole_moment + hirshfeld_dipole(i_coord, i_atom)**2.d0
    enddo
    dipole_moment = sqrt(dipole_moment) * bohr
    write(info_str,'(2X,A,F15.8)') "|   Hirshfeld dipole moment : ", &
      dipole_moment
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A,3(F15.8,2X))') "|   Hirshfeld second moments: ", &
      hirshfeld_quadrupole(1, i_atom) * bohr**2.d0, hirshfeld_quadrupole(2,i_atom) * bohr**2.d0, &
      hirshfeld_quadrupole(4, i_atom) * bohr**2.d0
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,24X,3(F15.8,2X))') "|     ", &
      hirshfeld_quadrupole(2, i_atom) * bohr**2.d0, hirshfeld_quadrupole(3,i_atom) * bohr**2.d0, &
      hirshfeld_quadrupole(5, i_atom) * bohr**2.d0
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,24X,3(F15.8,2X))') "|     ", &
      hirshfeld_quadrupole(4, i_atom) * bohr**2.d0, hirshfeld_quadrupole(5,i_atom) * bohr**2.d0, &
      hirshfeld_quadrupole(6, i_atom) * bohr**2.d0
    call localorb_info( info_str, use_unit,'(A)',OL_norm)

    write(info_str,'(2X,A)') "----------------------------------------------------------------------"
    call localorb_info( info_str, use_unit,'(A)',OL_norm)
  enddo

  call localorb_info ('')

  enddo ! iterations

  ! finally deallocate everything, well, whatever was allocated initially.
  if (allocated(hirshfeld_charge_iter))      deallocate ( hirshfeld_charge_iter )
  if (allocated(prev_hirshfeld_charge)) deallocate ( prev_hirshfeld_charge )
  if (allocated(free_int))              deallocate ( free_int )
  if (allocated(hirsh_int))             deallocate ( hirsh_int )
  if (allocated(hirshfeld_spin_moment)) deallocate ( hirshfeld_spin_moment )
  if (allocated(hirshfeld_dipole))      deallocate ( hirshfeld_dipole )
  if (allocated(hirshfeld_quadrupole))  deallocate ( hirshfeld_quadrupole )
  if (allocated(free_int))              deallocate ( free_int )
  if (allocated(hirsh_int))             deallocate ( hirsh_int )
  if (allocated(reference_rho_spl))     deallocate ( reference_rho_spl )
  if (allocated(reference_rho_spl1))     deallocate ( reference_rho_spl1 )

  ! If output was requested under all circumstances, must reduce the
  ! output priority for this subroutine. Must reset at end of subroutine!
  if (out_hirshfeld_always) then
     if (output_level.eq.'MD_light') output_priority = output_priority_old
  end if

end subroutine hirshfeld_analysis_iterative
!******	
