!****s* FHI-aims/get_lda_free_atom_densities
!  NAME
!   get_lda_free_atom_densities
!  SYNOPSIS

     subroutine get_lda_free_atom_densities &
     ( density_spl, density_deriv_spl, fill_derivative )

!  PURPOSE
!  Subroutine get_lda_free_atom_densities provides the free-atom density in the
!  local density approximation, non spin-polarized. This is only used for partition tables, which enter the
!  integration weights in all integrals and in the Hartree potential. The purpose is
!  to use one and the same set of integration weights in all integrations, independent
!  of the XC functional that is used in the actual self-consistent electronic structure
!  calculation later. 
!
!  this used to be a wrapper around Martin Fuchs' sr_atom.f
!  It is now a procedure that calls atomic solvers, but it still bears the 
!  marks of its original purpose.
!
!  We do switch according to the specified relativistic treatment ...
!
!  USES

      use dimensions
      use runtime_choices
      use grids
      use species_data
      use spline
      use localorb_io
      use constants
      use synchronize_mpi, only: sync_atomic_solver_outputs
      use psi_at_nucleus_mod, only: psi_at_nucleus_atomic
      use mpi_tasks, only: aims_stop, check_allocation
      implicit none

!  ARGUMENTS
!  INPUTS
     logical :: fill_derivative
!    fill_derivative - logical variable that specifies whether we touch the
!                      density_deriv_spl array, or not
!  OUTPUT
     real*8 :: density_spl(n_max_spline,n_max_grid, n_species)
     real*8 :: density_deriv_spl(n_max_spline,n_max_grid, n_species)
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
!    Release version, FHI-aims (2009).
!  SOURCE
!

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
!              max_shells INDEX FOR ANYTHING, IT IS A HARTREE_FO!K RELI!!
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

      real*8 nuclear_charge
      real*8 ::  dr_coef, alpha

      integer n_core
      integer n_valence

      integer, dimension(:), allocatable :: n_shell
      integer, dimension(:), allocatable :: l_shell
      integer, dimension(:), allocatable :: m_shell
      real*8, dimension(:), allocatable :: occ_shell
      real*8, dimension(:), allocatable :: eigenval
      integer, dimension(:), allocatable :: i_grid_max
      logical, dimension(:), allocatable :: core_type

      real*8 d_eig_sum

      real*8, dimension(:), allocatable :: r_peak

      real*8, dimension(:), allocatable :: ini_pot
      real*8, dimension(:,:), allocatable :: atom_pot
      real*8, dimension(:), allocatable :: atom_pot_es
      real*8, dimension(:,:), allocatable :: atom_pot_xc
      real*8, dimension(:,:), allocatable :: density
      real*8, dimension(:,:), allocatable :: density_deriv
      real*8, dimension(:,:), allocatable :: density_2nd_deriv
      real*8, dimension(:,:), allocatable :: core_density
      real*8, dimension(:,:), allocatable :: core_dd
      real*8, dimension(:,:), allocatable :: core_d2d

      real*8, dimension(:,:), allocatable :: wave
      real*8, dimension(:,:), allocatable :: wave_deriv

      real*8 :: outermost_charge

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
      character*120 :: info_str

      real*8 :: multipole_radius_old

      real*8, dimension(:,:), allocatable :: kinetic

!  counters

      integer i_species
      integer i_grid
      integer i_channels
      integer i_function, info

!  begin work

      write (info_str,'(2X,A)') &
        "| Obtaining LDA free-atom density for the partition table(s)."
      call localorb_info(info_str,use_unit,'(A)')      

!  allocations - needed because some of these arrays may collide with the
!  stack size on small machines if they were automatic arrays

      allocate(n_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'n_shell                       ') 

      allocate(l_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'l_shell                       ') 

      allocate(m_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'm_shell                       ') 

      allocate(occ_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'occ_shell                     ') 

      allocate(eigenval(n_max_ind_fns),stat=info)
      call check_allocation(info, 'eigenval                      ') 

      allocate(i_grid_max(n_max_ind_fns),stat=info)
      call check_allocation(info, 'i_grid_max                    ') 

      allocate(r_peak(n_max_ind_fns),stat=info)
      call check_allocation(info, 'r_peak                        ') 

      allocate(core_type(n_max_ind_fns),stat=info)
      call check_allocation(info, 'core_type                     ') 

      allocate(ini_pot(n_max_grid),stat=info)
      call check_allocation(info, 'ini_pot                       ') 

      allocate(atom_pot(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'atom_pot                      ') 

      allocate(atom_pot_es(n_max_grid),stat=info)
      call check_allocation(info, 'atom_pot_es                   ') 

      allocate(atom_pot_xc(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'atom_pot_xc                   ') 

      allocate(density(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'density                       ') 

      allocate(density_deriv(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'density_deriv                 ') 

      allocate(density_2nd_deriv(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'density_2nd_deriv             ') 

      allocate(core_density(n_max_grid,n_channels),stat=info)
      call check_allocation(info, 'core_density                  ') 

      allocate(core_dd(n_max_grid,n_channels),stat=info)
      call check_allocation(info, 'core_dd                       ') 

      allocate(core_d2d(n_max_grid,n_channels),stat=info)
      call check_allocation(info, 'core_d2d                      ') 
                                                                                       
      allocate(wave(n_max_grid, n_max_ind_fns),stat=info)
      call check_allocation(info, 'wave                          ') 

      allocate(e_sic(n_max_ind_fns),stat=info)
      call check_allocation(info, 'e_sic                         ') 

      allocate(v_orb(n_max_grid,n_max_ind_fns),stat=info)
      call check_allocation(info, 'v_orb                         ') 

      allocate(sv_kli(n_max_ind_fns),stat=info)
      call check_allocation(info, 'sv_kli                        ') 

      allocate(svm(n_max_ind_fns),stat=info)
      call check_allocation(info, 'svm                           ') 

      allocate(wave_deriv(n_max_grid, n_max_ind_fns),stat=info)
      call check_allocation(info, 'wave_deriv                    ') 

      allocate(kinetic(n_max_grid, n_max_ind_fns), stat=info)
      call check_allocation(info, 'kinetic                       ') 

!  input parameters for sr_atom()

      ! We always use the Perdew-Wang 1991 LDA for the free-atom densities produced here.
      i_exch = 8

      if (fill_derivative) then
        i_grad = 1
      else
        i_grad = 0
      end if

      ! output density varies with relativistic treatment for now, this is 
      ! not necessarily ideal ... should simply always use Koelling-Harmon here.
      !
      ! This will be the next step, after debugging.
      if(flag_KH_core_states)then
         i_rel = 10
      elseif (flag_rel.eq.0) then
        i_rel=2
      else if ((flag_rel.eq.REL_zora).or.(flag_rel.eq.REL_atomic_zora)) &
            then
        i_rel=9  
      else if (flag_rel.eq.REL_own.or.flag_rel==REL_KOLNING_HARMON)then
        i_rel = 10
      end if
      
      !  more input parameters for sr_atom()

      sc_max_iter = 100
      t_core = .true.
      n_core = 0
      t_in_pot = .false.
      t_sic = .false.
      t_kli = .false.
      t_mgga = .false.

      do i_species = 1,n_species,1
!        if(species_pseudoized(i_species)) cycle

        !  initialize also the unused core density parts
        do i_grid = 1, n_max_grid, 1
          do i_channels = 1,n_channels, 1
            core_density(i_grid,1) = 0.d0
            core_dd(i_grid,1) = 0.d0
            core_d2d(i_grid,1) = 0.d0
          enddo
        enddo

        nuclear_charge = species_z(i_species)

        !  filled states only

!       add total number of (n,l) shells
!       store l(shell), n(shell) in separate arrays
!       FIXME - When we switch array places i_species <-> i_function,
!       this ugly crudge can go entirely.
        n_valence = n_atomic(i_species)
        do i_function = 1, n_atomic(i_species), 1

            n_shell (i_function) = atomic_n(i_species, i_function)
            l_shell (i_function) = atomic_l(i_species,i_function)
            m_shell (i_function) = 1
            occ_shell (i_function) = &
                 valence_occ &
                 (n_shell(i_function),l_shell(i_function),i_species)

            core_type(i_function) = core_fn(i_species,i_function)
        enddo
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

        do i_grid = 1, n_grid(i_species), 1
          atom_pot(i_grid, 1) = ini_pot(i_grid)
          atom_pot_xc(i_grid, 1) = 0.d0
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
            n_max_grid, n_max_ind_fns, n_channels, &
            n_grid(i_species), r_grid(1,i_species),  &
            cut_free_atom(i_species), cutoff_type(i_species), &
            free_r_cut(i_species),&
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

          ! Note: Because we want PW91-LDA quantities here, we hard-core flag_xc = 8 here
          call atom_sphere_wrapper &
          ( n_grid(i_species), r_grid(1,i_species), species_name(i_species), 8, hybrid_coeff, hse_omega, &
            nuclear_charge, n_core, n_valence, n_max_ind_fns, n_shell, l_shell, m_shell, occ_shell, n_max_grid, & 
            n_channels, 1, &
            cut_free_atom(i_species), cutoff_type(i_species), free_r_cut(i_species), scale_cutoff(i_species), &
                 w_cutoff(i_species), & ! Cutoff potential
            eigenval, wave, wave_deriv, density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
            kinetic, psi_at_nucleus_atomic(:n_valence,i_species)) ! Output
          call sync_atomic_solver_outputs( n_max_ind_fns, n_max_grid, n_channels, 1, &
               eigenval, wave, wave_deriv, kinetic, &
               density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es )
        else
          call aims_stop('The atomic solver specified by the variable atomic_solver was not recognized, exiting.','func')  
        end if

!       now copy the density onto storage arrays

!       density has an extra index which used to be shell-dependent (Hartree-Fock?),
!       then indicated spin up / down (Martin Fuchs), and is unused here. 

!       spline rho
        call cubic_spline &
        ( density(1,1), n_grid(i_species),  &
          density_spl(1,1,i_species) )

        if (fill_derivative) then
!         spline drho_dr
          call cubic_spline &
          ( density_deriv(1,1), n_grid(i_species), &
            density_deriv_spl(1,1,i_species) )
        end if

!  Finally: Find the radius outside of which all partition tables
!  (via the density) are effectively zero for this atom.
!
!  Notice that we have _already_ done this once, in get_free_atoms.f90, but
!  there a different functional was used. Are we sure that we now integrate everything correctly?

        ! Store multipole radius of free atoms, using the XC specified in control.in
        multipole_radius_old = multipole_radius_free(i_species)

        i_grid = n_grid(i_species)
        do while ( (dabs(density(i_grid, 1)).le. &
             multipole_radius_free_threshold) .and.  &
                    (i_grid.ge.1) )
          i_grid = i_grid-1
        enddo
 
        ! Set the outermost radius to the innermost point where the 
        ! free-atom charge density is still zero.
        if (i_grid.eq.n_grid(i_species)) then
          outermost_charge = r_grid( i_grid, i_species )
        else
          outermost_charge = r_grid( i_grid+1, i_species )
        end if

        multipole_radius_free(i_species) = max(outermost_charge, multipole_radius_old)
        multipole_radius_free_sq(i_species) = multipole_radius_free(i_species)**2

        ! write test output to make sure we understand the impact of this
        ! somewhat ad-hoc free-atom radius modification ...
        write (info_str,'(2X,3A)') &
        "| Species ", trim(species_name(i_species)), ":"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(2X,A,F15.8,A)') &
        "|   Outermost radius of density used for partition table      : ", outermost_charge*bohr, " AA."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(2X,A,F15.8,A)') &
        "|   Outermost radius of free-atom density for XC in control.in: ", outermost_charge*bohr, " AA."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(2X,A)') &
        "|   We use the max. of either value to determine the extent of all summations / integrations later."
        call localorb_info(info_str,use_unit,'(A)')

        ! end loop over species        
      enddo

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
      deallocate(atom_pot_es)
      deallocate(ini_pot)
      deallocate(core_type)
      deallocate(r_peak)
      deallocate(i_grid_max)
      deallocate(eigenval)
      deallocate(occ_shell)
      deallocate(m_shell)
      deallocate(l_shell)
      deallocate(n_shell)

      return
    end subroutine get_lda_free_atom_densities
!******
