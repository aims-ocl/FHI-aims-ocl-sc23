!****s* FHI-aims/get_free_atoms_polarized
!  NAME
!   get_free_atoms_polarized
!  SYNOPSIS

subroutine get_free_atoms_polarized(  )

! PURPOSE
!  Subroutine get_free_atom provides atomic input data for each species
!
!  this used to be a wrapper around Martin Fuchs' sr_atom.f
!  It is now a procedure that calls atomic solvers, but it still bears the 
!  marks of its original purpose.
!
!  We use a non-relativistic version for now. For a proper relativistic treatment, 
!  must have a radial relativistic Hamiltonian which can be extended directly to 
!  the full (non-radial-symmetric) potential. This is not the case at present.
! spin-polarized version (R.Gehrke)
! USES
  use constants,       only : bohr, hartree, pi4_inv
  use dimensions,      only : n_wave_max, l_wave_max, n_spin, n_channels, &
                              n_ini_type, n_max_grid, n_max_ind_fns, &
                              use_density_gradient
  use runtime_choices, only : atomic_solver, atomic_solver_atom_sphere, &
                              atomic_solver_sratom, flag_KH_core_states, &
                              flag_rel, flag_xc, hse_omega, hybrid_coeff, & 
                              REL_atomic_zora, REL_kolning_harmon, REL_own
  use grids,           only : n_grid, r_grid
  use species_data,    only : species_z, species_name, n_atomic, atomic_n, &
                              atomic_l, core_fn, n_atomic, r_cutoff, &
                              cutoff_type, cutoff_type, scale_cutoff, w_cutoff
  use free_atoms,      only : initial_pot_es, initial_rho, initial_drho_dr, &
                              initial_rho_spl, initial_pot_es_spl, &
                              initial_drho_dr_spl, get_atomic_occ_numbers
  use geometry,        only : type_species, type_charge, type_kind, type_moment
  use spline,          only : cubic_spline
  use xc_library,      only : xc__pw_lda, xc_sr__pw_lda, xc__pbe, xc_sr__pbe, &
                              xc__pbe0, is_dfauto, normalize_flag_xc
  use localorb_io,     only : use_unit, localorb_info
  use mpi_tasks,       only : aims_stop
  use synchronize_mpi, only : sync_atomic_solver_outputs
  use psi_at_nucleus_mod, only: psi_at_nucleus_atomic
  implicit none
! INPUT
!  none
!  OUTPUT
!  none
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
      

!  imported variables

!  local variables

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
!              FIXME: SAME AS atom_pot - THE max_shells INDEX IS COMPLETELY USELESS!
!  core_density, core_dd, core_d2d : UNUSED - core density, density derivative, and second 
!              derivative, for frozen core calculations (never done here)
!  wave :      s.-c. atomic wave functions     
!  t_in_pot : if true, sratom() assumes an input non-selfconsistent potential - always false for now.
!  t_sic :    if true, sratom() does self-interaction correction - this is never the case. 
!  e_sic :    array of self-energies for SIC - not used
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
  
  real*8, dimension(:,:), allocatable :: ini_pot
  real*8, dimension(:,:), allocatable :: atom_pot
  real*8, dimension(:,:), allocatable :: atom_pot_xc
  real*8, dimension(:,:), allocatable :: atom_pot_es
  real*8, dimension(:,:), allocatable :: density
  real*8, dimension(:,:), allocatable :: density_deriv
  real*8, dimension(:,:), allocatable :: density_2nd_deriv
  real*8, dimension(:,:), allocatable :: core_density
  real*8, dimension(:,:), allocatable :: core_dd
  real*8, dimension(:,:), allocatable :: core_d2d
  
  real*8, dimension(:,:), allocatable :: wave
  real*8, dimension(:,:), allocatable :: wave_deriv
  real*8, dimension(:,:), allocatable :: kinetic
  
  logical t_in_pot 
  logical t_sic
  real*8, dimension(:), allocatable :: e_sic
  real*8, dimension(:,:), allocatable :: v_orb
  
  logical t_kli
  real*8, dimension(:), allocatable :: sv_kli
  real*8, dimension(:), allocatable :: svm
      
  logical t_mgga

!  other aux variables

      ! dimension for spin-polarized case
  integer n_max_ind_fns_pol
  integer i_function_start

  character*100 :: info_str

  real*8, dimension(1:n_wave_max, 0:l_wave_max, n_spin) :: atomic_occ_numbers

  integer :: comp_species
  real*8  :: dummy_moment
  integer :: dummy_kind
  logical :: cut_free_atom_to_pass
  real*8  :: free_r_cut_to_pass

!  counters

  integer i_grid
  integer i_channel
  integer i_function
  integer i_type

  ! Rundong: for relativistic case, I currently don't use this subroutine
  if(flag_rel.eq.7 .or. flag_rel.eq.8)return

!  begin work

  n_max_ind_fns_pol = n_spin * n_max_ind_fns

  write (info_str,'(2X,A,A)') &
       "Creating spin-polarized density ", &
       "from free atoms."
  
  call localorb_info('',use_unit)
  call localorb_info(info_str,use_unit,'(A)')
  
  !  allocations - needed because some of these arrays may collide with the
  !  stack size on small machines if they were automatic arrays
  
  allocate(n_shell(n_max_ind_fns_pol))
  allocate(l_shell(n_max_ind_fns_pol))
  allocate(m_shell(n_max_ind_fns_pol))
  allocate(occ_shell(n_max_ind_fns_pol))
  allocate(eigenval(n_max_ind_fns_pol))
  allocate(i_grid_max(n_max_ind_fns_pol))
  allocate(r_peak(n_max_ind_fns_pol))
  allocate(core_type(n_max_ind_fns_pol))
  
  allocate(ini_pot(n_max_grid, n_spin))
  allocate(atom_pot(n_max_grid, n_spin))
  allocate(atom_pot_xc(n_max_grid, n_spin))
  allocate(atom_pot_es(n_max_grid, n_spin))
  allocate(density(n_max_grid, n_spin))
  allocate(density_deriv(n_max_grid, n_spin))
  allocate(density_2nd_deriv(n_max_grid, n_spin))
  allocate(core_density(n_max_grid,n_spin))
  allocate(core_dd(n_max_grid,n_spin))
  allocate(core_d2d(n_max_grid,n_spin))
  
  allocate(wave(n_max_grid, n_max_ind_fns_pol))
 
  allocate(e_sic(n_max_ind_fns_pol))
  allocate(v_orb(n_max_grid,n_max_ind_fns_pol))
  allocate(sv_kli(n_max_ind_fns_pol))
  allocate(svm(n_max_ind_fns_pol))
  allocate(wave_deriv(n_max_grid, n_max_ind_fns_pol))
      
  allocate(kinetic(n_max_grid, n_max_ind_fns_pol))

  !  input parameters for sr_atom(), from control.in 

  i_exch = flag_xc
  !  for Hartree-Fock calculation --Ren

  if(flag_xc .eq. 0) then
     i_exch = 8
  else if(flag_xc .eq. 1) then
    !   Hybrid-PBE0 calculation, using PBE instead as the atomic solver.
     i_exch = 6
  end if
! HSE calculation or lc-wpbeh, using PBE
  if (flag_xc .eq. 7 .or. flag_xc.eq.23) then
    i_exch = 6 
  endif
  if(flag_xc .eq. 10) then
!  Hybrid-B3LYP calculation, using BLYP instead as the atomic solver.
        i_exch = 9
  endif
! VWN-LDA approximation, using VWN as initial guess.
! Unfortunately it is number 7 in vexcor.f
! I hope there will be no mixing with HSE
      if (flag_xc.eq.15) then
        i_exch = 7
      endif
! revPBE 
  if(flag_xc .eq. 12) then
        i_exch = 12
  endif 
!revPBE_vdw.   SAG
 if(flag_xc.eq.14) then
        i_exch = 12
  endif 

! RPBE
  if(flag_xc .eq. 11) then
        i_exch = 11
  endif 
! PBEsol or PBEsol0 calculation, using PBE as the atomic solver
  if(flag_xc .eq. 17 .or. flag_xc.eq.13) then
     i_exch = 6
  end if
  !R48PBE
  if(flag_xc.eq.18) then
     i_exch = 6
  end if
!PBE_vdw. SAG
  if(flag_xc.eq.5) then
     i_exch = 6
  end if

! FIXME : A case should be added in external/vexcor.f in order to
! calculate the AM05 energies and potentials (using subroutine
! external/am05.f) and use that as atomic solver!
! AM05 calculation, using PW-LDA instead as the atomic solver
  if (flag_xc.eq.20) then
     i_exch = 8
  end if

! Caveat to catch all Meta-GGAs for now.
! We'll solve calculate the free atom for them all as PW91-LD
! similar to the approach above used for AM05. AJL
  if (flag_xc.ge.25) then
     i_exch = 8
  endif

!xPBE  IYZ
  if (flag_xc.eq.22) then
     i_exch = 8
  endif

! Caveat to catch LibXC 
! We'll solve calculate the free atom for them all as PW91-LD
! similar to the approach above used for AM05. AJL
  if (flag_xc.le.-10 .or. flag_xc.eq.22) then
     i_exch = 8
  endif

  if (is_dfauto(flag_xc)) then
      select case (normalize_flag_xc(flag_xc))
      case (xc__pw_lda)
          i_exch = xc_sr__pw_lda
      case (xc__pbe, xc__pbe0)
          i_exch = xc_sr__pbe
      case default
          i_exch = xc_sr__pw_lda
      end select
  end if

! Decide whether gradients are needed or not, depending on xc scheme
  if ((i_exch.eq.3).or.(i_exch.eq.8).or.(i_exch.eq.7)) then
     i_grad = 0
  else if ((i_exch.eq.4).or.(i_exch.eq.6).or.(i_exch.eq.9) &
            .or.(i_exch.eq.11).or.(i_exch.eq.12)) then
     i_grad = 1
  end if
  
  if (use_density_gradient) then
     i_grad = 1
  end if
  
  if(flag_KH_core_states)then
     i_rel = 10
  elseif (flag_rel.eq.0) then
     i_rel = 2
  else if ((flag_rel.eq.1).or.(flag_rel.eq.REL_atomic_zora)) then
    !TJ Changed to ZORA
     i_rel = 9  
  else if (flag_rel.eq.REL_own.or.flag_rel==REL_KOLNING_HARMON) then
     i_rel = 10
  end if
  
  if (n_spin .eq. 2) then
     ! spin-polarized case
     ! remap to spin-polarized flags which are different in the atomic solver - VERY MESSY!!!

     if (i_exch .eq. 6) then
        ! PBE
        i_exch = -9
!@@edu>
     else if (i_exch .eq. 21) then 
        ! PBEint
        i_exch = -21
!@@edu<
     else if (i_exch .eq. 4) then
        ! PW91_gga
        i_exch = -6
     else if (i_exch .eq. 8) then
        ! PW-LDA
        i_exch = -3
     else if (i_exch .eq. 9) then
        ! BLYP
        i_exch = -85
     else if (i_exch .eq. 7) then
        ! VWN5
        i_exch = -15
     else if (i_exch .eq. 16) then
        ! VWN5
        i_exch = -16
     else if (i_exch .eq. 12) then
        ! revPBE
        i_exch = -12
     else if (i_exch .eq. 11) then
        ! revPBE
        i_exch = -11
     else if (i_exch .eq. 3) then
        ! PZ-LDA
        ! not supported
        write(info_str,'(1X,A,A)') "*** WARNING: PZ-LDA currently not implemented ", &
             "in the atomic solver for spin-polarized systems!"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(1X,A)') "*** Initial density from Perdew-Wang 1991 LDA instead. "
        call localorb_info(info_str,use_unit,'(A)')
        i_exch = -3  ! VB: This is the Perdew-Wang LDA, not the Perdew-Zunger one!
     else 
        write(info_str,'(2X,A,A)') "WARNING: Serious inconsistency in ", &
             "get_free_atoms_polarized!"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)') "* Aborting. "
        call localorb_info(info_str,use_unit,'(A)')
        stop
     end if

  end if
  
  !  more input parameters for sr_atom()
  !  FIXME: Some of these should be adjustable in control.in !

  sc_max_iter = 100
  t_core = .true.
  n_core = 0
  t_in_pot = .false.
  t_sic = .false.
  t_kli = .false.
  t_mgga = .false.
  
  do i_type = 1, n_ini_type, 1

     if (n_spin .eq. 2) then
        ! spin-polarized case
        ! since the brilliant atomic solver changes i_exch to abs(i_exch)
        ! (so it completely ignores that i_exch is an INPUT-parameter,
        !  that's why fortran90 with intent(in)-statements helps a lot...)
        ! I need to set it to the spin-polarized value again and again...
        ! R. Gehrke
        i_exch = - abs(i_exch)
     end if
     
     core_density = 0.d0
     core_dd      = 0.d0
     core_d2d     = 0.d0
     
     comp_species = type_species(i_type)
     
     if (n_spin .eq. 2) then
       
        ! spin polarized case
        call get_atomic_occ_numbers &
             (species_z(comp_species) - type_charge(i_type), &
             type_moment(i_type), type_kind(i_type), comp_species, &
             atomic_occ_numbers, .false.) 
        
     else
        
        dummy_moment = 0.d0
        dummy_kind = 1

        ! non-polarized case
        call get_atomic_occ_numbers &
             (species_z(comp_species) - type_charge(i_type), &
             dummy_moment, dummy_kind, comp_species, & 
             atomic_occ_numbers, .false.) 
        
     end if
     
     call localorb_info('',use_unit)
     write(info_str,'(2X,A,A)') "Species: ", species_name(comp_species)
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A,F10.6)') "Charge: ", type_charge(i_type)
     call localorb_info(info_str,use_unit,'(A)')
     if (n_spin .eq. 2) then
        if (type_kind(i_type) .eq. 1) then
           write(info_str,'(2X,A)') "applying Hund's rules ... "
        else
           write(info_str,'(2X,A,F10.6)') "Moment: ", type_moment(i_type)
        end if
        call localorb_info(info_str,use_unit,'(A)')
     end if
     
     nuclear_charge = species_z(comp_species)
     
     ! filled states only

     ! add total number of (n,l) shells
     ! store l(shell), n(shell) in separate arrays

     n_valence = n_atomic(comp_species)
     
     ! spin-up
     do i_function = 1, n_atomic(comp_species), 1
        
        n_shell (i_function) = atomic_n(comp_species, i_function)
        l_shell (i_function) = atomic_l(comp_species, i_function)
        m_shell (i_function) = 1
        occ_shell (i_function) = atomic_occ_numbers &
             (n_shell(i_function), l_shell(i_function), 1)
        core_type(i_function) = core_fn(comp_species,i_function)
     end do

     if (n_spin .eq. 2) then
        ! spin-down
        do i_function = n_atomic(comp_species) + 1, 2 * n_valence, 1
        
           n_shell (i_function) = n_shell (i_function - n_atomic(comp_species))
           l_shell (i_function) = l_shell (i_function - n_atomic(comp_species))
           m_shell (i_function) = 2
           occ_shell (i_function) = atomic_occ_numbers &
                (n_shell(i_function), l_shell(i_function), 2)
        core_type(i_function) = core_fn(comp_species,i_function- n_atomic(comp_species))
        end do
     end if

     ! in the case that there is no electronic charge at all, don't start
     ! the atomic solver

     if (nuclear_charge .ne. type_charge(i_type)) then

        do i_channel = 1, n_spin, 1

           i_function_start = (i_channel - 1) * n_atomic(comp_species) + 1
           call atomini &
                (n_valence, n_grid(comp_species), nuclear_charge, &
                n_shell(i_function_start), occ_shell(i_function_start), &
                r_grid(1,comp_species), ini_pot(1, i_channel), eigenval(i_function_start) &
                )
     
           ! FIXME: This would be a superfluous copy, if atom_pot were not grossly overdimensioned
           do i_grid = 1, n_grid(comp_species), 1
              atom_pot(i_grid, i_channel) = ini_pot(i_grid, i_channel)
              atom_pot_xc(i_grid, i_channel) = 0.d0
              atom_pot_es(i_grid, i_channel) = 0.d0
           end do
        end do

        ! Set cutoff potential for initial, spin-polarized free atom density. 
        ! In principle, a cutoff potential should always apply here.
        ! There is no reason for the atomic density in the initialization
        ! to exceed the cutoff radius used for the basis functions in any way.
        cut_free_atom_to_pass = .true.
        free_r_cut_to_pass = r_cutoff(comp_species)
        if ((species_z(comp_species) .eq. 1) .and. &
             ((n_spin .eq. 2) .and. ((type_moment(i_type) .eq. 1) .or. (type_kind(i_type) .eq. 1))) ) then
           ! well known pathological case where a cutoff-potential is necessary
           ! limit free-atom cutoff to 5.0 AA if not already done.
           free_r_cut_to_pass = min (5.d0/bohr, free_r_cut_to_pass )
        end if

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
               n_max_grid, n_max_ind_fns_pol, n_spin, &
               n_grid(comp_species), r_grid(1,comp_species), &
               cut_free_atom_to_pass, cutoff_type(comp_species), &
               free_r_cut_to_pass, &
               scale_cutoff(comp_species), w_cutoff(comp_species), &
               r_peak, atom_pot, atom_pot_es, &
               atom_pot_xc, density, density_deriv, &
               density_2nd_deriv, core_density, core_dd, core_d2d, &
               wave, t_in_pot, t_sic, e_sic, v_orb, t_kli, sv_kli, &
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
               ( n_grid(comp_species), r_grid(1,comp_species), species_name(comp_species), flag_xc, hybrid_coeff, hse_omega, &
                 nuclear_charge, n_core, n_valence, n_max_ind_fns_pol, n_shell, l_shell, m_shell, occ_shell, n_max_grid, & 
                 n_channels, n_spin, &
                 cut_free_atom_to_pass, cutoff_type(comp_species), free_r_cut_to_pass, scale_cutoff(comp_species), &
                          w_cutoff(comp_species), & ! Cutoff potential
                 eigenval, wave, wave_deriv, density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
                 kinetic, psi_at_nucleus_atomic(:n_atomic(comp_species),comp_species)) ! Output
          call sync_atomic_solver_outputs( n_max_ind_fns_pol, n_max_grid, n_channels, n_spin, &
               eigenval, wave, wave_deriv, kinetic, &
               density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
               psi_at_nucleus_atomic(:n_atomic(comp_species),comp_species))
        else
          call aims_stop('The atomic solver specified by the variable atomic_solver was not recognized, exiting.','func')  
        end if

!       now copy all data from atom onto storage arrays

!       atom_pot and density have an extra index which used to be shell-dependent (Hatree-Fock?),
!       then indicated spin up / down (Martin Fuchs), and is unused here. 

        do i_grid = 1, n_grid(comp_species), 1
          initial_pot_es(i_grid, i_type) = atom_pot_es(i_grid, 1)
        end do

        do i_channel = 1, n_spin, 1
           do i_grid = 1, n_grid(comp_species), 1
              initial_rho    (i_grid, i_type, i_channel) = pi4_inv * density      (i_grid, i_channel)
              ! This version breaks for the case of an lda-calculation if it's compiled with ifort 9.1.036 and -O3
              ! The strange thing is, if the line with initial_drho_dr is commented out it works, though this line
              ! is NEVER called in the case of an LDA-calculation !!!
              ! I guess it is an optimizer bug. (no segfault with -O0)
              ! With XL Fortran 10.1 (on the thp's) this version works as well.
              
              ! if (use_density_gradient) then
              !   initial_drho_dr(i_grid, i_type, i_channel) = pi4_inv * density_deriv(i_grid, i_channel)
              ! end if
              
           end do
           ! That way it works with ifort 9.1.036 and -O3 as well.
           ! R. Gehrke
           if (use_density_gradient) then
              do i_grid = 1, n_grid(comp_species), 1
                 initial_drho_dr(i_grid, i_type, i_channel) = pi4_inv * density_deriv(i_grid, i_channel)
              end do
           end if
           ! spline rho
           call cubic_spline &
                (initial_rho(1, i_type, i_channel), n_grid(comp_species), &
                initial_rho_spl(1, 1, i_type, i_channel) )
           if (use_density_gradient) then
              ! spline drho_dr
              call cubic_spline &
                   (initial_drho_dr(1, i_type, i_channel), n_grid(comp_species), &
                   initial_drho_dr_spl(1, 1, i_type, i_channel) )
           end if
        end do
        
        ! spline potentials
        call cubic_spline &
             (initial_pot_es(1,i_type), n_grid(comp_species), &
             initial_pot_es_spl(1,1,i_type) )
        
        call localorb_info('',use_unit)
        write(info_str,'(2X,A)') &
             "List of occupied orbitals and eigenvalues: "
        call localorb_info(info_str,use_unit,'(A)')
        
        write(info_str,'(4X,A,4X,A,6X,A,4X,A,4X,A)') &
             "n", "l", "energy [Ha]", "energy [eV]","occupation"
        call localorb_info(info_str,use_unit,'(A)')
        
        write(info_str,'(2X,A)') "Spin-up "
        call localorb_info(info_str,use_unit,'(A)')
        do i_function = 1, n_atomic(comp_species), 1
           write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.8)') &
                n_shell(i_function), l_shell(i_function), &
                eigenval(i_function), &
                eigenval(i_function) * hartree, occ_shell(i_function)
           call localorb_info(info_str,use_unit,'(A)')
        end do
        
        if (n_spin .eq. 2) then
           write(info_str,'(2X,A)') "Spin-down "
           call localorb_info(info_str,use_unit,'(A)')
           do i_function = n_atomic(comp_species) + 1, 2*n_atomic(comp_species), 1
              write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.8)') &
                   n_shell(i_function), l_shell(i_function), &
                   eigenval(i_function), &
                   eigenval(i_function) * hartree, occ_shell(i_function)
              call localorb_info(info_str,use_unit,'(A)')
           end do
        end if
        
        call localorb_info('',use_unit)
        
     else
        write(info_str,'(2X,A)') "No electronic charge for this initialization-type. "
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(2X,A)') "Skipping atomic solver. "
        call localorb_info(info_str,use_unit,'(A)')
     end if
     
  end do

  deallocate(n_shell)
  deallocate(l_shell)
  deallocate(m_shell)
  deallocate(occ_shell)
  deallocate(eigenval)
  deallocate(i_grid_max)
  deallocate(r_peak)
  deallocate(core_type)
  
  deallocate(ini_pot)
  deallocate(atom_pot)
  deallocate(atom_pot_xc)
  deallocate(density)
  deallocate(density_deriv)
  deallocate(density_2nd_deriv)
  deallocate(core_density)
  deallocate(core_dd)
  deallocate(core_d2d)
  
  deallocate(wave)
  
  deallocate(e_sic)
  deallocate(v_orb)
  deallocate(sv_kli)
  deallocate(svm)
  deallocate(wave_deriv)
  
  return

end subroutine get_free_atoms_polarized
!******  
