!****s* FHI-aims/get_free_atoms
!  NAME
!   get_free_atoms
!  SYNOPSIS

 subroutine get_free_atoms (  )

!  PURPOSE
!  Subroutine get_free_atom provides atomic input data for each species
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
      use constants,       only : hartree, pi4
      use dimensions,      only : n_max_grid, n_species, n_channels, n_max_ind_fns, &
                                  use_basis_gradients, use_density_gradient, &
                                  use_partition_deriv, use_relativistic_basis
      use runtime_choices, only : add_embedding_grids, analytic_potential_average, &
                                  atomic_solver_atom_sphere, atomic_solver_sratom, &
                                  atomic_solver, flag_atomic_xc, flag_KH_core_states, &
                                  flag_rel, flag_xc, hse_omega, hybrid_coeff, &
                                  multipole_radius_free_threshold, REL_at_zora_spinor, &
                                  REL_atomic_zora, REL_KOLNING_HARMON, REL_x2c, &
                                  REL_4c_dks, REL_zora, REL_zora_spinor, REL_own
      use grids,           only : r_grid, n_grid, r_grid_inc
      use spline,          only : cubic_spline
      use species_data,    only : multipole_radius_free, multipole_radius_free_sq, &
                                  species_z, species_name, n_atomic, atomic_n, &
                                  atomic_l, atomic_k, valence_occ, core_fn, cut_free_atom, &
                                  cutoff_type, free_r_cut, scale_cutoff, w_cutoff, &
                                  species_pseudoized, no_basis
      use rel_x2c_mod,     only : atom_occ_shell, free_den_diff, free_den_diff_spl, free_drho_dr_diff_spl
      use free_atoms,      only : free_potential, free_pot_es, free_rho, free_wave, free_wave_small, &
                                  free_kinetic, free_wave_eigenval, free_wave_deriv, free_wave_small_deriv, &
                                  renormalized_free_rho_spl, renormalized_free_drho_dr_spl, &
                                  free_pot_es_at_zero,  free_atom_average_es_pot, &
                                  free_rho_spl, free_drho_dr_spl, free_atom_average_es_vol, &
                                  free_pot_es_spl, free_potential_spl, free_wave_spl, &
                                  free_kinetic_spl, free_wave_deriv_spl
      use localorb_io,     only : OL_norm, use_unit, localorb_info
      use xc_library,      only : xc__pw_lda, xc_sr__pw_lda, xc__pbe, xc_sr__pbe, &
                                  xc__pbe0, is_dfauto, normalize_flag_xc
      use synchronize_mpi, only : sync_atomic_solver_outputs
      use mpi_tasks,       only : aims_stop, aims_stop_coll, check_allocation
      use psi_at_nucleus_mod, only: psi_at_nucleus_atomic
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

!  local variables

      real*8 drho_dr(n_max_grid, n_species)
      real*8 d2rho_dr2(n_max_grid, n_species)
      real*8 drho_dr_diff(n_max_grid, n_species) ! for fully relativistic cases

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

      real*8 nuclear_charge
      real*8 ::  dr_coef, alpha

      integer n_core
      integer n_valence

      integer, dimension(:), allocatable :: n_shell
      integer, dimension(:), allocatable :: l_shell
      integer, dimension(:), allocatable :: k_shell
      integer, dimension(:), allocatable :: m_shell
      real*8, dimension(:), allocatable :: occ_shell
      real*8, dimension(:), allocatable :: eigenval
      integer, dimension(:), allocatable :: i_grid_max
      logical, dimension(:), allocatable :: core_type

      real*8 d_eig_sum

      real*8, dimension(:), allocatable :: r_peak

      real*8, dimension(:), allocatable :: ini_pot
      real*8, dimension(:,:), allocatable :: atom_pot
      real*8, dimension(:,:), allocatable :: atom_pot_es
      real*8, dimension(:,:), allocatable :: atom_pot_xc
      real*8, dimension(:,:), allocatable :: density
      real*8, dimension(:,:), allocatable :: density_deriv
      real*8, dimension(:,:), allocatable :: density_2nd_deriv
      real*8, dimension(:,:), allocatable :: density_diff
      real*8, dimension(:,:), allocatable :: density_diff_deriv
      real*8, dimension(:,:), allocatable :: core_density
      real*8, dimension(:,:), allocatable :: core_dd
      real*8, dimension(:,:), allocatable :: core_d2d

      real*8, dimension(:,:), allocatable :: wave
      real*8, dimension(:,:), allocatable :: kinetic
      real*8, dimension(:,:), allocatable :: small_wave

      real*8 :: outermost_charge

      logical t_in_pot 
      logical t_sic
      real*8, dimension(:), allocatable :: e_sic
      real*8, dimension(:,:), allocatable :: v_orb
      real*8, dimension(:,:), allocatable :: spline_coef

      logical t_kli
      real*8, dimension(:), allocatable :: sv_kli
      real*8, dimension(:), allocatable :: svm

      logical t_mgga

      real*8, dimension(:,:), allocatable :: wave_deriv
      real*8, dimension(:,:), allocatable :: small_deriv
      

!  other aux variables

      integer n_max

      character*100 :: info_str

!  counters

      integer i_species
      integer i_grid
      integer i_channels
      integer i_function, info
                        ! test:
                        real*8 :: int_log_mesh
                        real*8 :: scalar

!  begin work

      write (info_str,'(2X,A,A)') &
        "Creating wave function, potential, and density ", &
        "for free atoms."


      call localorb_info('',use_unit)
      call localorb_info(info_str,use_unit,'(A)')

!      write(use_unit,*) 
!      write(use_unit,*) " !reating wave function, potential, and density ",
!     +  "Fo free atoms."

!  allocations - needed because some of these arrays may collide with the
!  stack size on small machines if they were automatic arrays

      allocate(n_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'n_shell                       ') 

      allocate(l_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'l_shell                       ') 
      
      allocate(k_shell(n_max_ind_fns),stat=info)
      call check_allocation(info, 'k_shell                       ') 

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

      allocate(atom_pot_es(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'atom_pot_es                   ') 

      allocate(atom_pot_xc(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'atom_pot_xc                   ') 

      allocate(density(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'density                       ') 

      allocate(density_diff(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'density_diff                  ') 

      allocate(density_diff_deriv(n_max_grid, n_channels),stat=info)
      call check_allocation(info, 'density_diff_deriv            ') 

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
      
      allocate(kinetic(n_max_grid, n_max_ind_fns), stat=info)
      call check_allocation(info, 'kinetic                       ') 

      allocate(small_wave(n_max_grid, n_max_ind_fns),stat=info)
      call check_allocation(info, 'small_wave                    ')

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

      allocate(small_deriv(n_max_grid, n_max_ind_fns),stat=info)
      call check_allocation(info, 'small_deriv                   ')

      allocate(spline_coef(4, n_max_grid),stat=info)
      call check_allocation(info, 'spline_coef                   ')

!  input parameters for sr_atom(), from control.in 

      ! blanket setting for XC, if basis fn XC is the same as eventual self-consistent XC
      i_exch = flag_xc


      ! Set exchange-correlation option for sratom
      if (flag_atomic_xc.eq.12) then
         ! This means that the atomic exchange flag was specifically set in control.in
         ! This choice is therefore essential
         i_exch = flag_atomic_xc

      else 

         !  for Hartree-Fock calculation --Ren
         if(flag_xc .eq. 0) then
            ! default is set in read_control.f90 - as of this writing (August 2017), the default is 8, that is, pw-lda.
            i_exch = flag_atomic_xc
         elseif(flag_xc .eq. 1) then
         !  Hybrid-PBE0 calculation, using PBE instead as the atomic solver.
            i_exch = 6
         endif

         if(flag_xc .eq. 10) then
         !  Hybrid-B3LYP calculation, using BLYP instead as the atomic solver.
            i_exch = 9
         endif

         ! HSE calculation or LC-wPBEh, using PBE
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

         !PBE_vdw.   SAG
         if (flag_xc.eq.5) then
            i_exch = 6
         end if

         !R48PBE
         if (flag_xc.eq.18) then
            i_exch = 6
         endif

         !xPBE    IYZ
         if (flag_xc.eq.22) then
            i_exch = 8
         endif


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

         ! Caveat to catch all Meta-GGAs
         ! We'll solve calculate the free atom for them all as PW91-LDA
         ! similar to the approach above used for AM05. AJL
         if (flag_xc.ge.25) then
            i_exch = 8 
         endif

         ! Caveat to catc LibXC
         ! We'll solve calculate the free atom for them all as PW91-LDA
         ! similar to the approach above used for AM05. AJL
         if (flag_xc.le.-10) then
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
      end if ! Set exchange-correlation option for atomic solver


! LDA's, no gradient derivatives needed   
      if ((i_exch.eq.3).or.(i_exch.eq.8).or.(i_exch.eq.7)) then
        i_grad = 0
!     PBE and its relatives or BLYP, should calculate the derivatives for free-atom
!@@edu>
!      else if ((i_exch.eq.6).or.(i_exch.eq.9).or.(i_exch.eq.14).or.(i_exch.eq.15)) then
      else if ((i_exch.eq.4).or.(i_exch.eq.6).or.(i_exch.eq.9) &
           .or.(i_exch.eq.14).or.(i_exch.eq.15).or.(i_exch.eq.21)) then
!@@edu<
        i_grad = 1
      end if

      if (use_density_gradient) then
        i_grad = 1
      end if

     !---------shanghui add for free_drho_dr_spl for gradient_partition----
      if (use_partition_deriv) then 
        i_grad = 1
      end if
     !---------shanghui end add for free_drho_dr_spl for gradient_partition----
      if(flag_KH_core_states)then
         i_rel = 10
      elseif (flag_rel.eq.0) then
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
      t_kli = .false.
      t_mgga = .false.

      if(flag_xc.eq.0) then
       if(flag_atomic_xc.eq.12) then
          t_kli = .true.
       endif
      endif

!  initialize grids
!  initialize also the core density parts, currently unused.
!  FIXME: This should really be done f90 style!

      do i_species = 1, n_species, 1
        do i_grid = 1, n_max_grid, 1
          free_potential(i_grid, i_species)   = 0.0d0          
          free_pot_es   (i_grid, i_species)   = 0.0d0          
          free_rho(i_grid, i_species) = 0.0d0 
          drho_dr(i_grid, i_species) = 0.0d0          
          d2rho_dr2(i_grid, i_species) = 0.0d0          
          drho_dr_diff(i_grid, i_species) = 0.0d0
          do i_channels = 1,n_channels, 1
            core_density(i_grid,1) = 0.d0
            core_dd(i_grid,1) = 0.d0
            core_d2d(i_grid,1) = 0.d0
          enddo
          do i_function = 1, n_max_ind_fns, 1
            free_wave(i_grid,i_species,i_function) = 0.d0
            free_kinetic(i_grid,i_species,i_function) = 0.d0
            if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) free_wave_small(i_grid,i_species,i_function) = 0.d0
          enddo 
        enddo
      enddo

      do i_species = 1,n_species,1

        write(info_str,'(2X,A,A)') "Species: ", species_name(i_species)
        call localorb_info('',use_unit)
        call localorb_info(info_str,use_unit,'(A)')

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
          if(use_relativistic_basis) then
              k_shell (i_function) = atomic_k(i_species,i_function)
          end if
          m_shell (i_function) = 1

          if(use_relativistic_basis) then
            if (l_shell(i_function).gt.0) then
                if (valence_occ(n_shell(i_function),l_shell(i_function),i_species).lt.1.d-8) then
                    occ_shell(i_function) = valence_occ(n_shell(i_function),l_shell(i_function),i_species)
                ! Check whether the shell is closed
                elseif ( dabs( valence_occ(n_shell(i_function),l_shell(i_function),i_species) - 2*(2*l_shell(i_function)+1) ) .lt. 1.d-8 ) then
                    if(k_shell(i_function).gt.0) then
                        occ_shell(i_function) = 2*l_shell(i_function)
                    else
                        occ_shell(i_function) = 2*(l_shell(i_function)+1)
                    end if
                ! If the shell is open, divide electrons correctly to two split spin orbitals. 
                ! One scheme is: all electrons should firstly go to the larger-kappa orbital, then the smaller-kappa orbital.
                ! Here we adopt the most popular average occupation scheme to keep the atomic electron density spherical.
                else
                    if(k_shell(i_function).gt.0) then
                        occ_shell(i_function) = dble(2*l_shell(i_function))/dble(2*(2*l_shell(i_function)+1)) * valence_occ(n_shell(i_function),l_shell(i_function),i_species)
                    else
                        occ_shell(i_function) = dble(2*(l_shell(i_function)+1))/dble(2*(2*l_shell(i_function)+1)) * valence_occ(n_shell(i_function),l_shell(i_function),i_species)
                    end if
               !else
               !    if( dabs( valence_occ(n_shell(i_function),l_shell(i_function),i_species) - 2*(l_shell(i_function)) ) .lt. 1.d-8 ) then                        
               !        if(k_shell(i_function).gt.0) then
               !            occ_shell(i_function) = valence_occ(n_shell(i_function),l_shell(i_function),i_species)
               !        else
               !            occ_shell(i_function) = 0.0
               !        end if
               !    else
               !        if(k_shell(i_function).gt.0) then
               !            occ_shell(i_function) = 2*l_shell(i_function)
               !        else
               !            occ_shell(i_function) = valence_occ(n_shell(i_function),l_shell(i_function),i_species) - 2*l_shell(i_function)
               !        end if
               !    end if
                end if 
            else ! (l_shell(i_function).le.0), for s orbital, there is no splitting due to spin
                occ_shell(i_function) = valence_occ(n_shell(i_function),l_shell(i_function),i_species)
            end if
          else ! non-rel
              occ_shell(i_function) = valence_occ(n_shell(i_function),l_shell(i_function),i_species)
          end if

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

!       FIXME: This would be a superfluous copy, if atom_pot were not grossly overdimensioned
        do i_grid = 1, n_grid(i_species), 1
          atom_pot(i_grid, 1) = ini_pot(i_grid)
          atom_pot_xc(i_grid, 1) = 0.d0
          atom_pot_es(i_grid, 1) = 0.d0
        enddo

!       Select correct basis generator for relativistic calculations
!        
        if (flag_rel.eq.REL_at_zora_spinor .or. flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks .or. flag_rel.eq.REL_zora_spinor) then 
!          if we take account spin-orbit coupling or use X2C, use basis set from the dirac equation
            call rel_atoms(i_species, nuclear_charge, n_channels, n_atomic(i_species), n_shell, k_shell, l_shell, occ_shell, &
            n_grid(i_species), r_grid(1,i_species), .false., cutoff_type(i_species), free_r_cut(i_species),&
            scale_cutoff(i_species), w_cutoff(i_species),1,i_exch, &
!           outputs  
            atom_pot, atom_pot_es, density, density_deriv, density_2nd_deriv, density_diff, &
            wave, wave_deriv, small_wave, small_deriv, eigenval)      

            atom_occ_shell(:,i_species) = occ_shell

            call fderiv(1,n_grid(i_species),r_grid(1,i_species),density_diff,density_diff_deriv,spline_coef)
        else 
!       otherwise use old atomic solver
!
             ! WPH: Here we run the chosen atomic solver to generate candidate basis functions
             ! Should you wish to add in another solver, the values extracted from these solvers are:
             !      eigenval, wave, atom_pot, atom_pot_es, density, density_deriv, density_2nd_deriv
             ! defined on the radial grid specified by r_grid (note that the l_channel index in these 
             ! variables is legacy code that is not needed)
             ! 
             ! If you wish to bypass the recalculation of basis functions during get_species_basis_fns,
             ! (read the comments for the keyword for use_atomic_solver_fns_as_minimal_basis for more 
             ! details,) you will also need to provide the following variables
             !      wave_deriv, kinetic
             ! yourself.
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

               call atom_sphere_wrapper &
               ( n_grid(i_species), r_grid(1,i_species), species_name(i_species), flag_xc, hybrid_coeff, hse_omega, &
                 nuclear_charge, n_core, n_valence, n_max_ind_fns, n_shell, l_shell, m_shell, occ_shell, n_max_grid, & 
                 n_channels, 1, &
                 cut_free_atom(i_species), cutoff_type(i_species), free_r_cut(i_species), scale_cutoff(i_species), &
                          w_cutoff(i_species), & ! Cutoff potential
                 eigenval, wave, wave_deriv, density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
                 kinetic, psi_at_nucleus_atomic(:n_valence,i_species)) ! Output
               call sync_atomic_solver_outputs( n_max_ind_fns, n_max_grid, n_channels, 1, &
                    eigenval, wave, wave_deriv, kinetic, &
                    density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, psi_at_nucleus_atomic(:n_valence,i_species))
             else
               call aims_stop('The atomic solver specified by the variable atomic_solver was not recognized, exiting.','func')  
            end if
        end if 


!       now copy all data from atom onto storage arrays

!       atom_pot and density have an extra index which used to be shell-dependent (Hatree-Fock?),
!       then indicated spin up / down (Martin Fuchs), and is unused here. 

        do i_function = 1, n_atomic(i_species), 1
          free_wave_eigenval(i_species,i_function) = &
            eigenval(i_function)
        enddo

        ! Note:  not all atomic solvers will be outputting meaningful data for
        ! these variables!  Be careful before using them in your own code.
        do i_grid = 1, n_grid(i_species), 1
          free_potential(i_grid, i_species)   = atom_pot(i_grid, 1) 
          free_pot_es(i_grid, i_species)   = atom_pot_es(i_grid, 1) 
          free_rho(i_grid, i_species) = density(i_grid, 1)
          drho_dr(i_grid, i_species) = density_deriv(i_grid, 1)
          d2rho_dr2(i_grid, i_species) = density_2nd_deriv(i_grid, 1)
          if(flag_rel.eq.REL_4c_dks)then
             free_den_diff(i_grid, i_species) = density_diff(i_grid, 1)
             drho_dr_diff(i_grid, i_species) = density_diff_deriv(i_grid, 1)
          endif
          do i_function = 1, n_atomic(i_species), 1
            free_wave (i_grid, i_species, i_function) = wave(i_grid, i_function)
            free_kinetic (i_grid, i_species, i_function) = kinetic(i_grid, i_function)
            if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
              free_wave_small (i_grid, i_species, i_function) = small_wave(i_grid, i_function)
              free_wave_deriv (i_grid, i_species, i_function) = wave_deriv(i_grid, i_function)
              free_wave_small_deriv (i_grid, i_species, i_function) = small_deriv(i_grid, i_function)
            endif
            if (use_basis_gradients) then
              free_wave_deriv (i_grid, i_species, i_function) = wave_deriv(i_grid, i_function)
            end if
          enddo
        enddo

!       spline free waves
        do i_function = 1,n_atomic(i_species),1
          call cubic_spline &
          ( free_wave(1,i_species,i_function),n_grid(i_species), &
            free_wave_spl(1,1,i_function,i_species) &
          )
          call cubic_spline &
          ( free_kinetic(1,i_species,i_function),n_grid(i_species), &
            free_kinetic_spl(1,1,i_function,i_species) &
          )
          if (use_basis_gradients) then
            call cubic_spline &
            ( free_wave_deriv(1,i_species,i_function),n_grid(i_species), &
              free_wave_deriv_spl(1,1,i_function,i_species) &
            )
          end if
        enddo

!       spline potentials
        call cubic_spline &
        ( free_potential(1,i_species), n_grid(i_species), &
          free_potential_spl(1,1,i_species) ) 
        call cubic_spline &
        ( free_pot_es(1,i_species), n_grid(i_species), &
          free_pot_es_spl(1,1,i_species) )

!       spline rho
        call cubic_spline &
        ( free_rho(1,i_species), n_grid(i_species),  &
          free_rho_spl(1,1,i_species) )
        ! Additional copy for use in code parts that concern electrostatics.
        ! For the benefit of improved electrostatics, we may renormalize
        ! the overall free-atom density slightly to keep its norm on the
        ! 3D integration grid.
        renormalized_free_rho_spl(:,:,i_species) = free_rho_spl(:,:,i_species)

        if(flag_rel.eq.REL_4c_dks)then
          call cubic_spline &
          ( free_den_diff(1,i_species), n_grid(i_species), &
          free_den_diff_spl(1,1,i_species) )
          if (use_density_gradient) then ! spline drho_dr_diff
             call cubic_spline( drho_dr_diff(1,i_species), n_grid(i_species), &
                  free_drho_dr_diff_spl(1,1,i_species) )
          endif
        endif
        

        if (use_density_gradient) then
!         spline drho_dr
          call cubic_spline &
          ( drho_dr(1,i_species), n_grid(i_species), &
            free_drho_dr_spl(1,1,i_species) )
          ! see comment above
          renormalized_free_drho_dr_spl(:,:,i_species) = free_drho_dr_spl(:,:,i_species)
        end if

    !---------shanghui add for free_drho_dr_spl for gradient_partition----
       if (use_partition_deriv) then !shanghui add 
          call cubic_spline &
          ( drho_dr(1,i_species), n_grid(i_species), &
            free_drho_dr_spl(1,1,i_species) )
        end if
     !---------shanghui end add for free_drho_dr_spl for gradient_partition----



        call localorb_info('',use_unit)
        write(info_str,'(2X,A)') "List of occupied orbitals and eigenvalues: "  
        call localorb_info(info_str,use_unit,'(A)')

        if (flag_rel.eq.REL_at_zora_spinor.or.flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks.or.flag_rel.eq.REL_zora_spinor) then 
            
            write(info_str,'(2X,A3,2X,A3,2X,A3,2X,A15,2X,A15,A15)') "n", "l", "k", "occ", "energy [Ha]", "energy [eV]"
            call localorb_info(info_str,use_unit,'(A)')

            do i_function = 1, n_atomic(i_species), 1
                write(info_str,'(2X,I3,2X,I3,2X,I3,2X,F15.4,2X,F15.6,F15.4)')&
                    atomic_n(i_species,i_function), &
                    atomic_l(i_species,i_function), &
                    atomic_k(i_species,i_function), &
                    occ_shell(i_function), &
                    free_wave_eigenval(i_species, i_function), &
                    free_wave_eigenval(i_species, i_function)*hartree
                call localorb_info(info_str,use_unit,'(A)')
            enddo
        else
            write(info_str,'(2X,A3,2X,A3,2X,A15,2X,A15,A15)') "n", "l", "occ", "energy [Ha]", "energy [eV]"
            call localorb_info(info_str,use_unit,'(A)')

            do i_function = 1, n_atomic(i_species), 1
                write(info_str,'(2X,I3,2X,I3,2X,F15.4,2X,F15.6,F15.4)') &
                    atomic_n(i_species,i_function), &
                    atomic_l(i_species,i_function), &
                    occ_shell(i_function), &
                    free_wave_eigenval(i_species, i_function), &
                    free_wave_eigenval(i_species, i_function)*hartree
                call localorb_info(info_str,use_unit,'(A)')
            enddo
        end if
        call localorb_info('',use_unit)

        ! Finally: Find the radius outside of which all free-atom quantities
        ! (namely the density) are effectively zero for this atom.
        !
        i_grid = n_grid(i_species)
        do while ( (dabs(free_rho(i_grid, i_species)).le. &
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

        multipole_radius_free(i_species) = outermost_charge
        multipole_radius_free_sq(i_species) = outermost_charge**2

        ! Electrostatic potential limit at position of nucleus "zero limit".
        free_pot_es_at_zero(i_species) = 0.d0
        alpha = log(r_grid_inc(i_species))

        do  i_grid = 1,  n_grid(i_species),1 
           
          dr_coef  = r_grid(i_grid,i_species)**2 * alpha

          free_pot_es_at_zero(i_species)=free_pot_es_at_zero(i_species) &
                + free_rho( i_grid,i_species)  *  dr_coef
        end do

        if (analytic_potential_average) then
           ! use free_pot_es as a basis to determine the average electrostatic potential
           ! of each atom in the volume within which it is non-zero
           free_atom_average_es_pot(i_species) = 0.d0
           free_atom_average_es_vol(i_species) = 0.d0
           if (.not.species_pseudoized(i_species)) then
              do i_grid = 1, n_grid(i_species), 1
                 if (free_pot_es(i_grid,i_species).ne.0.d0) then
                    dr_coef  = r_grid(i_grid,i_species)**3 * alpha ! radial integration weight on logarithmic grid
                    free_atom_average_es_pot(i_species) = free_atom_average_es_pot(i_species) + &
                      dr_coef * free_pot_es(i_grid,i_species)
                    free_atom_average_es_vol(i_species) = free_atom_average_es_vol(i_species) + &
                      dr_coef 
                 end if
              enddo
              free_atom_average_es_pot(i_species) = free_atom_average_es_pot(i_species) * pi4
              free_atom_average_es_vol(i_species) = free_atom_average_es_vol(i_species) * pi4
           end if
        end if

!FIXME: DB: 22.11.: this is a first trial&error approach to figure out which arrays needs to be 
!                   set to zero, if (pseudoized(i_species))

        if (species_pseudoized(i_species)) then

           free_pot_es(:,i_species) = 0.d0 
           free_pot_es_spl(:,:,i_species) = 0.d0
           free_pot_es_at_zero(i_species) = 0.d0 

           free_rho(:,i_species) = 0.d0
 ! when no grids are wanted around pseudocores, we set their free density properties to zero
 !    this will result in a zero partition function
           if(.not.add_embedding_grids) then
              free_rho_spl(:,:,i_species) = 0.d0
              if(allocated(free_drho_dr_spl)) then
                free_drho_dr_spl(:,:,i_species) = 0.d0
              end if
           endif

           drho_dr(:,i_species) = 0.d0
           renormalized_free_rho_spl(:,:,i_species) = 0.d0

           if(allocated(renormalized_free_drho_dr_spl)) then
             renormalized_free_drho_dr_spl(:,:,i_species) = 0.d0
           end if

           free_potential(:,i_species) = 0.d0 
           free_potential_spl(:,:,i_species) = 0.d0
           
           species_z(i_species) = 0.d0

        end if

!DB: 092612: the free_pot_es_spl leads to spurious hellman-feynmann forces for every empty atom.
!            this cannot be done here, since it concerns specific atoms and not species.
!            should be done by simply cycling over empty atoms
! only for debugging
      if (no_basis(i_species).or.species_pseudoized(i_species)) then
          free_pot_es_spl(:,:,i_species) = 0.d0
!          free_rho_spl(:,:,i_species) = 0.d0
           free_potential(:,i_species) = 0.d0
           free_potential_spl(:,:,i_species) = 0.d0

       endif
!    100112: this is now done directly in sum_up_whole_potential_p1()
    
      enddo


      if (analytic_potential_average) then
         call localorb_info('',use_unit,'(A)', OL_norm)
         write (info_str,'(2X,A)') &
           "Analytic average potentials were requested. Free-atom quantities:"
         call localorb_info(info_str,use_unit,'(A)', OL_norm)
         write (info_str,'(2X,A)') &
           "Species          Avg. potential (a.u.)      Volume (a.u.)"
         call localorb_info(info_str,use_unit,'(A)', OL_norm)

         do i_species = 1, n_species, 1
            write (info_str,'(2X,I8,15X,F15.8,4X,F15.8)') &
              i_species, free_atom_average_es_pot(i_species), free_atom_average_es_vol(i_species)
            call localorb_info(info_str,use_unit,'(A)', OL_norm)
         enddo
         call localorb_info('',use_unit,'(A)', OL_norm)

      end if

!     deallocate local variables

      deallocate(wave_deriv)
      deallocate(small_deriv)
      deallocate(spline_coef)
      deallocate(svm)
      deallocate(sv_kli)
      deallocate(v_orb)
      deallocate(e_sic)
      deallocate(wave)
      deallocate(small_wave)
      deallocate(core_d2d)
      deallocate(core_dd)
      deallocate(core_density)
      deallocate(density_2nd_deriv)
      deallocate(density_deriv)
      deallocate(density)
      deallocate(density_diff)
      deallocate(density_diff_deriv)
      deallocate(atom_pot_xc)
      deallocate(atom_pot)
      deallocate(ini_pot)
      deallocate(core_type)
      deallocate(r_peak)
      deallocate(i_grid_max)
      deallocate(eigenval)
      deallocate(occ_shell)
      deallocate(m_shell)
      deallocate(l_shell)
      deallocate(n_shell)

 end subroutine get_free_atoms

