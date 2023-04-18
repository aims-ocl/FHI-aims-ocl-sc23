!  FIXME: Subroutine get_free_ion is really a cleaned-up version of get_free_atoms,
!         only for one atom, not for many.
!  FIXME: With very little generalisation, it can be used directly in get_free_atoms:
!         must rename / clean up species_z, ion_n_max, ion_occ, n_grid, r_grid
!         EXCEPT WE DO NOT WANT TO PASS ON the electrostatic pot, XC pot, density/derivatives!
!

!****s* FHI-aims/get_free_ion
!  NAME
!   get_free_ion
!  SYNOPSIS

 subroutine get_free_ion ( i_species, &
   free_ion_pot, n_free_ionic, free_ion_wave, &
   free_ion_wave_deriv, free_ion_eigenval, free_ion_kinetic, free_ion_n, free_ion_l )

!  PURPOSE
!    Subroutine get_free_ion solves the self-consistent Kohn-Sham equations
!    for a free ion, much like get_free_atom.
! USES
   use constants,       only : hartree
   use dimensions,      only : n_max_grid, n_max_ind_fns, n_channels
   use runtime_choices, only : atomic_solver, atomic_solver_sratom, &
                               atomic_solver_atom_sphere, flag_KH_core_states, flag_rel, &
                               flag_xc, hse_omega, hybrid_coeff, REL_atomic_zora, &
                               REL_kolning_harmon, REL_own, REL_zora, REL_x2c, REL_4c_dks
   use grids,           only : n_grid, r_grid
   use species_data,    only : species_name, species_z, ion_occ, core_n_max, ion_n_max, &
                               l_shell_max, cutoff_type, scale_cutoff, w_cutoff
   use xc_library,      only : xc__pw_lda, xc_sr__pw_lda, xc__pbe, xc_sr__pbe, &
                               xc__pbe0, is_dfauto, normalize_flag_xc
   use localorb_io,     only : use_unit, localorb_info
   use mpi_tasks,       only : aims_stop
   use synchronize_mpi, only : sync_atomic_solver_outputs
   use psi_at_nucleus_mod, only: psi_at_nucleus_ionic
   implicit none
! ARGUMENTS

   integer :: i_species
   real*8  :: free_ion_pot(n_max_grid)
   integer :: n_free_ionic
   real*8  :: free_ion_wave(n_max_grid,n_max_ind_fns)
   real*8  :: free_ion_wave_deriv(n_max_grid,n_max_ind_fns)
   real*8  :: small_wave(n_max_grid,n_max_ind_fns) ! small component part (for relativistic case)
   real*8  :: small_deriv(n_max_grid,n_max_ind_fns)
   real*8  :: free_ion_eigenval(n_max_ind_fns)
   real*8  :: free_ion_kinetic(n_max_grid,n_max_ind_fns)
   integer :: free_ion_n(n_max_ind_fns)
   integer :: free_ion_l(n_max_ind_fns)
   integer :: free_ion_k(n_max_ind_fns)

! INPUTS 
! o i_species -- the atom species index
!
! OUTPUTS
! o free_ion_pot -- radial potential of the free ion
! o n_free_ionic -- number of free ion functions
! o free_ion_wave -- wave functions
! o free_ion_eigenval -- corresponding eigenvalue
! o free_ion_n -- quantum number n for free ion functions
! o free_ion_l -- quantum number l for free ion functions
!
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

   ! This should be exported at some point. For now assume that
   ! all ionic basis functions are valence functions.
   logical core_type(n_max_ind_fns)

!  nuclear_charge : well, the nuclear charge ...
!  i_rel    : Masked version of flag_rel; i_rel = 1 <-> non-rel., i_rel=2 <-> Koelling-Harmon
!  i_grad   : if i_grad=1, calculate free-atom density derivatives (1st and 2nd)
!  sc_max_iter : Maximum number of self-consistency iterations for free atom
!  t_core   : Here, always true - we want to calculate core wavefunctions!
!  n_core   : number of frozen core shells in a pseudopot'l calculation?? This is not needed here!
!  t_in_pot : if true, sratom() assumes an input non-selfconsistent potential - always false for now.
!  t_sic :    if true, sratom() does self-interaction correction - this is never the case.
!  t_kli :    kli is not used
!  t_mgga :   meta-gga used? currently never.
!  n_shell  : 1D list of main quantum numbers of occupied shells
!  l_shell  : 1D list of angular momenta of occupied shells
!  m_shell  : This value seems to label the nature of each state - spin up (m=1) or spin down (m=2)
!             for each shell
!  occ_shell : occupation number of each shell
!  n_valence: number of "valence" (n,l) shells in the calculation
!  ini_pot :   initial radial potential for each individual atom.
!  atom_pot :  s.-c. atomic potential FIXME: THIS HAS ONE INDEX TOO MANY - WE DO NOT NEED THE
!              max_shells INDEX FOR ANYTHING, IT IS A HARTREE_FOCK RELIC!
!  i_grid_max: integer index of maximum grid point in integration for each shell
!  d_eig_sum : FIXME This is an apparently unnecessary variable which monitors the
!              change of the eigenvalue sum between successive scf iterations
!  r_peak    : FIXME: Unused. Outermost peak position of each individual wave function.
!  density, density_deriv, density_2nd_deriv : radial atomic density and its derivatives.
!              FIXME: SAME AS atom_pot - THE max_shells INDEX IS COMPLETELY USELESS!
!  core_density, core_dd, core_d2d : UNUSED - core density, density derivative, and second
!              derivative, for frozen core calculations (never done here)
!  wave :      s.-c. atomic wave functions
!  e_sic :    array of self-energies for SIC - not used
!  v_orb :    orbital-depedent potential(?) for sic - not used
!  sv_kli :   kli is not used
!  svm   :    kli is not used
!  wave_deriv : radial derivative of each wave function - might be useful after all.

   integer :: i_exch
   real*8 nuclear_charge
   integer i_rel
   integer i_grad
   integer sc_max_iter
   logical t_core
   integer n_core
   logical t_in_pot
   logical t_sic
   logical t_kli
   logical t_mgga
   integer m_shell(n_max_ind_fns)
   real*8  occ_shell(n_max_ind_fns)
   integer n_valence
   real*8 ini_pot(n_max_grid)
   real*8 atom_pot(n_max_grid, n_channels)
   real*8 atom_pot_es(n_max_grid)
   real*8 atom_pot_xc(n_max_grid, n_channels)
   integer i_grid_max(n_max_ind_fns)
   real*8 d_eig_sum
   real*8 r_peak(n_max_ind_fns)
   real*8 density(n_max_grid, n_channels)
   real*8 density_diff(n_max_grid, n_channels)
   real*8 density_deriv(n_max_grid, n_channels)
   real*8 density_2nd_deriv(n_max_grid, n_channels)
   real*8 core_density(n_max_grid,n_channels)
   real*8 core_dd(n_max_grid,n_channels)
   real*8 core_d2d(n_max_grid,n_channels)
   real*8 e_sic(n_max_ind_fns)
   real*8 v_orb(n_max_grid,n_max_ind_fns)
   real*8 sv_kli(n_channels)
   real*8 svm(n_channels)
   real*8 wave_deriv(n_max_grid, n_max_ind_fns)

   ! cut_free_ion and ion_r_cut are placeholder variables until their final use is decided.
   ! IF we were ever to initialize the charge density of the first iteration with a
   ! free ion anywhere (e.g. a Li+ ion in a place where we know there is a Li+),
   ! THEN we would need finite values for these radii - probably the same as for
   ! the free atom.
   ! If not, the only possible effect of these variables could be to modify the shape
   ! of ionic basis functions, which we do not want.
   logical cut_free_ion
   real*8 :: ion_r_cut


   character*100 :: info_str

!  counters

   integer i_grid, i_channel, i_l, i_shell, i_fn

!  functions

!  begin work

   i_exch = flag_XC

!  for Hartree-Fock calculation
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
!  HSE calculation or lc-wpbeh, using PBE
   if (flag_xc .eq.7.or.flag_xc.eq.23) then
     i_exch = 6
   endif
! PBEsol  and PBEsol0 calculation, using PBE as atomic solver
   if (flag_xc .eq. 17 .or. flag_xc .eq. 13) then
      i_exch = 6
   end if
! RPBE, revPBE calculation, using PBE instead as the
! atomic solver.  Or revpbe_vdw! SAG. 
   if (flag_xc.eq.11.or.flag_xc.eq.12.or.flag_xc.eq.14) then
      i_exch = 6
   end if
   !48PBE
   if (flag_xc.eq.18) then
      i_exch = 6
   end if
!PBE_vdw. SAG
   if (flag_xc.eq.5) then
      i_exch = 6
   end if
!xPBE     IYZ
   if (flag_xc.ge.22) then
      i_exch = 8
   endif
!VWN
   if (flag_xc.eq.15) then
      i_exch = 7
   endif


! FIXME : A case should be added in external/vexcor.f in order to
! calculate the AM05 energies and potentials (using subroutine
! external/am05.f) and use that as atomic solver!
! AM05 calculation, using PBE instead as the atomic solver
   if (flag_xc.eq.20) then
      i_exch = 6
   end if

! Caveat to catch all Meta-GGAs for now.
! We'll solve calculate the free atom for them all as PW91-LD
! similar to the approach above used for AM05. AJL
   if (flag_xc.ge.25) then
      i_exch = 8
   endif

! Caveat to catch all LibXC
! We'll solve calculate the free atom for them all as PW91-LD
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

!  prepare input for sr_atom()
   nuclear_charge = species_z(i_species)
   if ((i_exch.eq.3).or.(i_exch.eq.8)) then
     i_grad = 0
   else if ((i_exch.eq.6).or.(i_exch.eq.9)) then
     i_grad = 1
   end if
   if(flag_KH_core_states)then
      i_rel = 10
   elseif (flag_rel.eq.0) then
     i_rel=2
   else if (flag_rel.eq.REL_zora.or.flag_rel.eq.REL_atomic_zora) then
     i_rel=9
   else if (flag_rel.eq.REL_own.or.flag_rel==REL_KOLNING_HARMON) then
     i_rel = 10
   end if

!  input parameter for sr_atom()
!  FIXME not configured in control.in
   sc_max_iter = 100
   t_core = .true.
   n_core = 0
   t_in_pot = .false.
   t_sic = .false.
   t_kli = .false.
   t_mgga = .false.

   ! initialize handling of free-ion cutoff - see explanatory comment above ...
   cut_free_ion = .false.
   ion_r_cut = 0.d0

!  initialize auxiliary arrays for sratom()

   do i_grid = 1, n_max_grid, 1
     do i_channel = 1, n_channels,1
       density(i_grid, i_channel) = 0.d0
       density_deriv(i_grid, i_channel) = 0.d0
       density_2nd_deriv(i_grid, i_channel) = 0.d0
       core_density(i_grid,i_channel) = 0.d0
       core_dd(i_grid,i_channel) = 0.d0
       core_d2d(i_grid,i_channel) = 0.d0
     enddo
   enddo

!  index occupied shells
   i_fn = 0
   core_type = .false.
   do i_l = 0, l_shell_max(i_species), 1
   do i_shell = i_l+1, ion_n_max(i_l,i_species)
     if (ion_occ(i_shell, i_l, i_species).gt.0.) then
       if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then ! fully-relativistic case
          if(i_l.gt.0) then
             i_fn = i_fn+1
             free_ion_n(i_fn) = i_shell
             free_ion_l(i_fn) = i_l
             free_ion_k(i_fn) = i_l 
             i_fn= i_fn+1
             free_ion_n(i_fn) = i_shell
             free_ion_l(i_fn) = i_l
             free_ion_k(i_fn) = -(i_l+1)
          else
             i_fn = i_fn+1
             free_ion_n(i_fn) = i_shell
             free_ion_l(i_fn) = i_l
             free_ion_k(i_fn) = -(i_l+1)
          end if
       else ! non-rel
             i_fn = i_fn+1
             free_ion_n(i_fn) = i_shell
             free_ion_l(i_fn) = i_l
             m_shell(i_fn) = 1
             occ_shell(i_fn) = ion_occ(i_shell, i_l, i_species)
       endif
       if (i_shell.le.core_n_max(i_l,i_species)) then
         core_type(i_fn) = .true.
       end if
     end if
   enddo
   enddo
   n_free_ionic = i_fn

   ! Generate the occupation table for relativistic cases, where the orbitals are split.
   if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
     do i_fn = 1, n_free_ionic
       if (free_ion_l(i_fn).gt.0) then
           if (ion_occ(free_ion_n(i_fn),free_ion_l(i_fn),i_species).lt.1.d-8) then
               occ_shell(i_fn) = ion_occ(free_ion_n(i_fn),free_ion_l(i_fn),i_species)
           ! Check whether the shell is closed
           elseif ( dabs( ion_occ(free_ion_n(i_fn),free_ion_l(i_fn),i_species) - 2*(2*free_ion_l(i_fn)+1) ) &
               .lt. 1.d-8 ) then
               if(free_ion_k(i_fn).gt.0) then
                   occ_shell(i_fn) = 2*free_ion_l(i_fn)
               else
                   occ_shell(i_fn) = 2*(free_ion_l(i_fn)+1)
               end if
           ! If the shell is open, divide electrons correctly to two split spin orbitals. 
           ! One scheme is: all electrons should firstly go to the larger-kappa orbital,
           ! then the smaller-kappa orbital. Here we adopt the most popular average occupation
           ! scheme to keep the atomic electron density spherical.
           else
               if(free_ion_k(i_fn).gt.0) then
                   occ_shell(i_fn) = dble(2*free_ion_l(i_fn))/dble(2*(2*free_ion_l(i_fn)+1)) &
                                     * ion_occ(free_ion_n(i_fn),free_ion_l(i_fn),i_species)
               else
                   occ_shell(i_fn) = dble(2*(free_ion_l(i_fn)+1))/dble(2*(2*free_ion_l(i_fn)+1)) &
                                     * ion_occ(free_ion_n(i_fn),free_ion_l(i_fn),i_species)
               end if
           end if 
       else ! (free_ion_l(i_fn).le.0), for s orbital, there is no splitting due to spin
           occ_shell(i_fn) = ion_occ(free_ion_n(i_fn),free_ion_l(i_fn),i_species)
       end if

      !write(use_unit,"('i_fn:',i3,5x,'n:',i3,5x,'l:',i3,5x,'k:',i3,5x,'occu:',f7.4)"),i_fn,free_ion_n(i_fn),free_ion_l(i_fn),free_ion_k(i_fn),occ_shell(i_fn)
     enddo
  !else
  !  do i_fn = 1, n_free_ionic
  !    write(use_unit,"('i_fn:',i3,5x,'n:',i3,5x,'l:',i3,5x,'occu:',f7.4)"),i_fn,free_ion_n(i_fn),free_ion_l(i_fn),occ_shell(i_fn)
  !  enddo
   end if

!  initialize output
   do i_grid = 1, n_max_grid, 1
     free_ion_pot(i_grid)   = 0.0d0
     do i_fn = 1, n_free_ionic, 1
       free_ion_wave (i_grid, i_fn) = 0.0d0
       free_ion_kinetic (i_grid, i_fn) = 0.0d0
     enddo
   enddo

!  initialize atomic density by Thomas-Fermi model
   call atomini (n_free_ionic, n_grid(i_species), nuclear_charge, &
    free_ion_n, occ_shell, r_grid(1,i_species), ini_pot, free_ion_eigenval )

!  FIXME: This would be a superfluous copy, if atom_pot were not grossly overdimensioned
!  initialize also other potentials
   do i_grid = 1, n_grid(i_species), 1
     atom_pot(i_grid, 1) = ini_pot(i_grid)
     atom_pot_es(i_grid) = 0.0d0
     atom_pot_xc(i_grid, 1) = 0.0d0
   enddo

  ! Here we run the chosen atomic solver to generate candidate basis functions
  ! Should you wish to add in another solver, the values extracted from these 
  ! solvers are:
  !      eigenval, wave, atom_pot, atom_pot_es, density, density_deriv, density_2nd_deriv
  ! defined on the radial grid specified by r_grid (note that the l_channel
  ! index in these variables is legacy code that is not needed)
   if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then 
     ! relativistic case:
       call rel_atoms(i_species,nuclear_charge,n_channels,n_free_ionic,free_ion_n,free_ion_k,free_ion_l,&
       occ_shell,n_grid(i_species),r_grid(1,i_species),cut_free_ion,cutoff_type(i_species),ion_r_cut,&
       scale_cutoff(i_species),w_cutoff(i_species),i_grad,i_exch, &
!      outputs  
       atom_pot, atom_pot_es, density, density_deriv, density_2nd_deriv, density_diff, &
       free_ion_wave, free_ion_wave_deriv, small_wave, small_deriv, free_ion_eigenval)      
       
   elseif (atomic_solver .eq. ATOMIC_SOLVER_SRATOM) then
!    obtain self-consistent ion
     call sratom &
     ( sc_max_iter, i_rel, i_grad, i_exch, t_core, nuclear_charge, &
       n_core, n_free_ionic, free_ion_n, free_ion_l, m_shell, &
       occ_shell, free_ion_eigenval, i_grid_max, d_eig_sum, &
       n_max_grid, n_max_ind_fns, n_channels, &
       n_grid(i_species), r_grid(1,i_species), &
       cut_free_ion, cutoff_type(i_species), ion_r_cut, &
       scale_cutoff(i_species), w_cutoff(i_species), &
       r_peak, atom_pot, &
       atom_pot_es, atom_pot_xc, density, density_deriv, &
       density_2nd_deriv, core_density, core_dd, core_d2d, &
       free_ion_wave, t_in_pot, t_sic, e_sic, v_orb, t_kli, sv_kli, &
       svm, t_mgga, free_ion_wave_deriv, core_type &
     )
   else if (atomic_solver .eq. ATOMIC_SOLVER_ATOM_SPHERE) then
     ! atom_sphere_wrapper's argument list has been deliberately constructed to be an example of a generic 
     ! interface between aims and a radial atomic solver
     ! No modules are used within this wrapper, other than for system-level tasks like aims_stop.
     ! What you see here is an exhaustive list of what properties of the atom are needed from aims, what variables
     ! need to be exported back to aims to construct the basis sets, and nothing more (e.g. all work matrices reside 
     ! within the wrapper.) 

     ! Since the ionic potential must be meaningful (as it is used to generate a new set of basis functions which 
     ! may include unoccupied quantum numbers not found in the set of occupied states output by get_free_ion,) 
     ! always use the semi-local functional that sratom would have (i.e., pass in i_exch instead of flag_xc), so 
     ! that atom_pot is meaningful and can be integrated later
     !
     ! Idle thought:  is it possible to have atom_sphere calculate the unoccupied states needed by specifying their
     ! quantum numbers and assigning them an occupation of zero?  Do we even need it?
     call atom_sphere_wrapper &
       ( n_grid(i_species), r_grid(1,i_species), species_name(i_species), i_exch, hybrid_coeff, hse_omega, &
         nuclear_charge, n_core, n_free_ionic, n_max_ind_fns, free_ion_n, free_ion_l, m_shell, occ_shell, n_max_grid, & 
         n_channels, 1, &
         cut_free_ion, cutoff_type(i_species), ion_r_cut, scale_cutoff(i_species), &
              w_cutoff(i_species), & ! Cutoff potential
         free_ion_eigenval, free_ion_wave, free_ion_wave_deriv, density, density_deriv, density_2nd_deriv, &
         atom_pot, atom_pot_es, free_ion_kinetic, psi_at_nucleus_ionic(:n_valence,i_species)) ! Output
     call sync_atomic_solver_outputs( n_max_ind_fns, n_max_grid, n_channels, 1, &
          free_ion_eigenval, free_ion_wave, free_ion_wave_deriv, free_ion_kinetic, &
          density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es )
   else
     call aims_stop('The atomic solver specified by the variable atomic_solver was not recognized, exiting.','func')  
   end if

!test
!      write(use_unit,*) "After sratom."
!test end

!  now copy all data from atom onto storage arrays
   do i_grid = 1, n_grid(i_species), 1
     free_ion_pot(i_grid) = atom_pot(i_grid,1)
   enddo

!  output free ion data
   call localorb_info('',use_unit)

   write(info_str,'(2X,A)') "List of free ionic orbitals and eigenvalues: "
   call localorb_info(info_str,use_unit,'(A)')

   if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) then
     write(info_str,'(4X,A,4X,A,4X,A,6X,A,4X,A)') "n", "l", "k", "energy [Ha]", "energy [eV]"
     call localorb_info(info_str,use_unit,'(A)')

     do i_fn = 1, n_free_ionic, 1

       write(info_str,'(2X,I3,2X,I3,2X,I3,2X,F15.6,F15.4)') free_ion_n(i_fn), free_ion_l(i_fn), free_ion_k(i_fn), &
            free_ion_eigenval(i_fn), free_ion_eigenval(i_fn)*hartree
       call localorb_info(info_str,use_unit,'(A)')

     enddo
     call localorb_info('',use_unit)
   else
     write(info_str,'(4X,A,4X,A,6X,A,4X,A)') "n", "l", "energy [Ha]", "energy [eV]"
     call localorb_info(info_str,use_unit,'(A)')

     do i_fn = 1, n_free_ionic, 1

       write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4)') free_ion_n(i_fn), free_ion_l(i_fn), &
            free_ion_eigenval(i_fn), free_ion_eigenval(i_fn)*hartree
       call localorb_info(info_str,use_unit,'(A)')

     enddo
     call localorb_info('',use_unit)
   endif

!  that's all folks

   return
 end subroutine get_free_ion

