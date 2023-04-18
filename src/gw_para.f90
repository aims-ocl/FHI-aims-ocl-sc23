!****h* FHI-aims/gw_para
!  NAME
!    gw_para
!  SYNOPSIS

      module gw_para

!  PURPOSE
!  Defining the parameters used in GW calculations

!  USES

      use dimensions
      use contour_def_gw_types, only: cd_environment_type
      implicit none

!  INPUT
!  o none
!  OUTPUT
!  o none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!
!  SOURCE

!  n_homo   :  HOMO level
!  n_freq   :  number of frequency point for self energy
!  n_full_freq : number of frequency for the screened coulomb interaction,
!                typically n_full_freq = 2*n_freq
!  omegamax    : maximum frequency for the self energy
!  omega_grid(n_freq) : GL frequency grid for self energy
!  taumax      : maximum time for the self energy
!  womega(n_freq) :  weight of the frequency grid
!  omega_full_grid(n_full_freq) : GL frequency grid for the screened
!              coulomb interaction ranging from 0 to infinity
!  womega_full : weight for the full frequency grid
!
!  n_max_par :   number of parameters for the analytical continuation
!  wfitrange :   the frequency range in the analytical continuation
!  n_frozen_core : if (n_frozen_core .eq. 0) no frozen core shell
!                  if (n_frozen_core .gt. 0) n_frozen_core represents the number 
!                of the valence shells under which the shells are frozen
!                  if (n_frozen_core .lt. 0) n_frozen_core represents the number 
!                of the frozen core shells starting from the lowest shell (n=0) 
!  n_low_state : lowest KS state of the self energy will be calculated
!  n_high_state : highest KS state of the self energy will be calculated
!  max_n_prodbas(n_species) :  maximum principal shell used in the product
!                        basis construction
!  max_l_prodbas(n_species) :  maximum angular momentum of the product basis
!  anacon_type     : integer, anacon_type = 0, two-pole fitting,
!                                   = 1, Pade approximation
!  anacon_type_defined : logical, only true if anacon_type was specified in control.in
!  out_self_energy      : logical, if true, print out the self energy
!  out_gw_regular_kgrid : logical, if true, evaluate the self energy for a periodic system 
!                         on a uniform regular k grid
!  flag_coh     : logical, if true, use the coulomb hole accelerator
!  flag_qpe     : logical, if true, calculate quasi-particle energy in this run
!  flag_ene     : logical, if true, calculate RPA correlation energy
!  weighted_continuation  : logical, if true, do weighted analytical continuation

!  This is the default parameter setting for two-pole fitting, but now we changed to
!  the default to Pade approximation (anacon_type=1)
!      integer ::  n_freq = 40 
!      integer ::  n_full_freq = 80
!      integer ::  n_full_time = 80
!      integer ::  n_max_par =5
!      integer ::  n_low_state = 1
!      integer ::  n_high_state = 0
!      integer ::  state_to_print = 1
!      integer ::  anacon_type = 0
!      integer ::  n_iteration = 1 

! Default setting for Pade approximation
      integer ::  n_freq = 100
      integer ::  n_full_freq = 100
      integer ::  n_full_time = 100
      integer ::  n_max_par =16
      integer ::  n_low_state = 1
      integer ::  n_high_state = 0
      integer ::  state_to_print = 1
      integer ::  anacon_type = 1
      integer ::  n_iteration = 1
      integer ::  state_to_shift = -100
      logical :: anacon_type_defined = .false.
! Using the modified Gauss-Legendre grid as the default
! Thanks to Alexandre Tkatchenko for suggesting this. -- X. Ren
      integer ::  freq_grid_type = 1

      real(kind=8), dimension(2)  ::  print_range 
      real*8  ::  omegamax = 10.d0
      real*8  ::  taumax   = 7.d0

      logical :: out_self_energy  = .false.
      logical :: out_gw_regular_kgrid  = .false.
      logical :: flag_coh = .false.
      logical :: flag_qpe = .false.
      logical :: qpe_multi_solu = .false.
      logical :: gw_zshot = .false.
      logical :: use_hedin_shift = .false.

      real*8, allocatable ::   omega_grid(:)
      real*8, allocatable ::   omega_full_grid(:)
      real*8, allocatable ::   womega(:)
      real*8, allocatable ::   womega_full(:)
      real*8, allocatable ::   hedin_shift(:,:)
      real*8, allocatable ::   z_value(:,:)
      complex*16, allocatable ::   sigma_par(:,:,:)
      complex*16, allocatable ::   sigma_par_p0(:,:,:,:)
      complex*16, allocatable ::   dielec_func_imagfreq(:)

      type(cd_environment_type) :: gw_cd

      contains

         subroutine allocate_gw ()
! allocate arrays

         if(freq_grid_type.eq.0) then
            n_freq=n_full_freq/2
         elseif(freq_grid_type.eq.1) then
!            n_freq=n_full_freq/2
!            n_freq=40
            n_freq=n_full_freq
         elseif(freq_grid_type.eq.2) then
            n_freq=n_full_freq
         endif 

         if(.not.allocated(omega_grid)) then
           allocate (omega_grid(n_freq))
         endif
         if(.not.allocated(womega)) then
           allocate (womega(n_freq))
         endif
         if(.not.allocated(omega_full_grid)) then
           allocate (omega_full_grid(n_full_freq))
         endif
         if(.not.allocated(womega_full)) then
           allocate (womega_full(n_full_freq))
         endif
         if(.not.allocated(hedin_shift)) then
           allocate (hedin_shift(n_states,n_spin))
           hedin_shift(:,:) = 0.0d0
         endif
         if(.not.allocated(z_value).and.gw_zshot) then
           allocate (z_value(n_states,n_spin))
           z_value(:,:) = 0.0d0
         endif

         end subroutine allocate_gw

! deallocate things
         subroutine deallocate_gw ()

         if(allocated (omega_grid)) then
           deallocate (omega_grid)
         endif
         if(allocated (womega)) then
           deallocate (womega)
         endif
         if(allocated (omega_full_grid)) then
           deallocate (omega_full_grid)
         endif
         if(allocated (womega_full)) then
           deallocate (womega_full)
         endif
         if(allocated(hedin_shift)) then
           deallocate(hedin_shift)
         endif
         if(allocated(z_value)) then
           deallocate(z_value)
         endif

         end subroutine deallocate_gw

! **************************************************************************************************
!> brief set defaults for contour deformation
!  o cd --contour deformation environment
! **************************************************************************************************
    subroutine set_defaults_cd(cd)

      type(cd_environment_type)                          :: cd

      integer                                            :: i_spin   
 
      cd%iterative = .true.
      cd%restart = 'none'
      cd%full_cmplx_sigma = .false.
      cd%try_zshot = .true.
      cd%calc_single_state_spec = .false.

      cd%eta = 0.001d0
      cd%contour_def_offset = 0.002d0 
      cd%n_iter_sc = 20
      
      if(.not.allocated(cd%contour_def_start)) then
        allocate(cd%contour_def_start(n_spin))
      endif
      if(.not.allocated(cd%contour_def_end)) then
        allocate(cd%contour_def_end(n_spin))
      endif
      if(.not.allocated(cd%spin_channel)) then
        allocate(cd%spin_channel(n_spin))
      endif

      if(.not.allocated(cd%sc_env)) then
        allocate(cd%sc_env(n_spin))
      endif

      do i_spin=1,n_spin
        cd%sc_env(i_spin)%n_occ  = 5
        cd%sc_env(i_spin)%n_virt = 5
        cd%sc_env(i_spin)%reiterate = .false.
        cd%spin_channel(i_spin) = i_spin
      enddo
     
    end subroutine set_defaults_cd

! **************************************************************************************************
!> brief allocate QP_adjust environment
!  o cd --contour deformation environment
! **************************************************************************************************
    subroutine allocate_qp_adjust(cd)

      type(cd_environment_type)                          :: cd

      integer                                            :: n

      n = cd%QP_adjust%n_qps
      if(.not.allocated(cd%QP_adjust%state)) then
         allocate(cd%QP_adjust%state(n))
      endif
      if(.not.allocated(cd%QP_adjust%values)) then
         allocate(cd%QP_adjust%values(n))
      endif

    end subroutine allocate_qp_adjust


      end module gw_para
!******
