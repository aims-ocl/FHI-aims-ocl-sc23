!****h* FHI-aims/molecular_dynamics
!  NAME
!    molecular_dynamics - MD capacities of FHI-aims
!  SYNOPSIS

module molecular_dynamics

!  PURPOSE
!    This module takes care of everything related to MD and requires only energies and forces
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
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
  use localorb_io,only:use_unit
  implicit none

  ! thermostat internal variables
  real*8 :: s_NP
  real*8 :: s_NP_last
  real*8 :: s_NP_half
  real*8 :: s_dot_NP_half
  real*8 :: s_dot_NP_last
  real*8 :: tsystem
  real*8 :: tsystem_last
!  BDP thermostat internal variables
  !real*8 :: BDP_hh0
  !real*8 :: BDP_phase_space_factor
  real*8 :: BDP_conint
  real*8 :: BDP_consv ! still useful?
  real*8 :: BDP_factor
! here vars for NH chain GNH, pi_half, eta_half, eta_last, pi_last
  real*8, dimension(:,:), allocatable :: v_half, v_last, r_last!, v_comparison
  integer, dimension(:), allocatable :: counter
  integer*8, dimension(7), private :: RNG_a, RNG_c, RNG_m
  integer, private :: n_RNG
  logical :: use_MD_max_steps
  logical :: MD_RNG_firstcall
  real*8  :: MD_H0
  real*8  :: MD_Epot_last
  logical :: MD_successful_restart_read
  ! variables for 4th order integrator
  integer :: MD_max_high_order_index
  real*8, dimension(:), private, allocatable :: c_high_order_coeff, d_high_order_coeff
  integer :: MD_high_order_i
  integer :: MD_schedule_step
  integer :: n_atoms_MD
  logical, dimension(:), allocatable :: constrain_MD
! variables and matrices for GLE thermostat
  real*8, allocatable ::  gS(:,:), gT(:,:), gp(:,:,:), ngp(:,:,:), gA(:,:), gC(:,:)
  real*8  langham, langham_last, langham_half
  ! Internal RNG variables
  logical, private :: random_gauss_saved
  real*8, private :: random_gauss_r1, random_gauss_r2
  logical, private :: MB_sample_box_mueller_saved
  real*8, private :: MB_sample_box_mueller_r1, MB_sample_box_mueller_r2
  integer*8, private :: random_own_Y
  integer, private :: gasdev_iset
  real*8, private :: gasdev_gset
contains

!******
!------------------------------------------------------------------------------
!****s* molecular_dynamics/allocate_MD
!  NAME
!    allocate_MD
!  SYNOPSIS

subroutine allocate_MD

!  PURPOSE
!    allocation of MD
!  USES
  use dimensions
  use runtime_choices
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  if (.not.allocated(r_last)) then
     allocate(r_last(3,n_atoms))
     r_last(:,:) = 0d0
  end if
  if (.not.allocated(v_half)) then
     allocate(v_half(3,n_atoms))
     v_half(:,:) = 0d0
  end if
  if (.not.allocated(v_last)) then
     allocate(v_last(3,n_atoms))
     v_last(:,:) = 0d0
  end if
!  if (.not.allocated(v_comparison)) then
!     allocate(v_comparison(3,n_atoms))
!     v_comparison(:,:) = 0d0
!  end if
  if (.not.allocated(counter)) then
     allocate(counter(n_atoms))
     counter(:) = 0
  end if
  if (.not.allocated(constrain_MD)) then
     allocate(constrain_MD(n_atoms))
     constrain_MD(:) = .false.
  end if
! allocate the various settings for MD schedule
  if (MD_use_schedule) then
     if (.not.allocated(MD_schedule_ensemble      )) allocate(MD_schedule_ensemble      (MD_segments))
     if (.not.allocated(MD_schedule_temperature   )) allocate(MD_schedule_temperature   (MD_segments))
     if (.not.allocated(MD_schedule_time          )) allocate(MD_schedule_time          (MD_segments))
     if (.not.allocated(MD_schedule_tau_berendsen )) allocate(MD_schedule_tau_berendsen (MD_segments))
     if (.not.allocated(MD_schedule_Q             )) allocate(MD_schedule_Q             (MD_segments))
     if (.not.allocated(MD_schedule_damping_factor)) allocate(MD_schedule_damping_factor(MD_segments))
     if (.not.allocated(MD_schedule_nu_andersen   )) allocate(MD_schedule_nu_andersen   (MD_segments))
     if (.not.allocated(MD_schedule_tau_BDP       )) allocate(MD_schedule_tau_BDP       (MD_segments))
     if (.not.allocated(MD_schedule_random_BDP    )) allocate(MD_schedule_random_BDP    (MD_segments))
     ! initialize to some neutral value
     MD_schedule_ensemble      (:) = "NVE"
     MD_schedule_temperature   (:) = 0d0
     MD_schedule_time          (:) = 0d0
     MD_schedule_tau_berendsen (:) = 1d0
     MD_schedule_Q             (:) = 1d0
     MD_schedule_damping_factor(:) = 1d0
     MD_schedule_nu_andersen   (:) = 1d0
     MD_schedule_tau_BDP       (:) = 1d0
     MD_schedule_random_BDP    (:) = .true.
  end if
  ! set degrees of freedom for thermostats
  MD_g_DOF = 3*n_atoms

end subroutine allocate_MD
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/clean_MD
!  NAME
!    clean_MD - deallocation of molecular dynamics
!  SYNOPSIS
subroutine clean_MD
!  PURPOSE
!    deallocation of local variables
!  USES
   use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(r_last))         deallocate(r_last)
  if (allocated(v_half))         deallocate(v_half)
  if (allocated(v_last))         deallocate(v_last)
  if (allocated(counter))        deallocate(counter)
  if (allocated(MD_schedule_ensemble      )) deallocate(MD_schedule_ensemble      )
  if (allocated(MD_schedule_temperature   )) deallocate(MD_schedule_temperature   )
  if (allocated(MD_schedule_time          )) deallocate(MD_schedule_time          )
  if (allocated(MD_schedule_tau_berendsen )) deallocate(MD_schedule_tau_berendsen )
  if (allocated(MD_schedule_Q             )) deallocate(MD_schedule_Q             )
  if (allocated(MD_schedule_damping_factor)) deallocate(MD_schedule_damping_factor)
  if (allocated(MD_schedule_nu_andersen   )) deallocate(MD_schedule_nu_andersen   )
  if (allocated(MD_schedule_tau_BDP       )) deallocate(MD_schedule_tau_BDP       )
  if (allocated(MD_schedule_random_BDP    )) deallocate(MD_schedule_random_BDP    )
  if (allocated(c_high_order_coeff        )) deallocate(c_high_order_coeff        )
  if (allocated(d_high_order_coeff        )) deallocate(d_high_order_coeff        )
  if (allocated(constrain_MD              )) deallocate(constrain_MD              )
  if (MD_ensemble .eq. 'GLE_thermostat') then
      call clean_gle()
  end if
end subroutine clean_MD
!******
subroutine change_MD_schedule_step
  use runtime_choices
  use dimensions
  use species_data
  use geometry
  use localorb_io
  implicit none
  character*120 :: info_str
  integer :: i_atom

  ! determine if one of the steps in the molecular dynamics schedule is over and we need to switch to the next one ... 
  if ((tsystem.ge.MD_time).and.(MD_use_schedule).and.(MD_schedule_step.lt.MD_segments)) then
     ! Nose thermostats propagate momenta, others propagate velocities.
     ! make sure we change only velocities, i.e. take out mass factos if necessary
     if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
     	do i_atom = 1, n_atoms
        	v_half(:,i_atom) = v_half(:,i_atom)/species_m(species(i_atom))
		v_last(:,i_atom) = v_last(:,i_atom)/species_m(species(i_atom))
     	end do
     end if	


     ! set all MD variables to the current step
     MD_schedule_step   = MD_schedule_step + 1 
     MD_ensemble        = MD_schedule_ensemble      (MD_schedule_step)
     ! CC: Ensure that a consistent number of steps is perfomed in each segment
     !     MD_time            = MD_schedule_time          (MD_schedule_step) + MD_time 
     MD_time            = MD_schedule_time          (MD_schedule_step) + MD_time + MD_tstep*0.99999999d0
     MD_temperature     = MD_schedule_temperature   (MD_schedule_step)
     MD_tau_berendsen   = MD_schedule_tau_berendsen (MD_schedule_step)
     MD_Q_NP            = MD_schedule_Q             (MD_schedule_step)
     NVE_damping_factor = MD_schedule_damping_factor(MD_schedule_step)
     MD_nu_andersen     = MD_schedule_nu_andersen   (MD_schedule_step)
     MD_tau_BDP         = MD_schedule_tau_BDP       (MD_schedule_step)
     !MD_random_BDP      = MD_schedule_random_BDP    (MD_schedule_step)


     call localorb_info(' ',use_unit,'(A)')
     write(info_str,'(2X,A,I3,A)') 'Molecular dynamics: Reached end of step ',MD_schedule_step-1,' in molecular dynamics schedule' 
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A,I3)') '| Starting schedule step number : ', MD_schedule_step
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,2A)')   '| new ensemble                  : ', MD_ensemble 
     call localorb_info(info_str,use_unit,'(A)')
     if ((MD_ensemble.eq.'NVE').or.(MD_ensemble.eq.'NVE_4th_order')) then
     else if (MD_ensemble.eq.'NVE_damped') then
        write(info_str,'(2X,A,E14.6)')   '| damping factor                : ', NVE_damping_factor
        call localorb_info(info_str,use_unit,'(A)')        
     else if (MD_ensemble.eq.'NVT_berendsen') then
        write(info_str,'(2X,A,E14.6)')   '| molecular dynamics temperature: ', MD_temperature
        call localorb_info(info_str,use_unit,'(A)')        
        write(info_str,'(2X,A,E14.6)')   '| Berendsen relaxation time     : ', MD_tau_berendsen
        call localorb_info(info_str,use_unit,'(A)')        
     else if (MD_ensemble.eq.'NVT_andersen') then
        write(info_str,'(2X,A,E14.6)')   '| molecular dynamics temperature: ', MD_temperature
        call localorb_info(info_str,use_unit,'(A)')        
        write(info_str,'(2X,A,E14.6)')   '| Andersen collision frequency  : ', MD_nu_andersen
        call localorb_info(info_str,use_unit,'(A)')        
     else if (MD_ensemble.eq.'NVT_parrinello') then
        write(info_str,'(2X,A,E14.6)')   '| molecular dynamics temperature      : ', MD_temperature
        call localorb_info(info_str,use_unit,'(A)')        
        write(info_str,'(2X,A,E14.6)')   '| Bussi-Donadio-Parrinello relax. time: ', MD_tau_BDP
        call localorb_info(info_str,use_unit,'(A)')        
     else if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
        write(info_str,'(2X,A,E14.6)')   '| molecular dynamics temperature: ', MD_temperature
        call localorb_info(info_str,use_unit,'(A)')        
        write(info_str,'(2X,A,E14.6)')   '| thermostat mass               : ', MD_Q_NP
        call localorb_info(info_str,use_unit,'(A)')        
     end if
     call localorb_info(' ',use_unit,'(A)')

     !CC: Thermodynamic integration variables, if used:
     if(use_thermodynamic_integration .or. use_reffree_AS) then
       call TDI_change_schedule_step(MD_schedule_step)
     end if

     ! again, if we are in a nose-thermostat NOW, need to throw in masses to be able to propagate momenta in the future. 
     if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
     	do i_atom = 1, n_atoms
        	v_half(:,i_atom) = species_m(species(i_atom))*v_half(:,i_atom)
		v_last(:,i_atom) = species_m(species(i_atom))*v_last(:,i_atom)
     	end do
     end if	

    
  end if
end subroutine change_MD_schedule_step
!****s* molecular_dynamics/write_MD_restart_binary
!  NAME
!    write_MD_restart_binary
!  SYNOPSIS
subroutine write_MD_restart_binary ![OBSOLETE]
!  PURPOSE
!    keep restart information (velocity, forces, thermostat parameters & such)
!    at the end of every time step; potentially this can be extended to the end of
!    every nth time step if really necessary 
!    [OBSOLETE] A new write_MD_restart function is provided that uses a human
!    readable, cross-architecture-compatible format
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use timing
  use species_data
  use localorb_io
  implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
! Here is the format for these restart files: 
!   (1) integer :: number of atoms   --  (merely for consistency checking
!   (2) real*8, dimension(3,n_atoms) :: current positions
!   (3) real*8, dimension(3,n_atoms) :: v_half
!   (4) real*8                       :: s
!   (5) real*8                       :: phalf
! this current format is the minimal thing that would work for all currently
! implemented thermostats while it still allows switching between various
! thermostats in different runs (or changing the external parameters for 
! example). If the need for any other numbers arises, please add them to THE
! END of this file and initialize s=1; phalf=0 unless you absolutely need these
! for your integrator.
  real*8 :: s_write, s_last_write
  real*8, dimension(3,n_atoms) :: v_half_write, v_last_write
  integer i_atom  
  if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
     do i_atom = 1, n_atoms
        v_half_write(:,i_atom) = v_half(:,i_atom)/species_m(species(i_atom))
	v_last_write(:,i_atom) = v_last(:,i_atom)/species_m(species(i_atom))
     end do
  else
     do i_atom = 1, n_atoms
        v_half_write(:,i_atom) = v_half(:,i_atom)
	v_last_write(:,i_atom) = v_last(:,i_atom)
     end do
  end if	!this does not work if inside myid .eq.0
 ! only write from the zero task ... 
  if (myid.eq.0) then
     s_write      = s_NP
     s_last_write = s_NP_last
     if (MD_ensemble.eq.'NVT_nose-hoover') then
        s_write      = dexp(s_write)
        s_last_write = dexp(s_last_write)
     end if

     open(file = MD_restart_file, unit = 88, status = 'unknown', form = 'unformatted', action='write')
     write(88) n_atoms               
     write(88) coords                
     write(88) v_half_write          
     write(88) v_last_write          
     write(88) r_last                
     write(88) s_write               
     write(88) s_NP_half             
     write(88) s_last_write          
     write(88) s_dot_NP_half         
     write(88) s_dot_NP_last         
     write(88) tsystem               
     write(88) tsystem_last          
     write(88) MD_H0                 
     write(88) MD_Epot_last          
     write(88) BDP_conint          
     write(88) MD_stepcount          
     write(88) MD_high_order_i       
     write(88) MD_force_evaluations  
     close(unit=88)
  end if
end subroutine write_MD_restart_binary
!******

!****s* molecular_dynamics/initialize_MD
!  NAME
!    initialize_MD
!  SYNOPSIS
subroutine initialize_MD(Epot)
!  PURPOSE
!    initialization routine and first MD step
!  USES
  use runtime_choices
  use physics
  use geometry
  use constants
  use species_data
  use dimensions
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use timing
  use relaxation
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o Epot - potential energy of current configuration
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none

  integer :: i_atom, i_coord
  real*8  :: T_inst, Epot, Ekin, shat, factor, w0, w1, w2, w3
  real*8, dimension(3,n_atoms) :: MD_forces
  real*8, dimension(n_atoms)   :: masses
  real*8, dimension(n_atoms)   :: charges 
  real*8, dimension(n_atoms)   :: px,py,pz 
  real*8, dimension(n_atoms)   :: fx,fy,fz 
  real*8 :: plumed_energy,plumed_eunit 
  integer :: cell_type
  logical :: read_restart_succes
  character*100 :: info_str

  call initialize_RNG ()

!===============================
! higher-order integrators with magic coefficients. For the time being, simply 
! uncomment the integrator you would like to use until we settle on one particular integrator 
! for each task in question, THEN clean up the code ... 
!
!!$  ! 4th order integrator according to Eqn 2.11 of Yoshida, Phys. Lett. A, v150, p262 (1990)
!!$  ! These magically defined coefficients give the proper operator splitting, see original paper for details. 
!!$  ! need to specify the maximal number of coefficients (i.e. momentum updates) explicitly
!!$  MD_max_high_order_index =  4
!!$  if (.not.allocated(c_high_order_coeff)) allocate(c_high_order_coeff(MD_max_high_order_index))
!!$  if (.not.allocated(d_high_order_coeff)) allocate(d_high_order_coeff(MD_max_high_order_index))
!!$  ! assign constant magic coefficients to 4th order symplectic integrator:
!!$  c_high_order_coeff(:) = 0d0
!!$  d_high_order_coeff(:) = 0d0
!!$  ! coefficients for momentum update - four updates per time step 
!!$  c_high_order_coeff(1)   =  0.67560359597982881702d0  ! = 1/[2(2-2^(1/3))]
!!$  c_high_order_coeff(2)   = -0.17560359597982881702d0  ! = (1-2^(1/3))*c1
!!$  c_high_order_coeff(3)   = -0.17560359597982881702d0  ! = c2
!!$  c_high_order_coeff(4)   =  0.67560359597982881702d0  ! = c1
!!$  ! coefficients for position update - three force calculations per time step
!!$  d_high_order_coeff(1)   =  1.35120719195965763405d0  ! = 2 c1
!!$  d_high_order_coeff(2)   = -1.70241438391931526810d0  ! = -2^(1/3) d1
!!$  d_high_order_coeff(3)   =  1.35120719195965763405d0  ! = d1
!
!--------------------------------------
!
  s_NP = 1.0d0
  s_NP_last = 1.0d0
  s_NP_half = 1.0d0
  s_dot_NP_half = 0.0d0
  s_dot_NP_last = 0.0d0
  tsystem = 0.0d0
  tsystem_last = 0.0d0
  BDP_conint = 0.0d0
  BDP_consv = 0.0d0
  BDP_factor = 1.0d0
  n_RNG = 7
  use_MD_max_steps = .false.
  MD_RNG_firstcall = .true.
  MD_H0 = 0.0d0
  MD_Epot_last = 0.0d0
  MD_successful_restart_read = .false.
  MD_high_order_i = 0

  ! the 4th order version used by Ishida et al. CPL 282, v115, integrator SI4
  ! ATTENTION: this integrator has one extra force evaluation per time step compared to the 
  ! fourth order integrator borrowed from the Yoshida paper. 
  MD_max_high_order_index  = 5
  if (.not.allocated(c_high_order_coeff)) allocate(c_high_order_coeff(MD_max_high_order_index))
  if (.not.allocated(d_high_order_coeff)) allocate(d_high_order_coeff(MD_max_high_order_index))
  ! assign constant magic coefficients to 4th order symplectic integrator:
  c_high_order_coeff(:) = 0d0
  d_high_order_coeff(:) = 0d0
  ! coefficients for momentum update
  c_high_order_coeff(1) =  0.061758858135626325d0
  c_high_order_coeff(2) =  0.338978026553643355d0
  c_high_order_coeff(3) =  0.614791307175577566d0
  c_high_order_coeff(4) = -0.140548014659373380d0
  c_high_order_coeff(5) =  0.125019822794526133d0
  ! coefficients for position update 
  d_high_order_coeff(1) =  0.205177661542286386d0
  d_high_order_coeff(2) =  0.403021281604214587d0
  d_high_order_coeff(3) = -0.120920876338914008d0
  d_high_order_coeff(4) =  0.512721933192413035d0
!
!----------------------------------------
!
!!$  ! sixth order integrator, proper Yoshida splitting, 
!!$  ! according to Table 1 (and below) of Yoshida Phys. Lett. A, v150, p262 (1990)
!!$  ! there are three different possible solutions to the equations making up the integrator ...
!!$  ! NOTE: Yoshida also describes (and gives solutions for) an eigth-order integrator, which can 
!!$  ! be set up by simply copying down the solutions from his paper, page 267 and extending the 
!!$  ! equations for w0 as well as adding a few coefficients to the established pattern. 
!!$  MD_max_high_order_index  = 8
!!$  if (.not.allocated(c_high_order_coeff)) allocate(c_high_order_coeff(MD_max_high_order_index))
!!$  if (.not.allocated(d_high_order_coeff)) allocate(d_high_order_coeff(MD_max_high_order_index))
!!$  ! solution A - this was the only "good" solution, the other two were much worse at the same computational cost
!!$  w1 = -0.117767998417887d1
!!$  w2 =  0.235573213359357d0
!!$  w3 =  0.784513610477560d0
!!$  ! this integrator is obtained with Yoshida's m=3, thus giving the actual integration
!!$  ! coefficients with the right eqn on p.267:
!!$  w0 = 1d0-2d0*(w1+w2+w3)
!!$  c_high_order_coeff(:) = 0d0
!!$  d_high_order_coeff(:) = 0d0
!!$  c_high_order_coeff(1) = w3/2d0
!!$  c_high_order_coeff(2) = (w2+w3)/2d0
!!$  c_high_order_coeff(3) = (w1+w2)/2d0
!!$  c_high_order_coeff(4) = (w0+w1)/2d0
!!$  c_high_order_coeff(5) = (w0+w1)/2d0
!!$  c_high_order_coeff(6) = (w1+w2)/2d0
!!$  c_high_order_coeff(7) = (w2+w3)/2d0 
!!$  c_high_order_coeff(8) = w3/2d0
!!$  d_high_order_coeff(1) = w3
!!$  d_high_order_coeff(2) = w2
!!$  d_high_order_coeff(3) = w1
!!$  d_high_order_coeff(4) = w0
!!$  d_high_order_coeff(5) = w1
!!$  d_high_order_coeff(6) = w2
!!$  d_high_order_coeff(7) = w3
  
!===== End integrator listing ========
  ! Check for correct settings:
  if ((MB_clean_rotations).AND.(n_periodic .gt. 0)) then
     write(info_str,'(2X,A)') '*** WARNING: Rotational cleaning should not be used for periodic simulations.'
     call localorb_info(info_str)
  end if
  if ( (MD_ensemble.eq.'NVE').and.(.not.MB_clean_rotations).AND. &
       (n_periodic .eq. 0).AND.(.not.use_relaxation_constraints)) then
     write(info_str,'(2X,A)') '*** WARNING: Rotational cleaning should be used for non-periodic simulations.'
     call localorb_info(info_str)
  end if
  if ((MB_clean_rotations).AND.(use_relaxation_constraints)) then 
     write(info_str,'(2X,A)') 'WARNING: Rotational cleaning can not be used together with constrained atoms.'
     call localorb_info(info_str)
  end if

  ! initialize the force cleaning routines, akin to initialize_relaxation
  if (remove_unitary_force_components .eq. 2) then
     call initialize_cleaning_forces()
  end if
 ! initialize MD schedule if present: simply set all the necessary variables for the first step in the schedule
  if (MD_use_schedule) then
     ! FIXME: Intelligent output about the initialization!!
     MD_schedule_step   = 1 
     MD_ensemble        = MD_schedule_ensemble      (MD_schedule_step)
     MD_time            = MD_schedule_time          (MD_schedule_step) 
     MD_temperature     = MD_schedule_temperature   (MD_schedule_step)
     MD_tau_berendsen   = MD_schedule_tau_berendsen (MD_schedule_step)
     MD_Q_NP            = MD_schedule_Q             (MD_schedule_step)
     NVE_damping_factor = MD_schedule_damping_factor(MD_schedule_step)
     MD_nu_andersen     = MD_schedule_nu_andersen   (MD_schedule_step)
     MD_tau_BDP         = MD_schedule_tau_BDP       (MD_schedule_step)
     !MD_random_BDP      = MD_schedule_random_BDP    (MD_schedule_step)
  end if

  ! find out number of mobile atoms, if some are restricted 
  if (use_relaxation_constraints) then 
     n_atoms_MD = 0
     constrain_MD(:) = constrain_relaxation(:)
     do i_atom = 1, n_atoms
        if (.not.constrain_MD(i_atom)) then
           n_atoms_MD = n_atoms_MD + 1 
        end if        
     end do
     ! reset number of degrees of freedom
     MD_g_DOF = 3*n_atoms_MD
  else
     n_atoms_MD = n_atoms
  end if

  !------------------------------------------------------------------------------------------------------
  ! initialize velocities from MB distribution if requested
  if (MB_velocity_initialization) then
     write(info_str,'(2X,A)') 'Initializing velocities for molecular dynamics using Maxwell-Boltzmann distribution'
     call localorb_info(info_str)
     ! initialize on thread 0
     if (myid.eq.0) then
        if (.not.use_relaxation_constraints) then 
           do i_atom = 1, n_atoms
              do i_coord = 1, 3
                 call MB_sample_box_mueller(MD_init_temperature,species_m(species(i_atom)),v_half(i_coord,i_atom))
              end do
           end do
        else
           do i_atom = 1, n_atoms
              if (.not.constrain_MD(i_atom)) then
                 do i_coord = 1, 3
                    call MB_sample_box_mueller(MD_init_temperature,species_m(species(i_atom)),v_half(i_coord,i_atom))
                 end do
              end if
           end do
        end if
       ! set exact temperature
       ! first take out all net translations, also do that with constraints - but only if the user wants it.... (cluster case)
       ! CC: Take away translations in periodic case 
       if ( (MB_clean_rotations) .or. (n_periodic .gt. 0)) call clean_velocities(v_half,.true.)  
       call calculate_kinetic_energy(v_half,T_inst)           ! then calculate kinetic energy
       T_inst = 2d0*T_inst/(3d0*dble(n_atoms_MD)*boltzmann_kB)   ! turn this into a temperature
       factor = sqrt(MD_init_temperature/T_inst)              ! calculate velocity scaling factor
       v_half(:,:) = factor*v_half(:,:)                       ! scale all velocities
     end if
     ! distribute to all other threads
     call broadcast_MD_velocities(v_half,0)
  else
     ! CC: Take away translations in periodic case 
     if ( (MB_clean_rotations) .or. (n_periodic .gt. 0)) call clean_velocities(v_half)  
  end if
   ! Nosé Thermostats need to work with generalized momenta - must change from velocities!!!
  if (MD_ensemble.eq.'NVT_nose-poincare') then
     do i_atom = 1, n_atoms
        v_half(:,i_atom) = species_m(species(i_atom))*v_half(:,i_atom)
     end do     
  else if (MD_ensemble.eq.'NVT_nose-hoover') then
     s_NP_last = 0d0
     s_NP_half = 0d0
     s_NP      = 0d0
     do i_atom = 1, n_atoms
        v_half(:,i_atom) = species_m(species(i_atom))*v_half(:,i_atom)
     end do
  end if
  v_last(:,:)  = v_half(:,:)
  r_last(:,:)  = coords(:,:)
  MD_Epot_last = Epot
  MD_stepcount = 1   
  MD_high_order_i = 0
  MD_schedule_step = 1

   !initialize PLUMED
    if(plumed_plugin) then
	MD_forces(:,:) = total_forces(:,:)
	if(myid.eq.0) then
		do i_atom = 1, n_atoms
			masses(i_atom) = species_m(species(i_atom))
                        charges(i_atom) = 0.d0
                        px(i_atom)=coords(1,i_atom)
                        py(i_atom)=coords(2,i_atom)
                        pz(i_atom)=coords(3,i_atom)
                        fx(i_atom)=MD_forces(1,i_atom)
                        fy(i_atom)=MD_forces(2,i_atom)
                        fz(i_atom)=MD_forces(3,i_atom)
		end do
		!call init_metadyn_(n_atoms, MD_tstep, masses, trim(adjustl(plumed_file))//char(0))
                ! 1.3 interface
                if ( n_periodic .ge.1) then   
                    cell_type=1 
                else
                    cell_type=0 
                endif 
                plumed_eunit=1.d0 
                call init_metadyn_(n_atoms, MD_tstep, masses, charges, &
                     cell_type, plumed_eunit, &
                     trim(adjustl(plumed_file))//char(0))
                !write(use_unit,*) "fhi atoms MD_tstep plumed_file",n_atoms, MD_tstep,plumed_file
		!write(use_unit,*) "fhi masses", masses
   		!PLUMED force added
		!call meta_force_calculation_(n_periodic,lattice_vector,coords,MD_forces,MD_stepcount)
                ! 1.3 interface
                  
                call meta_force_calculation_( lattice_vector, MD_stepcount , px, py, pz , fx, fy, fz,  plumed_energy)

        	do i_atom = 1, n_atoms
                        MD_forces(1,i_atom)=fx(i_atom)
                        MD_forces(2,i_atom)=fy(i_atom)
                        MD_forces(3,i_atom)=fz(i_atom)
		end do
	        
                
		!write(use_unit,*) "fhi n_periodic, MD_stepcount", n_periodic, MD_stepcount
		!write(use_unit,*) "fhi coords", coords
   		!write(use_unit,*) "fhi forces", MD_forces
	endif
    	call broadcast_MD_velocities(MD_forces,0)
        MD_forces(:,:) = MD_forces(:,:)/MD_KE_factor
    else	
       ! change forces to MD integrator units of [amu bohr/ps^2]
       MD_forces(:,:) = total_forces(:,:)/MD_KE_factor
       ! delete those forces that should be zero from relaxation constraints: 
       if (use_relaxation_constraints) then 
          do i_atom = 1, n_atoms
             if (constrain_MD(i_atom)) MD_forces(:,i_atom) = 0d0
          end do
       end if
    endif
    ! end init PLUMED and first call
  
  !------------------------------------------------------------------------------------------------------
  ! initialize MD counters
  tsystem      = 0d0
  if (MD_maxsteps.gt.0) use_MD_max_steps = .true.
  
  !------------------------------------------------------------------------------------------------------


  ! first integration step:
  if (MD_ensemble.eq.'NVT_nose-poincare') then
!     MD_H0        = Epot + 1.5d0*dble(n_atoms)*boltzmann_kB*MD_temperature 
     call calculate_kinetic_energy_nose(v_last(:,:),1d0,Ekin)
     MD_H0        = Epot + Ekin !Benedikt: here we compute H0 with s=1 and sdot=0 
!     call initial_step_NP99_GLA(MD_tstep,MD_Q_NP,coords,v_half,s_NP,s_dot_NP_half,v_last, &
!          s_dot_NP_last,MD_forces,Epot,MD_H0,MD_temperature, tsystem_last)  
     call initial_step_NP01_GLA(MD_tstep,MD_Q_NP,coords,v_half,s_NP,s_dot_NP_half,v_last, &
          s_dot_NP_last,MD_forces)
     tsystem = tsystem + MD_tstep
  else if (MD_ensemble.eq.'NVT_nose-hoover') then
     call initial_step_NH96_GLA(MD_tstep,MD_Q_NP,coords,v_half,v_last, &
                 s_NP_half,s_dot_NP_half,shat,MD_forces,MD_temperature)
     tsystem = tsystem + MD_tstep
  else if (MD_ensemble.eq.'NVE_4th_order') then 
     call initial_fourth_order_NVE_step(MD_tstep,coords,v_half,MD_forces)
  else if ((MD_ensemble.eq.'NVE').or. & 
           (MD_ensemble.eq.'NVT_andersen').or.(MD_ensemble.eq.'NVT_berendsen').or. &
	    (MD_ensemble.eq.'NVT_damped')) then
     tsystem = MD_tstep
     !  integrator 1
     call initial_step_velocity_verlet(coords,v_last,v_half,MD_forces,MD_tstep)
  else if (MD_ensemble.eq.'NVT_parrinello') then
     tsystem = MD_tstep
     call initial_step_velocity_verlet(coords,v_last,v_half,MD_forces,MD_tstep)
  else if (MD_ensemble .eq. 'GLE_thermostat') then
     tsystem = MD_tstep
       !initialize and build colored noise propagators
     call gle_init(MD_tstep/2d0)
     call initial_step_velocity_verlet(coords,v_last,v_half,MD_forces,MD_tstep)
  end if
  call localorb_info(' ')
  !stop
end subroutine initialize_MD
!******
  
!------------------------------------------------------------------------------
!****s* molecular_dynamics/MD_step
!  NAME
!    MD_step
!  SYNOPSIS
subroutine MD_step(Epot)
!  PURPOSE
!    driver routine for single MD step
!  USES
  use runtime_choices
  use geometry
  use species_data
  use constants
  use physics
  use dimensions
  use timing
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o Epot - potential energy of current configuration
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none

  real*8 :: Epot, T_inst, Ekin, shat1, shat2
  integer :: i_coord, i_atom
  real*8, dimension(3,n_atoms) :: MD_forces
  real*8, dimension(n_atoms)   :: px,py,pz
  real*8, dimension(n_atoms)   :: fx,fy,fz
  real*8 :: plumed_energy
  character*100 :: info_str


  !PLUMED force added
  if(plumed_plugin) then
	MD_forces(:,:) = total_forces(:,:)
       
	if(myid.eq.0) then
   		!write(use_unit,*) "fhi forces before", MD_forces
     		!call meta_force_calculation_(n_periodic,lattice_vector,coords,MD_forces,MD_stepcount)
                ! 1.3 interface
                do i_atom = 1, n_atoms
                        px(i_atom)=coords(1,i_atom)
                        py(i_atom)=coords(2,i_atom)
                        pz(i_atom)=coords(3,i_atom)
                        fx(i_atom)=MD_forces(1,i_atom)
                        fy(i_atom)=MD_forces(2,i_atom)
                        fz(i_atom)=MD_forces(3,i_atom)
                end do
                call meta_force_calculation_( lattice_vector, MD_stepcount , px, py, pz , fx, fy, fz,  plumed_energy)
		!write(use_unit,*) "fhi n_periodic, MD_stepcount", n_periodic, MD_stepcount
		!write(use_unit,*) "fhi coords", coords
        	do i_atom = 1, n_atoms
                        MD_forces(1,i_atom)=fx(i_atom)
                        MD_forces(2,i_atom)=fy(i_atom)
                        MD_forces(3,i_atom)=fz(i_atom)
		end do
	
  	endif	
	call broadcast_MD_velocities(MD_forces,0)
	MD_forces(:,:) = MD_forces(:,:)/MD_KE_factor	
  else
     ! change forces to MD integrator units of [amu bohr/ps^2]
     MD_forces(:,:) = total_forces(:,:)/MD_KE_factor

     ! delete any forces on atoms that should not move:
     if (use_relaxation_constraints) then 
        do i_atom = 1, n_atoms
           if (constrain_MD(i_atom)) MD_forces(:,i_atom) = 0d0
        end do
     end if
  endif	

  if ((MD_ensemble.eq.'NVE').or.(MD_ensemble.eq.'NVT_berendsen').or.(MD_ensemble.eq.'NVE_damped')&
        .or.(MD_ensemble.eq.'NVT_andersen').or.(MD_ensemble.eq.'NVT_parrinello').or.(MD_ensemble.eq.'GLE_thermostat')) then
     tsystem      = dble(MD_stepcount)*MD_tstep
     MD_stepcount = MD_stepcount + 1
     r_last(:,:)  = coords(:,:)
     tsystem_last = tsystem     
     MD_Epot_last = Epot
     !  integrator for Newtonian(-like) MD
     call velocity_verlet(coords,v_last,v_half,MD_forces,MD_tstep,Epot)
  else if (MD_ensemble.eq.'NVE_4th_order') then 
     call integrate_fourth_order_NVE(MD_tstep,coords,r_last,v_half,v_last,MD_forces)
     if (MD_high_order_i.eq.1) then
        tsystem      = tsystem + MD_tstep
        tsystem_last = tsystem
        MD_stepcount = MD_stepcount + 1 
        MD_Epot_last = Epot
     end if
  else if (MD_ensemble.eq.'NVT_nose-poincare') then 
     ! keep variables for eventual output once the velocity has been calculated
     s_NP_last    = s_NP
     MD_stepcount   = MD_stepcount + 1
     r_last(:,:)    = coords(:,:)
     tsystem_last   = tsystem     
     MD_Epot_last   = Epot
     ! MD integration step
!     call NP99_generalized_leap_frog(MD_tstep,MD_Q_NP,coords,v_half,s_NP, &
!          s_dot_NP_half,v_last,s_dot_NP_last,MD_forces,Epot,MD_H0,        &
!          MD_temperature,tsystem_last)  
     call NP01_generalized_leap_frog(MD_tstep,MD_Q_NP,coords,v_half,   &
          s_NP,s_NP_last,s_dot_NP_half,v_last,s_dot_NP_last,MD_forces, &
          Epot,MD_Epot_last,MD_H0,MD_temperature)
     ! calculate system time (rather than simulation time)
     tsystem = tsystem + MD_tstep
  else if (MD_ensemble.eq.'NVT_nose-hoover') then
     MD_stepcount   = MD_stepcount + 1
     r_last(:,:)    = coords(:,:)
     tsystem_last   = tsystem     
     MD_Epot_last   = Epot
     call NH96_generalized_leap_frog(MD_tstep,MD_Q_NP,coords,v_half,v_last, &
          s_dot_NP_half,s_NP_half,MD_forces,MD_temperature,shat1,shat2,     &
          s_NP_last,s_dot_NP_last)  
     tsystem = tsystem + MD_tstep
  end if

end subroutine MD_step
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/momentum_rescale_berendsen
!  NAME
!    momentum_rescale_berendsen - initialization of molecular dynamics
!  SYNOPSIS
subroutine momentum_rescale_berendsen(v,T_inst,T_system, MD_tstep, tau)
!  PURPOSE
!    cheap and simple implementation of the Berendsen method for thermal equilibration
!  USES
  use dimensions
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o v - current velocity
!    o T_inst - instantaneous temperature
!    o T_system - desired system temperature
!    o tau - relaxation time
!  OUTPUT
!    o v - corrected velocities
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  real*8, dimension(3,n_atoms) :: v
  real*8 ::T_inst, T_system, tau, prefactor, MD_tstep
  prefactor = sqrt(1d0+(MD_tstep/tau)*(T_system/T_inst-1d0))
  v(:,:) = prefactor*v(:,:)
end subroutine momentum_rescale_berendsen
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/momentum_rescale_NVE_damped
!  NAME
!    momentum_rescale_NVE_damped
!  SYNOPSIS
subroutine momentum_rescale_NVE_damped(v)
!  PURPOSE
!    something simple to do a damped-MD relaxation: simply reduce the velocities by a certain factor
!  USES
  use dimensions
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o v - current velocities
!  OUTPUT
!    o v - damped velocities
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  real*8, dimension(3,n_atoms) :: v
  v(:,:) = NVE_damping_factor*v(:,:)
end subroutine momentum_rescale_NVE_damped
!******

!-----------------------------------------------------------------------------------------------------------------------
!
! Integrators for Newtonian MD
!
!-----------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!****s* molecular_dynamics/velocity_verlet
!  NAME
!    velocity_verlet
!  SYNOPSIS
subroutine velocity_verlet(r,v_last,v_half,f,deltat,uu)
!  PURPOSE
!    integrator for regular Newtonian mechanics
! (1)  calculate v_last(t) from input data (f,m,v_half)
! (2)  calculate deltar(t) from v_last and forces
! (2a) check bounding box properties, to be removed in actual MD
! (3)  calculate r_now(t)
! (4)  calculate v_half(t+deltat/2)
!  USES
  use dimensions
  use geometry
  use species_data
  use runtime_choices
  use constants
  use timing
  use localorb_io
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o r -- initial positions
!    o v_half -- initial velocities
!    o f -- atomic forces
!    o deltat -- time step
!  OUTPUTS
!    o v_last -- velocities at last full time step
!    o r -- updated atomic positions
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  real*8, dimension(3,n_atoms) :: r, v_last, v_half, f, velocity_comparison
  real*8 :: deltat, Ekin, T_inst
  real*8 :: Ekin0_BDP, uu !,MDB_uu !BDP thermostat 
  integer :: i_atom, i_coord
  character*120 :: info_str

  ! (1) Verlet full tstep velocity update: calculate v_last - this will be the velocity _at_ the present time step
  !     (for which we just computed the forces)
  do i_atom = 1, n_atoms
     v_last(:,i_atom)  = v_half(:,i_atom) + f(:,i_atom)*deltat/(2d0*species_m(species(i_atom)))
  enddo
  
  ! (1a) The GLE propagation is properly Trotter factorized, so the full tstep prop. goes here 
  if (MD_ensemble.eq.'GLE_thermostat') then
     call calculate_kinetic_energy(v_last,Ekin)
     langham=langham+Ekin
     call gle_step(v_last)     
     call calculate_kinetic_energy(v_last,Ekin)
     langham=langham-Ekin
     langham_last=langham ! this is the "correct" value of the langevin conserved Hamiltonian!
  end if
 
 
  ! (2) Now that we know the propagated velocities, calculate the effect of
  !     any thermostats that act directly on the velocities
  ! NOTE: at least for BDP one should factorize properly, otherwise the conserved
  ! quantitiy is calculated in a shifted place of the verlet scheme - for NH the same should hold.
  if (MD_ensemble.eq.'NVT_berendsen') then
    call calculate_kinetic_energy(v_last,Ekin)
    T_inst = 2d0*Ekin/(3*dble(n_atoms_MD)*boltzmann_kB)
    call  momentum_rescale_berendsen(v_last,T_inst,MD_temperature, MD_tstep, MD_tau_berendsen)
  else if (MD_ensemble.eq.'NVE_damped') then
    call momentum_rescale_NVE_damped(v_last)
  else if (MD_ensemble.eq.'NVT_andersen') then
    call stochastic_momentum_NVT_andersen(v_last,MD_temperature,MD_tstep,MD_nu_andersen,counter)
  else if (MD_ensemble.eq.'NVT_parrinello') then
    call calculate_kinetic_energy(v_last,Ekin0_BDP)
    T_inst = 2d0*Ekin0_BDP/(3*dble(n_atoms_MD)*boltzmann_kB) ! in K
    call  newscale_BDP(v_last, T_inst, MD_tau_BDP)
  end if

  
  !if (MD_ensemble.eq.'GLE_thermostat') then
  ! langham_half=langham  
  !endif

  ! (3) If requested, remove any residual rotations and translations from the actual velocities.
  !     With numerically exact forces AND FOR NVE and an infinitesimally small time step, no such 
  !     components would arise in the first place, but as a result of
  !     (a) the finite integration grids in DFT
  !     (b) the finite time integration step in MD
  !     energy would still be distributed into all degrees of freedom (i.e., also translations and rotations)
  ! NOTE: if there is a thermostat, it should either really distribute energy in ALL DOF's including translation
  ! and rotation, or it should simply not act on these modes - the way we do it now is not 100% correct.

  if ( (MB_clean_rotations) .or. (n_periodic .gt. 0) ) call clean_velocities(v_last)   ! remove rotations and translations from velocities
 

  if (MD_ensemble.eq.'NVT_parrinello') then ! conserved quantity
    call calculate_kinetic_energy(v_last,Ekin)
    !Ekin0_BDP = Ekin0_BDP - Ekin
    !BDP_hh0 = BDP_hh0 + Ekin0_BDP
    BDP_conint    = BDP_conint + Ekin0_BDP - Ekin
    BDP_consv     = BDP_consv + BDP_conint      ! still needed?
  end if

  ! (3) now here we start calculating v_half
  v_half=v_last
  
  ! (3a) The GLE propagation is properly Trotter factorized, so the half tstep prop. goes here 
  if (MD_ensemble.eq.'GLE_thermostat') then
     call calculate_kinetic_energy(v_half,Ekin)
     langham=langham+Ekin
     call gle_step(v_half)     
     call calculate_kinetic_energy(v_half,Ekin)
     langham=langham-Ekin                                                                                     
  end if


  ! (4) Verlet half tstep velocity update: calculate v_half(t+deltat/2) - this is the velocity at the present time step plus 1/2 delta t,        
  !     where the velocity change is computed based on the forces at the present time step                    
  do i_atom = 1, n_atoms                                                                                      
     v_half(:,i_atom) = v_half(:,i_atom) + f(:,i_atom)*deltat/(2d0*species_m(species(i_atom)))                
  end do
  
  ! MR note: in principle cleaning could be performed also here in order to get molecules that really don't rotate  
  
  ! (5) calculate r new, at (present time step plus delta t) - equivalent to updating it with the current v_half * deltat.
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*(v_half(:,i_atom))! + f(:,i_atom)*deltat/(2d0*species_m(species(i_atom))))
     if (n_periodic.gt.0) then
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end if
  enddo

  ! (6) with the new positions, the new forces are calculated and one can proceed to calculate v_last, velocity at full next tstep.
end subroutine velocity_verlet
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/initial_step_velocity_verlet
!  NAME
!    intial_step_velocity_verlet
!  SYNOPSIS
subroutine initial_step_velocity_verlet(r,v_last,v_half,f,deltat)
!  PURPOSE
!    initial half-step integrator for Newtonian mechanics
!  USES
  use dimensions
  use geometry
  use runtime_choices
  use species_data
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  real*8, dimension(3,n_atoms) :: r, v_last, v_half, f
  real*8 :: deltat, Ekin
  integer :: i_atom, i_coord

  v_half = v_last
! initialize GLE propagation here
  if (MD_ensemble.eq.'GLE_thermostat') then
     call calculate_kinetic_energy(v_half,Ekin)
     langham=langham+Ekin
     call gle_step(v_half)
     call calculate_kinetic_energy(v_half,Ekin)
     langham=langham-Ekin
  end if

  do i_atom = 1, n_atoms
     v_half(:,i_atom) = v_half(:,i_atom) + f(:,i_atom)*deltat/(2d0*species_m(species(i_atom)))
  end do

  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*v_half(:,i_atom)
     if (n_periodic.gt.0) then
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end if
  end do
end subroutine initial_step_velocity_verlet
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/integrate_fourth_order_NVE
!  NAME
!    integrate_fourth_order_NVE
!  SYNOPSIS
subroutine integrate_fourth_order_NVE(deltat,r,r_last,v_half,v_last,forces)
!  PURPOSE
!    common integrator for fourth order symplectic NVE scheme, to be used in highly accurate trajectories
!  USES
  use dimensions
  use geometry
  use species_data
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o deltat - time step
!    o r      - current positions
!    o v_half - current velocities
!    o forces - forces from the last known coordinates
!  OUTPUT
!    o r_last - positions at full time step
!    o v_last - velocities at full time step
!    o r      - current positions
!    o v_half - current velocities
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r, r_last, v_last, v_half, forces
  real*8 :: deltat

  MD_high_order_i = MD_high_order_i + 1
  if (MD_high_order_i.eq.1) then
     do i_atom = 1, n_atoms 
        ! finish last step: calculate velocity and coordinates at end of time step for proper output
        r_last(:,i_atom) = r(:,i_atom)
        v_last(:,i_atom) = v_half(:,i_atom) + & 
                           c_high_order_coeff(MD_max_high_order_index)*deltat*forces(:,i_atom)/species_m(species(i_atom))
     end do 
     if ( (MB_clean_rotations) .or. (n_periodic .gt. 0)) call clean_velocities(v_last)   ! remove rotations and translations from velocities
     do i_atom = 1, n_atoms 
        ! start new step ... 
        v_half(:,i_atom) = v_last(:,i_atom) + c_high_order_coeff(1)*deltat*forces(:,i_atom)/species_m(species(i_atom))
        r     (:,i_atom) = r     (:,i_atom) + d_high_order_coeff(1)*deltat*v_half(:,i_atom)
        if (n_periodic.gt.0) then
          if (use_thermodynamic_integration) then 
            call TDI_map_to_center_cell(r(:,i_atom),i_atom)
          else
            call map_to_center_cell(r(:,i_atom))
          end if
        end if
     end do
  else
     ! (i = 2,3,4)
     do i_atom = 1, n_atoms
        v_half(:,i_atom) = v_half(:,i_atom) & 
                           + c_high_order_coeff(MD_high_order_i)*deltat*forces(:,i_atom)/species_m(species(i_atom))       
        r     (:,i_atom) = r     (:,i_atom) + d_high_order_coeff(MD_high_order_i)*deltat*v_half(:,i_atom)
        if (n_periodic.gt.0) then
          if (use_thermodynamic_integration) then 
            call TDI_map_to_center_cell(r(:,i_atom),i_atom)
          else
            call map_to_center_cell(r(:,i_atom))
          end if
        end if
     end do
  end if

  ! new cycle on next step?
  if (MD_high_order_i.eq.(MD_max_high_order_index-1)) MD_high_order_i = 0

end subroutine integrate_fourth_order_NVE
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/initial_fourth_order_NVE_step
!  NAME
!    initial_fourth_order_NVE_step
!  SYNOPSIS
subroutine initial_fourth_order_NVE_step(deltat,r,v_half,forces)
!  PURPOSE
!    first step & initialization for fourth order NVE integrator
!  USES
  use dimensions 
  use geometry
  use species_data
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o deltat - time step
!    o r      - current positions
!    o v_half - current velocities
!    o forces - forces from the last known coordinates
!  OUTPUT
!    o r      - current positions
!    o v_half - current velocities
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r, r_last, v_last, v_half, forces
  real*8 :: deltat
  MD_high_order_i = 1
  do i_atom = 1, n_atoms
     v_half(:,i_atom) = v_half(:,i_atom) + c_high_order_coeff(1)*deltat*forces(:,i_atom)/species_m(species(i_atom))
     r     (:,i_atom) = r     (:,i_atom) + d_high_order_coeff(1)*deltat*v_half(:,i_atom)
     if (n_periodic.gt.0) then
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end if
  end do
end subroutine initial_fourth_order_NVE_step
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/NP99_generalized_leap_frog
!  NAME
!    NP99_generalized_leap_frog - GLA for Nosé-Poincaré MD
!  SYNOPSIS
subroutine NP99_generalized_leap_frog(deltat,Q,r,p_half,s,sdot_half,p_last,sdot_last, &
     f,Epot,H0,temperature,time)  
!  PURPOSE
!    one integrator for Nosé-Poincaré MD - superceded by NP01_generalized_leap_frog
!  USES
  use constants
  use dimensions
  use geometry
  use species_data
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o deltat - time step
!    o r - atomic positions
!    o p_half - momenta
!    o s - thermostat variable
!    o sdot_half - thermostat momentum
!    o f - DFT forces
!    o Epot - total energy of configuration
!    o H0 - Nosé-Poincaré integration constant
!    o temperature - desired T
!    o time - simulation time
!  OUTPUT
!    o p_half - the updated momenta at half-time step
!    o p_last - the momenta at the last full time step, for output
!    o sdot_last - thermostat momentum at last full time step
!    o sdot_half - updated thermostat momentum
!    o s - updated thermostat variable
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r,p_half,p_last,f
  real*8 :: s, sdot_half, sdot_last, Epot, Ekin, deltat, H0, temperature, threeNKTlnsp1, C, buf, sold, Q, time

  threeNKTlnsp1 = dble(MD_g_DOF)*boltzmann_kB*temperature*(1d0+log(s))
  p_last(:,:)   = p_half(:,:) + deltat*s*f(:,:)/2d0
  call calculate_kinetic_energy_nose(p_half,s,Ekin)
  sdot_last     = sdot_half + deltat*(Ekin - Epot - threeNKTlnsp1 + H0)/(2d0*MD_KE_factor) - deltat*sdot_half*sdot_half/(4d0*Q)
  p_half(:,:)   = p_half(:,:) + deltat*s*f(:,:)
  call calculate_kinetic_energy_nose(p_half,s,Ekin)
  C             = deltat*(threeNKTlnsp1-Ekin+Epot-H0)/(2d0*MD_KE_factor)-sdot_last
  sdot_half     = -2d0*C/(1d0+sqrt(1d0-(deltat*C/Q)))
  buf           = deltat*sdot_half/(2d0*Q)
  sold          = s
  s             = s*(1d0+buf)/(1d0-buf)
  buf           = deltat*(s+sold)/(2d0*s*sold)
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + buf*p_half(:,i_atom)/species_m(species(i_atom))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms 
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end do
  end if
end subroutine NP99_generalized_leap_frog
!******


!------------------------------------------------------------------------------
!****s* molecular_dynamics/initial_step_NP99_GLA
!  NAME
!    initial_step_NP99_GLA
!  SYNOPSIS
subroutine initial_step_NP99_GLA(deltat,Q,r,p_half,s,sdot_half,&
                               p_last,sdot_last,f,Epot,H0,temperature,time)  
!  PURPOSE
!    first MD step for a Nosé-Poincaré thermostat
!  USES
  use constants
  use dimensions
  use species_data
  use geometry
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    deltat - time step
!    Q      - thermostat mass
!    r      - atomic positions
!    p_half - atomic momenta at half-time
!    s      - thermostat value
!    sdot_half - its momentum at half-time
!    p_last - momenta at last full time step
!    sdot_last - momenta at last full time step
!    f      - forces
!    Epot   - total energy
!    temperature
!    time   - initial time 
!  OUTPUT
!    r      - updated positions
!    p_half - updated momenta
!    s      - updated thermostat variable
!    sdot_half - updated thermostat momentum
!    time   - final time
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r,p_half,p_last,f
  real*8 :: s, sdot_half, sdot_last, Epot, Ekin, deltat, H0, temperature, threeNKTlnsp1, C, buf, sold, Q, time
  threeNKTlnsp1 = dble(MD_g_DOF)*boltzmann_kB*temperature*1d0
  p_last(:,:)   = p_half(:,:)
  s             = 1d0
  sdot_last     = 0d0
  p_half(:,:)   = p_half(:,:) + deltat*s*f(:,:)
  call calculate_kinetic_energy_nose(p_half,s,Ekin)
  C             = deltat*(threeNKTlnsp1-Ekin+Epot-H0)/(2d0*MD_KE_factor)-sdot_last
  sdot_half     = -2d0*C/(1d0+sqrt(1d0-(deltat*C/Q)))
  buf           = deltat*sdot_half/(2d0*Q)
  sold          = s
  s             = s*(1d0+buf)/(1d0-buf)
  buf           = deltat*(s+sold)/(2d0*s*sold)
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + buf*p_half(:,i_atom)/species_m(species(i_atom))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms 
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end do
  end if
end subroutine initial_step_NP99_GLA
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/NP01_generalized_leap_frog
!  NAME
!    NP01_generalized_leap_frog - GLA for Nosé-Poincaré MD
!  SYNOPSIS
subroutine NP01_generalized_leap_frog(deltat,Q,r,p_half,s,s_last,&
     sdot_half,p_last,sdot_last,f,Epot,Epot_last,H0,temperature)
!  PURPOSE
!    integrator for Nosé-Poincaré MD
!  USES
  use constants
  use dimensions
  use geometry
  use species_data
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o deltat -- simulation(!) time step
!    o Q -- thermostat mass
!    o r -- current atomic positions
!    o p_half -- velocities at t-dt/2
!    o s -- current thermostat var
!    o sdot_half -- thermostat momentum at t-dt/2
!    o f -- forces for the current time step
!    o Epot -- total energy for current time step
!    o Epot_last -- potential energy of last time step
!    o H0 --  NP integration constant
!    o temperature
!  OUTPUT
!    o r -- updated atomic positions
!    o p_half -- updated velocities at t+dt/2
!    o s_last -- thermostat variable at last full time step
!    o p_last -- velocities for the last time step 
!    o sdot_last -- thermostat momentum for last time step
!    o s -- updated thermostat variable
!    o sdot_half -- updated thermostat momentum at t+dt/2
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r,p_half,p_last,f
  real*8 :: s, sdot_half, sdot_last, Epot, Ekin, deltat, H0, temperature, threeNKTlnsp1
  real*8 :: Epot_last, buf, s_last, Q, time, sdot_doublestar
  threeNKTlnsp1   = dble(MD_g_DOF)*boltzmann_kB*temperature*(1d0+log(s))
  p_last(:,:)     = p_half(:,:) + deltat*s*f(:,:)/2d0
  call calculate_kinetic_energy_nose(p_half,s,Ekin)
  sdot_doublestar = sdot_last/(1d0+sdot_last*deltat/(4d0*Q))+deltat*(Ekin-(Epot+Epot_last)/2d0+H0-threeNKTlnsp1)/MD_KE_factor
  buf             = 1d0+sdot_doublestar*deltat/(4d0*Q)
  s_last          = s*buf*buf
  sdot_last       = sdot_doublestar/buf
  buf             = 1d0 + sdot_last*deltat/(4d0*Q)
  s               = s_last*buf*buf
  p_half(:,:)     = p_last(:,:) + s*f(:,:)*deltat/2d0
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*p_half(:,i_atom)/(s*species_m(species(i_atom)))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms 
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end do
  end if
  Epot_last = Epot
end subroutine NP01_generalized_leap_frog
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/initial_step_NP01_GLA
!  NAME
!    initial_step_NP01_GLA
!  SYNOPSIS
subroutine initial_step_NP01_GLA(deltat,Q,r,p_half,s,sdot_half,p_last,sdot_last,f)
!  PURPOSE
!    first MD step for a Nosé-Poincaré thermostat with NP01 integrator
!  USES
  use constants
  use dimensions
  use species_data
  use geometry
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    deltat - time step
!    Q      - thermostat mass
!    r      - atomic positions
!    p_half - atomic momenta at half-time
!    s      - thermostat value
!    sdot_half - its momentum at half-time
!    p_last - momenta at last full time step
!    sdot_last - momenta at last full time step
!    f      - forces
!    Epot   - total energy
!    temperature
!    time   - initial time 
!  OUTPUT
!    r      - updated positions
!    p_half - updated momenta
!    s      - updated thermostat variable
!    sdot_half - updated thermostat momentum
!    time   - final time
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r,p_half,p_last,f
  real*8 :: s, sdot_half, sdot_last, Q, deltat
  s           = 1d0
  sdot_last   = 0d0
  sdot_half   = 0d0
  p_half(:,:) = p_last(:,:) + deltat*f(:,:)/2d0
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*p_half(:,i_atom)/(s*species_m(species(i_atom)))
  end do
  if (n_periodic.gt.0) then
    do i_atom = 1, n_atoms 
      if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
      else
         call map_to_center_cell(r(:,i_atom))
      end if
    end do
  end if 
end subroutine initial_step_NP01_GLA
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/NH96_generalized_leap_frog
!  NAME
!    NH96_generalized_leap_frog
!  SYNOPSIS
!     call NH96_generalized_leap_frog(MD_tstep,MD_Q_NP,coords,v_half,v_last, &
!          s_dot_NP_half,s_NP_half,MD_forces,MD_temperature,shat1,shat2,     &
!          s_NP_last,s_dot_NP_last)  
subroutine NH96_generalized_leap_frog(deltat,Q,r,p_half,p_last,pi_half,&
     eta_half,f,T,shat1,shat2,eta_last,pi_last)  
!  PURPOSE
!    integrator for Nosé-Hoover MD
!  USES
  use constants
  use dimensions
  use geometry
  use species_data
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!   o deltat -- time step
!   o Q -- thermostat mass
!   o r -- initial positions
!   o p_half -- initial momenta
!   o pi_half -- thermostat momentum
!   o eta_half -- thermostat value 
!   o f -- total forces
!   o T -- temperature
!   o pi_last -- thermostat value at last full time step
!  OUTPUT
!   o eta_last -- thermostat momentum at last full time step
!   o r -- initial positions
!   o p_half -- initial momenta
!   o p_last -- momenta at last full time step
!   o pi_half -- thermostat momentum
!   o eta_half -- thermostat value 
!   o shat1 -- Nosé s
!   o shat2 -- actually need 2 of those ... 
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r,p_half,p_last,f, pstar, pkin
  real*8 :: eta, eta_half, eta_last, Epot, Ekin, deltat, T, gkT
  real*8 :: etaold, Q, time, shat1, shat2, pidoublestar, pi_half, pi_last, pistar
  gkT          = dble(MD_g_DOF)*boltzmann_kB*T
  pstar(:,:)   = p_half(:,:) + deltat*f(:,:)/2d0
  call calculate_kinetic_energy_nose(pstar,1d0,Ekin)
  pidoublestar = pi_half + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  shat1        = exp(-pidoublestar*deltat/(2d0*Q))
  eta_last     = eta_half + pidoublestar*deltat/(2d0*Q)
  p_last(:,:)  = pstar(:,:)*shat1
  call calculate_kinetic_energy_nose(p_last,1d0,Ekin)
  pi_last      = pidoublestar + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  call calculate_kinetic_energy_nose(p_last,1d0,Ekin)
  pistar       = pi_last + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  shat2        = exp(-pistar*deltat/(2d0*Q))
  pkin(:,:)    = p_last(:,:) * shat2
  call calculate_kinetic_energy_nose(pkin,1d0,Ekin)
  pi_half      = pistar + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  eta_half     = eta_last + pistar*deltat/(2d0*Q)
  p_half(:,:)  = p_last(:,:)*shat2+f(:,:)*deltat/2d0

!  MB_clean_rotations = .false.
  !if (MB_clean_rotations) call clean_velocities_NH(p_half)   ! remove rotations and translations from velocities, hopefully	
  if ( (MB_clean_rotations) .or. (n_periodic .gt. 0)) call clean_velocities_NH(p_half) ! remove rotations and translations from velocities, hopefully  

  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*p_half(:,i_atom)/species_m(species(i_atom))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms 
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end do
  end if
end subroutine NH96_generalized_leap_frog
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/initial_step_NH96_GLA
!  NAME
!    initial_step_NH96_GLA - GLA first step for Nosé-Hoover MD
!  SYNOPSIS
subroutine initial_step_NH96_GLA(deltat,Q,r,p_half,p_last,eta_half,pi_half,&
     shat,f,T)
!  PURPOSE
!    first MD step for a Nosé-Poincaré thermostat
!  USES
  use constants
  use dimensions
  use species_data
  use geometry
  use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    deltat - time step
!    Q      - thermostat mass
!    r      - atomic positions
!    p_half - atomic momenta at half-time
!    s      - thermostat value
!    sdot_half - its momentum at half-time
!    p_last - momenta at last full time step
!    sdot_last - momenta at last full time step
!    f      - forces
!    Epot   - total energy
!    temperature
!    time   - initial time 
!  OUTPUT
!    r      - updated positions
!    p_half - updated momenta
!    s      - updated thermostat variable
!    sdot_half - updated thermostat momentum
!    time   - final time
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
! SOURCE
  implicit none
  integer :: i_atom
  real*8, dimension(3,n_atoms) :: r,p_half,p_last,f, pkin
  real*8 :: Ekin, deltat, T, Q, time, eta_half, pi_half, shat, pistar,gkT, buf
  p_last(:,:) = p_half(:,:)
  gkT         = dble(MD_g_DOF)*boltzmann_kB*T 
  call calculate_kinetic_energy_nose(p_half,1d0,Ekin)
  pistar      = deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  eta_half    = pistar*deltat/(2d0*Q)
  shat        = exp(-eta_half)
  pkin(:,:)   = p_half(:,:)*shat
  call calculate_kinetic_energy_nose(pkin,1d0,Ekin)
  pi_half     = pistar + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  p_half(:,:) = p_half(:,:)*shat + f(:,:)*deltat/2d0
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*p_half(:,i_atom)/species_m(species(i_atom))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms 
       if (use_thermodynamic_integration) then 
         call TDI_map_to_center_cell(r(:,i_atom),i_atom)
       else
         call map_to_center_cell(r(:,i_atom))
       end if
     end do
  end if
end subroutine initial_step_NH96_GLA
!******
!------------------------------------------------------------------------------
!****s* molecular_dynamics/calculate_kinetic_energy_nose
!  NAME
!    calculate_kinetic_energy_nose 
!  SYNOPSIS
subroutine calculate_kinetic_energy_nose(p,s,Ekin)
!  PURPOSE
!    kinetic energy with Nose-Hamiltonian
!  USES
  use dimensions
  use constants
  use geometry 
  use species_data
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    p - momenta
!    s - thermostat variable
!  OUTPUT
!    Ekin - kinetic energy
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom 
  real*8, dimension(3,n_atoms) :: p
  real*8 :: s,Ekin
  Ekin = 0d0
  do i_atom = 1, n_atoms
     Ekin = Ekin + (p(1,i_atom)*p(1,i_atom)+p(2,i_atom)*p(2,i_atom)+p(3,i_atom)*p(3,i_atom))/species_m(species(i_atom))
  end do
  Ekin = MD_KE_factor*Ekin/(2d0*s*s)
end subroutine calculate_kinetic_energy_nose
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/calculate_kinetic_energy
!  NAME
!    calculate_kinetic_energy
!  SYNOPSIS
subroutine calculate_kinetic_energy(v,Ekin,start_i,end_i)
!  PURPOSE
!    kinetic energy
!  USES
  use dimensions
  use constants
  use species_data
  use geometry
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    v - velocities
!  OUTPUT
!    Ekin - kinetic energy
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom
  integer, optional :: start_i,end_i
  real*8, dimension(3,n_atoms) :: v
  real*8 :: Ekin
  Ekin = 0d0
  if ((present(start_i) ).and.(present(end_i) )) then
    do i_atom = start_i,end_i
       Ekin = Ekin+species_m(species(i_atom))*(v(1,i_atom)*v(1,i_atom)+v(2,i_atom)*v(2,i_atom)+v(3,i_atom)*v(3,i_atom))
    end do
  else
    do i_atom = 1, n_atoms
       Ekin = Ekin+species_m(species(i_atom))*(v(1,i_atom)*v(1,i_atom)+v(2,i_atom)*v(2,i_atom)+v(3,i_atom)*v(3,i_atom))
    end do
  end if
  Ekin = MD_KE_factor*Ekin/2d0
end subroutine calculate_kinetic_energy
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/clean_velocities
!  NAME
!    clean_velocities
!  SYNOPSIS
subroutine clean_velocities(v, clean_despite_constraints)
!  PURPOSE
!    clean center-of-mass motion from the velocities
!  USES
  use dimensions
  use geometry
  use species_data
  use runtime_choices
  use localorb_io
  use relaxation
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    v - initial velocities
!    clean_despite_constraints - if .true., the routine will clean the
!                                velocities regardless of any constrained atoms
!  OUTPUT
!    v - cleaned velocities
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  logical, optional :: clean_despite_constraints
  integer :: i_atom, lwork, info, i_coord, i_coord_p1, i_coord_p2, i_coord_2
  real*8, dimension(3,n_atoms) :: v, momenta, relative_coords
  real*8, dimension(3) :: center_of_mass, moments_of_inertia, p_total
  real*8, dimension(3,n_atoms,3) :: translation_vectors, rotation_vectors
  real*8, dimension(:), allocatable :: work
  real*8, dimension(3,3) :: tensor_of_inertia
  real*8 :: inv_norm, norm, rotation_component, translation_component, dnrm2, ddot
  real*8 :: E_kin_in, E_kin_out, velocity_scale
  character*100 :: info_str
  logical :: run_cleaning

  ! run routine in the case of constrained atoms? - not unless it is requested explicitly!
  run_cleaning = .true.
  if (use_relaxation_constraints) then 
     if (present(clean_despite_constraints)) then 
        run_cleaning = clean_despite_constraints
     else
        run_cleaning = .false.
     end if
  end if

  if (run_cleaning) then 
     if(output_level .eq. 'MD_light') output_priority = 2
     ! get initial kinetic energy
     call calculate_kinetic_energy(v,E_kin_in)
     if (.not.MB_clean_rotations) then
        write(info_str,'(2X,A)') 'Cleaning translations from velocities'
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        ! simple hack to remove all translations from the momenta
        p_total = 0d0
        do i_atom = 1, n_atoms
           p_total(:) = p_total(:) + species_m(species(i_atom))*v(:,i_atom)
        end do
        p_total(:) = p_total(:)/dble(n_atoms_MD)
        do i_atom = 1, n_atoms
           if (.not.constrain_MD(i_atom)) then  ! prevent constrained atoms from acquiring velocities ... 
              v(:,i_atom) = v(:,i_atom) - p_total(:)/species_m(species(i_atom))
           end if
        end do
     else 
        ! apply fully fledged sayvetz conditions to the momenta in order to
        ! project out rotations and translations at the same time. 
        ! In principle, this is the same as (and copied from) the routines
        ! initialize_cleaning_forces and remove_translation_and_rotation
        ! in the module translation, but they had to be modified slightly 
        ! to work for momenta
        write(info_str,'(2X,A)') 'Cleaning translations and rotations from velocities'
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        ! initialize the cleaning
        allocate(work(1))
        tensor_of_inertia = 1.d0
        lwork = -1
        call dsyev('V', 'U', 3, tensor_of_inertia, 3, moments_of_inertia, work, lwork, info)
        lwork = work(1)
        deallocate(work)
        allocate(work(lwork))
        translation_vectors = 0.d0
        inv_norm = 1.d0 / dsqrt(dble(n_atoms))  
        do i_atom = 1, n_atoms, 1
           do i_coord = 1, 3, 1
              translation_vectors(i_coord, i_atom, i_coord) = inv_norm
           end do
        end do        
        ! calculate momentum for each atom
        do i_atom = 1, n_atoms
           momenta(:,i_atom) = v(:,i_atom)*species_m(species(i_atom))
        end do
        ! calculate relative coordinates & center of mass
        center_of_mass = 0d0
        norm           = 0d0
        do i_atom = 1, n_atoms, 1
           center_of_mass(:) = center_of_mass(:) + coords(:, i_atom)
        end do
        center_of_mass(:) = center_of_mass(:) / dble(n_atoms)
        do i_atom = 1, n_atoms, 1
           relative_coords(:, i_atom) = coords(:, i_atom) - center_of_mass(:)
        end do
        
        ! calculate total angular momentum and total momentum; use same routines as for the relaxation
        call get_net_force(momenta,p_total)
        write(info_str,'(2X,A,3E14.4E3)') '| Net linear  momentum before cleaning: ',(p_total(i_coord),i_coord = 1, 3)     
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        call get_net_torque(momenta,relative_coords,p_total)
        write(info_str,'(2X,A,3E14.4E3)') '| Net angular momentum before cleaning: ',(p_total(i_coord),i_coord = 1, 3)
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        
        ! now, get moments of inertia
        tensor_of_inertia = 0.d0
        do i_atom = 1, n_atoms, 1
           do i_coord = 1, 3, 1
              i_coord_p1 = mod(i_coord    , 3) + 1
              i_coord_p2 = mod(i_coord + 1, 3) + 1
              ! diagonal elements
              tensor_of_inertia(i_coord, i_coord) = tensor_of_inertia(i_coord, i_coord) &
                   + (relative_coords(i_coord_p1, i_atom) * relative_coords(i_coord_p1, i_atom)) + &
                   (relative_coords(i_coord_p2, i_atom) * relative_coords(i_coord_p2, i_atom))
              ! non-diagonal elements
              do i_coord_2 = i_coord + 1, 3, 1
                 tensor_of_inertia(i_coord, i_coord_2) = tensor_of_inertia(i_coord, i_coord_2) &
                      - relative_coords(i_coord, i_atom) * relative_coords(i_coord_2, i_atom) 
              end do
           end do
        end do
        ! diagonalize it 
        call dsyev('V', 'U', 3, tensor_of_inertia, 3, moments_of_inertia, work, lwork, info)
        ! generate vectors of rotation according to gaussian
        do i_coord = 1, 3, 1
           i_coord_p1 = mod(i_coord    , 3) + 1
           i_coord_p2 = mod(i_coord + 1, 3) + 1
           do i_atom = 1, n_atoms, 1
              do i_coord_2 = 1, 3, 1
                 rotation_vectors(i_coord_2, i_atom, i_coord) = &
                      (ddot(3, relative_coords(:,i_atom), 1, tensor_of_inertia(:        , i_coord_p1), 1) * &
                      tensor_of_inertia(i_coord_2, i_coord_p2) - &
                      ddot(3, relative_coords(:,i_atom), 1, tensor_of_inertia(:        , i_coord_p2), 1) * &
                      tensor_of_inertia(i_coord_2, i_coord_p1))
                 if (i_coord .eq. 2) then
                    rotation_vectors(i_coord_2, i_atom, i_coord) = - rotation_vectors(i_coord_2, i_atom, i_coord)
                 end if
              end do
           end do
           ! normalize them
           norm = dnrm2(3*n_atoms, rotation_vectors(:, :, i_coord), 1)
           if (norm .gt. 0) then
              inv_norm = 1.d0 / norm
              rotation_vectors(:, :, i_coord) = rotation_vectors(:, :, i_coord) * inv_norm
              ! project out rotation
              rotation_component = ddot(3*n_atoms, rotation_vectors(1,1,i_coord), 1, momenta, 1)
              momenta(:,:) = momenta(:,:) - rotation_component * rotation_vectors(:,:,i_coord)           
           end if
        end do
        ! project out translation
        do i_coord = 1, 3, 1
           translation_component = ddot(3*n_atoms, translation_vectors(1,1,i_coord), 1, momenta, 1)
           momenta(:,:) = momenta(:,:) - translation_component * translation_vectors(:,:,i_coord)
        end do
        
        ! calculate total angular momentum
        call get_net_force(momenta,p_total)
        write(info_str,'(2X,A,3E14.4E3)') '| Net linear  momentum after cleaning:  ',(p_total(i_coord),i_coord = 1, 3)     
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        call get_net_torque(momenta,relative_coords,p_total)
        write(info_str,'(2X,A,3E14.4E3)') '| Net angular momentum after cleaning:  ',(p_total(i_coord),i_coord = 1, 3)
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        
        ! calculate velocities from the momenta
        do i_atom = 1, n_atoms
           v(:,i_atom) = momenta(:,i_atom)/species_m(species(i_atom))
        end do
        
        deallocate(work)
     end if
     ! After removing artificial rotations and translations, we have now removed
     ! a small amount of kinetic energy from the system (from specific degrees of freedom).
     ! We here choose to restore that "lost" kinetic energy by an a-posteriori
     ! rescaling of all velocities.
     ! Note that the associated scaling factor should go to 1 as we reduce the time step to zero,
     ! and possibly as we increase our integration grid density to infinity.
     call calculate_kinetic_energy(v,E_kin_out)
     if(E_kin_out==0) then
        ! catch a (very unlikely!) division by zero
        velocity_scale = 1.
     else
        velocity_scale = sqrt (E_kin_in/E_kin_out)
     endif
     write(info_str,'(2X,A,F12.8)') & 
          '| Velocity rescaling factor to preserve kinetic energy before and after cleaning: ', velocity_scale
     call localorb_info(info_str, use_unit,'(A)',OL_norm)
     
     ! actually rescale velocity
     v = v * velocity_scale
     if(output_level .eq. 'MD_light') output_priority = 1
  end if
end subroutine clean_velocities
!******


!------------------------------------------------------------------------------
!****s* molecular_dynamics/clean_velocities
!  NAME
!    clean_velocities
!  SYNOPSIS
subroutine clean_velocities_NH(v)
!  PURPOSE
!    clean center-of-mass motion from the velocities for NH thermostat
!  USES
  use dimensions
  use geometry
  use species_data
  use runtime_choices
  use localorb_io
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    v - initial velocities
!  OUTPUT
!    v - cleaned velocities
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom, lwork, info, i_coord, i_coord_p1, i_coord_p2, i_coord_2
  real*8, dimension(3,n_atoms) :: v, momenta, relative_coords
  real*8, dimension(3) :: center_of_mass, moments_of_inertia, p_total
  real*8, dimension(3,n_atoms,3) :: translation_vectors, rotation_vectors
  real*8, dimension(:), allocatable :: work
  real*8, dimension(3,3) :: tensor_of_inertia
  real*8 :: inv_norm, norm, rotation_component, translation_component, dnrm2, ddot
  real*8 :: E_kin_in, E_kin_out, velocity_scale
  character*100 :: info_str
  ! get initial kinetic energy
  call calculate_kinetic_energy_nose(v,1d0,E_kin_in) 
  if (.not.MB_clean_rotations) then
     write(info_str,'(2X,A)') 'Cleaning translations from velocities'
     call localorb_info(info_str)
     ! simple hack to remove all translations from the momenta
     p_total = 0d0
     do i_atom = 1, n_atoms
        p_total(:) = p_total(:) + v(:,i_atom)
     end do
     p_total(:) = p_total(:)/dble(n_atoms)
     do i_atom = 1, n_atoms
        v(:,i_atom) = v(:,i_atom) - p_total(:)
     end do
  else 
     ! apply fully fledged sayvetz conditions to the momenta in order to
     ! project out rotations and translations at the same time. 
     ! In principle, this is the same as (and copied from) the routines
     ! initialize_cleaning_forces and remove_translation_and_rotation
     ! in the module translation, but they had to be modified slightly 
     ! to work for momenta
     write(info_str,'(2X,A)') 'Cleaning translations and rotations from velocities'
     call localorb_info(info_str)
     ! initialize the cleaning
     allocate(work(1))
     tensor_of_inertia = 1.d0
     lwork = -1
     call dsyev('V', 'U', 3, tensor_of_inertia, 3, moments_of_inertia, work, lwork, info)
     lwork = work(1)
     deallocate(work)
     allocate(work(lwork))
     translation_vectors = 0.d0
     inv_norm = 1.d0 / dsqrt(dble(n_atoms))  
     do i_atom = 1, n_atoms, 1
        do i_coord = 1, 3, 1
           translation_vectors(i_coord, i_atom, i_coord) = inv_norm
        end do
     end do
     ! calculate momentum for each atom
     do i_atom = 1, n_atoms
        momenta(:,i_atom) = v(:,i_atom)!*species_m(species(i_atom))
     end do
     ! calculate relative coordinates & center of mass
     center_of_mass = 0d0
     norm           = 0d0
     do i_atom = 1, n_atoms, 1
        center_of_mass(:) = center_of_mass(:) + coords(:, i_atom)
     end do
     center_of_mass(:) = center_of_mass(:) / dble(n_atoms)
     do i_atom = 1, n_atoms, 1
        relative_coords(:, i_atom) = coords(:, i_atom) - center_of_mass(:)
     end do
     
     ! calculate total angular momentum and total momentum; use same routines as for the relaxation
     call get_net_force(v,p_total)
     write(info_str,'(2X,A,3E14.4E3)') '| Net linear  momentum before cleaning: ',(p_total(i_coord),i_coord = 1, 3)     
     call localorb_info(info_str)
     call get_net_torque(momenta,relative_coords,p_total)
     write(info_str,'(2X,A,3E14.4E3)') '| Net angular momentum before cleaning: ',(p_total(i_coord),i_coord = 1, 3)
     call localorb_info(info_str)

     ! now, get moments of inertia
     tensor_of_inertia = 0.d0
     do i_atom = 1, n_atoms, 1
        do i_coord = 1, 3, 1
           i_coord_p1 = mod(i_coord    , 3) + 1
           i_coord_p2 = mod(i_coord + 1, 3) + 1
           ! diagonal elements
           tensor_of_inertia(i_coord, i_coord) = tensor_of_inertia(i_coord, i_coord) &
                + (relative_coords(i_coord_p1, i_atom) * relative_coords(i_coord_p1, i_atom)) + &
                  (relative_coords(i_coord_p2, i_atom) * relative_coords(i_coord_p2, i_atom))
           ! non-diagonal elements
           do i_coord_2 = i_coord + 1, 3, 1
              tensor_of_inertia(i_coord, i_coord_2) = tensor_of_inertia(i_coord, i_coord_2) &
                   - relative_coords(i_coord, i_atom) * relative_coords(i_coord_2, i_atom) 
           end do
        end do
     end do
     ! diagonalize it 
     call dsyev('V', 'U', 3, tensor_of_inertia, 3, moments_of_inertia, work, lwork, info)
     ! generate vectors of rotation according to gaussian
     do i_coord = 1, 3, 1
        i_coord_p1 = mod(i_coord    , 3) + 1
        i_coord_p2 = mod(i_coord + 1, 3) + 1
        do i_atom = 1, n_atoms, 1
           do i_coord_2 = 1, 3, 1
              rotation_vectors(i_coord_2, i_atom, i_coord) = &
                   (ddot(3, relative_coords(:,i_atom), 1, tensor_of_inertia(:        , i_coord_p1), 1) * &
                                                          tensor_of_inertia(i_coord_2, i_coord_p2) - &
                    ddot(3, relative_coords(:,i_atom), 1, tensor_of_inertia(:        , i_coord_p2), 1) * &
                                                          tensor_of_inertia(i_coord_2, i_coord_p1))
              if (i_coord .eq. 2) then
                 rotation_vectors(i_coord_2, i_atom, i_coord) = - rotation_vectors(i_coord_2, i_atom, i_coord)
              end if
           end do
        end do
        ! normalize them
        norm = dnrm2(3*n_atoms, rotation_vectors(:, :, i_coord), 1)
        if (norm .gt. 0) then
           inv_norm = 1.d0 / norm
           rotation_vectors(:, :, i_coord) = rotation_vectors(:, :, i_coord) * inv_norm
           ! project out rotation
           rotation_component = ddot(3*n_atoms, rotation_vectors(1,1,i_coord), 1, momenta, 1)
           momenta(:,:) = momenta(:,:) - rotation_component * rotation_vectors(:,:,i_coord)           
        end if
     end do
     
     ! project out translation
     do i_coord = 1, 3, 1
        translation_component = ddot(3*n_atoms, translation_vectors(1,1,i_coord), 1, momenta, 1)
        momenta(:,:) = momenta(:,:) - translation_component * translation_vectors(:,:,i_coord)
     end do
     
     ! calculate total angular momentum
     call get_net_force(momenta,p_total)
     write(info_str,'(2X,A,3E14.4E3)') '| Net linear  momentum after cleaning:  ',(p_total(i_coord),i_coord = 1, 3)     
     call localorb_info(info_str)
     call get_net_torque(momenta,relative_coords,p_total)
     write(info_str,'(2X,A,3E14.4E3)') '| Net angular momentum after cleaning:  ',(p_total(i_coord),i_coord = 1, 3)
     call localorb_info(info_str)

     ! calculate velocities from the momenta
     do i_atom = 1, n_atoms
        v(:,i_atom) = momenta(:,i_atom)!/species_m(species(i_atom))
     end do     
     
     deallocate(work)
  end if

  ! After removing artificial rotations and translations, we have now removed
  ! a small amount of kinetic energy from the system (from specific degrees of freedom).
  ! We here choose to restore that "lost" kinetic energy by an a-posteriori
  ! rescaling of all velocities.
  ! Note that the associated scaling factor should go to 1 as we reduce the time step to zero,
  ! and possibly as we increase our integration grid density to infinity.
  call calculate_kinetic_energy_nose(v,1d0,E_kin_out)
  if(E_kin_out==0) then
     ! catch a (very unlikely!) division by zero
     velocity_scale = 1.
  else
     velocity_scale = sqrt (E_kin_in/E_kin_out)
  endif
  write(info_str,'(2X,A,F12.8)') & 
  '| Velocity rescaling factor to preserve kinetic energy before and after cleaning: ', velocity_scale
  call localorb_info(info_str)

  ! actually rescale velocity
  v = v * velocity_scale

end subroutine clean_velocities_NH
!******



!------------------------------------------------------------------------------
!****s* molecular_dynamics/clean_momenta
!  NAME
!    clean_momenta
!  SYNOPSIS
subroutine clean_momenta(n_atoms, p)
!  PURPOSE
!    clean center-of-mass motion from the momenta
!  USES
!    none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    n_atoms - number of atoms
!    p       - initial momenta
!  OUTPUT
!    p       - cleaned momenta
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: n_atoms, i_atom
  real*8, dimension(3,n_atoms) :: p
  real*8, dimension(3) :: p_total
  p_total = 0d0
  do i_atom = 1, n_atoms
     p_total(:) = p_total(:) + p(:,i_atom)
  end do
  p_total(:) = p_total(:)/dble(n_atoms)
  do i_atom = 1, n_atoms
     p(:,i_atom) = p(:,i_atom) - p_total(:)
  end do
end subroutine clean_momenta
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/stochastic_momentum_NVT_andersen
!  NAME
!    stochastic_momentum_NVT_andersen
!  SYNOPSIS

subroutine stochastic_momentum_NVT_andersen(velocity,temperature,deltat,nu_andersen,counter)
!  PURPOSE
!    selects a stochastic new momentum from a Maxwell-Boltzmann distribution
!  USES
  use runtime_choices
  use species_data
  use dimensions
  use geometry
  use mpi_tasks
  use synchronize_mpi
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    v_last
!  OUTPUT
!    new v_last
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  integer :: i_atom, i_coord
  integer, dimension(n_atoms) :: counter
  real*8, dimension(3,n_atoms) :: velocity !, v_comparison
  real*8 :: temperature, deltat, nu_andersen
  real*8 :: random_number
  if (myid.eq.0) then
     ! only apply random thermostat on thread number 0
     counter = 0
     do i_atom = 1, n_atoms
        if (.not.constrain_MD(i_atom)) then
           random_number = random_own()
           if (random_number .lt. deltat*nu_andersen) then
              do i_coord = 1, 3
                 call MB_sample_box_mueller(temperature,species_m(species(i_atom)),velocity(i_coord,i_atom))
              end do
              counter(i_atom) = 1
           end if
        end if
     end do
  end if
  ! broadcast new velocity to all other threads
  call broadcast_MD_velocities(velocity,0)
end subroutine stochastic_momentum_NVT_andersen
!******


!------------------------------------------------------------------------------
!****s* molecular_dynamics/newscale_BDP
!  NAME
!    newscale_BDP
!  SYNOPSIS
subroutine newscale_BDP(velocity,kki,tau)
! velocity velocity to be rescaled
! kki      istantaneous temperature, in units of the external temperature
! tau     relaxation time, in units of the timestep (i.e. in units of 'how often this routine is called')
!         negative value means no rescale
!  PURPOSE
!    calculate the velocity scaling factor according to Bussi Donadio Parrinello thermostat, new simpler version
!  USES
!  AUTHOR
!    adapted from a Davide Donadio source
!  HISTORY
!    Release version, FHI-aims (2009).
!  INPUTS
!    kki,tau
!  OUTPUT
!    velocity
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  use runtime_choices
  use dimensions
  use timing
  use localorb_io
  use mpi_tasks
  use synchronize_mpi
!  use random_module

  implicit none
  real*8, dimension(3,n_atoms) :: velocity 
  real*8,  intent(in)  :: kki
  !real*8,  intent(in)  :: uu
  real*8,  intent(in)  :: tau
  character*120 :: info_str

  real*8 :: kk0	! actual value during the evolution
  real*8 :: kkf	! final value for the kinetic energy
  !real*8 :: factor

  if (myid.eq.0) then
      if(kki.gt.1.d-6) then
          kkf=0.0d0
          kkf=resamplekin(kki,tau)
          BDP_factor=sqrt(kkf/kki)
      else
        BDP_factor=1.0
        kkf=kki
      end if
      velocity(:,:) = BDP_factor*velocity(:,:)
  end if
  ! broadcast new velocity to all other threads
  call broadcast_MD_velocities(velocity,0)

!     conserved quantity
!     in practise, I subtract the change in total kinetic energy due to the thermostat
!      conint = conint - (sumke-sumke0)
!      consv  = consv + conint
!     write(info_str,'(X,a,i10,f18.9,f18.9,f18.9)') 'tine factor', &
!     MD_stepcount,factor,kki,kkf
!          call localorb_info(info_str)

end subroutine newscale_BDP
!******

!------------------------------------------------------------------------------
!****s* molecular_dynamics/MB_sample_box_mueller
!  NAME
!    MB_sample_box_mueller
!  SYNOPSIS
subroutine MB_sample_box_mueller(T,m,v)
!  PURPOSE
!    Gaussian distribution from uniform distribution of random numbers
!  USES
  use constants
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    seed - RNG seed
!    T    - temperature
!    m    - mass of particle
!  OUTPUT
!    v    - single velocity component
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  real*8 :: T,m,v
  if (.not.MB_sample_box_mueller_saved) then
     MB_sample_box_mueller_saved = .true.
     MB_sample_box_mueller_r1    = random_own()
     MB_sample_box_mueller_r2    = random_own()
     v     = sqrt(boltzmann_kB*T/(MD_KE_factor*m))*sqrt(-2d0*dlog(MB_sample_box_mueller_r1))*cos(2d0*pi*MB_sample_box_mueller_r2)
  else
     MB_sample_box_mueller_saved = .false.
     v     = sqrt(boltzmann_kB*T/(MD_KE_factor*m))*sqrt(-2d0*dlog(MB_sample_box_mueller_r1))*sin(2d0*pi*MB_sample_box_mueller_r2)
  end if
end subroutine MB_sample_box_mueller
!******

!****s* molecular_dynamics/random_own
!  NAME
!    random_own
!  SYNOPSIS
real*8  function random_own()
!  PURPOSE
!    really cheap random number generator
!  USES
  use runtime_choices
  use mpi_tasks
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    seed - initial variable to seed the RNG
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
  implicit none
! SOURCE
  integer :: RNG_num
  real*8  :: rnum
  logical :: user_seeded
  character*8  :: cdate
  character*10 :: ctime
  ! hardwired warning in case someone ever tries again to call the
  ! rng without synchronizing between different tasks
  if (myid.ne.0) then
    write(use_unit,'(1X,A,I8,A)') "* Warning: Random number generator invoked from task myid = ", myid, "."
    write(use_unit,'(1X,A)')      "* To guarantee that all results are properly synchronized between different threads, "
    write(use_unit,'(1X,A)')      "* random_own() must never be called from any other task than myid = 0 . "
    write(use_unit,'(1X,A)')      "* This appears to be a programming error - please correct."
    stop
  end if
  if (MD_RNG_firstcall) then
     ! seed was already set when rng was initialized
     random_own_Y = seed
     MD_RNG_firstcall = .false.
  end if
  if (.not.MD_RNG_seeded) then
     ! if a seed was not specified, we randomize the rng variables used also
     call date_and_time(cdate, ctime)
     read(ctime,'(F10.3)') rnum
     RNG_num = mod(int(rnum),7)+1
  else
     ! if a seed was specified, we keep the rng variables fixed also
     ! this ensures always the same random number sequence for the same seed
     RNG_num = 7
  end if
  random_own_Y = mod(RNG_a(RNG_num)*random_own_Y+RNG_c(RNG_num),RNG_m(RNG_num))
  random_own   = dble(random_own_Y)/dble(RNG_m(RNG_num))
end function random_own
!******


!------------------------------------------------------------------------------
!****s* molecular_dynamics/initialize_RNG
!  NAME
!    initialize_RNG
!  SYNOPSIS
subroutine initialize_RNG
!  PURPOSE
!    initialization of random number generator
!  USES
  use runtime_choices
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none
  character*8  :: cdate
  character*10 :: ctime
  real*8       :: rng_seed
  integer*8    :: temp, i
  character*120 :: info_str
  ! some compilers handle 8 byte integers but can not interpret them
  ! as constants when you type out the digits
  ! Fix by calculating the RNG initialization from "scratch"
  temp = 1
  do i = 1, 31
     temp = temp * 2
  end do
  RNG_m(:) = temp*2 ! 4294967296
  RNG_m(7) = temp-1 ! 2147483647
  RNG_a(1) = 1664525
  RNG_a(2) = 22695477
  RNG_a(3) = 69069
  RNG_a(4) = 1103515245
  RNG_a(5) = 134775813
  RNG_a(6) = 214013
  RNG_a(7) = 16807
  RNG_c(1) = 1013904223
  RNG_c(2) = 1
  RNG_c(3) = 5
  RNG_c(4) = 12345
  RNG_c(5) = 1
  RNG_c(6) = 2531011
  RNG_c(7) = 0
  ! seed RNG from system time, this might be required for the momentum-initialization 
  if (.not.MD_RNG_seeded) then
     ! we must have one and the same seed on all processes, therefore initialize on process 0 only.
     ! random_own is protected against calls on other processes anyway but one never knows.
     if (myid.eq.0) then
        ! initialize random number generator on myid=0
        call date_and_time(cdate, ctime)
        read(ctime,*) rng_seed
        seed = int(rng_seed)
     else
       ! other processes
       seed = 0 ! for sync via allreduce below
     end if
     call sync_integer(seed)
     write(info_str,'(2X,A,I24)') & 
     '| Initial seed for random number generator from system time: ', seed
     call localorb_info(info_str)
  end if
  if (seed.lt.0) seed = - seed
  if (seed.eq.0) seed = seed + 1  

  random_gauss_saved = .false.
  MB_sample_box_mueller_saved = .false.
  random_own_Y = 0
  gasdev_iset = 0
end subroutine initialize_RNG
!******

!****s* molecular_dynamics/gasdev
!  NAME
!    gasdev
!  SYNOPSIS
double precision function gasdev()
!  PURPOSE
!    extracts a random number from a gaussian, for BDP thermostat
!    it does more or less the same thing as MB_sample_box_mueller()
!  USES
!  AUTHOR
!    adapted from Davide Donadio
!  HISTORY
!    Release version, FHI-aims (2009).
!  INPUTS
!    none
!  OUTPUT
!    random number, normal distributed
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
  use runtime_choices
!  use random_module
  implicit none
  real*8 fac,rsq,v1,v2
  if(gasdev_iset==0) then
!     if(MD_random_BDP) then
1       v1=2.*random_own()-1.0d0
        v2=2.*random_own()-1.0d0
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
!      else
! 2       v1=2.*grnd()-1.0d0
!         v2=2.*grnd()-1.0d0
!         rsq=v1**2+v2**2
!         if(rsq.ge.1..or.rsq.eq.0.)goto 2
!      end if
     fac=sqrt(-2.*log(rsq)/rsq)
     gasdev_gset=v1*fac
     gasdev=v2*fac
     gasdev_iset=1
  else
     gasdev=gasdev_gset
     gasdev_iset=0
  end if
end function gasdev
!******

!****s* molecular_dynamics/gasdev
!  NAME
!    gasdev
!  SYNOPSIS
double precision function gamdev(ia)
!  PURPOSE
!    extracts a random number from a gamma distribution, for BDP thermostat
!  USES
!  AUTHOR
!    adapted from Davide Donadio
!  HISTORY
!    Release version, FHI-aims (2009).
!  INPUTS
!    # degrees of freedom
!  OUTPUT
!    random number gamma distributed
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
  use runtime_choices
!  use random_module
      implicit none
      integer, intent(in) :: ia
      integer j
      real*8 am,e,s,v1,v2,x,y
!      if(ia.lt.1)pause 'bad argument in gamdev'
      if(ia.lt.6)then
         x=1.
         do j=1,ia
            x=x*random_own()
         enddo
         x=-log(x)
      else
!          if(MD_random_BDP) then
1           v1=2.*random_own()-1.
            v2=2.*random_own()-1.
            if(v1**2+v2**2.gt.1.)goto 1
            y=v2/v1
            am=ia-1
            s=sqrt(2.*am+1.)
            x=s*y+am
            if(x.le.0.)goto 1
            e=(1.+y**2)*exp(am*log(x/am)-s*y)
!         if(grnd().gt.e)goto 1
            if(random_own().gt.e)goto 1
!          else
! 2           v1=2.*grnd()-1.
!             v2=2.*grnd()-1.
!             if(v1**2+v2**2.gt.1.)goto 2
!             y=v2/v1
!             am=ia-1
!             s=sqrt(2.*am+1.)
!             x=s*y+am
!             if(x.le.0.)goto 2
!             e=(1.+y**2)*exp(am*log(x/am)-s*y)
!             if(grnd().gt.e)goto 2
!          end if
      end if
      gamdev=x
end function gamdev
!******

!****s* molecular_dynamics/sumnoises
!  NAME
!    sumnoises
!  SYNOPSIS
double precision function sumnoises(nn)
!  PURPOSE
!    extracts a random number from a gamma distribution, for BDP thermostat
!  USES
!  AUTHOR
!    adapted from Davide Donadio
!  HISTORY
!    Release version, FHI-aims (2009).
!  INPUTS
!    # degrees of MD freedom -1
!  OUTPUT
!    random number gamma distributed
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
        implicit none
        integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
        !real*8, external :: gamdev,gasdev
        if(modulo(nn,2)==0) then
          sumnoises=2.0*gamdev(nn/2)
        else
          sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
        end if
end function sumnoises

!******


!****s* molecular_dynamics/resamplekin
!  NAME
!    resamplekin
!  SYNOPSIS
function resamplekin(kk,taut)
!  PURPOSE
!    gives the new kinetic energy (in unit of K)
!  USES
  use runtime_choices
!  AUTHOR
!    adapted from Davide Donadio
!  HISTORY
!    Release version, FHI-aims (2009).
!  INPUTS
!    current istantaneous temperature, relaxation time
!  OUTPUT
!    new istantaneous temperature
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
        implicit none
        real*8 :: resamplekin
        real*8, intent(in) :: kk,taut

        real*8 :: fact,rr
        !real*8, external :: sumnoises,gasdev
        
        if(taut>0.1) then
          fact=exp(-1.0/taut)
        else
          fact=0.0
        end if
        rr = gasdev()

        resamplekin = kk + (1.0-fact)* &
                       (MD_temperature*(sumnoises(MD_g_DOF-1)+rr**2)/MD_g_DOF-kk) &
                   + 2.0*rr*sqrt(kk*MD_temperature/MD_g_DOF*(1.0-fact)*fact)

end function resamplekin

!****s* molecular_dynamics/write_MD_restart
!  NAME
!    write_MD_restart
!  SYNOPSIS
subroutine write_MD_restart
!  PURPOSE
!    keep restart information (velocity, forces, thermostat parameters & such)
!    at the end of every time step; potentially this can be extended to the end of
!    every nth time step if really necessary 
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use timing
  use species_data
  use generate_aims_uuid, only: write_aims_uuid
  use localorb_io
  implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
! Here is the format for these restart files: 
!   (1) integer :: number of atoms   --  (merely for consistency checking
!   (2) real*8, dimension(3,n_atoms) :: current positions
!   (3) real*8, dimension(3,n_atoms) :: v_half
!   (4) real*8                       :: s
!   (5) real*8                       :: phalf
! this current format is the minimal thing that would work for all currently
! implemented thermostats while it still allows switching between various
! thermostats in different runs (or changing the external parameters for 
! example). If the need for any other numbers arises, please add them to THE
! END of this file and initialize s=1; phalf=0 unless you absolutely need these
! for your integrator.
  real*8 :: s_write, s_last_write
  real*8, dimension(3,n_atoms) :: v_half_write, v_last_write
  integer i_atom, i_s
  character(LEN=80) :: uuid_str
  if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
     do i_atom = 1, n_atoms
        v_half_write(:,i_atom) = v_half(:,i_atom)/species_m(species(i_atom))
	v_last_write(:,i_atom) = v_last(:,i_atom)/species_m(species(i_atom))
     end do
  else
     do i_atom = 1, n_atoms
        v_half_write(:,i_atom) = v_half(:,i_atom)
	v_last_write(:,i_atom) = v_last(:,i_atom)
     end do
  end if	!this does not work if inside myid .eq.0
 ! only write from the zero task ... 
  if (myid.eq.0) then
     s_write      = s_NP
     s_last_write = s_NP_last
     if (MD_ensemble.eq.'NVT_nose-hoover') then
        s_write      = dexp(s_write)
        s_last_write = dexp(s_last_write)
     end if

     call write_aims_uuid(uuid_str)
     open(file = MD_restart_file, unit = 88, status = 'replace', action='write')
     write(88,'(A)') '############### aims_MD_restart.dat #######################'
     write(88,'(A)') '# This file contains all information that is necessary to #'
     write(88,'(A)') '# continue a MD run. The data in this file should not be  #'
     write(88,'(A)') '# altered by the user, unless you REALLY know what you    #'
     write(88,'(A)') '# are doing.                                              #'
     write(88,'(A)') '###########################################################'
     write(88,'(A,2X,A)') '#', uuid_str
     write(88,'(A,I6)')        'n_atoms              ',n_atoms
     do i_atom=1, n_atoms
       write(88,'(A,3E30.20)') 'coords               ', coords(:,i_atom)
     end do
     do i_atom=1, n_atoms
       write(88,'(A,3E30.20)') 'v_half               ',v_half_write(:,i_atom)
     end do
     do i_atom=1, n_atoms
       write(88,'(A,3E30.20)') 'v_last               ',v_last_write(:,i_atom)
     end do
     do i_atom=1, n_atoms
       write(88,'(A,3E30.20)') 'r_last               ',r_last(:,i_atom)
     end do
     write(88,'(A,E30.20)')    's_write              ', s_write
     write(88,'(A,E30.20)')    's_NP_half            ', s_NP_half
     write(88,'(A,E30.20)')    's_last_write         ', s_last_write
     write(88,'(A,E30.20)')    's_dot_NP_half        ', s_dot_NP_half
     write(88,'(A,E30.20)')    's_dot_NP_last        ', s_dot_NP_last
     write(88,'(A,E30.20)')    'tsystem_half         ', tsystem
     write(88,'(A,E30.20)')    'tsystem_last         ', tsystem_last
     write(88,'(A,E30.20)')    'MD_H0                ', MD_H0
     write(88,'(A,E30.20)')    'MD_Epot_last         ', MD_Epot_last
     write(88,'(A,E30.20)')    'BDP_psuedo-Hamilt.   ', BDP_conint
     write(88,'(A,I6)')        'MD_stepcount         ', MD_stepcount
     write(88,'(A,I6)')        'MD_high_order_i      ', MD_high_order_i
     write(88,'(A,I6)')        'MD_force_evaluations ', MD_force_evaluations
     if (MD_ensemble.eq.'GLE_thermostat') then
        write(88,'(A,E30.20)') 'GLE_pseudo', langham_last
        write(88,'(A,I6)') 'MD_gle_ns', MD_gle_ns
        do i_atom=1, n_atoms
            do i_s=1, MD_gle_ns+1
                write(88,*) 'gp', gp(:,i_atom,i_s)
            enddo
        enddo
!        do i_atom=1, n_atoms
!            do i_s=1, MD_gle_ns+1
!                write(88,*) 'ngp', ngp(:,i_atom,i_s)
!            enddo
!        enddo
        do i_s=1, MD_gle_ns+1
                write(88,'(A,20E30.20)') 'gS', gS(i_s,:)
        enddo
        do i_s=1, MD_gle_ns+1
                write(88,'(A,20E30.20)') 'gT', gT(i_s,:)
        enddo 
     endif
     close(unit=88)
  end if
end subroutine write_MD_restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! new !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GLE thermostat routines
! *********************************************************
! * code is licensed under GPLv3 [www.gnu.org]            *
! * please consider citing the relevant papers (listed    *
! * below) if you use GLE in your simulations.            *
! *                                                       *
! * e-mail me at michele dot ceriotti at gmail dot com    *
! *********************************************************
! Routines slightly modified by MR for compatibility with FHI-aims -- 2012 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MC simple & clean gaussian PRNG

real*8 function random_gauss()
  use constants
  implicit none
  if (.not.random_gauss_saved) then
     random_gauss_saved = .true.
     random_gauss_r1    = random_own()
     random_gauss_r2    = random_own()
     random_gauss     = sqrt(-2d0*dlog(random_gauss_r1))*cos(2d0*pi*random_gauss_r2)
  else
     random_gauss_saved = .false.
     random_gauss     = sqrt(-2d0*dlog(random_gauss_r1))*sin(2d0*pi*random_gauss_r2)
  end if
end function random_gauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_gle(MD_temperature, MD_gle_ns)
   use constants, only: boltzmann_kB, MD_KE_factor 
   use dimensions
   use localorb_io
   use mpi_tasks
   implicit none
!   real *8, allocatable :: A_mat(:,:), C_mat(:,:)
!   real*8, intent(in):: wopt
    real*8, intent(in) :: MD_temperature
    integer, intent(in) :: MD_gle_ns
    character*160 :: info_str
    character*40 :: desc_str
    character*400 :: inputline
    character(*), parameter :: func = 'read_gle'
    integer :: i_code, ia, ic, i

    i_code = 0
    ia = 0
    ic = 0

   !allocate A and C matrices, which need to be inputed for the thermostat
    allocate(gA(MD_gle_ns+1,MD_gle_ns+1))
    allocate(gC(MD_gle_ns+1,MD_gle_ns+1))
    
    gA(:,:)=0.d0
    gC(:,:)=0.d0
    
    ! allocation check
    if ((.not.allocated(gA)).and.(.not.allocated(gC))) then
      if (myid.eq.0) then
         write(0,*) "** Matrices not correctly allocated/initialized, bye bye :("
      end if
      stop
    end if
   
    lineloop: do
        read(7,'(A)',iostat=i_code) inputline
        if(i_code<0) then ! EOF
            backspace(7)
            exit lineloop
        elseif(i_code>0) then
            ! handle error at main level
            backspace(7)
            return
        endif
!       if (verbatim) call localorb_info(inputline,use_unit,'(2X,A)')
        read(inputline,*,iostat=i_code) desc_str
        if (i_code.ne.0) then
            cycle  ! empty line
        elseif (desc_str(1:1).eq."#") then
            cycle  ! comment
        elseif (desc_str.eq."MD_gle_A") then
        !reads A (in units of 1/ps)
            ia=ia+1
! MR: here must introduce a error message if ia > MD_gle_ns+1
            read(inputline,*,end=88,err=99) desc_str, gA(ia,:)
            ! DEBUG
            ! write(use_unit,*) 'My A matrix', ia, gA(ia,:)
            ! END DEBUG
        elseif (desc_str.eq."MD_gle_C") then
! MR: here must introduce a error message if ia > MD_gle_ns+1
            ic=ic+1
            read(inputline,*,end=88,err=99) desc_str, gC(ic,:)
            ! DEBUG
            ! write(use_unit,*) 'My C matrix', ic, gC(ic,:)
            ! END DEBUG
        else
!          must have reached end of matrices
            backspace(7)
            exit lineloop
        endif
    enddo lineloop

!   now verify input data   

    if (ia.eq.0) then
       if (myid .eq. 0) then
          write(info_str,'(X,A)') "** Error: could not read any A matrix!"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(X,A)') "** Please specify at this input matrix for this thermostat"
          call localorb_info(info_str,use_unit,'(A)')
       end if
       call aims_stop_coll('', func)
    end if

    if (ic.eq.0) then
        if (myid .eq. 0) then
          write(info_str,'(X,A)') "| No C matrix was specified"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(X,A)') "| Using canonical-sampling, C_p=kT"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(X,A)') "| Running molecular dynamics in the NVT ensemble"
          call localorb_info(info_str,use_unit,'(A)')
       end if
       gC=0.
       do i=1,MD_gle_ns+1
          gC(i,i)=boltzmann_kB*MD_temperature/(MD_KE_factor) ! this factor transforms hartree to internal units
       enddo
    else
       if (myid .eq. 0) then
          write(info_str,'(X,A)') "| C matrix was specified"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(X,A)') "| Reading specialized Cp matrix"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(X,A)') "| Detailed balance is (probably) broken"
          call localorb_info(info_str,use_unit,'(A)')
       end if
    ! Perform unit conversion of the C matrix to internal units
       gC(:,:)=gC(:,:)*boltzmann_kB/MD_KE_factor
    endif
   
   
!reads in matrices

!    open(121,file='GLE-A',status='OLD',iostat=ios)
!    if (ios.ne.0) then
!       if (myid .eq. 0) then
!          write(info_str,'(X,A)') "** Error: could not read GLE-A file!"
!          call localorb_info(info_str,use_unit,'(A)')
!          write(info_str,'(X,A)') "** Please specify at least 'GLE-A' file for this thermostat"
!          call localorb_info(info_str,use_unit,'(A)')
!       end if
!       call aims_stop_coll('', func)
!    end if



! read matrix from separate files for now
!   do i=1,MD_gle_ns+1
!       read(121,*) gA(i,:)
!    enddo
!    close(121)

    ! reads C (in K), or init to kT
!    open(121,file='GLE-C',status='OLD',iostat=ios)
!    if (ios.ne.0) then
!       if (myid .eq. 0) then
!          write(info_str,'(X,A)') "| No C matrix was specified"
!          call localorb_info(info_str,use_unit,'(A)')
!          write(info_str,'(X,A)') "| Using canonical-sampling, C_p=kT"
!          call localorb_info(info_str,use_unit,'(A)')
!          write(info_str,'(X,A)') "| Running molecular dynamics in the NVT ensemble"
!          call localorb_info(info_str,use_unit,'(A)')
!       end if
!       gC=0.
!       do i=1,MD_gle_ns+1
!          gC(i,i)=boltzmann_kB*MD_temperature/(MD_KE_factor) ! this factor transforms hartree to internal units
!       enddo
!    else
!       if (myid .eq. 0) then
!          write(info_str,'(X,A)') "| C matrix was specified"
!          call localorb_info(info_str,use_unit,'(A)')
!          write(info_str,'(X,A)') "| Reading specialized Cp matrix"
!          call localorb_info(info_str,use_unit,'(A)')
!          write(info_str,'(X,A)') "| Detailed balance is (probably) broken"
!          call localorb_info(info_str,use_unit,'(A)')
!       end if
!
!       do i=1,MD_gle_ns+1
!          read(121,*) gC(i,:)
!       enddo
!       ! Perform unit conversion of the C matrix to internal units
!       gC(:,:)=gC(:,:)*boltzmann_kB/MD_KE_factor
!    end if
!    close(121)

    return

    88 continue
        if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in' (missing arguments)"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
        end if
        call aims_stop_coll('', func)

    99 continue
        if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'"
            write(use_unit,*) "line: '"//trim(inputline)//"'"
        end if
        call aims_stop_coll('', func)

 end subroutine read_gle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! initialize gle_init
  subroutine gle_init(dt)
    use synchronize_mpi
    use mpi_tasks
    use runtime_choices
    use dimensions
    use localorb_io
    use species_data
    implicit none
    real*8, intent(in)  :: dt!, wopt
!    integer, intent(inout) :: irnd
    real *8, allocatable :: gr(:)
    character*160 :: info_str
    integer i, j, k, h, i_atom, i_coord
    
    if (myid.eq.0) then
      write(info_str,'(X,A)') " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(X,A)') " | Initialization of GLE thermostat."
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(X,A)')" | Please cite the relevant works among:"
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(X,2A)')" | M. Ceriotti, G. Bussi and M. Parrinello ",&
      "Phy. Rev. Lett. 102, 020601 (2009)"
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(X,2A)')" | M. Ceriotti, G. Bussi and M. Parrinello ",&
      "Phy. Rev. Lett. 103, 030603 (2009)"
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(X,2A)')" | M. Ceriotti, G. Bussi and M. Parrinello ",&
      "J. Chem. Theory Comput. 6, 1170 (2010)"
      call localorb_info(info_str,use_unit,'(A)')
    end if
 !! perform allocation checks here !!

    !allocate all extra matrices everything we need
    call allocate_gle() !allocates general propagator quantities that will be needed until the end of the simulation
    allocate(gr(MD_gle_ns+1)) ! allocates vector with random numbers
    

    ! the deterministic part of the propagator is obtained in a second
    ! T matrix from equation after 2.27 in M. Ceriotti's thesis
    call matrix_exp(-dt*gA, MD_gle_ns+1,15,15,gT)
    ! the stochastic part is just as easy. we use gA as a temporary array
    ! Also from Cholesky decomposition of eq. after 2.27 in M. Ceriotti's thesis
    gA=gC-matmul(gT,matmul(gC,transpose(gT)))
    call cholesky_own(gA, gS, MD_gle_ns+1)


    ! then, we must initialize the auxiliary vectors. we keep general - as we might be 
    ! using non-diagonal C to break detailed balance - and we use cholesky decomposition
    ! of C. again, since one would like to initialize correctly the velocities in 
    ! case of generic C, we use an extra slot for gp for the physical momentum, as we 
    ! could then use it to initialize the momentum in the calling code
    gA=gC   
    call cholesky_own(gA, gC, MD_gle_ns+1)
! MR:from now on gC is already "square rooted"

! This part has to be done only in thread 0, otherwise there is too much randomness involved

   if (myid.eq.0) then
     do i_coord=1, 3    
       do i_atom=1,n_atoms
         do i=1,MD_gle_ns+1
           gr(i)=random_gauss()
         enddo
         gp(i_coord,i_atom,:)=matmul(gC,gr)
       end do
     end do
   endif

! MR: Should I take the square root to maintain momenta dimensions? 
   do i=1,MD_gle_ns+1
     call broadcast_MD_velocities(gp(:,:,i),0)
   end do

   deallocate(gr)
   
   langham=0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
   langham_last=0.d0

  end subroutine gle_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! the GLE propagator. 
  ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
  ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
  ! gp(j,2:ns+1) the current values of additional momenta. 
  ! the matrix multiplies are performed on the whole array at once, and the new momentum
  ! is passed back to the caller, while the new s's are kept stored in gp.
  ! please note that one can avoid the double conversion between mass-scaled and actual
  ! momentum/velocity (as used by the calling code) by scaling with the mass the white-noise
  ! random numbers. 
  subroutine gle_step(velocities)
    use mpi_tasks
    use synchronize_mpi
    use runtime_choices
    use dimensions
    use species_data
    use geometry, only: species
    implicit none
    real *8, intent(inout) :: velocities(:,:)
    integer i, j, k, i_atom, i_dim, ndim
!    real*8 mfac, totm

    ndim=n_atoms*3
    do i_dim=1, 3
      do i_atom=1,n_atoms
        gp(i_dim,i_atom,1)=dsqrt(species_m(species(i_atom)))*velocities(i_dim,i_atom)   ! go to mass-scaled momenta and load into gp array
      enddo
    enddo

! Calculates deterministic part of the propagator T.(p,s), from eq. 2.27 of M. Ceriotti's thesis
! ... yes, but for each _dimension_ separately? must write down
!    do i_atom=1, n_atoms
!    do i_dim=1, 3
!      call dgemm('n','t', 3,ns+1,ns+1,1.0d0,gp(:,i_atom,:),3,gT,ns+1,0.0d0,ngp(:,i_atom,:),3)
!      call dgemm('n','t',n_atoms,ns+1,ns+1,1.0d0,gp(i_dim,:,:),n_atoms,gT,ns+1,0.0d0,ngp(i_dim,:,:),n_atoms)
!    end do
!    end do

    ! we have 1 possibility in 8 to have done this right
    call dgemm('n','t',ndim,MD_gle_ns+1,MD_gle_ns+1,1.0d0,gp,ndim,gT,MD_gle_ns+1,0.0d0,ngp,ndim)
!    call dgemm('n','n',ns+1,ndim,ns+1,1.0d0,gT,ns+1,gp,ns+1,0.0d0,ngp,ndim)
!    ngp=transpose(matmul(gT,transpose(gp))

    !now, must compute random part. 
    !first, fill up gp of random n
    !MR: but only in thread 0, otherwise it gets alllll messed up
    if (myid .eq. 0) then
      do i_dim=1, 3
        do i_atom=1,n_atoms
          do i=1,MD_gle_ns+1
            gp(i_dim,i_atom,i)=random_gauss() 
          end do
        end do
      end do
    end if

    do i=1,MD_gle_ns+1
      call broadcast_MD_velocities(gp(:,:,i),0)
    end do
 
    call dgemm('n','t',ndim,MD_gle_ns+1,MD_gle_ns+1,1.0d0,gp,ndim,gS,MD_gle_ns+1,1.0d0,ngp,ndim)
    gp = ngp
!    gp=transpose(matmul(gS,transpose(gp))+ngp
 

      
    do i_dim=1, 3
      do i_atom=1,n_atoms
        velocities(i_dim,i_atom)=gp(i_dim,i_atom,1)/dsqrt(species_m(species(i_atom)))
      end do
    end do


    
  end subroutine gle_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! matrix exponential by scale & square.
  ! one can diagonalize with lapack, but it's not worth it, as 
  ! we call this routine only once!      
  subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    real*8, intent(in)   :: M(n,n)
    real*8, intent(out)   :: EM(n,n)
    real *8 :: tc(j+1), tmp(n,n), SM(n,n)
    integer p, i
    
    tc(1)=1
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo
    
    !scale
    SM=M*(1./2.**k)
    EM=0.
    do i=1,n
       EM(i,i)=tc(j+1)
    enddo
    
    !taylor exp of scaled matrix
    do p=j,1,-1
       EM=matmul(SM,EM);
       do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       enddo
    enddo
    
    !square
    do p=1,k
       EM=matmul(EM,EM)
    enddo
  end subroutine matrix_exp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! brute-force "stabilized" cholesky decomposition.
  ! in practice, we compute LDL^T decomposition, and force
  ! to zero negative eigenvalues.
 subroutine cholesky_own(SST, S, n)
    integer, intent(in)  :: n
    real*8, intent(in)   :: SST(n,n)
    real*8, intent(out)   :: S(n,n)
    real *8 :: L(n,n), D(n,n) 
    integer i,j,k
    S=0.
    L=0.
    D=0.
    do i=1,n
       L(i,i)=1.0
       do j=1,i-1
          L(i,j)=SST(i,j);
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k,k)
          enddo
          if (D(j,j).ne. 0.0) then
            L(i,j)=L(i,j)/D(j,j) 
          else
            write(use_unit,*) "Warning: zero eigenvalue in LDL^T decomposition."
            L(i,j)=0.
          end if
       enddo
       D(i,i)=SST(i,i)
       do k=1,i-1
          D(i,i)=D(i,i)-L(i,k)**2*D(k,k)
       end do
    enddo
    do i=1,n
       if ( D(i,i).ge. 0.0d0 ) then
         D(i,i)=sqrt(D(i,i))
       else
         write(use_unit,*) "Warning: negative eigenvalue (",D(i,i),")in LDL^T decomposition."
         D(i,i)=0.0
       end if
    end do
    S=matmul(L,D)
end subroutine cholesky_own
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine allocate_gle()

    use runtime_choices
    use dimensions
    
    if(.not.allocated(gS)) then
      allocate(gS(MD_gle_ns+1,MD_gle_ns+1))
      gS(:,:)=0.0
    end if
    if(.not.allocated(gT)) then
      allocate(gT(MD_gle_ns+1,MD_gle_ns+1))
      gT(:,:)=0.0
    end if
    if(.not.allocated(gp)) then
      allocate(gp(3,n_atoms,MD_gle_ns+1))
      gp(:,:,:)=0.0
    end if
    if(.not.allocated(ngp)) then
      allocate(ngp(3,n_atoms,MD_gle_ns+1))
      ngp(:,:,:)=0.0
    end if
    
end subroutine allocate_gle


subroutine clean_gle()
  
    if (allocated(gS      )) deallocate(gS      )
    if (allocated(gT      )) deallocate(gT      )
    if (allocated(gp      )) deallocate(gp      )
    if (allocated(ngp     )) deallocate(ngp      )
    ! deallocate temporary initialization matrices   
    if (allocated(gA     )) deallocate(gA)
    if (allocated(gC     )) deallocate(gC)
 
end subroutine clean_gle
  


!******
end module molecular_dynamics
