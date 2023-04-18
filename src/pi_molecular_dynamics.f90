!****h* FHI-aims/pi_molecular_dyn
!  NAME
!    pi_molecular_dynamics - PIMD capacities of FHI-aims
!  SYNOPSIS

module pi_molecular_dynamics

!  PURPOSE
!    This module takes care of everything related to PIMD and requires only energies and forces
!  AUTHOR
!    Xin-Zheng Li, based on the module molecular_dynamics from the FHI-aims team
!  HISTORY
!    Release version, FHI-aims (2012).
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
!    Release version, FHI-aims (2012).
!  SOURCE
  implicit none

  ! thermostat internal variables and starting values
  real*8 :: s_NP_beads          = 1d0
  real*8 :: s_NP_beads_last     = 1d0
  real*8 :: s_NP_beads_half     = 1d0
  real*8 :: s_dot_NP_beads_half = 0d0
  real*8 :: s_dot_NP_beads_last = 0d0
  real*8 :: tsystem_beads       = 0d0
  real*8 :: tsystem_beads_last  = 0d0

  real*8, dimension(:,:,:), allocatable :: v_beads_half, v_beads_last, r_beads_last!, v_comparison
  real*8, dimension(:,:,:), allocatable :: v_beads_nm_half, v_beads_nm_last, r_beads_nm_last!, v_comparison
  real*8, dimension(:,:,:), allocatable :: v_beads_nm_half_orig, v_beads_nm_last_orig, v_beads_orig

  real*8, dimension(:,:), allocatable :: vb_half, vb_last, rb_last
  real*8, dimension(:), allocatable :: free_energy_beads, lamda_nm

  real*8, dimension(:,:,:), allocatable :: total_forces_beads
  real*8, dimension(:,:,:), allocatable :: total_forces_beads_nm
  real*8, dimension(:,:,:), allocatable :: forces_lv_beads 

  real*8, dimension(:,:,:), allocatable :: spring_forces_beads
  real*8, dimension(:,:,:), allocatable :: spring_forces_beads_nm

  real*8, dimension(:,:,:), allocatable :: coords_beads  
  real*8, dimension(:,:,:), allocatable :: frac_coords_beads
  real*8, dimension(:,:,:), allocatable :: coords_beads_nm
  real*8, dimension(:,:,:), allocatable :: frac_coords_beads_nm

  real*8, dimension(:,:), allocatable :: transform_u
  real*8, dimension(:,:), allocatable :: matrix_a
  real*8, dimension(:,:), allocatable :: matrix_lamda
  real*8, dimension(:,:), allocatable :: matrix_au

  integer*8, dimension(:,:), allocatable :: counter_beads
  integer*8, dimension(7), private :: RNGB_a, RNGB_c, RNGB_m

  logical :: use_PIMD_max_steps = .false.
  logical :: PIMD_RNG_firstcall = .true.
  real*8  :: PIMD_H0        = 0d0
  real*8  :: PIMD_Epot_last = 0d0
  logical :: PIMD_successful_restart_read = .false.
!  integer :: PIMD_schedule_step
  integer :: n_atoms_PIMD 
!  logical, dimension(:), allocatable :: constrain_PIMD

contains

!******	
!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/allocate_PIMD
!  NAME
!    allocate_PIMD
!  SYNOPSIS

subroutine allocate_PIMD

!  PURPOSE
!    allocation of PIMD
!  USES

  use dimensions
  use runtime_choices
  use geometry
  use physics

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2012).
!  SOURCE
  implicit none

  integer :: p

  if (.not.allocated(r_beads_last)) then
    allocate(r_beads_last(3,n_beads,n_atoms))
    r_beads_last(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_half)) then
    allocate(v_beads_half(3,n_beads,n_atoms))
    v_beads_half(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_last)) then
    allocate(v_beads_last(3,n_beads,n_atoms))
    v_beads_last(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_orig)) then
    allocate(v_beads_orig(3,n_beads,n_atoms))
    v_beads_orig(:,:,:) = 0d0
  end if 

  if (.not.allocated(r_beads_nm_last)) then
    allocate(r_beads_nm_last(3,n_beads,n_atoms))
    r_beads_nm_last(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_nm_half)) then
    allocate(v_beads_nm_half(3,n_beads,n_atoms))
    v_beads_nm_half(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_nm_last)) then
    allocate(v_beads_nm_last(3,n_beads,n_atoms))
    v_beads_nm_last(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_nm_half_orig)) then
    allocate(v_beads_nm_half_orig(3,n_beads,n_atoms))
    v_beads_nm_half_orig(:,:,:) = 0d0
  end if
  if (.not.allocated(v_beads_nm_last_orig)) then
    allocate(v_beads_nm_last_orig(3,n_beads,n_atoms))
    v_beads_nm_last_orig(:,:,:) = 0d0
  end if
  if (.not.allocated(rb_last)) then
    allocate(rb_last(3,n_atoms))
    rb_last(:,:) = 0d0
  end if
  if (.not.allocated(vb_half)) then
    allocate(vb_half(3,n_atoms))
    vb_half(:,:) = 0d0
  end if
  if (.not.allocated(vb_last)) then
    allocate(vb_last(3,n_atoms))
    vb_last(:,:) = 0d0
  end if
  if (.not.allocated(counter_beads)) then
    allocate(counter_beads(n_beads,n_atoms))
    counter_beads(:,:) = 0
  end if

  if (.not.allocated(transform_u)) then
    allocate(transform_u(n_beads,n_beads))
    transform_u(:,:) = 0
    allocate(matrix_a(n_beads,n_beads))
    allocate(matrix_lamda(n_beads,n_beads))
    allocate(matrix_au(n_beads,n_beads))
    matrix_a = 0.0d0
    matrix_lamda = 0.0d0
    matrix_au = 0.0d0
  end if

  if(.not.allocated(lamda_nm)) then
    allocate(lamda_nm(n_beads))
    lamda_nm(:) = 0.0d0
  end if

!  if (.not.allocated(constrain_PIMD)) then
!     allocate(constrain_PIMD(n_atoms))
!     constrain_PIMD(:) = .false.
!  end if
  if (.not.allocated(free_energy_beads)) then
    allocate(free_energy_beads(n_beads))
    free_energy_beads(:) = 0d0
  end if

! allocate the spatial coordinates for PIMD

  if (.not.allocated(coords_beads)) then
    allocate(coords_beads(3,n_beads,n_atoms))
    coords_beads(:,:,:) = 0d0
  end if

  if (.not.allocated(frac_coords_beads)) then
    allocate(frac_coords_beads(3,n_beads,n_atoms))
    frac_coords_beads = 0d0
  end if

  if (.not.allocated(coords_beads_nm)) then
    allocate(coords_beads_nm(3,n_beads,n_atoms))
    coords_beads_nm(:,:,:) = 0d0
  end if

  if (.not.allocated(frac_coords_beads_nm)) then
    allocate(frac_coords_beads_nm(3,n_beads,n_atoms))
    frac_coords_beads_nm = 0d0
  end if

!  allocate the forces for PIMD

  if (.not.allocated(total_forces_beads)) then
    allocate(total_forces_beads(3,n_beads,n_atoms))
    total_forces_beads = 0.d0
  end if

  if (.not.allocated(total_forces_beads_nm)) then
    allocate(total_forces_beads_nm(3,n_beads,n_atoms))
    total_forces_beads_nm = 0.d0
  end if

  if (.not.allocated(forces_lv_beads)) then
    allocate(forces_lv_beads(3,n_beads,n_atoms))
  end if

  if (.not.allocated(spring_forces_beads)) then
    allocate(spring_forces_beads(3,n_beads,n_atoms))
    spring_forces_beads = 0d0
  end if

  if (.not.allocated(spring_forces_beads_nm)) then
    allocate(spring_forces_beads_nm(3,n_beads,n_atoms))
    spring_forces_beads_nm = 0d0
  end if

  PIMD_g_DOF = 3*n_atoms

!! allocate the various settings for MD schedule
!  if (PIMD_use_schedule) then
!     if (.not.allocated(PIMD_schedule_ensemble      )) allocate(PIMD_schedule_ensemble      (PIMD_segments))
!     if (.not.allocated(PIMD_schedule_temperature   )) allocate(PIMD_schedule_temperature   (PIMD_segments))
!     if (.not.allocated(PIMD_schedule_time          )) allocate(PIMD_schedule_time          (PIMD_segments))
!     if (.not.allocated(PIMD_schedule_Q             )) allocate(PIMD_schedule_Q             (PIMD_segments))
!     ! initialize to some neutral value
!     PIMD_schedule_ensemble      (:) = "NVT"
!     PIMD_schedule_temperature   (:) = 0d0
!     PIMD_schedule_time          (:) = 0d0
!     PIMD_schedule_Q             (:) = 1d0
!  end if
!  ! set degrees of freedom for thermostats 
end subroutine allocate_PIMD
!******

!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/clean_PIMD
!  NAME
!    clean_PIMD - deallocation of path-integral molecular dynamics
!  SYNOPSIS
subroutine clean_PIMD
!  PURPOSE
!    deallocation of local variables
!  USES
   use runtime_choices
   use geometry
   use physics
!  AUTHOR
!    FHI-aims
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(r_beads_last))         deallocate(r_beads_last)
  if (allocated(v_beads_half))         deallocate(v_beads_half)
  if (allocated(v_beads_last))         deallocate(v_beads_last)
  if (allocated(v_beads_orig))         deallocate(v_beads_orig)
  if (allocated(r_beads_nm_last))      deallocate(r_beads_nm_last)
  if (allocated(v_beads_nm_half))      deallocate(v_beads_nm_half)
  if (allocated(v_beads_nm_last))      deallocate(v_beads_nm_last)
  if (allocated(v_beads_nm_half_orig)) deallocate(v_beads_nm_half_orig)
  if (allocated(v_beads_nm_last_orig)) deallocate(v_beads_nm_last_orig)
  if (allocated(rb_last))              deallocate(rb_last)
  if (allocated(vb_half))              deallocate(vb_half)
  if (allocated(vb_last))              deallocate(vb_last)
  if (allocated(counter_beads))        deallocate(counter_beads)
  if (allocated(coords_beads))         deallocate(coords_beads)
  if (allocated(frac_coords_beads))    deallocate(frac_coords_beads)
  if (allocated(total_forces_beads))   deallocate(total_forces_beads)
  if (allocated(total_forces_beads_nm)) deallocate(total_forces_beads_nm)
  if (allocated(forces_lv_beads))      deallocate(forces_lv_beads)
  if (allocated(free_energy_beads))    deallocate(free_energy_beads)
  if (allocated(spring_forces_beads))  deallocate(spring_forces_beads)
!  if (allocated(PIMD_schedule_ensemble      )) deallocate(PIMD_schedule_ensemble      )
!  if (allocated(PIMD_schedule_temperature   )) deallocate(PIMD_schedule_temperature   )
!  if (allocated(PIMD_schedule_time          )) deallocate(PIMD_schedule_time          )
!  if (allocated(PIMD_schedule_Q             )) deallocate(PIMD_schedule_Q             )
!  if (allocated(constrain_PIMD              )) deallocate(constrain_PIMD              )

end subroutine clean_PIMD
!******
!!subroutine change_PIMD_schedule_step
!  use runtime_choices
!  use dimensions
!  use species_data
!  use geometry
!  use localorb_io
!  implicit none
!  character*120 :: info_str
!  integer :: i_atom
!
!  ! determine if one of the steps in the molecular dynamics schedule is over and we need to switch to the next one ... 
!  if ((tsystem_beads.ge.PIMD_time).and.(PIMD_use_schedule).and.(PIMD_schedule_step.lt.PIMD_segments)) then
!     ! set all PIMD variables to the current step
!     PIMD_schedule_step   = PIMD_schedule_step + 1 
!     PIMD_ensemble        = PIMD_schedule_ensemble      (PIMD_schedule_step)
!     ! CC: Ensure that a consistent number of steps is perfomed in each segment
!     !     MD_time            = MD_schedule_time          (MD_schedule_step) + MD_time 
!     PIMD_time            = PIMD_schedule_time          (PIMD_schedule_step) + PIMD_time + PIMD_tstep*0.99999999d0
!     PIMD_temperature     = PIMD_schedule_temperature   (PIMD_schedule_step)
!     PIMD_Q_NP            = PIMD_schedule_Q             (PIMD_schedule_step)
!  end if
!end subroutine change_PIMD_schedule_step
!
!****s* pi_molecular_dynamics/initialize_PIMD
!  NAME
!    initialize_PIMD
!  SYNOPSIS
subroutine initialize_PIMD(Epot)
!  PURPOSE
!    initialization routine and first PIMD step
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

  integer :: i_atom, i_bead, i_coord, j_bead, p, nm 
  real*8  :: T_inst, Epot, Ekin, shat, factor, w0, w1, w2, w3, vel
  real*8, dimension(3,n_atoms) :: BeadMD_forces
  real*8, dimension(n_atoms)   :: masses
  logical :: read_restart_succes
  character*100 :: info_str

  call initialize_RNG_pi ()

  call initialize_nm ()

 ! initialize MD schedule if present: simply set all the necessary variables for the first step in the schedule
!  if (PIMD_use_schedule) then
!     ! FIXME: Intelligent output about the initialization!!
!     PIMD_schedule_step   = 1 
!     PIMD_ensemble        = PIMD_schedule_ensemble      (PIMD_schedule_step)
!     PIMD_time            = PIMD_schedule_time          (PIMD_schedule_step) 
!     PIMD_temperature     = PIMD_schedule_temperature   (PIMD_schedule_step)
!     PIMD_Q_NP            = PIMD_schedule_Q             (PIMD_schedule_step)
!     !MD_random_BDP      = MD_schedule_random_BDP    (MD_schedule_step)
!  end if

   PIMD_init_temperature = PIMD_temperature
  ! find out number of mobile atoms, if some are restricted 
!  if (use_relaxation_constraints) then 
!     n_atoms_PIMD = 0
!     constrain_PIMD(:) = constrain_relaxation(:)
!     do i_atom = 1, n_atoms
!        if (.not.constrain_PIMD(i_atom)) then
!           n_atoms_PIMD = n_atoms_PIMD + 1 
!        end if        
!     end do
!     ! reset number of degrees of freedom
!     PIMD_g_DOF = 3*n_atoms_PIMD
!  else
   n_atoms_PIMD = n_atoms
   PIMD_g_DOF = 3*n_atoms_PIMD
!  end if

  if (use_nm) then

!  when nm.eq.1, do propagation in the normal-mode coordinates

    coords_beads_nm(1:3,1:n_beads,1:n_atoms) = 0.0d0

    if(.not.PIMD_initialconf_from_restart) then
      do i_atom = 1, n_atoms, 1
        do i_bead = 1, n_beads, 1
          do j_bead = 1, n_beads, 1
            do i_coord = 1, 3, 1
              coords_beads_nm(i_coord,i_bead,i_atom) = coords_beads_nm(i_coord,i_bead,i_atom) + 1/sqrt(1.0d0*n_beads) * &
   &                                                   transform_u(i_bead,j_bead) * coords_beads(i_coord,j_bead,i_atom)
            end do
          end do
        end do
      end do
      r_beads_nm_last(:,:,:) = coords_beads_nm(:,:,:)
    else
      r_beads_nm_last(1:3,1:n_beads,1:n_atoms) = 0.0d0
      do i_atom = 1, n_atoms, 1
        do i_bead = 1, n_beads, 1
          do j_bead = 1, n_beads, 1
            do i_coord = 1, 3, 1
              coords_beads_nm(i_coord,i_bead,i_atom) = coords_beads_nm(i_coord,i_bead,i_atom) + 1/sqrt(1.0d0*n_beads) * &
   &                                                   transform_u(i_bead,j_bead) * coords_beads(i_coord,j_bead,i_atom)
              r_beads_nm_last(i_coord,i_bead,i_atom) = r_beads_nm_last(i_coord,i_bead,i_atom) + 1/sqrt(1.0d0*n_beads) * &
   &                                                   transform_u(i_bead,j_bead) * r_beads_last(i_coord,j_bead,i_atom)
            end do
          end do
        end do
      end do
    end if

    if(.not.PIMD_initialconf_from_restart) then
      v_beads_nm_half(1:3,1:n_beads,1:n_atoms) = 0.0d0
      !------------------------------------------------------------------------------------------------------
      ! initialize velocities from MB distribution if requested
      MB_velocity_initialization = .true.
      if (MB_velocity_initialization) then
        do i_atom = 1, n_atoms, 1
          do i_bead = 1, n_beads, 1
            do i_coord = 1, 3, 1
              if(use_pimd) then
                if(i_bead.ne.1) then
                  call MB_sample_box_mueller_pi(PIMD_init_temperature,species_m(species(i_atom))*lamda_nm(i_bead),&
   &                                            v_beads_nm_half(i_coord,i_bead,i_atom))
                else
                  call MB_sample_box_mueller_pi(PIMD_init_temperature,species_m(species(i_atom)),&
   &                                            v_beads_nm_half(i_coord,i_bead,i_atom))
                end if
              else if(use_cmd) then
                if(i_bead.ne.1) then
                  call MB_sample_box_mueller_pi(PIMD_init_temperature,species_m(species(i_atom))*lamda_nm(i_bead)/6.4d+1,&
   &                                            v_beads_nm_half(i_coord,i_bead,i_atom))
                else
                  call MB_sample_box_mueller_pi(PIMD_init_temperature,species_m(species(i_atom)),&
   &                                            v_beads_nm_half(i_coord,i_bead,i_atom))
                end if
              else if(use_rpmd) then
                call MB_sample_box_mueller_pi(PIMD_init_temperature,species_m(species(i_atom)),&
   &                                         v_beads_nm_half(i_coord,i_bead,i_atom))
              else
                write(stderr,*) 'initialize_PIMD error, please use PIMD.or.CMD.or.RPMD'
                stop
              end if
            end do
          end do
        enddo
      end if
      v_beads_nm_last(:,:,:) = v_beads_nm_half(:,:,:)
    else 
      v_beads_nm_half(1:3,1:n_beads,1:n_atoms) = 0.0d0
      v_beads_nm_last(1:3,1:n_beads,1:n_atoms) = 0.0d0
      do i_atom = 1, n_atoms, 1
        do i_bead = 1, n_beads, 1
          do j_bead = 1, n_beads, 1
            do i_coord = 1, 3, 1
              v_beads_nm_half(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom) + 1/sqrt(1.0d0*n_beads) * &
   &                                                 transform_u(i_bead,j_bead) * v_beads_half(i_coord,j_bead,i_atom)
              v_beads_nm_last(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom) + 1/sqrt(1.0d0*n_beads) * &
   &                                                 transform_u(i_bead,j_bead) * v_beads_last(i_coord,j_bead,i_atom)
            end do
          end do
        end do
      end do
    end if

  ! Nosé Thermostats need to work with generalized momenta - must change from velocities!!!

    if (PIMD_ensemble.eq.'NVT_nose-hoover') then
       s_NP_beads_last = 0d0
       s_NP_beads_half = 0d0
       s_NP_beads      = 0d0
       do i_atom = 1, n_atoms, 1
         do i_bead = 1, n_beads, 1
           do i_coord = 1, 3, 1
             if(use_pimd) then
               if(i_bead.eq.1) then
                 v_beads_nm_half(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * &
    &                                                 v_beads_nm_half(i_coord,i_bead,i_atom)
                 v_beads_nm_last(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * &
    &                                                 v_beads_nm_last(i_coord,i_bead,i_atom)
               else
                 v_beads_nm_half(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * lamda_nm(i_bead) * &
    &                                                 v_beads_nm_half(i_coord,i_bead,i_atom)
                 v_beads_nm_last(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * lamda_nm(i_bead) * &
    &                                                 v_beads_nm_last(i_coord,i_bead,i_atom)
               end if
             else if(use_cmd) then
               if(i_bead.eq.1) then
                 v_beads_nm_half(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * &
    &                                                 v_beads_nm_half(i_coord,i_bead,i_atom)
                 v_beads_nm_last(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * &
    &                                                 v_beads_nm_last(i_coord,i_bead,i_atom)
               else
                 v_beads_nm_half(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * lamda_nm(i_bead) * &
    &                                                 v_beads_nm_half(i_coord,i_bead,i_atom)/6.4d+1
                 v_beads_nm_last(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * lamda_nm(i_bead) * &
    &                                                 v_beads_nm_last(i_coord,i_bead,i_atom)/6.4d+1
               end if
             else if(use_rpmd) then
               v_beads_nm_half(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * &
    &                                                 v_beads_nm_half(i_coord,i_bead,i_atom)
               v_beads_nm_last(i_coord,i_bead,i_atom) = species_m(species(i_atom)) * &
    &                                                 v_beads_nm_last(i_coord,i_bead,i_atom)
             else
               write(stderr,*) 'Must use PIMD.or.CMD.or.RPMD'
               stop
             end if
           end do
         end do
       end do     
    else
       write(stderr,*) 'only Nose-Hoover supplied by PIMD at the moment'
       stop
    end if

    PIMD_Epot_last = Epot
!   PIMD_schedule_step = 1

  !------------------------------------------------------------------------------------------------------
  ! initialize PIMD counters
    tsystem_beads      = 0d0
    PIMD_maxsteps = int(PIMD_time/PIMD_tstep)
    use_PIMD_max_steps = .true.
  
  !-------------------------------------------------------------------------------------------------

    total_forces_beads_nm = 0.0d0
    spring_forces_beads_nm = 0.0d0

  ! first integration step:
    do i_atom = 1, n_atoms, 1
      do i_bead = 1, n_beads, 1
        do j_bead = 1, n_beads, 1
          do i_coord = 1, 3, 1
            total_forces_beads_nm(i_coord,i_bead,i_atom) = total_forces_beads_nm(i_coord,i_bead,i_atom) + &
   &                                 1/sqrt(1.0d0*n_beads) * &
   &                                 transform_u(i_bead,j_bead) * total_forces_beads(i_coord,j_bead,i_atom)
          end do
        end do
        do i_coord = 1, 3, 1
          spring_forces_beads_nm(i_coord,i_bead,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                              mass_atomic_unit * species_m(species(i_atom)) * lamda_nm(i_bead) *  &
   &                                              coords_beads_nm(i_coord,i_bead,i_atom)
        end do
      end do
    end do

    do i_bead = 1, n_beads, 1
      if (PIMD_ensemble.eq.'NVT_nose-hoover') then
        BeadMD_forces(:,:) = 0.0d0
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            BeadMD_forces(i_coord,i_atom) = (total_forces_beads_nm(i_coord,i_bead,i_atom) +  &
   &                                        spring_forces_beads_nm(i_coord,i_bead,i_atom))/MD_KE_factor
            rb_last(i_coord,i_atom) = r_beads_nm_last(i_coord,i_bead,i_atom)
            vb_last(i_coord,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)
            vb_half(i_coord,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)
          end do
        end do
        call initial_step_NH96_GLA_pi(PIMD_tstep,PIMD_Q_NP,rb_last,vb_half,vb_last, &
   &                 s_NP_beads_half,s_dot_NP_beads_half,shat,BeadMD_forces,PIMD_temperature)
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            r_beads_nm_last(i_coord,i_bead,i_atom) = rb_last(i_coord,i_atom)
            coords_beads_nm(i_coord,i_bead,i_atom) = rb_last(i_coord,i_atom)
            v_beads_nm_last(i_coord,i_bead,i_atom) = vb_last(i_coord,i_atom)
            v_beads_nm_half(i_coord,i_bead,i_atom) = vb_half(i_coord,i_atom)
          end do
        end do
      else 
        write(stderr,*) 'only Nose-Hoover is supplied in PIMD now'
        stop
      end if
    end do

    v_beads_nm_half_orig(:,:,:) = 0.0d0
    v_beads_nm_last_orig(:,:,:) = 0.0d0

    do i_bead = 1, n_beads, 1
      do i_atom = 1, n_atoms, 1
        do i_coord = 1, 3, 1
          if(use_pimd) then
            if(i_bead.eq.1) then
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
            else
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                    (species_m(species(i_atom)) * lamda_nm(i_bead))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                    (species_m(species(i_atom)) * lamda_nm(i_bead))
            end if
          else if(use_cmd) then
            if(i_bead.eq.1) then
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
            else
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)*6.4d+1/  &
    &                                                  (species_m(species(i_atom)) * lamda_nm(i_bead))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)*6.4d+1/  &
    &                                                  (species_m(species(i_atom)) * lamda_nm(i_bead))
            end if
          else if(use_rpmd) then 
            v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                species_m(species(i_atom))
            v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                species_m(species(i_atom))
          else
            write(stderr,*) 'Must use PIMD.or.CMD.or.RPMD'
            stop
          end if
        end do
      end do
    end do

    coords_beads = 0.0d0
    r_beads_last = 0.0d0
    v_beads_half = 0.0d0
    v_beads_last = 0.0d0
    do j_bead = 1, n_beads, 1
      do i_atom = 1, n_atoms, 1
        do i_bead = 1, n_beads, 1
          do i_coord = 1, 3, 1
            coords_beads(i_coord,j_bead,i_atom) = coords_beads(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * coords_beads_nm(i_coord,i_bead,i_atom)
            v_beads_half(i_coord,j_bead,i_atom) = v_beads_half(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * v_beads_nm_half_orig(i_coord,i_bead,i_atom)
            v_beads_last(i_coord,j_bead,i_atom) = v_beads_last(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * v_beads_nm_last_orig(i_coord,i_bead,i_atom)
            r_beads_last(i_coord,j_bead,i_atom) = r_beads_last(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * r_beads_nm_last(i_coord,i_bead,i_atom)
          end do
        end do
      end do
    end do

  else  !  use nm or not

    if(.not.PIMD_initialconf_from_restart) then
      !---------------------------------------------------------------------------------------------------
      ! initialize velocities from MB distribution if requested
      MB_velocity_initialization = .true.
      if (MB_velocity_initialization) then
        do i_atom = 1, n_atoms, 1
          do j_bead = 1, n_beads, 1
            do i_coord = 1, 3, 1
              call MB_sample_box_mueller_pi(PIMD_init_temperature,species_m(species(i_atom))/(1.0d0*n_beads), &
   &                                        v_beads_half(i_coord,j_bead,i_atom))
            end do
          end do
        end do
      end if
      v_beads_last(:,:,:) = v_beads_half(:,:,:)
      r_beads_last(:,:,:) = coords_beads(:,:,:)
    end if
    ! Below are the lines for propogation in the Cartesian coordinates 
    if (PIMD_ensemble.eq.'NVT_nose-hoover') then
      s_NP_beads_last = 0d0
      s_NP_beads_half = 0d0
      s_NP_beads      = 0d0
      do i_atom = 1, n_atoms, 1
        do j_bead = 1, n_beads, 1
          do i_coord = 1, 3, 1
            v_beads_half(i_coord,j_bead,i_atom) = (species_m(species(i_atom))/(1.0d0*n_beads)) *   &
    &                                              v_beads_half(i_coord,j_bead,i_atom)
            v_beads_last(i_coord,j_bead,i_atom) = (species_m(species(i_atom))/(1.0d0*n_beads)) *   &
    &                                              v_beads_last(i_coord,j_bead,i_atom)
          end do
        end do
      end do
    else
      write(stderr,*) 'only Nose-Hoover supplied by PIMD at the moment'
      stop
    end if

    spring_forces_beads(:,:,:) = 0.0d0
    do i_atom = 1, n_atoms, 1
      do j_bead = 2, n_beads - 1, 1
        do i_coord = 1, 3, 1
          spring_forces_beads(i_coord,j_bead,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                           mass_atomic_unit * species_m(species(i_atom)) *  &
   &                                           (coords_beads(i_coord,j_bead+1,i_atom) + coords_beads(i_coord,j_bead-1,i_atom) -  &
   &                                           2.0d0 * coords_beads(i_coord,j_bead,i_atom))
        end do
      end do
      do i_coord = 1, 3, 1
        spring_forces_beads(i_coord,1,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                           mass_atomic_unit * species_m(species(i_atom)) *  &
   &                                           (coords_beads(i_coord,n_beads,i_atom) + coords_beads(i_coord,2,i_atom) -  &
   &                                           2.0d0 * coords_beads(i_coord,1,i_atom))  
        spring_forces_beads(i_coord,n_beads,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                           mass_atomic_unit * species_m(species(i_atom)) *  &
   &                                           (coords_beads(i_coord,n_beads-1,i_atom) + coords_beads(i_coord,1,i_atom) -  &
   &                                           2.0d0 * coords_beads(i_coord,n_beads,i_atom))
      end do
    end do

    do j_bead = 1, n_beads, 1
      if (PIMD_ensemble.eq.'NVT_nose-hoover') then
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            BeadMD_forces(i_coord,i_atom) = (total_forces_beads(i_coord,j_bead,i_atom)+  & 
   &                                         spring_forces_beads(i_coord,j_bead,i_atom))/MD_KE_factor
            rb_last(i_coord,i_atom) = r_beads_last(i_coord,j_bead,i_atom)
            vb_last(i_coord,i_atom) = v_beads_last(i_coord,j_bead,i_atom)
            vb_half(i_coord,i_atom) = v_beads_half(i_coord,j_bead,i_atom)
          end do
        end do

        call initial_step_NH96_GLA_pi(PIMD_tstep,PIMD_Q_NP,rb_last,vb_half,vb_last, &
  &                  s_NP_beads_half,s_dot_NP_beads_half,shat,BeadMD_forces,PIMD_temperature)
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            r_beads_last(i_coord,j_bead,i_atom) = rb_last(i_coord,i_atom)
            coords_beads(i_coord,j_bead,i_atom) = rb_last(i_coord,i_atom)
            v_beads_last(i_coord,j_bead,i_atom) = vb_last(i_coord,i_atom)
            v_beads_half(i_coord,j_bead,i_atom) = vb_half(i_coord,i_atom)
          end do
        end do
      else
        write(stderr,*) 'only Nose-Hoover is supplied in PIMD now'
        stop
      end if
    end do

  end if
 
  tsystem_beads_last = tsystem_beads
  tsystem_beads = tsystem_beads + PIMD_tstep
  PIMD_stepcount = PIMD_stepcount + 1
  PIMD_force_evaluations = PIMD_force_evaluations +  1
    !stop
end subroutine initialize_PIMD
!******
  
!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/PIMD_step
!  NAME
!    PIMD_step
!  SYNOPSIS
subroutine PIMD_step(Epot)
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
  integer :: i_coord, i_atom, i_bead, p, j_bead, nm
  real*8, dimension(3,n_atoms) :: BeadMD_forces
  character*100 :: info_str

  PIMD_Epot_last = Epot

 
  if (use_nm) then  
    ! p = (n_beads - 2)/2

    total_forces_beads_nm(:,:,:) = 0.0d0

    ! first integration step:
    do i_atom = 1, n_atoms, 1
      do i_bead = 1, n_beads, 1
        do j_bead = 1, n_beads, 1
          do i_coord = 1, 3, 1
            total_forces_beads_nm(i_coord,i_bead,i_atom) = total_forces_beads_nm(i_coord,i_bead,i_atom) + &
   &                                 1/sqrt(1.0d0*n_beads) * &
   &                                 transform_u(i_bead,j_bead) * total_forces_beads(i_coord,j_bead,i_atom)
          end do
        end do
        do i_coord = 1, 3, 1
          spring_forces_beads_nm(i_coord,i_bead,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                              mass_atomic_unit * species_m(species(i_atom)) * lamda_nm(i_bead) *    &
   &                                              coords_beads_nm(i_coord,i_bead,i_atom)
        end do
      end do
    end do
  
    do i_bead = 1, n_beads, 1
      if (PIMD_ensemble.eq.'NVT_nose-hoover') then
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            BeadMD_forces(i_coord,i_atom) = (total_forces_beads_nm(i_coord,i_bead,i_atom) +  &
   &                                spring_forces_beads_nm(i_coord,i_bead,i_atom))/MD_KE_factor
            rb_last(i_coord,i_atom) = r_beads_nm_last(i_coord,i_bead,i_atom)
            vb_last(i_coord,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)
            vb_half(i_coord,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)
          end do
        end do
        call NH96_generalized_leap_frog_pi(PIMD_tstep,PIMD_Q_NP,rb_last,vb_half,vb_last, &
   &             s_dot_NP_beads_half,s_NP_beads_half,BeadMD_forces,PIMD_temperature,shat1,shat2,     &
                 s_NP_beads_last,s_dot_NP_beads_last)
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            r_beads_nm_last(i_coord,i_bead,i_atom) = rb_last(i_coord,i_atom)
            coords_beads_nm(i_coord,i_bead,i_atom) = rb_last(i_coord,i_atom)
            v_beads_nm_last(i_coord,i_bead,i_atom) = vb_last(i_coord,i_atom)
            v_beads_nm_half(i_coord,i_bead,i_atom) = vb_half(i_coord,i_atom)
          end do
        end do
      else
        write(stderr,*) 'only Nose-Hoover is supplied in PIMD now'
        stop
      end if
    end do

    v_beads_nm_half_orig(:,:,:) = 0.0d0
    v_beads_nm_last_orig(:,:,:) = 0.0d0

    do i_bead = 1, n_beads, 1
      do i_atom = 1, n_atoms, 1
        do i_coord = 1, 3, 1
          if(use_pimd) then
            if(i_bead.eq.1) then
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
            else
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                    (species_m(species(i_atom)) * lamda_nm(i_bead))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                    (species_m(species(i_atom)) * lamda_nm(i_bead))
            end if
          else if(use_cmd) then
            if(i_bead.eq.1) then
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                    species_m(species(i_atom))
            else
              v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)*6.4d+1/  &
    &                                                  (species_m(species(i_atom)) * lamda_nm(i_bead))
              v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)*6.4d+1/  &
    &                                                  (species_m(species(i_atom)) * lamda_nm(i_bead))
            end if
          else if(use_rpmd) then
            v_beads_nm_half_orig(i_coord,i_bead,i_atom) = v_beads_nm_half(i_coord,i_bead,i_atom)/  &
    &                                                species_m(species(i_atom))
            v_beads_nm_last_orig(i_coord,i_bead,i_atom) = v_beads_nm_last(i_coord,i_bead,i_atom)/  &
    &                                                species_m(species(i_atom))
          else
            write(stderr,*) 'Must use PIMD.or.CMD.or.RPMD'
            stop
          end if
        end do
      end do
    end do

    coords_beads = 0.0d0
    r_beads_last = 0.0d0
    v_beads_half = 0.0d0
    v_beads_last = 0.0d0
    do j_bead = 1, n_beads, 1
      do i_atom = 1, n_atoms, 1
        do i_bead = 1, n_beads, 1
          do i_coord = 1, 3, 1
            coords_beads(i_coord,j_bead,i_atom) = coords_beads(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * coords_beads_nm(i_coord,i_bead,i_atom)
            v_beads_half(i_coord,j_bead,i_atom) = v_beads_half(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * v_beads_nm_half_orig(i_coord,i_bead,i_atom)
            v_beads_last(i_coord,j_bead,i_atom) = v_beads_last(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * v_beads_nm_last_orig(i_coord,i_bead,i_atom)
            r_beads_last(i_coord,j_bead,i_atom) = r_beads_last(i_coord,j_bead,i_atom) + sqrt(1.0d0*n_beads) * &
 &                            transform_u(i_bead,j_bead) * r_beads_nm_last(i_coord,i_bead,i_atom)
          end do
        end do
      end do
    end do

  else  ! use nm or not

    !   in below we propagate in the Cardesian coordiantes
    spring_forces_beads(:,:,:) = 0.0d0
    do i_atom = 1, n_atoms, 1
      do j_bead = 2, n_beads - 1, 1
        do i_coord = 1, 3, 1
          spring_forces_beads(i_coord,j_bead,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                             mass_atomic_unit * species_m(species(i_atom)) *  & 
   &                                             (coords_beads(i_coord,j_bead+1,i_atom) + coords_beads(i_coord,j_bead-1,i_atom) -  &
   &                                             2.0d0 * coords_beads(i_coord,j_bead,i_atom))
        end do
      end do
      do i_coord = 1, 3, 1
        spring_forces_beads(i_coord,1,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                             mass_atomic_unit * species_m(species(i_atom)) *  & 
   &                                             (coords_beads(i_coord,n_beads,i_atom) + coords_beads(i_coord,2,i_atom) -  &
   &                                             2.0d0 * coords_beads(i_coord,1,i_atom))  
        spring_forces_beads(i_coord,n_beads,i_atom) = - n_beads * (boltzmann_kB * PIMD_temperature )**2 * &
   &                                             mass_atomic_unit * species_m(species(i_atom)) *  &
   &                                             (coords_beads(i_coord,n_beads-1,i_atom) + coords_beads(i_coord,1,i_atom) -  &
   &                                             2.0d0 * coords_beads(i_coord,n_beads,i_atom))
      end do
    end do

    do j_bead = 1, n_beads, 1
      if (PIMD_ensemble.eq.'NVT_nose-hoover') then
        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            BeadMD_forces(i_coord,i_atom) = (total_forces_beads(i_coord,j_bead,i_atom)+  &
   &                                        spring_forces_beads(i_coord,j_bead,i_atom))/MD_KE_factor
            rb_last(i_coord,i_atom) = r_beads_last(i_coord,j_bead,i_atom)
            vb_last(i_coord,i_atom) = v_beads_last(i_coord,j_bead,i_atom)
            vb_half(i_coord,i_atom) = v_beads_half(i_coord,j_bead,i_atom)
          end do
        end do

        call NH96_generalized_leap_frog_pi(PIMD_tstep,PIMD_Q_NP,rb_last,vb_half,vb_last, &
   &           s_dot_NP_beads_half,s_NP_beads_half,BeadMD_forces,PIMD_temperature,shat1,shat2,     &
   &           s_NP_beads_last,s_dot_NP_beads_last)

        do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
            r_beads_last(i_coord,j_bead,i_atom) = rb_last(i_coord,i_atom)
            coords_beads(i_coord,j_bead,i_atom) = rb_last(i_coord,i_atom)
            v_beads_last(i_coord,j_bead,i_atom) = vb_last(i_coord,i_atom)
            v_beads_half(i_coord,j_bead,i_atom) = vb_half(i_coord,i_atom)
          end do
        end do
      else
        write(stderr,*) 'only Nose-Hoover is supplied in PIMD now'
        stop
      end if
    end do

  end if

  PIMD_stepcount = PIMD_stepcount + 1

  tsystem_beads_last   = tsystem_beads

  tsystem_beads = tsystem_beads + PIMD_tstep

end subroutine PIMD_step
!******

!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/NH96_generalized_leap_frog_pi
!  NAME
!    NH96_generalized_leap_frog
!  SYNOPSIS
!     call NH96_generalized_leap_frog(MD_tstep,MD_Q_NP,coords,v_half,v_last, &
!          s_dot_NP_half,s_NP_half,MD_forces,MD_temperature,shat1,shat2,     &
!          s_NP_last,s_dot_NP_last)  
subroutine NH96_generalized_leap_frog_pi(deltat,Q,r,p_half,p_last,pi_half,&
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
  gkT          = dble(PIMD_g_DOF)*boltzmann_kB*T
  pstar(:,:)   = p_half(:,:) + deltat*f(:,:)/2d0
  call calculate_kinetic_energy_nose_pi(pstar,1d0,Ekin)
  pidoublestar = pi_half + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  shat1        = exp(-pidoublestar*deltat/(2d0*Q))
  eta_last     = eta_half + pidoublestar*deltat/(2d0*Q)
  p_last(:,:)  = pstar(:,:)*shat1
  call calculate_kinetic_energy_nose_pi(p_last,1d0,Ekin)
  pi_last      = pidoublestar + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  call calculate_kinetic_energy_nose_pi(p_last,1d0,Ekin)
  pistar       = pi_last + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  shat2        = exp(-pistar*deltat/(2d0*Q))
  pkin(:,:)    = p_last(:,:) * shat2
  call calculate_kinetic_energy_nose_pi(pkin,1d0,Ekin)
  pi_half      = pistar + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  eta_half     = eta_last + pistar*deltat/(2d0*Q)
  p_half(:,:)  = p_last(:,:)*shat2+f(:,:)*deltat/2d0

!  MB_clean_rotations = .false.
  if (MB_clean_rotations) call clean_velocities_NH_pi(p_half)   ! remove rotations and translations from velocities, hopefully	

  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*p_half(:,i_atom)/species_m(species(i_atom))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms 
       call map_to_center_cell(r(:,i_atom))
     end do
  end if
end subroutine NH96_generalized_leap_frog_pi
!******
!------------------------------------------------------------------------------
!****s* molecular_dynamics/initial_step_NH96_GLA_pi
!  NAME
!    initial_step_NH96_GLA - GLA first step for Nosé-Hoover MD
!  SYNOPSIS
subroutine initial_step_NH96_GLA_pi(deltat,Q,r,p_half,p_last,eta_half,pi_half,&
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
  gkT         = dble(PIMD_g_DOF)*boltzmann_kB*T
  call calculate_kinetic_energy_nose_pi(p_half,1d0,Ekin)
  pistar      = deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  eta_half    = pistar*deltat/(2d0*Q)
  shat        = exp(-eta_half)
  pkin(:,:)   = p_half(:,:)*shat
  call calculate_kinetic_energy_nose_pi(pkin,1d0,Ekin)
  pi_half     = pistar + deltat*(2d0*Ekin-gkT)/(4d0*MD_KE_factor)
  p_half(:,:) = p_half(:,:)*shat + f(:,:)*deltat/2d0
  do i_atom = 1, n_atoms
     r(:,i_atom) = r(:,i_atom) + deltat*p_half(:,i_atom)/species_m(species(i_atom))
  end do
  if (n_periodic.gt.0) then
     do i_atom = 1, n_atoms
       call map_to_center_cell(r(:,i_atom))
     end do
  end if
end subroutine initial_step_NH96_GLA_pi

!******
!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/calculate_kinetic_energy_nose_pi
!  NAME
!    calculate_kinetic_energy_nose_pi 
!  SYNOPSIS
subroutine calculate_kinetic_energy_nose_pi(p,s,Ekin)
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
end subroutine calculate_kinetic_energy_nose_pi

!******

!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/calculate_kinetic_energy_pi
!  NAME
!    calculate_kinetic_energy_pi
!  SYNOPSIS
subroutine calculate_kinetic_energy_pi(v,Ekin,start_i,end_i)
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
end subroutine calculate_kinetic_energy_pi

!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/clean_velocities_pi
!  NAME
!    clean_velocities
!  SYNOPSIS
subroutine clean_velocities_pi(v, clean_despite_constraints)
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
     if(output_level .eq. 'PIMD_light') output_priority = 2
     ! get initial kinetic energy
     call calculate_kinetic_energy_pi(v,E_kin_in)
     if (.not.MB_clean_rotations) then
        write(info_str,'(2X,A)') 'Cleaning translations from velocities'
        call localorb_info(info_str, use_unit,'(A)',OL_norm)
        ! simple hack to remove all translations from the momenta
        p_total = 0d0
        do i_atom = 1, n_atoms
           p_total(:) = p_total(:) + species_m(species(i_atom))*v(:,i_atom)
        end do
        p_total(:) = p_total(:)/dble(n_atoms_PIMD)
        do i_atom = 1, n_atoms_PIMD
!           if (.not.constrain_PIMD(i_atom)) then  ! prevent constrained atoms from acquiring velocities ... 
            v(:,i_atom) = v(:,i_atom) - p_total(:)/species_m(species(i_atom))
!           end if
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
        do i_atom = 1, n_atoms, 1
           do i_coord = 1, 3, 1
              translation_vectors(i_coord, i_atom, i_coord) = 1.d0
           end do
        end do
        do i_coord =1, 3, 1
           inv_norm = 1.d0 / dnrm2(3*n_atoms, translation_vectors(:, :, i_coord), 1)
           translation_vectors(:, :, i_coord) = translation_vectors(:, :, i_coord) * inv_norm
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
     call calculate_kinetic_energy_pi(v,E_kin_out)
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
     if(output_level .eq. 'PIMD_light') output_priority = 1
  end if
end subroutine clean_velocities_pi

!******


!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/clean_velocities_pi
!  NAME
!    clean_velocities_pi
!  SYNOPSIS
subroutine clean_velocities_NH_pi(v)
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
  call calculate_kinetic_energy_nose_pi(v,1d0,E_kin_in)
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
     do i_atom = 1, n_atoms, 1
        do i_coord = 1, 3, 1
           translation_vectors(i_coord, i_atom, i_coord) = 1.d0
        end do
     end do
     do i_coord =1, 3, 1
        inv_norm = 1.d0 / dnrm2(3*n_atoms, translation_vectors(:, :, i_coord), 1)
        translation_vectors(:, :, i_coord) = translation_vectors(:, :, i_coord) * inv_norm
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
  call calculate_kinetic_energy_nose_pi(v,1d0,E_kin_out)
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

end subroutine clean_velocities_NH_pi
!******


!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/clean_momenta_pi
!  NAME
!    clean_momenta_pi
!  SYNOPSIS
subroutine clean_momenta_pi(n_atoms, p)
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
end subroutine clean_momenta_pi
!!******

!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/MB_sample_box_mueller_pi
!  NAME
!    MB_sample_box_mueller_pi
!  SYNOPSIS
subroutine MB_sample_box_mueller_pi(T,m,v)
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
  real*8 :: T,m,v,vmax,rand1,rand2
  !integer :: seed !Why declaring this seed here?
  logical, save :: saved = .false.
  real*8, save :: r1, r2
  if (.not.saved) then
     saved = .true.
     call random_number(r1) 
     call random_number(r2)
     v     = sqrt(boltzmann_kB*T/(MD_KE_factor*m))*sqrt(-2d0*dlog(r1))*cos(2d0*pi*r2)
  else
     saved = .false.
     v     = sqrt(boltzmann_kB*T/(MD_KE_factor*m))*sqrt(-2d0*dlog(r1))*sin(2d0*pi*r2)
  end if
end subroutine MB_sample_box_mueller_pi
!******

!****s* pi_molecular_dynamics/random_own_pi
!  NAME
!    random_own_path_in
!  SYNOPSIS
real*8  function random_own_pi()
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
  integer :: RNGB_num
  real*8  :: rnum
  logical :: user_seeded
  integer*8,save :: Y = 0
  character*8  :: cdate
  character*10 :: ctime
  ! hardwired warning in case someone ever tries again to call the
  ! rng without synchronizing between different tasks
  if (myid.ne.0) then
    write(stderr,'(1X,A,I8,A)') "* Warning: Random number generator invoked from task myid = ", myid, "."
    write(stderr,'(1X,A)')      "* To guarantee that all results are properly synchronized between different threads, "
    write(stderr,'(1X,A)')      "* random_own() must never be called from any other task than myid = 0 . "
    write(stderr,'(1X,A)')      "* This appears to be a programming error - please correct."
    stop
  end if
  if (PIMD_RNG_firstcall) then
     ! seed was already set when rng was initialized
     Y = seed
     PIMD_RNG_firstcall = .false.
  end if
  if (.not.PIMD_RNG_seeded) then
     ! if a seed was not specified, we randomize the rng variables used also
     call date_and_time(cdate, ctime)
     read(ctime,'(F10.3)') rnum
     RNGB_num = mod(int(rnum),7)+1
  else
     ! if a seed was specified, we keep the rng variables fixed also
     ! this ensures always the same random number sequence for the same seed
     RNGB_num = 7
  end if
  Y       = mod(RNGB_a(RNGB_num)*Y+RNGB_c(RNGB_num),RNGB_m(RNGB_num))
  random_own_pi = dble(Y)/dble(RNGB_m(RNGB_num))
end function random_own_pi
!******


!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/initialize_RNG_pi
!  NAME
!    initialize_RNG_pi
!  SYNOPSIS
subroutine initialize_RNG_pi
!  PURPOSE
!    initialization of random number generator
!  USES
  use runtime_choices
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
  ! some compilers handle 8 byte integers but can not interpret them
  ! as constants when you type out the digits
  ! Fix by calculating the RNG initialization from "scratch"
  temp = 1
  do i = 1, 31
     temp = temp * 2
  end do
  RNGB_m(:) = temp*2 ! 4294967296
  RNGB_m(7) = temp-1 ! 2147483647
  RNGB_a(1) = 1664525
  RNGB_a(2) = 22695477
  RNGB_a(3) = 69069
  RNGB_a(4) = 1103515245
  RNGB_a(5) = 134775813
  RNGB_a(6) = 214013
  RNGB_a(7) = 16807
  RNGB_c(1) = 1013904223
  RNGB_c(2) = 1
  RNGB_c(3) = 5
  RNGB_c(4) = 12345
  RNGB_c(5) = 1
  RNGB_c(6) = 2531011
  RNGB_c(7) = 0
  ! seed RNG from system time, this might be required for the momentum-initialization 
  if (.not.PIMD_RNG_seeded) then
     ! initialize random number generator
     call date_and_time(cdate, ctime)
     read(ctime,*) rng_seed
     seed = int(rng_seed)
  end if
  if (seed.lt.0) seed = - seed
  if (seed.eq.0) seed = seed + 1  
end subroutine initialize_RNG_pi
!******

!****s* molecular_dynamics/write_PIMD_restart
!  NAME
!    write_PIMD_restart
!  SYNOPSIS
subroutine write_PIMD_restart
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
  real*8 :: s_beads_write, s_beads_last_write
  real*8, dimension(3,n_beads,n_atoms) :: v_beads_half_write, v_beads_last_write
  integer i_atom, i_bead

   if(use_nm) then
     do i_atom = 1, n_atoms
       do i_bead = 1, n_beads
         v_beads_half_write(:,i_bead,i_atom) = v_beads_half(:,i_bead,i_atom)
         v_beads_last_write(:,i_bead,i_atom) = v_beads_last(:,i_bead,i_atom)
       end do
     end do
   else
     do i_atom = 1, n_atoms
       do i_bead = 1, n_beads
         v_beads_half_write(:,i_bead,i_atom) = n_beads * v_beads_half(:,i_bead,i_atom)/species_m(species(i_atom))
         v_beads_last_write(:,i_bead,i_atom) = n_beads * v_beads_last(:,i_bead,i_atom)/species_m(species(i_atom))
       end do
     end do
   end if

  ! only write from the zero task ... 

  if (myid.eq.0) then
     s_beads_write      = s_NP_beads
     s_beads_last_write = s_NP_beads_last
     if (PIMD_ensemble.eq.'NVT_nose-hoover') then
        s_beads_write      = dexp(s_beads_write)
        s_beads_last_write = dexp(s_beads_last_write)
     end if

     open(file = PIMD_restart_file, unit = 89, status = 'replace', action='write')
     write(89,'(A)') '############### aims_PIMD_restart.dat #######################'
     write(89,'(A)') '# This file contains all information that is necessary to #'
     write(89,'(A)') '# continue a PIMD run. The data in this file should not   #'
     write(89,'(A)') '# be altered by the user, unless you REALLY know what you #'
     write(89,'(A)') '# are doing.                                              #'
     write(89,'(A)') '###########################################################'
     write(89,'(A,I6)')        'n_atoms              ',n_atoms
     write(89,'(A,I6)')        'n_beads              ',n_beads
     do i_atom = 1, n_atoms
       write(89,'(A,I6)')      'atom number          ',i_atom
       do i_bead=1, n_beads
         write(89,'(A,3E30.20)') 'coords_beads       ', coords_beads(:,i_bead,i_atom)
       end do
       do i_bead=1, n_beads
         write(89,'(A,3E30.20)') 'v_beads_half       ',v_beads_half_write(:,i_bead,i_atom)
       end do
       do i_bead=1, n_beads
         write(89,'(A,3E30.20)') 'v_beads_last       ',v_beads_last_write(:,i_bead,i_atom)
       end do
       do i_bead=1, n_beads
         write(89,'(A,3E30.20)') 'r_beads_last       ',r_beads_last(:,i_bead,i_atom)
       end do
     enddo
     write(89,'(A,E30.20)')    's_NP_beads            ', s_beads_write
     write(89,'(A,E30.20)')    's_NP_beads_half       ', s_NP_beads_half
     write(89,'(A,E30.20)')    's_NP_beads_last       ', s_beads_last_write
     write(89,'(A,E30.20)')    's_dot_NP_beads_half   ', s_dot_NP_beads_half
     write(89,'(A,E30.20)')    's_dot_NP_beads_last   ', s_dot_NP_beads_last
     write(89,'(A,E30.20)')    'tsystem_beads         ', tsystem_beads
     write(89,'(A,E30.20)')    'tsystem_beads_last    ', tsystem_beads_last
     write(89,'(A,E30.20)')    'PIMD_H0               ', PIMD_H0
     write(89,'(A,E30.20)')    'PIMD_Epot_last        ', PIMD_Epot_last
     write(89,'(A,I6)')        'PIMD_stepcount        ', PIMD_stepcount
     write(89,'(A,I6)')        'PIMD_force_evaluations', PIMD_force_evaluations
     close(unit=89)
  end if
end subroutine write_PIMD_restart

!------------------------------------------------------------------------------
!****s* pi_molecular_dynamics/initialize_nm
!  NAME
!    initialize_nm
!  SYNOPSIS
subroutine initialize_nm
!  PURPOSE
!    initialize the normal mode matrices and arrays
!  USES
  use dimensions
  use constants
!  AUTHOR
!    XZL
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
  integer*8    :: i_bead, j_bead, p, k_bead, lwork, info1
  real*8, dimension(:), allocatable :: work, w
  real*8, dimension(:,:), allocatable :: matrix_1, matrix_2

  allocate(matrix_1(1:n_beads,1:n_beads))
  allocate(matrix_2(1:n_beads,1:n_beads))

  matrix_1 = 0.0d0
  matrix_1(1,1) = 2.0d0 * n_beads
  matrix_1(1,2) = - n_beads
  matrix_1(1,n_beads) = - n_beads
  matrix_1(n_beads,n_beads) = 2.0d0 * n_beads
  matrix_1(n_beads,n_beads-1) = - n_beads
  matrix_1(n_beads,1) = - n_beads
  do i_bead = 2, n_beads-1
    matrix_1(i_bead, i_bead-1) = - n_beads
    matrix_1(i_bead, i_bead) = 2.0d0*n_beads
    matrix_1(i_bead, i_bead+1) = - n_beads
  end do

  matrix_2 = matrix_1

  lwork = 3*n_beads + 1
  allocate(work(1:lwork))
  allocate(w(1:n_beads))

  call dsyev('V', 'U', n_beads, matrix_2, n_beads, w, work, lwork, info1)   

  do i_bead = 1, n_beads
    do j_bead = 1, n_beads
      transform_u(i_bead,j_bead) = matrix_2(j_bead,i_bead)
    end do
  enddo
  
  lamda_nm(1:n_beads) = w(1:n_beads)

  do j_bead = 2, n_beads-1
    matrix_a(j_bead,j_bead-1) = -1.0d0
    matrix_a(j_bead,j_bead) = 2.0d0
    matrix_a(j_bead,j_bead+1) = -1.0d0
  end do   
  matrix_a(1,1) = 2.0d0
  matrix_a(1,2) = -1.0d0
  matrix_a(1,n_beads) = -1.0d0

  matrix_a(n_beads,1) = -1.0d0
  matrix_a(n_beads,n_beads-1) = -1.0d0
  matrix_a(n_beads,n_beads) = 2.0d0

  do j_bead = 1, n_beads
    do i_bead = 1, n_beads
      do k_bead = 1, n_beads
        matrix_au(j_bead,i_bead) = matrix_au(j_bead,i_bead) + matrix_a(j_bead, k_bead) * &
     &                             transform_u(i_bead,k_bead)
      end do 
    end do
  end do

  do k_bead = 1, n_beads
    do i_bead = 1, n_beads
      do j_bead = 1, n_beads
        matrix_lamda(k_bead,i_bead) = matrix_lamda(k_bead,i_bead) + n_beads * transform_u(k_bead,j_bead) * &
     &                                matrix_au(j_bead,i_bead)
      end do
    end do
  end do

end subroutine initialize_nm

end module pi_molecular_dynamics
