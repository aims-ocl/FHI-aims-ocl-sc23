!****h* FHI-aims/contour_def_gw_types
!  NAME
!    types for contour deformation in single-shot GW
!  SYNOPSIS

module contour_def_gw_types

   implicit none

   private

   public :: cd_environment_type, spectral_environment_type, qp_adjust_type
   public :: deallocate_cd_environment
 
! **************************************************************************************************
  type complex_type
     real(kind=8), dimension(:,:), allocatable           :: re
     real(kind=8), dimension(:,:), allocatable           :: im
     complex(kind=8), dimension(:,:), allocatable        :: complx
  end type complex_type

! **************************************************************************************************
  type spectral_environment_type
     integer                                             :: npoints_global
     integer                                             :: npoints_local
     integer, dimension(:), allocatable                  :: freq_index
     real(kind=8), dimension(:), allocatable             :: freqs
     real(kind=8), dimension(:), allocatable             :: spectrum                 
  end type spectral_environment_type

! **************************************************************************************************
  type sc_environment_type
     integer                                             :: n_occ          ! # explicit occ
     integer                                             :: n_virt         ! # explicit virt
     logical                                             :: reiterate
     logical                                             :: last_step
     logical, dimension(:), allocatable                  :: qp_converged        ! sc loop convergence
     logical, dimension(:), allocatable                  :: converged_last_step ! sc loop convergence
  end type sc_environment_type

! **************************************************************************************************
   type qp_adjust_type
      integer                                            :: n_qps
      integer, dimension(:), allocatable                 :: state
      real(kind=8), dimension(:), allocatable            :: values
   end type qp_adjust_type

! **************************************************************************************************

   type cd_environment_type
     character(len=12)                                   :: sctype
     character(len=40)                                   :: restart
     logical                                             :: iterative
     logical                                             :: full_cmplx_sigma
     logical                                             :: self_consistent
     logical                                             :: try_zshot
     logical                                             :: calc_single_state_spec
     real(kind=8)                                        :: contour_def_offset
     real(kind=8)                                        :: eta
     real(kind=8)                                        :: max_freq_spec
     real(kind=8)                                        :: min_freq_spec
     real(kind=8)                                        :: spec_resolution
     integer                                             :: state_spec
     integer                                             :: n_iter_sc   ! # of iterations sc loop
     integer                                             :: num_levels_tot
     integer                                             :: state_to_shift ! state used for deltaE
     integer                                             :: num_residues
     integer, dimension(:), allocatable                  :: spin_channel
     integer, dimension(:), allocatable                  :: contour_def_start
     integer, dimension(:), allocatable                  :: contour_def_end
     integer, dimension(:,:), allocatable                :: corrected_levels
     integer, dimension(:,:), allocatable                :: index_cd_levels
     integer, dimension(:), allocatable                  :: num_levels
     type(qp_adjust_type)                                :: qp_adjust                               
     type(complex_type), allocatable                     :: self_energy
     real(kind=8), dimension(:), allocatable             :: real_freq     
     real(kind=8), dimension(:,:), allocatable           :: shift  
     integer, dimension(:), allocatable                  :: residue_from_freq   
     type(spectral_environment_type), dimension(:),&
      allocatable                                        :: spectral_env
     type(sc_environment_type), dimension(:),&
      allocatable                                        :: sc_env
   end type cd_environment_type


! **************************************************************************************************

contains

! **************************************************************************************************
!> brief deallocates contour deformation environment
!  o gw_cd --contour deformation environment
! **************************************************************************************************

  
   subroutine deallocate_cd_environment(gw_cd)

      type(cd_environment_type)                          :: gw_cd

      integer                                            :: i

      if(allocated(gw_cd%self_energy)) then
        if(allocated(gw_cd%self_energy%re)) then
           deallocate(gw_cd%self_energy%re)
        endif
        if(allocated(gw_cd%self_energy%im)) then
           deallocate(gw_cd%self_energy%im)
        endif
        if(allocated(gw_cd%self_energy%complx)) then
           deallocate(gw_cd%self_energy%complx)
        endif
        deallocate(gw_cd%self_energy)
      endif
      if(allocated(gw_cd%real_freq)) then
         deallocate(gw_cd%real_freq) 
      endif
      if(allocated(gw_cd%residue_from_freq)) then
         deallocate(gw_cd%residue_from_freq) 
      endif
      if(allocated(gw_cd%spectral_env)) then
        do i = 1, size(gw_cd%spectral_env)
           if(allocated(gw_cd%spectral_env(i)%freq_index)) then
             deallocate(gw_cd%spectral_env(i)%freq_index)
           endif
           if(allocated(gw_cd%spectral_env(i)%freqs)) then
             deallocate(gw_cd%spectral_env(i)%freqs)
           endif
           if(allocated(gw_cd%spectral_env(i)%spectrum)) then
             deallocate(gw_cd%spectral_env(i)%spectrum)
           endif
        enddo
        deallocate(gw_cd%spectral_env)
      endif
     if(allocated(gw_cd%contour_def_start)) then
       deallocate(gw_cd%contour_def_start)
     endif
     if(allocated(gw_cd%contour_def_end)) then
       deallocate(gw_cd%contour_def_end)
     endif
     if(allocated(gw_cd%spin_channel)) then
       deallocate(gw_cd%spin_channel)
     endif
     if(allocated(gw_cd%num_levels)) then
       deallocate(gw_cd%num_levels)
     endif
     if(allocated(gw_cd%corrected_levels)) then
       deallocate(gw_cd%corrected_levels)
     endif
     if(allocated(gw_cd%index_cd_levels)) then
       deallocate(gw_cd%index_cd_levels)
     endif
     if(allocated(gw_cd%sc_env)) then
        do i = 1, size(gw_cd%sc_env)
           if(allocated(gw_cd%sc_env(i)%qp_converged)) then
             deallocate(gw_cd%sc_env(i)%qp_converged)
           endif
           if(allocated(gw_cd%sc_env(i)%converged_last_step)) then
             deallocate(gw_cd%sc_env(i)%converged_last_step)
           endif
        enddo
        deallocate(gw_cd%sc_env)
     endif
     if(allocated(gw_cd%shift)) then
       deallocate(gw_cd%shift)
     endif
     if(allocated(gw_cd%qp_adjust%state)) then
        deallocate(gw_cd%qp_adjust%state)
     endif
     if(allocated(gw_cd%qp_adjust%values)) then
        deallocate(gw_cd%qp_adjust%values)
     endif
           
   end subroutine deallocate_cd_environment
   
end module contour_def_gw_types
