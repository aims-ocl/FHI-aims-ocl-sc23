      module scgw_allocations 

      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use gw_para
      use constants
      use mpi_tasks
      use synchronize_mpi
      use numerical_utilities 
      use scgw_grid
      use poles_fit
      use gt
 

            implicit none

      !allocatable quantities
      real*8,     dimension(:,:,:,:), allocatable :: ovlp_3KS
      real*8,     dimension(:,:,:,:), allocatable :: green_fn_time
      real*8,     dimension(:,:,:,:), allocatable :: green_fn_time_ens !for ensemble calculations 
      real*8,     dimension(:,:,:,:), allocatable :: green_fn_time_dft
      real*8,     dimension(:,:,:,:), allocatable :: green_fn_time_FT
      real*8,     dimension(:,:,:,:), allocatable :: new_green_fn_time
      real*8,     dimension(:,:,:),   allocatable :: green_0_t_0
      real*8,     dimension(:,:),     allocatable :: green_tmp
      real*8,     dimension(:,:),     allocatable :: aux_ovlp
      real*8,     dimension(:,:,:),   allocatable :: polar_green_time
      real*8,     dimension(:,:,:),   allocatable :: polar_green_freq
      real*8,     dimension(:,:),     allocatable :: aux_polar_freq
      real*8,     dimension(:,:,:),   allocatable :: screened_coul_int_freq
      real*8,     dimension(:,:,:),   allocatable :: screened_coul_int_time
      real*8,     dimension(:,:,:,:), allocatable :: self_energy_time
      real*8,     dimension(:,:,:),   allocatable :: xc_matr
      real*8,     dimension(:,:,:),   allocatable :: exchange_self_energy
      real*8,     dimension(:,:,:),   allocatable :: exchange_self_energy0
      real*8,     dimension(:,:),     allocatable :: hartree_pot
      real*8,     dimension(:,:),     allocatable :: hartree_pot0
      real*8,     dimension(:,:),   allocatable :: hartree_dft
      real*8,     dimension(:,:),     allocatable :: aux_overlap_matrix 

      complex*16, dimension(:,:,:,:), allocatable :: self_energy_freq
      complex*16, dimension(:,:,:),   allocatable :: self_energy_omega
      complex*16, dimension(:,:,:,:), allocatable :: green_fn_freq
      complex*16, dimension(:,:,:,:), allocatable :: green_fn_freq_ens !for ensemble calculations 
      complex*16, dimension(:,:,:,:), allocatable :: inv_green_fn_freq

      ! timing
      real*8 fourier, fourier1
      real*8 green_fn_timing, green_fn_timing1
      real*8 time_chi_1, time_chi_0     
      real*8 time_W, time_W1
      real*8 time_self_0, time_self_1 
      real*8 time_inversion_1,time_inversion_0
      real*8 time_3KS_intg
      real*8 time_xc_intg
      real*8 time_transform
      real*8 tot_time_sigma
      real*8 time_analy_continue
      real*8 time_qp_energy
      real*8 rtime
      real*8 time11, time22

      ! other constants
      integer      max_reiteration
      integer      number_reiterations
      real*8       threshold_green
      real*8       average_error, max_error
      real*8       fraction_p 
      !real*8       gev(n_basis)
      real*8       new_chem_pot
      real*8       freq
      real*8       interval
      real*8       n_particles
      real*8       alpha 
      real*8       n_part
      real*8       rpa_c_energy
      real*8       rpa_c_integrand
      !real*8       error_on_G (n_spin)
      real*8       g1
      complex*16   x
      complex*16   selfe
      complex*16   integral_vchi
      character*5  name_of_quantity
      character*15 filename
      character*2  iter
      logical      output
      logical      scgw_converged
      logical      print_spectrum
      logical      get_ensemble_G 
      logical      update_W
      logical      compare_chi_0
      logical      test_FT 

      ! counters
      integer i_state
      integer i_basis
      integer j_basis
      integer i_task
      integer i_spin
      integer i_index
      integer i_basbas
      integer j_basbas
      integer i_freq
      integer i_tau,i
      integer matr_el1, matr_el2
      integer mpierr

      character(*), parameter :: func = 'self_consistent_gw'


      end module scgw_allocations 
