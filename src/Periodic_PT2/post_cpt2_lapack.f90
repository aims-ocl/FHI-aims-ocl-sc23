!****s* FHI-aims/post_scf_correlation_treatment_p0
!  NAME
!   post_scf_correlation_treatment_p0
!  SYNOPSIS

      subroutine post_cpt2_lapack()

!  PURPOSE
!  This subroutine evaluate the necessary matrix elements (3-center overlap,
!  coulomb matrix) used for accurate correlation energy calculations beyond DFT
!  (e.g. MP2, RPA, RPA+, etc).
!
!  USES

      use dimensions
      use species_data
      use runtime_choices
      use pbc_lists, only : k_point_list
      use prodbas
      use hartree_fock
      use hartree_fock_p0
      use physics
      use gw_para
      use lvl_triples
      use tight_binding_auxmat
      use timing
      use mpi_tasks
      use sbt_overlap_aims
      use calculate_fock_matrix_p0, only : cleanup_fock_matrix_calculations, &
          init_fock_matrix_calculations, evaluate_exchange_matr_realspace_p0
      use pbc_lists, only: n_cells
      implicit none

!  ARGUMENTS

!  INPUT
!    none
!  OUTPUT
!    none
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

!  local variables
      character*20 filename
      logical :: need_o3fn
      real*8 :: time_prodbas_add, clock_time_prodbas_add
      real*8 :: rpa_tot_en
      real*8 :: pt2_tot_en
      real*8 :: en_se, en_rse
      real*8 :: exx_ene, d_exx_ene(3) !SVL dummy variables, not actually used anywhere but needed for hf_realspace call 

      integer, allocatable :: n_homo_k(:,:)
      integer :: n_low_state_polar
      integer :: n_k_points_min

!   GW self-energy
      complex*16, allocatable :: gw_selfenergy(:,:,:,:)
      complex*16, allocatable :: gw_selfe_band(:,:,:,:)
!   exact-exhcange term
      real*8, allocatable :: exact_x_kspace(:,:,:)
      real*8, allocatable :: xc_kspace(:,:,:)
      real*8, allocatable :: qp_energy(:,:,:)
      complex*16, allocatable :: coulomb_gamma(:,:)

!   the XC part of the DFT hamiltonian in realspace
      real*8, allocatable :: xc_realspace(:,:)
      real*8, allocatable :: x_realspace(:,:)
      real*8, allocatable :: c_realspace(:,:)
! Counter
      integer i_spin
      integer i_state
      integer i_k_point, i_k_point_local
      integer i_band
      integer i_k_point_index
      integer i_prodbas_1, i_prodbas_2

      integer :: info
      character*150 :: info_str
      character(*), parameter :: &
          func='post_scf_correlation_treatment_p0'

      call get_timestamps(time_prodbas_add, clock_time_prodbas_add)

!  n_irk_points, irk_point_list, irk_point_mapping, inv_irk_point_mapping,
!  irk_point_included (all defined in pbc_lists.f90) will be determined here
      call determine_irreducible_k_grid ()

      if(real_eigenvectors) then
        allocate(KS_eigenvector_irk(n_basis,n_states,n_spin,n_irk_points_task),stat=info) 
        call check_allocation(info, 'KS_eigenvector_irk', func)
        allocate(KS_eigenvector_complex_irk(1,1,1,1),stat=info) 
        call check_allocation(info, 'KS_eigenvector_complex_irk', func)
      else
        allocate(KS_eigenvector_irk(1,1,1,1),stat=info) 
        call check_allocation(info, 'KS_eigenvector_irk', func)
        allocate(KS_eigenvector_complex_irk(n_basis,n_states,n_spin,n_irk_points_task),stat=info) 
        call check_allocation(info, 'KS_eigenvector_complex_irk', func)
      endif

      call distribute_irreducible_eigenvectors &
           ( KS_eigenvector, KS_eigenvector_complex, &
             KS_eigenvector_irk, KS_eigenvector_complex_irk )

      ! --- Ensure ovlp_3fn / coeff_3fn_ten&coulomb_matr_lvl

      ! for non-hartree-fock self-consistent calculations, we need to
      ! construct the product (auxiliary) basis functions, and evaluate
      ! the 3-center overlap and coulomb matrix elements here.
      if (.not.use_hf_kspace) then

        if (.not.allocated(n_homo)) then
         allocate (n_homo(n_spin))
          n_homo(:)=0
        endif

        call initialize_prodbas()

! estimate the largest number of real-space unit cells locally
        if(mod(n_cells, n_tasks).eq.0) then
            n_cells_task = n_cells/n_tasks
        else
            n_cells_task = n_cells/n_tasks+1
        endif

        if(.not.use_hf_realspace)then
           call allocate_hartree_fock_p0()
           if(.not.use_hf_kspace_with_rpa)then
              call cleanup_fock_matrix_calculations
              call init_fock_matrix_calculations(.false., .false.)
              call evaluate_exchange_matr_realspace_p0(KS_eigenvector,KS_eigenvector_complex,occ_numbers, &
                   exx_ene, d_exx_ene, .false., .false.)
              call cleanup_fock_matrix_calculations
           endif
        endif

        n_k_points_min=max(n_k_points_task,1)
        allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_min),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_recip1', func)
        allocate(lvl_tricoeff_recip2(max_n_basbas_sp,n_basis,n_basis),stat=info) 
        call check_allocation(info, 'lvl_tricoeff_recip2', func)
        allocate(kq_point_list(n_k_points,n_k_points),stat=info) 
        call check_allocation(info, 'kq_point_list', func)
        allocate(kpq_point_list(n_k_points,n_k_points),stat=info) 
        call check_allocation(info, 'kpq_point_list', func)
        allocate(coulomb_matr_recip(n_basbas,n_basbas,n_k_points_min),stat=info) 
        call check_allocation(info, 'coulomb_matr_recip', func)
        allocate(coulomb_gamma(n_basbas,n_basbas),stat=info) 
        call check_allocation(info, 'coulomb_gamma', func)

        call initialize_lvl_triples(OVLP_TYPE_COULOMB)
        call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)


        call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)
        call deallocate_tb_auxmat()

        if(use_hse .and. hse_omega_hf /= 0.d0 .and. (.not. use_gw_and_hse)&
            .and. (.not. use_dftpt2_and_hse)) then
           call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
        else
           call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
           !call initialize_tb_auxmat(1, OVLP_TYPE_CUT_ANALYTIC)
!          call initialize_periodic_tb_auxmat(1, 1.d0)
        endif

        call determine_k_minus_q_list(kq_point_list,kpq_point_list)
        call get_coulomb_matr_recip(coulomb_matr_recip,1)


        call evaluate_exchange_matr_kspace_p0 &
             (KS_eigenvalue,KS_eigenvector,KS_eigenvector_complex,occ_numbers)
!<ADDED CODE>
        call deallocate_tb_auxmat()
!</ADDED CODE>

       do i_k_point = 1, n_k_points, 1
        if(myid.ne.mod(i_k_point, n_tasks)) cycle
        i_k_point_local =(i_k_point-1)/n_tasks + 1 

        if(i_k_point .eq. 1 .and. all(abs(k_point_list(i_k_point,:)).lt.1.e-2)) then
          coulomb_gamma(:,:)= coulomb_matr_recip(:,:,i_k_point_local)        
        endif
      enddo

      call deallocate_tb_auxmat()
      call initialize_periodic_tb_auxmat(1, 1.d0)

      call get_coulomb_matr_recip(coulomb_matr_recip,1)
      do i_k_point = 1, n_k_points, 1
       if(myid.ne.mod(i_k_point, n_tasks)) cycle
       i_k_point_local =(i_k_point-1)/n_tasks + 1 
       if(i_k_point .eq. 1 .and. all(abs(k_point_list(i_k_point,:)).lt.1.e-2)) then
                coulomb_matr_recip(:,:, i_k_point_local) = &       
                coulomb_gamma(:, :)

       endif
      enddo

      endif

      call deallocate_tb_auxmat()
      call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)

       if(use_mp2 .or. use_os_mp2 .or. use_dftpt2) then  ! add by igor.
           call pt2_calculation(pt2_tot_en)
       endif

      end subroutine post_cpt2_lapack
!***************
