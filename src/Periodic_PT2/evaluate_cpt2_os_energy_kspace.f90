!****s* FHI-aims/evaluate_cpt2_os_energy_kspace
!  NAME
!   evaluate_cpt2_os_energy_kspace
!  SYNOPSIS

      subroutine evaluate_cpt2_os_energy_kspace &
           ( occ_numbers, &
            pt2_c_energy_os&
           )
      !KS_eigenvalue, &
      !KS_eigenvector, &
      !KS_eigenvector_complex, &

!  PURPOSE
!  Subroutine evaluate_cpt2_os_energy_kspace evaluates the correlation
!  energy at the PT2 level, which is the simplest post-SCF correlation.
!

! USES
      use dimensions
      use prodbas
      use pbc_lists
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      !
      use hartree_fock_p0
      use runtime_choices
      use restart_pt2
      use physics, only: KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue
      use aims_memory_tracking, only : aims_allocate, aims_deallocate
      use localorb_io, only: localorb_info, use_unit
      use geometry, only: species

      implicit none

! ARGUMENTS 

      !integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
!      integer :: n_lumo(n_spin)
!      integer :: n_homo(n_spin)

      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: e_diff
      complex*16  :: E_mp2_tmp_os, E_mp2_tmp_ss
      !real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      !real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      !complex*16  :: KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: E_mp2_a, E_mp2_b

!     output
      real*8  :: pt2_c_energy_os
      ! dump
      real*8  :: pt2_c_energy, pt2_c_energy_ss

! INPUTS
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarisability calcua            ltions
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate. In the present case, n_high_state >= n_homo
!            should be fine. 
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  KS_eigenvector_complex -- complex array,
!            the complex eigenvector of the single-particle calculation,
!            used when "real_eigenvectors == .false."
!           
!
! OUTPUT
! o  pt2_c_energy -- real number, the calculated PT2 correlation energy
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

      real*8  :: pt2_c_energy_local
      real*8  :: pt2_c_energy_start
      real*8  :: pt2_c_energy_local_os
      real*8  :: pt2_c_energy_start_os
      real*8  :: pt2_c_energy_local_ss
      real*8  :: pt2_c_energy_start_ss
      

!    local timing
      real*8  time_pt2(8)
      real*8  time_distribution(6)
      real*8  time_debug(6)

      ! for parallelization
      integer :: i_index, j_index, n_index, m_index
      real*8  :: pt2_term
      real*8  :: pt2_term_os
      real*8  :: pt2_term_ss

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)


      integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
      integer, allocatable :: n_lumo_k(:,:) ! # of non-fully occupied states at a given k-point q
      complex*16, allocatable :: lvl_tricoeff_sum(:,:) ! sum of LVL triple coefficient at k and q points
      complex*16, allocatable :: lvl_tricoeff_sum_all(:,:,:) ! lvl_tricoeff_sum for all babas values
      complex*16, allocatable :: coeff_tmp(:,:) ! temporary triple coefficient, one basis index is transformed to
                                                   ! (occupied) molecular orbital index 
      complex*16, allocatable :: coeff_tmp_2(:,:) ! temporary triple coefficient, two basis indices are transformed to molecular orbital indices.
      complex*16, dimension(:,:), allocatable :: coulomb_tmp_qk
      complex*16, dimension(:,:), allocatable :: coulomb_tmp_qpk
      complex*16, allocatable :: lvl_tricoeff_recip_q(:,:,:)  ! LVL triple coefficients in the q point
      complex*16, allocatable :: lvl_tricoeff_recip_k(:,:,:)  ! LVL triple coefficients in the k point 
      complex*16, allocatable :: lvl_tricoeff_recip_qp(:,:,:) ! LVL triple coefficients in the q' point
      complex*16, allocatable :: lvl_tricoeff_recip_kp(:,:,:) ! LVL triple coefficients in the k' point 
      complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector at a given q point in the BZ
      complex*16, allocatable :: KS_eigenvector_k(:,:,:) ! KS eigenvector at a given k point in the BZ
      complex*16, allocatable :: KS_eigenvector_qp(:,:,:) ! KS eigenvector at a given q point in the BZ
      complex*16, allocatable :: KS_eigenvector_kp(:,:,:) ! KS eigenvector at a given k point in the BZ
      complex*16, allocatable :: lvl_tricoeff_kq(:,:,:,:)   ! Full LVL triple coefficient, one basis index is transformed to
                                                            ! (occupied) molecular orbital index 
      complex*16, allocatable :: lvl_tricoeff_kpqp(:,:,:,:) ! Full LVL triple coefficient, one basis index is transformed to
                                                            ! (occupied) molecular orbital index 
      complex*16, allocatable :: lvl_tricoeff_kpq(:,:,:,:)  ! Full LVL triple coefficient, one basis index is transformed to
                                                            ! (occupied) molecular orbital index 
      complex*16, allocatable :: lvl_tricoeff_kqp(:,:,:,:)  ! Full LVL triple coefficient, one basis index is transformed to
                                                            ! (occupied) molecular orbital index 

      complex*16, allocatable :: lvl_tricoeff_tmp(:,:) ! temporary triple coefficients multiplied by a factor
      complex*16, allocatable :: lvl_tricoeff_tmp_2(:,:) ! temporary triple coefficients


      character(*), parameter :: func = 'evaluate_cpt2_os_energy_kspace.f90'
      integer :: info, mpierr
      character*150 :: info_str

      !==================================================================================
      !
      ! Argument block for parallism
      !
      !==================================================================================
      ! tabulate the required arrays for a batch of tasks (n_tasks in a batch)
      real*8  :: weight_total, weight_total_array(n_tasks)
      integer :: cpt2_para(6,n_tasks), i_tasks(6)
      integer :: i_p_index, i_task, i_task_real
      integer :: cpt2_total_grid
      integer :: cpt2_max_per_thread
      ! additional cost for the parallelism communication
      complex*16, dimension(:,:), allocatable :: coulomb_tmp_qk_para
      complex*16, dimension(:,:), allocatable :: coulomb_tmp_qpk_para
      complex*16, allocatable :: lvl_tricoeff_recip_q_para(:,:,:)  ! LVL triple coefficients in the q point
      complex*16, allocatable :: lvl_tricoeff_recip_k_para(:,:,:)  ! LVL triple coefficients in the k point 
      complex*16, allocatable :: lvl_tricoeff_recip_qp_para(:,:,:) ! LVL triple coefficients in the q' point
      complex*16, allocatable :: lvl_tricoeff_recip_kp_para(:,:,:) ! LVL triple coefficients in the k' point 
      complex*16, allocatable :: KS_eigenvector_q_para(:,:,:) ! KS eigenvector at a given q point in the BZ
      complex*16, allocatable :: KS_eigenvector_k_para(:,:,:) ! KS eigenvector at a given k point in the BZ
      complex*16, allocatable :: KS_eigenvector_qp_para(:,:,:) ! KS eigenvector at a given q point in the BZ
      complex*16, allocatable :: KS_eigenvector_kp_para(:,:,:) ! KS eigenvector at a given k point in the BZ
      integer :: size_coulomb_tmp
      integer :: size_lvl_tricoeff_recip
      integer :: size_KS_eigenvector
      ! request for unblocking MPI communications 
      integer, dimension(MPI_STATUS_SIZE) :: request_o3fn_k  
      integer, dimension(MPI_STATUS_SIZE) :: request_o3fn_kp 
      integer, dimension(MPI_STATUS_SIZE) :: request_o3fn_q  
      integer, dimension(MPI_STATUS_SIZE) :: request_o3fn_qp 
      integer, dimension(MPI_STATUS_SIZE) :: request_eigen_k 
      integer, dimension(MPI_STATUS_SIZE) :: request_eigen_kp
      integer, dimension(MPI_STATUS_SIZE) :: request_eigen_q 
      integer, dimension(MPI_STATUS_SIZE) :: request_eigen_qp
      integer, dimension(MPI_STATUS_SIZE) :: request_c_qk    
      integer, dimension(MPI_STATUS_SIZE) :: request_c_qpk   
      integer, dimension(MPI_STATUS_SIZE) :: status
      ! duplicate required arrays for a vast parallism
      integer :: cpt2_duplicate
      integer :: i_duplicate
      !==================================================================================
      

!     timing

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

!     counters

      ! Three loops of k-mesh: 1) q-k (qk); 2) k; 3) q'(qp)
      integer :: i_qk_point
      integer :: i_qk_point_local
      integer :: i_k_point
      integer :: i_k_point_local
      integer :: i_qp_point
      integer :: i_qp_point_local
      ! We also need other three independent k-mesh arguments
      ! 1) q; 2) k'(kp); 3) q'-k (qpk)
      integer :: i_q_point
      integer :: i_q_point_local
      integer :: i_kp_point
      integer :: i_kp_point_local
      integer :: i_qpk_point
      integer :: i_qpk_point_local
      ! index for basis func, MO(occ, vir), spin, prodbas
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: a_state
      integer :: b_state
      integer :: n_state
      integer :: m_state
      integer :: i_spin
      integer :: i_prodbas_1
      integer :: i_prodbas_2
      integer :: i_atom_1, i_species_1, bboff_1, n_spbb_1
      integer :: i_atom_2, i_species_2, bboff_2, n_spbb_2


      integer :: n_unocc_max
      integer :: n_lumo(n_spin)
      integer :: n_lumo_min
      integer :: n_lumo_eff
      integer :: n_unocc_eff


      complex*16, dimension(n_basbas,1) :: tmp_pt2

      !=====================================================
      ! Numerical test for the squre of the matrix
      !=====================================================
      !complex*16, dimension(2,2) :: test_squre, test_trans
      !real*8, dimension(2) :: test_eigen
      !real*8  :: ev_sqrt
      !integer :: n_nonsingular

      !test_squre(1,1) = (6.1d0,0.0d0)
      !test_squre(1,2) = (5.65685d0,2.26274d0)
      !test_squre(2,1) = (5.65685d0,-2.26274d0)
      !test_squre(2,2) = (6.1d0,0.0d0)
      !if (myid .eq. 0) then
      !    test_squre = -test_squre
      !    call diagonalize_auxmat_lapack_complex(2,test_squre,safe_minimum,0.0d0,&
      !        n_nonsingular, test_eigen, test_trans,'')
      !    test_eigen = - test_eigen
      !    write(use_unit,*) test_eigen
      !    do i_basis_1 = 1, n_nonsingular
      !       ev_sqrt = sqrt(test_eigen(i_basis_1))
      !       do i_basis_2 = 1, 2, 1
      !          test_trans(i_basis_2, i_basis_1) = &
      !          test_trans(i_basis_2,i_basis_1) * ev_sqrt**(0.5)
      !       enddo
      !    enddo
      !    test_squre = (0.d0,0.d0)
      !    call zgemm('N', 'C', 2, 2, n_nonsingular, &
      !    &          (1.0d0,0.d0), test_trans(:,1:n_nonsingular), 2, &
      !    &                 test_trans(:,1:n_nonsingular), 2, &
      !    &           (0.d0,0.d0), test_squre, 2)
      !    write(use_unit,*) test_squre
      !endif
      !=====================================================
      !complex*16, dimension(n_basbas,n_basbas) :: test_squre, test_trans
      !real*8, dimension(n_basbas) :: test_eigen
      !real*8  :: ev_sqrt
      !integer :: n_nonsingular



!     begin work

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-----------------------------------------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the periodic osPT2 correlation energy  ... "
      endif

      if(flag_KS_eigenfunc_conjg) then
         KS_eigenvector_complex = conjg(KS_eigenvector_complex)
      endif
!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin
          do a_state = 1, n_states
           if (abs(occ_numbers(a_state,i_spin,i_k_point)-dble(2/n_spin)) &
                           .lt.1.d-8) then
             n_first(i_spin)= a_state + 1
           endif
          enddo
          if(n_first(i_spin) .gt. n_states) then
           n_first(i_spin) = n_states
          endif
        enddo
      enddo

      allocate(n_homo_k(n_k_points,n_spin),stat=info) 
      call check_allocation(info, 'n_homo_k', func)
      allocate(n_lumo_k(n_k_points,n_spin),stat=info) 
      call check_allocation(info, 'n_unocc_k', func)

      n_homo_k(:,:) = 0
      n_lumo_k(:,:)= n_states
      do i_spin = 1, n_spin, 1
        do i_k_point = 1, n_k_points, 1
          do a_state = 1, n_states
            if(occ_numbers(a_state,i_spin,i_k_point) .gt. 1.e-12) then
             n_homo_k(i_k_point,i_spin) = a_state
            endif
          enddo
         
          do a_state = n_states, 1, -1
            if(occ_numbers(a_state,i_spin,i_k_point) .lt. 1.d0) then
             n_lumo_k(i_k_point,i_spin) = a_state 
            endif
          enddo
!          if(myid.eq.0) then
!            write(use_unit,*) i_spin, i_k_point, n_homo_k(i_k_point, i_spin), n_lumo_k(i_k_point,i_spin)
!          endif
        enddo
      enddo
      n_homo(1)=maxval(n_homo_k(:,1), n_k_points)
      n_homo(n_spin)=maxval(n_homo_k(:,n_spin), n_k_points)
      n_homo_max = max(n_homo(1), n_homo(n_spin)) 

      n_lumo(1)=minval(n_lumo_k(:,1), n_k_points)
      n_lumo(n_spin)=minval(n_lumo_k(:,n_spin), n_k_points)
      n_lumo_min = min(n_lumo(1), n_lumo(n_spin)) 
      n_unocc_max = n_states - n_lumo_min + 1

      ! coulomb matrices for q-k and q'-k
      allocate(coulomb_tmp_qk(n_basbas, n_basbas),stat=info)
      call check_allocation(info, 'coulomb_tmp_qk                  ')
      allocate(coulomb_tmp_qpk(n_basbas, n_basbas),stat=info)
      call check_allocation(info, 'coulomb_tmp_qpk                 ')
      ! tricoeff for four k-points, including k, q, k', and q'
      allocate(lvl_tricoeff_recip_k(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_k', func)
      allocate(lvl_tricoeff_recip_q(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_q', func)
      allocate(lvl_tricoeff_recip_kp(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_kp', func)
      allocate(lvl_tricoeff_recip_qp(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_qp', func)
      ! KS eigenvector matrices for four k-points, including k, q, k', and q'
      allocate(KS_eigenvector_k(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_k', func)
      allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_q', func)
      allocate(KS_eigenvector_kp(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_kp', func)
      allocate(KS_eigenvector_qp(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_qp', func)
      ! molecular-orbital tricoeff matrices for four crossing of k-points including
      ! M(kq), M(k'q'), M(kq'), and M(k'q) 
      allocate(lvl_tricoeff_kq(n_basbas,n_states,n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_kq', func)
      allocate(lvl_tricoeff_kpqp(n_basbas,n_states,n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_kpqp', func)
      allocate(lvl_tricoeff_kqp(n_basbas,n_states,n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_kqp', func)
      allocate(lvl_tricoeff_kpq(n_basbas,n_states,n_homo_max,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_kpq', func)
      ! temp matrices
      allocate(lvl_tricoeff_sum(n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_sum', func)
      allocate(lvl_tricoeff_sum_all(n_basbas,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_sum_all', func)
      allocate(coeff_tmp(n_basis,n_homo_max),stat=info) 
      call check_allocation(info, 'coeff_tmp', func)
      allocate(coeff_tmp_2(n_states,n_homo_max),stat=info) 
      call check_allocation(info, 'coeff_tmp_2', func)
      allocate(lvl_tricoeff_tmp(n_basbas,n_unocc_max),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_tmp', func)
      allocate(lvl_tricoeff_tmp_2(n_basbas,n_unocc_max),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_tmp_2', func)

      ! For parallelization use
      ! coulomb matrices for q-k and q'-k
      allocate(coulomb_tmp_qk_para(n_basbas, n_basbas),stat=info)
      call check_allocation(info, 'coulomb_tmp_qk_para      ')
      allocate(coulomb_tmp_qpk_para(n_basbas, n_basbas),stat=info)
      call check_allocation(info, 'coulomb_tmp_qpk_para     ')
      size_coulomb_tmp = size(coulomb_tmp_qk_para)
      ! tricoeff for four k-points, including k, q, k', and q'
      allocate(lvl_tricoeff_recip_k_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_k_para', func)
      allocate(lvl_tricoeff_recip_q_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_q_para', func)
      allocate(lvl_tricoeff_recip_kp_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_kp_para', func)
      allocate(lvl_tricoeff_recip_qp_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_qp_para', func)
      size_lvl_tricoeff_recip = size(lvl_tricoeff_recip_k_para)
      ! KS eigenvector matrices for four k-points, including k, q, k', and q'
      allocate(KS_eigenvector_k_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_k_para', func)
      allocate(KS_eigenvector_q_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_q_para', func)
      allocate(KS_eigenvector_kp_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_kp_para', func)
      allocate(KS_eigenvector_qp_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_qp_para', func)
      size_KS_eigenvector = size(KS_eigenvector_k_para)



      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(n_low_state)
      else
          n_low_state = 1
      endif

      call get_timestamps(time_pt2(1), time_pt2(2) )

      ! =====================================================
      ! calculate the square root of coulomb_matr_recip
      ! and replace the original array "coulomb_matr_recip
      do i_k_point = 1, n_k_points, 1
        if (myid .eq. mod(i_k_point, n_tasks)) then
            i_k_point_local = (i_k_point-1)/n_tasks + 1
            call power_auxmat_lapack_complex &
                (coulomb_matr_recip(:,:,i_k_point_local), 0.5d0, '')
        endif
      enddo

      call get_timestamps(time_pt2(3),time_pt2(4))
      if (myid .eq. 0) then
          write(use_unit, "(2X,'Timing for V^{0.5} of coulomb_matr_recip     :',1X,F12.3,' s',4X,F12.3,' s')") &
              time_pt2(3)-time_pt2(1), &
              time_pt2(4)-time_pt2(2)
      endif
      ! =====================================================
      ! Duplicating the required arrays to free threads,
      ! which can improve the communication efficiency 
      ! in a vast parallelization
      !
      ! Sketch map:
      ! 
      ! <--n_k_points-->                             <--A1-->
      ! |-----1st------|------2nd-----|------3rd-----|------|
      ! |---------------------------------------------------|
      ! <----------------------n_tasks---------------------->
      !
      ! A1      : The rest threads = mod(n_tasks, n_k_points)
      !
      ! Purpose : Copy the related arrays in the 1st block
      !           of threads to the other blocks.
      !           i.e. "2nd", "3rd", and "A1" in this map
      !           Then the PT2 calculation in each thread only
      !           communicate the required arrays within the
      !           relevant block.
      !           For the tasks in "A1" block, the communication
      !           happens in "3rd" and "A1" blocks.
      cpt2_duplicate = (n_tasks - 1)/n_k_points + 1 ! in the map, cpt2_duplicate = 4
      do i_task = n_k_points + 1, n_tasks, 1
        ! determine the index of the target thread
        if (i_task .eq. n_tasks) then
            i_task_real = 0
        else
            i_task_real = i_task
        endif
        ! determine the index of the original thread
        i_tasks(1) = mod(i_task, n_k_points) 
        if (i_tasks(1) .eq. 0) i_tasks(1) = n_k_points
        ! send the required arrays from the orginal thread to the target
        if (myid .eq. i_tasks(1)) then
            call mpi_send(n_k_points_task, 1, MPI_INTEGER, i_task_real, 0, mpi_comm_global, info)
            call send_complex_vector(lvl_tricoeff_recip1, size(lvl_tricoeff_recip1),i_task_real)
            call send_complex_vector(coulomb_matr_recip, size(coulomb_matr_recip),i_task_real)
            if(real_eigenvectors) then
               call send_real_vector(KS_eigenvector, size(KS_eigenvector),i_task_real)
            else
               call send_complex_vector(KS_eigenvector_complex, size(KS_eigenvector_complex),i_task_real)
            endif
        endif
        if (myid .eq. i_task_real) then
            call mpi_recv(n_k_points_task, 1, MPI_INTEGER, i_tasks(1), 0, mpi_comm_global, status, info)
            if (allocated(lvl_tricoeff_recip1)) then
                deallocate(lvl_tricoeff_recip1)
            endif
            allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_task),stat=info)
            if (allocated(coulomb_matr_recip)) then
                deallocate(coulomb_matr_recip)
            endif
            allocate(coulomb_matr_recip(n_basbas, n_basbas,n_k_points_task),stat=info)
            if(real_eigenvectors) then
                if (allocated(KS_eigenvector)) then
                    call aims_deallocate(KS_eigenvector, "KS_eigenvector" )
                endif
                call aims_allocate(KS_eigenvector,n_basis,n_states,n_spin,n_k_points_task,"KS_eigenvector")
            else
                if (allocated(KS_eigenvector_complex)) then
                    call aims_deallocate(KS_eigenvector_complex,"KS_eigenvector_complex")
                endif
               call aims_allocate(KS_eigenvector_complex,n_basis,n_states,n_spin,n_k_points_task,"KS_eigenvector_complex")
            endif
            call receive_complex_vector(lvl_tricoeff_recip1, size(lvl_tricoeff_recip1),i_tasks(1))
            call receive_complex_vector(coulomb_matr_recip, size(coulomb_matr_recip),i_tasks(1))
            if(real_eigenvectors) then
                call receive_real_vector(KS_eigenvector, size(KS_eigenvector),i_tasks(1))
            else
                call receive_complex_vector(KS_eigenvector_complex, size(KS_eigenvector_complex),i_tasks(1))
            endif
        endif
      enddo
      ! =====================================================


      cpt2_total_grid = n_irk_points * n_k_points * n_k_points
      cpt2_max_per_thread = (cpt2_total_grid-1) / n_tasks + 1

      time_distribution       = 0.d0
      call get_timestamps(time_pt2(1), time_pt2(2))

      ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
      ! pt2_c_energy_local :: the PT2 correlation in each thread
      ! pt2_c_energy       :: the final PT2 correlation
      ! os :: opposite-spin pair
      ! ss :: parallel-spin pair
      pt2_c_energy_start      = 0.d0 
      pt2_c_energy_local      = 0.d0 
      pt2_c_energy            = 0.d0 
      pt2_c_energy_start_os   = 0.d0 
      pt2_c_energy_local_os   = 0.d0 
      pt2_c_energy_os         = 0.d0 
      pt2_c_energy_start_ss   = 0.d0 
      pt2_c_energy_local_ss   = 0.d0 
      pt2_c_energy_ss         = 0.d0 
      ! =====================================================
      ! reading the restart info., Igor
      ! the story behind the following five lines is:
      ! 1) in the 0th thread, loading the pt2_c_energy_start for the restart file, 
      !       if restart_pt2_read = .true.
      ! 2) broadcasting pt2_c_energy_start to all threads from the 0th thread
      ! 3) broadcasting pt2_finish to all threads from the 0th thread
      pt2_finish              = 0
      current_pt2             = 0.0d0
      current_pt2_os          = 0.0d0
      current_pt2_ss          = 0.0d0
      call read_restart_pt2_info()
      pt2_c_energy_start      = current_pt2    ! only make sence for myid.eq.0 if restart_pt2_read = .true.
      pt2_c_energy_start_os   = current_pt2_os
      pt2_c_energy_start_ss   = current_pt2_ss
      call bcast_real(pt2_c_energy_start,0)
      call bcast_real(pt2_c_energy_start_os,0)
      call bcast_real(pt2_c_energy_start_ss,0)
      call mp_bcast_int(pt2_finish)
      ! =====================================================

      do i_p_index = pt2_finish+1, cpt2_max_per_thread, 1

        call get_timestamps(time_pt2(3),time_pt2(4))
        call get_timestamps(time_distribution(3),time_distribution(4))
        !
        ! tabulate the memory distribution for this round
        !
        call initialize_cpt2_para(cpt2_para, weight_total_array, i_p_index)

        !
        ! Here, we start the distribution from the thread with largest index,
        !
        do i_task = n_tasks, 1, -1
        

            if(cpt2_para(1,i_task) .eq. 0) cycle

        
            !i_k_point   = cpt2_para(1,i_task)
            !i_kp_point  = cpt2_para(2,i_task)
            !i_q_point   = cpt2_para(3,i_task)
            !i_qp_point  = cpt2_para(4,i_task)
            !i_qk_point  = cpt2_para(5,i_task)
            !i_qpk_point = cpt2_para(6,i_task)

            !===============================================================
            ! The PT2 calculation in each task only communicate the 
            ! required arrays within the relevant block.
            ! For the tasks in "A1" block, the communication
            ! happens between the last second and "A1" blocks.
            !===============================================================

            ! i_duplicate+1 : determine the index of the task block
            i_duplicate = (i_task - 1)/n_k_points

            do i_index = 1, 6, 1

              !! MARK here by igor
              !i_task_real = mod(cpt2_para(i_index,i_task), n_tasks)
              !if (i_task_real .eq. 0 .and. n_tasks .gt. n_k_points) i_task_real = n_k_points
              i_tasks(i_index) = mod(cpt2_para(i_index,i_task), n_tasks) + i_duplicate * n_k_points

              if (i_tasks(i_index) .gt. n_tasks) then ! if the task out of "A1" block, using the last second block
                  i_tasks(i_index) = i_tasks(i_index) - n_k_points
              else if (i_tasks(i_index) .eq. n_tasks) then ! if the task is "n_tasks", using myid=0
                  i_tasks(i_index) = 0
              endif

            enddo
            !if (myid .eq. 0) then
            !    write(use_unit,'(6I6)') cpt2_para
            !    write(use_unit,'(6I6)') i_tasks
            !endif
            !===============================================================

            i_k_point_local   = (cpt2_para(1,i_task)-1)/n_tasks + 1
            i_kp_point_local  = (cpt2_para(2,i_task)-1)/n_tasks + 1
            i_q_point_local   = (cpt2_para(3,i_task)-1)/n_tasks + 1
            i_qp_point_local  = (cpt2_para(4,i_task)-1)/n_tasks + 1
            i_qk_point_local  = (cpt2_para(5,i_task)-1)/n_tasks + 1
            i_qpk_point_local = (cpt2_para(6,i_task)-1)/n_tasks + 1

            ! loading arrays 
            if (myid.eq.i_tasks(1)) then ! for i_k_point
                lvl_tricoeff_recip_k_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_k_point_local)
                if(real_eigenvectors) then
                   KS_eigenvector_k_para(:,:,:)=KS_eigenvector(:,:,:,i_k_point_local)
                else
                   KS_eigenvector_k_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_k_point_local)
                endif
            endif
            if (myid.eq.i_tasks(2)) then ! for i_kp_point
                lvl_tricoeff_recip_kp_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_kp_point_local)
                if(real_eigenvectors) then
                   KS_eigenvector_kp_para(:,:,:)=KS_eigenvector(:,:,:,i_kp_point_local)
                else
                   KS_eigenvector_kp_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_kp_point_local)
                endif
            endif
            if (myid.eq.i_tasks(3)) then ! for i_q_point
                lvl_tricoeff_recip_q_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_q_point_local)
                if(real_eigenvectors) then
                   KS_eigenvector_q_para(:,:,:)=KS_eigenvector(:,:,:,i_q_point_local)
                else
                   KS_eigenvector_q_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_q_point_local)
                endif
            endif
            if (myid.eq.i_tasks(4)) then ! for i_qp_point
                lvl_tricoeff_recip_qp_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_qp_point_local)
                if(real_eigenvectors) then
                   KS_eigenvector_qp_para(:,:,:)=KS_eigenvector(:,:,:,i_qp_point_local)
                else
                   KS_eigenvector_qp_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_qp_point_local)
                endif
            endif
            if (myid.eq.i_tasks(5)) then ! for i_qk_point
                coulomb_tmp_qk_para(:,:) = coulomb_matr_recip(:,:,i_qk_point_local)
            endif
            if (myid.eq.i_tasks(6)) then ! for i_qpk_point
                coulomb_tmp_qpk_para(:,:) = coulomb_matr_recip(:,:,i_qpk_point_local)
            endif

            ! distribute relevant arrays to the target "i_task"
            if (i_task .eq. n_tasks) then
                i_task_real = 0
            else
                i_task_real = i_task
            endif

            if (myid.ne.i_task_real) then
              if (myid .eq. i_tasks(1)) then ! for i_k_point
                call isend_complex_vector(lvl_tricoeff_recip_k_para,&
                                            size_lvl_tricoeff_recip,i_task_real, &
                                            request_o3fn_k)
                call isend_complex_vector(KS_eigenvector_k_para,&
                                            size_KS_eigenvector,i_task_real, &
                                            request_eigen_k)
              endif
              if (myid .eq. i_tasks(2)) then ! for i_kp_point
                call isend_complex_vector(lvl_tricoeff_recip_kp_para,&
                                            size_lvl_tricoeff_recip,i_task_real, &
                                            request_o3fn_kp)
                call isend_complex_vector(KS_eigenvector_kp_para,&
                                            size_KS_eigenvector,i_task_real, &
                                            request_eigen_kp)
              endif
              if (myid .eq. i_tasks(3)) then ! for i_q_point
                call isend_complex_vector(lvl_tricoeff_recip_q_para,&
                                            size_lvl_tricoeff_recip,i_task_real, &
                                            request_o3fn_q)
                call isend_complex_vector(KS_eigenvector_q_para,&
                                            size_KS_eigenvector,i_task_real, &
                                            request_eigen_q)
              endif
              if (myid .eq. i_tasks(4)) then ! for i_qp_point
                call isend_complex_vector(lvl_tricoeff_recip_qp_para,&
                                            size_lvl_tricoeff_recip,i_task_real, &
                                            request_o3fn_qp)
                call isend_complex_vector(KS_eigenvector_qp_para,&
                                            size_KS_eigenvector,i_task_real, &
                                            request_eigen_qp)
              endif
              if (myid .eq. i_tasks(5)) then ! for i_qk_point
                call isend_complex_vector(coulomb_tmp_qk_para,&
                                            size_coulomb_tmp,i_task_real, &
                                            request_c_qk)
              endif
              if (myid .eq. i_tasks(6)) then ! for i_qpk_point
                call isend_complex_vector(coulomb_tmp_qpk_para,&
                                            size_coulomb_tmp,i_task_real, &
                                            request_c_qpk)
              endif
              ! now wait until the communications complete in related threads
              if (myid .eq. i_tasks(1)) then ! for i_k_point
                call isend_wait(request_o3fn_k)
                call isend_wait(request_eigen_k)
              endif
              if (myid .eq. i_tasks(2)) then ! for i_kp_point
                call isend_wait(request_o3fn_kp)
                call isend_wait(request_eigen_kp)
              endif
              if (myid .eq. i_tasks(3)) then ! for i_q_point
                call isend_wait(request_o3fn_q)
                call isend_wait(request_eigen_q)
              endif
              if (myid .eq. i_tasks(4)) then ! for i_qp_point
                call isend_wait(request_o3fn_qp)
                call isend_wait(request_eigen_qp)
              endif
              if (myid .eq. i_tasks(5)) then ! for i_qk_point
                call isend_wait(request_c_qk)
              endif
              if (myid .eq. i_tasks(6)) then ! for i_qpk_point
                call isend_wait(request_c_qpk)
              endif
          else ! if the thread relates to this task
              if (myid .ne. i_tasks(1)) then ! for i_k_point
                call receive_complex_vector(lvl_tricoeff_recip_k_para,&
                                            size_lvl_tricoeff_recip,i_tasks(1))
                call receive_complex_vector(KS_eigenvector_k_para,&
                                            size_KS_eigenvector,i_tasks(1))
              endif
              if (myid .ne. i_tasks(2)) then ! for i_kp_point
                call receive_complex_vector(lvl_tricoeff_recip_kp_para,&
                                            size_lvl_tricoeff_recip,i_tasks(2))
                call receive_complex_vector(KS_eigenvector_kp_para,&
                                            size_KS_eigenvector,i_tasks(2))
              endif
              if (myid .ne. i_tasks(3)) then ! for i_q_point
                call receive_complex_vector(lvl_tricoeff_recip_q_para,&
                                            size_lvl_tricoeff_recip,i_tasks(3))
                call receive_complex_vector(KS_eigenvector_q_para,&
                                            size_KS_eigenvector,i_tasks(3))
              endif
              if (myid .ne. i_tasks(4)) then ! for i_qp_point
                call receive_complex_vector(lvl_tricoeff_recip_qp_para,&
                                            size_lvl_tricoeff_recip,i_tasks(4))
                call receive_complex_vector(KS_eigenvector_qp_para,&
                                            size_KS_eigenvector,i_tasks(4))
              endif
              if (myid .ne. i_tasks(5)) then ! for i_qk_point
                call receive_complex_vector(coulomb_tmp_qk_para,&
                                            size_coulomb_tmp,i_tasks(5))
              endif
              if (myid .ne. i_tasks(6)) then ! for i_qpk_point
                call receive_complex_vector(coulomb_tmp_qpk_para,&
                                            size_coulomb_tmp,i_tasks(6))
              endif

              ! Finalize the memory distribution
              lvl_tricoeff_recip_k(:,:,:)  = lvl_tricoeff_recip_k_para(:,:,:)
              lvl_tricoeff_recip_kp(:,:,:) = lvl_tricoeff_recip_kp_para(:,:,:)
              lvl_tricoeff_recip_q(:,:,:)  = lvl_tricoeff_recip_q_para(:,:,:)
              lvl_tricoeff_recip_qp(:,:,:) = lvl_tricoeff_recip_qp_para(:,:,:)
              KS_eigenvector_k(:,:,:)      = KS_eigenvector_k_para(:,:,:)
              KS_eigenvector_kp(:,:,:)     = KS_eigenvector_kp_para(:,:,:)
              KS_eigenvector_q(:,:,:)      = KS_eigenvector_q_para(:,:,:)
              KS_eigenvector_qp(:,:,:)     = KS_eigenvector_qp_para(:,:,:)
              coulomb_tmp_qk(:,:)          = coulomb_tmp_qk_para(:,:)
              coulomb_tmp_qpk(:,:)         = coulomb_tmp_qpk_para(:,:)
              weight_total                 = weight_total_array(i_task)
              i_k_point                    = cpt2_para(1,i_task)
              i_kp_point                   = cpt2_para(2,i_task)
              i_q_point                    = cpt2_para(3,i_task)
              i_qp_point                   = cpt2_para(4,i_task)
              i_qk_point                   = cpt2_para(5,i_task)
              i_qpk_point                  = cpt2_para(6,i_task)
            endif

            !call mpi_barrier(mpi_comm_global,info)

        enddo ! i_task

        !call mpi_barrier(mpi_comm_global,info)

        !============
        ! Analysis the timing for the memory distribution.
        !============
        call get_timestamps(time_distribution(5), time_distribution(6))
        time_distribution(3) = time_distribution(5) - time_distribution(3)
        time_distribution(4) = time_distribution(6) - time_distribution(4)
        time_distribution(1) = time_distribution(1) + time_distribution(3)
        time_distribution(2) = time_distribution(2) + time_distribution(4)
        call sync_timing(time_distribution(3))
        call sync_timing(time_distribution(4))

        if (myid .eq. 0) then
            i_task_real = n_tasks
        else
            i_task_real = myid
        endif

        if (cpt2_para(1,i_task_real) .ne. 0) then

            call get_timestamps(time_debug(1),time_debug(2))

            ! ==================
            ! Inverse square root of the coulomb matrix
            !     coulomb_matr_qk
            !     coulomb_matr_qpk
            ! ==================
            ! Version 01
            ! =================
            !call power_auxmat_lapack_complex(coulomb_tmp_qk, 0.5d0, '')
            !call power_auxmat_lapack_complex(coulomb_tmp_qpk, 0.5d0, '')
            !coulomb_tmp_qk = conjg(coulomb_tmp_qk)
            !coulomb_tmp_qpk = conjg(coulomb_tmp_qpk)
            ! ===================
            ! Version 02
            ! ==================
            !call power_genmat_lapack_complex(n_basbas, coulomb_tmp_qk, 0.5d0,&
            !    safe_minimum, 0.0d0, '')
            !call power_genmat_lapack_complex(n_basbas, coulomb_tmp_qpk, 0.5d0,&
            !    safe_minimum, 0.0d0, '')

            !call get_timestamps(time_debug(3),time_debug(4))
            !if (myid .eq. 1) then
            !    write(use_unit, '(A35,2X,2F16.3)') "Timing for V^{0.5} of qk and qpk", &
            !        time_debug(3)-time_debug(1), &
            !        time_debug(4)-time_debug(2)
            !endif

            ! ==================
            ! for M(k,q)
            ! ==================
            !call get_timestamps(time_debug(1),time_debug(2))
            lvl_tricoeff_sum_all     = (0.d0,0.d0)
            do i_basis_1 = 1, n_basis, 1
              do i_basis_2 = 1, n_basis, 1
                i_atom_1 = Cbasis_to_atom(i_basis_1)
                i_species_1 = species(i_atom_1)
                bboff_1 = atom2basbas_off(i_atom_1)
                n_spbb_1 = sp2n_basbas_sp(i_species_1)
    
                i_atom_2 = Cbasis_to_atom(i_basis_2)
                i_species_2 = species(i_atom_2)
                bboff_2 = atom2basbas_off(i_atom_2)
                n_spbb_2 = sp2n_basbas_sp(i_species_2)
    
                ! in lvl_tricoeff_recip_q the q vector is associated with the second basis i_basis_2,
                ! for convenince we change this in lvl_tricoeff_sum_all
                lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,i_basis_2,i_basis_1) = &
                lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,i_basis_2,i_basis_1) + &
                  conjg(lvl_tricoeff_recip_q(1:n_spbb_1, i_basis_1, i_basis_2)) + &
                  lvl_tricoeff_recip2(1:n_spbb_1, i_basis_2, i_basis_1)
    
                lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) = &
                lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) + &
                  lvl_tricoeff_recip_k(1:n_spbb_2, i_basis_2, i_basis_1) + &
                  lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)
             enddo
            enddo

            ! multiply the 3-center overlap matrix with the squre root
            ! of the complex coulomb matrix
            call get_v_multi_ovlp3fn_complex(coulomb_tmp_qk, lvl_tricoeff_sum_all,1)
            
            lvl_tricoeff_kq(:,:,:,:) = (0.d0,0.d0)
            do i_prodbas_1 = 1, n_basbas, 1
              do i_basis_1 = 1, n_basis, 1
                do i_basis_2 = 1, n_basis, 1
                  lvl_tricoeff_sum(i_basis_2, i_basis_1) = lvl_tricoeff_sum_all(i_prodbas_1, i_basis_2, i_basis_1)
                !  end loop over i_basis_2
                enddo
              !  end loop over i_basis_1
              enddo
    
              do i_spin = 1, n_spin, 1
                coeff_tmp(:,:) = (0.d0,0.d0)
                call zgemm('N', 'N', n_basis, n_homo_max, n_basis, (1.d0,0.d0), &
                           lvl_tricoeff_sum, n_basis, &
                           KS_eigenvector_k(1,1,i_spin), n_basis, (0.d0,0.d0), &
                           coeff_tmp, n_basis)
    
                coeff_tmp_2(:,:) = (0.d0,0.d0)
                call zgemm('C', 'N', n_states, n_homo_max, n_basis, (1.d0,0.d0), &
                           KS_eigenvector_q(1,1,i_spin), n_basis,  &
                           coeff_tmp, n_basis, (0.d0,0.d0), &
                           coeff_tmp_2, n_states)

                lvl_tricoeff_kq(i_prodbas_1,:,:,i_spin) = coeff_tmp_2(:,:)
              ! end loop over i_spin
              enddo
            ! end loop over i_prodbas_1
            enddo
            !call get_timestamps(time_debug(3),time_debug(4))
            !if (myid .eq. 1) then
            !    write(use_unit, '(2X,A35,2X,2F16.3)') "Timing for M(k,q) generation:", &
            !        time_debug(3)-time_debug(1), &
            !        time_debug(4)-time_debug(2)
            !endif
            ! ==================
            ! for M(kp,qp)
            ! ==================
            lvl_tricoeff_sum_all = (0.d0,0.d0)
            do i_basis_1 = 1, n_basis, 1
              do i_basis_2 = 1, n_basis, 1
                i_atom_1 = Cbasis_to_atom(i_basis_1)
                i_species_1 = species(i_atom_1)
                bboff_1 = atom2basbas_off(i_atom_1)
                n_spbb_1 = sp2n_basbas_sp(i_species_1)
    
                i_atom_2 = Cbasis_to_atom(i_basis_2)
                i_species_2 = species(i_atom_2)
                bboff_2 = atom2basbas_off(i_atom_2)
                n_spbb_2 = sp2n_basbas_sp(i_species_2)
    
                ! in lvl_tricoeff_recip_q the q vector is associated with the second basis i_basis_2,
                ! for convenince we change this in lvl_tricoeff_sum_all
                lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,i_basis_2,i_basis_1) = &
                lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,i_basis_2,i_basis_1) + &
                  conjg(lvl_tricoeff_recip_qp(1:n_spbb_1, i_basis_1, i_basis_2)) + &
                  lvl_tricoeff_recip2(1:n_spbb_1, i_basis_2, i_basis_1)
    
                lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) = &
                lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) + &
                  lvl_tricoeff_recip_kp(1:n_spbb_2, i_basis_2, i_basis_1) + &
                  lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)
              enddo
            enddo

            ! multiply the 3-center overlap matrix with the inverse squre root
            ! of the coulomb matrix
            call get_v_multi_ovlp3fn_complex(coulomb_tmp_qk, lvl_tricoeff_sum_all,2)
            
            lvl_tricoeff_kpqp(:,:,:,:) = (0.d0,0.d0)
            do i_prodbas_1 = 1, n_basbas, 1
              do i_basis_1 = 1, n_basis, 1
                do i_basis_2 = 1, n_basis, 1
                  lvl_tricoeff_sum(i_basis_2, i_basis_1) = lvl_tricoeff_sum_all(i_prodbas_1, i_basis_2, i_basis_1)
                !  end loop over i_basis_2
                enddo
              !  end loop over i_basis_1
              enddo
    
              do i_spin = 1, n_spin, 1
                coeff_tmp(:,:) = (0.d0,0.d0)
                call zgemm('N', 'N', n_basis, n_homo_max, n_basis, (1.d0,0.d0), &
                           lvl_tricoeff_sum, n_basis, &
                           KS_eigenvector_kp(1,1,i_spin), n_basis, (0.d0,0.d0), &
                           coeff_tmp, n_basis)
    
                coeff_tmp_2(:,:) = (0.d0,0.d0)
                call zgemm('C', 'N', n_states, n_homo_max, n_basis, (1.d0,0.d0), &
                           KS_eigenvector_qp(1,1,i_spin), n_basis,  &
                           coeff_tmp, n_basis, (0.d0,0.d0), &
                           coeff_tmp_2, n_states)
    
                lvl_tricoeff_kpqp(i_prodbas_1,:,:,i_spin) = coeff_tmp_2(:,:)
              ! end loop over i_spin
              enddo
            ! end loop over i_prodbas_1
            enddo

            !call get_timestamps(time_debug(1),time_debug(2))
            ! Now do PT2 calculation
            pt2_term    = 0.0d0
            pt2_term_os = 0.0d0
            pt2_term_ss = 0.0d0
            do a_state = n_low_state, n_homo_k(i_k_point, 1), 1
                do b_state = n_low_state, n_homo_k(i_kp_point,1), 1
                    do n_state = n_lumo_k(i_q_point,1), n_states_k(i_q_point), 1
                        do m_state = n_lumo_k(i_qp_point,1), n_states_k(i_qp_point),1
                            e_diff = KS_eigenvalue(a_state,1,i_k_point)  &
                                   + KS_eigenvalue(b_state,1,i_kp_point) &
                                   - KS_eigenvalue(n_state,1,i_q_point)  &
                                   - KS_eigenvalue(m_state,1,i_qp_point)
                            if (abs(e_diff).lt.1e-6)then
                               write(use_unit,'(10X,A)') &
                                    "****************************************"
                               write(use_unit,'(10X,2A)') "| Warning :", &
                                 " too close to degeneracy"
                               write(use_unit,'(10X,A)') &
                                    "****************************************"
                            endif
                            E_mp2_a = 0.0d0
                            !E_mp2_b = 0.0d0
                            do i_prodbas_1 = 1, n_basbas, 1
                                E_mp2_a = E_mp2_a + &
                                  lvl_tricoeff_kq(i_prodbas_1,n_state,a_state,1) * &
                                  lvl_tricoeff_kpqp(i_prodbas_1,m_state,b_state,1)
                                !E_mp2_b = E_mp2_b + &
                                !  lvl_tricoeff_kqp(i_prodbas_1,m_state,a_state,1) * &
                                !  lvl_tricoeff_kpq(i_prodbas_1,n_state,b_state,1)
                            enddo  ! i_prodbas_1
                            E_mp2_tmp_os = conjg(E_mp2_a)*E_mp2_a
                            E_mp2_tmp_os = E_mp2_tmp_os / e_diff

                            !E_mp2_tmp_ss = conjg(E_mp2_a)*(E_mp2_a - E_mp2_b)
                            !E_mp2_tmp_ss = E_mp2_tmp_ss / e_diff

                            !if (i_p_index .eq. 7 .and. myid .eq. 0) then
                            !    write(use_unit,'(4I2,f10.3,f12.6,f12.6,f12.6)') &
                            !    a_state, b_state, n_state, m_state,& 
                            !    e_diff, real(E_mp2_a), real(E_mp2_b), real(E_mp2_tmp)
                            !endif

                            pt2_term_os = pt2_term_os + real(E_mp2_tmp_os)
                            !pt2_term_ss = pt2_term_ss + real(E_mp2_tmp_ss)
                            !pt2_term = pt2_term + pt2_term_ss + pt2_term_os

                        enddo ! m_state for the virtual states of i_qp_point
                    enddo ! n_state for the virtual states of i_kp_point
                enddo ! b_state for the occupied states of i_q_point
            enddo ! a_state for the occupied states of i_k_point

            !write(use_unit,*) pt2_term_os
            pt2_c_energy_local_os = pt2_c_energy_local_os + pt2_term_os * weight_total
            !pt2_c_energy_local_ss = pt2_c_energy_local_ss + pt2_term_ss * weight_total
            !pt2_c_energy_local    = pt2_c_energy_local_os + pt2_c_energy_local_ss

        endif
        !call get_timestamps(time_debug(3),time_debug(4))
        !if (myid .eq. 1) then
        !    write(use_unit, '(2X,A35,2X,2F16.3)') "Timing for PT2 evaluation", &
        !        time_debug(3)-time_debug(1), &
        !        time_debug(4)-time_debug(2)
        !endif

        ! Reflash the restarting file, Igor
        ! pt2_c_energy_local :: the PT2 correlation in each thread
        ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
        ! pt2_c_energy       :: the collected PT2 correlation up to now
        if (mod(i_p_index, 10)-1 .eq. 0) then
            !pt2_c_energy    = pt2_c_energy_local
            !pt2_c_energy_ss = pt2_c_energy_local_ss
            pt2_c_energy_os = pt2_c_energy_local_os
            !call sync_real_number(pt2_c_energy)
            !call sync_real_number(pt2_c_energy_ss)
            call sync_real_number(pt2_c_energy_os)
            pt2_finish      = i_p_index
            current_pt2_os  = pt2_c_energy_os + pt2_c_energy_start_os
            current_pt2_ss  = 0.0d0
            current_pt2     = current_pt2_os 
            call write_restart_pt2_info()
            write(info_str, "(2X,'The progress of PT2 corr. : ',F8.3,'% (',F16.8,' Ha)')") &
                 & real(i_p_index)/real(cpt2_max_per_thread)*100d0, &
                 current_pt2
            call localorb_info(info_str)
        endif

        ! calculate the timing
        call get_timestamps(time_pt2(5),time_pt2(6))
        time_pt2(3) = time_pt2(5) - time_pt2(3)
        time_pt2(4) = time_pt2(6) - time_pt2(4)
        time_pt2(7) = time_pt2(3)
        time_pt2(8) = time_pt2(4)
        call sync_timing(time_pt2(3))
        call sync_timing(time_pt2(4))
        call sync_real_number(time_pt2(7))
        call sync_real_number(time_pt2(8))
        time_pt2(7) = time_pt2(7)/n_tasks
        time_pt2(8) = time_pt2(8)/n_tasks

        if (mod(i_p_index, 10)-1 .eq. 0) then
            write(info_str, "('PT2 corr. in the ',I5,' round :',F12.6)") &
                 & i_p_index
            call output_timeheader('2X', info_str)
            call output_timer('Total cost', time_pt2(3:4))
            call output_timer('Communication', time_distribution(3:4))
            write(info_str, "(2X,'Total time in average                        :',1X,F12.3,' s',4X,F12.3,' s')") &
                 & time_pt2(7), time_pt2(8)
            call localorb_info(info_str)
            call localorb_info('')
        endif

      enddo! i_p_index

      call mpi_barrier(mpi_comm_global,info)
      ! Note again:
      ! pt2_c_energy_local :: the PT2 correlation in each thread
      ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
      ! pt2_c_energy       :: the final total PT2 correlation
      !call sync_real_number(pt2_c_energy_local)
      call sync_real_number(pt2_c_energy_local_os)
      !call sync_real_number(pt2_c_energy_local_ss)
      !pt2_c_energy    = pt2_c_energy_start + pt2_c_energy_local
      pt2_c_energy_os = pt2_c_energy_start_os + pt2_c_energy_local_os
      pt2_c_energy_ss = 0.0d0
      pt2_c_energy    = pt2_c_energy_os

      call get_timestamps(time_pt2(5),time_pt2(6))
      time_pt2(3) = time_pt2(5) - time_pt2(1)
      time_pt2(4) = time_pt2(6) - time_pt2(2)
      call sync_timing(time_pt2(3))
      call sync_timing(time_pt2(4))
      call sync_timing(time_distribution(1))
      call sync_timing(time_distribution(2))

      call localorb_info('')
      write(info_str,*)"----------------------------------------------------", &
                  "-------------------------"
      call localorb_info(info_str)
      write(info_str,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
          " PT2 correlation (os)   :", pt2_c_energy_os, "Ha,", &
           pt2_c_energy_os*hartree, "eV"
      call localorb_info(info_str)
      call localorb_info('')
      call output_timeheader('2X', 'PT2 corr. calculation')
      call output_timer('Total cost', time_pt2(3:4))
      call output_timer('Communication', time_distribution(1:2))



      if (allocated (coulomb_tmp_qk)) then
        deallocate (coulomb_tmp_qk)
      endif
      if (allocated (coulomb_tmp_qpk)) then
        deallocate (coulomb_tmp_qpk)
      endif

      return

      end subroutine evaluate_cpt2_os_energy_kspace
