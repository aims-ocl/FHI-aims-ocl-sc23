module evaluate_cpt2_energy_skpace_mod

  !  USES
  use dimensions
  use prodbas, only : atom2basbas_off, sp2n_basbas_sp, max_n_basbas_sp, basbas_atom
  use hartree_fock, only : lvl_tricoeff_mod_r, n_homo_max, n_homo, kq_point_list
  use hartree_fock_p0
  use mpi_tasks
  use synchronize_mpi
  use pbc_lists, only : Cbasis_to_atom
  use runtime_choices
  use timing
  use constants
  use geometry, only : species
  use crpa_blacs
  use exchange_trico
  use exchange_ev
  implicit none

  save
  private
  public evaluate_cpt2_energy_kspace_v02

  !MODULE VARIABLES
  integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
  integer, allocatable :: n_lumo_k(:,:) ! # of non-fully occupied states at a given k-point q
  complex*16, dimension(:,:,:,:), allocatable, target :: lvl_tricoeff_k, lvl_tricoeff_q ! LVL triple coefficients in k space
  complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_k_r, lvl_tricoeff_recip_q_r ! LVL triple coefficients in k space

  character(*), parameter :: func = 'evaluate_polarisability_kspace.f90'

  integer, dimension(:), allocatable:: lb_atom, ub_atom

  logical:: use_adjoint=.true.
  logical:: doublereal

  contains
!****s* FHI-aims/evaluate_cpt2_energy_kspace_v02
!  NAME
!   evaluate_cpt2_energy_kspace_v02
!  SYNOPSIS

   subroutine evaluate_cpt2_energy_kspace_v02 &
           ( occ_numbers, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            pt2_c_energy, &
            pt2_c_energy_os, &
            pt2_c_energy_ss &
           )

!  PURPOSE
!  Subroutine evaluate_cpt2_energy_kspace_v02 evaluates the correlation
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

      use crpa_blacs
      use exchange_trico
      use exchange_ev

      implicit none

! ARGUMENTS 

      !integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state
!      integer :: n_lumo(n_spin)
!      integer :: n_homo(n_spin)

      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: e_diff
      !real*8  :: E_mp2_tmp
      complex*16  :: E_mp2_tmp
      !real*8  :: omega_full(n_full_freq)
      !real*8  :: womega_full(n_full_freq)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: E_mp2_a, E_mp2_b

!     output
      real*8  :: pt2_c_energy
      real*8  :: pt2_c_energy_os
      real*8  :: pt2_c_energy_ss

      ! mpi related variable
      integer:: win_tri, win_ev

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

      real*8  :: weight_total, weight_total_array(n_tasks)
      integer :: cpt2_para(6,n_tasks), i_tasks(6)
      integer :: i_p_index, i_task, i_task_real
      integer :: cpt2_total_grid
      integer :: cpt2_max_per_thread
      real*8  :: pt2_c_energy_local
      real*8  :: pt2_c_energy_start
      

!    local timing
      real*8  time_pt2(8)
      real*8  time_distribution(6)

      ! for parallelization
      integer :: i_index, j_index, n_index, m_index
      real*8  :: pt2_term

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
      complex*16, allocatable, target :: lvl_tricoeff_recip_q(:,:,:,:)  ! LVL triple coefficients in the q point
      complex*16, allocatable, target :: lvl_tricoeff_recip_k(:,:,:,:)  ! LVL triple coefficients in the k point 
      complex*16, allocatable, target :: lvl_tricoeff_recip_qp(:,:,:,:) ! LVL triple coefficients in the q' point
      complex*16, allocatable, target :: lvl_tricoeff_recip_kp(:,:,:,:) ! LVL triple coefficients in the k' point 
      complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_q_r  
      complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_k_r  
      complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_qp_r 
      complex*16, dimension(:,:,:,:), pointer :: lvl_tricoeff_recip_kp_r 
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

      character(*), parameter :: func = 'evaluate_cpt2_energy_kspace_v02.f90'
      integer :: info, mpierr
      character*150 :: info_str

!     timing

!     parameters of the fitting tails
!       real*8  s1, s2, omega_1, omega_2
!       real*8  alpha, beta, a, b

!     counters

      ! Three loops of k-mesh: 1) q-k (qk); 2) k; 3) q'(qp)
      integer :: i_qk_point
      integer :: i_qk_point_local
      integer :: i_irqk_point
      integer :: i_irqk_point_local
      integer :: i_k_point
      integer :: i_k_point_local
      integer :: i_irk_point
      integer :: i_irk_point_local
      integer :: i_qp_point
      integer :: i_qp_point_local
      integer :: i_irqp_point
      integer :: i_irqp_point_local
      ! We also need other three independent k-mesh arguments
      ! 1) q; 2) k'(kp); 3) q'-k (qpk)
      integer :: i_q_point
      integer :: i_q_point_local
      integer :: i_irq_point
      integer :: i_irq_point_local
      integer :: i_kp_point
      integer :: i_kp_point_local
      integer :: i_irkp_point
      integer :: i_irkp_point_local
      integer :: i_qpk_point
      integer :: i_qpk_point_local
      integer :: i_irqpk_point
      integer :: i_irqpk_point_local
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

!     begin work

      !if (myid .eq. 1) then
      !    write(use_unit,*) "Igor debug k_weights :", k_weights
      !    write(use_unit,*) "Igor debug irk_weight :", irk_weight
      !endif

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"  -----------------------------------------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the periodic PT2 correlation energy  ... "
        write(use_unit,'(2X,A)') &
              "Parallel in terms of BLACS  ... "
      endif

      !if(flag_KS_eigenfunc_conjg) then
      !   KS_eigenvector_complex = conjg(KS_eigenvector_complex)
      !endif
      if (real_eigenvectors) then
         call init_access_ev_real(n_k_points,n_k_points_task,KS_eigenvector,win_ev)
      else
         call init_access_ev_complex(n_k_points,n_k_points_task,KS_eigenvector_complex,win_ev)
      end if

      call sinit_access_trico(n_k_points,n_ks_points_task,lvl_tricoeff_mod_r,win_tri)
      
      !determine the k and q required
      !i_irkq_point=n_tasks_irkq*(i_irkq_point_local-1) + myid_irkq + 1

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
      allocate(lvl_tricoeff_recip_k(lbb_row:ubb_row,max_n_basis_sp,n_basis,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_k', func)
      lvl_tricoeff_recip_k_r => lvl_tricoeff_recip_k
      allocate(lvl_tricoeff_recip_q(lbb_row:ubb_row,max_n_basis_sp,n_basis,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_q', func)
      lvl_tricoeff_recip_q_r => lvl_tricoeff_recip_q
      allocate(lvl_tricoeff_recip_kp(lbb_row:ubb_row,max_n_basis_sp,n_basis,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_kp', func)
      lvl_tricoeff_recip_kp_r => lvl_tricoeff_recip_kp
      allocate(lvl_tricoeff_recip_qp(lbb_row:ubb_row,max_n_basis_sp,n_basis,n_spin),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_qp', func)
      lvl_tricoeff_recip_qp_r => lvl_tricoeff_recip_qp
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
      ! tricoeff for four k-points, including k, q, k', and q'
      allocate(lvl_tricoeff_recip_k_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_k_para', func)
      allocate(lvl_tricoeff_recip_q_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_q_para', func)
      allocate(lvl_tricoeff_recip_kp_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_kp_para', func)
      allocate(lvl_tricoeff_recip_qp_para(max_n_basbas_sp,n_basis,n_basis),stat=info) 
      call check_allocation(info, 'lvl_tricoeff_recip_qp_para', func)
      ! KS eigenvector matrices for four k-points, including k, q, k', and q'
      allocate(KS_eigenvector_k_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_k_para', func)
      allocate(KS_eigenvector_q_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_q_para', func)
      allocate(KS_eigenvector_kp_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_kp_para', func)
      allocate(KS_eigenvector_qp_para(n_basis,n_states,n_spin),stat=info) 
      call check_allocation(info, 'KS_eigenvector_qp_para', func)



      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(n_low_state)
      else
          n_low_state = 1
      endif

      call get_timestamps(time_pt2(1), time_pt2(2) )

      cpt2_total_grid = n_irk_points * n_k_points * n_k_points
      cpt2_max_per_thread = (cpt2_total_grid-1) / n_tasks + 1

      time_distribution       = 0.d0
      call get_timestamps(time_pt2(1), time_pt2(2))

      ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
      ! pt2_c_energy_local :: the PT2 correlation in each thread
      ! pt2_c_energy       :: the final PT2 correlation
      pt2_c_energy_start      = 0.d0 
      pt2_c_energy_local      = 0.d0 
      pt2_c_energy            = 0.d0 
      ! =====================================================
      ! reading the restart info., Igor
      ! the story behind the following five lines is:
      ! 1) in the 0th thread, loading the pt2_c_energy_start for the restart file, 
      !       if restart_pt2_read = .true.
      ! 2) broadcasting pt2_c_energy_start to all threads from the 0th thread
      ! 3) broadcasting pt2_finish to all threads from the 0th thread
      pt2_finish              = 0
      current_pt2             = 0.0d0
      call read_restart_pt2_info()
      pt2_c_energy_start      = current_pt2 ! only make sence for myid.eq.0 if restart_pt2_read = .true.
      call bcast_real(pt2_c_energy_start,0)
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
        ! Here, we consider the distribution from the thread with biggest index,
        !
        do i_task = n_tasks, 1, -1

            if(cpt2_para(1,i_task) .eq. 0) cycle

            !i_k_point   = cpt2_para(1,i_task)
            !i_kp_point  = cpt2_para(2,i_task)
            !i_q_point   = cpt2_para(3,i_task)
            !i_qp_point  = cpt2_para(4,i_task)
            !i_qk_point  = cpt2_para(5,i_task)
            !i_qpk_point = cpt2_para(6,i_task)

            i_tasks(1)  = mod(cpt2_para(1,i_task), n_tasks)
            i_tasks(2)  = mod(cpt2_para(2,i_task), n_tasks)
            i_tasks(3)  = mod(cpt2_para(3,i_task), n_tasks)
            i_tasks(4)  = mod(cpt2_para(4,i_task), n_tasks)
            i_tasks(5)  = mod(cpt2_para(5,i_task), n_tasks)
            i_tasks(6)  = mod(cpt2_para(6,i_task), n_tasks)

            i_k_point_local   = (cpt2_para(1,i_task)-1)/n_tasks + 1
            i_kp_point_local  = (cpt2_para(2,i_task)-1)/n_tasks + 1
            i_q_point_local   = (cpt2_para(3,i_task)-1)/n_tasks + 1
            i_qp_point_local  = (cpt2_para(4,i_task)-1)/n_tasks + 1
            i_qk_point_local  = (cpt2_para(5,i_task)-1)/n_tasks + 1
            i_qpk_point_local = (cpt2_para(6,i_task)-1)/n_tasks + 1

            write(use_unit,*) "igor debug 1"
            weight_total                 = weight_total_array(i_task)
            i_k_point                    = cpt2_para(1,i_task)
            call saccess_tricoeff(n_spin,win_tri,i_k_point, &
                 lvl_tricoeff_recip_k_r)
            call access_ev(win_ev,k_point_loc(:,i_k_point),KS_eigenvector_k)
            i_kp_point                   = cpt2_para(2,i_task)
            call saccess_tricoeff(n_spin,win_tri,i_kp_point, &
                 lvl_tricoeff_recip_kp_r)
            call access_ev(win_ev,k_point_loc(:,i_kp_point),KS_eigenvector_kp)
            i_q_point                    = cpt2_para(3,i_task)
            call saccess_tricoeff(n_spin,win_tri,i_q_point, &
                 lvl_tricoeff_recip_q_r)
            call access_ev(win_ev,k_point_loc(:,i_q_point),KS_eigenvector_q)
            i_qp_point                   = cpt2_para(4,i_task)
            call saccess_tricoeff(n_spin,win_tri,i_qp_point, &
                 lvl_tricoeff_recip_qp_r)
            call access_ev(win_ev,k_point_loc(:,i_qp_point),KS_eigenvector_qp)
            i_qk_point                   = cpt2_para(5,i_task)
            i_qpk_point                  = cpt2_para(6,i_task)
            coulomb_tmp_qk = coulomb_matr_blacs(:,:,i_qk_point)
            coulomb_tmp_qpk = coulomb_matr_blacs(:,:,i_qpk_point)
            write(use_unit,*) "igor debug 2"


            ! loading arrays 
            !if (myid.eq.i_tasks(1)) then ! for i_k_point
            !    lvl_tricoeff_recip_k_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_k_point_local)
            !    if(real_eigenvectors) then
            !       KS_eigenvector_k_para(:,:,:)=KS_eigenvector(:,:,:,i_k_point_local)
            !    else
            !       KS_eigenvector_k_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_k_point_local)
            !    endif
            !endif
            !write(use_unit,*) "igor debug 2"
            !if (myid.eq.i_tasks(2)) then ! for i_kp_point
            !    lvl_tricoeff_recip_kp_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_kp_point_local)
            !    if(real_eigenvectors) then
            !       KS_eigenvector_kp_para(:,:,:)=KS_eigenvector(:,:,:,i_kp_point_local)
            !    else
            !       KS_eigenvector_kp_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_kp_point_local)
            !    endif
            !endif
            !write(use_unit,*) "igor debug 3"
            !if (myid.eq.i_tasks(3)) then ! for i_q_point
            !    lvl_tricoeff_recip_q_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_q_point_local)
            !    if(real_eigenvectors) then
            !       KS_eigenvector_q_para(:,:,:)=KS_eigenvector(:,:,:,i_q_point_local)
            !    else
            !       KS_eigenvector_q_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_q_point_local)
            !    endif
            !endif
            !write(use_unit,*) "igor debug 4"
            !if (myid.eq.i_tasks(4)) then ! for i_qp_point
            !    lvl_tricoeff_recip_qp_para(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_qp_point_local)
            !    if(real_eigenvectors) then
            !       KS_eigenvector_qp_para(:,:,:)=KS_eigenvector(:,:,:,i_qp_point_local)
            !    else
            !       KS_eigenvector_qp_para(:,:,:)=KS_eigenvector_complex(:,:,:,i_qp_point_local)
            !    endif
            !endif
            !write(use_unit,*) "igor debug 5"
            !if (myid.eq.i_tasks(5)) then ! for i_qk_point
            !    coulomb_tmp_qk_para(:,:) = coulomb_matr_recip(:,:,i_qk_point_local)
            !endif
            !if (myid.eq.i_tasks(6)) then ! for i_qpk_point
            !    coulomb_tmp_qpk_para(:,:) = coulomb_matr_recip(:,:,i_qpk_point_local)
            !endif

            !! distribute relevant arrays to the target "i_task"
            !if (i_task .eq. n_tasks) then
            !    i_task_real = 0
            !else
            !    i_task_real = i_task
            !endif

            !if (myid.ne.i_task_real) then
            !  if (myid .eq. i_tasks(1)) then ! for i_k_point
            !    call send_complex_vector(lvl_tricoeff_recip_k_para,&
            !                                size(lvl_tricoeff_recip_k_para),i_task_real)
            !    call send_complex_vector(KS_eigenvector_k_para,&
            !                                size(KS_eigenvector_k_para),i_task_real)
            !  endif
            !  if (myid .eq. i_tasks(2)) then ! for i_kp_point
            !    call send_complex_vector(lvl_tricoeff_recip_kp_para,&
            !                                size(lvl_tricoeff_recip_kp_para),i_task_real)
            !    call send_complex_vector(KS_eigenvector_kp_para,&
            !                                size(KS_eigenvector_kp_para),i_task_real)
            !  endif
            !  if (myid .eq. i_tasks(3)) then ! for i_q_point
            !    call send_complex_vector(lvl_tricoeff_recip_q_para,&
            !                                size(lvl_tricoeff_recip_q_para),i_task_real)
            !    call send_complex_vector(KS_eigenvector_q_para,&
            !                                size(KS_eigenvector_q_para),i_task_real)
            !  endif
            !  if (myid .eq. i_tasks(4)) then ! for i_qp_point
            !    call send_complex_vector(lvl_tricoeff_recip_qp_para,&
            !                                size(lvl_tricoeff_recip_qp_para),i_task_real)
            !    call send_complex_vector(KS_eigenvector_qp_para,&
            !                                size(KS_eigenvector_qp_para),i_task_real)
            !  endif
            !  if (myid .eq. i_tasks(5)) then ! for i_qk_point
            !    call send_complex_vector(coulomb_tmp_qk_para,&
            !                                size(coulomb_tmp_qk_para),i_task_real)
            !  endif
            !  if (myid .eq. i_tasks(6)) then ! for i_qpk_point
            !    call send_complex_vector(coulomb_tmp_qpk_para,&
            !                                size(coulomb_tmp_qpk_para),i_task_real)
            !  endif
            !else
            !  if (myid .ne. i_tasks(1)) then ! for i_k_point
            !    call receive_complex_vector(lvl_tricoeff_recip_k_para,&
            !                                size(lvl_tricoeff_recip_k_para),i_tasks(1))
            !    call receive_complex_vector(KS_eigenvector_k_para,&
            !                                size(KS_eigenvector_k_para),i_tasks(1))
            !  endif
            !  if (myid .ne. i_tasks(2)) then ! for i_kp_point
            !    call receive_complex_vector(lvl_tricoeff_recip_kp_para,&
            !                                size(lvl_tricoeff_recip_kp_para),i_tasks(2))
            !    call receive_complex_vector(KS_eigenvector_kp_para,&
            !                                size(KS_eigenvector_kp_para),i_tasks(2))
            !  endif
            !  if (myid .ne. i_tasks(3)) then ! for i_q_point
            !    call receive_complex_vector(lvl_tricoeff_recip_q_para,&
            !                                size(lvl_tricoeff_recip_q_para),i_tasks(3))
            !    call receive_complex_vector(KS_eigenvector_q_para,&
            !                                size(KS_eigenvector_q_para),i_tasks(3))
            !  endif
            !  if (myid .ne. i_tasks(4)) then ! for i_qp_point
            !    call receive_complex_vector(lvl_tricoeff_recip_qp_para,&
            !                                size(lvl_tricoeff_recip_qp_para),i_tasks(4))
            !    call receive_complex_vector(KS_eigenvector_qp_para,&
            !                                size(KS_eigenvector_qp_para),i_tasks(4))
            !  endif
            !  if (myid .ne. i_tasks(5)) then ! for i_qk_point
            !    call receive_complex_vector(coulomb_tmp_qk_para,&
            !                                size(coulomb_tmp_qk_para),i_tasks(5))
            !  endif
            !  if (myid .ne. i_tasks(6)) then ! for i_qpk_point
            !    call receive_complex_vector(coulomb_tmp_qpk_para,&
            !                                size(coulomb_tmp_qpk_para),i_tasks(6))
            !  endif

            !  ! Finalize the memory distribution
            !  lvl_tricoeff_recip_k(:,:,:)  = lvl_tricoeff_recip_k_para(:,:,:)
            !  lvl_tricoeff_recip_kp(:,:,:) = lvl_tricoeff_recip_kp_para(:,:,:)
            !  lvl_tricoeff_recip_q(:,:,:)  = lvl_tricoeff_recip_q_para(:,:,:)
            !  lvl_tricoeff_recip_qp(:,:,:) = lvl_tricoeff_recip_qp_para(:,:,:)
            !  KS_eigenvector_k(:,:,:)      = KS_eigenvector_k_para(:,:,:)
            !  KS_eigenvector_kp(:,:,:)     = KS_eigenvector_kp_para(:,:,:)
            !  KS_eigenvector_q(:,:,:)      = KS_eigenvector_q_para(:,:,:)
            !  KS_eigenvector_qp(:,:,:)     = KS_eigenvector_qp_para(:,:,:)
            !  coulomb_tmp_qk(:,:)          = coulomb_tmp_qk_para(:,:)
            !  coulomb_tmp_qpk(:,:)         = coulomb_tmp_qpk_para(:,:)
            !  weight_total                 = weight_total_array(i_task)
            !  i_k_point                    = cpt2_para(1,i_task)
            !  i_kp_point                   = cpt2_para(2,i_task)
            !  i_q_point                    = cpt2_para(3,i_task)
            !  i_qp_point                   = cpt2_para(4,i_task)
            !  i_qk_point                   = cpt2_para(5,i_task)
            !  i_qpk_point                  = cpt2_para(6,i_task)
            !endif

            !!call mpi_barrier(mpi_comm_global,info)

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

        !write(use_unit,*) "igor debug", myid, cpt2_para(1,i_task_real),&
        !cpt2_para(2,i_task_real),cpt2_para(3,i_task_real),cpt2_para(4,i_task_real)

        if (cpt2_para(1,i_task_real) .ne. 0) then
            ! ==================
            ! for M(k,q)
            ! ==================
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

            ! ==================
            ! for M(k,qp)
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
                  lvl_tricoeff_recip_k(1:n_spbb_2, i_basis_2, i_basis_1) + &
                  lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)
              enddo
            enddo
            
            lvl_tricoeff_kqp(:,:,:,:) = (0.d0,0.d0)
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
                           KS_eigenvector_qp(1,1,i_spin), n_basis,  &
                           coeff_tmp, n_basis, (0.d0,0.d0), &
                           coeff_tmp_2, n_states)
    
                lvl_tricoeff_kqp(i_prodbas_1,:,:,i_spin) = coeff_tmp_2(:,:)
              ! end loop over i_spin
              enddo
            ! end loop over i_prodbas_1
            enddo

            ! ==================
            ! for M(kp,q)
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
                  conjg(lvl_tricoeff_recip_q(1:n_spbb_1, i_basis_1, i_basis_2)) + &
                  lvl_tricoeff_recip2(1:n_spbb_1, i_basis_2, i_basis_1)
    
                lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) = &
                lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,i_basis_2,i_basis_1) + &
                  lvl_tricoeff_recip_kp(1:n_spbb_2, i_basis_2, i_basis_1) + &
                  lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)
              enddo
            enddo
            
            lvl_tricoeff_kpq(:,:,:,:) = (0.d0,0.d0)
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
                           KS_eigenvector_q(1,1,i_spin), n_basis,  &
                           coeff_tmp, n_basis, (0.d0,0.d0), &
                           coeff_tmp_2, n_states)
    
                lvl_tricoeff_kpq(i_prodbas_1,:,:,i_spin) = coeff_tmp_2(:,:)
              ! end loop over i_spin
              enddo
            ! end loop over i_prodbas_1
            enddo

            ! Now do PT2 calculation
            pt2_term = 0.0d0
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
                            E_mp2_b = 0.0d0
                            do i_prodbas_1 = 1, n_basbas, 1
                                do i_prodbas_2 = 1, n_basbas, 1
                                    E_mp2_a = E_mp2_a + &
                                      coulomb_tmp_qk(i_prodbas_2,i_prodbas_1) * &
                                      lvl_tricoeff_kq(i_prodbas_2,n_state,a_state,1) * &
                                      lvl_tricoeff_kpqp(i_prodbas_1,m_state,b_state,1)
                                    E_mp2_b = E_mp2_b + &
                                      coulomb_tmp_qpk(i_prodbas_2,i_prodbas_1) * &
                                      lvl_tricoeff_kqp(i_prodbas_2,m_state,a_state,1) * &
                                      lvl_tricoeff_kpq(i_prodbas_1,n_state,b_state,1)
                                enddo ! i_prodbas_2
                            enddo  ! i_prodbas_1
                            E_mp2_tmp = conjg(E_mp2_a)*(2.0d0*E_mp2_a - E_mp2_b)
                            E_mp2_tmp = E_mp2_tmp / e_diff

                            !if (i_p_index .eq. 1 .and. myid .eq. 1) then
                            !    write(use_unit,'(f10.3,f12.6,f12.6,f12.6)') &
                            !    e_diff, real(E_mp2_a), real(E_mp2_b), real(E_mp2_tmp)
                            !endif

                            pt2_term = pt2_term + real(E_mp2_tmp)

                        enddo ! m_state for the virtual states of i_qp_point
                    enddo ! n_state for the virtual states of i_kp_point
                enddo ! b_state for the occupied states of i_q_point
            enddo ! a_state for the occupied states of i_k_point

            pt2_c_energy_local = pt2_c_energy_local + pt2_term * weight_total

        endif

        ! Reflash the restarting file, Igor
        ! pt2_c_energy_local :: the PT2 correlation in each thread
        ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
        ! pt2_c_energy       :: the collected PT2 correlation up to now
        !if (mod(i_p_index, 10) .eq. 0) then
            pt2_c_energy = pt2_c_energy_local
            call sync_real_number(pt2_c_energy)
            pt2_finish   = i_p_index
            current_pt2  = pt2_c_energy + pt2_c_energy_start
            call write_restart_pt2_info()
        !endif

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
        write(info_str, "('PT2 corr. in the ',I5,' round',F12.6,' Ha')") &
             & i_p_index, pt2_c_energy
        call output_timeheader('2X', info_str)
        call output_timer('Total cost', time_pt2(3:4))
        call output_timer('Communication', time_distribution(3:4))
        write(info_str, "(2X,'Total time in average :',F12.3,F12.3)") &
             & time_pt2(7), time_pt2(8)
        call localorb_info(info_str)
        write(info_str, "(2X,'The progress of PT2 corr. : ',F8.3,'%')") &
             & real(i_p_index)/real(cpt2_max_per_thread)*100d0
        call localorb_info(info_str)
        call localorb_info('')

      enddo! i_p_index

      call mpi_barrier(mpi_comm_global,info)
      ! Note again:
      ! pt2_c_energy_local :: the PT2 correlation in each thread
      ! pt2_c_energy_start :: unfinished value inheriting from the previous calculation
      ! pt2_c_energy       :: the final total PT2 correlation
      call sync_real_number(pt2_c_energy_local)
      pt2_c_energy = pt2_c_energy_start + pt2_c_energy_local

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
          " PT2 correlation energy :", pt2_c_energy, "Ha,", &
           pt2_c_energy*hartree, "eV"
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
      end subroutine evaluate_cpt2_energy_kspace_v02

  subroutine compute_lvl_kq_complex_cpt2(i_spin, &
       lvl_tricoeff_kq, KS_eigenvector_q, KS_eigenvector_k, lvl_tricoeff_recip_k, &
       lvl_tricoeff_recip_q, lbb,ubb, state_pairs, n_pairs)
  
    implicit none
    integer :: lbb, ubb, i_state, lb1, ub1
    integer , intent(in) :: i_spin
    complex*16, dimension(lbb:ubb,max_n_basis_sp,n_states,n_spin), intent(in) :: lvl_tricoeff_recip_k, lvl_tricoeff_recip_q  
    complex*16, intent(inout) :: lvl_tricoeff_kq(lbb:ubb,n_pairs) ! Full LVL triple coefficient, one basis index is transformed to
  
    complex*16, dimension(n_basis,n_states,n_spin), intent(in) :: KS_eigenvector_q, KS_eigenvector_k
  
    integer:: i_atom, i_species, lb, ub, i_state_2, brange, lbb_atom, ubb_atom, ind, n_pairs, i_pair
    integer, dimension(2,n_pairs):: state_pairs
  
    call perfon('polind')
  
  !  do i_state_2=lb1,ub1
    do i_pair=1, n_pairs
       i_state=state_pairs(1,i_pair)
       i_state_2=state_pairs(2,i_pair)
       do i_atom=basbas_atom(lbb),basbas_atom(ubb)        
          !range in basis dimension
          lb=lb_atom(i_atom)
          ub=ub_atom(i_atom)
          brange=ub_atom(i_atom)-lb_atom(i_atom)+1
          
          !range in basbas dimension
          i_species = species(i_atom)
          lbb_atom=max(lbb,atom2basbas_off(i_atom)+1)
          ubb_atom=min(ubb,atom2basbas_off(i_atom)+sp2n_basbas_sp(i_species))
  
          lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)=0.
          do ind=1,brange
             lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)=lvl_tricoeff_kq(lbb_atom:ubb_atom,i_pair)+&
                  conjg(KS_eigenvector_q(lb+ind-1,i_state_2,i_spin)) * lvl_tricoeff_recip_k(lbb_atom:ubb_atom,ind,i_state,i_spin) +&
                  KS_eigenvector_k(lb+ind-1,i_state,i_spin) * conjg(lvl_tricoeff_recip_q(lbb_atom:ubb_atom,ind,i_state_2,i_spin))
          end do
  
       end do
    end do
  
    call perfoff
  
  end subroutine compute_lvl_kq_complex_cpt2

end module evaluate_cpt2_energy_skpace_mod
