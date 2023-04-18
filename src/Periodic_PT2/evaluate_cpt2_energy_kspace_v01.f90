!****s* FHI-aims/evaluate_cpt2_energy_kspace
!  NAME
!   evaluate_cpt2_energy_kspace
!  SYNOPSIS

      subroutine evaluate_cpt2_energy_kspace &
           ( occ_numbers, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            pt2_c_energy, &
            pt2_c_energy_os, &
            pt2_c_energy_ss &
           )

!  PURPOSE
!  Subroutine evaluate_cpt2_energy_kspace evaluates the correlation
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

      real*8, dimension(:), allocatable :: pt2_c_kgrid
      real*8  weight_total
      

!    local timing
      real*8  temp_time_pt2
      real*8  temp_clock_time_pt2
      real*8  rtime_1, clock_rtime_1
      real*8  rtime_2, clock_rtime_2
      real*8  rtime_3, clock_rtime_3

      ! for parallelization
      integer :: i_index, j_index, n_index, m_index
      real*8  :: pt2_term, pt2_term_2

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied

      integer :: n_first(n_spin)
      integer :: ipiv(n_basbas)


      integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
      integer, allocatable :: n_lumo_k(:,:) ! # of non-fully occupied states at a given k-point q
      complex*16, allocatable :: lvl_tricoeff_sum(:,:) ! sum of LVL triple coefficient at k and q points
      complex*16, allocatable :: lvl_tricoeff_sum_all(:,:,:) ! lvl_tricoeff_sum for all babas values
      complex*16, allocatable :: coeff_tmp(:,:) ! temporary triple coefficient, one basis index is transformed to
                                                   ! (occupied) molecular orbital index 
      complex*16, allocatable :: coeff_tmp_2(:,:) ! temporary triple coefficient, two basis indices are transformed to molecular orbital indices.
      !complex*16, allocatable :: coeff_tmp_3(:,:) 
      !complex*16, allocatable :: coeff_tmp_4(:,:) 

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

      character(*), parameter :: func = 'evaluate_cpt2_energy_kspace.f90'
      integer :: info, mpierr

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

      integer :: i_task_1
      integer :: i_task_2
      integer :: i_task_3
      integer :: i_task_4
      integer :: i_task_5
      integer :: id_recv
      integer :: id_send

      integer :: n_unocc_max
      integer :: n_lumo(n_spin)
      integer :: n_lumo_min
      integer :: n_lumo_eff
      integer :: n_unocc_eff

!     begin work

      if (myid .eq. 1) then
          write(use_unit,*) "Igor debug k_weights :", k_weights
          write(use_unit,*) "Igor debug irk_weight :", irk_weight
      endif

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-----------------------------------------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the periodic PT2 correlation energy  ... "
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

! work array

      !allocate(pt2_c_kgrid(n_irk_points),stat=i_index)

!      do i_k_point = 1, n_k_points, 1
!         do i_prodbas_1 = 1, n_basbas, 1
!           write(use_unit,'(3I4,4f18.8)') i_k_point, i_prodbas_1, basbas_l(i_prodbas_1), &
!                 coulomb_matr_recip(i_prodbas_1,i_prodbas_1,i_k_point) , &
!!                multipole_basbas_fn(basbas_fn(i_prodbas_1))
!         enddo
!      enddo

      !time_pt2_corr = 0.d0
      !clock_time_pt2_corr = 0.d0
      !time_polar = 0.d0
      !clock_time_polar = 0.d0

      pt2_c_energy = 0.d0 
      !pt2_c_kgrid = 0.d0 

      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(n_low_state)
      else
          n_low_state = 1
      endif

      call get_timestamps(temp_time_pt2, temp_clock_time_pt2 )

! n_k_points_original is the number of k points in the original mesh
      ! i_qk_point : delta k = q - k
      do i_qk_point = 1, n_k_points, 1
         i_qk_point_local   = (i_qk_point-1)/n_tasks + 1
         i_irqk_point       = irk_point_mapping(i_qk_point)
         if(.not. irk_point_included(i_qk_point) ) cycle
         !i_irqk_point_local = (i_irqk_point-1)/n_tasks + 1
         ! broadcast the two-center coulomb matrix at i_qk_point
         if(myid.eq.mod(i_qk_point,n_tasks)) then
           coulomb_tmp_qk(:,:) = coulomb_matr_recip(:,:,i_qk_point_local)
         endif
         call mpi_bcast(coulomb_tmp_qk,n_basbas*n_basbas, &
                        MPI_COMPLEX16, mod(i_qk_point,n_tasks), mpi_comm_global, mpierr)

         do i_k_point = 1, n_k_points, 1
             i_k_point_local   = (i_k_point-1)/n_tasks + 1
             !if(.not. irk_point_included(i_k_point) ) cycle
             i_irk_point = irk_point_mapping(i_k_point)
             i_irk_point_local   = (i_irk_point-1)/n_tasks + 1
             ! Determine the q grid
             i_q_point = kpq_point_list(i_k_point, i_qk_point)
             i_q_point_local   = (i_q_point-1)/n_tasks + 1
             !if(.not. irk_point_included(i_q_point) ) cycle
             i_irq_point = irk_point_mapping(i_q_point)
             i_irq_point_local   = (i_irq_point-1)/n_tasks + 1

             ! loading the eigenvectors of k and q grids
             if(myid.eq.mod(i_k_point,n_tasks)) then
                  lvl_tricoeff_recip_k(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_k_point_local)
                  if(real_eigenvectors) then
                     KS_eigenvector_k(:,:,:)=KS_eigenvector(:,:,:,i_k_point_local)
                  else
                     KS_eigenvector_k(:,:,:)=KS_eigenvector_complex(:,:,:,i_k_point_local)
                  endif
             endif
             if(myid.eq.mod(i_q_point,n_tasks)) then
                  lvl_tricoeff_recip_q(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_q_point_local)
                  if(real_eigenvectors) then
                     KS_eigenvector_q(:,:,:)=KS_eigenvector(:,:,:,i_q_point_local)
                  else
                     KS_eigenvector_q(:,:,:)=KS_eigenvector_complex(:,:,:,i_q_point_local)
                  endif
             endif

             call mpi_bcast(lvl_tricoeff_recip_k,max_n_basbas_sp*n_basis*n_basis, &
                                     MPI_COMPLEX16, mod(i_k_point,n_tasks), mpi_comm_global, mpierr)
             call mpi_bcast(KS_eigenvector_k,n_basis*n_states*n_spin, &
                                     MPI_COMPLEX16, mod(i_k_point,n_tasks), mpi_comm_global, mpierr)
             call mpi_bcast(lvl_tricoeff_recip_q,max_n_basbas_sp*n_basis*n_basis, &
                                     MPI_COMPLEX16, mod(i_q_point,n_tasks), mpi_comm_global, mpierr)
             call mpi_bcast(KS_eigenvector_q,n_basis*n_states*n_spin, &
                                     MPI_COMPLEX16, mod(i_q_point,n_tasks), mpi_comm_global, mpierr)

             ! ==================
             ! for M(k,q)
             ! ==================
             lvl_tricoeff_sum_all     = (0.d0,0.d0)
             !i_index = 0
             do i_basis_1 = 1, n_basis, 1
               do i_basis_2 = 1, n_basis, 1
                 !i_index = i_index + 1
                 !if (myid .eq. MOD(i_index,n_tasks)) then
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
                 !endif
              enddo
            enddo
            !call sync_vector_complex(lvl_tricoeff_sum_all,size(lvl_tricoeff_sum_all))
            
            call get_timestamps(rtime_1, clock_rtime_1)
            lvl_tricoeff_kq(:,:,:,:) = (0.d0,0.d0)
            do i_prodbas_1 = 1, n_basbas, 1
              if (myid .eq. MOD(i_prodbas_1,n_tasks)) then
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
              endif ! mpi, myid .eq. MOD(i_index,n_tasks)
            ! end loop over i_prodbas_1
            enddo
            call get_timestamps(rtime_2, clock_rtime_2)
            call sync_vector_complex(lvl_tricoeff_kq, size(lvl_tricoeff_kq))
            call get_timestamps(rtime_3, clock_rtime_3)
            if (myid .eq. 0) then
                write(use_unit,'(F16.8,F16.8)') rtime_2 - rtime_1, clock_rtime_2 - clock_rtime_1 
                write(use_unit,'(F16.8,F16.8)') rtime_3 - rtime_2, clock_rtime_3 - clock_rtime_2 
            endif

            do i_qp_point = 1 , n_k_points, 1
                i_qp_point_local   = (i_qp_point-1)/n_tasks + 1
                !if(.not. irk_point_included(i_qp_point) ) cycle
                i_irqp_point = irk_point_mapping(i_qp_point)
                i_irqp_point_local   = (i_irqp_point-1)/n_tasks + 1
                ! Determine the k' grid
                i_kp_point = kpq_point_list(i_qp_point,i_qk_point)
                i_kp_point_local   = (i_kp_point-1)/n_tasks + 1
                !if(.not. irk_point_included(i_kp_point) ) cycle
                i_irkp_point = irk_point_mapping(i_kp_point)
                i_irkp_point_local   = (i_irkp_point-1)/n_tasks + 1
                ! Now we have the k-mesh weight as: w_k*w_q*w_k'*w_q'
                !weight_total = k_weights(i_qk_point)*k_weights(i_k_point) &
                !         * k_weights(i_qp_point)
                weight_total = irk_weight(i_irqk_point)*k_weights(i_k_point) &
                         * k_weights(i_qp_point)
                !weight_total = irk_weight(i_irqk_point)*irk_weight(i_irk_point) &
                !         * irk_weight(i_irqp_point)
                ! Determine the q'-k grid
                i_qpk_point = kq_point_list(i_qp_point,i_k_point)
                i_qpk_point_local = (i_qpk_point-1)/n_tasks + 1
                !if(.not. irk_point_included(i_qpk_point) ) cycle
                i_irqpk_point = irk_point_mapping(i_qpk_point)
                i_irqpk_point_local = (i_irqpk_point-1)/n_tasks + 1
                ! ==========================================
                ! For the parallel implementation
                ! ==========================================
                ! load the two-center coulomb matrix at i_qpk_point
                ! and send the two-center coulomb matrix at i_qpk_point to the
                ! right process
                if(myid.eq.mod(i_qpk_point,n_tasks)) then
                  coulomb_tmp_qpk(:,:) = coulomb_matr_recip(:,:,i_qpk_point_local)
                endif
                call mpi_bcast(coulomb_tmp_qpk,n_basbas*n_basbas, &
                               MPI_COMPLEX16, mod(i_qpk_point,n_tasks), mpi_comm_global, mpierr)
                ! communicate the q' grid information
                if(myid.eq.mod(i_qp_point, n_tasks)) then
                     lvl_tricoeff_recip_qp(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_qp_point_local)
                     if(real_eigenvectors) then
                        KS_eigenvector_qp(:,:,:)=KS_eigenvector(:,:,:,i_qp_point_local)
                     else
                        KS_eigenvector_qp(:,:,:)=KS_eigenvector_complex(:,:,:,i_qp_point_local)
                     endif
                endif
                call mpi_bcast(lvl_tricoeff_recip_qp,max_n_basbas_sp*n_basis*n_basis, &
                                        MPI_COMPLEX16, mod(i_qp_point,n_tasks), mpi_comm_global, mpierr)
                call mpi_bcast(KS_eigenvector_qp,n_basis*n_states*n_spin, &
                                        MPI_COMPLEX16, mod(i_qp_point,n_tasks), mpi_comm_global, mpierr)
                ! communicate the k' grid information
                if(myid.eq.mod(i_kp_point, n_tasks)) then
                     lvl_tricoeff_recip_kp(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_kp_point_local)
                     if(real_eigenvectors) then
                        KS_eigenvector_kp(:,:,:)=KS_eigenvector(:,:,:,i_kp_point_local)
                     else
                        KS_eigenvector_kp(:,:,:)=KS_eigenvector_complex(:,:,:,i_kp_point_local)
                     endif
                endif
                call mpi_bcast(lvl_tricoeff_recip_kp,max_n_basbas_sp*n_basis*n_basis, &
                                        MPI_COMPLEX16, mod(i_kp_point,n_tasks), mpi_comm_global, mpierr)
                call mpi_bcast(KS_eigenvector_kp,n_basis*n_states*n_spin, &
                                        MPI_COMPLEX16, mod(i_kp_point,n_tasks), mpi_comm_global, mpierr)
                ! ==================
                ! for M(kp,qp)
                ! ==================
                lvl_tricoeff_sum_all = (0.d0,0.d0)
                !i_index = 0
                do i_basis_1 = 1, n_basis, 1
                  do i_basis_2 = 1, n_basis, 1
                    !i_index = i_index+1
                    !if (myid .eq. MOD(i_index,n_tasks)) then
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
                    !endif
                  enddo
                enddo
                !call sync_vector_complex(lvl_tricoeff_sum_all,size(lvl_tricoeff_sum_all))
           
                lvl_tricoeff_kpqp(:,:,:,:) = (0.d0,0.d0)
                do i_prodbas_1 = 1, n_basbas, 1
                  if (myid .eq. MOD(i_prodbas_1, n_tasks)) then
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
                  endif ! for mpi, myid .eq. MOD(i_index, n_tasks)
                ! end loop over i_prodbas_1
                enddo
                call sync_vector_complex(lvl_tricoeff_kpqp,size(lvl_tricoeff_kpqp))

                ! ==================
                ! for M(k,qp)
                ! ==================
                lvl_tricoeff_sum_all = (0.d0,0.d0)
                !i_index = 0
                do i_basis_1 = 1, n_basis, 1
                  do i_basis_2 = 1, n_basis, 1
                    !i_index = i_index + 1
                    !if (myid .eq. MOD(i_index, n_tasks)) then
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
                    !endif
                  enddo
                enddo
                !call sync_vector_complex(lvl_tricoeff_sum_all, size(lvl_tricoeff_sum_all))
           
                lvl_tricoeff_kqp(:,:,:,:) = (0.d0,0.d0)
                do i_prodbas_1 = 1, n_basbas, 1
                  if (myid .eq. MOD(i_prodbas_1, n_tasks)) then
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
                  endif ! for mpi, myid.eq.MOD(i_index,n_tasks)
                ! end loop over i_prodbas_1
                enddo
                call sync_vector_complex(lvl_tricoeff_kqp,size(lvl_tricoeff_kqp))

                ! ==================
                ! for M(kp,q)
                ! ==================
                lvl_tricoeff_sum_all = (0.d0,0.d0)
                !i_index = 0
                do i_basis_1 = 1, n_basis, 1
                  do i_basis_2 = 1, n_basis, 1
                    !i_index = i_index+1
                    !if (myid.eq.MOD(i_index,n_tasks)) then
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
                    !endif
                  enddo
                enddo
                !call sync_vector_complex(lvl_tricoeff_sum_all,size(lvl_tricoeff_sum_all))
           
                lvl_tricoeff_kpq(:,:,:,:) = (0.d0,0.d0)
                do i_prodbas_1 = 1, n_basbas, 1
                  if (myid .eq. MOD(i_prodbas_1,n_tasks)) then
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
                  endif ! for mpi, myid.eq.MOD(i_index,n_tasks)
                ! end loop over i_prodbas_1
                enddo
                call sync_vector_complex(lvl_tricoeff_kpq,size(lvl_tricoeff_kpq))
                ! Now do PT2 calculation
                i_index = 0
                pt2_term = 0.0d0
                do a_state = n_low_state, n_homo_k(i_k_point, 1), 1
                    do b_state = n_low_state, n_homo_k(i_kp_point,1), 1
                        do n_state = n_lumo_k(i_q_point,1), n_states_k(i_q_point), 1
                            do m_state = n_lumo_k(i_qp_point,1), n_states_k(i_qp_point),1
                                i_index = i_index + 1
                                if (myid .eq. MOD(i_index,n_tasks)) then
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
                                    !if (a_state .eq. b_state .and. &
                                    !    i_k_point .eq. i_kp_point) then
                                    !    E_mp2_tmp = E_mp2_tmp * 2
                                    !endif

                                    pt2_term = pt2_term + real(E_mp2_tmp)

                                    !if (i_k_point .eq. 1 .and. &
                                    !    i_kp_point .eq. 1 .and. &
                                    !    i_q_point .eq. 1 .and. &
                                    !    i_qp_point .eq. 1) &
                                    !    write(use_unit,'(f10.3,f12.6,f12.6,f12.6)') &
                                    !     e_diff, real(E_mp2_a), real(E_mp2_b), real(E_mp2_tmp)


                                endif ! if (myid .eq. MOD(i_index,n_tasks) then
                                
                            enddo ! m_state for the virtual states of i_qp_point
                        enddo ! n_state for the virtual states of i_kp_point
                    enddo ! b_state for the occupied states of i_q_point
                enddo ! a_state for the occupied states of i_k_point


                if (use_mpi) then
                    pt2_term_2 = 0.0d0
                    call MPI_ALLREDUCE(pt2_term, pt2_term_2, 1, &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM, mpi_comm_global, mpierr)
                endif

                pt2_c_energy = pt2_c_energy + pt2_term_2 * weight_total

                if (myid .eq. 0) then
                    write(use_unit,*) &
                        "k, k', q, q' :", i_k_point, i_kp_point, i_q_point, i_qp_point
                    write(use_unit,'(A12,f12.6,A35,f12.6,A3)') &
                        'weight_total', weight_total, 'pt2_c_energy for this point : ', pt2_term_2*weight_total, "Ha"
                endif

            enddo ! for i_qp_point
         enddo ! for i_k_point
      enddo ! for i_qk_point

      !if (use_mpi) then
      !    call MPI_ALLREDUCE(pt2_c_energy, pt2_c_energy, 1, &
      !        MPI_DOUBLE_PRECISION, &
      !        MPI_SUM, mpi_comm_global, mpierr)
      !endif

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            " PT2 correlation energy :", pt2_c_energy, "Ha,", &
             pt2_c_energy*hartree, "eV"
!        write(use_unit,*)"----------------------------------------------------", &
!                 "-------------------------"
        write(use_unit,*)
      endif

      if (allocated (coulomb_tmp_qk)) then
        deallocate (coulomb_tmp_qk)
      endif
      if (allocated (coulomb_tmp_qpk)) then
        deallocate (coulomb_tmp_qpk)
      endif

      return

      end subroutine evaluate_cpt2_energy_kspace
