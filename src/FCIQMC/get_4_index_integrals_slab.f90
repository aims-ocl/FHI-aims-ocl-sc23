
!  NAME get_4_index_integrals_slab
!   
!  SYNOPSIS

subroutine get_4_index_integrals_slab ()
!( omega_n, KS_eigenvalue, KS_eigenvector, &
!  KS_eigenvector_complex, &
!  occ_numbers, polar_kspace)

  !  PURPOSE
  !  Subroutine get_4_index_integrals_slab.f90 evaluates the exact-exchange part of
  !  Hartree-Fock hamiltonian in a periodic system. The algorithm used here are based
  !  on the reciprocal space and localized resolution of identity (RI-LVL)
  !
  !  USES

  use dimensions
  use prodbas
  use hartree_fock
  use hartree_fock_p0
  use mpi_tasks
  use synchronize_mpi
  use pbc_lists
  use runtime_choices
  use timing
  use constants
  use fciqmc_module
  implicit none

  !  ARGUMENTS

  !real*8 :: omega_n
  !real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
  !real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
  !complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
  !real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
  !complex*16, dimension(n_basbas,n_basbas,n_irk_points_task) :: polar_kspace

  !  INPUTS
  !  o  omega_n -- one frequency point (in the imaginary axis)
  !  o  occ_numbers -- real array,
  !       the occupation number of the electrons for each eigenstate and each spin
  !  o  KS_eigenvalue -- real array,
  !            the eigenvalue of the single-particle (KS/HF) self-consistent calculation
  !  o  KS_eigenvector -- real array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  o  KS_eigenvector_complex -- complex array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  OUTPUTS
  !  o  polar_kspace -- the polarisability matrix in k-space
  !  the exact exchange matrix (the "hf_exchange_matr_complex" in the source code) is defined in MODULE
  !  hartree_fock_p0
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

! FIXME: need to change the data structure of lvl_tricoeff_kq_2
! FIXME: need to distribute the memory storage of polar_kspace

  integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
  integer, allocatable :: n_lumo_k(:,:) ! # of non-fully occupied states at a given k-point q
  complex*16, allocatable :: lvl_tricoeff_recip_q(:,:,:) ! LVL triple coefficients in k space
  complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector at a given q point in the BZ
  complex*16, allocatable :: KS_eigenvector_k(:,:,:) ! KS eigenvector at a given k point in the BZ
  complex*16, allocatable :: lvl_tricoeff_sum(:,:) ! sum of LVL triple coefficient at k and q points
  complex*16, allocatable :: lvl_tricoeff_sum_all(:,:,:) ! lvl_tricoeff_sum for all babas values
  complex*16, allocatable :: coeff_tmp(:,:) ! temporary triple coefficient, one basis index is transformed to
                                               ! (occupied) molecular orbital index 
  complex*16, allocatable :: coeff_tmp_2(:,:) ! temporary triple coefficient, two basis indices are transformed to molecular orbital indices.
  complex*16, allocatable :: lvl_tricoeff_kq(:,:,:,:) ! Full LVL triple coefficient, one basis index is transformed to
                                                    ! (occupied) molecular orbital index 
  complex*16, dimension(:,:), allocatable :: coulomb_tmp
  complex*16, allocatable :: lvl_tricoeff_tmp(:,:) ! temporary triple coefficients multiplied by a factor
  complex*16, allocatable :: lvl_tricoeff_tmp_2(:,:) ! temporary triple coefficients
!  complex*16, allocatable :: exchange_matr_tmp(:,:) ! temparary exchange matrix per k_point
!  real*8, allocatable :: coulomb_matr_local(:,:) ! temparary exchange matrix per k_point
  real*8 :: e_diff
  real*8 :: occ_diff
  real*8 :: e_diff1, e_diff2
  real*8 :: occ_diff1, occ_diff2

  integer :: info, mpierr
  character(*), parameter :: func = 'get_4_index_integrals_slab.f90'
  character*150 :: info_str

!  timing info
   real*8   ::  time_polar_1, clock_time_polar_1
   real*8   ::  tot_time_polar_1, tot_clock_time_polar_1
   real*8   ::  time_polar_2, clock_time_polar_2
   real*8   ::  tot_time_polar_2, tot_clock_time_polar_2
   real*8   ::  time_polar_3, clock_time_polar_3
   real*8   ::  tot_time_polar_3, tot_clock_time_polar_3

!  integer :: n_homo_max
  integer :: n_unocc_max
  integer :: n_lumo(n_spin)
  integer :: n_lumo_min
  integer :: n_lumo_eff
  integer :: n_unocc_eff
!  real*8 :: k_minus_q(3), k_minus_q_lattvec(3)
! counter 
  integer i_cell
  integer i_cell_local
  integer i_k_point
  integer i_q_point
  integer i_kq_point
  integer i_irkq_point
  integer i_k_point_local
  integer i_irk_point
  integer i_irk_point_local
  integer i_q_point_local
  integer i_kq_point_local
  integer i_irkq_point_local
  integer i_basis_1
  integer i_basis_2
  integer i_prodbas_1
  integer i_prodbas_2
  integer i_state
  integer i_state_1
  integer i_state_2
  integer i_state_3
  integer i_state_4
  integer i_spin
  integer id_root
  integer i_task
  integer id_send
  integer i_send_kq
  integer id_recv
  integer i_atom_1, i_species_1, bboff_1, n_spbb_1
  integer i_atom_2, i_species_2, bboff_2, n_spbb_2

  character(len=128) :: file_name_1,file_name_2
  integer            :: unit_fci_1, unit_fci_2
  complex*16         :: integral_4ks_tmp
  real*8             :: exx_energy

  unit_fci_1  = 666 + myid
  unit_fci_2  = 6666 + myid

  allocate(n_homo_k(n_k_points,n_spin),stat=info) 
  call check_allocation(info, 'n_homo_k', func)
  allocate(n_lumo_k(n_k_points,n_spin),stat=info) 
  call check_allocation(info, 'n_unocc_k', func)

  n_homo_k(:,:) = 0
  n_lumo_k(:,:)= n_states
  do i_spin = 1, n_spin, 1
    do i_k_point = 1, n_k_points, 1
      do i_state = 1, n_states
        if(occ_numbers(i_state,i_spin,i_k_point) .gt. 1.e-12) then
         n_homo_k(i_k_point,i_spin) = i_state
        endif
      enddo
     
      do i_state = n_states, 1, -1
        if(occ_numbers(i_state,i_spin,i_k_point) .lt. 1.d0) then
         n_lumo_k(i_k_point,i_spin) = i_state 
        endif
      enddo
!      if(myid.eq.0) then
!        write(use_unit,*) i_spin, i_k_point, n_homo_k(i_k_point, i_spin), n_lumo_k(i_k_point,i_spin)
!      endif
    enddo
  enddo

  if (myid .eq. 0) then
      write(use_unit,*) "igor debug homo              : ",n_homo_k(:,1)
      write(use_unit,*) "igor debug k_weights         : ",k_weights(:)
      write(use_unit,*) "igor debug irk_weight        : ",irk_weight(:)
      write(use_unit,*) "igor debug irk_point_mapping : ", irk_point_mapping(:)
      write(use_unit,*) "igor debug irk_point_included: ", irk_point_included(:)
      do i_k_point = 1, n_k_points, 1
      !    write(use_unit,*) "igor debug occ_numbers in ",i_k_point," :", &
      !        occ_numbers(:,1,i_k_point)
          write(use_unit,*) "igor debug kq_point_list     :", &
              kq_point_list(:,i_k_point)
      enddo
  endif

  n_homo(1)=maxval(n_homo_k(:,1), n_k_points)
  n_homo(n_spin)=maxval(n_homo_k(:,n_spin), n_k_points)
  n_homo_max = max(n_homo(1), n_homo(n_spin)) 

  n_lumo(1)=minval(n_lumo_k(:,1), n_k_points)
  n_lumo(n_spin)=minval(n_lumo_k(:,n_spin), n_k_points)
  n_lumo_min = min(n_lumo(1), n_lumo(n_spin)) 
  n_unocc_max = n_states - n_lumo_min + 1

  allocate(coulomb_tmp(n_basbas,n_basbas),stat=info) 
  call check_allocation(info, 'coulomb_tmp', func)
  allocate(lvl_tricoeff_recip_q(max_n_basbas_sp,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_recip_q', func)
  allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
  call check_allocation(info, 'KS_eigenvector_q', func)
  allocate(KS_eigenvector_k(n_basis,n_states,n_spin),stat=info) 
  call check_allocation(info, 'KS_eigenvector_k', func)
  allocate(lvl_tricoeff_sum(n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_sum', func)
  allocate(lvl_tricoeff_sum_all(n_basbas,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_sum_sll', func)
  allocate(coeff_tmp(n_basis,n_basis),stat=info) 
  call check_allocation(info, 'coeff_tmp', func)
  allocate(coeff_tmp_2(n_states,n_basis),stat=info) 
  call check_allocation(info, 'coeff_tmp_2', func)
  allocate(lvl_tricoeff_kq(n_basbas,n_states,n_states,n_spin),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_kq', func)
  allocate(lvl_tricoeff_tmp(n_basbas,n_states),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_tmp', func)
  allocate(lvl_tricoeff_tmp_2(n_basbas,n_states),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_tmp_2', func)


!  write(info_str,'(2X,A)') "Evaluating the polarisability matrix in k space ..."
!  call localorb_info(info_str)

!  tot_time_polar_1 = 0.d0
!  tot_clock_time_polar_1 = 0.d0
!  tot_time_polar_2 = 0.d0
!  tot_clock_time_polar_2 = 0.d0
!  tot_time_polar_3 = 0.d0
!  tot_clock_time_polar_3 = 0.d0
  
  exx_energy = 0.0d0

  do i_q_point = 1, n_k_points, 1
    id_root = mod(i_q_point,n_tasks)
    i_q_point_local = (i_q_point-1)/n_tasks + 1
    if(myid .eq. id_root) then
        lvl_tricoeff_recip_q(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_q_point_local)
        if(real_eigenvectors) then
            KS_eigenvector_q(:,:,:)=KS_eigenvector(:,:,:,i_q_point_local)
        else
            KS_eigenvector_q(:,:,:)=KS_eigenvector_complex(:,:,:,i_q_point_local)
        endif
    endif

    call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!

    call mpi_bcast(lvl_tricoeff_recip_q,max_n_basbas_sp*n_basis*n_basis, &
                   MPI_COMPLEX16, id_root, mpi_comm_global, mpierr)

    call mpi_bcast(KS_eigenvector_q,n_basis*n_states*n_spin, &
                   MPI_COMPLEX16, id_root, mpi_comm_global, mpierr)

    !KS_eigenvector_q = conjg(KS_eigenvector_q)
   
    do i_k_point = 1, n_k_points, 1
      !call get_timestamps(time_polar_1, clock_time_polar_1)
      !write(use_unit,*) "i_k_point", k_point_list(i_k_point,1:3)
      i_irk_point = irk_point_mapping(i_k_point)
      i_task = mod(i_k_point,n_tasks)
      i_k_point_local = (i_k_point-1)/n_tasks + 1 

      i_kq_point = kq_point_list(i_k_point,i_q_point)
      i_kq_point_local = (i_kq_point-1)/n_tasks + 1
      ! only handle the irreducible k points
      if(.not. irk_point_included(i_kq_point) ) cycle

      if(myid.eq.i_task) then
          i_irkq_point = irk_point_mapping(i_kq_point)
          i_irkq_point_local = (i_irkq_point-1)/n_tasks + 1
          ! distribute the coulomb integrations of aux basis sets
          i_send_kq = mod(i_irkq_point, n_tasks)
          !id_recv = mod(i_irk_point, n_tasks)
          if(myid.eq.i_send_kq) then
              coulomb_tmp(:,:) = coulomb_matr_recip(:,:,i_irkq_point_local)
          endif
          if(i_send_kq.ne.i_task) then
              if(myid.eq.i_send_kq) then
                  call send_complex_vector(coulomb_tmp,n_basbas*n_basbas,i_task)
              elseif (myid.eq.i_task) then
                  call receive_complex_vector(coulomb_tmp,n_basbas*n_basbas,i_send_kq)
              endif
          endif

          call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!

          if(real_eigenvectors) then
              KS_eigenvector_k(:,:,:)=KS_eigenvector(:,:,:,i_k_point_local)
          else
              KS_eigenvector_k(:,:,:)=KS_eigenvector_complex(:,:,:,i_k_point_local)
          endif
          !write(use_unit,*) i_k_point, i_q_point, i_kq_point

          ! Set lvl_tricoeff_sum_all
          ! See evaluate_exchange_matr_kspace_p0.f90 for a more efficient way doing this!

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
                lvl_tricoeff_recip1(1:n_spbb_2, i_basis_2, i_basis_1, i_k_point_local) + &
                lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)

            enddo
          enddo
       
          lvl_tricoeff_kq(:,:,:,:) = (0.d0,0.d0)

          !call get_times(time_polar_1, clock_time_polar_1)
          !tot_time_polar_1 = tot_time_polar_1 + time_polar_1
          !tot_clock_time_polar_1 = tot_clock_time_polar_1 + clock_time_polar_1

          !call get_timestamps(time_polar_2, clock_time_polar_2)
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
              call zgemm('N', 'N', n_basis, n_basis, n_basis, (1.d0,0.d0), &
                         lvl_tricoeff_sum, n_basis, &
                         KS_eigenvector_k(1,1,i_spin), n_basis, (0.d0,0.d0), &
                         coeff_tmp, n_basis)

              coeff_tmp_2(:,:) = (0.d0,0.d0)
              call zgemm('C', 'N', n_states, n_basis, n_basis, (1.d0,0.d0), &
                         KS_eigenvector_q(1,1,i_spin), n_basis,  &
                         coeff_tmp, n_basis, (0.d0,0.d0), &
                         coeff_tmp_2, n_states)

              lvl_tricoeff_kq(i_prodbas_1,:,:,i_spin) = coeff_tmp_2(:,:)
              !if(i_prodbas_1.eq.1 ) then
              !   write(use_unit,'(50f18.8)') lvl_tricoeff_kq(i_prodbas_1,:,1,i_spin)
              !endif
            ! end loop over i_spin
            enddo
          ! end loop over i_prodbas_1
          enddo

          !call get_times(time_polar_2, clock_time_polar_2)
          !tot_time_polar_2 = tot_time_polar_2 + time_polar_2
          !tot_clock_time_polar_2 = tot_clock_time_polar_2 + clock_time_polar_2

      !endif (myid.eq.i_task) 
      endif

      ! mathematics behind: (\psi_{ik}\psi_{jp}|\psi_{mk}\psi_{np}) =
      !                        \sum_{\mu,\nv} C_{i,j}^{\mu}(k,q) * C_{m,n}^{\nv*}(k,q)
      !                        * (P_{\mu}^{k-q}|P_{\nv}^{q-k})

100 format('FCIDMP-4-index-',I1,'-',I1,'-',I1,'-',I1)
101 format('FCIDMP-4-index-',I1,'-',I1,'-',I1,'-',I1,'-real')
102 format('FCIDMP-4-index-',I1,'-',I1,'-',I1,'-',I1,'-imag')
103 FORMAT(1X,E23.16,4I4)
104 FORMAT(1X,4I4)
105 FORMAT(1X,G23.16,4I4)

      !write(file_name_1,100) i_q_point, i_k_point, i_k_point, i_q_point
      !open(unit_fci_1,File=trim(file_name_1),STATUS='NEW',FORM='FORMATTED')
      write(file_name_1,101) i_q_point, i_k_point, i_k_point, i_q_point
      write(file_name_2,102) i_q_point, i_k_point, i_k_point, i_q_point 
      open(unit_fci_1,File=trim(file_name_1),STATUS='NEW',FORM='FORMATTED')
      open(unit_fci_2,File=trim(file_name_2),STATUS='NEW',FORM='FORMATTED')
      write(unit_fci_1,104) i_q_point, i_k_point, i_k_point, i_q_point
      write(unit_fci_1,'(2X,I4)') n_states
      write(unit_fci_2,104) i_q_point, i_k_point, i_k_point, i_q_point 
      write(unit_fci_2,'(2X,I4)') n_states
      do i_state_1 = 1, n_states, 1
        do i_state_2 = 1, n_states, 1
          do i_state_3 = 1, n_states, 1
            do i_state_4 = 1, n_states, 1
              integral_4ks_tmp = (0.0d0,0.0d0)
              do i_prodbas_1 = 1, n_basbas, 1
                do i_prodbas_2 = 1, n_basbas, 1
                  integral_4ks_tmp = integral_4ks_tmp + &
                      lvl_tricoeff_kq(i_prodbas_1,i_state_1,i_state_2,1) * &
                      coulomb_tmp(i_prodbas_1, i_prodbas_2) * &
                      conjg(lvl_tricoeff_kq(i_prodbas_2,i_state_4,i_state_3,1))
                ! end of i_prodbas_1
                enddo
              ! end of i_prodbas_2
              enddo
              if (abs(integral_4ks_tmp) > Threshold_Int) then
                  !write(unit_fci_1,105) integral_4ks_tmp, &
                  !    i_state_1, i_state_2, i_state_3, i_state_4
                  write(unit_fci_1,103) real(integral_4ks_tmp), &
                      i_state_1, i_state_2, i_state_3, i_state_4
                  write(unit_fci_2,103) aimag(integral_4ks_tmp), &
                      i_state_1, i_state_2, i_state_3, i_state_4
              endif
              if (i_state_1 .eq. i_state_4 .and. i_state_1 .le. n_homo_k(i_q_point,1) .and. &
                  i_state_2 .eq. i_state_3 .and. i_state_2 .le. n_homo_k(i_k_point,1)) then
                  exx_energy = exx_energy + &
                      abs(integral_4ks_tmp) * &
                      k_weights(i_irkq_point) * k_weights(i_q_point)
                  !occ_numbers(i_state,i_spin,i_k_point)*k_p_weights(i_k_point)
              endif
            ! end of i_state_4
            enddo
          ! end of i_state_3
          enddo
        ! end of i_state_2
        enddo
      ! end of i_state_1
      enddo
      close(unit_fci_1)
      close(unit_fci_2)

      !call get_times(time_polar_3, clock_time_polar_3)
      !tot_time_polar_3 = tot_time_polar_3 + time_polar_3
      !tot_clock_time_polar_3 = tot_clock_time_polar_3 + clock_time_polar_3

     ! end loop over i_k_point
     enddo
   ! end loop i_q_point
   enddo

   call sync_real_number(exx_energy)
   if (myid .eq. 0) then
       exx_energy = - exx_energy
       write(use_unit,'(2X,A,F16.8,A)') 'Exx energy :',exx_energy,' a.u.'
   endif

   deallocate(n_homo_k)
   deallocate(n_lumo_k)
   deallocate(coulomb_tmp)
   deallocate(lvl_tricoeff_recip_q)
   deallocate(KS_eigenvector_q)
   deallocate(KS_eigenvector_k)
   deallocate(lvl_tricoeff_sum)
   deallocate(lvl_tricoeff_sum_all)
   deallocate(coeff_tmp)
   deallocate(coeff_tmp_2)
   deallocate(lvl_tricoeff_kq)
   deallocate(lvl_tricoeff_tmp)
   deallocate(lvl_tricoeff_tmp_2)

   return

end subroutine get_4_index_integrals_slab
!---------------------------------------------------------------------
!******
