!****s* FHI-aims/evaluate_exchange_matr_kspace_p0.f90
!  NAME evaluate_exchange_matr_kspace_p0
!   
!  SYNOPSIS

subroutine evaluate_exchange_matr_kspace_p0 &
(KS_eigenvalue,KS_eigenvector,KS_eigenvector_complex,occ_numbers)

  !  PURPOSE
  !  Subroutine evaluate_exchange_matr_kspace_p0 evaluates the exact-exchange part of
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
  use constants
  use basis
  use localorb_io, only: localorb_info, use_unit
  use geometry, only: species
  implicit none

  !  ARGUMENTS

  real*8, dimension(n_states,n_spin,n_k_points_task) :: KS_eigenvalue
  real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
  complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
  real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers

  !  INPUTS
  !  o  occ_numbers -- real array,
  !       the occupation number of the electrons for each eigenstate and each spin
  !  o  KS_eigenvector -- real array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  o  KS_eigenvector_complex -- complex array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  OUTPUTS
  !  none
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

  complex*16, allocatable :: lvl_tricoeff_recip_tmp(:,:,:) ! LVL triple coefficients in k space
  complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector a a given (romote) k point
  complex*16, allocatable :: lvl_tricoeff_sum1(:,:,:) ! sum of LVL triple coefficient at k and q points (atom 1)
  complex*16, allocatable :: lvl_tricoeff_sum2(:,:,:) ! sum of LVL triple coefficient at k and q points (atom 2)
  complex*16, allocatable :: lvl_tricoeff_MO(:,:,:,:) ! Full LVL triple coefficient, one basis index is transformed to
                                                      ! (occupied) molecular orbital index 
  complex*16, allocatable :: tmp_MO(:,:,:)
!  complex*16, allocatable :: coulomb_times_tricoeff(:,:) ! Coulomb matrix multiplies the LVL triple coefficient
  complex*16:: coulomb_times_tricoeff(n_basbas,n_basis) ! Coulomb matrix multiplies the LVL triple coefficient
  complex*16, allocatable :: exchange_matr_tmp(:,:) ! temparary exchange matrix per k_point
  complex*16, allocatable :: coulomb_matr_recip_tmp(:,:) !
  complex*16, allocatable :: ex_matr_complex(:,:) !
  real*8, allocatable :: ex_matr_real(:,:) !

  real*8  :: ex_vect_real
  complex*16  :: ex_vect_complex
  real*8  :: en_exx

  integer :: info, mpierr
  character(*), parameter :: func = 'evaluate_exchange_matr_kspace_p0'
  character*150 :: info_str

  integer :: max_n_homo
  integer :: n_homo_k(n_k_points,n_spin)
! counter 
  integer i_k_point
  integer i_q_point
  integer i_kq_point
  integer i_k_point_local
  integer i_q_point_local
  integer i_state
  integer i_state_1
  integer i_basis_1
  integer i_spin
  integer id_root
  integer i_task
  integer, allocatable :: k_points_at_task(:,:)
  integer :: i_req(0:n_tasks-1)
  integer :: i_atom_1, i_atom_2
  integer :: i_species_1, i_species_2
  integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
  integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
  integer :: i_1, i_2

print*,'basbas',n_basis,n_basbas,n_states
!stop
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  if(use_scalapack) call aims_stop('*** Periodic EXX in k-space must not be used if use_scalapack is in effect!')
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  n_homo_k(:,:) = 0
  do i_spin = 1, n_spin, 1
    do i_k_point = 1, n_k_points, 1
      do i_state = 1, n_states
        if(occ_numbers(i_state,i_spin,i_k_point) .gt. 1.e-12) then
         n_homo_k(i_k_point,i_spin) = i_state
        endif
      enddo
 !     if(myid.eq.0) then
 !       write(use_unit,*) i_spin, i_k_point, n_homo_k(i_k_point, i_spin), n_homo(i_spin)
 !     endif
    enddo
  enddo
  max_n_homo = max(maxval(n_homo_k(:,1), n_k_points), &
                   maxval(n_homo_k(:,n_spin), n_k_points))
 ! write(use_unit,*)"max_n_homo:", max_n_homo

!  max_n_homo = max(n_homo(1),n_homo(n_spin))

  allocate(lvl_tricoeff_recip_tmp(max_n_basbas_sp,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_recip_tmp', func)
  allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
  call check_allocation(info, 'KS_eigenvector_q', func)
  allocate(lvl_tricoeff_MO(n_basbas,n_basis,max_n_homo,n_spin),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_MO', func)
!  allocate(coulomb_times_tricoeff(n_basbas,n_basis),stat=info) 
!  call check_allocation(info, 'coulomb_times_tricoeff', func)
  allocate(exchange_matr_tmp(n_basis,n_basis),stat=info) 
  call check_allocation(info, 'exchange_matr_tmp', func)
  allocate(coulomb_matr_recip_tmp(n_basbas,n_basbas),stat=info) 
  call check_allocation(info, 'coulomb_matr_recip_tmp', func)

  if(real_eigenvectors) then
    allocate(ex_matr_real(n_basis,n_states),stat=info)
    call check_allocation(info, 'ex_matr_real', func)
    allocate(ex_matr_complex(1,1),stat=info)
  else
    allocate(ex_matr_complex(n_basis,n_states),stat=info)
    call check_allocation(info, 'ex_matr_complex', func)
    allocate(ex_matr_real(1,1),stat=info)
  endif

  write(info_str,'(2X,A)') "Constructing the exchange matrix in k space ..."
  call localorb_info(info_str)

  ! k_points_at_task: which task contains which k-point

  allocate(k_points_at_task(0:n_tasks-1,(n_k_points-1)/n_tasks+1))
  k_points_at_task(:,:) = 0
  do i_k_point = 1, n_k_points, 1
     i_k_point_local = (i_k_point-1)/n_tasks + 1 
     k_points_at_task(mod(i_k_point,n_tasks),i_k_point_local) = i_k_point
  enddo

  i_req(:) = MPI_REQUEST_NULL

  if(real_eigenvectors) then
    hf_exchange_matr_real(:,:,:,:) = 0.d0
  else
    hf_exchange_matr_complex(:,:,:,:) = (0.d0,0.d0)
  endif
  en_exx = 0.d0
  do i_q_point = 1, n_k_points, 1
!     if(myid.eq.1) print*,i_q_point
     id_root = mod(i_q_point,n_tasks)
     i_q_point_local = (i_q_point-1)/n_tasks + 1
     if(myid .eq. id_root) then

        lvl_tricoeff_recip_tmp(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_q_point_local)

        if(real_eigenvectors) then
           KS_eigenvector_q(:,:,:)=dble(KS_eigenvector(:,:,:,i_q_point_local))
        else
           KS_eigenvector_q(:,:,:)=KS_eigenvector_complex(:,:,:,i_q_point_local)
        endif
     endif

     call mpi_bcast(lvl_tricoeff_recip_tmp,max_n_basbas_sp*n_basis*n_basis, &
                             MPI_COMPLEX16, id_root, mpi_comm_global, mpierr)

     call mpi_bcast(KS_eigenvector_q,n_basis*n_states*n_spin, &
                             MPI_COMPLEX16, id_root, mpi_comm_global, mpierr)
     KS_eigenvector_q = conjg(KS_eigenvector_q)

     do i_k_point_local = 1, (n_k_points-1)/n_tasks + 1

        ! Send our coulomb_matr_recip to the task(s) which actually needs it
        ! Normally, only one send should be necessary

        do i_task = 0, n_tasks-1
           i_k_point = k_points_at_task(i_task, i_k_point_local)
           if(i_k_point == 0) cycle
           i_kq_point = kq_point_list(i_k_point,i_q_point)
           if(mod(i_kq_point,n_tasks) == myid) then
             ! Tasks i_task needs i_kq_point from me
             call mpi_isend(coulomb_matr_recip(1,1,(i_kq_point-1)/n_tasks+1), n_basbas*n_basbas, MPI_COMPLEX16, &
                            i_task, 111, mpi_comm_global, i_req(i_task), mpierr)
           endif
        enddo

        ! Our k-point in this turn:
        i_k_point = k_points_at_task(myid, i_k_point_local)

        if(i_k_point == 0) then
           call mpi_waitall(n_tasks, i_req, MPI_STATUSES_IGNORE, mpierr)
           cycle
        endif

        i_kq_point = kq_point_list(i_k_point,i_q_point)

        call mpi_recv(coulomb_matr_recip_tmp, n_basbas*n_basbas, MPI_COMPLEX16, &
                     mod(i_kq_point, n_tasks), 111, mpi_comm_global, MPI_STATUS_IGNORE, mpierr)

        call mpi_waitall(n_tasks, i_req, MPI_STATUSES_IGNORE, mpierr)


        lvl_tricoeff_MO = 0

        do i_atom_1 = 1, n_atoms, 1
           i_species_1 = species(i_atom_1)
           basis_off_1 = atom2basis_off(i_atom_1)
           n_sp_basis_1 = sp2n_basis_sp(i_species_1)
           bboff_1 = atom2basbas_off(i_atom_1)
           n_spbb_1 = sp2n_basbas_sp(i_species_1)

           allocate(lvl_tricoeff_sum1(n_spbb_1,n_sp_basis_1,n_basis),stat=info)
           call check_allocation(info, 'lvl_tricoeff_sum1', func)

           do i_1 = 1, n_sp_basis_1
              do i_2 = 1, n_basis
                 lvl_tricoeff_sum1(1:n_spbb_1,i_1,i_2) = &
                    conjg(lvl_tricoeff_recip_tmp(1:n_spbb_1, basis_off_1+i_1, i_2)) + &
                    lvl_tricoeff_recip2(1:n_spbb_1, i_2, basis_off_1+i_1)
              enddo
           enddo

           ! Multiply lvl_tricoeff_sum1 with KS_eigenvector_q
           ! Please note that the multipication uses the first 2 dimensions of
           ! lvl_tricoeff_sum1/tmp_MO as 1 logical dimension, thus these 2 arrays
           ! must be allocated to the exact size!

           allocate(tmp_MO(n_spbb_1,n_sp_basis_1,max_n_homo),stat=info)
           call check_allocation(info, 'tmp_MO', func)

           do i_spin = 1, n_spin
              call zgemm('N', 'N', n_spbb_1*n_sp_basis_1, max_n_homo, n_basis, (1.d0,0.d0), &
                       lvl_tricoeff_sum1, n_spbb_1*n_sp_basis_1, &
                       KS_eigenvector_q(1,1,i_spin), n_basis, (0.d0,0.d0), &
                       tmp_MO, n_spbb_1*n_sp_basis_1)
              lvl_tricoeff_MO(bboff_1+1:bboff_1+n_spbb_1,basis_off_1+1:basis_off_1+n_sp_basis_1,:,i_spin) = tmp_MO
           enddo
           deallocate(tmp_MO)

           deallocate(lvl_tricoeff_sum1)
        enddo

        do i_atom_2 = 1, n_atoms, 1
           i_species_2 = species(i_atom_2)
           basis_off_2 = atom2basis_off(i_atom_2)
           n_sp_basis_2 = sp2n_basis_sp(i_species_2)
           bboff_2 = atom2basbas_off(i_atom_2)
           n_spbb_2 = sp2n_basbas_sp(i_species_2)

           allocate(lvl_tricoeff_sum2(n_spbb_2,n_basis,n_sp_basis_2),stat=info)
           call check_allocation(info, 'lvl_tricoeff_sum2', func)

           do i_1 = 1, n_basis
              do i_2 = 1, n_sp_basis_2
                 lvl_tricoeff_sum2(1:n_spbb_2,i_1,i_2) = &
                    lvl_tricoeff_recip1(1:n_spbb_2, basis_off_2+i_2, i_1, i_k_point_local) + &
                    lvl_tricoeff_recip2(1:n_spbb_2, i_1, basis_off_2+i_2)
              enddo
           enddo

           ! Multiply lvl_tricoeff_sum2 with KS_eigenvector_q
           ! See remark above about matrix multiplication

           allocate(tmp_MO(n_spbb_2,n_basis,max_n_homo))
           call check_allocation(info, 'tmp_MO', func)

           do i_spin = 1, n_spin
              call zgemm('N', 'N', n_spbb_2*n_basis, max_n_homo, n_sp_basis_2, (1.d0,0.d0), &
                       lvl_tricoeff_sum2, n_spbb_2*n_basis, &
                       KS_eigenvector_q(basis_off_2+1,1,i_spin), n_basis, (0.d0,0.d0), &
                       tmp_MO, n_spbb_2*n_basis)
              lvl_tricoeff_MO(bboff_2+1:bboff_2+n_spbb_2,:,:,i_spin) = &
              lvl_tricoeff_MO(bboff_2+1:bboff_2+n_spbb_2,:,:,i_spin) + tmp_MO
           enddo
           deallocate(tmp_MO)

           deallocate(lvl_tricoeff_sum2)
        enddo

        do i_spin = 1, n_spin, 1
!           !$OMP PARALLEL DO private(i_state,coulomb_times_tricoeff,exchange_matr_tmp)
           do i_state = 1, max_n_homo, 1
                call perfon('eemk')

                call zgemm('N', 'N', n_basbas, n_basis, n_basbas, (1.d0,0.d0), &
                        coulomb_matr_recip_tmp, n_basbas, &
                        lvl_tricoeff_MO(1,1,i_state,i_spin), n_basbas, (0.d0,0.d0), &
                        coulomb_times_tricoeff, n_basbas)

                call zgemm('C', 'N', n_basis, n_basis, n_basbas, (1.d0,0.d0), &
                         lvl_tricoeff_MO(1,1,i_state,i_spin), n_basbas, &
                         coulomb_times_tricoeff, n_basbas, (0.d0,0.d0), &
                         exchange_matr_tmp,n_basis)
                 call perfoff

              if(real_eigenvectors) then
                 hf_exchange_matr_real(:,:,i_k_point_local,i_spin) = &
                   hf_exchange_matr_real(:,:,i_k_point_local,i_spin) + &
                   exchange_matr_tmp(:,:)*k_weights(i_q_point)* & 
                   occ_numbers(i_state,i_spin,i_q_point)*dble(n_spin)/2.d0

                   call dgemm('N', 'N', n_basis, n_states, n_basis, 1.d0, &
                              real(exchange_matr_tmp), n_basis, &
                              KS_eigenvector(:,:,i_spin,i_k_point_local), n_basis, 0.d0, &
                              ex_matr_real,n_basis)

                   do i_state_1 = 1, max_n_homo, 1
                      ex_vect_real = 0.d0
                      do i_basis_1 = 1, n_basis, 1
                        ex_vect_real = ex_vect_real + &
                                       KS_eigenvector(i_basis_1,i_state_1,i_spin,i_k_point_local) &
                                     * ex_matr_real(i_basis_1,i_state_1)
                      enddo
                      if(KS_eigenvalue(i_state,i_spin,i_q_point) .gt.  &
                          KS_eigenvalue(i_state_1,i_spin,i_k_point) ) then
                         en_exx = en_exx + ex_vect_real * occ_numbers(i_state,i_spin,i_q_point) * &
                                    k_weights(i_q_point)*k_weights(i_k_point) * dble(n_spin)
                               

                      elseif(KS_eigenvalue(i_state,i_spin,i_q_point) .eq.  KS_eigenvalue(i_state_1,i_spin,i_k_point))  then
                            en_exx = en_exx + ex_vect_real * occ_numbers(i_state,i_spin,i_q_point) * &
                                    k_weights(i_q_point)*k_weights(i_k_point) *dble(n_spin)/2.d0
                      endif
                                               
                   enddo

              else
                 hf_exchange_matr_complex(:,:,i_k_point_local,i_spin) = &
                   hf_exchange_matr_complex(:,:,i_k_point_local,i_spin) + &
                   exchange_matr_tmp(:,:)*k_weights(i_q_point)* & 
                   occ_numbers(i_state,i_spin,i_q_point)*dble(n_spin)/2.d0

                   call zgemm('N', 'N', n_basis, n_states, n_basis, (1.d0,0.d0), &
                              exchange_matr_tmp, n_basis, &
                              KS_eigenvector_complex(:,:,i_spin,i_k_point_local), n_basis, (0.d0,0.d0), &
                              ex_matr_complex,n_basis)

                   do i_state_1 = 1, max_n_homo, 1
                      ex_vect_complex = (0.d0,0.d0)
                      do i_basis_1 = 1, n_basis, 1
                        ex_vect_complex = ex_vect_complex + &
                                 conjg(KS_eigenvector_complex(i_basis_1,i_state_1,i_spin,i_k_point_local)) &
                                 * ex_matr_complex(i_basis_1,i_state_1) 
                      enddo
                      if(KS_eigenvalue(i_state,i_spin,i_q_point) .gt. &
                          KS_eigenvalue(i_state_1,i_spin,i_k_point) ) then
                         en_exx = en_exx + ex_vect_complex * occ_numbers(i_state,i_spin,i_q_point) * &
                                    k_weights(i_q_point)*k_weights(i_k_point) * dble(n_spin)

                      elseif(KS_eigenvalue(i_state,i_spin,i_q_point) ==   &
                                 KS_eigenvalue(i_state_1,i_spin,i_k_point) ) then
                         en_exx = en_exx + ex_vect_complex * occ_numbers(i_state,i_spin,i_q_point) * &
                                    k_weights(i_q_point)*k_weights(i_k_point) * dble(n_spin)/2.d0 
                      endif
                                               
                   enddo
              endif

! end loop over i_state
           enddo
!           !$OMP END PARALLEL DO
! end loop over i_spin
        enddo
! end loop over i_k_point
     enddo
           
! end loop i_q_point
  enddo
  call sync_real_number(en_exx)

  if(myid.eq.0) then
   write(use_unit,*) "EN_EXX:", en_exx, en_exx*hartree
  endif

  deallocate(lvl_tricoeff_recip_tmp)
  deallocate(KS_eigenvector_q)
  deallocate(lvl_tricoeff_MO)
!  deallocate(coulomb_times_tricoeff)
  deallocate(exchange_matr_tmp)
  deallocate(coulomb_matr_recip_tmp)
  deallocate(ex_matr_real)
  deallocate(ex_matr_complex)

  deallocate(k_points_at_task)

  return

end subroutine evaluate_exchange_matr_kspace_p0
!---------------------------------------------------------------------
!******
