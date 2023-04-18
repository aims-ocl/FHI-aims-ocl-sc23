
!****s* FHI-aims/evaluate_exchange_matr_kspace_single_kpoint_p0.f90
!  NAME evaluate_exchange_matr_kspace_single_kpoint_p0
!   
!  SYNOPSIS

subroutine evaluate_exchange_matr_kspace_single_kpoint_p0 &
     (k_point,KS_egnv,KS_egnv_complex,real_eigen,occ_numbers,&
     q_weights,lvl_tricoeff_recip1_new,lvl_tricoeff_recip2_new,fock_m)

  !  PURPOSE
  !  Subroutine evaluate_exchange_matr_kspace_single_kpoint_p0 evaluates the exact-exchange part of
  !  Hartree-Fock hamiltonian in a periodic system for a set of k-points on a band segment. 
  !  The algorithm used here is based
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
  use localorb_io
  use constants
  use basis
  use geometry, only: species
  implicit none

  !  ARGUMENTS

  integer :: k_point
  real*8, dimension(:,:,:,:) :: KS_egnv
  complex*16, dimension(:,:,:,:) :: KS_egnv_complex
  logical :: real_eigen
  real*8, dimension(:,:,:) :: occ_numbers
  real*8, dimension(:) :: q_weights
  complex*16, intent(IN)  :: lvl_tricoeff_recip1_new(:,:,:,:)
  real*8, intent(IN)      :: lvl_tricoeff_recip2_new(:,:,:)
  complex*16, intent(OUT) :: fock_m(:,:)
  

  !  INPUTS
  !
  !  o  k_point -- number of the k-point within n_k_points at which the exchange matrix is calculated
  !  o  occ_numbers -- real array,
  !       the occupation number of the electrons for each eigenstate and each spin
  !  o  q_weights -- weights of q-points
  !  o  KS_egnv -- real array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  o  KS_egnv_complex -- complex array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  o  real_eigen -- true if eigenvectors are real, false otherwise
  !  o  lvl_tricoeff_recip1_new and lvl_tricoeff_recip2_new -- parts of the lvl triple
  !            coeffs calculated at the set of n_k_points
  !
  !  OUTPUTS
  !
  !  o  fock_m -- the Fock matrix at k-point k_point
  !  
  !  
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
  complex*16, allocatable :: coulomb_times_tricoeff(:,:) ! Coulomb matrix multiplies the LVL triple coefficient
  complex*16, allocatable :: exchange_matr_tmp(:,:,:) ! temparary exchange matrix per k_point
!  complex*16, allocatable :: coulomb_matr_recip_tmp(:,:) !
  integer :: info, mpierr
  character(*), parameter :: func = 'evaluate_exchange_matr_kspace_single_kpoint_p0'
  character*150 :: info_str

  integer :: max_n_homo
  integer :: n_homo_k(n_q_points,n_spin)
! counter 
  integer i_k_point
  integer i_q_point
  integer i_kq_point
  integer i_k_point_local
  integer i_q_point_local
  integer i_state
  integer i_spin
  integer id_root
  integer i_task
  integer k_point_task
  integer, allocatable :: k_points_at_task(:,:)
  integer :: i_req(0:n_tasks-1)
  integer :: i_atom_1, i_atom_2
  integer :: i_species_1, i_species_2
  integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
  integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
  integer :: i_1, i_2

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  if(use_scalapack) call aims_stop('*** Periodic EXX in k-space must not be used if use_scalapack is in effect!')
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------


  n_homo_k(:,:) = 0
  do i_spin = 1, n_spin, 1
    do i_k_point = 1, n_q_points, 1
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
  max_n_homo = max(maxval(n_homo_k(:,1), n_q_points), &
                   maxval(n_homo_k(:,n_spin), n_q_points))
 ! write(use_unit,*)"max_n_homo:", max_n_homo

!  max_n_homo = max(n_homo(1),n_homo(n_spin))

  allocate(lvl_tricoeff_recip_tmp(max_n_basbas_sp,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_recip_tmp', func)
  allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
  call check_allocation(info, 'KS_eigenvector_q', func)
  allocate(lvl_tricoeff_MO(n_basbas,n_basis,max_n_homo,n_spin),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_MO', func)
  allocate(coulomb_times_tricoeff(n_basbas,n_basis),stat=info) 
  call check_allocation(info, 'coulomb_times_tricoeff', func)
  allocate(exchange_matr_tmp(n_basis,n_basis,n_spin),stat=info) 
  call check_allocation(info, 'exchange_matr_tmp', func)
!  allocate(coulomb_matr_recip_tmp(n_basbas,n_basbas),stat=info) 
!  call check_allocation(info, 'coulomb_matr_recip_tmp', func)

  if(output_priority .le. OL_norm) then
     write(info_str,'(2X,A)') "Constructing the exchange matrix in k space ..."
     call localorb_info(info_str)
  endif

  ! k_points_at_task: which task contains which k-point

  allocate(k_points_at_task(0:n_tasks-1,(n_q_points-1)/n_tasks+1))
  k_points_at_task(:,:) = 0
  do i_k_point = 1, n_q_points, 1
     i_k_point_local = (i_k_point-1)/n_tasks + 1 
     k_points_at_task(mod(i_k_point,n_tasks),i_k_point_local) = i_k_point
  enddo

  k_point_task = mod(k_point,n_tasks)
  i_k_point_local = (k_point-1)/n_tasks + 1
  if (myid .eq. k_point_task) then
     lvl_tricoeff_recip_tmp(:,:,:)=lvl_tricoeff_recip1_new(:,:,:,i_k_point_local)
  endif
  call mpi_bcast(lvl_tricoeff_recip_tmp,max_n_basbas_sp*n_basis*n_basis, &
       MPI_COMPLEX16, k_point_task, mpi_comm_global, mpierr)

  i_req(:) = MPI_REQUEST_NULL


  exchange_matr_tmp(:,:,:) = (0.d0,0.d0)

  do i_q_point_local = 1, n_q_points_task, 1

     i_q_point = k_points_at_task(myid,i_q_point_local)

     if(real_eigen) then
        KS_eigenvector_q(:,:,:)=KS_egnv(:,:,:,i_q_point_local)
     else
        KS_eigenvector_q(:,:,:)=KS_egnv_complex(:,:,:,i_q_point_local)
     endif

     KS_eigenvector_q = conjg(KS_eigenvector_q)

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
                   conjg(lvl_tricoeff_recip1(1:n_spbb_1, basis_off_1+i_1, i_2, i_q_point_local)) + &
                   lvl_tricoeff_recip2_new(1:n_spbb_1, i_2, basis_off_1+i_1)
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
                   lvl_tricoeff_recip_tmp(1:n_spbb_2, basis_off_2+i_2, i_1) + &
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
        do i_state = 1, max_n_homo, 1
           call zgemm('N', 'N', n_basbas, n_basis, n_basbas, (1.d0,0.d0), &
             coulomb_matr_recip(:,:,i_q_point_local), n_basbas, &
             lvl_tricoeff_MO(1,1,i_state,i_spin), n_basbas, (0.d0,0.d0), &
             coulomb_times_tricoeff, n_basbas)
           
           call zgemm('C', 'N', n_basis, n_basis, n_basbas, q_weights(i_q_point)* &
             occ_numbers(i_state,i_spin,i_q_point)*dble(n_spin)/2.d0*(1.d0,0.d0), &
             lvl_tricoeff_MO(1,1,i_state,i_spin), n_basbas, &
             coulomb_times_tricoeff, n_basbas, (1.d0,0.d0), &
             exchange_matr_tmp(:,:,i_spin),n_basis)
           
           ! end loop over i_state
        enddo
        ! end loop over i_spin
     enddo
     ! end loop i_q_point
  enddo
  
  fock_m = (0d0,0d0)
  do i_spin = 1, n_spin
     i_task = 0
     do i_k_point = 1, n_basis
        do i_q_point = 1, i_k_point
           i_task = i_task + 1
           fock_m(i_task,i_spin) = exchange_matr_tmp(i_q_point,i_k_point,i_spin)
        enddo
     enddo
  enddo
  call sync_vector_complex(fock_m,n_basis*(n_basis+1)/2*n_spin)

  deallocate(lvl_tricoeff_recip_tmp)
  deallocate(KS_eigenvector_q)
  deallocate(lvl_tricoeff_MO)
  deallocate(coulomb_times_tricoeff)
  deallocate(exchange_matr_tmp)

  deallocate(k_points_at_task)

  return

end subroutine evaluate_exchange_matr_kspace_single_kpoint_p0


module ex_mat_ksk_p0
  use dimensions
  use crpa_blacs
  use exchange_trico
  use exchange_ev
  implicit none

  integer, dimension(:), allocatable:: lb_atom, ub_atom
  
contains
  
  !---------------------------------------------------------------------
  !******
  !****s* FHI-aims/evaluate_exchange_matr_kspace_single_kpoint_p0.f90
  !  NAME evaluate_exchange_matr_kspace_single_kpoint_p0
  !   
  !  SYNOPSIS
  
  subroutine my_evaluate_exchange_matr_kspace_single_kpoint_p0 &
       (n_k_points_band, n_ks_points_band_task,k_point,KS_egnv,KS_egnv_complex,real_eigen,occ_numbers,&
       q_weights,lvl_tricoeff_recip_r_k,fock_m)
    
    !  PURPOSE
    !  Subroutine evaluate_exchange_matr_kspace_single_kpoint_p0 evaluates the exact-exchange part of
    !  Hartree-Fock hamiltonian in a periodic system for a set of k-points on a band segment. 
    !  The algorithm used here is based
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
    use localorb_io
    use constants
    use crpa_blacs
    use lvl_tricoeff
    implicit none
    
    !  ARGUMENTS
    
    integer :: n_k_points_band, n_ks_points_band_task, k_point
    real*8, dimension(1,1,1,1) :: KS_egnv
    complex*16, dimension(n_basis,n_states,n_spin,n_q_points_task) :: KS_egnv_complex
    logical :: real_eigen
    real*8, dimension(n_states,n_spin,n_q_points) :: occ_numbers
    real*8, dimension(n_q_points) :: q_weights
    complex*16, intent(IN)  :: lvl_tricoeff_recip_r_k(lbb_row:ubb_row,max_n_basis_sp,n_basis)
    
    complex*16, intent(OUT) :: fock_m(n_basis*(n_basis+1)/2,n_spin)
    
    
    !  INPUTS
    !
    !  o  k_point -- number of the k-point within n_k_points at which the exchange matrix is calculated
    !  o  occ_numbers -- real array,
    !       the occupation number of the electrons for each eigenstate and each spin
    !  o  q_weights -- weights of q-points
    !  o  KS_egnv -- real array,
    !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
    !  o  KS_egnv_complex -- complex array,
    !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
    !  o  real_eigen -- true if eigenvectors are real, false otherwise
    !  o  lvl_tricoeff_recip1_new and lvl_tricoeff_recip2_new -- parts of the lvl triple
    !            coeffs calculated at the set of n_k_points
    !
    !  OUTPUTS
    !
    !  o  fock_m -- the Fock matrix at k-point k_point
    !  
    !  
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
    
    complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector a a given (romote) k point
    
    complex*16, dimension(:,:,:), allocatable, target:: lvl_tricoeff_MO_r, lvl_tricoeff_MO_c_arr 
    complex*16, dimension(:,:,:), pointer:: lvl_tricoeff_MO_c
    
    complex*16, allocatable :: coulomb_times_tricoeff(:,:) ! Coulomb matrix multiplies the LVL triple coefficient
    complex*16, allocatable :: exchange_matr_tmp(:,:,:) ! temparary exchange matrix per k_point

    integer :: info, mpierr
    character(*), parameter :: func = 'evaluate_exchange_matr_kspace_p0'
    character*150 :: info_str
    
    integer :: max_n_homo
    integer :: n_homo_k(n_q_points,n_spin)
    ! counter 
    integer i_k_point
    integer i_q_point
    integer i_q_point_local
    integer max_q_points_task
    integer i_state
    integer i_spin
    integer i_task

    integer :: win_ev, win_tri_k
    integer :: count
    integer:: oldatom, i_basis, i_atom

!    integer:: mirror
!    integer:: status(MPI_STATUS_SIZE)
    integer:: blockrest, blocksize, n_blocks, n_thisblock, i_block, i_state_bl


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !the use_scalapack flag has no effect here, this always uses lapack (with the according k parallelization)

    allocate(lb_atom(n_atoms),ub_atom(n_atoms))
    oldatom=-1
    do i_basis = 1, n_basis
       i_atom = Cbasis_to_atom(i_basis)
       if(i_atom.ne.oldatom) then
          lb_atom(i_atom)=i_basis
          if (i_atom.gt.1) ub_atom(i_atom-1)=i_basis-1
          oldatom=i_atom
       end if
    end do
    ub_atom(n_atoms)=n_basis
    
    n_homo_k(:,:) = 0
    do i_spin = 1, n_spin, 1
       do i_q_point = 1, n_q_points, 1
          do i_state = 1, n_states
             if(occ_numbers(i_state,i_spin,i_q_point) .gt. 1.e-12) then
                n_homo_k(i_q_point,i_spin) = i_state
             endif
          enddo
       enddo
    enddo
    max_n_homo = max(maxval(n_homo_k(:,1), n_q_points), &
         maxval(n_homo_k(:,n_spin), n_q_points))


    !limit the memory consumption by the lvl_tricoeff_MO_r arrays to 800 MB per task
    blocksize=min(50000000/(3*bb_bl_row*n_basis),max_n_homo)
    !    blocksize=4

    n_blocks=max_n_homo/blocksize
    blockrest=mod(max_n_homo,blocksize)
    if (blockrest.gt.0) n_blocks=n_blocks+1
!    if (myid.eq.0) print*,'ojo emblocks',blocksize,max_n_homo, n_blocks


    ! write(use_unit,*)"max_n_homo:", max_n_homo

    !  max_n_homo = max(n_homo(1),n_homo(n_spin))

    allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
    call check_allocation(info, 'KS_eigenvector_q', func)
    allocate(lvl_tricoeff_MO_r(lbb_row:ubb_row,n_basis,blocksize),stat=info) 
    call check_allocation(info, 'lvl_tricoeff_MO_r', func)

    !use new array for lvl_tricoeff_MO_c for the tasks where it is not identical with lvl_tricoeff_MO_r
    if (myid_row.eq.myid_col) then
       lvl_tricoeff_MO_c => lvl_tricoeff_MO_r
    else
       allocate(lvl_tricoeff_MO_c_arr(lbb_col:ubb_col,n_basis,blocksize))
       lvl_tricoeff_MO_c => lvl_tricoeff_MO_c_arr
    end if

    allocate(coulomb_times_tricoeff(n_bb_row,n_basis),stat=info) 
    call check_allocation(info, 'coulomb_times_tricoeff', func)
    allocate(exchange_matr_tmp(n_basis,n_basis,n_spin),stat=info) 
    call check_allocation(info, 'exchange_matr_tmp', func)

    
    if(output_priority .le. OL_norm) then
       write(info_str,'(2X,A)') "Constructing the exchange matrix in k space ..."
       call localorb_info(info_str)
    endif
    
    
    if (real_eigenvectors) then
       call init_access_ev_real(n_q_points,n_q_points_task,KS_egnv,win_ev)
    else
       call init_access_ev_complex(n_q_points,n_q_points_task,KS_egnv_complex,win_ev)
    end if

    !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the irkq loop
    if (mod(n_q_points,n_tasks_irkq).eq.0) then    
       max_q_points_task=n_q_points/n_tasks_irkq
    else
       max_q_points_task=n_q_points/n_tasks_irkq+1
    end if
    
    exchange_matr_tmp(:,:,:) = (0.d0,0.d0)
    
    
    do i_q_point_local = 1, max_q_points_task
       i_q_point=(i_q_point_local-1)*n_tasks_irkq+myid_irkq+1

       if (irkblacs_member.and.(i_q_point.le.n_q_points)) then
          call access_ev(win_ev,k_point_loc(:,i_q_point),KS_eigenvector_q)

          call perfon('eemks')             
          do i_spin = 1, n_spin
             do i_block = 1, n_blocks
                n_thisblock=min(blocksize,max_n_homo-(i_block-1)*blocksize)
                if(myid_col.eq.0) then
                   call perfon('eemco')
                   do i_state_bl = 1, n_thisblock
                      i_state=(i_block-1)*blocksize+i_state_bl
                      call comp_tricoeffsum(lvl_tricoeff_recip_r_k, lvl_tricoeff_mod_r(:,:,i_state,i_spin,i_q_point_local), &
                           KS_eigenvector_q(:,i_state,i_spin), lvl_tricoeff_MO_r(:,:,i_state_bl), lbb_row, ubb_row)
                   end do
                   call perfoff
                end if
                if(n_tasks_bl.gt.1) then
                   count=n_bb_row*n_basis*n_thisblock
                   call mpi_bcast(lvl_tricoeff_MO_r,count,MPI_DOUBLE_COMPLEX,0,comm_blacs_col,mpierr)
                   count=n_bb_col*n_basis*n_thisblock
                   call mpi_bcast(lvl_tricoeff_MO_c,count,MPI_DOUBLE_COMPLEX,myid_col,comm_blacs_row,mpierr)
                end if
                
                call perfon('eemmu')
                do i_state_bl = 1, n_thisblock
                   i_state=(i_block-1)*blocksize+i_state_bl
                   call mult_with_tricoeffsum(q_weights(i_q_point), occ_numbers(i_state,i_spin,i_q_point), &
                        coulomb_matr_blacs(:,:,i_q_point_local), lvl_tricoeff_MO_r(:,:,i_state_bl), lvl_tricoeff_MO_c(:,:,i_state_bl), &
                        coulomb_times_tricoeff, exchange_matr_tmp(:,:,i_spin))
                end do
                call perfoff                
             end do
          end do
          call perfoff
          
       else
          call sync_ev(win_ev)
       end if
    enddo
    
    call finalize_access_ev(win_ev)
    
    !reduction over blacsdim to complete matrix multiplicaton combined with the sum over all k points
    count=n_basis*n_basis*n_spin
    call mpi_allreduce(MPI_IN_PLACE, exchange_matr_tmp, count, &
         MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, mpierr)

!    if(myid.eq.0) print*,'eemk',sum(abs(exchange_matr_tmp))
!stop
    fock_m = (0d0,0d0)
    do i_spin = 1, n_spin
       i_task = 0
       do i_k_point = 1, n_basis
          do i_q_point = 1, i_k_point
             i_task = i_task + 1
             fock_m(i_task,i_spin) = exchange_matr_tmp(i_q_point,i_k_point,i_spin)
          enddo
       enddo
    enddo
    
    deallocate(KS_eigenvector_q)
    deallocate(lvl_tricoeff_MO_r)
    if (myid_row.ne.myid_col) deallocate(lvl_tricoeff_MO_c_arr)
    deallocate(coulomb_times_tricoeff)
    deallocate(exchange_matr_tmp)
    deallocate(lb_atom, ub_atom)
    
    return
    
  end subroutine my_evaluate_exchange_matr_kspace_single_kpoint_p0
  !---------------------------------------------------------------------
  !******
  
  subroutine mult_with_tricoeffsum(q_weight, occ_number, coulomb_matr_blacs, lvl_tricoeff_MO_r, lvl_tricoeff_MO_c, &
       coulomb_times_tricoeff, exchange_matr_tmp)
    use dimensions
    use crpa_blacs
    
    implicit none
    real*8:: occ_number, q_weight
    complex*16, dimension(lbb_row:ubb_row,n_basis):: lvl_tricoeff_MO_r
    complex*16, dimension(lbb_col:ubb_col,n_basis):: lvl_tricoeff_MO_c
    complex*16, dimension(lbb_row:ubb_row,lbb_col:ubb_col):: coulomb_matr_blacs
    complex*16, dimension(n_bb_row,n_basis):: coulomb_times_tricoeff 
    complex*16, dimension(n_basis,n_basis):: exchange_matr_tmp
    
    !this is the local reduction, the matrix multiplication is only complete after mpi_allreduce in comm_blacs
!    call perfon('mu1')
    call zgemm('N', 'N', n_bb_row, n_basis, n_bb_col, (1.d0,0.d0), &
         coulomb_matr_blacs, n_bb_row, &
         lvl_tricoeff_MO_c(lbb_col,1), n_bb_col, (0.d0,0.d0), &
         coulomb_times_tricoeff, n_bb_row)
!    call perfoff
!    call perfon('mu2')
    call zgemm('C', 'N', n_basis, n_basis, n_bb_row, q_weight* &
         occ_number*dble(n_spin)/2.d0*(1.d0,0.d0), &
         lvl_tricoeff_MO_r(lbb_row,1), n_bb_row, &
         coulomb_times_tricoeff, n_bb_row, (1.d0,0.d0), &
         exchange_matr_tmp(:,:),n_basis)
!    call perfoff
  end subroutine mult_with_tricoeffsum
  
  
  subroutine comp_tricoeffsum(lvl_tricoeff_recip_k, lvl_tricoeff_recip_q, KS_eigenvector_q, lvl_tricoeff_MO, lbb, ubb)
    use dimensions
    use crpa_blacs
    use prodbas
    use basis
    use geometry
    
    implicit none
    
    integer:: lbb, ubb
    complex*16, dimension(n_basis):: KS_eigenvector_q
    complex*16, dimension(lbb:ubb,n_basis):: lvl_tricoeff_MO
    complex*16, dimension(lbb:ubb,max_n_basis_sp,n_basis), intent(in) :: lvl_tricoeff_recip_k
    complex*16, dimension(lbb:ubb,max_n_basis_sp), intent(in) :: lvl_tricoeff_recip_q  
    
    complex*16, dimension(n_basis):: evec_conjg
    integer:: i_atom, lb, ub, i_state_2, brange, i_prodbas, llbb, lubb, bboff, i_species, n_spbb, i_state, i_basis
    complex*16:: one

    evec_conjg=conjg(KS_eigenvector_q)
    one=1.

    lvl_tricoeff_MO=0.
    do i_atom = basbas_atom(lbb), basbas_atom(ubb)
       !basbas ranges for this atom
       bboff = atom2basbas_off(i_atom)
       i_species = species(i_atom)
       n_spbb = sp2n_basbas_sp(i_species)

       llbb=max(lbb,bboff+1)
       lubb=min(ubb,bboff+n_spbb)

       lb=lb_atom(i_atom)
       ub=ub_atom(i_atom)
       brange=ub_atom(i_atom)-lb_atom(i_atom)+1

       lvl_tricoeff_MO(llbb:lubb, lb:ub) = conjg(lvl_tricoeff_recip_q(llbb:lubb, 1:brange))
!call perfon('mu3')
       do i_basis=1,n_basis
          call zgemv('N',lubb-llbb+1, brange, one, &
          lvl_tricoeff_recip_k(llbb, 1, i_basis), ubb-lbb+1, &
          evec_conjg(lb), 1, one, &
          lvl_tricoeff_MO(llbb,i_basis), 1)
       end do
!call perfoff
    end do

  end subroutine comp_tricoeffsum

end module ex_mat_ksk_p0
