subroutine get_coulelement_lvl_v0(KS_eigen, KS_eigenvector,&
                                  KS_eigenvector_complex, &
                                  do_bcast)



!  PURPOSE
!  The subroutine gets the coulomb matrix in the product basis and transforms it
!  to the KS-basis for (k,k'.k-q,k'+q).
!
!  USES

      use synchronize_mpi
      use pbc_lists
      use constants
      use dimensions
      use runtime_choices
      use scalapack_wrapper, ONLY: mxld, mxcol
      use prodbas
      use physics, only: n_electrons
      use lvl_triples
      use tight_binding_auxmat
      use timing
      use synchronize_mpi_basic
      use mpi_tasks
      use sparse_tensor
      use hdf5_tools, only: HID_T, open_hdf5, close_hdf5, &
          open_coulelement_lvl_v0_1, out_coulelement_lvl_v0_1
      use calculate_mommat_base, only: get_state_minmax_k
      use localorb_io, only: localorb_info
      use geometry, only: species
      use basis
      implicit none 

  !  ARGUMENTS
  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
  complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: &
                                                        KS_eigenvector_complex
  logical :: do_bcast
  !  INPUTS
  !  o KS_eigen -- KS eigenvalues
  !  o  KS_eigenvector -- real array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent 
  !            calculation
  !  o  KS_eigenvector_complex -- complex array, the eigenvector of the 
  !                       single-particle (KS/HF) self-consistent calculation 
  !  o  do_bcast -- Switch between mpi_bcast and mpi_send/mpi_receive
  !  OUTPUTS
  !  hopefuly V_ak,bk',ck-q,dk'+q in our KS basis
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

  INTEGER(HID_T) :: file_id                  ! File identifier
  INTEGER(HID_T) :: plist_id                 ! Property list identifier
!     periodic Hartree Fock in kspace
  integer, allocatable :: k_minus_q_points_list(:,:)
  integer, allocatable :: kstrich_plus_q_points_list(:,:)
  complex*16, allocatable :: lvl_tricoeff_recip1(:,:,:,:)
  real*8, allocatable :: lvl_tricoeff_recip2(:,:,:)
  complex*16, allocatable :: coulomb_matr_recip(:,:,:)

  complex*16, allocatable :: lvl_tricoeff_recip_k(:,:,:) 
                                           ! LVL triple coefficients in k space
  complex*16, allocatable :: lvl_tricoeff_recip_k_minus_q(:,:,:) 
                                           ! LVL triple coefficients in k space
  complex*16, allocatable :: lvl_tricoeff_recip_kstrich(:,:,:) 
                                           ! LVL triple coefficients in k space
  complex*16, allocatable :: lvl_tricoeff_recip_kstrich_plus_q(:,:,:) 
                                          ! LVL triple coefficients in k space
  complex*16, allocatable :: KS_eigenvector_k_minus_q(:,:,:) 
  complex*16, allocatable :: KS_eigenvector_k(:,:,:) 
  complex*16, allocatable :: KS_eigenvector_kstrich(:,:,:) 
  complex*16, allocatable :: KS_eigenvector_kstrich_plus_q(:,:,:) 

  complex*16, allocatable :: lvl_tricoeff_sum1(:,:,:) 
                      ! sum of LVL triple coefficient at k and q points (atom 1)
  complex*16, allocatable :: lvl_tricoeff_sum2(:,:,:) 
                      ! sum of LVL triple coefficient at k and q points (atom 2)
  complex*16, allocatable :: lvl_tricoeff_l_tmp(:,:,:,:) 
  complex*16, allocatable :: lvl_tricoeff_r_tmp(:,:,:,:) 
  complex*16, allocatable :: lvl_tricoeff_l(:,:,:,:) 
               ! Full LVL triple coefficient, one basis index is transformed to
               ! (occupied) molecular orbital index 
  complex*16, allocatable :: lvl_tricoeff_r(:,:,:) 
               ! Full LVL triple coefficient, one basis index is transformed to
               ! (occupied) molecular orbital index 
  complex*16, allocatable :: tmp_MO(:,:,:)
  complex*16, allocatable :: coulomb_times_tricoeff(:,:) 
               ! Coulomb matrix multiplies the LVL triple coefficient
  complex*16, allocatable :: coulomb_matr_recip_tmp(:,:)

  complex*16, allocatable :: coul_matr(:,:,:,:)
  integer :: info, mpierr
  character(*), parameter :: func = 'calculate_Coulmat_p0'
  character*150 :: info_str

! counter 
  integer :: n_cells_task
  integer :: i_q_point, i_k_point, i_kstrich_point
  integer :: i_q_point_local, i_k_point_local, i_kstrich_point_local, &
             i_k_minus_q_point_local,  i_kstrich_plus_q_point_local
  integer :: i_spin, i_state, j_state
  integer :: id_root, id_rootk, id_rootkstrich, id_rootk_minus_q, &
             id_rootkstrich_plus_q
  integer :: i_task
  integer, allocatable :: k_points_at_task(:,:)
  integer :: i_req_ev1(0:n_tasks-1)
  integer :: i_req_ev2(0:n_tasks-1)
  integer :: i_req_lvl1(0:n_tasks-1)
  integer :: i_req_lvl2(0:n_tasks-1)
  integer :: i_atom_1, i_atom_2
  integer :: i_species_1, i_species_2
  integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
  integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
  integer :: i_1, i_2

  character(50) :: qstring
  integer::  n_state_min_in, n_state_max_in
  integer::  size_element, size_KS_vec
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  if(use_scalapack) call aims_stop('*** Periodic EXX in k-space must not be used if use_scalapack is in effect!')
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)	
  size_element = ((n_state_max_in-n_state_min_in+1)+1)*&
                              (n_state_max_in-n_state_min_in+1)/2
  size_KS_vec  = n_state_max_in-n_state_min_in+1
  allocate(coul_matr(size_KS_vec,size_KS_vec,size_KS_vec,size_KS_vec),stat=info) 
  call check_allocation(info, 'coul_matr', func)
  call open_hdf5('coulmat.h5', file_id, plist_id)
  call open_coulelement_lvl_v0_1("Coulomb_matrix", file_id, plist_id,& 
                                n_k_points, size_KS_vec)
  call mpi_barrier(mpi_comm_global,info)

  call get_timestamps(time_prodbas_total, clock_time_prodbas_total)
      ! --- Initialize

  ! First, clean HF quantities which depend on pairs
  ! (might change with relaxation)
  call cleanup_basbas()

  ! Construct product basis & basis pairs
  call initialize_prodbas()

  ! estimate the largest number of real-space unit cells locally
  if(mod(n_cells, n_tasks).eq.0) then
      n_cells_task = n_cells/n_tasks
  else
      n_cells_task = n_cells/n_tasks+1
  endif

  call initialize_lvl_triples(OVLP_TYPE_COULOMB)
  call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)

  allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_task)&
                               ,stat=info) 
  call check_allocation(info, 'lvl_tricoeff_recip1', func)
  allocate(lvl_tricoeff_recip2(max_n_basbas_sp,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_recip2', func)
  call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,&
                               lvl_tricoeff_recip2)

  call cleanup_lvl_triples()
  call deallocate_tb_auxmat()

!!call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
  call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
!!call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
!!call initialize_periodic_tb_auxmat(1, 1.d0)
  allocate(k_minus_q_points_list(n_k_points,n_k_points),stat=info) 
  call check_allocation(info, 'k_minus_q_points_list', func)
  allocate(kstrich_plus_q_points_list(n_k_points,n_k_points),stat=info) 
  call check_allocation(info, 'kstrich_plus_q_points_list', func)
  call determine_k_minus_q_list(k_minus_q_points_list,&
                                           kstrich_plus_q_points_list)

  allocate(coulomb_matr_recip(n_basbas,n_basbas,n_k_points_task),stat=info) 
  call check_allocation(info, 'coulomb_matr_recip', func)
  call get_coulomb_matr_recip(coulomb_matr_recip,1)

  call deallocate_tb_auxmat()

  call get_times(time_prodbas_total, clock_time_prodbas_total)

  write(info_str,'(2X,A)') "Constructing the Coulomb matrix elements in (k,k',q) space ..."
  call localorb_info(info_str)

  ! k_points_at_task: which task contains which k-point
  allocate(k_points_at_task(0:n_tasks-1,(n_k_points-1)/n_tasks+1))
  k_points_at_task(:,:) = 0
  do i_k_point = 1, n_k_points, 1
    i_k_point_local = (i_k_point-1)/n_tasks + 1 
    k_points_at_task(mod(i_k_point,n_tasks),i_k_point_local) = i_k_point
  enddo

  i_req_ev1(:) = MPI_REQUEST_NULL
  i_req_ev2(:) = MPI_REQUEST_NULL
  i_req_lvl1(:) = MPI_REQUEST_NULL
  i_req_lvl2(:) = MPI_REQUEST_NULL

  do i_q_point = 1, n_k_points, 1
     i_q_point_local =(i_q_point-1)/n_tasks + 1 
!     i_q_point = k_points_at_task(myid, i_q_point_local)
     write( qstring, '(i10)' )  i_q_point
     write(info_str,'(2X,A,A)') "Performing calcultions for q point: ",qstring
     call localorb_info(info_str)
     id_root = mod(i_q_point,n_tasks)
     allocate(coulomb_matr_recip_tmp(n_basbas,n_basbas),stat=info) 
     call check_allocation(info, 'coulomb_matr_recip_tmp', func)
     if(myid == id_root) then
	! Tasks i_task needs i_q_point from me
	coulomb_matr_recip_tmp(1:n_basbas,1:n_basbas)=coulomb_matr_recip(&
                                        1:n_basbas,1:n_basbas,i_q_point_local)
     endif
     call mpi_bcast(coulomb_matr_recip_tmp,n_basbas*n_basbas, MPI_COMPLEX16, &
                                             id_root, mpi_comm_global, mpierr)
!     do i_k_point_local = 1, (n_k_points-1)/n_tasks + 1
!     do i_k_point = 1, n_k_points, 1
     do i_k_point = 1, n_k_points, 1

	id_rootk = mod(i_k_point,n_tasks)
	id_rootk_minus_q = mod(k_minus_q_points_list(i_k_point,i_q_point)&
                                                                     ,n_tasks)
	i_k_point_local = (i_k_point-1)/n_tasks + 1 
	i_k_minus_q_point_local = (k_minus_q_points_list(i_k_point,&
                                                      i_q_point)-1)/n_tasks + 1 
	allocate(lvl_tricoeff_recip_k(max_n_basbas_sp,n_basis,n_basis)&
                 ,stat=info) 
	call check_allocation(info, 'lvl_tricoeff_recip_k', func)
	allocate(lvl_tricoeff_recip_k_minus_q(max_n_basbas_sp,n_basis,n_basis)&
                 ,stat=info) 
	call check_allocation(info, 'lvl_tricoeff_recip_k_minus_q', func)
	allocate(KS_eigenvector_k(n_basis,size_KS_vec,n_spin),stat=info) 
	call check_allocation(info, 'KS_eigenvector_k', func)
	allocate(KS_eigenvector_k_minus_q(n_basis,size_KS_vec,n_spin),stat=info) 
	call check_allocation(info, 'KS_eigenvector_k_minus_q', func)
	if(myid .eq. id_rootk) then
	    lvl_tricoeff_recip_k(:,:,:)=lvl_tricoeff_recip1(:,:,:,&
                                                            i_k_point_local)
	if(real_eigenvectors) then
	  KS_eigenvector_k(:,:,:)=KS_eigenvector(:,&
                    n_state_min_in:n_state_max_in,:,i_k_point_local)*(1.d0,0.d0)
	else
	  KS_eigenvector_k(:,:,:)=KS_eigenvector_complex(:,&
                                n_state_min_in:n_state_max_in,:,i_k_point_local)
	endif
	endif   
	if(myid .eq. id_rootk_minus_q) then
	    lvl_tricoeff_recip_k_minus_q(:,:,:)=lvl_tricoeff_recip1(:,:,:,&
                                                      i_k_minus_q_point_local)
	if(real_eigenvectors) then
	  KS_eigenvector_k_minus_q(:,:,:)=KS_eigenvector(:,&
                                            n_state_min_in:n_state_max_in,:,&
                                            i_k_minus_q_point_local)*(1.d0,0.d0)
	else
	  KS_eigenvector_k_minus_q(:,:,:)=KS_eigenvector_complex(:,&
                       n_state_min_in:n_state_max_in,:, i_k_minus_q_point_local)
	endif
	endif 
        if (do_bcast)then
	  call mpi_bcast(lvl_tricoeff_recip_k,max_n_basbas_sp*n_basis*n_basis,&
			    MPI_COMPLEX16, id_rootk, mpi_comm_global, mpierr)
	  call mpi_bcast(lvl_tricoeff_recip_k_minus_q,&
                     max_n_basbas_sp*n_basis*n_basis, &
		     MPI_COMPLEX16, id_rootk_minus_q, mpi_comm_global, mpierr)
	  call mpi_bcast(KS_eigenvector_k,n_basis*size_KS_vec*n_spin, &
			    MPI_COMPLEX16, id_rootk, mpi_comm_global, mpierr)
	  call mpi_bcast(KS_eigenvector_k_minus_q,n_basis*size_KS_vec*n_spin, &
		    MPI_COMPLEX16, id_rootk_minus_q, mpi_comm_global, mpierr)
        else
	    if (id_rootk .ne. id_rootk_minus_q) then
		if(myid.eq.id_rootk) then
                  !call mpi_isend(lvl_tricoeff_recip_k,max_n_basbas_sp*n_basis*&
                  !               n_basis, MPI_COMPLEX16, id_rootk_minus_q,&
                  !               111, mpi_comm_global, &
                  !               i_req_lvl1(id_rootk_minus_q), mpierr)
                  !call mpi_isend(KS_eigenvector_k,n_basis*n_spin, &
                  !               MPI_COMPLEX16, id_rootk_minus_q, 112, &
                  !               mpi_comm_global, i_req_ev1(id_rootk_minus_q),&
                  !               mpierr)   
                  call mpi_recv(lvl_tricoeff_recip_k_minus_q, &
                       max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
                       id_rootk_minus_q, 113, mpi_comm_global, &
                       MPI_STATUS_IGNORE, mpierr)
	          call mpi_wait(i_req_lvl1(id_rootk), MPI_STATUSES_IGNORE, &
                       mpierr) 
                  call mpi_recv(KS_eigenvector_k_minus_q,&
                       n_basis*size_KS_vec*n_spin, &
                       MPI_COMPLEX16, id_rootk_minus_q, 114, mpi_comm_global, &
                       MPI_STATUS_IGNORE, mpierr)
	          call mpi_wait(i_req_ev1(id_rootk), MPI_STATUSES_IGNORE, &
                       mpierr) 
		endif
		if(myid.eq.id_rootk_minus_q) then
                  call mpi_isend(lvl_tricoeff_recip_k_minus_q, &
                       max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
                       id_rootk, 113, mpi_comm_global, i_req_lvl1(id_rootk), &
                       mpierr) 
                  call mpi_isend(KS_eigenvector_k_minus_q,&
                       n_basis*size_KS_vec*n_spin, &
                       MPI_COMPLEX16, id_rootk, 114, mpi_comm_global, &
                       i_req_ev1(id_rootk), mpierr)
                  !call mpi_recv(lvl_tricoeff_recip_k, &
                  !     max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
                  !     id_rootk, 111, mpi_comm_global, MPI_STATUS_IGNORE, &
                  !     mpierr)
	          !call mpi_wait(i_req_lvl1(id_rootk_minus_q), &
                  !     MPI_STATUSES_IGNORE, mpierr) 
                  !call mpi_recv(KS_eigenvector_k,n_basis*n_spin,&
                  !     MPI_COMPLEX16, id_rootk, 112, mpi_comm_global, &
                  !     MPI_STATUS_IGNORE, mpierr)
	          !call mpi_wait(i_req_ev1(id_rootk_minus_q), &
                  !     MPI_STATUSES_IGNORE, mpierr)                  
		endif  
	    endif
        endif
	KS_eigenvector_k = conjg(KS_eigenvector_k)
	allocate(lvl_tricoeff_l(n_basbas,size_KS_vec,size_KS_vec,n_spin),&
                                stat=info) 
	call check_allocation(info, 'lvl_tricoeff_l', func)
        if(myid.eq.id_rootk) then
	  allocate(lvl_tricoeff_l_tmp(n_basbas,n_basis,size_KS_vec,n_spin),&
                                      stat=info) 
	  call check_allocation(info, 'lvl_tricoeff_l_tmp', func)

	  lvl_tricoeff_l_tmp = (0.d0,0.d0)

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
		  conjg(lvl_tricoeff_recip_k(1:n_spbb_1,basis_off_1+i_1, i_2))+&
		  lvl_tricoeff_recip2(1:n_spbb_1, i_2, basis_off_1+i_1)
		enddo
	    enddo

	    ! Multiply lvl_tricoeff_sum1 with KS_eigenvector_q
	    ! Please note that the multipication uses the first 2 dimensions of
	    ! lvl_tricoeff_sum1/tmp_MO as 1 logical dimension, thus these 2 
	    ! arrays must be allocated to the exact size!

	    allocate(tmp_MO(n_spbb_1,n_sp_basis_1,size_KS_vec),stat=info)
	    call check_allocation(info, 'tmp_MO', func)

	    do i_spin = 1, n_spin
		call zgemm('N', 'N', n_spbb_1*n_sp_basis_1,size_KS_vec,n_basis,&
		       (1.d0,0.d0),lvl_tricoeff_sum1,n_spbb_1*n_sp_basis_1, &
		       KS_eigenvector_k(1,1,i_spin), n_basis, (0.d0,0.d0), &
		       tmp_MO, n_spbb_1*n_sp_basis_1)
		lvl_tricoeff_l_tmp(bboff_1+1:bboff_1+n_spbb_1,basis_off_1+& 
				  1:basis_off_1+n_sp_basis_1,:,i_spin) = tmp_MO
	    enddo
	    deallocate(tmp_MO)
	    deallocate(lvl_tricoeff_sum1)
	  enddo
	  do i_spin = 1, n_spin
	    do i_state = 1, size_KS_vec, 1
	    enddo
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
		 lvl_tricoeff_recip_k_minus_q(1:n_spbb_2,basis_off_2+i_2, i_1)+&
		 lvl_tricoeff_recip2(1:n_spbb_2, i_1, basis_off_2+i_2)
		enddo
	    enddo

	    ! Multiply lvl_tricoeff_sum2 with KS_eigenvector_q
	    ! See remark above about matrix multiplication

	    allocate(tmp_MO(n_spbb_2,n_basis,size_KS_vec))
	    call check_allocation(info, 'tmp_MO', func)

	    do i_spin = 1, n_spin
		call zgemm('N', 'N', n_spbb_2*n_basis,size_KS_vec,n_sp_basis_2,&
		       (1.d0,0.d0), lvl_tricoeff_sum2, n_spbb_2*n_basis, &
		       KS_eigenvector_k(basis_off_2+1,1,i_spin), n_basis, &
		       (0.d0,0.d0), tmp_MO, n_spbb_2*n_basis)
		lvl_tricoeff_l_tmp(bboff_2+1:bboff_2+n_spbb_2,:,:,i_spin) = &
		lvl_tricoeff_l_tmp(bboff_2+1:bboff_2+n_spbb_2,:,:,i_spin) + &
                tmp_MO
	    enddo
	    deallocate(tmp_MO)
	    deallocate(lvl_tricoeff_sum2)
	  enddo
	  lvl_tricoeff_l=0
	  do i_spin = 1, n_spin
	    do i_state = 1, size_KS_vec, 1
	      call zgemm('N', 'N', n_basbas, size_KS_vec, n_basis, (1.d0,0.d0),&
	       lvl_tricoeff_l_tmp(1:n_basbas,1:n_basis,i_state,i_spin), &
               n_basbas, KS_eigenvector_k_minus_q(1,1,i_spin), n_basis, &
               (0.d0,0.d0), lvl_tricoeff_l(1:n_basbas,i_state, &
               1:size_KS_vec,i_spin), n_basbas)
	    enddo
	  enddo
  !        call mpi_bcast(lvl_tricoeff_l,n_basbas*n_spin, &
  !			MPI_COMPLEX16, id_rootk, mpi_comm_global, mpierr)
	  deallocate(lvl_tricoeff_l_tmp)
        endif
	if (allocated(lvl_tricoeff_recip_k)) deallocate(lvl_tricoeff_recip_k)
	if (allocated(lvl_tricoeff_recip_k_minus_q)) &
                                     deallocate(lvl_tricoeff_recip_k_minus_q)
	if (allocated(KS_eigenvector_k)) deallocate(KS_eigenvector_k)
	if (allocated(KS_eigenvector_k_minus_q)) &
                                         deallocate(KS_eigenvector_k_minus_q)
        do i_kstrich_point = 1,n_k_points, 1
!        do i_kstrich_point_local = 1,(n_k_points-1)/n_tasks + 1, 1

	    ! Send our coulomb_matr_recip to the task(s) which actually needs it
	    ! Normally, only one send should be necessary

	    id_rootkstrich = mod(i_kstrich_point,n_tasks)
	    id_rootkstrich_plus_q = mod(kstrich_plus_q_points_list(&
                                            i_kstrich_point,i_q_point),n_tasks)
	    i_kstrich_point_local = (i_kstrich_point-1)/n_tasks + 1 
	    i_kstrich_plus_q_point_local = (kstrich_plus_q_points_list(&
                                      i_kstrich_point,i_q_point)-1)/n_tasks + 1 

	    allocate(lvl_tricoeff_recip_kstrich(max_n_basbas_sp,n_basis,&
                                                            n_basis),stat=info) 
	    call check_allocation(info, 'lvl_tricoeff_recip_kstrich', func)
	    allocate(lvl_tricoeff_recip_kstrich_plus_q(max_n_basbas_sp,&
                                                    n_basis,n_basis),stat=info) 
	    call check_allocation(info, 'lvl_tricoeff_recip_kstrich_plus_q',&
                                  func)
	    allocate(KS_eigenvector_kstrich(n_basis,size_KS_vec,n_spin),&
                     stat=info) 
	    call check_allocation(info, 'KS_eigenvector_kstrich', func)
	    allocate(KS_eigenvector_kstrich_plus_q(n_basis,size_KS_vec,n_spin),&
                     stat=info) 
	    call check_allocation(info, 'KS_eigenvector_kstrich_plus_q', func)

	    if(myid .eq. id_rootkstrich) then
		lvl_tricoeff_recip_kstrich(:,:,:)=&
                                lvl_tricoeff_recip1(:,:,:,i_kstrich_point_local)
	    if(real_eigenvectors) then
	      KS_eigenvector_kstrich(:,:,:)=KS_eigenvector(:,&
              n_state_min_in:n_state_max_in,:,i_kstrich_point_local)*(1.d0,0.d0)
	    else
	      KS_eigenvector_kstrich(:,:,:)=KS_eigenvector_complex(:,&
                          n_state_min_in:n_state_max_in,:,i_kstrich_point_local)
	    endif
	    endif   
	    if(myid .eq. id_rootkstrich_plus_q) then
		lvl_tricoeff_recip_kstrich_plus_q(:,:,:)=&
                         lvl_tricoeff_recip1(:,:,:,i_kstrich_plus_q_point_local)
	    if(real_eigenvectors) then
	      KS_eigenvector_kstrich_plus_q(:,1:size_KS_vec,:)=&
                      KS_eigenvector(:,n_state_min_in:n_state_max_in,:,&
                                     i_kstrich_plus_q_point_local)*(1.d0,0.d0)
	    else
	      KS_eigenvector_kstrich_plus_q(1:n_basis,1:size_KS_vec,1:n_spin)=&
                      KS_eigenvector_complex(1:n_basis,&
                      n_state_min_in:n_state_max_in,1:n_spin,&
                      i_kstrich_plus_q_point_local)
	    endif
	    endif 
	    if (do_bcast)then
		call mpi_bcast(lvl_tricoeff_recip_kstrich,&
                               max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
                               id_rootkstrich, mpi_comm_global, mpierr)
		call mpi_bcast(lvl_tricoeff_recip_kstrich_plus_q,&
                               max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16,&
                               id_rootkstrich_plus_q, mpi_comm_global, mpierr)

		call mpi_bcast(KS_eigenvector_kstrich,&
                        n_basis*size_KS_vec*n_spin, &
			MPI_COMPLEX16, id_rootkstrich, mpi_comm_global, mpierr)
		call mpi_bcast(KS_eigenvector_kstrich_plus_q,&
                        n_basis*size_KS_vec*n_spin, MPI_COMPLEX16, &
                        id_rootkstrich_plus_q, mpi_comm_global, mpierr)
	    else
		if (id_rootkstrich .ne. id_rootkstrich_plus_q) then
		    if (myid.eq.id_rootkstrich) then
		      !call mpi_isend(lvl_tricoeff_recip_kstrich,&
                      !     max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
		      !     id_rootkstrich_plus_q, 115, mpi_comm_global, &
                      !     i_req_lvl2(id_rootkstrich_plus_q), mpierr)
                      !call mpi_isend(KS_eigenvector_kstrich,n_basis*n_spin, &
                      !     MPI_COMPLEX16, id_rootkstrich_plus_q, 116, &
                      !     mpi_comm_global, i_req_ev2(id_rootkstrich_plus_q), &
                      !     mpierr)
		      call mpi_recv(lvl_tricoeff_recip_kstrich_plus_q, &
                           max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
			   id_rootkstrich_plus_q, 117, mpi_comm_global, &
                           MPI_STATUS_IGNORE, mpierr)
	              call mpi_wait(i_req_lvl2(id_rootkstrich), &
                           MPI_STATUSES_IGNORE, mpierr) 
		      call mpi_recv(KS_eigenvector_kstrich_plus_q,&
                           n_basis*size_KS_vec*n_spin, MPI_COMPLEX16, &
			   id_rootkstrich_plus_q, 118, mpi_comm_global, &
                           MPI_STATUS_IGNORE, mpierr)
	              call mpi_wait(i_req_ev2(id_rootkstrich), &
                           MPI_STATUSES_IGNORE, mpierr) 
		    elseif(myid.eq.id_rootkstrich_plus_q) then
		      call mpi_isend(lvl_tricoeff_recip_kstrich_plus_q, &
                           max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
			   id_rootkstrich, 117, mpi_comm_global, &
                           i_req_lvl2(id_rootkstrich), mpierr)
                      call mpi_isend(KS_eigenvector_kstrich_plus_q,&
                           n_basis*size_KS_vec*n_spin, MPI_COMPLEX16, &
                           id_rootkstrich, 118, mpi_comm_global, &
                           i_req_ev2(id_rootkstrich), mpierr)
		      !call mpi_recv(lvl_tricoeff_recip_kstrich, &
                      !     max_n_basbas_sp*n_basis*n_basis, MPI_COMPLEX16, &
		      !     id_rootkstrich, 115, mpi_comm_global, &
                      !     MPI_STATUS_IGNORE, mpierr)
	              !call mpi_wait(i_req_lvl2(id_rootkstrich_plus_q), &
                      !     MPI_STATUSES_IGNORE, mpierr) 
		      !call mpi_recv(KS_eigenvector_kstrich,n_basis*n_spin, &
                      !     MPI_COMPLEX16, id_rootkstrich, 116, &
                      !     mpi_comm_global, MPI_STATUS_IGNORE, mpierr)
	              !call mpi_wait(i_req_ev2(id_rootkstrich_plus_q), &
                      !     MPI_STATUSES_IGNORE, mpierr) 
		    endif
		endif
	    endif
            if (myid.eq.id_rootkstrich) then
	      KS_eigenvector_kstrich = conjg(KS_eigenvector_kstrich)

	      allocate(lvl_tricoeff_r_tmp(n_basbas,n_basis,size_KS_vec,n_spin),&
                       stat=info) 
	      call check_allocation(info, 'lvl_tricoeff_r_tmp', func)

	      lvl_tricoeff_r_tmp = 0

	      do i_atom_1 = 1, n_atoms, 1
		i_species_1 = species(i_atom_1)
		basis_off_1 = atom2basis_off(i_atom_1)
		n_sp_basis_1 = sp2n_basis_sp(i_species_1)
		bboff_1 = atom2basbas_off(i_atom_1)
		n_spbb_1 = sp2n_basbas_sp(i_species_1)

		allocate(lvl_tricoeff_sum1(n_spbb_1,n_sp_basis_1,n_basis),&
                         stat=info)
		call check_allocation(info, 'lvl_tricoeff_sum1', func)

		do i_1 = 1, n_sp_basis_1
		    do i_2 = 1, n_basis
		      lvl_tricoeff_sum1(1:n_spbb_1,i_1,i_2) = &
			  conjg(lvl_tricoeff_recip_kstrich(1:n_spbb_1, &
                          basis_off_1+i_1, i_2)) + &
			  lvl_tricoeff_recip2(1:n_spbb_1, i_2, basis_off_1+i_1)
		    enddo
		enddo

		! Multiply lvl_tricoeff_sum1 with KS_eigenvector_q
		! Please note that the multipication uses the first 2 dimensions
		! of lvl_tricoeff_sum1/tmp_MO as 1 logical dimension, thus these 
		! 2 arrays must be allocated to the exact size!

		allocate(tmp_MO(n_spbb_1,n_sp_basis_1,size_KS_vec),stat=info)
		call check_allocation(info, 'tmp_MO', func)

		do i_spin = 1, n_spin
		    call zgemm('N', 'N', n_spbb_1*n_sp_basis_1, size_KS_vec, &
                           n_basis, (1.d0,0.d0), &
			   lvl_tricoeff_sum1, n_spbb_1*n_sp_basis_1, &
			   KS_eigenvector_kstrich(1,1,i_spin), n_basis, &
                           (0.d0,0.d0), tmp_MO, n_spbb_1*n_sp_basis_1)
		    lvl_tricoeff_r_tmp(bboff_1+1:bboff_1+n_spbb_1,&
                                       basis_off_1+1:basis_off_1+n_sp_basis_1,:&
                                       ,i_spin) = tmp_MO
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

		allocate(lvl_tricoeff_sum2(n_spbb_2,n_basis,n_sp_basis_2),&
                         stat=info)
		call check_allocation(info, 'lvl_tricoeff_sum2', func)

		do i_1 = 1, n_basis
		    do i_2 = 1, n_sp_basis_2
		      lvl_tricoeff_sum2(1:n_spbb_2,i_1,i_2) = &
			  lvl_tricoeff_recip_kstrich_plus_q(1:n_spbb_2, &
                          basis_off_2+i_2, i_1) + &
			  lvl_tricoeff_recip2(1:n_spbb_2, i_1, basis_off_2+i_2)
		    enddo
		enddo

		! Multiply lvl_tricoeff_sum2 with KS_eigenvector_q
		! See remark above about matrix multiplication

		allocate(tmp_MO(n_spbb_2,n_basis,size_KS_vec))
		call check_allocation(info, 'tmp_MO', func)

		do i_spin = 1, n_spin
		  call zgemm('N', 'N', n_spbb_2*n_basis, size_KS_vec, &
                         n_sp_basis_2, (1.d0,0.d0), &
			 lvl_tricoeff_sum2, n_spbb_2*n_basis, &
			 KS_eigenvector_kstrich(basis_off_2+1,1,i_spin), &
                         n_basis, (0.d0,0.d0), &
			 tmp_MO, n_spbb_2*n_basis)
  		  lvl_tricoeff_r_tmp(bboff_2+1:bboff_2+n_spbb_2,:,:,i_spin) = &
  		  lvl_tricoeff_r_tmp(bboff_2+1:bboff_2+n_spbb_2,:,:,i_spin) + &
                                                                          tmp_MO
		enddo
		deallocate(tmp_MO)
		deallocate(lvl_tricoeff_sum2)

	      enddo
            endif
	    allocate(coulomb_times_tricoeff(n_basbas,size_KS_vec),stat=info) 
	    call check_allocation(info, 'coulomb_times_tricoeff', func)
	    allocate(lvl_tricoeff_r(n_basbas,size_KS_vec,n_spin),stat=info) 
	    call check_allocation(info, 'lvl_tricoeff_r', func)
	    lvl_tricoeff_r=0
	    if (id_rootk .ne. id_rootkstrich) then
	      if(myid.eq.id_rootk) then
                  call mpi_isend(lvl_tricoeff_l,&
                       max_n_basbas_sp*size_KS_vec*size_KS_vec*n_spin, &
                       MPI_COMPLEX16, id_rootkstrich, 111, mpi_comm_global, &
                       i_req_lvl2(id_rootkstrich), mpierr) 
              endif
	      if(myid.eq.id_rootkstrich) then
                call mpi_recv(lvl_tricoeff_l,&
                     max_n_basbas_sp*size_KS_vec*size_KS_vec*n_spin, &
                     MPI_COMPLEX16, id_rootk, 111, mpi_comm_global, &
                     MPI_STATUS_IGNORE, mpierr)
	        call mpi_wait(i_req_lvl2(id_rootkstrich), MPI_STATUSES_IGNORE, &
                     mpierr) 
              endif
	    endif
            if(myid.eq.id_rootkstrich) then

	      do i_spin = 1, n_spin, 1
                 do j_state = 1, size_KS_vec, 1
		    call zgemm('N', 'N', n_basbas, size_KS_vec, n_basis, &
                      (1.d0,0.d0), lvl_tricoeff_r_tmp(1:n_basbas,1:n_basis,&
                      j_state,i_spin), n_basbas, &
                      KS_eigenvector_kstrich_plus_q(1,1,i_spin), n_basis, &
                      (0.d0,0.d0), lvl_tricoeff_r(1:n_basbas,1:size_KS_vec,&
                      i_spin), n_basbas)
		    call zgemm('N', 'N', n_basbas, size_KS_vec, n_basbas, &
                      (1.d0,0.d0), coulomb_matr_recip_tmp(1:n_basbas,&
                      1:n_basbas), n_basbas, lvl_tricoeff_r(1:n_basbas,&
                      1:size_KS_vec,i_spin), n_basbas, (0.d0,0.d0), &
		      coulomb_times_tricoeff(1:n_basbas,1:size_KS_vec), &
                      n_basbas)
                    do i_state = 1, size_KS_vec, 1
		      call zgemm('C', 'N', size_KS_vec, size_KS_vec, n_basbas, &
                        (1.d0,0.d0), lvl_tricoeff_l(1:n_basbas,i_state,&
                        1:size_KS_vec,i_spin), n_basbas, &
                        coulomb_times_tricoeff(1:n_basbas,1:size_KS_vec), &
                        n_basbas, (0.d0,0.d0), coul_matr(i_state,&
                        1:size_KS_vec,j_state,1:size_KS_vec),size_KS_vec)
                    enddo
                  enddo
	      enddo
	      call out_coulelement_lvl_v0_1(coul_matr,i_q_point,i_k_point,&
                   i_kstrich_point,"Coulomb_matrix", file_id, plist_id,&
                   size_KS_vec)
            endif
            if (allocated(coulomb_times_tricoeff)) &
                deallocate(coulomb_times_tricoeff)
            if (allocated(lvl_tricoeff_r)) deallocate(lvl_tricoeff_r)

	    if (allocated(lvl_tricoeff_r_tmp)) deallocate(lvl_tricoeff_r_tmp)
            if (allocated(lvl_tricoeff_recip_kstrich)) &
                deallocate(lvl_tricoeff_recip_kstrich) 
            if (allocated(lvl_tricoeff_recip_kstrich_plus_q)) &
                deallocate(lvl_tricoeff_recip_kstrich_plus_q) 
            if (allocated(KS_eigenvector_kstrich)) &
                deallocate(KS_eigenvector_kstrich) 
            if (allocated(KS_eigenvector_kstrich_plus_q)) &
                deallocate(KS_eigenvector_kstrich_plus_q) 
        enddo   
        if (allocated(lvl_tricoeff_l)) deallocate(lvl_tricoeff_l) 
     enddo
     if (allocated(coulomb_matr_recip_tmp)) deallocate(coulomb_matr_recip_tmp) 

  enddo
  if (allocated(lvl_tricoeff_recip1)) deallocate(lvl_tricoeff_recip1) 
  if (allocated(lvl_tricoeff_recip2)) deallocate(lvl_tricoeff_recip2)
  if (allocated(coulomb_matr_recip)) deallocate(coulomb_matr_recip)
  if (allocated(k_minus_q_points_list)) deallocate(k_minus_q_points_list)
  if (allocated(kstrich_plus_q_points_list)) &
      deallocate(kstrich_plus_q_points_list)
  if (allocated(coul_matr)) deallocate(coul_matr)
  call mpi_barrier(mpi_comm_global,info)
  call close_hdf5(file_id,plist_id)

endsubroutine get_coulelement_lvl_v0
