subroutine get_coulelement_ovl(partition_tab_std_in, basis_l_max_in,&
                         KS_eigen, KS_eigenvector_in,KS_eigenvector_complex_in&
                         )
!  PURPOSE
!
!  Wrapper function for calculating and outputting of the fourier components 
!  of the coulomb matrix
!  USES
      use synchronize_mpi
      use pbc_lists
      use constants
      use dimensions
      use runtime_choices
      use timing
      use synchronize_mpi_basic
      use mpi_tasks
      use calculate_mommat_base, only: get_state_minmax_k
      use calculate_coulmat_ovl
      use hdf5_tools, only: HID_T, open_hdf5, open_coulelement_kart, &
          out_k_points_kart, out_coulelement_kart, close_hdf5
      use localorb_io, only: use_unit
      implicit none 
!  ARGUMENTS
      real*8, target, dimension(n_full_points), INTENT(IN) :: &
                                                          partition_tab_std_in
      integer, INTENT(IN) :: basis_l_max_in (n_species)
      real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
      real*8, dimension(n_basis,n_states,n_spin,n_k_points_task), &
                                         INTENT(IN) :: KS_eigenvector_in
      complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task), &
                                         INTENT(IN) :: KS_eigenvector_complex_in
!  INPUTS
!   o partition_tab -- Partition tab
!   o l_shell_max
!   o KS_eigen -- KS eigenvalues
!   o KS_vec/KS_vec_complex -- KS coefficients
!  OUTPUT
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

      complex*16, allocatable :: coulelement_tmp(:)
      complex*16, allocatable :: KS_eigenvector_complex_k(:,:,:) 
      complex*16, allocatable :: KS_eigenvector_complex_k_strich(:,:,:) 
      real*8, allocatable :: KS_eigenvector_k(:,:,:) 
      real*8, allocatable :: KS_eigenvector_k_strich(:,:,:) 
      real*8, allocatable :: q_point_list(:,:) 
      integer :: n_q_point

      INTEGER(HID_T) :: file_id                  ! File identifier
      INTEGER(HID_T) :: plist_id                 ! Property list identifier

  !  counters
      integer::  i_q_point, i_k_point, i_k_strich_point
      integer::  id_root_q, id_rootk, id_rootk_strich
      integer::  i_k_point_local, i_k_strich_point_local, i_q_point_local
      integer::  i_state, j_state, num 
      integer::  n_state_min_in, n_state_max_in
      integer::  size_element, size_KS_vec
      integer :: info, mpierr
      character(*), parameter :: func = 'get_coulelement_ovl'
      character*150 :: info_str

      integer :: i_req(0:n_tasks-1)

      call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)	
      if(read_q_points)then
          call check_q_point_list(n_q_point)
	  allocate(q_point_list(n_q_point,3),stat=info) 
	  call check_allocation(info, 'q_point_list', func)
	  call get_q_point_list(q_point_list, n_q_point)
      else
	  allocate(q_point_list(n_k_points,3),stat=info) 
	  call check_allocation(info, 'q_point_list', func)
          q_point_list=k_point_list
          n_q_point=n_k_points
      endif
      call open_hdf5('coulmat.h5', file_id, plist_id)
      call open_coulelement_kart("Coulomb_matrix", file_id, plist_id,&
                            n_state_min_in, n_state_max_in, n_q_point)
      call out_k_points_kart(file_id,plist_id)
      call mpi_barrier(mpi_comm_global,info)
      size_element = ((n_state_max_in-n_state_min_in+1)+1)*&
                              (n_state_max_in-n_state_min_in+1)/2
      size_KS_vec  = n_state_max_in-n_state_min_in+1
      if(real_eigenvectors) then
	  allocate(KS_eigenvector_k(n_basis,size_KS_vec,n_spin),stat=info) 
	  call check_allocation(info, 'KS_eigenvector_k', func)
	  allocate(KS_eigenvector_k_strich(n_basis,size_KS_vec,n_spin),&
                                                                 stat=info) 
	  call check_allocation(info, 'KS_eigenvector_k_strich', func)
      else
	  allocate(KS_eigenvector_complex_k(n_basis,size_KS_vec,n_spin),&
                                                                 stat=info) 
	  call check_allocation(info, 'KS_eigenvector_complex_k', func)
	  allocate(KS_eigenvector_complex_k_strich(n_basis,size_KS_vec,n_spin)&
                                                               , stat=info) 
	  call check_allocation(info, 'KS_eigenvector_complex_k_strich', func)
      endif
      allocate(coulelement_tmp(size_element),stat=info)
      do i_q_point = 1, n_q_point, 1
          if (myid==0) write(use_unit,*) 'q-point: ', i_q_point
          call allocate_coulmat_w
	  coulmat_full_w_k=(0.0, 0.0)
	  call calculate_coulmat ( partition_tab_std_in, basis_l_max_in, &
                                   coulmat_full_w_k,q_point_list(i_q_point,1:3))  
	  call allocate_coulelement(n_state_min_in, n_state_max_in)
	  do i_k_point = 1, n_k_points, 1
	    id_rootk = mod(i_k_point,n_tasks)
	    i_k_point_local = (i_k_point-1)/n_tasks + 1 
	    coulelement_k = (0.0,0.0)
            i_req(:) = MPI_REQUEST_NULL
            coulelement_tmp(:) = (0.0, 0.0)
            do i_k_strich_point = 1, n_k_points, 1
		id_rootk_strich = mod(i_k_strich_point,n_tasks)
		i_k_strich_point_local = (i_k_strich_point-1)/n_tasks + 1
	
		if(myid .eq. id_rootk.and. myid <= n_k_points) then
		  if(real_eigenvectors) then
		    KS_eigenvector_k(:,:,:)=&
                                KS_eigenvector_complex_in(:,&
                                n_state_min_in:n_state_max_in,:,i_k_point_local)
		    call mpi_bcast(KS_eigenvector_k,&
                                   n_basis*size_KS_vec*n_spin,MPI_COMPLEX16, &
                                   id_rootk, mpi_comm_global, mpierr)
		  else
		    KS_eigenvector_complex_k(:,:,:)=&
                                KS_eigenvector_complex_in(:,&
                                n_state_min_in:n_state_max_in,:,i_k_point_local)
		    call mpi_isend(KS_eigenvector_complex_k,&
                                   n_basis*size_KS_vec*n_spin,MPI_COMPLEX16, &
                                   id_rootk_strich, 111, mpi_comm_global, &
                                   i_req(id_rootk_strich), mpierr)
		  endif
		endif 
		if(myid .eq. id_rootk_strich.and. myid <= n_k_points) then
		  if(real_eigenvectors) then
		    KS_eigenvector_k_strich(:,:,:)=&
                         KS_eigenvector_in(:,&
                         n_state_min_in:n_state_max_in,:,i_k_strich_point_local)
		  else
		    KS_eigenvector_complex_k_strich(:,:,:)=&
                         KS_eigenvector_complex_in(:,&
                         n_state_min_in:n_state_max_in,:,i_k_strich_point_local)
		  endif
		    call mpi_recv(KS_eigenvector_complex_k,&
                         n_basis*size_KS_vec*n_spin, MPI_COMPLEX16, id_rootk,&
                                 111,mpi_comm_global, MPI_STATUS_IGNORE, mpierr)
		    call mpi_wait( i_req(id_rootk), MPI_STATUSES_IGNORE, mpierr)
		    call calc_coulelement(coulelement_tmp, &
                                          coulmat_full_w_k(i_k_point, &
                                                          i_k_strich_point,:),&
                                          KS_eigenvector_k ,&
                                          KS_eigenvector_complex_k,&
		                          KS_eigenvector_k_strich, &
                                          KS_eigenvector_complex_k_strich, &
                                          n_state_min_in, n_state_max_in)
                    call out_coulelement_kart(coulelement_tmp, i_q_point, &
                    i_k_point, i_k_strich_point,"Coulomb_matrix",file_id, &
                    plist_id, n_state_min_in, n_state_max_in)
		endif
	    enddo
	  enddo
          call clean_coulmat_w
	  call clean_coulelement
          call mpi_barrier(mpi_comm_global,info)
	enddo
        if (allocated(q_point_list))deallocate(q_point_list)
        if (allocated(KS_eigenvector_k))deallocate(KS_eigenvector_k)
        if (allocated(KS_eigenvector_k_strich))&
                deallocate(KS_eigenvector_k_strich)
        if (allocated(KS_eigenvector_complex_k))&
                deallocate(KS_eigenvector_complex_k)
        if (allocated(KS_eigenvector_complex_k_strich))&
                deallocate(KS_eigenvector_complex_k_strich)
        call mpi_barrier(mpi_comm_global,info)
        call close_hdf5(file_id,plist_id)
end subroutine get_coulelement_ovl
