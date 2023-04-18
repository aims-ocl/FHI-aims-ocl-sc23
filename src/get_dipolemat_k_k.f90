subroutine get_dipolemat_k_k(KS_eigen, KS_eigenvector_in,&
                             KS_eigenvector_complex_in, occ_numbers, &
                             chemical_potential,partition_tab_std_in, &
                             basis_l_max_in, do_k_k)
!  PURPOSE
!
!  Wrapper function for calculating and outputting the dipolematrix depending 
!  on k and k'. 3 methods are available, that differ in the maximum matrix
!  size.
!  do_k_k == 1
!  max = n_basis*(n_basis+1)/2)*n_k_ponts*n_k_points
!  do_k_k == 2
!  max = n_basis*(n_basis+1)/2)*n_k_ponts
!  do_k_k == 3
!  max = n_basis*(n_basis+1)/2)
!  Memory needs reduce from 1 to 3, but more (repetitive) calculations are 
!  necessary
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
      use calculate_dipolemat_k_k
      use hdf5_tools, only: HID_T, open_hdf5, open_coulelement_kart, &
          out_bands_kart, out_k_points_kart, outmetafile, close_hdf5, &
          out_coulelement_kart
      use localorb_io, only: localorb_info, use_unit
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
      real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
      real*8, INTENT(IN) :: chemical_potential
      integer, INTENT(IN) :: do_k_k
!  INPUTS
!   o partition_tab -- Partition tab
!   o l_shell_max
!   o KS_eigen -- KS eigenvalues
!   o KS_vec/KS_vec_complex -- KS coefficients
!   o occ_numbers -- Occupation numbers
!   o chemical_potential -- E_F
!   o  do_k_k -- Switch to choose the method to calculate the dipolematrix(k,k')
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

      complex*16, allocatable :: KS_eigenvector_complex_k(:,:,:) 
      complex*16, allocatable :: KS_eigenvector_complex_k_strich(:,:,:) 
      real*8, allocatable :: KS_eigenvector_k(:,:,:) 
      real*8, allocatable :: KS_eigenvector_k_strich(:,:,:) 

      INTEGER(HID_T) :: file_id                  ! File identifier
      INTEGER(HID_T) :: plist_id                 ! Property list identifier

  !  counters
      integer::  i_coord
      integer::  i_k_point, i_k_strich_point
      integer::   id_rootk, id_rootk_strich
      integer::  i_k_point_local, i_k_strich_point_local
      integer::  i_state, j_state, num 
      integer::  n_state_min_in, n_state_max_in
      integer::  size_element, size_KS_vec
      integer :: info, mpierr
      character(*), parameter :: func = 'get_dipelement_k_k'
      character*150 :: info_str

      integer :: i_req(0:n_tasks-1)

 
      write(info_str,'(6X,A,1X,I4)') "Momentum Matrix (k,k') post processing starts"
      call localorb_info ( info_str ,use_unit)

      call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)	

      if (do_k_k==1) then
	  call open_hdf5('dipmat_k_k.h5', file_id, plist_id)
	  call open_coulelement_kart("Dipole_matrix_k_k", file_id, plist_id,&
				n_state_min_in, n_state_max_in, 3)
	  call out_bands_kart(KS_eigen,file_id,plist_id,n_state_min_in,&
			      n_state_max_in)
	  call out_k_points_kart(file_id,plist_id)
	  call outmetafile(KS_eigen,occ_numbers,chemical_potential,file_id,&
			  plist_id,n_state_min_in,n_state_max_in)
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
	      allocate(KS_eigenvector_complex_k_strich(n_basis,size_KS_vec,&
							n_spin), stat=info) 
	      call check_allocation(info,'KS_eigenvector_complex_k_strich',func)
	  endif
	  do i_coord=1, 3, 1
	      call allocate_dipmat_k_k
	      dipmat_full_k_k=(0.0, 0.0)
	      call calculate_dipmat_k_k( partition_tab_std_in, basis_l_max_in, &
				      dipmat_full_k_k,i_coord)  
	      call allocate_dipelement_k_k(n_state_min_in, n_state_max_in)
	      do i_k_point = 1, n_k_points, 1
		id_rootk = mod(i_k_point,n_tasks)
		i_k_point_local = (i_k_point-1)/n_tasks + 1 
		i_req(:) = MPI_REQUEST_NULL
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
				      n_basis*size_KS_vec*n_spin,MPI_COMPLEX16,&
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
			call mpi_wait( i_req(id_rootk), &
                                       MPI_STATUSES_IGNORE, mpierr)
			call calc_dipelement_k_k(dipelement_k_k, &
					      dipmat_full_k_k(i_k_point, &
							   i_k_strich_point,:),&
					      KS_eigenvector_k ,&
					      KS_eigenvector_complex_k,&
					      KS_eigenvector_k_strich, &
					      KS_eigenvector_complex_k_strich, &
					      n_state_min_in, n_state_max_in)
			call out_coulelement_kart(dipelement_k_k, i_coord, &
			   i_k_point, i_k_strich_point,"Dipole_matrix_k_k",&
			   file_id, plist_id, n_state_min_in, n_state_max_in)
		    endif
		enddo
	      enddo
	      call clean_dipmat_k_k
	      call clean_dipelement_k_k
	    enddo
	    if (allocated(KS_eigenvector_k))deallocate(KS_eigenvector_k)
	    if (allocated(KS_eigenvector_k_strich))&
		    deallocate(KS_eigenvector_k_strich)
	    if (allocated(KS_eigenvector_complex_k))&
		    deallocate(KS_eigenvector_complex_k)
	    if (allocated(KS_eigenvector_complex_k_strich))&
		    deallocate(KS_eigenvector_complex_k_strich)
	    call mpi_barrier(mpi_comm_global,info)
	    call close_hdf5(file_id,plist_id)
      elseif(do_k_k==2) then
	  call open_hdf5('dipmat_k.h5', file_id, plist_id)
	  call open_coulelement_kart("Dipole_matrix_k_k", file_id, plist_id,&
				n_state_min_in, n_state_max_in, 3)
	  call out_bands_kart(KS_eigen,file_id,plist_id,n_state_min_in,&
			      n_state_max_in)
	  call out_k_points_kart(file_id,plist_id)
	  call outmetafile(KS_eigen,occ_numbers,chemical_potential,file_id,&
			  plist_id,n_state_min_in,n_state_max_in)
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
	      allocate(KS_eigenvector_complex_k_strich(n_basis,size_KS_vec,&
							 n_spin), stat=info) 
	      call check_allocation(info,'KS_eigenvector_complex_k_strich',func)
	  endif
	  do i_coord=1, 3, 1
	      if (myid==0) then 
		write(use_unit,*) 'i_coord: ', i_coord 
	      endif
	      call allocate_dipmat_k()
	      call allocate_dipelement_k(size_element)
	      do i_k_point = 1, n_k_points, 1
		if (myid==0) then
		  write(use_unit,*) 'k_point: ', i_k_point
		endif
	        call calculate_dipmat_k( partition_tab_std_in, basis_l_max_in, &
				      dipmat_full_k,i_coord,i_k_point) 
		id_rootk = mod(i_k_point,n_tasks)
		i_k_point_local = (i_k_point-1)/n_tasks + 1 
		dipelement_k = (0.0,0.0)
		i_req(:) = MPI_REQUEST_NULL
		do i_k_strich_point = 1, n_k_points, 1
		    id_rootk_strich = mod(i_k_strich_point,n_tasks)
		    i_k_strich_point_local = (i_k_strich_point-1)/n_tasks + 1
	    
		    if(myid .eq. id_rootk.and. myid <= n_k_points) then
		      if(real_eigenvectors) then
			KS_eigenvector_k(:,:,:)=&
				KS_eigenvector_complex_in(:,&
			        n_state_min_in:n_state_max_in,:,i_k_point_local)
			call mpi_bcast(KS_eigenvector_k,&
				      n_basis*size_KS_vec*n_spin,MPI_COMPLEX16,&
				      id_rootk, mpi_comm_global, mpierr)
		      else
			KS_eigenvector_complex_k(:,:,:)=&
				    KS_eigenvector_complex_in(:,&
				    n_state_min_in:n_state_max_in,:,i_k_point_local)
			call mpi_isend(KS_eigenvector_complex_k,&
				      n_basis*size_KS_vec*n_spin,MPI_COMPLEX16,&
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
			call mpi_wait( i_req(id_rootk), &
                                       MPI_STATUSES_IGNORE, mpierr)
			call calc_dipelement_k_k(dipelement_k, &
					     dipmat_full_k(i_k_strich_point,:),&
					     KS_eigenvector_k ,&
					     KS_eigenvector_complex_k,&
					     KS_eigenvector_k_strich, &
					     KS_eigenvector_complex_k_strich, &
					     n_state_min_in, n_state_max_in)

			call out_coulelement_kart(dipelement_k, i_coord, &
			 i_k_point, i_k_strich_point,"Dipole_matrix_k_k",&
			 file_id, plist_id, n_state_min_in, n_state_max_in)
		    endif
		enddo
	      enddo
	      call clean_dipmat_k
	      call clean_dipelement_k
	  enddo
	  if (allocated(KS_eigenvector_k))deallocate(KS_eigenvector_k)
	  if (allocated(KS_eigenvector_k_strich))&
	                deallocate(KS_eigenvector_k_strich)
	  if (allocated(KS_eigenvector_complex_k))&
	                deallocate(KS_eigenvector_complex_k)
	  if (allocated(KS_eigenvector_complex_k_strich))&
	                deallocate(KS_eigenvector_complex_k_strich)
	  call mpi_barrier(mpi_comm_global,info)
	    call close_hdf5(file_id,plist_id)
      elseif(do_k_k==3) then
	  call open_hdf5('dipmat_k.h5', file_id, plist_id)
	  call open_coulelement_kart("Dipole_matrix_k_k", file_id, plist_id,&
				n_state_min_in, n_state_max_in, 3)
	  call out_bands_kart(KS_eigen,file_id,plist_id,n_state_min_in,&
			      n_state_max_in)
	  call out_k_points_kart(file_id,plist_id)
	  call outmetafile(KS_eigen,occ_numbers,chemical_potential,file_id,&
			  plist_id,n_state_min_in,n_state_max_in)
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
	      allocate(KS_eigenvector_complex_k_strich(n_basis,size_KS_vec,&
							     n_spin), stat=info) 
	      call check_allocation(info, 'KS_eigenvector_complex_k_strich',&
                                    func)
	  endif
	  do i_coord=1, 3, 1
	      if (myid==0) then 
		write(use_unit,*) 'i_coord: ', i_coord 
	      endif
	      call allocate_dipmat_one_k()
	      call allocate_dipelement_k(size_element)
	      do i_k_point = 1, n_k_points, 1
		if (myid==0) then
		  write(use_unit,*) 'k_point: ', i_k_point
		endif
		id_rootk = mod(i_k_point,n_tasks)
		i_k_point_local = (i_k_point-1)/n_tasks + 1 
		dipelement_k = (0.0,0.0)
		i_req(:) = MPI_REQUEST_NULL
		do i_k_strich_point = 1, n_k_points, 1
		    if (myid==0) then
		      write(use_unit,*) 'k_strich_point: ', i_k_strich_point
		    endif
	            call calculate_dipmat_one_k( partition_tab_std_in, &
                         basis_l_max_in, dipmat_full_one_k,i_coord,i_k_point,&
                         i_k_strich_point) 
		    id_rootk_strich = mod(i_k_strich_point,n_tasks)
		    i_k_strich_point_local = (i_k_strich_point-1)/n_tasks + 1
	    
		    if(myid .eq. id_rootk.and. myid <= n_k_points) then
		      if(real_eigenvectors) then
			KS_eigenvector_k(:,:,:)=&
				    KS_eigenvector_complex_in(:,&
				    n_state_min_in:n_state_max_in,:,&
                                    i_k_point_local)
			call mpi_bcast(KS_eigenvector_k,&
				      n_basis*size_KS_vec*n_spin,MPI_COMPLEX16,&
				      id_rootk, mpi_comm_global, mpierr)
		      else
			KS_eigenvector_complex_k(:,:,:)=&
				    KS_eigenvector_complex_in(:,&
				    n_state_min_in:n_state_max_in,:,&
                                    i_k_point_local)
			call mpi_isend(KS_eigenvector_complex_k,&
				      n_basis*size_KS_vec*n_spin,MPI_COMPLEX16,&
				      id_rootk_strich, 111, mpi_comm_global, &
				      i_req(id_rootk_strich), mpierr)
		      endif
		    endif 
		    if(myid .eq. id_rootk_strich.and. myid <= n_k_points) then
		      if(real_eigenvectors) then
			KS_eigenvector_k_strich(:,:,:)=&
			    KS_eigenvector_in(:,&
			    n_state_min_in:n_state_max_in,:,&
                            i_k_strich_point_local)
		      else
			KS_eigenvector_complex_k_strich(:,:,:)=&
			    KS_eigenvector_complex_in(:,&
			    n_state_min_in:n_state_max_in,:,&
                            i_k_strich_point_local)
		      endif
			call mpi_recv(KS_eigenvector_complex_k,&
			   n_basis*size_KS_vec*n_spin, MPI_COMPLEX16, id_rootk,&
		           111,mpi_comm_global, MPI_STATUS_IGNORE, mpierr)
			call mpi_wait( i_req(id_rootk), MPI_STATUSES_IGNORE, &
                                       mpierr)
			call calc_dipelement_k_k(dipelement_k, &
					      dipmat_full_one_k(:),&
					      KS_eigenvector_k ,&
					      KS_eigenvector_complex_k,&
					      KS_eigenvector_k_strich, &
					      KS_eigenvector_complex_k_strich, &
					      n_state_min_in, n_state_max_in)

			call out_coulelement_kart(dipelement_k, i_coord, &
			 i_k_point, i_k_strich_point,"Dipole_matrix_k_k",&
			 file_id, plist_id, n_state_min_in, n_state_max_in)
		    endif
		enddo
	      enddo
	      call clean_dipmat_one_k
	      call clean_dipelement_k
	  enddo
	  if (allocated(KS_eigenvector_k))deallocate(KS_eigenvector_k)
	  if (allocated(KS_eigenvector_k_strich))&
	                deallocate(KS_eigenvector_k_strich)
	  if (allocated(KS_eigenvector_complex_k))&
	                deallocate(KS_eigenvector_complex_k)
	  if (allocated(KS_eigenvector_complex_k_strich))&
	                deallocate(KS_eigenvector_complex_k_strich)
	  call mpi_barrier(mpi_comm_global,info)
	    call close_hdf5(file_id,plist_id)
      endif
end subroutine get_dipolemat_k_k
