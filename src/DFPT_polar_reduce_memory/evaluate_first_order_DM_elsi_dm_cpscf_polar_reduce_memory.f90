!****s* FHI-aims/evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory
!  NAME
!    evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory(   &
           first_order_H, first_order_density_matrix, n_matrix_size)
               

!  PURPOSE
!    calculate the first-order DM for polarizability 

!  shanghui,2013.12.12
!  USES

  use dimensions
  use runtime_choices
  use physics, only: overlap_matrix, hamiltonian 

  use synchronize_mpi_basic, only: sync_vector
  use scalapack_wrapper, only : construct_hamiltonian_real_for_elsi_scalapack, & 
                                construct_overlap_real_for_elsi_scalapack,     &
                                set_full_matrix_real, &
                                get_first_order_dm_polar_reduce_memory_for_elsi_scalapack, &
                                mxld, mxcol, &
                                ! wyj
                                ovlp, ham_stored, &
                                get_set_full_local_matrix_scalapack_cpscf, &
                                get_set_full_local_matrix_scalapack, &
                                print_ham_cpscf, &
                                first_order_ham_polar_reduce_memory_scalapack

  use aims_memory_tracking, only: aims_allocate, aims_deallocate 
  use elsi, only:   elsi_dm_real_cpscf
  use elsi_wrapper, only: eh_scf

  ! wyj
  use load_balancing
  use localorb_io
  use mpi_tasks

!  ARGUMENTS

  implicit none
  

  real*8, dimension( n_matrix_size, n_spin), intent(IN) :: first_order_H
  real*8, dimension( n_matrix_size, n_spin), intent(INOUT) :: first_order_density_matrix
  ! wyj add
  integer :: n_matrix_size

!  INPUTS
!   o first_order_H -- Kohn-Sham eigenvectors (real format)
!   o first_order_density_matrix -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_density_matrix

   real*8, dimension(:,:,:), allocatable   :: ham_local
   real*8, dimension(:,:),   allocatable   :: ovlp_local
   real*8, dimension(:,:,:), allocatable   :: dm
   real*8, dimension(:,:,:), allocatable   :: first_order_ham
   real*8, dimension(:,:),   allocatable   :: first_order_ovlp
   real*8, dimension(:,:,:), allocatable   :: first_order_dm
   real*8  energy
   integer i_spin
   real*8 time0, time_work, time_all, time2, time_work2
   character*200 :: info_str


   if (use_local_index .and. use_load_balancing) then
       write(info_str,'(2X,A)') "Evaluating first-order-DM matrix with &
           ELSI(local_index and load_balancing)"
       call localorb_info(info_str, use_unit,'(A)',OL_norm)

       call aims_allocate(dm,  mxld, mxcol, n_spin, "+dm")
       !call aims_allocate(first_order_ham,  mxld, mxcol, n_spin, "+first_order_ham")
       call aims_allocate(first_order_ovlp, mxld, mxcol, "+first_order_ovlp")
       !call aims_allocate(first_order_dm,  mxld, mxcol, n_spin, "+first_order_dm")



       !---(1) from sparse to scalapack----------------
       !call construct_overlap_real_for_elsi_scalapack(overlap_matrix, ovlp_local)
       !call construct_hamiltonian_real_for_elsi_scalapack(hamiltonian, ham_local)

       first_order_ovlp = 0.0d0
       !call construct_hamiltonian_real_for_elsi_scalapack(first_order_H, first_order_ham)
       !first_order_ham = first_order_ham_polar_reduce_memory_scalapack

       !---(2) from scalapack to elsi sparse--------
       do i_spin = 1, n_spin   
       time0 = mpi_wtime()
       call get_set_full_local_matrix_scalapack(hamiltonian, 1, i_spin)
       if(myid .eq. 0) print*, "myid=", myid, &
            " evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory:get_set_full_local_matrix_scalapack=", &
            mpi_wtime() - time0

       time0 = mpi_wtime()
       call set_full_matrix_real(ovlp)
       if(myid .eq. 0) print*, "myid=", myid, &
            " evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory:set_full_matrix_real=", &
            mpi_wtime() - time0

       !wyj:TODO debug
       !print *, 'print_ovlp'
       !call print_ham_cpscf(ovlp)

       time0 = mpi_wtime()
       call set_full_matrix_real(ham_stored(:,:,i_spin))
       if(myid .eq. 0) print*, "myid=", myid, &
            " evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory:set_full_matrix_real=", &
            mpi_wtime() - time0
       !wyj:TODO debug
       !print *, 'print_ham'
       !call print_ham_cpscf(ham_stored)

       time0 = mpi_wtime()
       call set_full_matrix_real(first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin))
       if(myid .eq. 0) print*, "myid=", myid, &
            " evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory:set_full_matrix_real=", &
            mpi_wtime() - time0
       !print *, 'print_ham_polar_re'
       !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

       !---(3) get first_order_DM----------
       !call elsi_dm_real_cpscf(eh_scf,ham_stored(:,:,i_spin),ovlp,dm(:,:,i_spin), energy, & 
       !    first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin), first_order_ovlp, first_order_dm(:,:,i_spin)  )

       time0 = mpi_wtime()
       call elsi_dm_real_cpscf(eh_scf,ham_stored(:,:,i_spin),ovlp,dm(:,:,i_spin), energy, & 
           first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin), first_order_ovlp, first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin)  )
       if(myid .eq. 0) print*, "myid=", myid, &
            " evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory:elsi_dm_real_cpscf=", &
            mpi_wtime() - time0
       !---(4) scalapack matrix to global matrix---------
       !call get_first_order_dm_polar_reduce_memory_for_elsi_scalapack( first_order_dm(:,:,i_spin), & 
       !     first_order_density_matrix(:,i_spin))
       !first_order_ham_polar_reduce_memory_scalapack = first_order_ham

       time0 = mpi_wtime()
       call get_set_full_local_matrix_scalapack_cpscf(first_order_density_matrix, 2, i_spin)
       if(myid .eq. 0) print*, "myid=", myid, &
            " evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory:get_set_full_local_matrix_scalapack_cpscf=", &
            mpi_wtime() - time0
       enddo ! i_spin

      !---(5) sync the first_order_DM------------
      !call sync_vector(first_order_density_matrix,n_hamiltonian_matrix_size*n_spin) 

      !if(allocated(ham_local))  call aims_deallocate(ham_local, "+ham_local")
      !if(allocated(ovlp_local)) call aims_deallocate(ovlp_local, "+ovlp_local")
      if(allocated(dm))   call aims_deallocate(dm, "+dm")
      !if(allocated(first_order_ham))  call aims_deallocate(first_order_ham, "+first_order_ham")
      if(allocated(first_order_ovlp)) call aims_deallocate(first_order_ovlp, "+first_order_ovlp")
      !if(allocated(first_order_dm))   call aims_deallocate(first_order_dm, "+first_order_dm")
  else
      write(info_str,'(2X,A)') "Evaluating first-order-DM matrix" 
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      call aims_allocate(ham_local,  mxld, mxcol, n_spin, "+ham_local")
      call aims_allocate(ovlp_local, mxld, mxcol, "+ovlp_local")
      call aims_allocate(dm,  mxld, mxcol, n_spin, "+dm")
      call aims_allocate(first_order_ham,  mxld, mxcol, n_spin, "+first_order_ham")
      call aims_allocate(first_order_ovlp, mxld, mxcol, "+first_order_ovlp")
      call aims_allocate(first_order_dm,  mxld, mxcol, n_spin, "+first_order_dm")


      !---(1) from sparse to scalapack----------------
      call construct_overlap_real_for_elsi_scalapack(overlap_matrix, ovlp_local)
      call construct_hamiltonian_real_for_elsi_scalapack(hamiltonian, ham_local)

      first_order_ovlp = 0.0d0
      call construct_hamiltonian_real_for_elsi_scalapack(first_order_H, first_order_ham)

      !---(2) from scalapack to elsi sparse--------
      do i_spin = 1, n_spin   
      call set_full_matrix_real(ovlp_local)
      !wyj:TODO debug
      !print *, 'print_ovlp'
      !call print_ham_cpscf(ovlp_local)

      call set_full_matrix_real(ham_local(:,:,i_spin))
      !wyj:TODO debug
      !print *, 'print_ham'
      !call print_ham_cpscf(ham_local)

      call set_full_matrix_real(first_order_ham(:,:,i_spin))
       !wyj:TODO debug
       !print *, 'print_ham_polar_re'
       !call print_ham_cpscf(first_order_ham)

      !---(3) get first_order_DM----------
      call elsi_dm_real_cpscf(eh_scf,ham_local(:,:,i_spin),ovlp_local,dm(:,:,i_spin), energy, & 
          first_order_ham(:,:,i_spin), first_order_ovlp, first_order_dm(:,:,i_spin)  )

      !---(4) scalapack matrix to global matrix---------
      call get_first_order_dm_polar_reduce_memory_for_elsi_scalapack( first_order_dm(:,:,i_spin), & 
          first_order_density_matrix(:,i_spin))
      enddo ! i_spin

      !---(5) sync the first_order_DM------------
      call sync_vector(first_order_density_matrix,n_hamiltonian_matrix_size*n_spin) 

      if(allocated(ham_local))  call aims_deallocate(ham_local, "+ham_local")
      if(allocated(ovlp_local)) call aims_deallocate(ovlp_local, "+ovlp_local")
      if(allocated(dm))   call aims_deallocate(dm, "+dm")
      if(allocated(first_order_ham))  call aims_deallocate(first_order_ham, "+first_order_ham")
      if(allocated(first_order_ovlp)) call aims_deallocate(first_order_ovlp, "+first_order_ovlp")
      if(allocated(first_order_dm))   call aims_deallocate(first_order_dm, "+first_order_dm")

  endif


end subroutine evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory
!******
