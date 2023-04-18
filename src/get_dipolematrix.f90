!****s* FHI-aims/get_angular_grid
!  NAME
!    get_angular_grid
!  SYNOPSIS

subroutine get_dipolematrix &
     (KS_eigen, KS_vec, KS_vec_complex, occ_numbers, chemical_potential, &
      partition_tab, l_shell_max)

!  PURPOSE
!
!  Wrapper function for calculating and outputting of Dipolematrix elements
!
!  USES
  use calculate_mommat_base
  use dimensions
  use runtime_choices
  use localorb_io
  use hdf5_tools, only: HID_T, open_hdf5, open_hdf5_kart, out_occ_kart, &
      out_k_points_kart, out_bands_kart, outmetafile, out_moment_kart, &
      close_hdf5
  use mpi_utilities
  use boys, only: boys_centers, boys_transform_subspace, boys_sub_flags
  use synchronize_mpi_basic, only: bcast_complex_vector
!  ARGUMENTS

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN)::  KS_vec
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, INTENT(IN) :: chemical_potential
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 
!  real*8 :: start, finish ! Time debug

!  INPUTS
!   o KS_eigen
!   o KS_vec/KS_vec_complex
!   o occ_numbers
!   o chemical_potential
!   o partition_tab
!   o l_shell_max
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

  INTEGER(HID_T) :: file_id                  ! File identifier
  INTEGER(HID_T) :: plist_id                 ! Property list identifier

  character*128 :: info_str
  integer :: info
  integer :: n_state_min_in
  integer :: n_state_max_in
  !  counters

  integer :: i_k
  integer :: new_k_point
  integer :: lenmom
  !  begin work

!    call cpu_time(start)

    write(info_str,'(6X,A,1X,I4)') "Dipole Matrix post processing starts"
    call localorb_info ( info_str )
    call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)
    if(k_point_mom==0)then
        call open_hdf5('dipmat.h5', file_id, plist_id)
        call open_hdf5_kart("_Dipolematrix_", file_id, plist_id,n_state_min_in,&
                            n_state_max_in)
        call mpi_barrier(mpi_comm_global,info)
        call out_bands_kart(KS_eigen,file_id,plist_id,n_state_min_in,&
                            n_state_max_in)
        call out_occ_kart(occ_numbers,file_id,plist_id,n_state_min_in,&
                            n_state_max_in)
        call out_k_points_kart(file_id,plist_id)
        call outmetafile(KS_eigen,occ_numbers,chemical_potential,file_id,&
                         plist_id,n_state_min_in,n_state_max_in)
   endif
   call allocate_mommat()
   call allocate_mommat_k()
   !if (k_point_mom==0)then
   !   call allocate_moment_one(n_state_min_in, n_state_max_in)
   !else
      call allocate_moment_one(n_state_min_in, n_state_max_in)
      call allocate_moment_two(n_state_min_in, n_state_max_in)
      call allocate_moment_three(n_state_min_in, n_state_max_in)
   !endif

!   call cpu_time(finish)
!   write(use_unit,'(A, 1X, F6.3)') "Time alloc = ", finish-start
!   call cpu_time(start)

   call calculate_dipmat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
                              mommat_full_oned_low,1)

!   call cpu_time(finish)
!   write(use_unit,'(A, 1X, F6.3)') "Time 1 = ", finish-start

   lenmom = (((n_state_max_in-n_state_min_in+1)+1)* &
           &  (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)

   i_k = 0
   do new_k_point = 1,n_k_points
      if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
        ! k-point stored on current mpi task
        i_k = i_k+1 !new_k_point
        if (new_k_point==k_point_mom.or.k_point_mom==0)then
          call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
                                  mommat_full_w_complex_up, new_k_point, &
                                  work_ovl_mom )
          call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
                                  mommat_full_w_complex_low, new_k_point, &
                                  work_ovl_mom )
          call calc_moment_p0(moment_one,mommat_full_w_up, mommat_full_w_low, &
                           mommat_full_w_complex_up, mommat_full_w_complex_low,&
                           KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
                           new_k_point, 1, n_state_min_in, n_state_max_in)
          if (k_point_mom==0) then
            call out_moment_kart(moment_one,KS_eigen,new_k_point,&
                                  "_Dipolematrix_", file_id, plist_id,1,&
                                  n_state_min_in,n_state_max_in)
          endif
        endif
      endif
      if (apply_boys_flag) then
        call bcast_complex_vector(moment_one, lenmom, MOD(new_k_point, n_tasks))
      end if
   enddo
   call calculate_dipmat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, &
                              mommat_full_oned_low,2)

!   call cpu_time(start)   
!   call cpu_time(finish)
!   write(use_unit,'(A, 1X, F6.3)') "Time 2 = ", finish-start
!   call cpu_time(start)

   i_k = 0
   do new_k_point = 1,n_k_points
     if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
       ! k-point stored on current mpi task
       i_k = i_k+1 !new_k_point
       if (new_k_point==k_point_mom.or.k_point_mom==0)then
         call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
                                  mommat_full_w_complex_up, new_k_point, &
                                  work_ovl_mom )
         call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
                                  mommat_full_w_complex_low, new_k_point, &
                                  work_ovl_mom )
         call calc_moment_p0(moment_two,mommat_full_w_up, mommat_full_w_low, &
                           mommat_full_w_complex_up, mommat_full_w_complex_low,&
                            KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
                            new_k_point, 2, n_state_min_in, n_state_max_in)
         if (k_point_mom==0)then
           call out_moment_kart(moment_two,KS_eigen,new_k_point,&
                                 "_Dipolematrix_", file_id, plist_id,2, &
                                  n_state_min_in,n_state_max_in)
         endif
       endif
     endif
     if (apply_boys_flag) then
        call bcast_complex_vector(moment_two, lenmom, MOD(new_k_point, n_tasks))
     end if
   enddo
   call calculate_dipmat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
                               mommat_full_oned_low,3)
!   call cpu_time(finish)
!   write(use_unit,'(A, 1X, F6.3)') "Time 3 = ", finish-start
!   call cpu_time(start)

    i_k = 0
    do new_k_point = 1,n_k_points
      if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
        ! k-point stored on current mpi task
        i_k = i_k+1 !new_k_point
        if (new_k_point==k_point_mom.or.k_point_mom==0)then
          call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
                                  mommat_full_w_complex_up, new_k_point, &
                                  work_ovl_mom )
          call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
                                  mommat_full_w_complex_low, new_k_point, &
                                  work_ovl_mom )
          call calc_moment_p0(moment_three,mommat_full_w_up, mommat_full_w_low,&
                           mommat_full_w_complex_up, mommat_full_w_complex_low,&
                           KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
                           new_k_point, 3, n_state_min_in, n_state_max_in)
          if (k_point_mom==0)then 
            call out_moment_kart(moment_three,KS_eigen,new_k_point,&
                                 "_Dipolematrix_", file_id, plist_id,3, &
                                 n_state_min_in,n_state_max_in)
          endif
        endif
      endif
      if (apply_boys_flag) then
         call bcast_complex_vector(moment_three, lenmom, MOD(new_k_point, n_tasks))
      end if
    enddo


    ! GSM: This is purely a debug print to check the dipole matrix elements ----
    !
    ! (complex values of upper triangular matrix)
    !
    ! Do not use these values for computations, they are truncated.
    ! This is purely for sanity checking..
    !
    ! open(unit=120, file='dipole_matrix_elements.csv',ACTION='WRITE')
    ! write(120, '(2X, A)') '# dipx dipy dipz '
    ! do i_dip=1, (((n_state_max_in-n_state_min_in+1)+1)* &
    ! & (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2)
    !   write(120,'(2X,3(2X,F12.8,:,SP,F12.8,"j"))')  moment_one(i_dip), &
    !   &                              moment_two(i_dip), &
    !   &                              moment_three(i_dip)
    ! enddo
    ! close(unit=120)
    !
    ! --------------------------------------------------------------------------

! Here all necessary ingredients for the calculation of the localized
! Boys orbitals have been calculated. Now it is just optimizing the
! functional from the external routine and outputting it. Boys orbitals
! are for non periodic systems only (k_point_mom=1). Implementation of
! periodic Wannier orbitals shouldn't be so hard... Also, all occupied
! orbitals should be taken into account! How to always consider all of them?
    if (flag_out_boys) then
       if (k_point_mom/=0.and.myid==modulo(k_point_mom, n_tasks))then
          write(info_str,'(6X,A,1X,I4)') 'Boys centers and rotation matrix calculation'
          call localorb_info ( info_str )
          call boys_centers(moment_one, moment_two, moment_three, n_state_max_in) 
       endif
    else if (apply_boys_flag) then
          write(info_str,'(6X,A,1X,I4)') 'Calculating and applying Boys matrix ...'
          call localorb_info ( info_str )
          call boys_transform_subspace(moment_one, moment_two, moment_three, n_state_max_in)
          if (maxval(boys_sub_flags) .ne. 1) then
              write(info_str,'(6X,A,1X,I4)') 'Deactivating Boys-calculation ...'
              call localorb_info ( info_str )
              apply_boys_flag = .false.
          end if
    else ! just output the momentum matrix
       if (k_point_mom/=0.and.myid==modulo(k_point_mom, n_tasks))then 
           call out_moment(moment_one,moment_two,moment_three, KS_eigen, &
                       k_point_mom,n_state_min_in,n_state_max_in)
       endif
    endif

!   call cpu_time(finish)
!   write(use_unit,'(A, 1X, F6.3)') "Time boys = ", finish-start

    call clean_mommat()
    call clean_mommat_final()
    if(k_point_mom==0)then
      call mpi_barrier(mpi_comm_global,info)
      call close_hdf5(file_id,plist_id)
    endif
    write(info_str,'(6X,A,1X,I4)') "Dipole Matrix post processing finished"
    call localorb_info ( info_str )
 
end subroutine get_dipolematrix
!******	
