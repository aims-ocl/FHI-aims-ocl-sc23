module redistribute_cpt2_kspace
  use dimensions
  use physics
  use hartree_fock
  use cpt2_blacs
  use lvl_tricoeff_cpt2
  implicit none

  save
  private
  public redistribute_cpt2_blacs_complex, redistribute_cpt2_blacs_real, deallocate_cpt2_blacs

  contains

  subroutine redistribute_cpt2_blacs_complex()
     use aims_memory_tracking, only : aims_deallocate
     use runtime_choices, only: real_eigenvectors

     implicit none
     integer :: info, mpierr
     character*150 :: info_str
     character(*), parameter :: func='redistribute_cpt2_kspace'

     ! for mpi windows
     integer :: win_tri_inblacs, win_ev, win_tri
     integer(kind=MPI_ADDRESS_KIND):: nbytes, offset
     integer:: dcsize, count, count_1, count_2, count_3, k_task
     ! for redistribution of KS_eigenvectors
     real(kind=8), allocatable :: KS_eigenvectors_occ_tmp(:,:,:,:)
     real(kind=8), allocatable :: KS_eigenvectors_unocc_tmp(:,:,:,:)
     complex(kind=8), dimension(:,:,:), pointer :: lvl_tricoeff_p
     integer :: max_k_points_task

     !counter
     integer :: i_spin
     integer :: i_k_point
     integer :: i_k_point_local

     !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the kq loop
     if (mod(n_k_points,n_tasks_kq_cpt2).eq.0) then    
        max_k_points_task=n_k_points/n_tasks_kq_cpt2
     else
        max_k_points_task=n_k_points/n_tasks_kq_cpt2+1
     end if

     ! ==========================================================================
     ! now re-distribute lvl_tricoeff in terms of (un)occupied orbitals for PT2
     ! ==========================================================================

     if (myid_col_cpt2.eq.0) then
       allocate(lvl_tricoeff_occ(lbb_row:ubb_row,max_n_basis_sp,n_low_state_cpt2:n_homo_max,n_spin,n_ks_points_task),stat=info)
       do i_k_point_local = 1, n_ks_points_task
         !i_k_point = n_tasks_kq_cpt2*(i_k_point_local-1)+myid_kq_cpt2+1
         !k_point_loc_cpt2(1,i_k_point)=myid
         !k_point_loc_cpt2(2,i_k_point)=i_k_point_local
         !k_point_loc_cpt2(3,i_k_point)=myid_bl_cpt2
         do i_spin = 1, n_spin, 1
           lvl_tricoeff_occ(lbb_row:ubb_row,:,n_low_state_cpt2:n_homo_max,i_spin,i_k_point_local) = &
               lvl_tricoeff_mod_r(lbb_row:ubb_row,:,n_low_state_cpt2:n_homo_max,i_spin,i_k_point_local)
         end do
       end do
     end if

     allocate(lvl_tricoeff_unocc(lbb_row:ubb_row,max_n_basis_sp,lpb_col:upb_col,n_spin,n_ks_points_task),stat=info)
     
     dcsize=16
     if (myid_kq_cpt2.lt.n_k_points .and. (size(lvl_tricoeff_mod_r) .ne. 0) .and. myid_col_cpt2 .eq. 0) then
        nbytes=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*int(n_states,8)&
            *int(n_spin,8)*int(n_ks_points_task,8)*int(dcsize,8), MPI_ADDRESS_KIND)
     else
        nbytes=int(0,MPI_ADDRESS_KIND)
     end if

     call mpi_win_create(lvl_tricoeff_mod_r,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri,mpierr)
     call mpi_win_fence(0,win_tri,mpierr)

     count_1 = n_bb_row*max_n_basis_sp*n_pb_col
     count_2 = n_bb_row*max_n_basis_sp*n_states*n_spin*n_ks_points_task

     do i_k_point_local = 1, max_k_points_task, 1
       if (i_k_point_local .le. n_ks_points_task) then
         do i_spin = 1, n_spin, 1
           k_task=myid_kq_cpt2*n_tasks_bl_cpt2+myid_bl_cpt2-myid_col_cpt2
           offset = n_bb_row*max_n_basis_sp*n_states*n_spin*(i_k_point_local-1) + &
               n_bb_row*max_n_basis_sp*n_states*(i_spin-1) + &
               n_bb_row*max_n_basis_sp*(lpb_col+n_lumo_min-2)
           call mpi_get(lvl_tricoeff_unocc(:,:,:,i_spin,i_k_point_local), count_1, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_1, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
           call mpi_win_fence(0,win_tri,mpierr)
         end do
       else
         do i_spin = 1, n_spin, 1
           call mpi_win_fence(0,win_tri,mpierr)
         end do
       end if
     end do

     call mpi_win_fence(0,win_tri,mpierr)

     call mpi_win_free(win_tri,mpierr)

     deallocate(lvl_tricoeff_mod_r)

     ! ==========================================================================
     ! now re-distribute KS_eigenvectors in terms of (un)occupied orbitals for PT2
     ! ==========================================================================
     if (myid_col_cpt2.eq.0) then
       allocate(KS_eigenvectors_occ(n_basis,n_low_state_cpt2:n_homo_max,n_spin,n_ks_points_task),stat=info)
     end if
     allocate(KS_eigenvectors_unocc(n_basis,lpb_col:upb_col,n_spin,n_ks_points_task),stat=info)
     !  creat mpi_win for KS_eigenvectors
     if (real_eigenvectors) then
       dcsize=8
       nbytes=int(int(n_basis,8)*int(n_states,8)*int(n_spin,8)*&
           int(n_k_points_task,8)*int(dcsize,8), MPI_ADDRESS_KIND)
       call mpi_win_create(KS_eigenvector,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
       if (myid_col_cpt2.eq.0) then
         allocate(KS_eigenvectors_occ_tmp(n_basis,n_low_state_cpt2:n_homo_max,n_spin,n_ks_points_task),stat=info)
       end if
       allocate(KS_eigenvectors_unocc_tmp(n_basis,lpb_col:upb_col,n_spin,n_ks_points_task),stat=info)
     else 
       if(flag_KS_eigenfunc_conjg) then
          KS_eigenvector_complex = conjg(KS_eigenvector_complex)
       endif
       dcsize=16
       nbytes=int(int(n_basis,8)*int(n_states,8)*int(n_spin,8)*&
           int(n_k_points_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       call mpi_win_create(KS_eigenvector_complex,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
     end if

     call mpi_win_fence(0,win_ev,mpierr)

     count_1=n_basis*n_states
     count_2=n_basis*n_occ_max
     count_3=n_basis*n_pb_col
     do i_k_point_local = 1, max_k_points_task, 1
       if (i_k_point_local .le. n_ks_points_task) then
         i_k_point = n_tasks_kq_cpt2*(i_k_point_local-1)+myid_kq_cpt2+1
         if (myid_col_cpt2.eq.0) then
           do i_spin = 1, n_spin, 1
             k_task=k_point_loc(1,i_k_point)
             offset=(k_point_loc(2,i_k_point)-1)*count_1*n_spin + &
                 count_1*(i_spin-1) + n_basis*(n_low_state_cpt2-1)
             if (real_eigenvectors) then
               call mpi_get(KS_eigenvectors_occ_tmp(:,:,i_spin,i_k_point_local), &
                   count_2, MPI_DOUBLE_PRECISION, k_task, &
                   offset, count_2, MPI_DOUBLE_PRECISION, win_ev, mpierr)
             else
               call mpi_get(KS_eigenvectors_occ(:,:,i_spin,i_k_point_local), &
                   count_2, MPI_DOUBLE_COMPLEX, k_task, &
                  offset, count_2, MPI_DOUBLE_COMPLEX, win_ev, mpierr)
             end if
             call mpi_win_fence(0,win_ev,mpierr)
           end do
         else
           do i_spin = 1, n_spin, 1
             call mpi_win_fence(0,win_ev,mpierr)
           end do
         end if
         do i_spin = 1, n_spin, 1
           k_task=k_point_loc(1,i_k_point)
           offset=(k_point_loc(2,i_k_point)-1)*count_1*n_spin + &
               count_1*(i_spin-1) + n_basis*(lpb_col+n_lumo_min-2)
           if (real_eigenvectors) then
             call mpi_get(KS_eigenvectors_unocc_tmp(:,:,i_spin,i_k_point_local), &
                 count_3, MPI_DOUBLE_PRECISION, k_task, &
                offset, count_3, MPI_DOUBLE_PRECISION, win_ev, mpierr)
           else
             call mpi_get(KS_eigenvectors_unocc(:,:,i_spin,i_k_point_local), &
                 count_3, MPI_DOUBLE_COMPLEX, k_task, &
                offset, count_3, MPI_DOUBLE_COMPLEX, win_ev, mpierr)
           end if
           call mpi_win_fence(0,win_ev,mpierr)
         end do
       else
         do i_spin = 1, n_spin, 1
           call mpi_win_fence(0,win_ev,mpierr)
         end do
         do i_spin = 1, n_spin, 1
           call mpi_win_fence(0,win_ev,mpierr)
         end do
       end if
     end do

     call mpi_win_fence(0,win_ev,mpierr)

     call mpi_win_free(win_ev,mpierr)

     if (real_eigenvectors) then
         if (myid_col_cpt2 .eq. 0) then
             KS_eigenvectors_occ = KS_eigenvectors_occ_tmp
         end if
         KS_eigenvectors_unocc = KS_eigenvectors_unocc_tmp
     end if

     if (real_eigenvectors) then
       if (myid_col_cpt2.eq.0) then
         deallocate(KS_eigenvectors_occ_tmp, KS_eigenvectors_unocc_tmp)
       else
         deallocate(KS_eigenvectors_unocc_tmp)
       end if
         call aims_deallocate(KS_eigenvector, "KS_eigenvector")
     else
         call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex")
     end if

   end subroutine redistribute_cpt2_blacs_complex

  subroutine redistribute_cpt2_blacs_real()

     implicit none
     integer :: info, mpierr
     character*150 :: info_str
     character(*), parameter :: func='redistribute_cpt2_kspace'

     ! for mpi windows
     integer :: win_tri_inblacs, win_ev, win_tri
     integer(kind=MPI_ADDRESS_KIND):: nbytes, offset
     integer:: dcsize, count, count_1, count_2, count_3, k_task
     ! for redistribution of KS_eigenvectors
     complex(kind=8), allocatable :: lvl_tricoeff_unocc_tmp(:,:,:)

     integer :: max_k_points_task

     !counter
     integer :: i_spin
     integer :: i_k_point
     integer :: i_k_point_local

     !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the kq loop
     if (mod(n_k_points,n_tasks_kq_cpt2).eq.0) then    
        max_k_points_task=n_k_points/n_tasks_kq_cpt2
     else
        max_k_points_task=n_k_points/n_tasks_kq_cpt2+1
     end if

     ! ==========================================================================
     ! now re-distribute lvl_tricoeff in terms of (un)occupied orbitals for PT2
     ! ==========================================================================

     if (myid_col_cpt2.eq.0) then
       allocate(lvl_tricoeff_occ_real(lbb_row:ubb_row,&
           max_n_basis_sp,n_low_state_cpt2:n_homo_max,n_spin,n_ks_points_task),stat=info)
       do i_k_point_local = 1, n_ks_points_task
         do i_spin = 1, n_spin, 1
           lvl_tricoeff_occ_real(lbb_row:ubb_row,:,n_low_state_cpt2:n_homo_max,i_spin,i_k_point_local) = &
               real(lvl_tricoeff_mod_r(lbb_row:ubb_row,:,n_low_state_cpt2:n_homo_max,i_spin,i_k_point_local))
         end do
       end do
     end if

     allocate(lvl_tricoeff_unocc_real(lbb_row:ubb_row,&
              max_n_basis_sp,lpb_col:upb_col,n_spin,n_ks_points_task),stat=info)
     allocate(lvl_tricoeff_unocc_tmp(lbb_row:ubb_row,&
              max_n_basis_sp,lpb_col:upb_col),stat=info)
     
     dcsize=16
     if ((myid_kq_cpt2.lt.n_k_points) .and. &
         (size(lvl_tricoeff_mod_r) .ne. 0) .and. &
         (myid_col_cpt2 .eq. 0)) then
        nbytes=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*int(n_states,8)*&
            int(n_spin,8)*int(n_ks_points_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
     else
        nbytes=int(0,MPI_ADDRESS_KIND)
     end if

     call mpi_win_create(lvl_tricoeff_mod_r,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri,mpierr)
     call mpi_win_fence(0,win_tri,mpierr)

     count_1 = n_bb_row*max_n_basis_sp*n_pb_col
     count_2 = n_bb_row*max_n_basis_sp*n_states*n_spin*n_ks_points_task

     do i_k_point_local = 1, max_k_points_task, 1
       if (i_k_point_local .le. n_ks_points_task) then
         do i_spin = 1, n_spin, 1
           k_task=myid_kq_cpt2*n_tasks_bl_cpt2+myid_bl_cpt2-myid_col_cpt2
           offset = n_bb_row*max_n_basis_sp*n_states*n_spin*(i_k_point_local-1) + &
               n_bb_row*max_n_basis_sp*n_states*(i_spin-1) + &
               n_bb_row*max_n_basis_sp*(lpb_col+n_lumo_min-2)
           call mpi_get(lvl_tricoeff_unocc_tmp(:,:,:), &
               count_1, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_1, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
           call mpi_win_fence(0,win_tri,mpierr)
           lvl_tricoeff_unocc_real(:,:,:,i_spin,i_k_point_local) = &
               real(lvl_tricoeff_unocc_tmp(:,:,:))
         end do
       else
         do i_spin = 1, n_spin, 1
           call mpi_win_fence(0,win_tri,mpierr)
         end do
       end if
     end do

     call mpi_win_fence(0,win_tri,mpierr)

     call mpi_win_free(win_tri,mpierr)

     deallocate(lvl_tricoeff_mod_r)

     deallocate(lvl_tricoeff_unocc_tmp)

     ! ==========================================================================
     ! now re-distribute KS_eigenvectors in terms of (un)occupied orbitals for PT2
     ! ==========================================================================
     if (myid_col_cpt2.eq.0) then
       allocate(KS_eigenvectors_occ_real(n_basis,n_low_state_cpt2:n_homo_max,n_spin,n_ks_points_task),stat=info)
     end if
     allocate(KS_eigenvectors_unocc_real(n_basis,lpb_col:upb_col,n_spin,n_ks_points_task),stat=info)
     !  creat mpi_win for KS_eigenvectors

     dcsize=8
     nbytes=int(int(n_basis,8)*int(n_states,8)*int(n_spin,8)*&
         int(n_k_points_task,8)*int(dcsize,8), MPI_ADDRESS_KIND)

     call mpi_win_create(KS_eigenvector,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
     call mpi_win_fence(0,win_ev,mpierr)

     count_1=n_basis*n_states
     count_2=n_basis*n_occ_max
     count_3=n_basis*n_pb_col
     do i_k_point_local = 1, max_k_points_task, 1
       if (i_k_point_local .le. n_ks_points_task) then
         i_k_point = n_tasks_kq_cpt2*(i_k_point_local-1)+myid_kq_cpt2+1
         if (myid_col_cpt2.eq.0) then
           do i_spin = 1, n_spin, 1
             k_task=k_point_loc(1,i_k_point)
             offset=(k_point_loc(2,i_k_point)-1)*count_1*n_spin + &
                 count_1*(i_spin-1) + n_basis*(n_low_state_cpt2-1)
             call mpi_get(KS_eigenvectors_occ_real(:,:,i_spin,i_k_point_local), &
                 count_2, MPI_DOUBLE_PRECISION, k_task, &
                 offset, count_2, MPI_DOUBLE_PRECISION, win_ev, mpierr)
             call mpi_win_fence(0,win_ev,mpierr)
           end do
         else
           do i_spin = 1, n_spin, 1
             call mpi_win_fence(0,win_ev,mpierr)
           end do
         end if
         do i_spin = 1, n_spin, 1
           k_task=k_point_loc(1,i_k_point)
           offset=(k_point_loc(2,i_k_point)-1)*count_1*n_spin + &
               count_1*(i_spin-1) + n_basis*(lpb_col+n_lumo_min-2)
           call mpi_get(KS_eigenvectors_unocc_real(:,:,i_spin,i_k_point_local), &
               count_3, MPI_DOUBLE_PRECISION, k_task, &
              offset, count_3, MPI_DOUBLE_PRECISION, win_ev, mpierr)
           call mpi_win_fence(0,win_ev,mpierr)
         end do
       else
         do i_spin = 1, n_spin, 1
           call mpi_win_fence(0,win_ev,mpierr)
         end do
         do i_spin = 1, n_spin, 1
           call mpi_win_fence(0,win_ev,mpierr)
         end do
       end if
     end do

     call mpi_win_fence(0,win_ev,mpierr)

     call mpi_win_free(win_ev,mpierr)

     deallocate(KS_eigenvector)

     ! ==========================================================================
     ! now change the complex coulomb matrix to real one
     ! ==========================================================================

     allocate(coulomb_matr_blacs_real(lbb_row:ubb_row,lbb_col:ubb_col,n_kq_points_task),stat=info)
     call check_allocation(info, 'coulomb_matr_blacs_real            ')
     coulomb_matr_blacs_real(:,:,:) = real(coulomb_matr_blacs(:,:,:))

     deallocate(coulomb_matr_blacs)

   end subroutine redistribute_cpt2_blacs_real

   subroutine deallocate_cpt2_blacs()
       implicit none

       if (allocated(KS_eigenvectors_occ)) deallocate(KS_eigenvectors_occ)
       if (allocated(KS_eigenvectors_unocc)) deallocate(KS_eigenvectors_unocc)

       if (allocated(lvl_tricoeff_occ)) deallocate(lvl_tricoeff_occ)
       if (allocated(lvl_tricoeff_unocc)) deallocate(lvl_tricoeff_unocc)

       if (allocated(coulomb_matr_blacs)) deallocate(coulomb_matr_blacs)

       if (allocated(KS_eigenvectors_occ_real)) deallocate(KS_eigenvectors_occ_real)
       if (allocated(KS_eigenvectors_unocc_real)) deallocate(KS_eigenvectors_unocc_real)

       if (allocated(lvl_tricoeff_occ_real)) deallocate(lvl_tricoeff_occ_real)
       if (allocated(lvl_tricoeff_unocc_real)) deallocate(lvl_tricoeff_unocc_real)

       if (allocated(coulomb_matr_blacs_real)) deallocate(coulomb_matr_blacs_real)

   end subroutine deallocate_cpt2_blacs


end module redistribute_cpt2_kspace
