
module exchange_trico_cpt2
  use cpt2_blacs
  use basis, only: max_n_basis_sp
  implicit none

  complex(kind=8), dimension(:,:,:,:,:),pointer:: lvl_tricoeff_occ_p => null()
  complex(kind=8), dimension(:,:,:,:,:),pointer:: lvl_tricoeff_unocc_p => null()

  real(kind=8), dimension(:,:,:,:,:),pointer:: lvl_tricoeff_occ_p_real => null()
  real(kind=8), dimension(:,:,:,:,:),pointer:: lvl_tricoeff_unocc_p_real => null()

contains
  subroutine sinit_access_trico_cpt2(n_k_points_special,n_k_points_special_task,&
                                     lvl_tricoeff_occ,win_tri_occ, &
                                     lvl_tricoeff_unocc,win_tri_unocc)
    integer:: n_k_points_special, n_k_points_special_task
    complex(kind=8), dimension(lbb_row:,1:,n_low_state_cpt2:,1:,1:), target:: lvl_tricoeff_occ
    complex(kind=8), dimension(lbb_row:,1:,lpb_col:,1:,1:), target:: lvl_tricoeff_unocc
    integer:: win_tri_occ, win_tri_unocc
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    !avoid copies if no kq parallelization
    if(n_tasks_kq_cpt2.eq.1) then
      if (myid_col_cpt2 .eq. 0) lvl_tricoeff_occ_p => lvl_tricoeff_occ
      lvl_tricoeff_unocc_p => lvl_tricoeff_unocc
    else
       dcsize=16
       !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for lvl_tricoeff_mod_r_k > 2 GiB)
       if (myid_kq_cpt2.lt.n_k_points_special .and. myid_col_cpt2.eq.0) then
          nbytes=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*int(n_occ_max,8)*&
              int(n_spin,8)*int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       else
          nbytes=int(0,MPI_ADDRESS_KIND)
       end if
       call mpi_win_create(lvl_tricoeff_occ,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri_occ,mpierr)
       call mpi_win_fence(0,win_tri_occ,mpierr)


       nbytes=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*int(n_pb_col,8)*&
           int(n_spin,8)*int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       call mpi_win_create(lvl_tricoeff_unocc,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri_unocc,mpierr)
       call mpi_win_fence(0,win_tri_unocc,mpierr)

    end if

  end subroutine sinit_access_trico_cpt2

  subroutine sfinalize_access_trico_cpt2(win_tri_occ,win_tri_unocc)
    integer:: win_tri_occ, win_tri_unocc
    integer:: mpierr
    
    if(n_tasks_kq_cpt2.gt.1) then
        call mpi_win_free(win_tri_occ,mpierr)
        call mpi_win_free(win_tri_unocc,mpierr)
    end if

  end subroutine sfinalize_access_trico_cpt2
  
  subroutine saccess_4_tricoeffs_cpt2(spins, win_tri_occ, win_tri_unocc,   &
                               lvl_tricoeff_curr_k, lvl_tricoeff_curr_q,   &
                               lvl_tricoeff_curr_kp, lvl_tricoeff_curr_qp, &
                               kq_pair_new,kq_pair_old)
    integer:: spins,win_tri_occ, win_tri_unocc
    complex(kind=8), dimension(:,:,:,:), pointer:: lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp
    complex(kind=8), dimension(:,:,:,:), pointer:: lvl_tricoeff_curr_q, lvl_tricoeff_curr_qp
    integer,dimension(6):: kq_pair_new,kq_pair_old

    integer(kind=MPI_ADDRESS_KIND):: offset, count_1, count_2
    integer:: k_task
    integer:: mpierr

    if ((kq_pair_new(1) .eq. 0) .or. (kq_pair_new(2) .eq. 0) .or. &
        (kq_pair_new(3) .eq. 0) .or. (kq_pair_new(4) .eq. 0)) then
      call ssync_trico_cpt2(win_tri_occ, win_tri_unocc)
      return
    end if

    count_1=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*&
        int(n_occ_max,8)*int(spins,8),MPI_ADDRESS_KIND)
    count_2=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*&
        int(n_pb_col,8)*int(spins,8),MPI_ADDRESS_KIND)

    call perfon('ex_tri')
    if(n_tasks_kq_cpt2.gt.1) then
      if (myid_col_cpt2.eq.0) then
        ! for lvl_tricoeff_curr_k
        if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
          if (kq_pair_new(1) .eq. kq_pair_old(3)) then
            lvl_tricoeff_curr_k = lvl_tricoeff_curr_kp
          else
            k_task=mod(kq_pair_new(1)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
            offset=int(int(((kq_pair_new(1)-1)/n_tasks_kq_cpt2),8)*count_1,MPI_ADDRESS_KIND)
            call mpi_get(lvl_tricoeff_curr_k, count_1, MPI_DOUBLE_COMPLEX, k_task, &
                 offset, count_1, MPI_DOUBLE_COMPLEX, win_tri_occ, mpierr)
          end if
        end if
        ! for lvl_tricoeff_curr_kp
        if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
          k_task=mod(kq_pair_new(3)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(3)-1)/n_tasks_kq_cpt2),8)*count_1,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_kp, count_1, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_1, MPI_DOUBLE_COMPLEX, win_tri_occ, mpierr)
        end if
      end if
      call ssync_trico_cpt2_single(win_tri_occ)
      ! for lvl_tricoeff_curr_q
      if (.not. (kq_pair_new(2) .eq. kq_pair_old(2))) then
        if (kq_pair_new(2) .eq. kq_pair_old(4)) then
          lvl_tricoeff_curr_q = lvl_tricoeff_curr_qp
        else
          k_task=mod(kq_pair_new(2)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(2)-1)/n_tasks_kq_cpt2),8)*count_2,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_q, count_2, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_2, MPI_DOUBLE_COMPLEX, win_tri_unocc, mpierr)
        end if
      end if
      ! for lvl_tricoeff_curr_qp
      if (.not. (kq_pair_new(4) .eq. kq_pair_old(4))) then
          k_task=mod(kq_pair_new(4)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(4)-1)/n_tasks_kq_cpt2),8)*count_2,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_qp, count_2, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_2, MPI_DOUBLE_COMPLEX, win_tri_unocc, mpierr)
      end if
      call ssync_trico_cpt2_single(win_tri_unocc)
    else
      if (myid_col_cpt2.eq.0) then
        lvl_tricoeff_curr_k  => lvl_tricoeff_occ_p(:,:,:,:,kq_pair_new(1))
        lvl_tricoeff_curr_kp => lvl_tricoeff_occ_p(:,:,:,:,kq_pair_new(3))
      end if
      lvl_tricoeff_curr_q  => lvl_tricoeff_unocc_p(:,:,:,:,kq_pair_new(2))
      lvl_tricoeff_curr_qp => lvl_tricoeff_unocc_p(:,:,:,:,kq_pair_new(4))
    end if

    !call ssync_trico_cpt2(win_tri_occ,win_tri_unocc)

    if (n_tasks_bl_cpt2.gt.1) then
      call mpi_bcast(lvl_tricoeff_curr_k, size(lvl_tricoeff_curr_k),   &
                     MPI_DOUBLE_COMPLEX,0,comm_blacs_col_cpt2,mpierr)
      call mpi_bcast(lvl_tricoeff_curr_kp, size(lvl_tricoeff_curr_kp), &
                     MPI_DOUBLE_COMPLEX,0,comm_blacs_col_cpt2,mpierr)
    end if
  
    call perfoff

  end subroutine saccess_4_tricoeffs_cpt2

  subroutine ssync_trico_cpt2_single(win_tri)
    integer:: win_tri
    integer:: mpierr

    if(n_tasks_kq_cpt2.gt.1) then
       call perfon('ex_tri')      
       call mpi_win_fence(0,win_tri,mpierr)
       call perfoff
    end if
  end subroutine ssync_trico_cpt2_single

  subroutine ssync_trico_cpt2(win_tri_occ, win_tri_unocc)
    integer:: win_tri_occ, win_tri_unocc
    integer:: mpierr

    if(n_tasks_kq_cpt2.gt.1) then
       call perfon('ex_tri')      
       call mpi_win_fence(0,win_tri_occ,mpierr)
       call mpi_win_fence(0,win_tri_unocc,mpierr)
       call perfoff
    end if
  end subroutine ssync_trico_cpt2

  subroutine saccess_4_tricoeffs_cpt2_complex(spins, win_tri_occ, win_tri_unocc,   &
                               lvl_tricoeff_curr_k, lvl_tricoeff_curr_q,   &
                               lvl_tricoeff_curr_kp, lvl_tricoeff_curr_qp, &
                               kq_pair_new,kq_pair_old,curr_run)
    integer:: spins,win_tri_occ, win_tri_unocc
    complex(kind=8), dimension(:,:,:,:), &
        pointer:: lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp
    complex(kind=8), dimension(:,:,:,:), &
        pointer:: lvl_tricoeff_curr_q, lvl_tricoeff_curr_qp
    integer,dimension(6):: kq_pair_new,kq_pair_old
    logical :: curr_run

    integer(kind=MPI_ADDRESS_KIND):: offset, count_1 ,count_2
    integer:: k_task    
    integer:: mpierr

    if (.not. curr_run)  return

    count_1=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*&
        int(n_occ_max,8)*int(spins,8),MPI_ADDRESS_KIND)
    count_2=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*&
        int(n_pb_col,8)*int(spins,8),MPI_ADDRESS_KIND)

    call perfon('ex_tri')
    if(n_tasks_kq_cpt2.gt.1) then
      if (myid_col_cpt2.eq.0) then
        ! for lvl_tricoeff_curr_k
        if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
          if (kq_pair_new(1) .eq. kq_pair_old(3)) then
            lvl_tricoeff_curr_k = lvl_tricoeff_curr_kp
          else
            k_task=mod(kq_pair_new(1)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
            offset=int(int(((kq_pair_new(1)-1)/n_tasks_kq_cpt2),8)*count_1,MPI_ADDRESS_KIND)
            call mpi_get(lvl_tricoeff_curr_k, count_1, MPI_DOUBLE_COMPLEX, k_task, &
                 offset, count_1, MPI_DOUBLE_COMPLEX, win_tri_occ, mpierr)
          end if
        end if
        ! for lvl_tricoeff_curr_kp
        if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
          k_task=mod(kq_pair_new(3)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(3)-1)/n_tasks_kq_cpt2),8)*count_1,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_kp, count_1, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_1, MPI_DOUBLE_COMPLEX, win_tri_occ, mpierr)
        end if
      end if
      ! for lvl_tricoeff_curr_q
      if (.not. (kq_pair_new(2) .eq. kq_pair_old(2))) then
        if (kq_pair_new(2) .eq. kq_pair_old(4)) then
          lvl_tricoeff_curr_q = lvl_tricoeff_curr_qp
        else
          k_task=mod(kq_pair_new(2)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(2)-1)/n_tasks_kq_cpt2),8)*count_2,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_q, count_2, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_2, MPI_DOUBLE_COMPLEX, win_tri_unocc, mpierr)
        end if
      end if
      ! for lvl_tricoeff_curr_qp
      if (.not. (kq_pair_new(4) .eq. kq_pair_old(4))) then
        k_task=mod(kq_pair_new(4)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
        offset=int(int(((kq_pair_new(4)-1)/n_tasks_kq_cpt2),8)*count_2,MPI_ADDRESS_KIND)
        call mpi_get(lvl_tricoeff_curr_qp, count_2, MPI_DOUBLE_COMPLEX, k_task, &
             offset, count_2, MPI_DOUBLE_COMPLEX, win_tri_unocc, mpierr)
      end if
    else
      if (myid_col_cpt2.eq.0) then
        lvl_tricoeff_curr_k  => lvl_tricoeff_occ_p(:,:,:,:,kq_pair_new(1))
        lvl_tricoeff_curr_kp => lvl_tricoeff_occ_p(:,:,:,:,kq_pair_new(3))
      end if
      lvl_tricoeff_curr_q  => lvl_tricoeff_unocc_p(:,:,:,:,kq_pair_new(2))
      lvl_tricoeff_curr_qp => lvl_tricoeff_unocc_p(:,:,:,:,kq_pair_new(4))
    end if

    call perfoff

  end subroutine saccess_4_tricoeffs_cpt2_complex

  subroutine ssync_trico_cpt2_complex(win_tri_occ, win_tri_unocc,   &
                               lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp, &
                               kq_pair_new,kq_pair_old,&
                               curr_run)
    integer:: win_tri_occ, win_tri_unocc
    complex(kind=8), dimension(:,:,:,:), &
        pointer:: lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp
    integer,dimension(6):: kq_pair_new,kq_pair_old
    logical :: curr_run
    integer:: mpierr

    !if ((kq_pair_new(1) .eq. 0) .or. (kq_pair_new(2) .eq. 0) .or. &
    !    (kq_pair_new(3) .eq. 0) .or. (kq_pair_new(4) .eq. 0)) then
    if (.not. curr_run) then
      call ssync_trico_cpt2(win_tri_occ, win_tri_unocc)
      return
    end if

    if(n_tasks_kq_cpt2.gt.1) then
       call perfon('ex_tri')      
       call mpi_win_fence(0,win_tri_occ,mpierr)
       call mpi_win_fence(0,win_tri_unocc,mpierr)
       call perfoff
    end if
    if (n_tasks_bl_cpt2.gt.1) then
      if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
        call mpi_bcast(lvl_tricoeff_curr_k, size(lvl_tricoeff_curr_k),   &
                       MPI_DOUBLE_COMPLEX,0,comm_blacs_col_cpt2,mpierr)
      end if
      if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
        call mpi_bcast(lvl_tricoeff_curr_kp, size(lvl_tricoeff_curr_kp), &
                       MPI_DOUBLE_COMPLEX,0,comm_blacs_col_cpt2,mpierr)
      end if
    end if
  end subroutine ssync_trico_cpt2_complex

  subroutine sinit_access_trico_real(n_k_points_special,n_k_points_special_task,&
                                     lvl_tricoeff_occ_real,win_tri_occ, &
                                     lvl_tricoeff_unocc_real,win_tri_unocc)
    integer:: n_k_points_special, n_k_points_special_task
    real(kind=8), dimension(lbb_row:,1:,n_low_state_cpt2:,1:,1:), &
        target:: lvl_tricoeff_occ_real
    real(kind=8), dimension(lbb_row:,1:,lpb_col:,1:,1:), &
        target:: lvl_tricoeff_unocc_real
    integer:: win_tri_occ, win_tri_unocc
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    !avoid copies if no kq parallelization
    if(n_tasks_kq_cpt2.eq.1) then
      if (myid_col_cpt2 .eq. 0) lvl_tricoeff_occ_p_real => lvl_tricoeff_occ_real
      lvl_tricoeff_unocc_p_real => lvl_tricoeff_unocc_real
    else
       dcsize=8
       !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for lvl_tricoeff_mod_r_k > 2 GiB)
       if (myid_kq_cpt2.lt.n_k_points_special .and. myid_col_cpt2.eq.0) then
          nbytes=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*int(n_occ_max,8)*&
              int(n_spin,8)*int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       else
          nbytes=int(0,MPI_ADDRESS_KIND)
       end if
       call mpi_win_create(lvl_tricoeff_occ_real,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri_occ,mpierr)
       call mpi_win_fence(0,win_tri_occ,mpierr)

       nbytes=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*int(n_pb_col,8)*&
           int(n_spin,8)*int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       call mpi_win_create(lvl_tricoeff_unocc_real,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri_unocc,mpierr)
       call mpi_win_fence(0,win_tri_unocc,mpierr)

    end if

  end subroutine sinit_access_trico_real

  subroutine saccess_4_tricoeffs_cpt2_real(spins, win_tri_occ, win_tri_unocc,   &
                               lvl_tricoeff_curr_k, lvl_tricoeff_curr_q,   &
                               lvl_tricoeff_curr_kp, lvl_tricoeff_curr_qp, &
                               kq_pair_new,kq_pair_old,curr_run)
    integer:: spins,win_tri_occ, win_tri_unocc
    real(kind=8), dimension(:,:,:,:), &
        pointer:: lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp
    real(kind=8), dimension(:,:,:,:), &
        pointer:: lvl_tricoeff_curr_q, lvl_tricoeff_curr_qp
    integer,dimension(6):: kq_pair_new,kq_pair_old
    logical :: curr_run

    integer(kind=MPI_ADDRESS_KIND):: offset, count_1 ,count_2
    integer:: k_task
    integer:: mpierr

    !if ((kq_pair_new(1) .eq. 0) .or. (kq_pair_new(2) .eq. 0) .or. &
    !    (kq_pair_new(3) .eq. 0) .or. (kq_pair_new(4) .eq. 0)) then
      !call ssync_trico_cpt2(win_tri_occ, win_tri_unocc)
    !end if
    if (.not. curr_run)  return

    count_1=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*&
        int(n_occ_max,8)*int(spins,8),MPI_ADDRESS_KIND)
    count_2=int(int(n_bb_row,8)*int(max_n_basis_sp,8)*&
        int(n_pb_col,8)*int(spins,8),MPI_ADDRESS_KIND)

    call perfon('ex_tri')
    if(n_tasks_kq_cpt2.gt.1) then
      if (myid_col_cpt2.eq.0) then
        ! for lvl_tricoeff_curr_k
        if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
          if (kq_pair_new(1) .eq. kq_pair_old(3)) then
            lvl_tricoeff_curr_k = lvl_tricoeff_curr_kp
          else
            k_task=mod(kq_pair_new(1)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
            offset=int(int(((kq_pair_new(1)-1)/n_tasks_kq_cpt2),8)*count_1,MPI_ADDRESS_KIND)
            call mpi_get(lvl_tricoeff_curr_k, count_1, MPI_DOUBLE_PRECISION, k_task, &
                 offset, count_1, MPI_DOUBLE_PRECISION, win_tri_occ, mpierr)
          end if
        end if
        ! for lvl_tricoeff_curr_kp
        if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
          k_task=mod(kq_pair_new(3)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(3)-1)/n_tasks_kq_cpt2),8)*count_1,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_kp, count_1, MPI_DOUBLE_PRECISION, k_task, &
               offset, count_1, MPI_DOUBLE_PRECISION, win_tri_occ, mpierr)
        end if
      end if
      ! for lvl_tricoeff_curr_q
      if (.not. (kq_pair_new(2) .eq. kq_pair_old(2))) then
        if (kq_pair_new(2) .eq. kq_pair_old(4)) then
          lvl_tricoeff_curr_q = lvl_tricoeff_curr_qp
        else
          k_task=mod(kq_pair_new(2)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(2)-1)/n_tasks_kq_cpt2),8)*count_2,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_q, count_2, MPI_DOUBLE_PRECISION, k_task, &
               offset, count_2, MPI_DOUBLE_PRECISION, win_tri_unocc, mpierr)
        end if
      end if
      ! for lvl_tricoeff_curr_qp
      if (.not. (kq_pair_new(4) .eq. kq_pair_old(4))) then
          k_task=mod(kq_pair_new(4)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(4)-1)/n_tasks_kq_cpt2),8)*count_2,MPI_ADDRESS_KIND)
          call mpi_get(lvl_tricoeff_curr_qp, count_2, MPI_DOUBLE_PRECISION, k_task, &
               offset, count_2, MPI_DOUBLE_PRECISION, win_tri_unocc, mpierr)
      end if
    else
      if (myid_col_cpt2.eq.0) then
        lvl_tricoeff_curr_k  => lvl_tricoeff_occ_p_real(:,:,:,:,kq_pair_new(1))
        lvl_tricoeff_curr_kp => lvl_tricoeff_occ_p_real(:,:,:,:,kq_pair_new(3))
      end if
      lvl_tricoeff_curr_q  => lvl_tricoeff_unocc_p_real(:,:,:,:,kq_pair_new(2))
      lvl_tricoeff_curr_qp => lvl_tricoeff_unocc_p_real(:,:,:,:,kq_pair_new(4))
    end if

    call perfoff

  end subroutine saccess_4_tricoeffs_cpt2_real

  subroutine ssync_trico_cpt2_real(win_tri_occ, win_tri_unocc,   &
                               lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp, &
                               kq_pair_new,kq_pair_old,&
                               curr_run)
    integer:: win_tri_occ, win_tri_unocc
    real(kind=8), dimension(:,:,:,:), &
        pointer:: lvl_tricoeff_curr_k, lvl_tricoeff_curr_kp
    integer,dimension(6):: kq_pair_new,kq_pair_old
    logical :: curr_run
    integer:: mpierr

    if (.not. curr_run) then
      call ssync_trico_cpt2(win_tri_occ, win_tri_unocc)
      return
    end if

    if(n_tasks_kq_cpt2.gt.1) then
       call perfon('ex_tri')      
       call mpi_win_fence(0,win_tri_occ,mpierr)
       call mpi_win_fence(0,win_tri_unocc,mpierr)
       call perfoff
    end if
    if (n_tasks_bl_cpt2.gt.1) then
      if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
        call mpi_bcast(lvl_tricoeff_curr_k, size(lvl_tricoeff_curr_k),   &
                       MPI_DOUBLE_PRECISION,0,comm_blacs_col_cpt2,mpierr)
      end if
      if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
        call mpi_bcast(lvl_tricoeff_curr_kp, size(lvl_tricoeff_curr_kp), &
                       MPI_DOUBLE_PRECISION,0,comm_blacs_col_cpt2,mpierr)
      end if
    end if
  end subroutine ssync_trico_cpt2_real

  subroutine ssync_trico_cpt2_gamma_only(lvl_tricoeff_curr_k,curr_run)
    real(kind=8), dimension(:,:,:,:), pointer:: lvl_tricoeff_curr_k
    logical :: curr_run
    integer:: mpierr

    if (.not. curr_run) return

    if (n_tasks_bl_cpt2.gt.1) then
       call mpi_bcast(lvl_tricoeff_curr_k, size(lvl_tricoeff_curr_k),   &
                      MPI_DOUBLE_PRECISION,0,comm_blacs_col_cpt2,mpierr)
    end if
  end subroutine ssync_trico_cpt2_gamma_only

end module exchange_trico_cpt2
