
module exchange_coulomb_cpt2
  use cpt2_blacs
  use basis, only: max_n_basis_sp
  implicit none
  complex(kind=8), dimension(:,:,:),pointer:: coulomb_matr_arr_cpt2 => null()
  real(kind=8), dimension(:,:,:),pointer:: coulomb_matr_arr_real => null()

contains
  subroutine sinit_access_coulomb_cpt2(n_k_points_special,n_k_points_special_task,coulomb_matr,win_coul)
    integer:: n_k_points_special, n_k_points_special_task
    complex(kind=8), dimension(lbb_row:,lbb_col:,1:), target:: coulomb_matr
    integer:: win_coul
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    !avoid copies if no kq parallelization
    if(n_tasks_kq_cpt2.eq.1) then
       coulomb_matr_arr_cpt2 => coulomb_matr
    else
       dcsize=16
       !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for coulomb_matr > 2 GiB)
       if (myid_kq_cpt2.lt.n_k_points_special .and. (size(coulomb_matr) .ne. 0)) then
          nbytes=int(int(n_bb_row,8)*int(n_bb_col,8)*&
              int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       else
          nbytes=int(0,MPI_ADDRESS_KIND)
       end if
       
       call mpi_win_create(coulomb_matr,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_coul,mpierr)
       call mpi_win_fence(0,win_coul,mpierr)
    end if

  end subroutine sinit_access_coulomb_cpt2

  subroutine sfinalize_access_coulomb_cpt2(win_coul)
    integer:: win_coul
    integer:: mpierr
    
    if(n_tasks_kq_cpt2.gt.1) call mpi_win_free(win_coul,mpierr)

  end subroutine sfinalize_access_coulomb_cpt2
  
  subroutine saccess_coulomb_cpt2(spins, win_coul, i_k_point, coulomb_matr_curr_r)
    integer:: i_k_point, spins
    integer:: win_coul
    complex(kind=8), dimension(:,:), pointer:: coulomb_matr_curr_r

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    if(n_tasks_kq_cpt2.gt.1) then
       call perfon('ex_coul')
       k_task=mod(i_k_point-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
       offset=int(int(((i_k_point-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*int(n_bb_col,8),&
           MPI_ADDRESS_KIND)
       count=n_bb_row*n_bb_col
       call mpi_get(coulomb_matr_curr_r, count, MPI_DOUBLE_COMPLEX, k_task, &
            offset, count, MPI_DOUBLE_COMPLEX, win_coul, mpierr)

       call sync_coulomb_cpt2(win_coul)

       call perfoff
    else
       coulomb_matr_curr_r => coulomb_matr_arr_cpt2(:,:,i_k_point)
    end if
  end subroutine saccess_coulomb_cpt2

  subroutine saccess_2_coulomb_cpt2(spins, win_coul, &
          i_k_point_1, i_k_point_1_old,coulomb_matr_curr_r_1, &
          i_k_point_2, i_k_point_2_old,coulomb_matr_curr_r_2 &
          )
    integer:: spins, i_k_point_1, i_k_point_1_old, i_k_point_2, i_k_point_2_old
    integer:: win_coul
    complex(kind=8), dimension(:,:), pointer:: coulomb_matr_curr_r_1,coulomb_matr_curr_r_2

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    if (i_k_point_1.eq.0 .or. i_k_point_2.eq.0) then
       call ssync_coulomb_cpt2(win_coul)
       return
    end if

    count=n_bb_row*n_bb_col

    if(n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_coul')
      ! for coulomb_matr_curr_r_1
      if (.not. (i_k_point_1.eq.i_k_point_1_old)) then
        if (i_k_point_1.eq.i_k_point_2_old) then
          coulomb_matr_curr_r_1=coulomb_matr_curr_r_2
        else
          k_task=mod(i_k_point_1-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((i_k_point_1-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*&
              int(n_bb_col,8), MPI_ADDRESS_KIND)
          call mpi_get(coulomb_matr_curr_r_1, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_coul, mpierr)
        end if
      end if
      ! for coulomb_matr_curr_r_2
      if (.not. (i_k_point_2.eq.i_k_point_2_old)) then
        k_task=mod(i_k_point_2-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
        offset=int(int(((i_k_point_2-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*&
            int(n_bb_col,8),MPI_ADDRESS_KIND)
        call mpi_get(coulomb_matr_curr_r_2, count, MPI_DOUBLE_COMPLEX, k_task, &
             offset, count, MPI_DOUBLE_COMPLEX, win_coul, mpierr)
      end if

      call sync_coulomb_cpt2(win_coul)

      call perfoff
    else
       coulomb_matr_curr_r_1 => coulomb_matr_arr_cpt2(:,:,i_k_point_1)
       coulomb_matr_curr_r_2 => coulomb_matr_arr_cpt2(:,:,i_k_point_2)
    end if
  end subroutine saccess_2_coulomb_cpt2

  subroutine sync_coulomb_cpt2(win_coul)
    integer:: win_coul
    integer:: mpierr

    call perfon('ex_coul')      
    call mpi_win_fence(0,win_coul,mpierr)
    call perfoff

  end subroutine sync_coulomb_cpt2
  
  subroutine ssync_coulomb_cpt2(win_coul)
    integer:: win_coul
    integer:: mpierr

    if(n_tasks_kq_cpt2.gt.1) then
       call perfon('ex_coul')      
       call mpi_win_fence(0,win_coul,mpierr)
       call perfoff
    end if
  end subroutine ssync_coulomb_cpt2

  subroutine saccess_2_coulomb_cpt2_complex(spins, win_coul, &
          i_k_point_1, i_k_point_1_old,coulomb_matr_curr_r_1, &
          i_k_point_2, i_k_point_2_old,coulomb_matr_curr_r_2 &
          )
    integer:: spins, i_k_point_1, i_k_point_1_old, i_k_point_2, i_k_point_2_old
    integer:: win_coul
    complex(kind=8), dimension(:,:), pointer:: coulomb_matr_curr_r_1,coulomb_matr_curr_r_2

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    if (i_k_point_1.eq.0 .or. i_k_point_2.eq.0) then
       !call ssync_coulomb_cpt2(win_coul)
       return
    end if

    count=n_bb_row*n_bb_col

    if(n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_coul')
      ! for coulomb_matr_curr_r_1
      if (.not. (i_k_point_1.eq.i_k_point_1_old)) then
        if (i_k_point_1.eq.i_k_point_2_old) then
          coulomb_matr_curr_r_1=coulomb_matr_curr_r_2
        else
          k_task=mod(i_k_point_1-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((i_k_point_1-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*&
              int(n_bb_col,8),MPI_ADDRESS_KIND)
          call mpi_get(coulomb_matr_curr_r_1, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_coul, mpierr)
        end if
      end if
      ! for coulomb_matr_curr_r_2
      if (.not. (i_k_point_2.eq.i_k_point_2_old)) then
        k_task=mod(i_k_point_2-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
        offset=int(int(((i_k_point_2-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*&
            int(n_bb_col,8),MPI_ADDRESS_KIND)
        call mpi_get(coulomb_matr_curr_r_2, count, MPI_DOUBLE_COMPLEX, k_task, &
             offset, count, MPI_DOUBLE_COMPLEX, win_coul, mpierr)
      end if

      call perfoff
    else
       coulomb_matr_curr_r_1 => coulomb_matr_arr_cpt2(:,:,i_k_point_1)
       coulomb_matr_curr_r_2 => coulomb_matr_arr_cpt2(:,:,i_k_point_2)
    end if
  end subroutine saccess_2_coulomb_cpt2_complex

  subroutine sinit_access_coulomb_real(n_k_points_special,n_k_points_special_task,coulomb_matr,win_coul)
    integer:: n_k_points_special, n_k_points_special_task
    real(kind=8), dimension(lbb_row:,lbb_col:,1:), target:: coulomb_matr
    integer:: win_coul
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    !avoid copies if no kq parallelization
    if(n_tasks_kq_cpt2.eq.1) then
       coulomb_matr_arr_real => coulomb_matr
    else
       dcsize=8
       !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for coulomb_matr > 2 GiB)
       if (myid_kq_cpt2.lt.n_k_points_special .and. (size(coulomb_matr) .ne. 0)) then
          nbytes=int(int(n_bb_row,8)*int(n_bb_col,8)*&
              int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
       else
          nbytes=int(0,MPI_ADDRESS_KIND)
       end if
       
       call mpi_win_create(coulomb_matr,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_coul,mpierr)
       call mpi_win_fence(0,win_coul,mpierr)
    end if

  end subroutine sinit_access_coulomb_real

  subroutine saccess_2_coulomb_cpt2_real(spins, win_coul, &
          i_k_point_1, i_k_point_1_old,coulomb_matr_curr_r_1, &
          i_k_point_2, i_k_point_2_old,coulomb_matr_curr_r_2 &
          )
    integer:: spins, i_k_point_1, i_k_point_1_old, i_k_point_2, i_k_point_2_old
    integer:: win_coul
    real(kind=8), dimension(:,:), pointer:: coulomb_matr_curr_r_1,coulomb_matr_curr_r_2

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    if (i_k_point_1.eq.0 .or. i_k_point_2.eq.0) then
       !call ssync_coulomb_cpt2(win_coul)
       return
    end if

    count=n_bb_row*n_bb_col

    if(n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_coul')
      ! for coulomb_matr_curr_r_1
      if (.not. (i_k_point_1.eq.i_k_point_1_old)) then
        if (i_k_point_1.eq.i_k_point_2_old) then
          coulomb_matr_curr_r_1=coulomb_matr_curr_r_2
        else
          k_task=mod(i_k_point_1-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((i_k_point_1-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*&
              int(n_bb_col,8),MPI_ADDRESS_KIND)
          call mpi_get(coulomb_matr_curr_r_1, count, MPI_DOUBLE_PRECISION, k_task, &
               offset, count, MPI_DOUBLE_PRECISION, win_coul, mpierr)
        end if
      end if
      ! for coulomb_matr_curr_r_2
      if (.not. (i_k_point_2.eq.i_k_point_2_old)) then
        k_task=mod(i_k_point_2-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
        offset=int(int(((i_k_point_2-1)/n_tasks_kq_cpt2),8)*int(n_bb_row,8)*&
            int(n_bb_col,8),MPI_ADDRESS_KIND)
        call mpi_get(coulomb_matr_curr_r_2, count, MPI_DOUBLE_PRECISION, k_task, &
             offset, count, MPI_DOUBLE_PRECISION, win_coul, mpierr)
      end if

      call perfoff
    else
       coulomb_matr_curr_r_1 => coulomb_matr_arr_real(:,:,i_k_point_1)
       coulomb_matr_curr_r_2 => coulomb_matr_arr_real(:,:,i_k_point_2)
    end if
  end subroutine saccess_2_coulomb_cpt2_real

end module exchange_coulomb_cpt2
