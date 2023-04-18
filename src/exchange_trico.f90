
module exchange_trico
  use crpa_blacs
  use basis, only: max_n_basis_sp
  implicit none
  complex*16, dimension(:,:,:,:,:),pointer:: lvl_tricoeff_arr
  
contains
  subroutine init_access_trico(n_k_points_special,n_k_points_special_task,lvl_tricoeff_mod_r_k,win_tri)
    integer:: n_k_points_special, n_k_points_special_task
    complex*16, dimension(lbb_row:,1:,1:,1:,1:):: lvl_tricoeff_mod_r_k
    integer:: win_tri
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr


    dcsize=16
    !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for lvl_tricoeff_mod_r_k > 2 GiB)
    if (myid_irkq.lt.n_k_points_special .and. (size(lvl_tricoeff_mod_r_k) .ne. 0)) then
       nbytes=int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*n_spin*n_k_points_special_task*dcsize
    else
       nbytes=0
    end if
    
    call mpi_win_create(lvl_tricoeff_mod_r_k,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri,mpierr)
    call mpi_win_fence(0,win_tri,mpierr)

  end subroutine init_access_trico

  subroutine finalize_access_trico(win_tri)
    integer:: win_tri
    integer:: mpierr

    call mpi_win_free(win_tri,mpierr)

  end subroutine finalize_access_trico

  subroutine sinit_access_trico(n_k_points_special,n_k_points_special_task,lvl_tricoeff_mod_r_k,win_tri)
    integer:: n_k_points_special, n_k_points_special_task
    complex*16, dimension(lbb_row:,1:,1:,1:,1:), target:: lvl_tricoeff_mod_r_k
    integer:: win_tri
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    !avoid copies if no irkq parallelization
    if(n_tasks_irkq.eq.1) then
       if(myid_col.eq.0) lvl_tricoeff_arr => lvl_tricoeff_mod_r_k
    else
       dcsize=16
       !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for lvl_tricoeff_mod_r_k > 2 GiB)
       if (myid_irkq.lt.n_k_points_special .and. (size(lvl_tricoeff_mod_r_k) .ne. 0)) then
          nbytes=int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*n_spin*n_k_points_special_task*dcsize
       else
          nbytes=0
       end if
       
       call mpi_win_create(lvl_tricoeff_mod_r_k,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_tri,mpierr)
       call mpi_win_fence(0,win_tri,mpierr)
    end if

  end subroutine sinit_access_trico

  subroutine sfinalize_access_trico(win_tri)
    integer:: win_tri
    integer:: mpierr
    
    if(n_tasks_irkq.gt.1) call mpi_win_free(win_tri,mpierr)

  end subroutine sfinalize_access_trico
  
  subroutine access_tricoeff(spins, win_tri, i_k_point, lvl_tricoeff_curr_r, lvl_tricoeff_curr_c)
    integer:: i_k_point, spins
    integer:: win_tri
    complex*16, dimension(n_bb_row,max_n_basis_sp,n_states,spins):: lvl_tricoeff_curr_r
    complex*16, dimension(n_bb_col,max_n_basis_sp,n_states,spins):: lvl_tricoeff_curr_c

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    call perfon('ex_tri')
    if(myid_col.eq.0) then
       k_task=mod(i_k_point-1,n_tasks_irkq)*n_tasks_bl+myid_bl
       offset=((i_k_point-1)/n_tasks_irkq)*int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
       count=n_bb_row*max_n_basis_sp*n_states*spins
       call mpi_get(lvl_tricoeff_curr_r, count, MPI_DOUBLE_COMPLEX, k_task, &
            offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
    end if

    if (n_tasks_bl.gt.1) then
       if(myid_row.eq.0) then
          k_task=mod(i_k_point-1,n_tasks_irkq)*n_tasks_bl+myid_row + myid_col*n_tasks_row
          offset=((i_k_point-1)/n_tasks_irkq)*int(n_bb_col,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
          count=n_bb_col*max_n_basis_sp*n_states*spins
          call mpi_get(lvl_tricoeff_curr_c, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
       end if
    end if

    call sync_trico(win_tri)

    if (n_tasks_bl.gt.1) then
       count=n_bb_row*max_n_basis_sp*n_states*spins
       call mpi_bcast(lvl_tricoeff_curr_r, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col, mpierr)

       count=n_bb_col*max_n_basis_sp*n_states*spins
       call mpi_bcast(lvl_tricoeff_curr_c, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_row, mpierr)
    end if
    call perfoff

  end subroutine access_tricoeff

  subroutine saccess_tricoeff(spins, win_tri, i_k_point, lvl_tricoeff_curr_r)
    integer:: i_k_point, spins
    integer:: win_tri
    complex*16, dimension(:,:,:,:), pointer:: lvl_tricoeff_curr_r

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    if(n_tasks_irkq.gt.1) then
       call perfon('ex_tri')
       if(myid_col.eq.0) then
          k_task=mod(i_k_point-1,n_tasks_irkq)*n_tasks_bl+myid_bl
          offset=((i_k_point-1)/n_tasks_irkq)*int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
          count=n_bb_row*max_n_basis_sp*n_states*spins
          call mpi_get(lvl_tricoeff_curr_r, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
       end if

       call sync_trico(win_tri)

       call perfoff
    else
       if (myid_col.eq.0) lvl_tricoeff_curr_r => lvl_tricoeff_arr(:,:,:,:,i_k_point)
    end if
  end subroutine saccess_tricoeff
  
  subroutine scopy_tricoeff(lvl_tricoeff_curr_in, lvl_tricoeff_curr_out)
    complex*16, dimension(:,:,:,:), pointer:: lvl_tricoeff_curr_in, lvl_tricoeff_curr_out

    if (myid_col.eq.0) then
       if(n_tasks_irkq.gt.1) then
          lvl_tricoeff_curr_out=lvl_tricoeff_curr_in
       else
          lvl_tricoeff_curr_out => lvl_tricoeff_curr_in
       end if
    end if

  end subroutine scopy_tricoeff
  

  subroutine access_2_tricoeffs(spins, win_tri, i_k_point_1, lvl_tricoeff_curr_r_1, lvl_tricoeff_curr_c_1, &
       i_k_point_2, lvl_tricoeff_curr_r_2, lvl_tricoeff_curr_c_2)
    integer:: i_k_point_1, i_k_point_2, spins
    integer:: win_tri
    complex*16, dimension(n_bb_row,max_n_basis_sp,n_states,spins):: lvl_tricoeff_curr_r_1, lvl_tricoeff_curr_r_2
    complex*16, dimension(n_bb_col,max_n_basis_sp,n_states,spins):: lvl_tricoeff_curr_c_1, lvl_tricoeff_curr_c_2

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    call perfon('ex_tri')
    if(myid_col.eq.0) then
       count=n_bb_row*max_n_basis_sp*n_states*spins

       k_task=mod(i_k_point_1-1,n_tasks_irkq)*n_tasks_bl+myid_bl
       offset=((i_k_point_1-1)/n_tasks_irkq)*int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
       call mpi_get(lvl_tricoeff_curr_r_1, count, MPI_DOUBLE_COMPLEX, k_task, &
            offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
       
       k_task=mod(i_k_point_2-1,n_tasks_irkq)*n_tasks_bl+myid_bl
       offset=((i_k_point_2-1)/n_tasks_irkq)*int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
       call mpi_get(lvl_tricoeff_curr_r_2, count, MPI_DOUBLE_COMPLEX, k_task, &
            offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
    end if
    if (n_tasks_bl.gt.1) then
       if(myid_row.eq.0) then
          count=n_bb_col*max_n_basis_sp*n_states*spins

          k_task=mod(i_k_point_1-1,n_tasks_irkq)*n_tasks_bl + myid_row + myid_col*n_tasks_row
          offset=((i_k_point_1-1)/n_tasks_irkq)*int(n_bb_col,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
          call mpi_get(lvl_tricoeff_curr_c_1, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
          
          k_task=mod(i_k_point_2-1,n_tasks_irkq)*n_tasks_bl + myid_row + myid_col*n_tasks_row
          offset=((i_k_point_2-1)/n_tasks_irkq)*int(n_bb_col,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
          call mpi_get(lvl_tricoeff_curr_c_2, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
       end if
    end if

    call sync_trico(win_tri)
    
    if (n_tasks_bl.gt.1) then
       count=n_bb_row*max_n_basis_sp*n_states*spins
       call mpi_bcast(lvl_tricoeff_curr_r_1, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col, mpierr)
       call mpi_bcast(lvl_tricoeff_curr_r_2, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col, mpierr)

       count=n_bb_col*max_n_basis_sp*n_states*spins
       call mpi_bcast(lvl_tricoeff_curr_c_1, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_row, mpierr)
       call mpi_bcast(lvl_tricoeff_curr_c_2, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_row, mpierr)
    end if
    call perfoff

  end subroutine access_2_tricoeffs

  subroutine saccess_2_tricoeffs(spins, win_tri, i_k_point_1, lvl_tricoeff_curr_r_1, &
       i_k_point_2, lvl_tricoeff_curr_r_2)
    integer:: i_k_point_1, i_k_point_2, spins
    integer:: win_tri
    complex*16, dimension(:,:,:,:), pointer:: lvl_tricoeff_curr_r_1, lvl_tricoeff_curr_r_2

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    if(n_tasks_irkq.gt.1) then
       call perfon('ex_tri')
       if(myid_col.eq.0) then
          count=n_bb_row*max_n_basis_sp*n_states*spins
          
          k_task=mod(i_k_point_1-1,n_tasks_irkq)*n_tasks_bl+myid_bl
          offset=((i_k_point_1-1)/n_tasks_irkq)*int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
          call mpi_get(lvl_tricoeff_curr_r_1, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
          
          k_task=mod(i_k_point_2-1,n_tasks_irkq)*n_tasks_bl+myid_bl
          offset=((i_k_point_2-1)/n_tasks_irkq)*int(n_bb_row,MPI_ADDRESS_KIND)*max_n_basis_sp*n_states*spins
          call mpi_get(lvl_tricoeff_curr_r_2, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_tri, mpierr)
       end if
       
       call sync_trico(win_tri)
       call perfoff
    else
       if (myid_col.eq.0) then
       lvl_tricoeff_curr_r_1 => lvl_tricoeff_arr(:,:,:,:,i_k_point_1)
       lvl_tricoeff_curr_r_2 => lvl_tricoeff_arr(:,:,:,:,i_k_point_2)
    end if
    end if
  end subroutine saccess_2_tricoeffs

  subroutine sync_trico(win_tri)
    integer:: win_tri
    integer:: mpierr

    call perfon('ex_tri')      
    call mpi_win_fence(0,win_tri,mpierr)
    call perfoff

  end subroutine sync_trico

  subroutine ssync_trico(win_tri)
    integer:: win_tri
    integer:: mpierr

    if(n_tasks_irkq.gt.1) then
       call perfon('ex_tri')      
       call mpi_win_fence(0,win_tri,mpierr)
       call perfoff
    end if
  end subroutine ssync_trico
end module exchange_trico
