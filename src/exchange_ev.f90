
module exchange_ev
  use crpa_blacs
  implicit none

  logical:: real_ev

contains
  subroutine init_access_ev_real(n_k_points_special,n_k_points_special_task,eigenvectors,win_ev)
    integer:: n_k_points_special, n_k_points_special_task
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_special_task):: eigenvectors
    integer:: win_ev
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    real_ev=.true.

    dcsize=8
    
    if (mod(myid+n_tasks-1,n_tasks).lt.n_k_points_special) then
       nbytes=int(n_basis,MPI_ADDRESS_KIND)*n_states*n_spin*n_k_points_special_task*dcsize
    else
       nbytes=0
    end if
    


    call mpi_win_create(eigenvectors,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
    call mpi_win_fence(0,win_ev,mpierr)
    
  end subroutine init_access_ev_real
  
  subroutine init_access_ev_complex(n_k_points_special,n_k_points_special_task,eigenvectors_complex,win_ev)
    integer:: n_k_points_special, n_k_points_special_task
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_special_task):: eigenvectors_complex
    integer:: win_ev

    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    real_ev=.false.
    dcsize=16
    
    if (mod(myid+n_tasks-1,n_tasks).lt.n_k_points_special) then
       nbytes=int(n_basis,MPI_ADDRESS_KIND)*n_states*n_spin*n_k_points_special_task*dcsize
    else
       nbytes=0
    end if

    call mpi_win_create(eigenvectors_complex,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev,mpierr)
    call mpi_win_fence(0,win_ev,mpierr)
    
  end subroutine init_access_ev_complex

  subroutine finalize_access_ev(win_ev)
    integer:: win_ev
    integer:: mpierr

    call mpi_win_free(win_ev,mpierr)

  end subroutine finalize_access_ev
  
  subroutine access_ev(win_ev, k_point_loc, ev_curr)
    integer, dimension(2):: k_point_loc
    integer:: win_ev
    complex*16, dimension(n_basis,n_states,n_spin):: ev_curr

    real*8, dimension(:,:,:), allocatable :: ev_real
    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    call perfon('ex_ev')
    count=n_basis*n_states*n_spin
    if(myid_col.eq.0) then
       if(real_ev) allocate(ev_real(n_basis,n_states,n_spin))

       k_task=k_point_loc(1)
       offset=(k_point_loc(2)-1)*int(count,MPI_ADDRESS_KIND)
       if(real_ev) then
          call mpi_get(ev_real, count, MPI_DOUBLE_PRECISION, k_task, &
               offset, count, MPI_DOUBLE_PRECISION, win_ev, mpierr)
       else
          call mpi_get(ev_curr, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_ev, mpierr)
       end if
    end if

    call sync_ev(win_ev)
    if(real_ev.and.(myid_col.eq.0)) then
       ev_curr=ev_real
       deallocate(ev_real)
    end if

    if (n_tasks_bl.gt.1) call mpi_bcast(ev_curr, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col, mpierr)
    call perfoff

  end subroutine access_ev
  

  subroutine access_2_evs(win_ev, k_point_loc_1, ev_curr_1, k_point_loc_2, ev_curr_2)

    integer, dimension(2):: k_point_loc_1, k_point_loc_2
    integer:: win_ev
    complex*16, dimension(n_basis,n_states,n_spin):: ev_curr_1, ev_curr_2

    real*8, dimension(:,:,:), allocatable :: ev_real_1, ev_real_2
    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count, k_task
    integer:: mpierr

    call perfon('ex_ev')
    count=n_basis*n_states*n_spin

    if(myid_col.eq.0) then
       if(real_ev) allocate(ev_real_1(n_basis,n_states,n_spin),ev_real_2(n_basis,n_states,n_spin))

       k_task=k_point_loc_1(1)
       offset=(k_point_loc_1(2)-1)*int(count,MPI_ADDRESS_KIND)
       if(real_ev) then
          call mpi_get(ev_real_1, count, MPI_DOUBLE_PRECISION, k_task, &
               offset, count, MPI_DOUBLE_PRECISION, win_ev, mpierr)
       else
          call mpi_get(ev_curr_1, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_ev, mpierr)
       end if

       k_task=k_point_loc_2(1)
       offset=(k_point_loc_2(2)-1)*int(count,MPI_ADDRESS_KIND)
       if(real_ev) then
          call mpi_get(ev_real_2, count, MPI_DOUBLE_PRECISION, k_task, &
               offset, count, MPI_DOUBLE_PRECISION, win_ev, mpierr)
       else
          call mpi_get(ev_curr_2, count, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count, MPI_DOUBLE_COMPLEX, win_ev, mpierr)
       end if
    end if

    call sync_ev(win_ev)
    if(real_ev.and.(myid_col.eq.0)) then
       ev_curr_1=ev_real_1    
       ev_curr_2=ev_real_2    
       deallocate(ev_real_1,ev_real_2)    
    end if

    if (n_tasks_bl.gt.1) then
       call mpi_bcast(ev_curr_1, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs, mpierr)
       call mpi_bcast(ev_curr_2, count, MPI_DOUBLE_COMPLEX, 0, comm_blacs, mpierr)
    end if
    call perfoff

  end subroutine access_2_evs
  

  subroutine sync_ev(win_ev)
    integer:: win_ev
    integer:: mpierr

    call mpi_win_fence(0,win_ev,mpierr)

  end subroutine sync_ev

end module exchange_ev
