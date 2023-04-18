
module exchange_ev_cpt2
  use cpt2_blacs
  use hartree_fock, only : n_homo_max

  implicit none

  complex(kind=8), dimension(:,:,:,:),pointer:: eigenvectors_occ_p => null()
  complex(kind=8), dimension(:,:,:,:),pointer:: eigenvectors_unocc_p => null()

  real(kind=8), dimension(:,:,:,:),pointer:: eigenvectors_occ_p_real => null()
  real(kind=8), dimension(:,:,:,:),pointer:: eigenvectors_unocc_p_real => null()

!  logical:: real_ev_cpt2

contains
  
  subroutine sinit_access_ev_cpt2(n_k_points_special,n_k_points_special_task, &
                                  eigenvectors_occ,win_ev_occ, &
                                  eigenvectors_unocc,win_ev_unocc)
    integer:: n_k_points_special, n_k_points_special_task
    complex(kind=8), dimension(1:,n_low_state_cpt2:,1:,1:), target:: eigenvectors_occ
    complex(kind=8), dimension(1:,lpb_col:,1:,1:), target:: eigenvectors_unocc
    integer:: win_ev_occ, win_ev_unocc
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    if (n_tasks_kq_cpt2.eq.1) then
      if (myid_col_cpt2 .eq. 0) eigenvectors_occ_p => eigenvectors_occ
      eigenvectors_unocc_p => eigenvectors_unocc
    else
      dcsize=16
      !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for lvl_tricoeff_mod_r_k > 2 GiB)
      if (myid_kq_cpt2.lt.n_k_points_special .and. myid_col_cpt2.eq.0) then
         nbytes=int(int(n_basis,8)*int(n_occ_max,8)*int(n_spin,8)*&
             int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
      else
         nbytes=int(0,MPI_ADDRESS_KIND)
      end if
      call mpi_win_create(eigenvectors_occ,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev_occ,mpierr)
      call mpi_win_fence(0,win_ev_occ,mpierr)

      nbytes=n_basis*n_pb_col*n_spin*n_k_points_special_task*dcsize
      call mpi_win_create(eigenvectors_unocc,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev_unocc,mpierr)
      call mpi_win_fence(0,win_ev_unocc,mpierr)
    end if
    
  end subroutine sinit_access_ev_cpt2

  subroutine sfinalize_access_ev_cpt2(win_ev_occ,win_ev_unocc)
    integer:: win_ev_occ, win_ev_unocc
    integer:: mpierr

    if (n_tasks_kq_cpt2.gt.1) then
      call mpi_win_free(win_ev_occ,mpierr)
      call mpi_win_free(win_ev_unocc,mpierr)
    end if

  end subroutine sfinalize_access_ev_cpt2
  
  subroutine saccess_4_evs_cpt2(spins,win_ev_occ,win_ev_unocc, &
                               ev_curr_k, ev_curr_q, ev_curr_kp, ev_curr_qp, &
                               kq_pair_new, kq_pair_old &
                              )

    integer:: spins, win_ev_occ,win_ev_unocc
    complex(kind=8), dimension(:,:,:), pointer:: ev_curr_k, ev_curr_kp
    complex(kind=8), dimension(:,:,:), pointer:: ev_curr_q, ev_curr_qp
    integer, dimension(6):: kq_pair_old, kq_pair_new

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count_1, count_2, k_task
    integer:: mpierr

    call perfon('ex_ev')
    count_1=n_basis*n_occ_max*n_spin
    count_2=n_basis*n_pb_col*n_spin

    if ((kq_pair_new(1) .eq. 0) .or. (kq_pair_new(2) .eq. 0) .or. &
        (kq_pair_new(3) .eq. 0) .or. (kq_pair_new(4) .eq. 0)) then
      call ssync_ev_cpt2(win_ev_occ,win_ev_unocc)
      return
    end if

    if (n_tasks_kq_cpt2.gt.1) then
      if(myid_col_cpt2.eq.0) then
        ! for ev_curr_k
        if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
          if (kq_pair_new(1) .eq. kq_pair_old(3)) then
            ev_curr_k = ev_curr_kp
          else
            k_task=mod(kq_pair_new(1)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
            offset=int(int(((kq_pair_new(1)-1)/n_tasks_kq_cpt2),8)*int(count_1,8),MPI_ADDRESS_KIND)
            call mpi_get(ev_curr_k, count_1, MPI_DOUBLE_COMPLEX, k_task, &
                 offset, count_1, MPI_DOUBLE_COMPLEX, win_ev_occ, mpierr)
          end if
        end if
        ! for ev_curr_kp
        if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
          k_task=mod(kq_pair_new(3)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(3)-1)/n_tasks_kq_cpt2),8)*int(count_1,8),MPI_ADDRESS_KIND)
          call mpi_get(ev_curr_kp, count_1, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_1, MPI_DOUBLE_COMPLEX, win_ev_occ, mpierr)
        end if
      end if
      call ssync_ev_cpt2_single(win_ev_occ)
      ! for ev_curr_q
      if (.not. (kq_pair_new(2) .eq. kq_pair_old(2))) then
        if (kq_pair_new(2) .eq. kq_pair_old(4)) then
          ev_curr_q = ev_curr_qp
        else
          k_task=mod(kq_pair_new(2)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(2)-1)/n_tasks_kq_cpt2),8)*int(count_2,8),MPI_ADDRESS_KIND)
          call mpi_get(ev_curr_q, count_2, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_2, MPI_DOUBLE_COMPLEX, win_ev_unocc, mpierr)
        end if
      end if
      ! for ev_curr_4
      if (.not. (kq_pair_new(4) .eq. kq_pair_old(4))) then
        k_task=mod(kq_pair_new(4)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(4)-1)/n_tasks_kq_cpt2),8)*int(count_2,8),MPI_ADDRESS_KIND)
        call mpi_get(ev_curr_qp, count_2, MPI_DOUBLE_COMPLEX, k_task, &
             offset, count_2, MPI_DOUBLE_COMPLEX, win_ev_unocc, mpierr)
      end if

      call ssync_ev_cpt2_single(win_ev_unocc)

    else
      if (myid_col_cpt2.eq.0) then
        ev_curr_k  => eigenvectors_occ_p(:,:,:,kq_pair_new(1))
        ev_curr_kp => eigenvectors_occ_p(:,:,:,kq_pair_new(3))
      end if
      ev_curr_q  => eigenvectors_unocc_p(:,:,:,kq_pair_new(2))
      ev_curr_qp => eigenvectors_unocc_p(:,:,:,kq_pair_new(4))
    end if

    if (n_tasks_bl_cpt2.gt.1) then
       call mpi_bcast(ev_curr_k,  count_1, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col_cpt2, mpierr)
       call mpi_bcast(ev_curr_kp, count_1, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col_cpt2, mpierr)
    end if
    call perfoff

  end subroutine saccess_4_evs_cpt2
  

  subroutine ssync_ev_cpt2(win_ev_occ,win_ev_unocc)
    integer:: win_ev_occ, win_ev_unocc
    integer:: mpierr

    if (n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_ev')
      call mpi_win_fence(0,win_ev_occ,mpierr)
      call mpi_win_fence(0,win_ev_unocc,mpierr)
      call perfoff
    end if

  end subroutine ssync_ev_cpt2

  subroutine ssync_ev_cpt2_single(win_ev)
    integer:: win_ev
    integer:: mpierr

    if (n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_ev')
      call mpi_win_fence(0,win_ev,mpierr)
      call perfoff
    end if
  end subroutine ssync_ev_cpt2_single

  subroutine saccess_4_evs_cpt2_complex(spins,win_ev_occ,win_ev_unocc, &
                               ev_curr_k, ev_curr_q, ev_curr_kp, ev_curr_qp, &
                               kq_pair_new, kq_pair_old, curr_run &
                              )

    integer:: spins, win_ev_occ,win_ev_unocc
    complex(kind=8), dimension(:,:,:), pointer:: ev_curr_k, ev_curr_kp
    complex(kind=8), dimension(:,:,:), pointer:: ev_curr_q, ev_curr_qp
    integer, dimension(6):: kq_pair_old, kq_pair_new
    logical :: curr_run

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count_1, count_2, k_task
    integer:: mpierr

    call perfon('ex_ev')
    count_1=n_basis*n_occ_max*n_spin
    count_2=n_basis*n_pb_col*n_spin

    !if ((kq_pair_new(1) .eq. 0) .or. (kq_pair_new(2) .eq. 0) .or. &
    !    (kq_pair_new(3) .eq. 0) .or. (kq_pair_new(4) .eq. 0)) then
    !  call ssync_ev_cpt2(win_ev_occ,win_ev_unocc)
    !  return
    !end if
    if (.not. curr_run) return

    if (n_tasks_kq_cpt2.gt.1) then
      if(myid_col_cpt2.eq.0) then
        ! for ev_curr_k
        if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
          if (kq_pair_new(1) .eq. kq_pair_old(3)) then
            ev_curr_k = ev_curr_kp
          else
            k_task=mod(kq_pair_new(1)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
            offset=int(int(((kq_pair_new(1)-1)/n_tasks_kq_cpt2),8)*int(count_1,8),MPI_ADDRESS_KIND)
            call mpi_get(ev_curr_k, count_1, MPI_DOUBLE_COMPLEX, k_task, &
                 offset, count_1, MPI_DOUBLE_COMPLEX, win_ev_occ, mpierr)
          end if
        end if
        ! for ev_curr_kp
        if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
          k_task=mod(kq_pair_new(3)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(3)-1)/n_tasks_kq_cpt2),8)*int(count_1,8),MPI_ADDRESS_KIND)
          call mpi_get(ev_curr_kp, count_1, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_1, MPI_DOUBLE_COMPLEX, win_ev_occ, mpierr)
        end if
      end if
      !call ssync_ev_cpt2_single(win_ev_occ)
      ! for ev_curr_q
      if (.not. (kq_pair_new(2) .eq. kq_pair_old(2))) then
        if (kq_pair_new(2) .eq. kq_pair_old(4)) then
          ev_curr_q = ev_curr_qp
        else
          k_task=mod(kq_pair_new(2)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(2)-1)/n_tasks_kq_cpt2),8)*int(count_2,8),MPI_ADDRESS_KIND)
          call mpi_get(ev_curr_q, count_2, MPI_DOUBLE_COMPLEX, k_task, &
               offset, count_2, MPI_DOUBLE_COMPLEX, win_ev_unocc, mpierr)
        end if
      end if
      ! for ev_curr_4
      if (.not. (kq_pair_new(4) .eq. kq_pair_old(4))) then
        k_task=mod(kq_pair_new(4)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
        offset=int(int(((kq_pair_new(4)-1)/n_tasks_kq_cpt2),8)*int(count_2,8),MPI_ADDRESS_KIND)
        call mpi_get(ev_curr_qp, count_2, MPI_DOUBLE_COMPLEX, k_task, &
             offset, count_2, MPI_DOUBLE_COMPLEX, win_ev_unocc, mpierr)
      end if

      !call ssync_ev_cpt2_single(win_ev_unocc)

    else
      if (myid_col_cpt2.eq.0) then
        ev_curr_k  => eigenvectors_occ_p(:,:,:,kq_pair_new(1))
        ev_curr_kp => eigenvectors_occ_p(:,:,:,kq_pair_new(3))
      end if
      ev_curr_q  => eigenvectors_unocc_p(:,:,:,kq_pair_new(2))
      ev_curr_qp => eigenvectors_unocc_p(:,:,:,kq_pair_new(4))
    end if

    !if (n_tasks_bl_cpt2.gt.1) then
    !   call mpi_bcast(ev_curr_k,  count_1, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col_cpt2, mpierr)
    !   call mpi_bcast(ev_curr_kp, count_1, MPI_DOUBLE_COMPLEX, 0, comm_blacs_col_cpt2, mpierr)
    !end if
    call perfoff

  end subroutine saccess_4_evs_cpt2_complex

  subroutine ssync_ev_cpt2_complex(win_ev_occ,win_ev_unocc,&
                               ev_curr_k, ev_curr_kp, &
                               kq_pair_new, kq_pair_old, &
                               curr_run)
    integer:: win_ev_occ, win_ev_unocc
    complex(kind=8), dimension(:,:,:), pointer:: ev_curr_k, ev_curr_kp
    integer, dimension(6):: kq_pair_old, kq_pair_new
    integer:: mpierr
    logical :: curr_run

    if (.not. curr_run) then
      call ssync_ev_cpt2(win_ev_occ,win_ev_unocc)
      return
    end if

    if (n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_ev')
      call mpi_win_fence(0,win_ev_occ,mpierr)
      call mpi_win_fence(0,win_ev_unocc,mpierr)
      call perfoff
    end if

    if (n_tasks_bl_cpt2.gt.1) then
      if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
        call mpi_bcast(ev_curr_k,  size(ev_curr_k),  MPI_DOUBLE_COMPLEX, &
            0, comm_blacs_col_cpt2, mpierr)
      end if
      if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
        call mpi_bcast(ev_curr_kp, size(ev_curr_kp), MPI_DOUBLE_COMPLEX, &
            0, comm_blacs_col_cpt2, mpierr)
      end if
    end if

  end subroutine ssync_ev_cpt2_complex

  subroutine sinit_access_ev_real(n_k_points_special,n_k_points_special_task, &
                                  eigenvectors_occ_real,win_ev_occ, &
                                  eigenvectors_unocc_real,win_ev_unocc)
    integer:: n_k_points_special, n_k_points_special_task
    real(kind=8), dimension(1:,n_low_state_cpt2:,1:,1:), target:: eigenvectors_occ_real
    real(kind=8), dimension(1:,lpb_col:,1:,1:), target:: eigenvectors_unocc_real
    integer:: win_ev_occ, win_ev_unocc
    
    integer(kind=MPI_ADDRESS_KIND):: nbytes
    integer:: dcsize
    integer:: mpierr

    if (n_tasks_kq_cpt2.eq.1) then
      if (myid_col_cpt2 .eq. 0) eigenvectors_occ_p_real => eigenvectors_occ_real
      eigenvectors_unocc_p_real => eigenvectors_unocc_real
    else
      dcsize=8
      !MPI_ADDRESS_KIND is just an integer for some MPIs and may be too small (for lvl_tricoeff_mod_r_k > 2 GiB)
      if (myid_kq_cpt2.lt.n_k_points_special .and. myid_col_cpt2.eq.0) then
         nbytes=int(int(n_basis,8)*int(n_occ_max,8)*int(n_spin,8)*&
             int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
      else
         nbytes=int(0,MPI_ADDRESS_KIND)
      end if
      call mpi_win_create(eigenvectors_occ_real,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev_occ,mpierr)
      call mpi_win_fence(0,win_ev_occ,mpierr)

      nbytes=int(int(n_basis,8)*int(n_pb_col,8)*int(n_spin,8)*&
          int(n_k_points_special_task,8)*int(dcsize,8),MPI_ADDRESS_KIND)
      call mpi_win_create(eigenvectors_unocc_real,nbytes,dcsize,MPI_INFO_NULL,mpi_comm_world,win_ev_unocc,mpierr)
      call mpi_win_fence(0,win_ev_unocc,mpierr)
    end if
    
  end subroutine sinit_access_ev_real

  subroutine saccess_4_evs_cpt2_real(spins,win_ev_occ,win_ev_unocc, &
                               ev_curr_k, ev_curr_q, ev_curr_kp, ev_curr_qp, &
                               kq_pair_new, kq_pair_old, curr_run &
                              )

    integer:: spins, win_ev_occ,win_ev_unocc
    real(kind=8), dimension(:,:,:), pointer:: ev_curr_k, ev_curr_kp
    real(kind=8), dimension(:,:,:), pointer:: ev_curr_q, ev_curr_qp
    integer, dimension(6):: kq_pair_old, kq_pair_new
    logical :: curr_run

    integer(kind=MPI_ADDRESS_KIND):: offset
    integer:: count_1, count_2, k_task
    integer:: mpierr

    call perfon('ex_ev')
    count_1=n_basis*n_occ_max*n_spin
    count_2=n_basis*n_pb_col*n_spin

    if (.not. curr_run) return

    if (n_tasks_kq_cpt2.gt.1) then
      if(myid_col_cpt2.eq.0) then
        ! for ev_curr_k
        if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
          if (kq_pair_new(1) .eq. kq_pair_old(3)) then
            ev_curr_k = ev_curr_kp
          else
            k_task=mod(kq_pair_new(1)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
            offset=int(int(((kq_pair_new(1)-1)/n_tasks_kq_cpt2),8)*int(count_1,8),MPI_ADDRESS_KIND)
            call mpi_get(ev_curr_k, count_1, MPI_DOUBLE_PRECISION, k_task, &
                 offset, count_1, MPI_DOUBLE_PRECISION, win_ev_occ, mpierr)
          end if
        end if
        ! for ev_curr_kp
        if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
          k_task=mod(kq_pair_new(3)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(3)-1)/n_tasks_kq_cpt2),8)*int(count_1,8),MPI_ADDRESS_KIND)
          call mpi_get(ev_curr_kp, count_1, MPI_DOUBLE_PRECISION, k_task, &
               offset, count_1, MPI_DOUBLE_PRECISION, win_ev_occ, mpierr)
        end if
      end if
      ! for ev_curr_q
      if (.not. (kq_pair_new(2) .eq. kq_pair_old(2))) then
        if (kq_pair_new(2) .eq. kq_pair_old(4)) then
          ev_curr_q = ev_curr_qp
        else
          k_task=mod(kq_pair_new(2)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
          offset=int(int(((kq_pair_new(2)-1)/n_tasks_kq_cpt2),8)*int(count_2,8),MPI_ADDRESS_KIND)
          call mpi_get(ev_curr_q, count_2, MPI_DOUBLE_PRECISION, k_task, &
               offset, count_2, MPI_DOUBLE_PRECISION, win_ev_unocc, mpierr)
        end if
      end if
      ! for ev_curr_4
      if (.not. (kq_pair_new(4) .eq. kq_pair_old(4))) then
        k_task=mod(kq_pair_new(4)-1,n_tasks_kq_cpt2)*n_tasks_bl_cpt2+myid_bl_cpt2
        offset=int(int(((kq_pair_new(4)-1)/n_tasks_kq_cpt2),8)*int(count_2,8),MPI_ADDRESS_KIND)
        call mpi_get(ev_curr_qp, count_2, MPI_DOUBLE_PRECISION, k_task, &
             offset, count_2, MPI_DOUBLE_PRECISION, win_ev_unocc, mpierr)
      end if

    else
      if (myid_col_cpt2.eq.0) then
        ev_curr_k  => eigenvectors_occ_p_real(:,:,:,kq_pair_new(1))
        ev_curr_kp => eigenvectors_occ_p_real(:,:,:,kq_pair_new(3))
      end if
      ev_curr_q  => eigenvectors_unocc_p_real(:,:,:,kq_pair_new(2))
      ev_curr_qp => eigenvectors_unocc_p_real(:,:,:,kq_pair_new(4))
    end if

    call perfoff

  end subroutine saccess_4_evs_cpt2_real

  subroutine ssync_ev_cpt2_real(win_ev_occ,win_ev_unocc,&
                               ev_curr_k, ev_curr_kp, &
                               kq_pair_new, kq_pair_old, &
                               curr_run)
    integer:: win_ev_occ, win_ev_unocc
    real(kind=8), dimension(:,:,:), pointer:: ev_curr_k, ev_curr_kp
    integer, dimension(6):: kq_pair_old, kq_pair_new
    integer:: mpierr
    logical :: curr_run

    if (.not. curr_run) then
      call ssync_ev_cpt2(win_ev_occ,win_ev_unocc)
      return
    end if

    if (n_tasks_kq_cpt2.gt.1) then
      call perfon('ex_ev')
      call mpi_win_fence(0,win_ev_occ,mpierr)
      call mpi_win_fence(0,win_ev_unocc,mpierr)
      call perfoff
    end if

    if (n_tasks_bl_cpt2.gt.1) then
       if (.not. (kq_pair_new(1) .eq. kq_pair_old(1))) then
         call mpi_bcast(ev_curr_k,  size(ev_curr_k),  &
             MPI_DOUBLE_PRECISION, 0, comm_blacs_col_cpt2, mpierr)
       end if
       if (.not. (kq_pair_new(3) .eq. kq_pair_old(3))) then
         call mpi_bcast(ev_curr_kp, size(ev_curr_kp), &
             MPI_DOUBLE_PRECISION, 0, comm_blacs_col_cpt2, mpierr)
       end if
    end if

  end subroutine ssync_ev_cpt2_real

  subroutine ssync_ev_cpt2_gamma_only(ev_curr_k, curr_run)
    real(kind=8), dimension(:,:,:), pointer:: ev_curr_k
    integer:: mpierr
    logical :: curr_run

    if (.not. curr_run) then
      return
    end if

    if (n_tasks_bl_cpt2.gt.1) then
       call mpi_bcast(ev_curr_k,  size(ev_curr_k),  &
           MPI_DOUBLE_PRECISION, 0, comm_blacs_col_cpt2, mpierr)
    end if
  end subroutine ssync_ev_cpt2_gamma_only

end module exchange_ev_cpt2
