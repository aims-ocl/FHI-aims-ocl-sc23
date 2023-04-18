!****s* FHI-aims/fciqmc/evaluate_w_vector_1a
!  NAME
!    evaluate_w_vector_1a
!  SYNOPSIS 

      subroutine evaluate_w_vector_1a()

!  PURPOSE
!
!  The '''evaluate_w_vector_1a''' evaluate the contributions of single-excitation
!  in the alpha spin space for w_vect
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Jan Kloppenburg
!  SOURCE
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use synchronize_mpi
  use fciqmc_module

    implicit none
    ! temp variables
    integer, dimension(:,:), allocatable :: occ_num_1
    ! temp indices
    integer :: i_state, j_state, a_state, b_state, c_state, k_state
    integer :: i_ci, i_ci_local, j_ci, k_ci
    integer :: errnum
    !character*128 :: info_str
    real*8  :: w_1, w_2, w_3, w_4, w_5
    real*8  :: V_s0,state_sign
    real*8  :: ddot
    integer :: tmp_spin(2)
    integer :: substitude(2,2), occupy(2,2)
    integer :: index_start(4)
    integer :: i_task, j_task
    integer :: mpierr, info

    real*8, allocatable, dimension(:) :: c_vect_tmp

    allocate(c_vect_tmp(n_configuration_table(1)),stat=errnum)
    call check_allocation(errnum, 'c_vect_tmp in CI')

    tmp_spin = 1
    if (n_spin .eq. 2) tmp_spin(2) = 2
    
    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 in CI')
    endif

    !write(use_unit,'(2X,A3,6(A16))') 'Ind', 'V_s0', 'w_1', 'w_2', 'w_3', 'w_4', 'w_5'

    ! ==============================================================
    ! now construct matrix elements between one-electron substutions
    ! ==============================================================
    ! First for the substitution in the alpha space
    ! --------------------------------------------------------------

    do i_task = 0, n_tasks-1, 1
      ! broadcast c_vect in i_task to all processors
      !write(use_unit,*) "igor debug",i_task, myid
      if (i_task .eq. myid) then
          c_vect_tmp(1:n_configuration_table(i_task+1)) = &
              c_vect(:)
      endif
      call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
      call mpi_bcast(c_vect_tmp,n_configuration_table(i_task+1), &
                     MPI_REAL8, i_task, mpi_comm_global, mpierr)

      i_ci        = memory_table(1,1,1,myid+1)-1
      i_ci_local  = memory_table(1,2,1,myid+1)-1
      index_start = index_table(:,1,myid+1)
      do i_state = index_start(1), n_elec(1), 1
        do a_state = index_start(2), n_states, 1

          i_ci       = i_ci + 1

          if (i_ci .gt. memory_table(2,1,1,myid+1)) goto 666

          i_ci_local = i_ci_local + 1

          ! ===================================================
          ! By definition, V_s0 is zero for singlet exitations
          ! ===================================================
          if (i_task .eq. 0) then
          !    w_vect(i_ci_local) = V_s0 * c_0
              w_vect(i_ci_local) = 0.0d0
          endif
          ! For HF orbitals, w_1 and w_2 are always zero
          w_1 = 0.0d0
          w_2 = 0.0d0
          w_3 = 0.0d0
          ! the first term in w_3
          do j_state = i_start_ci, n_elec(1), 1
            do b_state = n_elec(1)+1, n_states, 1
              substitude(1,1) = j_state
              occupy(1,1)     = b_state

              call seek_index_v03_1(substitude,occupy,j_ci,j_task)
              if (j_task .eq. i_task) then
                  w_3 = w_3 + (&
                      ddot(n_basbas,&
                      ovlp_3ks(:,j_state,i_state,tmp_spin(1)),1,&
                      ovlp_3ks(:,a_state,b_state,tmp_spin(1)),1)&
                      - &
                      ddot(n_basbas,&
                      ovlp_3ks(:,j_state,b_state,tmp_spin(1)),1,&
                      ovlp_3ks(:,a_state,i_state,tmp_spin(1)),1)&
                      ) * c_vect_tmp(j_ci)
              endif

            enddo
          enddo
          ! the second term in w_3
          do j_state = i_start_ci, n_elec(2), 1
            do b_state = n_elec(2)+1, n_states, 1
              substitude(1,2) = j_state
              occupy(1,2)     = b_state

              call seek_index_v03_2(substitude,occupy,j_ci,j_task)
              if (j_task .eq. i_task) then
                  w_3 = w_3 - (&
                      ddot(n_basbas,&
                      ovlp_3ks(:,j_state,b_state,tmp_spin(2)),1,&
                      ovlp_3ks(:,a_state,i_state,tmp_spin(1)),1)&
                      ) * c_vect_tmp(j_ci)
              endif

            enddo
          enddo
          w_3 = - w_3
          w_4 = 0.0d0
          ! the first term in w_4 j, b, c>b
          do j_state = i_start_ci, n_elec(1), 1 
            if (j_state.ne.i_state) then
              do b_state = n_elec(1)+1, n_states, 1
                do c_state = b_state+1, n_states, 1
                  substitude(1,1) = i_state
                  substitude(2,1) = j_state
                  occupy(1,1)     = b_state
                  occupy(2,1)     = c_state

                  call ci_coefficient_sign(substitude,occupy,state_sign,1)
                  call seek_index_v03_4(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_4 = w_4 - (&
                          ddot(n_basbas,&
                          ovlp_3ks(:,j_state,b_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,a_state,c_state,tmp_spin(1)),1)&
                          - &
                          ddot(n_basbas,&
                          ovlp_3ks(:,j_state,c_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,a_state,b_state,tmp_spin(1)),1)&
                          ) * c_vect_tmp(j_ci) * state_sign
                  endif

                enddo
              enddo
            end if
          enddo ! the first term in w_4

          ! the second term in w_4
          ! j in O_{\beta}, b in V_{\alpha}, c in V_{\beta}
          do j_state = i_start_ci, n_elec(2), 1 
            do b_state = n_elec(1)+1, n_states, 1
              do c_state = n_elec(2)+1, n_states, 1
                substitude(1,1) = i_state
                substitude(1,2) = j_state
                occupy(1,1)     = b_state
                occupy(1,2)     = c_state

                call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                if (j_task .eq. i_task) then
                    w_4 = w_4 + (&
                        ddot(n_basbas,&
                        ovlp_3ks(:,j_state,c_state,tmp_spin(2)),1,&
                        ovlp_3ks(:,a_state,b_state,tmp_spin(1)),1)&
                        ) * c_vect_tmp(j_ci)
                endif

              enddo
            enddo
          enddo ! the third term in w_4
          w_5 = 0.0d0
          ! there are two terms in w_5
          ! the first term
          !! 2) j, k>j (j,k in O_{\alpha}),  b (b in V_{\alpha})
          do j_state = i_start_ci, n_elec(1), 1 
            do k_state = j_state+1, n_elec(1), 1
              do b_state = n_elec(1)+1, n_states, 1
                if (b_state.ne.a_state) then
                  substitude(1,1) = j_state
                  substitude(2,1) = k_state
                  occupy(1,1)     = a_state
                  occupy(2,1)     = b_state

                  call ci_coefficient_sign(substitude,occupy,state_sign,1)
                  call seek_index_v03_4(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_5 = w_5 - (&
                          ddot(n_basbas,&
                          ovlp_3ks(:,j_state,i_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,k_state,b_state,tmp_spin(1)),1)&
                          - &
                          ddot(n_basbas,&
                          ovlp_3ks(:,j_state,b_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,k_state,i_state,tmp_spin(1)),1)&
                          ) * c_vect_tmp(j_ci) * state_sign
                  endif

                end if
              enddo
            enddo
          enddo ! the first term in w_5
          ! the second one
          ! j in O_{\alpha}, k in O_{\beta}, b in V_{\beta})
          do j_state = i_start_ci, n_elec(1), 1 
            do k_state = i_start_ci, n_elec(2), 1
              do b_state = n_elec(2)+1, n_states, 1
                substitude(1,1) = j_state
                substitude(1,2) = k_state
                occupy(1,1)     = a_state
                occupy(1,2)     = b_state

                call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                if (j_task .eq. i_task) then
                    w_5 = w_5 - (&
                        ddot(n_basbas,&
                        ovlp_3ks(:,j_state,i_state,tmp_spin(1)),1,&
                        ovlp_3ks(:,k_state,b_state,tmp_spin(2)),1)&
                        ) * c_vect_tmp(j_ci)
                endif

              enddo
            enddo
          enddo ! the second term in w_5
          !write(use_unit,'(2X,I3,6(F16.8))') i_ci+1, V_s0, w_1, w_2, w_3, w_4, w_5
          !w_vect(i_ci) = V_s0 *c_0 + w_1 + w_2 + w_3 + w_4 + w_5
          w_vect(i_ci_local) = w_vect(i_ci_local) + w_1 + w_2 + w_3 + w_4 + w_5
          !w_vect(i_ci_local) = w_1 + w_2 + w_3 + w_4 + w_5
        enddo
        index_start(2) = n_elec(1)+1
      enddo  ! one-electron excitations in the alpha spin space

666   continue
    enddo

    if (allocated(occ_num_1)) then
        deallocate(occ_num_1)
    endif
    if (allocated(c_vect_tmp)) then
        deallocate(c_vect_tmp)
    endif

end subroutine evaluate_w_vector_1a

