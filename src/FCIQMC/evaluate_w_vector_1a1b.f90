!****s* FHI-aims/fciqmc/evaluate_w_vector_1a1b
!  NAME
!    evaluate_w_vector_1a1b
!  SYNOPSIS 

      subroutine evaluate_w_vector_1a1b()

!  PURPOSE
!
!  Calculate w_vector on-the-fly for a given c_vect
!  '''evaluate_w_vector_1a1b''' generate the w_vect of double excitations in both spin spaces
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
  use dimensions
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
    integer :: i_state, j_state, a_state, b_state
    integer :: c_state, d_state, k_state, l_state
    integer :: i_ci, i_ci_local, j_ci, k_ci, j_ci_start, j_ci_end
    integer :: errnum
    !character*128 :: info_str
    real*8  :: w_1, w_2, w_3, w_4, w_5
    real*8  :: V_s0,state_sign
    real*8  :: ddot
    integer :: tmp_spin(2)
    integer :: substitude(2,2), occupy(2,2)
    integer :: index_start(4)
    integer :: i_task, j_task, j_task_start, j_task_end
    integer :: mpierr, info,wnum

    real*8, allocatable, dimension(:) :: c_vect_tmp

    !allocate(c_vect_tmp(maxval(n_configuration_table,n_tasks)),stat=errnum)
    allocate(c_vect_tmp(n_configuration_table(1)),stat=errnum)
    call check_allocation(errnum, 'c_vect_tmp in CI')


    tmp_spin = 1
    if (n_spin .eq. 2) tmp_spin(2) = 2
    
    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 (1D) in CI')
    endif

    do i_task = 0, n_tasks-1, 1

      c_vect_tmp=0.0d0

      ! broadcast c_vect in i_task to all processors
      if (i_task .eq. myid) then
          !write(use_unit,*) "shen debug 1", i_task,n_configuration_table(i_task+1)
          c_vect_tmp(1:n_configuration_table(i_task+1)) = &
              c_vect(:)


      end if

      call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
      call mpi_bcast(c_vect_tmp,n_configuration_table(i_task+1), &
                     MPI_REAL8, i_task, mpi_comm_global, mpierr)

      !write(use_unit,*) "igor debug 1a1b_1",i_task, myid
      !i_ci = n_configuration_1a + n_configuration_1b
      !do i_state = i_start_ci, n_elec(1), 1
      !  do a_state = n_elec(1)+1, n_states, 1
      !    do j_state = i_start_ci, n_elec(2), 1
      !      do b_state = n_elec(2)+1, n_states, 1
      i_ci        = memory_table(1,1,3,myid+1)-1
      i_ci_local  = memory_table(1,2,3,myid+1)-1
      index_start = index_table(:,3,myid+1)
      do i_state = index_start(1), n_elec(1), 1
        do a_state = index_start(2), n_states, 1
          do j_state = index_start(3), n_elec(2), 1
            do b_state = index_start(4), n_states, 1


              i_ci_local = i_ci_local + 1
              if (i_ci_local .gt. memory_table(2,2,3,myid+1)) goto 666

              if (i_task .eq. 0) then
                  ! be careful. For seek_index_v01, should bring them out
                  occ_num_1 = occ_num_0
                  occ_num_1(i_state,1) = 0 ! occ_alpha
                  occ_num_1(a_state,1) = 1 ! vir_alpha
                  occ_num_1(j_state,2) = 0 ! occ_beta
                  occ_num_1(b_state,2) = 1 ! vir_beta
                  ! only calculate V_s0 one time
                  call evaluate_Hamiltonian_matrix_element(occ_num_1,occ_num_0,V_s0)
                  ! NOTICE :: initialize w_vect(i_ci) in the 0th processor
                  w_vect(i_ci_local) = V_s0 * c_0
                  !write(use_unit,*) "igor debug V_s0", V_s0, w_vect(i_ci), i_ci
              endif

              w_1 = 0.0d0

              do c_state = n_elec(1)+1, n_states, 1
                do d_state = n_elec(2)+1, n_states, 1
                  substitude(1,1) = i_state
                  substitude(1,2) = j_state
                  occupy(1,1)     = c_state
                  occupy(1,2)     = d_state

                  call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_1 = w_1 + &
                          ddot(n_basbas,&
                          ovlp_3ks(:,a_state,c_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,b_state,d_state,tmp_spin(2)),1)&
                          * c_vect_tmp(j_ci)
                  endif

                enddo
              enddo

              w_2 = 0.0d0
              do k_state = i_start_ci, n_elec(1), 1
                do l_state = i_start_ci, n_elec(2), 1
                  substitude(1,1) = k_state
                  substitude(1,2) = l_state
                  occupy(1,1)     = a_state
                  occupy(1,2)     = b_state

                  call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_2 = w_2 + &
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,i_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,l_state,j_state,tmp_spin(2)),1)&
                          * c_vect_tmp(j_ci)
                  endif

                enddo
              enddo

              w_3 = 0.0d0
              !--------------------------------------
              ! for the first term in w_3
              !--------------------------------------
              ! for the first term in the first term
              ! both exciate and occupy in the alpha spin space
              do k_state = i_start_ci, n_elec(1), 1
                do c_state = n_elec(1)+1, n_states, 1
                  if ((k_state.ne.i_state).and.(a_state.ne.c_state)) then
                    substitude(1,1) = i_state
                    substitude(2,1) = k_state
                    occupy(1,1)     = a_state
                    occupy(2,1)     = c_state

                    call ci_coefficient_sign(substitude,occupy,state_sign,1)
                    call seek_index_v03_4(substitude,occupy,j_ci,j_task)

                    if (j_task .eq. i_task) then
                        w_3 = w_3 - &
                            ddot(n_basbas,&
                            ovlp_3ks(:,k_state,c_state,tmp_spin(1)),1,&
                            ovlp_3ks(:,b_state,j_state,tmp_spin(2)),1)&
                            * c_vect_tmp(j_ci) * state_sign
                    endif
                  end if
                enddo
              enddo

              ! for the second term in the first term
              do k_state = i_start_ci, n_elec(2), 1
                do c_state = n_elec(2)+1, n_states, 1
                  substitude(1,1) = i_state
                  substitude(1,2) = k_state
                  occupy(1,1)     = a_state
                  occupy(1,2)     = c_state

                  call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_3 = w_3 + (&
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,j_state,tmp_spin(2)),1,&
                          ovlp_3ks(:,b_state,c_state,tmp_spin(2)),1)&
                          - &
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,c_state,tmp_spin(2)),1,&
                          ovlp_3ks(:,b_state,j_state,tmp_spin(2)),1)&
                          )* c_vect_tmp(j_ci)
                  endif

                enddo
              enddo
              !--------------------------------------
              ! for the second term in w_3
              !--------------------------------------
              do k_state = i_start_ci, n_elec(2), 1
                do c_state = n_elec(1)+1, n_states, 1
                  substitude(1,1) = i_state
                  substitude(1,2) = k_state
                  occupy(1,1)     = c_state
                  occupy(1,2)     = b_state

                  call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_3 = w_3 + &
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,j_state,tmp_spin(2)),1,&
                          ovlp_3ks(:,a_state,c_state,tmp_spin(1)),1)&
                          * c_vect_tmp(j_ci)
                  endif

                enddo
              enddo
              !--------------------------------------
              ! for the third term in w_3
              !--------------------------------------
              do k_state = i_start_ci, n_elec(1), 1
                do c_state = n_elec(2)+1, n_states, 1
                  substitude(1,1) = k_state
                  substitude(1,2) = j_state
                  occupy(1,1)     = a_state
                  occupy(1,2)     = c_state

                  call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_3 = w_3 + &
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,i_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,b_state,c_state,tmp_spin(2)),1)&
                          * c_vect_tmp(j_ci)
                  endif

                enddo
              enddo
              !--------------------------------------
              ! for the fourth term in w_3
              !--------------------------------------
              ! for the first term in the fourth term
              do k_state = i_start_ci, n_elec(1), 1
                do c_state = n_elec(1)+1, n_states, 1
                  substitude(1,1) = k_state
                  substitude(1,2) = j_state
                  occupy(1,1)     = c_state
                  occupy(1,2)     = b_state

                  call seek_index_v03_3(substitude,occupy,j_ci,j_task)
                  if (j_task .eq. i_task) then
                      w_3 = w_3 + (&
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,i_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,a_state,c_state,tmp_spin(1)),1)&
                          - &
                          ddot(n_basbas,&
                          ovlp_3ks(:,k_state,c_state,tmp_spin(1)),1,&
                          ovlp_3ks(:,a_state,i_state,tmp_spin(1)),1)&
                          )* c_vect_tmp(j_ci)
                  endif

                enddo
              enddo
              ! exciate and occupy both in the beta spin space
              do k_state = i_start_ci, n_elec(2), 1
                do c_state = n_elec(2)+1, n_states, 1
                  if ((k_state.ne.j_state).and.(c_state.ne.b_state)) then
                    substitude(1,2) = k_state
                    substitude(2,2) = j_state
                    occupy(1,2)     = c_state
                    occupy(2,2)     = b_state

                    call ci_coefficient_sign(substitude,occupy,state_sign,2)
                    call seek_index_v03_5(substitude,occupy,j_ci,j_task)
                    if (j_task .eq. i_task) then
                        w_3 = w_3 - &
                            ddot(n_basbas,&
                            ovlp_3ks(:,k_state,c_state,tmp_spin(2)),1,&
                            ovlp_3ks(:,a_state,i_state,tmp_spin(1)),1)&
                            * c_vect_tmp(j_ci) * state_sign
                    endif
                  end if
                enddo
              enddo

              w_3 = - w_3
              w_4 = 0.0d0
              w_5 = 0.0d0
              if (flag_single_excitation) then
                  ! There are two terms in w_4
                  ! the first one goes over V_{\alpha}
                  do c_state = n_elec(1)+1, n_states, 1
                    substitude(1,1) = i_state
                    occupy(1,1)     = c_state

                    call seek_index_v03_1(substitude,occupy,j_ci,j_task)
                    if (j_task .eq. i_task) then
                        w_4 = w_4 + &
                            ddot(n_basbas,&
                            ovlp_3ks(:,a_state,c_state,tmp_spin(1)),1,&
                            ovlp_3ks(:,b_state,j_state,tmp_spin(2)),1)&
                            * c_vect_tmp(j_ci)
                    endif

                  enddo

                  ! the second one goes over V_{\beta}
                  do c_state = n_elec(2)+1, n_states, 1
                    substitude(1,2) = j_state
                    occupy(1,2)     = c_state

                    call seek_index_v03_2(substitude,occupy,j_ci,j_task)
                    if (j_task .eq. i_task) then
                        w_4 = w_4 + &
                            ddot(n_basbas,&
                            ovlp_3ks(:,a_state,i_state,tmp_spin(1)),1,&
                            ovlp_3ks(:,b_state,c_state,tmp_spin(2)),1)&
                            * c_vect_tmp(j_ci)
                    endif
                  enddo

                  ! There are two terms in w_5
                  ! the first one goes over O_{\alpha}
                  do k_state = i_start_ci, n_elec(1), 1
                    substitude(1,1) = k_state
                    occupy(1,1)     = a_state

                    call seek_index_v03_1(substitude,occupy,j_ci,j_task)
                    if (j_task .eq. i_task) then
                        w_5 = w_5 - &
                            ddot(n_basbas,&
                            ovlp_3ks(:,k_state,i_state,tmp_spin(1)),1,&
                            ovlp_3ks(:,b_state,j_state,tmp_spin(2)),1)&
                            * c_vect_tmp(j_ci)
                    endif

                  enddo
                  ! the second one goes over O_{\beta}
                  do k_state = i_start_ci, n_elec(2), 1
                    substitude(1,2) = k_state
                    occupy(1,2)     = b_state

                    call seek_index_v03_2(substitude,occupy,j_ci,j_task)
                    if (j_task .eq. i_task) then
                        w_5 = w_5 - &
                            ddot(n_basbas,&
                            ovlp_3ks(:,k_state,j_state,tmp_spin(2)),1,&
                            ovlp_3ks(:,a_state,i_state,tmp_spin(1)),1)&
                            * c_vect_tmp(j_ci)
                    endif

                  enddo

              endif
              !write(use_unit,'(2X,I3,6(F16.8))') i_ci,V_s0, w_1, w_2, w_3, w_4, w_5
              !w_vect(i_ci) = V_s0 * c_0 + w_1 + &
              !    w_2 + w_3 + w_4 + w_5
              ! For parallel implementation

              w_vect(i_ci_local) = w_vect(i_ci_local) + w_1 + &
                  w_2 + w_3 + w_4 + w_5



            enddo
            index_start(4) = n_elec(2)+1
          enddo
          index_start(3) = i_start_ci 
        enddo
        index_start(2) = n_elec(1)+1 
      enddo ! two-electron excitations in both spin spaces

666   continue
    enddo

    if (allocated(occ_num_1)) then
        deallocate(occ_num_1)
    endif
    if (allocated(c_vect_tmp)) then
        deallocate(c_vect_tmp)
    endif

end subroutine evaluate_w_vector_1a1b
