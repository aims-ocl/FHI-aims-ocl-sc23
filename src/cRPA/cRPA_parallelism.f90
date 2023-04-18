!Functions to communicate
MODULE cRPA_parallelism
    USE mpi_tasks
    USE cRPA_storage
    USE cRPA_view
    USE cRPA_parallelism_storage
    USE cRPA_parallelism_planing
    IMPLICIT NONE

    ABSTRACT INTERFACE
        SUBROUTINE on_recv_template(data_recv)
            IMPORT message_matrix_recv
            type(message_matrix_recv), INTENT(IN) :: data_recv
        END SUBROUTINE on_recv_template

        SUBROUTINE on_can_send_template()

        END SUBROUTINE on_can_send_template
    END INTERFACE

    INTEGER :: COMM_MATRIX_DISTRIBUTION_CPU_NEW
    PARAMETER(COMM_MATRIX_DISTRIBUTION_CPU_NEW=10000)

    INTEGER :: COMM_MATRIX_DISTRIBUTION_DATA
    PARAMETER(COMM_MATRIX_DISTRIBUTION_DATA=10001)

    INTEGER :: COMM_MATRIX_DISTRIBUTION_CPU_FINISHED
    PARAMETER(COMM_MATRIX_DISTRIBUTION_CPU_FINISHED=10002)

    INTEGER :: COMM_MATRIX_DISTRIBUTION_SHUTDOWN
    PARAMETER(COMM_MATRIX_DISTRIBUTION_SHUTDOWN=10003)

     LOGICAL :: is_a_send_pending
    INTEGER :: mpi_request_last_send

    LOGICAL :: is_a_read_pending
    INTEGER :: mpi_request_last_recv

    logical :: active_communication
    logical, allocatable, DIMENSION(:) :: active_cpus

    REAL*8, ALLOCATABLE, DIMENSION(:) :: buffer_send, buffer_recv, data
    INTEGER :: buffer_send_len, buffer_recv_len

    type(custom_communicator) :: comm_basisfuncs, &
                                 comm_exp_coeffs, &
                                 comm_mo_coeffs

    CONTAINS

    SUBROUTINE communicate_resp_func_for_one_k_point(n_tasks,worklist)
        INTEGER,INTENT(IN) :: n_tasks
        type(basbas_distribution), DIMENSION(n_tasks), INTENT(IN), TARGET :: worklist

        type(message_matrix_send), ALLOCATABLE, DIMENSION(:) :: data_send
        INTEGER :: i_k, n_my_basbas, i_spin,i_row, i_packet
        INTEGER, DIMENSION(n_k_points) :: k_point_to_cpu_list
        INTEGER, ALLOCATABLE,DIMENSION(:) :: packet_to_cpu_list

        PROCEDURE(on_recv_template),POINTER :: do_on_recv
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: data_as_real_array

        n_my_basbas = worklist(myid+1)%n_columns
!TODO: check alloc
        ALLOCATE(data_send(n_k_points*worklist(myid+1)%n_rows))
        ALLOCATE(packet_to_cpu_list(n_k_points*worklist(myid+1)%n_rows))
        ALLOCATE(data_as_real_array(n_my_basbas,2,1))

        resp_func_k_point_worklist => worklist

        CALL get_k_point_to_cpu_list(n_k_points, k_point_to_cpu_list)

        do i_k = 1, n_k_points
            do i_row = 1, worklist(myid+1)%n_rows

                data_as_real_array(:,1,1) = &
                                  dble(response_func%one_omega_in_k(:,i_row,i_k))
                data_as_real_array(:,2,1) = &
                                  dimag(response_func%one_omega_in_k(:,i_row,i_k))

! if(ANY(isnan(data_as_real_array(:,:,1)))) write(use_unit,*) "on send nan",data_as_real_array(:,:,1)

                i_packet = (i_k-1)*worklist(myid+1)%n_rows+i_row
                CALL create_message_matrix_send_from_matrix(n_my_basbas, 2, 1,&
                                                            data_as_real_array(:,:,:) ,&
                                                            data_send(i_packet))
                data_send(i_packet)%tag = i_k
                data_send(i_packet)%offset = worklist(myid+1)%n_offset_rows + i_row

                packet_to_cpu_list(i_packet) = k_point_to_cpu_list(i_k)

            enddo
        enddo

        do_on_recv => handle_resp_func_one_k_gathering

        CALL communicate_three_dimensional_matrices(n_k_points*worklist(myid+1)%n_rows, &
                                                    data_send, &
                                                    packet_to_cpu_list, &
                                                    do_on_recv)

        DEALLOCATE(packet_to_cpu_list)
        DEALLOCATE(data_as_real_array)

        NULLIFY(resp_func_k_point_worklist)
        CALL print_chksum_resp_func_my_k_points()
    END SUBROUTINE communicate_resp_func_for_one_k_point

    SUBROUTINE map_k_point_to_my_k_points(n_k_points,i_k_point,i_my_k_point)
        INTEGER, INTENT(IN) :: n_k_points, i_k_point
        INTEGER, INTENT(OUT) :: i_my_k_point

        INTEGER :: i
        INTEGER, DIMENSION(n_k_points) :: k_point_to_cpu_list

        CALL get_k_point_to_cpu_list(n_k_points, k_point_to_cpu_list)

        i_my_k_point = 0
        do i= 1, n_k_points
            if(k_point_to_cpu_list(i) == myid) then! + 1) then
                i_my_k_point = i_my_k_point + 1

                if(i == i_k_point) then
                    return
                endif

            endif
        enddo

        i_my_k_point = -1

    END SUBROUTINE map_k_point_to_my_k_points

!TODO refactor
    SUBROUTINE map_my_k_point_to_k_point(n_k_points,i_my_k_point, i_k_point)
        INTEGER, INTENT(IN) :: n_k_points, i_my_k_point
        INTEGER, INTENT(OUT) :: i_k_point


        INTEGER :: i, i_res

        do i= 1, n_k_points
            CALL map_k_point_to_my_k_points(n_k_points, i, i_res)
            if(i_res == i_my_k_point) then
                i_k_point = 1
                return
            endif
        enddo

        i_k_point = -1


    END SUBROUTINE map_my_k_point_to_k_point

    SUBROUTINE get_k_point_to_cpu_list(n_k_points,k_point_to_cpu_list)
        INTEGER, INTENT(IN) :: n_k_points
        INTEGER,DIMENSION(n_k_points), INTENT(INOUT) :: k_point_to_cpu_list

        INTEGER :: i_k

        do i_k = 1, n_k_points
           k_point_to_cpu_list(i_k) = MOD(i_k, n_tasks) !+ 1
        enddo

    END SUBROUTINE get_k_point_to_cpu_list

    SUBROUTINE handle_resp_func_one_k_gathering(data_recv)
        type(message_matrix_recv), INTENT(IN) :: data_recv

        COMPLEX*16 :: basbas_pair
        INTEGER :: i,i_my_k_point, source_cpu
        INTEGER :: i_basbas1, i_basbas2

        CALL map_k_point_to_my_k_points(n_k_points, data_recv%tag, i_my_k_point)

        source_cpu = data_recv%source+1
!CALL write_debug("my_k_point"//num2str(i_my_k_point)//num2str(data_recv%tag))
!CALL write_debug("dims"//num2str(resp_func_k_point_worklist(source_cpu)%dim_row)// &
!                 num2str(data_recv%dim1) // &
!                 num2str(data_recv%dim2) // &
!                 num2str(data_recv%dim3))

        do i = 1, resp_func_k_point_worklist(source_cpu)%n_columns
            i_basbas1 = data_recv%offset
            i_basbas2 = i

! if(ANY(isnan(data_recv%data(i,:,1)))) write(use_unit,*) "on recv nan",data_recv%data(i,:,1)

            basbas_pair = DCMPLX(data_recv%data(i,1,1), 0.d0)
            basbas_pair = basbas_pair + DCMPLX(0.d0, &
                                              data_recv%data(i,2,1))

            response_func%one_omega_my_k_points(i_basbas1, &
                                                i_basbas2, &
                                                i_my_k_point) = basbas_pair

        enddo

    END SUBROUTINE handle_resp_func_one_k_gathering

    SUBROUTINE communicate_three_dimensional_matrices(n_matrices, mat, map_to_cpu_list, &
                                                      do_on_recv)

        INTEGER, INTENT(IN) :: n_matrices
        type(message_matrix_send), DIMENSION(n_matrices), INTENT(INOUT) :: mat
        INTEGER, DIMENSION(n_matrices), INTENT(IN) :: map_to_cpu_list
        PROCEDURE(on_recv_template),POINTER :: do_on_recv

        INTEGER :: res_dim(3), mpi_result, dest_cpu,send_index, mpi_request

        INTEGER :: status(MPI_STATUS_SIZE), message_type_data, my_data_len, data_len
        LOGICAL :: can_send


        CALL set_res_dim(n_matrices,mat)

        my_data_len = product(res_dim)

        CALL MPI_ALLREDUCE(my_data_len, data_len, 1, &
                           MPI_INTEGER, MPI_MAX, mpi_comm_global , mpi_result)


        CALL init_non_blocking_communication(data_len+30,data_len+30)

!        CALL write_debug(" start")

            !            CALL MPI_BCAST (i,1,MPI_INTEGER, &
            !                            run_cpu-1,mpi_comm_global,mpi_result)
            !                    write(use_unit,*) myid,"bdown",res_dim, buffer_send_len, buffer_recv_len

        send_index=1

        do while (active_communication)

           CALL is_send_ready(can_send)
           if(can_send .AND. send_index<=n_matrices) then
               dest_cpu = map_to_cpu_list(send_index)

!TODO debug setting!
!.OR. (dest_cpu==myid .AND. n_tasks ==1)
               !if(dest_cpu/=myid) then
!                   write(use_unit,*) "Send index", send_index
                   CALL set_send_message_3d_mat(mat(send_index))
                   CALL send_message_3d_matrix(dest_cpu,mat(send_index),mpi_result)
               !endif

               if(send_index>1) then
                   CALL free_message_matrix_send(mat(send_index-1))
               endif

               send_index = send_index+1
           endif

           CALL is_send_ready(can_send)
           if(can_send .AND. send_index==n_matrices+1) then
!                write(use_unit,*) myid, "SENDS SHUTDOWN"
                CALL send_message_cpu_finished_to_root(mpi_result)

                CALL free_message_matrix_send(mat(send_index-1))

                send_index = send_index+1
           endif

           if(send_index>n_matrices) then
           !                  write(use_unit,*) ,myid, "sleeps"
            !                  CALL sleep(1)
           endif

           CALL read_message_if_any(do_on_recv,mpi_result)

      end do

!    write(use_unit,*) myid,"SHUTTED DOWN"

    CALL finalize_non_blocking_communication(mpi_result)

!    write(use_unit,*) myid,"FINISHED"

        CONTAINS

            SUBROUTINE set_res_dim(n_matrices, mat)
                INTEGER, INTENT(IN) :: n_matrices
                type(message_matrix_send), DIMENSION(n_matrices), INTENT(IN) :: mat

                res_dim(1) = maxval(mat(1:n_matrices)%dim1)
                res_dim(2) = maxval(mat(1:n_matrices)%dim2)
                res_dim(3) = maxval(mat(1:n_matrices)%dim3)
            END SUBROUTINE

            SUBROUTINE set_send_message_3d_mat(data_type)
                type(message_matrix_send), INTENT(INOUT) :: data_type

!                data_type%source = i
                data_type%message_id=COMM_MATRIX_DISTRIBUTION_DATA
                data_type%sources_index = send_index
!                data_type%dim1 = res_dim(1)
!                data_type%dim2 = res_dim(2)
!                data_type%dim3 = res_dim(3)
!                data_type%data(1:data_len) = mat(i)%data(1:data_len)
            END SUBROUTINE set_send_message_3d_mat

    END SUBROUTINE communicate_three_dimensional_matrices

    SUBROUTINE init_non_blocking_communication(max_buffer_send_len,max_buffer_recv_len)
        INTEGER, INTENT(IN) :: max_buffer_send_len,max_buffer_recv_len

        INTEGER :: mpi_result

        buffer_send_len = max_buffer_send_len + COMM_MATRIX_HEADER
        buffer_recv_len = max_buffer_recv_len + COMM_MATRIX_HEADER

        is_a_send_pending = .FALSE.
        is_a_read_pending = .FALSE.

        ALLOCATE(buffer_send(buffer_send_len))
        ALLOCATE(buffer_recv(buffer_recv_len))

        buffer_send=0
        buffer_recv=0

!        CALL write_stdout("Starting nonblocking communication")
!        CALL write_stdout("buffersizes: send:"//num2str(buffer_send_len) &
!                          // " - read:" // num2str(buffer_recv_len))

        CALL MPI_BARRIER(mpi_comm_global,mpi_result)

        ALLOCATE(active_cpus(n_tasks))
        active_cpus = .TRUE.
        active_communication = .TRUE.


    END SUBROUTINE init_non_blocking_communication

    SUBROUTINE finalize_non_blocking_communication(mpi_result)
        INTEGER,INTENT(OUT) :: mpi_result

!        CALL write_stdout("Finalizing non blocking communication")

        if(active_communication) then
            stop "cannot finalize communication - it is still activ"
        endif

        DEALLOCATE(buffer_send)
        DEALLOCATE(buffer_recv)

        DEALLOCATE(active_cpus)

        CALL MPI_BARRIER(mpi_comm_global,mpi_result)
    END SUBROUTINE finalize_non_blocking_communication

    SUBROUTINE is_send_ready(flag)

        LOGICAL, INTENT(INOUT) :: flag
        INTEGER :: mpi_status(MPI_STATUS_SIZE), mpi_result
        LOGICAL :: test_state

        if(is_a_send_pending) then
            CALL MPI_TEST(mpi_request_last_send,test_state, mpi_status,mpi_result)

            is_a_send_pending = .NOT. test_state
        endif


        if(is_a_send_pending) then
            flag = .FALSE.
        else
            flag = .TRUE.
        endif

    END SUBROUTINE

    SUBROUTINE send_message_3d_matrix(dest_cpu,matrix_3d,mpi_result)

       INTEGER, INTENT(IN) :: dest_cpu
       type(message_matrix_send),INTENT(IN) :: matrix_3d
       INTEGER, INTENT(OUT) :: mpi_result

       INTEGER :: message_type, data_len

       if(is_a_send_pending) then
          write(use_unit,*) myid, " sent not ready"
          stop
          return
       endif

       buffer_send(1)=real(myid)
       buffer_send(2)=real(matrix_3d%message_id)
       buffer_send(3)=real(matrix_3d%sources_index)
       buffer_send(4)=real(matrix_3d%offset)
       buffer_send(5)=real(matrix_3d%tag)
       buffer_send(6)=real(matrix_3d%dim1)
       buffer_send(7)=real(matrix_3d%dim2)
       buffer_send(8)=real(matrix_3d%dim3)
       buffer_send(9)=real(matrix_3d%i_cell1)
       buffer_send(10)=real(matrix_3d%i_cell2)
       buffer_send(11)=real(matrix_3d%i_cell3)


       data_len = matrix_3d%dim1 * matrix_3d%dim2 * matrix_3d%dim3

       buffer_send(12:data_len+11) = matrix_3d%data(1:data_len)

!       CALL write_debug("SENDS TO " // int2str(dest_cpu) &
!                        // "MSGID:" // int2str(matrix_3d%message_id) &
!                        // "TAG:" // int2str(matrix_3d%tag))


       call MPI_ISEND(buffer_send, data_len+11, MPI_REAL8, dest_cpu, &
                     COMM_MATRIX_DISTRIBUTION_DATA, mpi_comm_global, &
                     mpi_request_last_send, mpi_result)

!                    matrix_3d%i_cell1,matrix_3d%i_cell2,matrix_3d%i_cell3

       is_a_send_pending = .TRUE.

    END SUBROUTINE send_message_3d_matrix

    SUBROUTINE send_message(dest_cpu,message_id,mpi_result)
       INTEGER, INTENT(IN) :: dest_cpu, message_id
       INTEGER, INTENT(OUT) :: mpi_result

        CALL send_message_with_tag(dest_cpu,message_id,0, mpi_result)

    END SUBROUTINE send_message

    SUBROUTINE send_message_with_tag(dest_cpu, message_id, &
                                     tag,mpi_result)
       INTEGER, INTENT(IN) :: dest_cpu, message_id,tag
       INTEGER, INTENT(OUT) :: mpi_result

       type(message_matrix_send) :: message

       CALL create_message_matrix_send(1,1,1,message)

       message%source=myid
       message%message_id = message_id
       message%tag = tag

       CALL send_message_3d_matrix(dest_cpu,message,mpi_result)

       CALL free_message_matrix_send(message)

    END SUBROUTINE send_message_with_tag

    SUBROUTINE send_message_cpu_finished_to_root(mpi_result)
        INTEGER, INTENT(OUT) :: mpi_result
        LOGICAL ::flag

        flag = .TRUE.

  !      if(myid /= 0) then
            if(is_a_send_pending) then
                do while (flag)
                    CALL is_send_ready(flag)
                end do
            endif
            CALL send_message(0, COMM_MATRIX_DISTRIBUTION_CPU_FINISHED,mpi_result)
  !      else
  !          CALL handle_cpu_finished(0,mpi_result)
  !      endif
    END SUBROUTINE send_message_cpu_finished_to_root

    SUBROUTINE has_message(flag,mpi_result)
        LOGICAL,INTENT(INOUT) :: flag
        INTEGER,INTENT(INOUT) :: mpi_result

        INTEGER :: message_type,status(MPI_STATUS_SIZE)

        message_type = COMM_MATRIX_DISTRIBUTION_DATA

        !CALL mpi_iprobe(MPI_ANY_SOURCE, message_type, &
        !                 mpi_comm_global, flag, status, mpi_result)

    END SUBROUTINE

    SUBROUTINE is_read_ready(flag)

        LOGICAL, INTENT(INOUT) :: flag
        INTEGER :: mpi_status(MPI_STATUS_SIZE), mpi_result
        LOGICAL :: test_state

        if(is_a_read_pending) then
            CALL MPI_TEST(mpi_request_last_send,test_state, mpi_status,mpi_result)

            is_a_read_pending = .NOT. test_state
        endif


        if(is_a_read_pending) then
            flag = .FALSE.
        else
            flag = .TRUE.
        endif

    END SUBROUTINE

    SUBROUTINE read_message_if_any(do_on_recv,mpi_result)
       PROCEDURE(on_recv_template),POINTER :: do_on_recv
       INTEGER, INTENT(OUT) :: mpi_result

       INTEGER :: mpi_status(MPI_STATUS_SIZE)

       LOGICAL :: is_read_complete


       LOGICAL :: flag

        if(.NOT. is_a_read_pending) then
            CALL has_message(flag,mpi_result)

            if(flag) then
                CALL read_message(do_on_recv,mpi_result)
            endif
            return
        endif

        if(is_a_read_pending) then
           CALL MPI_TEST(mpi_request_last_recv,is_read_complete,mpi_status,mpi_result)
           if(is_read_complete) then
              CALL read_message(do_on_recv,mpi_result)
           endif
        endif

    END SUBROUTINE read_message_if_any

    SUBROUTINE read_message(do_on_recv,mpi_result)
       PROCEDURE(on_recv_template),POINTER :: do_on_recv
       INTEGER, INTENT(OUT) :: mpi_result

       type(message_matrix_recv) :: data_recv
       INTEGER :: data_len, message_type, mpi_request, mpi_status(MPI_STATUS_SIZE), &
                  dim1,dim2,dim3

       LOGICAL :: is_read_complete

       message_type = COMM_MATRIX_DISTRIBUTION_DATA

!       write(use_unit,*) myid, "start irecv"

        if(.NOT. is_a_read_pending) then

            call MPI_IRECV(buffer_recv, buffer_recv_len, MPI_REAL8, &
                           MPI_ANY_SOURCE, message_type, mpi_comm_global, &
                           mpi_request_last_recv,mpi_result)

            is_a_read_pending = .TRUE.
            return
        end if

!        CALL MPI_TEST(mpi_request_last_recv,is_read_complete,mpi_status,mpi_result)
!        write(use_unit,*) myid, "stop irecv"

        is_read_complete = .TRUE.

        if(is_read_complete) then
           CALL MPI_WAIT(mpi_request_last_recv,mpi_status,mpi_result)

           dim1=int(buffer_recv(6))
           dim2=int(buffer_recv(7))
           dim3=int(buffer_recv(8))


           CALL create_message_matrix_recv(dim1,dim2,dim3,data_recv)

           data_recv%source=int(buffer_recv(1))
           data_recv%message_id=int(buffer_recv(2))
           data_recv%sources_index=int(buffer_recv(3))
           data_recv%offset=int(buffer_recv(4))
           data_recv%tag=int(buffer_recv(5))
           data_recv%dim1=dim1
           data_recv%dim2=dim2
           data_recv%dim3=dim3
           data_recv%i_cell1=int(buffer_recv(9))
           data_recv%i_cell2=int(buffer_recv(10))
           data_recv%i_cell3=int(buffer_recv(11))

!        write(use_unit,*) myid," RECVED FROM ", data_recv%source, "MSGID:", data_recv%message_id , "TAG:", data_recv%tag, &
!                                                                     data_recv%i_cell1,data_recv%i_cell2,data_recv%i_cell3


           data_len = dim1 * dim2 * dim3
!write(use_unit,*) myid,"reshaping"
           data_recv%data(1:dim1,1:dim2,1:dim3)=reshape(buffer_recv(12:data_len+11),(/dim1, dim2, dim3/))
!write(use_unit,*) myid,"reshape done"

           if(data_recv%message_id == COMM_MATRIX_DISTRIBUTION_DATA) then
                CALL do_on_recv(data_recv)
           endif

           if(data_recv%message_id == COMM_MATRIX_DISTRIBUTION_CPU_FINISHED) then
              CALL handle_cpu_finished(data_recv%source,mpi_result)
           endif

           if(data_recv%message_id == COMM_MATRIX_DISTRIBUTION_SHUTDOWN) then
               active_communication = .FALSE.
           endif

           is_a_read_pending = .FALSE.
           CALL free_message_matrix_recv(data_recv)
       end if
    END SUBROUTINE

    SUBROUTINE handle_cpu_finished(source,mpi_result)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(OUT) :: mpi_result


        active_cpus(source+1) = .FALSE.

!        write(use_unit,*) "GOT SHUTDOWN from ",source ," - ", active_cpus


        if(ALL(active_cpus.EQV. .FALSE.)) then
!            CALL write_stdout("-- SEND ALL SHUTDOWN --")
            CALL send_shutdown(mpi_result)
            active_communication = .FALSE.
        endif
     END SUBROUTINE handle_cpu_finished


    INTEGER FUNCTION get_bvk_idx(i_cell1,i_cell2,i_cell3)
        INTEGER, INTENT(IN) :: i_cell1,i_cell2,i_cell3
        INTEGER :: i

        do i=1, n_cells_bvk
            if(cell_index_bvk(i,1) == i_cell1 &
               .AND. &
               cell_index_bvk(i,2) == i_cell2 &
               .AND. &
               cell_index_bvk(i,3) == i_cell3) then

               get_bvk_idx = i
               return
            endif
        enddo

        stop "could not find bvk idx"
    END FUNCTION get_bvk_idx

    SUBROUTINE send_shutdown(mpi_result)
       INTEGER, INTENT(OUT) :: mpi_result
       INTEGER :: dest_cpu
       LOGICAL :: can_send

       do dest_cpu=1,n_tasks-1
          if(myid/=dest_cpu) then
             do
                CALL is_send_ready(can_send)
                if(can_send) then
                  exit
                endif
             enddo

 !            write(use_unit,*) myid,"send shutdown to ",dest_cpu
             CALL send_message(dest_cpu,COMM_MATRIX_DISTRIBUTION_SHUTDOWN,mpi_result)
          endif
       end do

    END SUBROUTINE

    SUBROUTINE setup_communication_mo_coeffs(this, n_k_points_task, my_mo_coeffs)
        type(mo_coefficients), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: n_k_points_task
        DOUBLE COMPlEX, DIMENSION(this%n_basis, this%n_states, &
                                  this%n_spin, n_k_points_task), INTENT(IN) :: my_mo_coeffs


        INTEGER :: i_k_point, n_my_k_points
        INTEGER :: i_k_point_my, id_of_k_point, datalen, mpi_result
        DOUBLE COMPlEX, DIMENSION(this%n_basis, this%n_states, &
                                  this%n_spin) :: tmp

        do i_k_point = 1,this%n_k_points
            this%k_point_to_id(i_k_point) = mod(i_k_point,comm_mo_coeffs%n_tasks)
            this%k_point_to_k_point_local(i_k_point) = (i_k_point-1)/comm_mo_coeffs%n_tasks + 1
        enddo

        n_my_k_points = COUNT(this%k_point_to_id==comm_mo_coeffs%myid)

        CALL allocate_for_my_k_points(this,n_my_k_points)

        CALL MPI_BARRIER(mpi_comm_global, mpi_result)

        datalen = this%n_states*this%n_basis*this%n_spin
        do i_k_point = 1,this%n_k_points
!        CALL write_debug("mo comm"//num2str(i_k_point)//&
!                         num2str(size(mo_coeffs%c(:,:,:,i_k_point))) &
!                         //num2str(datalen))
            id_of_k_point = MOD(i_k_point, n_tasks)
            i_k_point_my = (i_k_point-1)/n_tasks + 1

            if(myid == id_of_k_point) then
                tmp(:,:,:) = my_mo_coeffs(:,:,:,i_k_point_my)
            endif

            CALL MPI_BCAST(tmp, &
                           datalen, &
                           MPI_COMPLEX16, &
                           id_of_k_point, &
                           mpi_comm_global, &
                           mpi_result)

            if(this%k_point_to_id(i_k_point)==comm_mo_coeffs%myid) then
                this%c(:,:,:,this%k_point_to_k_point_local(i_k_point)) = tmp(:,:,:)
            endif
        enddo

    END SUBROUTINE setup_communication_mo_coeffs

    SUBROUTINE communicate_mo_coeffs(this, i_k_point, mo_coeff)
        type(mo_coefficients), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: i_k_point
        DOUBLE COMPlEX, DIMENSION(this%n_basis, this%n_states, &
                                  this%n_spin), INTENT(INOUT) :: mo_coeff
        INTEGER :: mpi_result

        if(this%k_point_to_id(i_k_point) == comm_mo_coeffs%myid) then
            mo_coeff(:,:,:) = this%c(:,:,:,this%k_point_to_k_point_local(i_k_point))
        endif

        CALL MPI_BCAST(mo_coeff, size(mo_coeff), &
                       MPI_COMPLEX16, &
                       this%k_point_to_id(i_k_point), &
                       comm_mo_coeffs%mpi_comm, &
                       mpi_result)

    END SUBROUTINE communicate_mo_coeffs

    SUBROUTINE communicate_exp_coeffs_integration_cell(max_n_basis_sp,max_n_basbas_sp, &
                                                       i_id_basbas_owner,i_cell, &
                                                       int_dist, &
                                                       data_send, &
                                                       data_recv)
        INTEGER,INTENT(IN) :: max_n_basis_sp,max_n_basbas_sp, i_id_basbas_owner,i_cell
        type(exp_coeff_calculation_distribution) :: int_dist
        REAL*8, DIMENSION(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2), INTENT(IN) :: data_send
        REAL*8, DIMENSION(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2), INTENT(OUT) :: data_recv


        INTEGER :: i_size, mpi_result

        if(int_dist%i_cell_to_id(i_cell) == i_id_basbas_owner) then
            data_recv(:,:,:,:) = data_send(:,:,:,:)
            return
        endif

        i_size = (max_n_basis_sp**2)*max_n_basbas_sp*2

        if(int_dist%i_cell_to_id(i_cell) == comm_exp_coeffs%myid) then
            CALL MPI_SEND(data_send,i_size,MPI_REAL8, &
                          i_id_basbas_owner, 0, mpi_result)
        endif

        if(comm_exp_coeffs%myid == i_id_basbas_owner) then
            CALL MPI_RECV(data_recv,i_size,MPI_REAL8, &
                          int_dist%i_cell_to_id(i_cell), 0, mpi_result)
        endif
    END SUBROUTINE communicate_exp_coeffs_integration_cell

!This function is very slow! It could be significantly improved by reducing the amount
!of communication calls. The amount of bytes sent is small.
    SUBROUTINE communicate_exp_coeffs_headers(this)
        type(expansion_coefficients) :: this
        INTEGER :: i_basbas, i_id,i_pos, &
                   max_compressed_size, &
                   max_alloc_size, &
                   max_row_extend, &
                   max_column_extend, &
                   n_basbas_my, &
                   i_cur_size, &
                   total_compressed_size, &
                   mpi_result

        REAL*8 :: d_cur_size, d_total_size, d_ratio

        this%basbas_to_id(:)=-1
        do i_basbas = 1,this%n_basbas
            if(this%coeff(i_basbas)%is_stored_locally) then
                this%basbas_to_id(i_basbas) = comm_exp_coeffs%myid
            endif
        enddo

        CALL communicate_array(this%n_basbas, this%basbas_to_id)
!if(comm_exp_coeffs%myid==0) write(use_unit,*) "basbas_to_id", this%basbas_to_id
        do i_basbas = 1,this%n_basbas
            CALL communicate_header_data(this%n_cells+1, &
                                         this%coeff(i_basbas)%cell_start)
            CALL communicate_header_data(this%n_cells+1, &
                                         this%coeff(i_basbas)%cell_end)
            CALL communicate_header_data(this%n_cells+1, &
                                         this%coeff(i_basbas)%dictionary_start)
            CALL communicate_header_data(this%n_cells+1, &
                                         this%coeff(i_basbas)%dictionary_end)

            CALL communicate_header_value(this%coeff(i_basbas)%last_stored_cell)
            CALL communicate_header_value(this%coeff(i_basbas)%dictionary_length)
            CALL communicate_header_value(this%coeff(i_basbas)%dictionary_pos)
            CALL communicate_header_value(this%coeff(i_basbas)%c_compressed_size)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_first_filled_row_onsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_last_filled_row_onsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_extend_filled_row_onsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_first_filled_column_onsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_last_filled_column_onsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_extend_filled_column_onsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_first_filled_row_offsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_last_filled_row_offsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_extend_filled_row_offsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_first_filled_column_offsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_last_filled_column_offsite)
            CALL communicate_header_value(this%coeff(i_basbas)%all_cells_extend_filled_column_offsite)


            if(this%basbas_to_id(i_basbas) /= comm_exp_coeffs%myid) then
                DEALLOCATE(this%coeff(i_basbas)%dictionary)
                ALLOCATE(this%coeff(i_basbas)%dictionary(this%coeff(i_basbas)%dictionary_length))
            endif

            CALL communicate_header_data(this%coeff(i_basbas)%dictionary_length, &
                                         this%coeff(i_basbas)%dictionary)

        enddo

        max_alloc_size = maxval(this%coeff(:)%c_compressed_size)
        CALL write_stdout("Expansion coeffs max alloc size:"//num2str(max_alloc_size))

        total_compressed_size = 0
        max_compressed_size = 0
        do i_basbas = 1, this%n_basbas
            i_cur_size = this%coeff(i_basbas)%cell_end(this%n_cells+1)
            total_compressed_size = total_compressed_size + i_cur_size
            max_compressed_size = max(max_compressed_size, &
                                      i_cur_size)
        enddo

        d_cur_size = to_megabyte(dble(max_compressed_size))
        d_total_size = to_megabyte(dble((this%n_cells+1)*n_basis**2))
        CALL write_stdout("Expansion coeffs max compressed size in Mb:"//num2str(d_cur_size)// &
                          "of"//num2str(d_total_size))

        d_cur_size = to_megabyte(dble(total_compressed_size))
        d_total_size = to_megabyte(dble(this%n_basbas*(this%n_cells+1)*n_basis**2))
        CALL write_stdout("Expansion coeffs total compressed size in Mb:"//num2str(d_cur_size)// &
                          "of"//num2str(d_total_size))


        max_row_extend = maxval(this%coeff(:)%all_cells_extend_filled_row_onsite)
        CALL write_stdout("Expansion coeffs max row onsite extend:"//num2str(max_row_extend)//"of"//num2str(n_basis))
        this%n_max_all_cells_filled_rows=max_row_extend

        max_column_extend = maxval(this%coeff(:)%all_cells_extend_filled_column_onsite)
        CALL write_stdout("Expansion coeffs max column onsite extend:"//num2str(max_column_extend)//"of"//num2str(n_basis))


        max_row_extend = maxval(this%coeff(:)%all_cells_extend_filled_row_offsite)
        CALL write_stdout("Expansion coeffs max row offsite extend:"//num2str(max_row_extend)//"of"//num2str(n_basis))

        max_column_extend = maxval(this%coeff(:)%all_cells_extend_filled_column_offsite)
        CALL write_stdout("Expansion coeffs max column offsite extend:"//num2str(max_column_extend)//"of"//num2str(n_basis))


        this%n_max_all_cells_filled_rows=max(this%n_max_all_cells_filled_rows,&
                                             max_column_extend)


        this%n_basbas_local_max = 0
        do i_id = 0, comm_exp_coeffs%n_tasks-1
            this%n_basbas_local_max =max(this%n_basbas_local_max, &
                                         COUNT(this%basbas_to_id==i_id))
        enddo
        this%n_basbas_local = COUNT(this%basbas_to_id==comm_exp_coeffs%myid)

        CALL write_stdout("Expansion coeffs max per task:"//num2str(this%n_basbas_local_max))


        ALLOCATE(this%basbas_local_by_id(0:comm_exp_coeffs%n_tasks-1, &
                                         this%n_basbas_local_max))

        this%basbas_local_by_id(:,:) = EXP_COEFFICIENT_UNASSIGNED

        do i_id = 0,comm_exp_coeffs%n_tasks-1
            i_pos = 1
            do i_basbas = 1, this%n_basbas
                if(this%basbas_to_id(i_basbas) == i_id) then
                    this%basbas_local_by_id(i_id,i_pos) = i_basbas
                    i_pos = i_pos + 1
                endif
            enddo
        enddo

!TODO: make more strict
        ALLOCATE(this%requested_coeff_compressed(max_compressed_size * this%n_basbas_local_max))


!write(use_unit,*) "ec myid", comm_exp_coeffs%myid, &
!           this%basbas_local_by_id(comm_exp_coeffs%myid,:)
return
       do i_basbas = 1,this%n_basbas
       if(myid==0) write(use_unit,*) "row extends", i_basbas, &
            this%coeff(i_basbas)%all_cells_first_filled_row_onsite, &
            this%coeff(i_basbas)%all_cells_last_filled_row_onsite, &
            this%coeff(i_basbas)%all_cells_extend_filled_row_onsite, &
            this%coeff(i_basbas)%all_cells_first_filled_row_offsite,&
            this%coeff(i_basbas)%all_cells_last_filled_row_offsite,&
            this%coeff(i_basbas)%all_cells_extend_filled_row_offsite, &
            this%coeff(i_basbas)%all_cells_first_filled_column_offsite,&
            this%coeff(i_basbas)%all_cells_last_filled_column_offsite,&
            this%coeff(i_basbas)%all_cells_extend_filled_column_offsite
       enddo
       CALL mpi_barrier(mpi_comm_global,mpi_result)
    CONTAINS

        SUBROUTINE communicate_array(n_array, array)
            INTEGER, INTENT(IN) :: n_array
            INTEGER, DIMENSION(n_array), INTENT(INOUT) :: array

            INTEGER, DIMENSION(n_array) :: array_buffer
            INTEGER :: mpi_result

            array_buffer(:) = array(:)
            CALL MPI_ALLREDUCE(array_buffer, &
                               array, &
                               n_array,&
                               MPI_INTEGER, &
                               MPI_MAX, &
                               comm_exp_coeffs%mpi_comm, &
                               mpi_result)


        END SUBROUTINE communicate_array

        SUBROUTINE communicate_header_data(n_array, array)
            INTEGER, INTENT(IN) :: n_array
            INTEGER, DIMENSION(n_array), INTENT(INOUT) :: array
            INTEGER :: mpi_result

            CALL MPI_BCAST (array,&
                            n_array, &
                            MPI_INTEGER, &
                            this%basbas_to_id(i_basbas), &
                            comm_exp_coeffs%mpi_comm, &
                            mpi_result)


        END SUBROUTINE communicate_header_data

       SUBROUTINE communicate_header_value(value)
            INTEGER, INTENT(INOUT) :: value
            INTEGER :: mpi_result

            CALL MPI_BCAST (value,&
                            1, &
                            MPI_INTEGER, &
                            this%basbas_to_id(i_basbas), &
                            comm_exp_coeffs%mpi_comm, &
                            mpi_result)


        END SUBROUTINE communicate_header_value

        REAL*8 FUNCTION to_megabyte(d_doubles)
            REAL*8, INTENT(IN) :: d_doubles

            to_megabyte = (d_doubles* 8.d0) / (1024.d0**2 )
        END FUNCTION

    END SUBROUTINE communicate_exp_coeffs_headers

    SUBROUTINE communicate_exp_coeffs_compressed(this, i_basbas)
        type(expansion_coefficients) :: this
        INTEGER,INTENT(IN) :: i_basbas
        INTEGER :: mpi_result,i_size

        i_size = this%coeff(i_basbas)%cell_end(this%n_cells+1)

        if(this%basbas_to_id(i_basbas) == comm_exp_coeffs%myid) then
            this%requested_coeff_compressed(1:i_size) = &
                                      this%coeff(i_basbas)%c_compressed(1:i_size)
        endif

        CALL MPI_BCAST (this%requested_coeff_compressed,&
                        i_size, &
                        MPI_REAL8, &
                        this%basbas_to_id(i_basbas), &
                        comm_exp_coeffs%mpi_comm, &
                        mpi_result)


    END SUBROUTINE communicate_exp_coeffs_compressed

    SUBROUTINE communicate_exp_coeffs_compressed_all(this, i_id_master)
        type(expansion_coefficients) :: this
        INTEGER,INTENT(IN) :: i_id_master
        INTEGER :: mpi_result,i_size, i_pos, i_basbas


        i_pos = 1
        do i_basbas = 1,this%n_basbas

              if(this%basbas_to_id(i_basbas)/=i_id_master) CYCLE

              i_size = this%coeff(i_basbas)%cell_end(this%n_cells+1)

              if(comm_exp_coeffs%myid==i_id_master) then
                  this%requested_coeff_compressed(i_pos:i_size+i_pos-1) = &
                                                this%coeff(i_basbas)%c_compressed(1:i_size)
              endif
              i_pos = i_pos + i_size
        enddo

        if(comm_exp_coeffs%myid==i_id_master) then
           this%requested_id = EXP_COEFFICIENT_NO_REQUESTED
        endif

        CALL MPI_BCAST (this%requested_coeff_compressed,&
                        i_pos-1, &
                        MPI_REAL8, &
                        i_id_master, &
                        comm_exp_coeffs%mpi_comm, &
                        mpi_result)

        if(comm_exp_coeffs%myid/=i_id_master) then
           i_pos = 1

           do i_basbas = 1,this%n_basbas
              if(this%basbas_to_id(i_basbas)/=i_id_master) CYCLE
              if(i_basbas==EXP_COEFFICIENT_UNASSIGNED) CYCLE

              i_size = this%coeff(i_basbas)%cell_end(this%n_cells+1)
              ALLOCATE(this%coeff(i_basbas)%c_compressed(i_size))

              this%coeff(i_basbas)%c_compressed(1:i_size) = &
                           this%requested_coeff_compressed(i_pos:i_size+i_pos-1)

              i_pos = i_pos + i_size
           enddo

           this%requested_id = i_id_master
        endif


    END SUBROUTINE communicate_exp_coeffs_compressed_all


    SUBROUTINE communicate_M_basbas(this, i_basbas,M)
        type(expansion_coefficients) :: this
        INTEGER,INTENT(IN) :: i_basbas
        DOUBLE COMPLEX, DIMENSION(n_states, n_states,n_spin),INTENT(INOUT) :: M
        INTEGER :: mpi_result,i_size

        i_size = (n_states**2) * n_spin

        CALL MPI_BCAST (M,&
                        i_size, &
                        MPI_COMPLEX16, &
                        this%basbas_to_id(i_basbas), &
                        comm_exp_coeffs%mpi_comm, &
                        mpi_result)

    END SUBROUTINE communicate_M_basbas

!    SUBROUTINE communicate_M_basbas_nonblocking_start(this, i_basbas,M,mpi_request)
!        type(expansion_coefficients) :: this
!        INTEGER,INTENT(IN) :: i_basbas
!        DOUBLE COMPLEX, DIMENSION(n_states, n_states,n_spin),INTENT(INOUT) :: M
!        INTEGER, INTENT(OUT) :: mpi_request
!        INTEGER :: mpi_result,i_size
!
!        i_size = (n_states**2) * n_spin
!
!        CALL MPI_IBCAST (M,&
!                        i_size, &
!                        MPI_COMPLEX16, &
!                        this%basbas_to_id(i_basbas), &
!                        comm_exp_coeffs%mpi_comm, &
!                        mpi_request, &
!                        mpi_result)
!
!    END SUBROUTINE communicate_M_basbas_nonblocking_start
!
!    SUBROUTINE communicate_M_basbas_nonblocking_finish(mpi_request)
!        INTEGER,INTENT(IN) :: mpi_request
!        INTEGER :: mpi_status,mpi_result
!
!        CALL MPI_WAIT(mpi_request,mpi_status,mpi_result)
!    END SUBROUTINE communicate_M_basbas_nonblocking_finish


    SUBROUTINE communicate_M_basbas_multiple(this, i_id_owner,n_count,M)
        type(expansion_coefficients) :: this
        INTEGER,INTENT(IN) :: i_id_owner,n_count
        DOUBLE COMPLEX, DIMENSION(n_states, n_states,n_spin,n_count),INTENT(INOUT) :: M
        INTEGER :: mpi_result,i_size

        i_size = (n_states**2) * n_spin * n_count

        CALL MPI_BCAST (M,&
                        i_size, &
                        MPI_COMPLEX16, &
                        i_id_owner, &
                        comm_exp_coeffs%mpi_comm, &
                        mpi_result)

    END SUBROUTINE communicate_M_basbas_multiple

    SUBROUTINE register_custom_communicator(a_custom_communicator)
        type(custom_communicator), INTENT(INOUT) :: a_custom_communicator

        INTEGER :: ierror, &
                   group_total, &
                   group_custom_communicator, &
                   comm_custom_communicator, &
                   rank_custom_communicator, &
                   n_tasks_custom_communicator

        call MPI_COMM_GROUP(MPI_COMM_WORLD, group_total, ierror)
        call MPI_GROUP_INCL(group_total, &
                            a_custom_communicator%n_ranks_of_comm_global, &
                            a_custom_communicator%ranks_of_comm_global, &
                            group_custom_communicator, &
                            ierror)

       call MPI_COMM_CREATE(MPI_COMM_WORLD, &
                            group_custom_communicator, &
                            comm_custom_communicator, &
                            ierror)

       call MPI_COMM_SIZE(comm_custom_communicator, &
                          n_tasks_custom_communicator, &
                          ierror)

       call MPI_COMM_RANK(comm_custom_communicator, &
                          rank_custom_communicator, &
                          ierror)


       CALL assign_custom_communicator(group_custom_communicator, &
                                       comm_custom_communicator, &
                                       n_tasks_custom_communicator, &
                                       rank_custom_communicator, &
                                       a_custom_communicator)
    END SUBROUTINE register_custom_communicator

END MODULE cRPA_parallelism
