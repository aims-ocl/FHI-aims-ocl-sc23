MODULE cRPA_parallelism_test
    USE cRPA_flow_realspace_vector
    USE cRPA_parallelism
    implicit none
CONTAINS

    SUBROUTINE communicate_three_dimensional_matrices_test()
        INTEGER, ALLOCATABLE, DIMENSION(:) :: map_to_cpu_list
        type(message_matrix_send), ALLOCATABLE, DIMENSION(:) :: mat
        INTEGER :: dim1(n_tasks),dim2(n_tasks),dim3(n_tasks), dim_index(n_tasks)
        INTEGER :: i,j,k,l
        PROCEDURE(on_recv_template),POINTER :: do_on_recv


        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: aux_mat

!        CALL initialize_mpi()

        do i=1, n_tasks
            dim1(i) = 100-i-10
            dim2(i) = 100+i
            dim3(i) = 100-i

            dim_index(i) = 40-i
        enddo

        write(use_unit,*) dim1
        write(use_unit,*) dim2
        write(use_unit,*) dim3
        write(use_unit,*) dim_index

!TODO check allocation
        ALLOCATE(mat(dim_index(myid+1)))
        ALLOCATE(aux_mat(dim1(myid+1),dim2(myid+1),dim3(myid+1)))
        ALLOCATE(map_to_cpu_list(dim_index(myid+1)))

        do i = 1, dim_index(myid+1)
            map_to_cpu_list(i) = mod(i,n_tasks)

            do j=1,dim1(myid+1)
                do k=1, dim2(myid+1)
                    do l=1, dim3(myid+1)
                        aux_mat(j,k,l) = myid+j+k+l
                    enddo
                end do
            end do

            CALL create_message_matrix_send_from_matrix(dim1(myid+1),dim2(myid+1),dim3(myid+1),aux_mat,mat(i))
        enddo


!        write(use_unit,*), map_to_cpu_list

        do_on_recv => communicate_three_dimensional_matrices_test_on_recv
        CALL communicate_three_dimensional_matrices(dim_index(myid+1),mat, map_to_cpu_list,do_on_recv)

! validate

!            write(use_unit,*) "", myid,": from ", res(i)%source,"  ",res(i)%data(1)
        write(use_unit,*) myid, "all passed"


!        do i=1, dim_index(myid+1)
!            CALL free_message_matrix_send(mat(i))
!        enddo

        DEALLOCATE(map_to_cpu_list)
        DEALLOCATE(aux_mat)
        DEALLOCATE(mat)

     END SUBROUTINE communicate_three_dimensional_matrices_test

    SUBROUTINE communicate_three_dimensional_matrices_test_on_recv(data_recv)
        type(message_matrix_recv), INTENT(IN) :: data_recv

        integer ::j,k,l,dim1,dim2,dim3

        dim1 = data_recv%dim1
        dim2 = data_recv%dim2
        dim3 = data_recv%dim3

        do j=1,dim1
            do k=1,dim2
                do l=1,dim3
                    if(data_recv%data(j,k,l)/=data_recv%source+j+k+l) then
                        stop "!communicate_three_dimensional_matrices_test failed! (validation)"
                    endif
                enddo
            enddo
        enddo
    END SUBROUTINE communicate_three_dimensional_matrices_test_on_recv

    SUBROUTINE generate_my_basbas_contribution_test()
        type(basbas_distribution) :: worklist
        INTEGER, ALLOCATABLE, DIMENSION(:) :: basbas_to_cpu

        INTEGER :: num_cpu,num_aux

        print*, "Starting generate_my_basbas_contribution_test"

        do num_cpu = 60,60!,480, 50
            do num_aux = num_cpu, 5000, 500
                ALLOCATE(basbas_to_cpu(num_aux))
!                CALL generate_my_basbas_contribution(num_aux,num_cpu,worklist,basbas_to_cpu)
                DEALLOCATE(basbas_to_cpu)
            enddo
        enddo
    END SUBROUTINE generate_my_basbas_contribution_test

!TODO add 2d
    SUBROUTINE distribute_basbas_test()
        INTEGER :: num_aux,num_cpu,i,j
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_matrix

        print*, "Starting distribute_aux_test"

        do num_cpu = 60,380, 70
            do num_aux = num_cpu, 3000, 1000
                ALLOCATE (test_matrix(num_aux,num_aux))

                test_matrix(:,:) = 1
                CALL distribute_basbas_1D_test_matrix(num_aux,num_cpu,test_matrix)

                DEALLOCATE (test_matrix)
           end do
        end do

        print*, "Ended distribute_aux_test - all passed"

    END SUBROUTINE distribute_basbas_test

    SUBROUTINE distribute_basbas_1D_test_matrix(num_aux,num_cpu,test_matrix)
        INTEGER, INTENT(IN) :: num_aux,num_cpu
        INTEGER, DIMENSION(num_aux,num_aux), INTENT(INOUT) :: test_matrix

        INTEGER i,j
        type(basbas_distribution), ALLOCATABLE, DIMENSION(:) :: worklist,worklist_total
        type(point) :: element

        print*, "Testing distribute_aux with ", num_cpu, " CPUs and ", num_aux," num_aux"
        CALL distribute_basbas(num_aux,1,num_cpu,worklist,worklist_total)

        do i = 1, num_cpu
!            CALL print_matrixrow(worklist(i))

            do j = 1, worklist(i)%dim_row
                element = worklist(i)%element(j)
                if(.NOT. is_point_zero(element)) then
                    test_matrix(element%x,element%y) = test_matrix(element%x,element%y)-1
                endif
            end do
        end do

        CALL distribute_basbas_test_validate_test_matrix(num_aux,test_matrix)

        DEALLOCATE(worklist)
        DEALLOCATE(worklist_total)
!        if(ANY(basbas_to_cpu==-1)) then
!            write(use_unit,*) "Failed to test basbas to cpu!!"
!            stop
!        endif

    END SUBROUTINE distribute_basbas_1D_test_matrix

    SUBROUTINE distribute_basbas_test_validate_test_matrix(num_aux,test_matrix)
        INTEGER, INTENT(IN) :: num_aux
        INTEGER, DIMENSION(num_aux,num_aux), INTENT(IN) :: test_matrix

        INTEGER :: i,j

        if(ANY(test_matrix/=0)) then
            do i = 1,num_aux
!                print*, test_matrix(i,:)
            end do

            do i = 1,num_aux
                do j = 1,num_aux
                    if(test_matrix(i,j)/=0) then
                        !print*, "i,j :",i,j,test_matrix(i,j)
                    end if
                end do
            end do


            print*, "Failed!!!"
            stop
        end if
    END SUBROUTINE distribute_basbas_test_validate_test_matrix

    SUBROUTINE distribute_basisfunctions_test()
        INTEGER :: n_cpus,n_basisfuncs, i
        type(basisfunc_distribution), ALLOCATABLE,DIMENSION(:) :: dist
        INTEGER, ALLOCATABLE, DIMENSION(:) :: test_funcs
        type(realspace_vectors) :: rvecs
        REAL*8, DIMENSION(3,3) :: lattice_vector

        lattice_vector = reshape((/ 1,0,0 ,0,1,0, 0,0,1/), (/3,3/))

        CALL create_all_realspace_vectors_sorted(n_cells_bvk, cell_index_bvk,lattice_vector,rvecs)


        write(use_unit,*) "Testing distribute_basisfunctions"

!TODO check allocation
        do n_cpus=1,1010,12
            do n_basisfuncs = 100,10100,500


                ALLOCATE(dist(n_cpus))
                ALLOCATE(test_funcs(n_basisfuncs))

                write(use_unit,*) "Testing distribute basisfunctions: n_cpus=",n_cpus, "n_basis", n_basisfuncs


                test_funcs = 1
                CALL distribute_basisfunctions(n_cpus,n_basisfuncs,rvecs,dist)


                do i = 1, n_cpus
!                    write(use_unit,*), i,dist(i)%basisfunc_start,dist(i)%basisfunc_end
                    test_funcs(dist(i)%basisfunc_start: &
                               dist(i)%basisfunc_end) = test_funcs(dist(i)%basisfunc_start : &
                                                                    dist(i)%basisfunc_end) - 1
                end do


                if(ANY(test_funcs /= 0)) then
                    stop "failed to test distribute_basisfuncs"
                endif

                DEALLOCATE(test_funcs)
                DEALLOCATE(dist)
            enddo
        enddo

        write(use_unit,*) "Testing distribute basisfuncs: ok"
    END SUBROUTINE distribute_basisfunctions_test

    SUBROUTINE distribute_basbas_row_blocks_test()

        INTEGER :: i,num_aux, num_blocks
        type(basbas_distribution),ALLOCATABLE, DIMENSION(:) :: basbas_dist
        INTEGER, ALLOCATABLE, DIMENSION(:) :: test_vector

        write(use_unit,*) "Starting test: distribute_basbas_row_block_test"

        do i=1,40,10
            num_aux = i**2
            write(use_unit,*) "num_aux", num_aux
            do num_blocks=1,num_aux,ceiling(real(num_aux)/10d0)
                CALL distribute_basbas_row_blocks(num_aux, num_blocks,basbas_dist)
                ALLOCATE(test_vector(num_aux))
                test_vector(:)=1

                CALL distribute_basbas_row_blocks_test_validate()

                DEALLOCATE(basbas_dist)
                DEALLOCATE(test_vector)
            enddo
        enddo

        write(use_unit,*) "Passed test: distribute_basbas_row_block_test"

    CONTAINS
        SUBROUTINE distribute_basbas_row_blocks_test_validate()
            INTEGER :: i_block,i_row_start,i_row_end
            do i_block = 1,num_blocks
                i_row_start = basbas_dist(i_block)%n_offset_rows + 1
                i_row_end = (i_row_start -1) + basbas_dist(i_block)%n_rows

                test_vector(i_row_start:i_row_end) = test_vector(i_row_start:i_row_end)-1

                if(ANY(test_vector<0)) then
                  write(use_unit,*) "block ",i_block, " caused negative test vector"
                  write(use_unit,*) "num_aux", num_aux
                  write(use_unit,*) "num_blocks", num_blocks
                  write(use_unit,*) "i_row_start ",i_row_start
                  write(use_unit,*) "i_row_end ",i_row_end
                  stop
                endif
            enddo

            if(ANY(test_vector/=0)) then
              write(use_unit,*) "num_aux", num_aux
              write(use_unit,*) "num_blocks", num_blocks
              stop "failed! test_vector is not zero vector!"
            endif

        END SUBROUTINE distribute_basbas_row_blocks_test_validate
    END SUBROUTINE distribute_basbas_row_blocks_test

!test is NOT perfect!
    SUBROUTINE get_basbas_row_blocks_test()

       INTEGER :: n_tasks_global,myid_global, n_blocks

       INTEGER :: n_rank_block
       INTEGER,ALLOCATABLE, DIMENSION(:) :: rank_block, test_vector

       write(use_unit,*) "Starting test: get_basbas_row_blocks_test"

       do n_tasks_global=1,200,12
         write(use_unit,*) "n_tasks_global", n_tasks_global
         do n_blocks=1,n_tasks_global,10

            ALLOCATE(test_vector(n_tasks_global))
            test_vector(:)=-100

            write(use_unit,*) "n_tasks_global", n_tasks_global,"n_blocks", n_blocks
            do myid_global=0,n_tasks_global-1
                 CALL get_cpus_basbas_row_block(n_tasks_global,myid_global, &
                                                n_blocks, &
                                                n_rank_block,rank_block)

                 CALL get_cpus_basbas_row_block_test_update()

                 DEALLOCATE(rank_block)
             enddo

             CALL get_cpus_basbas_row_block_test_validate()
             DEALLOCATE(test_vector)
          enddo
       enddo

       write(use_unit,*) "Passed test: get_basbas_row_blocks_test"
    CONTAINS
        SUBROUTINE get_cpus_basbas_row_block_test_update()
            INTEGER :: i

            do i=1,n_rank_block
                if(test_vector(rank_block(i)+1)==-100) then
                    test_vector(rank_block(i)+1) = n_rank_block -1
                else
                    test_vector(rank_block(i)+1) = test_vector(rank_block(i)+1) - 1
                endif
            enddo
        END SUBROUTINE get_cpus_basbas_row_block_test_update

        SUBROUTINE get_cpus_basbas_row_block_test_validate()
            if(ANY(test_vector/=0)) then
                write(use_unit,*) "test_vector"
                write(use_unit,*) test_vector
                stop "failed"
            endif
        END SUBROUTINE get_cpus_basbas_row_block_test_validate
    END SUBROUTINE get_basbas_row_blocks_test
END MODULE cRPA_parallelism_test
