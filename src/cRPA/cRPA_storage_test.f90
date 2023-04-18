MODULE cRPA_storage_test
    USE basis
    USE geometry
    USE pbc_lists
    USE prodbas
    USE dimensions
    USE cRPA_storage
    USE cRPA_parallelism
    IMPLICIT NONE

CONTAINS
    SUBROUTINE exp_coeffs_storage_test()

        REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: coeff_matrix
        INTEGER :: i_basbas,i_species1,i_species2,n_Rvec, pos
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Rvec

        max_n_basis_sp = 50
        max_n_basbas_sp=5 * max_n_basis_sp


        write(use_unit,*) "Testing storage engine - plain"
        stop "this test is broken"

        n_Rvec = 16

        ALLOCATE(Rvec(3,n_Rvec))
        ALLOCATE(coeff_matrix(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2,n_Rvec))

!        CALL init_exp_coeff_storage(max_n_basbas_sp, (/ 4,7,11 /))

        do i_species1 = 1,10,2
            do i_species2 = 2,10,2

                CALL get_pairwise_coeff_3fn(i_species1,i_species2,n_Rvec,Rvec, &
                                            coeff_matrix,coeff_matrix,.FALSE.)

                if(.NOT. is_exp_coeff_stored(exp_coeffs_local,max_n_basbas_sp)) then
                    write(use_unit,*) "failed to test is_exp_coeff_stored"
                    stop
                endif

                pos = find_exp_coeff_by_basbas(exp_coeffs_local,max_n_basbas_sp)

                do i_basbas=1, max_n_basbas_sp
                    CALL drop_exp_coeff(exp_coeffs_local,i_basbas)
                enddo

            enddo
        enddo

        CALL finalize_exp_coeff_storage(exp_coeffs_local)

        DEALLOCATE(coeff_matrix)
        DEALLOCATE(Rvec)

        write(use_unit,*) "Test storage engine: ok"
    END SUBROUTINE exp_coeffs_storage_test

    SUBROUTINE greens_functions_storage_test()
        type(realspace_vectors) :: vectors, loc_vectors
        REAL*8, DIMENSION(3,3) :: lattice_vector

        type(basisfunc_distribution) :: dist(2)
        type(greens_functions) :: a_green_func

        DOUBLE COMPLEX, DIMENSION(100,100) :: test,test2

        lattice_vector = reshape((/ 1,0,0 ,0,1,0, 0,0,1/), (/3,3/))
        CALL create_all_realspace_vectors_sorted(n_cells_bvk, cell_index_bvk,lattice_vector,vectors)
        CALL get_n_smallest_realspace_vectors(5,n_cells_bvk, cell_index_bvk,lattice_vector,loc_vectors)

        write(use_unit,*) "Start testing greens function storage"

        CALL distribute_basisfunctions(2,n_basis,vectors, dist)
        CALL init_green_func_storage(a_green_func,vectors,1)

        dist%basisfunc_extend = 100
        test = 1
        CALL store_green_func(a_green_func,vectors%vecs(1),100,100,test)
        CALL store_green_func(a_green_func,vectors%vecs(3),100,100,test)
        CALL store_green_func(a_green_func,vectors%vecs(2),100,100,test)


        CALL store_green_func(a_green_func,vectors%vecs(4),100,100,test)

        CALL finalize_green_func_storage(a_green_func)
        CALL free_all_lattice_vectors(vectors)
        CALL free_all_lattice_vectors(loc_vectors)

    END SUBROUTINE greens_functions_storage_test

    SUBROUTINE queue_test()
        type(queue) :: aqueue
        INTEGER :: n_length, n_extend

        write(use_unit,*) "Testing queue functionalities"

        do n_extend = 1,5,2
            do n_length = 2,1002,200

                write(use_unit,*) "Testing queue (length,extend)", n_length, n_extend

                CALL create_queue(n_length, n_extend, aqueue)
                CALL queue_test_reliability(aqueue, n_extend)
                CALL free_queue(aqueue)

            enddo
        enddo

        write(use_unit,*) "Testing queue functionalities: ok"
    END SUBROUTINE queue_test

    SUBROUTINE queue_test_reliability(aqueue,n_extend)
        type(queue), INTENT(INOUT) :: aqueue
        INTEGER, INTENT(IN) :: n_extend

        INTEGER :: i,j,k,l,data(n_extend)

        !enqueue
        do i=1, aqueue%queue_length
            FORALL(j=1:n_extend) data(j) = i*j
            CALL enqueue(aqueue,n_extend, data)
        enddo

        if(can_enqueue(aqueue)) then
            stop "failed! queue should be full."
        endif

        if(.NOT. can_dequeue(aqueue)) then
            stop "failed! queue has data."
        endif

        !dequeue
        do i=1,aqueue%queue_length
            if(.NOT. can_dequeue(aqueue)) then
              stop "failed! queue has data."
            endif

            CALL dequeue(aqueue,data)
            FORALL(j=1:n_extend) data(j) = data(j) - i*j

            if(.NOT. ALL(data==0)) then
               stop "failed! sequential queuing and dequeuing."
            endif

        enddo

        if(can_dequeue(aqueue)) then
            stop "failed! queue should be empty."
        endif

        if(.NOT. can_enqueue(aqueue)) then
            stop "failed! queue should have space."
        endif

        if(aqueue%queue_pos /=1) then
            stop "failed! pos_queue is not equal to 1."
        endif

        !enqueue & dequeue
        k=1
        l=1
        do i=1, aqueue%queue_length
            do j=i,aqueue%queue_length
                data = k
                k=k+1
                CALL enqueue(aqueue,n_extend, data)
            enddo

            do j=i,aqueue%queue_length-1
                CALL dequeue(aqueue,data)
                if(.NOT. ALL(data==l)) then
                   stop "failed! parallel queuing and dequeuing."
                endif
                l=l+1
            enddo
        enddo

        do i=1, aqueue%queue_length
           if(.NOT. can_dequeue(aqueue)) then
              stop "failed! queue has data."
           endif

           CALL dequeue(aqueue,data)
             if(.NOT. ALL(data==l)) then
               stop "failed! parallel queuing and dequeuing in remainer."
             endif
             l=l+1
       enddo
    END SUBROUTINE queue_test_reliability

    SUBROUTINE symmetry_greensfunc_test(greensfunc)
        type(greens_functions) :: greensfunc

        REAL*8 :: max_value, &
                  diff_transpose, &
                  max_diff_transpose(greensfunc%rvecs%n_vecs), &
                  diff_inversion, &
                  max_diff_inversion

        INTEGER :: i_vec, &
                   i_spin, &
                   i_vec_inversion(greensfunc%rvecs%n_vecs), &
                   mpi_result

        CALL write_stdout("Testing symmetry: greensfunc")

        max_value = maxval(abs(greensfunc%g(:,:,:,:)))
        CALL write_stdout("Max value:"//num2str(max_value))

        CALL get_realspace_vectors_index_inversion(greensfunc%rvecs, &
                                                   i_vec_inversion)


!        max_diff_transpose(:) = 0.d0
!        do i_vec = 1, greensfunc%rvecs%n_vecs
!            do i_spin = 1, greensfunc%n_spin
!                diff_transpose = maxval(abs(greensfunc%g(:,:,i_vec,i_spin) -  &
!                                            transpose(greensfunc%g(:,:,i_vec,i_spin))))
!                max_diff_transpose(i_vec) = max(max_diff_transpose(i_vec), diff_transpose)
!            enddo
!
!            CALL write_stdout("Max diff transpose:"//num2str(i_vec)//num2str(i_vec_inversion(i_vec)) &
!                                //num2str(max_diff_transpose(i_vec)))
!        enddo
!
!        max_diff_inversion = 0.d0
!        do i_vec = 1, greensfunc%rvecs%n_vecs
!            do i_spin = 1, greensfunc%n_spin
!                diff_inversion = maxval(abs(greensfunc%g(:,:,i_vec,i_spin) +  &
!                                            greensfunc%g(:,:,i_vec_inversion(i_vec),i_spin)))
!                max_diff_inversion = max(max_diff_inversion, diff_inversion)
!            enddo
!        enddo
!
!        CALL write_stdout("Max diff inversion:"//num2str(max_diff_inversion))

        max_diff_transpose(:) = 0.d0
        do i_vec = 1, greensfunc%rvecs%n_vecs
            do i_spin = 1, greensfunc%n_spin
                diff_transpose = maxval(abs(greensfunc%g(:,:,i_vec,i_spin) -  &
                                            transpose(greensfunc%g(:,:,i_vec_inversion(i_vec),i_spin))))
                max_diff_transpose(i_vec) = max(max_diff_transpose(i_vec), diff_transpose)
            enddo

            if(max_diff_transpose(i_vec) > 1.d-8) then
            CALL write_stdout(">>>Max diff transpose/inversion:"//num2str(i_vec)//num2str(i_vec_inversion(i_vec))&
                              //num2str(max_diff_transpose(i_vec)))
            endif
        enddo



        CALL write_stdout("Passed symmetry test")
    END SUBROUTINE symmetry_greensfunc_test

    SUBROUTINE symmetry_exp_coeffs_test(exp_coeffs, exp_coeffs_offsite)
        type(expansion_coefficients) :: exp_coeffs, exp_coeffs_offsite

        REAL*8 :: max_value, &
                  diff_transpose, &
                  max_diff_transpose(exp_coeffs%rvecs%n_vecs), &
                  diff_inversion, &
                  max_diff_inversion

        INTEGER :: i_vec, &
                   i_basbas, &
                   i_vec_inversion(exp_coeffs%rvecs%n_vecs), &
                   mpi_result


        CALL MPI_BARRIER(mpi_comm_global,mpi_result)

        CALL write_stdout("Testing symmetry: exp_coeffs")

        if(exp_coeffs%basis_dist%basisfunc_start &
             /= &
           exp_coeffs%basis_dist%my_basis_start) then

           CALL write_stdout("Cannot test: exp_coeffs is not square")
           CALL write_stdout("passed")
        endif


!        max_diff_transpose(:) = 0.d0
!        do i_vec = 1, exp_coeffs%rvecs%n_vecs
!
!            do i_basbas = 1, n_basbas
!                diff_transpose = maxval(abs(exp_coeffs%coeff(i_basbas)%c(:,:,i_vec) -  &
!                                            transpose(exp_coeffs_offsite%coeff(i_basbas)%c(:,:,i_vec))))
!                max_diff_transpose(i_vec) = max(max_diff_transpose(i_vec), diff_transpose)
!            enddo
!
!            CALL write_stdout("ec Max diff transpose:"//num2str(i_vec)//num2str(max_diff_transpose(i_vec)))
!        enddo
!
!
!
        CALL get_realspace_vectors_index_inversion(exp_coeffs%rvecs, &
                                                   i_vec_inversion)
!
!
!        max_diff_inversion = 0.d0
!        do i_vec = 1, exp_coeffs%rvecs%n_vecs
!            do i_basbas = 1, n_basbas
!                diff_inversion = maxval(abs(exp_coeffs%coeff(i_basbas)%c(:,:,i_vec) -  &
!                                            exp_coeffs%coeff(i_basbas)%c(:,:,i_vec_inversion(i_vec))))
!                max_diff_inversion = max(max_diff_inversion, diff_inversion)
!            enddo
!        enddo
!
!        CALL write_stdout("ec Max diff inversion:"//num2str(max_diff_inversion))

        max_diff_transpose(:) = 0.d0
        do i_vec = 1, exp_coeffs%rvecs%n_vecs
            do i_basbas = 1, n_basbas
!                diff_transpose = maxval(abs(exp_coeffs%coeff(i_basbas)%c(:,:,i_vec) -  &
!                                            transpose(exp_coeffs_offsite%coeff(i_basbas)%c(:,:,i_vec_inversion(i_vec)))))
                max_diff_transpose(i_vec) = max(max_diff_transpose(i_vec), diff_transpose)
            enddo

            if(max_diff_transpose(i_vec) > 1.d-8) then
                CALL write_stdout("ec Max diff transpose/inversion:"//num2str(i_vec)//num2str(max_diff_transpose(i_vec)))
            endif
        enddo


        CALL write_stdout("Passed symmetry test")

        CALL MPI_BARRIER(mpi_comm_global,mpi_result)

    END SUBROUTINE symmetry_exp_coeffs_test

END MODULE cRPA_storage_test
