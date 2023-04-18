!Storage of greensfunction g_0
MODULE cRPA_storage_greensfunc
    USE dimensions


    USE cRPA_view
    USE cRPA_parallelism_storage
    USE cRPA_flow_realspace_vector
    USE cRPA_storage_queue
    USE cRPA_calculation_integration
    USE cRPA_calculation_fouriertrans
    USE cRPA_flow_adaptive_grid
    USE cRPA_storage_norms
    IMPLICIT NONE

    type greens_functions

        !(atom,atom,store_id)
        DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: g
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: g_norm_row, &
                                               g_norm_column
        type(realspace_vectors) :: rvecs
        INTEGER :: n_rvecs, n_spin
        LOGICAL is_greens_func_storage_inited
    end type

    type(greens_functions) :: green_funcs_local, &
                              green_funcs_local_minus_tau

CONTAINS

    SUBROUTINE init_green_func_storage(this,all_vectors,n_spin)
        type(greens_functions), INTENT(INOUT) :: this
        type(realspace_vectors), INTENT(IN) :: all_vectors
        INTEGER,INTENT(IN) :: n_spin

        INTEGER :: n_my_basis_size

        if(this%is_greens_func_storage_inited) then
            write(use_unit,*) "Greens function storage already inited!"
            stop
        end if

!        CALL write_stdout("Initializing greens function storage")

!TODO check allocation

        n_my_basis_size = n_basis
        this%n_rvecs = all_vectors%n_vecs
        this%n_spin = n_spin

        ALLOCATE(this%g(n_my_basis_size,n_basis,this%n_rvecs,n_spin))
        ALLOCATE(this%g_norm_row(this%n_rvecs,this%n_spin))
        ALLOCATE(this%g_norm_column(this%n_rvecs,this%n_spin))

        CALL clone_realspace_vectors(all_vectors, this%rvecs)

        this%is_greens_func_storage_inited = .TRUE.

    END SUBROUTINE init_green_func_storage

    SUBROUTINE store_green_func(this,rvec,dim_my_basis,dim_basis,matrix_green_func)
        type(greens_functions), INTENT(INOUT) :: this
        type(realspace_vector), INTENT(IN) :: rvec
        INTEGER, INTENT(IN) :: dim_my_basis,dim_basis
        DOUBLE COMPLEX, DIMENSION(dim_my_basis,dim_basis,this%n_spin), INTENT(IN) :: matrix_green_func

        INTEGER :: pos,i_spin

        pos = rvec%i_cell_bvk

        this%g(:,:,pos,:)=matrix_green_func(:,:,:)

        do i_spin = 1,this%n_spin
            CALL get_norm_max_row(dim_my_basis,n_basis, abs(this%g(:,:,pos,i_spin)), &
                                  this%g_norm_row(pos,i_spin))

            CALL get_norm_max_column(dim_my_basis,n_basis, abs(this%g(:,:,pos,i_spin)), &
                                     this%g_norm_column(pos,i_spin))
        enddo
    END SUBROUTINE store_green_func

    SUBROUTINE set_green_func_norms(this)
        type(greens_functions), INTENT(INOUT) :: this

        INTEGER :: i_vec, n_my_basis, i_spin

        n_my_basis = n_basis

        do i_spin = 1, this%n_spin
            do i_vec = 1, this%rvecs%n_vecs
                CALL get_norm_max_row(n_my_basis,n_basis, abs(this%g(:,:,i_vec,i_spin)), &
                                      this%g_norm_row(i_vec,i_spin))

                CALL get_norm_max_column(n_my_basis,n_basis, abs(this%g(:,:,i_vec,i_spin)), &
                                         this%g_norm_column(i_vec,i_spin))
            enddo
        enddo
    END SUBROUTINE set_green_func_norms

    SUBROUTINE get_green_func_norm_mean(this, norm_mean)
        type(greens_functions), INTENT(INOUT) :: this
        REAL*8, INTENT(OUT) :: norm_mean

        INTEGER :: i_vec, i_spin, mpi_result
        REAL*8 :: sum_norm, my_norm_mean

        my_norm_mean = 0.d0
        do i_spin = 1,this%n_spin
            do i_vec = 1, this%rvecs%n_vecs
                CALL get_norm_sum_column(ubound(this%g,1), &
                                         ubound(this%g,2), &
                                         abs(this%g(:,:,i_vec,i_spin)), my_norm_mean)
                my_norm_mean = my_norm_mean +&
                               sum_norm
            enddo
        enddo

        CALL MPI_ALLREDUCE(my_norm_mean, norm_mean, 1, &
                           MPI_REAL8, MPI_SUM, mpi_comm_global, mpi_result)


        norm_mean = norm_mean / dble(n_basbas)

    END SUBROUTINE get_green_func_norm_mean

    SUBROUTINE get_green_func_norm_max_column(this, norm_max_row)
        type(greens_functions), INTENT(INOUT) :: this
        REAL*8, INTENT(OUT) :: norm_max_row

        INTEGER :: i_vec, i_spin, mpi_result
        REAL*8 :: max_norm, my_norm_max_row

        my_norm_max_row = 0.d0
        do i_spin = 1,this%n_spin
            do i_vec = 1, this%rvecs%n_vecs
    !            CALL get_norm_max_column(ubound(this%g,1), &
    !                                     ubound(this%g,2), &
    !                                     this%g(:,:,i_vec), max_norm)
                max_norm = chksum_2d_matrix(abs(this%g(:,:,i_vec,i_spin)))
                my_norm_max_row = my_norm_max_row +&
                                  max_norm
            enddo
        enddo

        CALL MPI_ALLREDUCE(my_norm_max_row, norm_max_row, 1, &
                           MPI_REAL8, MPI_SUM, mpi_comm_global, mpi_result)


        norm_max_row = norm_max_row / dble(this%rvecs%n_vecs)

    END SUBROUTINE get_green_func_norm_max_column

    REAL*8 FUNCTION chksum_green_func(this)
        type(greens_functions), INTENT(INOUT) :: this
        REAL*8 :: chksum
        INTEGER :: i_spin

        chksum = 0.d0
        do i_spin = 1,1!this%n_spin
            chksum = chksum + &
                     chksum_3d_matrix(abs(this%g(:,:,:,i_spin)))
        enddo

        chksum_green_func = chksum
    END FUNCTION chksum_green_func


    SUBROUTINE print_green_func(this)
        type(greens_functions), INTENT(INOUT) :: this

        REAL*8 :: my_min_norm_column, &
                  my_min_norm_row, &
                  my_max_norm_column, &
                  my_max_norm_row, &
                  min_norm_column, &
                  min_norm_row, &
                  max_norm_column, &
                  max_norm_row, &
                  mean_norm_column


        my_max_norm_column = maxval(this%g_norm_column(:,:))
        my_min_norm_column = minval(this%g_norm_column(:,:))
        my_max_norm_row = maxval(this%g_norm_row(:,:))
        my_min_norm_row = minval(this%g_norm_row(:,:))

        CALL get_max_over_all_cpus(my_max_norm_column,max_norm_column)
        CALL get_min_over_all_cpus(my_min_norm_column,min_norm_column)

        CALL get_max_over_all_cpus(my_max_norm_row,max_norm_row)
        CALL get_min_over_all_cpus(my_min_norm_row,min_norm_row)

        CALL get_green_func_norm_mean(this,mean_norm_column)

        CALL write_stdout("Greens function info")
        CALL write_stdout("Norm column (max,min,mean)"//num2str(max_norm_column) &
                                                 //num2str(min_norm_column)&
                                                 //num2str(mean_norm_column))
        CALL write_stdout("Row column (max,min)"//num2str(max_norm_row) &
                                                 //num2str(min_norm_row))
    END SUBROUTINE print_green_func

    SUBROUTINE debug_print_green_func(this)
        type(greens_functions), INTENT(INOUT) :: this

        INTEGER :: i,j,i_spin,i_vec

        do i_spin =1, n_spin
            do i_vec = 1, this%rvecs%n_vecs
                do j =1, ubound(this%g,2)
                    do i =1, ubound(this%g,1)
                        write(use_unit,*) myid,i_spin,i_vec,i,j, this%g(i,j,i_vec,i_spin)
                    enddo
                enddo
            enddo
        enddo

        CALL write_stdout("Total chksum:"//num2str(chksum_green_func(this)))
    END SUBROUTINE debug_print_green_func

    SUBROUTINE finalize_green_func_storage(this)
        type(greens_functions), INTENT(INOUT) :: this

        if(this%is_greens_func_storage_inited) then

            DEALLOCATE(this%g)
            DEALLOCATE(this%g_norm_row)
            DEALLOCATE(this%g_norm_column)

            CALL free_all_lattice_vectors(this%rvecs)

            this%is_greens_func_storage_inited = .FALSE.
        end if
    END SUBROUTINE finalize_green_func_storage

END MODULE cRPA_storage_greensfunc
