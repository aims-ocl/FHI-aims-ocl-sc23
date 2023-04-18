!Implements basic matrix norms
MODULE cRPA_storage_norms
    USE cRPA_view
    implicit none

    real*8, external :: dasum

CONTAINS

    SUBROUTINE get_norm_max_row(dim_rows, dim_columns, mat, norm_max_row)
        INTEGER, INTENT(IN) :: dim_rows, dim_columns
        REAL*8,DIMENSION(dim_rows,dim_columns), INTENT(IN) :: mat
        REAL*8, INTENT(OUT) :: norm_max_row

        CALL get_norm_max_column(dim_columns,dim_rows,transpose(mat),norm_max_row)

    END SUBROUTINE get_norm_max_row

    SUBROUTINE get_norm_max_column(dim_rows, dim_columns, mat, norm_max_column)
        INTEGER, INTENT(IN) :: dim_rows, dim_columns
        REAL*8,DIMENSION(dim_rows,dim_columns), INTENT(IN) :: mat
        REAL*8, INTENT(OUT) :: norm_max_column

        INTEGER :: i_column

        norm_max_column = 0
        do i_column = 1,dim_columns
           norm_max_column = max(norm_max_column, &
                                 !sum(abs(mat(:,i_column))))
                                 dasum(dim_rows,mat(:,i_column),1))
        enddo

    END SUBROUTINE get_norm_max_column

    SUBROUTINE get_norm_sum_row(dim_rows, dim_columns, mat, norm_sum_row)
        INTEGER, INTENT(IN) :: dim_rows, dim_columns
        REAL*8,DIMENSION(dim_rows,dim_columns), INTENT(IN) :: mat
        REAL*8, INTENT(OUT) :: norm_sum_row

        CALL get_norm_sum_column(dim_columns,dim_rows,transpose(mat),norm_sum_row)

    END SUBROUTINE get_norm_sum_row

    SUBROUTINE get_norm_sum_column(dim_rows, dim_columns, mat, norm_sum_column)
        INTEGER, INTENT(IN) :: dim_rows, dim_columns
        REAL*8,DIMENSION(dim_rows,dim_columns), INTENT(IN) :: mat
        REAL*8, INTENT(OUT) :: norm_sum_column

        INTEGER :: i_column

        norm_sum_column = 0
        do i_column = 1,dim_columns
           norm_sum_column = norm_sum_column + &
                             dasum(dim_rows,mat(:,i_column),1)
        enddo

    END SUBROUTINE get_norm_sum_column


    SUBROUTINE get_min_over_all_cpus(my_min, total_min)
       REAL*8,INTENT(IN) :: my_min
       REAL*8,INTENT(OUT):: total_min

       INTEGER :: mpi_result

       CALL MPI_ALLREDUCE(my_min, total_min, 1, &
                          MPI_REAL8, MPI_MIN, mpi_comm_global , mpi_result)
    END SUBROUTINE get_min_over_all_cpus

    SUBROUTINE get_max_over_all_cpus(my_max, total_max)
        REAL*8,INTENT(IN) :: my_max
        REAL*8,INTENT(OUT):: total_max

        INTEGER :: mpi_result

        CALL MPI_ALLREDUCE(my_max, total_max, 1, &
                           MPI_REAL8, MPI_MAX, mpi_comm_global , mpi_result)

    END SUBROUTINE get_max_over_all_cpus


    REAL*8 FUNCTION chksum_2d_matrix(matrix)
        REAL*8,DIMENSION(:,:), INTENT(IN) :: matrix

        INTEGER :: n_dim,i,j

        n_dim=ubound(matrix,1)* &
              ubound(matrix,2)


!        CALL write_debug("n_dim:"//num2str(n_dim))

        chksum_2d_matrix=0.d0
        do j=1, ubound(matrix,2)
!                CALL write_debug("stored"//num2str(response_func%element(1,j,k)))

           do i=1, ubound(matrix,1)
                   chksum_2d_matrix = chksum_2d_matrix + abs(matrix(i,j))
           enddo
        enddo
    END FUNCTION chksum_2d_matrix



    REAL*8 FUNCTION chksum_3d_matrix(matrix)
        REAL*8,DIMENSION(:,:,:), INTENT(IN) :: matrix

        INTEGER :: n_dim,i,j,k

        n_dim=ubound(matrix,1)* &
              ubound(matrix,2)* &
              ubound(matrix,3)


!        CALL write_debug("n_dim:"//num2str(n_dim))

        chksum_3d_matrix=0.d0

        do k=1, ubound(matrix,3)
           do j=1, ubound(matrix,2)
!                CALL write_debug("stored"//num2str(response_func%element(1,j,k)))

            do i=1, ubound(matrix,1)
                   chksum_3d_matrix = chksum_3d_matrix + &
                                              abs(matrix(i,j,k))

              !if(isnan(abs(matrix(i,j,k)))) then
              !  CALL write_debug("nan"//num2str(i)//num2str(j)//num2str(k)//num2str(matrix(i,j,k)))
              !endif
            enddo
           enddo
        enddo
    END FUNCTION chksum_3d_matrix

END MODULE cRPA_storage_norms
