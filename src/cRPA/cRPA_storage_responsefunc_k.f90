!Storage of response function if calculated in k-space
MODULE cRPA_storage_responsefunc_k
    USE cRPA_parallelism_storage

    type result_last_calculation
        INTEGER :: i_k_point_start
        INTEGER :: i_k_point_end
        INTEGER :: n_frequencies
        REAL*8, ALLOCATABLE, DIMENSION(:) :: frequencies
        DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: data
        type(basbas_distribution) :: basbas_dist
    end type result_last_calculation

CONTAINS

    SUBROUTINE init_result_last_calculation(this, i_k_point_start, n_points , &
                                     n_frequencies, frequencies, &
                                     work_dist)
       type(result_last_calculation) :: this
       INTEGER, INTENT(IN) :: i_k_point_start, n_points, n_frequencies
       REAL*8, DIMENSION(n_frequencies), INTENT(IN) :: frequencies
       type(basbas_distribution) :: work_dist
       INTEGER :: allocation_info
       CHARACTER(*), PARAMETER :: funcname = 'init_result_last_calculation'


       this%i_k_point_start = i_k_point_start
       this%i_k_point_end = i_k_point_start + n_points - 1
       this%n_frequencies = n_frequencies
       ALLOCATE(this%frequencies(n_frequencies), stat=allocation_info)
       call check_allocation(allocation_info, 'frequencies', funcname)

       this%frequencies(:) = frequencies(:)
       this%basbas_dist = work_dist
       ALLOCATE(this%data(n_frequencies, n_points, &
                          this%basbas_dist%n_columns, &
                          this%basbas_dist%n_rows), stat=allocation_info)
       call check_allocation(allocation_info, 'data', funcname)
       this%data(:,:,:,:) = 0.d0
    END SUBROUTINE init_result_last_calculation

    LOGICAL FUNCTION has_result_last_calculation(this,i_k_point, frequency)
       type(result_last_calculation) :: this
       INTEGER,INTENT(IN) :: i_k_point
       REAL*8, INTENT(IN) :: frequency

       INTEGER :: i_freq

       if(this%i_k_point_start <= i_k_point &
          .AND. &
          i_k_point <= this%i_k_point_end) then

          do i_freq = 1,this%n_frequencies
              if(abs(this%frequencies(i_freq)-frequency) < 1.d-10) then
                  has_result_last_calculation = .TRUE.
                  return
              endif
          enddo
       endif

       has_result_last_calculation = .FALSE.
    END FUNCTION

    SUBROUTINE get_result_last_calculation_index(this,i_k_point, frequency, &
                                                 i_k_point_index, &
                                                 i_frequency_index)
       type(result_last_calculation) :: this
       INTEGER,INTENT(IN) :: i_k_point
       REAL*8, INTENT(IN) :: frequency
       INTEGER,INTENT(OUT) :: i_k_point_index
       INTEGER, INTENT(OUT) :: i_frequency_index

       INTEGER :: i_freq

       i_k_point_index = i_k_point - this%i_k_point_start +1

       do i_freq = 1,this%n_frequencies
         if(abs(this%frequencies(i_freq)-frequency) < 1.d-10) then
                i_frequency_index = i_freq
            EXIT
         endif
       enddo

    END SUBROUTINE get_result_last_calculation_index

    SUBROUTINE get_result(this, i_k_point,frequency, n_basbas, data_out)
        type(result_last_calculation) :: this
        INTEGER, INTENT(IN) :: i_k_point
        REAL*8, INTENT(IN) :: frequency
        INTEGER, INTENT(IN) :: n_basbas
        DOUBLE COMPLEX, DIMENSION(n_basbas,n_basbas),INTENT(OUT) :: data_out

        INTEGER :: i_k_point_index, i_freq_index, mpierr, i_row
        DOUBLE COMPLEX, DIMENSION(n_basbas,n_basbas) :: my_part

        CALL get_result_last_calculation_index(this, i_k_point,frequency, &
                                               i_k_point_index, i_freq_index)
        my_part(:,:) = DCMPLX(0.d0,0.d0)
        do i_row = 1,this%basbas_dist%n_rows
            my_part(this%basbas_dist%n_offset_rows + i_row,:) = &
                                this%data(i_freq_index,i_k_point_index,:,i_row)
        enddo

        CALL MPI_ALLREDUCE(my_part, &
                           data_out, size(my_part), &
                           MPI_COMPLEX16, MPI_SUM, mpi_comm_global , mpierr)
    END SUBROUTINE get_result

    SUBROUTINE free_result_last_calculation(this)
       type(result_last_calculation) :: this
       DEALLOCATE(this%frequencies)
       DEALLOCATE(this%data)
    END SUBROUTINE free_result_last_calculation


END MODULE cRPA_storage_responsefunc_k
