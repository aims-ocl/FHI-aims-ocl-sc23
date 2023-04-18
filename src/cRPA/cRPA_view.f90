!Functions for output to console or file
MODULE cRPA_view
    USE mpi_tasks
    use localorb_io, only : use_unit

    IMPLICIT NONE

    INTEGER :: ROOT_MESSAGE
    PARAMETER(ROOT_MESSAGE = -1)


    INTERFACE num2str
        MODULE PROCEDURE int2str_array
        MODULE PROCEDURE int2str_single
        MODULE PROCEDURE real2str_single
    END INTERFACE num2str


CONTAINS
    subroutine write_stdout(message, format)

      implicit none

      character(len=*) :: message
      character(len=*), optional :: format

      CALL write_stdout_if_myid(ROOT_MESSAGE, message,format)

    end subroutine write_stdout

    subroutine write_debug(message, format)

      implicit none

      INTEGER :: mpi_result
      character(len=*) :: message
      character(len=*), optional :: format

      CALL write_stdout_if_myid(myid, message,format)
    end subroutine write_debug


    subroutine write_stdout_as_task(message, format)

      implicit none

      character(len=*) :: message
      character(len=*), optional :: format

      CALL write_stdout_if_myid(myid, message,format)

    end subroutine write_stdout_as_task


    subroutine write_stdout_if_myid(tasks_id, message, format)

      implicit none

      INTEGER, INTENT(IN) :: tasks_id
      character(len=*) :: message
      character(len=*), optional :: format

      ! VB: Corrected construct below which crashes the pgi compiler.
      !     The original code may have been technically legal but is
      !     risky because it performs a sophisticated dance where a
      !     simple write statement sufficed.

      !     (I think the problem may have been be a function calling a function in 
      !     a write statement but this is just a guess. I do not
      !     have access to the PGI compiler.)

      if (tasks_id == ROOT_MESSAGE .AND. myid == 0) then
         if (present(format)) then
            write(use_unit, format) trim(message)
         else
            write(use_unit, *) trim(message)
         end if
      end if

      if (tasks_id /= ROOT_MESSAGE .AND. myid == tasks_id) then

         ! format will be ignored in this case. It looks like ironically
         ! the code below ignored a few other conventions, those for the
         ! user-visible output of FHI-aims
         !if (present(format)) then
         !   write(use_unit, ) trim("Task "//num2str(myid)//" :"//message)
         !else
         !   write(use_unit, *) trim("Task "//num2str(myid)//" :"//message)
         !end if

         write(use_unit,'(2X,A,I8,A,A)') "Task ", myid, ":v", trim(message)

      end if

    end subroutine write_stdout_if_myid

    CHARACTER(len=256) FUNCTION int2str_array(input_int)
        INTEGER, INTENT(IN), DIMENSION(:) :: input_int
        write(int2str_array,*) input_int
        int2str_array = trim(int2str_array)
    END FUNCTION

    CHARACTER(len=20) FUNCTION int2str_single(input_int)
        INTEGER, INTENT(IN) :: input_int

        write(int2str_single,*) input_int
        int2str_single = trim(int2str_single)
    END FUNCTION


    CHARACTER(len=40) FUNCTION real2str_single(input_real)
        REAL*8, INTENT(IN) :: input_real

        write(real2str_single,*) input_real
    END FUNCTION

    SUBROUTINE remove_all_spaces(a_str, new_str)
        CHARACTER(len=*),INTENT(IN) :: a_str

        CHARACTER(len=256) :: new_str
        INTEGER :: i,ls1,ls2

        ls1 = len_trim(a_str)
        ls2 = 0
        do i = 1,ls1
           if(a_str(i:i).ne.' ') then
              ls2 = ls2 + 1
              new_str(ls2:ls2) = a_str(i:i)
           endif
        enddo

        do i=ls2+1, 256
            new_str(i:i) =' '
        enddo
    END SUBROUTINE

    SUBROUTINE write_matrix_array_to_file(n_rows,n_columns, n_dim, &
                                          mat, x, matname)
        INTEGER, INTENT(IN) :: n_rows,n_columns, n_dim
        REAL*8, DIMENSION(n_rows,n_columns,n_dim), INTENT(IN) :: mat
        REAL*8, DIMENSION(n_dim), INTENT(IN) :: x
        CHARACTER(len=*), INTENT(IN) :: matname

        INTEGER :: i_row,i_column, i_dim
        CHARACTER(len=256) :: filename, filename_trimmed

        do i_row=1, n_rows
            do i_column=1, n_columns
                filename =   trim(matname) &
                           //trim("_")//trim(num2str(i_row)) &
                           //trim("_")//trim(num2str(i_column))&
                           //trim(".dat")

                CALL remove_all_spaces(filename,filename_trimmed)


                CALL write_two_arrays_to_file(n_dim,x, mat(i_row,i_column,:), filename_trimmed)
            enddo
        enddo

    END SUBROUTINE write_matrix_array_to_file

    SUBROUTINE write_two_arrays_to_file(n_dim, array1,array2, filename)
        INTEGER, INTENT(IN) :: n_dim
        REAL*8, DIMENSION(n_dim), INTENT(IN) :: array1, &
                                                array2
        CHARACTER(len=*), INTENT(IN) :: filename

        integer :: stat,fd,i

        fd = 67
        open(fd, file=filename, iostat=stat)
!        write(use_unit,*) filename
        if(stat == 0) then
            do i = 1,n_dim
                write(fd,*) array1(i),array2(i)
!                write(use_unit,*) array1(i),array2(i)
            enddo
        else
            stop "Could not create file"
        end if

        close(fd)

    END SUBROUTINE write_two_arrays_to_file

!    CHARACTER(len=80) FUNCTION int2str(input_int)
!        INTEGER, INTENT(IN) :: input_int
!
!        write(int2str,*) input_int
!    END FUNCTION


END MODULE cRPA_view
