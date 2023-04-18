MODULE ipc

    CONTAINS

    FUNCTION ipc_start_transaction(transaction_name) RESULT(is_transaction_active)
        CHARACTER(len=*) :: transaction_name
        LOGICAL :: is_transaction_active

        is_transaction_active = .FALSE.
    END FUNCTION

    SUBROUTINE ipc_write_string(string_to_write)
        CHARACTER(LEN=*) :: string_to_write
    END SUBROUTINE

    SUBROUTINE ipc_end_session()
        return
    END SUBROUTINE

    subroutine ipc_read_array(a, b)
      character(*), intent(in) :: a
      real*8 :: b(*)
    end subroutine ipc_read_array

    subroutine ipc_write_array(a, b)
      character(*), intent(in) :: a
      real*8 :: b(*)
    end subroutine ipc_write_array
END MODULE ipc

