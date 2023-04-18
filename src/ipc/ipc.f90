! Fortran convenience interface for the IPC code.
! The real implememention is found in the IPC C source files.
MODULE ipc
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    PRIVATE

    ! Make fortran ipc functions public.
    PUBLIC :: ipc_start_aims,        &
              ipc_start_transaction, &
              ipc_write_string,      &
              ipc_read_string,       &
              ipc_write_array,       & 
              ipc_read_array,        &
              ipc_end_session

!--------------------------------------------------------------------------------------------------
! Import all necessary IPC functions
!--------------------------------------------------------------------------------------------------
    INTERFACE
        INTEGER (C_INT) FUNCTION startTransaction(cptr_transactionname) BIND(C, NAME="startTransaction")
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            TYPE(C_PTR), INTENT(IN), VALUE :: cptr_transactionname
        END FUNCTION
    END INTERFACE    

    INTERFACE
        TYPE(C_PTR) FUNCTION createArrayPacketDouble(i_dimensions, i_dimension_lengths, data_double) BIND(C, NAME="createArrayPacketDouble")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                INTEGER(C_INT), VALUE :: i_dimensions
                TYPE(C_PTR), VALUE :: i_dimension_lengths
                TYPE(C_PTR), VALUE :: data_double
        END FUNCTION
    END INTERFACE

    INTERFACE
          TYPE(C_PTR) FUNCTION createArrayPacketInteger(i_dimensions, i_dimension_lengths, data_integer) BIND(C, NAME="createArrayPacketInteger")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                INTEGER(C_INT), VALUE :: i_dimensions
                TYPE(C_PTR), VALUE :: i_dimension_lengths
                TYPE(C_PTR), VALUE :: data_integer
          END FUNCTION
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION isArrayPacketDoubles(packet) BIND(C, NAME="isArrayPacketDoubles")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet                
          END FUNCTION
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION isArrayPacketIntegers(packet) BIND(C, NAME="isArrayPacketIntegers")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet
          END FUNCTION
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION isArrayPacketString(packet) BIND(C, NAME="isArrayPacketString")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet
          END FUNCTION
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION arrayPacketDimensionSize(packet, i_dimension) BIND(C, NAME="arrayPacketDimensionSize")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet
                INTEGER(C_INT),VALUE :: i_dimension
          END FUNCTION
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION copyArrayPacketData(packet, data_destination) BIND(C, NAME="copyArrayPacketData")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet
                TYPE(C_PTR), VALUE :: data_destination
          END FUNCTION
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION copyArrayPacketDataToInt(packet, data_destination) BIND(C, NAME="copyArrayPacketDataToInt")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet
                TYPE(C_PTR), VALUE :: data_destination
          END FUNCTION
    END INTERFACE

    INTERFACE
         SUBROUTINE freeArrayPacket(packet) BIND(C, NAME="freeArrayPacket")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: packet
          END SUBROUTINE
    END INTERFACE

    INTERFACE
         SUBROUTINE writeArray(array_name, packet) BIND(C, NAME="writeArray")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: array_name
                TYPE(C_PTR),VALUE :: packet
          END SUBROUTINE
    END INTERFACE

    INTERFACE
          TYPE(C_PTR) FUNCTION readArray(array_name) BIND(C, NAME="readArray")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: array_name
          END FUNCTION
    END INTERFACE

    INTERFACE
         SUBROUTINE writeString(string_to_write) BIND(C, NAME="writeString")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR),VALUE :: string_to_write
          END SUBROUTINE
    END INTERFACE

    INTERFACE
          INTEGER(C_INT) FUNCTION readString(max_length_string, read_string) BIND(C, NAME="readString")
                USE, INTRINSIC :: ISO_C_BINDING
                IMPLICIT NONE
                INTEGER(C_INT), VALUE :: max_length_string
                TYPE(C_PTR), VALUE :: read_string
          END FUNCTION
    END INTERFACE

    INTERFACE
          SUBROUTINE endSession() BIND(C, NAME="endSession")
                USE, INTRINSIC :: ISO_C_BINDING

          END SUBROUTINE
    END INTERFACE

!--------------------------------------------------------------------------------------------------
! Define fortran interface for the usage of ipc_write_array and ipc_read_array
!--------------------------------------------------------------------------------------------------
    ! Declare interface for write array.
    INTERFACE ipc_write_array
        MODULE PROCEDURE write_array_double_d1
        MODULE PROCEDURE write_array_double_d2
        MODULE PROCEDURE write_array_double_d3
        MODULE PROCEDURE write_array_double_d4
        MODULE PROCEDURE write_array_double_d5

        MODULE PROCEDURE write_array_integer_d1
        MODULE PROCEDURE write_array_integer_d2
        MODULE PROCEDURE write_array_integer_d3
        MODULE PROCEDURE write_array_integer_d4
        MODULE PROCEDURE write_array_integer_d5
    END INTERFACE

    ! Declare interface for write array.
    INTERFACE ipc_read_array
        MODULE PROCEDURE read_array_double_d1
        MODULE PROCEDURE read_array_double_d2
        MODULE PROCEDURE read_array_double_d3
        MODULE PROCEDURE read_array_double_d4
        MODULE PROCEDURE read_array_double_d5

        MODULE PROCEDURE read_array_integer_d1
        MODULE PROCEDURE read_array_integer_d2
        MODULE PROCEDURE read_array_integer_d3
        MODULE PROCEDURE read_array_integer_d4
        MODULE PROCEDURE read_array_integer_d5
    END INTERFACE

    CONTAINS

!--------------------------------------------------------------------------------------------------
! Log and start
!--------------------------------------------------------------------------------------------------
    SUBROUTINE open_log_file(filename, file_unit) 
         CHARACTER(LEN=*),INTENT(IN) :: filename
         INTEGER, INTENT(INOUT) :: file_unit
 
         LOGICAL :: is_file_open 

         file_unit = 49

         open(unit=file_unit, file=filename, action='write')

         inquire(unit=file_unit,opened=is_file_open)
         if (.not. is_file_open) then
             write(use_unit,*) "ERROR: Could not open log file for writing."
             file_unit = -1
         end if
    END SUBROUTINE open_log_file

    SUBROUTINE close_log_file(file_unit)
        INTEGER :: file_unit

        close(file_unit)
    END SUBROUTINE close_log_file

    SUBROUTINE ipc_start_aims(is_master_process, use_mpi, mpi_comm_global, filename_length, cptr_filename) BIND(C)
         USE, INTRINSIC :: ISO_C_BINDING
         INTEGER(C_INT), VALUE :: is_master_process
         INTEGER(C_INT), VALUE :: use_mpi
         INTEGER(C_INT), VALUE :: mpi_comm_global
         INTEGER(C_INT), VALUE :: filename_length
         TYPE(C_PTR), INTENT(IN), VALUE :: cptr_filename

         CHARACTER(KIND=C_CHAR, LEN=filename_length), POINTER :: filename
         CHARACTER(KIND=C_CHAR, LEN=256) :: f_filename
         INTEGER :: file_unit

         !CALL C_F_POINTER(cptr_filename, filename, [filename_length])
         CALL C_F_POINTER(cptr_filename, filename)

         f_filename(1:filename_length+1) = trim(filename(1:filename_length+1))

         if(is_master_process == 1) then
             CALL open_log_file(f_filename, file_unit)
         else
             ! Slaves write to stdout.
             file_unit = 6 
         endif

         CALL aims( mpi_comm_global, file_unit , use_mpi )

         if(file_unit /= 6) then
             CALL close_log_file(file_unit)
         endif
    END SUBROUTINE ipc_start_aims

!--------------------------------------------------------------------------------------------------
! Fortran <---> C conversion
!--------------------------------------------------------------------------------------------------
    SUBROUTINE clone_and_add_zero_termination(input, clone)
        CHARACTER(LEN=*) :: input, clone

        INTEGER :: input_length

        input_length = LEN_TRIM(input) + 1

        clone(:) = input(:)
        clone(input_length:input_length) = CHAR(0)
    END SUBROUTINE

!--------------------------------------------------------------------------------------------------
! write_array implementations
!--------------------------------------------------------------------------------------------------
    SUBROUTINE write_array_packet(array_name, packet)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        TYPE(C_PTR) :: packet

        CHARACTER(LEN=256), target :: c_array_name

        CALL clone_and_add_zero_termination(array_name, c_array_name)
        CALL writeArray(C_LOC(c_array_name),packet)
    END SUBROUTINE

    SUBROUTINE write_array_packet_double(array_name,i_dimension, i_dimension_length, array_double)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, INTENT(IN), target :: i_dimension, i_dimension_length(i_dimension)
        TYPE(C_PTR) :: array_double
        TYPE(C_PTR) :: packet

        INTEGER ::  array_size
        CHARACTER(LEN=256) :: c_array_name

        array_size = product(i_dimension_length(:))
 
        packet = createArrayPacketDouble(i_dimension, &
                                         C_LOC(i_dimension_length), &
                                         array_double)

        CALL write_array_packet(array_name, packet)
        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE write_array_packet_integer(array_name,i_dimension, i_dimension_length, array_integer)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, INTENT(IN), target :: i_dimension, i_dimension_length(i_dimension)
        TYPE(C_PTR) :: array_integer
        TYPE(C_PTR) :: packet

        INTEGER ::  array_size
        CHARACTER(LEN=256) :: c_array_name

        array_size = product(i_dimension_length(:))

        packet = createArrayPacketInteger(i_dimension, &
                                         C_LOC(i_dimension_length), &
                                         array_integer)

        CALL write_array_packet(array_name, packet)
        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE write_array_double_d1(array_name, array_double_d1)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, DIMENSION(:), target :: array_double_d1

        INTEGER :: i_dimension, i_dimension_length(1)
        
        i_dimension = 1
        i_dimension_length(1) = size(array_double_d1,1)

        CALL write_array_packet_double(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_double_d1))
    END SUBROUTINE

    SUBROUTINE write_array_double_d2(array_name, array_double_d2)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, DIMENSION(:,:), target :: array_double_d2

        INTEGER :: i_dimension, i_dimension_length(2)

        i_dimension = 2
        i_dimension_length(1) = size(array_double_d2,1)
        i_dimension_length(2) = size(array_double_d2,2)

        CALL write_array_packet_double(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_double_d2))
        
    END SUBROUTINE

    SUBROUTINE write_array_double_d3(array_name, array_double_d3)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, DIMENSION(:,:,:), target :: array_double_d3

        INTEGER :: i_dimension, i_dimension_length(3)

        i_dimension = 3
        i_dimension_length(1) = size(array_double_d3,1)
        i_dimension_length(2) = size(array_double_d3,2)
        i_dimension_length(3) = size(array_double_d3,3)

        CALL write_array_packet_double(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_double_d3))
    END SUBROUTINE

    SUBROUTINE write_array_double_d4(array_name, array_double_d4)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, DIMENSION(:,:,:,:), target :: array_double_d4

        INTEGER :: i_dimension, i_dimension_length(4)
        i_dimension = 4

        i_dimension_length(1) = size(array_double_d4,1)
        i_dimension_length(2) = size(array_double_d4,2)
        i_dimension_length(3) = size(array_double_d4,3)
        i_dimension_length(4) = size(array_double_d4,4)

        CALL write_array_packet_double(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_double_d4))
    END SUBROUTINE

    SUBROUTINE write_array_double_d5(array_name, array_double_d5)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, DIMENSION(:,:,:,:,:), target :: array_double_d5

        INTEGER :: i_dimension, i_dimension_length(5)
        i_dimension = 5

        i_dimension_length(1) = size(array_double_d5,1)
        i_dimension_length(2) = size(array_double_d5,2)
        i_dimension_length(3) = size(array_double_d5,3)
        i_dimension_length(4) = size(array_double_d5,4)
        i_dimension_length(5) = size(array_double_d5,5)

        CALL write_array_packet_double(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_double_d5))
    END SUBROUTINE

    SUBROUTINE write_array_integer_d1(array_name, array_integer_d1)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, DIMENSION(:), target :: array_integer_d1

        INTEGER :: i_dimension, i_dimension_length(1)
        
        i_dimension = 1
        i_dimension_length(1) = size(array_integer_d1,1)

        CALL write_array_packet_integer(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_integer_d1))
    END SUBROUTINE

    SUBROUTINE write_array_integer_d2(array_name, array_integer_d2)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, DIMENSION(:,:), target :: array_integer_d2

        INTEGER :: i_dimension, i_dimension_length(2)

        i_dimension = 2
        i_dimension_length(1) = size(array_integer_d2,1)
        i_dimension_length(2) = size(array_integer_d2,2)

        CALL write_array_packet_integer(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_integer_d2))
    END SUBROUTINE

    SUBROUTINE write_array_integer_d3(array_name, array_integer_d3)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, DIMENSION(:,:,:), target :: array_integer_d3

        INTEGER :: i_dimension, i_dimension_length(3)

        i_dimension = 3
        i_dimension_length(1) = size(array_integer_d3,1)
        i_dimension_length(2) = size(array_integer_d3,2)
        i_dimension_length(3) = size(array_integer_d3,3)

        CALL write_array_packet_integer(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_integer_d3))
    END SUBROUTINE

    SUBROUTINE write_array_integer_d4(array_name, array_integer_d4)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, DIMENSION(:,:,:,:), target :: array_integer_d4

        INTEGER :: i_dimension, i_dimension_length(4)
        i_dimension = 4

        i_dimension_length(1) = size(array_integer_d4,1)
        i_dimension_length(2) = size(array_integer_d4,2)
        i_dimension_length(3) = size(array_integer_d4,3)
        i_dimension_length(4) = size(array_integer_d4,4)

        CALL write_array_packet_integer(array_name, i_dimension, i_dimension_length, &
                                       C_LOC(array_integer_d4))
    END SUBROUTINE

    SUBROUTINE write_array_integer_d5(array_name, array_integer_d5)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, DIMENSION(:,:,:,:,:), target :: array_integer_d5

        INTEGER :: i_dimension, i_dimension_length(5)
        i_dimension = 5

        i_dimension_length(1) = size(array_integer_d5,1)
        i_dimension_length(2) = size(array_integer_d5,2)
        i_dimension_length(3) = size(array_integer_d5,3)
        i_dimension_length(4) = size(array_integer_d5,4)
        i_dimension_length(5) = size(array_integer_d5,5)

        CALL write_array_packet_integer(array_name, i_dimension, i_dimension_length, &
                                        C_LOC(array_integer_d5))
    END SUBROUTINE

!--------------------------------------------------------------------------------------------------
! read_array implementations
!--------------------------------------------------------------------------------------------------
    SUBROUTINE read_array_packet(array_name, packet)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        TYPE(C_PTR) :: packet

        CHARACTER(LEN=256), target :: c_array_name

        CALL clone_and_add_zero_termination(array_name, c_array_name)
        packet = readArray(C_LOC(c_array_name))
    END SUBROUTINE

    SUBROUTINE read_array_double_d1(array_name, array_double_d1)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), target :: array_double_d1

        TYPE(C_PTR) :: packet
        INTEGER :: i_read

        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_double_d1(arrayPacketDimensionSize(packet, 1)))

        i_read = copyArrayPacketData(packet, C_LOC(array_double_d1))

        CALL freeArrayPacket(packet)

    END SUBROUTINE

    SUBROUTINE read_array_double_d2(array_name, array_double_d2)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:), target :: array_double_d2

        TYPE(C_PTR) :: packet
        INTEGER :: i_read

        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_double_d2(arrayPacketDimensionSize(packet, 1),&
                                 arrayPacketDimensionSize(packet, 2)))

        i_read = copyArrayPacketData(packet, C_LOC(array_double_d2))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_double_d3(array_name, array_double_d3)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:), target :: array_double_d3

        TYPE(C_PTR) :: packet
        INTEGER :: i_read
        
        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_double_d3(arrayPacketDimensionSize(packet, 1),&
                                 arrayPacketDimensionSize(packet, 2),&
                                 arrayPacketDimensionSize(packet, 3)))

        i_read = copyArrayPacketData(packet, C_LOC(array_double_d3))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_double_d4(array_name, array_double_d4)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:,:), target :: array_double_d4
        TYPE(C_PTR) :: packet
        INTEGER :: i_read
        
        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_double_d4(arrayPacketDimensionSize(packet, 1),&
                                 arrayPacketDimensionSize(packet, 2),&
                                 arrayPacketDimensionSize(packet, 3),&
                                 arrayPacketDimensionSize(packet, 4)))

        i_read = copyArrayPacketData(packet, C_LOC(array_double_d4))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_double_d5(array_name, array_double_d5)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:,:,:), target :: array_double_d5
        TYPE(C_PTR) :: packet
        INTEGER :: i_read
        
        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_double_d5(arrayPacketDimensionSize(packet, 1),&
                                 arrayPacketDimensionSize(packet, 2),&
                                 arrayPacketDimensionSize(packet, 3),&
                                 arrayPacketDimensionSize(packet, 4),&
                                 arrayPacketDimensionSize(packet, 5)))

        i_read = copyArrayPacketData(packet, C_LOC(array_double_d5))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_integer_d1(array_name, array_integer_d1)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, ALLOCATABLE, DIMENSION(:), target :: array_integer_d1

        TYPE(C_PTR) :: packet
        INTEGER :: i_read

        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_integer_d1(arrayPacketDimensionSize(packet, 1)))

        i_read = copyArrayPacketDataToInt(packet, C_LOC(array_integer_d1))

        CALL freeArrayPacket(packet)

    END SUBROUTINE

    SUBROUTINE read_array_integer_d2(array_name, array_integer_d2)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, ALLOCATABLE,DIMENSION(:,:), target :: array_integer_d2

        TYPE(C_PTR) :: packet
        INTEGER :: i_read

        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_integer_d2(arrayPacketDimensionSize(packet, 1),&
                                  arrayPacketDimensionSize(packet, 2)))

        i_read = copyArrayPacketDataToInt(packet, C_LOC(array_integer_d2))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_integer_d3(array_name, array_integer_d3)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, ALLOCATABLE,DIMENSION(:,:,:), target :: array_integer_d3

        TYPE(C_PTR) :: packet
        INTEGER :: i_read
        
        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_integer_d3(arrayPacketDimensionSize(packet, 1),&
                                  arrayPacketDimensionSize(packet, 2),&
                                  arrayPacketDimensionSize(packet, 3)))

        i_read = copyArrayPacketDataToInt(packet, C_LOC(array_integer_d3))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_integer_d4(array_name, array_integer_d4)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:), target :: array_integer_d4
        TYPE(C_PTR) :: packet
        INTEGER :: i_read
        
        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_integer_d4(arrayPacketDimensionSize(packet, 1),&
                                  arrayPacketDimensionSize(packet, 2),&
                                  arrayPacketDimensionSize(packet, 3),&
                                  arrayPacketDimensionSize(packet, 4)))

        i_read = copyArrayPacketDataToInt(packet, C_LOC(array_integer_d4))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

    SUBROUTINE read_array_integer_d5(array_name, array_integer_d5)
        CHARACTER(LEN=*), INTENT(IN) :: array_name
        INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:,:), target :: array_integer_d5
        TYPE(C_PTR) :: packet
        INTEGER :: i_read
        
        CALL read_array_packet(array_name, packet)

        ALLOCATE(array_integer_d5(arrayPacketDimensionSize(packet, 1),&
                                  arrayPacketDimensionSize(packet, 2),&
                                  arrayPacketDimensionSize(packet, 3),&
                                  arrayPacketDimensionSize(packet, 4),&
                                  arrayPacketDimensionSize(packet, 5)))

        i_read = copyArrayPacketDataToInt(packet, C_LOC(array_integer_d5))

        CALL freeArrayPacket(packet)
    END SUBROUTINE

!--------------------------------------------------------------------------------------------------
! Write and read string
!--------------------------------------------------------------------------------------------------
    SUBROUTINE ipc_write_string(string_to_write)
        CHARACTER(LEN=*) :: string_to_write
        CHARACTER(LEN=1024), target :: c_string_to_write

        CALL clone_and_add_zero_termination(string_to_write, c_string_to_write)

        CALL writeString(C_LOC(c_string_to_write))
    END SUBROUTINE

    SUBROUTINE ipc_read_string(max_size, string)
        INTEGER :: max_size
        CHARACTER(LEN=max_size), target :: string

        INTEGER :: characters_read, i_run

        characters_read = readString(max_size, C_LOC(string))

        do i_run=characters_read+1, max_size
            string(i_run+1:i_run+1) = ''
        enddo
    END SUBROUTINE

!--------------------------------------------------------------------------------------------------
! Transaction and session
!--------------------------------------------------------------------------------------------------
    FUNCTION ipc_start_transaction(transaction_name) RESULT(is_transaction_active)
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(len=*) :: transaction_name
        LOGICAL :: is_transaction_active
        TYPE(C_PTR) :: stringPtr
        INTEGER :: reply
        CHARACTER(LEN=256), target :: c_transaction_name

        INTEGER :: TRANSACTION_UNHANDELED = 2001
        INTEGER :: TRANSACTION_REJECTED = 2002
        INTEGER :: TRANSACTION_ACTIVE = 2003

        CALL clone_and_add_zero_termination(transaction_name, c_transaction_name)

        stringPtr = C_LOC(c_transaction_name)
        reply = startTransaction(stringPtr)

        is_transaction_active = .FALSE.

        if(reply == 1) then
            is_transaction_active = .TRUE.
        endif
    END FUNCTION

    SUBROUTINE ipc_end_session()
        CALL endSession()
    END SUBROUTINE

!--------------------------------------------------------------------------------------------------
! Unittests (invoked from outside peer)
!--------------------------------------------------------------------------------------------------
    SUBROUTINE mirrorTest() BIND(C, NAME="mirrorTest")
        TYPE(C_PTR) :: packet
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: test_array
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_array_integer

        INTEGER :: i_index
        CHARACTER(LEN=1024), target :: test_array_name

        write(use_unit,*) "FHI-aims: MIRROR test starts."

        if(ipc_start_transaction("MIRROR_TEST")) then
            CALL clone_and_add_zero_termination("TEST_ARRAY", test_array_name)

            packet = readArray(C_LOC(test_array_name))

            write(use_unit,*) "FHI-aims: Got packet. Going to send it back."
            CALL writeArray(C_LOC(test_array_name), packet)

            write(use_unit,*) "FHI-aims: Sent packet back."
            CALL freeArrayPacket(packet)
        endif

        CALL ipc_end_session()
    END SUBROUTINE mirrorTest

    SUBROUTINE testWrite_array() BIND(C, NAME="testWrite_array")
        TYPE(C_PTR) :: packet
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: test_array
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_array_integer
        INTEGER :: i,j

        write(use_unit,*) "FHI-aims: write_array test starts."

        if(ipc_start_transaction("WRITE_ARRAY_TEST")) then
            allocate(test_array(10,30))
            do j=1, 30
               do i=1, 10
                    test_array(i,j) = dble(j) * 10001.d-2 + dble(i)
               enddo
            enddo

            call ipc_write_array("test_array_dbl", test_array)

            allocate(test_array_integer(10,30))
            do j=1, 30
               do i=1, 10
                    test_array_integer(i,j) = j * 100 + i
               enddo
            enddo
            call ipc_write_array("test_array_int", test_array_integer)

        endif

        CALL ipc_end_session()
    END SUBROUTINE testWrite_array

    SUBROUTINE testRead_array() BIND(C, NAME="testRead_array")
        TYPE(C_PTR) :: packet
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: test_array
        DOUBLE PRECISION :: diff
        INTEGER :: diff_int
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: test_array_integer
        INTEGER :: i,j
        write(use_unit,*) "FHI-aims: read_array test starts."

        if(ipc_start_transaction("READ_ARRAY_TEST")) then
            call ipc_read_array("test_array_dbl", test_array)
            call ipc_read_array("test_array_int", test_array_integer)

            do j=1, 30
               do i=1, 10
                    diff = test_array(i,j) - (dble(j) * 10001.d-2 + dble(i))
                    if(abs(diff) > 1.d-6) then
                        write(use_unit,*) "Found differnece:" ,i,j, diff, test_array(i,j)
                        stop "ABORT"
                    endif

                    diff_int = test_array_integer(i,j) - (j * 100 + i)
                    if(abs(diff_int) > 0) then
                        write(use_unit,*) "Found differnece:" ,i,j, diff_int, test_array_integer(i,j)
                        stop "ABORT"
                    endif
               enddo
            enddo

            deallocate(test_array)
            deallocate(test_array_integer)
        endif

        CALL ipc_end_session()
    END SUBROUTINE testRead_array

    SUBROUTINE testWrite_string() BIND(C, NAME="testWrite_string")
        write(use_unit,*) "FHI-aims: write_string test starts."

        if(ipc_start_transaction("WRITE_STRING_TEST")) then
           CALL ipc_write_string("write test string")   
        endif
    END SUBROUTINE testWrite_string

    SUBROUTINE testRead_string() BIND(C, NAME="testRead_string")
        CHARACTER(len=80) :: string

        write(use_unit,*) "FHI-aims: read_string test starts."

        if(ipc_start_transaction("READ_STRING_TEST")) then
            call ipc_read_string(80, string)

            if(.NOT.(string(1:16).EQ."read test string")) then
               write(use_unit,*) "Unexpectedly read string", string
               stop "ABORT"
            endif
        endif
    END SUBROUTINE testRead_string
END MODULE ipc
