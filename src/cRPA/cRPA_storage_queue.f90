!Implments a simple queue
MODULE cRPA_storage_queue
    IMPLICIT NONE

    type queue
        INTEGER :: queue_length
        INTEGER :: queue_extend
        LOGICAL :: is_new_queue
        INTEGER :: queue_pos
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: data
    end type queue

CONTAINS
    SUBROUTINE create_queue(length,extend,aqueue)
        INTEGER, INTENT(IN) :: length,extend
        type(queue), INTENT(OUT) :: aqueue

        if(length == 1) then
            stop "queue of length 1 forbidden"
        endif

        ALLOCATE(aqueue%data(2,length,extend))

!TODO: check alloc

        aqueue%queue_length = length
        aqueue%queue_extend = extend

        aqueue%queue_pos = 1
        aqueue%is_new_queue = .TRUE.

    END SUBROUTINE create_queue

    SUBROUTINE enqueue(aqueue,n_datalen, data)
         type(queue), INTENT(INOUT) :: aqueue
         INTEGER, INTENT(IN) :: n_datalen
         INTEGER, DIMENSION(n_datalen), INTENT(IN) :: data

         if(.NOT. can_enqueue(aqueue)) then
            stop "error in enqueue: queue is too small. can not override data!"
         endif

         aqueue%data(1,aqueue%queue_pos,1:aqueue%queue_extend) = data(:)

         aqueue%queue_pos = aqueue%queue_pos + 1

         if(aqueue%is_new_queue) then
            aqueue%is_new_queue = .FALSE.
         endif

    END SUBROUTINE enqueue

    LOGICAL FUNCTION can_enqueue(aqueue)
         type(queue), INTENT(IN) :: aqueue

         if(aqueue%queue_pos > aqueue%queue_length) then
            can_enqueue = .FALSE.
            return
         endif

         can_enqueue = .TRUE.
    END FUNCTION can_enqueue

    SUBROUTINE dequeue(aqueue,data)
         type(queue), INTENT(INOUT) :: aqueue
         INTEGER, DIMENSION(aqueue%queue_extend) :: data

         if(.NOT. can_dequeue(aqueue)) then
            stop "error in dequeue: there is no data to dequeue!"
         endif


         data(1:aqueue%queue_extend) = aqueue%data(1,1,:)

         aqueue%data(2,1:aqueue%queue_pos - 2,:) = aqueue%data(1,2:aqueue%queue_pos-1,:)
         aqueue%queue_pos = aqueue%queue_pos - 1
         aqueue%data(1,1:aqueue%queue_pos-1,:) = aqueue%data(2,1:aqueue%queue_pos-1,:)

    END SUBROUTINE dequeue


    LOGICAL FUNCTION can_dequeue(aqueue)
         type(queue), INTENT(IN) :: aqueue

         if(aqueue%queue_pos == 1) then
            can_dequeue = .FALSE.
            return
         endif

         can_dequeue = .TRUE.
    END FUNCTION can_dequeue

    SUBROUTINE shrink_queue(aqueue)
        type(queue), INTENT(INOUT) :: aqueue

        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: tmp

        ALLOCATE(tmp(2,aqueue%queue_length,aqueue%queue_extend))
!TODO check alloc
        tmp(:,:,:) = aqueue%data(:,:,:)

        aqueue%queue_length = aqueue%queue_pos - 1
        DEALLOCATE(aqueue%data)
        ALLOCATE(aqueue%data(2,aqueue%queue_length,aqueue%queue_extend))
!TODO check alloc

        aqueue%data(:,:,:) = tmp(:,1:aqueue%queue_length,1:aqueue%queue_extend)

        DEALLOCATE(tmp)
    END SUBROUTINE


    SUBROUTINE free_queue(aqueue)
        type(queue), INTENT(INOUT) :: aqueue

        DEALLOCATE(aqueue%data)

        aqueue%queue_length = 0
        aqueue%queue_extend = 0
    END SUBROUTINE free_queue
END MODULE cRPA_storage_queue
