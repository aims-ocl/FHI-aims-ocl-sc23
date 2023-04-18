
MODULE cRPA_parallelism_planing
    USE cRPA_view
    USE cRPA_flow_realspace_vector
    USE cRPA_parallelism_storage
    USE cRPA_storage
    IMPLICIT NONE

    type exp_coeff_calculation_distribution
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: i_cell_by_id, &
                                                i_atoms_to_id
        INTEGER, ALLOCATABLE, DIMENSION(:) :: i_cell_to_id, &
                                              i_atom_pair_my_to_atom1, &
                                              i_atom_pair_my_to_atom2
        INTEGER :: n_atoms, n_atom_pairs_my, &
                   n_cells, n_cells_my
        INTEGER :: n_cell_tasks_max
    end type

CONTAINS
    SUBROUTINE create_exp_coeff_calc_dist(n_atoms,n_cells,a_comm, a_dist)
        INTEGER, INTENT(IN) :: n_atoms,n_cells
        type(custom_communicator) :: a_comm
        type(exp_coeff_calculation_distribution) :: a_dist

        INTEGER :: i_atom_1,i_atom_2,i_cell, i_pos, i_id

        a_dist%n_atoms = n_atoms
        a_dist%n_cells = n_cells
        a_dist%n_cell_tasks_max = 1 + n_cells / a_comm%n_tasks

        ALLOCATE(a_dist%i_atoms_to_id(n_atoms,n_atoms))
        ALLOCATE(a_dist%i_atom_pair_my_to_atom1(n_atoms**2))
        ALLOCATE(a_dist%i_atom_pair_my_to_atom2(n_atoms**2))

        a_dist%n_atom_pairs_my = 0
        do i_atom_1 = 1,a_dist%n_atoms
            do i_atom_2 = 1,a_dist%n_atoms

                i_id = mod( (i_atom_1-1) * n_atoms + i_atom_2,a_comm%n_tasks)

                a_dist%i_atoms_to_id(i_atom_1,i_atom_2) = i_id

                if(i_id ==a_comm%myid) then
                    a_dist%n_atom_pairs_my = a_dist%n_atom_pairs_my +1
                    a_dist%i_atom_pair_my_to_atom1(a_dist%n_atom_pairs_my) = i_atom_1
                    a_dist%i_atom_pair_my_to_atom2(a_dist%n_atom_pairs_my) = i_atom_2
                endif
            enddo
        enddo


        ALLOCATE(a_dist%i_cell_by_id(0:a_comm%n_tasks-1, &
                                      a_dist%n_cell_tasks_max))
        ALLOCATE(a_dist%i_cell_to_id(n_cells))

        a_dist%i_cell_by_id(:,:) = -1
        i_pos = 1
        do i_cell = 1,a_dist%n_cells
           i_id = mod(i_cell,a_comm%n_tasks)

           a_dist%i_cell_to_id(i_cell) = i_id
           a_dist%i_cell_by_id(i_id,i_pos) = i_cell

           if(i_id == a_comm%n_tasks-1) i_pos = i_pos +1
        enddo

        a_dist%n_cells_my = COUNT(a_dist%i_cell_by_id(a_comm%myid,:)/=-1)
    END SUBROUTINE create_exp_coeff_calc_dist

    SUBROUTINE free_exp_coeff_calc_dist(a_dist)
        type(exp_coeff_calculation_distribution) :: a_dist

        DEALLOCATE(a_dist%i_cell_by_id)
        DEALLOCATE(a_dist%i_cell_to_id)

        DEALLOCATE(a_dist%i_atom_pair_my_to_atom1)
        DEALLOCATE(a_dist%i_atom_pair_my_to_atom2)
    END SUBROUTINE free_exp_coeff_calc_dist
!exp coeffs
    SUBROUTINE distribute_exp_coeff_integration(n_species, n_vectors,n_max_vectors_per_call, &
                                                num_cpu, worklist)
        INTEGER, INTENT(IN) :: n_species,n_vectors,n_max_vectors_per_call,num_cpu
        type(exp_coeff_int_distribution), DIMENSION(num_cpu), INTENT(OUT) :: worklist

        type(exp_coeff_int_task), ALLOCATABLE, DIMENSION(:) :: worktasks

        INTEGER :: i,j,k,run_cpu, num_elements, n_vectors_per_call, &
                   offset_k, offset_task, offset_i_j, diff, &
                   num_combinations, minimal_chunksize,rest_of_partion

        if(n_species == 0 &
           .OR. n_vectors == 0 &
           .OR. n_max_vectors_per_call == 0 &
           .OR. num_cpu == 0) then

           CALL write_debug("species,vectors or cpu is zero - this should not happen")
           stop
        endif

        num_combinations = ( n_species*(n_species+1) ) / 2
        num_combinations = num_combinations * n_vectors

        CALL factorize_integer_product(num_cpu,num_combinations, &
                                       minimal_chunksize,rest_of_partion)

        ALLOCATE(worktasks(num_combinations))
        if(.NOT. ALLOCATED(worktasks)) then
            write(use_unit,*) "error in dist int"
            stop
        endif


        i=1
        j=1
        offset_k = 1
        offset_i_j = n_vectors
        do run_cpu = 1, num_cpu

            CALL get_num_elements(num_elements)


            offset_task=0

            do while (num_elements>0)
                n_vectors_per_call = min(n_max_vectors_per_call, num_elements)
                offset_task = offset_task + 1
!call write_debug("cpu"//int2str(run_cpu)//"num_ele:"//int2str(num_elements))
                if(offset_i_j >= n_vectors_per_call) then
                   CALL set_task(min(n_vectors_per_call,n_vectors), &
                                 worktasks(offset_task))
                   CALL increase_i_j_if_possible()
                else
                   CALL set_task(min(offset_i_j,n_vectors), &
                                 worktasks(offset_task))
                   CALL increase_i_j_if_possible()
                end if
            end do

            CALL copy_worktasks_to_worklist(offset_task)
        end do


        CONTAINS
            SUBROUTINE increase_i_j_if_possible()
                if(offset_i_j /= 0) return

                if(j<n_species) then
                    j=j+1
                elseif(j==n_species) then
                    i=i+1
                    j=i
                end if

                offset_k = 1
                offset_i_j = n_vectors
            END SUBROUTINE increase_i_j_if_possible

            SUBROUTINE set_task(amount_vectors,task)
                INTEGER, INTENT(IN) :: amount_vectors
                type(exp_coeff_int_task), INTENT(INOUT) :: task

!                write(use_unit,*) "setting",i,j,offset_k,amount_vectors+offset_k-1
                task = exp_coeff_int_task(i,j,offset_k,amount_vectors+offset_k-1)
                offset_i_j = offset_i_j - amount_vectors
                offset_k = offset_k + amount_vectors
                num_elements = num_elements-amount_vectors
            END SUBROUTINE set_task

            SUBROUTINE get_num_elements(num_elements)
                INTEGER, INTENT(INOUT) :: num_elements

                num_elements = minimal_chunksize

                if(run_cpu <= rest_of_partion) then
                    num_elements = num_elements+1
                endif
!CALL write_debug(int2str(num_elements))
            END SUBROUTINE get_num_elements

            SUBROUTINE copy_worktasks_to_worklist(worktasks_length)
                INTEGER, INTENT(IN) :: worktasks_length
                INTEGER :: r

                ALLOCATE(worklist(run_cpu)%tasks(worktasks_length))
!TODO check allocation

                do r = 1, worktasks_length
                    worklist(run_cpu)%tasks(r) = exp_coeff_int_task( &
                            worktasks(r)%species1,worktasks(r)%species2, &
                            worktasks(r)%min_vector,worktasks(r)%max_vector)
                end do
            END SUBROUTINE copy_worktasks_to_worklist

    END SUBROUTINE distribute_exp_coeff_integration

    SUBROUTINE create_integration_distribution(num_cpu,worklist)
        INTEGER, INTENT(IN) :: num_cpu
        type(exp_coeff_int_distribution), ALLOCATABLE, DIMENSION(:) &
                                        , INTENT(OUT) :: worklist

!TODO check allocation
        ALLOCATE(worklist(num_cpu))

    END SUBROUTINE create_integration_distribution

    SUBROUTINE get_max_tasks_of_integration_distribution(int_dist,max_tasks)
        type(exp_coeff_int_distribution),INTENT(IN) :: int_dist
        INTEGER, INTENT(OUT) :: max_tasks

!        INTEGER, DIMENSION(ubound(int_dist,1)) :: num_tasks
!        INTEGER :: i

!        do i = myid+1,myid+1 !1, ubound(int_dist,1)
!            num_tasks(i) =ubound(int_dist(i)%tasks,1)
!        end do
        max_tasks = ubound(int_dist%tasks,1) ! maxval(num_tasks)

    END SUBROUTINE get_max_tasks_of_integration_distribution

    SUBROUTINE free_integration_distribution(num_cpu,worklist)
        INTEGER, INTENT(IN) :: num_cpu
        type(exp_coeff_int_distribution), DIMENSION(num_cpu) &
                                        ,INTENT(INOUT) :: worklist
        INTEGER :: i

        do i = 1, ubound(worklist,1)
            DEALLOCATE(worklist(i)%tasks)
        end do

    END SUBROUTINE free_integration_distribution

!basbas
    SUBROUTINE get_basbas_row_block_number_by_myid(n_tasks_global,myid_global, &
                                                   n_blocks, &
                                                   i_block)
        INTEGER,INTENT(IN) :: n_tasks_global,myid_global, n_blocks
        INTEGER,INTENT(OUT) :: i_block

        INTEGER :: i_block_start, n_rank_block

        CALL get_basbas_row_block_shape_by_myid(n_tasks_global,myid_global, &
                                                n_blocks, &
                                                i_block,i_block_start, n_rank_block)

    END SUBROUTINE get_basbas_row_block_number_by_myid

    SUBROUTINE get_basbas_row_block_extend_by_block(i_block,n_tasks_global,myid_global, &
                                                    n_blocks, n_rank_block)
        INTEGER, INTENT(IN) :: i_block,n_tasks_global,myid_global, n_blocks
        INTEGER, INTENT(OUT) :: n_rank_block

        INTEGER :: rest_block, n_block_mean

        CALL factorize_integer_product(n_blocks,n_tasks_global,n_block_mean,rest_block)

        if(i_block-1 < rest_block) then
          n_rank_block = n_block_mean + 1
        else
          n_rank_block = n_block_mean
        endif


    END SUBROUTINE get_basbas_row_block_extend_by_block

    SUBROUTINE get_basbas_row_block_shape_by_myid(n_tasks_global,myid_global, &
                                                  n_blocks, &
                                                  i_block,i_block_start, n_rank_block)


        INTEGER,INTENT(IN) :: n_tasks_global,myid_global, n_blocks
        INTEGER,INTENT(OUT) :: i_block, i_block_start, n_rank_block

        n_rank_block = 0
        i_block_start = 0
        do i_block = 1,n_blocks
           i_block_start = i_block_start + n_rank_block

           CALL get_basbas_row_block_extend_by_block(i_block,n_tasks_global, &
                                                     myid_global, n_blocks, &
                                                     n_rank_block)

           if(i_block_start+n_rank_block -1 >= myid_global) then
             exit
           endif
        enddo

    END SUBROUTINE get_basbas_row_block_shape_by_myid

    SUBROUTINE get_cpus_basbas_row_block(n_tasks_global,myid_global, &
                                         n_blocks, &
                                         n_rank_block,rank_block)
        INTEGER,INTENT(IN) :: n_tasks_global,myid_global, n_blocks

        INTEGER,INTENT(OUT) :: n_rank_block
        INTEGER,ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: rank_block

        INTEGER :: i_block, n_block_mean, rest_block, rest_id, n_block_offset, &
                   i, i_block_start


        CALL get_basbas_row_block_shape_by_myid(n_tasks_global,myid_global, &
                                                 n_blocks, &
                                                 i_block, i_block_start, n_rank_block)
!TODO check alloc
        ALLOCATE(rank_block(n_rank_block))
!write(use_unit,*) myid_global,"block",i_block,"block_start",i_block_start,"block_mean", n_block_mean,"n_rank_block", n_rank_block
        FORALL(i = 0:n_rank_block-1) rank_block(i+1)=i_block_start+i
!write(use_unit,*) rank_block
    END SUBROUTINE get_cpus_basbas_row_block

!basbas row blocks
    SUBROUTINE distribute_basbas_row_blocks(num_aux,num_blocks, basbas_dist)
        INTEGER, INTENT(IN) :: num_aux, num_blocks
        type(basbas_distribution),ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: basbas_dist

!TODO check alloc
        ALLOCATE(basbas_dist(num_blocks))
        CALL create_basbas_distribution_simple(num_aux,0, num_blocks,basbas_dist)

    END SUBROUTINE distribute_basbas_row_blocks



!basbas rows within block
    SUBROUTINE distribute_basbas(num_aux,num_blocks,num_cpu,worklist, worklist_total)
        INTEGER, INTENT(IN) :: num_aux, &
                               num_blocks, &
                               num_cpu

        type(basbas_distribution),ALLOCATABLE, DIMENSION(:), INTENT(OUT) ::  &
                                                                            worklist, &
                                                                            worklist_total
        type(basbas_distribution),ALLOCATABLE, DIMENSION(:) :: basbas_block_dist
        INTEGER :: i,i_block,n_rank_block, i_start, i_end

!        CALL write_stdout("Distributing "//num2str(num_aux)//" auxilliary functions over" &
!                          //num2str(num_cpu)//" cpus")


!TODO check alloc
        ALLOCATE(basbas_block_dist(num_blocks))
        CALL distribute_basbas_row_blocks(num_aux,num_blocks,basbas_block_dist)

!TODO check alloc
        CALL get_basbas_row_block_number_by_myid(num_cpu,myid, &
                                                 num_blocks, &
                                                 i_block)

        ALLOCATE(worklist_total(num_cpu))

        i_start = 1
        do i=1,num_blocks
            CALL get_basbas_row_block_extend_by_block(i,num_cpu, &
                                                      myid, num_blocks, &
                                                      n_rank_block)

            i_end=i_start-1+n_rank_block
            CALL create_basbas_distribution_for_block(n_rank_block,basbas_block_dist(i),&
                                                      worklist_total(i_start:i_end) )
!            CALL print_basbas_distributions(n_rank_block,worklist_total(i_start:i_end))
            i_start=i_end+1


            if(i==i_block) then
 !               CALL write_debug("Taking basbas row block:"//num2str(i_block))
                ALLOCATE(worklist(n_rank_block))
                CALL create_basbas_distribution_for_block(n_rank_block,basbas_block_dist(i_block),worklist)
            endif
        enddo

!        CALL print_basbas_distributions(num_cpu,worklist_total)

    END SUBROUTINE distribute_basbas

    SUBROUTINE create_basbas_distribution_for_block(num_cpu, basbas_block_dist, &
                                                    basbas_dist)
       INTEGER,INTENT(IN) :: num_cpu
       type(basbas_distribution), INTENT(IN) :: basbas_block_dist
       type(basbas_distribution),DIMENSION(num_cpu) , INTENT(OUT) :: basbas_dist

       CALL create_basbas_distribution_simple(basbas_block_dist%n_rows, &
                                              basbas_block_dist%n_offset_rows, &
                                              num_cpu, basbas_dist)
    END SUBROUTINE create_basbas_distribution_for_block

    SUBROUTINE create_basbas_distribution_simple(num_aux,num_aux_offset, num_parts,worklist)
        INTEGER, INTENT(in) :: num_aux,num_aux_offset, num_parts
        type(basbas_distribution),DIMENSION(num_parts) , INTENT(OUT) :: worklist

        INTEGER :: i_part, n_row,rest, n_elements, my_n_row, i,j, row_offset

        CALL factorize_integer_product(num_parts,num_aux,n_row,rest)

        row_offset = num_aux_offset
        do i_part = 1, num_parts
            if( (num_parts-i_part) < rest) then
                my_n_row = n_row + 1
            else
                my_n_row = n_row
            endif

            n_elements = num_aux*my_n_row
            CALL create_basbas_distribution(n_elements,worklist(i_part))

!TODO refactor
            worklist(i_part)%n_rows = my_n_row
            worklist(i_part)%n_columns = n_basbas
            worklist(i_part)%n_offset_rows = row_offset
            worklist(i_part)%n_offset_columns = 1

            do i = 1,my_n_row
                do j=1, num_aux
                    worklist(i_part)%element((i-1)*num_aux+j) = point(i+row_offset,j)
                enddo
            enddo

           row_offset = row_offset + my_n_row
        enddo
    END SUBROUTINE create_basbas_distribution_simple

!basis functions
    SUBROUTINE distribute_basisfunctions(n_cpu, n_basisfuncs, rvecs, distribution)
        INTEGER, INTENT(IN) :: n_cpu, n_basisfuncs
        type(realspace_vectors), INTENT(IN) :: rvecs
        type(basisfunc_distribution), DIMENSION(n_cpu),INTENT(OUT) :: distribution

        INTEGER :: i, mean_size, rest, n_offset, n_end, atom_start, atom_end

        CALL write_stdout("Creating distribution strategy for " // num2str(n_basisfuncs) &
                          // "basisfuncs over " // num2str(n_cpu) // " CPUs")


        CALL factorize_integer_product(n_cpu,n_basisfuncs, mean_size, rest)

        n_offset = 1
        n_end = 0
        do i = 1, n_cpu
            n_end = n_end + mean_size

            if(i> n_cpu - rest) then
                n_end = n_end + 1
            endif

            CALL create_basisfunc_distribution(n_offset,n_end, rvecs, distribution(i))
            n_offset = n_end + 1
        end do

        if(myid == 0) then
            do i=1,n_cpu
                CALL write_stdout("Basis dist for task "//num2str(i))
                CALL print_basisfunc_distribution(distribution(i))
            enddo
        endif
    END SUBROUTINE distribute_basisfunctions

!auxilliary functions

    !f1*f2+r=p
    SUBROUTINE factorize_integer_product(factor1,prod,factor2,rest)
        INTEGER, INTENT(IN) :: factor1,prod
        INTEGER, INTENT(OUT) :: factor2,rest

        factor2 = INT(FLOOR( FLOAT(prod) / FLOAT(factor1) ))
        rest = prod - factor1 * factor2

    END SUBROUTINE factorize_integer_product


END MODULE cRPA_parallelism_planing
