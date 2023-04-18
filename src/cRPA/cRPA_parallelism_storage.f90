!Manages storages of parallelism related tasks (planning, buffering of communication)
MODULE crpa_parallelism_storage
    USE physics
    USE basis
    USE geometry
    USE pbc_lists
    USE prodbas
    USE dimensions
    USE cRPA_view
    USE cRPA_flow_realspace_vector

    implicit none
    INTEGER ::COMM_MATRIX_HEADER
    PARAMETER( COMM_MATRIX_HEADER = 16)

    type message_matrix_send
        INTEGER :: source
        INTEGER :: message_id
        INTEGER :: sources_index
        INTEGER :: offset
        INTEGER :: tag
        INTEGER :: i_cell1
        INTEGER :: i_cell2
        INTEGER :: i_cell3
        INTEGER :: dim1
        INTEGER :: dim2
        INTEGER :: dim3
        REAL*8, ALLOCATABLE,DIMENSION(:) :: data
    end type

    type message_matrix_recv
        INTEGER :: source
        INTEGER :: message_id
        INTEGER :: sources_index
        INTEGER :: offset
        INTEGER :: tag
        INTEGER :: i_cell1
        INTEGER :: i_cell2
        INTEGER :: i_cell3
        INTEGER :: dim1
        INTEGER :: dim2
        INTEGER :: dim3
        REAL*8, ALLOCATABLE,DIMENSION(:,:,:) :: data
    end type

    type exp_coeff_int_task
        INTEGER :: species1, species2, &
                   min_vector, max_vector
    end type exp_coeff_int_task

    type exp_coeff_int_distribution
        type(exp_coeff_int_task), ALLOCATABLE, DIMENSION(:) :: tasks
    end type exp_coeff_int_distribution

    type realspace_vector_distribution
        type(realspace_vectors) :: my_vecs
        INTEGER :: n_my_vecs
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: vec_coeffs_to_cpu
    end type realspace_vector_distribution

    type basisfunc_distribution_atom
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: overlaping_atoms_connection_vector
        INTEGER, ALLOCATABLE, DIMENSION(:) :: overlaping_atoms
        INTEGER, ALLOCATABLE, DIMENSION(:) :: overlaping_atoms_supercell
        INTEGER :: part_of_atom_start
        INTEGER :: part_of_atom_end
        INTEGER :: my_basis_atom_start
        INTEGER :: my_basis_atom_end
        INTEGER :: n_overlaping_atoms
    end type basisfunc_distribution_atom

    type basisfunc_distribution
        type(basisfunc_distribution_atom), ALLOCATABLE, DIMENSION(:) :: atoms
        type(realspace_vectors) :: rvecs
        INTEGER :: basisfunc_start
        INTEGER :: basisfunc_end
        INTEGER :: basisfunc_extend
        INTEGER :: atom_start
        INTEGER :: atom_end
        INTEGER :: my_basis_start
        INTEGER :: my_basis_end
    end type basisfunc_distribution

    type point
        INTEGER :: x,y
    end type point

    type matrixblock_of_point
        type(point), ALLOCATABLE, DIMENSION(:,:) :: element
        INTEGER :: dim_row, dim_column
    end type matrixblock_of_point

    type basbas_distribution
        type(point), ALLOCATABLE, DIMENSION(:) :: element
        INTEGER :: dim_row
        INTEGER :: n_rows
        INTEGER :: n_columns
        INTEGER :: n_offset_rows
        INTEGER :: n_offset_columns
    end type basbas_distribution

    type custom_communicator
        INTEGER, ALLOCATABLE,DIMENSION(:) :: ranks_of_comm_global
        INTEGER :: n_ranks_of_comm_global
        INTEGER :: mpi_group
        INTEGER :: mpi_comm
        INTEGER :: n_tasks
        INTEGER :: myid
    end type custom_communicator

    INTEGER :: CUSTOM_COMMUNICATOR_UNASSIGNED
    PARAMETER(CUSTOM_COMMUNICATOR_UNASSIGNED = -1)

!storage of runtime information
    type(exp_coeff_int_distribution),POINTER, DIMENSION(:) :: cur_exp_coeff_int_dist
    type(basbas_distribution),POINTER, DIMENSION(:) :: resp_func_k_point_worklist

    INTEGER :: cur_exp_coeff_int_dist_size

    LOGICAL :: is_response_func_buffer_inited = .FALSE.

CONTAINS


    SUBROUTINE create_message_matrix_send(dim1,dim2,dim3,amessage)
        INTEGER, INTENT(IN) :: dim1,dim2,dim3
        type(message_matrix_send), INTENT(INOUT) :: amessage

        amessage%dim1=dim1
        amessage%dim2=dim2
        amessage%dim3=dim3
        ALLOCATE(amessage%data(dim1*dim2*dim3))

        if(.NOT. allocated(amessage%data)) then
            stop "failed to create message 3d matrix"
        endif

        amessage%data=0

        amessage%tag=0
        amessage%i_cell1=0
        amessage%i_cell2=0
        amessage%i_cell3=0
    END SUBROUTINE create_message_matrix_send

    SUBROUTINE create_message_matrix_send_from_matrix(dim1,dim2,dim3,mat, message_mat)
        INTEGER, INTENT(IN) :: dim1,dim2,dim3
        REAL*8, DIMENSION(dim1,dim2,dim3), INTENT(IN) :: mat
        type(message_matrix_send), INTENT(OUT) :: message_mat

        CALL create_message_matrix_send(dim1,dim2,dim3,message_mat)

        message_mat%data = reshape(mat(:,:,:), (/dim1*dim2*dim3/))

    END SUBROUTINE create_message_matrix_send_from_matrix

    SUBROUTINE create_message_matrix_send_array(n_dim,array)
        INTEGER, INTENT(IN) :: n_dim
        type(message_matrix_send), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: array

        ALLOCATE(array(n_dim))
        if(.NOT. allocated(array)) then
            stop "cannot allocate message_matrix_send_array"
        endif
!TODO check allocation
    END SUBROUTINE create_message_matrix_send_array

    SUBROUTINE free_message_matrix_send(amessage)
        type(message_matrix_send), INTENT(INOUT) :: amessage

        DEALLOCATE(amessage%data)
    END SUBROUTINE free_message_matrix_send

    SUBROUTINE free_message_matrix_send_array(array)
        type(message_matrix_send), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: array
!TODO: check for allcoted elements
        DEALLOCATE(array)
    END SUBROUTINE free_message_matrix_send_array


    SUBROUTINE copy_message_matrix_send(sour,dest)
        type(message_matrix_send), INTENT(IN) :: sour
        type(message_matrix_send), INTENT(INOUT) :: dest

        dest%source = sour%source
        dest%sources_index = sour%sources_index
        dest%dim1=sour%dim1
        dest%dim2=sour%dim2
        dest%dim2=sour%dim3
        dest%data=sour%data
    END SUBROUTINE copy_message_matrix_send

    SUBROUTINE create_message_matrix_recv(dim1,dim2,dim3,amessage)
        INTEGER, INTENT(IN) :: dim1,dim2,dim3
        type(message_matrix_recv), INTENT(INOUT) :: amessage

        amessage%dim1=dim1
        amessage%dim2=dim2
        amessage%dim3=dim3
        ALLOCATE(amessage%data(dim1,dim2,dim3))

        if(.NOT. allocated(amessage%data)) then
            stop "failed to create message 3d matrix - recv"
        endif

        amessage%data=0
    END SUBROUTINE create_message_matrix_recv

    SUBROUTINE clone_message_matrix_recv(source,dest)
        type(message_matrix_recv), INTENT(IN) :: source
        type(message_matrix_recv), INTENT(INOUT) :: dest

        CALL create_message_matrix_recv(source%dim1,source%dim2,source%dim3, dest)

        dest%source = source%source
        dest%sources_index = source%sources_index
        dest%dim1=source%dim1
        dest%dim2=source%dim2
        dest%dim2=source%dim3
        dest%data=source%data
    END SUBROUTINE clone_message_matrix_recv


    SUBROUTINE free_message_matrix_recv(amessage)
        type(message_matrix_recv), INTENT(INOUT) :: amessage

        if(allocated(amessage%data)) then
            DEALLOCATE(amessage%data)
        endif
    END SUBROUTINE free_message_matrix_recv

!exp coeff distribution
    SUBROUTINE set_exp_coeff_distribution(n_dim,adistribution)
        INTEGER, INTENT(IN) :: n_dim
        type(exp_coeff_int_distribution), DIMENSION(n_dim), INTENT(IN),TARGET :: adistribution

        cur_exp_coeff_int_dist => adistribution
        cur_exp_coeff_int_dist_size = n_dim

    END SUBROUTINE set_exp_coeff_distribution

    SUBROUTINE get_exp_coeff_distribution_by_cpu(n_cpu,adistribution)
        INTEGER, INTENT(IN) :: n_cpu
        type(exp_coeff_int_distribution), INTENT(OUT) :: adistribution

        adistribution = cur_exp_coeff_int_dist(n_cpu+1)
    END SUBROUTINE get_exp_coeff_distribution_by_cpu

    SUBROUTINE release_exp_coeff_distribution()
        NULLIFY (cur_exp_coeff_int_dist)
        cur_exp_coeff_int_dist_size = 0
    END SUBROUTINE release_exp_coeff_distribution

!basis func distribution
    SUBROUTINE create_basisfunc_distribution(n_start, n_end, rvecs, adistribution)
        INTEGER, INTENT(IN) :: n_start, n_end
        type(realspace_vectors), INTENT(IN) :: rvecs
        type(basisfunc_distribution), INTENT(OUT) :: adistribution

        adistribution%basisfunc_start = n_start
        adistribution%basisfunc_end = n_end
        adistribution%basisfunc_extend = n_end - n_start + 1

        adistribution%my_basis_start = 1
        adistribution%my_basis_end = adistribution%basisfunc_extend


        CALL clone_realspace_vectors(rvecs,adistribution%rvecs)

        CALL set_start_and_end_atom_for_basisfunc_distribution(adistribution)
!Might be memory intensive - comment: security reason
!        CALL set_basisfunc_distribution_atoms(adistribution)

    END SUBROUTINE create_basisfunc_distribution

    SUBROUTINE set_start_and_end_atom_for_basisfunc_distribution(dist)

        type(basisfunc_distribution), INTENT(INOUT) :: dist
        INTEGER :: atom_start, atom_end, i_atom

        atom_start = -1
        atom_end = -1
        do i_atom = 1,n_atoms
            if(atom_start == -1) then
                if(atom2basis_off(i_atom) >= dist%basisfunc_start) then
                    atom_start = i_atom -1
                endif
            endif

            if(atom_end== -1) then
                if(atom2basis_off(i_atom) >= dist%basisfunc_end) then
                    atom_end = i_atom-1
                endif
            endif
        enddo

        if(atom_start == -1) atom_start = n_atoms
        if(atom_end == -1) atom_end = n_atoms

        dist%atom_start = atom_start
        dist%atom_end = atom_end

    END SUBROUTINE set_start_and_end_atom_for_basisfunc_distribution

    SUBROUTINE set_basisfunc_distribution_atoms(dist)
        type(basisfunc_distribution), INTENT(INOUT) :: dist

        INTEGER :: atom_start, atom_end, i_atom

        atom_start = dist%atom_start
        atom_end = dist%atom_end

!TODO check alloc
        ALLOCATE(dist%atoms(atom_start:atom_end))

        if(atom_start == atom_end) then
            CALL set_basisfunc_distribution_atoms_one_atom(dist,atom_start)
        endif

        if(atom_start/=atom_end) then
            CALL set_basisfunc_distribution_atoms_multiple_atom(dist, &
                                                                atom_start,atom_end)
        endif

    END SUBROUTINE set_basisfunc_distribution_atoms

    SUBROUTINE set_basisfunc_distribution_atoms_one_atom(dist,i_atom)
        type(basisfunc_distribution), INTENT(INOUT) :: dist
        INTEGER, INTENT(IN) :: i_atom

        INTEGER :: part_of_atom_start, part_of_atom_end

        part_of_atom_start = dist%basisfunc_start - atom2basis_off(i_atom)
        part_of_atom_end = dist%basisfunc_end - atom2basis_off(i_atom)

        CALL set_basisfunc_distribution_atoms_overlapping_atom(dist,part_of_atom_start, &
                                                               part_of_atom_end, &
                                                               i_atom)

    END SUBROUTINE set_basisfunc_distribution_atoms_one_atom

    SUBROUTINE set_basisfunc_distribution_atoms_multiple_atom(dist,atom_start,atom_end)
        type(basisfunc_distribution), INTENT(INOUT) :: dist
        INTEGER, INTENT(IN) :: atom_start, atom_end

        INTEGER :: i_atom, part_of_atom_start, part_of_atom_end

        part_of_atom_start = dist%basisfunc_start-atom2basis_off(atom_start)
        part_of_atom_end = atom2basis_off(atom_start+1)-atom2basis_off(atom_start)

        CALL set_basisfunc_distribution_atoms_overlapping_atom(dist,part_of_atom_start, &
                                                               part_of_atom_end, &
                                                               atom_start)


        do i_atom = atom_start+1, atom_end-1
            part_of_atom_start = 1
            part_of_atom_end = atom2basis_off(i_atom+1)-atom2basis_off(i_atom)

            CALL set_basisfunc_distribution_atoms_overlapping_atom(dist,part_of_atom_start, &
                                                                   part_of_atom_end, &
                                                                   i_atom)
        enddo

        part_of_atom_start = 1
        part_of_atom_end = dist%basisfunc_end - atom2basis_off(atom_end)

        CALL set_basisfunc_distribution_atoms_overlapping_atom(dist,part_of_atom_start, &
                                                               part_of_atom_end, &
                                                               atom_end)
    END SUBROUTINE set_basisfunc_distribution_atoms_multiple_atom

    SUBROUTINE set_basisfunc_distribution_atoms_overlapping_atom(dist, &
                                                                 part_of_atom_start, &
                                                                 part_of_atom_end, &
                                                                 i_atom)
        type(basisfunc_distribution), INTENT(INOUT) :: dist
        INTEGER, INTENT(IN) :: part_of_atom_start, part_of_atom_end, &
                               i_atom

        INTEGER :: run_atom, run_vector, n_overlaping_atoms
        INTEGER, DIMENSION(n_atoms*dist%rvecs%n_vecs) :: overlaping_atoms
        INTEGER, DIMENSION(n_atoms*dist%rvecs%n_vecs) :: overlaping_atoms_supercell
        REAL*8,  DIMENSION(3,n_atoms*dist%rvecs%n_vecs) :: overlaping_atoms_connection_vector
        REAL*8 :: sum_radius
        REAL*8, DIMENSION(3) :: connection_vector

        n_overlaping_atoms = 0
        do run_vector = 1, dist%rvecs%n_vecs
            do run_atom = 1,n_atoms

!               if(is_zero_vector(dist%rvecs%vecs(run_vector)) &
!                  .AND. &
!                  run_atom == i_atom) then
!                  CYCLE
!                endif

               connection_vector = -1.d0 * coords(:,i_atom) + &
                                   coords(:,run_atom) + dist%rvecs%vecs(run_vector)%vector

               sum_radius = atom_radius(species(i_atom)) &
                            + atom_radius(species(run_atom))

!CALL write_debug("atom radi and sum:"// num2str(atom_radius(species(i_atom))) //&
!                 num2str(atom_radius(species(run_atom))) //&
!                 num2str(sum_radius))

!TODO TURN ON
!                if(vector_length(connection_vector) <= sum_radius * 1000) then
               if(vector_length(connection_vector) <= sum_radius) then
                  n_overlaping_atoms = n_overlaping_atoms + 1
                  overlaping_atoms(n_overlaping_atoms) = run_atom
                  overlaping_atoms_supercell(n_overlaping_atoms) = run_vector
                  overlaping_atoms_connection_vector(:,n_overlaping_atoms) = &
                                                                    connection_vector
               endif
            enddo
        enddo

        CALL create_basisfunc_distribution_atom(atom2basis_off(i_atom)-dist%basisfunc_start +1, &
                                                part_of_atom_start,part_of_atom_end, &
                                                n_overlaping_atoms, &
                                                overlaping_atoms,overlaping_atoms_supercell, &
                                                overlaping_atoms_connection_vector, &
                                                dist%atoms(i_atom))

    CONTAINS

        REAL*8 FUNCTION vector_length(vec)
            REAL*8,DIMENSION(3), INTENT(IN) :: vec
            vector_length = sqrt(sum(vec(1:3)**2))
        END FUNCTION vector_length

    END SUBROUTINE set_basisfunc_distribution_atoms_overlapping_atom

    SUBROUTINE free_basisfunc_distribution(adistribution)
        type(basisfunc_distribution), INTENT(INOUT) :: adistribution

        INTEGER :: i

        do i=adistribution%atom_start, adistribution%atom_end
            CALL free_basisfunc_distribution_atom(adistribution%atoms(i))
        enddo

        CALL free_all_lattice_vectors(adistribution%rvecs)

    END SUBROUTINE free_basisfunc_distribution

    SUBROUTINE print_basisfunc_distribution(adist)
        type(basisfunc_distribution) :: adist
        INTEGER :: i,i_atom,n_active_supercells


        CALL write_stdout("Basis start                   "// &
                          "end                  " //&
                          "extend            " // &
                          "Atom start            " // &
                          "end")

        CALL write_stdout(num2str(adist%basisfunc_start)// &
                          num2str(adist%basisfunc_end)// &
                          num2str(adist%basisfunc_extend) //&
                          num2str(adist%atom_start)// &
                          num2str(adist%atom_end))

        CALL write_stdout("Having overlapping atoms")
RETURN !debug
        do i_atom = adist%atom_start, adist%atom_end
            n_active_supercells = 0
            do i=1,adist%rvecs%n_vecs
                if(ANY(adist%atoms(i_atom)%overlaping_atoms_supercell == i)) then
                    n_active_supercells = n_active_supercells +1
                endif
            enddo

            CALL write_stdout("   For atom"//num2str(i_atom)//" overlaps with "&
                              //num2str(adist%atoms(i_atom)%n_overlaping_atoms) &
                              //" atoms in "//num2str(n_active_supercells)//"supercells")
            CALL write_stdout("       "//&
                              "Species basis start    "//&
                              "end          "//&
                              "my basis start          "// &
                              "end")

            CALL write_stdout(num2str(adist%atoms(i_atom)%part_of_atom_start)//&
                              num2str(adist%atoms(i_atom)%part_of_atom_end)//&
                              num2str(adist%atoms(i_atom)%my_basis_atom_start)//&
                              num2str(adist%atoms(i_atom)%my_basis_atom_end))

        end do

        CALL write_stdout("")
        CALL write_stdout("")
    END SUBROUTINE print_basisfunc_distribution

    SUBROUTINE create_basisfunc_distribution_atom(atom_offset, &
                                                  part_of_atom_start,part_of_atom_end, &
                                                  n_overlaping_atoms, &
                                                  overlaping_atoms,overlaping_atoms_supercell, &
                                                  overlaping_atoms_connection_vector, &
                                                  an_atom)

        INTEGER,INTENT(IN) :: atom_offset, &
                              part_of_atom_start, &
                              part_of_atom_end, &
                              n_overlaping_atoms

        INTEGER, DIMENSION(n_overlaping_atoms), INTENT(IN) :: overlaping_atoms
        INTEGER, DIMENSION(n_overlaping_atoms), INTENT(IN) :: overlaping_atoms_supercell
        REAL*8, DIMENSION(3,n_overlaping_atoms) :: overlaping_atoms_connection_vector

        type(basisfunc_distribution_atom), INTENT(OUT) :: an_atom

        INTEGER :: i

!TODO check alloc
        ALLOCATE(an_atom%overlaping_atoms(n_overlaping_atoms))
        ALLOCATE(an_atom%overlaping_atoms_supercell(n_overlaping_atoms))
        ALLOCATE(an_atom%overlaping_atoms_connection_vector(3,n_overlaping_atoms))

        an_atom%part_of_atom_start = part_of_atom_start
        an_atom%part_of_atom_end = part_of_atom_end

        an_atom%my_basis_atom_start =  atom_offset + part_of_atom_start
        an_atom%my_basis_atom_end = atom_offset + part_of_atom_end


        an_atom%n_overlaping_atoms = n_overlaping_atoms
        an_atom%overlaping_atoms(:) = overlaping_atoms(:)
        an_atom%overlaping_atoms_supercell(:) = overlaping_atoms_supercell(:)
        an_atom%overlaping_atoms_connection_vector(:,:) = overlaping_atoms_connection_vector(:,:)

    END SUBROUTINE create_basisfunc_distribution_atom

    SUBROUTINE free_basisfunc_distribution_atom(an_atom)
       type(basisfunc_distribution_atom), INTENT(INOUT) :: an_atom

       DEALLOCATE(an_atom%overlaping_atoms)
       DEALLOCATE(an_atom%overlaping_atoms_supercell)
       DEALLOCATE(an_atom%overlaping_atoms_connection_vector)
    END SUBROUTINE free_basisfunc_distribution_atom

!green's function distribution

    SUBROUTINE create_realspace_vector_distribution(my_vectors, &
                                               n_supercells, &
                                               vector_coefficients_to_cpu, &
                                               adistribution)

       type(realspace_vectors)  :: my_vectors
       INTEGER, DIMENSION(3,2), INTENT(IN) :: n_supercells
       INTEGER, DIMENSION(n_supercells(1,1):n_supercells(1,2), &
                          n_supercells(2,1):n_supercells(2,2), &
                          n_supercells(3,1):n_supercells(3,2)) :: vector_coefficients_to_cpu
       type(realspace_vector_distribution), INTENT(OUT) :: adistribution

       INTEGER :: i,n_my_vectors

       n_my_vectors = my_vectors%n_vecs

!TODO check alloc
       ALLOCATE(adistribution%vec_coeffs_to_cpu(n_supercells(1,1):n_supercells(1,2), &
                                                n_supercells(2,1):n_supercells(2,2), &
                                                n_supercells(3,1):n_supercells(3,2)))

       CALL clone_realspace_vectors(my_vectors,adistribution%my_vecs)
       adistribution%n_my_vecs = n_my_vectors

       adistribution%vec_coeffs_to_cpu(:,:,:) = vector_coefficients_to_cpu(:,:,:)

    END SUBROUTINE create_realspace_vector_distribution

    SUBROUTINE free_realspace_vector_distribution(adistribution)
       type(realspace_vector_distribution), INTENT(INOUT) :: adistribution

       CALL free_all_lattice_vectors(adistribution%my_vecs)
       DEALLOCATE(adistribution%vec_coeffs_to_cpu)
    END SUBROUTINE free_realspace_vector_distribution

!custom communicator function

    SUBROUTINE create_custom_communicator(n_ranks_of_comm_global, ranks_of_comm_global, &
                                          a_custom_communicator)
       INTEGER, INTENT(IN) :: n_ranks_of_comm_global
       INTEGER, DIMENSION(:),INTENT(IN) :: ranks_of_comm_global
       type(custom_communicator), INTENT(OUT) :: a_custom_communicator

!TODO check alloc
       ALLOCATE(a_custom_communicator%ranks_of_comm_global(n_ranks_of_comm_global))
       a_custom_communicator%n_ranks_of_comm_global = n_ranks_of_comm_global

       a_custom_communicator%ranks_of_comm_global(:) = ranks_of_comm_global(1:n_ranks_of_comm_global)

       a_custom_communicator%mpi_group = CUSTOM_COMMUNICATOR_UNASSIGNED
       a_custom_communicator%mpi_comm = CUSTOM_COMMUNICATOR_UNASSIGNED
       a_custom_communicator%n_tasks = CUSTOM_COMMUNICATOR_UNASSIGNED
       a_custom_communicator%myid = CUSTOM_COMMUNICATOR_UNASSIGNED
    END SUBROUTINE create_custom_communicator

    SUBROUTINE assign_custom_communicator(group_custom_communicator, &
                                          comm_custom_communicator, &
                                          n_tasks_custom_communicator, &
                                          rank_custom_communicator, &
                                          a_custom_communicator)
       INTEGER, INTENT(IN) :: group_custom_communicator, &
                              comm_custom_communicator, &
                              rank_custom_communicator, &
                              n_tasks_custom_communicator

       type(custom_communicator), INTENT(INOUT) :: a_custom_communicator

       a_custom_communicator%mpi_group = group_custom_communicator
       a_custom_communicator%mpi_comm = comm_custom_communicator
       a_custom_communicator%n_tasks = n_tasks_custom_communicator
       a_custom_communicator%myid = rank_custom_communicator

    END SUBROUTINE assign_custom_communicator

    SUBROUTINE print_custom_communicator(a_custom_communicator)
       type(custom_communicator), INTENT(IN) :: a_custom_communicator

       CALL write_debug("Custom communicator")
       CALL write_debug("group:"//num2str(a_custom_communicator%mpi_group))
       CALL write_debug("comm:"//num2str(a_custom_communicator%mpi_comm))
       CALL write_debug("n_tasks:"//num2str(a_custom_communicator%n_tasks))
       CALL write_debug("myid:"//num2str(a_custom_communicator%myid))
       CALL write_debug("rank global:"//num2str(a_custom_communicator%ranks_of_comm_global))
    END SUBROUTINE print_custom_communicator

    SUBROUTINE free_custom_communicator(a_custom_communicator)
       type(custom_communicator), INTENT(INOUT) :: a_custom_communicator

       DEALLOCATE(a_custom_communicator%ranks_of_comm_global)
    END SUBROUTINE free_custom_communicator

!integer array functions

    logical FUNCTION has_array_this_element(inarray,element)
        INTEGER, DIMENSION(:), INTENT(IN) :: inarray
        INTEGER, INTENT(IN) :: element

        has_array_this_element = (ANY(inarray==element))

    END FUNCTION has_array_this_element

    SUBROUTINE replace_first_zero_with_element(inarray,arraysize,element)
        INTEGER, INTENT(IN) :: arraysize
        INTEGER, DIMENSION(arraysize), INTENT(INOUT) :: inarray
        INTEGER, INTENT(IN) :: element

        INTEGER :: i

        do i = 1, arraysize
            if(inarray(i)==0) then
                inarray(i)=element
                return
            end if
        end do
    END SUBROUTINE replace_first_zero_with_element

    SUBROUTINE add_to_set_if_not_included(inarray,arraysize,element)
        INTEGER, INTENT(IN) :: arraysize
        INTEGER, DIMENSION(arraysize), INTENT(INOUT) :: inarray
        INTEGER, INTENT(IN) :: element

        if(.NOT. has_array_this_element(inarray,element)) then
           CALL replace_first_zero_with_element(inarray,arraysize, element)
        endif


    END SUBROUTINE add_to_set_if_not_included

    INTEGER FUNCTION size_of_set(inarray,arraysize)
        INTEGER, INTENT(IN) :: arraysize
        INTEGER, DIMENSION(arraysize), INTENT(IN) :: inarray

        INTEGER :: i, num_elements

        num_elements = 0

        do i = 1, arraysize
            if(inarray(i) /=0) then
                num_elements=num_elements+1
            endif
        end do

        size_of_set = num_elements

    END FUNCTION size_of_set

!point functions
    LOGICAL FUNCTION is_point_equal(first_point,second_point)
        type(point), INTENT(IN) :: first_point,second_point

        if(first_point%x == second_point%x .AND. first_point%y == second_point%y) then
            is_point_equal = .true.
            return
        endif

        is_point_equal= .false.
    END FUNCTION


    LOGICAL FUNCTION is_point_zero(apoint)
        type(point), INTENT(IN) :: apoint

        is_point_zero= is_point_equal(apoint,point(0,0))
    END FUNCTION

    LOGICAL FUNCTION is_coordinate_in_point(coordinate,apoint)
        INTEGER, INTENT(IN) :: coordinate
        type(point), INTENT(IN) :: apoint

        if(apoint%x == coordinate .OR. apoint%y == coordinate) then
            is_coordinate_in_point = .true.
            return
        endif

        is_coordinate_in_point= .false.
    END FUNCTION

!matrix functions

    SUBROUTINE create_matrixblock(dim_row,dim_column,amatrixblock)
        INTEGER, INTENT(IN) :: dim_row,dim_column

        type(matrixblock_of_point), INTENT(OUT) :: amatrixblock

        !local variables
        INTEGER :: i,j

        amatrixblock%dim_row = dim_row
        amatrixblock%dim_column = dim_column
        allocate( amatrixblock%element(dim_row,dim_column) )

        if(.NOT. ALLOCATED(amatrixblock%element)) then
            write(use_unit,*) "Failed to allocate memory"
            stop
        endif

        ! init with zeros
        FORALL(i=1:dim_row, j=1:dim_column) &
                amatrixblock%element(i,j) = point(0,0)

    END SUBROUTINE create_matrixblock

    SUBROUTINE free_matrixblock(amatrixblock)
        type(matrixblock_of_point), INTENT(INOUT) :: amatrixblock

        deallocate(amatrixblock%element)

    END SUBROUTINE free_matrixblock

    INTEGER FUNCTION count_of_nonzero_elements(amatrixblock)
        type(matrixblock_of_point), INTENT(IN) :: amatrixblock

        INTEGER :: i,j,num_of_nonzeros

        num_of_nonzeros = 0

        do i = 1, amatrixblock%dim_row
            do j = 1, amatrixblock%dim_column
                if (.NOT. is_point_zero(amatrixblock%element(i,j))) then
                    num_of_nonzeros=num_of_nonzeros+1
                end if
            end do
        end do

        count_of_nonzero_elements = num_of_nonzeros
    END FUNCTION count_of_nonzero_elements

    SUBROUTINE copy_matrix_to_row (amatrixblock,offset_in_row,amatrixrow)
        type(matrixblock_of_point), INTENT(IN) :: amatrixblock
        INTEGER, INTENT(IN) :: offset_in_row
        type(basbas_distribution), INTENT(INOUT) :: amatrixrow


        !local variable
        INTEGER :: i,j,pos_in_array
        pos_in_array = offset_in_row+1

        do j = 1, amatrixblock%dim_column
            do i = 1, amatrixblock%dim_row
                if(.NOT. is_point_zero(amatrixblock%element(i,j)) ) then
                    amatrixrow%element(pos_in_array) = amatrixblock%element(i,j)
                    pos_in_array = pos_in_array+1
                end if
            end do
        end do

    END SUBROUTINE copy_matrix_to_row

    SUBROUTINE extract_row_from_matrix(num_row,amatrixblock,amatrixrow)
        INTEGER, INTENT(IN) :: num_row
        type(matrixblock_of_point), INTENT(IN) :: amatrixblock
        type(basbas_distribution), INTENT(OUT) :: amatrixrow

        INTEGER :: i

        CALL create_basbas_distribution(amatrixblock%dim_column,amatrixrow)

        do i = 1, amatrixblock%dim_column
            amatrixrow%element(i) =point(amatrixblock%element(num_row,i)%x, &
                                         amatrixblock%element(num_row,i)%y)
        end do


    END SUBROUTINE extract_row_from_matrix

    SUBROUTINE extract_minimal_row_from_matrix(num_row,amatrixblock,amatrixrow)
        INTEGER, INTENT(IN) :: num_row
        type(matrixblock_of_point), INTENT(IN) :: amatrixblock
        type(basbas_distribution), INTENT(OUT) :: amatrixrow

        type(basbas_distribution) :: afull_matrixrow

        CALL extract_row_from_matrix(num_row,amatrixblock,afull_matrixrow)
        CALL minimize_matrixrow(afull_matrixrow,amatrixrow)
        CALL free_matrixrow(afull_matrixrow)

    END SUBROUTINE extract_minimal_row_from_matrix


    SUBROUTINE print_matrixblock(ablock)
        type(matrixblock_of_point), INTENT(IN) :: ablock

        INTEGER :: i,j
        INTEGER, DIMENSION(ablock%dim_row*ablock%dim_column) ::  distinct_elements

        write(use_unit,*) "Matrix content row,column,i,j"

        do j = 1, ablock%dim_row
            do i = 1, ablock%dim_column
                write(use_unit,*) i,j,ablock%element(i,j)%x, ablock%element(i,j)%y
            end do
        end do

    END SUBROUTINE print_matrixblock

! row functions

    SUBROUTINE create_basbas_distribution(dim_row,amatrixrow)
        INTEGER, INTENT(IN) :: dim_row

        type(basbas_distribution), INTENT(OUT) :: amatrixrow

        !local variables
        INTEGER :: i

        amatrixrow%dim_row = dim_row
        allocate( amatrixrow%element(dim_row) )

        if(.NOT. ALLOCATED(amatrixrow%element)) then
            write(use_unit,*) "Failed to allocate memory"
            stop
        endif

        ! init with zeros
        FORALL(i=1:dim_row) &
                amatrixrow%element(i) = point(0,0)

    END SUBROUTINE create_basbas_distribution

    SUBROUTINE copy_matrixrow(source,dest)
        type(basbas_distribution), INTENT(IN) :: source
        type(basbas_distribution), INTENT(INOUT) :: dest

        INTEGER :: i

        if(ALLOCATED(dest%element)) THEN
            CALL free_matrixrow(dest)
        endif

        CALL create_basbas_distribution(source%dim_row,dest)

        do i = 1, source%dim_row
            dest%element(i) = point(source%element(i)%x, &
                                    source%element(i)%y)
        enddo

        dest%n_rows = source%n_rows
        dest%n_columns = source%n_columns
        dest%n_offset_rows = source%n_offset_rows
        dest%n_offset_columns = source%n_offset_columns

    END SUBROUTINE copy_matrixrow

    SUBROUTINE free_matrixrow(amatrixblock)
        type(basbas_distribution), INTENT(INOUT) :: amatrixblock

        deallocate(amatrixblock%element)

    END SUBROUTINE free_matrixrow

    INTEGER FUNCTION row_count_of_nonzero_elements(amatrixrow)
        type(basbas_distribution), INTENT(IN) :: amatrixrow

        INTEGER :: i,num_of_nonzeros

        num_of_nonzeros = 0

        do i = 1, amatrixrow%dim_row
            if (.NOT. is_point_zero(amatrixrow%element(i))) then
               num_of_nonzeros=num_of_nonzeros+1
            end if
        end do

        row_count_of_nonzero_elements = num_of_nonzeros
    END FUNCTION row_count_of_nonzero_elements


    SUBROUTINE minimize_matrixrow(amatrixrow,aminimal_row)
        type(basbas_distribution), INTENT(IN) :: amatrixrow
        type(basbas_distribution), INTENT(OUT) :: aminimal_row

        INTEGER :: i,offset

        CALL create_basbas_distribution(row_count_of_nonzero_elements(amatrixrow), aminimal_row)

        offset = 1

        do i = 1, amatrixrow%dim_row
            if( .NOT. is_point_zero(amatrixrow%element(i)) ) then
                aminimal_row%element(offset)=point(amatrixrow%element(i)%x , &
                                                   amatrixrow%element(i)%y)

                offset = offset + 1
            end if
        end do

    END SUBROUTINE minimize_matrixrow

    LOGICAL FUNCTION has_row_element(arow,element)
        type(basbas_distribution),INTENT(IN) :: arow
        type(point), INTENT(IN) :: element

        INTEGER :: i

        do i = 1, arow%dim_row
            if(is_point_equal(arow%element(i),element)) then
                has_row_element = .true.
                return
            end if
        end do

        has_row_element = .false.
    END FUNCTION

    SUBROUTINE get_distinct_field_elements(amatrixrow,distinct_elements, max_elements)
        INTEGER, INTENT(IN) :: max_elements
        type(basbas_distribution), INTENT(IN) :: amatrixrow
        INTEGER, DIMENSION(max_elements), INTENT(OUT) :: distinct_elements

        CALL get_distinct_field_elements_with_offset(amatrixrow,distinct_elements, &
                                                     max_elements,1,amatrixrow%dim_row)
    END SUBROUTINE

    SUBROUTINE get_distinct_field_elements_with_offset(amatrixrow,distinct_elements, max_elements,offset,length)
        INTEGER, INTENT(IN) :: max_elements
        type(basbas_distribution), INTENT(IN) :: amatrixrow
        INTEGER, DIMENSION(max_elements), INTENT(OUT) :: distinct_elements
        INTEGER, INTENT(IN) :: offset,length

        !local variables
        INTEGER :: i,j,an_element
        distinct_elements(:)=0

        if(offset + length > amatrixrow%dim_row) then
            stop "get_distinct_field_elements_with_offset will exceed dim row"
        endif

        do i = offset, offset + length !amatrixrow%dim_row
                an_element = amatrixrow%element(i)%x
                CALL add_to_set_if_not_included(distinct_elements,max_elements, an_element)

                an_element = amatrixrow%element(i)%y
                CALL add_to_set_if_not_included(distinct_elements,max_elements, an_element)
        end do
    END SUBROUTINE get_distinct_field_elements_with_offset

    SUBROUTINE print_basbas_distributions(n_dists,dists)
        INTEGER, INTENT(IN) :: n_dists
        type(basbas_distribution),DIMENSION(n_dists), INTENT(IN) :: dists

        INTEGER :: i

        if(myid/=0) return

        write(use_unit,*) "Basbas distribution"

        write(use_unit,*) "                          ", &
                   "start        stop        extend"

        do i = 1,n_dists

            write(use_unit,*) "Task ", i, ": ", &
                       "", &
                       dists(i)%n_offset_rows+1, &
                       "", &
                       dists(i)%n_offset_rows+dists(i)%n_rows, &
                       "", &
                       dists(i)%n_rows

        enddo
    END SUBROUTINE print_basbas_distributions

    SUBROUTINE print_basbas_distribution(dist)
        type(basbas_distribution), INTENT(IN) :: dist

        write(use_unit,*) "             ", &
                   "start: ", dist%n_offset_rows+1,"stop: ",dist%n_offset_rows+dist%n_rows, &
                   "extend: ", dist%n_rows

    END SUBROUTINE print_basbas_distribution


END MODULE crpa_parallelism_storage
