MODULE cRPA_flow_test
    USE cRPA_flow_realspace_vector

    implicit none

CONTAINS

    SUBROUTINE add_vector_to_periodic_vector_test()
        type(realspace_vectors) :: vectors
        REAL*8, DIMENSION(3,3) :: lattice_vector
        INTEGER :: i,j,d, vec_index
        type(realspace_vector) :: vec_first, vec_second, vec_sum

        write(use_unit,*) "Start testing realspace vector distribution"

        lattice_vector = reshape((/ 0.5d0,0.0d0,0.0d0 ,0.0d0,0.5d0,0.0d0, 0.0d0,0.0d0,0.5d0/), (/3,3/))

        CALL create_all_realspace_vectors_sorted(n_cells_bvk, cell_index_bvk,lattice_vector,vectors)

        do d=1,6

            CALL find_vector_in_vector_array_by_coeffs(0,0,0, &
                                                       vectors, &
                                                       vec_index)
            vec_first = vectors%vecs(vec_index)

            SELECT CASE (d)
                CASE (1)
                   CALL find_vector_in_vector_array_by_coeffs(1,0,0, &
                                                              vectors, &
                                                             vec_index)
                CASE (2)
                   CALL find_vector_in_vector_array_by_coeffs(0,1,0, &
                                                              vectors, &
                                                             vec_index)
                CASE (3)
                   CALL find_vector_in_vector_array_by_coeffs(0,0,-1, &
                                                              vectors, &
                                                             vec_index)
                CASE (4)
                   CALL find_vector_in_vector_array_by_coeffs(1,0,1, &
                                                              vectors, &
                                                             vec_index)
                CASE (5)
                   CALL find_vector_in_vector_array_by_coeffs(1,1,1, &
                                                              vectors, &
                                                             vec_index)

                CASE (6)
                   CALL find_vector_in_vector_array_by_coeffs(-1,1,-1, &
                                                              vectors, &
                                                             vec_index)
            END SELECT

            vec_second = vectors%vecs(vec_index)

            do i=1,10

                CALL add_vector_to_periodic_vector(vec_second,vec_first, &
                                                   vectors, &
                                                   vec_sum)

                CALL print_vector(vec_sum)
                vec_first = vec_sum
            enddo
         enddo

        write(use_unit,*) "testing realspace vector distribution: ok"
    END SUBROUTINE add_vector_to_periodic_vector_test

    SUBROUTINE get_point_inversion_of_periodic_vector_test()
        type(realspace_vectors) :: vectors,periodic_vectors
        REAL*8, DIMENSION(3,3) :: lattice_vector
        INTEGER :: i,j,d, vec_index
        type(realspace_vector) :: vec, vec_inversion

        write(use_unit,*) "Start testing realspace vector distribution"


        lattice_vector = reshape((/ 0.5d0,0.0d0,0.0d0 ,0.0d0,0.5d0,0.0d0, 0.0d0,0.0d0,0.5d0/), (/3,3/))

        CALL create_all_realspace_vectors_sorted(n_cells_bvk, cell_index_bvk,lattice_vector,periodic_vectors)

        lattice_vector = reshape((/ 0.5d0,0.0d0,0.0d0 ,0.0d0,0.5d0,0.0d0, 0.0d0,0.0d0,0.5d0/), (/3,3/))

        CALL create_all_realspace_vectors_sorted(n_cells_bvk, cell_index_bvk,lattice_vector,vectors)


        do d=0,6
            SELECT CASE (d)
                CASE (0)
                   CALL find_vector_in_vector_array_by_coeffs(1,0,0, &
                                                              vectors, &
                                                             vec_index)

                CASE (1)
                   CALL find_vector_in_vector_array_by_coeffs(5,0,0, &
                                                              vectors, &
                                                             vec_index)
                CASE (2)
                   CALL find_vector_in_vector_array_by_coeffs(0,1,0, &
                                                              vectors, &
                                                             vec_index)
                CASE (3)
                   CALL find_vector_in_vector_array_by_coeffs(0,0,-1, &
                                                              vectors, &
                                                             vec_index)
                CASE (4)
                   CALL find_vector_in_vector_array_by_coeffs(-4,0,1, &
                                                              vectors, &
                                                             vec_index)
                CASE (5)
                   CALL find_vector_in_vector_array_by_coeffs(1,7,1, &
                                                              vectors, &
                                                             vec_index)

                CASE (6)
                   CALL find_vector_in_vector_array_by_coeffs(-6,2,-4, &
                                                              vectors, &
                                                             vec_index)
            END SELECT

            vec = vectors%vecs(vec_index)
            CALL get_point_inversion_of_periodic_vector(vec,periodic_vectors,vec_inversion)


            CALL write_debug("Normal vec")
            CALL print_vector(vec)
            CALL write_debug("Inversion vec")
            CALL print_vector(vec_inversion)
        enddo

    END SUBROUTINE get_point_inversion_of_periodic_vector_test

    SUBROUTINE get_realspace_vectors_index_inversion_test()
        type(realspace_vectors) :: rvecs
        INTEGER, ALLOCATABLE, DIMENSION(:) :: index_list
        INTEGER :: i_vec,i_vec_inversion, i_vec_double_inversion
        type(realspace_vector) :: vec_cur, vec_inversion, vec_inversion_by_list, &
                                  vec_double_inversion

        write(use_unit,*) "Testing get_realspace_vectors_index_inversion"

        CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk, lattice_vector, &
                                          rvecs)

        ALLOCATE(index_list(rvecs%n_vecs))
        CALL get_realspace_vectors_index_inversion(rvecs, index_list)

        do i_vec = 1, rvecs%n_vecs
            vec_cur = rvecs%vecs(i_vec)
            CALL get_point_inversion_of_periodic_vector(vec_cur, rvecs, vec_inversion)
            vec_inversion_by_list = rvecs%vecs(index_list(i_vec))

            if(.NOT. are_equal_vectors(vec_inversion, vec_inversion_by_list)) then
                write(use_unit,*) "vec_cur"
                CALL print_vector(vec_cur)
                write(use_unit,*) "vec_inversion_by_list"
                CALL print_vector(vec_inversion_by_list)
                write(use_unit,*) "vec_inversion"
                CALL print_vector(vec_inversion)

                stop "FAILED! vec inversion and vec_inversion by list are unequal!"
            endif

            CALL find_vector_in_vector_array(vec_inversion_by_list, rvecs, &
                                             i_vec_inversion)
            i_vec_double_inversion = index_list(i_vec_inversion)
            vec_double_inversion = rvecs%vecs(i_vec_double_inversion)

            if(.NOT. are_equal_vectors(vec_cur, vec_double_inversion)) then
                write(use_unit,*) "vec_cur"
                CALL print_vector(vec_cur)
                write(use_unit,*) "vec_double_inversion"
                CALL print_vector(vec_double_inversion)
                stop "FAILED! vec double inversion and vec_cur are unequal!"
            endif
        enddo

        DEALLOCATE(index_list)

        write(use_unit,*) "Test get_realspace_vectors_index_inversion: ok"
    END SUBROUTINE get_realspace_vectors_index_inversion_test

    SUBROUTINE get_realspace_vectors_index_all_possible_differences_test()
        type(realspace_vectors) :: rvecs
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: index_list
        INTEGER :: i_vec,i_vec_subst, i_vec_diff
        type(realspace_vector) :: vec_cur, vec_subst, vec_diff, &
                                  vec_diff_by_list, &
                                  vec_zero

        write(use_unit,*) "Testing get_realspace_vectors_index_all_possible_differences"

        CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk, lattice_vector, &
                                          rvecs)

        ALLOCATE(index_list(rvecs%n_vecs,rvecs%n_vecs))
        CALL get_realspace_vectors_index_all_possible_differences(rvecs, index_list)

        do i_vec = 1, rvecs%n_vecs
            do i_vec_subst = 1, rvecs%n_vecs
                vec_cur = rvecs%vecs(i_vec)
                CALL get_point_inversion_of_periodic_vector(rvecs%vecs(i_vec_subst), &
                                                            rvecs, vec_subst)
                CALL add_vector_to_periodic_vector(vec_subst,vec_cur,rvecs, &
                                                   vec_diff)

                vec_diff_by_list = rvecs%vecs(index_list(i_vec,i_vec_subst))

                if(.NOT. are_equal_vectors(vec_diff, vec_diff_by_list)) then
                    write(use_unit,*) "vec_cur"
                    CALL print_vector(vec_cur)
                    write(use_unit,*) "vec_subst"
                    CALL print_vector(vec_subst)
                    write(use_unit,*) "vec_diff_by_list"
                    CALL print_vector(vec_diff_by_list)
                    write(use_unit,*) "vec_diff"
                    CALL print_vector(vec_diff)

                    stop "FAILED! vec diff and vec_diff by list are unequal!"
                endif

                CALL find_vector_in_vector_array(vec_diff_by_list, rvecs, &
                                                 i_vec_diff)

                vec_zero = rvecs%vecs(index_list(i_vec_diff,i_vec_diff))

                if(.NOT. is_zero_vector(vec_zero)) then
                    write(use_unit,*) "vec_diff"
                    CALL print_vector(vec_diff)
                    write(use_unit,*) "vec_zero"
                    CALL print_vector(vec_zero)
                    stop "FAILED! vec zero is not zero!"
                endif
            enddo
        enddo

        DEALLOCATE(index_list)

        write(use_unit,*) "Test get_realspace_vectors_index_all_possible_differences: ok"
    END SUBROUTINE get_realspace_vectors_index_all_possible_differences_test
END MODULE cRPA_flow_test
