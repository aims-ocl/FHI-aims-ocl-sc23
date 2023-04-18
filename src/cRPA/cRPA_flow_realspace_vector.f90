!Implements realspace lattice vectors (cf. BvK translation vectors)
MODULE cRPA_flow_realspace_vector
    USE pbc_lists
    USE cRPA_view
    implicit none

    type realspace_vector
        REAL*8, DIMENSION(3) :: vector
        REAL*8 :: norm
        INTEGER :: i_cell1, i_cell2,i_cell3
        INTEGER :: i_cell_bvk
    end type realspace_vector

    type realspace_vectors
        type(realspace_vector),ALLOCATABLE, DIMENSION(:) :: vecs
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: coeff_index
        INTEGER, DIMENSION(3,2) :: n_supercells
        INTEGER :: n_vecs
        INTEGER :: sorting_order
    end type

    INTEGER :: SORTING_ORDER_AIMS_BVK
    PARAMETER(SORTING_ORDER_AIMS_BVK=1)

    INTEGER :: SORTING_ORDER_BY_NORM
    PARAMETER(SORTING_ORDER_BY_NORM=2)

CONTAINS

!realspace vector functions
    SUBROUTINE create_all_realspace_vectors_sorted(n_cells_bvk, cell_index_bvk, &
                                                   lattice_vector,all_vecs_sorted)
        INTEGER, INTENT(IN) :: n_cells_bvk
        INTEGER, DIMENSION(n_cells_bvk,3), INTENT(IN) :: cell_index_bvk
        REAL*8, DIMENSION(3,3), INTENT(IN) :: lattice_vector
        type(realspace_vectors), INTENT(OUT) :: all_vecs_sorted

        CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk, &
                                          lattice_vector,all_vecs_sorted)
        CALL sort_realspace_vectors_by_norm(all_vecs_sorted)
        CALL build_coeff_index(all_vecs_sorted)
        CALL build_supercell_entry(all_vecs_sorted)
    END SUBROUTINE create_all_realspace_vectors_sorted

    SUBROUTINE create_all_realspace_vectors(n_cells_bvk, cell_index_bvk, &
                                            lattice_vector,all_vecs)
        INTEGER, INTENT(IN) :: n_cells_bvk
        INTEGER, DIMENSION(n_cells_bvk,3), INTENT(IN) :: cell_index_bvk
        REAL*8, DIMENSION(3,3), INTENT(IN) :: lattice_vector
        type(realspace_vectors), INTENT(OUT) :: all_vecs

        INTEGER :: i,j,k,all_vec_index,m

        REAL*8, DIMENSION(3) :: avector
        ALLOCATE(all_vecs%vecs(n_cells_bvk))
!TODO check alloc
        all_vec_index = 1

        do m = 1, n_cells_bvk
          i=cell_index_bvk(m,1)
          j=cell_index_bvk(m,2)
          k=cell_index_bvk(m,3)


!        do i=-n_supercells(1),n_supercells(1)
!            do j=-n_supercells(2),n_supercells(2)
!                do k=-n_supercells(3),n_supercells(3)
!
                    avector = MATMUL(lattice_vector,(/i,j,k/))

                    all_vecs%vecs(all_vec_index)%vector = avector
                    all_vecs%vecs(all_vec_index)%norm = sqrt(dot_product(avector, avector))

                    all_vecs%vecs(all_vec_index)%i_cell1 = i
                    all_vecs%vecs(all_vec_index)%i_cell2 = j
                    all_vecs%vecs(all_vec_index)%i_cell3 = k
                    all_vecs%vecs(all_vec_index)%i_cell_bvk = m !cell_index_bvk(i,j,k)

                    all_vec_index=all_vec_index+1
!                enddo
!            enddo
        enddo

        all_vecs%n_vecs = all_vec_index - 1
        all_vecs%sorting_order = SORTING_ORDER_AIMS_BVK
        CALL build_coeff_index(all_vecs)
        CALL build_supercell_entry(all_vecs)
    END SUBROUTINE create_all_realspace_vectors

    SUBROUTINE build_coeff_index(rvecs)
        type(realspace_vectors), INTENT(INOUT) :: rvecs

        INTEGER,DIMENSION(3,2) :: n_supercells
        INTEGER :: i

!TODO check allocation
        if(ALLOCATED(rvecs%coeff_index)) then
            DEALLOCATE(rvecs%coeff_index)
        endif

        CALL get_supercells_spanned_by_vec_array(rvecs, n_supercells)

        ALLOCATE(rvecs%coeff_index(n_supercells(1,1):n_supercells(1,2), &
                                   n_supercells(2,1):n_supercells(2,2), &
                                   n_supercells(3,1):n_supercells(3,2)))

        rvecs%coeff_index(:,:,:) = -1

        do i = 1, rvecs%n_vecs
            rvecs%coeff_index(rvecs%vecs(i)%i_cell1, &
                              rvecs%vecs(i)%i_cell2, &
                              rvecs%vecs(i)%i_cell3) = i
        enddo

    END SUBROUTINE build_coeff_index

    SUBROUTINE build_supercell_entry(rvecs)
        type(realspace_vectors), INTENT(INOUT) :: rvecs

        CALL get_supercells_spanned_by_vec_array(rvecs,rvecs%n_supercells)
    END SUBROUTINE build_supercell_entry

    SUBROUTINE sort_realspace_vectors_by_norm(rvecs)
        type(realspace_vectors), INTENT(INOUT) :: rvecs

        type(realspace_vector) :: tmp
        INTEGER :: n,i
!bubblesort
        do n=rvecs%n_vecs, 1,-1
            do i=1,n-1
                if (rvecs%vecs(i)%norm > rvecs%vecs(i+1)%norm) then
                   CALL copy_realspace_vector(rvecs%vecs(i+1),tmp)
                   CALL copy_realspace_vector(rvecs%vecs(i),rvecs%vecs(i+1))
                   CALL copy_realspace_vector(tmp,rvecs%vecs(i))
                endif
            enddo
        enddo

        rvecs%sorting_order = SORTING_ORDER_BY_NORM
    END SUBROUTINE sort_realspace_vectors_by_norm

    SUBROUTINE get_max_radi_of_n_first_vectors(n_vecs,vecs,n,max_radi)
        INTEGER,INTENT(IN) :: n,n_vecs
        type(realspace_vector), DIMENSION(n_vecs), INTENT(IN) :: vecs
        REAL*8, INTENT(OUT) :: max_radi

        INTEGER :: i
        REAL*8 :: min_value

        min_value=0.d0

        if(n>=n_vecs) then
            max_radi = maxval(vecs(:)%norm,1)
        elseif(n<n_vecs) then
            do i=1,n
                min_value = minval(vecs(:)%norm,1, vecs(:)%norm > min_value)
            enddo

            max_radi = min_value
        endif

    END SUBROUTINE

    SUBROUTINE clone_realspace_vectors(source,clone)
        type(realspace_vectors), INTENT(IN) :: source
        type(realspace_vectors), INTENT(OUT) :: clone


        INTEGER :: i

!TODO check alloc
        ALLOCATE(clone%vecs(source%n_vecs))

        do i=1,source%n_vecs
            CALL copy_realspace_vector(source%vecs(i),clone%vecs(i))
        enddo

        clone%n_vecs = source%n_vecs
        CALL build_coeff_index(clone)
        CALL build_supercell_entry(clone)
    END SUBROUTINE clone_realspace_vectors

    SUBROUTINE copy_realspace_vector(source,dest)
        type(realspace_vector),INTENT(IN) :: source
        type(realspace_vector),INTENT(INOUT) :: dest

        dest%vector = source%vector
        dest%norm = source%norm

        dest%i_cell1 = source%i_cell1
        dest%i_cell2 = source%i_cell2
        dest%i_cell3 = source%i_cell3
        dest%i_cell_bvk = source%i_cell_bvk

    END SUBROUTINE copy_realspace_vector

    SUBROUTINE get_n_smallest_realspace_vectors(n,n_cells_bvk, cell_index_bvk, &
                                                lattice_vector,smallest_vectors)
        INTEGER, INTENT(IN) :: n,n_cells_bvk
        INTEGER, DIMENSION(n_cells_bvk,3), INTENT(IN) :: cell_index_bvk
        REAL*8, DIMENSION(3,3) :: lattice_vector
        type(realspace_vectors), INTENT(OUT) :: smallest_vectors

        type(realspace_vectors) :: all_vecs
        REAL*8 :: max_radius
        INTEGER :: n_vecs, i, vecs_index

        CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk,lattice_vector,all_vecs)
        n_vecs = all_vecs%n_vecs

        CALL sort_realspace_vectors_by_norm(all_vecs)
        ALLOCATE(smallest_vectors%vecs(min(n,n_vecs)))
!TODO check alloc
!        CALL get_max_radi_of_n_first_vectors(n_vecs,all_vecs,n,max_radius)

        do i = 1, min(n,n_vecs)
            CALL copy_realspace_vector(all_vecs%vecs(i), smallest_vectors%vecs(i))
        enddo

        CALL free_all_lattice_vectors(all_vecs)

    END SUBROUTINE get_n_smallest_realspace_vectors

    SUBROUTINE find_vector_in_vector_array(vector, vector_array, pos)
        type(realspace_vector), INTENT(IN) :: vector
        type(realspace_vectors), INTENT(IN) :: vector_array
        INTEGER, INTENT(OUT) :: pos


        CALL find_vector_in_vector_array_by_coeffs(vector%i_cell1, &
                                                   vector%i_cell2, &
                                                   vector%i_cell3, &
                                                   vector_array, pos)

!        INTEGER :: i
!
!        do i=1,vector_array%n_vecs
!            if(are_equal_vectors(vector_array%vecs(i),vector)) then
!                pos = i
!                return
!            endif
!        enddo
!
!        pos = -1

    END SUBROUTINE find_vector_in_vector_array

    SUBROUTINE find_vector_in_vector_array_by_coeffs(i_cell1,i_cell2,i_cell3, &
                                                     vector_array, pos)
        INTEGER, INTENT(IN) :: i_cell1, i_cell2, i_cell3
        type(realspace_vectors), INTENT(IN) :: vector_array
        INTEGER, INTENT(OUT) :: pos

        pos = vector_array%coeff_index(i_cell1,i_cell2,i_cell3)

if(pos/=vector_array%vecs(pos)%i_cell_bvk ) then
CALL write_debug( "pos unequal to i_cell_bvk"//num2str(pos)//num2str(vector_array%vecs(pos)%i_cell_bvk))
stop
endif

    END SUBROUTINE find_vector_in_vector_array_by_coeffs

    SUBROUTINE get_supercells_spanned_by_vec_array(rvecs,n_supercells)
        type(realspace_vectors), INTENT(IN) :: rvecs
        INTEGER, DIMENSION(3,2), INTENT(OUT) :: n_supercells

        n_supercells(1,1) = minval(rvecs%vecs(:)%i_cell1)
        n_supercells(1,2) = maxval(rvecs%vecs(:)%i_cell1)

        n_supercells(2,1) = minval(rvecs%vecs(:)%i_cell2)
        n_supercells(2,2) = maxval(rvecs%vecs(:)%i_cell2)

        n_supercells(3,1) = minval(rvecs%vecs(:)%i_cell3)
        n_supercells(3,2) = maxval(rvecs%vecs(:)%i_cell3)

    END SUBROUTINE

    LOGICAL FUNCTION are_equal_vectors(vector1,vector2)
        type(realspace_vector), INTENT(IN) :: vector1, vector2

        if( vector1%i_cell1 == vector2%i_cell1 .AND. &
            vector1%i_cell2 == vector2%i_cell2 .AND. &
            vector1%i_cell3 == vector2%i_cell3) THEN

            are_equal_vectors = .TRUE.
            return
        endif

       are_equal_vectors = .FALSE.

    END FUNCTION are_equal_vectors

    LOGICAL FUNCTION is_zero_vector(avector)
       type(realspace_vector), INTENT(IN) :: avector

        if( avector%i_cell1 == 0 .AND. &
            avector%i_cell2 == 0 .AND. &
            avector%i_cell3 == 0) THEN

            is_zero_vector = .TRUE.
            return
        endif

       is_zero_vector = .FALSE.

    END FUNCTION is_zero_vector

    SUBROUTINE free_all_lattice_vectors(all_vecs)
        type(realspace_vectors), INTENT(INOUT) :: all_vecs

        DEALLOCATE(all_vecs%vecs)
        DEALLOCATE(all_vecs%coeff_index)
    END SUBROUTINE free_all_lattice_vectors

    SUBROUTINE create_occupation_set_from_realspace_vectors(n_supercells,vecs,occupation_set)
        INTEGER, DIMENSION(3), INTENT(IN) :: n_supercells
        type(realspace_vectors), INTENT(IN) :: vecs

        LOGICAL, DIMENSION(-n_supercells(1):n_supercells(1), &
                           -n_supercells(2):n_supercells(2), &
                           -n_supercells(3):n_supercells(3))&
               , INTENT(OUT) :: occupation_set

        INTEGER :: i

        occupation_set(:,:,:) = .FALSE.

        do i=1, vecs%n_vecs
            occupation_set(vecs%vecs(i)%i_cell1, &
                           vecs%vecs(i)%i_cell2, &
                           vecs%vecs(i)%i_cell3) = .TRUE.
        enddo

    END SUBROUTINE create_occupation_set_from_realspace_vectors

    SUBROUTINE add_vector_to_periodic_vector(vector, periodic_vector, &
                                             all_periodic_vectors, &
                                             sum_vector)

    type(realspace_vector), INTENT(IN) :: vector, periodic_vector
    type(realspace_vectors), INTENT(IN) :: all_periodic_vectors
    type(realspace_vector), INTENT(OUT) :: sum_vector

    INTEGER, DIMENSION(3,2) :: n_supercell
    INTEGER, DIMENSION(3) :: n_supercell_extend
    INTEGER :: i_cell1,i_cell2,i_cell3, sum_vector_pos

!    CALL get_supercells_spanned_by_vec_array(all_periodic_vectors, &
!                                             n_supercell)

    n_supercell = all_periodic_vectors%n_supercells

    !sum vectors
    i_cell1 = vector%i_cell1 + periodic_vector%i_cell1
    i_cell2 = vector%i_cell2 + periodic_vector%i_cell2
    i_cell3 = vector%i_cell3 + periodic_vector%i_cell3

    !do mapping
    n_supercell_extend = n_supercell(:,2) - n_supercell(:,1) +1

    CALL map_dimension_into_supercell(n_supercell(1,:),n_supercell_extend(1), &
                                      i_cell1)

    CALL map_dimension_into_supercell(n_supercell(2,:),n_supercell_extend(2), &
                                      i_cell2)

    CALL map_dimension_into_supercell(n_supercell(3,:),n_supercell_extend(3), &
                                      i_cell3)

    CALL find_vector_in_vector_array_by_coeffs(i_cell1,i_cell2,i_cell3, &
                                               all_periodic_vectors, &
                                               sum_vector_pos)

    CALL copy_realspace_vector(all_periodic_vectors%vecs(sum_vector_pos), sum_vector)

    CONTAINS
        SUBROUTINE map_dimension_into_supercell(supercell_extend, supercell_width,i_cell)
            INTEGER,DIMENSION(2), INTENT(IN) :: supercell_extend
            INTEGER, INTENT(IN) :: supercell_width
            INTEGER, INTENT(INOUT) :: i_cell

            if(i_cell > supercell_extend(2)) then
                do while (i_cell > supercell_extend(2))
                    i_cell = i_cell - supercell_width
                end do
            endif

            if(i_cell < supercell_extend(1)) then
                do while (i_cell < supercell_extend(1))
                    i_cell = i_cell + supercell_width
                end do
            endif

        END SUBROUTINE map_dimension_into_supercell


    END SUBROUTINE add_vector_to_periodic_vector

    SUBROUTINE get_point_inversion_of_periodic_vector(vec,periodic_vectors,inversion_vec)
        type(realspace_vector), INTENT(IN) :: vec
        type(realspace_vectors), INTENT(IN) :: periodic_vectors
        type(realspace_vector), INTENT(OUT) :: inversion_vec

        type(realspace_vector) :: minus_vec

        minus_vec%i_cell1 = -2 * vec%i_cell1
        minus_vec%i_cell2 = -2 * vec%i_cell2
        minus_vec%i_cell3 = -2 * vec%i_cell3

        CALL add_vector_to_periodic_vector(minus_vec,vec,periodic_vectors, &
                                           inversion_vec)

    END SUBROUTINE get_point_inversion_of_periodic_vector

    SUBROUTINE print_vector(avector)
        type(realspace_vector) :: avector

        write(use_unit,*) "icell1", avector%i_cell1
        write(use_unit,*) "icell2", avector%i_cell2
        write(use_unit,*) "icell3", avector%i_cell3
        write(use_unit,*) "coords:", avector%vector
        write(use_unit,*) "Norm", avector%norm

    END SUBROUTINE print_vector

    SUBROUTINE convert_realspace_vecs_to_array_vecs(realspace_vecs,arrayvecs)
        type(realspace_vectors), INTENT(IN) :: realspace_vecs
        REAL*8,DIMENSION(3,realspace_vecs%n_vecs), INTENT(OUT) :: arrayvecs

        INTEGER :: i

        do i=1,realspace_vecs%n_vecs
            arrayvecs(:,i) = realspace_vecs%vecs(i)%vector(:)
        enddo
    END SUBROUTINE convert_realspace_vecs_to_array_vecs

    SUBROUTINE get_realspace_vectors_index_inversion(rvecs, index_list)
        type(realspace_vectors), INTENT(IN) :: rvecs
        INTEGER, DIMENSION(rvecs%n_vecs), INTENT(OUT) :: index_list

        type(realspace_vector) :: vec_cur, vec_inversion
        INTEGER :: i_vec

        do i_vec = 1,rvecs%n_vecs
            vec_cur = rvecs%vecs(i_vec)
            CALL get_point_inversion_of_periodic_vector(vec_cur, &
                                                        rvecs, &
                                                        vec_inversion)
            CALL find_vector_in_vector_array(vec_inversion, rvecs, &
                                             index_list(i_vec))
        enddo
    END SUBROUTINE get_realspace_vectors_index_inversion

    SUBROUTINE get_realspace_vectors_index_all_possible_differences(rvecs, index_list)
        type(realspace_vectors), INTENT(IN) :: rvecs
        !vec_1 - vec_2
        INTEGER, DIMENSION(rvecs%n_vecs,rvecs%n_vecs), INTENT(OUT) :: index_list

        type(realspace_vector) :: vec_cur, vec_subst, vec_diff
        INTEGER :: i_vec, i_vec_subst

        do i_vec = 1, rvecs%n_vecs
            vec_cur = rvecs%vecs(i_vec)
            do i_vec_subst = 1, rvecs%n_vecs

                CALL get_point_inversion_of_periodic_vector(rvecs%vecs(i_vec_subst), &
                                                            rvecs, &
                                                            vec_subst)
                CALL add_vector_to_periodic_vector(vec_subst,vec_cur, rvecs, &
                                                   vec_diff)

                CALL find_vector_in_vector_array(vec_diff, rvecs, &
                                                 index_list(i_vec,i_vec_subst))
            enddo
        enddo

    END SUBROUTINE get_realspace_vectors_index_all_possible_differences
END MODULE cRPA_flow_realspace_vector
