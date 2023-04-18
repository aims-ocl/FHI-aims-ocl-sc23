!Calculates/creates the auxiliary exp basis for least square fit and ft
module cRPA_calculation_expbasis

    USE numerical_utilities
    USE cRPA_calculation_integration
    USE cRPA_calculation_fouriertrans

    IMPLICIT NONE


    type exp_basis
        INTEGER :: n
        REAL*8,ALLOCATABLE,DIMENSION(:) :: exponents, &
                                           evaluated_points, &
                                           int_weights
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: transform_matrix, &
                                               evaluated_at_points, &
                                               orthonormal_rep, &
                                               evaluated_at_points_onb
        LOGICAL :: is_ft_prepared
        INTEGER :: n_evaluated_points
    end type exp_basis
CONTAINS

    SUBROUTINE create_exp_basis(n,exponents, an_exp_basis)
        INTEGER, INTENT(IN) :: n
        REAL*8, DIMENSION(n) :: exponents
        type(exp_basis) :: an_exp_basis


!TODO check alloc
        ALLOCATE(an_exp_basis%exponents(n))
        an_exp_basis%n = n
        an_exp_basis%exponents(:) = exponents(:)
        an_exp_basis%is_ft_prepared = .FALSE.

        CALL bubblesort_in_place(an_exp_basis%n, &
                                 an_exp_basis%exponents(:))

        CALL set_orthonormal_coefficients(an_exp_basis)

    END SUBROUTINE create_exp_basis

    SUBROUTINE set_orthonormal_coefficients(an_exp_basis)
        type(exp_basis) :: an_exp_basis
        INTEGER :: i,j, i_basis_normal, i_basis_new
        REAL*8, DIMENSION(an_exp_basis%n) :: tau, work
        REAL*8, DIMENSION(an_exp_basis%n,an_exp_basis%n) :: T, Tinv, Q,id
        INTEGER, DIMENSION(an_exp_basis%n) :: ipiv
        INTEGER :: lda, lwork, info, n

        ALLOCATE(an_exp_basis%orthonormal_rep(an_exp_basis%n, &
                                              an_exp_basis%n))

        n = an_exp_basis%n

        an_exp_basis%orthonormal_rep(:,:) = 0.d0

        do i_basis_normal = 1,an_exp_basis%n
            an_exp_basis%orthonormal_rep(i_basis_normal,i_basis_normal) = &
                            1.d0
        enddo
return
        do i_basis_normal = 1,an_exp_basis%n
            an_exp_basis%orthonormal_rep(i_basis_normal,:) = &
                  an_exp_basis%orthonormal_rep(i_basis_normal,:)/norm_orthogonal_vec(an_exp_basis,i_basis_normal)

            do i_basis_new= i_basis_normal+1, an_exp_basis%n
                an_exp_basis%orthonormal_rep(i_basis_new,:) = &
                            an_exp_basis%orthonormal_rep(i_basis_new,:) &
                            - &
                            project_along_vec(an_exp_basis,i_basis_normal,i_basis_new) * &
                            an_exp_basis%orthonormal_rep(i_basis_normal,:)
            enddo
        enddo

    END SUBROUTINE set_orthonormal_coefficients

    REAL*8 FUNCTION scalprod_exp(an_exp_basis,i_normal_1, i_normal_2)
        type(exp_basis) :: an_exp_basis
        INTEGER, INTENT(IN) :: i_normal_1,i_normal_2
        scalprod_exp = -1.d0 / &
                       (an_exp_basis%exponents(i_normal_1) + &
                        an_exp_basis%exponents(i_normal_2))
    END FUNCTION

     REAL*8 FUNCTION scalprod_orthogonal_vec(an_exp_basis,i_orthogonal1,i_orthogonal2)
         type(exp_basis) :: an_exp_basis
         INTEGER, INTENT(IN) :: i_orthogonal1,i_orthogonal2
         REAL*8, DIMENSION(an_exp_basis%n**2) :: part
         INTEGER :: i,j

         do i = 1, an_exp_basis%n
             do j = 1, an_exp_basis%n
                 part((i-1)*an_exp_basis%n+j) = an_exp_basis%orthonormal_rep(i_orthogonal1,j) &
                                 * an_exp_basis%orthonormal_rep(i_orthogonal2,i) &
                                 * scalprod_exp(an_exp_basis,i,j)
            enddo
         enddo

          scalprod_orthogonal_vec = sort_and_sum_vector(an_exp_basis%n**2, part)
      END FUNCTION


     REAL*8 FUNCTION norm_orthogonal_vec(an_exp_basis,i_orthogonal)
         type(exp_basis) :: an_exp_basis
         INTEGER, INTENT(IN) :: i_orthogonal
         REAL*8, DIMENSION(i_orthogonal**2) :: norm_part
         INTEGER :: i,j

         norm_orthogonal_vec = sqrt(scalprod_orthogonal_vec(an_exp_basis,i_orthogonal,i_orthogonal))
      END FUNCTION

      REAL*8 FUNCTION project_along_vec(an_exp_basis,i_direction, i_vec)
        type(exp_basis) :: an_exp_basis
        INTEGER,INTENT(IN) :: i_direction, i_vec

        project_along_vec = scalprod_orthogonal_vec(an_exp_basis,i_direction, i_vec) &
                              /&
                            norm_orthogonal_vec(an_exp_basis,i_direction)
      END FUNCTION

    SUBROUTINE evaluate_exp_basis(x,an_exp_basis,h_x)
        REAL*8,INTENT(IN) :: x
        type(exp_basis) :: an_exp_basis
        REAL*8, DIMENSION(an_exp_basis%n),INTENT(OUT) :: h_x

        h_x(:) = exp(an_exp_basis%exponents(:) * x)

    END SUBROUTINE evaluate_exp_basis

    SUBROUTINE evaluate_exp_basis_on_points(i,x,an_exp_basis,h_x)
        INTEGER, INTENT(IN) :: i
        REAL*8,DIMENSION(i), INTENT(IN) :: x
        type(exp_basis) :: an_exp_basis
        REAL*8, DIMENSION(i,an_exp_basis%n),INTENT(OUT) :: h_x

        INTEGER :: j

        do j=1,i
            CALL evaluate_exp_basis(x(j),an_exp_basis, &
                                              h_x(j,:))
        enddo
    END SUBROUTINE evaluate_exp_basis_on_points

    SUBROUTINE evaluate_onb_basis(x,an_exp_basis,h_x)
        REAL*8,INTENT(IN) :: x
        type(exp_basis) :: an_exp_basis
        REAL*8, DIMENSION(an_exp_basis%n),INTENT(OUT) :: h_x

        INTEGER :: i
        REAL*8, DIMENSION(an_exp_basis%n) :: parts

        do i=1,an_exp_basis%n
            parts(:) = an_exp_basis%orthonormal_rep(i,:) * &
                         exp(an_exp_basis%exponents(:) * x)

            h_x(i) = sort_and_sum_vector(an_exp_basis%n, parts(:))
        enddo

!write(use_unit,*) "eval_onb", x, h_x(1:min(an_exp_basis%n,3)), h_x(an_exp_basis%n)
    END SUBROUTINE evaluate_onb_basis

    SUBROUTINE evaluate_onb_basis_on_points(i,x,an_exp_basis,h_x)
        INTEGER, INTENT(IN) :: i
        REAL*8,DIMENSION(i), INTENT(IN) :: x
        type(exp_basis) :: an_exp_basis
        REAL*8, DIMENSION(i,an_exp_basis%n),INTENT(OUT) :: h_x

        INTEGER :: j

        do j=1,i
            CALL evaluate_onb_basis(x(j),an_exp_basis, &
                                              h_x(j,:))
        enddo
    END SUBROUTINE evaluate_onb_basis_on_points


    SUBROUTINE get_exp_interpolation_matrix(n,x_n, an_exp_basis, A)
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN), DIMENSION(n) :: x_n
        type(exp_basis) :: an_exp_basis
        REAL*8, INTENT(OUT), DIMENSION(an_exp_basis%n, n) :: A

        REAL*8, DIMENSION(n) :: work
        REAL*8, DIMENSION(an_exp_basis%n, n) :: tmp
        INTEGER, DIMENSION(n) :: ipiv
        INTEGER :: lda, lwork, info

        if(an_exp_basis%n == n) then
            CALL evaluate_exp_basis_on_points(n, x_n, &
                                              an_exp_basis, A)

            lda = n
            lwork = n
            A(:,:)=transpose(A(:,:))
            call dgetrf( n, n, A, lda, ipiv, info )
            call dgetri( n, A, lda, ipiv, work, lwork, info )
        else
            CALL pseudo_inverse('cRPA get_exp_interpolation matrix', &
                           an_exp_basis%n, &
                           n, &
                           A, &
                           tmp, &
                           0.d-7)
            A(:,:) = tmp(:,:)
        endif

    END SUBROUTINE get_exp_interpolation_matrix

    SUBROUTINE prepare_fouriertransform_by_exp_basis(n_points, points, an_exp_basis, int_weights)
       INTEGER, INTENT(IN) :: n_points
       REAL*8, DIMENSION(n_points), INTENT(IN) :: points
       type(exp_basis) :: an_exp_basis
       REAL*8, DIMENSION(n_points), INTENT(IN) :: int_weights


       if(an_exp_basis%is_ft_prepared) then
          stop "Already prepared"
       endif

       ALLOCATE(an_exp_basis%evaluated_points(n_points))
       an_exp_basis%evaluated_points(:) = points
       an_exp_basis%n_evaluated_points = n_points


       ALLOCATE(an_exp_basis%evaluated_at_points(an_exp_basis%n,n_points))
       ALLOCATE(an_exp_basis%evaluated_at_points_onb(an_exp_basis%n,n_points))

       CALL evaluate_exp_basis_on_points(n_points, points, &
                                         an_exp_basis, &
                                         an_exp_basis%evaluated_at_points)

       CALL evaluate_onb_basis_on_points(n_points, points, &
                                         an_exp_basis, &
                                         an_exp_basis%evaluated_at_points_onb)

       ALLOCATE(an_exp_basis%int_weights(n_points))
       an_exp_basis%int_weights(:) = int_weights(:)

       an_exp_basis%is_ft_prepared = .TRUE.

    END SUBROUTINE

    SUBROUTINE fouriertransform_by_exp_basis(n_points, points, &
                                             n_elements, data, &
                                             omega, &
                                             an_exp_basis, &
                                             res)
       INTEGER, INTENT(IN) :: n_points
       REAL*8, DIMENSION(n_points), INTENT(IN) :: points
       INTEGER, INTENT(IN) :: n_elements
       REAL*8, DIMENSION(n_elements, n_points), INTENT(IN) :: data
       REAL*8, INTENT(IN) :: omega
       type(exp_basis) :: an_exp_basis
       REAL*8, DIMENSION(n_elements), INTENT(OUT) :: res


       REAL*8, DIMENSION(an_exp_basis%n) :: rep_vector, &
                                            ft_of_basis
       INTEGER :: i_element,i_part,i_parts, i_start_points, i_end_points, i_extend_points, &
                  i_start_basis, i_end_basis, i_extend_basis

       if(.NOT. an_exp_basis%is_ft_prepared) then
            stop "Prepare ft first!!!"
       endif

       ft_of_basis(:) = abs(an_exp_basis%exponents(:)) &
                            / &
                        ( (omega*1.d0)**2 + abs(an_exp_basis%exponents(:))**2)
       res(:) = 0.d0
       i_parts = 1
       do i_part = 1,i_parts
          i_start_points = 1 + ((dble(i_part)-1.d0)/dble(i_parts)) * n_points
          i_end_points = (dble(i_part)/dble(i_parts)) * n_points
          i_extend_points = n_points !i_end_points - i_start_points + 1
 
          i_start_basis = 1 + ((dble(i_part)-1.d0)/dble(i_parts)) * an_exp_basis%n
          i_end_basis = (dble(i_part)/dble(i_parts)) * an_exp_basis%n
          i_extend_basis = an_exp_basis%n !i_end_basis - i_start_basis + 1

!          write(use_unit,*) i_part, i_start_points,i_end_points, i_extend_points
!          write(use_unit,*) i_part, i_start_basis, i_end_basis, i_extend_basis
       do i_element = 1,n_elements

          CALL solve_LSQ('cRPA resp one omega', &
                         i_extend_points, &
                         i_extend_basis, &
                         an_exp_basis%evaluated_at_points, & !(i_start_points:i_end_points, i_start_basis:i_end_basis), &
                         data(i_element,:), & !i_start_points:i_end_points), &
                         rep_vector)

!if(myid==0 .AND. i_element == 1) write(use_unit,*) rep_vector(1:i_extend_basis)
          res(i_element) = res(i_element) + &
                           sum_stable(i_extend_basis,rep_vector(:) * ft_of_basis(:))
       enddo
       enddo

    END SUBROUTINE fouriertransform_by_exp_basis

    SUBROUTINE fouriertransform_by_onb_exp_basis(n_points, points, &
                                             n_elements, data, &
                                             omega, &
                                             an_exp_basis, &
                                             res)
       INTEGER, INTENT(IN) :: n_points
       REAL*8, DIMENSION(n_points), INTENT(IN) :: points
       INTEGER, INTENT(IN) :: n_elements
       REAL*8, DIMENSION(n_elements, n_points), INTENT(IN) :: data
       REAL*8, INTENT(IN) :: omega
       type(exp_basis) :: an_exp_basis
       REAL*8, DIMENSION(n_elements), INTENT(OUT) :: res


       REAL*8, DIMENSION(an_exp_basis%n) :: rep_vector, &
                                            ft_of_basis_exp, &
                                            ft_of_basis_onb, &
                                            exp_integrals, &
                                            tmp, tmp2
       REAL*8, DIMENSION(n_points)::tmp_res


       INTEGER :: i_element, i

       if(.NOT. an_exp_basis%is_ft_prepared) then
            stop "Prepare ft first!!!"
       endif

       ft_of_basis_exp(:) = abs(an_exp_basis%exponents(:)) &
                            / &
                           ( (omega*1.d0)**2 + abs(an_exp_basis%exponents(:))**2)

       do i=1,an_exp_basis%n
          tmp(:) = an_exp_basis%orthonormal_rep(i,:) &
                     * ft_of_basis_exp(:)

          ft_of_basis_onb(i) = sort_and_sum_vector(an_exp_basis%n, tmp)
       enddo



       res(:) = 0.d0

       do i_element = 1,n_elements

          CALL solve_LSQ('cRPA resp one omega', &
                         n_points, &
                         an_exp_basis%n, &
                         an_exp_basis%evaluated_at_points_onb, &
                         data(i_element,:), &
                         rep_vector)
          res(i_element) = sum_stable(an_exp_basis%n,rep_vector(:) * ft_of_basis_onb(:))

!          do i=1,an_exp_basis%n
!              tmp_res(:) = data(i_element,:) &
!                           *exp(points(:) * an_exp_basis%exponents(i)) &
!                           * an_exp_basis%int_weights(:)
!              exp_integrals(i) = sum_stable(n_points,tmp_res)
!          enddo
! tmp2 = rep_vector
!          do i=1,an_exp_basis%n
!              tmp(:) = an_exp_basis%orthonormal_rep(i,:) * exp_integrals(:)
!              rep_vector(i) = sum_stable(n_points,tmp)
!          enddo

!          tmp_res(:) = data(i_element,:) * an_exp_basis%int_weights(:)
!          tmp_res(:) = an_exp_basis%evaluated_at_points_onb(1,:) * an_exp_basis%int_weights(:)

!          tmp(:) = rep_vector

!          do i=1,an_exp_basis%n
!            tmp_res(:) = an_exp_basis%evaluated_at_points_onb(i,:) * &
!                         data(i_element,:) * an_exp_basis%int_weights(:)
!            rep_vector(i) = sort_and_sum_vector(n_points,tmp_res)
!          enddo

!          rep_vector = MATMUL(an_exp_basis%evaluated_at_points_onb,tmp_res)

!write(use_unit,*) "rep_vector", tmp2(1:3), rep_vector(1:3)


!          tmp(1) = sum_stable(an_exp_basis%n,rep_vector(:) * ft_of_basis_onb(:))

!          write(use_unit,*) "res", res(i_element), tmp(1)
!          res(i_element) = tmp(1)

       enddo

    END SUBROUTINE fouriertransform_by_onb_exp_basis

    SUBROUTINE fouriertransform_by_integration(n_points, points, &
                                             n_elements, data, &
                                             omega, &
                                             an_exp_basis, &
                                             res)
       INTEGER, INTENT(IN) :: n_points
       REAL*8, DIMENSION(n_points), INTENT(IN) :: points
       INTEGER, INTENT(IN) :: n_elements
       REAL*8, DIMENSION(n_elements, n_points), INTENT(IN) :: data
       REAL*8, INTENT(IN) :: omega
       type(exp_basis) :: an_exp_basis
       REAL*8, DIMENSION(n_elements), INTENT(OUT) :: res

       COMPLEX*16,DIMENSION(n_elements) :: tmp_res
       type(integration_grid) :: int_grid

       if(.NOT. an_exp_basis%is_ft_prepared) then
            stop "Prepare ft first!!!"
       endif

       CALL create_integration_grid(n_points, points, an_exp_basis%int_Weights, &
                                    int_Grid)

       CALL fouriertransform_x_to_t(int_grid, &
                                    n_elements, data, &
                                    omega, tmp_res)

       res(:) = dble(tmp_res(:))

       CALL free_integration_grid(int_grid)
    END SUBROUTINE fouriertransform_by_integration


    SUBROUTINE print_exp_basis(an_exp_basis)
        type(exp_basis) :: an_exp_basis

        INTEGER :: i,j

        CALL write_stdout("Printing exponents of this basis")

        do i=1,an_exp_basis%n
            CALL write_stdout(num2str(i)//num2str(an_exp_basis%exponents(i)))
        enddo

return

        CALL write_stdout("Printing orthonormal coeffs in this basis")
        do i=1,min(8,an_exp_basis%n)
            do j=1,min(8,an_exp_basis%n)
                CALL write_stdout(num2str(i)//num2str(j) &
                                  //num2str(an_exp_basis%orthonormal_rep(i,j)))
            enddo
        enddo

        CALL write_stdout("Printing norms in this basis")
        do i=1,min(8,an_exp_basis%n)
             CALL write_stdout(num2str(i) &
                               //num2str(norm_orthogonal_vec(an_exp_basis,i)))
        enddo


        CALL write_stdout("Printing scalprods in this basis")
        do i=1,min(8,an_exp_basis%n)
            do j=1,min(8,an_exp_basis%n)
                CALL write_stdout(num2str(i)//num2str(j) &
                                  //num2str(scalprod_orthogonal_vec(an_exp_basis,i,j)))
            enddo
        enddo


    END SUBROUTINE print_exp_basis

    SUBROUTINE free_exp_basis(an_exp_basis)
        type(exp_basis) :: an_exp_basis

        DEALLOCATE(an_exp_basis%exponents)
        DEALLOCATE(an_exp_basis%orthonormal_rep)

        if(an_exp_basis%is_ft_prepared) then
            DEALLOCATE(an_exp_basis%evaluated_points)
            DEALLOCATE(an_exp_basis%evaluated_at_points)
            DEALLOCATE(an_exp_basis%evaluated_at_points_onb)
            DEALLOCATE(an_exp_basis%int_weights)

            an_exp_basis%is_ft_prepared = .FALSE.
        endif
    END SUBROUTINE free_exp_basis

    SUBROUTINE get_at_most_n_distinct_energies(n_energies,energies, &
                                               min_energy, max_energy, &
                                               n, &
                                               n_mean_energies, mean_energies)
        INTEGER, INTENT(IN) :: n_energies
        REAL*8, INTENT(IN), DIMENSION(n_energies) :: energies
        REAL*8, INTENT(IN) :: min_energy,max_energy
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: n_mean_energies
        REAL*8, INTENT(OUT), DIMENSION(:), ALLOCATABLE :: mean_energies

        INTEGER :: i, i_pos
        REAL*8 :: last_added
        REAL*8, DIMENSION(n) :: final_energies
        REAL*8, DIMENSION(n_energies) :: tmp_energies
        REAL*8 :: class_r

        CALL write_stdout("Looking in "//num2str(n_energies)//"energies for good energies")
!        CALL get_n_good_mean_energies(n_energies, energies, min(n,n_energies), tmp_energies)

        tmp_energies(:) = energies(:)
        CALL bubblesort_in_place(n_energies,tmp_energies)
        CALL write_stdout("Energies sorted")

        final_energies(1) =0.d0! tmp_energies(1)
        last_added = 1.d0 !final_energies(1)
!---- change to i_pos =2
        i_pos = 1

        class_r = 0.001d0

        do i = 1,n_energies
            if( tmp_energies(i) < min_energy) then
               cycle
            endif

            if( tmp_energies(i) > max_energy) then
               cycle
            endif

            if(abs(tmp_energies(i)) < class_r) then
               class_r = class_r / 10.d0
            endif

            if (abs(tmp_energies(i)-last_added)<class_r) then
               CYCLE
            endif

            final_energies(i_pos) = tmp_energies(i)
            last_added = final_energies(i_pos)
            i_pos = i_pos + 1

            if (i_pos>n) exit
        enddo

        n_mean_energies = i_pos-1
        ALLOCATE(mean_energies(n_mean_energies))
        mean_energies = final_energies(1:n_mean_energies)

        if(n_mean_energies < n) then
            CALL write_stdout("Had to condense energies: "//num2str(n_mean_energies) &
                              //"<"//num2str(n))
        else
            CALL write_stdout("Choose"//num2str(n_mean_energies) &
                              //"energies")
        endif

    END SUBROUTINE

    SUBROUTINE get_n_good_mean_energies(n_energies,energies,n, mean_energies)
        INTEGER, INTENT(IN) :: n_energies
        REAL*8, INTENT(IN), DIMENSION(n_energies) :: energies
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(OUT), DIMENSION(n) :: mean_energies

        REAL*8, DIMENSION(n_energies) :: sorted_energies
        INTEGER :: i_categorie, i_energies, i_start,i_end
        type(integration_intervals) :: int_intervals
        REAL*8 :: int_start,int_end
        INTEGER :: i

        sorted_energies(:) = energies(:)

        if(n>n_energies) then
            stop "There must be at least n energies!"
        endif

        CALL bubblesort_in_place(n_energies, sorted_energies)

        CALL create_integration_intervals_logarithmic(0.d0, &
                                                      maxval(abs(sorted_energies(:))), &
                                                      n, int_intervals)

        do i_categorie = 1,n

            int_start = -int_intervals%interval(i_categorie)%start
            int_end = -int_intervals%interval(i_categorie)%end
!write(use_unit,*) int_start,int_end
            i_start = -1
            i_end = -1

            do i =1, n_energies-1
                if(sorted_energies(i+1)>int_start .AND. i_end <= 0) i_end = i
                if(sorted_energies(i)>int_end .AND. i_start <= 0) i_start = i
            enddo

            if(i_start == -1) i_start = 1
            if(i_end == -1) i_end = n_energies

            i_energies = i_end - i_start

!write(use_unit,*) i_categorie,i_energies, i_start,i_end
            mean_energies(i_categorie) = SUM(sorted_energies(i_start:i_end)) &
                                            / &
                                         dble(i_energies)

        enddo

    END SUBROUTINE get_n_good_mean_energies

    SUBROUTINE bubblesort_in_place(n_data,data)
        INTEGER,INTENT(IN) :: n_data
        REAL*8, INTENT(INOUT), DIMENSION(n_data) :: data

        INTEGER :: n,i
        REAL*8 :: tmp

!bubblesort
        do n=n_data, 1,-1
            do i=1,n-1
                if (data(i) > data(i+1)) then
                   tmp = data(i+1)
                   data(i+1) = data(i)
                   data(i) = tmp
                endif
            enddo
        enddo
    END SUBROUTINE

    SUBROUTINE get_minkowski_difference(n_vec1,vec1,n_vec2,vec2, vec_diff, i_vec_size)
        INTEGER, INTENT(IN) :: n_vec1
        REAL*8, DIMENSION(n_vec1),INTENT(IN) :: vec1
        INTEGER, INTENT(IN) :: n_vec2
        REAL*8, DIMENSION(n_vec2),INTENT(IN) :: vec2
        REAL*8, DIMENSION(n_vec1*n_vec2), INTENT(OUT) :: vec_diff
        INTEGER,INTENT(OUT) :: i_vec_size

        REAL*8 :: cur_diff
        INTEGER :: i_vec1, i_vec2, i_pos

        vec_diff(1) = 0.d0
        i_pos = 2
        do i_vec1 = 1,n_vec1
            do i_vec2 = 1,n_vec2
                cur_diff = vec1(i_vec1) - vec2(i_vec2)

                if(cur_diff<0.d0) then
                    vec_diff(i_pos) = cur_diff
                    i_pos = i_pos +1
                endif
            enddo
        enddo

        i_vec_size = i_pos - 1
    END SUBROUTINE
end module cRPA_calculation_expbasis
