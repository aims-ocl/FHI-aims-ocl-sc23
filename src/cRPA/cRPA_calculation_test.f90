MODULE cRPA_calculation_test
    USE cRPA_calculation

    IMPLICIT NONE

CONTAINS

    SUBROUTINE calculate_imaginary_greens_func_comp_eigenvecs_test()
        REAL*8 :: tau
        REAL*8, dimension(n_states,n_spin,n_k_points)  :: KS_eigenvalue
        COMPLEX*16, dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector_complex
        REAL*8,  dimension(n_states, n_spin, n_k_points) :: occ_numbers
        COMPLEX*16,allocatable :: greens_func_complex(:,:,:,:)

        KS_eigenvalue = (1,1)
        KS_eigenvector_complex = 2.0d0
        occ_numbers = 1
        tau = 0.0d0


        ALLOCATE(greens_func_complex(n_basis,n_basis,n_k_points_task,n_spin))

        write(use_unit,*) "Testing calculate_imaginary_greens_func_comp_eigenvecs"
        CALL calculate_imaginary_greens_func_comp_eigenvecs(tau,KS_eigenvalue, &
                                                            KS_eigenvector_complex, &
                                                            occ_numbers, &
                                                            greens_func_complex)
        write(use_unit,*) "Test calculate_imaginary_greens_func_comp_eigenvecs: ok"

        DEALLOCATE(greens_func_complex)

    END SUBROUTINE calculate_imaginary_greens_func_comp_eigenvecs_test

    SUBROUTINE fouriertransform_comp_greensfunction_k_to_R_test()
      REAL*8 :: tau
      REAL*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      COMPLEX*16, dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector_complex
      REAL*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
      COMPLEX*16, allocatable :: greens_func_complex(:,:,:,:)
      INTEGER :: rest
      type(realspace_vectors) :: rvecs
      type(basisfunc_distribution), DIMENSION(n_tasks) :: dist
      type(greens_functions) :: a_greens_func

       KS_eigenvalue = (1,1)
       KS_eigenvector_complex = 2.0d0
       occ_numbers = 1
       tau = 0.0d0

       ALLOCATE(greens_func_complex(n_basis,n_basis,n_k_points_task,n_spin))

       write(use_unit,*) "Testing fouriertransform_comp_greensfunction_k_to_R_test"

       CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk,lattice_vector, rvecs)
       CALL distribute_basisfunctions(n_tasks,n_basis,rvecs, dist)

       CALL init_green_func_storage(a_greens_func,rvecs,1)

       CALL fouriertransform_comp_greensfunction_k_to_R(greens_func_complex, &
                                                        rvecs, &
                                                        a_greens_func)

       CALL free_all_lattice_vectors(rvecs)
       DEALLOCATE(greens_func_complex)
       CALL finalize_green_func_storage(a_greens_func)

       write(use_unit,*) "Test fouriertransform_comp_greensfunction_k_to_R_test: ok"

    END SUBROUTINE fouriertransform_comp_greensfunction_k_to_R_test

    SUBROUTINE get_n_good_mean_energies_test()
        INTEGER :: n_energies, n, i
        REAL*8, ALLOCATABLE,DIMENSION(:) :: energies, mean_energies


        write(use_unit,*) "Testing get_n_good_mean_energies"


        n_energies = 10

        ALLOCATE(energies(n_energies))

        do i=1, n_energies
            energies(i) = -dble(i)/2.d0
        enddo


        do n = 1, n_energies, 3

            ALLOCATE(mean_energies(n))

            CALL get_n_good_mean_energies(n_energies,energies, &
                                          n,mean_energies)

            write(use_unit,*) mean_energies

            DEALLOCATE(mean_energies)
        enddo
        DEALLOCATE(energies)


        write(use_unit,*) "Passed testing get_n_good_mean_energies"


    END SUBROUTINE get_n_good_mean_energies_test

    SUBROUTINE exp_basis_test()
        type(exp_basis) :: an_exp_basis

        REAL*8, ALLOCATABLE, DIMENSION(:) :: exponents,e_x, x_n, y_n
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: A
        REAL*8 :: x, max_diff,tmp
        INTEGER :: n_exponents,n,i

        write(use_unit,*) "Testing exp_basis"


        do n_exponents = 1, 5
            ALLOCATE(exponents(n_exponents))

            do i = 1, n_exponents
                exponents(i) = -dble(i)/2.d0
            enddo

            CALL create_exp_basis(n_exponents, exponents, an_exp_basis)

            CALL print_exp_basis(an_exp_basis)
!
!            x=1.d0
!            ALLOCATE(e_x(n_exponents))
!            CALL evaluate_exp_basis(x,an_exp_basis, e_x)
!            write(use_unit,*) "x",x,"e_x",e_x
!            DEALLOCATE(e_x)

            n = n_exponents
            ALLOCATE(x_n(n))
            ALLOCATE(y_n(n))
            ALLOCATE(A(n,n))

            do i = 1,n
                x_n(i) = dble(i)/dble(n)
            enddo

            CALL get_exp_interpolation_matrix(n,x_n,an_exp_basis,A)

            do i=1,n
               y_n(i) = SUM(exp(exponents(:)*x_n(i)))
            enddo

            ALLOCATE(e_x(1))
            CALL prepare_fouriertransform_by_exp_basis(n,x_n, an_exp_basis,x_n)
            CALL fouriertransform_by_exp_basis(n,x_n,1,y_n,10.d0,an_exp_basis,e_x)
            write(use_unit,*) "f_t", e_x
            DEALLOCATE(e_x)


            x_n = MATMUL(A,y_n)
            max_diff = maxval(abs(x_n-1))

            write(use_unit,*) "max_diff rep of y_n :", n,max_diff

            if(max_diff > 1.d-4) then
                write(use_unit,*) x_n
  !              stop "FAILED to expand exponentials"
            endif


            DEALLOCATE(y_n)
            DEALLOCATE(x_n)
            DEALLOCATE(A)


            CALL free_exp_basis(an_exp_basis)
            DEALLOCATE(exponents)
        enddo

        write(use_unit,*) "Passed testing exp_basis"

    END SUBROUTINE exp_basis_test

    SUBROUTINE calc_lvl_tricoeff_test()
        DOUBLE COMPLEX :: lvl_tricoeff_k(n_basis,n_basis,max_n_basbas_sp)
        INTEGER :: i_k_point, i_basbas

        CALL write_stdout("Testing calc_lvl_tricoeff")


if(n_tasks == 1) then
        do i_k_point = 1,n_k_points
            do i_basbas = 1,max_n_basbas_sp !,max_n_basbas_sp
                write(use_unit,*) "Testing basbas", i_basbas
!                CALL fouriertransform_exp_coeff(exp_coeffs_local,i_basbas,1)
!                write(use_unit,*) "Max diff",i_basbas,maxval(abs(lvl_tricoeff_recip1(i_basbas,:,:,i_k_point)-&
!                                                        exp_coeffs_local%ft_coeff(:,:,i_k_point)))
            enddo
        enddo
endif

        CALL write_stdout("Passed test calc_lvl_tricoeff")
    END SUBROUTINE calc_lvl_tricoeff_test
END MODULE cRPA_calculation_test
