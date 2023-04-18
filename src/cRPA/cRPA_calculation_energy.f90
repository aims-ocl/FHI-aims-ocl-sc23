!Calculates the RPA correlation energy
MODULE cRPA_calculation_energy
      use dimensions
      use prodbas
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      USE cRPA_calculation_integration
      USE cRPA_calculation_fouriertrans
      use evaluate_polarisability_kspace_mod
      USE gw_para

      implicit none

      ABSTRACT INTERFACE
           SUBROUTINE get_resp_func_for_omega_template(n_basbas, i_k_point, &
                                                      frequency, &
                                                      data_out)
              INTEGER, INTENT(IN) :: n_basbas,i_k_point
              REAL*8, INTENT(IN) :: frequency
              DOUBLE COMPLEX, DIMENSION(n_basbas, n_basbas), INTENT(OUT) :: data_out
            END SUBROUTINE get_resp_func_for_omega_template
      END INTERFACE

      INTEGER :: n_frequencies

      REAL*8, ALLOCATABLE, DIMENSION(:)  :: omega_points
      REAL*8, ALLOCATABLE, DIMENSION(:)  :: omega_weights

      complex*16, dimension(:,:,:,:), allocatable :: polar_Rspace

      INTEGER :: LOOP_ORDER_K_FREQ
      PARAMETER(LOOP_ORDER_K_FREQ = 1)

      INTEGER :: LOOP_ORDER_FREQ_K
      PARAMETER(LOOP_ORDER_FREQ_K = 2)

CONTAINS
    SUBROUTINE set_energy_omega_grid(grid_omega)
        type(integration_grid), INTENT(IN) :: grid_omega
        n_frequencies = grid_omega%n_points
!TODO check alloc
        ALLOCATE(omega_points(n_frequencies))
        ALLOCATE(omega_weights(n_frequencies))

        omega_points(:) = grid_omega%points(:)
        omega_weights(:) = grid_omega%weights(:)

!        CALL print_integration_grid(grid_omega)

    END SUBROUTINE set_energy_omega_grid

    SUBROUTINE unset_energy_omega_grid()
        DEALLOCATE(omega_points)
        DEALLOCATE(omega_weights)
    END SUBROUTINE unset_energy_omega_grid

    DOUBLE COMPLEX FUNCTION matrix_trace(n_dim,mat)
        INTEGER, INTENT(IN) :: n_dim
        DOUBLE COMPLEX,DIMENSION(n_dim,n_dim), INTENT(IN) :: mat
        INTEGER :: i

        matrix_trace = DCMPLX(0.d0,0.d0)
        do i=1, n_dim
           matrix_trace=matrix_trace + mat(i,i)
        enddo
    END FUNCTION matrix_trace

    DOUBLE COMPLEX FUNCTION matrix_determinante_by_LU(n_dim,mat)
        INTEGER, INTENT(IN) :: n_dim
        DOUBLE COMPLEX,DIMENSION(n_dim,n_dim), INTENT(INOUT) :: mat
        INTEGER :: i,ipiv(n_dim),info

        CALL  zgetrf (n_dim,n_dim,mat,n_basbas,ipiv,info)

        matrix_determinante_by_LU = DCMPLX(1.d0,0.d0)
        do i = 1, n_dim
             matrix_determinante_by_LU = matrix_determinante_by_LU &
                                         * mat(i,i)
        enddo

    END FUNCTION matrix_determinante_by_LU

    SUBROUTINE build_in_place_identity_minus_matrix(n_dim,mat)
        INTEGER, INTENT(IN) :: n_dim
        DOUBLE COMPLEX,DIMENSION(n_dim,n_dim), INTENT(INOUT) :: mat
        INTEGER :: i,j

        do j=1,n_dim
            do i=1,n_dim
                if(i==j) then
                    mat(i,j) = DCMPLX(1.d0,0.d0) -mat(i,j)
                else
                    mat(i,j) = -mat(i,j)
                endif
            enddo
        enddo
    END SUBROUTINE build_in_place_identity_minus_matrix

!Copied and modified from evaluate_crpa_energy_kspace.f90
    SUBROUTINE get_crpa_energy &
              (get_resp_func_for_omega,i_loop_order, rpa_c_energy_out)

      PROCEDURE(get_resp_func_for_omega_template),POINTER,INTENT(INOUT) :: &
                                                         get_resp_func_for_omega
      INTEGER, INTENT(IN) :: i_loop_order

      REAL*8, INTENT(OUT) :: rpa_c_energy_out

      DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: polar_kspace, &
                                                     aux_mat

      DOUBLE COMPLEX :: det_v_times_polar, &
                        trace_v_times_polar
      REAL*8  rpa_c_integrand, &
              rpa_c_energy, &
              rpa_c_contribution(n_frequencies), &
              temp_time_rpa, temp_clock_time_rpa, &
              temp_time_resp, temp_clock_time_resp, &
              time_in_resp, time_clock_in_resp

      INTEGER :: i_loop1, i_loop1_end, &
                 i_loop2, i_loop2_end, &
                 i_freq, &
                 i_k_point, &
                 i_k_point_local,&
                 allocation_info

      CALL get_timestamps(temp_time_rpa, temp_clock_time_rpa )

      ALLOCATE(polar_kspace(n_basbas, n_basbas),stat=allocation_info)
      CALL check_allocation(allocation_info, 'polar_kspace                    ')

      time_in_resp = 0.d0
      time_clock_in_resp = 0.d0

      rpa_c_energy = 0.d0
      rpa_c_contribution(:) = 0.d0
      rpa_c_integrand = 0.d0

      if(i_loop_order == LOOP_ORDER_K_FREQ) then
         i_loop1_end = n_k_points
         i_loop2_end = n_frequencies
      elseif(i_loop_order == LOOP_ORDER_FREQ_K) then
         i_loop1_end = n_frequencies
         i_loop2_end = n_k_points
      endif


      do i_loop1 = 1, i_loop1_end
         do i_loop2 = 1, i_loop2_end

            if(i_loop_order == LOOP_ORDER_K_FREQ) then
              i_k_point = i_loop1
              i_freq = i_loop2
            elseif(i_loop_order == LOOP_ORDER_FREQ_K) then
              i_freq = i_loop1
              i_k_point = i_loop2
            endif

            CALL get_timestamps(temp_time_resp, temp_clock_time_resp)

            CALL get_resp_func_for_omega(n_basbas,i_k_point,omega_points(i_freq), &
                                         polar_kspace)

            CALL get_times(temp_time_resp, temp_clock_time_resp)
            time_in_resp = time_in_resp + temp_time_resp
            time_clock_in_resp = time_clock_in_resp + temp_clock_time_resp


            if(myid .ne. mod(i_k_point, n_tasks)) cycle
            i_k_point_local = (i_k_point-1)/n_tasks + 1


            ALLOCATE(aux_mat(n_basbas, n_basbas),stat=allocation_info)
            CALL check_allocation(allocation_info, 'aux_mat')


            CALL zgemm('T', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                  coulomb_matr_recip(1,1,i_k_point_local), n_basbas, &
                  polar_kspace(1,1), n_basbas, (0.d0, 0.d0), &
                  aux_mat, n_basbas)

            trace_v_times_polar = matrix_trace(n_basbas,aux_mat)

            CALL build_in_place_identity_minus_matrix(n_basbas,aux_mat)
            det_v_times_polar = matrix_determinante_by_LU(n_basbas,aux_mat)

            DEALLOCATE(aux_mat)


            rpa_c_integrand = log (abs(real(det_v_times_polar))) + &
                                   real(trace_v_times_polar)

            rpa_c_contribution(i_freq) =  rpa_c_contribution(i_freq) + &
                                          rpa_c_integrand * omega_weights(i_freq) * &
                                          k_weights(i_k_point) /(2.d0*pi)

            rpa_c_energy = rpa_c_energy + &
                           rpa_c_integrand * omega_weights(i_freq) * &
                           k_weights(i_k_point)

          enddo
      enddo

      CALL sync_vector(rpa_c_contribution, n_frequencies, mpi_comm_global)

      if(myid == 0) then
          do i_freq=1, n_frequencies
            write(use_unit,*) "Contribution for",omega_points(i_freq), &
                                          rpa_c_contribution(i_freq)
          enddo
      endif

      call sync_real_number(rpa_c_energy)

      rpa_c_energy=rpa_c_energy/(2.d0*pi)


      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            " rsRPA correlation energy :", rpa_c_energy, "Ha,", &
             rpa_c_energy*hartree, "eV"

        write(use_unit,*)"----------------------------------------------------", &
                 "-------------------------"
        write(use_unit,*)
      endif

      if (allocated (polar_kspace)) then
        DEALLOCATE (polar_kspace)
      endif

      CALL get_times(temp_time_rpa, temp_clock_time_rpa)


      CALL output_times("2X","Overall time in get_resp_func_for_omega",time_in_resp, time_clock_in_resp)
      CALL output_times("2X","Overall time in get_crpa_energy",temp_time_rpa, temp_clock_time_rpa)

      rpa_c_energy_out = rpa_c_energy

END SUBROUTINE get_crpa_energy

!used only for debugging
!calculate chi_0 in k for later transform to realspace
subroutine calculate_kspace_polarisability &
           ( n_low_state, occ_numbers, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, rvecs)

      integer :: n_low_state
      real*8  :: occ_numbers(n_states,n_spin)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)
      type(realspace_vectors), INTENT(IN) :: rvecs



      integer :: ipiv(n_basbas)
      integer :: info

      complex*16, dimension(:,:,:), allocatable :: polar_kspace

      integer :: i_freq
      integer :: i_index
      INTEGER :: i_column

      allocate(polar_kspace(n_basbas, n_basbas, n_k_points),stat=i_index)
      call check_allocation(i_index, 'polar_kspace                    ')

      allocate(polar_Rspace(n_basbas, n_basbas, n_k_points, n_frequencies),stat=i_index)
      call check_allocation(i_index, 'polar_Rspace                    ')


! work array

      do i_freq = 1, n_frequencies
         call write_stdout("Frequency "//num2str(i_freq)//num2str(n_frequencies))
         call  evaluate_polarisability_kspace &
              ( 1, omega_points(i_freq), &
! from Xinguo
!              (  n_low_state, omega_points(i_freq), &

                 KS_eigenvalue, KS_eigenvector, &
                 KS_eigenvector_complex, &
                 occ_numbers, polar_kspace )

          do i_column = 1, n_basbas
              CALL discrete_fouriertransform_k_to_R(rvecs, n_basbas, &
                                                    polar_kspace(:,i_column,:), &
                                                    polar_Rspace(:,i_column,:,i_freq))
          enddo

      enddo


      DEALLOCATE(polar_kspace)
END SUBROUTINE

SUBROUTINE get_kspace_polarisability_for_tau(tau, rvecs, polar_R_for_tau)
    REAL*8,INTENT(IN) :: tau
    type(realspace_vectors), INTENT(IN) :: rvecs
    COMPLEX*16,DIMENSION(n_basbas,n_basbas, rvecs%n_vecs), INTENT(OUT) :: polar_R_for_tau

    INTEGER :: i_column, i_vec
    type(integration_grid) :: int_grid

    CALL create_integration_grid(n_frequencies, omega_points, omega_weights, &
                                 int_grid)

    do i_vec = 1, rvecs%n_vecs
        do i_column = 1, n_basbas
            CALL fouriertransform_t_to_x(int_grid,n_basbas, &
                                         polar_Rspace(:,i_column,i_vec,:), &
                                         tau, &
                                         polar_R_for_tau(:,i_column,i_vec))
        enddo
    enddo

    polar_R_for_tau(:,:,:) = 2.d0 * polar_R_for_tau

    CALL free_integration_grid(int_grid)
END SUBROUTINE get_kspace_polarisability_for_tau

END MODULE cRPA_calculation_energy
