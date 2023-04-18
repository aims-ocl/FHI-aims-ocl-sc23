!Calculates the greensfunction g_0 of the noninteracting reference system
MODULE cRPA_calculation_greensfunc
    USE cRPA_view
    USE cRPA_storage
    USE cRPA_parallelism

    IMPLICIT NONE

CONTAINS

!calculate green's funtion
!notice here that tau>0.d0 calculates the advaced and tau < 0.d0 calculates the
!retarded greens function
    SUBROUTINE calculate_local_greensfunc(rvecs, cur_tau)
        type(realspace_vectors), INTENT(IN) :: rvecs
        REAL*8, INTENT(IN) :: cur_tau

        CALL calculate_greensfunc(rvecs, cur_tau, &
                                  green_funcs_local)



!      CALL debug_print_green_func(green_funcs_local)
!      stop


        CALL calculate_greensfunc(rvecs, -1.d0*cur_tau, &
                                  green_funcs_local_minus_tau)

    END SUBROUTINE calculate_local_greensfunc

    SUBROUTINE calculate_greensfunc(rvecs,cur_tau, &
                                    greensfunc)
        type(realspace_vectors), INTENT(IN) :: rvecs
        REAL*8, INTENT(IN) :: cur_tau
        type(greens_functions), INTENT(INOUT) :: greensfunc

        REAL*8, allocatable :: greens_func_real(:,:,:,:)
        DOUBLE COMPLEX, allocatable :: greens_func_complex(:,:,:,:)

        if(real_eigenvectors) then
            stop
!            ALLOCATE(greens_func_real(n_basis,n_basis,n_k_points_task,n_spin))
!
!            CALL calculate_imaginary_greens_func_real_eigenvecs(cur_tau,KS_eigenvalue, &
!                                                                KS_eigenvector, &
!                                                                occ_numbers, &
!                                                                greens_func_real)
!
!            CALL fouriertransform_real_greensfunction_k_to_R(greens_func_real, &
!                                                             rvecs, &
!                                                             my_basis_dist, &
!                                                             real_or_img_part, &
!                                                             greensfunc)

            DEALLOCATE(greens_func_real)
        else
            ALLOCATE(greens_func_complex(n_basis,n_basis,n_k_points_task,n_spin))

            CALL calculate_imaginary_greens_func_comp_eigenvecs(cur_tau,KS_eigenvalue, &
                                                                KS_eigenvector_complex, &
                                                                occ_numbers, &
                                                                greens_func_complex)

            CALL fouriertransform_comp_greensfunction_k_to_R(greens_func_complex, &
                                                             rvecs, &
                                                             greensfunc)

            DEALLOCATE(greens_func_complex)
        endif

        CALL set_green_func_norms(greensfunc)

!        CALL print_green_func(greensfunc)

    END SUBROUTINE calculate_greensfunc

!TODO refactor
    SUBROUTINE score_greensfunc(tau, norm)
        REAL*8, INTENT(IN) :: tau
        REAL*8, INTENT(OUT) :: norm

        REAL*8 :: energy_fermi(2)
        INTEGER, DIMENSION(n_spin) :: i_k_point_of_highest_state
        INTEGER, DIMENSION(n_k_points,n_spin) :: n_state_highest_occupied
        REAL*8 :: decay_constant, cur_decay_constant

        INTEGER :: i_k, i_spin, i_state, n_state_start, n_state_end

        CALL write_stdout("Scoring greens function for tau="//num2str(tau))


        CALL get_highest_occupied_state(occ_numbers, n_state_highest_occupied, &
                                        i_k_point_of_highest_state)


        CALL get_energy_fermi(occ_numbers,KS_eigenvalue,energy_fermi)

        cur_decay_constant = 1.d6

        norm = 0.d0
        do i_spin = 1, n_spin
            do i_k = 1,n_k_points

                if(tau<=0.d0) then
                    n_state_start = 1
                    n_state_end = n_state_highest_occupied(i_k,i_spin)
                else
                    n_state_start = n_state_highest_occupied(i_k,i_spin) +1
                    n_state_end = n_states
                endif

                do i_state = n_state_start, n_state_end
                    cur_decay_constant = KS_eigenvalue(i_state,i_spin,i_k) &
                                         - energy_fermi(i_spin)

                    if(cur_decay_constant > 0.d0) then
                        norm = norm + exp(-tau*cur_decay_constant)
                        decay_constant = min(decay_constant, &
                                             cur_decay_constant)
                    endif
                enddo
            enddo
        enddo


!        write(use_unit,*) "Scored tau ",  tau, "with ", norm!, "using ", decay_constant

    END SUBROUTINE score_greensfunc

    SUBROUTINE calculate_imaginary_greens_func_real_eigenvecs(tau,KS_eigenvalue, &
                                                              KS_eigenvector, occ_numbers, &
                                                              greens_func_real)
        REAL*8, INTENT(IN) :: tau
        REAL*8, INTENT(IN) :: KS_eigenvalue(n_states,n_spin,n_k_points)
        REAL*8, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN):: KS_eigenvector
        REAL*8, dimension(n_states, n_spin, n_k_points), INTENT(IN) :: occ_numbers
        REAL*8, allocatable :: greens_func_real(:,:,:,:)


        REAL*8, allocatable :: tmp_r(:,:)
        REAL*8 :: cur_occ_number,cur_energy_diff,cur_exp_value, max_occ_number


        INTEGER :: i_state, i_spin, i_k, i_k_point

        ALLOCATE(tmp_r(n_basis,n_states))
!TODO check allocation

        max_occ_number = 3.d0 - dble(n_spin)

        i_k = 0
        do i_k_point = 1, n_k_points
          if(myid == MOD(i_k_point, n_tasks)) then
            i_k = i_k + 1
            do i_spin = 1, n_spin
                tmp_r(:,:) = 0.d0

                do i_state = 1, n_states

                    if(tau<=0.d0) then
                       cur_occ_number = occ_numbers(i_state,i_spin,i_k_point)
                    else
                       cur_occ_number = max_occ_number-occ_numbers(i_state,i_spin,i_k_point)
                    endif

                    cur_energy_diff = (KS_eigenvalue(i_state,i_spin,i_k_point))! - energy_fermi(i_spin))
                    cur_exp_value = exp(-1.d0 * tau * cur_energy_diff)

                    if(cur_occ_number > 1.d-5) then

                       tmp_r(:,i_state) = KS_eigenvector(:,i_state,i_spin,i_k) &
                                          * cur_occ_number &
                                          * cur_exp_value
                    endif
                enddo

                call dgemm('N','T',n_basis,n_basis,n_states,1.d0,KS_eigenvector(1,1,i_spin,i_k),ubound(KS_eigenvector,1), &
                           tmp_r,ubound(tmp_r,1),0.d0,greens_func_real(1,1,i_k,i_spin),ubound(greens_func_real,1))
            enddo
          endif
        enddo

    deallocate(tmp_r)

    END SUBROUTINE calculate_imaginary_greens_func_real_eigenvecs

    SUBROUTINE calculate_imaginary_greens_func_comp_eigenvecs(tau,KS_eigenvalue, &
                                                              KS_eigenvector_complex, occ_numbers, &
                                                              greens_func_complex)
        REAL*8, INTENT(IN) :: tau
        REAL*8, INTENT(IN)  :: KS_eigenvalue(n_states,n_spin,n_k_points)
        COMPLEX*16, INTENT(IN), dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector_complex
        REAL*8, INTENT(IN), dimension(n_states, n_spin, n_k_points) :: occ_numbers
        COMPLEX*16, allocatable :: greens_func_complex(:,:,:,:)

        COMPLEX*16, allocatable :: tmp_c(:,:)


        INTEGER :: i_state, i_spin, i_k, i_k_point, &
                   n_state_start, n_state_end,i
        REAL*8,DIMENSION(n_spin) :: energy_fermi

        INTEGER, DIMENSION(n_spin) :: i_k_point_of_highest_state
        INTEGER, DIMENSION(n_k_points,n_spin) :: n_state_highest_occupied
        REAL*8 :: cur_occ_number,cur_energy_diff,cur_exp_value, max_occ_number

!        CALL write_stdout("Calculating greens function for tau="//num2str(tau))

        ALLOCATE(tmp_c(n_basis,n_states))
!TODO check allocation

        CALL get_highest_occupied_state(occ_numbers, n_state_highest_occupied, &
                                        i_k_point_of_highest_state)
!        CALL write_stdout("k-point of highest state:"//num2str(i_k_point_of_highest_state(1)))

        CALL get_energy_fermi(occ_numbers,KS_eigenvalue,energy_fermi)
!        CALL write_stdout("Fermi KS _energy_:"//num2str(energy_fermi(1))//num2str(energy_fermi(n_spin)))


        max_occ_number = 3.d0 - dble(n_spin)

        i_k = 0
        do i_k_point = 1, n_k_points
!write(use_unit,*) i_k_point,k_phase_exx(i_k_point,:)

          if(myid == MOD(i_k_point, n_tasks)) then
            i_k = i_k + 1
            do i_spin = 1, n_spin
                tmp_c(:,:) = 0.d0

                if(tau<=0.d0) then
                    n_state_start = 1
                    n_state_end = n_states ! n_state_highest_occupied(i_k_point,i_spin)
                else
                    n_state_start = 1 !n_state_highest_occupied(i_k_point,i_spin)+1
                    n_state_end = n_states
                endif
!CALL write_debug("n_state_start:"//num2str(n_state_start))
!CALL write_debug("n_state_end:"//num2str(n_state_end))

                do i_state = n_state_start, n_state_end

!write(use_unit,*) "kS_eigenstate",i_k_point,i_state,i_spin, KS_eigenvector_complex(1:2,i_state,i_spin,i_k)

                    if(tau<=0.d0) then
                        cur_occ_number = occ_numbers(i_state,i_spin,i_k_point)
!if(myid==1) write(use_unit,*) "occu", i_state,i_spin,i_k_point,cur_occ_number
                    else
                       cur_occ_number = max_occ_number-occ_numbers(i_state,i_spin,i_k_point)
!if(myid==1) write(use_unit,*) "unoccu", i_state,i_spin,i_k_point,cur_occ_number
                    endif

                    cur_energy_diff = (KS_eigenvalue(i_state,i_spin,i_k_point) - energy_fermi(i_spin))
                    cur_exp_value = exp(-1.d0 * tau * cur_energy_diff)

                    if(cur_occ_number > 1.d-5) then

                       tmp_c(:,i_state) = KS_eigenvector_complex(:,i_state,i_spin,i_k) &
                                          * DCMPLX(cur_occ_number * &
                                                  cur_exp_value,   &
                                                  0.d0)
                    endif
!write(use_unit,*) "cur_energy_diff", cur_energy_diff, exp(-1.d0 * tau * cur_energy_diff)
!write(use_unit,*) "kS_eigenstate",i_k_point,i_state,i_spin, KS_eigenvector_complex(:,i_state,i_spin,i_k)
!write(use_unit,*) "cur_occ ",cur_occ_number
!write(use_unit,*) "KS_eigenvalue", i_spin,i_k_point,i_state, KS_eigenvalue(i_state,i_spin,i_k_point)
!write(use_unit,*) "tmp_c", tmp_c(:,i_state)
                enddo

                call zgemm('N','C', &
                           n_basis,n_basis,n_states, &
                           DCMPLX(1.d0,0.d0), &
                           KS_eigenvector_complex(:,:,i_spin,i_k), &
                           ubound(KS_eigenvector_complex,1), &
                           tmp_c,ubound(tmp_c,1), &
                           DCMPLX(0.d0,0.d0),&
                           greens_func_complex(:,:,i_k,i_spin),ubound(greens_func_complex,1))


!do i = 1, n_basis
!write(use_unit,*) "greens_func_complex",tau,i_state,dble(greens_func_complex(:,i,i_k,i_spin))
!enddo
!write(use_unit,*) "green", maxval(abs(&
!            dble(greens_func_complex(:,:,i_k,i_spin)- &
!            transpose(dble(greens_func_complex(:,:,i_k,i_spin))))))
            enddo
          endif
        enddo

    deallocate(tmp_c)

    END SUBROUTINE calculate_imaginary_greens_func_comp_eigenvecs

    SUBROUTINE get_highest_occupied_state(occ_numbers, &
                                          n_state_highest_occupied, &
                                          k_point_of_highest_state)
        REAL*8, INTENT(IN), dimension(n_states, n_spin, n_k_points) :: occ_numbers
        INTEGER, DIMENSION(n_k_points,n_spin), INTENT(OUT) :: n_state_highest_occupied
        INTEGER, DIMENSION(n_spin), INTENT(OUT) :: k_point_of_highest_state

        INTEGER, DIMENSION(n_spin) :: i_highest_state
        INTEGER :: i_spin, i_k_point, i_state

        i_highest_state(:) = 0
        k_point_of_highest_state(:) = 0
        n_state_highest_occupied(:,:) = 0
        do i_spin = 1, n_spin
            do i_k_point = 1, n_k_points
                do i_state = 1, n_states
                    if(occ_numbers(i_state,i_spin,i_k_point) .gt. 1.e-12) then
                        n_state_highest_occupied(i_k_point,i_spin) = i_state

                        if(i_state > i_highest_state(i_spin)) then
                            i_highest_state(i_spin) = i_state
                            k_point_of_highest_state(i_spin) = i_k_point
                        endif
                    endif
                enddo
 !     if(myid.eq.0) then
 !       write(use_unit,*) i_spin, i_k_point, n_homo_k(i_k_point, i_spin), n_homo(i_spin)
 !     endif
            enddo
        enddo

        !n_homo(1)=maxval(n_homo_k(:,1), n_k_points)
        !n_homo(n_spin)=maxval(n_homo_k(:,n_spin), n_k_points)
        !max_n_homo = max(n_homo(1), n_homo(n_spin))


    END SUBROUTINE get_highest_occupied_state

    SUBROUTINE get_energy_fermi(occ_numbers,KS_eigenvalue,energy_fermi)
        REAL*8, INTENT(IN)  :: KS_eigenvalue(n_states,n_spin,n_k_points)
        REAL*8, INTENT(IN), dimension(n_states, n_spin, n_k_points) :: occ_numbers
        REAL*8,DIMENSION(n_spin), INTENT(OUT) :: energy_fermi


        INTEGER, DIMENSION(n_k_points,n_spin) :: n_state_highest_occupied
        INTEGER, DIMENSION(n_spin) :: i_k_point_of_highest_state
        INTEGER :: i_spin, i_k_point, i_state
        REAL*8 :: cur_max_energy

        CALL get_highest_occupied_state(occ_numbers,n_state_highest_occupied, &
                                        i_k_point_of_highest_state)

        do i_spin = 1, n_spin
            energy_fermi(i_spin)=minval(KS_eigenvalue(:,i_spin,:))
            do i_k_point = 1,n_k_points
                do i_state=1, n_states
        !            cur_max_energy = maxval(KS_eigenvalue(1:n_state_highest_occupied(i_k_point,i_spin),&
         !                                                           i_spin,i_k_point))

                     if(occ_numbers(i_state,i_spin,i_k_point) >= 0.99d0) then
                     cur_max_energy = KS_eigenvalue(i_state, i_spin,i_k_point)
                energy_fermi(i_spin) = max(energy_fermi(i_spin), &
                                           cur_max_energy)
                     endif
                enddo
            enddo
        enddo

!        energy_fermi(:) = maxval(energy_fermi(:))
        energy_fermi(:) = chemical_potential_spin(:)
    END SUBROUTINE get_energy_fermi


    SUBROUTINE get_possible_energy_differences(occ_numbers,KS_eigenvalue, &
                                               min_energy_diff, &
                                               n_energy_diff,energy_diff)
        REAL*8, INTENT(IN)  :: KS_eigenvalue(n_states,n_spin,n_k_points)
        REAL*8, INTENT(IN), dimension(n_states, n_spin, n_k_points) :: occ_numbers
        REAL*8, INTENT(IN) :: min_energy_diff
        INTEGER, INTENT(OUT) :: n_energy_diff
        REAL*8,DIMENSION(:),ALLOCATABLE, INTENT(OUT) :: energy_diff


        INTEGER, DIMENSION(n_k_points,n_spin) :: n_state_highest_occupied
        INTEGER, DIMENSION(n_spin) :: i_k_point_of_highest_state
        INTEGER :: i_spin, i_k_point, i_q_point, i_state, i_state_snd
        REAL*8 :: tmp_swap,tmp_candidate

        REAL*8, DIMENSION(size(KS_eigenvalue)**2) :: tmp_energy_diffs
        INTEGER :: i_pos,i_sort

        tmp_energy_diffs(:) = 0.d0
        i_pos = 1

        do i_spin = 1, n_spin
            do i_k_point = 1,n_k_points
                do i_q_point = 1,n_k_points
                    do i_state=1, n_states
                       do i_state_snd = 1,n_states
                         if(abs(occ_numbers(i_state,i_spin,i_k_point) * &
                            (1.d0-occ_numbers(i_state_snd,i_spin,i_q_point))-1.d0) < 1.d-1) then

                            tmp_candidate = KS_eigenvalue(i_state,i_spin, i_k_point) &
                                            -&
                                            KS_eigenvalue(i_state_snd,i_spin, i_q_point)
                            if(tmp_candidate>=0.d0) cycle

                            if(i_pos == 1) then
                                tmp_energy_diffs(1) = tmp_candidate
                                i_pos = i_pos + 1
                            elseif(.NOT. has_energy(tmp_candidate)) then
                                tmp_energy_diffs(i_pos) = tmp_candidate
                                i_pos = i_pos + 1

!                                do i_sort=i_pos-1, i_pos
!                                    if(tmp_energy_diffs(i_sort)>tmp_energy_diffs(i_sort-1)) then
!                                       tmp_swap = tmp_energy_diffs(i_sort)
!                                       tmp_energy_diffs(i_sort) = tmp_energy_diffs(i_sort-1)
!                                       tmp_energy_diffs(i_sort-1) = tmp_swap
!                                    endif
!                                enddo

                            endif
                         endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

        n_energy_diff = i_pos - 1
        ALLOCATE(energy_diff(n_energy_diff))
        energy_diff(:) = tmp_energy_diffs(1:n_energy_diff)

write(use_unit,*) "energy_diff",n_energy_diff, energy_diff(1:10)

    CONTAINS
        LOGICAL FUNCTION has_energy(test_energy)
            REAL*8, INTENT(IN) :: test_energy

            INTEGER :: i

            do i=1,i_pos-1
                if(abs(tmp_energy_diffs(i)-test_energy) < 1.d-10) then
                    has_energy = .TRUE.
                    return
                endif
            enddo

            has_energy = .FALSE.
        END FUNCTION
    END SUBROUTINE get_possible_energy_differences


    SUBROUTINE fouriertransform_real_greensfunction_k_to_R (greens_func_real, &
                                                            rvecs, &
                                                            distribution, &
                                                            real_or_img_part, &
                                                            greensfunc)
      REAL*8, allocatable :: greens_func_real(:,:,:,:)
      type(realspace_vectors) :: rvecs
      type(basisfunc_distribution) :: distribution
      INTEGER,INTENT(IN) :: real_or_img_part
      type(greens_functions), INTENT(INOUT) :: greensfunc


      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dm_tmp
      INTEGER :: i_cell,i_spin, i_k,i_k_point

      ALLOCATE(dm_tmp(n_basis,n_basis,n_spin))
!TODO check alloc

      do i_cell = 1, n_cells_bvk
        dm_tmp(:,:,:) = 0.d0

        do i_spin = 1, n_spin
          i_k = 0
          do i_k_point = 1, n_k_points
            if(myid == MOD(i_k_point, n_tasks)) then
              i_k = i_k + 1
              dm_tmp(:,:,i_spin) = dm_tmp(:,:,i_spin) + &
                           greens_func_real(:,:,i_k,i_spin) *  &
                           dble (k_phase_exx(i_cell,i_k_point)*k_weights(i_k_point))
            endif
          enddo
        enddo

        call sync_vector(dm_tmp, size(dm_tmp))


        CALL store_green_func(greensfunc, rvecs%vecs(i_cell), &
                              n_basis, &
                              n_basis, &
                              dm_tmp(:,:,:)*DCMPLX(1.d0,0.d0))

      enddo

      DEALLOCATE(dm_tmp)

    END SUBROUTINE fouriertransform_real_greensfunction_k_to_R

    SUBROUTINE fouriertransform_comp_greensfunction_k_to_R (greens_func_complex, &
                                                            rvecs, &
                                                            greensfunc)

      DOUBLE COMPLEX, ALLOCATABLE,INTENT(IN) :: greens_func_complex(:,:,:,:)
      type(realspace_vectors), INTENT(IN) :: rvecs
      type(greens_functions), INTENT(INOUT) :: greensfunc


      DOUBLE COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: dm_tmp
      INTEGER :: i,j,k,i_cell,i_spin, i_k,i_k_point

!      CALL write_stdout("Calculate fouriertransformation of greensfunction (k -> R)")

      ALLOCATE(dm_tmp(n_basis,n_basis,n_spin))
!TODO check alloc
      do i_cell = 1, n_cells_bvk
        dm_tmp(:,:,:) = 0.d0

        do i_spin = 1, n_spin
          i_k = 0
          do i_k_point = 1, n_k_points
            if(myid == MOD(i_k_point, n_tasks)) then
              i_k = i_k + 1
                dm_tmp(:,:,i_spin) = dm_tmp(:,:,i_spin) + &
                              greens_func_complex(:,:,i_k,i_spin) * &
                              conjg(k_phase_exx(i_cell,i_k_point))*k_weights(i_k_point)
            endif
          enddo
        enddo

        dm_tmp(:,:,:) = dm_tmp(:,:,:) !* k_weights(1)

        call sync_vector_complex(dm_tmp, size(dm_tmp))

      !   CALL prune_nan()

!CALL write_debug("Storing greens func "//num2str(i_cell) &
!                //num2str(rvecs(i_cell)%i_cell1)&
!                //num2str(rvecs(i_cell)%i_cell2)&
!                //num2str(rvecs(i_cell)%i_cell3))
!CALL write_debug("gf img real:"//num2str(chksum_2d_matrix(dble(dm_tmp(:,:)))))
!CALL write_debug("gf img part:"//num2str(chksum_2d_matrix(aimag(dm_tmp(:,:)))))

!        CALL to_eight_digest()
!        dm_tmp(:,:,:) = 0.d0

       CALL store_green_func(greensfunc, rvecs%vecs(i_cell), &
                             n_basis, &
                             n_basis, &
                             dm_tmp(:,:,:))

if(maxval(abs(dimag(dm_tmp(:,:,:)))) > 2 .or. &
   chksum_3d_matrix(dble(dm_tmp(:,:,:))) > 10) then
CALL write_stdout("Max diff green_c"//num2str(i_cell) &
                                    //num2str(maxval(abs(dimag(dm_tmp(:,:,:))))) &
                                    //num2str(chksum_3d_matrix(dble(dm_tmp(:,:,:)))))
endif

      enddo

      DEALLOCATE(dm_tmp)

    END SUBROUTINE fouriertransform_comp_greensfunction_k_to_R

    SUBROUTINE debug_print_ks_vector(KS_eigenvalue, &
                                     KS_eigenvector_complex)
        REAL*8, INTENT(IN)  :: KS_eigenvalue(n_states,n_spin,n_k_points)
        COMPLEX*16, INTENT(IN), dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector_complex


        INTEGER :: i_state, i_spin, i_k, i_k_point


        i_k = 0
        do i_k_point = 1, n_k_points
!write(use_unit,*) i_k_point,k_phase_exx(i_k_point,:)

          if(myid == MOD(i_k_point, n_tasks)) then
             i_k = (i_k_point-1)/n_tasks + 1
!            i_k = i_k + 1
            do i_spin = 1, n_spin

                do i_state = 1, n_states

                    write(use_unit,*) "kS_eigenstate",i_k_point,i_state,i_spin, &
                                KS_eigenvector_complex(1:2,i_state,i_spin,i_k)

                enddo
            enddo
          endif
        enddo

        stop
    END SUBROUTINE

    LOGICAL FUNCTION is_energy_degenerate(KS_eigenvalue,occ_numbers, &
                                          i_state_cur, i_spin_cur,i_k_cur)
        REAL*8, INTENT(IN), DIMENSION(n_states,n_spin,n_k_points)  :: KS_eigenvalue
        REAL*8, INTENT(IN), DIMENSION(n_states, n_spin, n_k_points) :: occ_numbers
        INTEGER, INTENT(IN) :: i_state_cur,i_spin_cur,i_k_cur

        INTEGER :: i_state,i_k,i_spin
        REAL*8 :: energy_cur, occupation_cur

        energy_cur = KS_eigenvalue(i_state_cur,i_spin_cur,i_k_cur)
        occupation_cur = occ_numbers(i_state_cur,i_spin_cur,i_k_cur)

        do i_k=1,n_k_points
            do i_spin = 1,n_spin
                do i_state = 1,n_states
                    if(abs(KS_eigenvalue(i_state,i_spin,i_k)-energy_cur) < 1.d-8) then! &
!                       .AND. &
!                       abs(1.d0-(occ_numbers(i_state,i_spin,i_k)-occupation_cur))< 1.d-8) then

                        if(i_k /= i_k_cur .OR. &
                           i_spin /= i_spin_cur .OR. &
                           i_state /= i_state_cur) then
                             is_energy_degenerate = .TRUE.


!                             write(use_unit,*) occupation_cur, (occ_numbers(i_state,i_spin,i_k)-occupation_cur)


!                             write(use_unit,*) "DEGENERATE", i_state_cur,i_spin_cur,i_k_cur, &
!                                                      i_state,i_spin,i_k

                             return
                        endif
                    endif
                enddo
            enddo
        enddo

        is_energy_degenerate = .FALSE.

    END FUNCTION is_energy_degenerate

END MODULE cRPA_calculation_greensfunc
