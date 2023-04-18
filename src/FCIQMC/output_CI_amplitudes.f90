!****s* FHI-aims/RRS-PBC/output_CI_amplitudes
!  NAME
!   output_CI_amplitudes
!  SYNOPSIS

    subroutine output_CI_amplitudes_v02(scaled_factor,coeff_threshold)
!  USES
    use dimensions
    use runtime_choices
    use fciqmc_module
!  PURPOSE
!   The subroutine prints out the array '''p_array'''
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!                                                                  
!  SEE ALSO
!    
!  COPYRIGHT
!   
!  HISTORY
!    
!  SOURCE

    real*8  :: scaled_factor, coeff_threshold

    integer :: i_state, j_state, a_state, b_state
    integer :: i_ci

    real*8  :: tmp_c
    integer :: index_start(4)

    !if (myid .eq. 0) then
    !    write(use_unit,'(4X,A,4(A5),A21)') &
    !        '| Spin Case ', 'I','J','A','B','Value     '
    !endif


    if (flag_single_excitation) then
        ! ==============================================================
        ! now construct matrix elements between one-electron substutions
        ! ==============================================================
        ! First for the substitution in the alpha space
        ! --------------------------------------------------------------
        !i_ci = 0
        !do i_state = i_start_ci, n_elec(1), 1
        !  do a_state = n_elec(1)+1, n_states, 1
        i_ci = memory_table(1,2,1,myid+1)-1
        index_start = index_table(:,1,myid+1)
        do i_state = index_start(1), n_elec(1), 1
          do a_state = index_start(2), n_states, 1
            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,1,myid+1)) goto 666
            tmp_c = c_vect(i_ci)*scaled_factor
            !write(use_unit,*) "igor debug 1",i_ci, tmp_c, scaled_factor
            if (abs(tmp_c).ge.coeff_threshold) then
                write(use_unit,'(4X,A,I5,A5,I5,A5,E21.8)') '|    AA     ',&
                    i_state,'',a_state,'', tmp_c
            endif
          enddo
          index_start(2) = n_elec(1)+1
        enddo ! one-electron excitations in the alpha spin space
666     continue

        ! --------------------------------------------------------------
        ! Second for the substitution in the beta space
        ! --------------------------------------------------------------
        if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
            !i_ci = n_configuration_1a
            !do i_state = i_start_ci, n_elec(2), 1
            !  do a_state = n_elec(2)+1, n_states, 1
            i_ci = memory_table(1,2,2,myid+1)-1
            index_start = index_table(:,2,myid+1)
            do j_state = index_start(1), n_elec(2), 1
              do b_state = index_start(2), n_states, 1
                i_ci = i_ci + 1
                if (i_ci .gt. memory_table(2,2,2,myid+1)) goto 667
                tmp_c = c_vect(i_ci)*scaled_factor
                if (abs(tmp_c).ge.coeff_threshold) then
                    write(use_unit,'(4X,A,A5,I5,A5,I5,E21.8)') '|    BB     ',&
                        '',j_state,'',b_state,tmp_c
                endif
              enddo
              index_start(2) = n_elec(2)+1
            enddo ! one-electron excitations in the beta spin space
667         continue
        endif
    endif ! if (flag_single_excitation) then
    ! ==============================================================
    ! now construct matrix elements between two-electron substutions
    ! ==============================================================
    ! --------------------------------------------------------------
    ! First for the substitution in both spaces
    ! --------------------------------------------------------------
    !i_ci = n_configuration_1a + n_configuration_1btest.log.cisd.FC.v02
    !do i_state = i_start_ci, n_elec(1), 1
    !  do a_state = n_elec(1)+1, n_states, 1
    !    do j_state = i_start_ci, n_elec(2), 1
    !      do b_state = n_elec(2)+1, n_states, 1
    i_ci = memory_table(1,2,3,myid+1)-1
    index_start = index_table(:,3,myid+1)
    do i_state = index_start(1), n_elec(1), 1
      do a_state = index_start(2), n_states, 1
        do j_state = index_start(3), n_elec(2), 1
          do b_state = index_start(4), n_states, 1
            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,3,myid+1)) goto 668
            tmp_c = c_vect(i_ci)*scaled_factor
            !write(use_unit,*) "igor debug 2",i_ci, tmp_c
            if (abs(tmp_c).ge.coeff_threshold) then
               write(use_unit,'(4X,A,4(I5),E21.8)') '|   ABAB    ',&
                   i_state,a_state,j_state,b_state,tmp_c
           endif
          enddo
          index_start(4) = n_elec(2)+1
        enddo
        index_start(3) = i_start_ci 
      enddo
      index_start(2) = n_elec(1)+1 
    enddo ! two-electron excitations in both spin spaces
668 continue
    ! --------------------------------------------------------------
    ! Second for the substitution in the alpha spin space
    ! --------------------------------------------------------------
    !i_ci = n_configuration_1a + n_configuration_1b + &
    !       n_configuration_1a1b
    !do i_state = i_start_ci, n_elec(1), 1
    !  do j_state = i_state+1, n_elec(1), 1 
    !    do a_state = n_elec(1)+1, n_states, 1
    !      do b_state = a_state+1, n_states, 1
    i_ci = memory_table(1,2,4,myid+1)-1
    index_start = index_table(:,4,myid+1)
    do i_state = index_start(1), n_elec(1), 1
      do j_state = index_start(3), n_elec(1), 1 
        do a_state = index_start(2), n_states, 1
          do b_state = index_start(4), n_states, 1
            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,4,myid+1)) goto 669
            tmp_c = c_vect(i_ci)*scaled_factor
            if (abs(tmp_c).ge.coeff_threshold) then
                write(use_unit,'(4X,A,4(I5),E21.8)') '|   AAAA    ',&
                    i_state,a_state,j_state,b_state,tmp_c
            endif
          enddo
          if(a_state .eq. n_states) then
              index_start(4) = n_elec(1)+2
          else
              index_start(4) = a_state+2
          endif
        enddo
        index_start(2) = n_elec(1)+1 
      enddo
      index_start(3) = i_state+2
    enddo ! two-electron excitations in the alpha spin space
669 continue
    ! --------------------------------------------------------------
    ! Second for the substitution in the beta spin space
    ! --------------------------------------------------------------
    if (n_elec(1) .ne. n_elec(2) .or. n_spin .eq. 2) then
        !i_ci = n_configuration_1a + n_configuration_1b + &
        !       n_configuration_1a1b + n_configuration_1a
        !do i_state = i_start_ci, n_elec(2), 1
        !  do j_state = i_state+1, n_elec(2), 1 
        !    do a_state = n_elec(2)+1, n_states, 1
        !      do b_state = a_state+1, n_states, 1
        i_ci = memory_table(1,2,5,myid+1)-1
        index_start = index_table(:,5,myid+1)
        do i_state = index_start(1), n_elec(2), 1
          do j_state = index_start(3), n_elec(2), 1 
            do a_state = index_start(2), n_states, 1
              do b_state = index_start(4), n_states, 1
                i_ci = i_ci + 1
                if (i_ci .gt. memory_table(2,2,5,myid+1)) goto 670
                tmp_c = c_vect(i_ci)*scaled_factor
                if (abs(tmp_c).ge.coeff_threshold) then
                    write(use_unit,'(4X,A,4(I5),E21.8)') '|   BBBB    ',&
                        i_state,a_state,j_state,b_state,tmp_c
                endif
              enddo
              if(a_state .eq. n_states) then
                  index_start(4) = n_elec(2)+2
              else
                  index_start(4) = a_state+2
              endif
           enddo
           index_start(2) = n_elec(2)+1 
         enddo
         index_start(3) = i_state+2
        enddo ! two-electron excitations in the beta spin space
670     continue
    endif

end subroutine output_CI_amplitudes_v02

