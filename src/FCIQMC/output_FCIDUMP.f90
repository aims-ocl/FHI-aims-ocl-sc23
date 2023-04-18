!****s* FHI-aims/fciqmc/output_FCIDUMP
!  NAME
!    evaluate_cisd
!  SYNOPSIS 

      subroutine output_FCIDUMP()

!  PURPOSE
!
!  The '''output_FCIDUMP.f90''' evaluate w_vect(1)
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Jan Kloppenburg
!  SOURCE
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

    implicit none
    ! temp variables
    integer, dimension(:,:), allocatable :: occ_num_1
    ! temp indices
    integer :: i_state, j_state, k_state, l_state, i_spin, j_spin
    integer :: i_ci
    integer :: errnum
    integer :: unit_fci
    !character*128 :: info_str
    logical :: exit_flag
    character(len=128) :: string_form

    real*8  :: ddot
    real*8  :: integral_4ks_tmp

    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 (1D) in CI')
    endif

    occ_num_1(:,1)=1
    occ_num_1(:,2)=1

101 FORMAT(1X,E23.16,4I4)
102 FORMAT(1X,'&FCI NORB=',I3,',NELEC=',I3,',MS2=',I2,',')
103 FORMAT(2X,'ORBSYM=',12(I1,','))
104 FORMAT(2X,'ISYM=',I1)
105 FORMAT(2X,'ISYM=',I1,' UHF=.TRUE.')
    open(unit_fci,File='FCIDUMP',STATUS='NEW',FORM='FORMATTED')
    if (n_spin .eq. 1) then
        write(unit_fci,102) n_states, n_elec(1)+n_elec(2), n_elec(1)-n_elec(2)
        write(string_form,'(I3)') n_states
        string_form = "(2X,'ORBSYM=',"//trim(string_form)//"(I1,','))"
        !write(unit_fci,*) string_form
        write(unit_fci,trim(string_form)) occ_num_1(:,1)
        !write(unit_fci,103) occ_num_1(:,1)
        write(unit_fci,104) 0
    else
        write(unit_fci,102) n_states*2, n_elec(1)+n_elec(2), n_elec(1)-n_elec(2)
        write(string_form,'(I3)') n_states*2
        string_form = "(2X,'ORBSYM=',"//trim(string_form)//"(I1,','))"
        write(unit_fci,trim(string_form)) occ_num_1(:,:)
        write(unit_fci,105) 0
    endif
    write(unit_fci,*) '&END'
    if (n_spin==1) then
        do i_state = 1, n_states, 1
          do j_state = 1, i_state, 1
            do k_state = 1, n_states, 1
              do l_state = 1, k_state, 1
                integral_4ks_tmp = ddot(n_basbas, &
                    ovlp_3ks(:,l_state,k_state,1),1, &
                    ovlp_3ks(:,j_state,i_state,1),1)
                if (abs(integral_4ks_tmp).gt.Threshold_Int) then
                    write(unit_fci,101) integral_4ks_tmp, &
                        i_state, j_state, k_state, l_state
                endif
              enddo
            enddo
          enddo
        enddo
        do i_state = 1, n_states, 1
          do j_state = 1, i_state, 1
             if (abs(integral_2ks(j_state,i_state,1)).gt.&
                 Threshold_Int) then
                write(unit_fci,101) integral_2ks(&
                    j_state,i_state,1&
                    ), & 
                    i_state, j_state, 0, 0
            endif
          enddo
        enddo
    else
        do i_spin = 1, n_spin, 1
          do i_state = 1, n_states, 1
            do j_state = 1, n_states, 1
              do j_spin = 1, n_spin, 1
                do k_state = 1, n_states, 1
                  do l_state = 1, n_states, 1
                    integral_4ks_tmp = ddot(n_basbas, &
                        ovlp_3ks(:,l_state,k_state,j_spin),1, &
                        ovlp_3ks(:,j_state,i_state,i_spin),1)
                    if (abs(integral_4ks_tmp).gt.Threshold_Int) then
                        write(unit_fci,101) integral_4ks_tmp, &
                            (i_state-1)*2+i_spin, &
                            (j_state-1)*2+i_spin, &
                            (k_state-1)*2+j_spin, &
                            (l_state-1)*2+j_spin
                            !i_state*2-i_spin+1, &
                            !j_state*2-i_spin+1, &
                            !k_state*2-j_spin+1, &
                            !l_state*2-j_spin+1
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        do i_spin = 1, n_spin, 1
          do i_state = 1, n_states, 1
            do j_state = 1, i_state, 1
               if (abs(integral_2ks(j_state,i_state,i_spin)).gt.&
                   Threshold_Int) then
                  write(unit_fci,101) integral_2ks(&
                      j_state,i_state,i_spin), & 
                      (i_state-1)*2+i_spin, &
                      (j_state-1)*2+i_spin, &
                      0, 0
                      !i_state*2-i_spin+1, &
                      !j_state*2-i_spin+1, &
                      !0, 0
              endif
            enddo
          enddo
        enddo
        do i_spin = 1, n_spin, 1
          do i_state = 1, n_states, 1
            write(unit_fci,101) KS_eigenvalue(&
                i_state,i_spin,1), & 
                (i_state-1)*2+i_spin, &
                0,  0, 0
                !i_state*2-i_spin+1, &
                !0,  0, 0
          enddo
        enddo
        !call aims_stop('************   Spin UNrestricted calculations are not yet possible!')
    endif
    
    !
    write(unit_fci,101) en_ion_ion, 0, 0, 0, 0

    close(unit_fci)

    if (myid .eq. 1) then
        write(use_unit,'(2X,A)') 'Two- and four-integrals have been prepared in FCIDUMP for NECI'
    endif

end subroutine output_FCIDUMP

