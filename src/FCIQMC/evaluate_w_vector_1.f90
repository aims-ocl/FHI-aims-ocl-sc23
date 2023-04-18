!****s* FHI-aims/fciqmc/evaluate_w_vector_1
!  NAME
!    evaluate_w_vector_1
!  SYNOPSIS 

      subroutine evaluate_w_vector_1()

!  PURPOSE
!
!  The '''evaluate_w_vector_1''' evaluate w_0
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
  use synchronize_mpi
  use fciqmc_module

    implicit none
    ! temp variables
    integer, dimension(:,:), allocatable :: occ_num_1
    ! temp indices
    integer :: i_state, j_state, a_state, b_state
    integer :: i_ci
    integer :: errnum
    !character*128 :: info_str
    real*8  :: ddot
    integer :: tmp_spin(2)
    integer :: index_start(4)

    tmp_spin = 1
    if (n_spin .eq. 2) tmp_spin(2) = 2
    
    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 (1D) in CI')
    endif

    w_0 = 0.0d0

    ! ==============================================================
    ! now construct matrix elements between two-electron substutions
    ! ==============================================================
    ! --------------------------------------------------------------
    ! First for the substitution in both spaces
    ! --------------------------------------------------------------

    i_ci = memory_table(1,2,3,myid+1)-1
    index_start = index_table(:,3,myid+1)
    do i_state = index_start(1), n_elec(1), 1
      do a_state = index_start(2), n_states, 1
        do j_state = index_start(3), n_elec(2), 1
          do b_state = index_start(4), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,3,myid+1)) goto 668

            occ_num_1 = occ_num_0
            occ_num_1(i_state,1) = 0  !occ_alpha
            occ_num_1(j_state,2) = 0  !occ_beta
            occ_num_1(a_state,1) = 1  !vir_alpha
            occ_num_1(b_state,2) = 1  !vir_beta
            w_0 = w_0 + (&
                ddot(n_basbas,&
                ovlp_3ks(:,i_state,a_state,tmp_spin(1)),1,&
                ovlp_3ks(:,j_state,b_state,tmp_spin(2)),1)&
                )*c_vect(i_ci)
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
    i_ci = memory_table(1,2,4,myid+1)-1
    index_start = index_table(:,4,myid+1)
    do i_state = index_start(1), n_elec(1), 1
      do j_state = index_start(3), n_elec(1), 1 
        do a_state = index_start(2), n_states, 1
          do b_state = index_start(4), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,4,myid+1)) goto 669

            occ_num_1 = occ_num_0
            occ_num_1(i_state,1) = 0  !occ_alpha
            occ_num_1(j_state,1) = 0  !occ_alpha
            occ_num_1(a_state,1) = 1  !vir_alpha
            occ_num_1(b_state,1) = 1  !vir_alpha
            w_0 = w_0 + (&
                ddot(n_basbas,&
                ovlp_3ks(:,i_state,a_state,tmp_spin(1)),1,&
                ovlp_3ks(:,j_state,b_state,tmp_spin(1)),1) - &
                ddot(n_basbas,&
                ovlp_3ks(:,i_state,b_state,tmp_spin(1)),1,&
                ovlp_3ks(:,j_state,a_state,tmp_spin(1)),1)&
                )*c_vect(i_ci)
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
    i_ci = memory_table(1,2,5,myid+1)-1
    index_start = index_table(:,5,myid+1)
    ! index_table(:,...) are ordered by each pair of exciations
    do i_state = index_start(1), n_elec(2), 1
      do j_state = index_start(3), n_elec(2), 1 
        do a_state = index_start(2), n_states, 1
          do b_state = index_start(4), n_states, 1

            i_ci = i_ci + 1
            if (i_ci .gt. memory_table(2,2,5,myid+1)) goto 670

            occ_num_1 = occ_num_0
            occ_num_1(i_state,2) = 0  !occ_beta
            occ_num_1(j_state,2) = 0  !occ_beta
            occ_num_1(a_state,2) = 1  !vir_beta
            occ_num_1(b_state,2) = 1  !vir_beta
            w_0 = w_0 + (&
                ddot(n_basbas,&
                ovlp_3ks(:,i_state,a_state,tmp_spin(2)),1,&
                ovlp_3ks(:,j_state,b_state,tmp_spin(2)),1) - &
                ddot(n_basbas,&
                ovlp_3ks(:,i_state,b_state,tmp_spin(2)),1,&
                ovlp_3ks(:,j_state,a_state,tmp_spin(2)),1)&
                )*c_vect(i_ci)
          enddo
          if (a_state .eq. n_states) then
              index_start(4) = n_elec(2)+2
          else
              index_start(4) = a_state+2
          endif
        enddo
        index_start(2) = n_elec(2)+1 
      enddo
      index_start(3) = i_state+2
    enddo ! two-electron excitations in the beta spin space
670 continue

    call sync_real_number(w_0)

end subroutine evaluate_w_vector_1

