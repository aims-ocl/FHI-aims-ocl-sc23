!****s* FHI-aims/fciqmc/evaluate_w_vector
!  NAME
!    evaluate_w_vector
!  SYNOPSIS 

      subroutine evaluate_w_vector_1a1b()

!  PURPOSE
!
!  The '''evaluate_w_vector''' calculate w_vector on-the-fly for a given c_vect
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
  use dimensions
  use runtime_choices
  use physics
  use species_data
  use constants
  use gw_para
  use grids
  use geometry
  use mpi_utilities
  use synchronize_mpi
  use sbt_overlap_aims
  use localorb_io
  use basis
  use prodbas
  use grids
  use hartree_fock
  use timing
  use mpi_tasks
  use mpi_utilities
  use timing
  use fciqmc_module

    implicit none
    ! temp variables
    integer, dimension(:,:), allocatable :: occ_num_1, occ_num_2
    ! temp indices
    integer :: i_spin
    integer :: ii_state, jj_state, kk_state, ll_state
    integer :: i_state, j_state, a_state, b_state
    integer :: c_state, d_state, k_state, l_state
    integer :: i_ci, j_ci
    integer :: i_ci_start, j_ci_start
    integer :: errnum
    !character*128 :: info_str
    logical :: exit_flag = .false.
    real*8  :: w_1, w_2, w_3, w_4, w_5
    

    ! --------------------------------------------------------------
    ! First for the substitution in both spaces
    ! --------------------------------------------------------------
    i_ci = n_configuration_1a + n_configuration_1b
    j_ci = 0
    w_vect(1) = 0.0d0
    do i_state = 1, n_elec(1), 1
      do a_state = n_elec(1)+1, n_states, 1
        !j_ci_start = i_ci
        do j_state = 1, n_elec(2), 1
          do b_state = n_elec(2)+1, n_states, 1
            occ_num_1 = occ_num_0
            occ_num_1(i_state,1) = 0 ! occ_alpha
            occ_num_1(j_state,2) = 0 ! occ_beta
            occ_num_1(a_state,1) = 1 ! vir_alpha
            occ_num_1(b_state,2) = 1 ! vir_beta
            i_ci = i_ci + 1
            w_vect(1) = w_vect(1) + &
                integral_4ks(&
                i_state,a_state,1,&
                j_state,b_state,1) * c_vect(i_ci)
            !write(use_unit,'(2X,A,4(I3),F16.8)') "igor debug", &
            !    i_state, j_state, a_state, b_state, w_vect(1)
            w_1 = 0.0d0
            !j_ci = j_ci_start
            do c_state = n_elec(1)+1, n_states, 1
              do d_state = n_elec(2)+1, n_states, 1
                !j_ci = j_ci + 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(j_state,2) = 0 ! occ_beta
                occ_num_1(c_state,1) = 1 ! vir_alpha
                occ_num_1(d_state,2) = 1 ! vir_beta
                call seek_index(occ_num_1,3,j_ci)
                w_1 = w_1 + &
                    integral_4ks(&
                    a_state,c_state,1,&
                    b_state,d_state,1) * c_vect(j_ci)
              enddo
            enddo
            w_2 = 0.0d0
            do k_state = 1, n_elec(1), 1
              do l_state = 1, n_elec(2), 1
                occ_num_1 = occ_num_0
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(l_state,2) = 0 ! occ_beta
                occ_num_1(a_state,1) = 1 ! vir_alpha
                occ_num_1(b_state,2) = 1 ! vir_beta
                call seek_index(occ_num_1,3,j_ci)
                w_2 = w_2 + &
                    integral_4ks(&
                    k_state,i_state,1,&
                    l_state,j_state,1) * c_vect(j_ci)
                !write(use_unit,'(2X,A10,6(I3))') "igor debug", &
                !    i_state, a_state, j_state, b_state, i_ci
                !write(use_unit,'(12X,5(I3),F16.8)') k_state, a_state, &
                !    l_state, b_state, j_ci, w_2
              enddo
            enddo
            w_2 = w_2
            w_3 = 0.0d0
            !--------------------------------------
            ! for the first term in w_3
            !--------------------------------------
            ! for the first term in the first term
            ! exciate and occupy both in the alpha spin space
            ! due to antisymmetry, this term should be separated to four terms
            ! 1) k<i and c<a
            do k_state = 1, i_state-1, 1
              do c_state = n_elec(1)+1, a_state-1, 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(a_state,1) = 1 ! vir_alpha
                occ_num_1(c_state,1) = 1 ! vir_alpha
                call seek_index(occ_num_1,4,j_ci)
                w_3 = w_3 + (&
                    -integral_4ks(&
                    k_state,c_state,1, &
                    b_state,j_state,1) ) * c_vect(j_ci)
              enddo
            enddo
            ! 2) k<i and c>a
            do k_state = 1, i_state-1, 1
              do c_state = a_state+1, n_states, 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(c_state,1) = 1 ! vir_alpha
                occ_num_1(a_state,1) = 1 ! vir_alpha
                call seek_index(occ_num_1,4,j_ci)
                w_3 = w_3 + (&
                    integral_4ks(&
                    k_state,c_state,1, &
                    b_state,j_state,1) ) * c_vect(j_ci)
              enddo
            enddo
            ! 3) k>i and c<a
            do k_state = i_state+1, n_states, 1
              do c_state = n_elec(1)+1, a_state-1, 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(c_state,1) = 1 ! vir_alpha
                occ_num_1(a_state,1) = 1 ! vir_alpha
                call seek_index(occ_num_1,4,j_ci)
                w_3 = w_3 + (&
                    integral_4ks(&
                    k_state,c_state,1, &
                    b_state,j_state,1) ) * c_vect(j_ci)
              enddo
            enddo
            ! 4) k>i and c>a
            do k_state = i_state+1, n_states, 1
              do c_state = a_state+1, n_states, 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(c_state,1) = 1 ! vir_alpha
                occ_num_1(a_state,1) = 1 ! vir_alpha
                call seek_index(occ_num_1,4,j_ci)
                w_3 = w_3 + (&
                    -integral_4ks(&
                    k_state,c_state,1, &
                    b_state,j_state,1) ) * c_vect(j_ci)
              enddo
            enddo
            ! for the second term in the first term
            do k_state = 1, n_elec(2), 1
              do c_state = n_elec(2)+1, n_states, 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(k_state,2) = 0 ! occ_beta
                occ_num_1(a_state,1) = 1 ! vir_alpha
                occ_num_1(c_state,2) = 1 ! vir_beta
                call seek_index(occ_num_1,3,j_ci)
                w_3 = w_3 + (&
                    integral_4ks(&
                    k_state,j_state,1, &
                    b_state,c_state,1) - &
                    integral_4ks(&
                    k_state,c_state,1, &
                    b_state,j_state,1)) * c_vect(j_ci)
              enddo
            enddo
            !--------------------------------------
            ! for the second term in w_3
            !--------------------------------------
            do k_state = 1, n_elec(2), 1
              do c_state = n_elec(1)+1, n_states, 1
                occ_num_1 = occ_num_0
                occ_num_1(i_state,1) = 0 ! occ_alpha
                occ_num_1(k_state,2) = 0 ! occ_beta
                occ_num_1(c_state,1) = 1 ! vir_alpha
                occ_num_1(b_state,2) = 1 ! vir_beta
                call seek_index(occ_num_1,3,j_ci)
                w_3 = w_3 + (&
                    integral_4ks(&
                    k_state,j_state,1, &
                    a_state,c_state,1)) * c_vect(j_ci)
              enddo
            enddo
            !--------------------------------------
            ! for the third term in w_3
            !--------------------------------------
            do k_state = 1, n_elec(1), 1
              do c_state = n_elec(2)+1, n_states, 1
                occ_num_1 = occ_num_0
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(j_state,2) = 0 ! occ_beta
                occ_num_1(a_state,1) = 1 ! vir_alpha
                occ_num_1(c_state,2) = 1 ! vir_beta
                call seek_index(occ_num_1,3,j_ci)
                w_3 = w_3 + (&
                    integral_4ks(&
                    k_state,i_state,1, &
                    b_state,c_state,1)) * c_vect(j_ci)
              enddo
            enddo
            !--------------------------------------
            ! for the fourth term in w_3
            !--------------------------------------
            ! for the first term in the fourth term
            do k_state = 1, n_elec(1), 1
              do c_state = n_elec(1)+1, n_states, 1
                occ_num_1 = occ_num_0
                occ_num_1(k_state,1) = 0 ! occ_alpha
                occ_num_1(j_state,2) = 0 ! occ_beta
                occ_num_1(c_state,1) = 1 ! vir_alpha
                occ_num_1(b_state,2) = 1 ! vir_beta
                call seek_index(occ_num_1,3,j_ci)
                w_3 = w_3 + (&
                    integral_4ks(&
                    k_state,i_state,1, &
                    a_state,c_state,1) -&
                    integral_4ks(&
                    k_state,c_state,1, &
                    a_state,i_state,1) )* c_vect(j_ci)
              enddo
            enddo
            w_3 = - w_3
            w_4 = 0.0d0
            w_5 = 0.0d0
            !write(use_unit,'(2X,I3,4(F16.8))') i_ci+1, b_vect(i_ci),w_1, w_2, w_3
            w_vect(i_ci+1) = b_vect(i_ci)/Norm_A + w_1 + &
                w_2 + w_3 + w_4 + w_5
          enddo
        enddo
      enddo
    enddo ! two-electron excitations in both spin spaces
    !! --------------------------------------------------------------
    !! Second for the substitution in alpha spaces
    !! --------------------------------------------------------------
    !i_ci = n_configuration_1a + n_configuration_2a + n_configuration_1a1b
    !do i_state = 1, n_elec(1), 1
    !  do a_state = n_elec(1)+1, n_states, 1
    !    ! consider the combination
    !    ! (n_elec(1),2)=n_elec(1)!/(2!*(n_elec(1)-2)!
    !    !              =n_elec(1)*(n_elec(1)-1)/2
    !    do j_state = i_state+1, n_elec(1), 1 
    !      do b_state = j_state+1, n_states, 1
    !        occ_num_1 = occ_num_0
    !        occ_num_1(i_state,1) = 0
    !        occ_num_1(j_state,1) = 0
    !        occ_num_1(a_state,1) = 1
    !        occ_num_1(b_state,1) = 1
    !        i_ci = i_ci + 1
    !        w_vect(1) = w_vect(1) + &
    !            (&
    !            integral_4ks(&
    !            i_state,a_state,1,&
    !            j_state,b_state,1) -&
    !            integral_4ks(&
    !            i_state,b_state,1,&
    !            j_state,a_state,1) )&
    !            * c_vect(i_ci)
    !        w_1 = 0.0d0
    !        do c_state = n_elec(1)+1, n_states, 1
    !          do d_state = c_state+1, n_states, 1
    !            occ_num_1 = occ_num_0
    !            occ_num_1(i_state,1) = 0 ! occ_alpha
    !            occ_num_1(j_state,1) = 0 ! occ_alpha
    !            occ_num_1(c_state,1) = 1 ! vir_alpha
    !            occ_num_1(d_state,1) = 1 ! vir_alpha
    !            call seek_index(occ_num_1,4,j_ci)
    !            w_1 = w_1 + (&
    !                integral_4ks(&
    !                a_state,c_state,1,&
    !                b_state,d_state,1) - &
    !                integral_4ks(&
    !                a_state,d_state,1,&
    !                b_state,c_state,1) &
    !                )* c_vect(j_ci)
    !          enddo
    !        enddo
    !        w_2 = 0.0d0
    !        do k_state = 1, n_elec(1), 1
    !          do l_state = k_state+1, n_states, 1
    !            occ_num_1 = occ_num_0
    !            occ_num_1(k_state,1) = 0 ! occ_alpha
    !            occ_num_1(l_state,1) = 0 ! occ_beta
    !            occ_num_1(a_state,1) = 1 ! vir_alpha
    !            occ_num_1(b_state,1) = 1 ! vir_beta
    !            call seek_index(occ_num_1,4,j_ci)
    !            w_1 = w_1 + (&
    !                integral_4ks(&
    !                k_state,i_state,1,&
    !                l_state,j_state,1) - &
    !                integral_4ks(&
    !                k_state,j_state,1,&
    !                l_state,i_state,1) &
    !                )* c_vect(j_ci)
    !          enddo
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !enddo ! two-electron excitations in the alpha spin space

end subroutine evaluate_w_vector_1a1b
