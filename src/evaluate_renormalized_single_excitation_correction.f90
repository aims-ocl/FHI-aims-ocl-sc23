!****s*  FHI-aims/evaluate_renormalized_single_excitation_correction
!  NAME
!    evaluate_renormalized_single_excitation_correction
!  SYNOPSIS

     subroutine evaluate_renormalized_single_excitation_correction &
         ( n_KS_states,n_homo, &
           occ_numbers, ovlp_3KS, &
           KS_eigenvalue,KS_eigenvector, &
           xc_matr, en_rse_full &
          )

!  PURPOSE
!  Subroutine evaluate_renormalized_single_excitation_correction evaluates a summation of 
!  a selected type of single excitation contribution to the correlation energy. 
!  This is zero for a HF reference
!  state, but non-zero for a KS reference state.
!
!  There is a full renormalized single excitation contribution, namely, the non-diagonal
!  terms which was neglected before, is also included the following treatment.
!
!  It is given by \sum_{ia} |<i|V_HFx - V_KSxc|a>|^2 / (E_i - E_a)
!  where V_HFx is the HF exchange operator,  V_KSxc is the KS
!  exchange-correlation potential operator, and the summation of i, a 
!  runs over the occupied states and unoccupied states respectively.
! USES

      use runtime_choices, only: flag_xc, hybrid_coeff, safe_minimum
      use dimensions
      use runtime_choices
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  KS_eigenvector(n_basis,n_states,n_spin)
      real*8  ovlp_3KS(n_loc_prodbas,n_states,n_KS_states,n_spin)
      real*8  xc_matr(n_basis,n_basis,n_spin)
      real*8  en_rse_full


! INPUTS
!  o n_KS_states -- number of KS states included in the calculations
!  o occ_numbers -- the occupation number 
!  o KS_eigenvector -- Single-particle KS eigenvector
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!  o xc_matr -- the matrix form the KS exchange-correlation potential operator
!            -- within the atomic orbital basis
!
! OUTPUT      
!  o en_rse_full -- a summation of a ladder-type single excitation contribution to the 
!              correlation energy, now also the non-diagonal terms are included
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




!  local variables

      real*8  tmp_numerator
      real*8, dimension(:,:), allocatable :: HF_ex_matr  
      real*8, dimension(:,:), allocatable :: KS_xc_matr  
      real*8, dimension(:,:), allocatable :: aux_xc_matr
      real*8, dimension(:), allocatable :: HF_ev_occ
      real*8, dimension(:), allocatable :: HF_ev_unocc
      real*8, dimension(:,:), allocatable :: HF_ham_occ
      real*8, dimension(:,:), allocatable :: HF_ham_unocc
      real*8, dimension(:,:), allocatable :: trans_matr_occ
      real*8, dimension(:,:), allocatable :: trans_matr_unocc
      real*8, dimension(:,:), allocatable :: delta_v_num
      real*8, dimension(:,:), allocatable :: tmp_delta_v_num
      integer n_homo(n_spin)
      integer n_homo_max, n_homo_min
      integer n_occ, n_unocc
      integer n_nonsingular
      real*8 ddot
!     counters

      integer :: i_state, i_state_1
      integer :: j_state, j_state_1
      integer :: k_state
      integer :: i_spin


!     begin work

      if(flag_xc.eq.0) then
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A)') "HF reference state: single excitation contribution is zero"
        endif
        en_rse_full=0.d0
        return
      endif 

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the renormalized single excitation (full) correction ..."
      endif

!  determine the HOMO level
      do i_spin = 1, n_spin
       do i_state = 1, n_states
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         endif
       enddo
      enddo

      n_homo_max=max(n_homo(1),n_homo(n_spin))
      n_homo_min=min(n_homo(1),n_homo(n_spin))

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      allocate(HF_ex_matr(n_states, n_states))
      allocate(KS_xc_matr(n_states, n_states))
      allocate(aux_xc_matr(n_states, n_basis))

      en_rse_full = 0.d0
      do i_spin = 1, n_spin

        if(n_homo(i_spin).eq.0) cycle

        HF_ex_matr(:,:) = 0.d0
        do k_state = 1, n_homo(i_spin), 1

          do j_state = 1, n_states, 1
           do i_state = 1, n_states, 1

             HF_ex_matr(i_state, j_state) = &
              HF_ex_matr(i_state, j_state) - &
                 ddot( n_loc_prodbas, &
                      ovlp_3KS(1, i_state, k_state, i_spin), 1, &
                      ovlp_3KS(1, j_state, k_state, i_spin), 1 &
                     ) * &
                 occ_numbers(k_state,i_spin,1)*dble(n_spin)/2.d0
!       if(i_state.eq.1) then
!        write(use_unit,*) k_state, exchange_energy(i_state,i_spin),
!     +           "ovlp",   ovlp_3KS(:,k_state,i_state,i_spin),
!     +           "occ",    occ_numbers(k_state,i_spin,1)
!       endif
           enddo
         enddo
        enddo

        call sync_matrix(HF_ex_matr, n_states, n_states)
  
        KS_xc_matr(:,:) = 0.d0
        aux_xc_matr = 0.0d0
       
        call dgemm('T', 'N', n_states, n_basis, &
                  n_basis, 1.0d0, &
                  KS_eigenvector(1,1,i_spin), n_basis, &
                  xc_matr(1,1,i_spin), n_basis, 0.d0, &
                  aux_xc_matr(1,1), &
                  n_states)

!     second multiplication between eigenvector and temporary XC matrix.
        call dgemm('N', 'N', n_states, n_states, &
                   n_basis, 1.0d0, &
                   aux_xc_matr(1,1), n_states, &
                   KS_eigenvector(1,1,i_spin), &
                   n_basis, 0.d0, KS_xc_matr, &
                   n_states)

!  Build the Hartree-Fock Hamiltonian for the occupied and unoccupied blocks separately         
      n_occ=n_homo(i_spin)
      n_unocc=n_states-n_occ
      allocate(HF_ham_occ(n_occ,n_occ))
      allocate(HF_ham_unocc(n_unocc,n_unocc))
      allocate(trans_matr_occ(n_occ,n_occ))
      allocate(trans_matr_unocc(n_unocc,n_unocc))
      allocate(HF_ev_occ(n_occ))
      allocate(HF_ev_unocc(n_unocc))
      allocate(delta_v_num(n_occ,n_unocc))
      allocate(tmp_delta_v_num(n_occ,n_unocc))

       HF_ham_occ(:,:) = 0.d0
       do i_state = 1, n_occ, 1
         do j_state = 1, n_occ, 1
           HF_ham_occ(j_state,i_state) = &
             HF_ex_matr(j_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(j_state,i_state)
         enddo
         HF_ham_occ(i_state,i_state) = HF_ham_occ(i_state,i_state) + & 
         KS_eigenvalue(i_state,i_spin)
       enddo

       call diagonalize_auxmat_lapack &
          (  n_occ, HF_ham_occ, safe_minimum, &
             -1.d9, n_nonsingular, HF_ev_occ, trans_matr_occ, &
             "RSE-HF-Occ" &
           )
!       if(myid.eq.0) then
!         write(use_unit,*) "RSE-HF-Occ"
!         write(use_unit,*) "n_nonsingular", n_nonsingular
!         do i_state = 1, n_occ
!           write(use_unit,'(I4,3f18.8)') i_state, HF_ev_occ(i_state), KS_eigenvalue(i_state,i_spin)
!         enddo
!       endif
            
       HF_ham_unocc(:,:) = 0.d0
       do i_state = n_occ+1, n_states, 1
         i_state_1 = i_state - n_occ
         do j_state = n_occ+1, n_states, 1
           j_state_1 = j_state - n_occ
           HF_ham_unocc(j_state_1,i_state_1) = &
             HF_ex_matr(j_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(j_state,i_state)
         enddo
!         if(myid.eq.0) then
!           write(use_unit,'(I4,2f18.8)') i_state, HF_ham_unocc(i_state_1,i_state_1)
!         endif
          HF_ham_unocc(i_state_1,i_state_1) = HF_ham_unocc(i_state_1,i_state_1) + & 
          KS_eigenvalue(i_state,i_spin)
!         if(myid.eq.0) then
!           write(use_unit,'(I4,2f18.8)') i_state, HF_ham_unocc(i_state_1,i_state_1), &
!                    KS_eigenvalue(i_state,i_spin)
!         endif
      enddo

      call diagonalize_auxmat_lapack &
          (  n_unocc, HF_ham_unocc, safe_minimum, &
            -1.d9, n_nonsingular, HF_ev_unocc, trans_matr_unocc, &
             "RSE-HF-Unocc" &
           )
       
!       if(myid.eq.0) then
!         write(use_unit,*) "RSE-HF-Unocc"
!         write(use_unit,*) "n_nonsingular", n_nonsingular
!         do i_state = 1, n_unocc
!           write(use_unit,'(I4,2f18.8)') i_state+n_occ, HF_ev_unocc(i_state), &
!                    KS_eigenvalue(i_state+n_occ,i_spin)
!         enddo
!       endif

! Build the delta_v_ia matrix appearing in the numberator
        do j_state = n_occ+1, n_states, 1
          j_state_1 = j_state - n_occ
          do i_state = 1, n_occ, 1
            delta_v_num(i_state, j_state_1) = &
              HF_ex_matr(i_state,j_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,j_state)
          enddo
        enddo

! transforming the delta_v_ia matrix from KS basis into HF basis
        call dgemm('N','N', n_occ, n_unocc, n_unocc, 1.d0, &   
                   delta_v_num, n_occ, trans_matr_unocc, n_unocc, 0.d0, &
                   tmp_delta_v_num, n_occ)

        call dgemm('T','N', n_occ, n_unocc, n_occ, 1.d0, &
                   trans_matr_occ, n_occ, tmp_delta_v_num, n_occ, 0.d0, &
                   delta_v_num, n_occ)
                   
        do j_state = 1, n_unocc, 1
          j_state_1 = j_state + n_occ
          do i_state = 1, n_occ, 1

          en_rse_full = en_rse_full + &
               delta_v_num(i_state,j_state)*delta_v_num(i_state,j_state) &
               /(HF_ev_occ(i_state)-HF_ev_unocc(j_state)) &
               *occ_numbers(i_state,i_spin,1)

        enddo
       enddo

      deallocate(HF_ham_occ)
      deallocate(HF_ham_unocc)
      deallocate(trans_matr_occ)
      deallocate(trans_matr_unocc)
      deallocate(HF_ev_occ)
      deallocate(HF_ev_unocc)
      deallocate(delta_v_num)
      deallocate(tmp_delta_v_num)
!  end loop over i_spin
      enddo

      if(myid.eq.0) then
        write(use_unit,*) 
        write(use_unit,'(2X,2A,f13.8, A)') &
           "| Renormalized single excitation contribution (full)", &
           " to the correlation energy:", &
           en_rse_full*hartree, "  eV"
        write(use_unit,*) 
      endif

      if(allocated(HF_ex_matr))then
        deallocate(HF_ex_matr)
      endif
      if(allocated(KS_xc_matr))then
        deallocate(KS_xc_matr)
      endif
      if(allocated(aux_xc_matr))then
        deallocate(aux_xc_matr)
      endif

      return

      end subroutine evaluate_renormalized_single_excitation_correction
!---------------------------------------------------------------------
!******
     subroutine evaluate_renormalized_single_excitation_correction_2 &
         ( n_KS_states,n_homo, &
           occ_numbers, ovlp_3KS, &
           KS_eigenvalue,KS_eigenvector, &
           xc_matr, en_rse_full &
          )

!  PURPOSE
!  Subroutine evaluate_renormalized_single_excitation_correction evaluates the single excitation
!  contribution to the second order correlation energy. This is zero for a HF reference
!  state, but non-zero for a KS reference state.
!
!  It is given by \sum_{ia} |<i|V_HFx - V_KSxc|a>|^2 / (E_i - E_a)
!  where V_HFx is the HF exchange operator,  V_KSxc is the KS
!  exchange-correlation potential operator, and the summation of i, a 
!  runs over the occupied states and unoccupied states respectively.
! USES

      use runtime_choices, only: flag_xc, hybrid_coeff, safe_minimum
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      integer :: n_KS_states
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  KS_eigenvector(n_basis,n_states,n_spin)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  xc_matr(n_basis,n_basis,n_spin)
      real*8  en_rse_full


! INPUTS
!  o n_KS_states -- number of KS states included in the calculations
!  o occ_numbers -- the occupation number 
!  o KS_eigenvector -- Single-particle KS eigenvector
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!  o xc_matr -- the matrix form the KS exchange-correlation potential operator
!            -- within the atomic orbital basis
!
! OUTPUT      
!  o en_rse_full -- a summation of a ladder-type single excitation contribution to the 
!              correlation energy, now also the non-diagonal terms are included.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




!  local variables

      real*8  tmp_numerator
      real*8, dimension(:,:), allocatable :: HF_ex_matr  
      real*8, dimension(:,:), allocatable :: KS_xc_matr  
      real*8, dimension(:,:), allocatable :: aux_xc_matr
      real*8, dimension(:,:), allocatable :: tmp_o3KS
      real*8, dimension(:), allocatable :: HF_ev_occ
      real*8, dimension(:), allocatable :: HF_ev_unocc
      real*8, dimension(:,:), allocatable :: HF_ham_occ
      real*8, dimension(:,:), allocatable :: HF_ham_unocc
      real*8, dimension(:,:), allocatable :: trans_matr_occ
      real*8, dimension(:,:), allocatable :: trans_matr_unocc
      real*8, dimension(:,:), allocatable :: delta_v_num
      real*8, dimension(:,:), allocatable :: tmp_delta_v_num

      integer n_homo(n_spin)
      integer n_homo_max, n_homo_min
      integer n_occ, n_unocc
      integer n_nonsingular

      integer mpierr
      real*8 ddot
!     counters

      integer :: i_state, i_state_loc, i_state_1
      integer :: j_state, j_state_loc, j_state_1
      integer :: k_state, k_state_loc
      integer :: i_spin
      integer :: n_p1



!     begin work

      if(flag_xc.eq.0) then
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A)') "HF reference state: single excitation contribution is zero"
        endif
        en_rse_full=0.d0
        return
      endif 

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the renormalized single excitation (full) correction ..."
      endif

!  determine the HOMO level
      do i_spin = 1, n_spin
       do i_state = 1, n_states
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
         endif
       enddo
      enddo

      n_homo_max=max(n_homo(1),n_homo(n_spin))
      n_homo_min=min(n_homo(1),n_homo(n_spin))

!    rearrange the index of ovlp_3KS to get a two-dimensioal, working matrix

      allocate(HF_ex_matr(n_states,n_states))
      allocate(KS_xc_matr(n_states,n_states))
      allocate(aux_xc_matr(n_states,n_basis))

!      n_homo_cur_max = 0
!      do i_spin = 1, n_spin
!        do n_p1 = 0, np1_o3KS-1
!          n_homo_cur = COUNT(own_dim1_o3ks(1:n_homo(i_spin))==n_p1)
!          n_homo_cur_max = MAX(n_homo_cur_max,n_homo_cur)
!        enddo
!      enddo

      allocate(tmp_o3KS(n_basbas,ndim1_o3KS))

      en_rse_full = 0.d0
      do i_spin = 1, n_spin

        if(n_homo(i_spin).eq.0) cycle

        HF_ex_matr(:,:) = 0.d0

        do k_state = 1, n_homo(i_spin), 1
          if(own_dim2_o3ks(k_state) /= myp2_o3ks) cycle
          k_state_loc = loc_dim2_o3ks(k_state)

          do n_p1 = 0, np1_o3KS-1

!            n_homo_cur = COUNT(own_dim1_o3ks(1:n_homo(i_spin))==n_p1)

            if(myp1_o3KS==n_p1) &
              tmp_o3KS(:,:) = ovlp_3KS(:,:,k_state_loc,i_spin)
            call mpi_bcast(tmp_o3KS,n_basbas*ndim1_o3ks,MPI_REAL8,n_p1,mpi_comm_rows_aux_2d,mpierr)

            do i_state = 1, n_states, 1
              if(own_dim1_o3ks(i_state) /= n_p1) cycle
              i_state_loc = loc_dim1_o3ks(i_state)

!              do j_state = n_homo(i_spin)+1, n_states, 1
              do j_state = 1, n_states, 1
                if(own_dim1_o3ks(j_state) /= myp1_o3ks) cycle
                j_state_loc = loc_dim1_o3ks(j_state)

                HF_ex_matr(i_state, j_state) = &
                 HF_ex_matr(i_state, j_state) - &
                    ddot( n_basbas, &
                         tmp_o3KS(1, i_state_loc), 1, &
                         ovlp_3KS(1, j_state_loc, k_state_loc, i_spin), 1 &
                        ) * &
                    occ_numbers(k_state,i_spin,1)*dble(n_spin)/2.d0
              enddo
            enddo
          enddo
        enddo

        call sync_matrix(HF_ex_matr, n_states, n_states)
    
        KS_xc_matr(:,:) = 0.d0
        aux_xc_matr = 0.0d0
        call dgemm('T', 'N', n_states, n_basis, &
                    n_basis, 1.0d0, &
                    KS_eigenvector(1,1,i_spin), n_basis, &
                    xc_matr(1,1,i_spin), n_basis, 0.d0, &
                    aux_xc_matr(1,1), &
                    n_states)

!     second multiplication between eigenvector and temporary XC matrix.
        call dgemm('N', 'N', n_states, n_states, &
                   n_basis, 1.0d0, &
                   aux_xc_matr(1,1), n_states, &
                   KS_eigenvector(1,1,i_spin), &
                   n_basis, 0.d0, KS_xc_matr, &
                   n_states)

!  Build the Hartree-Fock Hamiltonian for the occupied and unoccupied blocks separately         
        n_occ=n_homo(i_spin)
        n_unocc=n_states-n_occ
        allocate(HF_ham_occ(n_occ,n_occ))
        allocate(HF_ham_unocc(n_unocc,n_unocc))
        allocate(trans_matr_occ(n_occ,n_occ))
        allocate(trans_matr_unocc(n_unocc,n_unocc))
        allocate(HF_ev_occ(n_occ))
        allocate(HF_ev_unocc(n_unocc))
        allocate(delta_v_num(n_occ,n_unocc))
        allocate(tmp_delta_v_num(n_occ,n_unocc))

        HF_ham_occ(:,:) = 0.d0
        do i_state = 1, n_occ, 1
          do j_state = 1, n_occ, 1
            HF_ham_occ(j_state,i_state) = &
              HF_ex_matr(j_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(j_state,i_state)
          enddo
          HF_ham_occ(i_state,i_state) = HF_ham_occ(i_state,i_state) + & 
          KS_eigenvalue(i_state,i_spin)
        enddo

        call diagonalize_auxmat_lapack &
          (  n_occ, HF_ham_occ, safe_minimum, &
             -1.d9, n_nonsingular, HF_ev_occ, trans_matr_occ, &
             "RSE-HF-Occ" &
           )


            
        HF_ham_unocc(:,:) = 0.d0
        do i_state = n_occ+1, n_states, 1
          i_state_1 = i_state - n_occ
          do j_state = n_occ+1, n_states, 1
            j_state_1 = j_state - n_occ
            HF_ham_unocc(j_state_1,i_state_1) = &
              HF_ex_matr(j_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(j_state,i_state)
          enddo
          HF_ham_unocc(i_state_1,i_state_1) = HF_ham_unocc(i_state_1,i_state_1) + & 
          KS_eigenvalue(i_state,i_spin)

        enddo

        call diagonalize_auxmat_lapack &
          (  n_unocc, HF_ham_unocc, safe_minimum, &
            -1.d9, n_nonsingular, HF_ev_unocc, trans_matr_unocc, &
             "RSE-HF-Unocc" &
           )
       
! Build the delta_v_ia matrix appearing in the numberator
        do j_state = n_occ+1, n_states, 1
          j_state_1 = j_state - n_occ
          do i_state = 1, n_occ, 1
            delta_v_num(i_state, j_state_1) = &
              HF_ex_matr(i_state,j_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,j_state)
          enddo
        enddo

! transforming the delta_v_ia matrix from KS basis into HF basis
        call dgemm('N','N', n_occ, n_unocc, n_unocc, 1.d0, &   
                   delta_v_num, n_occ, trans_matr_unocc, n_unocc, 0.d0, &
                   tmp_delta_v_num, n_occ)

        call dgemm('T','N', n_occ, n_unocc, n_occ, 1.d0, &
                   trans_matr_occ, n_occ, tmp_delta_v_num, n_occ, 0.d0, &
                   delta_v_num, n_occ)
                   
        do j_state = 1, n_unocc, 1
          j_state_1 = j_state + n_occ
          do i_state = 1, n_occ, 1

          en_rse_full = en_rse_full + &
               delta_v_num(i_state,j_state)*delta_v_num(i_state,j_state) &
               /(HF_ev_occ(i_state)-HF_ev_unocc(j_state)) &
               *occ_numbers(i_state,i_spin,1)

        enddo
       enddo

       deallocate(HF_ham_occ)
       deallocate(HF_ham_unocc)
       deallocate(trans_matr_occ)
       deallocate(trans_matr_unocc)
       deallocate(HF_ev_occ)
       deallocate(HF_ev_unocc)
       deallocate(delta_v_num)
       deallocate(tmp_delta_v_num)

!end loop over i_spin
      enddo

      if(myid.eq.0) then
        write(use_unit,*) 
        write(use_unit,'(2X,2A,f13.8, A)') &
           "| Renormalized single excitation contribution (full)", &
           " to the correlation energy:", &
           en_rse_full*hartree, "  eV"
        write(use_unit,*) 
      endif

      if(allocated(HF_ex_matr))then
        deallocate(HF_ex_matr)
      endif
      if(allocated(KS_xc_matr))then
        deallocate(KS_xc_matr)
      endif
      if(allocated(aux_xc_matr))then
        deallocate(aux_xc_matr)
      endif
      if(allocated(tmp_o3KS))then
        deallocate(tmp_o3KS)
      endif

      return

      end subroutine evaluate_renormalized_single_excitation_correction_2
!---------------------------------------------------------------------
!******
