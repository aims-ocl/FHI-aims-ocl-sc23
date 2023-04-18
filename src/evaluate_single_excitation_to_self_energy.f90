!****s*  FHI-aims/evaluate_single_excitation_to_self_energy
!  NAME
!    evaluate_single_excitation_to_self_energy
!  SYNOPSIS

     subroutine evaluate_single_excitation_to_self_energy &
         ( n_low_state,n_high_state, &
           n_homo,n_freq, omega, &
           occ_numbers, ovlp_3KS, &
           KS_eigenvalue,KS_eigenvector, &
           xc_matr, self_energy_freq &
          )

!  PURPOSE
!  Subroutine evaluate_single_excitation_to_self_energy evaluates the single excitation
!  contribution to the second order correlation energy. This is zero for a HF reference
!  state, but non-zero for a KS reference state.
!
!  It is given by \sum_{ia} |<i|V_HFx - V_KSxc|a>|^2 / (E_i - E_a)
!  where V_HFx is the HF exchange operator,  V_KSxc is the KS
!  exchange-correlation potential operator, and the summation of i, a 
!  runs over the occupied states and unoccupied states respectively.
! USES

      use runtime_choices, only: flag_xc, hybrid_coeff
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      integer :: n_low_state
      integer :: n_high_state
      integer :: n_freq
      integer :: n_homo(n_spin)
      real*8  omega(n_freq)
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  KS_eigenvector(n_basis,n_states,n_spin)
      real*8  ovlp_3KS(n_loc_prodbas,n_states,n_high_state,n_spin)
      real*8  xc_matr(n_basis,n_basis,n_spin)
      real*8  self_energy_freq(n_low_state:n_high_state,n_freq,n_spin)


! INPUTS
!  o n_low_state -- the lowest KS state included in the self-energy calculation, 
!                    =1 in the all-electron case
!  o n_high_state -- the highest KS state included in the self-energy calculation,
!                    >= n_homo
!  o n_homo       -- HOMO level
!  o n_freq       -- number of frequency points
!  o omega      -- the Gauss-Legendre frequency grid for the self-energy
!  o occ_numbers -- the occupation number 
!  o KS_eigenvector -- Single-particle KS eigenvector
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!  o xc_matr -- the matrix form the KS exchange-correlation potential operator
!            -- within the atomic orbital basis
!  o self_energy_freq -- self energy on the imaginary frequency axis before adding 
!  o               single excitation correction
!
! OUTPUT      
!  o self_energy_freq -- self energy on the imaginary frequency axis after adding 
!  o               single excitation correction
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
      integer n_homo_max, n_homo_min
      real*8 ddot
      real*8 :: qpe_se(n_low_state:n_high_state,n_spin)
      real*8 :: se_energy
!     counters

      integer :: i_state
      integer :: j_state
      integer :: k_state
      integer :: i_spin
      integer :: i_freq


!     begin work

      if(flag_xc.eq.0) then
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A)') "HF reference state: single excitation contribution is zero"
        endif
        return
      endif 

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the single excitation contribution to the self-energy ..."
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

      qpe_se(:,:) = 0.d0
      se_energy = 0.d0
      do i_spin = 1, n_spin
        HF_ex_matr(:,:) = 0.d0
        do j_state = 1, n_states, 1

         do i_state = 1, n_states, 1
          do k_state = 1, n_homo(i_spin), 1

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

      call sync_matrix( &
             HF_ex_matr, n_states, n_states)
  
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

      do i_freq = 1, n_freq, 1

        do j_state = 1, n_states, 1
!        do j_state = n_homo(i_spin)+1, n_states, 1
          do i_state = n_low_state, n_high_state, 1
!          do i_state = 1, n_homo(i_spin), 1
           
!          if( (i_state .gt. n_homo(i_spin) .and. j_state .gt. n_homo(i_spin)) .or. &
!              (i_state .le. n_homo(i_spin) .and. j_state .le. n_homo(i_spin)) ) cycle
!           if( (i_state .le. n_homo(i_spin) .and. j_state .le. n_homo(i_spin)) ) cycle
!           if( i_state .le. n_homo(i_spin) ) cycle

          tmp_numerator = (HF_ex_matr(i_state,j_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,j_state))*  &
          (HF_ex_matr(i_state,j_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,j_state))

         if(abs(KS_eigenvalue(i_state,i_spin) - KS_eigenvalue(j_state,i_spin)) .lt. 1.e-6) cycle
!          self_energy_freq(i_state,i_freq,i_spin) = &
!               self_energy_freq(i_state,i_freq,i_spin)  + &
!               tmp_numerator/( (0.d0,1.d0)*omega(i_freq) - KS_eigenvalue(j_state,i_spin) )
!               tmp_numerator/( KS_eigenvalue(i_state,i_spin) - KS_eigenvalue(j_state,i_spin) )
!               tmp_numerator/( (0.d0,1.d0)*omega(i_freq) - KS_eigenvalue(j_state,i_spin) + &
!                HF_ex_matr(i_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,i_state) ) 

          if(i_freq .eq. 1) then
!             qpe_se(i_state,i_spin) = qpe_se(i_state,i_spin)  + &
!               tmp_numerator/( KS_eigenvalue(i_state,i_spin) - KS_eigenvalue(j_state,i_spin) + &
!                HF_ex_matr(i_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,i_state) ) 
             if(abs(KS_eigenvalue(i_state,i_spin) - KS_eigenvalue(j_state,i_spin)) .lt. 1.e-6) cycle
             qpe_se(i_state,i_spin) = qpe_se(i_state,i_spin)  + &
               tmp_numerator/( KS_eigenvalue(i_state,i_spin) - KS_eigenvalue(j_state,i_spin)) 
             se_energy = se_energy + &
               tmp_numerator/( KS_eigenvalue(i_state,i_spin) - KS_eigenvalue(j_state,i_spin))

!               write(use_unit,'(3I6,2f12.6)')i_spin, i_state, j_state, &
!                  tmp_numerator/(KS_eigenvalue(i_state,i_spin)- KS_eigenvalue(j_state,i_spin)), &
!                  se_energy
!              if(i_state .eq. 7) then
!               write(use_unit,'(A, I4,6f12.6)')'HOMO', j_state,tmp_numerator,KS_eigenvalue(i_state,i_spin), &
!                KS_eigenvalue(j_state,i_spin), HF_ex_matr(i_state,i_state), &
!                KS_xc_matr(i_state,i_state), qpe_se(i_state,i_spin)
!              endif
!              if(i_state .eq. 8) then
!               write(use_unit,'(A, I4,6f12.6)')'LUMO', j_state,tmp_numerator,KS_eigenvalue(i_state,i_spin), &
!                KS_eigenvalue(j_state,i_spin), HF_ex_matr(i_state,i_state), &
!                KS_xc_matr(i_state,i_state), qpe_se(i_state,i_spin)
!              endif
          endif

        enddo
       enddo

! end of loop over i_freq
      enddo

!  end loop over i_spin
      enddo

      se_energy = 0.d0
      if(myid.eq.0) then
        write(use_unit,*) 
        write(use_unit,'(2X,A)') &
           "| Single excitation correction to quasiparticle energy level: "
        do i_spin = 1, n_spin, 1
          do i_state = n_low_state, n_high_state, 1
            se_energy = se_energy + qpe_se(i_state,i_spin)
            write(use_unit, '(2I6,2f13.4)') i_spin, i_state, &
               qpe_se(i_state,i_spin)*hartree, se_energy*hartree
         enddo
        enddo
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

      end subroutine evaluate_single_excitation_to_self_energy
!---------------------------------------------------------------------
!******
     subroutine evaluate_single_excitation_to_self_energy_2 &
         ( n_low_state,n_high_state, &
           n_homo,n_freq, omega, &
           occ_numbers, ovlp_3KS, &
           KS_eigenvalue,KS_eigenvector, &
           xc_matr, self_energy_freq &
          )

!  PURPOSE
!  Subroutine evaluate_single_excitation_to_self_energy evaluates the single excitation
!  contribution to the second order correlation energy. This is zero for a HF reference
!  state, but non-zero for a KS reference state.
!
!  It is given by \sum_{ia} |<i|V_HFx - V_KSxc|a>|^2 / (E_i - E_a)
!  where V_HFx is the HF exchange operator,  V_KSxc is the KS
!  exchange-correlation potential operator, and the summation of i, a 
!  runs over the occupied states and unoccupied states respectively.
! USES

      use runtime_choices, only: flag_xc, hybrid_coeff
      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use synchronize_mpi
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      integer :: n_low_state
      integer :: n_high_state
      integer :: n_homo(n_spin)
      integer :: n_freq
      real*8  omega(n_freq)
      real*8  occ_numbers(n_states,n_spin,n_k_points)
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  KS_eigenvector(n_basis,n_states,n_spin)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  xc_matr(n_basis,n_basis,n_spin)
      real*8  self_energy_freq(n_low_state:n_high_state,n_freq,n_spin)


! INPUTS
!  o n_low_state -- the lowest KS state included in the self-energy calculation, 
!                    =1 in the all-electron case
!  o n_high_state -- the highest KS state included in the self-energy calculation,
!                    >= n_homo
!  o n_homo   --    HOMO level
!  o n_freq   --    number of frequency points
!  o omega   --    the Gauss-Legendre frequency grid for the self-energy
!  o occ_numbers -- the occupation number 
!  o KS_eigenvector -- Single-particle KS eigenvector
!  o ovlp_3KS -- the transformed 3-function overlap integral involving
!                   two KS states and one auxiliary basis
!  o xc_matr -- the matrix form the KS exchange-correlation potential operator
!            -- within the atomic orbital basis
!  o self_energy_freq -- self energy on the imaginary frequency axis before adding 
!  o               single excitation correction
!
! OUTPUT      
!  o self_energy_freq -- self energy on the imaginary frequency axis after adding 
!  o               single excitation correction
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
      integer n_homo_max, n_homo_min
      integer n_homo_cur, n_homo_cur_max
      integer mpierr
      real*8 ddot
!     counters

      integer :: i_state, i_state_loc
      integer :: j_state, j_state_loc
      integer :: k_state, k_state_loc
      integer :: i_spin
      integer :: i_freq
      integer :: n_p1



!     begin work

      if(flag_xc.eq.0) then
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A)') "HF reference state: single excitation contribution is zero"
        endif
        return
      endif 

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,'(2X,A)')"Starts to calculate the single excitation to_self_energy ..."
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

      allocate(HF_ex_matr(n_homo_max,n_states))
      allocate(KS_xc_matr(n_homo_max,n_states))
      allocate(aux_xc_matr(n_homo_max,n_basis))

      n_homo_cur_max = 0
      do i_spin = 1, n_spin
        do n_p1 = 0, np1_o3KS-1
          n_homo_cur = COUNT(own_dim1_o3ks(1:n_homo(i_spin))==n_p1)
          n_homo_cur_max = MAX(n_homo_cur_max,n_homo_cur)
        enddo
      enddo

      allocate(tmp_o3KS(n_basbas,n_homo_cur_max))

      do i_spin = 1, n_spin
        HF_ex_matr(:,:) = 0.d0

        do k_state = 1, n_homo(i_spin), 1
          if(own_dim2_o3ks(k_state) /= myp2_o3ks) cycle
          k_state_loc = loc_dim2_o3ks(k_state)

          do n_p1 = 0, np1_o3KS-1

            n_homo_cur = COUNT(own_dim1_o3ks(1:n_homo(i_spin))==n_p1)

            if(myp1_o3KS==n_p1) &
              tmp_o3KS(:,1:n_homo_cur) = ovlp_3KS(:,1:n_homo_cur,k_state_loc,i_spin)
            call mpi_bcast(tmp_o3KS,n_basbas*n_homo_cur,MPI_REAL8,n_p1,mpi_comm_rows_aux_2d,mpierr)

            do i_state = 1, n_homo(i_spin), 1
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

        call sync_matrix(HF_ex_matr, n_homo_max, n_states)
    
        KS_xc_matr(:,:) = 0.d0
        aux_xc_matr = 0.0d0
        call dgemm('T', 'N', n_homo(i_spin), n_basis, &
                  n_basis, 1.0d0, &
                  KS_eigenvector(1,1,i_spin), n_basis, &
                  xc_matr(1,1,i_spin), n_basis, 0.d0, &
                  aux_xc_matr(1,1), &
                  n_homo_max)
         
!     second multiplication between eigenvector and temporary XC matrix.
        call dgemm('N', 'N', n_homo(i_spin), n_states, &
                   n_basis, 1.0d0, &
                   aux_xc_matr(1,1), n_homo_max, &
                   KS_eigenvector(1,1,i_spin), &
                   n_basis, 0.d0, KS_xc_matr, &
                   n_homo_max)




      do i_freq = 1, n_freq, 1

        do j_state = 1, n_states, 1
          do i_state = n_low_state, n_high_state, 1
           
          if( (i_state .gt. n_homo(i_spin) .and. j_state .gt. n_homo(i_spin)) .or. &
              (i_state .lt. n_homo(i_spin) .and. j_state .lt. n_homo(i_spin)) ) cycle

!          write(use_unit,'(3I4,2f13.6)') i_spin, i_state, j_state, &
!                 HF_ex_matr(i_state,j_state), &
!                 KS_xc_matr(i_state,j_state)
!                 HF_ex_matr(i_state,i_state),KS_xc_matr(i_state,i_state), &
!                 KS_eigenvalue(i_state,i_spin)-KS_eigenvalue(j_state,i_spin), &
!                 en_single_excit*hartree

          tmp_numerator = (HF_ex_matr(i_state,j_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,j_state))*  &
          (HF_ex_matr(i_state,j_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,j_state))

          self_energy_freq(i_state,i_freq,i_spin) = &
               self_energy_freq(i_state,i_freq,i_spin)  + &
               tmp_numerator/( (0.d0,1.d0)*omega(i_freq) - KS_eigenvalue(j_state,i_spin) + &
                HF_ex_matr(i_state,i_state)*(1.0-hybrid_coeff)-KS_xc_matr(i_state,i_state) ) &
               *occ_numbers(i_state,i_spin,1)
 
        enddo
       enddo
! end of loop over i_freq
      enddo
!  end loop over i_spin
      enddo


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

      end subroutine evaluate_single_excitation_to_self_energy_2
!---------------------------------------------------------------------
!******
