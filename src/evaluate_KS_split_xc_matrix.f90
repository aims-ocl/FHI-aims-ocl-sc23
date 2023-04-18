!****s*  FHI-aims/evaluate_KS_split_xc_matrix
!  NAME
!     evaluate_KS_split_xc_matrix
!  SYNOPSIS

      subroutine evaluate_KS_split_xc_matrix &
           (KS_eigenvector, &
            x_matr,c_matr, &
            x_KS_array, c_KS_array &
           )

!  PURPOSE
!  Subroutine evaluate_KS_kinetic_matrix evaluates the exchange-correlation
!  energy  matrix within 2 KS orbitals, produced by matrix multiplication of
!  the KS eigenvectors  and exchange-correlation energy
!  matrix within local basis.
!
!  Vxc_ij = \sum_mm c_mi c_nj  Vxc_mn
!
!  USES
      use dimensions
      use mpi_tasks

      implicit none

!  ARGUMENTS

      real*8 KS_eigenvector(n_basis,n_states,n_spin)
      real*8 x_matr(n_basis,n_basis,n_spin)
      real*8 c_matr(n_basis,n_basis,n_spin)
      real*8 x_KS_array(n_states,n_spin)
      real*8 c_KS_array(n_states,n_spin)

! INPUTS
!  o KS_eigenvector -- real array, the eigenvector from the single-particle calculation
!  o xc_matr -- the matrix elements of the exchange correlation potential
!    within the regular basis functions 
!          
! OUTPUT 
!  o xc_KS_matr -- the matrix elements of the exchange correlation potential
!    witin KS orbitals
!
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

!     auxiliary matrices for Level 3 Blas matrix multiplications


      real*8 aux_xc_matr(n_basis,n_states)


!     counters

      integer :: i_state

      integer :: i_basis_1
      integer :: i_spin

!     begin work

!    rearrange the index of ovlp_3fn to get a two-dimensioal matrix

!     first multiplication between eigenvector and basis kinetic matrix.
      x_KS_array(:,:)=0.d0
      c_KS_array(:,:)=0.d0
      do i_spin = 1, n_spin

        call dgemm('N', 'N', n_basis, n_states, &
              n_basis, 1.0d0, &
              x_matr(:,:,i_spin),  n_basis, &
              KS_eigenvector(:,:,i_spin), &
              n_basis, 0.d0, aux_xc_matr, &
              n_basis)

        do i_state = 1, n_states, 1
          do i_basis_1 =1, n_basis, 1
            x_KS_array(i_state,i_spin) = x_KS_array(i_state,i_spin) + &
            KS_eigenvector(i_basis_1,i_state,i_spin) * &
            aux_xc_matr(i_basis_1,i_state)
          enddo
        enddo

        call dgemm('N', 'N', n_basis, n_states, &
              n_basis, 1.0d0, &
              c_matr(:,:,i_spin),  n_basis, &
              KS_eigenvector(:,:,i_spin), &
              n_basis, 0.d0, aux_xc_matr, &
              n_basis)

        do i_state = 1, n_states, 1
          do i_basis_1 =1, n_basis, 1
            c_KS_array(i_state,i_spin) = c_KS_array(i_state,i_spin) + &
            KS_eigenvector(i_basis_1,i_state,i_spin) * &
            aux_xc_matr(i_basis_1,i_state)
          enddo
        enddo

      enddo
!      if(myid.eq.0) then
!        do i_spin = 1, n_spin, 1
!         do i_basis_1 = 1, n_basis, 1
!           write(use_unit,'(2I4,2f16.8)')i_spin,i_basis_1, x_KS_array(i_basis_1,i_spin), &
!                 c_KS_array(i_basis_1,i_spin)
!           write(use_unit,*)
!         enddo
!        enddo
!      endif

      return
      end subroutine evaluate_KS_split_xc_matrix
!---------------------------------------------------------------------
!******
