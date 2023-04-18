!****s*  FHI-aims/evaluate_KS_xc_matrix
!  NAME
!     evaluate_KS_xc_matrix
!  SYNOPSIS

      subroutine evaluate_KS_xc_matrix &
           (KS_eigenvector, &
            xc_matr,  xc_KS_matr &
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
      real*8 xc_matr(n_basis,n_basis,n_spin)
      real*8 xc_KS_matr(n_states,n_states,n_spin)

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
      integer :: i_spin
      integer :: i_basis_1

!     begin work

!    rearrange the index of ovlp_3fn to get a two-dimensioal matrix

!     first multiplication between eigenvector and basis kinetic matrix.
      do i_spin = 1, n_spin
        aux_xc_matr = 0.0d0
        call dgemm('N', 'N', n_basis, n_states, &
              n_basis, 1.0d0, &
              xc_matr(:,:,i_spin),  n_basis, &
              KS_eigenvector(:,:,i_spin), &
              n_basis, 0.d0, aux_xc_matr, &
              n_basis)


!     second multiplication between eigenvector and temporary kinetic matrix.
        call dgemm('T', 'N', n_states, n_states, &
             n_basis, 1.0d0, &
             KS_eigenvector(:,:,i_spin), n_basis, &
             aux_xc_matr, &
             n_basis, 0.d0, xc_KS_matr(:,:,i_spin), &
             n_states)
      enddo


!      if(myid.eq.0) then
!        do i_spin = 1, n_spin, 1
!         do i_basis_1 = 1, n_basis, 1
!           write(use_unit,'(2I4,f16.8)')i_spin,i_basis_1, &
!                xc_KS_matr(i_basis_1,i_basis_1,i_spin)
!           write(use_unit,*)
!         enddo
!        enddo
!      endif


      end subroutine evaluate_KS_xc_matrix
!---------------------------------------------------------------------
!******
