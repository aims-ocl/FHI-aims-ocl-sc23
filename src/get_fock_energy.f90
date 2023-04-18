!****s* FHI-aims/get_fock_energy
!  NAME
!   get_fock_hamiltonian
!  SYNOPSIS
      subroutine get_fock_energy &
           (alpha,KS_eigenvector,occ_numbers,fock_energy,energy_xc)

!  PURPOSE
!  Subroutine get_fock_energy evaluates the exact exchange energy
!  for Hatree-Fock (HF) and HF-based calculations.
!
!  USES

      use dimensions
      use hartree_fock

      implicit none

!  imported variables
!  now defind in modle hartree_fock
!   alpha = 1: Hartree-Fock; 0.25: PBE0

!  ARGUMENTS

      real*8 alpha
      real*8, dimension(n_basis,n_states,n_spin,n_k_points) ::  &
              KS_eigenvector
      real*8, dimension(n_states,n_spin,n_k_points) :: &
              occ_numbers

      real*8  fock_energy 
      real*8  energy_xc
!  INPUTS
!  o  alpha -- the mixing parameter for HF-based calculations, 1.0 for HF,
!         0.25 for PBE0, etc.
!  o  KS_eigenvector -- real array,
!         the eigenvector of the single-particle calculation
!  o  occ_numbers -- real array,
!         the occupation number of the electrons for each eigenstate and each spin
!  OUTPUTS
!  o  fock_energy -- real number, here the exact exchange energy
!  o  en_xc -- real number, the exchange-correlation energy, differs form fock_energy
!          in cases of hybrid functional calculations
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
      real*8 energy_fock

!     auxiliary matrices for Level 3 Blas matrix multiplications
      real*8   aux_fock_matr(n_basis,n_basis)


!     counters

      integer :: i_state
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_index
      integer :: i_spin

!     begin work

!     first multiplication between eigenvector and Fock matrix

      energy_fock=0.d0
      do i_spin = 1, n_spin
        do i_state = 1, n_homo(i_spin)

          do i_basis_1 = 1, n_basis
             do i_basis_2 = 1, n_basis
               energy_fock = energy_fock + &
                fock_matr(i_basis_1,i_basis_2,i_spin) * &
                KS_eigenvector(i_basis_1, i_state, i_spin, 1) * &
                KS_eigenvector(i_basis_2, i_state, i_spin, 1) * &
                occ_numbers(i_state,i_spin,1)
           enddo
          enddo

        enddo
      enddo

      fock_energy = -0.5d0*energy_fock
      energy_xc = energy_xc + 0.5d0*alpha*energy_fock


      end subroutine get_fock_energy
!---------------------------------------------------------------------
!******
