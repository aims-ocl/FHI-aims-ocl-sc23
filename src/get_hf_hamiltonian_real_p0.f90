!****s* FHI-aims/get_hf_hamiltonian_real_p0
!  NAME
!   get_hf_hamiltonian_real_p0
!  SYNOPSIS
      subroutine get_hf_hamiltonian_real_p0  &
       ( exchange_matr_real, hamiltonian_w, i_k) 

!  PURPOSE
!   get_hf_hamiltonian_p0 adds the exchange contribution to the hamiltonian (H_0) 
!   obtained from integrate_hamiltonian matrix to get the Hartree-Fock 
!   Hamiltonian (i.e. Fock matrix).
!
!  USES
      use dimensions 
      use basis
      use runtime_choices
      use prodbas
      use hartree_fock_p0
      use mpi_tasks

      implicit none

!  ARGUMENTS 

      real*8, dimension(n_basis,n_basis,n_k_points_task,n_spin) ::  &
             exchange_matr_real
      real*8, dimension(n_basis*(n_basis+1)/2,n_spin) :: hamiltonian_w
      integer :: i_k

!  INPUTS
!  o  exchange_matr_real -- real array,
!            the exchange matrix in case of real eigenvectors
!  o i_k -- the number of k-point on given task, for which the hamiltonian
!            matrix is updated
!  OUTPUTS
!  o  hamiltonian_w -- real array, HF hamiltonian in this case (i.e., the Fock matrix)
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

!     counters
      integer i_spin 
      integer i_state
      integer i_basis_1
      integer i_basis_2
      integer i_index

!  start to work
      
      do i_spin = 1, n_spin

         i_index = 0
         do i_basis_1 = 1, n_basis
            do i_basis_2 = 1, i_basis_1
               
               i_index = i_index + 1
               
               if (use_lc_wpbeh) then
                 hamiltonian_w(i_index, i_spin) = &
                    hamiltonian_w(i_index, i_spin) &
                    - exchange_matr_real( i_basis_2, i_basis_1, i_k, i_spin)
               else
                 hamiltonian_w(i_index, i_spin) = &
                    hamiltonian_w(i_index, i_spin) &
                    - hybrid_coeff * exchange_matr_real( i_basis_2, i_basis_1, i_k, i_spin)
               endif
            enddo
         enddo
         
      enddo

      return
    end  subroutine get_hf_hamiltonian_real_p0
!******
