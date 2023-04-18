!****s* FHI-aims/get_screx_hamiltonian
!  NAME
!   get_screx_hamiltonian
!  SYNOPSIS

      subroutine get_screx_hamiltonian &
        ( number_of_loops, KS_eigenvalue, KS_eigenvector, &
          n_electrons, occ_numbers, &
          hamiltonian, en_xc )
!  PURPOSE
!    get_screx_hamiltonian adds the screened exchange term to the hamiltonian (H_0)
!    obtained from integrate_hamiltonian matrix
!    It evaluates the screened (RPA screening at zero frequency) exchange matrix 
!    and if doing COHSEX calculation, (use_cohsex = .true.), it evaluates the
!    COHSEX matrix
!
!  USES

      use dimensions
      use basis
      use runtime_choices
      use prodbas
      use hartree_fock

      implicit none

!  ARGUMENTS  

      integer number_of_loops
      real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
      real*8, dimension(n_basis,n_states,n_spin,n_k_points) :: &
              KS_eigenvector
      real*8, dimension(n_states,n_spin,n_k_points) :: &
              occ_numbers
      real*8 :: n_electrons

      real*8  hamiltonian (n_basis*(n_basis+1)/2,n_spin)
      real*8  en_xc

!  INPUTS
!  o  number_of_loops -- integer number, the current self-consistent loop
!  o  KS_eigenvalue -- real array,
!        the eigenvalues of the single-particle calculation. For DFT calculation,
!        this is the KS eigenvalue, but for HF calculation, this is then the HF
!        eigenvalue
!  o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
!  o  occ_numbers -- real array,
!            the occupation number of the electrons for each eigenstate and each spin
!  o  n_electrons -- real number
!            the total number of electrons in the system
!  OUTPUTS
!  o  hamiltonian -- real number, HF hamiltonian in this case (i.e., the Fock matrix)
!  o  en_xc -- real number, the exchange-correlation energy, differs form fock_energy
!          in cases of hybrid functional calculations
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

!     local variables
      real*8 alpha, fock_energy
!   now defined in the module hartree_fock
!      real*8, dimension(:,:,:,:), allocatable ::  transformed_ovlp3fn

!     counters
      integer i_spin
      integer i_state
      integer i_basis_1
      integer i_basis_2
      integer i_index

!  start to work

       screx_matr(:,:,:) = 0.d0
!      if(.not. allocated(transformed_ovlp3fn)) then
!        allocate(transformed_ovlp3fn(
!     +           n_loc_prodbas,n_basis,n_states,n_spin))
!        transformed_ovlp3fn(:,:,:,:) = 0.d0
!      endif

      n_homo = n_states

      do i_state = 1, n_states
       do i_spin = 1, n_spin
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state

         endif
       enddo
      enddo

!      if (n_spin .eq. 2 .and. use_hf_multiplicity) then
!        n_homo (1) =  int((n_electrons+1)/2.d0) +
!     +                (hf_multiplicity - 1)/2
!        n_homo (2) =  int((n_electrons)/2.d0) -
!     +                (hf_multiplicity - 1)/2
!        occ_numbers(1:n_homo(1),1,1) = 1.d0
!        occ_numbers(n_homo(1)+1:n_states,1,1) = 0.d0
!        occ_numbers(1:n_homo(2),2,1) = 1.d0
!        occ_numbers(n_homo(2)+1:n_states,2,1) = 0.d0
!      endif

      n_homo_max = max(n_homo(1),n_homo(n_spin))

!     evaluate Fock matrix
      call evaluate_exchange_matr_v0(KS_eigenvector, &
                                       number_of_loops, &
                                      occ_numbers)

!     evaluate O_2bs1HF(1:n_basbas, 1:n_basis, 1:n_states)
      call evaluate_O_2bs1HF ( KS_eigenvector )

!     evaluate Fock matrix
      call evaluate_screx_matrix(KS_eigenvalue, &
                                occ_numbers,KS_eigenvector, &
                                O_2bs1HF)

      do i_spin = 1, n_spin

        i_index = 0
        do i_basis_1 = 1, n_basis
          do i_basis_2 = 1, i_basis_1

            i_index = i_index + 1
            hamiltonian(i_index, i_spin) = &
               hamiltonian(i_index, i_spin) &
               -  fock_matr( i_basis_2, i_basis_1, i_spin) &
               -  screx_matr( i_basis_2, i_basis_1, i_spin)
!            if(i_basis_1.eq.i_basis_2) then
!               write(use_unit,'(I4,3f18.8)')i_basis_1, hamiltonian(i_index, i_spin), fock_matr( i_basis_2, i_basis_1, i_spin), &
!                         screx_matr( i_basis_2, i_basis_1, i_spin)
!            endif
          enddo
        enddo

      enddo

      fock_matr(:,:,:) = fock_matr(:,:,:) + &
                2.d0 * screx_matr(:,:,:)

      call get_fock_energy &
          (hybrid_coeff,KS_eigenvector,occ_numbers,fock_energy,en_xc)

!      if(allocated(transformed_ovlp3fn)) then
!        deallocate(transformed_ovlp3fn)
!      endif

      return
      end  subroutine get_screx_hamiltonian
!******
