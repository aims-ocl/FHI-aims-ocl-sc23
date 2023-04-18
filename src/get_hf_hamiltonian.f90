!****s* FHI-aims/get_hf_hamiltonian
!  NAME
!   get_hf_hamiltonian
!  SYNOPSIS
      subroutine get_hf_hamiltonian  &
       ( number_of_loops, KS_eigenvalue, KS_eigenvector, &
         n_electrons, occ_numbers, &
         hamiltonian,fock_energy,en_xc ) 

!  PURPOSE
!   get_hf_hamiltonian adds the exchange contribution to the hamiltonian (H_0) 
!   obtained from integrate_hamiltonian matrix to get the Hartree-Fock 
!   Hamiltonian (i.e. Fock matrix).
!
!  USES
      use dimensions 
      use basis
      use runtime_choices
      use prodbas
      use hartree_fock
      use mpi_tasks
      use localorb_io
      use lc_wpbeh
      implicit none

!  ARGUMENTS 

      integer number_of_loops 
      real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
      real*8, dimension(n_basis,n_states,n_spin,n_k_points) ::  &
             KS_eigenvector
      real*8, dimension(n_states,n_spin,n_k_points) :: &
             occ_numbers
      real*8 :: n_electrons


      real*8  hamiltonian (n_basis*(n_basis+1)/2,n_spin)
      real*8  fock_energy
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

!     counters
      integer i_spin 
      integer i_state
      integer i_basis_1
      integer i_basis_2
      integer i_index
      integer :: info
      character*150 :: info_str
      character(*), parameter :: func = 'get_hf_hamiltonian'

!  start to work
       
      n_homo = 0
!      do i_state = 1, n_states
!       do i_spin = 1, n_spin
!         if(KS_eigenvalue(i_state, i_spin,1).le.chemical_potential) then
!           n_homo(i_spin)=i_state
!         endif
!       enddo
!      enddo

      do i_state = 1, n_states
       do i_spin = 1, n_spin
         if(occ_numbers(i_state,i_spin,1).gt.1.d-6) then
           n_homo(i_spin)=i_state
            
         endif
       enddo
      enddo

      if (n_spin .eq. 2 .and. use_hf_multiplicity) then
        n_homo (1) =  int((n_electrons+1)/2.d0) +  &
                     (hf_multiplicity - 1)/2
        n_homo (2) =  int((n_electrons)/2.d0) - &
                     (hf_multiplicity - 1)/2
        occ_numbers(1:n_homo(1),1,1) = 1.d0
        occ_numbers(n_homo(1)+1:n_states,1,1) = 0.d0
        occ_numbers(1:n_homo(2),2,1) = 1.d0
        occ_numbers(n_homo(2)+1:n_states,2,1) = 0.d0
      endif

      n_homo_max = max(n_homo(1),n_homo(n_spin))
      if (allocated(O_2bs1HF)) then
         if (n_homo_max > size(O_2bs1HF, 3)) then
            write(info_str, "(2X,A,': * O_2bs1HF needs to be reallocated')") &
            & trim(func)
            call localorb_info(info_str)
            deallocate(O_2bs1HF)
            allocate(O_2bs1HF(n_loc_prodbas, n_basis, n_homo_max, n_spin), &
            &        stat=info)
            call check_allocation(info, 'O_2bs1HF', func)
         end if
      end if
      

!     evaluate Fock matrix 
 	  if (sparse_o3fn) then
         select case (hf_version)
         case(HF_EIGEN)  ! Transformed overlap matrix
            if (RI_type == RI_LVL .or. RI_type == RI_LVL_2nd) then
               call evaluate_exchange_matr_LVL_eigen &
               &    (KS_eigenvector, number_of_loops, occ_numbers)
            else
               call aims_stop('Invalid RI_method', func)
            end if
         case(HF_DM)     ! Density matrix
            call aims_stop('LVL with density matrix not implemented', func)
         case default
            call aims_stop('Invalid hf_version', func)
         end select
      else
         select case (hf_version)
         case(HF_EIGEN)  ! Transformed overlap matrix
            call evaluate_exchange_matr_v1 &
            & (KS_eigenvector, number_of_loops, occ_numbers)
         case(HF_DM)     ! Density matrix
            call evaluate_exchange_matr_v0 &
            (KS_eigenvector, number_of_loops, occ_numbers)
         case default
            call aims_stop('Invalid hf_version', func)
         end select
      end if

      if (use_lc_wpbeh) then
      	call evaluate_2nd_fock_matr_for_lc_wpbeh(hamiltonian)
         call get_fock_energy_for_lc_wpbeh &
            (hybrid_coeff,KS_eigenvector,occ_numbers,fock_energy,en_xc)
      else
        if (packed_matrix_format /= PM_none) then
           call aims_stop('Invalid choice of packed_matrix_format for HF', func)
        end if
        do i_spin = 1, n_spin

          i_index = 0
          do i_basis_1 = 1, n_basis 
            do i_basis_2 = 1, i_basis_1 

              i_index = i_index + 1
              hamiltonian(i_index, i_spin) = &
                hamiltonian(i_index, i_spin) &
                - hybrid_coeff * fock_matr( i_basis_2, i_basis_1, i_spin)
!               write(use_unit,'(2I4,2f16.8)') i_basis_2, i_basis_1, fock_matr(i_basis_2,i_basis_1,i_spin), &
!                       hamiltonian(i_index,i_spin)
            enddo 
          enddo
    
        enddo
      
        call get_fock_energy &
            (hybrid_coeff,KS_eigenvector,occ_numbers,fock_energy,en_xc)
      end if
!! testing
!      write(use_unit,*) "fock_energy: ", fock_energy
!!      
      return
      end  subroutine get_hf_hamiltonian
!******
