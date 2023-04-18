!****s* FHI-aims/f_external_constraint
!  NAME
!    f_external_constraint
!  SYNOPSIS 
      subroutine f_external_constraint(n_arguments, arguments, &
           f_value,ifn, &
           hamiltonian, &
           overlap_matrix, KS_eigenvalue, KS_eigenvector, occ_numbers, &
           n_electrons, chemical_potential, constraint_proj, &
           electrons_in_region, &
           current_spin &
           )
!  PURPOSE
!    Provides the target function for the constrained DFT calculations as the
!    l2-norm of the difference between the actual number of electrons in the regions and
!    given target number of electrons in the regions. This routine is passed externally 
!    to the BFGS-optimizer.
!  USES
      use constraint
      use dimensions
      use localorb_io
      use constants

      implicit none
!  ARGUMENTS
      real*8 :: f_value
      integer :: n_arguments
      real*8, dimension(3*n_arguments) :: arguments
      integer :: ifn
      real*8 hamiltonian( n_basis*(n_basis+1)/2, n_spin )
      real*8 overlap_matrix( n_basis*(n_basis+1)/2 )
      real*8, dimension(n_states, n_spin) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential
      real*8, dimension(n_states, n_region, n_spin) :: constraint_proj
      real*8, dimension(n_region, n_spin) :: electrons_in_region
      integer :: current_spin
!  INPUTS
!    o n_arguments -- number of arguments to f
!    o arguments -- the array of arguments to f
!    o ifn -- function call counter
!    o hamiltonian -- the Hamilton matrix
!    o overlap_matrix -- the overlap matrix
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!    o occ_numbers -- occupation numbers
!    o n_electrons -- number of electrons in the spin channel
!    o chemical_potential -- the chemical potential
!    o constraint_proj -- projection to the constrained regions
!    o electrons_in_region -- number of electrons in the regions
!    o current_spin -- current spin channel
!  OUTPUT
!    o f_value -- value of the target function f
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


!     locals

      character*100 :: info_str

      real*8 :: avg_zero
      integer :: i_region, i_hamiltonian
      logical :: t_out = .false.

      ! dummy for output purposes only, constrained DFT does not support periodic boundary conditions (yet)
      integer :: i_k_point = 1

      real*8 :: avg_of_arguments

!     begin work

      if (constraint_debug) then
        call localorb_info('',use_unit,'(A)')
        write(info_str,'(2X,A,I5,A)') &
        "Constrained DFT: Eigenvector update number ", &
        ifn,"."
        call localorb_info(info_str,use_unit,'(A)')
        call localorb_info('',use_unit,'(A)')
      end if

      avg_of_arguments = 0.0d0

      ! FIXME EITHER I DO NOT UNDERSTAND THIS AVERAGING OR IT IS SIMPLY INCORRECT.

      do i_region = 1, n_arguments, 1
         avg_of_arguments = avg_of_arguments + arguments(i_region) * &
                            constraint_electrons(i_region,current_spin)
      enddo
      avg_of_arguments = avg_of_arguments / n_electrons

      do i_region = 1, n_active_regions - 1, 1
         constraint_potential(i_region,current_spin) = &
              arguments(i_region) - avg_of_arguments
      enddo

      constraint_potential(n_active_regions,current_spin) = &
           - avg_of_arguments

      !  Restore original hamiltonian for current spin component
      do i_hamiltonian=1, n_basis*(n_basis+1)/2, 1
            hamiltonian_work(i_hamiltonian,current_spin) = &
                 hamiltonian(i_hamiltonian,current_spin)
      enddo

      ! ... and add the current constraint potentials.
      call add_constraint_potentials_v2 &
           ( current_spin, overlap_matrix, hamiltonian_work &
           )

      if (constraint_debug) then
                   write(info_str,'(2X,A)') &
                        ''
                   write(info_str,'(2X,A)') &
                        'Region                    Constraint potential'
                   call localorb_info(info_str,use_unit,'(A)')
                   do i_region = 1, n_active_regions, 1
                     write(info_str,'(2X,I5,20X,F20.15)') &
                     i_region, &
                     constraint_potential(i_region,current_spin)*hartree
                     call localorb_info(info_str,use_unit,'(A)')
                   enddo
                   write(info_str,'(2X,A)') &
                        ''
      end if

      ! solve for new KS eigenstates and -energies ...
      ! note that i_k_point is here only a dummy, and used only for output purposes
      ! in improve_real_eigenfunctions
      call improve_real_eigenfunctions &
           ( overlap_matrix, hamiltonian_work(:,current_spin), &
           t_out, &
           KS_eigenvalue(:,current_spin), &
           KS_eigenvector(:,:,current_spin), i_k_point &
           )

      ! ... and obtain the updated global Fermi level and
      ! and global occupation numbers

      ! Notice n_electrons is here the number of electrons in the current spin channel
      call get_occupation_numbers_v2 &
           ( KS_eigenvalue(:,current_spin), n_electrons, &
           t_out, &
           occ_numbers(:,current_spin), chemical_potential &
           )

!test
!      write(use_unit,*) "occ_numbers: "
!      write(use_unit,*) occ_numbers
!      write(use_unit,*)
!test end

      ! update constraint projectors and count number of electrons in each region, each spin
      call get_electrons_per_region_v2 &
           ( current_spin, KS_eigenvector, overlap_matrix, occ_numbers, &
           constraint_proj, &
           electrons_in_region &
           )

      f_value = 0.0d0
      do i_region = 1, n_active_regions, 1
         f_value = f_value + &
              ( electrons_in_region(i_region,current_spin) - &
              constraint_electrons(i_region,current_spin) )**2.d0
      enddo

      f_value = dsqrt(f_value)

      end subroutine f_external_constraint
!******
