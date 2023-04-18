!****s* FHI-aims/unshift_eigenvalues
!  NAME
!    unshift_eigenvalues
!  SYNOPSIS

subroutine unshift_eigenvalues(n_electrons, KS_eigenvalue,KS_eigenvector, occ_numbers, chemical_potential, chemical_potential_spin)

!  PURPOSE
!  If spin constrain is used, this routine unshifts the eigenvalues after SC is achieved.
!
!  USES

      use dimensions
      use geometry
      use basis
      use runtime_choices
      use constraint
      use localorb_io
      use constants
      implicit none

!  ARGUMENTS

      real*8, dimension(n_states, n_spin,n_k_points) :: KS_eigenvalue
      real*8, dimension(n_basis, n_states, n_spin, n_k_points_task) ::  KS_eigenvector
      real*8, dimension(n_states, n_spin,n_k_points) :: occ_numbers
      real*8 :: n_electrons
      real*8 :: chemical_potential  
      real*8, dimension(n_spin) :: chemical_potential_spin

!  INPUTS
!   o KS_eigenvalue
!   o KS_eigenvector
!   o occ_numbers
!   -- others are not used, should be taken out
!    
!  OUTPUT
!   o Prints new eigenvalues.
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
!




      integer:: info


! local variables

      character*120 :: info_str

!  counters

      integer :: i_spin
      integer :: i_region
      integer :: i_state

      
!  begin work

      write(info_str,'(2X,A)') ''
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,A)') &
       "Unshifting Kohn-Sham eigenvalues."
      call localorb_info(info_str,use_unit,'(A)')
               
              do i_spin = 1, n_spin, 1
                 do i_state = 1, n_states, 1
                   KS_eigenvalue(i_state,i_spin,1) =   &
                   KS_eigenvalue(i_state,i_spin,1) - &
                   constraint_potential(1, i_spin)
                 enddo
                 chemical_potential_spin(i_spin) = chemical_potential - constraint_potential(1, i_spin)
              enddo


          if(output_level .ne. 'full') output_priority = 0
          call output_real_eigenfunctions &
              ( KS_eigenvalue, KS_eigenvector, occ_numbers  )
          if(output_level .ne. 'full') output_priority = 1
end subroutine unshift_eigenvalues
!******
