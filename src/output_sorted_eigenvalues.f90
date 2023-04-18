!****s* FHI-aims/output_KS_dos
!  NAME
!   output_KS_dos
!  SYNOPSIS

subroutine output_sorted_eigenvalues &
     ( KS_eigenvalue, occ_numbers, chemical_potential)

! PURPOSE
! The subroutine prints to file the Kohn-Sham density of states, derived from an advanced version of 
! output_mulliken.f90.
!
! Note that the present version, same as the partial density of states, does NOT include
! support scalapack. This means scalapack type of eigenvectors. If normal form of eigenvectors
! are created from scalapack routine then this this can be called.
!
!  USES

  use dimensions
  use constants
  use mpi_tasks
  use localorb_io
  use basis
  use geometry
  use species_data
  use runtime_choices
  use pbc_lists
  use synchronize_mpi
  implicit none

!  ARGUMENTS

  real*8, dimension(n_states, n_spin, n_k_points)               :: KS_eigenvalue 
  real*8, dimension(n_states, n_spin,n_k_points)                :: occ_numbers
  real*8:: chemical_potential
  real*8  :: n_electrons

!  INPUTS
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o occ_numbers -- occupations weights of the Kohn-Sham eigenstates
!   o n_electrons -- total number of electrons in the system
!  OUTPUT
!   o chemical_potential -- chemical potential
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

  ! local variables

  ! counters
  integer :: i_spin, i_state, i_counter, i_value
  integer :: n_values
  integer :: i_k_point
  character*100 :: info_str
  integer :: n, newn

  real*8, allocatable :: Sorted_Eigenvalues(:,:)
  real*8 :: swap_me(1,2)
  real*8 :: max_weight

  !!!VB: I removed the call to check_norm , for now.
  !      * For one thing, I do not see how occ_numbers are used in this routine.
  !      * Then, the Fermi level may be spin dependent if a spin constraint is used - the 
  !      results below are then wrong.
  !      * The occupation numbers must be updated whenever the KS eigenvalues change,
  !        right away, not here
  !      * Finally, if the KS eigenvalues changed, the Fermi level should itself be updated. The
  !        routine below only recomputes occupation numbers for a fixed Fermi level, the Fermi level
  !        does not change.
  !
  ! calculate position of Fermi level, among other things ... 
  ! i_counter = 0
  ! call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)

  n_values=n_states*n_k_points*n_spin
  if(.not. allocated (Sorted_Eigenvalues)) allocate(Sorted_Eigenvalues(n_values,2)) !first is the eigenvalue, second the maximum occupation 
  
  !First, populate the list unsorted
  i_value=0
!  write(use_unit,*) ' Debug: Populate'
  do i_state=1,n_states,1
   do i_k_point=1,n_k_points,1
     do i_spin=1,n_spin,1
       !calculate maximum occupation
       max_weight=k_weights(i_k_point)
       if (n_spin.eq.1) max_weight=max_weight*2 !If spin-unpolarized, state can hold 2 electrons instead of one
       
       i_value=i_value+1
       Sorted_Eigenvalues(i_value,1)=KS_eigenvalue(i_state,i_spin,i_k_point)
       Sorted_Eigenvalues(i_value,2)=max_weight

     enddo
   enddo
  enddo

   !Next,bubblesort. There are better ways, but this is efficient enough for now.
!   write(use_unit,*) 'Debug: Sort'
   n=n_values
   do while (n.gt.1)
      newn=1
      do i_value=1,n-1,1
         if (Sorted_Eigenvalues(i_value,1).gt.Sorted_Eigenvalues(i_value+1,1)) then
          !SWAP
            swap_me(1,:)=Sorted_Eigenvalues(i_value,:)
            Sorted_Eigenvalues(i_value,:)=Sorted_Eigenvalues(i_value+1,:)
            Sorted_Eigenvalues(i_value+1,:)=swap_me(1,:)
            newn=i_value+1
         endif
         n=newn
      enddo
    
   enddo

!   write(use_unit,*) 'Debug: Sort done, starting to write file'
   
  

  if(myid.eq.0)then
    open(unit=80,file='Sorted_Eigenvalues.dat' )
    write(80,fmt='(2X,A,2X,A)') 'Eigenvalue [eV]', 'Maximum Occupation [e]'

    do i_value=1,n_values,1
       write(80,fmt='(2X,F10.3,2X,F20.9)') &  !@OS: The final number has to be adjusted for your needs. F20.9 means that the output has 20 digits in total, with 9 printed after the comma. 
               (Sorted_Eigenvalues(i_value,1)-chemical_potential)*hartree, &
                Sorted_Eigenvalues(i_value,2)
    enddo

    close(80)
  end if

  if (allocated(Sorted_Eigenvalues)) deallocate(Sorted_Eigenvalues)

end subroutine output_sorted_eigenvalues
!******	
