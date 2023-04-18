!----------------------------------------------------------------------
! Subroutine get_constraint_fermi
!
! determines occupation numbers according to two available schemes:
! 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
! 2) fermi smearing
! 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
! R. Gehrke (2005)
!---------------------------------------------------------------------- 

subroutine get_constraint_fermi(KS_eigenvalue,constraint_electrons,constraint_proj,constraint_fermi)
  
  use dimensions
  use runtime_choices
  use localorb_io
  use constants


  implicit none
  
  
  !  imported variables

  !  input 
  
  real*8, dimension(n_states), intent(in) :: KS_eigenvalue
  real*8, dimension(n_states), intent(in) :: constraint_proj
  real*8, intent(in) :: constraint_electrons

  !  output
  
  real*8, intent(out) :: constraint_fermi

!  local variables
  
  real*8, dimension(n_states) :: occ_numbers
  real*8 :: chemical_potential

  real*8 :: degeneracy_threshold
  integer :: n_degenerate
  integer :: highest_full_state
  integer :: lowest_empty_state
  real*8 :: shared_electrons
  logical :: degenerate
  real*8  :: chemical_potential_l
  real*8  :: chemical_potential_r
  real*8 :: diff_l
  real*8 :: diff_r
  real*8 :: diff_electrons
  real*8 :: diff_electrons_thr
  real*8 :: lowest_eigenvalue
  real*8 :: highest_eigenvalue
  
  real*8 :: electron_count
  
  character*120 :: info_str

  !  counters
  
  integer :: i_state
  integer :: i_state_2
!  integer :: i_spin
  integer :: i_counter
  
!  write(use_unit,'(2X,A)') "Determining occupation numbers for Kohn-Sham eigenstates."

  i_counter = 0

!  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states
!  write(use_unit,*)"CONSTRAINT ELECTRONS: ",constraint_electrons

  lowest_eigenvalue = KS_eigenvalue(1)
  do i_state = 1, n_states, 1
     if (KS_eigenvalue(i_state) .lt. lowest_eigenvalue) then
        lowest_eigenvalue = KS_eigenvalue(i_state)
     end if
  end do

  highest_eigenvalue = KS_eigenvalue(n_states)
  do i_state = 1, n_states, 1
     if (KS_eigenvalue(i_state) .gt. highest_eigenvalue) then
        highest_eigenvalue = KS_eigenvalue(i_state)
     end if
  end do
 



!  initialize zeroin algorithm
!  to avoid problems if only one state is calculated (so for the H-atom)
!  do not initialize mit KS_eigenvalue(n_states), because it is then
!  identical to KS_eigenvalue(1)

  chemical_potential_r = highest_eigenvalue
  chemical_potential_l = lowest_eigenvalue
  call check_norm_constraint(chemical_potential_l, KS_eigenvalue, constraint_electrons, &
                             constraint_proj, occ_numbers, diff_l, i_counter)
  call check_norm_constraint(chemical_potential_r, KS_eigenvalue, constraint_electrons, & 
                             constraint_proj, occ_numbers, diff_r, i_counter)
  do while (diff_l * diff_r .gt. 0.d0)

!     interval for chemical potential still not found

     chemical_potential_l = chemical_potential_r
     chemical_potential_r = chemical_potential_r + 0.5 * dabs(highest_eigenvalue-lowest_eigenvalue)
     diff_l = diff_r
     call check_norm_constraint(chemical_potential_r, KS_eigenvalue, constraint_electrons, &
                                constraint_proj, occ_numbers, diff_r, i_counter)

!     write(use_unit,*) chemical_potential_l, chemical_potential_r
     !     If the right interval cannot be found, stop potentially
     !     infinite loop
     
     if (i_counter .gt. max_zeroin) then
        write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
              call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* Too many iterations needed to find good interval around E_F."
              call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* According to this check, you may not be able to find one."
              call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "Please look carefully at your settings."
              call localorb_info(info_str,use_unit,'(A)')
        stop
     end if
  end do
      
  !test
  !     write(use_unit,*) "Before zeroin."
  !     write(use_unit,*) "left: ", chemical_potential_l
  !         write(use_unit,*) "right: ", chemical_potential_r
  !         write(use_unit,*) "occ_numbers : ", occ_numbers
  !test end 
  
  call zeroin_constraint(chemical_potential_l, chemical_potential_r, diff_l, &
       diff_r, fermi_acc, occupation_acc, max_zeroin, &
       n_states, KS_eigenvalue, constraint_electrons, chemical_potential, &
       diff_electrons, constraint_proj, occ_numbers, i_counter)
      
  !     call check_norm one more time, just to make sure we really have the 
  !     occupation numbers for the final fermi level [this is critical for H]

  call check_norm_constraint( chemical_potential, KS_eigenvalue, constraint_electrons, &
                              constraint_proj, occ_numbers, diff_electrons, i_counter)

  !     and decrement i_counter which was falsely incremented ...

  i_counter = i_counter - 1
         
  !       If zeroin did not converge to the desired accuracy, issue
  !       warning and stop, for now; we can create nicer behavior after
  !       we gain more experience.

  if (i_counter .gt. max_zeroin) then
     write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
              call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') "* Too many iterations needed to find correct occupation of levels."
              call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') &
       "* Check variables occupation_acc, max_zeroin before continuing your calculation."
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* Status of get_occ_numbers when aborted: "
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
              call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* State #         Eigenvalue       Occ. number"
              call localorb_info(info_str,use_unit,'(A)')

!        write(use_unit,*) "* spin # ", i_spin
        write (info_str,*) "----------"
        do i_state = 1, n_states, 1
           write (info_str,*) i_state, KS_eigenvalue(i_state), occ_numbers(i_state)
        end do
     
     stop
  end if
  
  !     finally, occupation thresholding - if needed

  if (occupation_thr.gt.0.d0) then
     call threshold_occ_numbers_constraint( constraint_electrons, occ_numbers, diff_electrons_thr)
  end if
      
  constraint_fermi = chemical_potential

!  write(use_unit,*)"PROJ_OP, OCC_NUMBERS, OCC_NUMBERS*PROJ_OP, KS_EIGENVALUES"
!  do i_state = 1, n_states, 1
!     write(use_unit,*)"STATE:", i_state
!     write(use_unit,*) constraint_proj(i_state),occ_numbers(i_state),occ_numbers(i_state)*constraint_proj(i_state),KS_eigenvalue(i_state)
!  end do
!  write(use_unit,*)"CONSTRAINT FERMI:", constraint_fermi
!  write(use_unit, '(2X, A, E10.4)') "| Chemical potential (Fermi level) in eV                 : ", chemical_potential * hartree
!  write(use_unit, '(2X, A, E10.4)') "| Error in electron count due to remaining E_F inaccuracy: ", diff_electrons
!  if (occupation_thr.gt.0.d0) then
!     write(use_unit, '(2X, A, E10.4)') "| Error in electron count after thresholding : ", diff_electrons_thr
!  end if
!  write(use_unit,*)
  
end subroutine get_constraint_fermi
      
!-------------------------------------------------------------------------------
      
