!****s* FHI-aims/check_norm_constraint
!  NAME
!   check_norm_constraint
!  SYNOPSIS

subroutine check_norm_constraint(chemical_potential, KS_eigenvalue, constraint_electrons, &
                                 constraint_proj, occ_numbers, diff_electrons, i_counter)

!  PURPOSE
!   checks and returns the charge density norm (i.e. the electron count) for
!   a given Fermi level, and for whichever smearing type was chosen ...
!  USES

  use dimensions
  use runtime_choices
  use arch_specific
  use constants
  use localorb_io, only: use_unit
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: chemical_potential
  real*8, dimension(n_states), intent(in) :: KS_eigenvalue
  real*8, dimension(n_states), intent(in) :: constraint_proj  
  real*8, intent(in) :: constraint_electrons
  real*8, dimension(n_states), intent(out) :: occ_numbers
  real*8, intent(out) :: diff_electrons
  integer, intent(inout) :: i_counter 

!  INPUTS
!   o chemical_potential -- chemical potential
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o constraint_proj -- ?????????
!   o constraint_electrons -- number of nostrain electrons
!
!  OUTPUT
!   o occ_numbers -- occupation numbers
!   o diff_electrons -- ???????????
!   o i_counter -- ?????????????
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


  
  !  local variables
  real*8 :: temp_n_electrons
  real*8 :: one_over_ow
  real*8 :: H_odd
  real*8 :: H_even
  real*8 :: hermite_arg
  real*8 :: gauss_weight
  real*8 :: A

     
  !  counter
  integer :: i_state
  integer :: i_spin

  integer :: i_mp
  
  i_counter = i_counter + 1
  one_over_ow   = 1. / occupation_width
  temp_n_electrons   = 0.d0
  
  select case(occupation_type)
  case (0)
     !  gaussian smearing (0 <= occ_numbers )
!     do i_spin =1, n_spin, 1
        do i_state = 1, n_states, 1
            occ_numbers(i_state) = spin_degeneracy * 0.5d0 * &
           (1 - arch_erf((KS_eigenvalue(i_state) - chemical_potential) * one_over_ow))
           temp_n_electrons = temp_n_electrons + constraint_proj(i_state)*occ_numbers(i_state)
        end do
!     end do
     
  case (1)
     !  fermi smearing (0 <= occ_numbers )
!     do i_spin =1, n_spin, 1
        do i_state = 1, n_states, 1
           occ_numbers(i_state) = spin_degeneracy / &
           (1 + exp((KS_eigenvalue(i_state) - chemical_potential) * one_over_ow))
           temp_n_electrons = temp_n_electrons + occ_numbers(i_state) * constraint_proj(i_state)
        end do
!     end do
  case (2)
     ! Methfessel-Paxton smearing (0 <= occ_numbers <= pin_degeneracy)
!     do i_spin =1, n_spin, 1
        do i_state = 1, n_states, 1
           hermite_arg  = (KS_eigenvalue(i_state) - chemical_potential) * one_over_ow 
           gauss_weight = exp(- hermite_arg * hermite_arg)
           !     zero order contribution
           occ_numbers(i_state) = 0.5 * (1 - arch_erf(hermite_arg))
           if (n_methfessel_paxton .gt. 0) then
!              write(use_unit,*) A, H_even, H_odd, hermite_arg
              !     first order contribution
              A = - 0.25 * pisqrt_inv
              !     H_even = H_0 = 1
              H_even = 1.d0
              !     H_odd =  H_1 = 2 * x
              H_odd  = 2 * hermite_arg
              occ_numbers(i_state) = occ_numbers(i_state) + A * H_odd * gauss_weight
           end if
           if (n_methfessel_paxton .gt. 1) then
              do i_mp = 2, n_methfessel_paxton, 1
                 A = - 1. / dble(4 * i_mp) * A
                 H_even = 2 * hermite_arg * H_odd  - 2 *  i_mp      * H_even
                 H_odd  = 2 * hermite_arg * H_even - 2 * (i_mp + 1) * H_odd
                 occ_numbers(i_state) = occ_numbers(i_state) + A * H_odd * gauss_weight
              end do
           end if
!     factor 2 due to spin
!            occ_numbers(i_state,i_spin) = occ_numbers(i_state,i_spin) * spin_degeneracy
!            write(use_unit,*) "occ = ", occ_numbers(i_state,i_spin)
            temp_n_electrons = temp_n_electrons + occ_numbers(i_state) * constraint_proj(i_state)
         end do
!     end do
 
  case default
         write(use_unit,*) "* Unknown smearing type in subroutine check_norm_constraint."
         write(use_unit,*) "* Abort."
         stop

  end select

!  do i_spin = 1, n_spin, 1
    diff_electrons = temp_n_electrons - constraint_electrons
!  end do
  

end subroutine check_norm_constraint

!-------------------------------------------------------------------------------
!******
