!****s* FHI-aims/get_occupation_numbers_fsm
!  NAME
!    get_occupation_numbers_fsm
!  SYNOPSIS

subroutine get_occupation_numbers_fsm &
  ( KS_eigenvalue, occ_numbers, chemical_potential_spin )

!  PURPOSE
!
!  Based on a given Kohn-Sham eigenvalue spectrum, returns the occupation
!  numbers for the case of a fixed spin moment prescribed in control.in,
!  so that separate chemical potentials (Fermi levels) for each spin 
!  channels must be found.
!
!  Unlike the "other" get_occupation_numbers_* subroutines,
!  get_occupation_numbers_fsm is only a wrapper around (presently) 
!  get_occupation_numbers_single_channel , a subroutine that determines 
!  the occupation numbers only of a single spin channel (for fixed number
!  of electrons in that spin channel), and without any knowledge of 
!  of possible other spin channels at all.
!
!  Why a wrapper subroutine? The array dimension n_states can change along the
!  way in rare cases. With two spin channels, this will implicitly reshape 
!  the array KS_eigenvalues throughout the code - but if KS_eigenvalues were used 
!  directly from physics.f90, the old, formally larger (but factually unused) 
!  dimension would be used to shape the array instead. We thus need this wrapper to
!  declare the correct shape.
!
!  Why not just one get_occupation_numbers subroutine to bind them all? 
!  The trouble is that the chemical potentials are found iteratively, by guessing a value
!  and then counting to see whether the guessed chemical potential produces the right
!  number of electrons. So a sum of occupation numbers must be taken repeatedly, as the
!  chemical potential is adjusted. If the chemical potential is computed per spin channel,
!  the sum must be taken and iterated inside each spin channel. If, instead, the 
!  chemical potential is unique for all spin channels, the sum must be taken over 
!  all spin channels, and iterated outside the loop. So the loop structure of the associated
!  subroutines changes fundamentally.
!
!  USES

  use dimensions
  use runtime_choices
  use localorb_io
  use constants
  implicit none

!  ARGUMENTS

  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue

  real*8, dimension(n_states, n_spin, n_k_points), intent(out) :: occ_numbers
  real*8, dimension(n_spin), intent(out) :: chemical_potential_spin

!  INPUTS
!  o KS_eigenvalue -- Kohn-Sham eigenvalues
!
!  OUTPUT
!  o occ_numbers -- occupation weights of different KS states
!  o chemical_potential_spin -- separate chemical potentials for each spin channel
!  
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175-2196 (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2011).
!  SOURCE

!  local variables
   logical :: t_out

!  counter
  integer :: i_spin

  t_out = .true. ! print output in get_occupation_numbers_single_channel

  ! All we do is to determine the occupation numbers independently for each spin channels,
  ! but with the correct shape of the arrays KS_eigenvalues and occ_numbers declared above.
  do i_spin = 1, n_spin, 1
     call get_occupation_numbers_single_channel( & 
        KS_eigenvalue(1:n_states,i_spin,1:n_k_points), &
        fixed_spin_moment_electrons(i_spin),t_out,occ_numbers(1:n_states,i_spin,1:n_k_points), &
        chemical_potential_spin(i_spin),i_spin & 
        )
  end do

end subroutine get_occupation_numbers_fsm

!-------------------------------------------------------------------------------
!******	
