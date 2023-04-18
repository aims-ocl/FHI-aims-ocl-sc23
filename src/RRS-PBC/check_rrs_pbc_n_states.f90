!****s* FHI-aims/check_rrs_pbc_n_states
!  NAME
!    check_rrs_pbc_n_states
!  SYNOPSIS

subroutine check_rrs_pbc_n_states(n_states, n_spin, n_k_points, occ_numbers, t_out)

  !  PURPOSE
  !
  !    Check if n_states was large enough.  If there is (significant)
  !    occupation within the last state, it wasn't.
  !
  !  USES

  use localorb_io
  use dimensions, only: rrs_pbc_n_center_basis
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_states, n_spin, n_k_points
  real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
  logical, intent(IN) :: t_out

  !  INPUTS
  !    o n_states, n_spin, n_k_points -- Dimensions off occ_numbers
  !    o occ_numbers -- Occupation numbers
  !    o t_out -- Force output of highes occupation number
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  logical :: is_minimal_basis
  real*8 :: max_inv_occ
  character*150 :: info_str
  character(*), parameter :: func = 'check_rrs_pbc_n_states'

  max_inv_occ = maxval(occ_numbers(n_states, :,:))
  if (max_inv_occ > 1d-6) then
     is_minimal_basis = (n_states >= rrs_pbc_n_center_basis)
     if (is_minimal_basis) then
        write(info_str, "(2X,A,ES10.2,A)") &
        & 'Max occ number of highest state is', max_inv_occ, &
        & ' because we are using a minimal basis.'
        call localorb_info(info_str)
     else
        write(info_str, "(2X,A,ES10.2,A)") &
             & '*** Max occ number of highest state of ', max_inv_occ, &
             & ' for some k-piont is too high!'
        call localorb_info(info_str)
        write(info_str, "(2X,A)") &
        & '*** This means that the number of states kept after ', &
        & 'diagonalization is not large enough.'
        call localorb_info(info_str)
        write(info_str, "(2X,A)") &
        & '*** Make sure that empty_states (per atom!) matches your smearing.'
        call localorb_info(info_str)
        if (max_inv_occ > 1d-3) then
           call aims_stop('Number of states too low.  Stopping.', func)
        end if
     end if
  else if (max_inv_occ > 1d-12 .or. t_out) then
     write(info_str, "(2X,'| ',A,ES10.2)") &
     & 'Maximum occupation number of highest empty state:', max_inv_occ
     call localorb_info(info_str, use_unit, '(A)', OL_norm)
  end if

end subroutine check_rrs_pbc_n_states
!******
