!****s* FHI-aims/evaluate_KS_orbital_gradients_complex_p1
!  NAME
!   evaluate_KS_orbital_gradients_complex_p1
!  SYNOPSIS

subroutine evaluate_KS_orbital_gradients_complex_p1(n_points, &
     gradient_basis_wave, n_compute, &
     KS_ev_compute, max_occ_number, &
     KS_orbital_gradient, n_basis_list )

!  PURPOSE
!  Evaluates gradient of KS_orbitals in complex format
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_basis_list
  integer :: n_points
  integer :: max_occ_number
  integer :: n_compute

  complex*16, dimension(max_occ_number, n_compute) :: KS_ev_compute
  real*8, dimension(n_basis_list, 3, n_points) :: gradient_basis_wave

  complex*16, dimension(n_states*n_k_points_group,n_max_batch_size,3) :: KS_orbital_gradient

!  INPUTS
!  o n_basis_list -- total number of basis functions
!  o n_points -- number of grid poitns
!  o max_occ_number -- maximum number of occupated eigenstates
!  o n_compute -- number of relevant basis functions
!  o KS_ev_compute -- ????????
!  o gradient_basis_wave -- ?????????
!
!  OUTPUT
!  o KS_orbital_gradient --  gradient of Khan Sham orbitals
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


  !  local variables

  complex*16, dimension(n_compute,n_points) :: aux_wave_gradient

  !     counters

  integer :: i_compute
  integer :: i_coord
  integer :: i_point

  !  begin work

  !     We do a straightforward implementation of 
  !     grad(rho) = 2 * Sum_i [f_i * psi_i * grad(psi_i) ]
  !     where psi_i are the (occupied) Kohn-Sham orbitals.

  do i_coord = 1, 3, 1

     !        condense and reorganize wave_gradient (this is expensive!)
     do i_point = 1, n_points, 1
        do i_compute = 1, n_compute, 1

           aux_wave_gradient(i_compute,i_point) = &
                gradient_basis_wave(i_compute,i_coord,i_point)

        enddo
     enddo

    call zgemm('N','N', max_occ_number, n_points, &
       n_compute, (1.0d0,0.d0), KS_ev_compute, max_occ_number,&
       aux_wave_gradient, n_compute, (0.0d0,0.d0), &
       KS_orbital_gradient(1,1,i_coord), n_states*n_k_points_group)

  end do

  !  end work

end subroutine evaluate_KS_orbital_gradients_complex_p1
!---------------------------------------------------------------------
!******	
