!****s* FHI-aims/evaluate_KS_density_densmat
!  NAME
!   evaluate_KS_density_densmat
!  SYNOPSIS

subroutine evaluate_KS_density_densmat &
     ( n_points, wave, n_compute, &
     rho, n_basis_compute, n_max_basis_T, &
     density_matrix, work & 
     )

!  PURPOSE
!  Subroutine evaluates Khan Sham electron density using density matrix formalism.
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: n_basis_compute
  real*8, dimension(n_basis_compute, n_points) :: wave

  integer :: n_compute

  integer :: n_max_basis_T
  integer :: max_occ_number
  real*8, dimension(n_compute,n_compute) :: density_matrix

  real*8 :: rho(n_points)


! INPUTS
! o n_points -- number of grid points
! o n_basis_compute -- maximum number of relevant basis functions
! o wave -- basis functions
! o n_compute -- number of relevant basis functions
! o n_max_basis_T -- total number of basis functions
! o max_occ_number -- maximum number of states with non-zero occupation
! o density_matrix -- density matrix components for relevant basis function
!
! OUTPUT
! o rho -- electron density
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

  real*8, external :: ddot

  !     counters

  integer :: i_state
  integer :: i_compute1, i_compute2
  integer :: i_point

  real*8:: work(n_compute, n_points)

  !  begin work

!Safeguard if nothing to do:
  if (n_compute.eq.0) then
     rho(1:n_points) = 0.0
     return
  else
     call dsymm('L','U', n_compute, n_points, 1.0d0,  &
          density_matrix, n_compute, wave, n_basis_compute, &
         0.0d0, work, n_compute)

     do i_point = 1, n_points,1

        rho(i_point) = &
           dot_product(work(1:n_compute,i_point), wave(1:n_compute,i_point))

    end do
  endif

  !  end work

end subroutine evaluate_KS_density_densmat
!---------------------------------------------------------------------
!******	
