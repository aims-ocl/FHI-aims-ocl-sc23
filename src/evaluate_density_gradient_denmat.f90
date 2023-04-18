!****s* FHI-aims/evaluate_density_gradient_denmat
!  NAME
!   evaluate_density_gradient_denmat
!  SYNOPSIS

subroutine evaluate_density_gradient_denmat( & 
     n_points, gradient_basis_wave, n_compute, &
     density_matrix, rho_gradient, wave,  n_basis_list,work, vec_work )

!  PURPOSE
!  Evaluates density gradient for one grid batch using density matrix formalism.
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_basis_list
  integer :: n_points
  integer :: n_compute

  real*8, dimension(n_basis_list, 3, n_points) :: gradient_basis_wave

  real*8:: density_matrix(n_compute, n_compute)

  real*8, dimension(3, n_points) :: rho_gradient
  real*8, dimension(n_basis_list, n_points) :: wave
  real*8:: vec_work(n_compute, n_points)
  real*8:: work(n_compute, n_points)
!MEC_CB new workspace
  real*8:: work2(n_compute, 3, n_points)


! INPUTS
! o n_points -- number of grid points
! o gradient_basis_wave -- gradients of basis functions
! o n_compute -- number of relevant basis functions
! o density_matrix -- relevant numbers of density matrix
! o wave -- basis functions
! o n_basis_list -- 
! o work -- work memory
! o vec_work -- work memory
!
!  OUTPUT
! o rho_gradient -- gradient of electron density
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


  integer :: i_coord
  integer :: i_point

  !  begin work


!NEC_CB  do i_coord = 1, 3, 1
!NEC_CB
!NEC_CB     vec_work =  gradient_basis_wave(1:n_compute, i_coord, 1:n_points)

!NEC_CB Do all three dsymms at once
     call dsymm('L','U', n_compute, n_points*3, 2.0d0,  &
          density_matrix, n_compute, gradient_basis_wave, n_basis_list, &
          0.0d0, work2, n_compute)
!NEC_CB     call dsymm('L','U', n_compute, n_points, 2.0d0,  &
!NEC_CB          density_matrix, n_compute, vec_work, n_compute, &
!NEC_CB          0.0d0, work, n_compute)

  do i_coord = 1, 3, 1
     do i_point = 1, n_points,1

        rho_gradient(i_coord, i_point) = &
           dot_product(work2(1:n_compute, i_coord,i_point), &
                       wave(1:n_compute, i_point))

     end do
  end do

end subroutine evaluate_density_gradient_denmat
!---------------------------------------------------------------------
!******	 
