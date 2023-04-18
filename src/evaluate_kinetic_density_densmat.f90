!****s* FHI-aims/evaluate_kinetic_density_densmat
!  NAME
!    evaluate_kinetic_density_densmat
!  SYNOPSIS

subroutine evaluate_kinetic_density_denmat( & 
     n_points, gradient_basis_wave, n_compute, &
     density_matrix,kinetic_density_own, n_basis_list)

!  PURPOSE
!  Evaluates kinetic density for one grid batch using density matrix formalism.
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_basis_list
  integer :: n_points ! number of points in this batch
  integer :: n_compute ! number of non zero basis functions in this batch

  real*8, dimension(n_basis_list, 3, n_points) :: gradient_basis_wave ! input: the gradient of the basis
  real*8, dimension(n_compute, n_compute) :: density_matrix ! input: the part of the density matrix relevant for this batch

  real*8, dimension(n_points) :: kinetic_density_own ! output: the kinetic density for this batch

  real*8, dimension(3, n_points) :: kinetic_density_vec ! temporary variable that stores the kinetic density before summing up in the coordinates
  real*8:: vec_work(n_compute, n_points) ! temporary variable that stores the gradient_basis_wave for each coordinate
  real*8:: work(n_compute, n_points) ! temporary variable that stores the result of DM * grad phi


! INPUTS
! o gradient_basis_wave
! o density_matrix - see descriptions above
!  OUTPUT
! o temp_kinetic_density - see descriptions above
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

  kinetic_density_own(:)=0.d0

  do i_coord = 1, 3, 1

     vec_work =  gradient_basis_wave(1:n_compute, i_coord, 1:n_points)

     ! do DM * grad phi for all basis n

     call dsymm('L','U', n_compute, n_points, 1.0d0,  &
          density_matrix, n_compute, vec_work, n_compute, &
          0.0d0, work, n_compute)

     do i_point = 1, n_points,1
  
        ! now do (DM * grad phi) * grad phi for all basis m, for each point

        kinetic_density_vec(i_coord, i_point) = & 
          dot_product(work(1:n_compute, i_point), gradient_basis_wave(1:n_compute, i_coord, i_point))

       ! now sum over the coordinates, for each point

        kinetic_density_own(i_point) = kinetic_density_own(i_point)+kinetic_density_vec(i_coord, i_point)

     end do

  end do


end subroutine evaluate_kinetic_density_denmat
