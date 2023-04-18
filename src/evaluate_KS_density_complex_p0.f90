!****s* FHI-aims/evaluate_KS_density_complex_p0
!  NAME
!   evaluate_KS_density_complex_p0
!  SYNOPSIS

subroutine evaluate_KS_density_complex_p0 &
     ( n_points, wave, n_compute, i_basis, &
     KS_vec_times_occ_sqrt, KS_ev_compute, max_occ_number, &
     KS_orbital, rho, n_basis_compute, n_max_basis_T &
     )

!  PURPOSE
!  Subroutine takes the local geometry and basis functions at a given integration point, 
!  and evaluates the charge density. This is subroutine for complex format eigevectors.
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: n_basis_compute
  real*8, dimension(n_basis_compute, n_points) :: wave
  complex*16, dimension(n_basis_compute, n_points) :: wave_c

  integer :: n_compute
  integer :: i_basis(n_compute)

  integer :: n_max_basis_T
  integer :: max_occ_number
  complex*16, dimension(n_max_basis_T, max_occ_number) :: KS_vec_times_occ_sqrt

  complex*16, dimension(max_occ_number, n_points) :: KS_orbital
  real*8 :: rho(n_points)
  complex*16, dimension(max_occ_number, n_compute):: KS_ev_compute




! INPUTS
! o n_points -- number of integration points
! o n_basis_compute -- maximum number of relevant basis functions
! o wave -- one part of basis functions
! o wave_c -- other part of basis functions
! o n_compute -- number of relevant basis functions
! o i_basis -- list of relevant basis functions
! o n_max_basis_T -- total number of basis functions
! o max_occ_number -- maximum number of occupated states
! o KS_vec_times_occ_sqrt -- Kohn-Sham eigenvector multiplied with sqrt(occupation)
!  
! OUTPUT
! o KS_orbital -- Kohn-Sham orbitals
! o rho -- electron density
! o KS_ev_compute  -- relevant part of Kohn-Sham eigenvector multiplied with sqrt(occupation)
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
  integer :: i_compute
  integer :: i_point

  !      real*8, dimension(n_compute, n_points) :: wave_compute

  !  begin work

  !     FIXME: In the matrix-vector product below, we fail to take advantage of the benefits of
  !            compute() . This is an obvious place for improvement, but not by a simple if statement.
  !            Must condense wave() array itself already in evaluate_waves, to house only non-zero values, rest can be done via 
  !            index arrays.

  ! costs time ...      KS_orbital = 0.0d0

  ! must condense eigenvectors to only the needed ones ...

  do i_compute = 1, n_compute, 1
     do i_state = 1, max_occ_number, 1

        KS_ev_compute(i_state, i_compute) =  &
             KS_vec_times_occ_sqrt(i_basis(i_compute), i_state)

     enddo
  enddo


  wave_c = dcmplx(wave)


  call zgemm('N','N', max_occ_number, n_points, n_compute, (1.0d0,0d0),  &
       KS_ev_compute, max_occ_number, wave_c, n_basis_compute, &
       (0.0d0,0d0), KS_orbital, max_occ_number)

  do i_point = 1, n_points, 1


     rho(i_point) = dble(dot_product(KS_orbital(1:max_occ_number,i_point), &
          KS_orbital(1:max_occ_number,i_point)))

  enddo


end subroutine evaluate_KS_density_complex_p0
!---------------------------------------------------------------------
!******	
