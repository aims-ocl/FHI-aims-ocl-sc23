!****s* FHI-aims/evaluate_KS_orbital_density_complex
!  NAME
!    evaluate_KS_orbital_density_complex
!  SYNOPSIS

subroutine evaluate_KS_orbital_density_complex_p0 &
     (orbital_number, n_points, wave, n_compute, i_basis, &
     KS_vector, max_num_calc_orb,KS_orbital,rho, &
     n_basis_compute, n_max_basis_T)

!  PURPOSE
!  Subroutine evaluate_KS_orbital_density takes the local geometry and basis
!  functions at a given integration point and evaluates the charge density
!  for a single KS orbital
!
!  USES

  use dimensions
  use mpi_tasks, only: STDERR
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: n_basis_compute
  real*8, dimension(n_basis_compute, n_points) :: wave
  complex*16, dimension(n_basis_compute, n_points) :: wave_c

  integer :: n_compute
  integer :: i_basis(n_compute)
  integer :: orbital_number
  integer :: max_num_calc_orb
  integer :: n_max_basis_T

  complex*16, dimension(n_max_basis_T, n_states) :: KS_vector

  complex*16, dimension(max_num_calc_orb, n_points) :: KS_orbital
  complex*16, dimension(:,:),allocatable :: KS_ev_compute
  real*8 :: rho(n_points)


! INPUTS
!  o n_points -- number of grid points
!  o wave -- basis functions
!  o n_compute -- number of relavant basis functions
!  o i_basis -- list of relevant basis functions
!  o orbital_number -- ??????
!  o max_num_calc_orb -- ?????????
!  o KS_vector -- Kohn-Sham eigenvectors
!
!  OUTPUT
!  o KS_orbital -- Kohn-Sham orbitals
!  o KS_ev_compute -- ??????
!  o rho -- electron density
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


  !     local variables
  real*8, external :: ddot

  !     counters
  integer :: i_state, info
  integer :: i_compute
  integer :: i_point
  real*8, dimension(n_compute, n_points) :: wave_compute


  if(.not. allocated( KS_ev_compute))then
     allocate( KS_ev_compute(max_num_calc_orb,n_compute),stat=info)
     if(info/=0)then
        write(STDERR,*)'Error in allocation: KS_ev_compute'
        stop
     end if
  end if


  do i_compute = 1, n_compute, 1
     KS_ev_compute(1,i_compute) = KS_vector( &
          i_basis(i_compute), orbital_number)
  enddo

  wave_c = dcmplx(wave)

  call zgemm('N','N', max_num_calc_orb, n_points, n_compute, (1.0d0,0d0),  &
       KS_ev_compute, max_num_calc_orb, wave_c, n_basis_compute, &
       (0.0d0,0d0), KS_orbital, max_num_calc_orb)

  do i_point = 1, n_points, 1
     rho(i_point) = Real(KS_orbital(1,i_point))**2.d0+Aimag(KS_orbital(1,i_point))**2.d0
  enddo

  if(allocated( KS_ev_compute)) deallocate( KS_ev_compute )

  !     end work
end subroutine evaluate_KS_orbital_density_complex_p0
!******	
