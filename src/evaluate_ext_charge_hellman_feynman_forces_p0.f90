! Hellman-Feynman forces on external multipoles
! AT

subroutine evaluate_ext_charge_hellman_feynman_forces_p0(rho, partition_tab, ext_charge_forces)

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
 
  implicit none
  ! imported variables
  
  ! input
  real*8, dimension(n_spin, n_full_points) :: rho
  real*8, dimension(n_full_points) :: partition_tab

! local variable

  real*8 :: point_term
  real*8 :: atomic_term

  integer :: i_multipole
  integer :: i_coords
  integer :: i_point
  integer :: i_spin

  integer :: i_index
  integer :: i_my_batch, i_basis

  integer :: i_full_points

  real*8, dimension(3) :: direction
  real*8, dimension(3) :: coord_current
  real*8 :: distance_squared


  ! output
  real*8, dimension(3, n_multipoles) :: ext_charge_forces 


  i_full_points = 0



! from update_density_and_forces_orbital.f90
  do i_my_batch = 1, n_my_batches, 1

        i_basis = 0

        i_point = 0

        ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1

              !     get current integration point coordinate
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
              do i_spin = 1, n_spin, 1
                point_term = partition_tab(i_full_points) * rho(i_spin, i_full_points)

                do i_multipole = 1, n_multipoles, 1

                  distance_squared = 0.d0
                  do i_coords = 1, 3, 1
                     direction(i_coords) =   &
                          coord_current(i_coords) - multipole_coords(i_coords,i_multipole)
                     distance_squared = distance_squared +   &
                          direction(i_coords) * direction(i_coords)
                  end do

                  if (multipole_order(i_multipole) == 0) then
                    atomic_term = multipole_charge(i_multipole) * point_term &
                                  / (distance_squared ** 1.5d0)
                    do i_coords = 1, 3, 1
                         ext_charge_forces(i_coords, i_multipole) = ext_charge_forces(i_coords, i_multipole) + &
                                    atomic_term * direction(i_coords)
                    end do

                  end if
                end do ! multipole 
              end do   ! spin
           endif       ! part tab
        enddo          ! i_index
     !endif
  enddo

!DB: 08/14/13 commenting this before forces have not jet been synchronized at this point.
!

!  if (myid.eq.0) then
!    write (use_unit,'(2X,A)') "Total external charge forces (derivative of free energy) [a.u.]:"
!    do i_multipole = 1, n_multipoles, 1
!      if (multipole_order(i_multipole) == 0) then
!        write (use_unit,'(2X,A,I4,1X,3(E30.15,1X))') "|",i_multipole, ext_charge_forces(:,i_multipole) 
!      endif
!    enddo
!  endif

end subroutine evaluate_ext_charge_hellman_feynman_forces_p0
