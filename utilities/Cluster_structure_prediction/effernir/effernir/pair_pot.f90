! pair_potential is just a splined binding curve of a dimer
! to assess features of the pair-potential-community
! like pair-energy, angular-move, etc.. (see D.Wales)
! makes hopefully sense
!
! R.Gehrke (2007)
!

module pair_potential

  use spline
  use cluster

  implicit none

  real*8, dimension(:,:,:,:), allocatable :: pot_spl
  real*8, dimension(:,:,:), allocatable :: pot_par
  real*8, dimension(:,:), allocatable :: d_min
  real*8, dimension(:,:), allocatable :: d_max
  real*8, dimension(:,:), allocatable :: delta_d
  integer, dimension(:,:), allocatable  :: n_data
  integer :: n_max_data

  contains

    subroutine allocate()

      write (*,*) "Allocate memory for lennard jones data"
      if ((n_atoms .gt. 0) .and. (n_species .gt. 0)) then
         if (.not.allocated(pot_spl)) then
            allocate(pot_spl(4, n_max_data, n_species, n_species))
         end if
         if (.not.allocated(pot_par)) then
            allocate(pot_par(n_max_data, n_species, n_species))
         end if
         if (.not.allocated(n_data)) then
            allocate(n_data(n_species, n_species))
         end if
         if (.not.allocated(d_min)) then
            allocate(d_min(n_species, n_species))
         end if
         if (.not.allocated(d_max)) then
            allocate(d_max(n_species, n_species))
         end if
         if (.not.allocated(delta_d)) then
            allocate(delta_d(n_species, n_species))
         end if
      else
         write (*,*) "No atoms and species found yet."
         write(*,*) "* Aborting."
         stop
      end if
      
    end subroutine allocate

    subroutine deallocate()

      if (allocated(pot_spl)) then
         deallocate(pot_spl)
      end if
      if (allocated(pot_par)) then
         deallocate(pot_par)
      end if
      if (allocated(n_data)) then
         deallocate(n_data)
      end if
      if (allocated(d_min)) then
         deallocate(d_min)
      end if
      if (allocated(d_max)) then
         deallocate(d_max)
      end if
      if (allocated(delta_d)) then
         deallocate(delta_d)
      end if

    end subroutine deallocate

    subroutine initialize()

      ! local variables
      integer :: i_species
      integer :: i_species_2
      
      real*8 :: r_index
      real*8 :: d

      ! counter
      integer :: i_counter

      do i_species = 1, n_species, 1
         do i_species_2 = i_species, n_species, 1
            call cubic_spline(pot_par(1, i_species, i_species_2), n_data(i_species, i_species_2), pot_spl(1,1,i_species, i_species_2))
            d_max(i_species, i_species_2) = d_min(i_species, i_species_2) + (n_data(i_species, i_species_2) - 1) * delta_d(i_species, i_species_2) 
            if (i_species_2 .ne. i_species) then
               pot_spl(:,:,i_species_2, i_species) = pot_spl(:,:,i_species, i_species_2)
               d_max(i_species_2, i_species) = d_max(i_species, i_species_2)
            end if
!            do i_counter = 0, 100, 1
!               d = d_min(i_species, i_species_2) + i_counter * 0.01
!               r_index = (d - d_min(i_species, i_species_2)) / delta_d(i_species, i_species_2) + 1
!               write (6,*) "d= ", d, " ", val_spline(r_index, &
!                    pot_spl(1,1,i_species, i_species_2), n_data(i_species, i_species_2))
!            end do
         end do
      end do

    end subroutine initialize

end module pair_potential
