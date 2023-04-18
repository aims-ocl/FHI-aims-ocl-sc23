! module to read energies and forces from
! file "dft_data.dat" as long as provided,
! otherwise ask for new DFT-energy by aims via
! the status-flag
!
! R.Gehrke (2006)

module dft_file

  use cluster
  
  implicit none
  
  integer :: file_call
  integer :: max_file_call
  integer :: seed_uni
  integer :: seed_poiss
  
contains
  
  subroutine initialize_dft_file()

    open (30, FILE = "startindex.dat")
    read (30,*) max_file_call, seed_uni, seed_poiss
    close (30)

    if (seed_poiss .lt. 0) then
       write (6,*) "Negative seed for poisson-RNG senseless!"
       write (6,*) "Change sign."
       seed_poiss = - seed_poiss
    end if
    
    open (30, FILE = "dft_data.dat")
    file_call = 0
    
  end subroutine initialize_dft_file
  
  subroutine cleanup_dft_file()
    
    close (30)
    
  end subroutine cleanup_dft_file
  
  subroutine rewind_dft_file()

    rewind (30)
    file_call = 0

  end subroutine rewind_dft_file

  subroutine get_dft_data(coords, energy, forces, spin, status)
    
    ! imported variables
    
    ! output
    type (vector_array), intent(out) :: coords
    type (vector_array), intent(out) :: forces
    real*8, intent(out) :: energy
    integer, intent(out) :: status
    real*8, dimension(n_atoms), intent(out) :: spin
    
    ! counter
    integer :: i_atom
    
    file_call = file_call + 1

    if (file_call .le. max_file_call) then

       do i_atom = 1, n_atoms, 1
          read (30, *) coords%coords(i_atom)%x, coords%coords(i_atom)%y, coords%coords(i_atom)%z
       end do
       read (30, *) energy, status
       do i_atom = 1, n_atoms, 1
          read (30, *) forces%coords(i_atom)%x, forces%coords(i_atom)%y, forces%coords(i_atom)%z
       end do
       if (spin_polarized) then
          do i_atom = 1, n_atoms, 1
             read (30, *) spin(i_atom)
          end do
       end if
    else
       ! no forces in file anymore => call aims
       status = 1
    end if

  end subroutine get_dft_data
  
end module dft_file
