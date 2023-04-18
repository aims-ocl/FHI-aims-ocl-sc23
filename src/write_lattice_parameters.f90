!****f* FHI-aims/write_lattice_parameters
!*  NAME
!*    write_lattice_parameters
subroutine write_lattice_parameters(lattice_vector)
!*  SYNOPSIS
!*    call write_lattice_parameters(lattice_vector)
!*  FUNCTION
!*    Writes the crystallographic lattice parameters to standard output.
!*    Trivial functionality, but I've grown tired of re-calculating it
!*    in every new script I write.
!*  USES
  use dimensions
  use constants
  use localorb_io
  implicit none
!*  ARGUMENTS
!*    lattive_vector(3, n_periodic) - the lattice vectors
!*  INPUTS
!*    lattive_vector(3, n_periodic) - the lattice vectors
!*  OUTPUTS
!*    None
!*  AUTHORS
!*    William Huhn (Duke University)
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject
!*    the terms and conditions of the respective license agreement."
!*  SOURCE

  real*8, dimension(3,n_periodic), intent(in) :: lattice_vector

  character*300 :: info_str
  real*8 :: len_a, len_b, len_c
  real*8 :: ang_alpha, ang_beta, ang_gamma ! Turns out gamma is a Fortran keyword.  Who knew?

  if (n_periodic .lt. 1) then
    ! Non-periodic case, do nothing
    return
    
  else if (n_periodic .eq. 1) then
    len_a = sqrt(dot_product(lattice_vector(:,1),lattice_vector(:,1)))
    write(info_str, "(2X, A, F12.6)") "Lattice parameters for 1D lattice (in Angstroms) : ", &
         len_a*bohr
    call localorb_info( info_str )

  else if (n_periodic .eq. 2) then
    len_a = sqrt(dot_product(lattice_vector(:,1),lattice_vector(:,1)))
    len_b = sqrt(dot_product(lattice_vector(:,2),lattice_vector(:,2)))
    if (len_a .gt. 0 .and. len_b .gt. 0) then
      ang_gamma = acos(dot_product(lattice_vector(:,1),lattice_vector(:,2))/len_a/len_b) * 180. / pi 
    else
      ang_gamma = 0.
    end if
    write(info_str, "(2X, A, F12.6, F12.6)") "Lattice parameters for 2D lattice (in Angstroms) : ", &
         len_a*bohr, len_b*bohr
    call localorb_info( info_str )
    write(info_str, "(2X, A, F12.6)")        "Angle(s) between unit vectors (in degrees)       : ", &
         ang_gamma
    call localorb_info( info_str )

  else if (n_periodic .eq. 3) then
    len_a = sqrt(dot_product(lattice_vector(:,1),lattice_vector(:,1)))
    len_b = sqrt(dot_product(lattice_vector(:,2),lattice_vector(:,2)))
    len_c = sqrt(dot_product(lattice_vector(:,3),lattice_vector(:,3)))

    if (len_b .gt. 0 .and. len_c .gt. 0) then
      ang_alpha = acos(dot_product(lattice_vector(:,2),lattice_vector(:,3))/len_b/len_c) * 180. / pi
    else
      ang_alpha = 0.
    end if
    if (len_a .gt. 0 .and. len_c .gt. 0) then
      ang_beta  = acos(dot_product(lattice_vector(:,1),lattice_vector(:,3))/len_a/len_c) * 180. / pi
    else
      ang_beta  = 0.
    end if
    if (len_a .gt. 0 .and. len_b .gt. 0) then
      ang_gamma = acos(dot_product(lattice_vector(:,1),lattice_vector(:,2))/len_a/len_b) * 180. / pi
    else
      ang_gamma = 0.
    end if
    write(info_str, "(2X, A, F12.6, F12.6, F12.6)") "Lattice parameters for 3D lattice (in Angstroms) : ", &
         len_a*bohr, len_b*bohr, len_c*bohr
    call localorb_info( info_str )
    write(info_str, "(2X, A, F12.6, F12.6, F12.6)") "Angle(s) between unit vectors (in degrees)       : ", &
         ang_alpha, ang_beta, ang_gamma
    call localorb_info( info_str )

  else 
    ! Some kind of crazy 4+ dimenstional DFT is going on.
    ! Which actually sounds pretty cool.
    ! But I don't think there's any grant funding for that.
    ! So skip it.
    write(info_str, "(A)")  " Normally lattice parameters are output here, but you're working with", &
         " dimensionalities that exceed human understanding."
    call localorb_info( info_str )

  end if

end subroutine write_lattice_parameters
