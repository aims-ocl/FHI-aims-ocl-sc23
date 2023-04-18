!****f* FHI-aims/get_distance
!  NAME
!    get_distance
!  SYNOPSIS

real*8 function get_distance &
     ( coord_current, coord_other &
     )

!  PURPOSE
!  Function get_distance computes the distance between two 3D cartesian points
!
  implicit none
!  ARGUMENTS

    real*8 coord_current(3)
    real*8 coord_other(3)

!  INPUTS
!    o coord_current -- first point coordinates
!    o coord_other -- second point coordinates
!  OUTPUT
!    o get_distance -- distance between points
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


    !  counters

      integer i_coord

      !  begin work

!     get distance from i_atom_2 to current integration point
      get_distance = 0.d0
      do i_coord = 1,3,1

        get_distance = get_distance + &
          (coord_current(i_coord) - &
          coord_other(i_coord)   )**2.0d0

      enddo
      get_distance = sqrt(get_distance)

!  that's all folks

      return
    end function get_distance
!******
!----------------------------------------------------------------------
