!****s* FHI-aims/tab_trigonom_v2
!  NAME
!   tab_trigonom_v2
!  SYNOPSIS

      subroutine tab_trigonom_v2 &
      ( n_compute_atoms, dir_tab, trigonom_tab &
      )

!  PURPOSE
!  Subroutine tab_trigonom converts the directions between current integration point
!  and all atoms into trigonometric function of the spherical angles.
!
!  We obtain: * Trigonometric functions of spherical angles between current point and each atom,
!               cos(theta), sind(theta), cos(phi), sin(phi)
!
!  USES


      use dimensions
      use geometry
      use grids
      implicit none

!  ARGUMENTS

      integer :: n_compute_atoms
      real*8 dir_tab ( 3, n_atoms )
      real*8 trigonom_tab ( 4, n_atoms )

!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o dir_tab -- direction to atoms
!  OUTPUT
!    o trigonom_tab -- cos(theta), sind(theta), cos(phi), sin(phi)
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
!



!  local variables

      real*8 abmax, abcmax, ab, abc

!     counters

      integer :: i_atom
      integer :: i_coord

!  begin work

      do i_atom = 1, n_compute_atoms, 1

!       Convert unit vector (direction) into trigonometric functions
!       of spherical angles ...

!       begin with xy-plane vector phi

        abmax  = max( abs(dir_tab(1,i_atom)), abs(dir_tab(2,i_atom)) )

        if (abmax .gt. 0.0d+0) then

          ab = sqrt( dir_tab(1,i_atom)**2.d0 + dir_tab(2,i_atom)**2.d0 )

!         cos(phi)
          trigonom_tab(4,i_atom) = dir_tab(1,i_atom) / ab

!         sin(phi)
          trigonom_tab(3,i_atom) = dir_tab(2,i_atom) / ab

        else
!         we're on the z axis - phi not defined.

           trigonom_tab(4,i_atom) = 1.0d+0
           trigonom_tab(3,i_atom) = 0.0d+0
           ab = 0.d0

        end if

!       next, out-of-xy-plane vector theta

        abcmax = max( abmax, abs(dir_tab(3,i_atom)) )

        if (abcmax .gt. 0.0d+0) then

          abc = sqrt( ab**2.0d0 + dir_tab(3,i_atom)**2.d0 )

!         cos(theta)
          trigonom_tab(2,i_atom) = dir_tab(3,i_atom) / abc

!         sin(theta)
          trigonom_tab(1,i_atom) = ab / abc

        else
!         this piece of code should never trigger, but place us on the z axis anyway

          trigonom_tab(2,i_atom) = 1.0d+0
          trigonom_tab(1,i_atom) = 0.0d+0

        endif

      enddo

!  that's it

      end subroutine tab_trigonom_v2
!----------------------------------------------------------------------
!******
