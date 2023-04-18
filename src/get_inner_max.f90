!****f* FHI-aims/get_inner_max
!  NAME
!   get_inner_max
!  SYNOPSIS

real*8 function get_inner_max(wave,r_grid,n_grid)

!  PURPOSE
!  Function get_inner_max
!  returns the real x value of the innermost maximum of a wave
!  function on a logarithmic grid.
!
!  USES
!  ARGUMENTS

  integer n_grid
  real*8 wave(n_grid)
  real*8 r_grid(n_grid)

!  INPUTS
!  o wave -- basis function
!  o r_grid -- radial grid
!  o n_grid -- number of radial grid points
!
!  OUTPUT
!  o get_inner_max -- the innermost maximum of a wave function
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

  !     counters

  integer i_grid

  !  begin work

  i_grid = 1


  do while ( (i_grid.lt.n_grid)  .and. ( abs(wave(i_grid+1)) .ge. abs(wave(i_grid)) ))
     i_grid = i_grid+1
  enddo


  if (i_grid.eq.n_grid) then
     ! no maximum, i.e. no reasonable wave function - reject.
     get_inner_max = 0.
  else
     ! return r value of innermost maximum
     get_inner_max = r_grid(i_grid)
  end if

  !        write(use_unit,*) "| Inner max.: ", get_inner_max

  return
end function get_inner_max
!******	
