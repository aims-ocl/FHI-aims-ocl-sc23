!****s* FHI-aims/radial_fourier_hse
!  NAME
!    radial_fourier_hse
!  SYNOPSIS

real*8 function radial_fourier_hse(k, omega)

  !  PURPOSE
  !
  !    For a given k, return the factor to be multiplied to the
  !    radial function in k-space in order to obtain the Fourier transform
  !    of its HSE potential.
  !
  !  USES

  use constants
  use prodbas
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: k, omega

  !  INPUTS
  !    o k -- Magnitude of k vector
  !    o omega -- Screening parameter
  !  OUTPUTS
  !    o radial_fourier_hse -- Factor to be multiplied to function.
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  real*8 :: exp_const
  character(*), parameter :: func = 'radial_fourier_hse'

  exp_const = 0.25d0 / omega**2
  if (k**2 * exp_const < 1d-5) then
     ! Beginning of Taylor series...
      if (ovlp_type_bare_or_hse_coul .ne. OVLP_TYPE_LR) then
    		 radial_fourier_hse = 4*pi * (exp_const - exp_const**2/2 * k**2 &
     		 &                               + exp_const**3/6 * k**4)
      else
      	radial_fourier_hse = 4*pi * (1.d0 / k**2 - (exp_const - exp_const**2/2 * k**2 &
     		 &                               + exp_const**3/6 * k**4))
      end if
  else
  	 if (ovlp_type_bare_or_hse_coul .ne. OVLP_TYPE_LR) then
	     radial_fourier_hse = 4*pi / k**2 * (1.d0 - exp(- k**2 * exp_const))
	 else
	 	  radial_fourier_hse = 4*pi / k**2 * (exp(- k**2 * exp_const))
	 end if
  end if

end function radial_fourier_hse
!******
