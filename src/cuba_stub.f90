!****s* FHI-aims/cuba_stub
!  NAME
!   divonne
!  SYNOPSIS
subroutine divonne(ndim, ncomp, integrand, epsrel, epsabs, verbose, mineval, maxeval,	 &
        key1, key2, key3, maxpass, border, maxchisq, mindeviation,		 &
        ngiven, ldxgiven, tken, nextra, tken2,		 &
        nregions, neval, fail, integral, error, prob)

! PURPOSE
!   This is an Cuba stub. It ensures that AIMS can be 
!   compiled without cuba library (for Monte Carlo integration in vdw density functiona).
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
  
  use localorb_io

  implicit none

  integer          :: ndim, ncomp, neval, key1, key2, key3, maxpass, ngiven, ldxgiven, &
                      nregions, fail, verbose, tken, tken2
  double precision :: integral, error, prob, maxchisq, mindeviation, epsrel, epsabs, &
                      mineval, maxeval, integrand, border, nextra
  character*160    :: info_str

  write(info_str,'(1X,A,A)') '* You have called a cuba stub routine. ', &
       'Please check conrol.in and the libraries you linked.'
  call localorb_info(info_str, use_unit, '(A)')
  stop

  
end subroutine divonne
!******

subroutine lldivonne(ndim, ncomp, integrand, epsrel, epsabs, verbose, mineval, maxeval,	 &
        key1, key2, key3, maxpass, border, maxchisq, mindeviation,		 &
        ngiven, ldxgiven, tken, nextra, tken2,		 &
        nregions, neval, fail, integral, error, prob)

! PURPOSE
!   This is an Cuba stub. It ensures that AIMS can be 
!   compiled without cuba library (for Monte Carlo integration in vdw density functiona).
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
  
  use localorb_io

  implicit none

  integer          :: ndim, ncomp, neval, key1, key2, key3, maxpass, ngiven, ldxgiven, &
                      nregions, fail, verbose, tken, tken2
  double precision :: integral, error, prob, maxchisq, mindeviation, epsrel, epsabs, &
                      mineval, maxeval, integrand, border, nextra
  character*160    :: info_str

  write(info_str,'(1X,A,A)') '* You have called a cuba stub routine. ', &
       'Please check conrol.in and the libraries you linked.'
  call localorb_info(info_str, use_unit, '(A)')
  stop

  
end subroutine lldivonne
!******

