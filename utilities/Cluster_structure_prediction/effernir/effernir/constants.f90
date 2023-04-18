!  fundamental constants for use in joana

!  pi:   well, Pi ...
!  pi4:  4*Pi
!  pi4_inv: 1/(4 pi)
!  pisq3: Pi^2*3
!  third: one third
!  one_over_sqrt2 : might be better as a constant, saves a number of calls near the ylm fns
!  img_unit : complex imaginar unit i
!  const_rs : constant prefactor for Wigner-Seitz radius

      real*8 pi
      real*8 pi4
      real*8 pi4_inv
      real*8 pisq3
      real*8 third
      real*8 one_over_sqrt2
      complex*16 img_unit
      real*8 const_rs

      parameter ( pi      = 3.14159265d0 )
      parameter ( pi4     = 12.5663706143592d0 )
      parameter ( pi4_inv = 1. / pi4 )
      parameter ( pisq3   = .296088132032680740d2 )
      parameter ( third   = .333333333333333333d0 )
      parameter ( one_over_sqrt2 = 0.70710678d0 )
      parameter ( img_unit = (0.0d0, 1.0d0) )
      parameter ( const_rs = 1.91915829267751281d0 )

!  get constants from: http://physics.nist.gov/cuu/Constants/,
!                      (current as of 2002)
!  bohr: H radius in Hartrees
!        from http://physics.nist.gov/cuu/Constants/ (current as of 2002)
!  hartree: 1 Hartree in eV
!        from http://physics.nist.gov/cuu/Constants/energy.html (2004)
!  boltzmann_hartree: k_boltzmann in hartree / K
!        from http://physics.nist.gov/cgi-bin/cuu
!  boltzmann_eV: k_boltzmann in eV / K
!        from http://physics.nist.gov/cgi-bin/cuu
 
      real*8 :: bohr
      real*8 :: hartree
      real*8 :: boltzmann_hartree
      real*8 :: boltzmann_eV
      real*8 :: inv_bohr
      real*8 :: inv_hartree

      parameter ( bohr          = 0.52917721d0 ) 
      parameter ( hartree       = 27.2113845d0 )

      parameter ( boltzmann_hartree = 8.617343 * 3.67493245e-7)
      parameter ( boltzmann_eV = 8.617343e-5)
      
      parameter ( inv_bohr = 1.d0 / bohr)
      parameter ( inv_hartree = 1.d0 / hartree)
