!****h* FHI-aims/constants
!  NAME
!    constants
!  SYNOPSIS

      module constants

! PURPOSE
!
!  Fundamental constants:
!  o  pi                : Pi
!  o  pi4               :  4*Pi
!  o  pi4_inv           : 1/(4 Pi)
!  o  sqrt_pi           : sqrt(Pi)
!  o  pi2sqrt_inv       : 1 / sqrt(2 * Pi)
!  o  pisqrt_inv        : 1 / sqrt(Pi)
!  o  pisq3             : Pi^2*3
!  o  third             : 1/3
!  o  one_over_sqrt2    : 1 /sqrt(2)
!  o  euler_gamma       : Euler-Mascheroni constant (gamma ~ 0.57721)
!  o  img_unit          : complex imaginar unit i
!  o  const_rs          : constant prefactor for Wigner-Seitz radius
!  o  radians           : radians in degree
!  o  light_speed       : Speed of light in atomic units
!  o  light_speed_sq    :  sqrt(Speed of light in atomic units)
!  o  bohr              : H radius in Hartrees
!  o  hartree           : 1 Hartree in eV
!  o  hartree_over_bohr : Hartree/bohr
!  o  boltzmann_kB      : Boltzmann kb
!  o  MD_KE_factor      : molecular dynamics kinetic energy factor  
!                         (amu*bohr^2/ps^2)
!  o  NUM_ZERO          : numerical Zero
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




      implicit none


      real*8 pi
      real*8 pi4
      real*8 pi4_inv
      real*8 sqrt_pi
      real*8 pi2sqrt_inv
      real*8 pisqrt_inv
      real*8 pisq3
      real*8 third
      real*8 one_over_sqrt2
      real*8 euler_gamma
      complex*16 img_unit
      real*8 const_rs
      real*8 radians
      real*8 light_speed
      real*8 light_speed_sq
      real*8 bohr
      real*8 hartree
      real*8 hartree_over_bohr
      real*8 bohr_over_hartree
      real*8 boltzmann_kB
      real*8 MD_KE_factor
      real*8 MD_Q_conversion
      real*8 hbar_eV_ps
      real*8 giga_pascal
      real*8 mass_atomic_unit  ! XZL: added for PIMD
      real*8,dimension(0:200):: LMtable

      complex*16, public, parameter :: C_ONE      = (1.0d0, 0.0d0)
      complex*16, public, parameter :: C_IMAG_ONE = (0.0d0, 1.0d0)
      complex*16, public, parameter :: C_ZERO = (0.0d0, 0.0d0)
      real*8, public, parameter :: M_ZERO  =  0.0d0
      real*8, public, parameter :: M_HALF  =  0.5d0
      real*8, public, parameter :: M_ONE   =  1.0d0
      real*8, public, parameter :: M_TWO   =  2.0d0
      real*8, public, parameter :: M_THREE =  3.0d0
      real*8, public, parameter :: M_FOUR  =  4.0d0
      real*8, public, parameter :: M_FIVE  =  5.0d0
      real*8, public, parameter :: M_SIX   =  6.0d0
      real*8, public, parameter :: M_SEVEN =  7.0d0
      real*8, public, parameter :: M_EIGHT =  8.0d0
      real*8, public, parameter :: M_NINE  =  9.0d0
      real*8, public, parameter :: M_TEN   = 10.0d0
      real*8, public, parameter :: NUM_ZERO= 1e-20

      parameter ( pi      = 3.14159265358979323846d0 )
      parameter ( pi4     = 12.56637061435917295376d0 )
      parameter ( pi4_inv = 1. / pi4 )
      parameter ( sqrt_pi = 1.77245385090551602729d0 )
      parameter ( pisq3   = 29.60881320326807585607d0 )
      parameter ( third   = 0.333333333333333333d0 )
      parameter ( one_over_sqrt2 = 0.70710678118654752440d0 )
      parameter ( img_unit = (0.0d0,1.0d0) )
      parameter ( const_rs = 1.91915829267751281d0 )
      parameter ( pi2sqrt_inv = 0.39894228040143267794d0 )
      parameter ( pisqrt_inv =  0.56418958354775628695d0 )
      parameter ( euler_gamma = 0.57721566490153286060d0 )
      parameter ( radians = 360.d0 / (2.d0 * pi) )
      parameter ( mass_atomic_unit = 1.8226806d+3) ! XZL: added for PIMD

!  for scalar relativity always:     parameter ( light_speed = 100000000d0   )
      parameter ( light_speed = 137.0359895d0 )
      parameter ( light_speed_sq = light_speed*light_speed )

!  get constants from: http://physics.nist.gov/cuu/Constants/,
!                      (current as of 2002)
!  bohr: H radius in Hartrees
!        from http://physics.nist.gov/cuu/Constants/ (current as of 2002)
!  hartree: 1 Hartree in eV
!        from http://physics.nist.gov/cuu/Constants/energy.html (2004)


      parameter ( bohr    = 0.52917721d0 )
      parameter ( hartree = 27.2113845d0 )
      parameter ( hartree_over_bohr = 51.42206426d0 )
      parameter ( bohr_over_hartree = 0.019446905d0 )

      ! for MD, need Boltzmann constant (in Ha/K) from http://physics.nist.gov/cgi-bin/cuu/Value?tkev|search_for=boltzmann+constant
      parameter ( boltzmann_kB = 3.16681535d-6)

      ! this changes from kinetic energy in MD units (amu*bohr^2/ps^2) to Hartree ...
      parameter ( MD_KE_factor = 1.06657213d-6)

      ! this changes the molecular dynamics Q from [cm^-1] to amu bohr^2, via the formula 
      ! Q = MD_Q_conversion*(3 n_atoms T)/ omega^2 ; MD_Q_conversion = 1E-4 kB /(4 pi^2 amu c^2 bohr^2); all in their SI numerical values. 
      parameter ( MD_Q_conversion = 8.36818496904E1)

      ! h_bar (Planck constant over 2Pi in eV ps
      ! http://physics.nist.gov/cgi-bin/cuu/Value?hev|search_for=universal_in!
      parameter ( hbar_eV_ps = 6.58211899d-4)

      ! Conversion factor for pressure from eV/A**3 to GPa
      parameter ( giga_pascal = 160.2176565d0)


!Sqrts of integers. Used only in increment_ylm_deriv_forvdw.f90. SAG
  data LMtable /         0d0, &    
         1d0, &
         1.4142135623731d0, &
         1.73205080756888d0, &
         2d0, &
         2.23606797749979d0, &
         2.44948974278318d0, &
         2.64575131106459d0, &
         2.82842712474619d0, &
         3d0, &
         3.16227766016838d0, &
         3.3166247903554d0, &
         3.46410161513775d0, &
         3.60555127546399d0, &
         3.74165738677394d0, &
         3.87298334620742d0, &
         4d0, &
         4.12310562561766d0, &
         4.24264068711928d0, &
         4.35889894354067d0, &
         4.47213595499958d0, &
         4.58257569495584d0, &
         4.69041575982343d0, &
         4.79583152331272d0, &
         4.89897948556636d0, &
         5d0, &
         5.09901951359278d0, &
         5.19615242270663d0, &
         5.29150262212918d0, &
         5.3851648071345d0, &
         5.47722557505166d0, &
         5.56776436283002d0, &
         5.65685424949238d0, &
         5.74456264653803d0, &
         5.8309518948453d0, &
         5.91607978309962d0, &
         6d0, &
         6.08276253029822d0, &
         6.16441400296898d0, &
         6.2449979983984d0, &
         6.32455532033676d0, &
         6.40312423743285d0, &
         6.48074069840786d0, &
         6.557438524302d0, &
         6.6332495807108d0, &
         6.70820393249937d0, &
         6.78232998312527d0, &
         6.85565460040104d0, &
         6.92820323027551d0, &
         7d0, &
         7.07106781186548d0, &
         7.14142842854285d0, &
         7.21110255092798d0, &
         7.28010988928052d0, &
         7.34846922834953d0, &
         7.41619848709566d0, &
         7.48331477354788d0, &
         7.54983443527075d0, &
         7.61577310586391d0, &
         7.68114574786861d0, &
         7.74596669241483d0, &
         7.81024967590665d0, &
         7.87400787401181d0, &
         7.93725393319377d0, &
         8d0, &
         8.06225774829855d0, &
         8.12403840463596d0, &
         8.18535277187245d0, &
         8.24621125123532d0, &
         8.30662386291807d0, &
         8.36660026534076d0, &
         8.42614977317636d0, &
         8.48528137423857d0, &
         8.54400374531753d0, &
         8.60232526704263d0, &
         8.66025403784439d0, &
         8.71779788708135d0, &
         8.77496438739212d0, &
         8.83176086632785d0, &
         8.88819441731559d0, &
         8.94427190999916d0, &
         9d0, &
         9.05538513813742d0, &
         9.1104335791443d0, &
         9.16515138991168d0, &
         9.21954445729289d0, &
         9.2736184954957d0, &
         9.32737905308882d0, &
         9.38083151964686d0, &
         9.4339811320566d0, &
         9.48683298050514d0, &
         9.53939201416946d0, &
         9.59166304662544d0, &
         9.64365076099295d0, &
         9.69535971483266d0, &
         9.74679434480896d0, &
         9.79795897113271d0, &
         9.8488578017961d0, &
         9.89949493661167d0, &
         9.9498743710662d0, &
         10d0, &
         10.0498756211209d0, &
         10.0995049383621d0, &
         10.1488915650922d0, &
         10.1980390271856d0, &
         10.2469507659596d0, &
         10.295630140987d0, &
         10.3440804327886d0, &
         10.3923048454133d0, &
         10.4403065089106d0, &
         10.4880884817015d0, &
         10.5356537528527d0, &
         10.5830052442584d0, &
         10.6301458127346d0, &
         10.6770782520313d0, &
         10.7238052947636d0, &
         10.770329614269d0, &
         10.816653826392d0, &
         10.8627804912002d0, &
         10.9087121146357d0, &
         10.9544511501033d0, &
         11d0, &
         11.0453610171873d0, &
         11.0905365064094d0, &
         11.13552872566d0, &
         11.1803398874989d0, &
         11.2249721603218d0, &
         11.2694276695846d0, &
         11.3137084989848d0, &
         11.3578166916005d0, &
         11.4017542509914d0, &
         11.4455231422596d0, &
         11.4891252930761d0, &
         11.5325625946708d0, &
         11.5758369027902d0, &
         11.6189500386223d0, &
         11.6619037896906d0, &
         11.7046999107196d0, &
         11.7473401244707d0, &
         11.7898261225516d0, &
         11.8321595661992d0, &
         11.8743420870379d0, &
         11.916375287813d0, &
         11.9582607431014d0, &
         12d0, &
         12.0415945787923d0, &
         12.0830459735946d0, &
         12.1243556529821d0, &
         12.1655250605964d0, &
         12.2065556157337d0, &
         12.2474487139159d0, &
         12.2882057274445d0, &
         12.328828005938d0, &
         12.369316876853d0, &
         12.4096736459909d0, &
         12.4498995979887d0, &
         12.4899959967968d0, &
         12.5299640861417d0, &
         12.5698050899765d0, &
         12.6095202129185d0, &
         12.6491106406735d0, &
         12.6885775404495d0, &
         12.7279220613579d0, &
         12.7671453348037d0, &
         12.8062484748657d0, &
         12.8452325786651d0, &
         12.8840987267251d0, &
         12.9228479833201d0, &
         12.9614813968157d0, &
         13d0, &
         13.0384048104053d0, &
         13.076696830622d0, &
         13.114877048604d0, &
         13.1529464379659d0, &
         13.1909059582729d0, &
         13.228756555323d0, &
         13.2664991614216d0, &
         13.3041346956501d0, &
         13.3416640641263d0, &
         13.3790881602597d0, &
         13.4164078649987d0, &
         13.4536240470737d0, &
         13.490737563232d0, &
         13.5277492584687d0, &
         13.5646599662505d0, &
         13.6014705087354d0, &
         13.6381816969859d0, &
         13.6747943311773d0, &
         13.7113092008021d0, &
         13.7477270848675d0, &
         13.7840487520902d0, &
         13.8202749610853d0, &
         13.856406460551d0, &
         13.8924439894498d0, &
         13.9283882771841d0, &
         13.9642400437689d0, &
         14d0, &
         14.0356688476182d0, &
         14.0712472794703d0, &
         14.1067359796659d0, &
         14.142135623731d0 /













      end  module constants
!******
