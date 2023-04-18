C---------------------------------------------------------------------------
C  VB 2005, FHI: Modified version of original pbe code, to generate and
C  extract only the partial derivatives of the functional which are needed 
C  for matrix elements later:
C
C  i.e. we do not construct the local potential itself, because that would
C  would require the Hessian of the charge density.
C
C---------------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
c
      SUBROUTINE EXCH_HSE_derivs
     +     (hybrid_coef,hse_omega_pbe, rho, squared_grad_rho,
     +     en_density_x, x_density_deriv, x_gradient_deriv)
c
c---------------------------------------------------------------------
c
c  subroutine to get the short range pbe part   
c
c  SG - subroutine based on:      
c
c---------------------------------------------------------------------     
C
C  VB - Subroutine based on:
c
C----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C
C  INPUT: 
C
C  rho : DENSITY
C  squared_grad_rho: |grad(rho)|^2
C
C  OUTPUT:
C
C  en_density_x :     exchange energy density e(x)
C  x_density_deriv :  partial derivative d(rho*e)/d(rho)
C  x_gradient_deriv : partial derivative d(rho*e)/d[|grad(rho)|^2]
C
C  Derived quantities:
C
C  exunif: would-be exchange contribution of the uniform electron gas
C          for given density
C  scale:  scaling factor (2*KF*rho)^2 see below
C  SSQ:    squared scaled gradient [ ABS(GRAD rho)/(2*KF*rho) ]^2, 
C          where kf=(3 pi^2 rho)^(1/3)
C          dubbed s*s below
C  P0:     1 + ul * SSQ
C  FxPBE:  PBE enhancement factor (gradient correction)
C  Fs:     Partial derivative of FxPBE by squared scaled gradient
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submitted to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c       e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c   e_x[PBE]=e_x[unif]*FxPBE(s)
c   FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none
c
C     declare constants
c
      real*8 :: thrd, thrd2, thrd4, thrd8, ax, um, uk, ul,pi
c
      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c
c     INPUT 
c
      real*8 :: hybrid_coef
      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: hse_omega_pbe
c
c     OUTPUT
c
      real*8 :: en_density_x
      real*8 :: x_density_deriv
      real*8 :: x_gradient_deriv

c
c     VARIABLES THAT ARE NEEDE DURING CALCULRATION
c     
      real*8 :: exunif
      real*8 :: scale
      real*8 :: scale_hse
      real*8 :: en_density_x_pbe
      real*8 :: x_density_deriv_pbe
      real*8 :: x_gradient_deriv_pbe
      real*8 :: SSQ
      real*8 :: scaledgradient
c
c     HSE  stuff
c
      real*8 :: FxHSE
      real*8 :: dFxHSE_ds
      real*8 :: dFxHSE_drho
c
c     PBE stuff
c
      real*8 :: FxPBE
      real*8 :: dFxPBE_ds
      real*8 :: dFxPBE_drho
c
c     construct LDA exchange energy density
c
      exunif = AX*rho**THRD
c
c     now begin hse preparations
c
      scale_hse = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale_hse*(rho**thrd8)
      SSQ = squared_grad_rho / scale
c      
      scaledgradient =dsqrt(ssq)
c
c     first: get the pbe-relevant terms
c      
      call HSEFx_partial(rho,scaledgradient,0.0d0,
     +     FxPBE,dFxPBE_drho,dFxPBE_ds)
c         
      en_density_x_pbe = exunif*FxPBE
c      
      x_density_deriv_pbe = thrd4*en_density_x_pbe + exunif * rho * 
     +     dFxPBE_drho  -  thrd4 * Ax*dFxPBE_ds*dsqrt(squared_grad_rho)* 
     +     0.5d0  / (rho *(3.d0*pi*pi)**thrd)
c      
      x_gradient_deriv_pbe = Ax * dFxPBE_ds * 0.25d0 
     +     / (dsqrt(squared_grad_rho * scale_hse))
c
c     second: get the hse-relevant terms
c         
      call HSEFx_partial(rho,scaledgradient,hse_omega_pbe,
     +     FxHSE,dFxHSE_drho,dFxHSE_ds)
c
      en_density_x = exunif*FxHSE
c      
      x_density_deriv = thrd4*en_density_x + exunif * rho * 
     +     dFxHSE_drho 
     +     -  thrd4 * Ax * dFxHSE_ds * dsqrt(squared_grad_rho)* 
     +     0.5d0  / (rho *(3.d0*pi*pi)**thrd)
c
      x_gradient_deriv = Ax * dFxHSE_ds * 0.25d0 
     +     / (dsqrt(squared_grad_rho * scale_hse))
c
c
c    finally: add pbe- and hse-part together via hybrid_coefficient
c
      en_density_x = en_density_x_pbe - hybrid_coef * en_density_x
c
      x_density_deriv = x_density_deriv_pbe -hybrid_coef*x_density_deriv
c
      x_gradient_deriv=x_gradient_deriv_pbe-hybrid_coef*x_gradient_deriv
c     
      RETURN
      END
      
c#####################################################################
c---------------------------------------------------------------------

      SUBROUTINE exch_lc_wpbeh_derivs
     +     (hybrid_coef,lc_wpbeh_omega, rho, squared_grad_rho,
     +     en_density_x, x_density_deriv, x_gradient_deriv,
     +     lc_dielectric_constant)
c
c---------------------------------------------------------------------
c
c  subroutine to get the short range pbe part   
c
c  SG - subroutine based on:      
c
c---------------------------------------------------------------------     
C
C  VB - Subroutine based on:
c
C----------------------------------------------------------------------
C  
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C
C  INPUT: 
C
C  rho : DENSITY
C  squared_grad_rho: |grad(rho)|^2
C
C  OUTPUT:
C
C  en_density_x :     exchange energy density e(x)
C  x_density_deriv :  partial derivative d(rho*e)/d(rho)
C  x_gradient_deriv : partial derivative d(rho*e)/d[|grad(rho)|^2]
C
C  Derived quantities:
C
C  exunif: would-be exchange contribution of the uniform electron gas
C          for given density
C  scale:  scaling factor (2*KF*rho)^2 see below
C  SSQ:    squared scaled gradient [ ABS(GRAD rho)/(2*KF*rho) ]^2, 
C          where kf=(3 pi^2 rho)^(1/3)
C          dubbed s*s below
C  P0:     1 + ul * SSQ
C  FxPBE:  PBE enhancement factor (gradient correction)
C  Fs:     Partial derivative of FxPBE by squared scaled gradient
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c       e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c   e_x[PBE]=e_x[unif]*FxPBE(s)
c   FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none
c
C     declare constants
c
      real*8 :: thrd, thrd2, thrd4, thrd8, ax, um, uk, ul,pi
c
      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c
c     INPUT 
c
      real*8 :: hybrid_coef
      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: lc_wpbeh_omega
c     This was added 26.04.2016. Thats why it is not in order in
c     the argument list
      real*8 :: lc_dielectric_constant
c
c     OUTPUT
c
      real*8 :: en_density_x
      real*8 :: x_density_deriv
      real*8 :: x_gradient_deriv

c
c     VARIABLES THAT ARE NEEDE DURING CALCULRATION
c     
      real*8 :: exunif
      real*8 :: scale
      real*8 :: scale_hse
      real*8 :: en_density_x_hf
      real*8 :: x_density_deriv_hf
      real*8 :: x_gradient_deriv_hf
      real*8 :: en_density_x_pbe
      real*8 :: x_density_deriv_pbe
      real*8 :: x_gradient_deriv_pbe
      real*8 :: SSQ
      real*8 :: scaledgradient
c
c     SR  stuff
c
      real*8 :: FxSR
      real*8 :: dFxSR_ds
      real*8 :: dFxSR_drho
c
c     HF stuff
c
      real*8 :: FxHF
      real*8 :: dFxHF_ds
      real*8 :: dFxHF_drho
      
c     PBE stuff
c
      real*8 :: FxPBE
      real*8 :: dFxPBE_ds
      real*8 :: dFxPBE_drho
c
c     construct LDA exchange energy density
c
      exunif = AX*rho**THRD
c
c     now begin hse preparations
c
      scale_hse = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale_hse*(rho**thrd8)
      SSQ = squared_grad_rho / scale
c      
      scaledgradient =dsqrt(ssq)
c      
c     get the sr-relevant terms for wPBE (from HSE subroutine)
c     
c         
      call HSEFx_partial(rho,scaledgradient,lc_wpbeh_omega,
     +     FxSR,dFxSR_drho,dFxSR_ds)
c
      en_density_x = exunif*FxSR
c      
      x_density_deriv = thrd4*en_density_x + exunif * rho * 
     +     dFxSR_drho 
     +     -  thrd4 * Ax * dFxSR_ds * dsqrt(squared_grad_rho)* 
     +     0.5d0  / (rho *(3.d0*pi*pi)**thrd)
c
      x_gradient_deriv = Ax * dFxSR_ds * 0.25d0 
     +     / (dsqrt(squared_grad_rho * scale_hse))
c

      IF (lc_dielectric_constant .gt. 1.0d0) THEN
          call HSEFx_partial(rho,scaledgradient,0.0d0,
     +     FxPBE,dFxPBE_drho,dFxPBE_ds)
c         
          en_density_x_pbe = exunif*FxPBE
c      
          x_density_deriv_pbe = thrd4*en_density_x_pbe + exunif * rho * 
     +     dFxPBE_drho  -  thrd4 * Ax*dFxPBE_ds*dsqrt(squared_grad_rho)* 
     +     0.5d0  / (rho *(3.d0*pi*pi)**thrd)
c      
          x_gradient_deriv_pbe = Ax * dFxPBE_ds * 0.25d0 
     +     / (dsqrt(squared_grad_rho * scale_hse))
         
          en_density_x =  (1.d0/lc_dielectric_constant - 
     +     hybrid_coef)*
     +     en_density_x + (1.d0 - 1.d0/lc_dielectric_constant)*
     +     en_density_x_pbe
c
          x_density_deriv =  (1.d0/lc_dielectric_constant - 
     +     hybrid_coef)*
     +     x_density_deriv + (1.d0 - 1.d0/lc_dielectric_constant)*
     +     x_density_deriv_pbe
c
          x_gradient_deriv =  (1.d0/lc_dielectric_constant - 
     +     hybrid_coef)* 
     +     x_gradient_deriv + (1.d0 - 1.d0/lc_dielectric_constant)*
     +     x_gradient_deriv_pbe
      ELSE
c
c    finally: add together via hybrid_coefficient
c
      en_density_x =  (1.d0 - hybrid_coef) *  en_density_x
c
      x_density_deriv =  (1.d0 - hybrid_coef) * x_density_deriv
c
      x_gradient_deriv =  (1.d0 - hybrid_coef) * x_gradient_deriv
c     
      END IF
      RETURN
      END

c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPBE_derivs
     +( rho, squared_grad_rho,
     +  en_density_x, x_density_deriv, x_gradient_deriv )
C
C  VB - Subroutine based on:
c
C----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C
C  INPUT: 
C
C  rho : DENSITY
C  squared_grad_rho: |grad(rho)|^2
C
C  OUTPUT:
C
C  en_density_x :     exchange energy density e(x)
C  x_density_deriv :  partial derivative d(rho*e)/d(rho)
C  x_gradient_deriv : partial derivative d(rho*e)/d[|grad(rho)|^2]
C
C  Derived quantities:
C
C  exunif: would-be exchange contribution of the uniform electron gas
C          for given density
C  scale:  scaling factor (2*KF*rho)^2 see below
C  SSQ:    squared scaled gradient [ ABS(GRAD rho)/(2*KF*rho) ]^2, 
C          where kf=(3 pi^2 rho)^(1/3)
C          dubbed s*s below
C  P0:     1 + ul * SSQ
C  FxPBE:  PBE enhancement factor (gradient correction)
C  Fs:     Partial derivative of FxPBE by squared scaled gradient
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submitted to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c       e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c   e_x[PBE]=e_x[unif]*FxPBE(s)
c   FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none

C     declare constants

      real*8 :: thrd, thrd2, thrd4, thrd8, pi, ax, um, uk, ul

      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)

c     declare variables

      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: en_density_x
      real*8 :: x_density_deriv
      real*8 :: x_gradient_deriv

      real*8 :: exunif
      real*8 :: scale
      real*8 :: SSQ
      real*8 :: P0
      real*8 :: FxPBE
      real*8 :: Fs

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density

      exunif = AX*rho**THRD

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      scale = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale*(rho**thrd8)
C      write(*,*) "PPE:",rho,squared_grad_rho
      SSQ = squared_grad_rho / scale
      P0=1.d0+ul*SSQ
      FxPBE = 1d0+uk-uk/P0
      
      en_density_x = exunif*FxPBE
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE partial derivatives of rho*exunif(rho)*FxPBE(ssq) :
c  first, find derivative of Fx w.r.t ssq:
c  Fs = d FxPBE/ d(ssq)
      Fs=uk*ul/(P0*P0)
C  next, multiply with derivatives d(ssq)/d(rho) and d(ssq)/d(squared_grad_rho)
C  to obtain desired partial derivatives;
C  the definition of ssq gives
C  d(ssq)/d(rho) = - 8/3 ssq/rho
C  d(ssq)/d(squared_grad_rho) = 1/ (2*KF*rho)^2 = 1/scale

      x_density_deriv = thrd4*en_density_x - exunif*Fs*thrd8*ssq

      x_gradient_deriv = rho * exunif * Fs / scale

      RETURN
      END
      
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORPBE_derivs
     +( rho, squared_grad_rho, spin_pol, 
     +  en_density_c, c_density_deriv, c_gradient_deriv, c_spin_deriv
     +)
C
C  rho:               Density
C  squared_grad_rho : |grad(rho)|^2
C  spin_pol:          zeta = (rhoup-rhodn)/rho
C  en_density_c:      Correlation energy density e(x) = e_lda(x) + H(x)
C                     Really e = e( rho, squared_grad_rho, zeta )
C  c_density_deriv:   partial derivative d(rho*e)/d(rho)
C  c_gradient_deriv:  partial derivative d(rho*e)/d(|grad(rho)|^2)
C  c_spin_deriv:      partial derivative (1/rho)*d(rho*e)/d(zeta)
C
C                     Memo: c_spin_deriv formally contains an 1/rho part
C                           which is needed to form the XC potential.
C                           This means that the full XC potential could be obtained 
C                           simply by
C
C                           vcup = c_density_deriv - (zeta-1)*c_spin_deriv - 2 * div ( c_gradient_deriv * grad(rho) )
C                           vcdn = c_density_deriv - (zeta+1)*c_spin_deriv - 2 * div ( c_gradient_deriv * grad(rho) )
C
C  VB - Subroutine based on:
C
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  Original input (that's when the local potential _was_ computed)
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none

ctest
c      integer test
c      real*8 Q0,Q1,Q2
ctest end

C  constants go first

C  original constants from PBE code

c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at 
c          |zeta|=1.

      real*8 thrd, thrdm, thrd2, sixthm, thrd4

      real*8 gam, fzz, gamma, bet, delt, eta

      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)

C     constants from Martin Fuchs' wrapper

      real*8 crs, pi, pisq3
      parameter (crs=1.91915829267751281d0)
      parameter (pi=.314159265358979312d1)
      parameter (pisq3=.296088132032680740d2)

C  variables

C     input

      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: spin_pol

C     output

      real*8 :: en_density_c 
      real*8 :: c_density_deriv
      real*8 :: c_gradient_deriv
      real*8 :: c_spin_deriv 

C     local variables

C     fk: Fermi momentum k_f
C     rs: Seitz radius
C     rtrs: sqrt(rs)
C     scale: scaling factor for t^2

      real*8 :: fk
      real*8 :: rs
      real*8 :: rtrs
      real*8 :: scale

C     for uniform electron gas term:
c     EU=unpolarized LSD correlation energy
c     EURS=dEU/drs
c     EP=fully polarized LSD correlation energy
c     EPRS=dEP/drs
c     ALFM=-spin stiffness, [c](3).
c     ALFRSM=-dalpha/drs

C     alfc: spin stiffness -alfm
C     z4 :  zeta^4 (zeta=spin polarisation spin_pol)
c     F=spin-scaling factor from [c](9).
c     ecunif: correlation energy density of the uniform electron gas.
c     ECRS = dEc/drs [c](A2)
c     ECZET = dEc/dzeta [c](A3)
c     FZ = dF/dzeta [c](A4)

      real*8 :: eu, eurs
      real*8 :: ep, eprs
      real*8 :: alfm, alfrsm

      real*8 :: alfc, z4, f

      real*8 :: ecunif

      real*8 :: ECRS 
      real*8 :: FZ 
      real*8 :: ECZET

C     for correlation correction H( rho, squared_grad_rho, zeta)
c
c     these quantities are all intermediates of the orig PBE paper, but they are 
c     strangely renamed in the original PBE code. 
c
c     G=phi(zeta), given after [a](3)
c     G3=phi(zeta)**3
C     PON = -ecunif/(phi(zeta)**3 * gam*e^2/a0)
c     B=A of [a](8)
C     B_sq = B*B
C
C     T_sq : squared scaled gradient t of PBE paper
C
C     T4 : t^4
C     Q4 : 1+At^2 - see PBE paper
C     Q5 : 1 + A t^2 + A^2 t^4 - see PBE paper
C
C     ec_gradient_correction : gradient correction to correlation energy - H of PBE paper
C
      real*8 :: G, G3, PON, B
      real*8 :: B_sq

      real*8 :: T_sq
      real*8 :: T4
      real*8 :: Q4, Q5

      real*8 :: ec_gradient_correction

      real*8 :: G4,T6,RSTHRD,GZ
      real*8 :: FAC, BG, BEC, Q8, Q9
      real*8 :: hB, hRS, hZ, hT2

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c determine derived quantities needed later

      fk=(pisq3*rho)**thrd
      rs=crs/fk

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c construct ec, using [c](8)

      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,rtRS,ALFM,ALFRSM)

      ALFC = -ALFM
      Z4 = spin_pol**4
      F=((1.D0+spin_pol)**THRD4+(1.D0-spin_pol)**THRD4-2.D0)/GAM

ctest
c      if (test.eq.1) then
c        write(6,*) "rtrs = ", rtrs
c        write(6,*) "eu = ", eu
c        Q0 = -2.D0*0.0310907D0*(1.D0+0.21370D0*rtrs*rtrs)
c        Q1 = 2.D0*0.0310907D0*rtrs*(7.5957D0+rtrs*
c     +       (3.5876D0+rtrs*(1.6382D0+0.49294D0*rtrs)))
c
c        Q2 = DLOG(1.D0+1.D0/Q1)
c
c        write(6,*) "Q0 = ", Q0
c        write(6,*) "Q1 = ", Q1
c        write(6,*) "Q2(log) = ", Q2
c
c        Q2 = 1.D0/Q1
c        write(6,*) "Q2(taylor) = ", Q2
c
cc        GG = Q0*Q2
c      end if
ctest end

      ecunif = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ

c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  Next, calculate partial derivatives of LDA energy density ...

c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)

      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ

C     add uniform part of partial derivative of correlation energy by density
      c_density_deriv = ecunif - rs*ecrs*thrd

c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)

      FZ = THRD4*((1.D0+spin_pol)**THRD-(1.D0-spin_pol)**THRD)/GAM
      ECZET = 4.D0*(spin_pol**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)

c     add uniform part of partial derivative of correlation energy by zeta
      c_spin_deriv = eczet

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)

      G=((1.d0+spin_pol)**thrd2+(1.d0-spin_pol)**thrd2)/2.d0
      G3 = G**3
      PON=-ecunif/(G3*gamma)

ctest
c      if (test) then
c        write(6,*) "ecunif = ", ecunif
c        write(6,*) "PON    = ", PON
c      end if
ctest end

cvb - this is a fix for a potential floating point exception if dexp is numerically unstable
      if (dabs(pon).gt.1.d-8) then
        ! correct expression
        B = DELT/(DEXP(PON)-1.D0)
      else
        ! use Taylor limit
        B = PON
      end if
      B_sq = B*B

      scale = 16.d0*g*g*rho*rho*fk/pi

      T_sq = squared_grad_rho / scale

ctest
c      if (test) then
c        write(6,*) "B    = ", B
c        write(6,*) "B_sq = ", B_sq
c        write(6,*) "T_sq = ", T_sq
c      end if
ctest end

      T4 = T_sq*T_sq
      Q4 = 1.D0+B*T_sq
      Q5 = 1.D0+B*T_sq+B_sq*T4

ctest
c      if (test) then
c        write(6,*) "T4 = ", T4
c        write(6,*) "Q4 = ", Q4
c        write(6,*) "Q5 = ", Q5
c      end if
ctest end

C     this is H(rs,t,zeta)
      ec_gradient_correction = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T_sq/Q5)

      en_density_c = ecunif + ec_gradient_correction

c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE partial derivatives.
C
C We want:
C 1) d(rho*H)/d(rho)
C 2) (1/rho)*d(rho*H)/d(zeta)
C 3) d(rho*H)/d(|grad(rho)|^2)
C
C Steps:
C 1) d(rho*H)/d(rho) = H + rho*dH/d(rho)
C    dH/d(rho) = dH/d(r_s)*d(r_s)/d(rho) + dH/d(t^2)*d(t^2)/d(rho)
C    d(r_s)/d(rho)=-1/3*r_s/rho
C    d(t^2)/d(rho)=-7/3*t^2/rho
C    Hence:
C    d(rho*H) = H - 1/3*r_s*dH/d(r_s) - 7/3*t^2*dH/d(t^2)
C
C 2) (1/rho)*d(rho*H)/d(zeta) = dH/d(zeta) = rho * dH/d(phi) * d(phi)/d(zeta)
C
C    Notice that zeta appears here in two different functions, technically:
C    We want d[H(rho,|grad(rho)|^2,zeta]/d(zeta) from
C    d[h(r_s,zeta,t^2)] where t^2=t^2(phi(zeta))
C    Must apply chain rule to get
C    dH/d(zeta) = dh/d(zeta) + dh/d(t^2)*d(t^2)/d(phi)*d(phi)/d(zeta)
C
C 3) d(rho*H)/d(|grad(rho)|^2) = rho * d(H)/d(t^2) * d(t^2)/d(|grad(rho)|^2)
C    d(H)/d(t^2) = prefactor * Q9/Q8 below after lengthy derivation
C    d(t^2)/d(|grad(rho)|^2) = (2*k_s*phi*rho)^(-2)  [by definition of t]
C                            = 1/scale               [see above]
C

      G4 = G3*G
      T6 = T4*T_sq
      RSTHRD = RS/3.D0

c     this is d(phi)/d(zeta) in notation of paper
      GZ=(((1.d0+spin_pol)**2+eta)**sixthm-
     1((1.d0-spin_pol)**2+eta)**sixthm)/3.d0

      FAC = DELT/B+1.D0

c     this is dA/dphi in notation of paper
      BG = -3.D0*B_sq*ecunif*FAC/(BET*G4)

c     this is dA/d(ecunif) in notation of paper
      BEC = B_sq*FAC/(BET*G3)

      Q8 = Q5*Q5+DELT*Q4*Q5*T_sq
      Q9 = 1.D0+2.D0*B*T_sq

C     this is dH/dA in notation of paper
      hB = -BET*G3*B*T6*(2.D0+B*T_sq)/Q8

C     this is -1/3*r_s*dH/d(r_s) in notation of paper
      hRS = -RSTHRD*hB*BEC*ECRS

C     dh/d(zeta) in notation of paper
      hZ = 3.D0*GZ*ec_gradient_correction/G + hB*(BG*GZ+BEC*ECZET)

C     dH/d(t^2) in notation of paper
      hT2 = BET*G3*Q9/Q8

C     from here on, can deviate from orig version entirely - do not need second derivatives anywhere
C     d(rho*H) = H - 1/3*r_s*dH/d(r_s) - 7/3*t^2*dH/d(t^2)

      c_density_deriv = c_density_deriv + ec_gradient_correction
     + + hRS - 7.d0 * thrd * t_sq * hT2

C     add (1/rho)*dH/d(zeta)
C     See above; because of change of variables, apply chain rule
C     dH/d(zeta) = dh/d(zeta) + dh/d(t^2)*d(t^2)/d(phi)*d(phi)/d(zeta)
      
      ! sign was wrong!!! (R.Gehrke) 20.12.2007

!      c_spin_deriv = c_spin_deriv + hZ + 2.d0*hT2*t_sq/G*GZ
      c_spin_deriv = c_spin_deriv + hZ - 2.d0*hT2*t_sq/G*GZ

C     d(rho*H)/d(|grad(rho)|^2) = rho * d(H)/d(t^2) / scale

      c_gradient_deriv = rho * hT2 / scale

      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
c
c  VB 2005 : gcor2 is fully intact and correct; yet commented away 
c  because it is also present in pbe.f .
c
c      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
cc slimmed down version of GCOR used in PW91 routines, to interpolate
cc LSD correlation energy, as given by (10) of
cc J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
cc K. Burke, May 11, 1996.
c      IMPLICIT REAL*8 (A-H,O-Z)
c      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
c      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
c      Q2 = DLOG(1.D0+1.D0/Q1)
c      GG = Q0*Q2
c      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
c      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
c      RETURN
c      END
c--c----------------------------------------------------------------------

c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPBEsol_derivs
     +( rho, squared_grad_rho,
     +  en_density_x, x_density_deriv, x_gradient_deriv )
C
C  VB - Subroutine based on:
c
C----------------------------------------------------------------------
C  (PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)) AND
c
c  Modifications proposed in arXiv 0711.0156v2 (2008/02): As PBE, but
c  with parameters mu ("um" in code) = 0.1235 and beta ("bet") = 0.046
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C
C  INPUT: 
C
C  rho : DENSITY
C  squared_grad_rho: |grad(rho)|^2
C
C  OUTPUT:
C
C  en_density_x :     exchange energy density e(x)
C  x_density_deriv :  partial derivative d(rho*e)/d(rho)
C  x_gradient_deriv : partial derivative d(rho*e)/d[|grad(rho)|^2]
C
C  Derived quantities:
C
C  exunif: would-be exchange contribution of the uniform electron gas
C          for given density
C  scale:  scaling factor (2*KF*rho)^2 see below
C  SSQ:    squared scaled gradient [ ABS(GRAD rho)/(2*KF*rho) ]^2, 
C          where kf=(3 pi^2 rho)^(1/3)
C          dubbed s*s below
C  P0:     1 + ul * SSQ
C  FxPBE:  PBE enhancement factor (gradient correction)
C  Fs:     Partial derivative of FxPBE by squared scaled gradient
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submitted to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c       e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c   e_x[PBE]=e_x[unif]*FxPBE(s)
c   FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none

C     declare constants

      real*8 :: thrd, thrd2, thrd4, thrd8, pi, ax, um, uk, ul

      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.123456790123457d0,uk=0.8040d0,ul=um/uk)

c     declare variables

      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: en_density_x
      real*8 :: x_density_deriv
      real*8 :: x_gradient_deriv

      real*8 :: exunif
      real*8 :: scale
      real*8 :: SSQ
      real*8 :: P0
      real*8 :: FxPBE
      real*8 :: Fs

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density

      exunif = AX*rho**THRD

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      scale = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale*(rho**thrd8)
      SSQ = squared_grad_rho / scale
      P0=1.d0+ul*SSQ
      FxPBE = 1d0+uk-uk/P0
      en_density_x = exunif*FxPBE
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE partial derivatives of rho*exunif(rho)*FxPBE(ssq) :
c  first, find derivative of Fx w.r.t ssq:
c  Fs = d FxPBE/ d(ssq)
      Fs=uk*ul/(P0*P0)
C  next, multiply with derivatives d(ssq)/d(rho) and d(ssq)/d(squared_grad_rho)
C  to obtain desired partial derivatives;
C  the definition of ssq gives
C  d(ssq)/d(rho) = - 8/3 ssq/rho
C  d(ssq)/d(squared_grad_rho) = 1/ (2*KF*rho)^2 = 1/scale

      x_density_deriv = thrd4*en_density_x - exunif*Fs*thrd8*ssq

      x_gradient_deriv = rho * exunif * Fs / scale

      RETURN
      END

c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------

c######################################################################
c subroutine based on:
c
c PBEsol, see arXiv 0711.0156v2 (2008/02): As PBE, but
c with parameters mu ("um" in code) = 0.1235 and beta ("bet") = 0.046
c######################################################################

      SUBROUTINE CORPBEsol_derivs
     +( rho, squared_grad_rho, spin_pol, 
     +  en_density_c, c_density_deriv, c_gradient_deriv, c_spin_deriv
     +)
C
C  rho:               Density
C  squared_grad_rho : |grad(rho)|^2
C  spin_pol:          zeta = (rhoup-rhodn)/rho
C  en_density_c:      Correlation energy density e(x) = e_lda(x) + H(x)
C                     Really e = e( rho, squared_grad_rho, zeta )
C  c_density_deriv:   partial derivative d(rho*e)/d(rho)
C  c_gradient_deriv:  partial derivative d(rho*e)/d(|grad(rho)|^2)
C  c_spin_deriv:      partial derivative (1/rho)*d(rho*e)/d(zeta)
C
C                     Memo: c_spin_deriv formally contains an 1/rho part
C                           which is needed to form the XC potential.
C                           This means that the full XC potential could be obtained 
C                           simply by
C
C                           vcup = c_density_deriv - (zeta-1)*c_spin_deriv - 2 * div ( c_gradient_deriv * grad(rho) )
C                           vcdn = c_density_deriv - (zeta+1)*c_spin_deriv - 2 * div ( c_gradient_deriv * grad(rho) )
C
C  VB - Subroutine based on:
C
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  Original input (that's when the local potential _was_ computed)
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none

ctest
c      integer test
c      real*8 Q0,Q1,Q2
ctest end

C  constants go first

C  original constants from PBE code

c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at 
c          |zeta|=1.

      real*8 thrd, thrdm, thrd2, sixthm, thrd4

      real*8 gam, fzz, gamma, bet, delt, eta

      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.046d0,delt=bet/gamma)
      parameter(eta=1.d-12)

C     constants from Martin Fuchs' wrapper

      real*8 crs, pi, pisq3
      parameter (crs=1.91915829267751281d0)
      parameter (pi=.314159265358979312d1)
      parameter (pisq3=.296088132032680740d2)

C  variables

C     input

      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: spin_pol

C     output

      real*8 :: en_density_c 
      real*8 :: c_density_deriv
      real*8 :: c_gradient_deriv
      real*8 :: c_spin_deriv 

C     local variables

C     fk: Fermi momentum k_f
C     rs: Seitz radius
C     rtrs: sqrt(rs)
C     scale: scaling factor for t^2

      real*8 :: fk
      real*8 :: rs
      real*8 :: rtrs
      real*8 :: scale

C     for uniform electron gas term:
c     EU=unpolarized LSD correlation energy
c     EURS=dEU/drs
c     EP=fully polarized LSD correlation energy
c     EPRS=dEP/drs
c     ALFM=-spin stiffness, [c](3).
c     ALFRSM=-dalpha/drs

C     alfc: spin stiffness -alfm
C     z4 :  zeta^4 (zeta=spin polarisation spin_pol)
c     F=spin-scaling factor from [c](9).
c     ecunif: correlation energy density of the uniform electron gas.
c     ECRS = dEc/drs [c](A2)
c     ECZET = dEc/dzeta [c](A3)
c     FZ = dF/dzeta [c](A4)

      real*8 :: eu, eurs
      real*8 :: ep, eprs
      real*8 :: alfm, alfrsm

      real*8 :: alfc, z4, f

      real*8 :: ecunif

      real*8 :: ECRS 
      real*8 :: FZ 
      real*8 :: ECZET

C     for correlation correction H( rho, squared_grad_rho, zeta)
c
c     these quantities are all intermediates of the orig PBE paper, but they are 
c     strangely renamed in the original PBE code. 
c
c     G=phi(zeta), given after [a](3)
c     G3=phi(zeta)**3
C     PON = -ecunif/(phi(zeta)**3 * gam*e^2/a0)
c     B=A of [a](8)
C     B_sq = B*B
C
C     T_sq : squared scaled gradient t of PBE paper
C
C     T4 : t^4
C     Q4 : 1+At^2 - see PBE paper
C     Q5 : 1 + A t^2 + A^2 t^4 - see PBE paper
C
C     ec_gradient_correction : gradient correction to correlation energy - H of PBE paper
C
      real*8 :: G, G3, PON, B
      real*8 :: B_sq

      real*8 :: T_sq
      real*8 :: T4
      real*8 :: Q4, Q5

      real*8 :: ec_gradient_correction

      real*8 :: G4,T6,RSTHRD,GZ
      real*8 :: FAC, BG, BEC, Q8, Q9
      real*8 :: hB, hRS, hZ, hT2

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c determine derived quantities needed later

      fk=(pisq3*rho)**thrd
      rs=crs/fk

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c construct ec, using [c](8)

      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,rtRS,ALFM,ALFRSM)

      ALFC = -ALFM
      Z4 = spin_pol**4
      F=((1.D0+spin_pol)**THRD4+(1.D0-spin_pol)**THRD4-2.D0)/GAM

ctest
c      if (test.eq.1) then
c        write(6,*) "rtrs = ", rtrs
c        write(6,*) "eu = ", eu
c        Q0 = -2.D0*0.0310907D0*(1.D0+0.21370D0*rtrs*rtrs)
c        Q1 = 2.D0*0.0310907D0*rtrs*(7.5957D0+rtrs*
c     +       (3.5876D0+rtrs*(1.6382D0+0.49294D0*rtrs)))
c
c        Q2 = DLOG(1.D0+1.D0/Q1)
c
c        write(6,*) "Q0 = ", Q0
c        write(6,*) "Q1 = ", Q1
c        write(6,*) "Q2(log) = ", Q2
c
c        Q2 = 1.D0/Q1
c        write(6,*) "Q2(taylor) = ", Q2
c
cc        GG = Q0*Q2
c      end if
ctest end

      ecunif = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ

c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  Next, calculate partial derivatives of LDA energy density ...

c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)

      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ

C     add uniform part of partial derivative of correlation energy by density
      c_density_deriv = ecunif - rs*ecrs*thrd

c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)

      FZ = THRD4*((1.D0+spin_pol)**THRD-(1.D0-spin_pol)**THRD)/GAM
      ECZET = 4.D0*(spin_pol**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)

c     add uniform part of partial derivative of correlation energy by zeta
      c_spin_deriv = eczet

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)

      G=((1.d0+spin_pol)**thrd2+(1.d0-spin_pol)**thrd2)/2.d0
      G3 = G**3
      PON=-ecunif/(G3*gamma)

ctest
c      if (test) then
c        write(6,*) "ecunif = ", ecunif
c        write(6,*) "PON    = ", PON
c      end if
ctest end

cvb - this is a fix for a potential floating point exception if dexp is numerically unstable
      if (dabs(pon).gt.1.d-8) then
        ! correct expression
        B = DELT/(DEXP(PON)-1.D0)
      else
        ! use Taylor limit
        B = PON
      end if
      B_sq = B*B

      scale = 16.d0*g*g*rho*rho*fk/pi

      T_sq = squared_grad_rho / scale

ctest
c      if (test) then
c        write(6,*) "B    = ", B
c        write(6,*) "B_sq = ", B_sq
c        write(6,*) "T_sq = ", T_sq
c      end if
ctest end

      T4 = T_sq*T_sq
      Q4 = 1.D0+B*T_sq
      Q5 = 1.D0+B*T_sq+B_sq*T4

ctest
c      if (test) then
c        write(6,*) "T4 = ", T4
c        write(6,*) "Q4 = ", Q4
c        write(6,*) "Q5 = ", Q5
c      end if
ctest end

C     this is H(rs,t,zeta)
      ec_gradient_correction = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T_sq/Q5)

      en_density_c = ecunif + ec_gradient_correction

c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE partial derivatives.
C
C We want:
C 1) d(rho*H)/d(rho)
C 2) (1/rho)*d(rho*H)/d(zeta)
C 3) d(rho*H)/d(|grad(rho)|^2)
C
C Steps:
C 1) d(rho*H)/d(rho) = H + rho*dH/d(rho)
C    dH/d(rho) = dH/d(r_s)*d(r_s)/d(rho) + dH/d(t^2)*d(t^2)/d(rho)
C    d(r_s)/d(rho)=-1/3*r_s/rho
C    d(t^2)/d(rho)=-7/3*t^2/rho
C    Hence:
C    d(rho*H) = H - 1/3*r_s*dH/d(r_s) - 7/3*t^2*dH/d(t^2)
C
C 2) (1/rho)*d(rho*H)/d(zeta) = dH/d(zeta) = rho * dH/d(phi) * d(phi)/d(zeta)
C
C    Notice that zeta appears here in two different functions, technically:
C    We want d[H(rho,|grad(rho)|^2,zeta]/d(zeta) from
C    d[h(r_s,zeta,t^2)] where t^2=t^2(phi(zeta))
C    Must apply chain rule to get
C    dH/d(zeta) = dh/d(zeta) + dh/d(t^2)*d(t^2)/d(phi)*d(phi)/d(zeta)
C
C 3) d(rho*H)/d(|grad(rho)|^2) = rho * d(H)/d(t^2) * d(t^2)/d(|grad(rho)|^2)
C    d(H)/d(t^2) = prefactor * Q9/Q8 below after lengthy derivation
C    d(t^2)/d(|grad(rho)|^2) = (2*k_s*phi*rho)^(-2)  [by definition of t]
C                            = 1/scale               [see above]
C

      G4 = G3*G
      T6 = T4*T_sq
      RSTHRD = RS/3.D0

c     this is d(phi)/d(zeta) in notation of paper
      GZ=(((1.d0+spin_pol)**2+eta)**sixthm-
     1((1.d0-spin_pol)**2+eta)**sixthm)/3.d0

      FAC = DELT/B+1.D0

c     this is dA/dphi in notation of paper
      BG = -3.D0*B_sq*ecunif*FAC/(BET*G4)

c     this is dA/d(ecunif) in notation of paper
      BEC = B_sq*FAC/(BET*G3)

      Q8 = Q5*Q5+DELT*Q4*Q5*T_sq
      Q9 = 1.D0+2.D0*B*T_sq

C     this is dH/dA in notation of paper
      hB = -BET*G3*B*T6*(2.D0+B*T_sq)/Q8

C     this is -1/3*r_s*dH/d(r_s) in notation of paper
      hRS = -RSTHRD*hB*BEC*ECRS

C     dh/d(zeta) in notation of paper
      hZ = 3.D0*GZ*ec_gradient_correction/G + hB*(BG*GZ+BEC*ECZET)

C     dH/d(t^2) in notation of paper
      hT2 = BET*G3*Q9/Q8

C     from here on, can deviate from orig version entirely - do not need second derivatives anywhere
C     d(rho*H) = H - 1/3*r_s*dH/d(r_s) - 7/3*t^2*dH/d(t^2)

      c_density_deriv = c_density_deriv + ec_gradient_correction
     + + hRS - 7.d0 * thrd * t_sq * hT2

C     add (1/rho)*dH/d(zeta)
C     See above; because of change of variables, apply chain rule
C     dH/d(zeta) = dh/d(zeta) + dh/d(t^2)*d(t^2)/d(phi)*d(phi)/d(zeta)
      
      ! sign was wrong!!! (R.Gehrke) 20.12.2007

!      c_spin_deriv = c_spin_deriv + hZ + 2.d0*hT2*t_sq/G*GZ
      c_spin_deriv = c_spin_deriv + hZ - 2.d0*hT2*t_sq/G*GZ

C     d(rho*H)/d(|grad(rho)|^2) = rho * d(H)/d(t^2) / scale

      c_gradient_deriv = rho * hT2 / scale

      RETURN
      END
