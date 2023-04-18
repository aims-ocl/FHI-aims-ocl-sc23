      SUBROUTINE EXCHREVPBE_derivs
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
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
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
      parameter(um=0.2195149727645171d0,uk=1.2450d0,ul=um/uk)

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
