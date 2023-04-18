      SUBROUTINE EXCHPBEint_derivs
     +( rho, squared_grad_rho,
     +  en_density_x, x_density_deriv, x_gradient_deriv )
C
C  PBEint exchange - energy density and derivative
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
c [c]  E. Fabiano, L. A. Constantin, F. Della Sala, 
c       Phys. Rev. B 82, 113104 (2010). 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none

C     declare constants

      real*8 :: thrd, thrd2, thrd4, thrd8, pi, ax, um, uk,umsol
      real*8 :: thrd5,umint
      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0,thrd5=5.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)

c     declare variables
      common /maxss_comm/ amaxss,domaxss
      double precision amaxss
      logical domaxss

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

      real*8 :: upbe,uge,fx,dfx,t,tmp1,fden,k,aa,ab,s


      upbe = 0.2195d0
      uge =  10.d0/81.d0
      k = 0.804d0

      aa = 0.197d0
      ab= aa
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density

      exunif = AX*rho**THRD

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBEint enhancement factor
      scale = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale*(rho**thrd8)
      SSQ = squared_grad_rho / scale


      if (domaxss) then
      if (dsqrt(ssq).gt.amaxss) then
       amaxss=dsqrt(ssq)
       write(888,*) amaxss
      endif
      endif

       t=SSQ 
       s=dsqrt(ssq)
   
     
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE partial derivatives of rho*exunif(rho)*FxPBE(ssq) :
c  first, find derivative of Fx w.r.t ssq:
c  Fs = d FxPBE/ d(ssq)

      fx = 1 +k-k / (1 + (uge + (upbe - uge) * aa * s ** 2 / (1 + ab
     +* s ** 2)) * s ** 2 / k)


      FxPBE=fx
      en_density_x = exunif*FxPBE

      dfx = k/(1 +(uge + (upbe - uge) * aa * s ** 2 / (1 + ab * s **
     +2)) * s **2/k)**2 * ((2 * (upbe - uge) * aa * s / (1 + ab * s
     + ** 2) - 2*(upbe-uge) * aa * s ** 3 / (1 + ab * s ** 2) ** 2 *
     + ab) * s**2/k+2 * (uge + (upbe - uge) * aa * s ** 2 / (1 + a
     +b * s ** 2)) * s / k)

      fs=dfx/2.d0/s
      if(s.le.0.000001)dfx=uge
  
      x_density_deriv = thrd4*en_density_x - exunif*Fs*thrd8*ssq
      x_gradient_deriv = rho * exunif * Fs / scale

      RETURN
      END

C=====================================================================
C=====================================================================


      SUBROUTINE CORPBEint_derivs
     +( rho, squared_grad_rho, spin_pol, 
     +  en_density_c, c_density_deriv, c_gradient_deriv, c_spin_deriv
     +)
C
C PBE int correlation functional
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
c derived from the
c  Official PBE correlation code by K. Burke, May 14, 1996.
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

      real*8 thrd, thrdm, thrd2, sixthm, thrd4

      real*8 gam, fzz, gamma, bet, delt, eta

      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.052d0,delt=bet/gamma)
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

      ecunif = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ

c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  Next, calculate partial derivatives of LDA energy density ...

c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)

      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ

C     add uniform part of partial derivative of correlation energy by density
      c_density_deriv = ecunif - rs*ecrs*thrd

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


      T4 = T_sq*T_sq
      Q4 = 1.D0+B*T_sq
      Q5 = 1.D0+B*T_sq+B_sq*T4


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
c----------
