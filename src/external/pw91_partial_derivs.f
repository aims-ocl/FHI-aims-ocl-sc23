C---------------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPW91_derivs
     +( rho, squared_grad_rho,
     +  en_density_x, x_density_deriv, x_gradient_deriv )
C
C  VB - Subroutine based on:
c
C----------------------------------------------------------------------
C  PW91_gga EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
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
C
C  S  :    gradient [ ABS(GRAD rho)/(2*KF*rho) ], 
C          where kf=(3 pi^2 rho)^(1/3) 
C  S2:    squared scaled gradient [ ABS(GRAD rho)/(2*KF*rho) ]^2, 
C          where kf=(3 pi^2 rho)^(1/3)
C          dubbed s*s below
C  FxPW91:  PW91 enhancement factor (gradient correction)
C  Fs:     Partial derivative of FxPW91 by squared scaled gradient
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c Perdew and Wang Phys Rev. B 46 6671 (1992)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c      e_x[unif]=ax*rho^(4/3)  [LDA]
c      ax = -0.75*(3/pi)^(1/3)
c      e_x[pw91]=e_x[unif]*FxPW91(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT none

C     declare constants

      real*8 :: thrd, thrd2, thrd4, thrd8, pi, ax, um, uk, ul
      real*8 :: a1, a2 ,a3, a4, a, b1
      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
C     PW91
      parameter(a1=0.19645d0,a2=0.27430d0,a3=0.15084d0,
     +           a4=100.d0)
      parameter(a=7.7956d0,b1=0.004d0)

c     declare variables

      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: en_density_x
      real*8 :: x_density_deriv
      real*8 :: x_gradient_deriv

      real*8 :: exunif
      real*8 :: scale
      real*8 :: s
      real*8 :: FxPW91
      real*8 :: Fs,dFx
      real*8 :: s2,s3,s4
      real*8 :: p0,p1,p2,p3,p4
      real*8 :: p5,p6,p7
      
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density

      exunif = AX*rho**THRD

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PW91 enhancement factor
      scale = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale*(rho**thrd8)
      s2 = squared_grad_rho / scale
      s = dsqrt(s2)
      s3 = s2*s
      s4 = s2*s2
      p0 = 1.d0/dsqrt(1.d0+a*a*s2)
      p1 = dlog(a*s+1.d0/p0)
      p2 = dexp(-a4*s2)      
      
      p3 = 1.d0/(1.d0+a1*s*p1+b1*s4)
      p4 = 1.d0+a1*s*p1+(a2-a3*p2)*s2
      FxPW91 = p3*p4

      en_density_x = exunif*FxPW91
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE partial derivatives of rho*exunif(rho)*FxPW91(ssq) :
c  first, find derivative of Fx w.r.t s2:
c  dFx = d FxPW91/ d(s)

      p5 = b1*s2-(a2-a3*p2)
      p6 = a1*s*(p1+a*s*p0)
      p7 = 2.D0*(a2-a3*p2)+2.D0*a3*a4*s2*p2-4.D0*b1*s2*FxPW91
      dFx = p3*(p3*p5*p6+p7)*s
C  next, multiply with derivatives d(ssq)/d(rho) and d(ssq)/d(squared_grad_rho)
C  to obtain desired partial derivatives;
C  the definition of ssq gives
C  d(ssq)/d(rho) = - 8/3 ssq/rho
C  d(ssq)/d(squared_grad_rho) = 1/ (2*KF*rho)^2 = 1/scale
C
C  d(s)/d(s2)
      Fs = dFx/2.0d0/s

      x_density_deriv = thrd4*en_density_x - exunif*Fs*thrd8*s2

      x_gradient_deriv = rho * exunif * Fs / scale

      RETURN
      END      

c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORPW91_derivs
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
C  MR - Subroutine based on:
C
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c--C  pw91 CORRELATION, modified by K. Burke to put all arguments 
c--c  as variables in calling statement, rather than in common block
c--c  May, 1996.
c--C  INPUT RS: SEITZ RADIUS
c--C  INPUT ZET: RELATIVE SPIN POLARIZATION
c--C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
c--C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
c--C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
c--C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
c--C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c--C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
c
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
      real*8 xnu, cc0,cx,alf, c1,c2, c3,c4,c5,c6,a4
      real*8 gam, fzz, eta

      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))

      parameter(eta=1.d-12)

      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
      parameter(alf=0.09D0)
      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)

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
C      real*8 :: G, G3, PON, B
C      real*8 :: B_sq

C      real*8 :: T_sq
C      real*8 :: T4
C      real*8 :: Q4, Q5

      real*8 :: ec_gradient_correction

C      real*8 :: G4,T6,RSTHRD,GZ
C      real*8 :: FAC, BG, BEC, Q8, Q9
C      real*8 :: hB, hRS, hZ, hT2


c     PW91 code
      real*8 bet,delt
      real*8 G, G3, G4, PON, B, B2
      real*8 t2, t4, t6
      real*8 q4, q5, q6, q7, q8, q9
      real*8 cc, r0, r1, r2, r3,r4
      real*8 rs2, rs3
      real*8 coef
      real*8 ccrs,rsthrd
      real*8 gz, fac, bg, bec
      real*8 h0, h1
      real*8 h0b, h0rs, h1rs, h0z, h1z, h0t, h1t
      real*8 hrs, ht, hz
      
      
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
c PW91 correlation energy

      bet=xnu*cc0
      delt = 2.D0*alf/bet

      g=((1.d0+spin_pol)**thrd2+(1.d0-spin_pol)**thrd2)/2.d0
      g3 = g**3
      g4=g3*g
      PON = -delt*ecunif/(g3*bet)

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
      b2 = b*b
c      B_sq = B*B

      scale = 16.d0*g*g*rho*rho*fk/pi

C      T_sq = squared_grad_rho / scale
      t2 = squared_grad_rho / scale
ctest
c      if (test) then
c        write(6,*) "B    = ", B
c        write(6,*) "B_sq = ", B_sq
c        write(6,*) "T_sq = ", T_sq
c      end if
ctest end

      t4 = t2*t2
      t6 = t4*t2
      rs2 = rs*rs
      rs3 = rs2*rs
      q4 = 1.D0+b*t2
      q5 = 1.D0+b*t2+b2*t4
      q6 = c1+c2*rs+c3*rs2
      q7 = 1.D0+c4*rs+c5*rs2+c6*rs3
      cc = -cx + q6/q7
      r0 = 0.663436444d0*rs
      r1 = a4*r0*g4

      coef = cc-cc0-3.D0*cx/7.D0
      r2 = xnu*coef*g3
      r3 = DEXP(-r1*t2)
      h0 = g3*(bet/delt)*DLOG(1.D0+delt*q4*t2/q5)
      h1 = r3*r2*t2
      
ctest
c      if (test) then
c        write(6,*) "T4 = ", T4
c        write(6,*) "Q4 = ", Q4
c        write(6,*) "Q5 = ", Q5
c      end if
ctest end

C     this is H(rs,t,zeta)
      ec_gradient_correction = h0+h1
C      ec_gradient_correction = h0
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

      ccrs = (c2+2.*c3*rs)/q7 - q6*(c4+2.*c5*rs+3.*c6*rs2)/q7**2
      rsthrd = rs/3.D0
      r4 = rsthrd*ccrs/coef

      GZ=(((1.d0+spin_pol)**2+eta)**sixthm-
     1((1.d0-spin_pol)**2+eta)**sixthm)/3.d0
      fac = delt/b+1.D0

c     this is dA/dphi in notation of paper (PBE paper)
      bg = -3.D0*b2*ecunif*fac/(bet*g4)

c     this is dA/d(ecunif) in notation of paper
      bec = b2*fac/(bet*g3)

      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.D0+2.D0*b*t2

C     this is dH/dA in notation of paper (PBE comment)
      h0b = -bet*g3*b*t6*(2.D0+b*t2)/q8
      
C     this is -1/3*r_s*dH/d(r_s) in notation of paper (PBE comment)
      h0rs = -rsthrd*h0b*bec*ecrs
      h1rs = r3*r2*t2*(-r4+r1*t2/3.D0)

C     dh/d(zeta) in notation of paper (PBE comment)
      h0z = 3.D0*gz*h0/g + h0b*(bg*gz+bec*eczet)
      h1z = gz*r3*r2*t2*(3.D0-4.D0*r1*t2)/g

C     dH/d(t^2) in notation of paper
       h0t = bet*g3*q9/q8
       h1t = r3*r2*(1.D0-r1*t2)


       hrs=h0rs+h1rs
       ht=h0t+h1t
       hz=h0z+h1z

C     from here on, can deviate from orig version entirely - do not need second derivatives anywhere
C     d(rho*H) = H - 1/3*r_s*dH/d(r_s) - 7/3*t^2*dH/d(t^2)

      c_density_deriv = c_density_deriv + ec_gradient_correction
     + + hrs - 7.d0 * thrd * t2 * ht

C     add (1/rho)*dH/d(zeta)
C     See above; because of change of variables, apply chain rule
C     dH/d(zeta) = dh/d(zeta) + dh/d(t^2)*d(t^2)/d(phi)*d(phi)/d(zeta)
      
      ! sign was wrong!!! (R.Gehrke) 20.12.2007

!      c_spin_deriv = c_spin_deriv + hZ + 2.d0*hT2*t_sq/G*GZ
      c_spin_deriv = c_spin_deriv + hz - 2.d0*ht*t2/g*gz

C     d(rho*H)/d(|grad(rho)|^2) = rho * d(H)/d(t^2) / scale

      c_gradient_deriv = rho * ht / scale

      RETURN
      END
