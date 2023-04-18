! This is an edited version of CORPBE_derivs, taken from
! the file pbe_partial_derivs.f. In principle, we could re-write
! the subroutines in pbe_partial_derivs.f to make this data accessible
! mid-calculation, and I welcome anyone who is willing to do this,
! as all I have done is make the partial derivatives accessible from
! midway through the calculation of the potential, namely the derivative
! wrt to rho, grad_rho and zeta
!
! However, as I've spent far too long looking at this numbers by now,
! I'm going to halt any further work I do in this area and 
! concentrate now on getting the range-separated meta-gga functionals
! working as that was my original purpose for pursuing this.
!
! Everything I have changed from the original subroutine is commented with 
! double !! marks, so that transferability is possible. AJL

C      SUBROUTINE CORPBE_derivs
C     +( rho, squared_grad_rho, spin_pol, 
C     +  en_density_c, c_density_deriv, c_gradient_deriv, c_spin_deriv
C     +)


       SUBROUTINE PBEH0(rho,grad_rho,spin_pol,
     +                  ecunif,ecrs,eczet,ec_gradient_correction,
     +                  c_density_deriv,c_gradient_deriv,
     +                  c_spin_deriv)

C
C  rho:               Density
!!  grad_rho:         |grad(rho)|, used to calculate squared_grad_rho VVV
C  spin_pol:          zeta = (rhoup-rhodn)/rho
!! ecunif:            e_lda(x)
!! ecrs:              partial derivative d(ec)/d(rs)
!! eczet:             partial derivative d(ec)/d(zet)
!! ec_gradient_correction: H(x)
!! c_density_deriv:   partial derivative d(H)/d(rho) 
!! c_gradient_deriv:  partial derivative d(H)/d(|grad(rho)|)
!! c_spin_deriv:      partial derivative d(H)/d(zeta)
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
!! AJL: Additional variable
      real*8 :: grad_rho
!! AJL
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

!! AJL: Additional definition
      squared_grad_rho = grad_rho**2
!! AJL

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
!! AJL: Not needed
!!      c_density_deriv = ecunif - rs*ecrs*thrd
!! AJL

c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)

      FZ = THRD4*((1.D0+spin_pol)**THRD-(1.D0-spin_pol)**THRD)/GAM
      ECZET = 4.D0*(spin_pol**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)

c     add uniform part of partial derivative of correlation energy by zeta
!! AJL: Not needed
!!      c_spin_deriv = eczet
!! AJL

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
!      if ((dabs(pon).gt.1.d-8).and.(dabs(pon).lt.40.0d0)) then
      if (dabs(pon).gt.1.d-8) then
        ! correct expression
        B = DELT/(DEXP(PON)-1.D0)
!      else if (dabs(pon).ge.40.0d0) then
!        B = DELT/(-1.D0)
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

!! AJL: Not needed
!!      en_density_c = ecunif + ec_gradient_correction
!! AJL

c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE partial derivatives.
C
C We want:
C 1) d(rho*H)/d(rho)
C 2) (1/rho)*d(rho*H)/d(zeta)
C 3) d(rho*H)/d(|grad(rho)|^2)
!!
!! AJL: At the end of the subroutine we manipulate these to be:
!! 1) d(H)/d(rho)
!! 2) d(H)/d(zeta) [implicitly, this is what we have anyway
!! 3) d(H)/d(|grad(rho)|)
!! AJL
!!
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

!! AJL
!!      c_density_deriv = c_density_deriv + ec_gradient_correction
!!     + + hRS - 7.d0 * thrd * t_sq * hT2
       c_density_deriv = hRS - 7.d0 * thrd * t_sq * hT2

C     add (1/rho)*dH/d(zeta)
C     See above; because of change of variables, apply chain rule
C     dH/d(zeta) = dh/d(zeta) + dh/d(t^2)*d(t^2)/d(phi)*d(phi)/d(zeta)

      ! sign was wrong!!! (R.Gehrke) 20.12.2007

!      c_spin_deriv = c_spin_deriv + hZ + 2.d0*hT2*t_sq/G*GZ
!!      c_spin_deriv = c_spin_deriv + hZ - 2.d0*hT2*t_sq/G*GZ
!! AJL
      c_spin_deriv = hZ - 2.d0*hT2*t_sq/G*GZ
!! AJL

C     d(rho*H)/d(|grad(rho)|^2) = rho * d(H)/d(t^2) / scale

      c_gradient_deriv = rho * hT2 / scale

!! AJL
      c_density_deriv  = 1.0d0 * c_density_deriv / rho
      c_gradient_deriv = 2.0d0 * grad_rho * c_gradient_deriv / rho
      c_spin_deriv     = 1.0d0 * c_spin_deriv 
!! AJL

      RETURN
      END
