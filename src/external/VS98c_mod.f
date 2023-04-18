      SUBROUTINE  VS98c(Tol,F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy)
c     ******************************************************************
c     *                                                                *
c     *  VS98c evaluates the correlation part of VS98 and its          *
c     *  variants on a grid.                                           *
c     *  !!! Second derivatives are not available yet.                 *
c     *                                                                *
c     *  Ref:  T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109,   *
c     *        400 (1998).                                             *
c     *       Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125,           * 
c     *        194101 (2006).                                          *
c     *                                                                *
c     *       ijzy - 1 correlation functional in VS98                  *
c     *       ijzy - 2 correlation functional in M06-L                 *
c     *       ijzy - 3 correlation functional in M06-HF                *
c     *       ijzy - 4 correlation functional in M06                   *
c     *       ijzy - 5 correlation functional in M06-2X                *
c     *                                                                *
c     *  OUTPUT:                                                       *
c     *     F      - Functional values                                 *
c     *     D1F    - First derivatives with respect to RA, RB, GA, GB  *
c     *              TA, TB                                            *
c     *                                                                *
c     *  INPUT:                                                        *
c     *     RA,B   - Spin densities                                    *
c     *     D1RA,B - Spin density gradients                            *
c     *     TA,B   - Spin kinetic energy densities                     *
c     *     NGrid  - number of grids                                   *
c*                                                                      *
c*  Note:  You need to code up the Perdew-Wang uniform-gas correlation  *
c*          lsdac                                                       *
c*  YZ (12/08)                                                          *
c*                                                                      *
c************************************************************************
      IMPLICIT NONE

      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC
      REAL*8  F,D1F(7),RA,RB,
     $        D1RA(3),D1RB(3),TA,TB
      REAL*8 Pi, F6, F43, Pi34, F13, 
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ,
     &  P, EUEG,Denom, DenPA, DenPB, DenGA, DenGB,EUEGPA,EUEGPB
      REAL*8 PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA,ChiAP,
     &ChiAG,ZA,ZAP,ZAT
      REAL*8 PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB,ChiBP,
     &ChiBG,ZB,ZBP,ZBT
      REAL*8 ZAB, XAB, kab, xk, zk
      REAL*8 dgdx,dgdz,dgdPA,dgdGA,dgdTA,dgdPB,dgdGB,dgdTB
      REAL*8 gcab,cf,gab
      REAL*8 r7, r8, r9, r10, r11, r12

      
      REAL*8  Tol, DTol,F1, F2, F3, F4, COpp
      Data F1/1.0d0/,F2/2.0d0/,F3/3.0d0/,F4/4.0d0/ 
     $,gab/0.00304966d0/,cf/9.115599720d0/

      INTEGER i,j,ijzy

C     Global assignments
C      dRA = 1
C      dRB = 2
C      dTA = 3
C      dTB = 4
C      dGA = 5
C      dGB = 6
C     Global assignments
      dRA = 1
      dRB = 2
      dGA = 3
      dGB = 4
      dGC = 5
      dTA = 6
      dTB = 7
C Initializing all stupid variables
       PA=0.0
       GAA=0.0
       TauA=0.0
       FA=0.0
       FPA=0.0
       FGA=0.0
       FTA=0.0
       EUA=0.0
       EUEGA=0.0
       ChiA=0.0
       EUPA=0.0
       ChiAP=0.0
       ChiAG=0.0
       PB=0.0
       GBB=0.0
       TauB=0.0
       FB=0.0
       FPB=0.0 
       FGB=0.0
       FTB=0.0
       EUB=0.0
       EUEGB=0.0
       ChiB=0.0
       EUPB=0.0
       ChiBP=0.0
       ChiBG=0.0 
       ZA = 0.0
       ZAP = 0.0
       ZAT = 0.0
       ZB = 0.0
       ZBP = 0.0
       ZBT = 0.0
C     Parameters for VS98 
      if (ijzy.eq.1) then
              r7=   7.035010d-01
              r8=   7.694574d-03
              r9=   5.152765d-02
              r10=   3.394308d-05
              r11=  -1.269420d-03
              r12=   1.296118d-03
C     Parameters for M06-L
      elseif (ijzy.eq.2) then
              r7=      3.957626D-01
              r8=      -5.614546D-01
              r9=      1.403963D-02
              r10=     9.831442D-04
              r11=     -3.577176D-03
              r12=     0.000000D+00
C     Parameters for M06-HF
      elseif (ijzy.eq.3) then
              r7=    -6.746338D-01
              r8=    -1.534002D-01
              r9=    -9.021521D-02
              r10=   -1.292037D-03
              r11=   -2.352983D-04
              r12=   0.000000D+00

C     Parameters for M06
      elseif (ijzy.eq.4) then
               r7= -2.741539D+00
               r8= -6.720113D-01
               r9= -7.932688D-02
               r10=1.918681D-03
               r11=-2.032902D-03
               r12=0.000000D+00

C     Parameters for M06-2X
      elseif (ijzy.eq.5) then
              r7=  1.166404D-01
              r8=  -9.120847D-02
              r9=  -6.726189D-02
              r10= 6.720580D-05
              r11= 8.448011D-04
              r12= 0.000000D+00
      endif

      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3
      DTol = Tol 
       
c      DO i = 1,NGrid
         PA = RA
         IF (PA.gt.DTol.and.TA.gt.DTol) THEN
            GAA   =  D1RA(1)**2 + D1RA(2)**2 + D1RA(3)**2
C            IF (SQRT(GAA).gt.DTol) THEN ! AJL 
               TauA = TA

               Call vs98ss(PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,ZA,
     &                   ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,ijzy)
               F = F + FA
               D1F(dRA)  = D1F(dRA) + FPA
               D1F(dGA) = D1F(dGA) + FGA
               D1F(dTA)  = D1F(dTA) + FTA
C            ENDIF
         ENDIF
         PB = RB
         IF (PB.gt.DTol.and.TB.gt.DTol) THEN
            GBB   =  D1RB(1)**2 + D1RB(2)**2 + D1RB(3)**2
C            IF (SQRT(GBB).gt.DTol) THEN ! AJL
               TauB = TB

               Call vs98ss(PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,ZB,
     &                   ChiB,EUPB,ChiBP,ChiBG,ZBP,ZBT,ijzy)
               F = F + FB
               D1F(dRB)  = D1F(dRB) + FPB
               D1F(dGB) = D1F(dGB) + FGB
               D1F(dTB)  = D1F(dTB) + FTB
C            ENDIF
         ENDIF
         P=PA+PB
         IF (PB.gt.DTol.and.PA.gt.DTol) THEN
          RS = (Pi34/P) ** F13 
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ)
          EUEG = P*PotLC - EUA - EUB
          ZAB = ZA + ZB
          XAB = ChiA+ChiB
          kab = F1 + gab*(XAB+ZAB)
          xk = XAB/kab
          zk = ZAB/kab
       call gvt4(gcab,dgdx,dgdz,xk,zk,kab,gab,r7,r8,r9,r10,r11,r12)
          F = F + gcab*EUEG
          dgdPA = dgdx*ChiAP + dgdz*ZAP
          dgdGA = dgdx*ChiAG
          dgdTA = dgdz*ZAT
          dgdPB = dgdx*ChiBP + dgdz*ZBP
          dgdGB = dgdx*ChiBG
          dgdTB = dgdz*ZBT
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          D1F(dRA) = D1F(dRA) + (EUEGPA*gcab + EUEG*dgdPA)
          D1F(dRB) = D1F(dRB) + (EUEGPB*gcab + EUEG*dgdPB)
          D1F(dGA) = D1F(dGA) + EUEG*dgdGA 
          D1F(dGB) = D1F(dGB) + EUEG*dgdGB
          D1F(dTA)=  D1F(dTA) + EUEG*dgdTA
          D1F(dTB)=  D1F(dTB) + EUEG*dgdTB 
         ENDIF
c        ENDDO
        
        RETURN
        END  
                                        

      Subroutine vs98ss(PX,GX,TX,F,FP,FG,FT,EUEG,Z,Chi,EUEGP,
     &                   ChiP,ChiG,ZP,ZT,ijzy)
      Implicit none
C
C     Compute the same-spin part of the vs98 correlation functional for one grid
C     point and one spin-case.
C

      integer ijzy
      double precision r13, r14, r15, r16, r17, r18
      double precision PX, GX, TX, F, FP, FG, FT, DTol, Z, ZP, ZT
      double precision EUEG, Chi, EUEGP, ChiP, ChiG, cf, gcc
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision RS, D, RSP, PotLC, DX, DZ, dgdP, dgdG, dgdT
      double precision E,DP, DG, DT, rhoo, rho43, rho53, rho83
      double precision rrho, F4o3, rho13, kc, xk, zk, gc, dgdx, dgdz
      double precision d2LdSS, d2LdSZ, d2LdZZ, dLdS, dLdZ

      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  gcc/0.00515088d0/,cf/9.115599720d0/
 
 
      F4o3 = 4.0d0/3.0d0
C     Parameters for VS98 
      if (ijzy.eq.1) then
              r13=   3.270912d-01
              r14=  -3.228915d-02
              r15=  -2.942406d-02
              r16=   2.134222d-03
              r17=  -5.451559d-03
              r18=   1.577575d-02
C     Parameters for M06-L
      elseif (ijzy.eq.2) then
              r13=   4.650534D-01
              r14=   1.617589D-01
              r15=   1.833657D-01
              r16=   4.692100D-04
              r17=  -4.990573D-03
              r18=   0.000000D+00
C     Parameters for M06-HF
      elseif (ijzy.eq.3) then
              r13=   8.976746D-01
              r14=  -2.345830D-01
              r15=   2.368173D-01
              r16=  -9.913890D-04
              r17=  -1.146165D-02
              r18=   0.000000D+00
C     Parameters for M06
      elseif (ijzy.eq.4) then
               r13=  4.905945D-01
               r14= -1.437348D-01
               r15=  2.357824D-01
               r16=  1.871015D-03
               r17= -3.788963D-03
               r18=  0.000000D+00
C     Parameters for M06-2X
      elseif (ijzy.eq.5) then
              r13=  6.902145D-01
              r14=  9.847204D-02
              r15=  2.214797D-01
              r16= -1.968264D-03
              r17= -6.775479D-03
              r18=  0.000000D+00
      endif

      DTol =1.0d-10
      If(PX.le.DTol) then
        EUEG = Zero
        Chi = Zero
        EUEGP = Zero
        ChiP = Zero
        ChiG = Zero
        PX = Zero
        GX = Zero 
        TX = Zero
        F  = Zero
        FP = Zero
        FG = Zero
        FT = Zero
        Z  = Zero
        ZP = Zero
        ZT = Zero
      else
        Pi = F4*ATan(F1)
        Pi34 = F3 / (F4*Pi)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        rhoo = PX
c       I did this, assuming it is the same as in M06x... MR
        rrho = 1.0d0/rhoo 
        rho43 = rhoo**F4o3
        rho13 = rho43*rrho
        rho53 = rhoo**F53
        rho83 = rho53*rhoo
        
        RS = (Pi34/PX) ** F13
        Call lsdac(RS,F1,PotLC,dLdS,dLdZ)
        EUEG = PX*PotLC
        Chi = GX/rho83
        Z = (TX/rho53) - cf
        kc = F1 + gcc*(Chi + Z)
        xk = Chi/kc
        zk = Z/kc
        D = F1 - Chi/(F4*(Z + cf)) 
        call gvt4(gc,dgdx,dgdz,xk,zk,kc,gcc,r13,r14,r15,r16,r17,r18)
        E = D*EUEG*gc
        F = E 
c
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        ZP = -F53 * TX/rho83
        ZT =  F1/rho53
        DZ = Chi/(F4*(Z + cf)*(Z + cf)) 
        DX = -F1/(F4*(Z + cf))
        DP = DZ*ZP + DX*ChiP
        DG = DX*ChiG
        DT = DZ*ZT
        dgdP = dgdx*ChiP + dgdz*ZP
        dgdG = dgdx*ChiG 
        dgdT = dgdz*ZT
        EUEGP = PotLC + PX*dLdS*RSP
        FP = DP*EUEG*gc + D*EUEGP*gc + D*EUEG*dgdP
        FG = DG*EUEG*gc + D*EUEG*dgdG
        FT = DT*EUEG*gc + D*EUEG*dgdT
       Endif
       Return
       End


c lsdac starts here
c This seems to be a code duplication? AJL
c I'm going to turn it in to a wrapper around the lsda correlation
c which exists in gga91_sr.f as it is the same code and removes
c one source of possible bugs in my implementation of M06 potentials,
c though I think the problem is more likely small numbers....

      Subroutine lsdac(rs,zet,ec,ecrs,eczet)

      ! implicit real*8 (A-H,O-Z)

      ! Variables are as recieved, except vcup and vcdown, which are the discarded potential,
      ! and alfc which is the correlation stiffness

      real*8 rs, zet, ec, ecrs, eczet

      ! To be discarded

      real*8 vcup, vcdn, alfc

      call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)

      return
      end

Cc  uniform-gas correlation of perdew and wang 1991
Cc  input: seitz radius (rs), relative spin polarization (zet)
Cc  output: correlation energy per electron (ec), up- and down-spin
Cc     potentials (vcup,vcdn), derivatives of ec wrt rs (ecrs) & zet (eczet)
Cc  output: correlation contribution (alfc) to the spin stiffness
C      implicit real*8 (a-h,o-z)
C      data gam,fzz/0.5198421d0,1.709921d0/
C      data thrd,thrd4/0.333333333333d0,1.333333333333d0/
C      f = ((1.d0+zet)**thrd4+(1.d0-zet)**thrd4-2.d0)/gam
C      call gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,
C     1    0.49294d0,1.00d0,rs,eu,eurs)
C      call gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,
C     1    0.62517d0,1.00d0,rs,ep,eprs)
C      call gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,
C     1    0.49671d0,1.00d0,rs,alfm,alfrsm)
Cc  alfm is minus the spin stiffness alfc
C      alfc = -alfm
C      z4 = zet**4
C      ec = eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
Cc  energy done. now the potential:
C      ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
C      fz = thrd4*((1.d0+zet)**thrd-(1.d0-zet)**thrd)/gam
C      eczet = 4.d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu
C     1        -(1.d0-z4)*alfm/fzz)
C      comm = ec -rs*ecrs/3.d0-zet*eczet
Cc      vcup = comm + eczet
Cc      vcdn = comm - eczet
C
C      Return
C      End

