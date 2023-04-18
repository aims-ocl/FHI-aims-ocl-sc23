      Subroutine M06c(F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy)

c************************************************************************
c*                                                                      *
c*  CM06 evaluates the M05 correlation part of the M06 suite of         *
c*  functionals on the grids.                                           *
c*  !!! Second derivatives are not available yet.                       *
c*                                                                      *
c*  Ref: (a) Zhao, Y.  and Truhlar, D. G. J. Chem. Phys. 125,           *
c*    194101 (2006).                                                    *
c*       (b) Y. Zhao and D. G. Truhlar, J. Phys. Chem. A (2006),        *
c*    110(49),13126-13130.                                              *
c*                                                                      *
c*       ijzy - 1 M06-L                                                 *
c*       ijzy - 2 M06-HF                                                *
c*       ijzy - 3 M06                                                   *
c*       ijzy - 4 M06-2X                                                *
c*                                                                      *
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

c      INTEGER NGrid
      REAL*8  F,D1F(7),RA,RB,
     $        D1RA(3),D1RB(3),TA,TB
      REAL*8 Pi, F6, F43, Pi34, F13, 
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ,
     &  P, EUEG,Denom, DenPA, DenPB, DenGA, DenGB
      REAL*8 PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA,ChiAP,ChiAG
      REAL*8 PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB,ChiBP,ChiBG
      REAL*8 sopp0, sopp1,sopp2, sopp3, sopp4
      REAL*8 U, W, dUdChiA,dUdChiB,dUdPA,dUdPB,dUdGA,dUdGB,
     &dWdU,dWdPA,dWdPB, dWdGA,dWdGB,EUEGPA,EUEGPB
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC

      
      REAL*8  Tol, DTol,F1, F2, F3, F4, COpp
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,F3/3.0d0/,F4/4.0d0/ 


      INTEGER i,j,ijzy
C     Global assignments
      dRA = 1
      dRB = 2
      dGA = 3
      dGB = 4
      dGC = 5
      dTA = 6
      dTB = 7
      Tol = 1.d-10
C The variables initialization...
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
      if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
         sopp0= 6.042374D-01
         sopp1= 1.776783D+02
         sopp2= -2.513252D+02
         sopp3= 7.635173D+01
         sopp4= -1.255699D+01
      elseif (ijzy.eq.2) then
C     Parameters for M06-HF Correlation
         sopp0= 1.674634D+00
         sopp1= 5.732017D+01
         sopp2= 5.955416D+01
         sopp3= -2.311007D+02
         sopp4= 1.255199D+02
      elseif (ijzy.eq.3) then
C     Parameters for M06 Correlation
         sopp0= 3.741539D+00
         sopp1= 2.187098D+02
         sopp2= -4.531252D+02
         sopp3= 2.936479D+02
         sopp4= -6.287470D+01
      elseif (ijzy.eq.4) then
C     Parameters for M06-2X Correlation
         sopp0= 8.833596D-01
         sopp1= 3.357972D+01
         sopp2= -7.043548D+01
         sopp3= 4.978271D+01
         sopp4= -1.852891D+01
      endif
      DTol = Tol 
      CALL VS98c(Tol,F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy+1)
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

c      DO i = 1,NGrid
         PA = RA
         IF (PA.gt.DTol.and.TA.gt.DTol) THEN
            GAA   =  D1RA(1)**2 + D1RA(2)**2 + D1RA(3)**2
C            IF (SQRT(GAA).gt.DTol) THEN ! AJL
               TauA = TA

               call m06css(DTol,PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,
     &                   ChiA,EUPA,ChiAP,ChiAG,ijzy)
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

               call m06css(DTol,PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,
     &                   ChiB,EUPB,ChiBP,ChiBG,ijzy)
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
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
          F = F + EUEG*W
          dUdChiA =COpp/(F1 + COpp*(ChiA+ChiB))**2
          dUdChiB =COpp/(F1 + COpp*(ChiA+ChiB))**2
          dUdPA= dUdChiA*ChiAP
          dUdPB= dUdChiB*ChiBP
          dUdGA= dUdChiA*ChiAG
          dUdGB= dUdChiB*ChiBG
          dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
          dWdPA= dWdU*dUdPA
          dWdPB= dWdU*dUdPB
          dWdGA= dWdU*dUdGA
          dWdGB= dWdU*dUdGB
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          D1F(dRA) = D1F(dRA)+EUEGPA*W + EUEG*dWdPA 
          D1F(dGA) = D1F(dGA)+EUEG*dWdGA 
          D1F(dRB) = D1F(dRB)+EUEGPB*W + EUEG*dWdPB  
          D1F(dGB) = D1F(dGB)+EUEG*dWdGB
         ENDIF
c        ENDDO
        
        RETURN
        END 


      Subroutine m06css(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG,ijzy)
      Implicit none
C
C     Compute the same-spin part of the m06 correlation functional for one grid
C     point and one spin-case.
C
C
      integer ijzy
      REAL*8 PX, GX, TX, F, FP, FG, FT, DTol
      REAL*8 EUEG, Chi, EUEGP, ChiP, ChiG
      REAL*8 Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      REAL*8 ss, sss0,sss1, sss2, sss3, sss4, Css
      REAL*8 Pi, Pi34, F13, F23, F43, F53, F83, F113
      REAL*8 RS, FDUEG, D, Fscc, RSP, dFsccP, dFsccG
      REAL*8 E, W, U, dFsccT, dUdChi, dWdU, dWdP, dWdG 
      REAL*8 d2LdSS,d2LdSZ,d2LdZZ,PotLC,dLdS,dLdZ
     

      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  Css/0.06d0/
C
c      DTol=1.0D-7
      ss=1.0
      if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
         sss0=  5.349466D-01
         sss1=  5.396620D-01
         sss2=  -3.161217D+01
         sss3=  5.149592D+01
         sss4=  -2.919613D+01
      elseif (ijzy.eq.2) then
C     Parameters for M06-HF Correlation
         sss0=  1.023254D-01
         sss1=  -2.453783D+00
         sss2=  2.913180D+01
         sss3=  -3.494358D+01
         sss4=  2.315955D+01
      elseif (ijzy.eq.3) then
C     Parameters for M06 Correlation
         sss0=  5.094055D-01
         sss1=  -1.491085D+00
         sss2=  1.723922D+01
         sss3=  -3.859018D+01
         sss4=  2.845044D+01
      elseif (ijzy.eq.4) then
C     Parameters for M06-2X Correlation
         sss0=  3.097855D-01
         sss1=  -5.528642D+00
         sss2=  1.347420D+01
         sss3=  -3.213623D+01
         sss4=  2.846742D+01
      endif
      
      If ((PX.le.DTol))  then
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
      else
        Pi = F4*ATan(F1)
        Pi34 = F3 / (F4*Pi)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        FDUEG = (F3/F5)*(F6*Pi*Pi)**F23
        RS = (Pi34/PX) ** F13
        Call lsdac(RS,F1,PotLC,dLdS,dLdZ)
        EUEG = PX*PotLC
        D = TX - Pt25*GX/PX
C        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        U = Css*Chi/(F1 + Css*Chi)
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)))
        Fscc=D/TX
        E = Fscc*W*EUEG
        F = E*ss
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        dFsccP=Pt25*GX/(TX*PX**2)
        dFsccG=-Pt25/(TX*PX)
        dFsccT=Pt25*GX/(PX*TX**2)
        dUdChi=Css/((F1+Css*Chi)**2)
        dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4))
        dWdP=dWdU*dUdChi*ChiP
        dWdG=dWdU*dUdChi*ChiG 
        EUEGP = PotLC + PX*dLdS*RSP
        FP = ss*(dFsccP*W*EUEG 
     $                 + Fscc*dWdP*EUEG
     $                 + Fscc*W*EUEGP)
        FG = ss*(dFsccG*W*EUEG
     $                 + Fscc*dWdG*EUEG)

        FT = ss*(dFsccT*W*EUEG)
       Endif

       Return
       End


