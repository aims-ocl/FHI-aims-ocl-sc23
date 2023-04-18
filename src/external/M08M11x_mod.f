!      Subroutine M08M11x(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
      Subroutine M08M11x(F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy)
************************************************************************
*                                                                      *
*  M08M11x evaluates the exchange part of the M08 and M11 suite of     *
*  functionals on the grid.                                            *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  Ref: (a) Zhao, Y.  and Truhlar, D. G. JCTC, 2008, 4 , 1849          *
*       (b) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2011, 2, 2810  *
*       (c) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2012, 3, 117   *
*                                                                      *
*       ijzy - 1 M08-HX (a)                                            *
*       ijzy - 2 M08-SO (a)                                            *
*       ijzy - 3 M11 (b)                                               *
*       ijzy - 4 M11-L (c)                                             *
*                                                                      *
*  OUTPUT:                                                             *
*     F      - Functional values                                       *
*     D1F    - First derivatives with respect to RA, RB, GA, GB        *
*              TA, TB                                                  *
*                                                                      *
*  INPUT:                                                              *
*     RA,B   - Spin densities                                          *
*     D1RA,B - Spin density gradients                                  *
*     TA,B   - Spin kinetic energy densities                           *
*     NGrid  - number of grids                                         *
*                                                                      *
*  RP (11/11), YZ (12/08)                                              *
*                                                                      *
************************************************************************
      Implicit Real*8(A-H,O-Z)
!      INTEGER NGrid
!      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
!     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
      REAL*8  F,D1F(7),RA,RB,D1RA(3),D1RB(3),TA,TB
      Real*8 kapa,kapas,Mu,Mus
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGc
      Logical UseLC
      Save F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11
      DATA F0/0.0D+00/,  F1/1.0D+00/,  F2/2.0D+00/,
     $     F3/3.0D+00/,  F4/4.0D+00/,  F5/5.0D+00/,
     $     F6/6.0D+00/,  F7/7.0D+00/,  F8/8.0D+00/,
     $     F9/9.0D+00/,  F10/10.0D+00/,F11/11.0D+00/
      parameter( pi = 3.1415926535897932384626433832795d0 )
      dRA = 1
      dRB = 2
      dGA = 3
      dGB = 4
      dGC = 5
      dTA = 6
      dTB = 7
c
      DTol=1.0d-10
c
      F1o3 = F1/F3 
      F2o3 = F2/F3
      F3o5 = F3/F5
      F4o3 = F4/F3 
      F5o3 = F5/F3
      F48 = 48.0d0
      F81 = 81.0d0
      Ax = -(F3/F2) * (F4o3*Pi)**(-F1o3) 
C     RPBE parameters
      Mus = F10/F81
      kapas = 0.552d0
C     PBE parameters 
      Mu = 0.21951d0
      kapa = 0.804d0
C      
        ct00= 0D+00
        ct01= 0D+00
        ct02= 0D+00
        ct03= 0D+00
        ct04= 0D+00
        ct05= 0D+00
        ct06= 0D+00
        ct07= 0D+00
        ct08= 0D+00
        ct09= 0D+00
        ct10= 0D+00
        ct11= 0D+00
C
        dt00= 0D+00
        dt01= 0D+00
        dt02= 0D+00
        dt03= 0D+00
        dt04= 0D+00
        dt05= 0D+00
        dt06= 0D+00
        dt07= 0D+00
        dt08= 0D+00
        dt09= 0D+00
        dt10= 0D+00
        dt11= 0D+00
C
      if (ijzy.eq.1) then
C     Parameters for M08-HX
        at00=  1.3340172D+00
        at01= -9.4751087D+00
        at02= -1.2541893D+01
        at03=  9.1369974D+00
        at04=  3.4717204D+01
        at05=  5.8831807D+01
        at06=  7.1369574D+01
        at07=  2.3312961D+01
        at08=  4.8314679D+00
        at09= -6.5044167D+00
        at10= -1.4058265D+01
        at11=  1.2880570D+01

        bt00= -8.5631823D-01
        bt01=  9.2810354D+00
        bt02=  1.2260749D+01
        bt03= -5.5189665D+00
        bt04= -3.5534989D+01
        bt05= -8.2049996D+01
        bt06= -6.8586558D+01
        bt07=  3.6085694D+01
        bt08= -9.3740983D+00
        bt09= -5.9731688D+01
        bt10=  1.6587868D+01
        bt11=  1.3993203D+01
C
        UseLC=.False.
C
       elseif (ijzy.eq.2) then
C     Parameters for M08-SO
        at00= -3.4888428D-01
        at01= -5.8157416D+00
        at02=  3.7550810D+01
        at03=  6.3727406D+01
        at04= -5.3742313D+01
        at05= -9.8595529D+01
        at06=  1.6282216D+01
        at07=  1.7513468D+01
        at08= -6.7627553D+00
        at09=  1.1106658D+01
        at10=  1.5663545D+00
        at11=  8.7603470D+00

        bt00=  7.8098428D-01
        bt01=  5.4538178D+00
        bt02= -3.7853348D+01
        bt03= -6.2295080D+01
        bt04=  4.6713254D+01
        bt05=  8.7321376D+01
        bt06=  1.6053446D+01
        bt07=  2.0126920D+01
        bt08= -4.0343695D+01
        bt09= -5.8577565D+01
        bt10=  2.0890272D+01
        bt11=  1.0946903D+01
C
        UseLC=.False.
C
      elseif (ijzy.eq.3) then
C     Parameters for M11
        at00= -0.18399900D+00
        at01= -1.39046703D+01
        at02=  1.18206837D+01
        at03=  3.10098465D+01
        at04= -5.19625696D+01
        at05=  1.55750312D+01
        at06= -6.94775730D+00
        at07= -1.58465014D+02
        at08= -1.48447565D+00
        at09=  5.51042124D+01
        at10= -1.34714184D+01
        at11=  0.00000000D+00

        bt00=  0.75599900D+00
        bt01=  1.37137944D+01
        bt02= -1.27998304D+01
        bt03= -2.93428814D+01
        bt04=  5.91075674D+01
        bt05= -2.27604866D+01
        bt06= -1.02769340D+01
        bt07=  1.64752731D+02
        bt08=  1.85349258D+01
        bt09= -5.56825639D+01
        bt10=  7.47980859D+00
        bt11=  0.00000000D+00
C
        UseLC=.True.
        Emu =0.25D+00
C
      elseif (ijzy.eq.4) then
C     Parameters for M11-L
        at00=  8.121131D-01
        at01=  1.738124D+01
        at02=  1.154007D+00
        at03=  6.869556D+01
        at04=  1.016864D+02
        at05= -5.887467D+00
        at06=  4.517409D+01
        at07= -2.773149D+00
        at08= -2.617211D+01
        at09=  0.000000D+00
        at10=  0.000000D+00 
        at11=  0.000000D+00
C
        bt00=  1.878869D-01
        bt01= -1.653877D+01
        bt02=  6.755753D-01
        bt03= -7.567572D+01
        bt04= -1.040272D+02
        bt05=  1.831853D+01
        bt06= -5.573352D+01
        bt07= -3.520210D+00
        bt08=  3.724276D+01
        bt09=  0.000000D+00
        bt10=  0.000000D+00
        bt11=  0.000000D+00
C
        ct00= -4.386615D-01
        ct01= -1.214016D+02
        ct02= -1.393573D+02
        ct03= -2.046649D+00
        ct04=  2.804098D+01
        ct05= -1.312258D+01
        ct06= -6.361819D+00
        ct07= -8.055758D-01
        ct08=  3.736551D+00
        ct09=  0.000000D+00
        ct10=  0.000000D+00
        ct11=  0.000000D+00
C
        dt00=  1.438662D+00
        dt01=  1.209465D+02
        dt02=  1.328252D+02
        dt03=  1.296355D+01
        dt04=  5.854866D+00
        dt05= -3.378162D+00
        dt06= -4.423393D+01
        dt07=  6.844475D+00
        dt08=  1.949541D+01
        dt09=  0.000000D+00
        dt10=  0.000000D+00
        dt11=  0.000000D+00
C
        UseLC=.True.
        Emu =0.25D+00
C
      endif
c
!      DO i = 1,NGrid
         IF ((RA.gt.DTol).and.(TA.gt.DTol)) THEN
            rhoo  = RA
            rho43 = rhoo**F4o3
            rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c
            tauN = TA
            TauUEG=F3o5*((F6*pi*pi)**F2o3)*rho53
            TSIG =TAUUEG/TAUN
            Wsig =(Tsig - F1)/(Tsig + F1)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig1 =(at00    + at01*W1 + at02*W2 + at03*W3
     $            + at04*W4 + at05*W5 + at06*W6 + at07*W7
     $            + at08*W8 + at09*W9 + at10*W10+ at11*W11)
            Fsig2 =(bt00    + bt01*W1 + bt02*W2 + bt03*W3
     $            + bt04*W4 + bt05*W5 + bt06*W6 + bt07*W7
     $            + bt08*W8 + bt09*W9 + bt10*W10+ bt11*W11)
            Fsig3 =(ct00    + ct01*W1 + ct02*W2 + ct03*W3
     $            + ct04*W4 + ct05*W5 + ct06*W6 + ct07*W7
     $            + ct08*W8 + ct09*W9 + ct10*W10+ ct11*W11)
            Fsig4 =(dt00    + dt01*W1 + dt02*W2 + dt03*W3
     $            + dt04*W4 + dt05*W5 + dt06*W6 + dt07*W7
     $            + dt08*W8 + dt09*W9 + dt10*W10+ dt11*W11)


           Gam =  D1RA(1)**F2 + D1RA(2)**F2 + D1RA(3)**F2
           Gam12 = dsqrt(Gam)
           x = Gam12/rho43
           s = x/(F48*Pi*Pi)**F1o3
           y = s*s
           Deno = (F1 + Mu*y/kapa)
           fx1=F1+kapa*(F1-F1/Deno)
           fx2=F1+kapas*(F1-Exp(-Mus*y/kapas))
c
c lrclsda is a routine to calculate the short-range part of thr
c local spin density approximation (LSDA) to the exchange energy.
c This routine is provided within MFM at http://comp.chem.umn.edu/mfm/
c
           If(UseLC) then
             CALL LRCLSDA(EMU,RHOO,ElSR,PDUM)
             ElLR = Ax*Rho43-ElSR
           else
             ElSR = Ax*Rho43
             ElLR = F0
           endIf
           GGA1 = ElSR*fx1
           GGA2 = ElSR*fx2
           GGA3 = ElLR*fx1
           GGA4 = ElLR*fx2
C
          F=F+GGA1*Fsig1+GGA2*Fsig2
     $             +GGA3*Fsig3+GGA4*Fsig4
c
c     functional derivatives
c
            dydRho = -(F8/F3) * y/rhoo
            dydG = y/Gam
C
            dfx1dy = Mu*(F1/(Deno*Deno)) 
            dfx1dRho = dfx1dy*dydRho
            dfx1dG = dfx1dy*dydG
C
            dfx2dy = Mus*Exp(-Mus*y/kapas)
            dfx2dRho = dfx2dy*dydRho
            dfx2dG = dfx2dy*dydG

            dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
     $                    + F4*at04*W3 + F5*at05*W4
     $                    + F6*at06*W5 + F7*at07*W6
     $                    + F8*at08*W7 + F9*at09*W8
     $                    + F10*at10*W9+F11*at11*W10)
            dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
     $                    + F4*bt04*W3 + F5*bt05*W4
     $                    + F6*bt06*W5 + F7*bt07*W6
     $                    + F8*bt08*W7 + F9*bt09*W8
     $                    + F10*Bt10*W9+F11*Bt11*W10)
            dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
     $                    + F4*ct04*W3 + F5*ct05*W4
     $                    + F6*ct06*W5 + F7*ct07*W6
     $                    + F8*ct08*W7 + F9*ct09*W8
     $                    + F10*ct10*W9+F11*ct11*W10)
            dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
     $                    + F4*dt04*W3 + F5*dt05*W4
     $                    + F6*dt06*W5 + F7*dt07*W6
     $                    + F8*dt08*W7 + F9*dt09*W8
     $                    + F10*dt10*W9+F11*dt11*W10)

            dWdT = F2/((F1 + Tsig)**F2)
            dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
            dTdTau = -TauUEG/tauN**F2
C
           If(UseLC) then
             dElSRdR = PDUM
             dElLRdR = Ax*F4o3*Rho13-PDUM
           else
             dElSRdR=Ax*F4o3*Rho13
             dElLRdR=F0
           endIf  
            dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
            dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
            dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
            dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
            dF1dR = dF1dW*dWdT*dTdR
            dF1dTau=dF1dW*dWdT*dTdTau
            dF2dR = dF2dW*dWdT*dTdR
            dF2dTau=dF2dW*dWdT*dTdTau
            dF3dR = dF3dW*dWdT*dTdR
            dF3dTau=dF3dW*dWdT*dTdTau
            dF4dR = dF4dW*dWdT*dTdR
            dF4dTau=dF4dW*dWdT*dTdTau
c
            dGGA1dG = ElSR*dfx1dG
            dGGA2dG = ElSR*dfx2dG
            dGGA3dG = ElLR*dfx1dG
            dGGA4dG = ElLR*dfx2dG
C dF/dRhoa
            D1F(dRA) = D1F(dRA) +dGGA1dR*Fsig1 + GGA1*dF1dR
     $                              +dGGA2dR*Fsig2 + GGA2*dF2dR
     $                              +dGGA3dR*Fsig3 + GGA3*dF3dR
     $                              +dGGA4dR*Fsig4 + GGA4*dF4dR
C dF/dGammaaa
            D1F(dGA) = D1F(dGA) +dGGA1dG*Fsig1 + dGGA2dG*Fsig2
     $                              +dGGA3dG*Fsig3 + dGGA4dG*Fsig4
C dF/dTaua
            D1F(dTA) = D1F(dTA) +GGA1*dF1dTau + GGA2*dF2dTau
     $                              +GGA3*dF3dTau + GGA4*dF4dTau
         ENDIF
c
c beta component
c
         IF ((RB.gt.DTol).and.(TB.gt.DTol)) THEN
            rhoo  = RB
            rho43 = rhoo**F4o3
            rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c
            tauN = TB
            TauUEG=F3o5*((F6*pi*pi)**F2o3)*rho53
            TSIG =TAUUEG/TAUN
            Wsig =(Tsig - F1)/(Tsig + F1)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig1 =(at00    + at01*W1 + at02*W2 + at03*W3
     $            + at04*W4 + at05*W5 + at06*W6 + at07*W7
     $            + at08*W8 + at09*W9 + at10*W10+ at11*W11)
            Fsig2 =(bt00    + bt01*W1 + bt02*W2 + bt03*W3
     $            + bt04*W4 + bt05*W5 + bt06*W6 + bt07*W7
     $            + bt08*W8 + bt09*W9 + bt10*W10+ bt11*W11)
            Fsig3 =(ct00    + ct01*W1 + ct02*W2 + ct03*W3
     $            + ct04*W4 + ct05*W5 + ct06*W6 + ct07*W7
     $            + ct08*W8 + ct09*W9 + ct10*W10+ ct11*W11)
            Fsig4 =(dt00    + dt01*W1 + dt02*W2 + dt03*W3
     $            + dt04*W4 + dt05*W5 + dt06*W6 + dt07*W7
     $            + dt08*W8 + dt09*W9 + dt10*W10+ dt11*W11)

           Gam =  D1RB(1)**F2 + D1RB(2)**F2 + D1RB(3)**F2
           Gam12 = dsqrt(Gam)
           x = Gam12/rho43
           s = x/(F48*Pi*Pi)**F1o3
           y = s*s
           Deno = (F1 + Mu*y/kapa)
           fx1=F1+kapa*(F1-F1/Deno)
           fx2=F1+kapas*(F1-Exp(-Mus*y/kapas))
c
c lrclsda is a routine to calculate the short-range part of thr
c local spin density approximation (LSDA) to the exchange energy.
c This routine is provided within MFM at http://comp.chem.umn.edu/mfm/
c
           If(UseLC) then
             CALL LRCLSDA(EMU,RHOO,ElSR,PDUM)
             ElLR = Ax*Rho43-ElSR
           else
             ElSR = Ax*Rho43
             ElLR = F0
           endIf
           GGA1 = ElSR*fx1
           GGA2 = ElSR*fx2
           GGA3 = ElLR*fx1
           GGA4 = ElLR*fx2
C
          F=F+GGA1*Fsig1+GGA2*Fsig2
     $               +GGA3*Fsig3+GGA4*Fsig4
c
c     functional derivatives
c
            dydRho = -(F8/F3) * y/rhoo
            dydG = y/Gam
C
            dfx1dy = Mu*(F1/(Deno*Deno)) 
            dfx1dRho = dfx1dy*dydRho
            dfx1dG = dfx1dy*dydG
C
            dfx2dy = Mus*Exp(-Mus*y/kapas)
            dfx2dRho = dfx2dy*dydRho
            dfx2dG = dfx2dy*dydG
C
            dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
     $                    + F4*at04*W3 + F5*at05*W4
     $                    + F6*at06*W5 + F7*at07*W6
     $                    + F8*at08*W7 + F9*at09*W8
     $                    + F10*at10*W9+F11*at11*W10)
            dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
     $                    + F4*bt04*W3 + F5*bt05*W4
     $                    + F6*bt06*W5 + F7*bt07*W6
     $                    + F8*bt08*W7 + F9*bt09*W8
     $                    + F10*Bt10*W9+F11*Bt11*W10)
            dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
     $                    + F4*ct04*W3 + F5*ct05*W4
     $                    + F6*ct06*W5 + F7*ct07*W6
     $                    + F8*ct08*W7 + F9*ct09*W8
     $                    + F10*ct10*W9+F11*ct11*W10)
            dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
     $                    + F4*dt04*W3 + F5*dt05*W4
     $                    + F6*dt06*W5 + F7*dt07*W6
     $                    + F8*dt08*W7 + F9*dt09*W8
     $                    + F10*dt10*W9+F11*dt11*W10)

            dWdT = F2/((F1 + Tsig)**F2)
            dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
            dTdTau = -TauUEG/tauN**F2
C
           If(UseLC) then
             dElSRdR = PDUM
             dElLRdR = Ax*F4o3*Rho13-PDUM
           else
             dElSRdR=Ax*F4o3*Rho13
             dElLRdR=F0
           endIf  
            dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
            dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
            dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
            dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
            dF1dR = dF1dW*dWdT*dTdR
            dF1dTau=dF1dW*dWdT*dTdTau
            dF2dR = dF2dW*dWdT*dTdR
            dF2dTau=dF2dW*dWdT*dTdTau
            dF3dR = dF3dW*dWdT*dTdR
            dF3dTau=dF3dW*dWdT*dTdTau
            dF4dR = dF4dW*dWdT*dTdR
            dF4dTau=dF4dW*dWdT*dTdTau
c
            dGGA1dG = ElSR*dfx1dG
            dGGA2dG = ElSR*dfx2dG
            dGGA3dG = ElLR*dfx1dG
            dGGA4dG = ElLR*dfx2dG
C dF/dRhob
            D1F(dRB) = D1F(dRB) +dGGA1dR*Fsig1 + GGA1*dF1dR
     $                              +dGGA2dR*Fsig2 + GGA2*dF2dR
     $                              +dGGA3dR*Fsig3 + GGA3*dF3dR
     $                              +dGGA4dR*Fsig4 + GGA4*dF4dR
C dF/dGammabb
            D1F(dGB) = D1F(dGB) +dGGA1dG*Fsig1 + dGGA2dG*Fsig2
     $                              +dGGA3dG*Fsig3 + dGGA4dG*Fsig4
C dF/dTaub
            D1F(dTB) = D1F(dTB) +GGA1*dF1dTau + GGA2*dF2dTau
     $                              +GGA3*dF3dTau + GGA4*dF4dTau
         ENDIF
!      enddo
      Return
      End

