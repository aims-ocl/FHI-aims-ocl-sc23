!      Subroutine M08M11c(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
      Subroutine M08M11c(F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy)
************************************************************************
*                                                                      *
*  M08M11c evaluates the correlation part of the M08 and M11 suite of  *
*  functionals on the grid.                                            *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  OUTPUT:                                                             *
*     F      - Functional values                                       *
*     D1F    - First derivatives with respect to RA, RB, GA, GB        *
*              TA, TB                                                  *
*                                                                      *
*  INPUT:                                                              *
*                                                                      *
*       ijzy - 1 M08-HX                                                *
*       ijzy - 2 M08-SO                                                *
*       ijzy - 3 M11                                                   *
*       ijzy - 4 M11-L                                                 *
*       ijzy - 5 MN12-L                                                *
*       ijzy - 6 N12-SX                                                *
*                                                                      *
*     RA,B   - Spin densities                                          *
*     D1RA,B - Spin density gradients                                  *
*     TA,B   - Spin kinetic energy densities                           *
*     NGrid  - number of grids                                         *
*                                                                      *
*  RP (09/12), YZ (12/08)                                              *
*                                                                      *
************************************************************************
      Implicit Real*8(A-H,O-Z)
      Real*8 LSDA
!      INTEGER NGrid
!      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
!     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
      REAL*8  F,D1F(7),RA,RB,D1RA(3),D1RB(3),TA,TB
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC
      Save F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11
      DATA F1/1.0D+00/,  F2/2.0D+00/,  F3/3.0D+00/,
     $     F4/4.0D+00/,  F5/5.0D+00/,  F6/6.0D+00/,
     $     F7/7.0D+00/,  F8/8.0D+00/,  F9/9.0D+00/, 
     $     F10/10.0D+00/,F11/11.0D+00/
      parameter( pi = 3.1415926535897932384626433832795d0 )
      dRA = 1
      dRB = 2
      dGA = 3
      dGB = 4
      dGC = 5
      dTA = 6
      dTB = 7
C
      DTol   = 1.0d-10
C
      F1o3 = F1/F3
      F2o3 = F2/F3
      F5o3 = F5/F3
C
      Pi34   = F3/(F4*Pi)
C
      if (ijzy.eq.1) then
C     Parameters for M08-HX
        at0=    1.0000000D+00
        at1=   -4.0661387D-01
        at2=   -3.3232530D+00
        at3=    1.5540980D+00
        at4=    4.4248033D+01
        at5=   -8.4351930D+01
        at6=   -1.1955581D+02
        at7=    3.9147081D+02
        at8=    1.8363851D+02
        at9=   -6.3268223D+02
        at10=  -1.1297403D+02
        at11=   3.3629312D+02

        bt0=    1.3812334D+00
        bt1=   -2.4683806D+00
        bt2=   -1.1901501D+01
        bt3=   -5.4112667D+01
        bt4=    1.0055846D+01
        bt5=    1.4800687D+02
        bt6=    1.1561420D+02
        bt7=    2.5591815D+02
        bt8=    2.1320772D+02
        bt9=   -4.8412067D+02
        bt10=  -4.3430813D+02
        bt11=   5.6627964D+01
       elseif (ijzy.eq.2) then
C     Parameters for M08-SO 
        at0=    1.0000000D+00
        at1=    0.0000000D+00
        at2=   -3.9980886D+00
        at3=    1.2982340D+01
        at4=    1.0117507D+02
        at5=   -8.9541984D+01
        at6=   -3.5640242D+02
        at7=    2.0698803D+02
        at8=    4.6037780D+02
        at9=   -2.4510559D+02
        at10=  -1.9638425D+02
        at11=   1.1881459D+02

        bt0=    1.0000000D+00
        bt1=   -4.4117403D+00
        bt2=   -6.4128622D+00
        bt3=    4.7583635D+01
        bt4=    1.8630053D+02
        bt5=   -1.2800784D+02
        bt6=   -5.5385258D+02
        bt7=    1.3873727D+02
        bt8=    4.1646537D+02
        bt9=   -2.6626577D+02
        bt10=   5.6676300D+01
        bt11=   3.1673746D+02
       elseif (ijzy.eq.3) then
C     Parameters for M11
        at0=   1.0000000D+00
        at1=   0.0000000D+00
        at2=  -3.8933250D+00
        at3=  -2.1688455D+00
        at4=   9.3497200D+00
        at5=  -1.9845140D+01
        at6=   2.3455253D+00
        at7=   7.9246513D+01
        at8=   9.6042757D+00
        at9=  -6.7856719D+01
        at10= -9.1841067D+00
        at11=  0.0000000D+00

        bt0=   7.2239798D-01
        bt1=   4.3730564D-01
        bt2=  -1.6088809D+01
        bt3=  -6.5542437D+01
        bt4=   3.2057230D+01
        bt5=   1.8617888D+02
        bt6=   2.0483468D+01
        bt7=  -7.0853739D+01
        bt8=   4.4483915D+01
        bt9=  -9.4484747D+01
        bt10= -1.1459868D+02
        bt11=  0.0000000D+00
       elseif (ijzy.eq.4) then
C     Parameters for M11-L
        at0=   1.000000D+00
        at1=   0.000000D+00
        at2=   2.750880D+00
        at3=  -1.562287D+01
        at4=   9.363381D+00
        at5=   2.141024D+01
        at6=  -1.424975D+01
        at7=  -1.134712D+01
        at8=   1.022365D+01
        at9=   0.000000D+00
        at10=  0.000000D+00
        at11=  0.000000D+00
C
        bt0=   1.000000D+00
        bt1=  -9.082060D+00
        bt2=   6.134682D+00
        bt3=  -1.333216D+01
        bt4=  -1.464115D+01
        bt5=   1.713143D+01
        bt6=   2.480738D+00
        bt7=  -1.007036D+01
        bt8=  -1.117521D-01
        bt9=   0.000000D+00
        bt10=  0.000000D+00
        bt11=  0.000000D+00
       elseif (ijzy.eq.5) then
C     Parameters for MN12-L
        at00=  8.844610D-01
        at01= -2.202279D-01
        at02=  5.701372D+00
        at03= -2.562378D+00
        at04= -9.646827D-01
        at05=  1.982183D-01
        at06=  1.019976D+01
        at07=  9.789352D-01
        at08= -1.512722D+00
        at09=  0.000000D+00
        at10=  0.000000D+00
        at11=  0.000000D+00
C
        bt00=  5.323948D-01
        bt01= -5.831909D+00
        bt02=  3.882386D+00
        bt03=  5.878488D+00
        bt04=  1.493228D+01
        bt05= -1.374636D+01
        bt06= -8.492327D+00
        bt07= -2.486548D+00
        bt08= -1.822346D+01
        bt09=  0.000000D+00
        bt10=  0.000000D+00
        bt11=  0.000000D+00
       elseif (ijzy.eq.6) then
C     Parameters for MN12-SX
        at00=  7.171161D-01
        at01= -2.380914D+00
        at02=  5.793565D+00
        at03= -1.243624D+00
        at04=  1.364920D+01
        at05= -2.110812D+01
        at06= -1.598767D+01
        at07=  1.429208D+01
        at08=  6.149191D+00
        at09=  0.000000D+00
        at10=  0.000000D+00
        at11=  0.000000D+00
C
        bt00=  4.663699D-01
        bt01= -9.110685D+00
        bt02=  8.705051D+00
        bt03= -1.813949D+00
        bt04= -4.147211D-01
        bt05= -1.021527D+01
        bt06=  8.240270D-01
        bt07=  4.993815D+00
        bt08= -2.563930D+01
        bt09=  0.000000D+00
        bt10=  0.000000D+00
        bt11=  0.000000D+00
       endif
      
!      DO i = 1,NGrid
       RhoA = RA
       RhoB = RB
       Rho = RhoA + RhoB
       TauA = TA/F2
       TauB = TB/F2
       Tau = TauA + TauB

       If(Rho.gt.DTol.and.Tau.gt.DTol) then
        RS = (Pi34/Rho)**F1o3
        Zeta = (RhoA-RhoB)/Rho
        TauUEG=F3*(F3*Pi*Pi)**(F2o3)*Rho**(F5o3)/F10
        Tsig =TauUEG/Tau
        Wsig =(Tsig - F1)/(Tsig + F1)
        Fsig1=(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))

        Fsig2=(bt0 + Wsig*(bt1 + Wsig*(bt2 + Wsig*(bt3 + Wsig*(
     &            bt4 + Wsig*(bt5 + Wsig*(bt6 + Wsig*(bt7 + Wsig*(
     &            bt8 + Wsig*(bt9 + Wsig*(bt10+Wsig*bt11)))))))))))

        Y = (D1RA(1) + D1RB(1))**F2
     $      + (D1RA(2) + D1RB(2))**F2
     $      + (D1RA(3) + D1RB(3))**F2
        GRho = Sqrt(Y)
c       
c      lsdac is a subroutine to evaluate the Perdew-Wang-91 correlation functional 
c      local spin density approximation (LSDA) to the correlation energy of a uniform 
c      electron gas. (Phys. Rev. B 45, 13244 (1992)). Users should provid their own
c      for this LSDA correlation functional or they may find this routine on Kieron 
c      Burke's Web site at http://www.chem.uci.edu/~kieron/dftold2/pubs/PBE.asc
c
c        Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ)
c        LSDA = Rho*PotLC
c
c      PBEH0 is a subroutine to evaluate the H0 term in the PBE correlation functional
c      (Phys. Rev. Lett. 77, 3865 - 3868 (1996)) Users should provid their own
c      for this H0 subroutine or they may find this routine on Kieron
c      Burke's Web site at http://www.chem.uci.edu/~kieron/dftold2/pubs/PBE.asc
c
        Call PBEH0(Rho,GRho,Zeta,PotLC,dLdS,dLdZ,H,dHdR,dHdG,dHdZ)

c        Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ)

        LSDA = Rho*PotLC
        GGA = Rho*H 
        E1 = LSDA*Fsig1
        E2 = GGA*Fsig2
        F = F + E1 +E2
        
c
c     functional derivatives
c
         RSP = -RS/(F3*Rho)
         dZdA = (F1-Zeta)/Rho
         dZdB = (-F1-Zeta)/Rho
         dLdRA = dLdS*RSP + dLdZ*dZdA
         dLdRB = dLdS*RSP + dLdZ*dZdB
         dF1dW=( at1 + Wsig*(F2  *at2 + Wsig*(F3*at3 + Wsig*(
     &            F4 *at4 + Wsig*(F5 *at5 + Wsig*(F6  *at6 + Wsig*(
     &            F7*at7 + Wsig*(F8*at8 + Wsig*(F9 *at9 + Wsig*(
     &            F10  *at10+ Wsig*F11*at11))))))))))
         dF2dW=( bt1 + Wsig*(F2  *bt2 + Wsig*(F3*bt3 + Wsig*(
     &            F4 *bt4 + Wsig*(F5 *bt5 + Wsig*(F6  *bt6 + Wsig*(
     &            F7*bt7 + Wsig*(F8*bt8 + Wsig*(F9 *bt9 + Wsig*(
     &            F10  *bt10+ Wsig*F11*bt11))))))))))
         dWdT = F2/((F1 + Tsig)**F2)
         dTdR = Tsig*F5/(F3*Rho) 
         dTdTau = -Tsig/Tau
         dF1dR = dF1dW*dWdT*dTdR
         dF1dTau=dF1dW*dWdT*dTdTau
         dF2dR = dF2dW*dWdT*dTdR
         dF2dTau=dF2dW*dWdT*dTdTau
         dLDdRA = PotLC + Rho*dLdRA
         dLDdRB = PotLC + Rho*dLdRB
         dHdRA = dHdR + dHdZ*dZdA
         dHdRB = dHdR + dHdZ*dZdB
         dGRhodY = F1/(F2*GRho)
         dHdY = dHdG * dGRhodY
         dHdGA = dHdY
         dHdGC = dHdY*F2  
         dGGAdRA = H + Rho*dHdRA
         dGGAdRB = H + Rho*dHdRB
         dGGAdGA = Rho*dHdGA
         dGGAdGB = dGGAdGA
         dGGAdGC = Rho*dHdGC
C
         dE1dRA = dLDdRA*Fsig1 + LSDA*dF1dR
         dE1dRB = dLDdRB*Fsig1 + LSDA*dF1dR
         dE1dKA = LSDA*dF1dTau
         dE1dKB = dE1dKA
C
         dE2dRA = dGGAdRA*Fsig2 + GGA*dF2dR
         dE2dRB = dGGAdRB*Fsig2 + GGA*dF2dR 
         dE2dKA = GGA*dF2dTau
         dE2dKB = dE2dKA
         dE2dGA = dGGAdGA*Fsig2
         dE2dGB = dGGAdGB*Fsig2
         dE2dGC = dGGAdGC*Fsig2   

         D1F(dRA)=D1F(dRA) + dE1dRA + dE2dRA 
         D1F(dRB)=D1F(dRB) + dE1dRB + dE2dRB  
         D1F(dTA)=D1F(dTA) + (dE1dKA + dE2dKA)/F2
         D1F(dTB)=D1F(dTB) + (dE1dKB + dE2dKB)/F2
         D1F(dGA)=D1F(dGA) + dE2dGA
         D1F(dGB)=D1F(dGB) + dE2dGB
C     GC is the dot product of the vectors D1RA and D1RB
         D1F(dGC)=D1F(dGC) + dE2dGC 
       Endif
!      Enddo
      Return
      End

