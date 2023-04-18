 
c      Subroutine M06x(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
      Subroutine M06x(F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy)
************************************************************************
*                                                                      *
*  M06x evaluates the exchange part of the M06 suite of                *
*  functionals on the grid.                                            *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  Ref: (a) Zhao, Y.  and Truhlar, D. G. J. Chem. Phys. 125,           *
*    194101 (2006).                                                    *
*       (b) Y. Zhao and D. G. Truhlar, J. Phys. Chem. A (2006),        *
*    110(49),13126-13130.                                              *
*                                                                      *
*       ijzy - 1 M06-L                                                 *
*       ijzy - 2 M06-HF                                                *
*       ijzy - 3 M06                                                   *
*       ijzy - 4 M06-2X                                                *
*                                                                      *
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
*                                                                      *
*  YZ (12/08)                                                          *
*                                                                      *
************************************************************************
      Implicit Real*8(A-H,O-Z)
c      INTEGER NGrid
c      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
c     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
      REAL*8  F,D1F(7),RA,RB,D1RA(3),D1RB(3),TA,TB
      REAL*8 pi
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC
      integer i,j,ijzy
      REAL*8 at1, at2, at3, at4, at5, at6, at7, at8, at9
      REAL*8 at, at10, at11, C1, C2, fL, fNL, at0
      REAL*8 rrho, rho43, rho13, rhoo, rho53
      REAL*8 Gamma2, Gamma
      REAL*8 TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      REAL*8 W7, W8, W9, W10, W11, Fsig
      REAL*8 tauu,DTol
      REAL*8 F83, F23, F53, F1o3
      REAL*8 F1o4, F2o3, F3o2, F4o3, F4o9, F3o5
      REAL*8 Nine, F10, F11
      REAL*8  Ax,  x, x2, En, Ed, E, dE, dEn, dEd
      REAL*8 dFdW, dWdT, dTdR, dTdTau, dGGAdR, dFdR
      REAL*8 dFdTau, dGGAdG
      parameter( pi = 3.1415926535897932384626433832795d0 )
      parameter (F1o3=1.d0/3.d0, F1o4=1.d0/4.d0, F2o3=2.d0/3.d0, 
     &             F3o2=3.d0/2.d0)
      parameter (F4o3=4.d0/3.d0, F4o9=4.0d0/9.0d0, F3o5=3.d0/5.d0)
      parameter (F83=8.0d0/3.0d0, F23=2.0d0/3.d0, F53=5.d0/3.d0)
      parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &             Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &             Eight=8.0d0, Nine=9.0d0,F10=10.0d0, F11=11.d0)

      DTol=1.0d-10
C     Global assignments
      dRA = 1
      dRB = 2
      dGA = 3
      dGB = 4
      dGC = 5
      dTA = 6
      dTB = 7
      if (ijzy.eq.1) then
C     Parameters for M06-L
        at0=    3.987756D-01
        at1=    2.548219D-01
        at2=    3.923994D-01
        at3=    -2.103655D+00
        at4=    -6.302147D+00
        at5=    1.097615D+01
        at6=    3.097273D+01
        at7=    -2.318489D+01
        at8=    -5.673480D+01
        at9=    2.160364D+01
        at10=   3.421814D+01
        at11=   -9.049762D+00
       elseif (ijzy.eq.2) then
C     Parameters for M06-HF
        at0=    1.179732D-01
        at1=    -1.066708D+00
        at2=    -1.462405D-01
        at3=    7.481848D+00
        at4=    3.776679D+00
        at5=    -4.436118D+01
        at6=    -1.830962D+01
        at7=    1.003903D+02
        at8=    3.864360D+01
        at9=    -9.806018D+01
        at10=   -2.557716D+01
        at11=   3.590404D+01
       elseif (ijzy.eq.3) then
C     Parameters for M06
        at0=    5.877943D-01
        at1=    -1.371776D-01
        at2=    2.682367D-01
        at3=    -2.515898D+00
        at4=    -2.978892D+00
        at5=    8.710679D+00
        at6=    1.688195D+01
        at7=    -4.489724D+00
        at8=    -3.299983D+01
        at9=    -1.449050D+01
        at10=   2.043747D+01
        at11=   1.256504D+01
C The next command is originally commented,
C but it looks strange this way. MR
C Agreed! AJL
      elseif (ijzy.eq.4) then
C     Parameters for M06-2X
        at0=    4.600000D-01
        at1=    -2.206052D-01
        at2=    -9.431788D-02
        at3=    2.164494D+00
        at4=    -2.556466D+00
        at5=    -1.422133D+01
        at6=    1.555044D+01
        at7=    3.598078D+01
        at8=    -2.722754D+01
        at9=    -3.924093D+01
        at10=   1.522808D+01
        at11=   1.522227D+01
      endif
*
c      if (ijzy.LT.4) call VS98x(DTol,1.0d0,F,D1F,RA,RB,D1RA,D1RB,
c     &                            TA,TBxxx,NGrid,ijzy+1)
      if (ijzy.LT.4) call VS98x(DTol,1.0d0,F,D1F,RA,RB,D1RA,D1RB,
     &                            TA,TB,ijzy+1)      

      C1     = 3.36116D-03
      C2     = 4.49267D-03
      
      Ax = -F3o2*(F4o3*PI)**(-F1o3)
C      DO i = 1,NGrid
C         IF ((RA(i).gt.DTol).and.(TA(i).gt.DTol)) THEN
         IF (RA.gt.DTol.and.TA.gt.DTol) THEN
C            rhoo = RA(i)
            rhoo = RA
            rho43 = rhoo**F4o3  
            rrho = 1.0d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53
c
C            tauN = TA(i)
            tauN = TA
            tauu = tauN
            TauUEG=F3o5*((Six*pi*pi)**F2o3)*rho53
            Tsig =TauUEG/tauN
            Wsig =(Tsig-One)/(Tsig+One)
            Fsig=(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))
C         Gamma2 =  D1RA(i,1)**2 + D1RA(i,2)**2 + D1RA(i,3)**2
         Gamma2 =  D1RA(1)**2 + D1RA(2)**2 + D1RA(3)**2
         Gamma = dsqrt(Gamma2)
         x = Gamma/rho43
         x2 = x*x
         En = C1*x2
         Ed = One + C2*x2
         E  = -En/Ed
C         F(i)=F(i)+(Ax+E)*Fsig*rho43
c Interestingly when removing the array structure from above, unless
c one includes the spacing as below around "F = F +" then I seem to
c get erroneous values, ifort-mpi compiler on 4 AMD processor machine
c running GNU/Linux. Just a note for other editors. AJL
         F = F +(Ax+E)*Fsig*rho43
c
c     functional derivatives
c

         dEn   = Two*C1*x
         dEd   = Two*C2*x
         dE    = -(dEn*Ed-En*dEd)/(Ed*Ed)
         dFdW=( at1 + Wsig*(Two  *at2 + Wsig*(Three*at3 + Wsig*(
     &            Four *at4 + Wsig*(Five *at5 + Wsig*(Six  *at6 + Wsig*(
     &            Seven*at7 + Wsig*(Eight*at8 + Wsig*(Nine *at9 + Wsig*(
     &            F10  *at10+ Wsig*F11*at11))))))))))
         dWdT = Two/((One + Tsig)**2)
         dTdR = ((Six*PI*PI)**F2o3)*(rhoo**F2o3)/tauu
         dTdTau = -TauUEG/tauu**2
         dGGAdR = F4o3*rho13*(Ax+(E-x*dE))
         dFdR = dFdW*dWdT*dTdR
         dFdTau=dFdW*dWdT*dTdTau
         dGGAdG =(dE/(Two*Gamma))
*        dF/dRhoa
c            D1F(i,GAB) = 0.0d0
*        dF/dRhoa
c            D1F(i,dRA) = D1F(i,dRA) + dGGAdR*Fsig + (Ax+E)*rho43*dFdR
            D1F(dRA) = D1F(dRA) + dGGAdR*Fsig + (Ax+E)*rho43*dFdR
*        dF/dGammaaa
c            D1F(i,dGA) = D1F(i,dGA) + dGGAdG*Fsig 
            D1F(dGA) = D1F(dGA) + dGGAdG*Fsig
*        dF/dTaua
c            D1F(i,dTA) = D1F(i,dTA) + rho43*(Ax+E)*dFdTau
            D1F(dTA) = D1F(dTA) + rho43*(Ax+E)*dFdTau
         ENDIF
c
c beta component
c
c         IF ((RB(i).gt.DTol).and.(TB(i).gt.DTol)) THEN
         IF (RB.gt.DTol.and.TB.gt.DTol) THEN
c            rhoo = RB(i)
            rhoo = RB
            rho43 = rhoo**F4o3
            rrho = 1.0d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53
c
c            tauN = TB(i)
            tauN = TB
            tauu = tauN
            TauUEG=F3o5*((Six*pi*pi)**F2o3)*rho53
            Tsig =TauUEG/tauN
            Wsig =(Tsig-One)/(Tsig+One)
            Fsig=(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))
c         Gamma2 = D1RB(i,1)**2 + D1RB(i,2)**2 + D1RB(i,3)**2
         Gamma2 = D1RB(1)**2 + D1RB(2)**2 + D1RB(3)**2
         Gamma = dsqrt(Gamma2)
         x = Gamma/rho43
         x2 = x*x
         En = C1*x2
         Ed = One + C2*x2
         E  = -En/Ed
c         F(i)=F(i)+(Ax+E)*Fsig*rho43
         F = F + (Ax+E)*Fsig*rho43
c
c     functional derivatives
c

         dEn   = Two*C1*x
         dEd   = Two*C2*x
         dE    = -(dEn*Ed-En*dEd)/(Ed*Ed)

c I took out first parameter here but it isi just a guess based
C on the other component ! "at" by itself
c is never initialized. Must check the real derivations..... I am not
c worried about derivatives right now. Ugh. MR

c I've just compared this source to Truhlars webpage:
c http://comp.chem.umn.edu/mfm/index.html
c and can't see any missing parameters from these equations.
c I do also note the "at" is never used, as well as not being initialised
c so perhaps it is just a redundant variable. AJL

         dFdW=   (      at1 + Wsig*(Two  *at2 + Wsig*(Three*at3 + Wsig*(
     &            Four *at4 + Wsig*(Five *at5 + Wsig*(Six  *at6 + Wsig*(
     &            Seven*at7 + Wsig*(Eight*at8 + Wsig*(Nine *at9 + Wsig*(
     &            F10  *at10+ Wsig*F11*at11))))))))))
         dWdT = Two/((One + Tsig)**2)
         dTdR = ((Six*PI*PI)**F2o3)*(rhoo**F2o3)/tauu
         dTdTau = -TauUEG/tauu**2
         dGGAdR = F4o3*rho13*(Ax+(E-x*dE))
         dFdR = dFdW*dWdT*dTdR
         dFdTau=dFdW*dWdT*dTdTau
         dGGAdG =(dE/(Two*Gamma))
c
*        dF/dRhob
c           D1F(i,dRB) = D1F(i,dRB) + dGGAdR*Fsig+ (Ax+E)*rho43*dFdR
           D1F(dRB) = D1F(dRB) + dGGAdR*Fsig+ (Ax+E)*rho43*dFdR
*        dF/dGammaaa
c           D1F(i,dGB) = D1F(i,dGB) + dGGAdG*Fsig
           D1F(dGB) = D1F(dGB) + dGGAdG*Fsig
*        dF/dTaua
c            D1F(i,dTB) = D1F(i,dTB) + rho43*(Ax+E)*dFdTau
            D1F(dTB) = D1F(dTB) + rho43*(Ax+E)*dFdTau

       ENDIF
c      enddo
      Return
      End
