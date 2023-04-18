      SUBROUTINE CRPAPW(DUP,DDN,VCUP,VCDN,ECPERE)
!engel
!
!  $Id: crpapw.f,v 1.1 1996/10/23 11:28:02 engel Stab $
!
!  $Log: crpapw.f,v $
!  Revision 1.1  1996/10/23 11:28:02  engel
!  Initial revision
!
!
!  This subroutine calculates the spin-up and spin-down correlation
!  potentials (VCUP,VCDN --- in Hartree) as well as the correlation
!  energy per electron (ECPERE --- in Hartree) within the RPA/LSDA (for
!  given spin-densities DUP and DDN), using the form of Perdew-Wang
!  92. It is a slightly modified version of Perdew and Wang's original
!  LSDA correlation energy subroutine.
!
!  Input:  DUP,DDN = spin-up,spin-down densities (in Bohr**(-3))
!  Output: VCUP,VCDN = spin-up,dn correlation potentials (in Hartree)
!          ECPERE = correlation energy per electron (in Hartree)
!  Refs:          J.P.Perdew, Y.Wang, PRB45, 13244 (1992)
!  Written by:    J.P.Perdew, Y.Wang, E.Engel
!  Last revision: October 23, 1996
!engel
      IMPLICIT NONE
      REAL*8 DUP,DDN,VCUP,VCDN,ECPERE,PI,FPI,ONETRD,FORTRD,GAM,FZZ, &
             RS,ZETA,F,FZ,EU,EURS,EP,EPRS,ALFM,ALFRSM,ALFC,Z4,ECRS, &
             ECZETA,COMM
      PARAMETER( PI     = 3.141592653589793D0, &
                 FPI    = 4.D0*PI, &
                 ONETRD = 1.D0/3.D0, &
                 FORTRD = 4.D0/3.D0, &
                 GAM    = 0.5198421D0, &
                 FZZ    = 1.709921D0 )
      EXTERNAL GCOR

      IF( (DUP+DDN).LT.1.D-30 ) THEN
         RS = ( 3.D0 / (FPI*1.D-30) )**0.333333333333333333D0
         ZETA = ( DUP - DDN ) / 1.D-30
      ELSE
         RS = ( 3.D0 / (FPI*(DUP+DDN)) )**0.333333333333333333D0
         ZETA = ( DUP - DDN ) / ( DUP + DDN )
      END IF
      IF( ZETA.GE.1.D0 .OR. ZETA.LE.-1.D0 ) ZETA = SIGN( 1.D0, ZETA )

      F = ( (1.D0+ZETA)**FORTRD + (1.D0-ZETA)**FORTRD - 2.D0 ) / GAM
      FZ = FORTRD * ((1.D0+ZETA)**ONETRD-(1.D0-ZETA)**ONETRD) / GAM

      CALL GCOR(0.0310907D0,0.082477D0,5.1486D0,1.6483D0,0.23647D0, &
          0.20614D0,0.75D0,RS,EU,EURS)
      CALL GCOR(0.01554535D0,0.035374D0,6.4869D0,1.3083D0,0.15180D0, &
          0.082349D0,0.75D0,RS,EP,EPRS)
      CALL GCOR(0.0168869D0,0.028829D0,10.357D0,3.6231D0,0.47990D0, &
          0.12279D0,0.75D0,RS,ALFM,ALFRSM)
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZETA**4
      ECPERE = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      ECZETA = 4.D0*(ZETA**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = ECPERE -RS*ECRS/3.D0-ZETA*ECZETA
      VCUP = COMM + ECZETA
      VCDN = COMM - ECZETA
      RETURN
      END
!C--------------------------------------------------------------------------------
!      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
!Cengel
!C
!C  $Id: gcor.f,v 1.1 1996/10/23 11:30:48 engel Stab $
!C
!C  $Log: gcor.f,v $
!C  Revision 1.1  1996/10/23 11:30:48  engel
!C  Initial revision
!C
!C
!C  This subroutine evaluates the basic function which is used by
!C  Perdew and Wang as parametrization of the various components of
!C  their LSDA.
!C
!C  Refs:          J.P.Perdew, Y.Wang, PRB45, 13244 (1992) (-->Eq.10)
!C  Written by:    J.P.Perdew, Y.Wang, E.Engel
!C  Last revision: October 23, 1996
!Cengel
!      IMPLICIT REAL*8 (A-H,O-Z)
!      P1 = P + 1.D0
!      Q0 = -2.D0*A*(1.D0+A1*RS)
!      RS12 = SQRT(RS)
!      RS32 = RS12**3
!      RSP = RS**P
!      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
!      Q2 = LOG(1.D0+1.D0/Q1)
!      GG = Q0*Q2
!      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
!      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
!      RETURN
!      END
