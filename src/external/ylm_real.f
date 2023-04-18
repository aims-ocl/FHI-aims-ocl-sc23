C  VB, 23.11.04:
C  Subroutine ylm obtained from Ricardo Gomez-Abral, the original source is not
C  known to me. 
C  Given the absence of a formal copyright notice, I assume that the code is in the
C  public domain; however, I do not understand the signature "Institut fuer theoretische
C  Chemie - TU Vienna". 
C
C  VB, 07.02.05:
C  created this modified version which obtains the real Ylm only, i.e. does not
C  bake them together to a complex version as does the original ylm.f
C
C  JW, 21.06.10:
C  * Updated comments to reflect the actual choice of linear combinations.
C  * Add f2py directives.

      SUBROUTINE YLM_REAL(V,LMAX,Y)
      IMPLICIT NONE
Cf2py integer, intent(IN) :: LMAX
Cf2py double precision, intent(IN) :: V(3)
Cf2py real*8, intent(OUT) :: Y((LMAX+1)*(LMAX+1))
      INTEGER            LMAX
      DOUBLE PRECISION   V(3)
      real*8         Y((LMAX+1)*(LMAX+1))
!
!     ..................................................................
! 1.     PROGRAM UNIT 'YLM_REAL'
!           Calculates real-valued spherical harmonics
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Real and imaginary parts of
!           the spherical harmonics (Condon and Shortley convention)
!             Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!           for vector V (given in Cartesian coordinates)
!           are calculated. In the Condon Shortley convention the
!           spherical harmonics are defined as
!                             +------+ 
!                        m    |   1     m              im(Phi)
!           Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!                            \| 2(Pi)   l
!                  m      
!           where P (cos(Theta)) is the normalized Associated Legendre
!                  l
!           function. Thus,
!                                          m      *
!                            Y(l,-m) = (-1) Y(l,m)
!
!           This subroutine then calculates
!
!                             /    sqrt(2.) * real( Y(l, m) )   for m > 0
!                            |
!             Y^r(l, m) :=  <                       Y(l, 0)     for m = 0
!                            |
!                             \  - sqrt(2.) * imag( Y(l, m) )   for m < 0
!
! 3.     USAGE
!           DOUBLE PRECISION V(3), Y(5*5)
!           V(1) = ...
!           V(2) = ...
!           V(3) = ...
!           CALL YLM(V,4,Y)
!
!        ARGUMENT-DESCRIPTION
!           V      - DOUBLE PRECISION vector, dimension 3        (input)
!                    Must be given in Cartesian coordinates.
!                    Conversion of V to polar coordinates gives the
!                    angles Theta and Phi necessary for the calculation
!                    of the spherical harmonics.
!           LMAX   - INTEGER value                               (input)
!                    upper bound of L for which spherical harmonics
!                    will be calculated
!                    constraint:
!                       LMAX .GE. 0 (not checked)
!           Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!                    contains the calculated spherical harmonics
!                    Y(1)                   for L .EQ. 0 (M = 0)
!                    Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!                    ...
!                    Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!                                           for L .EQ. LMAX
!                                              (M = -L,...,L)
!                    constraint:
!                       Dimension of Y .GE. (LMAX+1)**2 (not checked)
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           Type COMPLEX*16 is used which does not conform to the
!           FORTRAN 77 standard.
!           Also the non-standard type conversion function DCMPLX()
!           is used which combines two double precision values into
!           one double complex value.
!
! 4.     REMARKS
!           none
!
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     ..................................................................
!
!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
cvb
      INTEGER            I22L, I24L2
cvb end
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13, PI
      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!
      PI = (4.0D+0)*ATAN(1.0D+0)
!
!        Y(0,0)
!
      YLLR = 1.0D+0/SQRT(4.0D+0*PI)
      YLLI = 0.0D+0
      Y(1) = YLLR
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      IF (LMAX .LE. 0) GOTO 999
!
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0D+0) THEN
         A = V(1)/ABMAX
         B = V(2)/ABMAX
         AB = SQRT(A*A+B*B)
         COSPH = A/AB
         SINPH = B/AB
      ELSE
         COSPH = 1.0D+0
         SINPH = 0.0D+0
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         AB = A*A + B*B
         ABC = SQRT(AB + C*C)
         COSTH = C/ABC
         SINTH = SQRT(AB)/ABC
      ELSE
         COSTH = 1.0D+0
         SINTH = 0.0D+0
      ENDIF
!
!        Y(1,0)
!
      Y(3) = SQRT(3.0D+0)*YLLR*COSTH
!
!        Y(1,1) ( = -DCONJG(Y(1,-1)))
!
      TEMP1 = -SQRT(3.0D+0)*YLLR*SINTH
      Y(4) = TEMP1*COSPH
      Y(2) = -TEMP1*SINPH
!
      DO 20 L = 2, LMAX
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         YL1L1R = Y(INDEX-1)
         YL1L1I = - MSIGN * Y(INDEX-2*L+1)
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
         YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
         YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
         Y(INDEX2) = YLLR
         Y(INDEX)  = MSIGN * YLLI
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = SQRT(DBLE(2*L+1))*COSTH
         YLL1R = TEMP2*YL1L1R
         YLL1I = TEMP2*YL1L1I
         Y(INDEX2) = YLL1R
         Y(INDEX)  = - MSIGN * YLL1I
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX - 4*L + 2
         I2L  = INDEX - 2*L
         I24L2 = INDEX2 - 4*L + 2
         I22L  = INDEX2 - 2*L
         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO 10 M = L - 2, 0, -1
!
!        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            YLMR = TEMP2*Y(I22L) + TEMP3*Y(I24L2)
            YLMI = TEMP2*Y(I2L) + TEMP3*Y(I4L2)
            Y(INDEX2) = YLMR
            Y(INDEX)  = YLMI
!
            INDEX2 = INDEX2 - 1
            INDEX  = INDEX  + 1
            I24L2   = I24L2   - 1
            I22L    = I22L    - 1
            I4L2   = I4L2   + 1
            I2L    = I2L    + 1
   10    CONTINUE
   20 CONTINUE

!
  999 RETURN
!
!        End of 'YLM'
!
      END
