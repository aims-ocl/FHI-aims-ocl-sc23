!---------------------------------------------------------------------------------
!  Subroutine increment_ylm_deriv calculates recursively the derivatives
!  d(Y)/d(theta) and 1/sin(theta)*d(Y)/d(phi), which are needed to compute
!  the [cartesian] gradient of the wave functions.
!
!  I know this could be done directly from the Y_lm, without recourse to an explicit
!  recursion formulae for the derivatives, BUT - whenever we hit the coordinate singularity
!  theta=0 (z axis), the 1/sin(theta) division kills the computation of gradients.
!
!  My best solution is to cancel out the singularity directly in the calculation of the
!  derivative.
!
!  See the original subroutine increment_ylm for a detailed explanation of the recursion
!  formalism, and also for a notice regarding the (unclear, but hopefully safe) copyright
!  of the original routine on which the code below is based.


      SUBROUTINE increment_ylm_deriv &
      ( SINTH, COSTH, SINPH, COSPH, LMIN, LMAX, Y, YTH, &
        YPH)

      IMPLICIT NONE

!     input

      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      INTEGER            LMIN
      INTEGER            LMAX

!     output
!     Y() is recomputed as well, even if it was calculated already earlier.
!     This does not cost anything as we need to follow through the entire recursion
!     anyhow.
!     YTH is dy/d(theta)
!     YPH is 1/sin(theta)*[dy/d(phi)]

      real*8             Y((LMAX+1)*(LMAX+1))
      real*8             YTH((LMAX+1)*(LMAX+1))
      real*8             YPH((LMAX+1)*(LMAX+1))

!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
!vb
      INTEGER            I22L, I24L2
!vb end
      DOUBLE PRECISION   D4LL1C, D4LL1S, D2L13, PI
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3, TEMP4, TEMP5
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YTHLLR, YTHLL1R, YTHL1L1R, YTHLMR
      DOUBLE PRECISION   YPHLLR, YPHLL1R, YPHL1L1R, YPHLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
      DOUBLE PRECISION   YTHLLI, YTHLL1I, YTHL1L1I, YTHLMI
      DOUBLE PRECISION   YPHLLI, YPHLL1I, YPHL1L1I, YPHLMI
!
      PI = (4.0D+0)*ATAN(1.0D+0)

      if (lmin.le.0) then
!
!        Y(0,0)
!
         YLLR = 1.0D+0/SQRT(4.0D+0*PI)
         YLLI = 0.0D+0
         Y(1) = YLLR

         YTH(1) = 0.d0
         YPH(1) = 0.d0
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      end if

      if ( (lmin.le.1).and.(lmax.ge.1)) then
!
!       Y(1,0)
!
        Y(3) = SQRT(3.0D+0)*YLLR*COSTH

        YTH(3) = - SQRT(3.0D+0)*YLLR*SINTH
        YPH(3) = 0.d0

!
!       Y(1,1) ( = -DCONJG(Y(1,-1)))
!
        TEMP1 = SQRT(3.0D+0)*YLLR

        TEMP2 = -TEMP1*SINTH
        Y(4) = TEMP2*COSPH
        Y(2) = -TEMP2*SINPH

        TEMP2 = -TEMP1*COSTH
        YTH(4) = TEMP2*COSPH
        YTH(2) = -TEMP2*SINPH

        YPH(4) = TEMP1*SINPH
        YPH(2) = TEMP1*COSPH

      end if

!     Now calculate remaining ylm's, derivatives in table

      DO 20 L = max( 2, lmin), LMAX, 1
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         YL1L1R = Y(INDEX-1)
         YL1L1I = - MSIGN * Y(INDEX-2*L+1)
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))
         TEMP2 = TEMP1*SINTH

         YLLR = (COSPH*YL1L1R - SINPH*YL1L1I)
         YLLI = (COSPH*YL1L1I + SINPH*YL1L1R)

         Y(INDEX2) = TEMP2 * YLLR
         Y(INDEX)  = TEMP2 * MSIGN * YLLI

!        YTHLL = YTH(L,L) = f(Y(L-1,L-1),YTH(L-1,L-1)) ... Formula 1
!
         YTHL1L1R = YTH(INDEX-1)
         YTHL1L1I = - MSIGN * YTH(INDEX-2*L+1)
         TEMP3 = TEMP1*COSTH

         YTHLLR = (COSPH*YTHL1L1R - SINPH*YTHL1L1I)
         YTHLLI = (COSPH*YTHL1L1I + SINPH*YTHL1L1R)

         YTH(INDEX2) = TEMP3 * YLLR + TEMP2 * YTHLLR
         YTH(INDEX)  = TEMP3 * MSIGN * YLLI + TEMP2 * MSIGN * YTHLLI

!        YPHLL = YPH(L,L) = f(Y(L-1,L-1),YPH(L-1,L-1)) ... Formula 1
!
         YPHL1L1R = YPH(INDEX-1)
         YPHL1L1I = - MSIGN * YPH(INDEX-2*L+1)

         YLLR = (-SINPH*YL1L1R - COSPH*YL1L1I)
         YLLI = (-SINPH*YL1L1I + COSPH*YL1L1R)

         YPHLLR = (COSPH*YPHL1L1R - SINPH*YPHL1L1I)
         YPHLLI = (COSPH*YPHL1L1I + SINPH*YPHL1L1R)

         YPH(INDEX2) = TEMP1 * YLLR + TEMP2 * YPHLLR
         YPH(INDEX)  = TEMP1 * MSIGN * YLLI + TEMP2 * MSIGN * YPHLLI
!
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP1 = SQRT(DBLE(2*L+1))
         TEMP2 = TEMP1*COSTH
         YLL1R = TEMP2*YL1L1R
         YLL1I = TEMP2*YL1L1I
         Y(INDEX2) = YLL1R
         Y(INDEX)  = - MSIGN * YLL1I
!
!        YTHLL1 = YTH(L,L-1) = f(Y(L-1,L-1),YTH(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1),YTH(L-2,L-1) in Formula 2 is zero)
!
         TEMP3 = -TEMP1*SINTH
         YLL1R = TEMP3*YL1L1R
         YLL1I = TEMP3*YL1L1I
         YTHLL1R = TEMP2*YTHL1L1R
         YTHLL1I = TEMP2*YTHL1L1I
         YTH(INDEX2) = YLL1R + YTHLL1R
         YTH(INDEX)  = - MSIGN * ( YLL1I + YTHLL1I )

!
!        YPHLL1 = YPH(L,L-1) = f(YPH(L-1,L-1)) ... Formula 2
!               (the coefficient for YPH(L-2,L-1) in Formula 2 is zero)
!
!         TEMP2 = TEMP2*SINTH
         YPHLL1R = TEMP2*YPHL1L1R
         YPHLL1I = TEMP2*YPHL1L1I
         YPH(INDEX2) = YPHLL1R
         YPH(INDEX)  = - MSIGN * YPHLL1I

         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX - 4*L + 2
         I2L  = INDEX - 2*L
         I24L2 = INDEX2 - 4*L + 2
         I22L  = INDEX2 - 2*L
         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
         D4LL1S = SINTH*SQRT(DBLE(4*L*L-1))
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
            TEMP4 = - D4LL1S*TEMP1
            YTHLMR = TEMP4*Y(I22L) + TEMP2*YTH(I22L) + TEMP3*YTH(I24L2)
            YTHLMI = TEMP4*Y(I2L) + TEMP2*YTH(I2L) + TEMP3*YTH(I4L2)
            YTH(INDEX2) = YTHLMR
            YTH(INDEX)  = YTHLMI
!
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            YPHLMR = TEMP2*YPH(I22L) + TEMP3*YPH(I24L2)
            YPHLMI = TEMP2*YPH(I2L) + TEMP3*YPH(I4L2)
            YPH(INDEX2) = YPHLMR
            YPH(INDEX)  = YPHLMI
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
!        End of 'YLM'
!
      END
