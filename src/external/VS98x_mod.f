      SUBROUTINE  VS98x(Tol,fac,F,D1F,RA,RB,D1RA,D1RB,TA,TB,ijzy)
c     ******************************************************************
c     *                                                                *
c     *  VS98x evaluates the exchange part of VS98 and its             *
c     *  variants on a grid.                                           *
c     *  !!! Second derivatives are not available yet.                 *
c     *                                                                *
c     *  Ref:  T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109,   *
c     *        400 (1998).                                             *
c     *       Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125,           * 
c     *        194101 (2006).                                          *
c     *                                                                *
c     *       ijzy - 1 exchange functional in VS98                     *
c     *       ijzy - 2 exchange functional in M06-L                    *
c     *       ijzy - 3 exchange functional in M06-HF                   *
c     *       ijzy - 4 exchange functional in M06                      *
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
c*  YZ (12/08)                                                          *
c*                                                                      *
c************************************************************************

      IMPLICIT NONE

      REAL*8  F,D1F(7),RA,RB,
     $        D1RA(3),D1RB(3),TA,TB
      REAL*8 pi,Tol, fac

c
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC
      integer i,j, ijzy
      REAL*8 rrho, rho43, rho13, rhoo, rho53, rho83
      REAL*8  Gamma
c
c     kinetic energy density or tau
c

      REAL*8 tauN,tauu,DTol
      REAL*8 Tiny, f13, f43, f53, f83, f113 
      REAL*8 gx, gg, x, z, kx,xk,zk
      REAL*8 One, Two, Three, Four, Five, Six, Seven, Eight
      REAL*8 Nine, F10, F11
      REAL*8 cf, Axlsda, r1, r2, r3, r4, r5, r6

c      functional derivatives below FFFFFFFFFFFF

       REAL*8 dxdr, dxdg, dzdr, dzdt, dgdx, dgdz

c      functional derivatives above FFFFFFFFFFFF

       parameter( pi = 3.1415926535897932384626433832795d0 )
         
       parameter (cf = 9.115599720d0, Axlsda = -0.9305257363491d0 )
       parameter (gg  = 0.00186726d0) 
       parameter (f13=1.d0/3.d0,f43=4.0d0/3.0d0,f53=5.0d0/3.0d0)
       parameter (f83=8.d0/3.0d0, F113=11.0d0/3.d0)
       parameter (One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, 
     &             Five=5.0d0,Six=6.0d0, Seven=7.0d0,
     &             Eight=8.0d0, Nine=9.0d0,F10=10.d0, F11=11.d0)

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

      if (ijzy.eq.1) then
c
c     Parameters for VS98
c
        r1=  -9.800683d-01
        r2=  -3.556788d-03
        r3=   6.250326d-03
        r4=  -2.354518d-05
        r5=  -1.282732d-04
        r6=   3.574822d-04
      elseif (ijzy.eq.2) then
c
c     Parameters for M06-L
c
        r1 =   6.012244D-01*Axlsda
        r2 =   4.748822D-03*Axlsda
        r3 =  -8.635108D-03*Axlsda
        r4 =  -9.308062D-06*Axlsda
        r5 =   4.482811D-05*Axlsda
        r6 =   0.000000D+00
      elseif (ijzy.eq.3) then
c
c     Parameters for M06-HF
c
        r1 =   -1.179732D-01*Axlsda
        r2 =   -2.500000D-03*Axlsda
        r3 =   -1.180065D-02*Axlsda
        r4 =   0.000000D+00
        r5 =   0.000000D+00
        r6 =   0.000000D+00
      elseif (ijzy.eq.4) then
c
c     Parameters for M06
c
        r1 =   1.422057D-01*Axlsda
        r2 =   7.370319D-04*Axlsda
        r3 =   -1.601373D-02*Axlsda
        r4 =   0.000000D+00
        r5 =   0.000000D+00
        r6 =   0.000000D+00
      endif
      
      DTol = Tol 

c      DO i = 1,NGrid
         IF (RA.gt.DTol.and.TA.gt.DTol) THEN
            rhoo = RA
            rho43 = rhoo**F43  
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53
            rho83 = rho53*rhoo
c
            tauu = TA
            Gamma =  D1RA(1)**2 + D1RA(2)**2 + D1RA(3)**2
            x = gamma/rho83
            dxdr = -f83*x*rrho
            dxdg = One/rho83
            z = tauu/rho53 - cf
            dzdr = -f53 * tauu/rho83
            dzdt = One/rho53
            kx = One + gg*x + gg*z
            xk = x/kx
            zk = z/kx
            call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,r1,r2,r3,r4,r5,r6)
            
            F = F + rho43*gx
c
c     functional derivatives 
c
            D1F(dRA) = D1F(dRA) + f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)    
            D1F(dGA) = D1F(dGA) + rho43*(dgdx*dxdg)
            D1F(dTA) = D1F(dTA) + rho43*(dgdz*dzdt)
         ENDIF
         IF (RB.gt.DTol.and.TB.gt.DTol) THEN
            rhoo = RB
            rho43 = rhoo**F43
            rrho = 1d0/rhoo       ! reciprocal of rho
            rho13 = rho43*rrho
            rho53 = rhoo**F53
            rho83 = rho53*rhoo
c
            tauu = TB
            Gamma =  D1RB(1)**2 + D1RB(2)**2 + D1RB(3)**2
            x = gamma/rho83
            dxdr = -f83*x*rrho
            dxdg = One/rho83
            z = tauu/rho53 - cf
            dzdr = -f53 * tauu/rho83
            dzdt = One/rho53
            kx = One + gg*x + gg*z
            xk = x/kx
            zk = z/kx
            call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,r1,r2,r3,r4,r5,r6)
            F = F + rho43*gx
c
c     functional derivatives
c
            D1F(dRB) = D1F(dRB) + f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
            D1F(dGB) = D1F(dGB) + rho43*(dgdx*dxdg)
            D1F(dTB) = D1F(dTB) + rho43*(dgdz*dzdt)
         ENDIF
c       ENDDO
 
       RETURN
       END

      Subroutine gvt4(gvt,dg_dx,dg_dz,xg,zg,gama,ct,a,b,c,d,e,f)
      Implicit none
c
c     Evaluate the GVT4 form in VS98 
c
c    some working variables
      REAL*8 gvt,dg_dx,dg_dz,xg,zg,gama,ct,a,b,c,d,e,f
      REAL*8 F1,F2,F3,g,g2,ct2 
      Data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/
C

      g=gama
      g2=gama*gama
      
      gvt =(a + b*xg + c*zg + d*xg*xg + e*zg*xg + f*zg*zg)/g
      dg_dx =(-a*ct+b*(F1-F2*ct*xg)-F2*c*zg*ct+d*(F2*xg-F3*xg*xg*ct)
     $  +e*(zg -F3*zg*xg*ct)-F3*f*zg*zg*ct )/g2
      dg_dz =(-a*ct -F2*b*xg*ct +c*(F1-F2*zg*ct)-F3*d*xg*xg*ct
     $  +e*(xg-F3*xg*zg*ct)+f*(F2*zg-F3*zg*zg*ct))/g2

      return
      end


