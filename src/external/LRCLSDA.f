      SUBROUTINE LRCLSDA(Emu,Rho,F,D1F)
c
c***********************************************
c                                               
c   INPUT:                                      
c      Emu - Value of mu (or omega)
c      Rho - Spin density                 
c                                               
c   OUTPUT:                                     
c      F      - Functional value               
c      D1F    - First derivative               
c                                               
c***********************************************
c
      IMPLICIT REAL*8 (a-h,o-z)
      Save F1, F2, F3, F4, F5, F6, F7, F8, F9
      DATA F1/1.0D+00/,F2/2.0D+00/,F3/3.0D+00/,F4/4.0D+00/,F5/5.0D+00/,
     $     F6/6.0D+00/,F7/7.0D+00/,F8/8.0D+00/,F9/9.0D+00/
C
      PARAMETER( PI = 3.1415926535897932384626433832795D+00 )
C
      F1o2 = F1 / F2
      F1o3 = F1 / F3
      F1o4 = F1 / F4
      F2o3 = F2 / F3
      F4o3 = F4 / F3
      F4o9 = F4 / F9
      F8o3 = F8 / F3
      PI12 = SQRT(Pi)
C
      AX   = -(F3/F2) * (F4o3*PI)**(-F1o3)
      Cmu  = (F6*Pi**F2)**F1o3   
C
      Rho13 = Rho**F1o3
      Rho43 = Rho**F4o3
c
      tmu  = Emu/(F2*Cmu*Rho13)
      tmu2 = tmu*tmu
      tmu3 = tmu*tmu2
      tmu4 = tmu*tmu3
c
      W    = Exp(-F1o4/tmu2)
      ERFV = Erf( F1o2/tmu)
      dtmudR = -F1o3*tmu / Rho
c
      Fsr = F1-F4o3*tmu*(-F6*tmu+F8*tmu3+W*
     $        (F4*tmu-F8*tmu3)+F2*PI12*ERFV)
      dFsrdtmu = F8o3*(F2*tmu*(F3-F8*tmu2+W*
     $          (-F1+F8*tmu2))-PI12*ERFV)
c
      F = Ax*Rho43*Fsr
      D1F = Ax*F4o3*Rho13*Fsr + Ax*Rho43*(dFsrdtmu*dtmudR)
C
      RETURN
      END


