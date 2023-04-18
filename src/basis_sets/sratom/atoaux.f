c Header:
c***********************************************************************
c auxiliary subroutines for atom calculation
c***********************************************************************
c adams extrapolation and interpolation formulas for
c outward and inward integration, abramowitz and
c stegun, p. 896
      real*8 function aeo(y,j,mx)
c
      real*8 y
      integer j
      integer mx
      include 'parameter.h'
c
      dimension y(mx)
      aeo=(4.1666666666667d-2)*(55.0d0*y(j)-59.0d0*y(j-1)
     & +37.0d0*y(j-2)-9.0d0*y(j-3))
      return
      end
c
c
      real*8 function aio(y,j,mx)
c
      real*8 y
      integer j
      integer mx
      include 'parameter.h'
c
      dimension y(mx)
      aio=(4.1666666666667d-2)*(9.0d0*y(j+1)+19.0d0*y(j)
     & -5.0d0*y(j-1)+y(j-2))
      return
      end
c
      real*8 function aei(y,j,mx)
c
      real*8 y
      integer j
      integer mx
      include 'parameter.h'
c
      dimension y(mx)
      aei=-(4.1666666666667d-2)*(55.0d0*y(j)-59.0d0*y(j+1)
     & +37.0d0*y(j+2)-9.0d0*y(j+3))
      return
      end
c
      real*8 function aii(y,j,mx)
c
      real*8 y
      integer j
      integer mx
      include 'parameter.h'
c
      dimension y(mx)
      aii=-(4.1666666666667d-2)*(9.0d0*y(j-1)+19.0d0*y(j)
     & -5.0d0*y(j+1)+y(j+2))
      return
      end
c
