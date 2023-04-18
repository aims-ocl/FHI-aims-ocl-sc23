c Header:
c**********************************************************************c
c  Becke exchange for a spin-unpolarized electronic system 
c
c  Gradient-corrected exchange energy based on
c     [A.D. Becke, J.Chem.Phys.96, 2155, 1992].
c  The LSDA energy functional, obtained as E{n+,n-}=(E{2n+}+E{2n-})/2,
c     and the functional derivative formula are given by 
c     [J.P. Perdew , PRB 33, 8800, 1986].
c     [J.P. Perdew , PRB 34, 7406, 1986].
c  see also [G.Ortiz ...,PRB 43, 6376 (1991)] eq. (A2)
c  
c  Hartree a.u.
c
c  constants    b is beta
c               c = 2*(6*pi*pi)^(1/3) = c1
c               ax = 3/4 (3/pi)^(1/3) = Cx           
c  input
c  d            density
c  s            abs(grad d)/(2kf*d)
c  fs           1 / s df / ds
c  output
c  en_density_x           exchange energy per electron
c  x_density_deriv        energy derivative with respect to density
c  x_gradient_deriv       energy derivative with respect to gradient
c**********************************************************************
c
      subroutine xbecke_derivs (rho, squared_grad_rho,
     +  en_density_x, x_density_deriv, x_gradient_deriv )
c
      implicit none

C     declare constants

      real*8 :: thrd, thrd2, thrd4, thrd8, pi, ax, c, b, bb
      real*8 :: x, y1, y0, ddi, dd1, g, s
      parameter(thrd=1.d0/3.d0,thrd2=2.d0/3.d0,thrd4=4.d0/3.d0,
     +          thrd8=8.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(c = .779555417944150792d1)
      parameter(b = .42d-2)
      parameter(bb = -.451357747124625192d-2)

c     declare variables

      real*8 :: rho
      real*8 :: squared_grad_rho
      real*8 :: en_density_x
      real*8 :: x_density_deriv
      real*8 :: x_gradient_deriv

      real*8 :: exunif
      real*8 :: scale
      real*8 :: SSQ
      real*8 :: f
      real*8 :: fs

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density

      exunif = ax*rho**thrd

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      scale = (3.d0*pi*pi)**thrd2
      scale = 4.d0*scale*(rho**thrd8)
      SSQ = squared_grad_rho / scale

      s = dsqrt(SSQ)
c      write(6,*) "s in own:", s
      x  = c*s
      y1 = 1.d0/dsqrt(1.d0+x*x)
      y0 = dlog(x+1.d0/y1)
c      y2 = -x*y1*y1*y1
      ddi= 1.d0/(1.d0 + 6.d0*b*x*y0)
      dd1= 6.d0*b*(y0+x*y1)
      g  = 1.d0 - 0.5d0*x*dd1*ddi
      fs = -2.d0*bb*c*c*ddi
c      g1 = -3.d0*b*(y0+x*(3.d0*y1+x*y2-dd1*dd1*ddi/(6.d0*b)))
c      fss= fs*c*(g1 - g*dd1)*ddi
      fs = fs*g


      f  = 1.d0 - bb*x*x*ddi
     
      en_density_x = exunif*f
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE partial derivatives of rho*exunif(rho)*FxPBE(ssq) :
c  first, find derivative of Fx w.r.t ssq:
c  Fs = d f/ d(ssq)
C  next, multiply with derivatives d(ssq)/d(rho) and d(ssq)/d(squared_grad_rho)
C  to obtain desired partial derivatives;
C  the definition of ssq gives
C  d(ssq)/d(rho) = - 8/3 ssq/rho
C  d(ssq)/d(squared_grad_rho) = 1/ (2*KF*rho)^2 = 1/scale

      x_density_deriv = thrd4*en_density_x - exunif*fs*thrd4*ssq

      x_gradient_deriv = rho*exunif*fs/ (2*scale)


      return
      end

