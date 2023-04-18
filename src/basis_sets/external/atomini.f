c Header:
c***********************************************************************
c
c spherical DFT atom
c initialize effective potential and eigenvalues
c here: tfapot() == Thomas-Fermi potential
c
c original version by D.R. Hamann, gncpp
c***********************************************************************
      subroutine atomini(imx,mmax,z,n,f,r,vi,e)
c
      implicit real*8 (a-h,o-z)
      integer imx, mmax, n
c
      dimension e(imx),f(imx),n(imx),r(mmax),vi(mmax)
      external tfapot

      integer i
c
      do 100 i=1,mmax
        vi(i) = tfapot(r(i),z) 
  100 continue
c
      sf = 0.d0
      do 200 i=1,imx
        sf = sf+f(i)
        zz = z-sf+1.d0
        e(i) = -.5d0*(zz/n(i))**2
        if(e(i) .gt. vi(mmax)) e(i)=2.0d0*vi(mmax)
  200 continue
c
      return
      end
c
c***********************************************************************
c generalized thomas fermi atomic potential
c
c...to an article of N. H. March ( "The Thomas-Fermi Approximation in
c Quantum Mechanics", Adv. Phys. 6, 1 (1957)). He has the formula,
c but it is not the result of his work. The original publication is:
c     R. Latter, Phys. Rev. 99, 510 (1955).
c He says that it's an analytic fit to an improved calculation of the
c potential distribution of a Thomas-Fermi atom without exchange first
c performed by Miranda (C. Miranda, Mem. Acc. Italia 5, 285 (1934)).
c                                 Alexander Seidl, TU Munich
c
c original version by D.R. Hamann, gncpp
c***********************************************************************
c
      double precision function tfapot(r,z)
c
      double precision b, dsqrt, r, t, x, xs, z
c
      b=(0.69395656d0/z)**.33333333d0
      x=r/b
      xs=dsqrt(x)
c
      t=z/(1.0d0+xs*(0.02747d0 - x*(0.1486d0 - 0.007298d0*x))
     &   + x*(1.243d0 + x*(0.2302d0 + 0.006944d0*x)))
c
      if(t .lt. 1.0d0) t=1.0d0
      tfapot=-t/r
      return
      end
ce
