c Header:
cc**********************************************************************
c  LDA exchange
c  Hartree a.u.
c***********************************************************************
      subroutine xlda(d,vx,ex)
c
      implicit real *8 (a-h,o-z)
      data ax,thd,thd4,eps/-.738558766382022406d0
     &,                     .333333333333333333d0
     &,                     .133333333333333333d1
     &,                    1.d-100/
c
      if(d .le. eps) then
        vx=0.d0
        ex=0.d0
      else
        ex=ax*d**thd
        vx=thd4*ex
      endif
c
      return
      end
ce
c***********************************************************************
c Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).
c Hartree a.u.
c see W.E. Pickett, Comp.Phys.Rep. 9, 115 (1989), pg. 187 : w2 = 7.79
c***********************************************************************
c
      subroutine wigner(rh,ex,fx,exc,fxc)
      implicit real*8 (a-h,o-z)
      data thrd,w1,w2,ax,pif,eps/.333333333333333333d0
     &,                         -.44d0
     &,                          .779d1
     &,                         -.738558766382022447d0
     &,                          .620350490899400087d0
     &,                         1.d-100/

c
      if(rh .lt. eps) then
        ex = 0.d0
        fx = 0.d0
        exc= 0.d0
        fxc= 0.d0
      else
c
        x  = rh**thrd
        ex = ax * x
        fx = 4.d0*thrd*ex
        rs = pif/x
        y  = 1.d0/(rs + w2)
        exc= ex + w1*y
        fxc= fx + w1*(4.d0*thrd*rs + w2)*y*y
      endif
c
      return
      end
ce
c**********************************************************************
c LDA - Ceperley - Alder exchange-correlation potential and energy
c as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
c Hartree a.u.
c original version by D.R. Hamann, gncpp
c**********************************************************************
c
      subroutine cepal(rh,ex,fx,exc,fxc)
      implicit real *8 (a-h,o-z)
      data eps/1.d-100/
c
      if(rh .lt. 0.23873241d0) then
c
        if(rh .le. eps) then
          fxc=0.d0
          exc=0.d0
          ex=0.d0
          fx=0.d0
          goto 100
        endif
c
        rs=0.62035049d0*rh**(-0.3333333333333333d0)
        sqrs=dsqrt(rs)
        den=1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
        exc=-0.4582d0/rs - 0.1423d0/den
        fxc=exc - rs*(0.15273333d0/rs**2
     &    + (0.02497128d0/sqrs + 0.01581427d0)/den**2)
      else
        rs=0.62035049d0*rh**(-0.3333333333333333d0)
        rsl=dlog(rs)
        exc=-0.4582d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs
     &    + 0.002d0*rs*rsl
        fxc=exc - rs*(0.15273333d0/rs**2
     &    + 0.01036667d0/rs - 0.003866667d0
     &    + 0.00066667d0*(1.0d0 + rsl))
      end if
c
c exchange-only energy and potential
      ex=-0.7385587664d0*rh**(.3333333333333333d0) 
      fx=4.d0/3.d0*ex
c
  100 return
      end
c
cc**********************************************************************
c LDA - scaled Wigner exchange-correlation
c from [Q. Zhao, R.G. Parr, PRA 46, R5320 (1992)] eqs. (5) and (7)
c
c Input   
c rh        density
c
c Output
c ex        exchange energy per electron
c fx           ""    potential
c           ... based on continuum-LDA (electron-gas like)
c
c exc       exchange-correlation energy per electron
c fxc          ""        ""      potential
c
c Hartree a.u.
c Martin Fuchs, FHI der MPG, Berlin, 01-1993
cc**********************************************************************
      subroutine wigscaled(rh,ex,fx,exc,fxc)
      implicit real *8 (a-h,o-z)
      intrinsic log
      data thrd,aa,bb,ax,eps/.333333333333333333d0
     &,                      .93222d0
     &,                      .947362d-2
     &,                      .738558766382022406d0
     &,                     1.d-100/
c
      if(rh .le. 0.d0) then
        rh = 0.d0
        ex = 0.d0
        fx = 0.d0
        exc= 0.d0
        fxc= 0.d0
      else
c
        x = bb*rh**thrd
        y = x/(x + 1.d0)
        z = -aa*x/bb
        xyln = x*log(y)
c
        exc = z * (1.d0 + xyln)
        fxc = thrd*z *(4.d0 + 5.d0*xyln + y)
c
c electron-gas based LDA exchange
c
        ex = ax*z/aa
        fx = 4.d0*thrd*ex
      endif
c      
      return
      end
c
