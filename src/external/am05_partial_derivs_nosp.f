c Header:
c**********************************************************************
c Armiento Mattsson am05 functional for exchange and correlation 
c Version: 3
c
c input
c   n        electron density [bohr**(-3)]
c   s        scaled abs of gradient of density
c   u        scaled grad n * grad | grad n |
c   t        scaled laplacian of density
c   pot      integer: 1 = calculate potential,
c              0 = don't calculate potential (then u and t 
c              are never touched)
c
c For exact definitions of s, u, and t, see PRB 33, 8800 (1986).
c
c output
c   ex       exchange energy per electron [hartree]
c   ec       correlation energy per electron [hartree]
c   vx       exchange potential [hartree]
c   vc       correlation potential [hartree]
c
c Citation request: when using this functional, please cite:
c "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005). " 
c
c Copyright (c) 2005, Rickard Armiento
c All rights reserved.
c
c Permission is hereby granted, free of charge, to any person obtaining 
c a copy of this software and associated documentation files (the 
c "Software"), to deal in the Software without restriction, including 
c without limitation the rights to use, copy, modify, merge, publish, 
c distribute, sublicense, and/or sell copies of the Software, and to 
c permit persons to whom the Software is furnished to do so, subject 
c to the following conditions:
c
c The above copyright notice and this permission notice shall be 
c included in all copies or substantial portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
c OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
c NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
c HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
c WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
c FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
c OTHER DEALINGS IN THE SOFTWARE.
c
c**********************************************************************

      subroutine am05_partial_derivs_nosp(n,squared_grad,ex,ec,dexdr,decdr,
     *           dexdgrad, decdgrad)

      implicit none

c     ** Input parameters
      real*8 n, squared_grad
c      integer pot

c     ** Output parameters
      real*8 ex, ec, dexdr, decdr, dexdgrad, decdgrad 

c     ** Constants and shorthand notation
      real*8 pi, g, a, s2, c
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (g = 0.8098d0, a = 2.804d0)
      parameter (c = 0.7168d0)

c     ** Local variables
      real*8 exlda, vxlda, eclda, vclda, X, Xsos, Xsossos, s
      real*8 Hx, Hxsos, Hxsossos, Hc, Hcsos, Hcsossos
      real*8 F, Fsos, s2Fsossos
      real*8 szsoz, mixder, scale
      real*8 denom, denomsos, sdenomsoss
      real*8 zfac, zosn, w

c     ** Cutoff 
c      if((n .le. 1.0d-16)) then
c         ex = 0.0d0
c         ec = 0.0d0
c         return
c      endif
c   ******************
c     Aims scale
c   ******************
      scale = (3.d0*pi*pi)**(2.d0/3.d0)
      scale = 4.d0*scale*(n**(8.d0/3.d0))
      s2 =  abs(squared_grad / scale)

      s = sqrt(s2)

c     *******************
c       LDA correlation
c     *******************
      call am05_ldapwc(n,eclda,vclda)

c     *******************
c        LDA exchange
c     *******************
      call am05_ldax(n,exlda,vxlda)
      
c     ********************
c          Exchange
c     ********************

c     ** Interpolation index
      X = 1.0d0/(1.0d0 + a*s2)
      
c     ** Airy LAA refinement function
      call am05_lambertw(s**(3.0d0/2.0d0)/sqrt(24.0d0),w)

c     ** am05_lambertw give back argument if it is < 1.0e-20
c     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
c     ** zosn = normalized z/s
      if (s < 1.0e-14) then
              zosn = 1.0d0
      else
              zosn = 24.0d0**(1.0d0/3.0d0)*w**(2.0d0/3.0d0)/s
      end if
      zfac = s2*(zosn*27.0d0/32.0d0/pi**2)**2

c     ** denom = denominator of Airy LAA refinement function
      denom = 1.0d0 + c*s2*zosn*(1.0d0 + zfac)**(1.0d0/4.0d0)
      F = (c*s2 + 1.0d0)/denom
      
c     ** Exchange refinement function
      Hx = X + (1.0d0 - X)*F

c     ** Exchange energy per particle, Ex = Integrate[n*ex]
      ex = exlda*Hx
      
c     ********************
c         Correlation
c     ********************
      
c     ** Correlation refinement function 
      Hc = X + g*(1.0d0 - X)
      
c     ** Correlation energy per particle, Ec = Integrate[n*ec]
      ec = eclda*Hc
      
c      if (pot .eq. 0) return
      
c     ***************************
c         Exchange potential
c     ***************************

c     ** Interpolation index derivatives: 1/s dX/ds, 1/s d/ds(1/s dX/ds)                   
      Xsos = -2.0d0*a*X**2
      Xsossos = 8.0d0*a**2*X**3

c     ** Airy LAA refinement function derivatives, 1/s dF/ds and s^2 1/s d/ds (1/s dF/ds) 
c     ** szsoz = s*(dz/ds)/z, mixder = szsoz + s^2*(d^2z/ds^2)/z
      szsoz = 1.0d0/(1.0d0 + w)
      mixder = (2.0d0 - w)/(2.0d0*(1.0d0 + w)**3)

c     ** denomsos = 1/s d(denom)/ds, sdenomsoss = s*d/ds(1/s d(denom)/ds))
      denomsos = c*zosn/(1.0d0 + zfac)**(3.0d0/4.0d0)*
     *     (1.0d0 + zfac + (1.0d0 + 3.0d0/2.0d0*zfac)*szsoz)
      sdenomsoss = c*zosn/(1.0d0 + zfac)**(7.0d0/4.0d0)*
     *     (-1.0d0 - zfac*(2.0d0 + zfac)
     *      + (1.0d0 + zfac/2.0d0*(5.0d0 + 3.0d0*zfac))*mixder
     *      + 3.0d0/2.0d0*zfac*(1.0d0 + zfac/2.0d0)*szsoz**2)

      Fsos = c/denom**2*(2.0d0 - zosn*
     *            ((1.0d0 - c*s2)*(1.0d0 + zfac)**(1.0d0/4.0d0) +
     *             (1.0d0 + c*s2)*(1.0d0 + 3.0d0/2.0d0*zfac)/
     *                      (1.0d0 + zfac)**(3.0d0/4.0d0)*szsoz))
      s2Fsossos = (-4.0d0*c*s2*denom*denomsos + (c*s2 + 1.0d0)*
     *                    (2.0d0*s2*denomsos**2 - denom*sdenomsoss))/
     *                      denom**3

c     ** Refinement function derivatives, 1/s dHx/ds, 1/s d/ds (1/s dHx/ds) 
c     ** We use that (1 - X) = a*X*s2
      Hxsos = (1.0d0 - X)*Fsos - (F - 1.0d0)*Xsos
c      Hxsossos = - 2.0d0*Fsos*Xsos + a*X*s2Fsossos -
c     *                  (F - 1.0d0)*Xsossos

c     ** vx formula for gradient dependent functional,
c      vx = vxlda*(Hx - s2*Hxsos) +
c     *     exlda*((4.0d0/3.0d0*s2-t)*Hxsos - 
c     *     (u-4.0d0/3.0d0*s**3)*s*Hxsossos)
      dexdr = 4.d0/3.d0*exlda*Hx - exlda*Hxsos*(4.d0/3.d0)*s2
      dexdgrad = n * exlda * 0.5d0 * Hxsos / scale 

c     *****************************
c        Correlation potential
c     *****************************
      
c     ** Correlation refinement function derivatives, dF/ds 
      Hcsos = Xsos*(1.0d0 - g)
c      Hcsossos = Xsossos*(1.0d0 - g)
      
c     ** vc formula for gradient dependent functional,
c     Generalized form of Eq. (24) in PRB 33, 8800 (1986) 
c      vc = vclda*(Hc - s2*Hcsos) +
c     *     eclda*((4.0d0/3.0d0*s2 - t)*Hcsos - 
c     *     (u - 4.0d0/3.0d0*s**3)*s*Hcsossos)

      decdr = vclda*Hc - eclda*Hcsos*(4.d0/3.d0)*s2
      decdgrad = n * eclda * 0.5d0 * Hcsos / scale

      return
      end

c     ******************************************
c       Local density approximation exchange
c
c       input
c       n        electron density [bohr**(-3)]
c
c       output
c       ex       exchange energy per electron [hartree]
c       vx       exchange potential [hartree]
c
c       Copyright (c) 2005, Rickard Armiento
c       All rights reserved.
c     ******************************************

      subroutine am05_ldax(n,ex,vx)
c     ** Input parameters
      real*8 n

c     ** Output parameters
      real*8 ex, vx

c     ** Constants
      real*8 pi
      parameter (pi = 3.141592653589793238462643383279502884197d0)

      vx = -(3.0d0*n/pi)**(1.0d0/3.0d0)
      ex = (3.0d0/4.0d0*vx)

      return
      end

c     ***********************************************
c       Local density approximation correlation
c
c       input
c       n        electron density [bohr**(-3)]
c
c       output
c       ec       correlation energy per electron [hartree]
c       vc       correlation potential [hartree]
c
c       As parameterized by Perdew Wang,
c         Phys. Rev. B 45, 13244 (1992) 
c       Based on monte carlo data by Ceperley Alder, 
c         Phys. Rev. Lett. 45, 566 (1980)
c
c       (Clean room implementation from paper)
c
c       Copyright (c) 2005, Rickard Armiento
c       All rights reserved.
c     ***********************************************
      subroutine am05_ldapwc(n,ec,vc)
      implicit none
c     ** Input parameters
      real*8 n

c     ** Output parameters
      real*8 ec, vc

c     ** Constants
      real*8 pi
      real*8 A0,a01,b01,b02,b03,b04
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (a01 = 0.21370d0)
      parameter (b01 = 7.5957d0)
      parameter (b02 = 3.5876d0)
      parameter (b03 = 1.6382d0)
      parameter (b04 = 0.49294d0)
c     ** Paper actually use this:
c      parameter (A0 = 0.031091d0)
c     ** But routines now "defacto standard" was distributed using:
      parameter (A0 = 0.0310907d0)

c     ** Local variables
      real*8 rsq
      real*8 Q0, Q1, Q1p, ecrs

      rsq = (3.0d0/(4.0d0*pi*n))**(1.0d0/6.0d0)

      ec = -2.0d0*A0*(1.0d0 + a01*rsq**2)* 
     *     log(1.0d0 + 1.0d0/
     *     (2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))))

      Q0 = -2.0d0*A0*(1.0d0 + a01*rsq**2)
      Q1 = 2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))
      Q1p = A0*(b01/rsq+2.0d0*b02+3.0d0*b03*rsq+4.0d0*b04*rsq**2)
      ecrs = -2.0d0*A0*a01*log(1.0d0 + 1.0d0/Q1)-Q0*Q1p/(Q1**2+Q1)

      vc=ec - rsq**2/3.0d0*ecrs
      
      end

c     ***********************************************
c       LambertW function. 
c
c       Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), 
c         Adv. in Comp. Math. 5(4):329-359. 
c       Implementation approach loosely inspired by the 
c       GNU Octave version by N. N. Schraudolph, but this 
c       implementation is only for real values and 
c       principal branch.
c
c       Copyright (c) 2005, Rickard Armiento
c       All rights reserved.
c     ***********************************************

      subroutine am05_lambertw(z,result)
      implicit none

c     input
      real*8 z
c     output
      real*8 result
c     local variables
      real*8 e,t,p    
      integer i

c     ** If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.0d-20) then
         result = z
         return
      endif

      e = exp(1.0d0)

c     ** Inital guess
      if( abs(z + 1.0d0/e) .gt. 1.45d0 ) then
c        ** Asymptotic expansion at 0 and Inf
         result = log(z)
         result = result - log(result)
      else
c        ** Series expansion about -1/e to first order
         result = 1.0d0*sqrt(2.0d0*e*z + 2.0d0) - 1.0d0
      endif

c     ** Find result through iteration
      do i=1,10
         p = exp(result)
         t = result*p - z
         if( result .ne. -1.0d0 ) then
            t = t/(p*(result + 1.0d0) - 
     *           0.5d0*(result + 2.0d0)*t/(result + 1.0d0))
         else
            t = 0.0d0
         endif
         result = result - t
         if(abs(t) < (2.48d0*1.0d-14)*(1.0d0 + abs(result))) then
            return
         endif
      enddo
c     ** This should never happen!;
      write(*,*) 'am05_lambertw: iteration limit reached.' 
      write(*,*) 'Should never happen: execution aborted.'
      stop
      end


