CVB---------------------------------------------------------------------
C
C  Old version of the spin-polarized version of Perdew / Zunger's 
C  1981 LDA interpolation for the XC functional.
C
C  I prefer this version because it gives explicit and direct expressions 
C  for Exc and Vxc. The alternative which I had was from the Daresbury
C  "DFT repository", but that requires a hugely bloated interface, as
C  well as consisting of entirely unreadable source code.
C
C  Cheers, VB
C
CVB end

        subroutine stvxc_spin(rho, zeta, vxc, dvxc,
     &                       exc)
        implicit none

        real*8   exc
        real*8   rho, zeta, vxc,
     &           dvxc
c-----------------------------------------------------------------------
c       Calculate the exchange-correlation energy and the exchange
c       correlation potentials for a spin-polarized electron gas within
c       the LSDA. Use the parametrization of the XC-energy of the homo-
c       geneous electron gas described by Perdew and Zunger in 
c       Phys. Rev. B 23, 5048 (1981).
c
c       Input:
c       rho     total charge density (a.u.)
c       zeta    spin polarization
c
c       Output:
c       vxc_vec         average XC-potential (a.u.)
c       dvxc_vec        XC-potential for up-spin
cmb                     electrons minus the spin-down XC potential (a.u.)
c       exc             exchange-correlation energy density in a.u. (Hartrees)
c
c       Version  1.0, Author: E. Pehlke, TU Muenchen, August 14, 1997.
c
c       modified (reduced to basics) by VB 2006
c
c-----------------------------------------------------------------------
c
c       Tabulate XC-parameters by Perdew and Zunger.
c       U: unpolarized (paramagnetic)
c       P: polarized   (ferromagnetic)
c
        real*8   A_U, B_U, C_U, D_U, b1_U, b2_U, g_U
        real*8   A_P, B_P, C_P, D_P, b1_P, b2_P, g_P
        parameter       (A_U=0.0311D0,  B_U=-0.0480D0)
        parameter       (A_P=0.01555D0, B_P=-0.0269D0)
        parameter       (C_U=0.0020D0,  D_U=-0.0116D0)
        parameter       (C_P=0.0007D0,  D_P=-0.0048D0)
        parameter       (b1_U=1.0529D0, b2_U=0.3334D0, g_U=-0.1423D0)
        parameter       (b1_P=1.3981D0, b2_P=0.2611D0, g_P=-0.0843D0)
c
c       Other constants.
c
        real*8  x4b3, x7b6, x1b3, x2b3, xn, pi, eps, one, zero, two
c       parameter       (x4b3=4.0D0/3.0D0, x7b6=7.0D0/6.0D0)
c       parameter       (x1b3=1.0D0/3.0D0, x2b3=2.0D0/3.0D0)
c       parameter       (xn=2.0D0**(4.0D0/3.0D0)-2.0D0)
        parameter       (pi=3.141592653589793238D0)
        parameter       (eps=1.D-13, one=1.0D0, zero=0.0D0, two=2.0D0)
c
        real*8   Cx_U, Cx_P
c       parameter       (Cx_U=3.0D0/(4.0D0*pi)*(9.0D0*pi/4.0D0)**x1b3)
c       parameter       (Cx_P=Cx_U*2**x1b3)
c
c
        real*8   rs, rrs, lrs, azeta, f, fprime, exc_U, exc_P,
     &           vxc_P, vxc_U
c
c
        exc = zero
        x4b3=4.0D0/3.0D0 
        x7b6=7.0D0/6.0D0
        x1b3=1.0D0/3.0D0
        x2b3=2.0D0/3.0D0
        xn=2.0D0**(4.0D0/3.0D0)-2.0D0
        Cx_U=3.0D0/(4.0D0*pi)*(9.0D0*pi/4.0D0)**x1b3
        Cx_P=Cx_U*2**x1b3
c
           if (rho.lt.-eps) then
              rho=0.0
           else if (rho.lt.eps) then
              vxc=zero
              dvxc=zero
           else
c
c             Calculate Wigners density parameter R_s:
c
              rs=(x4b3*pi*rho)**(-x1b3)
              rrs=sqrt(rs)
              lrs=log(rs)
c
c             fprime is the derivative of f(zeta) with respect to zeta.
c             Note that f(zeta) is an even function of zeta.
c
              azeta=abs(zeta)
              if (azeta.lt.eps) then
                 f=zero
                 fprime=zero
              else if (azeta.lt.one-eps) then
                 f=((one+zeta)**x4b3+(one-zeta)**x4b3-two)/xn
                 fprime=((one+zeta)**x1b3-(one-zeta)**x1b3)*x4b3/xn
              else
                 f=one
                 fprime=two**x1b3*x4b3/xn
                 if (zeta.lt.zero) fprime=-fprime
              endif
c
c             Calculate exchange energy and potentials for U and P:
c
              exc_U=-Cx_U/rs
              exc_P=-Cx_P/rs
              vxc_U=x4b3*exc_U
              vxc_P=x4b3*exc_P
c
c             Add correlation energy and potentials for U and P:
c
              if (rs.lt.one) then
                 exc_U=exc_U+A_U*lrs+B_U+C_U*rs*lrs+D_U*rs
                 exc_P=exc_P+A_P*lrs+B_P+C_P*rs*lrs+D_P*rs
                 vxc_U=vxc_U+A_U*lrs+(B_U-x1b3*A_U)+x2b3*C_U*rs*lrs+
     &                (x2b3*D_U-x1b3*C_U)*rs
                 vxc_P=vxc_P+A_P*lrs+(B_P-x1b3*A_P)+x2b3*C_P*rs*lrs+
     &                (x2b3*D_P-x1b3*C_P)*rs
              else
                 exc_U=exc_U+g_U/(one+b1_U*rrs+b2_U*rs)
                 exc_P=exc_P+g_P/(one+b1_P*rrs+b2_P*rs)
                 vxc_U=vxc_U+g_U*(one+x7b6*b1_U*rrs+x4b3*b2_U*rs)/
     &                (one+b1_U*rrs+b2_U*rs)**2
                 vxc_P=vxc_P+g_P*(one+x7b6*b1_P*rrs+x4b3*b2_P*rs)/
     &                (one+b1_P*rrs+b2_P*rs)**2
              endif

              exc=exc_U+f*(exc_P-exc_U)

              vxc=vxc_U+f*(vxc_P-vxc_U)
     &                   -zeta*fprime*(exc_P-exc_U)

              dvxc=2*fprime*(exc_P-exc_U)
           endif
c
        return
        end
