c Header:
c**********************************************************************
c Lee Yang Parr correlation energy functional
c one-dimensional densities only
c no provisions taken against division by zero
c 
c see e.g.
c C. Lee et al. Phys Rev B 37 785 (1988)
c
c Hartree a.u.
c               
c input
c tpot .... T  evaluate correlation energy and potential
c           F  evaluate energy only (a posteriori scheme)
c x ....... dummy
c dp ...... spin up density
c dp1 ..... grad(dp)
c dp2 ..... laplace(dp)
c dm ...... spin down density
c ...
c
c output
c ec ...... correlation energy per electron
c vcp0 .... correlation potential for spin up
c vcm0 .... correlation potential for spin down
c
c Martin Fuchs, FHI der MPG, 06-1995
c**********************************************************************
      subroutine corlyp(tpot,x,dp,dm,dp1,dm1,dp2,dm2,ec,vcp0,vcm0)
c
      implicit real*8 (a-h,o-z)
      logical  tpot
      parameter(aa=0.04918d0,bb=0.132d0,cc=0.2533d0,dd=0.349d0,
     &         c1=-4*aa,c2=dd,c3=2*bb,c4=cc,c5=4.55779986d0,
     &         c6=1.d0/72.d0,c7=1.d0/18.d0,c8=0.125d0,
     &         t13=1.d0/3.d0,t53=5*t13,t43=4*t13,t83=2*t43,
     &         t89=8.d0/9.d0,c9=t43+t89,zero=0.d0)
c
      d0= dp+ dm
      dxsq= 1.d0/(d0*d0)
      d1= dp1+ dm1
      d2= dp2+ dm2
      d0xt13= d0**(-t13)
      d0xt53= d0xt13*d0xt13/d0
      dpt53= dp**t53
      dmt53= dm**t53

c polarization factor
      z= c1*(dp*dm)*dxsq

c scaling function
      sc= 1.d0/(1.d0+ c2*d0xt13)
      h= c3*d0xt53*exp(-c4*d0xt13)

c kinetic energy density expansion
      ga= c5*(dp*dpt53+ dm*dmt53)

      gb= c6*(dp1*dp1-dp*dp2+ dm1*dm1-dm*dm2) 
     &   +c7*(dp*dp2+ dm*dm2)
     &   +c8*(d0*d2- d1*d1)

c calculate potential
      if(tpot) then

        gafp= t83*c5*dpt53
        gafm= t83*c5*dmt53

        scf= t13*c2*d0xt13/d0*sc*sc
        sc2= scf*(d2+ 2*(scf/sc- 2*t13/d0)*d1*d1)

        chf= t13*(c4*d0xt13 -5)/d0
        hf= chf*h
        hff= h*(chf**2+ t13*(5.d0-4*t13*c4*d0xt13)*dxsq)
        h2= (hf*d2+ hff*d1*d1)
  
        zfp= (c1*dm- 2*z*d0)*dxsq
        zfm= (c1*dp- 2*z*d0)*dxsq
        yz= z/c1
        yy1= dp*dm1+dm*dp1
        yz1= (yy1-2*yz*d1*d0)*dxsq
        yz2= (2*yz*d1*d1- 2*(yz1*d1 +yz*d2)*d0
     &      -2*d1*yy1/d0+ (dp*dm2+2*dp1*dm1+dm*dp2))*dxsq
        z1= c1*yz1
        z2= c1*yz2
  
        ya= sc*z*d0
        yafp= sc*(d0*zfp+z)+ z*d0*scf
        yafm= sc*(d0*zfm+z)+ z*d0*scf

        yb= sc*z*h
        ybfp= sc*(h*zfp+z*hf)+ z*h*scf
        ybfm= sc*(h*zfm+z*hf)+ z*h*scf
        yb1= sc*(h*z1+z*hf*d1)+ z*h*scf*d1
        yb2= (sc*hf+h*scf)*d1*z1+ h*sc*z2
     &     + (sc*z1+z*scf*d1)*hf*d1+ z*sc*h2
     &     + (z*hf*d1+h*z1)*scf*d1+ z*h*sc2
  
c collect contributions
        vcp0= yafp+ ybfp*(ga+gb)
     &      + yb*(gafp+2*c8*(c9*dp2+2*dm2))
     &      + yb1*2*c8*(c9*dp1+2*dm1)
     &      + yb2*c8*(t43*dp+dm)

        vcm0= yafm+ ybfm*(ga+gb)
     &      + yb*(gafm+2*c8*(c9*dm2+2*dp2))
     &      + yb1*2*c8*(c9*dm1+2*dp1)
     &      + yb2*c8*(t43*dm+dp)
  
      else

        vcp0 = zero
        vcm0 = zero

      endif
  
c correlation energy per electron
      ec= z*sc*(d0+ h*(ga+gb))/d0

      return
      end
c
