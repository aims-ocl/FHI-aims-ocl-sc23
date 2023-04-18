      subroutine cepvwn_gauss(dens,ceng,vc,npolrz)

      implicit real*8(a-h,o-z)
      integer npolrz
      dimension dens(3)
      dimension vc(2)

      real*8 :: ec, ceng
      dimension grmod(3),dlap(3),dldlp(3)
      real*8 :: rs,rx,zeta
      dimension dkf(2),sperd(2),faca(2),facb(2),
     &          facc(2),facd(2),fperd(2),sdfds(2),
     &          dsdfd(2),tperd(2),uperd(2),eperd(2)
      real*8 :: ecp,decp,vcp,ecfp,decfp,
     &          acp,dacp,vcfp
      
      dimension cap(4),cep(4),cef(4),cjp(4)
      
      data dperdm,dperdn,dperdb,dperdc,alalph,albeta
     $     / 0.06666666666667d0,0.06666666666667d0,
     $     14.d0,0.2d0,0.05d0,5.d0/
      data detol,zetol,ftyl,bkbeta/1.d-20,1.d-10,0.11d0,0.0042d0/
      data cap/-0.0337737,1.13107,13.0045,-0.0047584/
     &     cep/ 0.0621814,13.0720,42.7198,-0.409286/
     &    ,cef/ 0.0310907,20.1231,101.578,-0.743294/
     &    ,cjp/0.023266,7.389e-06,8.723,0.472/
      data third,fthrd,fthrd2,fdenom,pi,pi43,pi49
     $     / 0.33333333333333d0,1.33333333333333d0,5.1297633d0,
     $     0.51984204d0,3.14159265358979d0,
     $     4.18879020478638d0,1.63696380104572d0 /
      data twthrd,sevsix,fvthrd,ffifth,svthrd,zero,one,two,three,four,
     $     five,six,dnine,expcut
     $     / 0.66666666666667d0,1.16666666666667d0,1.66666666666667d0,
     $     0.8d0,2.33333333333333d0,
     $     0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,9.d0,80.d0 /
      data eps
     $     /1.d-100/
      integer i_spin
      
c----------------------------------------------------------------------
      
      fspin(z)=(one+z)**fthrd+(one-z)**fthrd
      fspol(z)=(one+z)**third
      fspol1(z)=(fspol(+z)-fspol(-z))*fthrd2/two
      fspol2(z)=(fspin(z)-two)/fdenom
      
      gvwn(x,a,b,c,x0)=a*log(x*x/(x-x0)**2)-a/(x0*x0+b*x0+c)*
     &     ((x0*x0+c)*log((x*x+b*x+c)/(x-x0)**2)+
     &     (x0*x0-c)*atan(sqrt(four*c-b*b)/(two*x+b))*two*b/
     &     sqrt(four*c-b*b))
      dgvwn(x,a,b,c,x0)=two*a*(one/x-(one+b/(x-x0))*x/(x*x+b*x+c))
c if total density is zero, nothing should be calculated
      if(dens(3) .le. eps) then
        ceng=0.d0
        vc=0.d0
      else

      if ( dens(3).gt.detol ) then
          rs = one/(pi43*dens(3))**third
      else
          rs = one/(pi43*detol*detol)**third
      end if
      
      
      if ( dens(3).gt.detol ) then
        zeta = (dens(1)-dens(2))/dens(3)
        if ( zeta.gt.+one ) zeta = +one
        if ( zeta.lt.-one ) zeta = -one
      else
        zeta = zero
      end if

      

c***********************************************************************
c                                                                      *
c  local correlation energy and potentials                             *
c                                                                      *
c***********************************************************************
        rx = sqrt(rs)
      
        ecp = gvwn(rx,cep(1),cep(2),cep(3),cep(4))/two
      
        decp = dgvwn(rx,cep(1),cep(2),cep(3),cep(4))/two
      
        vcp = ecp-decp*rx/six
      
        ec = ecp
      
      do i_spin=1, 2
          vc(i_spin)=vcp
      end do
      
c     
c  skip the following part if there is no spin polarization
c     
      if ( npolrz.eq.0 ) goto 3053
      
        ecfp = (gvwn(rx,cef(1),cef(2),cef(3),cef(4))/two)
     $       - ecp
      
        decfp = (dgvwn(rx,cef(1),cef(2),cef(3),cef(4))/two)
     $       - decp
      
c        acp = (gvwn(rx,cap(1),cap(2),cap(3),cap(4))/two)
c     $       * three/fthrd2
      
c        dacp=(dgvwn(rx,cap(1),cap(2),cap(3),cap(4))/two)
c     $       * three/fthrd2
      
        vcfp = fspol2(zeta)*(ecfp
     &       - (rx/six)*decfp)
      
        ec = ec + fspol2(zeta)*(ecfp)
      
      vc(1) = vcp + vcfp + (one-zeta)
     $       * ( fspol1(zeta)
     $       * ecfp )
      
      vc(2) = vcp + vcfp - (one+zeta)
     $       * ( fspol1(zeta)
     $       * ecfp)
      
 3053 continue
      
        ceng = ec
      
      
      end if

      return
      end

