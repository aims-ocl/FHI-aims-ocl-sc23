      subroutine vwn_c_derivs(rho, ec, dedcp, decz ,n_spin)

      implicit none
      integer :: n_spin, detol, i_spin
      real*8, dimension(n_spin) :: dens, rho
      real*8 :: dens_tot
      real*8 :: rs, rx, zeta
      real*8 :: ecp, ecfp, decfp, ec,
     &          acp,dacp, denmin, decp
c    output derivatives
      real*8 :: dedcp, decz 
c    constants
      real*8 :: third,fthrd,fthrd2,fdenom,pi,pi43,pi49,
     &         twthrd,sevsix,fvthrd,ffifth,svthrd,zero,one,two,three
     &         ,four,five,six,dnine
c functions ?!
      real*8 :: fspin, fspol, fspol1, fspol2, gvwn, dgvwn
c parameters for functions
      real*8 :: z, x, a, b, c, x0
      real*8 :: cap, cep, cef, cjp
      dimension cap(4),cep(4),cef(4),cjp(4)
      
      data detol/1.d-20/
      data cap/-0.0337737,1.13107,13.0045,-0.0047584/
     &    ,cep/ 0.0621814,3.72744,12.9352,-0.1049800/
     &    ,cef/ 0.0310907,7.06042,18.0578,-0.3250000/
     &    ,cjp/0.023266,7.389e-06,8.723,0.472/
      data third,fthrd,fthrd2,fdenom,pi,pi43,pi49
     $     / 0.33333333333333d0,1.33333333333333d0,5.1297633d0,
     $     0.51984204d0,3.14159265358979d0,
     $     4.18879020478638d0,1.63696380104572d0 /
      data twthrd,sevsix,fvthrd,ffifth,svthrd,zero,one,two,three,four,
     $     five,six,dnine
     $     / 0.66666666666667d0,1.16666666666667d0,1.66666666666667d0,
     $     0.8d0,2.33333333333333d0,
     $     0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,9.d0 /
      parameter (denmin=1.d-8)
      
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

      if (n_spin .eq. 1) then
        dens(1) = rho(1)
        dens(1) = max(denmin,dens(1))
        dens_tot = max( denmin, rho(1) )
        zeta = 0.d0
      else
        dens(1) = rho(1)
        dens(2) = rho(2)
        do i_spin=1,2
         dens(i_spin) = max (denmin,dens(i_spin))
        enddo
        dens_tot = max( denmin, rho(1)+rho(2) )
        zeta = (dens(1)-dens(2))/dens_tot
        if ( zeta.gt.+one ) zeta = +one
        if ( zeta.lt.-one ) zeta = -one
      endif
      
      
        if ( dens_tot.gt.detol ) then
          rs = one/(pi43*dens_tot)**third
        else
          rs = one/(pi43*detol*detol)**third
        end if
      
      
c***********************************************************************
c                                                                      *
c  local correlation energy and potentials                             *
c                                                                      *
c***********************************************************************
        rx = sqrt(rs)
      
        ecp = gvwn(rx,cep(1),cep(2),cep(3),cep(4))/two
c        write(*,*) rs, rx, cep(1)
      
        decp = dgvwn(rx,cep(1),cep(2),cep(3),cep(4))/two
      
c        vcp = ecp-decp*rx/six

        ec = ecp
        
        dedcp = ecp-decp*rx/six

        decz = 0.d0
      
c      do ispin=1, 2
c          vc(ispin)=vcp
c      end do
      
c     
c  skip the following part if there is no spin polarization
c     
      if ( n_spin.eq.2 ) then

        ecfp = (gvwn(rx,cef(1),cef(2),cef(3),cef(4))/two)
     $       - ecp
      
        decfp = (dgvwn(rx,cef(1),cef(2),cef(3),cef(4))/two)
     $       - decp
      
        acp = (gvwn(rx,cap(1),cap(2),cap(3),cap(4))/two)
     $       * three/fthrd2
      
        dacp=(dgvwn(rx,cap(1),cap(2),cap(3),cap(4))/two)
     $       * three/fthrd2
      
    
        ec = ec + fspol2(zeta)*(ecfp*(zeta**4)
     $       + acp*(one-zeta**4))

        dedcp = dedcp + fspol2(zeta)*(ecfp*(zeta**4)
     &       + acp*(one-(zeta**4))
     &       - (rx/six)*(decfp*(zeta**4)
     &       + dacp*(one-(zeta**4))))
        
        decz = ( fspol1(zeta)
     $       * (ecfp*(zeta**4) +acp*(one-(zeta**4)))
     $       + four*fspol2(zeta)*(ecfp*(zeta**3)
     $       - acp*(zeta**3)) )

      
c      vc(1) = vcp + vcfp + (one-zeta)
c     $       * ( fspol1(zeta)
c     $       * (ecfp*(zeta**4) +acp*(one-(zeta**4)))
c     $       + four*fspol2(zeta)*(ecfp*(zeta**3)
c     $       - acp*(zeta**3)) )
      

c      vc(2) = vcp + vcfp - (one+zeta)
c     $       * ( fspol1(zeta)
c     $       * (ecfp*(zeta**4) +acp*(one-(zeta**4)))
c     $       + four*fspol2(zeta)*(ecfp*(zeta**3)
c     $       - acp*(zeta**3)) )

      end if


c      do ispin=1, 2
c          vxc(ispin) = vc(ispin)
c      end do
      
      return
      end

