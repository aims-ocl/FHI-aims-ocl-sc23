      subroutine vwnrpa_c_derivs(rho, ec, dedcp, n_spin)

      implicit none
      integer :: n_spin, i_spin
      real*8, dimension(n_spin) :: dens, rho
      real*8 :: dens_tot
      real*8 :: rs, rx, zeta, decz
      real*8 :: ecp, ecfp, decfp, ec,
     &          denmin, decp, dedcptmp
c    output derivatives
      real*8, dimension(n_spin) :: dedcp 
c    constants
      real*8 :: third,fthrd,fthrd2,fdenom,pi,pi43,pi49,
     &         twthrd,sevsix,fvthrd,ffifth,svthrd,zero,one,two,three
     &         ,four,five,six,dnine
c functions ?!
      real*8 :: fspin, fspol, fspol1, fspol2, gvwn, dgvwn
c parameters for functions
      real*8 :: z, x, a, b, c, x0
      real*8 :: cep, cef, cjp
      dimension cep(4),cef(4),cjp(4)
      
      data cep/ 0.0621814,13.0720,42.7198,-0.409286/
     &    ,cef/ 0.0310907,20.1231,101.578,-0.743294/
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
      parameter (denmin=1.d-20)
      
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
        dens_tot =  rho(1)
        zeta = 0.d0
      else
        dens(1) = rho(1)
        dens(2) = rho(2)
        dens_tot = rho(1)+rho(2)
      endif

      if (dens_tot .lt. denmin) then
         ec = 0.d0
         dedcp = 0.d0
      else
        if(n_spin .eq. 2) then
        zeta = (dens(1)-dens(2))/dens_tot
         if ( zeta.gt.+one ) zeta = +one
         if ( zeta.lt.-one ) zeta = -one
        endif
      
      
C        if ( dens_tot.gt.detol ) then
          rs = one/(pi43*dens_tot)**third
C        else
C          rs = one/(pi43*detol*detol)**third
C        end if
      
      
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
        
        dedcptmp = ecp-decp*rx/six

        decz = 0.d0
      
      do i_spin=1, n_spin 
          dedcp(i_spin)=dedcptmp
      end do
      
c     
c  skip the following part if there is no spin polarization
c     
       if ( n_spin.eq.2 ) then

        ecfp = (gvwn(rx,cef(1),cef(2),cef(3),cef(4))/two)
     $       - ecp
      
        decfp = (dgvwn(rx,cef(1),cef(2),cef(3),cef(4))/two)
     $       - decp
      
 
        ec = ec + fspol2(zeta)*(ecfp)
   
        decz = fspol1(zeta)
     $       * (ecfp)

        dedcp(1) = dedcptmp + fspol2(zeta)*(ecfp
     &       - (rx/six)*(decfp))+(1.d0-zeta)*decz

        dedcp(2) = dedcptmp + fspol2(zeta)*(ecfp
     &       - (rx/six)*(decfp))-(1.d0+zeta)*decz
        

      

       end if
      endif

c      do ispin=1, 2
c          vxc(ispin) = vc(ispin)
c      end do
      
      return
      end

