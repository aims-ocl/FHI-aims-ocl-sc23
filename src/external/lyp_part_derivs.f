!-------------------------------------------------------
! LYP functional and partial derivatives
! Not dependent on the laplacian of the density
! References:
! Johnson, Gill and Pople - JCP 98. 5612 (1993)
! ------------------------------------------------------

      subroutine lyp_part_derivs(rho, rho_gradient, 
     +			   lypc, dlypdr, dlypdgr, n_spin)
      
      use constants
      implicit none
      integer n_spin
      double precision  rho(n_spin), rho_gradient(3,n_spin), lypc,
     +                  dlypdr(n_spin), dlypdgr(3,n_spin)


c Internal variables
      integer i_spin,i_coord
      double precision beta, thd, tthd, thrhlf, half, fothd,
     +                 d(2),gd(3,2), gdm(2),denmin,dt, 
     +                 g(2),x(2),a,b,c,dd,onzthd,gdmin, 	     
     +                 den,omega, domega, delta, ddelta,cf,
     +                 gam11, gam12, gam22, LYPa, LYPb1,
     +                 LYPb2,dLYP11,dLYP12,dLYP22,LYP,
     +                 dd1g11,dd1g12,dd1g22,dd2g12,dd2g11,dd2g22,
     +                 dLYPdd(2),dg11dd(3,2),dg22dd(3,2),
     +                 dg12dd(3,2),dLYPgd(3,2)
  
c numerical parameters 
      parameter ( thd = 1.d0/3.d0, tthd=2.d0/3.d0 )
      parameter ( thrhlf=1.5d0, half=0.5d0,
     +            fothd=4.d0/3.d0, onzthd=11.d0/3.d0)

c Constants for LYP functional (a.u.) 
      parameter(a=0.04918d0, b=0.132d0, c=0.2533d0, dd=0.349d0)

c Lower bounds of density and its gradient to avoid divisions by zero
      parameter (denmin=1.d-8)
      parameter (gdmin=1.d-8)

c Translate density and its gradient to new variables
      if (n_spin .eq. 1) then
        d(1) = half * rho(1)
        d(1) = max(denmin,d(1))
        d(2) = d(1)
        dt = max( denmin, rho(1) )
        do i_coord = 1,3
          gd(i_coord,1) = half * rho_gradient(i_coord,1)    
          gd(i_coord,2) = gd(i_coord,1)
        enddo 
      else
        d(1) = rho(1)
        d(2) = rho(2)
        do i_spin=1,2
         d(i_spin) = max (denmin,d(i_spin))
        enddo
        dt = max( denmin, rho(1)+rho(2) )  
        do i_coord = 1,3
          gd(i_coord,1) = rho_gradient(i_coord,1)
          gd(i_coord,2) = rho_gradient(i_coord,2)
        enddo
      endif

      gdm(1) = sqrt( gd(1,1)**2 + gd(2,1)**2 + gd(3,1)**2 )
      gdm(2) = sqrt( gd(1,2)**2 + gd(2,2)**2 + gd(3,2)**2 )
 
      do i_spin=1,2
       gdm(i_spin)= max(gdm(i_spin),gdmin)
      enddo

      den=1+dd*dt**(-thd)
      omega=dt**(-onzthd)*exp(-c*dt**(-thd))/den
      delta=c*dt**(-thd)+dd*dt**(-thd)/den
      cf=3.*(3*pi**2)**tthd/10.
      gam11=gdm(1)**2
      gam12=gd(1,1)*gd(1,2)+gd(2,1)*gd(2,2)+gd(3,1)*gd(3,2)
      gam22=gdm(2)**2
      LYPa=-4*a*d(1)*d(2)/(den*dt)
      LYPb1=2**onzthd*cf*a*b*omega*d(1)*d(2)
      LYPb2=d(1)**(8./3.)+d(2)**(8./3.)
      dLYP11=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)
     +*d(1)/dt)-d(2)**2)
      dLYP12=-a*b*omega*(d(1)*d(2)/9.*(47.-7.*delta)
     +-fothd*dt**2)
      dLYP22=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)*
     +d(2)/dt)-d(1)**2)

c    Density of lyp energy
      LYP=(LYPa-LYPb1*LYPb2+dLYP11*gam11+dLYP12*gam12
     ++dLYP22*gam22)/dt

c   Correlation energy derivatives
       domega=-thd*dt**(-fothd)*omega*(11.*dt**thd-c-dd/den)
       ddelta=thd*(dd**2*dt**(-5./3.)/den**2-delta/dt)

c   Second derivatives with respect to the density
       dd1g11=domega/omega*dLYP11-a*b*omega*(d(2)/9.*
     + (1.-3.*delta-2*(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     + ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2))

       dd1g12=domega/omega*dLYP12-a*b*omega*(d(2)/9.*
     + (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)

      dd1g22=domega/omega*dLYP22-a*b*omega*(1./9.*d(2)
     + *(1.-3.*delta-(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     + ((3.+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2)-2*d(1))

       
      dd2g22=domega/omega*dLYP22-a*b*omega*(d(1)/9.*
     + (1.-3.*delta-2*(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     + ((3+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2))
      
 
      dd2g12=domega/omega*dLYP12-a*b*omega*(d(1)/9.*
     + (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)
      
      dd2g11=domega/omega*dLYP11-a*b*omega*(1./9.*d(1)
     + *(1.-3.*delta-(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     + ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2)-2*d(2))


        dLYPdd(1)=-4*a/den*d(1)*d(2)/dt*
     + (thd*dd*dt**(-fothd)/den
     + +1./d(1)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     + (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(2)*(onzthd*
     + d(1)**(8./3.)+d(2)**(8./3.)))+dd1g11*gam11+
     + dd1g12*gam12+dd1g22*gam22


       dLYPdd(2)=-4*a/den*d(1)*d(2)/dt*(thd*dd*dt**(-fothd)/den
     + +1./d(2)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     + (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(1)*(onzthd*
     + d(2)**(8./3.)+d(1)**(8./3.)))+dd2g22*gam22+
     + dd2g12*gam12+dd2g11*gam11


c derivatives with respect to the density gradient

        do i_spin=1,2
          do i_coord=1,3
           dg11dd(i_coord,i_spin)=2*gd(i_coord,i_spin)
           dg22dd(i_coord,i_spin)=2*gd(i_coord,i_spin)
          enddo
        enddo
        do i_coord=1,3
          dLYPgd(i_coord,1)=dLYP11*dg11dd(i_coord,1)
     +     +dLYP12*gd(i_coord,2)
          dLYPgd(i_coord,2)=dLYP22*dg22dd(i_coord,2)
     +     +dLYP12*gd(i_coord,1)
        enddo

c assinging the values to the output variables

       lypc=LYP
       do i_spin=1,n_spin
        dlypdr(i_spin)=dLYPdd(i_spin) 
        do i_coord=1,3
         dlypdgr(i_coord,i_spin)=dLYPgd(i_coord,i_spin)/2
        enddo
       enddo
       
       return
       end 



