c**********************************************************************
c LDA - Ceperley - Alder correlation potential and energy
c as parameterized by John P. Perdew and Yue Wang  
c Phys. Rev. B 45, 13244–13249 (1992) 

c Hartree a.u.
c Now is for n_spin=1, so only eu is considered, 
c for n_spin=2, we need consider zeta.

c written by shanghui @ FHI, 2013.07.31
c**********************************************************************
c
      subroutine pw_lda(rh,exc,ex,ec,vxc,fxc)
      implicit none 
      real*8 rh,rs,eu,eurs,eursrs
      real*8 ex,ec,vx,vc,fx,fc,exc,vxc,fxc

      real*8, parameter:: 
     &  PI=3.14159265358979312d0 ,
     &  eps=1.d-100,
     &  crs=0.620350490899400087d0, ! crs=(3/4pi)^(1/3)  
      !define : 3/(4pi rs^3)=rh ===> so rs= (3/4pi)^(1/3)*rh^(-1/3)
     &  ax=-0.738558766382022406d0 ! ax=-3/4*(3/pi)^(1/3)
      !deine  : ex=ax*rh^(1/3)


      !-------begin PW-LDA parameters for n_spin=1------------------
      real*8, parameter:: 
     &   a=0.0310907d0 ,
     &  a1=0.21370d0 ,
     &  b1=7.5957d0 ,
     &  b2=3.5876d0 ,
     &  b3=1.6382d0 ,
     &  b4=0.49294d0
      !-------end PW-LDA parameters for n_spin=1--------------    


        if(rh .le. eps) then
          fxc=0.d0
          vxc=0.0d0
          exc=0.d0
          goto 100
        endif

      rs=crs*rh**(-0.3333333333333333d0)

      call gcor_shanghui(a,a1,b1,b2,b3,b4,1.00d0,rs,eu,eurs,eursrs)

        ex=ax*rh**(.3333333333333333d0) 
        ec=eu
        exc=ex+ec        

        vx=4.d0/3.d0*ex
        vc= ec-rs*eurs/3.d0 
        vxc=vx+vc

        fx=4.0d0/9.0d0*ax*rh**(-2.0d0/3.0d0) 
        fc=-rs**4.0d0 
        fc= fc * 4.0d0 * PI 
        fc= fc * ( 2.0d0*eurs/3.0d0 - rs*eursrs/3.d0 )/9.0d0
        fxc=fx+fc

  100 return
      end

      subroutine gcor_shanghui(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs,ggrsrs)
      implicit none
      real*8 a,a1,b1,b2,b3,b4,p,rs,gg,ggrs,ggrsrs
      real*8 p1,q0,q0_1,rs12,rs32,rsp,q1,q1_1,q1_2,q2
 
      p1 = p + 1.d0
      q0 = -2.d0*a*(1.d0+a1*rs)
      q0_1= -2.d0*a*a1

      rs12 = dsqrt(rs)
      rs32 = rs12**3
      rsp = rs**p
      q1 = 2.d0*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
      q2 = dlog(1.d0+1.d0/q1)
      gg = q0*q2

      q1_1 = a*(b1/rs12+2.d0*b2+3.d0*b3*rs12+2.d0*b4*p1*rsp)
      q1_2 = a*(-0.5d0*b1/rs32+1.5d0*b3/rs12+2.d0*b4*p*p1*rsp/rs )

      ggrs = -2.d0*a*a1*q2-q0*q1_1/(q1**2+q1)

      ggrsrs= 2.0d0*a*a1*q1_1/(q1**2+q1) - 
     . ( (q0_1*q1_1+q0*q1_2)*(q1**2+q1)-q0*q1_1**2*(2.0d0*q1+1.0d0) )/ 
     .       (q1**2+q1)**2  

      return
      end
