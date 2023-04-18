c**********************************************************************
c LDA - Ceperley - Alder exchange-correlation potential and energy
c as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
c Hartree a.u.
c now is for n_spin=1, need to extended to n_spin=2 in future.
c written by shanghui @ FHI, 2012.09.28
c**********************************************************************
c
      subroutine pz_lda(rh,exc,ex,ec,vxc,fxc)
      !----------------------------------
      !input:  rh 
      !output: exc,vxc,fxc
      !---------------------------------
      implicit none 
      real*8 rh,rs,sqrs,den,rsl
      real*8 ex,ec,vx,vc,fx,fc,exc,vxc,fxc

      real*8, parameter:: 
     &  PI=3.14159265358979312d0 ,
     &  eps=1.d-100,
     &  crs=0.620350490899400087d0, ! crs=(3/4pi)^(1/3)  
      !define : 3/(4pi rs^3)=rh ===> so rs= (3/4pi)^(1/3)*rh^(-1/3)
     &  ax=-0.738558766382022406d0 ! ax=-3/4*(3/pi)^(1/3)
      !deine  : ex=ax*rh^(1/3)


      !-------PZ-CA parameters for n_spin=1------------------
      real*8, parameter:: 
     &  gamma_U=-0.1423d0,
     &  beta_1_U=1.0529d0, beta_2_U=0.3334d0,
     &  A_U= 0.0311d0,
     &  B_U=-0.048d0,
     &  C_U= 0.002d0,
     &  D_U=-0.0116d0
      !-------end PZ-CA parameters for n_spin=1--------------    


        if(rh .le. eps) then
          fxc=0.d0
          vxc=0.0d0
          exc=0.d0
          goto 100
        endif

      rs=crs*rh**(-0.3333333333333333d0)

      if(rs .gt. 1.0d0) then

        sqrs=dsqrt(rs)
        den=1.0d0 + beta_1_U*sqrs + beta_2_U*rs

        ex=ax*rh**(.3333333333333333d0) 
        ec= gamma_U/den
        exc=ex+ec        

        vx=4.d0/3.d0*ex
        vc=ec*( 1+ (7.0d0/6.0d0)*beta_1_U*sqrs+  
     &      (4.0d0/3.0d0)*beta_2_U*rs )/den
        vxc=vx+vc

        fx=4.0d0/9.0d0*ax*rh**(-2.0d0/3.0d0) 
        fc=-ec*PI*rs**(7.0d0/2.0d0)
     &  *(-7.0d0*beta_1_U**(2.0d0)*sqrs -8.0d0*beta_2_U*sqrs  
     &  *(1.0d0+2.0d0*beta_2_U*rs)-beta_1_U*(5.0d0+21.0d0*beta_2_U*rs))
     &  /(27.0d0*den**(2.0d0))
        fxc=fx+fc

      else
        rsl=dlog(rs)
    
        ex= ax*rh**(.3333333333333333d0)
        ec= A_U*rsl + B_U + C_U*rs*rsl+ D_U*rs
        exc=ex+ec

        vx= 4.d0/3.d0*ex
        vc= A_U*rsl + (B_U - .3333333333333333d0*A_U) 
     &    + (2.0d0/3.0d0)*C_U*rs*rsl 
     &    + (1.0d0/3.0d0)*(2.0d0*D_U-C_U)*rs
        vxc=vx+vc

        fx=4.0d0/9.0d0*ax*rh**(-2.0d0/3.0d0)
        fc=-4.0d0*PI*rs**(4.0d0)*(C_U+2.0d0*D_U 
     &  +3.0d0*A_U/rs+2.0d0*C_U*rsl)/27.0d0
        fxc=fx+fc
 
      end if


  100 return
      end
