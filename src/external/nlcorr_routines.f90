module nlcorr_routines

!Routines written by Andris Gulans, Department of Applied Physics, COMP, Aalto University, Espoo, Finland. agl@cc.hut.fi 


implicit none


contains


  double precision function qq00(n,nabla)
    implicit none
    double precision n,nabla
    double precision kf,tt
    double precision pi,Zab,alpha
    parameter (pi=3.141592653589793116d0,Zab=-0.8491d0,alpha=4d0*pi/3d0)
    !double precision EcLDA
    !external EcLDA
    kf=(3d0*pi*pi*n)**0.333333333333d0
    tt=nabla**2
    qq00=kf*(1d0-Zab*tt/((6d0*n*kf)**2))-alpha*EcLDA(kf)
  end function qq00
  
  
  
  
  double precision function EcLDA(kf)
    implicit none
    double precision kf,rs
    double precision pi,thr
    parameter (pi=3.1415926535898d0,thr=1.9191582926775d0)
    double precision c0,c1,b1,b2,b3,b4,a1
    parameter (c0=0.031091d0,c1=0.046644d0,b3=1.6382d0,b4=0.49294d0,a1=0.21370d0)
    rs=thr/kf 
    b1=exp(-c1/(2*c0))/(2*c0)
    b2=2*c0*b1*b1
    EcLDA=-2*c0*(1+a1*rs)*log(1+1/(2*c0*(b1*sqrt(rs)+b2*rs+b3*sqrt(rs*rs*rs)+b4*rs*rs)))
  end function EcLDA
  
  
  
  
  subroutine Calc_q0_full(n,nabla,q0,dqdn,dqdg)
    !  use prec
    implicit none
    real*8 n,nabla(3),q0,dqdn,dqdg
    real*8 kf,tt
    real*8 pi,Zab,alpha,EcLDA,Egrad,rs,dEcLDA,denom,denom2,rs12
    parameter (pi=3.141592653589793116d0,Zab=-0.8491d0,alpha=4d0*pi/3d0)
    real*8 thr,c0,c1,b1,b2,b3,b4,a1
    parameter (c0=0.031091d0,c1=0.046644d0,b3=1.6382d0,b4=0.49294d0,a1=0.21370d0)
    parameter (thr=1.9191582926775d0)
    
    kf=(3d0*pi*pi*n)**0.333333333333d0
    tt=nabla(1)**2+nabla(2)**2+nabla(3)**2
    
    rs=thr/kf
    rs12=sqrt(rs)
    
    b1=exp(-c1/(2*c0))/(2*c0)
    b2=2*c0*b1*b1
    
    denom=2d0*c0*(b1*rs12+b2*rs+b3*rs*rs12+b4*rs*rs)
    denom2=2d0*c0*(0.5d0*b1/rs12+b2+1.5d0*b3*rs12+2d0*b4*rs)
    EcLDA=-2*c0*(1+a1*rs)*log(1d0+1d0/denom)
    dEcLDA=a1*EcLDA/(1d0+a1*rs)+2*c0*(1d0+a1*rs)/(denom*(denom+1))*denom2
    denom=Zab*tt/((6d0*n*kf)**2)
    Egrad=kf*(1d0-denom)
    
    !      answer=kf*(1_q-Zab*tt/((6_q*n*kf)**2))-alpha*EcLDA
    q0=Egrad-alpha*Eclda
    dqdn=(Egrad+8d0*kf*denom+alpha*dEcLDA*rs)/(3d0)
    !     dqdg=-2d0*kf*Zab*sqrt(tt)/((6d0*n*kf)**2)*n  !original dqdg: w/ respect to grad only.
    dqdg=-kf*Zab/((6d0*n*kf)**2)*n  !new
    
    
    
    !      q0=rs
    !      dqdn=-rs/(3_q*n)
  end subroutine Calc_q0_full
  
  








    end module nlcorr_routines
