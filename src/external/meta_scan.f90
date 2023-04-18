subroutine meta_scan_x(f,f_der,rhoa,rhob,grada,gradb,ta,tb,nfun)

implicit none


! input variables
real(8), intent(in) :: rhoa,rhob                  ! densities (up, dn)
real(8), dimension(3), intent(in) :: grada, gradb ! gradients (up, dn)
real(8), intent(in) :: ta, tb                     ! KE densities (up, dn)
integer, intent(in) :: nfun                       ! number of functional

! output variables
real(8), intent(out) :: f                         ! exch. energy density
real(8), dimension(7), intent(out) :: f_der       ! partials of fa & fb

! local variables
real(8) :: gaa, gbb,gab       ! square gradients for up & dn
real(8) :: fa, fb         ! double exchange energy densities for up & dn
real(8) :: farho, fbrho   ! partials of fa & fb wrt spin densities
real(8) :: fagrad, fbgrad ! partials of fa & fb wrt square gradients
real(8) :: fatau, fbtau   ! partials of fa & fb wrt KE densities

real(8), parameter :: threshold = 1.0d-10

!----------------------!

gaa = grada(1)*grada(1) + grada(2)*grada(2) + grada(3)*grada(3)
gbb = gradb(1)*gradb(1) + gradb(2)*gradb(2) + gradb(3)*gradb(3)
!gab = grada(1)*gradb(1) + grada(2)*gradb(2) + grada(3)*gradb(3)


fa = 0.0d0
fb = 0.0d0
farho = 0.0d0
fbrho = 0.0d0
fagrad = 0.0d0
fbgrad = 0.0d0
fatau = 0.0d0
fbtau = 0.0d0


if ( rhoa >  threshold ) then
    call meta_scan_exch(1,1.d0,rhoa,gaa,ta,fa,farho,fagrad,fatau)
endif


if ( rhob >  threshold ) then
    call meta_scan_exch(1,1.d0,rhob,gbb,tb,fb,fbrho,fbgrad,fbtau)
endif


f = 0.5d0*(fa+fb)
f_der(1) = farho
f_der(2) = fbrho
f_der(3) = fagrad
f_der(4) = fbgrad
f_der(5) = 0.0d0
f_der(6) = fatau
f_der(7) = fbtau

end subroutine meta_scan_x


!=====================================================================!
!=====================================================================!


subroutine meta_scan_exch(nfun,scal,rhoa,gaa,ta,f,frho,fgrad,ftau)

! This subroutine implements the exchange part of the SCAN  
!    (nfun=1)functional.
! It finds the needed exchange energy density and derivatives by spin 
!     scaling. 
!
! Ref. J. Sun, A, Ruzsinszky and J.P. Perdew,
!      Phys. Rev. Lett. 105, 036402 (2015).
! 
!
! Authors: A. Ruzsinszky, I.Y. Zhang, 
! Date: March 2016
! e-mail: aruzsinszky@temple.edu
!         zhang@fhi-berlin.mpg.de

! All expessions independently re-worked as a check on Jianwei Sun's 
!    coding of selfconsistent SCAN in Gaussian and VASP.

! The set of derivatives is the minimal set needed to construct the SCAN !  exchange potential, and may differ from the set of Furche and Perdew. !  Thus compatibility with calling programs should be checked.


implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: nfun    ! functional number
real(r8k), intent(in) :: scal       ! scaling factor
real(r8k), intent(in) :: rhoa       ! half electron density 
real(r8k), intent(in) :: gaa        ! square gradient of the density
real(r8k), intent(in) :: ta         ! kinetic energy density

! output variables
real(r8k), intent(out) :: f         ! energy density epsilon^{scan}_x = epsilon^{unif}*F_x
real(r8k), intent(out) :: frho      ! = df/d(rho)
real(r8k), intent(out) :: fgrad     ! = df/d(grad)
real(r8k), intent(out) :: ftau      ! = df/d(tau)

! local variables
real(r8k) :: rho     ! total electron density
real(r8k) :: grad    ! gradient**2 of the total density
real(r8k) :: tau     ! total kinetic energy density 

real(r8k) :: p
! = s^2 = square of reduced gradient
! dimensionless density gradient s= |\grad n| / [2(3 \pi^2)^{1/3} n^{4/3}]

real(r8k) :: alpha   ! DDSOS
real(r8k) :: dpdrho  ! = dp/d(rho)
real(r8k) :: dpdgrad ! = dp/d(grad)
real(r8k) :: dadrho  ! = d(alpha)/d(rho)
real(r8k) :: dadgrad ! = d(alpha)/d(grad)
real(r8k) :: dadtau  ! = d(alpha)/d(tau)

!real(r8k) :: ds
real(r8k) :: da,da2
real(r8k) :: y1,y2
real(r8k) :: xx      ! function x
real(r8k) :: dy2dp 
real(r8k) :: dy1dp 
real(r8k) :: dxxdp  ! = dx/d(p)
real(r8k) :: dxxda  ! = dx/d(alpha)
real(r8k) :: Fx,Fx1                    ! enhancement factor
real(r8k) :: dFxdp, dFxda 
real(r8k) :: dFxdrho, dFxdgrad, dFxdtau ! derivatives of Fx
real(r8k) :: ldax, dldaxdrho            ! LDA exchange & derivative
real(r8k) :: hx1                        ! part of enhancement factor
real(r8k) :: gx                         ! part of enhancement factor
real(r8k) :: sfx                        ! part of enhancement factor f_x(\alpha)
real(r8k) :: dhx0dp,dhx1dp, dhx1da
real(r8k) :: dy2da,dsfxda
real(r8k) :: dgxdp, dFx1dp, dFx1da
real(r8k) :: dfdrho, dfdgrad, dfdtau 

real(r8k) :: dirac_ep, theta_alpha_one,dirac_alpha_one ! approxition of unit function and dirac fuction with band width epsilon(dirac_ep)

real(r8k) :: alpha_2

! constants
real(r8k), parameter :: A = -0.738558766382d0       ! -(3/4)(3/pi)^(1/3)
real(r8k), parameter :: B =  0.026121173d0        ! 1/[4(3pi**2)**(2/3)] 
real(r8k), parameter :: C = 0.125d0                 ! 1/8
real(r8k), parameter :: D = 2.871233996d0       ! (3/10)(3*pi**2)**(2/3)
real(r8k), parameter :: uter = 0.3333333333333333d0 ! 1/3
real(r8k), parameter :: dter = 0.6666666666666667d0 ! 2/3
real(r8k), parameter :: qter = 1.3333333333333333d0 ! 4/3
real(r8k), parameter :: cter = 1.6666666666666667d0 ! 5/3 
real(r8k), parameter :: oter = 2.6666666666666667d0 ! 8/3
real(r8k), parameter :: unter = 3.666666666666667d0 ! 11/3


! for meta-scan
real(r8k), parameter :: one   = 1.000000000d0
real(r8k), parameter :: two   = 2.000000000d0
real(r8k), parameter :: four  = 4.000000000d0
real(r8k), parameter :: kappa = 0.065000000d0   ! kappa for meta-scan
real(r8k), parameter :: uak   = 0.123456790d0   ! 10/81
real(r8k), parameter :: b2    = 0.120830460d0   ! (5913/405000)^(1/2)
real(r8k), parameter :: b1    = 0.156632077d0   ! (511/13500)/(2*b2)
real(r8k), parameter :: b3    = 0.500000000d0
real(r8k), parameter :: b4    = 0.121831510d0 ! uak^2/k1-1606/18225-b1^2
real(r8k), parameter :: a1    = 4.947900000d0   ! fit Ex[HF] of H atom
real(r8k), parameter :: hx0   = 1.174000000d0   !
real(r8k), parameter :: c1x   = 0.667000000d0
real(r8k), parameter :: c2x   = 0.800000000d0
real(r8k), parameter :: dx    = 1.240000000d0
real(r8k), parameter :: ds    = 1.d0


real(r8k), parameter :: theta_threshold = 1.d-13

Integer(8) :: alpha_phy = -1

dirac_ep = 1.d-2

!------------------------------------------!

! build main variables
rho  = two*rhoa
!if (rho < 1.d-10) rho = 1.d-10
grad   = four*gaa
tau   = two*ta


p     = B*grad/rho**oter

!ds = (0.5d0**(cter)+ 0.5d0**(cter))/2.d0

alpha = (tau - C*grad/rho)/(D*rho**cter)*ds


call theta_a(dirac_ep,alpha-one,theta_alpha_one)
call theta_a(dirac_ep,alpha-one,dirac_alpha_one)
!if(dirac_alpha_one > 100.) then
!    write(*,*) alpha,theta_alpha_one,dirac_alpha_one
!    stop
!end if
dpdrho = (-oter)*B*grad/rho**unter
dpdgrad = B/rho**oter

dadrho = ((-cter)*(1.d0/D)*(tau/(rho**oter))+oter*C*(1.d0/D)*grad/(rho**unter))/ds


dadgrad = -C*(1.d0/D)/rho**oter/ds
dadtau = (1.d0/D)/rho**cter/ds


da = one-alpha
da2 = da*da
y1 = one+(b4*p/uak)*dexp(-b4*p/uak)
y2 = b1*p+b2*da*dexp(-b3*da2)

xx = uak*p*y1 +y2*y2
  
hx1 = one + kappa - kappa/(one + xx/kappa)
! gx
gx = one - dexp(-a1/p**(0.25d0))

! fx
if (alpha .lt. one-theta_threshold) then
    sfx = dexp(-c1x*alpha/da)
elseif (alpha .gt. one+theta_threshold) then
    sfx = -dx* dexp(c2x/da)
else
    sfx = 0.d0
endif


! == build the enhancement factor Fx==
Fx1 =(hx1+sfx*(hx0-hx1))
Fx = (hx1+sfx*(hx0-hx1))*gx


! == LDA exchange ==
ldax = A*(rho**qter)


! == build energy density==
f = ldax*Fx
!f = ldax

! The first derivatives


! Derivative elements of Fx, dhx1dp, dh1xda, dgxdp, dsfxda          

! dhx0dp
dhx0dp =0.0d0           
! dhx1dp
dy2dp = b1
dy1dp = dexp(-(b4)*p/uak)*((b4/uak) -((b4**2)*p/(uak**2)))
dxxdp = uak*(y1 + p*dy1dp)+two*y2*dy2dp
dhx1dp = dxxdp/(one+xx/kappa)**two

! dhx1da

!==================
! BS comment
! Add Date      : 2017.11.15
! modify Date   : 2017.11.15
!------------------
! error
! changed
dy2da = b2*(-one+two*b3*da2)*dexp(-b3*da2)
!------------------
!dy2da = -b2*(one+two*b3*da2)*dexp(-b3*da2)
!==================


dxxda= two*y2*dy2da
dhx1da=dxxda/(one+xx/kappa)**two
                                   
! dsfxda
dsfxda = 0.0d0


  if (alpha .lt. one-theta_threshold) then !<
    dsfxda = -c1x*dexp(-c1x*alpha/da)/da2
  else if (alpha .gt. one+theta_threshold) then !>
    dsfxda = -c2x*dx*dexp(c2x/da)/da2
  else
   dsfxda = 0.d0
  end if

!==================
! debug
!dsfxda = 0.d0

! dgxdp
dgxdp=0.0d0
 if (p .gt. 0.d0) then
  dgxdp = -(a1/4.0d0)*(dexp(-a1*(p**(-0.25d0)))/p**(1.25d0))
 endif
         
! Derivative elements of Fx without gx       
  
! dFx1dp
dFx1dp = dhx1dp+sfx*(-dhx1dp)
! dfx1da
dFx1da = (one-sfx)*dhx1da+dsfxda*(hx0-hx1)

! The full derivatives of the Fx enhancement factor

! dFx/dp
dFxdp = dFx1dp*gx + Fx1*dgxdp
! dFx/da
dFxda = dFx1da*gx
                       
!  Derivatives of Fx       
dFxdrho = dFxdp*dpdrho+dFxda*dadrho
dFxdgrad = dFxdp*dpdgrad+dFxda*dadgrad
dFxdtau = dFxda*dadtau

! Derivatives of the energy density

dldaxdrho = qter*A*(rho**uter)
dfdrho = dldaxdrho*Fx + ldax*dFxdrho
dfdgrad = 2.d0*ldax*dFxdgrad
dfdtau = ldax*dFxdtau


! For a spin-unpolarized density rho=2*rho_sigma 
! with kinetic-energy density tau=2*tau_sigma,
! f is the exchange energy density, which integrates
! over space to the exchange energy,
! dfdrho is the partial derivative of f with respect to rho
! dfdgrad is the partial derivative of f with respect to 
! the square of the density gradient,
! dfdtau is the partial derivative of f with respect to tau. 


! outputs
f = scal*f
frho   = scal*dfdrho
fgrad = scal*dfdgrad
ftau  = scal*dfdtau    
               
                    
end subroutine meta_scan_exch

!=====================================================================!
!=====================================================================!

subroutine meta_scan_c(f,f_der,rhoa,rhob,grada,gradb,ta,tb,nfun)

implicit none

! input variables
real(8), intent(in) :: rhoa, rhob                  ! densities
real(8), dimension(3), intent(in) :: grada, gradb ! density gradients
real(8), intent(in) :: ta, tb ! kinetic energy densities
integer, intent(in) :: nfun   ! number of the functional

! output variables
real(8), intent(out) :: f  ! energy density       
real(8), dimension(7), intent(out) :: f_der  ! functional derivatives

! local variables
real(8) :: gaa, gbb, gab  ! square gradients  
real(8) :: fa, fb         ! partial derivatives of f wrt rhoa & rhob
real(8) :: fgaa,fgab,fgbb ! partial derivative of f wrt square gradient
real(8) :: fta,ftb        ! partial derivative of f wrt tau  
!----------------------!

gaa = grada(1)*grada(1) + grada(2)*grada(2) + grada(3)*grada(3)
gbb = gradb(1)*gradb(1) + gradb(2)*gradb(2) + gradb(3)*gradb(3)
gab = grada(1)*gradb(1) + grada(2)*gradb(2) + grada(3)*gradb(3)
!grq = gaa+2d0*gab+gbb

!tau = ta + tb
!if( ta+tb > (gaa+2.d0*gab+gbb)/8.d0/(rhoa+rhob))  then
!if( dabs((rhoa-rhob)/(rhob+rhoa)) > 1.d-4) then
   call meta_scan_cor(1,1.d0,rhoa,rhob,gaa,gbb,gab,ta,tb,f,&
  &           fa,fb,fgaa,fgab,fgbb,fta,ftb)
!end if

!end if
f_der(1) = fa
f_der(2) = fb
f_der(3) = fgaa
f_der(5) = fgab
f_der(4) = fgbb
f_der(6) = fta
f_der(7) = ftb

end subroutine

!=====================================================================!
!=====================================================================!

subroutine meta_scan_cor(nfun,scal,rhoa,rhob,gaa,gbb,gab,ta,tb,f,&
       &           fa,fb,fgaa,fgab,fgbb,fta,ftb)

! This subroutine implements the correlation part of the meta-GGA SCAN 
!    functional.
!
! Ref. J. Sun, A. Ruzsinszky, J. P. Perdew, 
!      Phys. Rev. Lett. 105, 036402 (2015).

!
! Authors:  A. Ruzsinszky, I.Y. Zhang,
! Date: March 2016
! e-mail: aruzsinszky@temple.edu
!         zhang@fhi-berlin.mpg.de

! All expressions independently reworked as a check on Jianwei Sun's 
!   coding of selfconsistent SCAN in Gaussian and VASP.

! The set of derivatives is the minimal one needed to construct the SCAN
! correlation potential, and may differ from the set of Furche and Perdew.

implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: nfun   ! number of the functional
                                   !   = 1 standard meta-SCAN
                                   !   = 2 none
                                   !   = 3 none
real(r8k), intent(in) :: scal      ! scaling factor
real(r8k), intent(in) :: rhoa      ! alpha spin-density
real(r8k), intent(in) :: rhob      ! beta spin-density
real(r8k), intent(in) :: gaa       ! = |nabla rhoa|^2
real(r8k), intent(in) :: gab       ! = (nabla rhoa).(nabla rhob)
real(r8k), intent(in) :: gbb       ! = |nabla rhob|^2
real(r8k), intent(in) :: ta        ! kinetic energy density,a
real(r8k), intent(in) :: tb        ! kinetic energy density,b

! output variables
real(r8k), intent(out) :: f     ! energy density
real(r8k), intent(out) :: fa    ! df/d(rhoa)
real(r8k), intent(out) :: fb    ! df/d(rhob)
real(r8k), intent(out) :: fgaa  ! df/d(gab)
real(r8k), intent(out) :: fgab  ! df/d(gab)
real(r8k), intent(out) :: fgbb  ! df/d(gbb)
real(r8k), intent(out) :: fta   ! df/d(ta)
real(r8k), intent(out) :: ftb   ! df/d(tb)

! local variables
real(r8k) :: up,dn,dens ! spin-up, spin-down and total density
real(r8k) :: grq,gr      ! square gradient, gradient
real(r8k) :: p         ! gradient square for exchange (not a mistake...)
real(r8k) :: zeta       ! relative spin polarization
real(r8k) :: alpha      !  DDSOS
real(r8k) :: tau
real(r8k) :: rs          
real(r8k) :: dgrdgaa,dgrdgbb,dgrdgab ! derivatives of the gradient
real(r8k) :: dtaudta,dtaudtb ! = d(tau)/d(ta) , d(tau)/d(tb)
real(r8k) :: dpdgr,dadgr
real(r8k) :: ds, ddsdzeta          
real(r8k) :: drhodup, drhoddn  
real(r8k) :: drsdup, drsddn    
real(r8k) :: dzetadup, dzetaddn    
real(r8k) :: dpdup, dpddn    
real(r8k) :: dadup, daddn    
real(r8k) :: drhodgrq,drsdgrq,dpdgaa,dpdgbb,dpdgab,dadgaa,dadgbb,dadgab
real(r8k) :: drhodtau, drsdtau,dzetadtau,dpdta,dpdtb,dadtau,dadta,dadtb   
real(r8k) :: epsc,depscdrs, depscdzeta, depscdp,depscda 
real(r8k) :: depscdup,depscddn,depscdgaa,depscdgbb,depscdgab
real(r8k) :: depscdta,depscdtb



! constants
real(r8k), parameter :: one = 1.d0
real(r8k), parameter :: ethird = 2.666666666666667d0   ! 8/3
real(r8k), parameter :: tfp = 1.6666666666666667d0     ! 5/3
real(r8k), parameter :: B = 0.026121173d0      ! 1/(4*(3pi**2)**(2/3))
real(r8k), parameter :: C = 0.125d0                    ! 1/8
real(r8k), parameter :: D = 2.871233996d0      ! (3/10)(3pi**2)**(2/3)
real(r8k), parameter :: tfpar = 0.3d0                  ! 3/10
real(r8k), parameter :: dter = 0.66666666666666667d0   ! 2/3
real(r8k), parameter :: elthird = 3.6666666666666667d0 ! 11/3
real(r8k), parameter :: pi = 3.1415926535897932d0
real(r8k), parameter :: ex1=0.333333333333333333d0
real(r8k), parameter :: bb = 0.6203505d0             ! (3/(4*pi))**ex1
real(r8k), parameter :: fthird = 1.333333333333333d0   ! 4/3
!--------------------------------------------------------!



! compute relevant quantities

 up = rhoa
 dn = rhob
 if(up.le.1.d-10)up=1.d-10
 if(dn.le.1.d-10)dn=1.d-10
 dens = up+dn
! if(dens.le.1.d-6)dens=1.d-6

! zeta=min(max(zeta,-0.99999999999990d0),0.99999999999990d0)

 zeta = (up - dn)/(up + dn)
!if(dabs(zeta) -1. < 1.d-4 .and. max(up,dn) < 1.d-4) then
!    zeta = 0.d0
!end if

 zeta=min(max(zeta,-0.99999999999990d0),0.99999999999990d0)
! zeta=min(max(zeta,-1.d0),1.d0)
 grq = gaa+2.d0*gab+gbb
 gr = sqrt(grq)
 p = B*gr*gr*(dens**(-ethird))
 tau = ta + tb
! if( tau < 1.d-6) then
!    tau = 1.d-6
!end if
! z = (C*grq*(1.d0/dens))/tau
 ds = ((one + zeta)**(tfp)+ (one - zeta)**(tfp))/2.d0
 alpha = (tau - C*gr*gr*(1.d0/dens))/(D*dens**tfp*ds)
 if(alpha < 0) then
!
!    alpha = 1.
!    return
!!    write(*,*) "warning : need more grids"
!!    call sleep(1)
!
 end if

 rs = bb*dens**(-ex1)
! write(6,*) rs



! The derivatives of the basic quantities in correlation 
 drhodup = 1.d0
 drhoddn = 1.d0
 drsdup = (-ex1)*bb*dens**(-fthird)
 drsddn = drsdup
 dgrdgaa = 1.d0/2.d0/gr
 dgrdgbb = 1.d0/2.d0/gr
 dgrdgab = 1.d0/gr
 dzetadup = 2.d0*dn/dens**2
!==================
! BS comment
! Add Date      : 2017.11.15
! modify Date   : 2017.11.15
!------------------
! error
! changed
 dzetaddn = -2.d0*up/dens**2
!------------------
! dzetaddn = 2.d0*up/dens**2
!==================
 dtaudta = 1.d0
 dtaudtb = 1.d0
 ddsdzeta = tfp*((one + zeta)**dter - (one - zeta)**dter)/2.d0
 dpdup = -(ethird)*B*grq*dens**(-elthird)
 dpddn = dpdup
 dadup = (-tfp)*(1.d0/D)*tau*dens**(-ethird)/ds + &
   &  ethird*C*(1.d0/D)*grq*dens**(-elthird)/ds - & 
   &  (1.d0/D)*tau*dens**(-tfp)*(1.d0/(ds**2))*ddsdzeta*dzetadup &
   &  +C*(1.d0/D)*grq*dens**(-ethird)*(1.d0/(ds**2))*ddsdzeta*dzetadup 
 daddn = (-tfp)*(1.d0/D)*tau*dens**(-ethird)/ds + &
   &  ethird*C*(1.d0/D)*grq*dens**(-elthird)/ds - &
   &  (1.d0/D)*tau*dens**(-tfp)*(1.d0/(ds**2))*ddsdzeta*dzetaddn &
   &  +C*(1.d0/D)*grq*dens**(-ethird)*(1.d0/(ds**2))*ddsdzeta*dzetaddn 
 drhodgrq = 0.d0
 drsdgrq = 0.d0
 dpdgr = 2.d0*B*gr*(dens**(-ethird))
 dpdgaa = dpdgr*dgrdgaa
 dpdgab = dpdgr*dgrdgab
 dpdgbb = dpdgr*dgrdgbb
 dadgr = (-2.d0*C*gr*(1.d0/dens))/(D*dens**tfp*ds)
 dadgaa = dadgr*dgrdgaa
 dadgab = dadgr*dgrdgab
 dadgbb = dadgr*dgrdgbb
 drhodtau = 0.d0
 drsdtau = 0.d0
 dzetadtau = 0.d0
 dpdta = 0.d0
 dpdtb = 0.d0

 dadtau = (1.d0/D)*dens**(-tfp)*(1.d0/ds)
! if( dens < 1.d-4) then
! dadtau = 0.d0
! end if
 dadta = dadtau*dtaudta
 dadtb = dadtau*dtaudtb
 
! -- compute PBE-like part --

  call meta_scan_details(nfun,rs,zeta,p,alpha,epsc, &
  & depscdrs,depscdzeta,depscdp,depscda)
          
depscdup = depscdrs*drsdup + depscdzeta*dzetadup + &
  &  depscdp*dpdup + depscda*dadup
depscddn = depscdrs*drsddn + depscdzeta*dzetaddn + &
  &  depscdp*dpddn + depscda*daddn
!write(*,*) '-----------'
!write(*,*) depscdzeta*dzetaddn,depscdzeta,dzetaddn,depscddn

depscdgaa = depscdp*dpdgaa + depscda*dadgaa
depscdgab = depscdp*dpdgab + depscda*dadgab
depscdgbb = depscdp*dpdgbb + depscda*dadgbb
depscdta = depscda*dadta
depscdtb = depscda*dadtb

! -- compute the final quantities --
f = dens*epsc

fa  = epsc + dens*depscdup
fb  = epsc + dens*depscddn
!write(*,*) '-+----+----'
!write(*,*) tau,grq
!write(*,*) '-+----+----'
!write(*,*) epsc,dens*depscdup,dens*depscddn
fgaa = dens*depscdgaa
fgab = dens*depscdgab
fgbb = dens*depscdgbb
fta  = dens*depscdta
ftb  = dens*depscdtb



! f is the correlation energy density, which integrates over space to the !    correlation energy
! fa is the parial derivative of f with respect to the up-spin density
! fb is the partial derivative of f with respect to the down-spin density
! faa, fab, fbb are the partial derivatives of f with respect to gaa, gab, and gbb
! fta, ftb  partials wrt to the kinetic energy density

f = scal*f
fa  = scal*fa
fb  = scal*fb
fgaa = scal*fgaa
fgab = scal*fgab
fgbb = scal*fgbb
fta = scal*fta
ftb = scal*ftb

end subroutine meta_scan_cor

!=======================================================================!
!                                                                       !
!=======================================================================!

subroutine meta_scan_details(method,rs,zeta,p,alpha,epsc, &
   &                   depscdrs,depscdzeta,depscdp,depscda)

  implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: method
real(r8k), intent(in) :: rs
real(r8k), intent(in) :: zeta
real(r8k), intent(in) :: p
real(r8k), intent(in) :: alpha

! output variables
real(r8k), intent(out) :: epsc
real(r8k), intent(out) :: depscdrs,depscdzeta,depscdp, depscda 

! local variables
real(r8k) :: t, ts, ah, pon
real(r8k) :: phi,phi3,dphidzeta
real(r8k) :: dphidzetaphi !dphidzeta*phi
real(r8k) :: dtsdrs,dtsdzeta, dtsdp
real(r8k) :: A,alpha1,beta1,beta2,beta3,beta4,pp
real(r8k) :: G,dGdrs
real(r8k) :: epscunif0,depscunif0drs
real(r8k) :: epscunif1,depscunif1drs
real(r8k) :: alphac,dalphacdrs
real(r8k) :: ff,dffdzeta
real(r8k) :: epscunif,depscunifdrs,depscunifdzeta,depscunifdp
real(r8k) :: beta,betanum, betaden,dbetanum, dbetaden,dbetadrs
real(r8k) :: H1,dH1drs,dH1dzeta,dH1dp, dH1dts
real(r8k) :: H0,dH0drs,dH0dp
real(r8k) :: w1scan,gscan,Ascan
real(r8k) :: dw1scandrs, dw1scandzeta
real(r8k) :: dgscandrs, dgscandzeta, dgscandts
real(r8k) :: h1core, h0core
real(r8k) :: rshalf, factor1
real(r8k) :: epsc1,depsc1drs,depsc1dzeta, depsc1dp
real(r8k) :: eclda0scan, declda0scandrs
real(r8k) :: epsc0,depsc0drs, depsc0dzeta, depsc0dp
real(r8k) :: w0scan, ginftyscan,dxzeta,gczeta
real(r8k) :: ddxzetadzeta,dgczetadzeta,dw0scandrs,dginftyscandp
real(r8k) :: fc, dfcda

real(r8k) :: dirac_ep, theta_alpha_one,dirac_alpha_one ! approxition of unit function and dirac fuction with band width epsilon(dirac_ep)

! constants
real(r8k), parameter :: pi = 3.1415926535897932d0
real(r8k), parameter :: ex1=0.333333333333333333d0       ! 1/3
real(r8k), parameter :: dter = 0.66666666666666667d0     ! 2/3
real(r8k), parameter :: qter = 1.3333333333333333d0      ! 4/3  
real(r8k), parameter :: one = 1.0d0
real(r8k), parameter :: ethird = 2.6666666666666667d0    ! 8/3
real(r8k), parameter :: tfp = 1.666666666666666667d0     ! 5/3
real(r8k), parameter :: C = 0.125d0                      ! 1/8
real(r8k), parameter :: tfpar = 0.3d0                    ! 3/10
real(r8k), parameter :: elthird = 3.666666666666666667d0 ! 11/3
real(r8k), parameter :: b = 0.079577472               ! (3/(4*pi))**ex1
real(r8k), parameter :: betamb = 0.066725d0              ! for beta(rs)
real(r8k), parameter :: ct = 1.22772285d0                ! relating t, s (3 pi**2 /16)**(1./3.)
real(r8k), parameter :: gamma1 = 0.031091d0              ! for H1
real(r8k), parameter :: afac = 0.10d0                    ! for beta(rs)
real(r8k), parameter :: bfac = 0.1778d0                  ! for beta(rs)
real(r8k), parameter :: fthird = 1.3333333333333333d0    ! 4/3
real(r8k), parameter :: c1c = 0.64d0                     ! interpolation
real(r8k), parameter :: c2c = 1.5d0                      ! interpolation 
real(r8k), parameter :: dc = 0.7d0                       ! interpolation 
real(r8k), parameter :: b1c = 0.0285764d0                ! ec0 parameter
real(r8k), parameter :: b2c = 0.0889d0                   ! ec0 parameter
real(r8k), parameter :: b3c = 0.125541d0                 ! ec0 parameter
real(r8k), parameter :: ckaild = 0.12802585262626d0      ! in ginftyscan

real(r8k), parameter :: theta_threshold = 1.d-13

!-------------------------------------!

dirac_ep = 1.d-2
call theta_a(dirac_ep,alpha-1.d0,theta_alpha_one)
call eta(dirac_ep,alpha-1.d0,dirac_alpha_one)
if(dirac_alpha_one > 100.) then
    write(*,*) alpha,theta_alpha_one,dirac_alpha_one
    stop
end if
! compute relevant quantities for the PBE-like correlation


phi = ((1.d0+zeta)**dter)/2.d0 + ((1.d0-zeta)**dter)/2.d0
!if(phi < 1.d-4) then
!    phi = 1.d-4
!end if
t = ct*sqrt(p)/(sqrt(rs)*phi)
ts = t**2

!==================
! BS comment
! Add Date      : 2017.11.15
! modify Date   : 2017.11.15
!------------------
! error
! changed
if(dabs(1.d0 - zeta) < 1.d-4 ) then

    dphidzeta = ((1.d0+zeta)**(-ex1))/3.d0

elseif(dabs(1.d0 + zeta) < 1.d-4 ) then

    dphidzeta = - ((1.d0- zeta)**(-ex1))/3.d0

else

    dphidzeta = ((1.d0+zeta)**(-ex1))/3.d0 - ((1.d0- zeta)**(-ex1))/3.d0
end if
!------------------
!dphidzeta = ((1.d0+zeta)**(-ex1))/3.d0 - ((1.d0- zeta)**(-ex1))/3.d0
!==================

! -- compute LDA correlation with its first derivatives --
!==================
! BS comment
! Add Date      : 2017.08.02
! modify Date   : 2017.08.03
!------------------
! \epsilon^{LSDA1}_c is pw92 correlation energy
! J. P. Perdew and Y. Wang, Physical Review B 45, 13244 (1992),
! URL http://dx.doi.org/10.1103/physrevb.45.13244.
!==================
pp=1.d0
A=0.031091d0
alpha1=0.21370d0
beta1=7.5957d0
beta2=3.5876d0
beta3=1.6382d0
beta4=0.49294d0

G = -0.2D1*A * dble(1 + alpha1 * rs) * dlog(0.1D1 + 0.1D1 / A /&
  & (beta1*sqrt(dble(rs)) + dble(beta2 * rs) + dble(beta3 * &
  & rs**(0.3D1/0.2D1)) + dble(beta4 * rs ** (pp + 1))) / 0.2D1)
dGdrs = -0.2D1 * A * alpha1 * dlog(0.1D1 + 0.1D1 / A / (beta1 * sqrt(rs) +&
  &      beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + beta4 * &
  &      rs **(pp + 1)) / 0.2D1) + (0.1D1 + alpha1 * rs) / (beta1 * sqrt(rs) +&
  &      beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + &
  &      beta4 * rs ** (pp + 1))** 2 * (beta1 * rs ** (-0.1D1 / 0.2D1) / &
  &      0.2D1 + beta2 + 0.3D1 /0.2D1 * beta3 * sqrt(rs) + &
  &      beta4 * rs ** (pp + 1) * dble(pp + 1) / rs) / (0.1D1 + 0.1D1 / &
  &      A / (beta1 * sqrt(rs) + beta2 * rs + beta3 *&
  &      rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (pp + 1)) / 0.2D1)

epscunif0=G
depscunif0drs=dGdrs


pp=1.d0
A=0.015545d0
alpha1=0.20548d0
beta1=14.1189d0
beta2=6.1977d0
beta3=3.3662d0
beta4=0.62517d0
! =================
G = -0.2D1* A*dble(1+alpha1*rs) * dlog(0.1D1 + 0.1D1 / A / &
  & (beta1 * sqrt(dble(rs)) + dble(beta2 * rs) + &
  & dble(beta3*rs**(0.3D1/0.2D1))+dble(beta4* rs**(pp+1))) / 0.2D1)

dGdrs = -0.2D1*A*alpha1 * dlog(0.1D1 + 0.1D1 / A / (beta1 * &
  &     sqrt(rs) + beta2*rs + beta3*rs **(0.3D1/0.2D1) + &
  &     beta4*rs**(pp+1))/0.2D1) + (0.1D1+alpha1*rs) / &
  &     (beta1*sqrt(rs)+beta2*rs + beta3*rs**(0.3D1/0.2D1) + &
  &     beta4*rs**(pp+1)) ** 2 * (beta1*rs**(-0.1D1/0.2D1) / &
  &     0.2D1 + beta2 + 0.3D1/0.2D1*beta3*sqrt(rs) + beta4* &
  &     rs**(pp+1)* dble(pp+1) / rs) / (0.1D1 + 0.1D1 / A / &
  &     (beta1*sqrt(rs) + beta2*rs + beta3* &
  &     rs**(0.3D1/0.2D1) + beta4*rs**(pp + 1)) / 0.2D1)
epscunif1=G
depscunif1drs=dGdrs

pp=1.d0
A=0.016887d0
alpha1=0.11125d0
beta1=10.357d0
beta2=3.6231d0
beta3=0.88026d0
beta4=0.49671d0
G = -0.2D1*A*dble(1+alpha1*rs) * dlog(0.1D1 + 0.1D1 / A / &
  & (beta1*sqrt(dble(rs)) + dble(beta2*rs) + dble(beta3 * &
  & rs**(0.3D1/0.2D1)) + dble(beta4*rs**(pp+1))) / 0.2D1)
dGdrs = -0.2D1*A*alpha1*dlog(0.1D1 + 0.1D1 / A / (beta1 * &
  &    sqrt(rs) + beta2*rs + beta3*rs**(0.3D1/0.2D1) + beta4 * &
  &    rs**(pp+1))/0.2D1) + (0.1D1+alpha1*rs) / (beta1*sqrt(rs) + &
  &    beta2*rs + beta3*rs**(0.3D1/0.2D1) + beta4* &
  &    rs**(pp+1))**2 * (beta1*rs**(-0.1D1/0.2D1) / 0.2D1 + &
  &    beta2 + 0.3D1/0.2D1 * beta3*sqrt(rs) + &
  &    beta4*rs**(pp+1)*dble(pp+1) / rs) / (0.1D1 + 0.1D1 / A / &
  &    (beta1*sqrt(rs) + beta2*rs + beta3* &
  &    rs**(0.3D1/0.2D1) + beta4*rs**(pp + 1)) / 0.2D1)
alphac=-G
dalphacdrs=-dGdrs

ff = ((1+zeta)**(0.4D1/ 0.3D1) + (1-zeta)**(0.4D1/0.3D1) - 2) / &
  &  (2 * 2 ** (0.1D1 / 0.3D1) - 2)
dffdzeta = (0.4D1/0.3D1 * dble((1 + zeta) ** (0.1D1 / 0.3D1)) - &
  &  0.4D1/ 0.3D1 * dble((1 - zeta) ** (0.1D1 / 0.3D1))) / &
  &  dble(2 * 2 ** (0.1D1 / 0.3D1) - 2)

epscunif = epscunif0+alphac*ff*(1.d0-zeta**4)/1.709921d0+ &
  &         (epscunif1-epscunif0)*ff*zeta**4
depscunifdrs = depscunif0drs+dalphacdrs*ff*(1.d0-zeta**4)/1.709921d0+ &
  &         (depscunif1drs-depscunif0drs)*ff*zeta**4
depscunifdzeta = alphac/1.709921d0*(dffdzeta*(1.d0-zeta**4)- &
  &          ff*4.d0*zeta**3)+(epscunif1-epscunif0)*(dffdzeta*zeta**4+ &
  &          ff*4.d0*zeta**3)
depscunifdp = 0.d0

! -- compute GGA part --

! epsc1 for meta-SCAN
phi3 = phi**3.
pon =-epscunif/(phi3*gamma1)
w1scan = dexp(pon)- one
betanum = one + afac*rs
betaden = one + bfac*rs
beta = betamb*betanum/betaden
Ascan = beta/(gamma1*w1scan)
gscan = (1.0d0+4.0d0*Ascan*ts)**(-0.25d0)
h1core = one + w1scan*(one - gscan)
ah = gamma1*phi3
H1 = ah*dlog(h1core)
epsc1 = epscunif+H1

! The derivatives of epsc1
!******************************************************************
!******************************************************************
!==================
! BS comment
! Add Date      : 2017.11.14
! modify Date   : 2017.11.14
!------------------
! error
! changed
dtsdrs = ct**2.d0*p*(-1.d0/rs**2)*(1.d0/phi**2)
dtsdzeta = ct**2.d0*p*(1.d0/rs)*(-2.d0)*(1.d0/phi**3)*dphidzeta
dtsdp = ct**2.d0*(1.d0/rs)*(1.d0/phi**2)
!-------------------
!dtsdrs = ct*p*(-1.d0/rs**2)*(1.d0/phi**2)
!dtsdzeta = ct*p*(1.d0/rs)*(-2.d0)*(1.d0/phi**3)*dphidzeta
!dtsdp = ct*(1.d0/rs)*(1.d0/phi**2)
!==================


!epsc1 = epscunif + H1

! The rs dependence of beta
dbetadrs = betamb*(betaden*afac-betanum*bfac) &
   &         /betaden**2.
     
dw1scandrs =-(w1scan + one)*depscunifdrs/(gamma1*phi3)
dw1scandzeta =-(w1scan + one)/(gamma1*phi3) &
  &         *(depscunifdzeta-3.d0*epscunif*dphidzeta/phi)
!dw1scandzeta =-(w1scan + one)/(gamma1*phi3) &
!  &         *(depscunifdzeta-3.d0*epscunif*dphidzetaphi/phi**2)

!Ascan = beta(gamma1*w1scan)
        
! gscan = one/(one + 4.0d0*Ascan*ts)**0.250d0

!==================
! BS comment
! Add Date      : 2017.11.14
! modify Date   : 2017.11.14
!------------------
! error
! changed
dgscandrs = -(ts/(one+4.d0*Ascan*ts)**(1.25d0))*(dbetadrs/(gamma1*w1scan)-&
  &          beta*dw1scandrs/(gamma1*w1scan**2))
!-------------------
!dgscandrs = -(ts/(one+Ascan*ts)**(1.25d0))*(dbetadrs/(gamma1*w1scan)-&
!  &          beta*dw1scandrs/(gamma1*w1scan**2))
!==================

!==================
! BS comment
! Add Date      : 2017.11.14
! modify Date   : 2017.11.14
!------------------
! error
! changed
dgscandzeta = (ts/(one+4.d0*Ascan*ts)**(1.25d0))*beta*dw1scandzeta/ &
  &         (gamma1*w1scan**2)
!-------------------
!dgscandzeta = (ts/(one+Ascan*ts)**(1.25d0))*beta*dw1scandzeta/ &
!  &         (gamma1*w1scan**2)
!==================
dgscandts = -Ascan/(one + 4.0d0*Ascan*ts)**(1.250d0)

!==================
! BS comment
! Add Date      : 2017.11.14
! modify Date   : 2017.11.14
!------------------
! error
! changed
dH1drs = (ah/h1core)*(dw1scandrs*(1.d0-gscan) - w1scan*dgscandrs)
!-------------------
!dH1drs = (ah/h1core)*(dw1scandrs*(1.d0-gscan) + w1scan*dgscandrs)
!==================

dH1dzeta = (ah/h1core)*(dw1scandzeta*(1.d0-gscan)-w1scan*dgscandzeta)

dH1dts = (ah/h1core)*(-w1scan*dgscandts)

dH1drs = dH1drs + dH1dts*dtsdrs

!==================
! BS comment
! Add Date      : 2017.11.14
! modify Date   : 2017.11.14
!------------------
! error
! changed
dH1dzeta = dH1dzeta + dH1dts*dtsdzeta + ah*dlog(h1core)*3./phi*dphidzeta
!dH1dzeta = dH1dzeta + dH1dts*dtsdzeta + ah*dlog(h1core)*3./(phi**2)*dphidzetaphi
!-------------------
!dH1dzeta = dH1dzeta + dH1dts*dtsdzeta
!==================
dH1dp = dH1dts*dtsdp

!Output epsc1 and derivatives depsc1drs, depsc1dp,depscdzeta, & depscda 

depsc1drs = depscunifdrs + dH1drs 
depsc1dzeta = depscunifdzeta + dH1dzeta 
depsc1dp = dH1dp  


!*****************************************************************************  
! epscGGA0 
!*****************************************************************************
!*****************************************************************************

! p = grq/(38.28312d0*dens**(ethird))
ginftyscan = 1.0d0/(one+4.0d0*ckaild*p)**(0.25d0)
factor1=(one+b2c*rs**(0.5d0)+b3c*rs)
eclda0scan = -b1c/factor1
w0scan = dexp(-eclda0scan/b1c)-1.0d0
h0core = one + w0scan*(one - ginftyscan)
H0 = b1c*dlog(h0core)
dxzeta = ((one+zeta)**qter+(one-zeta)**qter)/2.0d0
gczeta = (3.3631d0-2.3631d0*dxzeta)*(one-zeta**12)
epsc0 = (eclda0scan + H0)*gczeta

!******************************************************************************
! The derivatives of epscGGA0

rshalf = sqrt(rs)
declda0scandrs =(b3c + 0.5d0*b2c/rshalf)*eclda0scan**2/b1c
      
! gc(zeta) and its derivatives

ddxzetadzeta=2.d0*((one+zeta)**ex1-(one-zeta)**ex1)/3.d0
dgczetadzeta=-2.3631d0*ddxzetadzeta*(one-zeta**12)- &
  & 12.d0*(3.3631d0-2.363d0*dxzeta)*zeta**11

! H0 and its derivatives
dw0scandrs =-(w0scan + one)*declda0scandrs/b1c 
dginftyscandp = -ckaild*(ginftyscan**5)

dH0drs = b1c*dw0scandrs*(one-ginftyscan)/h0core 
dH0dp = -b1c*w0scan*dginftyscandp/h0core


if (alpha .lt. 1.0d0-theta_threshold) then
    fc=dexp(-c1c*alpha/(one - alpha))
elseif (alpha .gt. 1.0d0 +theta_threshold) then
    fc=-dc*dexp(c2c/(one - alpha))
else
    fc=0.0d0
endif
!==================


!  The derivatives of the fc interpolation function



if (alpha .lt. one-theta_threshold) then
     dfcda = -c1c*fc/(alpha - one)**2
elseif(alpha .gt. one+theta_threshold) then
     dfcda = c2c*fc/(alpha - one)**2
else
     dfcda = 0.0d0
endif
!==================


! meta-SCAN correlation
epsc = epsc1+fc*(epsc0-epsc1)

! Output epsc0 AND ITS DERIVATIVES depsc0drs, depsc0dp,depsc0dzeta, and 
!     depsc0da 

depsc0drs = (declda0scandrs + dH0drs)*gczeta
depsc0dzeta = (eclda0scan+H0)*dgczetadzeta
depsc0dp = dH0dp*gczeta  

!*****************************************************************
! The derivatives of the corr. energy density wrt rs, zeta, p and alpha

 depscdrs = depsc1drs + fc*(depsc0drs - depsc1drs)
 depscdzeta = depsc1dzeta + fc*(depsc0dzeta -depsc1dzeta)
 depscdp = depsc1dp + fc*(depsc0dp -depsc1dp)
 depscda = dfcda*(epsc0-epsc1)
 
! -- compute final values --

end subroutine meta_scan_details

subroutine eta(ep,x,y)
!================
! auther : BS
! aim : to approximate Diract function delta(x)
! limit_{ep -> 0+} eta(x) = delta(x)
! we use Poisson function
! eta(x) = 1/\pi \frac{\epsilon}{\epsilon^2 + x^2}
!--------------------
! in :
! ep real(r8k) the width of eta function
! x  real(r8k) input position
! out :
! y  real(r8k) return the Poisson function value
!================

 integer, parameter :: r8k = selected_real_kind(15,307)

 ! in
 real(r8k),intent(in) ::  ep
 real(r8k),intent(in) ::  x

 !out
 real(r8k),intent(out) :: y


  if(ep <= 0.d0) then
    y = 0.d0
  else

    y = 1.d0/(3.1415926535897932d0) *(ep)/(ep**2.d0+x**2.d0)

  end if


end subroutine

subroutine theta_a(ep,x,y)
!================
! auther : BS
! aim : to approximate unit function \theta(x) = \int delta(x)
! limit_{ep -> 0+} theta_a(x) = theta(x)
! theta(0+) = 1,theta(0-) = 0,theta(0) = 0.5
!
! we use the integration of Poisson function
! theeta_a(x) = 1/\pi * actan(x/\epsilon)
!--------------------
! in :
! ep real(r8k) the width of eta function
! x  real(r8k) input position
! out :
! y  real(r8k) return the function value
!================
 integer, parameter :: r8k = selected_real_kind(15,307)
 ! in
 real(r8k),intent(in) ::  ep
 real(r8k),intent(in) ::  x

 !out
 real(r8k),intent(out) :: y

 real(r8k) :: exp_arg
 real(r8k) :: max_arg

  if(ep <= 0.d0) then
    y = 0.d0
  else
!    write(*,*) x,ep,x/ep,atan(x/ep)
!    y = 1.d0/(3.1415926535897932d0)* atan(x/ep)+0.5d0

    ! Must check for overflow
    exp_arg = -2.d0*x/ep
    max_arg = maxexponent(x)*log(2.d0)

    if (exp_arg < max_arg) then
      y = 1.d0/(1.d0+exp(-2.d0*x/ep))
    else
      y = 0.d0
    end if

  end if

end subroutine

