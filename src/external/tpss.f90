subroutine tpssx(f,f_der,rhoa,rhob,grada,gradb,ta,tb,nfun)

use constants, only : NUM_ZERO

implicit none

! input variables
real(8), intent(in) :: rhoa,rhob  ! densities
real(8), dimension(3), intent(in) :: grada, gradb ! density gradients
real(8), intent(in) :: ta, tb ! kinetic energy densities
integer, intent(in) :: nfun  ! number of the functional

! output variables
real(8), intent(out) :: f  ! energy density
real(8), dimension(7), intent(out) :: f_der  ! functional derivatives

! local variables
real(8) :: gaa, gbb  ! square gradients  
real(8) :: f1, f2    ! temporary energy densities
! real(8) :: dummy     ! dummy variable

!----------------------!

gaa = grada(1)*grada(1) + grada(2)*grada(2) + grada(3)*grada(3)
gbb = gradb(1)*gradb(1) + gradb(2)*gradb(2) + gradb(3)*gradb(3)


if (nfun .eq. 1) then  ! TPSS
 
   if (rhoa.gt.NUM_ZERO) then
     call tpss_exch(1,1.d0,rhoa,gaa,ta,f1,f_der(1),f_der(3),f_der(6))
   else
     f1 = 0.d0
     f_der(1) = 0.d0
     f_der(3) = 0.d0 
     f_der(6) = 0.d0
   end if 

   if (rhob.gt.NUM_ZERO) then
     call tpss_exch(1,1.d0,rhob,gbb,tb,f2,f_der(2),f_der(4),f_der(7))
   else
     f2 = 0.d0
     f_der(2) = 0.d0
     f_der(4) = 0.d0 
     f_der(7) = 0.d0
   end if 

   f = 0.5d0*(f1+f2)

end if

if (nfun.eq.2 .or. nfun.eq.3) then  ! revTPSS (also TPSSloc)

   if (rhoa.gt.NUM_ZERO) then
     call tpss_exch(2,1.d0,rhoa,gaa,ta,f1,f_der(1),f_der(3),f_der(6))
   else
     f1 = 0.d0
     f_der(1) = 0.d0
     f_der(3) = 0.d0 
     f_der(6) = 0.d0
   end if 

   if (rhob.gt.NUM_ZERO) then
     call tpss_exch(2,1.d0,rhob,gbb,tb,f2,f_der(2),f_der(4),f_der(7))
   else
     f2 = 0.d0
     f_der(2) = 0.d0
     f_der(4) = 0.d0 
     f_der(7) = 0.d0
   end if 

   f = 0.5d0*(f1+f2)

end if

! f_der = 0.d0 ! at present derivatives are not used!!

end subroutine tpssx

!=====================================================================!
!=====================================================================!

subroutine tpss_exch(nfun,scal,rho,gaa,ta,f,fa,fgaa,fta)

! This subroutine implements the exchange part of the TPPS [1] (nfun=1) and
! revTPSS [2] (nfun=2) functionals
!
! [1] J. Tao, J. P. Perdew, V. N. Staroverov, G. E. Scuseria,
!                                   Phys. Rev. Lett. 91, 146401 (2003).
! [2] J. P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin,  J. Sun,
!                                   Phys. Rev. Lett. 103, 026403 (2009).
!
! Author: E. Fabiano, L.A. Constantin
! Date: June 2012
! e-mail: eduardo.fabiano@nano.cnr.it
!         Lucian.Constantin@iit.it
!

implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: nfun   ! functional number
real(r8k), intent(in) :: scal      ! scaling factor
real(r8k), intent(in) :: rho       ! electron density (half of total density)
real(r8k), intent(in) :: gaa       ! square gradient of the density
real(r8k), intent(in) :: ta        ! kinetic energy density

! output variables
real(r8k), intent(out) :: f        ! energy density
real(r8k), intent(out) :: fa       ! = df/d(rho)
real(r8k), intent(out) :: fgaa     ! = df/d(gaa)
real(r8k), intent(out) :: fta      ! = df/d(tau)

! local variables
real(r8k) :: dens   ! total electron density
real(r8k) :: gra    ! gradeint of the total density
real(r8k) :: tau    ! kinetic energy density = 2*ta
real(r8k) :: p      ! = s^2 = |nabla dens|^2/[((2(3pi^2)^(1/3))^2)*dens**oter]
real(r8k) :: dpdrho ! = dp/d(rho)
real(r8k) :: dpdgaa ! = dp/d(gaa)
real(r8k) :: tw     ! von weizsaecker kinetic energy
real(r8k) :: dtwdrho ! = d(tw)/d(rho)
real(r8k) :: dtwdgaa ! = d(tw)/d(gaa)
real(r8k) :: z      ! = tw/tau
real(r8k) :: dzdrho ! = dz/d(rho)
real(r8k) :: dzdgaa ! =dz/d(gaa)
real(r8k) :: dzdtau ! =dz/d(tau)
real(r8k) :: alpha  ! function alpha
real(r8k) :: dalphadrho ! =d(alpha)/d(rho)
real(r8k) :: dalphadgaa ! =d(alpha)/d(gaa)
real(r8k) :: dalphadtau ! =d(alpha)/d(tau)
real(r8k) :: qb         ! function qb
real(r8k) :: dqbdrho    ! = d(qb)/d(rho)
real(r8k) :: dqbdgaa ! =d(qb)/d(gaa)
real(r8k) :: dqbdtau ! =d(qb)/d(tau)
real(r8k) :: tmp1,tmp2,tmp3,tmp4 ! scratch variables
real(r8k) :: y1,y2,y3,y4,y5,y6,y7
real(r8k) :: dy1drho,dy1dgaa,dy1dtau
real(r8k) :: dy2drho,dy2dgaa,dy2dtau
real(r8k) :: dy3drho,dy3dgaa,dy3dtau
real(r8k) :: dy4drho,dy4dgaa,dy4dtau
real(r8k) :: dy5drho,dy5dgaa,dy5dtau
real(r8k) :: dy6drho,dy6dgaa,dy6dtau
real(r8k) :: dy7drho,dy7dgaa,dy7dtau
real(r8k) :: xx    ! function x
real(r8k) :: dxxdrho  ! = dx/d(rho)
real(r8k) :: dxxdgaa  ! = dx/d(gaa)
real(r8k) :: dxxdtau  ! = dx/d(tau)
real(r8k) :: fx   ! enhancement factor
real(r8k) :: dfxdrho, dfxdgaa, dfxdtau ! derivatives of fx
real(r8k) :: ldax, dldaxdrho ! LDA exchange and its derivative
real(r8k) :: cc   ! constant c
real(r8k) :: ee   ! constant e
real(r8k) :: mu   ! mu

! constants
real(r8k), parameter :: A = -0.738558766382d0       ! -(3/4)(3/pi)^(1/3)
real(r8k), parameter :: B =  38.2831200025d0        ! (2(3pi^2)^(1/3))^2
real(r8k), parameter :: uott = 0.125d0 ! = 1/8
real(r8k), parameter :: uter = 0.3333333333333333d0 ! 1/3
real(r8k), parameter :: dter = 0.6666666666666667d0 ! 2/3
real(r8k), parameter :: qter = 1.3333333333333333d0 ! 4/3 
real(r8k), parameter :: oter = 2.6666666666666667d0 ! = 8/3
real(r8k), parameter :: unter = 3.666666666666667d0 ! 11/3
real(r8k), parameter :: tmez = 1.5d0  ! = 3/2
real(r8k), parameter :: c10_81 = 0.12345679012345679012d0 ! = 10/81
real(r8k), parameter :: c146_2025 = 0.07209876543209876543d0 ! = 146/2025
real(r8k), parameter :: c73_405 = 0.18024691358024691358d0 ! = 73/405
real(r8k), parameter :: bb = 0.40d0  ! constant b
real(r8k), parameter :: kappa = 0.804d0 ! kappa

!------------------------------------------!

! build main variables and their derivatives
dens = 2.d0*rho
gra = 4.d0*gaa
tau = 2.d0*ta
p = gra/(B*dens**oter)
dpdrho = -oter*gra/(B*dens**unter)
dpdgaa = 1.d0/(B*dens**oter)
tw = uott*gra/dens
dtwdrho = -tw/dens
dtwdgaa = uott/dens
z = tw/tau
dzdrho = dtwdrho/tau
dzdgaa = dtwdgaa/tau
dzdtau = -z/tau

! build function qb and its derivatives
alpha = (5.d0*p/3.d0)*((1.d0/z)-1.d0)
tmp1 = -((5.d0*p/3.d0)/(z*z))
dalphadrho = (5.d0*dpdrho/3.d0)*((1.d0/z)-1.d0) + tmp1*dzdrho
dalphadgaa = (5.d0*dpdgaa/3.d0)*((1.d0/z)-1.d0) + tmp1*dzdgaa
dalphadtau = tmp1*dzdtau 
qb = 0.45d0*(alpha-1)/dsqrt(1.d0+bb*alpha*(alpha-1.d0)) + dter*p
tmp1 = 0.225d0*(bb*(alpha-1.d0)+2.d0)/((bb*(alpha-1.d0)*alpha+1.d0)**tmez)
dqbdrho = tmp1*dalphadrho + dter*dpdrho
dqbdgaa = tmp1*dalphadgaa + dter*dpdgaa
dqbdtau = tmp1*dalphadtau


! == build x (y1+y2+y3+y4+y5+y6)/y7 ==

if (nfun .eq. 1) then
   cc = 1.59096d0  
   ee = 1.537d0 
   mu = 0.21951d0 
else if (nfun .eq. 2) then
   cc = 2.35204d0
   ee = 2.1677d0
   mu = 0.14d0
end if


! y1
if (nfun .eq. 1) then
   tmp1 = 1+z*z
   tmp2 = c10_81 + cc*z*z/(tmp1*tmp1)
   y1 = tmp2*p
   tmp3 = -2.d0*cc*z*(z*z-1.d0)/(tmp1*tmp1*tmp1)
   dy1drho = tmp3*dzdrho*p + tmp2*dpdrho
   dy1dgaa = tmp3*dzdgaa*p + tmp2*dpdgaa
   dy1dtau = tmp3*dzdtau*p
else if (nfun .eq. 2) then
   tmp1 = 1+z*z
   tmp2 = c10_81 + cc*z*z*z/(tmp1*tmp1)
   y1 = tmp2*p
   tmp3 = -cc*z*z*(z*z-3.d0)/(tmp1*tmp1*tmp1)
   dy1drho = tmp3*dzdrho*p + tmp2*dpdrho
   dy1dgaa = tmp3*dzdgaa*p + tmp2*dpdgaa
   dy1dtau = tmp3*dzdtau*p
end if


! y2
y2 = c146_2025*qb*qb
dy2drho = 2.d0*c146_2025*qb*dqbdrho
dy2dgaa = 2.d0*c146_2025*qb*dqbdgaa
dy2dtau = 2.d0*c146_2025*qb*dqbdtau

! y3
tmp1 = 0.5d0*(0.6d0*z)*(0.6d0*z)+0.5d0*p*p
y3 = -c73_405*qb*dsqrt(tmp1)
tmp2 = -c73_405*dsqrt(tmp1)
tmp3 = -c73_405*qb*0.5d0*p/dsqrt(tmp1)
tmp4 = -c73_405*qb*0.18d0*z/dsqrt(tmp1)
dy3drho = tmp2*dqbdrho + tmp3*dpdrho + tmp4*dzdrho
dy3dgaa = tmp2*dqbdgaa + tmp3*dpdgaa + tmp4*dzdgaa
dy3dtau = tmp2*dqbdtau + tmp4*dzdtau

! y4
y4 = c10_81*c10_81*p*p/kappa
dy4drho = 2.d0*c10_81*c10_81*p*dpdrho/kappa
dy4dgaa = 2.d0*c10_81*c10_81*p*dpdgaa/kappa

! y5
y5 = 2.d0*dsqrt(ee)*c10_81*0.36d0*z*z
dy5drho = 4.d0*dsqrt(ee)*c10_81*0.36d0*z*dzdrho
dy5dgaa = 4.d0*dsqrt(ee)*c10_81*0.36d0*z*dzdgaa
dy5dtau = 4.d0*dsqrt(ee)*c10_81*0.36d0*z*dzdtau

! y6
y6 = ee*mu*p*p*p
dy6drho = 3.d0*ee*mu*p*p*dpdrho
dy6dgaa = 3.d0*ee*mu*p*p*dpdgaa

! y7
tmp1 = 1.d0 +dsqrt(ee)*p
y7 = tmp1*tmp1
dy7drho = 2.d0*tmp1*dsqrt(ee)*dpdrho
dy7dgaa = 2.d0*tmp1*dsqrt(ee)*dpdgaa

!! x
tmp1 = y1 + y2 + y3 + y4 + y5 + y6
xx = tmp1/y7
tmp2 = dy1drho + dy2drho + dy3drho + dy4drho + dy5drho + dy6drho
dxxdrho = (tmp2*y7 - dy7drho*tmp1)/(y7*y7)
tmp2 = dy1dgaa + dy2dgaa + dy3dgaa + dy4dgaa + dy5dgaa + dy6dgaa
dxxdgaa = (tmp2*y7 - dy7dgaa*tmp1)/(y7*y7)
tmp2 = dy1dtau + dy2dtau + dy3dtau + dy5dtau
dxxdtau = tmp2*y7/(y7*y7)


! == build the enhancement factor Fx and its derivatives ==
fx = 1.d0 + kappa - kappa/(1.d0 + xx/kappa)
tmp1 = kappa*kappa/((kappa+xx)*(kappa+xx))
dfxdrho = tmp1*dxxdrho
dfxdgaa = tmp1*dxxdgaa
dfxdtau = tmp1*dxxdtau


! == LDA exchange ==
ldax = A*(dens**qter)
dldaxdrho = qter*A*(dens**uter)


! == build energy density and derivatives ==
f = ldax*fx
fa = dldaxdrho*fx + ldax*dfxdrho
fgaa = 2.d0*ldax*dfxdgaa
fta = ldax*dfxdtau
f = scal*f
fa = scal*fa
fgaa = scal*fgaa
fta = scal*fta


end subroutine tpss_exch


!=====================================================================!
!=====================================================================!

subroutine tpssc(f,f_der,rhoa,rhob,grada,gradb,ta,tb,nfun)

implicit none

! input variables
real(8), intent(in) :: rhoa,rhob  ! densities
real(8), dimension(3), intent(in) :: grada, gradb ! density gradients
real(8), intent(in) :: ta, tb ! kinetic energy densities
integer, intent(in) :: nfun  ! number of the functional

! output variables
real(8), intent(out) :: f  ! energy density
real(8), dimension(7), intent(out) :: f_der  ! functional derivatives

! local variables
real(8) :: gaa, gbb, gab  ! square gradients  
!real(8) :: dummy     ! dummy variable

!----------------------!

gaa = grada(1)*grada(1) + grada(2)*grada(2) + grada(3)*grada(3)
gbb = gradb(1)*gradb(1) + gradb(2)*gradb(2) + gradb(3)*gradb(3)
gab = grada(1)*gradb(1) + grada(2)*gradb(2) + grada(3)*gradb(3)

if (nfun .eq. 1) then  ! TPSS
   call tpss_cor(1,1.d0,rhoa,rhob,gaa,gab,gbb,ta,tb,f,&
     &           f_der(1),f_der(2),f_der(3),f_der(5),f_der(4),f_der(6),f_der(7))
end if

if (nfun .eq. 2) then  ! revTPSS
   call tpss_cor(2,1.d0,rhoa,rhob,gaa,gab,gbb,ta,tb,f,&
     &           f_der(1),f_der(2),f_der(3),f_der(5),f_der(4),f_der(6),f_der(7))
end if

if (nfun .eq. 3) then  ! TPSSloc
   call tpss_cor(3,1.d0,rhoa,rhob,gaa,gab,gbb,ta,tb,f,&
     &           f_der(1),f_der(2),f_der(3),f_der(5),f_der(4),f_der(6),f_der(7))
end if

!f_der = 0.d0 ! at present derivatives are not used!!

end subroutine

!=====================================================================!
!=====================================================================!

subroutine tpss_cor(nfun,scal,rhoa,rhob,gaa,gab,gbb,ta,tb,f,&
     &           fa,fb,fgaa,fgab,fgbb,fta,ftb)

! This subroutine implements the correlation part of the TPSS-like functionals
! TPSS [1] (nfun=1), revTPSS [2] (nfun=2), and TPSSloc [3] (nfun=3).
!
! [1] J. Tao, J. P. Perdew, V. N. Staroverov, G. E. Scuseria,
!                                   Phys. Rev. Lett. 91, 146401 (2003).
! [2] J. P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin,  J. Sun,
!                                   Phys. Rev. Lett. 103, 026403 (2009).
! [3] L. A. Constantin, E. Fabiano, F. Della Sala, (2012) submitted.
!
!
! Author: E. Fabiano, L.A. Constantin
! Date: June 2012
! e-mail: eduardo.fabiano@nano.cnr.it
!         Lucian.Constantin@iit.it
!

  use constants, only : NUM_ZERO

 implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: nfun   ! number of the functional
                                   !   = 1 TPSSc
                                   !   = 2 revTPSSc
                                   !   = 3 TPSSloc
real(r8k), intent(in) :: scal   ! scaling factor
real(r8k), intent(in) :: rhoa   ! alpha spin-density
real(r8k), intent(in) :: rhob   ! beta spin-density
real(r8k), intent(in) :: gaa    ! = |nabla rhoa|^2
real(r8k), intent(in) :: gab    ! = (nabla rhoa).(nabla rhob)
real(r8k), intent(in) :: gbb    ! = |nabla rhob|^2
real(r8k), intent(in) :: ta     ! alpha kinetic energy density
real(r8k), intent(in) :: tb     ! beta kinetic energy density

! output variables
real(r8k), intent(out) :: f    ! energy density
real(r8k), intent(out) :: fa   ! df/d(rhoa)
real(r8k), intent(out) :: fb   ! df/d(rhob)
real(r8k), intent(out) :: fgaa ! df/d(gaa)
real(r8k), intent(out) :: fgab ! df/d(gab)
real(r8k), intent(out) :: fgbb ! df/d(gbb)
real(r8k), intent(out) :: fta  ! df/d(ta)
real(r8k), intent(out) :: ftb  ! df/d(tb)

! local variables
real(r8k) :: up,dn,dens ! spin-up, spin-down and total density
real(r8k) :: grq, gr    ! square gradient and gradient norm
real(r8k) :: dgrdgaa,dgrdgbb,dgrdgab ! derivatives of the gradient
real(r8k) :: zeta       ! relative spin polarization
real(r8k) :: dzetadup,dzetaddn ! = d(zeta)/d(rhoa) , d(zeta)/d(rhob)
real(r8k) :: tau        ! total kinetic energy density
real(r8k) :: dtaudta,dtaudtb ! = d(tau)/d(ta) , d(tau)/d(tb)
real(r8k) :: tw         ! Von Weiszaecker kinetic energy
real(r8k) :: z          ! = tw/tau
real(r8k) :: dzdup,dzddn
real(r8k) :: dzdgr,dzdgaa,dzdgab,dzdgbb
real(r8k) :: dzdtau,dzdta,dzdtb
real(r8k) :: x1
real(r8k) :: eta,detadup,detaddn
real(r8k) :: detadgaa,detadgab,detadgbb
real(r8k) :: grzeta
real(r8k) :: dx1dup,dx1ddn,dx1dgaa,dx1dgab,dx1dgbb
real(r8k) :: dgetadup,dgetaddn,dgetadgaa,dgetadgab,dgetadgbb
real(r8k) :: epscGGA,depscGGAdup,depscGGAddn,depscGGAdgaa
real(r8k) :: depscGGAdgbb,depscGGAdgab,depscGGAdz
real(r8k) :: epsGGA2,depsGGAdn2,depsGGAdgaa2,depsGGAdz2
real(r8k) :: epsGGA3,depsGGAdn3,depsGGAdgaa3,depsGGAdz3
real(r8k) :: epscUPt,depscUPtdup,depscUPtddn
real(r8k) :: depscUPtdgaa,depscUPtdgbb,depscUPtdgab
real(r8k) :: depscUPtdz
real(r8k) :: epscDNt,depscDNtdup,depscDNtddn
real(r8k) :: depscDNtdgaa,depscDNtdgbb,depscDNtdgab
real(r8k) :: depscDNtdz
real(r8k) :: m,n,d,c1,c2,c3,c4
real(r8k) :: cz,dczdzeta,dczdeta,dczdz
real(r8k) :: cg
real(r8k) :: dczdup,dczddn
real(r8k) :: dczdgaa,dczdgab,dczdgbb
real(r8k) :: dczdta,dczdtb
real(r8k) :: y1,dy1dup,dy1ddn,dy1dgaa,dy1dgab,dy1dgbb,dy1dta,dy1dtb
real(r8k) :: y2,dy2dup,dy2ddn,dy2dgaa,dy2dgab,dy2dgbb,dy2dta,dy2dtb
real(r8k) :: y3,dy3dup,dy3ddn,dy3dgaa,dy3dgab,dy3dgbb,dy3dta,dy3dtb
real(r8k) :: y4,dy4dup,dy4ddn,dy4dgaa,dy4dgab,dy4dgbb,dy4dta,dy4dtb
real(r8k) :: y5,dy5dup,dy5ddn,dy5dgaa,dy5dgab,dy5dgbb,dy5dta,dy5dtb
real(r8k) :: epscREV,depscREVdup,depscREVddn
real(r8k) :: depscREVdgaa,depscREVdgab,depscREVdgbb
real(r8k) :: depscREVdta,depscREVdtb
real(r8k) :: epscMGGA,depscMGGAdup,depscMGGAddn
real(r8k) :: depscMGGAdgaa,depscMGGAdgbb,depscMGGAdgab
real(r8k) :: depscMGGAdta,depscMGGAdtb

real(r8k) :: ee1,ee2,ee3,ee4,ee5,ee6

real(r8k) :: threshold

! constants
real(r8k), parameter :: pi = 3.1415926535897932d0
real(r8k), parameter :: ex1=0.333333333333333333d0

!--------------------------------------------------------!

! VB: threshold determines the density in one spin channel below which we
!     expect a floating exception (division by zero). 1.d-10 works but is
!     unfortunately quite high. 1.d-20 did not work in my test using ifort.
  threshold = 1.d-16

! compute relevant quantities
 up=rhoa
 dn=rhob
 dens=up+dn
 ! AJL: Original threshold here was 1.d-18
 ! if(dens.le.1.d-18)dens=1.d-18
 ! if(up.le.1.d-18)up=1.d-18
 ! if(dn.le.1.d-18)dn=1.d-18
 ! Increased so that we can handle errors where one channel is empty (e.g. H)
 ! if (up.lt.threshold.or.dn.lt.threshold) return
 if(dens.le.threshold)dens=threshold
 if(up.le.threshold)up=threshold
 if(dn.le.threshold)dn=threshold
 grq = gaa+2d0*gab+gbb
 gr=dsqrt(grq)
 dgrdgaa=1.d0/2.d0/gr
 dgrdgbb=1.d0/2.d0/gr
 dgrdgab=1.d0/gr
 zeta=(up-dn)/dens
 dzetadup=0.1D1/(up + dn)-(up - dn)/(up + dn)**2
 dzetaddn=-0.1D1/(up + dn)-(up - dn)/(up + dn)**2
 tau=ta+tb
 dtaudta=1.d0
 dtaudtb=1.d0
 tw=gr*gr/8.d0/dens
 z=tw/tau
 if (z.ge.1.d0) z = 1.d0
 dzdup=-z/dens
 dzddn=-z/dens
 dzdgr=2.d0*gr/8.d0/dens/tau
 dzdgaa=dzdgr*dgrdgaa
 dzdgbb=dzdgr*dgrdgbb
 dzdgab=dzdgr*dgrdgab
 dzdtau=-z/tau
 dzdta=dzdtau*dtaudta
 dzdtb=dzdtau*dtaudtb

 x1=dsqrt(dabs(dn*dn*gaa+up*up*gbb-2.d0*up*dn*gab))
 if(x1.le.1.d-98)then
    eta=0.d0
    detadup=0.d0
    detaddn=0.d0
    detadgaa=0.d0
    detadgbb=0.d0
    detadgab=0.d0
 else
    grzeta=2.d0/dens/dens*x1
    eta=grzeta/2.d0/((3.d0*pi*pi*dens)**ex1)
    dx1dup=1.d0/2.d0/x1*(2.d0*up*gbb-2.d0*dn*gab)
    dx1ddn=1.d0/2.d0/x1*(2.d0*dn*gaa-2.d0*up*gab)
    dx1dgaa=1.d0/2.d0/x1*(dn*dn)
    dx1dgbb=1.d0/2.d0/x1*(up*up)
    dx1dgab=-1.d0/2.d0/x1*(2.d0*up*dn)
    dgetadup=2.d0/dens/dens*dx1dup-4.d0*x1/dens**3
    dgetaddn=2.d0/dens/dens*dx1ddn-4.d0*x1/dens**3
    dgetadgaa=2.d0/dens/dens*dx1dgaa
    dgetadgbb=2.d0/dens/dens*dx1dgbb
    dgetadgab=2.d0/dens/dens*dx1dgab
    detadup=dgetadup/2.d0/((3.d0*pi*pi*dens)**ex1)-eta/3.d0/dens
    detaddn=dgetaddn/2.d0/((3.d0*pi*pi*dens)**ex1)-eta/3.d0/dens
    detadgaa=dgetadgaa/2.d0/((3.d0*pi*pi*dens)**ex1)
    detadgbb=dgetadgbb/2.d0/((3.d0*pi*pi*dens)**ex1)
    detadgab=dgetadgab/2.d0/((3.d0*pi*pi*dens)**ex1)
 endif


! -- compute PBE-like part --

 call cGGA_u1(nfun,dens,zeta,grq,z,epscGGA,depscGGAdup,&
     & depscGGAddn,depscGGAdgaa,depscGGAdgbb,depscGGAdgab,&
     & depscGGAdz)

 ! If statement added to prevent empty-channel errors
 if (up.gt.threshold) then
   call cGGAspin_u1(nfun,up,gaa,z,epsGGA2,depsGGAdn2,depsGGAdgaa2,depsGGAdz2)
 else
   epsGGA2       = 0.d0
   depsGGAdn2    = 0.d0
   depsGGAdgaa2  = 0.d0
   depsGGAdz2    = 0.d0
 end if

 ! If statement added to prevent empty-channel errors
 if (dn.gt.threshold) then
   call cGGAspin_u1(nfun,dn,gbb,z,epsGGA3,depsGGAdn3,depsGGAdgaa3,depsGGAdz3)
 else
   epsGGA3       = 0.d0
   depsGGAdn3    = 0.d0
   depsGGAdgaa3  = 0.d0
   depsGGAdz3    = 0.d0
 end if

 if(epscGGA.ge.epsGGA2)then
    epscUPt=epscGGA
    depscUPtdup=depscGGAdup
    depscUPtddn=depscGGAddn
    depscUPtdgaa=depscGGAdgaa
    depscUPtdgbb=depscGGAdgbb
    depscUPtdgab=depscGGAdgab
    depscUPtdz=depscGGAdz
 else
    epscUPt=epsGGA2
    depscUPtdup=depsGGAdn2
    depscUPtddn=0.d0
    depscUPtdgaa=depsGGAdgaa2
    depscUPtdgbb=0.d0
    depscUPtdgab=0.d0
    depscUPtdz=depsGGAdz2
 endif

 if(epscGGA.ge.epsGGA3)then
    epscDNt=epscGGA
    depscDNtdup=depscGGAdup
    depscDNtddn=depscGGAddn
    depscDNtdgaa=depscGGAdgaa
    depscDNtdgbb=depscGGAdgbb
    depscDNtdgab=depscGGAdgab
    depscDNtdz=depscGGAdz
 else
    epscDNt=epsGGA3
    depscDNtdup=0.d0
    depscDNtddn=depsGGAdn3
    depscDNtdgaa=0.d0
    depscDNtdgbb=depsGGAdgaa3
    depscDNtdgab=0.d0
    depscDNtdz=depsGGAdz3
 endif


! -- compute the function C --

 if (nfun .eq. 1) then
    m=2
    n=3
    d=2.8d0
    c1=0.53d0
    c2=0.87d0
    c3=0.5d0
    c4=2.26d0
 else if (nfun .eq. 2) then
    m=2
    n=3
    d=2.8d0
    c1=0.59d0
    c2=0.9269d0
    c3=0.6225d0
    c4=2.1540d0
 else if (nfun .eq. 3) then
    m=2
    n=3
    d=4.5
    c1=0.35
    c2=0.87d0
    c3=0.5d0
    c4=2.26d0
 endif

 if(eta.le.1.d-28)then
    cz=c1+c2*zeta**2+c3*zeta**4+c4*zeta**6
    dczdzeta=2.d0*c2*zeta+4.d0*c3*zeta**3+6.d0*c4*zeta**5
    dczdeta=0.d0
    dczdz=0.d0
 else
    cz = dble(c1 + c2 * zeta ** 2 + c3 * zeta ** 4 + c4 * zeta ** 6) /&
  &      (0.1D1 + dble(eta ** 2 * ((1 + zeta) ** (-0.4D1 / 0.3D1) + &
  &      (1 - zeta) ** (-0.4D1 / 0.3D1))) / 0.2D1) ** 4

    cg = dble(2 * c2 * zeta + 4 * c3 * zeta ** 3 + 6 * c4 * zeta ** 5) /&
  &      (0.1D1 + dble(eta ** 2 * ((1 + zeta) ** (-0.4D1 / 0.3D1) + &
  &      (1 - zeta) ** (-0.4D1 / 0.3D1))) / 0.2D1) ** 4 - 0.2D1 * &
  &      dble(c1 + c2* zeta ** 2 + c3 * zeta ** 4 + c4 * zeta ** 6) /&
  &      (0.1D1 + dble(eta** 2 * ((1 + zeta) ** (-0.4D1 / 0.3D1) +&
  &      (1 - zeta) ** (-0.4D1 / 0.3D1))) / 0.2D1) ** 5 *&
  &      dble(eta ** 2) * (-0.4D1 / 0.3D1 * &
  &      dble((1 + zeta) ** (-0.7D1 / 0.3D1)) + 0.4D1 / 0.3D1 *&
  &      dble((1 - zeta) ** (-0.7D1 / 0.3D1)))
    dczdzeta=cg

    cg = -0.4D1 * dble(c1 + c2 * zeta ** 2 + c3 * zeta ** 4 + &
  &      c4 * zeta ** 6) / (0.1D1 + dble(eta ** 2 * &
  &      ((1 + zeta) ** (-0.4D1 / 0.3D1) + (1 - zeta) ** (-0.4D1 / 0.3D1)))/&
  &      0.2D1) ** 5 * dble(eta) * dble((1 + zeta) ** (-0.4D1 / 0.3D1) +&
  &      (1 - zeta) ** (-0.4D1 / 0.3D1))
    dczdeta=cg
    dczdz=0.d0
 endif

 dczdup=dczdzeta*dzetadup+dczdeta*detadup+dczdz*dzdup
 dczddn=dczdzeta*dzetaddn+dczdeta*detaddn+dczdz*dzddn
 dczdgaa=dczdeta*detadgaa+dczdz*dzdgaa
 dczdgbb=dczdeta*detadgbb+dczdz*dzdgbb
 dczdgab=dczdeta*detadgab+dczdz*dzdgab
 dczdta=dczdz*dzdta
 dczdtb=dczdz*dzdtb


! -- compute ecREV and its derivatives --
 y1=(1.d0+cz*z**m)
 dy1dup=dczdup*z**m+cz*m*z**(m-1)*dzdup
 dy1ddn=dczddn*z**m+cz*m*z**(m-1)*dzddn
 dy1dgaa=dczdgaa*z**m+cz*m*z**(m-1)*dzdgaa
 dy1dgbb=dczdgbb*z**m+cz*m*z**(m-1)*dzdgbb
 dy1dgab=dczdgab*z**m+cz*m*z**(m-1)*dzdgab
 dy1dta=dczdta*z**m + cz*m*z**(m-1)*dzdta
 dy1dtb=dczdtb*z**m + cz*m*z**(m-1)*dzdtb

 y2=1.d0+cz
 dy2dup=dczdup
 dy2ddn=dczddn
 dy2dgaa=dczdgaa
 dy2dgbb=dczdgbb
 dy2dgab=dczdgab
 dy2dta=dczdta
 dy2dtb=dczdtb

 y3=z**m
 dy3dup=m*z**(m-1)*dzdup
 dy3ddn=m*z**(m-1)*dzddn
 dy3dgaa=m*z**(m-1)*dzdgaa
 dy3dgbb=m*z**(m-1)*dzdgbb
 dy3dgab=m*z**(m-1)*dzdgab
 dy3dta=m*z**(m-1)*dzdta
 dy3dtb=m*z**(m-1)*dzdtb

 y4=up/dens*epscUPt+dn/dens*epscDNt
 dy4dup=1.d0/dens*epscUPt-up/dens/dens*epscUPt+&
  &     up/dens*depscUPtdup-dn/dens/dens*epscDNt+dn/dens*depscDNtdup
 dy4ddn=-up/dens/dens*epscUPt+up/dens*depscUPtddn+&
  &     1.d0/dens*epscDNt-dn/dens/dens*epscDNt+dn/dens*depscDNtddn
 dy4dgaa=up/dens*depscUPtdgaa+dn/dens*depscDNtdgaa
 dy4dgbb=up/dens*depscUPtdgbb+dn/dens*depscDNtdgbb
 dy4dgab=up/dens*depscUPtdgab+dn/dens*depscDNtdgab
 dy4dta=up/dens*depscUPtdz*dzdta+dn/dens*depscDNtdz*dzdta
 dy4dtb=up/dens*depscUPtdz*dzdtb+dn/dens*depscDNtdz*dzdtb

 epscREV=epscGGA*y1-y2*y3*y4
 depscREVdup = depscGGAdup*y1+epscGGA*dy1dup-&
  &            (dy2dup*y3*y4+y2*dy3dup*y4+y2*y3*dy4dup)
 depscREVddn = depscGGAddn*y1+epscGGA*dy1ddn-&
  &            (dy2ddn*y3*y4+y2*dy3ddn*y4+y2*y3*dy4ddn)
 depscREVdgaa = depscGGAdgaa*y1+epscGGA*dy1dgaa-&
  &             (dy2dgaa*y3*y4+y2*dy3dgaa*y4+y2*y3*dy4dgaa)
 depscREVdgbb = depscGGAdgbb*y1+epscGGA*dy1dgbb-&
  &             (dy2dgbb*y3*y4+y2*dy3dgbb*y4+y2*y3*dy4dgbb)
 depscREVdgab = depscGGAdgab*y1+epscGGA*dy1dgab-&
  &             (dy2dgab*y3*y4+y2*dy3dgab*y4+y2*y3*dy4dgab)
 depscREVdta = epscGGA*dy1dta+depscGGAdz*dzdta*y1-&
  &            dy2dta*y3*y4-y2*dy3dta*y4-y2*y3*dy4dta
 depscREVdtb = epscGGA*dy1dtb+depscGGAdz*dzdtb*y1-&
  &            dy2dtb*y3*y4-y2*dy3dtb*y4-y2*y3*dy4dtb


! -- compute ecMGGA and its derivatives
 y5=1.d0+d*epscREV*z**n
 dy5dup=d*depscREVdup*z**n+d*epscREV*n*z**(n-1)*dzdup
 dy5ddn=d*depscREVddn*z**n+d*epscREV*n*z**(n-1)*dzddn
 dy5dgaa=d*depscREVdgaa*z**n+d*epscREV*n*z**(n-1)*dzdgaa
 dy5dgbb=d*depscREVdgbb*z**n+d*epscREV*n*z**(n-1)*dzdgbb
 dy5dgab=d*depscREVdgab*z**n+d*epscREV*n*z**(n-1)*dzdgab
 dy5dta=d*depscREVdta*z**n+d*epscREV*n*z**(n-1)*dzdta
 dy5dtb=d*depscREVdtb*z**n+d*epscREV*n*z**(n-1)*dzdtb

 epscMGGA=epscREV*y5
 depscMGGAdup=depscREVdup*y5+epscREV*dy5dup
 depscMGGAddn=depscREVddn*y5+epscREV*dy5ddn
 depscMGGAdgaa=depscREVdgaa*y5+epscREV*dy5dgaa
 depscMGGAdgbb=depscREVdgbb*y5+epscREV*dy5dgbb
 depscMGGAdgab=depscREVdgab*y5+epscREV*dy5dgab
 depscMGGAdta=depscREVdta*y5+epscREV*dy5dta
 depscMGGAdtb=depscREVdtb*y5+epscREV*dy5dtb


! -- compute the final quantities --
 f=dens*epscMGGA
 fa=(epscMGGA+dens*depscMGGAdup)
 fb=(epscMGGA+dens*depscMGGAddn)
 fgaa=dens*depscMGGAdgaa
 fgbb=dens*depscMGGAdgbb
 fgab=dens*depscMGGAdgab
 fta=dens*depscMGGAdta
 ftb=dens*depscMGGAdtb

 f = scal*f
 fa = scal*fa
 fb = scal*fb
 fgaa = scal*fgaa
 fgab = scal*fgab
 fgbb = scal*fgbb
 fta = scal*fta
 ftb = scal*ftb

end subroutine tpss_cor

!=======================================================================!
!                                                                       !
!=======================================================================!

subroutine cGGA_u1(method,rho,zeta,grq,z,epsGGA,depsGGAdup,depsGGAddn,&
  &                depsGGAdgaa,depsGGAdgbb,depsGGAdgab,depsGGAdz)

!c if method=1, then PBE
!c if method=2, then GGA for revTPSS
!c if method=3, then PBEloc for TPSSloc

  implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: method
real(r8k), intent(in) :: rho
real(r8k), intent(in) :: zeta
real(r8k), intent(in) :: grq
real(r8k), intent(in) :: z

! output variables
real(r8k), intent(out) :: epsGGA
real(r8k), intent(out) :: depsGGAdup,depsGGAddn
real(r8k), intent(out) :: depsGGAdgaa,depsGGAdgbb,depsGGAdgab
real(r8k), intent(out) :: depsGGAdz

! local variables
real(r8k) :: up,dn,dens
real(r8k) :: gr,dgrdgaa,dgrdgab,dgrdgbb
real(r8k) :: dzetadup,dzetaddn
real(r8k) :: rs,drsdn,drsdup,drsddn
real(r8k) :: cg
real(r8k) :: phi,dphidzeta
real(r8k) :: akf,aks
real(r8k) :: t,dtdphi,dtdn,dtdup,dtddn
real(r8k) :: dtdgr,dtdgaa,dtdgab,dtdgbb
real(r8k) :: p,A,alpha1,beta1,beta2,beta3,beta4
real(r8k) :: G,dGdrs
real(r8k) :: epsc0,depsc0drs
real(r8k) :: epsc1,depsc1drs
real(r8k) :: alphac,dalphacdrs
real(r8k) :: ff,dffdzeta
real(r8k) :: epsc,depscdrs,depscdzeta
real(r8k) :: gamma1
real(r8k) :: beta,dbetadrs,dbetadz,dbetadt
real(r8k) :: bet1,c3,c4,bet2,dbet2drs
real(r8k) :: beta11,betaT,bx
real(r8k) :: ann
real(r8k) :: AA,dAAdphi,dAAdepsc,dAAdbeta
real(r8k) :: H,dHdAA,dHdt,dHdbeta
real(r8k) :: dHdrs,dHdzeta,dHdz

! constants
real(r8k), parameter :: pi = 3.1415926535897932d0
real(r8k), parameter :: ex1=0.333333333333333333d0

!-------------------------------------!

! compute relevant quantities

up=rho*(1.d0+zeta)/2.d0
dn=rho*(1.d0-zeta)/2.d0
dens=rho
gr=dsqrt(grq)
dgrdgaa=1.d0/2.d0/gr
dgrdgbb=1.d0/2.d0/gr
dgrdgab=1.d0/gr
dzetadup=0.1D1 / (up + dn) - (up - dn) / (up + dn) ** 2
dzetaddn=-0.1D1 / (up + dn) - (up - dn) / (up + dn) ** 2
rs=(3.d0/4.d0/pi/dens)**ex1
drsdn = -dble(3 ** (0.1D1 / 0.3D1)) * dble(2 ** (0.1D1 / 0.3D1)) * &
  &     0.3141592654D1 ** (-0.1D1 / 0.3D1) *&
  &     (0.1D1 / dens) ** (-0.2D1 /0.3D1) / dens ** 2 / 0.6D1
drsdup=drsdn
drsddn=drsdn
phi = dble((1 + zeta) ** (0.2D1 / 0.3D1)) / 0.2D1 + &
  &   dble((1 - zeta) ** (0.2D1 / 0.3D1)) / 0.2D1
cg = dble((1 + zeta) ** (-0.1D1 / 0.3D1)) / 0.3D1 - &
  &  dble((1 - zeta) ** (-0.1D1 / 0.3D1)) / 0.3D1
dphidzeta=cg

akf=(3.d0*pi*pi*dens)**(1.d0/3.d0)
aks=dsqrt(4.d0*akf/pi)
t=gr/2.d0/phi/aks/dens
dtdphi=-gr/2.d0/aks/dens/phi/phi
dtdn = -7.d0/6.d0*gr/2.d0/phi/dsqrt(4.d0/pi)/&
  &    ((3.d0*pi*pi)**(1.d0/6.d0))/(dens**(13.d0/6.d0))
dtdup=dtdn
dtddn=dtdn
dtdgr=1.d0/2.d0/phi/aks/dens
dtdgaa=dtdgr*dgrdgaa
dtdgbb=dtdgr*dgrdgbb
dtdgab=dtdgr*dgrdgab


! -- compute LDA correlation --
p=1.d0
A=0.031091d0
alpha1=0.21370d0
beta1=7.5957d0
beta2=3.5876d0
beta3=1.6382d0
beta4=0.49294d0
G = -0.2D1 * A * dble(1 + alpha1 * rs) * log(0.1D1 + 0.1D1 / A /&
  & (beta1 * sqrt(dble(rs)) + dble(beta2 * rs) + dble(beta3 * &
  & rs ** (0.3D1 / 0.2D1)) + dble(beta4 * rs ** (p + 1))) / 0.2D1)
dGdrs = -0.2D1 * A * alpha1 * log(0.1D1 + 0.1D1 / A / (beta1 * sqrt(rs) +&
  &      beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + beta4 * &
  &      rs **(p + 1)) / 0.2D1) + (0.1D1 + alpha1 * rs) / (beta1 * sqrt(rs) +&
  &      beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + &
  &      beta4 * rs ** (p + 1))** 2 * (beta1 * rs ** (-0.1D1 / 0.2D1) / &
  &      0.2D1 + beta2 + 0.3D1 /0.2D1 * beta3 * sqrt(rs) + &
  &      beta4 * rs ** (p + 1) * dble(p + 1) / rs) / (0.1D1 + 0.1D1 / &
  &      A / (beta1 * sqrt(rs) + beta2 * rs + beta3 *&
  &      rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (p + 1)) / 0.2D1)
epsc0=G
depsc0drs=dGdrs

p=1.d0
A=0.015545d0
alpha1=0.20548d0
beta1=14.1189d0
beta2=6.1977d0
beta3=3.3662d0
beta4=0.62517d0
G = -0.2D1 * A * dble(1 + alpha1 * rs) * log(0.1D1 + 0.1D1 / A / &
  & (beta1 * sqrt(dble(rs)) + dble(beta2 * rs) + &
  & dble(beta3 * rs ** (0.3D1 / 0.2D1)) + dble(beta4 * rs ** (p + 1))) / 0.2D1)
dGdrs = -0.2D1 * A * alpha1 * log(0.1D1 + 0.1D1 / A / (beta1 *&
  &     sqrt(rs) + beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) +&
  &     beta4 * rs **(p + 1)) / 0.2D1) + (0.1D1 + alpha1 * rs) /&
  &     (beta1 * sqrt(rs) + beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) +&
  &     beta4 * rs ** (p + 1)) ** 2 * (beta1 * rs ** (-0.1D1 / 0.2D1) /&
  &     0.2D1 + beta2 + 0.3D1 /0.2D1 * beta3 * sqrt(rs) + beta4 *&
  &     rs ** (p + 1) * dble(p + 1) / rs) / (0.1D1 + 0.1D1 / A / &
  &     (beta1 * sqrt(rs) + beta2 * rs + beta3 * &
  &     rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (p + 1)) / 0.2D1)
epsc1=G
depsc1drs=dGdrs

p=1.d0
A=0.016887d0
alpha1=0.11125d0
beta1=10.357d0
beta2=3.6231d0
beta3=0.88026d0
beta4=0.49671d0
G = -0.2D1 * A * dble(1 + alpha1 * rs) * log(0.1D1 + 0.1D1 / A / &
  & (beta1 * sqrt(dble(rs)) + dble(beta2 * rs) + dble(beta3 * &
  & rs ** (0.3D1 / 0.2D1)) + dble(beta4 * rs ** (p + 1))) / 0.2D1)
dGdrs = -0.2D1 * A * alpha1 * log(0.1D1 + 0.1D1 / A / (beta1 *&
  &    sqrt(rs) + beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + beta4 * &
  &    rs **(p + 1)) / 0.2D1) + (0.1D1 + alpha1 * rs) / (beta1 * sqrt(rs) +&
  &    beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + beta4 * &
  &    rs ** (p + 1))** 2 * (beta1 * rs ** (-0.1D1 / 0.2D1) / 0.2D1 + &
  &    beta2 + 0.3D1 /0.2D1 * beta3 * sqrt(rs) + &
  &    beta4 * rs ** (p + 1) * dble(p + 1) / rs) / (0.1D1 + 0.1D1 / A /&
  &    (beta1 * sqrt(rs) + beta2 * rs + beta3 *&
  &    rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (p + 1)) / 0.2D1)
alphac=-G
dalphacdrs=-dGdrs

ff = ((1 + zeta) ** (0.4D1 / 0.3D1) + (1 - zeta) ** (0.4D1 / 0.3D1) - 2) /&
  &  (2 * 2 ** (0.1D1 / 0.3D1) - 2)
cg = (0.4D1 / 0.3D1 * dble((1 + zeta) ** (0.1D1 / 0.3D1)) - &
  &  0.4D1/ 0.3D1 * dble((1 - zeta) ** (0.1D1 / 0.3D1))) / &
  &  dble(2 * 2 ** (0.1D1 / 0.3D1) - 2)
dffdzeta=cg

epsc = epsc0+alphac*ff*(1.d0-zeta**4)/1.709921d0+(epsc1-epsc0)*ff*zeta**4
depscdrs = depsc0drs+dalphacdrs*ff*(1.d0-zeta**4)/1.709921d0+&
  &       (depsc1drs-depsc0drs)*ff*zeta**4
depscdzeta = alphac/1.709921d0*(dffdzeta*(1.d0-zeta**4)- &
  &          ff*4.d0*zeta**3)+(epsc1-epsc0)*(dffdzeta*zeta**4+ff*4.d0*zeta**3)


! -- compute GGA part --

gamma1=0.031091d0
if (method .eq. 1) then
   beta=0.066725d0
   dbetadrs=0.d0
   dbetadz=0.d0
   dbetadt=0.d0
else if (method .eq. 2) then
   bet1=0.066725d0
   c3=0.1d0
   c4=0.1778d0
   bet2=(1.d0+c3*rs)/(1.d0+c4*rs)
   dbet2drs=c3/(1.d0+c4*rs)-(1.d0+c3*rs)/(1.d0+c4*rs)/(1.d0+c4*rs)*c4
   beta=bet1*bet2
   dbetadrs=bet1*dbet2drs
   dbetadz=0.d0
   dbetadt=0.d0
else if (method .eq. 3) then
   beta11=0.0375
   betaT=0.08
   bx=1.d0
   beta=beta11+betaT*t*t*(1.d0-dexp(-bx*rs**2))
   dbetadrs=betaT*t*t*2.d0*rs*dexp(-bx*rs**2)*bx
   dbetadt=betaT*2.d0*t*(1.d0-dexp(-bx*rs**2))
   dbetadz=0.d0
endif
ann=3.d0

AA = beta / gamma1 / (exp(-epsc / gamma1 / phi ** 3) - 0.1D1)
cg = -0.3D1 * beta / gamma1 ** 2 / (exp(-epsc/gamma1 / phi ** 3) -&
  &   0.1D1) ** 2 *epsc / phi ** 4 * exp(-epsc / gamma1 / phi ** 3)
dAAdphi=cg
cg = beta / gamma1 ** 2 / (exp(-epsc/ gamma1 / phi ** 3) - 0.1D1)** 2 /&
  &  phi ** 3 * exp(-epsc / gamma1 / phi ** 3)
dAAdepsc=cg
dAAdbeta=AA/beta

H = dble(gamma1) * phi ** ann * log(dble(1 + beta * t ** 2 /&
  & gamma1 * (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4)))
dHdAA = gamma1 * phi ** ann * (beta * t ** 4 / gamma1 / &
  &     (1 + AA *t ** 2 + AA ** 2 * t ** 4) - beta * t ** 2 / &
  &     gamma1 * (1 + AA * t** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4)**2 *&
  &     (t ** 2 + 2 *AA * t ** 4)) / (1 + beta * t ** 2 / gamma1 *&
  &     (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4))
dHdt = gamma1 * phi ** ann * (2 * beta * t / gamma1 * (1 + AA * t** 2) /&
  &    (1 + AA * t ** 2 + AA ** 2 * t ** 4) + &
  &    2 * beta * t ** 3 / gamma1 * AA/(1 + AA * t ** 2 + AA ** 2 * t ** 4) - &
  &    beta * t **2 / gamma1 * (1 + AA * t ** 2) /&
  &    (1 + AA * t ** 2 + AA ** 2 * t **4) ** 2 * (2 * AA * t + &
  &    4 * AA ** 2 * t ** 3)) / (1 + beta * t **2 / gamma1 *&
  &    (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4))
dHdbeta = phi**ann* t ** 2 * (1 + AA * t ** 2) / (1 + AA * t ** 2 + &
  &       AA ** 2 * t ** 4) / (1 + beta * t ** 2 / gamma1 * &
  &       (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4))
dHdt = dHdt+dHdbeta*dbetadt+dHdAA*dAAdbeta*dbetadt

dHdrs = dHdAA*dAAdepsc*depscdrs+dHdAA*dAAdbeta*dbetadrs+dHdbeta*dbetadrs
dHdzeta = ann*dphidzeta/phi*H+dHdAA*(dAAdphi*dphidzeta+ &
  &     dAAdepsc*depscdzeta)+dHdt*dtdphi*dphidzeta
dHdz = dHdAA*dAAdbeta*dbetadz+dHdbeta*dbetadz


! -- compute final values --

epsGGA = epsc+H
depsGGAdup = depscdrs*drsdup+depscdzeta*dzetadup+dHdrs*drsdup+&
  &          dHdzeta*dzetadup+dHdt*dtdup
depsGGAddn = depscdrs*drsddn+depscdzeta*dzetaddn+dHdrs*drsddn+&
  &          dHdzeta*dzetaddn+dHdt*dtddn
depsGGAdgaa = dHdt*dtdgaa
depsGGAdgbb = dHdt*dtdgbb
depsGGAdgab = dHdt*dtdgab
depsGGAdz = dHdz

end subroutine cGGA_u1


!=======================================================================!
!                                                                       !
!=======================================================================!

subroutine cGGAspin_u1(method,rho,gaa,z,epsGGA,depsGGAdn,depsGGAdgaa,depsGGAdz)

!c if method=1, then PBE
!c if method=2, then GGA for revTPSS
!c if method=3, then PBEloc

  implicit none

! select variables kinds
integer, parameter :: i4k = selected_int_kind(9) 
integer, parameter :: r8k = selected_real_kind(15,307)

! input variables
integer(i4k), intent(in) :: method
real(r8k), intent(in) :: rho
real(r8k), intent(in) :: gaa
real(r8k), intent(in) :: z

! output variables
real(r8k), intent(out) :: epsGGA
real(r8k), intent(out) :: depsGGAdn,depsGGAdgaa,depsGGAdz

! local variables
real(r8k) :: dens
real(r8k) :: gr,dgrdgaa
real(r8k) :: zeta
real(r8k) :: rs
real(r8k) :: drsdn
real(r8k) :: phi
real(r8k) :: akf,aks,t,dtdn,dtdgr,dtdgaa
real(r8k) :: p,A,alpha1,beta1,beta2,beta3,beta4
real(r8k) :: G,dGdrs
real(r8k) :: epsc1,depsc1drs
real(r8k) :: epsc,depscdrs
real(r8k) :: gamma1,beta,dbetadrs,dbetadz,dbetadt
real(r8k) :: bet1,c3,c4,bet2,dbet2drs
real(r8k) :: beta11,betaT,bx
real(r8k) :: AA,cg,dAAdepsc,dAAdbeta
real(r8k) :: H,dHdAA,dHdt,dHdbeta,dHdrs,dHdz
real(r8k) :: ann

! constants
real(r8k), parameter :: pi = 3.1415926535897932d0
real(r8k), parameter :: ex1=0.333333333333333333d0

!-------------------------------------!

! compute relevant quantities

dens=rho
gr=dsqrt(gaa)
dgrdgaa=1.d0/2.d0/gr
zeta=1.d0
rs=(3.d0/4.d0/pi/dens)**ex1
drsdn = -dble(3 ** (0.1D1 / 0.3D1)) * dble(2 ** (0.1D1 / 0.3D1)) * &
  &     0.3141592654D1 ** (-0.1D1 / 0.3D1) * &
  &     (0.1D1 / dens) ** (-0.2D1 / 0.3D1) / dens ** 2 / 0.6D1
phi = dble((1 + zeta) ** (0.2D1 / 0.3D1)) / 0.2D1 + &
  &   dble((1 - zeta) ** (0.2D1 / 0.3D1)) / 0.2D1

akf = (3.d0*pi*pi*dens)**(1.d0/3.d0)
aks = dsqrt(4.d0*akf/pi)
t = gr/2.d0/phi/aks/dens
dtdn = -7.d0/6.d0*gr/2.d0/phi/dsqrt(4.d0/pi)/&
  &   ((3.d0*pi*pi)**(1.d0/6.d0))/(dens**(13.d0/6.d0))
dtdgr=1.d0/2.d0/phi/aks/dens
dtdgaa=dtdgr*dgrdgaa


! -- compute LDA correlation --

p=1.d0
A=0.015545d0
alpha1=0.20548d0
beta1=14.1189d0
beta2=6.1977d0
beta3=3.3662d0
beta4=0.62517d0
G = -0.2D1 * A * dble(1 + alpha1 * rs) * log(0.1D1 + 0.1D1 / A /&
  & (beta1 * sqrt(dble(rs)) + dble(beta2 * rs) + dble(beta3 * &
  & rs ** (0.3D1 / 0.2D1)) + dble(beta4 * rs ** (p + 1))) / 0.2D1)
dGdrs = -0.2D1 * A * alpha1 * log(0.1D1 + 0.1D1 / A / (beta1 *& 
  &     sqrt(rs) + beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + &
  &     beta4 * rs **(p + 1)) / 0.2D1) + (0.1D1 + alpha1 * rs) /&
  &     (beta1 * sqrt(rs) + beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) +&
  &     beta4 * rs ** (p + 1))** 2 * (beta1 * rs ** (-0.1D1 / 0.2D1) /&
  &     0.2D1 + beta2 + 0.3D1 /0.2D1 * beta3 * sqrt(rs) + &
  &     beta4 * rs ** (p + 1) * dble(p + 1) / rs) / &
  &     (0.1D1 + 0.1D1 / A / (beta1 * sqrt(rs) + beta2 * rs + beta3 *&
  &     rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (p + 1)) / 0.2D1)
epsc1=G
depsc1drs=dGdrs

epsc = epsc1
depscdrs = depsc1drs


! -- compute the GGA part --

gamma1=0.031091d0
if (method .eq. 1) then
   beta=0.066725d0
   dbetadrs=0.d0
   dbetadz=0.d0
   dbetadt=0.d0
else if (method .eq. 2) then
   bet1=0.066725d0
   c3=0.1d0
   c4=0.1778d0
   bet2=(1.d0+c3*rs)/(1.d0+c4*rs)
   dbet2drs=c3/(1.d0+c4*rs)-(1.d0+c3*rs)/(1.d0+c4*rs)/(1.d0+c4*rs)*c4
   beta=bet1*bet2
   dbetadrs=bet1*dbet2drs
   dbetadz=0.d0
   dbetadt=0.d0
else if (method .eq. 3) then
   beta11=0.0375
   betaT=0.08
   bx=1.d0
   beta=beta11+betaT*t*t*(1.d0-dexp(-bx*rs**2))
   dbetadrs=betaT*t*t*2.d0*rs*dexp(-bx*rs**2)*bx
   dbetadt=betaT*2.d0*t*(1.d0-dexp(-bx*rs**2))
   dbetadz=0.d0
endif
ann=3.d0

AA = beta / gamma1 / (exp(-epsc / gamma1 / phi ** 3) - 0.1D1)
cg = beta / gamma1 ** 2 / (exp(-epsc/ gamma1 / phi ** 3) - 0.1D1)** 2 /&
  &  phi ** 3 * exp(-epsc / gamma1 / phi ** 3)
dAAdepsc=cg
dAAdbeta=AA/beta

H = dble(gamma1) * phi ** ann * log(dble(1 + beta * t ** 2 /&
  & gamma1 * (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4)))
dHdAA = gamma1 * phi ** ann * (beta * t ** 4 / gamma1 / &
  &    (1 + AA *t ** 2 + AA ** 2 * t ** 4) - beta * t ** 2 / gamma1 * &
  &    (1 + AA * t** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4) ** 2 * &
  &    (t ** 2 + 2 *AA * t ** 4)) / (1 + beta * t ** 2 / gamma1 * &
  &    (1 + AA * t ** 2) /(1 + AA * t ** 2 + AA ** 2 * t ** 4))
dHdt = gamma1 * phi ** ann * (2 * beta * t / gamma1 * (1 + AA * t** 2) /&
  &    (1 + AA * t ** 2 + AA ** 2 * t ** 4) + 2 * beta * t ** 3 /&
  &    gamma1 * AA / (1 + AA * t ** 2 + AA ** 2 * t ** 4) - beta * t **2 /&
  &    gamma1 * (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t **4)**2 *&
  &    (2 * AA * t + 4 * AA ** 2 * t ** 3)) / (1 + beta * t **2 / gamma1 *&
  &    (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4))
dHdbeta = phi**ann* t ** 2 * (1 + AA * t ** 2) / (1 + AA * t ** 2 + &
  &       AA ** 2 * t ** 4) / (1 + beta * t ** 2 / gamma1 * &
  &       (1 + AA * t ** 2) / (1 + AA * t ** 2 + AA ** 2 * t ** 4))
dHdt=dHdt+dHdbeta*dbetadt+dHdAA*dAAdbeta*dbetadt
dHdrs=dHdAA*dAAdepsc*depscdrs+dHdAA*dAAdbeta*dbetadrs+dHdbeta*dbetadrs
dHdz=dHdAA*dAAdbeta*dbetadz+dHdbeta*dbetadz


! -- compute final values --

epsGGA = epsc+H
depsGGAdn = depscdrs*drsdn+dHdrs*drsdn+dHdt*dtdn
depsGGAdgaa = dHdt*dtdgaa
depsGGAdz = dHdz

end subroutine cGGAspin_u1
