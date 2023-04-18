!****s*  FHI-aims/tf_ini
!  NAME
!    tf_ini
!  SYNOPSIS

subroutine tf_ini_trans(ntau,nomega,taumax,omegamax,tau,omega,wtau,womega,output)
!  PURPOSE
!
! initialises the frequency axis with modified Gauss-Legendre grids in the range
! (0.,omegamax) and (0.,taumax), respectively. Here omegamax for screened Coulomb
! potential and taumax for the self-energy are chosen to be the same, in contrast 
! to the "tf_ini" case.  
! Uses transformation of weights w = w_0 (1+w)/(1-w) to accelerate quadrature convergence.
!
! The abscissa and weights on the axis are chosen such that
! sum_i^N P_N(x_i) w_i is equal to int_x1^xN P_N(x) dx for a Legendre polynomial
! of order N (done in the routine GAULEG).
!

!  USES
   use mpi_tasks
   use localorb_io,only:use_unit
   implicit none
!  INPUTS
!   o ntau -- of points on the self-energy frequency axis
!   o nomega --  of points on the ccreened coulomb and poalrizability frequency axis.
!   o taumax --  cutoff on self-energy frequency axis
!   o omegamax -- cutoff on frequency axis
!   o output --  flag for control output
!
!  OUTPUT
!  o tau -- self-energy frequency grid
!  o omega -- omega grid
!  o wtau -- weights for self-energy frequency grid
!  o womega -- weights for omega grid
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!  ARGUMENTS
   integer,          intent(in)  :: ntau, nomega
   double precision, intent(in)  :: taumax, omegamax
   double precision, intent(out) :: tau(ntau), omega(nomega), wtau(ntau), womega(nomega)
   logical,          intent(in)  :: output
   
   double precision :: omega0
   integer          :: i


  if(myid.eq.0) then
  if(output) write(use_unit,'(/2x,2a)') 'Initialising transformed Gauss-Legendre time and ',&
                              'frequency grids'
  endif

! do some checks
  if(ntau<0)     call aims_stop ('ntau < 0','tf_ini')
  if(nomega<0)   call aims_stop ('nomega < 0','tf_ini')
  if(taumax<0)   call aims_stop ('taumax < 0','tf_ini')
  if(omegamax<0) call aims_stop ('omegamax < 0','tf_ini')

  omega0 = 0.5d0

! calculate Gauss-Legendre grids and weights for omega and tau
!  call gauleg(0,taumax,tau(1:ntau),wtau(1:ntau),ntau)
  call gauleg(-1.d0,1.d0, tau(1:ntau),wtau(1:ntau),ntau)
!  call gauleg(0d0,taumax,tau(1:ntau),wtau(1:ntau),ntau)
  call gauleg(-1.d0,1.d0, omega(1:nomega),womega(1:ntau),nomega)
  do i = 1, nomega
   womega(i) = womega(i) * 2.0d0*omega0/( (1.0d0-omega(i))**2 )
   omega(i) = omega0*(1.0d0+omega(i))/(1.0d0-omega(i))
  enddo
  do i = 1, ntau
   wtau(i) = wtau(i) * 2.0d0*omega0/( (1.0d0-tau(i))**2 )
   tau(i) = omega0*(1.0d0+tau(i))/(1.0d0-tau(i))
  enddo

!  call gauleg1(omegamax, omega(ntau:nomega),womega(ntau:nomega),nomega-ntau+1)
!  tau(ntau)        = 1.2*taumax  ! tail fit points
!  wtau(ntau)       = 0d0
!  omega(nomega)    = 1.2*omegamax
!  womega(nomega)   = 0d0

! control output
  if(output) then
   if(myid.eq.0) then
!    write(use_unit,'(2X,a,f7.3)') 'Frequency range for self energy...............:',taumax
!    write(use_unit,'(2X,a,f7.3)') 'Frequency range...............................:',omegamax
    write(use_unit,'(2X,a,i3)')   'Number of frequency points for self energy.....:',ntau
    write(use_unit,'(2X,a,i3)')   'Number of frequency points....................:',nomega
   endif
  endif

end subroutine
!******


subroutine tf_ini(ntau,nomega,taumax,omegamax,tau,omega,wtau,womega,output)

!  PURPOSE
! 
! initialises the time and frequency axis with Gauss-Legendre grids in the range
! (0..omegamax) and (0..taumax), respectively.
!
! The abscissa and weights on the axis are chosen such that
! sum_i^N P_N(x_i) w_i is equal to int_x1^xN P_N(x) dx for a Legendre polynomial
! of order N (done in the routine GAULEG).
!
! In order to utilise the analytical enhancements the last point on the axis
! is chosen at 1.2 * xmax (where xmax= taumax or ommegamax). This point is
! referred to as the 'tail fitting point'.
! Note that the tail fitting points don't contribute any wait to the integration
!
! Further note that the intervals are open, i.e. (0..taumax)

!  USES

   use mpi_tasks
   use localorb_io,only:use_unit
   implicit none

!  ARGUMENTS
   
   integer,          intent(in)  :: ntau, nomega  
   double precision, intent(in)  :: taumax, omegamax
   double precision, intent(out) :: tau(ntau), omega(nomega), wtau(ntau), womega(nomega) 
   logical,          intent(in)  :: output
   
!  INPUTS
!   o ntau -- of points on time axis
!   o nomega --  of points on frequency axis.
!   o taumax --  cutoff on time axis
!   o omegamax -- cutoff on frequency axis
!   o output --  flag for control output
!
!  OUTPUT
!  o tau -- time grid
!  o omega -- omega grid
!  o wtau -- weights for time grid
!  o womega -- weights for omega grid
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE






  if(myid.eq.0) then
  if(output) write(use_unit,'(/2x,2a)') 'Initialising Gauss-Legendre time and ',&
                              'frequency grids'
  endif 

! do some checks
  if(ntau<0)     call aims_stop ('ntau < 0','tf_ini')
  if(nomega<0)   call aims_stop ('nomega < 0','tf_ini')
  if(taumax<0)   call aims_stop ('taumax < 0','tf_ini')
  if(omegamax<0) call aims_stop ('omegamax < 0','tf_ini')

! calculate Gauss-Legendre grids and weights for omega and tau 
  call gauleg(0d0,taumax,tau(1:ntau),wtau(1:ntau),ntau)
  call gauleg(0d0,omegamax, omega(1:ntau),womega(1:ntau),ntau)
  call gauleg1(omegamax, omega(ntau:nomega),womega(ntau:nomega),nomega-ntau+1)
!  tau(ntau)        = 1.2*taumax  ! tail fit points
!  wtau(ntau)       = 0d0
!  omega(nomega)    = 1.2*omegamax
!  womega(nomega)   = 0d0

! control output
  if(output) then
   if(myid.eq.0) then
    write(use_unit,'(2X,a,f7.3)') 'Frequency range for self energy...............: ',taumax
    write(use_unit,'(2X,a,f7.3)') 'Frequency range...............................: ',omegamax
    write(use_unit,'(2X,a,i3)')   'Number of freuency points for self energy.....: ',ntau
    write(use_unit,'(2X,a,i3)')   'Number of frequency points....................: ',nomega
   endif
  endif  

end subroutine
!******

!-------------------------------------------------------------------------
! tf_ini
! 
! Gauss Legendre abscissa and weights (P.R),
! Given the lower and upper limits of integration x1 and x2, and given n,
! the routine returns arrays x(1:n) and w(1:n) of length n, containing the
! abcsissas and weights of the Gauss-Legendre n-point quadrature formula.
!-------------------------------------------------------------------------
subroutine gauleg(x1,x2,x,w,n)
use constants
implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(n),w(n)   ! abcsissa and weigths

! internal vars
double precision :: eps
Parameter (eps=3D-14) ! eps is the relative precision
integer :: i,j,m
double precision :: p1,p2,p3,pp,xl,xm,z,z1

! The roots are symmetric in the interval, so we only have to find half
! of them
  m=(n+1)/2
  xm=0.5D0*(x2+x1)
  xl=0.5D0*(x2-x1)
! loop over the desired roots
  do i=1,m
    z=cos(pi*(i-.25D0)/(n+.5D0))
!   starting with the above approximation to the ith root, we enter the 
!   main loop of refinement by Newton's method.
  1       continue
      p1=1D0
      p2=0D0
!     Loop up the recurrence relation to get the Legendre polynomial
!     evaluated at z.
      do j=1,n
	p3=p2
	p2=p1
	p1=((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
      enddo
!     p1 is now the desired Legendre polynomial. We next compute pp, its
!     derivative, by a standard rlation involving also p2, the polynomial
!     of one lower order
      pp=n*(z*p1-p2)/(z*z-1D0)
      z1=z
!     Newton's method
      z=z1-p1/pp	    
    if(abs(z-z1).gt.eps) goto 1
!   Scale the root to the desired interval, and put in its symmetric counter
!   part.
    x(i)=xm-xl*z
    x(n+1-i)=xm+xl*z
!   compute the weight and its symmetric counterpart
    w(i)=2D0*xl/((1D0-z*z)*pp*pp)
    w(n+1-i)=w(i)
  enddo                	 
  return
end

! this subroutine is slightly changed from the original one
! for the gaussian-legendre grid in the range (w0, infty).
subroutine gauleg1(x1,x,w,n)
use constants
implicit none

integer,          intent(in)  :: n           ! number of points
! x1 is the lower bound, the upper bound is the positive infinity by default 
double precision, intent(in)  :: x1          ! ranges
double precision, intent(out) :: x(n),w(n)   ! abcsissa and weigths

! internal vars
double precision :: eps
Parameter (eps=3D-14) ! eps is the relative precision
integer :: i,j,m
double precision :: p1,p2,p3,pp,xl,xm,z,z1

! The roots are symmetric in the interval, so we only have to find half
! of them
  m=(n+1)/2
!  xm=0.5D0*(x2+x1)
!  xl=0.5D0*(x2-x1)
! loop over the desired roots
  do i=1,m
    z=cos(pi*(i-.25D0)/(n+.5D0))
!   starting with the above approximation to the ith root, we enter the 
!   main loop of refinement by Newton's method.
  1       continue
      p1=1D0
      p2=0D0
!     Loop up the recurrence relation to get the Legendre polynomial
!     evaluated at z.
      do j=1,n
	p3=p2
	p2=p1
	p1=((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
      enddo
!     p1 is now the desired Legendre polynomial. We next compute pp, its
!     derivative, by a standard rlation involving also p2, the polynomial
!     of one lower order
      pp=n*(z*p1-p2)/(z*z-1D0)
      z1=z
!     Newton's method
      z=z1-p1/pp	    
    if(abs(z-z1).gt.eps) goto 1
!   Scale the root to the desired interval, and put in its symmetric counter
!   part.
    x(i)=2*x1/(z+1.D0)
    x(n+1-i)=2*x1/(1.D0-z)
!   compute the weight and its symmetric counterpart
    w(i)=4D0*x1/((1D0-z*z)*pp*pp)/(1.d0+z)/(1.d0+z)
    w(n+1-i)= 4D0*x1/((1D0-z*z)*pp*pp)/(1.d0-z)/(1.d0-z)
  enddo                	 
  return
end

