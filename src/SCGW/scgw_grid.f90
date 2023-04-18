      module scgw_grid 

      ! containes subroutines for the definition of the frequency/time grid
      ! the one currently used in scGW is the logarithmic ones, that works very for
      ! fourier transforms
      ! This module at a certain point should be merged with tf_ini.f90

      use dimensions 
      use gw_para
      implicit none
      integer ::  nomega
      integer ::  ntau
      real*8  ::  freqmax1
      real*8  ::  taumax1
      real*8, allocatable ::   omega(:)
      real*8, allocatable ::   womega1(:)
!      real*8, allocatable ::   womega(:)
      real*8, allocatable ::   tau(:)
      real*8, allocatable ::   wtau(:)
      integer , dimension(:,:) , allocatable :: map_index
      integer , dimension(:,:) , allocatable :: map_index_t
      integer l_freq
      real*8 t_0, w_0
!      integer :: n_loc_grid 
!      integer :: n_loc_grid_t 

      contains 
!----------------------------------------------------------------------      
       subroutine init_grid ()

          use gw_para
          use mpi_tasks

          implicit none 

          integer n_remain
                     

          nomega = n_full_freq
          ntau = n_full_time
          freqmax1 = omegamax
          taumax1 = taumax

          allocate (tau (-ntau : ntau))
          allocate (wtau (-ntau : ntau))
          allocate (omega (nomega))
          allocate (womega1 (nomega))

          call tf_ini_log (ntau,nomega, taumax1,freqmax1, &
            tau, omega, &
            wtau,womega1,.true.)

        if(.false. .and. myid == 0 )then
          open (44,file = 'prova_grid')
          do l_freq = 1, nomega, 1
            write(44,*) l_freq, omega(l_freq)
          enddo
          close(44)
        endif

        n_remain = MOD(nomega, n_tasks)
        if (n_remain.eq.0) then
          n_loc_grid = nomega / n_tasks
        else
          n_loc_grid = nomega / n_tasks + 1
        endif
       allocate(map_index(n_tasks, n_loc_grid))



       call distribute_grid (nomega, map_index,n_loc_grid)

!stop
       end subroutine init_grid 

       subroutine reset_scgw_grid ()
         nomega = 0
         ntau = 0
         freqmax1 = 0
         taumax1 = 0
         l_freq = 0
         t_0 = 0
         w_0 = 0

         if (allocated(omega)) deallocate(omega)
         if (allocated(womega1)) deallocate(womega1)
         if (allocated(tau)) deallocate(tau)
         if (allocated(wtau)) deallocate(wtau)
         if (allocated(map_index)) deallocate(map_index)
         if (allocated(map_index_t)) deallocate(map_index_t)
       end subroutine reset_scgw_grid


subroutine tf_ini_log(ntau,nomega,taumax,omegamax,tau,omega,wtau,womega,output)

      use localorB_io, only: use_unit
      use mpi_tasks
      use synchronize_mpi

   implicit none

!  ARGUMENTS
   
   integer,          intent(in)  :: ntau, nomega  
   double precision, intent(in)  :: taumax, omegamax
   double precision, intent(out) :: tau(-ntau:ntau), omega(nomega), wtau(-ntau:ntau), womega(nomega) 
   logical,          intent(in)  :: output
   integer i
   real*8 omega0 
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

  if(myid.eq.0)then
    if(output) write(use_unit,'(/2x,2a)') 'Initialising logarithmic time and ',&
                              'frequency grids'
  endif

! do some checks
  if(ntau<0)     stop 'ntau < 0'
  if(nomega<0)   stop 'nomega < 0'
  if(taumax<0)   stop 'taumax < 0'
  if(omegamax<0) stop 'omegamax < 0'

! calculate Gauss-Legendre grids and weights for omega and tau 

  call loggrid_time(0d0,taumax,tau(-ntau:ntau),wtau(-ntau:ntau),ntau)
  call loggrid_freq(0d0,omegamax, omega(1:nomega),womega(1:nomega),nomega)
!if(use_dmft_gw.or.use_dmft_pbe0)&
!omega(1)=omega(2)/5
!  omega0 = 0.5d0
!  call gauleg(-1.d0,1.d0, omega(1:nomega),womega(1:ntau),nomega)
!  do i = 1, nomega
!   womega(i) = womega(i) * 2.0d0*omega0/( (1.0d0-omega(i))**2 )
!   omega(i) = omega0*(1.0d0+omega(i))/(1.0d0-omega(i))
!   omega(i) = omega(i)-omega(1) !first point placed at zero
!  enddo
  
 
  if(myid.eq.0)then
    do i = 1, nomega
      write(use_unit,*) ' i_freq ',i,omega(i)
    enddo
  endif
!stop

!  call tf_ini_trans(n_freq,n_full_freq, omegamax,omegamax, &
!     omega_grid,omega_full_grid, &
!     womega,womega_full,.true.)


!  call gauleg(0d0,omegamax, omega(1:nomega/2),womega(1:nomega/2),nomega/2)
!  call gauleg1(omegamax, omega(nomega/2+1:nomega),womega(nomega/2+1:nomega),nomega/2+1)

!  call gauleg(-1.d0,1.d0, omega(1:nomega),womega(1:ntau),nomega)
!  do i = 1, nomega
!   womega(i) = womega(i) * 2.0d0*0.5d0/( (1.0d0-omega(i))**2 )
!   omega(i) = 0.5d0*(1.0d0+omega(i))/(1.0d0-omega(i))
!  enddo

! control output
  if(output) then
   if(myid.eq.0)then
    write(use_unit,'(2X,a,f7.1)') 'Time range....................................: ',taumax
    write(use_unit,'(2X,a,f7.1)') 'Frequency range...............................: ',omegamax
    write(use_unit,'(2X,a,i4)')   'Number of time points ........................: ',ntau
    write(use_unit,'(2X,a,i4)')   'Number of frequency points....................: ',nomega
    write(use_unit,*) " "
   endif  
  endif  

end subroutine

!-------------------------------------------------------------------------
! tf_ini
! 
! Gauss Legendre abscissa and weights (P.R),
! Given the lower and upper limits of integration x1 and x2, and given n,
! the routine returns arrays x(1:n) and w(1:n) of length n, containing the
! abcsissas and weights of the Gauss-Legendre n-point quadrature formula.
!-------------------------------------------------------------------------
subroutine loggrid_time(x1,x2,x,w,n)

implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(-n:n),w(-n:n)   ! abcsissa and weigths

! internal vars
real*8 h, a
!real*8 t_0
integer :: i,j,m

  t_0 = 0.01
  h = 1./real(n)*log ((x2-x1)/t_0)
  a = 1.05

  do i=1, n, 1
    x(i) = t_0 * (exp(i*h)- 1.)
    x(-i) = -x(i)
    w(i) = h * t_0 * exp(i*h) 
    w(-i) = w(i)
  enddo

  x(0)=0
  w(0)=t_0 * h 

  return
end subroutine loggrid_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine loggrid_freq(x1,x2,x,w,n)

implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(1:n),w(1:n)   ! abcsissa and weigths

! internal vars
real*8 h
!real*8 w_0 
integer :: i,j,m

  w_0 = 0.01!5
  h = 1./real(n)*log ((x2-x1)/w_0)
  
  !if(.not. use_scgw  .and. .not. use_scgw0)then
  if(.not. use_scgw  .and. .not. use_scgw0 .and. .not. use_dmft_gw &
     .and. .not. use_dmft_pbe0)then

    w_0 = 0.001
    h = 1./real(n)*log ((x2-x1)/w_0)
  endif


  do i=1, n, 1
    x(i) = w_0 * (exp((i-1)*h)-1.)
    w(i)  = h * w_0 * exp((i-1)*h)
    !if(.not. use_scgw  .and. .not. use_scgw0 )then
    if(.not. use_scgw  .and. .not. use_scgw0 .and. .not. use_dmft_gw &
     .and. .not. use_dmft_pbe0)then

       x(i) = w_0 * (exp((i-1)*h)-1.)
       w(i)  = h * w_0 * exp((i-1)*h)
    endif

  enddo

  return
end subroutine loggrid_freq

!******

      subroutine distribute_grid (n_points, map_index, n_loc_grid1)
      use localorb_io, only: use_unit
      use mpi_tasks

      implicit none

      integer n_points
      integer n_loc_grid1
      integer n_max_loc
      integer n_remain
      integer map_index ( n_tasks, n_loc_grid1 ) 

!counters
      integer i_index
      integer i_task
      integer i,n 

       n_remain = MOD(n_points, n_tasks)

       map_index = 0
       i_index = 0
       do i_task = 1, n_tasks

         if((n_remain.eq.0) .or. &
            i_task.le.n_remain ) then

            do i = 1, n_loc_grid1
             i_index = i_index + 1
             map_index( i_task, i) = i_index
            enddo

         else

            do i = 1, n_loc_grid1 -1
             i_index = i_index + 1
             map_index(i_task, i) = i_index
            enddo

         endif
       enddo

       if (i_index .ne. n_points) then
         write(use_unit,*) " * Error in the grid distribution. Stop! "
         stop
       endif

      end subroutine distribute_grid
 
!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------
!     OLD VERSIONS
!

subroutine homolog_freq(x1,x2,x,w,n)

use mpi_tasks

implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(1:n),w(1:n)   ! abcsissa and weigths

! internal vars
real*8 h, t_0, a, sum_v,dx
 
integer :: i,j,m, i_index, n_points

  t_0 = 0.01
  h = 1./real(n/5)*log ((x2-x1)/t_0/5.d0/8)
  i_index = 0
  sum_v = x1

  n_points = 5 
 
  do i=1, n/n_points, 1
    dx = t_0 * (exp((i)*h)-1.)
    !if(myid.eq.0) print *,dx
   do m = 1,n_points,1
    i_index = i_index + 1
    x(i_index)  = sum_v 
    w(i_index)  = 0.d0
    sum_v = sum_v + dx
   enddo
  enddo

  return
end subroutine homolog_freq


subroutine quadratic_grid (x1,x2,x,w,n)

use mpi_tasks

implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(1:n),w(1:n)   ! abcsissa and weigths

! internal vars
real*8 h, t_0, a, sum_v,dx

integer :: i,j,m, i_index

  a = x2 / n**3
!  t_0 = 2.d0
!  h = 1./real(n/5)*log ((x2-x1)/t_0/5.d0/8)
!  i_index = 0
!  sum_v = 0.d0
  do i=1, n, 1
    x(i)  = a* (i-1)**3 
    w(i)  = 3*a*(i-1)**2
   enddo

  return
end subroutine quadratic_grid
!-----------------------
!subroutine tf_ini_homo2(ntau,nomega,taumax,omegamax,tau,omega,wtau,womega,output)
subroutine tf_ini_homo2(n_freq,n_full_freq, omegamax, &
     omega_grid,omega_full_grid, &
     womega,womega_full,output)

      use localorb_io, only: use_unit
      use mpi_tasks
      use synchronize_mpi

   implicit none

!  ARGUMENTS
   integer,          intent(in)  :: n_freq, n_full_freq
   double precision, intent(in)  :: omegamax
   double precision, intent(out) :: omega_grid      (n_freq)
   double precision, intent(out) :: womega          (n_freq) 
   double precision, intent(out) :: omega_full_grid (n_full_freq) 
   double precision, intent(out) :: womega_full     (n_full_freq)
   logical,          intent(in)  :: output
   integer i_index,i
   integer n_full_freq1

  if(myid.eq.0)then
    if(output) write(use_unit,'(/2x,2a)') 'Initialising the homogeneous ',&
                              'frequency grid'
  endif

! do some checks
  if(n_freq<0)     stop 'n_freq < 0'
  if(n_full_freq<0)   stop 'n_full_freq < 0'
  if(omegamax<0) stop 'omegamax < 0'

  if(myid.eq.0)write(use_unit,*)n_freq, n_full_freq 
  call homo_grid_freq(0d0,omegamax, omega_grid     , womega,      n_freq)
  call homo_grid_freq(0d0,omegamax, omega_full_grid, womega_full, n_full_freq)

! control output
  if(output) then
   if(myid.eq.0)then
    write(use_unit,'(2X,a,f7.3)') 'Frequency range...............................: ',omegamax
    write(use_unit,'(2X,a,i3)')   'Number of freuency points for self energy.....: ',n_freq
    write(use_unit,'(2X,a,i3)')   'Number of frequency points....................: ',n_full_freq
    write(use_unit,*) " "
   endif
  endif

end subroutine



!****s*  FHI-aims/tf_ini
!  NAME
!    tf_ini
!  SYNOPSIS
subroutine tf_ini_homo(ntau,nomega,taumax,omegamax,tau,omega,wtau,womega,output)

      use localorb_io, only: use_unit
      use mpi_tasks
      use synchronize_mpi

   implicit none

!  ARGUMENTS
   integer,          intent(in)  :: ntau, nomega
   double precision, intent(in)  :: taumax, omegamax
   double precision, intent(out) :: tau(-ntau:ntau), omega(nomega), wtau(-ntau:ntau), womega(nomega)
   logical,          intent(in)  :: output
   integer i_index,i
   integer nomega1

  if(myid.eq.0)then
    if(output) write(use_unit,'(/2x,2a)') 'Initialising logarithmic time and ',&
                              'frequency grids'
  endif

! do some checks
  if(ntau<0)     stop 'ntau < 0'
  if(nomega<0)   stop 'nomega < 0'
  if(taumax<0)   stop 'taumax < 0'
  if(omegamax<0) stop 'omegamax < 0'

! calculate Gauss-Legendre grids and weights for omega and tau 

!  call loggrid_time(0d0,taumax,tau(-(ntau-1):ntau-1),wtau(-(ntau-1):ntau-1),ntau-1)
  call homo_grid_time(0d0,taumax,tau(-(ntau):(ntau)),wtau(-(ntau):(ntau)),ntau)
!  call loggrid_freq(0d0,omegamax, omega(1:nomega-3),womega(1:nomega-3),nomega-3)
  call homo_grid_freq(0d0,omegamax, omega(1:nomega),womega(1:nomega),nomega)
  nomega1 = 200
!   call homolog_freq(0.d0,omegamax, omega(1:nomega),womega(1:nomega),nomega)
!   call quadratic_grid(0.d0,omegamax, omega(1:nomega),womega(1:nomega),nomega)
!   call homo_grid_freq (100.d0,omegamax, omega(nomega/2+1:nomega),womega(nomega/2+1:nomega),nomega/2)
!  call homolog_freq(0d0,5.d0, omega(1:nomega1),womega(1:nomega1),nomega1)
!  call loggrid_freq (0d0,5.d0, omega(1:nomega1),womega(1:nomega1),nomega1)
!  call homo_grid_freq(omega(nomega1),omegamax, omega(nomega1+1:nomega),womega(nomega1+1:nomega),nomega-nomega1)
!  call homo_grid_freq(omega(nomega1),omegamax, omega(nomega1+1:nomega),womega(nomega1+1:nomega),nomega-nomega1)

!  call gauleg(0d0, omegamax, omega(1:nomega-1),womega(1:nomega-1),nomega-1)
!   call homo_grid_time(0d0,0.1d0, tau(-(ntau-1):(ntau-1)),wtau(-(ntau-1):(ntau-1)),ntau-1)
!   i = 0
!   do i_index =1, nomega, 10 
!    i = i+1
!    call homo_grid_freq(0.d0, 0.1d0, omega(1:10),womega(1:10),10)
!    call homo_grid_freq(0.1d0, 1.d0, omega(11:20),womega(11:20),10)
!    call homo_grid_freq(1.d0, 10.d0, omega(21:30),womega(21:30),10)
!    call homo_grid_freq(10.d0, 100.d0, omega(31:40),womega(31:40),10)
!    call homo_grid_freq(100.d0, 1000.d0, omega(41:100),womega(41:100),60)
    

     

!Q   enddo

!  tau(ntau) = omegamax*1.2
!  tau(-ntau) =- omegamax*1.2
!  wtau(ntau) = 0.d0
!    call homo_grid_freq(10.d0, 100.d0, omega(31:40),womega(31:40),10)
!  wtau(-ntau) = 0.d0

!  omega(nomega)   = omegamax * 1.2
!  omega(nomega-1) = omegamax * 1.1
!  omega(nomega-2) = omegamax * 1.15
!  womega(nomega)  = 0.d0
!  womega(nomega-1)= 0.d0
!  womega(nomega-2)= 0.d0
!  tau(ntau) = taumax*1.2
!  wtau(ntau) = 0.d0
!  tau(-ntau)= -taumax*1.2
!  wtau(-ntau)= 0.d0


! control output
  if(output) then
   if(myid.eq.0)then
    write(use_unit,'(2X,a,f7.3)') 'Time range....................................: ',taumax
    write(use_unit,'(2X,a,f7.3)') 'Frequency range...............................: ',omegamax
    write(use_unit,'(2X,a,i3)')   'Number of freuency points for self energy.....: ',ntau
    write(use_unit,'(2X,a,i3)')   'Number of frequency points....................: ',nomega
    write(use_unit,*) " "
   endif
  endif

end subroutine

!******
subroutine homo_grid_time(x1,x2,x,w,n)

implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(-n:n),w(-n:n)   ! abcsissa and weigths

! internal vars
real*8 h, t_0, a
integer :: i,j,m

  t_0 = x2 /n
  do i=1, n, 1
    x(i) = t_0*i
    x(-i) = -x(i)
  enddo
  w(:) = t_0

  return
end subroutine homo_grid_time
!******
subroutine homo_grid_freq(x1,x2,x,w,n)

implicit none

integer,          intent(in)  :: n           ! number of points
double precision, intent(in)  :: x1,x2       ! ranges
double precision, intent(out) :: x(1:n),w(1:n)   ! abcsissa and weigths

! internal vars
real*8 h, t_0, a
integer :: i,j,m

  t_0 = (x2-x1)/n
  x(1)=0.d0
  do i=2, n, 1
    x(i) = x1 + t_0*i
  enddo
  w(:) = t_0

  return
end subroutine homo_grid_freq



!****s*  FHI-aims/tf_ini
!  NAME
!    tf_ini
!  SYNOPSIS

subroutine tf_ini_scgw(ntau,nomega,taumax,omegamax,tau,omega,wtau,womega,output)

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

   use localorb_io, only: use_unit
   use mpi_tasks
   implicit none

!  ARGUMENTS
   
   integer,          intent(in)  :: ntau, nomega  
   double precision, intent(in)  :: taumax, omegamax
   double precision, intent(out) :: tau(-ntau:ntau), omega(nomega), wtau(-ntau:ntau), womega(nomega) 
   logical,          intent(in)  :: output

!counter 

   integer i   
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
  if(ntau<0)     stop 'ntau < 0'
  if(nomega<0)   stop 'nomega < 0'
  if(taumax<0)   stop 'taumax < 0'
  if(omegamax<0) stop 'omegamax < 0'

! calculate Gauss-Legendre grids and weights for omega and tau 

  call gauleg(0d0, omegamax, omega(1:nomega-3),womega(1:nomega-3),nomega-3)

!  call gauleg(0d0, omegamax, omega(1:nomega/2),womega(1:nomega/2),nomega/2)
!  call gauleg(omegamax, 5* omegamax, omega(nomega/2+1:nomega),womega(nomega/2+1:nomega),nomega/2)

!  call gauleg(0d0,taumax,tau(1:ntau-1),wtau(0:ntau-1),ntau-1)
  call gauleg(0d0,taumax,tau(1:ntau),wtau(1:ntau),ntau)

!  call gauleg(0d0,taumax,tau(1:ntau/2),wtau(1:ntau/2),ntau/2)
!  call gauleg(taumax,5*taumax,tau(ntau/2+1:ntau),wtau(ntau/2+1:ntau),ntau/2)
  
  womega(nomega-2)= 0
  womega(nomega-1)= 0
  womega(nomega)  = 0
!  womega(nomega-3)= 0

  omega(nomega-2) = 1.1 * omegamax
  omega(nomega-1) = 1.2 * omegamax
  omega(nomega)   = 1.3 * omegamax
 
!  wtau(nomega) = 0 
  tau (0) = 0
  wtau (0) = 0

  do i =1, ntau, 1
    tau(-i)= -tau(i)
    wtau(-i)=wtau(i)
  enddo 

  

  
!  call gauleg1(omegamax, omega(ntau:nomega),womega(ntau:nomega),nomega-ntau+1)
!  tau(ntau)        = 1.2*taumax  ! tail fit points
!  wtau(ntau)       = 0d0
!  omega(nomega)    = 1.2*omegamax
!  womega(nomega)   = 0d0

! control output
  if(output) then
   if(myid.eq.0) then
    write(use_unit,'(2X,a,f7.3)') 'Time range....................................: ',taumax
    write(use_unit,'(2X,a,f7.3)') 'Frequency range...............................: ',omegamax
    write(use_unit,'(2X,a,i3)')   'Number of freuency points for self energy.....: ',ntau
    write(use_unit,'(2X,a,i3)')   'Number of frequency points....................: ',nomega
   endif
  endif  

end subroutine
      end module scgw_grid
