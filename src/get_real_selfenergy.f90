!****s* FHI-aims/get_real_selfenergy
!  NAME
!     get_real_selfenergy
!  SYNOPSIS

      subroutine  get_real_selfenergy (flag_ana,n_freq,omega,x,npar,par,yexp)

!  PURPOSE
!     this subroutine gives the value of the fitting function evaluated 
!     at a given real argument x
!     the fitting function has a form of Pade approximation
!
!  USES     

      implicit none

!  ARGUMENTS

      integer :: flag_ana

      integer :: n_freq
      real*8  :: omega(n_freq)

      integer :: npar 
      complex*16  :: x
      complex*16 :: par(npar)

      complex*16  :: yexp

!  INPUTS
!  o flag_ana -- is an input integer number. 
!               If flag_ana = 0, use the expansion form
!                 P_N(x) = (a1 + a2*x + a3*x^2 + ...) /(1 + a4*x + a5*x^2 + ...)
!               If flag_ana = 1, use the expansion form
!                 P_N(x) = a1/(1 + a2(x-x[1])/(1 + a3(x-x[2])/(1 + ...)) )
!  o n_freq -- is an input integer number corresponding to the number of points
!               on the imagniary energy axis
!  o omega -- is an input real energy containing the energy points (Gauss-Legendre)
!               on the imagniary energy axis
!  o x -- is an complex input number. 
!  o npar -- is a positive integer input variable set to the number of parameters
!            in the fitting function
!  o par --  is a complex input array of length npar. On input par contains the 
!            values of the fitting parameter
!
!  OUTPUTS
!  o yexp --  is an output real number, correponding to the value of the
!             fitting function at x
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
!
!    Local variables

      complex*16 :: xtmp
      complex*16 :: gtmp
      complex*16 :: xdata(npar)
      complex*16 :: pade_numerator
      complex*16 :: pade_denominator

      integer :: n_step

!    counter
      integer i_par 
      integer i_dat

!    Start to work


! evaluate the value of the fitting function at abscisa x(i_dat) with given
! fitting parameters par(1:npar)

      xtmp = x

      if(flag_ana.eq.0) then

         pade_numerator = dcmplx(0.d0,0.d0)
         pade_denominator = dcmplx(0.d0,0.d0)
         if(mod(npar,2) .eq.0) then

            do i_par = npar/2, 2, -1 
              pade_numerator = (pade_numerator + par(i_par))*xtmp
            enddo

            do i_par = npar, npar/2+1, -1 
              pade_denominator = (pade_denominator + par(i_par))*xtmp
            enddo

         else

            do i_par = npar/2+1, 2, -1 
              pade_numerator = (pade_numerator + par(i_par))*xtmp
            enddo

            do i_par = npar, npar/2+2, -1 
              pade_denominator = (pade_denominator + par(i_par))*xtmp
            enddo

         endif

         pade_numerator = pade_numerator + par(1)
         pade_denominator = pade_denominator + dcmplx(1.d0,0.d0)

         yexp = pade_numerator/pade_denominator

      elseif(flag_ana.eq.1) then

         n_step = n_freq/(npar-1)
         i_dat = 1
         do i_par = 1, npar-1, 1
           xdata(i_par) = dcmplx(0.d0,omega(i_dat))
           i_dat = i_dat + n_step
         enddo
         xdata(npar) = dcmplx(0.d0,omega(n_freq))
         
         gtmp = dcmplx(1.d0,0.d0)
         do i_par = npar, 2, -1
           gtmp = 1.d0 + par(i_par)*(x-xdata(i_par-1))/gtmp 
         enddo
         yexp = par(1)/gtmp

      endif

      return

     end subroutine get_real_selfenergy
!***************
