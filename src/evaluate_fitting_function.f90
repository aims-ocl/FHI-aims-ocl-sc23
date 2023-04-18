!****s* FHI-aims/evaluate_fitting_function
!  NAME
!     evaluate_fitting_function
!  SYNOPSIS

      subroutine  evaluate_fitting_function (n, x, npar, par, yexp, dydp)

!  PURPOSE
!     this subroutine evaluate the fitting function at abscisas x(n) with
!     give fitting parameters par(npar), yielding the experimental fitting
!     value yexp(n), and its derivative with respect to the fitting parameters
!     dydp = d(yexp)/d(par)
!
!  USES

      implicit none

!  ARGUMENTS

      integer n 
      integer npar 
      complex*16  x(n)
      complex*16  yexp(n)

      complex*16  par(npar)
      complex*16  dydp(npar,n)

!  INPUTS
!
!  o n -- is a positive integer input variable set to be the number of data points
!
!  o x -- is a compelx input array of length n. On input x contains the abscisas 
!         of the data points
!
!  o npar -- is a positive integer input variable set to the number of parameters
!            in the fitting function
!
!  o par -- is a complex input array of length npar. On input par contains the 
!            values of the fitting parameter
!
!
!  OUTPUTS
!
!  o yexp -- is an output complex array of length n. On output yexp contains the 
!             initally the value of the fitting function evaluated at abscias x(n)
!
!  o dydp --  is an output complex array of dimension (n,npar). On output dydp 
!             contains the derivative of the fitting function yexp(n) with respect to
!             the fitting parameter par(npar)
!
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

!    Local variables

      complex*16 :: pade_numerator
      complex*16 :: pade_denominator

      integer :: info


!    counter
      integer i_par, j_par 
      integer i_dat

!    Start to work

      do i_dat = 1, n, 1

! evaluate the value of the fitting function at abscisa x(i_dat) with given
! fitting parameters par(1:npar)

        pade_numerator = dcmplx(0.d0,0.d0)
        pade_denominator = dcmplx(0.d0,0.d0)
        if(mod(npar,2) .eq.0) then
          do i_par = npar/2, 2, -1 
            pade_numerator = (pade_numerator + par(i_par))*x(i_dat)
          enddo
          do i_par = npar, npar/2+1, -1 
            pade_denominator = (pade_denominator + par(i_par))*x(i_dat)
          enddo
        else
          do i_par = npar/2+1, 2, -1 
            pade_numerator = (pade_numerator + par(i_par))*x(i_dat)
          enddo
          do i_par = npar, npar/2+2, -1 
            pade_denominator = (pade_denominator + par(i_par))*x(i_dat)
          enddo
        endif
        pade_numerator = pade_numerator + par(1)
        pade_denominator = pade_denominator + dcmplx(1.d0,0.d0)

        yexp(i_dat) = pade_numerator/pade_denominator

! evaluate the derivative of the fitting function at abscisa x(i_dat) with respect to
! fitting parameters par(1:npar)

        if(mod(npar,2).eq.0) then

          dydp(1,i_dat) = dcmplx(1.d0,0.d0)/pade_denominator
          dydp(npar/2+1,i_dat) = -yexp(i_dat)*x(i_dat)/pade_denominator
          do i_par = 2, npar/2, 1
            dydp(i_par,i_dat) = x(i_dat)*dydp(i_par-1,i_dat)
            dydp(i_par+npar/2,i_dat) = x(i_dat)*dydp(i_par+npar/2-1,i_dat)
          enddo

        else

          dydp(1,i_dat) = dcmplx(1.d0,0.d0)/pade_denominator
          dydp(npar/2+2,i_dat) = -yexp(i_dat)*x(i_dat)/pade_denominator
          do i_par = 2, npar/2, 1
            dydp(i_par,i_dat) = x(i_dat)*dydp(i_par-1,i_dat)
            dydp(i_par+npar/2+1,i_dat) = x(i_dat)*dydp(i_par+npar/2,i_dat)
          enddo
          dydp(npar/2+1,i_dat) = x(i_dat)*dydp(npar/2,i_dat)

        endif

!  end of loop over x(i_dat)
      enddo

     end subroutine evaluate_fitting_function
!******
!****s* FHI-aims/evaluate_fitting_function
!  NAME
!     evaluate_fitting_function
!  SYNOPSIS

      subroutine  evaluate_fitting_function_W (n, x, npar, par, yexp, dydp)

!  PURPOSE
!     this subroutine evaluate the fitting function at abscisas x(n) with
!     give fitting parameters par(npar), yielding the experimental fitting
!     value yexp(n), and its derivative with respect to the fitting parameters
!     dydp = d(yexp)/d(par)
!
!  USES

      implicit none

!  ARGUMENTS

      integer n 
      integer npar 
      real*8  x(n)
      real*8  yexp(n)

      real*8  par(npar)
      real*8  dydp(npar,n)

!  INPUTS
!
!  o n -- is a positive integer input variable set to be the number of data points
!
!  o x -- is a compelx input array of length n. On input x contains the abscisas 
!         of the data points
!
!  o npar -- is a positive integer input variable set to the number of parameters
!            in the fitting function
!
!  o par -- is a complex input array of length npar. On input par contains the 
!            values of the fitting parameter
!
!
!  OUTPUTS
!
!  o yexp -- is an output complex array of length n. On output yexp contains the 
!             initally the value of the fitting function evaluated at abscias x(n)
!
!  o dydp --  is an output complex array of dimension (n,npar). On output dydp 
!             contains the derivative of the fitting function yexp(n) with respect to
!             the fitting parameter par(npar)
!
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

!    Local variables

      real*8 :: pade_numerator
      real*8 :: pade_denominator

      integer :: info


!    counter
      integer i_par, j_par 
      integer i_dat

!    Start to work

      do i_dat = 1, n, 1

        yexp(i_dat) = par(1)/(par(2)+x(i_dat)) +&
                      par(3)/(par(4)+x(i_dat))

        dydp(1,i_dat) = 1.d0/(par(2)+x(i_dat))
        dydp(2,i_dat) = - par(1)*dydp(1,i_dat)**2
        dydp(3,i_dat) = 1.d0/(par(4)+x(i_dat))
        dydp(4,i_dat) = - par(3)*dydp(1,i_dat)**2

      enddo

      end subroutine evaluate_fitting_function_W
!******
