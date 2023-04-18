!****s* FHI-aims/initialize_fitting_parameter
!  NAME
!     initialize_fitting_parameter
!  SYNOPSIS

      subroutine  initialize_fitting_parameter_dmft (n, x, y, npar, par)
!
!  PURPOSE
!     this subroutine evaluate the initial fitting parameter using Pade
!     approximation P_N(x) = A_N(x)/B_N(x). Here N is the order of the 
!     Pade approximaton equivalent to number of parameters npar (i.e. N=npar). 
!
!  USES

      use localorb_io, only: use_unit
      implicit none 

!  ARGUMENTS

      integer n 
      integer npar 
      complex*16  x(n)
      complex*16  y(n)

      complex*16  par(npar)

!  INPUTS
!  o n -- is a positive integer input variable set to be the number of data points
!
!  o x -- is a compelx input array of length n. On input x contains the abscisas 
!         of the data points
!  o y -- is a complex input array of length n. On input contains is the values of the 
!         data points
!  o npar --  is a positive integer input variable set to the number of parameters
!          in the fitting function
!
!  OUTPUTS
!  o par -- is a output complex array of length npar. On output par contains the 
!          initally estimated fitting parameters.
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

      integer :: n_step 
      complex*16 :: coeff_matr(npar, npar)
      complex*16 :: xtmp(npar)

      integer :: ipiv(npar)
      integer :: info

      character*140 :: info_str

!    counter
      integer i_par, j_par 
      integer i_dat

!    Start to work

!    here we select npar many data points for determing the coefficents for
!    Pade approximation
!     
      n_step = n/npar

      i_dat = 1
      do i_par = 1, npar, 1
        xtmp(i_par) = x(i_dat) 
        par(i_par) = y(i_dat) 
        i_dat =  i_dat + n_step  
!        i_dat =  i_dat + 1
      enddo
         
!    Construct the coefficient matrix
!    here A(x) = par(1) + par(2)*x + par(3)*x**2 + ...
!         B(x) = 1 + par(N/2+1)*x + par(N/2+2)*x**2 + ...
!    note here N = npar 

     if(mod(npar,2) .eq. 0) then 
        do i_par = 1, npar, 1
          coeff_matr(i_par,1) = 1
          coeff_matr(i_par,npar/2+1) = -xtmp(i_par)*par(i_par) 
          do j_par = 2, npar/2, 1
            coeff_matr(i_par,j_par) = xtmp(i_par)*coeff_matr(i_par,j_par-1) 
            coeff_matr(i_par,j_par+npar/2) = xtmp(i_par)* &
                                  coeff_matr(i_par,j_par+npar/2-1) 
          enddo
        enddo
     else
        do i_par = 1, npar, 1
          coeff_matr(i_par,1) = 1
          coeff_matr(i_par,npar/2+2) = -xtmp(i_par)*par(i_par) 
          do j_par = 2, npar/2, 1
            coeff_matr(i_par,j_par) = xtmp(i_par)*coeff_matr(i_par,j_par-1) 
            coeff_matr(i_par,j_par+npar/2+1) = xtmp(i_par)* &
                                       coeff_matr(i_par,j_par+npar/2) 
          enddo
          coeff_matr(i_par,npar/2+1) = xtmp(i_par)*coeff_matr(i_par,npar/2) 
        enddo
     endif

!   Solve the complex linear equation 
     call zgesv(npar, 1, coeff_matr, npar, ipiv, par, npar, info)

     if(info.lt.0) then
        write(use_unit,'(2X,A,A)') &
       " Warning! Initialization of the fitting parameter has a problem:", &
       " error in solving the linear equation."
       write(use_unit,*) "info=", info
     endif

     end subroutine initialize_fitting_parameter_dmft
!***************
!***************
