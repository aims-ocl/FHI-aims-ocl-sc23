!****s* FHI-aims/get_pade_approx_para
!  NAME
!     get_pade_approx_para
!  SYNOPSIS

      subroutine  get_pade_approx_para (n, x, y, npar, par)

!  PURPOSE
!     this subroutine get the parameter of Pade approximation in the form of
!     continued fraction P_N(x) = a1/(1 + a2(x-x1)/(1 + a3(x-x2)/(1 + ... ))).
!     Here N is the order of the Pade approximaton equivalent to number of 
!     parameters npar (i.e. N=npar). 
!
!  USES     
      use localorb_io,only :use_unit
      implicit none 

!  ARGUMENTS

      integer n 
      integer npar 
      real*8  x(n)
      complex*16  y(n)

      complex*16  par(npar)

!  INPUTS
!  o n -- is a positive integer input variable set to be the number of data points
!  o x -- is a compelx input array of length n. On input x contains the abscisas 
!         of the data points
!  o y -- is a complex input array of length n. On input contains is the values of 
!         the data points
!  o npar -- is a positive integer input variable set to the number of parameters
!         in the fitting function
!
!  OUTPUTS
!  o par --  is a output complex array of length npar. On output par contains the 
!         initally estimated fitting parameters.
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

      integer :: n_step 
      complex*16 :: g_func(npar, npar)
      complex*16 :: xtmp(npar), ytmp(npar)

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
      if(npar.gt.n) then

        write(use_unit,*) &
            "* Warning! Pade approximation has more parameters than ", &
            "the data point. Set the number of parameters as ", &
            "available data points"
        npar=n

      endif

      n_step = n/(npar-1)

      i_dat = 1
      do i_par = 1, npar-1, 1
        xtmp(i_par) = dcmplx(0.d0,x(i_dat))
        ytmp(i_par) = y(i_dat) 
        i_dat = i_dat + n_step  
      enddo
      xtmp(npar)=dcmplx(0.d0,x(n))
      ytmp(npar)=y(n)
         
!    evaluate the parameter using the recursion relation  
!    g(n,x[n]) = (g(n-1,x[n-1])-g(n-1,x[n]))/((x[n-1]-x)g(n-1,x[n])),
!    and a_n = g(n,x[n]), g(1,x[n]) = y[n].
!    
!    this recursion relation really comes from the definition of
!    g(n,x), namely
!
!    P_N(x) = g(1,x) = g(1,x[1])/((1+g(2,x)(x-x[1])) 
!                    = g(1,x[1])/((1+g(2,x[2])(x-x[1])/(1+g(3,x)(x-x[2]))) 
!                    = ...
!   

     g_func(:,1) = ytmp(:)
     do i_par = 2, npar, 1
       do i_dat = i_par, npar, 1
         g_func(i_dat, i_par) = &
           (g_func(i_par-1,i_par-1)-g_func(i_dat,i_par-1))/ &
           ((xtmp(i_dat)-xtmp(i_par-1))*g_func(i_dat,i_par-1))
       enddo
     enddo

     do i_par = 1, npar, 1
       par(i_par) = g_func(i_par,i_par)
!       write(use_unit,'(I4,2f16.6)') i_par, par(i_par)
     enddo
!       write(use_unit,*)

     return

     end subroutine get_pade_approx_para
!************
