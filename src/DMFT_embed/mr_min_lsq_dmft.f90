!****s* FHI-aims/mr_min_lsq
!  NAME
!     mr_min_lsq
!  SYNOPSIS
      subroutine  mr_min_lsq_dmft (n, x, y, npar, par, i_k_point)

     use dimensions
     use localorb_io, only: use_unit
     use mpi_tasks 
!
!  PURPOSE
!     this subroutine produces the least square fitting to a data set (x,y),
!     using Levenberg-Marquardt algorithm.
!
!  USES     

      implicit none

!  ARGUMENTS

      integer n, i_k_point 
      integer npar 
      complex*16  x(n)
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
!  o par --  is a output real array of length npar. On output par contains the 
!         final fitting parameters.
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
     
      complex*16 :: yexp(n) 
      complex*16 :: dydp(npar,n) 
      complex*16 :: dy 
      complex*16 :: dydp1, dydp2 


      real*8 :: coeff_matr(2*npar,2*npar) 
      real*8 :: adapted_coeff_matr(2*npar,2*npar) 
      real*8 :: c_vect(2*npar) 
      real*8 :: delta_par(2*npar) 

      complex*16 :: par_exp(npar) 

      real*8  :: chisq
      real*8  :: chisq_old
      real*8  :: lamda

      integer :: ipiv(2*npar), info 

      logical :: converg
      logical :: accepted

!    Parameters 
       
      real*8  :: fit_acc
!      parameter (fit_acc = 1.d-15)

      integer  :: n_max_cycle
      parameter (n_max_cycle = 300)

!   counter
      integer i_dat
      integer i_par, j_par
      integer i_cycle

!    Start to work

! initialize the parameter to get started
!if(myid.eq.0) write(use_unit,*) 'k point nummer', i_k_point

      if(use_scgw .or. use_scgw0 )then 
         fit_acc = 1.d-10
      elseif(use_dmft_gw)then 
         fit_acc = 1.d-10
      else 
         fit_acc = 1.d-15
      endif

      call initialize_fitting_parameter_dmft  &
             (n, x, y, npar, par)   

      call evaluate_fitting_function  &
           (n, x, npar, par, yexp, dydp)   

! calculate the initial chi squre
      chisq_old = 0.d0
      do i_dat = 1, n
         dy = y(i_dat)-yexp(i_dat)
         chisq_old = chisq_old + conjg(dy)*dy
      enddo 
!      chisq_old = chisq_old/dble(n)
      
         
!      write(use_unit,*) 
!      write(use_unit,*) "  | initial chisq: ", chisq_old
!      write(use_unit,*) 

      lamda = 0.001
      i_cycle = 0
      converg = .false.
      do while (.not.converg)

         i_cycle = i_cycle + 1

!         write(use_unit,*) 
!         write(use_unit,*) "  | Cycle: ", i_cycle
!         write(use_unit,*) 

         coeff_matr(:,:) = 0.d0
         c_vect(:) = 0.d0
         do i_dat = 1, n, 1

!  construct the coefficient matrix and column vector for the linear equation
!   to be solved in Gauss-Newton method (or Levenberg-Marquardt method).

           dy = y(i_dat) - yexp(i_dat)
           do i_par = 1, 2*npar, 1

             if(i_par .le. npar) then
                dydp1 = dydp(i_par,i_dat)
             else
                dydp1 = dydp(i_par-npar,i_dat)*dcmplx(0.d0,1.d0)
             endif

             do j_par = 1, 2*npar, 1

               if(j_par .le. npar) then
                  dydp2 = dydp(j_par,i_dat)
               else
                  dydp2 = dydp(j_par-npar,i_dat)*dcmplx(0.d0,1.d0)
               endif

               coeff_matr(j_par,i_par) = coeff_matr(j_par,i_par) + &
                                         dreal( conjg(dydp1)*dydp2 )
             enddo

             c_vect(i_par) = c_vect(i_par) + dreal( conjg(dy)*dydp1 ) 

            enddo
          enddo

          accepted = .false. 
          do while (.not. accepted)

!  adapt the coefficient matrix (by enhancing the diagonal element) for 
!  the Levenberg-Marquardt method

            adapted_coeff_matr(:,:) = coeff_matr(:,:)
            do i_par = 1, 2*npar, 1
               adapted_coeff_matr(i_par,i_par) = coeff_matr(i_par,i_par) *  &
                                                 (1.d0+lamda)
            enddo

! solve the complex linear equation, on entry delta_par(:) = c_vect
! on exit, delta_par(:) gives  the change of the parameters
            delta_par(:) = c_vect(:)
            call dgesv(2*npar,1,adapted_coeff_matr,2*npar,ipiv,delta_par,2*npar,info)

! update the the parameters
            do i_par = 1, npar 
              par_exp(i_par) = par(i_par)+ dcmplx(delta_par(i_par),  &
                                                  delta_par(i_par+npar))
            enddo

! calculate the new value of the fitting function with updated parameters
            call evaluate_fitting_function  &
                 (n, x, npar, par_exp, yexp, dydp)   

! calculate the new chi square
            chisq = 0.d0
            do i_dat = 1, n
              dy = y(i_dat)-yexp(i_dat)
              chisq = chisq + conjg(dy)*dy
            enddo 
!            chisq = chisq/dble(n)

!            write(use_unit,*) "  | chisq =  ", chisq
!            write(use_unit,*) 

            if(chisq .gt. chisq_old) then
               lamda = lamda * 1.d1
            else
               accepted = .true.
               lamda = lamda * 1.d-1
               par(:) = par_exp(:)
            endif


! end of do while (accepted)
          enddo

          if(abs(chisq-chisq_old) .lt. fit_acc) then

            converg = .true.
!            write(use_unit,'(A,e16.3,A,I6,A)') "  Successful fitting with accuaray ", &
!                 chisq,  " is achieved after ", i_cycle,  " iterations."

          elseif (i_cycle .gt. n_max_cycle) then
            converg = .true.

             !write(use_unit,'(A,e16.3,A,I6,A)') "  Self-energy fitting with accuary  ", &



     
if(myid.eq.0) then           !write(use_unit,'(A,e16.3,A,I6,A)') "  Self-energy fitting with accuary  ", &
            write(use_unit,*) "  Self-energy fitting with accuary  ", &
                 fit_acc,  " is NOT achieved after ", i_cycle,  " iterations." ,"at k_point", i_k_point
            write(use_unit,'(A,e16.3)') " The accuracy achieved currently is ", &
                 abs(chisq-chisq_old)  
endif
          endif 
          chisq_old = chisq

      enddo

!      do  i_dat = 1, n, 1
!       write(use_unit,'(I6,6f16.6)') i_dat, x(i_dat), &
!       y(i_dat), yexp(i_dat)
!       write(use_unit,*)
!       do i_par = 1, npar
!        write(use_unit,'(I6,2f16.6)')i_par, par(i_par)
!       enddo
!      enddo
!      write(use_unit,*)
      end subroutine mr_min_lsq_dmft
