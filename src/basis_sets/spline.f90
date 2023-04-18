!****h* FHI-aims/spline
!  NAME
!   spline
!  SYNOPSIS

      module spline

!  PURPOSE
!  The module includes routines for the polynomial spline.
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


    implicit none

!  Private declarations for more efficient allocation

!     "vector" and "matrix" are the objects to be handled by LAPACK.

      integer, private :: spline_dimension

      real*8, dimension(:), allocatable, private :: vector
      real*8, dimension(:), allocatable, private :: matrix_diag
      real*8, dimension(:), allocatable, private :: matrix_upper
      real*8, dimension(:), allocatable, private :: matrix_lower

      contains
!******
!-----------------------------------------------------------------------------------
!****s* spline/cubic_spline
!  NAME
!    cubic_spline
!     Algorithm as provided by Eric W. Weisstein, http://mathworld.wolfram.com/CubicSpline.html
!     VB, 11/14/2004: My implementation is a proof-of-concept only, for two reasons.
!     (ii) xmgrace gives better fitting spline polynomials than the functions provided here.
!     This is likely due to the Endpoint boundary conditions, which enforce zero second derivatives.
!     Can someone confirm that?
!
!     We return an array spl_param, which contains four polynomial parameters for each input point.
!     Thus, the input function now has an analytical shape, defined piece-wise between each point.
!
!     Matching conditions:
!
!      * derivatives of all spline pieces match at the input points
!      * second derivatives match at all input points
!      * Endpoint conditions: Second(!) derivatives are zero at each end point.
!
!     imported variables
!
!     n_points is the number of points on which to-be-splined function is given
!     f_grid   is the actual function to be splined, given on an evenly spaces grid n_points
!     spl_param is the array of spline function parameters which approximate our actual, continuous function f
!  SYNOPSIS

      subroutine cubic_spline ( f_grid, n_points, spl_param )

!  PURPOSE
!     Subroutine cubic_spline splines a function f_grid, given on grid points 1,2,3,...,n_points.
!
!  USES
    use localorb_io, only : use_unit
    use mpi_tasks, only : stderr, check_allocation
    implicit none
!  ARGUMENTS

      integer n_points
      real*8 f_grid(n_points)
      real*8 spl_param (4, n_points)

!  INPUTS
!   o n_points -- number of points
!   o f_grid -- values of data
!
!  OUTPUT
!   o spl_param -- spline parameters
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE






 
!     local variables

      integer i_info

!     counters

      integer i_point

!     begin work

      if (n_points.le.0) then
         write(stderr,'(1X,A)') &
        "* Cannot spline a function of zero grid points."
        stop
      else if (n_points.eq.1) then
!       return a constant
        spl_param(1,1) = f_grid(1)
        spl_param(2,1) = 0.d0
        spl_param(3,1) = 0.d0
        spl_param(4,1) = 0.d0
      else
!       general setup for spline

!       check whether spline dimension changed since last call
        if (n_points.ne.spline_dimension) then
!         set up spline matrix from scratch

!         "vector" and "matrix" are the objects to be handled by LAPACK.

          if (allocated(vector)) then
            deallocate(vector)
          end if
          if (allocated(matrix_diag)) then
            deallocate(matrix_diag)
          end if
          if (allocated(matrix_upper)) then
            deallocate(matrix_upper)
          end if
          if (allocated(matrix_lower)) then
            deallocate(matrix_lower)
          end if

          spline_dimension = n_points

          allocate(vector(spline_dimension),stat=i_info)
          call check_allocation(i_info, 'vector                        ') 

          allocate(matrix_diag(spline_dimension),stat=i_info)
          call check_allocation(i_info, 'matrix_diag                   ') 

          allocate(matrix_upper(spline_dimension-1),stat=i_info)
          call check_allocation(i_info, 'matrix_upper                  ') 

          allocate(matrix_lower(spline_dimension-1),stat=i_info)
          call check_allocation(i_info, 'matrix_lower                  ') 


        end if

!       eq. for first grid point

        vector(1) = 3d0*(f_grid(2) - f_grid(1))
        matrix_diag(1) = 2d0
        matrix_upper(1) = 1d0
        matrix_lower(1) = 1d0

!       eqs. for other points

        do i_point = 2, n_points-1, 1
          vector(i_point) = 3d0*(f_grid(i_point+1) - f_grid(i_point-1))
          matrix_diag(i_point) = 4d0
          matrix_upper(i_point) = 1d0
          matrix_lower(i_point) = 1d0
        enddo

!       eq. for last grid point

        vector(n_points) = 3d0*(f_grid(n_points) - f_grid(n_points-1))
        matrix_diag(n_points) = 2d0
!       no matrix_upper or matrix_lower elements.

!       solve tridiagonal system to obtain all spline derivatives, using lapack

        call dgtsv &
          ( n_points, 1, matrix_lower, matrix_diag, matrix_upper, &
            vector, n_points, i_info )

        if (i_info.ne.0) then
          write(use_unit,'(1X,A)') &
           "* spline.f : A cubic spline failed - investigate!"
          stop
        end if

!     now calculate parameters for spline polynomials

        do i_point = 1, n_points-1, 1

          spl_param(1,i_point) = f_grid(i_point)
          spl_param(2,i_point) = vector(i_point)
          spl_param(3,i_point) = 3d0*(f_grid(i_point+1)-f_grid(i_point)) &
            - 2d0*vector(i_point) - vector(i_point+1)
          spl_param(4,i_point) = 2d0*(f_grid(i_point)-f_grid(i_point+1)) &
            + vector(i_point) + vector(i_point+1)

        enddo

!       spline_vector_v2 needs value and derivative at the last point
        spl_param(1,n_points) = f_grid(n_points)
        spl_param(2,n_points) = vector(n_points)

      end if

      end subroutine cubic_spline

!******
!-----------------------------------------------------------------------------------
!****s* spline/cubic_spline_v2
!  NAME
!    cubic_spline_v2
!  SYNOPSIS

      subroutine cubic_spline_v2 ( spl_param, n_l_dim, n_coeff, n_grid_dim,  n_points, n_vector)

!  PURPOSE
!    Calculates polynomial spline coefficients for given interpolation points.
!    Algorithm etc. is the same as in cubic_spline, the distribution of
!    input/output data, however, is different.
!

!  USES
    use mpi_tasks, only : stderr
    implicit none
!  ARGUMENTS

      integer :: n_l_dim
      integer :: n_coeff
      integer :: n_grid_dim
      integer :: n_points
      integer :: n_vector
      real*8  :: spl_param(n_l_dim,n_coeff,n_grid_dim)


!  INPUTS
!   o n_l_dim -- maximum dimension of spl_param
!   o n_coeff -- number of spline coefficients, must be 2 or 4
!   o n_grid_dim -- dimension of splined vector
!   o n_points -- number of grid points
!   o n_vector -- dimension of splined function
!
!  INPUT/OUTPUT
!   o spl_param -- sline parameters
!        On input, spl_param(1:n_vector,1,1:n_points) must contain the values
!        of the interpolation points.
!        On output, spl_param(1:n_vector,2:n_coeff,1:n_points) is set
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!     local variables

   real*8 d_inv(20)
   integer i

   if (n_points.le.0) then
      write(stderr,'(1X,A)') &
      "* Cannot spline a function of zero grid points."
      stop
   endif

   if (n_points.eq.1) then
!     return a constant (although this doesn't make much sense for n_coeff==2)
      spl_param(1:n_vector,2,1) = 0.d0
      if(n_coeff==4) then
         spl_param(1:n_vector,3,1) = 0.d0
         spl_param(1:n_vector,4,1) = 0.d0
      endif
   endif

   ! Calculate inverse of the diagonal elements of the LR decomposition of the
   ! system matrix - after 20 elements these are (numerically) constant

   d_inv(1) = 0.5
   do i=2,20
      d_inv(i) = 1./(4-d_inv(i-1))
   enddo

   ! Calculate right hand side in spl_param(1:n_vector,2,1:n_points)
   ! and at the same time do forward elimination

   spl_param(1:n_vector,2,1) = 3*(spl_param(1:n_vector,1,2) - spl_param(1:n_vector,1,1))

   do i= 2, n_points-1
      spl_param(1:n_vector,2,i) = 3*(spl_param(1:n_vector,1,i+1) - spl_param(1:n_vector,1,i-1)) &
                                - d_inv(MIN(i-1,20))*spl_param(1:n_vector,2,i-1)
   enddo

   ! The same for the last equation
   spl_param(1:n_vector,2,n_points) = 3*(spl_param(1:n_vector,1,n_points) - spl_param(1:n_vector,1,n_points-1)) &
                                    - d_inv(MIN(n_points-1,20))*spl_param(1:n_vector,2,n_points-1)
   ! Start of backward elimination
   spl_param(1:n_vector,2,n_points) = spl_param(1:n_vector,2,n_points)/(2-d_inv(MIN(n_points-1,20)))

   ! Backward elimination
   do i=n_points-1,1,-1
      spl_param(1:n_vector,2,i) = (spl_param(1:n_vector,2,i) - spl_param(1:n_vector,2,i+1))*d_inv(MIN(i,20))
      if(n_coeff==4) then
          spl_param(1:n_vector,3,i) = 3d0*(spl_param(1:n_vector,1,i+1)-spl_param(1:n_vector,1,i)) &
                                    - 2d0*spl_param(1:n_vector,2,i) - spl_param(1:n_vector,2,i+1)
          spl_param(1:n_vector,4,i) = 2d0*(spl_param(1:n_vector,1,i)-spl_param(1:n_vector,1,i+1)) &
                                    +     spl_param(1:n_vector,2,i) + spl_param(1:n_vector,2,i+1)
      endif
   enddo

end subroutine cubic_spline_v2

!****s* spline/val_spline
!  NAME
!    val_spline
!  SYNOPSIS

     real*8 function val_spline ( r_output, spl_param, n_points )

!  PURPOSE
!     Function val_spline
!     takes a set of spline parameters spl_param for a function, given on grid points
!     1, ... n_points, and produce the splined interpolation onto some intermediate point
!
!  ARGUMENTS

     implicit none
      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters
!
!  OUTPUT
!    o val_spline -- value of a splined function in asked distance.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!     local variables

      integer i_spl
      real*8 t, term
      integer i_term

!     begin work

      i_spl = int(r_output)

      if ( i_spl.lt.1 ) then
        i_spl = 1
      else if ( i_spl.gt.(n_points-1) ) then
        i_spl = n_points-1
      end if

!      write(use_unit,*) "i_spl: ", i_spl

      t = r_output - dble(i_spl)

!      write(use_unit,*) "t: ", t

      val_spline = spl_param(1,i_spl)

      term = 1.0d0

      do i_term = 2, 4,1
        term = term * t

        val_spline = val_spline + term*spl_param(i_term,i_spl)

      enddo

      end function val_spline
!******
!-----------------------------------------------------------------------
!****s* spline/val_spline_deriv
!  NAME
!   val_spline_deriv
!  SYNOPSIS

      real*8 function val_spline_deriv ( r_output, spl_param, n_points )

!  PURPOSE
!     takes a set of spline parameters spl_param for a function, given on grid points
!     1, ... n_points, and produce the splined interpolation of derivative onto some intermediate point
!
!  USES
!  ARGUMENTS
      use localorb_io, only: use_unit
      implicit none
      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters
!
!  OUTPUT
!    o val_spline_deriv -- derivative of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables

      integer i_spl
      real*8 t, term
      integer i_term

!     begin work

      i_spl = int(r_output)

      if ( i_spl.lt.1 ) then
        i_spl = 1
      else if ( i_spl.gt.(n_points-1) ) then
        i_spl = n_points-1
      end if

!      write(use_unit,*) "i_spl: ", i_spl

      t = r_output - dble(i_spl)
      if (r_output .lt. 1) then
         write(use_unit,*) "warning: bad integration grid!!!!"
         write(use_unit,*) "r_output= ", r_output
      end if
!      write(use_unit,*) "t: ", t
      val_spline_deriv = spl_param(2,i_spl)

      term = 1.0d0

      do i_term = 3, 4,1
        term = term * t

        val_spline_deriv = val_spline_deriv + &
             (i_term - 1)*term*spl_param(i_term,i_spl)

      enddo
      end function




!******
!-----------------------------------------------------------------------
!****s* spline/val_spline_2nd_deriv
!  NAME
!   val_spline_2nd_deriv
!  SYNOPSIS

      real*8 function val_spline_2nd_deriv ( r_output, spl_param, n_points )

!  PURPOSE
!   2nd_deriv 
!
!  USES
!  ARGUMENTS
      use localorb_io, only: use_unit
      implicit none
      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters
!
!  OUTPUT
!    o val_spline_2nd_deriv -- derivative of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables

      integer i_spl
      real*8 t
      integer i_term

!     begin work

      i_spl = int(r_output)

      if ( i_spl.lt.1 ) then
        i_spl = 1
      else if ( i_spl.gt.(n_points-1) ) then
        i_spl = n_points-1
      end if

!      write(use_unit,*) "i_spl: ", i_spl

      t = r_output - dble(i_spl)
      if (r_output .lt. 1) then
         write(use_unit,*) "warning: bad integration grid!!!!"
         write(use_unit,*) "r_output= ", r_output
      end if
!      write(use_unit,*) "t: ", t
      val_spline_2nd_deriv = 2.0d0*spl_param(3,i_spl) + 6.0d0*t*spl_param(4,i_spl)

      end function








!******
!----------------------------------------------------------------------
!****s* spline/spline_vector
!  NAME
!   spline_vector
!  SYNOPSIS

      subroutine spline_vector &
      ( r_output, spl_param, n_grid_dim, n_l_dim, &
        n_points, n_vector, out_result )

!  PURPOSE
!    Subroutine spline_vector produces splined values for a whole array
!    of functions at once; this version should vectorize much better than
!    separately calling val_spline n_vector times
!

!  USES
!  ARGUMENTS

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_l_dim
      integer:: n_points
      integer:: n_vector
      real*8:: spl_param(n_l_dim,4,n_grid_dim)
      real*8:: out_result(n_vector)


!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters
!   o n_grid_dim -- dimension of splined vector
!   o n_l_dim -- maximum dimension of spl_param
!   o n_vector -- dimension of splined function
!
!  OUTPUT
!   o out_result -- values of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE





!C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term





!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

!ctest
!                write(use_unit,*) "inside spline_vector"
!                write(use_unit,*) i_spl
!                do i_index = 1,n_vector,1
!                  do i_term = 1,4,1
!                write(use_unit,*)  i_index, i_term,
!     +          spl_param(i_index,i_term,i_spl)
!                  enddo
!                enddo
!                stop
!ctest end


      t = r_output - dble(i_spl)

      out_result(1:n_vector) = spl_param(1:n_vector,1,i_spl)
      term = 1.d0
      do i_term = 2, 4, 1
        term = term * t

        out_result(1:n_vector) = out_result(1:n_vector) + &
          term * spl_param(1:n_vector,i_term,i_spl)

      enddo

!      out_result(:) =        spl_param(:,1,i_spl) +&
!     +                t    * spl_param(:,2,i_spl) +&
!     +                t**2 * spl_param(:,3,i_spl) +&
!     +                t**3 * spl_param(:,4,i_spl)
! 
!     out_result = 0.0d0

!      out_result(:n_vector) =        spl_param(:n_vector,2,i_spl) +&
!                     2*t    * spl_param(:n_vector,3,i_spl) +&
!                     3*t**2 * spl_param(:n_vector,4,i_spl) 
! !     +                t**3 * spl_param(:,4,i_spl)

      end subroutine spline_vector

!******
!-----------------------------------------------------------------------------------
!****s* spline/spline_deriv_vector
!  NAME
!   spline_deriv_vector
!  SYNOPSIS

      subroutine spline_deriv_vector &
           ( r_output, spl_param, n_grid_dim, n_l_dim, &
           n_points, n_vector, out_result )

!  PURPOSE
!    Subroutine  produces splined vderivatives for a whole array
!    of functions at once; this version should vectorize much better than
!    separately calling val_spline n_vector times
!
!  USES
!  ARGUMENTS

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_l_dim
      integer:: n_points
      integer:: n_vector
      real*8:: spl_param(n_l_dim,4,n_grid_dim)
      real*8:: out_result(n_vector)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters
!   o n_grid_dim -- dimension of splined vector
!   o n_l_dim -- maximum dimension of spl_param
!   o n_vector -- dimension of splined function
!
!  OUTPUT
!   o out_result -- derivative of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!     C     local variables

      integer :: i_spl
      real*8:: t



!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      

      out_result(:) =            spl_param(:,2,i_spl) + &
                      2 * t    * spl_param(:,3,i_spl) + &
                      3 * t**2 * spl_param(:,4,i_spl)

      end subroutine spline_deriv_vector

!  NAME
!   spline_deriv_vector
!  SYNOPSIS

      subroutine spline_2nd_deriv_vector &
           ( r_output, spl_param, n_grid_dim, n_l_dim, &
           n_points, n_vector, out_result )

!  PURPOSE
!    Subroutine  produces splined vderivatives for a whole array
!    of functions at once; this version should vectorize much better than
!    separately calling val_spline n_vector times
!
!  USES
!  ARGUMENTS

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_l_dim
      integer:: n_points
      integer:: n_vector
      real*8:: spl_param(n_l_dim,4,n_grid_dim)
      real*8:: out_result(n_vector)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters
!   o n_grid_dim -- dimension of splined vector
!   o n_l_dim -- maximum dimension of spl_param
!   o n_vector -- dimension of splined function
!
!  OUTPUT
!   o out_result -- derivative of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!     C     local variables

      integer :: i_spl
      real*8:: t



!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      out_result(:)  = 2 * spl_param(:,3,i_spl) + &
                      6 * t * spl_param(:,4,i_spl)


      end subroutine spline_2nd_deriv_vector
!******
!----------------------------------------------------------------------
!****s* spline/spline_vector_v2
!  NAME
!   spline_vector_v2
!  SYNOPSIS

      subroutine spline_vector_v2 &
      ( r_output, spl_param, &
        n_l_dim, n_coeff, n_grid_dim, &
        n_points, n_vector, out_result )

!  PURPOSE
!    Like spline_vector but with the possibility to use 2 or 4 spline coefficients
!

!  USES
!  ARGUMENTS

      real*8  :: r_output
      integer :: n_l_dim
      integer :: n_coeff
      integer :: n_grid_dim
      integer :: n_points
      integer :: n_vector
      real*8  :: spl_param(n_l_dim,n_coeff,n_grid_dim)
      real*8  :: out_result(n_vector)


!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o spl_param -- sline parameters
!   o n_l_dim -- maximum dimension of spl_param
!   o n_coeff -- number of spline coefficients, must be 2 or 4
!   o n_grid_dim -- dimension of splined vector
!   o n_points -- number of grid points
!   o n_vector -- dimension of splined function
!
!  OUTPUT
!   o out_result -- values of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!  local variables

      integer :: i_spl
      real*8  :: t, t2, t3, ta, tb, tc, td

!  begin work

      i_spl = int(r_output)
      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      if(n_coeff == 4) then

         ! 4 coeff version, i.e. a normal polynomial evaluation

         t2 = t*t
         t3 = t*t2

         out_result(1:n_vector) = &
              spl_param(1:n_vector,1,i_spl)    &
            + spl_param(1:n_vector,2,i_spl)*t  &
            + spl_param(1:n_vector,3,i_spl)*t2 &
            + spl_param(1:n_vector,4,i_spl)*t3

      else

         ! 2 coeff version, i.e. a Hermite polynom evaluation
         ! For the formula below, see http://en.wikipedia.org/wiki/Cubic_Hermite_spline

         ta = (t-1)*(t-1)*(1+2*t)
         tb = (t-1)*(t-1)*t
         tc = t*t*(3-2*t)
         td = t*t*(t-1)

         out_result(1:n_vector) = &
              spl_param(1:n_vector,1,i_spl  )*ta &
            + spl_param(1:n_vector,2,i_spl  )*tb &
            + spl_param(1:n_vector,1,i_spl+1)*tc &
            + spl_param(1:n_vector,2,i_spl+1)*td

      endif

      end subroutine spline_vector_v2

!******
!----------------------------------------------------------------------
!****s* spline/spline_deriv_vector_v2
!  NAME
!   spline_deriv_vector_v2
!  SYNOPSIS

      subroutine spline_deriv_vector_v2 &
      ( r_output, spl_param, &
        n_l_dim, n_coeff, n_grid_dim, &
        n_points, n_vector, out_result )

!  PURPOSE
!    Like spline_deriv_vector but with the possibility to use 2 or 4 spline coefficients
!

!  USES
!  ARGUMENTS

      real*8  :: r_output
      integer :: n_l_dim
      integer :: n_coeff
      integer :: n_grid_dim
      integer :: n_points
      integer :: n_vector
      real*8  :: spl_param(n_l_dim,n_coeff,n_grid_dim)
      real*8  :: out_result(n_vector)


!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o spl_param -- sline parameters
!   o n_l_dim -- maximum dimension of spl_param
!   o n_coeff -- number of spline coefficients, must be 2 or 4
!   o n_grid_dim -- dimension of splined vector
!   o n_points -- number of grid points
!   o n_vector -- dimension of splined function
!
!  OUTPUT
!   o out_result -- values of a splined function in asked distance.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!  local variables

      integer :: i_spl
      real*8  :: t, ta, tb, tc, td

!  begin work

      i_spl = int(r_output)
      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      if(n_coeff == 4) then

         ! 4 coeff version, i.e. a normal polynomial evaluation

         out_result(1:n_vector) = &
              spl_param(1:n_vector,2,i_spl)      &
            + spl_param(1:n_vector,3,i_spl)*2*t  &
            + spl_param(1:n_vector,4,i_spl)*3*t*t

      else

         ! 2 coeff version, i.e. a Hermite polynom evaluation

         ta = 6*t*(t-1)
         tb = (t-1)*(3*t-1)
         tc = -ta
         td = (3*t-2)*t

         out_result(1:n_vector) = &
              spl_param(1:n_vector,1,i_spl  )*ta &
            + spl_param(1:n_vector,2,i_spl  )*tb &
            + spl_param(1:n_vector,1,i_spl+1)*tc &
            + spl_param(1:n_vector,2,i_spl+1)*td

      endif

      end subroutine spline_deriv_vector_v2

!******
!-----------------------------------------------------------------------------------
!****s* spline/spline_deriv
!  NAME
!   spline_deriv
!  SYNOPSIS

      real*8 function spline_deriv ( r_output, spl_param, n_points )

!  PURPOSE
!     takes a set of spline parameters spl_param for a function, given on grid points
!     1, ... n_points, and produce the splined interpolation of derivative onto some intermediate point
!
!  USES
!  ARGUMENTS

      implicit none
      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters      real*8 r_output
!
!  OUTPUT
!   o spline_deriv -- derivative of a splined function in asked distance
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!     local variables

      integer i_spl
      real*8 t !, term
      !integer i_term


!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      spline_deriv  =            spl_param(2,i_spl) + &
                      2 * t    * spl_param(3,i_spl) + &
                      3 * t**2 * spl_param(4,i_spl)

      end function spline_deriv
!******
!------------------------------------------------------------------------------------------
!****s* spline/spline_2nd_deriv
!  NAME
!   spline_2nd_deriv
!  SYNOPSIS

      real*8 function spline_2nd_deriv ( r_output, spl_param, n_points )

!  PURPOSE
!     takes a set of spline parameters spl_param for a function, given on grid points
!     1, ... n_points, and produce the splined interpolation of second derivative onto some intermediate point
!
!  USES
!  ARGUMENTS

      implicit none
      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_points -- number of grid points
!   o spl_param -- sline parameters      real*8 r_output
!
!  OUTPUT
!   o spline_2nd_deriv -- derivative of a splined function in asked distance
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!     local variables

      integer i_spl
      real*8 t !, term
      !integer i_term


!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      spline_2nd_deriv  = 2 * spl_param(3,i_spl) + &
                      6 * t * spl_param(4,i_spl)

      end function spline_2nd_deriv
!******
!------------------------------------------------------------------------------------------
!****s* spline/spline_vector_waves
!  NAME
!    spline_vector_waves
!  SYNOPSIS

      subroutine spline_vector_waves &
      ( r_output, spl_param, n_grid_dim, n_compute_fns, &
        spline_start, spline_end, n_spl_points, n_spline, out_wave )

!  PURPOSE
!     takes a set of spline parameters spl_param for a function, given on grid points
!     1, ... n_points, and produce the splined interpolation of function onto some intermediate point
!     This subroutine is tailored to the spline of the waves.
!  USES
!  ARGUMENTS

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_compute_fns
      integer :: spline_start
      integer :: spline_end
      integer:: n_spl_points
      integer:: n_spline
      real*8:: spl_param(n_compute_fns,4,n_grid_dim)
      real*8:: out_wave(n_spline)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_spl_points -- number of grid points
!   o spl_param -- spline parameters      real*8 r_output
!   o n_grid_dim -- number of grid points in spl.
!   o n_compute_fns -- number of non-zero fns functions 
!   o spline_start -- starting point of the wave spl param
!   o spline_end -- ending point of the wave spl param
!   o n_spline -- number of splined functions
! 
!  OUTPUT
!   o out_wave -- splined value of wave
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE








!C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term





!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_spl_points-1, i_spl)

      t = r_output - dble(i_spl)

      out_wave(1:n_spline) = spl_param(spline_start:spline_end,1,i_spl)
      term = 1.d0
      do i_term = 2, 4, 1
        term = term * t

        out_wave(1:n_spline) = out_wave(1:n_spline) + &
             term * spl_param(spline_start:spline_end,i_term,i_spl)

      enddo

      end subroutine spline_vector_waves
!******
!---------------------------------------------------------------------------
!****s* spline/spline_vector_waves_deriv
!  NAME
!   spline_vector_waves_deriv
!  SYNOPSIS

      subroutine spline_vector_waves_deriv &
      ( r_output, spl_param, n_grid_dim, n_compute_fns, &
        spline_start, spline_end, n_spl_points, n_spline, out_wave )

!  PURPOSE
!     takes a set of spline parameters spl_param for a function, given on grid points
!     1, ... n_points, and produce the splined interpolation of derivative onto some intermediate point
!     This subroutine is tailored to the spline of the waves.
!
!  USES
!  ARGUMENTS


      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_compute_fns
      integer :: spline_start
      integer :: spline_end
      integer:: n_spl_points
      integer:: n_spline
      real*8:: spl_param(n_compute_fns,4,n_grid_dim)
      real*8:: out_wave(n_spline)

!  INPUTS
!   o r_output -- distance in units of logarithmic grid
!   o n_spl_points -- number of grid points
!   o spl_param -- sline parameters      real*8 r_output
!   o n_grid_dim -- number of grid points in spl.
!   o n_compute_fns -- number of non-zero fns functions 
!   o spline_start -- starting point of the wave spl param
!   o spline_end -- ending point of the wave spl param
!   o n_spline -- number of splined functions
! 
!  OUTPUT
!   o out_wave -- splined value of wave derivative
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE







!C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term





!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_spl_points-1, i_spl)

      t = r_output - dble(i_spl)

      out_wave(1:n_spline) = spl_param(spline_start:spline_end,2,i_spl)
      term = 1.d0
      do i_term = 3, 4, 1
        term = term * t

        out_wave(1:n_spline) = out_wave(1:n_spline) + &
             (i_term - 1) * term * &
             spl_param(spline_start:spline_end,i_term,i_spl)

      enddo

      end subroutine spline_vector_waves_deriv
!******
!-----------------------------------------------------------------------------------
!****s* spline/hermite
!  NAME
!   hermite
!  SYNOPSIS

      subroutine hermite &
           ( f_grid, fp_grid, n_points, spl_param )


!  PURPOSE
!  Subroutine hermite produces hermite interpolation parameters for any function,
!  using the same output format as a spline
!
!  USES
    use mpi_tasks, only : stderr
    implicit none
!  ARGUMENTS

        integer n_points
        real*8 f_grid(n_points)
        real*8 fp_grid(n_points)
        real*8 spl_param (4, n_points)

!  INPUTS
!   o n_points -- number of points
!   o f_grid -- ????????????
!   o fp_grid -- ????????????
!
!  OUTPUT
!   o spl_param -- spline parameters
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!       local variables

!       counters

        integer i_point

!        begin work

        if (n_points.le.0) then
          write(stderr,'(1X,A)') &
          "* Interpolation error in Hermite interpolation: "
          write(stderr,'(1X,A)') &
          "* Cannot interpolate a function of zero grid points."
          stop
        else if (n_points.eq.1) then
          ! return a constant
          spl_param(1,1) = f_grid(1)
          spl_param(2,1) = 0.d0
          spl_param(3,1) = 0.d0
          spl_param(4,1) = 0.d0
        else
          ! general case
          do i_point = 1, n_points-1, 1
            spl_param(1,i_point) = f_grid(i_point)
            spl_param(2,i_point) = fp_grid(i_point)
            spl_param(3,i_point) = &
              - 3 * f_grid(i_point) + 3 * f_grid(i_point+1) &
              - 2 * fp_grid(i_point) - fp_grid(i_point+1)
            spl_param(4,i_point) = &
              2 * f_grid(i_point) - 2 * f_grid(i_point+1) &
              + fp_grid(i_point) + fp_grid(i_point+1)
          enddo

          ! boundary conditions for extrapolation beyond n_points:
          ! keep second, third derivative constant ...
          spl_param(1,n_points) = f_grid(i_point)
          spl_param(2,n_points) = fp_grid(i_point)
          spl_param(3,n_points) = spl_param(3,n_points-1)
          spl_param(4,n_points) = spl_param(4,n_points-1)

        end if

        end subroutine hermite

!******
!------------------------------------------------------------------------------------------
!****s* spline/cleanup_spline
!  NAME
!   cleanup_spline
!  SYNOPSIS

      subroutine cleanup_spline (  )

!  PURPOSE
!  The subroutine deallocates module variables.
!
!  USES
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE






      implicit none
      spline_dimension = 0
          if (allocated(vector)) then
            deallocate(vector)
          end if
          if (allocated(matrix_diag)) then
            deallocate(matrix_diag)
          end if
          if (allocated(matrix_upper)) then
            deallocate(matrix_upper)
          end if
          if (allocated(matrix_lower)) then
            deallocate(matrix_lower)
          end if

      end subroutine cleanup_spline
!******
!-----------------------------------------------------------------------------------
      end module spline
