      module spline

C  Private declarations for more efficient allocation

C     "vector" and "matrix" are the objects to be handled by LAPACK.

      integer, save, private :: spline_dimension = 0

      real*8, dimension(:), allocatable, save, private :: vector
      real*8, dimension(:), allocatable, save, private :: matrix_diag
      real*8, dimension(:), allocatable, save, private :: matrix_upper
      real*8, dimension(:), allocatable, save, private :: matrix_lower

      contains
C-----------------------------------------------------------------------------------
C     Subroutine cubic_spline splines a function f_grid, given on grid points 1,2,3,...,n_points.

      subroutine cubic_spline ( f_grid, n_points, spl_param )

      implicit none

C     Algorithm as provided by Eric W. Weisstein, http://mathworld.wolfram.com/CubicSpline.html
C     VB, 11/14/2004: My implementation is a proof-of-concept only, for two reasons.
C     (ii) xmgrace gives better fitting spline polynomials than the functions provided here.
C          This is likely due to the Endpoint boundary conditions, which enforce zero second derivatives.
c          Can someone confirm that?

C     We return an array spl_param, which contains four polynomial parameters for each input point.
C     Thus, the input function now has an analytical shape, defined piece-wise between each point.

C     Matching conditions: 

C     * derivatives of all spline pieces match at the input points
C     * second derivatives match at all input points
C     * Endpoint conditions: Second(!) derivatives are zero at each end point.

C     imported variables

C     n_points is the number of points on which to-be-splined function is given
C     f_grid   is the actual function to be splined, given on an evenly spaces grid n_points
C     spl_param is the array of spline function parameters which approximate our actual, continuous function f

C     input

      integer n_points
      real*8 f_grid(n_points)

C     output

      real*8 spl_param (4, n_points)

C     local variables

      integer i_info

C     counters

      integer i_point

C     begin work

      if (n_points.le.0) then
         write(6,'(1X,A)') 
     +  "* Cannot spline a function of zero grid points."
        stop
      else if (n_points.eq.1) then
C       return a constant
        spl_param(1,1) = f_grid(1)
        spl_param(2,1) = 0.d0
        spl_param(3,1) = 0.d0
        spl_param(4,1) = 0.d0
      else
C       general setup for spline

C       check whether spline dimension changed since last call
        if (n_points.ne.spline_dimension) then
C         set up spline matrix from scratch

C         "vector" and "matrix" are the objects to be handled by LAPACK.

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

          allocate(vector(spline_dimension))
          allocate(matrix_diag(spline_dimension))
          allocate(matrix_upper(spline_dimension-1))
          allocate(matrix_lower(spline_dimension-1))

        end if

C       eq. for first grid point

        vector(1) = 3d0*(f_grid(2) - f_grid(1))
        matrix_diag(1) = 2d0
        matrix_upper(1) = 1d0
        matrix_lower(1) = 1d0

C       eqs. for other points

        do i_point = 2, n_points-1, 1
          vector(i_point) = 3d0*(f_grid(i_point+1) - f_grid(i_point-1))
          matrix_diag(i_point) = 4d0
          matrix_upper(i_point) = 1d0
          matrix_lower(i_point) = 1d0
        enddo

C       eq. for last grid point

        vector(n_points) = 3d0*(f_grid(n_points) - f_grid(n_points-1))
        matrix_diag(n_points) = 2d0
C       no matrix_upper or matrix_lower elements.

C       solve tridiagonal system to obtain all spline derivatives, using lapack

        call dgtsv
     +    ( n_points, 1, matrix_lower, matrix_diag, matrix_upper, 
     +      vector, n_points, i_info )

        if (i_info.ne.0) then
          write (6,'(1X,A)') 
     +     "* spline.f : A cubic spline failed - investigate!"
          stop
        end if

C     now calculate parameters for spline polynomials

        do i_point = 1, n_points-1, 1

          spl_param(1,i_point) = f_grid(i_point)
          spl_param(2,i_point) = vector(i_point)
          spl_param(3,i_point) = 3d0*(f_grid(i_point+1)-f_grid(i_point))
     +      - 2d0*vector(i_point) - vector(i_point+1)
          spl_param(4,i_point) = 2d0*(f_grid(i_point)-f_grid(i_point+1))
     +      + vector(i_point) + vector(i_point+1)

        enddo

C       for beauty's sake only, val_spline does not use this point ...
        spl_param(1,n_points) = f_grid(n_points)

      end if

      end subroutine cubic_spline

C-----------------------------------------------------------------------------------
C     Function val_spline
C     takes a set of spline parameters spl_param for a function, given on grid points
C     1, ... n_points, and produce the splined interpolation onto some intermediate point

      real*8 function val_spline ( r_output, spl_param, n_points )

      implicit none

!     Input variables

      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

C     local variables

      integer i_spl
      real*8 t, term
      integer i_term

C     begin work

      i_spl = int(r_output)

      if ( i_spl.lt.1 ) then
        i_spl = 1
      else if ( i_spl.gt.(n_points-1) ) then
        i_spl = n_points-1
      end if

c      write(6,*) "i_spl: ", i_spl

      t = r_output - dble(i_spl)

c      write (6,*) "t: ", t

      val_spline = spl_param(1,i_spl) 

      term = 1.0d0

      do i_term = 2, 4,1
        term = term * t

        val_spline = val_spline + term*spl_param(i_term,i_spl) 

      enddo

      end function val_spline

!-----------------------------------------------------------------------

      real*8 function val_spline_deriv ( r_output, spl_param, n_points )

      implicit none

!     Input variables

      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

C     local variables

      integer i_spl
      real*8 t, term
      integer i_term

C     begin work

      i_spl = int(r_output)

      if ( i_spl.lt.1 ) then
        i_spl = 1
      else if ( i_spl.gt.(n_points-1) ) then
        i_spl = n_points-1
      end if

c      write(6,*) "i_spl: ", i_spl

      t = r_output - dble(i_spl)
      if (r_output .lt. 1) then
         write (6,*) "warning: bad integration grid!!!!"
         write (6,*) "r_output= ", r_output
      end if
c      write (6,*) "t: ", t
      val_spline_deriv = spl_param(2,i_spl) 

      term = 1.0d0

      do i_term = 3, 4,1
        term = term * t

        val_spline_deriv = val_spline_deriv + 
     +       (i_term - 1)*term*spl_param(i_term,i_spl) 

      enddo
      end function

!----------------------------------------------------------------------
C    Subroutine spline_vector produces splined values for a whole array
C    of functions at once; this version should vectorize much better than 
C    separately calling val_spline n_vector times 
C
      subroutine spline_vector
     +( r_output, spl_param, n_grid_dim, n_l_dim, 
     +  n_points, n_vector, out_result )

!      implicit none

!     Input variables

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_l_dim
      integer:: n_points
      integer:: n_vector
      real*8:: spl_param(n_l_dim,4,n_grid_dim)
      real*8:: out_result(n_vector)

!C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term

!test
      integer ::  i_index
!test end

!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

!ctest
c                write(6,*) "inside spline_vector"
c                write(6,*) i_spl
c                do i_index = 1,n_vector,1
c                  do i_term = 1,4,1
c                write(6,*)  i_index, i_term,
c     +          spl_param(i_index,i_term,i_spl)
c                  enddo
c                enddo
c                stop
!ctest end


      t = r_output - dble(i_spl)

      out_result(1:n_vector) = spl_param(1:n_vector,1,i_spl) 
      term = 1.d0
      do i_term = 2, 4, 1
        term = term * t

        out_result(1:n_vector) = out_result(1:n_vector) + 
     +    term * spl_param(1:n_vector,i_term,i_spl)

      enddo

c      out_result(:) =        spl_param(:,1,i_spl) + 
c     +                t    * spl_param(:,2,i_spl) + 
c     +                t**2 * spl_param(:,3,i_spl) + 
c     +                t**3 * spl_param(:,4,i_spl)

      end subroutine spline_vector
C-----------------------------------------------------------------------------------

      subroutine spline_deriv_vector
     +     ( r_output, spl_param, n_grid_dim, n_l_dim, 
     +     n_points, n_vector, out_result )
      
!     implicit none
      
!     Input variables
      
      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_l_dim
      integer:: n_points
      integer:: n_vector
      real*8:: spl_param(n_l_dim,4,n_grid_dim)
      real*8:: out_result(n_vector)
      
!     C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term

!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      out_result(:) =            spl_param(:,2,i_spl) + 
     +                2 * t    * spl_param(:,3,i_spl) + 
     +                3 * t**2 * spl_param(:,4,i_spl)

      end subroutine spline_deriv_vector
C-----------------------------------------------------------------------------------

      real*8 function spline_deriv ( r_output, spl_param, n_points )

      implicit none

!     Input variables

      real*8 r_output
      integer n_points
      real*8 spl_param(4,n_points)

C     local variables

      integer i_spl
      real*8 t !, term
      !integer i_term


!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_points-1, i_spl)

      t = r_output - dble(i_spl)

      spline_deriv  =            spl_param(2,i_spl) + 
     +                2 * t    * spl_param(3,i_spl) + 
     +                3 * t**2 * spl_param(4,i_spl)

      end function spline_deriv

!------------------------------------------------------------------------------------------


      subroutine spline_vector_waves
     +( r_output, spl_param, n_grid_dim, n_compute_fns, 
     +  spline_start, spline_end, n_spl_points, n_spline, out_wave )

!      implicit none

!     Input variables

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_compute_fns
      integer :: spline_start
      integer :: spline_end
      integer:: n_spl_points
      integer:: n_spline
      real*8:: spl_param(n_compute_fns,4,n_grid_dim)
      real*8:: out_wave(n_spline)

!C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term

!test
      integer ::  i_index
!test end

!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_spl_points-1, i_spl)

      t = r_output - dble(i_spl)

      out_wave(1:n_spline) = spl_param(spline_start:spline_end,1,i_spl) 
      term = 1.d0
      do i_term = 2, 4, 1
        term = term * t

        out_wave(1:n_spline) = out_wave(1:n_spline) + 
     +       term * spl_param(spline_start:spline_end,i_term,i_spl)

      enddo

      end subroutine spline_vector_waves

      subroutine spline_vector_waves_deriv
     +( r_output, spl_param, n_grid_dim, n_compute_fns, 
     +  spline_start, spline_end, n_spl_points, n_spline, out_wave )

!      implicit none

!     Input variables

      real*8:: r_output
      integer :: n_grid_dim
      integer :: n_compute_fns
      integer :: spline_start
      integer :: spline_end
      integer:: n_spl_points
      integer:: n_spline
      real*8:: spl_param(n_compute_fns,4,n_grid_dim)
      real*8:: out_wave(n_spline)

!C     local variables

      integer :: i_spl

      real*8:: t, term
      integer:: i_term

!test
      integer ::  i_index
!test end

!C     begin work

      i_spl = int(r_output)

      i_spl = max(1,i_spl)
      i_spl = min(n_spl_points-1, i_spl)

      t = r_output - dble(i_spl)

      out_wave(1:n_spline) = spl_param(spline_start:spline_end,2,i_spl) 
      term = 1.d0
      do i_term = 3, 4, 1
        term = term * t

        out_wave(1:n_spline) = out_wave(1:n_spline) + 
     +       (i_term - 1) * term * 
     +       spl_param(spline_start:spline_end,i_term,i_spl)

      enddo

      end subroutine spline_vector_waves_deriv
C-----------------------------------------------------------------------------------
!  Subroutine hermite produces hermite interpolation parameters for any function,
!  using the same output format as a spline

        subroutine hermite
     +   ( f_grid, fp_grid, n_points, spl_param )

!       input

        integer n_points
        real*8 f_grid(n_points)
        real*8 fp_grid(n_points)

!       output

        real*8 spl_param (4, n_points)

!       local variables

!       counters

        integer i_point

!        begin work

        if (n_points.le.0) then
          write(6,'(1X,A)') 
     +    "* Interpolation error in Hermite interpolation: "
          write(6,'(1X,A)') 
     +    "* Cannot interpolate a function of zero grid points."
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
            spl_param(3,i_point) = 
     +        - 3 * f_grid(i_point) + 3 * f_grid(i_point+1)
     +        - 2 * fp_grid(i_point) - fp_grid(i_point+1)
            spl_param(4,i_point) = 
     +        2 * f_grid(i_point) - 2 * f_grid(i_point+1)
     +        + fp_grid(i_point) + fp_grid(i_point+1)
          enddo
          
          ! boundary conditions for extrapolation beyond n_points:
          ! keep second, third derivative constant ...
          spl_param(1,n_points) = f_grid(i_point)
          spl_param(2,n_points) = fp_grid(i_point)
          spl_param(3,n_points) = spl_param(3,n_points-1)
          spl_param(4,n_points) = spl_param(4,n_points-1) 

        end if

        end subroutine hermite
!------------------------------------------------------------------------------------------
C  deallocate stuff

      subroutine cleanup_spline
     +  (  )

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
C-----------------------------------------------------------------------------------
      end module spline
