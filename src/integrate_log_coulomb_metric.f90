!----------------------------------------------------------------------
!
!  Integration functions / subroutines on the logarithmic mesh
!  * function int_log_coulomb_metric: Calculates Int (u_1(r)*v(r-r')*u_2(r')drdr')
!
!
!  Because we do not know the limiting behavior of the functions which we
!  integrate, we simply integrate sum(i) dr/di * u_1(i) * u_2(i)
!  FIXME: This is not a great way, this would work much better on the
!  radial integration mesh!
!  See also the better integration in dftseq(), which takes the limiting
!  behavior of wave functions into account explicitly
!
!----------------------------------------------------------------------
!
      real*8 function int_log_coulomb_metric &
      (i_l, wave_1, wave_2, n_grid, r_grid &
      )

      use  constants
      implicit none

!  imported variables

      integer i_l
      integer n_grid
      real*8 wave_1(n_grid)
      real*8 wave_2(n_grid)
      real*8 r_grid(n_grid)

! input
! o i_l    - angular momentum 
! o n_grid - logrithmic grid point number
! o r_grid - logrithmic grid point 
! o wave_1 - first radail wave
! o wave_2 - second radail wave
!
!  local variables

      real*8 alpha
      real*8 prefactor
      real*8 integral_zero_r
      real*8 integral_r_infinity
      real*8 integral_v_times_wave2(n_grid)

!  counters

      integer i_grid

!  begin work

!     assume logarithmic grid
      if (n_grid.gt.1) then

        alpha = log(r_grid(2)/r_grid(1))
        prefactor = alpha*pi4/(2*i_l + 1.d0)

        integral_zero_r = 0.d0
        do i_grid = 1, n_grid, 1
          integral_zero_r = integral_zero_r + &
            r_grid(i_grid)**2 * wave_2(i_grid)*r_grid(i_grid)**i_l

          integral_v_times_wave2(i_grid) =  &
              integral_zero_r * prefactor/r_grid(i_grid)**(i_l+1)
        enddo

        integral_r_infinity = 0.d0
        do i_grid = n_grid, 1, -1
          integral_v_times_wave2(i_grid)= integral_v_times_wave2(i_grid) + &
              integral_r_infinity * prefactor * r_grid(i_grid)**i_l

          integral_r_infinity = integral_r_infinity + &
            r_grid(i_grid)**2 * wave_2(i_grid)/r_grid(i_grid)**(i_l+1)
        enddo

        int_log_coulomb_metric = 0.
        do i_grid = 1, n_grid, 1
          int_log_coulomb_metric = int_log_coulomb_metric + &
          alpha * r_grid(i_grid)**2 * wave_1(i_grid) * &
             integral_v_times_wave2(i_grid)
        enddo
      else
        int_log_coulomb_metric = 1.
      end if

!  that's all folks

      return
      end

!----------------------------------------------------------------------
