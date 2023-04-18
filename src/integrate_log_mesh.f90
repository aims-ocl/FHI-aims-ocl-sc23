!----------------------------------------------------------------------
!
!  Integration functions / subroutines on the logarithmic mesh
!  * function int_log_mesh: Calculates Int (u_1(r) * u_2(r) dr)
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
      real*8 function int_log_mesh &
      ( wave_1, wave_2, n_grid, r_grid &
      )

      implicit none

!  imported variables

      integer n_grid
      real*8 wave_1(n_grid)
      real*8 wave_2(n_grid)
      real*8 r_grid(n_grid)

!  local variables

      real*8 alpha

!  counters

      integer i_grid

!  begin work

!     assume logarithmic grid
      if (n_grid.gt.1) then
        alpha = log(r_grid(2)/r_grid(1))
        int_log_mesh = 0.
        do i_grid = 1, n_grid, 1
          int_log_mesh = int_log_mesh + &
          alpha * r_grid(i_grid) * wave_1(i_grid) * wave_2 (i_grid)
        enddo
      else
        int_log_mesh = 1.
      end if

!  that's all folks

      return
      end

!----------------------------------------------------------------------
