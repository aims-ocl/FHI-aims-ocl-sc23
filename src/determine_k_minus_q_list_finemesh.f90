  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/determine_k_minus_q_list_finemesh
  !  NAME
  !    determine_k_minus_q_list_finemesh
  !  SYNOPSIS

  subroutine determine_k_minus_q_list_finemesh &
            (n_k_points_old, k_point_list_old, kq_point_list)

    !  PURPOSE
    !
    !    mapping the k_minus_q on a fine k grid back to the original grid.
    !   
    !
    !  USES

    use dimensions, only: n_k_points
    use pbc_lists
    use mpi_tasks, only: aims_stop
    use runtime_choices
    implicit none

    !  ARGUMENTS

       integer, intent(IN) :: n_k_points_old
       real*8, intent(IN) :: k_point_list_old(n_k_points_old,3)
       integer, intent(OUT) :: kq_point_list(n_k_points,n_k_points)

    !  INPUTS
    !    o n_k_points_old :: number of the k points in the orginal mesh
    !    o k_point_list_old :: original k grid
    !  OUTPUTS
    !    o kq_point_list :: maps the k-minus-q grid point back to the 1st Brillouin zone
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales

    integer :: i_kq_point
    integer :: i_k_point, i_q_point
    integer :: i_periodic

    real*8 :: k_minus_q(3)
    real*8 :: k_plus_q(3)
    real*8 :: diff(3)
    character(*), parameter :: func = 'determine_k_minus_q_list_finemesh'

! Note that now n_k_points = n_k_points_xyz(1)*n_k_points_xyz(2)*n_k_points_xyz(3)*enhancement_factor

    kq_point_list = 0
    do i_q_point = 1, n_k_points, 1
      do i_k_point = 1, n_k_points, 1

!         write(use_unit,*) i_q_point, i_k_point
         k_minus_q(:)=k_point_list(i_k_point,:)-k_point_list(i_q_point,:)

! convert the k-q point back to the 1st Brillouin zone
         do i_periodic=1, 3, 1
           if(k_minus_q(i_periodic).lt.-1d-10) then
             k_minus_q(i_periodic) = k_minus_q(i_periodic)+1.d0
           endif
         enddo
!         write(use_unit,*) "k_minus_q new:", k_minus_q(:)

! check which point in the 1st BZ that the k-q point corresponds to
         do i_kq_point = 1, n_k_points_old, 1
            diff(:) = k_minus_q(:) - k_point_list_old(i_kq_point,:)
            if(dot_product(diff,diff).lt.1.e-12) then
              kq_point_list(i_k_point,i_q_point) = i_kq_point
              exit
            endif
         enddo

!         write(use_unit,*) "i_kq_point", i_kq_point
         if(kq_point_list(i_k_point,i_q_point).lt.0) then
          call aims_stop(" Error in determining the k-q in fine k mesh, stop!", func)
         endif

      enddo
   enddo
  

  end subroutine determine_k_minus_q_list_finemesh
