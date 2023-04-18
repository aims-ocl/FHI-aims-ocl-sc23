  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/determine_k_minus_q_list
  !  NAME
  !    determine_k_minus_q_list
  !  SYNOPSIS

  subroutine determine_k_minus_q_list &
            (kq_point_list, kpq_point_list)

    !  PURPOSE
    !
    !    mapping the k_minus_q and k_plus_q reciprocal grid point back to the 1st BZ.
    !   
    !
    !  USES

    use pbc_lists
    use runtime_choices
    use localorb_io
    use mpi_tasks, only: aims_stop
    use dimensions, only: n_k_points, n_irk_points_task
    implicit none

    !  ARGUMENTS

       integer, intent(OUT) :: kq_point_list(n_k_points,n_k_points)
       integer, intent(OUT) :: kpq_point_list(n_k_points,n_k_points)

    !  INPUTS
    !    o all in pbc_list module
    !  OUTPUTS
    !    o kq_point_list :: maps the k-minus-q grid point back to the 1st Brillouin zone
    !    o kpq_point_list :: maps the k-plus-q grid point back to the 1st Brillouin zone
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
    real*8 :: k_diff_thr
    character*140 :: info_str
    character(*), parameter :: func = 'determine_k_minus_q_list'

    k_diff_thr=min(dot_product(k_points_offset,k_points_offset),1.d-5)+1.d-10

    kq_point_list = 0
    kpq_point_list = 0
    do i_q_point = 1, n_k_points, 1
      do i_k_point = 1, n_k_points, 1

!         write(use_unit,*) i_q_point, i_k_point
         k_minus_q(:)=k_point_list(i_k_point,:)-k_point_list(i_q_point,:)
         k_plus_q(:)=k_point_list(i_k_point,:)+k_point_list(i_q_point,:)

! convert the k-q point back to the 1st Brillouin zone
         do i_periodic=1, 3, 1
           if(k_minus_q(i_periodic).lt.-1d-10) then
             k_minus_q(i_periodic) = k_minus_q(i_periodic)+1.d0
           endif
           if(k_plus_q(i_periodic).ge.1.d0) then
             k_plus_q(i_periodic) = k_plus_q(i_periodic)-1.d0
           endif
         enddo
!         write(use_unit,*) "k_minus_q new:", k_minus_q(:)

! check which point in the 1st BZ that the k-q point corresponds to
         do i_kq_point = 1, n_k_points, 1
            diff(:) = k_minus_q(:) - k_point_list(i_kq_point,:)
            if(dot_product(diff,diff).lt.k_diff_thr) then
              kq_point_list(i_k_point,i_q_point) = i_kq_point
              exit
            endif
         enddo

         do i_kq_point = 1, n_k_points, 1
            diff(:) = k_plus_q(:) - k_point_list(i_kq_point,:)
            if(dot_product(diff,diff).lt.k_diff_thr) then
              kpq_point_list(i_k_point,i_q_point) = i_kq_point
              exit
              endif
         enddo
!         write(use_unit,*) "i_kq_point", i_kq_point,kq_point_list(i_k_point,i_q_point),kpq_point_list(i_k_point,i_q_point),i_k_point, i_q_point
         if(kq_point_list(i_k_point,i_q_point).le.0 .or. kpq_point_list(i_k_point,i_q_point).le.0 ) then

           if(n_k_points .eq. 1 .and. dot_product(k_point_list(1,:),k_point_list(1,:)).gt.1.e-12) then
              kq_point_list(1,1) = 1
              write(info_str, '(2X,A)') &
               & "Using single k point slightly away from the Gamma point."
              call localorb_info(info_str)  
 
           else
             call aims_stop(" Error in determining the k-q mesh, stop!", func)
           endif

         endif

      enddo
   enddo

!Dirty hack waiting for a final fix
   if (n_k_points.eq.1) then
      kq_point_list=1
      kpq_point_list=1
   end if

  end subroutine determine_k_minus_q_list
