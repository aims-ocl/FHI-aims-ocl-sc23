  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/determine_k_minus_q_mixed_grid
  !  NAME
  !    determine_k_minus_q_mixed_grid
  !  SYNOPSIS

  subroutine determine_k_minus_q_mixed_grid &
            (n_q_points, current_k_point, q_point_list, k_minus_q_list)

    !  PURPOSE
    !
    !  For a given k point, and a list of q point, determine the k-q point
    !  list 
    !   
    !
    !  USES
       
    implicit none

    !  ARGUMENTS

    integer :: n_q_points
    real*8, dimension(3) :: current_k_point
    real*8, dimension(n_q_points, 3) :: q_point_list
    real*8, dimension(n_q_points, 3) :: k_minus_q_list
    !  INPUTS
    !    o all in pbc_list module
    !  OUTPUTS
    !    o kq_point_list :: maps the k-q grid point back to the 1st Brillouin zone
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales

    real*8, allocatable :: k_minus_q_mixed(:,:)
    real*8 :: k_grid_point(3)
    real*8 :: k_minus_q(3)
    real*8 :: diff(3)
    character(*), parameter :: func = 'determine_k_minus_q_mixed_grid'

    integer :: n_mixed_points
    integer :: info
    integer :: n_special_k_points

!  counter
    integer :: i_kq_point
    integer :: i_k_point, i_q_point
    integer :: i_periodic
    integer :: i_band
    integer :: i_k_point_special

    n_special_k_points = 0
    do i_band = 1, n_plot_band, 1
      n_special_k_points = n_special_k_points +  n_points_in_band(i_band)
    enddo

    allocate(k_minus_q_mixed(3,n_q_points*n_special_k_points),stat=info)
    call check_allocation(info,'k_minus_q_mixed',func)


    n_mixed_points = 0
    i_k_point_special = 0
    do i_band =1, n_plot_band, 1

      n_k_points = n_points_in_band(i_band)
      do i_k_point = 1, n_k_points, 1

        i_k_point_special = i_k_point_special + 1
!FIXME: better to swap the index in "band_begin" and "band_end"
        k_grid_point(:) =  band_begin(i_band,:) +  real(i_k_point-1)/real(n_k_points-1) &
                *( band_end(i_band,:) -  band_begin(i_band,:))

        do i_q_point = 1, n_q_points, 1

          n_mixed_points = n_mixed_points + 1
!         write(use_unit,*) i_q_point, i_k_point
          k_minus_q(:)=k_grid_point(:)-k_point_list(i_q_point,:)
!         write(use_unit,*) "k_minus_q:", k_minus_q(:)

! convert the k-q point back to the 1st Brillouin zone
          do i_periodic=1, 3, 1
            if(k_minus_q(i_periodic).lt.-5d-1) then
              k_minus_q(i_periodic) = k_minus_q(i_periodic)+1.d0
            elseif(k_minus_q(i_periodic).gt.5d-1) then
              k_minus_q(i_periodic) = k_minus_q(i_periodic)-1.d0
            endif
          enddo
!         write(use_unit,*) "k_minus_q new:", k_minus_q(:)

! check which point in the 1st BZ that the k-q point corresponds to
          do i_kq_point = 1, n_mixed_points-1, 1
             diff(:) = k_minus_q(:) - k_point_list(i_kq_point,:)
             if(dot_product(diff,diff).lt.1.e-12) then
               n_mixed_points = n_mixed_points - 1 
             else
               k_minus_q_mixed(:,n_mixed_points) = k_minus_q(:)
             endif
          enddo

! end of loop over i_q_point
        enddo
! end of loop over i_k_point
      enddo
! end of loop over i_band
   enddo
   write(use_unit,*) "k_minus_q_points real : ", n_mixed_points
   write(use_unit,*) "k_minus_q_points total : ", n_special_k_points*n_q_points, n_special_k_points, n_q_points
  
   deallocate(k_minus_q_mixed)

  end subroutine determine_k_minus_q_mixed_grid
