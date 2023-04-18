  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/determine_irreducible_k_grid
  !  NAME
  !    determine_irreducible_k_grid
  !  SYNOPSIS

  subroutine determine_irreducible_k_grid ()

    !  PURPOSE
    !
    !    Determine the irreducible k point wedge according to certain point group
    !    symmetry, which is also determined within this routine.
    !   
    !
    !  USES

    use constants
    use dimensions, only: n_irk_points, n_irk_points_task, n_k_points, &
        n_basis_fns, n_centers
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use pbc_lists
    use runtime_choices
    use geometry, only: lattice_vector, recip_lattice_vector
    use localorb_io, only: use_unit
    implicit none

    !  ARGUMENTS

    !  INPUTS
    !    o all in pbc_list module
    !  OUTPUTS
    !    o n_irk_points :: # of irreducible k points (at the moment only the inversion symmetry is considered)
    !    o irk_point_list :: maps the k-q grid point back to the 1st Brillouin zone
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales

    real*8, allocatable :: k_point_included(:,:)
    logical :: included_already, rot_included_already
    logical, parameter :: use_full_symmetry1=.false.
    integer, parameter :: n_max_symmetries=48
    integer :: n_symmetries
    real*8, allocatable :: symmats(:,:,:)
    real*8, allocatable :: origins(:,:)
    real*8, allocatable :: k_rot_list(:,:,:)
    integer :: info

    ! counter
    integer :: i_k_point
    integer :: i_irk_point
    integer :: i_periodic
    integer :: i_sym

    real*8 :: diff(3)
    real*8 :: summation(3)
    real*8 :: k_vec(3), k_rot_vec(3)
    character(*), parameter :: func = 'determine_irreducible_k_grid'

    allocate(symmats(3,3,n_max_symmetries),stat=info)
    call check_allocation(info,'symmats',func)
    allocate(origins(3,n_max_symmetries),stat=info)
    call check_allocation(info,'origins',func)
    allocate(k_rot_list(3,n_max_symmetries,n_k_points),stat=info)
    call check_allocation(info,'k_rot_list',func)

    call get_system_symmetries(n_max_symmetries, n_symmetries, symmats, origins)
!    if(myid.eq.0) then
!      write(use_unit,*) "n_symmetries :", n_symmetries
!      do i_sym = 1, n_symmetries, 1
!        write(use_unit,*) "i_sym : ", i_sym
!        do i_periodic = 1, 3
!          write(use_unit,'(3f12.3)') symmats(:,i_periodic,i_sym) 
!        enddo
!        write(use_unit,*)
!        write(use_unit,'(3f12.3)') origins(:,i_sym) 
!        write(use_unit,*)
!      enddo
!    endif
    
    do i_k_point = 1, n_k_points, 1
      k_vec = matmul(recip_lattice_vector,k_point_list(i_k_point,:))
      do i_sym = 1, n_symmetries, 1
       k_rot_vec = matmul(symmats(:,:,i_sym),k_vec)/2.d0/pi
       k_rot_list(:,i_sym,i_k_point) = mod(matmul(k_rot_vec,lattice_vector),1.d0)
!       if(myid.eq.0) then
!        write(use_unit,'(I3,3f16.3)') i_sym, k_rot_list(:)
!      endif
      enddo
    enddo

    allocate(k_point_included(n_k_points,3),stat=info)
    call check_allocation(info,'k_point_included',func)
    n_irk_points=0
    do i_k_point = 1, n_k_points, 1
         included_already = .false.
         rot_included_already = .false.
         do i_irk_point = 1, n_irk_points, 1

            diff(:) = k_point_included(i_irk_point,:) - k_point_list(i_k_point,:)

            if (mod(abs(diff(1)),1.d0) .lt. 1.e-6 .and. &
                mod(abs(diff(2)),1.d0) .lt. 1.e-6 .and. &
                mod(abs(diff(3)),1.d0) .lt. 1.e-6 ) then
                 included_already = .true.
                 exit
            endif

            if (.not. included_already) then
!  this is a different k_point, and check if it is connected to other k points by inversion
               summation(:) = k_point_included(i_irk_point,:) + k_point_list(i_k_point,:)
 
               if (mod(abs(summation(1)),1.d0) .lt. 1.e-6 .and. &
                   mod(abs(summation(2)),1.d0) .lt. 1.e-6 .and. &
                   mod(abs(summation(3)),1.d0) .lt. 1.e-6 ) then
                     included_already = .true.
                     exit
               endif
               
            endif

            if(.not.included_already .and. use_full_symmetry) then
              do i_sym = 1, n_symmetries, 1
                 diff(:) = k_point_included(i_irk_point,:) - k_rot_list(:,i_sym,i_k_point)
                 if (mod(abs(diff(1)),1.d0) .lt. 1.e-6 .and. &
                     mod(abs(diff(2)),1.d0) .lt. 1.e-6 .and. &
                     mod(abs(diff(3)),1.d0) .lt. 1.e-6 ) then
                     rot_included_already = .true.
                     exit
                 endif
              enddo
            endif

            if(rot_included_already) then
              included_already = .true.
              exit
            endif
!  end of loop over i_irk_point
         enddo

         if(.not. included_already ) then
            n_irk_points = n_irk_points + 1
            k_point_included(n_irk_points,:) = k_point_list(i_k_point,:)
            irk_point_mapping(i_k_point) = n_irk_points
         else
            irk_point_mapping(i_k_point) = i_irk_point
         endif

!  end of loop over i_k_point
   enddo

   allocate(irk_point_list(n_irk_points,3))
   irk_point_list(:,:) = k_point_included(1:n_irk_points,:)

   n_irk_points_task = 0
   do i_irk_point = 1, n_irk_points, 1
       if(myid.eq.mod(i_irk_point,n_tasks) .and. myid .le. n_irk_points) then
          n_irk_points_task = n_irk_points_task + 1
       endif
   enddo
!   n_irk_points_task = max (n_irk_points_task,1)

   deallocate(k_point_included) 

! defined in "pbc_list.f90", determine the inverse mapping between the irreducible k grid and the full k grid
   allocate(inv_irk_point_mapping(n_irk_points),stat=info)
   call check_allocation(info,'inv_irk_point_mapping',func)

   do i_irk_point = 1, n_irk_points, 1
      do i_k_point = 1, n_k_points, 1
         if(i_irk_point  .eq. irk_point_mapping(i_k_point)) then
            inv_irk_point_mapping(i_irk_point) = i_k_point
            exit
         endif
      enddo
   enddo
   
! defined in "pbc_list.f90", tell if a given k point included in the irreducible set or not
   allocate(irk_point_included(n_k_points),stat=info)
   call check_allocation(info,'irk_point_included',func)

   do i_k_point = 1, n_k_points, 1
      irk_point_included(i_k_point) = .false.
      do i_irk_point = 1, n_irk_points, 1
         if (i_k_point .eq. inv_irk_point_mapping(i_irk_point) ) then
           irk_point_included(i_k_point) = .true.
         endif
      enddo
   enddo

   if(.not.allocated(irk_weight)) then
     allocate(irk_weight(n_irk_points))
   endif
   irk_weight(:)=0.d0
!    write(use_unit,*) "n_irk_points :", n_irk_points
   do i_k_point = 1, n_k_points, 1
     i_irk_point = irk_point_mapping(i_k_point)
     irk_weight(i_irk_point)=irk_weight(i_irk_point)+1
!     write(use_unit,*) i_k_point, i_irk_point, inv_irk_point_mapping(i_irk_point), irk_point_included(i_k_point)
   enddo
   irk_weight(:)=irk_weight(:)/n_k_points
   
   if(myid.eq.0) then
     write(use_unit,*) " Irreducible k points # ", n_irk_points
   endif
   if(allocated(symmats)) then
     deallocate(symmats)
   endif
   if(allocated(origins)) then
     deallocate(origins)
   endif
   if(allocated(k_rot_list)) then
     deallocate(k_rot_list)
   endif
   
  
   end subroutine determine_irreducible_k_grid
