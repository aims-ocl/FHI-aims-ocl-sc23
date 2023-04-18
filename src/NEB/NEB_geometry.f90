!
!
!  subroutines that do the geometry evaluation/update in the NEB
!  contains:
!  subroutine get_difference_magnitude
!  subroutine get_difference_image
!  subroutine BFGS_n_coords
!

!---------------------------------------------------------------------------------------------------
! subroutine to calculate the absolute values |r_i - r_i-1| for all images. 
! returns an array of magnitudes delta_R(0:n_images) for which 
! delta_R(0) = |r_start - r_1|
! delta_R(1) = |r_1 - r_2|
! ... 
subroutine get_difference_magnitude(n_atoms,n_images,start_coords,coords,end_coords,delta_R)
  implicit none
  integer :: n_atoms, n_images, i_coords, i_images
  real*8 :: ddot
  real*8, dimension(3*n_atoms) :: start_coords, end_coords, buf
  real*8, dimension(3*n_atoms,n_images) :: coords
  real*8, dimension(0:n_images) :: delta_R
  delta_R = 0d0
  buf(:) = start_coords(:) - coords(:,1)
  delta_R(0) = ddot(3*n_atoms,buf,1,buf,1)
  do i_images = 2, n_images
     buf(:)              = coords(:,i_images-1)-coords(:,i_images)
     delta_R(i_images-1) = ddot(3*n_atoms,buf,1,buf,1)
  end do
  buf(:)            = coords(:,n_images)-end_coords(:)
  delta_R(n_images) = ddot(3*n_atoms,buf,1,buf,1)
  delta_R(:) = sqrt(delta_R(:))
end subroutine get_difference_magnitude


!---------------------------------------------------------------------------------------------------
! calculates the absolute and RMS differences between the the stored and current i_image^th image
! helps when deciding whether or not the transition state has converged!
subroutine get_difference_image(n_atoms,n_images,i_image,coords,stored_coords,diff_max,diff_RMS)
  implicit none
  integer :: n_atoms, n_images, i_image, i_coords
  real*8 :: diff_max, diff_RMS
  real*8, dimension(3*n_atoms,n_images) :: coords, stored_coords
  diff_max = 0d0
  diff_RMS = 0d0
  do i_coords = 1, 3*n_atoms
     diff_max = max(diff_max,abs(coords(i_coords,i_image)-stored_coords(i_coords,i_image)))
     diff_RMS = diff_RMS + (coords(i_coords,i_image)-stored_coords(i_coords,i_image))**2d0
  end do
  diff_RMS = sqrt(diff_RMS)/dble(3*n_atoms)
end subroutine get_difference_image

!---------------------------------------------------------------------------------------------------
!   Broyden-Fletcher-Goldfarb-Shanno method, straight from Wikipedia
!   Line Search for the minimum along a given search direction may be done, 
!   but algorithm may also be run with a simple quadratic extrapolation (no line search)
!
!         - this version does not know anything about atoms, it only 'knows' about n_coords 
!           coordinates and does its thing within that space, whatever n_coords might describe.
!         - implements a line search for the initial convergence when using steepest-descent
subroutine BFGS_n_coords(n_coords,                  &
                         coords,                    &
                         total_forces,              &
                         ini_relaxation,            &
                         initial_hessian,           &
                         stored_forces,             &
                         stored_coords,             &
                         search_direction,          &
                         hessian,                   &
                         max_atomic_move,           &
                         min_line_step,             &
                         line_step_reduce,          &
                         line_step,                 &
                         object_function_exists,    &
                         object_function,           &
                         stored_object_function,    &
                         object_function_tolerance, &
                         line_search_on,            &
                         force_times_search_dir,    &
                         max_line_step_increase)
  implicit none
  ! imported variables
  integer n_coords
  real*8, dimension(n_coords) :: coords
  real*8, dimension(n_coords) :: total_forces
  real*8, dimension(n_coords) :: stored_forces
  real*8, dimension(n_coords) :: stored_coords
  real*8, dimension(n_coords) :: search_direction
  real*8, dimension(n_coords, n_coords) :: hessian, hessian_work
  real*8 :: max_atomic_move, object_function, stored_object_function, object_function_tolerance, max_line_step_increase
  logical :: ini_relaxation,initial_hessian, object_function_exists, line_search_on
  ! Local variables
  real*8, dimension(n_coords) :: force_difference
  real*8 :: forcediff_times_searchdir
  real*8, dimension(n_coords) :: Hs
  real*8 :: sHs, line_step
  real*8 :: force_times_search_dir, force_times_search_dir_old
  real*8 :: max_displacement
  real*8 :: limit_move
  real*8 :: min_line_step
  real*8 :: line_step_reduce
  real*8, external :: ddot
  ! equation solver workspace arrays etc
  integer, dimension(n_coords) :: ipiv
  integer :: lwork, info
  real*8, dimension(:), allocatable :: work
  
  ! counters
  integer :: i_coord
  integer :: i_atom
  integer :: i_mobile_atom
  integer :: i_hessian
  ! begin work
  
  write(use_unit,'(2X,A)')"Advancing geometry using BFGS."
  
  ! Determine whether relaxation step was successful, current
  ! coordinates, Hessian, search direction, and line step
  if (ini_relaxation) then
     
     ! This is the initial relaxation step - must store energy, forces etc.
     ! Store energy, coordinates, forces of present iteration for next iteration
     stored_forces          = total_forces
     stored_coords          = coords
     stored_object_function = object_function
     line_step              = 1d0
     search_direction(:)    = total_forces(:)            
     
     ! reinitialize hessian to be unity ...
     ! This means the first step will be a simple steepest descent.
     hessian = 0.d0
     do i_hessian = 1, n_coords
        hessian(i_hessian,i_hessian) = 1.d0
     enddo
     initial_hessian = .true.
     line_search_on  = .true.

  ! check if the last step was part 2 of a line search, if so update coordinates directly!
  else if (line_search_on) then

     write(use_unit,'(2X,A)') 'Using parabolic line search to determine next position:'
     ! this is part 1 of the derivative along the line step
     force_times_search_dir_old = force_times_search_dir

     ! calculate   forces.search_dir 
     force_times_search_dir = ddot(n_coords,total_forces,1,search_direction,1)

     write(use_unit,*) ' | force . search_direction (old) ', force_times_search_dir_old
     write(use_unit,*) ' | force . search_direction (new) ', force_times_search_dir
     write(use_unit,'(A,F10.6)') '  | line_step (old)                 ', line_step
     ! calculate new line step
     if ((max_line_step_increase.ge.1d0).and. &
         (force_times_search_dir/(force_times_search_dir_old-force_times_search_dir).gt.max_line_step_increase)) then
        write(use_unit,*) ' | The predicted new line step is too large: restricting to control file requirements'
        line_step = line_step*max_line_step_increase
     else
        line_step = line_step*force_times_search_dir/(force_times_search_dir_old-force_times_search_dir)     
     end if
     line_search_on = .false.    ! next step: no line search!
     
     write(use_unit,'(A,F10.6)') '  | line_step (new)                 ', line_step

  else if ((.not.object_function_exists).or. ( ((object_function - stored_object_function).lt.0.d0) .or. & 
       ( ((object_function - stored_object_function).lt.object_function_tolerance) .and. &
       initial_hessian .and. (line_step.lt.min_line_step)) )) then

!  else ! if (initial_hessian .and. (line_step.lt.min_line_step) ) then
     ! Step accepted if
     !   (a) no object function exists, i.e. for NEB
     !   (b) object function decreases
     !   (c) prevents noise at the end of a run, where the forces may point slightly upwards
     !       in energy simply due to numerical inaccuracies at the 10^-5 eV level ...
     ! update Hessian matrix using previous results
     ! compute force difference between present and previous geometry
     force_difference(:) = stored_forces(:) - total_forces(:)
     ! scalar product: force_difference * search_direction
     forcediff_times_searchdir = line_step * ddot(n_coords,force_difference,1,search_direction,1)
     ! Hessian times search_direction: Hs
     call dsymv('u',n_coords,line_step,hessian,n_coords,search_direction,1,0.d0,Hs,1)
     ! scalar product: search_direction * hessian * search_direction
     sHs = line_step*ddot(n_coords,search_direction,1,Hs,1)
     ! Now update the hessian ...
     call dsyr('u',n_coords,1.d0/forcediff_times_searchdir,force_difference,1,hessian,n_coords ) 
     call dsyr('u',n_coords,-1.d0/sHs,Hs,1,hessian,n_coords) 

     ! Store energy, coordinates, forces of present iteration for next iteration
     stored_forces          = total_forces
     stored_coords          = coords
     stored_object_function = object_function
     ! restore line_step
     line_step = 1.0
     initial_hessian = .false.
     ! Use current Hessian matrix and current forces to predict search direction s, 
     ! by inverting hessian * search_direction = - total_forces
     ! preload current forces into search_direction (after dsysv, this array will contain the actual search direction)
     ! note that search_direction is a 1D array for technical reasons but contains the search direction in
     ! the same order as a 2D array search_direction(i_coord,i_atom)
     search_direction(:) = total_forces(:)
     ! actually solve system of linear equations (lapack):
     ! first determine optimal workspace size with call to dsysv
     allocate(work(1))
     ipiv = 0
     lwork = -1
     call dsysv('u',n_coords,1,hessian_work,n_coords,ipiv,search_direction,n_coords,work,lwork,info) 
     lwork = work(1)
     deallocate(work)
     allocate(work(lwork))
     hessian_work = hessian
     call dsysv('u',n_coords,1,hessian_work,n_coords,ipiv,search_direction,n_coords,work,lwork,info) 
     ! check whether dgesv gave a reasonable result:
     if (info.lt.0) then
        write(use_unit,'(1X,A,I3,A)')"* Search direction: Argument number ",-info," to dsysv had an illegal value."
        stop
     else if (info.gt.0) then
        write(use_unit,'(1X,A,I3,A)')"* Search direction: Hessian factorisation ",info," is exactly zero."
        write(use_unit,'(1X,A)')"Your input Hessian may be singular - cannot proceed."
        stop     
     end if
     ! check whether the search direction makes any sense compared to the current forces
     ! scalar product: force_difference * search_direction / |search_direction|
     force_times_search_dir =  ddot(n_coords,total_forces,1,search_direction,1)

     ! if the forces point away from the current search direction, the Hessian can't be any good
     ! in that case, reinitialize the Hessian
     if (force_times_search_dir.lt.0.d0) then
        write(use_unit,'(2X,A,E14.6)')"BFGS: Forces * search direction is negative: ", force_times_search_dir
        write(use_unit,'(2X,A)')"Search direction should not point away from forces - reinitializing Hessian."
        line_step = 1.0
        search_direction(:) = total_forces(:)
        ! reinitialize hessian to be unity ...
        ! This means the first step will be a simple steepest descent.
        hessian = 0.d0
        do i_hessian = 1, n_coords 
           hessian(i_hessian,i_hessian) = 1.d0
        enddo
        initial_hessian = .true.
     end if
     line_search_on = .true.

!
!    FIXME: in the original AIMS code, this was the part in which we decided whether a step was accepted or 
!           not, based on (1) minimal line step (to avoid numerical noise)
!                         (2) energy increase (quit if too large)
!                         (3) something I am forgetting, check relaxation.f90 in aims source directory.
!           While this self-consistency check is, in principle, a very great idea, it was causing trouble in
!           the current NEB formalism. What criteria (if any?) are applicable here??? 
!
  else      
     ! Step was not accepted - revert to previous coordinates, 
     ! use same search direction as before but decrease line_step
     ! stored_object_function remains the same as the last accepted value should be used as reference
     coords = stored_coords
     total_forces = stored_forces
     line_step = line_step * line_step_reduce
     if ((line_step.ge.min_line_step).or.initial_hessian) then
        ! initial hessian: line step will be decreased until we find a geometry that 
        ! remains below the built-in noise threshold
        write(use_unit,'(2X,A)')"BFGS line minimization: Last trial step was too large."
        write(use_unit,'(2X,A,E14.6)')"Reverting to previous geometry - new trial step: ", line_step
     else
        write(use_unit,'(2X,A)')"BFGS: Extrapolation using current Hessian along current search direction"
        write(use_unit,'(2X,A)')"is not nearly quadratic - incorrect Hessian for current point?"
        write(use_unit,'(2X,A)')"Resetting search direction and Hessian to steepest descent."
        line_step = 1.0
        search_direction(:) = total_forces(:)
        ! reinitialize hessian to be unity ...
        ! This means the first step will be a simple steepest descent.
        hessian = 0.d0
        do i_hessian = 1, n_coords
           hessian(i_hessian,i_hessian) = 1.d0
        enddo
        initial_hessian = .true.
     end if
  end if

  ! now that we have a current search direction, make sure that the displacement remains below
  ! a threshold
  ! Check for overly large movements and reduce line_step if needed
  ! this can only happen after the search direction was changed, i.e. when line_step.eq.1)
  if ((line_step.eq.1.0d0).and.(line_search_on)) then
     ! find out max atomic displacement ... 
     max_displacement = abs(line_step*search_direction(1))
     do i_coord = 1, n_coords
        if ( abs(line_step*search_direction(i_coord)).gt.max_displacement ) then
           max_displacement = abs(line_step*search_direction(i_coord))
        end if
     end do
     ! check whether that's too big?
     if (max_displacement.gt.max_atomic_move) then
        line_step = max_atomic_move/max_displacement
        write(use_unit,'(2X,A,E14.6,A)') &
             "BFGS: Suggested max. displacement in a single coordinate = ", &
             max_displacement, " A."
        write(use_unit,'(2X,A,E14.6,A,E14.6,A)') &
             "Reducing line step to ", line_step, ", to avoid exceeding ", &
             max_atomic_move, " A."
     end if
  end if

  ! just in case this has not been calculated previously, we need this information for the line search!
  force_times_search_dir =  ddot(n_coords,total_forces    ,1,search_direction,1)

  ! We now have a new search direction; must update the actual coordinates
  coords(:) = coords(:) + line_step * search_direction(:)
end subroutine BFGS_n_coords
