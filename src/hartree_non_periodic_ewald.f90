!****h* FHI-aims/hartree_non_periodic_ewald
!
!  NAME
!    hartree_non_periodic_ewald
!
!  SYNOPSIS

module hartree_non_periodic_ewald

  !  PURPOSE
  !
  !    If the switch "use_hartree_non_periodic_ewald" is set to
  !    "true", then the Hartree term in the non-periodic case is
  !    decomposed into two parts according to Ewald's method. The
  !    present module contains the necessary data structures and
  !    subroutines for this task.
  !
  !    NOTE: THIS MODULE IS EXPERIMENTAL AT PRESENT!
  !
  !  AUTHOR
  !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
  !    Please note that any use of the FHI-aims-Software is subject to
  !    the terms and conditions of the respective license agreement.
  !
  !  SOURCE

  use constants ! XX
  implicit none

  !  Global variables of the module

  save

  real*8, parameter :: boinpm =  100 * bohr    ! XX Bohr radius in pm
  real*8, parameter :: HinmV  = 1000 * hartree ! XX Hartree energy in mV

  ! The following parameters determine the Cartesian grid. They are
  ! shared by the subroutines
  !  - initialize_hartree_non_periodic_ewald
  !  - calculate_extended_hartree_non_periodic_ewald
  ! while for
  !  - determine_bounding_box
  ! better encapsulation is achieved by using local variables and
  ! return arguments.

  real*8,                 private :: gridwidth
  real*8, dimension(3),   private :: origin
  real*8, dimension(3,3), private :: axes
  integer, dimension(3),  private :: maxind


  ! Ewald's compensating potential on the Cartesian grid

  real*8, dimension(:,:,:),   allocatable, private :: cartpot
  real*8, dimension(:,:,:,:), allocatable, private :: deriv2


  ! Interpolation methods

  integer, parameter :: const   = 1   ! constant interpolation
  integer, parameter :: trilin  = 2   ! trilinear interpolation
  integer, parameter :: numherm = 3   ! tricubic Hermite interpolation based
                                      !   on numerical derivatives
  integer, parameter :: ispline = 4   ! spline interpolation

  integer, parameter :: intpmeth = ispline






  contains






  !****s* hartree_non_periodic_ewald/initialize_hartree_non_periodic_ewald
  !
  !  NAME
  !    initialize_hartree_non_periodic_ewald
  !
  !  SYNOPSIS

  subroutine initialize_hartree_non_periodic_ewald

    !  PURPOSE
    !    Perform the initialization for the Ewald decomposition
    !  USES

    use constants
    use mpi_tasks
    use runtime_choices

    !  AUTHOR
    !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE

    implicit none

    integer :: info

    ! Start of initialize_hartree_non_periodic_ewald

!    if (myid == 0) then
!      print *
!      print *, 'initialize_hartree_non_periodic_ewald: gridwidth in pm:', gridwidth_in_pm
!      print *
!    end if

    gridwidth = gridwidth_in_pm / boinpm


    if (intpmeth == numherm) then   ! extraplanes = 1 for Hermite with numerical derivatives
      call determine_bounding_box( gridwidth, 1, origin, axes, maxind )
    else
      call determine_bounding_box( gridwidth, 0, origin, axes, maxind )
    end if


    if ( allocated(cartpot) )   deallocate(cartpot)
    allocate(   cartpot(     0:maxind(1), 0:maxind(2), 0:maxind(3) ),  stat=info   )
    call check_allocation( info, 'cartpot' )


    if (intpmeth == ispline) then
      if ( allocated(deriv2) )    deallocate(deriv2)
      allocate(   deriv2( 1:3, 0:maxind(1), 0:maxind(2), 0:maxind(3) ),  stat=info   )
      call check_allocation( info, 'deriv2' )
    end if


  end subroutine initialize_hartree_non_periodic_ewald
  !******







  !****s* hartree_non_periodic_ewald/determine_bounding_box
  !
  !  NAME
  !    determine_bounding_box
  !
  !  SYNOPSIS

  subroutine determine_bounding_box( gridwidth, extraplanes,  &
                                     origin, axes, maxind     )

    !  PURPOSE
    !
    !    Determine an arbitrarily oriented box with minimal volume
    !    that contains the molecule. At present, the bounding box is
    !    determined by considering the locations of the nuclei and
    !    adding afterwards some space to include the radial
    !    integration grids around the nuclei. In order to determine
    !    the bounding box for the nuclei, an approximate method is
    !    used based on principal component analysis, see:
    !
    !    - D. Dimitrov, C. Knauer, K. Kriegel, G. Rote:
    !      "Bounds on the quality of the PCA bounding boxes"
    !      Computational Geometry, vol. 42, 2009, p. 772-789
    !
    !    - M. Lahanas, T. Kemmerer, N. Milickovic, K. Karouzakis,
    !      D. Baltas, N. Zamboglou:
    !      "Optimized bounding boxes for three-dimensional treatment
    !       planning in brachytherapy"
    !      Medical Physics, vol. 27, 2000, p. 2333-2342
    !
    !    A method for obtaining the _true_ minimal-volume box is
    !    available, but this procedure is very complex, see:
    !
    !    - J. O'Rourke:
    !      "Finding minimal enclosing boxes"
    !      International Journal of Computer and Information Sciences,
    !      vol. 14, 1985, p. 183-199
    !
    !  USES

    use mpi_tasks
    use dimensions
    use localorb_io
    use grids
    use basis
    use geometry
    use numerical_utilities

    !  ARGUMENTS

    implicit none

    real*8,                 intent(in)  :: gridwidth
    integer,                intent(in)  :: extraplanes  ! XX documentation

    real*8, dimension(3),   intent(out) :: origin
    real*8, dimension(3,3), intent(out) :: axes
    integer, dimension(3),  intent(out) :: maxind

    ! Note that module-global variables with the same name exist, but
    ! it is sufficient for the present subroutine to use local
    ! variables and return arguments.

    !  INPUTS
    !    - gridwidth: mesh width of the Cartesian grid
    !
    !  OUTPUT
    !    - origin: point (0,0,0) of the bounding box
    !    - axes:   the three orthogonal axes aligned with the edges of the bounding box
    !    - maxind: highest index in the respective direction
    !
    !  AUTHOR
    !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE

    integer  ::  i, j, k    ! indices 1, 2, 3 for the original coordinate directions x, y, z
    integer  ::  at         ! index over atoms

    integer :: i_coord

    real*8, dimension(3)          :: mean
    real*8, dimension(3, n_atoms) :: relative_coords

    real*8, dimension(3,3)        :: cov
    real*8, dimension(3)          :: eigenval
    integer                       :: ev          ! index 1, 2, 3 for the eigenvectors  XX remove later
    integer                       :: dir         ! direction 1, 2, 3 of the new axes

    real*8, dimension(n_atoms)    :: proj
    real*8                        :: proj_min, proj_max, length
    real*8, dimension(3)          :: mincoor

    integer                       :: sp          ! index over species     XX
    real*8                        :: r_radial_outer_max  ! XX
    real*8                        :: outer_radius_max

   character*150 :: info_str

    ! Start of determine_bounding_box

!    if (myid == 0) then
!
      
     call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

     write(info_str,'(A)') &
       & '========== Bounding box determination for Ewald-like Hartree potential interpolation grid ============'
     call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

     call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

!
!      print *, 'positions of the nuclei:'
!      do at = 1, n_atoms
!        print *, coords(:,at)
!      end do
!      print *
!
!    end if


    ! 1. Determine the mean location and relative coordinates of the
    !    set of nuclei.

    call get_relative_coords( coords, mean, relative_coords )

 !   if (myid == 0) then
 !     print *
 !     print *, 'mean:'
 !     print *, mean
 !     print *
 !   end if


    ! 2. Calculate the covariance matrix between the scatter of the
    !    nuclei positions in the coordinate directions i and j. Only
    !    the upper triangle of cov is computed.

    do i = 1, 3
      do j = i, 3
        cov(i,j)  =  dot_product(  relative_coords(i,:), relative_coords(j,:)  )  /  n_atoms
      end do
    end do

!    if (myid == 0) then
!      print *
!      print *, 'cov:'
!      do i = 1, 3
!        print *, cov(i,:)
!      end do
!      print *
!    end if


    ! 3. Determine the eigenvectors of cov, which represent the
    !    principal components of this matrix. The edges of the
    !    bounding box will be aligned with the eigenvectors.

    call diagonalize_rmatrix( 3, cov, eigenval, .true. )

    !    The eigenvectors are now the columns of cov, that is,
    !    cov(:,ev) is eigenvector number ev.  The eigenvectors are
    !    assumed to be unit vectors.  XX normalize explicitly

!    if (myid == 0) then
!
!      print *
!      print *, 'eigenvalues:'
!      print *, eigenval
!      print *
!
!      print *
!      print *, 'eigenvectors:'
!      do ev = 1, 3
!        print *, cov(ev,:)
!      end do
!      print *
!
!    end if

    ! 4. Copy the eigenvectors in `cov' to the array `axes' that
    !    contains the axes. The axis vector in direction dir = 1, 2, 3
    !    is then given by axes(:,dir).

    axes = cov


    ! 5. Now, the size of the box has to be calculated. This is done
    !    for each direction `dir' (that is, for each edge of the box.)

    !  5a) First, determine the maximum of outer_radius


 !   if (myid == 0)  print *, 'outer_radius:', outer_radius
    outer_radius_max = maxval(outer_radius)
 !   if (myid == 0)  print *, 'outer_radius_max:', outer_radius_max

!    if (myid == 0)  print *
    r_radial_outer_max = 0
    do sp = 1, n_species
!      if (myid == 0) print *, 'r_radial( n_radial(sp), sp ):', r_radial( n_radial(sp), sp )
      r_radial_outer_max  =  max(  r_radial_outer_max,  r_radial(n_radial(sp),sp)  )
    end do
!    if (myid == 0) print *, 'r_radial_outer_max: ', r_radial_outer_max

    outer_radius_max = r_radial_outer_max  ! XX


    do dir = 1, 3  ! loop over directions

      ! 5b) Calculate the projection `proj' of all nuclei positions
      !     onto the direction `dir'.

      do at = 1, n_atoms
        proj(at) = dot_product( coords(:,at), axes(:,dir) )
      end do

      ! 5c) Determine the minimum and maximum of the projections since
      !     these have to be included at least in the bounding box.

      proj_min  =  minval(proj)
      proj_max  =  maxval(proj)
      ! print *, 'proj_min:', proj_min
      ! print *, 'proj_max:', proj_max

      ! 5d) The (minimum) length of the bounding box in the present
      !     direction is given by the follwing formula because on each
      !     side `outer_radius_max' has to be added to include the
      !     radial grid.  XX

      length  =  proj_max - proj_min + 2*outer_radius_max
      ! print *, 'length: ', length

      ! 5e) The final length of the bounding box is somewhat larger
      !     than `length' since it has to be an integer multiple of
      !     `gridwidth'. The maximum index of the grid in the
      !     considered direction is thus found by rounding up in the
      !     following manner:

      maxind(dir) = ceiling( length/gridwidth )

      ! 5f) Now, the length of the box in direction `dir' is
      !     maxind(dir)*gridwidth (>= length). Half of the extra
      !     length is added on each side of the box. The minimum
      !     coordinate (location of the face of the box) in direction
      !     `dir' is thus:

      mincoor(dir)  =  proj_min - outer_radius_max - 0.5*( maxind(dir)*gridwidth - length )

      ! 5g) Finally, extras planes (if requested) have do be added to
      !     each side of the box:

      mincoor(dir) =  mincoor(dir) -   extraplanes*gridwidth
      maxind(dir)  =  maxind(dir)  + 2*extraplanes

    end do ! dir


    ! 6. The 'origin' of the box, that is, the corner point with the
    !    smallest coordinates in the Cartesian `axes'-system,
    !    expressed in the basic coordiate system (x,y,z) is:

    origin = mincoor(1) * axes(:,1) +  &
             mincoor(2) * axes(:,2) +  &
             mincoor(3) * axes(:,3)

!  In principle, bounding box is finished. Write results.

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

      write(info_str,'(A,2X,F20.8,2X,F20.8,2X,F20.8)') &
        & 'Grid origin:                 ', origin(1), origin(2), origin(3)
      call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

      write(info_str,'(A,2X,I8,2X,I8,2X,I8)') &
        & 'Grid size in each direction: ', maxind(1), maxind(2), maxind(3) 
      call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)


      write(info_str,'(A,2X,I8)') &
        & 'Number of cubes in grid    : ', maxind(1) * maxind(2) * maxind(3) 
      call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

      write(info_str,'(A,2X,F20.8,2X,F20.8,2X,F20.8)') &
        & 'Total grid extension (pm)  : ', & 
        & maxind(1) * gridwidth * boinpm, &
        & maxind(2) * gridwidth * boinpm, &
        & maxind(3) * gridwidth * boinpm
      call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

      write(info_str,'(A)') &
        & 'Corner points in pm :'
      call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

      do i = 0, maxind(1),  maxind(1)
        do j = 0, maxind(2),  maxind(2)
          do k = 0, maxind(3),  maxind(3)

             write(info_str,'(F20.8,2X,F20.8,2X,F20.8)') &
             ( boinpm * (   &
                       origin(i_coord) + i * gridwidth * axes(i_coord,1)  &
                              + j * gridwidth * axes(i_coord,2)  &
                              + k * gridwidth * axes(i_coord,3)     ), &
             i_coord = 1,3,1 ) 
             call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

          end do
        end do
      end do

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)
      write(info_str,'(A)') &
        & '=========================================================='
      call localorb_info(info_str, use_unit,'(2X,A)',OL_norm)

      call localorb_info(" ", use_unit,'(2X,A)',OL_norm)

  end subroutine determine_bounding_box
  !******







  !****s* hartree_non_periodic_ewald/calculate_extended_hartree_non_periodic_ewald
  !
  !  NAME
  !    calculate_extended_hartree_non_periodic_ewald
  !
  !  SYNOPSIS

  subroutine calculate_extended_hartree_non_periodic_ewald(                           &
               l_hartree_max_far_distance, multipole_radius_sq, adap_outer_radius_sq  )

    !  PURPOSE
    !    Sum up the extended ('long-range') part of Ewald's decomposition of the
    !    Hartree potential on the equidistant Cartesian grid.

    !  USES

    !  use constants XX
    use mpi_tasks
    use geometry
    use species_data
    use localorb_io
    use hartree_potential_real_p0
    use timing  ! XX

    !  ARGUMENTS

    implicit none

    integer, dimension(n_atoms), intent(in) :: l_hartree_max_far_distance
    real*8,  dimension(n_atoms), intent(in) :: multipole_radius_sq, adap_outer_radius_sq

    !  AUTHOR
    !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE


    ! Local variables

    integer                               :: i, j, k,  at
    real*8, dimension(3)                  :: cartpoint
    real*8, dimension(:,:,:), allocatable :: cartpot_mpi
    integer                               :: info, mpierr
    real*8                                :: cpu_time_gridsum, clock_time_gridsum


    ! Start of calculate_extended_hartree_non_periodic_ewald

    call localorb_info & 
    & ("  | Summing up the extended ('long-range') part of Ewald's decomposition of the Hartree potential.", &
    &   use_unit,'(A)',OL_norm)

    call get_timestamps( cpu_time_gridsum, clock_time_gridsum )


    ! 1. Compute the extended Ewald part for my k-slice of the Cartesian grid.

    cartpot = 0

    do k = 0, maxind(3)

      if  (  k * n_tasks  >=   myid    * (maxind(3)+1)  .and.  &
             k * n_tasks  <   (myid+1) * (maxind(3)+1)           )  then

        do i = 0, maxind(1)
          do j = 0, maxind(2)   ! second longest loop as innermost loop

            cartpoint = origin + gridwidth * (   i * axes(:,1)  &
                                               + j * axes(:,2)  &
                                               + k * axes(:,3)     )

            do at = 1, n_atoms

              call far_distance_hartree_Fp_periodic_single_atom &
                   ( at, at, sqrt(  sum(  (coords(:,at) - cartpoint)**2 )   ),   &
                     l_hartree_max_far_distance, .true., .false.,                &
                     multipole_radius_sq(at), sqrt(adap_outer_radius_sq(at)),    &
                     non_peri_extd = .true.                                         )  ! XX forces

              ! XX inside (=.true.) is irrelevant

              call far_distance_real_hartree_potential_single_atom &
                   ( at, at, cartpot(i,j,k), l_hartree_max_far_distance, cartpoint )

              ! call far_distance_real_hartree_potential_single_atom_p2 &
              !      ( at, cartpot(i,j,k), l_hartree_max_far_distance(at), cartpoint )
              !    ! ( at, cartpot(i,j,k), 0                             , cartpoint )

            end do ! at

          end do ! j
        end do ! i

      end if  ! ( k * n_tasks ... )

    end do ! k



    ! 2. Distribute the k-slices among the processors.

    if (USE_MPI) then
       
       allocate(cartpot_mpi(0:maxind(1), 0:maxind(2), 0:maxind(3)), &
             stat=info   )
       call check_allocation( info, 'cartpot_mpi' )

       cartpot_mpi = 0

       call mpi_allreduce( cartpot, cartpot_mpi, size(cartpot), &
               mpi_double_precision, mpi_sum, mpi_comm_global, mpierr )

       cartpot = cartpot_mpi

       deallocate(cartpot_mpi)

    endif


    ! 3. Calculate the spline coefficients.

    if (intpmeth == ispline)  call calculate_2nd_deriv_from_splines



    call get_times( cpu_time_gridsum, clock_time_gridsum )
    call output_times( '2x', 'Cartesian grid summation', cpu_time_gridsum, clock_time_gridsum, OL_norm )

    ! call hartree_non_periodic_ewald_tests(l_hartree_max_far_distance)


  end subroutine calculate_extended_hartree_non_periodic_ewald
  !******







  !****s* hartree_non_periodic_ewald/interpolate_extended_hartree_non_periodic_ewald
  !
  !  NAME
  !    interpolate_extended_hartree_non_periodic_ewald
  !
  !  SYNOPSIS

  subroutine interpolate_extended_hartree_non_periodic_ewald( coor, potential, gradient )

    !  PURPOSE
    !    Calculate the potential and (if requested) also the gradient of the
    !    extended part of Ewald's decomposition. Potential and gradient are
    !    determined by interpolation on the even spaced grid. The target point
    !    is "coor".
    !      Note that the interpolated potential is _added_ to variable
    !    "potential", whereas the gradient is _returned_ (not added) in variable
    !    "gradient".

    !  ARGUMENTS

    implicit none

    real*8, dimension(3), intent(in)             :: coor
    real*8,               intent(inout)          :: potential
    real*8, dimension(3), intent(out),  optional :: gradient

    !  AUTHOR
    !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE


    ! Local variables

    integer               :: i, j, dir, graddir            ,k !XX remove k
    integer, dimension(3) :: ind
    real*8,  dimension(3) :: d, gradcart
    real*8                :: coordinate,  u1, u2,  v1, v2,  w1, w2,  interp_pot


    ! Start of interpolate_extended_hartree_non_periodic_ewald


    do dir = 1, 3

      coordinate  =  dot_product(  coor - origin,  axes(:,dir)  )
      ind(dir)    =  floor( coordinate / gridwidth )

      if (   intpmeth == numherm   .and.                          &
             (  ind(dir) < 1  .or.  ind(dir) > maxind(dir)-2  )   ) then
        print *, 'stop: to close to the faces of the box. dir =', dir, ' ind(dir) =', ind(dir)
        stop
      end if

      d(dir)      =  ( coordinate - ind(dir)*gridwidth ) / gridwidth

    end do



    select case (intpmeth)


      case (const)

        interp_pot  =  cartpot( ind(1),  ind(2),  ind(3) )


      case (trilin)

        u1  =  cartpot(ind(1),  ind(2),  ind(3)) * (1-d(3))  +  cartpot(ind(1),  ind(2),  ind(3)+1) * d(3)
        u2  =  cartpot(ind(1),  ind(2)+1,ind(3)) * (1-d(3))  +  cartpot(ind(1),  ind(2)+1,ind(3)+1) * d(3)
        v1  =  cartpot(ind(1)+1,ind(2),  ind(3)) * (1-d(3))  +  cartpot(ind(1)+1,ind(2),  ind(3)+1) * d(3)
        v2  =  cartpot(ind(1)+1,ind(2)+1,ind(3)) * (1-d(3))  +  cartpot(ind(1)+1,ind(2)+1,ind(3)+1) * d(3)

        w1  =  u1 * (1-d(2)) + u2 * d(2)
        w2  =  v1 * (1-d(2)) + v2 * d(2)

        interp_pot  =  w1 * (1-d(1)) +  w2 * d(1)


      case (numherm)

        ! 1. Calculate the potential.

        interp_pot = numhermite_interp(0)    ! graddir = 0 means no derivative


        ! 2. Calculate the gradient if requested.

        if ( present(gradient) ) then

          do graddir = 1, 3
            gradcart(graddir) = numhermite_interp(graddir)
          end do

        end if


      case (ispline)

        if ( .not. present(gradient) ) then

          call splint( ind, d, interp_pot )

        else

          ! XX begin tests

          ! do i = 0, maxind(1)
          !   do j = 0, maxind(2)
          !     do k = 0, maxind(3)

          !       cartpot(i,j,k) =  i + k**3

          !     end do
          !   end do
          ! end do

          ! call calculate_2nd_deriv_from_splines

          ! ind = (/ 10, 10, 10  /)
          ! d   = (/ 0.,  0., 0. /)

          ! print *, 'deriv2(:,ind(1),ind(2),ind(3)), deriv2(:,ind(1)+1,ind(2),ind(3))', &
          !           deriv2(:,ind(1),ind(2),ind(3)), deriv2(:,ind(1)+1,ind(2),ind(3))

          ! XX end tests

          call splint( ind, d, interp_pot, gradcart )

          ! XX begin tests
          ! print *, 'ind, d, interp_pot, gradcart:', ind, d, interp_pot, gradcart
          ! stop
          ! XX end tests

        end if


    end select  ! intpmeth



    potential = potential + interp_pot



    if ( present(gradient) ) then

      ! Transform the gradient from the Cartesian coordinate system (with unit
      ! mesh width) to the primary (x,y,z) coordinate system of aims.

      gradient(:)   =   (  gradcart(1) * axes(:,1) +  &
                           gradcart(2) * axes(:,2) +  &
                           gradcart(3) * axes(:,3)     )   /   gridwidth
    end if



    contains



    function numhermite_interp(graddir)

      integer, intent(in)          :: graddir
      real*8                       :: numhermite_interp

      real*8, dimension(4)         :: ipvx, ipvy, ipvz
      real*8, dimension(-1:2,-1:2) :: t
      real*8, dimension(-1:2)      :: u

      ! XX cite wiki; reorder x, y, z


      ! 1. Interpolate in the Cartesian z-direction

      ipvz = ipvec_choose( 3, graddir )
      do i = -1, 2
        do j = -1, 2
          t(i,j) = dot_product( ipvz,                                             &
                                cartpot( ind(1)+i, ind(2)+j, ind(3)-1:ind(3)+2 )   )
        end do
      end do


      ! 2. Interpolate in the Cartesian y-direction

      ipvy = ipvec_choose( 2, graddir )
      do i = -1, 2
        u(i)  =  dot_product(  ipvy,  t( i, -1:2 )  )
      end do


      ! 3. Interpolate in the Cartesian x-direction

      ipvx = ipvec_choose( 1, graddir )

      numhermite_interp = dot_product( ipvx, u )


    end function numhermite_interp



    function ipvec_choose( dir, graddir )

      integer, intent(in)   :: dir, graddir
      real*8,  dimension(4) :: ipvec_choose

      if ( dir == graddir ) then
        ipvec_choose  =  ipvecgrad(d(dir))
      else
        ipvec_choose  =  ipvec    (d(dir))
      end if

    end function ipvec_choose


    ! "ipvec" means "interpolation vector"

    function ipvec(d)
      real*8, intent(in)    :: d
      real*8, dimension(4)  :: ipvec
      ipvec(1)  =  0.5  *  d  *  (  (2-d)*d - 1  )
      ipvec(2)  =  0.5  *  (  d**2 * (3*d-5) + 2  )
      ipvec(3)  =  0.5  *  d  *  (  (4-3*d)*d + 1  )
      ipvec(4)  =  0.5  *  (d-1) * d**2
    end function ipvec


    function ipvecgrad(d)
      real*8, intent(in)    :: d
      real*8, dimension(4)  :: ipvecgrad
      ipvecgrad(1)   =   0.5  *  ( -3*d**2 +  4*d - 1  )
      ipvecgrad(2)   =   0.5  *  (  9*d**2 - 10*d      )
      ipvecgrad(3)   =   0.5  *  ( -9*d**2 +  8*d + 1  )
      ipvecgrad(4)   =   0.5  *  (  3*d**2 -  2*d      )
    end function ipvecgrad



  end subroutine interpolate_extended_hartree_non_periodic_ewald
  !******







  !****s* hartree_non_periodic_ewald/add_true_hartree_non_periodic_ewald
  !
  !  NAME
  !    add_true_hartree_non_periodic_ewald
  !
  !  SYNOPSIS

  subroutine add_true_hartree_non_periodic_ewald( coor, potential, l_hartree_max_far_distance, &
                                                  mindist )

    !  USES

    ! XX which necessary?
    use geometry
    use species_data
    use hartree_potential_real_p0
    use grids, only: r_grid_min

    implicit none

    integer, dimension(n_atoms) :: l_hartree_max_far_distance

    real*8, dimension(3), intent(in)    :: coor
    real*8,               intent(inout) :: potential
    logical                             :: mindist

    !  AUTHOR
    !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE

    integer  :: at
    real*8   :: dist_sq, true_pot

    ! Start of add_true_hartree_non_periodic_ewald

    true_pot = 0

    do at = 1, n_atoms

      dist_sq  =  sum(  (coords(:,at) - coor)**2  )

      if ( mindist .and. dist_sq < 1e-20 ) then
        dist_sq = r_grid_min(species(at))**2 + 1e-15
      end if

      ! call far_distance_hartree_Fp_periodic_single_atom &
      !      ( at, at, sqrt(dist_sq),                                         &
      !        l_hartree_max_far_distance, inside=.true., forces_on=.false.,  &
      !        non_peri_extd = .true.                                        )

      print *, 'stop: not adapted'  ! concerning the radii
      stop

      call far_distance_real_hartree_potential_single_atom &
           ( at, at, true_pot, l_hartree_max_far_distance, coor )

    end do ! at

    potential = potential + true_pot


  end subroutine add_true_hartree_non_periodic_ewald
  !******







  !****s* hartree_non_periodic_ewald/hartree_non_periodic_ewald_tests
  !
  !  NAME
  !    hartree_non_periodic_ewald_tests
  !
  !  SYNOPSIS

  subroutine hartree_non_periodic_ewald_tests(l_hartree_max_far_distance)

    !  USES

    !  use constants XX
    use geometry
    use species_data
    use localorb_io
    use hartree_potential_real_p0
    use grids, only: batches

    implicit none

    integer, dimension(n_atoms) :: l_hartree_max_far_distance

    !  AUTHOR
    !    Werner J"urgens, Fritz-Haber Institute of the Max-Planck-Society, Germany
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE

    integer                :: i, j, k,  at, b, p,  dir, sig
    integer, dimension(3)  :: ind
    real*8, dimension(3)   :: integr_point_coord, coor, d
    real*8                 :: coordinate,  interp_pot, true_pot

    ! Start of hartree_non_periodic_ewald_tests


    ! Write the Cartesian potential to `cartpot.cube'.

    if ( maxind(1) * maxind(2) * maxind(3) <= 200000 ) then
      open( 123, file = 'cartpot.cube' )
      write(123,*) 'potential'
      write(123,*) 'relative to:', maxval(abs(cartpot))
      write(123,'( "1  ", 3(1x,g12.6) )')  origin
      write(123,'( i0,3x, 3(1x,g12.6) )')  maxind(1)+1,  gridwidth * axes(:,1)
      write(123,'( i0,3x, 3(1x,g12.6) )')  maxind(2)+1,  gridwidth * axes(:,2)
      write(123,'( i0,3x, 3(1x,g12.6) )')  maxind(3)+1,  gridwidth * axes(:,3)
      write(123,*) '1  0  0 0 0'
      do i = 0, maxind(1)
        do j = 0, maxind(2)
          do k = 0, maxind(3)
            write(123,*)  cartpot(i,j,k) / maxval(abs(cartpot))
          end do
        end do
      end do
      close(123)
    end if


    ! Determine statistics for interpolation errors on the integration grid points.

    open( 123, file = 'interpdiff.dat' )
    write( 123, '("# no. of cubes: ", i0)' )  maxind(1) * maxind(2) * maxind(3)

    do b = 1, n_grid_batches
      do p = 1, batches(b) % size
        ! print *, 'b:', b, '  p:', p, ' coords:', batches(b) % points(p) % coords
        integr_point_coord = batches(b) % points(p) % coords

        at = batches(b) % points(p) % index_atom

        ! if (  sum(  ( integr_point_coord -           &
        !               coords( :, at )                &
        !             )**2                             &
        !           )                                  &
        !        <  atom_radius_sq(  species(at )  )   &
        !     ) then

          do dir = 1, 3
            coordinate  =  dot_product(  integr_point_coord - origin,  axes(:,dir)  )
            ind(dir)    =  floor( coordinate / gridwidth )
            ! d(dir)      =  ( coordinate - ind(dir)*gridwidth ) / gridwidth
          end do


          interp_pot = 0
          call interpolate_extended_hartree_non_periodic_ewald( integr_point_coord, interp_pot )

          true_pot = 0
          call add_true_hartree_non_periodic_ewald( integr_point_coord, true_pot, &
                                                    l_hartree_max_far_distance, mindist=.false. )

          sig = 0
          do i = 0, 1
            do j = 0, 1
              do k = 0, 1
                sig  =  sig  +  nint(   sign(  real( 1, kind(cartpot) ),                &
                                               cartpot( ind(1)+i, ind(2)+j, ind(3)+k )  &
                                            )                                           &
                                    )
              end do
            end do
          end do


          ! if ( b == 1 .and. p < 10 ) print *, 'b:', b, '  p:', p, ' coords:', integr_point_coord, 'ind:', ind, &
          !                                                         ' true_pot:', true_pot, ' interp_pot:', interp_pot
          write(123,'( i4,1x, i4,5x, 1p, 3(1x,g13.6), 5x, 3(1x,g13.6), 3x,i0 )')  &
                         b, p,  integr_point_coord,  &
                         true_pot*HinmV, interp_pot*HinmV,  (interp_pot - true_pot)*HinmV,  abs(sig)/8

          ! if ( abs(sig) < 8 ) then
          !   print *, 'true_pot:', true_pot,  'interp_pot:', interp_pot
          !   print *,  cartpot(ind(1),  ind(2),  ind(3)),  cartpot(ind(1),  ind(2),  ind(3)+1)
          !   print *,  cartpot(ind(1),  ind(2)+1,ind(3)),  cartpot(ind(1),  ind(2)+1,ind(3)+1)
          !   print *,  cartpot(ind(1)+1,ind(2),  ind(3)),  cartpot(ind(1)+1,ind(2),  ind(3)+1)
          !   print *,  cartpot(ind(1)+1,ind(2)+1,ind(3)),  cartpot(ind(1)+1,ind(2)+1,ind(3)+1)
          ! end if

        ! end if ! sum ...

      end do ! p
    end do ! b

    close(123)


    ! Interpolation along a line

    open( 123, file = 'linediff.dat' )
    write( 123, '( "# origin(3) = ", g13.6 )' ) origin(3) * boinpm
    write( 123, '( "# gridwidth = ", g13.6 )' ) gridwidth * boinpm
    write( 123, '( "# maxind(3) = ", i0    )' ) maxind(3)

    coor(1) = 0
    coor(2) = 0
    coor(3) = origin(3) + gridwidth + 1E-10

    do while ( coor(3) <= origin(3) + (maxind(3)-1)*gridwidth )

      interp_pot = 0
      call interpolate_extended_hartree_non_periodic_ewald( coor, interp_pot )

      true_pot = 0
      call add_true_hartree_non_periodic_ewald( coor, true_pot, &
                                                l_hartree_max_far_distance, mindist=.false. )

      write(123,'( 3(1x,g13.6) )')  coor(3)*boinpm, true_pot*HinmV, interp_pot*HinmV

      coor(3) = coor(3) + 0.001

    end do

    close(123)

    ! stop   ! XX

  end subroutine hartree_non_periodic_ewald_tests
  !******







  !****s* hartree_non_periodic_ewald/calculate_2nd_deriv_from_splines
  !
  !  NAME
  !    calculate_2nd_deriv_from_splines
  !
  !  SYNOPSIS

  subroutine calculate_2nd_deriv_from_splines

    !  USES
    use spline

    !  AUTHOR
    !    Andris Gulans and Werner J"urgens
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE

    implicit none

    integer             :: i, j, k,  l, maxdim
    real*8, allocatable :: ys(:), spl_param(:,:)


    ! setting up the splines
    maxdim = max( maxind(1), maxind(2), maxind(3) )   ! XX rename
    allocate( ys(0:maxdim), spl_param(4,0:maxdim) )   ! XX check allocation


    ! x-direction
    do k = 0, maxind(3)
      do j = 0, maxind(2)

        ys( 0:maxind(1) ) = cartpot(:,j,k)
        call cubic_spline( ys, maxind(1)+1, spl_param )

        do l = 0, maxind(1)
          deriv2(1,l,j,k)  =  spline_2nd_deriv( real(l+1,8), spl_param, maxind(1)+1 )
        end do

      enddo
    enddo


    ! y-direction
    do k = 0, maxind(3)
      do i = 0, maxind(1)

        ys( 0:maxind(2) ) = cartpot(i,:,k)
        call cubic_spline( ys, maxind(2)+1, spl_param )

        do l = 0, maxind(2)
          deriv2(2,i,l,k)  =  spline_2nd_deriv( real(l+1,8), spl_param, maxind(2)+1 )
        end do

      enddo
    enddo


    ! z-direction
    do j = 0, maxind(2)
      do i = 0, maxind(1)

        ys( 0:maxind(3) ) = cartpot(i,j,:)
        call cubic_spline( ys, maxind(3)+1, spl_param )

        do l = 0, maxind(3)
          deriv2(3,i,j,l)  =  spline_2nd_deriv( real(l+1,8), spl_param, maxind(3)+1 )
        end do

      enddo
    enddo


    deallocate(ys)

  end subroutine calculate_2nd_deriv_from_splines
  !******







  !****s* hartree_non_periodic_ewald/splint
  !
  !  NAME
  !    splint
  !
  !  SYNOPSIS

  subroutine splint( ind, d, answer, grad )

    !  AUTHOR
    !    Andris Gulans
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e. V.
    !    Please note that any use of the FHI-aims-Software is subject to
    !    the terms and conditions of the respective license agreement.
    !
    !  SOURCE

    ! XX adapt the following description:
    ! uses spline interpolation to evaluate a function between grid points in 3D
    ! initialisation above is necessary to use this routine
    ! GRIDC   - a derived type variable used in VASP, the only relevant fields here are GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ
    ! density - the table of function to interpolate
    !           density(0,:,:,:) - the function itself
    !           density(1,:,:,:) - the second derivative w.r.t. x
    !           density(2,:,:,:) - the second derivative w.r.t. y
    !           density(3,:,:,:) - the second derivative w.r.t. z
    ! answer  - the interpolated function value
    ! grad    - interpolated first derivatives of the function


    implicit none

    integer, dimension(3), intent(in)            :: ind
    real*8,  dimension(3), intent(in)            :: d
    real*8,                intent(out)           :: answer
    real*8,  dimension(3), intent(out), optional :: grad

    real*8  a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3
    real*8  e1,e2,e3, f1,f2,f3
    real*8  t1,t2,t3, u1,u2,u3
    real*8  a1a2,a1a3,a2a3, b1b2,b2b3,b1b3, a1b2,a1b3,a2b3, b1a2,b1a3,b2a3
    real*8  sixth
    parameter(sixth=0.16666666666666666666666666666667d0)

    b1 = d(1)
    a1 = 1-b1
    c1 = (a1**3-a1)*sixth
    d1 = (b1**3-b1)*sixth

    b2 = d(2)
    a2 = 1-b2
    c2 = (a2**3-a2)*sixth
    d2 = (b2**3-b2)*sixth

    b3 = d(3)
    a3 = 1-b3
    c3 = (a3**3-a3)*sixth
    d3 = (b3**3-b3)*sixth


    answer =        a1*a2*a3*cartpot (ind(1),ind(2),ind(3))+b1*a2*a3*cartpot (ind(1)+1,ind(2),ind(3))
    answer = answer+c1*a2*a3*deriv2(1,ind(1),ind(2),ind(3))+d1*a2*a3*deriv2(1,ind(1)+1,ind(2),ind(3))
    answer = answer+a1*c2*a3*deriv2(2,ind(1),ind(2),ind(3))+b1*c2*a3*deriv2(2,ind(1)+1,ind(2),ind(3))
    answer = answer+a1*a2*c3*deriv2(3,ind(1),ind(2),ind(3))+b1*a2*c3*deriv2(3,ind(1)+1,ind(2),ind(3))

    answer = answer+a1*b2*a3*cartpot (ind(1),ind(2)+1,ind(3))+b1*b2*a3*cartpot (ind(1)+1,ind(2)+1,ind(3))
    answer = answer+c1*b2*a3*deriv2(1,ind(1),ind(2)+1,ind(3))+d1*b2*a3*deriv2(1,ind(1)+1,ind(2)+1,ind(3))
    answer = answer+a1*d2*a3*deriv2(2,ind(1),ind(2)+1,ind(3))+b1*d2*a3*deriv2(2,ind(1)+1,ind(2)+1,ind(3))
    answer = answer+a1*b2*c3*deriv2(3,ind(1),ind(2)+1,ind(3))+b1*b2*c3*deriv2(3,ind(1)+1,ind(2)+1,ind(3))

    answer = answer+a1*a2*b3*cartpot (ind(1),ind(2),ind(3)+1)+b1*a2*b3*cartpot (ind(1)+1,ind(2),ind(3)+1)
    answer = answer+c1*a2*b3*deriv2(1,ind(1),ind(2),ind(3)+1)+d1*a2*b3*deriv2(1,ind(1)+1,ind(2),ind(3)+1)
    answer = answer+a1*c2*b3*deriv2(2,ind(1),ind(2),ind(3)+1)+b1*c2*b3*deriv2(2,ind(1)+1,ind(2),ind(3)+1)
    answer = answer+a1*a2*d3*deriv2(3,ind(1),ind(2),ind(3)+1)+b1*a2*d3*deriv2(3,ind(1)+1,ind(2),ind(3)+1)

    answer = answer+a1*b2*b3*cartpot (ind(1),ind(2)+1,ind(3)+1)+b1*b2*b3*cartpot (ind(1)+1,ind(2)+1,ind(3)+1)
    answer = answer+c1*b2*b3*deriv2(1,ind(1),ind(2)+1,ind(3)+1)+d1*b2*b3*deriv2(1,ind(1)+1,ind(2)+1,ind(3)+1)
    answer = answer+a1*d2*b3*deriv2(2,ind(1),ind(2)+1,ind(3)+1)+b1*d2*b3*deriv2(2,ind(1)+1,ind(2)+1,ind(3)+1)
    answer = answer+a1*b2*d3*deriv2(3,ind(1),ind(2)+1,ind(3)+1)+b1*b2*d3*deriv2(3,ind(1)+1,ind(2)+1,ind(3)+1)



    if ( present(grad) ) then


      e1 = (3*a1**2-1)*sixth
      f1 = (3*b1**2-1)*sixth

      e2 = (3*a2**2-1)*sixth
      f2 = (3*b2**2-1)*sixth

      e3 = (3*a3**2-1)*sixth
      f3 = (3*b3**2-1)*sixth


      t1 = (a1**2-1)*sixth
      t2 = (a2**2-1)*sixth
      t3 = (a3**2-1)*sixth

      u1 = (b1**2-1)*sixth
      u2 = (b2**2-1)*sixth
      u3 = (b3**2-1)*sixth


      a1a2 = a1*a2
      a1a3 = a1*a3
      a2a3 = a2*a3

      b1b2 = b1*b2
      b1b3 = b1*b3
      b2b3 = b2*b3

      a1b2 = a1*b2
      a1b3 = a1*b3
      a2b3 = a2*b3

      b1a2 = b1*a2
      b1a3 = b1*a3
      b2a3 = b2*a3


      grad(1) =           -a2a3*cartpot (ind(1),ind(2),ind(3))+   a2a3*cartpot (ind(1)+1,ind(2),ind(3))
      grad(1) = grad(1)-e1*a2a3*deriv2(1,ind(1),ind(2),ind(3))+f1*a2a3*deriv2(1,ind(1)+1,ind(2),ind(3))
      grad(1) = grad(1)-t2*a2a3*deriv2(2,ind(1),ind(2),ind(3))+t2*a2a3*deriv2(2,ind(1)+1,ind(2),ind(3))
      grad(1) = grad(1)-t3*a2a3*deriv2(3,ind(1),ind(2),ind(3))+t3*a2a3*deriv2(3,ind(1)+1,ind(2),ind(3))

      grad(2) =           -a1a3*cartpot (ind(1),ind(2),ind(3))-   b1a3*cartpot (ind(1)+1,ind(2),ind(3))
      grad(2) = grad(2)-t1*a1a3*deriv2(1,ind(1),ind(2),ind(3))-u1*b1a3*deriv2(1,ind(1)+1,ind(2),ind(3))
      grad(2) = grad(2)-e2*a1a3*deriv2(2,ind(1),ind(2),ind(3))-e2*b1a3*deriv2(2,ind(1)+1,ind(2),ind(3))
      grad(2) = grad(2)-t3*a1a3*deriv2(3,ind(1),ind(2),ind(3))-t3*b1a3*deriv2(3,ind(1)+1,ind(2),ind(3))

      grad(3) =           -a1a2*cartpot (ind(1),ind(2),ind(3))-   b1a2*cartpot (ind(1)+1,ind(2),ind(3))
      grad(3) = grad(3)-t1*a1a2*deriv2(1,ind(1),ind(2),ind(3))-u1*b1a2*deriv2(1,ind(1)+1,ind(2),ind(3))
      grad(3) = grad(3)-t2*a1a2*deriv2(2,ind(1),ind(2),ind(3))-t2*b1a2*deriv2(2,ind(1)+1,ind(2),ind(3))
      grad(3) = grad(3)-e3*a1a2*deriv2(3,ind(1),ind(2),ind(3))-e3*b1a2*deriv2(3,ind(1)+1,ind(2),ind(3))


      grad(1) = grad(1)-   b2a3*cartpot (ind(1),ind(2)+1,ind(3))+   b2a3*cartpot (ind(1)+1,ind(2)+1,ind(3))
      grad(1) = grad(1)-e1*b2a3*deriv2(1,ind(1),ind(2)+1,ind(3))+f1*b2a3*deriv2(1,ind(1)+1,ind(2)+1,ind(3))
      grad(1) = grad(1)-u2*b2a3*deriv2(2,ind(1),ind(2)+1,ind(3))+u2*b2a3*deriv2(2,ind(1)+1,ind(2)+1,ind(3))
      grad(1) = grad(1)-t3*b2a3*deriv2(3,ind(1),ind(2)+1,ind(3))+t3*b2a3*deriv2(3,ind(1)+1,ind(2)+1,ind(3))

      grad(2) = grad(2)+   a1a3*cartpot (ind(1),ind(2)+1,ind(3))+   b1a3*cartpot (ind(1)+1,ind(2)+1,ind(3))
      grad(2) = grad(2)+t1*a1a3*deriv2(1,ind(1),ind(2)+1,ind(3))+u1*b1a3*deriv2(1,ind(1)+1,ind(2)+1,ind(3))
      grad(2) = grad(2)+f2*a1a3*deriv2(2,ind(1),ind(2)+1,ind(3))+f2*b1a3*deriv2(2,ind(1)+1,ind(2)+1,ind(3))
      grad(2) = grad(2)+t3*a1a3*deriv2(3,ind(1),ind(2)+1,ind(3))+t3*b1a3*deriv2(3,ind(1)+1,ind(2)+1,ind(3))

      grad(3) = grad(3)-   a1b2*cartpot (ind(1),ind(2)+1,ind(3))-   b1b2*cartpot (ind(1)+1,ind(2)+1,ind(3))
      grad(3) = grad(3)-t1*a1b2*deriv2(1,ind(1),ind(2)+1,ind(3))-u1*b1b2*deriv2(1,ind(1)+1,ind(2)+1,ind(3))
      grad(3) = grad(3)-u2*a1b2*deriv2(2,ind(1),ind(2)+1,ind(3))-u2*b1b2*deriv2(2,ind(1)+1,ind(2)+1,ind(3))
      grad(3) = grad(3)-e3*a1b2*deriv2(3,ind(1),ind(2)+1,ind(3))-e3*b1b2*deriv2(3,ind(1)+1,ind(2)+1,ind(3))


      grad(1) = grad(1)-   a2b3*cartpot (ind(1),ind(2),ind(3)+1)+   a2b3*cartpot (ind(1)+1,ind(2),ind(3)+1)
      grad(1) = grad(1)-e1*a2b3*deriv2(1,ind(1),ind(2),ind(3)+1)+f1*a2b3*deriv2(1,ind(1)+1,ind(2),ind(3)+1)
      grad(1) = grad(1)-t2*a2b3*deriv2(2,ind(1),ind(2),ind(3)+1)+t2*a2b3*deriv2(2,ind(1)+1,ind(2),ind(3)+1)
      grad(1) = grad(1)-u3*a2b3*deriv2(3,ind(1),ind(2),ind(3)+1)+u3*a2b3*deriv2(3,ind(1)+1,ind(2),ind(3)+1)

      grad(2) = grad(2)-   a1b3*cartpot (ind(1),ind(2),ind(3)+1)-   b1b3*cartpot (ind(1)+1,ind(2),ind(3)+1)
      grad(2) = grad(2)-t1*a1b3*deriv2(1,ind(1),ind(2),ind(3)+1)-u1*b1b3*deriv2(1,ind(1)+1,ind(2),ind(3)+1)
      grad(2) = grad(2)-e2*a1b3*deriv2(2,ind(1),ind(2),ind(3)+1)-e2*b1b3*deriv2(2,ind(1)+1,ind(2),ind(3)+1)
      grad(2) = grad(2)-u3*a1b3*deriv2(3,ind(1),ind(2),ind(3)+1)-u3*b1b3*deriv2(3,ind(1)+1,ind(2),ind(3)+1)

      grad(3) = grad(3)+   a1a2*cartpot (ind(1),ind(2),ind(3)+1)+   b1a2*cartpot (ind(1)+1,ind(2),ind(3)+1)
      grad(3) = grad(3)+t1*a1a2*deriv2(1,ind(1),ind(2),ind(3)+1)+u1*b1a2*deriv2(1,ind(1)+1,ind(2),ind(3)+1)
      grad(3) = grad(3)+t2*a1a2*deriv2(2,ind(1),ind(2),ind(3)+1)+t2*b1a2*deriv2(2,ind(1)+1,ind(2),ind(3)+1)
      grad(3) = grad(3)+f3*a1a2*deriv2(3,ind(1),ind(2),ind(3)+1)+f3*b1a2*deriv2(3,ind(1)+1,ind(2),ind(3)+1)


      grad(1) = grad(1)-   b2b3*cartpot (ind(1),ind(2)+1,ind(3)+1)+   b2b3*cartpot (ind(1)+1,ind(2)+1,ind(3)+1)
      grad(1) = grad(1)-e1*b2b3*deriv2(1,ind(1),ind(2)+1,ind(3)+1)+f1*b2b3*deriv2(1,ind(1)+1,ind(2)+1,ind(3)+1)
      grad(1) = grad(1)-u2*b2b3*deriv2(2,ind(1),ind(2)+1,ind(3)+1)+u2*b2b3*deriv2(2,ind(1)+1,ind(2)+1,ind(3)+1)
      grad(1) = grad(1)-u3*b2b3*deriv2(3,ind(1),ind(2)+1,ind(3)+1)+u3*b2b3*deriv2(3,ind(1)+1,ind(2)+1,ind(3)+1)

      grad(2) = grad(2)+   a1b3*cartpot (ind(1),ind(2)+1,ind(3)+1)+   b1b3*cartpot (ind(1)+1,ind(2)+1,ind(3)+1)
      grad(2) = grad(2)+t1*a1b3*deriv2(1,ind(1),ind(2)+1,ind(3)+1)+u1*b1b3*deriv2(1,ind(1)+1,ind(2)+1,ind(3)+1)
      grad(2) = grad(2)+f2*a1b3*deriv2(2,ind(1),ind(2)+1,ind(3)+1)+f2*b1b3*deriv2(2,ind(1)+1,ind(2)+1,ind(3)+1)
      grad(2) = grad(2)+u3*a1b3*deriv2(3,ind(1),ind(2)+1,ind(3)+1)+u3*b1b3*deriv2(3,ind(1)+1,ind(2)+1,ind(3)+1)

      grad(3) = grad(3)+   a1b2*cartpot (ind(1),ind(2)+1,ind(3)+1)+   b1b2*cartpot (ind(1)+1,ind(2)+1,ind(3)+1)
      grad(3) = grad(3)+t1*a1b2*deriv2(1,ind(1),ind(2)+1,ind(3)+1)+u1*b1b2*deriv2(1,ind(1)+1,ind(2)+1,ind(3)+1)
      grad(3) = grad(3)+u2*a1b2*deriv2(2,ind(1),ind(2)+1,ind(3)+1)+u2*b1b2*deriv2(2,ind(1)+1,ind(2)+1,ind(3)+1)
      grad(3) = grad(3)+f3*a1b2*deriv2(3,ind(1),ind(2)+1,ind(3)+1)+f3*b1b2*deriv2(3,ind(1)+1,ind(2)+1,ind(3)+1)


    end if  ! present(grad)



  end subroutine splint
  !******




end module hartree_non_periodic_ewald
!******
