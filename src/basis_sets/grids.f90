!****h* FHI-aims/grids
!  NAME
!    grids
!  SYNOPSIS
      module grids
! PURPOSE
!  Module grids handles everything related to the real-space grids:
!  * the logarithmic radial grids for free-atom related quantities
!  * the radial integration grids
!  * the Lebedev-style angular grids
!
!  Subroutines:
!  * allocate_grids
!  * get_grids
!  * get_lebedev - historic and no longer used
!  * test_radial_grid
!  * get_logarithmic_grid
!  * cleanup_grids
!
!  Functions:
!  * invert_radial_grid
!  * invert_log_grid
!  * get_radial_weight
!
!  global variable declarations - exported to other program parts:

!  *   r_grid_min: for each species, the innermost and
!  *   r_grid_max: and outermost points of the logarithmic grid for wave
!  *   r_grid_inc: functions, and the increase factor
!                 !!! Notice that r_grid_min = r_grid_min/species_z !!!
!                 on (Lebedev) angular integration grids
!  *   n_grid    : number of radial grid points for atomic(!) potential, density,
!                 and radial wave functions for each species
!                 Please note that
!                    r_grid_max(i) ~ r_grid_min * r_grid_inc**(n_grid-1)
!                 is only approximately true.
!  *   r_grid    : Radial grid points for atomic(!) potential, density,
!                 radial wave functions - for each species

!  *   n_radial : for each species, number of points for radial integration grids
!  *   scale_radial : for each species, scale factor for radial integration grids
!  *   angular_acc : requested integration accuracy per integration shell
!  *   angular_limit : maximum number of angular integration points per shell for a given species
!  *   r_radial  : for each species, radial integration grid points
!  *   w_radial  : for each species, radial integration grid weights
   
!  *   n_angular : for each atom and each radial shell(!), actual number of grid points (Angular grids must be allocated per atom.)
!  *   r_angular : for each atom, angular integration grid points
!  *   w_angular : for each atom, angular integration grid weights
!
!  USES
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
      implicit none

      real*8,  dimension(:),       allocatable :: r_grid_min            ! r_grid_min(n_species)
      real*8,  dimension(:),       allocatable :: r_grid_max            ! r_grid_max(n_species)
      real*8,  dimension(:),       allocatable :: r_grid_inc            ! r_grid_inc(n_species)
      real*8,  dimension(:),       allocatable :: log_r_grid_inc        ! r_grid_inc(n_species)
      integer, dimension(:),       allocatable :: n_grid                ! n_grid(n_species)
      real*8,  dimension(:,:),     allocatable :: r_grid                ! r_grid(n_max_grid, n_species)
      integer, dimension(:),       allocatable :: n_radial              ! n_radial(n_species)
      real*8,  dimension(:),       allocatable :: scale_radial          ! scale_radial(n_species)
      integer, dimension(:),       allocatable :: angular_limit         ! angular_limit(n_species)
      integer, dimension(:),       allocatable :: angular_min           ! angular_min(n_species)
      real*8,  dimension(:),       allocatable :: angular_acc           ! angular_acc(n_species)
      real*8,  dimension(:,:),     allocatable :: r_radial              ! r_radial(n_max_radial, n_species)
      real*8,  dimension(:,:),     allocatable :: w_radial              ! w_radial(n_max_radial, n_species )
      integer, dimension(:,:),     allocatable :: n_angular             ! n_angular(n_max_radial, n_species)
      real*8,  dimension(:,:,:,:), allocatable :: r_angular             ! r_angular(3, n_max_angular, n_max_radial, n_species)
      real*8,  dimension(:,:,:),   allocatable :: w_angular             ! w_angular(n_max_angular, n_max_radial, n_species)
      integer,    dimension(:),     allocatable :: n_ang_lebedev         !for nlcorr. SAG
      real*8,  dimension(:,:,:),   allocatable :: r_ang_lebedev
      real*8,  dimension(:,:),     allocatable :: w_ang_lebedev
          
      integer, dimension(:),       allocatable :: n_angular_lebedev
      real*8,  dimension(:,:,:),   allocatable :: r_angular_lebedev
      real*8,  dimension(:,:),     allocatable :: w_angular_lebedev
      integer, dimension(:),       allocatable :: n_division_lebedev
      integer, dimension(:,:),     allocatable :: division_boundaries_lebedev
      integer, dimension(:,:),     allocatable :: fixed_grid_index
      integer, dimension(:,:),     allocatable :: lebedev_grid_index
      integer, dimension(:,:),     allocatable :: n_division            ! n_division(n_max_radial, n_species)
      integer, dimension(:,:,:),   allocatable :: division_boundaries   ! division_boundaries(n_max_angular_division+1,
      real*8,  dimension(:,:,:),   allocatable :: local_ylm_tab

      integer       :: n_max_lebedev
      integer       :: n_grids
      logical       :: grid_partitioned

      type grid_point
        real*8, dimension(3) :: coords
        integer :: index_atom
        integer :: index_radial
        integer :: index_angular
      end type grid_point

      type batch_of_points
        integer :: size
!        type(grid_point), dimension(:), allocatable :: points
        type(grid_point), pointer, dimension(:) :: points
        integer :: batch_n_compute
        integer, pointer, dimension(:) :: batch_i_basis
      end type batch_of_points

      type(batch_of_points), pointer, dimension(:) :: batches

!******

      contains
        subroutine set_grids_defaults( )
        implicit none

        grid_partitioned = .false.
        end subroutine set_grids_defaults

!---------------------------------------------------------------------------
!****s* grids/allocate_grids
!  NAME
!    allocate_grids
!  SYNOPSIS
        subroutine allocate_grids( )
!  PURPOSE
!  Subroutine allocate_grids allocates the necessary memory for all real-space grids.
!  This allocation is done once at the beginning of the code, since the grids
!  are needed everywhere.
!
!  In principle one could allocate everything in a much more fine-grained way - but only
!  if there is a clear need.
!  USES
        use dimensions,   only : n_max_angular, n_max_angular_division, n_max_grid, &
                                 n_max_radial, n_species
! ARGUMENTS
!  INPUT
!   none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
        implicit none

        n_max_lebedev = 50
        allocate ( r_grid_min(n_species) )
        allocate ( r_grid_max(n_species) )
        allocate ( r_grid_inc(n_species) )
        allocate ( log_r_grid_inc(n_species) )
        allocate ( n_grid(n_species) )
        allocate ( r_grid(n_max_grid, n_species) )
        allocate ( n_radial(n_species) )
        allocate ( scale_radial(n_species) )
        allocate ( angular_limit(n_species) )
        allocate ( angular_min(n_species) )
        allocate ( angular_acc(n_species) )
	  allocate ( r_radial(n_max_radial, n_species) )
	  allocate ( w_radial(n_max_radial, n_species) )
	  allocate ( n_angular(n_max_radial, n_species) )
	  allocate ( r_angular(3, n_max_angular, n_max_radial, n_species))
	  allocate ( w_angular(n_max_angular, n_max_radial, n_species) )
	  allocate ( n_division(n_max_radial, n_species) )
	  allocate ( division_boundaries(n_max_angular_division+1, n_max_radial, n_species) )
        allocate(n_angular_lebedev(n_max_lebedev))
        allocate(r_angular_lebedev(3,n_max_angular,n_max_lebedev))
        allocate(w_angular_lebedev(n_max_angular,n_max_lebedev))
        allocate(n_division_lebedev(n_max_lebedev))
        allocate ( division_boundaries_lebedev(n_max_angular_division+1,n_max_lebedev) )
	allocate(lebedev_grid_index(n_max_radial, n_species))
        
        
        
        allocate(n_ang_lebedev(15))  !SAG
        allocate(r_ang_lebedev(3,350,15))
        allocate(w_ang_lebedev(350,15))
        

        end subroutine allocate_grids
!******
!------------------------------------------------------------------------
!****s* grids/get_grids
!  NAME
!    get_grids
!  SYNOPSIS
        subroutine get_grids( out_grids )
!  PURPOSE
!    set up the grids for various species
!  USES
        use dimensions,  only : n_max_points_per_div, n_species, &
                                n_max_angular, use_vdw_method, &
                                use_vdw_post, use_nlcorr_post
        use mpi_tasks,   only : myid 
        use localorb_io, only : use_unit, localorb_info
        implicit none

!  ARGUMENTS
        logical :: out_grids
!  INPUTS
!   o out_grids -- print grids to files or not?
!  OUTPUT
!   none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).

!  SOURCE

        real*8 :: r_scaled
        character*15 :: out_file
        integer :: n_lebedev
        integer :: try_n_angular_lebedev
        real*8 :: try_r_angular_lebedev(3,n_max_angular)
        real*8 :: try_w_angular_lebedev(n_max_angular)
        integer :: i_species, i_radial
        integer :: i_lebedev, i_lebedev_copy, i_lebedev_prev, i_grid
        integer :: i_division
        integer :: i_leb, nal  !SAG
        real*8 :: ral(3,350)
        real*8 :: wal(350)

        call localorb_info( &
             "Setting up grids for atomic and cluster calculations.", &
             use_unit,'(2X,A)')

!     allocate and initialize array for fixed grid data
!        if (any(species_fixed_grid)) then
!           allocate(fixed_grid_index(n_max_radial, n_species))
!           fixed_grid_index = 0
!        end if

!       initialize logarithmic grid size

        do i_species = 1, n_species, 1
          n_grid(i_species) = 0
        enddo

!       set up grids species by species



	 do i_species = 1, n_species,1

! ----->  set up the radial grid

	    do i_radial = 1, n_radial (i_species)

      !           scale radial grid by n
		r_scaled = dble(i_radial/(n_radial(i_species)+1.0d0))

    !           get radial grid point
		r_radial (i_radial, i_species ) = -log(1.d0-r_scaled**2.d0)

		r_radial (i_radial, i_species) = &
		  scale_radial(i_species)*r_radial(i_radial, i_species)

    !           get radial grid weight dr/di
		w_radial (i_radial, i_species) = &
		(2.0d0/dble(i_radial)) / ( 1.0d0/(r_scaled**2.d0) - 1.0d0 )

		w_radial (i_radial, i_species) = &
		scale_radial(i_species)*w_radial(i_radial, i_species)

            enddo

!  Check if minimal logarithmic grid point is smaller than minimal
!  radial integration grid point; if not, adjust.

          if (r_grid_min(i_species).gt.r_radial(1,i_species)) then
             if (myid.eq.0) then
                write(use_unit,*)
                write(use_unit,*) "* Minimum radial grid point for ", &
                     "logarithmic grid of species ", i_species
                write(use_unit,*) "* , r_min = ", &
                     r_grid_min(i_species), ","
                write(use_unit,*) "* is chosen above the minimum radial", &
                     " integration grid point,"
                write(use_unit,*) "* r_min = ", r_radial(i_species,1), "."
             endif

             r_grid_min(i_species) = 0.5d0*r_radial(1,i_species)

             if (myid.eq.0) then
                write(use_unit,*) "* Setting r_grid_min to ", &
                     r_grid_min(i_species), "."
                write(use_unit,*)
             end if

          end if

          call get_logarithmic_grid &
          ( i_species )

          log_r_grid_inc(i_species) = log(r_grid_inc(i_species))

        enddo

        
!       Initialize angular grids separately.
!       Angular grids are also self-adapting, and reset in initialize_integrals.

        do i_species = 1, n_species, 1

           do i_radial = 1, n_radial(i_species)

              n_lebedev = &
                   lebedev_grid_floor(angular_limit(i_species))
              
                   
              call get_angular_grid &
                   ( n_lebedev, &
                   n_angular(i_radial,i_species), &
                   r_angular(1,1,i_radial,i_species), &
                   w_angular(1,i_radial,i_species) &
                   )
                   
             
              if (n_angular(i_radial,i_species).gt.n_max_angular) then
!             smallest grid is already too large, abort calculation
                 if (myid.eq.0) then
                    write(use_unit,'(1X,A,I5,A,I5,A)') &
                         "* Species ", i_species, " radial shell ", &
                         i_radial, ":"
                    write(use_unit,'(1X,A)') &
                         "* No suitable angular grid available. "
                    write(use_unit,'(1X,A,A)') &
                       "* Adjust the prescribed number of grid points ", &
                       "in control.in."
                 end if
                 stop

              end if

              call divide_angular_grid( &
                   n_angular(i_radial, i_species), &
                   r_angular(1,1,i_radial,i_species), &
                   w_angular(1,i_radial,i_species), &
                   n_division(i_radial, i_species), &
                   division_boundaries(1,i_radial,i_species), &
                   1.0d8 )
                   
           enddo

        enddo

        ! next, tabulate Lebedev grids themselved for the possibility of self-adapting
        ! grids later on. This option is not the default, but since tabulating the Lebedev 
        ! grids does not cost anything, we retain it here.

        ! initialize unused parts of auxiliary arrays to avoid suprious floating exceptions later
        try_r_angular_lebedev = 0.d0
        try_w_angular_lebedev = 0.d0

        n_lebedev = 0

        do i_species = 1, n_species, 1

           n_lebedev = MAX(n_lebedev, &
                lebedev_grid_floor(angular_limit(i_species)))
        enddo

        n_max_points_per_div = 0

        i_lebedev_prev = 0
        i_grid = 0
        do i_lebedev = 1, n_lebedev, 1

           i_lebedev_copy = i_lebedev

           call get_angular_grid &
                ( i_lebedev_copy, &
                try_n_angular_lebedev, &
                try_r_angular_lebedev, &
                try_w_angular_lebedev &
                )


           if(i_lebedev_copy.ne.i_lebedev_prev) then

              i_grid = i_grid + 1

              n_angular_lebedev(i_grid) = try_n_angular_lebedev
              r_angular_lebedev(:,:,i_grid) = try_r_angular_lebedev
              w_angular_lebedev(:,i_grid) = try_w_angular_lebedev

              call divide_angular_grid_p0 &
                   ( n_angular_lebedev(i_grid), &
                   r_angular_lebedev(1,1,i_grid), &
                   w_angular_lebedev(1,i_grid), &
                   n_division_lebedev(i_grid), &
                   division_boundaries_lebedev(1,i_grid), &
                   n_angular_lebedev(i_grid) &
                   )

              i_lebedev_prev = i_lebedev_copy

              do i_division = 1, n_division_lebedev(i_grid), 1
                n_max_points_per_div = &
                  max( n_max_points_per_div, &
                       division_boundaries_lebedev(i_division+1,i_grid) &
                      -division_boundaries_lebedev(i_division,i_grid) )
              enddo

           end if

        enddo

        n_grids = i_grid

        if (out_grids) then
!         output option for radial grid
          do i_species = 1, n_species,1
             if (myid.eq.0) then
                if (i_species.le.9) then
                   write(out_file,'(A9,I1,A4)') "rad_grid_", &
                        i_species,".dat"
                else if (i_species.le.99) then
                   write(out_file,'(A9,I2,A4)') "rad_grid_", &
                        i_species,".dat"
                end if
                open(50,file=out_file)
                do i_radial = 1, n_radial (i_species)
                   write(50,*) r_radial(i_radial,i_species), &
                        w_radial(i_radial,i_species)
                enddo
                close(50)
             end if
          enddo
       end if
       
       
       
       !SAG: make lebedev grids available for nlcorr
       !right now only one of these is actually used...
       !      try_r_angular_lebedev = 0.d0
       !      try_w_angular_lebedev = 0.d0
       !      try_n_angular_lebedev = 0.d0
       
       if(use_vdw_method.or.use_vdw_post.or.use_nlcorr_post)then
          do i_leb = 1, 15       
             ! call get_angular_grid(i_leb, try_n_angular_lebedev, &  !wouldn't work!
             !try_r_angular_lebedev, try_w_angular_lebedev )
             !n_ang_lebedev(i_leb) = try_n_angular_lebedev
             !r_ang_lebedev(:,:,i_leb) = try_r_angular_lebedev(:,1:350)
             !w_ang_lebedev(:,i_leb) = try_w_angular_lebedev(1:350)
             
              call get_leb_grids(i_leb, nal,ral,wal,350)
             n_ang_lebedev(i_leb) = nal
             r_ang_lebedev(:,:,i_leb) = ral
             w_ang_lebedev(:,i_leb) = wal
             
             
             
          end do
          
       endif






     end subroutine get_grids
!******
!---------------------------------------------------------------------------------------
!****s* grids/get_lebedev
!  NAME
!    get_atomic_occ_numbers
!  SYNOPSIS
        subroutine get_lebedev( n_ang, i_species, i_radial )
! PURPOSE
!  obtains the Lebedev integration grid points for given angular resolution.
!
!  This version is a legacy version and no longer used; we now use
!  subroutine get_angular_grids for the same purpose
!
!   USES
        use dimensions,  only : n_max_angular
        use mpi_tasks,   only : myid
        use localorb_io, only : use_unit, localorb_info
        implicit none
! ARGUMENTS

        integer :: n_ang
        integer :: i_species 
        integer :: i_radial

!  INPUTS
!    o n_ang     -- desired number of grid points in shell 
!    o i_species -- species number in question
!    o i_radial  -- radial shell in question 
!
!  OUTPUTS 
!    none
!  SOURCE
        real*8  :: x(n_max_angular)
        real*8  :: y(n_max_angular)
        real*8  :: z(n_max_angular)
        real*8  :: w(n_max_angular)
        integer :: n_leb
        integer :: i_ang


        if (n_ang.lt.6) then
           call localorb_info("* No suitable Lebedev grid found.", use_unit)
           stop
        else if (n_ang.lt.14) then
          n_ang = 6
          call LD0006(x,y,z,w,n_leb)
        else if (n_ang.lt.26) then
          n_ang = 14
          call LD0014(x,y,z,w,n_leb)
        else if (n_ang.lt.38) then
          n_ang = 26
          call LD0026(x,y,z,w,n_leb)
        else if (n_ang.lt.50) then
          n_ang = 38
          call LD0038(x,y,z,w,n_leb)
!patch        else if (n_ang.lt.74) then
        else if (n_ang.lt.86) then
          n_ang = 50
          call LD0050(x,y,z,w,n_leb)
!patch        else if (n_ang.lt.86) then
!patch          n_ang = 74
!patch          call LD0074(x,y,z,w,n_leb)
        else if (n_ang.lt.110) then
          n_ang = 86
          call LD0086(x,y,z,w,n_leb)
        else if (n_ang.lt.146) then
          n_ang = 110
          call LD0110(x,y,z,w,n_leb)
        else if (n_ang.lt.170) then
          n_ang = 146
          call LD0146(x,y,z,w,n_leb)
        else if (n_ang.lt.194) then
          n_ang = 170
          call LD0170(x,y,z,w,n_leb)
!patch        else if (n_ang.lt.230) then
        else if (n_ang.lt.302) then
          n_ang = 194
          call LD0194(x,y,z,w,n_leb)
!patch        else if (n_ang.lt.266) then
!patch          n_ang = 230
!patch          call LD0230(x,y,z,w,n_leb)
!patch        else if (n_ang.lt.302) then
!patch          n_ang = 266
!patch          call LD0266(x,y,z,w,n_leb)
        else if (n_ang.lt.350) then
          n_ang = 302
          call LD0302(x,y,z,w,n_leb)
        else if (n_ang.lt.434) then
          n_ang = 350
          call LD0350(x,y,z,w,n_leb)
        else if (n_ang.lt.590) then
          n_ang = 434
          call LD0434(x,y,z,w,n_leb)
        else if (n_ang.lt.770) then
          n_ang = 590
          call LD0590(x,y,z,w,n_leb)
        else if (n_ang.lt.974) then
          n_ang = 770
          call LD0770(x,y,z,w,n_leb)
        else if (n_ang.lt.1202) then
          n_ang = 974
          call LD0974(x,y,z,w,n_leb)
        else if (n_ang.lt.1454) then
          n_ang = 1202
          call LD1202(x,y,z,w,n_leb)
        else if (n_ang.lt.1730) then
          n_ang = 1454
          call LD1454(x,y,z,w,n_leb)
        else if (n_ang.lt.2030) then
          n_ang = 1730
          call LD1730(x,y,z,w,n_leb)
        else if (n_ang.lt.2354) then
          n_ang = 2030
          call LD2030(x,y,z,w,n_leb)
        else if (n_ang.lt.2702) then
          n_ang = 2354
          call LD2354(x,y,z,w,n_leb)
        else if (n_ang.lt.3074) then
          n_ang = 2702
          call LD2702(x,y,z,w,n_leb)
        else if (n_ang.lt.3470) then
          n_ang = 3074
          call LD3074(x,y,z,w,n_leb)
        else if (n_ang.lt.3890) then
          n_ang = 3470
          call LD3470(x,y,z,w,n_leb)
        else if (n_ang.lt.4334) then
          n_ang = 3890
          call LD3890(x,y,z,w,n_leb)
        else if (n_ang.lt.4802) then
          n_ang = 4334
          call LD4334(x,y,z,w,n_leb)
        else if (n_ang.lt.5294) then
          n_ang = 4802
          call LD4802(x,y,z,w,n_leb)
        else if (n_ang.lt.5810) then
          n_ang = 5294
          call LD5294(x,y,z,w,n_leb)
        else
         if (n_ang.gt.5810) then
           if (myid.eq.0) then
             write(use_unit,'(1X,A,A)') &
             "! Defaulting to maximum number of angular grid points: ", &
                     "5810."
           end if
           n_ang = 5810
         end if
         call LD5810(x,y,z,w,n_leb)
        end if

!test
!          write(use_unit,*) "Grid chosen: n_ang = ", n_ang
!test end

        do i_ang = 1, n_ang, 1
!test
!          write(use_unit,*) "i_ang   = ", i_ang
!          write(use_unit,*) "x(i_ang)= ", x(i_ang)
!          write(use_unit,*) "y(i_ang)= ", y(i_ang)
!          write(use_unit,*) "z(i_ang)= ", z(i_ang)
!          write(use_unit,*) "w(i_ang)= ", w(i_ang)
!test end
          r_angular (1, i_ang, i_radial, i_species) = x(i_ang)
          r_angular (2, i_ang, i_radial, i_species) = y(i_ang)
          r_angular (3, i_ang, i_radial, i_species) = z(i_ang)
          w_angular (i_ang, i_radial, i_species)    = w(i_ang)
!test
!          write(use_unit,*) "Next i_ang."
!test end

        enddo

!test
!          write(use_unit,*) "Leave get_lebedev."
!test end

        end subroutine get_lebedev
!******

!------------------------------------------------------------------------------
!****f* grids/lebedev_grid_ceil
!  NAME
!    lebedev_grid_ceil
!  SYNOPSIS
        integer function lebedev_grid_ceil ( n_ang )
! PURPOSE
!  function lebedev_grid_ceil obtains the next higher number
!  of the Lebedev grid for a given number
!  of angular integration grid points.
! USES

        use runtime_choices, only : force_lebedev
        use mpi_tasks,       only : myid
        use localorb_io,     only : use_unit
!  ARGUMENTS
        integer :: n_ang
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o n_ang - number of integration points

!  OUTPUTS
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
        if (force_lebedev.eq.2) then
            if (n_ang.le.0) then
                lebedev_grid_ceil = 0
            else if (n_ang.le.50) then
                lebedev_grid_ceil = 1
            else if (n_ang.le.62) then
                lebedev_grid_ceil = 2
            else if (n_ang.le.72) then
                lebedev_grid_ceil = 3
            else if (n_ang.le.86) then
                lebedev_grid_ceil = 4
            else if (n_ang.le.98) then
                lebedev_grid_ceil = 5
            else if (n_ang.le.114) then
                lebedev_grid_ceil = 6
            else if (n_ang.le.128) then
                lebedev_grid_ceil = 7
            else if (n_ang.le.146) then
                lebedev_grid_ceil = 8
            else if (n_ang.le.163) then
                lebedev_grid_ceil = 9
            else if (n_ang.le.182) then
                lebedev_grid_ceil = 10
            else if (n_ang.le.201) then
                lebedev_grid_ceil = 11
            else if (n_ang.le.222) then
                lebedev_grid_ceil = 12
            else if (n_ang.le.243) then
                lebedev_grid_ceil = 13
            else if (n_ang.le.266) then
                lebedev_grid_ceil = 14
            else if (n_ang.le.289) then
                lebedev_grid_ceil = 15
            else if (n_ang.le.314) then
                lebedev_grid_ceil = 16
            else if (n_ang.le.339) then
                lebedev_grid_ceil = 17
            else if (n_ang.le.366) then
                lebedev_grid_ceil = 18
            else if (n_ang.le.393) then
                lebedev_grid_ceil = 19
            else if (n_ang.le.422) then
                lebedev_grid_ceil = 20
            else if (n_ang.le.451) then
                lebedev_grid_ceil = 21
            else if (n_ang.le.482) then
                lebedev_grid_ceil = 22
            else if (n_ang.le.513) then
                lebedev_grid_ceil = 23
            else if (n_ang.le.546) then
                lebedev_grid_ceil = 24
            else if (n_ang.le.579) then
                lebedev_grid_ceil = 25
            else if (n_ang.le.614) then
                lebedev_grid_ceil = 26
            else if (n_ang.le.649) then
                lebedev_grid_ceil = 27
            else if (n_ang.le.686) then
                lebedev_grid_ceil = 28
            else if (n_ang.le.723) then
                lebedev_grid_ceil = 29
            else if (n_ang.le.762) then
                lebedev_grid_ceil = 30
            else if (n_ang.le.801) then
                lebedev_grid_ceil = 31
            else if (n_ang.le.842) then
                lebedev_grid_ceil = 32
            else if (n_ang.le.883) then
                lebedev_grid_ceil = 33
            else if (n_ang.le.926) then
                lebedev_grid_ceil = 34
            else if (n_ang.le.969) then
                lebedev_grid_ceil = 35
            else if (n_ang.le.1014) then
                lebedev_grid_ceil = 36
            else if (n_ang.le.1059) then
                lebedev_grid_ceil = 37
            else if (n_ang.le.1106) then
                lebedev_grid_ceil = 38
            else if (n_ang.le.1153) then
                lebedev_grid_ceil = 39
            else if (n_ang.le.1202) then
                lebedev_grid_ceil = 40
            else if (n_ang.le.1251) then
                lebedev_grid_ceil = 41
            else if (n_ang.le.1302) then
                lebedev_grid_ceil = 42
            else if (n_ang.le.1862) then
                lebedev_grid_ceil = 43
            else if (n_ang.le.3282) then
                lebedev_grid_ceil = 44
            else if (n_ang.le.5102) then
                lebedev_grid_ceil = 45
            else if (n_ang.le.7322) then
                lebedev_grid_ceil = 46
            else if (n_ang.le.9942) then
                lebedev_grid_ceil = 47
            else if (n_ang.le.12962) then
                lebedev_grid_ceil = 48
            else if (n_ang.le.16382) then
                lebedev_grid_ceil = 49
            else
                lebedev_grid_ceil = 49
            end if
        else if ((force_lebedev.eq.1).or.(n_ang.gt.1202)) then
          if (n_ang.le.0) then
            lebedev_grid_ceil = 0
          else if (n_ang.le.6) then
            lebedev_grid_ceil = 1
          else if (n_ang.le.14) then
            lebedev_grid_ceil = 2
          else if (n_ang.le.26) then
            lebedev_grid_ceil = 3
          else if (n_ang.le.38) then
            lebedev_grid_ceil = 4
          else if (n_ang.le.50) then
            lebedev_grid_ceil = 5
          else if (n_ang.le.78) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_ceil = 6
          else if (n_ang.le.86) then
            lebedev_grid_ceil = 7
          else if (n_ang.le.110) then
            lebedev_grid_ceil = 8
          else if (n_ang.le.146) then
            lebedev_grid_ceil = 9
          else if (n_ang.le.170) then
            lebedev_grid_ceil = 10
          else if (n_ang.le.194) then
            lebedev_grid_ceil = 11
          else if (n_ang.le.248) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_ceil = 12
          else if (n_ang.le.288) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_ceil = 13
          else if (n_ang.le.302) then
            lebedev_grid_ceil = 14
          else if (n_ang.le.350) then
            lebedev_grid_ceil = 15
          else if (n_ang.le.434) then
            lebedev_grid_ceil = 16
          else if (n_ang.le.590) then
            lebedev_grid_ceil = 17
          else if (n_ang.le.770) then
            lebedev_grid_ceil = 18
          else if (n_ang.le.974) then
            lebedev_grid_ceil = 19
          else if (n_ang.le.1202) then
            lebedev_grid_ceil = 20
          else if (n_ang.le.1454) then
            lebedev_grid_ceil = 21
          else if (n_ang.le.1730) then
            lebedev_grid_ceil = 22
          else if (n_ang.le.2030) then
            lebedev_grid_ceil = 23
          else if (n_ang.le.2354) then
            lebedev_grid_ceil = 24
          else if (n_ang.le.2702) then
            lebedev_grid_ceil = 25
          else if (n_ang.le.3074) then
            lebedev_grid_ceil = 26
          else if (n_ang.le.3470) then
            lebedev_grid_ceil = 27
          else if (n_ang.le.3890) then
            lebedev_grid_ceil = 28
          else if (n_ang.le.4334) then
            lebedev_grid_ceil = 29
          else if (n_ang.le.4802) then
            lebedev_grid_ceil = 30
          else if (n_ang.le.5294) then
            lebedev_grid_ceil = 31
          else if (n_ang.le.5810) then
            lebedev_grid_ceil = 32
          else
           if (n_ang.gt.5810) then
             if (myid.eq.0) then
               write(use_unit,'(1X,A,A)') &
               "! Defaulting to maximum number of angular grid points: ", &
               "5810."
             end if
             lebedev_grid_ceil = 32
           end if
          end if
        else
!         Delley style grids; not all are supplied
          if (n_ang.le.0) then
            lebedev_grid_ceil = 0
          else if (n_ang.le.6) then
            lebedev_grid_ceil = 1
          else if (n_ang.le.14) then
            lebedev_grid_ceil = 2
          else if (n_ang.le.26) then
            lebedev_grid_ceil = 3
          else if (n_ang.le.50) then
            lebedev_grid_ceil = 5
          else if (n_ang.le.110) then
            lebedev_grid_ceil = 8
          else if (n_ang.le.194) then
            lebedev_grid_ceil = 11
          else if (n_ang.le.302) then
            lebedev_grid_ceil = 14
          else if (n_ang.le.434) then
            lebedev_grid_ceil = 16
          else if (n_ang.le.590) then
            lebedev_grid_ceil = 17
          else if (n_ang.le.770) then
            lebedev_grid_ceil = 18
          else if (n_ang.le.974) then
            lebedev_grid_ceil = 19
          else if (n_ang.le.1202) then
            lebedev_grid_ceil = 20
          end if

        end if

        end function lebedev_grid_ceil
!******

!------------------------------------------------------------------------------
!****f* grids/grid_ceil
!  NAME
!    grid_ceil
!  SYNOPSIS
      integer function grid_ceil( n_ang )
!  PURPOSE
!    function grid_ceil obtains the next higher number
!    of the Lebedev grid for a given number
!    of angular integration grid points.
!
!  USES
        implicit none
! ARGUMENTS
        
        integer :: n_ang

!  INPUTS
!    o n_ang - number of integration points
!  OUTPUTS
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        integer i_grid

        if (n_ang.lt.0) then
           grid_ceil = 0

        else
           do i_grid = 1, n_grids, 1
              grid_ceil = i_grid
              if (n_angular_lebedev(i_grid).ge.n_ang) then
                 exit
              end if
           enddo
        end if

        end function grid_ceil
!******

!------------------------------------------------------------------------------
!****f* grids/grid_floor
!  NAME
!    grid_floor
!  SYNOPSIS
      integer function grid_floor( n_ang )
!  PURPOSE
!    function grid_floor obtains the next lower number
!    of the Lebedev grid for a given number
!    of angular integration grid points.
!  USES
        implicit none
! ARGUMENTS

        integer :: n_ang

!  INPUTS
!    o n_ang - number of integration points
!  OUTPUTS
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        integer :: i_grid

        if (n_ang.lt.n_angular_lebedev(1)) then
           grid_floor = 0

        else
           do i_grid = n_grids, 1, -1
              grid_floor = i_grid
              if (n_angular_lebedev(i_grid).le.n_ang) then
                 exit
              end if
           enddo
        end if

        end function grid_floor
!******

!------------------------------------------------------------------------------
!****f* grids/lebedev_grid_floor
!  NAME
!    lebedev_grid_floor
!  SYNOPSIS
        integer function lebedev_grid_floor ( n_ang )
!  PURPOSE
!    function lebedev_grid_ceil obtains the next lower number
!    of the Lebedev grid for a given number
!    of angular integration grid points.
!  USES

        use runtime_choices, only : force_lebedev
        use mpi_tasks,       only : myid
        use localorb_io,     only : use_unit
!  ARGUMENTS

        integer :: n_ang

!  INPUTS
!    o n_ang - number of integration points
!
!  OUTPUTS
!    none
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


        if (force_lebedev.eq.2) then
            if (n_ang.lt.50) then
                lebedev_grid_floor = 0
            else if (n_ang.lt.62) then
                lebedev_grid_floor = 1
            else if (n_ang.lt.72) then
                lebedev_grid_floor = 2
            else if (n_ang.lt.86) then
                lebedev_grid_floor = 3
            else if (n_ang.lt.98) then
                lebedev_grid_floor = 4
            else if (n_ang.lt.114) then
                lebedev_grid_floor = 5
            else if (n_ang.lt.128) then
                lebedev_grid_floor = 6
            else if (n_ang.lt.146) then
                lebedev_grid_floor = 7
            else if (n_ang.lt.163) then
                lebedev_grid_floor = 8
            else if (n_ang.lt.182) then
                lebedev_grid_floor = 9
            else if (n_ang.lt.201) then
                lebedev_grid_floor = 10
            else if (n_ang.lt.222) then
                lebedev_grid_floor = 11
            else if (n_ang.lt.243) then
                lebedev_grid_floor = 12
            else if (n_ang.lt.266) then
                lebedev_grid_floor = 13
            else if (n_ang.lt.289) then
                lebedev_grid_floor = 14
            else if (n_ang.lt.314) then
                lebedev_grid_floor = 15
            else if (n_ang.lt.339) then
                lebedev_grid_floor = 16
            else if (n_ang.lt.366) then
                lebedev_grid_floor = 17
            else if (n_ang.lt.393) then
                lebedev_grid_floor = 18
            else if (n_ang.lt.422) then
                lebedev_grid_floor = 19
            else if (n_ang.lt.451) then
                lebedev_grid_floor = 20
            else if (n_ang.lt.482) then
                lebedev_grid_floor = 21
            else if (n_ang.lt.513) then
                lebedev_grid_floor = 22
            else if (n_ang.lt.546) then
                lebedev_grid_floor = 23
            else if (n_ang.lt.579) then
                lebedev_grid_floor = 24
            else if (n_ang.lt.614) then
                lebedev_grid_floor = 25
            else if (n_ang.lt.649) then
                lebedev_grid_floor = 26
            else if (n_ang.lt.686) then
                lebedev_grid_floor = 27
            else if (n_ang.lt.723) then
                lebedev_grid_floor = 28
            else if (n_ang.lt.762) then
                lebedev_grid_floor = 29
            else if (n_ang.lt.801) then
                lebedev_grid_floor = 30
            else if (n_ang.lt.842) then
                lebedev_grid_floor = 31
            else if (n_ang.lt.883) then
                lebedev_grid_floor = 32
            else if (n_ang.lt.926) then
                lebedev_grid_floor = 33
            else if (n_ang.lt.969) then
                lebedev_grid_floor = 34
            else if (n_ang.lt.1014) then
                lebedev_grid_floor = 35
            else if (n_ang.lt.1059) then
                lebedev_grid_floor = 36
            else if (n_ang.lt.1106) then
                lebedev_grid_floor = 37
            else if (n_ang.lt.1153) then
                lebedev_grid_floor = 38
            else if (n_ang.lt.1202) then
                lebedev_grid_floor = 39
            else if (n_ang.lt.1251) then
                lebedev_grid_floor = 40
            else if (n_ang.lt.1302) then
                lebedev_grid_floor = 41
            else if (n_ang.lt.1862) then
                lebedev_grid_floor = 42
            else if (n_ang.lt.3282) then
                lebedev_grid_floor = 43
            else if (n_ang.lt.5102) then
                lebedev_grid_floor = 44
            else if (n_ang.lt.7322) then
                lebedev_grid_floor = 45
            else if (n_ang.lt.9942) then
                lebedev_grid_floor = 46
            else if (n_ang.lt.12962) then
                lebedev_grid_floor = 47
            else if (n_ang.lt.16382) then
                lebedev_grid_floor = 48
            else
                lebedev_grid_floor = 49
            end if
        else if ((force_lebedev.eq.1).or.(n_ang.ge.1454)) then
          if (n_ang.lt.6) then
            lebedev_grid_floor = 0
          else if (n_ang.lt.14) then
            lebedev_grid_floor = 1
          else if (n_ang.lt.26) then
            lebedev_grid_floor = 2
          else if (n_ang.lt.38) then
            lebedev_grid_floor = 3
          else if (n_ang.lt.50) then
            lebedev_grid_floor = 4
          else if (n_ang.lt.78) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_floor = 5
          else if (n_ang.lt.86) then
            lebedev_grid_floor = 6
          else if (n_ang.lt.110) then
            lebedev_grid_floor = 7
          else if (n_ang.lt.146) then
            lebedev_grid_floor = 8
          else if (n_ang.lt.170) then
            lebedev_grid_floor = 9
          else if (n_ang.lt.194) then
            lebedev_grid_floor = 10
          else if (n_ang.lt.248) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_floor = 11
          else if (n_ang.lt.288) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_floor = 12
          else if (n_ang.lt.302) then
            lebedev_grid_floor = 13
          else if (n_ang.lt.350) then
            lebedev_grid_floor = 14
          else if (n_ang.lt.434) then
            lebedev_grid_floor = 15
          else if (n_ang.lt.590) then
            lebedev_grid_floor = 16
          else if (n_ang.lt.770) then
            lebedev_grid_floor = 17
          else if (n_ang.lt.974) then
            lebedev_grid_floor = 18
          else if (n_ang.lt.1202) then
            lebedev_grid_floor = 19
          else if (n_ang.lt.1454) then
            lebedev_grid_floor = 20
          else if (n_ang.lt.1730) then
            lebedev_grid_floor = 21
          else if (n_ang.lt.2030) then
            lebedev_grid_floor = 22
          else if (n_ang.lt.2354) then
            lebedev_grid_floor = 23
          else if (n_ang.lt.2702) then
            lebedev_grid_floor = 24
          else if (n_ang.lt.3074) then
            lebedev_grid_floor = 25
          else if (n_ang.lt.3470) then
            lebedev_grid_floor = 26
          else if (n_ang.lt.3890) then
            lebedev_grid_floor = 27
          else if (n_ang.lt.4334) then
            lebedev_grid_floor = 28
          else if (n_ang.lt.4802) then
            lebedev_grid_floor = 29
          else if (n_ang.lt.5294) then
            lebedev_grid_floor = 30
          else if (n_ang.lt.5810) then
            lebedev_grid_floor = 31
          else if (n_ang.eq.5810) then
            lebedev_grid_floor = 32
          else
            if (n_ang.gt.5810) then
              if (myid.eq.0) then
                write(use_unit,'(1X,A,A)') &
                "! Defaulting to maximum number of angular grid points: ", &
                "5810."
              end if
              lebedev_grid_floor = 32
            end if
         end if
      else
!         Delley style grids; not all are supplied
          if (n_ang.lt.6) then
            lebedev_grid_floor = 0
          else if (n_ang.lt.14) then
            lebedev_grid_floor = 1
          else if (n_ang.lt.26) then
            lebedev_grid_floor = 2
          else if (n_ang.lt.50) then
            lebedev_grid_floor = 3
          else if (n_ang.lt.110) then
            lebedev_grid_floor = 5
          else if (n_ang.lt.194) then
            lebedev_grid_floor = 8
          else if (n_ang.lt.302) then
            lebedev_grid_floor = 11
          else if (n_ang.lt.434) then
            lebedev_grid_floor = 14
          else if (n_ang.lt.590) then
            lebedev_grid_floor = 16
          else if (n_ang.lt.770) then
            lebedev_grid_floor = 17
          else if (n_ang.lt.974) then
            lebedev_grid_floor = 18
          else if (n_ang.lt.1202) then
            lebedev_grid_floor = 19
          else
            lebedev_grid_floor = 20
          end if

        end if

        end function lebedev_grid_floor
!******
!------------------------------------------------------------------------------------
!****s* grids/get_logarithmic_grid
!  NAME
!    get_logarithmic_grid
!  SYNOPSIS
      subroutine get_logarithmic_grid ( i_species )
! PURPOSE
!    Subroutine get_logarithmic_grid sets up the logarithmic integration grids
!    for radial atomic wave functions (basis functions). Input data are specified
!    in control.in .
!    grid definition: 
!       r(i) = r_grid_min/z * exp[ log(r_grid_inc) * (i-1) ]
!  USES
         use mpi_tasks,   only : myid
         use localorb_io, only : use_unit
         implicit none
!  ARGUMENTS

         integer :: i_species

!  INPUTS
!    o i_species -- index of the species in question
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


        real*8 :: r_test
        integer :: i_grid


        i_grid = 0
        r_test = r_grid_min(i_species)

        do while ( r_test.le.r_grid_max(i_species) )

          i_grid = i_grid + 1

          r_grid(i_grid, i_species) = r_test
          r_test=r_test*r_grid_inc(i_species)

        end do

        if (i_grid.le.0) then
           if (myid.eq.0) then
              write(use_unit,*) "Illegal logarithmic grid for species ", &
                   i_species, ". r_min > r_max. Aborting."
           end if
           stop
        else
          n_grid(i_species) = i_grid
        end if

        end subroutine get_logarithmic_grid
!******
!-----------------------------------------------------------------------------------
!****s* grids/test_radial_grid
!  NAME
!    test_radial_grid
!  SYNOPSIS

      subroutine test_radial_grid (  )

!  PURPOSE
!    Debug subroutine, not normally used.
!    subroutine test_radial_grid performs two simple 1D test integrals with functions
!    whose characteristic width is 1 bohr: a step function, and an exponential function
!  USES
      use dimensions,   only : n_max_radial, n_species
      use mpi_tasks,    only : myid
      use localorb_io,  only : use_unit, localorb_info
      implicit none

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


!  local variables

      real*8 sum
      real*8 Radius
      real*8 funct(n_max_radial)

!  counters

      integer i_species, i_radial

!  begin tests


      call localorb_info('',use_unit)
      call localorb_info("* Testing radial integration grids:",use_unit)

      do i_species = 1, n_species, 1
        do i_radial = 1, n_radial(i_species), 1

           if (myid.eq.0) then
              open (50, file = "grid.dat")
              write (50,*) r_radial (i_radial, i_species), "1.0", &
                   w_radial (i_radial, i_species)
              close (50)
           end if

        enddo
      enddo

      do i_species = 1, n_species, 1

         if (myid.eq.0) then
            write(use_unit,*) "* Species ", i_species, ":"
         end if

! step integral
        Radius = 1.
        do i_radial = 1, n_radial(i_species), 1
          if (r_radial (i_radial,i_species).le.Radius ) then
            funct(i_radial) = 1.d0
        else
          funct(i_radial) = 0.d0
        end if
        enddo
        sum = 0.
        do i_radial = 1, n_radial(i_species), 1
          sum = sum + funct(i_radial) * w_radial(i_radial,i_species)

        enddo

        if (myid.eq.0) then
           write(use_unit,*) "* Step Integral: ", sum, ", should be 1."
        end if

! exponential integral
        do i_radial = 1, n_radial(i_species), 1
           funct(i_radial) = exp (- r_radial(i_radial,i_species) )
        enddo
        sum = 0.
        do i_radial = 1, n_radial(i_species), 1
          sum = sum + funct(i_radial) * w_radial(i_radial,i_species)
        enddo

        if (myid.eq.0) then
           write(use_unit,*) "* Exponential Integral: ", sum, ", should be 1."
        end if

      enddo

      write(use_unit,*)

      end subroutine test_radial_grid
!******

!****s* grids/cleanup_grids
!  NAME
!     cleanup_grids
!  SYNOPSIS
        subroutine cleanup_grids( )
!  PURPOSE
!     deallocation of all arrays pertaining to the grids. 
! USES
          use dimensions, only: use_vdw_method
          implicit none

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


        if ( allocated (r_grid_min) ) then
          deallocate ( r_grid_min )
        end if
        if ( allocated (r_grid_max) ) then
          deallocate ( r_grid_max )
        end if
        if ( allocated (r_grid_inc) ) then
          deallocate ( r_grid_inc )
        end if
        if ( allocated (log_r_grid_inc) ) then
          deallocate ( log_r_grid_inc )
        end if

        if ( allocated (n_grid ) ) then
          deallocate ( n_grid )
        end if
        if ( allocated (r_grid ) ) then
          deallocate ( r_grid )
        end if

        if ( allocated (n_radial ) ) then
          deallocate ( n_radial )
        end if
        if ( allocated (scale_radial ) ) then
          deallocate ( scale_radial )
        end if
        if ( allocated (n_angular ) ) then
          deallocate ( n_angular )
        end if
        if ( allocated (angular_limit ) ) then
          deallocate ( angular_limit )
        end if
        if ( allocated (angular_acc ) ) then
          deallocate ( angular_acc )
        end if
        if ( allocated (angular_min ) ) then
          deallocate ( angular_min )
        end if

        if ( allocated (r_radial ) ) then
          deallocate ( r_radial )
        end if
        if ( allocated (w_radial ) ) then
          deallocate ( w_radial )
        end if
        if ( allocated (r_angular ) ) then
          deallocate ( r_angular )
        end if
        if ( allocated (w_angular ) ) then
          deallocate ( w_angular )
        end if
        if ( allocated (n_angular ) ) then
           deallocate ( n_angular )
        end if
        if (allocated(local_ylm_tab)) then
          deallocate ( local_ylm_tab )
        end if
        if (allocated (n_division ) ) then
           deallocate (n_division )
        end if
        if (allocated (division_boundaries)) then
           deallocate (division_boundaries)
        end if

        if (allocated(n_angular_lebedev)) then
           deallocate(n_angular_lebedev)
        end if
        if (allocated(r_angular_lebedev)) then
           deallocate(r_angular_lebedev)
        end if
        if (allocated(w_angular_lebedev)) then
           deallocate(w_angular_lebedev)
        end if

        
        if (allocated(n_ang_lebedev)) then  !SAG
           deallocate(n_ang_lebedev)
        end if
        if (allocated(r_ang_lebedev)) then
           deallocate(r_ang_lebedev)
        end if
        if (allocated(w_ang_lebedev)) then
           deallocate(w_ang_lebedev)
        end if

        
        if (allocated(n_division_lebedev)) then
           deallocate(n_division_lebedev)
        end if
        if (allocated(division_boundaries_lebedev)) then
           deallocate(division_boundaries_lebedev)
        end if

        if (allocated(fixed_grid_index)) then
           deallocate(fixed_grid_index)
        end if

        if (allocated(lebedev_grid_index)) then
           deallocate(lebedev_grid_index)
        end if

        end subroutine cleanup_grids
!******

!****f* grids/invert_radial_grid
!  NAME
!    invert_radial_grid
!  SYNOPSIS

        real*8 function invert_radial_grid(r_current, n_scale, r_scale)

!  PURPOSE
!    Function invert_radial_grid takes a given real-space radius value
!    and returns a (fractional) number i(r) which indicates its position
!    between the grid points i_radial of the radial integration grid.
!
!    Grid chosen according to Eq. (2) of B. Delley review 1995.
!
!    Since 
!        r(i) = - C ln ( 1 - (i/(n+1))**2 ), we have
!
!        i(r) = (n+1) * ( 1 - exp(-r/C) )**(1/2)
!
!  USES
        implicit none

!  ARGUMENTS

        real*8  :: r_current
        integer :: n_scale
        real*8  :: r_scale

!  INPUTS
!    o r_current -- Radius for which i(r) is desired
!    o n_scale -- Total number of grid points on integration grid
!    o r_scale -- Radial scaling parameter of integration grid
!
! OUTPUTS
!    o invert_radial_grid -- radial grid index
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        invert_radial_grid = 1.0d0 - exp( - r_current/r_scale)
        invert_radial_grid = &
          dble(n_scale + 1) * sqrt(invert_radial_grid)
      end function invert_radial_grid
!******
!----------------------------------------------------------------------------------
!****f* grids/invert_log_grid
!  NAME
!    invert_log_grid
!  SYNOPSIS

     real*8 function invert_log_grid ( r_current,r_min,scale)

!  PURPOSE
!    Function invert_log_grid takes a given real-space radius value
!    and returns a (fractional) number i(r) which indicates its position
!    between the grid points i_grid of the logarithmic grid.
!
!    i(r) = 1 + [ln(r / r_grid_min)]/[ln(r_grid_inc)]
!
!  USES

        implicit none
!  ARGUMENTS

        real*8 :: r_current
        real*8 :: r_min
        real*8 :: scale

!  INPUTS
!    o r_current -- radius to be inverted
!    o r_min -- starting radius for log grid
!    o scale -- logarithmic scaling parameter
!
! OUTPUTS
!    o invert_log_grid -- logarithmic grid index
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



        invert_log_grid = log (r_current/r_min)
        invert_log_grid = 1.d0 + invert_log_grid/log(scale)
      end function invert_log_grid
!******
!----------------------------------------------------------------------------
!****s* grids/invert_log_grid_deriv
!  NAME
!    invert_log_grid_deriv
!  SYNOPSIS

!real*8 function invert_log_grid_deriv(r_current, r_min, scale)
     real*8 function invert_log_grid_deriv(r_current, scale)

!  PURPOSE
!    Derivative i'(r) of invert_log_grid.
!  USES

        implicit none

!  ARGUMENTS

        real*8, intent(IN) :: r_current
!       real*8, intent(IN) :: r_min
        real*8, intent(IN) :: scale

!  INPUTS
!    o r_current -- radius to be inverted
!    o r_min -- starting radius for log grid
!    o scale -- logarithmic scaling parameter
!
!  OUTPUTS
!    o invert_log_grid_deriv -- derivative of logarithmic grid index
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2011).
!  SOURCE

        character(*), parameter :: func = 'invert_log_grid_deriv'

!       invert_log_grid = 1.d0 + log(r_current/r_min) / log(scale)
        invert_log_grid_deriv = 1.d0 / r_current / log(scale)

      end function invert_log_grid_deriv
!******
!------------------------------------------------------------------------------
!****f* grids/invert_log_grid_p2
!  NAME
!    invert_log_grid_p2
!  SYNOPSIS
        real*8 function invert_log_grid_p2( r_current, i_species)
! PURPOSE
!   Function invert_log_grid_p2 takes a given real-space radius value
!   and returns a (fractional) number i(r) which indicates its position
!   between the grid points i_grid of the logarithmic grid.
!
!   i(r) = 1 + [ln(r / r_grid_min)]/[ln(r_grid_inc)]
!
!   This does not care about scales, it only takes a species number and looks up its own parameters.
!
!  ARGUMENTS

        real*8  :: r_current
        integer :: i_species

!  INPUTS
!    o r_current -- radius to be inverted
!    o i_species -- number of the species
!
! OUTPUTS
!    o invert_log_grid_p2 -- logarithmic grid index
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        invert_log_grid_p2 = log (r_current/r_grid_min(i_species))
        invert_log_grid_p2 = 1.d0 + invert_log_grid_p2/log_r_grid_inc(i_species)
        end function invert_log_grid_p2
!******
!-------------------------------------------------------------------------------
!****f* grids/get_radial_weight
!  NAME
!    get_radial_weight
!  SYNOPSIS

      real*8 function get_radial_weight(index, scale, n_max)

!  PURPOSE
!     gives di/dr !!
!
!  ARGUMENTS

      real*8  :: index
      real*8  :: scale
      integer :: n_max

!  INPUTS
!    o index -- shell index
!    o scale -- shell scale
!    o n_max -- number of shells

! OUTPUTS
!    o get_radial_weight -- di/dr for given radial shell
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      get_radial_weight = (dble(n_max + 1) / index) ** 2.d0
      get_radial_weight = index / (2.d0 * scale) * &
           (get_radial_weight - 1)
      end function get_radial_weight
!******
!-------------------------------------------------------------------------------
!****s* grids/get_local_ylm_tab
!  NAME
!    get_local_ylm_tab
!  SYNOPSIS
        subroutine get_local_ylm_tab ( l_hartree )
!  PURPOSE
!       tabulate Ylm function at each integration point around each species
!       for the benefit of update_hartree_potential later on
!  USES

        use dimensions, only : n_species, l_pot_max, n_max_angular
        implicit none

!  ARGUMENTS

       integer, dimension(n_species) :: l_hartree

!  INPUTS
!    o l_hartree -- max angular momentum shell for integrations
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


        real*8, dimension(3) :: dir_tab
        real*8, dimension(4) :: trigonom_tab
        integer :: i_species
        integer :: i_radial
        integer :: i_angular
        integer :: i_lebedev
        integer :: max_lebedev_index

        max_lebedev_index = -1
        do i_species = 1, n_species, 1
           do i_radial = 1, n_radial( i_species ), 1
              max_lebedev_index = &
                   MAX( lebedev_grid_index( i_radial, i_species ), &
                   max_lebedev_index )
           end do
        end do
        if (.not.allocated(local_ylm_tab)) then
           allocate ( local_ylm_tab ( (l_pot_max+1)**2, n_max_angular, &
                max_lebedev_index ) )
           local_ylm_tab = 0.0d0
        end if

        do i_species = 1, n_species, 1
           do i_radial = 1, n_radial( i_species ), 1
              i_lebedev = lebedev_grid_index( i_radial, i_species )
              do i_angular = 1, n_angular( i_radial, i_species ), 1
!                dist_tab = r_radial( i_radial, i_species )
                 dir_tab(:) = &
                      r_angular_lebedev( : , i_angular, i_lebedev )

!               compute trigonometric functions of spherical coordinate angles
!               of current integration point
                 call tab_local_trigonom &
                ( dir_tab, trigonom_tab &
                )

                call tab_local_ylm &
                ( trigonom_tab, l_hartree(i_species), l_pot_max, &
                  local_ylm_tab(:,i_angular,i_lebedev) )

              enddo
           enddo
        enddo

      end subroutine get_local_ylm_tab
!******
!-------------------------------------------------------------------------
!****s* grids/invert_log_grid_vector
!  NAME
!    invert_log_grid_vector
!  SYNOPSIS
        subroutine invert_log_grid_vector ( r_current,r_min,scale,n_points,out_value )
!  PURPOSE
!    Subroutine invert_log_grid_vector takes a set of real-space radiuses
!    and returns (fractional) numbers i(r) which indicates its position
!    between the grid points i_grid of the logarithmic grid.
!
!    i(r) = 1 + [ln(r / r_grid_min)]/[ln(r_grid_inc)]
!
!  USES
    implicit none
!  ARGUMENTS

        integer :: n_points
        real*8  :: r_current(n_points)
        real*8  :: r_min
        real*8  :: scale
        real*8 :: out_value(n_points)

!  INPUTS
!    o n_points -- number of points
!    o r_current -- radii to be inverted
!    o r_min -- log grid minimal radius
!    o scale -- log grid scale
!
!  OUTPUTS
!    o out_value -- inverted radii
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        out_value = log (r_current/r_min)
        out_value = 1.d0 + out_value/log(scale)
        end subroutine invert_log_grid_vector
!******
      end module grids
