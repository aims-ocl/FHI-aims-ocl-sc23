
!****h* FHI-aims/grids
!  NAME
!    grids
!  SYNOPSIS
      module esp_grids
! PURPOSE
!  Module esp_grids handles everything related to the real-space grids for 
!  the fitting of esp-charges (modified grids module):
!  * the logarithmic radial grids for free-atom related quantities
!  * the radial integration grids
!  * the Lebedev-style angular grids
!
!  Subroutines:
!  * allocate_grids_esp
!  * get_grids_esp
!  * get_lebedev_esp - historic and no longer used
!  * test_radial_grid_esp
!  * get_logarithmic_grid_esp
!  * cleanup_grids_esp
!
!  Functions:
!  * invert_radial_grid_esp
!  * invert_log_grid_esp
!  * get_radial_weight_esp
!
!  global variable declarations - exported to other program parts:

!  *   r_grid_min_esp: for each species, the innermost and
!  *   r_grid_max_esp : and outermost points of the logarithmic grid for wave
!  *   r_grid_inc: functions, and the increase factor
!                 !!! Notice that r_grid_min_esp = r_grid_min_esp/species_z !!!
!                 on (Lebedev) angular integration grids
!  *   n_grid_esp    : number of radial grid points for atomic(!) potential, density,
!                 and radial wave functions for each species
!                 Please note that
!                    r_grid_max_esp (i) ~ r_grid_min_esp * r_grid_inc**(n_grid_esp-1)
!                 is only approximately true.
!  *   r_grid_esp    : Radial grid points for atomic(!) potential, density,
!                 radial wave functions - for each species

!  *   n_radial_esp : for each species, number of points for radial integration grids
!  *   scale_radial_esp : for each species, scale factor for radial integration grids
!  *   angular_acc_esp : requested integration accuracy per integration shell
!  *   angular_limit_esp : maximum number of angular integration points per shell for a given species
!  *   r_radial_esp  : for each species, radial integration grid points
!  *   w_radial_esp  : for each species, radial integration grid weights
   
!  *   n_angular_esp : for each atom and each radial shell(!), actual number of grid points (Angular grids must be allocated per atom.)
!  *   r_angular_esp : for each atom, angular integration grid points
!  *   w_angular_esp : for each atom, angular integration grid weights
!
!  USES
      use constants,only:bohr
      use localorb_io,only:use_unit
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

      real*8,  dimension(:),       allocatable :: r_grid_min_esp            
      ! r_grid_min(n_species)
      real*8,  dimension(:),       allocatable :: r_grid_max_esp                        
      ! r_grid_max(n_species)
      real*8,  dimension(:),       allocatable :: r_grid_inc_esp                        
      ! r_grid_inc(n_species)
      real*8,  dimension(:),       allocatable :: log_r_grid_inc_esp                    
      ! r_grid_inc(n_species)
      integer, dimension(:),       allocatable :: n_grid_esp                            
      ! n_grid_esp(n_species)
      real*8,  dimension(:,:),     allocatable :: r_grid_esp                            
      ! r_grid_esp(n_max_grid, n_species)
      integer, dimension(:),       allocatable :: n_radial_esp                        
      ! n_radial_esp(n_species)
      real*8,  dimension(:),       allocatable :: scale_radial_esp                      
      ! scale_radial_esp(n_species)
      integer, dimension(:),       allocatable :: angular_limit_esp                     
      ! angular_limit_esp(n_species)
      integer, dimension(:),       allocatable :: angular_min_esp                       
      ! angular_min_esp(n_species)
      real*8,  dimension(:),       allocatable :: angular_acc_esp                       
      ! angular_acc_esp(n_species)
      real*8,  dimension(:,:),     allocatable :: r_radial_esp                          
      ! r_radial_esp(n_max_radial, n_species)
      real*8,  dimension(:,:),     allocatable :: w_radial_esp                          
      ! w_radial_esp(n_max_radial, n_species )
      integer, dimension(:,:),     allocatable :: n_angular_esp                         
      ! n_angular_esp(n_max_radial, n_species)
      real*8,  dimension(:,:,:,:), allocatable :: r_angular_esp                         
      ! r_angular_esp(3, n_max_angular, n_max_radial, n_species)
      real*8,  dimension(:,:,:),   allocatable :: w_angular_esp                         
      ! w_angular_esp(n_max_angular, n_max_radial, n_species)
      integer,    dimension(:),     allocatable :: n_ang_lebedev_esp                     
      !for nlcorr. SAG
      real*8,  dimension(:,:,:),   allocatable :: r_ang_lebedev_esp
      real*8,  dimension(:,:),     allocatable :: w_ang_lebedev_esp
          
      integer, dimension(:),       allocatable :: n_angular_lebedev_esp
      real*8,  dimension(:,:,:),   allocatable :: r_angular_lebedev_esp
      real*8,  dimension(:,:),     allocatable :: w_angular_lebedev_esp
      integer, dimension(:),       allocatable :: n_division_lebedev_esp
      integer, dimension(:,:),     allocatable :: &
                                        division_boundaries_lebedev_esp
      integer, dimension(:,:),     allocatable :: fixed_grid_index_esp
      integer, dimension(:,:),     allocatable :: lebedev_grid_index_esp
      integer, dimension(:,:),     allocatable :: n_division_esp                       
      ! n_division_esp(n_max_radial, n_species)
      integer, dimension(:,:,:),   allocatable :: division_boundaries_esp              
      ! division_boundaries_esp(n_max_angular_division+1,
      real*8,  dimension(:,:,:),   allocatable :: local_ylm_tab_esp

      integer       :: n_max_lebedev_esp
      integer       :: n_grids_esp
      logical       :: grid_partitioned_esp = .false.

      type grid_point_esp
        real*8, dimension(3) :: coords_esp
        integer :: index_atom_esp
        integer :: index_radial_esp
        integer :: index_angular_esp
      end type grid_point_esp

      type batch_of_points_esp
        integer :: size_esp
!        type(grid_point), dimension(:), allocatable :: points
        type(grid_point_esp), pointer, dimension(:) :: points_esp
        integer :: batch_n_compute_esp
        integer, pointer, dimension(:) :: batch_i_basis_esp
      end type batch_of_points_esp

      type(batch_of_points_esp), pointer, dimension(:) :: batches_esp

      integer :: n_grid_batches_esp
      integer :: n_points_in_batch_esp
      integer :: n_my_batches_esp
      integer :: n_max_batch_size_esp

      integer :: n_full_points_esp

      integer :: n_max_compute_atoms_esp
      integer :: n_max_compute_dens_esp
      integer :: n_max_compute_fns_dens_esp
      logical :: got_n_compute_maxes_esp = .false.
!******

      contains
!---------------------------------------------------------------------------
!****s* grids/allocate_grids
!  NAME
!    allocate_grids
!  SYNOPSIS
        subroutine allocate_grids_esp( )
!  PURPOSE
!  Subroutine allocate_grids allocates the necessary memory for all real-space 
!  grids. This allocation is done once at the beginning of the code, since the 
!  grids are needed everywhere.
!
!  In principle one could allocate everything in a much more fine-grained way - 
!  but only if there is a clear need.
!  USES
        use dimensions,only:n_species,n_max_radial,n_max_angular,&
                            n_max_angular_division,n_max_grid

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

        n_max_lebedev_esp = 50
        allocate ( r_grid_min_esp(n_species) )
        allocate ( r_grid_max_esp(n_species) )
        allocate ( r_grid_inc_esp(n_species) )
        allocate ( log_r_grid_inc_esp(n_species) )
        allocate ( n_grid_esp(n_species) )
        allocate ( r_grid_esp(n_max_grid, n_species) )
        allocate ( n_radial_esp(n_species) )
        allocate ( scale_radial_esp(n_species) )
        allocate ( angular_limit_esp(n_species) )
        allocate ( angular_min_esp(n_species) )
        allocate ( angular_acc_esp(n_species) )
        allocate ( r_radial_esp(n_max_radial, n_species) )
        allocate ( w_radial_esp(n_max_radial, n_species) )
        allocate ( n_angular_esp(n_max_radial, n_species) )
        allocate ( r_angular_esp(3, n_max_angular, n_max_radial, n_species))
        allocate ( w_angular_esp(n_max_angular, n_max_radial, n_species) )
        allocate ( n_division_esp(n_max_radial, n_species) )
        allocate ( division_boundaries_esp(n_max_angular_division+1, &
                   n_max_radial, n_species) )
        allocate(n_angular_lebedev_esp(n_max_lebedev_esp))
        allocate(r_angular_lebedev_esp(3,n_max_angular,n_max_lebedev_esp))
        allocate(w_angular_lebedev_esp(n_max_angular,n_max_lebedev_esp))
        allocate(n_division_lebedev_esp(n_max_lebedev_esp))
        allocate ( division_boundaries_lebedev_esp(n_max_angular_division+1,&
                   n_max_lebedev_esp) )
        allocate(lebedev_grid_index_esp(n_max_radial, n_species))
        
        
        allocate(n_ang_lebedev_esp(15))  !SAG
        allocate(r_ang_lebedev_esp(3,350,15))
        allocate(w_ang_lebedev_esp(350,15))
        


        end subroutine allocate_grids_esp
!******
!------------------------------------------------------------------------
!****s* esp_grids/get_grids_esp
!  NAME
!    get_grids_esp
!  SYNOPSIS
        subroutine get_grids_esp(radius_esp_min_new, radius_esp_max_new, &
                   out_grids,esp_log_grid_in )
!  PURPOSE
!    set up the grids for various species
!  USES
        use dimensions,only:n_max_angular,n_species,use_nlcorr_post,&
                            n_max_points_per_div,use_vdw_post,use_vdw_method
        use mpi_tasks,only:myid
        use localorb_io,only:localorb_info
        implicit none

!  ARGUMENTS
        real*8,  dimension(n_species) :: radius_esp_min_new
        real*8,  dimension(n_species) :: radius_esp_max_new
        logical :: out_grids
        logical :: esp_log_grid_in
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
	r_grid_max_esp(:) = 100.
	r_grid_inc_esp(:) = 1.0123
	angular_limit_esp(:) = n_max_angular
        n_points_in_batch_esp = 100
!     allocate and initialize array for fixed grid data
!        if (any(species_fixed_grid)) then
!           allocate(fixed_grid_index_esp(n_max_radial, n_species))
!           fixed_grid_index_esp = 0
!        end if

!       initialize logarithmic grid size
        r_grid_min_esp  = minval(radius_esp_min_new)
        do i_species = 1, n_species, 1
          if (esp_log_grid_in)then
            scale_radial_esp(i_species) =  -(radius_esp_max_new(i_species)-&
                                          radius_esp_min_new(i_species)) / &
            log( 1.d0 - &
                (n_radial_esp(i_species)/(1.d0+n_radial_esp(i_species)))**2.d0 )
          else
            scale_radial_esp(i_species) =  (radius_esp_max_new(i_species)-&
                                         radius_esp_min_new(i_species)) / &
            (n_radial_esp(i_species)/(1.d0+n_radial_esp(i_species)))
          endif
        enddo
        do i_species = 1, n_species, 1
          n_grid_esp(i_species) = 0
        enddo

!       set up grids species by species

        do i_species = 1, n_species,1

! ----->  set up the radial grid

          do i_radial = 1, n_radial_esp (i_species)

!           scale radial grid by n
            r_scaled = dble(i_radial/(n_radial_esp(i_species)+1.0d0))

!           get radial grid point
            if (esp_log_grid_in)then
              r_radial_esp (i_radial, i_species ) = -log(1.d0-r_scaled**2.d0)
            else
              r_radial_esp (i_radial, i_species ) = r_scaled
            endif

            r_radial_esp (i_radial, i_species) = radius_esp_min_new(i_species) +&
              scale_radial_esp(i_species)*r_radial_esp(i_radial, i_species)

!           get radial grid weight dr/di
            w_radial_esp (i_radial, i_species) = &
            (2.0d0/dble(i_radial)) / ( 1.0d0/(r_scaled**2.d0) - 1.0d0 )

            w_radial_esp (i_radial, i_species) = &
            scale_radial_esp(i_species)*w_radial_esp(i_radial, i_species)
          enddo

!  Check if minimal logarithmic grid point is smaller than minimal
!  radial integration grid point; if not, adjust.

          if (r_grid_min_esp(i_species).gt.r_radial_esp(1,i_species)) then
             if (myid.eq.0) then
                write(use_unit,*)
                write(use_unit,*) "* Minimum radial grid point for ", &
                     "logarithmic grid of species ", i_species
                write(use_unit,*) "* , r_min = ", &
                     r_grid_min_esp(i_species), ","
                write(use_unit,*) "* is chosen above the minimum radial", &
                     " integration grid point,"
                write(use_unit,*) "* r_min = ", r_radial_esp(i_species,1), "."
             end if

             r_grid_min_esp(i_species) = 0.5d0*r_radial_esp(1,i_species)

             if (myid.eq.0) then
                write(use_unit,*) "* Setting r_grid_min_esp to ", &
                     r_grid_min_esp(i_species), "."
                write(use_unit,*)
             end if

          end if

          call get_logarithmic_grid_esp &
          ( i_species )

          log_r_grid_inc_esp(i_species) = log(r_grid_inc_esp(i_species))

        enddo

!       Initialize angular grids separately.
!       Angular grids are also self-adapting, and reset in initialize_integrals.

        do i_species = 1, n_species, 1

           do i_radial = 1, n_radial_esp(i_species)

              n_lebedev = &
                   lebedev_grid_floor_esp(angular_limit_esp(i_species))
              call get_angular_grid &
                   ( n_lebedev, &
                   n_angular_esp(i_radial,i_species), &
                   r_angular_esp(1,1,i_radial,i_species), &
                   w_angular_esp(1,i_radial,i_species) &
                   )
              if (n_angular_esp(i_radial,i_species).gt.n_max_angular) then
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
                   n_angular_esp(i_radial, i_species), &
                   r_angular_esp(1,1,i_radial,i_species), &
                   w_angular_esp(1,i_radial,i_species), &
                   n_division_esp(i_radial, i_species), &
                   division_boundaries_esp(1,i_radial,i_species), &
                   1.0d8 )
           enddo

        enddo

        ! next, tabulate Lebedev grids themselved for the possibility of 
        ! self-adapting grids later on. This option is not the default, but 
        ! since tabulating the Lebedev  grids does not cost anything, we retain 
        ! it here.

        ! initialize unused parts of auxiliary arrays to avoid suprious floating 
        ! exceptions later
        try_r_angular_lebedev = 0.d0
        try_w_angular_lebedev = 0.d0

        n_lebedev = 0

        do i_species = 1, n_species, 1

           n_lebedev = MAX(n_lebedev, &
                lebedev_grid_floor_esp(angular_limit_esp(i_species)))
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

              n_angular_lebedev_esp(i_grid) = try_n_angular_lebedev
              r_angular_lebedev_esp(:,:,i_grid) = try_r_angular_lebedev
              w_angular_lebedev_esp(:,i_grid) = try_w_angular_lebedev

              call divide_angular_grid_p0 &
                   ( n_angular_lebedev_esp(i_grid), &
                   r_angular_lebedev_esp(1,1,i_grid), &
                   w_angular_lebedev_esp(1,i_grid), &
                   n_division_lebedev_esp(i_grid), &
                   division_boundaries_lebedev_esp(1,i_grid), &
                   n_angular_lebedev_esp(i_grid) &
                   )

              i_lebedev_prev = i_lebedev_copy

              do i_division = 1, n_division_lebedev_esp(i_grid), 1
                n_max_points_per_div = &
                  max( n_max_points_per_div, &
                       division_boundaries_lebedev_esp(i_division+1,i_grid) &
                      -division_boundaries_lebedev_esp(i_division,i_grid) )
              enddo

           end if

        enddo

        n_grids_esp = i_grid

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
                do i_radial = 1, n_radial_esp (i_species)
                   write(50,*) r_radial_esp(i_radial,i_species), &
                        w_radial_esp(i_radial,i_species)
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
             ! call get_angular_grid(i_leb, try_n_angular_lebedev, &  !wouldn't
             ! work!
             !try_r_angular_lebedev, try_w_angular_lebedev )
             !n_ang_lebedev_esp(i_leb) = try_n_angular_lebedev
             !r_ang_lebedev_esp(:,:,i_leb) = try_r_angular_lebedev(:,1:350)
             !w_ang_lebedev_esp(:,i_leb) = try_w_angular_lebedev(1:350)
             
              call get_leb_grids(i_leb, nal,ral,wal,350)
             n_ang_lebedev_esp(i_leb) = nal
             r_ang_lebedev_esp(:,:,i_leb) = ral
             w_ang_lebedev_esp(:,i_leb) = wal
             
             
             
          end do
          
       endif





     end subroutine get_grids_esp
!******
!-------------------------------------------------------------------------------
!****s* esp_grids/get_lebedev_esp
!  NAME
!    get_lebedev_esp
!  SYNOPSIS
        subroutine get_lebedev_esp( n_ang, i_species, i_radial )
! PURPOSE
!  obtains the Lebedev integration grid points for given angular resolution.
!
!  This version is a legacy version and no longer used; we now use
!  subroutine get_angular_grids for the same purpose
!
!   USES
        use dimensions,only:n_max_angular
        use mpi_tasks,only:myid
        use localorb_io,only:localorb_info

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
           call localorb_info("* No suitable Lebedev grid found.", &
           use_unit)
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
          r_angular_esp (1, i_ang, i_radial, i_species) = x(i_ang)
          r_angular_esp (2, i_ang, i_radial, i_species) = y(i_ang)
          r_angular_esp (3, i_ang, i_radial, i_species) = z(i_ang)
          w_angular_esp (i_ang, i_radial, i_species)    = w(i_ang)
!test
!          write(use_unit,*) "Next i_ang."
!test end

        enddo

!test
!          write(use_unit,*) "Leave get_lebedev_esp."
!test end

        end subroutine get_lebedev_esp
!******

!------------------------------------------------------------------------------
!****f* esp_grids/lebedev_grid_ceil_esp
!  NAME
!    lebedev_grid_ceil_esp
!  SYNOPSIS
        integer function lebedev_grid_ceil_esp ( n_ang )
! PURPOSE
!  function lebedev_grid_ceil_esp obtains the next higher number
!  of the Lebedev grid for a given number
!  of angular integration grid points.
! USES
        use dimensions,only:n_max_angular
        use mpi_tasks,only:myid
        use runtime_choices,only:force_lebedev 

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
                lebedev_grid_ceil_esp = 0
            else if (n_ang.le.50) then
                lebedev_grid_ceil_esp = 1
            else if (n_ang.le.62) then
                lebedev_grid_ceil_esp = 2
            else if (n_ang.le.72) then
                lebedev_grid_ceil_esp = 3
            else if (n_ang.le.86) then
                lebedev_grid_ceil_esp = 4
            else if (n_ang.le.98) then
                lebedev_grid_ceil_esp = 5
            else if (n_ang.le.114) then
                lebedev_grid_ceil_esp = 6
            else if (n_ang.le.128) then
                lebedev_grid_ceil_esp = 7
            else if (n_ang.le.146) then
                lebedev_grid_ceil_esp = 8
            else if (n_ang.le.163) then
                lebedev_grid_ceil_esp = 9
            else if (n_ang.le.182) then
                lebedev_grid_ceil_esp = 10
            else if (n_ang.le.201) then
                lebedev_grid_ceil_esp = 11
            else if (n_ang.le.222) then
                lebedev_grid_ceil_esp = 12
            else if (n_ang.le.243) then
                lebedev_grid_ceil_esp = 13
            else if (n_ang.le.266) then
                lebedev_grid_ceil_esp = 14
            else if (n_ang.le.289) then
                lebedev_grid_ceil_esp = 15
            else if (n_ang.le.314) then
                lebedev_grid_ceil_esp = 16
            else if (n_ang.le.339) then
                lebedev_grid_ceil_esp = 17
            else if (n_ang.le.366) then
                lebedev_grid_ceil_esp = 18
            else if (n_ang.le.393) then
                lebedev_grid_ceil_esp = 19
            else if (n_ang.le.422) then
                lebedev_grid_ceil_esp = 20
            else if (n_ang.le.451) then
                lebedev_grid_ceil_esp = 21
            else if (n_ang.le.482) then
                lebedev_grid_ceil_esp = 22
            else if (n_ang.le.513) then
                lebedev_grid_ceil_esp = 23
            else if (n_ang.le.546) then
                lebedev_grid_ceil_esp = 24
            else if (n_ang.le.579) then
                lebedev_grid_ceil_esp = 25
            else if (n_ang.le.614) then
                lebedev_grid_ceil_esp = 26
            else if (n_ang.le.649) then
                lebedev_grid_ceil_esp = 27
            else if (n_ang.le.686) then
                lebedev_grid_ceil_esp = 28
            else if (n_ang.le.723) then
                lebedev_grid_ceil_esp = 29
            else if (n_ang.le.762) then
                lebedev_grid_ceil_esp = 30
            else if (n_ang.le.801) then
                lebedev_grid_ceil_esp = 31
            else if (n_ang.le.842) then
                lebedev_grid_ceil_esp = 32
            else if (n_ang.le.883) then
                lebedev_grid_ceil_esp = 33
            else if (n_ang.le.926) then
                lebedev_grid_ceil_esp = 34
            else if (n_ang.le.969) then
                lebedev_grid_ceil_esp = 35
            else if (n_ang.le.1014) then
                lebedev_grid_ceil_esp = 36
            else if (n_ang.le.1059) then
                lebedev_grid_ceil_esp = 37
            else if (n_ang.le.1106) then
                lebedev_grid_ceil_esp = 38
            else if (n_ang.le.1153) then
                lebedev_grid_ceil_esp = 39
            else if (n_ang.le.1202) then
                lebedev_grid_ceil_esp = 40
            else if (n_ang.le.1251) then
                lebedev_grid_ceil_esp = 41
            else if (n_ang.le.1302) then
                lebedev_grid_ceil_esp = 42
            else if (n_ang.le.1862) then
                lebedev_grid_ceil_esp = 43
            else if (n_ang.le.3282) then
                lebedev_grid_ceil_esp = 44
            else if (n_ang.le.5102) then
                lebedev_grid_ceil_esp = 45
            else if (n_ang.le.7322) then
                lebedev_grid_ceil_esp = 46
            else if (n_ang.le.9942) then
                lebedev_grid_ceil_esp = 47
            else if (n_ang.le.12962) then
                lebedev_grid_ceil_esp = 48
            else if (n_ang.le.16382) then
                lebedev_grid_ceil_esp = 49
            else
                lebedev_grid_ceil_esp = 49
            end if
        else if ((force_lebedev.eq.1).or.(n_ang.gt.1202)) then
          if (n_ang.le.0) then
            lebedev_grid_ceil_esp = 0
          else if (n_ang.le.6) then
            lebedev_grid_ceil_esp = 1
          else if (n_ang.le.14) then
            lebedev_grid_ceil_esp = 2
          else if (n_ang.le.26) then
            lebedev_grid_ceil_esp = 3
          else if (n_ang.le.38) then
            lebedev_grid_ceil_esp = 4
          else if (n_ang.le.50) then
            lebedev_grid_ceil_esp = 5
          else if (n_ang.le.78) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_ceil_esp = 6
          else if (n_ang.le.86) then
            lebedev_grid_ceil_esp = 7
          else if (n_ang.le.110) then
            lebedev_grid_ceil_esp = 8
          else if (n_ang.le.146) then
            lebedev_grid_ceil_esp = 9
          else if (n_ang.le.170) then
            lebedev_grid_ceil_esp = 10
          else if (n_ang.le.194) then
            lebedev_grid_ceil_esp = 11
          else if (n_ang.le.248) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_ceil_esp = 12
          else if (n_ang.le.288) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_ceil_esp = 13
          else if (n_ang.le.302) then
            lebedev_grid_ceil_esp = 14
          else if (n_ang.le.350) then
            lebedev_grid_ceil_esp = 15
          else if (n_ang.le.434) then
            lebedev_grid_ceil_esp = 16
          else if (n_ang.le.590) then
            lebedev_grid_ceil_esp = 17
          else if (n_ang.le.770) then
            lebedev_grid_ceil_esp = 18
          else if (n_ang.le.974) then
            lebedev_grid_ceil_esp = 19
          else if (n_ang.le.1202) then
            lebedev_grid_ceil_esp = 20
          else if (n_ang.le.1454) then
            lebedev_grid_ceil_esp = 21
          else if (n_ang.le.1730) then
            lebedev_grid_ceil_esp = 22
          else if (n_ang.le.2030) then
            lebedev_grid_ceil_esp = 23
          else if (n_ang.le.2354) then
            lebedev_grid_ceil_esp = 24
          else if (n_ang.le.2702) then
            lebedev_grid_ceil_esp = 25
          else if (n_ang.le.3074) then
            lebedev_grid_ceil_esp = 26
          else if (n_ang.le.3470) then
            lebedev_grid_ceil_esp = 27
          else if (n_ang.le.3890) then
            lebedev_grid_ceil_esp = 28
          else if (n_ang.le.4334) then
            lebedev_grid_ceil_esp = 29
          else if (n_ang.le.4802) then
            lebedev_grid_ceil_esp = 30
          else if (n_ang.le.5294) then
            lebedev_grid_ceil_esp = 31
          else if (n_ang.le.5810) then
            lebedev_grid_ceil_esp = 32
          else
           if (n_ang.gt.5810) then
             if (myid.eq.0) then
               write(use_unit,'(1X,A,A)') &
               "! Defaulting to maximum number of angular grid points: ", &
               "5810."
             end if
             lebedev_grid_ceil_esp = 32
           end if
          end if
        else
!         Delley style grids; not all are supplied
          if (n_ang.le.0) then
            lebedev_grid_ceil_esp = 0
          else if (n_ang.le.6) then
            lebedev_grid_ceil_esp = 1
          else if (n_ang.le.14) then
            lebedev_grid_ceil_esp = 2
          else if (n_ang.le.26) then
            lebedev_grid_ceil_esp = 3
          else if (n_ang.le.50) then
            lebedev_grid_ceil_esp = 5
          else if (n_ang.le.110) then
            lebedev_grid_ceil_esp = 8
          else if (n_ang.le.194) then
            lebedev_grid_ceil_esp = 11
          else if (n_ang.le.302) then
            lebedev_grid_ceil_esp = 14
          else if (n_ang.le.434) then
            lebedev_grid_ceil_esp = 16
          else if (n_ang.le.590) then
            lebedev_grid_ceil_esp = 17
          else if (n_ang.le.770) then
            lebedev_grid_ceil_esp = 18
          else if (n_ang.le.974) then
            lebedev_grid_ceil_esp = 19
          else if (n_ang.le.1202) then
            lebedev_grid_ceil_esp = 20
          end if

        end if

        end function lebedev_grid_ceil_esp
!******

!------------------------------------------------------------------------------
!****f* esp_grids/grid_ceil_esp
!  NAME
!    grid_ceil_esp
!  SYNOPSIS
      integer function grid_ceil_esp( n_ang )
!  PURPOSE
!    function grid_ceil_esp obtains the next higher number
!    of the Lebedev grid for a given number
!    of angular integration grid points.
!
!  USES

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
           grid_ceil_esp = 0

        else
           do i_grid = 1, n_grids_esp, 1
              grid_ceil_esp = i_grid
              if (n_angular_lebedev_esp(i_grid).ge.n_ang) then
                 exit
              end if
           enddo
        end if

        end function grid_ceil_esp
!******

!------------------------------------------------------------------------------
!****f* grids/grid_floor_esp
!  NAME
!    grid_floor_esp
!  SYNOPSIS
      integer function grid_floor_esp( n_ang )
!  PURPOSE
!    function grid_floor_esp obtains the next lower number
!    of the Lebedev grid for a given number
!    of angular integration grid points.
!  USES

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

        if (n_ang.lt.n_angular_lebedev_esp(1)) then
           grid_floor_esp = 0

        else
           do i_grid = n_grids_esp, 1, -1
              grid_floor_esp = i_grid
              if (n_angular_lebedev_esp(i_grid).le.n_ang) then
                 exit
              end if
           enddo
        end if

        end function grid_floor_esp
!******

!------------------------------------------------------------------------------
!****f* grids/lebedev_grid_floor_esp
!  NAME
!    lebedev_grid_floor_esp
!  SYNOPSIS
        integer function lebedev_grid_floor_esp ( n_ang )
!  PURPOSE
!    function lebedev_grid_ceil_esp obtains the next lower number
!    of the Lebedev grid for a given number
!    of angular integration grid points.
!  USES

        use runtime_choices,only:force_lebedev
        use mpi_tasks,only:myid

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
                lebedev_grid_floor_esp = 0
            else if (n_ang.lt.62) then
                lebedev_grid_floor_esp = 1
            else if (n_ang.lt.72) then
                lebedev_grid_floor_esp = 2
            else if (n_ang.lt.86) then
                lebedev_grid_floor_esp = 3
            else if (n_ang.lt.98) then
                lebedev_grid_floor_esp = 4
            else if (n_ang.lt.114) then
                lebedev_grid_floor_esp = 5
            else if (n_ang.lt.128) then
                lebedev_grid_floor_esp = 6
            else if (n_ang.lt.146) then
                lebedev_grid_floor_esp = 7
            else if (n_ang.lt.163) then
                lebedev_grid_floor_esp = 8
            else if (n_ang.lt.182) then
                lebedev_grid_floor_esp = 9
            else if (n_ang.lt.201) then
                lebedev_grid_floor_esp = 10
            else if (n_ang.lt.222) then
                lebedev_grid_floor_esp = 11
            else if (n_ang.lt.243) then
                lebedev_grid_floor_esp = 12
            else if (n_ang.lt.266) then
                lebedev_grid_floor_esp = 13
            else if (n_ang.lt.289) then
                lebedev_grid_floor_esp = 14
            else if (n_ang.lt.314) then
                lebedev_grid_floor_esp = 15
            else if (n_ang.lt.339) then
                lebedev_grid_floor_esp = 16
            else if (n_ang.lt.366) then
                lebedev_grid_floor_esp = 17
            else if (n_ang.lt.393) then
                lebedev_grid_floor_esp = 18
            else if (n_ang.lt.422) then
                lebedev_grid_floor_esp = 19
            else if (n_ang.lt.451) then
                lebedev_grid_floor_esp = 20
            else if (n_ang.lt.482) then
                lebedev_grid_floor_esp = 21
            else if (n_ang.lt.513) then
                lebedev_grid_floor_esp = 22
            else if (n_ang.lt.546) then
                lebedev_grid_floor_esp = 23
            else if (n_ang.lt.579) then
                lebedev_grid_floor_esp = 24
            else if (n_ang.lt.614) then
                lebedev_grid_floor_esp = 25
            else if (n_ang.lt.649) then
                lebedev_grid_floor_esp = 26
            else if (n_ang.lt.686) then
                lebedev_grid_floor_esp = 27
            else if (n_ang.lt.723) then
                lebedev_grid_floor_esp = 28
            else if (n_ang.lt.762) then
                lebedev_grid_floor_esp = 29
            else if (n_ang.lt.801) then
                lebedev_grid_floor_esp = 30
            else if (n_ang.lt.842) then
                lebedev_grid_floor_esp = 31
            else if (n_ang.lt.883) then
                lebedev_grid_floor_esp = 32
            else if (n_ang.lt.926) then
                lebedev_grid_floor_esp = 33
            else if (n_ang.lt.969) then
                lebedev_grid_floor_esp = 34
            else if (n_ang.lt.1014) then
                lebedev_grid_floor_esp = 35
            else if (n_ang.lt.1059) then
                lebedev_grid_floor_esp = 36
            else if (n_ang.lt.1106) then
                lebedev_grid_floor_esp = 37
            else if (n_ang.lt.1153) then
                lebedev_grid_floor_esp = 38
            else if (n_ang.lt.1202) then
                lebedev_grid_floor_esp = 39
            else if (n_ang.lt.1251) then
                lebedev_grid_floor_esp = 40
            else if (n_ang.lt.1302) then
                lebedev_grid_floor_esp = 41
            else if (n_ang.lt.1862) then
                lebedev_grid_floor_esp = 42
            else if (n_ang.lt.3282) then
                lebedev_grid_floor_esp = 43
            else if (n_ang.lt.5102) then
                lebedev_grid_floor_esp = 44
            else if (n_ang.lt.7322) then
                lebedev_grid_floor_esp = 45
            else if (n_ang.lt.9942) then
                lebedev_grid_floor_esp = 46
            else if (n_ang.lt.12962) then
                lebedev_grid_floor_esp = 47
            else if (n_ang.lt.16382) then
                lebedev_grid_floor_esp = 48
            else
                lebedev_grid_floor_esp = 49
            end if
        else if ((force_lebedev.eq.1).or.(n_ang.ge.1454)) then
          if (n_ang.lt.6) then
            lebedev_grid_floor_esp = 0
          else if (n_ang.lt.14) then
            lebedev_grid_floor_esp = 1
          else if (n_ang.lt.26) then
            lebedev_grid_floor_esp = 2
          else if (n_ang.lt.38) then
            lebedev_grid_floor_esp = 3
          else if (n_ang.lt.50) then
            lebedev_grid_floor_esp = 4
          else if (n_ang.lt.78) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_floor_esp = 5
          else if (n_ang.lt.86) then
            lebedev_grid_floor_esp = 6
          else if (n_ang.lt.110) then
            lebedev_grid_floor_esp = 7
          else if (n_ang.lt.146) then
            lebedev_grid_floor_esp = 8
          else if (n_ang.lt.170) then
            lebedev_grid_floor_esp = 9
          else if (n_ang.lt.194) then
            lebedev_grid_floor_esp = 10
          else if (n_ang.lt.248) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_floor_esp = 11
          else if (n_ang.lt.288) then ! check external/HeoXu.f for explanation of this number
            lebedev_grid_floor_esp = 12
          else if (n_ang.lt.302) then
            lebedev_grid_floor_esp = 13
          else if (n_ang.lt.350) then
            lebedev_grid_floor_esp = 14
          else if (n_ang.lt.434) then
            lebedev_grid_floor_esp = 15
          else if (n_ang.lt.590) then
            lebedev_grid_floor_esp = 16
          else if (n_ang.lt.770) then
            lebedev_grid_floor_esp = 17
          else if (n_ang.lt.974) then
            lebedev_grid_floor_esp = 18
          else if (n_ang.lt.1202) then
            lebedev_grid_floor_esp = 19
          else if (n_ang.lt.1454) then
            lebedev_grid_floor_esp = 20
          else if (n_ang.lt.1730) then
            lebedev_grid_floor_esp = 21
          else if (n_ang.lt.2030) then
            lebedev_grid_floor_esp = 22
          else if (n_ang.lt.2354) then
            lebedev_grid_floor_esp = 23
          else if (n_ang.lt.2702) then
            lebedev_grid_floor_esp = 24
          else if (n_ang.lt.3074) then
            lebedev_grid_floor_esp = 25
          else if (n_ang.lt.3470) then
            lebedev_grid_floor_esp = 26
          else if (n_ang.lt.3890) then
            lebedev_grid_floor_esp = 27
          else if (n_ang.lt.4334) then
            lebedev_grid_floor_esp = 28
          else if (n_ang.lt.4802) then
            lebedev_grid_floor_esp = 29
          else if (n_ang.lt.5294) then
            lebedev_grid_floor_esp = 30
          else if (n_ang.lt.5810) then
            lebedev_grid_floor_esp = 31
          else if (n_ang.eq.5810) then
            lebedev_grid_floor_esp = 32
          else
            if (n_ang.gt.5810) then
              if (myid.eq.0) then
                write(use_unit,'(1X,A,A)') &
                "! Defaulting to maximum number of angular grid points: ", &
                "5810."
              end if
              lebedev_grid_floor_esp = 32
            end if
         end if
      else
!         Delley style grids; not all are supplied
          if (n_ang.lt.6) then
            lebedev_grid_floor_esp = 0
          else if (n_ang.lt.14) then
            lebedev_grid_floor_esp = 1
          else if (n_ang.lt.26) then
            lebedev_grid_floor_esp = 2
          else if (n_ang.lt.50) then
            lebedev_grid_floor_esp = 3
          else if (n_ang.lt.110) then
            lebedev_grid_floor_esp = 5
          else if (n_ang.lt.194) then
            lebedev_grid_floor_esp = 8
          else if (n_ang.lt.302) then
            lebedev_grid_floor_esp = 11
          else if (n_ang.lt.434) then
            lebedev_grid_floor_esp = 14
          else if (n_ang.lt.590) then
            lebedev_grid_floor_esp = 16
          else if (n_ang.lt.770) then
            lebedev_grid_floor_esp = 17
          else if (n_ang.lt.974) then
            lebedev_grid_floor_esp = 18
          else if (n_ang.lt.1202) then
            lebedev_grid_floor_esp = 19
          else
            lebedev_grid_floor_esp = 20
          end if

        end if

        end function lebedev_grid_floor_esp
!******
!-------------------------------------------------------------------------------
!****s* grids/get_logarithmic_grid_esp
!  NAME
!    get_logarithmic_grid_esp
!  SYNOPSIS
        subroutine get_logarithmic_grid_esp ( i_species )
! PURPOSE
!    Subroutine get_logarithmic_grid_esp sets up the logarithmic integration 
!    grids for radial atomic wave functions (basis functions). Input data are 
!    specified in control.in .
!    grid definition: 
!       r(i) = r_grid_min_esp/z * exp[ log(r_grid_inc_esp) * (i-1) ]
!  USES
        use mpi_tasks,only:myid
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
        r_test = r_grid_min_esp(i_species)

        do while ( r_test.le.r_grid_max_esp (i_species) )

          i_grid = i_grid + 1

          r_grid_esp(i_grid, i_species) = r_test
          r_test=r_test*r_grid_inc_esp(i_species)

        end do

        if (i_grid.le.0) then
           if (myid.eq.0) then
              write(use_unit,*) "Illegal logarithmic grid for species ", &
                   i_species, ". r_min > r_max. Aborting."
           end if
           stop
        else
          n_grid_esp(i_species) = i_grid
        end if

        end subroutine get_logarithmic_grid_esp
!******
!-------------------------------------------------------------------------------
!****s* grids/test_radial_grid_esp
!  NAME
!    test_radial_grid_esp
!  SYNOPSIS

      subroutine test_radial_grid_esp (  )

!  PURPOSE
!    Debug subroutine, not normally used.
!    subroutine test_radial_grid_esp performs two simple 1D test integrals with 
!    functions whose characteristic width is 1 bohr: a step function, and an 
!    exponential function
!  USES
      use dimensions,only:n_max_radial,n_species
      use mpi_tasks,only:myid
      use localorb_io,only:localorb_info
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
        do i_radial = 1, n_radial_esp(i_species), 1

           if (myid.eq.0) then
              open (50, file = "grid.dat")
              write (50,*) r_radial_esp (i_radial, i_species), "1.0", &
                   w_radial_esp (i_radial, i_species)
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
        do i_radial = 1, n_radial_esp(i_species), 1
          if (r_radial_esp (i_radial,i_species).le.Radius ) then
            funct(i_radial) = 1.d0
        else
          funct(i_radial) = 0.d0
        end if
        enddo
        sum = 0.
        do i_radial = 1, n_radial_esp(i_species), 1
          sum = sum + funct(i_radial) * w_radial_esp(i_radial,i_species)

        enddo

        if (myid.eq.0) then
           write(use_unit,*) "* Step Integral: ", sum, ", should be 1."
        end if

! exponential integral
        do i_radial = 1, n_radial_esp(i_species), 1
           funct(i_radial) = exp (- r_radial_esp(i_radial,i_species) )
        enddo
        sum = 0.
        do i_radial = 1, n_radial_esp(i_species), 1
          sum = sum + funct(i_radial) * w_radial_esp(i_radial,i_species)
        enddo

        if (myid.eq.0) then
           write(use_unit,*) "* Exponential Integral: ", sum, ", should be 1."
        end if

      enddo

      write(use_unit,*)

      end subroutine test_radial_grid_esp
!******

!****s* esp_grids/cleanup_grids_esp
!  NAME
!     cleanup_grids_esp
!  SYNOPSIS
        subroutine cleanup_grids_esp( )
!  PURPOSE
!     deallocation of all arrays pertaining to the grids. 
! USES

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


        if ( allocated (r_grid_min_esp) ) then
          deallocate ( r_grid_min_esp )
        end if
        if ( allocated (r_grid_max_esp ) ) then
          deallocate ( r_grid_max_esp  )
        end if
        if ( allocated (r_grid_inc_esp) ) then
          deallocate ( r_grid_inc_esp )
        end if
        if ( allocated (log_r_grid_inc_esp) ) then
          deallocate ( log_r_grid_inc_esp )
        end if

        if ( allocated (n_grid_esp ) ) then
          deallocate ( n_grid_esp )
        end if
        if ( allocated (r_grid_esp ) ) then
          deallocate ( r_grid_esp )
        end if

        if ( allocated (n_radial_esp ) ) then
          deallocate ( n_radial_esp )
        end if
        if ( allocated (scale_radial_esp ) ) then
          deallocate ( scale_radial_esp )
        end if
        if ( allocated (angular_limit_esp ) ) then
          deallocate ( angular_limit_esp )
        end if
        if ( allocated (angular_min_esp ) ) then
          deallocate ( angular_min_esp )
        end if
        if ( allocated (angular_acc_esp ) ) then
          deallocate ( angular_acc_esp )
        end if
        if ( allocated (r_radial_esp ) ) then
          deallocate ( r_radial_esp )
        end if
        if ( allocated (w_radial_esp ) ) then
          deallocate ( w_radial_esp )
        end if
        if ( allocated (n_angular_esp ) ) then
          deallocate ( n_angular_esp )
        end if
        if ( allocated (r_angular_esp ) ) then
          deallocate ( r_angular_esp )
        end if
        if ( allocated (w_angular_esp ) ) then
          deallocate ( w_angular_esp )
        end if
        if (allocated (n_division_esp ) ) then
           deallocate (n_division_esp )
        end if
        if (allocated (division_boundaries_esp)) then
           deallocate (division_boundaries_esp)
        end if
        if (allocated(n_angular_lebedev_esp)) then
           deallocate(n_angular_lebedev_esp)
        end if
        if (allocated(r_angular_lebedev_esp)) then
           deallocate(r_angular_lebedev_esp)
        end if
        if (allocated(w_angular_lebedev_esp)) then
           deallocate(w_angular_lebedev_esp)
        end if
        if (allocated(n_division_lebedev_esp)) then
           deallocate(n_division_lebedev_esp)
        end if
        if (allocated(division_boundaries_lebedev_esp)) then
           deallocate(division_boundaries_lebedev_esp)
        end if
        if (allocated(n_ang_lebedev_esp)) then  !SAG
           deallocate(n_ang_lebedev_esp)
        end if
        if (allocated(r_ang_lebedev_esp)) then
           deallocate(r_ang_lebedev_esp)
        end if
        if (allocated(w_ang_lebedev_esp)) then
           deallocate(w_ang_lebedev_esp)
        end if




        if (allocated(local_ylm_tab_esp)) then
          deallocate ( local_ylm_tab_esp )
        end if
        if (allocated(fixed_grid_index_esp)) then
           deallocate(fixed_grid_index_esp)
        end if

        if (allocated(lebedev_grid_index_esp)) then
           deallocate(lebedev_grid_index_esp)
        end if

        end subroutine cleanup_grids_esp
!******

!****f* grids/invert_radial_grid_esp
!  NAME
!    invert_radial_grid_esp
!  SYNOPSIS

        real*8 function invert_radial_grid_esp(r_current, n_scale, r_scale)

!  PURPOSE
!    Function invert_radial_grid_esp takes a given real-space radius value
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
!    o invert_radial_grid_esp -- radial grid index
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        invert_radial_grid_esp = 1.0d0 - exp( - r_current/r_scale)
        invert_radial_grid_esp = &
          dble(n_scale + 1) * sqrt(invert_radial_grid_esp)
      end function invert_radial_grid_esp
!******
!-------------------------------------------------------------------------------
!****f* esp_grids/invert_log_grid_esp
!  NAME
!    invert_log_grid_esp
!  SYNOPSIS

        real*8 function invert_log_grid_esp ( r_current,r_min,scale)

!  PURPOSE
!    Function invert_log_grid_esp takes a given real-space radius value
!    and returns a (fractional) number i(r) which indicates its position
!    between the grid points i_grid of the logarithmic grid.
!
!    i(r) = 1 + [ln(r / r_grid_min_esp)]/[ln(r_grid_inc)]
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
!    o invert_log_grid_esp -- logarithmic grid index
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



        invert_log_grid_esp = log (r_current/r_min)
        invert_log_grid_esp = 1.d0 + invert_log_grid_esp/log(scale)
        end function invert_log_grid_esp
!******
!------------------------------------------------------------------------------
!****f* wsp_grids/invert_log_grid_p2_esp
!  NAME
!    invert_log_grid_p2_esp
!  SYNOPSIS
        real*8 function invert_log_grid_p2_esp( r_current, i_species)
! PURPOSE
!   Function invert_log_grid_p2_esp takes a given real-space radius value
!   and returns a (fractional) number i(r) which indicates its position
!   between the grid points i_grid of the logarithmic grid.
!
!   i(r) = 1 + [ln(r / r_grid_min_esp)]/[ln(r_grid_inc)]
!
!   This does not care about scales, it only takes a species number and looks 
!   up its own parameters.
!
!  ARGUMENTS

        real*8  :: r_current
        integer :: i_species

!  INPUTS
!    o r_current -- radius to be inverted
!    o i_species -- number of the species
!
! OUTPUTS
!    o invert_log_grid_p2_esp -- logarithmic grid index
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        invert_log_grid_p2_esp = log (r_current/r_grid_min_esp(i_species))
        invert_log_grid_p2_esp = 1.d0 + invert_log_grid_p2_esp/&
        log_r_grid_inc_esp(i_species)
        end function invert_log_grid_p2_esp
!******
!-------------------------------------------------------------------------------
!****f* esp_grids/get_radial_weight_esp
!  NAME
!    get_radial_weight_esp
!  SYNOPSIS

      real*8 function get_radial_weight_esp(index, scale, n_max)

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
!    o get_radial_weight_esp -- di/dr for given radial shell
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      get_radial_weight_esp = (dble(n_max + 1) / index) ** 2.d0
      get_radial_weight_esp = index / (2.d0 * scale) * &
           (get_radial_weight_esp - 1)
      end function get_radial_weight_esp
!******
!-------------------------------------------------------------------------------
!****s* grids/get_local_ylm_tab_esp
!  NAME
!    get_local_ylm_tab_esp
!  SYNOPSIS
        subroutine get_local_ylm_tab_esp ( l_hartree )
!  PURPOSE
!       tabulate Ylm function at each integration point around each species
!       for the benefit of update_hartree_potential later on
!  USES

        use dimensions
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
           do i_radial = 1, n_radial_esp( i_species ), 1
              max_lebedev_index = &
                   MAX( lebedev_grid_index_esp( i_radial, i_species ), &
                   max_lebedev_index )
           end do
        end do
        if (.not.allocated(local_ylm_tab_esp)) then
           allocate ( local_ylm_tab_esp ( (l_pot_max+1)**2, n_max_angular, &
                max_lebedev_index ) )
           local_ylm_tab_esp = 0.0d0
        end if

        do i_species = 1, n_species, 1
           do i_radial = 1, n_radial_esp( i_species ), 1
              i_lebedev = lebedev_grid_index_esp( i_radial, i_species )
              do i_angular = 1, n_angular_esp( i_radial, i_species ), 1
!                dist_tab = r_radial_esp( i_radial, i_species )
                 dir_tab(:) = &
                      r_angular_lebedev_esp( : , i_angular, i_lebedev )

!               compute trigonometric functions of spherical coordinate angles
!               of current integration point
                 call tab_local_trigonom &
                ( dir_tab, trigonom_tab &
                )

                call tab_local_ylm &
                ( trigonom_tab, l_hartree(i_species), l_pot_max, &
                  local_ylm_tab_esp(:,i_angular,i_lebedev) )

              enddo
           enddo
        enddo

      end subroutine get_local_ylm_tab_esp
!******
!-------------------------------------------------------------------------
!****s* esp_grids/invert_log_grid_vector_esp
!  NAME
!    invert_log_grid_vector_esp
!  SYNOPSIS
        subroutine invert_log_grid_vector_esp ( r_current,r_min,scale,&
                   n_points,out_value )
!  PURPOSE
!    Subroutine invert_log_grid_vector_esp takes a set of real-space radiuses
!    and returns (fractional) numbers i(r) which indicates its position
!    between the grid points i_grid of the logarithmic grid.
!
!    i(r) = 1 + [ln(r / r_grid_min_esp)]/[ln(r_grid_inc)]
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
        end subroutine invert_log_grid_vector_esp
!******
!****s* FHI-aims/esp_get_n_compute_dens
!  NAME
!    esp_get_n_compute_dens
!  SYNOPSIS 
subroutine esp_get_n_compute_dens( partition_tab )
!  PURPOSE
!    Computes the threadwise maximal and average number of non-zero basis functions.
!  USES
  use dimensions,only: n_basis_fns, n_centers_basis_T, n_centers_basis_I,&
                       n_periodic, n_centers, n_centers_integrals
!  use esp_grids, only:  batches_esp,n_my_batches_esp, n_max_compute_dens_esp,&
!                        n_max_compute_fns_dens_esp,&
!                        n_max_compute_atoms_esp,got_n_compute_maxes_esp,&
!                        n_full_points_esp,n_max_batch_size_esp
  use pbc_lists,only:inv_centers_basis_integrals,centers_basis_integrals                     
  use runtime_choices, only: use_metis_batch_distribution,use_load_balancing
  use localorb_io, only: use_unit,OL_norm,OL_low,output_priority,&
                         localorb_info
  use mpi_tasks,only:myid
  use load_balancing,only:use_batch_permutation
  use synchronize_mpi,only:sync_find_max

  implicit none
!  ARGUMENTS
  real*8, dimension(n_full_points_esp) :: partition_tab
! the following parameters are needed for the automatic selection of charge density/force update


!  INPUTS
!    o partition_tab -- the partition tab
!  OUTPUT
!    n_max_compute_dens_esp, n_max_compute_fns_dens_esp, n_max_compute_atoms_esp
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
! SOURCE


 ! real*8, dimension( n_full_points_hamiltonian_integrals) :: hamiltonian_partition_tab
  
  ! locals
  integer :: i_my_batch, i_index
  integer :: i_point

  integer :: n_compute_c, n_compute_a
  integer :: i_full_points_A, i_full_points_C, i_full_points_2C
  integer :: i_full_points, i_full_points_2, i_full_points_3

  real*8 :: coord_current(3)
  real*8 :: dist_tab(n_centers_integrals, n_max_batch_size_esp)
  real*8 :: dist_tab_sq(n_centers_integrals, n_max_batch_size_esp)
  real*8 :: dir_tab(3,n_centers_integrals, n_max_batch_size_esp)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

  integer :: n_centers_max, n_points
  integer :: n_max

  integer,dimension(:),allocatable :: i_basis

  character*200 :: info_str
  integer :: mpierr

  logical, dimension(:), allocatable :: my_functions



  if (got_n_compute_maxes_esp) return

  write (info_str,'(2X,A)') "Obtaining max. number of non-zero basis functions in each batch (esp_get_n_compute_dens)."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)


  n_centers_max = MAX(n_centers_basis_I, n_centers_basis_T)

  allocate(i_basis(n_centers_max))

  n_max_compute_atoms_esp = 0

  n_max_compute_dens_esp = 0
  n_max_compute_fns_dens_esp = 0



  i_full_points_C = 0
  i_full_points_2C = 0
  i_full_points_A = 0

  i_basis_fns_inv = 0

  if (use_metis_batch_distribution) then
     allocate(my_functions(n_centers_basis_I))
     my_functions = .false.
  end if

  do i_my_batch = 1, n_my_batches_esp, 1

        n_compute_c = 0
        n_compute_a = 0
        i_basis = 0
     
        i_point = 0

        ! loop over one batch
        do i_index = 1, batches_esp(i_my_batch)%size_esp, 1

           i_full_points_2C = i_full_points_2C + 1
           
           if (partition_tab(i_full_points_2C).gt.0.d0) then

              i_point = i_point+1

              ! get current integration point coordinate
              coord_current(:) = batches_esp(i_my_batch) % points_esp(i_index) % coords_esp(:)
              
              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if
              
              ! compute atom-centered coordinates of current integration point,
              ! as viewed from all atoms
              call tab_atom_centered_coords_p0 &
                   ( coord_current,  &
                   dist_tab_sq(1,i_point),  &
                   dir_tab(1,1,i_point), &
                   n_centers_integrals, centers_basis_integrals )
              
              ! determine which basis functions are relevant at current integration point,
              ! and tabulate their indices
              
              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
              call prune_basis_p0X &
                   ( dist_tab_sq(1,i_point), &
                   n_compute_a, n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, i_point  )
!NEC_CB Using version which saves massively indirect addressed arrays

              
           end if
        enddo  ! end loop over one part of the angular division              

        n_points = i_point
        
        batches_esp(i_my_batch)%batch_n_compute_esp = n_compute_a
        allocate(batches_esp(i_my_batch)%batch_i_basis_esp(n_compute_a))
        batches_esp(i_my_batch)%batch_i_basis_esp = i_basis(1:n_compute_a)

        if (use_metis_batch_distribution) then
           my_functions(i_basis(1:n_compute_a)) = .true.
        end if

        n_max_compute_dens_esp = MAX(n_compute_c, n_max_compute_dens_esp)

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_a.gt.0) then

           i_point = 0

           ! loop over one division of the angular grid
           do i_index = 1, batches_esp(i_my_batch)%size_esp, 1

              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points_C = i_full_points_C + 1
              
              if (partition_tab(i_full_points_C).gt.0.d0) then
                 
                 i_point = i_point+1

                 n_compute_atoms = 0
                 n_compute_fns = 0
!                 i_basis_fns_inv = 0

                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                 ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                 ! without any copying and without doing any unnecessary operations. 
                 ! The price is that the interface is no longer explicit in terms of physical 
                 ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                 call prune_radial_basis_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dist_tab(1,i_point), &
                      dir_tab(1,1,i_point), &
                      n_compute_atoms, atom_index, atom_index_inv, &
                      n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                      i_atom_fns, spline_array_start, spline_array_end, &
                      n_centers_integrals, centers_basis_integrals)
                 
                 n_max_compute_fns_dens_esp = MAX(n_compute_fns, n_max_compute_fns_dens_esp)

                 n_max_compute_atoms_esp = MAX(n_compute_atoms, n_max_compute_atoms_esp)

              end if

           end do ! end loop over a batch
           
        else
          ! must increment grid counter in this case too, else we get an inconsistent
          ! partition_tab and all sorts of trouble

           i_full_points_C = i_full_points_C + batches_esp(i_my_batch)%size_esp

        end if ! end if (n_compute.gt.0)
        
     ! end if ! end distribution of tasks over threads
     
  end do ! end loop over batches

  if(use_load_balancing) then
     ! When load balancing is used, we need the global max numbers since batches are permuted
     n_max = n_max_compute_atoms_esp;        call sync_find_max(n_max, n_max_compute_atoms_esp)
     n_max = n_max_compute_dens_esp;     call sync_find_max(n_max, n_max_compute_dens_esp)
     n_max = n_max_compute_fns_dens_esp; call sync_find_max(n_max, n_max_compute_fns_dens_esp)
     write(info_str,'(2X,A,I8)') "| Maximal number of non-zero basis functions: ",n_max_compute_dens_esp
     call localorb_info(info_str)
  endif



  if( allocated(i_basis))then
     deallocate(i_basis)
  end if
  if (allocated(my_functions)) then
     deallocate(my_functions)
  end if
  
  
end subroutine esp_get_n_compute_dens
!******

      end module esp_grids
