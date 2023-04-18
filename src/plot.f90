!****h* FHI-aims/plot
!  NAME
!    plot 
!  SYNOPSIS

      module plot
! PURPOSE
!  Module plot handles graphical output options:
!
!  Gaussian-type cube files:
!  * total charge density
!  * delta charge density
!  * eigenstates
!
!  Subroutines:
!  * allocate_plot
!  * read_cube_parameters
!  * cleanup_plot
!
!  USES
      use localorb_io
      use dimensions
      use geometry
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

!     global variable declarations - exported to other program parts
!     All variables have the array-size (n_cubes), and the information is
!     stored seperately for each cube.
!
!     cube_type  :     specifies what kind of (3D) data is plotted in cube file # n_cube
!                  - "total_density"
!                  - "eigenstate"
!                  - "eigenstate_imag"
!                  - "eigenstate_density"
!                  - "long range potential"
!                  - "potential"
!                  - "elf"
!     cube_origin:     Cartesian origin of the voumentric output
!     cube_edge_unit:  Size of each voxel step
!     cube_edge_steps: Number of voxel steps for each edge 
!
!     cube_stm: 
!     cube_fermi_ref:
!     cube_index:      Identifier of the plotted cube_type
!                 - type="eigenstate" or "eigenstate_density" : index is the number of the eigenstate
!     TODO: cube_index is no longer needed, as it is redundant with cube_state,  and only kept for 
!           backwards compatibility with output_cubes_files.f90
!           At some point, both output_cube_files and output_cube_files_p1 should go, and this
!           keyword with it.
!     cube_state:      contains k-point and state of the cube file
!     cube_spin:       contains the factors to multiply for each spin.
!     cube_elf         specifies type of ELF calculation:
!                      0 - Savin et al.
!                      1 - Becke-Edgecombe for spin 1
!                      2 - Becke-Edgecombe for spin 2
!                      anything else - Kohout-Savin variant
!
!     cube_filename:   stores the name of the output file
!     cube_format:     governs the format of the output file
!
!     flag_origin:     true if user defined input has been detected for the cube origin
!     flag_edge:       true if user defined input has been detected for the cube edges
!     flag_divisor:    true if user defined input has been detected for the cube divisor
!
!     cube_divisor:    governs the distribution of cubes into minicubes, which are then parallalized
!
!     cube_type_needs_densmat:      true if the current cube type requires the density matrix
!     cube_type_needs_eigenvectors: true if the current cube type requires KS eigenvectors
!     cube_type_needs_soc_eigenvectors: true if the current cube type requires SOC eigenvectors

      character*20, dimension(:), allocatable      :: cube_type 
      
      !Geometry information
      real*8, dimension(:,:), allocatable          :: cube_origin
      real*8, dimension(:,:,:), allocatable        :: cube_edge_unit
      integer, dimension(:,:), allocatable         :: cube_edge_steps

      !Information about the content
      real*8, dimension(:,:), allocatable          :: cube_stm
      integer, dimension(:), allocatable           :: cube_elf
      logical, dimension(:), allocatable           :: cube_fermi_ref 
      integer, dimension(:), allocatable           :: cube_index
      integer, dimension(:,:), allocatable         :: cube_state
      integer, dimension(:,:), allocatable         :: cube_spin


      !output options
      character*100, dimension(:), allocatable     :: cube_filename
      character*100, dimension(:), allocatable     :: cube_inputname !Save input so that it can be changed later
      character*20, dimension(:), allocatable      :: cube_format !Currently supported options are "cube", "gOpenMol", and "xsf"
      logical, dimension(:), allocatable           :: cube_average

    
      !Checks if user inputs are made or default are needed
      logical, dimension(:), allocatable           :: flag_origin 
      logical, dimension(:), allocatable           :: flag_edge
      logical, dimension(:), allocatable           :: flag_divisor

      !Performance option - governs distribution in minicubes
      integer, dimension(:), allocatable           :: cube_divisor 

      !Flags governing interaction with remaining code
      logical, dimension(:), allocatable           :: cube_type_needs_densmat
      logical, dimension(:), allocatable           :: cube_type_needs_eigenvectors
      logical, dimension(:), allocatable           :: cube_type_needs_soc_eigenvectors

      character*200, private :: info_str

      contains
!******	
!------------------------------------------------------------------------
!****s* plot/allocate_plot
!  NAME
!    allocate_plot
!  SYNOPSIS
        subroutine allocate_plot ( )

! PURPOSE
!  Subroutine allocate_plot allocates the variables needed for plottable output.
!  This allocation is done once at the beginning of the code, since the input is read
!  from file control.in
! USES
        use dimensions
        use aims_memory_tracking, only : aims_allocate
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
        implicit none

!       begin work

!       check here in case we introduce other plottable output later ...
        if (use_cube_output) then
          ! AJL, Oct2018: aims_deallocate can't handle text
          allocate ( cube_type(n_cube) )
          call aims_allocate ( cube_index, n_cube, "cube_index" )
          call aims_allocate ( cube_fermi_ref, n_cube, "cube_fermi_ref")
          cube_fermi_ref(:) = .false.
          call aims_allocate ( cube_elf,n_cube, "cube_elf" )
          ! These can't be handled by aims_deallocate currently either
          allocate ( cube_origin(3,0:n_cube) ) !Include "0" as default
          allocate ( cube_edge_unit(3,3,0:n_cube) ) !Include "0" as default
          allocate ( cube_edge_steps(3,0:n_cube) ) !Include "0" as default
          allocate ( cube_divisor(0:n_cube)) !Include "0" as default
          call aims_allocate ( cube_spin,2,n_cube, "cube_spin" )
          call aims_allocate ( cube_state,2,n_cube, "cube_state" )
          call aims_allocate ( cube_stm,3,n_cube, "cube_stm" )
          allocate ( cube_filename(n_cube)) 
          allocate ( cube_inputname(n_cube)) 
          call aims_allocate ( flag_origin,n_cube, "flag_origin" )
          call aims_allocate (flag_edge,n_cube, "flag_edge" )
          allocate ( cube_format(n_cube)) 
          call aims_allocate ( flag_divisor,n_cube, "flag_advisor")
          call aims_allocate ( cube_type_needs_densmat,n_cube, "cube_type_needs_densmat")
          call aims_allocate ( cube_type_needs_eigenvectors,n_cube, "cube_type_needs_eigenvector")
          call aims_allocate ( cube_type_needs_soc_eigenvectors,n_cube, "cube_type_needs_soc_eigenvectors")
          call aims_allocate ( cube_average,n_cube, "cube_average" )
        end if

        end subroutine allocate_plot
!******	
!---------------------------------------------------------------------------
!****s* plot/read_cube_parameters
!  NAME
!    read_cube_parameters
!  SYNOPSIS
        subroutine read_cube_parameters ( i_cube )
!  PURPOSE
!    Subroutine read_cube_parameters reads the geometry parameters of an output
!    cube for requested cube output number i_cube
!  USES
        use constants
        use mpi_tasks, only: aims_stop, myid
        use runtime_choices
        implicit none

! ARGUMENTS

        integer :: i_cube

!  INPUTS
!    o i_cube - number of the cube file 
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE



!  local variables

!  desc_str : line identifier in control.in
!  i_code   : return status flag from read operation.

        character*20 desc_str
        character*132 inputline
        integer i_code

        integer :: i_edge
       character(*), parameter :: func = 'read_cube_parameters'

!     counters

      integer :: i_coord

!  begin work

!set defaults
      flag_origin(i_cube) = .false.
      flag_edge(i_cube) = .false.
      i_edge = 0
      cube_filename(i_cube) =""
      cube_inputname(i_cube) =""
      cube_format(i_cube)='cube'

      cube_elf(i_cube) = 0
      cube_state(1,i_cube)=1
      cube_state(2,i_cube)=1
      cube_spin(1,i_cube)=1
      cube_spin(2,i_cube)=1
      cube_divisor(i_cube)=10
      flag_divisor(i_cube)=.false.
      cube_average(i_cube)=.false.
        write(info_str,'(2X,A)') "| Cube output options: "
        call localorb_info(info_str)

!       read cube information
      if (flag_out_elec_real) then 
         cube_average(i_cube) = .true. 
      endif 
      lineloop: do

         read(7,'(A)',iostat=i_code) inputline
         if(i_code<0) then
            backspace(7)
            exit lineloop        ! end of file reached
         end if
         if(i_code>0) then
            write(use_unit,'(2X,A)') &
               "Error reading file 'control.in'..."
            stop
         end if

         read (inputline,*,iostat=i_code) desc_str

         if (i_code.ne.0) cycle lineloop

         if (desc_str(1:1) == "#") cycle lineloop

         if (desc_str.eq.'cube') then

            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'elf_type') then

               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                    cube_elf(i_cube)
               if (cube_elf(i_cube).eq.0) then
                  write(info_str,'(2X,A)') &
                       "| Savin et al. ELF variant will be calculated"
               else if (n_spin.eq.1) then
                  write(info_str,'(2X,A,A)') &
                       "| Only Savin et al. variant of ELF is avaliable for spin-unpolarized systems,",&
                       " defaulting to this variant"
                  cube_elf(i_cube) = 0
               else if (cube_elf(i_cube).eq.1.or.cube_elf(i_cube).eq.2) then
                  write(info_str,'(2X,A)') &
                       "| Becke-Edgecombe ELF variant requested for spin-channel 1"
               else
                  write(info_str,'(2X,A)') &
                       "| Kohout-Savin variant of ELF is requested"
               endif
               call localorb_info(info_str)

            else if (desc_str.eq.'state') then

! cube_state(1,i_cube) is spin
! cube_state(2,i_cube) is k_point

               if(n_periodic.eq.0) then
                  read(inputline,*,end=88,err=99) desc_str, desc_str, &
                       cube_state(1,i_cube)
                  write(info_str,'(2X,A,1X,I1.1)') &
                       "| Eigenstate spin: ", &
                       cube_state(1,i_cube)
                  cube_state(2,i_cube)=1
               else
                  read(inputline,*,end=88,err=99) desc_str, desc_str, &
                       cube_state(1,i_cube), cube_state(2,i_cube)

                  write(info_str,'(2X,A,1X,I1.1,1X,I4.4)') &
                       "| Eigenstate spin and k-point: ", &
                       cube_state(1,i_cube), cube_state(2,i_cube)
               endif
               call localorb_info(info_str)

            else if (desc_str.eq.'kpoint') then
                  read(inputline,*,end=88,err=99) desc_str, desc_str, &
                      cube_state(2,i_cube)
                  write(info_str,'(2X,A,1X,I1.1)') &
                       "| Eigenstate k-point: ", &
                       cube_state(2,i_cube)
                   call localorb_info(info_str)

            else if (desc_str.eq.'spinstate') then
! comment cl
!                 if (calculate_perturbative_soc) then
!                   call localorb_info("")
!                   write(info_str, "(2X,A)") '* You have specified a (collinear) spin channel when plotting cube files for & 
!                                             &a spin-orbit-coupled eigenvector.  This is not possible, as SOC is non-collinear. &
!                                             &Exiting.' 
!                   call aims_stop_coll(info_str, func)
!                 endif

                  read(inputline,*,end=88,err=99) desc_str, desc_str, &
                       cube_state(1,i_cube)
                  write(info_str,'(2X,A,1X,I1.1)') &
                       "| Eigenstate spin: ", &
                       cube_state(1,i_cube)
                   call localorb_info(info_str)

            else if (desc_str.eq.'spin'.or. desc_str.eq.'spinmask') then

! cube_spin(1,i_cube) is coefficient for i_spin = 1
! cube_spin(2,i_cube) is coefficient for i_spin = 2

               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                    cube_spin(1,i_cube), cube_spin(2,i_cube)

               write(info_str,'(2X,A,1X,I1.1,1X,I4.4)') &
                    "| Electron density spin mask: ", &
                    cube_spin(1,i_cube), cube_spin(2,i_cube)

               call localorb_info(info_str)

            else if (desc_str.eq.'origin') then

              read(inputline,*,end=88,err=99) desc_str, desc_str, &
              (cube_origin(i_coord,i_cube), i_coord=1,3,1)


              do i_coord = 1,3,1
                cube_origin(i_coord,i_cube) = &
                cube_origin(i_coord,i_cube)/bohr
              enddo

              flag_origin(i_cube) = .true.


            else if (desc_str.eq.'edge') then

              i_edge = i_edge+1
              if (i_edge.gt.3) then
                write(use_unit,'(1X,A)') &
                "* Too many cube edges specified - please correct!"
                stop
              end if

              read(inputline,*,end=88,err=99) desc_str, desc_str, &
              cube_edge_steps(i_edge,i_cube), &
              (cube_edge_unit(i_coord,i_edge,i_cube), i_coord=1,3,1)

              if (i_edge.eq.3) then
                flag_edge(i_cube) = .true.
              end if

            elseif (desc_str.eq.'filename') then
                read(inputline,*,end=88,err=99) desc_str, &
                    desc_str, cube_inputname(i_cube)
                    cube_filename(i_cube)=cube_inputname(i_cube)
        
            elseif (desc_str.eq.'format') then
                read(inputline,*,end=88,err=99) desc_str, &
                    desc_str, cube_format(i_cube)

            elseif (desc_str.eq.'divisor') then
               read(inputline,*,end=88,err=99) desc_str, &
                    desc_str, cube_divisor(i_cube)
               if (cube_divisor(i_cube).gt.45) then
                 call aims_stop("Cube divisor must be smaller than 45")
                endif
               flag_divisor(i_cube) = .true.
            else if (desc_str.eq.'gOpenMol_format') then

              read(inputline,*,end=88,err=99) desc_str,   desc_str, flag_gOpenMol_format

            elseif (desc_str.eq.'average') then 
                cube_average(i_cube) = .true.
            else

              write(use_unit,'(1X,A,A,A)') &
              "* Unknown cube specification ", desc_str, &
              ". Please correct."
              stop

            end if

          else
!           must have reached the end of the cube description.

            backspace(7)
            exit lineloop

          end if

        end do lineloop

!       check for output completeness

      if (myid.eq.0) then
        if (.not.flag_origin(i_cube)) then
!          write(use_unit,'(1X,A,A)') &
!          "* Cube origin was not provided following cube request.", &
!          " Please correct."
          write(use_unit,'(1X,A,A)') &
          "* Cube origin was not provided following cube request.", &
          " Will be using default value."
        end if

        if (.not.flag_edge(i_cube)) then
          if(i_edge.eq.0) then
               write(use_unit,'(1X,A,A)') &
               "* Cube edges were not provided following cube request.", &
               " Defaults will be used."
          elseif (i_edge.gt.0 .and. i_edge.le.3) then 
               write(use_unit,'(1X,A,A)') &
               "* Not all three cube edges were provided ", &
               "following cube request. Please correct."
               stop
          endif
        end if
      endif


!       Nonsense checks
       if (myid.eq.0) then
         if (flag_edge(i_cube) ) then  !Check only if no default will be taken
              if ( cube_edge_steps(i_edge,i_cube).le.0) then
                write(use_unit,'(1X,A,A)') &
                "* Zero or less grid steps are not possible. ", &
                "Please correct."
                stop
              end if
              if ( (cube_edge_unit(1,i_edge,i_cube).eq.0.d0) .and. &
                   (cube_edge_unit(2,i_edge,i_cube).eq.0.d0) .and. &
                   (cube_edge_unit(3,i_edge,i_cube).eq.0.d0) ) &
                then
                write(use_unit,'(1X,A,A)') &
                "* Zero-length grid units are not possible. ", &
                "Please correct."
                stop
              end if
          endif
         
          select case (cube_format(i_cube))
                 case ('cube') 
                    info_str=' Standard cube format requested.'
                    if (myid.eq.0) call localorb_info(info_str)
                 case ('gOpenMol') 
                    info_str=' gOpenMol cube format requested.'
                    if (myid.eq.0) call localorb_info(info_str)
                 case('xsf') 
                    info_str=' XCrysden format requested'
                    if (myid.eq.0) call localorb_info(info_str)
                 case default
                    info_str = ' Unknown cube format requested'
                    if (myid.eq.0) call aims_stop(info_str,func)
          end select
          if(flag_gOpenMol_format) then !Reason: Setting this keyword would give inconsistent file types. 
            write(unit=info_str,fmt='(A,A)')& 
                 'Keyword gOpenMol no longer supported.' ,  &
                 'Please use the keyword "cube format gOpenMol" instead'
            call aims_stop(info_str,func)
           
         endif
   
        if(cube_state(1,i_cube).gt.2.or.cube_state(1,i_cube).lt.1) then
             info_str = 'Cube file with illegal spin requested.' 
             call aims_stop(info_str,func)
        endif
        if(cube_state(1,i_cube).gt.n_spin) then
          info_str ='Cubefiles: Spin = 2 not possible for spin=none.' 
!          call aims_stop(info_str,func)
        endif


     ! the following nonsense check should really be here, but
     ! apparently, the number of k-points is not known yet. 
        if(cube_state(2,i_cube).gt.product(n_k_points_xyz)) then
            write(unit=info_str,fmt=('(A,I0,A,I0,A)')) &
             'Cube file output at k-point ', cube_state(2,i_cube),  &
             ' requested, but there are only ', product(n_k_points_xyz),  &
             ' k-points specified. Please correct.' 
            call aims_stop(info_str,func)
        elseif(cube_state(2,i_cube).lt.0) then
          info_str ='Cubefiles: Negative k-point requested.'
          call aims_stop(info_str,func)
        elseif(cube_state(2,i_cube)==0) then
          info_str = 'k-point 0 requested.'
          call localorb_info(info_str)
          info_str= 'Cube output will contain the sum over all  &
                  & k_points specified for scf'
          call localorb_info(info_str)
         endif
      endif

      if(cube_state(2,i_cube).ne.1) then
         write(info_str,*) 'WARNING: Off-Gamma eigenstate output requested'
         call localorb_info(info_str)
         write(info_str,*) 'There are currently doubts on the correctness of the results' 
         call localorb_info(info_str)
         write(info_str,*) 'Therefore, this feature will be deactivated for now' 
         call localorb_info(info_str)
         call aims_stop('Sorry for the inconvenience.')
      end if
      
      

!Block moved to when the cube-file is really created (so defaults are accounted for )
!              write(info_str,'(2X,A,I3,A,I6,A,3(1X,F11.4))') &
!              "| Grid step direction ", i_edge, ": ", &
!              cube_edge_steps(i_edge,i_cube), " steps of ", &
!              (cube_edge_unit(i_coord,i_edge,i_cube), i_coord=1,3,1)
!             call localorb_info(info_str)

!         Block moved to when the cube-file is really created
!              write(info_str,'(2X,A,3(1X,F11.4))') &
!              "| Origin: ", &
!              (cube_origin(i_coord,i_cube), i_coord=1,3,1)

!        Convert Angstrom input into atomic units
              do i_coord = 1,3,1
                do i_edge = 1,3,1
                cube_edge_unit(i_coord,i_edge,i_cube) = &
                cube_edge_unit(i_coord,i_edge,i_cube)/bohr
                enddo
              enddo

             call localorb_info(info_str)

        return

      88 continue
         write(use_unit,'(A)') "Syntax error reading 'output cube' block from 'control.in' (missing arguments)"
         write(use_unit,'(A)') "in line: '"//trim(inputline)//"'"
         stop

      99 continue
         write(use_unit,'(A)') "Syntax error reading 'output cube' block from 'control.in'"
         write(use_unit,'(A)') "in line: '"//trim(inputline)//"'"
         stop

      end subroutine read_cube_parameters
!******	
!--------------------------------------------------------------------------------
!****s* plot/cleanup_plot
!  NAME
!    cleanup_plot
!  SYNOPSIS
        subroutine cleanup_plot ( )
!  PURPOSE
!    Subroutine cleanup_plot deallocates everything to do with geometry
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  USES 
        use aims_memory_tracking, only : aims_deallocate
!  SOURCE

        implicit none

!  begin work

        if ( allocated (cube_type) )                   deallocate ( cube_type)
        if ( allocated (cube_index) )                  call aims_deallocate ( cube_index, "cube_index" )
        if ( allocated (cube_fermi_ref) )              call aims_deallocate ( cube_fermi_ref, "cube_fermi_ref" )
        if ( allocated (cube_origin) )                 deallocate (cube_origin )
        if ( allocated (cube_stm) )                    call aims_deallocate ( cube_stm, "cube_stm" )
        if ( allocated (cube_state) )                  call aims_deallocate ( cube_state, "cube_state" )
        if ( allocated (cube_spin) )                   call aims_deallocate ( cube_spin, "cube_spin" )
        if ( allocated (cube_edge_unit) )              deallocate ( cube_edge_unit )
        if ( allocated (cube_edge_steps) )             deallocate ( cube_edge_steps )
        if ( allocated (cube_elf) )                    call aims_deallocate (cube_elf, "cube_elf")
        if ( allocated (cube_filename) )               deallocate (cube_filename)
        if ( allocated (cube_inputname) )              deallocate (cube_inputname)
        if ( allocated (flag_origin) )                 call aims_deallocate ( flag_origin, "flag_origin" )
        if ( allocated (flag_edge) )                   call aims_deallocate ( flag_edge, "flag_edge" )
        if ( allocated (flag_divisor) )                call aims_deallocate (flag_divisor, "flag_divisor")
        if (allocated (cube_format))                   deallocate (cube_format)
        if (allocated (cube_divisor))                  deallocate (cube_divisor)
        if (allocated (cube_type_needs_densmat))       call aims_deallocate (cube_type_needs_densmat, "cube_type_needs_densmat")
        if (allocated (cube_type_needs_eigenvectors))  call aims_deallocate (cube_type_needs_eigenvectors, "cube_type_needs_eigenvectors")
        if (allocated (cube_type_needs_soc_eigenvectors))  call aims_deallocate (cube_type_needs_soc_eigenvectors, "cube_type_needs_soc_eigenvectors")
        if (allocated (cube_average) )                 call aims_deallocate(cube_average, "cube_average")

        end subroutine cleanup_plot
!******	
      end module plot
