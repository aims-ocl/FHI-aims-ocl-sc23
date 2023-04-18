!****h* FHI-aims/thermodynamic_integration
!  NAME
!    thermodynamic_integration - thermodynamic_integration of FHI-aims
!  SYNOPSIS

module thermodynamic_integration 

!  PURPOSE
!    This module takes care of all routines related to the thermodynamic_integration. It uses
!    the interface provided by classical_field.f90 to feed the required potentials to the MD code.
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  AUTHOR
!    Christian Carbogno,  Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Yet unpublished and unreleased
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE
  implicit none

save
!! Arrays containing ALL info for ALL Segments
! Lattice vectors of the equilibrium geometry
real*8,dimension(:,:,:)    ,allocatable :: TDI_lattice_vectors
! Basis vectors of the equilibrium geometry
real*8,dimension(:,:,:)    ,allocatable :: TDI_atoms          
! Force constant for the equilibrium geometry
real*8,dimension(:,:,:,:,:),allocatable :: TDI_force_constants

! Work arrays containing the info for a single segment
real*8,dimension(:,:)    ,allocatable :: TDI_Segment_lattice_vectors
real*8,dimension(:,:)    ,allocatable :: TDI_Segment_atoms          
real*8,dimension(:,:,:,:),allocatable :: TDI_Segment_force_constants

! Work arrays containing... 
! the deviation of the  actual geometry from the equilibrium geometry 
real*8,dimension(:,:)    ,allocatable :: TDI_atoms_DevFromEqGeo          
! the energy of the harmonic system
real*8                                :: TDI_atoms_energy
! the forces acting on the nuclei in the harmonic system
real*8,dimension(:,:)    ,allocatable :: TDI_atoms_forces
! the anharmonic energy, i.e. the DFT free energy 
real*8                                :: TDI_ANH_energy

! Work arrays containing the values read from control.in for the single segments
! Free energy of the harmonic system in equilibrium
real*8                                :: TDI_Segment_energy_offset  
! Starting value for lambda in this MD segment
real*8                                :: TDI_Segment_lambda_start  
! Stopping value for lambda in this MD segment
real*8                                :: TDI_Segment_lambda_end  

! work data
! Actual value of lambda in this time step
real*8                                :: TDI_lambda  
! Actual time derivative of lambda in this segment 
real*8                                :: TDI_dlambda_dt  
! Starting time t=0 for this segment
real*8                                :: TDI_time_offset
! Starting step t=0 for this segment
integer                               :: TDI_step_offset

!Statistics:
! Sum of the anharmonic contributions in this segment
real*8                                :: TDI_sum_anh_cont
! Sum of the anharmonic contributions in this segment (old value) (required for diff.)
real*8                                :: TDI_sum_anh_cont_old


contains

!******	
!------------------------------------------------------------------------------
!****s* thermodynamic_integration/allocate_TDI
!  NAME
!    allocate_TDI
!  SYNOPSIS

subroutine allocate_TDI

!  PURPOSE
!    allocation of all arrays required for the TDI. Safe initial settings are assigned, too.
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE
  implicit none

  if (.not.allocated(TDI_lambda_start    )) allocate(TDI_lambda_start   (MD_segments))
  if (.not.allocated(TDI_lambda_end      )) allocate(TDI_lambda_end     (MD_segments))
  if(use_thermodynamic_integration) then
     if (.not.allocated(TDI_QHA_file        )) allocate(TDI_QHA_file       (MD_segments))
     if (.not.allocated(TDI_QHA_file_exists )) allocate(TDI_QHA_file_exists(MD_segments))
     if (.not.allocated(TDI_QHA_E0          )) allocate(TDI_QHA_E0         (MD_segments))
  endif

  ! initialize to some neutral value
  TDI_lambda_start   (:) = 0d0  
  TDI_lambda_end     (:) = 1d0
  if(use_thermodynamic_integration) then
     TDI_QHA_file       (:) = ''
     TDI_QHA_file_exists(:) = .FALSE.
     TDI_QHA_E0         (:) = 0d0  
  endif

end subroutine allocate_TDI

!******	
!------------------------------------------------------------------------------
!****s* thermodynamic_integration/TDI_allocate_QHA_arrays()
!  NAME
!    TDI_allocate_QHA_arrays()
!  SYNOPSIS
subroutine TDI_allocate_QHA_arrays()
!  PURPOSE
!    allocation of all arrays required for the QHA parameters. Safe initial settings are assigned, too.
!  USES
  use dimensions
  use runtime_choices
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE
  implicit none
!******

  ! Arrays containing ALL info
  if(.not. allocated(TDI_lattice_vectors)) allocate(TDI_lattice_vectors(TDI_segments,3,n_periodic)       )
  if(.not. allocated(TDI_atoms          )) allocate(TDI_atoms          (TDI_segments,3,n_atoms)          )
  if(.not. allocated(TDI_force_constants)) allocate(TDI_force_constants(TDI_segments,3,n_atoms,3,n_atoms))
  ! Arrays containing the info for a single segment
  if(.not. allocated(TDI_Segment_lattice_vectors)) allocate(TDI_Segment_lattice_vectors(3,n_periodic)       )
  if(.not. allocated(TDI_Segment_atoms          )) allocate(TDI_Segment_atoms          (3,n_atoms)          )
  if(.not. allocated(TDI_Segment_force_constants)) allocate(TDI_Segment_force_constants(3,n_atoms,3,n_atoms))
  ! Array for computing the deviances from the eq. geom. and the respective forces
  if(.not. allocated(TDI_atoms_DevFromEqGeo)) allocate(TDI_atoms_DevFromEqGeo(3,n_atoms))
  if(.not. allocated(TDI_atoms_forces      )) allocate(TDI_atoms_forces      (3,n_atoms))

end subroutine TDI_allocate_QHA_arrays
!******

!******	
!------------------------------------------------------------------------------
!****s* thermodynamic_integration/deallocate_TDI
!  NAME
!    deallocate_TDI
!  SYNOPSIS

subroutine deallocate_TDI

!  PURPOSE
!    deallocation of all arrays required for the TDI. 
!  USES
  use dimensions
  use runtime_choices
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE
  implicit none
  if (allocated(TDI_lambda_start    )) deallocate(TDI_lambda_start   )
  if (allocated(TDI_lambda_end      )) deallocate(TDI_lambda_end     )
  if (allocated(TDI_QHA_file        )) deallocate(TDI_QHA_file       )
  if (allocated(TDI_QHA_file_exists )) deallocate(TDI_QHA_file_exists)
  if (allocated(TDI_QHA_E0          )) deallocate(TDI_QHA_E0         )

  if(allocated(TDI_lattice_vectors))         deallocate(TDI_lattice_vectors)
  if(allocated(TDI_atoms          ))         deallocate(TDI_atoms          )
  if(allocated(TDI_force_constants))         deallocate(TDI_force_constants)
  if(allocated(TDI_Segment_lattice_vectors)) deallocate(TDI_Segment_lattice_vectors)
  if(allocated(TDI_Segment_atoms          )) deallocate(TDI_Segment_atoms          )
  if(allocated(TDI_Segment_force_constants)) deallocate(TDI_Segment_force_constants)
  if(allocated(TDI_atoms_DevFromEqGeo))      deallocate(TDI_atoms_DevFromEqGeo)
  if(allocated(TDI_atoms_forces      ))      deallocate(TDI_atoms_forces      )

end subroutine deallocate_TDI
!******

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/read_TDI_QHA_file
!  NAME
!    read_TDI_QHA_file
!  SYNOPSIS
subroutine read_TDI_QHA_file
!  PURPOSE
!     read force constants and respective equilibrium geometries from the respective TDI_QHA_file
!     that is specified in control.in. Such a file is compulsory when using the thermodynamic 
!     integration. The different segments might use different force_constants, though.  
!     Also the integrity of the data, i.e., the "similarity" of the respective equilibrium geometry 
!     and the geometry specified in geometry.in is inspected.
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use localorb_io
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none

  ! local counters
  integer      :: TDI_schedule_step,i_atom

  ! local file i/o related variables 
  character*100 :: info_str
!  integer      :: i_code
!  character*20 :: desc_str

  ! local mean values
  real*8        :: MeanU,StdDevU
  

  if (.not. MD_initialconf_from_restart) then
     ! Give info on what I am doing
     call localorb_info("------------------------------------------------------------",use_unit,'(A)')
     call localorb_info("Parsing the force constant files for the thermodynamic_integration", use_unit,'(2X,A)')
     call localorb_info("",use_unit,'(A)')

     !! Pre-parsing and checking
     ! Iterate over MD_segemnts=TDI_segments
     do TDI_schedule_step = 1, TDI_segments

       ! Pre-Parse file and pre-parse integrity
       write(info_str,'(2X,A,A,A,I4)') "| Pre-parsing file ",trim(TDI_QHA_file(TDI_schedule_step))," for segment ",TDI_schedule_step
       call localorb_info(info_str,use_unit,'(A)')
       call parse_TDI_QHA_file(TDI_QHA_file(TDI_schedule_step))
       
     end do

     !! Allocate arrays for the thermodynamic integration
     call TDI_allocate_QHA_arrays()

     !! Really read the files now
     ! Iterate over MD_segemnts=TDI_segments
     do TDI_schedule_step = 1, TDI_segments

       !! Pre-Parse file and pre-parse integrity
       ! Infos to user
       call localorb_info('  ----------------------------------------------------------', use_unit,'(A)')
       write(info_str,'(2X,A,A,A,I4)') "| Now parsing file ",trim(TDI_QHA_file(TDI_schedule_step))," for segment ",TDI_schedule_step
       call localorb_info(info_str,use_unit,'(A)')
       ! NB: read_TDI_QHA_file writes into the TDI_Segment_* variables to avoid
       !     (de)allocation of multiple Intent(OUT) interfaces
       call read_values_TDI_QHA_file(TDI_QHA_file(TDI_schedule_step))

       call localorb_info('  | --------------------------------------------------------', use_unit,'(A)')
       call localorb_info('  | Equilibrium conditions read from input file:',use_unit,'(A)')
       call TDI_print_EQ_cond()
       call localorb_info('  | --------------------------------------------------------', use_unit,'(A)')
       
       ! Calculate elongations and potential and perform consistency checks:
       call localorb_info("| Calculating the quasi-harmonic potential for this geometry", use_unit,'(2X,A)')
       call update_QH_potential()

       ! Assign value to correct segment
       TDI_lattice_vectors(TDI_schedule_step,:,:)     = TDI_Segment_lattice_vectors(:,:)    
       TDI_atoms          (TDI_schedule_step,:,:)     = TDI_Segment_atoms          (:,:)    
       TDI_force_constants(TDI_schedule_step,:,:,:,:) = TDI_Segment_force_constants(:,:,:,:)
     end do

     ! Setup first run
     ! Due to the particular structure of the MD code, we have to setup the initial values
     ! for the first segment already here
     TDI_Segment_lattice_vectors(:,:)     = TDI_lattice_vectors(1,:,:)     
     TDI_Segment_atoms          (:,:)     = TDI_atoms          (1,:,:)     
     TDI_Segment_force_constants(:,:,:,:) = TDI_force_constants(1,:,:,:,:) 
     TDI_Segment_energy_offset            = TDI_QHA_E0         (1)
     TDI_Segment_lambda_start             = TDI_lambda_start   (1)
     TDI_Segment_lambda_end               = TDI_lambda_end     (1)
     TDI_lambda                           = TDI_lambda_start   (1)
     TDI_dlambda_dt                       = (TDI_Segment_lambda_end - TDI_Segment_lambda_start)/MD_schedule_time(1)
     TDI_time_offset                      = 0d0
     TDI_step_offset                      = 0d0

  else 
     !Only one segment allowed when restarting
     if (MD_segments.ne.1) then
       call localorb_info(" ********************************************************",use_unit,'(A)')
       call localorb_info(" * ERROR:                                                ",use_unit,'(A)')
       call localorb_info(" * When restarting the thermodynamic integration from    ",use_unit,'(A)')
       call localorb_info(" * a restart file, only one MD_segment can be specified. ",use_unit,'(A)')
       call localorb_info(" ********************************************************",use_unit,'(A)')
       call aims_stop
     end if

     !! Allocate arrays for the thermodynamic integration
     call TDI_allocate_QHA_arrays()
     
     !! Read TDI_segment values from file:
     call read_TDI_from_restart()

     call localorb_info('  | --------------------------------------------------------', use_unit,'(A)')
     call localorb_info('  | Equilibrium conditions read from restart file:',use_unit,'(A)')
     call TDI_print_EQ_cond()
     call localorb_info('  | --------------------------------------------------------', use_unit,'(A)')

     !! Assign these values to ALL segments
     !! FIXME At the moment, only one segment is allowed. Is there any need for more?
     do TDI_schedule_step = 1, TDI_segments
       TDI_lattice_vectors(TDI_schedule_step,:,:)     = TDI_Segment_lattice_vectors(:,:)    
       TDI_atoms          (TDI_schedule_step,:,:)     = TDI_Segment_atoms          (:,:)    
       TDI_force_constants(TDI_schedule_step,:,:,:,:) = TDI_Segment_force_constants(:,:,:,:)
     end do

     ! Calculate elongations and potential and perform consistency checks:
     call update_QH_potential()
  end if

  ! Give info to user
  if (.not. use_harmonic_pot_only ) then
    call localorb_info("  | ------------------------------------------------------------",use_unit,'(A)')
    write(info_str,'(2X,A)') "Setting up the thermodynamic integration for MD_segment 1"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A,E30.15,A)') "| Energy offset                 :",TDI_Segment_energy_offset*Hartree," eV"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A,E30.15)')   "| Starting with the lambda value:",TDI_Segment_lambda_start
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A,E30.15)')   "| Stopping   at the lambda value:",TDI_Segment_lambda_end
    call localorb_info(info_str,use_unit,'(A)')
    call localorb_info("  --------------------------------------------------------------",use_unit,'(A)')
  end if

end subroutine read_TDI_QHA_file

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/init_AS
!  NAME
!    init_AS
!  SYNOPSIS
subroutine init_AS
!  PURPOSE
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use localorb_io
!  AUTHOR
!    Christian Carbogno and LMG
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none

  ! local counters
  integer      :: TDI_schedule_step,i_atom

  ! local file i/o related variables 
  character*100 :: info_str
!  integer      :: i_code
!  character*20 :: desc_str

  ! local mean values
  real*8        :: MeanU,StdDevU
  

  if (.not. MD_initialconf_from_restart) then
     ! Setup first run
     ! Due to the particular structure of the MD code, we have to setup the initial values
     ! for the first segment already here
     TDI_Segment_lambda_start             = TDI_lambda_start   (1)
     TDI_Segment_lambda_end               = TDI_lambda_end     (1)
     TDI_lambda                           = TDI_lambda_start   (1)
     TDI_dlambda_dt                       = (TDI_Segment_lambda_end - TDI_Segment_lambda_start)/MD_schedule_time(1)
     TDI_time_offset                      = 0d0
     TDI_step_offset                      = 0d0

  else 
     !Only one segment allowed when restarting
     if (MD_segments.ne.1) then
       call localorb_info(" ********************************************************",use_unit,'(A)')
       call localorb_info(" * ERROR:                                                ",use_unit,'(A)')
       call localorb_info(" * When restarting the thermodynamic integration from    ",use_unit,'(A)')
       call localorb_info(" * a restart file, only one MD_segment can be specified. ",use_unit,'(A)')
       call localorb_info(" ********************************************************",use_unit,'(A)')
       call aims_stop
     end if

     !! Read TDI_segment values from file:
     call read_TDI_from_restart()
  end if

  ! Give info to user
  call localorb_info("  | ------------------------------------------------------------",use_unit,'(A)')
  write(info_str,'(2X,A)') "Setting up the thermodynamic integration for MD_segment 1"
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A,E30.15,A)') "| Energy offset                 :",TDI_Segment_energy_offset*Hartree," eV"
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A,E30.15)')   "| Stopping   at the lambda value:",TDI_Segment_lambda_end
  call localorb_info(info_str,use_unit,'(A)')
  call localorb_info("  --------------------------------------------------------------",use_unit,'(A)')

end subroutine init_AS


!------------------------------------------------------------------------------
!****s* thermodynamic_integration/parse_TDI_QHA_file
!  NAME
!    parse_TDI_QHA_file
!  SYNOPSIS
subroutine parse_TDI_QHA_file(single_TDI_QHA_file)
!  PURPOSE
!     parse_TDI_QHA_file to determine if number of lattice vectors, atoms and number of force constants 
!     is consistent with geometry.in 
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
  character*40,intent(in) :: single_TDI_QHA_file
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  ! local file i/o related variables 
  integer      :: i_code
  logical      :: eof
  character*20 :: desc_str
  
  ! local tag counters
  integer      :: TDI_QHA_lattice_vector
  integer      :: TDI_QHA_atom          
  integer      :: TDI_QHA_force_constant

  ! Open file & check if evereything works // NB: Existence of files already checked while reading control.in
  open (8, FILE=single_TDI_QHA_file)
  read (8,*,iostat=i_code) desc_str
  if (i_code.ne.0) then
    if (myid.eq.0) then
      write(stderr,*) " *** Invalid file ",trim(single_TDI_QHA_file)
      write(stderr,*) " *** specified for the thermodynamic integration."
      write(stderr,*) " *** Aborting."
    end if
    call aims_stop
  end if
  backspace(8)

  !Got until the end of the file and count the relevant tags
  TDI_QHA_lattice_vector = 0
  TDI_QHA_atom           = 0
  TDI_QHA_force_constant = 0
  eof                    = .false.
  desc_str = ""
  do while (.not.eof)
    
    if (desc_str(1:1).eq."#") then
      continue

    else if (desc_str.eq."lattice_vector") then
      TDI_QHA_lattice_vector = TDI_QHA_lattice_vector + 1

    else if (desc_str.eq."atom") then
      TDI_QHA_atom = TDI_QHA_atom + 1

    else if (desc_str.eq."force_constants") then
      TDI_QHA_force_constant = TDI_QHA_force_constant + 1
    end if

    ! Next line
    read (8,*,iostat=i_code) desc_str
    if (i_code.ne.0) then
      eof = .true.
    end if

  end do

  !Close file
  close(8)

  !Check number of lattice vectors
  if (TDI_QHA_lattice_vector .ne. n_periodic) then
     if (myid.eq.0) then
      write(stderr,*) "  *** Number of lattice vectors specified for the thermodynamic integration"
      write(stderr,*) "  *** is not consistent with geometry.in."
      write(stderr,*) "  *** Aborting."
     end if
     call aims_stop
  end if

  !Check number of atoms
  if (TDI_QHA_atom .ne. n_atoms) then
     if (myid.eq.0) then
      write(stderr,*) "  *** Number of atoms specified for the thermodynamic integration"
      write(stderr,*) "  *** is not consistent with geometry.in."
      write(stderr,*) "  *** Aborting."
     end if
     call aims_stop
  end if
 
  ! Force constants
  if (TDI_QHA_force_constant.ne.((n_atoms**2)*n_periodic)) then
     if (myid.eq.0) then
      write(stderr,*) "  *** Number of force_constant LINES specified for the thermodynamic integration"
      write(stderr,*) "  *** is not consistent with the number of atoms specified"
      write(stderr,*) "  *** Aborting."
     end if
     call aims_stop
  end if

end subroutine parse_TDI_QHA_file

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/read_values_TDI_QHA_file
!  NAME
!    read_values_TDI_QHA_file
!  SYNOPSIS
subroutine read_values_TDI_QHA_file(single_TDI_QHA_file)
!  PURPOSE
!     read_values_TDI_QHA_file reads all relevant data from single_TDI_QHA_file
!     and assigns it to the variables 
!     TDI_Segment_lattice_vectors(3,n_periodic)      )
!     TDI_Segment_atoms          (3,n_atoms)         )
!     TDI_Segment_force_constants(3,n_atoms,3,n_atoms)
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use species_data
  use localorb_io
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
  character*40,intent(in) :: single_TDI_QHA_file
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  ! local file i/o related variables 
  integer      :: i_code
  logical      :: eof
  character*20 :: desc_str

  ! local counters
  integer      :: i_atom,i_dim,j_atom,j_dim,i_latt
  
  ! local tag counters
  integer      :: TDI_QHA_lattice_vector
  integer      :: TDI_QHA_atom          
  integer      :: TDI_QHA_force_constant
  ! local tag counters as read from file
  integer      :: TDI_QHA_lattice_vector_fromfile
  integer      :: TDI_QHA_atom_fromfile          
  
  ! local species as read from file for integrity check
  character*20 :: QHA_species

  ! local coordinates counter
  integer   :: i_coord
  integer   :: FC_atom_1,FC_atom_1_coord,FC_atom_2
  integer   :: FC_atom_1_fromfile,FC_atom_1_coord_fromfile,FC_atom_2_fromfile

  ! local temp. fc
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms)  :: TDI_fc_temp
  integer*8     :: n_iter
  character*120 :: info_str

  ! Open file (Existence and usability already checked while parsing)
  open (8, FILE=single_TDI_QHA_file)

  ! Reset counters
  TDI_QHA_lattice_vector = 0
  TDI_QHA_atom           = 0
  TDI_QHA_force_constant = 0
  FC_atom_1       = 1
  FC_atom_1_coord = 1
  FC_atom_2       = 0
  ! Counters in files
  FC_atom_1_fromfile       = 1
  FC_atom_1_coord_fromfile = 1
  FC_atom_2_fromfile       = 0

  !Got until the end of the file and parse its text
  eof                    = .false.
  desc_str = ""
  do while (.not.eof)
    
    if (desc_str(1:1).eq."#") then
      continue

    ! Lattice vectors
    else if (desc_str.eq."lattice_vector") then
      backspace(8)
      TDI_QHA_lattice_vector = TDI_QHA_lattice_vector + 1
      read(8,*) desc_str,(TDI_Segment_lattice_vectors(i_coord,TDI_QHA_lattice_vector), i_coord=1,3,1),& 
                TDI_QHA_lattice_vector_fromfile 

      !Check for integrity
      if (TDI_QHA_lattice_vector .ne. TDI_QHA_lattice_vector_fromfile) then
        if (myid.eq.0) then
          write(stderr,*) "  *** Lattice vectors for the thermodynamic integration"
          write(stderr,*) "  *** are not in ascending order, see specification ",TDI_QHA_lattice_vector,"."
          write(stderr,*) "  *** Aborting."
        end if
        call aims_stop
      end if

    ! Atoms
    else if (desc_str.eq."atom") then
      backspace(8)
      TDI_QHA_atom = TDI_QHA_atom + 1
      read(8,*) desc_str,(TDI_Segment_atoms(i_coord,TDI_QHA_atom), i_coord=1,3,1), QHA_species, TDI_QHA_atom_fromfile 

      !Check for integrity
      if (TDI_QHA_atom .ne. TDI_QHA_atom_fromfile) then
        if (myid.eq.0) then
          write(stderr,*) "  *** Atoms specified for the thermodynamic integration"
          write(stderr,*) "  *** are not in ascending order, see atom line ",TDI_QHA_atom,"."
          write(stderr,*) "  *** Aborting."
        end if
        call aims_stop
      end if
      if (QHA_species .ne. species_name(species(TDI_QHA_atom))) then
        if (myid.eq.0) then
          write(stderr,*) "  *** Species specified for the thermodynamic integration"
          write(stderr,*) "  *** are not consistent with geometry.in, see atom line ",TDI_QHA_atom,"."
          write(stderr,*) "  *** Aborting."
        end if
        call aims_stop
      end if

    ! Force constants: 
    ! FIXME: force_constants MUST come after all lattice_vectro/atom specifications; this is not ensured yet.
    else if (desc_str.eq."force_constants") then

      ! Construct counters, fastest first
      TDI_QHA_force_constant = TDI_QHA_force_constant + 1
      if ( FC_atom_2 .lt. TDI_QHA_atom ) then
        FC_atom_2 = FC_atom_2 + 1
      else
        FC_atom_2 = 1
        if (FC_atom_1_coord .lt. TDI_QHA_lattice_vector) then
          FC_atom_1_coord = FC_atom_1_coord + 1   
        else
          FC_atom_1_coord = 1
          if (FC_atom_1 .lt. TDI_QHA_atom) then
            FC_atom_1 = FC_atom_1 +1
          else
            ! Should not happen due to consistency check before, but check anyway
            if (myid.eq.0) then
              write(stderr,*) "  *** ERROR: Too many force_constants."
              write(stderr,*) "  *** Aborting."
            end if
            call aims_stop
          end if
        end if
      end if

      ! Read input for force constants now
      backspace(8)
      read(8,*) desc_str,(TDI_Segment_force_constants(FC_atom_1_coord,FC_atom_1,i_coord,FC_atom_2),i_coord=1,3,1), &
               & FC_atom_1_fromfile, FC_atom_1_coord_fromfile, FC_atom_2_fromfile 
               

      ! Consistency checks:
      if (FC_atom_1_fromfile .ne. FC_atom_1) then
        if (myid.eq.0) then
          write(stderr,*) "  *** First index in the force constants is not consistent with given specification," 
          write(stderr,*) "  *** i.e.,       (",FC_atom_1_fromfile,"),",FC_atom_1_coord_fromfile,",",FC_atom_2_fromfile
          write(stderr,*) "  *** instead of, (",FC_atom_1,"),",FC_atom_1_coord,",",FC_atom_2,"."
          write(stderr,*) "  *** Are the force constants in the right order (last index fastest)?"
          write(stderr,*) "  *** Aborting."
        end if
        call aims_stop
      end if
      if (FC_atom_1_coord_fromfile .ne. FC_atom_1_coord) then
        if (myid.eq.0) then
          write(stderr,*) "  *** Second index in the force constants is not consistent with given specification," 
          write(stderr,*) "  *** i.e.,       ",FC_atom_1_fromfile,",(",FC_atom_1_coord_fromfile,"),",FC_atom_2_fromfile
          write(stderr,*) "  *** instead of, ",FC_atom_1,",(",FC_atom_1_coord,"),",FC_atom_2,"."
          write(stderr,*) "  *** Are the force constants in the right order (last index fastest)?"
          write(stderr,*) "  *** Aborting."
        end if
        call aims_stop
      end if
      if (FC_atom_2_fromfile .ne. FC_atom_2) then
        if (myid.eq.0) then
          write(stderr,*) "  *** Third index in the force constants is not consistent with given specification," 
          write(stderr,*) "  *** i.e.,       ",FC_atom_1_fromfile,",",FC_atom_1_coord_fromfile,",(",FC_atom_2_fromfile,")"
          write(stderr,*) "  *** instead of, ",FC_atom_1,",",FC_atom_1_coord,",(",FC_atom_2,")."
          write(stderr,*) "  *** Are the force constants in the right order (last index fastest)?"
          write(stderr,*) "  *** Aborting."
        end if
        call aims_stop
      end if
    end if

    ! Next line
    read (8,*,iostat=i_code) desc_str
    if (i_code.ne.0) then
      eof = .true.
    end if
  end do

  ! close file
  close(8)

  !! Get everything in the correct Atomic units
  ! Lattice vectors
  do i_latt=1,TDI_QHA_lattice_vector
    do i_dim=1,3
          TDI_Segment_lattice_vectors(i_dim,i_latt) = TDI_Segment_lattice_vectors(i_dim,i_latt)/bohr
    end do
  end do
  ! Atoms
  do i_atom=1,TDI_QHA_atom
    do i_dim=1,3
          TDI_Segment_atoms(i_dim,i_atom) = TDI_Segment_atoms(i_dim,i_atom)/bohr
    end do
  end do
  ! Force constants
  do i_atom=1,TDI_QHA_atom
    do i_dim=1,3
      do j_atom=1,TDI_QHA_atom
        do j_dim=1,3
          TDI_Segment_force_constants(i_dim,i_atom,j_dim,j_atom) = & 
            TDI_Segment_force_constants(i_dim,i_atom,j_dim,j_atom)*bohr**2/hartree
        end do
      end do
    end do
  end do
  
  ! Check if lattice_vectors are equal to the lattice_vectors defined 
  ! in geometry.in
  do i_latt=1,TDI_QHA_lattice_vector
    do i_dim=1,3
      if ( TDI_Segment_lattice_vectors(i_dim,i_latt) .ne. lattice_vector(i_dim,i_latt) ) then
        if (myid.eq.0) then
          write(stderr,*) " *** Component ",i_dim," of lattice vector",i_latt
          write(stderr,*) " *** is not consistent with geometry.in"
          write(stderr,*) " *** Aborting."
        end if
        call aims_stop
      end if
    end do
  end do

  ! !Close file
  ! close(8)

  ! Explicit symmetrization (translation & permutation) to avoid unit cell drifts:
  write(info_str,'(A)') "  |- Iteratively symmetrizing the force constants: "
  call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
  write(info_str,'(A,A)') "  |  | Iteration |  Max. Dev. in Permutation Symm. |", &
   & " Max. Dev. in Translational Symm. "
  call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
  write(info_str,'(A,I9,A,E21.8,A,E21.8,A)') "  |  | ", 0, " | ", & 
   & max_dev_from_permutation_symmetry(TDI_Segment_force_constants)," (eV/AA^2) | ", & 
   & max_dev_from_momentum_conservation(TDI_Segment_force_constants)," (eV/AA^2)"
  call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
  do n_iter=1,100
    TDI_Segment_force_constants = enforce_permutation_symmetry( TDI_Segment_force_constants )
    TDI_Segment_force_constants = enforce_translational_symmetry( TDI_Segment_force_constants )
    write(info_str,'(A,I9,A,E21.8,A,E21.8,A)') "  |  | ", n_iter, " | ", & 
     & max_dev_from_permutation_symmetry(TDI_Segment_force_constants)," (eV/AA^2) | ", &
     & max_dev_from_momentum_conservation(TDI_Segment_force_constants)," (eV/AA^2)"
    call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
    if ( ( max_dev_from_permutation_symmetry(TDI_Segment_force_constants) .lt. 1.0d-14 ) .and. & 
       & ( max_dev_from_momentum_conservation(TDI_Segment_force_constants) .lt. 1.0d-14 ) ) exit
  end do
  if ( n_iter .ge. 100 ) then
    write(info_str,'(A)') "  |- *** Symmetrization failed."
    call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
    call aims_stop
  end if

end subroutine read_values_TDI_QHA_file

function enforce_permutation_symmetry( fc )
  use dimensions
  implicit none
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms) :: enforce_permutation_symmetry
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms) :: fc
  !local
  integer*8 :: i_atom,j_atom, i_dim, j_dim

  do i_atom=1,n_atoms
    do i_dim=1,3
      do j_atom=1,n_atoms
        do j_dim=1,3
          enforce_permutation_symmetry(i_dim,i_atom,j_dim,j_atom) = &
             0.5d0 * ( fc(i_dim,i_atom,j_dim,j_atom) + fc(j_dim,j_atom,i_dim,i_atom) )
        end do
      end do
    end do
  end do
end function enforce_permutation_symmetry

function enforce_translational_symmetry( fc )
  use dimensions
  implicit none
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms) :: enforce_translational_symmetry
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms) :: fc
  !local
  integer*8 :: i_atom,i_dim, j_dim

  do i_atom=1,n_atoms
    do i_dim=1,3
      do j_dim=1,3
        enforce_translational_symmetry(i_dim,i_atom,j_dim,:) = &
          fc(i_dim,i_atom,j_dim,:) - sum(fc(i_dim,i_atom,j_dim,:))/dble(n_atoms)
      end do
    end do
  end do
end function enforce_translational_symmetry

function max_dev_from_permutation_symmetry( fc )
  use dimensions
  implicit none
  real*8 :: max_dev_from_permutation_symmetry
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms) :: fc
  !local
  real*8 :: dev
  integer*8 :: i_atom,j_atom, i_dim, j_dim

  max_dev_from_permutation_symmetry = 0.0d0
  do i_atom=1,n_atoms
    do i_dim=1,3
      do j_atom=i_atom,n_atoms
        do j_dim=1,3
          dev = abs(fc(i_dim,i_atom,j_dim,j_atom) - fc(j_dim,j_atom,i_dim,i_atom))
          if ( dev .gt. max_dev_from_permutation_symmetry ) max_dev_from_permutation_symmetry = dev
        end do
      end do
    end do
  end do
end function max_dev_from_permutation_symmetry

function max_dev_from_momentum_conservation( fc )
  use dimensions
  implicit none
  real*8 :: max_dev_from_momentum_conservation
  real*8,dimension(1:3,1:n_atoms,1:3,1:n_atoms) :: fc
  !local
  real*8 :: dev
  integer*8 :: i_atom,j_atom, i_dim, j_dim

  max_dev_from_momentum_conservation = 0.0d0
  do i_atom=1,n_atoms
    do i_dim=1,3
      do j_dim=1,3
        dev = abs( sum( fc(i_dim,i_atom,j_dim,:) ) )
        if ( dev .gt. max_dev_from_momentum_conservation ) max_dev_from_momentum_conservation = dev
      end do
    end do
  end do
end function max_dev_from_momentum_conservation

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/TDI_geometry_check
!  NAME
!    TDI_geometry_check
!  SYNOPSIS
subroutine  TDI_geometry_check(MeanValue,StandardDeviation)
!  PURPOSE
!     TDI_geometry_check calculates the mean elongation from equilibrium 
!      and its standard deviation
!     Again, all calculations refer to
!      TDI_Segment_lattice_vectors(3,n_periodic)      
!      TDI_Segment_atoms          (3,n_atoms)         
!      TDI_Segment_force_constants(3,n_atoms,3,n_atoms)
!      TDI_atoms_DevFromEqGeo     (3,n_atoms)
!     and the standard geometry
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use localorb_io
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!  OUTPUT
  real*8, intent(out) :: MeanValue,StandardDeviation
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  
  ! Store the absolute deviations from eq. geo.
  real*8, dimension(:), allocatable :: AbsDevFromEqGeo

  ! local counters
  integer :: i_atom,i_coord 

  ! local messages
  character*100 :: info_str

  ! Allocate 
  if(.not.allocated(AbsDevFromEqGeo)) allocate(AbsDevFromEqGeo(n_atoms))

  ! and ist absolute values
  do i_atom=1,n_atoms
    AbsDevFromEqGeo(i_atom) = sqrt(TDI_atoms_DevFromEqGeo(1,i_atom)**2d0 + &
                                &  TDI_atoms_DevFromEqGeo(2,i_atom)**2d0 + &
                                &  TDI_atoms_DevFromEqGeo(3,i_atom)**2d0)  
  end do
  
  ! Mean Value and Standard deviation
  call TDI_MeanValueAndAvg(AbsDevFromEqGeo,MeanValue,StandardDeviation)

  ! Consistency check
  do i_atom=1,n_atoms
    if((AbsDevFromEqGeo(i_atom).gt.(MeanValue+4d0*StandardDeviation)).or. &
      &(AbsDevFromEqGeo(i_atom).lt.(MeanValue-4d0*StandardDeviation))) then
      write(info_str,'(2X,A,I4,A)') "WARNING: Atom ",i_atom," is very far from equilibrium."
      call localorb_info(info_str, use_unit,'(2X,A)')
    end if
  end do

  ! Deallocate
  deallocate(AbsDevFromEqGeo)

end subroutine  TDI_geometry_check

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/TDI_forces_check
!  NAME
!    TDI_forces_check
!  SYNOPSIS
subroutine  TDI_forces_check(MeanValue,StandardDeviation)
!  PURPOSE
!     TDI_forces_check calculates the mean force and its 
!     standard deviation in the QHA
!     Again, all calculations refer to
!      TDI_Segment_lattice_vectors(3,n_periodic)      
!      TDI_Segment_atoms          (3,n_atoms)         
!      TDI_Segment_force_constants(3,n_atoms,3,n_atoms)
!      TDI_atoms_DevFromEqGeo     (3,n_atoms)
!      TDI_atoms_forces           (3,n_atoms)
!     and the standard geometry
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use localorb_io
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!  OUTPUT
  real*8, intent(out) :: MeanValue,StandardDeviation
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  
  ! Store the absolute deviations 
  real*8, dimension(:), allocatable :: AbsForce

  ! local counters
  integer :: i_atom,i_coord 

  ! local messages
  character*100 :: info_str

  ! Allocate 
  if(.not.allocated(AbsForce)) allocate(AbsForce(n_atoms))

  ! and ist absolute values
  do i_atom=1,n_atoms
    AbsForce(i_atom) = sqrt(TDI_atoms_forces(1,i_atom)**2d0 + &
                         &  TDI_atoms_forces(2,i_atom)**2d0 + &
                         &  TDI_atoms_forces(3,i_atom)**2d0)  
  end do

  ! Mean Value and Standard deviation
  call TDI_MeanValueAndAvg(AbsForce,MeanValue,StandardDeviation)

  ! Consistency check
  do i_atom=1,n_atoms
    if((AbsForce(i_atom).gt.(MeanValue+4d0*StandardDeviation)).or. &
      &(AbsForce(i_atom).lt.(MeanValue-4d0*StandardDeviation))) then
      write(info_str,'(2X,A,I4,A)') "WARNING: Strong forces on atom ",i_atom,", apparently far from equilibrium."
      call localorb_info(info_str, use_unit,'(2X,A)')
    end if
  end do

  ! Deallocate
  deallocate(AbsForce)

end subroutine  TDI_forces_check

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/TDI_elongations
!  NAME
!    TDI_elongations
!  SYNOPSIS
subroutine  TDI_elongations()
!  PURPOSE
!     TDI_elongations calculates the current elongation from equilibrium 
!     Again, all calculations refer to
!      TDI_Segment_lattice_vectors(3,n_periodic)      
!      TDI_Segment_atoms          (3,n_atoms)         
!      TDI_Segment_force_constants(3,n_atoms,3,n_atoms)
!      TDI_atoms_DevFromEqGeo     (3,n_atoms)
!     and the standard geometry
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use synchronize_mpi
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  ! local counters
  integer :: i_atom,i_coord 

  ! Compute elongations
  do i_atom=1,n_atoms

    !if(mod(i_atom-1,n_tasks) /= myid) cycle ! distribute effort over tasks

    do i_coord=1,3
      TDI_atoms_DevFromEqGeo(i_coord,i_atom) = coords(i_coord,i_atom) - TDI_Segment_atoms(i_coord,i_atom)
    end do

  end do
  
  !call sync_matrix(TDI_atoms_DevFromEqGeo,3,n_atoms)

end subroutine  TDI_elongations

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/TDI_compute_force_energy()
!  NAME
!    TDI_compute_forces_and_energy()
!  SYNOPSIS
subroutine  TDI_compute_force_energy()
!  PURPOSE
!     TDI_compute_force_energy() does what the name says :-)
!     Again, all calculations refer to
!      TDI_Segment_lattice_vectors(3,n_periodic)      
!      TDI_Segment_atoms          (3,n_atoms)         
!      TDI_Segment_force_constants(3,n_atoms,3,n_atoms)
!      TDI_atoms_DevFromEqGeo     (3,n_atoms)
!      TDI_atoms_forces           (3,n_atoms)
!      TDI_atoms_energy           
!  USES
  use mpi_tasks
  use dimensions
  use runtime_choices
  use geometry
  use synchronize_mpi
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  ! local counters
  integer :: i_atom,i_coord,j_atom,j_coord 

  ! Set distance from equilibirum
  call TDI_elongations()

  ! Reset values
  TDI_atoms_forces(:,:) = 0d0


  ! Compute forces
  ! First atom
  do i_atom=1,n_atoms

    !if(mod(i_atom-1,n_tasks) /= myid) cycle ! distribute effort over tasks

    ! its cartesian coordinate
    do i_coord=1,3
      ! Sum over all others
      do j_atom=1,n_atoms
        do j_coord=1,3
          !Normal term
          TDI_atoms_forces(i_coord,i_atom) = TDI_atoms_forces(i_coord,i_atom) - &
                  & TDI_Segment_force_constants(i_coord,i_atom,j_coord,j_atom)*TDI_atoms_DevFromEqGeo(j_coord,j_atom)
        end do
      end do
    end do
  end do
 
  ! Sync along MPI tasks
  !call sync_matrix(TDI_atoms_forces,3,n_atoms)

  ! Compute energy of the system
  TDI_atoms_energy = 0d0
  do i_atom=1,n_atoms
    do i_coord=1,3
      TDI_atoms_energy = TDI_atoms_energy - TDI_atoms_forces(i_coord,i_atom)*TDI_atoms_DevFromEqGeo(i_coord,i_atom)
    end do
  end do
  TDI_atoms_energy = 0.5d0*TDI_atoms_energy
  

end subroutine  TDI_compute_force_energy

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/update_QH_potential
!  NAME
!    update_QH_potential
!  SYNOPSIS
subroutine update_QH_potential
!  PURPOSE
!     Calculate elongations, QH potential and forces and add the values to the main energies
!     and forces 
!  USES
  use physics
  use dimensions
  use runtime_choices
  use localorb_io
  use mpi_tasks, only: aims_stop
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  implicit none

  ! local file i/o related variables 
  character*300 :: info_str

  ! Debug statistics
  real*8 :: MeanU,StdDevU
  real*8,dimension(1:3) :: sum_force
  integer :: i_cart

  ! Compute elongations
  call TDI_elongations()
  ! Check geometry now
  call TDI_geometry_check( MeanU, StdDevU )
  write(info_str,'(2X,A,E30.15,A)') "| Mean value              of the elongation from equilibrium  ", MeanU*bohr," A"
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A,E30.15,A)') "| Std. dev.               of the elongation from equilibrium  ", StdDevU*bohr," A"
  call localorb_info(info_str,use_unit,'(A)')

  ! Force and energy integrity checks
  call TDI_compute_force_energy()
  write(info_str,'(2X,A,E30.15,A)') "| Potential energy        in the quasi-harmonic approximation ", & 
       & (TDI_atoms_energy+TDI_Segment_energy_offset)*Hartree," eV"
  call localorb_info(info_str,use_unit,'(A)')

  !Consistency checks for the forces
  call TDI_forces_check( MeanU, StdDevU )
  write(info_str,'(2X,A,E30.15,A)') "| Mean force              in the quasi-harmonic approximation ", MeanU*Hartree/bohr," eV/A"
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A,E30.15,A)') "| Std. Dev. of the forces in the quasi-harmonic approximation ", StdDevU*Hartree/bohr," eV/A"
  call localorb_info(info_str,use_unit,'(A)')
  !! CC: Check for momentum conservation sum_i F_i = 0
  do i_cart=1,3,1
    sum_force(i_cart) = sum( TDI_atoms_forces(i_cart,:) )
  end do
  write(info_str,'(2X,2A,3E30.15,A)') "| Force on the center of mass in the ", &
     "quasi-harmonic approx.   ", sum_force*Hartree/bohr," eV/A"
  call localorb_info(info_str,use_unit,'(A)')
  if ( any( abs( sum_force ) .gt. 1d-8 ) ) then
    write(info_str,'(2X,A)') " *** SEVERE WARNING: Forces on center of mass way too large... Check your Hessian! "
    call localorb_info(info_str,use_unit,'(A)')
    call aims_stop()
  end if


end subroutine update_QH_potential

!------------------------------------------------------------------------------
!****s* thermodynamic_integration/output_TDI_statistics
!  NAME
!    output_TDI_statistics
!  SYNOPSIS
subroutine output_TDI_statistics(E_kin)
!  PURPOSE
!     Write stats to standard output file
!     and forces 
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
  use timing
  use molecular_dynamics
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
! Actual kinetic energy
  implicit none
  real*8,intent(in) :: E_kin
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  !local variables
  character*120 :: info_str

  if(use_thermodynamic_integration) then
     if (use_harmonic_pot_only) then
       write(info_str,'(2X,A)') 'Harmonic potential:'
     else
       write(info_str,'(2X,A)') 'Thermodynamic integration:'
     end if
  else if(use_reffree_AS) then
     write(info_str,'(2X,A)') 'Adiabatic switching:'
  end if
  call localorb_info(info_str,use_unit,'(A)')
  
  ! Statistics with respect to the TDI parameters
  write(info_str,'(2X,A)') 'Results for this single time step: ( with | without energy offset)'
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A)') 'Complete information for previous time-step:'
  write(info_str,'(2X,A,I10)')       "| Time step number                            : ", MD_stepcount-1
  call localorb_info(info_str,use_unit,'(A)')
  if (.not. use_harmonic_pot_only) then
    write(info_str,'(2X,A,E30.15)')    "| Value of lambda                             : ", TDI_lambda
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A,E30.15,A)')  "| Value of dlambda/dt                         : ", TDI_dlambda_dt," 1/ps"
    call localorb_info(info_str,use_unit,'(A)')
  end if

  ! Statistics with respect to the TDI energetics 
  if(use_thermodynamic_integration) then
     write(info_str,'(2X,A,E30.15,A,E30.15,A)')  "| Potential energy of the   harmonic system   : ", & 
          & (TDI_atoms_energy + TDI_Segment_energy_offset)*Hartree," eV | ", TDI_atoms_energy*Hartree," eV"
     call localorb_info(info_str,use_unit,'(A)')
     if (.not. use_harmonic_pot_only) then
       write(info_str,'(2X,A,E30.15,A,E30.15,A)') "| Potential energy of the anharmonic system   : ", &
          & TDI_ANH_energy*Hartree," eV | ", (TDI_ANH_energy - TDI_Segment_energy_offset)*Hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,E30.15,A,E30.15,A)')  & 
            "| Potential energy of the     hybrid system   : ", MD_Epot_last * hartree," eV | ", &
            & (MD_Epot_last-TDI_Segment_energy_offset) * hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,E30.15,A,E30.15,A)')  "| Total     energy of the     hybrid system   : ", & 
            & (E_kin+MD_Epot_last) * hartree," eV | ",(E_kin+MD_Epot_last-TDI_Segment_energy_offset)*Hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
     else 
       write(info_str,'(2X,A,E30.15,A,E30.15,A)')  "| Total     energy of the   harmonic system   : ", & 
            & (E_kin+MD_Epot_last) * hartree," eV | ",(E_kin+MD_Epot_last-TDI_Segment_energy_offset)*Hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
     end if
   
     ! "Interesting" values for TDI:
     if (.not. use_harmonic_pot_only) then
       write(info_str,'(2X,A)') "d/dlambda (Tot. En.  of the hybrid system)"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,E30.15,A)')           "|   Value in this step                        : ", &
            & (TDI_ANH_energy - TDI_atoms_energy - TDI_Segment_energy_offset)*Hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
       ! Statistics with respect to the energetics of the full run:
       TDI_sum_anh_cont_old    = TDI_sum_anh_cont      
       TDI_sum_anh_cont        = TDI_sum_anh_cont        + (TDI_ANH_energy - TDI_atoms_energy - TDI_Segment_energy_offset)
       write(info_str,'(2X,A,E30.15,A)')           "|   Thermod. Average in this segment          : ", & 
           & TDI_sum_anh_cont/real(MD_stepcount-TDI_step_offset)*Hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
       if ( (MD_stepcount-TDI_step_offset) .gt. 1) then
         write(info_str,'(2X,A,E30.15,A)')         "|   Change w.r.t. last step                   : ", & 
             & (TDI_sum_anh_cont/real(MD_stepcount-TDI_step_offset) - & 
                TDI_sum_anh_cont_old/real(MD_stepcount-TDI_step_offset-1))*Hartree," eV"
         call localorb_info(info_str,use_unit,'(A)')
       end if
     end if
  else if (use_reffree_AS) then
     write(info_str,'(2X,A,E30.15,A)') "| Potential energy of the anharmonic system   : ", &
        & TDI_ANH_energy*Hartree," eV | "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A,E30.15,A)')  "| Potential energy of the     scaled system   : ", MD_Epot_last * hartree," eV | "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A,E30.15,A)')  "| Total     energy of the     scaled system   : ", & 
          & (E_kin+MD_Epot_last) * hartree," eV | "
     call localorb_info(info_str,use_unit,'(A)')
   
     ! "Interesting" values for TDI:
     write(info_str,'(2X,A)') "d/dlambda (Tot. En.  of the hybrid system)"
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,'(2X,A,E30.15,A)')           "|   Value in this step                        : ", &
          & (TDI_ANH_energy )*Hartree," eV"
     call localorb_info(info_str,use_unit,'(A)')
     ! Statistics with respect to the energetics of the full run:
     TDI_sum_anh_cont_old    = TDI_sum_anh_cont      
     TDI_sum_anh_cont        = TDI_sum_anh_cont        + (TDI_ANH_energy )
     write(info_str,'(2X,A,E30.15,A)')           "|   Thermod. Average in this segment          : ", & 
         & TDI_sum_anh_cont/real(MD_stepcount-TDI_step_offset)*Hartree," eV"
     call localorb_info(info_str,use_unit,'(A)')
     if ( (MD_stepcount-TDI_step_offset) .gt. 1) then
       write(info_str,'(2X,A,E30.15,A)')         "|   Change w.r.t. last step                   : ", & 
           & (TDI_sum_anh_cont/real(MD_stepcount-TDI_step_offset)-& 
              TDI_sum_anh_cont_old/real(MD_stepcount-TDI_step_offset-1))*Hartree," eV"
       call localorb_info(info_str,use_unit,'(A)')
     end if
 
  end if

  ! Final values computed within this MD segment:
  if ( (tsystem.ge.MD_time).and.(.not.use_harmonic_pot_only) ) then

    write(info_str,'(A)') "  ----------------------------------------------------------"
    call localorb_info(info_str)
    if(use_thermodynamic_integration) then
       write(info_str,'(2X,A)') 'Final values for the thermodynamic integration in this segment: '
    else if(use_reffree_AS) then
       write(info_str,'(2X,A)') 'Final values for the adiabatic switching in this segment: '
    end if
    call localorb_info(info_str)
    write(info_str,'(2X,A,E20.10,A,E20.10,A,E20.10,A)') & 
       & "| Lambda was varied in between ",TDI_Segment_lambda_start," and ",& 
         TDI_Segment_lambda_end," in ",MD_time-TDI_time_offset," ps."
    call localorb_info(info_str)
    write(info_str,'(2X,A,E30.15,A)') & 
       & "| Thermod. Average of  d/dlambda (Tot. En.  of the hybrid system)  :", & 
       & TDI_sum_anh_cont/real(MD_stepcount-TDI_step_offset)*Hartree," eV"
    call localorb_info(info_str)
    TDI_step_offset=MD_stepcount
    TDI_sum_anh_cont        = 0d0
    TDI_sum_anh_cont_old    = 0d0
  end if

  write(info_str,'(A)') "------------------------------------------------------------"
  call localorb_info(info_str)

end subroutine output_TDI_statistics

!****s* thermodynamic_integration/TDI_MeanValueAndAvg
!  NAME
!    TDI_MeanValueAndAvg
!  SYNOPSIS
subroutine  TDI_MeanValueAndAvg(Array,Avg,StdDev)
!  PURPOSE
!     Calculate average and standard deviation for an array of data with the dimension n_atoms
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
  implicit none
!  AUTHOR
!    Christian Carbogno
!  HISTORY
!    Development version, FHI-aims (2010).
!  INPUTS
  real*8,dimension(n_atoms),intent(in)  :: Array
!  OUTPUT
  real*8, intent(out)                   :: Avg,StdDev
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  
  ! local counters
  integer :: i_atom    

  ! Compute mean value
  Avg         = 0d0
  do i_atom=1,n_atoms
    Avg = Avg + Array(i_atom)
  end do
  Avg = Avg/real(n_atoms)

  ! Compute standard deviation
  StdDev = 0d0
  do i_atom=1,n_atoms
    StdDev = StdDev + (Avg - Array(i_atom) )**2d0  
  end do
  if (n_atoms .eq. 1) then
    StdDev = 0d0
  else
    StdDev = sqrt( 1d0/real(n_atoms-1) * StdDev) 
  end if
end subroutine  TDI_MeanValueAndAvg
!******

subroutine read_TDI_from_restart()
!  PURPOSE
!   Reads all parameter of the thermodynamic integration from
!   the restart file
!
!  USES

  use dimensions
  use mpi_tasks, only: aims_stop
  use runtime_choices
  use timing
  use localorb_io, only: localorb_info
  implicit none

!  ARGUMENTS
!  INPUTS
! unit of the aims_MD_restart file
!  OUTPUT
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

  logical  :: file_exists
  character*120 :: info_str

  !Check for existence
  inquire(FILE="aims_TDI_restart.dat",EXIST=file_exists)
  if(.not. file_exists) then
     write(info_str,'(2X,2A)') 'Could not find TDI restart file: aims_TDI_restart.dat'
     call localorb_info(info_str)
     call aims_stop
  end if


  open(file = "aims_TDI_restart.dat", unit = 88, status = 'old', form = 'unformatted', action='read')
  if(use_thermodynamic_integration) then
    ! Work arrays containing the info for a single segment
    read(88) TDI_Segment_lattice_vectors
    read(88) TDI_Segment_atoms
    read(88) TDI_Segment_force_constants
    
    ! Free energy of the harmonic system in equilibrium
    read(88) TDI_Segment_energy_offset
  end if 

  ! Start/end value of lambda
  read(88) TDI_Segment_lambda_start 
  read(88) TDI_Segment_lambda_end   

  ! Actual value of lambda in this time step
  read(88) TDI_lambda
  ! Actual time derivative of lambda in this segment 
  read(88) TDI_dlambda_dt
  ! Starting time t=0 for this segment
  read(88) TDI_time_offset
  ! Starting step t=0 for this segment
  read(88) TDI_step_offset

  !Statistics:
  ! Sum of the anharmonic contributions in this segment
  read(88) TDI_sum_anh_cont
  ! Sum of the anharmonic contributions in this segment squared (old value)
  read(88) TDI_sum_anh_cont_old

  ! Finishing time for this segment
  ! FIXME: Should be in read/write_MD_restart, but is not written there.
  !        To ensure backwards compatibility, it is read/written here
  if (.not.MD_time_restart) then
    read(88) MD_time
  end if

  close(88)

end subroutine read_TDI_from_restart
!******

!****s* thermodynamic_integration/write_TDI_restart
!  NAME
!    TDI_write_MD_restart
!  SYNOPSIS

subroutine write_TDI_restart()

!  PURPOSE
!   Writes all parameter of the thermodynamic integration to 
!   the aims_TDI_restart.dat restart file
!
!  USES

  use dimensions
  use timing    
  use runtime_choices
  use mpi_tasks, only: myid
  use localorb_io, only: localorb_info, use_unit
  implicit none

!  ARGUMENTS
!  INPUTS
!  OUTPUT
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

  ! only write from the zero task ... 
  if (myid.eq.0) then
    open(file = "aims_TDI_restart.dat", unit = 88, status = 'unknown', form = 'unformatted', action='write')
    ! Work arrays containing the info for a single segment
    if(use_thermodynamic_integration) then
       write(88) TDI_Segment_lattice_vectors
       write(88) TDI_Segment_atoms
       write(88) TDI_Segment_force_constants
      
       ! Free energy of the harmonic system in equilibrium
       write(88) TDI_Segment_energy_offset
    end if

    ! Start/end value of lambda
    write(88) TDI_Segment_lambda_start 
    write(88) TDI_Segment_lambda_end   
   
    ! Actual value of lambda in this time step
    write(88) TDI_lambda
    ! Actual time derivative of lambda in this segment 
    write(88) TDI_dlambda_dt
    ! Starting time t=0 for this segment
    write(88) TDI_time_offset
    ! Starting step t=0 for this segment
    write(88) TDI_step_offset
   
    !Statistics:
    ! Sum of the anharmonic contributions in this segment
    write(88) TDI_sum_anh_cont
    ! Sum of the anharmonic contributions in this segment squared (old value)
    write(88) TDI_sum_anh_cont_old

    ! Finishing time for this segment
    ! FIXME: Should be in read/write_MD_restart, but is not written there.
    !        To ensure backwards compatibility, it is read/written here.
    write(88) MD_time

    close(88)
    
    if(use_thermodynamic_integration) then
       call localorb_info('------------------------------------------------------------', use_unit,'(2X,A)')
       call localorb_info('Equilibrium conditions written to restart file:',use_unit,'(2X,A)')
       call TDI_print_EQ_cond()
       call localorb_info('------------------------------------------------------------', use_unit,'(2X,A)')
    end if
  end if

end subroutine write_TDI_restart
!******
subroutine TDI_print_EQ_cond()
!  PURPOSE
!   Inform about the parameters read from
!   the restart file
!
!  USES

  use dimensions
  use timing    
  use runtime_choices
  use localorb_io, only : use_unit
  use mpi_tasks, only: myid
  implicit none

!  ARGUMENTS
!  INPUTS
!  OUTPUT
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

!!! Work arrays containing the info for a single segment
!!  read(iunit) TDI_Segment_lattice_vectors
!!  read(iunit) TDI_Segment_atoms
!!  read(iunit) TDI_Segment_force_constants
!!
!!! Free energy of the harmonic system in equilibrium
!!  read(iunit) TDI_Segment_energy_offset
!!
!!! Actual value of lambda in this time step
!!  read(iunit) TDI_lambda
!!! Actual time derivative of lambda in this segment 
!!  read(iunit) TDI_dlambda_dt
!!! Starting time t=0 for this segment
!!  read(iunit) TDI_time_offset
!!! Starting step t=0 for this segment
!!  read(iunit) TDI_step_offset
!!
!!!Statistics:
!!! Sum of the anharmonic contributions in this segment
!!  read(iunit) TDI_sum_anh_cont
!!! Sum of the anharmonic contributions in this segment (old value)
!!  read(iunit) TDI_sum_anh_cont_old
!!  integer :: iunit
    integer :: i_periodic,i_atom,i_coord

  if(myid.eq.0) then
    if (n_periodic.gt.0) then
       write(use_unit,'(2X,A)') "| Unit cell: "
       do i_periodic = 1, n_periodic, 1
          write(use_unit,'(2X,A,3(2X,F15.6))') "|", &
               (TDI_Segment_lattice_vectors(i_coord,i_periodic)*bohr, &
               i_coord=1,3,1)
       enddo
    else
       write(use_unit,'(2X,A)') "| No unit cell requested."
    end if

    write(use_unit,'(2X,A)') "| Atomic structure: "
    write(use_unit,'(2X,A1,2X,A,10X,A,12X,A,12X,A)') &
         "|","Atom ","x [A]","y [A]","z [A]"
    do i_atom = 1, n_occ_atoms, 1
       write(use_unit,'(2X,A1,I5,3(2X,F15.6))') &
            "|",i_atom, (TDI_Segment_atoms(i_coord,i_atom)*bohr, i_coord=1,3,1)
    enddo
  end if

end subroutine TDI_print_EQ_cond

end module thermodynamic_integration

