!****h* FHI-aims/MD_QH_init
!  NAME
!    Molecular Dynamics: setup positions and velocities
!	according to the quasi-harmonic approximation
!  SYNOPSIS

module MD_QH_init

!  PURPOSE
!    This module takes care of all routines related to the setup of the system according to the
!	quasiharmonic approximation. 
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
real*8, dimension(:,:,:),   allocatable:: MD_QH_latt_vec           ! MD_QH_eq_coord (segment,3,3)
real*8, dimension(:,:,:),   allocatable:: MD_QH_eq_coord           ! MD_QH_eq_coord (segment,3,n_atoms)
real*8, dimension(:,:),     allocatable:: MD_QH_eig_val            ! MD_QH_eig_val  (segment,3*n_atoms)
real*8, dimension(:,:,:,:), allocatable:: MD_QH_eig_vec            ! MD_QH_eig_val  (segment,3*n_atoms,3,n_atoms) "freq.,cartesian,atom"
real*8, dimension(:,:),     allocatable:: MD_QH_rnd_phase          ! MD_QH_rnd_phase(segment,3*n_atoms) "freq."
real*8, dimension(:,:),     allocatable:: MD_QH_rnd_MBamp          ! MD_QH_rnd_MBamp(segment,3*n_atoms) "freq."
integer                                :: MD_QH_n_acoustic = 3     ! Number of modes == 0

! CC: These are of overall switches that determine
! the boundary conditions at the various segments
logical :: MD_QH_allsamephase   = .true.
logical :: MD_QH_allsameMBamp   = .true.
logical :: exploit_localisation = .true.

contains

!******	
!------------------------------------------------------------------------------
!****s* MD_QH_init/allocate_MD_QH
!  NAME
!    allocate_MD_QH
!  SYNOPSIS
subroutine allocate_MD_QH
!  PURPOSE
!    allocation of all arrays required for the initialization of MD_QH_init
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

  if (.not.allocated(MD_QH_first_atom )) allocate(MD_QH_first_atom (MD_QH_init_segments))
  if (.not.allocated(MD_QH_last_atom  )) allocate(MD_QH_last_atom  (MD_QH_init_segments))
  if (.not.allocated(MD_QH_temperature)) allocate(MD_QH_temperature(MD_QH_init_segments))
  if (.not.allocated(MD_QH_file       )) allocate(MD_QH_file       (MD_QH_init_segments))

  ! initialize to some neutral value
  MD_QH_first_atom   (:) = -1   
  MD_QH_last_atom    (:) = -1 
  MD_QH_temperature  (:) = -1.0 
  MD_QH_file         (:) = ''

end subroutine allocate_MD_QH

!******	
!------------------------------------------------------------------------------
!****s* MD_QH_init/deallocate_MD_QH
!  NAME
!    deallocate_MD_QH
!  SYNOPSIS
subroutine deallocate_MD_QH
!  PURPOSE
!    deallocation of all arrays required for the initialization of MD_QH_init
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

  call deallocate_MD_QH_data()
 
  if (allocated(MD_QH_first_atom )) deallocate(MD_QH_first_atom )
  if (allocated(MD_QH_last_atom  )) deallocate(MD_QH_last_atom  )
  if (allocated(MD_QH_temperature)) deallocate(MD_QH_temperature)
  if (allocated(MD_QH_file       )) deallocate(MD_QH_file       )


end subroutine deallocate_MD_QH

!******	
!------------------------------------------------------------------------------
!****s* MD_QH_init/MD_QH_initialize
!  NAME
!    MD_QH_initialize
!  SYNOPSIS
subroutine MD_QH_initialize
!  PURPOSE
!    high level wrapper for MD_QH_init
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
  use geometry
  use species_data
  use molecular_dynamics
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
  character*160 :: info_str
  integer       :: i,i_atom,i_coord
  real*8        :: kineticE

  write(info_str,'(2X,A)') ""
  call localorb_info(info_str)

  write(info_str,'(2X,A)') "-------------------------------------------------------------------------------------------"
  call localorb_info(info_str)

  write(info_str,'(2X,A)') "Initializing atom positions and velocities on the basis of the quasi-harmonic approximation"
  call localorb_info(info_str)

  write(info_str,'(2X,A)') "| Allocating necessary arrays and creating random numbers."
  call localorb_info(info_str)
  call allocate_MD_QH_data()

  write(info_str,'(2X,A)') "| Reading input files containing the quasi-harmonic data:"
  call localorb_info(info_str)

  do i=1,MD_QH_init_segments 
    write(info_str,'(2X,A,A,A)') "| - Pre-parsing file ",trim(MD_QH_file(i)),"."
    call localorb_info(info_str)
    call parse_MD_QH_file(i)
    write(info_str,'(2X,A,A,A)') "| - Reading file ",trim(MD_QH_file(i)),"."
    call localorb_info(info_str)
    call read_MD_QH_file(i)
    write(info_str,'(2X,A)')     "| - Performing the necessary integrity checks"
    call localorb_info(info_str)
    call integrity_MD_QH_file(i)
  end do

  write(info_str,'(2X,A)') "| Initializing positions and velocities"
  call localorb_info(info_str)

  do i=1,MD_QH_init_segments 
    write(info_str,'(2X,A,I3,A,I5,A,I5)') & 
      "| - Segment ",i,": From atom number ",MD_QH_first_atom(i)," to atom number ",MD_QH_last_atom(i)
    call localorb_info(info_str)
    call setup_pos_and_vel_MD_QH(i)
    call calculate_kinetic_energy(v_half,kineticE,MD_QH_first_atom(i),MD_QH_last_atom(i))
    write(info_str,'(2X,A,F20.8,A)') "|   has kinetic energy:            ", &
       kineticE*hartree/dble(MD_QH_last_atom(i)-MD_QH_first_atom(i)+1), &
       " eV / atom "
    call localorb_info(info_str)
  end do
  write(info_str,'(2X,A)') "-------------------------------------------------------------------------------------------"
  call localorb_info(info_str)
  write(info_str,'(2X,A)') "| Cleaning spurious center of mass motion: "
  call localorb_info(info_str)
  call clean_center_of_mass_MD_QH(v_half)
  do i=1,MD_QH_init_segments 
    write(info_str,'(2X,A,I3,A,I5,A,I5)') &
       "| - Segment ",i,": From atom number ", MD_QH_first_atom(i), &
       " to atom number ",MD_QH_last_atom(i)
    call localorb_info(info_str)
    call calculate_kinetic_energy(v_half,kineticE,MD_QH_first_atom(i),MD_QH_last_atom(i))
    write(info_str,'(2X,A,F20.8,A)') &
       "|   has kinetic energy:            ", &
       kineticE * hartree / dble(MD_QH_last_atom(i)-MD_QH_first_atom(i)+1), &
       " eV / atom "
    call localorb_info(info_str)
  end do
  write(info_str,'(2X,A)') "-------------------------------------------------------------------------------------------"
  call localorb_info(info_str)
  write(info_str,'(2X,A)') "| New atomic structure and velocities: "
  call localorb_info(info_str)
  write(info_str,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') "|","Atom ","x [A]","y [A]","z [A]"
  call localorb_info(info_str)
  do i_atom = 1, n_atoms, 1
     write(info_str,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
          "|",i_atom, ": Species ", species_name(species(i_atom)), (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
     call localorb_info(info_str)
     write(info_str,'(2X,A,3(2X,F15.6))') &
          "|         velocity", (v_half(i_coord,i_atom)*bohr, i_coord=1,3,1)
     call localorb_info(info_str)
  end do
  write(info_str,'(2X,A)') "-------------------------------------------------------------------------------------------"
  call localorb_info(info_str)

end subroutine MD_QH_initialize
!******	

subroutine clean_center_of_mass_MD_QH(velocities)
use constants, only: bohr, hartree
use dimensions
use species_data   
use localorb_io
use geometry
use molecular_dynamics
implicit none
  real*8,dimension(1:3,n_atoms),intent(inout) :: velocities


  !local
  real*8                :: E_kin_in,E_kin_out,total_mass
  real*8,dimension(1:3) :: mom_in,mom_out
  integer               :: i_atom, i_cart
  character*200         :: info_str

  call calculate_kinetic_energy(velocities,E_kin_in)

  mom_in(:)  = 0.0d0
  total_mass = 0.0d0
  do i_atom=1,n_atoms
    do i_cart=1,3
      mom_in(i_cart) = mom_in(i_cart) + velocities(i_cart,i_atom)*species_m(species(i_atom))
    end do
    total_mass = total_mass + species_m(species(i_atom))
  end do

  do i_atom=1,n_atoms
    do i_cart=1,3
      velocities(i_cart,i_atom) = velocities(i_cart,i_atom) - mom_in(i_cart) / species_m(species(i_atom)) / dble(n_atoms)
    end do
  end do

  mom_out(:) = 0.d0
  do i_atom=1,n_atoms
    do i_cart=1,3
      mom_out(i_cart) = mom_out(i_cart) + velocities(i_cart,i_atom)*species_m(species(i_atom))
    end do
  end do

  call calculate_kinetic_energy(velocities,E_kin_out)
  do i_atom=1,n_atoms
    do i_cart=1,3
      velocities(i_cart,i_atom) = velocities(i_cart,i_atom) * sqrt(E_kin_in/E_kin_out)
    end do
  end do

  write(info_str,'(2X,A,3(F20.8,X))') "| - Translational momentum before cleaning: ", mom_in*bohr/total_mass
  call localorb_info(info_str)
  write(info_str,'(2X,A,3(F20.8,X))') "| - Translational momentum after  cleaning: ", mom_out*bohr/total_mass
  call localorb_info(info_str)
  write(info_str,'(2X,A,3(F20.8,X))') "| - Kinetic energy scaling factor         : ", sqrt(E_kin_in/E_kin_out)
  call localorb_info(info_str)
  call calculate_kinetic_energy(velocities,E_kin_out)
  write(info_str,'(2X,A,3(F20.8,X))') "| - Total kinetic energy   before cleaning: ", E_kin_in*hartree
  call localorb_info(info_str)
  write(info_str,'(2X,A,3(F20.8,X))') "| - Total kinetic energy   after  cleaning: ", E_kin_out*hartree
  call localorb_info(info_str)

end subroutine

!------------------------------------------------------------------------------
!****s* MD_QH_init/allocate_MD_QH_data
!  NAME
!    allocate_MD_QH_data
!  SYNOPSIS
subroutine allocate_MD_QH_data
!  PURPOSE
!    allocate all arrays required for the data contained in MD_QH_file
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
  use mpi_tasks
  use synchronize_mpi
  use molecular_dynamics
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
  integer :: i_atom,i_segment

  if (.not.allocated(MD_QH_latt_vec )  ) allocate(MD_QH_latt_vec (MD_QH_init_segments,3,3))
  if (.not.allocated(MD_QH_eq_coord )  ) allocate(MD_QH_eq_coord (MD_QH_init_segments,3,n_atoms))
  if (.not.allocated(MD_QH_eig_val  )  ) allocate(MD_QH_eig_val  (MD_QH_init_segments,3*n_atoms))
  if (.not.allocated(MD_QH_eig_vec  )  ) allocate(MD_QH_eig_vec  (MD_QH_init_segments,3*n_atoms,3,n_atoms))
  if (.not.allocated(MD_QH_rnd_phase)  ) allocate(MD_QH_rnd_phase(MD_QH_init_segments,3*n_atoms))
  if (.not.allocated(MD_QH_rnd_MBamp)  ) allocate(MD_QH_rnd_MBamp(MD_QH_init_segments,3*n_atoms))
  
 ! MD_QH_latt_vec (segment,3,3)
 ! MD_QH_eq_coord (segment,3,n_atoms)
 ! MD_QH_eig_val  (segment,3*n_atoms)
 ! MD_QH_eig_vec  (segment,3*n_atoms,3,n_atoms) "freq.,cartesian,atom"
 ! MD_QH_rnd_phase(segment,3*n_atoms) "freq."
 ! MD_QH_rnd_MBamp(segment,3*n_atoms) "freq."

  !Initialize arrays to dummy values
  MD_QH_latt_vec(:,:,:)   = 0.0d0
  MD_QH_eq_coord(:,:,:)   = 0.0d0
  MD_QH_eig_val (:,:)     = 0.0d0
  MD_QH_eig_vec (:,:,:,:) = 0.0d0
  
  !Generate random numbers on task 0:
  call initialize_RNG()
  if (myid.eq.0) then
    !call random_seed
    do i_atom=1,3*n_atoms
       ! call random_number(MD_QH_rnd_MBamp(1,i_atom))
       ! call random_number(MD_QH_rnd_phase(1,i_atom))
       MD_QH_rnd_MBamp(1,i_atom) = random_own()
       MD_QH_rnd_phase(1,i_atom) = random_own()
       do i_segment=2,MD_QH_init_segments
         ! Same random numbers for all segments?
         if (MD_QH_allsamephase) then
            MD_QH_rnd_phase(i_segment,i_atom) = MD_QH_rnd_phase(1,i_atom)
         else
            !call random_number(MD_QH_rnd_phase(i_segment,i_atom))
            MD_QH_rnd_phase(i_segment,i_atom) = random_own()
         end if
         if (MD_QH_allsameMBamp) then
            MD_QH_rnd_MBamp(i_segment,i_atom) = MD_QH_rnd_MBamp(1,i_atom)
         else
            !call random_number(MD_QH_rnd_MBamp(i_segment,i_atom))
            MD_QH_rnd_MBamp(i_segment,i_atom) = random_own()
         end if
       end do
    end do

  end if

end subroutine allocate_MD_QH_data
!******	

!------------------------------------------------------------------------------
!****s* MD_QH_init/deallocate_MD_QH_data
!  NAME
!    deallocate_MD_QH_data
!  SYNOPSIS
subroutine deallocate_MD_QH_data
!  PURPOSE
!    deallocate all arrays required for the data contained in MD_QH_file
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

  if ( allocated(MD_QH_latt_vec) ) deallocate(MD_QH_latt_vec)
  if ( allocated(MD_QH_eq_coord) ) deallocate(MD_QH_eq_coord)
  if ( allocated(MD_QH_eig_val ) ) deallocate(MD_QH_eig_val )
  if ( allocated(MD_QH_eig_vec ) ) deallocate(MD_QH_eig_vec )


end subroutine deallocate_MD_QH_data
!******	

!------------------------------------------------------------------------------
!****s* MD_QH_init/parse_MD_QH_file
!  NAME
!    parse_MD_QH_files
!  SYNOPSIS
subroutine parse_MD_QH_file(i_segment)
!  PURPOSE
!    parse all files specified in control.in for MD_QH_init
!    Just check for correct number of specifications before reading
!  USES
  use dimensions
  use localorb_io
  use mpi_tasks, only: aims_stop, myid
  use runtime_choices
  implicit none
!  ARGUMENTS
!  INPUTS
integer,intent(in) ::  i_segment
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE

  ! local file i/o related variables 
  integer      :: i_code
  logical      :: eof
  character*20 :: desc_str

  ! local counters
  integer :: latt_vec_counter
  integer :: eq_coord_counter 
  integer :: eig_val_counter
  integer :: eig_vec_counter
  integer :: eig_vec_atom_counter

  ! Reset counters
  latt_vec_counter      = 0
  eq_coord_counter      = 0
  eig_val_counter       = 0
  eig_vec_counter       = 0
  
  ! Open file (Existence and usability already checked while parsing)
  open (8, FILE=MD_QH_file(i_segment),status='old',action='read')

  !Go until the end of the file and parse its text
  eof                    = .false.
  desc_str = ""
  do while (.not.eof)

    ! Skip comments
    if (desc_str(1:1).eq."#") then
      continue

    ! Lattice vectors
    else if (desc_str.eq."lattice_vector") then
      latt_vec_counter = latt_vec_counter + 1

    ! Equilibirum positions
    else if (desc_str.eq."atom") then
      eq_coord_counter = eq_coord_counter + 1

    ! Eigenfrequencies 
    else if (desc_str.eq."frequency") then
      eig_val_counter = eig_val_counter + 1

    ! Eigenvectors     
    else if (desc_str.eq."eigenvector") then
      eig_vec_counter = eig_vec_counter + 1

    end if

    ! Next line
    read (8,*,iostat=i_code) desc_str
    if (i_code.ne.0) then
      eof = .true.
    end if
  end do

  ! Close file
  close(8)
  
  !Integrity check
  !Lattice vectors
  if (latt_vec_counter .ne. 3) then
    if (myid.eq.0) then
      write(use_unit,'(1X,A)')         "*** ERROR:"
      write(use_unit,'(1X,A,A,A)')     &
         "*** Number of lattice vectors specified in file ", &
         trim(MD_QH_file(i_segment))," not correct."
    end if
    call aims_stop ()
  end if 
  !Atoms
  if (eq_coord_counter .ne. n_atoms) then
    if (myid.eq.0) then
      write(use_unit,'(1X,A)')         "*** ERROR:"
      write(use_unit,'(1X,A,A,A)')     "*** Number of atoms specified in file ",trim(MD_QH_file(i_segment))," not correct."
    end if
    call aims_stop ()
  end if 
  !Eigenfrequencies
  if (eig_val_counter .ne. 3*n_atoms) then
    if (myid.eq.0) then
      write(use_unit,'(1X,A)')         "*** ERROR:"
      write(use_unit,'(1X,A,A,A)')     &
         "*** Number of eigenfrequencies specified in file ", &
         trim(MD_QH_file(i_segment))," not correct."
    end if
    call aims_stop ()
  end if 
  !Eigenvectors
  if (eig_vec_counter .ne. 3*n_atoms*n_atoms) then
    if (myid.eq.0) then
      write(use_unit,'(1X,A)')         "*** ERROR:"
      write(use_unit,'(1X,A,A,A)')     "*** Number of eigenvectors specified in file ",trim(MD_QH_file(i_segment))," not correct."
    end if
    call aims_stop ()
  end if 

end subroutine parse_MD_QH_file
!******	

!------------------------------------------------------------------------------
!****s* MD_QH_init/read_MD_QH_file
!  NAME
!    read_MD_QH_files
!  SYNOPSIS
subroutine read_MD_QH_file(i_segment)
!  PURPOSE
!    read all files specified in control.in for MD_QH_init
!  USES
  use constants, only: hbar_eV_ps
  use dimensions
  use runtime_choices
  use localorb_io
  use geometry    
  use species_data
  use mpi_tasks, only: aims_stop, myid, n_tasks
  implicit none
!  ARGUMENTS
!  INPUTS
integer,intent(in) ::  i_segment
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE

  ! local file i/o related variables 
  integer      :: i_code
  logical      :: eof
  character*20 :: desc_str
  integer      :: i_coord

  ! Temporary variables
  character*2  :: species_tmp
  integer      :: index_tmp
  integer      :: mode_tmp
  integer      :: atom_tmp

  ! local counters
  integer :: latt_vec_counter
  integer :: eq_coord_counter 
  integer :: eig_val_counter
  integer :: eig_vec_counter
  integer :: eig_vec_atom_counter

  ! Reset counters
  latt_vec_counter      = 0
  eq_coord_counter      = 0
  eig_val_counter       = 0
  eig_vec_counter       = 1
  eig_vec_atom_counter  = 0
  
  ! Open file (Existence and usability already checked while parsing)
  open (8, FILE=MD_QH_file(i_segment),status='old',action='read')

  !Go until the end of the file and parse its text
  eof                    = .false.
  desc_str = ""
  do while (.not.eof)

    ! Skip comments
    if (desc_str(1:1).eq."#") then
      continue

    ! Lattice vectors
    ! MD_QH_latt_vec (segment,3,3)
    else if (desc_str.eq."lattice_vector") then
      backspace(8)
      latt_vec_counter = latt_vec_counter + 1
      read(8,*) desc_str,(MD_QH_latt_vec(i_segment,i_coord,latt_vec_counter), i_coord=1,3,1)

    ! Equilibirum positions
    ! MD_QH_eq_coord (segment,3,n_atoms)
    else if (desc_str.eq."atom") then
      backspace(8)
      eq_coord_counter = eq_coord_counter + 1
      read(8,*) desc_str,(MD_QH_eq_coord(i_segment,i_coord,eq_coord_counter), i_coord=1,3,1),species_tmp,index_tmp

      !Check for correct ordering
      if (index_tmp .ne. eq_coord_counter) then
            if (myid.eq.0) then
              write(use_unit,'(1X,A)')      "*** ERROR:"
              write(use_unit,'(1X,A,A,A)')  "*** Atoms specified in file ",trim(MD_QH_file(i_segment))," are"
              write(use_unit,'(1X,A)')      "*** not in ascending order."
            end if
            call aims_stop ()
      end if
      !Check for correct ordering
      if (species_tmp .ne. species_name(species(eq_coord_counter))) then
            if (myid.eq.0) then
              write(use_unit,'(1X,A)')         "*** ERROR:"
              write(use_unit,'(1X,A,I5,A,A)')  &
                 "*** Species specified for atom ", eq_coord_counter, &
                 " in file ", trim(MD_QH_file(i_segment))
              write(use_unit,'(1X,A)')         "*** is not consistent with geometry.in."
            end if
            call aims_stop ()
      end if

    ! Eigenfrequencies 
    ! MD_QH_eig_val  (segment,3*n_atoms)
    else if (desc_str.eq."frequency") then
      backspace(8)
      eig_val_counter = eig_val_counter + 1
      read(8,*) desc_str,MD_QH_eig_val(i_segment,eig_val_counter)

    ! Eigenvectors     
    ! MD_QH_eig_vec  (segment,3*n_atoms,3,n_atoms) "freq.,cartesian,atom"
    else if (desc_str.eq."eigenvector") then
      backspace(8)
      eig_vec_atom_counter = eig_vec_atom_counter + 1
      if (eig_vec_atom_counter .gt. n_atoms) then 
        eig_vec_counter      = eig_vec_counter + 1
        eig_vec_atom_counter = 1
      end if
      
      read(8,*) desc_str,mode_tmp,atom_tmp,(MD_QH_eig_vec(i_segment,eig_vec_counter,i_coord,eig_vec_atom_counter), i_coord=1,3,1)

      !Check for correct ordering
      if ((mode_tmp .ne. eig_vec_counter).or.(atom_tmp .ne. eig_vec_atom_counter)) then
            if (myid.eq.0) then
              write(use_unit,'(1X,A)')      "*** ERROR:"
              write(use_unit,'(1X,A,I4,A,I4,A,I4,A,I4)') "*** Eigenvector for mode ",mode_tmp," /",eig_vec_counter, &
                                                    & " and atom ",atom_tmp," /",eig_vec_atom_counter 
              write(use_unit,'(1X,A,A,A)')  "*** is not given in the required order in file ",trim(MD_QH_file(i_segment)),"."
            end if
            call aims_stop ()
      end if

    !Number of 0 modes
    else if (desc_str.eq."n_acoustic") then
      backspace(8)
      read(8,*) desc_str,MD_QH_n_acoustic

    end if

    ! Next line
    read (8,*,iostat=i_code) desc_str
    if (i_code.ne.0) then
      eof = .true.
    end if
  end do

  ! Close file
  close(8)

  ! Get all quantities in the correct units:
  ! MD_QH_latt_vec (segment,3,3)
  do latt_vec_counter=1,3
   do i_coord=1,3
     MD_QH_latt_vec(i_segment,i_coord,latt_vec_counter) = MD_QH_latt_vec(i_segment,i_coord,latt_vec_counter) / bohr
   end do
  end do

  ! MD_QH_eq_coord (segment,3,n_atoms)
  do eq_coord_counter=1,n_atoms
   do i_coord=1,3
     MD_QH_eq_coord(i_segment,i_coord,eq_coord_counter) = MD_QH_eq_coord(i_segment,i_coord,eq_coord_counter) / bohr
   end do
  end do

  ! MD_QH_eig_val  (segment,3*n_atoms)
  do eig_val_counter=1,3*n_atoms
    MD_QH_eig_val(i_segment,eig_val_counter) = MD_QH_eig_val(i_segment,eig_val_counter)/1000.0d0/hbar_eV_ps
  end do

end subroutine read_MD_QH_file
!******	

!------------------------------------------------------------------------------
!****s* MD_QH_init/integrity_MD_QH_file
!  NAME
!    integrity_MD_QH_files
!  SYNOPSIS
subroutine integrity_MD_QH_file(i_segment)
!  PURPOSE
!    Integrity check for the data read from MD_QH_file
!  USES
  use dimensions
  use runtime_choices
  use localorb_io
  use geometry
  use mpi_tasks, only: aims_stop, myid
  implicit none
!  ARGUMENTS
!  INPUTS
integer,intent(in) ::  i_segment
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE

 !Local vars
 integer      :: i_coord,j_vector,k_freq
 real*8       :: norm

 ! Lattice vectors
 ! MD_QH_latt_vec (segment,3,3)
 do j_vector=1,3
   norm = 0.0d0
   do i_coord=1,3  
     norm = norm + ( MD_QH_latt_vec(i_segment,i_coord,j_vector)  - lattice_vector(i_coord,j_vector) )**2.0d0
   end do
   if ( sqrt(norm) .gt. 1.0d-8 ) then
     if (myid.eq.0) then
       write(use_unit,'(1X,A)')         "*** ERROR:"
       write(use_unit,'(1X,A,A)')       "*** Lattice vectors specified in file ",trim(MD_QH_file(i_segment))
       write(use_unit,'(1X,A)')         "*** are not consistent with geometry.in."
     end if
     call aims_stop ()
   end if
 end do
 
 ! Eq. positions
 ! MD_QH_eq_coord (segment,3,n_atoms)
 do j_vector=1,n_atoms
   norm = 0.0d0
   do i_coord=1,3  
     norm = norm + ( MD_QH_eq_coord(i_segment,i_coord,j_vector)  - coords(i_coord,j_vector) )**2.0d0
   end do
   if ( sqrt(norm) .gt. 1.0d-8 ) then
     if (myid.eq.0) then
       write(use_unit,'(1X,A)')         "*** ERROR:"
       write(use_unit,'(1X,A,I5,A,A)')  "*** Coordinates for atom ",j_vector," specified in file ",trim(MD_QH_file(i_segment))
       write(use_unit,'(1X,A)')         "*** are not consistent with geometry.in."
     end if
     call aims_stop ()
   end if
 end do

 ! MD_QH_eig_val  (segment,3*n_atoms)
 do k_freq=2,3*n_atoms
   if (( MD_QH_eig_val(i_segment,k_freq) - MD_QH_eig_val(i_segment,k_freq-1) ).lt. -1.0d-11) then
     if (myid.eq.0) then
       write(use_unit,'(1X,A)')         "*** ERROR:"
       write(use_unit,'(1X,A,A)')       "*** Eigenfrequencies specified in file ",trim(MD_QH_file(i_segment))
       write(use_unit,'(1X,A)')         "*** are not in ascending order."
     end if
     call aims_stop ()
   end if
 end do
 do k_freq=1,MD_QH_n_acoustic
   if ( abs(MD_QH_eig_val(i_segment,k_freq)) .gt. 1.0d-2) then
     if (myid.eq.0) then
       write(use_unit,'(1X,A)')         "*** ERROR:"
       write(use_unit,'(1X,A,A)')       "*** The acoustic eigenmodes specified in file ",trim(MD_QH_file(i_segment))
       write(use_unit,'(1X,A)')         "*** do not vanish at Gamma (Soft modes?)."
     end if
     call aims_stop ()
   end if
 end do

 ! MD_QH_eig_vec  (segment,3*n_atoms,3,n_atoms) "freq.,cartesian,atom"
 do k_freq=1,3*n_atoms
   norm = 0.0d0
   do j_vector=1,n_atoms
     do i_coord=1,3
       norm = norm + (MD_QH_eig_vec(i_segment,k_freq,i_coord,j_vector))**2.0d0
     end do
   end do
   if ( abs(sqrt(norm)-1.0d0) .gt. 1.0d-8) then
     if (myid.eq.0) then
       write(use_unit,'(1X,A)')         "*** ERROR:"
       write(use_unit,'(1X,A,I5,A,A)')  "*** Eigenvector for mode ",k_freq," in file ",trim(MD_QH_file(i_segment))
       write(use_unit,'(1X,A)')         "*** is not normalized."
     end if
     call aims_stop ()
   end if  
 end do


end subroutine integrity_MD_QH_file
!******	

!******	
!------------------------------------------------------------------------------
!****s* MD_QH_init/setup_pos_and_vel_MD_QH
!  NAME
!    setup_pos_and_vel_MD_QH
!  SYNOPSIS
subroutine setup_pos_and_vel_MD_QH(i_segment)
!  PURPOSE
!    Integrity check for the data read from MD_QH_file
!  USES
  use constants, only: boltzmann_kb, MD_KE_factor, pi
  use dimensions
  use runtime_choices
  use localorb_io
  use geometry
  use mpi_tasks
  use species_data
  use molecular_dynamics
  use synchronize_mpi
  implicit none
!  ARGUMENTS
!  INPUTS
integer,intent(in) ::  i_segment
!  OUTPUT
!    none
!  AUTHOR
!    Christian Carbogno, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2010).
!  SOURCE

 !Local vars
 integer                    :: i_coord,j_atom,k_freq
 real*8                     :: MB_energy
 real*8,dimension(3*n_atoms):: amplitude
 real*8                     :: loc_factor

 if (myid.eq.0) then

    !Set amplitudes accordingly
    do k_freq=1,MD_QH_n_acoustic
      amplitude(k_freq) = 0.0d0 ! No translation or rotation
    end do
    MB_energy = 0.0d0
    do k_freq=MD_QH_n_acoustic+1,3*n_atoms
      ! Amplitude
      amplitude(k_freq) = sqrt(2.0d0*boltzmann_kB*MD_QH_temperature(i_segment)/MD_KE_factor ) &
                        & / MD_QH_eig_val(i_segment,k_freq) * sqrt(-1.0d0*log(1.0d0 - MD_QH_rnd_MBamp(i_segment,k_freq))) 

      ! Localisation factor for this mode:
      if (exploit_localisation) then
        loc_factor = 0.0d0
        do j_atom=MD_QH_first_atom(i_segment),MD_QH_last_atom(i_segment)
            do i_coord=1,3
               loc_factor = loc_factor + MD_QH_eig_vec(i_segment,k_freq,i_coord,j_atom)**2.0d0
            end do
        end do  
        loc_factor =  loc_factor / ( 3.0d0 * ( MD_QH_last_atom(i_segment) - MD_QH_first_atom(i_segment) + 1)) * (3.0d0*n_atoms)
        MB_energy = MB_energy +  (MD_QH_eig_val(i_segment,k_freq))**2.0d0*(amplitude(k_freq))**2.0d0/4.0d0*loc_factor
      else
        MB_energy = MB_energy +  (MD_QH_eig_val(i_segment,k_freq))**2.0d0*(amplitude(k_freq))**2.0d0/4.0d0
      end if
    end do
    MB_energy   = (3.0d0*n_atoms)*boltzmann_kB*MD_QH_temperature(i_segment)/MD_KE_factor/2.0d0 / MB_energy

    !Generate displacements and velocities
    ! Reset to eq. positions:
    do j_atom=MD_QH_first_atom(i_segment),MD_QH_last_atom(i_segment)
        do i_coord=1,3
          coords(i_coord,j_atom) = MD_QH_eq_coord(i_segment,i_coord,j_atom)
          v_half(i_coord,j_atom) = 0.0d0
        end do
    end do  
  
    ! Iterate over modes
    do k_freq=MD_QH_n_acoustic+1,3*n_atoms
      ! Iterate over atoms
      do j_atom=MD_QH_first_atom(i_segment),MD_QH_last_atom(i_segment)
        ! All components
        do i_coord=1,3
          coords(i_coord,j_atom) = coords(i_coord,j_atom) + &
           & (1.0/sqrt(species_m(species(j_atom)))) * sqrt(MB_energy) * amplitude(k_freq) &
           & * MD_QH_eig_vec(i_segment,k_freq,i_coord,j_atom) * cos(2.0d0*pi*MD_QH_rnd_phase(i_segment,k_freq))
          v_half(i_coord,j_atom) = v_half(i_coord,j_atom) - &
           & (1.0/sqrt(species_m(species(j_atom)))) * sqrt(MB_energy) * amplitude(k_freq) &
           & * MD_QH_eig_vec(i_segment,k_freq,i_coord,j_atom) * sin(2.0d0*pi*MD_QH_rnd_phase(i_segment,k_freq)) &
           & * MD_QH_eig_val(i_segment,k_freq)
        end do
      end do
    end do
 end if
 
 call broadcast_MD_velocities(coords,0)
 call broadcast_MD_velocities(v_half,0)


end subroutine setup_pos_and_vel_MD_QH
!******	

end module MD_QH_init


