! symmetry checker for aims phonon calculations:
! parameters: 
! file with all the forces computed so far 
!     -> this is where we append one line if necessary, then the perl script can simply read it from there ... 
! geometry.in
! control.in
! number of current displacement

! program:
! - sets up supercell including displacement (I)
! - sets up all previous supercells (M)
! - does a number of symmetry operations on (I) and compares them against all (M)
!     o this first of all means checking if the lattice vectors make sense (are integer multiples)
!     o then it means checking if ALL the atoms are the same
! - if symmetry found, program prints forces to standard out from where they can easily be read by the perl parser
!   
program symmetry_check
  use MPI_ROUTINES
  implicit none

  character*100 :: arg, geometry_name, control_name, force_name, keyword, jobtype
  integer :: force_number, ioerr, n_atoms, sc_max, displacement_number
  integer :: n_supercells, n_atoms_supercell, info
  logical :: symmetry_found, check_atom_position_equivalence
  real*8  :: displacement, work(3), symmetry_thresh, buf, force_temp(3)

  ! counter variables:
  integer :: i_coord, i_atom, i_lattice, isc1, isc2, isc3, i_disp, i_species, i_atom_sc
  integer :: index, i_atom2, i_force, i_basis_1, i_basis_2, symmetries_found, i_symm

  ! variables for symmetry setup:
  integer :: perm(3), count, i_sign(3), i_sign_x, i_sign_y, i_sign_z, iswap1, iswap2, ibuf

  integer, parameter :: max_symmetry_count = 48 ! maximal number of transformations checked
  
  real*8, dimension(3,3  )                  :: lattice_vector, lattice_supercell, inverse_lattice_vector
  real*8, dimension(:,:  ),  allocatable    :: coords_in, coords_check, coords_check_2, forces_sum
  real*8, dimension(:,:,:),  allocatable    :: forces, coords
  real*8, dimension(:)    ,  allocatable    :: forces_line
  character*5, dimension(:), allocatable    :: species_in, species_supercell
  integer, dimension(3)                     :: supercell, ipivot
  real*8, dimension(3,3,max_symmetry_count) :: symmetry_transformation, inverse_transformation
  integer, dimension(:),     allocatable    :: correspondence
  integer, dimension(:,:),   allocatable    :: complete_correspondence
  real*8, dimension(:,:),    allocatable    :: offset
  integer, dimension(max_symmetry_count)    :: symmetry_index
  logical, dimension(:),     allocatable    :: atom_treated

  call initialize_mpi ( )

  ! only work on thread 0, the rest are only for the benefit of a number of queuing systems
  if (myid.eq.0) then

     ! evaluate runtime arguments:
     call getarg(1,arg)
     jobtype = trim(arg)
     call getarg(2,arg)
     force_name = trim(arg)
     call getarg(3,arg)
     geometry_name = trim(arg)
     call getarg(4,arg)
     control_name = trim(arg)
     call getarg(5,arg)
     read(arg,*) force_number
     if (jobtype.eq.'DFT_symmetrization') then
        displacement_number = force_number  ! this displacement is to be checked for symmetries
        force_number        = 1             ! this is the number of force points we will have throughout the code
        symmetries_found    = 0
     end if

     ! this is to keep track whether or not we are done yet
     symmetry_found  = .false.
     displacement    = 1d-3  ! default displacement
     symmetry_thresh = 1d-6  ! maximally allowed coordinate difference for atoms to be considered at the same position, keyword driven

     ! parse  control.in: need supercell, and displacement & symmetry threshold if specified
     open(unit = 20, file = control_name, status = "old")
     ioerr = 0
     do while (ioerr.eq.0)
        read(unit = 20, fmt = *, iostat = ioerr) keyword
        if (ioerr.eq.0) then
           if (keyword.eq.'phonon') then
              backspace(20)
              read(20,*) keyword, keyword
              if (keyword.eq.'supercell') then
                 backspace(20)
                 read(20,*) keyword, keyword, supercell(1), supercell(2), supercell(3)
              else if (keyword.eq.'displacement') then
                 backspace(20) 
                 read(20,*) keyword, keyword, displacement
              else if (keyword.eq.'symmetry_thresh') then
                 backspace(20)
                 read(20,*) keyword,keyword,symmetry_thresh
              end if
           end if
        end if
     end do
     close(20)

     ! parse geometry.in: need number of atoms & positions
     ! (1) count atoms
     n_atoms = 0
     open(unit = 20, file = geometry_name, status = "old")
     ioerr = 0
     do while (ioerr.eq.0) 
        read(unit = 20, fmt = *, iostat = ioerr) keyword
        if (ioerr.eq.0) then
           if (keyword.eq.'atom') then
              n_atoms = n_atoms + 1
           end if
        end if
     end do
     close(20)

     ! (2) allocate & read atoms and lattice vectors
     if (.not.allocated(coords_in)) then
        allocate(coords_in(3,n_atoms))
     end if
     if (.not.allocated(species_in)) then
        allocate(species_in(n_atoms))
     end if
     open(unit = 20, file = geometry_name, status = "old")
     ioerr = 0
     i_atom = 0
     i_lattice = 0
     do while (ioerr.eq.0) 
        read(unit = 20, fmt = *, iostat = ioerr) keyword
        if (ioerr.eq.0) then
           if (keyword.eq.'atom') then
              i_atom = i_atom + 1
              backspace(20)
              read(20,*) keyword, coords_in(1,i_atom), coords_in(2,i_atom), coords_in(3,i_atom), species_in(i_atom)
           else if (keyword.eq.'lattice_vector') then
              i_lattice = i_lattice + 1
              backspace(20) 
              read(20,*) keyword, lattice_vector(1,i_lattice), lattice_vector(2,i_lattice), lattice_vector(3,i_lattice)
           end if
        end if
     end do
     close(20)

     ! calculate overall lattice vector for extended supercell - this is to be used for map_to_center_cell only!
     lattice_supercell(:,:) = 0d0
     n_supercells = 1
     do i_coord = 1, 3
        lattice_supercell(:,i_coord) = lattice_vector(:,i_coord)*dble(supercell(i_coord))
        n_supercells = n_supercells * supercell(i_coord)
     end do

     n_atoms_supercell = n_atoms*n_supercells

     if (.not.allocated(species_supercell)) allocate(species_supercell(n_atoms_supercell))

     ! set species in entire extended supercell
     do i_atom = 1, n_atoms_supercell
        i_species = mod(i_atom,n_atoms)
        if (i_species.eq.0) i_species = n_atoms
        species_supercell(i_atom) = species_in(i_species)
     end do

     ! set up all extended supercells
     ! allocate necessary data
     if (.not.allocated(coords))                  allocate(coords(3,n_atoms_supercell,force_number))
     if (.not.allocated(forces))                  allocate(forces(3,n_atoms_supercell,force_number))
     if (.not.allocated(forces_line))             allocate(forces_line(n_atoms*3*supercell(1)*supercell(2)*supercell(3)))
     if (.not.allocated(forces_sum))              allocate(forces_sum(3,n_atoms_supercell))
     if (.not.allocated(coords_check))            allocate(coords_check  (3,n_atoms_supercell))
     if (.not.allocated(coords_check_2))          allocate(coords_check_2(3,n_atoms_supercell))
     if (.not.allocated(correspondence))          allocate(correspondence(n_atoms_supercell))
     if (.not.allocated(offset))                  allocate(offset(3,force_number))
     if (.not.allocated(complete_correspondence)) allocate(complete_correspondence(n_atoms_supercell,max_symmetry_count))
     if (.not.allocated(atom_treated))            allocate(atom_treated(n_atoms_supercell))

     coords(:,:,:)                = 0d0
     forces(:,:,:)                = 0d0
     forces_sum(:,:)              = 0d0
     offset(:,:)                  = 0d0
     complete_correspondence(:,:) = 0
     symmetry_index(:)            = 0
     atom_treated(:)              = .false.

     ! first set all atoms in original position within extended supercell, can be done for all force_numbers at once
     i_atom_sc = 0
     do isc1 = 1, supercell(1)
        do isc2 = 1, supercell(2)
           do isc3 = 1, supercell(3)
              do i_atom = 1, n_atoms
                 i_atom_sc = i_atom_sc + 1
                 do i_coord = 1, 3
                    coords(i_coord,i_atom_sc,:)=coords_in(i_coord,i_atom)+dble(isc1-1)*lattice_vector(i_coord,1)&
                                                                         +dble(isc2-1)*lattice_vector(i_coord,2)&
                                                                         +dble(isc3-1)*lattice_vector(i_coord,3)
                 end do
              end do
           end do
        end do
     end do
     
     ! forces & finite displacements are set in different manners according to jobtype:
     if (jobtype.eq.'new_displacement') then
        ! then go through all the necessary displacements & set atoms accordingly (in central unit cell, starting with atom 1 ... )
        i_force = 0
        do i_atom = 1, n_atoms
           do i_coord = 1, 3
              do i_disp = -1, 1, 2
                 i_force = i_force+1 
                 if (i_force.le.force_number) then
                    offset(:,i_force) = coords(:,i_atom,i_force)      ! remember offset for given force
                    coords(i_coord,i_atom,i_force) = coords(i_coord,i_atom,i_force) + dble(i_disp)*displacement
                 end if
              end do
           end do
        end do
        
        ! read forces from file force_name
        open(unit = 20, file = force_name, status = 'old')
        do i_force = 1, force_number - 1
           forces_line(:) = 0d0
           ! read entire line of forces
           read(20,*) (forces_line(i_coord), i_coord = 1, n_atoms*3*supercell(1)*supercell(2)*supercell(3))
           ! transform into properly arranged force array        
           i_atom_sc = 0
           index     = 0
           do isc1 = 1, supercell(1)
              do isc2 = 1, supercell(2)
                 do isc3 = 1, supercell(3)
                    do i_atom = 1, n_atoms
                       i_atom_sc = i_atom_sc + 1
                       do i_coord = 1, 3
                          index = index + 1
                          forces(i_coord,i_atom_sc,i_force) = forces_line(index)
                       end do
                    end do
                 end do
              end do
           end do
        end do  ! do i_force = 1, force_number - 1
        close(unit = 20)

     else if (jobtype.eq.'DFT_symmetrization') then

        ! update coords: setting the single proper displacement
        i_force = 0
        do i_atom = 1, n_atoms
           do i_coord = 1, 3
              do i_disp = -1, 1, 2
                 i_force = i_force+1 
                 if (i_force.eq.displacement_number) then
                    offset(:,1) = coords(:,i_atom,1)
                    coords(i_coord,i_atom,1) = coords(i_coord,i_atom,1) + dble(i_disp)*displacement
                 end if
              end do
           end do
        end do
 
        ! read forces from input file
        forces_line(:) = 0d0
        open(unit = 30, file = force_name, status = "old")
        read(30,*) (forces_line(i_coord), i_coord = 1, n_atoms*3*supercell(1)*supercell(2)*supercell(3))        
        close(unit = 30)

        ! transform into properly arranged force array        
        index = 0
        i_atom_sc = 0
        do isc1 = 1, supercell(1)
           do isc2 = 1, supercell(2)
              do isc3 = 1, supercell(3)
                 do i_atom = 1, n_atoms
                    i_atom_sc = i_atom_sc + 1
                    do i_coord = 1, 3
                       index = index + 1 
                       forces(i_coord,i_atom_sc,1) = forces_line(index)
                    end do
                 end do
              end do
           end do
        end do
     end if     ! picking forces and finite displacements according to job type

     ! set up inverse lattice matrix from supercell lattice
     inverse_lattice_vector(:,:) = lattice_supercell(:,:)
     call DGETRF(3, 3, inverse_lattice_vector, 3, ipivot, info )
     call DGETRI(3, inverse_lattice_vector, 3, ipivot, work, 3, info)     

     ! calculate all unitary transformation matrices that retain octahedral symmetry, 
     ! this is done by simply permuting all coordinates in all possible ways and with the directions +/-1
     perm(1) = 1
     perm(2) = 2
     perm(3) = 3
     count   = 0
     symmetry_transformation(:,:,:) = 0d0
     inverse_transformation (:,:,:) = 0d0
     do iswap1 = 1, 2
        do iswap2 = 1, 3
           do i_sign_x = 1, -1, -2
              i_sign(1) = i_sign_x
              do i_sign_y = 1, -1, -2
                 i_sign(2) = i_sign_y
                 do i_sign_z = 1, -1, -2
                    i_sign(3) = i_sign_z
                    count = count + 1 
                    do i_coord = 1, 3
                       symmetry_transformation(i_coord,perm(i_coord),count) = dble(i_sign(i_coord))
                       inverse_transformation (perm(i_coord),i_coord,count) = dble(i_sign(i_coord))
                    end do
                 end do
              end do
           end do
           ibuf    = perm(1)
           perm(1) = perm(2)
           perm(2) = perm(3)
           perm(3) = ibuf   
        end do
        ibuf    = perm(1)
        perm(1) = perm(2)
        perm(2) = ibuf
     end do
     
     ! loop through all structures that have been calculated so far 
     i_force = 0
     symmetry_found = .false.
     do while ((.not.symmetry_found).and.((i_force.lt.(force_number-1)).or.((i_force.eq.0).and.(jobtype.eq."DFT_symmetrization"))))
        i_force = i_force + 1 
        count = 0

        ! count enumerates the symmetry operation, loop over all symmetries
        do while ((.not.symmetry_found).and.(count.lt.max_symmetry_count))
           count = count + 1 

! DEBUG
!           write(use_unit,'(I4, A,$)') count, ' | '
           do i_atom = 1, n_atoms_supercell
              coords_check(:,i_atom) = coords(:,i_atom,i_force) - offset(:,i_force)   ! fetch coordinates to be transformed
              call matrix_times_vector(symmetry_transformation(:,:,count),coords_check(:,i_atom))
              coords_check_2(:,i_atom) = coords(:,i_atom,force_number) - offset(:,force_number)
           end do
           
           ! check through all pairs of atoms, see if they correspond, if yes, store number of atom in coords_check 
           ! at the position of atom in coords and remember in index array correspondence
           correspondence(:) = -1
           symmetry_found    = .true.  
           i_atom            = 0           
           do while(symmetry_found.and.(i_atom.lt.n_atoms_supercell))
              i_atom = i_atom + 1 
              do i_atom2 = 1, n_atoms_supercell 
                 if (check_atom_position_equivalence(coords_check(:,i_atom),coords_check_2(:,i_atom2), & ! check position equivalence ...
                      inverse_lattice_vector,symmetry_thresh).and. &
                      (species_supercell(i_atom).eq.species_supercell(i_atom2))) then     ! ... and atom type
                    correspondence(i_atom) = i_atom2
                 end if
              end do
              symmetry_found = symmetry_found .and. (correspondence(i_atom).gt.0)                    
! DEBUG
!              write(use_unit,'(I4,$)') correspondence(i_atom)
           end do
! DEBUG
!           write(use_unit,*)
! DEBUG: output those symmetry transformations that have a 1 on the (1,1) position
!           if (symmetry_transformation(1,1,count).eq.1d0) then
!              do i_coord = 1, 3
!                 write(use_unit,'(2X,A,3F10.6)') '| ',(symmetry_transformation(i_coord,ibuf,count),ibuf = 1,3)
!              end do
!           end if

           ! if there is the desire for a DFT symmetrization, one should transform ALL the forces
           ! and add them to forces_sum
           if ((jobtype.eq."DFT_symmetrization").and.symmetry_found) then
              symmetries_found = symmetries_found + 1
              
              ! remember symmetry & correspondence list 
              complete_correspondence(:,symmetries_found) = correspondence(:)
              symmetry_index(symmetries_found)            = count 
                            
              ! make sure that the loop runs through until the end, done with this particular symmetry
              symmetry_found = .false.
              
! DEBUG: output the symmetry just found
!              do i_coord = 1, 3
!                 write(use_unit,'(2X,A,3F10.6)') '| ',(symmetry_transformation(i_coord,ibuf,count),ibuf = 1,3)
!              end do

           end if
        end do
     end do
     
     ! final status output before launching into the DFT calculation if necessary
     if (symmetry_found.and.(jobtype.eq.'new_displacement')) then

        ! at this point, the applicable symmetry transformation is stored in 
        !     symmetry_transformation(:,:,count)
        !     coords(:,:,i_force)    [no longer needed here]
        !     forces(:,:,i_force)
        !     correspondence(:)      [index array to figure out which atom in the transformed structure 
        !                             corresponds to which atom in the untransformed structure]
        write(use_unit,'(2X,A,I4,A)') 'Found symmetry operation to match previously calculated displacement number', i_force,' : '
        do i_coord = 1, 3
           write(use_unit,'(2X,A,3F10.6)') '| ',(symmetry_transformation(i_coord,ibuf,count),ibuf = 1,3)
        end do
        write(use_unit,'(2X,A)') '| applying unitary transformation to find current forces:'

        ! apply symmetry operation to forces of the previously calculated cell, store in coords_check
        coords_check(:,:) = forces(:,:,i_force)
        do i_atom = 1, n_atoms_supercell      ! transform by symmetry operation count
           call matrix_times_vector(symmetry_transformation(:,:,count),coords_check(:,i_atom))
        end do

        ! set forces(:,:,force_number)
        do i_atom = 1, n_atoms_supercell
           forces(:,correspondence(i_atom),force_number) = coords_check(:,i_atom)
        end do

        ! the same structure used to read the forces into forces_line can be used to write them too - the final result of this program:
        i_atom_sc = 0
        do isc1 = 1, supercell(1)
           do isc2 = 1, supercell(2)
              do isc3 = 1, supercell(3)
                 do i_atom = 1, n_atoms
                    i_atom_sc = i_atom_sc + 1
                    do i_coord = 1, 3
                       write(use_unit,'(F20.15,2X,$)') forces(i_coord,i_atom_sc,force_number) 
                    end do
                 end do
              end do
           end do
        end do
        write(use_unit,*)

     else if (jobtype.eq.'DFT_symmetrization') then
        
        ! output number of supercells:
        if (symmetries_found.gt.1) then
           write(use_unit,'(2X,A,I4,A)') "Found ",symmetries_found," symmetry operations for current configuration and applied those to the DFT forces."
        else
           write(use_unit,'(2X,A)') "Found identity to be the only applicable symmetry transformation for current configuration."
        end if
        write(use_unit,*) "Symmetrized forces: "
        forces_sum(:,:) = 0d0

        ! treat all atoms in groups by adding up all the symmetric forces, adding them, and then returning the forces to their respective 
        ! atoms all the same. 
        do i_atom = 1, n_atoms_supercell
           if (.not.atom_treated(i_atom)) then
              ! add up the force on atom i_atom and all its symmetry equivalents into array work
              work(:) = 0d0
              do i_symm = 1, symmetries_found
                 count   = symmetry_index(i_symm)
                 i_atom2 = complete_correspondence(i_atom,i_symm)
                 force_temp(:) = forces(:,i_atom2,1)
                 call matrix_times_vector(inverse_transformation(:,:,count),force_temp)
                 work(:) = work(:) + force_temp(:)
              end do
              do i_symm = 1, symmetries_found
                 count = symmetry_index(i_symm)
                 i_atom2 = complete_correspondence(i_atom,i_symm)
                 force_temp(:) = work(:)
                 call matrix_times_vector(symmetry_transformation(:,:,count), force_temp)
                 forces_sum(:,i_atom2) = force_temp(:)
                 atom_treated(i_atom2) = .true.
              end do 
           end if
        end do
        

        ! normalize forces by number of symmetries found, there should be at least one (identity)           
        forces_sum(:,:) = forces_sum(:,:) / dble(symmetries_found)
        
        ! DEBUG
        index = 0
        i_atom_sc = 0
        do isc1 = 1, supercell(1)
           do isc2 = 1, supercell(2)
              do isc3 = 1, supercell(3)
                 do i_atom = 1, n_atoms
                    i_atom_sc = i_atom_sc + 1
                    write(use_unit,'(3(F21.16,2X),$)') (forces_sum(i_coord,i_atom_sc), i_coord = 1, 3)

! DEBUG: This eliminates the force symmetrization procedure completely by printing the input data again
!!$                    do i_coord = 1, 3
!!$                       index = index +  1                       
!!$                       write(use_unit,'(F20.15,2X,$)') forces_line(index)
!!$                    end do

                 end do
              end do
           end do
        end do

        write(use_unit,*)
     else         
        write(use_unit,*) 'No symmetry found in previously calculated displacements. '
     end if  !  (symmetry_found)
     
     ! deallocate all data
     if (allocated(coords_in ))              deallocate(coords_in)
     if (allocated(species_in))              deallocate(species_in)
     if (allocated(forces))                  deallocate(forces)
     if (allocated(coords))                  deallocate(coords)
     if (allocated(species_supercell))       deallocate(species_supercell)
     if (allocated(forces_line))             deallocate(forces_line)
     if (allocated(coords_check))            deallocate(coords_check)
     if (allocated(coords_check_2))          deallocate(coords_check_2)
     if (allocated(correspondence))          deallocate(correspondence)
     if (allocated(offset))                  deallocate(offset)
     if (allocated(complete_correspondence)) deallocate(complete_correspondence)
     if (allocated(atom_treated))            deallocate(atom_treated)
  end if ! (myid.eq.0)
  call finalize_mpi ()
end program symmetry_check

! routine to simply multiply a 3x3 matrix with a 3x3 vector and store the result in vector
subroutine matrix_times_vector(matrix,vector)
  implicit none
  real*8  :: matrix(3,3), vector(3), v(3)
  integer :: i_coord, i_coord2
  v(:) = 0d0
  do i_coord = 1, 3
     do i_coord2 = 1, 3
        v(i_coord) = v(i_coord) + matrix(i_coord,i_coord2)*vector(i_coord2)
     end do
  end do
  vector(:) = v(:)
end subroutine matrix_times_vector

! function that determines whether two atoms sit on periodically equal lattice sites
logical function check_atom_position_equivalence(coord1,coord2,inverse_lattice_vector, thresh)
  implicit none
  real*8 :: coord1(3), coord2(3), inverse_lattice_vector(3,3), length(3), thresh
  length(:) = coord1(:) - coord2(:)
  call matrix_times_vector(inverse_lattice_vector,length)
  check_atom_position_equivalence = &
       (dabs(length(1)-dble(nint(length(1)))).lt.thresh).and. &
       (dabs(length(2)-dble(nint(length(2)))).lt.thresh).and. &
       (dabs(length(3)-dble(nint(length(3)))).lt.thresh)
  return
end function check_atom_position_equivalence