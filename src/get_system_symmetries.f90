!****s* FHI-aims/get_system_symmetries
!  NAME
!    get_system_symmetries
!  SYNOPSIS

subroutine get_system_symmetries(n_max_symmetries, n_symmetries, &
&                                symmats, origins)

  !  PURPOSE
  !
  !     Get symmetry matrices which leave the system invariant
  !
  !  USES

  use bravais
  use runtime_choices
  use dimensions
  use geometry
  use localorb_io
  use mpi_tasks, only: aims_stop, check_allocation
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_max_symmetries
  integer, intent(OUT) :: n_symmetries
  real*8, intent(OUT) :: symmats(3, 3, n_max_symmetries)
  real*8, intent(OUT) :: origins(3, n_max_symmetries)

  !  INPUTS
  !    o n_max_symmetries -- Maximum number of possible symmtries (48 is safe)
  !    [global data]
  !  OUTPUTS
  !    o n_symmetries -- Actual number of symmetries
  !    o symmats(:,:, 1:n_symmetries) -- Symmetry matrices
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! maximal number of transformations checked (cubic + hexagonal)
  integer, parameter :: n_trial_transformation = 72
  ! maximal number of symmetry origin transformations (+atoms origins)
  integer, parameter :: max_symmetry_origin = 3
  real*8 :: trial_transformation(3, 3, n_trial_transformation)
  real*8 :: this_trans(3,3), this_origin(3), mapped(3), diff(3)
  real*8, allocatable :: origin(:,:)
  real*8, allocatable :: mapcoord(:,:), refcoord(:,:)
  integer, allocatable :: ref2map(:)
  integer :: n_origin
  integer :: i_coord, i_origin, i_atom, i_atom_ref, i_atom_map
  integer :: perm(3), iswap1, iswap2, ibuf
  integer :: i_sign_x, i_sign_y, i_sign_z, i_sign(3)
  real*8 :: sqrt_three_over_two
  integer :: i_trial_transformation, n_symmetries_uptonow
  logical :: symmetry_found, is_equiv
  character*140 :: info_str
  integer :: info
  character(*), parameter :: func = 'get_system_symmetries'

  call localorb_info("Checking for symmetries", use_unit, '(10X,A)')

  n_origin = max_symmetry_origin + n_atoms
  allocate(origin(3, n_origin), stat=info)
  call check_allocation(info, 'origin', func)
  allocate(mapcoord(3, n_atoms), refcoord(3, n_atoms), stat=info)
  call check_allocation(info, 'mapcoord, refcoord', func)
  allocate(ref2map(n_atoms), stat=info)
  call check_allocation(info, 'ref2map', func)

  ! --- Initialise the rotation origins:

  ! 1: same than origin at geometry.in
  ! 2: average of coordinates
  ! 3: average of max and min values
  ! 3-3+n_atoms: atom coordinates
  ! Here you can add more relevant origin places.
  origin = 0.0d0
  do i_coord = 1, 3
     do i_atom = 1, n_atoms
        origin(i_coord, 2) = origin(i_coord,2) + coords(i_coord,i_atom)/n_atoms
     end do
     origin(i_coord, 3) = 0.5d0*(maxval(coords(i_coord,:)) &
     &                           + minval(coords(i_coord,:)))
  end do
  do i_coord = 1, 3
     do i_atom = 1, n_atoms
        origin(i_coord,3+i_atom) = coords(i_coord,i_atom)
     end do
  end do

  ! --- Initialize the symmetry matrixes:

  ! 3x3 matrixes, which makes symmetry rotations and so on.
  ! Here you can add more symmetry operations.

  ! Calculate all unitary transformation matrices that retain octahedral
  ! symmetry, this is done by simply permuting all coordinates in all possible
  ! ways and with the directions +/-1
  perm(1) = 1
  perm(2) = 2
  perm(3) = 3
  i_trial_transformation = 0
  trial_transformation(:,:,:) = 0d0
  do iswap1 = 1, 2
     do iswap2 = 1, 3
        do i_sign_x = 1, -1, -2
           i_sign(1) = i_sign_x
           do i_sign_y = 1, -1, -2
              i_sign(2) = i_sign_y
              do i_sign_z = 1, -1, -2
                 i_sign(3) = i_sign_z
                 i_trial_transformation = i_trial_transformation + 1
                 do i_coord = 1, 3
                    trial_transformation(i_coord, perm(i_coord), &
                    &                    i_trial_transformation) &
                    & = dble(i_sign(i_coord))
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

  ! --- hexagonal symmetries:

  ! all missing n*60-deg rotations and rotoinversions Currently writing them
  ! such the c-axis is along one of the three coordinate axes, treating the
  ! resulting three orientations separately. Not sure that this is the most
  ! general way and always applicable, but it is pragmatic in the sense that
  ! normally people choose that convention too, so it should cover almost 100%
  ! of all users.
  sqrt_three_over_two = sqrt(3d0)/2d0
  ! (1) with the c-axis of the lattice in the z-direction:
  trial_transformation(3,3,49:52) =   1d0
  trial_transformation(3,3,53:56) =  -1d0
  trial_transformation(1,1,49)    =   0.5d0
  trial_transformation(1,2,49)    =   sqrt_three_over_two
  trial_transformation(2,1,49)    =  -sqrt_three_over_two
  trial_transformation(2,2,49)    =   0.5d0
  trial_transformation(1,1,50)    =  -0.5d0
  trial_transformation(1,2,50)    =  -sqrt_three_over_two
  trial_transformation(2,1,50)    =   sqrt_three_over_two
  trial_transformation(2,2,50)    =  -0.5d0
  trial_transformation(1,1,51)    =  -0.5d0
  trial_transformation(1,2,51)    =   sqrt_three_over_two
  trial_transformation(2,1,51)    =  -sqrt_three_over_two
  trial_transformation(2,2,51)    =  -0.5d0
  trial_transformation(1,1,52)    =   0.5d0
  trial_transformation(1,2,52)    =  -sqrt_three_over_two
  trial_transformation(2,1,52)    =   sqrt_three_over_two
  trial_transformation(2,2,52)    =   0.5d0
  trial_transformation(1:2,1:2,53:56) = trial_transformation(1:2,1:2,49:52)
  ! (2) with the c-axis of the lattice in the x-direction:
  trial_transformation(1,1,57:60) =   1d0
  trial_transformation(1,1,61:64) =  -1d0
  trial_transformation(2,2,57)    =   0.5d0
  trial_transformation(2,3,57)    =   sqrt_three_over_two
  trial_transformation(3,2,57)    =  -sqrt_three_over_two
  trial_transformation(3,3,57)    =   0.5d0
  trial_transformation(2,2,58)    =  -0.5d0
  trial_transformation(2,3,58)    =  -sqrt_three_over_two
  trial_transformation(3,2,58)    =   sqrt_three_over_two
  trial_transformation(3,3,58)    =  -0.5d0
  trial_transformation(2,2,59)    =  -0.5d0
  trial_transformation(2,3,59)    =   sqrt_three_over_two
  trial_transformation(3,2,59)    =  -sqrt_three_over_two
  trial_transformation(3,3,59)    =  -0.5d0
  trial_transformation(2,2,60)    =   0.5d0
  trial_transformation(2,3,60)    =  -sqrt_three_over_two
  trial_transformation(3,2,60)    =   sqrt_three_over_two
  trial_transformation(3,3,60)    =   0.5d0
  trial_transformation(2:3,2:3,61:64) = trial_transformation(2:3,2:3,57:60)
  ! (3) with the c-axis of the lattice in the y-direction:
  trial_transformation(2,2,65:68) =   1d0
  trial_transformation(1,1,65)    =   0.5d0
  trial_transformation(1,3,65)    =   sqrt_three_over_two
  trial_transformation(3,1,65)    =  -sqrt_three_over_two
  trial_transformation(3,3,65)    =   0.5d0
  trial_transformation(1,1,66)    =  -0.5d0
  trial_transformation(1,3,66)    =  -sqrt_three_over_two
  trial_transformation(3,1,66)    =   sqrt_three_over_two
  trial_transformation(3,3,66)    =  -0.5d0
  trial_transformation(1,1,67)    =  -0.5d0
  trial_transformation(1,3,67)    =   sqrt_three_over_two
  trial_transformation(3,1,67)    =  -sqrt_three_over_two
  trial_transformation(3,3,67)    =  -0.5d0
  trial_transformation(1,1,68)    =   0.5d0
  trial_transformation(1,3,68)    =  -sqrt_three_over_two
  trial_transformation(3,1,68)    =   sqrt_three_over_two
  trial_transformation(3,3,68)    =   0.5d0
  trial_transformation(:,:,69:72) =   trial_transformation(:,:,65:68)
  trial_transformation(2,2,69:72) =  -1d0

  ! --- Find out which symmetries the system has

  n_symmetries_uptonow = 0

  do i_trial_transformation = 1, n_trial_transformation
     this_trans = trial_transformation(:,:,i_trial_transformation)
     symmetry_found    = .true.  ! assume the best, plan for the worst

     ! --- Lattice symmetry

     ! Symmetry of lattice vectors == supercell
     ! for the system to have a symmetry, all the transformed reciprocal
     ! lattice vectors have to point to another lattice point, which maps back
     ! onto the origin
     do i_coord = 1, n_periodic
        mapped = matmul(this_trans, lattice_vector(:, i_coord))
        call check_if_lattice_point(n_periodic, lattice_vector, &
        &                           symmetry_thresh, mapped, is_equiv)
        symmetry_found = symmetry_found .and. is_equiv
     end do

     ! If the lattice vectors passed then check on the
     ! Symmetry of atoms in supercell
     if (symmetry_found) then
        ! The k-points do not care of origin, however the atoms do
        ! Test the all symmetry operations with different origins, and check
        ! whether or not the atoms can be transformed on top of one another
        ORIGIN_LOOP: do i_origin = 1, n_origin
           this_origin = origin(:, i_origin)
           symmetry_found = .true.
           do i_atom = 1, n_atoms
              ! fetch coordinates to be transformed
              refcoord(:, i_atom) = coords(:,i_atom) - this_origin
              mapcoord(:, i_atom) = matmul(this_trans, refcoord(:, i_atom))
           end do

           ! check through all pairs of atoms, see if they correspond, if yes,
           ! store number of atom in mapcoord at the position of atom in
           ! coords and remember in index array ref2map
           ref2map(:) = -1
           ATOM_REF_LOOP: do i_atom_ref = 1, n_atoms
              ATOM_MAP_LOOP: do i_atom_map = 1, n_atoms
                 if (species(i_atom_ref)== species(i_atom_map)) then
                    ! check for equivalence of the remapped lattice point ...
                    diff = mapcoord(:, i_atom_map) - refcoord(:, i_atom_ref)
                    call check_if_lattice_point(n_periodic, lattice_vector, &
                    &                           symmetry_thresh, diff, is_equiv)
                    if (is_equiv) then
                       ref2map(i_atom_ref) = i_atom_map
                    end if
                 end if
              end do ATOM_MAP_LOOP
              symmetry_found = symmetry_found .and. (ref2map(i_atom_ref) > 0)
              if (.not. symmetry_found) exit ATOM_REF_LOOP
           end do ATOM_REF_LOOP
           if(symmetry_found)then
              ! Found the symmetry, do not need to test other origins
              exit ORIGIN_LOOP
           end if
        end do ORIGIN_LOOP
     end if

     ! these are the things we need to remember: How many symmetry operations
     ! have been found as well as which ones (and which origin).
     if (symmetry_found) then
        n_symmetries_uptonow = n_symmetries_uptonow + 1
        if (n_symmetries_uptonow > n_max_symmetries) then
           call aims_stop('n_max_symmetries exceeded', func)
        end if
        if (i_origin > n_origin) call aims_stop('Inconsistent i_origin', func)
        symmats(:,:, n_symmetries_uptonow) = this_trans
        origins(:, n_symmetries_uptonow) = this_origin
        write(info_str, "('Found symmetry',3(3F5.2,' '),' around',3F5.2)") &
        & this_trans, this_origin
        call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
     end if
     if (i_trial_transformation == 1 .and. .not. symmetry_found) then
        call aims_stop('Identity missing as symmetry operation', func)
     end if
  end do
  n_symmetries = n_symmetries_uptonow
  write(info_str,'(2X,A,I2,A)') &
  & "Found ", n_symmetries, " symmetry operations for current configuration"
  call localorb_info(info_str,use_unit,'(A)')

  deallocate(origin, mapcoord, refcoord, ref2map)

end subroutine get_system_symmetries
!******
