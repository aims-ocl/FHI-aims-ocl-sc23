!****h* FHI-aims/bravais
!  NAME
!    bravais
!  SYNOPSIS

module bravais

  !  PURPOSE
  !
  !  USES

  implicit none

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

contains

! Moved to geometry.f90
!
!  !----------------------------------------------------------------------------
!  !****s* bravais/min_atomic_dist_old
!  !  NAME
!  !    min_atomic_dist_old
!  !  SYNOPSIS
!  
!  ! THIS VERSION DOES NOT WORK FOR ALL PERIODIC SYSTEMS!!
!  ! TAKE SI, CONVENTIONAL CELL, AND THE SECOND ATOM AT FRAC COORDS 0.3, 0.3, 0.3
!  ! VERSION IS KEPT HERE ONLY FOR ARCHIVAL PURPOSES, FOR NOW!
!
!  subroutine min_atomic_dist_old(lattice_vector, coords, dist, i_atom, j_atom)
!
!    !  PURPOSE
!    !
!    !    Return minimal distance between two atoms.
!    !
!    !  USES
!
!    use dimensions
!    implicit none
!
!    !  ARGUMENTS
!
!    real*8, intent(IN) :: lattice_vector(3, n_periodic)
!    real*8, intent(IN) :: coords(3,n_atoms)
!    real*8, intent(OUT) :: dist
!    integer, intent(OUT), optional :: i_atom, j_atom
!
!    !  INPUTS
!    !   o lattice_vector -- Bravais vectors
!    !   o coords -- Atomic coordinates
!    !  OUTPUTS
!    !   o dist -- Minimal distance between two atoms
!    !   o i_atom, j_atom (optional) -- Atoms with minimal distance.
!    !  AUTHOR
!    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!    !  HISTORY
!    !    Release version, FHI-aims (2011).
!    !  SOURCE
!
!    integer :: i, j, i_periodic
!    real*8 :: thissq, minsq
!    real*8 :: diff_vec(3), length(n_periodic)
!    real*8 :: map_to_center_cell_matrix(n_periodic, 3)
!    character(*), parameter :: func = 'min_atomic_dist'
!
!    ! --- Initialization
!
!    if (present(i_atom)) i_atom = 1
!    if (present(j_atom)) j_atom = 1
!    minsq = huge(minsq)
!
!    if (n_periodic > 0) then
!       ! Cannot use map_to_center_cell() because of custom lattice_vector.
!       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
!       &                                  map_to_center_cell_matrix)
!    end if
!
!    ! --- Distances to own images
!
!    ! JW: For badly chosen Bravais vectors, this would be incorrect.
!    do i_periodic = 1, n_periodic
!       minsq = min(minsq, sum(lattice_vector(:, i_periodic)**2))
!    end do
!
!    ! --- Distances from atoms to centers
!
!    ! THE FOLLOWING IS WHERE THIS GOES WRONG ... but even if we fixed it,
!    ! the loop over unit cells above would make sure we will miss any 
!    ! self-images of atoms in different unit cells.
!
!    do i = 1, n_atoms
!       do j = i+1, n_atoms
!          diff_vec = coords(:, j) - coords(:, i)
!          if (n_periodic > 0) then    ! Map to center cell
!             length = matmul(map_to_center_cell_matrix, diff_vec)
!             length = length - nint(length)
!             diff_vec = matmul(lattice_vector, length)
!          end if
!          thissq = sum(diff_vec**2)
!          if (thissq < minsq) then
!             if (present(i_atom)) i_atom = i
!             if (present(j_atom)) j_atom = j
!             minsq = thissq
!          end if
!       end do
!    end do
!    dist = sqrt(minsq)
!
!  end subroutine min_atomic_dist_old
!  !******
!  !----------------------------------------------------------------------------
!  !****s* bravais/min_atomic_dist
!  !  NAME
!  !    min_atomic_dist
!  !  SYNOPSIS
!  
!  subroutine min_atomic_dist(lattice_vector, coords, dist, out_atom_1, out_atom_2, & 
!                             cell_1, cell_2, cell_3, &
!                             dist_just_atoms, out_just_atoms_1, out_just_atoms_2, &
!                             cell_just_atoms_1, cell_just_atoms_2, cell_just_atoms_3)
!
!    !  PURPOSE
!    !
!    !    Return minimal distance between two atoms / grid centres
!    !
!    !  USES
!
!    use dimensions
!    use localorb_io, only: use_unit
!    ! AJL, Sept2017: Is it a problem if we import these arrays?
!    use species_data, only: species_pseudoized, no_basis
!!    use geometry, only: species, empty
!    implicit none
!
!    !  ARGUMENTS
!
!    real*8, intent(IN) :: lattice_vector(3, n_periodic)
!    real*8, intent(IN) :: coords(3,n_atoms)
!    real*8, intent(OUT) :: dist
!    ! AJL, Sept2017: These labels are deceptive. "atom" here actually means "all grid centres", and it must be 
!    ! remembered that this does not include multipole centres
!    integer, intent(OUT), optional :: out_atom_1, out_atom_2
!    integer, intent(OUT), optional :: cell_1, cell_2, cell_3
!    ! AJL, Sept2017. We need to differentiate between atom distances and all-centre distances
!    ! specifically so we can pick up input errors that aren't to do with geometries.
!    ! I don't want to change the default behaviour, which is to get min distance between all centres,
!    ! so am adding extra functionality to pick out just physically meaningful distances.
!    ! TODO: rename "out_atom_1" etc. to something more representative e.g. out_all_centres_1
!    ! The following is to store distances between specifically atoms and/or pseudocores:
!    real*8, intent(OUT), optional :: dist_just_atoms
!    integer, intent(OUT), optional :: out_just_atoms_1, out_just_atoms_2
!    integer, intent(OUT), optional :: cell_just_atoms_1, cell_just_atoms_2, cell_just_atoms_3
!
!    !  INPUTS
!    !   o lattice_vector -- Bravais vectors
!    !   o coords -- Atomic coordinates
!    !  OUTPUTS
!    !   o dist -- Minimal distance between two grid-centres (so includes empty sites)
!    !   o out_atom_1, out_atom_2 (optional) -- Grid-centres with minimal distance.
!    !   o cell_1, cell_2, cell_3 (optional) -- Unit cell images where minimum-distance grid-centres were found.
!    !   Added Sept2017
!    !   o dist_just_atoms (optional) -- Minimal distance between just physically meaningful centres (atoms and pseudocores)
!    !   o just_atoms_1, just_atoms_2 (optional) --Physically meaningful centres with minimal distance
!    !   o cell_just_atoms_1, cell_just_atoms_2, cell_just_atoms_3 (optional) -- Unit cell images where minimum distance
!    !   physically meaningful centres were found.
!    !  AUTHOR
!    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!    !  HISTORY
!    !    Release version, FHI-aims (2013).
!    !    Edited Sept2017 by AJL to differentiate between atomic centres and all-grid centres 
!    !  SOURCE
!
!    integer :: i_periodic
!    ! AJL, Sept2017: Again, technically "atom" actually means "all grid centres"
!    integer :: i_atom, j_atom
!    integer :: n_shell
!    integer :: n_1, n_2, n_3
!!    real*8 :: minsq, thissq
!    real*8 :: minsq_just_atoms
!
!    real*8 :: previous_minsq
!    real*8 :: previous_minsq_just_atoms
!    real*8 :: length(3)
!
!    integer :: safe_shell_convergence
!
!    real*8 :: diff_vec(3)
!    real*8 :: map_to_center_cell_matrix(n_periodic, 3)
!    real*8, dimension(:,:), allocatable :: coords_temp
!
!    character*150 :: info_str
!
!    character(*), parameter :: func = 'min_atomic_dist'
!
!    ! --- Initialization
!
!    ! Sure, the coords_temp array can be saved altogether in non-periodic
!    ! systems. If this ever causes memory problems, please eliminate coords_temp
!    ! if you can.
!    allocate(coords_temp(3,1:n_atoms))
!
!    ! AJL, Sept2017: TODO: rename these variables to something more representative.
!    if (present(out_atom_1)) out_atom_1 = 1
!    if (present(out_atom_2)) out_atom_2 = 1
!    if (present(cell_1)) cell_1 = 0
!    if (present(cell_2)) cell_2 = 0
!    if (present(cell_3)) cell_3 = 0
!
!    if (present(out_just_atoms_1)) out_just_atoms_1 = 1
!    if (present(out_just_atoms_2)) out_just_atoms_2 = 1
!    if (present(cell_1)) cell_just_atoms_1 = 0
!    if (present(cell_2)) cell_just_atoms_2 = 0
!    if (present(cell_3)) cell_just_atoms_3 = 0
!
!    minsq = huge(minsq)
!    minsq_just_atoms = huge(minsq_just_atoms)
!
!    if (n_periodic > 0) then
!       ! Cannot use map_to_center_cell() because of custom lattice_vector.
!       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
!       &                                  map_to_center_cell_matrix)
!
!       do i_atom = 1, n_atoms, 1
!          length = matmul (map_to_center_cell_matrix,coords(:,i_atom))
!          length = length - nint(length)
!          coords_temp(:,i_atom) = matmul (lattice_vector, length)
!       enddo
!
!    else
!       coords_temp = coords
!    end if
!
!    ! Nice idea but running over arbitrary dimensions systematically
!    ! is best done by a recursive algorithm ... or much more simply
!    ! by copying and simplifying our loops. Please implement if needed,
!    ! not hard to do.
!    if ( (n_periodic.ne.0) .and. (n_periodic.ne.3) ) then
!       write(info_str, '(A)') &
!       & 'Error: This subroutine currently ONLY supports zero-', &
!       & 'or three-dimensional periodic geometries.'
!       call aims_stop(info_str, func)
!    end if
!
!    ! --- Distances to own images
!    ! Initialize minsq. This may not strictly be necessary, but it puts us
!    ! on safe ground.
!    do i_periodic = 1, n_periodic
!       minsq = min(minsq, sum(lattice_vector(:, i_periodic)**2))
!    end do
!
!    ! For non-periodic systems the rest is trivial. Periodic systems are
!    ! more ... confounding.
!
!    ! The primary problem is that unit cell shapes can be truly awkward, with
!    ! very non-orthogonal angles between them. Thus, even in a primitive Bravais
!    ! lattice one can easily construct unit cells where the shortest connecting
!    ! lattice vector is given by (3*a_1 - 4*a_2) or similar.
!    ! Since lattice vectors are integer linear combinations of the 
!    ! unit cell vectors specified in geometry.in, many clever strategies to find
!    ! the shortest linear combinations anaytically will fail because they implicitly
!    ! rely on continuous coordinates somehow.
!    ! Thus, the best algorithm that I could think of is the trivial one - just 
!    ! enumerate unit cells in shell after shell until we are safely converged.
!
!    ! The second problem is that the atoms in the unit cell can be located 
!    ! anywhere they like. We can map them back into the Wigner-Seitz cell, 
!    ! but that will still leave corner cases for unpleasantly chosen unit 
!    ! cell geometries. Thus, we just iterate over all those atoms too.
!
!    ! I know we could clean the math for elegance. If it becomes necessary, please do.
!    ! Sometimes, however, it is better to have a safe solution, even if it is not
!    ! very aesthetic.
!
!    ! Comment recorded after the fact:
!    !
!    ! What would be even better is to first cycle only over possible lattice vectors,
!    ! shell by shell, as done below.
!    !
!    ! From among those, find the triple of lattice vectors R_1, R_2, R_3 that 
!    !
!    ! (a) yields the volume of only the primitive cell (not 2 or more, not zero)
!    !
!    ! (b) has the smallest sum of squares (R_1)^2 + (R_2)^2 + (R_3)^2 .
!    !
!    ! This would have to be the "most orthogonal" triple of lattice vectors.
!    !
!    ! This triple of lattice vectors would be the one to use for real-space
!    ! lists of atoms with roughly equal extent in all directions.
!    !
!    ! The above algorithm is not implemented, but if, for any reason, we see
!    ! sustained trouble with inconvenient unit cell choices in the future,
!    ! we should consider it.
!
!    n_shell = 0
!    safe_shell_convergence = 2 ! This means that we iterate over an extra
!                               ! shell of unit cells before we declare convergence.
!                               ! If this becomes a time problem for large structures
!                               ! there are certainly better ways. I (VB) just do not see
!                               ! that yet.
!    previous_minsq = huge(previous_minsq)
!    ! AJL, Sept2017 Not sure if we need to store information specific to just atoms, but adding it anyway to be safe.
!    previous_minsq_just_atoms = huge(previous_minsq_just_atoms)
!    do while (safe_shell_convergence.gt.0)
!
!       ! run over shells of unit cells
!       do n_1 = -n_shell, n_shell, 1
!          do n_2 = -n_shell, n_shell, 1
!             do n_3 = -n_shell, n_shell, 1
!
!                if ( (n_1.ne.n_shell) .and. (n_2.ne.n_shell) .and. (n_3.ne.n_shell) ) then
!                   ! we were here before - skip this combination of unit cells
!                   continue
!                else
!
!                   do i_atom = 1, n_atoms, 1
!
!                      if (n_periodic.ne.0) then
!                         length(:) = coords_temp(:,i_atom) + n_1 * lattice_vector(:,1) & 
!                                                           + n_2 * lattice_vector(:,2) & 
!                                                           + n_3 * lattice_vector(:,3) 
!                      else
!                         length(:) = coords_temp(:,i_atom)
!                      end if
!
!                      ! Can NOT do just the upper triangle here, in case of periodic images.
!                      do j_atom = 1, n_atoms, 1
!
!                         ! Evaluate all this EXCEPT if we are in the zeroth unit cell,
!                         ! where we must skip any self-images
!                         if (.not.( (i_atom.eq.j_atom) .and. & 
!!                                    (n_1.eq.0) .and. (n_2.eq.0) .and. (n_3.eq.0) )) then
!                            diff_vec = coords_temp(:, j_atom) - length(:)
!                            thissq = sum(diff_vec**2)
!                            if (thissq < minsq) then
!                               if (present(out_atom_1)) out_atom_1 = i_atom
!                               if (present(out_atom_2)) out_atom_2 = j_atom
!                               if (present(cell_1)) cell_1 = n_1
!                               if (present(cell_2)) cell_2 = n_2
!                               if (present(cell_3)) cell_3 = n_3
!                               minsq = thissq
!                            end if
!                            ! AJL/Sept2017
!                            ! Store here information only if the centres are physically meaningful
!                            ! Can this be simplified? I can't remember the difference betwenn empty and no_basis.
!                            ! For now, as is lifted from close_encounters in pbc_lists.f90
!!                            if ( (.not.(empty(i_atom).and.species_pseudoized(species(j_atom)))) .and. &
!!                                 (.not.(species_pseudoized(species(i_atom)).and.empty(j_atom))) .and. &
!!                                 (.not.(no_basis(species(i_atom)).and.(no_basis(species(j_atom))))) .and. & !) then
!!                                 (thissq < minsq_just_atoms) ) then
!                            if (thissq < minsq_just_atoms) then
!                               if (present(out_just_atoms_1)) out_just_atoms_1 = i_atom
!                               if (present(out_just_atoms_2)) out_just_atoms_2 = j_atom
!                               if (present(cell_just_atoms_1)) cell_just_atoms_1 = n_1
!                               if (present(cell_just_atoms_2)) cell_just_atoms_2 = n_2
!                               if (present(cell_just_atoms_3)) cell_just_atoms_3 = n_3
!                               minsq_just_atoms = thissq
!                            end if
!                         end if
!
!                      enddo
!                   enddo
!
!                end if
!
!             enddo ! n_3
!          enddo ! n_2
!       enddo ! n_1
!
!       if (n_periodic.eq.0) then
!          ! we do not need to run over any shells, just state that we are done.
!          safe_shell_convergence = 0
!       else
!         ! check shell convergence
!         if (minsq .lt. previous_minsq) then
!            ! we found something new - do another shell.
!            previous_minsq = minsq
!            previous_minsq_just_atoms = minsq_just_atoms
!            safe_shell_convergence = 2
!         else
!            safe_shell_convergence = safe_shell_convergence - 1
!         end if
!         previous_minsq = minsq
!         previous_minsq_just_atoms = minsq_just_atoms
!       end if
!
!        n_shell = n_shell + 1
!    enddo
!
!    ! Reset to last shell used, even if only needed for test output.
!    n_shell = n_shell - 1
!
!    ! now we are done - here is the minimum distance found ...
!    dist = sqrt(minsq)
!    ! AJL, Sept2017: And just for atomic (physical) centres
!    dist_just_atoms = sqrt(minsq_just_atoms)
!
!    deallocate(coords_temp)
!
!    ! test output only
!    !if (myid.eq.0) then
!    !
!    !  write (use_unit, *) " | Minimum distance found: ", dist*bohr
!    !  write (use_unit, *) " | Shells used:            ", n_shell - 1
!    !  write (use_unit, *) " | Closest atoms:          ", out_atom_1, out_atom_2
!    !  write (use_unit, *) " | Unit cell of atom 1:    ", cell_1, cell_2, cell_3
!    !
!    !end if
!    ! end test
!
!  end subroutine min_atomic_dist
!  !******
!  !----------------------------------------------------------------------------
!  !****s* bravais/min_multipole_distance
!  !  NAME
!  !    min_multipole_distance
!  !  SYNOPSIS
!  
!  subroutine min_multipole_dist & 
!  & ( lattice_vector, occ_coords, species, atom_radius_sq, multipole_coords, & 
!  &   empty_coords, dist, min_atom, min_multipole, offending_multipole )
! 
!    !  PURPOSE
!    !
!    !    Determine distance between each external embedding multipole and closest atom.
!    !    If too close to the nearest atom, report back to calling subroutine.
!    !
!    !  USES
!
!    use dimensions
!    implicit none
!
!    !  ARGUMENTS
!
!    real*8, intent(IN) :: lattice_vector(3, n_periodic)
!    real*8, intent(IN) :: occ_coords(3,n_occ_atoms)
!    integer, intent(IN) :: species(n_atoms)
!    real*8, intent(IN) :: atom_radius_sq(n_species)
!    real*8, intent(IN) :: multipole_coords(3,n_multipoles)
!    real*8, intent(IN) :: empty_coords(3,n_empty_atoms)
!    real*8, intent(OUT) :: dist
!    integer, intent(OUT) :: min_atom, min_multipole
!    integer, dimension(n_multipoles), intent(OUT) :: offending_multipole
!
!    !  INPUTS
!    !   o lattice_vector -- Bravais vectors
!    !   o occ_coords -- Occupied atomic coordinates
!    !   o species -- Species number of each atom
!    !   o atom_radius_sq -- square of the radius of the most extended basis function of each species
!    !   o multipole_coords -- Multipole coordinates
!    !   o empty_coords -- Empty site coordinates
!    !  OUTPUTS
!    !   o dist -- Minimal distance between two atoms
!    !   o i_atom, i_multipole -- Atom and multipole with minimal distance.
!    !   o offending_multipole -- List (number in order of appearance in geometry.in) of all
!    !                            multipoles that are too close to an atom.
!    !  AUTHOR
!    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!    !  HISTORY
!    !    Release version, FHI-aims (2011).
!!    !  SOURCE
!
!    integer :: i_multipole, i_atom, i_periodic
!    integer :: n_offending_multipoles
!    real*8 :: thissq, minsq, thissq_empty, minsq_empty
!    real*8 :: diff_vec(3), length(n_periodic)
!    real*8 :: map_to_center_cell_matrix(n_periodic, 3)
!
!    real*8 :: safety_margin  ! hardwired here for now but should become configurable later.
!                     
!    character(*), parameter :: func = 'min_atomic_dist'
!
!    ! --- Initialization
!
!    min_atom = 0
!    min_multipole = 0
!    n_offending_multipoles = 1       ! real value is one less 
!    offending_multipole = 0
!    minsq = huge(minsq)
!
!    safety_margin = 1.0  ! Value is chosen in units of bohr radii .
!
!    if (n_periodic > 0) then
!       ! Cannot use map_to_center_cell() because of custom lattice_vector.
!!       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
!       &                                  map_to_center_cell_matrix)
!    end if
!
!    ! --- Distances to own images
!
!    ! What we wish to find 
!
!
!    do i_periodic = 1, n_periodic
!       minsq = min(minsq, sum(lattice_vector(:, i_periodic)**2))
!    end do
!
!    ! --- Distances from multipoles to centers
!
!    do i_multipole = 1, n_multipoles
!  
!       ! these need resetting every loop
!       minsq_empty = huge(minsq_empty)
!       do i_periodic = 1, n_periodic
!          minsq_empty = min(minsq_empty, sum(lattice_vector(:, i_periodic)**2))
!       end do
!
!       ! check if there is an empty site on the multipole
!       do i_atom = 1, n_empty_atoms
!
!          ! collect distance
!          diff_vec = empty_coords(:, i_atom) - multipole_coords(:, i_multipole)
!          ! correct for periodicity
!          if (n_periodic > 0) then    ! Map to center cell
!             length = matmul(map_to_center_cell_matrix, diff_vec)
!             length = length - nint(length)
!             diff_vec = matmul(lattice_vector, length)
!          end if
!
!          thissq_empty = sum(diff_vec**2)
!
!          if (thissq_empty < minsq_empty) then
!             ! set shortest distance between multipole and empty sites
!             minsq_empty = thissq_empty
!          end if
!
!       end do 
!
!       if (minsq_empty.gt.0.0d0) then
!          ! check proximity to other atoms if no empty site is present
! 
!          do i_atom = 1, n_occ_atoms 
!             ! collect distance
!             diff_vec = occ_coords(:, i_atom) - multipole_coords(:, i_multipole)
!             ! correct for periodicity
!             if (n_periodic > 0) then    ! Map to center cell
!                length = matmul(map_to_center_cell_matrix, diff_vec)
!                length = length - nint(length)
!                diff_vec = matmul(lattice_vector, length)
!             end if
!          
!             thissq = sum(diff_vec**2)
!
!             ! ideally we would take the square root of these numbers,
!             ! then add the safety_margin and square again, as the relationship is non-linear,
!             ! but for now this comparison is safe enough, and simpler.
!
!             if (thissq .lt. ( atom_radius_sq(species(i_atom)) + safety_margin ) ) then
!                ! this multipole is too close to the quantum region and needs an empty site
!                ! with grids on top of it. Mark as an offending multipole.
!                offending_multipole(n_offending_multipoles) = i_multipole
!             end if   
!
!             ! returning the closest multipole is actually quite redundant,
!             ! but we will leave it here for now as we've already computed the distances.
!
!             if (thissq < minsq) then
!                ! this is checked irrespective of whether we are within the basis set radius
!                min_atom = i_atom
!                min_multipole = i_multipole
!                minsq = thissq
!             end if
!   
!          end do
!
!          if (offending_multipole(n_offending_multipoles).eq.i_multipole) then
!             ! increase n_offending_multipoles if this value has been added
!             n_offending_multipoles = n_offending_multipoles + 1
!          end if
!
!       end if
!
!    end do
!    dist = sqrt(minsq)
!
!  end subroutine min_multipole_dist
!
!

  !******
  !----------------------------------------------------------------------------
  !****s* bravais/get_map_to_center_cell_matrix
  !  NAME
  ! get_map_to_center_cell_matrix
  !  SYNOPSIS

  subroutine get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
  &                                        map_to_center_cell_matrix)

    !  PURPOSE
    !
    !  Updates map to center cell matrix. 
    !  Mind that :
    !     (1) map_to_center_cell_matrix == inverse of lattice vector matrix
    !         --> somewhat odd naming
    !     (2) in the case that the unit cell changes, the center cell matrix
    !         changes as well
    !
    !  USES

    use mpi_tasks, only: aims_stop
    use numerical_utilities, only: pseudo_inverse
    implicit none

    !  ARGUMENTS

    integer, intent(IN)  :: n_periodic
    real*8,  intent(IN)  :: lattice_vector(3, n_periodic)
    real*8,  intent(OUT) :: map_to_center_cell_matrix(n_periodic, 3)
    
    !  INPUTS
    !    o n_periodic -- Number of periodic dimensions
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    o map_to_center_cell_matrix -- inverse of lattice_vector
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank
    character(*), parameter :: func = 'get_map_to_center_cell_matrix'

    call pseudo_inverse(func, 3, n_periodic, lattice_vector, &
    &                   map_to_center_cell_matrix, 1d-10, rank)
    if(rank /= n_periodic) then
       call aims_stop('ERROR: lattice_vector is singular!')
    end if

  end subroutine get_map_to_center_cell_matrix
  !******
  !----------------------------------------------------------------------------
  !****s* get_cell_volume
  !  NAME
  !   get_cell_volume
  !  SYNOPSIS

  subroutine get_cell_volume(lattice_vector, cell_volume)

    !  PURPOSE
    !
    !       Calculates the supercell volume in the periodic systems
    !
    !  USES
    
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: lattice_vector(3, 3)
    real*8, intent(OUT) :: cell_volume

    !  INPUTS
    !    o lattice_vector -- Bravais lattice
    !  OUTPUT
    !    o cell_volume -- (Positive) volume of the unit cell.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    
    character*150 :: info_str
    character(*), parameter :: func = 'get_cell_volume'

    cell_volume = &
       abs(   lattice_vector(1,1) * lattice_vector(2,2) * lattice_vector(3,3)  &
            + lattice_vector(1,2) * lattice_vector(2,3) * lattice_vector(3,1)  &
            + lattice_vector(1,3) * lattice_vector(2,1) * lattice_vector(3,2)  &
            - lattice_vector(1,1) * lattice_vector(2,3) * lattice_vector(3,2)  &
            - lattice_vector(1,2) * lattice_vector(2,1) * lattice_vector(3,3)  &
            - lattice_vector(1,3) * lattice_vector(2,2) * lattice_vector(3,1)  &
          )

    if( cell_volume < 1d-10 ) then
       write(info_str, "(A,ES10.2)") &
       & 'Error: linear dependent lattice vectors; cell_volume =', cell_volume
       call aims_stop(info_str, func)
    end if

  end subroutine get_cell_volume
  !******
  !----------------------------------------------------------------------------
  !****s* get_length_of_lv
  !  NAME
  !       get_length_of_lv
  !  SYNOPSIS

  subroutine get_length_of_lv(n_periodic, lattice_vector, &
  &                           length_of_lattice_vector) 

    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(OUT) :: length_of_lattice_vector(n_periodic)

    !  INPUTS
    !    o n_periodic -- Number of periodic dimensions
    !    o lattice_vector -- Bravais lattice
    !  OUTPUT
    !    o length_of_lattice_vector -- Guess what
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_periodic


    do i_periodic = 1, n_periodic
       length_of_lattice_vector(i_periodic) = &
            sqrt(sum(lattice_vector(:,i_periodic)**2))
    end do

  end subroutine get_length_of_lv
  !******
  !----------------------------------------------------------------------------
  !****s* get_reciprocal_vectors
  !  NAME
  !       get_reciprocal_vectors
  !  SYNOPSIS

  subroutine get_reciprocal_vectors(n_periodic, lattice_vector, &
  &                                 recip_lattice_vector)

    !  PURPOSE
    !
    !   Calculates reciprocal lattice vectors for periodic systems from lattice
    !   vectors
    !
    !     b1 = 2p (a2 x a3 / a1 . (a2 x a3))
    !    
    !     b2 = 2p (a3 x a1 / a1 . (a2 x a3))
    !    
    !     b3 = 2p (a1 x a2 / a1 . (a2 x a3)) 
    !
    !  USES

    use constants, only: pi
    use mpi_tasks, only: aims_stop
    use numerical_utilities, only: pseudo_inverse
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(OUT) :: recip_lattice_vector(3, n_periodic)

    !  INPUTS
    !    o n_periodic -- Number of periodic dimensions
    !    o lattice_vector -- Bravais lattice
    !  OUTPUT
    !    o recip_lattice_vector -- Reciprocal lattice vectors
    !                              2*pi* transpose(lattice_vector^-1)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: inv_bra(n_periodic, 3)
    integer :: rank
    character(*), parameter :: func = 'get_reciprocal_vectors'

    call pseudo_inverse(func, 3, n_periodic, lattice_vector, inv_bra, &
    &                   1d-10, rank)
    if (rank /= n_periodic) call aims_stop('Invalid lattice_vector', func)
    recip_lattice_vector = 2*pi * transpose(inv_bra)

  end subroutine get_reciprocal_vectors
  !******
  !----------------------------------------------------------------------------
  !****s* bravais/get_n_supercells_pairs
  !  NAME
  !    get_n_supercells_pairs
  !  SYNOPSIS

  subroutine get_n_supercells_pairs(n_periodic, lattice_vector, &
  &                                 n_coords, coords, Rmax, n_supercells)

    !  PURPOSE
    !
    !     Figure out how many unit cells we must take into account to capture
    !     all pairs of atoms not separated further than Rmax.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    integer, intent(IN) :: n_coords
    real*8, intent(IN) :: coords(3, n_coords)
    real*8, intent(IN) :: Rmax
    integer, intent(OUT) :: n_supercells(n_periodic)

    !  INPUTS
    !    o n_periodic -- Number of periodic directions
    !                    [n_supercells(n_periodic+1:)==0]
    !    o lattice_vector(:,i) -- i-th lattice vector
    !    o n_coords -- number of atoms
    !    o coords -- atomic coordinates
    !    o Rmax -- maximum pair distance
    !  OUTPUTS
    !    o n_supercells -- number of cells for each Bravais direction
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Bmax, Tmax, vec(3)
    integer :: i_coord
    character(*), parameter :: func = 'get_n_supercells_pairs'

    ! Tmax is the largest distance of an coord to the origin
    Tmax = 0.d0
    do i_coord = 1, n_coords
       vec = coords(:, i_coord)
       Tmax = max(Tmax, sqrt(dot_product(vec, vec)))
    end do
    Bmax = Rmax + 2*Tmax

    call get_n_supercells(n_periodic, lattice_vector, Bmax, n_supercells)

  end subroutine get_n_supercells_pairs
  !******
  !----------------------------------------------------------------------------
  !****s* bravais/get_n_supercells
  !  NAME
  !    get_n_supercells
  !  SYNOPSIS

  subroutine get_n_supercells(n_periodic, lattice_vector, Rmax, n_supercells)

    !  PURPOSE
    !
    !     Figure out how many unit cells we must take into account to capture
    !     all Bravais vectors shorter than Rmax.
    !
    !  USES

    use constants, only: pi
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(IN) :: Rmax
    integer, intent(OUT) :: n_supercells(n_periodic)

    !  INPUTS
    !    o n_periodic -- Number of periodic directions
    !                    [n_supercells(n_periodic+1:)==0]
    !    o lattice_vector(:,i) -- i-th lattice vector
    !    o n_coords -- number of atoms
    !    o coords -- atomic coordinates
    !    o Rmax -- maximum pair distance
    !  OUTPUTS
    !    o n_supercells -- number of cells for each Bravais direction
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    [1] Richard M. Martin, "Electronic Structure: Basic Theory and
    !        Practical Methods", 1st ed., Cambridge University Press (2004).
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: rec(3, n_periodic), len_rec(n_periodic)
    integer :: i
    character(*), parameter :: func = 'get_n_supercells'

    n_supercells = 0
    if (n_periodic > 0) then
       call get_reciprocal_vectors(n_periodic, lattice_vector, rec)
       do i = 1, n_periodic
          len_rec(i) = sqrt(dot_product(rec(:, i), rec(:, i)))
       end do
       do i = 1, n_periodic
          ! See, e.g., [1, Chap. 4.2] or have a look at the routine
          ! "initialize_recip_hartree_potential" in hartree_potential_recip.f90
          ! where the following formula is derived - but for the case of points
          ! in _reciprocal_ space.
          n_supercells(i) = floor(len_rec(i) * Rmax / (2*pi) + 1d-5)
       end do
    end if

  end subroutine get_n_supercells
  !******
  !----------------------------------------------------------------------------
  !****s* bravais/get_neighboring_wigner_seitz
  !  NAME
  !    get_neighboring_wigner_seitz
  !  SYNOPSIS

  subroutine get_neighboring_wigner_seitz(n_periodic, lattice_vector, &
  &                                       max_n_ngh, n_ngh, a_ngh)

    !  PURPOSE
    !
    !    Get all Bravais vectors to neighboring Wigner-Seitz cells.  From these,
    !    the central Wigner-Seitz cell can be constructed.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    integer, intent(IN) :: max_n_ngh
    integer, intent(OUT) :: n_ngh
    integer, intent(OUT) :: a_ngh(n_periodic, max_n_ngh)

    !  INPUTS
    !    o n_periodic -- Number of periodic directions
    !    o lattice_vector(:,i) -- i-th lattice vector
    !    o max_n_ngh -- Guess on maximum number of neighbors
    !                   JW: I cannot imagine a case where 14 is not enough.
    !  OUTPUTS
    !    o n_ngh -- Number of neighbors
    !    o a_ngh -- Lattice indices of neighbors:
    !               Rvec(:, i_ngh) = matmul(lattice_vector, a_ngh(:, i_ngh)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i, n_supercells(3), a1, a2, a3, a123(3)
    integer :: n_cand, n_cand_uptonow, i_cand, j_cand
    integer :: n_ngh_uptonow
    real*8 :: Rmax, Rvec(3), R0vec(3)
    real*8, allocatable :: candvec(:,:)
    integer, allocatable :: canda123(:,:)
    integer :: info
    character(*), parameter :: func = 'get_neighboring_wigner_seitz'

    Rmax = 0.d0
    do i = 1, n_periodic
       Rmax = max(Rmax, 2*sqrt(sum(lattice_vector(:, i)**2)))
    end do
    n_supercells = 0   ! Important for n_periodic < 3.
    call get_n_supercells(n_periodic, lattice_vector, Rmax, n_supercells)

    ! --- Get candidates

    n_cand = product(2*n_supercells+1) - 1  ! no origin
    allocate(candvec(3, n_cand), canda123(3, n_cand), stat=info)
    call check_allocation(info, 'candvec', func)
    n_cand_uptonow = 0
    do a1 = -n_supercells(1), n_supercells(1)
       do a2 = -n_supercells(2), n_supercells(2)
          do a3 = -n_supercells(3), n_supercells(3)
             a123 = (/a1, a2, a3/)
             if (all(a123 == 0)) cycle
             n_cand_uptonow = n_cand_uptonow + 1
             canda123(:, n_cand_uptonow) = a123
             candvec(:, n_cand_uptonow) = matmul(lattice_vector, dble(a123))
          end do
       end do
    end do
    n_cand = n_cand_uptonow

    ! --- Check them

    n_ngh_uptonow = 0
    OUTER_L: do i_cand = 1, n_cand
       R0vec = candvec(:, i_cand)
       do j_cand = 1, n_cand
          if (i_cand == j_cand) cycle
          Rvec = candvec(:, j_cand)
          ! If the midpoint between the origin and R0vec is closer to
          ! Rvec than to R0vec, R0vec cannot be a neighbor.
          if (dot_product(R0vec, Rvec) > dot_product(Rvec, Rvec) - 1d-10) then
             cycle OUTER_L
          end if
       end do
       ! If the midpoint between the origin and R0vec is closer to R0vec
       ! than to any other Bravais vector, it is a neighbor.
       n_ngh_uptonow = n_ngh_uptonow + 1
       if (n_ngh_uptonow > max_n_ngh) then
          call aims_stop('max_n_ngh too small', func)
       end if
       a_ngh(:, n_ngh_uptonow) = canda123(:, i_cand)
    end do OUTER_L
    n_ngh = n_ngh_uptonow

    deallocate(candvec, canda123)

  end subroutine get_neighboring_wigner_seitz
  !******
  !----------------------------------------------------------------------------
  !****s* bravais/check_if_lattice_point
  !  NAME
  !    check_if_lattice_point
  !  SYNOPSIS

  subroutine check_if_lattice_point(n_periodic, lattice_vector, thres, &
  &                                 vec, is_lattice_point)

    !  PURPOSE
    !
    !     Check if given point is equivalent to lattice point.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    use numerical_utilities, only: solve_lsq
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_periodic
    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(IN) :: thres
    real*8, intent(IN) :: vec(3)
    logical, intent(OUT) :: is_lattice_point

    !  INPUTS
    !    o n_periodic -- Number of periodic directions
    !    o lattice_vector(:,i) -- i-th lattice vector
    !    o thres -- Threshold for accepting as equal
    !    o vec -- Vector to test
    !  OUTPUTS
    !    o is_lattice_point -- .true. if vec is integer
    !                          multiple of lattice_vector
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8 :: coeff(n_periodic), mapped_vec(3)
    integer :: rank
    character(*), parameter :: func = 'check_if_lattice_point'

    ! Least squares solution of Ax = b
    ! call solve_LSQ(func, M, N, A, b, x, rank, sing_thres)
    call solve_LSQ(func, 3, n_periodic, lattice_vector, vec, coeff, rank)
    if (rank /= n_periodic) then
       call aims_stop(func, 'Linear dependent lattice vectors')
    end if
    coeff = coeff - nint(coeff)
    mapped_vec = matmul(lattice_vector, coeff)
    is_lattice_point = all(abs(mapped_vec) < thres)

  end subroutine check_if_lattice_point
  !******
end module bravais
!******
