module sym_base

  !  PURPOSE
  !
  !  Base module for symmetry related routines
  !  USES
  !    none
  !  ARGUMENTS
  implicit none
  save
  
  ! Inverse symmetry operations 
  integer, dimension(:,:,:), allocatable :: rotations_reciprocal
  real*8, dimension(:,:), allocatable :: translations_reciprocal
  
  ! Mapping arrays (grid)
  integer, dimension(:,:),    allocatable :: map_atom
  integer, dimension(:,:),    allocatable :: map_atom_full
  integer, dimension(:),  allocatable :: map_centers
  integer, dimension(:),  allocatable :: map_centers_sym 
  real*8, dimension(:,:),    allocatable :: map_real
  integer, dimension(:,:),    allocatable :: map_grid

  ! Local partition_tab
  real*8, dimension(:, :), allocatable :: atom_atom_tab_sym
  real*8, dimension(:), allocatable :: atom_atom_dist_list_sym
  integer, dimension(:), allocatable :: atom_idx_A_sym
  integer, dimension(:), allocatable :: atom_idx_B_sym
  real*8, dimension(:),    allocatable :: min_atom_atom_tab_sym
  integer, dimension(:),   allocatable :: atom_atom_index_sym
  integer :: n_atom_atom_tab_sym
  real*8, dimension(:,:,:),    allocatable :: partition_tab_sym
  real*8, dimension(:,:,:),    allocatable :: hartree_partition_tab_sym 
  
  real*8, dimension(:),    allocatable ::  partition_tab_sym_batches
  real*8, dimension(:),    allocatable ::  hartree_partition_tab_sym_batches

  ! Y_lm transformation matrix
  complex*16, dimension(:,:,:,:), allocatable :: Delta
  complex*16, dimension(:,:,:,:), allocatable :: Delta_matrix_full
  integer, dimension(:), allocatable ::  rotations_at_task_map
  integer, dimension(:), allocatable ::  rotations_at_task

  real*8, dimension(:,:), allocatable ::  k_point_list_nosym
  integer, dimension(:), allocatable ::  k_points_at_task
  integer, dimension(:), allocatable ::  n_k_per_task
  integer, dimension(:,:), allocatable ::  k_per_task
  integer, dimension(:), allocatable ::  n_k_per_irr
  integer, dimension(:,:), allocatable ::  k_per_irr
  complex*16, dimension(:,:), allocatable ::  k_phase_nosym
  real*8, dimension(:,:,:), allocatable :: occ_numbers_nosym
  real*8, dimension(:,:,:), allocatable :: occ_numbers_sym  
  real*8, dimension(:), allocatable ::  k_weights_nosym
  complex*16,dimension(:,:),allocatable :: k_phase_base_nosym
  integer::  n_k_points_task_full
  complex*16, dimension(:,:), allocatable ::  k_phase_exx_nosym
  integer, dimension(:), allocatable ::    map_back
  integer :: n_rotations_at_task
  integer ::   n_kq_points, n_kq_points_task
  integer, dimension(:), allocatable :: k_points_coul
 contains
!****s* FHI-aims/symmat/map_centers_syms
!  NAME
!   map_centers_syms
!  SYNOPSIS
subroutine map_centers_syms(n_centers,coords_center,&
                            center_to_atom, num_symmetry,spg_shift, spg_rotations)

  !  PURPOSE
  !
  !    Map multiplicity of center to map_centers and save index of rotation matrix 
  !    in spg_rotations to map_centers_sym
  !
  !  USES
  use localorb_io, only: localorb_info, use_unit
  use geometry, only: lattice_vector, frac2cart, cart2frac
  use mpi_tasks, only: check_allocation
  use constants, only: bohr
  implicit none
  !  ARGUMENTS
  integer, intent(in) :: n_centers
  real*8, intent(in) :: coords_center(3,n_centers)
  integer, intent(in) :: center_to_atom(n_centers)
  integer, intent(in) :: num_symmetry
  integer, intent(in) :: spg_rotations(3,3,num_symmetry)
  real*8, intent(in)  :: spg_shift(3,num_symmetry)
  !  INPUTS
  ! o n_centers - number of centers
  ! o coords_center - center coordinates
  ! o center_to_atom - center to atom map
  ! o num_symmetry - elements in space group
  ! o spg_rotations - Rotations
  ! o spg_shift - Translations
  !  OUTPUTS
  ! o map_centers - multiplicity of center
  ! o map_centers_sym - index of rotation matrix tramsforming center
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
  ! Local Variables

  integer :: info, is, n_new_centers
  character(len=120)  :: info_str
  real*8, dimension(3) :: coord_current, coord_current_rot, trans,&
                          atom_pos, coord_current_rot_scaled
  integer :: new_centers, n_real_centers, i_center, j_center
  
  ! Allocations
  if (.not.allocated(map_centers))then
    allocate ( map_centers(n_centers) ,stat=info)
      call check_allocation(info, 'map_centers                       ')
  endif
  if (.not.allocated(map_centers_sym))then
    allocate ( map_centers_sym(n_centers) ,stat=info)
      call check_allocation(info, 'map_centers_sym                   ')
  endif
  
  ! Init counters
  n_real_centers=0
  map_centers = 0
  map_centers_sym = 1
  new_centers = 0
  
  ! Loop over centers
  do i_center = 1, n_centers, 1
     n_real_centers = n_real_centers + 1
     map_centers(i_center) = map_centers(i_center) +1
     coord_current(:) = coords_center( :, i_center) 
     atom_pos = coords_center( :,center_to_atom(i_center))
     
     ! Loop over symmetry operations
     do is=2, num_symmetry, 1 !skip identity
     
       ! Rotate center
       call frac2cart(lattice_vector,spg_shift(1:3,is),trans)
       coord_current_rot = coord_current - trans
       call cart2frac(lattice_vector, coord_current_rot, coord_current_rot_scaled)
       coord_current_rot_scaled = matmul(inv(dble(spg_rotations(1:3,1:3,is))),coord_current_rot_scaled)
       call frac2cart(lattice_vector, coord_current_rot_scaled, coord_current_rot)
       
       ! Finde equivalent center
       do j_center = 1, i_center-1, 1  
          
          if (sqrt(sum((coord_current_rot-coords_center( :, j_center))**2))<1e-4)then
             map_centers(i_center) = map_centers(i_center)-1
             map_centers(j_center) = map_centers(j_center)+1
             new_centers = new_centers + 1
             map_centers_sym(n_real_centers) = is
             exit
          endif
       end do
       if (sqrt(sum((coord_current_rot-coords_center(:,j_center))**2))<1e-4)then
         exit
       endif
     end do
  end do
  
  ! Some output
  n_new_centers = 0
  do i_center = 1, n_centers, 1
    if (count(map_centers == i_center) .ne. 0) then
      n_new_centers = n_new_centers + 1
    endif
  enddo
  !write(use_unit,*) map_centers
  write(info_str, '(2X,A,I8,A,I8)') &
  & '| centers reduced from: ', n_centers, ' to ', n_new_centers
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str, '(2X,A,I8,A,I8)') &
  & '| centers reduced from: ', n_centers, ' to ', count(map_centers .ne. 0)
  call localorb_info(info_str,use_unit,'(A)')

end subroutine map_centers_syms
!******
!****s* FHI-aims/symmat/map_atoms
!  NAME
!   map_atoms
!  SYNOPSIS
subroutine map_atoms()

  !  PURPOSE
  !
  !    Map between original and rotated/translated atom 
  !    for all operations in spg_rotations (num_symmetry) to map_atom.
  !
  !  USES
  use spglib_symmetry,only: spg_shift, spg_rotations
  use dimensions,only: n_atoms, n_centers
  use localorb_io, only: localorb_info
  use geometry, only: lattice_vector, frac2cart, cart2frac, species
  use mpi_tasks, only: check_allocation
  use pbc_lists,only: center_to_atom,coords_center
  use constants, only: bohr
  implicit none
  !  ARGUMENTS
  !  INPUTS
  !  OUTPUTS
  ! o map_atom -- Map between original and rotated/translated atom
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
  ! Local Variables


  real*8, dimension(3) :: coord_current, coord_current_rot, center_coord
  real*8, dimension(3,3) :: temp
  integer :: i_atom, i_center, is, info, s_a_current, s_c_current
  real*8 :: diff
  logical :: found
  
  if (.not.allocated(map_atom))then
      allocate ( map_atom(n_atoms,n_rotations_at_task) ,stat=info)
      call check_allocation(info, 'map_atom                  ')
  endif 
  ! Loop over atoms
  do i_atom = 1, n_atoms, 1
     call cart2frac(lattice_vector,coords_center( :, i_atom),coord_current) 
     ! Loop over symmetry operations
     do is=1, n_rotations_at_task, 1
       coord_current_rot = coord_current - spg_shift(1:3,rotations_at_task(is))
       temp = inv(dble(spg_rotations(1:3,1:3,rotations_at_task(is))))
       coord_current_rot = matmul(temp,coord_current_rot)
       
       map_atom(i_atom,is) = i_atom
       ! find equivalent center
       found = .false.
       do i_center = 1, n_centers, 1
         call cart2frac(lattice_vector,coords_center(:,i_center),center_coord)
         diff = sqrt(sum((coord_current_rot-center_coord)**2))
         s_a_current = species(i_atom)
         s_c_current = species(center_to_atom(i_center))
         if (diff<1e-4.and.(s_a_current.eq.s_c_current)) then
            map_atom(i_atom,is) = center_to_atom(i_center)
            found = .true.
            exit
         endif
       enddo
       
       ! debug write(use_unit,'(A)') 'Atom map'
       ! debug write(use_unit,'(A,I3,A,I4)') 'Symmetry op.: ',  rotations_at_tasks(is), ', Atom: ', i_atom
       ! debug write(use_unit,'(A,3E11.4,A,3E11.4)') 'Translation (frac, \AA): ',   &
       ! debug                                spg_shift(1:3,rotations_at_tasks(is)),', ', trans*bohr
       ! debug write(use_unit,'(A,L2)') 'Equivalent atom found: ', found
       ! debug write(use_unit,'(A,3E11.4,A,3E11.4)') 'Atom coord.: ', coord_current, ',&
       ! debug                                Rotated coord: ', coord_current_rot
       ! debug call cart2frac(lattice_vector, coords_center(1:3,&
       ! debug                map_atom(i_atom,is)),center_coord)
       ! debug write(use_unit,'(A,I4,A,3E11.4)') 'Equivalent atom: ', &
       ! debug                             center_to_atom(map_atom(i_atom,is)),&
       ! debug                            ', Coord.: ', center_coord
       ! debug center_coord = coords_center(:,map_atom(i_atom,is))
       ! debug call map_to_center_cell(center_coord)
       ! debug write(use_unit,*) center_coord,center_coord*bohr
     end do
  end do

end subroutine map_atoms

subroutine map_atoms_all()

  !  PURPOSE
  !
  !    Map between original and rotated/translated atom 
  !    for all operations in spg_rotations (num_symmetry) to map_atom.
  !
  !  USES
  use spglib_symmetry,only: spg_shift, spg_rotations, num_symmetry
  use dimensions,only: n_atoms,n_centers
  use localorb_io, only: localorb_info
  use geometry, only: lattice_vector, frac2cart, cart2frac, species
  use mpi_tasks, only: check_allocation
  use pbc_lists,only: center_to_atom,coords_center
  use constants, only: bohr
  implicit none
  !  ARGUMENTS
  !  INPUTS
  !  OUTPUTS
  ! o map_atom -- Map between original and rotated/translated atom
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
  ! Local Variables


  real*8, dimension(3) :: coord_current, coord_current_rot, center_coord
  real*8, dimension(3,3) :: temp
  integer :: i_atom, i_center, is, info, s_a_current, s_c_current
  real*8 :: diff
  logical :: found
  
  if (.not.allocated(map_atom_full))then
      allocate ( map_atom_full(n_atoms,num_symmetry) ,stat=info)
      call check_allocation(info, 'map_atom_full                 ')
  endif 
  ! Loop over atoms
  do i_atom = 1, n_atoms, 1
     call cart2frac(lattice_vector,coords_center( :, i_atom),coord_current) 
     ! Loop over symmetry operations
     do is=1, num_symmetry, 1
       coord_current_rot = coord_current - spg_shift(1:3,(is))
       temp = inv(dble(spg_rotations(1:3,1:3,(is))))
       coord_current_rot = matmul(temp,coord_current_rot)
       
       map_atom_full(i_atom,is) = i_atom
       ! find equivalent center
       found = .false.
       do i_center = 1, n_centers, 1
         call cart2frac(lattice_vector,coords_center(:,i_center),center_coord)
         diff = sqrt(sum((coord_current_rot-center_coord)**2))
         s_a_current = species(i_atom)
         s_c_current = species(center_to_atom(i_center))
         if (diff<1e-4.and.(s_a_current.eq.s_c_current)) then
            map_atom_full(i_atom,is) = center_to_atom(i_center)
            found = .true.
            exit
         endif
       enddo
       
       ! debug write(use_unit,'(A)') 'Atom map'
       ! debug write(use_unit,'(A,I3,A,I4)') 'Symmetry op.: ',  rotations_at_tasks(is), ', Atom: ', i_atom
       ! debug write(use_unit,'(A,3E11.4,A,3E11.4)') 'Translation (frac, \AA): ',   &
       ! debug                                spg_shift(1:3,rotations_at_tasks(is)),', ', trans*bohr
       ! debug write(use_unit,'(A,L2)') 'Equivalent atom found: ', found
       ! debug write(use_unit,'(A,3E11.4,A,3E11.4)') 'Atom coord.: ', coord_current, ',&
       ! debug                                Rotated coord: ', coord_current_rot
       ! debug call cart2frac(lattice_vector, coords_center(1:3,&
       ! debug                map_atom(i_atom,is)),center_coord)
       ! debug write(use_unit,'(A,I4,A,3E11.4)') 'Equivalent atom: ', &
       ! debug                             center_to_atom(map_atom(i_atom,is)),&
       ! debug                            ', Coord.: ', center_coord
       ! debug center_coord = coords_center(:,map_atom(i_atom,is))
       ! debug call map_to_center_cell(center_coord)
       ! debug write(use_unit,*) center_coord,center_coord*bohr
     end do
  end do

end subroutine map_atoms_all
!****s* FHI-aims/symmat/get_inv_sym
!  NAME
!   get_inv_sym
!  SYNOPSIS
subroutine get_inv_sym(num_symmetry, rotations, translations, is_time_reversal,&
                       num_inv_symmetry)
  !  PURPOSE
  !
  !    Invert rotation matrix in rotations.
  !
  !  USES
  use mpi_tasks, only: check_allocation
  implicit none
  !  ARGUMENTS
  integer, intent(in) :: num_symmetry
  integer, dimension(3,3,num_symmetry), intent(in) :: rotations
  real*8, dimension(3,num_symmetry), intent(in) :: translations
  integer, intent(in) :: is_time_reversal
  integer, intent(out) :: num_inv_symmetry
  !  INPUTS
  ! o num_symmetry -- number of elements in space group
  ! o rotations -- rotation matrices
  ! o translations -- translation vectors
  ! o is_time_reversal -- use time reversal
  !  OUTPUTS
  ! o num_inv_symmetry -- number of inverse operations
  ! o rotations_reciprocal
  ! o translations_reciprocal
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
  ! Local Variables

  integer :: info, i, j
  integer, dimension(:,:,:), allocatable :: rotations_work
  real*8, dimension(:,:), allocatable :: translations_work
  integer, dimension(:),  allocatable :: unique_rot
  integer,  dimension(3,3) :: inversion
  integer :: num_rot, num_int



  inversion(1,1:3) =  (/-1, 0, 0 /)
  inversion(2,1:3) =  (/ 0,-1, 0 /)
  inversion(3,1:3) =  (/ 0, 0,-1 /)

  if (.not.allocated(rotations_work))then
    if (is_time_reversal==1)then 
      allocate ( rotations_work(3,3,2*num_symmetry) ,stat=info)
    else
      allocate ( rotations_work(3,3,num_symmetry) ,stat=info)
    endif
      call check_allocation(info, 'rotations_work                 ')
  endif
  if (.not.allocated(translations_work))then
    if (is_time_reversal==1)then 
      allocate ( translations_work(3,2*num_symmetry) ,stat=info)
    else
      allocate ( translations_work(3,num_symmetry) ,stat=info)
    endif
      call check_allocation(info, 'translations_work                 ')
  endif
  if (.not.allocated(unique_rot))then
    if (is_time_reversal==1)then 
      allocate ( unique_rot(2*num_symmetry) ,stat=info)
    else
      allocate ( unique_rot(num_symmetry) ,stat=info)
    endif
      call check_allocation(info, 'unique_rot                 ')
  endif

  rotations_work = 0
  unique_rot = 0
  translations_work = 0d0



  do i = 1, num_symmetry, 1
    rotations_work(:,:,i) = transpose(rotations(:,:,i))
    translations_work(:,i) = translations(:,i)
    if (is_time_reversal==1)then 
      rotations_work(:,:,num_symmetry+i) = matmul(rotations_work(:,:,i),&
                         inversion)
      translations_work(:,num_symmetry+i) = translations(:,i)
    endif
  enddo


  num_rot = 1;
  if (is_time_reversal==1)then
    num_int = 2*num_symmetry
  else
    num_int = num_symmetry
  endif 

  do i = 1, num_int, 1
    unique_rot(i) = 1
  enddo
  do i =1, num_int, 1
    do j = 1, num_rot, 1
      if(all(rotations_work(:,:,unique_rot(j)).eq.rotations_work(:,:,i))) then
        cycle
      endif
    enddo
    unique_rot(num_rot) = i;
    num_rot = num_rot + 1
  enddo
  num_rot = num_rot - 1

  if (.not.allocated(rotations_reciprocal))then
      allocate ( rotations_reciprocal(3,3,num_rot) ,stat=info)
      call check_allocation(info, 'rotations_reciprocal             ')
  endif
  if (.not.allocated(translations_reciprocal))then
      allocate ( translations_reciprocal(3,num_rot) ,stat=info)
      call check_allocation(info, 'translations_reciprocal             ')
  endif
  do i = 1, num_rot, 1
      translations_reciprocal(:,i) = translations_work(:,unique_rot(i))
      rotations_reciprocal(:,:,i) = rotations_work(:,:,unique_rot(i))
  enddo
  num_inv_symmetry = num_rot

  if(allocated(rotations_work)) deallocate(rotations_work)
  if(allocated(translations_work)) deallocate(translations_work)
  if(allocated(unique_rot)) deallocate(unique_rot)
end subroutine get_inv_sym

!****s* FHI-aims/symmat/symmattoeuler
!  NAME
!   symmattoeuler
!  SYNOPSIS
subroutine symmattoeuler(symmat,alpha, beta, gamma)
  !  PURPOSE
  !
  !    Get euler angles from rotation matrix.
  !    see https://de.wikipedia.org/wiki/Eulersche_Winkel
  !    Given a rotation matrix
  !    \begin{align}
  !    M_{zyz} & = 
  !       \begin{pmatrix}
  !         \cos \gamma  & \sin \gamma  & 0 \\
  !        -\sin \gamma  & \cos \gamma  & 0 \\
  !           0          & 0            & 1
  !      \end{pmatrix}
  !      \begin{pmatrix}
  !         \cos \beta   & 0   & -\sin \beta \\
  !           0          & 1   &    0 \\
  !         \sin \beta   & 0   & \cos \beta
  !      \end{pmatrix}
  !      \begin{pmatrix}
  !         \cos \alpha  & \sin \alpha  & 0 \\
  !        -\sin \alpha  & \cos \alpha  & 0 \\
  !           0          & 0            & 1
  !      \end{pmatrix} \\
  !    \end{align}
  !    \begin{align}
  !    M_{zyz}^T & = 
  !     \begin{pmatrix}
  !       -\sin\alpha \sin\gamma  + \cos\alpha \cos\beta \cos\gamma
  !    &  -\sin\alpha \cos\gamma  - \cos\alpha \cos\beta \sin\gamma
  !    &   \cos\alpha \sin\beta \\
  !        \cos\alpha \sin\gamma  + \sin\alpha \cos\beta \cos\gamma
  !    &   \cos\alpha \cos\gamma  - \sin\alpha \cos\beta \sin\gamma
  !    &   \sin\alpha \sin\beta    \\
  !       -\sin\beta \cos\gamma
  !    &   \sin\beta \sin\gamma
  !    &   \cos \beta
  !    \end{pmatrix}
  !    \end{align}
  !   this routine determines the Euler angles, $(\alpha,\beta,\gamma)$. This
  !   corresponds to the Standard-y-Konvention (z, y′, z″), which involves the 
  !   following successive rotations of the coordinate system:
  !   \begin{itemize}
  !    \item[1.]{The $x'$-, $y'$-, $z'$-axes are rotated anticlockwise
  !     through an angle $\alpha$ about the z' axis}
  !    \item[2.]{The $x''$-, $y''$-, $z''$-axes are rotated anticlockwise
  !     through an angle $\beta$ about the $y'$ axis}
  !    \item[3.]{The $x'''$-, $y'''$-, $z'''$-axes are rotated anticlockwise
  !     through an angle $\gamma$ about the $x_3''$ axis}
  !   \end{itemize}
  !   The Euler angles are not necessarily unique for a given rotation matrix.  
  !
  !  USES
  use constants, only: pi
  use mpi_tasks, only: aims_stop

  implicit none
  !  ARGUMENTS
  real*8, dimension(3,3), intent(in) :: symmat
  real*8, intent(out) :: alpha, beta, gamma
  !  INPUTS
  !   o symmat - Rotation matrix
  !  OUTPUTS
  ! Euler angles
  !   o alpha
  !   o beta
  !   o gamma
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
  ! Local Variables


  real*8, parameter :: safe_minimum = 1.d-10
  real*8 :: det 

  det = symmat(1,2) * symmat(2,3) * symmat(3,1) &
      - symmat(1,3) * symmat(2,2) * symmat(3,1) &
      + symmat(1,3) * symmat(2,1) * symmat(3,2) &
      - symmat(1,1) * symmat(2,3) * symmat(3,2) &
      + symmat(1,1) * symmat(2,2) * symmat(3,3) &
      - symmat(1,2) * symmat(2,1) * symmat(3,3) 

  if ((det > 1.d0+safe_minimum) .or. (det < 1.d0-safe_minimum)) then
    call aims_stop("symmattoeuler: Symmetry matrix improper or not unitary")
  endif

  if ((abs(symmat(3, 1)) > safe_minimum) .or. (abs(symmat(3, 2)) > safe_minimum)) then
    
    alpha = atan2(symmat(3, 2), symmat(3, 1))
    if (alpha < 0d0) then
      alpha = alpha + 2d0*pi
    endif
    if (abs(symmat(3, 1)) > abs(symmat(3, 2))) then
      beta = atan2 (symmat(3,1)/cos(alpha), symmat(3,3))
    else
      beta = atan2 (symmat(3,2)/Sin(alpha), symmat(3,3))
    endif
    gamma = atan2(symmat(2,3),-symmat(1,3))
    if (gamma < 0d0) then
      gamma = gamma + 2d0*pi
    endif
  else
    alpha = atan2(symmat(1, 2), symmat(1, 1))
    if (alpha < 0d0) then
      alpha = alpha + 2d0*pi
    endif
    if (symmat(3,3) > 0d0) then
      beta = 0.d0
      gamma = 0.d0
    else
      beta = pi
      gamma = pi
    endif    

  endif

end subroutine symmattoeuler
!****s* FHI-aims/symmat/get_TVlmm
!  NAME
!   get_TVlmm
!  SYNOPSIS
subroutine get_TVlmm(TVlmm,alpha, beta, gamma,p)
  !  PURPOSE
  !
  !    Calculate the transformation matrix T(V,l,m,m') for Y_lm function from
  !    the euler angle with the recursive algorithm described in
  !    -Blanco, Florez, Bermejo, Journal of Molecular strucutre (TheoChem) 419 (1997) 19-27
  !    -Kudlicki, Rowicki, Gilski, Olwinowski, J. Appl. Cryst. (2005) 38, 501-504
  !    aka Wigner D-Matrix for complex spherical harmonics
  !
  !  USES
  use dimensions, only: l_wave_max
  use constants, only: pi
  implicit none

  
  !  ARGUMENTS
  real*8, intent(in) :: alpha, beta, gamma
  integer, intent(in) :: p
  complex*16 :: TVlmm(0:l_wave_max,-l_wave_max:l_wave_max,-l_wave_max:l_wave_max)
  !logical :: inversion
  !  INPUTS
  !  Euler angles
  !  o alpha
  !  o beta
  !  o gamma
  !  OUTPUTS
  !  o  TVlmm - Y_lm transformation matrix
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
  ! Local Variables

  integer :: l,m1,m2
  real*8  :: phase_ang
  ! Initialize
  do l =0, l_wave_max, 1
    do m1 = -l_wave_max, l_wave_max, 1
      do m2 = -l_wave_max, l_wave_max, 1
         TVlmm(l,m1,m2) = 10d0
      enddo
    enddo
  enddo

 
  do l = 0, l_wave_max , 1
    ! l = 0
    if (l==0) then
      TVlmm(l,0,0) = 1d0
    ! l = 1
    elseif (l==1) then
      TVlmm(l,0,0) = cos(beta)
      TVlmm(l,1,0) = -1d0/sqrt(2d0)*sin(beta)
      TVlmm(l,0,1) = (-1d0)*TVlmm(l,1,0)
      TVlmm(l,-1,0) = (-1d0)*TVlmm(l,1,0)
      TVlmm(l,0,-1) = TVlmm(l,1,0)
      TVlmm(l,1,-1) = sin(beta*0.5)**2
      TVlmm(l,-1,1) = TVlmm(l,1,-1)
      TVlmm(l,1,1) = 1d0/2d0*(1 + cos(beta))
      TVlmm(l,-1,-1) = TVlmm(l,1,1)
    ! l > 1
    else
      ! m1 = 0,..,l-2; m2 =-m1,..,m1 
      do m1 = 0, l-2, 1
        do m2 = -m1, m1, 1     
           TVlmm(l,m1,m2) = (1d0*(2*l**2-l)/(((l**2-m1**2)*(l**2-m2**2))**(0.5)))*&
                            ((TVlmm(1,0,0)-(1d0*(m1*m2)/(l*(l-1))))*TVlmm(l-1,m1,m2)&
                            -(((((l-1)**2-m1**2)*((l-1)**2-m2**2))**0.5)*1d0/((l-1)&
                            *(2*l-1)))*TVlmm(l-2,m1,m2))
           ! mirror properties
           TVlmm(l,m2,m1) = ((-1d0)**(m1-m2))*TVlmm(l,m1,m2)
           TVlmm(l,-m1,-m2) = ((-1d0)**(m1-m2))*TVlmm(l,m1,m2)
           TVlmm(l,-m2,-m1) = ((-1d0)**(m2-m1))*TVlmm(l,m2,m1)
        enddo
      enddo
      ! m1=l && m2=2 / m1 =l-1 $$ m2 =l-1
      TVlmm(l,l,l) = cos(beta*0.5)**(2d0*l)
      TVlmm(l,l-1,l-1) =  (l*cos(beta)-l+1)*(cos(beta*0.5)**(2d0*(l-1)))
      ! mirror properties
      TVlmm(l,-l,-l) = TVlmm(l,l,l)
      TVlmm(l,-(l-1),-(l-1)) = TVlmm(l,l-1,l-1)
      ! special case beta == pi
      if (beta.eq.pi)then
        TVlmm(l,l-1:l,:) = 0
        TVlmm(l,:,l-1:l) = 0
        TVlmm(l,:,-l:-l+1) = 0d0
        TVlmm(l,-l:-l+1,:) = 0d0
        TVlmm(l,l,-l) = 1d0
        TVlmm(l,-l,l) = 1d0
        TVlmm(l,l-1,-(l-1)) = (-1d0)
        TVlmm(l,-(l-1),l-1) = (-1d0)
      else
        ! m2 = l-1,..,-l
        do m2 = l-1, -l, -1
          TVlmm(l,l,m2) = -((dble(l+(m2+1))/dble(l-(m2+1)+1))**0.5)*TVlmm(l,l,m2+1)*tan(beta*0.5)
          TVlmm(l,m2,l) = ((-1d0)**(l-(m2)))*TVlmm(l,l,m2)
          ! mirror properties
          TVlmm(l,-l,-m2) = ((-1d0)**(l-(m2)))*TVlmm(l,l,m2)
          TVlmm(l,-m2,-l) = ((-1d0)**((m2-l)))*TVlmm(l,m2,l)
        enddo
        ! m2 = l-2,..,1-l
        do m2 = l-2, 1-l, -1
          ! special case beta =pi/2 or 3/2pi and m2==-1
          !if ((beta.eq.0.5*pi.or.beta.eq.1.5*pi).and.m2==-1)then
             TVlmm(l,l-1,m2)=(l*cos(beta)-m2)*((cos(beta*0.5))**(l-1+m2))*&
                             ((-sin(beta*0.5))**(l-1-m2))*&
                             sqrt(factn(2*l-1)/(factn(l+m2)*factn(l-m2)))
          !else
          !   TVlmm(l,l-1,m2) = (dble(-l*cos(beta)+(m2+1)-1)/dble(l*cos(beta)-(m2+1)))*&
          !                  ((dble(l+(m2+1))/dble(l-(m2+1)+1))**0.5)*TVlmm(l,l-1,m2+1)*tan(beta*0.5)
          !endif
          ! mirror properties
          TVlmm(l,m2,l-1) = ((-1d0)**((l-1)-(m2)))*TVlmm(l,l-1,m2)
          TVlmm(l,-(l-1),-m2) = ((-1d0)**((l-1)-(m2)))*TVlmm(l,l-1,m2)
          TVlmm(l,-m2,-(l-1)) = ((-1d0)**(m2-(l-1)))*TVlmm(l,m2,l-1)
        enddo             
      endif
    endif

      
  enddo
  ! debug write(use_unit,*)alpha, beta, gamma 
  ! m, m' dependency
  !debug write(use_unit,'(A)') '-----------------------------------------------------'
  !debug write(use_unit,'(A)') 'TVlmm'
  do l =0, l_wave_max, 1
    do m1 = -l, l, 1
      do m2 = -l, l, 1
          ! alpha and gamma dependency
          phase_ang=-dble(m1)*alpha-dble(m2)*gamma
          TVlmm(l,m1,m2) = TVlmm(l,m1,m2)*dcmplx(cos(phase_ang),sin(phase_ang))
          ! Parity
          if ((p .eq.-1) .and. (mod(l, 2) .ne. 0)) then
            TVlmm(l,m1,m2) = - TVlmm(l,m1,m2)
          endif
          !debug write(use_unit,'(A,I4,A,I4,A,I4,A,2F11.4)') 'l',l, ' m1',m1, ' m2', &
          !debug        m2, ' TVlmm:',real(TVlmm(l,m1,m2)),dimag(TVlmm(l,m1,m2))
      enddo
    enddo
  enddo
 !debug  write(use_unit,'(A)') '-----------------------------------------------------'
end subroutine get_TVlmm

!****s* FHI-aims/symmat/prepare_local_partition_tab
!  NAME
!   prepare_local_partition_tab
!  SYNOPSIS
subroutine prepare_local_partition_tab()

!  PURPOSE
!   initializaton and allocation for calculation of the local partiton tab
!
!
!  USES
!
      use dimensions,only:n_species,n_centers_basis_integrals
      use runtime_choices,only:partition_type,stratmann_a
      use pbc_lists,only:centers_basis_integrals
      use species_data,only:cut_free_atom,outer_partition_radius,free_r_cut,&
                            w_cutoff, multipole_radius_free
      use mpi_tasks,only:check_allocation
      use localorb_io,only:localorb_info,use_unit,OL_low,OL_high
      use constants,only:bohr
      implicit none

!  ARGUMENTS
! INPUTS
!   none

! OUTPUTS

!
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
!


!     local variables
      character*140 :: info_str
      integer :: info

!     counters

      integer :: i_species,i_atom,i_atom_2

      write (info_str,'(2X,A)') & 
        "Preparing partition"
      call localorb_info(info_str,use_unit,'(A)',OL_high)

      ! We begin by preparing the partition table of Stratmann and coworkers. This part
      ! became lengthy over time, please simply take it as one block 
      ! (the end of this part is also marked by a comment below)

      if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) &
            .or. (partition_type.eq.9) ) then
        ! All these are variants of the partition table by Stratmann and coworkers.
        ! Incremental improvements during development, keeping them does not hurt. 

        ! Determine what will be the outermost distance at which the stratmann partition
        ! table for a given atom should be non-zero

        ! outer_partition_radius is the outermost radius at which any integrand
        ! or density component will be partitioned to belong to a given atom.
        ! See the extended comment at the top of this subroutines for Equation
        ! numbers in the FHI-aims CPC and other applicable references.

        if (cut_free_atom(1)) then  
        ! there has to be at least one species. Referring to cut_free_atom(1) is an ugly hack, 
        ! but the output is explanatory only anyway.
            write (info_str,'(2X,A)') &
              "| initialize_grid_storage: Actual outermost partition radius vs. multipole_radius_free"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
            write (info_str,'(2X,A)') &
              "| (-- VB: in principle, multipole_radius_free should be larger, hence this output)"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
        end if

        do i_species = 1, n_species, 1
          ! If the Stratmann et al. partition table is used to determine the integration weights, then
          ! make sure that the outermost radius where the partition table elements for a given atom
          ! can be non-zero is bounded by the maximum extent of basis functions / densities associated
          ! with that atom.
          !
          ! See the comment at the top of this subroutine for references and equation numbers.
          !
          if (cut_free_atom(i_species)) then
            ! This is the case where the free-atom radius is bounded somehow - should be the normal
            ! case. If the free atom radius is not bounded, we have a problem for periodic systems -
            ! very extended integration regions must be accounted for.
            outer_partition_radius(i_species) = free_r_cut (i_species) + w_cutoff (i_species)
            ! ... and add some diagnostic output here.
            write (info_str,'(2X,A,I8,A,F30.15,A,F30.15,A)') &
              "| Species ", i_species, ": Confinement radius = ", outer_partition_radius(i_species)*bohr, &
              " AA, multipole_radius_free = ", multipole_radius_free(i_species)*bohr, " AA."
            call localorb_info(info_str,use_unit,'(A)',OL_low)
            !
            ! In the case of analytically tabulated Gaussian functions or hydrogenic orbitals,
            ! the extent of basis functions is not guarded by the confinement potential.
            ! In these cases (especially diffuse Gaussian functions) we may need to set a 
            ! larger outer boundary.
            !
            ! multipole_radius_free(i_species) is made sure to be larger than any radial function centered
            ! at this atom in subroutine shrink_fixed_basis_phi_thresh.f90 .
            !
            ! Hence, we employ multipole_radius_free here if it is larger. In the normal case,
            ! where both the free atom and the basis functions are subject to the same cutoff
            ! potential, the difference should be very small.
            !
            ! And yes, I purposely do not employ the 'max' function here. 
            if (multipole_radius_free(i_species).gt.outer_partition_radius(i_species)) then
               outer_partition_radius(i_species) = multipole_radius_free(i_species)
            end if
          else
            ! The free-atom radius is not bounded. Well, we must then use the outermost
            ! radius of the free atom density (very large), for better or for worse.
            outer_partition_radius(i_species) = multipole_radius_free(i_species)
          end if
          write (info_str,'(2X,A,I8,A,F30.15,A)') &
            "| Species ", i_species, ": outer_partition_radius set to ", outer_partition_radius(i_species)*bohr, &
            " AA ."
          call localorb_info(info_str,use_unit,'(A)')

        end do

      end if

      if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) ) then
        ! for stratmann partitioning scheme, make available only if necessary
        ! tabulate the interatomic distances for all relevant atoms
        allocate(atom_atom_tab_sym(n_centers_basis_integrals,n_centers_basis_integrals),stat=info)
        call check_allocation(info, 'atom_atom_tab_sym                 ')
      else
        ! Not needed, allocate dummy
        allocate(atom_atom_tab_sym(1,1))
      endif
      
      if (partition_type.eq.9) then
        ! for stratmann partitioning scheme, make available only if necessary
        ! tabulate the interatomic distances for all relevant atoms
        call tab_interatomic_distances_count_entries(n_centers_basis_integrals, centers_basis_integrals, n_atom_atom_tab_sym)

        allocate(atom_atom_dist_list_sym(n_atom_atom_tab_sym), stat=info)
        call check_allocation(info, 'atom_atom_dist_list_sym               ')

        allocate(atom_idx_A_sym(n_atom_atom_tab_sym), stat=info)
        call check_allocation(info, 'atom_idx_A_sym                    ')


        allocate(atom_idx_B_sym(n_centers_basis_integrals+1), stat=info)
        call check_allocation(info, 'atom_idx_B_sym               ')
      else
        ! Not needed, allocate dummy
        n_atom_atom_tab_sym = 1
        allocate(atom_atom_dist_list_sym(1))
        allocate(atom_idx_A_sym(1))
        allocate(atom_idx_B_sym(1))
      endif


      allocate(atom_atom_index_sym(n_centers_basis_integrals),stat=info)
      call check_allocation(info, 'atom_atom_index_sym               ')

      allocate(min_atom_atom_tab_sym(n_centers_basis_integrals),stat=info)
      call check_allocation(info, 'min_atom_atom_tab_sym             ')


      if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) ) then
        ! FIXME:
        ! (1) count distance^2 and figure out which atoms could potentially be important in each other's stratmann partition tab
        ! (2) calculate min_distance^2
        ! (3) only allocate as much data as is required, otherwise the arrays may get too big ...
        ! (4) tabulate interatomic distances ONLY for those atoms needed.
        ! (5) remember how many relevant atoms there are for each possible i_atom
        ! the resulting list should be enough to build the stratmann list from it...
        call tab_interatomic_distances(n_centers_basis_integrals, centers_basis_integrals, atom_atom_tab_sym)
        ! for each atom, determine distance to next neighbour:
        do i_atom = 1, n_centers_basis_integrals
          min_atom_atom_tab_sym(i_atom) = 1d100
          do i_atom_2 = 1, n_centers_basis_integrals
            if (i_atom.ne.i_atom_2) then
              min_atom_atom_tab_sym(i_atom) = min(min_atom_atom_tab_sym(i_atom),atom_atom_tab_sym(i_atom,i_atom_2))
            end if
          end do
        end do
        ! need this for comparison ...
        min_atom_atom_tab_sym(:) = (1d0-stratmann_a)*min_atom_atom_tab_sym(:)/2d0

      else if (partition_type.eq.9) then
        call tab_interatomic_distances_local(n_centers_basis_integrals, centers_basis_integrals, &
                                             n_atom_atom_tab_sym, atom_atom_dist_list_sym, atom_idx_A_sym, atom_idx_B_sym)
        ! for each atom, determine distance to next neighbour:
        do i_atom = 1, n_centers_basis_integrals
          min_atom_atom_tab_sym(i_atom) = 1d100
          do i_atom_2 = atom_idx_B_sym(i_atom),atom_idx_B_sym(i_atom+1)-1
            min_atom_atom_tab_sym(i_atom) = min(min_atom_atom_tab_sym(i_atom),atom_atom_dist_list_sym(i_atom_2))
          end do
        end do
        ! need this for comparison ...
        min_atom_atom_tab_sym(:) = (1d0-stratmann_a)*min_atom_atom_tab_sym(:)/2d0
      end if

      ! End all preparations for the Stratmann and coworkers partition table.
end subroutine prepare_local_partition_tab
!****s* FHI-aims/sym_base/clean_local_partition_tab
!  NAME
!   clean_local_partition_tab
!  SYNOPSIS
subroutine clean_local_partition_tab()

!  PURPOSE
!   Deallocations for calculation of local partition tabs
!
!
!  USES
!
      use localorb_io,only:localorb_info,use_unit,OL_high
      implicit none

!  ARGUMENTS
! INPUTS
!   none

! OUTPUTS

!
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
!


!     local variables
      character*140 :: info_str

!     counters


      ! clean up the odd allocated variable
      deallocate(min_atom_atom_tab_sym)
      deallocate(atom_atom_index_sym)
      deallocate(atom_atom_tab_sym)
      deallocate(atom_atom_dist_list_sym)
      deallocate(atom_idx_A_sym)
      deallocate(atom_idx_B_sym)

      write (info_str,'(2X,A)') & 
        "Clean partition"
      call localorb_info(info_str,use_unit,'(A)',OL_high)

end subroutine clean_local_partition_tab
!****s* FHI-aims/sym_base/local_partition_tab
!  NAME
!   local_partition_tab
!  SYNOPSIS
subroutine local_partition_tab( partition_tab, hartree_partition_tab, &
                                coord_current, current_atom, &
                                current_radial, current_angular&
    )

!  PURPOSE
!   Calculate the integration weight for one point on the grid
!   adapted from: initialize_grid_storage_p1
!
!  USES
!
      use dimensions,only:n_centers_basis_integrals,n_multipoles,n_pp_atoms,n_atoms
      use grids,only:invert_log_grid,r_grid_min,r_grid_inc,n_grid,w_angular,&
                     w_radial
      use runtime_choices,only:flag_hartree_partition_type,partition_type
      use geometry,only:multipole_coords,pp_coords,empty,species
      use spline,only:val_spline
      use free_atoms,only:hartree_partition_rho_spl,partition_rho_spl,free_rho_spl
      use pbc_lists,only:coords_center,centers_basis_integrals,species_center,&
                         center_to_atom
      use species_data,only:species_pseudoized,multipole_radius_free_sq
      use localorb_io,only:localorb_info,OL_low,OL_high
      implicit none

!  ARGUMENTS

      real*8 :: partition_tab
      real*8 :: dummy
      real*8 :: hartree_partition_tab
      real*8, dimension(3) :: coord_current
      integer :: current_atom, current_radial, current_angular
! INPUTS
! o coord_current - cartesian coordinate
! o current_atom - atom
! o current_radial - radial coordinate
! o current_angular - angular 
! OUTPUTS
! o partition_tab -- grid integration weight for one point 
! o hartree_partition_tab -- grid integration weight for hartree potential for one point 
!
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
!

      ! kept only for convenience, will not hurt production, may simplify
      ! future debugging efforts:
      integer :: write_out = 0 
!     local variables
      real*8 :: dir_current(3)
      real*8 :: dist_current, dist_current_sq, dens_current, r_temp

      real*8 :: dist_tab_sq(n_centers_basis_integrals)
      real*8 :: dir_tab(3,n_centers_basis_integrals)

      logical :: point_on_atom, points_on_multipole, point_on_pseudocore
      real*8 ::  dist_to_multipole_sq, dist_to_pseudocore_sq
      integer :: i_pp_atom,i_center_L,i_multipole, i_center

      integer :: n_compute_atoms
      integer :: n_compute_occ_atoms
      integer, dimension(n_centers_basis_integrals) :: center_index

      real*8 :: dist_tab(n_centers_basis_integrals)
      real*8 :: i_r(n_centers_basis_integrals)
      real*8 :: dir_tab_norm(3,n_centers_basis_integrals)

      integer :: i_center_L2
      real*8 :: aux_spl
      real*8, dimension(n_centers_basis_integrals) :: temp_free_rho

      integer, dimension(n_centers_basis_integrals) :: i_occ2i_compute
      real*8, dimension(n_centers_basis_integrals) :: temp_occ_rho

      real*8 :: hartree_partition_deriv
      real*8, dimension(3,n_atoms) :: partition_deriv_delley




          ! remember where this point is with respect to its own atom as well as its free atom density
          dir_current(:)  = coord_current(:)- &
                            coords_center(:,centers_basis_integrals(current_atom))
          dist_current_sq = dir_current(1)*dir_current(1) &
                          + dir_current(2)*dir_current(2) &
                          + dir_current(3)*dir_current(3)
          dist_current    = sqrt(dist_current_sq)
          dir_current(:)  = dir_current(:)/dist_current
          r_temp          = invert_log_grid &
                            ( dist_current, &
                            r_grid_min(species_center(centers_basis_integrals(current_atom))), &
                            r_grid_inc(species_center(centers_basis_integrals(current_atom))))
          dens_current    =  val_spline &
                ( r_temp, hartree_partition_rho_spl(1,1,species_center(centers_basis_integrals(current_atom))), &
                  n_grid(species_center(centers_basis_integrals(current_atom))) )

!              tabulate current integration point as it appears from spherical
!              coordinates centered at each atom

          call tab_atom_centered_coords_p0 &
              ( coord_current, &
              dist_tab_sq, &
              dir_tab, &
              n_centers_basis_integrals, centers_basis_integrals )

          ! To avoid a floating-point exception, check here whether this integration point happens to sit on
          ! another atom.
          ! VB: We could also check here whether we are inside the innermost logarithmic
          !     grid shell of an atom, as is done further down in evaluate_partition_tab.
          !     This would be even better, since we could then remove the same check from
          !     evaluate_partition altogether ...
          point_on_atom = .false.
          do i_center_L = 1, n_centers_basis_integrals, 1
            if ( dist_tab_sq(i_center_L).eq.0.d0) then
              point_on_atom = .true.
              exit ! exit do loop
            end if
          end do

          point_on_pseudocore = .false.
          do i_pp_atom = 1, n_pp_atoms, 1
            dist_to_pseudocore_sq = &
               (coord_current(1) - pp_coords(1, i_pp_atom))**2 &
              +(coord_current(2) - pp_coords(2, i_pp_atom))**2 &
              +(coord_current(3) - pp_coords(3, i_pp_atom))**2 
            if(dist_to_pseudocore_sq.eq.0.d0) then
              point_on_pseudocore = .true.
              exit ! exit do loop
            end if
          end do

          ! also check for multipole singularities
          points_on_multipole = .false.
          do i_multipole = 1, n_multipoles
            dist_to_multipole_sq = &
               (coord_current(1) - multipole_coords(1, i_multipole))**2 &
              +(coord_current(2) - multipole_coords(2, i_multipole))**2 &
              +(coord_current(3) - multipole_coords(3, i_multipole))**2 
            if(dist_to_multipole_sq.eq.0.d0) then
              points_on_multipole = .true.
              exit
            endif
          enddo

          if (.not. point_on_atom .and. .not.points_on_multipole .and. .not.point_on_pseudocore) then
            ! This is the normal case. If our grid point does not sit on an atom, we keep it for later use.

            ! For the following operations, we only need those atoms that have a non-zero
            ! free-atom charge density at the current point. We assemble that list here explicitly ...
            n_compute_atoms = 0
            n_compute_occ_atoms = 0

            do i_center_L = 1, n_centers_basis_integrals, 1
              i_center = centers_basis_integrals(i_center_L)

              if ((i_center.eq.current_atom).or. &
                  (dist_tab_sq(i_center_L).lt.multipole_radius_free_sq(species_center(i_center)) )) then
                  ! this center has a non-zero free-atom density, or it belongs to current integration point
                  n_compute_atoms                  = n_compute_atoms + 1
                  center_index(n_compute_atoms)    = i_center
                  dist_tab_sq(n_compute_atoms)     = dist_tab_sq(i_center_L)
                  dir_tab(:,n_compute_atoms)       = dir_tab(:,i_center_L)   ! This ensures consistent handling later, note n_compute_atoms <= i_center_L
                  atom_atom_index_sym(n_compute_atoms) = i_center_L              ! indexing for later use of atom_atom_tab_sym, which is NOT recomputed here for speed reasons
              end if

            end do


            call tab_global_geometry_p0 &
                ( dist_tab_sq,         &
                  dir_tab,             &
                  dist_tab,            &
                  i_r,                 &
                  dir_tab_norm,        &
                  n_compute_atoms,     &
                  center_index )

            ! calculate the free-atom density only for the (now) known atoms ...
            do i_center_L2 = 1, n_compute_atoms
              i_center = center_index(i_center_L2)
              if (i_center.eq.current_atom) i_center_L = i_center_L2 ! remember the center we are currently at!
              aux_spl = val_spline &
                  ( i_r(i_center_L2), hartree_partition_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )

              temp_free_rho(i_center_L2) = aux_spl
!DB101813
               if ((.not.empty(center_to_atom(i_center))).and.&
                  (.not.(species_pseudoized(species(center_to_atom(i_center)))))) then
!              if (.not.empty(center_to_atom(i_center))) then
                  n_compute_occ_atoms = n_compute_occ_atoms + 1
                  i_occ2i_compute(n_compute_occ_atoms) = i_center_L2
                  temp_occ_rho(n_compute_occ_atoms) = aux_spl
              end if

            end do


!               evaluate Hartree partition tab first
            call evaluate_partition_p2 &
            ( flag_hartree_partition_type,      &
              current_atom,                     &
              i_center_L,                       &
              dist_current,                     &
              dir_current,                      &
              dist_current_sq,                  &
              dens_current,                     &
              dist_tab,                         &
              i_r,                              &
              dir_tab_norm,                     &
              w_angular(current_angular, current_radial, species(current_atom)),  &
              hartree_partition_tab,   &
              dummy, &
              hartree_partition_deriv, &
              n_centers_basis_integrals,        &
              n_compute_atoms,                  &
              center_index,                     &
              temp_free_rho,                    &
              dist_tab_sq,                      &
              atom_atom_tab_sym,                    &
              atom_atom_index_sym,                  &
              min_atom_atom_tab_sym,                &
              n_atom_atom_tab_sym,                  &
              atom_atom_dist_list_sym,                 &
              atom_idx_A_sym,                       &
              atom_idx_B_sym,                       &
              write_out, & 
              partition_deriv_delley(1:3,1:n_atoms) )  
!           evaluate partition_tab
            if (partition_type.ne.flag_hartree_partition_type) then
              ! Partition functions for integration and Hartree potential are different

              ! re-evaluate pieces for the partition tab to be on the safe side
              if ( (partition_type .ne. 5) .and. (partition_type .ne. 7) &
                   .and. (partition_type .ne. 8) .and. (partition_type .ne. 9) ) then ! no density needed for stratmann
                dens_current    =  val_spline &
                  ( r_temp, partition_rho_spl(1,1,species_center(centers_basis_integrals(current_atom))), &
                    n_grid(species_center(centers_basis_integrals(current_atom))) )
                do i_center_L2 = 1, n_compute_atoms
                  i_center = center_index(i_center_L2)
                  if (i_center.eq.current_atom) i_center_L = i_center_L2 ! remember the center we are currently at!
                  aux_spl = val_spline &
                  ( i_r(i_center_L2), hartree_partition_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )

                  temp_free_rho(i_center_L2) = aux_spl

                end do

              end if


  
              call evaluate_partition_tab_p2  &
                ( current_atom,              &
                  i_center_L,                &
                  dist_current,              &
                  dist_current_sq,           &
                  dens_current,              &
                  dist_tab,                  &
                  i_r,                       &
                  w_radial(current_radial, species(current_atom)),  &
                  w_angular(current_angular, current_radial, species(current_atom)),  &
                  partition_tab,    &
                  dummy, &
                  n_centers_basis_integrals, &
                  n_compute_atoms,           &
                  center_index,              &
                  temp_free_rho,             &
                  dist_tab_sq,               &
                  atom_atom_tab_sym,             &
                  atom_atom_index_sym,           &
                  min_atom_atom_tab_sym,         &
                  n_atom_atom_tab_sym,           &
                  atom_atom_dist_list_sym,          &
                  atom_idx_A_sym,                &
                  atom_idx_B_sym)


            else
              ! Partition function for integration is the same as for Hartree potential
              ! and just needs some integration weights ...

! test
            write (80,'(2X,F15.8,2X,F15.8,2X,F15.8,2X,F15.8,2X, I3,A)')  coord_current(1), coord_current(2), coord_current(3), &
              hartree_partition_tab,current_atom, " coord, partition_tab"


              partition_tab = hartree_partition_tab * &
                w_radial(current_radial, species(current_atom)) * &
                dist_current_sq


            write (81,'(2X,F15.8,2X,F15.8,2X,F15.8,2X,F15.8,2X, I3,A)')  coord_current(1), coord_current(2), coord_current(3), &
              partition_tab,current_atom, " coord, partition_tab"

            end if


            ! Fallout from the partition_tab fix rho_r2_lda: More complexity
            if ( flag_hartree_partition_type.eq.6 ) then
              ! must retabulate densities of non-empty atoms for the initial electron density here;
              ! the previous tabulation was LDA, possibly not the correct XC.
              n_compute_occ_atoms = 0
              do i_center_L2 = 1, n_compute_atoms
                i_center = center_index(i_center_L2)
! DB 092812:  we should also exclude pseudocores at this place.
!             won't change the result, but should give a minimal speed-up.
               if ((.not.empty(center_to_atom(i_center))).and.&
                  (.not.(species_pseudoized(species(center_to_atom(i_center)))))) then

!                if (.not.empty(center_to_atom(i_center))) then

                  aux_spl = val_spline &
                  ( i_r(i_center_L2), free_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )

                  n_compute_occ_atoms = n_compute_occ_atoms + 1
                  i_occ2i_compute(n_compute_occ_atoms) = i_center_L2
                  temp_occ_rho(n_compute_occ_atoms) = aux_spl
                end if

              end do
            end if


          else
          ! this is the case where the current integration point sits exactly on top of another atom
          ! must never be considered in later dealings!

            hartree_partition_tab = 0.d0
            partition_tab = 0.d0

          end if







    end subroutine local_partition_tab
!******
!****s* FHI-aims/sym_base/count_points
!  NAME
!   count_points
!  SYNOPSIS
subroutine count_points(n_points_in_grid)

!  PURPOSE
!
!  Counts the number of points in the grid
!
!  USES
!
      use dimensions,only:n_atoms
      use grids,only:n_angular,n_radial
      use geometry,only:species

      implicit none

!  ARGUMENTS
      integer :: n_points_in_grid
! INPUTS
!   none
! OUTPUTS
!   o n_points_in_grid - total number of grid points
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
!

  integer :: i_radial,i_atom

  n_points_in_grid = 0
  do i_atom = 1, n_atoms, 1
     do i_radial = 1, n_radial(species(i_atom)), 1
        n_points_in_grid = n_points_in_grid + &
             n_angular( i_radial,species(i_atom) )
     end do
  end do

end subroutine count_points 
!****s* FHI-aims/sym_base/map_real_space
!  NAME
!   map_real_space
!  SYNOPSIS
subroutine map_real_space(num_symmetry,spg_rotations,spg_shift,n_points_in_grid)

!  PURPOSE
!
!  Get full partion tabs for to irreducible points (other points are zero)
!
!  USES
!
      use dimensions,only:n_atoms,n_periodic,n_max_radial,n_max_angular
      use grids,only:n_angular,n_radial,r_angular,r_radial
      use runtime_choices
      use geometry,only:species,lattice_vector, frac2cart, cart2frac
      use pbc_lists,only:coords_center
      use mpi_tasks, only: check_allocation
      use localorb_io,only:localorb_info,use_unit
      use synchronize_mpi_basic,only:sync_vector
      implicit none

!  ARGUMENTS
  integer, intent(in) :: num_symmetry
  integer, intent(in) :: spg_rotations(3,3,num_symmetry)
  real*8, intent(in)  :: spg_shift(3,num_symmetry)
  integer, intent(in) :: n_points_in_grid
! INPUTS
! o num_symmetry - number of elements in space group
! o spg_rotations - Rotations
! o spg_shift - Translations
! o n_points_in_grid - total number of grid points
! OUTPUTS
! o partition_tab_sym(i_atom,i_radial,i_angulal) - partion tab for irreducible points
! o hartree_partition_tab_sym(i_atom,i_radial,i_angulal)- partion tab for irreducible points
! o map_real(i_point) - cartesian coordinate and multiplicity of coordinate (reduced)
! o map_grid - radial grid index for reduced grid
!
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
!

  integer :: n_real_grid
  integer :: i_radial,i_atom,i_angular,is,i_grid
  integer :: info
  real*8, dimension(3) :: coord_current, coord_current_temp,coord_current_rot,&
                          trans, coord_current_rot_scaled
  real*8  :: partition_tab_local,hartree_partition_tab_local
  character*140 :: info_str

  ! Allocations
  if (.not.allocated(map_real))then
    allocate ( map_real(4,n_points_in_grid) ,stat=info)
      call check_allocation(info, 'map_real                       ')
  endif
  if (.not.allocated(map_grid))then
    allocate ( map_grid(3,n_points_in_grid) ,stat=info)
      call check_allocation(info, 'map_grid                       ')
  endif
  if (.not.allocated(partition_tab_sym))then
    allocate ( partition_tab_sym(n_atoms,n_max_radial,n_max_angular) ,stat=info)
      call check_allocation(info, 'partition_tab_sym               ')
  endif
  if (.not.allocated(hartree_partition_tab_sym))then
    allocate ( hartree_partition_tab_sym(n_atoms,n_max_radial,n_max_angular) ,stat=info)
      call check_allocation(info, 'hartree_partition_tab_sym               ')
  endif
  
  ! initialize counters
  n_real_grid = 0
  map_real = 0d0
  map_grid = 0

  ! initialize partition tab calculations
  call prepare_local_partition_tab()
  
  ! initialize partition tabs
  partition_tab_sym = 0d0
  hartree_partition_tab_sym = 0d0
  
  ! loop over grid 
  do i_atom = 1, n_atoms, 1
    do i_radial = 1, n_radial(species(i_atom)), 1
       do i_angular = 1, n_angular( i_radial,species(i_atom) )

          n_real_grid = n_real_grid + 1
          map_real(4,n_real_grid) = map_real(4,n_real_grid) +1
          
          ! cartesian coordinate
          coord_current(:) = coords_center( :,i_atom ) + &
                             r_angular(:, i_angular, i_radial, species(i_atom)) * &
                             r_radial(i_radial, species(i_atom))
          coord_current_temp = coord_current

          ! center unit cell
          if(n_periodic > 0)then
             call map_to_center_cell(coord_current_temp(1:3) )
          end if
          !call local_partition_tab( partition_tab_local, hartree_partition_tab_local, &
          !                       coord_current_temp, i_atom, &
          !                     i_radial, i_angular )
          !partition_tab_sym(i_atom,i_radial, i_angular) = partition_tab_local
          !hartree_partition_tab_sym(i_atom,i_radial, i_angular) = hartree_partition_tab_local
          
          ! Set original values for mapping array and partition tabs
          map_real(1:3,n_real_grid) = coord_current_temp(1:3)
          map_grid(1,n_real_grid) = i_atom
          map_grid(2,n_real_grid) = i_radial
          map_grid(3,n_real_grid) = i_angular
          
     
          call local_partition_tab( partition_tab_local, hartree_partition_tab_local, &
                               coord_current, i_atom, i_radial, i_angular )
          partition_tab_sym(i_atom,i_radial,i_angular) = partition_tab_local
          hartree_partition_tab_sym(i_atom,i_radial,i_angular) = hartree_partition_tab_local
          
          ! Loop over symmetry operations (spg_rotations,spg_shift)
          do is = 2, num_symmetry, 1
          
             ! rotate current coordinate 
             call frac2cart(lattice_vector,spg_shift(1:3,is),trans)
             coord_current_rot = coord_current_temp - trans
             !coord_current_rot = matmul(transpose(spg_rotations(1:3,1:3,is)),coord_current_rot)
             call cart2frac(lattice_vector, coord_current_rot, coord_current_rot_scaled)
             coord_current_rot_scaled = matmul(inv(dble(spg_rotations(1:3,1:3,is))),coord_current_rot_scaled)
             call frac2cart(lattice_vector, coord_current_rot_scaled, coord_current_rot)
             if(n_periodic > 0)then
               call map_to_center_cell(coord_current_rot(1:3) )
             end if

             ! Find match for rotated coordinate in grid
             do i_grid = 1, n_real_grid-1, 1  
                          
               if (sqrt(sum((coord_current_rot-map_real(1:3,i_grid))**2))<1e-6)then
                  partition_tab_sym(i_atom,i_radial,i_angular) = partition_tab_sym(i_atom,i_radial,i_angular)- partition_tab_local
                  hartree_partition_tab_sym(i_atom,i_radial,i_angular) = hartree_partition_tab_sym(i_atom,i_radial,i_angular)-hartree_partition_tab_local         
                  map_real(4,n_real_grid) = map_real(4,n_real_grid) - 1
                  map_real(1:3,n_real_grid) = map_real(1:3,i_grid)
                  map_real(4,i_grid) =  map_real(4,i_grid) +1
                  map_grid(1,n_real_grid) = map_grid(1,i_grid)
                  map_grid(2,n_real_grid) = map_grid(2,i_grid)
                  map_grid(3,n_real_grid) = map_grid(3,i_grid)

                  partition_tab_sym(map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid)) =&
                  partition_tab_sym(map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid)) + partition_tab_local
                  hartree_partition_tab_sym(map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid)) =&
                  hartree_partition_tab_sym(map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid)) + hartree_partition_tab_local
                  exit
               endif
             end do
             if (sqrt(sum((coord_current_rot-map_real(1:3,i_grid))**2))<1e-4)then
                exit
             endif
          end do

       end do
     end do
  end do
  !        write(use_unit,*) sum(partition_tab_sym)
  write(info_str, '(2X,A,I8,A,I8)') &
  & '| points reduced from: ', n_real_grid, ' to ', count(map_real(4,:) > 0)
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str, '(2X,A,I8,A,I8)') &
  & '| nonzero points in partition_tab_sym: ', count(partition_tab_sym > 0), '  hartree_partition_tab_sym', count(hartree_partition_tab_sym > 0)
  call localorb_info(info_str,use_unit,'(A)')
  deallocate(map_real)
  !deallocate(map_grid)
 call clean_local_partition_tab()
end subroutine map_real_space
!****s* FHI-aims/sym_base/set_partition_tab_sym
!  NAME
!   set_partition_tab_sym
!  SYNOPSIS
subroutine set_partition_tab_sym(partition_tab_std,hartree_partition_tab_std)

!  PURPOSE
!  Set the the local partition_tab for all points in batch

! To use this add before line 600 to initialize_grid_storage_p1.f90
!            partition_tab(i_point) = partition_tab_sym(current_atom,current_radial,current_angular)
!            hartree_partition_tab(i_point) = hartree_partition_tab_sym(current_atom,current_radial,current_angular)


!  USES
  use dimensions,only:n_full_points,n_my_batches
  use grids, only: batches,batch_of_points
  use synchronize_mpi_basic,only:sync_integer
  use localorb_io,only:localorb_info,use_unit

  implicit none
!  ARGUMENTS
  real*8, target, dimension(n_full_points) :: partition_tab_std
  real*8, target, dimension(n_full_points) :: hartree_partition_tab_std


!  INPUTS

!  OUTPUT
!    o partition_tab -- the partition tab
!    o hartree_partition_tab -- the hartree_partition tab
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


  integer :: current_atom, current_radial, current_angular

  integer :: num_hartree_partition,num_partition
  character*140 :: info_str

  real*8, pointer :: partition_tab(:)
  real*8, pointer :: hartree_partition_tab(:)
  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    hartree_partition_tab => hartree_partition_tab_std
  partition_tab = 0d0
  hartree_partition_tab = 0d0
  i_point = 0
  do i_my_batch = 1, n_my_batches_work, 1
     

        ! loop over one batch
        do i_index = 1, batches_work(i_my_batch)%size, 1

           
              i_point = i_point+1

              ! get current integration point coordinate
              current_atom    = batches_work(i_my_batch) % points(i_index) % index_atom
              current_radial  = batches_work(i_my_batch) % points(i_index) % index_radial
              current_angular = batches_work(i_my_batch) % points(i_index) % index_angular
              partition_tab(i_point) = partition_tab_sym(current_atom,current_radial,current_angular)
              hartree_partition_tab(i_point) =  hartree_partition_tab_sym(current_atom,current_radial,current_angular)
        enddo  ! end loop over one part of the angular division              

  end do ! end loop over batches

  num_partition = count(partition_tab > 0)
  num_hartree_partition = count(hartree_partition_tab > 0)

  call sync_integer(num_partition)
  call sync_integer(num_hartree_partition)

  write(info_str, '(2X,A,I8,A,I8)') &
  & '| nonzero points in partition_tab: ', num_partition, '  hartree_partition_tab', num_hartree_partition
  call localorb_info(info_str,use_unit,'(A)')



end subroutine set_partition_tab_sym
!****s* FHI-aims/sym_base/out_ks_coeff
!  NAME
!   out_ks_coeff
!  SYNOPSIS
subroutine out_ks_coeff(KS_eigenvector_complex)

!  PURPOSE
!    Print KS eigenvector 
!  USES
  use dimensions,only:n_k_points,n_basis,n_states,n_spin,n_k_points_task
  use localorb_io,only:localorb_info,localorb_allinfo
  use mpi_tasks,only:myid, n_tasks
  implicit none
!  ARGUMENTS
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task), intent(IN):: KS_eigenvector_complex
!  INPUTS
!    o KS_eigenvector,KS_eigenvector_complex
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
! SOURCE  
  integer :: i_k, i_k_point,  i_basis, i_state
  character*120 :: info_str
  open(11,file='KS_org.dat')
       i_k = 0
       do i_k_point = 1, n_k_points,1
          if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
             i_k = i_k + 1
             !write(info_str,'(A,E12.4,E12.4,E12.4 )') 'k_pointlist', k_point_list(i_k_point,1:3)
             !call localorb_allinfo(info_str)
             !if(i_k_point.eq.n_k_points)then
               do i_basis = 1, n_basis
                   i_state = 14
                   write(11,'(I10,I10,I10,I10,2F20.14,1F20.14)') i_k_point,  i_k_point, i_basis, 1, &
                       & KS_eigenvector_complex(i_basis, i_state,1,i_k),abs(KS_eigenvector_complex(i_basis, i_state,1,i_k)) 
               enddo
             !endif
          end if !myid == mod(k...
       end do ! i_k_point
  close(11)
endsubroutine out_ks_coeff
!****s* FHI-aims/sym_base/reduced_density_to_full_density
!  NAME
!   reduced_density_to_full_density
!  SYNOPSIS
subroutine reduced_density_to_full_density(rho)

!  PURPOSE
!    Transfrom the symmetry reduced density (grid) back to the full grid
!    -Firt we map from points on batches to radial coordinates
!    -Sync
!    -Transfromation
!    -Distribute on batches
!  USES
  use dimensions,only:n_atoms,n_max_radial,n_max_angular,n_full_points,n_my_batches
  use localorb_io,only:localorb_info
  use mpi_tasks,only:check_allocation
  use geometry, only: species
  use grids,only: batches,n_angular,n_radial
  use synchronize_mpi_basic,only:sync_vector
  implicit none
!  ARGUMENTS
  real*8,     dimension(n_full_points), intent(INOUT):: rho
!  INPUTS
!    o rho - density
!  OUTPUT
!    o rho - density
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
  integer :: i_my_batch,i_index,i_point
  integer :: info
  integer :: i_radial,i_atom,i_angular,i_grid
  integer :: current_atom,current_radial,current_angular
  real*8, dimension(:,:,:), allocatable :: rho_transfrom
  if (.not.allocated(rho_transfrom))then
    allocate ( rho_transfrom(n_atoms,n_max_radial,n_max_angular) ,stat=info)
      call check_allocation(info, 'rho_transfrom               ')
  endif

  rho_transfrom = 0d0

  i_point = 0

  do i_my_batch = 1, n_my_batches, 1

     do i_index = 1, batches(i_my_batch)%size, 1

        i_point = i_point + 1

        current_atom    = batches(i_my_batch) % points(i_index) % index_atom
        current_radial  = batches(i_my_batch) % points(i_index) % index_radial
        current_angular = batches(i_my_batch) % points(i_index) % index_angular

        rho_transfrom(current_atom,current_radial,current_angular) = rho(i_point)

     enddo
  enddo

  call sync_vector(rho_transfrom,n_atoms*n_max_radial*n_max_angular)

  i_grid = 0

  do i_atom = 1, n_atoms, 1
    do i_radial = 1, n_radial(species(i_atom)), 1
       do i_angular = 1, n_angular( i_radial,species(i_atom) )
         i_grid = i_grid + 1
         rho_transfrom(i_atom,i_radial,i_angular) = rho_transfrom(map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid))
        !if(myid==0)then
        !write(use_unit,*) i_atom, i_radial, i_angular, map_grid(1,i_grid),map_grid(2,i_grid),map_grid(3,i_grid)
        !endif
       enddo
    enddo
  enddo

  i_point = 0

  do i_my_batch = 1, n_my_batches, 1

     do i_index = 1, batches(i_my_batch)%size, 1

        i_point = i_point + 1

        current_atom    = batches(i_my_batch) % points(i_index) % index_atom
        current_radial  = batches(i_my_batch) % points(i_index) % index_radial
        current_angular = batches(i_my_batch) % points(i_index) % index_angular

        rho(i_point) = rho_transfrom(current_atom,current_radial,current_angular)
        !if(myid==1.and.current_atom==2.and.current_radial==42)then
        !write(use_unit,*) current_atom, current_radial, current_angular, rho(i_point),rho_transfrom(current_atom,current_radial,current_angular)
        !endif
     enddo
  enddo
if (allocated(rho_transfrom)) deallocate(rho_transfrom)
endsubroutine reduced_density_to_full_density
!****s* FHI-aims/sym_base/evaluate_k_densmat_sym
!  NAME
!   evaluate_k_densmat_sym
!  SYNOPSIS
  subroutine evaluate_k_densmat_sym(kdensmat, kdensmat_complex, occ_numbers, &
  &                             KS_eigenvector, KS_eigenvector_complex, &
  &                             i_spin, i_k_point, i_k)

    !  PURPOSE
    !    Calculate the k-dependend density matrix
    !    for one k and one spin only for real_eigenvectors.
    !    -> Adopted from accumulate_k_densmat in density_matrix_evaluation.f90
    !       For Symmetry reconstruction, reconstructed KS-coeff as Input
    !  USES
    use spglib_symmetry,only:map, map_sym
    use dimensions,only:n_spin,n_states,n_basis,n_k_points_task, n_k_points_nosym
    use runtime_choices,only:real_eigenvectors,PM_index,&
                             use_spg_full_Delta
    use mpi_tasks,only: check_allocation, aims_stop
    implicit none

    !  ARGUMENTS


    real*8, intent(OUT) :: kdensmat(n_basis, n_basis)
    complex*16, intent(OUT) :: kdensmat_complex(n_basis, n_basis)
    real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points_nosym)
    ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
    !        properly k-weighted (and are thus not the "true" occupation numbers.)
    !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
    real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)
    integer, intent(IN) :: i_spin, i_k_point, i_k

    !  INPUTS
    !    o occ_numbers -- occupation of states, with k-weighting already applied
    !    o KS_eigenvector, KS_eigenvector_complex -- eigencoefficients
    !    o i_spin -- spin component
    !    o i_k_point -- k-point index
    !    o i_k -- node-local k-point index
    !  OUTPUTS
    !    o kdensmat, kdensmat_complex -- density matrix of this k-point
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_state, n_max_occupied, info, i_bas, is, new_k_point, i_basis
    real*8 :: occu
    real*8, allocatable :: KS_scaled(:,:)
    complex*16, allocatable :: KS_scaled_cmplx(:,:)
    character(*), parameter :: func = 'evaluate_k_densmat'

    if (real_eigenvectors) then
       kdensmat = 0.d0
    else
       kdensmat_complex = (0.d0, 0.d0)
    end if
  
    is= rotations_at_task_map(map_sym(i_k_point))
    new_k_point = map(i_k_point)

    if (count(occ_numbers(:, i_spin, i_k_point) > 0d0) > n_states/2) then
       ! Most states are occupied; use collective operations.
       if (real_eigenvectors) then
          allocate(KS_scaled(n_basis, n_states), stat=info)
       else
          allocate(KS_scaled_cmplx(n_basis, n_states), stat=info)
       end if
       call check_allocation(info, 'evaluate_k_densmat:KS_scaled{_cmplx}')

       n_max_occupied = 1
       do i_state = 1, n_states, 1
          occu =  occ_numbers(i_state, i_spin, i_k_point)
          if (occu .gt. 0.d0) then
             n_max_occupied = i_state
          endif
          if (real_eigenvectors) then
            if(map_sym(i_k_point).eq.1)then
              KS_scaled(:, i_state) = KS_eigenvector(:, i_state, i_spin, i_k) * sqrt(occu)
            else
              if(use_spg_full_Delta)then
                call dgemv('N', n_basis,  n_basis,sqrt(occu), dble(Delta_matrix_full(:,:,is,i_k)), n_basis, &
            &          KS_eigenvector(:, i_state, i_spin, i_k), 1, 0.d0,&
                      KS_scaled(:, i_state), 1 )
              else
                if(real_eigenvectors)then
                  call aims_stop('On the fly calculation of Delta matrix for irr. k-point set not implemented.')
                else
                  call get_kvec_sym(rotations_at_task(is), new_k_point, map_atom(:,is), &
                    Delta(:,:,:,is), KS_eigenvector_complex(:, i_state, i_spin, i_k), &
                    KS_scaled_cmplx(:, i_state),i_state,i_k_point)
                  !debug if(i_state==14)then
                  !debug    do i_basis=1,n_basis
                  !debug      write(10,'(I10,I10,I10,I10,2F20.14,1F20.14)') i_k_point, new_k_point, i_basis, map_sym(i_k_point),&
                  !debug      & KS_scaled_cmplx(i_basis, i_state),abs(KS_scaled_cmplx(i_basis, i_state)) 
                  !debug    enddo
                  !debug  endif
                   KS_scaled(:, i_state)=sqrt(occu)*KS_scaled_cmplx(:, i_state)
                endif
              endif
           endif
          else
            if(map_sym(i_k_point).eq.1)then
               KS_scaled_cmplx(:, i_state) = sqrt(occu)*KS_eigenvector_complex(:, i_state, i_spin, i_k)
               !debug    if(i_state==14)then
               !debug       do i_basis=1,n_basis
               !debug         write(10,'(I10,I10,I10,I10,2F20.14,1F20.14)') i_k_point, new_k_point, i_basis, map_sym(i_k_point),&
               !debug         & KS_eigenvector_complex(i_basis, i_state, i_spin, i_k),abs(KS_eigenvector_complex(i_basis, i_state, i_spin, i_k)) 
               !debug       enddo 
               !debug     endif
                   
            else
              if(use_spg_full_Delta)then
                call zgemv('N', n_basis,  n_basis, dcmplx(sqrt(occu),0.d0), Delta_matrix_full(:,:,is,i_k), n_basis, &
            &          KS_eigenvector_complex(:, i_state, i_spin, i_k), 1, dcmplx(0.d0,0.d0),&
                      KS_scaled_cmplx(:, i_state), 1 )
              else
                call get_kvec_sym(rotations_at_task(is), new_k_point, map_atom(:,is), &
                  Delta(:,:,:,is), dcmplx(1.d0,0.d0)*KS_eigenvector_complex(:, i_state, i_spin, i_k), &
                  KS_scaled_cmplx(:, i_state),i_state,i_k_point)
                  !debug if(i_state==14)then
                  !debug    do i_basis=1,n_basis
                  !debug      write(10,'(I10,I10,I10,I10,2F20.14,1F20.14)') i_k_point, new_k_point, i_basis, map_sym(i_k_point),&
                  !debug      & KS_scaled_cmplx(i_basis, i_state),abs(KS_scaled_cmplx(i_basis, i_state)) 
                  !debug    enddo
                  !debug  endif
                   KS_scaled_cmplx(:, i_state)=dcmplx(sqrt(occu),0.d0)*KS_scaled_cmplx(:, i_state)
              endif
            endif
          end if
       end do
       if (real_eigenvectors) then
          call dsyrk('U', 'N', n_basis, n_max_occupied, 1.d0, KS_scaled, n_basis, &
          &          0.d0, kdensmat, n_basis)
       else
          call zherk('U', 'N', n_basis, n_max_occupied, 1.d0, KS_scaled_cmplx, n_basis, &
          &          0.d0, kdensmat_complex, n_basis)
       end if
    else
       ! Most states are unoccupied; use selective operations.
       do i_state = 1, n_states
          occu =  occ_numbers(i_state, i_spin, i_k_point)
          if (occu .gt. 0.d0) then
             if (real_eigenvectors) then
                if(map_sym(i_k_point).eq.1)then
                  call dsyr('U', n_basis, occu, &
                  &         KS_eigenvector(:, i_state, i_spin, i_k), 1, kdensmat, n_basis)
                else             
                  if(use_spg_full_Delta)then
                    call dgemv('N', n_basis,  n_basis, 1.d0, dble(Delta_matrix_full(:,:,is,i_k)), n_basis, &
                    &          KS_eigenvector(:, i_state, i_spin, i_k), 1, 0.d0,&
                            KS_scaled(:, i_state), 1 )
                  else
                    if(real_eigenvectors)then
                     call aims_stop('On the fly calculation of Delta matrix for irr. k-point set not implemented.')
                    else                 
                      call get_kvec_sym(rotations_at_task(is), new_k_point, map_atom(:,is), &
                        Delta(:,:,:,is), KS_eigenvector_complex(:, i_state, i_spin, i_k), &
                        KS_scaled_cmplx(:, i_state),i_state,i_k_point)
                    endif
                  endif
                  call dsyr('U', n_basis, occu, &
                  &         KS_scaled(:, i_state), 1, &
                  &         kdensmat, n_basis)
                endif
             else
                if(map_sym(i_k_point).eq.1)then
                  call zher('U', n_basis, occu, &
                  &         KS_eigenvector_complex(:, i_state, i_spin, i_k), 1, &
                  &         kdensmat_complex, n_basis)
                else             
                  if(use_spg_full_Delta)then
                    call zgemv('N', n_basis,  n_basis, dcmplx(1.d0,0.d0), Delta_matrix_full(:,:,is,i_k), n_basis, &
                    &          KS_eigenvector_complex(:, i_state, i_spin, i_k), 1, dcmplx(0.d0,0.d0),&
                            KS_scaled_cmplx(:, i_state), 1 )
                  else
                    call get_kvec_sym(rotations_at_task(is), new_k_point, map_atom(:,is), &
                        Delta(:,:,:,is), KS_eigenvector_complex(:, i_state, i_spin, i_k), &
                        KS_scaled_cmplx(:, i_state),i_state,i_k_point)
                  endif
                  call zher('U', n_basis, occu, &
                  &         KS_scaled_cmplx(:, i_state), 1, &
                  &         kdensmat_complex, n_basis)
                endif
             end if
          end if
       end do
    end if

    if (real_eigenvectors) then
       kdensmat = kdensmat + transpose(kdensmat)
       do i_bas = 1, n_basis
          kdensmat(i_bas, i_bas) = kdensmat(i_bas, i_bas)*0.5d0
       end do
    else
       kdensmat_complex = kdensmat_complex + transpose(dconjg(kdensmat_complex))
       do i_bas = 1, n_basis
          kdensmat_complex(i_bas, i_bas) = kdensmat_complex(i_bas, i_bas)*0.5d0
       end do
    end if

   if(allocated(KS_scaled)) deallocate(KS_scaled)
   if(allocated(KS_scaled_cmplx)) deallocate(KS_scaled_cmplx)

  end subroutine evaluate_k_densmat_sym
!****s* FHI-aims/sym_base/accumulate_k_densmat_sym
!  NAME
!   accumulate_k_densmat_sym
!  SYNOPSIS
subroutine accumulate_k_densmat_sym(density_matrix_sparse, density_matrix, force_packed, &
  &                               kdm, kdm_complex, k_phase)

    !  PURPOSE
    !    Accumulate the k-dependent density matrices to the real-space representation.
    !    -> Adopted from accumulate_k_densmat in density_matrix_evaluation.f90
    !       For Symmetry reconstruction full k-phase is required
    !  USES
    use dimensions,only:n_basis,n_hamiltonian_matrix_size,n_centers_basis_T
    use pbc_lists,only:index_hamiltonian, column_index_hamiltonian,&
                       n_cells_in_hamiltonian,Cbasis_to_basis,center_to_cell,&
                       Cbasis_to_center, n_cells
    use runtime_choices,only:real_eigenvectors,packed_matrix_format,PM_index,&
                             PM_none
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: density_matrix_sparse(n_hamiltonian_matrix_size)
    real*8, intent(INOUT) :: density_matrix(n_centers_basis_T, n_centers_basis_T)
    logical, intent(IN) :: force_packed
    real*8, intent(IN) :: kdm(n_basis, n_basis)
    complex*16, intent(IN) :: kdm_complex(n_basis, n_basis)
    complex*16, intent(IN) :: k_phase(n_cells)
    !  INPUTS
    !    o density_matrix{_sparse} -- real-space density matrix
    !    o force_packed -- use (lapack-)packed density_matrix even for PM_none
    !    o kdm{_cmplx} -- k-dependent density matrix for i_k_point
    !  OUTPUTS
    !    o density_matrix{_sparse} -- real-space density matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    ! Add these matrix elements with their phases to the corresponding
    ! real space matrix.

    real*8 :: add
    integer :: i_bas1, i_bas2, i_basT1, i_basT2
    integer :: i_cell1, i_cell2
    integer :: i_cell
    integer :: i_index
    complex*16 :: conj_phase
    character(*), parameter :: func = 'accumulate_k_densmat'

    select case(packed_matrix_format)
    case(PM_index)
       do i_bas2 = 1, n_basis
          do i_cell = 1,n_cells_in_hamiltonian-1
             conj_phase = conjg(k_phase(i_cell))
             if (index_hamiltonian(1,i_cell, i_bas2) > 0) then
                do i_index = index_hamiltonian(1, i_cell, i_bas2), &
                &            index_hamiltonian(2, i_cell, i_bas2)
                   i_bas1 =  column_index_hamiltonian(i_index)

                   if (real_eigenvectors) then
                      density_matrix_sparse(i_index) &
                      & = density_matrix_sparse(i_index) &
                      & + kdm(i_bas1, i_bas2) * dble(conj_phase)
                   else

                      density_matrix_sparse(i_index) &
                      & = density_matrix_sparse(i_index) &
                      & + dble(kdm_complex(i_bas1, i_bas2) * conj_phase)
                      
                   end if

                end do
             end if
          end do
       end do

    case(PM_none)

       i_index = 0
       do i_basT2 = 1,  n_centers_basis_T
          i_bas2 = Cbasis_to_basis(i_basT2)
          i_cell2 = center_to_cell(Cbasis_to_center(i_basT2))
          do i_basT1 = 1, n_centers_basis_T
             if (force_packed .and. i_basT1 > i_basT2) cycle
             i_bas1 = Cbasis_to_basis(i_basT1)
             i_cell1 = center_to_cell(Cbasis_to_center(i_basT1))
             i_index = i_index + 1

             if (real_eigenvectors) then
                add = kdm(i_bas1, i_bas2) &
                &     * dble(k_phase(i_cell1)) * dble(k_phase(i_cell2))
             else
                add = dble(kdm_complex(i_bas1, i_bas2) &
                &     * dconjg(k_phase(i_cell1)) * k_phase(i_cell2))
             end if
             if (force_packed) then
                density_matrix_sparse(i_index) = density_matrix_sparse(i_index) + add
             else
                density_matrix(i_basT1, i_basT2) = density_matrix(i_basT1, i_basT2) + add
             end if
          end do
       end do

    case default
       call aims_stop('Invalid packing type', func)

    end select

end subroutine accumulate_k_densmat_sym
!****s* FHI-aims/sym_base/reset_map
!  NAME
!   reset_map
!  SYNOPSIS
subroutine reset_map()

    !  PURPOSE
    !    Reset the mapping between old and new k-points
    !    mostly taken from pbc_list.f90
    !  USES

    use spglib_symmetry,only: map, num_symmetry, map_sym, map_inv
    use runtime_choices,only: k_points_offset, reconstruct_proper_only, &
                              n_k_points_xyz_nosym, use_k_phase_exx
    use dimensions, only: n_k_points,n_k_points_task, n_periodic, use_hf_kspace, &
                          use_hf_realspace, use_full_spectrum, n_k_points_nosym
    use mpi_tasks,only:myid, n_tasks, check_allocation
    use pbc_lists,only:cell_index, n_cells,n_cells_bvk
    use constants, only: pi
    use localorb_io,only:localorb_info
    implicit none

    !  ARGUMENTS

    real*8 :: r_x,r_y,r_z
    integer :: i_k, i_x, i_y, i_z, i_k_point, old_k_point, new_k_point, i_cell,&
               nr, is, iss, i_k_point_inv,n_k_full,i,max_k_number,i_cell_n,&
               i_cell_1,i_cell_2,i_cell_3
    integer :: info
    logical:: found 
    integer,dimension(:),allocatable :: map_sym_no_sym
    integer,dimension(:),allocatable :: map_sym_inv_sym
    integer,dimension(:),allocatable :: map_inv_sym

    integer,dimension(n_k_points_xyz_nosym(1),n_k_points_xyz_nosym(2),n_k_points_xyz_nosym(3)) :: k_number_inv
    
    n_k_points_nosym=product(n_k_points_xyz_nosym)
    if (.not.allocated(map_sym_no_sym))then
      allocate ( map_sym_no_sym(n_k_points_nosym) ,stat=info)
        call check_allocation(info, 'map_sym_no_sym             ')
    endif
    if (.not.allocated(map_sym_inv_sym))then
      allocate ( map_sym_inv_sym(n_k_points_nosym) ,stat=info)
        call check_allocation(info, 'map_sym_inv_sym             ')
    endif
    old_k_point = 0
    i_k_point_inv = 0
    i_k_point = 0
    map_sym_no_sym = 0
    map_sym_inv_sym = 0
    do i_x = 1, n_k_points_xyz_nosym(1)
      do i_y = 1, n_k_points_xyz_nosym(2)
          do i_z = 1, n_k_points_xyz_nosym(3)
            old_k_point = old_k_point + 1
            k_number_inv(i_x,i_y,i_z)= count(map_inv == old_k_point)
            if(count(map_inv == old_k_point)> 0)then
                i_k_point_inv = i_k_point_inv + 1
                map_sym_inv_sym(old_k_point) = i_k_point_inv
            endif
            if(count(map == old_k_point)> 0)then
                i_k_point = i_k_point + 1
                map_sym_no_sym(old_k_point) = i_k_point
            endif
          enddo
      enddo
    enddo
      i_k_point_inv = 0
      do old_k_point = 1,n_k_points_nosym
        if (map_sym_inv_sym(old_k_point).ne.0) then
          i_k_point_inv = i_k_point_inv+1
        endif
      enddo   
    ! map k-number without symmetry to k-number with symmetry
    if (reconstruct_proper_only) then
      if (allocated(map_inv)) deallocate(map_inv)
      allocate ( map_inv(count(k_number_inv.ne.0)) ,stat=info)
      call check_allocation(info, 'map_inv             ')
      allocate ( map_inv_sym(count(k_number_inv.ne.0)) ,stat=info)
      call check_allocation(info, 'map_inv_sym          ')
      i_k_point_inv = 0
      do old_k_point = 1,n_k_points_nosym
        if (map_sym_inv_sym(old_k_point).ne.0) then
          i_k_point_inv = i_k_point_inv+1
          map_inv(i_k_point_inv) = map_sym_no_sym(map(old_k_point))
          map_inv_sym(i_k_point_inv) = map_sym(old_k_point)
        endif
      enddo

      n_k_points_nosym = count(k_number_inv.ne.0)
      if (allocated(map)) deallocate(map)
      allocate ( map(n_k_points_nosym) ,stat=info)
      call check_allocation(info, 'map             ')
      if (allocated(map_sym)) deallocate(map_sym)
      allocate ( map_sym(n_k_points_nosym) ,stat=info)
      call check_allocation(info, 'map_sym             ')
      do old_k_point = 1,n_k_points_nosym
        map(old_k_point) = map_inv(old_k_point)
        map_sym(old_k_point) = map_inv_sym(old_k_point)
      enddo
    else
      k_number_inv = 1
      i_k=0
      do old_k_point = 1,n_k_points_nosym
        map(old_k_point) = map_sym_no_sym(map(old_k_point))
      enddo
    endif
    
    max_k_number = 0
    do i_k_point = 1,n_k_points
        max_k_number= max( max_k_number,count(map==i_k_point))
    enddo
    if (.not.allocated(k_point_list_nosym))then
      allocate ( k_point_list_nosym(n_k_points_nosym,3) ,stat=info)
        call check_allocation(info, 'k_point_list_nosym             ')
    endif

    if (.not.allocated( k_points_at_task))then
      allocate (  k_points_at_task(n_k_points) ,stat=info)
        call check_allocation(info, ' k_points_at_task             ')
    endif

    if (.not.allocated( k_phase_base_nosym))then
      allocate (  k_phase_base_nosym(3,n_k_points_nosym) ,stat=info)
        call check_allocation(info, ' k_phase_base_nosym            ')
    endif
    
    if (.not.allocated( k_phase_nosym))then
      allocate (  k_phase_nosym(n_cells,n_k_points_nosym) ,stat=info)
        call check_allocation(info, ' k_phase_nosym             ')
    endif
    !write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells*n_k_points_nosym*16)/2**10, ' KiB for k_phase_nosym.'
    if (.not.allocated( k_weights_nosym))then
      allocate (  k_weights_nosym(n_k_points_nosym) ,stat=info)
        call check_allocation(info, ' k_weights_nosym             ')
    endif

    if (.not.allocated( rotations_at_task))then
      allocate (  rotations_at_task(num_symmetry) ,stat=info)
        call check_allocation(info, ' rotations_at_task            ')
    endif
    if(use_hf_kspace)then
      if (.not.allocated( k_per_irr))then
	allocate (  k_per_irr(n_k_points,max_k_number) ,stat=info)
	  call check_allocation(info, ' k_per_irr            ')
      endif
      if (.not.allocated( n_k_per_irr))then
	allocate (  n_k_per_irr(n_k_points) ,stat=info)
	  call check_allocation(info, ' n_k_per_irr            ')
      endif    
    endif    
      if (.not.allocated( k_per_task))then
	allocate (  k_per_task(n_k_points_task,max_k_number) ,stat=info)
	  call check_allocation(info, ' k_per_task            ')
      endif          
      if (.not.allocated( n_k_per_task))then
	allocate (  n_k_per_task(n_k_points_task) ,stat=info)
	  call check_allocation(info, ' n_k_per_task            ')
      endif
    if((use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt. 0)).and.use_k_phase_exx)then
      if (.not.allocated( k_phase_exx_nosym))then
        allocate (  k_phase_exx_nosym(n_cells_bvk,n_k_points_nosym) ,stat=info)
          call check_allocation(info, ' k_phase_nosym             ')
      endif
      k_phase_exx_nosym = (1.d0,0.d0)
    endif      
      if (.not.allocated(map_back))then
          allocate ( map_back(n_k_points) ,stat=info)
          call check_allocation(info, 'map_back          ')
      endif

    ! We need to reconstruc k-point arrays without symmetry

    ! First where are my k-points at
    i_k = 0
    k_points_at_task = 0
    k_per_task = 0
    n_k_per_task=0
    n_k_points_task_full = 0
    do new_k_point = 1,n_k_points
      if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
          i_k = i_k+1
          k_points_at_task(new_k_point) = i_k
	    i=0
	    do n_k_full = 1, n_k_points_nosym, 1
	      if(map(n_k_full).eq.new_k_point)then
		i=i+1
		k_per_task(i_k,i) = n_k_full      
	      endif
	    enddo
	    n_k_per_task(i_k) = i
	    n_k_points_task_full = n_k_points_task_full + i
      endif
      if(use_hf_kspace)then
	i=0
	do n_k_full = 1, n_k_points_nosym, 1
	  if(map(n_k_full).eq.new_k_point)then
	    i=i+1
	    k_per_irr(new_k_point,i) = n_k_full      
	  endif
	enddo
	n_k_per_irr(new_k_point) = i   
      endif
    enddo
   !debug  write(use_unit,*) '--------------------------------------------'
   !debug  write(use_unit,*) '--------------------------------------------'
    !debug do i=1, n_k_points_task
   !debug    write(use_unit,'(I10)')n_k_per_task(i)
   !debug    do new_k_point = 1,n_k_per_task(i)
   !debug      write(use_unit,'(I10,I10,I10)') i, new_k_point, k_per_task(i,new_k_point)
   !debug    enddo
   !debug  enddo
   !debug  write(use_unit,*) '--------------------------------------------'
   !debug  write(use_unit,*) '--------------------------------------------'
    !debug write(use_unit,*) n_k_points_xyz_nosym
    ! original k_point_list, k_phase_base, map_sym
    i_k_point = 0
    old_k_point = 0
    do i_x = 1, n_k_points_xyz_nosym(1)
      do i_y = 1, n_k_points_xyz_nosym(2)
          do i_z = 1, n_k_points_xyz_nosym(3)
            if(k_number_inv(i_x,i_y,i_z).gt.0)then
            old_k_point = old_k_point + 1
            r_x = dble(i_x-1) / dble(n_k_points_xyz_nosym(1)) + k_points_offset(1)
            r_y = dble(i_y-1) / dble(n_k_points_xyz_nosym(2)) + k_points_offset(2)
            r_z = dble(i_z-1) / dble(n_k_points_xyz_nosym(3)) + k_points_offset(3)
            k_point_list_nosym(old_k_point,1) = r_x
            k_point_list_nosym(old_k_point,2) = r_y
            k_point_list_nosym(old_k_point,3) = r_z
            !write(use_unit,'(A,I4,A,I4,I4,I4)') 'i_k: ', old_k_point, ', k-coord: ', (i_x-1), (i_y-1), (i_z-1)
            k_phase_base_nosym(1, old_k_point) = exp((0.d0,1.d0) * 2*pi * r_x)
            k_phase_base_nosym(2, old_k_point) = exp((0.d0,1.d0) * 2*pi * r_y)
            k_phase_base_nosym(3, old_k_point) = exp((0.d0,1.d0) * 2*pi * r_z)
            k_weights_nosym(old_k_point) = dble(k_number_inv(i_x,i_y,i_z))
            
            if((use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt. 0)).and.use_k_phase_exx)then
              i_cell_n = 1
              do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
                do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                   do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2
                                     
                     if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then
                       i_cell_n = i_cell_n + 1
                       k_phase_exx_nosym(i_cell_n, old_k_point) &
                       & = k_phase_base_nosym(1, old_k_point)**i_cell_1 &
                       & * k_phase_base_nosym(2, old_k_point)**i_cell_2 &
                       & * k_phase_base_nosym(3, old_k_point)**i_cell_3
                     endif
                   enddo
                 enddo
               enddo
            endif
            endif
          enddo
      enddo
    enddo
    k_weights_nosym(:) = k_weights_nosym(:) / ( n_k_points_xyz_nosym(1)* n_k_points_xyz_nosym(2)* n_k_points_xyz_nosym(3))

    ! k_phase without symmetry
    do i_cell = 1, n_cells-1
       do i_k_point = 1, n_k_points_nosym
          k_phase_nosym(i_cell, i_k_point) &
          & = product(k_phase_base_nosym(:, i_k_point)**cell_index(i_cell,:))
       enddo
    end do
!do i_k_point = 1, n_k_points_nosym
!do i_cell = 1,n_cells-1
!write(use_unit,'(I5,I5,1F14.8,1F14.8)') i_k_point, i_cell,dble(k_phase_nosym(i_cell,i_k_point)),dimag(k_phase_nosym(i_cell,i_k_point))
!enddo
!enddo
    ! The last cell is not set in the above loops since it is a dummy cell
    k_phase_nosym(n_cells,:) = 0

    if(allocated(map_sym_no_sym)) deallocate(map_sym_no_sym)
    if(allocated(map_sym_inv_sym)) deallocate(map_sym_inv_sym)
    if(allocated(map_inv)) deallocate(map_inv)
    if(allocated(map_inv_sym)) deallocate(map_inv_sym)


    n_rotations_at_task = 0
    rotations_at_task = 0
    !do is = 1, num_symmetry, 1
      i_k = 0
      do new_k_point = 1,n_k_points
        if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
            i_k = i_k+1            
            do i_k_point = 1, n_k_points_nosym, 1
            found=.false.
            if (map(i_k_point)==new_k_point)then
              do nr= 1, n_rotations_at_task, 1
                if (rotations_at_task(nr)==map_sym(i_k_point))then
                  found=.true.
                  exit
                endif
              enddo
              if (.not.found)then
                n_rotations_at_task=n_rotations_at_task+1
                rotations_at_task(n_rotations_at_task)=map_sym(i_k_point)
              endif
            endif
            enddo
        endif
      enddo  
    !enddo

   if (.not.allocated( rotations_at_task_map))then
      allocate (  rotations_at_task_map(num_symmetry) ,stat=info)
        call check_allocation(info, ' rotations_at_task_map         ')
   endif

   rotations_at_task_map=1
   do iss = 1, n_rotations_at_task, 1
     do is = 1, num_symmetry, 1
       if(rotations_at_task(iss)==is)then
         rotations_at_task_map(is)=iss
         exit
       endif
     enddo
   enddo
   
   if((use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt. 0)).and.use_k_phase_exx)then  
     i_k=0
     do old_k_point=1, n_k_points_nosym
        if(map_sym(old_k_point).eq.1)then
          i_k=i_k+1
          map_back(i_k)=old_k_point
        endif
      enddo
    endif
    !debug write(use_unit,'(I5,I5)') myid, n_rotations_at_task
    !debug do is = 1, n_rotations_at_task, 1
    !debug   write(use_unit,'(I5,I5,I5,I5)') myid, is, rotations_at_task(is), rotations_at_task_map(rotations_at_task(is))
    !debug enddo
end subroutine reset_map
!****s* FHI-aims/sym_base/get_kvec_sym
!  NAME
!   get_kvec_sym
!  SYNOPSIS
subroutine get_kvec_sym(is,new_k_point, map_atom, Delta, KS_eigenvector_complex, KS_new,i_state,i_k_point)

    !  PURPOSE
    !    Reconstruct the full KS-coefficents (i_k_point) from reduced set (new_k_point)
    !  USES
    use spglib_symmetry, only: spg_rotations
    use dimensions, only: n_basis, l_wave_max, n_atoms
    use geometry, only: frac_coords, species
    use pbc_lists, only: Cbasis_to_atom, k_point_list
    use basis,only: basis_m, basis_l,basis_fn
    use localorb_io,only:localorb_info,localorb_allinfo,use_unit
    use constants, only: pi
    implicit none
    !  INPUTS
    !    o new_k_point - reduced k-point
    !    o KS_eigenvector{_complex} - Reduced set of KS-coefficents
    !  OUTPUTS
    !    o KS_new - Full KS-coefficents
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    integer, intent(IN) :: is,i_state,i_k_point
    integer, intent(IN) :: new_k_point
    integer, dimension(n_atoms), intent(IN) :: map_atom
    complex*16,intent(IN), dimension(0:l_wave_max,-l_wave_max:l_wave_max,-l_wave_max:l_wave_max) :: Delta
    complex*16,intent(IN), dimension(n_basis) :: KS_eigenvector_complex
    complex*16,intent(OUT), dimension(n_basis) :: KS_new

    integer :: i_species, i_basis, i_atom, i_fn, i_basis_l, i_basis_m,&
               j_species, j_basis, j_atom, j_fn, j_basis_l, j_basis_m

    real*8 ::  arg1, arg2
    complex*16 :: phase
    real*8 ::  k_point_rot(3)

    logical :: debug

    ! debug
    debug = .false.



    KS_new = (0.d0,0.d0)
    k_point_rot = matmul(transpose(inv(dble(spg_rotations(:, :, is)))),k_point_list(new_k_point,:))

    do i_basis = 1, n_basis
      i_atom = Cbasis_to_atom( i_basis )
      i_species = species(i_atom)
      i_fn = basis_fn(i_basis)
      i_basis_l = basis_l(i_basis)
      i_basis_m = basis_m(i_basis)

      do j_basis = 1, n_basis
        j_atom = Cbasis_to_atom( j_basis )
        j_species = species(j_atom)
        j_fn = basis_fn(j_basis)
        j_basis_l = basis_l(j_basis)
        j_basis_m = basis_m(j_basis)
        
        ! Transformationmatrix of YLM: Delta
        if (j_atom.eq.map_atom(i_atom).and.i_species.eq.j_species.and.i_basis_l.eq.j_basis_l.and.i_fn.eq.j_fn) then 

          ! Phase due to shift between basis functions
          arg1 = dot_product (k_point_list(new_k_point,:), frac_coords(:,j_atom))
          arg2 = dot_product (k_point_rot(:), frac_coords(:,i_atom))         
          phase = exp (dcmplx (0.d0, 2.0d0*pi*(arg2-arg1)))
          KS_new(i_basis)=KS_new(i_basis)+&
          phase*Delta(i_basis_l,i_basis_m,j_basis_m)*KS_eigenvector_complex(j_basis)
          ! debug if(i_k_point==3.and.i_state==14)then
          ! debug    write(use_unit,'(I5,I5,I5,I5,I5,2F20.14,2F20.14,2F20.14,2F20.14)') new_k_point,i_basis,&
          ! debug      & j_basis,j_atom, i_atom, KS_eigenvector_complex(j_basis),KS_new(i_basis),&
          ! debug     &  Delta(i_basis_l,i_basis_m,j_basis_m),phase
          ! debug endif
        endif
      enddo

    enddo


    if (debug) then
      do i_basis = 1, n_basis
          write(use_unit,'(A,I4,I4,I4,2F11.4)') 'i_k_point, is, i_basis, KS_new: ', new_k_point, is, i_basis, KS_new(i_basis)      
      enddo
    endif

end subroutine get_kvec_sym
!****s* FHI-aims/sym_base/get_density_rotation
!  NAME
!   get_density_rotation
!  SYNOPSIS
subroutine get_density_rotation(is, new_k_point, map_atom, Delta, Delta_matrix)

    !  PURPOSE
    !    Reconstruct the full KS-coefficents (i_k_point) from reduced set (new_k_point)
    !  USES
    use spglib_symmetry, only: spg_rotations
    use dimensions, only: n_basis, l_wave_max, n_atoms
    use geometry, only: frac_coords, species
    use pbc_lists, only: Cbasis_to_atom, k_point_list
    use basis,only: basis_m, basis_l,basis_fn
    use localorb_io,only:localorb_info,localorb_allinfo,use_unit
    use constants, only: pi
    implicit none
    !  INPUTS
    !    o new_k_point - reduced k-point
    !    o i_k_point - full k-point
    !    o KS_eigenvector{_complex} - Reduced set of KS-coefficents
    !  OUTPUTS
    !    o KS_new - Full KS-coefficents
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    integer, intent(IN) :: is, new_k_point
    integer, dimension(n_atoms), intent(IN) :: map_atom
    complex*16,intent(IN), dimension(0:l_wave_max,-l_wave_max:l_wave_max,-l_wave_max:l_wave_max) :: Delta
    complex*16,intent(OUT), dimension(n_basis,n_basis) :: Delta_matrix

    integer :: i_species, i_basis, i_atom, i_fn, i_basis_l, i_basis_m,&
               j_species, j_basis, j_atom, j_fn, j_basis_l, j_basis_m

    real*8 ::  arg1, arg2
    complex*16 :: phase
    real*8 ::  k_point_rot(3)

    logical :: debug

    ! debug
    debug = .false.



    Delta_matrix = 0.d0
    k_point_rot = matmul(transpose(inv(dble(spg_rotations(:, :, is)))),k_point_list(new_k_point,:))

    do i_basis = 1, n_basis
      i_atom = Cbasis_to_atom( i_basis )
      i_species = species(i_atom)
      i_fn = basis_fn(i_basis)
      i_basis_l = basis_l(i_basis)
      i_basis_m = basis_m(i_basis)

      do j_basis = 1, n_basis
        j_atom = Cbasis_to_atom( j_basis )
        j_species = species(j_atom)
        j_fn = basis_fn(j_basis)
        j_basis_l = basis_l(j_basis)
        j_basis_m = basis_m(j_basis)
        
        ! Transformationmatrix of YLM: Delta
        if (map_atom(i_atom).eq.j_atom.and.i_species.eq.j_species.and.i_basis_l.eq.j_basis_l.and.i_fn.eq.j_fn) then 

          ! Phase due to shift between basis functions
          arg1 = dot_product (k_point_list(new_k_point,:), frac_coords(:,j_atom))
          arg2 = dot_product (k_point_rot(:), frac_coords(:,i_atom))         
          phase = exp (dcmplx (0.d0, 2.0d0*pi*(arg2-arg1)))
          Delta_matrix(i_basis,j_basis)=phase*Delta(i_basis_l,i_basis_m,j_basis_m)

        endif
      enddo

    enddo


    if (debug) then
      do i_basis = 1, n_basis
        do j_basis = 1, n_basis
          write(use_unit,'(A,I4,I4,I4,I4,2F20.14)') 'i_k_point, is, i_basis, j_basis, Delta_matrix: ', new_k_point, is, i_basis, j_basis, Delta_matrix(i_basis,j_basis)      
        enddo
      enddo
    endif

end subroutine get_density_rotation
!****s* FHI-aims/sym_base/evaluate_densmat_sym
!  NAME
!   evaluate_densmat_sym
!  SYNOPSIS
subroutine evaluate_densmat_sym &
  ( KS_eigenvector, KS_eigenvector_complex, occ_numbers,  &
  density_matrix, density_matrix_sparse, i_spin, force_packed &
  ) 
    !  PURPOSE
    !    Adapted from evaluate_densmat (density_matrix_evaluation.f90)
    !    Evaluates the density matrix for symmetry reduced k-points
    !    
    !     -KS-coeffs are reconstructed from reduced set
    !     -k_phase and occ_numbers are adapted
    !  USES
    use spglib_symmetry, only: map, map_sym
    use localorb_io,only:localorb_info,use_unit, OL_norm
    use runtime_choices,only:real_eigenvectors,packed_matrix_format,PM_index,&
                             PM_none,use_scalapack,use_local_index,&
                             use_spg_full_Delta,use_spg_mv_mm
    use dimensions,only:n_k_points,n_spin,n_basis,n_states,n_k_points_task,&
                        n_k_points_nosym, n_hamiltonian_matrix_size,n_centers_basis_T
    use mpi_tasks,only:myid, n_tasks, check_allocation, aims_stop
    use load_balancing,only: batch_perm, use_batch_permutation, get_full_local_matrix
    use synchronize_mpi,only:sync_density_matrix_sparse, sync_density_matrix
    use scalapack_wrapper, only: construct_dm_scalapack, &
                                 get_sparse_matrix_scalapack, get_full_matrix_scalapack,&
                                 get_sparse_local_matrix_scalapack
    use density_matrix_evaluation, only:evaluate_k_densmat, accumulate_k_densmat
    use pbc_lists, only: k_weights, n_cells
    implicit none

    !  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
    ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
    !        properly k-weighted (and are thus not the "true" occupation numbers.)
    !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
    integer, intent(IN) :: i_spin
    real*8, dimension(n_centers_basis_T,n_centers_basis_T), intent(OUT) :: density_matrix
    ! when this routine is called, density_matrix_sparse has either the dimension
    ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
    ! so we declare it here as a 1D assumed size array
    real*8, intent(OUT) :: density_matrix_sparse(*)
    logical, intent(IN) :: force_packed

    !  INPUTS
    !   o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !   o occ_numbers -- occupation of states, with k-weighting already applied
    !   o i_spin -- spin index
    !   o force_packed -- use packed storage (Lapack-type in ..._sparse) even if
    !                   packed_matrix_format is PM_none
    !
    !  OUTPUT
    !   o density_matrix -- density matrix if non-packed matrix is in use
    !   o density_matrix_sparse -- density matrix if packed matrix is in use
    ! 
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
    !

    !     other local variables
    real*8, dimension(:,:), allocatable :: kdm
    complex*16, dimension(:,:), allocatable:: kdm_complex
    integer :: i_k_point, i_k, n, i ,j, k_task
    integer :: info
    character(*), parameter :: func = 'evaluate_densmat_sum'

    integer :: new_k_point

    ! Reconstruced KS-coefficent
    logical:: old
    complex*16,dimension(n_cells):: kp
    call localorb_info("Evaluating symmetry reduced density matrix", use_unit,'(2X,A)', OL_norm )


    !if (.not.allocated( KS_eigenvalue_nosym))then
    !  allocate (  KS_eigenvalue_nosym(n_states,n_spin,n_k_points_nosym) ,stat=info)
    !        call check_allocation(info, ' KS_eigenvalue_nosym      ')
    !endif
    !KS_eigenvalue_nosym = 0.d0
    !do i_k_point = 1, n_k_points_nosym, 1
    !  KS_eigenvalue_nosym(1:n_states,1:n_spin,i_k_point)=KS_eigenvalue(1:n_states,1:n_spin,map(i_k_point))
      !write(use_unit,'(I5,I5,1F14.8,1F14.8)') i_k_point,map(i_k_point),  KS_eigenvalue(6,1,map(i_k_point)),KS_eigenvalue_nosym(6,1,i_k_point)
    !enddo
    if(use_spg_mv_mm)then
      if (.not.allocated( occ_numbers_nosym))then
        allocate (  occ_numbers_nosym(n_states,n_spin,n_k_points_nosym) ,stat=info)
        call check_allocation(info, 'occ_numbers_nosym             ')
      endif
      occ_numbers_nosym = 0.d0
      do i_k = 1, n_k_points_nosym
        new_k_point = map(i_k)
        occ_numbers_nosym(:,:, i_k) = (occ_numbers(:,:,new_k_point)/k_weights(new_k_point)) * k_weights_nosym(i_k)
      enddo
    else
      if (.not.allocated( occ_numbers_sym))then
        allocate (  occ_numbers_sym(n_states,n_spin,n_k_points) ,stat=info)
        call check_allocation(info, 'occ_numbers_sym            ')
      endif
      do i_k = 1, n_k_points   
        occ_numbers_sym(:,:, i_k) = (occ_numbers(:,:,i_k)/k_weights(i_k))
      enddo
    endif
    !if (allocated(KS_eigenvalue_nosym)) deallocate(KS_eigenvalue_nosym)


    if(packed_matrix_format /= PM_none)then
       if(use_batch_permutation > 0) then
          n = batch_perm(use_batch_permutation)%n_local_matrix_size
          density_matrix_sparse(1:n) = 0.d0     
       else
          density_matrix_sparse(1:n_hamiltonian_matrix_size) = 0.d0     
       endif
    else if (force_packed) then
       density_matrix_sparse(1:n_hamiltonian_matrix_size) = 0.d0     
    else
       density_matrix = 0.d0
    end if

    if(use_scalapack)then

       call construct_dm_scalapack(occ_numbers_nosym, i_spin)

       select case (packed_matrix_format)
       case(PM_index)
          if(use_local_index) then
             if(use_batch_permutation > 0) then
                call get_full_local_matrix(density_matrix_sparse, i_spin)
             else
                call get_sparse_local_matrix_scalapack(density_matrix_sparse, i_spin)
             endif
          else
             call get_sparse_matrix_scalapack(density_matrix_sparse,i_spin)
             call sync_density_matrix_sparse(density_matrix_sparse)
          endif

       case(PM_none)
          if (force_packed) then
             call aims_stop('lapack-packed matrix not for scalapack', func)
          end if
          call get_full_matrix_scalapack( density_matrix, i_spin )
          call sync_density_matrix(density_matrix)           
       case default
          call aims_stop('Invalid packing!', func)
       end select

    else

       !------------ construct density matrix -------------------------------------------

       if (real_eigenvectors) then
          allocate(kdm(n_basis, n_basis), stat=info)
          call check_allocation(info, 'evaluate_densmat:kdm')
          ! dummy for the real case
          allocate(kdm_complex(1, 1), stat=info)
          call check_allocation(info, 'evaluate_densmat:kdm_complex')
       else
          allocate(kdm_complex(n_basis, n_basis), stat=info)
          call check_allocation(info, 'evaluate_densmat:kdm_complex')
          ! dummy for the complex case
          allocate(kdm(1, 1), stat=info)
          call check_allocation(info, 'evaluate_densmat:kdm')   
       end if  
       ! JW: If there is need for optimization, we could just calculate
       ! those density matrix elements which are actually needed.
       ! This would scale O(N^2) instead O(N^3).
  
 
       if(use_spg_mv_mm)then
          ! Here the density matrix is constructed from the reconstructed KS-coefficents
          do i_k_point = 1, n_k_points_nosym, 1
            ! We map the k-point to the reduced set of k-points
            new_k_point = map(i_k_point)
            if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
              i_k = k_points_at_task(new_k_point)
                  ! k-point stored on current mpi task
                  ! Build density (with adapted occupation numbers)
                  call evaluate_k_densmat_sym(kdm, kdm_complex, occ_numbers_nosym, &
                  &                       KS_eigenvector, KS_eigenvector_complex, &
                  &                       i_spin,i_k_point,i_k)
                  call accumulate_k_densmat_sym(density_matrix_sparse, density_matrix, force_packed, &
                  &                         kdm, kdm_complex, k_phase_nosym(:,i_k_point))
            endif
          end do ! i_k_point 
        else
          i_k = 0
          do i_k_point = 1,n_k_points
            if (myid ==  modulo(i_k_point, n_tasks) .and. myid <= n_k_points) then
                i_k = i_k+1 
                k_task = k_per_task(i_k,1)
                kp = k_phase_nosym(:,k_task)*dcmplx(k_weights_nosym(k_task),0.d0)
                call evaluate_k_densmat(kdm, kdm_complex, occ_numbers_sym, &
                      &                       KS_eigenvector, KS_eigenvector_complex, &
                      &                       i_spin,i_k_point,i_k)
                call accumulate_k_densmat_sym(density_matrix_sparse, density_matrix, &
                                force_packed, kdm, kdm_complex,kp)
                call get_reducible_KS_dens(kdm, kdm_complex, density_matrix_sparse, &
                             density_matrix,  force_packed, i_k, i_k_point)
            endif
          enddo! i_k_point
       endif
       
       if(packed_matrix_format /= PM_none .or. force_packed) then
          call sync_density_matrix_sparse(density_matrix_sparse)
       else
          call sync_density_matrix(density_matrix)
       end if

       !write(use_unit,*) '************zero_order_dm in density_matrix_evaluation.f90'
       !write(use_unit,'(40f20.15)') (density_matrix_sparse(i_k),i_k=1,n_hamiltonian_matrix_size-1)
       !debug do i = 1, n_hamiltonian_matrix_size, 1
       !debug   write(use_unit,'(A,I5,2F11.8)') 'DM:', i, density_matrix_sparse(i) 
       !debug enddo
       if (allocated(kdm)) deallocate(kdm)
       if (allocated(kdm_complex)) deallocate(kdm_complex)
       if (allocated(occ_numbers_sym)) deallocate(occ_numbers_sym)
       if (allocated(occ_numbers_nosym)) deallocate(occ_numbers_nosym)

    end if ! use_scalapack
end subroutine evaluate_densmat_sym

!****s* FHI-aims/sym_base/factn
!  NAME
!   factn
!  SYNOPSIS
Real (8) function factn (n)
!  PURPOSE
!
! 
! Returns the factorial n!, n should be less than 150.
!
!  USES
  use mpi_tasks,only: aims_stop
  implicit none
! ARGUMENTS
!   n : input (in,integer)

      integer, intent (IN) :: n
!  INPUTS
!  o n -- n
!
!  OUTPUT
!  o factn -- n!
!  
! ARGUMENTS
      integer :: i
      real (8) :: f1 (24)
      Data f1 / 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0, &
     & 40320.d0, 362880.d0, 3628800.d0, 39916800.d0, 479001600.d0, &
     & 6227020800.d0, 87178291200.d0, 1307674368000.d0, &
     & 20922789888000.d0, 355687428096000.d0, 6402373705728000.d0, &
     & 121645100408832000.d0, 2432902008176640000.d0, &
     & 51090942171709440000.d0, 1124000727777607680000.d0, &
     & 25852016738884976640000.d0, 620448401733239439360000.d0 /
     ! fast return if possible
      if (n .Eq. 0) then
         factn = 1.d0
         return
      end if
      if ((n .Ge. 1) .And. (n .Le. 24)) then
        factn = f1 (n)
        return
      end if
      if (n .Lt. 0) then
         call aims_stop("Error(factnm): n < 0 : ","factn")
      end if
      if (n .Gt. 150) then
         call aims_stop("Error(factnm): n out of range","factn")
      end if
      factn = f1 (24)
      do i = 25, n
        factn = factn * dble (i)
      enddo
      return
end function


!****s* FHI-aims/sym_base/inv
!  NAME
!   inv
!  SYNOPSIS
function inv(A) result(Ainv)
!  PURPOSE
!
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
!
!  USES
  use mpi_tasks,only: aims_stop
  implicit none
! ARGUMENTS
  real*8, dimension(:,:), intent(in) :: A
  real*8, dimension(size(A,1),size(A,2)) :: Ainv
!  INPUTS
!  o A -- Matrix
!
!  OUTPUT
!  o Ainv -- inverted Matrix
!  
! ARGUMENTS
  real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     call aims_stop('Matrix is numerically singular!','inv')
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     call aims_stop('Matrix inversion failed!','inv')
  end if
end function inv


!****s* FHI-aims/sym_base/get_occupation_numbers_sym
!  NAME
!   get_occupation_numbers_sym
!  SYNOPSIS
subroutine get_occupation_numbers_sym(KS_eigenvalue, n_electrons, &
                              t_out, occ_numbers, chemical_potential, k_weights, n_k_points)

!  PURPOSE
! Subroutine get_occupation_numbers for full set of k-points
!
! determines occupation numbers according to two available schemes:
! o 1) gaussian smearing according to Kresse et al., Comp. Mat. Sci. 6, 15 - 50 (1996)
! o 2) fermi smearing
! o 3) methfessel-paxton (methfessel et al., PRB 40, 6, 40 (1989))
!
!  USES

  use dimensions,only:n_states,n_spin
  use runtime_choices,only:fermi_acc,max_zeroin,occupation_acc,occupation_thr
  use localorb_io, only:OL_norm,use_unit,localorb_info
  use constants, only: hartree
  use mpi_tasks,only:check_allocation, aims_stop
  implicit none

!  ARGUMENTS

  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue
  real*8, intent(in) :: n_electrons

  logical :: t_out

  real*8, dimension(n_states, n_spin, n_k_points), intent(out) :: occ_numbers
  real*8, intent(out) :: chemical_potential
  real*8, dimension(n_k_points), intent(in) :: k_weights
  integer, intent(in) :: n_k_points
!  INPUTS
!  o KS_eigenvalue -- Kohn-Sham eigenvalues
!  o n_electrons -- number of electrons
!  o t_out -- is the information printed out or not ?
!
!  OUTPUT
!  o occ_numbers -- occupation weights of different KS states
!  o chemical_potential -- chemical potential
!  
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




  !  local variables
!  real*8 :: degeneracy_threshold
!  integer :: n_degenerate
!  integer :: highest_full_state
!  integer :: lowest_empty_state
!  real*8 :: shared_electrons
!  logical :: degenerate
  real*8 :: chemical_potential_l
  real*8 :: chemical_potential_r
  real*8 :: diff_l
  real*8 :: diff_r
  real*8 :: diff_electrons
  real*8 :: diff_electrons_thr
  real*8 :: lowest_eigenvalue
  real*8 :: highest_eigenvalue

!  real*8 :: electron_count

  character*120 :: info_str

  !  counters
  integer :: i_state
!  integer :: i_state_2
  integer :: i_spin
  integer :: i_k_points
  integer :: i_counter


  if (t_out) then
     write(info_str,'(2X,A)') &
          "Determining occupation numbers for Kohn-Sham eigenstates."
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if



  i_counter = 0
  !  write(use_unit,*) "get_occ_numbers.... n_spin=", n_spin, "n_states=", n_states, "n_k_points=", n_k_points

  lowest_eigenvalue = KS_eigenvalue(1,1,1)

  do i_k_points = 1, n_k_points,1
     do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
           if (KS_eigenvalue(i_state, i_spin,i_k_points) .lt. lowest_eigenvalue) then
              lowest_eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_points)
           end if
        end do
     end do
  end do

  highest_eigenvalue = KS_eigenvalue(n_states,1,1)
  do i_k_points = 1,n_k_points,1
     do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
           if (KS_eigenvalue(i_state, i_spin, i_k_points) .gt. highest_eigenvalue) then
              highest_eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_points)
           end if
        end do
     end do
  end do

  !  initialize zeroin algorithm
  !  to avoid problems if only one state is calculated (so for the H-atom)
  !  do not initialize mit KS_eigenvalue(n_states), because it is then
  !  identical to KS_eigenvalue(1)

  if (lowest_eigenvalue.ne.highest_eigenvalue) then
     chemical_potential_r = highest_eigenvalue
  else
     chemical_potential_r = 0.d0
  end if
  chemical_potential_l = lowest_eigenvalue

  call check_norm_sym(chemical_potential_l, KS_eigenvalue, n_electrons, occ_numbers, diff_l, i_counter, k_weights, n_k_points)
  call check_norm_sym(chemical_potential_r, KS_eigenvalue, n_electrons, occ_numbers, diff_r, i_counter, k_weights, n_k_points)

  do while (diff_l * diff_r .gt. 0.d0)
     !     interval for chemical potential still not found
     !     Must extend intervals both ways: 
     !     * For zero electrons, need Fermi level below lowest state
     !     * For all electrons in one channel, need Fermi level above highest state

     chemical_potential_l = chemical_potential_l - 0.5d0 * dabs(highest_eigenvalue-lowest_eigenvalue)
     chemical_potential_r = chemical_potential_r + 0.5d0 * dabs(highest_eigenvalue-lowest_eigenvalue)
     diff_l = diff_r

     call check_norm_sym(chemical_potential_l, KS_eigenvalue, n_electrons, occ_numbers, diff_l, i_counter, k_weights,n_k_points)
     call check_norm_sym(chemical_potential_r, KS_eigenvalue, n_electrons, occ_numbers, diff_r, i_counter, k_weights,n_k_points)

     !     If the right interval cannot be found, stop potentially
     !     infinite loop

     if (i_counter .gt. max_zeroin) then
        write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') "* Too many iterations needed to find good interval around E_F."
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
             "* According to this check, you may not be able to find one. "
        call localorb_info(info_str,use_unit,'(A)')
        write (info_str,'(1X,A)') &
             "Please look carefully at your settings."
        call localorb_info(info_str,use_unit,'(A)')
        call aims_stop("*** Error in finding electronic chemical potential", &
              "get_occupation_numbers_p0")
     end if
  end do


  call zeroin_sym(chemical_potential_l, chemical_potential_r, diff_l, diff_r, fermi_acc, occupation_acc, max_zeroin, &
       n_states, KS_eigenvalue, n_electrons, chemical_potential, diff_electrons, occ_numbers, i_counter,k_weights,n_k_points)


  !     call check_norm one more time, just to make sure we really have the 
  !     occupation numbers for the final fermi level [this is critical for H]
  call check_norm_sym( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter, k_weights, n_k_points)
  !     and decrement i_counter which was falsely incremented ...
  i_counter = i_counter - 1

  !       If zeroin did not converge to the desired accuracy, issue
  !       warning and stop, for now; we can create nicer behavior after
  !       we gain more experience.
  if (i_counter .gt. max_zeroin) then
     write (info_str,'(1X,A)') "* Determination of electronic chemical potential (Fermi level):"
     call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') "* Too many iterations needed to find correct occupation of levels."
     call localorb_info(info_str,use_unit,'(A)')
     write (info_str,'(1X,A)') "* Check variables occupation_acc, max_zeroin before continuing your calculation."
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* Status of get_occ_numbers when aborted: "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* "
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) "* State #         Eigenvalue       Occ. number"
     call localorb_info(info_str,use_unit,'(A)')
     do i_k_points = 1, n_k_points, 1
        do i_spin = 1, n_spin, 1
           write (info_str,*) "* spin # ", i_spin
           call localorb_info(info_str,use_unit,'(A)')
           write (info_str,*) "----------"
           call localorb_info(info_str,use_unit,'(A)')
           do i_state = 1, n_states, 1
              write (info_str,*) i_state, KS_eigenvalue(i_state, i_spin,i_k_points), occ_numbers(i_state, i_spin,i_k_points)
              call localorb_info(info_str,use_unit,'(A)')
           end do
        end do
     end do
     stop
  end if

  !     finally, occupation thresholding - if needed
  if (occupation_thr.gt.0.d0) then
     call threshold_occ_numbers( n_electrons, occ_numbers, diff_electrons_thr)
  end if

  if (t_out) then
     write (info_str, '(2X, A, E15.8)') "| Chemical potential (Fermi level) in eV                 : ", chemical_potential * hartree
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
     write (info_str, '(2X, A, E15.8)') "| Error in electron count due to remaining E_F inaccuracy: ", diff_electrons
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
     if (occupation_thr.gt.0.d0) then
        write (info_str, '(2X, A, E15.8)') "| Error in electron count after thresholding : ", diff_electrons_thr
        call localorb_info(info_str,use_unit,'(A)',OL_norm)
     end if
  end if

  !DB 04/26/12 .. has to get uncommented again .. however causes segfault
  call check_n_states(n_states, n_spin, n_k_points, occ_numbers, t_out)


  if (t_out) then
     write(info_str,*)
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if

end subroutine get_occupation_numbers_sym
!****s* FHI-aims/check_norm_sym
!  NAME
!    check_norm_p0_sym
!  SYNOPSIS
subroutine check_norm_sym(chemical_potential, KS_eigenvalue, n_electrons, &
                          occ_numbers, diff_electrons, i_counter,&
                          k_weights, n_k_points)

!  PURPOSE
!  checks and returns the charge density norm (i.e. the electron count) for
!  a given Fermi level, and for whichever smearing type was chosen ...
!
!  Thus, this routine also produces updated occupation numbers.
!
!  USES
  use dimensions,only:n_states, n_spin, spin_degeneracy
  use runtime_choices,only:n_methfessel_paxton,occupation_type,occupation_width
  use arch_specific,only:arch_erf
  use constants,only: pisqrt_inv
  use synchronize_mpi,only:sync_integer_vector
  use localorb_io,only:use_unit
  implicit none

!  ARGUMENTS

  !real*8, intent(in) :: chemical_potential
  integer, intent(in) :: n_k_points
  real*8 :: chemical_potential
  real*8, dimension(n_states, n_spin, n_k_points), intent(in) :: KS_eigenvalue
  real*8, intent(in) :: n_electrons
  real*8, dimension(n_states, n_spin,n_k_points), intent(out) :: occ_numbers
  real*8, intent(out) :: diff_electrons
  integer, intent(inout) :: i_counter 
  real*8, dimension(n_k_points), intent(in) :: k_weights


!  INPUTS
!   o chemical_potential -- chemical potential
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o n_electrons -- number of electrons
!    
!  OUTPUT
!   o occ_numbers -- occupations of states
!   o diff_electrons -- difference between formally expected number of electrons
!                       and actual electron count for present trial Fermi level
!   o i_counter -- Number of iterations we have taken so far to find the Fermi
!                  level
! 
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
!



  
  !  local variables
  real*8 :: temp_n_electrons
  real*8 :: one_over_ow
  real*8 :: H_odd
  real*8 :: H_even
  real*8 :: hermite_arg
  real*8 :: gauss_weight
  real*8 :: A
  real*8 :: max_exponent
  real*8 :: exp_arg

  !  counter
  integer :: i_state
  integer :: i_spin
  integer :: i_mp
  integer :: i_k_point



  one_over_ow   = 1.d0 / occupation_width
  !debug write(use_unit,*) "check_norm... n_spin=", n_spin, "n_states=", n_states
  !debug do i_k_point = 1, n_k_points,1
  !debug   do i_spin = 1, n_spin, 1
  !debug      do i_state = 1, n_states, 1
  !debug         write(use_unit,*) i_k_point, i_state, i_spin, KS_eigenvalue(i_state, i_spin, i_k_point)
  !debug         write(use_unit,*) (KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow
  !debug      end do
  !debug   end do
  !debug end do
  i_counter = i_counter + 1
  temp_n_electrons   = 0.d0
  select case(occupation_type)
  case (0)
     !  gaussian smearing (0 <= occ_numbers <= spin_degenracy)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              occ_numbers(i_state, i_spin, i_k_point) = spin_degeneracy *  0.5d0 * & 
                   (1.d0 - arch_erf((KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow))

              temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin, i_k_point) * k_weights(i_k_point)

           end do
        end do
     end do


  case (1)
     !  fermi smearing (0 <= occ_numbers <= spin_degeneracy)
     max_exponent = maxexponent(chemical_potential) * log(2.d0)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              exp_arg = (KS_eigenvalue(i_state, i_spin,i_k_point) - chemical_potential) * one_over_ow
              if (exp_arg .lt. max_exponent) then

                 occ_numbers(i_state, i_spin, i_k_point) = spin_degeneracy / (1.d0 + exp(exp_arg))

                 temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin, i_k_point) * k_weights(i_k_point)

              else
                 occ_numbers(i_state, i_spin, i_k_point) = 0.d0
              end if
           end do
        end do
     end do

  case (2)
     ! Methfessel-Paxton smearing (0 <= occ_numbers <= spin_degeneracy)

     do i_k_point = 1, n_k_points,1
        do i_spin = 1, n_spin, 1
           do i_state = 1, n_states, 1

              hermite_arg  = (KS_eigenvalue(i_state, i_spin, i_k_point) - chemical_potential) * one_over_ow 
              gauss_weight = exp(- hermite_arg * hermite_arg)

              !     zero order contribution
              occ_numbers(i_state, i_spin,i_k_point) = 0.5d0 * (1.d0 - arch_erf(hermite_arg))

              if (n_methfessel_paxton .gt. 0) then
             
                 !     first order contribution
                 A = - 0.25d0 * pisqrt_inv
                 !     H_even = H_0 = 1
                 H_even = 1.d0
                 !     H_odd =  H_1 = 2 * x
                 H_odd  = 2d0 * hermite_arg
                 occ_numbers(i_state, i_spin,i_k_point) = occ_numbers(i_state, i_spin, i_k_point) + A * H_odd * gauss_weight
              end if

              if (n_methfessel_paxton .gt. 1) then
                 do i_mp = 2, n_methfessel_paxton, 1

                    A = - 1.d0 / dble(4d0 * i_mp) * A
                    H_even = 2d0 * hermite_arg * H_odd  - 2d0 *  dble(i_mp)      * H_even
                    H_odd  = 2d0 * hermite_arg * H_even - 2d0 * (dble(i_mp) + 1d0) * H_odd
                    occ_numbers(i_state, i_spin, i_k_point) = occ_numbers(i_state, i_spin,i_k_point) + A * H_odd * gauss_weight
                 end do
              end if


           occ_numbers(i_state, i_spin,i_k_point) = occ_numbers(i_state, i_spin,i_k_point) &
                * spin_degeneracy 

           temp_n_electrons = temp_n_electrons + occ_numbers(i_state, i_spin,i_k_point) * k_weights(i_k_point)
        end do
     end do
  end do

  case (3)
    ! integer occupation (occ_numbers = 0, or spin_degenracy) test by igor
    if (n_k_points.gt.1) then
         write(use_unit,*) "* At present, the integer occupation scheme can not be used to the periodic model."
         write(use_unit,*) "* Abort."
         stop
    end if
    chemical_potential = 0.5d0*(KS_eigenvalue(1,1,1)+KS_eigenvalue(2,1,1))
    do i_k_point = 1, n_k_points,1
        temp_n_electrons = 0.d0
        i_state = 1
        diff_electrons = n_electrons
        do while (diff_electrons .gt. 1.0d-8)
           occ_numbers(i_state, 1, i_k_point) = spin_degeneracy
           temp_n_electrons = temp_n_electrons + occ_numbers(i_state, 1, i_k_point)
           diff_electrons = n_electrons - temp_n_electrons
           i_state = i_state + 1
        end do
        if (KS_eigenvalue(i_state-1,1,i_k_point).gt.chemical_potential) then
            chemical_potential = 0.5d0*(KS_eigenvalue(i_state-1,1,i_k_point)+KS_eigenvalue(i_state,1,i_k_point))
        end if
    end do

  case default
     write(use_unit,*) "* Unknown smearing type in subroutine check_norm."
     write(use_unit,*) "* Abort."
     stop
     
  end select


  diff_electrons = temp_n_electrons - n_electrons
  
end subroutine check_norm_sym

!   VB 2005:
!   Routine taken from http://www.netlib.org/go/
!
!   Copyright: My impression is that these routines come with a liberal (if nonexistant)
!   copyright provision, but netlib states that if in doubt, ask the original
!   author. In any case, it is publicly available for anyone to use.
!
!   Adapted by Ralf Gehrke 2005 to work for our purposes
!
!   Description from netlib
!   lib      ../go
!   for      Golden Oldies:  widely used,  but not in standard libraries.
!   #      Nominations welcome!
!   rel      excellent
!   age      old
!   editor      Eric Grosse
!   master      netlib.bell-labs.com
!
!  To get d1mach, mail netlib
!       send d1mach from core
!  VB: In my opinion, we already have the machine accuracy in question, but we can always get it
!      from netlib by downloading the dependencies also.
!
!  BB: Had to copy routine because is calls check_norm which imports n_k_points
!      -->calls check_norm_sym(...,n_k_points,...)
subroutine zeroin_sym(ax, bx, fax, fbx, tol, &
           occupation_acc, max_zeroin, &
           n_states, KS_eigenvalue, &
           n_electrons, chemical_potential, diff_electrons, &
           occ_numbers, i_counter, k_weights, n_k_points)



      use localorb_io, only: use_unit
      implicit none


! imported variables

! input
      real*8, intent(in) :: ax
      real*8, intent(in) :: bx
      real*8, intent(in) :: fax
      real*8, intent(in) :: fbx
      real*8, intent(in) :: tol

      real*8, intent(in) :: occupation_acc
      integer, intent(in) :: max_zeroin

      integer, intent(in) :: n_states
!      integer, intent(in) :: n_spin
      real*8, dimension(n_states), intent(in) :: &
           KS_eigenvalue
      real*8, intent(in) :: n_electrons

! output
      real*8, dimension(n_states), intent(out) :: &
           occ_numbers
!      real*8, dimension(n_states,n_spin,n_k_points), intent(out) ::
!     +     occ_numbers
!      real*8, dimension(:,:,:), intent(out) ::
!     +     occ_numbers
      real*8, intent(out) :: chemical_potential
      real*8, intent(out) :: diff_electrons
      integer, intent(inout) :: i_counter
      integer, intent(in) :: n_k_points
      real*8, dimension(n_k_points), intent(in) :: k_weights

!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
      real*8  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, &
           q, r, s
      real*8  dabs, d1mach

      eps = d1mach(4)
      tol1 = eps+1.0d0
!
      a  = ax
      b  = bx
      fa = fax
      fb = fbx


!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(use_unit,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' zeroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if (abs(fb) .le. occupation_acc) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge. &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
!  140 fb=f(b)
 140  continue
      call check_norm_sym(b, KS_eigenvalue, n_electrons, &
           occ_numbers, fb, i_counter, k_weights, n_k_points)
!test
!      write(use_unit,*) "i_counter = ", i_counter
!      write(use_unit,*) "n_electrons: ", n_electrons
!      write(use_unit,*) "occ_numbers: ", occ_numbers
!test end

      if (i_counter .gt. max_zeroin) then
!       force end despite non-convergence
        go to 150
      end if
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 chemical_potential = b
      diff_electrons     = fb
      return

end subroutine zeroin_sym
!****s* FHI-aims/sym_base/construct_C
!  NAME
!    construct_C
!  SYNOPSIS
subroutine construct_C(l_max, l, C)

    !  PURPOSE
    !  
    !  Construct the matrix for the transformation of the rotation matrix in 
    !  the basis of the complex!! spherical harmonics into the basis of the
    !  real sherical harmonics for on l.
    !  The matrix is constructed according to:
    !  Blanco, Florez, Bermejo, Journal of Molecular strucutre 
    !  (TheoChem) 419 (1997) 19-27
    !  The code was adapted from and test against aimsutil by C. Schober:
    !  https://gitlab.lrz.de/theochem/aimsutils
    !
    !  USES

    implicit none

    !  ARGUMENTS
    integer, intent(IN) :: l_max
    integer, intent(IN) :: l
    complex*16,dimension(-l_max:l_max,-l_max:l_max), intent(OUT) :: C 
    !  INPUTS
    !    o l_max -- maximum angular moment
    !    o l -- angular moment for which C is calculated
    !  OUTPUT
    !    o C -- Transfromation matrix
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
    !    Release version, FHI-aims (2016).
    !  SOURCE
    ! local variables

    complex*16 :: element

    integer :: m, m2
    
    C= dcmplx(0.d0,0.d0)

    do m = -l, l, 1
      do m2 =-l, l, 1
            if (abs(m).ne. abs(m2))then
                element = dcmplx(0.d0,0.d0)
            elseif ((m == 0) .and. (m2 == 0))then
                element = dcmplx(sqrt(2.d0),0.d0)

            elseif (abs(m) .gt. 0) then
                if ((abs(m).ne.m).and. (abs(m2) .ne. m2))then
                    ! upper left, m<0 + m2<0
                    element = dcmplx(0.d0,1.d0)
                elseif ((abs(m) == m) .and. (abs(m2) .ne. m2))then
                    ! lower left, m>0 + m2<0
                    element = dcmplx(1.d0,0.d0)
                elseif ((abs(m) == m) .and. (abs(m2) == m2))then
                    ! lower right, m>0 + m2>0
                    element = dcmplx((-1)**(l-(l-m)),0.d0)
                elseif (.not.(abs(m).eq.m) .and. (abs(m2).eq.m2))then
                    ! upper right, m<0 + m2 > 0
                    element = -dcmplx(0.d0,(-1.d0)**(l-(l-m)))
                endif
            endif
            C(m, m2) = element
       enddo
    enddo
            
    C = dcmplx((1.d0/sqrt(2.d0)),0.d0)*C


end subroutine construct_C
!****s* FHI-aims/sym_base/construct_C_full
!  NAME
!    construct_C_full
!  SYNOPSIS
subroutine construct_C_full(l_max, C_full)
    !  PURPOSE
    !  
    !  Construct the matrix for the transformation of the rotation matrix in 
    !  the basis of the complex!! spherical harmonics into the basis of the
    !  real sherical harmonics for all l uptp l_max.
    !
    !  USES

    implicit none

    !  ARGUMENTS
    integer, intent(IN) :: l_max
    complex*16,dimension(0:l_max,-l_max:l_max,-l_max:l_max) :: C_full
    !  INPUTS
    !    o l_max -- maximum angular moment
    !  OUTPUTS
    !    o C_full -- Transfromation matrix for all l
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    integer :: l
    !debug    integer :: m1, m2

    do l=0, l_max, 1
      call construct_C(l_max, l, C_full(l,:,:))
    enddo

    !debug    write(use_unit,'(A)') '---------------------------------------------------------------------'
    !debug    write(use_unit,'(A)') 'C Matrix'
    !debug    do l =0, l_max, 1
    !debug      do m1 = -l, l, 1
    !debug          do m2 = -l, l, 1
    !debug              write(use_unit,'(A,I4,A,I4,A,I4,A,2F11.4)') 'l',l, ' m1',m1, ' m2', m2, ' C:',&
    !debug                    real(C_full(l,m1,m2)), dimag(C_full(l,m1,m2))
    !debug          enddo
    !debug      enddo
    !debug    enddo
    !debug    write(use_unit,'(A)') '---------------------------------------------------------------------'

end subroutine construct_C_full
!****s* FHI-aims/sym_base/construct_Delta
!  NAME
!    construct_Delta
!  SYNOPSIS
subroutine construct_Delta(l_max,TVlmm,C_full,Delta)
    !  PURPOSE
    !  
    !  Transforms the Wigner D-Matrix to the basis of real spherical harmonics.
    !  Delta(l) = C_full^* x TVlmm x C_full^T
    !
    !  USES
    implicit none

    !  ARGUMENTS
    integer, intent(IN) :: l_max
    complex*16,dimension(0:l_max,-l_max:l_max,-l_max:l_max), intent(IN) :: C_full
    complex*16,dimension(0:l_max,-l_max:l_max,-l_max:l_max), intent(IN) :: TVlmm
    complex*16,dimension(0:l_max,-l_max:l_max,-l_max:l_max), intent(OUT) :: Delta
    !  INPUTS
    !    o l_max -- maximum angular moment
    !    o C_full -- Transfromation matrix for all l
    !    o TVlmm  -- Wigner D Matrix
    !  OUTPUTS
    !    o Delta  -- Rotation matrix in basis of real spherical harmonics   
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    integer :: l , m1, m2

    do l=0, l_max, 1
      Delta(l,:,:)=matmul(dconjg(C_full(l,:,:)),matmul(TVlmm(l,:,:),transpose(C_full(l,:,:))))
    enddo

         ! debug write(use_unit,'(A)') '---------------------------------------------------------------------'
         ! debug write(use_unit,'(A)') 'D Matrix'
         ! debug do l =0, l_max, 1
         ! debug   do m1 = -l, l, 1
         ! debug     do m2 = -l, l, 1
         ! debug         write(use_unit,'(A,I4,A,I4,A,I4,A,2F20.14)') 'l',l, ' m1',m1, ' m2', m2, ' Delta:', &
         ! debug                  dble(Delta(l,m1,m2)), dimag(Delta(l,m1,m2))
         ! debug     enddo
         ! debug   enddo
         ! debug enddo
         ! debug write(use_unit,'(A)') '---------------------------------------------------------------------'

end subroutine construct_Delta
!****s* FHI-aims/sym_base/apply_T
!  NAME
!    apply_T
!  SYNOPSIS
subroutine apply_T(l_max,C_full)
    !  PURPOSE
    !  
    !  Apply the FHI-aims sign convention for the real spherical harmonics to
    !  transformation matrix between Wigner-D and Delta matrix
    !
    !  USES
    use mpi_tasks,only: aims_stop
    implicit none

    !  ARGUMENTS
    integer, intent(IN) :: l_max
    complex*16,dimension(0:l_max,-l_max:l_max,-l_max:l_max) :: C_full
    !  INPUTS
    !    o l_max -- maximum angular moment
    !    o C_full -- Transfromation matrix for all l
    !  OUTPUTS
     !    o C_full -- Transfromation matrix for all l   
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    integer,dimension(0:12,-12:12) :: T
    integer,dimension(-l_max:l_max,-l_max:l_max) :: T_mult

    integer :: l, m1 , m2

    if (l_max.gt.12) then
      call aims_stop('Y_lm rotation matrix only implemented for l<=12.')
    endif
    ! hard coded
    T = 0
    T(0,:) =  (/0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0/) 
    T(1,:) =  (/0,0,0,0,0,0,0,0,0,0,0,1,1,-1,0,0,0,0,0,0,0,0,0,0,0/)
    T(2,:) =  (/0,0,0,0,0,0,0,0,0,0,1,1,1,-1,1,0,0,0,0,0,0,0,0,0,0/)
    T(3,:) =  (/0,0,0,0,0,0,0,0,0,1,1,1,1,-1,1,-1,0,0,0,0,0,0,0,0,0/)
    T(4,:) =  (/0,0,0,0,0,0,0,0,1,1,1,1,1,-1,1,-1,1,0,0,0,0,0,0,0,0/)
    T(5,:) =  (/0,0,0,0,0,0,0,1,1,1,1,1,1,-1,1,-1,1,-1,0,0,0,0,0,0,0/)
    T(6,:) =  (/0,0,0,0,0,0,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,0,0,0,0,0,0/)
    T(7,:) =  (/0,0,0,0,0,1,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,0,0,0,0,0/)
    T(8,:) =  (/0,0,0,0,1,1,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,0,0,0,0/)
    T(9,:) =  (/0,0,0,1,1,1,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,0,0,0/)
    T(10,:) = (/0,0,1,1,1,1,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,0,0/)
    T(11,:) = (/0,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,0/)
    T(12,:) = (/1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1/)


    do l = 0, l_max, 1
      T_mult = 0
      do m1 = -l ,l, 1
        T_mult(m1,m1)=T(l,m1)
      enddo
      C_full(l,:,:) = matmul(T_mult(:,:),C_full(l,:,:))
    enddo

        ! debug write(use_unit,'(A)') '---------------------------------------------------------------------'
        ! debug write(use_unit,'(A)') 'C Matrix with T'
        ! debug do l =0, l_max, 1
        ! debug   do m1 = -l, l, 1
        ! debug     do m2 = -l, l, 1
        ! debug         write(use_unit,'(A,I4,A,I4,A,I4,A,2F11.4)') 'l',l, ' m1',m1, ' m2', m2, ' C:',&
        ! debug                  real(C_full(l,m1,m2)), dimag(C_full(l,m1,m2))
        ! debug     enddo
        ! debug   enddo
        ! debug enddo
        ! debug write(use_unit,'(A)') '---------------------------------------------------------------------'

end subroutine apply_T

!****s* FHI-aims/sym_base/destroy_symmetry_arrays
!  NAME
!    destroy_symmetry_arrays
!  SYNOPSIS
subroutine destroy_symmetry_arrays()

  !  PURPOSE
  !
  !  Free Memory   
  !
  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
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

  ! Local Variables
  implicit none


  if (allocated(map_atom)) deallocate(map_atom)
  if (allocated(map_atom_full)) deallocate(map_atom_full)
  if (allocated(Delta)) deallocate(Delta)
  if (allocated(Delta_matrix_full)) deallocate ( Delta_matrix_full )
  if (allocated(rotations_at_task_map)) deallocate ( rotations_at_task_map )
  if (allocated(rotations_at_task)) deallocate ( rotations_at_task )
  if (allocated(k_point_list_nosym)) deallocate ( k_point_list_nosym )
  if (allocated(k_points_at_task)) deallocate ( k_points_at_task )
  if (allocated(k_phase_nosym)) deallocate ( k_phase_nosym )
  if (allocated(occ_numbers_nosym)) deallocate ( occ_numbers_nosym )
  if (allocated(k_weights_nosym)) deallocate ( k_weights_nosym )
  if (allocated(k_phase_base_nosym)) deallocate(k_phase_base_nosym)
  if (allocated( k_phase_exx_nosym)) deallocate(k_phase_exx_nosym)
  if (allocated(k_per_irr )) deallocate(k_per_irr)  
  if (allocated(n_k_per_irr )) deallocate(n_k_per_irr)
  if (allocated(k_per_task )) deallocate(k_per_task)
  if (allocated(n_k_per_task )) deallocate(n_k_per_task)
  if (allocated(k_phase_exx_nosym )) deallocate(k_phase_exx_nosym)
  if (allocated(map_back )) deallocate(map_back)
end subroutine destroy_symmetry_arrays

subroutine out_KSeigen()
    !  USES
    use dimensions,only: n_states,n_spin,n_k_points, n_k_points_nosym
    use runtime_choices,only: use_symmetry_reduced_spg
    use spglib_symmetry,only:map
    use mpi_tasks, only: check_allocation
    use physics, only: KS_eigenvalue
    !  ARGUMENTS
    !  INPUTS
    !    o l_max -- maximum angular moment
    !    o C_full -- Transfromation matrix for all l
    !  OUTPUTS
     !    o C_full -- Transfromation matrix for all l   
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    real*8,dimension(:,:,:), allocatable :: KS_eigenvalue_nosym
    integer:: i_k_point,i_spin,i_state
    integer:: info
    if ( use_symmetry_reduced_spg ) then 
      if (.not.allocated( KS_eigenvalue_nosym))then
        allocate (  KS_eigenvalue_nosym(n_states,n_spin,n_k_points_nosym) ,stat=info)
        call check_allocation(info, ' KS_eigenvalue_nosym      ')
      endif
      KS_eigenvalue_nosym = 0.d0
      do i_k_point = 1, n_k_points_nosym, 1
        KS_eigenvalue_nosym(1:n_states,1:n_spin,i_k_point)=KS_eigenvalue(1:n_states,1:n_spin,map(i_k_point))
      enddo
      open(67,file='eigen_irr.dat')
      i_spin=1
      do i_k_point = 1, n_k_points_nosym, 1
        do i_state=1,n_states
          write(67,'(I4,I4,1F20.14)') i_k_point, i_state, KS_eigenvalue_nosym(i_state,i_spin,i_k_point)
        enddo
      enddo
      close(67)
    else
      open(66,file='eigen_org.dat')
      i_spin=1
      do i_k_point = 1, n_k_points, 1
        do i_state=1,n_states
           write(66,'(I4,I4,1F20.14)') i_k_point, i_state, KS_eigenvalue(i_state,i_spin,i_k_point)
        enddo
      enddo
      close(66)    
    endif
endsubroutine

  !----------------------------------------------------------------------------
  !****s* sym_base/get_reducible_KS_dens
  !  NAME
  !    get_reducible_KS_dens
  !  SYNOPSIS
subroutine get_reducible_KS_dens(kdm, kdm_complex, density_matrix_sparse, &
                                 density_matrix,  force_packed , i_k, i_k_point)
    !  PURPOSE
    !
    !  Rotate unweighted irreducible density matrix to full set and sum over k   
    !
    !  USES
    use dimensions,only: n_basis, n_centers_basis_T, n_k_points_nosym
    use runtime_choices,only: real_eigenvectors,use_spg_full_Delta
    use spglib_symmetry,only:map_sym
    use mpi_tasks, only: check_allocation
    !  ARGUMENTS
    real*8, dimension(n_basis,n_basis) :: kdm
    complex*16, dimension(n_basis,n_basis):: kdm_complex   
    real*8, dimension(n_centers_basis_T,n_centers_basis_T), intent(INOUT) :: density_matrix
    ! when this routine is called, density_matrix_sparse has either the dimension
    ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
    ! so we declare it here as a 1D assumed size array
    real*8, intent(INOUT) :: density_matrix_sparse(*)
    logical, intent(IN) :: force_packed   
    integer, intent(IN) :: i_k 
    integer, intent(IN) :: i_k_point
    !  INPUTS
    !    o kdm, kdm_complex --  unweighted density at irreducible k-point
    !    o i_kpoint -- k-point in irreducible set
    !    o i_k -- task local k-point in irreducible set 
    !    o force_packed
    !  OUTPUTS
     !    o density_matrix, density_matrix_sparse -- density summed over k
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  ARGUMENTS
    real*8, dimension(:,:), allocatable :: kdm_work_L
    complex*16, dimension(:,:), allocatable:: kdm_complex_work_L
    real*8, dimension(:,:), allocatable :: kdm_work_R
    complex*16, dimension(:,:), allocatable:: kdm_complex_work_R
    complex*16, dimension(:,:), allocatable:: Delta_temp
    integer :: is, full_k_point, k_task, info
     if (real_eigenvectors) then
        allocate(kdm_work_L(n_basis, n_basis), stat=info)    
        call check_allocation(info, 'evaluate_densmat:kdm_work_L')
        allocate(kdm_work_R(n_basis, n_basis), stat=info)
        call check_allocation(info, 'evaluate_densmat:kdm_work_R')
        ! dummy for the real case 
        allocate(kdm_complex_work_L(1, 1), stat=info)          
        call check_allocation(info, 'evaluate_densmat:kdm_complex_work_L')
        allocate(kdm_complex_work_R(1, 1), stat=info)
        call check_allocation(info, 'evaluate_densmat:kdm_complex_work_R')
    else
        allocate(kdm_complex_work_L(n_basis, n_basis), stat=info)
        call check_allocation(info, 'evaluate_densmat:kdm_complex_work_L')
        allocate(kdm_complex_work_R(n_basis, n_basis), stat=info)
        call check_allocation(info, 'evaluate_densmat:kdm_complex_work_R')          
        ! dummy for the complex case
        allocate(kdm_work_L(1, 1), stat=info)
        call check_allocation(info, 'evaluate_densmat:kdm_work_L')
        allocate(kdm_work_R(1, 1), stat=info)
        call check_allocation(info, 'evaluate_densmat:kdm_work_R')          
    end if
    !if(.not.use_spg_full_Delta)then
        allocate(Delta_temp(n_basis, n_basis), stat=info)
        call check_allocation(info, 'evaluate_densmat:Delta_temp')   
    !endif   
    do full_k_point = 2, n_k_per_task(i_k)
       k_task = k_per_task(i_k,full_k_point)                  
       is = rotations_at_task_map(map_sym(k_task))
       if(use_spg_full_Delta)then
           Delta_temp = Delta_matrix_full(:,:, is,i_k)
       else
           call get_density_rotation(rotations_at_task(is), i_k_point, &
                        map_atom(:,is), Delta(:,:,:,is), Delta_temp)
       endif
       if (real_eigenvectors) then
         call dsymm ('R', 'U',  n_basis, n_basis,k_weights_nosym(k_task), &
                 kdm, n_basis,dble(Delta_temp),n_basis,0.d0,&
                 kdm_work_L, n_basis)  
         call dgemm ('N', 'T', n_basis, n_basis, n_basis, 1.d0, &
                      kdm_work_L, n_basis,dble(Delta_temp),n_basis,0.d0,&
                      kdm_work_R, n_basis)  
       else                    
          call zhemm ('R', 'U',  n_basis, n_basis, dcmplx(k_weights_nosym(k_task),0.d0), &
                kdm_complex, n_basis,Delta_temp,n_basis,(0.d0,0.d0),&
                kdm_complex_work_L, n_basis)  
          call zgemm ('N', 'C', n_basis, n_basis, n_basis, (1.d0,0.d0), &
                kdm_complex_work_L, n_basis,Delta_temp,n_basis,(0.d0,0.d0),&
                kdm_complex_work_R, n_basis)  
        endif
        call accumulate_k_densmat_sym(density_matrix_sparse, density_matrix, force_packed, &
                & kdm_work_R, kdm_complex_work_R, k_phase_nosym(:,k_task))
    enddo
    if (allocated(kdm_work_L)) deallocate(kdm_work_L)
    if (allocated(kdm_complex_work_L)) deallocate(kdm_complex_work_L)   
    if (allocated(kdm_work_R)) deallocate(kdm_work_R)
    if (allocated(kdm_complex_work_R)) deallocate(kdm_complex_work_R) 
    if (allocated(Delta_temp)) deallocate(Delta_temp)   
    
    
endsubroutine get_reducible_KS_dens

  !----------------------------------------------------------------------------
  !****s* sym_base/evaluate_densmat_hf_sym
  !  NAME
  !    evaluate_densmat_hf_sym
  !  SYNOPSIS

  subroutine evaluate_densmat_hf_sym(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

    !  PURPOSE
    !    Evaluates the density matrix at the k-points and stores it in
    !    hf_exchange_matr_real/complex as an intermediate location
    !    
    !    The full density matrix is reconstructed from the irreducible
    !    k-points. hf_exchange_matr_real/complex is deallocated and allocated
    !    for full number of k-points per task. 
    !  USES
    use dimensions,only: n_basis, n_k_points, n_states, n_spin, n_k_points_task,&
                         n_k_points_nosym
    use runtime_choices,only: real_eigenvectors,use_spg_full_Delta
    use spglib_symmetry,only:map_sym
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use hartree_fock_p0, only: hf_exchange_matr_real, hf_exchange_matr_complex
    implicit none

    !  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector_complex
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers

    !  INPUTS
    !   o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !   o occ_numbers -- occupation of states
    !
    !  OUTPUT
    !   only via hf_exchange_matr_real/complex
    !
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
    !

    ! local variables

    real*8, allocatable :: tmp_r(:,:)
    complex*16, allocatable :: tmp_c(:,:)
    complex*16, dimension(:,:), allocatable:: Delta_temp
    real*8, allocatable :: KS_scaled(:,:)
    complex*16, allocatable :: KS_scaled_cmplx(:,:)
    
    integer :: i_state, i_spin, i_k, i_k_point, info, i_k_sym, full_k_point, k_task, is, i_bas1, i_bas2

    if(real_eigenvectors)then
      if(allocated(hf_exchange_matr_real)) then
        deallocate(hf_exchange_matr_real)
        allocate(hf_exchange_matr_real(n_basis,n_basis,n_k_points_task_full,n_spin))
      endif
    else
      if(allocated(hf_exchange_matr_complex)) then
        deallocate(hf_exchange_matr_complex)
        allocate(hf_exchange_matr_complex(n_basis,n_basis,n_k_points_task_full,n_spin))
      endif
    endif
         
    if(real_eigenvectors)then
      allocate(tmp_r(n_basis,n_states),stat=info)
      allocate(KS_scaled(n_basis, n_states), stat=info)
    else
      allocate(tmp_c(n_basis,n_states),stat=info)
      allocate(KS_scaled_cmplx(n_basis, n_states), stat=info)
    endif
    allocate(Delta_temp(n_basis, n_basis), stat=info)
    call check_allocation(info, 'evaluate_densmat:Delta_temp') 
    
    do i_spin = 1, n_spin   
      i_k = 0
      i_k_sym = 0
      do i_k_point = 1, n_k_points
        if(myid == MOD(i_k_point, n_tasks)) then
            i_k = i_k + 1
            ! irreducible k-points per task
            if(real_eigenvectors)then
              do full_k_point = 1, n_k_per_task(i_k)
                ! All k-points from rotated from i_k
                i_k_sym = i_k_sym + 1
                k_task = k_per_task(i_k,full_k_point)                  
                if(map_sym(k_task).eq.1)then
                  ! Identity (irreducible), do standard
                  do i_state = 1, n_states
                    tmp_r(:,i_state) = KS_eigenvector(:,i_state,i_spin,i_k)*occ_numbers(i_state,i_spin,i_k_point)
                  enddo
                  call dgemm('N','T',n_basis,n_basis,n_states,1.d0,KS_eigenvector(1,1,i_spin,i_k),&
                             ubound(KS_eigenvector,1),tmp_r,ubound(tmp_r,1),0.d0,&
                             hf_exchange_matr_real(1,1,i_k,i_spin),ubound(hf_exchange_matr_real,1))
                else 
                  ! Rotations
                      is = rotations_at_task_map(map_sym(k_task))
                      if(use_spg_full_Delta)then
                          Delta_temp = Delta_matrix_full(:,:, is,i_k)      
                      else
                          call get_density_rotation(rotations_at_task(is), i_k_point, &
                                        map_atom(:,is), Delta(:,:,:,is), Delta_temp)
                      endif 
                      ! Rotate KS_eigenvector
                      do i_state = 1, n_states
                        call dgemv('N', n_basis,  n_basis, 1.d0, dble(Delta_temp), n_basis, &
                        &          KS_eigenvector(1, i_state, i_spin, i_k), 1, 0.d0,&
                                KS_scaled(1, i_state), 1 )
                        tmp_r(:,i_state) = KS_scaled(:,i_state)*occ_numbers(i_state,i_spin,i_k_point)
                      enddo
                      call dgemm('N','T',n_basis,n_basis,n_states,(1.d0,0.d0),KS_scaled, &
                              ubound(KS_scaled,1), &
                              tmp_r,ubound(tmp_r,1),(0.d0,0.d0),&
                              hf_exchange_matr_real(1,1,i_k_sym,i_spin),ubound(hf_exchange_matr_real,1))
                endif
              enddo
            else              
              do full_k_point = 1, n_k_per_task(i_k)
                i_k_sym = i_k_sym + 1
                k_task = k_per_task(i_k,full_k_point)    
                if(map_sym(k_task).eq.1)then
                    ! Identity (irreducible), do standard
                    do i_state = 1, n_states
                      tmp_c(:,i_state) = KS_eigenvector_complex(:,i_state, i_spin, i_k)*occ_numbers(i_state,i_spin,i_k_point)
                    enddo
                    call zgemm('N','C',n_basis,n_basis,n_states,(1.d0,0.d0),KS_eigenvector_complex(1,1,i_spin,i_k), &
                              ubound(KS_eigenvector_complex,1), &
                              tmp_c,ubound(tmp_c,1),(0.d0,0.d0),&
                              hf_exchange_matr_complex(1,1,i_k_sym,i_spin),ubound(hf_exchange_matr_complex,1))
                else 
                  ! Rotations
                    is = rotations_at_task_map(map_sym(k_task))
                    if(use_spg_full_Delta)then
                        Delta_temp = Delta_matrix_full(:,:, is,i_k)      
                    else
                        call get_density_rotation(rotations_at_task(is), i_k_point, &
                                      map_atom(:,is), Delta(:,:,:,is), Delta_temp)
                    endif
                    ! Rotate KS_eigenvector
                    do i_state = 1, n_states
                      call zgemv('N', n_basis,  n_basis, dcmplx(1.d0,0.d0), Delta_temp, n_basis, &
                          &       KS_eigenvector_complex(1,i_state, i_spin, i_k), 1, dcmplx(0.d0,0.d0),&
                                  KS_scaled_cmplx(1,i_state), 1 )
                      tmp_c(:,i_state) = KS_scaled_cmplx(:,i_state)*occ_numbers(i_state,i_spin,i_k_point)
                    enddo
                    call zgemm('N','C',n_basis,n_basis,n_states,(1.d0,0.d0),KS_scaled_cmplx, &
                            ubound(KS_scaled_cmplx,1), &
                            tmp_c,ubound(tmp_c,1),(0.d0,0.d0),&
                            hf_exchange_matr_complex(1:n_basis,1:n_basis,i_k_sym,i_spin),ubound(hf_exchange_matr_complex,1))
                endif  ! sym = identity    
              enddo ! n_k_per_task
            endif ! real_eigenvectors
        endif ! mod n_tasks
      enddo ! i_k_point
    enddo ! i_spin
    if(real_eigenvectors)then
      deallocate(tmp_r)
    else
      deallocate(tmp_c)
    endif
    deallocate(Delta_temp)
end subroutine evaluate_densmat_hf_sym
  !----------------------------------------------------------------------------
  !****s* sym_base/FT_densmat_sym
  !  NAME
  !    FT_densmat_sym
  !  SYNOPSIS

subroutine FT_densmat_sym(my_n_atoms,my_n_basis,my_atom_list,my_basis_off,dm_tmp,dm_cols)

    !  PURPOSE
    !    FT and add density matrix to dm_cols (for eack BvK cell and k-point)
    !    If not calculate previously the the full density is rotated from 
    !    the irreducible set of k-points for each cell.
    !  USES
    use dimensions,only: n_basis, n_k_points, n_states, n_spin, n_k_points_task,&
			n_atoms, n_k_points_nosym
    use runtime_choices,only: real_eigenvectors,use_spg_full_Delta, &
                              get_full_density_first, use_k_phase_exx, n_k_points_xyz_nosym
    use spglib_symmetry,only:map_sym
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use hartree_fock_p0, only: hf_exchange_matr_real, hf_exchange_matr_complex
    use basis, only: atom2basis_off, sp2n_basis_sp
    use geometry, only: species
    use pbc_lists, only: n_cells_bvk
    use synchronize_mpi_basic,only:sync_vector    
    implicit none

    !  ARGUMENTS
    real*8, dimension(n_basis,n_basis):: dm_tmp
    integer:: my_n_atoms
    integer:: my_n_basis
    integer, dimension(n_atoms):: my_atom_list
    integer, dimension(n_atoms):: my_basis_off
    real*8, dimension(n_basis,my_n_basis,n_cells_bvk,n_spin):: dm_cols
    !  INPUTS
    !   o my_n_atoms - node local atoms
    !   o my_n_basis - node local basis
    !   o my_atom_list - node local atom list
    !   o my_basis_off - node local basis
    !   o dm_tmp - workspace for density
    !
    !  OUTPUT
    !   o  dm_cols
    !
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
    !

    ! local variables
    real*8, dimension(:,:), allocatable :: kdm_work_L
    complex*16, dimension(:,:), allocatable:: kdm_complex_work_L
    real*8, dimension(:,:), allocatable :: kdm_work_R
    complex*16, dimension(:,:), allocatable:: kdm_complex_work_R
    complex*16, dimension(:,:), allocatable:: Delta_temp
    integer :: info, is, full_k_point, k_task, i_cell_1, i_cell_2, i_cell_3,&
              i_k_sym, i_bas1, i_bas2, ind, i_spin, i_cell, i_k, i_k_point, &
              i_my_atom, i_atom, i_basis_s, i_basis_e, my_s, my_e
    complex*16 :: k_phase_new
    ! Statement functions for offset of basis/basbas (only for better readability!!!)
    integer :: atom2basis_len
    atom2basis_len(i_atom)  = sp2n_basis_sp(species(i_atom))
    if(use_k_phase_exx)then
        if(.not.get_full_density_first)then
          if (real_eigenvectors) then
                allocate(kdm_work_L(n_basis, n_basis), stat=info)    
                call check_allocation(info, 'evaluate_densmat:kdm_work_L')
                allocate(kdm_work_R(n_basis, n_basis), stat=info)
                call check_allocation(info, 'evaluate_densmat:kdm_work_R')
          else
                allocate(kdm_complex_work_L(n_basis, n_basis), stat=info)
                call check_allocation(info, 'evaluate_densmat:kdm_complex_work_L')
                allocate(kdm_complex_work_R(n_basis, n_basis), stat=info)
                call check_allocation(info, 'evaluate_densmat:kdm_complex_work_R')               
          end if
          !if(.not.use_spg_full_Delta)then
              allocate(Delta_temp(n_basis, n_basis), stat=info)
              call check_allocation(info, 'evaluate_densmat:Delta_temp')   
          !endif   
          do i_cell = 1, n_cells_bvk          
                    do i_spin = 1, n_spin
                      dm_tmp(:,:) = 0.
                      i_k = 0
                      do i_k_point = 1, n_k_points,1
                        if(myid == MOD(i_k_point, n_tasks)) then
                          i_k = i_k + 1
                          k_task = k_per_task(i_k,1)
                          if(real_eigenvectors)then
                            dm_tmp(:,:) = dm_tmp(:,:) + &
                                hf_exchange_matr_real   (:,:,i_k,i_spin)*&
                                dble (k_phase_exx_nosym(i_cell,k_task)*k_weights_nosym(k_task))
                          else
                            dm_tmp(:,:) = dm_tmp(:,:) + &
                                dble(hf_exchange_matr_complex(:,:,i_k,i_spin)*&
                                conjg(k_phase_exx_nosym(i_cell,k_task)*k_weights_nosym(k_task)))
                          endif
                          do full_k_point = 2, n_k_per_task(i_k)
                                  k_task = k_per_task(i_k,full_k_point)                  
                                  is = rotations_at_task_map(map_sym(k_task))
                                  if(use_spg_full_Delta)then
                                      Delta_temp = Delta_matrix_full(:,:, is,i_k)
                                  else
                                      call get_density_rotation(rotations_at_task(is), i_k_point, &
                                                    map_atom(:,is), Delta(:,:,:,is), Delta_temp)
                                  endif
                                  if (real_eigenvectors) then
                                    call dsymm ('R', 'U',  n_basis, n_basis,1.d0, &
                                            hf_exchange_matr_real   (:,:,i_k,i_spin), &
                                            n_basis,dble(Delta_temp),n_basis,0.d0,&
                                            kdm_work_L, n_basis)  
                                    call dgemm ('N', 'T', n_basis, n_basis, n_basis, 1.d0, &
                                                  kdm_work_L, n_basis,dble(Delta_temp),n_basis,0.d0,&
                                                  kdm_work_R, n_basis)  
                                    dm_tmp(:,:) = dm_tmp(:,:) + &
                                        kdm_work_R*dble (k_phase_exx_nosym(i_cell,k_task)*k_weights_nosym(k_task))
                                  else                    
                                      call zhemm ('R', 'U',  n_basis, n_basis,(1.d0,0.d0), &
                                            hf_exchange_matr_complex(:,:,i_k,i_spin), n_basis,&
                                            Delta_temp,n_basis,(0.d0,0.d0),&
                                            kdm_complex_work_L, n_basis)  
                                      call zgemm ('N', 'C', n_basis, n_basis, n_basis, (1.d0,0.d0), &
                                            kdm_complex_work_L, n_basis,Delta_temp,n_basis,(0.d0,0.d0),&
                                            kdm_complex_work_R, n_basis)  
                                      dm_tmp(:,:) = dm_tmp(:,:) + &
                                        dble(kdm_complex_work_R*conjg(k_phase_exx_nosym(i_cell,k_task)*k_weights_nosym(k_task)))
                                  endif
                                enddo
                        endif
                      enddo
                      
                      call sync_vector(dm_tmp, size(dm_tmp))

                      do i_my_atom = 1, my_n_atoms
                        i_atom = my_atom_list(i_my_atom)
                        i_basis_s = atom2basis_off(i_atom) + 1
                        i_basis_e = atom2basis_off(i_atom) + atom2basis_len(i_atom)
                        my_s = my_basis_off(i_atom) + 1
                        my_e = my_basis_off(i_atom) + atom2basis_len(i_atom)
                        !SVL this is original version:
                        !dm_cols(:,my_s:my_e,i_cell,i_spin) = dm_tmp(:,i_basis_s:i_basis_s)
                        !New version:
                        dm_cols(:,my_s:my_e,i_cell,i_spin) = dm_tmp(:,i_basis_s:i_basis_e)
                      enddo

                    enddo  
          enddo
          if (allocated(kdm_work_L)) deallocate(kdm_work_L)
          if (allocated(kdm_complex_work_L)) deallocate(kdm_complex_work_L)   
          if (allocated(kdm_work_R)) deallocate(kdm_work_R)
          if (allocated(kdm_complex_work_R)) deallocate(kdm_complex_work_R) 
          if (allocated(Delta_temp)) deallocate(Delta_temp) 
        else
          do i_cell = 1, n_cells_bvk          
                    do i_spin = 1, n_spin
                      dm_tmp(:,:) = 0.
                      i_k = 0
                      i_k_sym = 0
                      do i_k_point = 1, n_k_points,1
                        if(myid == MOD(i_k_point, n_tasks)) then
                          i_k = i_k + 1
                          do full_k_point = 1, n_k_per_task(i_k)
                            i_k_sym = i_k_sym + 1
                            k_task = k_per_task(i_k,full_k_point)
                            if (real_eigenvectors) then
                              dm_tmp(:,:) = dm_tmp(:,:) + &
                                  hf_exchange_matr_real(:,:,i_k_sym,i_spin)*&
                                  dble (k_phase_exx_nosym(i_cell,k_task)*k_weights_nosym(k_task))
                            else                    
                                dm_tmp(:,:) = dm_tmp(:,:) + &
                                  dble(hf_exchange_matr_complex(:,:,i_k_sym,i_spin)*&
                                  conjg(k_phase_exx_nosym(i_cell,k_task)*k_weights_nosym(k_task)))
                            endif
                          enddo
                        endif
                      enddo
                      
                      call sync_vector(dm_tmp, size(dm_tmp))

                      do i_my_atom = 1, my_n_atoms
                        i_atom = my_atom_list(i_my_atom)
                        i_basis_s = atom2basis_off(i_atom) + 1
                        i_basis_e = atom2basis_off(i_atom) + atom2basis_len(i_atom)
                        my_s = my_basis_off(i_atom) + 1
                        my_e = my_basis_off(i_atom) + atom2basis_len(i_atom)
                        dm_cols(:,my_s:my_e,i_cell,i_spin) = dm_tmp(:,i_basis_s:i_basis_e)
                      enddo

                    enddo 
          enddo 
          if(real_eigenvectors)then
            if(allocated(hf_exchange_matr_real)) then
              deallocate(hf_exchange_matr_real)
              allocate(hf_exchange_matr_real(n_basis,n_basis,n_k_points_task,n_spin))
            endif
          else
            if(allocated(hf_exchange_matr_complex)) then
              deallocate(hf_exchange_matr_complex)
              allocate(hf_exchange_matr_complex(n_basis,n_basis,n_k_points_task,n_spin))
            endif
          endif        
        endif
    else
        if(.not.get_full_density_first)then
          if (real_eigenvectors) then
                allocate(kdm_work_L(n_basis, n_basis), stat=info)    
                call check_allocation(info, 'evaluate_densmat:kdm_work_L')
                allocate(kdm_work_R(n_basis, n_basis), stat=info)
                call check_allocation(info, 'evaluate_densmat:kdm_work_R')
          else
                allocate(kdm_complex_work_L(n_basis, n_basis), stat=info)
                call check_allocation(info, 'evaluate_densmat:kdm_complex_work_L')
                allocate(kdm_complex_work_R(n_basis, n_basis), stat=info)
                call check_allocation(info, 'evaluate_densmat:kdm_complex_work_R')               
          end if
          !if(.not.use_spg_full_Delta)then
              allocate(Delta_temp(n_basis, n_basis), stat=info)
              call check_allocation(info, 'evaluate_densmat:Delta_temp')   
          !endif   
          i_cell = 1          
          do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
            do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2              
                  if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then
                    i_cell = i_cell + 1
                    ind = i_cell
                  else
                    ind = 1
                  endif
                    do i_spin = 1, n_spin
                      dm_tmp(:,:) = 0.
                      i_k = 0
                      do i_k_point = 1, n_k_points,1
                        if(myid == MOD(i_k_point, n_tasks)) then
                          i_k = i_k + 1
                          k_task = k_per_task(i_k,1)
                          k_phase_new = k_phase_base_nosym(1, k_task)**i_cell_1 &
                            * k_phase_base_nosym(2, k_task)**i_cell_2 &
                            * k_phase_base_nosym(3, k_task)**i_cell_3
                          if(real_eigenvectors)then
                            dm_tmp(:,:) = dm_tmp(:,:) + &
                                hf_exchange_matr_real   (:,:,i_k,i_spin)*&
                                dble (k_phase_new*k_weights_nosym(k_task))
                          else
                            dm_tmp(:,:) = dm_tmp(:,:) + &
                                dble(hf_exchange_matr_complex(:,:,i_k,i_spin)*&
                                conjg(k_phase_new*k_weights_nosym(k_task)))
                          endif
                          do full_k_point = 2, n_k_per_task(i_k)
                                  k_task = k_per_task(i_k,full_k_point)                  
                                  is = rotations_at_task_map(map_sym(k_task))
                                  k_phase_new = k_phase_base_nosym(1, k_task)**i_cell_1 &
                                          * k_phase_base_nosym(2, k_task)**i_cell_2 &
                                          * k_phase_base_nosym(3, k_task)**i_cell_3
                                  if(use_spg_full_Delta)then
                                      Delta_temp = Delta_matrix_full(:,:, is,i_k)
                                  else
                                      call get_density_rotation(rotations_at_task(is), i_k_point, &
                                                    map_atom(:,is), Delta(:,:,:,is), Delta_temp)
                                  endif
                                  if (real_eigenvectors) then
                                    call dsymm ('R', 'U',  n_basis, n_basis,1.d0, &
                                            hf_exchange_matr_real   (:,:,i_k,i_spin), &
                                            n_basis,dble(Delta_temp),n_basis,0.d0,&
                                            kdm_work_L, n_basis)  
                                    call dgemm ('N', 'T', n_basis, n_basis, n_basis, 1.d0, &
                                                  kdm_work_L, n_basis,dble(Delta_temp),n_basis,0.d0,&
                                                  kdm_work_R, n_basis)  
                                    dm_tmp(:,:) = dm_tmp(:,:) + &
                                        kdm_work_R*dble (k_phase_new*k_weights_nosym(k_task))
                                  else                    
                                      call zhemm ('R', 'U',  n_basis, n_basis,(1.d0,0.d0), &
                                            hf_exchange_matr_complex(:,:,i_k,i_spin), n_basis,&
                                            Delta_temp,n_basis,(0.d0,0.d0),&
                                            kdm_complex_work_L, n_basis)  
                                      call zgemm ('N', 'C', n_basis, n_basis, n_basis, (1.d0,0.d0), &
                                            kdm_complex_work_L, n_basis,Delta_temp,n_basis,(0.d0,0.d0),&
                                            kdm_complex_work_R, n_basis)  
                                      dm_tmp(:,:) = dm_tmp(:,:) + &
                                        dble(kdm_complex_work_R*conjg(k_phase_new*k_weights_nosym(k_task)))
                                  endif
                                enddo
                        endif
                      enddo
                      
                      call sync_vector(dm_tmp, size(dm_tmp))

                      do i_my_atom = 1, my_n_atoms
                        i_atom = my_atom_list(i_my_atom)
                        i_basis_s = atom2basis_off(i_atom) + 1
                        i_basis_e = atom2basis_off(i_atom) + atom2basis_len(i_atom)
                        my_s = my_basis_off(i_atom) + 1
                        my_e = my_basis_off(i_atom) + atom2basis_len(i_atom)
                        !SVL this is original version:
                        !dm_cols(:,my_s:my_e,i_cell,i_spin) = dm_tmp(:,i_basis_s:i_basis_s)
                        !New version:
                        dm_cols(:,my_s:my_e,ind,i_spin) = dm_tmp(:,i_basis_s:i_basis_e)
                      enddo

                    enddo  
              enddo
            enddo
          enddo
          if (allocated(kdm_work_L)) deallocate(kdm_work_L)
          if (allocated(kdm_complex_work_L)) deallocate(kdm_complex_work_L)   
          if (allocated(kdm_work_R)) deallocate(kdm_work_R)
          if (allocated(kdm_complex_work_R)) deallocate(kdm_complex_work_R) 
          if (allocated(Delta_temp)) deallocate(Delta_temp) 
        else
          i_cell = 1          
          do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
            do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2              
                  if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then
                    i_cell = i_cell + 1
                    ind = i_cell
                  else
                    ind = 1
                  endif
                    do i_spin = 1, n_spin
                      dm_tmp(:,:) = 0.
                      i_k = 0
                      i_k_sym = 0
                      do i_k_point = 1, n_k_points,1
                        if(myid == MOD(i_k_point, n_tasks)) then
                          i_k = i_k + 1
                          do full_k_point = 1, n_k_per_task(i_k)
                            i_k_sym = i_k_sym + 1
                            k_task = k_per_task(i_k,full_k_point)
                            !write(use_unit,*) i_k_point, i_k_sym, k_task
                            k_phase_new = k_phase_base_nosym(1, k_task)**i_cell_1 &
                                    * k_phase_base_nosym(2, k_task)**i_cell_2 &
                                    * k_phase_base_nosym(3, k_task)**i_cell_3
                            if (real_eigenvectors) then
                              dm_tmp(:,:) = dm_tmp(:,:) + &
                                  hf_exchange_matr_real(:,:,i_k_sym,i_spin)*&
                                  dble (k_phase_new*k_weights_nosym(k_task))
                            else                    
                                dm_tmp(:,:) = dm_tmp(:,:) + &
                                  dble(hf_exchange_matr_complex(:,:,i_k_sym,i_spin)*&
                                  conjg(k_phase_new*k_weights_nosym(k_task)))
                            endif
                          enddo
                        endif
                      enddo
                      
                      call sync_vector(dm_tmp, size(dm_tmp))

                      do i_my_atom = 1, my_n_atoms
                        i_atom = my_atom_list(i_my_atom)
                        i_basis_s = atom2basis_off(i_atom) + 1
                        i_basis_e = atom2basis_off(i_atom) + atom2basis_len(i_atom)
                        my_s = my_basis_off(i_atom) + 1
                        my_e = my_basis_off(i_atom) + atom2basis_len(i_atom)
                        dm_cols(:,my_s:my_e,ind,i_spin) = dm_tmp(:,i_basis_s:i_basis_e)
                      enddo

                    enddo 
              enddo
            enddo
          enddo 
          if(real_eigenvectors)then
            if(allocated(hf_exchange_matr_real)) then
              deallocate(hf_exchange_matr_real)
              allocate(hf_exchange_matr_real(n_basis,n_basis,n_k_points_task,n_spin))
            endif
          else
            if(allocated(hf_exchange_matr_complex)) then
              deallocate(hf_exchange_matr_complex)
              allocate(hf_exchange_matr_complex(n_basis,n_basis,n_k_points_task,n_spin))
            endif
          endif        
        endif
    endif
  endsubroutine FT_densmat_sym
  !----------------------------------------------------------------------------
  !****s* sym_base/FT_hf_exchange_matr_sym
  !  NAME
  !    FT_hf_exchange_matr_sym
  !  SYNOPSIS
subroutine FT_hf_exchange_matr_sym(my_n_basis,my_basis_off,&
                                   fock_tmp,fock_matrix,fock_tmp_SR,fock_matrix_SR)

    !  PURPOSE
    !
    !  "Construct" routine for Fock matrix. Real space Fock matrix is Fourier
    !  transformed by summing over BvK cells with k_phase
    !  For symmetry we only use the irreducible k-point set and map those to
    !  the k_phase array for the full set
    !
    !  Adapted from calculate_fock_matrix_p0
    !  USES
    use dimensions,only: n_basis, n_k_points, n_states, n_spin ,n_atoms,&
                         use_lc_wpbeh, n_k_points_nosym
    use runtime_choices,only: real_eigenvectors, hybrid_coeff,use_scalapack, &
                              use_k_phase_exx, n_k_points_xyz_nosym
    use spglib_symmetry,only: map, map_sym
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use hartree_fock_p0, only: hf_exchange_matr_real, hf_exchange_matr_complex, &
     hf_exchange_matr_real_SR, hf_exchange_matr_complex_SR
    use basis, only: atom2basis_off, sp2n_basis_sp, max_n_basis_sp
    use geometry, only: species
    use pbc_lists, only: n_cells_bvk, k_phase_base
    use scalapack_wrapper, only: l_col, l_row
    use synchronize_mpi_basic,only:sync_vector 
    implicit none

    !  ARGUMENTS
    real*8, dimension(n_basis,max_n_basis_sp):: fock_tmp
    real*8, dimension(n_basis,max_n_basis_sp):: fock_tmp_SR
    integer:: my_n_basis
    integer, dimension(n_atoms):: my_basis_off
    real*8, dimension(n_basis,my_n_basis,n_cells_bvk,n_spin):: fock_matrix
    real*8, dimension(n_basis,my_n_basis,n_cells_bvk,n_spin):: fock_matrix_SR
    !  INPUTS
    !   o my_n_basis    - node local basis
    !   o my_basis_off  - node local basis offset
    !   o fock_tmp      - work space
    !   o fock_tmp_SR   - work space
    !   o fock_matrix   - Piece of Fock Matrix to add
    !   o fock_matrix_SR- Piece of Fock Matrix to add
    !
    !  OUTPUT
    !   o output via hf_exchange_matr{_complex}
    !
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
    !

    ! local variables
    real*8, dimension(:,:), allocatable :: kdm_work_L
    complex*16, dimension(:,:), allocatable:: kdm_complex_work_L
    real*8, dimension(:,:), allocatable :: kdm_work_R
    complex*16, dimension(:,:), allocatable:: kdm_complex_work_R
    complex*16, dimension(:,:), allocatable:: Delta_temp
    integer :: info, i_cell_1, i_cell_2, i_cell_3, ind, i_spin, i_k, i_k_point, &
              i_atom, my_s, my_e, i_cell_fock, my_k_point, i, i_basis, j, full_k_point
    complex*16 :: k_phase_new, cfact
    ! Statement functions for offset of basis/basbas (only for better readability!!!)
    integer :: atom2basis_len
    atom2basis_len(i_atom)  = sp2n_basis_sp(species(i_atom))

    
    if(use_k_phase_exx)then
      do i_atom = 1, n_atoms
        do i_cell_fock = 1, n_cells_bvk
          do i_spin = 1, n_spin
                    if(my_basis_off(i_atom) >= 0) then
                      my_s = my_basis_off(i_atom) + 1
                      my_e = my_basis_off(i_atom) + atom2basis_len(i_atom)
                      fock_tmp(:,1:atom2basis_len(i_atom)) = fock_matrix(:,my_s:my_e,i_cell_fock,i_spin)
                      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          fock_tmp_SR(:,1:atom2basis_len(i_atom)) = fock_matrix_SR(:,my_s:my_e,i_cell_fock,i_spin)
                      end if
                    else
                      fock_tmp(:,1:atom2basis_len(i_atom)) = 0
                      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          fock_tmp_SR(:,1:atom2basis_len(i_atom)) = 0
                      end if
                    endif

                    call sync_vector(fock_tmp, n_basis*atom2basis_len(i_atom))
                    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          call sync_vector(fock_tmp_SR, n_basis*atom2basis_len(i_atom))
                    end if

                    if(use_scalapack) then
                      full_k_point=map_back(i_k_point)
                      cfact = k_phase_exx_nosym(i_cell_fock,full_k_point)*0.5*n_spin

                      do i = 1, atom2basis_len(i_atom)
                        i_basis = atom2basis_off(i_atom) + i
                        if(l_col(i_basis) > 0) then
                          do j = 1, n_basis
                            if(l_row(j)>0) then
                              if(real_eigenvectors)then
                                hf_exchange_matr_real(l_row(j),l_col(i_basis),1,i_spin) = &
                                  hf_exchange_matr_real(l_row(j),l_col(i_basis),1,i_spin) + dble(fock_tmp(j,i)*cfact)
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_real_SR(l_row(j),l_col(i_basis),1,i_spin) = &
                                          hf_exchange_matr_real_SR(l_row(j),l_col(i_basis),1,i_spin) + dble(fock_tmp_SR(j,i)*cfact)
                                end if
                              else
                                hf_exchange_matr_complex(l_row(j),l_col(i_basis),1,i_spin) = &
                                  hf_exchange_matr_complex(l_row(j),l_col(i_basis),1,i_spin) + fock_tmp(j,i)*cfact
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_complex_SR(l_row(j),l_col(i_basis),1,i_spin) = &
                                          hf_exchange_matr_complex_SR(l_row(j),l_col(i_basis),1,i_spin) + fock_tmp_SR(j,i)*cfact
                                end if
                              endif
                            endif
                          enddo
                        endif

                        if(l_row(i_basis) > 0) then
                          do j = 1, n_basis
                            if(l_col(j)>0) then
                              if(real_eigenvectors)then
                                hf_exchange_matr_real(l_row(i_basis),l_col(j),1,i_spin) = &
                                  hf_exchange_matr_real(l_row(i_basis),l_col(j),1,i_spin) + dble(fock_tmp(j,i)*conjg(cfact))
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_real_SR(l_row(i_basis),l_col(j),1,i_spin) = &
                                          hf_exchange_matr_real_SR(l_row(i_basis),l_col(j),1,i_spin) + dble(fock_tmp_SR(j,i)*conjg(cfact))
                                end if
                              else
                                hf_exchange_matr_complex(l_row(i_basis),l_col(j),1,i_spin) = &
                                  hf_exchange_matr_complex(l_row(i_basis),l_col(j),1,i_spin) + fock_tmp(j,i)*conjg(cfact)
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_complex_SR(l_row(i_basis),l_col(j),1,i_spin) = &
                                          hf_exchange_matr_complex_SR(l_row(i_basis),l_col(j),1,i_spin) + fock_tmp_SR(j,i)*conjg(cfact)
                                end if
                              endif
                            endif
                          enddo
                        endif
                      enddo

                    else

                      i_k = 0
                      do i_k_point = 1, n_k_points
                        full_k_point=map_back(i_k_point)
                        if (myid == MOD(i_k_point, n_tasks)) then
                          i_k = i_k + 1
                          !full_k_point=k_per_task(i_k,1)
                          do i = 1, atom2basis_len(i_atom)
                            i_basis = atom2basis_off(i_atom) + i
                            if(real_eigenvectors)then
                              hf_exchange_matr_real(:,i_basis,i_k,i_spin) = &
                                hf_exchange_matr_real(:,i_basis,i_k,i_spin) &
                                + dble(fock_tmp(:,i)*k_phase_exx_nosym(i_cell_fock,full_k_point)*0.5*n_spin)
                              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_real_SR(:,i_basis,i_k,i_spin) = &
                                            hf_exchange_matr_real_SR(:,i_basis,i_k,i_spin) &
                                            + dble(fock_tmp_SR(:,i)*k_phase_exx_nosym(i_cell_fock,full_k_point)*0.5*n_spin)
                              end if
                            else
                              hf_exchange_matr_complex(:,i_basis,i_k,i_spin) = &
                                hf_exchange_matr_complex(:,i_basis,i_k,i_spin) &
                                + fock_tmp(:,i)*k_phase_exx_nosym(i_cell_fock,full_k_point)*0.5*n_spin
                              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                          hf_exchange_matr_complex_SR(:,i_basis,i_k,i_spin) = &
                                            hf_exchange_matr_complex_SR(:,i_basis,i_k,i_spin) &
                                            + fock_tmp_SR(:,i)*k_phase_exx_nosym(i_cell_fock,full_k_point)*0.5*n_spin
                              end if
                            endif
                          enddo
                        endif
                      enddo

                    endif ! use_scalapack

                  enddo ! i_spin
        enddo
      enddo  
    else
      do i_atom = 1, n_atoms
        i_cell_fock = 1
        do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
          do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
             do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2
                                     
                if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then   
                  i_cell_fock = i_cell_fock + 1
                  ind = i_cell_fock
                else
                  ind = 1
                endif
                do i_spin = 1, n_spin

                    if(my_basis_off(i_atom) >= 0) then
                      my_s = my_basis_off(i_atom) + 1
                      my_e = my_basis_off(i_atom) + atom2basis_len(i_atom)
                      fock_tmp(:,1:atom2basis_len(i_atom)) = fock_matrix(:,my_s:my_e,ind,i_spin)
                      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          fock_tmp_SR(:,1:atom2basis_len(i_atom)) = fock_matrix_SR(:,my_s:my_e,ind,i_spin)
                      end if
                    else
                      fock_tmp(:,1:atom2basis_len(i_atom)) = 0
                      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          fock_tmp_SR(:,1:atom2basis_len(i_atom)) = 0
                      end if
                    endif

                    call sync_vector(fock_tmp, n_basis*atom2basis_len(i_atom))
                    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          call sync_vector(fock_tmp_SR, n_basis*atom2basis_len(i_atom))
                    end if

                    if(use_scalapack) then

                      cfact = k_phase_base(1, my_k_point)**i_cell_1 &
                          * k_phase_base(2, my_k_point)**i_cell_2 &
                          * k_phase_base(3, my_k_point)**i_cell_3*0.5*n_spin

                      do i = 1, atom2basis_len(i_atom)
                        i_basis = atom2basis_off(i_atom) + i
                        if(l_col(i_basis) > 0) then
                          do j = 1, n_basis
                            if(l_row(j)>0) then
                              if(real_eigenvectors)then
                                hf_exchange_matr_real(l_row(j),l_col(i_basis),1,i_spin) = &
                                  hf_exchange_matr_real(l_row(j),l_col(i_basis),1,i_spin) + dble(fock_tmp(j,i)*cfact)
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_real_SR(l_row(j),l_col(i_basis),1,i_spin) = &
                                          hf_exchange_matr_real_SR(l_row(j),l_col(i_basis),1,i_spin) + dble(fock_tmp_SR(j,i)*cfact)
                                end if
                              else
                                hf_exchange_matr_complex(l_row(j),l_col(i_basis),1,i_spin) = &
                                  hf_exchange_matr_complex(l_row(j),l_col(i_basis),1,i_spin) + fock_tmp(j,i)*cfact
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_complex_SR(l_row(j),l_col(i_basis),1,i_spin) = &
                                          hf_exchange_matr_complex_SR(l_row(j),l_col(i_basis),1,i_spin) + fock_tmp_SR(j,i)*cfact
                                end if
                              endif
                            endif
                          enddo
                        endif

                        if(l_row(i_basis) > 0) then
                          do j = 1, n_basis
                            if(l_col(j)>0) then
                              if(real_eigenvectors)then
                                hf_exchange_matr_real(l_row(i_basis),l_col(j),1,i_spin) = &
                                  hf_exchange_matr_real(l_row(i_basis),l_col(j),1,i_spin) + dble(fock_tmp(j,i)*conjg(cfact))
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_real_SR(l_row(i_basis),l_col(j),1,i_spin) = &
                                          hf_exchange_matr_real_SR(l_row(i_basis),l_col(j),1,i_spin) + dble(fock_tmp_SR(j,i)*conjg(cfact))
                                end if
                              else
                                hf_exchange_matr_complex(l_row(i_basis),l_col(j),1,i_spin) = &
                                  hf_exchange_matr_complex(l_row(i_basis),l_col(j),1,i_spin) + fock_tmp(j,i)*conjg(cfact)
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_complex_SR(l_row(i_basis),l_col(j),1,i_spin) = &
                                          hf_exchange_matr_complex_SR(l_row(i_basis),l_col(j),1,i_spin) + fock_tmp_SR(j,i)*conjg(cfact)
                                end if
                              endif
                            endif
                          enddo
                        endif
                      enddo

                    else

                      i_k = 0
                      do i_k_point = 1, n_k_points
                        if (myid == MOD(i_k_point, n_tasks)) then
                          i_k = i_k + 1
                          k_phase_new = k_phase_base(1, i_k_point)**i_cell_1 &
                                      * k_phase_base(2, i_k_point)**i_cell_2 &
                                      * k_phase_base(3, i_k_point)**i_cell_3
                          do i = 1, atom2basis_len(i_atom)
                            i_basis = atom2basis_off(i_atom) + i
                            if(real_eigenvectors)then
                              hf_exchange_matr_real(:,i_basis,i_k,i_spin) = &
                                hf_exchange_matr_real(:,i_basis,i_k,i_spin) &
                                + dble(fock_tmp(:,i)*k_phase_new*0.5*n_spin)
                              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_real_SR(:,i_basis,i_k,i_spin) = &
                                            hf_exchange_matr_real_SR(:,i_basis,i_k,i_spin) &
                                            + dble(fock_tmp_SR(:,i)*k_phase_new*0.5*n_spin)
                              end if
                            else
                              hf_exchange_matr_complex(:,i_basis,i_k,i_spin) = &
                                hf_exchange_matr_complex(:,i_basis,i_k,i_spin) &
                                + fock_tmp(:,i)*k_phase_new*0.5*n_spin
                              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                          hf_exchange_matr_complex_SR(:,i_basis,i_k,i_spin) = &
                                            hf_exchange_matr_complex_SR(:,i_basis,i_k,i_spin) &
                                            + fock_tmp_SR(:,i)*k_phase_new*0.5*n_spin
                              end if
                            endif
                          enddo
                        endif
                      enddo

                    endif ! use_scalapack

                  enddo ! i_spin

              enddo
          enddo
        enddo      
      enddo
    endif
  endsubroutine FT_hf_exchange_matr_sym
  !----------------------------------------------------------------------------
  !****s* sym_base/FT_k_phase_sym
  !  NAME
  !    FT_k_phase_sym
  !  SYNOPSIS  
subroutine FT_k_phase_sym(Cvec, c_phases_fft)
 
    !  PURPOSE
    !   Fourier transformation of k-phase
    !     For symmetry we have to use full (inversion symmetry) k-point set
    !
    !   This version circumvents use of k_phase_exx
    !
    !  USES
    use dimensions,only:n_k_points_nosym, n_atoms
    use runtime_choices, only: n_k_points_xyz_nosym
    use geometry, only: recip_lattice_vector
    use pbc_lists, only: n_cells_bvk
    use mpi_tasks, only: myid, n_tasks, aims_stop
    use synchronize_mpi_basic,only:sync_vector_complex
    implicit none

    !  ARGUMENTS
    real*8, dimension(3):: Cvec
    complex*16 :: c_phases_fft(n_cells_bvk)
    !  INPUTS
    !   o Cvec - periodic (real space) replica
    !
    !  OUTPUT
    !   o c_phases_fft - Fourier transforme phase
    !
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
    !

    ! local variables
    integer :: info, i_cell_1, i_cell_2, i_cell_3, ind, i_spin, i_k_point, &
              i_cell
    real*8 :: phi, qvec(3)
    complex*16 :: c_phases_nosym(n_k_points_nosym), phase
              


                ! If reconstruct_proper_only .false. not set we do not reconstruct fully
                ! but till inversion symmetry (k_point_list_nosym and k-_weights_nosym
                ! differ here)
                do i_k_point = 1, n_k_points_nosym
                      qvec(:) = matmul(recip_lattice_vector, k_point_list_nosym(i_k_point,:))
                      phi = dot_product(qvec(:), Cvec)
                      c_phases_nosym(i_k_point) = cmplx(cos(phi), sin(phi), kind(0.d0))
                end do 
                
                c_phases_fft(:) = 0
                    i_cell = 1
                    ! Next three dos and if clause circumvent storing of k_phase_exx,
                    ! which is a monster (n_points x n_points) for large k-point sets.
                    ! But we do the three exponentials.
                    do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
                      do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                        do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2                                     
                          if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then
                            i_cell = i_cell + 1
                            ind = i_cell

                          else
                            ! Gamma point is number 1
                            ind = 1
                          endif ! Gamma point
                          do i_k_point = 1, n_k_points_nosym
                            phase =  conjg(k_phase_base_nosym(1, i_k_point)**i_cell_1 &
                                         * k_phase_base_nosym(2, i_k_point)**i_cell_2 &
                                       * k_phase_base_nosym(3, i_k_point)**i_cell_3)
                            c_phases_fft(ind) = c_phases_fft(ind) &
                            + c_phases_nosym(i_k_point) * phase * k_weights_nosym(i_k_point) 
			  enddo ! i_k_point
                        enddo ! i_cell_1
                      enddo ! i_cell_2
                    enddo ! i_cell_3
endsubroutine FT_k_phase_sym
  !----------------------------------------------------------------------------
  !****s* sym_base/FT_k_phase_sym
  !  NAME
  !    FT_k_phase_sym
  !  SYNOPSIS  
subroutine FT_k_phase_sym_dist(Cvec, c_phases_fft, my_atom_id, my_tasks_per_atom, &
			  mpi_comm_atom)
 
    !  PURPOSE
    !   Fourier transformation of k-phase
    !     For symmetry we have to use full (inversion symmetry) k-point set
    !
    !   This version circumvents use of k_phase_exx
    !
    !  USES
    use dimensions,only:n_k_points_nosym, n_atoms
    use runtime_choices, only: n_k_points_xyz_nosym
    use geometry, only: recip_lattice_vector
    use pbc_lists, only: n_cells_bvk
    use mpi_tasks, only: myid, n_tasks, aims_stop
    use synchronize_mpi_basic,only:sync_vector_complex
    implicit none

    !  ARGUMENTS
    real*8, dimension(3):: Cvec
    complex*16 :: c_phases_fft(n_cells_bvk)
    integer :: my_atom_id, my_tasks_per_atom, mpi_comm_atom
    !  INPUTS
    !   o Cvec - periodic (real space) replica
    !
    !  OUTPUT
    !   o c_phases_fft - Fourier transforme phase
    !
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
    !

    ! local variables
    integer :: info, i_cell_1, i_cell_2, i_cell_3, ind, i_spin, i_k_point, &
              i_cell
    real*8 :: phi, qvec(3)
    complex*16 :: c_phases_nosym(n_k_points_nosym), phase
              


                ! If reconstruct_proper_only .false. not set we do not reconstruct fully
                ! but till inversion symmetry (k_point_list_nosym and k-_weights_nosym
                ! differ here)
                do i_k_point = 1, n_k_points_nosym
                      qvec(:) = matmul(recip_lattice_vector, k_point_list_nosym(i_k_point,:))
                      phi = dot_product(qvec(:), Cvec)
                      c_phases_nosym(i_k_point) = cmplx(cos(phi), sin(phi), kind(0.d0))
                end do 
                
                c_phases_fft(:) = 0
                    i_cell = 1
                    ! Next three dos and if clause circumvent storing of k_phase_exx,
                    ! which is a monster (n_points x n_points) for large k-point sets.
                    ! But we do the three exponentials.
                    do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
                      do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                        do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2                                     
                          if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then
                            i_cell = i_cell + 1
                            ind = i_cell

                          else
                            ! Gamma point is number 1
                            ind = 1
                          endif ! Gamma point
                          if (my_atom_id == MOD(ind, my_tasks_per_atom)) then
                          !write(use_unit,*) myid, my_atom_id, my_tasks_per_atom
                          do i_k_point = 1, n_k_points_nosym
                            phase =  conjg(k_phase_base_nosym(1, i_k_point)**i_cell_1 &
                                         * k_phase_base_nosym(2, i_k_point)**i_cell_2 &
                                       * k_phase_base_nosym(3, i_k_point)**i_cell_3)
                            c_phases_fft(ind) = c_phases_fft(ind) &
                            + c_phases_nosym(i_k_point) * phase * k_weights_nosym(i_k_point) 
			  enddo ! i_k_point
			  endif
                        enddo ! i_cell_1
                      enddo ! i_cell_2
                    enddo ! i_cell_3
                    call sync_vector_complex(c_phases_fft,n_cells_bvk,mpi_comm_atom)
endsubroutine FT_k_phase_sym_dist
  !----------------------------------------------------------------------------
  !****s* sym_base/FT_k_phase_sym_exx
  !  NAME
  !    FT_k_phase_sym_exx
  !  SYNOPSIS  
subroutine FT_k_phase_sym_exx(Cvec, c_phases_fft)
 
    !  PURPOSE
    !   Fourier transformation of k-phase
    !     For symmetry we have to use full (inversion symmetry) k-point set
    !  USES
    use dimensions, only: n_k_points_nosym
    use runtime_choices, only: n_k_points_xyz_nosym
    use geometry, only: recip_lattice_vector
    use pbc_lists, only: n_cells_bvk, k_phase_exx
    use synchronize_mpi_basic,only:sync_vector_complex
    use mpi_tasks, only: myid, n_tasks, aims_stop, mpi_comm_global
    implicit none

    !  ARGUMENTS
    real*8, dimension(3):: Cvec
    complex*16 :: c_phases_fft(n_cells_bvk)
    

    !  INPUTS
    !   o Cvec - periodic (real space) replica
    !
    !  OUTPUT
    !   o c_phases_fft - Fourier transforme phase
    !
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
    !

    ! local variables
    integer :: info, i_spin, i_k_point, i_cell,mpierr
    real*8 :: phi, qvec(3)
    complex*16 :: c_phases_nosym(n_k_points_nosym)
              
                ! If reconstruct_proper_only .false. not set we do not reconstruct fully
                ! but till inversion symmetry (k_point_list_nosym and k-_weights_nosym
                ! differ here)
                do i_k_point = 1, n_k_points_nosym
                      qvec(:) = matmul(recip_lattice_vector, k_point_list_nosym(i_k_point,:))
                      phi = dot_product(qvec(:), Cvec)
                      c_phases_nosym(i_k_point) = cmplx(cos(phi), sin(phi), kind(0.d0))
                end do 
                
                c_phases_fft(:) = 0
                do i_k_point = 1, n_k_points_nosym

                    do i_cell = 1, n_cells_bvk
                          c_phases_fft(i_cell) = c_phases_fft(i_cell) &
                          + c_phases_nosym(i_k_point) * k_weights_nosym(i_k_point) &
                          * conjg(k_phase_exx_nosym(i_cell,i_k_point)) 
                    enddo ! i_cell
                enddo ! i_k_point  
endsubroutine FT_k_phase_sym_exx
  !----------------------------------------------------------------------------
  !****s* sym_base/FT_k_phase_sym_exx
  !  NAME
  !    FT_k_phase_sym_exx
  !  SYNOPSIS  
subroutine FT_k_phase_sym_exx_dist(Cvec, c_phases_fft, my_atom_id, my_tasks_per_atom, &
			  mpi_comm_atom)
 
    !  PURPOSE
    !   Fourier transformation of k-phase
    !     For symmetry we have to use full (inversion symmetry) k-point set
    !  USES
    use dimensions, only: n_k_points_nosym
    use runtime_choices, only: n_k_points_xyz_nosym
    use geometry, only: recip_lattice_vector
    use pbc_lists, only: n_cells_bvk, k_phase_exx
    use synchronize_mpi_basic,only:sync_vector_complex
    use mpi_tasks, only: myid, n_tasks, aims_stop, mpi_comm_global
    implicit none

    !  ARGUMENTS
    real*8, dimension(3):: Cvec
    complex*16 :: c_phases_fft(n_cells_bvk)
    integer :: my_atom_id, my_tasks_per_atom, mpi_comm_atom    

    !  INPUTS
    !   o Cvec - periodic (real space) replica
    !
    !  OUTPUT
    !   o c_phases_fft - Fourier transforme phase
    !
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
    !

    ! local variables
    integer :: info, i_spin, i_k_point, i_cell,mpierr
    real*8 :: phi, qvec(3)
    complex*16 :: c_phases_nosym(n_k_points_nosym)
              
                ! If reconstruct_proper_only .false. not set we do not reconstruct fully
                ! but till inversion symmetry (k_point_list_nosym and k-_weights_nosym
                ! differ here)
                do i_k_point = 1, n_k_points_nosym
                      qvec(:) = matmul(recip_lattice_vector, k_point_list_nosym(i_k_point,:))
                      phi = dot_product(qvec(:), Cvec)
                      c_phases_nosym(i_k_point) = cmplx(cos(phi), sin(phi), kind(0.d0))
                end do 
                
                c_phases_fft(:) = 0

		do i_cell = 1, n_cells_bvk
		    if (my_atom_id == MOD(i_cell, my_tasks_per_atom)) then
                      do i_k_point = 1, n_k_points_nosym                    
                          c_phases_fft(i_cell) = c_phases_fft(i_cell) &
                          + c_phases_nosym(i_k_point) * k_weights_nosym(i_k_point) &
                          * conjg(k_phase_exx_nosym(i_cell,i_k_point)) 
                      enddo ! i_k_point  
                    endif
                enddo ! i_cell
                call sync_vector_complex(c_phases_fft,n_cells_bvk,mpi_comm_atom)
endsubroutine FT_k_phase_sym_exx_dist
  !----------------------------------------------------------------------------
  !****s* sym_base/calculate_realspace_coulomb_sym
  !  NAME
  !    calculate_realspace_coulomb_sym
  !  SYNOPSIS 
subroutine calculate_realspace_coulomb_sym(i_species_1, i_species_2, Dvec, &
                                         auxmat, d_auxmat, AS_d_auxmat, nbb1, &
                                         nbb2_s, nbb2_e,my_cm_cell_num,init_calc_deriv,&
                                         AS_init_stress,AS_l_index,AS_m_index,&
                                         my_cm_cell_start,my_cm_cell_inc, my_atom_id,&
                                         my_tasks_per_atom, mpi_comm_atom)

    !  PURPOSE
    !
    !  Does the same as sum_up_auxmat_qvecs_from_realspace() but only for 1 atom/atom pair
    !  and immediatly transforms the resulting matrices back to realspace.
    !  The combination of calculation and immediate backtransformation can be made much
    !  more efficiently than doing both steps separate.
    !
    !  Description of the original sum_up_auxmat_qvecs_from_realspace:
    !
    !    Calculate auxiliary matrix for a given set of q-points by a
    !    real-space sum.  This will only work if the radial interaction
    !    splines of this modules are "short-ranged", i.e. if either
    !    (.not. have_mp_far) [bare overlap, HSE, cutCb] or if prepared for an
    !    Ewald-like treatment by initialize_periodic_tb_auxmat().
    !
    !  Adapted from calculate_fock_matrix_p0 for k-point symmetry
    !   Difference: call to FT_k_phase_sym
    !
    !  USES
    use dimensions, only: n_k_points, n_periodic
    use pbc_lists, only: n_cells_bvk
    use runtime_choices, only: AS_components, use_k_phase_exx
    use bravais, only: get_n_supercells
    use mpi_tasks, only: myid, check_allocation, aims_stop
    use geometry, only: lattice_vector
    use prodbas, only: sp2n_basbas_sp
    use tight_binding_auxmat, only: have_mp_far, iarange_species, fast_calculate_tb_auxmat
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2, nbb1, nbb2_s, nbb2_e,&
                           my_cm_cell_start, my_cm_cell_inc, my_cm_cell_num
    real*8, intent(IN)  :: Dvec(3)
    real*8, intent(OUT) :: auxmat(nbb1, nbb2_s:nbb2_e, my_cm_cell_num)
    real*8, intent(OUT) :: d_auxmat(nbb1, nbb2_s:nbb2_e, 3, my_cm_cell_num)
    real*8, intent(OUT) :: AS_d_auxmat(nbb1, nbb2_s:nbb2_e, AS_components, my_cm_cell_num)
    logical, intent(IN) :: init_calc_deriv, AS_init_stress
    integer, intent(IN) :: AS_l_index(1:9), AS_m_index(1:9)
    integer, intent(IN) :: my_atom_id, my_tasks_per_atom, mpi_comm_atom
    !  INPUTS
    !    o i_species_1, i_species_2 -- Species of the two atoms involved.
    !    o Dvec -- Essentially coords(:, i_atom_2) - coords(:, i_atom_1).
    !    o nbb1 -- first dimension of auxmat/d_auxmat
    !    o nbb2_s, nbb2_e - limits of second dimension of auxmat/d_auxmat
    !    o further inputs to decouple routine from module calculate_fock_matrix_p0
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !    o d_auxmat -- derivative of auxmat with respect to Dvec,
    !                  calculated only if init_calc_deriv is set.
    !    o AS_d_auxmat -- derivative of auxmat with respect to strain,
    !                     calculated only if AS_init_stress is .true.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: n_spbb_1, n_spbb_2
    integer :: n_supercells(3), a1, a2, a3, i_k_point, i_cell, n_cnt, &
               info, AS_index
    real*8 :: Rvec(3), Cvec(3), phi, Rmax, Cmax, qvec(3)
    real*8, allocatable :: aux(:,:), d_aux(:,:,:)
    complex*16 :: c_phases(n_k_points), c_phases_fft(n_cells_bvk)
    logical :: dist_cell_sum
    character(*), parameter :: func = 'calculate_realspace_coulomb'

    if(my_tasks_per_atom>1) then
      dist_cell_sum=.True.
    else
      dist_cell_sum=.False.
    endif
    if (have_mp_far) then
       ! For the bare Coulomb potential, the radial interaction splines
       ! have to be fudged according to the Ewald procedure and have_mp_far
       ! has to be set to .false., afterwards.
       call aims_stop('Cannot use realspace sum with multipoles', func)
    end if

    n_spbb_1 = sp2n_basbas_sp(i_species_1)
    n_spbb_2 = sp2n_basbas_sp(i_species_2)

    allocate(aux(n_spbb_1,n_spbb_2), stat=info)
    call check_allocation(info, 'aux', func)

    if(init_calc_deriv) then
       allocate(d_aux(n_spbb_1,n_spbb_2,3), stat=info)
    else
       allocate(d_aux(1,1,1), stat=info)
    endif
    call check_allocation(info, 'd_aux', func)

    auxmat(:,:,:) = 0.
    if(init_calc_deriv) d_auxmat(:,:,:,:) = 0.
    if(AS_init_stress) AS_d_auxmat(:,:,:,:) = 0.0d0


    Rmax = iarange_species(i_species_1, i_species_2) 
    Cmax = Rmax + sqrt(sum(Dvec**2))
    n_supercells = 0   ! Important for n_periodic < 3.

    call get_n_supercells(n_periodic, lattice_vector, Cmax, n_supercells)
    do a1 = -n_supercells(1), n_supercells(1)
       do a2 = -n_supercells(2), n_supercells(2)
          do a3 = -n_supercells(3), n_supercells(3)
             Cvec = matmul(lattice_vector, (/a1, a2, a3/))
             Rvec = Dvec + Cvec
             if (sum(Rvec**2) > Rmax**2) cycle

             call fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec, aux, d_aux, init_calc_deriv)
             
             if (dist_cell_sum)then
	      if(use_k_phase_exx)then
		call FT_k_phase_sym_exx_dist(Cvec, c_phases_fft, my_atom_id, my_tasks_per_atom, &
				    mpi_comm_atom)
	      else
		call FT_k_phase_sym_dist(Cvec, c_phases_fft, my_atom_id, my_tasks_per_atom, &
				    mpi_comm_atom)
	      endif
	     else
	      if(use_k_phase_exx)then
		call FT_k_phase_sym_exx(Cvec, c_phases_fft)
	      else
		call FT_k_phase_sym(Cvec, c_phases_fft)
	      endif
             endif
             n_cnt = 0
             do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc
               n_cnt = n_cnt + 1
               auxmat(1:n_spbb_1,nbb2_s:nbb2_e,n_cnt) = &
                 auxmat(1:n_spbb_1,nbb2_s:nbb2_e,n_cnt) + aux(1:n_spbb_1,nbb2_s:nbb2_e)*dble(c_phases_fft(i_cell))
               if(init_calc_deriv) &
                 d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,1:3,n_cnt) = &
                   d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,1:3,n_cnt) + d_aux(1:n_spbb_1,nbb2_s:nbb2_e,1:3)*dble(c_phases_fft(i_cell))

               !FK: Multiply nuclear distances and derivatives with respect to nuclear coordinates
               if (AS_init_stress) then
                 do AS_index = 1, AS_components, 1
                   AS_d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,AS_index,n_cnt) = &
                     AS_d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,AS_index,n_cnt) + &
                     d_aux(1:n_spbb_1,nbb2_s:nbb2_e,AS_l_index(AS_index)) * &
                     dble(c_phases_fft(i_cell)) * Rvec(AS_m_index(AS_index))
                 end do
               end if
             enddo
          end do
       end do
    end do
    
  end subroutine calculate_realspace_coulomb_sym
  
subroutine determine_k_minus_q_list_spg &
            (kq_point_list, kpq_point_list)

    !  PURPOSE
    !
    !    mapping the k_minus_q and k_plus_q reciprocal grid point back to the 1st BZ.
    !   
    !
    !  USES
    use dimensions, only:   n_k_points_nosym,n_k_points
    use runtime_choices,only:  k_points_offset
    use localorb_io,only: localorb_info
    use pbc_lists, only: k_point_list
    use spglib_symmetry,only: map
    use mpi_tasks, only: n_tasks, aims_stop
    implicit none

    !  ARGUMENTS

       integer, intent(OUT) :: kq_point_list(n_k_points,n_k_points)
       integer, intent(OUT) :: kpq_point_list(n_k_points,n_k_points)

    !  INPUTS
    !    o all in pbc_list module
    !  OUTPUTS
    !    o kq_point_list :: maps the k-minus-q grid point back to the 1st Brillouin zone
    !    o kpq_point_list :: maps the k-plus-q grid point back to the 1st Brillouin zone
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales

    integer :: i_kq_point
    integer :: i_k_point, i_q_point
    integer :: i_periodic

    real*8 :: k_minus_q(3)
    real*8 :: k_plus_q(3)
    real*8 :: diff(3)
    real*8 :: k_diff_thr
    character*140 :: info_str
    character(*), parameter :: func = 'determine_k_minus_q_list_spg'
    integer :: k_points_test(n_k_points_nosym)
    k_diff_thr=min(dot_product(k_points_offset,k_points_offset),1.d-5)+1.d-10

    kq_point_list = 0
    kpq_point_list = 0
    n_kq_points=0
    k_points_test = 0
    do i_q_point = 1, n_k_points, 1
      do i_k_point = 1, n_k_points, 1

!         write(use_unit,*) i_q_point, i_k_point
         k_minus_q(:)=k_point_list(i_k_point,:)-k_point_list(i_q_point,:)
         k_plus_q(:)=k_point_list(i_k_point,:)+k_point_list(i_q_point,:)

! convert the k-q point back to the 1st Brillouin zone
         do i_periodic=1, 3, 1
           if(k_minus_q(i_periodic).lt.-1d-10) then
             k_minus_q(i_periodic) = k_minus_q(i_periodic)+1.d0
           endif
           if(k_plus_q(i_periodic).ge.1.d0) then
             k_plus_q(i_periodic) = k_plus_q(i_periodic)-1.d0
           endif
         enddo
!         write(use_unit,*) "k_minus_q new:", k_minus_q(:)

! check which point in the 1st BZ that the k-q point corresponds to
         do i_kq_point = 1, n_k_points_nosym, 1
            diff(:) = k_minus_q(:) - k_point_list_nosym(i_kq_point,:)
            if(dot_product(diff,diff).lt.k_diff_thr) then
              kq_point_list(i_k_point,i_q_point) = i_kq_point
              if (k_points_test(i_kq_point)==0)then
              n_kq_points=n_kq_points+1
              k_points_test(i_kq_point) = n_kq_points
              endif
              exit
            endif
         enddo
         do i_kq_point = 1, n_k_points_nosym, 1
            diff(:) = k_plus_q(:) - k_point_list_nosym(i_kq_point,:)
            if(dot_product(diff,diff).lt.k_diff_thr) then
              kpq_point_list(i_k_point,i_q_point) = map(i_kq_point)
              exit
	    endif
         enddo   
!         write(use_unit,*) "i_kq_point", i_kq_point,kq_point_list(i_k_point,i_q_point),kpq_point_list(i_k_point,i_q_point),i_k_point, i_q_point
         if(kq_point_list(i_k_point,i_q_point).le.0 .or. kpq_point_list(i_k_point,i_q_point).le.0 ) then

           if(n_k_points .eq. 1 .and. dot_product(k_point_list(1,:),k_point_list(1,:)).gt.1.e-12) then
              kq_point_list(1,1) = 1
              write(info_str, '(2X,A)') &
               & "Using single k point slightly away from the Gamma point."
              call localorb_info(info_str)  
 
           else
             call aims_stop(" Error in determining the k-q mesh, stop!", func)
           endif

         endif

      enddo
   enddo
   n_kq_points = n_k_points_nosym
   allocate(k_points_coul(n_kq_points))
   do i_q_point = 1, n_k_points_nosym, 1
     !if(k_points_test(i_q_point).ne.0)then
       !k_points_coul(k_points_test(i_q_point))=i_q_point
       k_points_coul(i_q_point)=i_q_point
     !endif
   enddo 
   !do i_q_point = 1, n_k_points, 1
   !   do i_k_point = 1, n_k_points, 1
        !write(use_unit,*) kq_point_list(i_k_point,i_q_point), k_points_test(kq_point_list(i_k_point,i_q_point)),k_points_coul(k_points_test(kq_point_list(i_k_point,i_q_point)))
   !     kq_point_list(i_k_point,i_q_point)=k_points_test(kq_point_list(i_k_point,i_q_point))
   !   enddo
  ! enddo
   n_kq_points_task=(n_kq_points-1)/n_tasks+1
!Dirty hack waiting for a final fix
   if (n_k_points.eq.1) then
      kq_point_list=1
      kpq_point_list=1
   end if
  end subroutine determine_k_minus_q_list_spg

  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/get_coulomb_matr_recip
  !  NAME
  !    get_coulomb_matr_recip
  !  SYNOPSIS

  subroutine get_coulomb_matr_recip_spg &
            (coulomb_matr_recip, output_info)

    !  PURPOSE
    !
    !    Fourier transfrom the LVL triple coefficients from real space (on a Bravais lattice)
    !    to reciprocal space.
    !   
    !
    !  USES

    use dimensions
    use prodbas
    use pbc_lists
    use geometry
    use tight_binding_auxmat
    use localorb_io
    use mpi_tasks, only: check_allocation, myid, n_tasks
    implicit none

    !  ARGUMENTS

       complex*16, intent(INOUT) :: coulomb_matr_recip(n_basbas,n_basbas,n_kq_points_task)
       integer :: output_info

    !  INPUTS
    !    o none : all are in module pbc_list.f90 or pbc_lists.f90
    !  OUTPUTS
    !    o coulomb_matr_recip: Coulomb matrix in reciprocal sapce
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    ! Local variabales
    real*8, allocatable :: k_lattvec(:,:) ! k-minus-q lattice vectors in reciprocal space

    integer :: i_k_point, i_k_point_local
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_coulomb_matr_recip'

    if(output_info == 1) then
      write(info_str,'(2X,A)') "Calculating the Coulomb matrix in reciprocal space ..."
      call localorb_info(info_str)
    endif

   allocate(k_lattvec(3,n_kq_points_task),stat=info)
   call check_allocation(info, 'k_lattvect', func)

    k_lattvec(:,:) = 0.d0
    do i_k_point = 1, n_kq_points, 1

        if(myid.eq.mod(i_k_point,n_tasks)) then
         
          i_k_point_local = (i_k_point-1)/n_tasks + 1
          k_lattvec(1:3,i_k_point_local) = matmul(recip_lattice_vector,k_point_list_nosym(k_points_coul(i_k_point),1:3))
          ! k_minus_q_point_is a copy of k_point_list
        endif

    enddo

  
    call get_qspace_auxmat(n_kq_points_task,k_lattvec,coulomb_matr_recip)

    if(allocated(k_lattvec)) then
      deallocate(k_lattvec)
    endif

  end subroutine get_coulomb_matr_recip_spg
!****s* FHI-aims/evaluate_exchange_matr_kspace_p0.f90
!  NAME evaluate_exchange_matr_kspace_p0
!   
!  SYNOPSIS

subroutine evaluate_exchange_matr_kspace_spg &
(KS_eigenvalue,KS_eigenvector,KS_eigenvector_complex,occ_numbers)

  !  PURPOSE
  !  Subroutine evaluate_exchange_matr_kspace_p0 evaluates the exact-exchange part of
  !  Hartree-Fock hamiltonian in a periodic system. The algorithm used here are based
  !  on the reciprocal space and localized resolution of identity (RI-LVL)
  !
  !  USES

  use dimensions
  use prodbas
  use hartree_fock
  use hartree_fock_p0
  use mpi_tasks
  use synchronize_mpi
  use pbc_lists
  use runtime_choices
  use constants
  use spglib_symmetry, only:map, map_sym
  use basis
  use localorb_io, only: localorb_info, use_unit
  use geometry, only: species

  implicit none

  !  ARGUMENTS

  real*8, dimension(n_states,n_spin,n_k_points_task) :: KS_eigenvalue
  real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
  complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
  real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers

  !  INPUTS
  !  o  occ_numbers -- real array,
  !       the occupation number of the electrons for each eigenstate and each spin
  !  o  KS_eigenvector -- real array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  o  KS_eigenvector_complex -- complex array,
  !            the eigenvector of the single-particle (KS/HF) self-consistent calculation
  !  OUTPUTS
  !  none
  !  the exact exchange matrix (the "hf_exchange_matr_complex" in the source code) is defined in MODULE
  !  hartree_fock_p0
  !
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

  complex*16, allocatable :: lvl_tricoeff_recip_tmp(:,:,:) ! LVL triple coefficients in k space
  complex*16, allocatable :: KS_eigenvector_q(:,:,:) ! KS eigenvector a a given (romote) k point
  complex*16, allocatable :: KS_eigenvector_q_rotated(:)
  complex*16, allocatable :: lvl_tricoeff_sum1(:,:,:) ! sum of LVL triple coefficient at k and q points (atom 1)
  complex*16, allocatable :: lvl_tricoeff_sum2(:,:,:) ! sum of LVL triple coefficient at k and q points (atom 2)
  complex*16, allocatable :: lvl_tricoeff(:,:,:) ! Full LVL triple coefficient, one basis index is transformed to
  complex*16, allocatable :: lvl_tricoeff_rotated(:,:,:)   ! (occupied) molecular orbital index 
  complex*16, allocatable :: lvl_tricoeff_rotated_tmp(:,:,:)
  complex*16, allocatable :: lvl_tricoeff_KS(:,:)
  complex*16, allocatable :: lvl_tricoeff_KS_rotated(:,:)
  complex*16, allocatable :: tmp_MO(:,:,:)
!  complex*16, allocatable :: coulomb_times_tricoeff(:,:) ! Coulomb matrix multiplies the LVL triple coefficient
  complex*16:: coulomb_times_tricoeff(n_basbas,n_basis) ! Coulomb matrix multiplies the LVL triple coefficient
  complex*16, allocatable :: exchange_matr_tmp(:,:) ! temparary exchange matrix per k_point
  complex*16, allocatable :: coulomb_matr_recip_tmp(:,:) !
  complex*16, allocatable :: ex_matr_complex(:,:) !
  real*8, allocatable :: ex_matr_real(:,:) !
  complex*16, dimension(:,:), allocatable:: Delta_temp
  real*8  :: ex_vect_real
  complex*16  :: ex_vect_complex
  real*8  :: en_exx

  integer :: info, mpierr
  character(*), parameter :: func = 'evaluate_exchange_matr_kspace_p0'
  character*150 :: info_str

  integer :: max_n_homo
  integer :: n_homo_k(n_k_points,n_spin)
! counter 
  integer i_k_point
  integer i_q_point
  integer i_kq_point
  integer i_k_point_local
  integer i_q_point_local
  integer i_state
  integer i_state_1
  integer i_basis_1
  integer i_spin
  integer id_root
  integer i_task
  integer, allocatable :: k_points_at_task(:,:)
  integer :: i_req(0:n_tasks-1)
  integer :: i_atom_1, i_atom_2
  integer :: i_species_1, i_species_2
  integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
  integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
  integer :: i_1, i_2,i,j,k, l
  integer :: i_q_point_nosym, is, full_q_point, q_task, i_basbas, i_basis
!stop
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  if(use_scalapack) call aims_stop('*** Periodic EXX in k-space must not be used if use_scalapack is in effect!')
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  n_homo_k(:,:) = 0
  do i_spin = 1, n_spin, 1
    do i_k_point = 1, n_k_points, 1
      do i_state = 1, n_states
        if(occ_numbers(i_state,i_spin,i_k_point) .gt. 0.0) then
         n_homo_k(i_k_point,i_spin) = i_state
        endif
      enddo
 !     if(myid.eq.0) then
 !       write(use_unit,*) i_spin, i_k_point, n_homo_k(i_k_point, i_spin), n_homo(i_spin)
 !     endif
    enddo
  enddo
  max_n_homo = max(maxval(n_homo_k(:,1), n_k_points), &
                   maxval(n_homo_k(:,n_spin), n_k_points))
 ! write(use_unit,*)"max_n_homo:", max_n_homo

!  max_n_homo = max(n_homo(1),n_homo(n_spin))

  allocate(lvl_tricoeff_recip_tmp(max_n_basbas_sp,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_recip_tmp', func)
  allocate(KS_eigenvector_q(n_basis,n_states,n_spin),stat=info) 
  call check_allocation(info, 'KS_eigenvector_q', func)
  allocate(KS_eigenvector_q_rotated(n_basis),stat=info) 
  call check_allocation(info, 'KS_eigenvector_q_rotated', func)
  allocate(lvl_tricoeff(n_basbas,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff', func)
  allocate(lvl_tricoeff_rotated_tmp(n_basbas,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_rotated_tmp', func)
  allocate(lvl_tricoeff_rotated(n_basbas,n_basis,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_rotated', func)
  allocate(lvl_tricoeff_KS(n_basbas,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_KS', func)  
  allocate(lvl_tricoeff_KS_rotated(n_basbas,n_basis),stat=info) 
  call check_allocation(info, 'lvl_tricoeff_KS_rotated', func)  
!  allocate(coulomb_times_tricoeff(n_basbas,n_basis),stat=info) 
!  call check_allocation(info, 'coulomb_times_tricoeff', func)
  allocate(exchange_matr_tmp(n_basis,n_basis),stat=info) 
  call check_allocation(info, 'exchange_matr_tmp', func)
  allocate(coulomb_matr_recip_tmp(n_basbas,n_basbas),stat=info) 
  call check_allocation(info, 'coulomb_matr_recip_tmp', func)

  if(real_eigenvectors) then
    allocate(ex_matr_real(n_basis,n_states),stat=info)
    call check_allocation(info, 'ex_matr_real', func)
    allocate(ex_matr_complex(1,1),stat=info)
  else
    allocate(ex_matr_complex(n_basis,n_states),stat=info)
    call check_allocation(info, 'ex_matr_complex', func)
    allocate(ex_matr_real(1,1),stat=info)
  endif

  allocate(Delta_temp(n_basis, n_basis), stat=info)
  call check_allocation(info, 'Delta_temp')
              
  write(info_str,'(2X,A)') "Constructing the exchange matrix in k space ..."
  call localorb_info(info_str)

  ! k_points_at_task: which task contains which k-point

  allocate(k_points_at_task(0:n_tasks-1,(n_k_points-1)/n_tasks+1))
  k_points_at_task(:,:) = 0
  do i_k_point = 1, n_k_points, 1
     i_k_point_local = (i_k_point-1)/n_tasks + 1 
     k_points_at_task(mod(i_k_point,n_tasks),i_k_point_local) = i_k_point
  enddo

  i_req(:) = MPI_REQUEST_NULL

  if(real_eigenvectors) then
    hf_exchange_matr_real(:,:,:,:) = 0.d0
  else
    hf_exchange_matr_complex(:,:,:,:) = (0.d0,0.d0)
  endif
  en_exx = 0.d0
open(unit=8, file='lvl_tricoeff.dat',ACTION='WRITE') 
  do i_q_point = 1, n_k_points, 1
!     if(myid.eq.1) print*,i_q_point
     id_root = mod(i_q_point,n_tasks)
     i_q_point_local = (i_q_point-1)/n_tasks + 1
     

     if(myid .eq. id_root) then

        lvl_tricoeff_recip_tmp(:,:,:)=lvl_tricoeff_recip1(:,:,:,i_q_point_local)

        if(real_eigenvectors) then
           KS_eigenvector_q(:,:,:)=dble(KS_eigenvector(:,:,:,i_q_point_local))
        else
           KS_eigenvector_q(:,:,:)=KS_eigenvector_complex(:,:,:,i_q_point_local)
        endif
     endif

     call mpi_bcast(lvl_tricoeff_recip_tmp,max_n_basbas_sp*n_basis*n_basis, &
                             MPI_COMPLEX16, id_root, mpi_comm_global, mpierr)

     call mpi_bcast(KS_eigenvector_q,n_basis*n_states*n_spin, &
                             MPI_COMPLEX16, id_root, mpi_comm_global, mpierr)
     KS_eigenvector_q = conjg(KS_eigenvector_q)


     do i_k_point_local = 1, (n_k_points-1)/n_tasks + 1

        ! Send our coulomb_matr_recip to the task(s) which actually needs it
        ! Normally, only one send should be necessary

        do i_task = 0, n_tasks-1
           i_k_point = k_points_at_task(i_task, i_k_point_local)
           if(i_k_point == 0) cycle
           i_kq_point = kq_point_list(i_k_point,i_q_point)
           if(mod(i_kq_point,n_tasks) == myid) then
             ! Tasks i_task needs i_kq_point from me
             call mpi_isend(coulomb_matr_recip(1,1,(i_kq_point-1)/n_tasks+1), n_basbas*n_basbas, MPI_COMPLEX16, &
                            i_task, 111, mpi_comm_global, i_req(i_task), mpierr)
           endif
        enddo

        ! Our k-point in this turn:
        i_k_point = k_points_at_task(myid, i_k_point_local)

        if(i_k_point == 0) then
           call mpi_waitall(n_tasks, i_req, MPI_STATUSES_IGNORE, mpierr)
           cycle
        endif

        i_kq_point = kq_point_list(i_k_point,i_q_point)
!write(use_unit,*) i_kq_point, k_points_coul(i_kq_point)
        call mpi_recv(coulomb_matr_recip_tmp, n_basbas*n_basbas, MPI_COMPLEX16, &
                     mod(i_kq_point, n_tasks), 111, mpi_comm_global, MPI_STATUS_IGNORE, mpierr)

        call mpi_waitall(n_tasks, i_req, MPI_STATUSES_IGNORE, mpierr)

		   
 do full_q_point = 1, n_k_per_irr(i_q_point)
                   q_task = k_per_irr(i_q_point,full_q_point) 
                   is = map_sym(q_task)
                   if(use_spg_full_Delta)then
                     Delta_temp = Delta_matrix_full(:,:, is,i_q_point)
                   else
                     call get_density_rotation(is, i_q_point, &
			 map_atom(:,is), Delta(:,:,:,is), Delta_temp)
                   endif   
        lvl_tricoeff = 0                 
        do i_atom_1 = 1, n_atoms, 1
           i_species_1 = species(i_atom_1)
           basis_off_1 = atom2basis_off(i_atom_1)
           n_sp_basis_1 = sp2n_basis_sp(i_species_1)
           bboff_1 = atom2basbas_off(i_atom_1)
           n_spbb_1 = sp2n_basbas_sp(i_species_1)

           allocate(lvl_tricoeff_sum1(n_spbb_1,n_sp_basis_1,n_basis),stat=info)
           call check_allocation(info, 'lvl_tricoeff_sum1', func)

           do i_1 = 1, n_sp_basis_1
              do i_2 = 1, n_basis
                 lvl_tricoeff_sum1(1:n_spbb_1,i_1,i_2) = &
                    conjg(lvl_tricoeff_recip_tmp(1:n_spbb_1, basis_off_1+i_1, i_2)) + &
                    lvl_tricoeff_recip2(1:n_spbb_1, i_2, basis_off_1+i_1)
              enddo
           enddo

	   lvl_tricoeff(bboff_1+1:bboff_1+n_spbb_1,basis_off_1+1:basis_off_1+n_sp_basis_1,:) = &
	    lvl_tricoeff_sum1

           deallocate(lvl_tricoeff_sum1)
        enddo

        do i_atom_2 = 1, n_atoms, 1
           i_species_2 = species(i_atom_2)
           basis_off_2 = atom2basis_off(i_atom_2)
           n_sp_basis_2 = sp2n_basis_sp(i_species_2)
           bboff_2 = atom2basbas_off(i_atom_2)
           n_spbb_2 = sp2n_basbas_sp(i_species_2)

           allocate(lvl_tricoeff_sum2(n_spbb_2,n_basis,n_sp_basis_2),stat=info)
           call check_allocation(info, 'lvl_tricoeff_sum2', func)

           do i_1 = 1, n_basis
              do i_2 = 1, n_sp_basis_2
                 lvl_tricoeff_sum2(1:n_spbb_2,i_1,i_2) = &
                    lvl_tricoeff_recip1(1:n_spbb_2, basis_off_2+i_2, i_1, i_k_point_local) + &
                    lvl_tricoeff_recip2(1:n_spbb_2, i_1, basis_off_2+i_2)
              enddo
           enddo

	   lvl_tricoeff(bboff_2+1:bboff_2+n_spbb_2,:,basis_off_2+1:basis_off_2+n_sp_basis_2) = &
	   lvl_tricoeff(bboff_2+1:bboff_2+n_spbb_2,:,basis_off_2+1:basis_off_2+n_sp_basis_2) + lvl_tricoeff_sum2


           deallocate(lvl_tricoeff_sum2)
        enddo

		   !call zgemm ('N', 'T', n_basbas*n_basis, n_basis, n_basis, (1.d0,0.d0), &
                   !       lvl_tricoeff, n_basbas*n_basis,Delta_temp,n_basis,(0.d0,0.d0),&
                   !       lvl_tricoeff_rotated, n_basbas*n_basis)        
                    do i_basbas=1,max_n_basbas_sp, 1
                     do i_basis=1,n_basis, 1
                     	call zgemv('N', n_basis,  n_basis, dcmplx(1.d0,0.d0), Delta_temp, n_basis, &
                    &          lvl_tricoeff(i_basbas,i_basis,:), 1, dcmplx(0.d0,0.d0),&
                               lvl_tricoeff_rotated(i_basbas,i_basis,:), 1 )     
                     enddo
                     !do i_basis=1,n_basis, 1
                     !	call zgemv('N', n_basis,  n_basis, dcmplx(1.d0,0.d0), Delta_temp, n_basis, &
                    !&          lvl_tricoeff_recip2(i_basbas,i_basis,:), 1, dcmplx(0.d0,0.d0),&
                    !           lvl_tricoeff_rotated(i_basbas,i_basis,:), 1 )         
                     !enddo            
                   enddo  
    
        do i_spin = 1, n_spin, 1
!           !$OMP PARALLEL DO private(i_state,coulomb_times_tricoeff,exchange_matr_tmp)

             
           do i_state = 1, 7, 1
                                
	      call perfon('eemk')

	           call zgemv('N', n_basis,  n_basis, dcmplx(1.d0,0.d0), Delta_temp, n_basis, &
                    &          KS_eigenvector_q(:, i_state, i_spin), 1, dcmplx(0.d0,0.d0),&
                               KS_eigenvector_q_rotated, 1 )   
		   call zgemv('N', n_basbas*n_basis,  n_basis, dcmplx(1.d0,0.d0), lvl_tricoeff_rotated, n_basbas*n_basis, &
                    &          (KS_eigenvector_q_rotated), 1, dcmplx(0.d0,0.d0),lvl_tricoeff_KS_rotated, 1 ) 
		   !call zgemm ('N', 'N', n_basbas, n_basis, n_basis, (1.d0,0.d0), &
                    !      lvl_tricoeff_KS, n_basbas,Delta_temp,n_basis,(0.d0,0.d0),&
                    !      lvl_tricoeff_KS_rotated, n_basbas)
if(i_k_point==1)then                    
do i = 1, n_basbas,1
do j = 1, n_basis,1
!do l = 1, n_basis,1
write(8,'(5I5,2F15.10)') q_task,map_back(i_k_point),i,j,i_state,lvl_tricoeff_KS_rotated(i,j)
!enddo
enddo
enddo
endif

                 call zgemm('N', 'N', n_basbas, n_basis, n_basbas, (1.d0,0.d0), &
                         coulomb_matr_recip_tmp, n_basbas, &
                         lvl_tricoeff_KS_rotated, n_basbas, (0.d0,0.d0), &
                         coulomb_times_tricoeff, n_basbas)

                 call zgemm('C', 'N', n_basis, n_basis, n_basbas, (1.d0,0.d0), &
                          lvl_tricoeff_KS_rotated, n_basbas, &
                          coulomb_times_tricoeff, n_basbas, (0.d0,0.d0), &
                          exchange_matr_tmp,n_basis)
	      call perfoff	
                 

              if(real_eigenvectors) then
                 hf_exchange_matr_real(:,:,i_k_point_local,i_spin) = &
                   hf_exchange_matr_real(:,:,i_k_point_local,i_spin) + &
                   exchange_matr_tmp(:,:)*k_weights(i_q_point)* & 
                   occ_numbers(i_state,i_spin,i_q_point)*dble(n_spin)/2.d0

                   call dgemm('N', 'N', n_basis, n_states, n_basis, 1.d0, &
                              real(exchange_matr_tmp), n_basis, &
                              KS_eigenvector(:,:,i_spin,i_k_point_local), n_basis, 0.d0, &
                              ex_matr_real,n_basis)

                   do i_state_1 = 1, max_n_homo, 1
                      ex_vect_real = 0.d0
                      do i_basis_1 = 1, n_basis, 1
                        ex_vect_real = ex_vect_real + &
                                       KS_eigenvector(i_basis_1,i_state_1,i_spin,i_k_point_local) &
                                     * ex_matr_real(i_basis_1,i_state_1)
                      enddo
                      if(KS_eigenvalue(i_state,i_spin,i_q_point) .gt.  &
                          KS_eigenvalue(i_state_1,i_spin,i_k_point) ) then
                         en_exx = en_exx + ex_vect_real * occ_numbers(i_state,i_spin,i_q_point) * &
                                    k_weights(i_q_point)*k_weights(i_k_point) * dble(n_spin)
                               

                      elseif(KS_eigenvalue(i_state,i_spin,i_q_point) .eq.  KS_eigenvalue(i_state_1,i_spin,i_k_point))  then
                            en_exx = en_exx + ex_vect_real * occ_numbers(i_state,i_spin,i_q_point) * &
                                    k_weights(i_q_point)*k_weights(i_k_point) *dble(n_spin)/2.d0
                      endif
                                               
                   enddo

              else
               
                   hf_exchange_matr_complex(:,:,i_k_point_local,i_spin) = &
                      hf_exchange_matr_complex(:,:,i_k_point_local,i_spin) + &
                      exchange_matr_tmp(:,:)*k_weights_nosym(q_task)* & 
                      occ_numbers(i_state,i_spin,i_q_point)*dble(n_spin)/2.d0

                   call zgemm('N', 'N', n_basis, n_states, n_basis, (1.d0,0.d0), &
                              exchange_matr_tmp, n_basis, &
                              KS_eigenvector_complex(:,:,i_spin,i_k_point_local), n_basis, (0.d0,0.d0), &
                              ex_matr_complex,n_basis)

                   do i_state_1 = 1, max_n_homo, 1
                      ex_vect_complex = (0.d0,0.d0)
                      do i_basis_1 = 1, n_basis, 1
                        ex_vect_complex = ex_vect_complex + &
                                 conjg(KS_eigenvector_complex(i_basis_1,i_state_1,i_spin,i_k_point_local)) &
                                 * ex_matr_complex(i_basis_1,i_state_1) 
                      enddo
                      if(KS_eigenvalue(i_state,i_spin,map(q_task)) .gt. &
                          KS_eigenvalue(i_state_1,i_spin,i_k_point) ) then
                         en_exx = en_exx + ex_vect_complex * occ_numbers(i_state,i_spin,map(q_task)) * &
                                    k_weights_nosym(q_task)*k_weights(i_k_point) * dble(n_spin)

                      elseif(KS_eigenvalue(i_state,i_spin,map(q_task)) ==   &
                                 KS_eigenvalue(i_state_1,i_spin,i_k_point) ) then
                         en_exx = en_exx + ex_vect_complex * occ_numbers(i_state,i_spin,map(q_task)) * &
                                    k_weights_nosym(q_task)*k_weights(i_k_point) * dble(n_spin)/2.d0 
                      endif
                                               
                   enddo
              endif

! end loop over i_state
           enddo
!           !$OMP END PARALLEL DO
! end loop over i_spin
        enddo
! end loop over i_k_point
     enddo
enddo           
! end loop i_q_point
  enddo
  close(unit=8)
do i = 1, n_basis,1
do j = 1, n_basis,1
write(use_unit,'(2I5,2F15.10)') i,j,hf_exchange_matr_complex(i,j,2,1)
enddo
enddo

  call sync_real_number(en_exx)

  if(myid.eq.0) then
   write(use_unit,*) "EN_EXX:", en_exx, en_exx*hartree
  endif


  deallocate(lvl_tricoeff_recip_tmp)
  deallocate(KS_eigenvector_q)
  deallocate(KS_eigenvector_q_rotated)
  deallocate(lvl_tricoeff)
  deallocate(lvl_tricoeff_rotated)
  deallocate(lvl_tricoeff_KS_rotated)  
!  deallocate(coulomb_times_tricoeff)
  deallocate(exchange_matr_tmp)
  deallocate(coulomb_matr_recip_tmp)
  deallocate(ex_matr_real)
  deallocate(ex_matr_complex)

  deallocate(k_points_at_task)

  return

end subroutine evaluate_exchange_matr_kspace_spg
!---------------------------------------------------------------------
!******

  
subroutine symmetrize_forces_spg(forces)

  !  PURPOSE
  !
  !    Symmetrize the (total) forces by averaging over all symmetry operations.
  !
  !  USES
  use spglib_symmetry, only: spg_rotations, num_symmetry
  use dimensions, only: n_atoms
  use localorb_io, only: localorb_info, use_unit
  implicit none

  !  ARGUMENTS
  !    o forces ({x,y,z},n_atoms) forces of the system
  !  INPUTS
  real*8, dimension(3,n_atoms), intent(inout) :: forces
  !  OUTPUTS
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

  ! Local Variables
  real*8, dimension(3,n_atoms) :: forces_work
  real*8, dimension(3) :: work
  integer              :: i_atom, j_atom, is
  character(len=120)   :: info_str


  write(info_str,'(2X,A)') "| Symmetrizing total forces. "
  call localorb_info(info_str,use_unit,'(A)')
  forces_work = 0.d0
  do is = 1, num_symmetry, 1

    do i_atom = 1, n_atoms, 1
      
       j_atom = map_atom_full(i_atom,is)
       !temp = inv(dble(spg_rotations(1:3,1:3,is)))
       work = matmul(spg_rotations(:,:,is),forces(:,j_atom))
       forces_work(:,i_atom) =  forces_work(:,i_atom) + work
    enddo

  enddo

  forces = forces_work/num_symmetry

end subroutine symmetrize_forces_spg

subroutine symmetrize_stress_spg(stress)

  !  PURPOSE
  !
  !    Symmetrize the stress by averaging over all symmetry operations,
  !    similar to force symmetrization above
  !
  !  USES
  use spglib_symmetry, only: spg_rotations, num_symmetry
  use dimensions     , only: n_atoms
  use localorb_io    , only: localorb_info, use_unit
  implicit none

  !  ARGUMENTS
  !    o stress ({xx, yx, zx, xy, yy, zy, xz, yz, zz}) stress tensor of
  !    the system
  !  INPUTS
  real*8, dimension(3, 3), intent(inout) :: stress
  !  OUTPUTS
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

  ! Local Variables
  real*8, dimension(3, 3)     :: stress_work
  real*8, dimension(3, 3)     :: work1, work2, unit_mat, rot, rotT, rotTrot
  real*8                      :: count
  integer                     :: i_atom, j_atom, is, ii, jj
  character(len=120)          :: info_str
  logical                     :: isrot

  write(info_str,'(2X,A)') "| Symmetrizing stress "
  call localorb_info(info_str,use_unit,'(A)')
  stress_work = 0.d0

  ! define unit matrix
  unit_mat(1,1:3) =  (/ 1.d0, 0.d0, 0.d0 /)
  unit_mat(2,1:3) =  (/ 0.d0, 1.d0, 0.d0 /)
  unit_mat(3,1:3) =  (/ 0.d0, 0.d0, 1.d0 /)

  count = 0.d0 ! count the number of symmetry transformations applied
  do is = 1, num_symmetry, 1
     ! filter real rotations
     rot         = spg_rotations(:, :, is)
     rotT        = transpose(rot)
     rotTrot = matmul(rot, rotT)
     isrot = .true.
     if (maxval(abs(rotTrot - unit_mat)) > 1d-10) isrot = .false.
     if (isrot) then
        work1       = matmul(rot, stress)
        work2       = matmul(work1, rotT)
        stress_work = stress_work + work2
        count       = count + 1.d0
     end if
  enddo

  stress = stress_work / count
  !write(use_unit,*) stress
end subroutine symmetrize_stress_spg
end module sym_base
