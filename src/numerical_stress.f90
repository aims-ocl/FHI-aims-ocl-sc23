!****h* FHI-aims/numerical_stress
!  NAME
!   numercal_stress 
!  SYNOPSIS

module numerical_stress 

!  PURPOSE
!
!    This module contains all the routines for 
!    calculating the numerical stress
!
!    Because the numerical stress is based on finite differences
!    it's  calculated in predict_new_geometry.f90 

!  USES

  use localorb_io
  use dimensions
  use constants
  use geometry
  use pbc_lists
  use physics

!  ARGUMENTS

  implicit none
 
  real*8                               :: delta_numerical_stress, default_delta_numerical_stress =0.0001d0
  integer                              :: counter_numerical_stress = 0
  real*8, dimension(:,:), allocatable  :: original_coords 
  real*8, dimension(:,:), allocatable  :: original_lattice_vector
  real*8                               :: original_volume
  real*8, dimension(:,:), allocatable  :: original_periodic_unit_cell_translations 
  character*20                         :: original_output_level 
  logical                              :: original_force_flag 
  logical                              :: original_analytical_stress_flag
  logical                              :: original_skip_scf_flag
  real*8, dimension (3,3,0:1)          :: numerical_stress_tensor_components=0d0

  !VA: moved to physics.90
  !real*8,        dimension (3,3)       :: numerical_stress_tensor=0d0
  !real*8                               :: numerical_pressure



!  INPUTS
!    
!  OUTPUT
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
! 
!-----------------------------------------------------------------------------------------------------------------
! MATHEMATICAL BACKGROUND:
!
! The stress-tensor T _ij is obtained in two steps:
!
!      (1) Apply the  transformation matrix (1+epsilon) on the Hamiltonian. 
!          epsilon: 'strain matrix'
!
!      (2) Derive the Free energy F with respect to the strain matrix epsilon
!           
!                    T_ij=d/d epsilon_ij F
!
!
! In principle (1+epsilon) acts on the nuclear and electronic configuration space. 
! However since here the numerical stress is computed by finited displacements, it 
! effectively only acts on the nuclear configuration space (The electrons are then
! 'distributed automatically' by the scf).
!
! Here the parametriztion of the symmetric strain matrix epsilon is according to 
! REFERENCE:         K.Doll, Molecular Physics, Vol 108, 223-227 [doll]
! FK: We remove the factor 1/2 for the off-diagonal elements because otherwise they are 
!     divided two times by 2 which is incorrect. One division by 2 is correct.
!
!    (      x_1      1/2 x_6      1/2 x_5  )
!    (  1/2 x_6          x_2      1/2 x_4  )  =:  epsilon
!    (  1/2 x_5      1.2 x_4          x_3  )
!
! In particular the symmetry property of the stress tensor is exploited and only the lower triangle
! of the stress matrix is calculated.
!
!  The derivation of step (2) is made numerically by central finite displacements (2 displacements per component)
!  It is otained component wise.
!
!          (*)    First the geometry is distorted according to the given component. 
!          (**)   Run a full SCF cycle
!          (***)  Corresponding free energy build up the stress tensor
!
!  e.g. the 1,1 with grading 0 component reads  ('partial_strain_matrix')
!
!    (  1+delta     0       0  )    
!    (  1           1       0  ) *x  
!    (  1           0       1  )    
! 
!  where x represents coordinates and lattice (grading -1 would correspond to 1-delta)
!
! REMARKS:
!
! Delta is dimensionless!!
! => never ever convert it to "other" units
!-----------------------------------------------------------------------------------------------------------------

!  SOURCE

   contains


  !******
  !------------------------------------------------------------------------------
  !****s*  apply_strain_transformation_on_vector 
  !  NAME
  !   apply_strain_transformation_on_vector 
  !  SYNOPSIS

  subroutine apply_strain_transformation_on_vector  (column, row, grading, v_in, v_out)

  !  PURPOSE
  !    Apply partial strain transformation on a given vector 
  !
  !  USES
    implicit none
  
  !  ARGUMENTS
  integer,               intent(in)       ::  column, row, grading     ! grading should be 0 or 1
  real*8, dimension (3), intent(in)       ::  v_in
  real*8, dimension (3), intent(out)      ::  v_out
  real*8, dimension (3,3)                 ::  unity, partial_strain_matrix
  integer                                 ::  i,j

  !  INPUTS
  !    o column    -- component of the partial strain matrix
  !    o row       -- component of the partial strain matrix 
  !    o grading   -- grading for the partial strain matrix
  !    o v_in      -- vector on which partial strain matrix acts
  !  OUTPUTS
  !    o v_out     -- result after applying a partial strain on a given vector
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE


! fill unity
  do i=1,3,1
    do j=1,3,1
    unity(i,j)=0.0
    end do
    unity(i,i)=1.0
  end do
  
  partial_strain_matrix = 0.0
 
! Construct strain transformation for a given set of components
! Also convert the delta to bohr on the fly 
  do i=1,3,1
    do j=1,3,1

      ! Element matches to component of strain transformation
      if ((i==column).and.(j==row)) then  

         if (column==row) then  !diagonal elements
         partial_strain_matrix(i,j) = unity(i,j)+(-1d0)**grading *( delta_numerical_stress )
         
         else ! offdiagonal elements weighted by a factor of 0.5
         ! FK: We remove the factor 0.5, otherwise the off-diagonal elements are incorrect, see comment above.
         partial_strain_matrix(i,j) = unity(i,j)+(-1d0)**grading * (delta_numerical_stress )
         !partial_strain_matrix(i,j) = unity(i,j)+(-1d0)**grading* 0.5 * (delta_numerical_stress )

         end if

      ! Other elements are given by the unity  
      else 
         partial_strain_matrix(i,j) = unity(i,j)

      end if
    end do
  end do

 ! Apply the transformation on input vector and store in into output vector
 ! Blas: matrix on vector 

  call dgemv ('n', 3,3, 1d0, partial_strain_matrix, 3, v_in, 1, 0d0, v_out, 1)

 ! test begin
 !  write(use_unit,*) "grading", grading
 !  
 !  write(use_unit,*) "----- partial strain matrix ------------"  
 !  do i=1,3,1   ! print row-by-row
 !    write(use_unit,*) (partial_strain_matrix(i,j), j=1,3,1)
 !  end do
 !  write(use_unit,*) "----- end partial strain matrix ------------"  
 ! test debug

  end subroutine apply_strain_transformation_on_vector


  !******
  !------------------------------------------------------------------------------
  !****s*  remember_original_geometry 
  !  NAME
  !        remember_original_geometry  
  !  SYNOPSIS

  subroutine remember_original_geometry ()

  !  PURPOSE
  !    Stores current coordinate lattice geometry and volume
  !    
  !  USES
  !  ARGUMENTS
  implicit none 
  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
  
  original_coords                          = 0d0
  original_lattice_vector                  = 0d0 
  original_periodic_unit_cell_translations = 0d0 
  ! original_volume could be removed, since its re-calculated in initialize_bc_dependent_lists
  original_volume                          = 0d0

  original_coords          = coords
  original_lattice_vector  = lattice_vector      
  original_volume          = cell_volume 
  original_periodic_unit_cell_translations = periodic_unit_cell_translations 

  end subroutine remember_original_geometry


  !******
  !------------------------------------------------------------------------------
  !****s*  restore_zero_point_geometry  
  !  NAME
  !        restore_zero_point_geometry   
  !  SYNOPSIS

  subroutine restore_zero_point_geometry ()

  !  PURPOSE
  !    Restores original coordinate and lattice geometry  
  !    Calculates cell volume and reciprocial lattice vectors again
  !    
  !  USES
  !  ARGUMENTS
  implicit none

  integer :: i_atom
  character*80 :: info_str
 
  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  coords          = original_coords         
  lattice_vector  = original_lattice_vector 
  periodic_unit_cell_translations = original_periodic_unit_cell_translations 

  ! Recalculate conjugated quantities
  call initialize_bravais_quantities()

  end subroutine restore_zero_point_geometry


  !******
  !------------------------------------------------------------------------------
  !****s*  deallocate_numerical_stress  
  !  NAME
  !        deallocate_numerical_stress   
  !  SYNOPSIS

  subroutine deallocate_numerical_stress ()

  !  PURPOSE
  !    Deallocates arrays for calculating the numerical stress  
  !    Only those which store the zero point geometry:
  !    - original_coord 
  !    - original_lattice_vector
  !
  !  USES
  !  ARGUMENTS
  implicit none
  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
  
  if (allocated (original_coords         ))                    deallocate (original_coords)
  if (allocated (original_lattice_vector ))                    deallocate (original_lattice_vector)
  if (allocated (original_periodic_unit_cell_translations ))   deallocate (original_periodic_unit_cell_translations)

  end subroutine deallocate_numerical_stress


  !******
  !------------------------------------------------------------------------------
  !****s*  allocate_numerical_stress  
  !  NAME
  !        allocate_numerical_stress   
  !  SYNOPSIS
  
  subroutine allocate_numerical_stress ()
  
  !  PURPOSE
  !    Allocates arrays for calculating the numerical stress  
  !    Only those which store the zero point geometry:
  !    - original_coord 
  !    - original_lattice_vector
  !  USES
  !  ARGUMENTS
  implicit none
  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  if (.not. allocated (original_coords         ))                 allocate (original_coords(3,n_atoms))
  if (.not. allocated (original_lattice_vector))                  allocate (original_lattice_vector(3,n_periodic))
  if (.not. allocated (original_periodic_unit_cell_translations)) allocate (original_periodic_unit_cell_translations(3,n_atoms))

  end subroutine allocate_numerical_stress


  !******
  !------------------------------------------------------------------------------
  !****s*  get_indices_numerical_stress  
  !  NAME
  !        get_indices_numerical_stress    
  !  SYNOPSIS

  subroutine get_indices_numerical_stress (in_number, out_column, out_row, out_grading)

  !  PURPOSE
  !    Asigns to every valid counter for the stress {0 .. 11} the triplet {column, row, grading}
  !
  !  USES
  !  ARGUMENTS

  implicit none

  integer, intent (in)  :: in_number
  integer, intent (out) :: out_column, out_row, out_grading

  !  INPUTS
  !    o Counter for the stress tensor
  !        - in_number
  !  OUTPUTS
  !    o Index picture for the partial strain matrix
  !        - column
  !        - row
  !        - grading
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
  
  ! Asigning the indices by hand -- it's ugly!
  select case (in_number)

      case (0) 
          out_grading = 0
          out_column  = 1
          out_row     = 1

      case (1) 
          out_grading = 1
          out_column  = 1
          out_row     = 1

      case (2)
          out_grading = 0
          out_column  = 2
          out_row     = 1

      case (3) 
          out_grading = 1
          out_column  = 2
          out_row     = 1

      case (4) 
          out_grading = 0
          out_column  = 2
          out_row     = 2

      case (5)
          out_grading = 1
          out_column  = 2
          out_row     = 2

      case (6) 
          out_grading = 0
          out_column  = 3
          out_row     = 1

      case (7) 
          out_grading = 1
          out_column  = 3
          out_row     = 1

      case (8) 
          out_grading = 0
          out_column  = 3
          out_row     = 2

      case (9) 
          out_grading = 1
          out_column  = 3
          out_row     = 2

      case (10) 
          out_grading = 0
          out_column  = 3
          out_row     = 3

      case (11) 
          out_grading = 1
          out_column  = 3
          out_row     = 3

     case default 
         call localorb_info & 
("* Warning: something is seriously wrong with the counter for scf calculations for the numerical stress tensor",use_unit,'(2X,A)' )
         call localorb_info("* Don't trust the calculated stress tensor!",use_unit,'(2X,A)' )

  end select

  end subroutine get_indices_numerical_stress

  
  !******
  !------------------------------------------------------------------------------
  !****s*  apply_strain_transformation  
  !  NAME
  !        apply_strain_transformation   
  !  SYNOPSIS

  subroutine apply_strain_transformation (in_number)

  !  PURPOSE
  !  Construct for a given counter_number {0..11} all the corresponding distorted geometries
  !  (coordinates & lattice vectors) which come from application of the partial_strain_transformation.
  !  Also adapt quantities derived from these:
  !       - cell_volume
  !       - reciprocial_lattice_vector
  !  USES

  use species_data, only: species_name ! only for  debugging (print the stuff)

  !  ARGUMENTS

  implicit none

  integer, intent(in) :: in_number
  integer             :: column, row, grading
  integer             :: i_atom, i_coord, i_periodic
  character*120       :: info_str

  real*8, dimension(:,:), allocatable  :: coords_out
  real*8, dimension(:,:), allocatable  :: lattice_vector_out

  !  INPUTS
  !     o Counter for the stress tensor
  !
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  ! Allocate arrays in order to call 'apply_strain_transformation_on_vector'
  if (.not. allocated (coords_out         ))           allocate (coords_out(3,n_atoms))
  if (.not. allocated (lattice_vector_out ))           allocate (lattice_vector_out(3,n_periodic))


  ! Get the correct index picture
  call get_indices_numerical_stress (in_number, column, row, grading) 

  ! Apply strain transformation on lattice vectors
  do  i_periodic = 1, n_periodic, 1
  call apply_strain_transformation_on_vector &
       (column, row, grading, original_lattice_vector(:,i_periodic),lattice_vector_out(:,i_periodic))
  end do

  lattice_vector =  lattice_vector_out
  
  ! Apply strain transformation on the coords
  do  i_atom = 1, n_atoms, 1
  call apply_strain_transformation_on_vector  (column, row, grading, original_coords(:,i_atom),coords_out(:,i_atom) )
  end do
  
  coords=coords_out
  
  ! Update conjugated quantities
  call initialize_bravais_quantities()
   
  !------------------------
  ! Debug begin
  !------------------------

        ! print lattice
        ! do i_periodic = 1, n_periodic, 1
        !    write(info_str,'(2X,A1,A,3(2X,F15.6))') "|", "lattice_vector", &
        !         (lattice_vector(i_coord,i_periodic)*bohr ,i_coord=1,3,1)
        !    call localorb_info(info_str,use_unit,'(A)')
        ! enddo


        ! reciprocial print lattice
        ! do i_periodic = 1, n_periodic, 1
        !    write(info_str,'(2X,A1,A,3(2X,F10.6))') "|", "reciprocial lattice_vector ", &
        !         (recip_lattice_vector(i_coord,i_periodic)*bohr,i_coord=1,3,1)
        !    call localorb_info(info_str,use_unit,'(A)')
        ! enddo


        !  do i_atom = 1, n_atoms
        !     write(info_str,'(2X,A1,I5,A,A2,3(2X,F15.6))')    &
        !        & "|",i_atom, ": Species  ", species_name(species(i_atom)), (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
        !     call localorb_info(info_str,use_unit,'(A)')
        !  end do

  !------------------------
  ! Debug end 
  !------------------------

  if (allocated (coords_out         ))       deallocate (coords_out)
  if (allocated (lattice_vector_out ))       deallocate (lattice_vector_out)

  end subroutine apply_strain_transformation


  !******
  !------------------------------------------------------------------------------
  !****s* energy_contribution_to_numerical_stress   
  !  NAME
  !       energy_contribution_to_numerical_stress    
  !  SYNOPSIS

  subroutine  energy_contribution_to_numerical_stress (free_energy, counter)

  !  PURPOSE
  !   Writes for a given geometry (that corresponds to a counter
  !   for the numerical stress) the free energy into the array 
  !   'numerical_stress_tensor_components'
  !
  !  USES
  !  ARGUMENTS

  implicit none
  integer, intent(in)    :: counter
  real*8, intent(in)     :: free_energy
  integer                :: column, row, grading

  !  INPUTS
  !    o free_energy
  !    o Counter for the stress tensor 
  !   
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  ! Get the right index picture 
  call get_indices_numerical_stress (counter, column, row, grading)
  
  ! Write energy into the stress_tensor_components
  numerical_stress_tensor_components (column, row, grading) = free_energy  

  end subroutine  energy_contribution_to_numerical_stress


  !******
  !------------------------------------------------------------------------------
  !****s*  update_numerical_stress_tensor  
  !  NAME  update_numerical_stress_tensor
  !           
  !  SYNOPSIS

  subroutine  update_numerical_stress_tensor ()

  !  PURPOSE
  !      Calculates the numerical stress tensor and also fills the
  !      upper triangle of the stress tensor
  !
  !  USES
  !  ARGUMENTS

  implicit none
  integer :: grading
  integer :: column, row  ! counters

  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  ! The stress tensor is possibly not zero because
  ! of previous relaxation steps:
  numerical_stress_tensor = 0
 
  ! Here we only calculate the difference in the free energies
  ! The correct sign is established by the grading 
  
  do grading = 0,1,1
    numerical_stress_tensor (:,:) =     &
       numerical_stress_tensor (:,:) +  &
       (-1d0)** grading * numerical_stress_tensor_components (:,:, grading)
  end do

  ! Get prefactors and apply them to stress tensor 1/(2*delta*V)
  ! Factor of 1/(2*delta) comes from the definition of the numerical derivative 
  ! Need to take the zero point volume, which might not correspond to the current geomtry
  ! and is thus stored sperately at the beginning of the numerical stress calculation.
  ! The unit of the stress tensor is hartree/bohr^3

  numerical_stress_tensor= numerical_stress_tensor / (2*delta_numerical_stress*original_volume)

  ! Since the stress is only calculated for the lower triangle (including the diagonal) in the 
  ! stress matrix, obtain the missing entries by symmetrization (property of the stress tensor)
  ! By symmetrization we mean that we mirror the entries of the lower triangle to the upper one

  ! Symmetrization -- consider only off-diagonal elements:

  do column =2,3,1
     do row =1, column-1,1
        numerical_stress_tensor(row,column) = numerical_stress_tensor(column,row) 
     end do
  end do

  end subroutine update_numerical_stress_tensor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine could/should be removed 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !******
!  !---------------------------------------------------------------------------------------------------------------------------------
!  !****s* calculate_matrix_and_map_to_center_of_cell 
!  !  NAME
!  !       calculate_matrix_and_map_to_center_of_cell 
!  !  SYNOPSIS
!
!
!  subroutine   calculate_matrix_and_map_to_center_of_cell ()
!
!  !  PURPOSE
!  !  Calculate matrix for mapping coordinates to center of cell
!  !  and apply that transformation
!
!    implicit none
! 
!  integer ::  i_latt, i_atom
!  real*8, dimension(3):: work
!  integer,dimension(3):: ipivot
!  integer:: info, info_str
!
!  !  INPUTS
!  !    none
!  !  OUTPUT
!  !    none
!  !  AUTHOR
!  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  !  HISTORY
!  !    Release version, FHI-aims (2008).
!  !  SOURCE
!    ! These are needed in mapping to center cell in
!    ! map_to_center_cell subroutine
!
!    do i_latt = 1, n_periodic
!       map_to_center_cell_matrix(:,i_latt) = lattice_vector(:,i_latt)
!    end do
!
!    ! Get inverse matrix
!    call DGETRF(3, 3, map_to_center_cell_matrix, 3, ipivot, info )
!    if(info /= 0) then
!       write(use_unit,*) 'ERROR: lattice_vector is singular!'
!       stop
!    endif
!    call DGETRI(3, map_to_center_cell_matrix, 3, ipivot, work, 3, info)
!    if(info /= 0) then
!       write(use_unit,*) 'ERROR: lattice_vector is singular!'
!       stop
!    endif
!
!!    write(info_str,'(2X,A)') 'Mapping all atomic coordinates to central unit cell.'
!!    call localorb_info(info_str)
!    periodic_unit_cell_translations(:,:) = coords(:,:) + periodic_unit_cell_translations(:,:)
!
!    !CC: Take care in case of TDI:
!    if(use_thermodynamic_integration) then
!      do i_atom = 1,n_atoms
!         call TDI_map_to_center_cell(coords(1:3,i_atom),i_atom)
!      end do
!    else
!      do i_atom = 1,n_atoms
!         call map_to_center_cell(coords(1:3,i_atom))
!      end do
!    end if
!    periodic_unit_cell_translations(:,:) = periodic_unit_cell_translations(:,:) - coords(:,:)
!
!end   subroutine   calculate_matrix_and_map_to_center_of_cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !******
  !---------------------------------------------------------------------------------------------------------------------------------
  !**** NAME check_pressure 
  !          check_pressure 
  !  SYNOPSIS

  subroutine  check_pressure (stress_tensor) 

  !  PURPOSE
  !  Since Pressure is a hydrostatic quantity it only makes sense
  !  if the stress tensor is a multiple of the unit matrix
  !  This routine checks the stress tensor:
  !     - if the diagonal elements are the same up to a given threshold
  !     - the diagonal elements are zero up to the same threshold
  !  If this is not the case then there is a warning.

    implicit none
  real*8, intent(in)  :: stress_tensor(3,3)

  ! threshold for which the system is asumed to be isotropic
  real*8   :: hydrostatic_threshold = 0.001 ! in ev/A^3
  real*8   :: diagonal_stress_entries(3) 
  logical  :: diagonal_ok           = .true.
  logical  :: off_diagonal_ok       = .true.

  ! counters
  integer  :: i, j

  ! convert threshold:
  hydrostatic_threshold =  hydrostatic_threshold * (bohr**3)/hartree
  
  ! check off-diagonals
  do i=1,3,1
    do j=1,i-1,1
      if (abs (stress_tensor(i,j)) .gt. hydrostatic_threshold) off_diagonal_ok = .false.
    end do
  end do 

  ! check diagonals
  do i=1,3,1
    diagonal_stress_entries(i)= stress_tensor(i,i)
  end do
  if (abs (maxval(diagonal_stress_entries) - minval (diagonal_stress_entries)) .gt.  hydrostatic_threshold ) then
     diagonal_ok =.false.
  end if 
 
  ! Write warnings if necessary
  if ( (.not. diagonal_ok) .or. (.not. off_diagonal_ok)) then 
    call localorb_info('',use_unit)
    call localorb_info("* Warning: Stress tensor is anisotropic. Be aware that pressure is an isotropic quantity.",use_unit,'(1X,A)' )
  end if  

  end subroutine check_pressure

end module numerical_stress 
