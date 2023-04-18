!****h* FHI-aims/geometry
!  NAME
!    geometry
!  SYNOPSIS

module geometry

  !  PURPOSE
  !  Module geometry handles basic tasks related to the atomic geometry.
  !
  !  Subroutines:
  !  o allocate_geometry
  !  o cleanup_geometry
  !  o output_structure
  !  o output_bsse_stucture
  !  o cart2frac
  !  o frac2cart
  !  o initialize_bravais_quantities
  !  o all_atoms_to_cart
  !  o min_atomic_dist_old (retired)
  !  o min_atomic_dist
  !  o min_grid_centre_dist
  !  o min_multipole_dist 
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




  implicit none

  !     global variable declarations - exported to other program parts

  !     coords   : all atomic coordinates
  !     species  : species type for all atoms

  real*8, dimension(:,:), allocatable:: coords           ! coords(3,n_atoms)
  real*8, dimension(:,:), allocatable:: frac_coords      ! frac_coords(3,n_atoms)
  real*8, dimension(:,:), allocatable:: input_coords     ! input_coords(3,n_atoms) VA: for cell relaxation + fixed cartesian coords
  real*8, dimension(:,:), allocatable:: occ_coords       ! coords(3,n_occ_atoms)
  real*8, dimension(:,:), allocatable:: empty_coords     ! coords(3,n_empty_atoms)


  integer, dimension(:), allocatable:: species           ! species(n_atoms)
  integer, dimension(:), allocatable:: coord_basis_atoms ! coord_basis_atoms(n_atoms)

  logical, dimension(:),   allocatable :: empty          ! empty(n_atoms)

  real*8, dimension(:,:), allocatable:: lattice_vector
  real*8, dimension(:,:), allocatable:: recip_lattice_vector
  real*8, dimension(:),   allocatable:: length_of_lattice_vector
  real*8                             :: map_to_center_cell_matrix(3,3)
  real*8,dimension(:,:), allocatable :: periodic_unit_cell_translations ! internal to-center shift of every atom 
  real*8                             :: cell_volume
  ! Logical variable that memorizes if the unit cell shape changed after
  ! a geometry update.
  logical                            :: unit_cell_shape_changed

  !==================================
  !   For constant force on atoms
  !   e.g. model nano indentation
  !==================================
  real*8, dimension(:,:), allocatable:: external_forces     ! forces acting on atom (3,n_atoms)
  logical                            :: external_forces_on
  
  !==================================
  !      For RRS-PBC scheme
  !==================================
  ! rrs_pbc_lattice_vector           :: the lattice vectors for RRS-PBC projection
  ! rrs_pbc_inv_lattice_vector       :: the inverse matrix of lattice vectors
  ! rrs_pbc_recip_lattice_vector     :: the associated k vectors
  ! rrs_pbc_length_of_lattice_vector
  !
  real*8, dimension(3,3)             :: rrs_pbc_lattice_vector
  real*8, dimension(3,3)             :: rrs_pbc_inv_lattice_vector
  real*8, dimension(3,3)             :: rrs_pbc_recip_lattice_vector
  real*8, dimension(3)               :: rrs_pbc_length_of_lattice_vector
  real*8                             :: rrs_pbc_cell_volume
  real*8  ,dimension(:),   allocatable   :: initial_moment
  real*8,  dimension(:),   allocatable   :: initial_charge
  integer, dimension(:),   allocatable   :: kind_of_initial
  real*8 :: total_initial_charge

  real*8,  dimension(:), allocatable :: type_moment
  real*8,  dimension(:), allocatable :: type_charge
  integer, dimension(:), allocatable :: atom_type
  integer, dimension(:), allocatable :: type_species
  integer, dimension(:), allocatable :: type_kind

  real*8,dimension(:),  allocatable  :: homogeneous_field ! homogeneous_field(3)
  real*8 :: image_plane_position, image_plane_cutoff !Image potential variables
  real*8 ::  image_plane_scaling ! Image potential variables


  real*8, dimension(:,:),allocatable  :: multipole_coords  ! species(3,n_multipoles)
  real*8, dimension(:,:),allocatable  :: multipole_data    ! data(x,n_multipoles)
  real*8, dimension(:),  allocatable  :: multipole_charge  ! charge(n_multipoles)
  integer,dimension(:),  allocatable  :: multipole_order   ! order(n_multipoles)

  logical:: geometry_initialized

  !      logical:: real_eigenvalues

  ! Each hess_blocks(:,:, i_block) contains the block corresponding to
  ! (/i_atom, j_atom/) == hess_block_pos(:,i_block).  Do not forget the
  ! transposed.
  real*8,  dimension(:,:,:), allocatable :: hess_blocks            ! (3,3,n_hess_blocks)
  integer, dimension(:,:),   allocatable :: hess_block_pos         ! (2, n_hess_blocks)
  real*8,  dimension(:,:,:), allocatable :: hess_blocks_lv         ! (3,3,n_hess_blocks_lv)
  integer, dimension(:,:),   allocatable :: hess_block_pos_lv      ! (2, n_hess_blocks_lv)
  real*8,  dimension(:,:,:), allocatable :: hess_blocks_lv_atom    ! (3,3,n_hess_blocks_lv_atom)
  integer, dimension(:,:),   allocatable :: hess_block_pos_lv_atom ! (2, n_hess_blocks_lv_atom)

  !pseudopotentials

  integer, dimension(:), allocatable :: pp_species           ! pp_species(n_pp_atoms)
  real*8, dimension(:,:),allocatable :: pp_coords            ! pp_coords(3,n_pp_atoms)

  ! Whether the given atom is included in the magnetic response
  ! calculations
  logical, allocatable :: mr_atoms(:)
  ! If specified, these values of the magnetic moment and the nuclear
  ! spin for a given atom are taken instead of the default values from
  ! MR_nuclear_data.f90.
  real(8), allocatable :: mr_custom_mus(:), mr_custom_spins(:)
  ! If specified, the nuclear data corresponding to this isotope is
  ! used for the magnetic response calculations.
  integer, allocatable :: isotope_numbers(:)

contains
  !******
  !---------------------------------------------------------------------
  !****s* geometry/allocate_geometry
  !  NAME
  !   allocate_geometry
  !  SYNOPSIS

  subroutine set_geometry_defaults (  )
  implicit none

  unit_cell_shape_changed = .false.
  external_forces_on = .false.
  total_initial_charge = 0.0d0
  image_plane_scaling = 0.0d0
  geometry_initialized = .false.

  end subroutine set_geometry_defaults

  subroutine allocate_geometry (  )

    !  PURPOSE
    !  Subroutine allocate_grids allocates the necessary memory for all real-space grids. 
    !  This allocation is done once at the beginning of the code, since the grids
    !  are needed everywhere.
    !
    !  In principle one could allocate everything in a much more fine-grained way - but only 
    !  if there is a clear need. 
    !
    !    AJL: If you add a new allocation *PLEASE PLEASE* include a
    !    deallocation in the cleanup routine, otherwise problems are
    !    encountered when treating aims as a library
    !
    !  USES

    use dimensions, only : n_atoms, n_empty_atoms, n_hess_blocks, n_hess_blocks_lv, &
                           n_hess_blocks_lv_atom, n_multipoles, n_occ_atoms, n_periodic, &
                           n_pp_atoms, n_spin, relax_unit_cell, use_embedding_potential, &
                           use_embedding_pp 
    use mpi_tasks,  only : check_allocation
    implicit none

    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    integer:: info


    allocate ( coords(3,n_atoms) ,stat=info)
    call check_allocation(info, 'coords                        ')

    allocate ( frac_coords(3,n_atoms) ,stat=info)
    call check_allocation(info, 'frac_coords                   ')
    frac_coords = 0d0

    allocate ( species(n_atoms) ,stat=info)
    call check_allocation(info, 'species                       ')

    allocate ( coord_basis_atoms(n_atoms) ,stat=info)
    call check_allocation(info, 'coord_basis_atoms             ')
    coord_basis_atoms = -1

    allocate ( empty(n_atoms) ,stat=info)
    call check_allocation(info, 'gempty                        ')
    empty = .false.
    
    allocate ( external_forces(3,n_atoms) ,stat=info)
    call check_allocation(info, 'external_forces                  ')
    external_forces = 0d0

    if (n_periodic.gt.0) then
       allocate ( lattice_vector(3,n_periodic) ,stat=info)
       call check_allocation(info, 'lattice_vector                ')

       allocate ( length_of_lattice_vector(n_periodic) ,stat=info)
       call check_allocation(info, 'length_of_lattice_vector      ')

       allocate ( recip_lattice_vector(3,n_periodic),stat=info )
       call check_allocation(info, 'recip_lattice_vector          ')

       allocate ( periodic_unit_cell_translations(3,n_atoms),stat=info)
       call check_allocation(info, 'periodic_unit_cell_translation')
       periodic_unit_cell_translations(:,:) = 0d0
    
       if (relax_unit_cell > 0) then
          allocate ( input_coords(3,n_atoms) ,stat=info)
          call check_allocation(info, 'frac_coords                   ')
       end if
     
    else 
       ! dummy lattice_vector for cluster case.
       ! All 3x3 components are needed since they are used in some contexts.
       ! TODO: Check if this is also needed for n_periodic = 1 or 2
       allocate ( lattice_vector(3,3) ,stat=info)
       call check_allocation(info, 'lattice_vector                ')
       ! for safety reasons asign  ridiculous value
       lattice_vector = 0
       lattice_vector(1,1) = 1e9
       lattice_vector(2,2) = 1e9
       lattice_vector(3,3) = 1e9

    end if

    allocate ( initial_charge(n_atoms),stat=info)
    call check_allocation(info, 'initial_charge                ')

    allocate (type_charge(n_atoms),stat=info)
    call check_allocation(info, 'type_charge                   ')

    if (n_spin.gt.1) then
       allocate ( initial_moment(n_atoms),stat=info)
       call check_allocation(info, 'initial_moment                ')

       allocate (kind_of_initial(n_atoms),stat=info)
       call check_allocation(info, 'kind_of_initial               ')

       allocate (type_moment(n_atoms),stat=info)
       call check_allocation(info, 'type_moment                   ')

       allocate (type_kind(n_atoms),stat=info)
       call check_allocation(info, 'type_kind                     ')

    end if
    allocate (atom_type(n_atoms),stat=info)
    call check_allocation(info, 'atom_type                     ')

    allocate (type_species(n_atoms),stat=info)
    call check_allocation(info, 'type_species                  ')

    !       Nadia
    if (use_embedding_potential)  then
       allocate ( multipole_coords(3,n_multipoles) ,stat=info)
       call check_allocation(info, 'use_embedding_potential       ')

       allocate ( multipole_order(n_multipoles) ,stat=info)
       call check_allocation(info, 'multipole_order               ')

       allocate ( multipole_charge(n_multipoles) ,stat=info)
       call check_allocation(info, 'multipole_charge              ')

       allocate ( homogeneous_field(3) ,stat=info) 
       call check_allocation(info, 'homogeneous_field             ')
       homogeneous_field = 0.d0

       !       FIXME: 3 is hard-coded; fix for higher order moments
       allocate ( multipole_data(3,n_multipoles) ,stat=info) 
       call check_allocation(info, 'multipole_data                ')

    end if

    if (use_embedding_pp) then
      allocate(pp_species(n_pp_atoms) ,stat=info)
       call check_allocation(info, 'pp_species       ')

      allocate(pp_coords(3,n_pp_atoms) ,stat=info)
       call check_allocation(info, 'pp_coords       ')

    end if


    !       BSSE    
    allocate ( occ_coords (3, n_occ_atoms),stat=info)
    call check_allocation(info, 'occ_coords                    ')

    allocate ( empty_coords (3, n_empty_atoms),stat=info)
    call check_allocation(info, 'empty_coords                  ')

    ! Magnetic response
    allocate(mr_atoms(n_atoms), mr_custom_mus(n_atoms), &
         & mr_custom_spins(n_atoms), isotope_numbers(n_atoms), stat=info)
    call check_allocation(info, 'MR work variable', &
         & 'geometry::allocate_geometry')
    mr_atoms = .false.
    mr_custom_mus = 1d300 ! Defaults are in MR_nuclear_data.f90
    mr_custom_spins = -1d0
    isotope_numbers = -1d0
    
  end subroutine allocate_geometry
  !******
  !---------------------------------------------------------------------
  !****s* geometry/output_structure
  !  NAME
  !   output_structure
  !  SYNOPSIS

  subroutine output_structure

    !  PURPOSE
    !   Prints out the coordinates of the atoms. This is used in geometry relaxations.
    !
    !   Note that this routine also does some work that reaches beyond a mere output.
    !   If the unit cell changed along the way, we also _update_ the stored internal
    !   shifts of coordinates by multiples of the lattice vectors, stored in the array:
    !
    !   periodic_unit_cell_translations
    !
    !   On entry, these reflect the internal shift of each atom (compared to the user input geometry)
    !   in units of the OLD lattice vectors (fractional coordinates, prior to relaxation). 
    !
    !   Upon exit, they reflect the internal shift of each atom in units of the new lattice vectors
    !   (after relaxation).
    !
    !  USES

    use constants,       only : bohr
    use dimensions,      only : n_atoms, n_multipoles, n_periodic, n_pp_atoms, n_real_atoms
    use runtime_choices, only : add_embedding_grids, output_in_original_unit_cell, output_level, relax_geo
    use species_data,    only : species_pseudoized, species_name
    use localorb_io,     only : localorb_info, use_unit
    implicit none

    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    integer                      :: i_atom, i_coord, i_lv, i_multipole, i_pp_atom
    real*8, dimension(3,n_atoms) :: coords_temp
    real*8, dimension(3,n_atoms) :: frac_coords_temp
    real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_cart
    character*100                :: info_str

    ! use temporary coordinates (i.e. NOT the internal ones!) to output structure ... 
    coords_temp(:,:) = coords(:,:)

    !... because if there was a request to keep the original unit cell configuration, that should not mess up other things...
    if ((output_in_original_unit_cell).and.(n_periodic>0)) then

       ! This routine touches only locally stored information.
       call frac2cart (lattice_vector,                    &
                       periodic_unit_cell_translations,   &
                       periodic_unit_cell_translations_cart)

       coords_temp(:,:) = coords(:,:)+periodic_unit_cell_translations_cart(:,:)

       ! This statement here is not just there for output purposes. It actually modified the
       ! stored vectors "periodic_unit_cell_translations" to reflect fractional coordinates
       ! in units of the _new_ lattice vectors
       call cart2frac (lattice_vector,                       &
                       periodic_unit_cell_translations_cart, &
                       periodic_unit_cell_translations)
    end if

    write(info_str,'(25X,A,13X,A,13X,A,7X)') &
         "x [A]","y [A]","z [A]"
    call localorb_info( info_str )

    if (n_periodic > 0) then 
       do i_lv = 1, n_periodic, 1
          write(info_str,'(2X,A,3(2X,F16.8),2X)') &
               "lattice_vector ", &
               (lattice_vector(i_coord,i_lv)*bohr, i_coord=1,3,1)
          call localorb_info( info_str )
       enddo
          call localorb_info('',use_unit)
    end if

    do i_atom = 1, n_atoms, 1
       if(species_pseudoized(species(i_atom))) then
 
      write(info_str,'(7X,A,3(2X,F16.8),2X,A)') &
            "pseudocore", &
            (coords_temp(i_coord,i_atom)*bohr, i_coord=1,3,1), &
            species_name(species(i_atom))

       else if(empty(i_atom)) then 

      write(info_str,'(12X,A,3(2X,F16.8),2X,A)') &
            "empty", &
            (coords_temp(i_coord,i_atom)*bohr, i_coord=1,3,1), &
            species_name(species(i_atom))

       else
      write(info_str,'(12X,A,3(2X,F16.8),2X,A)') &
            "atom ", &
            (coords_temp(i_coord,i_atom)*bohr, i_coord=1,3,1), &
            species_name(species(i_atom))


       endif
       call localorb_info( info_str )
    enddo

    if((n_pp_atoms.gt.0).and.(.not.add_embedding_grids)) then
      do i_pp_atom = 1, n_pp_atoms
         write(info_str,'(7X,A,3(2X,F16.8),2X,A)') &
              "pseudocore", &
              (pp_coords(i_coord,i_pp_atom)*bohr, i_coord=1,3,1), &
              trim(species_name(species(n_real_atoms + i_pp_atom)))
       call localorb_info( info_str )
      enddo
    endif


!DB multipole output -----
    if (n_multipoles.gt.0) then
      do i_multipole = 1, n_multipoles, 1
        write(info_str,'(7X,A,3(2X,F16.8),2X,I4,2X,F10.3)') &
            "multipole ", &
            (multipole_coords(i_coord,i_multipole)*bohr, i_coord=1,3,1), &
             multipole_order(i_multipole), multipole_charge(i_multipole)
        call localorb_info( info_str )
      enddo
    endif
!DB -pseudoatoms output-------
!    if (use_embedding_pp) then
!    coords_temp(:,:) = pp_coords(:,:)
!    do i_pp_atom = 1, n_pp_atoms, 1
!       write(info_str,'(12X,A,3(2X,F16.8),2X,A)') &
!            "pseudocore ", &
!            (coords_temp(i_coord,i_pp_atom)*bohr, i_coord=1,3,1), &
!            pp_species_name(pp_species(i_pp_atom))
!       call localorb_info( info_str )
!    enddo
!    coords_temp(:,:) = coords(:,:)
!    endif
!DB----

   ! Output fractional coordinates 
   if ((n_periodic>0) .and. ((relax_geo > 0) .or. (output_level /= 'MD_light' ))) then
      call cart2frac(lattice_vector, coords_temp, frac_coords_temp)
      call localorb_info('',use_unit)
      write(info_str,'(2X,A)') "Fractional coordinates: "
      call localorb_info( info_str )
      write(info_str,'(25X,A,16X,A,16X,A,7X)') &
            "L1","L2","L3"
      call localorb_info( info_str )
      do i_atom = 1, n_atoms, 1
         write(info_str,'(7X,A,3(2X,F16.8),2X,A)') &
              "atom_frac ", &
              (frac_coords_temp(i_coord,i_atom), i_coord=1,3,1), &
              species_name(species(i_atom))
         call localorb_info( info_str )
      enddo
   endif



  end subroutine output_structure
  !******
  !---------------------------------------------------------------------
  !****s* geometry/cleanup_geometry
  !  NAME
  !   cleanup_geometry
  !  SYNOPSIS

  subroutine cleanup_geometry ()

    !  PURPOSE
    !    Subroutine cleanup_geometry deallocates everything to do with geometry
    !
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    implicit none



    if ( allocated (coords) ) then
       deallocate ( coords )
    end if
    if ( allocated (frac_coords) ) then
       deallocate ( frac_coords )
    end if
    if ( allocated (input_coords) ) then
       deallocate (input_coords )
    end if
    if ( allocated (species) ) then
       deallocate ( species )
    end if
    if ( allocated (coord_basis_atoms) ) then
       deallocate ( coord_basis_atoms )
    end if
    if ( allocated (empty) ) then
       deallocate ( empty )
    end if
    if ( allocated (occ_coords) ) then
       deallocate ( occ_coords )
    end if
    if ( allocated (empty_coords) ) then
       deallocate (empty_coords )
    end if
    if ( allocated (external_forces) ) then
       deallocate ( external_forces )
    end if
    if ( allocated (lattice_vector) ) then
       deallocate ( lattice_vector )
    end if
    if ( allocated (recip_lattice_vector) ) then
       deallocate ( recip_lattice_vector )
    end if
    if ( allocated (initial_moment) ) then
       deallocate ( initial_moment )
    end if
    if ( allocated (kind_of_initial) ) then
       deallocate ( kind_of_initial )
    end if
    if ( allocated (initial_charge) ) then
       deallocate ( initial_charge )
    end if
    if ( allocated (multipole_coords) ) then
       deallocate ( multipole_coords )
    end if
    if ( allocated (multipole_order) ) then
       deallocate ( multipole_order )
    end if
    if ( allocated (multipole_charge) ) then
       deallocate ( multipole_charge )
    end if

    if ( allocated (multipole_data) ) then
       deallocate ( multipole_data )
    end if
    if ( allocated (pp_coords) ) then
       deallocate ( pp_coords)
    end if
    if ( allocated (pp_species) ) then
       deallocate (pp_species)
    end if
    if ( allocated (homogeneous_field) ) then
       deallocate ( homogeneous_field )
    end if
    if (allocated(type_moment)) then
       deallocate(type_moment)
    end if
    if (allocated(type_kind)) then
       deallocate(type_kind)
    end if
    if (allocated(type_charge)) then
       deallocate(type_charge)
    end if
    if (allocated(atom_type)) then
       deallocate(atom_type)
    end if
    if (allocated(type_species)) then
       deallocate(type_species)
    end if

    if (allocated(periodic_unit_cell_translations)) then
       deallocate(periodic_unit_cell_translations)
    end if

    if (allocated(length_of_lattice_vector)) then
       deallocate(length_of_lattice_vector)
    end if

    ! Magnetic response
    if (allocated(mr_atoms))         deallocate(mr_atoms)
    if (allocated(mr_custom_mus))    deallocate(mr_custom_mus)
    if (allocated(mr_custom_spins))  deallocate(mr_custom_spins)
    if (allocated(isotope_numbers))  deallocate(isotope_numbers)


  end subroutine cleanup_geometry
  !******
  !------------------------------------------------------------------------------
    !******
  !---------------------------------------------------------------------
  subroutine output_bsse_structure (current_bsse_atom, n_electrons)

    !  PURPOSE
    !   Prints out the coordinates of the atoms(or ghosts) during calculations for
    !   atomization BSSE correction.
    !
    !  USES
    use constants,       only : bohr
    use dimensions,      only : calculate_all_eigenstates, n_atoms, n_empty_atoms, &
                                n_states, n_occ_atoms, n_species, use_full_spectrum, &
                                n_max_basis
    use runtime_choices, only : charge, n_empty
    use species_data,    only : species_name, species_z
    use mpi_tasks,       only : myid
    use localorb_io,     only : use_unit
    implicit none

    !  INPUTS
    !    current_bsse_atom
    !  OUTPUT
    !    n_electrons
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    logical eof
    integer :: i_code, i_species
    integer :: current_bsse_atom
    integer :: i_coord, i_atom,i_atom_2
    integer :: i_empty_atoms,i_empty_atoms_2
    real*8 :: n_electrons
    integer, dimension(n_atoms) :: species_temp_2
    real*8, dimension(3,n_atoms) :: coords_temp_2
    real*8 :: distance, BSSE_radius
    character*20 desc_str
    character*20 species_temp
     
!   initialization 
    n_electrons = 0.d0
    BSSE_radius=1000000.0
    i_empty_atoms = 0
    
! deallocate arrays used elsewhere before    
    
    if ( allocated (empty) ) then
       deallocate ( empty )
    end if
    if ( allocated (empty_coords) ) then
       deallocate (empty_coords )
    end if

! now reallocate them :
 
    allocate ( empty(n_atoms))  
    empty = .false.  
    allocate ( empty_coords (3, n_empty_atoms))

    
! read in the original geometry.in file to reset the coordinates, similar to read_geo
       
      eof = .false.

      i_atom = 0
        !  read in the original coords again.
        !  open input file

      open (8, FILE="geometry.in" )

   !  read first line

      read (8,*,iostat=i_code) desc_str

      do while (.not.eof)

        if (desc_str(1:1).eq."#") then

          continue

        else
        
           backspace(8)

          i_atom = i_atom+1

          read(8,*) desc_str, (coords(i_coord,i_atom), i_coord=1,3,1), &
            species_temp  

        !  transform coordinates into atomic units (bohr)

          do i_coord = 1,3,1
            coords(i_coord,i_atom) = coords(i_coord,i_atom)/bohr
          enddo

        !  read in species

          do i_species=1, n_species, 1
            if (species_temp.eq.species_name(i_species)) then  
              species(i_atom) = i_species  
            end if
          end do
            
        end if    

        read (8,*,iostat=i_code) desc_str

        if (i_code.ne.0) then
          eof = .true.
        end if

      end do

!  close geometry.in

      close(8)  
    
    ! use temporary coordinates
    coords_temp_2(:,:) = coords(:,:)
    species_temp_2(:)=species(:)
   

    do i_atom_2 = 1, n_atoms, 1
       if (i_atom_2==current_bsse_atom) then
          ! The occupied atom is now the first atom 
          ! mimics what read_geo does
           do i_coord = 1,3,1
              coords(i_coord,1)=&  
                    coords_temp_2(i_coord,i_atom_2)  
           end do
          ! also update species
           species(1)=species_temp_2(i_atom_2)
          ! count the total number of electrons=number of electron in the occupied atom  
           n_electrons = n_electrons + species_z(species(1))
       else 
           distance=0
           do i_coord=1,3,1
               distance=distance+(coords_temp_2(i_coord,current_bsse_atom) &
                  - coords_temp_2(i_coord,i_atom_2))**2
           enddo
           ! print the ghost atom co-ordinates
           if (distance < (BSSE_radius**2) ) then
           ! re_assign positions of the empty atoms: refer to read_geo
              i_empty_atoms = i_empty_atoms+1 
              i_empty_atoms_2 = i_empty_atoms+n_occ_atoms 
              do i_coord = 1,3,1
              empty_coords(i_coord,i_empty_atoms) = &   
                     coords_temp_2(i_coord,i_atom_2)  
              coords(i_coord,i_empty_atoms_2)= & 
                 empty_coords (i_coord,i_empty_atoms)  
              enddo
              species(i_empty_atoms_2)=species_temp_2(i_atom_2)
           endif           
           !  store the position of empty atoms in the logical array
           empty(i_empty_atoms+n_occ_atoms) = .true. 
       endif
    enddo

    ! print out the new coordinates
    if (myid.eq.0) then
      write(use_unit,'(2X,A)') "| Atomic structure: "
      write(use_unit,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Atom ","x [A]","y [A]","z [A]"
      do i_atom_2 = 1, n_occ_atoms, 1
         write(use_unit,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
           "|",i_atom_2, ": Species ", &
           species_name(species(i_atom_2)), &
           (coords(i_coord,i_atom_2)*bohr, i_coord=1,3,1)
     enddo

     if (n_occ_atoms.lt.n_atoms) then
         write(use_unit,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Ghost","x [A]","y [A]","z [A]"
         do i_atom_2 = n_occ_atoms+1, n_atoms, 1
          write(use_unit,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
              "|",i_atom_2, ": Species ", &
              species_name(species(i_atom_2)), &
              (coords(i_coord,i_atom_2)*bohr, i_coord=1,3,1)
         enddo
     endif
    endif

  !  determine minimum number of Kohn-Sham states to be computed during the scf cycle
    n_states = ceiling(real(n_electrons + charge)/2.) + ceiling(real(n_empty * n_atoms)/2.)
    if (use_full_spectrum.or.calculate_all_eigenstates) then
        n_states = n_max_basis
    endif

  end subroutine output_bsse_structure

  !******
  !------------------------------------------------------------------------------
  !****s* geometry/cart2frac
  !  NAME
  !  cart2frac 
  !  SYNOPSIS

  subroutine cart2frac (lv, cart_coords, frac_coords)

    !  PURPOSE
    !  cartesian --> fractional basis transformation 
    !  with respect to the current lattice vectors 
    !
    !  USES
    use dimensions, only : n_atoms, n_periodic
    use bravais,    only : get_map_to_center_cell_matrix
    implicit none

    !  ARGUMENTS
    !  INPUTS
    !     o lv          -- lattice vectors
    !     o cart_coords -- cartesian coords
    !  OUTPUTS
    !     o frac_coords -- fractional coords
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  Local arguments
    real*8 , dimension(3,n_atoms),   intent(in)  :: cart_coords
    real*8 , dimension(3,n_periodic),intent(in)  :: lv 
    real*8 , dimension(3,n_atoms),   intent(out) :: frac_coords
    real*8 , dimension(n_periodic,3)             :: inv_lv  ! inverse lattice vector matrix (== map_to_center_cell_matrix)

    if (n_periodic > 0) then

       call get_map_to_center_cell_matrix(n_periodic,   &
                                          lv,           &
                                          inv_lv)

       frac_coords = matmul(inv_lv, cart_coords)
    !   frac_coords = matmul(map_to_center_cell_matrix, cart_coords)
    end if 
    
  end subroutine cart2frac

  !******
  !----------------------------------------------------------------------------
  !****s* geometry/frac2cart
  !  NAME
  !  frac2cart
  !  SYNOPSIS

  subroutine frac2cart (lv, frac_coords, cart_coords)

    !  PURPOSE
    !  fractional --> cartesian basis transformation 
    !  with respect to the current lattice vectors 
    !
    !  USES
    use dimensions, only : n_atoms, n_periodic
    implicit none

    !  ARGUMENTS
    !  INPUTS
    !     o lv          -- lattice vectors
    !     o frac_coords -- fractional coords
    !  OUTPUTS
    !     o cart_coords -- cartesian coords
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    !  Local arguments
    real*8 , dimension(3,n_atoms),    intent(in)  :: frac_coords
    real*8 , dimension(3,n_periodic), intent(in)  :: lv 
    real*8 , dimension(3,n_atoms),    intent(out) :: cart_coords

    if (n_periodic > 0) then
      cart_coords=matmul(lv, frac_coords)   
    end if

  end subroutine frac2cart
  !******
  !----------------------------------------------------------------------------
  !****s* geometry/initialize_bravais_quantities
  !  NAME
  !    initialize_bravais_quantities
  !  SYNOPSIS

  subroutine initialize_bravais_quantities()

    !  PURPOSE
    !
    !    Initialize quantities which depend only on the Bravais lattice (like
    !    the reciprocal lattice, the cell volume, ...).
    !
    !    This function *must* be called every time lattice_vector is changed.
    !
    !  USES

    use dimensions, only : n_periodic 
    use bravais,    only : get_cell_volume, get_reciprocal_vectors, get_length_of_lv, &
                           get_map_to_center_cell_matrix
    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none [geometry.f90:lattice_vector]
    !  OUTPUTS
    !    none [geometry.f90:cell_volume,recip_lattice,lenght_of_...]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'initialize_bravais_quantities'

    if(n_periodic > 0)then
       call get_cell_volume(lattice_vector, cell_volume)
       call get_reciprocal_vectors(n_periodic, lattice_vector, &
       &                           recip_lattice_vector)
       call get_length_of_lv(n_periodic, lattice_vector, &
       &                     length_of_lattice_vector)
       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
       &                                  map_to_center_cell_matrix)
    end if

  end subroutine initialize_bravais_quantities
  !******
  !----------------------------------------------------------------------------
  !****s* geometry/output_bravais_quantities
  !  NAME
  !    output_bravais_quantities
  !  SYNOPSIS

  subroutine output_bravais_quantities()

    !  PURPOSE
    !
    !    Output quantities which are derived from the Bravais lattice
    !    (i.e. the reciprocal lattice and the unit cell volume).
    !
    !  USES
    use constants,   only : bohr
    use dimensions,  only : n_periodic 
    use localorb_io, only : use_unit, localorb_info
    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none [geometry.f90:lattice_vector]
    !  OUTPUTS
    !    none [geometry.f90:cell_volume,recip_lattice,lenght_of_...]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character*150 :: info_str
    character(*), parameter :: func = 'output_bravais_quantities'

    if(n_periodic > 0)then
       call localorb_info('  ', use_unit, '(A)')
       write(info_str,'(2X,A)') "Quantities derived from the lattice vectors:"
       call localorb_info(info_str,use_unit,'(A)')        
       write(info_str,'(2X,A,3F10.6)') &
       & '| Reciprocal lattice vector 1:', recip_lattice_vector(:,1)/bohr
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str, '(2X,A,3F10.6)') &
       & '| Reciprocal lattice vector 2:', recip_lattice_vector(:,2)/bohr
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,3F10.6)') &
       & '| Reciprocal lattice vector 3:', recip_lattice_vector(:,3)/bohr
       call localorb_info(info_str,use_unit,'(A)')

       write(info_str,'(2X,A,E15.6,2X,A)') & 
       & '| Unit cell volume                               :', &
       &                                   cell_volume * bohr**3, 'A^3'
       call localorb_info(info_str,use_unit,'(A)')
    end if

  end subroutine output_bravais_quantities
  !******
  !----------------------------------------------------------------------------
  !****s* geometry/all_atoms_to_cart
  !  NAME
  !    all_atoms_to_cart 
  !  SYNOPSIS

  subroutine all_atoms_to_cart()

    !  PURPOSE
    !
    !  USES
    use constants,   only : bohr
    use dimensions,  only : n_atoms
    use mpi_tasks,   only : aims_stop
    use localorb_io, only : use_unit, localorb_info
    implicit none
    
    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none [geometry.f90:]
    !  OUTPUTS
    !    none [geometry.f90:]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE
   
    real*8, dimension(3,n_atoms) :: coords_temp
    integer                      :: i_atom
    character*150 :: info_str

    coords_temp(:,:) = 0.

    ! 1. Transform all coordinates: frac -> cart
    call frac2cart (lattice_vector, coords, coords_temp)

    ! 2. For fractional coordinates, the components should 
    !    not be multiplied by 1/bohr (since lattice_vectors are already)
    !    as it is done in read_geo.f90
    !    --> revert this

    do i_atom = 1,n_atoms,1
       coords_temp(:,i_atom) = coords_temp(:,i_atom)*bohr
    end do

    ! 3. Set back those coordinates which were already in a cartesian basis
    do i_atom = 1, n_atoms , 1
       if (coord_basis_atoms(i_atom) == 0) then
           coords_temp(:,i_atom) = coords(:,i_atom)
       else if (coord_basis_atoms(i_atom) == 1) then
          cycle
       else
           write(info_str,'(1X,A)') &
           "* Coordinate basis for atoms unknown!"
           call localorb_info(info_str,use_unit,'(A)')
           call aims_stop ()
       end if
    end do 
 
    ! 4. Now all coordinates should be cartesian -> work with them 
    coords = coords_temp 

  end subroutine all_atoms_to_cart 
  !******
    !----------------------------------------------------------------------------
  !****s* geometry/min_atomic_dist_old
  !  NAME
  !    min_atomic_dist_old
  !  SYNOPSIS

  ! THIS VERSION DOES NOT WORK FOR ALL PERIODIC SYSTEMS!!
  ! TAKE SI, CONVENTIONAL CELL, AND THE SECOND ATOM AT FRAC COORDS 0.3, 0.3, 0.3
  ! VERSION IS KEPT HERE ONLY FOR ARCHIVAL PURPOSES, FOR NOW!

  subroutine min_atomic_dist_old(lattice_vector, coords, dist, i_atom, j_atom)

    !  PURPOSE
    !
    !    Return minimal distance between two atoms.
    !
    !  USES

    use dimensions
    use bravais, only: get_map_to_center_cell_matrix
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(IN) :: coords(3,n_atoms)
    real*8, intent(OUT) :: dist
    integer, intent(OUT), optional :: i_atom, j_atom

    !  INPUTS
    !   o lattice_vector -- Bravais vectors
    !   o coords -- Atomic coordinates
    !  OUTPUTS
    !   o dist -- Minimal distance between two atoms
    !   o i_atom, j_atom (optional) -- Atoms with minimal distance.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i, j, i_periodic
    real*8 :: thissq, minsq
    real*8 :: diff_vec(3), length(n_periodic)
    real*8 :: map_to_center_cell_matrix(n_periodic, 3)
    character(*), parameter :: func = 'min_atomic_dist'

    ! --- Initialization

    if (present(i_atom)) i_atom = 1
    if (present(j_atom)) j_atom = 1
    minsq = huge(minsq)

    if (n_periodic > 0) then
       ! Cannot use map_to_center_cell() because of custom lattice_vector.
       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
       &                                  map_to_center_cell_matrix)
    end if

    ! --- Distances to own images

    ! JW: For badly chosen Bravais vectors, this would be incorrect.
    do i_periodic = 1, n_periodic
       minsq = min(minsq, sum(lattice_vector(:, i_periodic)**2))
    end do

    ! --- Distances from atoms to centers

    ! THE FOLLOWING IS WHERE THIS GOES WRONG ... but even if we fixed it,
    ! the loop over unit cells above would make sure we will miss any
    ! self-images of atoms in different unit cells.

    do i = 1, n_atoms
       do j = i+1, n_atoms
          diff_vec = coords(:, j) - coords(:, i)
          if (n_periodic > 0) then    ! Map to center cell
             length = matmul(map_to_center_cell_matrix, diff_vec)
             length = length - nint(length)
             diff_vec = matmul(lattice_vector, length)
          end if
          thissq = sum(diff_vec**2)
          if (thissq < minsq) then
             if (present(i_atom)) i_atom = i
             if (present(j_atom)) j_atom = j
             minsq = thissq
          end if
       end do
    end do
    dist = sqrt(minsq)

  end subroutine min_atomic_dist_old
  !******
  !----------------------------------------------------------------------------
  !****s* geometry/min_atomic_dist
  !  NAME 
  !    min_atomic_dist
  !  SYNOPSIS
  !
  subroutine min_atomic_dist(lattice_vector, coords, dist)

    !  PURPOSE
    !
    !    Wrapper to min_grid_centre_dist that removes need for atom/grid/cell indices
    !
    !  USES

    use dimensions, only: n_periodic, n_atoms 
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(IN) :: coords(3,n_atoms)
    real*8, intent(OUT) :: dist

    !  INPUTS
    !   o lattice_vector -- Bravais vectors
    !   o coords -- Atomic coordinates
    !  OUTPUTS
    !   o dist -- Minimal distance between two grid-centres (so includes empty sites)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Created Sept2017 by AJL as a wrapper to min_grid_centre_dist
    !  SOURCE


    integer :: out_grid_1, out_grid_2, cell_grid_1, cell_grid_2, cell_grid_3
    real*8  :: dist_physical
    integer :: out_physical_1, out_physical_2, cell_physical_1, cell_physical_2, cell_physical_3

    call min_grid_centre_dist(lattice_vector, coords, dist, out_grid_1, out_grid_2, &
                                  cell_grid_1, cell_grid_2, cell_grid_3, &
                                  dist_physical, out_physical_1, out_physical_2, &
                                  cell_physical_1, cell_physical_2, cell_physical_3)

  end subroutine min_atomic_dist
  !******
  !----------------------------------------------------------------------------
  !****s* geometry/min_grid_centre_dist
  !  NAME
  !    min_grid_centre_dist
  !  SYNOPSIS

  ! AJL: Renamed to reflect that we have grid centres in the structure, not *just* atoms.
  ! However, still calculates distances between atom centres as well as grid centres.
  subroutine min_grid_centre_dist(lattice_vector, coords, dist, out_grid_1, out_grid_2, &
                                  cell_grid_1, cell_grid_2, cell_grid_3, &
                                  dist_physical, out_physical_1, out_physical_2, &
                                  cell_physical_1, cell_physical_2, cell_physical_3)

    !  PURPOSE
    !
    !    Return minimal distance between two atoms / grid centres
    !
    !  USES

    use dimensions
    use localorb_io,  only: use_unit
    use bravais, only: get_map_to_center_cell_matrix
    use species_data, only: species_pseudoized, no_basis
    use mpi_tasks,    only: aims_stop
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(IN) :: coords(3,n_atoms)
    real*8, intent(OUT) :: dist
    integer, intent(OUT) :: out_grid_1, out_grid_2
    integer, intent(OUT) :: cell_grid_1, cell_grid_2, cell_grid_3
    real*8, intent(OUT) :: dist_physical
    integer, intent(OUT) :: out_physical_1, out_physical_2
    integer, intent(OUT) :: cell_physical_1, cell_physical_2, cell_physical_3

    !  INPUTS
    !   o lattice_vector -- Bravais vectors
    !   o coords -- Atomic coordinates
    !  OUTPUTS
    !   o dist -- Minimal distance between two grid-centres (so includes empty sites)
    !   o out_grid_1, out_grid_2 -- Grid-centres with minimal distance.
    !   o cell_grid_1, cell_grid_2, cell_grid_3 -- Unit cell images where minimum-distance grid-centres were found.
    !   Added Sept2017
    !   o dist_physical -- Minimal distance between just physically meaningful centres (atoms and pseudocores)
    !   o out_physical_1, out_physical_2 -- Physically meaningful centres with minimal distance
    !   o cell_physical_1, cell_physical_2, cell_phyisical_3 -- Unit cell images where minimum distance
    !   physically meaningful centres were found.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2013).
    !    Edited Sept2017 by AJL to differentiate between atomic centres and all-grid centres
    !  SOURCE

    integer :: i_periodic
    integer :: i_atom, j_atom
    integer :: n_shell
    integer :: n_1, n_2, n_3
    real*8 :: minsq, thissq
    real*8 :: minsq_physical

    real*8 :: previous_minsq
    real*8 :: previous_minsq_physical
    real*8 :: length(3)

    integer :: safe_shell_convergence

    real*8 :: diff_vec(3)
    real*8 :: map_to_center_cell_matrix(n_periodic, 3)
    real*8, dimension(:,:), allocatable :: coords_temp

    character*150 :: info_str

    character(*), parameter :: func = 'min_grid_centre_dist'

    ! --- Initialization

    ! Sure, the coords_temp array can be saved altogether in non-periodic
    ! systems. If this ever causes memory problems, please eliminate coords_temp
    ! if you can.
    allocate(coords_temp(3,1:n_atoms))

    out_grid_1 = 1 
    out_grid_2 = 1
    cell_grid_1 = 0
    cell_grid_2 = 0
    cell_grid_3 = 0

    out_physical_1 = 1
    out_physical_2 = 1
    cell_physical_1 = 0
    cell_physical_2 = 0
    cell_physical_3 = 0

    minsq = huge(minsq)
    minsq_physical = huge(minsq_physical)

    if (n_periodic > 0) then
       ! Cannot use map_to_center_cell() because of custom lattice_vector.
       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
       &                                  map_to_center_cell_matrix)

       do i_atom = 1, n_atoms, 1
          length = matmul (map_to_center_cell_matrix,coords(:,i_atom))
          length = length - nint(length)
          coords_temp(:,i_atom) = matmul (lattice_vector, length)
       enddo

    else
       coords_temp = coords
    end if

    ! Nice idea but running over arbitrary dimensions systematically
    ! is best done by a recursive algorithm ... or much more simply
    ! by copying and simplifying our loops. Please implement if needed,
    ! not hard to do.
    if ( (n_periodic.ne.0) .and. (n_periodic.ne.3) ) then
       write(info_str, '(A)') &
       & 'Error: This subroutine currently ONLY supports zero-', &
       & 'or three-dimensional periodic geometries.'
       call aims_stop(info_str, func)
    end if

    ! --- Distances to own images
    ! Initialize minsq. This may not strictly be necessary, but it puts us
    ! on safe ground.
    do i_periodic = 1, n_periodic
       minsq = min(minsq, sum(lattice_vector(:, i_periodic)**2))
    end do

    ! For non-periodic systems the rest is trivial. Periodic systems are
    ! more ... confounding.

    ! The primary problem is that unit cell shapes can be truly awkward, with
    ! very non-orthogonal angles between them. Thus, even in a primitive Bravais
    ! lattice one can easily construct unit cells where the shortest connecting
    ! lattice vector is given by (3*a_1 - 4*a_2) or similar.
    ! Since lattice vectors are integer linear combinations of the
    ! unit cell vectors specified in geometry.in, many clever strategies to find
    ! the shortest linear combinations anaytically will fail because they implicitly
    ! rely on continuous coordinates somehow.
    ! Thus, the best algorithm that I could think of is the trivial one - just
    ! enumerate unit cells in shell after shell until we are safely converged.

    ! The second problem is that the atoms in the unit cell can be located
    ! anywhere they like. We can map them back into the Wigner-Seitz cell,
    ! but that will still leave corner cases for unpleasantly chosen unit
    ! cell geometries. Thus, we just iterate over all those atoms too.

    ! I know we could clean the math for elegance. If it becomes necessary, please do.
    ! Sometimes, however, it is better to have a safe solution, even if it is not
    ! very aesthetic.

    ! Comment recorded after the fact:
    !
    ! What would be even better is to first cycle only over possible lattice vectors,
    ! shell by shell, as done below.
    !
    ! From among those, find the triple of lattice vectors R_1, R_2, R_3 that
    !
    ! (a) yields the volume of only the primitive cell (not 2 or more, not zero)
    !
    ! (b) has the smallest sum of squares (R_1)^2 + (R_2)^2 + (R_3)^2 .
    !
    ! This would have to be the "most orthogonal" triple of lattice vectors.
    !
    ! This triple of lattice vectors would be the one to use for real-space
    ! lists of atoms with roughly equal extent in all directions.
    !
    ! The above algorithm is not implemented, but if, for any reason, we see
    ! sustained trouble with inconvenient unit cell choices in the future,
    ! we should consider it.

    n_shell = 0
    safe_shell_convergence = 2 ! This means that we iterate over an extra
                               ! shell of unit cells before we declare convergence.
                               ! If this becomes a time problem for large structures
                               ! there are certainly better ways. I (VB) just do not see
                               ! that yet.
    previous_minsq = huge(previous_minsq)
    previous_minsq_physical = huge(previous_minsq_physical)
    do while (safe_shell_convergence.gt.0)

       ! run over shells of unit cells
       do n_1 = -n_shell, n_shell, 1
          do n_2 = -n_shell, n_shell, 1
             do n_3 = -n_shell, n_shell, 1

                if ( (n_1.ne.n_shell) .and. (n_2.ne.n_shell) .and. (n_3.ne.n_shell) ) then
                   ! we were here before - skip this combination of unit cells
                   continue
                else

                   do i_atom = 1, n_atoms, 1

                      if (n_periodic.ne.0) then
                         length(:) = coords_temp(:,i_atom) + n_1 * lattice_vector(:,1) &
                                                           + n_2 * lattice_vector(:,2) &
                                                           + n_3 * lattice_vector(:,3)
                      else
                         length(:) = coords_temp(:,i_atom)
                      end if

                      ! Can NOT do just the upper triangle here, in case of periodic images.
                      do j_atom = 1, n_atoms, 1

                         ! Evaluate all this EXCEPT if we are in the zeroth unit cell,
                         ! where we must skip any self-images
                         if (.not.( (i_atom.eq.j_atom) .and. &
                                    (n_1.eq.0) .and. (n_2.eq.0) .and. (n_3.eq.0) )) then
                            diff_vec = coords_temp(:, j_atom) - length(:)
                            thissq = sum(diff_vec**2)
                            if (thissq < minsq) then
                               out_grid_1 = i_atom
                               out_grid_2 = j_atom
                               cell_grid_1 = n_1
                               cell_grid_2 = n_2
                               cell_grid_3 = n_3
                               minsq = thissq
                            end if
                            ! AJL/Sept2017
                            ! Store here information only if the centres are physically meaningful
                            ! Empty is not a subset of no_basis, as empty can have a basis and no_basis
                            ! can have a charge!
                            ! For now, as is lifted from close_encounters in pbc_lists.f90
                            if ( (.not.(empty(i_atom).and.species_pseudoized(species(j_atom)))) .and. &
                                 (.not.(species_pseudoized(species(i_atom)).and.empty(j_atom))) .and. &
                                 (.not.(no_basis(species(i_atom)).and.(no_basis(species(j_atom))))) .and. & !) then
                                 (thissq < minsq_physical) ) then
                               out_physical_1 = i_atom
                               out_physical_2 = j_atom
                               cell_physical_1 = n_1
                               cell_physical_2 = n_2
                               cell_physical_3 = n_3
                               minsq_physical = thissq
                            end if
                         end if

                      enddo
                   enddo

                end if

             enddo ! n_3
          enddo ! n_2
       enddo ! n_1

       if (n_periodic.eq.0) then
          ! we do not need to run over any shells, just state that we are done.
          safe_shell_convergence = 0
       else
         ! check shell convergence
         if (minsq .lt. previous_minsq) then
            ! we found something new - do another shell.
            previous_minsq = minsq
            previous_minsq_physical = minsq_physical
            safe_shell_convergence = 2
         else
            safe_shell_convergence = safe_shell_convergence - 1
         end if
         previous_minsq = minsq
         previous_minsq_physical = minsq_physical
       end if

        n_shell = n_shell + 1
    enddo

    ! Reset to last shell used, even if only needed for test output.
    n_shell = n_shell - 1

    ! now we are done - here is the minimum distance found ...
    dist = sqrt(minsq)
    ! AJL, Sept2017: And just for atomic (physical) centres
    dist_physical = sqrt(minsq_physical)

    deallocate(coords_temp)

    ! test output only
    !if (myid.eq.0) then
    !
    !  write (use_unit, *) " | Minimum distance found: ", dist*bohr
    !  write (use_unit, *) " | Shells used:            ", n_shell - 1
    !  write (use_unit, *) " | Closest atoms:          ", out_grid_1, out_grid_2
    !  write (use_unit, *) " | Unit cell of atom 1:    ", cell_grid_1, cell_grid_2, cell_grid_3
    !
    !end if
    ! end test

  end subroutine min_grid_centre_dist
  !******
  !----------------------------------------------------------------------------
  !****s* geometry/min_multipole_distance
  !  NAME
  !    min_multipole_distance
  !  SYNOPSIS

  subroutine min_multipole_dist &
  & ( lattice_vector, occ_coords, species, atom_radius_sq, multipole_coords, &
  &   empty_coords, dist, min_atom, min_multipole, offending_multipole )

    !  PURPOSE
    !
    !    Determine distance between each external embedding multipole and closest atom.
    !    If too close to the nearest atom, report back to calling subroutine.
    !
    !  USES

    use dimensions
    use bravais, only: get_map_to_center_cell_matrix
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: lattice_vector(3, n_periodic)
    real*8, intent(IN) :: occ_coords(3,n_occ_atoms)
    integer, intent(IN) :: species(n_atoms)
    real*8, intent(IN) :: atom_radius_sq(n_species)
    real*8, intent(IN) :: multipole_coords(3,n_multipoles)
    real*8, intent(IN) :: empty_coords(3,n_empty_atoms)
    real*8, intent(OUT) :: dist
    integer, intent(OUT) :: min_atom, min_multipole
    integer, dimension(n_multipoles), intent(OUT) :: offending_multipole

    !  INPUTS
    !   o lattice_vector -- Bravais vectors
    !   o occ_coords -- Occupied atomic coordinates
    !   o species -- Species number of each atom
    !   o atom_radius_sq -- square of the radius of the most extended basis function of each species
    !   o multipole_coords -- Multipole coordinates
    !   o empty_coords -- Empty site coordinates
    !  OUTPUTS
    !   o dist -- Minimal distance between two atoms
    !   o i_atom, i_multipole -- Atom and multipole with minimal distance.
    !   o offending_multipole -- List (number in order of appearance in geometry.in) of all
    !                            multipoles that are too close to an atom.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i_multipole, i_atom, i_periodic
    integer :: n_offending_multipoles
    real*8 :: thissq, minsq, thissq_empty, minsq_empty
    real*8 :: diff_vec(3), length(n_periodic)
    real*8 :: map_to_center_cell_matrix(n_periodic, 3)

    real*8 :: safety_margin  ! hardwired here for now but should become configurable later.

    character(*), parameter :: func = 'min_multipole_dist'

    ! --- Initialization

    min_atom = 0
    min_multipole = 0
    n_offending_multipoles = 1       ! real value is one less
    offending_multipole = 0
    minsq = huge(minsq)

    safety_margin = 1.0  ! Value is chosen in units of bohr radii .

    if (n_periodic > 0) then
       ! Cannot use map_to_center_cell() because of custom lattice_vector.
       call get_map_to_center_cell_matrix(n_periodic, lattice_vector, &
       &                                  map_to_center_cell_matrix)
    end if

    ! --- Distances to own images

    ! What we wish to find


    do i_periodic = 1, n_periodic
       minsq = min(minsq, sum(lattice_vector(:, i_periodic)**2))
    end do

    ! --- Distances from multipoles to centers

    do i_multipole = 1, n_multipoles

       ! these need resetting every loop
       minsq_empty = huge(minsq_empty)
       do i_periodic = 1, n_periodic
          minsq_empty = min(minsq_empty, sum(lattice_vector(:, i_periodic)**2))
       end do

       ! check if there is an empty site on the multipole
       do i_atom = 1, n_empty_atoms

          ! collect distance
          diff_vec = empty_coords(:, i_atom) - multipole_coords(:, i_multipole)
          ! correct for periodicity
          if (n_periodic > 0) then    ! Map to center cell
             length = matmul(map_to_center_cell_matrix, diff_vec)
             length = length - nint(length)
             diff_vec = matmul(lattice_vector, length)
          end if

          thissq_empty = sum(diff_vec**2)

          if (thissq_empty < minsq_empty) then
             ! set shortest distance between multipole and empty sites
             minsq_empty = thissq_empty
          end if

       end do

       if (minsq_empty.gt.0.0d0) then
          ! check proximity to other atoms if no empty site is present

          do i_atom = 1, n_occ_atoms
             ! collect distance
             diff_vec = occ_coords(:, i_atom) - multipole_coords(:, i_multipole)
             ! correct for periodicity
             if (n_periodic > 0) then    ! Map to center cell
                length = matmul(map_to_center_cell_matrix, diff_vec)
                length = length - nint(length)
                diff_vec = matmul(lattice_vector, length)
             end if

             thissq = sum(diff_vec**2)

             ! ideally we would take the square root of these numbers,
             ! then add the safety_margin and square again, as the relationship is non-linear,
             ! but for now this comparison is safe enough, and simpler.

             if (thissq .lt. ( atom_radius_sq(species(i_atom)) + safety_margin ) ) then
                ! this multipole is too close to the quantum region and needs an empty site
                ! with grids on top of it. Mark as an offending multipole.
                offending_multipole(n_offending_multipoles) = i_multipole
             end if

             ! returning the closest multipole is actually quite redundant,
             ! but we will leave it here for now as we've already computed the distances.

             if (thissq < minsq) then
                ! this is checked irrespective of whether we are within the basis set radius
                min_atom = i_atom
                min_multipole = i_multipole
                minsq = thissq
             end if

          end do

          if (offending_multipole(n_offending_multipoles).eq.i_multipole) then
             ! increase n_offending_multipoles if this value has been added
             n_offending_multipoles = n_offending_multipoles + 1
          end if

       end if

    end do
    dist = sqrt(minsq)

  end subroutine min_multipole_dist
  !**********
end module geometry
