! subroutine to create lookup-table atom_type to map
! atom onto initialization type for initial_rho
!
! R.Gehrke (2007)


subroutine create_ini_type_lookup( initial_moment, initial_charge, kind_of_initial, species, empty, &
                                   type_charge, type_moment, type_kind, type_species, atom_type)

  use dimensions
  use runtime_choices
  use localorb_io
  use species_data

  implicit none
  
  ! imported variables

  ! input 
  real*8,dimension(n_atoms) :: initial_moment
  real*8,dimension(n_atoms) :: initial_charge
  integer, dimension(n_atoms) :: kind_of_initial
  integer, dimension(n_atoms) :: species
  logical, dimension(n_atoms) :: empty

  ! output
  real*8, dimension(n_atoms) :: type_charge
  real*8, dimension(n_atoms) :: type_moment
  integer, dimension(n_atoms) :: atom_type 
  integer, dimension(n_atoms) :: type_species 
  integer, dimension(n_atoms) :: type_kind

  ! local variables
  logical :: found
  character*100 :: info_str

  ! counter
  integer :: i_atom
  integer :: i_type

  ! FIXME: n_ini_type should be determined in parse_geo in dimensions.f - 
  !        not here.

  n_ini_type = 0

  do i_atom = 1, n_atoms, 1
  if (.not.empty(i_atom)) then
     
     ! look, if this combination of moment, charge, species and kind does already exist
     if (spin_treatment.eq.1) then
        found = .false.
        do i_type = 1, n_ini_type, 1
           if ((initial_moment(i_atom) .eq. type_moment(i_type)) .and. (initial_charge(i_atom) .eq. type_charge(i_type)) &
                .and. (species(i_atom) .eq. type_species(i_type)) .and. (kind_of_initial(i_atom) .eq. type_kind(i_type))) then
              found = .true.
              atom_type(i_atom) = i_type
           end if
        end do
        if (.not.found) then
           n_ini_type = n_ini_type + 1
           atom_type(i_atom) = n_ini_type
           type_moment(n_ini_type) =  initial_moment(i_atom)
           type_charge(n_ini_type) =  initial_charge(i_atom)
           type_species(n_ini_type) = species(i_atom)
           type_kind(n_ini_type) = kind_of_initial(i_atom)
        end if
     else if (spin_treatment.eq.0) then
        found = .false.
        do i_type = 1, n_ini_type, 1
           if ((initial_charge(i_atom) .eq. type_charge(i_type)) .and. (species(i_atom) .eq. type_species(i_type))) then
              found = .true.
              atom_type(i_atom) = i_type
           end if
        end do
        if (.not.found) then
           n_ini_type = n_ini_type + 1
           atom_type(i_atom) = n_ini_type
           type_charge(n_ini_type) = initial_charge(i_atom)
           type_species(n_ini_type) = species(i_atom)
        end if
     end if

  end if  ! .not. empty(i_atom)
  end do

end subroutine create_ini_type_lookup
