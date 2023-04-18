!  subroutine parse_control to determine sizes of arrays

subroutine parse_control()

  use cluster
  use pair_potential, pp_allocate => allocate
  use basin_hopping, bh_minimizer => minimizer

  implicit none

  ! local variables

  logical :: eof
  character*30 :: desc_str
  integer :: i_code
  integer :: n_data_temp

  ! counters

  integer :: i_species
  integer :: i_atom
  integer :: i_adsorp

  ! begin work

  write(6,'(2X, A)') "Parsing control.in.opt ..."

  potential_flag = 'external' 
  use_pair_pot   = .false.
  n_max_data = 0

  open(7, file = "control.in.opt", status='OLD', iostat = i_code)

  eof = .false.

  i_species = 0
  i_adsorp = 0

  if (i_code.ne.0) then
     write(*,*) "* Input file control.in.opt not found."
     stop
  end if
  
  do while (.not. eof)
     read(7,*,iostat = i_code) desc_str
     if (i_code.ne.0) then
        eof = .true.
     else if (desc_str.eq."species") then
        i_species = i_species + 1

     else if (desc_str.eq."adsorption_atom") then
        i_adsorp = i_adsorp + 1

     else if (desc_str.eq."potential") then

        backspace(7)
        read (7,*) desc_str, potential_flag

        select case (potential_flag)

        case ('NN')
           bh_minimizer = 'magic'
           
        case default
           continue

        end select

     else if (desc_str.eq."spin_polarized") then
        
        backspace(7)
        read (7,*) desc_str, spin_polarized

     else if (desc_str.eq."charged") then
        
        backspace(7)
        read (7,*) desc_str, charged

     else if (desc_str.eq."pair_pot") then
        use_pair_pot = .true.
        backspace(7)
        read (7,*) desc_str, desc_str, desc_str, n_data_temp
        n_max_data = max(n_max_data, n_data_temp)
     end if
  end do

  n_adsorption_atom = i_adsorp

  if (i_species.gt.0) then
     n_species = i_species
  end if
     
  close(7)

  if (potential_flag .eq. 'external') then

     write(6,'(2X, A)') "Parsing control.in ..."
  
     open(7, file = "control.in", status='OLD', iostat = i_code)
     
     eof = .false.
     
     i_species = 0
     
     if (i_code.ne.0) then
        write(*,*) "* Input file control.in not found."
        stop
     end if
     
     do while (.not. eof)
        read(7,*,iostat = i_code) desc_str
        if (i_code.ne.0) then
           eof = .true.
        else if (desc_str.eq."species") then
           i_species = i_species + 1
        end if
     end do
     
     close(7)
     
     if (i_species.gt.0) then
        n_species = i_species
     else
        write(*,*) "! No species in control.in"
        stop
     end if
     
  end if
  
  write(6,'(2X, A)') "Parsing geometry.in.basic ..."

  open(7, file="geometry.in.basic", status='OLD', iostat = i_code)
  
  eof = .false.

  if (i_code.ne.0) then
     write(*,*) "* Input file geometry.in.basic not found."
     stop
  end if

  i_atom = 0

  do while (.not. eof)
     read(7,*,iostat = i_code) desc_str
     if (i_code.ne.0) then
        eof = .true.
     else if (desc_str.eq."atom") then
        i_atom = i_atom + 1
     end if
  end do

  if (i_atom.gt.0) then
     n_atoms = i_atom
  else
     write(*,*) "! No atoms in geometry.in.basic"
     stop
  end if

  write (6,'(1X,A,I6)') "n_species = ", n_species
  write (6,'(1X,A,I6)') "n_atoms   = ", n_atoms
  write (6,'(1X,A,I6)') "n_adsorption_atom   = ", n_adsorption_atom

  if (use_pair_pot) then
     call pp_allocate()
  end if

  allocate(adsorption_atom(n_adsorption_atom))
  allocate(ad_radius_reduction(n_adsorption_atom))

end subroutine parse_control

