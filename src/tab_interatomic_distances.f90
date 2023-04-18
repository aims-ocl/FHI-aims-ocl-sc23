!****s* FHI-aims/tab_interatomic_distances
!  NAME
!   tab_interatomic_distances
!  SYNOPSIS

subroutine tab_interatomic_distances(n_atom_list, atom_list, atom_atom_tab)

!  PURPOSE
!  calculate the interatomic distances between all pairs of n_atom_list and store in atom_atom_tab
!  The atom_list is some list that, n_atom_list long, that points to the proper indices in coords_center.
!
!  USES

  use localorb_io, only: localorb_info, use_unit, OL_norm
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer, intent(in) :: n_atom_list
  real*8, dimension(n_atom_list,n_atom_list) :: atom_atom_tab
  integer, dimension(n_atom_list) :: atom_list


!  INPUTS
!    o n_atom_list -- number of relevant atoms
!    o atom_list -- list of relevant atoms
!  OUTPUT
!    o atom_atom_tab -- interatomic distances for all n_atom_list atoms, must be indexed with atom_atom_index
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


  integer :: i1, i2, i_center_1, i_center_2, i_center_L1, i_center_L2
  integer, dimension(n_atom_list) :: center_index
  real*8, dimension(3) :: delta
  character*140 :: info_str

! start calculating the original atom_atom_tab matrix
  
  do i1 = 1, n_atom_list
     i_center_1  = atom_list(i1)
     do i2 = i1, n_atom_list

        i_center_2  = atom_list(i2)

        delta(1) = coords_center(1,i_center_1)-coords_center(1,i_center_2)
        delta(2) = coords_center(2,i_center_1)-coords_center(2,i_center_2)
        delta(3) = coords_center(3,i_center_1)-coords_center(3,i_center_2)

        atom_atom_tab(i1,i2) = sqrt(delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3))
        atom_atom_tab(i2,i1) = atom_atom_tab(i1,i2)

     end do
  end do
  write(info_str,'(2X,A,1X,I10,A,I10,1X,A)') &
            "| Original list of interatomic distances with ", n_atom_list ," x ",n_atom_list, "entries is created."

  call localorb_info(info_str,use_unit,'(A)',OL_norm)
end subroutine tab_interatomic_distances

! end
!******

!****s* FHI-aims/tab_interatomic_distances_count_entries
!  NAME
!   tab_interatomic_distances_count_entries
!  SYNOPSIS
subroutine tab_interatomic_distances_count_entries(n_atom_list, atom_list, n_atom_atom_tab)

!  PURPOSE
!  find size of the array atom_atom_tab
!  The atom_list is some list that, n_atom_list long, that points to the proper indices in coords_center.
!
!  USES

  use pbc_lists
  use dimensions, only: n_species
  use species_data
  use grids

  implicit none

!  ARGUMENTS

  integer, intent(in) :: n_atom_list
  integer, dimension(n_atom_list) :: atom_list
  integer, intent(out) :: n_atom_atom_tab

!  INPUTS
!    o n_atom_list -- number of relevant atoms
!    o atom_list -- list of relevant atoms
!  OUTPUT
!    o n_atom_atom_tab -- length of array of all relevant distances for all n_atom_list atoms
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

  integer :: i1, i2, i_center_1, i_center_2, i_center_L1, i_center_L2
  real*8, dimension(3) :: delta
  integer :: counter
  real*8 :: dist_sq

  integer :: i_species_1, i_species_2
  real*8, dimension(n_species,n_species) :: cut_radius_sq
  real*8 :: trial_radius

  ! This routine only counts the dimension of an array that is later allocated and then filled with life
  ! in tab_interatomic_distances_local below.
  ! The counting in the present routine and in tab_interatomic_distances_local must be exactly the
  ! same, or else we will run into problems.

  ! Define cutoff radius beyond which atom pairs will never be needed for OUR version of
  ! the Stratmann et al. partition table. Note that our version is bounded in a way that
  ! should ensure a zero partition table if an atom is further away than its outermost 
  ! prescribed radius ("outer_partition_radius") from the current grid point or
  ! if two atoms are further away than the sum of their outer_partition_radius values
  ! (whichever is greater).

  ! However the list of centers over which we run is determined by a different and global criterion.
  !  
  ! Therefore, the furthest distance over which two atoms can "interact" in the partition table
  ! is here set to sum the maximum radii possible - out of the outermost integration grid shell
  ! of each species and the outer_partition_radius of each species. This is theoretically too
  ! large but for most species defaults should not be a significant impediment.

  trial_radius = 0.
  do i_species_1 = 1, n_species, 1
        trial_radius = max(trial_radius,outer_partition_radius(i_species_1))
        trial_radius = max(trial_radius,r_radial(n_radial(i_species_1),i_species_1))
  enddo
  trial_radius = (2.d0*trial_radius)**2
  
  do i_species_1 = 1, n_species, 1
    do i_species_2 = 1, n_species, 1
      
      if (i_species_1.le.i_species_2) then
         cut_radius_sq(i_species_1,i_species_2) = trial_radius
      else
        ! we already calculated the sum of outer partition radius for species 1 and outer grid radius
        ! for species 2, now see whether the switched sum is larger.
        ! In any case, take the bigger one.
        if (trial_radius .gt. cut_radius_sq(i_species_2,i_species_1) ) then
           cut_radius_sq(i_species_1,i_species_2) = trial_radius
           cut_radius_sq(i_species_2,i_species_1) = trial_radius
        else
           cut_radius_sq(i_species_1,i_species_2) = cut_radius_sq(i_species_2,i_species_1)
        end if
      end if

    end do
  end do

  ! now count the number of atom pairs for which we will ever need the interatomic distance.

  n_atom_atom_tab = 0
  
  do i1 = 1, n_atom_list
     i_center_1  = atom_list(i1)
     do i2 = 1, n_atom_list

        if (i1 .ne. i2) then

           i_center_2  = atom_list(i2)

           delta(1) = coords_center(1,i_center_1)-coords_center(1,i_center_2)
           delta(2) = coords_center(2,i_center_1)-coords_center(2,i_center_2)
           delta(3) = coords_center(3,i_center_1)-coords_center(3,i_center_2)

           dist_sq =  delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3)

           if (dist_sq .le. cut_radius_sq (species_center(i1),species_center(i2)) ) then
              n_atom_atom_tab = n_atom_atom_tab + 1
           end if
        end if
     end do
  end do
end subroutine tab_interatomic_distances_count_entries

!******

!****s* FHI-aims/tab_interatomic_distances_local
!  NAME
!   tab_interatomic_distances_local
!  SYNOPSIS
subroutine tab_interatomic_distances_local(n_atom_list, atom_list, n_atom_atom_tab, &
                                           atom_atom_dist_list, atom_idx_A, atom_idx_B)
!  PURPOSE
!  calculate the interatomic distances between all pairs of n_atom_list and store in atom_atom_tab
!  The atom_list is some list that, n_atom_list long, that points to the proper indices in coords_center.
!
!  USES

  use pbc_lists
  use dimensions, only: n_species
  use species_data
  use grids
  use localorb_io, only: localorb_info

  implicit none

!  ARGUMENTS

  integer, intent(in) :: n_atom_list
  integer, dimension(n_atom_list) :: atom_list
  integer, intent(in) :: n_atom_atom_tab
  real*8, dimension(n_atom_atom_tab) :: atom_atom_dist_list
  integer, dimension(n_atom_atom_tab) :: atom_idx_A
  integer, dimension(n_atom_list+1) :: atom_idx_B

!  INPUTS
!    o n_atom_list -- number of relevant atoms
!    o atom_list -- list of relevant atoms
!    o n_atom_atom_tab -- length of array of all relevant distances for all n_atom_list atoms
!  OUTPUT
!    o atom_atom_dist_list -- interatomic distances for all n_atom_list atoms, must be indexed with atom_atom_index
!    o atom_idx_A       -- index for atom A in array atom_atom_dist_list
!    o atom_idx_B       -- array of indices pointing to the starting entry for atom B in array atom_atom_dist_list
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

  integer :: i1, i2, i_center_1, i_center_2, i_center_L1, i_center_L2
  real*8, dimension(3) :: delta
  integer :: counter
  real*8 :: memory_usage, memory_usage_trad
  integer :: i_species_1, i_species_2
  real*8, dimension(n_species,n_species) :: cut_radius_sq
  real*8 :: trial_radius
  real*8 :: dist_sq
  real*8 :: distance
  character*140 :: info_str

  ! Define cutoff radius beyond which atom pairs will never be needed for OUR version of
  ! the Stratmann et al. partition table. Note that our version is bounded in a way that
  ! should ensure a zero partition table if an atom is further away than its outermost 
  ! prescribed radius ("outer_partition_radius") from the current grid point.
  !
  ! Therefore, the furthest distance over which two atoms can "interact" in the partition table
  ! is the sum of the outermost integration grid shell of one species 
  ! and the outer_partition_radius of the other.

  ! However the list of centers over which we run is determined by a different and global criterion.
  !  
  ! Therefore, the furthest distance over which two atoms can "interact" in the partition table
  ! is here set to sum the maximum radii possible - out of the outermost integration grid shell
  ! of each species and the outer_partition_radius of each species. This is theoretically too
  ! large but for most species defaults should not be a significant impediment.

  trial_radius = 0.
  do i_species_1 = 1, n_species, 1
        trial_radius = max(trial_radius,outer_partition_radius(i_species_1))
        trial_radius = max(trial_radius,r_radial(n_radial(i_species_1),i_species_1))
  enddo
  trial_radius = (2.d0*trial_radius)**2
  
  do i_species_1 = 1, n_species, 1
    do i_species_2 = 1, n_species, 1

      if (i_species_1.le.i_species_2) then
         cut_radius_sq(i_species_1,i_species_2) = trial_radius
      else
        ! we already calculated the sum of outer partition radius for species 1 and outer grid radius
        ! for species 2, now see whether the switched sum is larger.
        ! In any case, take the bigger one.
        if (trial_radius .gt. cut_radius_sq(i_species_2,i_species_1) ) then
           cut_radius_sq(i_species_1,i_species_2) = trial_radius
           cut_radius_sq(i_species_2,i_species_1) = trial_radius
        else
           cut_radius_sq(i_species_1,i_species_2) = cut_radius_sq(i_species_2,i_species_1)
        end if
      end if

    end do
  end do

! LN: start calculating new atom_atom_dist_list
  counter = 1
  do i1 = 1, n_atom_list
     i_center_1 = atom_list(i1)
     atom_idx_B(i1) = counter
     do i2 = 1, n_atom_list

        if (i1 .ne. i2) then
           i_center_2  = atom_list(i2)

           delta(1) = coords_center(1,i_center_1)-coords_center(1,i_center_2)
           delta(2) = coords_center(2,i_center_1)-coords_center(2,i_center_2)
           delta(3) = coords_center(3,i_center_1)-coords_center(3,i_center_2)

           dist_sq = delta(1)*delta(1)+delta(2)*delta(2)+delta(3)*delta(3)

           if (dist_sq .le. cut_radius_sq (species_center(i1),species_center(i2)) ) then
              atom_atom_dist_list(counter) = sqrt(dist_sq)
              atom_idx_A(counter) = i2
              counter = counter + 1
           end if
        end if
     end do
  end do
  ! calculating the total memory usage of the sparse atom_atom_tab and the standard 2-dim matrix
  memory_usage = (n_atom_atom_tab*8.0 + (n_atom_list+1+n_atom_atom_tab)*4.0)/1000.0
  memory_usage_trad = (n_atom_list * n_atom_list * 8.0 ) / 1000.0

  write(info_str,'(2X,A,1X,F12.2,1X,A,1X,F12.2,1X,A)') &
     "| The sparse table of interatomic distances needs ", memory_usage , &
     "kbyte instead of", memory_usage_trad , "kbyte of memory."
  call localorb_info(info_str)

  if (memory_usage.gt.memory_usage_trad) then
     write(info_str,'(2X,A)') &
     "| Using the partition_type stratmann_smoother will reduce your memory usage."
     call localorb_info(info_str)
  end if

  atom_idx_B(n_atom_list+1) = counter
end subroutine tab_interatomic_distances_local
!******
