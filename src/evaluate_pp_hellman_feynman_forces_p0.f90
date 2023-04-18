! Hellman-Feynman forces on the pseudocores
! AT

subroutine evaluate_pp_hellman_feynman_forces_p0(rho, partition_tab, pseudocore_forces)

  use dimensions
  use runtime_choices
  use grids
  use spline
  use geometry
  use pseudodata
  use species_data
  use mpi_utilities
  use localorb_io
  use constants
  use synchronize_mpi

  implicit none
  ! imported variables
  
  ! input
  real*8, dimension(n_spin, n_full_points) :: rho
  real*8, dimension(n_full_points) :: partition_tab

! local variable

  real*8, dimension(3, n_pp_atoms) :: pp_ionic_forces 
  real*8, dimension(3, n_pp_atoms) :: pp_elec_forces 
!  real*8, dimension(3, n_atoms) :: pp_elec_forces 
  real*8, dimension(3, n_pp_atoms) :: pp_multipole_forces 

  real*8 :: point_term
  real*8 :: atomic_term

  integer :: i_multipole
  integer :: i_coord
  integer :: i_spin

  integer :: i_index
  integer :: i_my_batch

  integer :: i_full_points

  real*8, dimension(3) :: direction
  real*8, dimension(3) :: coord_current
  real*8 :: distance_squared, dist


  real*8 :: d_V_loc_dr
  real*8 :: conversion

  integer :: i_atom, i_pp_atom, i_pp_species

  real*8 :: Z, q, i_r_log

  character*120 :: info_str
  character*40 :: func = 'evaluate_pp_hellman_feynman_forces_p0'

  ! output
  real*8, dimension(3, n_atoms) :: pseudocore_forces 

  conversion = hartree / bohr  

  ! first we add up the interaction with all nulcei

  pp_ionic_forces = 0.d0
  pp_elec_forces = 0.d0
  pseudocore_forces = 0.d0


  do i_pp_atom = 1, n_pp_atoms  

      q = pp_charge(pp_species(i_pp_atom))

      do i_atom = 1, n_real_atoms, 1
          ! just make sure, that we don't hit a ghost atom
          if (empty(i_atom)) cycle
 
          Z = species_z(species(i_atom))
            
          direction(:) = coords(:,i_atom) - pp_coords(:,i_pp_atom)

          distance_squared = 0.d0 
          do i_coord = 1,3,1
              distance_squared = distance_squared + (direction(i_coord))**2
          enddo

          !the plus in the line below is right - I checked that
          pp_ionic_forces(:,i_pp_atom) = pp_ionic_forces(:,i_pp_atom) - &
            Z*q / (distance_squared)**1.5d0*direction(:)

      enddo

  enddo


  i_full_points = 0

  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              !     get current integration point coordinate
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

!              i_atom = batches(i_my_batch)%points(i_index)%index_atom


              do i_spin = 1, n_spin, 1
                point_term = partition_tab(i_full_points) * rho(i_spin, i_full_points)

                do i_pp_atom = 1,n_pp_atoms

                 i_pp_species = pp_species(i_pp_atom)



              ! now find distance of pseudocore to integration point

                 distance_squared = 0.d0
                 do i_coord = 1,3,1

                   direction(i_coord) = coord_current(i_coord) - pp_coords(i_coord, i_pp_atom) 

                   distance_squared = distance_squared + &
                      (direction(i_coord))**2

                 enddo

		 d_V_loc_dr = 0.d0

                 if(distance_squared.gt.(localpot_outer_radius(i_pp_species))**2 ) then
                 ! extrapolate the local potential according to Coulomb character
	     	    d_V_loc_dr = pp_charge(i_pp_species)/distance_squared**1.5d0


                 elseif (distance_squared.lt.(pp_r_grid_min(i_pp_species))**2) then
                    d_V_loc_dr = 0.d0
                 else

                    dist = sqrt(distance_squared)
                    i_r_log = invert_log_grid(dist, &
                         pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))


                    ! chain rule: df/dr = df/di *di/dr

                    !df/di:
		    d_V_loc_dr = &
                         val_spline_deriv( i_r_log, &
                         local_pseudopot_spl(1,1,i_pp_species),  &
                         n_points_pp_fn(i_pp_species))

                   !df/di *di/dr:
                    d_V_loc_dr = d_V_loc_dr/log(pp_r_grid_inc(i_pp_species))/dist

                    !here we have to scale with 1/r, since vec(direction) is not normalized
                    d_V_loc_dr = d_V_loc_dr/dist


                 endif
 
                 do i_coord = 1,3,1

!                 pp_elec_forces(i_coord, i_atom) = pp_elec_forces(i_coord, i_atom) + &
!                    point_term*d_V_loc_dr*direction(i_coord)

                 pp_elec_forces(i_coord, i_pp_atom) = pp_elec_forces(i_coord, i_pp_atom) + &
                    point_term*d_V_loc_dr*direction(i_coord)

!                 pp_multipole_forces(i_coord, i_pp_atom) = pp_multipole_forces(i_coord, i_pp_atom) + &
!                               + (rho(i_spin, i_full_points) -  0.5*rho_multipole(i_full_points)) &
!                               * d_V_loc_dr*direction(i_coord) &
!                              * partition_tab(i_full_points)


                 enddo

!                 do i_coord = 1,3,1
!                 pp_multipole_forces(i_coord, i_pp_atom) = pp_multipole_forces(i_coord, i_pp_atom) + &
!                               + (rho(i_spin, i_full_points) -  0.5*rho_multipole(i_full_points)) &
!                               * d_V_loc_dr*direction(i_coord) &
!                               * partition_tab(i_full_points)
!                enddo



              enddo !i_pp_atom
             enddo !i_spin
           endif

        enddo

  enddo


  call sync_pp_charge_forces(pp_elec_forces)
!  call sync_pp_charge_forces(pp_multipole_forces)

  if (myid.eq.0) then

    write(info_str,'(A)') ''
    call localorb_info(info_str)

!    write (use_unit,'(2X,A)') "Forces acting on the pseudocore (local potential contributions):"
    write (use_unit,'(5X,A)') "ionic forces [eV/A]:"
    do i_pp_atom = 1, n_pp_atoms, 1
        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_pp_atom, pp_ionic_forces(:,i_pp_atom)*conversion 
    enddo

!    write (use_unit,'(5X,A)') "multipole forces [eV/A]:"
!    do i_atom = 1, n_pp_atoms, 1
!        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_pp_atom, pp_multipole_forces(:,i_atom)*conversion
!    enddo

    write (use_unit,'(5X,A)') "electronic forces [eV/A]:"
    do i_pp_atom = 1, n_pp_atoms, 1
        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_pp_atom, pp_elec_forces(:,i_pp_atom)*conversion
    enddo

    write (use_unit,'(5X,A)') "nonlocal forces [eV/A]:"
    do i_pp_atom = 1, n_pp_atoms, 1
        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_pp_atom, pp_nonlocal_forces(:,i_pp_atom)*conversion 
    enddo

    write (use_unit,'(5X,A)') "nlcc forces [eV/A]:"
    do i_pp_atom = 1, n_pp_atoms, 1
        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_pp_atom, pp_nlcc_forces(:,i_pp_atom)*conversion 
    enddo

    write (use_unit,'(5X,A)') "combined forces on pseudocore [eV/A]:"
    do i_pp_atom = 1, n_pp_atoms
        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_pp_atom, &
           (pp_ionic_forces(:,i_pp_atom) + pp_elec_forces(:,i_pp_atom) &
          + pp_nonlocal_forces(:,i_pp_atom) + pp_nlcc_forces(:,i_pp_atom)) &
          * conversion
    enddo

    write (use_unit,'(5X,A)') "nonlocal forces on atom [eV/A]:"
    do i_atom = 1, n_real_atoms
        write (use_unit,'(2X,A,I4,1X,3(E25.10,1X))') "|",i_atom, &
             (pp_nonlocal_forces_2(:,i_atom))*conversion
    enddo

  endif

 ! finally we have to put the pp_forces on the global array pseudocore_forces (with dimensions n_atoms)

 ! now, only the pseudocore_forces of those pseudopotentials close to the QM-region
  ! is added to the atomic forces

  i_pp_atom = 0
  do i_atom = 1, n_atoms,1
    if(species_pseudoized(species(i_atom))) then
        i_pp_atom = pp_atom2atom(i_atom)
      do i_coord = 1,3,1
         pseudocore_forces(i_coord,i_atom) = &
             pp_ionic_forces(i_coord,i_pp_atom) + pp_elec_forces(i_coord,i_pp_atom) + &
             pp_nonlocal_forces(i_coord,i_pp_atom) + pp_nlcc_forces(i_coord,i_pp_atom)
      enddo
    else
      do i_coord = 1,3,1
         pseudocore_forces(i_coord,i_atom) = pp_nonlocal_forces_2(i_coord,i_atom)
      enddo
    endif
  enddo



  if (i_pp_atom.ne.n_pp_in_qm) then
      if (myid.eq.0) then
       write(info_str,'(2X,A)') 'OOPS. Misscounting of pseudocores detected when forces are evaluated! AIMS ABORTING'
       call aims_stop(info_str,func)
      endif
  endif

end subroutine evaluate_pp_hellman_feynman_forces_p0
