!  NAME
!   read_geo
!  SYNOPSIS
    subroutine read_geo ( n_electrons )
!  PURPOSE
!  Subroutine read_geo provides cluster geometry, and calculates
!  related information (e.g. total number of Kohn-Sham eigenstates)
!  USES
      use aims_memory_tracking, only: aims_allocate, aims_deallocate
      use dimensions
      use runtime_choices
      use synchronize_mpi
      use geometry
      use species_data
      use plot
      use constraint
      use relaxation
      use mpi_tasks
      use localorb_io
      use constants
      use free_atoms
      use molecular_dynamics
      use ll_vdwdf
      use bravais
      use pseudodata
      use friction
      use esp_charges,only:esp_atom_constraint
      use xml_write
      use pbc_lists, only: check_for_close_encounters
      use distributed_hessian, only: check_hess_file
      implicit none
!   ARGUMENTS
      real*8 :: n_electrons
!   INPUTS
!    none
!   OUTPUTS
!    o n_electrons -- number of electrons in structure
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
!  eof : end-of-file marker
!  i_code : iostatus flag
!  desc_str : line descriptor
!  species_temp: temporary placeholder for species_name
!  found : flag indicating whether each species is actually defined.
!  mult : multiplicity of a species in structure

      logical eof
      logical:: MD_restart_file_exists = .false.
      integer i_code
      character*40 desc_str
      character*20 species_temp
      logical found
      integer mult

      !VA: in order to distinguish the relaxation constraints
      !    betwenn lattice vectors and atoms
      logical :: lv_found_before =.false.
      logical :: individual_lv_components_fixed =.false.
      logical :: individual_atom_components_frac_fixed =.false.

      integer region_temp
      real*8 :: electronsum

      real*8 :: species_z_max

      logical, dimension(n_region) :: flag_active_region
      logical, dimension(n_atoms) :: flag_initial_moment

      integer :: n_z_vectors

      character*120 :: info_str

      real*8 :: empty_states      ! this variable will count an upper estimate of the number of
                                  ! empty states needed in a calculation, for not too large smearing

      logical :: homogeneous_field_set = .false. !Determines whether a field has been requested
      logical :: image_potential_set = .false. !Determines whether a image-potential like 1/4z-potential has been requested

      ! for RRS-PBC scheme
      logical :: rrs_pbc_lv_found_before =.false.
      integer :: rrs_pbc_i_periodic

!  counters

      integer i_coord, i_species, i_atom, i_function

      integer :: i_spin, i_region

      integer i_periodic

      integer i_cube

      integer :: i_empty_atoms
      integer :: i_empty_atoms_2
      integer :: i_hess_block
      integer :: i_hess_block_lv
      integer :: i_hess_block_lv_atom

      integer, parameter :: max_n_ngh = 14
      integer :: a_ngh(3, max_n_ngh), n_ngh, i_ngh
      real*8 :: vec(3), rmin, BvK_lattvec(3, n_periodic)

!      logical :: flag_set_vacuum_level = .false. ! has the vacuum level been set?
!     Nadia
      integer i_multipole

!DB
      logical :: pp_found
      integer :: i_pp_atom
      integer :: i_pp_species
      character*20 :: pp_species_temp

      character(*), parameter :: func = 'read_geo'

      ! igor for the frozen core testing
      integer :: n_low_state

      ! for easier error handling
      ! should be updated to same algorithm as used in read_control.f90 one day.
      character*320 inputline

      ! for output of fractional coordinates
      real*8, dimension(3,n_atoms) :: frac_coords_temp
      real*8, dimension(3,n_atoms) :: coords_temp

      ! MOL: for Symmetry-constrained relaxation
      integer   :: SCR_i_atom, SCR_i_periodic
      integer   :: SCR_splitind0, SCR_splitind1, SCR_splitind2, SCR_splitind3, SCR_splitind4
      logical   :: SCR_lv_found_before

!  begin work

      empty_states = 0

      call localorb_info('',use_unit)
      call localorb_info( &
      "------------------------------------------------------------", &
           use_unit,'(A)' )
      call localorb_info( &
           "Reading geometry description geometry.in.", use_unit,'(10X,A)')
      call localorb_info( &
      "------------------------------------------------------------", &
           use_unit,'(A)' )

      if (n_periodic .ne. 0) then
        use_constraint = .false.
      end if
      ! If we are using any constraint at all, must initialize all arrays
      if (use_constraint) then

        !  assign every atom to region number 1 by default
        do i_atom=1,n_atoms
          constraint_region(i_atom)=1
        enddo

        flag_active_region = .true.

      end if

      if (use_relaxation_constraints) then
        ! initialize to no constraints first
        constrain_relaxation = .false.
        constrain_coords = .false.
        constrain_lv_relaxation = .false.
        constrain_lv_components = .false.
      end if

      if (spin_treatment.eq.1) then
        flag_initial_moment = .false.
      end if

!  initialize

      eof = .false.

      n_electrons = 0.d0

      i_periodic = 0

      i_atom = 0
      i_empty_atoms = 0

!    Nadia
      i_multipole = 0

      i_pp_atom = 0
      use_embedding_pp = .false.

      do i_species = 1, n_species, 1
        atoms_in_structure(i_species) = 0
      enddo


      initial_charge = 0.d0
      total_initial_charge = 0.d0

      i_hess_block           = 0
      i_hess_block_lv        = 0
      i_hess_block_lv_atom   = 0

      empty = .false.

! MOL: initailize for symmetry-constrained relaxation
      SCR_i_atom            = 0
      SCR_i_periodic        = 0
      SCR_lv_found_before   = .false.

! Allocate hess_blocks if needed
      if (relax_mode /= RELAX_OFF) then
         if (n_hess_blocks > 0) then
            call aims_allocate(hess_blocks,3,3,n_hess_blocks,"hess_blocks")
            call aims_allocate(hess_block_pos,2,n_hess_blocks,"hess_block_pos")
         end if

         if (n_hess_blocks_lv > 0) then
            call aims_allocate(hess_blocks_lv,3,3,n_hess_blocks_lv,"hess_blocks_lv")
            call aims_allocate(hess_block_pos_lv,2,n_hess_blocks_lv,"hess_block_pos_lv")
         end if

         if (n_hess_blocks_lv_atom > 0) then
            call aims_allocate(hess_blocks_lv_atom,3,3,n_hess_blocks_lv_atom,"hess_blocks_lv_atom")
            call aims_allocate(hess_block_pos_lv_atom,2,n_hess_blocks_lv_atom,"hess_block_pos_lv_atom")
         end if
      end if

!  Enforce external field if requested

      if (use_atom_ref) then

        homogeneous_field(3)=1.d-3

        if (myid.eq.0) then

          write (use_unit,'(2X,A,A)') &
          "Attention! Symmetry-breaking settings ", &
          "requested in control.in."

          write (use_unit,'(4X,A,3E14.6)') &
          "Homogeneous external field: Initial settings ", &
          homogeneous_field

          if (ext_safe.ne.'safe') then
            write (use_unit,'(2X,A,A)') &
            "Since no safe settings were requested, the actual ", &
            "field may be overwritten in geometry.in ."
          end if

        end if

        homogeneous_field=homogeneous_field/hartree_over_bohr

      end if

!  open input file

      open (8, FILE="geometry.in" )

!  read first line

      read (8,*,iostat=i_code) desc_str

      if (i_code.ne.0) then

         if (myid.eq.0) then
            write(use_unit,*) "Invalid file geometry.in, or file not found."
            write(use_unit,*) "Aborting."
         end if

         call aims_stop('Invalid file or file not found while reading geometry.in',func)

      end if

      do while (.not.eof)

        if (desc_str(1:1).eq."#") then

          continue

        elseif (desc_str.eq.'verbatim_writeout') then
           continue
        else if (desc_str.eq."lattice_vector") then
          backspace(8)

         i_periodic = i_periodic+1
         if (i_periodic .gt. 3) then
            if (myid.eq.0) then
              write(use_unit,*) "* You have specified four or more lattice vectors in geometry.in."
              write(use_unit,*) "* Unfortunately, FHI-aims only supports 3+1D electronic structure theory."
              write(use_unit,*) "* Please remove one or more lattice vectors from your geometry.in file.  Exiting."
           end if
           call aims_stop_coll('Too many lattice vectors found while reading geometry.in',func)
         end if
         lv_found_before =.true.

          read(8,*,err=99) desc_str, &
            (lattice_vector(i_coord,i_periodic), i_coord=1,3,1)

          do i_coord = 1,3,1
            lattice_vector(i_coord,i_periodic) = &
            lattice_vector(i_coord,i_periodic) / bohr
          enddo

        else if (desc_str.eq."rrs_pbc_lattice_vector") then
          backspace(8)

         rrs_pbc_i_periodic = rrs_pbc_i_periodic+1
         rrs_pbc_lv_found_before =.true.
         if (lv_found_before .and. rrs_pbc_lv_found_before) then
           write(use_unit,*) &
               "rrs-pbc lattice vector could not be specified with the normal lattice vector simultaneously"
           write(use_unit,*) "Aborting."
           call aims_stop('rrs-pbc lattice vectors conflict found while reading geometry.in',func)
         endif

          read(8,*,err=99) desc_str, &
            (rrs_pbc_lattice_vector(i_coord,rrs_pbc_i_periodic), i_coord=1,3,1)

          do i_coord = 1,3,1
            rrs_pbc_lattice_vector(i_coord,rrs_pbc_i_periodic) = &
            rrs_pbc_lattice_vector(i_coord,rrs_pbc_i_periodic) / bohr
          enddo

        else if (desc_str.eq."homogeneous_field" ) then

         homogeneous_field_set = .true.

         ! BL: Rotations and translations are now possible:
         ! Do not clean forces

         final_forces_cleaned = .false.

         ! Read this only if no "safe" external potential was
         ! requested earlier in control.in
         if (ext_safe.ne.'safe') then

           backspace(8)
           read(8,*,err=99) desc_str, &
            (homogeneous_field(i_coord), i_coord=1,3,1)

           write(info_str,'(2X,A,3E14.6)') &
            'Homogeneous electric field specified. Values:', &
               homogeneous_field
           call localorb_info(info_str,use_unit,'(A)')

           do i_coord = 1,3,1
             homogeneous_field(i_coord) = &
             homogeneous_field(i_coord)/hartree_over_bohr
           enddo


             if(n_periodic.eq.3 .and. &
               (dabs(homogeneous_field(1)) .gt. 1d-12 .or. &
                dabs(homogeneous_field(2)) .gt. 1d-12 )) then
                  write(info_str,*) &
                        '***Only fields in z-direction supported  &
                     &   for periodic boundary conditions'
                  call aims_stop(info_str,func)
             endif
             !this gradient
             if (use_dipole_correction.and.n_periodic.eq.3) then !
               write(info_str,*) &
                   'Dipole correction not possible when applying fields'
!               call aims_stop(info_str, func)
             endif



         else

           if (myid.eq.0) then
             write (use_unit, '(2X,A,A)') &
             "'safe' symmetry-breaking was explicitly requested ", &
             "in control.in."
             write (use_unit, '(2X,A,A)') &
             "Ignoring explicit setting of homogeneous_field."
           end if

         endif

         elseif (desc_str.eq."image_potential") then

             image_potential_set = .true.
             backspace(8)
             read(8,*,err=99) desc_str, image_plane_position, image_plane_cutoff, image_plane_scaling
             !Rescale to atomic units
             image_plane_position=image_plane_position/bohr
             image_plane_cutoff=image_plane_cutoff/bohr

             write(use_unit,'(2X,A)') '1/4(z-z0) (image-potential-like) potential has been requested'
             write(use_unit,*) 'z0=', image_plane_position, 'Bohr. z>', image_plane_cutoff+image_plane_position, 'Bohr.'
             write(use_unit,'(2X,A,F15.5)') 'The image potential will be scaled by a factor of ', image_plane_scaling
             write(use_unit,*) '*** WARNING: THIS IS AN EXPERIMENTAL FEATURE'
             write(use_unit,*) '*** PROCEED ONLY IF YOU KNOW WHAT YOU ARE DOING'
             write(use_unit,*) '*** AND ON YOUR OWN RISK'


        else if (desc_str.eq."data") then

         if (myid.eq.0) then
           write(use_unit,*) "Read data sub-tag without corresponding multipole tag"
           write(use_unit,*) "Aborting."
         endif
         call aims_stop('"data" subtag found without multipole tag while reading geometry.in',func)

        else if (desc_str.eq."multipole") then
          backspace(8)

          i_multipole = i_multipole+1

          read(8,*,err=99) desc_str, &
            (multipole_coords(i_coord,i_multipole), i_coord=1,3,1), &
            multipole_order(i_multipole), multipole_charge(i_multipole)

!
          if (multipole_order(i_multipole).gt.0) then
              call read_multipole_data &
                (i_multipole,multipole_order(i_multipole), &
                multipole_data(:,i_multipole),3)
          end if

!        transform multipole coordinates into atomic units (bohr)
          do i_coord = 1,3,1
            multipole_coords(i_coord,i_multipole) = &
            multipole_coords(i_coord,i_multipole)/bohr
          enddo

!
        else if (desc_str.eq."pseudocore") then
          backspace(8)
          i_pp_atom = i_pp_atom + 1

          read(8,*,err=99) desc_str, &
            (pp_coords(i_coord,i_pp_atom), i_coord=1,3,1), &
            pp_species_temp

!        transform pseudocore coordinates into atomic units (bohr)
          do i_coord = 1,3,1
            pp_coords(i_coord,i_pp_atom) = &
            pp_coords(i_coord,i_pp_atom)/bohr
          enddo

          pp_found = .false.

          do i_species=1, n_species, 1
             if (species_pseudoized(i_species)) then
                i_pp_species = species2pp_species(i_species)

                if (pp_species_temp.eq.species_name(i_species)) then
                   pp_species(i_pp_atom) = i_pp_species
                   pp_found = .true.
                end if
             end if
          enddo


          if (.not.pp_found) then
             if (myid.eq.0) then
                write(use_unit,'(1X,A,A,A,A,A)') &
                     "* Pseudoized species ", pp_species_temp, ", listed in ", &
                     "geometry.in, is not described in input file ", &
                     "control.in."
                write(use_unit,*) "* List of ", n_pp_species, " species: "
             end if
             do i_species=1, n_species, 1
                if (species_pseudoized(i_species)) then
                   i_pp_species = species2pp_species(i_species)
                   if (myid.eq.0) then
                      write(use_unit,*) "* ",i_pp_species,":",species_name(i_species)
                   end if
                end if
            enddo
            if (myid.eq.0) then
               write(use_unit,*) "* Aborting."
            end if

            call aims_stop('Non-existent pseudo species definition found while reading geometry.in',func)
          end if
          ! here we arrive at the final desicion whether to use pseudopotentials or not
          use_embedding_pp = .true.

        else if ((desc_str.eq."atom") .or. (desc_str.eq."atom_frac")) then
          ! ATTENTION: Please leave this keyword ONLY lower case - i.e., aims should NOT accept "ATOM" instead of "atom".
          ! The reason is that "ATOM" is a keyword used in .pdb files (CAPS only), and some tools, like jmol, use
          ! "ATOM" as a criterion to auto-identify the .pdb format.
          !  If we want jmol to autodetect our geometry.in files correctly, we must avoid "ATOM" as a valid keyword.

          ! VA: treat fractional coordinates "atom_frac" and cartesian corrdinates "atom"
          !     on equal footing. Here we just remember the TYPE of coordinate basis
          !     with the array coord_basis_atoms. An overall transformation -> cartesian coordinates
          !     is done after the structure is read with the call "all_atoms_to_cart"


          backspace(8)

          i_atom = i_atom+1

          !VA: remember the basis of atomic coordinates
          !     -- cartesian: 	 	0
          !     -- fractional 	 	1

          if      (desc_str .eq. "atom")      then
             coord_basis_atoms(i_atom) = 0
          else if (desc_str .eq. "atom_frac") then
             ! make sure that we are periodic:
             if (n_periodic /= 3) then
                write(info_str,'(1X,A)') &
                "* Fractional coordinates atom_frac make only sense in 3D periodic calculations!"
                call localorb_info(info_str,use_unit,'(A)')
                call aims_stop_coll ()
             else
                coord_basis_atoms(i_atom) = 1
             end if
          end if

          lv_found_before =.false.

          read(8,*,err=99) desc_str, (coords(i_coord,i_atom), i_coord=1,3,1), &
            species_temp

!  transform coordinates into atomic units (bohr)

          do i_coord = 1,3,1
            coords(i_coord,i_atom) = coords(i_coord,i_atom)/bohr
          enddo

!  check whether we know this species

          found = .false.

          do i_species=1, n_species, 1
            if (species_temp.eq.species_name(i_species)) then
              species(i_atom) = i_species
              found = .true.
            end if
          enddo

          if (.not.found) then
             if (myid.eq.0) then
                write(use_unit,'(1X,A,A,A,A,A)') &
                     "* Species ", species_temp, ", listed in ", &
                     "geometry.in, is not described in input file ", &
                     "control.in."
                write(use_unit,*) "* List of ", n_species, " species: "
             end if
             do i_species = 1,n_species,1
                if (myid.eq.0) then
                   write(use_unit,*) "* ",i_species,":",species_name(i_species)
                end if
            enddo
            if (myid.eq.0) then
               write(use_unit,*) "* Aborting."
            end if

            call aims_stop('Non-existent species definition found while reading geometry.in',func)
          end if

          atoms_in_structure(species(i_atom)) = &
                atoms_in_structure(species(i_atom))+1

!          if(species_pseudoized(species(i_atom))) then
!            n_electrons = n_electrons + pp_charge(species2pp_species(species(i_atom)))
!          else
            n_electrons = n_electrons + species_z(species(i_atom))
!          endif
          ! add up an estimate of the  number of empty states needed for reasonable (not too large) smearing
          empty_states = empty_states + 2. + real((maxval( atomic_l(species(i_atom),1:n_atomic(species(i_atom))) )+1)**2)

! constant force on atom

        else if (desc_str .eq. "external_force") then
           backspace(8)

           read(8,*,err=99) desc_str, (external_forces(i_coord,i_atom), i_coord=1,3,1)
           external_forces_on = .true.

           if (myid.eq.0) write(use_unit,'(2X,A,I5)')'i_atom : ',  i_atom

           !Transform external_force to atomic units (eV/AA -> Ha/bohr)
           external_forces(1,i_atom) = external_forces(1,i_atom) * bohr_over_hartree
           external_forces(2,i_atom) = external_forces(2,i_atom) * bohr_over_hartree
           external_forces(3,i_atom) = external_forces(3,i_atom) * bohr_over_hartree


! Ghost atom -> Counterpoise Correction

        else if (desc_str.eq."empty") then
          backspace(8)

          i_empty_atoms = i_empty_atoms+1
          i_empty_atoms_2 = i_empty_atoms+n_occ_atoms
        read(8,*,err=99) desc_str, (empty_coords(i_coord,i_empty_atoms), &
                   i_coord=1,3,1), &
                   species_temp

!  transform coordinates into atomic units (bohr)

          do i_coord = 1,3,1
       empty_coords(i_coord,i_empty_atoms) = &
            empty_coords(i_coord,i_empty_atoms)/bohr
       coords(i_coord,i_empty_atoms_2)= &
          empty_coords (i_coord,i_empty_atoms)
          enddo

! The coordinate basis is cartesian

       coord_basis_atoms(i_empty_atoms_2) = 0

!  this is a basis site only, no atom:
          empty(i_empty_atoms+n_occ_atoms) = .true.

!  check whether we know this species

          found = .false.

          do i_species=1, n_species, 1

            if (species_temp.eq.species_name(i_species)) then
              species(i_empty_atoms_2) = i_species
              found = .true.
            end if

          enddo

          if (.not.found) then

             if (myid.eq.0) then
                write(use_unit,'(1X,A,A,A,A,A)') &
                     "* Species ", species_temp, ", listed in ", &
                     "geometry.in, is not described in input file ", &
                     "control.in."
                write(use_unit,*) "* List of ", n_species, " species: "
             end if

             do i_species = 1,n_species,1

                if (myid.eq.0) then
                   write(use_unit,*) "* ",i_species,":",species_name(i_species)
                end if

            enddo

            if (myid.eq.0) then
               write(use_unit,*) "* Aborting."
            end if

            call aims_stop('Non-existent species definition found while reading geometry.in',func)

          end if

          atoms_in_structure(species(i_empty_atoms_2)) = &
           atoms_in_structure(species(i_empty_atoms_2))+1



        else if (desc_str.eq."initial_moment") then

          ! initial moment for spin polarisation

          if (spin_treatment.ne.1) then

            write(info_str,'(1X,A,A)') &
            "* Initial moment specified in geometry.in, ", &
            "but spin-unpolarized calculation requested."
            call localorb_info(info_str,use_unit,'(A)')

            write(info_str,'(1X,A,A)') &
            "* Ignoring initial moment."
            call localorb_info(info_str,use_unit,'(A)')

          else if (i_atom.eq.0) then

            write(info_str,'(1X,A,A)') &
            "* Initial moment specified in geometry.in, ", &
            "but no atoms defined as yet."
            call localorb_info(info_str,use_unit,'(A)')

            write(info_str,'(1X,A,A)') &
            "* Ignoring initial moment."
            call localorb_info(info_str,use_unit,'(A)')

          else

            backspace(8)

            read(8,*,err=99)desc_str, initial_moment(i_atom)

            if (abs(initial_moment(i_atom)) &
            .gt.(species_z(species(i_atom))) ) then

               write(info_str,'(1X,A,A,A,A)') &
               "* Default initial moment per atom is larger than ", &
               "number of electrons on species ", &
               species_name(i_species),"."
               call localorb_info(info_str,use_unit,'(A)')

               write(info_str,*) "This makes no sense - please correct."
               call localorb_info(info_str,use_unit,'(A)')

              call aims_stop('Invalid (default) initial moment encountered while reading geometry.in',func)
            end if

            flag_initial_moment(i_atom) = .true.
            kind_of_initial(i_atom) = 2

          end if

       else if (desc_str.eq."initial_charge") then
          ! initial charge for charged systems

          if (i_atom.eq.0) then

            write(info_str,'(1X,A,A)') &
            "* Initial charge specified in geometry.in, ", &
            "but no atoms defined as yet."
            call localorb_info(info_str,use_unit,'(A)')

            write(info_str,'(1X,A,A)') &
            "* Ignoring initial charge."
            call localorb_info(info_str,use_unit,'(A)')

          else

            backspace(8)
            read(8,*,err=99)desc_str, initial_charge(i_atom)
            total_initial_charge = total_initial_charge + initial_charge(i_atom)

            if (abs(initial_charge(i_atom)) &
            .gt.(species_z(species(i_atom))) ) then

               write(info_str,'(1X,A,A,A,A)') &
               "* Initial charge requested on an atom is larger than ", &
               "number of electrons on species ", &
               species_name(i_species),"."
               call localorb_info(info_str,use_unit,'(A)')

               write(info_str,*) "This makes no sense - please correct."
               call localorb_info(info_str,use_unit,'(A)')

              call aims_stop('Invalid initial charge encountered while reading geometry.in',func)
            else if ((n_periodic.gt.0).and.(initial_charge(i_atom).ne.0.d0)) then
               ! Unfortunately, cannot use free ions in the superposition of
               ! atomic densities initialization for periodic systems yet.

               write(info_str,'(1X,A,A)') &
               "* Warning - an initialization of periodic systems with charged free ions ", &
               "(initial_charge) is not yet implemented."
               call localorb_info(info_str,use_unit,'(A)')

               write(info_str,'(1X,A,A)') &
               "* The simple underlying issue is that an Ewald sum would also be needed ", &
               "in the initialization, but that is not yet implemented - sorry."
               call localorb_info(info_str,use_unit,'(A)')

              call aims_stop('Initial charge not possible in periodic systems - found while reading geometry.in',func)
            end if

          end if

       else if (desc_str.eq."velocity") then
          ! initial velocity for molecular dynamics run. v_half is the main integration array for velocities
          if (.not.use_molecular_dynamics) then
             write(info_str,'(2X,2A)') '* WARNING: Ignoring velocity ', &
                  'specification in geometry.in without MD request.'
             call localorb_info(info_str,use_unit,'(A)')
          else
            backspace(8)
            read(8,*,err=99) desc_str, v_half(1,i_atom), v_half(2,i_atom), &
               v_half(3,i_atom)
            v_half(:,i_atom) = v_half(:,i_atom)/bohr
          end if

       else if (desc_str.eq."calculate_friction") then
          !nonadiabatic friction is included for this atom
          if (.not.use_friction) then
            write(info_str,*) &
            "No friction requested - ignoring calculate_friction."
            call localorb_info(info_str,use_unit,'(A)')
          else if ((i_atom.eq.0).and.(i_periodic.eq.0)) then
            write(info_str,*) &
            "Found friction request for atom, but no atoms or ", &
            "lattice vectors specified yet."
            call localorb_info(info_str,use_unit,'(A)')

          else if (use_friction) then
            backspace(8)
            read(8,*,iostat=i_code) &
              desc_str, desc_str
            if (desc_str.eq.".true.") then
              friction_n_active_atoms = friction_n_active_atoms +1
              friction_atoms_list(i_atom) = .true.
              write(info_str,'(2X,A,I6,A)') &
              "Found friction request for atom ", i_atom
              call localorb_info(info_str,use_unit,'(A)')
            else if (desc_str.eq.".false.") then
              friction_atoms_list(i_atom) = .false.
            endif
          endif

       else if (desc_str.eq."constrain_relaxation") then
          ! this atom may be constrained during relaxation

          if ((.not.use_geo_relaxation).and. &
             (.not.use_molecular_dynamics)) then
            write(info_str,*) &
            "No relaxation requested - ignoring relaxation constraint."
            call localorb_info(info_str,use_unit,'(A)')
          else if ((i_atom.eq.0).and.(i_periodic.eq.0)) then
            write(info_str,*) &
            "Found relaxation constraint, but no atoms or lattice vectors specified yet."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,*) &
            "Ignoring constraint."
            call localorb_info(info_str,use_unit,'(A)')
          else if (use_relaxation_constraints) then

             ! VA: atoms constrained
             if (.not. lv_found_before) then

               ! this is the case of a valid relaxation constraint
               backspace(8)
               read(8,*,iostat=i_code) &
                 desc_str, desc_str
               if (desc_str.eq.".true.") then
                  constrain_relaxation(i_atom) = .true.
                  constrain_coords(1:3,i_atom) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for atom ", i_atom, &
                  ": All coordinates fixed."
                  call localorb_info(info_str,use_unit,'(A)')
               else if (desc_str.eq.".false.") then
                  constrain_relaxation(i_atom) = .false.
                  constrain_coords(1:3,i_atom) = .false.
               else if ( (desc_str.eq."x") .or. (desc_str.eq."X") ) then
                  constrain_relaxation(i_atom) = .false.
                  constrain_coords(1,i_atom) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for atom ", i_atom, &
                  ": x coordinate fixed."
                   call localorb_info(info_str,use_unit,'(A)')
               else if ( (desc_str.eq."y") .or. (desc_str.eq."Y") ) then
                  constrain_relaxation(i_atom) = .false.
                  constrain_coords(2,i_atom) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for atom ", i_atom, &
                  ": y coordinate fixed."
                  call localorb_info(info_str,use_unit,'(A)')
               else if ( (desc_str.eq."z") .or. (desc_str.eq."Z") ) then
                  constrain_relaxation(i_atom) = .false.
                  constrain_coords(3,i_atom) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for atom ", i_atom, &
                  ": z coordinate fixed."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               if ( ALL(constrain_coords(1:3,i_atom)) ) then
                  ! if all coordinates are fixed, we mark the atom as constrained
                  constrain_relaxation(i_atom) = .true.
               end if

             else if ((lv_found_before) .and. (relax_unit_cell .ne. 0)) then   ! lattice vectors are contrained

               ! this is the case of a valid relaxation constraint
               backspace(8)
               read(8,*,iostat=i_code) &
                 desc_str, desc_str
               if (desc_str.eq.".true.") then
                  constrain_lv_relaxation(i_periodic) = .true.
                  constrain_lv_components(1:3,i_periodic) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for lattice vector ", i_periodic, &
                  ": All components fixed."
                  call localorb_info(info_str,use_unit,'(A)')
               else if (desc_str.eq.".false.") then
                  constrain_lv_relaxation(i_periodic) = .false.
                  constrain_lv_components(1:3,i_periodic) = .false.
               else if ( (desc_str.eq."x") .or. (desc_str.eq."X") ) then
                  constrain_lv_relaxation(i_periodic) = .false.
                  constrain_lv_components(1,i_periodic) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for lattice vector ", i_periodic, &
                  ": x coordinate fixed."
                   call localorb_info(info_str,use_unit,'(A)')
               else if ( (desc_str.eq."y") .or. (desc_str.eq."Y") ) then
                  constrain_lv_relaxation(i_periodic) = .false.
                  constrain_lv_components(2,i_periodic) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for lattice vector ", i_periodic, &
                  ": y coordinate fixed."
                  call localorb_info(info_str,use_unit,'(A)')
               else if ( (desc_str.eq."z") .or. (desc_str.eq."Z") ) then
                  constrain_lv_relaxation(i_periodic) = .false.
                  constrain_lv_components(3,i_periodic) = .true.
                  write(info_str,'(2X,A,I6,A)') &
                  "Found relaxation constraint for lattice vector ", i_periodic, &
                  ": z coordinate fixed."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               if ( ALL(constrain_lv_components(1:3,i_periodic)) ) then
                  ! if all components are fixed, we mark the lattice vector as constrained
                  constrain_lv_relaxation(i_periodic) = .true.
               end if

            end if  ! lattice vectors found
          end if
        else if (desc_str.eq."esp_constraint") then
          !contraints for ESP charge fitting
          if (esp_constraint.eq.1)then
            backspace(8)
            read(8,*,err=99) desc_str, &
            esp_atom_constraint(1,i_atom), esp_atom_constraint(2,i_atom),&
            esp_atom_constraint(3,i_atom)
          elseif (esp_constraint.eq.2)then
            backspace(8)
            read(8,*,err=99) desc_str, &
            esp_atom_constraint(1,i_atom), esp_atom_constraint(2,i_atom)
          else
            write(info_str,*) &
            "No esp constrains (method 1/2) requested - ignoring constraint."
            call localorb_info(info_str,use_unit,'(A)')
          endif
        else if (desc_str.eq."constraint_region") then
          ! locally constrained DFT

          backspace(8)
          read(8,*,err=99)desc_str,region_temp

          if (.not.use_constraint) then

            write(info_str,'(1X,A,A)') &
            "* Constraint region specified in geometry.in, ", &
            "but no constraints defined in control.in."
            call localorb_info(info_str,use_unit,'(A)')

            write(info_str,'(1X,A,A)') &
            "* Ignoring constraint."
            call localorb_info(info_str,use_unit,'(A)')

          else if (i_atom.eq.0) then

            write(info_str,'(1X,A,A)') &
            "* Constraint region specified in geometry.in, ", &
            "but no atoms defined as yet."
            call localorb_info(info_str,use_unit,'(A)')

            write(info_str,'(1X,A,A)') &
            "* Ignoring constraint."
            call localorb_info(info_str,use_unit,'(A)')

          else

            ! check if we know the region label
            do i_region = 1, n_region, 1
              if (region_temp.eq.constraint_region_label(i_region)) then
                 flag_active_region(i_region) = .true.
                 constraint_region(i_atom) = i_region
                 exit
              end if
            enddo

            if (.not.flag_active_region(i_region)) then
              write(info_str,'(1X,A,I5,A,I5,A)') &
                "* Atom ", i_atom, " assigned to subsystem ", &
                region_temp," for locally constrained DFT - "
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(1X,A,A)') &
                "* but no such subsystem defined in control.in! ", &
                "Please correct."
              call localorb_info(info_str,use_unit,'(A)')
            end if

          end if

        else if (desc_str.eq."hessian_block") then
           if (relax_mode /= RELAX_OFF) then
              backspace(8)

              i_hess_block = i_hess_block+1

              if (i_hess_block > n_hess_blocks) then
                 call aims_stop('n_hess_blocks error',func)
              end if

              read(8,*,err=99) desc_str, &
              & hess_block_pos(1,i_hess_block), &
              & hess_block_pos(2,i_hess_block), &
              & hess_blocks(:,:,i_hess_block)

              hess_blocks(:,:,i_hess_block) &
              & = hess_blocks(:,:,i_hess_block)*bohr**2/hartree
           end if

        else if (desc_str.eq."hessian_block_lv") then
           if (relax_mode /= RELAX_OFF) then
              backspace(8)

              i_hess_block_lv = i_hess_block_lv+1

              if (i_hess_block_lv > n_hess_blocks_lv) then
                 call aims_stop('n_hess_blocks_lv error',func)
              end if

              read(8,*,err=99) desc_str, &
              & hess_block_pos_lv(1,i_hess_block_lv), &
              & hess_block_pos_lv(2,i_hess_block_lv), &
              & hess_blocks_lv(:,:,i_hess_block_lv)

              hess_blocks_lv(:,:,i_hess_block_lv) &
              & = hess_blocks_lv(:,:,i_hess_block_lv)*bohr**2/hartree
           end if

        else if (desc_str.eq."hessian_block_lv_atom") then
           if (relax_mode /= RELAX_OFF) then
              backspace(8)

              i_hess_block_lv_atom = i_hess_block_lv_atom+1

              if (i_hess_block_lv_atom > n_hess_blocks_lv_atom) then
                 call aims_stop('n_hess_blocks_lv_atom error',func)
              end if

              read(8,*,err=99) desc_str, &
              & hess_block_pos_lv_atom(1,i_hess_block_lv_atom), &
              & hess_block_pos_lv_atom(2,i_hess_block_lv_atom), &
              & hess_blocks_lv_atom(:,:,i_hess_block_lv_atom)

              hess_blocks_lv_atom(:,:,i_hess_block_lv_atom) &
              & = hess_blocks_lv_atom(:,:,i_hess_block_lv_atom)*bohr**2/hartree
           end if

        else if (desc_str.eq."trust_radius") then
           ! This keyword is only here because it needs to be part of
           ! the automatically written relaxation restart file, geometry.in.next_step .
           ! Conceptually, this keyword belongs into control.in, so please consider it a very rare exception.
           backspace(8)
           read(8,*,err=99) desc_str, trust_radius

           if (trust_radius.lt.0.d0) then
              write(info_str,'(1X,A)') &
                "* Invalid trust radius value found - the default (max_atomic_move) will be used later."
              call localorb_info(info_str,use_unit,'(A)')
              trust_radius = 0.d0
           else if (trust_radius.gt.(max_atomic_move*bohr)) then
              write(info_str,'(1X,A,F15.8,A)') &
                "* Large trust radius value found - resetting to the current value of max_atomic_move, ", &
                max_atomic_move*bohr, " A."
              call localorb_info(info_str,use_unit,'(A)')
              trust_radius = max_atomic_move
           else
              trust_radius = trust_radius/bohr
           end if

         elseif (desc_str.eq.'set_vacuum_level') then
            backspace(8)
            read(8,*,err=99) desc_str, desc_str
            if (desc_str=='auto') then
               call determine_vacuum_level()
            else
            backspace(8)
            read(8,*,err=99) desc_str, vacuum_z_level
               vacuum_z_level = vacuum_z_level/bohr
            endif

            flag_set_vacuum_level = .true.

            if (myid.eq.0) then
               write(use_unit,'(2X,A,E14.6)')'Vacuum level in z-axis : ',  vacuum_z_level*bohr
            end if

         else if (desc_str.eq.'magnetic_response') then
            if (i_atom < 1) &
                 & call aims_stop_coll('You have specified a magnetic reponse &
                 &to be calculated before any atoms. Exiting.', func)
            mr_atoms(i_atom) = .true.

         else if (desc_str.eq.'magnetic_moment') then
            backspace(8)
            if(i_atom < 1) &
                 & call aims_stop_coll('You have specified a magnetic &
                 &moment before any atoms. Exiting.', func)
            read(8,*,err=99) desc_str, mr_custom_mus(i_atom)

         else if (desc_str.eq.'nuclear_spin') then
            backspace(8)
            if (i_atom < 1) &
                 & call aims_stop_coll('You have specified a spin before any &
                 &atoms. Exiting.', func)
            read(8,*,err=99) desc_str, mr_custom_spins(i_atom)

         else if (desc_str.eq.'isotope') then
            backspace(8)
            if (i_atom < 1) &
                 & call aims_stop_coll('You have specified an isotope before &
                 &any atoms. Exiting.', func)
            read(8,*,err=99) desc_str, isotope_numbers(i_atom)

        else if (desc_str == "hessian_file") then
            write(info_str,"(2X,A)") "Found keyword 'hessian_file',"//&
               " indicating that this is the continuation of a structure"//&
               " optimization."
            call localorb_info(info_str,use_unit,"(A)")

            call check_hess_file()

            hess_in_file = .true.
        ! MOL: symmtry-constrained relaxation
        else if (desc_str.eq.'symmetry_n_params') then
            ! was already read in in dimensions.f90
            continue
        else if (desc_str.eq.'symmetry_params') then
            backspace(8)
            read(8,*,err=99) desc_str, (SCR_params(i_coord), i_coord=1,SCR_n_params,1)
            do i_coord=1,SCR_n_params
              SCR_params(i_coord) = adjustl(SCR_params(i_coord))
            end do
        else if (desc_str.eq.'symmetry_lv') then
            backspace(8)
            SCR_i_periodic = SCR_i_periodic+1
            if (SCR_i_periodic .gt. 3) then
               if (myid.eq.0) then
                 write(use_unit,*) "* Do not give more symmetry_lv than 3", &
                                   " (especially not more than lattice_vector)."
               end if
               call aims_stop_coll('Too many symmetry_lv tags found while reading geometry.in',func)
            end if
            SCR_lv_found_before =.true.
            read(8,'(A)') inputline
            SCR_splitind0 = SCAN(inputline, "symmetry_lv")
            SCR_splitind1 = SCAN(inputline(SCR_splitind0+1:),' ') + SCR_splitind0
            SCR_splitind2 = SCAN(inputline(SCR_splitind1+1:),',') + SCR_splitind1
            SCR_splitind3 = SCAN(inputline(SCR_splitind2+1:),',') + SCR_splitind2
            SCR_splitind4 = SCAN(inputline(SCR_splitind3+1:),'#') + SCR_splitind3
            SCR_lv_str(1,SCR_i_periodic) = adjustl(inputline(SCR_splitind1+1:SCR_splitind2-1))
            SCR_lv_str(2,SCR_i_periodic) = adjustl(inputline(SCR_splitind2+1:SCR_splitind3-1))
            if (SCR_splitind4 > SCR_splitind3) then
              SCR_lv_str(3,SCR_i_periodic) = adjustl(inputline(SCR_splitind3+1:SCR_splitind4-1))
            else
              SCR_lv_str(3,SCR_i_periodic) = adjustl(inputline(SCR_splitind3+1:))
            end if
        else if (desc_str.eq.'symmetry_frac') then
            backspace(8)
            SCR_i_atom = SCR_i_atom+1
            SCR_lv_found_before =.false.
            read(8,'(A)') inputline
            SCR_splitind0 = SCAN(inputline, "symmetry_lv")
            SCR_splitind1 = SCAN(inputline(SCR_splitind0+1:),' ') + SCR_splitind0
            SCR_splitind2 = SCAN(inputline(SCR_splitind1+1:),',') + SCR_splitind1
            SCR_splitind3 = SCAN(inputline(SCR_splitind2+1:),',') + SCR_splitind2
            SCR_splitind4 = SCAN(inputline(SCR_splitind3+1:),'#') + SCR_splitind3
            SCR_coord_str(1,SCR_i_atom) = adjustl(inputline(SCR_splitind1+1:SCR_splitind2-1))
            SCR_coord_str(2,SCR_i_atom) = adjustl(inputline(SCR_splitind2+1:SCR_splitind3-1))
            if (SCR_splitind4 > SCR_splitind3) then
              SCR_coord_str(3,SCR_i_atom) = adjustl(inputline(SCR_splitind3+1:SCR_splitind4-1))
            else
              SCR_coord_str(3,SCR_i_atom) = adjustl(inputline(SCR_splitind3+1:))
            end if
        ! MOL end block for reading symmetry-constrained relaxation
        else
          ! we do not know the specified label - abort run
              write(info_str,'(1X,A,A,A)') &
                "* Unknown descriptor in file geometry.in:", &
                desc_str," - please correct."
              call aims_stop_coll(info_str, func)
        end if

!  next line

        read (8,*,iostat=i_code) desc_str

        if (i_code.ne.0) then
          eof = .true.
        end if

      end do

!  close geometry.in

      close(8)
!  now, we have to evaluate which of the pseudoized atoms should be treated as i_atom

      i_atom=i_atom+i_empty_atoms
!  redundant consistency check, can be thrown out later
      if (i_atom.ne.n_atoms) then

         if (myid.eq.0) then
            write(use_unit,*) "* Internal consistency error: n_atoms changed "
            write(use_unit,*) "* between allocation and actual reading of "
            write(use_unit,*) "* geometry information"
            write(use_unit,*) "* n_atoms (before) = ", n_atoms
            write(use_unit,*) "* n_atoms (after) = ", i_atom
         end if

         call aims_stop('Atom number changed while reading geometry.in',func)
      end if

      !VA: in case of fractional coordinate input,
      !     get overall cartesian coordinates
      call all_atoms_to_cart ()

      ! Unit cell relaxation + cartesian constraints:
      if (relax_unit_cell > 0) then
         input_coords = coords
      end if

      call initialize_bravais_quantities()

      if (use_friction) then
        call friction_coordinate_workaround()
      endif

      if (out_vdwdf) call cube_cell_dimension( )

      ! PATCH - Disable empty sites in periodic boundary conditions for now
      if ((n_empty_atoms.gt.0).and.(n_periodic.gt.0)) then
         if (myid.eq.0) then
            write(use_unit,'(1X,A,A)') &
                 "* Input file geometry.in: 'empty' sites are presently ", &
                 " disabled "
            write(use_unit,'(1X,A,A)') &
                 "* for periodic boundary conditions. This is a known ", &
                 " limitation, and under construction."
            write(use_unit,'(1X,A)') &
                 "* Sorry for the inconvenience."
         end if
         call aims_stop &
         ('empty sites not supported for periodic boundary conditions - encountered while reading geometry.in',func)
      end if

      !  adjust number of electrons by overall system charge
      if (n_electrons.gt.charge) then
        n_electrons = n_electrons-charge
      else

         if (myid.eq.0) then
            write(use_unit,'(1X,A,E14.6,A)') &
                 "* The required charge of the system, ", charge, &
                 " electrons, "
            write(use_unit,'(1X,A,A,E14.6,A)') &
                 "* exceeds the total number of electrons of ", &
                 "all neutral atoms, ", &
                 n_electrons, " electrons."
            write(use_unit,'(1X,A,A)') &
                 "* Please adjust input file control.in ", &
                 "for a reasonable calculation."
         end if

      end if

      ! check consistency between charge of the system and initial_charge
      if (dabs(total_initial_charge - charge) .gt. 1d-5) then
         write(info_str,'(1X,A,A)') &
              "* Warning: Initial charge and total", &
              " charge of the system are inconsistent."
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A,F10.4,A,F10.4)') &
              "* Selected charge in control.in: ", charge, &
              "* Sum of initial charges in geometry.in: ", total_initial_charge
         call localorb_info(info_str,use_unit,'(A)')
      end if

!     Check consistency of relativity treatment
      species_z_max = 0
      do i_atom = 1, n_atoms, 1
        species_z_max = max(species_z(species(i_atom)),species_z_max)
      enddo
      if (flag_rel.eq.REL_none) then
        ! check whether _no_ relativistic treatment is justified
        !
        ! After learning a bit more about the cohesive energy of Cu, and since there is
        ! no harm in running with atomic_zora, relativity will now be enforced also for the 3d
        ! elements.
        !
        if (species_z_max.gt.20) then
          ! Elements beyond Ca should not run without relativity by default.
          write(info_str,'(1X,A)') &
            "* Warning! You are trying to treat an element beyond Ca (Z=20) without "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* accounting for relativity. This may lead to physically incorrect results."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* If you know what you are doing, you may override this warning using the "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* flag override_relativity in control.in. "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* Else, please use the keyword 'relativistic atomic_zora scalar' (for relaxation) or "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* 'relativistic zora scalar 1e-12' as described in the documentation."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* Additional information about a change made 04/05/2010 (will be removed 2011):"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* Since there is no harm in using atomic_zora, we now mandate relativity also"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* for the 3d elements. Apologies to all those whom this affects for the first time."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* Of course, any older nonrelativistic total energies for these elements must never"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* be compared (added or subtracted) directly to atomic_zora results."
          call localorb_info(info_str,use_unit,'(A)')
          if (.not.override_relativity) then
            call aims_stop_coll('', func)
          end if
        else if (species_z_max.gt.10) then
          ! Elements between Na and Ca will yield plausible results without relativity, but
          ! might as well be treated relativistically. We will put a warning here.
          write(info_str,'(1X,A)') &
            "* Attention. You are trying to treat an element beyond Ne (Z=10) without "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* accounting for relativity. This will lead to slight deviations from the true "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* physics of the system. Check your results if in doubt. "
          call localorb_info(info_str,use_unit,'(A)')
        end if
      else if (flag_rel.eq.REL_not_set) then
        if (species_z_max.le.20) then
          write(info_str,'(1X,A)') &
            "* Warning! You did not set the 'relativistic' flag in input file 'control.in'."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* Since you are not attempting to treat any elements beyond Ca(Z=20), FHI-aims"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* will default to 'relativistic none'. "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* You must ensure that this is really what you want!!"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* For example, DO NOT compute total-energy differences between results from"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* different settings of the 'relativistic' flag."
          call localorb_info(info_str,use_unit,'(A)')
          flag_rel = REL_none
        else
          write(info_str,'(1X,A)') &
            "* Warning! You did not set the 'relativistic' flag in input file 'control.in'."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* Since you are attempting to treat elements beyond Ca(Z=20), you must make "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* an explicit choice in input file 'control.in' - please READ the documentation. "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* For example, DO NOT compute total-energy differences between results from"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* different settings of the 'relativistic' flag."
          call localorb_info(info_str,use_unit,'(A)')
          call aims_stop_coll('', func)
        end if
      end if

!     regardless if spin-polarized or not, hund's rules
!     is default in any case
!     this is necessary for get_atomic_occ_number to work for
!     the non-polarized case as well (R.Gehrke)

!     prepare spin treatment if needed
      if (spin_treatment.eq.1) then

!       initialize atoms with default moment
        do i_atom = 1, n_atoms, 1
          if (.not.flag_initial_moment(i_atom)) then
             if (flag_moment) then
                initial_moment(i_atom) = default_initial_moment
                kind_of_initial(i_atom) = 2
             else
                ! apply hund's rules
                kind_of_initial(i_atom) = 1
                initial_moment(i_atom) = 0.d0
             end if
          end if
        enddo

      end if

      ! create lookuptables for charge initialization
      if (use_initial_rho) then
        call create_ini_type_lookup(initial_moment, initial_charge, &
           kind_of_initial, species, empty, type_charge, type_moment, &
           type_kind, type_species, atom_type)
      end if

      ! Condense and verify settings for locally constrained DFT
      if (use_constraint) then
        n_active_regions = 0
        do i_region = 1, n_region, 1
          if (flag_active_region(i_region)) then
            n_active_regions = n_active_regions+1
          end if
        enddo
        if (n_region.ge.n_active_regions) then
          ! condense region information to only that relevant for n_active_regions
          ! no more empty regions ...
          n_active_regions = 0
          do i_region = 1, n_region, 1

            if (flag_active_region(i_region) .and. &
                 (.not.use_hf_multiplicity)) then

              n_active_regions = n_active_regions+1

              constraint_region_label(n_active_regions) = &
                constraint_region_label(i_region)
              constraint_region_number &
                (constraint_region_label(n_active_regions)) = &
                n_active_regions

              do i_spin = 1, n_spin, 1
                constraint_electrons(n_active_regions, i_spin) = &
                constraint_electrons(i_region, i_spin)
              enddo

              do i_atom = 1, n_atoms, 1
                if (constraint_region(i_atom).eq.i_region) then
                  constraint_region(i_atom) = n_active_regions
                end if
              enddo
! continuing the wrapper for the multiplicity flag, calculating how many electrons
! to put in each spin channel
            else if (flag_active_region(i_region) .and. &
                             use_hf_multiplicity) then
! Determine if multiplicity makes sense
            if (MOD(hf_multiplicity, 2) .eq. 0 .and. &
                MOD(n_electrons, 2.d0) .eq. 0) then
             write (use_unit,'(2X,A)') "This multiplicity with this number ", &
             "of electrons works only with fractional occ. numbers. ", &
             "Please use constraint_spin flag."
              call aims_stop &
              ('Electron number and multiplicity incompatible - encountered while reading geometry.in',func)
            endif

              n_active_regions = n_active_regions+1

              constraint_region_label(n_active_regions) = &
                constraint_region_label(i_region)
              constraint_region_number &
                (constraint_region_label(n_active_regions)) = &
                n_active_regions

                constraint_electrons(n_active_regions, 1) = &
                  real(int((n_electrons+1.d0)/2.d0 + &
                  (hf_multiplicity - 1.d0)/2.d0))
                constraint_electrons(n_active_regions, 2) = &
                  real(int((n_electrons+1.d0)/2.d0 - &
                  (hf_multiplicity - 1.d0)/2.d0))
              do i_atom = 1, n_atoms, 1
                if (constraint_region(i_atom).eq.i_region) then
                  constraint_region(i_atom) = n_active_regions
                end if
              enddo
            end if !multiplicity, etc
          enddo
        end if ! n_region, n_active_region

        ! check if number of electrons corresponds to electrons in constraint
        electronsum=0.0d0
        do i_region=1,n_active_regions, 1
          do i_spin = 1, n_spin, 1
            electronsum=electronsum + &
              constraint_electrons(i_region,i_spin)
          enddo
        enddo

        if ( (abs(electronsum-n_electrons).ne.0.d0).and.(.not.force_occupation_projector) ) then
          write(info_str,'(1X,A,A)') &
            "* Total number of electrons in structure ", &
            "does not equal the number of electrons "
          call localorb_info(info_str,use_unit,'(A)')

          write(info_str,'(1X,A,A)') &
            "* specified in subsystems for locally constrained DFT:"
          call localorb_info(info_str,use_unit,'(A)')

          write(info_str,'(1X,A,F10.5)') &
            "* Total number of electrons from geometry.in: ",n_electrons
          call localorb_info(info_str,use_unit,'(A)')

          write(info_str,'(1X,A,A,F10.5)') &
            "* Number of electrons assigned to subsystems ", &
            "in control.in: ", electronsum
          call localorb_info(info_str,use_unit,'(A)')

          write(info_str,'(1X,A)') "* Please correct. "
          call localorb_info(info_str,use_unit,'(A)')

          call aims_stop &
          ('Electrons in geometry not consistent with request in control.in - encountered while reading geometry.in',func)
        endif

        if (use_scalapack.and.(n_active_regions.gt.1)) then

          write(info_str,'(1X,A)') &
            "* Error: Spin constraint with more than one region used together with ScaLapack."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* This combination is not yet functional. At this time, please use the Lapack"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* eigenvalue solver instead."
          call localorb_info(info_str,use_unit,'(A)')

          call aims_stop('Spin constraint with scalapack unsupported - encountered while reading geometry.in',func)

        else if (use_density_matrix.and.use_forces.and.(n_active_regions.gt.1)) then

          write(info_str,'(1X,A)') &
            "* Error: Spin constraint with more than one region used together with density."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* matrix and forces. This combination is not yet functional. "
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(1X,A)') &
            "* At this time, please use the Lapack eigenvalue solver instead."
          call localorb_info(info_str,use_unit,'(A)')

          call aims_stop &
          ('Spin constraint not supported with density matrix update - encountered while reading geometry.in',func)

        end if

      ! end treatment of locally constrained DFT
      end if

!  determine maximum numbers of atom-centered basis functions
!  (before basis reduction in shrink_*_basis)

      n_max_basis_fns = 0
      n_max_basis = 0
      do i_species = 1, n_species, 1
        mult = 1

        ! VB: mult is a multiplier that assists in finding the total number of
        !     different radial functions in the entire calculation. For
        !     radial functions that adapt their shape in the course of a calculation,
        !     atom-dependent rather than species-dependent arrays would have been needed.
        !     This feature war supported in an early version of the code, but support has
        !     long since ceased. I am keeping the "mult" option here anyway in case it
        !     ever comes back. The original usage was:
        !
        ! if (.not.adaptive) then
!       !   radial functions are stored per species
        !   mult = 1
        ! else
!       !   radial functions are stored per atom
        !   mult = atoms_in_structure(i_species)
        ! end if

        if (include_min_basis(i_species)) then
          n_max_basis_fns = n_max_basis_fns &
          + n_atomic(i_species) * mult
          do i_function = 1, n_atomic(i_species),1
            n_max_basis = n_max_basis &
            + (2*atomic_l(i_species,i_function)+1) &
            * atoms_in_structure(i_species)
          enddo
        end if
        n_max_basis_fns = n_max_basis_fns &
        + n_ionic(i_species) * mult
        do i_function = 1, n_ionic(i_species),1
          n_max_basis = n_max_basis &
          + (2*ionic_l(i_species,i_function)+1) &
          * atoms_in_structure(i_species)
        enddo


        n_max_basis_fns = n_max_basis_fns &
        + n_conf(i_species) * mult
        do i_function = 1, n_conf(i_species),1
          n_max_basis = n_max_basis &
          + (2*conf_l(i_species,i_function)+1) &
          * atoms_in_structure(i_species)
        enddo


        n_max_basis_fns = n_max_basis_fns &
        + n_hydro(i_species) * mult
        do i_function = 1, n_hydro(i_species),1
          n_max_basis = n_max_basis &
          + (2*hydro_l(i_species,i_function)+1) &
          * atoms_in_structure(i_species)
        enddo


        n_max_basis_fns = n_max_basis_fns &
        + n_gaussian(i_species) * mult
        do i_function = 1, n_gaussian(i_species),1
          n_max_basis = n_max_basis &
          + (2*gaussian_l(i_species,i_function)+1) &
          * atoms_in_structure(i_species)
        enddo

        n_max_basis_fns = n_max_basis_fns + n_sto(i_species)*mult
        do i_function = 1, n_sto(i_species)
           n_max_basis = n_max_basis + &
                & (2*sto_l(i_species,i_function)+1)* &
                & atoms_in_structure(i_species)
        end do

      enddo

      ! Why not use charge for this? In any case, if anyone uses this in control.in,
      ! there is a clear warning there.
      if (force_n_electrons) then
          n_electrons = forced_n_electrons  ! kind of brutal. poor code.
      end if

      ! Mostly output after here, in principle no more changes of key variables?

      write (info_str,'(2X,A)') "Input structure read successfully. "
      call localorb_info(info_str,use_unit,'(A)')

      if (.not. force_n_electrons) then
        write (info_str,'(2X,A,I8,A,A,F14.3,A)') &
            "The structure contains ", n_atoms, " atoms, ", &
            " and a total of ", n_electrons, " electrons."
      else
        write (info_str,'(2X,A,I8,A,A,F14.3,A)') &
            "The structure contains ", n_atoms, " atoms, ", &
            " and a FORCED number of ", n_electrons, " electrons."
      end if
      call localorb_info(info_str,use_unit,'(A)')

      if (n_multipoles.gt.0) then
          write (info_str,'(2X,A,I8,A)') "It furthermore comprises ", &
          n_multipoles, " multipoles."
          call localorb_info(info_str,use_unit,'(A)')
      endif
      call localorb_info(' ',use_unit,'(A)')

      if (use_molecular_dynamics.and.MD_initialconf_from_restart) then
        ! Check if restart file exists already here!
        ! We only check on process myid=0, as the restart file is only read there (for now).
        inquire(FILE=MD_restart_file,EXIST=MD_restart_file_exists)
        call bcast_logical(MD_restart_file_exists, 0)
      endif

      ! When provided by the user, output tags identifying this calculations as belonging
      ! to a larger set of calculations
      if (allocated(calculation_set)) then
        write(info_str, "(2X,A,1X,A)") "This calculation is part of set:     ", calculation_set
        call localorb_info(info_str)
        if (allocated(calculation_subset)) then
          write(info_str, "(2X,A,1X,A)") "This calculation is part of subset : ", calculation_subset
          call localorb_info(info_str)
        end if
        write(info_str, "(A)") ""
        call localorb_info(info_str)
      else if (allocated(calculation_subset)) then
        write(info_str, "(2X,A,1X,A)") "This calculation is part of subset : ", calculation_subset
        call localorb_info(info_str)
        write(info_str, "(A)") ""
        call localorb_info(info_str)
      end if

      if (use_molecular_dynamics.and.MD_initialconf_from_restart.and.MD_restart_file_exists) then
          write(info_str,'(2X,A)') "Input geometry read from MD restart information."
          call localorb_info(info_str,use_unit,'(A)')
      else
         write(info_str,'(2X,A)') "Input geometry:"
         call localorb_info(info_str,use_unit,'(A)')

         if (n_periodic.gt.0) then
            write(info_str,'(2X,A)') "| Unit cell: "
            call localorb_info(info_str,use_unit,'(A)')
            do i_periodic = 1, n_periodic, 1
               write(info_str,'(2X,A,3(2X,F16.8))') "|", &
                    (lattice_vector(i_coord,i_periodic)*bohr, &
                    i_coord=1,3,1)
               call localorb_info(info_str,use_unit,'(A)')
            enddo

         else
            write(info_str,'(2X,A)') "| No unit cell requested."
            call localorb_info(info_str,use_unit,'(A)')
         end if

         write(info_str,'(2X,A)') "| Atomic structure: "
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Atom ","x [A]","y [A]","z [A]"
         call localorb_info(info_str,use_unit,'(A)')
         do i_atom = 1, n_occ_atoms, 1
!            if (.not.species_pseudoized(species(i_atom))) then
            write(info_str,'(2X,A1,I5,A,A6,3(2X,F16.8))') &
                 "|",i_atom, ": Species ", &
                 species_name(species(i_atom)), &
                 (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
            call localorb_info(info_str,use_unit,'(A)')
!            end if
         enddo

        if (n_periodic.gt.0) then
          write(info_str,'(2X,A)')
          call localorb_info(info_str,use_unit,'(A)')
          call write_lattice_parameters(lattice_vector)
          write(info_str,'(2X,A)')
          call localorb_info(info_str,use_unit,'(A)')
        end if

        call xml_elem('coords', coords, unit='a.u.')
        if (n_periodic > 0) then
            call xml_elem('lattice_vector', lattice_vector, unit='a.u.')
        end if
        call xml_open('elems')
        do i_atom = 1, n_atoms
            call xml_elem('atom', species_name(species(i_atom)))
        end do
        call xml_close()

        if (n_multipoles.gt.0) then
           do i_multipole = 1, n_multipoles, 1
              write(info_str,'(2X,A1,I5,A,3(2X,F16.8),2X,I4,2X,F10.3)') &
               "|",i_multipole, ": Multipole ", &
               (multipole_coords(i_coord,i_multipole)*bohr, i_coord=1,3,1), &
                multipole_order(i_multipole), multipole_charge(i_multipole)
              call localorb_info(info_str,use_unit,'(A)')
           enddo
        endif


         if (use_embedding_pp) then
            do i_pp_atom = 1, n_pp_atoms, 1
                write(info_str,'(2X,A1,I5,A,A9,F16.8,2(2X,F16.8))') &
                  "|",i_pp_atom, ": Pseudocore ", &
                  species_name(pp_species2species(pp_species(i_pp_atom))), &
                  (pp_coords(i_coord,i_pp_atom)*bohr, i_coord=1,3,1)
                call localorb_info(info_str,use_unit,'(A)')
            enddo
         else
          ! in order to definitely set this flag false when no pseudocores are present
           use_nonlinear_core = .false.
         endif

         if (n_occ_atoms.lt.n_atoms) then
           write(info_str,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Ghost","x [A]","y [A]","z [A]"
           call localorb_info(info_str,use_unit,'(A)')
           do i_atom = n_occ_atoms+1, n_atoms, 1
              write(info_str,'(2X,A1,I5,A,A2,3(2X,F16.8))') &
                   "|",i_atom, ": Species ", &
                   species_name(species(i_atom)), &
                   (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
              call localorb_info(info_str,use_unit,'(A)')
       ! DB: 10/14/13
       ! empty sites with no basis sets should be excluded from relaxation.
       ! Actually, shouldn't we exclude all empty sites? Asignment of constraints
       ! to ghost atoms in the geometry.in does not work.
!              if(no_basis(species(i_atom))) then
                  constrain_relaxation(i_atom) = .true.
!              endif
       ! DB
           enddo
         endif

         call check_for_close_encounters(.true.)

         ! must call this after check_for_close_encounters,
         ! since spglib will crash if two atoms are in the same location.
         call perform_symmetry_analysis ()

         if (n_periodic > 0) then
            call output_bravais_quantities()

            ! This location, at the latest, is where both cell_volume and n_atoms are known.
            ! Set range separation parameter for Ewald summation here if automatic determination requested.
            ! Ewald_radius is the range separation parameter in bohr .
            if (Ewald_radius_automatic) then
               call estimate_Ewald_radius(n_atoms, cell_volume, Ewald_radius)
            end if
         end if

         coords_temp(:,:) = coords(:,:)
         if ((n_periodic>0) .and. (output_level /= 'MD_light' )) then
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


         call localorb_info('  ', use_unit, '(A)')

         if (use_initial_rho) then
           if (spin_treatment.eq.1) then
             write (info_str,'(2X,A)') "Initial moments and charges:"
             call localorb_info(info_str,use_unit,'(A)')
             write (info_str,'(2X,A)') &
             "|  Atom   Moment         Charge        Species"
             call localorb_info(info_str,use_unit,'(A)')
           else
             write (info_str,'(2X,A)') "Initial charges:"
             call localorb_info(info_str,use_unit,'(A)')
             write (info_str,'(2X,A)') &
             "|  Atom   Charge        Species "
             call localorb_info(info_str,use_unit,'(A)')
           end if
           do i_atom = 1, n_atoms, 1
             if (spin_treatment.eq.1) then
               if (kind_of_initial(i_atom) .eq. 2) then
                 write (info_str,'(2X,A,I5,1X,E14.6,1X,E14.6,3X,A2)') &
                 "| ", i_atom, &
                 initial_moment(i_atom), initial_charge(i_atom), &
                 species_name(species(i_atom))
                 call localorb_info(info_str,use_unit,'(A)')
               else
                 write (info_str,'(2X,A,I5,1X,A14,1X,E14.6,3X,A2)') &
                 "| ", i_atom, &
                 " (Hund's rule)", &
                 initial_charge(i_atom), species_name(species(i_atom))
                 call localorb_info(info_str,use_unit,'(A)')
               end if
             else
               write (info_str,'(2X,A,I5,1X,E14.6,3X,A2)') "| ", i_atom, &
               initial_charge(i_atom), species_name(species(i_atom))
               call localorb_info(info_str,use_unit,'(A)')
             end if
           enddo
           call localorb_info('  ', use_unit, '(A)')
         end if
      end if !not MD restart

      ! Next, adjust the number of states. Note that the system charge is already included
      ! in n_electrons at this point.
      if (calculate_all_eigenstates.and..not.use_full_spectrum) then
        if (n_empty.ne.-1) then
          write(info_str,'(2X,A)') &
               "You have specified both the number of empty states per atom and a method"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
               "| that uses calculate_all_eigenstates.  The latter setting (including all"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
               "| empty states possible) will override the former."
          call localorb_info(info_str,use_unit,'(A)')
        elseif(fo_finalise) then
           write(info_str,'(2X,A)') &
               "You have specified FODFT and a method that uses calculate_all_eigenstates."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
               "| This is not supported, out of sheer laziness, and implementing it"
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
               "| is likely a few lines of code."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
               "| Exiting."
          call localorb_info(info_str,use_unit,'(A)')
          call aims_stop('FODFT with calculate_all_eigenstates not supported - encountered while reading geometry.in',func)
        end if
        n_states = n_max_basis
      elseif ((n_empty.ne.-1).and.(fo_finalise)) then
        n_states = ceiling(real(n_electrons)/2.) + n_empty
        ! Adjust number of states for H(2n-1)@D^-A^- FODFT-mode
        ! necessary for electron transfer in di-anions
        if ((occ_fo).and.(fo_type.eq.'elec')) then
            n_states = n_states + 1
        endif
      elseif (n_empty.ne.-1) then ! in this case, the value per atom has been set in control.in
        ! account for spin degeneracy or individual spin channels
        n_states = ceiling(real(n_electrons)/2.) + ceiling(real(n_empty * n_atoms)/2.)
        if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)&
           n_states = n_electrons + ceiling(real(n_empty*n_atoms))
      else

        if (use_full_spectrum) then
            if (.not. flag_frozen_virtual_postSCF) then
                ! Note: For some methods (MP2, GW, RPA, ...) the number of KS states should be equal
                ! to the exact number of basis functions (Hamiltonian size).
                ! Depending on the details of the per-atom basis, n_states may later be revisited in
                ! subroutines shrink_fixed_basis_*()
                write(info_str,'(2X,A)') &
                  "Number of empty states per atom not set in control.in ."
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') &
                  "| Since you are using a method that relies on the unoccupied spectrum "
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') &
                  "| (MP2,GW,RPA et al.), will use the full Hamiltonian size (see below) "
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') &
                  "| as the max. possible number of states (occupied plus empty)."
                call localorb_info(info_str,use_unit,'(A)')
                n_states = n_max_basis
            else
                ! Igor Note: the frozen virtual scheme is well tested for MP2
                ! and RPA
                write(info_str,'(2X,A)') &
                  "Number of empty states per atom not set in control.in ."
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') &
                  "| Warning: although you are using a method that relies on the unoccupied spectrum "
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') &
                  "| (MP2,GW,RPA et al.), a frozen virtual scheme is required, where "
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A,I2,A)') &
                  "| ",n_frozen_virtual_states, " highest virtual orbitals are not taken into account."
                call localorb_info(info_str,use_unit,'(A)')
                n_states = n_max_basis - n_frozen_virtual_states
            endif
        else
        ! In this case, the number of empty states needed for the s.c.f. cycle is determined
        ! by default as counted by empty_states above.

          write(info_str,'(2X,A)') &
            "Number of empty states per atom not set in control.in - providing a guess from actual geometry."
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A,I8)') &
            "| Total number of empty states used during s.c.f. cycle: ", ceiling(empty_states/2.)
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,A)') &
            "If you use a very high smearing, use empty_states (per atom!) in control.in to increase this value."
          call localorb_info(info_str,use_unit,'(A)')

          ! adding more empty states for dielectric calculations, 50 more states
          ! for each atoms (TZ-start: tag1)
          if (flag_out_dielectric .or. flag_out_absorption) then
             empty_states = empty_states + 50.0*n_atoms
          endif
          ! (TZ-end: tag1)

          n_states = ceiling(real(n_electrons)/2.) + ceiling(empty_states/2.)
          if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)&
             n_states = n_electrons + ceiling(empty_states)

        end if

        write(info_str,'(2X)')
        call localorb_info(info_str,use_unit,'(A)')

      end if

      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(n_low_state)
      endif


      if (myid.eq.0) then
         write(use_unit,'(2X,A)') "Structure-dependent array size parameters: "
         if (use_initial_rho) then
           write (use_unit,'(2X,A,I8)') &
              "| Number of distinct atom types in initial rho : ", &
              n_ini_type
         end if
         write (use_unit,'(2X,A,I8)') &
              "| Maximum number of distinct radial functions  : ", &
              n_max_basis_fns
         write (use_unit,'(2X,A,I8)') &
              "| Maximum number of basis functions            : ", &
              n_max_basis
         write (use_unit,'(2X,A,I8)') &
              "| Number of Kohn-Sham states (occupied + empty): ", &
              n_states
         if (flag_frozen_core_postSCF) then ! count the frozen core states
            write (use_unit,'(2X,A,I8)') &
              "| Number of core electrons for post-SCF        : ", &
              (n_low_state - 1) * 2
         end if
         write(use_unit,'(A)') &
      "------------------------------------------------------------"

!     end if (myid.eq.0)
      end if

      if (use_hf_multiplicity) then
        fixed_spin_moment_electrons(1) = (n_electrons + hf_multiplicity - 1)/2.d0
        fixed_spin_moment_electrons(2) = (n_electrons - hf_multiplicity + 1)/2.d0
        ! write(info_str,*) "N_up and N_down used in calculation: ",fixed_spin_moment_electrons(1),fixed_spin_moment_electrons(2)
        ! call localorb_info(info_str,use_unit,'(A)')
      else if (fixed_spin_moment.and..not.use_hf_multiplicity) then
        fixed_spin_moment_electrons(1) = (n_electrons + spin_moment) / 2.d0
        fixed_spin_moment_electrons(2) = (n_electrons - spin_moment) / 2.d0
      end if


!     Unrelated sanity check:
!     n_states now known; in case that graphical output is required
!     for eigenstates ("cube" style), check now whether any unreasonable
!     eigenstates were requested.
      if (out_cube_soc) then
        ! SOC has twice as many states, needs its own conditional
        do i_cube = 1, n_cube, 1
          if (cube_type(i_cube).eq.'eigenstate' .or. &
                cube_type(i_cube).eq. 'eigenstate_density') then
            if (cube_index(i_cube).gt.2*n_states) then
               write(use_unit,*) ""
               write(use_unit,'(1X,A,A,A,I4,A)') &
                    "* Sanity check: Unreasonable ", &
                    "cube style output ", &
                    "requested: Eigenstate ", cube_index(i_cube), &
                    " does not exist for SOC system - please correct control.in."
               call aims_stop('Wrong cube eigenstate request detected while reading geometry.in',func)
            end if
          end if
        enddo
      else if (out_cube) then
        do i_cube = 1, n_cube, 1
          if (cube_type(i_cube).eq.'eigenstate' .or. &
                cube_type(i_cube).eq. 'eigenstate_density') then
            if (cube_index(i_cube).gt.n_states) then
               write(use_unit,*) ""
               write(use_unit,'(1X,A,A,A,I4,A)') &
                    "* Sanity check: Unreasonable ", &
                    "cube style output ", &
                    "requested: Eigenstate ", cube_index(i_cube), &
                    " does not exist - please correct control.in."
               call aims_stop('Wrong cube eigenstate detected while reading geometry.in',func)
            end if
          end if
        enddo
      end if

!     One more check: If a vacuum level is set, and used, for a surface calculation (work function calculation or dipole correction),
!     we must make sure that this vacuum level is where it belongs - i.e., in the vacuum. Currently, only a surface
!     perpendicular to z is supported

      if (calculate_work_function .or. use_dipole_correction) then

        ! check whether we have only one lattice vector with a z component
        n_z_vectors = 0
        if (n_periodic.gt.0) then
          do i_periodic = 1, n_periodic, 1
            if ( (lattice_vector(3,i_periodic).ne.0.d0) ) then
              ! this vector has a z component
              n_z_vectors = n_z_vectors+1
            end if
          enddo
        end if

        if (n_z_vectors.ne.1) then
           if (myid.eq.0) then
              write(use_unit,'(1X,A)') &
                "* Error: You have requested a work function calculation or dipole correction for a surface."
              write(use_unit,'(1X,A)') &
                "* This is only supported if there is a periodic surface, where the surface is in the xy plane"
              write(use_unit,'(1X,A)') &
                "* (i.e., only one of your lattice vectors can have a z component)."
           end if
           call aims_stop('Dipole correction clashes with lattice vector definition - found while reading geometry.in',func)
        end if

        ! Next, we need to check that all atoms are far enough away from the vacuum_level.
        ! This is only done after all atoms have been remapped to the central unit cell.
        ! See module pbc_lists.f90 / subroutine check_geometry_for_sanity

      end if

      if (use_cutCb) then
         if (flag_cutCb_width /= CUTCB_FROM_EXPLICIT) then
            cutCb_width = exp(cutCb_width_factor * sbtgrid_lnrange/sbtgrid_N)
            write(info_str,'(2X,A,F16.12,A)') &
            & 'Setting Coulomb cut-width to', cutCb_width, ' * rcut.'
            call localorb_info(info_str)
         end if
         if (flag_cutCb_rcut /= CUTCB_FROM_EXPLICIT) then
            if (n_periodic > 0) then
               do i_periodic = 1, n_periodic
                  BvK_lattvec(:, i_periodic) = lattice_vector(:, i_periodic) &
                  &                            * n_k_points_xyz(i_periodic)
               end do
               call get_neighboring_wigner_seitz(n_periodic, BvK_lattvec, &
               &                                 max_n_ngh, n_ngh, a_ngh)
               rmin = huge(rmin)
               do i_ngh = 1, n_ngh
                  vec = matmul(BvK_lattvec, dble(a_ngh(:, i_ngh)))
                  rmin = min(rmin, sqrt(sum(vec**2)))
               end do
               cutCb_rcut = cutCb_rcut_factor * 0.5d0 * rmin
            else
               cutCb_rcut = cutCb_rcut_factor * 8.d0/bohr
            end if
            write(info_str,'(2X,A,F16.12,A)') &
            & 'Setting Coulomb cutting radius rcut to', cutCb_rcut*bohr, ' A.'
            call localorb_info(info_str)
         end if
      end if


      ! VA: Unit cell relaxation and fixed fractional atoms:
      !     make sure that whole atom is fixed and not individual components

      if (relax_unit_cell /= 0 ) then
         if (use_relaxation_constraints) then
           do i_atom =1,n_atoms,1
             if ( ( any (constrain_coords(:,i_atom)))     .and. &
                 (.not.(constrain_relaxation (i_atom)))  .and. &
                 (coord_basis_atoms(i_atom) == 1)) then
                individual_atom_components_frac_fixed = .true.
             end if
           end do
           if (individual_atom_components_frac_fixed) then
             write(info_str,'(2X,A)') &
             & '* Cell Relaxation and constrained fractional coordinates: '
             call localorb_info(info_str)
             write(info_str,'(2X,A)') &
             & '* Currently only works for atoms which are entirely fixed. Fixed fractional components would lead to '
             call localorb_info(info_str)
             write(info_str,'(2X,A)') &
             & '* non-linear constraints which are not implemented yet! '
             call localorb_info(info_str)
             call aims_stop_coll ()
           end if
         end if
       end if


      ! VA: Check if relax_unit_cell shape and indiviual components of lattice vectors
      !     are constrained --> this is probably not intended. Give a warning

      if (relax_unit_cell == 2) then

         ! Check if there is more than one component fixed, but not all for
         ! a given lattice vector

         do i_periodic =1,n_periodic,1
            if ( (any (constrain_lv_components(1:3,i_periodic)))  .and. &
                 (.not. constrain_lv_relaxation (i_periodic)))  then
               individual_lv_components_fixed = .true.
            end if
         end do

         if (individual_lv_components_fixed) then
            write(info_str,'(2X,A)') &
            & '** WARNING: You requested unit cell relaxation with fixed angles between lattice vectors '
            call localorb_info(info_str)
            write(info_str,'(2X,A)') &
            & '**          AND fixed individual components of a lattice vector. This is probably not   '
            call localorb_info(info_str)
            write(info_str,'(2X,A)') &
            & '**          what you want to do. Continue with caution!'
            call localorb_info(info_str)
        end if

      end if

!Sanity checks for vacuum level
      if (n_periodic.eq.3 .and. .not. flag_set_vacuum_level &
           .and. homogeneous_field_set) then
          write(info_str,*) 'Electric fields in periodic systems require setting vacuum level'
          call localorb_info(info_str)
          write(info_str,*)  'Vacuum level will be determined automatically.'
          call localorb_info(info_str)
          call determine_vacuum_level()
          flag_set_vacuum_level = .true.
!          call aims_stop(info_str,func)
      endif

      if (.not. flag_set_vacuum_level .and. (use_dipole_correction.or.calculate_work_function))then
          write(info_str,*) 'Error: The dipole correction needs vacuum z-level to be defined'
          call localorb_info(info_str)
          write(info_str,*)  'Vacuum level will be determined automatically.'
          call localorb_info(info_str)
          call determine_vacuum_level()
          flag_set_vacuum_level = .true.
!            write(use_unit,*) 'Specify set_vacuum_level in geometry.in !'
!         call aims_stop_coll('', func)
      end if
!      if (.not. flag_set_vacuum_level .and. calculate_work_function)then
!         if (myid.eq.0) then
!            write(use_unit,*) 'Error: The work function calculation needs vacuum z-level to be defined'
!            write(use_unit,*) 'Specify set_vacuum_level in geometry.in !'
!         end if
!         call aims_stop_coll('', func)
!      end if

      if (image_potential_set.and..not.use_dipole_correction) then
          write(info_str,*) 'Error: Image potential requires ', &
                          & 'to set dipole correction. Please correct.'
           call aims_stop(info_str)
      endif
! MR: Ideally the wrapper should write pass the geometries first. For now I simply warn
!     that this geometry.in will be ignored
      if (use_pimd_wrapper) then
       write(info_str,'(A)') ''
       call localorb_info ( info_str )
       write(info_str,'(A)') '************************** W A R N I N G *******************************'
       call localorb_info ( info_str )
       write(info_str,'(A)') '* You are using the PIMD wrapper. Specifications and positions         *'
       call localorb_info ( info_str )
       write(info_str,'(A)') '* in geometry.in will be IGNORED - all is received from the wrapper.   *'
       call localorb_info ( info_str )
       write(info_str,'(A)') '* Please make sure species are declared in the same order in           *'
       call localorb_info ( info_str )
       write(info_str,'(A)') '* geometry.in and wrapper input.                                       *'
       call localorb_info ( info_str )
       write(info_str,'(A)') '************************************************************************'
       call localorb_info ( info_str )
       write(info_str,'(A)') ''
       call localorb_info ( info_str )
      endif

      ! TARP: Symmetry-constrained relaxation needs to have non-cell centered mapped coordinates,
      !       allocating and initializing here to ensure that
      if ( (use_geo_relaxation) .and. (use_symm_const_geo_relaxation) ) then
        call allocate_SCR()
        call SCR_create_coords_matrices()
        call SCR_create_lv_matrices()
        call cart2frac(lattice_vector, coords, frac_coords_temp)
        call SCR_symmetrize(frac_coords_temp)
        call SCR_symmetrize_lv(lattice_vector)
      end if

      return

   99 continue
      write(use_unit,*) 'Syntax error while reading file "geometry.in".'
      backspace(8)
      read(8,'(A)',iostat=i_code) inputline
      if(i_code<0) then
         call aims_stop_coll("Apparent end-of-file when re-reading line in 'geometry.in'?", func)
      end if
      if(i_code>0) then
         call aims_stop_coll("Unknown error re-reading line in 'geometry.in'...", func)
      endif
      write (use_unit,'(A)') "line: '"//trim(inputline)//"'"
      call aims_stop('Error encountered while re-reading line in geometry.in',func)

    end subroutine read_geo
!******
