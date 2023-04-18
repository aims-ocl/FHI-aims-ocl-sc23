!****s* FHI-aims/obtain_initial_dimensions
!  NAME
!   obtain_initial_dimensions
!  SYNOPSIS

      subroutine obtain_initial_dimensions &
        ( )

!  PURPOSE
!  Subroutine obtain_initial_dimensions is a wrapper about all
!  initial allocation work.
!
!  USES
        use dimensions
        use grids
        use geometry
        use species_data
        use mixing
        use plot
        use constraint
        use relaxation
        use mpi_tasks
        use localorb_io
        use vdw_correction
        use ll_vdwdf
        use molecular_dynamics
        use boys,                       only : allocate_boys
        use force_occupation, only: allocate_force_occupation, &
                                    allocate_force_occupation_basis 
        ! force_occupation_basis, force_occupation_projector (in
        ! dimensions.f90)
        use runtime_choices
        use thermodynamic_integration
        use MD_QH_init
        use pseudodata
        use friction
        use esp_charges,                only : allocate_esp_out
      implicit none

!  ARGUMENTS
!  none
!
!  INPUTS
!  none
!  OUTPUTS
!  files
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals:
!    FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject
!   to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!---------------------------------------------------------------------
      CHARACTER(LEN=1024) :: info_str

!  begin work


      call localorb_info('')
        call localorb_info( &
             "Obtaining array dimensions for all initial allocations:", &
             use_unit,'(2X,A)')

!       parse control.in only for necessary dynamic array dimensions.
        call parse_control &
        ( )

!       parse geometry.in only for necessary dynamic array dimensions
        call parse_geometry &
        ( )

!       write out the basic array dimensions

        if (myid.eq.0) then
           write(info_str,*)
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A)') "Basic array size parameters: "
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Number of species                 : ",n_species
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Number of atoms                   : ",n_atoms
           call localorb_info(info_str, use_unit)

           if (n_periodic.gt.0) then
              write(info_str,'(2X,A,I8)') &
                   "| Number of lattice vectors         : ",n_periodic
              call localorb_info(info_str, use_unit)
           end if
           write(info_str,'(2X,A,I8)') &
                "| Max. basis fn. angular momentum   : ",l_wave_max
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Max. atomic/ionic basis occupied n: ",n_wave_max
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Max. number of basis fn. types    : ",n_basis_types
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Max. radial fns per species/type  : ",n_max_ind_fns
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Max. logarithmic grid size        : ",n_max_grid
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Max. radial integration grid size : ",n_max_radial
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Max. angular integration grid size: ",n_max_angular
           call localorb_info(info_str, use_unit)


           if (use_angular_division) then
              write(info_str,'(2X,A,I8)') &
                   "| Max. angular grid division number : ", &
                   n_max_angular_division
              call localorb_info(info_str, use_unit)
           end if
           write(info_str,'(2X,A,I8)') &
                "| Radial grid for Hartree potential : ",n_hartree_grid
           call localorb_info(info_str, use_unit)

           write(info_str,'(2X,A,I8)') &
                "| Number of spin channels           : ",n_spin
           call localorb_info(info_str, use_unit)

           if (use_embedding_potential) then
              write(info_str,'(2X,A,I8)') &
                   "| External embedding multipoles     : ",n_multipoles
              call localorb_info(info_str, use_unit)
           end if
           if (n_pp_atoms.gt.0) then
              write(info_str,'(2X,A,I8)') &
                   "| Number of pseudoized species      : ",n_pp_species
              call localorb_info(info_str, use_unit)
              
              write(info_str,'(2X,A,I8)') &
                   "| Number of pseudocores             : ",n_pp_atoms
              call localorb_info(info_str, use_unit)
           end if
           if (use_vdw_correction) then
              write(info_str,'(2X,A,I8)') &
                   "| Van der Waals pairs for empirical correction  : ", &
                     vdw_pairs
              call localorb_info(info_str, use_unit)
           end if
           if (use_constraint) then
              write(info_str,'(2X,A,I8)') &
                   "| Locally constrained regions       : ", n_region
              call localorb_info(info_str, use_unit)
           end if
           if (use_cube_output) then
              write(info_str,'(2X,A,I8)') &
                   "| Cube type output requested        : ", n_cube
              call localorb_info(info_str, use_unit)
           end if
           if (use_plus_u) then
              write(info_str,'(2X,2A,I8)') &
                 "| Max. number of (n,l) shells per", &
                 " species with '+U' treatment: ", n_max_shells_plus_u
              call localorb_info(info_str, use_unit)
           end if
        end if

!       allocate arrays for all grids
!        print*,"allocate_grids"
        call allocate_grids &
        ( )

!       allocate arrays for species information
!        print*,"allocate_species_data"
        call allocate_species_data &
        ( )
!       allocate arrays needed for pseudopot infrastructure
!        if(use_embedding_pp) then
!           print*,"allocate_psp_arrays"
           call allocate_pseudospecies_arrays ( )
!        end if

!       allocate arrays for geometry
!        print*,"allocate_geometry"
        call allocate_geometry &
        ( )

!       allocate mixing parameters
!        print *,"allocate_mixing"
        call allocate_mixing &
        ( )

        ! allocate relaxation constraints, which might also become eligible for MD
!        print *,"allocate_relaxation"
        call allocate_relaxation ( )

        call allocate_friction( )

        if (use_molecular_dynamics) then
           call allocate_MD ( )
        end if

        if (use_thermodynamic_integration .or. use_reffree_AS) then
           call allocate_TDI ( )
        end if

        if (use_MD_QH_init) then
           call allocate_MD_QH ( )
        end if

        if (use_constraint) then
           call allocate_constraint ( )
        end if

!       if output for plotting requested, allocate here
        if (use_cube_output) then
          call allocate_plot &
          ( )
        end if

!       if output for esp charges requested, allocate here
        if (use_esp) then
          call allocate_esp_out &
          ( )
        end if

!       Allocate arrays for vdw correction
        call allocate_vdw ( )

!       Allocate arrays for ll_vdw_functional correction
        if (use_ll_vdwdf) then
          call allocate_ll_vdw &
          ( )
        end if

!       Allocate arrays for forced occupation 
        if (force_occupation_projector) then
          call allocate_force_occupation
        end if

        if (force_occupation_basis) then
          call allocate_force_occupation_basis
        end if

        if (apply_boys_flag) then
          call allocate_boys
        end if


        end subroutine obtain_initial_dimensions
!---------------------------------------------------------------------
!******
