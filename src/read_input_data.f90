!****s* FHI-aims/read_input_data
!  NAME
!   read_input_data
!  SYNOPSIS 
    subroutine read_input_data
! PURPOSE
!  High-level wrapper around all treatment of input data:
!  * Setup of all array dimensions for input data
!  * Allocation of all initially needed arrays
!  * Reading of all data from input file control.in
!  * Reading of all data from input file geometry.in
! USES 
      use physics
      use dimensions
      use runtime_choices
      use applicable_citations
      use thermodynamic_integration
      use MD_QH_init
      use MR_global, only: process_magnetic_keywords
      use pi_molecular_dynamics
      use geometry
      use pseudodata
      use plumed_new_interface
! INPUT
!  none
! OUTPUT
!  none
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
      implicit none
      integer*8 :: i_atom, i_bead, i_coord
      real*8 :: x

      ! local variables

      ! begin work
  
      !  obtain necessary dimensions and allocate all initially available arrays

      call obtain_initial_dimensions ( )


      !  now read all actual input data
      !   read global control variables from input file control.in
      !   this includes all variables in module species.f, and all variables in module runtime_choices.f

      call read_control ( )

      !  read geometry from input file geometry.in
      !  notice that n_electrons is determined explicitly here, as it enters the 
      !  general module physics.f
      call read_geo ( n_electrons )

      !  check whether there are pseudocores present and close too atoms
      !   if yes, count them as atoms
      !call prepare_pseudocoregrids()

      ! Process the keywords for magnetic response.
      if (magnetic_response) call process_magnetic_keywords()

      ! Read the data required for MD_QH_init from the respective files 
      ! that are specified in control.in
      if (use_MD_QH_init) then
         call MD_QH_initialize ( )
      end if

      ! Read initial molecular dynamics configuration if necessary;
      ! that might change the positions & velocities initialized from 
      ! the control files as the run might already be a long way from the 
      ! start
      if (use_molecular_dynamics.and.MD_initialconf_from_restart) then
!         if (MD_use_schedule) then
!           ! FIXME: Intelligent output about the initialization!!
!           MD_schedule_step   = 1 
!           MD_ensemble        = MD_schedule_ensemble      (MD_schedule_step)
!           MD_time            = MD_schedule_time          (MD_schedule_step) 
!           MD_temperature     = MD_schedule_temperature   (MD_schedule_step)
!           MD_tau_berendsen   = MD_schedule_tau_berendsen (MD_schedule_step)
!           MD_Q_NP            = MD_schedule_Q             (MD_schedule_step)
!           NVE_damping_factor = MD_schedule_damping_factor(MD_schedule_step)
!           MD_nu_andersen     = MD_schedule_nu_andersen   (MD_schedule_step)
!           MD_tau_BDP         = MD_schedule_tau_BDP       (MD_schedule_step)
!           MD_random_BDP      = MD_schedule_random_BDP    (MD_schedule_step)		
!          end if
         if (MD_restart_binary) then
           call read_MD_restart_binary ( )
         else
           call read_MD_restart ( )
         end if
      end if

      ! Read all data relevant for the thermodynamic integration, i.e.,
      ! the force constants -- either from input or from restart file
      if (use_thermodynamic_integration) then
         call read_TDI_QHA_file ( )
      else if (use_reffree_AS) then
         call init_AS ( )
      end if

      ! end work

      ! initialize plumed           
      if(plumed_new_plugin) then   
        call plumed_new_init ( )  
      endif

      ! PIMD restart
      if (use_pathint) then
        call allocate_PIMD ()
        if(PIMD_initialconf_from_restart) then 
          call read_PIMD_restart ( )
        else
          do i_atom = 1, n_atoms
            coords_beads(:,1,i_atom) = coords(:,i_atom)
          end do
          do i_bead = 2, n_beads
            do i_atom = 1, n_atoms
              do i_coord = 1, 3
                call random_number(x)
                coords_beads(i_coord, i_bead, i_atom) = coords(i_coord,i_atom) + 0.05 * x
              end do
            end do
          end do
        end if
      end if

      ! After ALL input was read, set applicable references that
      ! can be said to be used for certain in the present run
      if (use_prodbas) then
         ! Xinguo Ren et al. should please be cited for any use of
         ! the two-electron Coulomb operator in FHI-aims
         call cite_reference("Ren_NJP_2012")
      end if

      if (use_embedding_pp.or.use_embedding_potential) then
         ! Daniel Berger et al. should please be cited for any use of
         ! QM/MM embedding infrastructure or the use of pseudopotentials
         call cite_reference("Berger_JCP_2014")
      end if
      
      if (solvent_method.eq.SOLVENT_MPB) then
         ! Stefan Ringe et al. should please be cited for any use of
         ! the MPBE solver infrastructure
         call cite_reference("Ringe_JCTC_2016")
      end if

      if (use_analytical_stress) then
         call cite_reference("FHI-aims-analytical-stress-tensor")
      end if

      if (RI_type .eq. RI_LVL) then
         call cite_reference("FHI-aims-lvl")
      end if

      if (fo_dft) then
         call cite_reference("FODFT")
      end if
 
      if (use_DFPT.or.use_DFPT_reduce_memory.or.use_DFPT_phonon.or.use_DFPT_phonon_reduce_memory) then
         call cite_reference("dfpt")
      end if

      if (use_DFPT_polarizability.or.use_DFPT_polar_reduce_memory.or.use_DFPT_dielectric) then
         call cite_reference("dfpt_electricfield")
      end if

      if (use_plus_u) then
         call cite_reference("dftpU")
      endif

      ! In case of conflicts between keywords, the calculation should
      ! stop now.
      call check_consistency_of_keywords()

    end subroutine read_input_data
!******
