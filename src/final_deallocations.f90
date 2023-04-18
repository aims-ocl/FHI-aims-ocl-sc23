!****s* FHI-aims/final_deallocations
!  NAME
!   final_deallocations
!  SYNOPSIS

  subroutine final_deallocations

!  PURPOSE
!  Collection of all cleanup routines from all modules, to make main.f look tidier ...
!
!  USES

  use physics
  use grids
  use species_data
  use free_atoms
  use basis
  use geometry
  use lapack_wrapper
  use scalapack_wrapper
  use spline
  use plot
  use mixing
  use constraint
  use cartesian_ylm
  use analytic_multipole_coefficients
  use hartree_potential_real_p0
  use pbc_lists
  use hartree_potential_recip
  use relaxation
  use wf_extrapolation
  use separate_core_states
  use precondition
  use SPE_solver
  use cg
  use molecular_dynamics
  use KH_core_states
  use vdw_correction
  use force_occupation, only: deallocate_force_occupation, &
                              deallocate_previous_eigenvector, &
                              deallocate_force_occupation_basis
  ! force_occupation_basis, force_occupation_projector (in dimensions.f90)
  use Hartree_F_p_functions   
  use ll_vdwdf
  use thermodynamic_integration
  use MD_QH_init
  use plus_u
  use esp_charges, only: deallocate_esp_out
  use pseudodata
  use friction
  use hartree_potential_storage
  use lpb_solver_utilities, only: deallocate_lpbsolver, solve_lpbe_only
  use mpb_solver_utilities, only: deallocate_mpbsolver, ln_cosh_function
  use runtime_choices, only: solvent_method, SOLVENT_MPE
  use mpe_interface, only: free_mpe_interface
  use spglib_symmetry, only: destroy_symmetry
  use mbd_std_wrapper, only: mbd_std_finalize
  use sym_base, only: destroy_symmetry_arrays
  use scalapack_soc, only : finalize_scalapack_soc
  use calculate_fock_matrix_p0, only : cleanup_fock_matrix_calculations
  use psi_at_nucleus_mod, only: cleanup_psi_at_nucleus
  use boys, only: deallocate_boys
  use elsi_wrapper, only: aims_elsi_finalize_scf
  use debugmanager, only: cleanup_debugmanager
  use applicable_citations, only: cleanup_citationlist
  use scgw_grid, only: reset_scgw_grid
  use load_balancing, only: reset_load_balancing
  use optics_module, only: reset_optics_module
  use localized_basbas, only: cleanup_localized_basbas
  use hartree_fock_p0, only: cleanup_hartree_fock_p0
  use hartree_fock, only: cleanup_hartree_fock
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
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

      call deallocate_grid_storage ( )

      call deallocate_physics ( )

      call cleanup_grids ( )
      
      call cleanup_species_data ( )

      call cleanup_free_atoms (  )

      call cleanup_species_waves (  )

      call cleanup_geometry ( )

      call cleanup_basis ( )

      ! AJL/April 2015
      ! Seems this is only called from initialize_prodbas previously?
      ! Adding it here as allocations in shrink_fixed_large_basis are
      ! never deallocated and we get misalignments in QM/MM.
      ! Please correct if a neater alternative exists
      call cleanup_ext_basis ( )

      call cleanup_lapack_wrapper ( )

      call cleanup_spline ( )

      call cleanup_plot ( )

      call cleanup_pulay ( )

      call cleanup_mixing ( )

      call cleanup_constraint ( )

      call cleanup_multipole_moments ( )

      call cleanup_cartesian_ylm ( )
      
      call cleanup_analytic_multipole_coefficients ( )

      call cleanup_hartree_potential_real ( )

      call deallocate_pbc_lists( )

      call deallocate_hartree_p_recip ( )

      call cleanup_relaxation ( )

      call cleanup_wf ( )

      call cleanup_relativity ( )

      call deallocate_separate_core_states ( )

      call deallocate_preconditioner ( )

      call deallocate_SPE_solver ( )

      call deallocate_lpbsolver ()
      
      call deallocate_mpbsolver ()

       ! CC: very ugly, see comments in function:
      ! Must revert to if call to avoid crashing 
      ! when skipping the SCF calls
      if (.not.skip_SCF) then
        call reset_grid_partitions ( )
      end if

      call cleanup_cg ( )

      call clean_MD ( ) 

      call cleanup_ll_vdw ( ) 

      ! call clean_kdens ( )

      call free_mpe_interface()

      call aims_elsi_finalize_scf()

      if (use_scalapack) then
         call finalize_scalapack()
      end if

      call reset_load_balancing()

      call reset_scgw_grid ()

      call reset_optics_module ()

      if(.not.use_scalapack.and.use_scalapack_DFPT_phonon) then
         call finalize_scalapack_DFPT_phonon()
      endif

      call deallocate_KH_core_states

      call deallocate_vdw

      if (use_mbd_std) then
          call mbd_std_finalize()
      end if

      if (use_libmbd) call mbd_calc%destroy()

      if (force_occupation_projector) then
        call deallocate_force_occupation()
        call deallocate_previous_eigenvector()
      end if

      if (force_occupation_basis) then
        call deallocate_force_occupation_basis()
      end if

      if (apply_boys_flag) then
        call deallocate_boys()
      end if

      if (use_plus_u) then
         call deallocate_plus_u
      endif

      call clean_F_p_functions ()

      call deallocate_TDI ()

      call deallocate_MD_QH ()

      call cleanup_hartree_fock_p0()

      call cleanup_hartree_fock()

      call cleanup_localized_basbas()

      call clear_lapack_overlap_matrix

      call cleanup_pseudodata

      call cleanup_hartree_potential_storage()

      call deallocate_esp_out()

      call deallocate_friction()

      if (use_spglib) then
        call destroy_symmetry ()
        !call destroy_symmats()
      endif

      if (use_symmetry_reduced_spg) then
        !call destroy_sym_maps ()
        call destroy_symmetry_arrays ()
      endif

      if (calculate_perturbative_soc.and.use_scalapack) then
        call finalize_scalapack_soc ()
      end if

      call cleanup_fock_matrix_calculations ()

      call cleanup_psi_at_nucleus()

      call cleanup_debugmanager()

  end subroutine final_deallocations
!******	
