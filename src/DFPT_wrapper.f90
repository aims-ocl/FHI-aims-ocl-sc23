
!****s* FHI-aims/DFPT_wrapper
!  NAME
!    DFPT_wrapper
!  SYNOPSIS

! For now this is a wrapper for the DFPT modules to keep the main file cleaner.
! Later this should move into a cleaned-up DFPT module itself.

    subroutine DFPT_wrapper(converged_cpscf,converged_scf)

    use applicable_citations
    use runtime_choices, only: restarting_scf
    use dimensions, only: use_DFPT, use_DFPT_reduce_memory, & 
                          use_DFPT_polarizability, use_DFPT_polar_reduce_memory, &
                          use_DFPT_dielectric, use_DFPT_phonon_gamma, use_DFPT_phonon, &
                          use_DFPT_phonon_reduce_memory, use_friction

    implicit none

    logical, intent(OUT) :: converged_cpscf
    logical, intent(IN) :: converged_scf

    if (.not. restarting_scf) then

      if(use_DFPT) then
         converged_cpscf=.false.
         call cpscf_solver( converged_cpscf )
      endif

      if(use_DFPT_reduce_memory) then
         converged_cpscf=.false.
         call cpscf_solver_reduce_memory( converged_cpscf )
      endif

      if(use_DFPT_polarizability) then
         converged_cpscf=.false.
         call cpscf_solver_polarizability( converged_cpscf )
      endif
 
      if(use_DFPT_polar_reduce_memory) then
         converged_cpscf=.false.
         call cpscf_solver_polar_reduce_memory( converged_cpscf )
      endif

      if(use_DFPT_dielectric) then
         converged_cpscf=.false.
         call cpscf_solver_dielectric( converged_cpscf )
      endif

      if(use_DFPT_phonon_gamma) then
         converged_cpscf=.false.
         call cpscf_solver_phonon_p0( converged_cpscf )
      endif

      if(use_DFPT_phonon) then
         converged_cpscf=.false.
         call cpscf_solver_phonon_p1( converged_cpscf )
      endif

      if(use_DFPT_phonon_reduce_memory) then
         converged_cpscf=.false.
         call cpscf_solver_phonon_reduce_memory( converged_cpscf )
      endif

      if(use_friction) then
         call calculate_nonadiabatic_friction(converged_scf)
         call cite_reference("friction")
      endif
    endif 

    end subroutine DFPT_wrapper
