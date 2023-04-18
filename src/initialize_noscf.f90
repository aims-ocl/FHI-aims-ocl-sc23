!****s* FHI-aims/initialize_noscf
!  NAME
!    initialize_scf
!  SYNOPSIS

    subroutine initialize_noscf ( converged )

!  PURPOSE
!  High-level wrapper to avoid former "first" scf loop
!
!  USES

      use localorb_io
      use dimensions
      use runtime_choices
      use physics
      use pbc_lists
      implicit none

!  ARGUMENTS

      logical :: converged

!  INPUTS
!    o none
!  OUTPUTS
!    o converged -- always .true. 
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
!

      character*100 :: info_str

      ! begin work
      write(info_str,'(A)') ''
      call localorb_info ( info_str )
      write(info_str,'(A)') '************************** W A R N I N G *******************************'
      call localorb_info ( info_str )
      if (use_pimd_wrapper) then
        write(info_str,'(A)') '* Skipping the SCF initialization for now - done inside wrapper      *'
        call localorb_info ( info_str )
      else
        write(info_str,'(A)') '* Skipping the SCF initialization.                                   *'
        call localorb_info ( info_str )
      endif
      write(info_str,'(A)') '************************************************************************'
      call localorb_info ( info_str )
      write(info_str,'(A)') ''
      call localorb_info ( info_str )
      ! FIXME: CC: One could dig deeper into this function calls to really just
      !            setup the required data
      ! MR: I am using this routine to skip the initiallization here, when using the code with pimd wrapper.
      !     probably this is not the cleanest way possible
      if (.not.use_pimd_wrapper) then
        call allocate_physics ( )
        call initialize_bc_dependent_lists()
        if( packed_matrix_format /= PM_none ) then 
           if(n_periodic .eq. 0) then
              call initialize_packed_matrix_cluster
           end if
        end if
        converged = .true.
      endif

    end subroutine initialize_noscf
!******
