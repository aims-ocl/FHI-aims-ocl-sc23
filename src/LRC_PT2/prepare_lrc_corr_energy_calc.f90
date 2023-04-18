!****s* FHI-aims/prepare_lrc_corr_energy_calc
!  NAME
!   prepare_lrc_corr_energy_calc
!  SYNOPSIS

      subroutine prepare_lrc_corr_energy_calc()

!  PURPOSE
!  This subroutine evaluate the necessary matrix elements (3-center overlap,
!  coulomb matrix) used for accurate correlation energy calculations beyond DFT
!  (e.g. MP2, RPA, RPA+, etc).
!
!  USES

      use dimensions
      use species_data
      use runtime_choices
      use prodbas
      use hartree_fock
      use timing
      use mpi_tasks
      use sbt_overlap_aims
      use localorb_io, only: OL_norm
      implicit none

!  ARGUMENTS

!  INPUT
!    none
!  OUTPUT
!    none
!
!  AUTHOR
!    Igor Ying Zhang, FHI, December 2016
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
      character*20 filename
      logical :: need_o3fn
      real*8 :: time_prodbas_add, clock_time_prodbas_add

! Counter
      integer i_basis_1

      integer :: info
      character*150 :: info_str
      character(*), parameter :: func = 'prepare_lrc_corr_energy_calc'

      call get_timestamps(time_prodbas_add, clock_time_prodbas_add)

      ! --- Always recalculate ovlp_3fn for lrc-PT2 ---
      ! At present, it supports RI_type=RI_V only

      lrc_pt2_started = .true.
      hse_omega_hf = lrc_pt2_omega
      use_logsbt_for_radial_hse_integration = .true.
      !write(use_unit,*) &
      !  "igor debug use_logsbt_for_radial_hse_integration", use_logsbt_for_radial_hse_integration

      call cleanup_hartree_fock()

      !call initialize_prodbas()
      
      call initialize_hartree_fock()

      ! --- Timing

      call get_times(time_prodbas_add, clock_time_prodbas_add, &
      &              time_prodbas_total, clock_time_prodbas_total)

      call output_timeheader('2X', 'End of correlation preparation',&
                             OL_norm)
      call ovlp3fn_timings('2X')

      end subroutine prepare_lrc_corr_energy_calc
!***************
