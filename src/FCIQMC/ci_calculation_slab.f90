!****s* FHI-aims/fciqmc/ci_calculation_slab
!  NAME
!   ci_calculation_slab
!  SYNOPSIS 

 subroutine ci_calculation_slab ( )

!  PURPOSE
!    Calculation of the correlation energy based on the configuration interaction algorithm
!    Frozen core approximation is available too. 
! USES

      use physics
      use mpi_tasks
      use runtime_choices
      use fciqmc_module

      implicit none

!  ARGUMENTS

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society.
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
    real*8 :: ddot
    ! variables
    real*8, dimension(:,:), allocatable :: inv_d_matrix
    real*8, dimension(:),   allocatable :: c_vect_prev
    real*8 :: dE_ci, dE_ci_prev, E_ci_prev
    real*8 :: E_corr_MP2, E_corr_MP3
    integer :: i_scf
    ! temp indices
    integer :: i_state, j_state, i_spin
    integer :: i_start
    integer :: a_state, b_state, i_ci
    integer :: errnum
    integer :: n_valence_a, n_valence_b
    !character*128 :: info_str

    if (myid.eq.0) then
        write(use_unit,'(2X,A)') "--------------------------------------------------------"
        write(use_unit,'(2X,A)') "|"
        write(use_unit,'(2X,A)') "| Configuration interaction (CI) calculation starts ... "
        write(use_unit,'(2X,A)') "|"
        write(use_unit,'(2X,A)') "--------------------------------------------------------"
    endif

    ! prepare two- and four-index integrals
    !call ci_initialization_slab()
    call get_4_index_integrals_slab()

    end subroutine ci_calculation_slab
!****** 

