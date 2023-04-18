!****s* FHI-aims/evaluate_KS_orbitals
!  NAME
!   evaluate_KS_orbitals
!  SYNOPSIS

      subroutine evaluate_KS_orbitals &
           ( n_points, wave, n_compute, i_basis, &
             KS_eigenvector, max_occ_number, &
             KS_orbital &
           )

!  PURPOSE
!  Subroutine evaluate_KS_orbitals computes KS orbitals for a set of
!  integration points.
!
!  USES

      use dimensions
      implicit none

!  ARGUMENTS

      integer :: n_points
      real*8, dimension(n_basis, n_points) :: wave

      integer :: n_compute
      integer :: i_basis(n_compute)

      real*8, dimension(n_basis, n_states) :: KS_eigenvector
      integer :: max_occ_number

      real*8, dimension(n_states, n_points) :: KS_orbital

!  INPUTS
!  o n_points -- number of grid points
!  o wave -- values of basis functions
!  o n_compute -- number of relevant basis functions
!  o i_basis -- list of relevant basis functions
!  o KS_eigenvector -- Kohn-Sham eigenvectors
!  o max_occ_number -- number of occupated states
!
!  OUTPUT
!  o KS_orbital -- Kohn-Sham orbitals
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



!  local variables

      real*8, dimension(n_compute, max_occ_number) :: &
           KS_ev_compute

!     counters

      integer i_state
      integer :: i_compute


!  begin work

!     FIXME: In the matrix-vector product below, we fail to take advantage of the benefits of
!            compute() . This is an obvious place for improvement, but not by a simple if statement.
!            Must condense wave() array itself already in evaluate_waves, to house only non-zero values, rest can be done via
!            index arrays.

! costs time ...      KS_orbital = 0.0d0

      ! must condense eigenvectors to only the needed ones ...
      do i_compute = 1, n_compute, 1
         do i_state = 1, max_occ_number, 1
            KS_ev_compute(i_compute, i_state) = &
            KS_eigenvector(i_basis(i_compute), i_state)
         enddo
      enddo

      call dgemm('T','N', max_occ_number, n_points, n_compute, 1.0d0, &
           KS_ev_compute, n_compute, wave, n_basis, &
           0.0d0, KS_orbital, n_states)

!  end work

      end subroutine evaluate_KS_orbitals
!---------------------------------------------------------------------
!******	
