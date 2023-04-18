!****s* FHI-aims/asym_inv_sqrt_of_auxmat_scalapack
!  NAME
!    asym_inv_sqrt_of_auxmat_scalapack
!  SYNOPSIS

subroutine asym_inv_sqrt_of_auxmat_scalapack(auxmat, name, transposed)

  !  PURPOSE
  !
  !    Get V^{-0.5} for a non-symmetric auxiliary (Coulomb) matrix
  !    (asymmetry due to grid integration) in the sense of
  !       V^- * (V+V^T)/2.
  !    where M^- is the generalized inverse of M.
  !
  !    The input matrix is assumed to be a 1d distributed version of
  !      V_{rs} = <r| (G^+ G) v |s>
  !    with product functions |r>, |s>, the grid-mapping operator G
  !    and the Coulomb operator v.  This matrix is to be applied to the
  !    ovlp_3fn matrix like:
  !      C_{abs} := \sum_r O_{abr} (V^- * (V+V^T)/2.)_{rs}
  !
  !  USES

  use dimensions
  use runtime_choices
  use prodbas
  use mpi_tasks
  use synchronize_mpi
  use scalapack_wrapper
  use localorb_io
  implicit none

  !  ARGUMENTS

  real*8, intent(INOUT) :: auxmat(n_basbas, n_loc_prodbas)
  character(*), intent(IN) :: name
  logical, intent(IN) :: transposed

  !  INPUTS
  !    o auxmat -- real array, e.g. bare coulomb matrix within
  !                the auxiliary basis
  !    o name -- name of matrix (for output only)
  !    o transposed -- if .true., transpose auxmat, therefore to be applied
  !                    to ovlp_3fn to the right (V*T instead of T*V).
  !  OUTPUTS
  !    o auxmat -- "square root" of input auxmat in the above sense
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  real*8, allocatable :: dist_auxmat(:,:)
  integer :: info
  character(*), parameter :: func = 'asym_inv_sqrt_of_auxmat_scalapack'

  ! Distribute to 2d (for efficiency)
  allocate(dist_auxmat(max_row_2d,max_col_2d), stat=info)
  call check_allocation(info, 'dist_auxmat', func)
  call dist_1d_2d(n_basbas, auxmat, ubound(auxmat,1), &
  &               dist_auxmat, ubound(dist_auxmat,1))

  ! Do work
  call asym_inv_sqrt_of_auxmat_scalapack_2d(dist_auxmat, name, transposed)

  ! Redistribute back to 1D
  call dist_2d_1d(n_basbas, dist_auxmat, ubound(dist_auxmat,1), &
  &               auxmat, ubound(auxmat,1))

  ! Tidy up
  deallocate(dist_auxmat)

end subroutine asym_inv_sqrt_of_auxmat_scalapack
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/asym_inv_sqrt_of_auxmat_scalapack_2d
!  NAME
!    asym_inv_sqrt_of_auxmat_scalapack_2d
!  SYNOPSIS

subroutine asym_inv_sqrt_of_auxmat_scalapack_2d(auxmat, name, transposed)

  !  PURPOSE
  !
  !    Same as asym_inv_sqrt_of_auxmat_scalapack() for already 2d distributed
  !    array.
  !
  !  USES

  use dimensions
  use runtime_choices
  use prodbas
  use mpi_tasks
  use synchronize_mpi
  use scalapack_wrapper
  use localorb_io
  implicit none

  !  ARGUMENTS

  real*8, intent(INOUT) :: auxmat(max_row_2d, max_col_2d)
  character(*), intent(IN) :: name
  logical, intent(IN) :: transposed

  !  INPUTS
  !    o auxmat -- real array, e.g. bare coulomb matrix within
  !                the auxiliary basis
  !    o name -- name of matrix (for output only)
  !    o transposed -- if .true., transpose auxmat, therefore to be applied
  !                    to ovlp_3fn to the right (V*T instead of T*V).
  !  OUTPUTS
  !    o auxmat -- "square root" of input auxmat in the above sense
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  real*8, allocatable :: aux_inverse(:,:), aux_sqrt(:,:)
  integer :: info
  character(*), parameter :: func = 'asym_inv_sqrt_of_auxmat_scalapack_2d'

  ! Distribute original matrix to 2d (for efficiency) and copy to both
  ! factors

  allocate(aux_inverse(max_row_2d,max_col_2d), stat=info)
  call check_allocation(info, 'aux_inverse', func)
  allocate(aux_sqrt(max_row_2d,max_col_2d), stat=info)
  call check_allocation(info, 'aux_inverse', func)
  aux_inverse = auxmat
  aux_sqrt = aux_inverse

  ! Get sqrt of 0.5*(V+V^T):
  call power_auxmat_scalapack_2d(aux_sqrt, 0.5d0, name)

  ! Get pseudo-inverse of V
  call gen_inv_auxmat_scalapack_2d(aux_inverse, name)

  ! Multiply
  if (transposed) then
     call localorb_info('Multiplying V^0.5 x V^+ (scalapack)', &
     &                  use_unit, '(2X,A)', OL_norm)
     call pdgemm('T', 'T', n_basbas, n_basbas, n_basbas, &
     &           1.d0, aux_sqrt,     1, 1, aux_sc_desc_2d, &
     &                 aux_inverse,  1, 1, aux_sc_desc_2d, &
     &           0.d0, auxmat,       1, 1, aux_sc_desc_2d)
  else
     call localorb_info('Multiplying V^+ x V^0.5 (scalapack)', &
     &                  use_unit, '(2X,A)', OL_norm)
     call pdgemm('N', 'N', n_basbas, n_basbas, n_basbas, &
     &           1.d0, aux_inverse,  1, 1, aux_sc_desc_2d, &
     &                 aux_sqrt,     1, 1, aux_sc_desc_2d, &
     &           0.d0, auxmat,       1, 1, aux_sc_desc_2d)
  end if
  deallocate(aux_inverse, aux_sqrt)

end subroutine asym_inv_sqrt_of_auxmat_scalapack_2d
!******
