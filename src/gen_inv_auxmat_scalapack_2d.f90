!****s* FHI-aims/gen_inv_auxmat_scalapack_2d
!  NAME
!    gen_inv_auxmat_scalapack_2d
!  SYNOPSIS

subroutine gen_inv_auxmat_scalapack_2d(auxmat, name)

  !  PURPOSE
  !    Get the generalized inverse of auxmat
  !  USES

  use dimensions
  use runtime_choices
  use prodbas
  use mpi_tasks
  use synchronize_mpi
  use scalapack_wrapper
  use scalapack_utils
  use localorb_io
  implicit none

  !  ARGUMENTS

  real*8, intent(INOUT) :: auxmat(max_row_2d,max_col_2d)
  character(*), intent(IN) :: name

  !  INPUTS
  !    o auxmat -- real array, e.g. bare coulomb matrix within
  !                the auxiliary basis
  !    o name -- name of matrix (for output only)
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

  real*8, allocatable :: S(:), U(:,:), VT(:,:)
  real*8, allocatable :: work(:)
  real*8 :: work_tmp(1)
  integer :: lwork, i_sv, i_col
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'gen_inv_auxmat_scalapack_2d'

  if (name /= '') then
     write(info_str, "(A,' ',A,' ',A,' (2d-scalapack)')") &
     & 'Getting pseudo-inverse of', trim(name), 'matrix'
     call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  end if

  allocate(S(n_basbas), stat=info)
  call check_allocation(info, 'S', func)
  allocate(U(max_row_2d, max_col_2d), stat=info)
  call check_allocation(info, 'U', func)
  allocate(VT(max_row_2d, max_col_2d), stat=info)
  call check_allocation(info, 'VT', func)

  ! Workspace query
  lwork = -1
  call pdgesvd('V', 'V', n_basbas, n_basbas, &
  &            auxmat, 1, 1, aux_sc_desc_2d, &
  &            S, &
  &            U, 1, 1, aux_sc_desc_2d, &
  &            VT, 1, 1, aux_sc_desc_2d, &
  &            work_tmp, lwork, info)
  if (info /= 0) call aims_stop('Workspace query failed for pdgesvd', func)
  lwork = nint(work_tmp(1))
  allocate(work(lwork), stat=info)
  call check_allocation(info, 'work', func)

  ! Singular value decomposition: auxmat = U * diag(S) * VT
  call pdgesvd('V', 'V', n_basbas, n_basbas, &
  &            auxmat, 1, 1, aux_sc_desc_2d, &
  &            S, &
  &            U, 1, 1, aux_sc_desc_2d, &
  &            VT, 1, 1, aux_sc_desc_2d, &
  &            work, lwork, info)

  ! (Pseudo-)invert S
  do i_sv = 1, n_basbas
     if (use_smooth_prodbas_threshold) then
        S(i_sv) = S(i_sv) / (S(i_sv)**2 + prodbas_threshold**2)
     else if (S(i_sv) > prodbas_threshold) then
        S(i_sv) = 1.d0 / S(i_sv)
     else
        S(i_sv) = 0.d0
     end if
  end do

  ! Up = U * diag(inv_S)
  do i_sv = 1, n_basbas
     i_col = sclpck_loc_ind(i_sv, mypcol_aux_2d, npcol_aux_2d, nb_aux_2d, 0)
     if (i_col > 0) then
        U(:, i_col) = U(:, i_col) * S(i_sv)
     end if
  end do
  
  ! inv_auxmat = VT^T * diag(inv_S) * U^T = VT^T * Up^T
  call pdgemm('T', 'T', n_basbas, n_basbas, n_basbas, &
  &           1.d0, VT,     1, 1, aux_sc_desc_2d, &
  &                 U,      1, 1, aux_sc_desc_2d, &
  &           0.d0, auxmat, 1, 1, aux_sc_desc_2d)

  deallocate(S, U, VT, work)

end subroutine gen_inv_auxmat_scalapack_2d
!******
