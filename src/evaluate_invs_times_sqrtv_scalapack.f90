!****s* FHI-aims/evaluate_invs_times_sqrtv_scalapack
!  NAME
!    evaluate_invs_times_sqrtv_scalapack
!  SYNOPSIS

subroutine evaluate_invs_times_sqrtv_scalapack(ovlp_prodbas, coulomb_matr)

  !  PURPOSE
  !
  !   calculate the inverse of overlap matrix S times square root of the
  !   Coulomb matrix V, namely S^(-1)*V^(1/2).
  !   This is the scalapack version, 
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

  real*8  coulomb_matr(n_basbas,n_loc_prodbas)
  real*8  ovlp_prodbas(n_basbas,n_loc_prodbas)

  !  INPUTS/OUTPTS
  !    o ovlp_prodbas -- on entry, stores the overlap matrix with the
  !                                auxiliary basis
  !                      on exit, the inverse of the overlap
  !    o coulomb_matr -- on entry, stores the Coulomb interaction matrix
  !                      on exit, stores the inverted overlap matrix times the
  !                               square root of the Coulomb matrix
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

  integer :: info
      integer   sc_desc_am(dlen_)
  real*8, dimension(:,:), allocatable ::  invs_times_sqrtv(:,:)

  ! Get V^0.5 and S^-1

  call power_auxmat_scalapack(coulomb_matr, 0.5d0, 'Coulomb')
  call power_auxmat_scalapack(ovlp_prodbas, -1.d0, 'overlap')

  call localorb_info('Multiplying S^-1 x V^0.5 (scalapack)', &
  &                  use_unit, '(2X,A)', OL_norm)

  ! Multiply

  call descinit(sc_desc_am, n_basbas, n_basbas, mb_aux, nb_aux, &
  &             0, 0, my_blacs_ctxt_aux, MAX(1,max_row), info)

  allocate(invs_times_sqrtv(n_basbas,n_loc_prodbas), stat=info)
  call check_allocation(info, 'invs_times_sqrtv')

  call pdgemm('N', 'N', n_basbas, n_basbas, n_basbas, &
  &           1.0d0, ovlp_prodbas(1,1), 1, 1, sc_desc_am, &
  &           coulomb_matr(1,1), 1, 1, sc_desc_am, &
  &           0.d0, invs_times_sqrtv, 1, 1, sc_desc_am)
  coulomb_matr(:,:) = invs_times_sqrtv(:,:)

  ! Tidy up

  deallocate(invs_times_sqrtv)

end subroutine evaluate_invs_times_sqrtv_scalapack

!************
