!****s* FHI-aims/symmetry_reduce_k_points
!  NAME
!    symmetry_reduce_k_points
!  SYNOPSIS

subroutine symmetry_reduce_k_points(n_symmetries, symmats, &
&                                   n_k_xyz, k_off, k_number)

  !  PURPOSE
  !
  !    Routine to symmetrize the k-point grid according to a set of symmetry
  !    operations.
  !
  !  USES

  use localorb_io
  use bravais, only: check_if_lattice_point
  use geometry, only: recip_lattice_vector
  use dimensions, only: n_periodic,ik2irred_map
  use runtime_choices, only: symmetry_thresh
  use numerical_utilities, only: solve_LEQ
  use synchronize_mpi_basic, only: sync_logical
  use mpi_tasks, only: check_allocation, aims_stop, aims_warn, use_mpi, SYNC_OR
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_symmetries
  real*8, intent(IN) :: symmats(3, 3, n_symmetries)
  integer, intent(IN) :: n_k_xyz(3)
  real*8, intent(IN) :: k_off(3)
  integer, intent(OUT) :: k_number(n_k_xyz(1), n_k_xyz(2), n_k_xyz(3))

  !  INPUTS
  !    o n_symmetries -- Number of symmetry operations
  !    o symmats -- Symmetry operations
  !    o n_k_xyz -- Number of k-points
  !  OUTPUTS
  !    o k_number -- contains the number of k-points that have been mapped
  !                  on to the original entries, this is required to
  !                  determine the proper k-weights in the end
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
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_k_points_old
  integer :: n_unstabilized
  integer :: i_k1, i_k2, i_k3, i_symmetry, iik(3), trans_ik(3)
  real*8 :: k_tran(3), k_cur(3), k_diff(3), k_ind(3)
  logical :: found, is_equiv, stabiliser
  integer :: n_irred_uptonow, i_irred
  integer, allocatable :: irred2ik(:,:)
  integer, allocatable :: ik2irred(:,:,:)
  character*150 :: info_str
  integer :: info
  character(*), parameter :: func = 'symmetry_reduce_k_points'
  logical :: unsymmetric_k_mesh

  ! Initialize every k-point as representing only itself.
  k_number = 1

  unsymmetric_k_mesh = .false.

  if (n_symmetries <= 1) then
     write(info_str,'(2X,A)') &
     & "Found only identity as symmetry transformation; no k-point reduction."
     call localorb_info(info_str,use_unit,'(A)')
     return
  end if

  n_k_points_old = product(n_k_xyz)
  allocate(irred2ik(3, n_k_points_old), stat=info)
  call check_allocation(info, 'irred2ik', func)
  allocate(ik2irred(n_k_xyz(1), n_k_xyz(2), n_k_xyz(3)), stat=info)
  call check_allocation(info, 'ik2irred', func)
  ! Initialize with invalid/empty
  ik2irred = 0
  irred2ik = 0

  n_irred_uptonow = 0
  stabiliser = .false.
  n_unstabilized = 0
  do i_k1 = 1, n_k_xyz(1)
     do i_k2 = 1, n_k_xyz(2)
        do i_k3 = 1, n_k_xyz(3)
           k_ind(1) = dble(i_k1-1) / dble(n_k_xyz(1)) +  k_off(1)
           k_ind(2) = dble(i_k2-1) / dble(n_k_xyz(2)) +  k_off(2)
           k_ind(3) = dble(i_k3-1) / dble(n_k_xyz(3)) +  k_off(3)

           ! Moving from relative to absolute coordinates:
           k_cur = matmul(recip_lattice_vector, k_ind)

           found = .false.
           SYMMETRY_LOOP: do i_symmetry = 1, n_symmetries
              ! Make symmetry operator to k-point:
              k_tran = matmul(symmats(:,:, i_symmetry), k_cur)
              call solve_LEQ(func, 3, recip_lattice_vector, k_tran, k_ind)
              k_ind = (k_ind - k_off) * dble(n_k_xyz)
              if (any(k_ind - nint(k_ind) > 100*symmetry_thresh)) then
                 call aims_warn('Invalid rotated k-point indices detected', &
                                func)
                 unsymmetric_k_mesh = .true.
                 exit SYMMETRY_LOOP
              end if
              trans_ik = modulo(nint(k_ind), n_k_xyz) + 1  ! 1, ..., n_k_xyz(:)

              i_irred = ik2irred(trans_ik(1), trans_ik(2), trans_ik(3))
              if (i_irred > 0) then
                 ! Found the same k-point.
                 iik = irred2ik(:, i_irred)   ! irreducible k-point indices
                 if (any(iik /= trans_ik)) then
                    call aims_stop('Mismatch of irreducible k-point indices: irred2ik', func)
                 end if
                 ! Move occupation:
                 k_number(iik(1), iik(2), iik(3)) &
                 & = k_number(iik(1), iik(2), iik(3)) &
                 & + k_number(i_k1, i_k2, i_k3)
                 k_number(i_k1, i_k2, i_k3) = 0
                 ik2irred_map(i_k1, i_k2, i_k3) = &
                 ik2irred(iik(1), iik(2), iik(3))
                 ! This is here, so that all the k-points would not be
                 ! at the same corner in k-space:
                 ! JW: Why not?
                 if (stabiliser) then
                    k_number(i_k1, i_k2, i_k3) &
                    & = k_number(iik(1), iik(2), iik(3))
                    k_number(iik(1), iik(2), iik(3)) = 0
                    irred2ik(1, i_irred) = i_k1
                    irred2ik(2, i_irred) = i_k2
                    irred2ik(3, i_irred) = i_k3
                    ik2irred(i_k1, i_k2, i_k3) = i_irred
                    !ik2irred_map(i_k1, i_k2, i_k3) = &
                    !i_irred
                    ik2irred(iik(1), iik(2), iik(3)) = 0
                    stabiliser = .false.
                 else
                    n_unstabilized = n_unstabilized + 1
                    if (n_unstabilized == 3) then
                       stabiliser = .true.
                       n_unstabilized = 0
                    end if
                 end if
                 found = .true.
                 exit SYMMETRY_LOOP
              end if
           end do SYMMETRY_LOOP

           if (.not. found) then
              n_irred_uptonow = n_irred_uptonow + 1
              irred2ik(1, n_irred_uptonow) = i_k1
              irred2ik(2, n_irred_uptonow) = i_k2
              irred2ik(3, n_irred_uptonow) = i_k3
              ik2irred(i_k1, i_k2, i_k3) = n_irred_uptonow
              ik2irred_map(i_k1, i_k2, i_k3) = n_irred_uptonow
           end if
        end do ! i_k3
     end do !i_k2
  end do ! i_k1

  if (use_mpi) then 
     call sync_logical(unsymmetric_k_mesh, SYNC_OR)
  endif

  if (unsymmetric_k_mesh) then
     write(info_str,'(1X,A)') &
        "*** Error: An unsymmetric k-mesh has been detected."
     call localorb_info(info_str)
     write(info_str,'(1X,A)') &
        "    This error can be circumvent by explicitly setting 'symmetry_reduced_k_grid' to false"
     call localorb_info(info_str)
     write(info_str,'(1X,A)') &
        "    or going to a symmetric k-grid, e.g. Gamma centered 'k_offset 0 0 0'"
     call localorb_info(info_str)
     write(info_str,'(1X,A)') &
        "    However, nonsymmetric k-grids are normally a bad idea. We hope you know what you are doing."
     call localorb_info(info_str)
     stop
   endif

  write(info_str, '(2X,A,I8,A,I8)') &
  & '| k-points reduced from: ', n_k_points_old, ' to ', n_irred_uptonow
  call localorb_info(info_str,use_unit,'(A)')
  
  deallocate(irred2ik,ik2irred)

end subroutine symmetry_reduce_k_points
!******
