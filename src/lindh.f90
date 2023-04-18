!****h* FHI-aims/lindh
!  NAME
!    lindh
!  SYNOPSIS

module lindh

  !  PURPOSE
  !
  !    Calculate the Lindh model Hessian [2].
  !
  !    Elements not needed by the relaxation algorithm are not calculated.
  !
  !  USES

  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use bravais, only: get_n_supercells_pairs
  use constants, only: pi, bohr, hartree, radians
  use dimensions, only: n_atoms, n_periodic
  use distributed_hessian, only: nrow_hess, ncol_hess, offset, hess_id
  use localorb_io, only: localorb_info, use_unit, OL_low
  use mpi_tasks, only: aims_stop
  use numerical_utilities, only: pseudo_inverse
  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    [1] Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !      Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !      "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !      Computer Physics Communications 180, 2175 (2009).
  !    [2] Roland Lindh u. a., "On the use of a Hessian model function in
  !      molecular geometry optimizations", Chemical Physics Letters 241,
  !      423-428 (1995).
  !    [3] Vebjorn Bakken und Trygve Helgaker, "The efficient optimization
  !      of molecular geometries using redundant internal coordinates", The
  !      Journal of Chemical Physics 117, 9160-9174 (2002).
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !    the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  private

  public :: get_lindh_hessian

  ! --- elemental parameters

  integer, parameter :: n_elements = 18
  integer, parameter :: element_row(n_elements) = (/ &
  &  1,    1,    2,    2,    2,    2,    2,    2,    2,    2, &
  &              3,    3,    3,    3,    3,    3,    3,    3 /)
  integer, parameter :: max_row = 3

  ! --- hessian parameters

  ! from [2]
  real*8, parameter :: k_bond    = 0.450d0 ! Hartree / Bohr**2
  real*8, parameter :: k_bending = 0.150d0 ! Hartree
  real*8, parameter :: k_torsion = 0.005d0 ! Hartree

  real*8, parameter :: damping_gamma(max_row, max_row) = reshape( &
  & (/ 1.0000d0, 0.3949d0, 0.3949d0, &
  &    0.3949d0, 0.2800d0, 0.2800d0, &
  &    0.3949d0, 0.2800d0, 0.2800d0  /), (/3, 3/))
  real*8, parameter :: damping_len(max_row, max_row) = reshape( &
  & (/ 1.35d0, 2.10d0, 2.53d0, &
  &    2.10d0, 2.87d0, 3.40d0, &
  &    2.53d0, 3.40d0, 3.40d0  /), (/3, 3/))

contains

  !----------------------------------------------------------------------------
  !****s* lindh/get_lindh_hessian
  !  NAME
  !    get_lindh_hessian
  !  SYNOPSIS

  subroutine get_lindh_hessian(coords, lattice_vector, thres, full_hessian)

    !  PURPOSE
    !
    !    Get Lindh Hessian matrix.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: coords(3, n_atoms)
    real*8, intent(IN) :: lattice_vector(3, 3)
    real*8, intent(IN) :: thres
    real*8, intent(OUT) :: full_hessian(nrow_hess, ncol_hess)

    !  INPUTS
    !    o coords -- atomic coords in unit cell
    !    o lattice_vector -- Bravais vectors
    !    o thres -- Threshold [neglect terms smaller than ~exp(-thres)]
    !               Try something like 15.d0 if in doubt.
    !  OUTPUTS
    !    o full_hessian -- Full Lindh model Hessian
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Rmax, Vmax
    integer :: n_supercells(3)
    integer :: n_max_neighbor
    real*8, parameter :: model_spacing = 1.d0  ! abohr
    integer :: i_atom_1, i_atom_2, i_atom_3, i_atom_4
    integer :: i_neighbor_2, i_neighbor_3, i_neighbor_4
    integer :: avec12(3), avec13(3), avec14(3), avec24(3)
    integer :: chain2atom(4), chain2a(3,4)
    real*8 :: damp12, damp23, damp13, damp34, damp14
    real*8 :: thres_bending, thres_torsion, prefac
    real*8 :: pos(3, 4)
    real*8 :: q, dq(3, 4), dq2(3, 3)
    integer, allocatable :: atom2n_neighbor(:)
    integer, allocatable :: neighbor2atom(:,:), neighbor2cell(:,:,:)
    real*8, allocatable :: neighbor2damp(:,:)
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lindh_hessian'

    call localorb_info('Setting up Lindh matrix', use_unit, "(2X,A)", OL_low)

    full_hessian = 0.d0

    ! Reduced thresholds for less pronounced and more expensive coordinates:
    ! Take scaling from the old Lindh.py (-> /bohr**2)
    thres_bending = thres + log(k_bending / (k_bond/bohr**2))
    thres_torsion = thres + log(k_torsion / (k_bond/bohr**2))

    write(info_str, "(A,F8.4,A,F8.4,A,F8.4,A)") &
    & 'Use threshold of', thres, ' for bonds,', &
    & thres_bending, ' for bendings, and ', thres_torsion, ' for torsions.'
    call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)

    ! --- Get neighbor list

    ! Set n_max_neighbor to a value corresponding to a cubic grid with
    ! a lattice constant of model_spacing
    Rmax = sqrt(thres / minval(damping_gamma))
    Vmax = 4.d0 * pi / 3.d0 * Rmax**3
    n_max_neighbor = ceiling(Vmax) / model_spacing**3
    ! In principle, one could use n_supercells=number_of_super_cells, I guess.
    ! But I found it to be 0 at this stage of calculations.  So it was easier
    ! (and somehow even cleaner) to write a subroutine to get values of
    ! n_supercells which really match the problem.
    n_supercells = 0  ! Important for n_periodic < 3 (e.g. == 0)
    call get_n_supercells_pairs(n_periodic, lattice_vector, n_atoms, coords, &
    &                           Rmax, n_supercells)
    write(info_str,"('Number of supercells:',3I5)") n_supercells
    call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)

    call aims_allocate(atom2n_neighbor, n_atoms, "atom2n_neighbor")
    call aims_allocate(neighbor2atom, n_max_neighbor, n_atoms, "neighbor2atom")
    call aims_allocate(neighbor2cell, 3, n_max_neighbor, n_atoms, "neighbor2cell")
    call aims_allocate(neighbor2damp, n_max_neighbor, n_atoms, "neighbor2damp")

    call get_neighbor_lists(thres, coords, lattice_vector, n_supercells, &
    &                       n_max_neighbor, atom2n_neighbor, &
    &                       neighbor2atom, neighbor2cell, neighbor2damp)

    ! --- Go over internal coordinates

    ATOM_1: do i_atom_1 = 1, n_atoms
       pos(:, 1) = coords(:, i_atom_1)
       chain2atom(1) = i_atom_1
       chain2a(:, 1) = (/0, 0, 0/)

       write(info_str,"(A,I5,A,I5,A)") &
       & 'Atom', i_atom_1, ' has', atom2n_neighbor(i_atom_1), ' neighbors.'
       call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_low)


       ATOM_2: do i_neighbor_2 = 1, atom2n_neighbor(i_atom_1)
          i_atom_2 = neighbor2atom(i_neighbor_2, i_atom_1)
          avec12 = neighbor2cell(:, i_neighbor_2, i_atom_1)
          damp12 = neighbor2damp(i_neighbor_2, i_atom_1)
          if (damp12 > thres) exit ATOM_2
          pos(:, 2) = coords(:, i_atom_2) + matmul(lattice_vector, avec12)
          chain2atom(2) = i_atom_2
          chain2a(:, 2) = avec12

          if (positive_pair(i_atom_1, i_atom_2, avec12)) then
             prefac = k_bond * exp(- damp12)
             call bond2dq(pos(:, 1:2), q, dq(:, 1:2))
             call add_to_hessian(full_hessian, prefac, 2, &
             &                   chain2atom(1:2), chain2a(:, 1:2), dq(:, 1:2))
             write(info_str, &
             &     "(A,3X,I5,' (',3I2,'):',F5.2,A,F6.2,A,F8.3,A)") &
             & 'Bond to', i_atom_2, avec12, q*bohr, ' A,   damped by', damp12,&
             & '; prefac =', prefac*hartree/bohr**2, ' eV/A^2.'
             call localorb_info(info_str, use_unit, "(2X,'| ',1X,A)", OL_low)
          end if

          ATOM_3: do i_neighbor_3 = 1, atom2n_neighbor(i_atom_2)
             i_atom_3 = neighbor2atom(i_neighbor_3, i_atom_2)
             avec13 = avec12 + neighbor2cell(:, i_neighbor_3, i_atom_2)
             damp23 = neighbor2damp(i_neighbor_3, i_atom_2)
             if (damp23 > thres) exit ATOM_3
             if (i_atom_1 == i_atom_3 .and. all(avec13 == 0)) cycle
             damp13 = damp12 + damp23
             pos(:, 3) = coords(:, i_atom_3) + matmul(lattice_vector, avec13)
             chain2atom(3) = i_atom_3
             chain2a(:, 3) = avec13

             if (damp13 < thres_bending .and. &
             &   positive_pair(i_atom_1, i_atom_3, avec13)) then
                ! For nearly linear bends, dq2 contains the derivative
                ! of an orthogonally defined bending angle.
                prefac = k_bending * exp(- damp13)
                call bending2dq(pos(:, 1:3), q, dq(:, 1:3), dq2)
                write(info_str, &
                &     "(A,I5,' (',3I2,'):',I5,A,F6.2,A,F8.3,A)") &
                & 'Bending to', i_atom_3, avec13, nint(q*radians), &
                & ' deg, damped by', damp13, &
                & '; prefac =', prefac*hartree, ' eV.'
                call add_to_hessian(full_hessian, prefac, 3, &
                &                chain2atom(1:3), chain2a(:, 1:3), dq(:, 1:3))
                if (any(dq2 /= 0.d0)) then    ! 0.d0 signals absence
                   write(info_str, "(A,A)") trim(info_str), ' [Nearly linear]'
                   call add_to_hessian(full_hessian, prefac, 3, &
                   &                   chain2atom(1:3), chain2a(:, 1:3), dq2)
                end if
                call localorb_info(info_str, use_unit, "(2X,'| ',2X,A)", OL_low)
             end if

             ATOM_4: do i_neighbor_4 = 1, atom2n_neighbor(i_atom_3)
                i_atom_4 = neighbor2atom(i_neighbor_4, i_atom_3)
                avec14 = avec13 + neighbor2cell(:, i_neighbor_4, i_atom_3)
                damp34 = neighbor2damp(i_neighbor_4, i_atom_3)
                if (damp34 > thres) exit ATOM_4
                if (i_atom_1 == i_atom_4 .and. all(avec14 == 0)) cycle
                avec24 = avec14 - avec12
                if (i_atom_2 == i_atom_4 .and. all(avec24 == 0)) cycle
                damp14 = damp13 + damp34
                pos(:, 4) = coords(:, i_atom_4) + matmul(lattice_vector, avec14)
                chain2atom(4) = i_atom_4
                chain2a(:, 4) = avec14

                if (damp14 < thres_torsion .and. &
                &   positive_pair(i_atom_1, i_atom_4, avec14)) then
                   prefac = k_torsion * exp(- damp14)
                   call torsion2dq(pos, q, dq)
                   write(info_str, &
                   &     "(A,I5,' (',3I2,'):',I5,A,F6.2,A,F8.3,A)") &
                   & 'Torsion to', i_atom_4, avec14, nint(q*radians), &
                   & ' deg, damped by', damp14, &
                   & '; prefac =', prefac*hartree, ' eV.'
                   if (any(dq /= 0.d0)) then    ! 0.d0 signals absence
                      call add_to_hessian(full_hessian, prefac, 4, &
                      &                   chain2atom, chain2a, dq)
                   else
                      write(info_str, "(A,A)") &
                      & trim(info_str), ' [Nearly linear]'
                   end if
                   call localorb_info(info_str, use_unit, &
                   &                  "(2X,'| ',3X,A)", OL_low)
                end if
             end do ATOM_4
          end do ATOM_3
       end do ATOM_2
    end do ATOM_1

    call aims_deallocate(atom2n_neighbor, "atom2n_neighbor")
    call aims_deallocate(neighbor2atom, "neighbor2atom")
    call aims_deallocate(neighbor2cell, "neighbor2cell")
    call aims_deallocate(neighbor2damp, "neighbor2damp")

  end subroutine get_lindh_hessian
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/bond2dq
  !  NAME
  !    bond2dq
  !  SYNOPSIS

  subroutine bond2dq(pos, q, dq)

    !  PURPOSE
    !
    !    Calculate the bond length q and its derivative dq wrt to the atomic
    !    positions.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: pos(3, 2)
    real*8, intent(OUT) :: q
    real*8, intent(OUT) :: dq(3, 2)

    !  INPUTS
    !    o pos -- pos(:,i): coordinates of i-th atom
    !  OUTPUTS
    !    o q -- bond length
    !    o dq -- dq / dpos
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Avec(3), Asq, dq_A(3)

    character(*), parameter :: func = 'bond2dq'

    Avec = pos(:, 2) - pos(:, 1)
    Asq = dot_product(Avec, Avec)
    if (Asq < 1d-10) call aims_stop('Zero-length bond', func)
    q = sqrt(Asq)
    dq_A = Avec / q

    dq(:, 1) = - dq_A   ! = dq/dAvec * dAvec / dpos1
    dq(:, 2) = + dq_A   ! = dq/dAvec * dAvec / dpos2

  end subroutine bond2dq
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/bending2dq
  !  NAME
  !    bending2dq
  !  SYNOPSIS

  subroutine bending2dq(pos, q, dq, dq2)

    !  PURPOSE
    !
    !    Calculate the bending angle q and its derivative dq wrt to the atomic
    !    positions.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: pos(3, 3)
    real*8, intent(OUT) :: q
    real*8, intent(OUT) :: dq(3, 3)
    real*8, intent(OUT) :: dq2(3, 3)

    !  INPUTS
    !    o pos -- pos(:,i): coordinates of i-th atom
    !  OUTPUTS
    !    o q -- bending angle
    !    o dq -- dq / dpos
    !    o dq2 -- if nonzero, derivative of orthogonal angle
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Avec(3), Bvec(3)
    real*8 :: Asq, Bsq, lambda_u, lambda_v
    real*8 :: uvec(3), vvec(3)
    real*8 :: u_dot_v, sinsq_uv, wsq
    real*8 :: dq_up(3), dq_vp(3)
    real*8 :: wvec(3), wvec_tmp(3)
    real*8 :: mydq(3, 3, 2)
    integer :: n_get, i_get

    character(*), parameter :: func = 'bending2dq'

    ! The notation of this subroutine is similar to the one in [3].

    Avec = pos(:, 2) - pos(:, 1)
    Asq = dot_product(Avec, Avec)
    Bvec = pos(:, 3) - pos(:, 2)
    Bsq = dot_product(Bvec, Bvec)
    if (min(Asq,Bsq) < 1d-10) call aims_stop('Zero-length bond', func)

    ! Get bond angle q
    lambda_u = sqrt(Asq)
    lambda_v = sqrt(Bsq)
    uvec = - Avec / lambda_u
    vvec =   Bvec / lambda_v
    u_dot_v = dot_product(uvec, vvec)
    call abs_into_one(u_dot_v)  ! u_dot_v = 1+eps -> 1
    q = acos(u_dot_v)

    ! Check for nearly linear angle
    sinsq_uv = abs(1.d0-u_dot_v**2)
    if (sinsq_uv < 0.01d0) then
       n_get = 2
    else
       n_get = 1
    end if

    ! Get derivative(s)
    mydq = 0.d0
    do i_get = 1, n_get
       if (i_get == 1) then
          call cross_product(uvec, vvec, wvec)
          wsq = dot_product(wvec, wvec)
          if (wsq < 1d-10) then
             call cross_product(uvec, (/1.d0, -1.d0, 1.d0/), wvec)
             wsq = dot_product(wvec, wvec)
             if (wsq < 1d-5) then
                call cross_product(uvec, (/-1.d0, 1.d0, 1.d0/), wvec)
                wsq = dot_product(wvec, wvec)
             end if
          end if
       else ! i_get == 2
          wvec_tmp = wvec
          call cross_product(uvec, wvec_tmp, wvec)
          wsq = dot_product(wvec, wvec)
       end if
       if (wsq < 1d-10) call aims_stop('wsq still invalid', func)

       wvec = wvec / sqrt(wsq)
       call cross_product(uvec, wvec, dq_up)
       dq_up = dq_up / lambda_u
       call cross_product(wvec, vvec, dq_vp)
       dq_vp = dq_vp / lambda_v

       mydq(:, 1, i_get) =   dq_up
       mydq(:, 2, i_get) = - dq_up - dq_vp
       mydq(:, 3, i_get) =         + dq_vp
    end do

    dq = mydq(:,:, 1)
    dq2 = mydq(:,:, 2)  ! zero in standard case of n_get==1

  end subroutine bending2dq
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/torsion2dq
  !  NAME
  !    torsion2dq
  !  SYNOPSIS

  subroutine torsion2dq(pos, q, dq)

    !  PURPOSE
    !
    !    Calculate the torsional angle q and its derivative dq wrt to the
    !    atomic positions.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: pos(3, 4)
    real*8, intent(OUT) :: q
    real*8, intent(OUT) :: dq(3, 4)

    !  INPUTS
    !    o pos -- pos(:,i): coordinates of i-th atom
    !  OUTPUTS
    !    o q -- torsional angle
    !    o dq -- dq / dpos
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Avec(3), Bvec(3), Cvec(3)
    real*8 :: Asq, Bsq, Csq, lambda_u, lambda_v, lambda_w
    real*8 :: uvec(3), vvec(3), wvec(3)
    real*8 :: ucw(3), vcw(3), ucv(3)
    real*8 :: cphiu, sphiu, cphiv, sphiv, cosq
    real*8 :: dq_up(3), dq_vp(3), dq_wp(3)

    character(*), parameter :: func = 'torsion2dq'

    ! The notation of this subroutine is similar to the one in [3].

    Avec = pos(:, 2) - pos(:, 1)
    Asq = dot_product(Avec, Avec)
    Bvec = pos(:, 3) - pos(:, 2)
    Bsq = dot_product(Bvec, Bvec)
    Cvec = pos(:, 4) - pos(:, 3)
    Csq = dot_product(Cvec, Cvec)
    if (min(Asq,Bsq,Csq) < 1d-10) call aims_stop('Zero-length bond', func)

    lambda_u = sqrt(Asq)
    lambda_v = sqrt(Csq)
    lambda_w = sqrt(Bsq)
    uvec = - Avec / lambda_u
    vvec =   Cvec / lambda_v
    wvec =   Bvec / lambda_w
    call cross_product(uvec, wvec, ucw)
    call cross_product(vvec, wvec, vcw)
    ! phi_u
    cphiu = dot_product(uvec, wvec)
    call abs_into_one(cphiu)
    sphiu = sqrt(1.d0 - cphiu*cphiu)
    ! phi_v
    cphiv = - dot_product(vvec, wvec)
    call abs_into_one(cphiv)
    sphiv = sqrt(1.d0 - cphiv*cphiv)

    if (sphiu**2 < 1d-6 .or. sphiv**2 < 1d-6) then
       ! Either uvec or vvec is parallel to wvec.
       !   => The torsional angle is ill-defined.
       ! In principle, torsional angles should be weighted by something like
       ! sin(phi_{u,v})**2 to get this smooth.  But for now, smoothness of the
       ! Hessian matrix is not an issue.
       q = 0.d0
       dq = 0.d0
    else
       cosq = dot_product(ucw, vcw) / (sphiu*sphiv)
       call abs_into_one(cosq)
       call cross_product(uvec, vvec, ucv)
       if (dot_product(ucv, wvec) > 0.d0) then
          q =   acos(cosq)
       else
          q = - acos(cosq)
       end if
       dq_up =   ucw / (lambda_u * sphiu*sphiu)
       dq_vp = - vcw / (lambda_v * sphiv*sphiv)
       dq_wp = - ucw * cphiu / (lambda_w * sphiu*sphiu) &
       &       - vcw * cphiv / (lambda_w * sphiv*sphiv)

       dq(:, 1) = dq_up
       dq(:, 2) = - dq_wp - dq_up
       dq(:, 3) = + dq_wp - dq_vp
       dq(:, 4) = dq_vp
    end if

  end subroutine torsion2dq
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/get_neighbor_lists
  !  NAME
  !    get_neighbor_lists
  !  SYNOPSIS

  subroutine get_neighbor_lists(thres, coords, lattice_vector, n_supercells, &
  &                             n_max_neighbor, atom2n_neighbor, &
  &                             neighbor2atom, neighbor2cell, neighbor2damp)

    !  PURPOSE
    !
    !    Get a mapping of the atoms to their neighbors.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: thres
    real*8, intent(IN) :: coords(3, n_atoms)
    real*8, intent(IN) :: lattice_vector(3, 3)
    integer, intent(IN) :: n_supercells(3)
    integer, intent(IN) :: n_max_neighbor
    integer, intent(OUT) :: atom2n_neighbor(n_atoms)
    integer, intent(OUT) :: neighbor2atom(n_max_neighbor, n_atoms)
    integer, intent(OUT) :: neighbor2cell(3, n_max_neighbor, n_atoms)
    real*8, intent(OUT) :: neighbor2damp(n_max_neighbor, n_atoms)

    !  INPUTS
    !    o thres -- maximum allowed damping
    !    o coords -- atomic coords in unit cell
    !    o lattice_vector -- Bravais vectors
    !    o n_max_neighbor -- first dimension of neighbor_list
    !  OUTPUTS
    !    o atom2n_neighbor -- number of neighbors per atom
    !    o neighbor2atom -- atom number of neighbor
    !    o neighbor2cell -- Bravais indices of neighbor
    !    o neighbor2damp -- damping of neighborship relation [sorted increasing]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_atom, j_atom, i_row, j_row
    integer :: a1, a2, a3, avec(3)
    integer :: n_neighbor_uptonow, i_old, i_new
    real*8 :: dist_vec(3), dist_sq, damp
    real*8 :: n2damp(n_max_neighbor)
    integer :: n2atom(n_max_neighbor), n2cell(3, n_max_neighbor)
    integer :: new2old(n_max_neighbor), old2new(n_max_neighbor)
    character(*), parameter :: func = 'get_neighbor_lists'

    I_ATOM_L: do i_atom = 1, n_atoms
       i_row = atom_to_row(i_atom)
       n_neighbor_uptonow = 0
       J_ATOM_L: do j_atom = 1, n_atoms
          j_row = atom_to_row(j_atom)
          do a1 = -n_supercells(1), n_supercells(1)
             avec(1) = a1
             do a2 = -n_supercells(2), n_supercells(2)
                avec(2) = a2
                do a3 = -n_supercells(3), n_supercells(3)
                   avec(3) = a3
                   if (i_atom == j_atom .and. all(avec == 0)) cycle
                   dist_vec = coords(:, j_atom) &
                   &          + matmul(lattice_vector, avec) &
                   &          - coords(:, i_atom)
                   dist_sq = dot_product(dist_vec, dist_vec)
                   damp = pair_to_damping(dist_sq, i_row, j_row)
                   if (damp < thres) then
                      if (n_neighbor_uptonow >= n_max_neighbor) then
                         call aims_stop('Too many neighbors', func)
                      end if
                      n_neighbor_uptonow = n_neighbor_uptonow + 1
                      n2atom(n_neighbor_uptonow) = j_atom
                      n2cell(:, n_neighbor_uptonow) = avec
                      n2damp(n_neighbor_uptonow) = damp
                   end if
                end do
             end do
          end do
       end do J_ATOM_L
       call insertionsort(n2damp, n_neighbor_uptonow, old2new, new2old)
       do i_new = 1, n_neighbor_uptonow
          neighbor2damp(i_new, i_atom) = n2damp(i_new)
          i_old = new2old(i_new)
          neighbor2atom(i_new, i_atom) = n2atom(i_old)
          neighbor2cell(:, i_new, i_atom) = n2cell(:, i_old)
       end do
       atom2n_neighbor(i_atom) = n_neighbor_uptonow
    end do I_ATOM_L

  end subroutine get_neighbor_lists
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/add_to_hessian
  !  NAME
  !    add_to_hessian
  !  SYNOPSIS

  subroutine add_to_hessian(full_hessian, prefac, n_chain, &
  &                         chain2atom, chain2a, dq)

    !  PURPOSE
    !
    !    Make a rank-1 update to the given Hessian:
    !
    !       H += prefac * dq * dq^T
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: full_hessian(nrow_hess, ncol_hess)
    real*8, intent(IN) :: prefac
    integer, intent(IN) :: n_chain
    integer, intent(IN) :: chain2atom(1:n_chain)
    integer, intent(IN) :: chain2a(3, 1:n_chain)
    real*8, intent(IN) :: dq(3, n_chain)

    !  INPUTS
    !    o full_hessian -- Hessian to be updated
    !    o prefac -- prefactor
    !    o n_chain -- number of affected atoms (bond: 2, bending:3, torsion: 4)
    !    o chain2atom -- affected atom numbers
    !    o chain2a -- Lattice vector indices of affected centers
    !    o dq -- actual vector (derivative of q wrt atomic coordinates)
    !  OUTPUTS
    !    o full_hessian -- updated Hessian
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_chain_1, i_chain_2, i_atom_1, i_atom_2, i_1, i_2
    integer :: i_row
    integer :: i_col
    integer :: j_row
    integer :: j_col
    real*8 :: vec_1(3), vec_2(3)
    character(*), parameter :: func = 'add_to_hessian'

    do i_chain_1 = 1, n_chain
       i_atom_1 = chain2atom(i_chain_1)
       if (i_atom_1 <= n_atoms) then
          vec_1 = dq(:, i_chain_1)
       else
          ! Change of lattice vector i means chain2a(i, i_chain_1)-fold change
          ! in center_coord, and thus:
          vec_1 = dq(:, i_chain_1) * chain2a(:, i_chain_1)
       end if
       do i_chain_2 = 1, n_chain
          i_atom_2 = chain2atom(i_chain_2)
          if (i_atom_2 <= n_atoms) then
             vec_2 = dq(:, i_chain_2)
          else
             vec_1 = dq(:, i_chain_2) * chain2a(:, i_chain_2)
          end if

          do i_1 = 1, 3
             do i_2 = 1, 3
                i_row = (i_atom_1 - 1) * 3 + i_1
                i_col = (i_atom_2 - 1) * 3 + i_2

                if (hess_id(i_row) > 0 .and. hess_id(i_col) /= 0) then
                   j_row = hess_id(i_row)

                   if (hess_id(i_col) > 0) then
                      j_col = hess_id(i_col)+offset
                   else
                      j_col = -hess_id(i_col)
                   end if

                   full_hessian(j_row, j_col) = full_hessian(j_row, j_col) &
                      & + prefac * vec_1(i_1) * vec_2(i_2)
                end if
             end do
          end do
       end do
    end do

  end subroutine add_to_hessian
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/positive_pair
  !  NAME
  !    positive_pair
  !  SYNOPSIS

  logical function positive_pair(i_atom_1, i_atom_2, avec12) result(is_positive)

    !  PURPOSE
    !
    !    Defines some (arbitrary, but consistent) "sign" to a given atom pair
    !    which differs for transposed pairs.  It is invalid to call this
    !    function with two identical atoms.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_atom_1, i_atom_2
    integer, intent(IN) :: avec12(3)

    !  INPUTS
    !    o i_atom_1, i_atom_2 -- atomic numbers
    !    o avec12 -- Bravais indices of second atom
    !  OUTPUTS
    !    o is_positive -- sign of pair
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i
    character(*), parameter :: func = 'positive_pair'

    if (i_atom_1 == i_atom_2 .and. all(avec12 == 0)) then
       call aims_stop('Called with self-referencing pair', func)
    end if

    if (i_atom_1 < i_atom_2) then
       is_positive = .true.
    else if (i_atom_1 > i_atom_2) then
       is_positive = .false.
    else
       is_positive = .false.   ! appease compiler
       do i = 1, 3
          if (avec12(i) > 0) then
             is_positive = .true.
             exit
          else if (avec12(i) < 0) then
             is_positive = .false.
             exit
          end if
       end do
    end if

  end function positive_pair
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/pair_to_damping
  !  NAME
  !    pair_to_damping
  !  SYNOPSIS

  real*8 function pair_to_damping(dist_sq, i_row, j_row) result(damp)

    !  PURPOSE
    !
    !    Get damping factor for the given pair of atoms.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: dist_sq
    integer, intent(IN) :: i_row, j_row

    !  INPUTS
    !    o dist_sq -- squared distance between the two atoms in question
    !    o i_row, j_row -- rows of
    !  OUTPUTS
    !    o damp -- damping factor [to be put into exp(- damp)]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: gamma, def_len_sq
    character(*), parameter :: func = 'pair_to_damping'

    gamma = damping_gamma(i_row, j_row)
    def_len_sq = damping_len(i_row, j_row)**2
    damp = gamma * (dist_sq - def_len_sq)

  end function pair_to_damping
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/atom_to_row
  !  NAME
  !    atom_to_row
  !  SYNOPSIS

  integer function atom_to_row(i_atom) result(i_row)

    !  PURPOSE
    !
    !    For a given atom number, return the row in the periodic table (ceiled
    !    by max_row).
    !
    !  USES

    use geometry, only: species
    use species_data, only: species_z
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_atom

    !  INPUTS
    !    o i_atom -- atom number
    !  OUTPUTS
    !    o i_row -- row in periodic table
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_z
    character(*), parameter :: func = 'atom_to_row'

    i_row = 0  ! Appease compiler

    ! Yes, fractional nuclear charges have been observed to occur in actual
    ! input files.
    i_z = nint(species_z(species(i_atom)))
    if (i_z < 0) then
       call aims_stop('Invalid atomic number', func)
    else if (i_z == 0) then
       i_row = 1  ! Should not move, but what should we do here?
    else if (1 <= i_z .and. i_z <= n_elements) then
       i_row = element_row(i_z)
    else
       i_row = max_row
    end if

  end function atom_to_row
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/cross_product
  !  NAME
  !    cross_product
  !  SYNOPSIS

  subroutine cross_product(vec_1, vec_2, vec_out)

    !  PURPOSE
    !    Calculate the cross product
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: vec_1(3), vec_2(3)
    real*8, intent(OUT) :: vec_out(3)

    !  INPUTS
    !    o vec_1, vec_2 -- R^3 vectors
    !  OUTPUTS
    !    o vec_out -- vec_1 x vec_2
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'cross_product'
    !-

    vec_out(1) = vec_1(2)*vec_2(3) - vec_1(3)*vec_2(2)
    vec_out(2) = vec_1(3)*vec_2(1) - vec_1(1)*vec_2(3)
    vec_out(3) = vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1)

  end subroutine cross_product
  !******
  !----------------------------------------------------------------------------
  !****s* lindh/abs_into_one
  !  NAME
  !    abs_into_one
  !  SYNOPSIS

  subroutine abs_into_one(x)

    !  PURPOSE
    !
    !    Make sure that abs(x) <= 1.d0 by cutting it towards -1.d0 or 1.d0.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: x

    !  INPUTS
    !    o x -- number roughly within [-1, 1]
    !  OUTPUTS
    !    o x -- number strictly within [-1, 1]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'abs_into_one'

    x = min(x,  1.d0)
    x = max(x, -1.d0)

  end subroutine abs_into_one
  !******
end module lindh
