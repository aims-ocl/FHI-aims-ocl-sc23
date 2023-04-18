!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  The core subroutine for calculating magnetic response
!!  properties. The result is put into a response tensor and printed
!!  in MR_main. See the documentation for theory and other details.
!!
!!  The structure of this subroutine is the following. There are two
!!  large blocks - for evaluating first-order (STEP 1 below) and
!!  second-order responses (STEP 2). If the response has the form
!!  2Re<Psi1|H1|Psi0>, e.g. paramagnetic shielding, the first-order
!!  wavefunctions (|Psi1>) are found first, followed by multiplication
!!  with the H1|Psi0> term. If the response has the form
!!  <Psi0|H2|Psi0>, e.g. diamagnetic shielding, the first part is
!!  skipped and the expectation between the unperturbed wavefunctions
!!  can be taken directly.
!!
!!  This subroutine repeatedly calls 'integrate', which integrates a
!!  given function over all grid points to produce the matrix elements
!!  between basis functions. Functions to be integrated are found in
!!  integrands.f90
!!
!!  First-order responses are calculated using density functional
!!  perturbation theory (main code in DFPT.f90).
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
!!  COMMENTS
!!
!!  We use the S.I. atomic units.
!!
subroutine MR_core(term_name, timer)

  use aims_memory_tracking, only: aims_allocate, push_current_memory, &
       & pop_memory_estimate, update_when_allocating
  use constants,            only: light_speed_sq, pi
  use DFPT,                 only: do_DFPT, cleanup_DFPT
  use dimensions,           only: n_spin
  use geometry,             only: coords
  use integration,          only: integrate
  use localorb_io,          only: localorb_multi
  use integrands
  use MR_global
  use physics,              only: KS_eigenvalue
  use tools,                only: antisymmetrize, mul, safe_deallocate, &
       & start_wall, stop_wall, str
  use runtime_choices,      only: flag_rel, REL_none
  use scalapack_wrapper,    only: eigenvec, mxld, mxcol, my_row, my_col

  implicit none

  character(*), intent(in) :: term_name
  ! Timer for the non-self-consistent part of the first-order
  ! Hamiltonian (the DFPT part has its own timers)
  real(dp), intent(in out) :: timer(4)

  real(dp), parameter :: alpha2 = 1/light_speed_sq
  real(dp), parameter :: alpha4 = 1/light_speed_sq**2

  character(*), parameter :: THIS_SUB = 'MR_core::'
  ! General purpose 2D block cyclic distributed matrix. Global
  ! dimensions: (n_basis,n_basis,n_dims), where n_dims is a compact
  ! index for directions, spins, and any other dimensions. For
  ! example, for the diamagnetic shielding tensor the number of
  ! dimensions is 9 if calc_full_tensor==.true. and 3 otherwise. For
  ! spin-polarized paramagnetic SO it is 6 (2 spins x 3 directions).
  real(dp), allocatable :: matrix_BC(:,:,:)
  ! H|Psi>, where H is some Hamiltonian (unperturbed, first-order,
  ! second-order, ...) and Psi are some wavefunctions (unperturbed or
  ! perturbed). Global dimensions:
  ! (n_basis,n_basis,n_spin,3,size(active_nuclei))
  real(dp), allocatable :: H_Psi(:,:,:,:,:)
  ! First-order wavefunctions. Global dimensions:
  ! (n_basis,n_basis,n_spin,3,size(active_nuclei)), where 2 and 3 are spins
  ! and spatial directions.
  real(dp), allocatable :: Psi1(:,:,:,:,:)
  ! First-order GIAO overlap matrix. Global dimensions: (n_basis,n_basis,3)
  real(dp), allocatable :: overlap1(:,:,:)
  ! Second-order GIAO overlap matrix. Global dimensions: (n_basis,n_basis,3)
  real(dp), allocatable :: overlap2(:,:,:)
  ! First-order eigenvalue matrix. Global dimensions: (n_basis,n_basis,n_spin,3)
  real(dp), allocatable :: epsilon1(:,:,:,:)
  ! Work array for a response tensor. Results are then copied from
  ! this to the correct array.
  real(dp) :: tmp_tensor(3,3)
  ! 9 if calc_full_tensor, 3 otherwise
  integer :: n_tensor_components
  ! Whether the current perturbation is purely imaginary. One of the
  ! consequences of imaginary perturbations is that there is no
  ! first-order density.
  logical :: imaginary_pert
  ! Number of spatial directions of Psi1 that are simultaneously
  ! processed. This is 1 for isotropic perturbations and 3 otherwise.
  integer :: n_directions
  ! Number of perturbations considered in the DFPT part. This is
  ! size(active_nuclei) for J-couplings and 1 otherwise (i.e., when
  ! B-field is the perturbation).
  integer :: n_psi1
  ! Whether to skip the DFPT part (e.g., when we are calculating only
  ! a diamagnetic term)
  logical :: skip_dfpt
  ! Whether to compute the first-order GIAO overlap matrix
  logical :: include_overlap1
  ! Whether to compute the second-order GIAO overlap matrix
  logical :: include_overlap2
  ! Whether to compute the first-order GIAO eigenvalue matrix
  logical :: include_epsilon1
  ! Whether to compute H|Psi> in STEP 1 and save it for STEP 2. This
  ! applies to the case where the same perturbation is used for
  ! |Psi1[j]> and H1|Psi1[i]> in 2Re<Psi1[j]|H1|Psi[i]> (e.g.,
  ! J-couplings). Also, from the way the loop is performed in STEP 2,
  ! the index j is always above 1 for J-couplings.
  logical :: save_H_psi
  ! For J-couplings, the first-order wavefunctions are required for
  ! all but the last atom in the list active_nuclei. This results from
  ! the way the loops are constructed in STEP 2. When computing
  ! <Psi1[i]|H1|Psi0[j]>, the index i never reaches the last atom.
  logical :: skip_last_dfpt
  ! Whether the perturbation contains the Sz operator. If so, this
  ! information must be passed to the DFPT cycle.
  logical :: spin_operator
  ! Whether to calculate the first-order Hartree potential
  logical :: calc_H1_Ha
  ! Counters
  integer :: i_spin, i_atom, j_atom, i_dir, j_dir, i_row, i_col, shift
  ! For matrix_BC, the aims_allocate subroutine needs to be
  ! inlined. Needed for circumventing PGI bug.
  integer(8) :: mem_matrix_BC
  integer :: info

  ! STEP 0 - Preparatory tasks
  imaginary_pert = any(term_name == [character(28) :: &
       & 'Paramagnetic spin-orbit', 'Paramagnetic shielding', &
       & 'Paramagnetic magnetizability'])
  if (calc_full_tensor) then
     n_tensor_components = 9
  else
     n_tensor_components = 3
  end if
  if (term_name == 'Fermi contact') then
     n_directions = 1
  else
     n_directions = 3
  end if
  if (any(term_name == [character(28) :: 'Fermi contact', &
       & 'Paramagnetic spin-orbit', 'Spin-dipole', 'Diamagnetic spin-orbit'])) &
       & then
     n_psi1 = size(active_nuclei)
  else
     n_psi1 = 1
  end if
  include_overlap1 = any(term_name == [character(28) :: &
       & 'Paramagnetic shielding', 'Paramagnetic magnetizability']) .and. &
       & .not. no_giao
  include_overlap2 = term_name == 'Diamagnetic magnetizability' .and. &
       & .not. no_giao
  include_epsilon1 = term_name == 'Paramagnetic magnetizability' .and. &
       & .not. no_giao
  save_H_Psi = any(term_name == [character(28) :: 'Paramagnetic spin-orbit', &
       & 'Fermi contact', 'Spin-dipole', 'Diamagnetic spin-orbit', &
       & 'Paramagnetic magnetizability'])
  skip_dfpt = any(term_name == [character(28) :: 'Diamagnetic spin-orbit', &
       & 'Diamagnetic magnetizability', 'Diamagnetic shielding', &
       & 'Electric field gradient'])
  skip_last_dfpt = any(term_name == [character(28) :: 'Fermi contact', &
       & 'Spin-dipole', 'Paramagnetic spin-orbit'])
  spin_operator = any(term_name == &
       & [character(28) :: 'Fermi contact', 'Spin-dipole'])
  calc_H1_Ha = n_spin == 2 .and. .not. imaginary_pert

  ! Allocate the main work variables
  if (calc_full_tensor) then
     ! See comments above for matrix_BC dimensions
     if (term_name == 'Diamagnetic magnetizability' .and. .not. no_giao) then
        ! Spin-polarized diamagnetic magnetizability is a special case.
        mem_matrix_BC = int(mxld,8)*int(mxcol,8)*int(n_spin*9,8)*int(8,8)
        allocate(matrix_BC(mxld, mxcol, n_spin*9), stat=info)
     else
        mem_matrix_BC = int(mxld,8)*int(mxcol,8)*int(9,8)*int(8,8)
        allocate(matrix_BC(mxld, mxcol, 9), stat=info)
     end if
  else
     mem_matrix_BC = int(mxld,8)*int(mxcol,8)*int(n_directions*n_spin,8)* &
          & int(8,8)
     allocate(matrix_BC(mxld, mxcol, n_directions*n_spin), stat=info)
  end if
  call update_when_allocating(info, THIS_SUB//"matrix_BC", mem_matrix_BC)
  if (.not. skip_dfpt) then
     if (skip_last_dfpt) then
        call aims_allocate(Psi1, mxld, mxcol, n_spin, n_directions, n_psi1-1, &
             & THIS_SUB//'Psi1')
        call aims_allocate(H_Psi, 1, mxld, 1, mxcol, 1, n_spin, 1, &
             & n_directions, 2, n_psi1, THIS_SUB//'H_Psi')
     else
        call aims_allocate(Psi1, mxld, mxcol, n_spin, n_directions, n_psi1, &
             & THIS_SUB//'Psi1')
        call aims_allocate(H_Psi, mxld, mxcol, n_spin, n_directions, n_psi1, &
             & THIS_SUB//'H_Psi')
     end if
  else
     call aims_allocate(H_Psi, 1, mxld, 1, mxcol, 1, n_spin, 1, 1, 1, n_psi1, &
          & THIS_SUB//'H_Psi')
  end if

  ! STEP 1 - Calculate the first-order response. The target quantities
  !          are Psi1 and H1_Psi0. The part has a simple structure -
  !          there is one loop over all atoms of interest. All spatial
  !          directions (1, 3, or some other number depending on the
  !          perturbation) are simultaneously processed. Most
  !          important pieces here are calls to 'integrate' to compute
  !          the non-self-consistent first-order matrix elements and
  !          the call to 'do_DFPT', which performs a self-consistent
  !          DFPT cycle.
  if (.not. skip_dfpt) call localorb_multi( &
       & 'Calculating first-order response (DFPT)', &
       & '---------------------------------------', &
       & format='(2x, a)')
  calc_psi1: do i_atom = 1, n_psi1
     if (skip_dfpt) exit calc_psi1
     ! STEP 1a - Compute the non-self-consistent matrix elements
     !           that enter as input for DFPT.
     if (.not. any(term_name == [character(28) :: 'Paramagnetic shielding', &
          & 'Paramagnetic magnetizability'])) &
          & call localorb_multi(trim(c_atoms(i_atom)), format='(2x, a)')
     select case(term_name)
     case('Fermi contact')
        if (flag_rel == REL_none) then
           call integrate_FC(active_nuclei(i_atom), matrix_BC, timer)
           matrix_BC(:,:,1) = 4*pi*alpha2/3*matrix_BC(:,:,1)
        else
           call integrate(INT_FERMI_CONTACT_SR, matrix_BC, timer=timer, &
                & points_for_distance=coords(:,active_nuclei(i_atom)))
           matrix_BC(:,:,1) = -2*alpha2/3*matrix_BC(:,:,1)
        end if
        if (n_spin == 2) matrix_BC(:,:,2) = -matrix_BC(:,:,1)
     case('Paramagnetic spin-orbit')
        if (flag_rel == REL_none) then
           call integrate(INT_PARAMAGNETIC_SO, matrix_BC, timer=timer, &
                & points_for_distance=coords(:,active_nuclei(i_atom)))
        else
           call integrate(INT_PARAMAGNETIC_SO_SR, matrix_BC, &
                & timer=timer, &
                & points_for_distance=coords(:,active_nuclei(i_atom)))
        end if
        call antisymmetrize(matrix_BC, 3)
        ! If n_spin == 2: 1,3,5 spin-up, 2,4,6 spin-down
        if (n_spin == 2) call rearrange_spins(matrix_BC, 1)
        matrix_BC = -alpha2*matrix_BC ! -i*alpha^2
     case('Spin-dipole')
        if (flag_rel == REL_none) then
           if (calc_spin_dipole_diag) then
              call integrate(INT_SPIN_DIPOLE_XX_YY_ZZ, matrix_BC, &
                   & timer=timer, &
                   & points_for_distance=coords(:,active_nuclei(i_atom)))
           else
              call integrate(INT_SPIN_DIPOLE_XY_XZ_YZ, matrix_BC, &
                   & timer=timer, &
                   & points_for_distance=coords(:,active_nuclei(i_atom)))
           end if
        else
           if (calc_spin_dipole_diag) then
              call integrate(INT_SPIN_DIPOLE_XX_YY_ZZ, matrix_BC, &
                   & timer=timer, &
                   & points_for_distance=coords(:,active_nuclei(i_atom)))
           else
              call integrate(INT_SPIN_DIPOLE_XY_XZ_YZ, matrix_BC, &
                   & timer=timer, &
                   & points_for_distance=coords(:,active_nuclei(i_atom)))
           end if
        end if
        if (n_spin == 2) then
           call rearrange_spins(matrix_BC, 1)
           matrix_BC(:,:,2:6:2) = -matrix_BC(:,:,2:6:2)
        end if
        matrix_BC = alpha2/2*matrix_BC ! alpha^2/2
     case('Paramagnetic shielding', 'Paramagnetic magnetizability')
        if (no_giao) then
           call integrate(INT_ORBITAL_ZEEMAN, matrix_BC, timer=timer, &
                & points_for_distance=gauge_origin)
           call antisymmetrize(matrix_BC, 3)
           if (n_spin == 2) call rearrange_spins(matrix_BC, 1)
        else
           call integrate(INT_GIAO_PARAMAGNETIC, matrix_BC, timer=timer)
           if (n_spin == 2) call rearrange_spins(matrix_BC, n_spin)
        end if
        matrix_BC = -0.5d0*matrix_BC ! -i/2
     end select
     ! Done with the main integrals. Next compute additional
     ! quantities such as the first-order overlap matrix if necessary.
     if (include_overlap1) then
        call aims_allocate(overlap1, mxld, mxcol, 3, THIS_SUB//'overlap1')
        call integrate(INT_GIAO_OVERLAP1, overlap1, timer=timer)
        overlap1(:,:,:3) = 0.5d0*overlap1(:,:,:3) ! i/2
     end if
     if (include_epsilon1) then
        call aims_allocate(epsilon1, mxld, mxcol, n_spin, n_directions, &
             & THIS_SUB//'epsilon1')
        do i_dir = 1, n_directions
           do i_spin = 1, n_spin
              ! epsilon1_mn = -(E_m+E_n)/2 <Psi0|S1|Psi0>_mn
              call start_wall(walltime_mul)
              call mul(overlap1(:,:,i_dir), eigenvec(:,:,i_spin), &
                   & H_Psi(:,:,i_spin,i_dir,1), N=n_occ_states(i_spin))
              call mul(eigenvec(:,:,i_spin), H_Psi(:,:,i_spin,i_dir,1), &
                   & epsilon1(:,:,i_spin,i_dir), transa='t', &
                   & N=n_occ_states(i_spin), M=n_occ_states(i_spin))
              call stop_wall(walltime_mul)
              i_row = my_max_occ_row(i_spin)
              do i_col = 1, my_max_occ_col(i_spin)
                 epsilon1(:i_row,i_col,i_spin,i_dir) = &
                      & -epsilon1(:i_row,i_col,i_spin,i_dir)* &
                      & (KS_eigenvalue(my_row(:i_row),i_spin,1) + &
                      & KS_eigenvalue(my_col(i_col),i_spin,1))/2
              end do
              ! epsilon1 += <Psi0|H1|Psi0>
              ! Note: it is assumed that matrix_BC was determined
              ! above by the call with INT_GIAO_PARAMAGNETIC.
              call start_wall(walltime_mul)
              call mul(matrix_BC(:,:,i_spin+n_spin*(i_dir-1)), &
                   & eigenvec(:,:,i_spin), H_Psi(:,:,i_spin,i_dir,1), &
                   & N=n_occ_states(i_spin))
              call mul(eigenvec(:,:,i_spin), H_Psi(:,:,i_spin,i_dir,1), &
                   & epsilon1(:,:,i_spin,i_dir), transa='t', &
                   & N=n_occ_states(i_spin), M=n_occ_states(i_spin), beta=1d0)
              call stop_wall(walltime_mul)
           end do
        end do
     end if
     ! STEP 1b - Compute H1_nsc|psi0> for STEP 2.
     if (save_H_Psi .and. i_atom >= lbound(H_Psi,5)) then
        ! Not computed for the first atom. See comments about
        ! save_H_Psi above.
        do i_dir = 1, n_directions
           do i_spin = 1, n_spin
              call start_wall(walltime_mul)
              call mul(matrix_BC(:,:,i_spin+n_spin*(i_dir-1)), &
                   & eigenvec(:,:,i_spin), H_Psi(:,:,i_spin,i_dir,i_atom), &
                   & N=n_occ_states(i_spin))
              call stop_wall(walltime_mul)
           end do
        end do
     end if
     ! STEP 1c - Compute first-order wavefunctions from DFPT.
     if (i_atom > size(Psi1,5)) then
        ! Calculation of the last Psi1 skipped in this case.
        call localorb_multi('Psi1 not needed', format='(6x, a)')
        exit calc_psi1
     end if
     ! Measure memory usage only for the first nucleus in the list. We
     ! don't want to contaminate max_dfpt_memory with allocations from
     ! elsewhere such as the integration routines above.
     if (i_atom == 1) call push_current_memory()
     call do_DFPT(matrix_BC, overlap1, imaginary_pert, n_directions, &
          & calc_H1_Ha, include_overlap1, spin_operator, Psi1(:,:,:,:,i_atom))
     ! max_dfpt_memory was initialized in MR_initialize
     if (i_atom == 1) &
          & max_dfpt_memory = max(pop_memory_estimate(),max_dfpt_memory)
     ! No loop over atoms for these perturbations
     if (any(term_name == [character(28) :: 'Paramagnetic shielding', &
          & 'Paramagnetic magnetizability'])) exit calc_psi1
  end do calc_psi1
  ! We don't need any of the DFPT work variables anymore.
  call cleanup_DFPT()
  if (.not. skip_dfpt) &
       & call localorb_multi('', 'Done with first-order response', '', &
       & format='(2x, a)')

  ! STEP 2 - Calculate the second order response. Structure of this
  !          part is the following. Two outer loops are over pairs of
  !          atoms of interest. The 'if_dia' conditional determines
  !          whether we are evaluating a diamagnetic (<Psi0|H2|Psi0>)
  !          or a paramagnetic (2Re<Psi1|H1|Psi0>) response.
  call localorb_multi( &
       & 'Calculating second-order response', &
       & '---------------------------------', &
       & format='(2x, a)')
  atom_1: do i_atom = 1, max(1,size(active_nuclei))
     atom_2: do j_atom = i_atom, max(1,size(active_nuclei))
        ! When calculating J-couplings, j_atom runs over i_atom+1,
        ! size(active_nuclei). When calculating shieldings, we only
        ! consider j_atom == i_atom (i_atom == j_atom == 1 for
        ! magnetizability) and exit the loop otherwise.
        select case(term_name)
        case('Paramagnetic magnetizability', 'Diamagnetic magnetizability')
           if(i_atom /= 1 .or. j_atom /= 1) cycle atom_1
        case('Paramagnetic shielding', 'Diamagnetic shielding', &
             & 'Electric field gradient')
           if (j_atom /= i_atom) cycle atom_2
           call localorb_multi(trim(c_atoms(i_atom)), format='(4x, a)')
        case default
           if (j_atom == i_atom) cycle atom_2
           if (term_name == 'Diamagnetic spin-orbit') &
                & call localorb_multi(trim(c_atoms(i_atom))//' - '// &
                & trim(c_atoms(j_atom)), format='(4x, a)')
        end select
        ! STEP 2a - Compute the |H2|Psi0> terms and the respective
        !           responses by taking <Psi0|H2|Psi0>.
        if_dia: if (skip_dfpt) then
           select case(term_name)
           case('Diamagnetic spin-orbit')
              if (flag_rel == REL_none) then
                 call integrate(INT_DIAMAGNETIC_SO, matrix_BC, &
                      & timer=timer, points_for_distance=[ &
                      & coords(:,active_nuclei(i_atom)), &
                      & coords(:,active_nuclei(j_atom))])
              else
                 call integrate(INT_DIAMAGNETIC_SO_SR, matrix_BC, &
                      & timer=timer, points_for_distance=[ &
                      & coords(:,active_nuclei(i_atom)), &
                      & coords(:,active_nuclei(j_atom))])
              end if
              matrix_BC = alpha4*matrix_BC! alpha^4
           case('Diamagnetic shielding')
              if (no_giao) then
                 call integrate(INT_DIAMAGNETIC_SHIELDING_STD, matrix_BC, &
                      & timer=timer, points_for_distance=[gauge_origin, &
                      & coords(:,active_nuclei(i_atom))])
              else
                 call integrate(INT_DIAMAGNETIC_SHIELDING, matrix_BC, &
                      & timer=timer, points_for_distance= &
                      & [coords(:,active_nuclei(i_atom))])
              end if
              matrix_BC = alpha2/2*matrix_BC ! alpha^2/2
           case('Electric field gradient')
              call integrate(INT_ELECTRIC_FIELD_GRADIENT, matrix_BC, &
                   & timer=timer, &
                   & points_for_distance=[coords(:,active_nuclei(i_atom))])
              !ToDO check the constants!
           case('Diamagnetic magnetizability')
              if (no_giao) then
                 call integrate(INT_DIAMAGNETIC_MAGNETIZABILITY_STD, &
                      & matrix_BC, timer=timer, &
                      & points_for_distance=[gauge_origin])
              else
                 call integrate(INT_DIAMAGNETIC_MAGNETIZABILITY, &
                      & matrix_BC, timer=timer)
              end if
              matrix_BC = matrix_BC/4
           end select
           if (include_overlap2 .and. .not. allocated(overlap2)) then
              call aims_allocate(overlap2, mxld, mxcol, n_tensor_components, &
                   & THIS_SUB//'overlap2')
              call integrate(INT_GIAO_OVERLAP2, overlap2, timer=timer)
              overlap2 = overlap2/4
           end if
           ! All diamagnetic integrals have been calculated. Now take
           ! the sums over occupied states and put the results into
           ! respective response tensors.
           tmp_tensor = 0d0
           do i_dir = 1, 3
              do j_dir = 1, 3
                 do i_spin = 1, n_spin
                    if (term_name == 'Diamagnetic magnetizability' .and. &
                         & .not. no_giao) then
                       ! Diamagnetic magnetizability is the only case
                       ! here whose matrix elements could depend on
                       ! spin.
                       shift = (i_spin-1)*n_tensor_components
                    else
                       shift = 0
                    end if
                    ! If n_tensor_components == 3, matrix elements
                    ! corresponding to the 3 directions are stored as
                    ! follows:
                    ! matrix_BC(:,:,1) - direction 1
                    ! matrix_BC(:,:,2) - direction 2
                    ! matrix_BC(:,:,3) - direction 3
                    if (n_tensor_components == 3) then
                       if (j_dir == i_dir) then
                          call start_wall(walltime_mul)
                          call mul(matrix_BC(:,:,j_dir+shift), &
                               & eigenvec(:,:,i_spin), H_Psi(:,:,1,1,1), &
                               & N=n_occ_states(i_spin))
                          call stop_wall(walltime_mul)
                       else
                          H_Psi(:,:,1,1,1) = 0d0
                       end if
                    else
                       ! If n_tensor_components == 9, then
                       ! matrix_BC(:,:,1) - 11 or xx
                       ! matrix_BC(:,:,2) - 21 or yx
                       ! matrix_BC(:,:,3) - 31 or zx
                       ! matrix_BC(:,:,4) - 12 or xy
                       ! matrix_BC(:,:,5) - 22 or yy
                       ! matrix_BC(:,:,6) - 32 or zy
                       ! matrix_BC(:,:,7) - 13 or xz
                       ! matrix_BC(:,:,8) - 23 or yz
                       ! matrix_BC(:,:,9) - 33 or zz
                       call start_wall(walltime_mul)
                       call mul(matrix_BC(:,:,i_dir+3*(j_dir-1)+shift), &
                            & eigenvec(:,:,i_spin), H_Psi(:,:,1,1,1), &
                            & N=n_occ_states(i_spin))
                       call stop_wall(walltime_mul)
                    end if
                    ! Add the <Psi0|H2|Psi0> term to the response tensor.
                    tmp_tensor(i_dir,j_dir) = tmp_tensor(i_dir,j_dir) + &
                         & sum(eigenvec(:,:my_max_occ_col(i_spin),i_spin)* &
                         & H_Psi(:,:my_max_occ_col(i_spin),1,1,1))
                    if (include_overlap2) then
                       ! If required, add the -<Psi0|S2|Psi0>epsilon term.
                       if (n_tensor_components == 3) then
                          if (j_dir == i_dir) then
                             call start_wall(walltime_mul)
                             call mul(overlap2(:,:,j_dir), &
                                  & eigenvec(:,:,i_spin), H_Psi(:,:,1,1,1), &
                                  & N=n_occ_states(i_spin))
                             call stop_wall(walltime_mul)
                          else
                             H_Psi(:,:,1,1,1) = 0d0
                          end if
                       else
                          call start_wall(walltime_mul)
                          call mul(overlap2(:,:,i_dir+3*(j_dir-1)), &
                               & eigenvec(:,:,i_spin), H_Psi(:,:,1,1,1), &
                               & N=n_occ_states(i_spin))
                          call stop_wall(walltime_mul)
                       end if
                       do i_col = 1, my_max_occ_col(i_spin)
                          H_Psi(:,i_col,1,1,1) = &
                               & -H_Psi(:,i_col,1,1,1)* &
                               & KS_eigenvalue(my_col(i_col),i_spin,1)
                       end do
                       tmp_tensor(i_dir,j_dir) = tmp_tensor(i_dir,j_dir) + &
                            & sum(eigenvec(:,:my_max_occ_col(i_spin),i_spin)*&
                            & H_Psi(:,:my_max_occ_col(i_spin),1,1,1))
                    end if
                 end do
              end do
           end do
           tmp_tensor = (3-n_spin)*tmp_tensor ! Spin degeneracy
           select case(term_name)
           case('Diamagnetic spin-orbit')
              DSO_tensor(:,:,i_atom,j_atom) = tmp_tensor
           case('Diamagnetic shielding')
              shield_dia_tensor(:,:,i_atom) = tmp_tensor
           case('Diamagnetic magnetizability')
              mag_dia_tensor = tmp_tensor
           case('Electric field gradient')
              EFG_tensor(:,:,i_atom) = tmp_tensor
           end select
        else ! if_dia
           ! STEP 2b - Compute the <Psi1|H1|Psi0> terms and the
           !           respective responses. If |Psi> was calculated
           !           using the same perturbation as H1 (e.g.,
           !           paramagnetic spin-orbit term of J-coupling), no
           !           extra integrals are computed in this section.
           select case(term_name)
           case('Paramagnetic shielding')
              ! For paramagnetic shielding, calculate the
              ! H1|Psi0> vectors first.
              if (flag_rel == REL_none) then
                 call integrate(INT_PARAMAGNETIC_SO, matrix_BC, &
                      & timer=timer, &
                      & points_for_distance=coords(:,active_nuclei(j_atom)))
              else
                 call integrate(INT_PARAMAGNETIC_SO_SR, matrix_BC, &
                      & timer=timer, &
                      & points_for_distance=coords(:,active_nuclei(j_atom)))
              end if
              call antisymmetrize(matrix_BC, 3)
              matrix_BC(:,:,1:3) = -alpha2*matrix_BC(:,:,1:3) ! -i*alpha^2
              do i_dir = 1, 3
                 do i_spin = 1, n_spin
                    call start_wall(walltime_mul)
                    call mul(matrix_BC(:,:,i_dir), eigenvec(:,:,i_spin), &
                         & H_Psi(:,:,i_spin,i_dir,1), N=n_occ_states(i_spin))
                    call stop_wall(walltime_mul)
                 end do
              end do
              call update_para_response(1, 1, shield_para_tensor(:,:,i_atom))
           case('Paramagnetic spin-orbit')
              call update_para_response(i_atom, j_atom, &
                   & PSO_tensor(:,:,i_atom,j_atom))
           case('Fermi contact')
              call update_para_response(i_atom, j_atom, &
                   & FC_tensor(:1,:1,i_atom,j_atom))
              FC_tensor(2,2,i_atom,j_atom) = FC_tensor(1,1,i_atom,j_atom)
              FC_tensor(3,3,i_atom,j_atom) = FC_tensor(1,1,i_atom,j_atom)
           case('Spin-dipole')
              tmp_tensor = 0d0
              call update_para_response(i_atom, j_atom, tmp_tensor)
              ! If calc_spin_dipole_diag then i_dir = [xx,yy,zz]
              ! else i_dir =[xy,xz,yz].
              do i_dir = 1, 3
                 if (calc_spin_dipole_diag) then
                    spin_dipole_tensor(i_dir,i_dir,i_atom,j_atom) = &
                         & tmp_tensor(i_dir,i_dir)
                 else
                    select case(i_dir)
                    case(1)
                       spin_dipole_tensor(1,1,i_atom,j_atom) = &
                            & spin_dipole_tensor(1,1,i_atom,j_atom) + &
                            & tmp_tensor(i_dir,i_dir)
                       spin_dipole_tensor(2,2,i_atom,j_atom) = &
                            & spin_dipole_tensor(2,2,i_atom,j_atom) + &
                            & tmp_tensor(i_dir,i_dir)
                    case(2)
                       spin_dipole_tensor(1,1,i_atom,j_atom) = &
                            & spin_dipole_tensor(1,1,i_atom,j_atom) + &
                            & tmp_tensor(i_dir,i_dir)
                       spin_dipole_tensor(3,3,i_atom,j_atom) = &
                            & spin_dipole_tensor(3,3,i_atom,j_atom) + &
                            & tmp_tensor(i_dir,i_dir)
                    case(3)
                       spin_dipole_tensor(2,2,i_atom,j_atom) = &
                            & spin_dipole_tensor(2,2,i_atom,j_atom) + &
                            & tmp_tensor(i_dir,i_dir)
                       spin_dipole_tensor(3,3,i_atom,j_atom) = &
                            & spin_dipole_tensor(3,3,i_atom,j_atom) + &
                            & tmp_tensor(i_dir,i_dir)
                    end select
                 end if
              end do
           case('Paramagnetic magnetizability')
              ! First add the 2Re<Psi1|H1|Psi0> term to the response tensor.
              call update_para_response(1, 1, mag_para_tensor)
              if (.not. no_giao) then
                 ! Next add the -2Re<Psi1|S1|Psi0>epsilon term.
                 do i_dir = 1, n_directions
                    do i_spin = 1, n_spin
                       call start_wall(walltime_mul)
                       call mul(overlap1(:,:,i_dir), eigenvec(:,:,i_spin), &
                            & H_Psi(:,:,i_spin,i_dir,1), &
                            & N=n_occ_states(i_spin))
                       call stop_wall(walltime_mul)
                       do i_col = 1, my_max_occ_col(i_spin)
                          H_Psi(:,i_col,i_spin,i_dir,1) = &
                               & -H_Psi(:,i_col,i_spin,i_dir,1)* &
                               & KS_eigenvalue(my_col(i_col),i_spin,1)
                       end do
                    end do
                 end do
                 call update_para_response(1, 1, mag_para_tensor)
                 ! Next add the -<Psi0|S1|Psi0>epsilon1 term.
                 do i_dir = 1, 3
                    do i_spin = 1, n_spin
                       call start_wall(walltime_mul)
                       call mul(overlap1(:,:,i_dir), eigenvec(:,:,i_spin), &
                            & H_Psi(:,:,i_spin,i_dir,1), &
                            & N=n_occ_states(i_spin))
                       call mul(eigenvec(:,:,i_spin), &
                            & H_Psi(:,:,i_spin,i_dir,1), matrix_BC(:,:,1), &
                            & transa='t', N=n_occ_states(i_spin))
                       call stop_wall(walltime_mul)
                       i_row = my_max_occ_row(i_spin)
                       i_col = my_max_occ_col(i_spin)
                       do j_dir = 1, 3
                          mag_para_tensor(i_dir,j_dir) = &
                               & mag_para_tensor(i_dir,j_dir) - &
                               & (3-n_spin)* & ! Spin degeneracy
                               & sum(matrix_BC(:i_row,:i_col,1)* &
                               & epsilon1(:i_row,:i_col,i_spin,j_dir))
                       end do
                    end do
                 end do
              end if
           end select
        end if if_dia
     end do atom_2
  end do atom_1
  call localorb_multi('', 'Done with second order response', &
       & '', format='(2x, a)')

  ! STEP 3 - deallocations
  call safe_deallocate(matrix_BC)
  call safe_deallocate(H_Psi)
  call safe_deallocate(Psi1)
  call safe_deallocate(overlap1)
  call safe_deallocate(overlap2)
  call safe_deallocate(epsilon1)

contains
  ! Change the ordering of indices: (:,:,3,2) -> (:,:,2,3), where 2
  ! and 3 correspond to spins and directions.
  subroutine rearrange_spins(matrix_BC, n_spins)
    real(dp), intent(in out) :: matrix_BC(:,:,:)
    integer, intent(in) :: n_spins
    real(dp), allocatable :: tmp(:,:)
    if (n_spins == 1) then
       ! Copy the values from spin channel 1 to spin channel 2.
       matrix_BC(:,:,6) = matrix_BC(:,:,3)
       matrix_BC(:,:,5) = matrix_BC(:,:,3)
       matrix_BC(:,:,4) = matrix_BC(:,:,2)
       matrix_BC(:,:,3) = matrix_BC(:,:,2)
       matrix_BC(:,:,2) = matrix_BC(:,:,1)
    else
       ! If n_spins == 2, the values for both spin channels were
       ! explicitly calculated.
       call aims_allocate(tmp, mxld,mxcol, THIS_SUB//'rearrange_spins::tmp')
       tmp = matrix_BC(:,:,2)
       matrix_BC(:,:,2) = matrix_BC(:,:,4)
       matrix_BC(:,:,4) = matrix_BC(:,:,5)
       matrix_BC(:,:,5) = matrix_BC(:,:,3)
       matrix_BC(:,:,3) = tmp
       call safe_deallocate(tmp)
    end if
  end subroutine rearrange_spins

  subroutine update_para_response(i_atom, j_atom, X)
    integer, intent(in) :: i_atom, j_atom
    real(dp), intent(in out) :: X(:,:)
    do i_dir = 1, size(X,1)
       do j_dir = 1, size(X,2)
          do i_spin = 1, n_spin
             ! Prefactor is spin degeneracy
             X(i_dir,j_dir) = X(i_dir,j_dir) + (3-n_spin)* &
                  & 2*sum(Psi1(:,:my_max_occ_col(i_spin),i_spin,i_dir,i_atom)* &
                  & H_Psi(:,:my_max_occ_col(i_spin),i_spin,j_dir,j_atom))
          end do
       end do
    end do
  end subroutine update_para_response
end subroutine MR_core
