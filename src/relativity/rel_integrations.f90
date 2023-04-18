
! This file contains subroutines relating to fully-relativistic integrations.
! -- Rundong Zhao, Sep. 2018 @ Durham, NC
 subroutine prune_basis_p2_rel ( dist_tab_sq, n_compute, i_basis_small,  n_basis_list, n_atom_list, inv_list )

!  PURPOSE
!  Reduces full basis to relevant functions only, for a given set of integration points
!
  use dimensions
  use basis
  use grids
  use geometry
  use pbc_lists
  implicit none

!  ARGUMENTS
  integer :: n_basis_list, n_atom_list
  integer :: inv_list(n_centers)
  real*8  :: dist_tab_sq(n_atom_list)
  integer :: n_compute
  integer :: i_basis_small(n_basis_list)

!  INPUTS
!    o n_basis_list -- number of basis functions (including periodic mirror images)
!    o n_atom_list -- number of atoms (including periodic mirror images)
!    o inv_list -- citing: atom center -> relative position in distance list
!    o dist_tab_sq -- (distance to atoms)**2
!   
!  OUTPUT
!    o n_compute_a == n_compute_c -- number of relevant basis functions
!    o i_basis_small -- list of relevant basis functions

  !     counters
  integer :: i_basis_1
  integer :: i_compute_center
  integer :: i_compute_1

  logical :: flag = .false.

  ! tabulate total wave function value for each basis function

  if (n_compute.eq.0) then
     ! this is the first point of the batch - simply list all non-zero wave functions

     i_compute_center = 0

     do i_basis_1 = 1, n_basis_list, 1

         !write(6,"('i_basis_1=',i3,3x,'dist_tab_sq=',f12.6,3x,'outer_radius_sq=',f12.6)") &
         ! i_basis_1, dist_tab_sq(inv_list(Cbasis_to_center_s(i_basis_1))),outer_radius_sq(basis_small_fn(Cbasis_to_basis_s(i_basis_1)))
        if (dist_tab_sq(inv_list(Cbasis_to_center_s(i_basis_1))) .le. &
             outer_radius_sq(basis_small_fn(Cbasis_to_basis_s(i_basis_1)))) then

           !  nonzero basis function - use it ... 
           i_compute_center = i_compute_center+1
           i_basis_small(i_compute_center) = i_basis_1
                      ! write(6,"('i_compute_center=',i3,3x,'i_basis_small=',i3)")i_compute_center, i_basis_1
        end if

     enddo

     n_compute = i_compute_center

  else
     ! this is not the first integration point of the batch - check whether non-zero
     ! wave functions are already there, else add them to the list

     i_compute_center = 0

     do i_basis_1 = 1, n_basis_list, 1

        if (i_basis_small(i_compute_center+1).eq.i_basis_1) then
           ! basis function already in list

           i_compute_center = i_compute_center+1

        else if (dist_tab_sq(inv_list(Cbasis_to_center_s(i_basis_1))) .le. &
             outer_radius_sq(basis_small_fn(Cbasis_to_basis_s(i_basis_1)))) then
           ! basis function not in list - add it to the list of nonzero functions

           i_compute_center = i_compute_center+1

           do i_compute_1 = n_compute, i_compute_center, -1
              i_basis_small(i_compute_1+1) = i_basis_small(i_compute_1)
           enddo

           i_basis_small(i_compute_center) = i_basis_1

           n_compute = n_compute+1

        end if

     enddo

  end if

 end subroutine prune_basis_p2_rel


 subroutine collect_batch_centers_rel ( n_compute, i_basis, n_basis_list, n_atom_list, inv_list, n_batch_centers, batch_center )

!  PURPOSE
!     For a given batch of integration points, we tabulate all integration centers
!     that are ever relevant in current batch, to remove any O(N^2) steps later
!     in prune_radial_basis ...

  use dimensions
  use basis
  use grids
  use geometry
  use pbc_lists
  use runtime_choices, only : flag_rel, REL_x2c, REL_4c_dks
  implicit none

!  ARGUMENTS
  integer :: n_basis_list, n_atom_list
  integer :: inv_list(n_centers)
  integer :: n_compute
  integer :: i_basis(n_compute)
  integer :: n_batch_centers
  integer :: batch_center(n_atom_list)

!  INPUTS
!   o n_basis_list -- total number of basis functions
!   o n_atom_list -- number of relevant atoms
!   o inv_list -- inverse citing list of atoms
!   o n_compute -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
! 
!  OUTPUT
!   o n_batch_centers -- number of batch centers 
!   o batch_center -- list of batch centers
 
!  local variables
! counters

  integer :: i_compute

  integer :: i_basis_1
  integer :: i_center

  integer :: i_batch_center
  integer :: i_batch_center_2

  logical :: found

  n_batch_centers = 0

  do i_compute = 1, n_compute, 1
    i_basis_1 = i_basis(i_compute)

    ! number of center associated with current basis function
    i_center = inv_list(Cbasis_to_center_s(i_basis_1))

!test
!       write(6,*) i_compute, i_center
!test end

    ! check if we already know of this center, else add it to the list

    found = .false.
    i_batch_center = 0
    do while ( (.not.found) .and. (i_batch_center.lt.n_batch_centers) )
      i_batch_center = i_batch_center+1

      if (batch_center(i_batch_center).eq.i_center) then
        ! atom already in list - no need to search further
        found = .true.

      else if (batch_center(i_batch_center).gt.i_center) then
        ! atom was apparently not in the list, since we create a list in
        ! order of increasing center index here ...

        n_batch_centers = n_batch_centers+1
        do i_batch_center_2 = n_batch_centers, i_batch_center+1, -1
          batch_center(i_batch_center_2) = batch_center(i_batch_center_2-1)
        enddo
        batch_center(i_batch_center) = i_center

        found = .true.

      end if

    enddo

    if (.not.found) then
      ! Atom was not yet in list and must have a higher index than any atom already in list

        n_batch_centers = n_batch_centers+1
        batch_center(n_batch_centers) = i_center
     
    end if

  enddo

 end subroutine collect_batch_centers_rel


 subroutine prune_radial_basis_rel ( dim_atoms, dim_fns, dist_tab_sq_s, dist_tab_s, dir_tab_s, &
  n_compute_atoms_s, atom_index_s, atom_index_inv_s, n_compute_fns_s, i_basis_fns_s, i_basis_fns_inv_s, &
  i_atom_fns_s, spline_array_start_s, spline_array_end_s, n_atom_list, atom_list, n_compute_small, i_basis_small, &
  n_batch_centers_s, batch_center_s, one_over_dist_tab_s, rad_index_small, wave_index_small, l_index_small, &
  l_count_small, fn_atom_small, n_zero_compute_s, zero_index_point_s )

!  PURPOSE
!
!     For a given integration point, determines and indexes:
!     *  those atoms from which at least one (radial) wave function is needed at this point 
!     *  those radial functions which are non-zero at this point (# n_compute_fns_s)
!     *  among _all_ wave functions relevant in current batch of points (# n_compute), those 
!         which are zero at current point. (# n_zero_compute)
!
!--------------------------------------------------------------------------------------------------------
! (Rundong Zhao) For a fully-relativistic basis set, subroutine prune_radial_basis_rel 
! (the counterpart of prune_radial_basis_p2) is used for generating small component basis.
! In this case, the large comp. indices are generated in subroutine prune_radial_basis_p2.
!
! The relativistic SCALAR basis set involves large comp. and small comp. parts, the number of which 
! are DIFFERENT (small comp. usually contains more basis, as interpreted below).
!
! We will firstly obtain the SCALAR integrations (for S, T, V, and W matrices) with SCALAR basis; then 
! assemble the SPINOR integrations with the SCALAR integrations. -- Will be done in other subroutines.
!-------------------------------------------------------------------------------------------------------- 

  use dimensions
  use basis
  use grids
  use geometry
  use pbc_lists
  use localorb_io, only: use_unit
  implicit none
    
!  ARGUMENTS

  integer,intent(in) :: dim_atoms  ! dimensions that bound n_compute_atoms_s
  integer,intent(in) :: dim_fns ! dimensions that bound n_compute_fns_s
  integer,intent(in):: n_atom_list ! total number of atom centers
  integer,intent(in):: atom_list(n_atom_list) ! list of atom centers

  real*8,intent(inout) :: dist_tab_sq_s(n_atom_list) ! (distance to atom centers)**2
  real*8,intent(inout) :: dir_tab_s(3,n_atom_list) ! input: direction to atoms; output: direction to relevant atoms (normalized)
  real*8,intent(out) :: dist_tab_s(n_atom_list) ! distance to relevant atoms

  integer,intent(in) :: n_batch_centers_s ! number of batch centers
  integer,intent(in) :: batch_center_s(n_atom_list) ! list of batch centers

  integer,intent(out) :: n_compute_atoms_s ! number of relevant atoms
  integer,intent(out) :: atom_index_s(n_atom_list) !  list of relevant atoms
  integer,intent(out) :: atom_index_inv_s(n_centers) !  inverse citing list of atom_index_s

! n_compute_fns_s: total number of radial functions relevant at these points
! note that each atom (not species) contributes to this number
! separately, so that in most cases n_basis_fns < n_compute_fns_s < n_basis
  integer,intent(out) :: n_compute_fns_s
  integer,intent(out) :: i_basis_fns_s(n_basis_fns*n_atom_list) ! list of relevant radial basis functions
  integer,intent(out) :: i_basis_fns_inv_s(n_basis_fns,n_centers) ! inverse citing list of i_basis_fns_s
! i_atom_fns_s -- in the list of non-zero atoms at the current point, the index
! of the atom at which radial function i_compute_fn is centered
  integer,intent(out) :: i_atom_fns_s(n_basis_fns*n_atom_list)

  integer,intent(out) :: spline_array_start_s(n_atom_list) ! starting point of splined data information
  integer,intent(out) :: spline_array_end_s(n_atom_list) ! ending point of splined data information

  integer,intent(in) :: n_compute_small ! number of relevant basis functions ( For indexing of full basis, to avoid any reindexing later)
  integer,intent(in) :: i_basis_small(n_compute_small) ! list of relevant basis functions ( For indexing of full basis, to avoid any reindexing later)

  real*8,intent(out)  :: one_over_dist_tab_s(n_atom_list) ! 1/(distance to atoms)

! rad_index_small -- in the list of n_compute_fns_s radial functions, 
! the END of the radial function list associated with atom i_compute_atom
  integer,intent(out) :: rad_index_small(dim_atoms)
! wave_index -- in the list of non-zero wave functions, the START of
! those wave functions associated with current radial function i_compute_fn
  integer,intent(out) :: wave_index_small(dim_fns)
! l_index -- in the list of ylm functions, the START of
! those ylm functions associated with current radial function i_compute_fn
  integer,intent(out) :: l_index_small(dim_fns)
! l_count -- in the list of ylm functions, the NUMBER of
! ylm functions associated with current radial function i_compute_fn, MINUS 1
  integer,intent(out) :: l_count_small(dim_fns)
  integer,intent(out) :: fn_atom_small(dim_fns)
  
  integer,intent(out) :: n_zero_compute_s ! number of  known zero basis functions at current point
  integer,intent(out) :: zero_index_point_s(n_compute_small) ! list of  known zero basis functions at current point

  ! local variables
  integer :: i_offset_spl
  integer :: l_aux, k_aux

  ! counters

  integer :: i_batch_center
  integer :: i_basis_1
  integer :: i_compute
  integer :: i_fn
  integer :: i_lm
  integer :: i_center, i_center_L
  integer :: i_atom_compute
  integer :: i_spline


  !     check for atoms which contribute basis functions at the current integration point at all
  !     VB: This is now reduced to the list of only the atoms that are active in current batch of
  !         integration points. Thus goes the last O(N^2) part of the integration.
  do i_batch_center = 1, n_batch_centers_s, 1

     i_center_L = batch_center_s(i_batch_center)

!test
     if (.not.(dist_tab_sq_s(i_center_L).gt.0.d0)) then
       write(use_unit,*) "small comp. zero: ", i_center_L, dist_tab_sq_s(i_center_L)
     end if
!test end

     i_center = atom_list(i_center_L)

     ! if this atom is relevant at this point
     if ( dist_tab_sq_s(i_center_L) .le. atom_radius_sq(species_center(i_center)) .and. (dist_tab_sq_s(i_center_L).gt.0.d0) ) then

        n_compute_atoms_s = n_compute_atoms_s + 1
        atom_index_s(n_compute_atoms_s) = i_center
        atom_index_inv_s(i_center) = n_compute_atoms_s


        dist_tab_sq_s(n_compute_atoms_s) = dist_tab_sq_s(i_center_L)
        dist_tab_s(n_compute_atoms_s) = sqrt(dist_tab_sq_s(n_compute_atoms_s))

        one_over_dist_tab_s(n_compute_atoms_s) = 1.d0/dist_tab_s(n_compute_atoms_s) 

        dir_tab_s(1:3,n_compute_atoms_s) = dir_tab_s(1:3,i_center_L) * one_over_dist_tab_s(n_compute_atoms_s) 
     else

        ! Some routines, like small component needs this one.
        atom_index_inv_s(i_center) = 0

     end if

     !     loop over atoms
  enddo


  !     next, check for radial basis functions
  do i_atom_compute = 1, n_compute_atoms_s, 1

     i_center = atom_index_s(i_atom_compute)

     ! This is the index where the splines for the present atom
     ! begin within the condensed array of all splined radial functions:
     i_offset_spl = basis_fn_start_spl(species_center(i_center))

     ! simply run through array of splined functions to determine 
     ! the beginning of the section of the spline array which is
     ! needed from the present atom at the present integration point
     do i_spline = 1, n_basis_fns, 1

        ! original radial function index associated with current spline section
        i_basis_1 = perm_basis_fns_spl(i_spline)

        !          if current basis function is associated with current atom
        if (basis_fn_atom(i_basis_1, center_to_atom(i_center))) then

           if ( dist_tab_sq_s(i_atom_compute).le. outer_radius_sq(i_basis_1)) then
              ! Found the start of the non-zero radial spline functions for current atom;
              ! store for future use.
              spline_array_start_s(i_atom_compute) = (i_spline-1) + 1

              exit

           end if
        end if
     enddo


     ! Now, determine separately the correspondence between radial functions
     ! used at the current point, and radial functions as they were originally stored
     ! in basis_fn(i_basis)
     do i_basis_1 = 1, n_basis_fns  ! n_basis_fns: number of distinct radial basis functions

        !          if current radial function is associated with current atom
        if (basis_fn_atom(i_basis_1, center_to_atom(i_center))) then

           !            if this function is relevant at this point
           if ( dist_tab_sq_s(i_atom_compute) .le. outer_radius_sq(i_basis_1) ) then

              n_compute_fns_s = n_compute_fns_s + 1

              i_basis_fns_s(n_compute_fns_s) = i_basis_1
              i_atom_fns_s(n_compute_fns_s) = i_center
              i_basis_fns_inv_s(i_basis_1,i_center) = n_compute_fns_s

              fn_atom_small(n_compute_fns_s) = i_atom_compute

           end if

           i_offset_spl = i_offset_spl + 1

        end if

        !        loop over radial functions
     enddo
     ! store the last fn associated with current atom in functions list above
     rad_index_small(i_atom_compute) = n_compute_fns_s

     ! The end of the splined array is simply the end of the block where
     ! the splined functions of the current species are stored, which is
     ! what we have after the loop over basis functions ...
     spline_array_end_s(i_atom_compute) = i_offset_spl - 1

     !     loop over atoms
  enddo


  ! Loop over all basis functions u(r)/r * Ylm(Omega) that are active in current 
  ! BATCH of integration points. Sort out those basis functions that are ZERO at
  ! current point, and index those that are non-zero for later use.
  ! NOTE that this sequence relies on the fact that all Ylm functions associated
  ! with a single radial function appear DIRECTLY after one another in the array of
  ! active basis functions.

  ! For relativistic large comp. SCALAR basis: the next 2*basis_l(i_basis(i_compute))+1 elements of wave will always be zero.
  ! (This was done in subroutine prune_radial_basis_p2, not here.)

  ! Take Neon (1s, 2s, 2p) for example:

  ! The nonrel basis set is: 1s, 2s, 2p_y, 2p_z, 2p_x. -- 5 in total.

  ! The rel large comp. SCALAR basis set is: 1s, 2s, 2p_{1/2}y, 2p_{1/2}z, 2p_{1/2}x, 2p_{3/2}y, 2p_{3/2}z, 2p_{3/2}x. -- 8 in total.

  ! However, for small comp.: the FINAL p_{1/2} spinor basis is SOLELY composed of s orbitals, 
  ! while the p_{3/2} spinor basis is SOLELY composed of d orbitals.
  ! So to d_{3/2} and d_{5/2}: composed of p and f, respectively. See Daoling Peng's PhD thesis at Peking University (Appendix B2, pp 114) for details.

  ! Eventually, the rel small comp. SCALAR basis set is: 
  ! 1p_y, 1p_z, 1p_x, 2p_y, 2p_z, 2p_x, 2s, 2d_{xy}, 2d_{yz}, 2d_{z^2}, 2d_{xz}, 2d_{x^2-y^2}. -- 12 in total.

  ! For both large and small comp. basis, the radial function shell is: 1s, 2s, 2p_{1/2}, 2p_{3/2}.


  ! For small component (for large comp. see subroutine prune_radial_basis_p2):
  n_zero_compute_s = 0
  wave_index_small(1:n_compute_fns_s) = 0
  i_compute = 1
  do while (i_compute .le. n_compute_small)
     i_basis_1 = i_basis_small(i_compute)
     ! radial function in list of i_compute_fns that corresponds to 
     ! current basis function
     i_fn = i_basis_fns_inv_s(basis_small_fn(Cbasis_to_basis_s(i_basis_1)), Cbasis_to_center_s(i_basis_1))
     l_aux = basis_small_l(Cbasis_to_basis_s(i_basis_small(i_compute))) 
       ! Note that basis_small_l is already upgraded/degraded for small component part, e.g.:
       ! for p_{1/2}, l_aux = 0; for p_{3/2}, l_aux = 2.
     k_aux = basis_small_k(Cbasis_to_basis_s(i_basis_small(i_compute)))
              ! write(6,"('i_compute=',i3,3x,'i_basis_small=',i3,3x,'l_aux=',i3,3x,'k_aux=',i3,3x,'i_fn=',i3)")i_compute,i_basis_small(i_compute),l_aux,k_aux,i_fn
     if ( i_fn .eq. 0 ) then ! That basis function was listed as zero before and does not exist in list of radial functions.
 
        do i_lm = 0, 2*l_aux, 1
           n_zero_compute_s = n_zero_compute_s + 1
           zero_index_point_s(n_zero_compute_s) = i_compute + i_lm
        enddo

     else

        if (wave_index_small(i_fn).eq.0) then
           ! first time that this function shows up in the list - collect all related indices here.
           wave_index_small(i_fn) = i_compute

           ! starting index for current l in ylm_tab
           l_index_small(i_fn) = l_aux**2 + 1

           ! number of Ylm fns (minus one) for current l
           l_count_small(i_fn) = 2*l_aux

           ! k value of the corresponding basis function
        end if

     end if
     i_compute = i_compute + 2*l_aux + 1
  enddo

 end subroutine prune_radial_basis_rel


 subroutine evaluate_VTW_psi_rel (iop, n_compute, n_compute_atoms, n_compute_fns, l_ylm_max, ylm_tab, one_over_dist_tab, &
   radial_wave, op_times_psi, local_potential_parts, kinetic_wave, rad_index, wave_index, &
   l_index, l_count, fn_atom, n_zero_compute, zero_index_point)

!  PURPOSE
!  Subroutine evaluates V,T times large component basis function, and V times small comp. basis function.

  use dimensions
  use constants
  use basis
  use grids
  use geometry
  use spline
  use runtime_choices
  use pbc_lists

  implicit none

!  ARGUMENTS
  integer,intent(in) :: iop ! Control parameter: iop=1, for V or W matrix integration; iop=2, for T matrix integration; iop=3, one.
                            ! when iop=4, we use kinetic basis to generate W matrix.
  integer,intent(in) :: n_compute_atoms ! Number of relevant atoms
  integer,intent(in) :: n_compute_fns ! Number of relevant radial functions

  integer,intent(in) :: l_ylm_max ! Maximum l index
  real*8,intent(in)  :: ylm_tab ( (l_ylm_max+1)**2, n_compute_atoms ) ! Y_lm functions
  real*8,intent(in)  :: one_over_dist_tab( n_compute_atoms ) ! 1/(distance to atoms)
  integer,intent(in) :: n_compute ! Number of relevant basis functions
  real*8,intent(in)  :: radial_wave(n_compute_fns) ! radial part of basis functions
  real*8,intent(in)  :: local_potential_parts ! total potentials

  real*8,intent(in)  :: kinetic_wave(*) ! kinetic part of the radian basis functions: kinetic_wave(n_compute_fns)

  integer,intent(in) :: rad_index(n_compute_atoms) ! in the list of n_compute_fns radial functions, the END of the radial function list associated with atom i_compute_atom
  integer,intent(in) :: wave_index(n_compute_fns) ! in the list of non-zero wave functions, the START of those wave functions associated with current radial function i_compute_fn
  integer,intent(in) :: l_index(n_compute_fns) ! in the list of ylm functions, the START of those ylm functions associated with current radial function i_compute_fn
  integer,intent(in) :: l_count(n_compute_fns) ! in the list of ylm functions, the NUMBER of ylm functions associated with current radial function i_compute_fn, MINUS 1
  integer,intent(in) :: fn_atom(n_compute_fns)

  integer,intent(in) :: n_zero_compute ! number of  known zero basis functions at current point
  integer,intent(in) :: zero_index_point(n_compute) ! list of  known zero basis functions at current point

  real*8,intent(out) :: op_times_psi(n_compute)
!  OUTPUT
!   V_times_psi -- V times large comp. basis function
!   T_times_psi -- T times large comp. basis function
!   W_times_psi -- V times small comp. basis function

!  local variables
  real*8 op(n_compute_fns) ! The operator. It can be V, T, or W, for different cases.

!  counters
  integer :: i_compute
  integer :: i_compute_point
  integer :: i_compute_fn
  integer :: i_compute_atom
  integer :: index_start
  integer :: index_end

  index_start = 1
  do i_compute_atom = 1, n_compute_atoms, 1

     index_end = rad_index(i_compute_atom)

     !T_plus_V ( index_start:index_end ) = local_potential_parts * radial_wave( index_start:index_end ) + kinetic_wave ( index_start:index_end )
     if(iop.eq.1)then
        op( index_start:index_end ) = local_potential_parts * radial_wave( index_start:index_end )
     elseif(iop.eq.2)then
        op( index_start:index_end ) = kinetic_wave ( index_start:index_end )
     elseif(iop.eq.3)then
        op( index_start:index_end ) = radial_wave( index_start:index_end )
     elseif(iop.eq.4)then
        op( index_start:index_end ) = local_potential_parts * kinetic_wave( index_start:index_end )
     endif
     index_start = index_end + 1

  enddo

  index_start = 1

  do i_compute_atom = 1, n_compute_atoms, 1

     index_end = rad_index(i_compute_atom)

     op ( index_start:index_end ) = op ( index_start:index_end ) * one_over_dist_tab(i_compute_atom)

     index_start = index_end+1

  enddo

  ! Now tabulate full wave function kinetic energy for each radial function

  ! first, the nonzero functions
  do i_compute_fn = 1, n_compute_fns, 1
    !call mul_vec_2 ( H_times_psi(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), T_plus_V(i_compute_fn) )
     call mul_vec_2 ( op_times_psi(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), op(i_compute_fn) )
  enddo

  ! then, the zero functions
  do i_compute_point = 1, n_zero_compute, 1
     i_compute = zero_index_point(i_compute_point)
    !H_times_psi(i_compute) = 0.0d0
     op_times_psi(i_compute) = 0.0d0
  enddo

 end subroutine evaluate_VTW_psi_rel


 subroutine update_full_matrix_rel ( n_compute, i_basis, ld_dirac, matrix_shell, matrix )

!  PURPOSE
!  Subroutine update_full_matrix adds a part of the integrals in a
!  matrix (overlap or Hamiltonian matrix)
!  (only for the n_compute basis functions that are nonzero at the
!  current integration shell) to the full matrix (which contains
!  all n_basis basis functions). The link between i_compute = 1 ... n_compute
!  and i_basis = 1 ... n_basis is provided by the index array 
!  i_basis(i_compute).
!
!  USES

  use dimensions
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop
  use pbc_lists
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer,intent(in) :: n_compute ! n_compute_c == n_compute_a -- number of relevant basis functions
  integer,intent(in) :: i_basis(n_compute) ! i_basis -- list of relevant basis functions
  integer,intent(in) :: ld_dirac ! dimension
  real*8,intent(in)  :: matrix_shell(n_compute,n_compute) ! hamiltonian / overlap_matrix of relevant basis functions
  real*8,intent(out) :: matrix(ld_dirac) ! date from matrix_shell is added here

  integer :: i_compute_1, i_compute_2, i_offset
  integer :: i_index_real
  integer :: i_offset_first_part
  integer :: i_cell_index, i_cell_1
  integer :: i_max_basis, i_min_basis
  integer :: i_one_part, i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)
  integer :: help
  integer:: i_cell_old

!NEC_CB
  integer::help1(n_compute), help2(n_compute)

  ! For fully-rel cases, we currently use PM_none for benchmark calculations. This will be further expanded to other PM_index.
  ! packed_matrix_format was also set to 0, in subroutine initialize_bc_dependent_lists in pbc_lists.f90 
  packed_matrix_format = PM_none
  select case(packed_matrix_format)

  case(PM_none)

     i_index_real = 0
     do i_compute_2 = 1, n_compute
        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2
        do i_compute_1 = 1,i_compute_2
           i_index_real = i_offset + i_basis(i_compute_1)
           matrix(i_index_real) = matrix(i_index_real) + matrix_shell(i_compute_1, i_compute_2)
        enddo
     enddo


     !!!!!!---------(Rundong) The following code does not work for now !---------------
  case(PM_index) !------------------------------------------------------------------------------------


     if(n_periodic == 0)then

        do i_compute_1 = 1, n_compute, 1
           
           i_start =  index_hamiltonian(1,1, i_basis(i_compute_1))
           i_end   =  index_hamiltonian(2,1, i_basis(i_compute_1))

           do i_compute_2 = 1,i_compute_1,1
              
              i_basis_2 = i_basis(i_compute_2)
              
              place: do i_place = i_start, i_end, 1
                 
                 if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
                    if (i_compute_2.le.i_compute_1) then
     !!!!!             matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                    else
     !!!!!             matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_1, i_compute_2)
                    end if

                    i_index_real = i_place
                    exit place 
                 
                 else if(column_index_hamiltonian( i_place) > i_basis_2)then
                    i_index_real = i_place
                    exit place  

                 end if
              end do place
              i_start = i_index_real
           end do
        end do
        
        
     else ! Periodic case----------------

!NEC_CB Do not read multi-indirectly addressed indices multiple times
        do i_compute_1 = 1, n_compute, 1
          help1(i_compute_1)=Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1)=center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do

        do i_compute_1 = 1, n_compute, 1

           i_basis_1 = help1(i_compute_1)!Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = help2(i_compute_1)!center_to_cell(Cbasis_to_center(i_basis(i_compute_1))

           offset_end = -1
           offset = -1

           do i_cell_1 = 1, n_cells 

             i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

              offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
              offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

           end do


           do i_compute_2 = 1, n_compute
              
              i_basis_2 = help1(i_compute_2)!Cbasis_to_basis(i_basis(i_compute_2))

              if(i_basis_2 <= i_basis_1)then

                 i_cell    =  help2(i_compute_2)!center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))


                 place_2: do i_place = offset(i_cell), offset_end(i_cell),1 
                 
                    if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
                       if (i_compute_2.le.i_compute_1) then
     !!!!!               matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_2, i_compute_1)
                       else
     !!!!                matrix(i_place) = matrix(i_place) + matrix_shell(i_compute_1, i_compute_2)
                       end if

                       exit place_2
                 
                    else if(column_index_hamiltonian( i_place) > i_basis_2)then
                       exit place_2  


                    end if
                 end do place_2

                 offset(i_cell) = i_place

              end if
           end do
        end do

     end if ! n_periodic == 0 - else


  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select

 end subroutine update_full_matrix_rel
!---------------------------------------------------------------------
!



