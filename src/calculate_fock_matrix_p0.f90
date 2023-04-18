!****h* FHI-aims/calculate_fock_matrix
!  NAME
!    calculate_fock_matrix
!  SYNOPSIS

module calculate_fock_matrix_p0

  !  PURPOSE
  !    Calculate Fock matrix
  !  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate, &
      update_when_allocating
  use dimensions
  use localorb_io, only: use_unit, OL_LOW, OL_NORM, OL_HIGH, localorb_info, &
      localorb_allinfo
  use mpi_tasks
  use pbc_lists
  use runtime_choices

  implicit none

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
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE


  !-----------------------------------------------------------------------------
  ! By default all data is private:

  PRIVATE

  ! Only the following routines are public:

  PUBLIC :: init_fock_matrix_calculations, &
            cleanup_fock_matrix_calculations, &
            evaluate_exchange_matr_realspace_p0

  !-----------------------------------------------------------------------------
  !
  ! Please configure the following parameters,
  ! some of them should eventually go to runtime_choices:

  ! Approx error in exchange matrix entries:

  real*8, public :: crit_val = 1.d-7

  ! Rows/columns in the coulomb matrices which are completely below coul_mat_threshold are left away!
  real*8, public :: coul_mat_threshold = 1.d-10

  !----------------------------------------------------------------------------

  ! Flag if derivatives are initialized
  logical :: init_calc_deriv = .false.

  !----------------------------------------------------------------------------

  ! Variables describing work distribution:

  integer :: my_n_atoms, my_n_basis
  integer, allocatable :: my_atom_list(:), my_basis_off(:)

  integer :: my_cm_bb_s, my_cm_bb_e
  integer :: my_cm_cell_start, my_cm_cell_inc, my_cm_cell_num

  integer :: mpi_comm_atom
  integer :: my_atom_id, my_tasks_per_atom

  integer, allocatable :: my_atom_pair_s(:), my_atom_pair_e(:)

  integer :: n_atoms2
  integer, allocatable :: vb_atom2atom(:)

  integer,allocatable :: atom2basis_len2(:), atom2basis_off2(:),atom2vb_basis_off(:),basis_atom2(:)

  !----------------------------------------------------------------------------

  ! ovlp3fn handling:
  ! Since the storage requirements for ovlp3fn may get rather big,
  ! we don't want to waste memory by allocating a big array with
  ! leading dimension max_n_basbas_sp.
  ! Therefore we use an array of arrays which are dimensioned just
  ! to fit the needed size.

  type matrix_2d
    real*8, allocatable :: m(:,:) ! ovlp3fn entries
    real*8, allocatable :: n(:)   ! norm of columns of m
  end type

  type(matrix_2d), allocatable :: ovlp3fn(:), d_ovlp3fn(:,:), AS_d_ovlp3fn(:,:)

  real*8, allocatable :: max_norm_ovlp3fn(:,:,:), max_norm_ovlp3fn_rcv(:,:,:)
  real*8, allocatable :: max_norm_d_ovlp3fn(:,:,:), max_norm_d_ovlp3fn_rcv(:,:,:)
  real*8, allocatable :: AS_max_norm_d_ovlp3fn(:,:,:), AS_max_norm_d_ovlp3fn_rcv(:,:,:)
  real*8, allocatable :: max_norm_ovlp3fn_per_latom(:)
  real*8, allocatable :: max_norm_d_ovlp3fn_per_latom(:)
  real*8, allocatable :: AS_max_norm_d_ovlp3fn_per_latom(:)

  real*8, allocatable :: max_abs_dm_cols_old(:,:,:)
  real*8, allocatable :: max_norm_tmp_old(:)
  !used for LRC-wPBEh since we need two ovlp3fn arrays to construct the fock_matr later
  type(matrix_2d), allocatable :: ovlp3fn_SR(:), d_ovlp3fn_SR(:,:), AS_d_ovlp3fn_SR(:,:)

  real*8, allocatable :: max_norm_ovlp3fn_SR(:,:,:), max_norm_ovlp3fn_rcv_SR(:,:,:)
  real*8, allocatable :: max_norm_d_ovlp3fn_SR(:,:,:), max_norm_d_ovlp3fn_rcv_SR(:,:,:)
  real*8, allocatable :: AS_max_norm_d_ovlp3fn_SR(:,:,:), AS_max_norm_d_ovlp3fn_rcv_SR(:,:,:)
  real*8, allocatable :: max_norm_ovlp3fn_per_latom_SR(:)
  real*8, allocatable :: max_norm_d_ovlp3fn_per_latom_SR(:)
  real*8, allocatable :: AS_max_norm_d_ovlp3fn_per_latom_SR(:)

  !----------------------------------------------------------------------------

  ! Lower and upper boundary of the density matrix cells

  integer :: lbnd_bvk_cell(3), ubnd_bvk_cell(3)

  ! Mapping of 3D index in density matrix to cell number in density matrix

  integer, allocatable :: bvk_cell_idx(:,:,:)
  integer, allocatable :: inv_cell_bvk(:)

  integer, allocatable :: add_cells(:,:)
  integer, allocatable :: sub_cells(:,:)

  !----------------------------------------------------------------------------

  integer :: n_atom_pairs
  integer :: my_n_atom_pairs

  integer, allocatable :: pair_offset(:,:,:)

  logical, allocatable :: pair_flag_bvk(:,:,:)
  logical, allocatable :: my_pair_flag_bvk(:,:,:)

  integer, allocatable :: pair_list(:,:)

  !-----------------------------------------------------------------------------

  ! Coulomb matrix handling:
  ! Every coulomb matrix is stored compressed, i.e. rows/columns which are
  ! completely below a threshold are left away.
  ! The remaining row/cols are stored in matrix form for being able to use
  ! library routines like DGEMM.
  ! The matrix may still contain zero entries, i.e. this is not a sparse
  ! storage format

  type coul_mat_t
    integer n_rows, n_cols ! Number of rows/columns stored
    integer, allocatable :: row_idx(:) ! Index of stored rows within total rows
    integer, allocatable :: col_idx(:) ! Index of stored cols within total cols
    real*8, allocatable :: mat(:,:)    ! Stored matrix elements
  end type

  type(coul_mat_t), allocatable :: coul_mat_store(:,:,:), d_coul_mat_store(:,:,:,:), AS_d_coul_mat_store(:,:,:,:)

  integer*8 :: n_coulmat_elems(3)

  real*8, allocatable :: coul_mat_norm(:,:,:), d_coul_mat_norm(:,:,:), AS_d_coul_mat_norm(:,:,:)

  !used for LRC-wPBEh since we need two ovlp3fn arrays to construct the fock_matr later
  type(coul_mat_t), allocatable :: coul_mat_store_SR(:,:,:), d_coul_mat_store_SR(:,:,:,:), AS_d_coul_mat_store_SR(:,:,:,:)

  real*8, allocatable :: coul_mat_norm_SR(:,:,:), d_coul_mat_norm_SR(:,:,:), AS_d_coul_mat_norm_SR(:,:,:)

  real*8 :: time_mult_coul_mat

  ! For output using localorb_io
  character*8192 :: info_str

  ! Only for debugging, should be removed afterwards
  real*8 :: Rvec_length, Min_Rvec_length

  ! FK: Flag if parts for analytical stress are initialized
  logical :: AS_init_stress = .false.
  ! FK: Index variables for analytical stress
  integer                 :: AS_index
  ! Map index -> (l_index, m_index)
  ! 1 2 3
  ! 7 4 5
  ! 8 9 6
  integer, dimension(1:9) :: AS_l_index = (/ 1,1,1,2,2,3,2,3,3 /)
  integer, dimension(1:9) :: AS_m_index = (/ 1,2,3,2,3,3,1,1,2 /)

  ! For the LC-wPBEh functional with a hybrid_xc_coeff other than zero, we need to
  ! calculate the ovlp3fn_SR as well. So we need to run the calculation twice.
  ! To make sure the right ovlp3fn is
  ! calculated, there needs to be a control variable, which is following:
  logical :: lc_wpbeh_lr_run = .false.

contains
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/init_fock_calculations
  !  NAME
  !    init_fock_calculations
  !  SYNOPSIS

  subroutine init_fock_matrix_calculations(opt_calc_deriv, AS_opt_stress)

    !  PURPOSE
    !
    !  Does all initial settings for calculate_fock_matrix_p0 which have to
    !  be done once per SCF cycle.
    !  Most of the work done here should go somewhere else (e.g. condense_basis_pairs_v2)
    !  I left it here for convenience.
    !
    !  USES

    use basis, only: max_n_basis_sp, max_n_basis_sp2 ,sp2n_basis_sp, atom2basis_off, atom_radius, outer_radius, basis_fn
    use lvl_triples, only: initialize_lvl_triples, cleanup_lvl_triples
    use geometry, only: species, coords, lattice_vector
    use prodbas, only: OVLP_TYPE_COULOMB, OVLP_TYPE_HSE, OVLP_TYPE_LR, &
        OVLP_TYPE_CUT, ovlp_type_bare_or_hse_coul, max_n_basbas_sp, &
        atom2basbas_off, sp2n_basbas_sp
    use species_data, only: no_basis, species_pseudoized
    use sym_base, only: calculate_realspace_coulomb_sym
    use synchronize_mpi_basic, only: sync_logical, sync_vector
    use tight_binding_auxmat, only: initialize_tb_auxmat, deallocate_tb_auxmat,&
        fast_calculate_tb_auxmat

    use species_data, only : r_cutoff, outer_partition_radius

    implicit none

    !  ARGUMENTS

    logical, intent(in) :: opt_calc_deriv
    logical, intent(in) :: AS_opt_stress

    !  INPUTS
    !    o opt_calc_deriv -- flag if gradients of ovlp and Coulomb matrices should be calculated
    !  OUTPUTS
    !    None
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE


    ! Maximum number of coeff_3fn arrays calculated in one call
    ! integer, parameter :: n_max_coeff_3fn = 16 (Now defined in dimension.f90, XR, 2017.1.7)

    integer i, j, k, d, jj, i1, i2, i3, n, i_cell, i_atom, i_atom_l, i_atom_r, i_my_atom, i_split
    integer mpierr, nbb1, nbb2
    real*8, allocatable :: aux_mat(:,:,:), d_aux_mat(:,:,:,:), AS_d_aux_mat(:,:,:,:)
    real*8 work, total_work

    integer ii, j_atom, n_bb_s, n_bb_e, n_cnt
    real*8 Rvec(3), rad_sum, cmem, omem
    real*8, allocatable :: Rvecs(:,:)
    real*8, allocatable :: coeff_3fn(:,:,:,:,:), d_coeff_3fn(:,:,:,:,:,:)
    real*8, allocatable :: aux(:,:)
    real*8, allocatable :: d_coul_mat_norm_3(:,:,:,:)
    real*8, allocatable :: AS_d_coul_mat_norm_components(:,:,:,:)
    real*8, allocatable :: weight(:)
    real*8, allocatable :: vb_atom_radius(:)
    logical, allocatable :: pair_flag_full(:,:,:,:,:)

    integer n_pairs, n_super, n_pairs_cur, i_start, i_atom_pair
    integer, allocatable :: cell_list(:), tasks_per_atom(:), task_of_atom(:)
    integer i_task, my_atom

    real*8 :: ttt0

    ! Logic to make sure we only create atom communicators if that is really needed!
    logical, save :: distribution_changed = .true.
    logical, save :: atom_comm_created = .false.
    integer, save :: my_previous_atom = -1

    integer :: n_real_atoms_DB, n_atoms_save

    integer :: my_cm_cell_num_old, my_cm_cell_start_old, my_cm_cell_inc_old
    logical :: my_iterate
    logical :: cm_calculated(n_atoms,n_atoms)

    ! indices of largest coul_mat_store elements (used to avoid mass allocations in LVL_fast)
    integer :: peakcol(3), peakrow(3)
    integer :: index_1, index_2, index_3, maxrow, maxcol, colcount

    integer :: info
    character(*), parameter :: func = 'init_fock_matrix_calculations'

    integer :: n_very_big, i_basis, i_nbb
    real*8  :: atom2basis_len_sq(n_atoms),atom2basis_len_sq_av

    integer :: split(n_atoms), remain, divide, basis_off, min_split
    real*8  :: load(n_atoms), av_load_per_task, var_load, tmpsplit(n_atoms)

    init_calc_deriv = opt_calc_deriv
    AS_init_stress  = AS_opt_stress


    write(info_str,'(A)') '  -----------------------------------------'
    call localorb_info ( info_str,use_unit )
    write(info_str,'(A)') '  --- Initializing Fock matrix calculations'
    call localorb_info ( info_str,use_unit )
    write(info_str,'(A)') '  -----------------------------------------'
    call localorb_info ( info_str,use_unit )

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Presetting of parameters
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    ! Output crit_val and coul_mat_threshold

    write(info_str,'(A,F15.12)') 'screening_threshold (crit_val) = ', crit_val
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write(info_str,'(A,F15.12)') 'coul_mat_threshold = ', coul_mat_threshold
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Initialize indices needed in calculations
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

! DB 10/02/13
    ! For QM/MM embedding one mights to put additional integration grids around
    ! embedding monopoles (emty sites without any basis functions).
    ! This routine has a problem when dealing with such empty atoms, e.g. dgemm
    ! is annoyed.
    ! The workaround here is that we change the number of atoms to those which are
    ! real atoms with basis functions (and change it back at the end of the routine).
    ! This is not the nicest way to do it, but at least we maintain this highly
    ! optimzed structure of this routine.
    ! The only drawback: this relies on correct ordering of atoms
    ! (1:n_real_atoms_DB,1:n_empty_atoms,1:n_pp_atoms)

    n_real_atoms_DB = 0
    do i_atom = 1, n_atoms
       ! AL/HJ 28/09/2018
       ! This previous definition was inconsistent with the if statement 7 lines
       ! below. Made need further assessment if anyone ever does pseudocores
       ! with basis functions on top, but for now this is bullet proof.
       !if (no_basis(species(i_atom))) cycle
       if (species_pseudoized(species(i_atom)).or.no_basis(species(i_atom))) cycle
       n_real_atoms_DB = n_real_atoms_DB +1
    enddo
    n_atoms_save = n_atoms
    n_atoms = n_real_atoms_DB
    ! check for correct ordering
    do i_atom = 1, n_atoms
      ! AL/HJ 28/09/2018
      ! Debug statement for figuring out error with pseudocore counting
      ! write(use_unit,*) 'DEBUG: atom_number ', i_atom
      if(species_pseudoized(species(i_atom)).or.no_basis(species(i_atom))) then
         call aims_stop('Internal: inconsistent ordering of atoms')
      endif
    enddo

    !--------------------------------------------------------------------------
    ! SK: Check for very big atoms, that is, atoms with many basis functions,
    !     and split them into fractions.
    !--------------------------------------------------------------------------

    if (flag_split_atoms) then

      write(info_str, *) 'Trying to make the workload more even by internally splitting atoms into fractions'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      if (.not. flag_split_batch) split_batch = 14.d0
      do i_atom = 1, n_atoms
       load(i_atom) = dble(atom2basis_len(i_atom))
       split(i_atom)= max(1,int(load(i_atom)/split_batch)) ! Splitting basis_fn/atom into batches of split_batch
       write(info_str, *) 'Atom: ',i_atom, ' Basis/atom: ', load(i_atom)
       call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      enddo

      min_split = minval(split(:))
      do i_atom = 1, n_atoms
        if (flag_split_min_val) then
          split(i_atom) = (split(i_atom)/ min_split) * split_min_val
        endif
        if (flag_split_max_val) then
          split(i_atom) = min(split_max_val,split(i_atom))
        endif
      enddo

    else

     split(:) = 1

    endif

    n_atoms2 = sum(split(:))
    if (.not. allocated(basis_atom2)) allocate(basis_atom2(n_basis))
    if (.not. allocated(vb_atom2atom)) allocate(vb_atom2atom(n_atoms2))
    if (.not. allocated(atom2basis_len2)) allocate(atom2basis_len2(n_atoms2))
    if (.not. allocated(atom2basis_off2)) allocate(atom2basis_off2(n_atoms2))
    if (.not. allocated(atom2vb_basis_off)) allocate(atom2vb_basis_off(n_atoms2))
    atom2vb_basis_off = 0
    atom2basis_len2 = 0

    n_atoms2 = 0

    do i_atom = 1, n_atoms
      remain = mod(atom2basis_len(i_atom),split(i_atom))
      divide = atom2basis_len(i_atom)/split(i_atom)
      basis_off = 0
      do i_split = 1, split(i_atom)
       n_atoms2 = n_atoms2 + 1
       atom2basis_len2(n_atoms2) = divide + min(i_split,remain) - min(i_split-1,remain)
       atom2basis_off2(n_atoms2) = atom2basis_off(i_atom) + basis_off
       atom2vb_basis_off(n_atoms2) = basis_off
       vb_atom2atom(n_atoms2) = i_atom
       basis_off = basis_off + atom2basis_len2(n_atoms2)
      enddo
    enddo

    if (n_atoms == n_atoms2) then
      write(info_str, *) 'No atoms are split internally'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    else
      write(info_str, *) 'There are now: ', n_atoms2, ' atoms instead of: ', n_atoms
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      write(info_str,'(3A10)') 'new', 'old', 'bas/atom'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      do i_atom = 1, n_atoms2
       write(info_str, '(3i10)') i_atom,vb_atom2atom(i_atom),atom2basis_len2(i_atom) !,atom2basis_off2(i_atom)
       call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      enddo

    endif

    do i_atom = 1, n_atoms2
     basis_atom2(atom2basis_off2(i_atom)+1:atom2basis_off2(i_atom)+atom2basis_len2(i_atom))=i_atom
    enddo

    max_n_basis_sp2 = maxval(atom2basis_len2)

    ! This code assumes that all basis functions are ordered according to atom ordering
    ! Just for safety, check this assumption here for the case that somebody
    ! should change the code in an incompatible way:

    n = 0
    do i_atom = 1, n_atoms2
      if(atom2basis_off2(i_atom) /= n) &
        call aims_stop('init_fock_matrix_calculations: atom/basis ordering is not like assumed')
      n = n + atom2basis_len2(i_atom)
    enddo



    ! BvK cells boundaries

    lbnd_bvk_cell(:) = -(n_k_points_xyz_nosym(:)-1)/2
    ubnd_bvk_cell(:) =   n_k_points_xyz_nosym(:)/2

    ! Mapping from 3D cell to BvK cell number

    allocate(bvk_cell_idx(lbnd_bvk_cell(1):ubnd_bvk_cell(1),lbnd_bvk_cell(2):ubnd_bvk_cell(2),lbnd_bvk_cell(3):ubnd_bvk_cell(3)), stat=info)
    call check_allocation(info, 'bvk_cell_idx', func)
    bvk_cell_idx(:,:,:) = 0
    do i = 1, n_cells_bvk
      bvk_cell_idx(cell_index_bvk(i,1),cell_index_bvk(i,2),cell_index_bvk(i,3)) = i
    enddo

    allocate(inv_cell_bvk(n_cells_bvk), stat=info)
    call check_allocation(info, 'inv_cell_bvk', func, n_cells_bvk)
    do i = 1, n_cells_bvk
      inv_cell_bvk(i) = get_bvk_cell_idx(-cell_index_bvk(i,:))
    enddo

    allocate(add_cells(n_cells_bvk,n_cells_bvk), stat=info)
    call check_allocation(info, 'add_cells', func, n_cells_bvk, n_cells_bvk)
    allocate(sub_cells(n_cells_bvk,n_cells_bvk), stat=info)
    call check_allocation(info, 'sub_cells', func, n_cells_bvk, n_cells_bvk)
    do i = 1, n_cells_bvk
    do j = 1, n_cells_bvk
      add_cells(i,j) = get_bvk_cell_idx(cell_index_bvk(i,:) + cell_index_bvk(j,:))
      sub_cells(i,j) = get_bvk_cell_idx(cell_index_bvk(i,:) - cell_index_bvk(j,:))
    enddo
    enddo

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Get the total atom pairs within all super cells
    ! and the atom pairs within BvK cells
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

     write(info_str, '(a,3i8)') '  Number of super-cells in X Y Z:',number_of_super_cells(:)
     call localorb_info ( info_str,use_unit,'(A)', OL_norm )


    ! Flag all atom pairs which contribute integrals to ovlp3fn
    ! Please note:
    ! We are using the criterion that the distance between atoms must be
    ! smaller than the sum of atom_radius in order to contribute.
    ! If this assumption has to be changed, it must be changed here (and only here!)
    !
    ! In this code, there is made usage of the fact that every pair with
    ! distance Rvec has a corresponding pair with distance -Rvec and only
    ! half of the integrals has to be stored.
    ! In order to avoid the slightest chance that a symmetric pair is missed
    ! due to numeric inaccuracy, we go only over half of the supercells
    ! and flag both pairs at once.

    if (.not. allocated(vb_atom_radius)) allocate(vb_atom_radius(n_atoms2))
    vb_atom_radius = 0.0
    do i_atom = 1, n_atoms2
     do i_basis = 1+atom2basis_off2(i_atom),atom2basis_len2(i_atom)+atom2basis_off2(i_atom)
      vb_atom_radius(i_atom) = max(vb_atom_radius(i_atom),outer_radius(basis_fn(i_basis)))
     enddo
    enddo

    allocate(pair_flag_full(n_atoms2, n_atoms2, &
                            -number_of_super_cells(1):number_of_super_cells(1), &
                            -number_of_super_cells(2):number_of_super_cells(2), &
                            -number_of_super_cells(3):number_of_super_cells(3)), &
                            stat=info)
    call check_allocation(info, 'pair_flag_full', func)
    pair_flag_full(:,:,:,:,:) = .false.

    do i1 = 0, number_of_super_cells(1)
      do i2 = -number_of_super_cells(2), number_of_super_cells(2)
        if(i1==0 .and. i2<0) cycle
        do i3 = -number_of_super_cells(3), number_of_super_cells(3)
          if(i1==0 .and. i2==0 .and. i3<0) cycle
          do j=1,n_atoms2
            do i=1,n_atoms2
              if(i1==0 .and. i2==0 .and. i3==0 .and. i>j) cycle !SK checked: is correct
              Rvec = coords(:,vb_atom2atom(j)) - coords(:,vb_atom2atom(i)) + matmul(lattice_vector, (/ i1, i2, i3 /))
              rad_sum = vb_atom_radius(i) + vb_atom_radius(j) !skchange
!              rad_sum = atom_radius(species(vb_atom2atom(i))) + atom_radius(species(vb_atom2atom(j))) !skchange
!              rad_sum = r_cutoff(species(vb_atom2atom(i))) + r_cutoff(species(vb_atom2atom(j)))
              if(sum(Rvec(:)**2) <= rad_sum**2) then
                pair_flag_full(i,j, i1, i2, i3) = .true.
                pair_flag_full(j,i,-i1,-i2,-i3) = .true.
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

    write(info_str, *)' Number of atoms pairs total: ',count(pair_flag_full)
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str, *) ' Atom pair distribution:'
    call localorb_info ( info_str,use_unit,'(A)', OL_low )

    do i1 = -number_of_super_cells(1), number_of_super_cells(1)
      do i2 = -number_of_super_cells(2), number_of_super_cells(2)
        do i3 = -number_of_super_cells(3), number_of_super_cells(3)
          if(count(pair_flag_full(:,:,i1,i2,i3)) == 0) cycle
            write(info_str, '(a,3i8,a,i8)')  '  Cell: ',i1,i2,i3,' atom pairs: ',count(pair_flag_full(:,:,i1,i2,i3))
            call localorb_info ( info_str,use_unit,'(A)', OL_low )
        enddo
      enddo
    enddo

    ! Now get the atom pairs within the BvK cells
    ! If the number of BvK cells is big enough (>= 2*number_of_super_cells+1 in every direction),
    ! these are the same as we got above.
    ! Otherwise, some pairs coincide into the same BvK cell

    allocate(pair_flag_bvk(n_atoms2, n_atoms2, n_cells_bvk), stat=info)
    call check_allocation(info, 'pair_flag_bvk', func, n_atoms2, n_atoms2, n_cells_bvk)
    pair_flag_bvk(:,:,:) = .false.

    do i1 = -number_of_super_cells(1), number_of_super_cells(1)
      do i2 = -number_of_super_cells(2), number_of_super_cells(2)
        do i3 = -number_of_super_cells(3), number_of_super_cells(3)
          do j=1,n_atoms2
            do i=1,n_atoms2
              if(pair_flag_full(i,j,i1,i2,i3)) then

                i_cell = get_bvk_cell_idx( (/ i1, i2, i3 /) )
                pair_flag_bvk(i,j,i_cell) = .true.
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

    n_atom_pairs = count(pair_flag_bvk)
    write(info_str, *)' Number of atoms pairs in BvK cells: ',n_atom_pairs
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Work distribution
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    write(info_str, *)' Initializing work distribution for fock matrix calculation'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

    allocate(my_pair_flag_bvk(n_atoms2, n_atoms2, n_cells_bvk), stat=info)
    call check_allocation(info, 'my_pair_flag_bvk', func, n_atoms2, n_atoms2, n_cells_bvk)
    my_pair_flag_bvk(:,:,:) = .false.

    allocate(weight(n_atoms2), stat=info)
    call check_allocation(info, 'weight', func, n_atoms2)
    allocate(task_of_atom(n_atoms2), tasks_per_atom(n_atoms2), stat=info)
    call check_allocation(info, 'task_of_atom', func, n_atoms2)
    call check_allocation(info, 'tasks_per_atom', func, n_atoms2)

    ! We have to set the weight per atom for the subdivision which should reflect
    ! the time (and memory) needed per atom.
    ! This is hard to guess, currently we use the number of basbas functions:

    do i_atom = 1, n_atoms2
      weight(i_atom) = atom2basbas_len(i_atom) !?SK or:sqrt(float(atom2basbas_len(i_atom)*atom2basis_len2(i_atom)))
    enddo

    call divide_atoms(weight, task_of_atom, tasks_per_atom)

    allocate(my_atom_list(n_atoms2), stat=info)
    call check_allocation(info, 'my_atom_list', func, n_atoms2)

    if(count(task_of_atom==myid)>1) then

      ! Current task has more than 1 atom, which are never shared

      my_n_atoms = 0
      do i_atom = 1, n_atoms2
        if(task_of_atom(i_atom) == myid) then
           my_n_atoms = my_n_atoms+1
           my_atom_list(my_n_atoms) = i_atom
        endif
      enddo
      if(my_n_atoms <= 1) call aims_stop('Internal: my_n_atoms <= 1')

      my_atom = -999999999 ! Not defined in this case
      my_atom_id = 0
      my_tasks_per_atom = 1

    else

      ! Current task has 1 atom (which is potentially shared)

      ! Search my atom - task_of_atom contains the start task for every atom

      my_atom = -1
      do i_atom = 1, n_atoms2
        if(myid>=task_of_atom(i_atom) .and. myid < task_of_atom(i_atom)+tasks_per_atom(i_atom)) then
          my_atom = i_atom
          exit
        endif
      enddo

      if(my_atom < 0) call aims_stop('Internal: my_atom not found')

      my_n_atoms = 1
      my_atom_list(1) = my_atom

      my_atom_id = myid - task_of_atom(my_atom)
      my_tasks_per_atom = tasks_per_atom(my_atom)

    endif

    deallocate(weight)
    deallocate(task_of_atom, tasks_per_atom)

    ! Check if the distribution changed. We need to check so that we
    ! only ever destroy and recreate the atom communicators if really needed.
    ! Some MPI libraries will collapse if too many communicators are created
    ! during a single run.
    !
    ! my_previous_atom is set to -1 at the outset, so my_atom should not
    ! ever be equal to -1 at this point.

    distribution_changed = (my_previous_atom.ne.my_atom)

    if (distribution_changed) then
       ! For now, this output is kept at "high" priority to ensure that we find
       ! any possible problems. Should be downgraded to OL_norm later.
       write(info_str,'(2X,a,i5,a)') 'Task ',myid,' : Coulomb matrix distribution changed.'
       call localorb_allinfo(info_str,use_unit,'(A)',OL_high)
    end if

    call sync_logical (distribution_changed,SYNC_OR)

    my_previous_atom = my_atom

    if (distribution_changed) then
       ! For now, this output is kept at "high" priority to ensure that we find
       ! any possible problems. Should be downgraded to OL_norm later.
       write(info_str,'(2X,a)') 'Overall Coulomb matrix distribution changed.'
       call localorb_info(info_str,use_unit,'(A)',OL_high)
    end if

    write(info_str,'(2X,a)') 'Coulomb matrix subdivision:'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    if(my_tasks_per_atom > 1) then

      ! Distribute atom pairs within tasks sharing the atom.

      i_atom_l = my_atom
      my_n_atom_pairs = count(pair_flag_bvk(i_atom_l,:,:)) ! atom_pairs for all tasks sharing atom, not really mine

      ! Sort pairs by distance and distribute the sorted list round robin on CPUs.
      ! I hope this will give a better load balancing since pairs with less distance
      ! are expensive and these are distributed equally on all CPUs.

      allocate(aux(3,my_n_atom_pairs), stat=info)
      call check_allocation(info, 'aux', func, 3, my_n_atom_pairs)
      n = 0
      do i_atom_r=1,n_atoms2
        do i_cell=1,n_cells_bvk
          if(pair_flag_bvk(i_atom_l,i_atom_r,i_cell)) then
            n = n+1
            aux(1,n) = i_atom_r
            aux(2,n) = i_cell
            Rvec = coords(:,vb_atom2atom(i_atom_r)) - coords(:,vb_atom2atom(i_atom_l)) + matmul(lattice_vector, cell_index_bvk(i_cell,:))
            aux(3,n) = sum(Rvec(:)**2)
          endif
        enddo
      enddo

      call heapsort_general(aux, 3, my_n_atom_pairs, 3)

      do n = 1, (my_n_atom_pairs-1)/my_tasks_per_atom+1
        if(mod(n,2) == 1) then
          i = (n-1)*my_tasks_per_atom + my_atom_id + 1
        else
          i = n*my_tasks_per_atom - my_atom_id
        endif
        if(i<=my_n_atom_pairs) then
          i_atom_r = aux(1,i)
          i_cell   = aux(2,i)
          my_pair_flag_bvk(i_atom_l,i_atom_r,i_cell) = .true.
        endif
      enddo

      deallocate(aux)

      ! Coulomb matrix distribution:
      if(n_cells_bvk<my_tasks_per_atom) then
        ! Divide the rows of the Coulomb matrix
        ! my_cm_bb_s/e are start and end row for my (one and only) atom
        my_cm_bb_s = (my_atom_id*atom2basbas_len(my_atom))/my_tasks_per_atom + 1
        my_cm_bb_e = ((my_atom_id+1)*atom2basbas_len(my_atom))/my_tasks_per_atom
        my_cm_cell_start = 1
        my_cm_cell_inc = 1
        my_cm_cell_num = n_cells_bvk
        write(info_str,'(2(a,i5),2(a,i10))') '  Task ',myid,' atom ',my_atom,' has Coulomb matrix rows ', &
                                             my_cm_bb_s,' to ',my_cm_bb_e
      else
        ! Divide the cells of the Coulomb matrix:
        ! It is important to do this round robin like so that no task has a chunk
        ! of contigous cells which may lead to load balancing problems.
        my_cm_bb_s = 1
        my_cm_bb_e = atom2basbas_len(my_atom)
        my_cm_cell_start = 1 + my_atom_id
        my_cm_cell_inc = my_tasks_per_atom
        my_cm_cell_num = (n_cells_bvk - my_cm_cell_start)/my_tasks_per_atom + 1
        write(info_str,'(2(a,i5),2(a,i10))') '  Task ',myid,' atom ',my_atom,' has ',my_cm_cell_num,&
                                             ' Coulomb matrix cells'
      endif

      ! Create communicator for tasks sharing one atom

      if ( distribution_changed ) then
         if (atom_comm_created) then
            ! Place a safeguard here to avoid ever freeing up the global communicator.
            ! In principle, this should never happen.
            if (mpi_comm_atom.ne.mpi_comm_global) then
               call mpi_comm_free(mpi_comm_atom,mpierr)
            end if
         end if
         call mpi_comm_split(mpi_comm_global, my_atom, myid, mpi_comm_atom, mpierr)
         atom_comm_created = .true.
      end if

    else

      ! Atom pair distribution

      do i_my_atom = 1, my_n_atoms
        i_atom_l = my_atom_list(i_my_atom)
        do i_atom_r=1,n_atoms2
          do i_cell=1,n_cells_bvk
            if(pair_flag_bvk(i_atom_l,i_atom_r,i_cell)) my_pair_flag_bvk(i_atom_l,i_atom_r,i_cell) = .true.
          enddo
        enddo
      enddo

      ! Coulomb matrix distribution:
      ! my_cm_bb_s/e make no sense (especially if the task has several atoms)
      ! and must never be used in this case
      my_cm_bb_s =  999999999 ! for safety
      my_cm_bb_e = -999999999 ! for safety
      my_cm_cell_start = 1
      my_cm_cell_inc = 1
      my_cm_cell_num = n_cells_bvk

      ! The MPI communicator must never be used in this case,
      ! however we must participate in mpi_comm_split.
      ! When providing a color of MPI_UNDEFINED we get a communicator MPI_COMM_NULL.

      if ( distribution_changed ) then
         if (atom_comm_created) then
            ! Place a safeguard here to avoid ever freeing up the global communicator.
            ! In principle, this should never happen.
            if (mpi_comm_atom.ne.mpi_comm_global) then
               call mpi_comm_free(mpi_comm_atom,mpierr)
            end if
         end if
         call mpi_comm_split(mpi_comm_global, MPI_UNDEFINED, myid, mpi_comm_atom, mpierr)
         atom_comm_created = .false.
      end if

      write(info_str,'(a,i5,a)') '  Task ',myid,' works on full atoms (no Coulomb matrix division)'

    endif
    call localorb_info ( info_str,use_unit,'(A)', OL_norm  )

    write(info_str,'(a,i6,a,i6,a,20i6)') '  Task ',myid,': my_n_atoms: ',my_n_atoms, &
                                         ' my_atom_list: ',my_atom_list(1:min(my_n_atoms,20))
    if(my_n_atoms>20) info_str = trim(info_str) // ' ... more follow'
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)

    ! get the number of basis functions for this task and the offset for every atom

    my_n_basis = 0
    allocate(my_basis_off(n_atoms2), stat=info)
    call check_allocation(info, 'my_basis_off', func, n_atoms2)
    my_basis_off(:) = -1 ! i.e. atom not owned by task
    do i_my_atom = 1, my_n_atoms
      i_atom = my_atom_list(i_my_atom)
      my_basis_off(i_atom) = my_n_basis
      my_n_basis = my_n_basis + atom2basis_len2(i_atom)
    enddo

    my_n_atom_pairs = count(my_pair_flag_bvk)

    write(info_str,'(2(a,i6))') '  Task ',myid,': my_n_atom_pairs: ',my_n_atom_pairs
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)

    !---------------------------------------------------------------------------
    ! Output the memory for ovlp3fn and Coulomb matrix BEFORE they are allocated
    ! so that the user has an estimate even if the allocate fails

    omem = 0
    do i_my_atom = 1, my_n_atoms
      i_atom_l = my_atom_list(i_my_atom)
      n = 0
      do i_cell=1,n_cells_bvk
        do i_atom_r=1,n_atoms2
          if(my_pair_flag_bvk(i_atom_l,i_atom_r,i_cell)) then
            n = n + atom2basis_len2(i_atom_l)*atom2basis_len2(i_atom_r)
          endif
        enddo
      enddo
      omem = omem + dble(atom2basbas_len(i_atom_l))*dble(n)
    enddo

    if(my_tasks_per_atom > 1) then
      cmem = dble(n_basbas)*dble(max(my_cm_bb_e-my_cm_bb_s+1,0))*dble(my_cm_cell_num)
    else
      cmem = 0
      do i_my_atom = 1, my_n_atoms
        i_atom = my_atom_list(i_my_atom)
        cmem = cmem + dble(n_basbas)*dble(atom2basbas_len(i_atom))*dble(n_cells_bvk)
      enddo
    endif

    if (init_calc_deriv) then
      if (AS_init_stress) then
        omem = (1+3+AS_components)*omem
        cmem = (1+3+AS_components)*cmem
      else
        omem = (1+3)*omem
        cmem = (1+3)*cmem
      end if
    end if

    write(info_str, *)' Memory needs for ovlp3fn and Coulomb matrix:'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str, *) ' For Coulomb matrix, this is the upper limit if no compression occurs'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    if (use_lc_wpbeh .and. hybrid_coeff .ne. 0.0d0) then
	    write(info_str,'(a,2(f12.3,a))') '  ovlp3fn: ',8*2*omem/1.d6,' MB, Coulomb: ',8*2*cmem/1.d6,' MB' !should be the upper limit for each array
	 else
	   write(info_str,'(a,2(f12.3,a))') '  ovlp3fn: ',8*omem/1.d6,' MB, Coulomb: ',8*cmem/1.d6,' MB'
	 end if
	 call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    !---------------------------------------------------------------------------

    ! Set the mapping of basis pairs to storage in ovlp3fn

    allocate(pair_offset(n_atoms2, n_atoms2, n_cells_bvk), stat=info)
    call check_allocation(info, 'pair_offset', func, n_atoms2, n_atoms2, n_cells_bvk)
    pair_offset = -1 ! i.e. pair not present

    ! Allocate the ovlp3fn array itself, the storage for each atom is allocated in the loop

    allocate(ovlp3fn(n_atoms2), stat=info)
    call check_allocation(info, 'ovlp3fn', func, n_atoms2)
    if(init_calc_deriv) then
      allocate(d_ovlp3fn(n_atoms2,3), stat=info)
      call check_allocation(info, 'd_ovlp3fn', func, n_atoms2, 3)
    endif
    if(AS_init_stress) then
      allocate(AS_d_ovlp3fn(n_atoms2,AS_components), stat=info)
      call check_allocation(info, 'AS_d_ovlp3fn', func, n_atoms2, AS_components)
    endif
    if(use_lc_wpbeh .and. hse_omega_hf /= 0.d0 .and. hybrid_coeff /= 0.d0) then
      allocate(ovlp3fn_SR(n_atoms2))
      call check_allocation(info, 'ovlp3fn_SR', func, n_atoms2)
      if(init_calc_deriv) then
        allocate(d_ovlp3fn_SR(n_atoms2,3), stat=info)
        call check_allocation(info, 'd_ovlp3fn_SR', func, n_atoms2, 3)
      endif
      if(AS_init_stress) then
        allocate(AS_d_ovlp3fn_SR(n_atoms2,AS_components), stat=info)
        call check_allocation(info, 'AS_d_ovlp3fn_SR', func, n_atoms2, AS_components)
      endif
    endif

    do i_my_atom = 1, my_n_atoms
      i_atom_l = my_atom_list(i_my_atom)
      n = 0
      do i_cell=1,n_cells_bvk
        do i_atom_r=1,n_atoms2
          if(my_pair_flag_bvk(i_atom_l,i_atom_r,i_cell)) then
            pair_offset(i_atom_l,i_atom_r,i_cell) = n
            n = n + atom2basis_len2(i_atom_l)*atom2basis_len2(i_atom_r)
          endif
        enddo
      enddo
      ! Allocate ovlp3fn for this atom to the needed size
      write(info_str, '(A,3i8)')'Size of ovlp3fn for atom',n,i_atom_l,atom2basbas_len(i_atom_l)
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )


      call aims_allocate(ovlp3fn(i_atom_l)%m,atom2basbas_len(i_atom_l),n,'ovlp3fn for atom')
      call aims_allocate(ovlp3fn(i_atom_l)%n,n,'norm ovlp3fn for atom')
      if(init_calc_deriv) then
        call aims_allocate(d_ovlp3fn(i_atom_l,1)%m,atom2basbas_len(i_atom_l),n,'d_ovlp3fn for atom')
        call aims_allocate(d_ovlp3fn(i_atom_l,2)%m,atom2basbas_len(i_atom_l),n,'d_ovlp3fn for atom')
        call aims_allocate(d_ovlp3fn(i_atom_l,3)%m,atom2basbas_len(i_atom_l),n,'d_ovlp3fn for atom')
        call aims_allocate(d_ovlp3fn(i_atom_l,1)%n,n,'norm d_ovlp3fn for atom')
        call aims_allocate(d_ovlp3fn(i_atom_l,2)%n,n,'norm d_ovlp3fn for atom')
        call aims_allocate(d_ovlp3fn(i_atom_l,3)%n,n,'norm d_ovlp3fn for atom')
      endif
      if(AS_init_stress) then
        do AS_index = 1, AS_components, 1
          call aims_allocate(AS_d_ovlp3fn(i_atom_l,AS_index)%m,atom2basbas_len(i_atom_l),n,'AS_d_ovlp3fn for atom')
          call aims_allocate(AS_d_ovlp3fn(i_atom_l,AS_index)%n,n,'norm AS_d_ovlp3fn for atom')
        enddo
      endif
      if(use_lc_wpbeh .and. hse_omega_hf /= 0.d0 .and. hybrid_coeff /= 0.d0) then
        call aims_allocate(ovlp3fn_SR(i_atom_l)%m,atom2basbas_len(i_atom_l),n,'ovlp3fn_SR for atom')
        call aims_allocate(ovlp3fn_SR(i_atom_l)%n,n,'norm ovlp3fn_SR for atom')
        if(init_calc_deriv) then
		    call aims_allocate(d_ovlp3fn_SR(i_atom_l,1)%m,atom2basbas_len(i_atom_l),n,'d_ovlp3fn_SR for atom')
		    call aims_allocate(d_ovlp3fn_SR(i_atom_l,2)%m,atom2basbas_len(i_atom_l),n,'d_ovlp3fn_SR for atom')
		    call aims_allocate(d_ovlp3fn_SR(i_atom_l,3)%m,atom2basbas_len(i_atom_l),n,'d_ovlp3fn_SR for atom')
		    call aims_allocate(d_ovlp3fn_SR(i_atom_l,1)%n,n,'norm d_ovlp3fn_SR for atom')
		    call aims_allocate(d_ovlp3fn_SR(i_atom_l,2)%n,n,'norm d_ovlp3fn_SR for atom')
		    call aims_allocate(d_ovlp3fn_SR(i_atom_l,3)%n,n,'norm d_ovlp3fn_SR for atom')
		  endif
		  if(AS_init_stress) then
		    do AS_index = 1, AS_components, 1
		      call aims_allocate(AS_d_ovlp3fn_SR(i_atom_l,AS_index)%m,atom2basbas_len(i_atom_l),n,'AS_d_ovlp3fn_SR for atom')
		      call aims_allocate(AS_d_ovlp3fn_SR(i_atom_l,AS_index)%n,n,'norm AS_d_ovlp3fn_SR for atom')
		    enddo
		  endif
      endif
    enddo

    !---------------------------------------------------------------------------

    ! Get a linear list of all pairs

    allocate(pair_list(3,my_n_atom_pairs), stat=info)
    call check_allocation(info, 'pair_list', func, 3, my_n_atom_pairs)
    allocate(my_atom_pair_s(my_n_atoms), stat=info)
    call check_allocation(info, 'my_atom_pair_s', func, my_n_atoms)
    allocate(my_atom_pair_e(my_n_atoms), stat=info)
    call check_allocation(info, 'my_atom_pair_e', func, my_n_atoms)

    n = 0
    do i_my_atom = 1, my_n_atoms
      i_atom_l = my_atom_list(i_my_atom)
      my_atom_pair_s(i_my_atom) = n+1
      do i_cell = 1, n_cells_bvk
        do i_atom_r = 1, n_atoms2
          if(my_pair_flag_bvk(i_atom_l,i_atom_r,i_cell)) then
            n = n+1
            pair_list(1, n) = i_atom_l
            pair_list(2, n) = i_atom_r
            pair_list(3, n) = i_cell
          endif
        enddo
      enddo
      my_atom_pair_e(i_my_atom) = n
    enddo

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Calculate ovlp3fn
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    call mpi_barrier(mpi_comm_global, mpierr)
    ttt0 = mpi_wtime()
    if(use_lc_wpbeh .and. hse_omega_hf /= 0.d0 .and. hybrid_coeff .eq. 0.d0) then !.and. hybrid_coeff .eq. 0.d0
    	lc_wpbeh_lr_run = .true.
    endif
123 if(use_lc_wpbeh .and. hse_omega_hf /= 0.d0 .and. lc_wpbeh_lr_run) then
      write(info_str,'(A)') '  --- Initializing LVL triples, type LR'
      call localorb_info ( info_str,use_unit )
      ovlp_type_bare_or_hse_coul = OVLP_TYPE_LR
      call initialize_lvl_triples(OVLP_TYPE_LR)
      write(info_str,'(A)') '  --- Initializing Coulomb matrix calculations, type LR'
      call localorb_info ( info_str,use_unit )
      call initialize_tb_auxmat(1, OVLP_TYPE_LR)
      ovlp_type_bare_or_hse_coul = OVLP_TYPE_HSE
	else if(use_hse .and. hse_omega_hf /= 0.d0) then
      write(info_str,'(A)') '  --- Initializing LVL triples, type HSE'
      call localorb_info ( info_str,use_unit )
      call initialize_lvl_triples(OVLP_TYPE_HSE)
      write(info_str,'(A)') '  --- Initializing Coulomb matrix calculations, type HSE'
      call localorb_info ( info_str,use_unit )
      call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
    else
      write(info_str,'(A)') '  --- Initializing LVL triples'
      call localorb_info ( info_str,use_unit )
      call initialize_lvl_triples(OVLP_TYPE_COULOMB)
      write(info_str,'(A)') '  --- Initializing Coulomb matrix calculations'
      call localorb_info ( info_str,use_unit )
      call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
    endif
    call mpi_barrier(mpi_comm_global, mpierr)

    write(info_str,'(A,F12.3)') '  Finished with initialize_tb_auxmat/lvl_triples, time: ', mpi_wtime()-ttt0
    call localorb_info ( info_str,use_unit )

    write(info_str,'(A)') '  --- Calculating OVLP3FN'
    call localorb_info ( info_str,use_unit )
    ttt0 = mpi_wtime()

    n_super = product(2*number_of_super_cells(1:3)+1)
    allocate(cell_list(n_super), stat=info)
    call check_allocation(info, 'cell_list', func, n_super)
    allocate(Rvecs(3,n_super))
    call check_allocation(info, 'Rvecs', func, 3, n_super)
    if(.not.flag_max_coeff_3fn) then
       do while(max_n_basis_sp*max_n_basis_sp*max_n_basbas_sp*MIN(n_super,n_max_coeff_3fn)*16 .gt. 5.e7 .and. n_max_coeff_3fn .gt. 1)
           n_max_coeff_3fn = n_max_coeff_3fn -1
       enddo
    endif
    call aims_allocate(coeff_3fn,max_n_basis_sp2,max_n_basis_sp2,max_n_basbas_sp,2,MIN(n_super,n_max_coeff_3fn),'+coeff_3fn')
    if(init_calc_deriv) then
      call aims_allocate(d_coeff_3fn,max_n_basis_sp2,max_n_basis_sp2,max_n_basbas_sp,2,3,MIN(n_super,n_max_coeff_3fn),'+d_coeff_3fn')
    else
      call aims_allocate(d_coeff_3fn,1,1,1,1,1,1,'d_coeff_3fn')
    endif

    do i_my_atom = 1, my_n_atoms
      i_atom_l = my_atom_list(i_my_atom)

      ovlp3fn(i_atom_l)%m(:,:) = 0. ! Entries are ADDED below
      if(init_calc_deriv) then
        d_ovlp3fn(i_atom_l,1)%m(:,:) = 0. ! dito
        d_ovlp3fn(i_atom_l,2)%m(:,:) = 0. ! dito
        d_ovlp3fn(i_atom_l,3)%m(:,:) = 0. ! dito
      endif
      if(AS_init_stress) then
        do AS_index = 1, AS_components, 1
          AS_d_ovlp3fn(i_atom_l,AS_index)%m(:,:) = 0.0d0
        enddo
      endif

      do i_atom_r=1,n_atoms2

        ! Get all relevant atom pairs for atoms (i_atom_l,i_atom_r)

        Min_Rvec_length = 100.d0
        n_pairs = 0
        do i1 = -number_of_super_cells(1), number_of_super_cells(1)
          do i2 = -number_of_super_cells(2), number_of_super_cells(2)
            do i3 = -number_of_super_cells(3), number_of_super_cells(3)
              if(pair_flag_full(i_atom_l,i_atom_r,i1,i2,i3)) then

                i_cell = get_bvk_cell_idx( (/ i1, i2, i3 /) )

                if(.not.my_pair_flag_bvk(i_atom_l,i_atom_r,i_cell)) cycle

                n_pairs = n_pairs+1

                Rvecs(:,n_pairs) = coords(:,vb_atom2atom(i_atom_r)) - coords(:,vb_atom2atom(i_atom_l)) + matmul(lattice_vector, (/ i1, i2, i3 /))
                if(i1 == 0 .and. i2 == 0 .and. i3 == 0) i_cell = 0 ! Need to flag the center cell specifically !!!
                cell_list(n_pairs) = i_cell
! XR: debugging
!                Rvec_length = sqrt(Rvecs(1,n_pairs)**2 + Rvecs(2,n_pairs)**2 + Rvecs(3,n_pairs)**2)
!                if(Min_Rvec_length.gt.Rvec_length) then
!                   Min_Rvec_length = Rvec_length
!                   write(use_unit,'(3I4,f16.8)') i1, i2, i3, Rvec_length
!                   write(use_unit,'(3I4,f16.8)') i_atom_r, i_atom_l
!                   write(use_unit,'(3f16.8)') coords(:,i_atom_r)
!                   write(use_unit,'(3f16.8)') coords(:,i_atom_l)
!                endif
              endif
            enddo
          enddo
        enddo

        ! Calculate coefficients for all pairs with a maximum of max_coeff_3fn per call to get_pairwise_coeff_3fn

        do i_start = 1, n_pairs, n_max_coeff_3fn

          n_pairs_cur = min(n_max_coeff_3fn, n_pairs-i_start+1)

          ! VB: The following subroutine is called often, and the ovlp3fn setup is still very significant.
          !
          !     A great deal of allocations and deallocations happen inside this routine.
          !     They, at least, should be moved to the outside.
          !
          !     Also, a number of subroutine calls in there do reshape arrays automatically.
          !     What is the time impact of those extra copies?
          !
          !     I am writing down these possible optimizations here in case there is ever
          !     time to profile and optimize this part of the code.


          call get_pairwise_coeff_3fn_vb(i_atom_l, i_atom_r, species(vb_atom2atom(i_atom_l)), species(vb_atom2atom(i_atom_r)), &
                                         atom2basis_len2, atom2vb_basis_off, n_atoms2, &
                                         n_pairs_cur, Rvecs(1,i_start), &
                                         coeff_3fn, d_coeff_3fn, init_calc_deriv)
          do n=1,n_pairs_cur

            i_cell = cell_list(n+i_start-1)
            if(i_cell==0) i_cell = 1 ! Real center cell
            nbb1 = atom2basbas_len(i_atom_l)

            k = pair_offset(i_atom_l,i_atom_r,i_cell)

            do jj = 1, atom2basis_len2(i_atom_r)
              do ii = 1, atom2basis_len2(i_atom_l)
                k = k+1
                if(cell_list(n+i_start-1)==0 .and. vb_atom2atom(i_atom_l)==vb_atom2atom(i_atom_r)) then
                  ovlp3fn(i_atom_l)%m(1:nbb1,k) = ovlp3fn(i_atom_l)%m(1:nbb1,k) + 0.5*coeff_3fn(ii,jj,1:nbb1,1,n)
                  ! Contribution to d_ovlp3fn is 0
                else
                  ovlp3fn(i_atom_l)%m(1:nbb1,k) = ovlp3fn(i_atom_l)%m(1:nbb1,k) + coeff_3fn(ii,jj,1:nbb1,1,n)
                  if(init_calc_deriv) then
                    d_ovlp3fn(i_atom_l,1)%m(1:nbb1,k) = d_ovlp3fn(i_atom_l,1)%m(1:nbb1,k) + d_coeff_3fn(ii,jj,1:nbb1,1,1,n)
                    d_ovlp3fn(i_atom_l,2)%m(1:nbb1,k) = d_ovlp3fn(i_atom_l,2)%m(1:nbb1,k) + d_coeff_3fn(ii,jj,1:nbb1,1,2,n)
                    d_ovlp3fn(i_atom_l,3)%m(1:nbb1,k) = d_ovlp3fn(i_atom_l,3)%m(1:nbb1,k) + d_coeff_3fn(ii,jj,1:nbb1,1,3,n)
                  endif

                  !FK: Multiply nuclear distances and derivatives with respect to nuclear coordinates
                  if (AS_init_stress) then
                    do AS_index = 1, AS_components, 1
                      AS_d_ovlp3fn(i_atom_l,AS_index)%m(1:nbb1,k) = AS_d_ovlp3fn(i_atom_l,AS_index)%m(1:nbb1,k) + &
                        d_coeff_3fn(ii,jj,1:nbb1,1,AS_l_index(AS_index),n) * Rvecs(AS_m_index(AS_index),i_start-1+n)
                    end do
                  end if
                endif
              enddo
            enddo

          enddo
        enddo

      enddo
    enddo
    call mpi_barrier(mpi_comm_global, mpierr)

    write(info_str,'(A,F12.3)') '  --- Done calculating OVLP3FN, time: ',mpi_wtime()-ttt0
    call localorb_info ( info_str,use_unit )

    call cleanup_lvl_triples()
    write(info_str,'(A,2I8)') 'Compare n_pairs with n_super',n_pairs,n_super
    call localorb_info ( info_str,use_unit )

    deallocate(cell_list)
    deallocate(Rvecs)
    call aims_deallocate( coeff_3fn,     "coeff_3fn" )
    call aims_deallocate( d_coeff_3fn, "d_coeff_3fn" )
    if (.not. use_lc_wpbeh .or. (use_lc_wpbeh .and. lc_wpbeh_lr_run)) then
    	deallocate(pair_flag_full)
    endif

    !---------------------------------------------------------------------------

    ! Get maximum norm of ovlp3fn for every atom pair

    ttt0 = mpi_wtime()
    allocate(max_norm_ovlp3fn(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
    call check_allocation(info, 'max_norm_ovlp3fn', func, n_atoms2, n_atoms2, n_cells_bvk)
    max_norm_ovlp3fn = 0
    allocate(max_norm_d_ovlp3fn(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
    call check_allocation(info, 'max_norm_d_ovlp3fn', func, n_atoms2, n_atoms2, n_cells_bvk)
    max_norm_d_ovlp3fn = 0
    if(AS_init_stress) then
      allocate(AS_max_norm_d_ovlp3fn(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
      call check_allocation(info, 'AS_max_norm_d_ovlp3fn', func, n_atoms2, n_atoms2, n_cells_bvk)
      AS_max_norm_d_ovlp3fn = 0.0d0
    endif

    if (.not.(use_mpi_in_place)) then
       allocate(max_norm_ovlp3fn_rcv(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
       call check_allocation(info, 'max_norm_ovlp3fn_rcv', func, n_atoms2, n_atoms2, n_cells_bvk)
       max_norm_ovlp3fn_rcv = 0
       allocate(max_norm_d_ovlp3fn_rcv(n_atoms2,n_atoms2,n_cells_bvk))
       call check_allocation(info, 'max_norm_d_ovlp3fn_rcv', func, n_atoms2, n_atoms2, n_cells_bvk)
       max_norm_d_ovlp3fn_rcv = 0
       if(AS_init_stress) then
         allocate(AS_max_norm_d_ovlp3fn_rcv(n_atoms2,n_atoms2,n_cells_bvk))
         call check_allocation(info, 'AS_max_norm_d_ovlp3fn_rcv', func, n_atoms2, n_atoms2, n_cells_bvk)
         AS_max_norm_d_ovlp3fn_rcv = 0.0d0
       endif
    endif

    do i_atom_pair = 1, my_n_atom_pairs

      i_atom_l = pair_list(1,i_atom_pair)
      i_atom_r = pair_list(2,i_atom_pair)
      i_cell   = pair_list(3,i_atom_pair)

      k = pair_offset(i_atom_l,i_atom_r,i_cell)
      n = atom2basis_len2(i_atom_l)*atom2basis_len2(i_atom_r)

      do i = k+1, k+n
        ovlp3fn(i_atom_l)%n(i) = sqrt(sum(ovlp3fn(i_atom_l)%m(:,i)**2))
        max_norm_ovlp3fn(i_atom_l,i_atom_r,i_cell) = &
          max(max_norm_ovlp3fn(i_atom_l,i_atom_r,i_cell),ovlp3fn(i_atom_l)%n(i))
      enddo

      if(init_calc_deriv) then
        do i = k+1, k+n
          do j = 1, 3
            d_ovlp3fn(i_atom_l,j)%n(i) = sqrt(sum(d_ovlp3fn(i_atom_l,j)%m(:,i)**2))
            max_norm_d_ovlp3fn(i_atom_l,i_atom_r,i_cell) = &
              max(max_norm_d_ovlp3fn(i_atom_l,i_atom_r,i_cell), d_ovlp3fn(i_atom_l,j)%n(i))
          enddo
        enddo
      endif

      if (AS_init_stress) then
        do i = k+1, k+n
          do AS_index = 1, AS_components, 1
            AS_d_ovlp3fn(i_atom_l,AS_index)%n(i) = sqrt(sum(AS_d_ovlp3fn(i_atom_l,AS_index)%m(:,i)**2))
            AS_max_norm_d_ovlp3fn(i_atom_l,i_atom_r,i_cell) = &
              max(AS_max_norm_d_ovlp3fn(i_atom_l,i_atom_r,i_cell), AS_d_ovlp3fn(i_atom_l,AS_index)%n(i))
          end do
        end do
      end if

    enddo

    if (use_mpi) then
       if (.not.(use_mpi_in_place)) then
          call mpi_allreduce(max_norm_ovlp3fn,max_norm_ovlp3fn_rcv,size(max_norm_ovlp3fn),&
                             mpi_real8,mpi_max,mpi_comm_global,mpierr)
          max_norm_ovlp3fn = max_norm_ovlp3fn_rcv
          call mpi_allreduce(max_norm_d_ovlp3fn,max_norm_d_ovlp3fn_rcv,size(max_norm_d_ovlp3fn),&
                             mpi_real8,mpi_max,mpi_comm_global,mpierr)
          max_norm_d_ovlp3fn = max_norm_d_ovlp3fn_rcv
          if (AS_init_stress) then
             call mpi_allreduce(AS_max_norm_d_ovlp3fn,AS_max_norm_d_ovlp3fn_rcv,size(AS_max_norm_d_ovlp3fn),&
                                mpi_real8,mpi_max,mpi_comm_global,mpierr)
             AS_max_norm_d_ovlp3fn = AS_max_norm_d_ovlp3fn_rcv
          end if
       else
          call mpi_allreduce(mpi_in_place,max_norm_ovlp3fn,size(max_norm_ovlp3fn),&
                             mpi_real8,mpi_max,mpi_comm_global,mpierr)
          call mpi_allreduce(mpi_in_place,max_norm_d_ovlp3fn,size(max_norm_d_ovlp3fn),&
                             mpi_real8,mpi_max,mpi_comm_global,mpierr)
          if (AS_init_stress) then
             call mpi_allreduce(mpi_in_place,AS_max_norm_d_ovlp3fn,size(AS_max_norm_d_ovlp3fn),&
                                mpi_real8,mpi_max,mpi_comm_global,mpierr)
          end if
       endif
    end if



    ! Warning placed here by VB.
    !
    ! I encountered a buggy MPI library that simply returned all values of max_norm_ovlp3fn to be zero
    ! AFTER the mpi_allreduce call.
    ! If this happens, the exx energy is zero below.
    ! I did not place a stop here, but I did place a warning message below after the exx_ene variable is
    ! finally computed. If this is a more common problem with MPI libraries, we should place a check right here.
    !
    ! Remark added by A.Marek
    !
    ! the above written faulty behaviour should be gone now, after no
    ! "MPI_IN_PLACE" is used anymore. It is known that the usage of MPI_IN_PLACE
    ! can lead, with some MPI implementations, to wrong results in the reduction
    ! However, a performance degradation might be possible. You can switch back
    ! to the old behaviour (if it works correctly)

    ! Get maximum norm of ovlp3fn per left atom

    allocate(max_norm_ovlp3fn_per_latom(n_atoms2), stat=info)
    call check_allocation(info, 'max_norm_ovlp3fn_per_latom', func, n_atoms2)
    max_norm_ovlp3fn_per_latom(:) = 0

    allocate(max_norm_d_ovlp3fn_per_latom(n_atoms2), stat=info)
    call check_allocation(info, 'max_norm_d_ovlp3fn_per_latom', func, n_atoms2)
    max_norm_d_ovlp3fn_per_latom(:) = 0

    if(AS_init_stress) then
      allocate(AS_max_norm_d_ovlp3fn_per_latom(n_atoms2), stat=info)
      call check_allocation(info, 'AS_max_norm_d_ovlp3fn_per_latom', func, n_atoms2)
      AS_max_norm_d_ovlp3fn_per_latom(:) = 0.0d0
    endif

    do i_atom_l = 1, n_atoms2
      max_norm_ovlp3fn_per_latom(i_atom_l)   = maxval(max_norm_ovlp3fn  (i_atom_l,:,:))
      max_norm_d_ovlp3fn_per_latom(i_atom_l) = maxval(max_norm_d_ovlp3fn(i_atom_l,:,:))
      if (AS_init_stress) then
        AS_max_norm_d_ovlp3fn_per_latom(i_atom_l) = maxval(AS_max_norm_d_ovlp3fn(i_atom_l,:,:))
      end if
    enddo

    write(info_str,'(A,F12.3)') '  --- Done getting norm of OVLP3FN, time: ', mpi_wtime()-ttt0
    call localorb_info ( info_str,use_unit )

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Calculate Coulomb matrices
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    write(info_str,'(A)') '  --- Calculating Coulomb matrices'
    call localorb_info ( info_str,use_unit )

    allocate(coul_mat_store(n_atoms, n_atoms, n_cells_bvk), stat=info)
    call check_allocation(info, 'coul_mat_store', func, n_atoms, n_atoms, n_cells_bvk)
    if(init_calc_deriv) then
      allocate(d_coul_mat_store(n_atoms, n_atoms, n_cells_bvk, 3), stat=info)
      call check_allocation(info, 'd_coul_mat_store', func, n_atoms, n_atoms, n_cells_bvk, 3)
    endif
    if(AS_init_stress) then
      allocate(AS_d_coul_mat_store(n_atoms, n_atoms, n_cells_bvk, AS_components), stat=info)
      call check_allocation(info, 'AS_d_coul_mat_store', func, n_atoms, n_atoms, n_cells_bvk, AS_components)
    endif
    n_coulmat_elems = 0

    allocate(coul_mat_norm(n_atoms,n_atoms,n_cells_bvk), stat=info)
    call check_allocation(info, 'coul_mat_norm', func, n_atoms, n_atoms, n_cells_bvk)
    coul_mat_norm = 0.
    allocate(d_coul_mat_norm(n_atoms,n_atoms,n_cells_bvk), stat=info)
    call check_allocation(info, 'd_coul_mat_norm', func, n_atoms, n_atoms, n_cells_bvk)
    d_coul_mat_norm = 0.
    allocate(d_coul_mat_norm_3(n_atoms,n_atoms,n_cells_bvk,3), stat=info)
    call check_allocation(info, 'd_coul_mat_norm_3', func, n_atoms, n_atoms, n_cells_bvk, 3)
    d_coul_mat_norm_3 = 0.

    if(AS_init_stress) then
      allocate(AS_d_coul_mat_norm(n_atoms,n_atoms,n_cells_bvk), stat=info)
      call check_allocation(info, 'AS_d_coul_mat_norm', func, n_atoms, n_atoms, n_cells_bvk)
      AS_d_coul_mat_norm = 0.0d0
      allocate(AS_d_coul_mat_norm_components(n_atoms,n_atoms,n_cells_bvk,AS_components), stat=info)
      call check_allocation(info, 'AS_d_coul_mat_norm_components', func, n_atoms, n_atoms, n_cells_bvk, AS_components)
      AS_d_coul_mat_norm_components = 0.0d0
    endif

    if(n_periodic==0) then
      call aims_allocate(aux_mat, max_n_basbas_sp, max_n_basbas_sp, 1, '+aux_mat')
      if(init_calc_deriv) then
        call aims_allocate(d_aux_mat, max_n_basbas_sp, max_n_basbas_sp, 3, 1, '+d_aux_mat')
      else
        call aims_allocate(d_aux_mat, 1, 1, 1, 1, 'd_aux_mat')
      endif
    endif

    ttt0 = mpi_wtime()

    if(n_periodic > 0) then

      ! For the periodic non-HSE case OVLP_TYPE_CUT is used instead of OVLP_TYPE_COULOMB,
      ! so reinitalize the auxmat stuff if necessary
      if(use_hse .and. hse_omega_hf /= 0.d0) then
        ! Nothing to do
        continue
      else
        call deallocate_tb_auxmat()
        call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
        ! NB: integrate_ovlp3fn_p0 was using OVLP_TYPE_COULOMB here
      endif

    endif

    cm_calculated = .false.
    do i_my_atom = 1, my_n_atoms
      j_atom = my_atom_list(i_my_atom)
      do i_atom = 1, n_atoms2
        nbb1 = atom2basbas_len(i_atom)
        nbb2 = atom2basbas_len(j_atom)
        if(.not. cm_calculated(vb_atom2atom(i_atom),vb_atom2atom(j_atom))) then
        if(my_tasks_per_atom > 1) then
         n_bb_s = my_cm_bb_s
         n_bb_e = my_cm_bb_e
        else
          n_bb_s = 1
          n_bb_e = atom2basbas_len(j_atom)
        endif

        Rvec(:) = coords(:,vb_atom2atom(j_atom)) - coords(:,vb_atom2atom(i_atom))
        if(n_periodic > 0) then
          ! For a large number of cells, the full aux_mat (max_n_basbas_sp**2 * n_cells_bvk entries)
          ! can get quite expensive, so we allocate only what is really needed
          call aims_allocate(aux_mat, 1, nbb1, n_bb_s, n_bb_e, 1, my_cm_cell_num, 'aux_mat')
          if(init_calc_deriv) then
           call aims_allocate(d_aux_mat, 1, nbb1, n_bb_s, n_bb_e, 1, 3, 1, my_cm_cell_num, 'd_aux_mat')
          else
           call aims_allocate(d_aux_mat, 1, 1, 1, 1, 'd_aux_mat')
          endif
          if(AS_init_stress) then
           call aims_allocate(AS_d_aux_mat, 1, nbb1, n_bb_s, n_bb_e, 1, AS_components, 1, my_cm_cell_num, 'AS_d_aux_mat')
          else
           call aims_allocate(AS_d_aux_mat, 1, 1, 1, 1, 'AS_d_aux_mat')
          endif
          if(use_symmetry_reduced_spg)then
            call calculate_realspace_coulomb_sym(species(vb_atom2atom(i_atom)), species(vb_atom2atom(j_atom)), Rvec, &
                                           aux_mat, d_aux_mat, AS_d_aux_mat, nbb1, n_bb_s, n_bb_e,&
                                           my_cm_cell_num,init_calc_deriv,&
                                           AS_init_stress,AS_l_index,AS_m_index,&
                                           my_cm_cell_start,my_cm_cell_inc,&
                                           my_atom_id, my_tasks_per_atom, &
			                   mpi_comm_atom)
	  else
            call calculate_realspace_coulomb(species(vb_atom2atom(i_atom)), species(vb_atom2atom(j_atom)), Rvec, &
                                           aux_mat, d_aux_mat, AS_d_aux_mat, nbb1, n_bb_s, n_bb_e)
	  endif
        else
          call fast_calculate_tb_auxmat(species(vb_atom2atom(i_atom)), species(vb_atom2atom(j_atom)), Rvec, &
                                        aux_mat(:,:,1), d_aux_mat(:,:,:,1), init_calc_deriv)
        endif

        n_cnt = 0
        do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc
          n_cnt = n_cnt + 1
          call compress_matrix(aux_mat(1:nbb1,n_bb_s:n_bb_e,n_cnt),coul_mat_store(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell))
          coul_mat_norm(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell) = sum(aux_mat(1:nbb1,n_bb_s:n_bb_e,n_cnt)**2)

          if(init_calc_deriv) then
            call compress_matrix(d_aux_mat(1:nbb1,n_bb_s:n_bb_e,1,n_cnt),d_coul_mat_store(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,1))
            call compress_matrix(d_aux_mat(1:nbb1,n_bb_s:n_bb_e,2,n_cnt),d_coul_mat_store(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,2))
            call compress_matrix(d_aux_mat(1:nbb1,n_bb_s:n_bb_e,3,n_cnt),d_coul_mat_store(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,3))
            d_coul_mat_norm_3(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,1) = sum(d_aux_mat(1:nbb1,n_bb_s:n_bb_e,1,n_cnt)**2)
            d_coul_mat_norm_3(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,2) = sum(d_aux_mat(1:nbb1,n_bb_s:n_bb_e,2,n_cnt)**2)
            d_coul_mat_norm_3(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,3) = sum(d_aux_mat(1:nbb1,n_bb_s:n_bb_e,3,n_cnt)**2)
          endif

          if (AS_init_stress) then
            do AS_index = 1, AS_components, 1
              call compress_matrix &
              ( AS_d_aux_mat(1:nbb1,n_bb_s:n_bb_e,AS_index,n_cnt),AS_d_coul_mat_store(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,AS_index))

              AS_d_coul_mat_norm_components(vb_atom2atom(i_atom),vb_atom2atom(j_atom),i_cell,AS_index) = &
              sum(AS_d_aux_mat(1:nbb1,n_bb_s:n_bb_e,AS_index,n_cnt)**2)
            end do
          end if
        enddo

        if(n_periodic > 0) then
          call aims_deallocate( aux_mat,           "aux_mat" )
          call aims_deallocate( d_aux_mat,       "d_aux_mat" )
          call aims_deallocate( AS_d_aux_mat, "AS_d_aux_mat" )
        endif
        cm_calculated(vb_atom2atom(i_atom),vb_atom2atom(j_atom)) = .true.
        endif
      enddo
    enddo


      !determine the maximum number of rows/cols to allocate the temp arrays only once per step
      !colcount = maxval(sp2n_basis_sp)
      !call sync_integer(colcount,mpi_comm_global)
      !maxrow = 0
      !maxcol = 0
      !do index_1 = lbound(coul_mat_store,1), ubound(coul_mat_store,1), 1
      !   do index_2 = lbound(coul_mat_store,2), ubound(coul_mat_store,2), 1
      !      do index_3 = lbound(coul_mat_store,3), ubound(coul_mat_store,3), 1
      !         if (coul_mat_store(index_1,index_2,index_3)%n_rows > maxrow) then
      !            maxrow = coul_mat_store(index_1,index_2,index_3)%n_rows
      !            peakrow(:) = (/ index_1,index_2,index_3 /)
      !         endif
      !         if (coul_mat_store(index_1,index_2,index_3)%n_cols > maxcol) then
      !            maxcol = coul_mat_store(index_1,index_2,index_3)%n_cols
      !            peakcol(:) = (/ index_1,index_2,index_3 /)
      !         endif
      !      enddo
      !   enddo
      !enddo
      !call mult_coul_mat_left(colcount,coul_mat_store(peakrow(1),peakrow(2),peakrow(3)), coul_mat_norm, 1, coul_mat_norm, 1, 1)
      !call mult_coul_mat_left(colcount,coul_mat_store(peakcol(1),peakcol(2),peakcol(3)), coul_mat_norm, 1, coul_mat_norm, 1, 2)
      !call mult_coul_mat_left(colcount,coul_mat_store(peakrow(1),peakrow(2),peakrow(3)), coul_mat_norm, 1, coul_mat_norm, 1, 3)

    if(n_periodic==0) then
      call aims_deallocate( aux_mat,     "aux_mat" )
      call aims_deallocate( d_aux_mat, "d_aux_mat" )
    endif

    call sync_vector(coul_mat_norm, size(coul_mat_norm))
    coul_mat_norm(:,:,:) = sqrt(coul_mat_norm(:,:,:))
    if(init_calc_deriv) then
      call sync_vector(d_coul_mat_norm_3, size(d_coul_mat_norm_3))
      d_coul_mat_norm(:,:,:) = max(d_coul_mat_norm_3(:,:,:,1),d_coul_mat_norm_3(:,:,:,2),d_coul_mat_norm_3(:,:,:,3))
      d_coul_mat_norm(:,:,:) = sqrt(d_coul_mat_norm(:,:,:))
    endif
    deallocate(d_coul_mat_norm_3)

    if (AS_init_stress) then
      call sync_vector(AS_d_coul_mat_norm_components, size(AS_d_coul_mat_norm_components))

      do AS_index = 1, AS_components, 1
        AS_d_coul_mat_norm(:,:,:) = max(AS_d_coul_mat_norm(:,:,:),AS_d_coul_mat_norm_components(:,:,:,AS_index))
      end do
      AS_d_coul_mat_norm(:,:,:) = sqrt(AS_d_coul_mat_norm(:,:,:))

      if (allocated(AS_d_coul_mat_norm_components)) deallocate(AS_d_coul_mat_norm_components)
    end if

    !The following Part stores the calculated results in new arrays and runs the code again to calculate ovlp3n and coulomb_matr
    ! but this time for the LR-Part
    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
      write(info_str,'(A)') '  --- Done for short range, now for long range'
      call localorb_info ( info_str,use_unit )
      lc_wpbeh_lr_run=.true.
      ! switch arrays to store SR Part and deallocate those which will be reallocated again above after the goto.
      ! This is not really nice coding, but i haven't figuered out another easy way of doing this without copying the complete code.
      ! I think this is the best compromise at the moment. Maybe the deallocation can be done in a separat routine at some point.
      ovlp3fn_SR=ovlp3fn
      if(init_calc_deriv) d_ovlp3fn_SR=d_ovlp3fn
      if(AS_init_stress) AS_d_ovlp3fn_SR=AS_d_ovlp3fn

      allocate(max_norm_ovlp3fn_SR(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
      call check_allocation(info, 'max_norm_ovlp3fn_SR', func, n_atoms2, n_atoms2, n_cells_bvk)
      max_norm_ovlp3fn_SR=max_norm_ovlp3fn
      if(allocated(max_norm_ovlp3fn)) deallocate(max_norm_ovlp3fn)

      allocate(max_norm_d_ovlp3fn_SR(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
      call check_allocation(info, 'max_norm_d_ovlp3fn_SR', func, n_atoms2, n_atoms2, n_cells_bvk)
      max_norm_d_ovlp3fn_SR=max_norm_d_ovlp3fn
      if(allocated(max_norm_d_ovlp3fn)) deallocate(max_norm_d_ovlp3fn)

      if (AS_init_stress) then
     		allocate(AS_max_norm_d_ovlp3fn_SR(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
                call check_allocation(info, 'AS_max_norm_d_ovlp3fn_SR', func, n_atoms2, n_atoms2, n_cells_bvk)
      	AS_max_norm_d_ovlp3fn_SR=AS_max_norm_d_ovlp3fn
      	if(allocated(AS_max_norm_d_ovlp3fn)) deallocate(AS_max_norm_d_ovlp3fn)
      end if

      allocate(max_norm_ovlp3fn_per_latom_SR(n_atoms2), stat=info)
      call check_allocation(info, 'max_norm_ovlp3fn_per_latom_SR', func, n_atoms2)
      max_norm_ovlp3fn_per_latom_SR=max_norm_ovlp3fn_per_latom
      if(allocated(max_norm_ovlp3fn_per_latom)) deallocate(max_norm_ovlp3fn_per_latom)

      allocate(max_norm_d_ovlp3fn_per_latom_SR(n_atoms2), stat=info)
      call check_allocation(info, 'max_norm_d_ovlp3fn_per_latom_SR', func, n_atoms2)
      max_norm_d_ovlp3fn_per_latom_SR=max_norm_d_ovlp3fn_per_latom
      if(allocated(max_norm_d_ovlp3fn_per_latom)) deallocate(max_norm_d_ovlp3fn_per_latom)

      if (AS_init_stress) then
     		allocate(AS_max_norm_d_ovlp3fn_per_latom(n_atoms2), stat=info)
                call check_allocation(info, 'AS_max_norm_d_ovlp3fn_per_latom', func, n_atoms2)
      	AS_max_norm_d_ovlp3fn_per_latom_SR=AS_max_norm_d_ovlp3fn_per_latom
      	if(allocated(AS_max_norm_d_ovlp3fn_per_latom)) deallocate(AS_max_norm_d_ovlp3fn_per_latom)
      end if

		if (.not.(use_mpi_in_place)) then
			allocate(max_norm_ovlp3fn_rcv_SR(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
                        call check_allocation(info, 'max_norm_ovlp3fn_rcv_SR', func, n_atoms2, n_atoms2, n_cells_bvk)
			max_norm_ovlp3fn_rcv_SR=max_norm_ovlp3fn_rcv
		   if(allocated(max_norm_ovlp3fn_rcv)) deallocate(max_norm_ovlp3fn_rcv)

		   allocate(max_norm_d_ovlp3fn_rcv_SR(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
                   call check_allocation(info, 'max_norm_d_ovlp3fn_rcv_SR', func, n_atoms2, n_atoms2, n_cells_bvk)
			max_norm_d_ovlp3fn_rcv_SR=max_norm_d_ovlp3fn_rcv
	   	if(allocated(max_norm_d_ovlp3fn_rcv)) deallocate(max_norm_d_ovlp3fn_rcv)

		   if (AS_init_stress) then
		   	allocate(AS_max_norm_d_ovlp3fn_rcv_SR(n_atoms2,n_atoms2,n_cells_bvk), stat=info)
                        call check_allocation(info, 'AS_max_norm_d_ovlp3fn_rcv_SR', func, n_atoms2, n_atoms2, n_cells_bvk)
		   	AS_max_norm_d_ovlp3fn_rcv_SR=AS_max_norm_d_ovlp3fn_rcv
		   	if(allocated(AS_max_norm_d_ovlp3fn_rcv)) deallocate(AS_max_norm_d_ovlp3fn_rcv)
		   end if
		endif

      allocate(coul_mat_store_SR(n_atoms, n_atoms, n_cells_bvk), stat=info)
      call check_allocation(info, 'coul_mat_store_SR', func, n_atoms, n_atoms, n_cells_bvk)
      coul_mat_store_SR=coul_mat_store
      do k = lbound(coul_mat_store,3), ubound(coul_mat_store,3)
        do j = lbound(coul_mat_store,2), ubound(coul_mat_store,2)
          do i = lbound(coul_mat_store,1), ubound(coul_mat_store,1)
            if(allocated(coul_mat_store(i,j,k)%row_idx)) call aims_deallocate( coul_mat_store(i,j,k)%row_idx, 'coul_mat_store%row_idx')
            if(allocated(coul_mat_store(i,j,k)%col_idx)) call aims_deallocate( coul_mat_store(i,j,k)%col_idx, 'coul_mat_store%col_idx')
            if(allocated(coul_mat_store(i,j,k)%mat    )) call aims_deallocate( coul_mat_store(i,j,k)%mat,         'coul_mat_store%mat')
          enddo
        enddo
      enddo
      deallocate(coul_mat_store)

      if(init_calc_deriv) then
        allocate(d_coul_mat_store_SR(n_atoms, n_atoms, n_cells_bvk, 3), stat=info)
        call check_allocation(info, 'd_coul_mat_store_SR', func, n_atoms, n_atoms, n_cells_bvk, 3)
        d_coul_mat_store_SR=d_coul_mat_store
        if(allocated(d_coul_mat_store)) then
           do d = 1, 3
           do k = lbound(d_coul_mat_store,3), ubound(d_coul_mat_store,3)
           do j = lbound(d_coul_mat_store,2), ubound(d_coul_mat_store,2)
           do i = lbound(d_coul_mat_store,1), ubound(d_coul_mat_store,1)
             if(allocated(d_coul_mat_store(i,j,k,d)%row_idx)) &
                  call aims_deallocate( d_coul_mat_store(i,j,k,d)%row_idx, "d_coul_mat_store%row_idx" )
             if(allocated(d_coul_mat_store(i,j,k,d)%col_idx)) &
                  call aims_deallocate( d_coul_mat_store(i,j,k,d)%col_idx, "d_coul_mat_store%col_idx" )
             if(allocated(d_coul_mat_store(i,j,k,d)%mat    )) &
                  call aims_deallocate( d_coul_mat_store(i,j,k,d)%mat ,        "d_coul_mat_store%mat" )
           enddo
           enddo
           enddo
           enddo
           deallocate(d_coul_mat_store)
        endif
      end if

      if (AS_init_stress) then
        allocate(AS_d_coul_mat_store_SR(n_atoms, n_atoms, n_cells_bvk, AS_components), stat=info)
        call check_allocation(info, 'AS_d_coul_mat_store_SR', func, n_atoms, n_atoms, n_cells_bvk, AS_components)
        AS_d_coul_mat_store_SR=AS_d_coul_mat_store
        if(allocated(AS_d_coul_mat_store)) then
          do d = 1, AS_components, 1
          do k = lbound(AS_d_coul_mat_store,3), ubound(AS_d_coul_mat_store,3)
          do j = lbound(AS_d_coul_mat_store,2), ubound(AS_d_coul_mat_store,2)
          do i = lbound(AS_d_coul_mat_store,1), ubound(AS_d_coul_mat_store,1)
            if(allocated(AS_d_coul_mat_store(i,j,k,d)%row_idx)) &
                  call aims_deallocate( AS_d_coul_mat_store(i,j,k,d)%row_idx, "AS_d_coul_mat_store%row_idx" )
            if(allocated(AS_d_coul_mat_store(i,j,k,d)%col_idx)) &
                  call aims_deallocate( AS_d_coul_mat_store(i,j,k,d)%col_idx, "AS_d_coul_mat_store%col_idx" )
            if(allocated(AS_d_coul_mat_store(i,j,k,d)%mat    )) &
                  call aims_deallocate( AS_d_coul_mat_store(i,j,k,d)%mat,         "AS_d_coul_mat_store%mat" )
          enddo
          enddo
          enddo
          enddo
          deallocate(AS_d_coul_mat_store)
        endif
      end if

		 allocate(coul_mat_norm_SR(n_atoms,n_atoms,n_cells_bvk), stat=info)
                 call check_allocation(info, 'coul_mat_norm_SR', func, n_atoms, n_atoms, n_cells_bvk)
		 coul_mat_norm_SR=coul_mat_norm
		 if(allocated(coul_mat_norm)) deallocate(coul_mat_norm)

		 allocate(d_coul_mat_norm_SR(n_atoms,n_atoms,n_cells_bvk), stat=info)
                 call check_allocation(info, 'd_coul_mat_norm_SR', func, n_atoms, n_atoms, n_cells_bvk)
		 d_coul_mat_norm_SR=d_coul_mat_norm
		 if(allocated(d_coul_mat_norm)) deallocate(d_coul_mat_norm)

		 if (AS_init_stress) then
		 	allocate(AS_d_coul_mat_norm_SR(n_atoms,n_atoms,n_cells_bvk), stat=info)
                        call check_allocation(info, 'AS_d_coul_mat_norm_SR', func, n_atoms, n_atoms, n_cells_bvk)
			AS_d_coul_mat_norm_SR=AS_d_coul_mat_norm
			if(allocated(AS_d_coul_mat_norm)) deallocate(AS_d_coul_mat_norm)
		 end if

		 call deallocate_tb_auxmat()
		 call mpi_barrier(mpi_comm_global, mpierr)
		 ! now go back up and do the routine again but this time the LR will be stored in the arrays.
       goto 123
    endif

    write(info_str,'(A,F12.3)') '  Time to calculate Coulomb Matrices   : ', mpi_wtime()-ttt0
    call localorb_info ( info_str,use_unit )
    write(info_str, *) ' Number of Coulomb matrix elements:'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

    write(info_str,'(a,i5,3(a,i14))') '  | Task ',myid,' total: ',n_coulmat_elems(1), &
                         ', > threshold: ',n_coulmat_elems(2), &
                         ', stored: ',n_coulmat_elems(3)
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)

    write(info_str,'(A)') '  --- Done init_fock_matrix_calculations'
    call localorb_info ( info_str,use_unit )

    call deallocate_tb_auxmat()

    call aims_allocate(max_norm_tmp_old,n_cells_bvk,'max_norm_tmp_old')
    max_norm_tmp_old = 0.d0
    call aims_allocate(max_abs_dm_cols_old,n_atoms2,n_atoms2,n_cells_bvk,'max_abs_dm_cols_old')
    max_abs_dm_cols_old = 0.d0

! DB 10/02/13
    ! as promised above: we have to change back the global meaning of n_atoms
    n_atoms = n_atoms_save

    if(allocated(coeff_3fn))    call aims_deallocate( coeff_3fn,       "coeff_3fn" )
    if(allocated(d_coeff_3fn))  call aims_deallocate( d_coeff_3fn,   "d_coeff_3fn" )
    if(allocated(aux_mat))      call aims_deallocate( aux_mat,           "aux_mat" )
    if(allocated(d_aux_mat))    call aims_deallocate( d_aux_mat,       "d_aux_mat" )
    if(allocated(AS_d_aux_mat)) call aims_deallocate( AS_d_aux_mat, "AS_d_aux_mat" )

  contains

    integer function atom2basis_len(i_atom)
        integer :: i_atom
        atom2basis_len = sp2n_basis_sp(species(i_atom))
    end function
    integer function atom2basbas_len(i_atom)
        integer :: i_atom
        atom2basbas_len = sp2n_basbas_sp(species(vb_atom2atom(i_atom)))
    end function
    integer function atom2basbas_len1(i_atom)
        integer :: i_atom
        atom2basbas_len1 = sp2n_basbas_sp(species(i_atom))
    end function

  end subroutine init_fock_matrix_calculations
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/divide_atoms
  !  NAME
  !    divide_atoms
  !  SYNOPSIS

  subroutine divide_atoms(weight, task_of_atom, tasks_per_atom)

    !  PURPOSE
    !
    !  Divides the atoms among tasks so that the maximum weight gets minimal.
    !  The division can not be done in an arbitrary way:
    !  - either a task gets several atoms, but these cannot be shared with other tasks
    !    i.e. this task has only full atoms
    !  - or a atom may be divided among several tasks, but then these task can
    !    have only this, and really only this, atom
    !
    !  USES

    implicit none

    real*8, intent(in)   :: weight(n_atoms2)
    integer, intent(out) :: task_of_atom(n_atoms2), tasks_per_atom(n_atoms2)

    integer :: n_big, n_small, list_big(n_atoms2), list_small(n_atoms2), i_atom
    integer :: task_of_small_atom(n_atoms2), tasks_per_big_atom(n_atoms2), i_task
    real*8  :: aux(2,n_atoms2), weight_sum, weight_max, weight_big(n_atoms2), weight_small(n_atoms2)
    real*8  :: weight_max_big, weight_max_small


    weight_sum = sum(weight(1:n_atoms2))

    write(info_str,*) ' -----------------------------------------------------'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str,*) ' Dividing atoms among tasks for exact exchange energy:'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    do i_atom = 1, n_atoms2
      write(info_str, '(a,i5,a,g12.5)') '  Atom: ',i_atom,' Weight: ',weight(i_atom)
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    enddo
    write(info_str,*) ''
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str, '(a,g12.5)') '  Average weight per task: ',weight_sum/n_tasks
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str,*) ''
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

    ! Divide atoms into "big" (above average) and "small" (below average) atoms

    n_big = 0
    n_small = 0

    do i_atom = 1, n_atoms2
      if(weight(i_atom) > weight_sum/n_tasks) then
        n_big = n_big + 1
        list_big(n_big) = i_atom
        weight_big(n_big) = weight(i_atom)
      else
        n_small = n_small + 1
        list_small(n_small) = i_atom
        weight_small(n_small) = weight(i_atom)
      endif
    enddo

    write(info_str, '(a,i5)') '  Number of "big" (above average) atoms  : ',n_big
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str, '(a,i5)')'  Number of "small" (below average) atoms: ',n_small
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

    if(n_small > 1) then
      ! Sort the small atoms with descending weight
      do i_atom = 1, n_small
        aux(1,i_atom) = -weight_small(i_atom)
        aux(2,i_atom) = list_small(i_atom)
      enddo

      call heapsort_general(aux,2,n_small,1)
      do i_atom = 1, n_small
        weight_small(i_atom) = -aux(1,i_atom)
        list_small(i_atom) = aux(2,i_atom)
      enddo
    endif


    task_of_atom(:) = 0
    tasks_per_atom(:) = 0

    if(n_big == 0) then

      ! All atoms have less work than the average work per task.
      ! The number of tasks must be <= number of atoms, no atom is divided.

      write(info_str,*) ''
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      write(info_str, '(a)') '  All atoms have below average work per task'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      write(info_str, '(a)')'  No atom will be divided among tasks, atoms will be sorted to get best load distribution'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      if(n_tasks > n_atoms2) call aims_stop('divide_atoms INTERNAL ERROR: n_tasks > n_atoms')
      call assign_full_atoms(n_atoms2, n_tasks, weight, weight_max, task_of_atom)
      tasks_per_atom(:) = 1

    else if(n_big < n_atoms2) then

      ! Some atoms have more, some less work than the average work per task.
      ! The "big" atoms are cancidates for subdivision, the "small" atoms are not divided.
      ! The number of "big" atoms must be less than the number of tasks.

      write(info_str,*) ''
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      write(info_str, '(a)')'  There are "big" and "small" atoms'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      write(info_str, '(a,i4,a)')'  Reserving ',n_big,' tasks for "big" atoms'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      if(n_big >= n_tasks) call aims_stop('divide_atoms INTERNAL ERROR: n_big >= n_tasks')

      ! Each big atom gets at least 1 task
      tasks_per_big_atom(1:n_big) = 1.

      ! Get the max weight for the small atoms
      call assign_full_atoms(n_small, n_tasks-n_big, weight_small, weight_max_small, task_of_small_atom)
      ! Get the max weight for the big atoms (currently 1 task per big atom!)
      weight_max_big = maxval(weight_big(1:n_big))

      write(info_str, '(2(a,g12.5))')'  Max. weight for big atoms: ',weight_max_big,' small atoms: ',weight_max_small
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      weight_max = max(weight_max_big,weight_max_small)

      ! The remaining tasks are divided between big and small atoms

      do i_task = n_big+1, n_tasks-1

        write(info_str, '(a,i0,a)')'  Trying ',i_task,' tasks for "big" atoms:'
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )

        ! Check what happens if we assign i_task to the big atoms:

        ! Get the max weight for the small atoms
        call assign_full_atoms(n_small, n_tasks-i_task, weight_small, weight_max_small, task_of_small_atom)

        ! Get the max weight for the big atoms
        i_atom = maxloc(weight_big(1:n_big)/tasks_per_big_atom(1:n_big),1)
        tasks_per_big_atom(i_atom) = tasks_per_big_atom(i_atom)+1
        weight_max_big = maxval(weight_big(1:n_big)/tasks_per_big_atom(1:n_big))

        write(info_str, '(2(a,g12.5))')'  Max. weight for big atoms: ',weight_max_big,' small atoms: ',weight_max_small
        call localorb_info ( info_str,use_unit,'(A)', OL_norm )

        if(weight_max_small >= weight_max_big) then
          ! Stop redistribution, check if this step or the last one is better:
          if(max(weight_max_small,weight_max_big) >= weight_max) then
            ! Last one was better, undo task assignment from above:
            tasks_per_big_atom(i_atom) = tasks_per_big_atom(i_atom)-1
            write(info_str, '(a)')'  Last Distribution was better!'
            call localorb_info ( info_str,use_unit,'(A)', OL_norm )
          else
            write(info_str, '(a,i0,a,i0,a)')'  => Atom ',list_big(i_atom),' gets ',tasks_per_big_atom(i_atom),' tasks'
            call localorb_info ( info_str,use_unit,'(A)', OL_norm )
          endif
          exit ! we are done
        else
          ! It pays off to assign i_task to the big atoms, try next one
          weight_max = max(weight_max_small,weight_max_big)
          write(info_str, '(a,i0,a,i0,a)')'  => Atom ',list_big(i_atom),' gets ',tasks_per_big_atom(i_atom),' tasks'
          call localorb_info ( info_str,use_unit,'(A)', OL_norm )
        endif
      enddo

      ! Number of tasks for small atoms:
      i_task = n_tasks - sum(tasks_per_big_atom(1:n_big))
      if(i_task < 1) call aims_stop('divide_atoms INTERNAL ERROR: i_task < 1')

      write(info_str, '(a,i0,a,i0,a)')'  Using ',i_task,' tasks for "small" atoms ',n_tasks-i_task,', for "big" atoms'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      ! The small atoms go first and are finally assigned here
      call assign_full_atoms(n_small, i_task, weight_small, weight_max_small, task_of_small_atom)
      do i_atom = 1, n_small
        task_of_atom(list_small(i_atom)) = task_of_small_atom(i_atom)
        tasks_per_atom(list_small(i_atom)) = 1
      enddo

      ! Now the big atoms
      do i_atom = 1, n_big
        task_of_atom(list_big(i_atom)) = i_task
        tasks_per_atom(list_big(i_atom)) = tasks_per_big_atom(i_atom)
        i_task = i_task + tasks_per_big_atom(i_atom)
      enddo
      if(i_task/=n_tasks) call aims_stop('divide_atoms INTERNAL ERROR: sum tasks_per_atom')

    else

      write(info_str, '(a)')'  All atoms have above average work per task'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
      write(info_str, '(a)')'  Atoms will be divided among tasks to get best load distribution'
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      ! All atoms have more weight than the average work per task.
      ! The number of tasks must be => number of atoms, all atoms are cancidates for subdivision.
      if(n_tasks < n_atoms2) call aims_stop('divide_atoms INTERNAL ERROR: n_tasks < n_atoms')

      ! Get the tasks per each atom - start at 1 task per atom
      tasks_per_atom(1:n_atoms2) = 1

      ! Distribute the remaining tasks in a way that the maximum work is minimized
      do i_task = n_atoms2+1, n_tasks
        i_atom = maxloc(weight(:)/tasks_per_atom(:),1)
        tasks_per_atom(i_atom) = tasks_per_atom(i_atom) + 1
      enddo

      i_task = 0
      do i_atom = 1, n_atoms2
        task_of_atom(i_atom) = i_task
        i_task = i_task + tasks_per_atom(i_atom)
      enddo
      if(i_task/=n_tasks) call aims_stop('divide_atoms INTERNAL ERROR: sum tasks_per_atom')

    endif

    write(info_str, *)''
    call localorb_info ( info_str,use_unit )
    write(info_str, *)' Atom division among tasks:'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    do i_atom = 1, n_atoms2
      write(info_str, '(3(a,i5))')'  Atom ',i_atom,': (initial) task: ',task_of_atom(i_atom),&
            ', # of tasks: ',tasks_per_atom(i_atom)
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    enddo
    write(info_str, *)''
    call localorb_info ( info_str,use_unit )
    write(info_str, *)' -----------------------------------------------------'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

  end subroutine
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/assign_full_atoms
  !  NAME
  !    assign_full_atoms
  !  SYNOPSIS

  subroutine assign_full_atoms(n_atoms_full, n_tasks_full, weight, weight_max, task_of_atom)

    !  PURPOSE
    !
    !  Assigns full (i.e. undivided) atoms to tasks so that the maximum weight gets minimal.
    !  The atoms should be sorted by descending weights.
    !
    !  RJ: I hope this is the optimal way to do it
    !
    !  USES

    implicit none
    integer :: n_atoms_full, n_tasks_full, task_of_atom(n_atoms_full)
    real*8  :: weight(n_atoms_full), weight_max

    integer :: i_atom, i_task
    real*8 weight_per_task(n_tasks_full)

    weight_per_task = 0

    do i_atom = 1, n_atoms_full
      ! The task with the minimum weight gets the next atom
      i_task = minloc(weight_per_task(:),1)
      task_of_atom(i_atom) = i_task-1 ! Task numbers start at 0 (not 1!)
      weight_per_task(i_task) = weight_per_task(i_task) + weight(i_atom)
    enddo

    weight_max = maxval(weight_per_task)

  end subroutine
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/cleanup_fock_calculations
  !  NAME
  !    cleanup_fock_calculations
  !  SYNOPSIS

  subroutine cleanup_fock_matrix_calculations

    !  PURPOSE
    !
    !  Deallocates all stuff allocated in init_fock_matrix_calculations.
    !  Completely deallocates all stored coulomb matrices.
    !  Releases shared memory.
    !  Resets all other stuff
    !
    !  USES

    implicit none

    !  ARGUMENTS


    !  INPUTS
    !    None
    !  OUTPUTS
    !    None
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i, j, k, d

    if(allocated(max_norm_tmp_old))    call aims_deallocate( max_norm_tmp_old, "max_norm_tmp_old" )
    if(allocated(max_abs_dm_cols_old)) call aims_deallocate( max_abs_dm_cols_old, 'max_abs_dm_cols_old')

    if(allocated(my_atom_pair_s)) deallocate(my_atom_pair_s)
    if(allocated(my_atom_pair_e)) deallocate(my_atom_pair_e)

    if(allocated(my_atom_list)) deallocate(my_atom_list)
    if(allocated(my_basis_off)) deallocate(my_basis_off)

    if(allocated(ovlp3fn)) then
      do i = lbound(ovlp3fn,1), ubound(ovlp3fn,1)
        if(allocated(ovlp3fn(i)%m)) call aims_deallocate( ovlp3fn(i)%m, "ovlp3fn(i)%m" )
        if(allocated(ovlp3fn(i)%n)) call aims_deallocate( ovlp3fn(i)%n, "ovlp3fn(i)%n" )
      enddo
      deallocate(ovlp3fn)
    endif

    if(allocated(d_ovlp3fn)) then
      do d = 1, 3
      do i = lbound(d_ovlp3fn,1), ubound(d_ovlp3fn,1)
        if(allocated(d_ovlp3fn(i,d)%m)) call aims_deallocate( d_ovlp3fn(i,d)%m, "d_ovlp3fn(i,d)%m" )
        if(allocated(d_ovlp3fn(i,d)%n)) call aims_deallocate( d_ovlp3fn(i,d)%n, "d_ovlp3fn(i,d)%n" )
      enddo
      enddo
      deallocate(d_ovlp3fn)
    endif

    if(allocated(AS_d_ovlp3fn)) then
      do d = 1, AS_components, 1
      do i = lbound(AS_d_ovlp3fn,1), ubound(AS_d_ovlp3fn,1)
        if(allocated(AS_d_ovlp3fn(i,d)%m)) call aims_deallocate( AS_d_ovlp3fn(i,d)%m, "AS_d_ovlp3fn(i,d)%m" )
        if(allocated(AS_d_ovlp3fn(i,d)%n)) call aims_deallocate( AS_d_ovlp3fn(i,d)%n, "AS_d_ovlp3fn(i,d)%n" )
      enddo
      enddo
      deallocate(AS_d_ovlp3fn)
    endif

    if(allocated(ovlp3fn_SR)) then
      do i = lbound(ovlp3fn_SR,1), ubound(ovlp3fn,1)
        if(allocated(ovlp3fn_SR(i)%m)) call aims_deallocate( ovlp3fn_SR(i)%m, "ovlp3fn_SR(i)%m" )
        if(allocated(ovlp3fn_SR(i)%n)) call aims_deallocate( ovlp3fn_SR(i)%n, "ovlp3fn_SR(i)%n" )
      enddo
      deallocate(ovlp3fn_SR)
    endif

    if(allocated(d_ovlp3fn_SR)) then
      do d = 1, 3
      do i = lbound(d_ovlp3fn_SR,1), ubound(d_ovlp3fn_SR,1)
        if(allocated(d_ovlp3fn_SR(i,d)%m)) call aims_deallocate( d_ovlp3fn_SR(i,d)%m, "d_ovlp3fn_SR(i,d)%m" )
        if(allocated(d_ovlp3fn_SR(i,d)%n)) call aims_deallocate( d_ovlp3fn_SR(i,d)%n, "d_ovlp3fn_SR(i,d)%n" )
      enddo
      enddo
      deallocate(d_ovlp3fn_SR)
    endif

    if(allocated(AS_d_ovlp3fn_SR)) then
      do d = 1, AS_components, 1
      do i = lbound(AS_d_ovlp3fn_SR,1), ubound(AS_d_ovlp3fn_SR,1)
        if(allocated(AS_d_ovlp3fn_SR(i,d)%m)) call aims_deallocate( AS_d_ovlp3fn_SR(i,d)%m, "AS_d_ovlp3fn_SR(i,d)%m" )
        if(allocated(AS_d_ovlp3fn_SR(i,d)%n)) call aims_deallocate( AS_d_ovlp3fn_SR(i,d)%n, "AS_d_ovlp3fn_SR(i,d)%n" )
      enddo
      enddo
      deallocate(AS_d_ovlp3fn_SR)
    endif

    if(allocated(max_norm_ovlp3fn)) deallocate(max_norm_ovlp3fn)
    if(allocated(max_norm_d_ovlp3fn)) deallocate(max_norm_d_ovlp3fn)
    if(allocated(AS_max_norm_d_ovlp3fn)) deallocate(AS_max_norm_d_ovlp3fn)
    if(allocated(max_norm_ovlp3fn_per_latom)) deallocate(max_norm_ovlp3fn_per_latom)
    if(allocated(max_norm_d_ovlp3fn_per_latom)) deallocate(max_norm_d_ovlp3fn_per_latom)
    if(allocated(AS_max_norm_d_ovlp3fn_per_latom)) deallocate(AS_max_norm_d_ovlp3fn_per_latom)

    if (.not.(use_mpi_in_place)) then
       if(allocated(max_norm_ovlp3fn_rcv)) deallocate(max_norm_ovlp3fn_rcv)
       if(allocated(max_norm_d_ovlp3fn_rcv)) deallocate(max_norm_d_ovlp3fn_rcv)
       if(allocated(AS_max_norm_d_ovlp3fn_rcv)) deallocate(AS_max_norm_d_ovlp3fn_rcv)
    endif

    if(allocated(max_norm_ovlp3fn_SR)) deallocate(max_norm_ovlp3fn_SR)
    if(allocated(max_norm_d_ovlp3fn_SR)) deallocate(max_norm_d_ovlp3fn_SR)
    if(allocated(AS_max_norm_d_ovlp3fn_SR)) deallocate(AS_max_norm_d_ovlp3fn_SR)
    if(allocated(max_norm_ovlp3fn_per_latom_SR)) deallocate(max_norm_ovlp3fn_per_latom_SR)
    if(allocated(max_norm_d_ovlp3fn_per_latom_SR)) deallocate(max_norm_d_ovlp3fn_per_latom_SR)
    if(allocated(AS_max_norm_d_ovlp3fn_per_latom_SR)) deallocate(AS_max_norm_d_ovlp3fn_per_latom_SR)

    if (.not.(use_mpi_in_place)) then
       if(allocated(max_norm_ovlp3fn_rcv_SR)) deallocate(max_norm_ovlp3fn_rcv_SR)
       if(allocated(max_norm_d_ovlp3fn_rcv_SR)) deallocate(max_norm_d_ovlp3fn_rcv_SR)
       if(allocated(AS_max_norm_d_ovlp3fn_rcv_SR)) deallocate(AS_max_norm_d_ovlp3fn_rcv_SR)
    endif

    if(allocated(bvk_cell_idx)) deallocate(bvk_cell_idx)
    if(allocated(inv_cell_bvk)) deallocate(inv_cell_bvk)
    if(allocated(add_cells)) deallocate(add_cells)
    if(allocated(sub_cells)) deallocate(sub_cells)

    if(allocated(pair_offset)) deallocate(pair_offset)
    if(allocated(pair_flag_bvk)) deallocate(pair_flag_bvk)
    if(allocated(my_pair_flag_bvk)) deallocate(my_pair_flag_bvk)
    if(allocated(pair_list)) deallocate(pair_list)

    if(allocated(coul_mat_norm)) deallocate(coul_mat_norm)
    if(allocated(d_coul_mat_norm)) deallocate(d_coul_mat_norm)
    if(allocated(AS_d_coul_mat_norm)) deallocate(AS_d_coul_mat_norm)

    if(allocated(coul_mat_store)) then
      do k = lbound(coul_mat_store,3), ubound(coul_mat_store,3)
      do j = lbound(coul_mat_store,2), ubound(coul_mat_store,2)
      do i = lbound(coul_mat_store,1), ubound(coul_mat_store,1)
        if(allocated(coul_mat_store(i,j,k)%row_idx)) &
             call aims_deallocate( coul_mat_store(i,j,k)%row_idx, "coul_mat_store%row_idx" )
        if(allocated(coul_mat_store(i,j,k)%col_idx)) &
             call aims_deallocate( coul_mat_store(i,j,k)%col_idx, "coul_mat_store%col_idx" )
        if(allocated(coul_mat_store(i,j,k)%mat    )) &
             call aims_deallocate( coul_mat_store(i,j,k)%mat,         "coul_mat_store%mat" )
      enddo
      enddo
      enddo
      deallocate(coul_mat_store)
    endif

    if(allocated(d_coul_mat_store)) then
      do d = 1, 3
      do k = lbound(d_coul_mat_store,3), ubound(d_coul_mat_store,3)
      do j = lbound(d_coul_mat_store,2), ubound(d_coul_mat_store,2)
      do i = lbound(d_coul_mat_store,1), ubound(d_coul_mat_store,1)
        if(allocated(d_coul_mat_store(i,j,k,d)%row_idx)) &
             call aims_deallocate( d_coul_mat_store(i,j,k,d)%row_idx, "d_coul_mat_store%row_idx" )
        if(allocated(d_coul_mat_store(i,j,k,d)%col_idx)) &
             call aims_deallocate( d_coul_mat_store(i,j,k,d)%col_idx, "d_coul_mat_store%col_idx" )
        if(allocated(d_coul_mat_store(i,j,k,d)%mat    )) &
             call aims_deallocate( d_coul_mat_store(i,j,k,d)%mat,         "d_coul_mat_store%mat" )
      enddo
      enddo
      enddo
      enddo
      deallocate(d_coul_mat_store)
    endif

    if(allocated(AS_d_coul_mat_store)) then
      do d = 1, AS_components, 1
      do k = lbound(AS_d_coul_mat_store,3), ubound(AS_d_coul_mat_store,3)
      do j = lbound(AS_d_coul_mat_store,2), ubound(AS_d_coul_mat_store,2)
      do i = lbound(AS_d_coul_mat_store,1), ubound(AS_d_coul_mat_store,1)
        if(allocated(AS_d_coul_mat_store(i,j,k,d)%row_idx)) &
             call aims_deallocate( AS_d_coul_mat_store(i,j,k,d)%row_idx, "AS_d_coul_mat_store%row_idx" )
        if(allocated(AS_d_coul_mat_store(i,j,k,d)%col_idx)) &
             call aims_deallocate( AS_d_coul_mat_store(i,j,k,d)%col_idx, "AS_d_coul_mat_store%col_idx" )
        if(allocated(AS_d_coul_mat_store(i,j,k,d)%mat    )) &
             call aims_deallocate( AS_d_coul_mat_store(i,j,k,d)%mat,         "AS_d_coul_mat_store%mat"  )
      enddo
      enddo
      enddo
      enddo
      deallocate(AS_d_coul_mat_store)
    endif

    if(allocated(coul_mat_norm_SR)) deallocate(coul_mat_norm_SR)
    if(allocated(d_coul_mat_norm_SR)) deallocate(d_coul_mat_norm_SR)
    if(allocated(AS_d_coul_mat_norm_SR)) deallocate(AS_d_coul_mat_norm_SR)

    if(allocated(coul_mat_store_SR)) then
      do k = lbound(coul_mat_store_SR,3), ubound(coul_mat_store_SR,3)
      do j = lbound(coul_mat_store_SR,2), ubound(coul_mat_store_SR,2)
      do i = lbound(coul_mat_store_SR,1), ubound(coul_mat_store_SR,1)
        if(allocated(coul_mat_store_SR(i,j,k)%row_idx)) deallocate(coul_mat_store_SR(i,j,k)%row_idx)
        if(allocated(coul_mat_store_SR(i,j,k)%col_idx)) deallocate(coul_mat_store_SR(i,j,k)%col_idx)
        if(allocated(coul_mat_store_SR(i,j,k)%mat    )) deallocate(coul_mat_store_SR(i,j,k)%mat    )
      enddo
      enddo
      enddo
      deallocate(coul_mat_store_SR)
    endif

    if(allocated(d_coul_mat_store_SR)) then
      do d = 1, 3
      do k = lbound(d_coul_mat_store_SR,3), ubound(d_coul_mat_store_SR,3)
      do j = lbound(d_coul_mat_store_SR,2), ubound(d_coul_mat_store_SR,2)
      do i = lbound(d_coul_mat_store_SR,1), ubound(d_coul_mat_store_SR,1)
        if(allocated(d_coul_mat_store_SR(i,j,k,d)%row_idx)) deallocate(d_coul_mat_store_SR(i,j,k,d)%row_idx)
        if(allocated(d_coul_mat_store_SR(i,j,k,d)%col_idx)) deallocate(d_coul_mat_store_SR(i,j,k,d)%col_idx)
        if(allocated(d_coul_mat_store_SR(i,j,k,d)%mat    )) deallocate(d_coul_mat_store_SR(i,j,k,d)%mat    )
      enddo
      enddo
      enddo
      enddo
      deallocate(d_coul_mat_store_SR)
    endif

    if(allocated(AS_d_coul_mat_store_SR)) then
      do d = 1, AS_components, 1
      do k = lbound(AS_d_coul_mat_store_SR,3), ubound(AS_d_coul_mat_store_SR,3)
      do j = lbound(AS_d_coul_mat_store_SR,2), ubound(AS_d_coul_mat_store_SR,2)
      do i = lbound(AS_d_coul_mat_store_SR,1), ubound(AS_d_coul_mat_store_SR,1)
        if(allocated(AS_d_coul_mat_store_SR(i,j,k,d)%row_idx)) deallocate(AS_d_coul_mat_store_SR(i,j,k,d)%row_idx)
        if(allocated(AS_d_coul_mat_store_SR(i,j,k,d)%col_idx)) deallocate(AS_d_coul_mat_store_SR(i,j,k,d)%col_idx)
        if(allocated(AS_d_coul_mat_store_SR(i,j,k,d)%mat    )) deallocate(AS_d_coul_mat_store_SR(i,j,k,d)%mat    )
      enddo
      enddo
      enddo
      enddo
      deallocate(AS_d_coul_mat_store_SR)
    endif

    init_calc_deriv = .false.
    AS_init_stress  = .false.

  end subroutine cleanup_fock_matrix_calculations
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/calculate_realspace_coulomb
  !  NAME
  !    calculate_realspace_coulomb
  !  SYNOPSIS

  subroutine calculate_realspace_coulomb(i_species_1, i_species_2, Dvec, &
                                         auxmat, d_auxmat, AS_d_auxmat, nbb1, nbb2_s, nbb2_e)

    !  PURPOSE
    !
    !  Does the same as sum_up_auxmat_qvecs_from_realspace() but only for 1 atom/atom pair
    !  and immediatly transforms the resulting matrices back to realspace.
    !  The combination of calculation and immediate backtransformation can be made much
    !  more efficiently than doing both steps separate.
    !
    !  Description of the original sum_up_auxmat_qvecs_from_realspace:
    !
    !    Calculate auxiliary matrix for a given set of q-points by a
    !    real-space sum.  This will only work if the radial interaction
    !    splines of this modules are "short-ranged", i.e. if either
    !    (.not. have_mp_far) [bare overlap, HSE, cutCb] or if prepared for an
    !    Ewald-like treatment by initialize_periodic_tb_auxmat().
    !
    !  USES

    use bravais, only: get_n_supercells
    use geometry, only: recip_lattice_vector, lattice_vector
    use prodbas, only: sp2n_basbas_sp
    use tight_binding_auxmat, only: have_mp_far, iarange_species, &
        fast_calculate_tb_auxmat
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2, nbb1, nbb2_s, nbb2_e
    real*8, intent(IN)  :: Dvec(3)
    real*8, intent(OUT) :: auxmat(nbb1, nbb2_s:nbb2_e, my_cm_cell_num)
    real*8, intent(OUT) :: d_auxmat(nbb1, nbb2_s:nbb2_e, 3, my_cm_cell_num)
    real*8, intent(OUT) :: AS_d_auxmat(nbb1, nbb2_s:nbb2_e, AS_components, my_cm_cell_num)


    !  INPUTS
    !    o i_species_1, i_species_2 -- Species of the two atoms involved.
    !    o Dvec -- Essentially coords(:, i_atom_2) - coords(:, i_atom_1).
    !    o nbb1 -- first dimension of auxmat/d_auxmat
    !    o nbb2_s, nbb2_e - limits of second dimension of auxmat/d_auxmat
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !    o d_auxmat -- derivative of auxmat with respect to Dvec,
    !                  calculated only if init_calc_deriv is set.
    !    o AS_d_auxmat -- derivative of auxmat with respect to strain,
    !                     calculated only if AS_init_stress is .true.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: n_spbb_1, n_spbb_2
    integer :: n_supercells(3), a1, a2, a3, i_k_point, i_cell, n_cnt, i_cell_1, &
               info
    real*8 :: Rvec(3), Cvec(3), phi, Rmax, Cmax, qvec(3)
    real*8, allocatable :: aux(:,:), d_aux(:,:,:)
    complex*16 :: c_phases(n_k_points), c_phases_fft(n_cells_bvk)

    character(*), parameter :: func = 'calculate_realspace_coulomb'

    if (have_mp_far) then
       ! For the bare Coulomb potential, the radial interaction splines
       ! have to be fudged according to the Ewald procedure and have_mp_far
       ! has to be set to .false., afterwards.
       call aims_stop('Cannot use realspace sum with multipoles', func)
    end if

    n_spbb_1 = sp2n_basbas_sp(i_species_1)
    n_spbb_2 = sp2n_basbas_sp(i_species_2)

    allocate(aux(n_spbb_1,n_spbb_2), stat=info)
    call check_allocation(info, 'aux', func)

    if(init_calc_deriv) then
       allocate(d_aux(n_spbb_1,n_spbb_2,3), stat=info)
    else
       allocate(d_aux(1,1,1), stat=info)
    endif
    call check_allocation(info, 'd_aux', func)

    auxmat(:,:,:) = 0.
    if(init_calc_deriv) d_auxmat(:,:,:,:) = 0.
    if(AS_init_stress) AS_d_auxmat(:,:,:,:) = 0.0d0


    Rmax = iarange_species(i_species_1, i_species_2)
    Cmax = Rmax + sqrt(sum(Dvec**2))
    n_supercells = 0   ! Important for n_periodic < 3.
    call get_n_supercells(n_periodic, lattice_vector, Cmax, n_supercells)
    do a1 = -n_supercells(1), n_supercells(1)
       do a2 = -n_supercells(2), n_supercells(2)
          do a3 = -n_supercells(3), n_supercells(3)
             Cvec = matmul(lattice_vector, (/a1, a2, a3/))
             Rvec = Dvec + Cvec
             if (sum(Rvec**2) > Rmax**2) cycle

             call fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec, aux, d_aux, init_calc_deriv)

            do i_k_point = 1, n_k_points
              qvec(:) = matmul(recip_lattice_vector, k_point_list(i_k_point,:))
              phi = dot_product(qvec(:), Cvec)
              c_phases(i_k_point) = cmplx(cos(phi), sin(phi), kind(0.d0))
            end do

            c_phases_fft(:) = 0
            do i_k_point = 1, n_k_points
                do i_cell = 1, n_cells_bvk
                  c_phases_fft(i_cell) = c_phases_fft(i_cell) &
                      + c_phases(i_k_point) * k_weights(i_k_point) &
                      * conjg(k_phase_exx(i_cell,i_k_point))
                enddo
             enddo

             n_cnt = 0
             do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc
               n_cnt = n_cnt + 1
               auxmat(1:n_spbb_1,nbb2_s:nbb2_e,n_cnt) = &
                 auxmat(1:n_spbb_1,nbb2_s:nbb2_e,n_cnt) + aux(1:n_spbb_1,nbb2_s:nbb2_e)*dble(c_phases_fft(i_cell))
               if(init_calc_deriv) &
                 d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,1:3,n_cnt) = &
                   d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,1:3,n_cnt) + d_aux(1:n_spbb_1,nbb2_s:nbb2_e,1:3)*dble(c_phases_fft(i_cell))

               !FK: Multiply nuclear distances and derivatives with respect to nuclear coordinates
               if (AS_init_stress) then
                 do AS_index = 1, AS_components, 1
                   AS_d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,AS_index,n_cnt) = &
                     AS_d_auxmat(1:n_spbb_1,nbb2_s:nbb2_e,AS_index,n_cnt) + &
                     d_aux(1:n_spbb_1,nbb2_s:nbb2_e,AS_l_index(AS_index)) * &
                     dble(c_phases_fft(i_cell)) * Rvec(AS_m_index(AS_index))
                 end do
               end if
             enddo
          end do
       end do
    end do

  end subroutine calculate_realspace_coulomb
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/evaluate_exchange_matr_realspace_p0
  !  NAME
  !    evaluate_exchange_matr_realspace_p0
  !  SYNOPSIS

  subroutine evaluate_exchange_matr_realspace_p0(KS_eigenvector,KS_eigenvector_complex,occ_numbers, &
                                                 exx_ene, d_exx_ene, calc_deriv, AS_stress_on, density_matrix_real_in, density_matrix_complex_in)

    !  PURPOSE
    !     This subroutine calculates Fock matrix

    !
    !  USES

    use analytical_stress, only: as_exx_stress_local, as_sync_exx_stress
    use basis, only: max_n_basis_sp, sp2n_basis_sp, atom2basis_off, basis_atom, max_n_basis_sp2
    use constants, only: HARTREE_OVER_BOHR
    use geometry, only: species
    use hartree_fock_p0, only: hf_exchange_matr_real, hf_exchange_matr_complex,&
        hf_exchange_matr_real_SR, hf_exchange_matr_complex_SR
    use prodbas, only: max_n_basbas_sp, atom2basbas_off, sp2n_basbas_sp
    use restart_elsi, only: elsi_restart_scalapack
    use scalapack_wrapper, only: dm_scalapack_real => ham, &
        dm_scalapack_complex => ham_complex, l_col, l_row, my_k_point, &
        construct_dm_scalapack, set_full_matrix_real, set_full_matrix_complex
    use species_data, only: no_basis, species_pseudoized
    use sym_base, only: evaluate_densmat_hf_sym, FT_densmat_sym, &
        FT_hf_exchange_matr_sym
    use synchronize_mpi, only: sync_vector, sync_integer_vector, &
        sync_real_number
    use timing_core, only: get_times, get_timestamps
    use timing, only: tot_time_matrix_io, tot_clock_time_matrix_io

    implicit none

    !  ARGUMENTS
    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8  :: exx_ene, d_exx_ene(3,n_atoms)
    logical, intent(in) :: calc_deriv
    logical, intent(in) :: AS_stress_on
    ! If present, these will be used instead of constructing the
    ! density matrix below.
    real*8, intent(in), optional :: density_matrix_real_in(:,:,:)
    complex*16, intent(in), optional :: density_matrix_complex_in(:,:,:)

    !  OUTPUTS
    !  o exx_ene - Exact-exchange enery
    !  o d_exx_ene - Derivative of exx_ene
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
	 real*8  :: exx_ene_SR, d_exx_ene_SR(3,n_atoms2)

    real*8, allocatable :: fock_matrix(:,:,:,:)
    real*8, allocatable :: fock_matrix_SR(:,:,:,:)
    real*8, allocatable :: fock_matrix_row(:,:,:)
    real*8, allocatable :: fock_tmp(:,:)
    real*8, allocatable :: fock_tmp_SR(:,:)
    real*8, allocatable :: dm_cols(:,:,:,:)
    real*8, allocatable :: dm_row(:,:,:)
    real*8, allocatable :: dm_tmp(:,:)
    real*8, allocatable :: max_abs_dm_row(:,:)
    real*8, allocatable :: max_abs_dm_cols(:,:,:)
    real*8, allocatable :: ovlp3fn_rcv(:,:,:)
    real*8, allocatable :: norm_ovlp3fn_rcv(:,:)
    real*8, allocatable :: max_norm_ovlp3fn_rcv(:,:)
    real*8, allocatable :: dm_aux(:,:), dm_x_o3fn_aux(:,:)
    real*8, allocatable :: tmp_prod(:,:), tmp1(:,:,:,:), tmp2(:,:,:,:), tmpx(:,:,:,:)
    real*8, allocatable :: max_norm_tmp(:), max_norm_tmp2(:), max_fact_tmp1(:,:), AS_max_fact_tmp1(:)
    real*8, allocatable :: fock_matrix_tmp(:,:), aux(:,:), tmp_atom(:,:)
    real*8, allocatable :: ovl_fact11(:,:), ovl_fact11_rcv(:,:)
    real*8, allocatable :: ovl_fact12(:,:), ovl_fact12_rcv(:,:)
    real*8, allocatable :: ovl_fact22(:), ovl_fact22_rcv(:)
    real*8, allocatable :: ovl_fact22_SR(:), ovl_fact22_rcv_SR(:)
    real*8, allocatable :: ovl_fact11_d(:,:), ovl_fact11_d_rcv(:,:)
    real*8, allocatable :: ovl_fact12_d(:,:), ovl_fact12_d_rcv(:,:)
    real*8, allocatable :: ovl_fact22_d(:), ovl_fact22_d_rcv(:)
    real*8, allocatable :: ovl_fact22_d_SR(:), ovl_fact22_d_rcv_SR(:)
    real*8, allocatable :: AS_ovl_fact11_d(:,:), AS_ovl_fact11_d_rcv(:,:)
    real*8, allocatable :: AS_ovl_fact12_d(:,:), AS_ovl_fact12_d_rcv(:,:)
    real*8, allocatable :: AS_ovl_fact22_d(:), AS_ovl_fact22_d_rcv(:)
    real*8, allocatable :: AS_ovl_fact22_d_SR(:), AS_ovl_fact22_d_rcv_SR(:)
    real*8 :: x1, x2, s, max_fact1, max_fact2, AS_max_fact3, max_abs_dm, max_abs_dm_p1
    real*8 :: max_fact2_rcv, AS_max_fact3_rcv
    real*8 :: max_norm_dxo, fact, max_abs_dm_row_tot
    real*8 :: crit_val_deriv
    real*8 :: AS_crit_val_deriv
    real*8 :: cm_norm
    real*8 :: AS_cm_norm

    complex*16 :: cfact

    integer :: mpierr
    integer :: i_spin, i_grad
    integer :: my_i_basis
    integer :: i_cell_bvk
    integer :: i, j, i_k, i_k_point, n, i_atom_pair, i_row, n_off
    integer :: i_cell, i_cell_1, i_cell_2, i_cell_dm, i_cell_fock, i_cell_cm1, i_cell_cm2
    integer :: i_cell_p1, i_cell_m1
    integer :: i_basis, i_basis_l1, i_basis_r1, i_basis_l2, i_basis_r2
    integer :: i_pair_1, i_pair_2
    integer :: i_basis_s, i_basis_e, my_s, my_e
    integer :: i_atom_l1, i_atom_r1, i_my_atom
    integer :: nbb, work_count
    integer :: i_atom, na_off, na_len
    integer :: i_atom_l2, i_atom_r2, i_atom_r1_prev
    integer :: nl2_off, nl2_len, nl2_bb_off, nl2_bb_len
    integer :: nr2_off, nr2_len
    integer :: nl1_off, nl1_len
    integer :: n_ovlp3fn_rcv, n_ovlp3fn_rcv_1, n_ovlp3fn_rcv_2, no3fn_off, no3fn_len, n_bb_s, n_bb_e
    integer :: n_cells_used
    integer :: max_ovlp3fn_rcv
    integer :: max_pair_1, n_pair_1
    integer :: nr1_off, nr1_len, ir1, il1, ir2, il2
    integer :: idx_cells(n_cells_bvk)
    integer :: n_cols1, n_cols2, n_cols, AS_n_cols
    integer, allocatable :: pair_1_list(:,:)
    integer, allocatable :: idx_ovlp3fn_rcv(:,:)

    logical :: mult_cm_r1(n_cells_bvk), mult_d_cm_r1(n_cells_bvk), AS_mult_d_cm_r1(n_cells_bvk)
    logical :: cell_used(n_cells_bvk, n_spin), cell_tmpx_used(n_cells_bvk, n_spin)
    logical :: cell_used_rcv(n_cells_bvk, n_spin)

!DB
!BL: n_real_atoms is a global variable in dimensions.f90
!    renamed it to n_real_atoms_DB to overcome compilation error
!    to DB: Is this ok?
    integer :: n_real_atoms_DB, n_atoms_save

    integer :: info
    ! Needed for circumventing PGI bug
    integer*8 :: mem_tmpx

    ! Statement functions for offset of basis/basbas (only for better readability!!!)
    integer :: atom2basis_len, atom2basbas_len
    atom2basis_len(i_atom)  = sp2n_basis_sp(species(i_atom))
    atom2basbas_len(i_atom) = sp2n_basbas_sp(species(vb_atom2atom(i_atom)))

    real*8 :: t_read
    real*8 :: tmp

! WPH:  To discover why the "Other" times are so long on Theta, I've added
!       comments to the code to indicate when a timing block is beginning and
!       ending.  The subdivisions beyond what is already in the code (for
!       example, the #3 in "Other #3") are strictly my own.
!
!       The blocks are:
!       - times(1):  DM_x_o3fn  - matrix multiply, density matrix times ovlp3fn
!       - times(2):  CM_x_o3fn  - matrix multiply, coulomb matrix times ovlp3fn
!                                 That is, the collective timing for the body of
!                                 mult_coul_mat_left, which is called exclusively
!                                 within the "Products" block.  This timing is
!                                 subtracted off from times(3) at the end, so it
!                                 is an independent timing in the final output.
!       - times(3):  Products   - Roughly, all DGEMM and DGEMV calls that do not
!                                 belong to the previous two categories.
!       - times(4):  Sync       - The time spent in code blocks buttressed by
!                                 time_sync_start and time_sync_end subroutine
!                                 calls.  As you might have guessed from the
!                                 name, these code blocks are dominated by
!                                 MPI_BCAST/MPI_ALLREDUCE calls.  This timing does
!                                 not include the time spent in the MPI_BARRIER
!                                 subroutine call in time_sync_start (see next).
!       - times(5):  Imbalance  - The time spent at MPI_BARRIER throughout the
!                                 code.  Includes the MPI_BARRIERs in
!                                 time_sync_start as well as one explicit call
!                                 within the body of the main subroutine.
!       - times(6):  Pre/Post   - What it says on the tin, pre- and post-processing.
!                                 Contains (almost) everything outside the main
!                                 do loop.
!       - times(19): Other      - The placeholder for everything that did not occur
!                                 in one of the preceeding timing blocks.
!       - times(20): Total calc - The total time spent calculating the Fock matrix.
!
!       ttt0 is the temporary variable for timing the current block, and ttts is the
!       temporary variable for timing the total Fock matrix calculation.  tttx
!       is not used.

! WPH:  Variables to time the various "Other" blocks in the code
real*8 ttt_other, times_other(100)

! WPH:  As an aside, the following comment is from 2012.
! The following is for timing etc - to be removed later
real*8 ttt0, ttts, tttx, times(100), times_rcv(100)
times(:) = 0
times_rcv(:) = 0
times_other(:) =  0
call mpi_barrier(mpi_comm_global,mpierr)
ttts = mpi_wtime()
ttt_other = ttts
! Begin timing,       Total
! Begin timing block, Other #1
    if(calc_deriv .and. .not. init_calc_deriv) &
      call aims_stop('Requested derivatives of exact exchange energy are not initialized')

    if(AS_stress_on .and. .not. AS_init_stress) &
      call aims_stop('Requested strain derivatives of exact exchange energy are not initialized')

    ! The critical value for screening gradients is kept in an extra variable
    ! since it is possible that the gradients are needed to a higher accuracy.
    ! For now, we just set it to crit_val

    crit_val_deriv = crit_val
    AS_crit_val_deriv = crit_val

    !---------------------------------------------------------------------------

    ! begin work

    write(info_str,'(2X,A)') 'Calculating non-local Hartree-Fock exchange by two-center RI (RI-LVL).'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm  )
    write(info_str,'(A,F15.12)') 'screening_threshold (crit_val) = ', crit_val
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )


! DB 10/02/13
    ! For QM/MM embedding one mights to put additional integration grids around
    ! embedding monopoles (emty sites without any basis functions).
    ! This routine has a problem when dealing with such empty atoms, e.g. dgemm
    ! is annoyed.
    ! The workaround here is that we change the number of atoms to those which are
    ! real atoms with basis functions (and change it back at the end of the routine).
    ! This is not the nicest way to do it, but at least we maintain this highly
    ! optimzed structure of this routine.
    ! The only drawback: this relies on correct ordering of atoms
    ! (1:n_real_atoms_DB,1:n_empty_atoms,1:n_pp_atoms)

    n_real_atoms_DB = 0
    do i_atom = 1, n_atoms
       ! HJ 28/09/2018
       ! Changed this if statement to be consistent with the one 7 lines down.
       if (species_pseudoized(species(i_atom)).or.no_basis(species(i_atom))) cycle
       n_real_atoms_DB = n_real_atoms_DB + 1
    enddo
    n_atoms_save = n_atoms
    n_atoms = n_real_atoms_DB
    ! check for correct ordering
    do i_atom = 1, n_atoms
      if(species_pseudoized(species(i_atom)).or.no_basis(species(i_atom))) then
         call aims_stop('Internal: inconsistent ordering of atoms')
      endif
    enddo


    !---------------------------------------------------------------------------
    ! Get max number of columns in ovlp3fn_rcv, this is needed for allocations below
! End timing block,   Other #1
! Begin timing block, Pre/Post #1
ttt0 = mpi_wtime()
times_other(1) = times_other(1) + ttt0 - ttt_other
    max_ovlp3fn_rcv = 0

    do i_atom_r1 = 1, n_atoms2
      n_ovlp3fn_rcv = 0
      do i_cell_1 = 1, n_cells_bvk
        do i_atom_l1 = 1, n_atoms2
          nl1_len = atom2basis_len2(i_atom_l1)
          if(pair_flag_bvk(i_atom_l1,i_atom_r1,i_cell_1)) n_ovlp3fn_rcv = n_ovlp3fn_rcv + nl1_len
        enddo
      enddo
      max_ovlp3fn_rcv = max(max_ovlp3fn_rcv,n_ovlp3fn_rcv)
    enddo

    ! Max number of atom pairs for all values of i_basis_r1

    max_pair_1 = 0
    do i_atom_r1 = 1, n_atoms2
      max_pair_1 = max(max_pair_1,count(pair_flag_bvk(:,i_atom_r1,:)))
    enddo

    !---------------------------------------------------------------------------

    ! Allocations

    call aims_allocate(dm_aux,max_n_basis_sp2,n_cells_bvk,'dm_aux')
    call aims_allocate(dm_x_o3fn_aux,max_n_basbas_sp,n_cells_bvk,'dm_x_o3fn_aux')

    call aims_allocate(dm_cols,n_basis,my_n_basis,n_cells_bvk,n_spin,'+dm_cols')
    call aims_allocate(dm_row,n_basis,n_cells_bvk,n_spin,'dm_row')

    call aims_allocate(ovlp3fn_rcv,max_n_basbas_sp,max_ovlp3fn_rcv,2,'+ovlp3fn_rcv')
    call aims_allocate(norm_ovlp3fn_rcv,max_ovlp3fn_rcv,2,'norm_ovlp3fn_rcv')
    call aims_allocate(idx_ovlp3fn_rcv,max_ovlp3fn_rcv,2,'idx_ovlp3fn_rcv')

    call aims_allocate(max_abs_dm_row,n_atoms2,n_cells_bvk,'max_abs_dm_row')
    call aims_allocate(max_abs_dm_cols,n_atoms2,n_atoms2,n_cells_bvk,'max_abs_dm_cols')
    call aims_allocate(max_fact_tmp1,n_cells_bvk,2,'max_fact_tmp1')
    call aims_allocate(AS_max_fact_tmp1,n_cells_bvk,'AS_max_fact_tmp1')

    call aims_allocate(max_norm_ovlp3fn_rcv,max_pair_1,2,'max_norm_ovlp3fn_rcv')

    call aims_allocate(fock_matrix,n_basis,my_n_basis,n_cells_bvk,n_spin,'+fock_matrix')
    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
    	call aims_allocate(fock_matrix_SR,n_basis,my_n_basis,n_cells_bvk,n_spin,'+fock_matrix_SR')
    end if
    call aims_allocate(fock_matrix_row,n_basis,n_cells_bvk,n_spin,'+fock_matrix_row')

    call aims_allocate(tmp_prod,max_n_basbas_sp,max_n_basis_sp2,'tmp_prod')
    call aims_allocate(tmp_atom,max_n_basis_sp2,max_n_basis_sp2,'tmp_atom')

    call aims_allocate(max_norm_tmp,n_cells_bvk,'max_norm_tmp')
    call aims_allocate(max_norm_tmp2,n_cells_bvk,'max_norm_tmp2')
    call aims_allocate(fock_matrix_tmp,max_n_basis_sp2,n_cells_bvk,'fock_matrix_tmp')
    call aims_allocate(aux,2,max(n_cells_bvk,max_n_basis_sp2),'aux')

    call aims_allocate(pair_1_list,max_pair_1,6,'pair_1_list')

    call aims_allocate(ovl_fact22, n_atoms2, 'ovl_fact22')
    call aims_allocate(ovl_fact12, n_atoms2, n_cells_bvk, 'ovl_fact12')
    call aims_allocate(ovl_fact11, n_atoms2, n_cells_bvk, 'ovl_fact11')

    call aims_allocate(ovl_fact22_d, n_atoms2, 'ovl_fact22_d')
    call aims_allocate(ovl_fact12_d, n_atoms2, n_cells_bvk, 'ovl_fact12_')
    call aims_allocate(ovl_fact11_d, n_atoms2, n_cells_bvk, 'ovl_fact11_d')

    call aims_allocate(AS_ovl_fact22_d, n_atoms2, 'AS_ovl_fact22_d')
    call aims_allocate(AS_ovl_fact12_d, n_atoms2, n_cells_bvk, 'AS_ovl_fact12_d')
    call aims_allocate(AS_ovl_fact11_d, n_atoms2, n_cells_bvk, 'AS_ovl_fact11_d')

	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
    	call aims_allocate(ovl_fact22_SR, n_atoms2, 'ovl_fact22_SR')
    	call aims_allocate(ovl_fact22_d_SR, n_atoms2, 'ovl_fact22_d_SR')
    	call aims_allocate(AS_ovl_fact22_d_SR, n_atoms2, 'AS_ovl_fact22_d_SR')
    end if

    if (.not.(use_mpi_in_place)) then
       call aims_allocate(ovl_fact22_rcv, n_atoms2, 'ovl_fact22_rcv')
       call aims_allocate(ovl_fact12_rcv, n_atoms2, n_cells_bvk, 'ovl_fact12_rcv')
       call aims_allocate(ovl_fact11_rcv, n_atoms2, n_cells_bvk, 'ovl_fact11_rcv')
       call aims_allocate(ovl_fact22_d_rcv, n_atoms2, 'ovl_fact22_d_rcv')
       call aims_allocate(ovl_fact12_d_rcv, n_atoms2, n_cells_bvk, 'ovl_fact12_d_rcv')
       call aims_allocate(ovl_fact11_d_rcv, n_atoms2, n_cells_bvk, 'ovl_fact11_d_rcv')
       call aims_allocate(AS_ovl_fact22_d_rcv, n_atoms2, 'AS_ovl_fact22_d_rcv')
       call aims_allocate(AS_ovl_fact12_d_rcv, n_atoms2, n_cells_bvk, 'AS_ovl_fact12_d_rcv')
       call aims_allocate(AS_ovl_fact11_d_rcv, n_atoms2, n_cells_bvk, 'AS_ovl_fact11_d_rcv')

       if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
       		call aims_allocate(ovl_fact22_rcv_SR, n_atoms2, 'ovl_fact22_rcv_SR')
       		call aims_allocate(ovl_fact22_d_rcv_SR, n_atoms2, 'ovl_fact22_d_rcv_SR')
       		call aims_allocate(AS_ovl_fact22_d_rcv_SR, n_atoms2, 'AS_ovl_fact22_d_rcv_SR')
       end if
    endif

    !---------------------------------------------------------------------------

    ! Calculate density matrix, do a Fourier transform and store it in dm_cols.

    dm_cols = 0.

    if(use_scalapack) then

      ! Calculate density matrix using scalapack, it is stored in dm_scalapack_real/complex

      if (present(density_matrix_real_in) .or. &
           & present(density_matrix_complex_in)) then
         if (real_eigenvectors) then
            dm_scalapack_real = density_matrix_real_in
         else
            dm_scalapack_complex = density_matrix_complex_in
         end if
      else
         do i_spin = 1,n_spin
            if(elsi_read_dm) then
               call get_timestamps(t_read,tmp)

               if(real_eigenvectors) then
                  call elsi_restart_scalapack(i_spin,&
                       dm_scalapack_real(:,:,i_spin))

                  dm_scalapack_real(:,:,i_spin) = &
                     dm_scalapack_real(:,:,i_spin)/k_weights(my_k_point)
               else
                  call elsi_restart_scalapack(i_spin,&
                       dm_scalapack_complex(:,:,i_spin))

                  dm_scalapack_complex(:,:,i_spin) = &
                     dm_scalapack_complex(:,:,i_spin)/k_weights(my_k_point)
               endif

               call get_times(t_read,tmp,tot_time_matrix_io,&
                    tot_clock_time_matrix_io)

               write(info_str,"(2X,A)") "Finished reading density matrices"//&
                  " from file"
               call localorb_info(info_str,use_unit)
               write(info_str,"(2X,A,F10.3,A)") "| Time : ",t_read," s"
               call localorb_info(info_str,use_unit)
               write(info_str,"(A)") ""
               call localorb_info(info_str,use_unit)

               elsi_read_dm_done = .true.
            else
               call construct_dm_scalapack(occ_numbers,i_spin)

               ! Only upper half is ready, we need both!
               if(real_eigenvectors) then
                  call set_full_matrix_real(dm_scalapack_real(:,:,i_spin))
               else
                  call set_full_matrix_complex(dm_scalapack_complex(:,:,i_spin))
               endif
            endif
         enddo
      endif

      ! Fourier transform

      call aims_allocate(dm_tmp,n_basis,n_cells_bvk,'dm_tmp')

      do i_spin = 1, n_spin

        do i_atom = 1, n_atoms2
          i_basis_s = atom2basis_off2(i_atom) + 1
          i_basis_e = atom2basis_off2(i_atom) + atom2basis_len2(i_atom)
          do i_basis = i_basis_s, i_basis_e

            dm_tmp(:,:) = 0.

            if(l_col(i_basis) > 0) then
              do i_row = 1, n_basis
                if(l_row(i_row)>0) then
                  if(real_eigenvectors)then
                    dm_tmp(i_row,:) = dm_scalapack_real(l_row(i_row),l_col(i_basis),i_spin) &
                                      * dble (k_phase_exx(:,my_k_point))*k_weights(my_k_point)
                  else
                    dm_tmp(i_row,:) = dble(dm_scalapack_complex(l_row(i_row),l_col(i_basis),i_spin)&
                                      * conjg(k_phase_exx(:,my_k_point))*k_weights(my_k_point))
                  endif
                endif
              enddo
            endif

            call sync_vector(dm_tmp, size(dm_tmp))

            if(my_basis_off(i_atom) >= 0)  &
               dm_cols(:,i_basis-atom2basis_off2(i_atom)+my_basis_off(i_atom),:,i_spin) = dm_tmp(:,:)

          end do
        end do
      end do

      call aims_deallocate( dm_tmp, "dm_tmp" )

    else

      if (present(density_matrix_real_in) .or. &
           & present(density_matrix_complex_in)) &
           & call aims_stop("optional density matrix argument to &
           &evaluate_exchange_matr_realspace_p0 only works with use_scalapack")
      if(elsi_read_dm) then
         call read_densmat_hf()
      else
         ! Calculate density matrix and store it in
         ! hf_exchange_matr_real/complex as intermediate storage
         if(use_symmetry_reduced_spg .and. get_full_density_first) then
            call evaluate_densmat_hf_sym(KS_eigenvector,KS_eigenvector_complex,&
                    occ_numbers)
         else
            call evaluate_densmat_hf(KS_eigenvector,KS_eigenvector_complex,&
                    occ_numbers)
         endif
      endif

      ! Fourier transform

      call aims_allocate(dm_tmp,n_basis,n_basis,'dm_tmp')

      if(use_symmetry_reduced_spg)then
        call FT_densmat_sym(my_n_atoms,my_n_basis,my_atom_list,my_basis_off,dm_tmp,dm_cols)
      else
        do i_cell = 1, n_cells_bvk
          do i_spin = 1, n_spin

            dm_tmp(:,:) = 0.

            i_k = 0
            do i_k_point = 1, n_k_points,1
              if(myid == MOD(i_k_point, n_tasks)) then
                i_k = i_k + 1
                if(real_eigenvectors)then
                  dm_tmp(:,:) = dm_tmp(:,:) + &
                      hf_exchange_matr_real   (:,:,i_k,i_spin)*dble (k_phase_exx(i_cell,i_k_point)*k_weights(i_k_point))
                else
                  dm_tmp(:,:) = dm_tmp(:,:) + &
                  dble(hf_exchange_matr_complex(:,:,i_k,i_spin)*conjg(k_phase_exx(i_cell,i_k_point)*k_weights(i_k_point)))
                endif
              endif
            enddo

            call sync_vector(dm_tmp, size(dm_tmp))

            do i_my_atom = 1, my_n_atoms
              i_atom = my_atom_list(i_my_atom)
              i_basis_s = atom2basis_off2(i_atom) + 1
              i_basis_e = atom2basis_off2(i_atom) + atom2basis_len2(i_atom)
              my_s = my_basis_off(i_atom) + 1
              my_e = my_basis_off(i_atom) + atom2basis_len2(i_atom)
              !SVL this is original version:
              !dm_cols(:,my_s:my_e,i_cell,i_spin) = dm_tmp(:,i_basis_s:i_basis_s)
              !New version:
              dm_cols(:,my_s:my_e,i_cell,i_spin) = dm_tmp(:,i_basis_s:i_basis_e)
            enddo

          enddo
        enddo
      endif

      call aims_deallocate( dm_tmp, "dm_tmp" )

    endif
    !---------------------------------------------------------------------------

    ! Get max absolute value of dm_cols per atom/cell (for my atoms)

    max_abs_dm_cols = 0
    do i_spin = 1,n_spin
      do i_cell = 1,n_cells_bvk
        do i_atom_r1 = 1, n_atoms2
          nr1_off    = atom2basis_off2(i_atom_r1)
          nr1_len    = atom2basis_len2(i_atom_r1)
          do i_my_atom = 1, my_n_atoms
            i_atom_r2  = my_atom_list(i_my_atom)
            nr2_off    = my_basis_off(i_atom_r2)
            nr2_len    = atom2basis_len2(i_atom_r2)
            max_abs_dm_cols(i_atom_r1,i_atom_r2,i_cell) = max(max_abs_dm_cols(i_atom_r1,i_atom_r2,i_cell), &
              maxval(abs(dm_cols(nr1_off+1:nr1_off+nr1_len,nr2_off+1:nr2_off+nr2_len,i_cell,i_spin))))
          enddo
        enddo
      enddo
    enddo

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Algorithm used below for the HF Matrix:
    !
    ! The Fock Matrix consists of 4 different parts comming from the 2 different
    ! parts of each overlap integral (called OVL1 and OVL2 below):
    !
    ! HF11(l1,l2,icf) = OVL1(:,l1,r1,ic1) * CM(:,:,l1,l2,icf        ) * OVL1(:,l2,r2,ic2)
    !                                     * DM(    r1,r2,icf-ic1+ic2)
    !
    ! HF12(l1,l2,icf) = OVL1(:,l1,r1,ic1) * CM(:,:,l1,r2,icf+ic2    ) * OVL2(:,l2,r2,ic2)
    !                                     * DM(    r1,r2,icf-ic1+ic2)
    !
    ! HF21(l1,l2,icf) = OVL2(:,l1,r1,ic1) * CM(:,:,r1,l2,icf-ic1    ) * OVL1(:,l2,r2,ic2)
    !                                     * DM(    r1,r2,icf-ic1+ic2)
    !
    ! HF22(l1,l2,icf) = OVL2(:,l1,r1,ic1) * CM(:,:,r1,r2,icf-ic1+ic2) * OVL2(:,l2,r2,ic2)
    !                                     * DM(    r1,r2,icf-ic1+ic2)
    ! (all to be summed over r1, r2)
    !
    ! ATTENTION: The ':' in the above notation (indicating basbas functions) are left away below.
    ! They have to be considered implicitly!
    !
    ! Please note that HF12 and HF21 deliver a symmetrical result and thus
    ! only one of both needs to be calculated.
    !
    ! Please note also that OVL2(l,r,ic) = OVL1(r,l,-ic) and thus only OVL1 is stored.
    ! Using this for HF12/HF22 gives:
    !
    ! HF12(l1,l2,icf) = OVL2(r1,l1,-ic1) * CM(l1,r2,icf+ic2    ) * OVL1(r2,l2,-ic2)
    !                                    * DM(r1,r2,icf-ic1+ic2)
    !
    ! HF22(l1,l2,icf) = OVL1(r1,l1,-ic1) * CM(r1,r2,icf-ic1+ic2) * OVL1(r2,l2,-ic2)
    !                                    * DM(r1,r2,icf-ic1+ic2)
    !
    ! Exchanging l1 <-> r1, l2 <-> r2, ic1 <-> -ic1, ic2 <-> -ic2 in the
    ! above equation yields
    !
    ! HF12(r1,r2,icf) = OVL2(l1,r1,ic1) * CM(r1,l2,icf-ic2    ) * OVL1(l2,r2,ic2)
    !                                   * DM(l1,l2,icf+ic1-ic2)
    !
    ! HF22(r1,r2,icf) = OVL1(l1,r1,ic1) * CM(l1,l2,icf+ic1-ic2) * OVL1(l2,r2,ic2)
    !                                   * DM(l1,l2,icf+ic1-ic2)
    ! (to be summed over l1, l2)
    !
    ! Together with HF11 from above
    !
    ! HF11(l1,l2,icf) = OVL1(l1,r1,ic1) * CM(l1,l2,icf        ) * OVL1(l2,r2,ic2)
    !                                   * DM(r1,r2,icf-ic1+ic2)
    !
    ! this is what is actually calculated.
    ! The OVL terms at the left are broadcast in every iteration, the OVL terms
    ! at the right are locally stored.
    !
    ! This looks complicated, but it has some advantages when taking into account that
    ! - the data distriution is over l2 (ie. the number of elements in l2 direction is "small" locally)
    ! - the outermost loop is over r1
    ! - the right instance of OVL1 is kept locally
    !
    ! These advantages are:
    ! * only OVL1 (not OVL2) appears on the right (which is locally stored)
    ! * OVL1 on the right is needed only for local l2
    ! * the Coulomb matrices are needed only for local l2
    ! * HF11 is needed only for local l2 (and needs no exchange)
    ! * same for DM(l1,l2,.)
    ! * HF12/HF22 must be exchanged (they are going through all indices),
    !   but only once during the sweep over r1
    ! * same for DM(r1,r2,.)
    !
    ! Here is the pseudocode for what is done below.
    !
    ! Please note:
    ! - All OVL's and all TMP's have an implicit additional index (for basbas function)
    ! - CM(l1,l2,ic) is in reality CM(:,:,l1,l2,ic) (for basbas functions)
    ! - Multiplying OVL with CM is in reality a matrix-matrix product
    ! - Multiplying OVL with TMP's or two TMP's is a matrix-matrix product
    !
    ! DO r1 = 1, n_basis
    !
    !   Broadcast DM(r1,:,:)   - called DM_ROW(:,:) below
    !   Broadcast OVL1/2(:,r1,:) - called OVL_RCV1/2(:,:) below
    !
    !   FORALL l2a IN my atoms (the distribution of OVL and CM is over l2)
    !
    !     l2b = Number of basis functions of l2a
    !
    !     TMP1(1:l2b,ic) = SUM over r2 ( OVL1(1:l2b,r2,ic2)*DM_ROW(r2,ic+ic2) )
    !
    !     DO ic  = 1, n_cells_bvk
    !     DO l1a = 1, n_atoms ! The loop nest l1a/ic1 is collapsed in the code!
    !     DO ic1 = 1, n_cells_bvk
    !
    !       l1b = Number of basis functions of l1a
    !
    !       TMPP(1:l1b) = OVL_RCV1(1:l1b,ic1) * CM(l1,l2,ic) ! This can be used for HF11 and HF22
    !
    !       HF11(1:l1b,1:l2b,ic) = HF11(1:l1b,1:l2b,ic) + TMPP(1:l1b) * TMP1(1:l2b,ic-ic1)
    !
    !       TMP2(1:l2b,ic-ic1) = TMP2(1:l2b,ic-ic1) + TMPP(1:l1b) * DM(1:l1b,1:l2b,ic)
    !
    !       TMPX(1:l2b,ic-ic1) = TMPX(1:l2b,ic-ic1) + OVL_RCV2(1:l1b,ic1)*DM(1:l1b,1:l2b,ic)
    !
    !     ENDDO
    !     ENDDO
    !     ENDDO
    !
    !     DO ic = 1, n_cells_bvk
    !       TMP1(1:l2b,ic) = TMPX(1:l2b,ic) * CM(r1,l2,ic)
    !     ENDDO
    !
    !     DO r2a = 1, n_atoms     ! The loop nest r2a/ic2 is collapsed in the code!
    !     DO ic2 = 1, n_cells_bvk
    !     DO ic  = 1, n_cells_bvk
    !
    !       r2b = Number of basis functions of r1a
    !
    !       ! The following is done together in the code since there is no need to
    !       ! distiguish HF12 and HF22
    !       HF12(r1,1:r2b,ic+ic2) = SUM over l2b ( TMP1(l2b,ic) * OVL1(l2b,1:r2b,ic2) )
    !       HF22(r1,1:r2b,ic+ic2) = SUM over l2b ( TMP2(l2b,ic) * OVL1(l2b,1:r2b,ic2) )
    !
    !     ENDDO
    !     ENDDO
    !     ENDDO
    !
    !   ENDDO ! l2a
    ! ENDDO ! r1
    !
    ! This code has the following remarkable features:
    !
    ! - it scales with n_atoms^3 - no loop nest involving atoms or basis functions is deeper than 3
    ! - The product OVL_RCV1(1:l1b,ic1) * CM(l1,l2,ic) is used twice, the second product
    !   is needed much less times (time for it can be neglected)
    !
    !
    ! When unfolding the above code, you get:
    !
    ! TMP1(l2,ic) = OVL1(l2,r2,ic2)*DM(r1,r2,ic+ic2), summed over r2
    !
    ! TMPP(l1) = OVL1(l1,r1,ic1)*CM(l1,l2,ic)
    !
    ! HF11(l1,l2,ic) = TMPP(l1)*TMP1(l2,ic-ic1), summed over r1 in the outermost loop
    !                = OVL1(l1,r1,ic1)*CM(l1,l2,ic)*OVL1(l2,r2,ic2)*DM(r1,r2,ic-ic1+ic2), summed over r2 + r1
    !
    ! TMP2(l2,ic-ic1) = TMPP(l1)*DM(l1,l2,ic), summed over l1
    !                 = OVL1(l1,r1,ic1)*CM(l1,l2,ic)*DM(l1,l2,ic), summed over l1
    !
    ! TMPX(l2,ic-ic1) = OVL2(r1,l1,ic1)*DM(l1,l2,ic), summed over l1
    !
    ! TMP1(l2,ic) = TMPX(l2,ic)*CM(r1,l2,ic)
    !             = OVL2(r1,l1,ic1)*DM(l1,l2,ic+ic1)*CM(r1,l2,ic), summed over l1
    !
    ! HF12(r1,r2,ic+ic2) = TMP1(l2,ic)*OVL1(l2,r2,ic2), summed over l2
    !                    = OVL2(r1,l1,ic1)*DM(l1,l2,ic+ic1)*CM(r1,l2,ic)*OVL1(l2,r2,ic2), summed over l1 + l2
    ! Substituting ic by ic-ic2:
    ! HF12(r1,r2,ic) = OVL2(r1,l1,ic1)*DM(l1,l2,ic+ic1-ic2)*CM(r1,l2,ic-ic2)*OVL1(l2,r2,ic2), summed over l1 + l2
    !
    ! HF22(r1,r2,ic+ic2) = TMP2(l2,ic)*OVL1(l2,r2,ic2), summed over l2
    !                    = OVL1(l1,r1,ic1)*CM(l1,l2,ic+ic1)*DM(l1,l2,ic+ic1)*OVL1(l2,r2,ic2), summed over l1 + l2
    ! Substituting ic by ic-ic2:
    ! HF22(r1,r2,ic) = OVL1(l1,r1,ic1)*CM(l1,l2,ic+ic1-ic2)*DM(l1,l2,ic+ic1-ic2)*OVL1(l2,r2,ic2), summed over l1 + l2
    !
    ! i.e. exactly the definitions for HF from above!
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Exchange energy gradients:
    !
    ! exx_ene = 0.5*HF(l1,l2,ic)*DM(l1,l2,ic) summed over l1, l2, ic
    !
    ! It can be easily shown that the contributions from HF11 and HF22 are identical
    ! (as are the contributions from HF21 and HF12) so it is sufficient to set
    !
    ! exx_ene = (HF11 + HF12)*DM   or   exx_ene = (HF22 + HF12)*DM (we will use both!)
    !
    ! HF has the form OVL*CM*OVL*DM and since the gradient of DM is 0 (???):
    !
    ! d_HF = (d_OVL*CM*OVL + OVL*d_CM*OVL + OVL*CM*d_OVL)*DM
    !
    ! Again, the 1st and 3rd term yield identical results, so that we will calculate
    ! d_HF = (OVL*d_CM*OVL + 2*OVL*CM*d_OVL)*DM
    !
    ! For the first term OVL*d_CM*OVL we use exx_ene = (HF11 + HF12)*DM, which ends up in
    ! the two contributions
    !
    !     OVL1(l1,r1,ic1) * d_CM(l1,l2,icf        ) * OVL1(l2,r2,ic2)
    !                     *   DM(r1,r2,icf-ic1+ic2)
    !                     *   DM(l1,l2,icf)
    !   + OVL2(l1,r1,ic1) * d_CM(r1,l2,icf-ic2    ) * OVL1(l2,r2,ic2)
    !                     *   DM(l1,l2,icf+ic1-ic2)
    !                     *   DM(r1,r2,icf)
    !
    ! For the first part, we sum up
    !
    !     OVL1(l1,r1,ic1) * d_CM(l1,l2,icf        ) * OVL1(l2,r2,ic2)
    !                     *   DM(r1,r2,icf-ic1+ic2)
    !
    ! exactly as we do for HF11 (re-using TMP1 for OVL1(l2,r2,ic2)*DM(r1,r2,icf-ic1+ic2))
    ! but using d_CM instead of CM.
    ! The remaining multiplication with DM(l1,l2,icf) is trivial.
    !
    ! For the second part, we sum up
    !
    !     OVL2(l1,r1,ic1) * d_CM(r1,l2,icf-ic2    )
    !                     *   DM(l1,l2,icf+ic1-ic2)
    !
    ! exactly as we do for HF12 (re-using TMPX for OVL2(l1,r1,ic1)*DM(l1,l2,icf+ic1-ic2))
    ! but using d_CM instead of CM.
    ! For the remaining
    !
    !     OVL1(l2,r2,ic2)*DM(r1,r2,icf)
    !
    ! we can use the precalculated TMP1 (with some index magic for the cells)
    !
    ! For the second Term 2*OVL*CM*d_OVL we use exx_ene = (HF22 + HF12)*DM
    ! which ends up in
    !
    !
    !     OVL2(l1,r1,ic1) * CM(r1,l2,icf-ic2    ) * d_OVL1(l2,r2,ic2)
    !                     * DM(l1,l2,icf+ic1-ic2)
    !                     * DM(r1,r2,icf)
    !
    !   + OVL1(l1,r1,ic1) * CM(l1,l2,icf+ic1-ic2) * d_OVL1(l2,r2,ic2)
    !                     * DM(l1,l2,icf+ic1-ic2)
    !                     * DM(r1,r2,icf)
    !
    ! Here, the contributions
    !
    !     OVL2(l1,r1,ic1) * CM(r1,l2,icf-ic2    )
    !                     * DM(l1,l2,icf+ic1-ic2)
    !   + OVL1(l1,r1,ic1) * CM(l1,l2,icf+ic1-ic2)
    !                     * DM(l1,l2,icf+ic1-ic2)
    !
    ! are already done when calculating HF12 and HF22, so we have only to
    ! multiply with d_OVL1 (instead of OVL1 in the case of HF12/HF22)
    ! and then do a trivial multiplication with DM(r1,r2,icf).
    !
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    ! A lot of code fragments have been copied right below their original line but changed slighty by some variables.
    ! This needs to be done for the LR-wPBEh functional in which two fock_matr need to be calculated. So similar to the initalisation above
    ! a goto is used to run the routines again. This is not really nice coding, but again the easiest without coping everything. This also
    ! saves some time (in theory and not tested) as some arrays only need to be populated once.
    ! The fragments normaly start with the if statement: if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then

    fock_matrix  = 0.
    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
    	fock_matrix_SR = 0.
    	exx_ene_SR = 0.
    	if(calc_deriv) d_exx_ene_SR = 0.
    end if
    exx_ene = 0.
    if(calc_deriv) d_exx_ene = 0.
    if(AS_stress_on) AS_EXX_stress_local(1:9) = 0.0d0
    time_mult_coul_mat = 0.
! End timing block,   Pre/Post #1
ttt_other = mpi_wtime()
times(6) = times(6) + ttt_other - ttt0
! Begin timing block, Other #2

    ! Get maximum factor with which  OVL1(l1,r1,ic1) will be multiplied to get HF22
    !
    ! From above:
    ! HF22(r1,r2,icf) = OVL1(l1,r1,ic1) * CM(l1,l2,icf+ic1-ic2) * OVL1(l2,r2,ic2)
    !                                   * DM(l1,l2,icf+ic1-ic2)
    ! The factor is independent of r1 and can be calculated outside the i_basis_r1 loop

    ovl_fact22(:) = 0   ! Factor for Fock matrix
    ovl_fact22_d(:) = 0 ! Factor for gradients
    AS_ovl_fact22_d(:) = 0.0d0 ! Factor for analytical stress

    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
    	ovl_fact22_SR(:) = 0   ! Factor for Fock matrix
    	ovl_fact22_d_SR(:) = 0 ! Factor for gradients
    	AS_ovl_fact22_d_SR(:) = 0.0d0 ! Factor for analytical stress
    end if

    if (.not.(use_mpi_in_place)) then
       ovl_fact22_rcv(:) = 0
       ovl_fact22_d_rcv(:) = 0
       AS_ovl_fact22_d_rcv(:) = 0.0d0

       if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
       		ovl_fact22_rcv_SR(:) = 0
      	 	ovl_fact22_d_rcv_SR(:) = 0
       		AS_ovl_fact22_d_rcv_SR(:) = 0.0d0
       end if
    end if


    do i_atom_l1 = 1, n_atoms2
      do i_my_atom = 1, my_n_atoms
        i_atom_l2 = my_atom_list(i_my_atom)
        s = 0.
        do i_cell = 1, n_cells_bvk
          s = max(s,coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell))
        enddo
        ovl_fact22(i_atom_l1) = max(ovl_fact22(i_atom_l1),s*max_norm_ovlp3fn_per_latom(i_atom_l2))
        if(calc_deriv) &
          ovl_fact22_d(i_atom_l1) = max(ovl_fact22_d(i_atom_l1),s*max_norm_d_ovlp3fn_per_latom(i_atom_l2))
        if(AS_stress_on) &
          AS_ovl_fact22_d(i_atom_l1) = max(AS_ovl_fact22_d(i_atom_l1),s*AS_max_norm_d_ovlp3fn_per_latom(i_atom_l2))
        ! ovl_fact22 must still be multiplied with MAXVAL(DM(r1,:,:)) (inside i_basis_r1 loop)
      enddo
    enddo

    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
        do i_atom_l1 = 1, n_atoms2
		  do i_my_atom = 1, my_n_atoms
		    i_atom_l2 = my_atom_list(i_my_atom)
		    s = 0.
		    do i_cell = 1, n_cells_bvk
		      s = max(s,coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell))
		    enddo
		    ovl_fact22_SR(i_atom_l1) = max(ovl_fact22_SR(i_atom_l1),s*max_norm_ovlp3fn_per_latom_SR(i_atom_l2))
		    if(calc_deriv) &
		      ovl_fact22_d_SR(i_atom_l1) = max(ovl_fact22_d_SR(i_atom_l1),s*max_norm_d_ovlp3fn_per_latom_SR(i_atom_l2))
		    if(AS_stress_on) &
		      AS_ovl_fact22_d_SR(i_atom_l1) = max(AS_ovl_fact22_d_SR(i_atom_l1),s*AS_max_norm_d_ovlp3fn_per_latom_SR(i_atom_l2))
		    ! ovl_fact22 must still be multiplied with MAXVAL(DM(r1,:,:)) (inside i_basis_r1 loop)
		  enddo
		enddo
    end if

    if (use_mpi) then
       if (.not.(use_mpi_in_place)) then
          call mpi_allreduce(ovl_fact22,ovl_fact22_rcv,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
          ovl_fact22 = ovl_fact22_rcv
          if(calc_deriv) then
             call mpi_allreduce(ovl_fact22_d,ovl_fact22_d_rcv,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
             ovl_fact22_d = ovl_fact22_d_rcv
          endif
          if(AS_stress_on) then
             call mpi_allreduce(AS_ovl_fact22_d,AS_ovl_fact22_d_rcv,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
             AS_ovl_fact22_d = AS_ovl_fact22_d_rcv
          endif
       else
          call mpi_allreduce(mpi_in_place,ovl_fact22,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
          if(calc_deriv) then
             call mpi_allreduce(mpi_in_place,ovl_fact22_d,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
          endif
          if(AS_stress_on) then
             call mpi_allreduce(mpi_in_place,AS_ovl_fact22_d,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
          endif
       end if

       if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
		   if (.not.(use_mpi_in_place)) then
		      call mpi_allreduce(ovl_fact22_SR,ovl_fact22_rcv_SR,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
		      ovl_fact22_SR = ovl_fact22_rcv_SR
		      if(calc_deriv) then
		         call mpi_allreduce(ovl_fact22_d_SR,ovl_fact22_d_rcv_SR,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
		         ovl_fact22_d_SR = ovl_fact22_d_rcv_SR
		      endif
		      if(AS_stress_on) then
		         call mpi_allreduce(AS_ovl_fact22_d_SR,AS_ovl_fact22_d_rcv_SR,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
		         AS_ovl_fact22_d_SR = AS_ovl_fact22_d_rcv_SR
		      endif
		   else
		      call mpi_allreduce(mpi_in_place,ovl_fact22_SR,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
		      if(calc_deriv) then
		         call mpi_allreduce(mpi_in_place,ovl_fact22_d_SR,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
		      endif
		      if(AS_stress_on) then
		         call mpi_allreduce(mpi_in_place,AS_ovl_fact22_d_SR,n_atoms2,mpi_real8,mpi_max,mpi_comm_global,mpierr)
		      endif
		   end if
       end if
    end if


    !---------------------------------------------------------------------------
    ! Start of outermost loop over i_basis_r1
    !---------------------------------------------------------------------------

    i_atom_r1_prev = -1 ! no previous atom
	write(info_str,'(A,I8)') 'Exact exchange progress report - outermost basis function loop i_basis_r1: 1 ..', n_basis
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )

! End timing block,   Other #2
times_other(2) = times_other(2) + mpi_wtime() - ttt_other

    do i_basis_r1 = 1, n_basis
! Begin timing block, Other #3
ttt_other = mpi_wtime()

      if ( mod(i_basis_r1,10).eq.0) then
        write(info_str,'(I8)') i_basis_r1
        call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      else if ( mod(i_basis_r1,10).eq.1) then
        write(info_str,'(A,I8)') '| ', i_basis_r1
        call localorb_info ( info_str,use_unit,'(2X,A,$)', OL_norm  )
      else
        write(info_str,'(I8)') i_basis_r1
        call localorb_info ( info_str,use_unit,'(2X,A,$)', OL_norm  )
      end if

      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
		lc_wpbeh_lr_run=.false.
	  else
		lc_wpbeh_lr_run=.true.
	  end if

124   i_atom_r1  = basis_atom2(i_basis_r1)

      nr1_off    = atom2basis_off2(i_atom_r1)
      nr1_len    = atom2basis_len2(i_atom_r1)

      ! Number of atom pairs for this i_basis_r1

      n_pair_1 = count(pair_flag_bvk(:,i_atom_r1,:))

      ! Get row of density matrix
      ! Since the density matrix is symmetric, we can use a column of the inverse cell

      dm_row = 0
      if(my_atom_id == 0 .and. my_basis_off(i_atom_r1) >= 0) then
        do i_cell = 1, n_cells_bvk
          i_basis = i_basis_r1-atom2basis_off2(i_atom_r1)+my_basis_off(i_atom_r1)
          dm_row(:,i_cell,:) = dm_cols(:,i_basis,inv_cell_bvk(i_cell),:)
        enddo
      endif
! End timing block,   Other #3
times_other(3) = times_other(3) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #1
      call time_sync_start
      ! We could use a broadcast instead of sync_vector, but this is not time critical here
      call sync_vector(dm_row,size(dm_row))
      call time_sync_end
! End timing block,   Sync/Imbalance #1
! Begin timing block, Other #4
ttt_other = mpi_wtime()

      ! Calculate maximum absolute value of dm_row for all atoms.

      max_abs_dm_row = 0
      do i_atom = 1, n_atoms2
        na_off = atom2basis_off2(i_atom)
        na_len = atom2basis_len2(i_atom)
        do i_cell_dm = 1, n_cells_bvk
          max_abs_dm_row(i_atom,i_cell_dm) = maxval(abs(dm_row(na_off+1:na_off+na_len,i_cell_dm,:)))
        enddo
      enddo

      max_abs_dm_row_tot = maxval(max_abs_dm_row(:,:))

      ! Get maximum factor with which  OVL2(l1,r1,ic1) will be multiplied to get HF12
      !
      ! From above:
      ! HF12(r1,r2,icf) = OVL2(l1,r1,ic1) * CM(r1,l2,icf-ic2    ) * OVL1(l2,r2,ic2)
      !                                   * DM(l1,l2,icf+ic1-ic2)
      ! This factor depends only on i_atom_r1, not on i_basis_r1 and needs to be calculated
      ! only once for every atom

      if(i_atom_r1_prev .ne. i_atom_r1) then
        ovl_fact12(:,:) = 0   ! Factor for Fock matrix
        if(calc_deriv) ovl_fact12_d(:,:) = 0 ! Factor for gradients
        if(AS_stress_on) AS_ovl_fact12_d(:,:) = 0.0d0 ! Factor for analytical stress

        if (.not.(use_mpi_in_place)) then
           ovl_fact12_rcv(:,:) = 0
           ovl_fact12_d_rcv(:,:) = 0
           AS_ovl_fact12_d_rcv(:,:) = 0.0d0
        end if

        do i_atom_l1 = 1, n_atoms2
          do i_cell_1 = 1, n_cells_bvk
            do i_my_atom = 1, my_n_atoms
              i_atom_l2 = my_atom_list(i_my_atom)
              s = 0.
              do i_cell = 1, n_cells_bvk
                i_cell_p1 = add_cells(i_cell,i_cell_1)
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
                	s = max(s,coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1))
                else
	                s = max(s,coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1))
	            end if
              enddo
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          ovl_fact12(i_atom_l1,i_cell_1) = &
		            max(ovl_fact12(i_atom_l1,i_cell_1),s*max_norm_ovlp3fn_per_latom_SR(i_atom_l2))
              else
		  		  ovl_fact12(i_atom_l1,i_cell_1) = &
		            max(ovl_fact12(i_atom_l1,i_cell_1),s*max_norm_ovlp3fn_per_latom(i_atom_l2))
              end if

              if(calc_deriv) then
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(ovl_fact12_d(i_atom_l1,i_cell_1),s*max_norm_d_ovlp3fn_per_latom_SR(i_atom_l2))
		        else
		        	ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(ovl_fact12_d(i_atom_l1,i_cell_1),s*max_norm_d_ovlp3fn_per_latom(i_atom_l2))
		        end if
                s = 0.
                do i_cell = 1, n_cells_bvk
                  i_cell_p1 = add_cells(i_cell,i_cell_1)
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
	                  s = max(s,d_coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1))
	              else
	              	  s = max(s,d_coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1))
	              end if
                enddo
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(ovl_fact12_d(i_atom_l1,i_cell_1),s*max_norm_ovlp3fn_per_latom_SR(i_atom_l2))
                else
		            ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(ovl_fact12_d(i_atom_l1,i_cell_1),s*max_norm_ovlp3fn_per_latom(i_atom_l2))
                end if
              endif

              if (AS_stress_on) then
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            AS_ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(AS_ovl_fact12_d(i_atom_l1,i_cell_1),s*AS_max_norm_d_ovlp3fn_per_latom_SR(i_atom_l2))
                else
                	AS_ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(AS_ovl_fact12_d(i_atom_l1,i_cell_1),s*AS_max_norm_d_ovlp3fn_per_latom(i_atom_l2))
                end if
                s = 0.0d0
                do i_cell = 1, n_cells_bvk, 1
                  i_cell_p1 = add_cells(i_cell,i_cell_1)
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
                  	s = max(s,AS_d_coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1))
                  else
                  	s = max(s,AS_d_coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1))
                  end if
                end do
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            AS_ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(AS_ovl_fact12_d(i_atom_l1,i_cell_1),s*max_norm_ovlp3fn_per_latom_SR(i_atom_l2))
                else
		        	AS_ovl_fact12_d(i_atom_l1,i_cell_1) = &
		              max(AS_ovl_fact12_d(i_atom_l1,i_cell_1),s*max_norm_ovlp3fn_per_latom(i_atom_l2))
                end if
              end if
            enddo
          enddo
        enddo

        ! ovl_fact12_d still needs to be multiplied by max_abs_dm_row_tot
        if (calc_deriv) ovl_fact12_d(:,:) = ovl_fact12_d(:,:)*max_abs_dm_row_tot
        if (AS_stress_on) then
          AS_ovl_fact12_d(:,:) = AS_ovl_fact12_d(:,:)*max_abs_dm_row_tot
        end if

! End timing block,   Other #4 (if fork was taken)
times_other(4) = times_other(4) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #2
        call time_sync_start

        if (use_mpi) then
           if (.not.(use_mpi_in_place)) then
              call mpi_allreduce(ovl_fact12,ovl_fact12_rcv,n_atoms2*n_cells_bvk,&
                                 mpi_real8,mpi_max,mpi_comm_global,mpierr)
              ovl_fact12 = ovl_fact12_rcv
              if(calc_deriv) then
                 call mpi_allreduce(ovl_fact12_d,ovl_fact12_d_rcv,n_atoms2*n_cells_bvk,&
                                    mpi_real8,mpi_max,mpi_comm_global,mpierr)
                 ovl_fact12_d = ovl_fact12_d_rcv
              endif
              if(AS_stress_on) then
                 call mpi_allreduce(AS_ovl_fact12_d,AS_ovl_fact12_d_rcv,n_atoms2*n_cells_bvk,&
                                    mpi_real8,mpi_max,mpi_comm_global,mpierr)
                 AS_ovl_fact12_d = AS_ovl_fact12_d_rcv
              endif
           else
              call mpi_allreduce(mpi_in_place,ovl_fact12,n_atoms2*n_cells_bvk,&
                                 mpi_real8,mpi_max,mpi_comm_global,mpierr)
              if(calc_deriv) then
                 call mpi_allreduce(mpi_in_place,ovl_fact12_d,n_atoms2*n_cells_bvk,&
                                    mpi_real8,mpi_max,mpi_comm_global,mpierr)
              endif
              if(AS_stress_on) then
                 call mpi_allreduce(mpi_in_place,AS_ovl_fact12_d,n_atoms2*n_cells_bvk,&
                                    mpi_real8,mpi_max,mpi_comm_global,mpierr)
              endif
           end if
        end if

        call time_sync_end
! End timing block,   Sync/Imbalance #2

        i_atom_r1_prev = i_atom_r1
      else
! End timing block,   Other #4 (if fork was not taken)
times_other(4) = times_other(4) + mpi_wtime() - ttt_other
      endif
! Begin timing block, Other #5
ttt_other = mpi_wtime()

      ! Get maximum factor with which  OVL1(l1,r1,ic1) will be multiplied to get HF11
      !
      ! From above:
      ! HF11(l1,l2,icf) = OVL1(l1,r1,ic1) * CM(l1,l2,icf        ) * OVL1(l2,r2,ic2)
      !                                   * DM(r1,r2,icf-ic1+ic2)

      ovl_fact11(:,:) = 0   ! Factor for Fock matrix
      if(calc_deriv) ovl_fact11_d(:,:) = 0 ! Factor for gradients
      if(AS_stress_on) AS_ovl_fact11_d(:,:) = 0.0d0 ! Factor for analytical stress

      if (.not.(use_mpi_in_place)) then
         ovl_fact11_rcv(:,:) = 0
         ovl_fact11_d_rcv(:,:) = 0
         AS_ovl_fact11_d_rcv(:,:) = 0.0d0
      endif

      do i_my_atom = 1, my_n_atoms
        i_atom_l2 = my_atom_list(i_my_atom)
        max_norm_tmp(:) = 0
        do i_atom_pair = my_atom_pair_s(i_my_atom), my_atom_pair_e(i_my_atom)
          i_atom_r2 = pair_list(2,i_atom_pair)
          i_cell_2  = pair_list(3,i_atom_pair)
          do i_cell = 1, n_cells_bvk
            i_cell_dm = add_cells(i_cell,i_cell_2)
            if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		        max_norm_tmp(i_cell) = max(max_norm_tmp(i_cell),max_abs_dm_row(i_atom_r2,i_cell_dm) * &
                                                        max_norm_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2))
            else
                max_norm_tmp(i_cell) = max(max_norm_tmp(i_cell),max_abs_dm_row(i_atom_r2,i_cell_dm) * &
                                                        max_norm_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2))
            end if
          enddo
        enddo
        if (any(max_norm_tmp>max_norm_tmp_old*5.d0) .or. &
            (any(max_abs_dm_cols(:,i_atom_l2,:)>max_abs_dm_cols_old(:,i_atom_l2,:)*5.d0) .and. (calc_deriv .or. AS_stress_on))) then
        do i_atom_l1 = 1, n_atoms2
          do i_cell_1 = 1, n_cells_bvk
            do i_cell = 1, n_cells_bvk
              i_cell_p1 = add_cells(i_cell,i_cell_1)
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          ovl_fact11(i_atom_l1,i_cell_1) = max(ovl_fact11(i_atom_l1,i_cell_1), &
		                                               coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_tmp(i_cell))
              else
		          ovl_fact11(i_atom_l1,i_cell_1) = max(ovl_fact11(i_atom_l1,i_cell_1), &
		                                               coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_tmp(i_cell))
              end if
              if(calc_deriv) then
              	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            ovl_fact11_d(i_atom_l1,i_cell_1) = max(ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                                   d_coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)* &
		                                                   max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)* &
		                                                   max_norm_tmp(i_cell))
                else
		            ovl_fact11_d(i_atom_l1,i_cell_1) = max(ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                                   d_coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)* &
		                                                   max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)* &
		                                                   max_norm_tmp(i_cell))
                end if
              end if
              if(AS_stress_on) then
              	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            AS_ovl_fact11_d(i_atom_l1,i_cell_1) = max(AS_ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                                   AS_d_coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)* &
		                                                   max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)* &
		                                                   max_norm_tmp(i_cell))
                else
		            AS_ovl_fact11_d(i_atom_l1,i_cell_1) = max(AS_ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                                   AS_d_coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)* &
		                                                   max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)* &
		                                                   max_norm_tmp(i_cell))
                end if
              end if
            enddo
          enddo
        enddo
        max_norm_tmp_old = max_norm_tmp
        max_abs_dm_cols_old(:,i_atom_l2,:) = max_abs_dm_cols(:,i_atom_l2,:)
        end if
      enddo

! End timing block,   Other #5
times_other(5) = times_other(5) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #3
      call time_sync_start

      if (use_mpi) then
         if (.not.(use_mpi_in_place)) then
            call mpi_allreduce(ovl_fact11,ovl_fact11_rcv,n_atoms2*n_cells_bvk,&
                               mpi_real8,mpi_max,mpi_comm_global,mpierr)
            ovl_fact11 = ovl_fact11_rcv
            if(calc_deriv) then
               call mpi_allreduce(ovl_fact11_d,ovl_fact11_d_rcv,n_atoms2*n_cells_bvk,&
                                  mpi_real8,mpi_max,mpi_comm_global,mpierr)
               ovl_fact11_d = ovl_fact11_d_rcv
            endif
            if(AS_stress_on) then
               call mpi_allreduce(AS_ovl_fact11_d,AS_ovl_fact11_d_rcv,n_atoms2*n_cells_bvk,&
                                  mpi_real8,mpi_max,mpi_comm_global,mpierr)
               AS_ovl_fact11_d = AS_ovl_fact11_d_rcv
            endif
         else
            call mpi_allreduce(mpi_in_place,ovl_fact11,n_atoms2*n_cells_bvk,&
                               mpi_real8,mpi_max,mpi_comm_global,mpierr)

            if(calc_deriv) then
               call mpi_allreduce(mpi_in_place,ovl_fact11_d,n_atoms2*n_cells_bvk,&
                                  mpi_real8,mpi_max,mpi_comm_global,mpierr)
            endif
            if(AS_stress_on) then
               call mpi_allreduce(mpi_in_place,AS_ovl_fact11_d,n_atoms2*n_cells_bvk,&
                                  mpi_real8,mpi_max,mpi_comm_global,mpierr)
            endif
         end if
      end if

      call time_sync_end
! End timing block,   Sync/Imbalance #3
ttt_other = mpi_wtime()
! Begin timing block, Other #6


      ! Count the number of basis pairs which are needed in OVL1(l1,r1,ic1) and OVL2(l1,r1,ic1)
      ! for the current i_basis_r1 based upon the maximum factors with which they will be
      ! multiplied below.
      !
      ! Set up a list for the atom pairs (in OVL1/OVL2) containing
      ! pair_1_list(:,1) : i_atom_l1 (i.e. left atom of pair, right atom is always i_atom_r1)
      ! pair_1_list(:,2) : i_cell_1
      ! pair_1_list(:,3) : Number of basis pairs in OVL1
      ! pair_1_list(:,4) : Number of basis pairs in OVL2
      ! pair_1_list(:,5) : Offset in ovlp3fn_rcv for OVL1
      ! pair_1_list(:,6) : Offset in ovlp3fn_rcv for OVL2

      i_pair_1 = 0
      pair_1_list(:,:) = 0

      do i_cell_1 = 1, n_cells_bvk
        do i_atom_l1 = 1, n_atoms2

          nl1_off = atom2basis_off2(i_atom_l1)
          nl1_len = atom2basis_len2(i_atom_l1)

          if(pair_flag_bvk(i_atom_l1,i_atom_r1,i_cell_1)) then

            i_pair_1 = i_pair_1 + 1
            pair_1_list(i_pair_1,1) = i_atom_l1
            pair_1_list(i_pair_1,2) = i_cell_1

            ! Number of basis pairs in OVL1
            if(my_pair_flag_bvk(i_atom_l1,i_atom_r1,i_cell_1)) then
              n = pair_offset(i_atom_l1,i_atom_r1,i_cell_1) + (i_basis_r1-nr1_off-1)*nl1_len + 1
              do il1 = 1, nl1_len
              	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            if(ovlp3fn_SR(i_atom_l1)%n(n)*max(ovl_fact11(i_atom_l1,i_cell_1),ovl_fact22_SR(i_atom_l1)) > crit_val .or. &
		               (calc_deriv .and. (ovlp3fn_SR(i_atom_l1)%n(n)*max(ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                           ovl_fact22_d_SR(i_atom_l1)*max_abs_dm_row_tot) > crit_val_deriv)) .or. &
		               (AS_stress_on .and. (ovlp3fn_SR(i_atom_l1)%n(n)*max(AS_ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                           AS_ovl_fact22_d_SR(i_atom_l1)*max_abs_dm_row_tot) > AS_crit_val_deriv))) then
		              pair_1_list(i_pair_1,3) = pair_1_list(i_pair_1,3) + 1
		            endif
		        else
		        	if(ovlp3fn(i_atom_l1)%n(n)*max(ovl_fact11(i_atom_l1,i_cell_1),ovl_fact22(i_atom_l1)) > crit_val .or. &
		               (calc_deriv .and. (ovlp3fn(i_atom_l1)%n(n)*max(ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                           ovl_fact22_d(i_atom_l1)*max_abs_dm_row_tot) > crit_val_deriv)) .or. &
		               (AS_stress_on .and. (ovlp3fn(i_atom_l1)%n(n)*max(AS_ovl_fact11_d(i_atom_l1,i_cell_1), &
		                                           AS_ovl_fact22_d(i_atom_l1)*max_abs_dm_row_tot) > AS_crit_val_deriv))) then
		              pair_1_list(i_pair_1,3) = pair_1_list(i_pair_1,3) + 1
		            endif
		        end if
                n = n + 1
              enddo
            endif

            ! Number of basis pairs in OVL2
            if(my_pair_flag_bvk(i_atom_r1,i_atom_l1,inv_cell_bvk(i_cell_1))) then
              n = pair_offset(i_atom_r1,i_atom_l1,inv_cell_bvk(i_cell_1)) + i_basis_r1-nr1_off
              do il1 = 1, nl1_len
              	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            if(ovlp3fn_SR(i_atom_r1)%n(n)*ovl_fact12(i_atom_l1,i_cell_1) > crit_val .or. &
		               (calc_deriv .and. (ovlp3fn_SR(i_atom_r1)%n(n)*ovl_fact12_d(i_atom_l1,i_cell_1) > crit_val_deriv)) .or. &
		               (AS_stress_on .and. (ovlp3fn_SR(i_atom_r1)%n(n)*AS_ovl_fact12_d(i_atom_l1,i_cell_1) > AS_crit_val_deriv))) then
		              pair_1_list(i_pair_1,4) = pair_1_list(i_pair_1,4) + 1
		            endif
		        else
		        	if(ovlp3fn(i_atom_r1)%n(n)*ovl_fact12(i_atom_l1,i_cell_1) > crit_val .or. &
		               (calc_deriv .and. (ovlp3fn(i_atom_r1)%n(n)*ovl_fact12_d(i_atom_l1,i_cell_1) > crit_val_deriv)) .or. &
		               (AS_stress_on .and. (ovlp3fn(i_atom_r1)%n(n)*AS_ovl_fact12_d(i_atom_l1,i_cell_1) > AS_crit_val_deriv))) then
		              pair_1_list(i_pair_1,4) = pair_1_list(i_pair_1,4) + 1
		            endif
		        end if
                n = n + nr1_len
              enddo
            endif
          endif

        enddo
      enddo

! End timing block,   Other #6
times_other(6) = times_other(6) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #4
      call time_sync_start
      call sync_integer_vector(pair_1_list(1,3), n_pair_1)
      call sync_integer_vector(pair_1_list(1,4), n_pair_1)
      call time_sync_end
! End timing block,   Sync/Imbalance #4
ttt_other = mpi_wtime()
! Begin timing block, Other #7

      ! Get offsets and total number of basis pairs

      n_ovlp3fn_rcv_1 = 0
      n_ovlp3fn_rcv_2 = 0
      do i_pair_1 = 1, n_pair_1
        pair_1_list(i_pair_1,5) = n_ovlp3fn_rcv_1 ! Offset in received ovlp3fn
        pair_1_list(i_pair_1,6) = n_ovlp3fn_rcv_2 ! Offset in received ovlp3fn
        n_ovlp3fn_rcv_1 = n_ovlp3fn_rcv_1 + pair_1_list(i_pair_1,3)
        n_ovlp3fn_rcv_2 = n_ovlp3fn_rcv_2 + pair_1_list(i_pair_1,4)
      enddo

      ! Gather and broadcast the needed ovlp3fn entries to all
      ! Since a number of irregular broadcasts may get expensive,
      ! every task puts its own entries to ovlp3fn_rcv
      ! and then a sync_vector is done.

      ! The contributions for every atom pair are sorted with descending norm,
      ! only the number of basis pairs calculated above is entered into ovlp3fn_rcv.

      i_pair_1 = 0
      ovlp3fn_rcv(:,:,:) = 0.     ! Basis pairs
      norm_ovlp3fn_rcv(:,:) = 0.  ! Norm of basis pair
      idx_ovlp3fn_rcv(:,:) = 0    ! Left basis number for basis pair

      do i_cell_1 = 1, n_cells_bvk
        do i_atom_l1 = 1, n_atoms2

          nl1_off = atom2basis_off2(i_atom_l1)
          nl1_len = atom2basis_len2(i_atom_l1)

          if(pair_flag_bvk(i_atom_l1,i_atom_r1,i_cell_1)) then

            i_pair_1 = i_pair_1 + 1

            if(my_pair_flag_bvk(i_atom_l1,i_atom_r1,i_cell_1)) then
              n_off = pair_offset(i_atom_l1,i_atom_r1,i_cell_1) + (i_basis_r1-nr1_off-1)*nl1_len
              do i=1,nl1_len
              	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
                	aux(1,i) = -ovlp3fn_SR(i_atom_l1)%n(n_off+i)
                else
                	aux(1,i) = -ovlp3fn(i_atom_l1)%n(n_off+i)
                end if
                aux(2,i) = i
              enddo
              call heapsort_general(aux, 2, nl1_len, 1)
              j = pair_1_list(i_pair_1,5)
              nbb = atom2basbas_len(i_atom_l1)
              do i=1,pair_1_list(i_pair_1,3)
                j = j+1
                il1 = int(aux(2,i))
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            ovlp3fn_rcv(1:nbb,j,1) = ovlp3fn_SR(i_atom_l1)%m(1:nbb,n_off+il1)
		            norm_ovlp3fn_rcv (j,1) = ovlp3fn_SR(i_atom_l1)%n(n_off+il1)
                else
                	ovlp3fn_rcv(1:nbb,j,1) = ovlp3fn(i_atom_l1)%m(1:nbb,n_off+il1)
		            norm_ovlp3fn_rcv (j,1) = ovlp3fn(i_atom_l1)%n(n_off+il1)
                end if
                idx_ovlp3fn_rcv  (j,1) = nl1_off + il1
              enddo
            endif
            if(my_pair_flag_bvk(i_atom_r1,i_atom_l1,inv_cell_bvk(i_cell_1))) then
              n_off = pair_offset(i_atom_r1,i_atom_l1,inv_cell_bvk(i_cell_1)) + i_basis_r1-nr1_off-1
              do i=1,nl1_len
              	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
                	aux(1,i) = -ovlp3fn_SR(i_atom_r1)%n(n_off+(i-1)*nr1_len+1)
                else
                	aux(1,i) = -ovlp3fn(i_atom_r1)%n(n_off+(i-1)*nr1_len+1)
                end if
                aux(2,i) = i
              enddo
              call heapsort_general(aux, 2, nl1_len, 1)
              j = pair_1_list(i_pair_1,6)
              nbb = atom2basbas_len(i_atom_r1)
              do i=1,pair_1_list(i_pair_1,4)
                j = j+1
                il1 = int(aux(2,i))
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            ovlp3fn_rcv(1:nbb,j,2) = ovlp3fn_SR(i_atom_r1)%m(1:nbb,n_off+(il1-1)*nr1_len+1)
		            norm_ovlp3fn_rcv (j,2) = ovlp3fn_SR(i_atom_r1)%n(n_off+(il1-1)*nr1_len+1)
                else
                	ovlp3fn_rcv(1:nbb,j,2) = ovlp3fn(i_atom_r1)%m(1:nbb,n_off+(il1-1)*nr1_len+1)
		            norm_ovlp3fn_rcv (j,2) = ovlp3fn(i_atom_r1)%n(n_off+(il1-1)*nr1_len+1)
                end if
                idx_ovlp3fn_rcv  (j,2) = nl1_off + il1
              enddo
            endif
          endif

        enddo
      enddo

! End timing block,   Other #7
times_other(7) = times_other(7) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #5
      call time_sync_start
      call sync_vector(ovlp3fn_rcv(1,1,1),size(ovlp3fn_rcv,1)*n_ovlp3fn_rcv_1)
      call sync_vector(ovlp3fn_rcv(1,1,2),size(ovlp3fn_rcv,1)*n_ovlp3fn_rcv_2)
      call sync_vector(norm_ovlp3fn_rcv(1,1),n_ovlp3fn_rcv_1)
      call sync_vector(norm_ovlp3fn_rcv(1,2),n_ovlp3fn_rcv_2)
      call sync_integer_vector(idx_ovlp3fn_rcv(1,1),n_ovlp3fn_rcv_1)
      call sync_integer_vector(idx_ovlp3fn_rcv(1,2),n_ovlp3fn_rcv_2)
      call time_sync_end
! End timing block,   Sync/Imbalance #5
ttt_other = mpi_wtime()
! Begin timing block, Other #8


      ! Get the maximum norm of basis pairs per atom pair.
      ! Since basis pairs are sorted in descending order, we just need the first norm.

      do i_pair_1 = 1, n_pair_1

        no3fn_len = pair_1_list(i_pair_1,3)
        no3fn_off = pair_1_list(i_pair_1,5)
        if(no3fn_len==0) then
          max_norm_ovlp3fn_rcv(i_pair_1,1) = 0.
        else
          max_norm_ovlp3fn_rcv(i_pair_1,1) = norm_ovlp3fn_rcv(no3fn_off+1,1)
        endif

        no3fn_len = pair_1_list(i_pair_1,4)
        no3fn_off = pair_1_list(i_pair_1,6)
        if(no3fn_len==0) then
          max_norm_ovlp3fn_rcv(i_pair_1,2) = 0.
        else
          max_norm_ovlp3fn_rcv(i_pair_1,2) = norm_ovlp3fn_rcv(no3fn_off+1,2)
        endif

      enddo


      fock_matrix_row(:,:,:) = 0

! End timing block,   Other #8
times_other(8) = times_other(8) + mpi_wtime() - ttt_other

      do i_my_atom = 1, my_n_atoms
! End timing block,   Other #9
ttt_other = mpi_wtime()

        i_atom_l2 = my_atom_list(i_my_atom)

        nl2_bb_off = atom2basbas_off(vb_atom2atom(i_atom_l2))
        nl2_bb_len = atom2basbas_len(i_atom_l2)
        nl2_off    = atom2basis_off2(i_atom_l2)
        nl2_len    = atom2basis_len2(i_atom_l2)

        ! tmp1/tmp2 must be allocated to the exact size

        call aims_allocate(tmp1,nl2_bb_len,nl2_len,n_cells_bvk,n_spin,'tmp1')
        call aims_allocate(tmp2,nl2_bb_len,nl2_len,n_cells_bvk,n_spin,'tmp2')
        nbb = atom2basbas_len(i_atom_r1)
        mem_tmpx = int(nbb,kind=8) * int(nl2_len,kind=8) * &
             int(n_cells_bvk,kind=8) * int(n_spin,kind=8) * int(8,kind=8)
        allocate(tmpx(nbb, nl2_len, n_cells_bvk, n_spin), stat=info)
        call update_when_allocating(info, "tmpx", mem_tmpx)

        ! Get maximum norm of the factors with which tmp1 (calculated below) will be multiplied.
        ! This serves as a criterion, which entries have to be taken when calculating tmp1 below.

        max_fact_tmp1(:,:) = 0
        AS_max_fact_tmp1(:) = 0.0d0

        do i_cell = 1, n_cells_bvk
          do i_pair_1 = 1, n_pair_1

            i_atom_l1 = pair_1_list(i_pair_1,1)
            i_cell_1  = pair_1_list(i_pair_1,2)

            i_cell_p1 = add_cells(i_cell,i_cell_1)

            ! Factor for fock_matrix
            if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
	            fact = coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_ovlp3fn_rcv(i_pair_1,1)
	        else
	        	fact = coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_ovlp3fn_rcv(i_pair_1,1)
	        end if
            max_fact_tmp1(i_cell,1) = max(max_fact_tmp1(i_cell,1), fact)

            if(calc_deriv) then
              ! Factors for d_exx_ene
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          fact = max(d_coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_ovlp3fn_rcv(i_pair_1,1), &
		                     d_coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell   )*max_norm_ovlp3fn_rcv(i_pair_1,2))
              else
		          fact = max(d_coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_ovlp3fn_rcv(i_pair_1,1), &
				                 d_coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell   )*max_norm_ovlp3fn_rcv(i_pair_1,2))
              end if
              fact = fact*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)
              max_fact_tmp1(i_cell,2) = max(max_fact_tmp1(i_cell,2), fact)
            endif

            if(AS_stress_on) then
              ! Factors for analytical stress
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          fact = max(AS_d_coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_ovlp3fn_rcv(i_pair_1,1), &
		                     AS_d_coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell   )*max_norm_ovlp3fn_rcv(i_pair_1,2))
              else
		          fact = max(AS_d_coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell_p1)*max_norm_ovlp3fn_rcv(i_pair_1,1), &
		                     AS_d_coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell   )*max_norm_ovlp3fn_rcv(i_pair_1,2))
              end if
              fact = fact*max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)
              AS_max_fact_tmp1(i_cell) = max(AS_max_fact_tmp1(i_cell), fact)
            endif
          enddo
        enddo

        ! TMP1(1:l2b,ic) = SUM over r2 ( OVL1(1:l2b,r2,ic2)*DM_ROW(r2,ic+ic2) )
        ! ----------------------------------------------------------------------

        cell_used(:,:) = .false.

        if (.not.(use_mpi_in_place)) then
           cell_used_rcv(:,:) = .false.
        endif

        tmp1 = 0

! End timing block,   Other #9
ttt0 = mpi_wtime()
times_other(9) = times_other(9) + ttt0 - ttt_other
! Begin timing block, DM_x_o3fn
        do i_atom_pair = my_atom_pair_s(i_my_atom), my_atom_pair_e(i_my_atom)

          if(i_atom_l2 /= pair_list(1,i_atom_pair)) call aims_stop('Inconsistent pair_list')
          i_atom_r2 = pair_list(2,i_atom_pair)
          i_cell_2  = pair_list(3,i_atom_pair)

          nr2_off    = atom2basis_off2(i_atom_r2)
          nr2_len    = atom2basis_len2(i_atom_r2)

          do i_spin = 1, n_spin

            ! Find those cells of the density matrix for the actual (i_atom_r2,i_basis_r1)
            ! which are relevant for calculations and store them in dm_aux

            n_cells_used = 0
            do i_cell = 1, n_cells_bvk

              ! Density matrix cell used for dm_x_o3fn1:
              i_cell_dm = add_cells(i_cell,i_cell_2)
			  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
              	max_norm_dxo = max_abs_dm_row(i_atom_r2,i_cell_dm)*max_norm_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2)
              else
                max_norm_dxo = max_abs_dm_row(i_atom_r2,i_cell_dm)*max_norm_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2)
              end if
              if(max_norm_dxo*max_fact_tmp1(i_cell,1) > crit_val .or. &
                 max_norm_dxo*max_fact_tmp1(i_cell,2) > crit_val_deriv .or. &
                 max_norm_dxo*AS_max_fact_tmp1(i_cell) > AS_crit_val_deriv) then
                n_cells_used = n_cells_used + 1
                idx_cells(n_cells_used) = i_cell
                dm_aux(1:nr2_len,n_cells_used) = dm_row(nr2_off+1:nr2_off+nr2_len,i_cell_dm,i_spin)
                cell_used(i_cell,i_spin) = .true.
              endif
            enddo

            ! Multiply o3fn for the current i_atom_pair with dm_aux (i.e. the parts of the density matrix surviving screening)

            do i_basis_l2 = nl2_off+1, nl2_off+nl2_len
              i_pair_2 = pair_offset(i_atom_l2,i_atom_r2,i_cell_2) + i_basis_l2-nl2_off
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          call dgemm('N','N',nl2_bb_len,n_cells_used,nr2_len,1.d0, &
		                  ovlp3fn_SR(i_atom_l2)%m(1,i_pair_2), nl2_bb_len*nl2_len, &
		                  dm_aux,size(dm_aux,1), &
		                  0.d0,dm_x_o3fn_aux,size(dm_x_o3fn_aux,1))
                 else
		          call dgemm('N','N',nl2_bb_len,n_cells_used,nr2_len,1.d0, &
		                  ovlp3fn(i_atom_l2)%m(1,i_pair_2), nl2_bb_len*nl2_len, &
		                  dm_aux,size(dm_aux,1), &
		                  0.d0,dm_x_o3fn_aux,size(dm_x_o3fn_aux,1))
              end if
              do i_cell = 1, n_cells_used
                tmp1(1:nl2_bb_len,i_basis_l2-nl2_off,idx_cells(i_cell),i_spin) = &
                tmp1(1:nl2_bb_len,i_basis_l2-nl2_off,idx_cells(i_cell),i_spin) + &
                  dm_x_o3fn_aux(1:nl2_bb_len,i_cell)
              enddo
            enddo
          enddo
        enddo
! End timing block,   DM_x_o3fn
times(1) = times(1) + mpi_wtime()-ttt0
! Begin timing block, Sync/Imbalance #6
        if(my_tasks_per_atom > 1) then
          call time_sync_start(mpi_comm_atom)
          ! NB: sync_logical_vector doesn't support a communicator, so we nust use mpi_allreduce
          if (.not.(use_mpi_in_place)) then
             call mpi_allreduce(cell_used,cell_used_rcv,n_spin*n_cells_bvk,MPI_LOGICAL,MPI_LOR,mpi_comm_atom,mpierr)
             cell_used = cell_used_rcv
          else
             call mpi_allreduce(mpi_in_place,cell_used,n_spin*n_cells_bvk,MPI_LOGICAL,MPI_LOR,mpi_comm_atom,mpierr)
          endif

          call sync_partial(tmp1, cell_used)
          call time_sync_end
        endif

!SK
!        write(info_str,*) tmp1(1,:,:,:)
!        write(unit=use_unit,fmt="(A,I8,/,A)") 'Basis: ',i_basis_r1,trim(info_str)
! End timing block,   Sync/Imbalance #6
! Begin timing block, Other #10
ttt_other = mpi_wtime()

        ! Get max norm of dm*ovlp3fn

        max_norm_tmp(:) = 0
        do i_cell = 1, n_cells_bvk
          do i_spin = 1, n_spin
            if(cell_used(i_cell,i_spin)) &
              max_norm_tmp(i_cell) = max(max_norm_tmp(i_cell),max_column_norm(tmp1(:,:,i_cell,i_spin)))
          enddo
        enddo

        ! Get maximum norm of the factors with which tmp2 (calculated below) will be multiplied.

		if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
        	max_fact1 = max_norm_ovlp3fn_per_latom_SR(i_atom_l2)
        else
        	max_fact1 = max_norm_ovlp3fn_per_latom(i_atom_l2)
        end if
        max_fact2 = 0
        AS_max_fact3 = 0.0d0

        if (.not.(use_mpi_in_place)) then
           max_fact2_rcv = 0
           AS_max_fact3_rcv = 0.0d0
        endif

        if(calc_deriv .or. AS_stress_on) then
          do i_atom_pair = my_atom_pair_s(i_my_atom), my_atom_pair_e(i_my_atom)
            i_atom_r2 = pair_list(2,i_atom_pair)
            i_cell_2  = pair_list(3,i_atom_pair)
            do i_cell = 1, n_cells_bvk
              i_cell_fock = add_cells(i_cell,i_cell_2)
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
              	max_fact2 = max(max_fact2, max_norm_d_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2)*max_abs_dm_row(i_atom_r2,i_cell_fock))
              else
              	max_fact2 = max(max_fact2, max_norm_d_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2)*max_abs_dm_row(i_atom_r2,i_cell_fock))
              end if
              if(AS_stress_on) then
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            AS_max_fact3 = max(AS_max_fact3, &
		              AS_max_norm_d_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2)*max_abs_dm_row(i_atom_r2,i_cell_fock))
                else
		            AS_max_fact3 = max(AS_max_fact3, &
		              AS_max_norm_d_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2)*max_abs_dm_row(i_atom_r2,i_cell_fock))
                end if
              end if
            enddo
          enddo
          ! max_fact must be identical for all tasks working on one atom
! End timing block,   Other #10 (if fork taken)
times_other(10) = times_other(10) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #7
          if(my_tasks_per_atom > 1 .and. calc_deriv) then
            call time_sync_start(mpi_comm_atom)

            if (.not.(use_mpi_in_place)) then
               call mpi_allreduce(max_fact2,max_fact2_rcv,1,mpi_real8,mpi_max,mpi_comm_atom,mpierr)
               max_fact2 = max_fact2_rcv
               if(AS_stress_on) then
                 call mpi_allreduce(AS_max_fact3,AS_max_fact3_rcv,1,mpi_real8,mpi_max,mpi_comm_atom,mpierr)
                 AS_max_fact3 = AS_max_fact3_rcv
               endif
            else
               call mpi_allreduce(mpi_in_place,max_fact2,1,mpi_real8,mpi_max,mpi_comm_atom,mpierr)
               if(AS_stress_on) then
                 call mpi_allreduce(mpi_in_place,AS_max_fact3,1,mpi_real8,mpi_max,mpi_comm_atom,mpierr)
               endif
            endif

            call time_sync_end
          endif
! End timing block,   Sync/Imbalance #7
        else
! End timing block,   Other #10 (if fork not taken)
times_other(10) = times_other(10) + mpi_wtime() - ttt_other
        endif
! Begin timing block, Other #11
ttt_other = mpi_wtime()

        ! Calculate HF11, TMP2, TMPX
        ! See the above pseudocode for details!

        if(my_tasks_per_atom > 1) then
          n_bb_s = my_cm_bb_s
          n_bb_e = my_cm_bb_e
        else
          n_bb_s = 1
          n_bb_e = atom2basbas_len(i_atom_l2)
        endif

        if(n_bb_s > n_bb_e) call aims_stop('Illegal settings for my_cm_bb_s/e')

        nl2_off = my_basis_off(i_atom_l2) ! Offset in fock_matrix and dm_cols, nl2_off is only needed for that
        nl2_len = atom2basis_len2(i_atom_l2)
        if(nl2_off<0) call aims_stop('Illegal setting in my_basis_off')

        tmp2 = 0
        cell_used(:,:) = .false.

        tmpx = 0
        cell_tmpx_used(:,:) = .false.
        mult_cm_r1(:) = .false.
        mult_d_cm_r1(:) = .false.
        AS_mult_d_cm_r1(:) = .false.

        work_count = 0

! End timing block,   Other #11
ttt0 = mpi_wtime()
times_other(11) = times_other(11) + ttt0 - ttt_other
! Begin timing block, Products #1
        do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc

          do i_pair_1 = 1, n_pair_1

            i_atom_l1 = pair_1_list(i_pair_1,1)
            i_cell_1  = pair_1_list(i_pair_1,2)

            i_cell_m1 = sub_cells(i_cell,i_cell_1)
            i_cell_p1 = add_cells(i_cell,i_cell_1)

            nl1_off = atom2basis_off2(i_atom_l1)
            nl1_len = atom2basis_len2(i_atom_l1)

            max_abs_dm    = max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell)
            max_abs_dm_p1 = max_abs_dm_cols(i_atom_l1,i_atom_l2,i_cell_p1)

            ! For dealing with ovlp3fn_rcv(:,:,1):
            no3fn_len = pair_1_list(i_pair_1,3)
            no3fn_off = pair_1_list(i_pair_1,5)

			if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
            	cm_norm = coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)
            else
            	cm_norm = coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)
            end if

            ! Get the number of columns of ovlp3fn_rcv needed for HF11
            if(i_cell <= inv_cell_bvk(i_cell)) then
              do i=1,no3fn_len
                if(norm_ovlp3fn_rcv(no3fn_off+i,1)*cm_norm*max_norm_tmp(i_cell_m1) <= crit_val) exit
              enddo
              n_cols1 = i-1
            else
              n_cols1 = 0
            endif

            ! Get the number of columns of ovlp3fn_rcv needed for HF22 and gradient (part 1)
            do i=1,no3fn_len
              if(norm_ovlp3fn_rcv(no3fn_off+i,1)*cm_norm*max_abs_dm*max_fact1 <= crit_val .and. &
                 norm_ovlp3fn_rcv(no3fn_off+i,1)*cm_norm*max_abs_dm*max_fact2 <= crit_val_deriv .and. &
                 norm_ovlp3fn_rcv(no3fn_off+i,1)*cm_norm*max_abs_dm*AS_max_fact3 <= AS_crit_val_deriv) exit
            enddo
            n_cols2 = i-1

            n_cols = max(n_cols1,n_cols2)

            if(n_cols > 0) then

              ! TMPP(1:l1b) = OVL_RCV1(1:l1b,ic1) * CM(l1,l2,ic)
              ! ------------------------------------------------
              tmp_prod(1:n_bb_e-n_bb_s+1,1:n_cols) = 0
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          call mult_coul_mat_left(n_cols, coul_mat_store_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell), &
		                                  ovlp3fn_rcv(1,no3fn_off+1,1), size(ovlp3fn_rcv,1), &
		                                  tmp_prod, size(tmp_prod,1))
              else
		          call mult_coul_mat_left(n_cols, coul_mat_store(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell), &
		                                  ovlp3fn_rcv(1,no3fn_off+1,1), size(ovlp3fn_rcv,1), &
		                                  tmp_prod, size(tmp_prod,1))
              end if

              if(n_cols1 > 0) then
                ! We make use of the fact that we need to calculate only one of two
                ! symmetric cells - although it is questionable if that buys us much
                if(i_cell == inv_cell_bvk(i_cell)) then
                  s = 1.d0
                else
                  s = 2.d0
                endif

                ! HF11(1:l1b,1:l2b,ic) = HF11(1:l1b,1:l2b,ic) + TMPP(1:l1b) * TMP1(1:l2b,ic-ic1)
                ! ------------------------------------------------------------------------------
                do i_spin = 1, n_spin
                  call dgemm('T','N',n_cols1,nl2_len,n_bb_e-n_bb_s+1,s, &
                           tmp_prod, size(tmp_prod,1), &
                           tmp1(n_bb_s,1,i_cell_m1,i_spin), size(tmp1,1), &
                           0.d0,tmp_atom,size(tmp_atom,1))
                  do i = 1, nl2_len
                  do j = 1, n_cols1
                  	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		                fock_matrix_SR(idx_ovlp3fn_rcv(no3fn_off+j,1),nl2_off+i,i_cell,i_spin) = &
		                fock_matrix_SR(idx_ovlp3fn_rcv(no3fn_off+j,1),nl2_off+i,i_cell,i_spin) + tmp_atom(j,i)
                    else
		                fock_matrix(idx_ovlp3fn_rcv(no3fn_off+j,1),nl2_off+i,i_cell,i_spin) = &
		                fock_matrix(idx_ovlp3fn_rcv(no3fn_off+j,1),nl2_off+i,i_cell,i_spin) + tmp_atom(j,i)
                    end if
                  enddo
                  enddo
                enddo
              endif

              if(n_cols2 > 0) then
                ! TMP2(1:l2b,ic-ic1) = TMP2(1:l2b,ic-ic1) + TMPP(1:l1b) * DM(1:l1b,1:l2b,ic)
                ! --------------------------------------------------------------------------
                do i_spin = 1, n_spin
                  do i = 1, nl2_len
                  do j = 1, n_cols2
                    tmp_atom(j,i) = dm_cols(idx_ovlp3fn_rcv(no3fn_off+j,1),nl2_off+i,i_cell,i_spin)
                  enddo
                  enddo
                  call dgemm('N','N',n_bb_e-n_bb_s+1,nl2_len,n_cols2,1.d0, &
                           tmp_prod,size(tmp_prod,1),tmp_atom,size(tmp_atom,1), &
                           1.d0,tmp2(n_bb_s,1,i_cell_m1,i_spin),size(tmp2,1))
                  cell_used(i_cell_m1,i_spin) = .true.
                enddo
              endif

            endif

            ! Get the number of columns of ovlp3fn_rcv needed for gradient (part 2)
            if(calc_deriv .and. (i_atom_l2 /= i_atom_l1)) then
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
              	cm_norm = d_coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)
              else
              	cm_norm = d_coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)
              end if
              do i=1,no3fn_len
                if(norm_ovlp3fn_rcv(no3fn_off+i,1)*cm_norm*max_norm_tmp(i_cell_m1)*max_abs_dm <= crit_val_deriv) exit
              enddo
              n_cols = i-1
            else
              n_cols = 0
            endif

            if(n_cols > 0) then
              do i_grad = 1, 3
                tmp_prod(1:n_bb_e-n_bb_s+1,1:n_cols) = 0
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            call mult_coul_mat_left(n_cols, d_coul_mat_store_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell,i_grad), &
		                                    ovlp3fn_rcv(1,no3fn_off+1,1), size(ovlp3fn_rcv,1), &
		                                    tmp_prod, size(tmp_prod,1))
                else
		            call mult_coul_mat_left(n_cols, d_coul_mat_store(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell,i_grad), &
		                                    ovlp3fn_rcv(1,no3fn_off+1,1), size(ovlp3fn_rcv,1), &
		                                    tmp_prod, size(tmp_prod,1))
                end if
                do i_spin = 1, n_spin
                  call dgemm('T','N',n_cols,nl2_len,n_bb_e-n_bb_s+1,1.d0, &
                           tmp_prod, size(tmp_prod,1), &
                           tmp1(n_bb_s,1,i_cell_m1,i_spin), size(tmp1,1), &
                           0.d0,tmp_atom,size(tmp_atom,1))
                  s = sum(tmp_atom(1:n_cols,1:nl2_len)*  &
                          dm_cols(idx_ovlp3fn_rcv(no3fn_off+1:no3fn_off+n_cols,1), &
                                  nl2_off+1:nl2_off+nl2_len,i_cell,i_spin))
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		               d_exx_ene_SR(i_grad,i_atom_l2) = d_exx_ene_SR(i_grad,i_atom_l2) + 2*s
		               d_exx_ene_SR(i_grad,i_atom_l1) = d_exx_ene_SR(i_grad,i_atom_l1) - 2*s
                  else
		               d_exx_ene(i_grad,vb_atom2atom(i_atom_l2)) = d_exx_ene(i_grad,vb_atom2atom(i_atom_l2)) + 2*s
		               d_exx_ene(i_grad,vb_atom2atom(i_atom_l1)) = d_exx_ene(i_grad,vb_atom2atom(i_atom_l1)) - 2*s
                  end if
                enddo
              enddo
            endif

            if(AS_stress_on) then
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
             	AS_cm_norm = AS_d_coul_mat_norm_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)
              else
              	AS_cm_norm = AS_d_coul_mat_norm(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell)
              end if
              do i=1,no3fn_len
                if(norm_ovlp3fn_rcv(no3fn_off+i,1)*AS_cm_norm*max_norm_tmp(i_cell_m1)*max_abs_dm <= AS_crit_val_deriv) exit
              enddo
              AS_n_cols = i-1
            else
              AS_n_cols = 0
            endif

            if(AS_n_cols > 0) then
              do AS_index = 1, AS_components, 1
                tmp_prod(1:n_bb_e-n_bb_s+1,1:AS_n_cols) = 0
                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		            call mult_coul_mat_left(AS_n_cols, AS_d_coul_mat_store_SR(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell,AS_index), &
		                                    ovlp3fn_rcv(1,no3fn_off+1,1), size(ovlp3fn_rcv,1), &
		                                    tmp_prod, size(tmp_prod,1))
                else
		            call mult_coul_mat_left(AS_n_cols, AS_d_coul_mat_store(vb_atom2atom(i_atom_l1),vb_atom2atom(i_atom_l2),i_cell,AS_index), &
		                                    ovlp3fn_rcv(1,no3fn_off+1,1), size(ovlp3fn_rcv,1), &
		                                    tmp_prod, size(tmp_prod,1))
                end if
                do i_spin = 1, n_spin
                  call dgemm('T','N',AS_n_cols,nl2_len,n_bb_e-n_bb_s+1,1.d0, &
                          tmp_prod, size(tmp_prod,1), &
                          tmp1(n_bb_s,1,i_cell_m1,i_spin), size(tmp1,1), &
                          0.d0,tmp_atom,size(tmp_atom,1))
                  s = sum(tmp_atom(1:AS_n_cols,1:nl2_len)*  &
                          dm_cols(idx_ovlp3fn_rcv(no3fn_off+1:no3fn_off+AS_n_cols,1), &
                                  nl2_off+1:nl2_off+nl2_len,i_cell,i_spin))
                  AS_EXX_stress_local(AS_index) = AS_EXX_stress_local(AS_index) + 2*s
                enddo
              enddo
            endif

            ! For dealing with ovlp3fn_rcv(:,:,2):
            no3fn_len = pair_1_list(i_pair_1,4)
            no3fn_off = pair_1_list(i_pair_1,6)

			if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
            	cm_norm = coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)
            else
            	cm_norm = coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)
            end if

            ! Get the number of columns of ovlp3fn_rcv needed for HF12
            do i=1,no3fn_len
              if(norm_ovlp3fn_rcv(no3fn_off+i,2)*cm_norm*max_abs_dm_p1*max_fact1 <= crit_val .and. &
                 norm_ovlp3fn_rcv(no3fn_off+i,2)*cm_norm*max_abs_dm_p1*max_fact2 <= crit_val_deriv .and. &
                 norm_ovlp3fn_rcv(no3fn_off+i,2)*cm_norm*max_abs_dm_p1*AS_max_fact3 <= AS_crit_val_deriv) exit
            enddo
            n_cols1 = i-1

            ! Get the number of columns of ovlp3fn_rcv needed for gradient (part 3)
            if(calc_deriv) then
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
              	cm_norm = d_coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)
              else
              	cm_norm = d_coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)
              end if
              do i=1,no3fn_len
                if(norm_ovlp3fn_rcv(no3fn_off+i,2)*cm_norm*max_abs_dm_p1*max_norm_tmp(i_cell) <= crit_val_deriv) exit
              enddo
              n_cols2 = i-1
            else
              n_cols2 = 0
            endif

            ! Get the number of columns of ovlp3fn_rcv needed for stress (part 3)
            if(AS_stress_on) then
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
              	AS_cm_norm = AS_d_coul_mat_norm_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)
              else
              	AS_cm_norm = AS_d_coul_mat_norm(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell)
              end if
              do i=1,no3fn_len
                if(norm_ovlp3fn_rcv(no3fn_off+i,2)*AS_cm_norm*max_abs_dm_p1*max_norm_tmp(i_cell) <= AS_crit_val_deriv) exit
              enddo
              AS_n_cols = i-1
            else
              AS_n_cols = 0
            endif

            n_cols = max(n_cols1, n_cols2)

            if(n_cols > 0) then
              do i_spin = 1, n_spin
                ! Please note: If the Coulomb matrix cells are not subdivided when my_tasks_per_atom>1,
                ! there would be done redundant work here - all tasks of one atom do the same.
                ! Therefore only one processor does the work here (and we have an additional sync below)
                work_count = work_count+1
                cell_tmpx_used(i_cell,i_spin) = .true.
                if(my_tasks_per_atom>1 .and. my_cm_cell_inc==1 .and. mod(work_count-1,my_tasks_per_atom)/=my_atom_id) cycle

                ! TMPX(1:l2b,ic-ic1) = TMPX(1:l2b,ic-ic1) + OVL_RCV2(1:l1b,ic1)*DM(1:l1b,1:l2b,ic)
                ! --------------------------------------------------------------------------------
                do i = 1, nl2_len
                do j = 1, n_cols
                  tmp_atom(j,i) = dm_cols(idx_ovlp3fn_rcv(no3fn_off+j,2),nl2_off+i,i_cell_p1,i_spin)
                enddo
                enddo
                nbb = atom2basbas_len(i_atom_r1)
                call dgemm('N','N',nbb,nl2_len,n_cols,1.d0, &
                        ovlp3fn_rcv(1,no3fn_off+1,2), size(ovlp3fn_rcv,1), &
                        tmp_atom, size(tmp_atom,1), &
                        1.d0,tmpx(1,1,i_cell,i_spin),size(tmpx,1))
              enddo
              if(n_cols1 > 0)   mult_cm_r1(i_cell)      = .true.
              if(n_cols2 > 0)   mult_d_cm_r1(i_cell)    = .true.
              if(AS_n_cols > 0) AS_mult_d_cm_r1(i_cell) = .true.
            endif
          enddo

        enddo

        if(my_tasks_per_atom>1 .and. my_cm_cell_inc==1) then
          ! We have to sync, see above!
! End timing block,   Products #1
times(3) = times(3) + mpi_wtime()-ttt0
! Begin timing block, Sync/Imbalance #8
          call time_sync_start(mpi_comm_atom)
          ! NB: cell_tmpx_used is already in sync here
          call sync_partial(tmpx, cell_tmpx_used)

          call time_sync_end
! End timing block,   Sync/Imbalance #8
ttt0 = mpi_wtime()
! Begin timing block, Products #2
        endif

!SK
!        write(info_str,*) tmpx(1,:,:,:)
!        write(unit=use_unit,fmt="(A,I8,/,A)") 'tmpxBasis: ',i_basis_r1,trim(info_str)

        ! We have to do the gradients first since we want to re-use tmp1 below
        if(calc_deriv .and. (i_atom_l2 /= i_atom_r1)) then
          do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc
            if(mult_d_cm_r1(i_cell)) then
              do i_spin = 1, n_spin
                do i_grad = 1,3
                  tmp_prod(1:n_bb_e-n_bb_s+1,1:nl2_len) = 0
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		              call mult_coul_mat_left(nl2_len, d_coul_mat_store_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell,i_grad), &
		                                      tmpx(1,1,i_cell,i_spin),size(tmpx,1), &
		                                      tmp_prod, size(tmp_prod,1))
                  else
		              call mult_coul_mat_left(nl2_len, d_coul_mat_store(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell,i_grad), &
		                                      tmpx(1,1,i_cell,i_spin),size(tmpx,1), &
		                                      tmp_prod, size(tmp_prod,1))
                  end if
                  do i = 1,nl2_len
                    s = dot_product(tmp_prod(1:n_bb_e-n_bb_s+1,i),tmp1(n_bb_s:n_bb_e,i,i_cell,i_spin))
                    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		                 d_exx_ene_SR(i_grad,i_atom_l2) = d_exx_ene_SR(i_grad,i_atom_l2) + 2*s
		                 d_exx_ene_SR(i_grad,i_atom_r1) = d_exx_ene_SR(i_grad,i_atom_r1) - 2*s
                    else
		                 d_exx_ene(i_grad,vb_atom2atom(i_atom_l2)) = d_exx_ene(i_grad,vb_atom2atom(i_atom_l2)) + 2*s
		                 d_exx_ene(i_grad,vb_atom2atom(i_atom_r1)) = d_exx_ene(i_grad,vb_atom2atom(i_atom_r1)) - 2*s
                    end if
                  enddo
                enddo
              enddo
            endif
          enddo
        endif

        if(AS_stress_on) then
          do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc
            if(AS_mult_d_cm_r1(i_cell)) then
              do i_spin = 1, n_spin
                do AS_index = 1, AS_components, 1
                  tmp_prod(1:n_bb_e-n_bb_s+1,1:nl2_len) = 0
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		              call mult_coul_mat_left(nl2_len, AS_d_coul_mat_store_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell,AS_index), &
		                                      tmpx(1,1,i_cell,i_spin),size(tmpx,1), &
		                                      tmp_prod, size(tmp_prod,1))
                  else
		              call mult_coul_mat_left(nl2_len, AS_d_coul_mat_store(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell,AS_index), &
		                                      tmpx(1,1,i_cell,i_spin),size(tmpx,1), &
		                                      tmp_prod, size(tmp_prod,1))
                  end if
                  do i = 1,nl2_len
                    s = dot_product(tmp_prod(1:n_bb_e-n_bb_s+1,i),tmp1(n_bb_s:n_bb_e,i,i_cell,i_spin))
                    AS_EXX_stress_local(AS_index) = AS_EXX_stress_local(AS_index) + 2*s
                  enddo
                enddo
              enddo
            endif
          enddo
        endif

        if (calc_deriv .or. AS_stress_on) tmp1 = 0 ! can be re-used from here on

        do i_cell = my_cm_cell_start, n_cells_bvk, my_cm_cell_inc
          if(mult_cm_r1(i_cell)) then
            do i_spin = 1, n_spin

              ! TMP1(1:l2b,ic) = TMPX(1:l2b,ic) * CM(r1,l2,ic)
              ! ----------------------------------------------
              tmp_prod(1:n_bb_e-n_bb_s+1,1:nl2_len) = 0
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          call mult_coul_mat_left(nl2_len, coul_mat_store_SR(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell), &
		                                  tmpx(1,1,i_cell,i_spin),size(tmpx,1), &
		                                  tmp_prod, size(tmp_prod,1))
              else
		          call mult_coul_mat_left(nl2_len, coul_mat_store(vb_atom2atom(i_atom_r1),vb_atom2atom(i_atom_l2),i_cell), &
		                                  tmpx(1,1,i_cell,i_spin),size(tmpx,1), &
		                                  tmp_prod, size(tmp_prod,1))
              end if
              if (calc_deriv .or. AS_stress_on) then
                tmp1(n_bb_s:n_bb_e,1:nl2_len,i_cell,i_spin) = tmp_prod(1:n_bb_e-n_bb_s+1,1:nl2_len)
              else
                tmp2(n_bb_s:n_bb_e,1:nl2_len,i_cell,i_spin) =&
                    tmp2(n_bb_s:n_bb_e,1:nl2_len,i_cell,i_spin) + 2*tmp_prod(1:n_bb_e-n_bb_s+1,1:nl2_len)
              endif
              cell_used(i_cell,i_spin) = .true.
            enddo
          endif

        enddo
! End timing block,   Products #2
ttt_other = mpi_wtime()
times(3) = times(3) + ttt_other - ttt0
! Begin timing block, Other #12

        ! tmp1 contains the contributions to HF12
        ! tmp2 contains the contributions to HF22
        !
        ! For the Fock matrix we need   tmp2 + 2*tmp1 since HF21 is not calculated and HF12 is used twice.
        ! For the derivatives we need 2*tmp2 + 2*tmp1 since the contribution from HF11 is the same
        ! as from HF22 and the contribution from HF11 is also not calculated.

        call aims_deallocate(tmpx, "tmpx")
        ! WPH:  I'm manually allocating and checking tmpx here, as the allocate
        !       infrastructure was returning arrays with incorrect bounds in
        !       the first dimension when compiled with PGI compilers on Titan.
        !       Hopefully this is just a compiler quirk and not something
        !       deeper...
        ! WPH (11 February 2018):  Well here I am again, 1.5 years later.
        !       Updating this module to use the aims memory tracking module,
        !       which was based off the allocate infrastructure causing the bug
        !       back in July 2016, caused the bug to reappear again.  Only this
        !       time, it's occuring on timewarp (our development cluster at
        !       Duke) using PGI 17.4.  So it's *definitely* a compiler quirk
        !       that simply refuses to go away.  Inlining the offending
        !       subroutine fixes the problem.
        ! SK: In the remaining part we can save one array if
        !     calc_deriv .or. AS_stress_on==.false. by changing tmpx <--> tmp2
        !     This saves up to 10% time in case of large tmpx,tmp2,tmp1.
        ! WPH (14 October 2018): Well here I am again, 2.25 years later,
        !       unregressing this PGI bug that keeps regressing.  PGI 17.10 also
        !       experiences this bug, by the way.
        if (calc_deriv .or. AS_stress_on) then
          mem_tmpx = int(nl2_bb_len,kind=8) * int(nl2_len,kind=8) * &
             int(n_cells_bvk,kind=8) * int(n_spin,kind=8) * int(8,kind=8)
          allocate(tmpx(nl2_bb_len, nl2_len, n_cells_bvk, n_spin), stat=info)
          call update_when_allocating(info, "tmpx", mem_tmpx)

          tmpx = 2*tmp2 + 2*tmp1 ! for the derivatives
          tmp2 =   tmp2 + 2*tmp1 ! for the fock matrix
        endif

! End timing block,   Other #12
times_other(12) = times_other(12) + mpi_wtime() - ttt_other
! Begin timing block, Sync/Imbalance #9
        if(my_tasks_per_atom > 1) then
          call time_sync_start(mpi_comm_atom)

          if (.not.(use_mpi_in_place)) then
             call mpi_allreduce(cell_used,cell_used_rcv,n_spin*n_cells_bvk,MPI_LOGICAL,MPI_LOR,mpi_comm_atom,mpierr)
             cell_used = cell_used_rcv
          else
             call mpi_allreduce(mpi_in_place,cell_used,n_spin*n_cells_bvk,MPI_LOGICAL,MPI_LOR,mpi_comm_atom,mpierr)
          endif

          call sync_partial(tmp2, cell_used)
          if(calc_deriv .or. AS_stress_on) call sync_partial(tmpx, cell_used)
          call time_sync_end
        endif
! End timing block,   Sync/Imbalance #9
ttt_other = mpi_wtime()
! Begin timing block, Other #13

        ! Get maximum norm of tmp2, sort with descending norm and store in tmp1

        max_norm_tmp2(:) = 0
        do i_cell = 1, n_cells_bvk
          do i_spin = 1, n_spin
            if(cell_used(i_cell,i_spin)) &
              max_norm_tmp2(i_cell) = max(max_norm_tmp2(i_cell), max_column_norm(tmp2(:,:,i_cell,i_spin)))
          enddo
          aux(1,i_cell) = -max_norm_tmp2(i_cell)
          aux(2,i_cell) = i_cell
        enddo

        call heapsort_general(aux, 2, n_cells_bvk, 1)

        do i_cell = 1, n_cells_bvk
          max_norm_tmp2(i_cell) = -aux(1,i_cell)
          idx_cells(i_cell) = aux(2,i_cell)
          tmp1(:,:,i_cell,:) = tmp2(:,:,idx_cells(i_cell),:)
        enddo

        if(calc_deriv .or. AS_stress_on) then
          ! Get maximum norm of tmpx, presorting doesn't help for the derivatives
          max_norm_tmp(:) = 0
          do i_cell = 1, n_cells_bvk
            do i_spin = 1, n_spin
              if(cell_used(i_cell,i_spin)) &
                max_norm_tmp(i_cell) = max(max_norm_tmp(i_cell), max_column_norm(tmpx(:,:,i_cell,i_spin)))
            enddo
          enddo
        endif

        ! Please note: We could use the fact that only half of the fock cells are needed
        ! but then we would lose the possiblity to use a presorted tmp1

! End timing block,   Other #13
ttt0 = mpi_wtime()
times_other(13) = times_other(13) + ttt0 - ttt_other
! Begin timing block, Products #3
        do i_atom_pair = my_atom_pair_s(i_my_atom), my_atom_pair_e(i_my_atom)

          i_atom_r2 = pair_list(2,i_atom_pair)
          i_cell_2  = pair_list(3,i_atom_pair)

          nr2_off    = atom2basis_off2(i_atom_r2)
          nr2_len    = atom2basis_len2(i_atom_r2)

          ! Get the number of cells needed in tmp1
          do i_cell = 1, n_cells_bvk
          	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
            	if(max_norm_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2)*max_norm_tmp2(i_cell) <= crit_val) exit
            else
            	if(max_norm_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2)*max_norm_tmp2(i_cell) <= crit_val) exit
            end if
          enddo
          n_cells_used = i_cell-1

          ! HF12(r1,1:r2b,ic+ic2) = SUM over l2b ( TMP1(l2b,ic) * OVL1(l2b,1:r2b,ic2) )
          ! HF22(r1,1:r2b,ic+ic2) = SUM over l2b ( TMP2(l2b,ic) * OVL1(l2b,1:r2b,ic2) )
          ! ---------------------------------------------------------------------------
          ! tmp1 contains 2*TMP1+TMP2 from the above pseudocode

          if(n_cells_used>0) then
            do i_spin = 1, n_spin
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          call dgemm('T','N',nr2_len,n_cells_used,nl2_bb_len*nl2_len,1.d0, &
		                  ovlp3fn_SR(i_atom_l2)%m(1,pair_offset(i_atom_l2,i_atom_r2,i_cell_2)+1), nl2_bb_len*nl2_len, &
		                  tmp1(1,1,1,i_spin), nl2_bb_len*nl2_len, &
		                  0.d0,fock_matrix_tmp,size(fock_matrix_tmp,1))
              else
		          call dgemm('T','N',nr2_len,n_cells_used,nl2_bb_len*nl2_len,1.d0, &
		                  ovlp3fn(i_atom_l2)%m(1,pair_offset(i_atom_l2,i_atom_r2,i_cell_2)+1), nl2_bb_len*nl2_len, &
		                  tmp1(1,1,1,i_spin), nl2_bb_len*nl2_len, &
		                  0.d0,fock_matrix_tmp,size(fock_matrix_tmp,1))
              end if

              ! Insert fock_matrix_tmp into fock_matrix_row in the correct order
              do i_cell = 1, n_cells_used
                i_cell_fock = add_cells(idx_cells(i_cell),i_cell_2)
                do ir2 = 1, nr2_len
                  fock_matrix_row(nr2_off+ir2,i_cell_fock,i_spin) = &
                  fock_matrix_row(nr2_off+ir2,i_cell_fock,i_spin) + &
                    fock_matrix_tmp(ir2,i_cell)
                enddo
              enddo
            enddo
          endif

          if(calc_deriv .and. (i_atom_l2 /= i_atom_r2)) then
            do i_cell = 1, n_cells_bvk
              i_cell_fock = add_cells(i_cell,i_cell_2)
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          if(max_norm_tmp(i_cell)*max_norm_d_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2)* &
		             max_abs_dm_row(i_atom_r2,i_cell_fock) <= crit_val_deriv) cycle
              else
		          if(max_norm_tmp(i_cell)*max_norm_d_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2)* &
		             max_abs_dm_row(i_atom_r2,i_cell_fock) <= crit_val_deriv) cycle
              end if
              do i_grad = 1, 3
                do i_spin = 1, n_spin
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		              call dgemv('T',nl2_bb_len*nl2_len,nr2_len,1.d0, &
		                         d_ovlp3fn_SR(i_atom_l2,i_grad)%m(1,pair_offset(i_atom_l2,i_atom_r2,i_cell_2)+1), &
		                         nl2_bb_len*nl2_len, &
		                         tmpx(1,1,i_cell,i_spin), 1, &
		                         0.d0, fock_matrix_tmp, 1)
                  else
		              call dgemv('T',nl2_bb_len*nl2_len,nr2_len,1.d0, &
		                         d_ovlp3fn(i_atom_l2,i_grad)%m(1,pair_offset(i_atom_l2,i_atom_r2,i_cell_2)+1), &
		                         nl2_bb_len*nl2_len, &
		                         tmpx(1,1,i_cell,i_spin), 1, &
		                         0.d0, fock_matrix_tmp, 1)
                  end if

                  ! Add to energy gradient
                  do ir2 = 1, nr2_len
                    s = fock_matrix_tmp(ir2,1)*dm_row(nr2_off+ir2,i_cell_fock,i_spin)
                    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		                 d_exx_ene_SR(i_grad,i_atom_l2) = d_exx_ene_SR(i_grad,i_atom_l2) - 2*s
		                 d_exx_ene_SR(i_grad,i_atom_r2) = d_exx_ene_SR(i_grad,i_atom_r2) + 2*s
                    else
		                 d_exx_ene(i_grad,vb_atom2atom(i_atom_l2)) = d_exx_ene(i_grad,vb_atom2atom(i_atom_l2)) - 2*s
		                 d_exx_ene(i_grad,vb_atom2atom(i_atom_r2)) = d_exx_ene(i_grad,vb_atom2atom(i_atom_r2)) + 2*s
                    end if
                  enddo
                enddo
              enddo
            enddo
          endif

          if(AS_stress_on) then
            do i_cell = 1, n_cells_bvk
              i_cell_fock = add_cells(i_cell,i_cell_2)
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		          if(max_norm_tmp(i_cell)*AS_max_norm_d_ovlp3fn_SR(i_atom_l2,i_atom_r2,i_cell_2)* &
		             max_abs_dm_row(i_atom_r2,i_cell_fock) <= AS_crit_val_deriv) cycle
              else
		          if(max_norm_tmp(i_cell)*AS_max_norm_d_ovlp3fn(i_atom_l2,i_atom_r2,i_cell_2)* &
		             max_abs_dm_row(i_atom_r2,i_cell_fock) <= AS_crit_val_deriv) cycle
              end if
              do AS_index = 1, AS_components, 1
                do i_spin = 1, n_spin
                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		              call dgemv('T',nl2_bb_len*nl2_len,nr2_len,1.d0, &
		                         AS_d_ovlp3fn_SR(i_atom_l2,AS_index)%m(1,pair_offset(i_atom_l2,i_atom_r2,i_cell_2)+1), &
		                         nl2_bb_len*nl2_len, &
		                         tmpx(1,1,i_cell,i_spin), 1, &
		                         0.d0, fock_matrix_tmp, 1)
                  else
		              call dgemv('T',nl2_bb_len*nl2_len,nr2_len,1.d0, &
		                         AS_d_ovlp3fn(i_atom_l2,AS_index)%m(1,pair_offset(i_atom_l2,i_atom_r2,i_cell_2)+1), &
		                         nl2_bb_len*nl2_len, &
		                         tmpx(1,1,i_cell,i_spin), 1, &
		                         0.d0, fock_matrix_tmp, 1)
                  end if

                  do ir2 = 1, nr2_len
                    s = fock_matrix_tmp(ir2,1)*dm_row(nr2_off+ir2,i_cell_fock,i_spin)
                    AS_EXX_stress_local(AS_index) = AS_EXX_stress_local(AS_index) + 2*s
                  enddo
                enddo
              enddo
            enddo
          endif

        enddo
! End timing block,   Products #3
ttt_other = mpi_wtime()
times(3) = times(3) + ttt_other - ttt0
! I have my reasons.
! Begin timing block, Other #14
        call aims_deallocate( tmp1, "tmp1" )
        call aims_deallocate( tmp2, "tmp2" )
        if (calc_deriv .or. AS_stress_on) call aims_deallocate( tmpx, "tmpx" )
! End timing block,   Other #14
times_other(14) = times_other(14) + mpi_wtime() - ttt_other
      enddo ! i_my_atom

      ! Add fock_matrix_row to fock_matrix

! Begin timing block, Sync/Imbalance #10
      call time_sync_start
      call sync_vector(fock_matrix_row,size(fock_matrix_row))
      call time_sync_end
! End timing block,   Sync/Imbalance #10

! Begin timing block, Other #15
ttt_other = mpi_wtime()
      ! Since the fock_matrix matrix is symmetric, we can use a column of the inverse cell
      if(my_atom_id == 0 .and. my_basis_off(i_atom_r1) >= 0) then
        do i_cell = 1, n_cells_bvk
          i_basis = i_basis_r1-atom2basis_off2(i_atom_r1)+my_basis_off(i_atom_r1)
          if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
		      fock_matrix_SR(:,i_basis,i_cell,:) = &
		      fock_matrix_SR(:,i_basis,i_cell,:) + fock_matrix_row(:,inv_cell_bvk(i_cell),:)
          else
		      fock_matrix(:,i_basis,i_cell,:) = &
		      fock_matrix(:,i_basis,i_cell,:) + fock_matrix_row(:,inv_cell_bvk(i_cell),:)
          end if
        enddo
      endif
      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0 .and. .not. lc_wpbeh_lr_run) then
      	lc_wpbeh_lr_run = .true.
      	goto 124
      end if
times_other(15) = times_other(15) + mpi_wtime() - ttt_other
! End timing block,   Other #15
    enddo

! Begin timing block, Other #16
ttt_other = mpi_wtime()
    ! Add final line break to test output in case it is needed.
    if ( mod(n_basis,10).ne.0) then
        write(info_str,'(A)') ' '
        call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    end if

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

! End timing block,   Other #16
ttt0 = mpi_wtime()
times_other(16) = times_other(16) + ttt0 - ttt_other
! Begin timing block, Imbalance #11
      call mpi_barrier(mpi_comm_global, mpierr)
times(5) = times(5) + mpi_wtime()-ttt0
! End timing block,   Imbalance #11
! Begin timing block, Pre/Post #2
ttt0 = mpi_wtime()
!SK
!    do i_basis = 1, n_basis
!    write(info_str,*) fock_matrix(i_basis,:,:,:)
!    write(unit=use_unit,fmt="(A,I8,/,A)") 'fock_matrix: ',i_basis,trim(info_str)
!    enddo
! BL Optimization on ALCF BG/Q MIRA produces wrong results
!    exx_ene = sum(fock_matrix(1:n_basis,1:my_n_basis,1:n_cells_bvk,1:n_spin) &
!                    * dm_cols(1:n_basis,1:my_n_basis,1:n_cells_bvk,1:n_spin))

    exx_ene = 0d0
    exx_ene_SR = 0d0
    do i_spin = 1, n_spin, 1
      do i_cell_bvk = 1, n_cells_bvk, 1
        do my_i_basis = 1, my_n_basis, 1
          do i_basis = 1, n_basis, 1
             exx_ene = exx_ene &
                     + fock_matrix(i_basis, my_i_basis, i_cell_bvk, i_spin) &
                     * dm_cols(i_basis, my_i_basis, i_cell_bvk, i_spin)
	         if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
			      exx_ene_SR = exx_ene_SR &
			             + fock_matrix_SR(i_basis, my_i_basis, i_cell_bvk, i_spin) &
			             * dm_cols(i_basis, my_i_basis, i_cell_bvk, i_spin)
	         end if
          end do
        end do
      end do
    end do

    call sync_real_number(exx_ene)
    call sync_real_number(exx_ene_SR)

    exx_ene = exx_ene*0.25*n_spin
    exx_ene_SR = exx_ene_SR*0.25*n_spin

    write(info_str,'(A,F20.12,A)') 'Full exact exchange energy: ', &
      exx_ene*hartree,' eV.'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )

    if (exx_ene.eq.0.d0) then
       ! Put in a warning. If the EXX energy is zero here, something must have gone wrong.
       write(info_str,'(A)') '*** Error. This value should not normally be zero.'
       call localorb_info ( info_str,use_unit,'(1X,A)', OL_norm  )
       write(info_str,'(A)') '*** Unless you know of a good reason, this may indicate a bug in your MPI library.'
       call localorb_info ( info_str,use_unit,'(1X,A)', OL_norm  )
    end if

    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
 	 	 write(info_str,'(A,F20.12,A)') 'Full exact exchange energy short range: ', &
		   exx_ene_SR*hartree,' eV.'
		 call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    end if

    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
		 if (exx_ene_SR.eq.0.d0) then
		    ! Put in a warning. If the EXX energy is zero here, something must have gone wrong.
		    write(info_str,'(A)') '*** Error. This value should not normally be zero.'
		    call localorb_info ( info_str,use_unit,'(1X,A)', OL_norm  )
		    write(info_str,'(A)') '*** Unless you know of a good reason, this may indicate a bug in your MPI library.'
		    call localorb_info ( info_str,use_unit,'(1X,A)', OL_norm  )
		 end if
    end if

    if(calc_deriv) then
      call sync_vector(d_exx_ene,size(d_exx_ene))
      d_exx_ene = d_exx_ene*0.25*n_spin

      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
      	call sync_vector(d_exx_ene_SR,size(d_exx_ene_SR))
      	d_exx_ene_SR = d_exx_ene_SR*0.25*n_spin
      end if

      !FK: Output of EXX gradients not really necessary, they will show up in the forces summary annyway.
      !    Therefore, changed output priority to low.
      write(info_str,'(A)') 'Full exact exchange energy gradients (eV/A) : '
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_low  )

      do i_atom_r1 = 1, n_atoms, 1
        write(info_str,'(3g25.15)') &
          ( d_exx_ene(i_grad,i_atom_r1)*hartree_over_bohr ,i_grad=1,3,1)
        call localorb_info ( info_str,use_unit,'(2X,A)', OL_low  )
      enddo

      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
      	write(info_str,'(A)') 'Full exact exchange energy gradients short range (eV/A) : '
		   call localorb_info ( info_str,use_unit,'(2X,A)', OL_low  )

		   do i_atom_r1 = 1, n_atoms2, 1
		     write(info_str,'(3g25.15)') &
		       ( d_exx_ene_SR(i_grad,i_atom_r1)*hartree_over_bohr ,i_grad=1,3,1)
		     call localorb_info ( info_str,use_unit,'(2X,A)', OL_low  )
		   enddo
      end if
    endif

    if(AS_stress_on) then
      AS_EXX_stress_local(1:9) = AS_EXX_stress_local(1:9)*0.25*n_spin
      call AS_sync_EXX_stress()
    end if

    ! Set hf_exchange_matr_real/complex

    if(real_eigenvectors)then
      hf_exchange_matr_real(:,:,:,:) = 0.
    else
      hf_exchange_matr_complex(:,:,:,:) = 0.
    endif

    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
		if(real_eigenvectors)then
		  hf_exchange_matr_real_SR(:,:,:,:) = 0.
		else
		  hf_exchange_matr_complex_SR(:,:,:,:) = 0.
		endif
    end if

    call aims_allocate(fock_tmp,n_basis,max_n_basis_sp2,'fock_tmp')
	if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
		call aims_allocate(fock_tmp_SR,n_basis,max_n_basis_sp2,'fock_tmp_SR')
	end if

    if(use_symmetry_reduced_spg)then
      call FT_hf_exchange_matr_sym(my_n_basis,my_basis_off,&
                   fock_tmp,fock_matrix,fock_tmp_SR,fock_matrix_SR)
    else
      do i_atom = 1, n_atoms2

        do i_cell_fock = 1, n_cells_bvk

          do i_spin = 1, n_spin

            if(my_basis_off(i_atom) >= 0) then
              my_s = my_basis_off(i_atom) + 1
              my_e = my_basis_off(i_atom) + atom2basis_len2(i_atom)
              fock_tmp(:,1:atom2basis_len2(i_atom)) = fock_matrix(:,my_s:my_e,i_cell_fock,i_spin)
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                  fock_tmp_SR(:,1:atom2basis_len2(i_atom)) = fock_matrix_SR(:,my_s:my_e,i_cell_fock,i_spin)
              end if
            else
              fock_tmp(:,1:atom2basis_len2(i_atom)) = 0
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                  fock_tmp_SR(:,1:atom2basis_len2(i_atom)) = 0
              end if
            endif

            call sync_vector(fock_tmp, n_basis*atom2basis_len2(i_atom))
            if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                  call sync_vector(fock_tmp_SR, n_basis*atom2basis_len2(i_atom))
            end if

            if(use_scalapack) then

              cfact = k_phase_exx(i_cell_fock,my_k_point)*0.5*n_spin

              do i = 1, atom2basis_len2(i_atom)
                i_basis = atom2basis_off2(i_atom) + i
                if(l_col(i_basis) > 0) then
                  do j = 1, n_basis
                    if(l_row(j)>0) then
                      if(real_eigenvectors)then
                        hf_exchange_matr_real(l_row(j),l_col(i_basis),1,i_spin) = &
                          hf_exchange_matr_real(l_row(j),l_col(i_basis),1,i_spin) + dble(fock_tmp(j,i)*cfact)
                        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          hf_exchange_matr_real_SR(l_row(j),l_col(i_basis),1,i_spin) = &
                                  hf_exchange_matr_real_SR(l_row(j),l_col(i_basis),1,i_spin) + dble(fock_tmp_SR(j,i)*cfact)
                        end if
                      else
                        hf_exchange_matr_complex(l_row(j),l_col(i_basis),1,i_spin) = &
                          hf_exchange_matr_complex(l_row(j),l_col(i_basis),1,i_spin) + fock_tmp(j,i)*cfact
                        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          hf_exchange_matr_complex_SR(l_row(j),l_col(i_basis),1,i_spin) = &
                                  hf_exchange_matr_complex_SR(l_row(j),l_col(i_basis),1,i_spin) + fock_tmp_SR(j,i)*cfact
                        end if
                      endif
                    endif
                  enddo
                endif

                if(l_row(i_basis) > 0) then
                  do j = 1, n_basis
                    if(l_col(j)>0) then
                      if(real_eigenvectors)then
                        hf_exchange_matr_real(l_row(i_basis),l_col(j),1,i_spin) = &
                          hf_exchange_matr_real(l_row(i_basis),l_col(j),1,i_spin) + dble(fock_tmp(j,i)*conjg(cfact))
                        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          hf_exchange_matr_real_SR(l_row(i_basis),l_col(j),1,i_spin) = &
                                  hf_exchange_matr_real_SR(l_row(i_basis),l_col(j),1,i_spin) + dble(fock_tmp_SR(j,i)*conjg(cfact))
                        end if
                      else
                        hf_exchange_matr_complex(l_row(i_basis),l_col(j),1,i_spin) = &
                          hf_exchange_matr_complex(l_row(i_basis),l_col(j),1,i_spin) + fock_tmp(j,i)*conjg(cfact)
                        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          hf_exchange_matr_complex_SR(l_row(i_basis),l_col(j),1,i_spin) = &
                                  hf_exchange_matr_complex_SR(l_row(i_basis),l_col(j),1,i_spin) + fock_tmp_SR(j,i)*conjg(cfact)
                        end if
                      endif
                    endif
                  enddo
                endif
              enddo

            else

              i_k = 0
              do i_k_point = 1, n_k_points
                if (myid == MOD(i_k_point, n_tasks)) then
                  i_k = i_k + 1
                  do i = 1, atom2basis_len2(i_atom)
                    i_basis = atom2basis_off2(i_atom) + i
                    if(real_eigenvectors)then
                      hf_exchange_matr_real(:,i_basis,i_k,i_spin) = &
                        hf_exchange_matr_real(:,i_basis,i_k,i_spin) &
                        + dble(fock_tmp(:,i)*k_phase_exx(i_cell_fock,i_k_point)*0.5*n_spin)
                      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                          hf_exchange_matr_real_SR(:,i_basis,i_k,i_spin) = &
                                    hf_exchange_matr_real_SR(:,i_basis,i_k,i_spin) &
                                    + dble(fock_tmp_SR(:,i)*k_phase_exx(i_cell_fock,i_k_point)*0.5*n_spin)
                      end if
                    else
                      hf_exchange_matr_complex(:,i_basis,i_k,i_spin) = &
                        hf_exchange_matr_complex(:,i_basis,i_k,i_spin) &
                        + fock_tmp(:,i)*k_phase_exx(i_cell_fock,i_k_point)*0.5*n_spin
                      if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  hf_exchange_matr_complex_SR(:,i_basis,i_k,i_spin) = &
                                    hf_exchange_matr_complex_SR(:,i_basis,i_k,i_spin) &
                                    + fock_tmp_SR(:,i)*k_phase_exx(i_cell_fock,i_k_point)*0.5*n_spin
                      end if
                    endif
                  enddo
                endif
              enddo

            endif ! use_scalapack

          enddo ! i_spin
        enddo
      enddo
    endif
    ! Add transposed matrix for making use of symmetries

    if(use_scalapack) then
      ! TODO: Transpose with pdtran, remove ugly transposition above
      if(real_eigenvectors) then
        hf_exchange_matr_real = 0.5*hf_exchange_matr_real
        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
        	hf_exchange_matr_real_SR = 0.5*hf_exchange_matr_real_SR
        end if
      else
        hf_exchange_matr_complex = 0.5*hf_exchange_matr_complex
        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
        	hf_exchange_matr_complex_SR = 0.5*hf_exchange_matr_complex_SR
        end if
      endif
    else
      do i_spin = 1, n_spin
        i_k = 0
        do i_k_point = 1, n_k_points
          if (myid == MOD(i_k_point, n_tasks)) then
            i_k = i_k + 1
            if(real_eigenvectors) then
              hf_exchange_matr_real(1:n_basis,1:n_basis,i_k,i_spin) = &
                0.5*(hf_exchange_matr_real(1:n_basis,1:n_basis,i_k,i_spin) + &
                     transpose(hf_exchange_matr_real(1:n_basis,1:n_basis,i_k,i_spin)))
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
		          hf_exchange_matr_real_SR(1:n_basis,1:n_basis,i_k,i_spin) = &
		            0.5*(hf_exchange_matr_real_SR(1:n_basis,1:n_basis,i_k,i_spin) + &
		                 transpose(hf_exchange_matr_real_SR(1:n_basis,1:n_basis,i_k,i_spin)))
              end if
            else
              hf_exchange_matr_complex(1:n_basis,1:n_basis,i_k,i_spin) = &
                0.5*(hf_exchange_matr_complex(1:n_basis,1:n_basis,i_k,i_spin) + &
                     transpose(conjg(hf_exchange_matr_complex(1:n_basis,1:n_basis,i_k,i_spin))))
              if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
              	hf_exchange_matr_complex_SR(1:n_basis,1:n_basis,i_k,i_spin) = &
		            0.5*(hf_exchange_matr_complex_SR(1:n_basis,1:n_basis,i_k,i_spin) + &
		                 transpose(conjg(hf_exchange_matr_complex_SR(1:n_basis,1:n_basis,i_k,i_spin))))
              end if
            endif
          endif
        enddo
      enddo
      if(n_periodic == 0 .and. n_tasks > 1) then
        ! in this case only task 1 (not 0!) has set the exchange matrix but the hamiltonian
        ! is used on all tasks, therefore broadcast the exchange matrix to all:
        if(myid /= 1) hf_exchange_matr_real = 0
        call sync_vector(hf_exchange_matr_real, size(hf_exchange_matr_real))
        if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
        	if(myid /= 1) hf_exchange_matr_real_SR = 0
        	call sync_vector(hf_exchange_matr_real_SR, size(hf_exchange_matr_real_SR))
        end if
      endif
    endif

    call aims_deallocate( fock_tmp, "fock_tmp" )
    if (allocated(fock_tmp_SR)) call aims_deallocate( fock_tmp_SR, "fock_tmp_SR" )
times(6) = times(6) + mpi_wtime()-ttt0
! End timing block,   Pre/Post #2
! End timing,         Total
times(20) = mpi_wtime()-ttts

! The following is the timing for the body of the subroutine mult_coul_mat_left
times(2) = time_mult_coul_mat
! See comment at beginning of timings
times(3) = times(3) - time_mult_coul_mat


  !-----------------------------------------------------------------------------

    call aims_deallocate( dm_aux,                       "dm_aux" )
    call aims_deallocate( dm_x_o3fn_aux,         "dm_x_o3fn_aux" )

    call aims_deallocate( dm_row,                       "dm_row" )
    call aims_deallocate( dm_cols,                     "dm_cols" )

    call aims_deallocate( ovlp3fn_rcv,             "ovlp3fn_rcv" )
    call aims_deallocate( norm_ovlp3fn_rcv,   "norm_ovlp3fn_rcv" )
    call aims_deallocate( idx_ovlp3fn_rcv,     "idx_ovlp3fn_rcv" )

    call aims_deallocate( max_abs_dm_row,       "max_abs_dm_row" )
    call aims_deallocate( max_abs_dm_cols,     "max_abs_dm_cols" )
    call aims_deallocate( max_fact_tmp1,         "max_fact_tmp1" )
    call aims_deallocate( AS_max_fact_tmp1,   "AS_max_fact_tmp1" )
    call aims_deallocate( max_norm_ovlp3fn_rcv, "max_norm_ovlp3fn_rcv" )

    call aims_deallocate( fock_matrix,         "fock_matrix" )
    if(allocated(fock_matrix_SR)) call aims_deallocate( fock_matrix_SR, "fock_matrix_SR" )
    call aims_deallocate( fock_matrix_row, "fock_matrix_row" )
    call aims_deallocate( tmp_prod,               "tmp_prod" )
    call aims_deallocate( tmp_atom,               "tmp_atom" )
    call aims_deallocate( max_norm_tmp,       "max_norm_tmp" )
    call aims_deallocate( max_norm_tmp2,     "max_norm_tmp2" )
    call aims_deallocate( fock_matrix_tmp, "fock_matrix_tmp" )
    call aims_deallocate( aux,                         "aux" )

    call aims_deallocate( pair_1_list,         "pair_1_list" )

    call aims_deallocate( ovl_fact22, "ovl_fact22" )
    if(allocated(ovl_fact22_SR)) call aims_deallocate( ovl_fact22_SR, "ovl_fact22_SR" )
    call aims_deallocate( ovl_fact12, "ovl_fact12" )
    call aims_deallocate( ovl_fact11, "ovl_fact11" )

    call aims_deallocate( ovl_fact22_d, "ovl_fact22_d" )
    if(allocated(ovl_fact22_d_SR)) call aims_deallocate( ovl_fact22_d_SR, "ovl_fact22_d_SR" )
    call aims_deallocate( ovl_fact12_d, "ovl_fact12_d" )
    call aims_deallocate( ovl_fact11_d, "ovl_fact11_d" )

    call aims_deallocate( AS_ovl_fact22_d, "AS_ovl_fact22_d" )
    if(allocated(AS_ovl_fact22_d_SR)) call aims_deallocate( AS_ovl_fact22_d_SR, "AS_ovl_fact22_d_SR" )
    call aims_deallocate( AS_ovl_fact12_d, "AS_ovl_fact12_d" )
    call aims_deallocate( AS_ovl_fact11_d, "AS_ovl_fact11_d" )

    if (.not.(use_mpi_in_place)) then
       call aims_deallocate( ovl_fact22_rcv, "ovl_fact22_rcv" )
       if(allocated(ovl_fact22_rcv_SR)) call aims_deallocate( ovl_fact22_rcv_SR, "ovl_fact22_rcv_SR" )
       call aims_deallocate( ovl_fact12_rcv, "ovl_fact12_rcv" )
       call aims_deallocate( ovl_fact11_rcv, "ovl_fact11_rcv" )
       call aims_deallocate( ovl_fact22_d_rcv, "ovl_fact22_d_rcv" )
       if(allocated(ovl_fact22_d_SR)) call aims_deallocate( ovl_fact22_d_SR, "ovl_fact22_d_SR" )
       call aims_deallocate( ovl_fact12_d_rcv, "ovl_fact12_d_rcv" )
       call aims_deallocate( ovl_fact11_d_rcv, "ovl_fact11_d_rcv" )
       call aims_deallocate( AS_ovl_fact22_d_rcv, "AS_ovl_fact22_d_rcv" )
       if(allocated(AS_ovl_fact22_d_SR)) call aims_deallocate( AS_ovl_fact22_d_SR, "AS_ovl_fact22_d_SR" )
       call aims_deallocate( AS_ovl_fact12_d_rcv, "AS_ovl_fact12_d_rcv" )
       call aims_deallocate( AS_ovl_fact11_d_rcv, "AS_ovl_fact11_d_rcv" )
    endif

  !-----------------------------------------------------------------------------

    times(19) = times(20) - sum(times(1:18))
    times_other(19) = sum(times_other(1:16))
    times_other(20) = times(19) ! Will be overwritten by the MPI_ALLREDUCE

    write(info_str,'(a,7X,8a12)') &
    '  Times  ',   '   DM_x_o3fn','   CM_x_o3fn','    Products',&
    '        Sync','   Imbalance','    Pre/Post','       Other',&
    '  Total calc'
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,i5,a,20f12.3)')'  | Times',myid,': ',times(1:6),times(19:20)
    call localorb_allinfo(info_str,use_unit,'(A)', OL_norm)

    !BL: debug output
    !do i_atom = 1, n_atoms, 1
    !  write(info_str,'(2(a,i5),a,f12.3)')'  | Check ovl3fn' ,myid, ' Atom = ', i_atom, ': ',SUM(ovlp3fn(i_atom)%n)
    !  call localorb_allinfo(info_str,use_unit,'(A)')
    !end do


    if (use_mpi) then
       if (.not.(use_mpi_in_place)) then
          call mpi_allreduce(times,times_rcv,size(times),MPI_REAL8,MPI_SUM,mpi_comm_global,mpierr)
          times = times_rcv
       else
          call mpi_allreduce(mpi_in_place,times,size(times),MPI_REAL8,MPI_SUM,mpi_comm_global,mpierr)
       endif
    end if


    !SK: output hf_exchange_matr_real
!    do i_spin = 1, n_spin
!     do i_k = 1, n_k_points
!      do i = 1, n_basis
!       if(real_eigenvectors) then
!        write(info_str,*) hf_exchange_matr_real(1:n_basis,i,i_k,i_spin)
!        call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
!       else
!        write(info_str,*) hf_exchange_matr_complex(1:n_basis,i,i_k,i_spin)
!        call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
!       endif
!      enddo
!     enddo
!    enddo

    write(info_str,'(a,20f12.3)') '  | Times  sum: ',times(1:6),times(19:20)
    call localorb_info(info_str,use_unit,'(A)', OL_norm)

    write(info_str,'(a)')
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(2X,a)') "Decomposition of Other based on code block (only on MPI task 0):"
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,7X,8a12)') &
    '  Times  ',   '    Block #1','    Block #2','    Block #3',&
    '    Block #4','    Block #5','    Block #6','    Block #7',&
    '    Block #8'
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,i5,a,8f12.3)')'  | Times',myid,': ',times_other(1:8)
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a)')
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,7X,8a12)') &
    '  Times  ',   '    Block #9','   Block #10','   Block #11',&
    '   Block #12','   Block #13','   Block #14','   Block #15',&
    '   Block #16'
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,i5,a,8f12.3)')'  | Times',myid,': ',times_other(9:16)
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a)')
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,7X,2a12)') &
    '  Times  ',   '   Block Sum','    Expected'
    call localorb_info(info_str,use_unit,'(A)', OL_norm)
    write(info_str,'(a,i5,a,2f12.3)')'  | Times',myid,': ', times_other(19), times_other(20)
    call localorb_info(info_str,use_unit,'(A)', OL_norm)

! DB 10/02/13
    ! as promised above: we have to change back the global meaning of n_atoms
    n_atoms = n_atoms_save

    return

  contains

    subroutine sync_partial(vec, used)

      real*8, intent(inout) :: vec(:,:,:,:)
      logical, intent(in) :: used(:,:)

      integer :: i_spin, i_cell, n

      do i_spin = 1, n_spin
        n = 0
        ! Go forward in cells, thus not overwriting anything
        do i_cell = 1, n_cells_bvk
          if(used(i_cell,i_spin)) then
            n = n+1
            if(n/=i_cell) vec(:,:,n,i_spin) = vec(:,:,i_cell,i_spin)
          endif
        enddo
        if(n>0) call sync_vector(vec(:,:,:,i_spin),size(vec,1)*size(vec,2)*n,mpi_comm_atom)
        ! Go backward in cells, thus not overwriting anything
        do i_cell = n_cells_bvk, 1, -1
          if(used(i_cell,i_spin)) then
            if(n/=i_cell) vec(:,:,i_cell,i_spin) = vec(:,:,n,i_spin)
            n = n-1
          else
            vec(:,:,i_cell,i_spin) = 0
          endif
        enddo
      enddo

    end subroutine sync_partial

    ! Routines for timing syncs including a barrier at the begin
    ! Please note that these modify ttt0 and times()

    subroutine time_sync_start(opt_comm)
      integer, intent(in), optional :: opt_comm
      ttt0 = mpi_wtime()
      if(present(opt_comm)) then
        call mpi_barrier(opt_comm, mpierr)
      else
        call mpi_barrier(mpi_comm_global, mpierr)
      endif
      times(5) = times(5) + mpi_wtime()-ttt0
      ttt0 = mpi_wtime() ! to be used in time_sync_end
    end subroutine time_sync_start

    subroutine time_sync_end
      times(4) = times(4) + mpi_wtime()-ttt0 ! ttt0 must be set in time_sync_start
    end subroutine time_sync_end



  end subroutine evaluate_exchange_matr_realspace_p0

  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/get_bvk_cell_idx
  !  NAME
  !    get_bvk_cell_idx
  !  SYNOPSIS

  integer function get_bvk_cell_idx(idx)

    !  PURPOSE
    !
    !    Returns the BvK cell number for an arbitrary index (which may be outside BvK range)
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(in) :: idx(3)

    !  INPUTS
    !    idx -- index to be mapped
    !  OUTPUTS
    !    function result: BvK cell number
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: idx_local(3)

    idx_local(:) = idx(:)

    ! Map into BvK cell range
    idx_local(:) = idx_local(:) - lbnd_bvk_cell(:)
    idx_local(:) = mod(idx_local(:),n_k_points_xyz_nosym(:))
    where(idx_local(:) < 0) idx_local(:) = idx_local(:) + n_k_points_xyz_nosym(:)
    idx_local(:) = idx_local(:) + lbnd_bvk_cell(:)

    ! Safety check only, may be removed later
    if(any(idx_local(:) < lbnd_bvk_cell(:)) .or. any(idx_local(:) > ubnd_bvk_cell(:))) call aims_stop('get_bvk_cell_idx')

    ! Retrieve index of BvK cell
    get_bvk_cell_idx = bvk_cell_idx(idx_local(1),idx_local(2),idx_local(3))

    end function
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix_p0/max_column_norm
  !  NAME
  !    max_column_norm
  !  SYNOPSIS

  real*8 function max_column_norm(matrix)

    !  PURPOSE
    !
    !    Returns the maximum Euclidean norm of the columns of a matrix
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(in) :: matrix(:,:)

    !  INPUTS
    !    matrix -- matrix for which max column norm is calculated
    !  OUTPUTS
    !    function result: maximum Euclidean norm of the columns
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer i

    max_column_norm = 0

    do i=1, ubound(matrix,2)
      max_column_norm = max(max_column_norm,sum(matrix(:,i)**2))
    enddo

   max_column_norm = sqrt(max_column_norm)

  end function

  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix/evaluate_densmat_hf
  !  NAME
  !    evaluate_densmat_hf
  !  SYNOPSIS

  subroutine evaluate_densmat_hf(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

    !  PURPOSE
    !    Evaluates the density matrix at the k-points and stores it in
    !    hf_exchange_matr_real/complex as an intermediate location
    !  USES

    use hartree_fock_p0, only: hf_exchange_matr_real, hf_exchange_matr_complex
    implicit none

    !  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task):: KS_eigenvector_complex
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers

    !  INPUTS
    !   o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !   o occ_numbers -- occupation of states
    !
    !  OUTPUT
    !   only via hf_exchange_matr_real/complex
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
    !

    ! local variables

    real*8, allocatable :: tmp_r(:,:)
    complex*16, allocatable :: tmp_c(:,:)

    integer :: i_state, i_spin, i_k, i_k_point


    if(real_eigenvectors)then
      call aims_allocate(tmp_r,n_basis,n_states,'tmp_r')
    else
      call aims_allocate(tmp_c,n_basis,n_states,'tmp_c')
    endif

    i_k = 0
    do i_k_point = 1, n_k_points
      if(myid == MOD(i_k_point, n_tasks)) then
        i_k = i_k + 1
        do i_spin = 1, n_spin
          if(real_eigenvectors)then
            do i_state = 1, n_states
              tmp_r(:,i_state) = KS_eigenvector(:,i_state,i_spin,i_k)*occ_numbers(i_state,i_spin,i_k_point)
            enddo
            call dgemm('N','T',n_basis,n_basis,n_states,1.d0,KS_eigenvector(1,1,i_spin,i_k),ubound(KS_eigenvector,1), &
                    tmp_r,ubound(tmp_r,1),0.d0,hf_exchange_matr_real(1,1,i_k,i_spin),ubound(hf_exchange_matr_real,1))
          else
            do i_state = 1, n_states
              tmp_c(:,i_state) = KS_eigenvector_complex(:,i_state,i_spin,i_k)*occ_numbers(i_state,i_spin,i_k_point)
            enddo
            call zgemm('N','C',n_basis,n_basis,n_states,(1.d0,0.d0),KS_eigenvector_complex(1,1,i_spin,i_k), &
                    ubound(KS_eigenvector_complex,1), &
                    tmp_c,ubound(tmp_c,1),(0.d0,0.d0),&
                    hf_exchange_matr_complex(1,1,i_k,i_spin),ubound(hf_exchange_matr_complex,1))
          endif
        enddo
      endif
    enddo

    if(real_eigenvectors)then
      call aims_deallocate( tmp_r, "tmp_r" )
    else
      call aims_deallocate( tmp_c, "tmp_c" )
    endif

  end subroutine evaluate_densmat_hf
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix/read_densmat_hf
  !  NAME
  !    read_densmat_hf
  !  SYNOPSIS

  subroutine read_densmat_hf()

    !  PURPOSE
    !    Reads the density matrix at the k-points from disk and stores it in
    !    hf_exchange_matr_real/complex as an intermediate location
    !  USES
    use hartree_fock_p0, only: hf_exchange_matr_real, hf_exchange_matr_complex
    use restart_elsi, only: elsi_restart_lapack
    use timing_core, only: get_times, get_timestamps
    use timing, only: tot_time_matrix_io, tot_clock_time_matrix_io
    implicit none

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  SOURCE

    integer :: i_spin
    integer :: i_k
    integer :: i_k_point
    integer :: i_row
    integer :: i_col
    real*8 :: t_read
    real*8 :: tmp

    do i_spin = 1,n_spin
       call get_timestamps(t_read,tmp)

       i_k = 0

       do i_k_point = 1,n_k_points
          if(myid == mod(i_k_point,n_tasks)) then
             i_k = i_k+1

             if(real_eigenvectors) then
                call elsi_restart_lapack(i_spin,i_k_point,&
                     hf_exchange_matr_real(:,:,i_k,i_spin))

                hf_exchange_matr_real(:,:,i_k,i_spin) = &
                   hf_exchange_matr_real(:,:,i_k,i_spin)/k_weights(i_k_point)
             else
                call elsi_restart_lapack(i_spin,i_k_point,&
                     hf_exchange_matr_complex(:,:,i_k,i_spin))

                hf_exchange_matr_complex(:,:,i_k,i_spin) = &
                   hf_exchange_matr_complex(:,:,i_k,i_spin)/k_weights(i_k_point)
             endif
          endif
       enddo

       call get_times(t_read,tmp,tot_time_matrix_io,tot_clock_time_matrix_io)

       write(info_str,"(2X,A)") "Finished reading density matrices from file"
       call localorb_info(info_str,use_unit)
       write(info_str,"(2X,A,F10.3,A)") "| Time : ",t_read," s"
       call localorb_info(info_str,use_unit)
       write(info_str,"(A)") ""
       call localorb_info(info_str,use_unit)
    enddo

    ! Populate lower half
    if(real_eigenvectors) then
       do i_col = 1,n_basis-1
          do i_row = i_col+1,n_basis
             hf_exchange_matr_real(i_row,i_col,:,:) = &
                hf_exchange_matr_real(i_col,i_row,:,:)
          enddo
       enddo
    else
       do i_col = 1,n_basis-1
          do i_row = i_col+1,n_basis
             hf_exchange_matr_complex(i_row,i_col,:,:) = &
                hf_exchange_matr_complex(i_col,i_row,:,:)
          enddo
       enddo
    endif

  end subroutine read_densmat_hf
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix/compress_matrix
  !  NAME
  !    compress_matrix
  !  SYNOPSIS

  subroutine compress_matrix(matrix, cm)

    !  PURPOSE
    !    Compresses a Coulomb matrix by removing rows/columns
    !    which are completely below coul_mat_threshold and
    !    stores the result in type(coul_mat_t)
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(in) :: matrix(:,:)
    type (coul_mat_t), intent(inout) :: cm

    !  INPUTS
    !   o matrix -- matrix to be compressed
    !
    !  OUTPUT
    !   o cm -- compressed matrix of type coul_mat_t
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
    !

    integer n_rows, n_cols, i, j
    logical, allocatable :: need_row(:), need_col(:)

    integer :: info
    character(*), parameter :: func = 'compress_matrix'

    allocate(need_row(ubound(matrix,1)), stat=info)
    call check_allocation(info, 'need_row', func)
    allocate(need_col(ubound(matrix,2)), stat=info)
    call check_allocation(info, 'need_col', func)

    ! Total number of uncompressed coulomb matrix entries
    n_coulmat_elems(1) = n_coulmat_elems(1) + ubound(matrix,1)*ubound(matrix,2)

    ! Get rows/columns which have at least 1 element above coul_mat_threshold

    need_row(:) = .false.
    need_col(:) = .false.
    do i=1,ubound(matrix,2)
      do j=1,ubound(matrix,1)
        if(abs(matrix(j,i)) > coul_mat_threshold) then
          need_row(j) = .true.
          need_col(i) = .true.
          n_coulmat_elems(2) = n_coulmat_elems(2) + 1 ! Number of elems above threshold
        endif
      enddo
    enddo

    n_rows = count(need_row(:))
    n_cols = count(need_col(:))

    cm%n_rows = n_rows
    cm%n_cols = n_cols

    call aims_allocate( cm%row_idx, n_rows, 'Coulomb rows' )
    call aims_allocate( cm%col_idx, n_cols, 'Coulomb cols' )

    n_rows = 0
    n_cols = 0
    do i=1,ubound(matrix,1)
      if(need_row(i)) then
        n_rows = n_rows+1
        cm%row_idx(n_rows) = i
      endif
    enddo
    do i=1,ubound(matrix,2)
      if(need_col(i)) then
        n_cols = n_cols+1
        cm%col_idx(n_cols) = i
      endif
    enddo

    call aims_allocate(cm%mat,n_rows,n_cols,'Coulomb submatrix')

    do i=1,n_cols
    do j=1,n_rows
      cm%mat(j,i) = matrix(cm%row_idx(j),cm%col_idx(i))
    enddo
    enddo

    deallocate(need_row)
    deallocate(need_col)

    n_coulmat_elems(3) = n_coulmat_elems(3) + n_rows*n_cols ! Stored elems

  end subroutine compress_matrix
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix/mult_coul_mat_left
  !  NAME
  !    mult_coul_mat_left
  !  SYNOPSIS

  subroutine mult_coul_mat_left(ncols,cm,o3fn,ld_o3fn,cm_x_o3fn,ld_cm_x_o3fn, &
                  mode)

    !  PURPOSE
    !    Multiplies a compressed Coulomb matrix with ovlp3fn.
    !    The multiplication is done at the left, i.e. ovlp3fn^T * coulmat
    !    ATTENTION: The result is added to cm_x_o3fn, i.e. this must be set
    !    to zero before calling this routine
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(in) :: ncols,ld_o3fn,ld_cm_x_o3fn
    integer, intent(in), optional :: mode
    type (coul_mat_t), intent(in) :: cm
    real*8, intent(in)  :: o3fn(ld_o3fn,*)
    real*8, intent(inout) :: cm_x_o3fn(ld_cm_x_o3fn,*)

    !  INPUTS
    !   o ncols -- Number of columns in o3fn
    !   o cm -- compressed coulomb matrix
    !   o o3fn -- ovlp3fn, matrix to multiply with Coulomb matrix
    !   o ld_o3fn -- leading dimension of o3fn
    !   o ld_cm_x_o3fn -- leading dimension of cm_x_o3fn
    !   o mode -- operation mode:
    !                   0 - normal operation
    !                   1 - peak row count storage
    !                   2 - peak col count storage
    !                   3 - array allocation
    !
    !  OUTPUT
    !   o cm_x_o3fn -- result of multiply
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
    !
    integer, save :: rows, cols, col2
    integer :: info
    integer :: i, j
    !real*8, save, allocatable :: tmp_o3fn(:,:), tmp_res(:,:)
    real*8, allocatable :: tmp_o3fn(:,:), tmp_res(:,:)
    real*8 ttt0

    !if (present(mode)) then
    !  if (mode.eq.1) then
    !     rows = cm%n_rows
    !     col2 = ncols
    !  elseif (mode.eq.2.) then
    !     cols = cm%n_cols
    !  elseif (mode.eq.3) then
    !     if (allocated(tmp_o3fn)) deallocate(tmp_o3fn)
    !     if (allocated(tmp_res)) deallocate(tmp_res)
    !     allocate(tmp_o3fn(rows,col2), tmp_res(cols,col2))
    !  endif
    !  return
    !endif

    if(cm%n_rows==0 .or. cm%n_cols==0) return ! Matrix is zero

    ! AJL/Feb2014 (a.logsdail@ucl.ac.uk)
    ! Commented out until I reproduce the error
    ! This next section presents errors on my Intel machine, where I get:
    ! Subscript #1 of the array TMP_O3FN has value 234 which is greater than the upper bound of 233
    ! or such like on some processors as regular occurence, called from line 3352 of this file.
    ! It seems the performance code (mode) above prevents a re-allocation of tmp_o3fn

    ! My only fool-proof work around is to check the matrix is big enough
    ! and deallocate/reallocate if it is not. Alternatives welcomed!

    ! if (size(tmp_o3fn,1).lt.cm%n_rows.or.size(tmp_o3fn,2).lt.ncols) then
       !write(use_unit,*) size(tmp_o3fn,1), cm%n_rows, size(tmp_o3fn,2), ncols
       !write(use_unit,*) 'REDEFINING TMP_O3FN'
    !   deallocate(tmp_o3fn)
       ! Assume the allocation is the same as the comments below...
    !   allocate(tmp_o3fn(cm%n_rows,ncols))
    ! endif
    ! AJL/END

    ! Begin timing block, CM_x_o3fn
    ttt0 = mpi_wtime()

    ! AJL April 2015
    ! I keep getting problems with overflow on these arrays.
    ! Added checks to make sure the process is explicit and correct
    ! AJL
    allocate(tmp_o3fn(cm%n_rows,ncols),stat=info)
    call check_allocation(info, 'tmp_o3fn                      ')
    allocate(tmp_res(cm%n_cols,ncols),stat=info)
    call check_allocation(info, 'tmp_res                       ')

    do i=1,ncols
      do j=1,cm%n_rows
        tmp_o3fn(j,i) = o3fn(cm%row_idx(j),i)
      enddo
    enddo

    ! For performance reasons we calculate the transposed result coulmat^T * ovlp3fn
    call dgemm('T','N', cm%n_cols, ncols, cm%n_rows, 1.d0, cm%mat, ubound(cm%mat,1), &
            tmp_o3fn, ubound(tmp_o3fn,1), 0.d0, tmp_res, ubound(tmp_res,1))

    do i=1,ncols
      do j=1,cm%n_cols
        cm_x_o3fn(cm%col_idx(j),i) = cm_x_o3fn(cm%col_idx(j),i) + tmp_res(j,i)
      enddo
    enddo

    deallocate(tmp_o3fn)
    deallocate(tmp_res)

    time_mult_coul_mat = time_mult_coul_mat + mpi_wtime()-ttt0
    ! End timing block,   CM_x_o3fn

  end subroutine mult_coul_mat_left
  !******
  !----------------------------------------------------------------------------
  !****s* calculate_fock_matrix/mult_coul_mat_right
  !  NAME
  !    mult_coul_mat_right
  !  SYNOPSIS

  subroutine mult_coul_mat_right(ncols,cm,o3fn,ld_o3fn,cm_x_o3fn,ld_cm_x_o3fn)

    !  PURPOSE
    !    Multiplies a compressed Coulomb matrix with ovlp3fn.
    !    The multiplication is done at the right, i.e. coulmat * ovlp3fn
    !    ATTENTION: The result is added to cm_x_o3fn, i.e. this must be set
    !    to zero before calling this routine
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(in) :: ncols,ld_o3fn,ld_cm_x_o3fn
    type (coul_mat_t), intent(in) :: cm
    real*8, intent(in)  :: o3fn(ld_o3fn,*)
    real*8, intent(inout) :: cm_x_o3fn(ld_cm_x_o3fn,*)

    !  INPUTS
    !   o ncols -- Number of columns in o3fn
    !   o cm -- compressed coulomb matrix
    !   o o3fn -- ovlp3fn, matrix to multiply with Coulomb matrix
    !   o ld_o3fn -- leading dimension of o3fn
    !   o ld_cm_x_o3fn -- leading dimension of cm_x_o3fn
    !
    !  OUTPUT
    !   o cm_x_o3fn -- result of multiply
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
    !

    integer :: info
    integer :: i, j
    real*8, allocatable :: tmp_o3fn(:,:), tmp_res(:,:)


    if(cm%n_rows==0 .or. cm%n_cols==0) return ! Matrix is zero

    ! AJL April 2015
    ! I keep getting problems with overflow on these arrays.
    ! Added checks to make sure the process is explicit and correct
    ! AJL
    allocate(tmp_o3fn(cm%n_rows,ncols),stat=info)
    call check_allocation(info, 'tmp_o3fn                      ')
    allocate(tmp_res(cm%n_cols,ncols),stat=info)
    call check_allocation(info, 'tmp_res                       ')
    !allocate(tmp_o3fn(cm%n_cols,ncols), tmp_res(cm%n_rows,ncols))

    do i=1,ncols
      do j=1,cm%n_cols
        tmp_o3fn(j,i) = o3fn(cm%col_idx(j),i)
      enddo
    enddo

    call dgemm('N','N', cm%n_rows, ncols, cm%n_cols, 1.d0, cm%mat, ubound(cm%mat,1), &
            tmp_o3fn, ubound(tmp_o3fn,1), 0.d0, tmp_res, ubound(tmp_res,1))

    do i=1,ncols
      do j=1,cm%n_rows
        cm_x_o3fn(cm%row_idx(j),i) = cm_x_o3fn(cm%row_idx(j),i) + tmp_res(j,i)
      enddo
    enddo


    !deallocate(tmp_o3fn, tmp_res)
    deallocate(tmp_o3fn)
    deallocate(tmp_res)

  end subroutine mult_coul_mat_right

  ! Here follows a very simple (but inefficient) implementation of the
  ! exact exchange energy calculations including gradients.
  ! It can be used to test the rather complicated evaluate_exchange_matr_realspace
  ! since the code below is straight forward.
  ! If used with more than a few atoms/cells, however, it will take forever ...

  subroutine calc_exx_simple(dens, fock, exx_ene, d_exx_ene)

    use basis, only: atom2basis_off, sp2n_basis_sp
    use geometry, only: species
    use prodbas, only: max_n_basbas_sp, atom2basbas_off, sp2n_basbas_sp
    implicit none

    real*8, intent(in)  :: dens(:,:,:,:)
    real*8, intent(out) :: fock(:,:,:,:)
    real*8, intent(out) :: exx_ene
    real*8, intent(out) :: d_exx_ene(3,n_atoms)

    integer :: i_atom_r, i_atom_l, i_atom, i_atom_r1, i_atom_l1, i_atom_r2, i_atom_l2
    integer :: i_cell, i_cell_1, i_cell_2, i_cell_fock, i_spin
    integer :: i_cell11, i_cell12, i_cell21, i_cell22
    integer :: i, k, ii, jj, ii1, jj1, ii2, jj2, nbb1, nbb2, nbbx, irf, icf, ird, icd
    real*8 :: Rvecs(3,1), s, f, dummy(1,1,1,1,1,1)

    real*8 :: o3fn1(max_n_basbas_sp), o3fn2(max_n_basbas_sp)
    real*8 :: d_o3fn1(max_n_basbas_sp,3), d_o3fn2(max_n_basbas_sp,3)

    real*8, allocatable :: cxo1(:,:,:), cxo2(:,:,:), d_cxo1(:,:,:,:), d_cxo2(:,:,:,:)

    logical :: calc_deriv = .true.

    integer :: info
    character(*), parameter :: func = 'calc_exx_simple'


    ! Statement functions for offset of basis/basbas (only for better readability!!!)
    integer :: atom2basis_len, atom2basbas_len
    atom2basis_len(i_atom)  = sp2n_basis_sp(species(i_atom))
    atom2basbas_len(i_atom) = sp2n_basbas_sp(species(i_atom))

    allocate(cxo1(max_n_basbas_sp,n_atoms,n_cells_bvk), stat=info)
    call check_allocation(info, 'cxo1', func, max_n_basbas_sp, n_atoms, n_cells_bvk)
    allocate(cxo2(max_n_basbas_sp,n_atoms,n_cells_bvk), stat=info)
    call check_allocation(info, 'cxo2', func, max_n_basbas_sp, n_atoms, n_cells_bvk)
    if(calc_deriv) then
      allocate(d_cxo1(max_n_basbas_sp,n_atoms,n_cells_bvk,3), stat=info)
      call check_allocation(info, 'd_cxo1', func, max_n_basbas_sp, n_atoms, n_cells_bvk, 3)
      allocate(d_cxo2(max_n_basbas_sp,n_atoms,n_cells_bvk,3), stat=info)
      call check_allocation(info, 'd_cxo2', func, max_n_basbas_sp, n_atoms, n_cells_bvk, 3)
    endif

    fock(1:n_basis,1:n_basis,1:n_cells_bvk,1:n_spin) = 0.
    exx_ene = 0.
    if(calc_deriv) d_exx_ene(:,:) = 0.

    do i_cell_1  = 1, n_cells_bvk
    do i_atom_r1 = 1, n_atoms
    do i_atom_l1 = 1, n_atoms

      write(info_str, *)'Working on pair ',i_atom_l1,i_atom_r1,i_cell_1
      call localorb_info ( info_str,use_unit,'(A)', OL_norm )

      if(.not.pair_flag_bvk(i_atom_l1,i_atom_r1,i_cell_1)) cycle

      do jj1 = 1, atom2basis_len(i_atom_r1)
      do ii1 = 1, atom2basis_len(i_atom_l1)

        nbb1 = atom2basbas_len(i_atom_l1)
        nbb2 = atom2basbas_len(i_atom_r1)

        i_cell = i_cell_1
        k = pair_offset(i_atom_l1,i_atom_r1,i_cell) + (jj1-1)*atom2basis_len(i_atom_l1) + ii1
        o3fn1(1:nbb1) = ovlp3fn(i_atom_l1)%m(1:nbb1,k)

        i_cell = inv_cell_bvk(i_cell_1)
        k = pair_offset(i_atom_r1,i_atom_l1,i_cell) + (ii1-1)*atom2basis_len(i_atom_r1) + jj1
        o3fn2(1:nbb2) = ovlp3fn(i_atom_r1)%m(1:nbb2,k)

        cxo1(:,:,:) = 0
        cxo2(:,:,:) = 0
        d_cxo1(:,:,:,:) = 0
        d_cxo2(:,:,:,:) = 0
        do i_atom = 1, n_atoms
        do i_cell = 1, n_cells_bvk
          nbbx = atom2basbas_len(i_atom)
          call mult_coul_mat_left(1,coul_mat_store(i_atom_l1,i_atom,i_cell),o3fn1,max_n_basbas_sp, &
                                  cxo1(1,i_atom,i_cell),max_n_basbas_sp)
          call mult_coul_mat_left(1,coul_mat_store(i_atom_r1,i_atom,i_cell),o3fn2,max_n_basbas_sp, &
                                  cxo2(1,i_atom,i_cell),max_n_basbas_sp)
          if(calc_deriv) then
            do i=1,3
              call mult_coul_mat_left(1,d_coul_mat_store(i_atom_l1,i_atom,i_cell,i),o3fn1,max_n_basbas_sp, &
                                      d_cxo1(1,i_atom,i_cell,i),max_n_basbas_sp)
              call mult_coul_mat_left(1,d_coul_mat_store(i_atom_r1,i_atom,i_cell,i),o3fn2,max_n_basbas_sp, &
                                      d_cxo2(1,i_atom,i_cell,i),max_n_basbas_sp)
            enddo
          endif
        enddo
        enddo

        do i_cell_fock = 1, n_cells_bvk

          do i_cell_2  = 1, n_cells_bvk
          do i_atom_r2 = 1, n_atoms
          do i_atom_l2 = 1, n_atoms

            if(.not.pair_flag_bvk(i_atom_l2,i_atom_r2,i_cell_2)) cycle

            nbb1 = atom2basbas_len(i_atom_l2)
            nbb2 = atom2basbas_len(i_atom_r2)
            do jj2 = 1, atom2basis_len(i_atom_r2)
            do ii2 = 1, atom2basis_len(i_atom_l2)

              i_cell = i_cell_2
              k = pair_offset(i_atom_l2,i_atom_r2,i_cell) + (jj2-1)*atom2basis_len(i_atom_l2) + ii2
              o3fn1(1:nbb1) = ovlp3fn(i_atom_l2)%m(1:nbb1,k)
              if(calc_deriv) then
                d_o3fn1(1:nbb1,1) = d_ovlp3fn(i_atom_l2,1)%m(1:nbb1,k)
                d_o3fn1(1:nbb1,2) = d_ovlp3fn(i_atom_l2,2)%m(1:nbb1,k)
                d_o3fn1(1:nbb1,3) = d_ovlp3fn(i_atom_l2,3)%m(1:nbb1,k)
              endif

              i_cell = inv_cell_bvk(i_cell_2)
              k = pair_offset(i_atom_r2,i_atom_l2,i_cell) + (ii2-1)*atom2basis_len(i_atom_r2) + jj2
              o3fn2(1:nbb2) = ovlp3fn(i_atom_r2)%m(1:nbb2,k)
              if(calc_deriv) then
                d_o3fn2(1:nbb2,1) = d_ovlp3fn(i_atom_r2,1)%m(1:nbb2,k)
                d_o3fn2(1:nbb2,2) = d_ovlp3fn(i_atom_r2,2)%m(1:nbb2,k)
                d_o3fn2(1:nbb2,3) = d_ovlp3fn(i_atom_r2,3)%m(1:nbb2,k)
              endif

              i_cell11 = i_cell_fock
              i_cell12 = get_bvk_cell_idx(cell_index_bvk(i_cell_fock,:) + cell_index_bvk(i_cell_2,:))
              i_cell21 = get_bvk_cell_idx(cell_index_bvk(i_cell_fock,:) - cell_index_bvk(i_cell_1,:))
              i_cell22 = get_bvk_cell_idx(cell_index_bvk(i_cell_fock,:) + cell_index_bvk(i_cell_2,:) - cell_index_bvk(i_cell_1,:))

              s = dot_product(cxo1(1:nbb1,i_atom_l2,i_cell11),o3fn1(1:nbb1)) &
                + dot_product(cxo1(1:nbb2,i_atom_r2,i_cell12),o3fn2(1:nbb2)) &
                + dot_product(cxo2(1:nbb1,i_atom_l2,i_cell21),o3fn1(1:nbb1)) &
                + dot_product(cxo2(1:nbb2,i_atom_r2,i_cell22),o3fn2(1:nbb2))

              irf = atom2basis_off(i_atom_l1) + ii1
              icf = atom2basis_off(i_atom_l2) + ii2

              ird = atom2basis_off(i_atom_r1) + jj1
              icd = atom2basis_off(i_atom_r2) + jj2

              i_cell = get_bvk_cell_idx(cell_index_bvk(i_cell_fock,:) + cell_index_bvk(i_cell_2,:) - cell_index_bvk(i_cell_1,:))
              do i_spin = 1, n_spin
                fock(irf,icf,i_cell_fock,i_spin) = fock(irf,icf,i_cell_fock,i_spin) + s*dens(ird,icd,i_cell,i_spin)
                f = dens(ird,icd,i_cell,i_spin)*dens(irf,icf,i_cell_fock,i_spin)
                exx_ene = exx_ene + s*f

                if(calc_deriv) then
                  do i=1,3
                    s = dot_product(cxo1(1:nbb1,i_atom_l2,i_cell11),d_o3fn1(1:nbb1,i)) &
                      + dot_product(cxo2(1:nbb1,i_atom_l2,i_cell21),d_o3fn1(1:nbb1,i)) &
                      - dot_product(cxo1(1:nbb2,i_atom_r2,i_cell12),d_o3fn2(1:nbb2,i)) &
                      - dot_product(cxo2(1:nbb2,i_atom_r2,i_cell22),d_o3fn2(1:nbb2,i))
                    d_exx_ene(i,i_atom_r2) = d_exx_ene(i,i_atom_r2) + 2*s*f
                    d_exx_ene(i,i_atom_l2) = d_exx_ene(i,i_atom_l2) - 2*s*f
                    s = dot_product(d_cxo1(1:nbb1,i_atom_l2,i_cell11,i),o3fn1(1:nbb1))
                    d_exx_ene(i,i_atom_l1) = d_exx_ene(i,i_atom_l1) - s*f
                    d_exx_ene(i,i_atom_l2) = d_exx_ene(i,i_atom_l2) + s*f
                    s = dot_product(d_cxo1(1:nbb2,i_atom_r2,i_cell12,i),o3fn2(1:nbb2))
                    d_exx_ene(i,i_atom_l1) = d_exx_ene(i,i_atom_l1) - s*f
                    d_exx_ene(i,i_atom_r2) = d_exx_ene(i,i_atom_r2) + s*f
                    s = dot_product(d_cxo2(1:nbb1,i_atom_l2,i_cell21,i),o3fn1(1:nbb1))
                    d_exx_ene(i,i_atom_r1) = d_exx_ene(i,i_atom_r1) - s*f
                    d_exx_ene(i,i_atom_l2) = d_exx_ene(i,i_atom_l2) + s*f
                    s = dot_product(d_cxo2(1:nbb2,i_atom_r2,i_cell22,i),o3fn2(1:nbb2))
                    d_exx_ene(i,i_atom_r1) = d_exx_ene(i,i_atom_r1) - s*f
                    d_exx_ene(i,i_atom_r2) = d_exx_ene(i,i_atom_r2) + s*f
                  enddo
                endif
              enddo

            enddo ! ii2
            enddo ! jj2
          enddo ! i_atom_l2
          enddo ! i_atom_r2
          enddo ! i_cell_2
        enddo ! i_cell_fock
      enddo ! ii1
      enddo ! jj1
    enddo ! i_atom_l1
    enddo ! i_atom_r1
    enddo ! i_cell_1

    exx_ene = 0.25 * exx_ene * n_spin
    if(calc_deriv) d_exx_ene(:,:) = 0.25*d_exx_ene(:,:) * n_spin

  end subroutine calc_exx_simple

end module calculate_fock_matrix_p0
