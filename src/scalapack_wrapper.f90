! WPH:  The module is misnamed:  it is not a true ScaLAPACK wrapper, but rather it deals
!       almost exclusively with hamiltonian and matrices with similar packings and
!       dimensions (to be precise, leading dimensions.)
!       For using ScaLAPACK/ELPA/ELSI on generic matrices, which have a k-point dependence
!       like hamiltonian, but have different dimensions, please see the subroutines in
!       scalapack_generic_wrapper.f90.

!****h* FHI-aims/scalapack_wrapper
!  NAME
!    scalapack_wrapper -- wrappers for ScaLAPACK library
!  SYNOPSIS
module scalapack_wrapper
!  PURPOSE
!    This module provides the interface between aims and ScaLAPACK
!  USES
  ! WPH: The following use statements should not be here, but these variables
  !      are so widely used throughout scalapack_wrapper that it will take a
  !      significant effort to move them out of global use.  On the other hand,
  !      they're also widely used throughout aims, so leaving them here
  !      shouldn't cause any noticeable problems.
  use dimensions, only: n_basis, n_states, n_spin, n_k_points, &
      n_hamiltonian_matrix_size
  use runtime_choices, only: packed_matrix_format, real_eigenvectors, PM_index
  use localorb_io, only: use_unit
  implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
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

  private

  ! scalapack and blacs static information
  integer, public :: my_scalapack_comm_all
  integer, public :: my_scalapack_comm_work
  integer, public :: my_scalapack_id
  integer, private :: my_scalapack_task_size
  integer, public :: my_blacs_ctxt
  integer, public :: nprow, npcol
  integer, public :: myprow, mypcol
  integer, public :: mxld, mxcol ! dimension of distributed matrices

  integer, public :: my_k_point
  integer, public :: n_scalapack_tasks
  integer, private :: n_tasks_per_host, n_tasks_per_kpoint

  ! parameters for the blacs grid
  integer, parameter, public :: dlen_ = 9
  integer, parameter, public :: csrc = 0
  integer, parameter, public :: rsrc = 0
  integer, public :: mb = 64
  integer, public :: nb = 64

  ! We need only one descriptor, since all matrices have the same shape
  ! and distribution over the processor grid
  integer, dimension(dlen_), public :: sc_desc

  ! storage for real arrays
  real*8, dimension(:,:),   allocatable, public :: ovlp
  real*8, dimension(:,:),   allocatable, public :: ovlp_stored
  real*8, dimension(:,:,:), allocatable, public :: ham
  ! wyj add for save ham, used in elsi_dm_cpscf
  real*8, dimension(:,:,:), allocatable, public :: ham_stored
  real*8, dimension(:,:,:), allocatable, public :: eigenvec
  real*8, dimension(:,:,:), allocatable, private :: eigenvec_stored
  real*8, dimension(:,:,:), allocatable, private :: eigenvec_untrafo

!-------------begin variables for DFPT_phonon_reduced_memory--------------------------
  real*8, dimension(:,:),   allocatable, private :: first_order_ovlp_scalapack
  real*8, dimension(:,:,:),   allocatable, public :: first_order_ham_scalapack
  real*8, dimension(:,:),   allocatable, public :: first_order_U_scalapack
  real*8, dimension(:,:),   allocatable, private :: first_order_edm_scalapack
!-------------end variables for DFPT_phonon_reduced_memory--------------------------

!-------------begin variables for DFPT_polarizability--------------------------
  real*8, dimension(:,:,:,:),   allocatable, private :: first_order_ham_polar_scalapack
  real*8, dimension(:,:,:,:),   allocatable, private :: first_order_U_polar_scalapack
!-------------end variables for DFPT_polarizability--------------------------

!-------------begin variables for DFPT_polar_reduce_memory--------------------------
  !real*8, dimension(:,:,:),   allocatable, private :: first_order_ham_polar_reduce_memory_scalapack
  real*8, dimension(:,:,:),   allocatable, public :: first_order_ham_polar_reduce_memory_scalapack
  !real*8, dimension(:,:,:),   allocatable, private :: first_order_U_polar_reduce_memory_scalapack
  real*8, dimension(:,:,:),   allocatable, public :: first_order_U_polar_reduce_memory_scalapack
!-------------end variables for DFPT_polar_reduce_memory--------------------------

!-------------begin variables for DFPT_dielectric--------------------------
  real*8, dimension(:,:),   allocatable, private :: momentum_matrix_scalapack
  !real*8, dimension(:,:,:),   allocatable, private :: Omega_MO_scalapack
  ! wyj: TODO, just for debug
  real*8, dimension(:,:,:),   allocatable, public:: Omega_MO_scalapack
!-------------end variables for DFPT_dielectric--------------------------

  logical, private :: factor_overlap = .false.
  logical, public :: use_ovlp_trafo = .false.
  logical, public :: full_ovlp_ready = .false.
  integer, public :: n_nonsing_ovlp

  ! Set the following variable to true if ovlp_trafo should
  ! also be used for nonsingular matrices:
  logical, private :: force_use_ovlp_trafo = .false.

  ! Index arrays for distributed matrices

  ! l_row and l_col map global positions to local positions
  ! They have 0 at positions which are not on the local processor

  integer, dimension(:), allocatable, public :: l_row
  integer, dimension(:), allocatable, public :: l_col

  ! my_row and my_col are the counterparts of l_row, l_col and map global positions to local.
  ! n_my_rows, n_my_cols are the exact number of local rows/cols
  ! n_my_rows <= mxld, n_my_cols <= mxcol since mxld, mxcol may be bigger than needed

  integer, public :: n_my_rows, n_my_cols
  integer, allocatable, public :: my_row(:), my_col(:)

  ! storage for complex arrays
  complex*16, dimension(:,:),   allocatable, public :: ovlp_complex
  complex*16, dimension(:,:),   allocatable, public :: ovlp_complex_stored
  complex*16, dimension(:,:,:), allocatable, public :: ham_complex
  ! wyj add
  complex*16, dimension(:,:,:), allocatable, public :: ham_complex_stored
  complex*16, dimension(:,:,:), allocatable, public :: eigenvec_complex
  complex*16, dimension(:,:,:), allocatable, private :: eigenvec_complex_stored

!-------------begin variables for DFPT_phonon_reduced_memory--------------------------
  complex*16, dimension(:,:),   allocatable, private :: first_order_ovlp_complex_scalapack
  complex*16, dimension(:,:,:),   allocatable, public :: first_order_ham_complex_scalapack
  complex*16, dimension(:,:),   allocatable, public :: first_order_U_complex_scalapack
  complex*16, dimension(:,:),   allocatable, private :: first_order_edm_complex_scalapack
!-------------end variables for DFPT_phonon_reduced_memory--------------------------

!-------------begin variables for DFPT_dielectric--------------------------
  complex*16, dimension(:,:),   allocatable, private :: momentum_matrix_complex_scalapack
  complex*16, dimension(:,:,:),   allocatable, private :: Omega_MO_complex_scalapack
!-------------end variables for DFPT_dielectric--------------------------


  ! Global information about process mapping
  integer, allocatable, private :: gl_prow(:), gl_pcol(:), gl_nprow(:), gl_npcol(:)
  integer, private :: max_mxcol

  ! workspace for scalapack eigenvalue solver

  integer, private :: liwork, lrwork, lcwork, len_scalapack_work

  real*8, allocatable, private :: scalapack_work(:)

  ! restart with scalapack
  logical, private :: read_restart(2) = .true.
  integer, private :: last_restart_saving = 0

!-------------begin variables for DFPT_phonon--------------------------
  integer, private :: my_scalapack_id_DFPT_phonon
  integer, private :: my_blacs_ctxt_DFPT_phonon
  integer, private :: nprow_DFPT_phonon, npcol_DFPT_phonon
  integer, private :: myprow_DFPT_phonon, mypcol_DFPT_phonon
  integer, public :: mxld_DFPT_phonon, mxcol_DFPT_phonon ! dimension of distributed matrices

  ! parameters for the blacs grid
  integer, parameter, private :: dlen_DFPT_phonon = 9
  integer, parameter, private :: csrc_DFPT_phonon = 0
  integer, parameter, private :: rsrc_DFPT_phonon = 0
  integer, private :: mb_DFPT_phonon = 64
  integer, private :: nb_DFPT_phonon = 64

  ! shanghui add the second descriptor
  integer, dimension(dlen_), public :: sc_desc_DFPT_phonon

  ! Index arrays for distributed matrices
  ! global positions ===> local positions
  integer, dimension(:), allocatable, private :: l_row_DFPT_phonon
  integer, dimension(:), allocatable, private :: l_col_DFPT_phonon

  ! local positions  ===> global positions
  ! n_my_rows <= mxld, n_my_cols <= mxcol since mxld, mxcol may be bigger than needed
  integer, private :: n_my_rows_DFPT_phonon, n_my_cols_DFPT_phonon
  integer, allocatable, private :: my_row_DFPT_phonon(:), my_col_DFPT_phonon(:)

  integer, private :: liwork_DFPT_phonon, lrwork_DFPT_phonon, &
                   lcwork_DFPT_phonon, len_scalapack_work_DFPT_phonon

  real*8, allocatable, private :: scalapack_work_DFPT_phonon(:)

  logical, private :: factor_overlap_DFPT_phonon = .false.
  logical, private :: ovlp_singular_DFPT_phonon = .false.  !use_ovlp_trafo = .false. shanghui just change name here
  integer, private :: n_nonsing_ovlp_DFPT_phonon


  ! storage for real arrays
  real*8, dimension(:,:),   allocatable, private :: ovlp_supercell_scalapack
  real*8, dimension(:,:,:), allocatable, private :: ham_supercell_scalapack
  real*8, dimension(:,:,:), allocatable, private :: eigenvec_supercell_scalapack
  real*8, dimension(:,:), allocatable, private :: eigenvalues_supercell_scalapack
  real*8, dimension(:,:),   allocatable, private :: first_order_ovlp_supercell_scalapack
  real*8, dimension(:,:), allocatable, private :: first_order_ham_supercell_scalapack
  real*8, dimension(:,:),   allocatable, private :: first_order_U_supercell_scalapack
  real*8, dimension(:,:),   allocatable, private :: first_order_edm_supercell_scalapack
!-------------end variables for DFPT_phonon----------------------------



!-------------------------------------------------------------------------------

  ! Variables for communication when use_local_index is set

  ! The following variables are used for communicating SPARSE local matrices
  ! to/from distributed scalapack arrays

  integer, allocatable, private :: send_idx(:)       ! index into local matrix for elements to send
  integer, allocatable, private :: send_idx_count(:) ! how many elements go to every PE
  integer, allocatable, private :: send_idx_displ(:) ! displacements in send_idx array for every PE
  integer, private :: send_idx_count_tot             ! total number of elements to send

  ! Please note: The send_ccc... arrays are only needed when building the index tables
  !              and immediatly deallocated thereafter
  integer, allocatable, private :: send_ccc(:,:)     ! array coding Cell/Column/Count of idx
  integer, allocatable, private :: send_ccc_count(:) ! how many entries go to every PE
  integer, allocatable, private :: send_ccc_displ(:) ! displacements in send_ccc array for every PE
  integer, private :: send_ccc_count_tot             ! total number of entries to send

  integer, allocatable, private :: recv_row(:)       ! row numbers of received matrix elements,
                                            ! column and cell are coded with recv_ccc
  integer, allocatable, private :: recv_row_count(:) ! how many elements we get from every PE
  integer, allocatable, private :: recv_row_displ(:) ! displacements in recv_row array for every PE
  integer, private :: recv_row_count_tot             ! total number of elements to receive

  integer, allocatable, private :: recv_ccc(:,:)     ! array coding Cell/Column/Count of rows
  integer, allocatable, private :: recv_ccc_count(:) ! how many entries we get from every PE
  integer, allocatable, private :: recv_ccc_displ(:) ! displacements in recv_ccc array for every PE
  integer, private :: recv_ccc_count_tot             ! total number of entries to receive

  ! The following variables are used for communicating FULL local matrices
  ! to/from distributed scalapack arrays.
  ! Please note that in this case no big index arrays are needed
  ! (like send_idx/recv_row above)

  integer, public :: n_basis_local        ! Number of local basis functions
  integer, public :: n_local_matrix_size  ! = n_basis_local*(n_basis_local+1)/2
  integer, allocatable, public :: i_basis_local(:) ! List of local basis functions

  integer, allocatable, public :: basis_row(:) ! List of local basis functions of ALL tasks (relevant rows only)
  integer, allocatable, public :: basis_col(:) ! List of local basis functions of ALL tasks (relevant cols only)
  integer, allocatable, public :: basis_row_limit(:) ! Limits between tasks in basis_row
  integer, allocatable, public :: basis_col_limit(:) ! Limits between tasks in basis_col

  integer, allocatable, public :: send_mat_count(:), send_mat_displ(:)
  integer, allocatable, public :: recv_mat_count(:), recv_mat_displ(:)
  integer, public :: send_mat_count_tot, recv_mat_count_tot

!-------------------------------------------------------------------------------

  ! Description of k-points (grid size and processors involved)
  type k_point_description
    integer :: nprow  ! number of rows for k-point
    integer :: npcol  ! number of cols for k-point
    integer, pointer :: global_id(:,:) ! global ID for every participant of k-point
  end type k_point_description

  type(k_point_description), allocatable, private :: k_point_desc(:)

!-------------------------------------------------------------------------------

  integer, private :: num_ovlp = 0
  integer, private :: num_ham = 0

  ! Variables for ELPA solver
  integer, public :: mpi_comm_rows, mpi_comm_cols
  integer, public :: solver_method_used = 0 ! 0 for ELPA solver not yet determined,
                                            ! 1 for ELPA 1-stage, 2 for ELPA 2-stage

  ! Variables for the Landauer transport:
  complex*16, dimension(:,:), allocatable, private :: green_work, green, gamma, work_gamma
  integer, dimension(:), allocatable, private:: green_ipiv

  ! List of public subroutines
  public :: initialize_scalapack
  public :: reinitialize_scalapack
  public :: scalapack_err_exit
  public :: setup_scalapack_rmatrix
  public :: scalapack_output_global_matrix
  public :: setup_scalapack_full_rmatrix
  public :: setup_scalapack_full_zmatrix
  public :: get_scalapack_global_rmatrix
  public :: get_scalapack_global_zmatrix
  public :: scalapack_pdsyev
  public :: setup_hamiltonian_scalapack
  public :: setup_overlap_scalapack
  public :: construct_hamiltonian_scalapack
  public :: construct_hamiltonian_like_matrix_scalapack
  public :: construct_hamiltonian_like_matrix_zero_diag_scalapack
  public :: construct_overlap_scalapack
  public :: save_overlap_scalapack
  public :: save_ham_scalapack ! add by wyj for CPSCF
  public :: set_sparse_local_ovlp_scalapack
  public :: set_sparse_local_ham_scalapack
  public :: get_sparse_local_matrix_scalapack
  public :: get_set_full_local_matrix_scalapack
  public :: get_set_full_local_matrix_scalapack_cpscf ! add by wyj for CPSCF
  public :: print_ham_cpscf ! add by wyj for CPSCF
  public :: print_sparse_to_dense_global_index_cpscf ! add by wyj for CPSCF
  public :: evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf ! add by wyj for CPSCF
  public :: set_full_local_matrix_scalapack_generic
  public :: init_comm_full_local_matrix_scalapack
  public :: set_sparse_local_matrix_scalapack_generic
  public :: init_comm_sparse_local_matrix_scalapack
  public :: solve_evp_scalapack
  public :: set_full_matrix_real
  public :: solve_evp_scalapack_complex
  public :: set_full_matrix_complex
  public :: extrapolate_dm_scalapack
  public :: normalize_eigenvectors_scalapack
  public :: normalize_eigenvectors_scalapack_real
  public :: normalize_eigenvectors_scalapack_complex
  public :: orthonormalize_eigenvectors_scalapack_real
  public :: check_ev_orthogonality_real
  public :: orthonormalize_eigenvectors_scalapack_complex
  public :: construct_dm_scalapack
  public :: construct_mulliken_decomp_scalapack
  public :: construct_lowdin_decomp_scalapack
  public :: get_sparse_matrix_scalapack
  public :: get_full_matrix_scalapack
  public :: evaluate_scaled_zora_tra_scalapack
  public :: store_eigenvectors_scalapack
  public :: load_eigenvectors_scalapack
  public :: deallocate_scalapack
  public :: finalize_scalapack
  public :: collect_eigenvectors_scalapack
  public :: collect_eigenvectors_scalapack_complex
  public :: spread_eigenvectors_scalapack
  public :: sync_single_eigenvec_scalapack
  public :: sync_single_eigenvec_scalapack_complex
  public :: restart_scalapack_read
  public :: restart_scalapack_write
  public :: construct_ham_and_ovl_transport_scalapack
  public :: construct_greenfunction_scalapack
  public :: solve_greens_functions
  public :: add_self_energy_to_greenfunction_scalapack
  public :: transport_proj_weight
  public :: initialize_scalapack_DFPT_phonon
  public :: finalize_scalapack_DFPT_phonon
  public :: construct_overlap_supercell_scalapack
  public :: construct_hamiltonian_supercell_scalapack
  public :: solve_evp_supercell_scalapack
  public :: construct_first_order_overlap_supercell_scalapack
  public :: construct_first_order_hamiltonian_supercell_scalapack
  public :: evaluate_first_order_U_supercell_scalapack
  public :: construct_first_order_dm_supercell_scalapack
  public :: get_first_order_dm_sparse_matrix_from_supercell_scalapack
  public :: construct_first_order_edm_supercell_scalapack
  public :: get_first_order_edm_sparse_matrix_from_supercell_scalapack
  public :: construct_first_order_overlap_scalapack
  public :: get_first_order_overlap_sparse_matrix_scalapack
  public :: construct_first_order_hamiltonian_scalapack
  public :: get_first_order_hamiltonian_sparse_matrix_scalapack
  public :: construct_first_order_dm_scalapack
  public :: evaluate_first_order_U_scalapack
  public :: get_first_order_dm_complex_sparse_matrix_scalapack
  public :: construct_first_order_edm_scalapack
  public :: get_first_order_edm_complex_sparse_matrix_scalapack
  public :: construct_first_order_dm_polar_scalapack       ! polar  
  public :: get_first_order_dm_polar_scalapack             ! polar
  public :: construct_first_order_ham_polar_scalapack      ! polar
  public :: evaluate_first_order_U_polar_scalapack         ! polar 
  public :: construct_first_order_dm_polar_reduce_memory_scalapack    ! polar_reduce_memory
  public :: get_first_order_dm_polar_reduce_memory_scalapack          ! polar_reduce_memory
  public :: get_first_order_dm_polar_reduce_memory_for_elsi_scalapack ! polar_reduce_memory
  public :: construct_first_order_ham_polar_reduce_memory_scalapack   ! polar_reduce_memory
  public :: evaluate_first_order_U_polar_reduce_memory_scalapack      ! polar_reduce_memory
  public :: construct_momentum_matrix_dielectric_scalapack
  public :: construct_momentum_matrix_dielectric_for_elsi_scalapack
  public :: construct_first_order_hamiltonian_dielectric_scalapack
  public :: construct_first_order_hamiltonian_dielectric_for_elsi_scalapack
  public :: construct_first_order_dm_dielectric_scalapack
  public :: evaluate_first_order_U_dielectric_scalapack
  public :: get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
  public :: get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack
  public :: set_full_matrix_complex_L_to_U
  public :: set_full_matrix_real_L_to_U
  public :: construct_hamiltonian_real_for_elsi_scalapack
  public :: construct_overlap_real_for_elsi_scalapack
  public :: construct_first_order_hamiltonian_polar_for_elsi_scalapack
  public :: get_first_order_dm_polar_for_elsi_scalapack
!-------------------------------------------------------------------------------
  ! List of PRIVATE (= module internal) routines
  private :: scalapack_output_local_matrix
  private :: get_set_sparse_local_matrix_scalapack
  private :: set_local_index_comm_arrays
  private :: output_matrix_real
  private :: construct_overlap_like_matrix_scalapack_die
  private :: diagonalize_overlap_scalapack_real
  private :: diagonalize_overlap_scalapack_complex
  private :: orthonormalize_eigenvectors_scalapack_real_GS
  private :: orthonormalize_eigenvectors_scalapack_complex_GS
  private :: check_ev_orthogonality_complex
  private :: construct_dm_Pulay_forces_scalapack
  private :: construct_forces_dm_scalapack
  private :: remove_stored_eigenvectors_scalapack
  private :: collect_generic_eigenvectors_scalapack
  private :: deallocate_scalapack_DFPT_phonon
!-----------------------------------------------------------------------------------

!******
contains
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/initialize_scalapack
!  NAME
!    initialize_scalapack
!  SYNOPSIS
  subroutine initialize_scalapack()
!  PURPOSE
!    Initialize the ScaLAPACK environment.
!  USES
    use aims_memory_tracking, only: aims_allocate
    use mpi_utilities, only: get_my_processor
    use pbc_lists
    use synchronize_mpi_basic, only: sync_integer_vector, sync_find_max, &
        sync_int_vector
    use localorb_io
    use mpi_tasks
    use dimensions
    use runtime_choices
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    ScaLAPACK communicators, block sizes and local storage arrays plur BLACS grids
!    and the local indexing arrays are set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


    integer :: i_task
    integer :: max_npcol, max_nprow
    integer :: i, lc, lr, n
    integer :: np0, nq0, trilwmin, lwormtr
    integer :: block_size, info
    integer :: mpierr

    character(LEN=MPI_MAX_PROCESSOR_NAME) :: my_proc_name
    character(LEN=MPI_MAX_PROCESSOR_NAME), allocatable :: all_proc_names(:)
    integer :: my_proc_name_len
    integer :: my_master_proc
    integer :: my_host_rank, n_hosts
    integer :: i_k_point, k_dim(2)
    integer, allocatable :: master_proc_tasks(:)
    logical :: default_kpoint_parallelism

    integer, external :: numroc

    character*200 :: info_str


    ! Check which solver should be used, this should go into a config file switch

    if(.not. use_elsi) then
       if(use_elpa) then
          if(myid==0) write(use_unit, *) ' Using modified Scalapack Eigenvalue solver (ELPA)'
       else
          if(myid==0) write(use_unit, *) ' Using standard Scalapack Eigenvalue solver'
       endif
    endif

    default_kpoint_parallelism = .true.

    ! Check if there are exactly as many hosts (SMP-nodes) as k-points.
    ! In this case the tasks will split in a way that each SMP node works on 1 k-point.
    ! This will be done only if
    ! - there is more than 1 k-point
    ! - the assumption that there are as many SMP-nodes as k-points is reasonable
    !   We assume an SMP node has at most 512 MPI tasks in order to avoid
    !   large useless searches in the case of a huge number of tasks.

    ! JM:
    ! introduced max_tasks_per_smp_node with previously hard coded value of 512 as default
    ! this keyword might not yet have an appropriate name for what it does (or could do)
    ! and might be dropped in the future when this if statement might disappear altogether?
    ! there do not seem be any O(n_tasks^2) operations here,
    ! so what is the performance bottleneck when n_tasks is large(r)?!


    ! Get my processor name
    ! Matti:
    ! I moved this out from the if(n_k_points ... ) as the my_proc_name is always needed.
    ! If not done trim(my_proc_name) might crash later in code.

    call get_my_processor(my_proc_name, my_proc_name_len)
    my_proc_name(my_proc_name_len+1:) = ' ' ! for safety: pad with blanks

    if(n_k_points > 1 .and. n_tasks <= max_tasks_per_smp_node*n_k_points) then

       write(info_str, '(2X,A)') 'Checking for possible mapping of k-points to SMP-nodes'
       call localorb_info(info_str, use_unit, '(A)')

       ! Every task get a list of all processor names

       allocate(all_proc_names(0:n_tasks-1))

       call mpi_allgather(my_proc_name,   MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                          all_proc_names, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                          mpi_comm_global, mpierr)

       ! Search my master proc, i.e. the first proc with the same proc name as mine

       do i_task = 0, n_tasks-1
          if(all_proc_names(i_task) == my_proc_name) then
             my_master_proc = i_task
             exit
          endif
       enddo

       deallocate(all_proc_names) ! potentially huge and not needed any more ...

       ! generate a list where all master procs are flagged

       allocate(master_proc_tasks(0:n_tasks-1))
       master_proc_tasks(:) = 0
       master_proc_tasks(my_master_proc) = 1

       ! When syncing this list, we get also the number of tasks per master proc

       call sync_integer_vector(master_proc_tasks, n_tasks)

       ! get the total number of master procs as well as the rank of my master proc

       n_hosts = 0
       my_host_rank = 0
       do i_task = 0, n_tasks-1
          if(master_proc_tasks(i_task) > 0) then
             n_hosts = n_hosts+1
             if(i_task==my_master_proc) my_host_rank = n_hosts ! my_host_rank starts at 1, not 0
          endif
       enddo

       ! Safety check only:
       if(my_host_rank==0) call aims_stop('Internal error getting master procs')

       write(info_str, '(2X,A,I0,A)') 'Running on ',n_hosts,' SMP-nodes'
       call localorb_info(info_str, use_unit, '(A)')

       if(n_hosts == n_k_points) then

          ! There are exactly as many smp nodes as k-points

          write(info_str, '(2X,A)') &
             'Mapping of k-points to SMP-nodes is possible:'
          call localorb_info(info_str, use_unit, '(A)')
          write(info_str, '(2X,A,I0,A,I0,A)') &
             '| ', n_k_points, ' k-points  :  ', n_hosts, ' SMP-nodes)'

          my_k_point = my_host_rank

          n_scalapack_tasks = n_tasks

          default_kpoint_parallelism = .false.

          i_k_point = 0
          do i_task = 0, n_tasks-1
             if(master_proc_tasks(i_task) > 0) then
                i_k_point=i_k_point+1
                k_point_loc(1,i_k_point)=i_task
                k_point_loc(2,i_k_point)=1
             endif
          enddo

       else

          ! JM: There are less or more smp nodes than k-points

          write(info_str, '(2X,A)') 'Mapping of k-points to SMP-nodes is *NOT* possible:'
          call localorb_info(info_str, use_unit, '(A)')
          write(info_str, '(2X,A,I0,A,I0,A)') &
             '| ', n_k_points, ' k-points  : ', n_hosts, ' SMP-nodes'
          call localorb_info(info_str, use_unit, '(A)')

          if(restrict_kpoint_to_smp_node) then

             write(info_str, '(2X,A)') 'Restricting each k-point to a single SMP node.'
             call localorb_info(info_str, use_unit, '(A)')

             ! JM: implicitly assuming equal number of tasks per host here
             n_tasks_per_host = n_tasks / n_hosts

             if(n_hosts < n_k_points) then

                ! JM: (n_hosts < n_k_points) yields (n_tasks_per_host > n_tasks_per_kpoint)
                do n_tasks_per_kpoint = n_tasks / n_k_points, 1, -1
                   ! JM: this always succeeds for n_tasks_per_kpoint == 1
                   if (MOD(n_tasks_per_host,n_tasks_per_kpoint)==0) exit
                enddo

             else

                ! JM: (n_hosts > n_k_points) should always be true here

                n_tasks_per_kpoint = n_tasks_per_host

             endif

             n_scalapack_tasks = n_k_points * n_tasks_per_kpoint

             if(myid < n_scalapack_tasks) then
                my_k_point = myid*n_k_points/n_scalapack_tasks + 1
             else
                my_k_point = MPI_UNDEFINED
             endif

             default_kpoint_parallelism = .false.

             do i_k_point = 1, n_k_points
                k_point_loc(1,i_k_point)=(i_k_point-1)*n_tasks_per_kpoint
                k_point_loc(2,i_k_point)=1
             enddo

          else

             write(info_str, '(1X,A,A)') '* The calculation will run, ', &
                                           'but some performance degradation might result.'
             call localorb_info(info_str, use_unit, '(A)')
             write(info_str, '(1X,A,A)') '* Check the number of actual k-points above, ', &
                                           'and consider using a number'
             call localorb_info(info_str, use_unit, '(A)')
             write(info_str, '(1X,A,A)') '* compute nodes that can be divided cleanly ', &
                                           'by the number of k points'
             call localorb_info(info_str, use_unit, '(A)')
             write(info_str, '(1X,A)') '* (e.g., two k-points per node, two nodes per k-point, etc.).'
             call localorb_info(info_str, use_unit, '(A)')

          endif

       endif

       deallocate(master_proc_tasks)

    endif

    if(default_kpoint_parallelism) then
       n_scalapack_tasks = n_tasks

       ! for the non-periodic case, this message makes no sense.
       if (n_periodic.gt.0) then
         write(info_str, '(2X,A)') &
             'Using simple linear distribution as default method for k-point parallelism.'
         call localorb_info(info_str, use_unit, '(A)')
       end if
       my_k_point = myid*n_k_points/n_tasks + 1

       do i_k_point=1,n_k_points
          k_point_loc(1,i_k_point)=ceiling((i_k_point-1.)*n_tasks/n_k_points)
          k_point_loc(2,i_k_point)=1
       end do

       ! Quick check that (default) k-point setting is correct:
       if(my_k_point<1 .or. my_k_point>n_k_points) then
          call aims_stop('Internal error setting k-points using default method.')
       endif
    endif

    write(info_str, '(2X,A,I0,A)') &
          '* Using ', n_scalapack_tasks, ' tasks for Scalapack Eigenvalue solver.'
    call localorb_info(info_str, use_unit, '(A)')

    write(info_str, '(2X,a)') 'Detailed listing of tasks and assigned k-points:'
    call localorb_info(info_str, use_unit, '(A)')
    if (n_periodic.eq.0) then
      write(info_str, '(2X,a)') '(for non-periodic systems, the "k-point" denotes an internal label only)'
      call localorb_info(info_str, use_unit, '(A)')
    end if
    if (my_k_point == MPI_UNDEFINED) then
       write(info_str, '(3X,a,i5,a,a)') &
             'Task ',myid,' k-point NONE on ',trim(my_proc_name)
    else
       write(info_str, '(3X,a,i5,a,i5,a,a)') &
             'Task ',myid,' k-point ',my_k_point,' on ',trim(my_proc_name)
    endif
    call localorb_allinfo(info_str, use_unit, '(A)')

    if(myid < n_scalapack_tasks) then
       call MPI_Comm_split(mpi_comm_global, my_k_point, myid, my_scalapack_comm_all, mpierr )
       ! JM: leads to MPI errors on nodes which do not participate in Scalapack k-point parallelism
       call MPI_Comm_size( my_scalapack_comm_all, my_scalapack_task_size, mpierr)
       call MPI_Comm_rank( my_scalapack_comm_all, my_scalapack_id, mpierr)
    else
       ! JM: Are these reasonable defaults for tasks idling in Scalapack operations?
       my_scalapack_comm_work = MPI_ERR_COMM
       my_scalapack_task_size = 0
       my_scalapack_id = MPI_UNDEFINED
    endif

    ! JM: status of restrict_kpoint_to_smp_node above
    ! - code seems to succeed until here
    ! - apparently other parts cannot cope with not participating in Scalapack operations

    ! divide the BLACS grid into rows and columns for each task
    do npcol = NINT(sqrt(dble(my_scalapack_task_size))), 2, -1
       if (MOD(my_scalapack_task_size,npcol)==0) exit
    enddo
    ! at the end of the above loop, my_scalapack_task_size is always divisible by npcol
    nprow = my_scalapack_task_size/npcol ! always succeeds without remainder

    if (my_scalapack_id==0)  then
       if (n_periodic.gt.0) then
         write(info_str, '(5(a,i6))') '  K-point:',my_k_point,' Tasks:',my_scalapack_task_size, &
               ' split into ',nprow,' X ',npcol,' BLACS grid'
       else
         write(info_str, '(4(a,i6))') '  Tasks:',my_scalapack_task_size, &
               ' split into ',nprow,' X ',npcol,' BLACS grid'
       end if
       call localorb_info(info_str, use_unit, '(A)')
    end if

    ! If the number of working processors is smaller the total number
    ! of processors for the k-point we have to split the communicator again

    ! RJ: This is currently not necessary since always: npcol*nprow == my_scalapack_task_size
    ! We leave the code below here since it might make sense to use a smaller number
    ! of tasks than possible if this leads to a better divison (e.g. 4 x 4 = 16 instead of 17)

    if(npcol*nprow < my_scalapack_task_size) then

       if(my_scalapack_id<npcol*nprow) then
          n = 1
       else
          n = MPI_UNDEFINED
       endif
       call MPI_Comm_split( my_scalapack_comm_all, n, &
         my_scalapack_id, my_scalapack_comm_work, mpierr )

    else

       my_scalapack_comm_work = my_scalapack_comm_all

    endif

    ! initialize the BLACS grid
    if(my_scalapack_id<npcol*nprow) then
       my_blacs_ctxt = my_scalapack_comm_work
       call BLACS_Gridinit( my_blacs_ctxt, 'R', nprow, npcol )
       call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, myprow, mypcol )
    else
       myprow = -1
       mypcol = -1
    endif

    ! Allocate and set global process mapping info

    allocate(gl_prow(0:n_tasks-1),stat=info)
    call check_allocation(info, 'gl_prow                       ')

    allocate(gl_pcol(0:n_tasks-1),stat=info)
    call check_allocation(info, 'gl_pcol                       ')

    allocate(gl_nprow(0:n_tasks-1),stat=info)
    call check_allocation(info, 'gl_nprow                      ')

    allocate(gl_npcol(0:n_tasks-1),stat=info)
    call check_allocation(info, 'gl_npcol                      ')


    gl_prow = 0
    gl_pcol = 0
    gl_nprow = 0
    gl_npcol = 0

    gl_prow(myid) = myprow
    gl_pcol(myid) = mypcol
    gl_nprow(myid) = nprow
    gl_npcol(myid) = npcol

    call sync_integer_vector(gl_prow,n_tasks)
    call sync_integer_vector(gl_pcol,n_tasks)
    call sync_integer_vector(gl_nprow,n_tasks)
    call sync_integer_vector(gl_npcol,n_tasks)

    ! Calculate blocksize based on n_basis and nprow/npcol
    ! We want that every processor owns a part of the matrix,
    ! i.e. mxld>0 and mxcol>0 everywhere.
    ! This is needed to catch a "bug" in ScaLAPACK:
    ! If a processor doesn't own a part of a matrix, the results are not
    ! distributed to this one
    ! Theoretically, nprow/npcol can be different for different k-points,
    ! so we have to get the global maximum.

    call sync_find_max(nprow,max_nprow)
    call sync_find_max(npcol,max_npcol)

    write(info_str, *) ' Calculating block size based on n_basis = ',n_basis, &
                       ' max_nprow = ',max_nprow,' max_npcol = ',max_npcol
    call localorb_info(info_str, use_unit, '(A)')

    block_size = 1 ! Minimum permitted size
    if(block_size*MAX(max_nprow,max_npcol) > n_basis) then
       write(info_str, *) 'ERROR: n_basis = ',n_basis,' too small for this processor grid'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    endif

    ! Increase blocksize to maximum possible size or 64

    do while (2*block_size*MAX(max_nprow,max_npcol) <= n_basis &
          .and. block_size<64)
       block_size = 2*block_size
    end do

    ! use_elpa works best with SMALL blocksizes!
    if(use_elpa .and. block_size>16) block_size = 16

    ! If block_size is defined in control.in file, check it here.
    if(scalapack_block_size > 0) then ! Set by user.
       if(scalapack_block_size*(MAX(max_nprow,max_npcol)-1) .ge. n_basis) then
          ! See comments above!
          write(info_str, *) &
             'ERROR: User defined block size ',scalapack_block_size,&
             ' too large for this process grid and n_basis.'
          call localorb_info(info_str, use_unit, '(A)')
          call aims_stop
       else
          ! (Hopefully) reasonable value is accepted here.
          block_size = scalapack_block_size
          write(info_str, *) &
             ' Use block size defined in control.in file.'
          call localorb_info(info_str, use_unit, '(A)')
       endif
    endif

    nb = block_size
    mb = block_size

    write(info_str, *) ' ScaLAPACK block size set to: ',block_size
    call localorb_info(info_str, use_unit, '(A)')

    ! initialize the Scalapack descriptor

    if(my_scalapack_id<npcol*nprow) then

       mxld = numroc( n_basis, mb, myprow, rsrc, nprow )
       mxcol = numroc( n_basis, nb, mypcol, csrc, npcol )

! RJ: If mxld/mxcol are too small, they *might* trigger an error in the
! Intel/Scalapack implementation, so set them to a at least 64:

! BL: This is potentially dangerous in case the local
!     dimensions are smaller than n_basis: Horror for parallel writing based on
!     patterns

!       if(mxld  < 64) mxld  = 64
!       if(mxcol < 64) mxcol = 64

       call descinit( sc_desc, n_basis, n_basis, mb, nb, rsrc, csrc, &
            my_blacs_ctxt, MAX(1,mxld), info )

       ! Safety check only, the following should never happen
       if(mxld<=0 .or. mxcol<=0) then
          write(use_unit,*) 'ERROR Task #',myid,' mxld= ',mxld,' mxcol= ',mxcol
          call mpi_abort(mpi_comm_global,1,mpierr)
       endif

    else
       mxld = 1
       mxcol = 1
    endif

    ! allocate and set index arrays

    allocate(l_row(n_basis),stat=info)
    call check_allocation(info, 'l_row                         ')

    allocate(l_col(n_basis),stat=info)
    call check_allocation(info, 'l_col                         ')

    ! Mapping of global rows/cols to local

    l_row(:) = 0
    l_col(:) = 0

    ! ATTENTION: The following code assumes rsrc==0 and csrc==0 !!!!
    ! For processors outside the working set, l_row/l_col will stay
    ! completely at 0

    lr = 0 ! local row counter
    lc = 0 ! local column counter

    do i = 1, n_basis

      if( MOD((i-1)/mb,nprow) == myprow) then
        ! row i is on local processor
        lr = lr+1
        l_row(i) = lr
      endif

      if( MOD((i-1)/nb,npcol) == mypcol) then
        ! column i is on local processor
        lc = lc+1
        l_col(i) = lc
      endif

    enddo

    ! Mapping of local rows/cols to global

    n_my_rows = lr
    n_my_cols = lc
    allocate(my_row(n_my_rows))
    allocate(my_col(n_my_cols))
    lr = 0
    lc = 0

    do i = 1, n_basis
       if(l_row(i)>0) then; lr = lr+1; my_row(lr) = i; endif
       if(l_col(i)>0) then; lc = lc+1; my_col(lc) = i; endif
    enddo

    ! Set up description of all k_points

    allocate(k_point_desc(n_k_points))

    do i_k_point = 1, n_k_points

       ! Get number of rows/cols for current k point

       k_dim(:) = 0
       if(i_k_point==my_k_point .and. myprow==0 .and. mypcol==0) then
          k_dim(1) = nprow
          k_dim(2) = npcol
       endif
       call sync_int_vector(k_dim, 2)

       k_point_desc(i_k_point)%nprow = k_dim(1)
       k_point_desc(i_k_point)%npcol = k_dim(2)

       ! Get global id for all rows/cols of this k-point

       allocate(k_point_desc(i_k_point)%global_id(0:k_dim(1)-1,0:k_dim(2)-1))
       k_point_desc(i_k_point)%global_id(:,:) = 0
       if(i_k_point==my_k_point) k_point_desc(i_k_point)%global_id(myprow, mypcol) = myid
       call sync_int_vector(k_point_desc(i_k_point)%global_id, k_dim(1)*k_dim(2))

    enddo


    ! Allocate scalapack arrays
    if(real_eigenvectors)then

       call aims_allocate(ovlp, mxld, mxcol, "+ovlp")
       call aims_allocate(ham, mxld, mxcol, n_spin, "+ham")
       call aims_allocate(eigenvec, mxld, mxcol, n_spin, "+eigenvec")
       call aims_allocate(eigenvec_complex, 1, 1, 1, "eigenvec_complex")

       if (use_cg) allocate(eigenvec_untrafo(mxld,mxcol,n_spin))


       if(use_DFPT_phonon_reduce_memory) then
         allocate(first_order_ovlp_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_ovlp_scalapack                  ')
         allocate(first_order_ham_scalapack(mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_ham_scalapack                  ')
         allocate(first_order_U_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_U_scalapack                  ')
         allocate(first_order_edm_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_edm_scalapack                  ')

         first_order_ovlp_scalapack = 0.0d0
         first_order_ham_scalapack = 0.0d0
         first_order_U_scalapack = 0.0d0
         first_order_edm_scalapack = 0.0d0
       endif

       if(use_DFPT_polarizability) then
         allocate(first_order_ham_polar_scalapack(3,mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_ham_polar_scalapack                  ')
         allocate(first_order_U_polar_scalapack(3,mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_U_polar_scalapack                  ')

         first_order_ham_polar_scalapack = 0.0d0
         first_order_U_polar_scalapack = 0.0d0
       endif

       if(use_DFPT_polar_reduce_memory) then
         allocate(first_order_ham_polar_reduce_memory_scalapack(mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_ham_polar_reduce_memory_scalapack                  ')
         allocate(first_order_U_polar_reduce_memory_scalapack(mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_U_polar_reduce_memory_scalapack                  ')

         first_order_ham_polar_reduce_memory_scalapack = 0.0d0
         first_order_U_polar_reduce_memory_scalapack = 0.0d0
       endif

       if(use_DFPT_dielectric) then
         allocate(momentum_matrix_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'momentum_matrix_scalapack                  ')
         allocate(Omega_MO_scalapack(mxld,mxcol,3),stat=info)
         call check_allocation(info, 'Omega_MO_scalapack                  ')
         allocate(first_order_ham_scalapack(mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_ham_scalapack                  ')
         allocate(first_order_U_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_U_scalapack                  ')

         momentum_matrix_scalapack = 0.0d0
         Omega_MO_scalapack = 0.0d0
         first_order_ham_scalapack = 0.0d0
         first_order_U_scalapack = 0.0d0
       endif


       ! Safety only:
       ovlp     = 0
       ham      = 0
       eigenvec = 0
       if (use_cg) eigenvec_untrafo = 0.0d0

    else

       call aims_allocate(ovlp_complex, mxld ,mxcol, "+ovlp_complex")
       call aims_allocate(ham_complex, mxld, mxcol, n_spin, "+ham_complex")

       call aims_allocate(eigenvec, 1, 1, 1, "eigenvec")
       call aims_allocate(eigenvec_complex, mxld, mxcol, n_spin, "+eigenvec_complex")

       if(use_DFPT_phonon_reduce_memory) then
         allocate(first_order_ovlp_complex_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_ovlp_complex_scalapack                  ')
         allocate(first_order_ham_complex_scalapack(mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_ham_complex_scalapack                  ')
         allocate(first_order_U_complex_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_U_complex_scalapack                  ')
         allocate(first_order_edm_complex_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_edm_complex_scalapack                  ')

         first_order_ovlp_complex_scalapack = 0.0d0
         first_order_ham_complex_scalapack = 0.0d0
         first_order_U_complex_scalapack = 0.0d0
         first_order_edm_complex_scalapack = 0.0d0
       endif

       if(use_DFPT_dielectric) then
         allocate(momentum_matrix_complex_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'momentum_matrix_complex_scalapack                  ')
         allocate(Omega_MO_complex_scalapack(mxld,mxcol,3),stat=info)
         call check_allocation(info, 'Omega_MO_complex_scalapack                  ')
         allocate(first_order_ham_complex_scalapack(mxld,mxcol,n_spin),stat=info)
         call check_allocation(info, 'first_order_ham_complex_scalapack                  ')
         allocate(first_order_U_complex_scalapack(mxld,mxcol),stat=info)
         call check_allocation(info, 'first_order_U_complex_scalapack                  ')

         momentum_matrix_complex_scalapack = 0.0d0
         Omega_MO_complex_scalapack = 0.0d0
         first_order_ham_complex_scalapack = 0.0d0
         first_order_U_complex_scalapack = 0.0d0

       endif

       ! Safety only:
       ovlp_complex     = 0
       ham_complex      = 0
       eigenvec_complex = 0

    endif


    ! Calculate workspace needed for eigenvalue solver

    if(real_eigenvectors)then

      lcwork = 0 ! no complex workspace needed

      np0 = NUMROC( MAX(n_basis,nb,2), nb, 0, 0, nprow )
      nq0 = NUMROC( MAX(n_basis,nb,2), nb, 0, 0, npcol )
      TRILWMIN = 3*n_basis + MAX( NB*( NP0+1 ), 3*NB )
      lwormtr = MAX( (NB*(NB-1))/2, (np0 + nq0)*NB + 2*NB*NB)
      lrwork = MAX( 1+6*n_basis+2*NP0*NQ0, TRILWMIN, lwormtr ) + 2*n_basis

      liwork = MAX(7*n_basis + 8*NPCOL + 2, n_basis + 2*NB + 2*npcol)

      if(use_elpa) then
         lrwork = 1
         liwork = 1
      endif

      len_scalapack_work = lrwork ! Total workspace (real numbers)

      write(info_str, *) ' Required Scalapack workspace - INTEGER: ',liwork, &
                         ' REAL:  ',lrwork
      call localorb_info(info_str, use_unit, '(A)')

    else

      np0 = NUMROC( MAX(n_basis,nb,2), nb, 0, 0, nprow )
      nq0 = NUMROC( MAX(n_basis,nb,2), nb, 0, 0, npcol )
      lcwork = n_basis + ( NP0 + NQ0 + NB ) * NB

      lrwork = 1 + 9*n_basis + 3*NP0*NQ0

      liwork = 7*n_basis + 8*NPCOL + 2

      if(use_elpa) then
         lcwork = 1
         lrwork = 1
         liwork = 1
      endif

      len_scalapack_work = 2*lcwork + lrwork ! Total workspace (real numbers)

      write(info_str, *) ' Required Scalapack workspace - INTEGER: ',liwork, &
                         ' REAL:  ',lrwork,' COMPLEX: ',lcwork
      call localorb_info(info_str, use_unit, '(A)')

    endif

    ! ELPA needs MPI communicators for communicating between rows/cols
    ! of the processor grid.
    ! Please note about the nomenclature:
    ! mpi_comm_cols is used for communication between processor columns
    ! (where processors talking to each other have the same row number).
    ! Analogous for mpi_comm_rows.

    if(my_scalapack_id < npcol*nprow) then
      call mpi_comm_split( my_scalapack_comm_work, myprow, my_scalapack_id, mpi_comm_cols, mpierr)
      call mpi_comm_split( my_scalapack_comm_work, mypcol, my_scalapack_id, mpi_comm_rows, mpierr)
    else
      mpi_comm_cols = mpi_comm_null
      mpi_comm_rows = mpi_comm_null
    endif

  end subroutine initialize_scalapack





!******
!-------------------------------------------------------------------------------
!****s* scalapack_wrapper/reinitialize_scalapack
!  NAME
!    reinitialize_scalapack
!  SYNOPSIS

  subroutine reinitialize_scalapack()

!  PURPOSE
!    Reinitialize the ScaLAPACK environment after relaxation step.
!  USES
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!   o all arrays for local index communication are deallocated
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    ! Deallocate all arrays for local index communication

    ! For sparse matrix communication:

    if(allocated(send_idx))       deallocate(send_idx)
    if(allocated(send_idx_count)) deallocate(send_idx_count)
    if(allocated(send_idx_displ)) deallocate(send_idx_displ)
    send_idx_count_tot = 0 ! safety only

    ! The send_ccc... arrays are already deallocated

    if(allocated(recv_row))       deallocate(recv_row)
    if(allocated(recv_row_count)) deallocate(recv_row_count)
    if(allocated(recv_row_displ)) deallocate(recv_row_displ)
    recv_row_count_tot = 0 ! safety only

    if(allocated(recv_ccc))       deallocate(recv_ccc)
    if(allocated(recv_ccc_count)) deallocate(recv_ccc_count)
    if(allocated(recv_ccc_displ)) deallocate(recv_ccc_displ)
    recv_ccc_count_tot = 0 ! safety only

    ! For full matrix communication:
    ! Please note that this is not really necessary since
    ! init_comm_full_local_matrix_scalapack deallocates
    ! arrays also when called a second time.

    if(allocated(basis_row)) deallocate(basis_row)
    if(allocated(basis_col)) deallocate(basis_col)
    if(allocated(basis_row_limit)) deallocate(basis_row_limit)
    if(allocated(basis_col_limit)) deallocate(basis_col_limit)
    if(allocated(send_mat_count)) deallocate(send_mat_count)
    if(allocated(send_mat_displ)) deallocate(send_mat_displ)
    if(allocated(recv_mat_count)) deallocate(recv_mat_count)
    if(allocated(recv_mat_displ)) deallocate(recv_mat_displ)
    send_mat_count_tot = 0
    recv_mat_count_tot = 0

    n_basis_local = 0
    n_local_matrix_size = 0
    if(allocated(i_basis_local)) deallocate(i_basis_local)

    full_ovlp_ready = .false.

  end subroutine reinitialize_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/scalapack_err_exit
!  NAME
!    scalapack_err_exit
!  SYNOPSIS
  subroutine scalapack_err_exit(info, name)
!  PURPOSE
!    Exits ScaLAPACK after an error.
!  USES
    use localorb_io
    use mpi_tasks
    implicit none
!  ARGUMENTS
    integer :: info
    character*(*) :: name
!  INPUTS
!    o info -- ScaLAPACK error code
!    o name -- name of the routine where error occurred
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    character*100 :: info_str
!    integer :: mpierr

    write (info_str,'(2X,A,I5,A,A)') 'Error ',info,' in ',name
    call localorb_info(info_str,use_unit,'(A)')

    call aims_stop

!    call MPI_Finalize(mpierr)
!    stop
  end subroutine scalapack_err_exit

  !-----------------------------------------------------------------------------------
  subroutine setup_scalapack_rmatrix( input_matrix, scalapack_matrix )
    real*8 :: input_matrix(:,:)
    real*8 :: scalapack_matrix(:,:)

    integer :: i_index, i_col, i_row, lr, lc

    scalapack_matrix(:,:) = 0

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    ! store the matrices
    i_index = 0
    do i_col = 1, n_basis, 1
      lc = l_col(i_col) ! local column number
      if(lc>0) then
        do i_row = 1, i_col, 1
          lr = l_row(i_row) ! local row number
          if(lr>0) then
!            scalapack_matrix(lr,lc) = input_matrix(i_index+i_row)
            scalapack_matrix(lr,lc) = input_matrix(i_row, i_col)
          endif
        enddo
      endif
      i_index = i_index + i_col
    end do

  end subroutine setup_scalapack_rmatrix

  !-----------------------------------------------------------------------------------
  subroutine scalapack_output_local_matrix(local_matrix, file_prefix)
    use localorb_io, only: use_unit
    use mpi_tasks, only: myid
    implicit none

    real*8 :: local_matrix(:,:)
    character(len=*) :: file_prefix

    integer :: index_1, index_2
    character(len=50) :: filename

    write(filename, '(A,A,I4.4,A)') file_prefix, '.', myid, '.out'
    write(use_unit,*) filename, mxld, mxcol
    open(50+myid, file=filename)
    do index_1 = 1, mxld, 1
      write(50+myid, '(5000E20.11)') (local_matrix(index_1,index_2), index_2 = 1, mxcol)
    end do
    close(50+myid)

  end subroutine scalapack_output_local_matrix

  !-----------------------------------------------------------------------------------
  subroutine scalapack_output_global_matrix(nmax, global_matrix, file_prefix)
    use mpi_tasks, only: myid

    real*8 :: global_matrix(:,:)
    integer :: nmax
    character(len=*) :: file_prefix

    integer :: index_1, index_2
    character(len=200) :: filename

    if (myid.ne.0) return

    write(filename, '(A,A)') trim(file_prefix), '.global.out'
    open(51, file=filename)
    do index_1 = 1, nmax, 1
      write(51, '(5000E20.11)') (global_matrix(index_1,index_2), index_2 = 1, nmax)
    end do
    close(51)

  end subroutine scalapack_output_global_matrix

  !-----------------------------------------------------------------------------------
  subroutine setup_scalapack_full_rmatrix( input_matrix, scalapack_matrix )
    real*8 :: input_matrix(:,:)
    real*8 :: scalapack_matrix(:,:)

    integer :: i_index, i_col, i_row, lr, lc

    scalapack_matrix(:,:) = 0

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    ! store the matrices
    i_index = 0
    do i_col = 1, n_basis, 1
      lc = l_col(i_col) ! local column number
      if(lc>0) then
        do i_row = 1, i_col, 1
          lr = l_row(i_row) ! local row number
          if(lr>0) then
!            scalapack_matrix(lr,lc) = input_matrix(i_index+i_row)
            scalapack_matrix(lr,lc) = input_matrix(i_row, i_col)
          endif
        enddo
      endif
      i_index = i_index + i_col
    end do

    call set_full_matrix_real(scalapack_matrix)

  end subroutine setup_scalapack_full_rmatrix

  !-----------------------------------------------------------------------------------
  subroutine setup_scalapack_full_zmatrix( input_matrix, scalapack_matrix )
    complex*16 :: input_matrix(:,:)
    complex*16 :: scalapack_matrix(:,:)

    integer :: i_col, i_row, lr, lc

    scalapack_matrix(:,:) = 0

    ! store the matrices
    do i_col = 1, n_basis, 1
      lc = l_col(i_col) ! local column number
      if(lc>0) then
        do i_row = 1, n_basis, 1
          lr = l_row(i_row) ! local row number
          if(lr>0) then
            scalapack_matrix(lr,lc) = input_matrix(i_row, i_col)
          endif
        enddo
      endif
    end do

  end subroutine setup_scalapack_full_zmatrix

  !-----------------------------------------------------------------------------------
  subroutine get_scalapack_global_rmatrix( local_matrix, global_matrix )
    use dimensions, only: n_hamiltonian_matrix_size_no_symmetry
    use synchronize_mpi_basic, only: sync_vector
    implicit none

    real*8 :: local_matrix(:,:)
    real*8 :: global_matrix(:,:)

    integer :: i_col, i_row

    global_matrix(:,:) = 0
    do i_col = 1, n_basis
      if(l_col(i_col)==0) cycle
      do i_row = 1, n_basis
        if(l_row(i_row)>0) then
          global_matrix(i_row,i_col) = local_matrix(l_row(i_row),l_col(i_col))
        endif
      end do
    end do
    call sync_vector(global_matrix, n_basis*n_basis, my_scalapack_comm_all)
  end subroutine get_scalapack_global_rmatrix

  !-----------------------------------------------------------------------------------
  subroutine get_scalapack_global_zmatrix( local_matrix, global_matrix )
    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none

    complex*16 :: local_matrix(:,:)
    complex*16 :: global_matrix(:,:)

    integer :: i_col, i_row

    global_matrix(:,:) = 0
    do i_col = 1, n_basis
      if(l_col(i_col)==0) cycle
      do i_row = 1, n_basis
        if(l_row(i_row)>0) then
          global_matrix(i_row,i_col) = local_matrix(l_row(i_row),l_col(i_col))
        endif
      end do
    end do
    call sync_vector_complex(global_matrix, n_basis*n_basis, &
          my_scalapack_comm_all)
  end subroutine get_scalapack_global_zmatrix


  !-----------------------------------------------------------------------------------
  subroutine scalapack_pdsyev( hamiltonian_predictor_step_scalapack,&
        predictor_eigenvalues, hamiltonian_predictor_ev)
    real*8 :: hamiltonian_predictor_step_scalapack(:,:)
    real*8 :: hamiltonian_predictor_ev(:,:)
    real*8 :: predictor_eigenvalues(:)

    integer :: lwork, info, i_col, i_row
    real*8, dimension(:), allocatable :: work

    allocate(work(1))

    call pdsyev( 'V', 'U', n_basis, hamiltonian_predictor_step_scalapack, &
         1, 1, sc_desc, predictor_eigenvalues, hamiltonian_predictor_ev, 1, 1, sc_desc, work, -1, info )
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    call pdsyev( 'V', 'U', n_basis, hamiltonian_predictor_step_scalapack, &
         1, 1, sc_desc, predictor_eigenvalues, hamiltonian_predictor_ev, 1, 1, sc_desc, WORK, lwork, info )

    deallocate(work)
  end subroutine scalapack_pdsyev


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/setup_hamiltoninan_scalapack
!  NAME
!    setup_hamiltoninan_scalapack
!  SYNOPSIS
  subroutine setup_hamiltonian_scalapack( hamiltonian_w, one_spin_only )
!  PURPOSE
!    Sets the Hamiltonian in the ScaLAPACK array.
!  USES
    implicit none
!  ARGUMENTS
    real*8:: hamiltonian_w   ( n_hamiltonian_matrix_size, n_spin )
    logical, optional :: one_spin_only
!  INPUTS
!    o hamiltonian_w -- the Hamiltonian
!    o one_spin_only -- if true only spin channel one is set
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: n_spin_max, i_spin, i_index, i_col, i_row, lr, lc

    n_spin_max = n_spin
    if(present(one_spin_only)) then
       if(one_spin_only) n_spin_max = 1
    endif

    ham(:,:,:) = 0

    ! Attention: Only the upper half of the matrix is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    do i_spin = 1, n_spin_max
       i_index = 0
       do i_col = 1, n_basis
          lc = l_col(i_col) ! local column number
          if(lc>0) then
             do i_row = 1, i_col, 1
                lr = l_row(i_row) ! local row number
                if(lr>0) then
                   ham (lr,lc,i_spin) = hamiltonian_w(i_index+i_row,i_spin)
                endif
             enddo
          endif
          i_index = i_index + i_col
       end do
    end do

  end subroutine setup_hamiltonian_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/setup_overlap_scalapack
!  NAME
!    setup_overlap_scalapack
!  SYNOPSIS
  subroutine setup_overlap_scalapack( overlap_matrix_w )
!  PURPOSE
!    Sets the overlap matrix in the ScaLAPACK array.
!  USES
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8:: overlap_matrix_w( n_hamiltonian_matrix_size )
!  INPUTS
!    o overlap_matrix_w -- overlap matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ovlp is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: i_index, i_col, i_row, lr, lc

    ovlp(:,:) = 0

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    ! store the matrices
    i_index = 0
    do i_col = 1, n_basis, 1
       lc = l_col(i_col) ! local column number
       if(lc>0) then
          do i_row = 1, i_col, 1
             lr = l_row(i_row) ! local row number
             if(lr>0) then
                ovlp(lr,lc) = overlap_matrix_w(i_index+i_row)
             endif
          enddo
       endif
       i_index = i_index + i_col
    end do

    factor_overlap = .TRUE.
    use_ovlp_trafo = .FALSE.
    full_ovlp_ready = .false.
    n_nonsing_ovlp = n_basis

    if (use_wf_extrapolation) call wf_save_overlap(mxld, mxcol, ovlp)

  end subroutine setup_overlap_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_hamiltonian_scalapack
!  NAME
!    construct_hamiltonian_scalapack
!  SYNOPSIS
  subroutine construct_hamiltonian_scalapack( hamiltonian )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use runtime_choices
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis

    implicit none
!  ARGUMENTS
    real*8:: hamiltonian   ( n_hamiltonian_matrix_size, n_spin )
!  INPUTS
!    o hamiltonian -- the Hamilton matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_scalapack + use_local_index")

    if(real_eigenvectors)then
       ham(:,:,:) = 0.
    else
       ham_complex(:,:,:) = 0.
    end if

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         ham (lr,lc,i_spin) = ham (lr,lc,i_spin) &
                                            + dble(k_phase(i_cell,my_k_point)) * hamiltonian(idx,i_spin)
                      else ! complex eigenvectors
                         ham_complex (lr,lc,i_spin) = ham_complex (lr,lc,i_spin) &
                                                    + k_phase(i_cell,my_k_point) * hamiltonian(idx,i_spin)
                      end if ! real_eigenvectors

                   end do

                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_hamiltonian_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format
    end do

  end subroutine construct_hamiltonian_scalapack


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_hamiltonian_real_for_elsi_scalapack
!  NAME
!    construct_hamiltonian_real_for_elsi_scalapack
!  SYNOPSIS
  subroutine construct_hamiltonian_real_for_elsi_scalapack( matrix, mat )
!  PURPOSE
!    This subroutine give the scalapack version (mat) for global sparse matrix (matrix) 
!  USES
    use runtime_choices
    use dimensions, only: n_periodic
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS
    real*8,     intent(in)  :: matrix         ( n_hamiltonian_matrix_size, n_spin )
    real*8,     intent(out) :: mat            ( mxld, mxcol, n_spin )
!  INPUTS
!    o matrix          -- The analogue to hamiltonian
!  OUTPUT
!    o mat             -- The analogue to ham
!    Only the upper halves of mat are set on exit.
!  AUTHOR
!    Forked by shanghui
!  HISTORY
!    Added in 2019, forked from construct_hamiltonian_scalapack
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_like_matrix_scalapack must not be used when use_local_index is set! (apparantly)

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_like_matrix_scalapack + use_local_index")

     mat(:,:,:) = 0.0d0

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                         mat (lr,lc,i_spin) = mat (lr,lc,i_spin) &
                                            + dble(k_phase(i_cell,my_k_point)) * matrix(idx,i_spin)

                   end do

                end if
             end do
          end do ! i_cell

       case (PM_none) !---------------------------------------------------------

          if (n_periodic .gt. 0) then
             write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support periodic non-packed matrices.'
             call aims_stop
          end if

          idx = 0
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number

             do i_row = 1, i_col

                   lr = l_row(i_row) ! local row number

                   idx = idx + 1
                   if(lc==0 .or. lr==0) cycle   ! skip if not local

                      mat (lr,lc,i_spin) = matrix (idx,i_spin)

             end do

          end do

       case default

          write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support the selected matrix packing.'
          call aims_stop

       end select ! packed_matrix_format
     enddo
  end subroutine construct_hamiltonian_real_for_elsi_scalapack



!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_hamiltonian_polar_for_elsi_scalapack
!  NAME
!    construct_first_order_hamiltonian_polar_for_elsi_scalapack
!  SYNOPSIS
  subroutine construct_first_order_hamiltonian_polar_for_elsi_scalapack( matrix, mat )
!  PURPOSE
!    This subroutine give the scalapack version (mat) for global sparse matrix (matrix) 
!  USES
    use runtime_choices
    use dimensions, only: n_periodic
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS
    real*8,     intent(in)  :: matrix         ( n_basis, n_basis, n_spin )
    real*8,     intent(out) :: mat            ( mxld, mxcol, n_spin )
!  INPUTS
!    o matrix          -- The analogue to hamiltonian
!  OUTPUT
!    o mat             -- The analogue to ham
!    Only the upper halves of mat are set on exit.
!  AUTHOR
!    Forked by shanghui
!  HISTORY
!    Added in 2019, forked from construct_hamiltonian_scalapack
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_like_matrix_scalapack must not be used when use_local_index is set! (apparantly)

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_like_matrix_scalapack + use_local_index")

     mat(:,:,:) = 0.0d0

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    do i_spin = 1, n_spin

          idx = 0
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number

             do i_row = 1, i_col

                   lr = l_row(i_row) ! local row number

                   if(lc==0 .or. lr==0) cycle   ! skip if not local

                      mat (lr,lc,i_spin) = matrix (i_row, i_col, i_spin)

             end do

          end do

     enddo ! i_spin
  end subroutine construct_first_order_hamiltonian_polar_for_elsi_scalapack


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_polar_for_elsi_scalapack
!  NAME
!    get_first_order_dm_polar_for_elsi_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_polar_for_elsi_scalapack( mat, matrix )
!  PURPOSE
!    Reconstructs the global matrix from ScaLAPACK
!  USES
    implicit none
!  ARGUMENTS
    real*8,     intent(in) :: mat  ( mxld, mxcol )
    real*8,     intent(out)  :: matrix  ( n_basis, n_basis )

!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer ::  lr, lc, i_col, i_row, i_coord

    character*200 :: info_str

    matrix(:,:)= 0.0d0

       do i_col = 1, n_basis
          lc = l_col(i_col) ! local column number
          if(lc>0) then
             do i_row = 1, n_basis
                lr = l_row(i_row) ! local row number
                if(lr>0) then
                   matrix(i_row,i_col) = mat(lr,lc)
                endif
             enddo ! i_row
          endif
       end do ! i_col


  end subroutine get_first_order_dm_polar_for_elsi_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_polar_reduce_memory_for_elsi_scalapack
!  NAME
!    get_first_order_dm_polar_reduce_memory_for_elsi_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_polar_reduce_memory_for_elsi_scalapack( mat, matrix_sparse )
!  PURPOSE
!    Reconstructs the global matrix from ScaLAPACK : the real part of
!    get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
!  USES
   use pbc_lists, only:  position_in_hamiltonian, n_cells_in_hamiltonian,  &
                         column_index_hamiltonian,index_hamiltonian
    implicit none
!  ARGUMENTS
    real*8,     intent(in)   :: mat( mxld, mxcol )
    real*8,     intent(out)  :: matrix_sparse( n_hamiltonian_matrix_size )
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o first_order_density_matrix_sparse -- set to the contents of first_order_ham_polar_scalapack
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_spin
    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    matrix_sparse= 0.0d0
     

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) ! 
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

                matrix_sparse(i_index) = mat(lr,lc)

          end do
       end do
    end do


  end subroutine get_first_order_dm_polar_reduce_memory_for_elsi_scalapack



!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_overlap_real_for_elsi_scalapack
!  NAME
!    construct_overlap_real_for_elsi_scalapack
!  SYNOPSIS
  subroutine construct_overlap_real_for_elsi_scalapack( matrix, mat )
!  PURPOSE
!    This subroutine give the scalapack version (mat) for global sparse matrix (matrix), no spin for overlap
!  USES
    use runtime_choices
    use dimensions, only: n_periodic
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS
    real*8,     intent(in)  :: matrix         ( n_hamiltonian_matrix_size )
    real*8,     intent(out) :: mat            ( mxld, mxcol )
!  INPUTS
!    o matrix          -- The analogue to hamiltonian
!  OUTPUT
!    o mat             -- The analogue to ham
!    Only the upper halves of mat are set on exit.
!  AUTHOR
!    Forked by shanghui
!  HISTORY
!    Added in 2019, forked from construct_hamiltonian_scalapack
!  SOURCE


    integer::  i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_like_matrix_scalapack must not be used when use_local_index is set! (apparantly)

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_like_matrix_scalapack + use_local_index")

     mat(:,:) = 0.0d0

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!


       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                         mat (lr,lc) = mat (lr,lc) &
                                            + dble(k_phase(i_cell,my_k_point)) * matrix(idx)

                   end do

                end if
             end do
          end do ! i_cell

       case (PM_none) !---------------------------------------------------------

          if (n_periodic .gt. 0) then
             write(use_unit,*) 'Error: construct_overlap_real_for_elsi_scalapack does not support periodic non-packed matrices.'
             call aims_stop
          end if

          idx = 0
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number

             do i_row = 1, i_col

                   lr = l_row(i_row) ! local row number

                   idx = idx + 1
                   if(lc==0 .or. lr==0) cycle   ! skip if not local

                      mat (lr,lc) = matrix (idx)

             end do

          end do

       case default

          write(use_unit,*) 'Error: construct_overlap_real_for_elsi_scalapack does not support the selected matrix packing.'
          call aims_stop

       end select ! packed_matrix_format
  end subroutine construct_overlap_real_for_elsi_scalapack



!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_hamiltonian_like_matrix_scalapack
!  NAME
!    construct_hamiltonian_like_matrix_scalapack
!  SYNOPSIS
  subroutine construct_hamiltonian_like_matrix_scalapack( matrix, mat, mat_complex )
!  PURPOSE
!    This subroutine is identical to construct_hamiltonian_scalapack save one key
!    difference; it passes its outputs as variables instead of modifying module
!    variables, allowing this subroutine to be used for generic matrices with the same
!    dimensions and packing as the Hamiltonian (for example, matrix elements for weak
!    perturbations to the Hamiltonian.)
!  USES
    use runtime_choices
    use dimensions, only: n_periodic
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS
    real*8,     intent(in)  :: matrix         ( n_hamiltonian_matrix_size, n_spin )
    real*8,     intent(out) :: mat            ( mxld, mxcol, n_spin )
    complex*16, intent(out) :: mat_complex    ( mxld, mxcol, n_spin )
!  INPUTS
!    o matrix          -- The analogue to hamiltonian
!  OUTPUT
!    o mat             -- The analogue to ham
!    o mat_complex     -- The analogue to ham_complex
!    Only the upper halves of mat(_complex) are set on exit.
!  AUTHOR
!    Forked by William Huhn
!  HISTORY
!    Added in 2016, forked from construct_hamiltonian_scalapack
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_like_matrix_scalapack must not be used when use_local_index is set! (apparantly)

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_like_matrix_scalapack + use_local_index")

    if(real_eigenvectors)then
       mat(:,:,:) = 0.0d0
    else
       mat_complex(:,:,:) = (0.0d0,0.0d0)
    end if

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         mat (lr,lc,i_spin) = mat (lr,lc,i_spin) &
                                            + dble(k_phase(i_cell,my_k_point)) * matrix(idx,i_spin)
                      else ! complex eigenvectors
                         mat_complex (lr,lc,i_spin) = mat_complex (lr,lc,i_spin) &
                                                    + k_phase(i_cell,my_k_point) * matrix(idx,i_spin)
                      end if ! real_eigenvectors

                   end do

                end if
             end do
          end do ! i_cell

       case (PM_none) !---------------------------------------------------------

          if (n_periodic .gt. 0) then
             write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support periodic non-packed matrices.'
             call aims_stop
          end if

          idx = 0
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number

             do i_row = 1, i_col

                   lr = l_row(i_row) ! local row number

                   idx = idx + 1
                   if(lc==0 .or. lr==0) cycle   ! skip if not local

                   if(real_eigenvectors)then
                      mat (lr,lc,i_spin) = matrix (idx,i_spin)
                   else ! complex eigenvectors
                      mat_complex (lr,lc,i_spin) = matrix (idx,i_spin)
                   end if ! real_eigenvectors

             end do

          end do

       case default

          write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support the selected matrix packing.'
          call aims_stop

       end select ! packed_matrix_format
    end do
  end subroutine construct_hamiltonian_like_matrix_scalapack


  subroutine construct_hamiltonian_like_matrix_zero_diag_scalapack( matrix, mat, mat_complex )
!  PURPOSE
!    This subroutine is identical to construct_hamiltonian_scalapack save two key
!    differences; it passes its outputs as variables instead of modifying module
!    variables, allowing this subroutine to be used for generic matrices with the same
!    dimensions and packing as the Hamiltonian (for example, matrix elements for weak
!    perturbations to the Hamiltonian) AND it zeros the diagonal.
!  USES
    use runtime_choices
    use dimensions, only: n_periodic
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS
    real*8,     intent(in)  :: matrix         ( n_hamiltonian_matrix_size, n_spin )
    real*8,     intent(out) :: mat            ( mxld, mxcol, n_spin )
    complex*16, intent(out) :: mat_complex    ( mxld, mxcol, n_spin )
!  INPUTS
!    o matrix          -- The analogue to hamiltonian
!  OUTPUT
!    o mat             -- The analogue to ham
!    o mat_complex     -- The analogue to ham_complex
!    Only the upper halves of mat(_complex) are set on exit.
!  AUTHOR
!    Forked by William Huhn
!  HISTORY
!    Added in 2016, forked from construct_hamiltonian_scalapack
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_like_matrix_scalapack must not be used when use_local_index is set! (apparantly)

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_like_matrix_scalapack + use_local_index")

    if(real_eigenvectors)then
       mat(:,:,:) = 0.
    else
       mat_complex(:,:,:) = (0.0d0, 0.0d0)
    end if

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         if (i_col == i_row) then
                            mat (lr,lc,i_spin) = 0.0d0
                         else
                         mat (lr,lc,i_spin) = mat (lr,lc,i_spin) &
                                            + dble(k_phase(i_cell,my_k_point)) * matrix(idx,i_spin)
                         endif
                      else ! complex eigenvectors
                         if (i_col == i_row) then
                            mat_complex (lr,lc,i_spin) = (0.0d0,0.0d0)
                         else
                         mat_complex (lr,lc,i_spin) = mat_complex (lr,lc,i_spin) &
                                                    + k_phase(i_cell,my_k_point) * matrix(idx,i_spin)
                         endif
                      end if ! real_eigenvectors

                   end do

                end if
             end do
          end do ! i_cell

       case (PM_none) !---------------------------------------------------------

          if (n_periodic .gt. 0) then
             write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support periodic non-packed matrices.'
             call aims_stop
          end if

          idx = 0
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number

             do i_row = 1, i_col

                   lr = l_row(i_row) ! local row number

                   idx = idx + 1
                   if(lc==0 .or. lr==0) cycle   ! skip if not local

                   if(real_eigenvectors)then
                      mat (lr,lc,i_spin) = matrix (idx,i_spin)
                   else ! complex eigenvectors
                      mat_complex (lr,lc,i_spin) = matrix (idx,i_spin)
                   end if ! real_eigenvectors

             end do

          end do

       case default

          write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support the selected matrix packing.'
          call aims_stop

       end select ! packed_matrix_format
    end do
  end subroutine construct_hamiltonian_like_matrix_zero_diag_scalapack



 !******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_hamiltonian_like_matrix_scalapack
!  NAME
!    construct_hamiltonian_like_matrix_scalapack
!  SYNOPSIS
  subroutine construct_overlap_like_matrix_scalapack_die( matrix_up, matrix_low, mat, mat_complex )
!  PURPOSE
!    This subroutine is identical to construct_hamiltonian_scalapack save one key
!    difference; it passes its outputs as variables instead of modifying module
!    variables, allowing this subroutine to be used for generic matrices with the same
!    dimensions and packing as the Hamiltonian (for example, matrix elements for weak
!    perturbations to the Hamiltonian.)
!  USES
    use runtime_choices
    use dimensions, only: n_periodic
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS
    real*8,     intent(in)  :: matrix_up         ( n_hamiltonian_matrix_size )
    real*8,     intent(in)  :: matrix_low        ( n_hamiltonian_matrix_size )
    real*8,     intent(out) :: mat            ( mxld, mxcol )
    complex*16, intent(out) :: mat_complex    ( mxld, mxcol )
!  INPUTS
!    o matrix          -- The analogue to hamiltonian
!  OUTPUT
!    o mat             -- The analogue to ham
!    o mat_complex     -- The analogue to ham_complex
!    Only the upper halves of mat(_complex) are set on exit.
!  AUTHOR
!    Forked by William Huhn
!  HISTORY
!    Added in 2016, forked from construct_hamiltonian_scalapack
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_like_matrix_scalapack must not be used when use_local_index is set! (apparantly)

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_like_matrix_scalapack + use_local_index")

    if(real_eigenvectors)then
       mat(:,:) = 0.
    else
       mat_complex(:,:) = 0.
    end if

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!


       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                      if(real_eigenvectors)then

                         if (i_row .eq. i_col) then
                           mat (lr,lc) = mat (lr,lc) &
                                            + dble(k_phase(i_cell,my_k_point)) * matrix_up(idx)
                         else
                           mat (lr,lc) = mat (lr,lc) &
                                            + dble(k_phase(i_cell,my_k_point))* matrix_up(idx)
                           mat (lc,lr) = mat (lc,lr) &
                                            + dble(k_phase(i_cell,my_k_point))*matrix_low(idx)
                        endif

                      else ! complex eigenvectors
                         if (i_row .eq. i_col) then
                            mat_complex (lr,lc) = mat_complex (lr,lc) &
                                                    + k_phase(i_cell,my_k_point) * matrix_up(idx)
                         else
                            mat_complex (lr,lc) = mat_complex (lr,lc) &
                                                    + k_phase(i_cell,my_k_point) * matrix_up(idx)
                            mat_complex (lc,lr) = mat_complex (lc,lr) &
                                                    +conjg(k_phase(i_cell,my_k_point)* matrix_low(idx))
                         endif

                      end if ! real_eigenvectors

                   end do

                end if
             end do
          end do ! i_cell

       case (PM_none) !---------------------------------------------------------

          if (n_periodic .gt. 0) then
             write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support periodic non-packed matrices.'
             call aims_stop
          end if

          idx = 0
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number

             do i_row = 1, i_col

                   lr = l_row(i_row) ! local row number

                   idx = idx + 1
                   if(lc==0 .or. lr==0) cycle   ! skip if not local

                   if(real_eigenvectors)then
                      mat (lr,lc) = matrix_up (idx)
                   else ! complex eigenvectors
                      mat_complex (lr,lc) = matrix_low (idx)
                   end if ! real_eigenvectors

             end do

          end do

       case default

          write(use_unit,*) 'Error: construct_hamiltonian_like_matrix_scalapack does not support the selected matrix packing.'
          call aims_stop

       end select ! packed_matrix_format
  end subroutine construct_overlap_like_matrix_scalapack_die
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_overlap_scalapack
!  NAME
!    construct_overlap_scalapack
!  SYNOPSIS
  subroutine construct_overlap_scalapack( overlap_matrix )
!  PURPOSE
!    Sets the overlap matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use runtime_choices
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis

    implicit none
!  ARGUMENTS
    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
!  INPUTS
!    o overlap_matrix -- the overlap matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_cell, i_col, i_row, lr, lc, idx

    ! construct_overlap_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_overlap_scalapack + use_local_index")


    if(real_eigenvectors)then
       ovlp = 0.
    else
       ovlp_complex = 0.
    end if

    ! Attention: Only the upper half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'U'!

    select case(packed_matrix_format)

    case(PM_index) !------------------------------------------------

       do i_cell = 1, n_cells_in_hamiltonian-1
          do i_col = 1, n_basis

             lc = l_col(i_col) ! local column number
             if(lc==0) cycle   ! skip if not local

             if(index_hamiltonian(1,i_cell,i_col) > 0) then

                do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                   i_row = column_index_hamiltonian(idx)
                   lr = l_row(i_row) ! local row number
                   if(lr==0) cycle   ! skip if not local

                   if(real_eigenvectors)then
                      ovlp(lr,lc) = ovlp(lr,lc) + dble(k_phase(i_cell,my_k_point)) * overlap_matrix(idx)
                   else ! complex eigenvectors
                      ovlp_complex(lr,lc) = ovlp_complex(lr,lc) + k_phase(i_cell,my_k_point) * overlap_matrix(idx)
                   end if ! real_eigenvectors

                end do

             end if
          end do
       end do ! i_cell

    case default !---------------------------------------------------------

       write(use_unit,*) 'Error: construct_overlap_scalapack does not support non-packed matrices.'
       call aims_stop

    end select ! packed_matrix_format

    factor_overlap = .TRUE.
    use_ovlp_trafo = .FALSE.
    full_ovlp_ready = .false.
    n_nonsing_ovlp = n_basis

    if (use_wf_extrapolation) then
       if (real_eigenvectors) then
          call wf_save_overlap(mxld, mxcol, ovlp)
       else
          call wf_save_overlap_cmplx(mxld, mxcol, ovlp_complex)
       end if
    end if
  end subroutine construct_overlap_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/save_overlap_scalapack
!  NAME
!    save_overlap_scalapack
!  SYNOPSIS
  subroutine save_overlap_scalapack
!  PURPOSE
!    Saves the overlap matrix in the case it is needed for postprocessing
!  USES
    use aims_memory_tracking, only: aims_allocate
    implicit none

    if(real_eigenvectors)then
       if(.not. allocated(ovlp_stored)) then
          call aims_allocate(ovlp_stored,mxld,mxcol,"ovlp_stored")
       end if

       ovlp_stored = ovlp

       call set_full_matrix_real(ovlp_stored)
    else
       if(.not. allocated(ovlp_complex_stored)) then
          call aims_allocate(ovlp_complex_stored,mxld,mxcol,&
               "ovlp_complex_stored")
       end if

       ovlp_complex_stored = ovlp_complex

       call set_full_matrix_complex(ovlp_complex_stored)
    end if

  end subroutine save_overlap_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/save_ham_scalapack
!  NAME
!    save_ham_scalapack
!  SYNOPSIS
  subroutine save_ham_scalapack
!  PURPOSE
!    Saves the ham matrix in the case it is needed for
!    elsi_dm_polar_reduce_memory_cpscf [wyj add]
!  USES
    use aims_memory_tracking, only: aims_allocate
    implicit none

    if(real_eigenvectors)then
       if(.not. allocated(ham_stored)) then
          call aims_allocate(ham_stored,mxld,mxcol,n_spin,"ham_stored")
       end if

       ham_stored = ham

       !call set_full_matrix_real(ham_stored)
    else
       if(.not. allocated(ham_complex_stored)) then
          call aims_allocate(ham_complex_stored,mxld,mxcol,n_spin,&
               "ham_complex_stored")
       end if

       ham_complex_stored = ham_complex

       !call set_full_matrix_complex(ham_complex_stored)
    end if

  end subroutine save_ham_scalapack
!******
!-----------------------------------------------------------------------------------
! The following routines are only wrapper functions which should be used
! instead of the get_set_... versions
!-----------------------------------------------------------------------------------
  subroutine set_sparse_local_ovlp_scalapack(matrix)
    use runtime_choices
    implicit none

    real*8:: matrix(n_hamiltonian_matrix_size)

    call get_set_sparse_local_matrix_scalapack(matrix, 0, 1)

    if (use_wf_extrapolation) then
       if (real_eigenvectors) then
          call wf_save_overlap(mxld, mxcol, ovlp)
       else
          call wf_save_overlap_cmplx(mxld, mxcol, ovlp_complex)
       end if
    end if

  end subroutine set_sparse_local_ovlp_scalapack
!-----------------------------------------------------------------------------------
  subroutine set_sparse_local_ham_scalapack(matrix)
    implicit none

    real*8:: matrix(n_hamiltonian_matrix_size,n_spin)

    integer i_spin

    do i_spin = 1, n_spin
      call get_set_sparse_local_matrix_scalapack(matrix(1,i_spin), 1, i_spin)
    enddo

  end subroutine set_sparse_local_ham_scalapack
!-----------------------------------------------------------------------------------
  subroutine get_sparse_local_matrix_scalapack(matrix, i_spin)
    implicit none

    real*8:: matrix(n_hamiltonian_matrix_size)
    integer i_spin

    call get_set_sparse_local_matrix_scalapack(matrix, 2, i_spin)

  end subroutine get_sparse_local_matrix_scalapack
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_set_full_local_matrix_scalapack
!  NAME
!    get_set_full_local_matrix_scalapack
!  SYNOPSIS
  subroutine get_set_full_local_matrix_scalapack( matrix, which, i_spin )
!  PURPOSE
!    Gets or sets a ScaLAPACK matrix in local-index mode (working with full local matrices)
!    which = 0: Set ovlp from local overlap matrix (in matrix)
!    which = 1: Set ham from local hamiltonian (in matrix)
!    which = 2: Set matrix from scalapack ham (ham must be set to something usefull before call)
!
!    WPH:  To understand this subroutine, it is recommend that you first read
!    and absorb the init_comm_full_local_matrix_scalapack subroutine, which sets
!    up the communication for this subroutine.  It is in that subroutine where
!    most of the actual logic is performed to determine which chunks of data
!    goes where and the necessary data structures are defined and set up. This
!    subroutine then performs the heavy-and-stupid lifting of shuffling the data
!    and multiplying by k-point phases when needed.
!  USES
    use dimensions, only: n_periodic
    use localorb_io, only: localorb_info
    use pbc_lists, only: cbasis_to_basis, center_to_cell, cbasis_to_center, &
        position_in_hamiltonian, n_cells_in_hamiltonian, k_phase
    use mpi_tasks, only: n_tasks, mpi_wtime, mpi_logical, mpi_land, &
        mpi_comm_global, mpi_real8, myid, mpi_status_ignore, aims_stop
    use runtime_choices, only: use_alltoall, output_level
    implicit none
!  ARGUMENTS
    real*8 :: matrix(n_local_matrix_size)
    integer :: which, i_spin
!  INPUTS
!    o matrix -- the matrix to set/get
!    o which -- the operation to perform
!    o i_spin -- the spin channel
!  OUTPUT
!    if which = 0/1 upper half of the ScaLAPACK array ham/ovlp is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i, j, ip, i_k_point, k_nprow, k_npcol, ip_c, i_cell_r, i_cell_c
    integer :: i_cell, i_row, i_col, i_diff, i_send, i_recv, i_cnt_mat
    integer :: jj, n_off, n_rows_local, j_row_local(n_basis_local)
    integer :: lc, lr, istat1, istat2, mpierr
    integer :: send_mat_off(0:n_tasks-1), ip_r(n_basis_local)
    real*8, allocatable :: matrix_send(:), matrix_recv(:)
    logical alloc_success, use_alltoall_really
    real*8 ttt0

    character(*), parameter :: func = "get_set_full_local_matrix_scalapack"

    ttt0 = mpi_wtime()

    ! Reset the matrix to be set, check parameter which

    if(which==0) then
       if(real_eigenvectors)then
          ovlp(:,:) = 0
       else
          ovlp_complex(:,:) = 0
       end if
    else if(which==1) then
       if(real_eigenvectors)then
          ham(:,:,i_spin) = 0
       else
          ham_complex(:,:,i_spin) = 0
       end if
    else if(which==2) then
       matrix(:) = 0
    else
       call aims_stop("Illegal parameter for 'which'", func)
    endif

    use_alltoall_really = use_alltoall

    ! When use_alltoall is set we check if we can allocate the memory needed
    ! and if this doesn't succeed we fall back to sendrecv.

    ! The reason why we don't do that always this way but use a config switch is the fact that
    ! the system might be configured with swap space. In this case the allocate may succeed
    ! but the machine will start swapping - which is absolutely nothing we want!

    if(use_alltoall) then

      ! Try to allocate large matrices for using mpi_alltoallv ...
      allocate(matrix_send(send_mat_count_tot), stat=istat1)
      allocate(matrix_recv(recv_mat_count_tot), stat=istat2)

      ! ... check if allocation succeeded on ALL procs, otherways use mpi_sendrecv

      alloc_success = (istat1==0 .and. istat2==0)
      call mpi_allreduce(alloc_success, use_alltoall_really, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_global, mpierr)

      if(.not.use_alltoall_really) then
        ! fall back to sendrecv
        if(allocated(matrix_send)) deallocate(matrix_send)
        if(allocated(matrix_recv)) deallocate(matrix_recv)
        call localorb_info('  *** Not enough memory for using mpi_alltoall, falling back to mpi_sendrecv')
      endif

    endif

    if(.not.use_alltoall_really) then
       allocate(matrix_send(maxval(send_mat_count)))
       allocate(matrix_recv(maxval(recv_mat_count)))
    endif


    if(which<=1 .and. use_alltoall_really) then
      ! Send the complete local matrix to the receivers using 1 mpi_alltoallv call
      ! This is the most effective way but needs a large amount of memory

      ! Put matrix into matrix_send ...
      send_mat_off(:) = send_mat_displ(:)
      do i_k_point = 1, n_k_points

        k_nprow = k_point_desc(i_k_point)%nprow
        k_npcol = k_point_desc(i_k_point)%npcol

        do j=1,n_basis_local
          i_row = i_basis_local(j)
          ip_r(j) = mod((Cbasis_to_basis(i_row)-1)/nb, k_nprow) ! Processor row owning i_row
        enddo

        do i=1,n_basis_local
          i_col = i_basis_local(i)
          ip_c = mod((Cbasis_to_basis(i_col)-1)/nb, k_npcol) ! Processor column owning i_col
          do j=1,n_basis_local
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              ip = k_point_desc(i_k_point)%global_id(ip_r(j),ip_c)
              send_mat_off(ip) = send_mat_off(ip)+1
              if(j<=i) then
                matrix_send(send_mat_off(ip)) = matrix(i*(i-1)/2+j)
              else
                matrix_send(send_mat_off(ip)) = matrix(j*(j-1)/2+i)
              endif
            endif
          enddo
        enddo

      enddo

      ! ... and send it away
      call mpi_alltoallv(matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
                         matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
                         mpi_comm_global, mpierr)
    endif

    if(which>1 .and. use_alltoall_really) matrix_recv(:) = 0

    ! Insert data from local matrix into scalapack matrix (which<=1)
    ! or gather data from scalapack matrix for putting into local matrix (which>1)
    ! This is done for every remote task separatly for the case that
    ! mpi_alltoallv cannot be used

    do i_diff = 0, n_tasks-1

      i_send = mod(myid+i_diff,n_tasks)         ! Task to which we send data
      i_recv = mod(myid+n_tasks-i_diff,n_tasks) ! Task from which we get data

      if(.not.use_alltoall_really) then

        ! Gather all local rows going to proc i_send
        ! (for which > 1 this is needed below)
        n_rows_local = 0
        do j=1,n_basis_local
          i_row = i_basis_local(j)
          if(mod((Cbasis_to_basis(i_row)-1)/nb, gl_nprow(i_send)) == gl_prow(i_send)) then
            n_rows_local = n_rows_local+1
            j_row_local(n_rows_local) = j
          endif
        enddo

        if(which<=1) then

          n_off = 0
          do i=1,n_basis_local
            i_col = i_basis_local(i)
            if(mod((Cbasis_to_basis(i_col)-1)/nb, gl_npcol(i_send)) /= gl_pcol(i_send)) cycle
            do jj=1,n_rows_local
              j = j_row_local(jj)
              i_row = i_basis_local(j)
              if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
                n_off = n_off+1
                if(j<=i) then
                  matrix_send(n_off) = matrix(i*(i-1)/2+j)
                else
                  matrix_send(n_off) = matrix(j*(j-1)/2+i)
                endif
              endif
            enddo
          enddo

          ! Gather and send data for remote task i_send,
          ! receive corresponding data from task i_recv
          call mpi_sendrecv(matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
                            matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
                            mpi_comm_global, mpi_status_ignore, mpierr)

        else

          matrix_recv(1:recv_mat_count(i_recv)) = 0

        endif

      endif

      if(use_alltoall_really) then
        i_cnt_mat = recv_mat_displ(i_recv) ! Counter in matrix_recv
      else
        i_cnt_mat = 0 ! Counter in matrix_recv, reset for every task
      endif

      if(n_periodic==0) then

        ! We could also use the code for the periodic case here,
        ! but the code below is much more simple (and hopefully faster!)
        do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

          i_col = basis_col(j)
          lc = l_col(i_col)

          do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

            i_row = basis_row(i)
            if(i_row > i_col) exit ! done with this column
            lr = l_row(i_row)
            i_cnt_mat = i_cnt_mat+1

            if(which==0) then
              ovlp(lr,lc) = ovlp(lr,lc) + matrix_recv(i_cnt_mat)
            else if(which==1) then
              ham(lr,lc,i_spin) = ham(lr,lc,i_spin) + matrix_recv(i_cnt_mat)
            else
              matrix_recv(i_cnt_mat) = ham(lr,lc,i_spin)
            endif
          enddo
        enddo

      else ! periodic case

        do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

          i_col = basis_col(j)
          i_cell_c = center_to_cell(Cbasis_to_center(i_col))
          lc = l_col(Cbasis_to_basis(i_col))

          do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

            i_row = basis_row(i)
            if(Cbasis_to_basis(i_row) > Cbasis_to_basis(i_col)) cycle

            i_cell_r = center_to_cell(Cbasis_to_center(i_row))
            lr = l_row(Cbasis_to_basis(i_row))

            i_cell = position_in_hamiltonian(i_cell_c, i_cell_r) ! Attention: position_in_hamiltonian is not symmetric!

            i_cnt_mat = i_cnt_mat+1

            if(i_cell == n_cells_in_hamiltonian) cycle

            if(which==0) then
              if(real_eigenvectors)then
                ovlp(lr,lc) = ovlp(lr,lc) + dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
              else
                ovlp_complex(lr,lc) = ovlp_complex(lr,lc) + k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
              end if
            else if(which==1) then
              if(real_eigenvectors)then
                ham(lr,lc,i_spin) = ham(lr,lc,i_spin) + dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
              else
                ham_complex(lr,lc,i_spin) = ham_complex(lr,lc,i_spin) + k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
              end if
            else
              if(real_eigenvectors)then
                matrix_recv(i_cnt_mat) = ham(lr,lc,i_spin)*dble(k_phase(i_cell,my_k_point))
              else
                matrix_recv(i_cnt_mat) = dble(ham_complex(lr,lc,i_spin)*dconjg(k_phase(i_cell,my_k_point)))
              end if
            endif
          enddo
        enddo

      endif

      if(which>1 .and. .not.use_alltoall_really) then

        ! Send matrix_recv immediatly back to owner of local matrix
        call mpi_sendrecv(matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
                          matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
                          mpi_comm_global, mpi_status_ignore, mpierr)

        n_off = 0
        do i=1,n_basis_local
          i_col = i_basis_local(i)
          if(mod((Cbasis_to_basis(i_col)-1)/nb, gl_npcol(i_send)) /= gl_pcol(i_send)) cycle
          do jj=1,n_rows_local
            j = j_row_local(jj)
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              n_off = n_off+1
              ! some of the elements we got are duplicates, we must be careful
              ! not to add them twice!!!!
              if(Cbasis_to_basis(i_row)==Cbasis_to_basis(i_col) .and. j>i) cycle
              if(j<=i) then
                matrix(i*(i-1)/2+j) = matrix(i*(i-1)/2+j) + matrix_send(n_off)
              else
                matrix(j*(j-1)/2+i) = matrix(j*(j-1)/2+i) + matrix_send(n_off)
              endif
            endif
          enddo
        enddo

      endif

    enddo

    if(which>1 .and. use_alltoall_really) then
      ! Send the matrix gathered in matrix_recv to owners of local matrix with mpi_alltoallv
      ! This is the most effective way but needs a large amount of memory
      call mpi_alltoallv(matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
                         matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
                         mpi_comm_global, mpierr)

      ! Insert received matrix
      send_mat_off(:) = send_mat_displ(:)
      do i_k_point = 1, n_k_points

        k_nprow = k_point_desc(i_k_point)%nprow
        k_npcol = k_point_desc(i_k_point)%npcol

        do j=1,n_basis_local
          i_row = i_basis_local(j)
          ip_r(j) = mod((Cbasis_to_basis(i_row)-1)/nb, k_nprow) ! Processor row owning i_row
        enddo

        do i=1,n_basis_local
          i_col = i_basis_local(i)
          ip_c = mod((Cbasis_to_basis(i_col)-1)/nb, k_npcol) ! Processor column owning i_col
          do j=1,n_basis_local
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              ip = k_point_desc(i_k_point)%global_id(ip_r(j),ip_c)
              send_mat_off(ip) = send_mat_off(ip)+1
              ! some of the elements we got are duplicates, we must be careful
              ! not to add them twice!!!!
              if(Cbasis_to_basis(i_row)==Cbasis_to_basis(i_col) .and. j>i) cycle
              if(j<=i) then
                matrix(i*(i-1)/2+j) = matrix(i*(i-1)/2+j) + matrix_send(send_mat_off(ip))
              else
                matrix(j*(j-1)/2+i) = matrix(j*(j-1)/2+i) + matrix_send(send_mat_off(ip))
              endif
            endif
          enddo
        enddo

      enddo
    endif

    deallocate(matrix_send)
    deallocate(matrix_recv)

    if(which==0) then
       factor_overlap = .TRUE.
       use_ovlp_trafo = .FALSE.
       full_ovlp_ready = .false.
       n_nonsing_ovlp = n_basis
    endif

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
       write(use_unit,"(2X,A,F13.6,A)") "| Time get_set_full_local_matrix_scalapack:",mpi_wtime()-ttt0," s"

  end subroutine get_set_full_local_matrix_scalapack
!****s* scalapack_wrapper/set_full_local_matrix_scalapack_generic
!  NAME
!    set_full_local_matrix_scalapack_generic
!  SYNOPSIS
  subroutine set_full_local_matrix_scalapack_generic( matrix, loc_mat, loc_mat_complex )
!  PURPOSE
!    Sets a ScaLAPACK matrix in local-index mode (working with full local matrices)
!    This is the relevant code path when using both local index and load balancing
!    This subroutine is a fork of get_set_full_local_matrix_scalapack (which == 1);
!    the difference is that the desired ScaLAPACK matrix is passed as an argument
!    rather than pulled from a module variable, allowing for reuse by generic
!    matrices sharing the same packing (and ScaLAPACK descriptor) as the Hamiltonian
!    or overlap matrices
!    We deliberately do not set the incoming local matrices to zero, to allow
!    for accumulation
!  USES
    use dimensions, only: n_periodic
    use localorb_io, only: localorb_info
    use mpi_tasks
    use pbc_lists
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8,     intent(in)    :: matrix(n_local_matrix_size)
    real*8,     intent(inout) :: loc_mat(mxld,mxcol)
    complex*16, intent(inout) :: loc_mat_complex(mxld,mxcol)
!  INPUTS
!    o matrix -- the matrix, in local indexing
!  OUTPUT
!    The upper half of the appropriate ScaLAPACK array are set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    September 2017 - Forked off of get_set_full_local_matrix_scalapack by
!                     William Huhn (Duke University)
!  SOURCE

    integer :: i, j, ip, i_k_point, k_nprow, k_npcol, ip_c, i_cell_r, i_cell_c
    integer :: i_cell, i_row, i_col, i_diff, i_send, i_recv, i_cnt_mat
    integer :: jj, n_off, n_rows_local, j_row_local(n_basis_local)
    integer :: lc, lr, istat1, istat2, mpierr
    integer :: send_mat_off(0:n_tasks-1), ip_r(n_basis_local)
    real*8, allocatable :: matrix_send(:), matrix_recv(:)
    logical alloc_success, use_alltoall_really
real*8 ttt0
ttt0 = mpi_wtime()

    use_alltoall_really = use_alltoall

    ! When use_alltoall is set we check if we can allocate the memory needed
    ! and if this doesn't succeed we fall back to sendrecv.

    ! The reason why we don't do that always this way but use a config switch is the fact that
    ! the system might be configured with swap space. In this case the allocate may succeed
    ! but the machine will start swapping - which is absolutely nothing we want!

    if(use_alltoall) then

      ! Try to allocate large matrices for using mpi_alltoallv ...
      allocate(matrix_send(send_mat_count_tot), stat=istat1)
      allocate(matrix_recv(recv_mat_count_tot), stat=istat2)

      ! ... check if allocation succeeded on ALL procs, otherways use mpi_sendrecv

      alloc_success = (istat1==0 .and. istat2==0)
      call mpi_allreduce(alloc_success, use_alltoall_really, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_global, mpierr)

      if(.not.use_alltoall_really) then
        ! fall back to sendrecv
        if(allocated(matrix_send)) deallocate(matrix_send)
        if(allocated(matrix_recv)) deallocate(matrix_recv)
        call localorb_info('  *** Not enough memory for using mpi_alltoall, falling back to mpi_sendrecv')
      endif

    endif

    if(.not.use_alltoall_really) then
       allocate(matrix_send(maxval(send_mat_count)))
       allocate(matrix_recv(maxval(recv_mat_count)))
    endif


    if(use_alltoall_really) then
      ! Send the complete local matrix to the receivers using 1 mpi_alltoallv call
      ! This is the most effective way but needs a large amount of memory

      ! Put matrix into matrix_send ...
      send_mat_off(:) = send_mat_displ(:)
      do i_k_point = 1, n_k_points

        k_nprow = k_point_desc(i_k_point)%nprow
        k_npcol = k_point_desc(i_k_point)%npcol

        do j=1,n_basis_local
          i_row = i_basis_local(j)
          ip_r(j) = mod((Cbasis_to_basis(i_row)-1)/nb, k_nprow) ! Processor row owning i_row
        enddo

        do i=1,n_basis_local
          i_col = i_basis_local(i)
          ip_c = mod((Cbasis_to_basis(i_col)-1)/nb, k_npcol) ! Processor column owning i_col
          do j=1,n_basis_local
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              ip = k_point_desc(i_k_point)%global_id(ip_r(j),ip_c)
              send_mat_off(ip) = send_mat_off(ip)+1
              if(j<=i) then
                matrix_send(send_mat_off(ip)) = matrix(i*(i-1)/2+j)
              else
                matrix_send(send_mat_off(ip)) = matrix(j*(j-1)/2+i)
              endif
            endif
          enddo
        enddo

      enddo

      ! ... and send it away
      call mpi_alltoallv(matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
                         matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
                         mpi_comm_global, mpierr)
    endif

    ! Insert data from local matrix into scalapack matrix
    ! This is done for every remote task separatly for the case that
    ! mpi_alltoallv cannot be used

    do i_diff = 0, n_tasks-1

      i_send = mod(myid+i_diff,n_tasks)         ! Task to which we send data
      i_recv = mod(myid+n_tasks-i_diff,n_tasks) ! Task from which we get data

      if(.not.use_alltoall_really) then

        ! Gather all local rows going to proc i_send
        n_rows_local = 0
        do j=1,n_basis_local
          i_row = i_basis_local(j)
          if(mod((Cbasis_to_basis(i_row)-1)/nb, gl_nprow(i_send)) == gl_prow(i_send)) then
            n_rows_local = n_rows_local+1
            j_row_local(n_rows_local) = j
          endif
        enddo

        n_off = 0
        do i=1,n_basis_local
          i_col = i_basis_local(i)
          if(mod((Cbasis_to_basis(i_col)-1)/nb, gl_npcol(i_send)) /= gl_pcol(i_send)) cycle
          do jj=1,n_rows_local
            j = j_row_local(jj)
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              n_off = n_off+1
              if(j<=i) then
                matrix_send(n_off) = matrix(i*(i-1)/2+j)
              else
                matrix_send(n_off) = matrix(j*(j-1)/2+i)
              endif
            endif
          enddo
        enddo

        ! Gather and send data for remote task i_send,
        ! receive corresponding data from task i_recv
        call mpi_sendrecv(matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
                          matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
                          mpi_comm_global, mpi_status_ignore, mpierr)
      endif

      if(use_alltoall_really) then
        i_cnt_mat = recv_mat_displ(i_recv) ! Counter in matrix_recv
      else
        i_cnt_mat = 0 ! Counter in matrix_recv, reset for every task
      endif

      if(n_periodic==0) then

        ! We could also use the code for the periodic case here,
        ! but the code below is much more simple (and hopefully faster!)
        do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

          i_col = basis_col(j)
          lc = l_col(i_col)

          do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

            i_row = basis_row(i)
            if(i_row > i_col) exit ! done with this column
            lr = l_row(i_row)
            i_cnt_mat = i_cnt_mat+1

            if (real_eigenvectors) then
              loc_mat(lr,lc) = loc_mat(lr,lc) + matrix_recv(i_cnt_mat)
            else
              ! This code path is rare but can happen (for example, SOC
              ! for non-periodic systes)
              loc_mat_complex(lr,lc) = loc_mat_complex(lr,lc) &
                   + (1.0d0,0.0d0) * matrix_recv(i_cnt_mat)
            end if
          enddo
        enddo

      else ! periodic case

        do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

          i_col = basis_col(j)
          i_cell_c = center_to_cell(Cbasis_to_center(i_col))
          lc = l_col(Cbasis_to_basis(i_col))

          do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

            i_row = basis_row(i)
            if(Cbasis_to_basis(i_row) > Cbasis_to_basis(i_col)) cycle

            i_cell_r = center_to_cell(Cbasis_to_center(i_row))
            lr = l_row(Cbasis_to_basis(i_row))

            i_cell = position_in_hamiltonian(i_cell_c, i_cell_r) ! Attention: position_in_hamiltonian is not symmetric!

            i_cnt_mat = i_cnt_mat+1

            if(i_cell == n_cells_in_hamiltonian) cycle

            if(real_eigenvectors)then
              loc_mat(lr,lc) = loc_mat(lr,lc) + &
                   dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
            else
              loc_mat_complex(lr,lc) = loc_mat_complex(lr,lc) + &
                   k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
            end if
          enddo
        enddo

      endif

    enddo

    deallocate(matrix_send)
    deallocate(matrix_recv)

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
       write(use_unit,"(2X,A,F13.6,A)") "| Time set_full_local_matrix_scalapack_generic:",mpi_wtime()-ttt0," s"

  end subroutine set_full_local_matrix_scalapack_generic
!******

!****s* scalapack_wrapper/get_set_full_local_matrix_scalapack_cpscf
!  NAME
!    get_set_full_local_matrix_scalapack_cpscf
!  SYNOPSIS
  subroutine get_set_full_local_matrix_scalapack_cpscf( matrix, which, i_spin )
   !  PURPOSE
   !    Gets or sets a ScaLAPACK matrix in local-index mode (working with full local matrices)
   !    which = 0: Set ovlp from local overlap matrix (in matrix)
   !    which = 1: Set first_order_ham_polar_reduce_memory_scalapack from local hamiltonian (in matrix)
   !    which = 2: Set matrix from scalapack first_order_ham_polar_reduce_memory_scalapack (first_order_ham_polar_reduce_memory_scalapack must be set to something usefull before call)
   !
   !    WPH:  To understand this subroutine, it is recommend that you first read
   !    and absorb the init_comm_full_local_matrix_scalapack subroutine, which sets
   !    up the communication for this subroutine.  It is in that subroutine where
   !    most of the actual logic is performed to determine which chunks of data
   !    goes where and the necessary data structures are defined and set up. This
   !    subroutine then performs the heavy-and-stupid lifting of shuffling the data
   !    and multiplying by k-point phases when needed.
   !  USES
   use dimensions, only: n_periodic
   use localorb_io, only: localorb_info
   use pbc_lists, only: cbasis_to_basis, center_to_cell, cbasis_to_center, &
       position_in_hamiltonian, n_cells_in_hamiltonian, k_phase
   use mpi_tasks, only: n_tasks, mpi_wtime, mpi_logical, mpi_land, &
       mpi_comm_global, mpi_real8, myid, mpi_status_ignore, aims_stop
   use runtime_choices, only: use_alltoall, output_level
   implicit none
   !  ARGUMENTS
   real*8 :: matrix(n_local_matrix_size)
   integer :: which, i_spin
   !  INPUTS
   !    o matrix -- the matrix to set/get
   !    o which -- the operation to perform
   !    o i_spin -- the spin channel
   !  OUTPUT
   !    if which = 0/1 upper half of the ScaLAPACK array first_order_ham_polar_reduce_memory_scalapack/ovlp is set on exit
   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2008).
   !  SOURCE

   integer :: i, j, ip, i_k_point, k_nprow, k_npcol, ip_c, i_cell_r, i_cell_c
   integer :: i_cell, i_row, i_col, i_diff, i_send, i_recv, i_cnt_mat
   integer :: jj, n_off, n_rows_local, j_row_local(n_basis_local)
   integer :: lc, lr, istat1, istat2, mpierr
   integer :: send_mat_off(0:n_tasks-1), ip_r(n_basis_local)
   real*8, allocatable :: matrix_send(:), matrix_recv(:)
   logical alloc_success, use_alltoall_really
   real*8 ttt0

   character(*), parameter :: func = "get_set_full_local_matrix_scalapack_cpscf"

   ttt0 = mpi_wtime()

   ! Reset the matrix to be set, check parameter which

   if(which==0) then
       if(real_eigenvectors)then
           ovlp(:,:) = 0
       else
           ovlp_complex(:,:) = 0
       end if
   else if(which==1) then
       if(real_eigenvectors)then
           first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin) = 0
       else
           call aims_stop("cpscf poalr_reduce_memory error: ham_complex")
           first_order_ham_complex_scalapack(:,:,i_spin) = 0
       end if
   else if(which==2) then
       matrix(:) = 0
   else
       call aims_stop("Illegal parameter for 'which'", func)
   endif

   use_alltoall_really = use_alltoall

   ! When use_alltoall is set we check if we can allocate the memory needed
   ! and if this doesn't succeed we fall back to sendrecv.

   ! The reason why we don't do that always this way but use a config switch is the fact that
   ! the system might be configured with swap space. In this case the allocate may succeed
   ! but the machine will start swapping - which is absolutely nothing we want!

   if(use_alltoall) then

       ! Try to allocate large matrices for using mpi_alltoallv ...
       allocate(matrix_send(send_mat_count_tot), stat=istat1)
       allocate(matrix_recv(recv_mat_count_tot), stat=istat2)

       ! ... check if allocation succeeded on ALL procs, otherways use mpi_sendrecv

       alloc_success = (istat1==0 .and. istat2==0)
       call mpi_allreduce(alloc_success, use_alltoall_really, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_global, mpierr)

       if(.not.use_alltoall_really) then
           ! fall back to sendrecv
           if(allocated(matrix_send)) deallocate(matrix_send)
           if(allocated(matrix_recv)) deallocate(matrix_recv)
           call localorb_info('  *** Not enough memory for using mpi_alltoall, falling back to mpi_sendrecv')
       endif

   endif

   if(.not.use_alltoall_really) then
       allocate(matrix_send(maxval(send_mat_count)))
       allocate(matrix_recv(maxval(recv_mat_count)))
   endif


   if(which<=1 .and. use_alltoall_really) then
       ! Send the complete local matrix to the receivers using 1 mpi_alltoallv call
       ! This is the most effective way but needs a large amount of memory

       ! Put matrix into matrix_send ...
       send_mat_off(:) = send_mat_displ(:)
       do i_k_point = 1, n_k_points

       k_nprow = k_point_desc(i_k_point)%nprow
       k_npcol = k_point_desc(i_k_point)%npcol

       do j=1,n_basis_local
       i_row = i_basis_local(j)
       ip_r(j) = mod((Cbasis_to_basis(i_row)-1)/nb, k_nprow) ! Processor row owning i_row
       enddo

       do i=1,n_basis_local
       i_col = i_basis_local(i)
       ip_c = mod((Cbasis_to_basis(i_col)-1)/nb, k_npcol) ! Processor column owning i_col
       do j=1,n_basis_local
       i_row = i_basis_local(j)
       if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
           ip = k_point_desc(i_k_point)%global_id(ip_r(j),ip_c)
           send_mat_off(ip) = send_mat_off(ip)+1
           if(j<=i) then
               matrix_send(send_mat_off(ip)) = matrix(i*(i-1)/2+j)
           else
               matrix_send(send_mat_off(ip)) = matrix(j*(j-1)/2+i)
           endif
       endif
       enddo
       enddo

       enddo

       ! ... and send it away
       call mpi_alltoallv(matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
           matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
           mpi_comm_global, mpierr)
   endif

   if(which>1 .and. use_alltoall_really) matrix_recv(:) = 0

   ! Insert data from local matrix into scalapack matrix (which<=1)
   ! or gather data from scalapack matrix for putting into local matrix (which>1)
   ! This is done for every remote task separatly for the case that
   ! mpi_alltoallv cannot be used

   do i_diff = 0, n_tasks-1

   i_send = mod(myid+i_diff,n_tasks)         ! Task to which we send data
   i_recv = mod(myid+n_tasks-i_diff,n_tasks) ! Task from which we get data

   if(.not.use_alltoall_really) then

       ! Gather all local rows going to proc i_send
       ! (for which > 1 this is needed below)
       n_rows_local = 0
       do j=1,n_basis_local
       i_row = i_basis_local(j)
       if(mod((Cbasis_to_basis(i_row)-1)/nb, gl_nprow(i_send)) == gl_prow(i_send)) then
           n_rows_local = n_rows_local+1
           j_row_local(n_rows_local) = j
       endif
       enddo

       if(which<=1) then

           n_off = 0
           do i=1,n_basis_local
           i_col = i_basis_local(i)
           if(mod((Cbasis_to_basis(i_col)-1)/nb, gl_npcol(i_send)) /= gl_pcol(i_send)) cycle
           do jj=1,n_rows_local
           j = j_row_local(jj)
           i_row = i_basis_local(j)
           if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
               n_off = n_off+1
               if(j<=i) then
                   matrix_send(n_off) = matrix(i*(i-1)/2+j)
               else
                   matrix_send(n_off) = matrix(j*(j-1)/2+i)
               endif
           endif
           enddo
           enddo

           ! Gather and send data for remote task i_send,
           ! receive corresponding data from task i_recv
           call mpi_sendrecv(matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
               matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
               mpi_comm_global, mpi_status_ignore, mpierr)

       else

           matrix_recv(1:recv_mat_count(i_recv)) = 0

       endif

   endif

   if(use_alltoall_really) then
       i_cnt_mat = recv_mat_displ(i_recv) ! Counter in matrix_recv
   else
       i_cnt_mat = 0 ! Counter in matrix_recv, reset for every task
   endif

   if(n_periodic==0) then

       ! We could also use the code for the periodic case here,
       ! but the code below is much more simple (and hopefully faster!)
       do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

       i_col = basis_col(j)
       lc = l_col(i_col)

       do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

       i_row = basis_row(i)
       if(i_row > i_col) exit ! done with this column
       lr = l_row(i_row)
       i_cnt_mat = i_cnt_mat+1

       if(which==0) then
           ovlp(lr,lc) = ovlp(lr,lc) + matrix_recv(i_cnt_mat)
       else if(which==1) then
           first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin) = first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin) + matrix_recv(i_cnt_mat)
       else
           matrix_recv(i_cnt_mat) = first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin)
       endif
       enddo
       enddo

   else ! periodic case

       do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

       i_col = basis_col(j)
       i_cell_c = center_to_cell(Cbasis_to_center(i_col))
       lc = l_col(Cbasis_to_basis(i_col))

       do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

       i_row = basis_row(i)
       if(Cbasis_to_basis(i_row) > Cbasis_to_basis(i_col)) cycle

       i_cell_r = center_to_cell(Cbasis_to_center(i_row))
       lr = l_row(Cbasis_to_basis(i_row))

       i_cell = position_in_hamiltonian(i_cell_c, i_cell_r) ! Attention: position_in_hamiltonian is not symmetric!

       i_cnt_mat = i_cnt_mat+1

       if(i_cell == n_cells_in_hamiltonian) cycle

       if(which==0) then
           if(real_eigenvectors)then
               ovlp(lr,lc) = ovlp(lr,lc) + dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
           else
               ovlp_complex(lr,lc) = ovlp_complex(lr,lc) + k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
           end if
       else if(which==1) then
           if(real_eigenvectors)then
               first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin) = first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin) + dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
           else
               first_order_ham_complex_scalapack(lr,lc,i_spin) = first_order_ham_complex_scalapack(lr,lc,i_spin) + k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
           end if
       else
           if(real_eigenvectors)then
               matrix_recv(i_cnt_mat) = first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin)*dble(k_phase(i_cell,my_k_point))
           else
               matrix_recv(i_cnt_mat) = dble(first_order_ham_complex_scalapack(lr,lc,i_spin)*dconjg(k_phase(i_cell,my_k_point)))
           end if
       endif
       enddo
       enddo

   endif

   if(which>1 .and. .not.use_alltoall_really) then

       ! Send matrix_recv immediatly back to owner of local matrix
       call mpi_sendrecv(matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
           matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
           mpi_comm_global, mpi_status_ignore, mpierr)

       n_off = 0
       do i=1,n_basis_local
       i_col = i_basis_local(i)
       if(mod((Cbasis_to_basis(i_col)-1)/nb, gl_npcol(i_send)) /= gl_pcol(i_send)) cycle
       do jj=1,n_rows_local
       j = j_row_local(jj)
       i_row = i_basis_local(j)
       if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
           n_off = n_off+1
           ! some of the elements we got are duplicates, we must be careful
           ! not to add them twice!!!!
           if(Cbasis_to_basis(i_row)==Cbasis_to_basis(i_col) .and. j>i) cycle
           if(j<=i) then
               matrix(i*(i-1)/2+j) = matrix(i*(i-1)/2+j) + matrix_send(n_off)
           else
               matrix(j*(j-1)/2+i) = matrix(j*(j-1)/2+i) + matrix_send(n_off)
           endif
       endif
       enddo
       enddo

   endif

   enddo

   if(which>1 .and. use_alltoall_really) then
       ! Send the matrix gathered in matrix_recv to owners of local matrix with mpi_alltoallv
       ! This is the most effective way but needs a large amount of memory
       call mpi_alltoallv(matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
           matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
           mpi_comm_global, mpierr)

       ! Insert received matrix
       send_mat_off(:) = send_mat_displ(:)
       do i_k_point = 1, n_k_points

       k_nprow = k_point_desc(i_k_point)%nprow
       k_npcol = k_point_desc(i_k_point)%npcol

       do j=1,n_basis_local
       i_row = i_basis_local(j)
       ip_r(j) = mod((Cbasis_to_basis(i_row)-1)/nb, k_nprow) ! Processor row owning i_row
       enddo

       do i=1,n_basis_local
       i_col = i_basis_local(i)
       ip_c = mod((Cbasis_to_basis(i_col)-1)/nb, k_npcol) ! Processor column owning i_col
       do j=1,n_basis_local
       i_row = i_basis_local(j)
       if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
           ip = k_point_desc(i_k_point)%global_id(ip_r(j),ip_c)
           send_mat_off(ip) = send_mat_off(ip)+1
           ! some of the elements we got are duplicates, we must be careful
           ! not to add them twice!!!!
           if(Cbasis_to_basis(i_row)==Cbasis_to_basis(i_col) .and. j>i) cycle
           if(j<=i) then
               matrix(i*(i-1)/2+j) = matrix(i*(i-1)/2+j) + matrix_send(send_mat_off(ip))
           else
               matrix(j*(j-1)/2+i) = matrix(j*(j-1)/2+i) + matrix_send(send_mat_off(ip))
           endif
       endif
       enddo
       enddo

       enddo
   endif

   deallocate(matrix_send)
   deallocate(matrix_recv)

   if(which==0) then
       factor_overlap = .TRUE.
       use_ovlp_trafo = .FALSE.
       full_ovlp_ready = .false.
       n_nonsing_ovlp = n_basis
   endif

   if((myid==0).and.(output_level .ne. 'MD_light' )) &
       write(use_unit,"(2X,A,F13.6,A)") "| Time get_set_full_local_matrix_scalapack_cpscf:",mpi_wtime()-ttt0," s"

end subroutine get_set_full_local_matrix_scalapack_cpscf


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/init_comm_full_local_matrix_scalapack
!  NAME
!    init_comm_full_local_matrix_scalapack
!  SYNOPSIS
  subroutine init_comm_full_local_matrix_scalapack(n_basis_local_full, i_basis_local_full)
!  PURPOSE
!    Initializes the communication for get_set_full_local_matrix_scalapack
!  USES
    use dimensions, only: n_centers_basis_T
    use mpi_tasks, only: n_tasks, myid, mpi_comm_global, mpi_integer, aims_stop
    use pbc_lists, only: Cbasis_to_basis
    use synchronize_mpi_basic, only: sync_integer_vector
    implicit none
!  ARGUMENTS

    integer, intent(in) :: n_basis_local_full
    integer, intent(in) :: i_basis_local_full(n_basis_local_full)

!  INPUTS
!    o n_basis_local_full -- number of local basis functions for current task
!    o i_basis_local_full -- array with local basis functions for current task
!  OUTPUT
!    all index arrays for local_index communication are set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    ! Variables used for setting up receiving data
    integer :: irow, icol, lr, lc

    ! Variables used for setting up sending data
    integer :: k_nprow, k_npcol, i_k_point, ip_r, ip_c, i_col, i_row

    ! Generic multipurpose counters and dummy parameters
    integer :: i, j, ip, ncnt, mpierr

    ! The dimension of the local matrix for a given MPI task
    integer :: n_basis_all(0:n_tasks)

    ! The total number of copies of basis elements that the current MPI task
    ! will receive from other MPI tasks (and itself)
    integer :: ncnt_row, ncnt_col

    ! The number of MPI tasks whose local matrix contains a given basis element
    integer :: i_basis(n_centers_basis_T)

    character(*), parameter :: func = "init_comm_full_local_matrix_scalapack"

    ! Set up the dimensions and basis functions contained in the local matrix
    ! for *this* MPI tasks (a.k.a. myid)
    n_basis_local = n_basis_local_full
    if(allocated(i_basis_local)) deallocate(i_basis_local)
    allocate(i_basis_local(n_basis_local))
    i_basis_local(:) = i_basis_local_full(:)
    n_local_matrix_size = n_basis_local*(n_basis_local+1)/2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               Receiving Data                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! WPH: We here set up the communication where myid obtains local matrix
    ! elements from other MPI tasks so that it may construct its BLACS-formatted
    ! matrix.
    !
    ! This requires some communication from other MPI tasks to construct various
    ! indexing arrays as, by design, myid has no method to determine the details
    ! of the local matrices for other MPI tasks without asking them.

    ! Get total number of tasks having each basis function

    i_basis(:) = 0
    i_basis(i_basis_local(1:n_basis_local)) = 1
    call sync_integer_vector(i_basis,n_centers_basis_T)

    ! Count how much basis functions will go to my task

    ncnt_row = 0
    ncnt_col = 0
    do i=1,n_centers_basis_T
      lr = l_row(Cbasis_to_basis(i))
      lc = l_col(Cbasis_to_basis(i))
      if(lr>0) ncnt_row = ncnt_row+i_basis(i)
      if(lc>0) ncnt_col = ncnt_col+i_basis(i)
    enddo

    ! Allocate arrays for remote basis functions

    if(allocated(basis_row)) deallocate(basis_row)
    if(allocated(basis_col)) deallocate(basis_col)
    allocate(basis_row(ncnt_row))
    allocate(basis_col(ncnt_col))

    if(.not.allocated(basis_col_limit)) allocate(basis_col_limit(0:n_tasks))
    if(.not.allocated(basis_row_limit)) allocate(basis_row_limit(0:n_tasks))
    basis_col_limit(0) = 0
    basis_row_limit(0) = 0

    ! Get all remote basis functions

    ! First, tell every MPI task the dimensions of the local matrices on every
    ! other MPI task
    call mpi_allgather(n_basis_local,1,MPI_INTEGER,n_basis_all,1,MPI_INTEGER,mpi_comm_global,mpierr)

    ! Next, we loop through each MPI task, broadcasting the basis functions of
    ! its local matrices to every other MPI task
    do ip=0,n_tasks-1
      if(ip==myid) i_basis(1:n_basis_local) = i_basis_local(:)
      call mpi_bcast(i_basis,n_basis_all(ip),mpi_integer,ip,mpi_comm_global,mpierr)
      irow = 0
      icol = 0
      ! We check to see if *this* MPI task (a.k.a myid) contains the same basis
      ! function in its BLACS matrix as the broadcasting MPI task (a.k.a.
      ! ip) does in its local matrix.  (We'll call this the "current MPI task
      ! pair" from now on.)
      do i=1,n_basis_all(ip)
        lr = l_row(Cbasis_to_basis(i_basis(i)))
        lc = l_col(Cbasis_to_basis(i_basis(i)))
        ! When they match, update the number of matches (irow/icol) for the
        ! current MPI task pair and push the index of the basis function to the
        ! top of basis_row/basis_col
        if(lr>0) then
          irow = irow+1
          basis_row(basis_row_limit(ip)+irow) = i_basis(i)
        endif
        if(lc>0) then
          icol = icol+1
          basis_col(basis_col_limit(ip)+icol) = i_basis(i)
        endif
      enddo
      ! Update the offset for the next broadcasting MPI task based on the number
      ! of matches for the current MPI task pair
      basis_row_limit(ip+1) = basis_row_limit(ip) + irow
      basis_col_limit(ip+1) = basis_col_limit(ip) + icol
    enddo

    ! Plausibility checks
    if(basis_row_limit(n_tasks) /= ncnt_row) call aims_stop('INTERNAL ERROR: basis_row_limit', func)
    if(basis_col_limit(n_tasks) /= ncnt_col) call aims_stop('INTERNAL ERROR: basis_col_limit', func)

    ! set up count/displ arrays for matrix send/recv

    if(.not.allocated(send_mat_count)) allocate(send_mat_count(0:n_tasks-1))
    if(.not.allocated(send_mat_displ)) allocate(send_mat_displ(0:n_tasks-1))
    if(.not.allocated(recv_mat_count)) allocate(recv_mat_count(0:n_tasks-1))
    if(.not.allocated(recv_mat_displ)) allocate(recv_mat_displ(0:n_tasks-1))

    !---- what we recv:

    ! Knowing now the indices of the basis elements that will be communicated
    ! between the MPI task pairs, we'll now loop through broadcasting MPI tasks
    ! again to determine the number of matrix elements that will be communicated
    recv_mat_displ(0) = 0
    do ip=0,n_tasks-1
      ncnt = 0
      ! Loop over all pairs of basis elements involved in communication, only
      ! keeping those that lie in the upper triangular portion of the matrix
      do i=basis_col_limit(ip)+1,basis_col_limit(ip+1)
        i_col = basis_col(i)
        do j=basis_row_limit(ip)+1,basis_row_limit(ip+1)
           i_row = basis_row(j)
           if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) ncnt = ncnt+1
           ! wyj: TODO lower triangular
           !if(Cbasis_to_basis(i_row) >= Cbasis_to_basis(i_col)) ncnt = ncnt+1
        enddo
      enddo
      ! Update the number of matrix elements that will be received by myid from
      ! the current MPI task pair as well as the offset for the next MPI task
      ! pair
      recv_mat_count(ip) = ncnt
      if(ip>0) recv_mat_displ(ip) = recv_mat_displ(ip-1) + recv_mat_count(ip-1)
    enddo
    recv_mat_count_tot = recv_mat_displ(n_tasks-1)+recv_mat_count(n_tasks-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                Sending Data                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! WPH: We here set up the communication where myid sends its local matrix
    ! elements to other MPI tasks so that they may construct their
    ! BLACS-formatted matrix.
    !
    ! This section requires no communication with other MPI tasks due to the
    ! deterministic setup of the BLACS matrix format.  myid already knows what
    ! local matrix elements it has, and it can easily surmise the rows and
    ! columns contained in the BLACS matrices of other MPI tasks based on the
    ! k-point descriptors, giving myid all the information needed to know which
    ! MPI tasks it will be communicating with in advance.
    !
    ! This section is independent of the "Receiving Data" section and could
    ! easily be moved to its own subroutine.

    !---- what we send

    send_mat_count(:) = 0
    do i_k_point = 1, n_k_points

       k_nprow = k_point_desc(i_k_point)%nprow
       k_npcol = k_point_desc(i_k_point)%npcol

       ! For each k-point, we loop over every matrix element contained in the
       ! local matrix on myid.  Because we are using the dense BLACS format,
       ! there must exist one and only one MPI task which contains that matrix
       ! element at the given k-point in its BLACS matrix.
       do i=1,n_basis_local
         i_col = i_basis_local(i)
         ip_c = mod((Cbasis_to_basis(i_col)-1)/nb, k_npcol) ! Processor column owning i_col
         do j=1,n_basis_local
           i_row = i_basis_local(j)
           ip_r = mod((Cbasis_to_basis(i_row)-1)/nb, k_nprow) ! Processor row owning i_row
           ! As always, we only consider the upper triangular matrix elements
           if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
           ! wyj: TODO lower triangular
           !if(Cbasis_to_basis(i_row) >= Cbasis_to_basis(i_col)) then
             ip = k_point_desc(i_k_point)%global_id(ip_r,ip_c)
             send_mat_count(ip) = send_mat_count(ip)+1
           endif
         enddo
       enddo

    enddo

    send_mat_displ(0) = 0
    do ip=1,n_tasks-1
       send_mat_displ(ip) = send_mat_displ(ip-1) + send_mat_count(ip-1)
    enddo
    send_mat_count_tot = send_mat_displ(n_tasks-1)+send_mat_count(n_tasks-1)

  end subroutine init_comm_full_local_matrix_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_set_sparse_local_matrix_scalapack
!  NAME
!    get_set_sparse_local_matrix_scalapack
!  SYNOPSIS
  subroutine get_set_sparse_local_matrix_scalapack( matrix, which, i_spin )
!  PURPOSE
!    Gets or sets a ScaLAPACK matrix in local-index mode (working with sparse local matrices)
!    which = 0: Set ovlp from local overlap matrix (in matrix)
!    which = 1: Set ham from local hamiltonian (in matrix)
!    which = 2: Set matrix from scalapack ham (ham must be set to something usefull before call)
!  USES
    use runtime_choices
    use basis
    use localorb_io, only: localorb_info
    use mpi_tasks
    use pbc_lists
    use geometry
    implicit none
!  ARGUMENTS
    real*8:: matrix(n_hamiltonian_matrix_size)
    integer :: which, i_spin
!  INPUTS
!    o matrix -- the matrix to set/get
!    o which -- the operation to perform
!    o i_spin -- the spin channel
!  OUTPUT
!    if which = 0/1 upper half of the ScaLAPACK array ham/ovlp is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cnt, i_cell, i_col, i_ccc, i_diff, i_send, i_recv, i_cnt_idx, i_cnt_mat
    integer :: lc, lr, idx, n_rows, istat1, istat2
    integer mpierr
    real*8, allocatable :: matrix_send(:), matrix_recv(:)
    logical alloc_success, use_alltoall_really
real*8 ttt0
ttt0 = mpi_wtime()


    ! Reset the matrix to be set, check parameter which

    if(which==0) then
       if(real_eigenvectors)then
          ovlp(:,:) = 0
       else
          ovlp_complex(:,:) = 0
       end if
    else if(which==1) then
       if(real_eigenvectors)then
          ham(:,:,i_spin) = 0
       else
          ham_complex(:,:,i_spin) = 0
       end if
    else if(which==2) then
       matrix(1:n_hamiltonian_matrix_size) = 0
    else
       write(use_unit,*) 'Illegal parameter ''which'' in get_set_sparse_local_matrix_scalapack: ',which
       call aims_stop
    endif

    use_alltoall_really = use_alltoall

    ! When use_alltoall is set we check if we can allocate the memory needed
    ! and if this doesn't succeed we fall back to sendrecv.

    ! The reason why we don't do that always this way but use a config switch is the fact that
    ! the system might be configured with swap space. In this case the allocate may succeed
    ! but the machine will start swapping - which is absolutely nothing we want!

    if(use_alltoall) then

      ! Try to allocate large matrices for using mpi_alltoallv ...
      allocate(matrix_send(send_idx_count_tot), stat=istat1)
      allocate(matrix_recv(recv_row_count_tot), stat=istat2)

      ! ... check if allocation succeeded on ALL procs, otherways use mpi_sendrecv

      alloc_success = (istat1==0 .and. istat2==0)
      call mpi_allreduce(alloc_success, use_alltoall_really, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_global, mpierr)

      if(.not.use_alltoall_really) then
        ! fall back to sendrecv
        if(allocated(matrix_send)) deallocate(matrix_send)
        if(allocated(matrix_recv)) deallocate(matrix_recv)
        call localorb_info('  *** Not enough memory for using mpi_alltoall, falling back to mpi_sendrecv')
      endif

    endif

    if(.not.use_alltoall_really) then
       allocate(matrix_send(maxval(send_idx_count)))
       allocate(matrix_recv(maxval(recv_row_count)))
    endif

    if(which<=1 .and. use_alltoall_really) then
       ! Send the complete local matrix to the receivers using 1 mpi_alltoallv call
       ! This is the most effective way but needs a large amount of memory
       do i_cnt = 1, send_idx_count_tot
          matrix_send(i_cnt) = matrix(send_idx(i_cnt))
       enddo
       call mpi_alltoallv(matrix_send, send_idx_count, send_idx_displ, MPI_REAL8, &
                          matrix_recv, recv_row_count, recv_row_displ, MPI_REAL8, &
                          mpi_comm_global, mpierr)
    endif

    ! Insert data from local matrix into scalapack matrix (which<=1)
    ! or gather data from scalapack matrix for putting into local matrix (which>1)
    ! This is done for every remote task separatly for the case that
    ! mpi_alltoallv cannot be used

    do i_diff = 0, n_tasks-1

       i_send = mod(myid+i_diff,n_tasks)         ! Task to which we send data
       i_recv = mod(myid+n_tasks-i_diff,n_tasks) ! Task from which we get data

       if(which<=1 .and. .not.use_alltoall_really) then
          ! Gather and send data for remote task i_send,
          ! receive corresponding data from task i_recv
          do i_cnt = 1, send_idx_count(i_send)
             matrix_send(i_cnt) = matrix(send_idx(i_cnt+send_idx_displ(i_send)))
          enddo
          call mpi_sendrecv(matrix_send, send_idx_count(i_send), MPI_REAL8, i_send, 111, &
                            matrix_recv, recv_row_count(i_recv), MPI_REAL8, i_recv, 111, &
                            mpi_comm_global, mpi_status_ignore, mpierr)
       endif

       i_cnt_idx = recv_row_displ(i_recv) ! Counter in recv_row

       if(use_alltoall_really) then
          i_cnt_mat = recv_row_displ(i_recv) ! Counter in matrix_recv
       else
          i_cnt_mat = 0 ! Counter in matrix_recv, reset for every task
       endif

       do i_ccc = recv_ccc_displ(i_recv)+1, recv_ccc_displ(i_recv)+recv_ccc_count(i_recv)

          i_cell = recv_ccc(1,i_ccc)
          i_col  = recv_ccc(2,i_ccc)
          n_rows = recv_ccc(3,i_ccc)

          lc = l_col(i_col) ! is always > 0 (ie. local)

          do idx = 1, n_rows

             i_cnt_idx = i_cnt_idx+1
             lr = l_row(recv_row(i_cnt_idx)) ! is always > 0 (ie. local)

             i_cnt_mat = i_cnt_mat+1
             if(which==0) then
                if(real_eigenvectors)then
                   ovlp(lr,lc) = ovlp(lr,lc) + dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
                else
                   ovlp_complex(lr,lc) = ovlp_complex(lr,lc) + k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
                end if
             else if(which==1) then
                if(real_eigenvectors)then
                   ham(lr,lc,i_spin) = ham(lr,lc,i_spin) + dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
                else
                   ham_complex(lr,lc,i_spin) = ham_complex(lr,lc,i_spin) + &
                     k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
                end if
             else
                if(real_eigenvectors)then
                   matrix_recv(i_cnt_mat) = ham(lr,lc,i_spin)*dble(k_phase(i_cell,my_k_point))
                else
                   matrix_recv(i_cnt_mat) = dble(ham_complex(lr,lc,i_spin)*dconjg(k_phase(i_cell,my_k_point)))
                end if
             endif
          enddo
       enddo

       if(which>1 .and. .not.use_alltoall_really) then
          ! Send matrix_recv immediatly back to owner of local matrix
          call mpi_sendrecv(matrix_recv, recv_row_count(i_recv), MPI_REAL8, i_recv, 111, &
                            matrix_send, send_idx_count(i_send), MPI_REAL8, i_send, 111, &
                            mpi_comm_global, mpi_status_ignore, mpierr)
          do i_cnt = 1, send_idx_count(i_send)
             matrix(send_idx(i_cnt+send_idx_displ(i_send))) = &
                matrix(send_idx(i_cnt+send_idx_displ(i_send))) + matrix_send(i_cnt)
          enddo
       endif

    enddo

    if(which>1 .and. use_alltoall_really) then
       ! Send the matrix gathered in matrix_recv to owners of local matrix with mpi_alltoallv
       ! This is the most effective way but needs a large amount of memory
       call mpi_alltoallv(matrix_recv, recv_row_count, recv_row_displ, MPI_REAL8, &
                          matrix_send, send_idx_count, send_idx_displ, MPI_REAL8, &
                          mpi_comm_global, mpierr)
       ! Add matrix_send (which we actually received in this case) into local matrix
       do i_cnt = 1, send_idx_count_tot
          matrix(send_idx(i_cnt)) = matrix(send_idx(i_cnt)) + matrix_send(i_cnt)
       enddo
    endif

    deallocate(matrix_send)
    deallocate(matrix_recv)

    if(which==0) then
       factor_overlap = .TRUE.
       use_ovlp_trafo = .FALSE.
       full_ovlp_ready = .false.
       n_nonsing_ovlp = n_basis
    endif

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
      write(use_unit,"(2X,A,F13.6,A)") "| Time get_set_sparse_local_matrix_scalapack:",mpi_wtime()-ttt0," s"

  end subroutine get_set_sparse_local_matrix_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_sparse_local_matrix_scalapack_generic
!  NAME
!    set_sparse_local_matrix_scalapack_generic
!  SYNOPSIS
  subroutine set_sparse_local_matrix_scalapack_generic( matrix, loc_mat, loc_mat_complex )
!  PURPOSE
!    Sets a ScaLAPACK matrix in local-index mode (working with sparse local matrices)
!    This is the relevant code path when using local index but no load balancing
!    This subroutine is a fork of get_set_sparse_local_matrix_scalapack (which == 1);
!    the difference is that the desired ScaLAPACK matrix is passed as an argument
!    rather than pulled from a module variable, allowing for reuse by generic
!    matrices sharing the same packing (and ScaLAPACK descriptor) as the Hamiltonian
!    or overlap matrices
!    We deliberately do not set the incoming local matrices to zero, to allow
!    for accumulation
!  USES
    use localorb_io, only: localorb_info
    use mpi_tasks
    use pbc_lists
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8,      intent(in)    :: matrix(n_hamiltonian_matrix_size)
    real*8,      intent(inout) :: loc_mat(mxld,mxcol)
    complex*16,  intent(inout) :: loc_mat_complex(mxld,mxcol)
!  INPUTS
!    o matrix -- the matrix, in local indexing
!  OUTPUT
!    The upper half of the appropriate ScaLAPACK array are set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    September 2017 - Forked off of get_set_sparse_local_matrix_scalapack by
!                     William Huhn (Duke University)
!  SOURCE

    integer :: i_cnt, i_cell, i_col, i_ccc, i_diff, i_send, i_recv, i_cnt_idx, i_cnt_mat
    integer :: lc, lr, idx, n_rows, istat1, istat2
    integer mpierr
    real*8, allocatable :: matrix_send(:), matrix_recv(:)
    logical alloc_success, use_alltoall_really
real*8 ttt0
ttt0 = mpi_wtime()

    ! Reset the matrix to be set, check parameter which

    use_alltoall_really = use_alltoall

    ! When use_alltoall is set we check if we can allocate the memory needed
    ! and if this doesn't succeed we fall back to sendrecv.

    ! The reason why we don't do that always this way but use a config switch is the fact that
    ! the system might be configured with swap space. In this case the allocate may succeed
    ! but the machine will start swapping - which is absolutely nothing we want!

    if(use_alltoall) then

      ! Try to allocate large matrices for using mpi_alltoallv ...
      allocate(matrix_send(send_idx_count_tot), stat=istat1)
      allocate(matrix_recv(recv_row_count_tot), stat=istat2)

      ! ... check if allocation succeeded on ALL procs, otherways use mpi_sendrecv

      alloc_success = (istat1==0 .and. istat2==0)
      call mpi_allreduce(alloc_success, use_alltoall_really, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_global, mpierr)

      if(.not.use_alltoall_really) then
        ! fall back to sendrecv
        if(allocated(matrix_send)) deallocate(matrix_send)
        if(allocated(matrix_recv)) deallocate(matrix_recv)
        call localorb_info('  *** Not enough memory for using mpi_alltoall, falling back to mpi_sendrecv')
      endif

    endif

    if(.not.use_alltoall_really) then
       allocate(matrix_send(maxval(send_idx_count)))
       allocate(matrix_recv(maxval(recv_row_count)))
    endif

    if(use_alltoall_really) then
       ! Send the complete local matrix to the receivers using 1 mpi_alltoallv call
       ! This is the most effective way but needs a large amount of memory
       do i_cnt = 1, send_idx_count_tot
          matrix_send(i_cnt) = matrix(send_idx(i_cnt))
       enddo
       call mpi_alltoallv(matrix_send, send_idx_count, send_idx_displ, MPI_REAL8, &
                          matrix_recv, recv_row_count, recv_row_displ, MPI_REAL8, &
                          mpi_comm_global, mpierr)
    endif

    ! Insert data from local matrix into scalapack matrix
    ! This is done for every remote task separatly for the case that
    ! mpi_alltoallv cannot be used

    do i_diff = 0, n_tasks-1

       i_send = mod(myid+i_diff,n_tasks)         ! Task to which we send data
       i_recv = mod(myid+n_tasks-i_diff,n_tasks) ! Task from which we get data

       if(.not.use_alltoall_really) then
          ! Gather and send data for remote task i_send,
          ! receive corresponding data from task i_recv
          do i_cnt = 1, send_idx_count(i_send)
             matrix_send(i_cnt) = matrix(send_idx(i_cnt+send_idx_displ(i_send)))
          enddo
          call mpi_sendrecv(matrix_send, send_idx_count(i_send), MPI_REAL8, i_send, 111, &
                            matrix_recv, recv_row_count(i_recv), MPI_REAL8, i_recv, 111, &
                            mpi_comm_global, mpi_status_ignore, mpierr)
       endif

       i_cnt_idx = recv_row_displ(i_recv) ! Counter in recv_row

       if(use_alltoall_really) then
          i_cnt_mat = recv_row_displ(i_recv) ! Counter in matrix_recv
       else
          i_cnt_mat = 0 ! Counter in matrix_recv, reset for every task
       endif

       do i_ccc = recv_ccc_displ(i_recv)+1, recv_ccc_displ(i_recv)+recv_ccc_count(i_recv)

          i_cell = recv_ccc(1,i_ccc)
          i_col  = recv_ccc(2,i_ccc)
          n_rows = recv_ccc(3,i_ccc)

          lc = l_col(i_col) ! is always > 0 (ie. local)

          do idx = 1, n_rows

             i_cnt_idx = i_cnt_idx+1
             lr = l_row(recv_row(i_cnt_idx)) ! is always > 0 (ie. local)

             i_cnt_mat = i_cnt_mat+1
             if(real_eigenvectors)then
               loc_mat(lr,lc) = loc_mat(lr,lc) + &
                    dble(k_phase(i_cell,my_k_point)) * matrix_recv(i_cnt_mat)
             else
               loc_mat_complex(lr,lc) = loc_mat_complex(lr,lc) + &
                    k_phase(i_cell,my_k_point) * matrix_recv(i_cnt_mat)
             end if
          enddo
       enddo
    enddo

    deallocate(matrix_send)
    deallocate(matrix_recv)

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
      write(use_unit,"(2X,A,F13.6,A)") "| Time set_sparse_local_matrix_scalapack_generic:",mpi_wtime()-ttt0," s"

  end subroutine set_sparse_local_matrix_scalapack_generic
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_local_index_comm_arrays
!  NAME
!    set_comm_desc
!  SYNOPSIS
  subroutine set_local_index_comm_arrays(get_size)
!  PURPOSE
!    Computes the size of / sets the index arrays needed for local index communication
!  USES
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
   logical get_size
!  INPUTS
!    o get_size -- specifies if sizes should be computed or the arrays should be set
!  OUTPUT
!    if get_size is TRUE:
!      send_idx_count(:) and send_ccc_count(:) will be set,
!      send_idx(:) and send_ccc(:,:) will not be touched
!    if get_size is FALSE:
!      send_idx(:) and send_ccc(:,:) will be set,
!      send_idx_count(:) and send_ccc_count(:) must be correctly set on entry
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_k_point, i_cell, i_col, idx, k_nprow, k_npcol, ip_r, ip_c, ip
    integer, allocatable :: i_count(:), i_idx_cnt(:), i_ccc_cnt(:)


    if(get_size) then
       send_idx_count(:) = 0
       send_ccc_count(:) = 0
    else
       ! Allocate i_idx_cnt/i_ccc_cnt and let them point to the initial entry
       ! in send_idx/send_ccc for every task
       allocate(i_idx_cnt(0:n_tasks-1))
       allocate(i_ccc_cnt(0:n_tasks-1))
       i_idx_cnt(0) = 0
       i_ccc_cnt(0) = 0
       do ip=1,n_tasks-1
          i_idx_cnt(ip) = i_idx_cnt(ip-1) + send_idx_count(ip-1)
          i_ccc_cnt(ip) = i_ccc_cnt(ip-1) + send_ccc_count(ip-1)
       enddo
    endif

    do i_k_point = 1, n_k_points

       ! Number of rows/cols for current k point
       k_nprow = k_point_desc(i_k_point)%nprow
       k_npcol = k_point_desc(i_k_point)%npcol

       allocate(i_count(0:k_nprow-1))

       do i_cell = 1, n_cells_in_hamiltonian-1
          do i_col = 1, n_basis

             ip_c = mod((i_col-1)/nb, k_npcol) ! Processor column owning i_col

             if(index_hamiltonian(1,i_cell,i_col) > 0) then
                ! Count how many elements of the hamiltonian of current cell/col go to which processor row
                i_count(:) = 0
                do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)
                   ip_r = mod( (column_index_hamiltonian(idx)-1)/mb, k_nprow ) ! Processor row for this entry
                   i_count(ip_r) = i_count(ip_r) + 1
                   if(.not.get_size) then
                      ! Set corresponding entry in send_idx
                      ip = k_point_desc(i_k_point)%global_id(ip_r,ip_c)
                      i_idx_cnt(ip) = i_idx_cnt(ip)+1
                      send_idx(i_idx_cnt(ip)) = idx
                   endif
                end do
                do ip_r = 0, k_nprow-1
                   if(i_count(ip_r) > 0) then
                      ip = k_point_desc(i_k_point)%global_id(ip_r,ip_c)
                      if(get_size) then
                         ! Just count entries in send_idx / send_ccc
                         send_idx_count(ip) = send_idx_count(ip) + i_count(ip_r)
                         send_ccc_count(ip) = send_ccc_count(ip) + 1
                      else
                         ! Set corresponding entry in send_ccc (send_idx has been set above)
                         i_ccc_cnt(ip) = i_ccc_cnt(ip) + 1
                         send_ccc(1,i_ccc_cnt(ip)) = i_cell
                         send_ccc(2,i_ccc_cnt(ip)) = i_col
                         send_ccc(3,i_ccc_cnt(ip)) = i_count(ip_r)
                      endif
                   endif
                enddo
             end if

          end do ! i_col
       end do ! i_cell

       deallocate(i_count)

    enddo ! i_k_point

    if(.not.get_size) then
       deallocate(i_idx_cnt)
       deallocate(i_ccc_cnt)
    endif

  end subroutine set_local_index_comm_arrays
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/init_comm_sparse_local_matrix_scalapack
!  NAME
!    init_comm_sparse_local_matrix_scalapack
!  SYNOPSIS
  subroutine init_comm_sparse_local_matrix_scalapack
!  PURPOSE
!    Initializes the communication for get_set_sparse_local_matrix_scalapack
!  USES
    use mpi_tasks
    use pbc_lists
    use runtime_choices
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    all index arrays for local_index communication are set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer i, ip, mpierr



    ! Allocate all count and displacement arrays

    allocate(send_idx_count(0:n_tasks-1))
    allocate(send_idx_displ(0:n_tasks-1))

    allocate(send_ccc_count(0:n_tasks-1))
    allocate(send_ccc_displ(0:n_tasks-1))

    allocate(recv_row_count(0:n_tasks-1))
    allocate(recv_row_displ(0:n_tasks-1))

    allocate(recv_ccc_count(0:n_tasks-1))
    allocate(recv_ccc_displ(0:n_tasks-1))


    ! Count how many elements of the hamiltonian go to every processor
    ! and the number of entries in the send_ccc_desc array for every processor
    ! i.e. calculate send_idx_count and send_ccc_count

    call set_local_index_comm_arrays(get_size=.true.)

    ! Calculate displacements in the arrays to be send

    send_idx_displ(0) = 0
    send_ccc_displ(0) = 0

    do ip = 1, n_tasks-1
       send_idx_displ(ip) = send_idx_displ(ip-1) + send_idx_count(ip-1)
       send_ccc_displ(ip) = send_ccc_displ(ip-1) + send_ccc_count(ip-1)
    enddo

    ! Allocate arrays for sending index and ccc

    send_idx_count_tot = sum(send_idx_count(:))
    allocate(send_idx(send_idx_count_tot))

    send_ccc_count_tot = sum(send_ccc_count(:))
    allocate(send_ccc(3,send_ccc_count_tot))

    ! Calculate send_idx and send_ccc

    call set_local_index_comm_arrays(get_size=.false.)

    ! Tell the receivers how many points will be sent in send_idx/send_ccc

    call mpi_alltoall(send_idx_count, 1, MPI_INTEGER, recv_row_count, 1, MPI_INTEGER, mpi_comm_global, mpierr)
    call mpi_alltoall(send_ccc_count, 1, MPI_INTEGER, recv_ccc_count, 1, MPI_INTEGER, mpi_comm_global, mpierr)

    ! Calculate displacements in the arrays to be received

    recv_row_displ(0) = 0
    recv_ccc_displ(0) = 0

    do ip = 1, n_tasks-1
       recv_row_displ(ip) = recv_row_displ(ip-1) + recv_row_count(ip-1)
       recv_ccc_displ(ip) = recv_ccc_displ(ip-1) + recv_ccc_count(ip-1)
    enddo

    ! Allocate arrays for receiving row numbers and ccc

    recv_row_count_tot = sum(recv_row_count(:))
    allocate(recv_row(recv_row_count_tot))
    recv_ccc_count_tot = sum(recv_ccc_count(:))
    allocate(recv_ccc(3,recv_ccc_count_tot))

    ! Please note: The receiver does not need the index where to find every matrix element
    ! (as stored in send_idx) but the row number of every element.
    ! Since send_idx is rather big, we do not allocate a new array for getting the rows
    ! but over-write send_idx and recalculate it later

    do i=1, send_idx_count_tot
       send_idx(i) = column_index_hamiltonian(send_idx(i))
    enddo

    ! Send the row numbers to the receivers

    call mpi_alltoallv(send_idx, send_idx_count, send_idx_displ, MPI_INTEGER, &
                       recv_row, recv_row_count, recv_row_displ, MPI_INTEGER, &
                       mpi_comm_global, mpierr)

    ! Recalculate send_idx

    call set_local_index_comm_arrays(get_size=.false.)

    ! Send the CCC arrray entries to the receivers

    send_ccc_count(:) = send_ccc_count(:)*3
    recv_ccc_count(:) = recv_ccc_count(:)*3
    send_ccc_displ(:) = send_ccc_displ(:)*3
    recv_ccc_displ(:) = recv_ccc_displ(:)*3

    call mpi_alltoallv(send_ccc, send_ccc_count, send_ccc_displ, MPI_INTEGER, &
                       recv_ccc, recv_ccc_count, recv_ccc_displ, MPI_INTEGER, &
                       mpi_comm_global, mpierr)

    recv_ccc_count(:) = recv_ccc_count(:)/3
    recv_ccc_displ(:) = recv_ccc_displ(:)/3


    ! The send_ccc arrays are not needed any more and immediately deallocated here

    deallocate(send_ccc)
    deallocate(send_ccc_count)
    deallocate(send_ccc_displ)

  end subroutine init_comm_sparse_local_matrix_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/solve_evp_scalapack
!  NAME
!    solve_evp_scalapack
!  SYNOPSIS
  subroutine solve_evp_scalapack( KS_eigenvalue, KS_eigenvector, i_spin )
!  PURPOSE
!    Solves the eigenvalue problem with ScaLAPACK
!  USES
    use dimensions, only: n_periodic, n_states_k
    use timing, only: number_of_loops
    use physics, only: ev_sum_change, ev_sum
    use applicable_citations, only: cite_reference
    use cg_scalapack, only: initialize_cg_prec_scalapack, &
        initialize_cg_scalapack, cg_solver_scalapack
    use elpa1_2013
    use elpa2_2013
    use localorb_io
    use mpi_tasks
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states) :: KS_eigenvalue
    real*8, dimension(n_basis, n_states) :: KS_eigenvector
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector -- Kohn-Sham eigenvectors
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8, dimension(n_basis) :: aux_eigenvalue
    integer, dimension(:), allocatable :: iwork

    integer :: info, nwork, n
    integer :: i_row, i_col
    integer, dimension(dlen_) :: ns_desc

    integer :: mpierr
    logical :: use_cg_this_cycle, cg_converged, prec_done = .false.

    character*100 :: info_str
    character*256 :: filename

    integer, save :: n_calls=0
    real*8, save :: time_1stage, time_2stage

    real*8, allocatable :: tmp(:,:)
!-------------------- for accuracy testing only:
    real*8 :: err, errmax


    write (info_str,'(2X,A,A)') &
       "Solving real symmetric generalised eigenvalue problem ", &
       "by standard ScaLAPACK."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    allocate(tmp(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp')
    tmp = 0


    n_calls = n_calls+1

    if(my_scalapack_id<npcol*nprow) then ! The main work is done only on the working set

!!       call check_ev_orthogonality_real()

       use_cg_this_cycle = .false.
       if ( use_cg.and.(number_of_loops > 1).and. &
            (ABS(ev_sum_change*hartree) < lopcg_start_tol).and. &
            (number_of_loops > initial_ev_solutions) ) use_cg_this_cycle = .true.

       if (use_cg_this_cycle .and. (.not.prec_done)) then
          call initialize_cg_prec_scalapack( ham(1,1,1), KS_eigenvalue )
          prec_done = .true.
       end if

       ! Factor the overlap matrix if necessary

       if(factor_overlap) then

          num_ovlp = num_ovlp+1
          num_ham = 0

          if (  scan( out_mat_par, 'o' )  >  0  ) then
             write(filename,'("ovlp-K",i4.4,"-O",i4.4,".out")') my_k_point,num_ovlp
             call output_matrix_real(filename,ovlp)
          endif

          ! When using the iterative solver we need to save a copy of several variables (ouch!)
          if (use_cg) call initialize_cg_scalapack( ovlp, sc_desc, mxld, mxcol, npcol, nprow, nb, mb, &
               rsrc, csrc, my_blacs_ctxt, &
               my_scalapack_id, my_scalapack_comm_work, mpi_comm_rows, mpi_comm_cols, l_col, l_row, dlen_ )

          ! Check if the overlap matrix is singular and use the scaled eigenvectors
          ! for transformation to generalized eigenvalue problem in this case,
          ! else proceed to Cholesky-decomposition.
          ! If this check is not wanted, just delete the call to diagonalize_overlap_scalapack_real.
          ! The routine will then always use Cholesky.

          call diagonalize_overlap_scalapack_real

          if(.not. use_ovlp_trafo) then

             ! Factor the overlap matrix into U**T * U using Cholesky decomposition
             ! The upper triangular matrix U is saved in ovlp

             write (info_str,'(2X,A)') 'Factoring overlap matrix with Cholesky decomposition'
             call localorb_info(info_str,use_unit,'(A)', OL_norm)

             if(use_elpa) then

                ! For some reason PDPOTRF sometimes fails for small blocksizes
                ! We use therefore own routines

                call cholesky_real_2013(n_basis, ovlp, mxld, nb, &
                      mpi_comm_rows, mpi_comm_cols)
                call invert_trm_real_2013(n_basis, ovlp, mxld, nb,&
                      mpi_comm_rows, mpi_comm_cols)

             else

                CALL PDPOTRF('U', n_basis, ovlp, 1, 1, sc_desc, info )
                if(info /= 0) call scalapack_err_exit(info,"PDPOTRF")

                ! Calculate the inverse of U and save it in ovlp

                tmp(:,:) = 0
                do n = 1, n_basis
                   if(l_row(n)>0 .and. l_col(n)>0) tmp(l_row(n),l_col(n)) = 1.
                enddo
                do n=1,n_basis,2*nb*npcol
                   nwork = 2*nb*npcol
                   if(n+nwork-1>n_basis) nwork = n_basis-n+1
                   CALL PDTRSM('L', 'U', 'N', 'N', n+nwork-1, nwork, 1.d0, &
                        ovlp, 1, 1, sc_desc, tmp, 1, n, sc_desc)
                enddo
                ovlp(:,:) = tmp(:,:)

             endif
          end if

          factor_overlap = .false.

       endif

       ! Store the k-dependent n_nonsing_ovlp in n_states_k(:)
       ! n_states is no longer reduced here, but possibly in solve_KS_eigen.f90
       n_states_k(my_k_point) = n_nonsing_ovlp

       num_ham = num_ham+1

       ! working  on the full Ham matrix is faster than working on the upper triangle
       ! therefore we set the lower triangle from the upper triangle here

       call set_full_matrix_real(ham(1,1,i_spin))

       if (  scan( out_mat_par, 'h' )  >  0  ) then
          write(filename,'("ham-K",i4.4,"-O",i4.4,"-H",i4.4,"-S",i1,".out")') my_k_point,num_ovlp,num_ham,i_spin
          call output_matrix_real(filename,ham(1,1,i_spin))
       endif

       ! Transform problem to standard eigenvalue problem

       if(use_ovlp_trafo) then

          ! Compute: ovlp_trafo**T * Ham * ovlp_trafo

          write (info_str,'(2X,A,I8)') &
               'Using scaled eigenvectors for transformation, n_nonsing_ovlp=', &
               n_nonsing_ovlp
          call localorb_info(info_str,use_unit,'(A)')

          ! Step 1: tmp = Ham * ovlp_trafo

          call pdgemm('N','N',n_basis,n_nonsing_ovlp,n_basis,1.d0,ham(1,1,i_spin),1,1,sc_desc, &
              ovlp,1,1,sc_desc,0.d0,tmp,1,1,sc_desc)

          ! Step 2: Ham = ovlp_trafo**T * tmp

          call pdgemm('T','N',n_nonsing_ovlp,n_nonsing_ovlp,n_basis,1.d0,ovlp,1,1,sc_desc, &
               tmp,1,1,sc_desc,0.d0,ham(1,1,i_spin),1,1,sc_desc)

       elseif (use_cg_this_cycle) then

          call set_full_matrix_real(ham(1,1,i_spin))
          call cg_solver_scalapack( n_nonsing_ovlp, n_states, KS_eigenvalue, &
               ham(1,1,i_spin), eigenvec_untrafo(1,1,i_spin), i_spin, cg_converged)

       end if

       if ((.not.use_cg_this_cycle) .or. (.not.cg_converged)) then

          if (use_cg_this_cycle.and.(.not.cg_converged)) then
             write (info_str,'(2X,A)') &
                  'LOPCG failed to converge. Switching to ScaLAPACK'
             call localorb_info(info_str,use_unit,'(A)')
          end  if

          if (.not.use_ovlp_trafo) then

             ! Compute: U**-T * Ham * U**-1

             write (info_str,'(2X,A)') 'Using Cholesky factor for transformation'
             call localorb_info(info_str,use_unit,'(A)', OL_norm)

             ! The tricky part here is that U**-1 is an upper triangular matrix and
             ! that only one triangle of the result is needed (since it is symmetric).
             ! So doing this in a blocked fashion like below saves half of the computations.
             ! Using pdgemm in this way is still much faster than using pdsymm

             ! Step 1: tmp = Ham * U**-1, only the upper triangle of blocks is needed.
             ! For some reason, it is faster (on the IBM regatta) to calculate
             ! the transpose U**-T * Ham and to transpose back the result.

             tmp = 0 ! For safety, must not contain garbage, especially not NaN !!!!!

             if(use_elpa) then
                call mult_at_b_real_2013('U','L',n_basis,n_basis,ovlp,mxld,&
                     ham(1,1,i_spin),mxld,nb,mpi_comm_rows,mpi_comm_cols,tmp,&
                     mxld)

             else
                do n=1,n_basis,2*nb*npcol
                   nwork = 2*nb*npcol
                   if(n+nwork-1>n_basis) nwork = n_basis-n+1
                   call pdgemm('T','N',nwork,n+nwork-1,n+nwork-1,1.d0,ovlp,1,n,sc_desc, &
                        ham(1,1,i_spin),1,1,sc_desc,0.d0,tmp,n,1,sc_desc)
                enddo
             endif

             call pdtran(n_basis,n_basis,1.d0,tmp,1,1,sc_desc,0.d0,ham(1,1,i_spin),1,1,sc_desc)
             tmp(:,:) = ham(:,:,i_spin)

             ! Step 2: ham = U**-T * tmp, only the upper triangle of blocks is needed.

             if(use_elpa) then
                call mult_at_b_real_2013('U','U',n_basis,n_basis,ovlp,mxld,tmp,&
                     mxld,nb,mpi_comm_rows,mpi_comm_cols,ham(1,1,i_spin),mxld)
             else
                do n=1,n_basis,2*nb*npcol
                   nwork = 2*nb*npcol
                   if(n+nwork-1>n_basis) nwork = n_basis-n+1
                   call pdgemm('T','N',nwork,n_basis-n+1,n+nwork-1,1.d0,ovlp,1,n,sc_desc, &
                        tmp,1,n,sc_desc,0.d0,ham(1,1,i_spin),n,n,sc_desc)
                enddo
             endif
             ! Please note: Lower part of ham contains nonsense and must not be used!

          end if

          if (  scan( out_mat_par, 's' )  >  0  ) then
             write(filename,'("sys_mat-K",i4.4,"-O",i4.4,"-H",i4.4,"-S",i1,".out")') my_k_point,num_ovlp,num_ham,i_spin
             call output_matrix_real(filename,ham(1,1,i_spin))
          endif

          ! Solve eigenvalue problem with special form of PDSYEVD

          allocate(iwork(liwork),stat=info)
          call check_allocation(info, 'iwork                         ')

          allocate(scalapack_work(len_scalapack_work),stat=info)
          call check_allocation(info, 'scalapack_work                 ')

          ! Store ham (for accuracy testing only)
          eigenvec(:,:,i_spin) = ham(:,:,i_spin)


          ! PDSYEVD might need a descriptor exactly fitting the problem dimensions
          ! (PZHEEVD needs it)
          call descinit( ns_desc, n_nonsing_ovlp, n_nonsing_ovlp, mb, nb, rsrc, csrc, &
               my_blacs_ctxt, MAX(1,mxld), info )

          if(use_elpa) then

             ! solver_method is set in runtime_choices
             if(solver_method == 0) then ! decide during runtime
                if(n_nonsing_ovlp<256) then
                   solver_method_used = 1 ! always 1-stage for samll systems
                else if(n_calls == 1) then
                   solver_method_used = 1 ! try 1-stage
                else if(n_calls == 2) then
                   solver_method_used = 2 ! try 2-stage
                endif
                ! for n_calls>2, solver_method_used is already set
             else ! solver_method is fixed
                solver_method_used = solver_method
             endif

             ! The ELPA solver needs a full matrix, so fill the lower triangle from the upper

             call set_full_matrix_real(ham(1,1,i_spin))
             call MPI_BARRIER(my_scalapack_comm_work, mpierr) ! Just to get correct timings

             if(solver_method_used /= 2) then
                ! use 1-stage solver
                call localorb_info("  Using 1-stage real ELPA solver",use_unit,'(A)',OL_norm)
                call cite_reference("ELPA")
                call solve_evp_real_2013(n_nonsing_ovlp, n_states, ham(:,:,i_spin), &
                                  mxld, aux_eigenvalue, tmp, mxld, nb, mpi_comm_rows, &
                                  mpi_comm_cols, my_scalapack_comm_work)
                time_1stage = time_evp_fwd + time_evp_back
             else
                ! use 2-stage solver
                call localorb_info("  Using 2-stage real ELPA solver",use_unit,'(A)',OL_norm)
                call cite_reference("ELPA")
                call cite_reference("ELPA-2stage")
                call solve_evp_real_2stage_2013(n_nonsing_ovlp, n_states, &
                     ham(:,:,i_spin), mxld, aux_eigenvalue, &
                     tmp, mxld, nb, mpi_comm_rows, mpi_comm_cols, &
                     my_scalapack_comm_work)
                time_2stage = time_evp_fwd + time_evp_back
             endif

             if(solver_method == 0 .and. n_nonsing_ovlp>=256 .and. n_calls == 2) then

                ! The proc with my_scalapack_id==0 decides which method to use
                if(my_scalapack_id == 0) then
                   if(time_2stage<time_1stage) then
                      solver_method_used = 2
                   else
                      solver_method_used = 1
                   endif

                   write(info_str, '(a,i5,a,f10.3,a,f10.3,a,i1,a)') '  K-Point:',my_k_point, &
                        ', Trafo times: 1-stage=',time_1stage,' s, 2-stage=',time_2stage, &
                        ' s ==> Using ',solver_method_used,'-stage solver'
                   call localorb_info(info_str,use_unit,'(A)',OL_norm)
                endif

                call MPI_BCAST( solver_method_used, 1, MPI_INTEGER, 0, my_scalapack_comm_work, mpierr)

             endif

          end if

       endif

       KS_eigenvalue(1:n_states) = aux_eigenvalue(1:n_states)

!-       Accuracy testing of solution - normally disabled

      if(.false.) then ! DISABLED
         ham(:,:,i_spin) = eigenvec(:,:,i_spin)
         call pdsymm('L','U',n_basis,n_states,1.d0,ham(1,1,i_spin),1,1,sc_desc, &
                 tmp,1,1,sc_desc,0.d0,eigenvec(1,1,i_spin),1,1,sc_desc)

         err = 0
         do i_col=1,n_states
            if(l_col(i_col)==0) cycle
            do i_row=1,n_basis
               if(l_row(i_row)>0) then
                  err = max(err,abs(tmp(l_row(i_row),l_col(i_col))*aux_eigenvalue(i_col) &
                                    -eigenvec(l_row(i_row),l_col(i_col),i_spin)))
               endif
            enddo
         enddo

         call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,my_scalapack_comm_work,info)
         if(my_scalapack_id==0) write (use_unit,*) 'MAX Eigensolution ERR:',errmax
      endif

      if (allocated(iwork)) deallocate(iwork)
      if (allocated(scalapack_work)) deallocate(scalapack_work)

      ! Backtransform eigenvectors to the original problem

      if(use_ovlp_trafo) then
         call pdgemm('N','N',n_basis,n_states,n_nonsing_ovlp,1.d0,ovlp,1,1,sc_desc, &
              tmp,1,1,sc_desc,0.d0,eigenvec(1,1,i_spin),1,1,sc_desc)
      elseif ((.not.use_cg_this_cycle).or.(.not.cg_converged)) then
         if(use_elpa) then
            call pdtran(n_basis,n_basis,1.d0,ovlp,1,1,sc_desc,0.d0,ham(1,1,i_spin),1,1,sc_desc)

            call mult_at_b_real_2013('L','N',n_basis,n_states,&
                 ham(1,1,i_spin),mxld,tmp,mxld,nb,mpi_comm_rows,&
                 mpi_comm_cols,eigenvec(1,1,i_spin),mxld)
         else
            do n=1,n_basis,2*mb*nprow
               nwork = 2*mb*nprow
               if(n+nwork-1>n_basis) nwork = n_basis-n+1
               call pdgemm('N','N',nwork,n_states,n_basis-n+1,1.d0,ovlp,n,n,sc_desc, &
                    tmp,n,1,sc_desc,0.d0,eigenvec(1,1,i_spin),n,1,sc_desc)
            enddo
         endif
      end if

!!      call check_ev_orthogonality_real()
      if (use_cg_this_cycle .and. cg_converged) then
         eigenvec(:,:,i_spin) = eigenvec_untrafo(:,:,i_spin)
      elseif (use_cg) then
         eigenvec_untrafo(:,:,i_spin) = eigenvec(:,:,i_spin)
      end if
    endif ! work only on working set

    ! if needed, get the eigenvectors back to all threads

    if (collect_eigenvectors) then
       call collect_eigenvectors_scalapack(KS_eigenvector, i_spin)
    end if

    if (n_periodic > 0 .or. packed_matrix_format /= PM_none) then
       if (my_scalapack_id /= 0) then
          KS_eigenvalue = 0.0d0
       end if
    else
       ! The procs not in working set don't have the eigenvalues yet
       call MPI_BCAST( KS_eigenvalue, n_states, MPI_REAL8, 0, &
                       my_scalapack_comm_all, mpierr)
    end if

    deallocate(tmp)

  end subroutine solve_evp_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_full_matrix_real
!  NAME
!    set_full_matrix_real
!  SYNOPSIS
  subroutine set_full_matrix_real( mat )
!  PURPOSE
!    Sets the lower half of a distributed matrix from the upper half
!  USES
    implicit none
!  ARGUMENTS
    real*8, dimension(mxld, mxcol) :: mat
!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_col, i_row
    real*8, allocatable :: tmp2(:,:)

    ! Allocate tmp2 bigger than necessary to catch overwrites in pdtran

    allocate(tmp2(mxld,mxcol+2*nb)) ! no idea whats really needed

    tmp2 = 1.d99 ! fill with guard value

    call pdtran(n_basis,n_basis,1.d0,mat,1,1,sc_desc,0.d0,tmp2,1,1,sc_desc)

    ! Check if guard value has been overwritten

!    do i_col = mxcol+1, ubound(tmp2,2)
!       if(any(tmp2(:,i_col) /= 1.d99)) then
!         print *,'***Process:',myid,'PDTRAN corrupted column:',i_col
!         print *,'   Most likely, this warning is caused by a bug in Intels MKL library.'
!         print *,'   Check aimsclub for information, and if necessary, report the problem.'
!
!         call aims_stop('Stopping execution.', 'set_full_matrix_real')
!       end if
!    enddo

    do i_col=1,n_basis-1
       if(l_col(i_col)==0) cycle
       do i_row=i_col+1,n_basis
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp2(l_row(i_row),l_col(i_col))
       enddo
    enddo

    deallocate(tmp2)

  end  subroutine set_full_matrix_real
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/output_matrix_real
!  NAME
!    output_matrix_real
!  SYNOPSIS
  subroutine output_matrix_real(filename, mat)
!  PURPOSE
!    Outputs a matrix to a file
!  USES
    use mpi_tasks
    use runtime_choices
    use synchronize_mpi_basic, only: sync_vector
    implicit none
!  ARGUMENTS
    character(len=*), intent(in) :: filename
    real*8, dimension(mxld, mxcol), intent(in) :: mat
!  INPUTS
!    o filename -- name of the file to be written
!    o mat -- matrix to be written, only upper half needs to be set
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_col, i_row, iunit
    logical od
    real*8, allocatable :: matcol(:)


    allocate(matcol(n_basis))

    if(my_scalapack_id == 0) then

       do iunit = 10, 99
          inquire(iunit, opened=od)
          if(.not.od) exit
       enddo

       if(iunit > 99) call aims_stop('No free unit found')

       if ( out_mat_par_format == 'bin' ) then
          open(iunit,file=filename,form='UNFORMATTED',action='WRITE',status='REPLACE')
          write(iunit) n_basis
       else
          open(iunit,file=filename,form='FORMATTED',action='WRITE',status='REPLACE')
          write(iunit,'(i12)') n_basis
       endif

    endif

    do i_col = 1, n_basis

       matcol(:) = 0

       if(l_col(i_col)>0) then
          do i_row = 1, n_basis
             if(l_row(i_row)>0) matcol(i_row) = mat(l_row(i_row),l_col(i_col))
          enddo
       endif

       call sync_vector(matcol, n_basis, mpi_comm=my_scalapack_comm_work)

       if(my_scalapack_id == 0) then
          if ( out_mat_par_format == 'bin' ) then
             ! Binary
             write(iunit) matcol(1:i_col)
          else
             write(iunit,'(g25.16)') matcol(1:i_col)
          endif
       endif

    enddo

    if(my_scalapack_id == 0) close(iunit)

    deallocate(matcol)

  end  subroutine output_matrix_real
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/solve_evp_scalapack_complex
!  NAME
!    solve_evp_scalapack_complex
!  SYNOPSIS
  subroutine solve_evp_scalapack_complex( KS_eigenvalue, KS_eigenvector_complex, i_spin )
!  PURPOSE
!    Solves the complex eigenvalue problem with ScaLAPACK
!  USES
    use applicable_citations, only: cite_reference
    use dimensions, only: n_states_k
    use elpa1_2013
    use elpa2_2013
    use localorb_io
    use mpi_tasks
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states) :: KS_eigenvalue
    complex*16, dimension(n_basis, n_states) :: KS_eigenvector_complex
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!    o KS_eigenvector_complex -- complex Kohn-Sham eigenvectors
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    real*8, dimension(n_basis) :: aux_eigenvalue
    integer, dimension(:), allocatable :: iwork

    integer :: info, nwork, n, mpierr, i_row, i_col
    integer, dimension(dlen_) :: ns_desc

    character*100 :: info_str

    complex*16, allocatable :: tmp_complex(:,:)

    integer, save :: n_calls=0, solver_method_used=0
    real*8, save :: time_1stage, time_2stage

    write (info_str,'(2X,A,A)') &
         "Solving hermitian generalised eigenvalue problem ", &
         "by standard ScaLAPACK."
    call localorb_info(info_str,use_unit,'(A)')

    n_calls = n_calls+1

    allocate(tmp_complex(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_complex')
    tmp_complex = 0.


    if(my_scalapack_id<npcol*nprow) then ! The main work is done only on the working set

      ! Factor the overlap matrix if necessary

      if(factor_overlap) then

        ! Check if the overlap matrix is singular and use the scaled eigenvectors
        ! for transformation to generalized eigenvalue problem in this case,
        ! else proceed to Cholesky-decomposition.
        ! If this check is not wanted, just delete the call to diagonalize_overlap_scalapack_complex.
        ! The routine will then always use Cholesky.

        call diagonalize_overlap_scalapack_complex
        if(.not.use_ovlp_trafo) then

          ! Factor the overlap matrix into U**H * U using Cholesky decomposition
          ! The upper triangular matrix U is saved in ovlp_complex

          write (info_str,'(2X,A)') 'Factoring overlap matrix with Cholesky decomposition'
          call localorb_info(info_str,use_unit,'(A)', OL_norm)

          if(use_elpa) then

            ! For some reason PDPOTRF sometimes fails for small blocksizes
            ! We use therefore own routines
            call cholesky_complex_2013(n_basis, ovlp_complex, mxld, nb,&
                  mpi_comm_rows, mpi_comm_cols)
            call invert_trm_complex_2013(n_basis, ovlp_complex, mxld, nb,&
                  mpi_comm_rows, mpi_comm_cols)

          else

            CALL PZPOTRF('U', n_basis, ovlp_complex, 1, 1, sc_desc, info )
            if(info /= 0) call scalapack_err_exit(info,"PZPOTRF")

            ! Calculate the inverse of U and save it in ovlp_complex

            tmp_complex(:,:) = 0
            do n = 1, n_basis
               if(l_row(n)>0 .and. l_col(n)>0) tmp_complex(l_row(n),l_col(n)) = 1.
            enddo
            do n=1,n_basis,2*nb*npcol
               nwork = 2*nb*npcol
               if(n+nwork-1>n_basis) nwork = n_basis-n+1
               CALL PZTRSM('L', 'U', 'N', 'N', n+nwork-1, nwork, (1.d0,0.d0), &
                    ovlp_complex, 1, 1, sc_desc, tmp_complex, 1, n, sc_desc)
            enddo
            ovlp_complex(:,:) = tmp_complex(:,:)
          endif

        end if

        factor_overlap = .false.

      endif

      ! Store the k-dependent n_nonsing_ovlp in n_states_k(:)
      ! n_states is no longer reduced here, but possibly in solve_KS_eigen.f90
      n_states_k(my_k_point) = n_nonsing_ovlp

      ! working  on the full Ham matrix is faster than working on the upper triangle
      ! therefore we set the lower triangle from the upper triangle here

      call set_full_matrix_complex(ham_complex(1,1,i_spin))

      ! Transform problem to standard eigenvalue problem

      if(use_ovlp_trafo) then

         ! Compute: ovlp_trafo**H * Ham * ovlp_trafo

         write (info_str,'(2X,A,I8)') &
              'Using scaled eigenvectors for transformation, n_nonsing_ovlp=', &
              n_nonsing_ovlp
         call localorb_info(info_str,use_unit,'(A)')

         ! Step 1: tmp = Ham * ovlp_trafo

         call pzgemm('N','N',n_basis,n_nonsing_ovlp,n_basis,(1.d0,0.d0),ham_complex(1,1,i_spin),1,1,sc_desc, &
              ovlp_complex,1,1,sc_desc,(0.d0,0.d0),tmp_complex,1,1,sc_desc)

         ! Step 2: Ham = ovlp_trafo**H * tmp

         call pzgemm('C','N',n_nonsing_ovlp,n_nonsing_ovlp,n_basis,(1.d0,0.d0),ovlp_complex,1,1,sc_desc, &
                      tmp_complex,1,1,sc_desc,(0.d0,0.d0),ham_complex(1,1,i_spin),1,1,sc_desc)
      else

         ! Compute: U**-H * Ham * U**-1

         write (info_str,'(2X,A)') 'Using Cholesky factor for transformation'
         call localorb_info(info_str,use_unit,'(A)')

         ! The tricky part here is that U**-1 is an upper triangular matrix and
         ! that only one triangle of the result is needed (since it is symmetric).
         ! So doing this in a blocked fashion like below saves half of the computations.
         ! Using pzgemm in this way is still much faster than using pzhemm

         ! Step 1: tmp = Ham * U**-1, only the upper triangle of blocks is needed.
         ! For some reason, it is faster (on the IBM regatta) to calculate
         ! the transpose U**-H * Ham and to transpose back the result.

         tmp_complex = 0 ! For safety, must not contain garbage, especially not NaN !!!!!

         if(use_elpa) then
            call mult_ah_b_complex_2013('U','L',n_basis,n_basis,ovlp_complex,&
                 mxld,ham_complex(:,:,i_spin),mxld,nb,mpi_comm_rows,&
                 mpi_comm_cols,tmp_complex,mxld)
         else
            do n=1,n_basis,2*nb*npcol
               nwork = 2*nb*npcol
               if(n+nwork-1>n_basis) nwork = n_basis-n+1
               call pzgemm('C','N',nwork,n+nwork-1,n+nwork-1,(1.d0,0.d0),ovlp_complex,1,n,sc_desc, &
                    ham_complex(1,1,i_spin),1,1,sc_desc,(0.d0,0.d0),tmp_complex,n,1,sc_desc)
            enddo
         endif

         call pztranc(n_basis,n_basis,(1.d0,0.d0),tmp_complex,1,1,sc_desc,(0.d0,0.d0),ham_complex(1,1,i_spin),1,1,sc_desc)
         tmp_complex(:,:) = ham_complex(:,:,i_spin)

         ! Step 2: ham = U**-H * tmp, only the upper triangle of blocks is needed.

         if(use_elpa) then
            call mult_ah_b_complex_2013('U','U',n_basis,n_basis,ovlp_complex,&
                 mxld,tmp_complex,mxld,nb,mpi_comm_rows,mpi_comm_cols,&
                 ham_complex(:,:,i_spin),mxld)
         else
            do n=1,n_basis,2*nb*npcol
               nwork = 2*nb*npcol
               if(n+nwork-1>n_basis) nwork = n_basis-n+1
               call pzgemm('C','N',nwork,n_basis-n+1,n+nwork-1,(1.d0,0.d0),ovlp_complex,1,n,sc_desc, &
                    tmp_complex,1,n,sc_desc,(0.d0,0.d0),ham_complex(1,1,i_spin),n,n,sc_desc)
            enddo
         endif

         ! Please note: Lower part of ham contains nonsense and must not be used!
      end if


      if(use_elpa) then

         if(solver_method == 0) then ! decide during runtime
            if(n_nonsing_ovlp<256) then
               solver_method_used = 1 ! always 1-stage for samll systems
            elseif(n_calls == 1) then
               solver_method_used = 1 ! try 1-stage
            elseif(n_calls == 2) then
               solver_method_used = 2 ! try 2-stage
            endif
            ! for n_calls>2, solver_method_used is already set
         else ! solver_method is fixed
            solver_method_used = solver_method
         endif

         ! The ELPA solver needs a full matrix, so fill the lower triangle from the upper

         call set_full_matrix_complex(ham_complex(1,1,i_spin))

         if(solver_method_used /= 2) then
            ! use 1-stage solver
            call localorb_info("  Using 1-stage complex ELPA solver", use_unit, '(A)', OL_norm)
            call cite_reference("ELPA")
            call solve_evp_complex_2013(n_nonsing_ovlp, n_states, &
                  ham_complex(:,:,i_spin), mxld, aux_eigenvalue, &
                  tmp_complex, mxld, nb, mpi_comm_rows, mpi_comm_cols)
            time_1stage = time_evp_fwd + time_evp_back
         else
            ! use 2-stage solver
            call localorb_info("  Using 2-stage complex ELPA solver", use_unit, '(A)', OL_norm)
            call cite_reference("ELPA")
            call cite_reference("ELPA-2stage")

            call solve_evp_complex_2stage_2013(n_nonsing_ovlp, n_states,&
                 ham_complex(:,:,i_spin), mxld, aux_eigenvalue, &
                 tmp_complex, mxld, nb, mpi_comm_rows, mpi_comm_cols,&
                 my_scalapack_comm_work)

            time_2stage = time_evp_fwd + time_evp_back
         endif

         if(solver_method == 0 .and. n_nonsing_ovlp>=256 .and. n_calls == 2) then

            ! The proc with my_scalapack_id==0 decides which method to use
            if(my_scalapack_id == 0) then
               if(time_2stage<time_1stage) then
                  solver_method_used = 2
               else
                  solver_method_used = 1
               endif

               write(info_str, '(a,i5,a,f10.3,a,f10.3,a,i1,a)') '  K-Point:',my_k_point, &
                     ', Trafo times: 1-stage=',time_1stage,', 2-stage=',time_2stage, &
                     ' ==> Using ',solver_method_used,'-stage solver'
               call localorb_info(info_str,use_unit,'(A)')
            endif

            call MPI_BCAST( solver_method_used, 1, MPI_INTEGER, 0, my_scalapack_comm_work, mpierr)

         endif

      else

         ! Solve eigenvalue problem with special form of PZHEEVD

         allocate(iwork(liwork),stat=info)
         call check_allocation(info, 'iwork                         ')

         allocate(scalapack_work(len_scalapack_work),stat=info)
         call check_allocation(info, 'scalapack_work                 ')

         ! PZHEEVD needs a descriptor exactly fitting the problem dimensions!
         call descinit( ns_desc, n_nonsing_ovlp, n_nonsing_ovlp, mb, nb, rsrc, csrc, &
              my_blacs_ctxt, MAX(1,mxld), info )
         CALL PZHEEVD_X('V', 'U', n_nonsing_ovlp, ham_complex(1,1,i_spin), 1, 1, ns_desc, &
              aux_eigenvalue, tmp_complex, 1, 1, ns_desc, &
              n_states, scalapack_work, lcwork, scalapack_work(2*lcwork), lrwork, iwork, liwork, info )
         if(info /= 0) call scalapack_err_exit(info,"PZHEEVD_X")

         deallocate(iwork)
         deallocate(scalapack_work)

      endif

      KS_eigenvalue(1:n_states) = aux_eigenvalue(1:n_states)

      ! Backtransform eigenvectors to the original problem

      if(use_ovlp_trafo) then
         call pzgemm('N','N',n_basis,n_states,n_nonsing_ovlp,(1.d0,0.d0),ovlp_complex,1,1,sc_desc, &
              tmp_complex,1,1,sc_desc,(0.d0,0.d0),eigenvec_complex(1,1,i_spin),1,1,sc_desc)
      else
         if(use_elpa) then
            call pztranc(n_basis, n_basis, (1.d0,0.d0), ovlp_complex, 1, 1, sc_desc, &
            (0.d0,0.d0), ham_complex(1,1,i_spin), 1, 1, sc_desc)

            call mult_ah_b_complex_2013('L','N',n_basis,n_states,&
                 ham_complex(:,:,i_spin),mxld,tmp_complex,mxld,nb,&
                 mpi_comm_rows,mpi_comm_cols,eigenvec_complex(:,:,i_spin),mxld)
         else
            do n=1,n_basis,2*mb*nprow
               nwork = 2*mb*nprow
               if(n+nwork-1>n_basis) nwork = n_basis-n+1
               call pzgemm('N','N',nwork,n_states,n_basis-n+1,(1.d0,0.d0),ovlp_complex,n,n,sc_desc, &
                        tmp_complex,n,1,sc_desc,(0.d0,0.d0),eigenvec_complex(1,1,i_spin),n,1,sc_desc)
            enddo
         endif
      end if

    endif ! work only on working set

    ! if required, get the eigenvectors back to all threads

    if (collect_eigenvectors) then
       call collect_eigenvectors_scalapack_complex(KS_eigenvector_complex, i_spin)
    end if

    ! zero the eigenvalues for all other threads than the master of each k-point
    ! this way the fortcoming sync is easier to perform
    if (my_scalapack_id /= 0) then
       KS_eigenvalue = 0.0d0
    end if

    deallocate(tmp_complex)

  end subroutine solve_evp_scalapack_complex
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_full_matrix_complex
!  NAME
!    set_full_matrix_complex
!  SYNOPSIS
  subroutine set_full_matrix_complex( mat )
!  PURPOSE
!    Sets the lower half of a distributed matrix from the upper half
!  USES
    use mpi_tasks, only: check_allocation
    implicit none
!  ARGUMENTS
    complex*16, dimension(mxld, mxcol) :: mat
!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_col, i_row
    integer info
    complex*16, allocatable :: tmp2(:,:)

    ! Allocate tmp2 bigger than necessary to catch overwrites in pdtran

    allocate(tmp2(mxld,mxcol+2*nb),stat=info) ! no idea whats really needed
    call check_allocation(info, 'tmp2')
    ! This routine is called only from the working set, so no need to check here

    call pztranc(n_basis,n_basis,(1.d0,0.d0),mat,1,1,sc_desc,(0.d0,0.d0),tmp2,1,1,sc_desc)

    do i_col=1,n_basis-1
       if(l_col(i_col)==0) cycle
       do i_row=i_col+1,n_basis
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp2(l_row(i_row),l_col(i_col))
       enddo
    enddo

    ! For safety: Make diagonal real

    do i_col=1,n_basis
       if(l_col(i_col)==0 .or. l_row(i_col)==0) cycle
       mat(l_row(i_col),l_col(i_col)) = dble(mat(l_row(i_col),l_col(i_col)))
    enddo

    deallocate(tmp2)

  end  subroutine set_full_matrix_complex
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/diagonalize_overlap_scalapack_real
!  NAME
!    diagonalize_overlap_scalapack_real
!  SYNOPSIS
  subroutine diagonalize_overlap_scalapack_real
!  PURPOSE
!    Diagonalizes the overlap matrix with ScaLAPACK
!  USES
    use dimensions
    use elpa2_2013
    use localorb_io
    use mpi_tasks
    use runtime_choices
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    if singularity in the overlap matrix is found the array ovlp is
!    overwritten by the transformation/projection array
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    real*8, dimension(:), allocatable :: eigenvalues
    integer, dimension(:), allocatable :: iwork

    real*8, allocatable :: tmp(:,:)

    integer :: info, i_row, i_col
    real*8 :: ev_sqrt
    character*100 :: info_str


    ! This routine is called only from the working set, so no need to check here

    write (info_str,'(2X,A,A)') &
         "Transforming overlap matrix with ScaLAPACK ", &
         "and checking for singularities."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    allocate(eigenvalues(n_basis),stat=info)
    call check_allocation(info, 'eigenvalues                   ')

    allocate(tmp(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp')


    ! Since the eigenvalue calculations will destroy the overlap matrix
    ! we copy the ovlp(:,:) to tmp(:,:) for the calculations
    ! The nonsingular eigenvalues must be the first ones, so calculate the eigenvalues of -ovlp!

    if(use_elpa) then

       eigenvec(:,:,1) =  -ovlp(:,:)
       call set_full_matrix_real(eigenvec) ! uses tmp internally
       tmp(:,:) = eigenvec(:,:,1)
       call solve_evp_real_2stage_2013(n_basis, n_basis, tmp, mxld, &
            eigenvalues, eigenvec(:,:,1), mxld, nb, mpi_comm_rows, &
            mpi_comm_cols, my_scalapack_comm_work)

    else

       tmp(:,:) = -ovlp(:,:)

       allocate(iwork(liwork),stat=info)
       call check_allocation(info, 'iwork                         ')

       allocate(scalapack_work(len_scalapack_work),stat=info)
       call check_allocation(info, 'scalapack_work                 ')

       call PDSYEVD('V', 'U', n_basis, tmp, 1, 1, sc_desc, eigenvalues, &
                    eigenvec, 1, 1, sc_desc, scalapack_work, lrwork, iwork, liwork, info)
       if(info /= 0) call scalapack_err_exit(info,"PDSYEVD")

       deallocate(iwork)
       deallocate(scalapack_work)

    endif


    ! We have calculated the eigenvalues of -ovlp, so invert sign of eigenvalues
    eigenvalues(1:n_basis) = -eigenvalues(1:n_basis)

    ! Get the number of nonsingular eigenvalues


    do i_row = 1, n_basis
       if(eigenvalues(i_row) < basis_threshold) exit
    end do
    n_nonsing_ovlp = i_row-1

    if (n_nonsing_ovlp < n_basis .or. force_use_ovlp_trafo) then
       if (my_scalapack_id == 0) then
         if (n_k_points.eq.1) then
          write(use_unit,'(2X,A,I8,A)') "Overlap matrix is singular:"
         else
          write(use_unit,'(2X,A,I8,A)') "k-point ", my_k_point, ": Overlap matrix is singular:"
         end if
         write(use_unit,'(2X,A,E13.6)') "| Lowest eigenvalue: ", eigenvalues(n_basis)
         write(use_unit,'(2X,A,I8,A,I8,A)') "| Using only ", n_nonsing_ovlp, &
               " out of a possible ", n_basis, " specified basis functions."
       end if

         if (.not.override_illconditioning) then
            call stop_illconditioning()
         end if

       ! Overlap matrix is not needed any more, over-write it with scaled eigenvectors
       do i_col = 1, n_nonsing_ovlp
          ev_sqrt = sqrt(eigenvalues(i_col))
          if(l_col(i_col)==0) cycle
          ovlp(:,l_col(i_col)) = eigenvec(:,l_col(i_col),1)/ev_sqrt
       end do
       use_ovlp_trafo = .true.
    else
       if ((my_scalapack_id == 0).and.(output_priority.le.OL_norm)) then
          ! don't print for MD_light, might interfere with other simultaneous outputs by other tasks...
         if (n_k_points.eq.1) then
          write(use_unit,'(2X,A)') &
               ": Overlap matrix is nonsingular, using Cholesky decomposition"
         else
          write(use_unit,'(2X,A,I8,A)') "k-point ", my_k_point, &
               ": Overlap matrix is nonsingular, using Cholesky decomposition"
         end if
          write(use_unit,'(2X,A,E13.6)') "| Lowest eigenvalue: ", eigenvalues(n_basis)
          ! Add a warning if needed
          if (eigenvalues(n_basis).le.1d-5) then
            write(use_unit,'(1X,A)') "* Warning - overlap matrix is near-singular!"
            write(use_unit,'(1X,A)') "* Consider using a larger value of basis_threshold to ensure numerical stability!"

            if (.not.override_illconditioning) then
              call stop_illconditioning()
            end if
          end if

       end if
       use_ovlp_trafo = .false.
   end if

   deallocate(eigenvalues)
   deallocate(tmp)

  end subroutine diagonalize_overlap_scalapack_real
!******
!-----------------------------------------------------------------------------------
!******
!****s* scalapack_wrapper/diagonalize_overlap_scalapack_complex
!  NAME
!    diagonalize_overlap_scalapack_complex
!  SYNOPSIS
  subroutine diagonalize_overlap_scalapack_complex
!  PURPOSE
!    Diagonalizes the complex overlap matrix with ScaLAPACK
!  USES
    use dimensions
    use elpa2_2013
    use localorb_io
    use mpi_tasks
    use runtime_choices
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    if singularity in the overlap matrix is found the array ovlp_complex is
!    overwritten by the transformation/projection array
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    real*8, dimension(:), allocatable :: eigenvalues
    integer, dimension(:), allocatable :: iwork

    complex*16, allocatable :: tmp_complex(:,:)

    integer :: info, i_row, i_col
    real*8 :: ev_sqrt
    character*100 :: info_str


    ! This routine is called only from the working set, so no need to check here

    write (info_str,'(2X,A,A)') &
         "Transforming overlap matrix with ScaLAPACK ", &
         "and checking for singularities."
    call localorb_info(info_str,use_unit,'(A)')

    allocate(eigenvalues(n_basis),stat=info)
    call check_allocation(info, 'eigenvalues                   ')

    allocate(tmp_complex(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_complex')


    ! Since the eigenvalue calculations will destroy the overlap matrix
    ! we copy the ovlp_complex(:,:) to tmp_complex(:,:) for the calculations.
    ! The nonsingular eigenvalues must be the first ones, so calculate the eigenvalues of -ovlp_complex!

    if(use_elpa) then
       call MPI_BARRIER(my_scalapack_comm_work, info) ! Just to get correct timings

       eigenvec_complex(:,:,1) =  -ovlp_complex(:,:)
       call set_full_matrix_complex(eigenvec_complex) ! uses tmp_complex internally
       tmp_complex(:,:) = eigenvec_complex(:,:,1)

       call solve_evp_complex_2stage_2013(n_basis, n_basis, tmp_complex, mxld, &
            eigenvalues, eigenvec_complex(:,:,1), mxld, nb, mpi_comm_rows, &
            mpi_comm_cols, my_scalapack_comm_work)
       call MPI_BARRIER(my_scalapack_comm_work, info) ! Just to get correct timings

    else

       tmp_complex(:,:) = -ovlp_complex(:,:)

       allocate(iwork(liwork),stat=info)
       call check_allocation(info, 'iwork                         ')

       allocate(scalapack_work(len_scalapack_work),stat=info)
       call check_allocation(info, 'scalapack_work                 ')

       CALL PZHEEVD('V', 'U', n_basis, tmp_complex, 1, 1, sc_desc, &
            eigenvalues, eigenvec_complex, 1, 1, sc_desc, &
            scalapack_work, lcwork, scalapack_work(2*lcwork), lrwork, iwork, liwork, info )
       if(info /= 0) call scalapack_err_exit(info,"PZHEEVD")

       deallocate(iwork)
       deallocate(scalapack_work)

    endif

    ! We have calculated the eigenvalues of -ovlp, so invert sign of eigenvalues
    eigenvalues(1:n_basis) = -eigenvalues(1:n_basis)

    ! Get the number of nonsingular eigenvalues

    do i_row = 1, n_basis
      if(eigenvalues(i_row) < basis_threshold) exit
    end do
    n_nonsing_ovlp = i_row-1

    if (n_nonsing_ovlp < n_basis .or. force_use_ovlp_trafo) then
      if (my_scalapack_id == 0) then
        write(use_unit,'(2X,A,I8,A)') "k-point ", my_k_point, ": Overlap matrix is singular:"
        write(use_unit,'(2X,A,E13.6)') "| Lowest eigenvalue: ", eigenvalues(n_basis)
        write(use_unit,'(2X,A,I8,A,I8,A)') "| Using only ", n_nonsing_ovlp, &
              " out of a possible ", n_basis, " specified basis functions."
      end if

         if (.not.override_illconditioning) then
            call stop_illconditioning()
         end if

      ! Overlap matrix is not needed any more, over-write it with scaled eigenvectors
      do i_col = 1, n_nonsing_ovlp
        ev_sqrt = sqrt(eigenvalues(i_col))
        if(l_col(i_col)==0) cycle
        ovlp_complex(:,l_col(i_col)) = eigenvec_complex(:,l_col(i_col),1)/ev_sqrt
      end do
      use_ovlp_trafo = .true.
    else
       if ((myid == 0).and.(output_priority.le.OL_norm)) then
          write(use_unit,'(2X,A)') &
               "k-point 1: Overlap matrix is nonsingular, using Cholesky decomposition"
          write(use_unit,'(2X,A,E13.6)') "| Lowest eigenvalue: ", eigenvalues(n_basis)
       end if
       if (my_scalapack_id == 0) then
         if (eigenvalues(n_basis).le.1d-5) then
            ! issue a warning for too low eigenvalues anyway!
            write(use_unit,'(1X,A,I8,A)') "* Warning - overlap matrix is near-singular on k-point ", my_k_point, " ."
            write(use_unit,'(1X,A,E13.6)') "* Lowest eigenvalue: ", eigenvalues(n_basis)
            write(use_unit,'(1X,A)') "* Consider using a larger value of basis_threshold to ensure numerical stability!"

            if (.not.override_illconditioning) then
              call stop_illconditioning()
            end if
         end if
       end if
       use_ovlp_trafo = .false.
    end if

    deallocate(eigenvalues)
    deallocate(tmp_complex)

  end subroutine diagonalize_overlap_scalapack_complex
!******

!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/extrapolate_dm_scalapack
!  NAME
!    extrapolate_dm_scalapack
!  SYNOPSIS
  subroutine extrapolate_dm_scalapack
!  PURPOSE
!    Extrapolates the density matrix for a new overlap matrix.
!  USE
  use elsi_wrapper, only: eh_scf,aims_elsi_extrapolate_dm
  use runtime_choices, only: orthonormalize_evs
!  SOURCE
  implicit none

  integer :: i_spin

  if(.not. orthonormalize_evs) then
     return
  end if

  ! Work must be done only on the working set
  if(my_scalapack_id < npcol*nprow) then
     if(real_eigenvectors) then
        if(.not. full_ovlp_ready) then
           call set_full_matrix_real(ovlp)
        end if

        do i_spin = 1,n_spin
           call aims_elsi_extrapolate_dm(eh_scf,ovlp,ham(:,:,i_spin))
        end do
     else ! Complex eigenvectors
        if(.not. full_ovlp_ready) then
           call set_full_matrix_complex(ovlp_complex)
        end if

        do i_spin = 1,n_spin
           call aims_elsi_extrapolate_dm(eh_scf,ovlp_complex,&
                ham_complex(:,:,i_spin))
        end do
     end if
  end if

  full_ovlp_ready = .true.

end subroutine extrapolate_dm_scalapack
!******

!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/normalize_eigenvectors_scalapack
!  NAME
!    normalize_eigenvectors_scalapack
!  SYNOPSIS
  subroutine normalize_eigenvectors_scalapack
!  PURPOSE
!    Normalizes the eigenvectors with ScaLAPACK
!  USES
    use localorb_io
    use runtime_choices
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    eigevectors in arrays eigenvec and eigenvec_stored or
!    eigenvec_complex and eigenvec_complex_stored are normalized to unity
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    character*100 :: info_str

    if (.not. orthonormalize_evs) return

    write (info_str,'(2X,A,A)') "Normalizing  ScaLAPACK eigenvectors"
    call localorb_info(info_str,use_unit,'(A)')


    if (real_eigenvectors) then
       call  normalize_eigenvectors_scalapack_real(eigenvec)
       if(allocated(eigenvec_stored)) then
         ! VB: It is not clear to me why this needs to be orthonormalized at all, ever.
         !     The "if" is here to catch the situation in MD, where the stored eigenvectors do
         !     not even exist.
         ! JW: The overlap has changed.  Please note that this is no orthogonalization,
         !     just a normalization
         call  normalize_eigenvectors_scalapack_real(eigenvec_stored)
       end if
    else
       call  normalize_eigenvectors_scalapack_complex(eigenvec_complex)
       if(allocated(eigenvec_complex_stored)) then
         call  normalize_eigenvectors_scalapack_complex(eigenvec_complex_stored)
       end if
   end if


  end subroutine normalize_eigenvectors_scalapack
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/normalize_eigenvectors_scalapack_real
!  NAME
!    normalize_eigenvectors_scalapack_real
!  SYNOPSIS
subroutine normalize_eigenvectors_scalapack_real(eigenvector)
!  PURPOSE
!    Normalizes the eigenvectors with ScaLAPACK, real version
!  USES
  use localorb_io
  use mpi_tasks
  use runtime_choices
  use synchronize_mpi_basic, only: sync_vector
  implicit none
!  ARGUMENTS
  real*8:: eigenvector(mxld, mxcol,n_spin)
!  INPUTS
!    o eigenvector -- the ScaLAPACK eigenvectors
!  OUTPUT
!    eigevectors in array eigenvector are normalized to unity
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  character*100 :: info_str

  real*8, allocatable, dimension(:) :: r_ii
  real*8, allocatable, dimension(:,:) :: work2
  real*8, allocatable, dimension(:,:) :: work
!  real*8, allocatable, dimension(:,:) :: matrix_temp
  integer:: i_spin,  i_col, info


  if (.not. orthonormalize_evs) return

  ! The work in this routine must be done only on the working set
  if(my_scalapack_id>=npcol*nprow) return

  allocate(work(mxld, mxcol),stat=info)
  call check_allocation(info, 'work                          ')


  allocate(work2(mxld, mxcol),stat=info)
  call check_allocation(info, 'work2                         ')

  allocate(r_ii(n_states),stat=info)
  call check_allocation(info, 'r_ii                          ')


  do i_spin = 1, n_spin, 1

     work2 = ovlp

     ! complete lower part of ovlp in work2
     call set_full_matrix_real(work2)

     ! work = work2 . eigenvector
     call PdGEMM( 'N', 'N', n_basis , n_states, n_basis, 1.d0, work2, 1, 1, sc_desc, &
         eigenvector(1,1,i_spin), 1, 1, sc_desc, 0.d0, work, 1, 1, sc_desc )

     ! work2 = eigenvector^T . work
     call PdGEMM( 'T', 'N', n_states , n_states, n_basis, 1.d0, &
     eigenvector(1,1,i_spin), 1, 1, sc_desc, work, 1, 1, sc_desc,  0.d0, work2, 1, 1, sc_desc )

     r_ii = 0.d0

     do i_col = 1, n_states
        if(l_col(i_col) > 0 .and. l_row(i_col) > 0)then

           if(work2(l_row(i_col), l_col(i_col)) <= 0.) then
              write(info_str, '(a,g25.15)') '*** Warning normalize_eigenvectors_scalapack_real: norm^2 = ', &
                    work2(l_row(i_col), l_col(i_col))
              call localorb_info(info_str,use_unit,'(A)')
              r_ii(i_col) = 1.
           else
              r_ii(i_col) =  sqrt(work2(l_row(i_col), l_col(i_col)))
           endif

        end if
     end do

     call sync_vector(r_ii, n_states, my_scalapack_comm_work)

     do i_col = 1, n_states
        if(l_col(i_col) > 0)then

           eigenvector(:,l_col(i_col),i_spin) = eigenvector(:,l_col(i_col),i_spin)/r_ii(i_col)

        end if
     end do
  end do

  deallocate(work)
  deallocate(work2)
  deallocate(r_ii)


end subroutine normalize_eigenvectors_scalapack_real
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/normalize_eigenvectors_scalapack_complex
!  NAME
!    normalize_eigenvectors_scalapack_complex
!  SYNOPSIS
subroutine normalize_eigenvectors_scalapack_complex(eigenvector)
!  PURPOSE
!    Normalizes the eigenvectors with ScaLAPACK, real version
!  USES
  use localorb_io
  use mpi_tasks
  use synchronize_mpi_basic, only: sync_vector
  use runtime_choices
  implicit none
!  ARGUMENTS
  complex*16:: eigenvector(mxld, mxcol,n_spin)
!  INPUTS
!    o eigenvector -- the ScaLAPACK eigenvectors
!  OUTPUT
!    eigevectors in array eigenvector are normalized to unity
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  character*100 :: info_str

  real*8, allocatable, dimension(:) :: r_ii
  complex*16, allocatable, dimension(:,:) :: work2
  complex*16, allocatable, dimension(:,:) :: work
!  real*8, allocatable, dimension(:,:) :: matrix_temp
  integer:: i_spin,  i_col, info


  if (.not. orthonormalize_evs) return

  ! The work in this routine must be done only on the working set
  if(my_scalapack_id>=npcol*nprow) return

  allocate(work(mxld, mxcol),stat=info)
  call check_allocation(info, 'work                          ')

  allocate(work2(mxld, mxcol),stat=info)
  call check_allocation(info, 'work2                         ')

  allocate(r_ii(n_states),stat=info)
  call check_allocation(info, 'r_ii                          ')

  do i_spin = 1, n_spin, 1

     work2 = ovlp_complex

     ! complete lower part of ovlp_complex in work2
     call set_full_matrix_complex(work2)

     call PzGEMM( 'N', 'N', n_basis , n_states, n_basis, (1.d0,0.d0), work2, 1, 1, sc_desc, &
         eigenvector(1,1,i_spin), 1, 1, sc_desc, (0.d0,0.d0), work, 1, 1, sc_desc )


     call PzGEMM( 'C', 'N', n_states , n_states, n_basis, (1.d0,0.d0), &
          eigenvector(1,1,i_spin), 1, 1, sc_desc, work, 1, 1, sc_desc,  (0.d0,0.d0), work2, 1, 1, sc_desc )


     r_ii = 0.d0


     do i_col = 1, n_states
        if(l_col(i_col) > 0 .and. l_row(i_col) > 0)then

           if(dble(work2(l_row(i_col), l_col(i_col))) <= 0.) then
              write(info_str, '(a,g25.15)') '*** Warning normalize_eigenvectors_scalapack_complex: norm^2 = ', &
                    dble(work2(l_row(i_col), l_col(i_col)))
              call localorb_info(info_str,use_unit,'(A)')
              r_ii(i_col) = 1.
           else
              r_ii(i_col) =  sqrt(dble(work2(l_row(i_col), l_col(i_col))))
           endif

        end if
     end do

     call sync_vector(r_ii, n_states, my_scalapack_comm_work)

     do i_col = 1, n_states
        if(l_col(i_col) > 0)then

           eigenvector(:,l_col(i_col),i_spin) = eigenvector(:,l_col(i_col),i_spin)/r_ii(i_col)

        end if
     end do
  end do



  deallocate(work)
  deallocate(work2)
   deallocate(r_ii)


end subroutine normalize_eigenvectors_scalapack_complex
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/orthonormalize_eigenvectors_scalapack_real
!  NAME
!    orthonormalize_eigenvectors_scalapack_real
!  SYNOPSIS
subroutine orthonormalize_eigenvectors_scalapack_real(KS_eigenvector)
!  PURPOSE
!    Orthonormalizes the eigenvectors with ScaLAPACK, real version
!    This routine uses a modified Gram-Schmidt orthogonalization
!    It needs two additional work matrices but is fast
!    because most work is done in matrix-matrix-products
!
!    This routine must be called after the overlap matrix is set
!    but before it is factored !!!!
!
!  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use elsi_wrapper, only: eh_scf, aims_elsi_orthonormalize_ev
  use localorb_io, only: localorb_info, use_unit, OL_norm
  use runtime_choices, only: use_elsi, frozen_core_scf, collect_eigenvectors
  implicit none
!  ARGUMENTS
  real*8 :: KS_eigenvector(n_basis,n_states,n_spin)
!  INPUTS
!    -
!  OUTPUT
!    KS_eigenvector as defined in physics.f90 . Note that if scalapack is used,
!    there can only be one k-point per task anyway. So, the n_k_points_task
!    dimension of physics.f90 is ignored here as it is implicitly 1.
!    KS_eigenvector is only provided if collect_eigenvectors is requested.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  real*8, allocatable :: prod(:)
  real*8, allocatable :: ovlp_ev(:,:)
  real*8, allocatable :: work(:,:)
  real*8 :: dotprod
  real*8 :: fact
  integer :: i_spin
  integer :: i_col
  integer :: n_block
  integer :: i_done
  character*200 :: info_str

  if(use_elsi .and. .not. frozen_core_scf) then
     do i_spin = 1,n_spin
        ! The main work is done only on the working set
        if(my_scalapack_id < npcol*nprow) then
           call aims_elsi_orthonormalize_ev(eh_scf,ovlp,eigenvec(:,:,i_spin))
        end if

        ! if needed, get the eigenvectors back to all threads
        if (collect_eigenvectors) then
           call collect_eigenvectors_scalapack(KS_eigenvector(:,:,i_spin),&
                i_spin)
        end if
     end do
  else
     write(info_str,'(2X,A)') "Orthonormalizing eigenvectors"
     call localorb_info(info_str,use_unit,'(A)',OL_norm)

     n_block = npcol*nb

     call aims_allocate(prod,mxld,"prod")
     call aims_allocate(ovlp_ev,mxld,mxcol,"ovlp_ev")
     call aims_allocate(work,mxld,mxcol,"work")

     do i_spin = 1,n_spin
        ! The main work is done only on the working set
        if(my_scalapack_id < npcol*nprow) then
           ! ovlp_ev = ovlp * eigenvec
           call pdsymm('L','U',n_basis,n_states,1.d0,ovlp,1,1,sc_desc,&
                eigenvec(1,1,i_spin),1,1,sc_desc,0.d0,ovlp_ev,1,1,sc_desc)

           ! Number of vectors which are orthogonalized against all others
           i_done = 0

           do i_col = 1,n_states
              if(i_col > i_done+1) then
                 ! Build the matrix-dot product of eigenvec(i_col) with all
                 ! eigenvectors against which eigenvec(i_col) is not yet
                 ! orthogonalized, these are the ones from i_done+1 to i_col-1
                 call pdgemv('T',n_basis,i_col-1-i_done,1.d0,ovlp_ev,1,&
                      i_done+1,sc_desc,eigenvec(1,1,i_spin),1,i_col,sc_desc,1,&
                      0.d0,prod,1,1,sc_desc,1)

                 ! Orthogonalize eigenvec against the others, keep
                 ! ovlp * eigenvec up to date
                 ! eigenvec(i_col) = eigenvec(i_col)
                 !                   - eigenvec(i_done+1:i_col-1) * prod
                 call pdgemv('N',n_basis,i_col-1-i_done,-1.d0,&
                      eigenvec(1,1,i_spin),1,i_done+1,sc_desc,prod,1,1,sc_desc,&
                      1,1.d0,eigenvec(1,1,i_spin),1,i_col,sc_desc,1)

                 ! The same for ovlp * eigenvec
                 call pdgemv('N',n_basis,i_col-1-i_done,-1.d0,ovlp_ev,1,&
                      i_done+1,sc_desc,prod,1,1,sc_desc,1,1.d0,ovlp_ev,1,i_col,&
                      sc_desc,1)
              end if

              ! Now eigenvec(i_col) is orthogonalized against all previous ones,
              ! normalize it
              ! dotprod = ovlp_ev(i_col)**T * eigenvec(i_col)
              dotprod = 0.d0

              call pddot(n_basis,dotprod,ovlp_ev,1,i_col,sc_desc,1,&
                   eigenvec(1,1,i_spin),1,i_col,sc_desc,1)

              if(l_col(i_col) > 0) then
                 ! eigenvec(i_col) = eigenvec(i_col)/sqrt(dotprod)
                 if(dotprod > 0) then
                    fact = 1.d0/sqrt(dotprod)
                 else
                    fact = 1.d0 ! for safety only, should never happen!
                 end if

                 eigenvec(:,l_col(i_col),i_spin) &
                    = eigenvec(:,l_col(i_col),i_spin)*fact
                 ovlp_ev(:,l_col(i_col)) = ovlp_ev(:,l_col(i_col))*fact
              end if

              ! If i_col-i_done reaches block size, orthogonalize the vectors
              ! i_done+1 .. i_col against all following ones. This can be done
              ! with matrix-matrix operations
              if(i_col-i_done == n_block .and. i_col < n_states) then
                 ! Build the matrix-dot product of eigenvec(i_done+1..i_col)
                 ! with eigenvec(i_col+1 .. n_states)
                 call pdgemm('T','N',n_block,n_states-i_col,n_basis,1.d0,&
                      ovlp_ev,1,i_done+1,sc_desc,eigenvec(1,1,i_spin),1,&
                      i_col+1,sc_desc,0.d0,work,1,i_col+1,sc_desc)

                 ! Orthogonalize
                 call pdgemm('N','N',n_basis,n_states-i_col,n_block,-1.d0,&
                      eigenvec(1,1,i_spin),1,i_done+1,sc_desc,work,1,i_col+1,&
                      sc_desc,1.d0,eigenvec(1,1,i_spin),1,i_col+1,sc_desc)

                 i_done = i_done+n_block
              end if
           end do ! i_col
        end if ! work only on working set

        ! if needed, get the eigenvectors back to all threads
        if (collect_eigenvectors) then
           call collect_eigenvectors_scalapack(KS_eigenvector(:,:,i_spin),&
                i_spin)
        end if
     end do ! i_spin

     call aims_deallocate(prod,"prod")
     call aims_deallocate(ovlp_ev,"ovlp_ev")
     call aims_deallocate(work,"work")
  end if ! use_elsi

end subroutine orthonormalize_eigenvectors_scalapack_real
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/orthonormalize_eigenvectors_scalapack_real
!  NAME
!    orthonormalize_eigenvectors_scalapack_real_GS
!  SYNOPSIS
subroutine orthonormalize_eigenvectors_scalapack_real_GS()
!  PURPOSE
!    Orthonormalizes the eigenvectors with ScaLAPACK, real version
!    This routine uses the classical Gram-Schmidt orthogonalization
!    It needs no additional memory (only the 2 1d arrays) but is slow
!    because it only uses matrix vector products
!  USES
  use localorb_io
  use mpi_tasks
  implicit none
!  ARGUMENTS
!    -
!  INPUTS
!    -
!  OUTPUT
!    -
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  real*8, allocatable, dimension(:) :: work1, work2
  real*8 :: dotprod, fact

  integer:: i_spin, i_col, info

    character*200 :: info_str



  ! The work in this routine must be done only on the working set
  if(my_scalapack_id>=npcol*nprow) return

  allocate(work1(mxld),stat=info)
  call check_allocation(info, 'work1                         ')

  allocate(work2(mxld),stat=info)
  call check_allocation(info, 'work2                         ')

  write(info_str,'(2X,A)') "Orthonormalizing eigenvectors"
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  do i_spin = 1, n_spin, 1

     do i_col = 1, n_states

        if(i_col > 1) then

           ! Orthogonalize against eigenvectors 1 .. i_col-1

           ! work1 = ovlp * eigenvec(i_col)

           call PDSYMV('U', n_basis, 1.0d0, ovlp, 1, 1, sc_desc, &
                       eigenvec(1,1,i_spin), 1, i_col, sc_desc, 1, &
                       0.0d0, work1, 1, 1, sc_desc, 1)

           ! work2 = eigenvec(1:i_col-1)**T * work1

           call PDGEMV('T', n_basis, i_col-1, 1.0d0, eigenvec(1,1,i_spin), 1, 1, sc_desc, &
                       work1, 1, 1, sc_desc, 1, &
                       0.0d0, work2, 1, 1, sc_desc, 1)

           ! eigenvec(i_col) = eigenvec(i_col) - eigenvec(1:i_col-1) * work2

           call PDGEMV('N', n_basis, i_col-1, -1.0d0, eigenvec(1,1,i_spin), 1, 1, sc_desc, &
                       work2, 1, 1, sc_desc, 1, &
                       1.0d0, eigenvec(1,1,i_spin), 1, i_col, sc_desc, 1)
        endif

        ! Normalize eigenvec(i_col)

        ! work1 = ovlp * eigenvec(i_col)

        call PDSYMV('U', n_basis, 1.0d0, ovlp, 1, 1, sc_desc, &
                    eigenvec(1,1,i_spin), 1, i_col, sc_desc, 1, &
                    0.0d0, work1, 1, 1, sc_desc, 1)

        ! dotprod = work1**T * eigenvec(i_col)

        dotprod = 0
        call PDDOT(n_basis, dotprod, work1, 1, 1, sc_desc, 1, eigenvec(1,1,i_spin), 1, i_col, sc_desc, 1)

        if(l_col(i_col) > 0) then

           ! eigenvec(i_col) = eigenvec(i_col)/sqrt(dotprod)

           if(dotprod>0) then
              fact = 1.d0/sqrt(dotprod)
           else
              fact = 1.d0 ! for safety only, should never happen!
           endif

           eigenvec(:,l_col(i_col),i_spin) = eigenvec(:,l_col(i_col),i_spin)*fact

        endif

     enddo ! i_col

  enddo ! i_spin


  ! TEST only, may be removed later:
  ! call check_ev_orthogonality_real

  deallocate(work1)
  deallocate(work2)

end subroutine orthonormalize_eigenvectors_scalapack_real_GS
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/check_ev_orthogonality_real
!  NAME
!    check_ev_orthogonality_real
!  SYNOPSIS
subroutine check_ev_orthogonality_real()
!  PURPOSE
!    Checks the orthogonality of the eigenvectors
!    This is a test routine only!!!!!!!!!!!!
!  USES
    use localorb_io
    use mpi_tasks, only: myid
    implicit none

   integer i_spin, i_col, i_row
   real*8 :: amax, dmax
   real*8, allocatable :: work1(:,:), work2(:,:)

   allocate(work1(mxld,mxcol))
   allocate(work2(mxld,mxcol))

   ! The work in this routine must be done only on the working set
   if(my_scalapack_id>=npcol*nprow) return

   i_spin = 1

   call PDSYMM('L', 'U', n_basis, n_states, 1.0d0, ovlp, 1, 1, sc_desc, &
               eigenvec(1,1,i_spin), 1, 1, sc_desc, &
               0.0d0, work1, 1, 1, sc_desc)

   call PDGEMM('T', 'N', n_states, n_states, n_basis, &
               1.0d0, eigenvec(1,1,i_spin), 1, 1, sc_desc, &
               work1, 1, 1, sc_desc, &
               0.0d0, work2, 1, 1, sc_desc)

   dmax = 0.
   amax = 0.

   do i_col = 1, n_states
   do i_row = 1, n_states
      if(l_row(i_row) > 0 .and. l_col(i_col)>0) then
         if(i_row==i_col) then
            ! diagonal element
            dmax = max(dmax,abs(1.0-work2(l_row(i_row),l_col(i_col))))
         else
            ! off diagonal element
            amax = max(amax,abs(work2(l_row(i_row),l_col(i_col))))
         endif
      endif
   enddo
   enddo

   write(use_unit,'(A,I4,A,F10.5,F10.5)') 'check_ev_orthogonality_real, ID: ',myid,' Errors: ',amax, dmax
   deallocate(work1)
   deallocate(work2)

end subroutine check_ev_orthogonality_real
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/orthonormalize_eigenvectors_scalapack_complex
!  NAME
!    orthonormalize_eigenvectors_scalapack_complex
!  SYNOPSIS
subroutine orthonormalize_eigenvectors_scalapack_complex(KS_eigenvector_complex)
!  PURPOSE
!    Orthonormalizes the eigenvectors with ScaLAPACK, complex version
!    This routine uses a modified Gram-Schmidt orthogonalization
!    It needs two additional work matrices but is fast
!    because most work is done in matrix-matrix-products
!
!    This routine must be called after the overlap matrix is set
!    but before it is factored !!!!
!
!  USES
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use elsi_wrapper, only: eh_scf, aims_elsi_orthonormalize_ev
  use localorb_io, only: localorb_info, use_unit, OL_norm
  use runtime_choices, only: use_elsi, frozen_core_scf, collect_eigenvectors
  implicit none
!  ARGUMENTS
  complex*16 :: KS_eigenvector_complex(n_basis,n_states,n_spin)
!  INPUTS
!    -
!  OUTPUT
!    KS_eigenvector_complex as defined in physics.f90 . Note that if scalapack
!    is used, there can only be one k-point per task anyway. So, the
!    n_k_points_task dimension of physics.f90 is ignored here as it is
!    implicitly 1. KS_eigenvector_complex is only provided if
!    collect_eigenvectors is requested.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

  complex*16, allocatable :: prod(:)
  complex*16, allocatable :: ovlp_ev(:,:)
  complex*16, allocatable :: work(:,:)
  complex*16 :: dotprod
  real*8 :: fact
  integer :: i_spin
  integer :: i_col
  integer :: n_block
  integer :: i_done
  character*200 :: info_str

  if(use_elsi .and. .not. frozen_core_scf) then
     do i_spin = 1,n_spin
        ! The main work is done only on the working set
        if(my_scalapack_id < npcol*nprow) then
           call aims_elsi_orthonormalize_ev(eh_scf,ovlp_complex,&
                eigenvec_complex(:,:,i_spin))
        end if

        ! if needed, get the eigenvectors back to all threads
        if (collect_eigenvectors) then
           call collect_eigenvectors_scalapack_complex(&
                KS_eigenvector_complex(:,:,i_spin),i_spin)
        end if
     end do
  else
     write(info_str,'(2X,A)') "Orthonormalizing eigenvectors"
     call localorb_info(info_str,use_unit,'(A)',OL_norm)

     n_block = npcol*nb

     call aims_allocate(prod,mxld,"prod")
     call aims_allocate(ovlp_ev,mxld,mxcol,"ovlp_ev")
     call aims_allocate(work,mxld,mxcol,"work")

     do i_spin = 1,n_spin
        ! The main work is done only on the working set
        if(my_scalapack_id < npcol*nprow) then
           ! ovlp_ev = ovlp * eigenvec
           call pzhemm('L','U',n_basis,n_states,(1.d0,0.d0),ovlp_complex,1,1,&
                sc_desc,eigenvec_complex(1,1,i_spin),1,1,sc_desc,(0.d0,0.d0),&
                ovlp_ev,1,1,sc_desc)

           ! Number of vectors which are orthogonalized against all others
           i_done = 0

           do i_col = 1,n_states
              if(i_col > i_done+1) then
                 ! Build the matrix-dot product of eigenvec(i_col) with all
                 ! eigenvectors against which eigenvec(i_col) is not yet
                 ! orthogonalized, these are the ones from i_done+1 to i_col-1
                 call pzgemv('C',n_basis,i_col-1-i_done,(1.d0,0.d0),ovlp_ev,1,&
                      i_done+1,sc_desc,eigenvec_complex(1,1,i_spin),1,i_col,&
                      sc_desc,1,(0.d0,0.d0),prod,1,1,sc_desc,1)

                 ! Orthogonalize eigenvec against the others, keep
                 ! ovlp * eigenvec up to date
                 ! eigenvec(i_col) = eigenvec(i_col)
                 !                   - eigenvec(i_done+1:i_col-1) * prod
                 call pzgemv('N',n_basis,i_col-1-i_done,(-1.d0,0.d0),&
                      eigenvec_complex(1,1,i_spin),1,i_done+1,sc_desc,prod,1,1,&
                      sc_desc,1,(1.d0,0.d0),eigenvec_complex(1,1,i_spin),1,&
                      i_col,sc_desc,1)

                 ! The same for ovlp * eigenvec
                 call pzgemv('N',n_basis,i_col-1-i_done,(-1.0d0,0.0d0),ovlp_ev,&
                      1,i_done+1,sc_desc,prod,1,1,sc_desc,1,(1.0d0,0.0d0),&
                      ovlp_ev,1,i_col,sc_desc,1)
              end if

              ! Now eigenvec(i_col) is orthogonalized against all previous ones,
              ! normalize it
              ! dotprod = ovlp_ev(i_col)**T * eigenvec(i_col)
              dotprod = 0
              call pzdotc(n_basis,dotprod,ovlp_ev,1,i_col,sc_desc,1,&
                   eigenvec_complex(1,1,i_spin),1,i_col,sc_desc,1)

              if(l_col(i_col) > 0) then
                 ! eigenvec(i_col) = eigenvec(i_col)/sqrt(dotprod)
                 if(dble(dotprod) > 0) then
                    fact = 1.d0/sqrt(dble(dotprod))
                 else
                    fact = 1.d0 ! for safety only, should never happen!
                 end if

                 eigenvec_complex(:,l_col(i_col),i_spin) &
                    = eigenvec_complex(:,l_col(i_col),i_spin)*fact
                 ovlp_ev(:,l_col(i_col)) = ovlp_ev(:,l_col(i_col))*fact
              end if

              ! If i_col-i_done reaches block size, orthogonalize the vectors
              ! i_done+1 .. i_col against all following ones. This can be done
              ! with matrix-matrix operations
              if(i_col-i_done == n_block .and. i_col < n_states) then
                 ! Build the matrix-dot product of eigenvec(i_done+1..i_col)
                 ! with eigenvec(i_col+1 .. n_states)
                 call pzgemm('C','N',n_block,n_states-i_col,n_basis,&
                      (1.d0,0.d0),ovlp_ev,1,i_done+1,sc_desc,&
                      eigenvec_complex(1,1,i_spin),1,i_col+1,sc_desc,&
                      (0.d0,0.d0),work,1,i_col+1,sc_desc)

                 ! Orthogonalize
                 call pzgemm('N','N',n_basis,n_states-i_col,n_block,&
                      (-1.d0,0.d0),eigenvec_complex(1,1,i_spin),1,i_done+1,&
                      sc_desc,work,1,i_col+1,sc_desc,(1.d0,0.d0),&
                      eigenvec_complex(1,1,i_spin),1,i_col+1,sc_desc)

                 i_done = i_done+n_block
              end if
           end do ! i_col
        end if ! work only on working set

        ! if needed, get the eigenvectors back to all threads
        if (collect_eigenvectors) then
           call collect_eigenvectors_scalapack_complex(&
                KS_eigenvector_complex(:,:,i_spin),i_spin)
        end if
     end do ! i_spin

     call aims_deallocate(prod,"prod")
     call aims_deallocate(ovlp_ev,"ovlp_ev")
     call aims_deallocate(work,"work")
  end if ! use_elsi

end subroutine orthonormalize_eigenvectors_scalapack_complex
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/orthonormalize_eigenvectors_scalapack_complex_GS
!  NAME
!    orthonormalize_eigenvectors_scalapack_complex
!  SYNOPSIS
subroutine orthonormalize_eigenvectors_scalapack_complex_GS()
!  PURPOSE
!    Orthonormalizes the eigenvectors with ScaLAPACK, complex version
!    This routine uses the classical Gram-Schmidt orthogonalization
!    It needs no additional memory (only the 2 1d arrays) but is slow
!    because it only uses matrix vector products
!  USES
  use localorb_io
  use mpi_tasks
  implicit none
!  ARGUMENTS
!    -
!  INPUTS
!    -
!  OUTPUT
!    -
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  complex*16, allocatable, dimension(:) :: work1, work2
  complex*16 :: dotprod
  real*8:: fact

  integer:: i_spin, i_col, info

    character*200 :: info_str


  write(info_str,'(2X,A)') "Orthonormalizing eigenvectors"
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  ! The work in this routine must be done only on the working set
  if(my_scalapack_id>=npcol*nprow) return


  allocate(work1(mxld),stat=info)
  call check_allocation(info, 'work1                         ')

  allocate(work2(mxld),stat=info)
  call check_allocation(info, 'work2                         ')


  do i_spin = 1, n_spin, 1

     do i_col = 1, n_states

        if(i_col > 1) then

           ! Orthogonalize against eigenvectors 1 .. i_col-1

           ! work1 = ovlp_complex * eigenvec_complex(i_col)

           call PZHEMV('U', n_basis, (1.0d0,0.0d0), ovlp_complex, 1, 1, sc_desc, &
                       eigenvec_complex(1,1,i_spin), 1, i_col, sc_desc, 1, &
                       (0.0d0,0.0d0), work1, 1, 1, sc_desc, 1)

           ! work2 = eigenvec_complex(1:i_col-1)**T * work1

           call PZGEMV('C', n_basis, i_col-1, (1.0d0,0.0d0), eigenvec_complex(1,1,i_spin), 1, 1, sc_desc, &
                       work1, 1, 1, sc_desc, 1, &
                       (0.0d0,0.0d0), work2, 1, 1, sc_desc, 1)

           ! eigenvec_complex(i_col) = eigenvec_complex(i_col) - eigenvec_complex(1:i_col-1) * work2

           call PZGEMV('N', n_basis, i_col-1, (-1.0d0,0.0d0), eigenvec_complex(1,1,i_spin), 1, 1, sc_desc, &
                       work2, 1, 1, sc_desc, 1, &
                       (1.0d0,0.0d0), eigenvec_complex(1,1,i_spin), 1, i_col, sc_desc, 1)
        endif

        ! Normalize eigenvec_complex(i_col)

        ! work1 = ovlp_complex * eigenvec_complex(i_col)

        call PZHEMV('U', n_basis, (1.0d0,0.0d0), ovlp_complex, 1, 1, sc_desc, &
                    eigenvec_complex(1,1,i_spin), 1, i_col, sc_desc, 1, &
                    (0.0d0,0.0d0), work1, 1, 1, sc_desc, 1)

        ! dotprod = work1**T * eigenvec_complex(i_col)

        dotprod = 0
        call PZDOTC(n_basis, dotprod, work1, 1, 1, sc_desc, 1, eigenvec_complex(1,1,i_spin), 1, i_col, sc_desc, 1)

        if(l_col(i_col) > 0) then

           ! eigenvec_complex(i_col) = eigenvec_complex(i_col)/sqrt(dotprod)

           if(dble(dotprod)>0) then
              fact = 1.d0/sqrt(dble(dotprod))
           else
              fact = 1.d0 ! for safety only, should never happen!
           endif

           eigenvec_complex(:,l_col(i_col),i_spin) = eigenvec_complex(:,l_col(i_col),i_spin)*fact

        endif

     enddo ! i_col

  enddo ! i_spin

  ! TEST only, may be removed later:
  ! call check_ev_orthogonality_complex

  deallocate(work1)
  deallocate(work2)

end subroutine orthonormalize_eigenvectors_scalapack_complex_GS
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/check_ev_orthogonality_complex
!  NAME
!    check_ev_orthogonality_complex
!  SYNOPSIS
subroutine check_ev_orthogonality_complex()
!  PURPOSE
!    Checks the orthogonality of the eigenvectors
!    This is a test routine only!!!!!!!!!!!!
!  USES
    use localorb_io
    use mpi_tasks, only: myid
    implicit none

   integer i_spin, i_col, i_row
   real*8 :: amax, dmax
   complex*16, allocatable :: work1(:,:), work2(:,:)

   ! The work in this routine must be done only on the working set
   if(my_scalapack_id>=npcol*nprow) return

   allocate(work1(mxld,mxcol))
   allocate(work2(mxld,mxcol))

   i_spin = 1

   call PZHEMM('L', 'U', n_basis, n_states, (1.0d0,0.0d0), ovlp_complex, 1, 1, sc_desc, &
               eigenvec_complex(1,1,i_spin), 1, 1, sc_desc, &
               (0.0d0,0.0d0), work1, 1, 1, sc_desc)

   call PZGEMM('C', 'N', n_states, n_states, n_basis, &
               (1.0d0,0.0d0), eigenvec_complex(1,1,i_spin), 1, 1, sc_desc, &
               work1, 1, 1, sc_desc, &
               (0.0d0,0.0d0), work2, 1, 1, sc_desc)

   dmax = 0.
   amax = 0.

   do i_col = 1, n_states
   do i_row = 1, n_states
      if(l_row(i_row) > 0 .and. l_col(i_col)>0) then
         if(i_row==i_col) then
            ! diagonal element
            dmax = max(dmax,abs(1.0-work2(l_row(i_row),l_col(i_col))))
         else
            ! off diagonal element
            amax = max(amax,abs(work2(l_row(i_row),l_col(i_col))))
         endif
      endif
   enddo
   enddo

   write(use_unit,'(A,I4,A,F10.5,F10.5)') 'check_ev_orthogonality_complex, ID: ',myid,' Errors: ',amax, dmax
   deallocate(work1)
   deallocate(work2)

end subroutine check_ev_orthogonality_complex
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_dm_scalapack
!  NAME
!    construct_dm_scalapack
!  SYNOPSIS
  subroutine construct_dm_scalapack(occ_numbers,i_spin)
!  PURPOSE
!    Construct the density matrix in ScaLAPACK
!  USES
    use aims_memory_tracking, only: aims_allocate,aims_deallocate
    use mpi_tasks
    use pbc_lists
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8, intent(in) :: occ_numbers(n_states,n_spin,n_k_points)
    integer, intent(in) :: i_spin
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8 :: factor(n_states)
    integer :: max_occ_number
    integer :: i_state
    integer :: info

    real*8, allocatable :: tmp(:,:)
    complex*16, allocatable :: tmp_complex(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id >= npcol*nprow) return

    ! If use_elsi_dm, density matrix should be computed and stored elsewhere.
    if(use_elsi_dm) return

    ! We use ham/ham_complex as storage area for the density matrix
    if(occupation_type /= 2) then ! Not Methfessel-Paxton
       factor(:) = 0.d0
       max_occ_number = 0

       do i_state = 1,n_states
          if(occ_numbers(i_state,i_spin,my_k_point) > 0.d0) then
             factor(i_state) = sqrt(occ_numbers(i_state,i_spin,my_k_point))
             max_occ_number = i_state
          end if
       end do

       if(real_eigenvectors) then
          allocate(tmp(mxld,mxcol),stat=info)
          call check_allocation(info,"tmp")

          ham(:,:,i_spin) = 0.d0
          tmp(:,:) = eigenvec(:,:,i_spin)

          do i_state = 1,n_states
             if(factor(i_state) > 0.d0) then
                if(l_col(i_state) > 0) then
                   tmp(:,l_col(i_state)) = tmp(:,l_col(i_state))*factor(i_state)
                end if
             else if(l_col(i_state) /= 0) then
                tmp(:,l_col(i_state)) = 0.d0
             end if
          end do

          call pdsyrk("U","N",n_basis,max_occ_number,1.d0,tmp,1,1,sc_desc,0.d0,&
               ham(:,:,i_spin),1,1,sc_desc)

          deallocate(tmp)
       else ! Not real eigenvectors
          allocate(tmp_complex(mxld,mxcol),stat=info)
          call check_allocation(info,"tmp_complex")

          ham_complex(:,:,i_spin) = 0.d0
          tmp_complex(:,:) = eigenvec_complex(:,:,i_spin)

          do i_state = 1,n_states
             if(factor(i_state) > 0.d0) then
                if(l_col(i_state) > 0) then
                   tmp_complex(:,l_col(i_state)) = tmp_complex(:,l_col(i_state))*factor(i_state)
                end if
             else if(l_col(i_state) /= 0) then
                tmp_complex(:,l_col(i_state)) = 0.d0
             end if
          end do

          call pzherk("U","N",n_basis,max_occ_number,(1.d0,0.d0),tmp_complex,1,&
               1,sc_desc,(0.d0,0.d0),ham_complex(:,:,i_spin),1,1,sc_desc)

          deallocate(tmp_complex)
       end if
    else ! Methfessel-Paxton
       ! dsyrk/zherk doesn't work for MP because of negative occupation numbers
       if(real_eigenvectors) then
          allocate(tmp(mxld,mxcol),stat=info)
          call check_allocation(info,"tmp")

          ham(:,:,i_spin) = 0.d0
          tmp(:,:) = eigenvec(:,:,i_spin)

          do i_state = 1,n_states
             if(l_col(i_state) > 0) then
                tmp(:,l_col(i_state)) = tmp(:,l_col(i_state))*occ_numbers(i_state,i_spin,my_k_point)
             end if
          end do

          call pdgemm("N","T",n_basis,n_basis,n_states,1.d0,tmp,1,1,sc_desc,&
               eigenvec(:,:,i_spin),1,1,sc_desc,0.d0,ham(:,:,i_spin),1,1,&
               sc_desc)

          deallocate(tmp)
       else ! Not real eigenvectors
          allocate(tmp_complex(mxld,mxcol),stat=info)
          call check_allocation(info,"tmp_complex")

          ham_complex(:,:,i_spin) = 0.d0
          tmp_complex(:,:) = eigenvec_complex(:,:,i_spin)

          do i_state = 1,n_states
             if(l_col(i_state) > 0) then
                tmp_complex(:,l_col(i_state)) = tmp_complex(:,l_col(i_state))*occ_numbers(i_state,i_spin,my_k_point)
             end if
          end do

          call pzgemm("N","C",n_basis,n_basis,n_states,(1.d0,0.d0),tmp_complex,&
               1,1,sc_desc,eigenvec_complex(:,:,i_spin),1,1,sc_desc,&
               (0.d0,0.d0),ham_complex(:,:,i_spin),1,1,sc_desc)

          deallocate(tmp_complex)
       end if
    end if

  end subroutine construct_dm_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_mulliken_decomp_scalapack
!  NAME
!    construct_mulliken_decomp_scalapack
!  SYNOPSIS
  subroutine construct_mulliken_decomp_scalapack(mulliken_decomp)
    !  PURPOSE
    !    The subroutine constructs mulliken_decomp when Scalapack is used.
    !    Please note that before calling this routine the overlap matrix must be saved
    !    (before decomposition) by calling save_overlap_scalapack
    !
    !  USES
    use dimensions, only: l_wave_max, n_atoms
    use mpi_tasks
    use pbc_lists
    use basis, only: basis_l
    implicit none
    !  ARGUMENTS
    real*8, dimension( 0:l_wave_max, n_atoms, n_states, n_spin ) :: mulliken_decomp
    !  INPUTS
    !    o none - stored overlap and density matrix are used
    !  OUTPUT
    !    o mulliken_decomp -- Mulliken matrix
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer:: i_spin, i_state, i_basis, info
    real*8:: mul_temp
    real*8, allocatable :: tmp(:,:)
    complex*16, allocatable :: tmp_complex(:,:)

    mulliken_decomp = 0.d0

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    if (real_eigenvectors) then
      allocate(tmp(mxld, mxcol),stat=info)
      call check_allocation(info, 'tmp')
    else
      allocate(tmp_complex(mxld, mxcol),stat=info)
      call check_allocation(info, 'tmp_complex')
    endif

    do i_spin = 1, n_spin

      if (real_eigenvectors) then
        call pdgemm('N','N',n_basis,n_states,n_basis,1.d0,ovlp_stored,1,1,sc_desc, &
                    eigenvec(1,1,i_spin),1,1,sc_desc,0.d0,tmp,1,1,sc_desc)
      else
        call pzgemm('N','N',n_basis,n_states,n_basis,(1.d0,0.d0),ovlp_complex_stored,1,1,sc_desc, &
                    eigenvec_complex(1,1,i_spin),1,1,sc_desc,(0.d0,0.d0),tmp_complex,1,1,sc_desc)
      endif

      do i_state = 1, n_states
        if(l_col(i_state) == 0) cycle
        do i_basis = 1, n_basis
          if(l_row(i_basis) == 0) cycle
            if(real_eigenvectors)then
              mul_temp = eigenvec(l_row(i_basis),l_col(i_state),i_spin) * &
                         tmp     (l_row(i_basis),l_col(i_state))
            else
              mul_temp = dble(conjg(eigenvec_complex(l_row(i_basis),l_col(i_state),i_spin)) * &
                                    tmp_complex     (l_row(i_basis),l_col(i_state))  )
            endif

            mulliken_decomp(basis_l(i_basis), Cbasis_to_atom(i_basis), i_state, i_spin) = &
              mulliken_decomp(basis_l(i_basis), Cbasis_to_atom(i_basis), i_state, i_spin) + mul_temp

        enddo
      enddo
    enddo

    if (real_eigenvectors) then
      deallocate(tmp)
    else
      deallocate(tmp_complex)
    endif

  end subroutine construct_mulliken_decomp_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_lowdin_decomp_scalapack
!  NAME
!    construct_lowdin_decomp_scalapack
!  SYNOPSIS
  subroutine construct_lowdin_decomp_scalapack(lowdin_decomp)
    !  PURPOSE
    !    The subroutine constructs lowdin_decomp when Scalapack is used.
    !    Please note that before calling this routine the overlap matrix must be saved
    !    (before decomposition) by calling save_overlap_scalapack
    !
    !  USES
    use dimensions, only: l_wave_max, n_atoms
    use mpi_tasks
    use pbc_lists
    use basis, only: basis_l
    use elpa2_2013
    use runtime_choices
    implicit none
    !  ARGUMENTS
    real*8, dimension( 0:l_wave_max, n_atoms, n_states, n_spin ) :: lowdin_decomp
    !  INPUTS
    !    o none - stored overlap and density matrix are used
    !  OUTPUT
    !    o lowdin_decomp -- lowdin matrix
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer:: i_spin, i_state, i_basis, info
    real*8, allocatable :: tmp_ovlp(:,:), tmp_ev(:,:), eigenvalues(:), rwork(:)
    integer, allocatable :: iwork(:)

    lowdin_decomp = 0.d0

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    if (real_eigenvectors) then
      allocate(tmp_ovlp(mxld, mxcol),stat=info)
      call check_allocation(info, 'tmp_ovlp')
      allocate(tmp_ev(mxld, mxcol),stat=info)
      call check_allocation(info, 'tmp_ev')
      allocate(eigenvalues(n_basis),stat=info)
      call check_allocation(info, 'eigenvalues')
    else
      call aims_stop('Loewdin analysis requires real eigenvectors')
    endif

    ! Get eigenvalues/eigenvectors of ovlp_stored
    ! The nonsingular eigenvalues must be the first ones, so calculate the eigenvalues of -ovlp_stored!

    tmp_ovlp = -ovlp_stored
    tmp_ev = 0.d0

    if(use_elpa) then
       call solve_evp_real_2stage_2013(n_basis, n_basis, tmp_ovlp, mxld, &
            eigenvalues, tmp_ev, mxld, nb, mpi_comm_rows, mpi_comm_cols, &
            my_scalapack_comm_work)

    else

       allocate(iwork(liwork),stat=info)
       call check_allocation(info, 'iwork')

       allocate(rwork(lrwork),stat=info)
       call check_allocation(info, 'rwork')

       call PDSYEVD('V', 'U', n_basis, tmp_ovlp, 1, 1, sc_desc, eigenvalues, &
                    tmp_ev, 1, 1, sc_desc, rwork, lrwork, iwork, liwork, info)
       if(info /= 0) call scalapack_err_exit(info,"PDSYEVD")

       deallocate(iwork)
       deallocate(rwork)

    endif

    ! We have calculated the eigenvalues of -ovlp, so invert sign of eigenvalues
    eigenvalues(1:n_basis) = -eigenvalues(1:n_basis)

    ! Get the number of nonsingular eigenvalues

    do i_basis = 1, n_basis
       if(eigenvalues(i_basis) < basis_threshold) exit
    end do
    n_nonsing_ovlp = i_basis-1

    ! Multiply columns of tmp_ev with sqrt(sqrt(eigenvalues(:))
    do i_basis = 1, n_nonsing_ovlp
      if(l_col(i_basis) > 0) tmp_ev(:,l_col(i_basis)) = tmp_ev(:,l_col(i_basis)) * sqrt(sqrt(eigenvalues(i_basis)))
    enddo

    ! Get sqrt(ovlp_stored)
    call pdgemm('N','T',n_basis,n_basis,n_nonsing_ovlp,1.d0,tmp_ev,1,1,sc_desc, &
                tmp_ev,1,1,sc_desc,0.d0,tmp_ovlp,1,1,sc_desc)


    do i_spin = 1, n_spin

      call pdgemm('N','N',n_basis,n_states,n_basis,1.d0,tmp_ovlp,1,1,sc_desc, &
                  eigenvec(1,1,i_spin),1,1,sc_desc,0.d0,tmp_ev,1,1,sc_desc)

      do i_state = 1, n_states
        if(l_col(i_state) == 0) cycle
        do i_basis = 1, n_basis
          if(l_row(i_basis) == 0) cycle

            lowdin_decomp(basis_l(i_basis), Cbasis_to_atom(i_basis), i_state, i_spin) = &
              lowdin_decomp(basis_l(i_basis), Cbasis_to_atom(i_basis), i_state, i_spin) + &
              tmp_ev(l_row(i_basis),l_col(i_state))**2

        enddo
      enddo
    enddo

    deallocate(tmp_ovlp)
    deallocate(tmp_ev)
    deallocate(eigenvalues)

  end subroutine construct_lowdin_decomp_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_dm_Pulay_forces_scalapack
!  NAME
!    construct_dm_Pulay_forces_scalapack
!  SYNOPSIS
  subroutine construct_dm_Pulay_forces_scalapack(occ_numbers, eigenvalues, evaluate_grad_psi_psi)
!  PURPOSE
!    The subroutine constructs density matrix type of matrix for Pulay force calculations
!    from scalapack type of eigenvectors.
!  USES
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
    logical :: evaluate_grad_psi_psi
!  INPUTS
!    o occ_number -- occupation numbers
!    o eigenvalues -- Kohn-Sham eigenvalues
!    o evaluate_grad_psi_psi -- are we calculating the first or second part of the Pulay forces.
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



    real*8, dimension(n_states) :: factor
    integer :: max_occ_number, i_state, i_spin, info

    real*8, allocatable :: tmp(:,:)
    complex*16, allocatable :: tmp_complex(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return


    ! we use ham/ham_complex as storage area for the density matrix




    do i_spin = 1, n_spin

      factor(:) = 0.
      max_occ_number = 0

      if(evaluate_grad_psi_psi)then

         do i_state = 1, n_states
            if(occ_numbers(i_state, i_spin, my_k_point) > 0) then
               factor(i_state) = sqrt(- 2.d0 * occ_numbers(i_state, i_spin, my_k_point) &
                    * eigenvalues(i_state, i_spin, my_k_point))
               max_occ_number = i_state
            endif
         enddo

      else

         do i_state = 1, n_states
            if(occ_numbers(i_state, i_spin, my_k_point) > 0) then
               factor(i_state) = sqrt(2*occ_numbers(i_state, i_spin, my_k_point))
               max_occ_number = i_state
            endif
         enddo
      end if




      if (real_eigenvectors) then

        allocate(tmp(mxld, mxcol),stat=info)
        call check_allocation(info, 'tmp')

        ham(:,:,i_spin) = 0.
        tmp(:,:) = eigenvec(:,:,i_spin)

        do i_state = 1, n_states
          if (factor(i_state) > 0.0d0) then
            if(l_col(i_state)>0) &
              tmp(:,l_col(i_state)) = tmp(:,l_col(i_state)) * factor(i_state)
          end if
        end do

        call pdsyrk('U', 'N', n_basis, max_occ_number, 1.0d0, tmp, 1, 1, sc_desc, &
             0.0d0, ham(1,1,i_spin), 1, 1, sc_desc )

        deallocate(tmp)

      else

        allocate(tmp_complex(mxld, mxcol),stat=info)
        call check_allocation(info, 'tmp_complex')

        ham_complex(:,:,i_spin) = 0.
        tmp_complex(:,:) = eigenvec_complex(:,:,i_spin)

        do i_state = 1, n_states
          if (factor(i_state) > 0.0d0) then
            if(l_col(i_state)>0) &
              tmp_complex(:,l_col(i_state)) = tmp_complex(:,l_col(i_state)) * factor(i_state)
          end if
        end do

        call pzherk('U', 'N', n_basis, max_occ_number, (1.0d0, 0.0d0), tmp_complex, 1, 1, sc_desc, &
             (0.0d0, 0.0d0), ham_complex(1,1,i_spin), 1, 1, sc_desc )

        deallocate(tmp_complex)

      endif

    enddo ! i_spin

  end subroutine construct_dm_Pulay_forces_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_sparse_matrix_scalapack
!  NAME
!    get_sparse_matrix_scalapack
!  SYNOPSIS
   subroutine get_sparse_matrix_scalapack( matrix_sparse, i_spin )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8 :: matrix_sparse(n_hamiltonian_matrix_size)
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = 0d0

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas2 = 1, n_basis

          lc = l_col(i_bas2) ! local column number
          if(lc==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas2) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas2), index_hamiltonian(2,i_cell,i_bas2)

             i_bas1 = column_index_hamiltonian(i_index)
             lr = l_row(i_bas1) ! local row number
             if (lr==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   ham(lr,lc,i_spin)*dble(k_phase(i_cell,my_k_point))
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble(ham_complex(lr,lc,i_spin)*dconjg(k_phase(i_cell,my_k_point)))

             endif
          end do
       end do
    end do


  end subroutine get_sparse_matrix_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_full_matrix_scalapack
!  NAME
!    get_full_matrix_scalapack
!  SYNOPSIS
   subroutine get_full_matrix_scalapack( matrix, i_spin )
!  PURPOSE
!    Gets a full matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    use dimensions, only: n_centers_basis_T, n_periodic
    use runtime_choices, only: PM_none
    implicit none
!  ARGUMENTS
    real*8 :: matrix(n_centers_basis_T, n_centers_basis_T)
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    !integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str
    integer :: i, j

    if(n_periodic> 0)then
       write(info_str, '(A)') '* ERROR: Periodic systems plus ScaLapack REQUIRE packed matrices.'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    if (packed_matrix_format /= PM_none) then
       write(info_str, '(A)') '* ERROR: get_full_matrix_scalapack works only for packed_matrix_format == PM_none'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix = 0.d0
    do i = 1, n_centers_basis_T
      do j = 1, n_centers_basis_T
        if(l_row(j)>0 .and. l_col(i)>0) matrix(j,i) = ham(l_row(j),l_col(i),i_spin)
      enddo
    enddo

  end subroutine get_full_matrix_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_forces_dm_scalapack
!  NAME
!    construct_forces_dm_scalapack
!  SYNOPSIS
  subroutine construct_forces_dm_scalapack(occ_numbers, eigenvalues, evaluate_grad_psi_psi)
!  PURPOSE
!    The subroutine constructs density matrix from scalapack type of eigenvectors.
!  USES
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
    logical:: evaluate_grad_psi_psi
!  INPUTS
!    o occ_number -- occupation numbers
!    o eigenvalues -- Kohn-Sham eigenvalues
!    o evaluate_grad_psi_psi -- are we calculating the first or second part of the Pulay forces
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    real*8, dimension(n_states) :: factor
    integer :: max_occ_number, i_state, i_spin

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return


    ! we use ham/ham_complex as storage area for the density matrix

    do i_spin = 1, n_spin

      factor(:) = 0.
      max_occ_number = 0

      do i_state = 1, n_states
        if(occ_numbers(i_state, i_spin, my_k_point) > 0) then
          if (evaluate_grad_psi_psi) then
             factor(i_state) = sqrt(abs(2*occ_numbers(i_state, i_spin, my_k_point) * &
                                        eigenvalues(i_state, i_spin, my_k_point) ) )
          else
             factor(i_state) = sqrt(2*occ_numbers(i_state, i_spin, my_k_point) )
          endif
          max_occ_number = i_state
        endif
      enddo

      if (real_eigenvectors) then

        ham(:,:,i_spin) = 0.

        do i_state = 1, n_states
          if (factor(i_state) > 0.0d0) then
            if(l_col(i_state)>0) &
              eigenvec(:,l_col(i_state),i_spin) = eigenvec(:,l_col(i_state),i_spin) * factor(i_state)
          end if
        end do

        call pdsyrk('U', 'N', n_basis, max_occ_number, 1.0d0, eigenvec(1,1,i_spin), 1, 1, sc_desc, &
             0.0d0, ham(1,1,i_spin), 1, 1, sc_desc )

        ! We need eigenvectors later in the geometry relaxation, forces and scaled zora

        do i_state = 1, n_states
          if (factor(i_state) > 0.0d0) then
            if(l_col(i_state)>0) &
              eigenvec(:,l_col(i_state),i_spin) = eigenvec(:,l_col(i_state),i_spin) / factor(i_state)
          end if
        end do

      else

        ham_complex(:,:,i_spin) = 0.

        do i_state = 1, n_states
          if (factor(i_state) > 0.0d0) then
            if(l_col(i_state)>0) &
              eigenvec_complex(:,l_col(i_state),i_spin) = eigenvec_complex(:,l_col(i_state),i_spin) * factor(i_state)
          end if
        end do

        call pzherk('U', 'N', n_basis, max_occ_number, (1.0d0, 0.0d0), eigenvec_complex(1,1,i_spin), 1, 1, sc_desc, &
             (0.0d0, 0.0d0), ham_complex(1,1,i_spin), 1, 1, sc_desc )

        ! We need eigenvectors later in the geometry relaxation, forces and scaled zora

        do i_state = 1, n_states
          if (factor(i_state) > 0.0d0) then
            if(l_col(i_state)>0) &
              eigenvec_complex(:,l_col(i_state),i_spin) = eigenvec_complex(:,l_col(i_state),i_spin) / factor(i_state)
          end if
        end do

      endif

    enddo ! i_spin

  end subroutine construct_forces_dm_scalapack
!******
!-------------------------------------------------------------------
!****s* scalapack_wrapper/evaluate_scaled_zora_tra_scalapack
!  NAME
!    evaluate_scaled_zora_tra_scalapack
!  SYNOPSIS
  subroutine evaluate_scaled_zora_tra_scalapack( &
       KS_eigenvalue, SZ_hamiltonian )
!  PURPOSE
!    The subroutine evaluates scaled zora transform using scalapack type of eigenvectors.
!  USES
    use dimensions, only: n_periodic
    use grids
    use geometry
    use basis
    use mpi_utilities
    use constants
    use synchronize_mpi_basic, only: sync_vector, sync_vector_scalapac
    use synchronize_mpi, only: sync_eigenvalues
    use localorb_io
    use runtime_choices, only: use_local_index, PM_none
    implicit none
!  ARGUMENTS
    real*8:: SZ_hamiltonian(n_hamiltonian_matrix_size, n_spin)
    real*8, dimension(n_states, n_spin, n_k_points) :: KS_eigenvalue
!  INPUTS
!    o SZ_hamiltonian -- SZ integrals.
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!  OUTPUT
!     o KS_eigenvalue -- Kohn-Sham eigenvalues after sz transformation
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    real*8:: scaled_zora_integral, scaled_zora_integral_v(n_states, n_spin)

    real*8 :: vec(n_basis)

    integer :: i_k_point, i_spin, i_state, i_compute_1, i_compute_2, info

    real*8, allocatable :: tmp(:,:)
    complex*16, allocatable :: tmp_complex(:,:)

    character*100 :: info_str

    ! Statement function for calculating index in an upper triangular compressed matrix

    integer index_b, row, col
    index_b(row,col) = MIN(row,col) + (MAX(row,col)*(MAX(row,col)-1))/2


    write(info_str,'(2X,A,A)') &
         "Evaluating scaled zora correction for each Kohn-Sham eigenstate."
    call localorb_info(info_str,use_unit,'(A)')


    if(n_periodic .eq. 0 .and. packed_matrix_format == PM_none)then

       i_k_point = 1

       do i_spin = 1, n_spin, 1

          do i_state = 1, n_states

             scaled_zora_integral = 0.d0

             ! Get eigenvector i_state on every processor

             vec(:) = 0
             do i_compute_1 = 1, n_basis
                if(l_col(i_state)>0 .and. l_row(i_compute_1)>0) then
                   vec(i_compute_1) = eigenvec(l_row(i_compute_1),l_col(i_state),i_spin)
                endif
             enddo
             call sync_vector_scalapac( vec, n_basis, my_scalapack_comm_all  )

             do i_compute_2 = 1,n_basis
                do i_compute_1 = 1,n_basis
                   scaled_zora_integral =  scaled_zora_integral + &
                          vec(i_compute_1) * vec(i_compute_2) &
                          * SZ_hamiltonian(index_b(i_compute_1,i_compute_2),i_spin)
                end do
             end do

             if(my_scalapack_id == 0)then
                KS_eigenvalue(i_state,i_spin,i_k_point) =  KS_eigenvalue(i_state,i_spin,i_k_point)/ &
                     (1+ scaled_zora_integral)
             else
                KS_eigenvalue(i_state,i_spin,i_k_point) =  0.d0
             end if

          end do
       end do

    elseif(real_eigenvectors)then

       allocate(tmp(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp')
       tmp = 0

       do i_k_point = 1, n_k_points

          if (my_k_point == i_k_point) then

             if(use_local_index) then
                call set_sparse_local_ham_scalapack( SZ_hamiltonian )
             else
                call construct_hamiltonian_scalapack( SZ_hamiltonian )
             endif

             scaled_zora_integral_v(:,:) = 0.

             do i_spin = 1, n_spin

                if(my_scalapack_id<npcol*nprow) then
                   ! tmp(:,:) = ham(:,:)*eigenvec(:,1:n_states)
                   call PDSYMM('L','U',n_basis,n_states,1.d0,ham(1,1,i_spin),1,1,sc_desc, &
                               eigenvec(1,1,i_spin),1,1,sc_desc,0.d0,tmp,1,1,sc_desc)
                endif

                do i_state = 1, n_states

                   if(l_col(i_state) == 0) cycle

                   ! scaled_zora_integral = eigenvec(:,i_state)**T * tmp(:,i_state)
                   ! contribution from local matrix parts
                   scaled_zora_integral_v(i_state,i_spin) = &
                      dot_product(eigenvec(:,l_col(i_state),i_spin),tmp(:,l_col(i_state)))
                enddo
             enddo

             ! Sum up contributions from all matrix parts
             call sync_vector(scaled_zora_integral_v,n_states*n_spin,my_scalapack_comm_all)

             if(my_scalapack_id == 0)then
                   KS_eigenvalue(:,:,i_k_point) =  KS_eigenvalue(:,:,i_k_point)/ &
                        (1+ scaled_zora_integral_v(:,:))
             else
                KS_eigenvalue(1:n_states,1:n_spin,i_k_point) =  0.d0

             end if
          else

             KS_eigenvalue(1:n_states,1:n_spin,i_k_point) =  0.d0

          end if
       end do

       deallocate(tmp)

    else ! complex eigenvectors

       allocate(tmp_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_complex')
       tmp_complex = 0

       do i_k_point = 1, n_k_points

          if (my_k_point == i_k_point) then

             if(use_local_index) then
                call set_sparse_local_ham_scalapack( SZ_hamiltonian )
             else
                call construct_hamiltonian_scalapack( SZ_hamiltonian )
             endif

             scaled_zora_integral_v(:,:) = 0.

             do i_spin = 1, n_spin

                if(my_scalapack_id<npcol*nprow) then
                   ! tmp(:,:) = ham(:,:)*eigenvec(:,1:n_states)
                   call PZHEMM('L','U',n_basis,n_states,(1.d0,0.d0),ham_complex(1,1,i_spin),1,1,sc_desc, &
                               eigenvec_complex(1,1,i_spin),1,1,sc_desc,(0.d0,0.d0),tmp_complex,1,1,sc_desc)
                endif

                do i_state = 1, n_states

                   if(l_col(i_state) == 0) cycle

                   ! scaled_zora_integral = eigenvec(:,i_state)**H * tmp(:,i_state)
                   ! contribution from local matrix parts
                   scaled_zora_integral_v(i_state,i_spin) = &
                      dble(dot_product(eigenvec_complex(:,l_col(i_state),i_spin),tmp_complex(:,l_col(i_state))))
                enddo
             enddo

             ! Sum up contributions from all matrix parts
             call sync_vector(scaled_zora_integral_v,n_states*n_spin,my_scalapack_comm_all)

             if(my_scalapack_id == 0)then
                KS_eigenvalue(:,:,i_k_point) =  KS_eigenvalue(:,:,i_k_point)/ &
                           (1+ scaled_zora_integral_v(:,:))
             else
                KS_eigenvalue(1:n_states,1:n_spin,i_k_point) =  0.d0
             end if
          else
             KS_eigenvalue(1:n_states,1:n_spin,i_k_point) =  0.d0
          end if
       end do

       deallocate(tmp_complex)

    end if


!    if (n_periodic.ne.0) then

     ! This sync call would not be necessary if we just computed KS_eigenvectors
     ! on all threads above, or propoerly parallelized the calculation.
     ! Left here as a minimal fix ... but revisit if needed.
     call sync_eigenvalues( KS_eigenvalue)

!    end if

  end subroutine evaluate_scaled_zora_tra_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/store_eigenvectors_scalapack
!  NAME
!    store_eigenvectors_scalapack
!  SYNOPSIS
      subroutine store_eigenvectors_scalapack
!  PURPOSE
!    Stores the current eigenvectors.
!  USES
        use mpi_tasks
        implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    the array eigenvec_stored/eigenvec_stored_complex is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        integer:: info

        if(real_eigenvectors)then
           if (.not.allocated(eigenvec_stored)) then
              allocate(eigenvec_stored(mxld,mxcol,n_spin),stat=info)
              call check_allocation(info, 'eigenvec_stored               ')
           end if

           eigenvec_stored = eigenvec

        else
           if (.not.allocated(eigenvec_complex_stored)) then
              allocate(eigenvec_complex_stored(mxld,mxcol,n_spin),stat=info)
              call check_allocation(info, 'eigenvec_complex_stored       ')
           end if

           eigenvec_complex_stored = eigenvec_complex

        end if


      end subroutine store_eigenvectors_scalapack
!******
!------------------------------------------------------------------------------------------
!****s* scalapack_wrapper/load_eigenvectors_scalapack
!  NAME
!    load_eigenvectors_scalapack
!  SYNOPSIS
      subroutine load_eigenvectors_scalapack
!  PURPOSE
!    Loads the eigenvectors from storage.
!  USES
        implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    the array eigenvec/eigenvec_complex is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



        if(real_eigenvectors)then

           eigenvec = eigenvec_stored

        else

           eigenvec_complex = eigenvec_complex_stored

        end if


      end subroutine load_eigenvectors_scalapack
!******
!------------------------------------------------------------------------------------------
!****s* scalapack_wrapper/remove_stored_eigenvectors_scalapack
!  NAME
!    remove_stored_eigenvectors_scalapack
!  SYNOPSIS
      subroutine remove_stored_eigenvectors_scalapack
!  PURPOSE
!    Removes the arrays for eigenvector storage.
!  USES
        implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    the array eigenvec_stored/eigenvec_complex_stored is deallocated on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


        if (allocated(eigenvec_stored)) then
           deallocate(eigenvec_stored)
        end if


        if (allocated(eigenvec_complex_stored)) then
           deallocate(eigenvec_complex_stored)
        end if


      end subroutine remove_stored_eigenvectors_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/deallocate_scalapack
!  NAME
!    deallocate_scalapack
!  SYNOPSIS
  subroutine deallocate_scalapack()
!  PURPOSE
!    Deallocates ScaLAPACK arrays.
!  USES
    use aims_memory_tracking, only: aims_deallocate
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    ScaLAPACK arrays are deallocated on exit.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE

    integer :: i_k_point

    call reinitialize_scalapack ()

    if (allocated(l_row))    deallocate(l_row)
    if (allocated(l_col))    deallocate(l_col)

    if (allocated(my_row))    deallocate(my_row)
    if (allocated(my_col))    deallocate(my_col)
    if (allocated(k_point_desc)) then
      do i_k_point = 1, size(k_point_desc)
        if (associated(k_point_desc(i_k_point)%global_id)) then
          deallocate(k_point_desc(i_k_point)%global_id)
          nullify(k_point_desc(i_k_point)%global_id)
        end if
      enddo
      deallocate(k_point_desc)
    end if

    if (allocated(ovlp))     call aims_deallocate(ovlp, "ovlp")
    if (allocated(ovlp_stored))     call aims_deallocate(ovlp_stored, "ovlp_stored")
    if (allocated(ham))      call aims_deallocate(ham, "ham")
    ! wyj add
    if (allocated(ham_stored))      call aims_deallocate(ham_stored, "ham_stored")
    if (allocated(first_order_ovlp_scalapack))     deallocate(first_order_ovlp_scalapack)
    if (allocated(first_order_ham_scalapack))      deallocate(first_order_ham_scalapack)
    if (allocated(first_order_U_scalapack))        deallocate(first_order_U_scalapack)
    if (allocated(first_order_edm_scalapack))     deallocate(first_order_edm_scalapack)
    if (allocated(momentum_matrix_scalapack))      deallocate(momentum_matrix_scalapack)
    if (allocated(Omega_MO_scalapack))      deallocate(Omega_MO_scalapack)
    if (allocated(eigenvec)) call aims_deallocate(eigenvec, "eigenvec")
    if (allocated(eigenvec_stored)) deallocate(eigenvec_stored)

    if (allocated(first_order_ham_polar_scalapack))      deallocate(first_order_ham_polar_scalapack)
    if (allocated(first_order_U_polar_scalapack))        deallocate(first_order_U_polar_scalapack)

    if (allocated(first_order_ham_polar_reduce_memory_scalapack))   deallocate(first_order_ham_polar_reduce_memory_scalapack)
    if (allocated(first_order_U_polar_reduce_memory_scalapack))   deallocate(first_order_U_polar_reduce_memory_scalapack)


    if (allocated(ovlp_complex))     call aims_deallocate(ovlp_complex, "ovlp_complex")
    if (allocated(ovlp_complex_stored))     call aims_deallocate(ovlp_complex_stored, "ovlp_complex_stored")
    if (allocated(ham_complex))      call aims_deallocate(ham_complex, "ham_complex")
    ! wyj add
    if (allocated(ham_complex_stored))      call aims_deallocate(ham_complex_stored, "ham_complex_stored")
    if (allocated(first_order_ovlp_complex_scalapack))     deallocate(first_order_ovlp_complex_scalapack)
    if (allocated(first_order_ham_complex_scalapack))      deallocate(first_order_ham_complex_scalapack)
    if (allocated(first_order_U_complex_scalapack))      deallocate(first_order_U_complex_scalapack)
    if (allocated(first_order_edm_complex_scalapack))     deallocate(first_order_edm_complex_scalapack)
    if (allocated(momentum_matrix_complex_scalapack))      deallocate(momentum_matrix_complex_scalapack)
    if (allocated(Omega_MO_complex_scalapack))      deallocate(Omega_MO_complex_scalapack)
    if (allocated(eigenvec_complex)) call aims_deallocate(eigenvec_complex, "eigenvec_complex")
    if (allocated(eigenvec_complex_stored)) deallocate(eigenvec_complex_stored)


    if (allocated(gl_prow)) deallocate(gl_prow)
    if( allocated(gl_pcol)) deallocate(gl_pcol)
    if( allocated(gl_pcol)) deallocate(gl_pcol)
    if(allocated(gl_nprow)) deallocate(gl_nprow)
    if(allocated(gl_npcol)) deallocate(gl_npcol)


  end subroutine deallocate_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/finalize_scalapack
!  NAME
!    finalize_scalapack
!  SYNOPSIS
  subroutine finalize_scalapack
!  PURPOSE
!    Removes the BLACS grid.
!  USES
  use runtime_choices
  implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    if(my_scalapack_id<npcol*nprow) then
      call BLACS_Gridexit(my_blacs_ctxt)

!AJL/Feb2014: Debug
!      write(use_unit,*) 'Exiting BLACS. Processor: ', my_scalapack_id
!
! For QM/MM, where repeat calls are made to FHI-aims, the use of BLACS_Exit
! is a problem as it seems the Scalapack environment cannot be reinitialised
! once finalised.
!
! E.g. if we force KS_method lapack, irrespective of whether Scalapack is compiled in,
! then everything works OK because Scalapack is neither initialised or finalised.
!
! But if we use scalapack then the environment is initialised, finalised and then falls
! over when you try to reinitialise as it thinks no processors are available to work
! with. Any suggestions welcomed, else we have to hack this for QM/MM...
!
! ...which is not ideal. I have not edited anything here - just written my comments.
!
! Update: AJL Feb/2016
! After discussion with Volker is has been concluded that, if we cannot repeatedly make
! a call to BLACS_Exit, then this subroutine call should be moved *outside* of the
! main aims() call and any calling software will need to pay attention and close
! the BLACS environment down in an approriate manner.
!
! In practical terms for FHI-aim users, this subroutine call has been moved to aims.f90
!
!      call BLACS_Exit(1)
!
    endif

! It seems this was not called previously.
    call deallocate_scalapack
!AJL/Feb2014

  end subroutine finalize_scalapack
!******
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! Obsolete routines
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!  subroutine print_ominais_vektori
!
!    if(myid==0) write(use_unit,*) 'EEEp', eigenvec_complex(1,1,1)
!
!  end subroutine print_ominais_vektori
!
!-----------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !****s* scalapack_wrapper/collect_eigenvectors_scalapack
  !  NAME
  !    collect_eigenvectors_scalapack
  !  SYNOPSIS

  subroutine collect_eigenvectors_scalapack(KS_eigenvector, i_spin)

    !  PURPOSE
    !    Copy/collect the contents of eigenvec into KS_eigenvector.
    !  USES

    use synchronize_mpi_basic, only: sync_vector
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_spin
    real*8, intent(OUT) :: KS_eigenvector(n_basis, n_states)

    !  INPUTS
    !    o i_spin -- spin channel to collect
    !  OUTPUTS
    !    o KS_eigenvector -- updated one-particle eigencoefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_col, i_row

    KS_eigenvector = 0.d0
    do i_col = 1, n_states
       if(l_col(i_col)==0) cycle
       do i_row = 1, n_basis
          if(l_row(i_row)>0) then
             KS_eigenvector(i_row, i_col) &
             & = eigenvec(l_row(i_row), l_col(i_col), i_spin)
          end if
       end do
    end do
    call sync_vector(KS_eigenvector, n_basis*n_states, my_scalapack_comm_all)


  end subroutine collect_eigenvectors_scalapack
  !******
  !------------------------------------------------------------------------------
  !****s* scalapack_wrapper/collect_eigenvectors_scalapack_complex
  !  NAME
  !    collect_eigenvectors_scalapack_complex
  !  SYNOPSIS

  subroutine collect_eigenvectors_scalapack_complex(KS_eigenvector_complex, i_spin)

    !  PURPOSE
    !    Copy/collect the contents of eigenvec_complex into
    !    KS_eigenvector_complex.
    !  USES

    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_spin
    complex*16, intent(OUT) :: KS_eigenvector_complex(n_basis, n_states)

    !  INPUTS
    !    o i_spin -- spin channel to collect
    !  OUTPUTS
    !    o KS_eigenvector_complex -- updated one-particle eigencoefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_col, i_row

    KS_eigenvector_complex = (0.d0, 0.d0)
    do i_col = 1, n_states
       if(l_col(i_col)==0) cycle
       do i_row = 1, n_basis
          if(l_row(i_row)>0) then
             KS_eigenvector_complex(i_row, i_col) &
             & = eigenvec_complex(l_row(i_row), l_col(i_col), i_spin)
          end if
       end do
    end do
    call sync_vector_complex(KS_eigenvector_complex, n_basis*n_states, &
    &                        my_scalapack_comm_all)

  end subroutine collect_eigenvectors_scalapack_complex
  !******
  !------------------------------------------------------------------------------
  !****s* scalapack_wrapper/spread_eigenvectors_scalapack
  !  NAME
  !    spread_eigenvectors_scalapack
  !  SYNOPSIS

  subroutine collect_generic_eigenvectors_scalapack(matrix_eigenvec, n_rank, eigenvector_out, i_spin)

    !  PURPOSE
    !    Copy/collect the contents of distributed matrix_eigenvec (e.g. after solve_evp_real) into eigenvector_out.
    !  USES

    use synchronize_mpi_basic, only: sync_vector
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_spin
    integer, intent(IN) :: n_rank
    real*8, intent(IN) :: matrix_eigenvec(mxld, mxcol, n_spin)
    real*8, intent(OUT) :: eigenvector_out(n_rank, n_rank)

    !  INPUTS
    !    o i_spin -- spin channel to collect
    !    o n_rank -- rank of the matrix
    !    o matrix_eigenvec -- the distributed scalapack array to collect
    !  OUTPUTS
    !    o eigenvector_out -- collected eigenvector from matrix_eigenvec
    !  AUTHOR
    !    Christoph Schober
    !  HISTORY
    !    First version (Sept. 2014)
    !  SOURCE

    integer :: i_col, i_row

    eigenvector_out = 0.d0
    do i_col = 1, n_rank
       if(l_col(i_col)==0) cycle
       do i_row = 1, n_rank
          if(l_row(i_row)>0) then
             eigenvector_out(i_row, i_col) &
             & = matrix_eigenvec(l_row(i_row), l_col(i_col), i_spin)
          end if
       end do
    end do
    call sync_vector(eigenvector_out, n_rank*n_rank, my_scalapack_comm_all)


  end subroutine collect_generic_eigenvectors_scalapack
  subroutine spread_eigenvectors_scalapack(KS_eigenvector, KS_eigenvector_complex)

    !  PURPOSE
    !    Copy the contents of KS_eigenvector{_complex} into
    !    eigenvec{_complex}.
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) ::KS_eigenvector(n_basis, n_states, n_spin, 1)
    complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, &
    &                                                n_spin, 1)

    !  INPUTS
    !    o KS_eigenvector{_complex} -- one-particle eigencoefficients to be spread
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_spin, i_col, i_row

    do i_spin = 1, n_spin
       do i_col = 1, n_states
          if(l_col(i_col)==0) cycle
          do i_row = 1, n_basis
             if(l_row(i_row)>0) then
                if (real_eigenvectors) then
                   eigenvec(l_row(i_row), l_col(i_col), i_spin) &
                   & = KS_eigenvector(i_row, i_col, i_spin, 1)
                else
                   eigenvec_complex(l_row(i_row), l_col(i_col), i_spin) &
                   & = KS_eigenvector_complex(i_row, i_col, i_spin, 1)
                end if
             end if
          end do
       end do
    end do

  end subroutine spread_eigenvectors_scalapack
  !******
  !------------------------------------------------------------------------------
  !****s* scalapack_wrapper/sync_single_eigenvec_scalapack
  !  NAME
  !    sync_single_eigenvec_scalapack
  !  SYNOPSIS

  subroutine sync_single_eigenvec_scalapack(KS_vec, i_state, i_spin, i_k_point)

    !  PURPOSE
    !    Sync single state into KS_vec to /all/ procs.
    !  USES

    use synchronize_mpi_basic, only: sync_vector
    implicit none

    !  ARGUMENTS

    real*8, intent(OUT) :: KS_vec(n_basis)
    integer, intent(IN) :: i_state
    integer, intent(IN) :: i_spin
    integer, intent(IN) :: i_k_point

    !  INPUTS
    !    o i_state -- state number to collect
    !    o i_spin -- spin channel to collect
    !    o i_k_point -- k-point number to collect
    !  OUTPUTS
    !    o KS_eigenvector -- updated one-particle eigencoefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_col, i_row

    KS_vec = 0.d0

    if (i_k_point == my_k_point) then
       i_col = i_state
       if (l_col(i_col) /= 0) then
          do i_row = 1, n_basis
             if(l_row(i_row)>0) then
                KS_vec(i_row) = eigenvec(l_row(i_row), l_col(i_col), i_spin)
             end if
          end do
       end if
    end if
    call sync_vector(KS_vec, n_basis)

  end subroutine sync_single_eigenvec_scalapack
  !******
  !------------------------------------------------------------------------------
  !****s* scalapack_wrapper/sync_single_eigenvec_scalapack_complex
  !  NAME
  !    sync_single_eigenvec_scalapack_complex
  !  SYNOPSIS

  subroutine sync_single_eigenvec_scalapack_complex(KS_vec_complex, &
  &                                                 i_state, i_spin, i_k_point)

    !  PURPOSE
    !    Sync single state into KS_vec_complex to /all/ procs.
    !  USES

    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none

    !  ARGUMENTS

    complex*16, intent(OUT) :: KS_vec_complex(n_basis)
    integer, intent(IN) :: i_state
    integer, intent(IN) :: i_spin
    integer, intent(IN) :: i_k_point

    !  INPUTS
    !    o i_state -- state number to collect
    !    o i_spin -- spin channel to collect
    !    o i_k_point -- k-point number to collect
    !  OUTPUTS
    !    o KS_eigenvector -- updated one-particle eigencoefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_col, i_row

    KS_vec_complex = (0.d0, 0.d0)

    if (i_k_point == my_k_point) then
       i_col = i_state
       if (l_col(i_col) /= 0) then
          do i_row = 1, n_basis
             if(l_row(i_row)>0) then
                KS_vec_complex(i_row) &
                & = eigenvec_complex(l_row(i_row), l_col(i_col), i_spin)
             end if
          end do
       end if
    end if
    call sync_vector_complex(KS_vec_complex, n_basis)

  end subroutine sync_single_eigenvec_scalapack_complex
  !******

!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/restart_scalapack_read
!  NAME
!    restart_scalapack_read
!  SYNOPSIS
  subroutine restart_scalapack_read( density_matrix_sparse, i_spin, is_found)
!  PURPOSE
!    Reads  density matrix from external storage.
!  USES
    use localorb_io
    use mpi_tasks, only: myid
    use pbc_lists
    use runtime_choices
    use synchronize_mpi, only: bcast_logical, sync_density_matrix_sparse
    implicit none
!  ARGUMENTS
    real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
    integer::  i_spin
    logical, intent(OUT) :: is_found
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o density_matrix_sparse -- density matrix
!    o is_found -- whether restart file has been found and read
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    character*150 :: info_str, name_temp
    logical:: exists
    integer :: i_index

    last_restart_saving = 0
    is_found = .false.

    if(.not. read_restart(i_spin)) return


    if(n_spin == 1)then
       write(name_temp,'(A)') trim(restart_read_file)
    else
       write(info_str,'(A,A,I1)') trim(restart_read_file),'.',i_spin
       write(name_temp,'(A,A,I1)')trim(info_str)
    end if


    if(myid == 0)then
       inquire(FILE=name_temp,EXIST=exists)
    end if
    call bcast_logical(exists, 0)

    if(.not. exists)then
       read_restart(i_spin) = .false.
       return
    end  if

    write(info_str,'(2X,2A)') 'Reading scalapack restart information from file ', trim(name_temp)
    call localorb_info(info_str,use_unit,'(A)')

    if(myid==0)then
       open(file = name_temp,  unit = 7, status = 'old', form = 'unformatted')

       do i_index = 1, n_hamiltonian_matrix_size

          read(7)  density_matrix_sparse(i_index)

       end do

       close(7)
    else

       density_matrix_sparse = 0.d0

    end if

    call sync_density_matrix_sparse(density_matrix_sparse)
    read_restart(i_spin) = .false.
    is_found = .true.

!    write(use_unit,*) 'DM R', density_matrix_sparse(1),  density_matrix_sparse(2)

  end subroutine restart_scalapack_read
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/restart_scalapack_write
!  NAME
!    restart_scalapack_write
!  SYNOPSIS
  subroutine restart_scalapack_write( density_matrix_sparse,i_spin)
!  PURPOSE
!    Writes density matrix from external storage.
!  USES
    use localorb_io
    use mpi_tasks, only: myid
    use pbc_lists
    use runtime_choices
    implicit none
!  ARGUMENTS
    real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
    integer:: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!    o density_matrix_sparse -- density matrix
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    character*150 :: info_str, name_temp
    integer :: i_index





    if(last_restart_saving < restart_save_iterations) then

       if(i_spin ==n_spin ) last_restart_saving = last_restart_saving + 1

    else

       if(i_spin ==n_spin) last_restart_saving = 0


       if(n_spin == 1)then
          write(name_temp,'(A)') trim(restart_write_file)
       else
          write(info_str,'(A,A,I1)') trim(restart_write_file),'.',i_spin
          write(name_temp,'(A)')trim(info_str)
       end if


       write(info_str,'(2X,2A)') 'Writing scalapack restart information to file ', trim(name_temp)
       call localorb_info(info_str,use_unit,'(A)',OL_norm)

       if(myid == 0)then
          open(file = name_temp, unit = 88, status = 'replace', form = 'unformatted',action='write')

          do i_index = 1, n_hamiltonian_matrix_size

             write(88)  density_matrix_sparse(i_index)

          end do

          close(88)
       end if
    end if

!    write(use_unit,*) 'DM W', density_matrix_sparse(1),   density_matrix_sparse(2)

  end subroutine restart_scalapack_write
!******
!---------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! Transport routines
!-------------------------------------------------------------------------------
!******
!-------------------------------------------------------------------------
  !****s* scalapack_wrapper/construct_hamiltonian_scalapack
  !  NAME
  !    construct_hamiltonian_scalapack
  !  SYNOPSIS
  subroutine construct_ham_and_ovl_transport_scalapack( hamiltonian, overlap_matrix )
    !  PURPOSE
    !    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage.
    !  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io
    use runtime_choices
    implicit none
    !  ARGUMENTS
    real*8:: hamiltonian   ( n_hamiltonian_matrix_size, n_spin )
    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
    !  INPUTS
    !    o hamiltonian -- the Hamilton matrix
    !  OUTPUT
    !    upper half of the ScaLAPACK array ham is set on exit
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    if(myid==0) write(use_unit,*) 'Constructing Hamiltonian and Overlap matrix scalapack'
    if(myid==0) write(use_unit,*) ' format for transport calculation'

    ! If a local index is used, call special routine

    ! Ppaula


    if(use_local_index) then

       write(use_unit,*) 'use local index not implemented'
       stop
!       do i_spin = 1,n_spin
!          call get_set_sparse_local_matrix_scalapack(hamiltonian(1,i_spin),1,i_spin)
!       enddo
!       call get_set_sparse_local_matrix_scalapack(overlap_matrix,0,1)
!       return
    endif

    if(n_k_points > 1)then
       call  deallocate_scalapack
      ! call  finalize_scalapack
       n_k_points = 1
       call initialize_scalapack
    end if


    ! kohta
    if(real_eigenvectors)then
       ham(:,:,:) = 0.
       ovlp = 0.d0
    else
       ham_complex(:,:,:) = 0.
       ovlp_complex = 0.0
    end if


    select case(packed_matrix_format)

    case(PM_index) !------------------------------------------------

       do i_cell = 1, n_cells_in_hamiltonian-1

          if(cell_index(i_cell,1) == 0 .and. cell_index(i_cell,2) == 0 .and. cell_index(i_cell,3) == 0)then

             do i_col = 1, n_basis

                lc = l_col(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_row(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         do i_spin = 1, n_spin
                            ham (lr,lc,i_spin) = ham (lr,lc,i_spin) &
                                 +  hamiltonian(idx,i_spin)
                         end do
                         ovlp(lr,lc) =    ovlp(lr,lc) + overlap_matrix(idx)


                      else ! complex eigenvectors
                         do i_spin = 1, n_spin
                            ham_complex (lr,lc,i_spin) = ham_complex (lr,lc,i_spin) &
                                 + hamiltonian(idx,i_spin)
                         end do

                         ovlp_complex(lr,lc) =  ovlp_complex(lr,lc) + overlap_matrix(idx)
                      end if ! real_eigenvectors

                   end do

                end if
             end do
          end if
       end do ! i_cell


       do i_cell = 1, n_cells_in_hamiltonian-1

          if(cell_index(i_cell,1) == 0 .and. cell_index(i_cell,2) == 0 .and. cell_index(i_cell,3) == 0)then

             do i_col = 1, n_basis

                lc = l_row(i_col) ! local column number
                if(lc==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_col) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)

                      i_row = column_index_hamiltonian(idx)
                      lr = l_col(i_row) ! local row number
                      if(lr==0) cycle   ! skip if not local

                      if(i_col /= i_row)then

                         if(real_eigenvectors)then
                            do i_spin = 1, n_spin
                               ham (lc,lr,i_spin) = ham (lc,lr,i_spin) &
                                    +  hamiltonian(idx,i_spin)
                            end do
                            ovlp(lc,lr) =    ovlp(lc,lr) + overlap_matrix(idx)


                         else ! complex eigenvectors
                            do i_spin = 1, n_spin
                               ham_complex (lc,lr,i_spin) = ham_complex (lc,lr,i_spin) &
                                    + hamiltonian(idx,i_spin)
                            end do

                            ovlp_complex(lc,lr) =  ovlp_complex(lc,lr) + overlap_matrix(idx)
                         end if ! real_eigenvectors
                      end if

                   end do

                end if
             end do
          end if
       end do ! i_cell




    case(PM_none) !---------------------------------------------------------

       write(use_unit,*) 'Error: construct_hamiltonian_scalapack does not support non-packed matrices.'
       call aims_stop

    end select ! packed_matrix_format






  end subroutine construct_ham_and_ovl_transport_scalapack

!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_greenfunction_scalapack
!  NAME
!    construct_overlap_scalapack
!  SYNOPSIS
  subroutine construct_greenfunction_scalapack( energy )
!  PURPOSE
!    Sets the overlap matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use pbc_lists
    use geometry
    use basis
    use localorb_io
    use mpi_tasks
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    real*8:: energy
!  INPUTS
!    o overlap_matrix -- the overlap matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_cell


    IF(.not.allocated(green_work))then
       allocate(green_work(mxld,mxcol),stat=i_cell)
       call check_allocation(i_cell, 'green_work                   ')
    end IF

    if(use_local_index) then
       write(use_unit,*) 'Error use_local_index not implemented to transport!'
       stop
    endif

!    green_work = 0.d0

    if(real_eigenvectors)then
       green_work =  energy * ovlp - ham(:,:,1)
    else
       green_work =  energy * ovlp_complex - ham_complex(:,:,1)
    end if



  end subroutine construct_greenfunction_scalapack
!******
!******
!-----------------------------------------------------------------------------------
  subroutine solve_greens_functions( greenL, &  !green_13, green_14, green_34, &
        i_L1_1,  i_L1_2, i_L2_1, i_L2_2, i_L3_1, i_L3_2, i_L4_1, i_L4_2, &
        basis_L1, basis_L2, basis_L3, basis_L4, n_lead_basis, n_leads , work,n_task_energy_steps,  i_task,&
        receive_12,  receive_13,  receive_14,   receive_34 )

!lineaari

    use mpi_tasks
    use synchronize_mpi_basic
    implicit none

    integer:: info, i_task, n_task_energy_steps
    integer:: receive_12(n_task_energy_steps), receive_13(n_task_energy_steps), &
              receive_14(n_task_energy_steps), receive_34(n_task_energy_steps)
    integer:: jj, ii,s, sy, L, dim
    integer::  basis_L1, basis_L2, basis_L3, basis_L4, n_lead_basis
    integer::  i_L1_1,  i_L1_2, i_L2_1, i_L2_2, i_L3_1, i_L3_2, i_L4_1, i_L4_2, n_leads, lr, lc
    complex*16, dimension(n_lead_basis,n_lead_basis):: work
    complex*16, dimension(basis_L1,basis_L1)::  greenL
!    complex*16, dimension(:,:):: green_13, green_14, green_34


    IF(.not.allocated(green))then
       allocate(green(mxld,mxcol),stat=info)
       call check_allocation(info, 'green                        ')
       allocate(green_ipiv(nb + n_basis),stat=info)
       call check_allocation(info, 'green_ipiv                   ')

    end IF


    green = 0.d0
    s = 0

    do L = i_L1_1, i_L1_2
       s = s+1

       lr = l_row(L)     ! local row number
       if(lr==0) cycle   ! skip if not local

       lc = l_col(s)     ! local row number
       if(lc==0) cycle   ! skip if not local

       green(lr,lc) = 1.d0
    end do

    if(n_leads > 2)then
       do L = i_L3_1, i_L3_2
          s = s+1

          lr = l_row(L)     ! local row number
          if(lr==0) cycle   ! skip if not local

          lc = l_col(s)     ! local row number
          if(lc==0) cycle   ! skip if not local

          green(lr,lc) = 1.d0
       end do
    end if

    dim = s

    call PZGETRF(  n_basis, n_basis, green_work, 1, 1, sc_desc, green_ipiv, info )

    call PZGETRS( 'N',n_basis, dim, green_work, 1, 1 , sc_desc, green_ipiv, green, &
         1, 1, sc_desc, info )


    ! Tunneling lead 1 -> lead 2
    s = 0
    work = 0.d0
    do ii = i_L1_1,  i_L1_2
       s = s+1

       lc = l_col(s)    ! local column number
       if(lc==0) cycle   ! skip if not local

       sy = 0
       do jj =  i_L2_1,  i_L2_2
          sy = sy +1

          lr = l_row(jj)    ! local row number
          if(lr==0) cycle   ! skip if not local

          work(jj-i_L2_1+1,ii-i_L1_1+1) = green(lr,lc)
       end do
    end do

    !    write(use_unit,*) 'ulos'
    !    write(use_unit,*) work(:,1)
    !    stop

!    do ii =  1,n_task_energy_steps
!       write(use_unit,*) ii, 'col', receive_12(i_task)
       call collect_vector_complex(work, greenL,  n_lead_basis**2, receive_12(i_task))
!    end do
!       write(use_unit,*) 'ohi'

!    call  sync_vector_complex(work, n_lead_basis**2 )

!    if(myid == i_task)then
!       green_12 = work
!    end if

    if(n_leads >= 3)then

       ! Tunneling lead 1 -> lead 3
       s = 0
       work = 0.d0
       do ii = i_L1_1,  i_L1_2
          s = s+1

          lc = l_col(s)    ! local column number
          if(lc==0) cycle   ! skip if not local

          sy = 0
          do jj =  i_L3_1,  i_L3_2
             sy = sy +1

             lr = l_row(jj)    ! local row number
             if(lr==0) cycle   ! skip if not local

             work(jj-i_L3_1+1,ii-i_L1_1+1) = green(lr,lc)
          end do
       end do


!    do ii =  1,n_task_energy_steps
 !      write(use_unit,*) ii, 'col', receive_12(ii)
       call collect_vector_complex(work, greenL,  n_lead_basis**2, receive_13(i_task))
!    end do

!       call collect_vector_complex(work, green_13,  n_lead_basis**2, i_task)
!       call  sync_vector_complex(work, n_lead_basis**2 )

!       if(myid == i_task)then
!          green_13 = work
!       end if
    end if


    if(n_leads >= 4)then

       ! Tunneling lead 1 -> lead 4
       s = 0
       work = 0.d0
       do ii = i_L1_1,  i_L1_2
          s = s+1

          lc = l_col(s)    ! local column number
          if(lc==0) cycle   ! skip if not local

          sy = 0
          do jj =  i_L4_1,  i_L4_2
             sy = sy +1

             lr = l_row(jj)    ! local row number
             if(lr==0) cycle   ! skip if not local

             work(jj-i_L4_1+1,ii-i_L1_1+1) = green(lr,lc)
          end do
       end do

!    do ii =  1,n_task_energy_steps
        call collect_vector_complex(work, greenL,  n_lead_basis**2, receive_14(i_task))
!    end do

!       call collect_vector_complex(work, green_14,  n_lead_basis**2, i_task)
!       call  sync_vector_complex(work, n_lead_basis**2 )

!       if(myid == i_task)then
!          green_14 = work
!       end if


       ! Tunneling lead 3 -> lead 4
       work = 0.d0

       do ii = i_L3_1,  i_L3_2
          s = s+1

          lc = l_col(s)    ! local column number
          if(lc==0) cycle   ! skip if not local

          sy = 0
          do jj =  i_L4_1,  i_L4_2
             sy = sy +1

             lr = l_row(jj)    ! local row number
             if(lr==0) cycle   ! skip if not local

             work(jj-i_L4_1+1,ii-i_L3_1+1) = green(lr,lc)

          end do
       end do

!    do ii =  1,n_task_energy_steps
        call collect_vector_complex(work, greenL,  n_lead_basis**2, receive_34(i_task))
!    end do

!       call collect_vector_complex(work, green_34,  n_lead_basis**2, i_task)
!       call  sync_vector_complex(work, n_lead_basis**2 )

!       if(myid == i_task)then
!          green_34 = work
!       end if

    end if

  end subroutine solve_greens_functions


!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/add_self_energy_to_greenfunction_scalapack
!  NAME
!    add_self_energy_to_greenfunction_scalapack
!  SYNOPSIS
  subroutine add_self_energy_to_greenfunction_scalapack( self_energy, i1, i2, dim )
!  PURPOSE
!    Sets the overlap matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use pbc_lists
    use geometry
    use basis

    implicit none
!  ARGUMENTS
    integer::  i1, i2, dim
    complex*16:: self_energy(dim,dim)


!  INPUTS

!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!Paula

    integer:: lr, lc, ii, jj


    do ii = i1, i2, 1

       lc = l_col(ii)    ! local column number
       if(lc==0) cycle   ! skip if not local

       do jj = i1, i2, 1

          lr = l_row(jj)    ! local row number
          if(lr==0) cycle   ! skip if not local


          green_work(lr,lc) = green_work(lr,lc) - self_energy(ii-i1+1,jj-i1+1)

       end do
    end do





  end subroutine add_self_energy_to_greenfunction_scalapack
!******
!******
!-----------------------------------------------------------------------------------
  subroutine transport_proj_weight(  min_eigenvalue, min_states, overlap_matrix, KS_eigenvalue )

    use dimensions, only: n_atoms
    use pbc_lists
    use synchronize_mpi_basic, only: sync_vector, sync_vector_complex
    use runtime_choices, only: use_local_index
    implicit none

    real*8:: projected_weight_r(n_atoms),  min_eigenvalue(n_atoms,n_k_points)
    real*8:: min_states(n_atoms,n_k_points)
    integer :: i_state, i_spin, i_k_point, i_atom, i_basis_1, i_basis_2, i_row, i_cell, i_index, i_size

    real*8 :: KS_eigenvec(n_basis), mul_temp
    complex*16:: KS_eigenvec_complex(n_basis)
    real*8:: overlap_matrix( n_hamiltonian_matrix_size )
    real*8,     dimension(n_states, n_spin,n_k_points) :: KS_eigenvalue


!    write(use_unit,*) 'transport_find_local_min_occupated_state scalapack'






    i_k_point = my_k_point

    do i_spin = 1, n_spin
       do i_state = 1, n_states


          if(real_eigenvectors)then
             KS_eigenvec(:) = 0
             if(l_col(i_state)>0)then
                do i_row = 1, n_basis
                   if(l_row(i_row)>0) then
                      KS_eigenvec(i_row) = eigenvec(l_row(i_row),l_col(i_state),i_spin)
                   endif
                end do
             end if
             call sync_vector( KS_eigenvec, n_basis, my_scalapack_comm_all )

             projected_weight_r = 0.d0



             do i_cell = 1,n_cells_in_hamiltonian-1

                do i_basis_2 = 1, n_basis

                   if((my_scalapack_id == mod(i_basis_2,my_scalapack_task_size )) .or. use_local_index )then


                      if( index_hamiltonian(1,i_cell, i_basis_2) > 0 ) then

                         i_index = index_hamiltonian(1,i_cell, i_basis_2)-1

                         do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                            i_index = i_index + 1
                            i_basis_1 =  column_index_hamiltonian(i_index)

                            mul_temp = overlap_matrix(i_index) * KS_eigenvec(i_basis_1) &
                                 * KS_eigenvec(i_basis_2) * dble(k_phase(i_cell,i_k_point))

                            ! 1st pass over all matrix elements
                            projected_weight_r(Cbasis_to_atom(i_basis_1)) = projected_weight_r(Cbasis_to_atom(i_basis_1)) &
                                 + mul_temp

                            ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                            if (i_basis_1.ne.i_basis_2) then

                               projected_weight_r(Cbasis_to_atom(i_basis_2)) = projected_weight_r(Cbasis_to_atom(i_basis_2)) &
                                    + mul_temp

                            end if

                         end do
                      end if
                   end if
                end do
             end do ! i_cell

          else  ! Complex eigenvectors


             KS_eigenvec_complex(:) = 0
             if(l_col(i_state)>0)then
                do i_row = 1, n_basis
                   if(l_row(i_row)>0) then
                      KS_eigenvec_complex(i_row) = eigenvec_complex(l_row(i_row),l_col(i_state),i_spin)
                   endif
                end do
             end if
             call sync_vector_complex( KS_eigenvec_complex, n_basis, my_scalapack_comm_all )


             projected_weight_r = 0.d0



             do i_cell = 1,n_cells_in_hamiltonian-1

                do i_basis_2 = 1, n_basis

                   if((my_scalapack_id == mod(i_basis_2,my_scalapack_task_size )) .or. use_local_index )then


                      if( index_hamiltonian(1,i_cell, i_basis_2) > 0 ) then

                         i_index = index_hamiltonian(1,i_cell, i_basis_2)-1

                         do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)

                            i_index = i_index + 1
                            i_basis_1 =  column_index_hamiltonian(i_index)


!                            mul_temp = dble(overlap_matrix(i_index) * KS_eigenvec_complex(Cbasis_to_basis(i_basis_1)) &
!                                 *dconjg( KS_eigenvec_complex(Cbasis_to_basis(i_basis_2))) * dconjg(k_phase(i_cell,i_k_point)))

                            mul_temp = dble(overlap_matrix(i_index) * KS_eigenvec_complex(i_basis_1) &
                                 *dconjg( KS_eigenvec_complex(i_basis_2)) * dconjg(k_phase(i_cell,i_k_point)))

                            ! 1st pass over all matrix elements

                            projected_weight_r(Cbasis_to_atom(i_basis_1)) = projected_weight_r(Cbasis_to_atom(i_basis_1)) &
                                 + mul_temp

                            ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                            if (i_basis_1.ne.i_basis_2) then

                               mul_temp = dble(overlap_matrix(i_index) * dconjg( KS_eigenvec_complex(i_basis_1)) &
                                    * ( KS_eigenvec_complex(i_basis_2)) * (k_phase(i_cell,i_k_point)))

                               projected_weight_r(Cbasis_to_atom(i_basis_2)) = projected_weight_r(Cbasis_to_atom(i_basis_2)) &
                                    + mul_temp

                            end if
                         end do
                      end if
                   end if
                end do ! i_basis_2
             end do ! i_cell
          end if ! real_eigenvectors


          call sync_vector( projected_weight_r,  n_atoms , my_scalapack_comm_all )


          if(my_scalapack_id==0)then


             projected_weight_r = abs(projected_weight_r)*k_weights(i_k_point)

             do i_atom = 1, n_atoms

                if( projected_weight_r(i_atom)  > 1e-4)then


                   if( min_states(i_atom,i_k_point) < 5e-5) then

                      min_states(i_atom,i_k_point) = projected_weight_r(i_atom)
                      min_eigenvalue(i_atom,i_k_point) = KS_eigenvalue(i_state,i_spin,i_k_point)*projected_weight_r(i_atom)

                   else if( abs(min_eigenvalue(i_atom,i_k_point)/ min_states(i_atom,i_k_point) &
                        -  KS_eigenvalue(i_state,i_spin,i_k_point)) < 2)then

                      min_states(i_atom,i_k_point) =  min_states(i_atom,i_k_point) + projected_weight_r(i_atom)
                      min_eigenvalue(i_atom,i_k_point) =  min_eigenvalue(i_atom,i_k_point) &
                           +  KS_eigenvalue(i_state,i_spin,i_k_point)*projected_weight_r(i_atom)

                   else if( min_eigenvalue(i_atom,i_k_point)/ min_states(i_atom,i_k_point) &
                        >  KS_eigenvalue(i_state,i_spin,i_k_point)) then
                      min_states(i_atom,i_k_point) = projected_weight_r(i_atom)
                      min_eigenvalue(i_atom,i_k_point) = KS_eigenvalue(i_state,i_spin,i_k_point)*projected_weight_r(i_atom)
                   end if

                end if

             end do
          end if
       end do ! i_state
    end do ! i_spin


  end subroutine transport_proj_weight

!-----------------------------------------------------------------------------------

! End transport section of scalapack wrapper as provided by Paula



!=========================================================================================
!=========================begin for scalapack used in DFPT_phonon=========================
!=========================================================================================

!----------------lapack version------------  |---------scalapack version-------------------
!(0)                                         |  initialize_scalapack_DFPT_phonon
!(1)   trans_first_order_sparse_to_matrix    |  construct_first_order_overlap_supercell_scalapack
!(2)   trans_first_order_sparse_to_matrix    |  construct_first_order_hamiltonian_supercell_scalapack
!(3)   evaluate_first_order_DM_supercell_p1  |  construct_first_order_dm_supercell_scalapack
!(4)   add_matrix_to_sparse                  |  get_first_order_dm_sparse_matrix_from_supercell_scalapack
!(5)   evaluate_first_order_DM_supercell_p1  |  evaluate_first_order_U_supercell_scalapack
!(6)   evaluate_first_order_EDM_supercell_p1 |  construct_first_order_edm_supercell_scalapack
!(7)   add_matrix_to_sparse                  |  get_first_order_edm_sparse_matrix_from_supercell_scalapack

!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/initialize_scalapack_DFPT_phonon
!  NAME
!    initialize_scalapack_DFPT_phonon
!  SYNOPSIS
  subroutine initialize_scalapack_DFPT_phonon()
!  PURPOSE
!    Initialize the ScaLAPACK environment for DFPT_phonon.
!  USES
    use dimensions, only: n_basis_sc_DFPT
    use localorb_io
    use mpi_tasks
    use mpi_utilities, only: get_my_processor
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    ScaLAPACK communicators, block sizes and local storage arrays plur BLACS grids
!    and the local indexing arrays are set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


    integer :: i, j, lc_DFPT_phonon, lr_DFPT_phonon
    integer :: np0, nq0, trilwmin, lwormtr
    integer :: block_size_DFPT_phonon, info
    integer :: mpierr
    character(LEN=MPI_MAX_PROCESSOR_NAME) :: my_proc_name
    integer :: my_proc_name_len

    integer, external :: numroc

    character*200 :: info_str

    !----(1)  deallocate_scalapack using n_basis.----------------
    ! note:   we can also just add DFPT_phonon part but not deallocate.
    !call deallocate_scalapack()

   !-----(2) get new BLACS_Grid for DFPT_phonon
    call get_my_processor(my_proc_name, my_proc_name_len)
    my_proc_name(my_proc_name_len+1:) = ' ' ! for safety: pad with blanks

    n_scalapack_tasks = n_tasks
    write(info_str, '(2X,A,I0,A)') &
          '* Using ', n_scalapack_tasks, ' tasks for scalapack DFPT_phonon.'
    call localorb_info(info_str, use_unit, '(A)')

    call MPI_Comm_size( mpi_comm_global, n_scalapack_tasks, mpierr)
    call MPI_Comm_rank( mpi_comm_global, my_scalapack_id_DFPT_phonon, mpierr) ! the same as myid

    ! divide the BLACS grid into rows and columns for each task
    do npcol_DFPT_phonon = NINT(sqrt(dble(n_scalapack_tasks))), 2, -1
       if (MOD(n_scalapack_tasks,npcol_DFPT_phonon)==0) exit
    enddo
    nprow_DFPT_phonon = n_scalapack_tasks/npcol_DFPT_phonon ! always succeeds without remainder

    if (my_scalapack_id_DFPT_phonon==0)  then
         write(info_str, '(4(a,i6))') ' DFPT_phonon: Tasks:',n_scalapack_tasks, &
               ' split into ',nprow_DFPT_phonon,' X ',npcol_DFPT_phonon,' BLACS grid'
       call localorb_info(info_str, use_unit, '(A)')
    end if

   ! initialize the BLACS grid
    if(my_scalapack_id_DFPT_phonon<npcol_DFPT_phonon*nprow_DFPT_phonon) then
      my_blacs_ctxt_DFPT_phonon = mpi_comm_global
      call BLACS_Gridinit( my_blacs_ctxt_DFPT_phonon, 'R', nprow_DFPT_phonon, npcol_DFPT_phonon )
      call BLACS_Gridinfo( my_blacs_ctxt_DFPT_phonon, nprow_DFPT_phonon, npcol_DFPT_phonon, &
                           myprow_DFPT_phonon, mypcol_DFPT_phonon )
    else
      myprow_DFPT_phonon = -1
      mypcol_DFPT_phonon = -1
    endif

    write(info_str, '(a)') '-------------shanghui begin debug-----------------'
       call localorb_info(info_str, use_unit, '(A)')
    write(use_unit,*) '|my_scalapack_grid:',my_scalapack_id_DFPT_phonon, myprow_DFPT_phonon, mypcol_DFPT_phonon
    write(info_str, '(a)') '-------------shanghui end debug-----------------'
       call localorb_info(info_str, use_unit, '(A)')


    !--------(3) n_basis_sc_DFPT relateded reallocate
    ! Allocate and set global process mapping info

    write(info_str, *) ' Calculating blocksize based on n_basis_sc_DFPT = ',n_basis_sc_DFPT, &
                       ' max_nprow = ',nprow_DFPT_phonon,' max_npcol = ',npcol_DFPT_phonon
    call localorb_info(info_str, use_unit, '(A)')

!---------------------------begin block_size_DFPT_phonon------------------------------------------------
    block_size_DFPT_phonon = 1 ! Minimum permitted size
    if(block_size_DFPT_phonon*MAX(nprow_DFPT_phonon,npcol_DFPT_phonon) > n_basis_sc_DFPT) then
       write(info_str, *) 'ERROR: n_basis_sc_DFPT = ',n_basis_sc_DFPT,' too small for this processor grid'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    endif
    ! Increase blocksize to maximum possible size or 64

!================begin shanghui note=======================
!  n_block_DFPT_phonon is always set to 1 at present, becaue,
!  have tested n_block_DFPT_phonon = 1 (246 s) is faster than  2 (254 s),
!  you can uncomment the following if you need large n_block_DFPT_phonon
!================end shanghui note=======================
!!!    do while (2*block_size_DFPT_phonon*MAX(nprow_DFPT_phonon,npcol_DFPT_phonon) &
!!!             <= n_basis_sc_DFPT .and. block_size_DFPT_phonon<64)
!!!       block_size_DFPT_phonon = 2*block_size_DFPT_phonon
!!!    end do
!---------------------------begin block_size_DFPT_phonon------------------------------------------------

    nb_DFPT_phonon = block_size_DFPT_phonon
    mb_DFPT_phonon = block_size_DFPT_phonon
    write(info_str, *) 'DFPT_phonon: Scalapack blocksize set to: ',block_size_DFPT_phonon
    call localorb_info(info_str, use_unit, '(A)')

    ! initialize the Scalapack descriptor

    if(my_scalapack_id_DFPT_phonon<npcol_DFPT_phonon*nprow_DFPT_phonon) then

       mxld_DFPT_phonon =  numroc( n_basis_sc_DFPT, mb_DFPT_phonon, myprow_DFPT_phonon, &
                                   rsrc_DFPT_phonon, nprow_DFPT_phonon )
       mxcol_DFPT_phonon = numroc( n_basis_sc_DFPT, nb_DFPT_phonon, mypcol_DFPT_phonon, &
                                   csrc_DFPT_phonon, npcol_DFPT_phonon )

! RJ: If mxld/mxcol are too small, they *might* trigger an error in the
! Intel/Scalapack implementation, so set them to a at least 64:

       if(mxld_DFPT_phonon  < 64) mxld_DFPT_phonon  = 64
       if(mxcol_DFPT_phonon < 64) mxcol_DFPT_phonon = 64

       call descinit( sc_desc_DFPT_phonon, n_basis_sc_DFPT, n_basis_sc_DFPT, &
                      mb_DFPT_phonon, nb_DFPT_phonon, rsrc_DFPT_phonon, csrc_DFPT_phonon, &
                      my_blacs_ctxt_DFPT_phonon, MAX(1,mxld_DFPT_phonon), info )

       ! Safety check only, the following should never happen
       if(mxld_DFPT_phonon<=0 .or. mxcol_DFPT_phonon<=0) then
          write(use_unit,*) 'ERROR Task #', myid,' mxld_DFPT_phonon= ', mxld_DFPT_phonon, &
                            ' mxcol_DFPT_phonon= ', mxcol_DFPT_phonon
          call mpi_abort(mpi_comm_global,1,mpierr)
       endif

    else
       mxld_DFPT_phonon = 1
       mxcol_DFPT_phonon = 1
    endif



    ! allocate and set index arrays

    allocate(l_row_DFPT_phonon(n_basis_sc_DFPT),stat=info)
    call check_allocation(info, 'l_row_DFPT_phonon                         ')
    allocate(l_col_DFPT_phonon(n_basis_sc_DFPT),stat=info)
    call check_allocation(info, 'l_col_DFPT_phonon                         ')

    ! Mapping of global rows/cols to local

    l_row_DFPT_phonon(:) = 0
    l_col_DFPT_phonon(:) = 0

    ! ATTENTION: The following code assumes rsrc==0 and csrc==0 !!!!
    ! For processors outside the working set, l_row/l_col will stay
    ! completely at 0

    lr_DFPT_phonon = 0 ! local row counter
    lc_DFPT_phonon = 0 ! local column counter


    do i = 1, n_basis_sc_DFPT
      if( MOD((i-1)/mb_DFPT_phonon,nprow_DFPT_phonon) == myprow_DFPT_phonon) then
        ! row i is on local processor
        lr_DFPT_phonon = lr_DFPT_phonon + 1
        l_row_DFPT_phonon(i) = lr_DFPT_phonon
      endif

      if( MOD((i-1)/nb_DFPT_phonon,npcol_DFPT_phonon) == mypcol_DFPT_phonon) then
        ! column i is on local processor
        lc_DFPT_phonon = lc_DFPT_phonon + 1
        l_col_DFPT_phonon(i) = lc_DFPT_phonon
      endif
    enddo



    ! Mapping of local rows/cols to global

    n_my_rows_DFPT_phonon = lr_DFPT_phonon
    n_my_cols_DFPT_phonon = lc_DFPT_phonon
    allocate(my_row_DFPT_phonon(n_my_rows_DFPT_phonon))
    allocate(my_col_DFPT_phonon(n_my_cols_DFPT_phonon))
    lr_DFPT_phonon = 0
    lc_DFPT_phonon = 0

    do i = 1, n_basis_sc_DFPT
       if(l_row_DFPT_phonon(i)>0) then
          lr_DFPT_phonon = lr_DFPT_phonon + 1
          my_row_DFPT_phonon(lr_DFPT_phonon) = i
       endif
       if(l_col_DFPT_phonon(i)>0) then
          lc_DFPT_phonon = lc_DFPT_phonon + 1
          my_col_DFPT_phonon(lc_DFPT_phonon) = i
       endif
    enddo


    ! Allocate scalapack arrays
    allocate(ovlp_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon),stat=info)
    call check_allocation(info, 'ovlp_supercell_scalapack             ')
    allocate(ham_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon,n_spin),stat=info)
    call check_allocation(info, 'ham_supercell_scalapack             ')
    allocate(eigenvec_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon,n_spin),stat=info)
    call check_allocation(info, 'eigenvec_supercell_scalapack     ')
    allocate(eigenvalues_supercell_scalapack(n_basis_sc_DFPT, n_spin),stat=info)
    call check_allocation(info, 'eigenvalues_supercell_scalapack     ')

    allocate(first_order_ovlp_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon),stat=info)
    call check_allocation(info, 'first_order_ovlp_supercell_scalapack             ')
    allocate(first_order_ham_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon),stat=info)
    call check_allocation(info, 'first_order_ham_supercell_scalapack             ')
    allocate(first_order_U_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon),stat=info)
    call check_allocation(info, 'first_order_U_supercell_scalapack             ')
    allocate(first_order_edm_supercell_scalapack(mxld_DFPT_phonon,mxcol_DFPT_phonon),stat=info)
    call check_allocation(info, 'first_order_edm_supercell_scalapack             ')
      ! Safety only:
    ovlp_supercell_scalapack     = 0
    ham_supercell_scalapack      = 0
    eigenvec_supercell_scalapack = 0
    eigenvalues_supercell_scalapack = 0

    first_order_ovlp_supercell_scalapack     = 0
    first_order_ham_supercell_scalapack      = 0
    first_order_U_supercell_scalapack        = 0
    first_order_edm_supercell_scalapack      = 0



    ! Calculate workspace needed for eigenvalue solver
    np0 = NUMROC( MAX(n_basis_sc_DFPT,nb_DFPT_phonon,2), nb_DFPT_phonon, 0, 0, nprow_DFPT_phonon )
    nq0 = NUMROC( MAX(n_basis_sc_DFPT,nb_DFPT_phonon,2), nb_DFPT_phonon, 0, 0, npcol_DFPT_phonon )
    TRILWMIN = 3*n_basis_sc_DFPT+ MAX( nb_DFPT_phonon*( NP0+1 ), 3*nb_DFPT_phonon )
    lwormtr = MAX( (nb_DFPT_phonon*(nb_DFPT_phonon-1))/2, &
                   (np0 + nq0)*nb_DFPT_phonon + 2*nb_DFPT_phonon*nb_DFPT_phonon)
    lrwork_DFPT_phonon = MAX( 1+6*n_basis_sc_DFPT+2*NP0*NQ0, TRILWMIN, lwormtr ) + 2*n_basis_sc_DFPT

    liwork_DFPT_phonon = MAX(7*n_basis_sc_DFPT + 8*npcol_DFPT_phonon + 2, n_basis_sc_DFPT + &
                                     2*nb_DFPT_phonon + 2*npcol_DFPT_phonon)
    lcwork_DFPT_phonon = 0 ! no complex workspace needed

    len_scalapack_work_DFPT_phonon = lrwork_DFPT_phonon ! Total workspace (real numbers)
    write(info_str, *) ' Required Scalapack workspace - INTEGER: ',liwork_DFPT_phonon, &
                       ' REAL:  ',lrwork_DFPT_phonon
    call localorb_info(info_str, use_unit, '(A)')

  end subroutine initialize_scalapack_DFPT_phonon
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/deallocate_scalapack_DFPT_phonon
!  NAME
!    deallocate_scalapack_DFPT_phonon
!  SYNOPSIS
  subroutine deallocate_scalapack_DFPT_phonon()
!  PURPOSE
!    Deallocates ScaLAPACK arrays used in DFPT_phonon.
!  USES
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    ScaLAPACK arrays are deallocated on exit.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE

    if (allocated(l_row_DFPT_phonon))    deallocate(l_row_DFPT_phonon)
    if (allocated(l_col_DFPT_phonon))    deallocate(l_col_DFPT_phonon)

    if (allocated(my_row_DFPT_phonon))    deallocate(my_row_DFPT_phonon)
    if (allocated(my_col_DFPT_phonon))    deallocate(my_col_DFPT_phonon)

    if (allocated(ovlp_supercell_scalapack))     deallocate(ovlp_supercell_scalapack)
    if (allocated(ham_supercell_scalapack))      deallocate(ham_supercell_scalapack)
    if (allocated(first_order_ovlp_supercell_scalapack))     deallocate(first_order_ovlp_supercell_scalapack)
    if (allocated(first_order_ham_supercell_scalapack))      deallocate(first_order_ham_supercell_scalapack)
    if (allocated(first_order_U_supercell_scalapack))        deallocate(first_order_U_supercell_scalapack)
    if (allocated(first_order_edm_supercell_scalapack))     deallocate(first_order_edm_supercell_scalapack)
    if (allocated(eigenvec_supercell_scalapack)) deallocate(eigenvec_supercell_scalapack)
    if (allocated(eigenvalues_supercell_scalapack)) deallocate(eigenvalues_supercell_scalapack)

  end subroutine deallocate_scalapack_DFPT_phonon
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/finalize_scalapack_DFPT_phonon
!  NAME
!    finalize_scalapack_DFPT_phonon
!  SYNOPSIS
  subroutine finalize_scalapack_DFPT_phonon
!  PURPOSE
!    Removes the BLACS grid used in DFPT_phonon.
!  USES
  implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    if(my_scalapack_id_DFPT_phonon<npcol_DFPT_phonon*nprow_DFPT_phonon) then
      call BLACS_Gridexit(my_blacs_ctxt_DFPT_phonon)
      call BLACS_Exit(1)
    endif

    call deallocate_scalapack_DFPT_phonon

  end subroutine finalize_scalapack_DFPT_phonon


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_overlap_supercell_scalapack
!  NAME
!    construct_overlap_supercell_scalapack
!  SYNOPSIS
  subroutine construct_overlap_supercell_scalapack( overlap_matrix )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage using n_basis_sc_DFPT.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use runtime_choices, only: use_local_index

    implicit none
!  ARGUMENTS
    real*8:: overlap_matrix  ( n_hamiltonian_matrix_size )
!  INPUTS
!    o overlap_matrix -- the overlap matrix
!  OUTPUT
!    lower half of the ScaLAPACK array overlap_matrix is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: i_spin, i_cell_in_hamiltonian, iuo, io, juo, jo, i_index
    integer :: lr, lc
    integer :: i_cell_in_sc_DFPT, i_cell_trans, j_cell_trans, io_trans, jo_trans

    ! construct_overlap_supercell_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_overlap_supercell_scalapack + use_local_index")

    !--------in DFPT_phonon, the ham(n_basis_sc_DFPT,n_basis_sc_DFPT) is always real.
     ovlp_supercell_scalapack(:,:) = 0.

    ! Attention: Only the lower half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'L'!


       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell_in_hamiltonian = 1, n_cells_in_hamiltonian-1
          do iuo = 1, n_basis

             i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian)
             io  = cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

             if(index_hamiltonian(1,i_cell_in_hamiltonian,iuo) > 0) then
               do i_index = index_hamiltonian(1,i_cell_in_hamiltonian,iuo),  &
                            index_hamiltonian(2,i_cell_in_hamiltonian,iuo)

                juo = column_index_hamiltonian(i_index)
                lr  = l_row_DFPT_phonon(io)   ! local row number
                lc  = l_col_DFPT_phonon(juo)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                   ovlp_supercell_scalapack(lr,lc) = overlap_matrix(i_index)
                endif

                !--------begin the symmetric part-------
                lr  = l_row_DFPT_phonon(juo)   ! local row number
                lc  = l_col_DFPT_phonon(io)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                   ovlp_supercell_scalapack(lr,lc) = overlap_matrix(i_index)
                endif
                !--------end the symmetric part-------

                  do j_cell_trans = 2,n_cells_in_sc_DFPT
                     ! jo_trans = juo +j_cell
                     jo_trans = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)
                     i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans)
                     ! io_trans = iuo + i_cell_in_sc_DFPT + j_cell_in_sc_DFPT
                     io_trans = cell_and_basis_to_basis_sc_DFPT( i_cell_trans,iuo)

                     lr  = l_row_DFPT_phonon(io_trans)
                     lc  = l_col_DFPT_phonon(jo_trans)
                     if(lr.ne.0.and.lc.ne.0) then
                        ovlp_supercell_scalapack(lr,lc) = overlap_matrix(i_index)
                     endif

                     !--------begin the symmetric part-------
                     lr  = l_row_DFPT_phonon(jo_trans)   ! local row number
                     lc  = l_col_DFPT_phonon(io_trans)  ! local column number
                     if(lr.ne.0.and.lc.ne.0) then
                        ovlp_supercell_scalapack(lr,lc) = overlap_matrix(i_index)
                     endif
                     !--------end the symmetric part-------

                  enddo
               end do
             end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_overlap_supercell_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format

   factor_overlap_DFPT_phonon = .true.

  end subroutine construct_overlap_supercell_scalapack

!--------------for debug, need to be remove-----------
! subroutine get_overlap_sparse_matrix_from_supercell_scalapack( matrix_sparse)
! !  PURPOSE
! !    Gets a sparse matrix array from ScaLAPACK
! !  USES
!   use pbc_lists
!   implicit none
! !  ARGUMENTS
!   real*8 :: matrix_sparse(n_hamiltonian_matrix_size)
!   integer :: i_spin
!
!   integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc
!   integer :: i_cell_in_sc_DFPT, io
!
!   character*200 :: info_str
!
!   matrix_sparse = 0.0d0
!
!   do i_cell = 1,n_cells_in_hamiltonian-1
!   do i_bas1 = 1, n_basis
!
!      i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell)
!      io  = cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,i_bas1)
!      lr = l_row_DFPT_phonon(io)
!
!         if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries
!         do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)
!
!            i_bas2 = column_index_hamiltonian(i_index)
!            lc = l_col_DFPT_phonon(i_bas2) ! local row number
!
!               if(lr.ne.0.and.lc.ne.0) then
!               ! debug for ovlp
!               !matrix_sparse(i_index) = ovlp_supercell_scalapack(lr,lc)
!               ! debug for ham
!               matrix_sparse(i_index) = ham_supercell_scalapack(lr,lc,1)
!               endif
!
!         end do
!      end do
!   end do
!
! call sync_density_matrix_sparse(matrix_sparse)
!
! end subroutine get_overlap_sparse_matrix_from_supercell_scalapack



!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_hamiltonian_supercell_scalapack
!  NAME
!    construct_hamiltonian_supercell_scalapack
!  SYNOPSIS
  subroutine construct_hamiltonian_supercell_scalapack( hamiltonian )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage using n_basis_sc_DFPT.
!  USES
    use pbc_lists
    use geometry
    use basis
    use mpi_tasks
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    real*8:: hamiltonian   ( n_hamiltonian_matrix_size, n_spin )
!  INPUTS
!    o hamiltonian -- the Hamilton matrix
!  OUTPUT
!    lower half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: i_spin, i_cell_in_hamiltonian, iuo, io, juo, jo, i_index
    integer :: lr, lc
    integer :: i_cell_in_sc_DFPT, i_cell_trans, j_cell_trans, io_trans, jo_trans

    ! construct_hamiltonian_supercell_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_supercell_scalapack + use_local_index")

    !--------in DFPT_phonon, the ham(n_basis_sc_DFPT,n_basis_sc_DFPT) is always real.
     ham_supercell_scalapack(:,:,:) = 0.

    ! Attention: Only the lower half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'L'!

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell_in_hamiltonian = 1, n_cells_in_hamiltonian-1
          do iuo = 1, n_basis

             i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian)
             io  = cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

             if(index_hamiltonian(1,i_cell_in_hamiltonian,iuo) > 0) then
               do i_index = index_hamiltonian(1,i_cell_in_hamiltonian,iuo),  &
                            index_hamiltonian(2,i_cell_in_hamiltonian,iuo)

                juo = column_index_hamiltonian(i_index)
                lc  = l_col_DFPT_phonon(juo)  ! local column number
                lr  = l_row_DFPT_phonon(io)   ! local row number
                if(lr.ne.0.and.lc.ne.0) then
                  ham_supercell_scalapack(lr,lc,i_spin) = hamiltonian(i_index,i_spin)
                endif

                !--------begin the symmetric part-------
                lr  = l_row_DFPT_phonon(juo)   ! local row number
                lc  = l_col_DFPT_phonon(io)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                   ham_supercell_scalapack(lr,lc,i_spin) = hamiltonian(i_index,i_spin)
                endif
                !--------end the symmetric part-------

                  do j_cell_trans = 2,n_cells_in_sc_DFPT
                     ! jo_trans = juo +j_cell
                     jo_trans = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)
                     i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans)
                     ! io_trans = iuo + i_cell_in_sc_DFPT + j_cell_in_sc_DFPT
                     io_trans = cell_and_basis_to_basis_sc_DFPT( i_cell_trans,iuo)

                     lr  = l_row_DFPT_phonon(io_trans)
                     lc  = l_col_DFPT_phonon(jo_trans)
                     if(lr.ne.0.and.lc.ne.0) then
                      ham_supercell_scalapack(lr,lc,i_spin) = hamiltonian(i_index,i_spin)
                     endif

                     !--------begin the symmetric part-------
                     lr  = l_row_DFPT_phonon(jo_trans)   ! local row number
                     lc  = l_col_DFPT_phonon(io_trans)  ! local column number
                     if(lr.ne.0.and.lc.ne.0) then
                      ham_supercell_scalapack(lr,lc,i_spin) = hamiltonian(i_index,i_spin)
                     endif
                     !--------end the symmetric part-------

                  enddo
               end do
             end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_hamiltonian_supercell_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format
    end do ! i_spin

  end subroutine construct_hamiltonian_supercell_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/solve_evp_supercell_scalapack
!  NAME
!    solve_evp_supercell_scalapack
!  SYNOPSIS
  subroutine solve_evp_supercell_scalapack( )
!  PURPOSE
!    Solve evp in the ScaLAPACK array.
!  USES
    use dimensions, only: n_basis_sc_DFPT
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use localorb_io
    implicit none
!  ARGUMENTS
!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!  ARGUMENTS
!  INPUTS
!  OUTPUT
!    o
!    o
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: nwork, n
    integer :: i_row, i_col
    integer, dimension(dlen_DFPT_phonon) :: ns_desc_DFPT_phonon

    integer :: mpierr

    character*100 :: info_str
    character*256 :: filename

    real*8, allocatable :: tmp(:,:)
    integer :: i_spin

    integer,     parameter :: ibtype = 1
    character*1, parameter :: jobz = 'V'
    character*1, parameter :: range = 'A'
    character*1, parameter :: uplo = 'L'
    real*8      abstol_tiny
    real*8      abstol
    integer :: m, nz
    integer :: ifail(n_basis_sc_DFPT)
    integer, allocatable ::  iclustr(:)
    real*8 , allocatable ::  gap(:)
    integer :: info

    real*8 , allocatable ::  work(:)
    integer, allocatable :: iwork(:)
    integer lwork, liwrok

    !==========shanghui begin for debug========
    !external  PDLAPRNT
    !real*8 , allocatable ::  work(:)
    !==========shanghui end for debug==========

    write (info_str,'(2X,A,A)') &
         "| DFPT_phonon : Solving real symmetric generalised eigenvalue problem ", &
         "by standard ScaLAPACK."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)


!------begain scalapack version-----------------------
    allocate(iclustr(2*npcol_DFPT_phonon*nprow_DFPT_phonon))
    allocate(gap(npcol_DFPT_phonon*nprow_DFPT_phonon))

    do i_spin = 1, n_spin


      abstol=tiny(abstol_tiny)
!      If jobz = 'V', then the amount of workspace required to guarantee that all eigenvectors are computed is:
!For PDSYGVX, lwork ≥ 3(n+ja-1)+2n +max(5nn, (np0)(mq0)+2(nb)(nb)) + iceil(neig, (nprow)(npcol))(nn)
      allocate(work(1))
      allocate(iwork(1))
      CALL pdsygvx(ibtype, 'V', 'A', 'L', n_basis_sc_DFPT,  &
           ham_supercell_scalapack(:,:,i_spin), 1, 1, sc_desc_DFPT_phonon, &
           ovlp_supercell_scalapack(:,:), 1, 1, sc_desc_DFPT_phonon, &
           0.0d0, 0.0d0, 0.00d0, 0.00d0, abstol, m, nz, eigenvalues_supercell_scalapack(:,i_spin), &
           -1.0d0, eigenvec_supercell_scalapack(:,:,i_spin), &
           1, 1, sc_desc_DFPT_phonon, work, -1, iwork, -1, &
           ifail, iclustr, gap, info)

      lwork = work(1) !lwork = 100000
      deallocate(work)
      allocate(work(lwork))

      liwork = iwork(1) ! liwork = 5000
      deallocate(iwork)
      allocate(iwork(liwork))
      CALL pdsygvx(ibtype, 'V', 'A', 'L', n_basis_sc_DFPT,  &
           ham_supercell_scalapack(:,:,i_spin), 1, 1, sc_desc_DFPT_phonon, &
           ovlp_supercell_scalapack(:,:), 1, 1, sc_desc_DFPT_phonon, &
           0.0d0, 0.0d0, 0.00d0, 0.00d0, abstol, m, nz, eigenvalues_supercell_scalapack(:,i_spin), &
           -1.0d0, eigenvec_supercell_scalapack(:,:,i_spin), &
           1, 1, sc_desc_DFPT_phonon, work, lwork, iwork, liwork, &
           ifail, iclustr, gap, info)

         if(myid.eq.0) then
         write(use_unit,*) 'scalapack eigenvalues:', eigenvalues_supercell_scalapack(:,i_spin)
         endif

         deallocate(work)
         deallocate(iwork)

     enddo ! i_spin

     ! !======shanghui begin debug===============
     ! --------to make sure eigen solver is right----------
     !
     ! allocate(work(1))
     ! CALL PDLAPRNT( n_basis_sc_DFPT, n_basis_sc_DFPT, eigenvec_supercell_scalapack(:,:,1), &
     !                1, 1, sc_desc_DFPT_phonon, 0, 0, 'C0', 6, work )
     ! deallocate(work)

     !  !(1) print CT H C
     !  tmp_1 = 0.0d0
     !  call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
     !              1.0d0, eigenvec_supercell_scalapack(:,:,1),  1, 1, sc_desc_DFPT_phonon,  &  ! alpha, a, ia, ja, desc_a
     !                     ham_supercell_scalapack(:,:,1), 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
     !              0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c
     !  tmp_H1 = 0.0d0
     !  call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
     !              1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
     !                     eigenvec_supercell_scalapack(:,:,1), 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
     !              0.0d0, tmp_H1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
     ! allocate(work(1))
     ! CALL PDLAPRNT( n_basis_sc_DFPT, n_basis_sc_DFPT, tmp_H1, &
     !                1, 1, sc_desc_DFPT_phonon, 0, 0, 'CT H0 C', 6, work )
     ! deallocate(work)


     !  !(2) print CT S C
     !  tmp_1 = 0.0d0
     !  call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
     !              1.0d0, eigenvec_supercell_scalapack(:,:,1),  1, 1, sc_desc_DFPT_phonon,  &  ! alpha, a, ia, ja, desc_a
     !                     ovlp_supercell_scalapack(:,:), 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
     !              0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c
     !  tmp_S1 = 0.0d0
     !  call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
     !              1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
     !                     eigenvec_supercell_scalapack(:,:,1), 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
     !              0.0d0, tmp_S1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
     ! allocate(work(1))
     ! CALL PDLAPRNT( n_basis_sc_DFPT, n_basis_sc_DFPT, tmp_S1, &
     !                1, 1, sc_desc_DFPT_phonon, 0, 0, 'CT S0 C', 6, work )
     ! deallocate(work)


     ! deallocate(tmp_S1)
     ! deallocate(tmp_H1)
     ! deallocate(tmp_1)

     ! !======shanghui end debug===============


     !=================shanghui begin debug===============
     !---------using scalapack test code------------------
     !
     !  !-------begin print whole matrix----------------
     !  allocate(work(1))
     !  CALL PDLAPRNT( n_basis_sc_DFPT, n_basis_sc_DFPT, ovlp_supercell_scalapack, &
     !                 1, 1, sc_desc_DFPT_phonon, 0, 0, 'S1_supercell', 6, work )
     !  CALL PDLAPRNT( n_basis_sc_DFPT, n_basis_sc_DFPT, ham_supercell_scalapack(:,:,i_spin), &
     !                 1, 1, sc_desc_DFPT_phonon, 0, 0, 'H1_supercell', 6, work )
     !  deallocate(work)
     !  !-------end print whole matrix----------------
     ! CALL PDLAMODHILB( n_basis_sc_DFPT, ham_supercell_scalapack(:,:,i_spin), 1, 1, &
     !                   sc_desc_DFPT_phonon, info )
     ! CALL PDLATOEPLITZ( n_basis_sc_DFPT, ovlp_supercell_scalapack, 1, 1, &
     !                    sc_desc_DFPT_phonon, info )
     !
     ! lwork = 100000
     ! allocate(work(lwork))
     ! liwork = 5000
     ! allocate(iwork(liwork))
     !
     ! CALL PDSYGVX( ibtype, 'V', 'A', 'U', n_basis_sc_DFPT, &
     !               ham_supercell_scalapack(:,:,i_spin), 1, 1, sc_desc_DFPT_phonon, &
     !               ovlp_supercell_scalapack(:,:), 1, 1, sc_desc_DFPT_phonon, &
     !               0.0d0, 0.0d0, 13, -13, -1.0d0, m, nz,  eigenvalues_supercell_scalapack(:,i_spin),   &
     !               -1.0d0, eigenvec_supercell_scalapack(:,:,i_spin), &
     !                1, 1, sc_desc_DFPT_phonon, work, lwork, iwork, liwork,   &
     !               IFAIL, ICLUSTR, GAP, INFO )
     !    if(myid.eq.0) then
     !      do i_row = 1, n_basis_sc_DFPT
     !      write(use_unit,*) ' W(', i_row,')=', eigenvalues_supercell_scalapack(i_row,1)
     !      enddo
     !    endif
      !-------begin shanghui test code--------------



     deallocate(iclustr)
     deallocate(gap)
!------end scalapack version-----------------------


   end subroutine solve_evp_supercell_scalapack


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_overlap_supercell_scalapack
!  NAME
!    construct_first_order_overlap_supercell_scalapack
!  SYNOPSIS
  subroutine construct_first_order_overlap_supercell_scalapack( k_coord,k_atom,first_order_S_sparse )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage using n_basis_sc_DFPT.
!  USES
    use dimensions, only: n_centers_in_sc_DFPT
    use runtime_choices
    use pbc_lists
    use geometry
    use basis
    use mpi_tasks
    implicit none
!  ARGUMENTS
    integer, intent(in) :: k_coord
    integer, intent(in) :: k_atom
    real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size), intent(in) :: first_order_S_sparse
!  INPUTS
!    o overlap_matrix -- the overlap matrix
!  OUTPUT
!    upper half of the ScaLAPACK array overlap_matrix is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: i_spin, i_cell_in_hamiltonian, iuo, io, juo, jo, i_index
    integer :: lr, lc
    integer :: i_cell_in_sc_DFPT, i_cell_trans, j_cell_trans, io_trans, jo_trans
    integer :: k_cell_trans, k_center_trans

    ! construct_first_order_overlap_supercell_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_first_order_overlap_supercell_scalapack + use_local_index")

    !--------in DFPT_phonon, the ham(n_basis_sc_DFPT,n_basis_sc_DFPT) is always real.
     first_order_ovlp_supercell_scalapack(:,:) = 0.

    ! Attention: Only the lower half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'L'!


       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do iuo = 1, n_basis
          do i_cell_in_hamiltonian = 1, n_cells_in_hamiltonian-1

             i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian)
             io  = cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

             !------------------(1) j_cell = 1 ----------------------------
             !                                                            !
             !       d M_full (u R1, v)          d M_sparse (u R1, v)     !
             !    ----------------------- =    -----------------------    !
             !        d R_K                      d R_K                    !
             !                                                            !

             if(index_hamiltonian(1,i_cell_in_hamiltonian,iuo) > 0) then
             do i_index =  index_hamiltonian(1,i_cell_in_hamiltonian,iuo),  &
                           index_hamiltonian(2,i_cell_in_hamiltonian,iuo)

                juo = column_index_hamiltonian(i_index)
                lr  = l_row_DFPT_phonon(io)   ! local row number
                lc  = l_col_DFPT_phonon(juo)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                  first_order_ovlp_supercell_scalapack(lr,lc) = &
                  first_order_S_sparse(k_coord, k_atom, i_index)
                endif

               !--------begin the symmetric part-------
                lr  = l_row_DFPT_phonon(juo)   ! local row number
                lc  = l_col_DFPT_phonon(io)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                  first_order_ovlp_supercell_scalapack(lr,lc) = &
                  first_order_S_sparse(k_coord, k_atom, i_index)
                endif

               !--------end the symmetric part-------

            !------------------(2) j_cell > 1 ----------------------------
            !                                                            !
            !    d M_full (u R1+R2, v+R2)      d M_sparse (u R1, v)      !
            !    ----------------------- =    -----------------------    !
            !        d R_K [k_center]       d R_K-R2 [k_center_trans]    !
            !                                                            !

                 do j_cell_trans = 2,n_cells_in_sc_DFPT

                    ! jo_trans = juo +j_cell
                    jo_trans = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)

                    ! io_trans = iuo + i_cell + j_cell
                    i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans)
                    io_trans = cell_and_basis_to_basis_sc_DFPT(i_cell_trans,iuo)

                    ! k_center_trans = k_atom_cell - j_cell
                    k_cell_trans   = cell_diff_sc_DFPT(1,j_cell_trans)
                    k_center_trans = cell_and_atom_to_center_sc_DFPT(k_cell_trans,k_atom)

                    lr  = l_row_DFPT_phonon(io_trans)
                    lc  = l_col_DFPT_phonon(jo_trans)
                    if(lr.ne.0.and.lc.ne.0) then
                      first_order_ovlp_supercell_scalapack(lr,lc) = &
                      first_order_S_sparse(k_coord, k_center_trans, i_index)
                    endif

                   !--------begin the symmetric part-------
                    lr  = l_row_DFPT_phonon(jo_trans)
                    lc  = l_col_DFPT_phonon(io_trans)
                    if(lr.ne.0.and.lc.ne.0) then
                      first_order_ovlp_supercell_scalapack(lr,lc) = &
                      first_order_S_sparse(k_coord, k_center_trans, i_index)
                    endif
                   !--------end the symmetric part-------

                 enddo

             end do !i_index
             end if

          end do ! i_cell_in_hamiltonian
          end do ! iuo

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_first_order_overlap_supercell does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format


  end subroutine construct_first_order_overlap_supercell_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_hamiltonian_supercell_scalapack
!  NAME
!    construct_first_order_hamiltonian_supercell_scalapack
!  SYNOPSIS
  subroutine construct_first_order_hamiltonian_supercell_scalapack( k_coord,k_atom,first_order_H_sparse )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage using n_basis_sc_DFPT.
!  USES
    use dimensions, only: n_centers_in_sc_DFPT
    use runtime_choices
    use pbc_lists
    use geometry
    use basis
    use mpi_tasks
    implicit none
!  ARGUMENTS
    integer, intent(in) :: k_coord
    integer, intent(in) :: k_atom
    real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size), intent(in) :: first_order_H_sparse

!  INPUTS
!    o overlap_matrix -- the overlap matrix
!  OUTPUT
!    upper half of the ScaLAPACK array overlap_matrix is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: i_spin, i_cell_in_hamiltonian, iuo, io, juo, jo, i_index
    integer :: lr, lc
    integer :: i_cell_in_sc_DFPT, i_cell_trans, j_cell_trans, io_trans, jo_trans
    integer :: k_cell_trans, k_center_trans

    ! construct_first_order_overlap_supercell_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_first_order_overlap_supercell_scalapack + use_local_index")

    !--------in DFPT_phonon, the ham(n_basis_sc_DFPT,n_basis_sc_DFPT) is always real.
     first_order_ham_supercell_scalapack(:,:) = 0.

    ! Attention: Only the lower half of the matrices is set in the code below,
    ! so the parameter UPLO to the Scalapack routines must always be 'L'!

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do iuo = 1, n_basis
          do i_cell_in_hamiltonian = 1, n_cells_in_hamiltonian-1

             i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian)

             io  = cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

             !------------------(1) j_cell = 1 ----------------------------
             !                                                            !
             !       d M_full (u R1, v)          d M_sparse (u R1, v)     !
             !    ----------------------- =    -----------------------    !
             !        d R_K                      d R_K                    !
             !                                                            !


             if(index_hamiltonian(1,i_cell_in_hamiltonian,iuo) > 0) then
             do i_index =  index_hamiltonian(1,i_cell_in_hamiltonian,iuo),  &
                           index_hamiltonian(2,i_cell_in_hamiltonian,iuo)

                juo = column_index_hamiltonian(i_index)
                lr  = l_row_DFPT_phonon(io)   ! local row number
                lc  = l_col_DFPT_phonon(juo)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                  first_order_ham_supercell_scalapack(lr,lc) = &
                  first_order_H_sparse(k_coord, k_atom, i_index)
                endif

               !--------begin the symmetric part-------
                lr  = l_row_DFPT_phonon(juo)   ! local row number
                lc  = l_col_DFPT_phonon(io)  ! local column number
                if(lr.ne.0.and.lc.ne.0) then
                  first_order_ham_supercell_scalapack(lr,lc) = &
                  first_order_H_sparse(k_coord, k_atom, i_index)
                endif

               !--------end the symmetric part-------

            !------------------(2) j_cell > 1 ----------------------------
            !                                                            !
            !    d M_full (u R1+R2, v+R2)      d M_sparse (u R1, v)      !
            !    ----------------------- =    -----------------------    !
            !        d R_K [k_center]       d R_K-R2 [k_center_trans]    !
            !                                                            !

                 do j_cell_trans = 2,n_cells_in_sc_DFPT

                    ! jo_trans = juo +j_cell
                    jo_trans = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)

                    ! io_trans = iuo + i_cell + j_cell
                    i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans)
                    io_trans = cell_and_basis_to_basis_sc_DFPT(i_cell_trans,iuo)

                    ! k_center_trans = k_atom_cell - j_cell
                    k_cell_trans   = cell_diff_sc_DFPT(1,j_cell_trans)
                    k_center_trans = cell_and_atom_to_center_sc_DFPT(k_cell_trans,k_atom)

                    lr  = l_row_DFPT_phonon(io_trans)
                    lc  = l_col_DFPT_phonon(jo_trans)
                    if(lr.ne.0.and.lc.ne.0) then
                      first_order_ham_supercell_scalapack(lr,lc) = &
                      first_order_H_sparse(k_coord, k_center_trans, i_index)
                    endif

                   !--------begin the symmetric part-------
                    lr  = l_row_DFPT_phonon(jo_trans)
                    lc  = l_col_DFPT_phonon(io_trans)
                    if(lr.ne.0.and.lc.ne.0) then
                      first_order_ham_supercell_scalapack(lr,lc) = &
                      first_order_H_sparse(k_coord, k_center_trans, i_index)
                    endif
                   !--------end the symmetric part-------

                 enddo

             end do !i_index
             end if

          end do ! i_cell_in_hamiltonian
          end do ! iuo

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_first_order_hamiltonian_supercell_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format

  end subroutine construct_first_order_hamiltonian_supercell_scalapack
!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/evaluate_first_order_U_supercell_scalapack
!  NAME
!    evaluate_first_order_U_supercell_scalapack
!  SYNOPSIS
  subroutine evaluate_first_order_U_supercell_scalapack()
!  PURPOSE
!    evaluate first_order_U
!    here C^+ = C*T
!    U1_pq(my_k_point) =  (C^+ S1 C E - C^+ H1 C)_pq/(E_pp-E_qq)
!
!
!  USES
    use dimensions, only: n_basis_sc_DFPT
    use mpi_tasks
    use pbc_lists
    use physics, only : n_electrons
    implicit none
!  ARGUMENTS
!  INPUTS
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, j_state, info
    integer :: i_local_row, i_local_col, i_global_row, i_global_col
    integer :: i_spin

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_eig(:,:)
    real*8, allocatable :: tmp_S1(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)





    ! The work in this routine must be done only on the working set
    if(my_scalapack_id_DFPT_phonon>=npcol_DFPT_phonon*nprow_DFPT_phonon) return

    max_occ_number =  NINT(n_electrons * n_cells_in_sc_DFPT / 2.0d0)  ! n_basis_sc_DFPT / 2 ! todo : extend

    i_spin = 1 ! to be extended

   !------real-(1).prepare C, C_eig, S^(1), H^(1)----------------------
       allocate(tmp_C(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_eig(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C_eig')
       allocate(tmp_S1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_S1')
       allocate(tmp_H1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_H1')
       allocate(tmp_1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_1')

       tmp_C(:,:) = eigenvec_supercell_scalapack(:,:,i_spin)
       tmp_S1(:,:) = first_order_ovlp_supercell_scalapack(:,:)
       tmp_H1(:,:) = first_order_ham_supercell_scalapack(:,:)
       tmp_C_eig(:,:) = 0.0d0


       do i_state = 1, n_basis_sc_DFPT
          if(l_col_DFPT_phonon(i_state)>0) then
             tmp_C_eig(:,l_col_DFPT_phonon(i_state)) = tmp_C(:,l_col_DFPT_phonon(i_state)) * &
              eigenvalues_supercell_scalapack(i_state, i_spin)
          endif
       end do

    !------real-(2) C^+ S1 C E ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c

        tmp_S1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_eig, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_S1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c


    !------real-(3) C^+ H1 C ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1 , 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
        tmp_H1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_H1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c


    !------real-(4) U1_pq = (C^+ S1 C E - C^+ H1 C)_pq/(E_pp-E_qq) -----------------
        first_order_U_supercell_scalapack = 0.0d0

        do  i_local_row = 1, n_my_rows_DFPT_phonon
        do  i_local_col = 1, n_my_cols_DFPT_phonon

            i_global_row = my_row_DFPT_phonon(i_local_row)
            i_global_col = my_col_DFPT_phonon(i_local_col)

            if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

            first_order_U_supercell_scalapack(i_local_row,i_local_col) =  &
           ( tmp_S1(i_local_row,i_local_col) - tmp_H1(i_local_row,i_local_col)) / &
           ( eigenvalues_supercell_scalapack(i_global_row, i_spin) - &
             eigenvalues_supercell_scalapack(i_global_col, i_spin))

            endif

            if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

            first_order_U_supercell_scalapack(i_local_row,i_local_col) =  &
           ( tmp_S1(i_local_row,i_local_col) - tmp_H1(i_local_row,i_local_col)) / &
           ( eigenvalues_supercell_scalapack(i_global_row, i_spin) - &
             eigenvalues_supercell_scalapack(i_global_col, i_spin))

            endif

       enddo
       enddo


       deallocate(tmp_C)
       deallocate(tmp_C_eig)
       deallocate(tmp_S1)
       deallocate(tmp_H1)
       deallocate(tmp_1)


  end subroutine evaluate_first_order_U_supercell_scalapack


!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_dm_supercell_scalapack
!  NAME
!    construct_first_order_dm_supercell_scalapack
!  SYNOPSIS
  subroutine construct_first_order_dm_supercell_scalapack()
!  PURPOSE
!    Construct first_order density matrix supercell in ScaLAPACK
!    here C+ = C*T
!    DM1(my_k_point) = C*occ_number* (-C^+ S C) * C
!                      + (C U) C^+*occ_number + C*occ_number (C U)^+
!
!  USES
    use dimensions, only: n_basis_sc_DFPT
    use mpi_tasks
    use pbc_lists
    use physics, only : n_electrons
    implicit none
!  ARGUMENTS
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col
    integer :: i_spin

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_S1(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_2(:,:)
    real*8, allocatable :: tmp_C1(:,:)


    ! The work in this routine must be done only on the working set
    if(my_scalapack_id_DFPT_phonon>=npcol_DFPT_phonon*nprow_DFPT_phonon) return

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix

    max_occ_number =  NINT(n_electrons * n_cells_in_sc_DFPT / 2.0d0)  ! n_basis_sc_DFPT / 2 ! todo : extend

    if(mod(nint(n_electrons),2).ne.0) then
      call aims_stop('n_electrons is odd, please check', 'construct_first_order_dm_supercell_scalapack')
    endif

    i_spin = 1 ! to be extended

    !------real-(1).prepare C, C_occ, S^(1) ----------------------
       allocate(tmp_C(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_occ(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C_occ')
       allocate(tmp_S1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_S1')

       allocate(tmp_1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_1')
       allocate(tmp_2(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_2')
       allocate(tmp_C1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C1')

       first_order_ham_supercell_scalapack(:,:) = 0d0

       tmp_C(:,:) = eigenvec_supercell_scalapack(:,:,i_spin)
       tmp_S1(:,:) = - first_order_ovlp_supercell_scalapack(:,:)


    !------real-(2).DM^(1)_oo ----------------------
    !--------- = C*occ_number* (-C^+ S C) * C-------------

        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c

        tmp_2 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c

       tmp_C_occ = 0.0d0
       do i_state = 1, max_occ_number
          if(l_col_DFPT_phonon(i_state)>0) then
             tmp_C_occ(:,l_col_DFPT_phonon(i_state)) = tmp_C(:,l_col_DFPT_phonon(i_state)) * &
             2.0d0
          endif
       end do

       do i_state = max_occ_number+1, n_basis_sc_DFPT
          if(l_col_DFPT_phonon(i_state)>0) then
             tmp_C_occ(:,l_col_DFPT_phonon(i_state)) = 0.0d0
          endif
       end do

        tmp_1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_2, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c

        first_order_ham_supercell_scalapack = 0.0d0
        call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, first_order_ham_supercell_scalapack(:,:), 1, 1, sc_desc_DFPT_phonon)      ! beta,  c, ic, jc, desc_c


   !------real-(3).DM^(1)_(ov+vo) ----------------------
   !    = (C U) C^+*occ_number + C*occ_number (C U)^+
       tmp_C1 = 0.0d0
       call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                          first_order_U_supercell_scalapack, 1, 1, sc_desc_DFPT_phonon,   &  ! b, ib, jb, desc_b
                   0.0d0, tmp_C1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
       tmp_1 = 0.0d0
       call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C1,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                          tmp_C_occ, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                   0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
       first_order_ham_supercell_scalapack(:,:) = &
       first_order_ham_supercell_scalapack(:,:) + tmp_1(:,:)

       tmp_2 = 0.0d0
       call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C_occ,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                          tmp_C1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                   0.0d0, tmp_2,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c
       first_order_ham_supercell_scalapack(:,:) = &
       first_order_ham_supercell_scalapack(:,:) + tmp_2(:,:)


       deallocate(tmp_C)
       deallocate(tmp_C_occ)
       deallocate(tmp_S1)
       deallocate(tmp_1)
       deallocate(tmp_2)

  end subroutine construct_first_order_dm_supercell_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_sparse_matrix_from_supercell_scalapack
!  NAME
!    get_first_order_dm_sparse_matrix_from_supercell_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_sparse_matrix_from_supercell_scalapack( k_coord, k_atom, M_sparse )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use dimensions, only: n_centers_in_sc_DFPT
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    integer, intent(in) :: k_coord, k_atom
    real*8, intent(out) :: M_sparse(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size)
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, iuo, juo, io, jo, i_index, lr, lc
    integer :: j_cell_trans, jo_trans, i_cell_in_sc_DFPT, i_cell_trans, io_trans
    integer :: k_cell, k_center

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_first_order_dm_sparse_matrix_from_supercell_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if


    do iuo = 1, n_basis
    do i_cell = 1, n_cells_in_hamiltonian-1

       i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell)

       if (index_hamiltonian(1,i_cell, iuo) > 0) then
       do i_index = index_hamiltonian(1, i_cell, iuo), &
              &     index_hamiltonian(2, i_cell, iuo)

       !------------------(1) j_cell = 1 ----------------------------
       !                                                            !
       !      d M_sparse (u R1, v)           d M_full (u R1, v)     !
       !   -----------------------   +=   -----------------------   !
       !      d R_K                           d R_K                 !
       !                                                            !

        juo =  column_index_hamiltonian(i_index)
        io  =  cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

        lr = l_row_DFPT_phonon(io)
        lc = l_col_DFPT_phonon(juo)

        if(lr.ne.0.and.lc.ne.0) then
         M_sparse(k_coord, k_atom, i_index) = M_sparse(k_coord, k_atom, i_index) + &
             first_order_ham_supercell_scalapack(lr,lc)
        endif

       !------------------(2) j_cell > 1 ----------------------------
       !                                                            !
       !      d M_sparse (u R1, v)         d M_full (u R1+R2, v R2) !
       !   -----------------------   =   -----------------------    !
       !      d R_K-R2 [k_center]             d R_K [k_atom]        !
       !                                                            !
        do j_cell_trans = 2,n_cells_in_sc_DFPT                             ! j_cell_trans = R2

           jo_trans     = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)
           i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans) ! i_cell_trans = R1 + R2
           io_trans     = cell_and_basis_to_basis_sc_DFPT( i_cell_trans,iuo)

           lr = l_row_DFPT_phonon(io_trans)
           lc = l_col_DFPT_phonon(jo_trans)

           if(lr.ne.0.and.lc.ne.0) then
           k_cell       = cell_diff_sc_DFPT(1,j_cell_trans)                ! k_cell = 0 - R2
           k_center     = cell_and_atom_to_center_sc_DFPT( k_cell, k_atom)
           M_sparse(k_coord, k_center,i_index) = M_sparse(k_coord, k_center,i_index) + &
               first_order_ham_supercell_scalapack(lr,lc)
           endif

        enddo

       end do ! i_index
       end if

   end do ! i_cell
   end do ! i_basis

  end subroutine get_first_order_dm_sparse_matrix_from_supercell_scalapack

!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_edm_supercell_scalapack
!  NAME
!    construct_first_order_edm_supercell_scalapack
!  SYNOPSIS
  subroutine construct_first_order_edm_supercell_scalapack()
!  PURPOSE
!    Construct first_order density matrix supercell in ScaLAPACK
!    here C+ = C*T
!    EDM1(my_k_point) = C*occ_number* ( C^+H1C-(Ei+Ej)(C^+ S1 C) ) * C
!                      + (C U1) C^+*occ_number E + C*occ_number E (C U1)^+
!
!  USES
    use dimensions, only: n_basis_sc_DFPT
    use mpi_tasks
    use pbc_lists
    use physics, only : n_electrons
    implicit none
!  ARGUMENTS
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col
    integer :: i_spin

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_eig(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_C_occ_eig(:,:)
    real*8, allocatable :: tmp_S1(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_2(:,:)
    real*8, allocatable :: tmp_C1(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id_DFPT_phonon>=npcol_DFPT_phonon*nprow_DFPT_phonon) return

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix

    max_occ_number =  NINT(n_electrons * n_cells_in_sc_DFPT / 2.0d0)  ! n_basis_sc_DFPT / 2 ! todo : extend

    if(mod(nint(n_electrons),2).ne.0) then
      call aims_stop('n_electrons is odd, please check', 'construct_first_order_dm_supercell_scalapack')
    endif

    i_spin = 1 ! to be extended

       allocate(tmp_C(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_occ(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C_occ')
       allocate(tmp_C_eig(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C_eig')
       allocate(tmp_C_occ_eig(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C_occ_eig')

       allocate(tmp_S1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_S1')
       allocate(tmp_H1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_H1')

       allocate(tmp_1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_1')
       allocate(tmp_2(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_2')
       allocate(tmp_C1(mxld_DFPT_phonon, mxcol_DFPT_phonon),stat=info)
       call check_allocation(info, 'tmp_C1')

       first_order_edm_supercell_scalapack(:,:) = 0d0

       tmp_C(:,:) = eigenvec_supercell_scalapack(:,:,i_spin)
       tmp_S1(:,:) = first_order_ovlp_supercell_scalapack(:,:)
       tmp_H1(:,:) =  first_order_ham_supercell_scalapack(:,:)

       do i_state = 1, n_basis_sc_DFPT
          if(l_col_DFPT_phonon(i_state)>0) then
             tmp_C_eig(:,l_col_DFPT_phonon(i_state)) = tmp_C(:,l_col_DFPT_phonon(i_state)) * &
              eigenvalues_supercell_scalapack(i_state, i_spin)
          endif
       end do

    !------real-(1).EDM^(1) with E1 ----------------------
    !--------- = C*occ_number* ( C^+H1C-(C^+ S1 C)E ) * C-------------
    !---------here E1_{pp} = ( C^+ H1 C- C^+ S1 C E )_{pp}

        !-----real C^+ H1 C-----------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c
        tmp_H1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_H1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
        tmp_2(:,:) = tmp_H1(:,:)


        !-----real C^+ S1 C E-----------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c
        tmp_S1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_eig, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_S1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
        tmp_2(:,:) = tmp_2(:,:) - tmp_S1


       !-------C*occ_number* ( C^+H1C-(C^+ S1 C)E ) * C----------------
       do i_state = 1, max_occ_number
          if(l_col_DFPT_phonon(i_state)>0) then
             tmp_C_occ(:,l_col_DFPT_phonon(i_state)) = tmp_C(:,l_col_DFPT_phonon(i_state)) * &
             2.0d0
             tmp_C_occ_eig(:,l_col_DFPT_phonon(i_state)) = tmp_C_eig(:,l_col_DFPT_phonon(i_state)) * &
             2.0d0
          endif
       enddo
       do i_state = max_occ_number+1, n_basis_sc_DFPT
          if(l_col_DFPT_phonon(i_state)>0) then
             tmp_C_occ(:,l_col_DFPT_phonon(i_state)) = 0.0d0
             tmp_C_occ_eig(:,l_col_DFPT_phonon(i_state)) = 0.0d0
          endif
       enddo


        tmp_1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_2, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c

        first_order_edm_supercell_scalapack = 0.0d0
        call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, first_order_edm_supercell_scalapack(:,:), 1, 1, sc_desc_DFPT_phonon) ! beta,  c, ic, jc, desc_c

     !------real-(2) EDM^(1) with (E C1 C^+ + E C C1^+)_oo   ----------------------
     !-------------= C*occ_number*E  ( -(C^+ S1 C) ) * C-----------
        tmp_S1(:,:) = first_order_ovlp_supercell_scalapack(:,:)

        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c

        tmp_S1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                     1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                            tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                     0.0d0, tmp_S1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c

        tmp_1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C_occ_eig, 1, 1, sc_desc_DFPT_phonon, &  ! alpha, a, ia, ja, desc_a
                          tmp_S1, 1, 1, sc_desc_DFPT_phonon,     &  !        b, ib, jb, desc_b
                   0.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c

        tmp_2 = 0.0d0
        call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                  1.0d0, tmp_1, 1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                         tmp_C, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                  0.0d0, tmp_2, 1, 1, sc_desc_DFPT_phonon)     ! beta,  c, ic, jc, desc_c
        first_order_edm_supercell_scalapack(:,:) = first_order_edm_supercell_scalapack(:,:) &
                                                   - tmp_2(:,:)


    !------real-(3).EDM^(1)_(ov+vo) ----------------------
    !    =  (C U1) C^+*occ_number E + C*occ_number E (C U1)^+
        tmp_C1 = 0.0d0
        call pdgemm("N","N",n_basis_sc_DFPT, n_basis_sc_DFPT, n_basis_sc_DFPT, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
              first_order_U_supercell_scalapack, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
        tmp_1 = 0.0d0
        call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ_eig, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc_DFPT_phonon)        ! beta,  c, ic, jc, desc_c
        first_order_edm_supercell_scalapack(:,:) = first_order_edm_supercell_scalapack(:,:) + tmp_1(:,:)

        tmp_2 = 0.0d0
        call pdgemm("N","T",n_basis_sc_DFPT, n_basis_sc_DFPT, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ_eig,  1, 1, sc_desc_DFPT_phonon,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1, 1, 1, sc_desc_DFPT_phonon,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2,  1, 1, sc_desc_DFPT_phonon)         ! beta,  c, ic, jc, desc_c
        first_order_edm_supercell_scalapack(:,:) = first_order_edm_supercell_scalapack(:,:) + tmp_2(:,:)


       deallocate(tmp_C)
       deallocate(tmp_C_occ)
       deallocate(tmp_C_eig)
       deallocate(tmp_C_occ_eig)
       deallocate(tmp_S1)
       deallocate(tmp_H1)
       deallocate(tmp_1)
       deallocate(tmp_2)

  end subroutine construct_first_order_edm_supercell_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_edm_sparse_matrix_from_supercell_scalapack
!  NAME
!    get_first_order_edm_sparse_matrix_from_supercell_scalapack
!  SYNOPSIS
   subroutine get_first_order_edm_sparse_matrix_from_supercell_scalapack( k_coord, k_atom, M_sparse )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use dimensions, only: n_centers_in_sc_DFPT
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    integer, intent(in) :: k_coord, k_atom
    real*8,intent(out) :: M_sparse(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size)
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, iuo, juo, io, jo, i_index, lr, lc
    integer :: j_cell_trans, jo_trans, i_cell_trans, i_cell_in_sc_DFPT, io_trans
    integer :: k_cell, k_center


    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_first_order_dm_sparse_matrix_from_supercell_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if


    do iuo = 1, n_basis
    do i_cell = 1, n_cells_in_hamiltonian-1

       i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell)

       if (index_hamiltonian(1,i_cell, iuo) > 0) then
       do i_index = index_hamiltonian(1, i_cell, iuo), &
              &     index_hamiltonian(2, i_cell, iuo)

       !------------------(1) j_cell = 1 ----------------------------
       !                                                            !
       !      d M_sparse (u R1, v)           d M_full (u R1, v)     !
       !   -----------------------   +=   -----------------------   !
       !      d R_K                           d R_K                 !
       !                                                            !

        juo =  column_index_hamiltonian(i_index)
        io  =  cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

        lr = l_row_DFPT_phonon(io)
        lc = l_col_DFPT_phonon(juo)

        if(lr.ne.0.and.lc.ne.0) then
         M_sparse(k_coord, k_atom, i_index) = M_sparse(k_coord, k_atom, i_index) + &
             first_order_edm_supercell_scalapack(lr,lc)
        endif

       !------------------(2) j_cell > 1 ----------------------------
       !                                                            !
       !      d M_sparse (u R1, v)         d M_full (u R1+R2, v R2) !
       !   -----------------------   =   -----------------------    !
       !      d R_K-R2 [k_center]             d R_K [k_atom]        !
       !                                                            !
        do j_cell_trans = 2,n_cells_in_sc_DFPT                             ! j_cell_trans = R2

           jo_trans     = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)
           i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans) ! i_cell_trans = R1 + R2
           io_trans     = cell_and_basis_to_basis_sc_DFPT( i_cell_trans,iuo)

           lr = l_row_DFPT_phonon(io_trans)
           lc = l_col_DFPT_phonon(jo_trans)

           if(lr.ne.0.and.lc.ne.0) then
           k_cell       = cell_diff_sc_DFPT(1,j_cell_trans)                ! k_cell = 0 - R2
           k_center     = cell_and_atom_to_center_sc_DFPT( k_cell, k_atom)
           M_sparse(k_coord, k_center,i_index) = M_sparse(k_coord, k_center,i_index) + &
               first_order_edm_supercell_scalapack(lr,lc)
           endif

        enddo

       end do ! i_index
       end if

   end do ! i_cell
   end do ! i_basis


  end subroutine get_first_order_edm_sparse_matrix_from_supercell_scalapack
!=========================================================================================
!=========================end for scalapack used in DFPT_phonon=========================
!=========================================================================================


!=========================================================================================
!=========================begin for scalapack used in DFPT_phonon_reduced_memory=========
!=========================================================================================
!----------------lapack version-------------------------  |---------scalapack version-------------------
!(1) construct_first_order_matrix_phonon_reduce_memory    |  construct_first_order_overlap_scalapack
!(2) construct_first_order_matrix_phonon_reduce_memory    |  construct_first_order_hamiltonian_scalapack
!(3) evaluate_first_order_DM_phonon_reduce_memory         |  construct_first_order_dm_scalapack
!(4) evaluate_first_order_DM_phonon_reduce_memory         |  get_first_order_dm_complex_sparse_matrix_scalapack
!(5) evaluate_first_order_U_phonon_reduce_memory          |  evaluate_first_order_U_scalapack
!(6) evaluate_first_order_EDM_phonon_reduce_memory        |  construct_first_order_edm_scalapack
!(7) evaluate_first_order_EDM_phonon_reduce_memory        |  get_first_order_edm_complex_sparse_matrix_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_overlap_scalapack
!  NAME
!    construct_first_order_overlap_scalapack
!  SYNOPSIS
  subroutine construct_first_order_overlap_scalapack( first_order_S_sparse )
!  PURPOSE
!    Sets the overlap matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use runtime_choices, only: use_local_index

    implicit none
!  ARGUMENTS
    complex*16, dimension(n_hamiltonian_matrix_size), intent(in) :: first_order_S_sparse
!  INPUTS
!    o first_order_S_sparse -- the first order overlap matrix at i_q_point, i_atom,i_coord
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_cell, i_col, i_row, lr, lc, idx


    ! construct_overlap_scalapack must not be used when use_local_index is set!
    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_first_order_overlap_scalapack + use_local_index")


    !---------we put the inintial of U1 here, because, every pertubation need one S1 calcualtion, before DM1.
    ! the lapack verion, the U1 inital is put before DM1 calcualtion.

    if(real_eigenvectors)then
       first_order_ovlp_scalapack = 0.
       first_order_U_scalapack = 0.0d0
    else
       first_order_ovlp_complex_scalapack = (0.d0,0.d0)
       first_order_U_complex_scalapack = 0.0d0
    end if

    ! Attention: Only the lower half of the matrices is set here for S1

    select case(packed_matrix_format)

    case(PM_index) !------------------------------------------------

       do i_cell = 1, n_cells_in_hamiltonian-1
          do i_row = 1, n_basis

             lr = l_row(i_row) !
             if(lr==0) cycle   ! skip if not local

             if(index_hamiltonian(1,i_cell,i_row) > 0) then

                do idx = index_hamiltonian(1,i_cell,i_row),index_hamiltonian(2,i_cell,i_row)

                   i_col = column_index_hamiltonian(idx)
                   lc = l_col(i_col) !
                   if(lc==0) cycle   ! skip if not local

                   if(real_eigenvectors)then
                      first_order_ovlp_scalapack(lr,lc) = first_order_ovlp_scalapack(lr,lc) +  &
                      dble(k_phase(i_cell,my_k_point)) * dble(first_order_S_sparse(idx))
                   else
                      first_order_ovlp_complex_scalapack(lr,lc) = first_order_ovlp_complex_scalapack(lr,lc) + &
                      k_phase(i_cell,my_k_point) * first_order_S_sparse(idx)
                   end if

                end do

             end if
          end do
       end do ! i_cell

    case default !---------------------------------------------------------

       write(use_unit,*) 'Error: construct_first_order_overlap_scalapack does not support non-packed matrices.'
       call aims_stop

    end select ! packed_matrix_format

  end subroutine construct_first_order_overlap_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_overlap_sparse_matrix_scalapack
!  NAME
!    get_first_order_overlap_sparse_matrix_scalapack
!  SYNOPSIS
   subroutine get_first_order_overlap_sparse_matrix_scalapack( matrix_sparse)
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none
!  ARGUMENTS
    complex*16 :: matrix_sparse(n_hamiltonian_matrix_size)
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = 0.0d0

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) !
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   first_order_ovlp_scalapack(lr,lc)*dble(k_phase(i_cell,my_k_point))*k_weights(my_k_point)
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble(first_order_ovlp_complex_scalapack(lr,lc)*dconjg(k_phase(i_cell,my_k_point))) &
                      *k_weights(my_k_point)

             endif
          end do
       end do
    end do

    call sync_vector_complex(matrix_sparse, n_hamiltonian_matrix_size)

  end subroutine get_first_order_overlap_sparse_matrix_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_hamiltonian_scalapack
!  NAME
!    construct_first_order_hamiltonian_scalapack
!  SYNOPSIS
  subroutine construct_first_order_hamiltonian_scalapack( first_order_H_sparse )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    complex*16, dimension(n_hamiltonian_matrix_size, n_spin), intent(in) :: first_order_H_sparse
!  INPUTS
!    o first_order_H_sparse -- the first order Hamilton matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_scalapack + use_local_index")

    if(real_eigenvectors)then
       first_order_ham_scalapack(:,:,:) = 0.
    else
       first_order_ham_complex_scalapack(:,:,:) = 0.
    end if

    ! Attention: Only the lower half of the matrices is set here for H1

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_row = 1, n_basis
                lr = l_row(i_row) ! we do not exchange row,col as H0.
                if(lr==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_row) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_row),index_hamiltonian(2,i_cell,i_row)

                      i_col = column_index_hamiltonian(idx)
                      lc = l_col(i_col) ! we do not exchange row,col as H0.
                      if(lc==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         first_order_ham_scalapack(lr,lc,i_spin) = first_order_ham_scalapack (lr,lc,i_spin) &
                                + dble(k_phase(i_cell,my_k_point) * first_order_H_sparse(idx,i_spin))
                      else
                         first_order_ham_complex_scalapack (lr,lc,i_spin) =  &
                         first_order_ham_complex_scalapack (lr,lc,i_spin) &
                                + k_phase(i_cell,my_k_point) * first_order_H_sparse(idx,i_spin)
                      end if

                   end do

                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_hamiltonian_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format
    end do

  end subroutine construct_first_order_hamiltonian_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_hamiltonian_sparse_matrix_scalapack
!  NAME
!  SYNOPSIS
   subroutine get_first_order_hamiltonian_sparse_matrix_scalapack( matrix_sparse)
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none
!  ARGUMENTS
    complex*16 :: matrix_sparse(n_hamiltonian_matrix_size)
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = 0.0d0

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) !
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   first_order_ham_scalapack(lr,lc,1)*dble(k_phase(i_cell,my_k_point))*k_weights(my_k_point)
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble(first_order_ham_complex_scalapack(lr,lc,1)*dconjg(k_phase(i_cell,my_k_point))) &
                      *k_weights(my_k_point)

             endif
          end do
       end do
    end do

    call sync_vector_complex(matrix_sparse, n_hamiltonian_matrix_size)

  end subroutine get_first_order_hamiltonian_sparse_matrix_scalapack

!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_dm_scalapack
!  NAME
!    construct_first_order_dm_scalapack
!  SYNOPSIS
  subroutine construct_first_order_dm_scalapack(occ_numbers, i_spin)
!  PURPOSE
!    Construct first_order density matrix at my_k_point in ScaLAPACK
!    here C+ = C*T
!    DM1(my_k_point) = C*occ_number* (-C^+ S C) * C
!                      + (C U) C^+*occ_number + C*occ_number (C U)^+
!
!  USES
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    integer, intent(IN) :: i_spin
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_S1(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_2(:,:)
    real*8, allocatable :: tmp_C1(:,:)

    complex*16, allocatable :: tmp_C_complex(:,:)
    complex*16, allocatable :: tmp_C_occ_complex(:,:)
    complex*16, allocatable :: tmp_S1_complex(:,:)
    complex*16, allocatable :: tmp_1_complex(:,:)
    complex*16, allocatable :: tmp_2_complex(:,:)
    complex*16, allocatable :: tmp_C1_complex(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix
    max_occ_number = 0

    !Note:  Here occ_number have already contained with k_weights, but do not have in lapack version.
    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
          max_occ_number = i_state
       endif
    enddo

    if (real_eigenvectors) then

    !------real-(1).prepare C, C_occ, S^(1) ----------------------
       allocate(tmp_C(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_occ(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ')
       allocate(tmp_S1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_S1')

       allocate(tmp_1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1')
       allocate(tmp_2(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_2')
       allocate(tmp_C1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C1')

       first_order_ham_scalapack(:,:,i_spin) = 0d0

       tmp_C(:,:) = eigenvec(:,:,i_spin)
       tmp_S1(:,:) = - first_order_ovlp_scalapack(:,:)
       call set_full_matrix_real_L_to_U(tmp_S1)


    !------real-(2).DM^(1)_oo ----------------------
    !--------- = C*occ_number* (-C^+ S C) * C-------------

     !    call pdsyrk('U', 'N', n_basis, max_occ_number, !uplo, trans, n, k
     !                1.0d0, tmp, 1, 1, sc_desc, &       !alpha, a, ia, ja, desc_a
     !                0.0d0, ham(1,1,i_spin), 1, 1, sc_desc ) !  beta, c, ic, jc, desc_c

        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        tmp_2 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

       do i_state = 1, n_states
          if (occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
             if(l_col(i_state)>0) then
               tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
               occ_numbers(i_state, i_spin, my_k_point)
             endif
          elseif(l_col(i_state).ne.0) then
            tmp_C_occ(:,l_col(i_state)) = 0d0
          endif
       end do

        tmp_1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_2, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1, 1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        first_order_ham_scalapack = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, first_order_ham_scalapack(:,:,i_spin), 1, 1, sc_desc)      ! beta,  c, ic, jc, desc_c

    !------real-(3).DM^(1)_(ov+vo) ----------------------
    !    = (C U) C^+*occ_number + C*occ_number (C U)^+
        tmp_C1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           first_order_U_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_1 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        first_order_ham_scalapack(:,:,i_spin) = first_order_ham_scalapack(:,:,i_spin) + tmp_1(:,:)

        tmp_2 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        first_order_ham_scalapack(:,:,i_spin) = first_order_ham_scalapack(:,:,i_spin) + tmp_2(:,:)


       deallocate(tmp_C)
       deallocate(tmp_C_occ)
       deallocate(tmp_S1)
       deallocate(tmp_1)
       deallocate(tmp_2)

    else

    !------complex-(1).prepare C, C_occ, S^(1) ----------------------
       allocate(tmp_C_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       allocate(tmp_C_occ_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ_complex')
       allocate(tmp_S1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_S1_complex')

       allocate(tmp_1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1_complex')
       allocate(tmp_2_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_2_complex')
       allocate(tmp_C1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C1_complex')

       first_order_ham_complex_scalapack(:,:,i_spin) = 0d0

       ! Remeber that C_aims = C^*, so first we need to make C=C_aims^*
       tmp_C_complex(:,:) = dconjg(eigenvec_complex(:,:,i_spin))

       tmp_S1_complex(:,:) = - first_order_ovlp_complex_scalapack(:,:)
       call set_full_matrix_complex_L_to_U(tmp_S1_complex)

    !------complex-(2).DM^(1)_oo ----------------------
    !--------- = C*occ_number* (-C^+ S C) * C-------------

     !    call pdsyrk('U', 'N', n_basis, max_occ_number, !uplo, trans, n, k
     !                1.0d0, tmp, 1, 1, sc_desc, &       !alpha, a, ia, ja, desc_a
     !                0.0d0, ham(1,1,i_spin), 1, 1, sc_desc ) !  beta, c, ic, jc, desc_c

        !--------(-C^+ S C)------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_2_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis,         &  ! transa, transb, m, n, k
                     1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                            tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                     0.0d0, tmp_2_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


       do i_state = 1, n_states
          if (occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
             if(l_col(i_state)>0) then
               tmp_C_occ_complex(:,l_col(i_state)) = tmp_C_complex(:,l_col(i_state)) * &
                                             occ_numbers(i_state, i_spin, my_k_point)
             endif
          elseif(l_col(i_state).ne.0) then
               tmp_C_occ_complex(:,l_col(i_state)) = 0d0
          endif
       end do

!-----------------begin a DM0 benchmark test---------------------------------------
!       tmp_2_complex = 0
!       call pzgemm("N","C",n_basis, n_basis, n_basis,         &  ! transa, transb, m, n, k
!                     1.0d0, tmp_C_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
!                            tmp_C_occ_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
!                     0.0d0, tmp_2_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
!        first_order_ham_complex_scalapack(:,:,i_spin) = tmp_2_complex(:,:) !!!
!-----------------end a DM0 benchmark test---------------------------------------



      !---------C*occ_number* (-C^+ S C) * C-----------
      tmp_1_complex = 0.0d0
      call pzgemm("N","N",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C_occ_complex, 1, 1, sc_desc, &  ! alpha, a, ia, ja, desc_a
                          tmp_2_complex, 1, 1, sc_desc,     &  !        b, ib, jb, desc_b
                   0.0d0, tmp_1_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

      first_order_ham_complex_scalapack = 0.0d0
      call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                  1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                         tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                  0.0d0, first_order_ham_complex_scalapack(:,:,i_spin), 1, 1, sc_desc)     ! beta,  c, ic, jc, desc_c

!-------------begin shanghui's debug tool for scalapack version---------------
        print *, "wyj debug DM real"
       do  i_local_row = 1, n_my_rows
       do  i_local_col = 1, n_my_cols
           i_global_row = my_row(i_local_row)
           i_global_col = my_col(i_local_col)
           if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.1) then
           !write(use_unit,*) 'C(C+S1C)C+ (11):', first_order_ham_complex_scalapack(i_local_row,i_local_col,i_spin)
           write(use_unit,*) 'C(C+S1C)C+ (11):', first_order_ham_scalapack(i_local_row,i_local_col,i_spin)
           endif
           if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.2) then
           !write(use_unit,*) 'C(C+S1C)C+ (12):', first_order_ham_complex_scalapack(i_local_row,i_local_col,i_spin)
           write(use_unit,*) 'C(C+S1C)C+ (12):', first_order_ham_scalapack(i_local_row,i_local_col,i_spin)
           endif
       enddo
       enddo
!-------------end shanghui's debug tool for scalapack version---------------

!!     if(my_k_point.eq.1) then
!!     do i_state = 1, n_states
!!     do j_state = 1, n_states
!!
!!        if(l_row(i_state).ne.0.and.l_col(j_state).ne.0) then
!!        write(use_unit,*) 'U', i_state, j_state, first_order_U_complex_scalapack(l_row(i_state),l_col(j_state))
!!        endif
!!
!!     enddo
!!     enddo
!!     endif

    !------complex-(3).DM^(1)_(ov+vo) ----------------------
    !    = (C U) C^+*occ_number + C*occ_number (C U)^+

        tmp_C1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           first_order_U_complex_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1_complex,  1, 1, sc_desc)       ! beta,  c, ic, jc, desc_c

        tmp_1_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)            ! beta,  c, ic, jc, desc_c
        first_order_ham_complex_scalapack(:,:,i_spin) = first_order_ham_complex_scalapack(:,:,i_spin) &
                                                      + tmp_1_complex(:,:)

        tmp_2_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2_complex,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        first_order_ham_complex_scalapack(:,:,i_spin) = first_order_ham_complex_scalapack(:,:,i_spin) &
                                                      + tmp_2_complex(:,:)
!-------------begin shanghui's debug tool for scalapack version---------------
        print *, "wyj debug DM complex"
      do  i_local_row = 1, n_my_rows
      do  i_local_col = 1, n_my_cols
          i_global_row = my_row(i_local_row)
          i_global_col = my_col(i_local_col)
          if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.1) then
          write(use_unit,*) 'U1 (11):', first_order_U_complex_scalapack(i_local_row,i_local_col)
          write(use_unit,*) 'C1 C+ (11):', tmp_1_complex(i_local_row,i_local_col)
          endif
          if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.2) then
          write(use_unit,*) 'U1 (12):', first_order_U_complex_scalapack(i_local_row,i_local_col)
          write(use_unit,*) 'C1 C+ (12):', tmp_1_complex(i_local_row,i_local_col)
          endif
      enddo
      enddo
!-------------end shanghui's debug tool for scalapack version---------------



!!     if(my_k_point.eq.1) then
!!     do i_state = 1, n_states
!!     do j_state = 1, n_states
!!
!!        if(l_row(i_state).ne.0.and.l_col(j_state).ne.0) then
!!        write(use_unit,*) 'C^+UC', i_state, j_state,  &
!!                    tmp_1_complex(l_row(i_state),l_col(j_state))+tmp_2_complex(l_row(i_state),l_col(j_state))
!!        endif
!!
!!     enddo
!!     enddo
!!     endif


       deallocate(tmp_C_complex)
       deallocate(tmp_C_occ_complex)
       deallocate(tmp_S1_complex)
       deallocate(tmp_1_complex)
       deallocate(tmp_2_complex)


    endif

  end subroutine construct_first_order_dm_scalapack

  subroutine print_ham_cpscf(matrix_scalapack)
      ! use
      use dimensions
      use mpi_tasks, only: aims_stop, myid
      ! arguments
      real*8 matrix_scalapack(mxld, mxcol, n_spin)
      ! local variables
      integer :: i_local_row, i_local_col, i_global_row, i_global_col
      integer :: i_spin=1

      print *, 'myid=', myid, 'BLACS=', n_my_rows, 'X', n_my_cols
      do  i_local_row = 1, n_my_rows
      do  i_local_col = 1, n_my_cols
      i_global_row = my_row(i_local_row)
      i_global_col = my_col(i_local_col)
      !if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.1) then
      !    !write(use_unit,*) 'C(C+S1C)C+ (11):', first_order_ham_complex_scalapack(i_local_row,i_local_col,i_spin)
      !    !write(use_unit,*) 'wyj_ham: C(C+S1C)C+ (11):', first_order_ham_scalapack(i_local_row,i_local_col,i_spin)
      !    write(use_unit,*) 'wyj_matrix_scalapack: 1= ', matrix_scalapack(i_local_row,i_local_col,i_spin)
      !endif
      !if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.2) then
      !    !write(use_unit,*) 'C(C+S1C)C+ (12):', first_order_ham_complex_scalapack(i_local_row,i_local_col,i_spin)
      !    !write(use_unit,*) 'wyj_ham: C(C+S1C)C+ (12):', first_order_ham_scalapack(i_local_row,i_local_col,i_spin)
      !    write(use_unit,*) 'wyj_matrix_scalapack: 2= ', matrix_scalapack(i_local_row,i_local_col,i_spin)
      !endif

      write(use_unit,*), 'myid=', myid, 'wyj_matrix_scalapack: (', i_global_row, i_global_col, ')=', matrix_scalapack(i_local_row,i_local_col,i_spin)
      enddo
      enddo
      !-------------end shanghui's debug tool for scalapack version---------------
  end subroutine print_ham_cpscf

!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/evaluate_first_order_U_scalapack
!  NAME
!    evaluate_first_order_U_scalapack
!  SYNOPSIS
  subroutine evaluate_first_order_U_scalapack(occ_numbers, eigenvalues)
!  PURPOSE
!    evaluate first_order_U
!    here C^+ = C*T
!    U1_pq(my_k_point) =  (C^+ S1 C E - C^+ H1 C)_pq/(E_pp-E_qq)
!
!
!  USES
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
!  INPUTS
!    o occ_numbers -- the occupation numbers with k_weights
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, j_state, info, i_spin
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_eig(:,:)
    real*8, allocatable :: tmp_S1(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)

    complex*16, allocatable :: tmp_C_complex(:,:)
    complex*16, allocatable :: tmp_C_eig_complex(:,:)
    complex*16, allocatable :: tmp_S1_complex(:,:)
    complex*16, allocatable :: tmp_H1_complex(:,:)
    complex*16, allocatable :: tmp_1_complex(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    do i_spin = 1, n_spin

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix
    max_occ_number = 0

    !Note:  Here occ_number have already contained with k_weights, but do not have in lapack version.
    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
          max_occ_number = i_state
       endif
    enddo

    if (real_eigenvectors) then

    !------real-(1).prepare C, C_eig, S^(1), H^(1)----------------------
       allocate(tmp_C(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_eig(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_eig')
       allocate(tmp_S1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_S1')
       allocate(tmp_H1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_H1')
       allocate(tmp_1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1')

       tmp_C(:,:) = eigenvec(:,:,i_spin)
       tmp_S1(:,:) = first_order_ovlp_scalapack(:,:)
       call set_full_matrix_real_L_to_U(tmp_S1)

       tmp_H1(:,:) = first_order_ham_scalapack(:,:, i_spin)
       call set_full_matrix_real_L_to_U(tmp_H1)


       do i_state = 1, n_states
          if(l_col(i_state)>0) then
             tmp_C_eig(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
              eigenvalues(i_state, i_spin, my_k_point)
          endif
       end do

    !------real-(2) C^+ S1 C E ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        tmp_S1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_eig, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_S1, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------real-(3) C^+ H1 C ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1 , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_H1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_H1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------real-(4) U1_pq = (C^+ S1 C E - C^+ H1 C)_pq/(E_pp-E_qq) -----------------
        first_order_U_scalapack = 0.0d0

        do  i_local_row = 1, n_my_rows
        do  i_local_col = 1, n_my_cols

            i_global_row = my_row(i_local_row)
            i_global_col = my_col(i_local_col)

            if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

            first_order_U_scalapack(i_local_row,i_local_col) =  &
           ( tmp_S1(i_local_row,i_local_col) - tmp_H1(i_local_row,i_local_col)) / &
           ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

           endif

           if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

            first_order_U_scalapack(i_local_row,i_local_col) =  &
           ( tmp_S1(i_local_row,i_local_col) - tmp_H1(i_local_row,i_local_col)) / &
           ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

           endif

       enddo
       enddo

!
!       do i_state = 1, max_occ_number
!       do j_state = max_occ_number+1, n_states
!
!          if(l_row(i_state)>0.and.l_col(j_state)>0) then
!            first_order_U_scalapack(l_row(i_state),l_col(j_state)) = &
!          ( tmp_S1(l_row(i_state),l_col(j_state))-tmp_H1(l_row(i_state),l_col(j_state)) ) / &
!          ( eigenvalues(i_state, i_spin, my_k_point) - eigenvalues(j_state, i_spin, my_k_point))
!          endif
!
!          if(l_row(j_state)>0.and.l_col(i_state)>0) then
!            first_order_U_scalapack(l_row(j_state),l_col(i_state)) = &
!          ( tmp_S1(l_row(j_state),l_col(i_state))-tmp_H1(l_row(j_state),l_col(i_state)) ) / &
!          ( eigenvalues(j_state, i_spin, my_k_point) - eigenvalues(i_state, i_spin, my_k_point))
!          endif
!
!       enddo
!       enddo
!

       deallocate(tmp_C)
       deallocate(tmp_C_eig)
       deallocate(tmp_S1)
       deallocate(tmp_H1)
       deallocate(tmp_1)

    else

    !------complex-(1).prepare C, C_eig, S^(1), H^(1)----------------------
       allocate(tmp_C_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       allocate(tmp_C_eig_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_eig_complex')
       allocate(tmp_S1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_S1_complex')
       allocate(tmp_H1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_H1_complex')
       allocate(tmp_1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1_complex')

       tmp_C_complex(:,:) = dconjg(eigenvec_complex(:,:,i_spin))

       tmp_S1_complex(:,:) = first_order_ovlp_complex_scalapack(:,:)
       call set_full_matrix_complex_L_to_U(tmp_S1_complex)

       tmp_H1_complex(:,:) = first_order_ham_complex_scalapack(:,:,i_spin)
       call set_full_matrix_complex_L_to_U(tmp_H1_complex)


       do i_state = 1, n_states
          if(l_col(i_state)>0) then
             tmp_C_eig_complex(:,l_col(i_state)) = tmp_C_complex(:,l_col(i_state)) * &
              eigenvalues(i_state, i_spin, my_k_point)
          endif
       end do

    !------complex-(2) C^+ S1 C E ----------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        tmp_S1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_eig_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_S1_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------complex-(3) C^+ H1 C ----------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1_complex , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_H1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_H1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------real-(4) U1_pq = (C^+ S1 C E - C^+ H1 C)_pq/(E_pp-E_qq) -----------------
        first_order_U_complex_scalapack = 0.0d0

        do  i_local_row = 1, n_my_rows
        do  i_local_col = 1, n_my_cols

            i_global_row = my_row(i_local_row)
            i_global_col = my_col(i_local_col)

            if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

            first_order_U_complex_scalapack(i_local_row,i_local_col) =  &
           ( tmp_S1_complex(i_local_row,i_local_col) - tmp_H1_complex(i_local_row,i_local_col)) / &
           ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

!!               if(my_k_point.eq.1.and.i_global_row.eq.2.and.i_global_col.eq.1) then
!!               write(use_unit,*) '(S1)_{21}', tmp_S1_complex(i_local_row,i_local_col),tmp_S1_complex(i_local_col,i_local_row)
!!               write(use_unit,*) '(H1)_{21}', tmp_H1_complex(i_local_row,i_local_col),tmp_H1_complex(i_local_col,i_local_row)
!!               write(use_unit,*) 'E_2,E_1',eigenvalues(i_global_row, i_spin, my_k_point),eigenvalues(i_global_col, i_spin, my_k_point)
!!               write(use_unit,*) 'U1_{21}',first_order_U_complex_scalapack(i_local_row,i_local_col)
!!               endif

           endif

           if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

            first_order_U_complex_scalapack(i_local_row,i_local_col) =  &
           ( tmp_S1_complex(i_local_row,i_local_col) - tmp_H1_complex(i_local_row,i_local_col)) / &
           ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

           endif

       enddo
       enddo

!!!        do i_state = 1, max_occ_number
!!!        do j_state = max_occ_number+1, n_states
!!!
!!!         !!!  if(l_row(i_state)>0.and.l_col(j_state)>0) then
!!!         !!!    first_order_U_complex_scalapack(l_row(i_state),l_col(j_state)) = &
!!!         !!!  ( tmp_S1_complex(l_row(i_state),l_col(j_state))-tmp_H1_complex(l_row(i_state),l_col(j_state)) ) / &
!!!         !!!  ( eigenvalues(i_state, i_spin, my_k_point) - eigenvalues(j_state, i_spin, my_k_point))
!!!         !!!  endif
!!!
!!!           if(l_row(j_state)>0.and.l_col(i_state)>0) then
!!!             first_order_U_complex_scalapack(l_row(j_state),l_col(i_state)) = &
!!!           ( tmp_S1_complex(l_row(j_state),l_col(i_state))-tmp_H1_complex(l_row(j_state),l_col(i_state)) ) / &
!!!           ( eigenvalues(j_state, i_spin, my_k_point) - eigenvalues(i_state, i_spin, my_k_point))
!!!           endif
!!!
!!!        enddo
!!!        enddo


       deallocate(tmp_C_complex)
       deallocate(tmp_C_eig_complex)
       deallocate(tmp_S1_complex)
       deallocate(tmp_H1_complex)
       deallocate(tmp_1_complex)

    endif

   enddo ! i_spin

  end subroutine evaluate_first_order_U_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_complex_sparse_matrix_scalapack
!  NAME
!    get_first_order_dm_complex_sparse_matrix_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_complex_sparse_matrix_scalapack( matrix_sparse, i_spin )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    complex*16 :: matrix_sparse(n_hamiltonian_matrix_size)
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = (0.0d0,00d0)

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) !
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   first_order_ham_scalapack(lr,lc,i_spin)*dble(k_phase(i_cell,my_k_point))
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble( first_order_ham_complex_scalapack(lr,lc,i_spin)   &
                   * dconjg(k_phase(i_cell,my_k_point)) )


             endif
          end do
       end do
    end do

  end subroutine get_first_order_dm_complex_sparse_matrix_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_full_matrix_complex_L_to_U
!  NAME
!    set_full_matrix_complex
!  SYNOPSIS
  subroutine set_full_matrix_complex_L_to_U( mat )
!****s* scalapack_wrapper/set_full_matrix_complex_L_to_U
!  NAMEoutine set_full_matrix_complex_L_to_U( mat )
!  PURPOSE
!    Sets the Upper half of a distributed matrix from the Lower half
!  USES
    implicit none
!  ARGUMENTS
    complex*16, dimension(mxld, mxcol) :: mat
!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_col, i_row
    complex*16, allocatable :: tmp2(:,:)

    ! Allocate tmp2 bigger than necessary to catch overwrites in pdtran

    allocate(tmp2(mxld,mxcol+2*nb)) ! no idea whats really needed

    ! This routine is called only from the working set, so no need to check here
    call pztranc(n_basis,n_basis,(1.d0,0.d0),mat,1,1,sc_desc,(0.d0,0.d0),tmp2,1,1,sc_desc)

    do i_col=1,n_basis
       if(l_col(i_col)==0) cycle
       do i_row=1,i_col
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp2(l_row(i_row),l_col(i_col))
       enddo
    enddo

    ! For safety: Make diagonal real
    do i_col=1,n_basis
       if(l_col(i_col)==0 .or. l_row(i_col)==0) cycle
       mat(l_row(i_col),l_col(i_col)) = dble(mat(l_row(i_col),l_col(i_col)))
    enddo

    deallocate(tmp2)

  end  subroutine set_full_matrix_complex_L_to_U

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_full_matrix_real_L_to_U
!  NAME
!    set_full_matrix_complex
!  SYNOPSIS
  subroutine set_full_matrix_real_L_to_U( mat )
!****s* scalapack_wrapper/set_full_matrix_real_L_to_U
!  NAMEoutine set_full_matrix_real_L_to_U( mat )
!  PURPOSE
!    Sets the Upper half of a distributed matrix from the Lower half
!  USES
    implicit none
!  ARGUMENTS
    real*8, dimension(mxld, mxcol) :: mat
!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_col, i_row
    real*8, allocatable :: tmp2(:,:)

    ! Allocate tmp2 bigger than necessary to catch overwrites in pdtran

    allocate(tmp2(mxld,mxcol+2*nb)) ! no idea whats really needed

    ! This routine is called only from the working set, so no need to check here
    call pdtran(n_basis,n_basis,1.0d0,mat,1,1,sc_desc,0.d0,tmp2,1,1,sc_desc)

    do i_col=1,n_basis
       if(l_col(i_col)==0) cycle
       do i_row=1,i_col
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp2(l_row(i_row),l_col(i_col))
       enddo
    enddo

    deallocate(tmp2)

  end  subroutine set_full_matrix_real_L_to_U


!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_edm_scalapack
!  NAME
!    construct_first_order_edm_scalapack
!  SYNOPSIS
  subroutine construct_first_order_edm_scalapack(occ_numbers, eigenvalues, i_spin)
!  PURPOSE
!    Construct first_order energy density matrix at my_k_point in ScaLAPACK
!    here C+ = C*T
!    EDM1(my_k_point) = C*occ_number* ( C^+H1C-(Ei+Ej)(C^+ S1 C) ) * C
!                      + (C U1) C^+*occ_number E + C*occ_number E (C U1)^+
!
!  USES
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
    integer, intent(IN) :: i_spin
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_C_eig(:,:)
    real*8, allocatable :: tmp_C_occ_eig(:,:)
    real*8, allocatable :: tmp_S1(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_2(:,:)
    real*8, allocatable :: tmp_C1(:,:)

    complex*16, allocatable :: tmp_C_complex(:,:)
    complex*16, allocatable :: tmp_C_occ_complex(:,:)
    complex*16, allocatable :: tmp_C_eig_complex(:,:)
    complex*16, allocatable :: tmp_C_occ_eig_complex(:,:)
    complex*16, allocatable :: tmp_C_occ_eig1_complex(:,:)
    complex*16, allocatable :: tmp_S1_complex(:,:)
    complex*16, allocatable :: tmp_H1_complex(:,:)
    complex*16, allocatable :: tmp_1_complex(:,:)
    complex*16, allocatable :: tmp_2_complex(:,:)
    complex*16, allocatable :: tmp_C1_complex(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix
    max_occ_number = 0

    !Note:  Here occ_number have already contained with k_weights, but do not have in lapack version.
    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
          max_occ_number = i_state
       endif
    enddo

    if (real_eigenvectors) then

       allocate(tmp_C(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_occ(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ')
       allocate(tmp_C_eig(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_eig')
       allocate(tmp_C_occ_eig(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ_eig')

       allocate(tmp_S1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_S1')
       allocate(tmp_H1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_H1')

       allocate(tmp_1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1')
       allocate(tmp_2(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_2')
       allocate(tmp_C1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C1')

       first_order_edm_scalapack(:,:) = 0d0

       tmp_C(:,:) = eigenvec(:,:,i_spin)

       tmp_S1(:,:) = first_order_ovlp_scalapack(:,:)
       call set_full_matrix_real_L_to_U(tmp_S1)

       tmp_H1(:,:) =  first_order_ham_scalapack(:,:,i_spin)
       call set_full_matrix_real_L_to_U(tmp_H1)

       do i_state = 1, n_states
          if(l_col(i_state)>0) then
             tmp_C_eig(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
              eigenvalues(i_state, i_spin, my_k_point)
          endif
       end do

    !------real-(1).EDM^(1) with E1 ----------------------
    !--------- = C*occ_number* ( C^+H1C-(C^+ S1 C)E ) * C-------------
    !---------here E1_{pp} = ( C^+ H1 C- C^+ S1 C E )_{pp}

        !-----real C^+ H1 C-----------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        tmp_H1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_H1, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_2(:,:) = tmp_H1(:,:)


        !-----real C^+ S1 C E-----------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        tmp_S1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_eig, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_S1, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_2(:,:) = tmp_2(:,:) - tmp_S1


       do i_state = 1, n_states
          if (occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
             if(l_col(i_state)>0) then
               tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
               occ_numbers(i_state, i_spin, my_k_point)
               tmp_C_occ_eig(:,l_col(i_state)) = tmp_C_eig(:,l_col(i_state)) * &
               occ_numbers(i_state, i_spin, my_k_point)
             endif
          elseif(l_col(i_state).ne.0) then
            tmp_C_occ(:,l_col(i_state)) = 0d0
            tmp_C_occ_eig(:,l_col(i_state)) = 0d0
          endif
       end do

        tmp_1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_2, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1, 1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        first_order_edm_scalapack = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, first_order_edm_scalapack(:,:), 1, 1, sc_desc)      ! beta,  c, ic, jc, desc_c

     !------real-(2) EDM^(1) with (E C1 C^+ + E C C1^+)_oo   ----------------------
     !-------------= C*occ_number*E  ( -(C^+ S1 C) ) * C-----------
        tmp_S1(:,:) = first_order_ovlp_scalapack(:,:)
        call set_full_matrix_real_L_to_U(tmp_S1)

        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_S1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis,         &  ! transa, transb, m, n, k
                     1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                            tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                     0.0d0, tmp_S1, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C_occ_eig_complex, 1, 1, sc_desc, &  ! alpha, a, ia, ja, desc_a
                          tmp_S1, 1, 1, sc_desc,     &  !        b, ib, jb, desc_b
                   0.0d0, tmp_1, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_2 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                  1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                         tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                  0.0d0, tmp_2, 1, 1, sc_desc)     ! beta,  c, ic, jc, desc_c
        first_order_edm_scalapack(:,:) = first_order_edm_scalapack(:,:) &
                                                   - tmp_2(:,:)


    !------real-(3).EDM^(1)_(ov+vo) ----------------------
    !    =  (C U1) C^+*occ_number E + C*occ_number E (C U1)^+
        tmp_C1 = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           first_order_U_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_1 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ_eig, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        first_order_edm_scalapack(:,:) = first_order_edm_scalapack(:,:) + tmp_1(:,:)

        tmp_2 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ_eig,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        first_order_edm_scalapack(:,:) = first_order_edm_scalapack(:,:) + tmp_2(:,:)


       deallocate(tmp_C)
       deallocate(tmp_C_occ)
       deallocate(tmp_C_eig)
       deallocate(tmp_C_occ_eig)
       deallocate(tmp_S1)
       deallocate(tmp_H1)
       deallocate(tmp_1)
       deallocate(tmp_2)

    else

       allocate(tmp_C_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       allocate(tmp_C_occ_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ_complex')
       allocate(tmp_C_eig_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_eig_complex')
       allocate(tmp_C_occ_eig_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ_eig_complex')
       allocate(tmp_C_occ_eig1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ_eig_complex')

       allocate(tmp_S1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_S1_complex')
       allocate(tmp_H1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_H1_complex')

       allocate(tmp_1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1_complex')
       allocate(tmp_2_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_2_complex')
       allocate(tmp_C1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C1_complex')

       first_order_edm_complex_scalapack(:,:) = 0d0


       tmp_C_complex(:,:) = dconjg(eigenvec_complex(:,:,i_spin))

       tmp_S1_complex(:,:) = first_order_ovlp_complex_scalapack(:,:)
       call set_full_matrix_complex_L_to_U(tmp_S1_complex)

       tmp_H1_complex(:,:) = first_order_ham_complex_scalapack(:,:,i_spin)
       call set_full_matrix_complex_L_to_U(tmp_H1_complex)


       do i_state = 1, n_states
          if(l_col(i_state)>0) then
             tmp_C_eig_complex(:,l_col(i_state)) = tmp_C_complex(:,l_col(i_state)) * &
              eigenvalues(i_state, i_spin, my_k_point)
          endif
       end do




    !------complex-(1).EDM^(1) with E1 ----------------------
    !--------- = C*occ_number* ( C^+H1C-(C^+ S1 C)E ) * C-------------
    !---------here E1_{pp} = ( C^+ H1 C- C^+ S1 C E )_{pp}

       !------complex C^+ H1 C  ----------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_H1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis,         &  ! transa, transb, m, n, k
                     1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                            tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                     0.0d0, tmp_H1_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_2_complex(:,:) = tmp_H1_complex(:,:)

       !------complex C^+ S1 C E  ----------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_S1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis,         &  ! transa, transb, m, n, k
                     1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                            tmp_C_eig_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                     0.0d0, tmp_S1_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_2_complex(:,:) = tmp_2_complex(:,:) - tmp_S1_complex(:,:)


        do i_state = 1, n_states
          if (occ_numbers(i_state, i_spin, my_k_point) > 0d0) then
             if(l_col(i_state)>0) then
               tmp_C_occ_complex(:,l_col(i_state)) = tmp_C_complex(:,l_col(i_state)) * &
                                             occ_numbers(i_state, i_spin, my_k_point)

               tmp_C_occ_eig_complex(:,l_col(i_state)) = tmp_C_eig_complex(:,l_col(i_state)) * &
                                             occ_numbers(i_state, i_spin, my_k_point)

               tmp_C_occ_eig1_complex(:,l_col(i_state)) = tmp_C_complex(:,l_col(i_state)) * &
                  tmp_2_complex(l_col(i_state),l_col(i_state)) * occ_numbers(i_state, i_spin, my_k_point)
             endif
          elseif(l_col(i_state).ne.0) then
               tmp_C_occ_complex(:,l_col(i_state)) = 0d0
               tmp_C_occ_eig_complex(:,l_col(i_state)) = 0d0
               tmp_C_occ_eig1_complex(:,l_col(i_state)) = 0d0
          endif
        end do

!       first_order_edm_complex_scalapack = 0.0d0
!       call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
!                  1.0d0, tmp_C_occ_eig1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
!                         tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
!                  0.0d0, first_order_edm_complex_scalapack(:,:), 1, 1, sc_desc)     ! beta,  c, ic, jc, desc_c
       tmp_1_complex = 0.0d0
       call pzgemm("N","N",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                  1.0d0, tmp_C_occ_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                         tmp_2_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                  0.0d0, tmp_1_complex, 1, 1, sc_desc)     ! beta,  c, ic, jc, desc_c

       first_order_edm_complex_scalapack = 0.0d0
       call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                  1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                         tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                  0.0d0, first_order_edm_complex_scalapack(:,:), 1, 1, sc_desc)     ! beta,  c, ic, jc, desc_c


     !------complex-(2) EDM^(1) with (E C1 C^+ + E C C1^+)_oo   ----------------------
     !---------C*occ_number*E  ( -(C^+ S1 C) ) * C-----------
        tmp_S1_complex(:,:) = first_order_ovlp_complex_scalapack(:,:)
        call set_full_matrix_complex_L_to_U(tmp_S1_complex)

        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_S1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_S1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis,         &  ! transa, transb, m, n, k
                     1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                            tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                     0.0d0, tmp_S1_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                   1.0d0, tmp_C_occ_eig_complex, 1, 1, sc_desc, &  ! alpha, a, ia, ja, desc_a
                          tmp_S1_complex, 1, 1, sc_desc,     &  !        b, ib, jb, desc_b
                   0.0d0, tmp_1_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_2_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                  1.0d0, tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                         tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                  0.0d0, tmp_2_complex, 1, 1, sc_desc)     ! beta,  c, ic, jc, desc_c
        first_order_edm_complex_scalapack(:,:) = first_order_edm_complex_scalapack(:,:) &
                                                      - tmp_2_complex(:,:)


     !------complex-(3) EDM^(1) with (E C1 C^+  +  E C C1^+)_{ov+vo}   ----------------------
     !    = (C U) C^+*occ_number E + C*occ_number E (C U)^+

        tmp_C1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           first_order_U_complex_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1_complex,  1, 1, sc_desc)       ! beta,  c, ic, jc, desc_c

        tmp_1_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ_eig_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1_complex,  1, 1, sc_desc)            ! beta,  c, ic, jc, desc_c
        first_order_edm_complex_scalapack(:,:) = first_order_edm_complex_scalapack(:,:) &
                                                      + tmp_1_complex(:,:)

        tmp_2_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ_eig_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2_complex,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        first_order_edm_complex_scalapack(:,:) = first_order_edm_complex_scalapack(:,:) &
                                                      + tmp_2_complex(:,:)



       deallocate(tmp_C_complex)
       deallocate(tmp_C_occ_complex)
       deallocate(tmp_C_eig_complex)
       deallocate(tmp_C_occ_eig_complex)
       deallocate(tmp_C_occ_eig1_complex)
       deallocate(tmp_S1_complex)
       deallocate(tmp_H1_complex)
       deallocate(tmp_1_complex)
       deallocate(tmp_2_complex)


    endif

  end subroutine construct_first_order_edm_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_edm_complex_sparse_matrix_scalapack
!  NAME
!    get_first_order_edm_complex_sparse_matrix_scalapack
!  SYNOPSIS
   subroutine get_first_order_edm_complex_sparse_matrix_scalapack( matrix_sparse, i_spin )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    complex*16 :: matrix_sparse(n_hamiltonian_matrix_size)
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = (0.0d0,00d0)

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) !
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   first_order_edm_scalapack(lr,lc)*dble(k_phase(i_cell,my_k_point))
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble( first_order_edm_complex_scalapack(lr,lc)   &
                   * dconjg(k_phase(i_cell,my_k_point)) )


             endif
          end do
       end do
    end do

  end subroutine get_first_order_edm_complex_sparse_matrix_scalapack
!=========================================================================================
!=========================end for scalapack used in DFPT_phonon_reduced_memory===========
!=========================================================================================


!=========================================================================================
!=========================begin scalapack part used in DFPT_polarizability=========
!=========================================================================================

! Nath: For the moment this part of the code uses the first-order density matrix as a dense matrix,
! and not as a sparse matrix. For big molecules and/or heavy atoms, it may be necessary to build
! a sparse matrix.

!----------------lapack version-------------------------  |---------scalapack version-------------------

!(1) evaluate_first_order_DM_polarizability               |(1)  construct_first_order_dm_polar_scalapack
!    evaluate_first_order_DM_polarizability               |(2)  get_first_order_dm_polar_scalapack
!    evaluate_first_order_U_polarizability                |(3)  construct_first_order_ham_polar_scalapack
!(2) evaluate_first_order_U_polarizability                |(4)  evaluate_first_order_U_polar_scalapack


!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_dm_polar_scalapack
!  NAME
!    construct_first_order_dm_polar_scalapack
!  SYNOPSIS
  subroutine construct_first_order_dm_polar_scalapack(occ_numbers, i_spin)
!  PURPOSE
!    Construct first-order density matrix in ScaLAPACK
!    DM1 = (C U) C^T*occ_number + C*occ_number (C U)^T
!
!  USES
    ! n_states_k contains the actual number of states. Careful: the variable
    ! n_states (no _k) in Scalapack is given the same value as n_states for
    ! communication reasons. That's why there seems to be more eigenstates in
    ! the Scalapack version than in the lapack version, but these extra
    ! eigenvectors are artificial and should not be taken into account when
    ! summing over unoccupied states.
    use dimensions, only: n_states_k
    use mpi_tasks
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    integer, intent(IN) :: i_spin
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham_polar is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col
    integer :: i_coord, i_basis, j_basis

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_C1(:,:)
    real*8, allocatable :: tmp_U(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    ! We use first_order_ham_polar as a storage area for the first-order density matrix
    max_occ_number = 0

    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

    !------Prepare C, C_occ ----------------------
    allocate(tmp_C(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C')
    allocate(tmp_C_occ(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C_occ')
    allocate(tmp_1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_1')
    allocate(tmp_C1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C1')
    allocate(tmp_U(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_U')

    tmp_C(:,:) = eigenvec(:,:,i_spin)

    do i_state = 1, n_states
       if (occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          if(l_col(i_state)>0) then
            tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
            occ_numbers(i_state, i_spin, 1)
          endif
       elseif(l_col(i_state).ne.0) then
         tmp_C_occ(:,l_col(i_state)) = 0d0
       endif
    end do

    do i_coord = 1, 3

    !------ DM^(1)_(ov+vo) ----------------------
    !    = (C U) C^T*occ_number + C*occ_number (C U)^T

        tmp_C1 = 0.0d0
        first_order_ham_polar_scalapack(i_coord,:,:,i_spin) = 0.0d0
        tmp_U(:,:) = first_order_U_polar_scalapack(i_coord,:,:,i_spin)

        call pdgemm("N","N",n_basis, n_states_k, n_states_k, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_U, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_1 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        first_order_ham_polar_scalapack(i_coord,:,:,i_spin) = first_order_ham_polar_scalapack(i_coord,:,:,i_spin) + tmp_1(:,:)

        tmp_1 = 0.0d0

        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        first_order_ham_polar_scalapack(i_coord,:,:,i_spin) = first_order_ham_polar_scalapack(i_coord,:,:,i_spin) + tmp_1(:,:)

    end do ! i_coord

    deallocate(tmp_C)
    deallocate(tmp_C_occ)
    deallocate(tmp_1)
    deallocate(tmp_U)
    deallocate(tmp_C1)


  end subroutine construct_first_order_dm_polar_scalapack



!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_polar_scalapack
!  NAME
!    get_first_order_dm_polar_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_polar_scalapack( first_order_density_matrix_polar, i_spin )
!  PURPOSE
!    Reconstructs the global matrix from ScaLAPACK
!  USES
    implicit none
!  ARGUMENTS
    real*8 :: first_order_density_matrix_polar(3,n_basis,n_basis,n_spin)
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o first_order_density_matrix_polar -- set to the contents of first_order_ham_polar_scalapack
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer ::  lr, lc, i_col, i_row, i_coord

    character*200 :: info_str

    first_order_density_matrix_polar(:,:,:,i_spin)= 0.0d0

      do i_coord = 1, 3
       do i_col = 1, n_basis
          lc = l_col(i_col) ! local column number
          if(lc>0) then
             do i_row = 1, n_basis
                lr = l_row(i_row) ! local row number
                if(lr>0) then
                   first_order_density_matrix_polar(i_coord,i_row,i_col,i_spin) = first_order_ham_polar_scalapack(i_coord,lr,lc,i_spin)
                endif
             enddo ! i_row
          endif
       end do ! i_col
      end do ! i_coord


  end subroutine get_first_order_dm_polar_scalapack


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_ham_polar_scalapack
!  NAME
!    construct_first_order_ham_polar_scalapack
!  SYNOPSIS
  subroutine construct_first_order_ham_polar_scalapack( first_order_H_dense )
!  PURPOSE
!    Constructs the local first-order hamiltonian (ScaLAPACK format) from the global one
!  USES
    implicit none
!  ARGUMENTS
    real*8:: first_order_H_dense (3, n_basis, n_basis, n_spin )
!  INPUTS
!    o first_order_H_dense -- the first-order hamiltonian as a 2-dimensional array (=matrix)
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: n_spin_max, i_spin, i_index, i_col, i_row, lr, lc, i_coord

    first_order_ham_polar_scalapack(:,:,:,:) = 0

    ! Note: Here we directly construct the full matrix, not just the lower or upper part, so there is no need
    !       to symmetrize it later.

    do i_spin = 1, n_spin
      do i_coord = 1, 3
       do i_row = 1, n_basis
          lr = l_row(i_row) ! local row number
          if(lr>0) then
             do i_col = 1, n_basis
                lc = l_col(i_col) ! local column number
                if(lc>0) then
                   first_order_ham_polar_scalapack(i_coord,lr,lc,i_spin) = first_order_H_dense(i_coord,i_row,i_col,i_spin)
                endif
             enddo ! i_col
          endif
       end do ! i_row
      end do ! i_coord
    end do ! i_spin

  end subroutine construct_first_order_ham_polar_scalapack

!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/evaluate_first_order_U_polar_scalapack
!  NAME
!    evaluate_first_order_U_polar_scalapack
!  SYNOPSIS
  subroutine evaluate_first_order_U_polar_scalapack(occ_numbers, eigenvalues)
!  PURPOSE
!    evaluate first_order_U
!    U1_pq =  (-C^T H1 C)_pq/(E_pp-E_qq)
!
!  USES
    use mpi_tasks
    ! n_states_k contains the actual number of states. Careful: the variable
    ! n_states (no _k) in Scalapack is given the same value as n_states for
    ! communication reasons. That's why there seems to be more eigenstates in
    ! the Scalapack version than in the lapack version, but these extra
    ! eigenvectors are artificial and should not be taken into account when
    ! summing over unoccupied states.
    use dimensions, only: n_states_k
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
!  INPUTS
!    o occ_numbers -- the occupation numbers
!  OUTPUT
!    o the array ham_polar is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, j_state, info, i_spin
    integer :: i_local_row, i_local_col, i_global_row, i_global_col
    integer :: i_coord

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)

! -- Debug begin--
!    external  PDLAPRNT
!    real*8 , allocatable ::  work(:)
! -- Debug end--


    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    do i_spin = 1, n_spin

    max_occ_number = 0

    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

 ! Nath: For the polarizability (molecules), we have no k-points, so only real eigenvectors

    !------real-(1).prepare C, H^(1)----------------------

    allocate(tmp_C(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C')
    allocate(tmp_H1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_H1')
    allocate(tmp_1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_1')

    tmp_C(:,:) = eigenvec(:,:,i_spin)

    do i_coord = 1, 3

       ! At this points first_order_ham_polar_scalapack contains the first-order
       ! hamiltonian in scalapack format
       tmp_H1(:,:) = first_order_ham_polar_scalapack(i_coord,:,:,i_spin)
       ! call set_full_matrix_real_L_to_U(tmp_H1) ! Nath: not needed anymore since the matrix is already full


    !------ C^T H1 C ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_states_k, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1 , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
       tmp_H1 = 0.0d0
       call pdgemm("N","N",n_states_k, n_states_k, n_basis, &  ! transa, transb, m, n, k
                   1.0d0, tmp_1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                          tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                   0.0d0, tmp_H1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------ U1_pq = (- C^T H1 C)_pq/(E_pp-E_qq) -----------------
    ! Nath: Actually we only need to get the virt-occ and/or occ-virt part of U

       first_order_U_polar_scalapack(i_coord,:,:,i_spin) = 0.0d0

       do i_local_row = 1, n_my_rows
       do i_local_col = 1, n_my_cols

          i_global_row = my_row(i_local_row)
          i_global_col = my_col(i_local_col)

          if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

          first_order_U_polar_scalapack(i_coord,i_local_row,i_local_col,i_spin) =  &
          ( -tmp_H1(i_local_row,i_local_col)) / &
          ( eigenvalues(i_global_row, i_spin, 1) - eigenvalues(i_global_col, i_spin, 1) )

          endif

          if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

          first_order_U_polar_scalapack(i_coord,i_local_row,i_local_col,i_spin) =  &
          ( -tmp_H1(i_local_row,i_local_col)) / &
          ( eigenvalues(i_global_row, i_spin, 1) - eigenvalues(i_global_col, i_spin, 1) )

          endif

       enddo
       enddo


    end do ! i_coord

    deallocate(tmp_C)
    deallocate(tmp_H1)
    deallocate(tmp_1)


    enddo ! i_spin


! -- Debug begin--
 !     allocate(work(1))
 !     CALL PDLAPRNT( n_basis, n_basis, first_order_U_polar_scalapack(1,:,:,1), &
 !                    1, 1, sc_desc, 0, 0, 'U1_polar', 6, work )
 !     deallocate(work)
! -- Debug end--


  end subroutine evaluate_first_order_U_polar_scalapack



!! ----Nath: Debug >> needs to be removed ----
!!******
!!-----------------------------------------------------------------------------------
!!****s* scalapack_wrapper/get_first_order_U_polar_scalapack
!!  NAME
!!    get_first_order_U_polar_scalapack
!!  SYNOPSIS
!   subroutine get_first_order_U_polar_scalapack( first_order_U_polar, i_spin )
!!  PURPOSE
!!    Gets a matrix with true dimensions from ScaLAPACK
!!  USES
!    !use pbc_lists
!    implicit none
!!  ARGUMENTS
!    real*8 :: first_order_U_polar(3,n_basis,n_basis,n_spin)
!    integer :: i_spin
!!  INPUTS
!!    o i_spin -- the spin channel
!!  OUTPUT
!!    o matrix_dense -- set to the contents of ham_polar
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!
!    integer ::  lr, lc, i_col, i_row, i_coord
!
!    character*200 :: info_str
!
!!    if (packed_matrix_format /= PM_index) then
!!       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
!!       call localorb_info(info_str, use_unit, '(A)')
!!       call aims_stop
!!    end if
!
!    first_order_U_polar(:,:,:,i_spin)= 0.0d0
!
!      do i_coord = 1, 3
!       do i_col = 1, n_basis
!          lc = l_col(i_col) ! local column number
!          if(lc>0) then
!             !do i_row = 1, i_col, 1
!             do i_row = 1, n_basis
!                lr = l_row(i_row) ! local row number
!                if(lr>0) then
!                   first_order_U_polar(i_coord,i_row,i_col,i_spin) = first_order_U_polar_scalapack(i_coord,lr,lc,i_spin)
!                endif
!             enddo
!          endif
!       end do
!      end do
!
!    call sync_vector(first_order_U_polar,3*n_basis*n_basis*n_spin)
!
!  end subroutine get_first_order_U_polar_scalapack


!=========================================================================================
!=========================end scalapack part used in DFPT_polarizability===========
!=========================================================================================




!=========================================================================================
!=========================begin scalapack part used in DFPT_polar_reduce_memory=========
!=========================================================================================

! shanghui: This part of the code uses the first-order density matrix as a sparse matrix. 
! It is useful for big molecules and/or heavy atoms. 

!----------------lapack version-------------------------  |---------scalapack version-------------------

!(1) evaluate_first_order_DM_polarizability               |(1)  construct_first_order_dm_polar_reduce_memory_scalapack
!    evaluate_first_order_DM_polarizability               |(2)  get_first_order_dm_polar_reduce_memory_scalapack
!    evaluate_first_order_U_polarizability                |(3)  construct_first_order_ham_polar_reduce_memory_scalapack
!(2) evaluate_first_order_U_polarizability                |(4)  evaluate_first_order_U_polar_reduce_memory_scalapack


!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_dm_polar_reduce_memory_scalapack
!  NAME
!    construct_first_order_dm_polar_reduce_memory_scalapack
!  SYNOPSIS
  subroutine construct_first_order_dm_polar_reduce_memory_scalapack(occ_numbers)
!  PURPOSE
!    Construct first-order density matrix in ScaLAPACK
!    DM1 = (C U) C^T*occ_number + C*occ_number (C U)^T
!
!  USES
    ! n_states_k contains the actual number of states. Careful: the variable
    ! n_states (no _k) in Scalapack is given the same value as n_states for
    ! communication reasons. That's why there seems to be more eigenstates in
    ! the Scalapack version than in the lapack version, but these extra
    ! eigenvectors are artificial and should not be taken into account when
    ! summing over unoccupied states.
    use dimensions, only: n_states_k
    use mpi_tasks
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham_polar is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col
    integer :: i_coord, i_basis, j_basis
    integer :: i_spin

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_C1(:,:)
    real*8, allocatable :: tmp_U(:,:)
    real*8 :: time1

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return


    do i_spin = 1, n_spin
    ! We use first_order_ham_polar as a storage area for the first-order density matrix
    max_occ_number = 0

    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

    !------Prepare C, C_occ ----------------------
    allocate(tmp_C(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C')
    allocate(tmp_C_occ(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C_occ')
    allocate(tmp_1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_1')
    allocate(tmp_C1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C1')
    allocate(tmp_U(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_U')

    tmp_C(:,:) = eigenvec(:,:,i_spin)

    ! wyj: TODO debug
    !tmp_C = 1.0d0

    do i_state = 1, n_states
       if (occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          if(l_col(i_state)>0) then
            tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
            occ_numbers(i_state, i_spin, 1)
          endif
       elseif(l_col(i_state).ne.0) then
         tmp_C_occ(:,l_col(i_state)) = 0d0
       endif
    end do

    !------ DM^(1)_(ov+vo) ----------------------
    !    = (C U) C^T*occ_number + C*occ_number (C U)^T

        tmp_C1 = 0.0d0
        first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin) = 0.0d0
        tmp_U(:,:) = first_order_U_polar_reduce_memory_scalapack(:,:,i_spin)


        if(myid .eq. 0) then
            print*, "pdgemm, c_dm1, myid=", myid, ", N, N, n_basis=", n_basis, ", n_states_k=", n_states_k, ", n_states_k=", n_states_k
        endif
        time1 = mpi_wtime()
        call pdgemm("N","N",n_basis, n_states_k, n_states_k, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_U, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
      !   call hippdgemm("N","N",n_basis, n_states_k, n_states_k, &  ! transa, transb, m, n, k
      !               1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
      !                      tmp_U, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
      !               0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_1 = 0.0d0
        if(myid .eq. 0) then
            print*, "pdgemm, c_dm1, myid=", myid, ", time=", mpi_wtime() - time1
            print*, "pdgemm, c_dm2, myid=", myid, ", N, T, n_basis=", n_basis, ", n_basis=", n_basis, ", max_occ_number=", max_occ_number
        endif
        time1 = mpi_wtime()
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
      !   call hippdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
      !               1.0d0, tmp_C1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
      !                      tmp_C_occ, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
      !               0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin) = & 
        first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin) + tmp_1(:,:)

        tmp_1 = 0.0d0
        if(myid .eq. 0) then
            print*, "pdgemm, c_dm2, myid=", myid, ", time=", mpi_wtime() - time1
            print*, "pdgemm, c_dm3, myid=", myid, ", N, T, n_basis=", n_basis, ", n_basis=", n_basis, ", max_occ_number=", max_occ_number
        endif
        time1 = mpi_wtime()
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
      !   call hippdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
      !               1.0d0, tmp_C_occ,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
      !                      tmp_C1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
      !               0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        if(myid .eq. 0) then
            print*, "pdgemm, c_dm3, myid=", myid, ", time=", mpi_wtime() - time1
        endif
        first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin) = &
        first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin) + tmp_1(:,:)

    enddo ! i_spin

    deallocate(tmp_C)
    deallocate(tmp_C_occ)
    deallocate(tmp_1)
    deallocate(tmp_U)
    deallocate(tmp_C1)


  end subroutine construct_first_order_dm_polar_reduce_memory_scalapack



!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_polar_reduce_memory_scalapack
!  NAME
!    get_first_order_dm_polar_reduece_memory_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_polar_reduce_memory_scalapack( first_order_density_matrix_sparse )
!  PURPOSE
!    Reconstructs the global matrix from ScaLAPACK : the real part of
!    get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
!  USES
   use pbc_lists, only:  position_in_hamiltonian, n_cells_in_hamiltonian,  &
                         column_index_hamiltonian,index_hamiltonian
    implicit none
!  ARGUMENTS
    real*8 :: first_order_density_matrix_sparse(n_hamiltonian_matrix_size,n_spin)
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o first_order_density_matrix_sparse -- set to the contents of first_order_ham_polar_scalapack
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_spin
    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    first_order_density_matrix_sparse(:,:)= 0.0d0
     
    do i_spin = 1, n_spin 

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) ! 
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

                first_order_density_matrix_sparse(i_index,i_spin) = first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin)

          end do
       end do
    end do

    enddo ! i_spin

  end subroutine get_first_order_dm_polar_reduce_memory_scalapack


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_ham_polar_reduce_memory_scalapack
!  NAME
!    construct_first_order_ham_polar_reduce_memory_scalapack
!  SYNOPSIS
  subroutine construct_first_order_ham_polar_reduce_memory_scalapack( first_order_H_sparse )
!  PURPOSE
!    Constructs the local first-order hamiltonian (ScaLAPACK format) from the global one
!  USES
   use pbc_lists, only:  n_cells_in_hamiltonian,  &
                         column_index_hamiltonian,index_hamiltonian
   use mpi_tasks, only:  aims_stop
   use localorb_io, only: use_unit

    implicit none
!  ARGUMENTS
    real*8:: first_order_H_sparse (n_hamiltonian_matrix_size, n_spin )
!  INPUTS
!    o first_order_H_sparse -- the first-order hamiltonian 
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: n_spin_max, i_spin, i_index, i_col, i_row, lr, lc, i_cell, idx


    first_order_ham_polar_reduce_memory_scalapack(:,:,:) = 0.

    ! Attention: Only the lower half of the matrices is set here for H1

    do i_spin = 1, n_spin
       select case(packed_matrix_format)
       case(PM_index) !------------------------------------------------
          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_row = 1, n_basis

                lr = l_row(i_row) ! we do not exchange row,col as H0.
                if(lr==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_row) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_row),index_hamiltonian(2,i_cell,i_row)

                      i_col = column_index_hamiltonian(idx)
                      lc = l_col(i_col) ! we do not exchange row,col as H0.
                      if(lc==0) cycle   ! skip if not local

                         first_order_ham_polar_reduce_memory_scalapack(lr,lc,i_spin) =  & 
                         first_order_H_sparse(idx,i_spin)

                   end do
                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------
          write(use_unit,*) 'Error: construct_first_order_hamiltonian_polar_reduece_memor_scalapack does not support non-packed matrices.'
          call aims_stop
       end select ! packed_matrix_format
    end do  !i_spin

  end subroutine construct_first_order_ham_polar_reduce_memory_scalapack

!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/evaluate_first_order_U_polar_reduce_memory_scalapack
!  NAME
!    evaluate_first_order_U_polar_reduce_memory_scalapack
!  SYNOPSIS
  subroutine evaluate_first_order_U_polar_reduce_memory_scalapack(occ_numbers, eigenvalues)
!  PURPOSE
!    evaluate first_order_U
!    U1_pq =  (-C^T H1 C)_pq/(E_pp-E_qq)
!
!  USES
    use mpi_tasks
    ! n_states_k contains the actual number of states. Careful: the variable
    ! n_states (no _k) in Scalapack is given the same value as n_states for
    ! communication reasons. That's why there seems to be more eigenstates in
    ! the Scalapack version than in the lapack version, but these extra
    ! eigenvectors are artificial and should not be taken into account when
    ! summing over unoccupied states.
    use dimensions, only: n_states_k
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
!  INPUTS
!    o occ_numbers -- the occupation numbers
!  OUTPUT
!    o the array ham_polar is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, j_state, info, i_spin
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)

! -- Debug begin--
!    external  PDLAPRNT
!    real*8 , allocatable ::  work(:)
! -- Debug end--


    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    do i_spin = 1, n_spin

    max_occ_number = 0

    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

 ! Nath: For the polarizability (molecules), we have no k-points, so only real eigenvectors

    !------real-(1).prepare C, H^(1)----------------------

    allocate(tmp_C(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C')
    allocate(tmp_H1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_H1')
    allocate(tmp_1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_1')

    tmp_C(:,:) = eigenvec(:,:,i_spin)


       ! At this points first_order_ham_polar_reduce_meory_scalapack contains the first-order
       ! hamiltonian in scalapack format
       tmp_H1(:,:) = first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin)
       call set_full_matrix_real_L_to_U(tmp_H1)


    !------ C^T H1 C ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_states_k, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1 , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
       tmp_H1 = 0.0d0
       call pdgemm("N","N",n_states_k, n_states_k, n_basis, &  ! transa, transb, m, n, k
                   1.0d0, tmp_1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                          tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                   0.0d0, tmp_H1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------ U1_pq = (- C^T H1 C)_pq/(E_pp-E_qq) -----------------
    ! Nath: Actually we only need to get the virt-occ and/or occ-virt part of U

       first_order_U_polar_reduce_memory_scalapack(:,:,i_spin) = 0.0d0

       do i_local_row = 1, n_my_rows
       do i_local_col = 1, n_my_cols

          i_global_row = my_row(i_local_row)
          i_global_col = my_col(i_local_col)

          if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

          first_order_U_polar_reduce_memory_scalapack(i_local_row,i_local_col,i_spin) =  &
          ( -tmp_H1(i_local_row,i_local_col)) / &
          ( eigenvalues(i_global_row, i_spin, 1) - eigenvalues(i_global_col, i_spin, 1) )

          endif

          if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

          first_order_U_polar_reduce_memory_scalapack(i_local_row,i_local_col,i_spin) =  &
          ( -tmp_H1(i_local_row,i_local_col)) / &
          ( eigenvalues(i_global_row, i_spin, 1) - eigenvalues(i_global_col, i_spin, 1) )

          endif

       enddo
       enddo



    deallocate(tmp_C)
    deallocate(tmp_H1)
    deallocate(tmp_1)


    enddo ! i_spin


  end subroutine evaluate_first_order_U_polar_reduce_memory_scalapack

!******
!-----------------------------------------------------------
!****scalapack_wrapper/evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf
!  NAME
!    evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf
!  SYNOPSIS
  subroutine evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf(occ_numbers, eigenvalues, which)
!  PURPOSE
!    evaluate first_order_U
!    U1_pq =  (-C^T H1 C)_pq/(E_pp-E_qq)
!
!  USES
    use mpi_tasks
    ! n_states_k contains the actual number of states. Careful: the variable
    ! n_states (no _k) in Scalapack is given the same value as n_states for
    ! communication reasons. That's why there seems to be more eigenstates in
    ! the Scalapack version than in the lapack version, but these extra
    ! eigenvectors are artificial and should not be taken into account when
    ! summing over unoccupied states.
    use dimensions, only: n_states_k
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
    integer :: which
!  INPUTS
!    o occ_numbers -- the occupation numbers
!  OUTPUT
!    o the array ham_polar is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, j_state, info, i_spin
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8 :: time1

! -- Debug begin--
!    external  PDLAPRNT
!    real*8 , allocatable ::  work(:)
! -- Debug end--


    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    do i_spin = 1, n_spin

    max_occ_number = 0

    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, 1) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

 ! Nath: For the polarizability (molecules), we have no k-points, so only real eigenvectors

    !------real-(1).prepare C, H^(1)----------------------

    allocate(tmp_C(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_C')
    allocate(tmp_H1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_H1')
    allocate(tmp_1(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp_1')

    tmp_C(:,:) = eigenvec(:,:,i_spin)

    ! wyj: TODO debug
    ! eigenvec has precision problem, 10^-13
    !print *, myid, 'eigenvec='
    !call print_ham_cpscf(eigenvec)
    !tmp_C(:,:) = 1.0d0
    ! wyj


       ! At this points first_order_ham_polar_reduce_meory_scalapack contains the first-order
       ! hamiltonian in scalapack format
       tmp_H1(:,:) = first_order_ham_polar_reduce_memory_scalapack(:,:,i_spin)
       !call set_full_matrix_real_L_to_U(tmp_H1)
       if (which == 1) then
           ! global_index
           call set_full_matrix_real_L_to_U(tmp_H1)
       else
           ! local_index
           call set_full_matrix_real(tmp_H1)
       endif
    ! wyj: TODO debug
    !print *, myid, 'L+U='
    !call print_ham_cpscf(tmp_H1)


    !------ C^T H1 C ----------------------
        tmp_1 = 0.0d0
        if(myid .eq. 0) then
            print*, "pdgemm, U1, myid=", myid, ", T, N, n_states_k=", n_states_k, ", n_basis=", n_basis, ", n_basis=", n_basis
        endif
        time1 = mpi_wtime() 
        call pdgemm("T","N",n_states_k, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1 , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
      !   call hippdgemm("T","N",n_states_k, n_basis, n_basis, &  ! transa, transb, m, n, k
      !               1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
      !                      tmp_H1 , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
      !               0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
    ! wyj: TODO debug
    !print *, myid, 'C*H1=tmp_1='
    !call print_ham_cpscf(tmp_1)
       tmp_H1 = 0.0d0
        if(myid .eq. 0) then
            print*, "pdgemm, U1, myid=", myid, ", time=", mpi_wtime() - time1
            print*, "pdgemm, U2, myid=", myid, ", N, T, n_states_k=", n_states_k, ", n_states_k=", n_states_k, ", n_basis=", n_basis
        endif
       time1 = mpi_wtime() 
       call pdgemm("N","N",n_states_k, n_states_k, n_basis, &  ! transa, transb, m, n, k
                   1.0d0, tmp_1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                          tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                   0.0d0, tmp_H1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
      !  call hippdgemm("N","N",n_states_k, n_states_k, n_basis, &  ! transa, transb, m, n, k
      !              1.0d0, tmp_1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
      !                     tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
      !              0.0d0, tmp_H1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        if(myid .eq. 0) then
            print*, "pdgemm, U2, myid=", myid, ", time=", mpi_wtime() - time1
        endif
    ! wyj: TODO debug
    !print *, myid, '1*C=tmp_H1='
    !call print_ham_cpscf(tmp_H1)

    !------ U1_pq = (- C^T H1 C)_pq/(E_pp-E_qq) -----------------
    ! Nath: Actually we only need to get the virt-occ and/or occ-virt part of U

       first_order_U_polar_reduce_memory_scalapack(:,:,i_spin) = 0.0d0

       do i_local_row = 1, n_my_rows
       do i_local_col = 1, n_my_cols

          i_global_row = my_row(i_local_row)
          i_global_col = my_col(i_local_col)

          if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

          first_order_U_polar_reduce_memory_scalapack(i_local_row,i_local_col,i_spin) =  &
          ( -tmp_H1(i_local_row,i_local_col)) / &
          ( eigenvalues(i_global_row, i_spin, 1) - eigenvalues(i_global_col, i_spin, 1) )

          !first_order_U_polar_reduce_memory_scalapack(i_local_row,i_local_col,i_spin) =  &
          !( -tmp_H1(i_local_row,i_local_col)) / &
          ! 1 
        !print *, myid, 'u_eigenvalues(', i_global_row, i_spin, '1)', eigenvalues(i_global_row, i_spin, 1)
        !print *, myid, 'u_eigenvalues(', i_global_col, i_spin, '1)', eigenvalues(i_global_col, i_spin, 1)

          endif

          if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

          first_order_U_polar_reduce_memory_scalapack(i_local_row,i_local_col,i_spin) =  &
          ( -tmp_H1(i_local_row,i_local_col)) / &
          ( eigenvalues(i_global_row, i_spin, 1) - eigenvalues(i_global_col, i_spin, 1) )

          !first_order_U_polar_reduce_memory_scalapack(i_local_row,i_local_col,i_spin) =  &
          !( -tmp_H1(i_local_row,i_local_col)) / &
          !  1
        !print *, myid, 'u_eigenvalues(', i_global_row, i_spin, '1)', eigenvalues(i_global_row, i_spin, 1)
        !print *, myid, 'u_eigenvalues(', i_global_col, i_spin, '1)', eigenvalues(i_global_col, i_spin, 1)

          endif

       enddo
       enddo



    deallocate(tmp_C)
    deallocate(tmp_H1)
    deallocate(tmp_1)


    enddo ! i_spin


  end subroutine evaluate_first_order_U_polar_reduce_memory_scalapack_cpscf
!=========================================================================================
!=========================end scalapack part used in DFPT_polar_reduce_memory===========
!=========================================================================================






!=====================================================================================
!=========================begin for scalapack used in DFPT_dielectric=========
!=========================================================================================
!----------------lapack version-----------------|---------scalapack version-------------------
!(1) construct_matrix_complex                   |  construct_first_order_hamiltonian_dielectric_scalapack
!(2) evaluate_first_order_DM_dielectric         |  construct_first_order_dm_dielectric_scalapack
!(2) evaluate_first_order_DM_dielectric         |  get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
!(3) evaluate_first_order_U_dielectric          |  evaluate_first_order_U_dielectric_scalapack


!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_momentum_matrix_dielectric_scalapack
!  NAME
!    construct_momentum_matrix_dielectric_scalapack
!  SYNOPSIS
  subroutine construct_momentum_matrix_dielectric_scalapack( momentum_matrix_sparse,j_coord )
!  PURPOSE
!    Sets the momentum_matrix_sparse in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use dimensions, only: n_hamiltonian_matrix_size_no_symmetry,n_states_k
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    real*8, dimension(n_hamiltonian_matrix_size_no_symmetry), intent(in) :: momentum_matrix_sparse
    integer, intent(in) :: j_coord
!  INPUTS
!    o momentum_matrix_sparse
!  OUTPUT
!    the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_cell, i_col, i_row, lr, lc, idx, i_spin, info, i_basis,i_state, j_state
    real*8, allocatable :: tmp_MM(:,:)
    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_1(:,:)
    complex*16, allocatable :: tmp_MM_complex(:,:)
    complex*16, allocatable :: tmp_C_complex(:,:)
    complex*16, allocatable :: tmp_1_complex(:,:)

    ! construct_hamiltonian_scalapack must not be used when use_local_index is set!

    !if(use_local_index) call aims_stop("INTERNAL ERROR: construct_momentum_matrix_dielectric_scalapack + use_local_index")

    if(real_eigenvectors)then
       momentum_matrix_scalapack(:,:) = 0.
    else
       momentum_matrix_complex_scalapack(:,:) = 0.
    end if



       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_row = 1, n_basis
                lr = l_row(i_row) ! we do not exchange row,col as H0.
                if(lr==0) cycle   ! skip if not local

                if(index_hamiltonian_no_symmetry(1,i_cell,i_row) > 0) then

                   do idx = index_hamiltonian_no_symmetry(1,i_cell,i_row), &
                            index_hamiltonian_no_symmetry(2,i_cell,i_row)

                      i_col = column_index_hamiltonian_no_symmetry(idx)
                      lc = l_col(i_col) ! we do not exchange row,col as H0.
                      if(lc==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         momentum_matrix_scalapack(lr,lc) = momentum_matrix_scalapack(lr,lc) &
                                + dble(k_phase(i_cell,my_k_point) * momentum_matrix_sparse(idx))
                      else
                         momentum_matrix_complex_scalapack (lr,lc) =  &
                         momentum_matrix_complex_scalapack (lr,lc) &
                                + k_phase(i_cell,my_k_point) * momentum_matrix_sparse(idx)
                      end if

                   end do

                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_momentum_matrix_dielectric_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format


   do i_spin = 1, n_spin

    if (real_eigenvectors) then
       allocate(tmp_C(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       !allocate(tmp_MM(mxld, mxcol),stat=info)
       !call check_allocation(info, 'tmp_S1')
       allocate(tmp_1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1')
       tmp_C(:,:) = eigenvec(:,:,i_spin)
       !tmp_MM(:,:) = momentum_matrix_scalapack(:,:)
    !------real-(2) C^+ MM C  ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           momentum_matrix_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        !tmp_MM = 0.0d0
        call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, Omega_MO_scalapack(:,:,j_coord), 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

       deallocate(tmp_C)
       !deallocate(tmp_MM)
       deallocate(tmp_1)
    else
       allocate(tmp_C_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       !allocate(tmp_MM_complex(mxld, mxcol),stat=info)
       !call check_allocation(info, 'tmp_MM_complex')
       allocate(tmp_1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1_complex')
       tmp_C_complex(:,:) = dconjg(eigenvec_complex(:,:,i_spin))
       !tmp_MM_complex(:,:) = momentum_matrix_complex_scalapack(:,:)
    !--!----complex-(2) C^+ MM C  ----------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_states_k(my_k_point), n_basis, n_basis, &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           momentum_matrix_complex_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), tmp_1_complex,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        !tmp_MM_complex = 0.0d0
        call pzgemm("N","N",n_states_k(my_k_point), n_states_k(my_k_point), n_basis, &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), Omega_MO_complex_scalapack(:,:,j_coord), 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

       !if(my_k_point.eq.1) then
       ! print*, "Omega_MO_complex_scalapack", sum(real(Omega_MO_complex_scalapack(:,:,3))),j_coord
       !endif

       deallocate(tmp_C_complex)
       !deallocate(tmp_MM_complex)
       deallocate(tmp_1_complex)

    endif

   enddo !i_spin

  end subroutine construct_momentum_matrix_dielectric_scalapack

!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_momentum_matrix_dielectric_for_elsi_scalapack
!  NAME
!    construct_momentum_matrix_dielectric_for_elsi_scalapack
!  SYNOPSIS
  subroutine construct_momentum_matrix_dielectric_for_elsi_scalapack( momentum_matrix_sparse, mm_scalapack)
!  PURPOSE
!    Sets the momentum_matrix_sparse in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use dimensions, only: n_hamiltonian_matrix_size_no_symmetry,n_states_k
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    real*8, dimension(n_hamiltonian_matrix_size_no_symmetry), intent(in) :: momentum_matrix_sparse
    real*8, dimension(mxld,mxcol), intent(inout) :: mm_scalapack
!  INPUTS
!    o momentum_matrix_sparse
!  OUTPUT
!    the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_cell, i_col, i_row, lr, lc, idx, i_spin, info, i_basis,i_state, j_state

    ! construct_hamiltonian_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_momentum_matrix_dielectric_for_elsi_scalapack + use_local_index")

    if(real_eigenvectors)then
       momentum_matrix_scalapack(:,:) = 0.
       mm_scalapack(:,:) = 0.0d0
    else
       momentum_matrix_complex_scalapack(:,:) = 0.
    end if



       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_row = 1, n_basis
                lr = l_row(i_row) ! we do not exchange row,col as H0.
                if(lr==0) cycle   ! skip if not local

                if(index_hamiltonian_no_symmetry(1,i_cell,i_row) > 0) then

                   do idx = index_hamiltonian_no_symmetry(1,i_cell,i_row), &
                            index_hamiltonian_no_symmetry(2,i_cell,i_row)

                      i_col = column_index_hamiltonian_no_symmetry(idx)
                      lc = l_col(i_col) ! we do not exchange row,col as H0.
                      if(lc==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         momentum_matrix_scalapack(lr,lc) = momentum_matrix_scalapack(lr,lc) &
                                + dble(k_phase(i_cell,my_k_point) * momentum_matrix_sparse(idx))
                         mm_scalapack(lr,lc) = momentum_matrix_scalapack(lr,lc)
                      else
                         momentum_matrix_complex_scalapack (lr,lc) =  &
                         momentum_matrix_complex_scalapack (lr,lc) &
                                + k_phase(i_cell,my_k_point) * momentum_matrix_sparse(idx)
                      end if

                   end do

                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_momentum_matrix_dielectric_for_elsi_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format

  end subroutine construct_momentum_matrix_dielectric_for_elsi_scalapack
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_hamiltonian_dielectric_scalapack
!  NAME
!    construct_first_order_hamiltonian_dielectric_scalapack
!  SYNOPSIS
  subroutine construct_first_order_hamiltonian_dielectric_scalapack( first_order_H_sparse )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    real*8, dimension(n_hamiltonian_matrix_size), intent(in) :: first_order_H_sparse
!  INPUTS
!    o first_order_H_sparse -- the first order Hamilton matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_scalapack must not be used when use_local_index is set!

    !if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_scalapack + use_local_index")

    if(real_eigenvectors)then
       first_order_ham_scalapack(:,:,:) = 0.
    else
       first_order_ham_complex_scalapack(:,:,:) = 0.
    end if

    ! Attention: Only the lower half of the matrices is set here for H1

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_row = 1, n_basis
                lr = l_row(i_row) ! we do not exchange row,col as H0.
                if(lr==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_row) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_row),index_hamiltonian(2,i_cell,i_row)

                      i_col = column_index_hamiltonian(idx)
                      lc = l_col(i_col) ! we do not exchange row,col as H0.
                      if(lc==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         first_order_ham_scalapack(lr,lc,i_spin) = first_order_ham_scalapack (lr,lc,i_spin) &
                                + dble(k_phase(i_cell,my_k_point) * first_order_H_sparse(idx))
                      else
                         first_order_ham_complex_scalapack (lr,lc,i_spin) =  &
                         first_order_ham_complex_scalapack (lr,lc,i_spin) &
                                + k_phase(i_cell,my_k_point) * first_order_H_sparse(idx)
                      end if

                   end do

                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_first_order_hamiltonian_dielectric_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format
    end do

  end subroutine construct_first_order_hamiltonian_dielectric_scalapack


!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_hamiltonian_dielectric_for_elsi_scalapack
!  NAME
!    construct_first_order_hamiltonian_dielectric_for_elsi_scalapack
!  SYNOPSIS
  subroutine construct_first_order_hamiltonian_dielectric_for_elsi_scalapack( first_order_H_sparse, mat_scalapack )
!  PURPOSE
!    Sets the Hamilton matrix in the ScaLAPACK array form a sparse matrix storage.
!  USES
    use mpi_tasks
    use pbc_lists
    use geometry
    use basis
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS
    real*8, dimension(n_hamiltonian_matrix_size), intent(in) :: first_order_H_sparse
    real*8, dimension(mxld,mxcol), intent(inout)  :: mat_scalapack
    
!  INPUTS
!    o first_order_H_sparse -- the first order Hamilton matrix
!  OUTPUT
!    upper half of the ScaLAPACK array ham is set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer:: i_spin, i_cell, i_col, i_row, lr, lc, idx

    ! construct_hamiltonian_scalapack must not be used when use_local_index is set!

    if(use_local_index) call aims_stop("INTERNAL ERROR: construct_hamiltonian_scalapack + use_local_index")

    if(real_eigenvectors)then
       first_order_ham_scalapack(:,:,:) = 0.
       mat_scalapack(:,:) = 0.0d0
    else
       first_order_ham_complex_scalapack(:,:,:) = 0.
    end if

    ! Attention: Only the lower half of the matrices is set here for H1

    do i_spin = 1, n_spin

       select case(packed_matrix_format)

       case(PM_index) !------------------------------------------------

          do i_cell = 1, n_cells_in_hamiltonian-1
             do i_row = 1, n_basis
                lr = l_row(i_row) ! we do not exchange row,col as H0.
                if(lr==0) cycle   ! skip if not local

                if(index_hamiltonian(1,i_cell,i_row) > 0) then

                   do idx = index_hamiltonian(1,i_cell,i_row),index_hamiltonian(2,i_cell,i_row)

                      i_col = column_index_hamiltonian(idx)
                      lc = l_col(i_col) ! we do not exchange row,col as H0.
                      if(lc==0) cycle   ! skip if not local

                      if(real_eigenvectors)then
                         first_order_ham_scalapack(lr,lc,i_spin) = first_order_ham_scalapack (lr,lc,i_spin) &
                                + dble(k_phase(i_cell,my_k_point) * first_order_H_sparse(idx))
                         mat_scalapack(lr,lc) = first_order_ham_scalapack(lr,lc,i_spin)
                      else
                         first_order_ham_complex_scalapack (lr,lc,i_spin) =  &
                         first_order_ham_complex_scalapack (lr,lc,i_spin) &
                                + k_phase(i_cell,my_k_point) * first_order_H_sparse(idx)
                      end if

                   end do

                end if
             end do
          end do ! i_cell

       case default !---------------------------------------------------------

          write(use_unit,*) 'Error: construct_first_order_hamiltonian_dielectric_for_elsi_scalapack does not support non-packed matrices.'
          call aims_stop

       end select ! packed_matrix_format
    end do

  end subroutine construct_first_order_hamiltonian_dielectric_for_elsi_scalapack



!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/construct_first_order_dm_dielectric_scalapack
!  NAME
!    construct_first_order_dm_dielectric_scalapack
!  SYNOPSIS
  subroutine construct_first_order_dm_dielectric_scalapack(occ_numbers, i_spin)
!  PURPOSE
!    Construct first_order density matrix at my_k_point in ScaLAPACK
!    here C+ = C*T
!    DM1(my_k_point) = C*occ_number* (-C^+ S C) * C
!                      + (C U) C^+*occ_number + C*occ_number (C U)^+
!
!  USES
    use mpi_tasks
    use pbc_lists
    use runtime_choices
    use dimensions, only: n_states_k
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    integer, intent(IN) :: i_spin
!  INPUTS
!    o occ_numbers -- the occupation numbers
!    o i_spin -- spin component
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, info , j_state
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    real*8, allocatable :: tmp_C_occ(:,:)
    real*8, allocatable :: tmp_1(:,:)
    real*8, allocatable :: tmp_2(:,:)
    real*8, allocatable :: tmp_C1(:,:)

    complex*16, allocatable :: tmp_C_complex(:,:)
    complex*16, allocatable :: tmp_C_occ_complex(:,:)
    complex*16, allocatable :: tmp_1_complex(:,:)
    complex*16, allocatable :: tmp_2_complex(:,:)
    complex*16, allocatable :: tmp_C1_complex(:,:)

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix
    max_occ_number = 0

    !Note:  Here occ_number have already contained with k_weights, but do not have in lapack version.
    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, my_k_point) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

    if (real_eigenvectors) then

    !------real-(1).prepare C, C_occ, S^(1) ----------------------
       allocate(tmp_C(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C')
       allocate(tmp_C_occ(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ')

       allocate(tmp_1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1')
       allocate(tmp_2(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_2')
       allocate(tmp_C1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C1')

       first_order_ham_scalapack(:,:,i_spin) = 0d0

       tmp_C(:,:) = eigenvec(:,:,i_spin)
       do i_state = 1, n_states
          if (occ_numbers(i_state, i_spin, my_k_point) > 1.e-6) then
             if(l_col(i_state)>0) then
       ! Nath: double check if the occupation numbers should enter this routine when DFPT_width is called. It seems that yes.
       ! If not, careful that the occ_numbers here are different than the ones in evaluate_first_order_DM (lapack version)
               tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
               occ_numbers(i_state, i_spin, my_k_point)
             endif
          elseif(l_col(i_state).ne.0) then
            tmp_C_occ(:,l_col(i_state)) = 0d0
          endif
       end do

       print *, 'wyj_tmp_C(*,1)=', 'myid=', myid, tmp_C(:,1)

    !------real-(3).DM^(1)_(ov+vo) ----------------------
    !    = (C U) C^+*occ_number + C*occ_number (C U)^+
        tmp_C1 = 0.0d0
        !call pdgemm("N","N",n_basis, n_states_k(my_k_point), n_states_k(my_k_point), &  ! transa, transb, m, n, k
        !            1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
        !                   first_order_U_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
        !            0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        ! wyj:TODO: debug here
        call pdgemm("N","T",n_basis, n_states_k(my_k_point), n_states_k(my_k_point), &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_C1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

        tmp_1 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        first_order_ham_scalapack(:,:,i_spin) = first_order_ham_scalapack(:,:,i_spin) + tmp_1(:,:)

        tmp_2 = 0.0d0
        call pdgemm("N","T",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C_occ,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_2,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        first_order_ham_scalapack(:,:,i_spin) = first_order_ham_scalapack(:,:,i_spin) + tmp_2(:,:)

       print *, 'wyj_tmp_C1(*,1)=', 'myid=', myid, tmp_C1(:,1)
       print *, 'wyj_tmp_C_occ(*,1)=', 'myid=', myid, tmp_C_occ(:,1)
       print *, 'wyj_tmp_1(*,1)=', 'myid=', myid, tmp_1(:,1)
       print *, 'wyj_tmp_2(*,1)=', 'myid=', myid, tmp_2(:,1)
!-------------begin shanghui's debug tool for scalapack version---------------
        print *, "wyj debug DM real"
       do  i_local_row = 1, n_my_rows
       do  i_local_col = 1, n_my_cols
           i_global_row = my_row(i_local_row)
           i_global_col = my_col(i_local_col)
           if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.1) then
           write(use_unit,*) 'real: C(C+S1C)C+ (11):', first_order_ham_scalapack(i_local_row,i_local_col,i_spin)
           endif
           if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.2) then
           write(use_unit,*) 'real: C(C+S1C)C+ (12):', first_order_ham_scalapack(i_local_row,i_local_col,i_spin)
           endif
       enddo
       enddo
!-------------end shanghui's debug tool for scalapack version---------------

       deallocate(tmp_C)
       deallocate(tmp_C_occ)
       deallocate(tmp_1)
       deallocate(tmp_2)

    else

    !------complex-(1).prepare C, C_occ, S^(1) ----------------------
       allocate(tmp_C_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       allocate(tmp_C_occ_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_occ_complex')

       allocate(tmp_1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1_complex')
       allocate(tmp_2_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_2_complex')
       allocate(tmp_C1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C1_complex')



       first_order_ham_complex_scalapack(:,:,i_spin) = 0d0

              !if (DFPT_width .eq. 0.0) then ! No metal
              ! tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state)) * &
              ! occ_numbers(i_state, i_spin, my_k_point)
              !else
              ! tmp_C_occ(:,l_col(i_state)) = tmp_C(:,l_col(i_state))
              !endif ! metal

       ! Remeber that C_aims = C^*, so first we need to make C=C_aims^*
       tmp_C_complex(:,:) = dconjg(eigenvec_complex(:,:,i_spin))
       do i_state = 1, n_states
          if (occ_numbers(i_state, i_spin, my_k_point) > 1.e-6) then
             if(l_col(i_state)>0) then
               tmp_C_occ_complex(:,l_col(i_state)) = tmp_C_complex(:,l_col(i_state)) * &
                                             occ_numbers(i_state, i_spin, my_k_point)
             endif
          elseif(l_col(i_state).ne.0) then
               tmp_C_occ_complex(:,l_col(i_state)) = 0d0
          endif
       end do




    !------complex-(3).DM^(1)_(ov+vo) ----------------------
    !    = (C U) C^+*occ_number + C*occ_number (C U)^+
        tmp_C1_complex = 0.0d0
        call pzgemm("N","N",n_basis, n_states_k(my_k_point), n_states_k(my_k_point), &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           first_order_U_complex_scalapack, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), tmp_C1_complex,  1, 1, sc_desc)       ! beta,  c, ic, jc, desc_c

        tmp_1_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_C1_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_occ_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), tmp_1_complex,  1, 1, sc_desc)            ! beta,  c, ic, jc, desc_c
        first_order_ham_complex_scalapack(:,:,i_spin) = first_order_ham_complex_scalapack(:,:,i_spin) &
                                                      + tmp_1_complex(:,:)

        tmp_2_complex = 0.0d0
        call pzgemm("N","C",n_basis, n_basis, max_occ_number, &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_C_occ_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C1_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), tmp_2_complex,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c
        first_order_ham_complex_scalapack(:,:,i_spin) = first_order_ham_complex_scalapack(:,:,i_spin) &
                                                      + tmp_2_complex(:,:)

!-------------begin shanghui's debug tool for scalapack version---------------
        print *, "wyj debug DM real"
       do  i_local_row = 1, n_my_rows
       do  i_local_col = 1, n_my_cols
           i_global_row = my_row(i_local_row)
           i_global_col = my_col(i_local_col)
           if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.1) then
           write(use_unit,*) 'C(C+S1C)C+ (11):', first_order_ham_complex_scalapack(i_local_row,i_local_col,i_spin)
           endif
           if(my_k_point.eq.1.and.i_global_row.eq.1.and.i_global_col.eq.2) then
           write(use_unit,*) 'C(C+S1C)C+ (12):', first_order_ham_complex_scalapack(i_local_row,i_local_col,i_spin)
           endif
       enddo
       enddo
!-------------end shanghui's debug tool for scalapack version---------------
       !if(my_k_point.eq.1) then
       !do i_state = 1, n_states
       !do j_state = 1, n_states

       !   if(l_row(i_state).ne.0.and.l_col(j_state).ne.0) then
       !   write(use_unit,*) 'C^+UC', i_state, j_state,  &
       !               tmp_1_complex(l_row(i_state),l_col(j_state))+tmp_2_complex(l_row(i_state),l_col(j_state))
       !   endif

       !enddo
       !enddo
       !endif


       deallocate(tmp_C_complex)
       deallocate(tmp_C_occ_complex)
       deallocate(tmp_1_complex)
       deallocate(tmp_2_complex)


    endif

  end subroutine construct_first_order_dm_dielectric_scalapack

!******
!-----------------------------------------------------------
!****s* scalapack_wrapper/evaluate_first_order_U_dielectric_scalapack
!  NAME
!    evaluate_first_order_U_dielectric_scalapack
!  SYNOPSIS
  subroutine evaluate_first_order_U_dielectric_scalapack(occ_numbers, eigenvalues,dP_dE,j_coord)
!  PURPOSE
!    evaluate first_order_U
!    here C^+ = C*T
!    U1_pq(my_k_point) =  - (C^+ MM C)_pq /(E_pp-E_qq)^2  - (C^+ H1 C)_pq/(E_pp-E_qq)
!
!
!  USES
    use mpi_tasks
    use pbc_lists
    use physics, only : chemical_potential
    use runtime_choices
    use dimensions, only: n_states_k
    implicit none
!  ARGUMENTS
    real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
    real*8, dimension(n_states, n_spin, n_k_points) :: eigenvalues
    complex*16 :: dP_dE(3,3)
    integer, intent(in) :: j_coord
!  INPUTS
!    o occ_numbers -- the occupation numbers with k_weights
!  OUTPUT
!    o the array ham/ham_complex is overwritten by the density matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: max_occ_number, i_state, j_state, info, i_spin, i_coord
    integer :: i_local_row, i_local_col, i_global_row, i_global_col

    real*8, allocatable :: tmp_C(:,:)
    !real*8, allocatable :: tmp_MM(:,:)
    real*8, allocatable :: tmp_H1(:,:)
    real*8, allocatable :: tmp_1(:,:)

    complex*16, allocatable :: tmp_C_complex(:,:)
    !complex*16, allocatable :: tmp_MM_complex(:,:)
    complex*16, allocatable :: tmp_H1_complex(:,:)
    complex*16, allocatable :: tmp_1_complex(:,:)
    real*8 :: theta_focc, theta_fvirt, theta_virtocc, lim ! integral of the smearing function

    ! The work in this routine must be done only on the working set
    if(my_scalapack_id>=npcol*nprow) return

    do i_spin = 1, n_spin

    ! we use first_order_ham/first_order_ham_complex as storage area for first order density matrix
    max_occ_number = 0

    !Note:  Here occ_number have already contained with k_weights, but do not have in lapack version.
    do i_state = 1, n_states
       if(occ_numbers(i_state, i_spin, my_k_point) > 1.e-6) then
          max_occ_number = i_state
       endif
    enddo

    if (real_eigenvectors) then

    !------real-(1).prepare C, C_eig, S^(1), H^(1)----------------------
       allocate(tmp_C(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C')
       !allocate(tmp_MM(mxld, mxcol),stat=info)
       !call check_allocation(info, 'tmp_S1')
       allocate(tmp_H1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_H1')
       allocate(tmp_1(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1')

       tmp_C(:,:) = eigenvec(:,:,i_spin)
       !tmp_MM(:,:) = momentum_matrix_scalapack(:,:)
       tmp_H1(:,:) = first_order_ham_scalapack(:,:, i_spin)
       call set_full_matrix_real_L_to_U(tmp_H1)


    !------real-(2) C^+ MM C  ----------------------
        !tmp_1 = 0.0d0
        !call pdgemm("T","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
        !            1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
        !                   tmp_MM, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
        !            0.0d0, tmp_1,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        !tmp_MM = 0.0d0
        !call pdgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
        !            1.0d0, tmp_1, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
        !                   tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
        !            0.0d0, tmp_MM, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------real-(3) C^+ H1 C ----------------------
        tmp_1 = 0.0d0
        call pdgemm("T","N",n_states_k(my_k_point), n_basis, n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_C,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1 , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_H1 = 0.0d0
        call pdgemm("N","N",n_states_k(my_k_point), n_states_k(my_k_point), n_basis, &  ! transa, transb, m, n, k
                    1.0d0, tmp_1,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    0.0d0, tmp_H1,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------real-(4) U1_pq  -----------------
        first_order_U_scalapack = 0.0d0

        do  i_local_row = 1, n_my_rows
        do  i_local_col = 1, n_my_cols

            i_global_row = my_row(i_local_row)
            i_global_col = my_col(i_local_col)

            if (DFPT_width .eq. 0.0) then

              if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

                first_order_U_scalapack(i_local_row,i_local_col) =  &
                !- tmp_MM(i_local_row,i_local_col) /              &
                - Omega_MO_scalapack(i_local_row,i_local_col,j_coord) / &
                ( eigenvalues(i_global_row, i_spin, my_k_point) - &
                 eigenvalues(i_global_col, i_spin, my_k_point))**2  &
                - tmp_H1(i_local_row,i_local_col) / &
                ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

                  print *, 'wyj_first_order_U_scalapack(', i_global_row, i_global_col, ')=', &
                      first_order_U_scalapack(i_local_row, i_local_col)
                  print *, 'eigenvalues(', i_global_col, ')=', eigenvalues(i_global_col, i_spin, my_k_point)
                  print *, 'eigenvalues(', i_global_row, ')=', eigenvalues(i_global_row, i_spin, my_k_point)
                  print *, 'Omega_MO_scalapack(', i_global_row, i_global_col, ')=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
                  print *, 'tmp_H1(', i_global_row, i_global_col, ')=', tmp_H1(i_local_row, i_local_col)

              !if (my_k_point == 1 .and. i_global_row == max_occ_number + 1 .and. i_global_col == max_occ_number) then
              !    print *, 'wyj_first_order_U_scalapack(1,1)=', &
              !        first_order_U_scalapack(i_local_row, i_local_col)
              !    print *, 'eigenvalues(col1,x,x)=', eigenvalues(i_global_col, i_spin, my_k_point)
              !    print *, 'Omega_MO_scalapack(1,1,x)=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
              !    print *, 'tmp_H1(1,1)=', tmp_H1(i_local_row, i_local_col)
              !endif

              !if (my_k_point == 1 .and. i_global_row == max_occ_number + 2 .and. i_global_col == max_occ_number) then
              !    print *, 'wyj_first_order_U_scalapack(1,2)=', &
              !        first_order_U_scalapack(i_local_row, i_local_col)
              !    print *, 'eigenvalues(col2,x,x)=', eigenvalues(i_global_col, i_spin, my_k_point)
              !    print *, 'Omega_MO_scalapack(1,2,x)=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
              !    print *, 'tmp_H1(1,2)=', tmp_H1(i_local_row, i_local_col)
              !endif

              !if (my_k_point == 1 .and. i_global_row == max_occ_number + 3 .and. i_global_col == max_occ_number) then
              !    print *, 'wyj_first_order_U_scalapack(1,3)=', &
              !        first_order_U_scalapack(i_local_row, i_local_col)
              !    print *, 'eigenvalues(col3,x,x)=', eigenvalues(i_global_col, i_spin, my_k_point)
              !    print *, 'Omega_MO_scalapack(1,3,x)=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
              !    print *, 'tmp_H1(1,3)=', tmp_H1(i_local_row, i_local_col)
              !endif
              endif

              if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

                first_order_U_scalapack(i_local_row,i_local_col) =  &
                !- tmp_MM(i_local_row,i_local_col) /              &
                - Omega_MO_scalapack(i_local_row,i_local_col,j_coord) /              &
                ( eigenvalues(i_global_row, i_spin, my_k_point) - &
                 eigenvalues(i_global_col, i_spin, my_k_point))**2  &
                - tmp_H1(i_local_row,i_local_col) / &
                ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))
                  print *, 'wyj_first_order_U_scalapack(', i_global_row, i_global_col, ')=', &
                      first_order_U_scalapack(i_local_row, i_local_col)
                  print *, 'eigenvalues(', i_global_col, ')=', eigenvalues(i_global_col, i_spin, my_k_point)
                  print *, 'eigenvalues(', i_global_row, ')=', eigenvalues(i_global_row, i_spin, my_k_point)
                  print *, 'Omega_MO_scalapack(', i_global_row, i_global_col, ')=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
                  print *, 'tmp_H1(', i_global_row, i_global_col, ')=', tmp_H1(i_local_row, i_local_col)

              !if (my_k_point == 1 .and. i_global_row == max_occ_number .and. i_global_col == max_occ_number + 1) then
              !    print *, 'wyj_first_order_U_scalapack(1,1)=', &
              !        first_order_U_scalapack(i_local_row, i_local_col)
              !    print *, 'eigenvalues(col1,x,x)=', eigenvalues(i_global_col, i_spin, my_k_point)
              !    print *, 'Omega_MO_scalapack(1,1,x)=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
              !    print *, 'tmp_H1(1,1)=', tmp_H1(i_local_row, i_local_col)
              !endif

              !if (my_k_point == 1 .and. i_global_row == max_occ_number .and. i_global_col == max_occ_number + 2) then
              !    print *, 'wyj_first_order_U_scalapack(1,2)=', &
              !        first_order_U_scalapack(i_local_row, i_local_col)
              !    print *, 'eigenvalues(col2,x,x)=', eigenvalues(i_global_col, i_spin, my_k_point)
              !    print *, 'Omega_MO_scalapack(1,2,x)=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
              !    print *, 'tmp_H1(1,2)=', tmp_H1(i_local_row, i_local_col)
              !endif

              !if (my_k_point == 1 .and. i_global_row == max_occ_number .and. i_global_col == max_occ_number + 3) then
              !    print *, 'wyj_first_order_U_scalapack(1,3)=', &
              !        first_order_U_scalapack(i_local_row, i_local_col)
              !    print *, 'eigenvalues(col3,x,x)=', eigenvalues(i_global_col, i_spin, my_k_point)
              !    print *, 'Omega_MO_scalapack(1,3,x)=', Omega_MO_scalapack(i_local_row, i_local_col, j_coord)
              !    print *, 'tmp_H1(1,3)=', tmp_H1(i_local_row, i_local_col)
              !endif
              endif


            else ! Metals/Tiny energy differences
              ! See 'Lattice dynamics of metals from density-functional perturbation theory', Stefano de Gironcoli,
              ! Phys. Rev.B51, 6773(R), (doi.org/10.1103/PhysRevB.51.6773)

              if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then
                theta_focc = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_col,1,my_k_point))/(2.*DFPT_width)) )
                theta_fvirt = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_row,1,my_k_point))/(2.*DFPT_width)) )
                theta_virtocc = 0.5*(1.+tanh((eigenvalues(i_global_row,1,my_k_point)-eigenvalues(i_global_col,1,my_k_point))/(2.*DFPT_width)) )
                lim = -0.5*DFPT_width*1./(1+cosh((chemical_potential-eigenvalues(i_global_col,1,my_k_point))/DFPT_width))

                if (abs(eigenvalues(i_global_col,1,my_k_point)-eigenvalues(i_global_row,1,my_k_point)) .gt. 5.e-2 ) then

                  first_order_U_scalapack(i_local_row,i_local_col) =  &
                  !- tmp_MM(i_local_row,i_local_col) /              &
                  - Omega_MO_scalapack(i_local_row,i_local_col,j_coord)*theta_virtocc*(theta_focc-theta_fvirt) / &
                  ( eigenvalues(i_global_row, i_spin, my_k_point) - &
                   eigenvalues(i_global_col, i_spin, my_k_point))**2  &
                  - tmp_H1(i_local_row,i_local_col)*theta_virtocc*(theta_focc-theta_fvirt) / &
                  ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))
                else ! Switch to the limit when the difference between 2 eigenvalues is too small
                  first_order_U_scalapack( i_local_row,i_local_col)=         &
                  -Omega_MO_scalapack(i_local_row,i_local_col,j_coord)*lim&
                  -tmp_H1(i_local_row,i_local_col)*lim
                endif !limit

              endif ! unocc/occ

              if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then
                theta_focc = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_row,1,my_k_point))/(2.*DFPT_width)) )
                theta_fvirt = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_col,1,my_k_point))/(2.*DFPT_width)) )
                theta_virtocc = 0.5*(1.+tanh((eigenvalues(i_global_col,1,my_k_point)-eigenvalues(i_global_row,1,my_k_point))/(2.*DFPT_width)) )
                lim = -0.5*DFPT_width*1./(1+cosh((chemical_potential-eigenvalues(i_global_row,1,my_k_point))/DFPT_width))

                if (abs(eigenvalues(i_global_col,1,my_k_point)-eigenvalues(i_global_row,1,my_k_point)) .gt. 5.e-2 ) then
                  first_order_U_scalapack(i_local_row,i_local_col) =  &
                  !- tmp_MM(i_local_row,i_local_col) /              &
                  - Omega_MO_scalapack(i_local_row,i_local_col,j_coord) /              &
                  ( eigenvalues(i_global_row, i_spin, my_k_point) - &
                   eigenvalues(i_global_col, i_spin, my_k_point))**2  &
                  - tmp_H1(i_local_row,i_local_col) / &
                  ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

                else ! Switch to the limit when the difference between 2 eigenvalues is too small
                  first_order_U_scalapack( i_local_row,i_local_col)=         &
                  -Omega_MO_scalapack(i_local_row,i_local_col,j_coord)*lim&
                  -tmp_H1(i_local_row,i_local_col)*lim

                endif ! limit
              endif !occ/unocc

            endif ! DFPT_width

       enddo
       enddo


!---------------------begin test for diagonal dielectric constant-----------------
        dP_dE(:,j_coord) = (0.0d0,0.0d0)
        do  i_local_row = 1, n_my_rows
        do  i_local_col = 1, n_my_cols

            i_global_row = my_row(i_local_row)
            i_global_col = my_col(i_local_col)

            if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

            do i_coord=1,3
            dP_dE(i_coord,j_coord) = dP_dE(i_coord,j_coord) - &
            !4.0d0*k_weights(my_k_point)*tmp_MM(i_local_row,i_local_col) / &
            4.0d0*k_weights(my_k_point)*Omega_MO_scalapack(i_local_row,i_local_col,i_coord) / &
           ( eigenvalues(i_global_row, i_spin, my_k_point) -   &
             eigenvalues(i_global_col, i_spin, my_k_point) ) * &
           (-first_order_U_scalapack(i_local_row,i_local_col))
            enddo

           endif

       enddo
       enddo

!---------------------end test for diagonal dielectric constant-----------------


       deallocate(tmp_C)
       !deallocate(tmp_MM)
       deallocate(tmp_H1)
       deallocate(tmp_1)

    else

        print *, 'wyj_U: complex part*****'
    !------complex-(1).prepare C, C_eig, S^(1), H^(1)----------------------
       allocate(tmp_C_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_C_complex')
       !allocate(tmp_MM_complex(mxld, mxcol),stat=info)
       !call check_allocation(info, 'tmp_MM_complex')
       allocate(tmp_H1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_H1_complex')
       allocate(tmp_1_complex(mxld, mxcol),stat=info)
       call check_allocation(info, 'tmp_1_complex')

       tmp_C_complex(:,:) = dconjg(eigenvec_complex(:,:,i_spin))
       !tmp_MM_complex(:,:) = momentum_matrix_complex_scalapack(:,:)
       tmp_H1_complex(:,:) = first_order_ham_complex_scalapack(:,:,i_spin)
       call set_full_matrix_complex_L_to_U(tmp_H1_complex)

    !------complex-(2) C^+ MM C  ----------------------
        !tmp_1_complex = 0.0d0
        !call pzgemm("C","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
        !            (1.0d0,0.0d0), tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
        !                   tmp_MM_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
        !            (0.0d0,0.0d0), tmp_1_complex,  1, 1, sc_desc)         ! beta,  c, ic, jc, desc_c

        !tmp_MM_complex = 0.0d0
        !call pzgemm("N","N",n_basis, n_basis, n_basis, &  ! transa, transb, m, n, k
        !            (1.0d0,0.0d0), tmp_1_complex, 1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
        !                   tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
        !            (0.0d0,0.0d0), tmp_MM_complex, 1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c

       !if(my_k_point.eq.1) then
       ! print*, "tmp_MM_complex_U", sum(real(tmp_MM_complex(:,:)))
       !endif

    !------complex-(3) C^+ H1 C ----------------------
        tmp_1_complex = 0.0d0
        call pzgemm("C","N",n_states_k(my_k_point), n_basis, n_basis, &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_C_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_H1_complex , 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), tmp_1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c
        tmp_H1_complex = 0.0d0
        call pzgemm("N","N",n_states_k(my_k_point), n_states_k(my_k_point), n_basis, &  ! transa, transb, m, n, k
                    (1.0d0,0.0d0), tmp_1_complex,  1, 1, sc_desc,      &  ! alpha, a, ia, ja, desc_a
                           tmp_C_complex, 1, 1, sc_desc,      &  !        b, ib, jb, desc_b
                    (0.0d0,0.0d0), tmp_H1_complex,  1, 1, sc_desc)        ! beta,  c, ic, jc, desc_c


    !------real-(4) U1_pq = (C^+ S1 C E - C^+ H1 C)_pq/(E_pp-E_qq) -----------------
        first_order_U_complex_scalapack = 0.0d0


        do  i_local_row = 1, n_my_rows
        do  i_local_col = 1, n_my_cols

            i_global_row = my_row(i_local_row)
            i_global_col = my_col(i_local_col)

            if (DFPT_width .eq. 0.0) then

             if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then

               first_order_U_complex_scalapack(i_local_row,i_local_col) =  &
               !- tmp_MM_complex(i_local_row,i_local_col) / &
               - Omega_MO_complex_scalapack(i_local_row,i_local_col,j_coord) / &
               ( eigenvalues(i_global_row, i_spin, my_k_point) - &
               eigenvalues(i_global_col, i_spin, my_k_point))**2 &
               - tmp_H1_complex(i_local_row,i_local_col) / &
               ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

             endif

             if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

               first_order_U_complex_scalapack(i_local_row,i_local_col) =  &
               !- tmp_MM_complex(i_local_row,i_local_col) / &
               - Omega_MO_complex_scalapack(i_local_row,i_local_col,j_coord) / &
               ( eigenvalues(i_global_row, i_spin, my_k_point) - &
               eigenvalues(i_global_col, i_spin, my_k_point))**2 &
               - tmp_H1_complex(i_local_row,i_local_col) / &
               ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))
             endif

            else ! Metals/Tiny energy differences
              ! See 'Lattice dynamics of metals from density-functional perturbation theory', Stefano de Gironcoli,
              ! Phys. Rev.B51, 6773(R), (doi.org/10.1103/PhysRevB.51.6773)

              if(i_global_row.gt.max_occ_number.and.i_global_col.le.max_occ_number) then
                theta_focc = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_col,1,my_k_point))/(2.*DFPT_width)) )
                theta_fvirt = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_row,1,my_k_point))/(2.*DFPT_width)) )
                theta_virtocc = 0.5*(1.+tanh((eigenvalues(i_global_row,1,my_k_point)-eigenvalues(i_global_col,1,my_k_point))/(2.*DFPT_width)) )
                lim = -0.5*DFPT_width*1./(1+cosh((chemical_potential-eigenvalues(i_global_col,1,my_k_point))/DFPT_width))

                if (abs(eigenvalues(i_global_col,1,my_k_point)-eigenvalues(i_global_row,1,my_k_point)) .gt. 5.e-2 ) then

                 first_order_U_complex_scalapack(i_local_row,i_local_col) =  &
                 !- tmp_MM_complex(i_local_row,i_local_col) / &
                 - Omega_MO_complex_scalapack(i_local_row,i_local_col,j_coord)*theta_virtocc*(theta_focc-theta_fvirt) / &
                 ( eigenvalues(i_global_row, i_spin, my_k_point) - &
                 eigenvalues(i_global_col, i_spin, my_k_point))**2 &
                 - tmp_H1_complex(i_local_row,i_local_col)*theta_virtocc*(theta_focc-theta_fvirt) / &
                 ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

                else ! Switch to the limit when the difference between 2 eigenvalues is too small
                  first_order_U_complex_scalapack( i_local_row,i_local_col)=         &
                  -Omega_MO_complex_scalapack(i_local_row,i_local_col,j_coord)*lim&
                  -tmp_H1_complex(i_local_row,i_local_col)*lim

                endif ! limit
              endif !unocc/occ

              if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

                theta_focc = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_row,1,my_k_point))/(2.*DFPT_width)) )
                theta_fvirt = 0.5*(1.+tanh((chemical_potential-eigenvalues(i_global_col,1,my_k_point))/(2.*DFPT_width)) )
                theta_virtocc = 0.5*(1.+tanh((eigenvalues(i_global_col,1,my_k_point)-eigenvalues(i_global_row,1,my_k_point))/(2.*DFPT_width)) )
                lim = -0.5*DFPT_width*1./(1+cosh((chemical_potential-eigenvalues(i_global_row,1,my_k_point))/DFPT_width))

                if (abs(eigenvalues(i_global_col,1,my_k_point)-eigenvalues(i_global_row,1,my_k_point)) .gt. 5.e-2 ) then

                 first_order_U_complex_scalapack(i_local_row,i_local_col) =  &
                 !- tmp_MM_complex(i_local_row,i_local_col) / &
                 - Omega_MO_complex_scalapack(i_local_row,i_local_col,j_coord)*theta_virtocc*(theta_focc-theta_fvirt) / &
                 ( eigenvalues(i_global_row, i_spin, my_k_point) - &
                 eigenvalues(i_global_col, i_spin, my_k_point))**2 &
                 - tmp_H1_complex(i_local_row,i_local_col)*theta_virtocc*(theta_focc-theta_fvirt) / &
                 ( eigenvalues(i_global_row, i_spin, my_k_point) - eigenvalues(i_global_col, i_spin, my_k_point))

                else ! Switch to the limit when the difference between 2 eigenvalues is too small
                  first_order_U_complex_scalapack( i_local_row,i_local_col)=         &
                  -Omega_MO_complex_scalapack(i_local_row,i_local_col,j_coord)*lim&
                  -tmp_H1_complex(i_local_row,i_local_col)*lim

                endif ! limit

              endif !occ/unocc

            endif ! width

       enddo
       enddo

!---------------------begin test for diagonal dielectric constant-----------------
        dP_dE(:,j_coord) = (0.0d0,0.0d0)
        do  i_local_row = 1, n_my_rows
        do  i_local_col = 1, n_my_cols

            i_global_row = my_row(i_local_row)
            i_global_col = my_col(i_local_col)

            if(i_global_row.le.max_occ_number.and.i_global_col.gt.max_occ_number) then

            do i_coord=1,3
            ! dP_dE = -4 *k_weights* Omega_ij * U_ji where U_ij=-(U_ji)^*
            dP_dE(i_coord,j_coord) = dP_dE(i_coord,j_coord) - &
            !4.0d0*k_weights(my_k_point)*tmp_MM_complex(i_local_row,i_local_col) / &
            4.0d0*k_weights(my_k_point)*Omega_MO_complex_scalapack(i_local_row,i_local_col,i_coord) / &
           ( eigenvalues(i_global_row, i_spin, my_k_point) -   &
             eigenvalues(i_global_col, i_spin, my_k_point) ) * &
           (-dconjg(first_order_U_complex_scalapack(i_local_row,i_local_col)))
            enddo

           endif

       enddo
       enddo
!---------------------end test for diagonal dielectric constant-----------------

       deallocate(tmp_C_complex)
       !deallocate(tmp_MM_complex)
       deallocate(tmp_H1_complex)
       deallocate(tmp_1_complex)

    endif

   enddo ! i_spin

  end subroutine evaluate_first_order_U_dielectric_scalapack

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
!  NAME
!    get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_complex_sparse_matrix_dielectric_scalapack( matrix_sparse, i_spin )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8 :: matrix_sparse(n_hamiltonian_matrix_size)
    integer :: i_spin
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = 0.0d0

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) !
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   first_order_ham_scalapack(lr,lc,i_spin)*dble(k_phase(i_cell,my_k_point))
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble( first_order_ham_complex_scalapack(lr,lc,i_spin)   &
                   * dconjg(k_phase(i_cell,my_k_point)) )


             endif
          end do
       end do
    end do

  end subroutine get_first_order_dm_complex_sparse_matrix_dielectric_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack
!  NAME
!    get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack
!  SYNOPSIS
   subroutine get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack(mat, matrix_sparse )
!  PURPOSE
!    Gets a sparse matrix array from ScaLAPACK
!  USES
    use localorb_io
    use mpi_tasks
    use pbc_lists
    implicit none
!  ARGUMENTS
    real*8, dimension(mxld,mxcol), intent(in)  :: mat
    real*8, dimension(n_hamiltonian_matrix_size), intent(out) :: matrix_sparse
    
!  INPUTS
!    o i_spin -- the spin channel
!  OUTPUT
!    o matrix_sparse -- set to the contents of ham/ham_complex
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_cell, i_bas1, i_bas2, i_index, lr, lc

    character*200 :: info_str

    if (packed_matrix_format /= PM_index) then
       write(info_str, '(A)') '* ERROR: get_sparse_matrix_scalapack works only for packed_matrix_format == PM_index'
       call localorb_info(info_str, use_unit, '(A)')
       call aims_stop
    end if

    matrix_sparse = 0.0d0

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_bas1 = 1, n_basis
          lr = l_row(i_bas1) !
          if(lr==0) cycle    ! skip if not local

          if( index_hamiltonian(1,i_cell,i_bas1) == 0 ) cycle ! no entries

          do i_index = index_hamiltonian(1,i_cell,i_bas1), index_hamiltonian(2,i_cell,i_bas1)

             i_bas2 = column_index_hamiltonian(i_index)
             lc = l_col(i_bas2) !
             if (lc==0) cycle

             if (real_eigenvectors) then
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   mat(lr,lc)*dble(k_phase(i_cell,my_k_point))
             else
                matrix_sparse(i_index) = &
                   matrix_sparse(i_index) + &
                   dble( mat(lr,lc)   &
                   * dconjg(k_phase(i_cell,my_k_point)) )


             endif
          end do
       end do
    end do

  end subroutine get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack
!=========================================================================================
!=========================end for scalapack used in DFPT_dielectric=======================
!=========================================================================================


!-----------------------------------------------------------------------------------
!****s* FHI-aims/print_sparse_to_dense_global_index_cpscf
!  NAME
!    print_sparse_to_dense_global_index_cpscf
!  SYNOPSIS

subroutine print_sparse_to_dense_global_index_cpscf(sparse_matrix)
!  PURPOSE
!    The subroutine change the density matrix components belongs to non-zero basis functions
!    to density_matrix_con.
!
!  USES

  use dimensions
  use pbc_lists
  use mpi_tasks, only: myid
  implicit none

!  ARGUMENTS

  real*8 :: sparse_matrix(*)

!  INPUTS
!    o sparse_matrix -- total density matrix (packed matrix format)
!
!  OUTPUT
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
! SOURCE

  integer :: i_basis,j_basis, i_start, i_end, i_place

  do i_basis = 1, n_basis
 
     i_start = index_hamiltonian(1,1, i_basis)
     i_end   = index_hamiltonian(2,1, i_basis)

  do j_basis= 1, n_basis
 
     do i_place = i_start, i_end, 1

     if( column_index_hamiltonian(i_place) == j_basis)then
         if (j_basis == 1) then
             print *, myid, 'global_H/DM', '(', i_basis, j_basis, ')', sparse_matrix(i_place)
             !if (i_basis .ne. j_basis) then
             !    print *, myid, 'H/DM', '(', j_basis, i_basis, ')', sparse_matrix(i_place)
             !endif
         endif
     endif

     enddo ! i_place

  enddo ! j_basis
  enddo ! i_basis 

end subroutine print_sparse_to_dense_global_index_cpscf 
!******

end module scalapack_wrapper

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! Below are two modified ScalaPack routines, PDSYEVD_X and PZHEEVD_X
!
! These are declared OUTSIDE the module because otherways the usage
! of a REAL work array as a COMPLEX argument is not possible.
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

      SUBROUTINE PDSYEVD_X(JOBZ,UPLO,N,A,IA,JA,DESCA,W,Z,IZ,JZ,DESCZ,   &
     &                     N_EV,WORK,LWORK,IWORK,LIWORK,INFO)
!
! -------------------------------------------------------------------------
! PDSYEVD_X works like PDSYEVD but has an additional parameter N_EV
! which gives the number of eigenvectors which are backtransformed.
! If not all eigenvectors are needed, this should save some computing time.
! The last eigenvectors (beyond N_EV) are not valid on return!
! -------------------------------------------------------------------------
!
!  -- ScaLAPACK routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     March 14, 2000
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N, N_EV
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), W( * ), WORK( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  PDSYEVD computes  all the eigenvalues and eigenvectors
!  of a real symmetric matrix A by calling the recommended sequence
!  of ScaLAPACK routines.
!
!  In its present form, PDSYEVD assumes a homogeneous system and makes
!  no checks for consistency of the eigenvalues or eigenvectors across
!  the different processes.  Because of this, it is possible that a
!  heterogeneous system may return incorrect results without any error
!  messages.
!
!  Arguments
!  =========
!
!     NP = the number of rows local to a given process.
!     NQ = the number of columns local to a given process.
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;     (NOT IMPLEMENTED YET)
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (global input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (global input) INTEGER
!          The number of rows and columns to be operated on, i.e. the
!          order of the distributed submatrix sub( A ). N >= 0.
!
!  A       (local input/workspace) block cyclic DOUBLE PRECISION array,
!          global dimension (N, N), local dimension ( LLD_A,
!          LOCc(JA+N-1) )
!          On entry, the symmetric matrix A.  If UPLO = 'U', only the
!          upper triangular part of A is used to define the elements of
!          the symmetric matrix.  If UPLO = 'L', only the lower
!          triangular part of A is used to define the elements of the
!          symmetric matrix.
!          On exit, the lower triangle (if UPLO='L') or the upper
!          triangle (if UPLO='U') of A, including the diagonal, is
!          destroyed.
!
!  IA      (global input) INTEGER
!          A's global row index, which points to the beginning of the
!          submatrix which is to be operated on.
!
!  JA      (global input) INTEGER
!          A's global column index, which points to the beginning of
!          the submatrix which is to be operated on.
!
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix A.
!
!  W       (global output) DOUBLE PRECISION array, dimension (N)
!          If INFO=0, the eigenvalues in ascending order.
!
!  Z       (local output) DOUBLE PRECISION array,
!          global dimension (N, N),
!          local dimension ( LLD_Z, LOCc(JZ+N-1) )
!          Z contains the orthonormal eigenvectors
!          of the symmetric matrix A.
!
!  IZ      (global input) INTEGER
!          Z's global row index, which points to the beginning of the
!          submatrix which is to be operated on.
!
!  JZ      (global input) INTEGER
!          Z's global column index, which points to the beginning of
!          the submatrix which is to be operated on.
!
!  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix Z.
!          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
!
!  WORK    (local workspace/output) DOUBLE PRECISION array,
!          dimension (LWORK)
!          On output, WORK(1) returns the workspace required.
!
!  LWORK   (local input) INTEGER
!          LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
!          TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
!          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
!          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
!
!          If LWORK = -1, the LWORK is global input and a workspace
!          query is assumed; the routine only calculates the minimum
!          size for the WORK array.  The required workspace is returned
!          as the first element of WORK and no error message is issued
!          by PXERBLA.
!
!  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
!          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          LIWORK = 7*N + 8*NPCOL + 2
!
!  INFO    (global output) INTEGER
!          = 0:  successful exit
!          < 0:  If the i-th argument is an array and the j-entry had
!                an illegal value, then INFO = -(i*100+j), if the i-th
!                argument is a scalar and had an illegal value, then
!                INFO = -i.
!          > 0:  The algorithm failed to compute the INFO/(N+1) th
!                eigenvalue while working on the submatrix lying in
!                global rows and columns mod(INFO,N+1).
!
!  Alignment requirements
!  ======================
!
!  The distributed submatrices sub( A ), sub( Z ) must verify
!  some alignment properties, namely the following expression
!  should be true:
!  ( MB_A.EQ.NB_A.EQ.MB_Z.EQ.NB_Z .AND. IROFFA.EQ.ICOFFA .AND.
!    IROFFA.EQ.0 .AND.IROFFA.EQ.IROFFZ. AND. IAROW.EQ.IZROW)
!    with IROFFA = MOD( IA-1, MB_A )
!     and ICOFFA = MOD( JA-1, NB_A ).
!
!  Further Details
!  ======= =======
!
!  Contributed by Francoise Tisseur, University of Manchester.
!
!  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
!              Conquer Algorithm for the Symmetric Eigenvalue Problem
!              on Distributed Memory Architectures",
!              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
!              (see also LAPACK Working Note 132)
!                http://www.netlib.org/lapack/lawns/lawn132.ps
!
!  =====================================================================
!
!     .. Parameters ..
!
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_, &
     &                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
     &                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
     &                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IACOL, IAROW, ICOFFA, ICOFFZ, ICTXT, IINFO,    &
     &                   INDD, INDE, INDE2, INDTAU, INDWORK, INDWORK2,  &
     &                   IROFFA, IROFFZ, ISCALE, LIWMIN, LLWORK,        &
     &                   LLWORK2, LWMIN, MYCOL, MYROW, NB, NP, NPCOL,   &
     &                   NPROW, NQ, OFFSET, TRILWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,  &
     &                   SMLNUM
!     ..
!     .. Local Arrays ..
!     ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANSY
      EXTERNAL           LSAME, INDXG2P, NUMROC, PDLAMCH, PDLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DSCAL, PCHK1MAT,      &
     &                   PDLARED1D, PDLASCL, PDLASET, PDORMTR, PDSTEDC, &
     &                   PDSYTRD, PXERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MIN, MOD, SQRT
!     ..
!     .. Executable Statements ..
!       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*  &
     &    RSRC_.LT.0 )RETURN
!
!     Quick return
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     Test the input arguments.
!
      CALL BLACS_GRIDINFO( DESCZ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IZ, JZ, DESCZ, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            NB = DESCA( NB_ )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFZ = MOD( IZ-1, DESCZ( MB_ ) )
            ICOFFZ = MOD( JZ-1, DESCZ( NB_ ) )
            IAROW = INDXG2P( IA, NB, MYROW, DESCA( RSRC_ ), NPROW )
            IACOL = INDXG2P( JA, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
!
            LQUERY = ( LWORK.EQ.-1 )
            TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
            LWMIN = MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
            LIWMIN = 7*N + 8*NPCOL + 2
            WORK( 1 ) = DBLE( LWMIN )
            IWORK( 1 ) = LIWMIN
            IF( .NOT.LSAME( JOBZ, 'V' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( IROFFA.NE.ICOFFA .OR. ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( IROFFA.NE.IROFFZ .OR. ICOFFA.NE.ICOFFZ ) THEN
               INFO = -10
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 1200+M_ )
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( DESCZ( MB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 1200+NB_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 1200+MB_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 1200+CTXT_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 1200+RSRC_ )
            ELSE IF( DESCA( CSRC_ ).NE.DESCZ( CSRC_ ) ) THEN
               INFO = -( 1200+CSRC_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -14
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -16
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, 2, IDUM1, IDUM2,  &
     &                  INFO )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Set up pointers into the WORK array
!
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDE2 = INDD + N
      INDWORK = INDE2 + N
      LLWORK = LWORK - INDWORK + 1
      INDWORK2 = INDD
      LLWORK2 = LWORK - INDWORK2 + 1
!
!     Scale matrix to allowable range, if necessary.
!
      ISCALE = 0
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
      ANRM = PDLANSY( 'M', UPLO, N, A, IA, JA, DESCA, WORK( INDWORK ) )
!
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
!
      IF( ISCALE.EQ.1 ) THEN
         CALL PDLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
!
!     Reduce symmetric matrix to tridiagonal form.
!
!
      CALL PDSYTRD( UPLO, N, A, IA, JA, DESCA, WORK( INDD ),            &
     &              WORK( INDE2 ), WORK( INDTAU ), WORK( INDWORK ),     &
     &              LLWORK, IINFO )
!
!     Copy the values of D, E to all processes.
!
      CALL PDLARED1D( N, IA, JA, DESCA, WORK( INDD ), W,                &
     &                WORK( INDWORK ), LLWORK )
!
      CALL PDLARED1D( N, IA, JA, DESCA, WORK( INDE2 ), WORK( INDE ),    &
     &                WORK( INDWORK ), LLWORK )
!
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
!
      IF( UPPER ) THEN
         OFFSET = 1
      ELSE
         OFFSET = 0
      END IF

      CALL PDSTEDC( 'I', N, W, WORK( INDE+OFFSET ), Z, IZ, JZ, DESCZ,   &
     &              WORK( INDWORK2 ), LLWORK2, IWORK, LIWORK, INFO )
!
      CALL PDORMTR( 'L', UPLO, 'N', N, N_EV, A, IA, JA, DESCA,          &
     &              WORK( INDTAU ), Z, IZ, JZ, DESCZ, WORK( INDWORK2 ), &
     &              LLWORK2, IINFO )
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
!
      RETURN
!
!     End of PDSYEVD
!
      END SUBROUTINE PDSYEVD_X

!-----------------------------------------------------------------------------------

      SUBROUTINE PZHEEVD_X(JOBZ,UPLO,N,A,IA,JA,DESCA,W,Z,IZ,JZ,DESCZ,   &
     &                    N_EV, WORK, LWORK, RWORK, LRWORK, IWORK,      &
     &                    LIWORK, INFO )
!
! -------------------------------------------------------------------------
! PZHEEVD_X works like PZHEEVD but has an additional parameter N_EV
! which gives the number of eigenvectors which are backtransformed.
! If not all eigenvectors are needed, this should save some computing time.
! The last eigenvectors (beyond N_EV) are not valid on return!
! -------------------------------------------------------------------------
!
!  -- ScaLAPACK routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     March 25, 2002
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA,INFO,IZ,JA,JZ,LIWORK,LRWORK,LWORK,N,N_EV
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( * ), WORK( * ), Z( * )
!
!
!  Purpose
!  =======
!
!  PZHEEVD computes all the eigenvalues and eigenvectors of a Hermitian
!  matrix A by using a divide and conquer algorithm.
!
!  Arguments
!  =========
!
!     NP = the number of rows local to a given process.
!     NQ = the number of columns local to a given process.
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;    (NOT IMPLEMENTED YET)
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (global input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (global input) INTEGER
!          The number of rows and columns of the matrix A.  N >= 0.
!
!  A       (local input/workspace) block cyclic COMPLEX*16 array,
!          global dimension (N, N), local dimension ( LLD_A,
!          LOCc(JA+N-1) )
!
!          On entry, the symmetric matrix A.  If UPLO = 'U', only the
!          upper triangular part of A is used to define the elements of
!          the symmetric matrix.  If UPLO = 'L', only the lower
!          triangular part of A is used to define the elements of the
!          symmetric matrix.
!
!          On exit, the lower triangle (if UPLO='L') or the upper
!          triangle (if UPLO='U') of A, including the diagonal, is
!          destroyed.
!
!  IA      (global input) INTEGER
!          A's global row index, which points to the beginning of the
!          submatrix which is to be operated on.
!
!  JA      (global input) INTEGER
!          A's global column index, which points to the beginning of
!          the submatrix which is to be operated on.
!
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix A.
!          If DESCA( CTXT_ ) is incorrect, PZHEEV cannot guarantee
!          correct error reporting.
!
!  W       (global output) DOUBLE PRECISION array, dimension (N)
!          If INFO=0, the eigenvalues in ascending order.
!
!  Z       (local output) COMPLEX*16 array,
!          global dimension (N, N),
!          local dimension ( LLD_Z, LOCc(JZ+N-1) )
!          Z contains the orthonormal eigenvectors of the matrix A.
!
!  IZ      (global input) INTEGER
!          Z's global row index, which points to the beginning of the
!          submatrix which is to be operated on.
!
!  JZ      (global input) INTEGER
!          Z's global column index, which points to the beginning of
!          the submatrix which is to be operated on.
!
!  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix Z.
!          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
!
!  WORK    (local workspace/output) COMPLEX*16 array,
!          dimension (LWORK)
!          On output, WORK(1) returns the workspace needed for the
!          computation.
!
!  LWORK   (local input) INTEGER
!          If eigenvectors are requested:
!            LWORK = N + ( NP0 + MQ0 + NB ) * NB,
!          with  NP0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPROW )
!                MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
!
!          If LWORK = -1, then LWORK is global input and a workspace
!          query is assumed; the routine calculates the size for all
!          work arrays. Each of these values is returned in the first
!          entry of the corresponding work array, and no error message
!          is issued by PXERBLA.
!
!  RWORK   (local workspace/output) DOUBLE PRECISION array,
!          dimension (LRWORK)
!          On output RWORK(1) returns the real workspace needed to
!          guarantee completion.  If the input parameters are incorrect,
!          RWORK(1) may also be incorrect.
!
!  LRWORK  (local input) INTEGER
!          Size of RWORK array.
!          LRWORK >= 1 + 9*N + 3*NP*NQ,
!          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
!          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
!
!  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
!          On output IWORK(1) returns the integer workspace needed.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          LIWORK = 7*N + 8*NPCOL + 2
!
!  INFO    (global output) INTEGER
!          = 0:  successful exit
!          < 0:  If the i-th argument is an array and the j-entry had
!                an illegal value, then INFO = -(i*100+j), if the i-th
!                argument is a scalar and had an illegal value, then
!                INFO = -i.
!          > 0:  If INFO = 1 through N, the i(th) eigenvalue did not
!                converge in PDLAED3.
!
!  Alignment requirements
!  ======================
!
!  The distributed submatrices sub( A ), sub( Z ) must verify
!  some alignment properties, namely the following expression
!  should be true:
!  ( MB_A.EQ.NB_A.EQ.MB_Z.EQ.NB_Z .AND. IROFFA.EQ.ICOFFA .AND.
!    IROFFA.EQ.0 .AND.IROFFA.EQ.IROFFZ. AND. IAROW.EQ.IZROW)
!    with IROFFA = MOD( IA-1, MB_A )
!     and ICOFFA = MOD( JA-1, NB_A ).
!
!  Further Details
!  ======= =======
!
!  Contributed by Francoise Tisseur, University of Manchester.
!
!  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
!              Conquer Algorithm for the Symmetric Eigenvalue Problem
!              on Distributed Memory Architectures",
!              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
!              (see also LAPACK Working Note 132)
!                http://www.netlib.org/lapack/lawns/lawn132.ps
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_, &
     &                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
     &                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
     &                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION               ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),                  &
     &                   CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY
      INTEGER            CSRC_A, I, IACOL, IAROW, ICOFFA, IINFO, IIZ,   &
     &                   INDD, INDE, INDE2, INDRWORK, INDTAU, INDWORK,  &
     &                   INDZ, IPR, IPZ, IROFFA, IROFFZ, ISCALE, IZCOL, &
     &                   IZROW, J, JJZ, LDR, LDZ, LIWMIN, LLRWORK,      &
     &                   LLWORK, LRWMIN, LWMIN, MB_A, MYCOL, MYROW, NB, &
     &                   NB_A, NN, NP0, NPCOL, NPROW, NQ, NQ0, OFFSET,  &
     &                   RSRC_A
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,  &
     &                   SMLNUM
!     ..
!     .. Local Arrays ..
      INTEGER            DESCRZ( 9 ), IDUM1( 2 ), IDUM2( 2 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC
      DOUBLE PRECISION   PZLANHE, PDLAMCH
      EXTERNAL           LSAME, INDXG2L, INDXG2P, NUMROC, PZLANHE,      &
     &                   PDLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCINIT, INFOG2L,    &
     &                   PZELGET, PZHETRD, PCHK2MAT, PZLASCL, PZLASET,  &
     &                   PZUNMTR, PDLARED1D, PDLASET, PDSTEDC, PXERBLA, &
     &                   DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, ICHAR, MAX, MIN, MOD, DBLE, SQRT
!     ..
!     .. Executable Statements ..
!       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*  &
     &    RSRC_.LT.0 )RETURN
!
      INFO = 0
!
!     Quick return
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     Test the input arguments.
!
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         CALL CHK1MAT( N, 2, N, 2, IZ, JZ, DESCZ, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            LOWER = LSAME( UPLO, 'L' )
            NB_A = DESCA( NB_ )
            MB_A = DESCA( MB_ )
            NB = NB_A
            RSRC_A = DESCA( RSRC_ )
            CSRC_A = DESCA( CSRC_ )
            IROFFA = MOD( IA-1, MB_A )
            ICOFFA = MOD( JA-1, NB_A )
            IAROW = INDXG2P( IA, NB_A, MYROW, RSRC_A, NPROW )
            IACOL = INDXG2P( JA, MB_A, MYCOL, CSRC_A, NPCOL )
            NP0 = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ0 = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
            IROFFZ = MOD( IZ-1, MB_A )
            CALL INFOG2L( IZ, JZ, DESCZ, NPROW, NPCOL, MYROW, MYCOL,    &
     &                    IIZ, JJZ, IZROW, IZCOL )
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
!
!           Compute the total amount of space needed
!
            NN = MAX( N, NB, 2 )
            NQ = NUMROC( NN, NB, 0, 0, NPCOL )
            LWMIN = N + ( NP0+NQ+NB )*NB
            LRWMIN = 1 + 9*N + 3*NP0*NQ0
            LIWMIN = 7*N + 8*NPCOL + 2
            WORK( 1 ) = DCMPLX( LWMIN )
            RWORK( 1 ) = DBLE( LRWMIN )
            IWORK( 1 ) = LIWMIN
            IF( .NOT.LSAME( JOBZ, 'V' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
               INFO = -14
            ELSE IF( LRWORK.LT.LRWMIN .AND. LRWORK.NE.-1 ) THEN
               INFO = -16
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -4
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( IROFFA.NE.IROFFZ ) THEN
               INFO = -10
            ELSE IF( IAROW.NE.IZROW ) THEN
               INFO = -10
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 1200+M_ )
            ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
               INFO = -( 1200+N_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 1200+MB_ )
            ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 1200+NB_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 1200+RSRC_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 1200+CTXT_ )
            END IF
         END IF
         IF( LOWER ) THEN
            IDUM1( 1 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'U' )
         END IF
         IDUM2( 1 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IZ,   &
     &                  JZ, DESCZ, 11, 2, IDUM1, IDUM2, INFO )
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PZHEEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Get machine constants.
!
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
!
!     Set up pointers into the WORK array
!
      INDTAU = 1
      INDWORK = INDTAU + N
      LLWORK = LWORK - INDWORK + 1
!
!     Set up pointers into the RWORK array
!
      INDE = 1
      INDD = INDE + N
      INDE2 = INDD + N
      INDRWORK = INDE2 + N
      LLRWORK = LRWORK - INDRWORK + 1
!
!     Scale matrix to allowable range, if necessary.
!
      ISCALE = 0
!
      ANRM = PZLANHE( 'M', UPLO, N, A, IA, JA, DESCA,                   &
     &       RWORK( INDRWORK ) )
!
!
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
!
      IF( ISCALE.EQ.1 ) THEN
         CALL PZLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
!
!     Reduce Hermitian matrix to tridiagonal form.
!
      CALL PZHETRD( UPLO, N, A, IA, JA, DESCA, RWORK( INDD ),           &
     &              RWORK( INDE2 ), WORK( INDTAU ), WORK( INDWORK ),    &
     &              LLWORK, IINFO )
!
!     Copy the values of D, E to all processes
!
!     Here PxLARED1D is used to redistribute the tridiagonal matrix.
!     PxLARED1D, however, doesn't yet workMx Mawith arbritary matrix
!     distributions so we have PxELGET as a backup.
!
      OFFSET = 0
      IF( IA.EQ.1 .AND. JA.EQ.1 .AND. RSRC_A.EQ.0 .AND. CSRC_A.EQ.0 )   &
     &     THEN
         CALL PDLARED1D( N, IA, JA, DESCA, RWORK( INDD ), W,            &
     &                   RWORK( INDRWORK ), LLRWORK )
!
         CALL PDLARED1D( N, IA, JA, DESCA, RWORK( INDE2 ),              &
     &                   RWORK( INDE ), RWORK( INDRWORK ), LLRWORK )
         IF( .NOT.LOWER )                                               &
     &      OFFSET = 1
      ELSE
         DO 10 I = 1, N
            CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA-1, I+JA-1, &
     &                    DESCA )
            W( I ) = DBLE( WORK( INDWORK ) )
   10    CONTINUE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 I = 1, N - 1
               CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA-1, I+JA,&
     &                       DESCA )
               RWORK( INDE+I-1 ) = DBLE( WORK( INDWORK ) )
   20       CONTINUE
         ELSE
            DO 30 I = 1, N - 1
               CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA, I+JA-1,&
     &                       DESCA )
               RWORK( INDE+I-1 ) = DBLE( WORK( INDWORK ) )
   30       CONTINUE
         END IF
      END IF
!
!     Call PDSTEDC to compute eigenvalues and eigenvectors.
!
      INDZ = INDE + N
      INDRWORK = INDZ + NP0*NQ0
      LLRWORK = LRWORK - INDRWORK + 1
      LDR = MAX( 1, NP0 )
      CALL DESCINIT( DESCRZ, DESCZ( M_ ), DESCZ( N_ ), DESCZ( MB_ ),    &
     &               DESCZ( NB_ ), DESCZ( RSRC_ ), DESCZ( CSRC_ ),      &
     &               DESCZ( CTXT_ ), LDR, INFO )
      CALL PZLASET( 'Full', N, N, CZERO, CONE, Z, IZ, JZ, DESCZ )
      CALL PDLASET( 'Full', N, N, ZERO, ONE, RWORK( INDZ ), 1, 1,       &
     &              DESCRZ )
      CALL PDSTEDC( 'I', N, W, RWORK( INDE+OFFSET ), RWORK( INDZ ), IZ, &
     &              JZ, DESCRZ, RWORK( INDRWORK ), LLRWORK, IWORK,      &
     &              LIWORK, IINFO )
!
      LDZ = DESCZ( LLD_ )
      LDR = DESCRZ( LLD_ )
      IIZ = INDXG2L( IZ, NB, MYROW, MYROW, NPROW )
      JJZ = INDXG2L( JZ, NB, MYCOL, MYCOL, NPCOL )
      IPZ = IIZ + ( JJZ-1 )*LDZ
      IPR = INDZ - 1 + IIZ + ( JJZ-1 )*LDR
      DO 50 J = 0, NQ0 - 1
         DO 40 I = 0, NP0 - 1
            Z( IPZ+I+J*LDZ ) = RWORK( IPR+I+J*LDR )
   40    CONTINUE
   50 CONTINUE
!
!     Z = Q * Z
!
      CALL PZUNMTR( 'L', UPLO, 'N', N, N_EV, A, IA, JA, DESCA,          &
     &              WORK( INDTAU ), Z, IZ, JZ, DESCZ, WORK( INDWORK ),  &
     &              LLWORK, IINFO )
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
!
      WORK( 1 ) = DCMPLX( LWMIN )
      RWORK( 1 ) = DBLE( LRWMIN )
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of PZHEEVD
!
      END SUBROUTINE PZHEEVD_X

