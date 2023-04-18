!***** FHI-aims/prodbas
!  NAME
!    prodbas
!  SYNOPSIS

module prodbas

! PURPOSE
!  Module prodbas contains the information about the auxiliary (product) basis.
!  This auxiliary are used to expand the quantities in the HF, MP2, GW calculation,
!  including polarisability, bare and screened coulomb interaction, and self energy .
!
!  This module containing subroutines:
!  * allocate_basbas
!  * prodas_task_distribution
!  * allocate_basis_pairs
!  * cleanup_basbas
!

!  global variable declarations:
!  *   n_basbas        : total number of product basis functions
!  *   n_loc_prodbas   : number of auxiliary basis functions on this thread.
!  *   n_basbas_supercell   : number of product basis functions for periodic systems
!  *   map_prodbas     : Maps the auxiliary basis on each thread to the original auxiliary basis
!                        Size: (1:n_max_loc_prodbas,1:n_task)
!  *   n_basbas_fns    : number of different radial functions of auxiary basis (< n_basbas)
!  *   basbas_wave_spl : Basis functions, cubic splined on the logarithmic grid
!                         [e - v_basis(r)] * u_basis(r) which is also: u"(r) + [l*(l+1)/r^2]*u(r)
!                        Size: (n_max_spline(==4), n_max_grid, n_basbas_fns)
!    The following index fields are of size (n_basbas):
!  *   basbas_atom : atom to which a basis fn corresponds
!  *   basbas_l    : angular momentum quantum number of given basis function
!  *   basbas_m    : magnetic quantum number of given basis function
!  *   basbas_fn   : The radial basis function number which corresponds to
!                   basis function i_basbas [There are more basis functions (at,type,n,l,m)
!                   than there are necessary radial functions (at,type, n,l).]:
!                   1 <= basbas_fn(i_basbas) <= n_basbas_fns.
!  *  basbasfn_species : i_basbas_fn -> i_species
!  *  basbasfn_l       : i_basbas_fn -> i_l
!
!  *   charge_radius_basbas_fn(i_basbas_fn) :
!                       Radius at which the radial part has a
!                       magnitude of wave_threshold (~ 1d-6).
!            The overlap of two product basis functions can be neglected if
!            the separation is larger than the sum of their charge_radii.
!  *   field_radius_basbas_fn(i_basbas_fn) :
!                       Radius at which the Coulomb potential has a magnitude
!                       of wave_threshold [or wave_threshold**2 for finite
!                       multipole].
!            The Coulomb interaction of two product basis functions can be
!            neglected if the separation is larger than any of the two
!            distinct sums of charge_radius and outer_radius.  That is,
!            separated product basis functions only Coulomb interact if both
!            of them have finite multipole moments!
!  *   multipole_basbas_fn(i_basbas_fn) :
!                       Multipole moment of a given product basis function.
!
!  *  Some comments by RJ about the 2D distribution of ovlp_3KS:
!        The dimension is now ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
!        i.e. the second and third dimension is distributed onto a 2D grid.
!        The actual distribution is a simple cyclic distribution (not block
!        cyclic) [...]  but the distribution shouldn't matter at all since I
!        tried to hide it using the arrays:
!           own_dim1_o3ks(:)
!           loc_dim1_o3ks(:)
!           own_dim2_o3ks(:)
!           loc_dim2_o3ks(:)
!        own_dim1_o3ks(i) has the "owner" of the global index i of the first
!        distributed dimension, loc_dim1_o3ks(i) has the local index for this
!        global index.
!        Same for own_dim2_o3ks/loc_dim2_o3ks.
!        So the global elements ovlp_3KS(:,i,j,x) are located on processor
!        with coords (own_dim1_o3ks(i),own_dim2_o3ks(j)) and are stored there
!        in ovlp_3KS(:,loc_dim1_o3ks(i),loc_dim2_o3ks(j),x)
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!
!  SOURCE

      implicit none


      integer, dimension(:,:), allocatable :: map_prodbas

      real*8, dimension(:,:,:), allocatable :: basbas_wave_spl

      integer, dimension(:), allocatable :: basbas_atom
      integer, dimension(:), allocatable :: basbas_l
      integer, dimension(:), allocatable :: basbas_m
      integer, dimension(:), allocatable :: basbas_fn

      integer, dimension(:), allocatable :: basbasfn_species
      integer, dimension(:), allocatable :: basbasfn_l

      real*8, dimension(:), allocatable :: field_radius_basbas_fn
      real*8, dimension(:), allocatable :: charge_radius_basbas_fn
      real*8, dimension(:), allocatable :: multipole_basbas_fn

      ! i_loc_prodbas -> field_radius
      real*8, dimension(:), allocatable :: outer_radius_prodbas

      ! JW: Arrays to loop basis functions by their angular momentum.  See
      ! basis.f90 for documentation.
      integer :: max_basbas_L
      integer :: max_n_basbas_fnLsp
      integer, allocatable :: Lsp2n_basbas_fnLsp(:,:)
      integer, allocatable :: Lsp2basbas_fn(:,:,:)
      integer, allocatable :: Lsp2basbas_sp(:,:,:)
      integer, allocatable :: atom2basbas_off(:)
      integer, allocatable :: sp2n_basbas_sp(:)
      integer :: max_n_basbas_sp


      integer, parameter :: OVLP_TYPE_OVERLAP = 1
      integer, parameter :: OVLP_TYPE_COULOMB = 2
      integer, parameter :: OVLP_TYPE_HSE = 3
      integer, parameter :: OVLP_TYPE_CUT = 4
      integer, parameter :: OVLP_TYPE_CUT_ANALYTIC = 5
      integer, parameter :: OVLP_TYPE_ERS = 6
      integer, parameter :: OVLP_TYPE_LR = 7
      character*7, parameter :: OVLP_TYPE_NAMES(7) = (/'overlap', 'Coulomb', &
      &                                                'HSE-pot', 'cutCoul', &
      &                                                'Alavi-C', 'ERS-pot', &
      																 'LRC-pot' /)
      ! Either OVLP_TYPE_COULOMB or OVLP_TYPE_HSE.
      integer :: ovlp_type_bare_or_hse_coul
      ! Either OVLP_TYPE_COULOMB or OVLP_TYPE_ERS
      integer :: ovlp_type_bare_or_ers_coul

      integer :: n_max_basis_atom
      integer :: n_basis_pairs
      integer, dimension(:), allocatable ::  n_nghr_basis   ! See condense_basis_pairs()
      integer, dimension(:,:), allocatable :: basis_nghr    ! for documentation.
      integer, dimension(:), allocatable ::  n_nghr_basis_original  ! atom bsse 
      integer, dimension(:,:), allocatable :: basis_nghr_original   ! atom bsse
      integer, dimension(:,:), allocatable :: atom_pairs

      real*8, dimension(:,:,:), allocatable :: v_times_radialbasbas_spl

!   Variables to ensure we never call prodbas_tasks_distribution more than necessary
!   e.g. during a relaxation. (New MPI communicators are created there, and they
!   should not normally be structure dependent - but MPI libraries tend to give up
!   if one creates too many MPI communicators.)

      logical :: prodbas_tasks_distributed
      integer :: previous_n_basbas
      integer :: previous_prodbas_nb
      integer :: previous_n_states

      ! ... and similar for prodbas_tasks_distribution_pbc
      logical :: prodbas_tasks_pbc_distributed
      integer :: previous_n_basbas_supercell

!   scalapack for auxiliary basis
      integer :: my_blacs_ctxt_aux
      integer :: nprow_aux, npcol_aux
      integer :: myprow_aux, mypcol_aux
      integer :: mb_aux, nb_aux
      integer :: max_row, max_col

!   scalapack for auxiliary basis, 2D task distribution
      integer :: my_blacs_ctxt_aux_2d
      integer :: nprow_aux_2d, npcol_aux_2d
      integer :: myprow_aux_2d, mypcol_aux_2d
      integer :: nb_aux_2d
      integer :: max_row_2d, max_col_2d

!   mapping of 2D processor coords to global ids
      integer, dimension(:,:), allocatable :: global_id

!   communicators for processor rows/columns
      integer :: mpi_comm_rows_aux_2d
      integer :: mpi_comm_cols_aux_2d

      ! This is only a fallback variable in case the prodbas_tasks_distribution
      ! is called more than once during a run
      logical, private :: prodbas_distributed

!   Scalapack descriptors for 1D and 2D distributed matrices

      integer :: aux_sc_desc(9)
      integer :: aux_sc_desc_2d(9)

!   grid for 2D distribution of ovlp_3KS matrix

      integer :: np1_o3KS
      integer :: np2_o3KS
      integer :: myp1_o3KS
      integer :: myp2_o3KS

!   dimensions of 2D distributed ovlp_3KS; set when they are actually known

      integer :: ndim1_o3KS
      integer :: ndim2_o3KS

!   owner and local index for distributed ovlp_3KS

      integer, allocatable :: own_dim1_o3ks(:)
      integer, allocatable :: loc_dim1_o3ks(:)
      integer, allocatable :: own_dim2_o3ks(:)
      integer, allocatable :: loc_dim2_o3ks(:)

!   MPI communicators for 2D distributed ovlp_3KS

      integer :: mpi_comm_o3ks_1 ! communicator for all procs with identical myp1_o3KS
      integer :: mpi_comm_o3ks_2 ! communicator for all procs with identical myp2_o3KS

      ! This is only a fallback variable in case the prodbas_tasks_distribution
      ! is called more than once during a run
      logical, private :: o3ks_distributed

 ! SVL scalapack for auxiliary basis for entire supercell
      integer :: sc_my_blacs_ctxt_aux
      integer :: sc_nprow_aux, sc_npcol_aux
      integer :: sc_myprow_aux, sc_mypcol_aux
      integer :: sc_mb_aux, sc_nb_aux
      integer :: sc_max_row, sc_max_col


!******

      contains
!---------------------------------------------------------------------
!****s* prodbas/allocate_basbas
!  NAME
!   allocate_basbas
!  SYNOPSIS

        subroutine allocate_basbas &
        ( )

!  PURPOSE
!  Subroutine allocate_basis allocates all basis-related arrays, after
!  the relevant dimensions n_basis, n_basis_fn are known.
!  USES

        use runtime_choices
        use dimensions
        use mpi_tasks, only: check_allocation
        implicit none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

        integer:: info
        character(*), parameter :: func = 'allocate_basbas'

        if(.not.allocated(basbas_wave_spl)) then
          allocate &
          ( basbas_wave_spl ( &
                  n_max_spline,n_max_grid,n_basbas_fns),stat=info )
          call check_allocation(info, 'basbas_wave_spl               ')
        endif

        if(.not.allocated(v_times_radialbasbas_spl)) then
          allocate &
          ( v_times_radialbasbas_spl ( &
                  n_max_spline,n_hartree_grid,n_basbas_fns),stat=info )
          call check_allocation(info, 'v_times_radialbasbas_spl               ')
        endif

        if(.not.allocated(basbas_atom)) then
          allocate (basbas_atom(n_basbas_supercell),stat=info)
          call check_allocation(info, 'basbas_atom                   ')
        endif
        if(.not.allocated(basbas_l)) then
          allocate (basbas_l(n_basbas_supercell),stat=info)
          call check_allocation(info, 'basbas_l                      ')
        endif
        if(.not.allocated(basbas_m)) then
          allocate (basbas_m(n_basbas_supercell),stat=info)
          call check_allocation(info, 'basbas_m                      ')
        endif
        if(.not.allocated(basbas_fn)) then
          allocate (basbas_fn(n_basbas_supercell),stat=info)
          call check_allocation(info, 'basbas_fn                     ')
        endif

        if(.not.allocated(basbasfn_species)) then
           allocate (basbasfn_species(n_basbas_fns),stat=info)
           call check_allocation(info, 'basbasfn_species', func)
        endif
        if(.not.allocated(basbasfn_l)) then
           allocate (basbasfn_l(n_basbas_fns),stat=info)
           call check_allocation(info, 'basbasfn_l', func)
        endif

        if(.not.allocated(charge_radius_basbas_fn)) then
           allocate(charge_radius_basbas_fn(n_basbas_fns), stat=info)
           call check_allocation(info, 'charge_radius_basbas_fn', func)
        end if
        if(.not.allocated(field_radius_basbas_fn)) then
           allocate(field_radius_basbas_fn(n_basbas_fns), stat=info)
           call check_allocation(info, 'field_radius_basbas_fn', func)
        end if
        if(.not.allocated(multipole_basbas_fn)) then
           allocate(multipole_basbas_fn(n_basbas_fns), stat=info)
           call check_allocation(info, 'multipole_basbas_fn', func)
        end if

        if (.not. allocated(Lsp2n_basbas_fnLsp)) then
           allocate(Lsp2n_basbas_fnLsp(0:max_basbas_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2n_basbas_fnLsp', func)
        end if
        if (.not. allocated(Lsp2basbas_fn)) then
           allocate(Lsp2basbas_fn(max_n_basbas_fnLsp, 0:max_basbas_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2basbas_fn', func)
        end if
        if (.not. allocated(Lsp2basbas_sp)) then
           allocate(Lsp2basbas_sp(max_n_basbas_fnLsp, 0:max_basbas_L, n_species), stat=info)
           call check_allocation(info, 'Lsp2basbas_sp', func)
        end if
        if (.not. allocated(atom2basbas_off)) then
           allocate(atom2basbas_off(n_atoms), stat=info)
           call check_allocation(info, 'atom2basbas_off', func)
        end if
        if (.not. allocated(sp2n_basbas_sp)) then
           allocate(sp2n_basbas_sp(n_species), stat=info)
           call check_allocation(info, 'sp2n_basbas_sp', func)
        end if

        if (use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse &
            ) then
           ovlp_type_bare_or_hse_coul = OVLP_TYPE_HSE
        else if (lrc_pt2_started) then
           ovlp_type_bare_or_hse_coul = OVLP_TYPE_LR
        else
           ovlp_type_bare_or_hse_coul = OVLP_TYPE_COULOMB
        end if

        if (use_ers) then
           ovlp_type_bare_or_ers_coul = OVLP_TYPE_ERS
        elseif (use_erfc) then
           ovlp_type_bare_or_ers_coul = OVLP_TYPE_HSE
        end if

        end subroutine allocate_basbas


!******
  !----------------------------------------------------------------------------
  !****s* prodbas/initialize_prodbas
  !  NAME
  !    initialize_prodbas
  !  SYNOPSIS
 
  subroutine initialize_prodbas()

    !  PURPOSE
    !    Initializes the product basis and the basis pair arrays.
    !  USES

    use runtime_choices
    use dimensions
    use basis, only : cleanup_ext_basis
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'initialize_prodbas'
    

!    if(myid.eq.0) write(use_unit,*) " Initialization flags:   " 
!    if(myid.eq.0) write(use_unit,*) " flag_auxil_basis ", flag_auxil_basis
!    if(myid.eq.0) write(use_unit,*) " PRODBAS_OPT      ", PRODBAS_OPT
!    if(myid.eq.0) write(use_unit,*) " PRODBAS_FULL     ", PRODBAS_FULL
!    if(myid.eq.0) write(use_unit,*) " PRODBAS_SVD      ", PRODBAS_SVD 

    select case (flag_auxil_basis)
    case(PRODBAS_OPT)
       call shrink_opt_auxil_basis()
    case(PRODBAS_FULL)
!       if(use_aux_basis) then
       call shrink_full_auxil_basis_aux()
!       call shrink_full_auxil_basis()
!       else 
!           call shrink_full_auxil_basis()
!       end if
! 

!  VB:  The following cleanup has been disabled since the s.c.f. cycle may now restart
!       itself (reinitialize_scf) even for RI-V. If cleanup_ext_basis is performed,
!       the code will crash in reinitialize_scf.
!       However, the fact that this added memory is there throughout the scf cycle is 
!       not great. It adds about 5 MB per atom and task in terms of memory usage.
!       It would be better to always deallocate this memory and simply reallocate
!       and reconstruct it when needed. This amounts to a renewed call to the appropriate
!       shrink_fixed_basis_...* routines, though (I did not test what happens in that case).
!       This modification would be good both for RI-V and for RI-LVL.

!       Basbas basis has been constructed and we can free memory by deallocating 
!       extended basis set as it is not used anymore
!        if(.not. (use_hse .and. (use_gw .or. use_dftpt2))) call cleanup_ext_basis()


!
!  
    case(PRODBAS_SVD)
       call shrink_svd_auxil_basis()
    case default
       call aims_stop("This choice of auxiliary basis is NOT allowed!", func)
    end select

    call allocate_basis_pairs()
    call condense_basis_pairs()

  end subroutine initialize_prodbas
  !******
!--------------------------------------------------------------------
!****s* prodbas/prodbas_tasks_distribution_pbc
!  NAME
!   prodbas_tasks_distribution_pbc
!  SYNOPSIS

      subroutine prodbas_tasks_distribution_pbc &
                 ( )

!  PURPOSE
!  the auxiliary(product) basis is equally distributed among the available threads.
!  n_loc_prodbas_supercell:  the largest number of the amount of auxiliary basis on each thread.
!     map_product(1:n_loc_prodbas_supercell,1:n_task)
!  maps the auxiliary basis on each thread to the original auxiliary basis
!

!  USES

      use dimensions
      use localorb_io, only: localorb_allinfo, localorb_info, use_unit
      use mpi_tasks
      use runtime_choices
      use synchronize_mpi
      implicit none


!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer n_tmp
!     local counter
      integer :: i_task
      integer :: i_basis, i_loc_prodbas, i_basbas, i_basbas_fn
      integer :: i_index
      integer :: n_remain_prodbas

      character*150 :: info_str
      integer, external :: numroc
      character(*), parameter :: func = 'prodbas_tasks_distribution_pbc'

!     begin tasks

! test
!      use_scalapack = .true.
! end test
      if( use_scalapack ) then

!  a special choice of the processor grid, not most efficient, but is proper
!  for our task

        sc_nprow_aux = 1
        sc_npcol_aux = n_tasks

        sc_mb_aux = 16
        sc_nb_aux = 16

        if(sc_nb_aux*sc_npcol_aux .gt. n_basbas_supercell ) then
          sc_nb_aux = n_basbas_supercell/sc_npcol_aux
          sc_mb_aux = sc_nb_aux
        endif

! If this is not the first time we are distributing the tasks, the BLACS context below exists already. 
! We should free it first to avoid creating plenty of excess MPI communicators.
        if (prodbas_tasks_pbc_distributed) then
           call BLACS_Gridexit ( sc_my_blacs_ctxt_aux )
        end if

! get the context handle for auxiliary basis distribution
        call BLACS_GET (0,0,sc_my_blacs_ctxt_aux)
        call BLACS_Gridinit( sc_my_blacs_ctxt_aux, 'R', sc_nprow_aux, sc_npcol_aux )
        call BLACS_Gridinfo( sc_my_blacs_ctxt_aux, sc_nprow_aux, sc_npcol_aux, &
                             sc_myprow_aux, sc_mypcol_aux )
!        myprow_aux = 0
!        mypcol_aux = myid
        sc_max_row = numroc( n_basbas_supercell, sc_mb_aux, sc_myprow_aux, 0, sc_nprow_aux )
        sc_max_col = numroc( n_basbas_supercell, sc_nb_aux, sc_mypcol_aux, 0, sc_npcol_aux )

        n_tmp = sc_npcol_aux*sc_nb_aux
        n_tmp = mod(n_basbas_supercell, sc_npcol_aux*sc_nb_aux)

        if(n_tmp .gt. sc_nb_aux) n_tmp = sc_nb_aux

        n_loc_prodbas_supercell = sc_max_col

        n_max_loc_prodbas_supercell = n_basbas_supercell/(sc_npcol_aux*sc_nb_aux)*sc_nb_aux + n_tmp

        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A,2I6)') &
             "| Maximal number of product basis functions per thread :", &
              n_loc_prodbas_supercell, n_loc_prodbas_supercell
          write(use_unit,*)
!          write(use_unit,'(2X,A,f12.3,A)') &
!             "| Minimal requirement for computer memories :", &
!              (dble(n_max_loc_prodbas_supercell*n_tasks)*dble(n_basis*(n_basis+1)/2) &
!               +2.0*n_basbas_supercell*n_loc_prodbas_supercell*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3, &
!               " Gbs"
        endif
        write(info_str,"(2X,'| ',A,I6,' : ',I6)") &
        & "Number of local auxiliary basis funcs on node", &
        & myid, n_loc_prodbas_supercell
        call localorb_allinfo(info_str)

      else

         if(n_tasks.ge.n_basbas_supercell)then
            if(myid.le.n_basbas_supercell-1)then
               n_loc_prodbas_supercell = 1
            else
               n_loc_prodbas_supercell = 0
            endif
            n_max_loc_prodbas_supercell = 1
            n_remain_prodbas = 0
         else
            n_remain_prodbas = MOD(n_basbas_supercell,n_tasks)

            n_loc_prodbas_supercell = n_basbas_supercell/n_tasks
            n_max_loc_prodbas_supercell = n_loc_prodbas_supercell

            if (n_remain_prodbas.gt.0) then
               if (myid.le.n_remain_prodbas-1) n_loc_prodbas_supercell = n_loc_prodbas_supercell + 1
               n_max_loc_prodbas_supercell = n_max_loc_prodbas_supercell + 1
            endif
         endif

        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A,I6)') &
             "| Number of product basis functions per thread :", &
              n_loc_prodbas_supercell

!          write(use_unit,'(2X,A,f16.3,A)') &
!             "| Minimal requirement for computer memories :", &
!              (dble(n_loc_prodbas_supercell*n_tasks)*dble(n_basis*(n_basis+1)/2) &
!               +n_basbas_supercell*n_loc_prodbas_supercell*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3, &
!              " Gbs"
        endif

! end if scalapack
      endif

! SVL
      if (.not.allocated(n_prodbas_per_proc)) then
         allocate (n_prodbas_per_proc(n_tasks),stat=i_index)
         call check_allocation(i_index, 'n_prodbas_per_proc        ')
      endif
      n_prodbas_per_proc = 0
      n_prodbas_per_proc(myid+1) = n_loc_prodbas_supercell
      call sync_integer_vector(n_prodbas_per_proc,n_tasks)

! test
!      use_scalapack = .false.
! end test

      end subroutine prodbas_tasks_distribution_pbc
!******
!--------------------------------------------------------------------
!****s* prodbas/prodbas_tasks_distribution
!  NAME
!   prodbas_tasks_distribution
!  SYNOPSIS

      subroutine prodbas_tasks_distribution &
                 ( )

!  PURPOSE
!  the auxiliary(product) basis is equally distributed among the available threads.
!     map_prodbas(1:n_max_loc_prodbas,1:n_task)
!  maps the auxiliary basis on each thread to the original auxiliary basis
!

!  USES

      use dimensions
      use mpi_tasks
      use localorb_io, only: localorb_info, localorb_allinfo, use_unit, OL_norm
      use runtime_choices
      use synchronize_mpi
      use scalapack_utils
      implicit none
      
      
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer n_tmp 
!     local counter
      integer :: i_task
      integer :: i_basis
      integer :: i_index
      integer :: n_remain_prodbas
      integer :: info
      integer :: mpierr
      integer :: i_loc_prodbas, i_basbas, i_basbas_fn
      character*150 :: info_str
      integer, external :: numroc
      character(*), parameter :: func = 'prodbas_tasks_distribution'

!     begin tasks

      if( use_scalapack ) then
!  a special choice of the processor grid, not most efficient, but is proper
!  for our task

        nprow_aux = 1
        npcol_aux = n_tasks

        if (prodbas_nb > 0) then
           mb_aux = prodbas_nb
           nb_aux = prodbas_nb
        else
           mb_aux = 16
           nb_aux = 16
        end if

        if(nb_aux*npcol_aux .gt. n_basbas ) then
           if (prodbas_nb > 0) then
              ! If I remember correctly, it would not even work. -- JW
              call localorb_info('*** prodbas_nb too large; reducing.')
           end if
           nb_aux = n_basbas/npcol_aux
           mb_aux = nb_aux

           if(nb_aux == 0) then
              write(use_unit,*) "You have chosen a value for npcol_aux (", npcol_aux, ") which at the time  &
              &of the writing of this error message is equivalent to the total number of tasks  (", &
              n_tasks, ").  The value of npcol_aux is greater than the number of product basis elements (", & 
              n_basbas, ").  To avoid this error message, reduce the number of tasks spawned below &
              &the number of basis elements."
              call aims_stop('System is to small for scalapack (1d) - use fewer CPUs?')
           end if
        endif

! If this is not the first time we are distributing the tasks, the BLACS context below exists already. 
! We should free it first to avoid creating plenty of excess MPI communicators.
        if (prodbas_tasks_distributed) then
           call BLACS_Gridexit ( my_blacs_ctxt_aux )
        end if

! get the context handle for auxiliary basis distribution
        call BLACS_GET (0,0,my_blacs_ctxt_aux)
        call BLACS_Gridinit( my_blacs_ctxt_aux, 'R', nprow_aux, npcol_aux )
        call BLACS_Gridinfo( my_blacs_ctxt_aux, nprow_aux, npcol_aux, &
                             myprow_aux, mypcol_aux )
        max_row = numroc( n_basbas, mb_aux, myprow_aux, 0, nprow_aux )
        max_col = numroc( n_basbas, nb_aux, mypcol_aux, 0, npcol_aux )

        n_tmp = npcol_aux*nb_aux
        n_tmp = mod(n_basbas, npcol_aux*nb_aux)

        if(n_tmp .gt. nb_aux) n_tmp = nb_aux

        n_loc_prodbas = max_col

        n_max_loc_prodbas = n_basbas/(npcol_aux*nb_aux)*nb_aux + n_tmp

        ! If we get here and map_prodbas is already allocate, the task distribution
        ! must have changed. Deallocate and reallocate.
        if (allocated(map_prodbas)) then
           deallocate (map_prodbas)
        endif
        allocate(map_prodbas(n_max_loc_prodbas,n_tasks),stat=i_index)
        call check_allocation(i_index, 'map_prodbas                   ')

        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A,2I6)') &
             "| Maximal number of product basis functions per thread :", &
              n_loc_prodbas, n_loc_prodbas
          write(use_unit,*)
          if(.not.use_lc_wpbeh .or. hybrid_coeff .eq. 0.0d0) then
           write(use_unit,'(2X,A,f12.3,A)') &
              "| Minimal total memory requirement for three-center overlap :", &
               (dble(n_max_loc_prodbas*n_tasks)*dble(n_basis*(n_basis+1)/2) &
                +2.0*n_basbas*n_loc_prodbas*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3, &
                " Gbs"
          else if (use_lc_wpbeh .and. hybrid_coeff .ne. 0.0d0) then
           write(use_unit,'(2X,A,f12.3,A)') &
              "| Minimal total memory requirement for three-center overlap :", &
               ((dble(n_max_loc_prodbas*n_tasks)*dble(n_basis*(n_basis+1)/2) &
                +2.0*n_basbas*n_loc_prodbas*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3)*2, &
                " Gbs"
          end if
        endif
        write(info_str,'(2X,A,2I6)') "| Number of local auxiliary basis functions on task: ", &
        & myid, n_loc_prodbas
        call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)

        i_index = 0
        map_prodbas(:,:) = 0
        do i_basis = 1, n_basbas

          if( MOD((i_basis-1)/nb_aux,npcol_aux) .eq. mypcol_aux) then
           i_index = i_index + 1
           map_prodbas(i_index,mypcol_aux+1) = i_basis
          endif
        enddo

        call sync_int_vector(map_prodbas, n_max_loc_prodbas*npcol_aux)

        call descinit( aux_sc_desc, n_basbas, n_basbas, mb_aux, nb_aux, &
                       0, 0, &
                       my_blacs_ctxt_aux, MAX(1,max_row), info )

! Calculate 2D task distribution and create BLACS context

        do nprow_aux_2d = INT(SQRT(DBLE(n_tasks)+0.5)),1,-1
          if(MOD(n_tasks,nprow_aux_2d) == 0) exit
        enddo
        npcol_aux_2d = n_tasks / nprow_aux_2d

        if(myid==0) write(use_unit,'(2x,a,i6,a,i6)') &
          '| 2D aux task distribution: ',nprow_aux_2d,' X ',npcol_aux_2d

! If this is not the first time we are distributing the tasks, the BLACS context below exists already. 
! We should free it first to avoid creating plenty of excess MPI communicators.
        if (prodbas_tasks_distributed) then
           call BLACS_Gridexit ( my_blacs_ctxt_aux_2d )
        end if

        call BLACS_GET (0,0,my_blacs_ctxt_aux_2d)
        call BLACS_Gridinit( my_blacs_ctxt_aux_2d, 'R', nprow_aux_2d, npcol_aux_2d )
        call BLACS_Gridinfo( my_blacs_ctxt_aux_2d, nprow_aux_2d, npcol_aux_2d, &
                             myprow_aux_2d, mypcol_aux_2d )

        ! We are about to create new column/row communicators.
        ! If these exist, we must have been here before. In that case, 
        ! free up the old ones first.
        if (prodbas_distributed) then
           ! Place a safeguard here to avoid ever freeing up the global communicator.
           ! In principle, this should never happen.
           if (mpi_comm_cols_aux_2d.ne.mpi_comm_global) then
              call mpi_comm_free(mpi_comm_cols_aux_2d,mpierr)
           end if
           if (mpi_comm_rows_aux_2d.ne.mpi_comm_global) then
              call mpi_comm_free(mpi_comm_rows_aux_2d,mpierr)
           end if
           prodbas_distributed = .false.
        end if

        ! Column/Row communicators
        if (.not.prodbas_distributed) then
           call mpi_comm_split(mpi_comm_global,myprow_aux_2d,myid,mpi_comm_cols_aux_2d,mpierr)
           call mpi_comm_split(mpi_comm_global,mypcol_aux_2d,myid,mpi_comm_rows_aux_2d,mpierr)
           prodbas_distributed = .true.
        end if

        nb_aux_2d = 64

        do while (nb_aux_2d * (nprow_aux_2d - 1) >= n_basbas .or. &
                  nb_aux_2d * (npcol_aux_2d - 1) >= n_basbas)
           ! There is not enough cake for anyone.
           nb_aux_2d = nb_aux_2d / 2
           if (nb_aux_2d == 0) then
              call aims_stop('System too small for scalapack (2d) - for now, use fewer CPUs?')
           end if
        end do

        write(info_str, "(2X, '| 2D block size (nb_aux_2d):',I7)") nb_aux_2d
        call localorb_info(info_str)

        max_row_2d = numroc( n_basbas, nb_aux_2d, myprow_aux_2d, 0, nprow_aux_2d )
        max_col_2d = numroc( n_basbas, nb_aux_2d, mypcol_aux_2d, 0, npcol_aux_2d )
        if (n_tasks <= 16) then
           write(info_str,"(2X,'| ',A,I5,' x',I5,' gets',I7,' x',I7,' ',A)")&
           & 'Task at', myprow_aux_2d, mypcol_aux_2d, max_row_2d, max_col_2d, &
           & 'entries of 2D aux'
           call localorb_allinfo(info_str)
        end if


        call descinit( aux_sc_desc_2d, n_basbas, n_basbas, nb_aux_2d, nb_aux_2d, &
                       0, 0, &
                       my_blacs_ctxt_aux_2d, MAX(1,max_row_2d), info )

        ! get mapping of 2D coords (row,col) to global ids

        if (allocated(global_id)) then
           ! Make sure that the dimensions are still right, else deallocate and reallocate
           if (ubound(global_id, 1) /= nprow_aux_2d-1) then
              deallocate(global_id)
           else if (ubound(global_id, 2) /= npcol_aux_2d-1) then
              deallocate(global_id)
           end if
        end if

        if (.not.allocated(global_id)) then
           allocate(global_id(0:nprow_aux_2d-1,0:npcol_aux_2d-1),stat=i_index)
           call check_allocation(i_index, 'global_id                     ')
        end if

        global_id = 0
        global_id(myprow_aux_2d,mypcol_aux_2d) = myid
        call sync_integer_vector(global_id,nprow_aux_2d*npcol_aux_2d)

        ! Select 2D distribution of ovlp_3KS matrix
        ! Here we choose the same processor grid as for scalapack

        np1_o3KS = nprow_aux_2d
        np2_o3KS = npcol_aux_2d
        myp1_o3KS = myprow_aux_2d
        myp2_o3KS = mypcol_aux_2d

!       Set owner and local index for distributed ovlp_3KS

        ! only reallocate if dimensions changed
        if (allocated(own_dim1_o3ks)) then 
           if (ubound(own_dim1_o3ks, 1) .ne. n_states) then
              deallocate(own_dim1_o3ks)
           end if
        end if
        if (allocated(loc_dim1_o3ks)) then 
           if (ubound(loc_dim1_o3ks, 1) .ne. n_states) then
              deallocate(loc_dim1_o3ks)
           end if
        end if
        if (allocated(own_dim2_o3ks)) then 
              if (ubound(own_dim2_o3ks, 1) .ne. n_states) then
           deallocate(own_dim2_o3ks)
           end if
        end if
        if (allocated(loc_dim2_o3ks)) then 
           if (ubound(loc_dim2_o3ks, 1) .ne. n_states) then
              deallocate(loc_dim2_o3ks)
           end if
        end if

        if (.not. allocated(own_dim1_o3ks)) allocate(own_dim1_o3ks(n_states))
        if (.not. allocated(loc_dim1_o3ks)) allocate(loc_dim1_o3ks(n_states))
        if (.not. allocated(own_dim2_o3ks)) allocate(own_dim2_o3ks(n_states))
        if (.not. allocated(loc_dim2_o3ks)) allocate(loc_dim2_o3ks(n_states))

        do i_index = 1, n_states
          own_dim1_o3ks(i_index) = MOD(i_index-1, np1_o3KS)
          loc_dim1_o3ks(i_index) = (i_index-1)/np1_o3KS + 1
          own_dim2_o3ks(i_index) = MOD(i_index-1, np2_o3KS)
          loc_dim2_o3ks(i_index) = (i_index-1)/np2_o3KS + 1
        enddo

!       Set MPI communicators:
!       These are the same as myprow_aux_2d/mypcol_aux_2d.
!       We use our own if the strict correspondence between 2D grid for matrices
!       and for ovlp_3KS should change in the future.

        if (o3ks_distributed) then
           ! Place a safeguard here to avoid ever freeing up the global communicator.
           ! In principle, this should never happen.
           if (mpi_comm_o3ks_1.ne.mpi_comm_global) then
              call mpi_comm_free(mpi_comm_o3ks_1,mpierr)
           end if
           if (mpi_comm_o3ks_2.ne.mpi_comm_global) then
              call mpi_comm_free(mpi_comm_o3ks_2,mpierr)
           end if
           o3ks_distributed = .false.
        end if

        ! o3KS communicators
        if (.not.o3ks_distributed) then
          call mpi_comm_split(mpi_comm_global,myp1_o3KS,myp2_o3KS,mpi_comm_o3ks_1,mpierr)
          call mpi_comm_split(mpi_comm_global,myp2_o3KS,myp1_o3KS,mpi_comm_o3ks_2,mpierr)
           o3ks_distributed = .true.
        end if


      else

         n_remain_prodbas = MOD(n_basbas,n_tasks)
            
         if (n_remain_prodbas.eq.0.and.n_basbas.ge.n_tasks) then
            n_loc_prodbas = n_basbas/n_tasks
         else
            n_loc_prodbas = n_basbas/n_tasks + 1
         endif
         n_max_loc_prodbas = n_loc_prodbas
         
        if(myid.eq.0) then
          write(use_unit,*)
          write(use_unit,'(2X,A,I6)') &
             "| Number of product basis functions per thread :", &
              n_loc_prodbas
          if (.not.use_lc_wpbeh) then
           write(use_unit,'(2X,A,f16.3,A)') &
              "| Minimal requirement for computer memory :", &
               (dble(n_loc_prodbas*n_tasks)*dble(n_basis*(n_basis+1)/2) &
                +n_basbas*n_loc_prodbas*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3, &
               " Gbs"
          else if (use_lc_wpbeh .and. .not.use_gw) then
           write(use_unit,'(2X,A,f12.3,A)') &
              "| Minimal total memory requirement for three-center overlap :", &
               ((dble(n_loc_prodbas*n_tasks)*dble(n_basis*(n_basis+1)/2) &
                +n_basbas*n_loc_prodbas*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3)*2, &
                " Gbs"
          else
           write(use_unit,'(2X,A,f12.3,A)') &
              "| Minimal total memory requirement for three-center overlap :", &
               ((dble(n_loc_prodbas*n_tasks)*dble(n_basis*(n_basis+1)/2) &
                +n_basbas*n_loc_prodbas*n_tasks)*8.d0/1.024d3/1.024d3/1.024d3)*3, &
                " Gbs"
          end if
        endif

        ! If we get here and map_prodbas is already allocate, the task distribution
        ! must have changed. Deallocate and reallocate.
        if (allocated(map_prodbas)) then
           deallocate (map_prodbas)
        endif
        if (.not. allocated(map_prodbas)) then
           allocate(map_prodbas(n_loc_prodbas,n_tasks),stat=i_index)
           call check_allocation(i_index, 'map_prodbas                   ')
        endif

        map_prodbas = 0
        i_index = 0
        do i_task = 1, n_tasks
           do i_basis = 1, n_loc_prodbas
              i_index = i_index + 1
              if (i_index > n_basbas) exit
              map_prodbas(i_basis, i_task) = i_index
           end do
        end do

        if(myid.eq.n_basbas/n_loc_prodbas)then
           n_loc_prodbas = n_basbas - n_loc_prodbas*(n_basbas/n_loc_prodbas)
        elseif(myid.gt.n_basbas/n_loc_prodbas)then
           n_loc_prodbas = 0
        endif

! end if scalapack
      endif

      if (allocated(outer_radius_prodbas) .and. (ubound(outer_radius_prodbas, 1).ne.n_loc_prodbas) ) then
         deallocate (outer_radius_prodbas)
      endif
      if (.not.allocated(outer_radius_prodbas)) then
        allocate (outer_radius_prodbas(n_loc_prodbas),stat=i_index)
        call check_allocation(i_index, 'outer_radius_prodbas          ')
      endif

      do i_loc_prodbas = 1, n_loc_prodbas
         i_basbas = map_prodbas(i_loc_prodbas, myid+1)
         if (i_basbas > 0) then
            i_basbas_fn = basbas_fn(i_basbas)
            outer_radius_prodbas(i_loc_prodbas) = &
            & field_radius_basbas_fn(i_basbas_fn)
         end if
      end do

      end subroutine prodbas_tasks_distribution
!******
  !----------------------------------------------------------------------------
  !****s* prodbas/get_basbas_to_rowcol
  !  NAME
  !    get_basbas_to_rowcol
  !  SYNOPSIS

  subroutine get_basbas_to_rowcol(distr_2d, basbas2row, basbas2col)

    !  PURPOSE
    !    Provide index mapping arrays for auxiliary basis matrices.
    !  USES

    use runtime_choices
    use scalapack_utils
    use dimensions
    use mpi_tasks, only: myid
    implicit none

    !  ARGUMENTS

    logical, intent(IN) :: distr_2d
    integer, intent(OUT) :: basbas2row(n_basbas)
    integer, intent(OUT) :: basbas2col(n_basbas)

    !  INPUTS
    !    o distr_2d -- for (use_scalapack), decide if 1D or 2D
    !  OUTPUTS
    !    o basbas2row, basbas2col -- global to local indices
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_basbas, i_loc_prodbas
    character(*), parameter :: func = 'get_basbas_to_rowcol'

    if (use_scalapack) then
       if (distr_2d) then
          do i_basbas = 1, n_basbas
             basbas2row(i_basbas) &
             & = sclpck_loc_ind(i_basbas, myprow_aux_2d, nprow_aux_2d, &
             &                  nb_aux_2d, LOCIND_ZERO)
             basbas2col(i_basbas) &
             & = sclpck_loc_ind(i_basbas, mypcol_aux_2d, npcol_aux_2d, &
             &                  nb_aux_2d, LOCIND_ZERO)
          end do
       else
          do i_basbas = 1, n_basbas
             basbas2row(i_basbas) &
             & = sclpck_loc_ind(i_basbas, myprow_aux, nprow_aux, &
             &                  mb_aux, LOCIND_ZERO)
             basbas2col(i_basbas) &
             & = sclpck_loc_ind(i_basbas, mypcol_aux, npcol_aux, &
             &                  nb_aux, LOCIND_ZERO)
          end do
       end if
    else
       do i_basbas = 1, n_basbas
          basbas2row(i_basbas) = i_basbas
       end do
       basbas2col = 0
       do i_loc_prodbas = 1, n_loc_prodbas
          i_basbas = map_prodbas(i_loc_prodbas, myid+1)
          if (i_basbas > 0) then
             basbas2col(i_basbas) = i_loc_prodbas
          end if
       end do
    end if

  end subroutine get_basbas_to_rowcol
  !******
!------------------------------------------------------------------------------
!****s* prodbas/allocate_basis_pairs
!  NAME
!    allocate_basis_pairs
!  SYNOPSIS

      subroutine allocate_basis_pairs &
                ( )
!  PURPOSE
!  Allocates basis pairs.
!
!  USES
         use dimensions
         use mpi_tasks, only: check_allocation
         implicit none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


         integer:: info


         if (.not.allocated(n_nghr_basis)) then
            allocate (n_nghr_basis(n_basis),stat=info)
            call check_allocation(info, 'n_nghr_basis                  ')
         endif

         if (.not.allocated(basis_nghr)) then
            allocate (basis_nghr(n_centers_basis_T,n_basis),stat=info)
            call check_allocation(info, 'basis_nghr                    ')
         endif

         if (.not.allocated(n_nghr_basis_original)) then
            allocate (n_nghr_basis_original(n_basis),stat=info)
            call check_allocation(info, 'n_nghr_basis_original         ')
         endif

         if (.not.allocated(basis_nghr_original)) then
            allocate (basis_nghr_original(n_centers_basis_T,n_basis),stat=info)
            call check_allocation(info, 'basis_nghr_original           ')
         endif
         if (.not.allocated(atom_pairs)) then
            allocate (atom_pairs(n_atoms,n_atoms),stat=info)
            call check_allocation(info, 'atom_pairs                    ')
         endif

        end subroutine allocate_basis_pairs
!******

!----------------------------------------------------------------------
!****s* prodbas/cleanup_basbas
!  NAME
!    cleanup_basbas
!  SYNOPSIS

        subroutine cleanup_basbas &
        ( )

!  PURPOSE
!  Subroutine cleanup_basis deallocates all basis-related arrays
!
!  USES
         use dimensions
         use pbc_lists, only: n_cells_pairs

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!        if (allocated(basbas_wave_spl)) then
!          deallocate(basbas_wave_spl)
!        end if
!        if (allocated(basis_aug_deriv_spl)) then
!          deallocate(basis_aug_deriv_spl)
!        end if
!       if (allocated(outer_radius_prodbas)) then
!          deallocate(outer_radius_prodbas)
!        end if
!
!        if (allocated(basbas_atom)) then
!          deallocate(basbas_atom)
!        end if
!        if (allocated(basbas_l)) then
!          deallocate(basbas_l)
!        end if
!        if (allocated(basbas_m)) then
!          deallocate(basbas_m)
!        end if
!        if (allocated(basbas_fn)) then
!          deallocate(basbas_fn)
!        end if
!        if (allocated(map_prodbas)) then
!          deallocate(map_prodbas)
!        end if
        prodbas_tasks_distributed = .false.
        previous_n_basbas = -1
        previous_prodbas_nb = -1
        previous_n_states = -1

        prodbas_tasks_pbc_distributed = .false.
        previous_n_basbas_supercell = -1

        prodbas_distributed = .false.

        o3ks_distributed = .false.
!        ! JW: For whatever reason the above are not deallocated, it probably
!        !     also applies to basbasfn_species, basbasfn_l,
!        !     field_radius_basbas_fn, charge_radius_babbas_fn,
!        !     multipole_basbas_fn, ... which I also do not deallocate.
        if (allocated(n_nghr_basis)) then
           deallocate (n_nghr_basis)
        endif
        if (allocated(basis_nghr)) then
           deallocate (basis_nghr)
        endif
        if (allocated(atom_pairs)) then
           deallocate (atom_pairs)
        endif
        if (allocated(n_nghr_basis_original)) then ! atom bsse
           deallocate (n_nghr_basis_original)
        endif
        if (allocated(basis_nghr_original)) then ! atom bsse
           deallocate (basis_nghr_original)
        endif

        ! FK: For relaxation it is CRITICAL to deallocate arrays that depend on the size of the product basis
        !     since this size changes with changing geometry!
        !     For simplicity and safety reasons I deallocate every array that is allocated in this module
        !     if not done above.
        if (allocated(basbas_wave_spl)) then
           deallocate (basbas_wave_spl)
        endif
        if (allocated(v_times_radialbasbas_spl)) then
           deallocate (v_times_radialbasbas_spl)
        endif
        if (allocated(basbas_atom)) then
           deallocate (basbas_atom)
        endif
        if (allocated(basbas_l)) then
           deallocate (basbas_l)
        endif
        if (allocated(basbas_m)) then
           deallocate (basbas_m)
        endif
        if (allocated(basbas_fn)) then
           deallocate (basbas_fn)
        endif
        if (allocated(basbasfn_species)) then
           deallocate (basbasfn_species)
        endif
        if (allocated(basbasfn_l)) then
           deallocate (basbasfn_l)
        endif
        if (allocated(charge_radius_basbas_fn)) then
           deallocate (charge_radius_basbas_fn)
        endif
        if (allocated(field_radius_basbas_fn)) then
           deallocate (field_radius_basbas_fn)
        endif
        if (allocated(multipole_basbas_fn)) then
           deallocate (multipole_basbas_fn)
        endif
        if (allocated(Lsp2n_basbas_fnLsp)) then
           deallocate (Lsp2n_basbas_fnLsp)
        endif
        if (allocated(Lsp2basbas_fn)) then
           deallocate (Lsp2basbas_fn)
        endif
        if (allocated(Lsp2basbas_sp)) then
           deallocate (Lsp2basbas_sp)
        endif
        if (allocated(atom2basbas_off)) then
           deallocate (atom2basbas_off)
        endif
        if (allocated(sp2n_basbas_sp)) then
           deallocate (sp2n_basbas_sp)
        endif

        end subroutine cleanup_basbas
        


!******
!---------------------------------------------------------------------
      end module prodbas
